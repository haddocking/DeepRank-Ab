#!/usr/bin/env python3
"""
neuralnet.py

Training harness for DeepRank-AB. Builds/loads HDF5 datasets, instantiates a
user-provided GNN, and handles train/valid/test loops, metrics, plotting, and
checkpointing. Supports regression or classification, optional community
clustering, and an EMA evaluation model.
"""

import sys
import os
import time
import warnings
import torch
import numpy as np
import h5py
import inspect
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from tools.FocalLoss import BMCLoss
from DataSet import HDF5DataSet, DivideDataSet, PreCluster
from Metrics import Metrics
from torch.optim.lr_scheduler import StepLR, ReduceLROnPlateau
from copy import deepcopy

# deprecation noise from deps
warnings.filterwarnings("ignore", message="TypedStorage is deprecated")
warnings.filterwarnings("ignore", message="'dropout_adj' is deprecated")


# Simple EMA wrapper 

class EMA:
    """
    Exponential moving average wrapper.

    - constructor: EMA(model, beta, update_after_step, update_every)
    - property: ema_model (smoothed copy)
    - methods: update(), state_dict(), load_state_dict()
    - attribute: update_after_step (step warmup)
    """
    def __init__(self, model, beta=0.999, update_after_step=0, update_every=1):
        self.ema_model = deepcopy(model).eval()
        for p in self.ema_model.parameters():
            p.requires_grad_(False)

        self.model = model
        self.beta = float(beta)
        self.update_after_step = int(update_after_step)
        self.update_every = int(update_every)
        self.step = 0

    @torch.no_grad()
    def update(self):
        """Update EMA weights after optimizer.step()."""
        self.step += 1
        if self.step <= self.update_after_step:
            return
        if (self.step % self.update_every) != 0:
            return

        model_params = dict(self.model.named_parameters())
        ema_params   = dict(self.ema_model.named_parameters())
        for name, p_ema in ema_params.items():
            p_src = model_params[name]
            if p_ema.dtype.is_floating_point:
                p_ema.data.mul_(self.beta).add_(p_src.detach().data, alpha=1.0 - self.beta)
            else:
                p_ema.data.copy_(p_src.data)

        model_buffers = dict(self.model.named_buffers())
        ema_buffers   = dict(self.ema_model.named_buffers())
        for name, b_ema in ema_buffers.items():
            b_src = model_buffers[name]
            if b_ema.dtype.is_floating_point:
                b_ema.data.mul_(self.beta).add_(b_src.detach().data, alpha=1.0 - self.beta)
            else:
                b_ema.data.copy_(b_src.data)

    def state_dict(self):
        return {
            "ema_model": self.ema_model.state_dict(),
            "beta": self.beta,
            "update_after_step": self.update_after_step,
            "update_every": self.update_every,
            "step": self.step,
        }

    def load_state_dict(self, state):
        self.beta = float(state.get("beta", self.beta))
        self.update_after_step = int(state.get("update_after_step", self.update_after_step))
        self.update_every = int(state.get("update_every", self.update_every))
        self.step = int(state.get("step", 0))
        self.ema_model.load_state_dict(state["ema_model"], strict=False)



# Trainer

class NeuralNet(object):
    @staticmethod
    def _fmt(v, spec=".4f"):
        import math
        if v is None:
            return "nan"
        try:
            if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
                return "nan"
            return format(v, spec)
        except Exception:
            return "nan"

    def __init__(
        self,
        database,
        Net,
        node_feature=["type", "polarity", "bsa"],
        edge_feature=["dist"],
        target="irmsd",
        lr=0.01,
        scheduler_type=None,        # "plateau" | "step" | None
        scheduler_params=None,
        batch_size=32,
        percent=[1.0, 0.0],
        database_eval=None,
        index=None,
        index_eval=None,
        class_weights=None,
        task=None,
        classes=[0, 1],
        threshold=None,
        device_name=None,
        num_workers=1,
        pretrained_model=None,
        shuffle=True,
        outdir="./",
        cluster_nodes="mcl",
        transform_sigmoid=False,
        # EMA
        use_ema: bool = True,
        ema_decay: float = 0.999,
        ema_warmup_steps: int = 0,
        ema_use_for_eval: bool = True,
    ):
        """Configure training, data, model, and (optionally) EMA."""
        sig = inspect.signature(Net.__init__)
        n_params = len(sig.parameters) - 1
        if n_params == 3:
            self.use_equivariant = False
        elif n_params == 5:
            self.use_equivariant = True
        else:
            raise ValueError(f"Unsupported signature for {Net.__name__}: expected 3 or 5 parameters")

        # args to attributes
        self.num_workers = num_workers
        self.device_name = device_name

        # EMA
        self.ema = None
        self.use_ema = use_ema
        self.ema_decay = ema_decay
        self.ema_warmup_steps = ema_warmup_steps
        self.ema_use_for_eval = ema_use_for_eval
        self.global_step = 0

        # LR/scheduler
        self.lr = lr
        self.scheduler = None
        self.scheduler_type = scheduler_type
        self.scheduler_params = scheduler_params or {}

        # logs
        self.train_out, self.train_y = [], []
        self.valid_out, self.valid_y = [], []
        self.test_out,  self.test_y  = [], []

        self.train_acc = []
        self.train_loss = []
        self.valid_acc = []
        self.valid_loss = []
        self.train_auc = []
        self.valid_auc = []
        self.train_recall = []
        self.valid_recall = []
        self.train_mse = []
        self.valid_mse = []
        self.train_rmse = []
        self.valid_rmse = []
        self.train_r2 = []
        self.valid_r2 = []
        self.train_prec = []
        self.valid_prec = []
        self.train_f1 = []
        self.valid_f1 = []

        # ---- case 1: fresh training ----
        if pretrained_model is None:
            # mirror constructor args onto self (except handled ones)
            for k, v in dict(locals()).items():
                if k not in ["self", "database", "Net", "database_eval"]:
                    self.__setattr__(k, v)

            # infer task if missing
            if self.task is None:
                if self.target in ["irmsd", "lrmsd", "fnat", "dockQ"]:
                    self.task = "reg"
                elif self.target in ["bin_class", "capri_classes"]:
                    self.task = "class"
                else:
                    raise ValueError(
                        "Missing task. Set task='class' or task='reg'."
                    )

            if self.task == "class" and self.threshold is None:
                print(f"threshold set to {self.classes[1]}")
                self.threshold = self.classes[1]
            if self.task == "reg" and self.threshold is None:
                print("threshold set to 0.23")
                self.threshold = 0.23

            self.load_model(database, Net, database_eval)

        # ---- case 2: load checkpoint ----
        else:
            self.load_params(pretrained_model)
            self.outdir = outdir
            self.load_pretrained_model(database, Net)

    def _net_for_eval(self):
        """Pick model for eval/test (EMA if enabled)."""
        use_ema_now = (
            self.use_ema
            and self.ema_use_for_eval
            and getattr(self, "ema", None) is not None
        )
        return self.ema.ema_model if use_ema_now else self.model

    # -------------------- checkpoint path: data & model --------------------
    def load_pretrained_model(self, database, Net):
        """Load datasets and restore a pretrained model + optimizer (+EMA if present)."""
        test_dataset = HDF5DataSet(
            name="Test",
            root="./",
            database=database,
            node_feature=self.node_feature,
            edge_feature=self.edge_feature,
            target=self.target,
            clustering_method=self.cluster_nodes,
        )
        #PreCluster(test_dataset, method=self.cluster_nodes)
        self.test_loader = DataLoader(test_dataset, num_workers=self.num_workers)


        self.put_model_to_device(test_dataset, Net)
        self.set_loss()

        self.optimizer = torch.optim.Adam([
            {'params': self.model.parameters(), 'lr': self.lr},
            {'params': self.loss.noise_sigma,   'lr': 1e-4, 'name': 'noise_sigma'}
        ])

        if self.scheduler_type == "step":
            self.scheduler = StepLR(self.optimizer, **self.scheduler_params)
        elif self.scheduler_type == "plateau":
            self.scheduler = ReduceLROnPlateau(self.optimizer, mode='min', **self.scheduler_params)
        else:
            self.scheduler = None

        self.optimizer.load_state_dict(self.opt_loaded_state_dict)
        self.model.load_state_dict(self.model_load_state_dict)

        if self.use_ema:
            if not getattr(self, "ema", None):
                self.ema = EMA(
                    self.model,
                    beta=self.ema_decay,
                    update_after_step=max(0, int(self.ema_warmup_steps or 0)),
                    update_every=1,
                )
            if getattr(self, "ema_loaded_state_dict", None) is not None:
                self.ema.load_state_dict(self.ema_loaded_state_dict)
            else:
                self.ema.ema_model.load_state_dict(self.model.state_dict(), strict=False)
        if self.use_ema and self.ema is not None:
            self.ema.update_after_step = int(max(0, self.ema_warmup_steps))

    # -------------------- fresh training path: data & model --------------------
    def load_model(self, database, Net, database_eval):
        """Build datasets/loaders and create the model."""
        dataset = HDF5DataSet(
            name="Train",
            root="./",
            database=database,
            index=self.index,
            node_feature=self.node_feature,
            edge_feature=self.edge_feature,
            target=self.target,
            clustering_method=self.cluster_nodes,
        )

        if self.cluster_nodes is not None:
            if self.cluster_nodes in ("mcl", "louvain"):
                print("loading clusters")
            else:
                raise ValueError(
                    "Invalid node clustering method. Set 'mcl', 'louvain', or None."
                )

        train_dataset, valid_dataset = DivideDataSet(dataset, percent=self.percent)

        #PreCluster(train_dataset, method=self.cluster_nodes)
        self.train_loader = DataLoader(
            train_dataset,
            batch_size=self.batch_size,
            shuffle=self.shuffle,
            num_workers=self.num_workers,
        )


        if self.percent[1] > 0.0:
            self.valid_loader = DataLoader(
                valid_dataset,
                batch_size=self.batch_size,
                shuffle=False,
                num_workers=self.num_workers,
            )


        if database_eval is not None:
            print("loading independent validation set")
            valid_dataset = HDF5DataSet(
                name="Evaluation",
                root="./",
                database=database_eval,
                index=self.index,
                node_feature=self.node_feature,
                edge_feature=self.edge_feature,
                target=self.target,
                clustering_method=self.cluster_nodes,
            )
            if self.cluster_nodes in ("mcl", "louvain"):
                print("loading clusters for evaluation")
                #PreCluster(valid_dataset, method=self.cluster_nodes)

            self.valid_loader = DataLoader(
                valid_dataset,
                num_workers=self.num_workers,
                batch_size=self.batch_size,
                shuffle=False,
            )
            print("independent validation set loaded")
        else:
            print("no independent validation set")

        self.put_model_to_device(dataset, Net)
        self.set_loss()

        self.optimizer = torch.optim.Adam([
            {'params': self.model.parameters(), 'lr': self.lr},
            {'params': self.loss.noise_sigma,   'lr': 1e-4, 'name': 'noise_sigma'}
        ])

        print("optimizer param groups")
        for i, g in enumerate(self.optimizer.param_groups):
            name = g.get('name', f'group_{i}')
            print(f"group {i}: {name}, lr {g['lr']}, nparams {len(g['params'])}")

        if self.scheduler_type == "step":
            self.scheduler = StepLR(self.optimizer, **self.scheduler_params)
        elif self.scheduler_type == "plateau":
            self.scheduler = ReduceLROnPlateau(self.optimizer, mode='min', **self.scheduler_params)
        else:
            self.scheduler = None

        if self.use_ema:
            steps_per_epoch = len(self.train_loader)
            if self.ema_warmup_steps is None or self.ema_warmup_steps <= 0:
                self.ema_warmup_steps = 2 * steps_per_epoch
            print(f"ema warmup steps {self.ema_warmup_steps}")

        if self.use_ema and self.ema is not None:
            self.ema.update_after_step = int(max(0, self.ema_warmup_steps))

    def put_model_to_device(self, dataset, Net):
        """Infer input dims from a sample and build the network on the chosen device."""
        self.device = (torch.device(self.device_name)
                       if self.device_name
                       else torch.device('cuda' if torch.cuda.is_available() else 'cpu'))

        example = dataset.get(0)

        if self.use_equivariant:
            node_s, node_v = example.x
            edge_s, edge_v = example.edge_attr
            in_ns, in_nv = node_s.shape[1], node_v.shape[1]
            in_es, in_ev = edge_s.shape[1], edge_v.shape[1]

            if self.task == "reg":
                out_dim = 1
            else:
                self.classes_to_idx = {c: i for i, c in enumerate(self.classes)}
                self.idx_to_classes = {i: c for i, c in enumerate(self.classes)}
                out_dim = len(self.classes)

            self.model = Net(in_ns, in_nv, in_es, in_ev, out_dim=out_dim).to(self.device)
        else:
            in_feats   = example.x.shape[1]
            edge_feats = example.edge_attr.shape[1]


            if self.task == "reg":
                out_dim = 1
            else:
                self.classes_to_idx = {c: i for i, c in enumerate(self.classes)}
                self.idx_to_classes = {i: c for i, c in enumerate(self.classes)}
                out_dim = len(self.classes)

            self.model = Net(in_feats, out_dim, edge_feats).to(self.device)

        if self.use_ema:
            self.ema = EMA(
                self.model,
                beta=self.ema_decay,
                update_after_step=max(0, int(self.ema_warmup_steps or 0)),
                update_every=1,
            )

    def set_loss(self):
        """Set loss: BMC for regression, CrossEntropy for classification (optional class weights)."""
        if self.task == "reg":
            self.loss = BMCLoss(init_noise_sigma=0.5, learn_noise=True, detach_scale=False)

        elif self.task == "class":
            self.weights = None
            if self.class_weights is True:
                targets_all = []
                for batch in self.train_loader:
                    targets_all.append(batch.y)
                targets_all = torch.cat(targets_all).squeeze().tolist()
                self.weights = torch.tensor(
                    [targets_all.count(i) for i in self.classes], dtype=torch.float32
                )
                print(f"class counts {self.weights}")
                self.weights = 1.0 / self.weights
                self.weights = self.weights / self.weights.sum()
                print(f"class weights {self.weights}")
            self.loss = nn.CrossEntropyLoss(weight=self.weights, reduction="mean")

    # -------------------- pretty console printing helpers --------------------
    def _print_epoch_block(self, phase, epoch, vals, dt_seconds):
        """
        Pretty-print a per-epoch summary line for TRAIN/VALID in the requested style.
        """
        # Build metric line
        line = (
            f"Epoch [{epoch:04d}] {phase} "
            f"loss {self._fmt(vals['loss'], '.6e')} | "
            f"acc {self._fmt(vals['acc'])} | "
            f"recall {self._fmt(vals['recall'])} | "
            f"auc {self._fmt(vals['auc'])} | "
            f"mse {self._fmt(vals['mse'], '.4e')} | "
            f"rmse {self._fmt(vals['rmse'], '.4e')} | "
            f"r2 {self._fmt(vals['r2'])} | "
            f"prec {self._fmt(vals['prec'])} | "
            f"f1 {self._fmt(vals['f1'])}"
        )
        print(line)
        print(f" time {dt_seconds:.2f}s")

    def train(
        self,
        nepoch=1,
        validate=False,
        save_model="best",
        hdf5="train_data.hdf5",
        save_epoch="intermediate",
        save_every=5,
        early_stopping_patience=None,
        min_delta=1e-4,
        topk=1
    ):
        """Train with optional validation/early stopping; write per-epoch HDF5 logs.
        Also prints compact per-epoch TRAIN/VALID summaries with timings.
        """
        best_loss = float("inf")
        best_mse = float("inf")
        self.best_epoch = None
        topk_mse = []
        epochs_no_improve = 0
        fname = self.update_name(hdf5, self.outdir)

        with h5py.File(fname, "w") as self.f5:
            self.nepoch = nepoch
            self.data = {}

            # Header line (matches requested style)
            print(f"[2] Training for {nepoch} epochs…")

            for epoch in range(1, nepoch + 1):
                # ---------------- train ----------------
                self.model.train()

                total_loss = 0.0
                total_examples = 0
                train_out, train_raw, train_y, train_mol = [], [], [], []

                t0 = time.time()
                for batch in self.train_loader:
                    batch = batch.to(self.device)
                    bsize = batch.y.size(0)

                    self.optimizer.zero_grad()
                    pred = self.model(batch)
                    pred, batch.y = self.format_output(pred, batch.y)
                    batch.y = batch.y.to(self.device)

                    loss = self.loss(pred, batch.y)
                    loss_val = loss.item()
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
                    self.optimizer.step()

                    # EMA
                    self.global_step += 1
                    if self.use_ema and self.ema is not None:
                        self.ema.update()

                    total_loss += loss_val * bsize
                    total_examples += bsize

                    if self.task == "class":
                        probs = F.softmax(pred.cpu().detach(), dim=1)
                        train_raw += probs.tolist()
                        preds = np.argmax(probs, axis=1)
                    else:
                        preds = pred.cpu().detach().reshape(-1)
                        train_raw += preds.tolist()

                    train_out += preds.tolist()
                    train_y += batch.y.tolist()
                    train_mol += batch["mol"]

                train_dt = time.time() - t0
                avg_train_loss = total_loss / total_examples if total_examples > 0 else 0.0
                self.train_loss.append(avg_train_loss)
                self.train_out, self.train_y = train_out, train_y

                train_metrics = self.get_metrics("train", self.threshold)
                train_vals = {
                    "loss": avg_train_loss,
                    "acc": train_metrics.accuracy,
                    "recall": train_metrics.sensitivity,
                    "auc": train_metrics.auc(),
                    "mse": train_metrics.mean_squared_error,
                    "rmse": train_metrics.root_mean_squared_error,
                    "r2": train_metrics.r2_score,
                    "prec": train_metrics.precision,
                }
                print(train_vals)
                if train_vals["prec"] is not None and train_vals["recall"] is not None:
                    p, r = train_vals["prec"], train_vals["recall"]
                    train_vals["f1"] = 2 * p * r / (p + r) if (p + r) > 0 else float("nan")
                else:
                    train_vals["f1"] = float("nan")

                train_payload = train_vals.copy()
                train_payload["outputs"]     = train_out
                train_payload["raw_outputs"] = train_raw
                train_payload["targets"]     = train_y
                train_payload["mol"]         = train_mol
                self.data["train"] = train_payload

                for key in ["acc", "auc", "recall", "prec", "f1"]:
                    getattr(self, f"train_{key}").append(train_vals[key])
                for key in ["mse", "rmse", "r2"]:
                    getattr(self, f"train_{key}").append(train_vals[key])

                # Pretty TRAIN print
                self._print_epoch_block("TRAIN", epoch, train_vals, train_dt)

                # ---------------- valid / early stop ----------------
                if validate and hasattr(self, "valid_loader") and self.valid_loader is not None:
                    vt0 = time.time()
                    _out, _y, avg_val_loss, _payload = self.eval(self.valid_loader)
                    val_dt = time.time() - vt0

                    self.valid_out, self.valid_y = _out, _y

                    val_metrics = self.get_metrics("eval", self.threshold)
                    val_vals = {
                        "loss": avg_val_loss,
                        "acc": val_metrics.accuracy,
                        "recall": val_metrics.sensitivity,
                        "auc": val_metrics.auc(),
                        "mse": val_metrics.mean_squared_error,
                        "rmse": val_metrics.root_mean_squared_error,
                        "r2": val_metrics.r2_score,
                        "prec": val_metrics.precision,
                    }

                    if val_vals["prec"] is not None and val_vals["recall"] is not None:
                        p, r = val_vals["prec"], val_vals["recall"]
                        val_vals["f1"] = 2 * p * r / (p + r) if (p + r) > 0 else float("nan")
                    else:
                        val_vals["f1"] = float("nan")

                    self.valid_loss.append(avg_val_loss)
                    for key in ["acc", "auc", "recall", "prec", "f1"]:
                        getattr(self, f"valid_{key}").append(val_vals[key])
                    for key in ["mse", "rmse", "r2"]:
                        getattr(self, f"valid_{key}").append(val_vals[key])

                    # top-K by MSE save
                    if save_model == "topk_mse":
                        cur_mse = val_vals["mse"]
                        if cur_mse is not None and np.isfinite(cur_mse):
                            qualifies = (len(topk_mse) < topk) or (cur_mse < max(topk_mse, key=lambda t: t[0])[0])
                            if qualifies:
                                filename = (
                                    f"top_mse_ep{epoch}_mse{cur_mse:.4e}"
                                    f"_t{self.task}_y{self.target}_b{self.batch_size}"
                                    f"_e{nepoch}_lr{self.lr:.0e}.pth.tar"
                                )
                                outpath = os.path.join(self.outdir, filename)
                                self.save_model(filename=outpath)
                                topk_mse.append((cur_mse, epoch, outpath))
                                topk_mse.sort(key=lambda t: t[0])
                                while len(topk_mse) > topk:
                                    worst_mse, worst_ep, worst_path = topk_mse.pop()
                                    try:
                                        os.remove(worst_path)
                                        # Match requested style
                                        print(f"Removed {worst_path}")
                                    except OSError:
                                        pass
                                # no extra print here; keep console lean

                    # Best-MSE banner and early stopping book-keeping
                    if (val_vals["mse"] is not None and np.isfinite(val_vals["mse"])
                            and (val_vals["mse"] < best_mse - min_delta)):
                        prev = float("inf") if not np.isfinite(best_mse) else best_mse
                        print(f"→ [BEST MSE] Epoch {epoch}: mse {val_vals['mse']:.4e} < prev_best {prev if prev==float('inf') else f'{prev:.4e}'}")
                        best_mse = val_vals["mse"]
                        self.best_epoch = epoch
                        if early_stopping_patience is not None:
                            epochs_no_improve = 0
                    else:
                        if early_stopping_patience is not None:
                            epochs_no_improve += 1
                            # keep the classic “no improvement” trace to help debugging
                            print(f"no improvement {epochs_no_improve}/{early_stopping_patience}")
                            if epochs_no_improve >= early_stopping_patience:
                                print(f"early stopping at epoch {epoch}")
                                # print final VALID block before break
                                self._print_epoch_block("VALID", epoch, val_vals, val_dt)
                                # step schedulers appropriately before break if needed
                                if self.scheduler_type == "plateau" and self.scheduler:
                                    self.scheduler.step(avg_val_loss)
                                break

                    if self.scheduler_type == "plateau" and self.scheduler:
                        self.scheduler.step(avg_val_loss)

                    # Merge payload and store
                    val_payload = val_vals.copy()
                    val_payload["outputs"]     = _payload["outputs"]
                    val_payload["raw_outputs"] = _payload["raw_outputs"]
                    val_payload["targets"]     = _payload["targets"]
                    val_payload["mol"]         = _payload["mol"]
                    self.data["eval"] = val_payload

                    # Pretty VALID print
                    self._print_epoch_block("VALID", epoch, val_vals, val_dt)

                if self.scheduler_type == "step" and self.scheduler:
                    self.scheduler.step()

                # Explicit LR line
                print(f"Epoch {epoch}: lr = {self.optimizer.param_groups[0]['lr']:.2e}")

                # Per-epoch HDF5 export cadence
                if save_epoch == "all" or epoch == nepoch or (
                    save_epoch == "intermediate" and epoch % save_every == 0
                ):
                    self._export_epoch_hdf5(epoch, self.data)

        if save_model == "best" and self.best_epoch is not None:
            filename = (
                f"best_t{self.task}_y{self.target}"
                f"_b{self.batch_size}_e{nepoch}"
                f"_lr{self.lr:.0e}"
                f"_ep{self.best_epoch}.pth.tar"
            )
            outpath = os.path.join(self.outdir, filename)
            self.save_model(filename=outpath)
            print(f"saved best model from epoch {self.best_epoch}")

        if save_model == "last":
            last_fn = f"t{self.task}_y{self.target}_b{self.batch_size}_e{nepoch}_lr{self.lr:.0e}.pth.tar"
            self.save_model(filename=os.path.join(self.outdir, last_fn))

        if save_model == "topk_mse":
            print(f"saved top {topk} by val mse:")
            for rank, (mse, ep, path) in enumerate(topk_mse, 1):
                print(f"#{rank}: epoch {ep}, mse {mse:.4e} -> {path}")

    def test(self, database_test=None, threshold=0.23, hdf5="test_data.hdf5"):
        """Evaluate on a test set and export results to HDF5."""
        self.data = {}
        fname = self.update_name(hdf5, self.outdir)
        with h5py.File(fname, "w") as self.f5:
            if database_test is not None:
                test_dataset = HDF5DataSet(
                    name="Test",
                    root="./",
                    database=database_test,
                    node_feature=self.node_feature,
                    edge_feature=self.edge_feature,
                    target=self.target,
                    clustering_method=self.cluster_nodes,
                )
                print("test set loaded")
                #PreCluster(test_dataset, method=self.cluster_nodes)
                self.test_loader = DataLoader(test_dataset, num_workers=self.num_workers, batch_size=self.batch_size)

            net = self._net_for_eval()
            net.eval()

            total_loss = 0.0
            total_examples = 0
            out, raw_outputs, y = [], [], []
            data = {"outputs": [], "raw_outputs": [], "targets": [], "mol": [], "loss": 0}
            for batch in self.test_loader:
                batch = batch.to(self.device)
                bsize = batch.y.size(0)
                #bszie = self.batch_size

                pred = net(batch)
                pred, batch.y = self.format_output(pred, batch.y)

                batch.y = batch.y.to(self.device)
                loss_val = self.loss(pred, batch.y).item()
                total_loss += loss_val * bsize
                total_examples += bsize

                if self.task == "class":
                    probs = F.softmax(pred.cpu().detach(), dim=1)
                    raw_outputs += probs.tolist()
                    preds = np.argmax(probs, axis=1)
                else:
                    preds = pred.cpu().detach().reshape(-1)
                    raw_outputs += preds.tolist()

                out += preds.tolist()
                y += batch.y.tolist()
                data["mol"] += batch["mol"]

            avg_test_loss = total_loss / total_examples if total_examples > 0 else 0.0
            self.test_loss = avg_test_loss
            data["loss"] = avg_test_loss

            if self.task == "class":
                data["targets"] = [self.idx_to_classes[i] for i in y]
                data["outputs"] = [self.idx_to_classes[i] for i in out]
            else:
                data["targets"] = y
                data["outputs"] = out

            data["raw_outputs"] = raw_outputs
            self.test_out = out
            self.test_y = y
            self.test_acc = None if not y else self.get_metrics("test", threshold).accuracy

            self.data["test"] = data
            self._export_epoch_hdf5(0, self.data)

            test_metrics = self.get_metrics("test", threshold)
            data["acc"]    = test_metrics.accuracy
            data["recall"] = test_metrics.sensitivity
            data["auc"]    = test_metrics.auc()
            data["mse"]    = test_metrics.mean_squared_error
            data["rmse"]   = test_metrics.root_mean_squared_error
            data["r2"]     = test_metrics.r2_score
            data["prec"]   = test_metrics.precision
            if data["prec"] is not None and data["recall"] is not None:
                p, r = data["prec"], data["recall"]
                data["f1"] = 2 * p * r / (p + r) if (p + r) > 0 else float("nan")
            else:
                data["f1"] = float("nan")

    def predict(self, database_test=None, hdf5="predictions.hdf5"):
        """
        Run inference on a dataset and export raw predictions to HDF5.
        Does not require targets (y) to exist.
        """
        self.data = {}
        fname = self.update_name(hdf5, self.outdir)
        with h5py.File(fname, "w") as self.f5:
            if database_test is not None:
                pred_dataset = HDF5DataSet(
                    name="Predict",
                    root="./",
                    database=database_test,
                    node_feature=self.node_feature,
                    edge_feature=self.edge_feature,
                    target=self.target,  # still pass, but unused if missing
                    clustering_method=self.cluster_nodes,
                )
                self.pred_loader = DataLoader(
                    pred_dataset,
                    num_workers=self.num_workers,
                    batch_size=self.batch_size,
                    shuffle=False
                )

            net = self._net_for_eval()
            net.eval()

            out, raw_outputs, mols = [], [], []

            for batch in self.pred_loader:
                batch = batch.to(self.device)

                # Some HDF5s might lack targets (y=None)
                # so we don't touch batch.y at all
                with torch.no_grad():
                    pred = net(batch)

                if self.task == "class":
                    probs = F.softmax(pred.cpu(), dim=1)
                    raw_outputs += probs.tolist()
                    preds = np.argmax(probs, axis=1)
                else:
                    preds = pred.cpu().reshape(-1)
                    raw_outputs += preds.tolist()

                out += preds.tolist()
                mols += batch["mol"]

            # Save predictions to HDF5
            data = {
                "outputs": out,
                "raw_outputs": raw_outputs,
                "mol": mols,
            }
            self.data["pred"] = data
            self._export_epoch_hdf5(0, self.data)

            return out, mols

    def eval(self, loader):
        """Evaluate on a loader; return predictions, targets, avg loss, and payload."""
        net = self._net_for_eval()
        net.eval()

        loss_func = self.loss
        total_loss = 0.0
        total_examples = 0

        out, raw_outputs, y = [], [], []
        data = {"outputs": [], "raw_outputs": [], "targets": [], "mol": [], "loss": 0}

        for data_batch in loader:
            data_batch = data_batch.to(self.device)
            bsize = data_batch.y.size(0)

            pred = net(data_batch)
            pred, data_batch.y = self.format_output(pred, data_batch.y)

            if data_batch.y is not None:
                data_batch.y = data_batch.y.to(self.device)
                batch_loss = loss_func(pred, data_batch.y).item()
                total_loss += batch_loss * bsize
                total_examples += bsize
                y += data_batch.y.tolist()

            if self.task == "class":
                probs = F.softmax(pred.cpu().detach(), dim=1)
                raw_outputs += probs.tolist()
                preds = np.argmax(probs, axis=1)
            else:
                preds = pred.cpu().detach().reshape(-1)
                raw_outputs += preds.tolist()

            out += preds.tolist()
            data["mol"] += data_batch["mol"]

        if self.task == "class":
            data["targets"] = [self.idx_to_classes[i] for i in y]
            data["outputs"] = [self.idx_to_classes[i] for i in out]
        else:
            data["targets"] = y
            data["outputs"] = out

        data["raw_outputs"] = raw_outputs
        avg_loss = total_loss / total_examples if total_examples > 0 else 0.0
        data["loss"] = avg_loss
        return out, y, avg_loss, data

    def _epoch(self, epoch):
        """Single training epoch; returns predictions, targets, running loss, and payload."""
        running_loss = 0
        out, raw_outputs, y = [], [], []
        data = {"outputs": [], "raw_outputs": [], "targets": [], "mol": [], "loss": 0}

        for data_batch in self.train_loader:
            data_batch = data_batch.to(self.device)
            self.optimizer.zero_grad()
            pred = self.model(data_batch)
            pred, data_batch.y = self.format_output(pred, data_batch.y)

            pred = pred.to(self.device)
            data_batch.y = data_batch.y.to(self.device)

            loss = self.loss(pred, data_batch.y)
            running_loss += loss.detach().item()
            loss.backward()
            self.optimizer.step()

            try:
                y += data_batch.y.tolist()
            except ValueError:
                print("provide target values (y) for the training set")

            if self.task == "class":
                pred = F.softmax(pred, dim=1).cpu().detach()
                raw_outputs += pred.tolist()
                pred = np.argmax(pred, axis=1)
            else:
                pred = pred.cpu().detach().reshape(-1)
                raw_outputs += pred.tolist()

            out += pred.tolist()
            data["mol"] += data_batch["mol"]

        if self.task == "class":
            data["targets"] += [self.idx_to_classes[x] for x in y]
            data["outputs"] += [self.idx_to_classes[x] for x in out]
        else:
            data["targets"] += y
            data["outputs"] += out

        data["raw_outputs"] += raw_outputs
        data["loss"] = running_loss
        return out, y, running_loss, data

    def get_metrics(self, data="eval", threshold=0.23, binary=True):
        """Compute metrics for a given stage."""
        if self.task == "class":
            threshold = self.classes_to_idx[threshold]

        if data == "eval":
            pred, y = self.valid_out, self.valid_y
        elif data == "train":
            pred, y = self.train_out, self.train_y
        elif data == "test":
            pred, y = self.test_out, self.test_y

        # --- Sanitize predictions and targets ---
        pred = np.array(pred, dtype=np.float32)
        y = np.array(y, dtype=np.float32)

        return Metrics(pred, y, self.target, threshold, binary)


    def compute_class_weights(self):
        """Compute normalized class weights from the training loader."""
        targets_all = []
        for batch in self.train_loader:
            targets_all.append(batch.y)
        targets_all = torch.cat(targets_all).squeeze().tolist()
        weights = torch.tensor(
            [targets_all.count(i) for i in self.classes], dtype=torch.float32
        )
        print(f"class counts {weights}")
        weights = 1.0 / weights
        weights = weights / weights.sum()
        print(f"class weights {weights}")
        return weights

    @staticmethod
    def print_epoch_data(stage, epoch, loss, acc, time_):
        """Short per-epoch print (legacy helper; kept for compatibility)."""
        acc_str = "None" if acc is None else "%1.4e" % acc
        print(f"epoch {epoch:04d} {stage} loss {loss:.3e} acc {acc_str} time {time_:.2e}s")

    def format_output(self, pred, target=None):
        """Format network outputs by task."""
        if self.task == "class":
            if target is not None:
                target = torch.tensor([self.classes_to_idx[int(x)] for x in target])
        elif self.transform_sigmoid is True:
            pred = torch.sigmoid(pred.reshape(-1))
        else:
            pred = pred.reshape(-1)
        return pred, target

    @staticmethod
    def update_name(hdf5, outdir):
        """Return a non-clashing HDF5 path in outdir."""
        fname = os.path.join(outdir, hdf5)
        count = 0
        hdf5_name = hdf5.split(".")[0]
        while os.path.exists(fname):
            count += 1
            hdf5 = f"{hdf5_name}_{count:03d}.hdf5"
            fname = os.path.join(outdir, hdf5)
        return fname

    def plot_loss(self, name=""):
        """Plot train/valid loss vs epoch and save to loss_epoch{name}.png."""
        nepoch = self.nepoch
        train_loss = self.train_loss
        valid_loss = self.valid_loss
        import matplotlib.pyplot as plt

        if len(valid_loss) > 1:
            plt.plot(range(1, nepoch + 1), valid_loss, c="red", label="valid")
        if len(train_loss) > 1:
            plt.plot(range(1, nepoch + 1), train_loss, c="blue", label="train")
            plt.title("Loss per epoch")
            plt.xlabel("Epoch")
            plt.ylabel("Loss")
            plt.legend()
            plt.savefig(f"loss_epoch{name}.png")
            plt.close()

    def plot_acc(self, name=""):
        """Plot train/valid accuracy vs epoch and save to acc_epoch{name}.png."""
        nepoch = self.nepoch
        train_acc = self.train_acc
        valid_acc = self.valid_acc
        import matplotlib.pyplot as plt

        if len(valid_acc) > 1:
            plt.plot(range(1, nepoch + 1), valid_acc, c="red", label="valid")
        if len(train_acc) > 1:
            plt.plot(range(1, nepoch + 1), train_acc, c="blue", label="train")
            plt.title("Accuracy per epoch")
            plt.xlabel("Epoch")
            plt.ylabel("Accuracy")
            plt.legend()
            plt.savefig(f"acc_epoch{name}.png")
            plt.close()

    def plot_hit_rate(self, data="eval", threshold=4, mode="percentage", name=""):
        """Plot hit-rate curve and save to hitrate{name}.png."""
        import matplotlib.pyplot as plt
        try:
            hitrate = self.get_metrics(data, threshold).hitrate()
            nb_models = len(hitrate)
            X = range(1, nb_models + 1)
            if mode == "percentage":
                hitrate /= hitrate.sum()
            plt.plot(X, hitrate, c="blue", label="train")
            plt.title("Hit rate")
            plt.xlabel("Rank")
            plt.ylabel("Hit rate")
            plt.legend()
            plt.savefig(f"hitrate{name}.png")
            plt.close()
        except Exception:
            print(f"no hit-rate plot for {self.task}")

    def plot_scatter(self):
        """Scatter plot of train/valid predictions vs truth."""
        import matplotlib.pyplot as plt
        net = self._net_for_eval()
        net.eval()

        pred, truth = {"train": [], "valid": []}, {"train": [], "valid": []}
        for data in self.train_loader:
            data = data.to(self.device)
            truth["train"] += data.y.tolist()
            pred["train"]  += net(data).reshape(-1).tolist()
        for data in self.valid_loader:
            data = data.to(self.device)
            truth["valid"] += data.y.tolist()
            pred["valid"]  += net(data).reshape(-1).tolist()

        plt.scatter(truth["train"], pred["train"], c="blue")
        plt.scatter(truth["valid"], pred["valid"], c="red")
        plt.show()

    def save_model(self, filename="model.pth.tar"):
        """Save model, optimizer, and (if present) EMA state."""
        state = {
            "model": self.model.state_dict(),
            "optimizer": self.optimizer.state_dict(),
            "node": self.node_feature,
            "edge": self.edge_feature,
            "target": self.target,
            "task": self.task,
            "classes": self.classes,
            "class_weight": self.class_weights,
            "batch_size": self.batch_size,
            "percent": self.percent,
            "lr": self.lr,
            "index": self.index,
            "shuffle": self.shuffle,
            "threshold": self.threshold,
            "cluster_nodes": self.cluster_nodes,
            "transform_sigmoid": self.transform_sigmoid,
            "ema_decay": getattr(self, "ema_decay", None),
            "ema_warmup_steps": getattr(self, "ema_warmup_steps", None),
            "ema_active": bool(getattr(self, "use_ema", False) and getattr(self, "ema", None) is not None),
            "ema_state": (self.ema.state_dict()
                          if getattr(self, "use_ema", False) and getattr(self, "ema", None) is not None
                          else None),
        }
        torch.save(state, filename)

    def load_params(self, filename):
        """Load hyper-params and states from a checkpoint (model/optimizer/EMA)."""
        self.device = torch.device(self.device_name)
        state = torch.load(filename, map_location=torch.device(self.device))

        self.node_feature = state["node"]
        self.edge_feature = state["edge"]
        self.target = state["target"]
        self.batch_size = state["batch_size"]
        self.percent = state["percent"]
        self.lr = state["lr"]
        self.index = state["index"]
        self.class_weights = state["class_weight"]
        self.task = state["task"]
        self.classes = state["classes"]
        self.threshold = state["threshold"]
        self.shuffle = state["shuffle"]
        self.cluster_nodes = state["cluster_nodes"]
        self.transform_sigmoid = state["transform_sigmoid"]

        self.opt_loaded_state_dict = state["optimizer"]
        self.model_load_state_dict = state["model"]

        self.ema_decay = state.get("ema_decay", getattr(self, "ema_decay", 0.999))
        self.ema_warmup_steps = state.get("ema_warmup_steps", getattr(self, "ema_warmup_steps", 0))
        self.use_ema = state.get("ema_active", getattr(self, "use_ema", False))
        self.ema_loaded_state_dict = state.get("ema_state", None)

    def _export_epoch_hdf5(self, epoch, data):
        """Write per-epoch data to HDF5 (train/valid/test groups) safely."""
        grp_name = f"epoch_{epoch:04d}"
        grp = self.f5.create_group(grp_name)

        # store general attributes
        grp.attrs["task"] = self.task
        grp.attrs["target"] = self.target
        grp.attrs["batch_size"] = self.batch_size

        for pass_type, pass_data in data.items():
            if pass_data is None or not isinstance(pass_data, dict):
                continue  # skip empty or invalid entries

            sg = grp.create_group(pass_type)

            for data_name, data_value in pass_data.items():
                if data_value is None:
                    continue  # skip empty data

                # handle string scalars
                if isinstance(data_value, str):
                    string_dt = h5py.string_dtype(encoding='utf-8')
                    sg.create_dataset(data_name, data=np.array(data_value, dtype=object), dtype=string_dt)
                    continue

                # convert everything else to numpy array
                try:
                    arr = np.array(data_value)
                except Exception as e:
                    print(f"Skipping {data_name} in {pass_type}: cannot convert to array ({e})")
                    continue

                if arr.size == 0:
                    continue  # skip empty arrays

                # handle string arrays
                if arr.dtype.kind in ('U', 'S', 'O'):  # unicode, bytes, or object
                    arr_list = arr.tolist()  # convert to Python list of strings
                    string_dt = h5py.string_dtype(encoding='utf-8')
                    sg.create_dataset(data_name, data=arr_list, dtype=string_dt)
                else:
                    # numeric arrays
                    sg.create_dataset(data_name, data=arr)



