#!/usr/bin/env python3
"""
DeepRank-Ab Inference Pipeline
(robust against Biopython chain/residue ID collisions)
"""

import argparse
import hashlib
import json
import logging
import os
import shutil
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

import h5py
import pandas as pd
import requests
import torch
from Bio.PDB import PDBParser, PDBIO
from esm import FastaBatchedDataset, pretrained
import numpy as np
from anarci import anarci

# Add project root
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

BASE = Path(__file__).resolve().parents[1]  # DeepRank-Ab/
hmmscan_path = BASE / "src" / "tools" / "ANARCI" 


# Local modules
from src.DataSet import HDF5DataSet, PreCluster
from src.GraphGenMP import GraphHDF5
from src.tools.annotate import annotate_folder_one_by_one_mp
from src.NeuralNet_focal_EMA import NeuralNet
from src.EGNN import egnn


# ===============================================================
# CONFIGURATION
# ===============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent

MODEL_PATH = (
    ROOT_DIR
    / "src"
    / "weights"
    / "af"
    / "top_mse_ep37_mse2.3192e-02_treg_ydockq_b128_e100_lr2e-03.pth.tar"
)
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

TOKS_PER_BATCH = 4096
REPR_LAYERS = [33]
TRUNCATION_SEQ_LENGTH = 2500
INCLUDE = ["mean", "per_tok"]

MAX_cores = 50
BATCH_SIZE = 64
NUM_WORKERS = 96

TARGET = "dockq"
TASK = "reg"
THRESHOLD = 0.23

NODE_FEATURES = ["atom_type", "polarity", "bsa", "region", "embedding"]
EDGE_FEATURES = ["voro_area", "covalent", "vdw", "orientation"]

ESM_MODEL = "esm2_t33_650M_UR50D"
EXPECTED_CHECKSUMS = [
    "ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c",
    "8ffe6edbd4173dc8d45c2cd5cb27d43aad77ec26b4c768200c58ae1f96693575",
]

# Logging
log = logging.getLogger("DeepRankAB")
log.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(" [%(levelname)s] %(message)s"))
log.addHandler(handler)


# ===============================================================
# HELPERS
# ===============================================================

def _norm_chain(x: Optional[str]) -> Optional[str]:
    """Normalize chain ID args. Treat '-', 'none', 'null' etc. as missing."""
    if x is None:
        return None
    s = str(x).strip()
    if s.lower() in {"", "none", "null", "-", "na"}:
        return None
    return s


def clamp_cores(requested: int, task_count: int) -> int:
    """
    Clamp number of cores to avoid creating more processes than tasks.
    Prevents nproc from exceeding CPU count or task count.
    """
    avail = os.cpu_count() or 1
    return max(1, min(requested, avail, task_count))


# ===============================================================
# WORKSPACE SETUP
# ===============================================================

def setup_workspace(identificator: str) -> Path:
    workspace = Path.cwd() / identificator
    workspace.mkdir(parents=True, exist_ok=True)
    log.info(f"Workspace → {workspace}")
    return workspace


# ===============================================================
# CHAIN AUTO-DETECTION
# ===============================================================
def detect_chains_anarci(pdb_file: Path) -> Tuple[str, Optional[str], str]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)
    model = next(structure.get_models())
    mapping = three_to_one()

    heavy_chains: List[str] = []
    light_chains: List[str] = []
    antigen_chains: List[str] = []

    for chain in model.get_chains():
        seq = "".join(mapping.get(r.get_resname(), "X") for r in chain.get_residues())

        if len(seq) < 10:
            antigen_chains.append(chain.id)
            continue

        results, _, details = anarci(
            [(chain.id, seq)],
            scheme="imgt",
            assign_germline=True,
            hmmerpath=hmmscan_path,
        )

        hit = results[0]      # hits for this sequence
        det = details[0]      # detail dicts for this sequence


        if hit is None or all(h is None for h in hit):
            antigen_chains.append(chain.id)
        else:
            # details[0] is a list of dicts, one per domain hit
            chain_type = det[1][0][-1]
            if chain_type == "H":
                heavy_chains.append(chain.id)
            elif chain_type in ("L", "K"):
                light_chains.append(chain.id)
            else:
                antigen_chains.append(chain.id)

        log.info(f"Chain {chain.id}: ANARCI chain_type = {chain_type} if hit and det else 'no hit'")


    log.info(
        f"ANARCI detected → heavy={heavy_chains}, light={light_chains}, "
        f"antigen={antigen_chains}"
    )

    if not heavy_chains:
        raise ValueError(
            f"No heavy chain detected by ANARCI in {pdb_file.name}. "
            f"Chains present: {[c.id for c in model.get_chains()]}. "
            f"Pass --heavy_chain_id manually if needed."
        )
    if not antigen_chains:
        raise ValueError(
            f"No antigen chain detected in {pdb_file.name} (all chains matched as "
            f"antibody). Pass --antigen_chain_id manually if needed."
        )

    if len(heavy_chains) > 1:
        log.warning(
            f"Multiple heavy chains detected {heavy_chains}, using '{heavy_chains[0]}'"
        )
    if len(light_chains) > 1:
        log.warning(
            f"Multiple light chains detected {light_chains}, using '{light_chains[0]}'"
        )
    if len(antigen_chains) > 1:
        log.warning(
            f"Multiple antigen chains detected {antigen_chains}, "
            f"using '{antigen_chains[0]}'"
        )

    return (
        heavy_chains[0],
        light_chains[0] if light_chains else None,
        antigen_chains[0],
    )


# ===============================================================
# PDB PROCESSING
# ===============================================================

def split_input_pdb(pdb_file: Path, out_dir: Path) -> List[Path]:
    """
    Split an ensemble PDB into separate model PDB files
    while preserving REMARK model names.

    Example REMARK:
        REMARK  4      MODEL     0 FROM Target319/055a74d834b0b51a855178676f0f16f2.pdb

    Output:
        055a74d834b0b51a855178676f0f16f2.pdb
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    out_dir.mkdir(exist_ok=True)

    # map: model index -> original model name
    model_remarks = {}
    with open(pdb_file, "r") as f:
        for line in f:
            if "REMARK" in line and "MODEL" in line:
                parts = line.split()
                try:
                    ensemble_model_id = int(parts[3])
                    original_name = parts[-1].split("/")[-1]
                    original_name = original_name.replace(".pdb", "")
                    model_remarks[ensemble_model_id] = {
                        "remark": line.strip(),
                        "name": original_name,
                    }
                except Exception as e:
                    print(f"[WARN] failed parsing line:\n{line}\n{e}")

    io = PDBIO()
    saved = []

    for model in structure:
        model_id = model.id
        if model_id in model_remarks:
            model_name = model_remarks[model_id]["name"]
        else:
            model_name = f"model_{model_id}"
        out = out_dir / f"{model_name}.pdb"
        io.set_structure(model)
        io.save(str(out))

        if model_id in model_remarks:
            with open(out, "r") as f:
                pdb_text = f.read()
            with open(out, "w") as f:
                f.write(model_remarks[model_id]["remark"] + "\n")
                f.write(pdb_text)
        saved.append(out)

    print(f"Split PDB into {len(saved)} model(s)")
    return saved


def three_to_one() -> dict:
    return {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }


def get_chain_sequence(structure, chain_id: Optional[str]) -> str:
    if not chain_id:
        return ""
    try:
        chain = structure[0][chain_id]
    except KeyError:
        return ""

    mapping = three_to_one()
    seq = []
    for res in chain:
        seq.append(mapping.get(res.get_resname(), "X"))
    return "".join(seq)


def create_merged_pdb(
    pdb_file, heavy_chain_id, light_chain_id, antigen_chain_id, output_pdb
):
    """
    Collision-proof merged PDB builder.

    Output contains:
      - antibody (heavy + optional light) in chain 'A'
      - antigen in chain 'B'
    """
    parser = PDBParser(QUIET=True)
    s0 = parser.get_structure("model", pdb_file)
    m0 = next(s0.get_models())

    HID = _norm_chain(heavy_chain_id)
    LID = _norm_chain(light_chain_id)
    AID = _norm_chain(antigen_chain_id)

    chain_ids = {c.id for c in m0.get_chains()}
    if HID not in chain_ids:
        raise ValueError(f"Heavy chain '{HID}' not found (have {sorted(chain_ids)})")
    if AID not in chain_ids:
        raise ValueError(f"Antigen chain '{AID}' not found (have {sorted(chain_ids)})")
    if LID is not None and LID not in chain_ids:
        raise ValueError(f"Light chain '{LID}' not found (have {sorted(chain_ids)})")

    chH = m0[HID]
    chL = m0[LID] if LID is not None else None
    chAg = m0[AID]

    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain

    s = Structure("merged")
    m = Model(0)
    s.add(m)

    chAb = Chain("A")
    chBg = Chain("B")
    m.add(chAb)
    m.add(chBg)

    def _sorted_res(chain_obj):
        residues = list(chain_obj.get_residues())

        def key(r):
            hetflag, resseq, icode = r.id
            icode_key = "" if icode == " " else icode
            return (hetflag.strip(), int(resseq), icode_key)

        return sorted(residues, key=key)

    def _copy_into(dst_chain, src_chain, start_resseq=1):
        cur = start_resseq
        for r0 in _sorted_res(src_chain):
            r = r0.copy()
            hetflag, _, _icode = r.id
            r.id = (hetflag, cur, " ")
            dst_chain.add(r)
            cur += 1
        return cur

    next_res = 1
    next_res = _copy_into(chAb, chH, start_resseq=next_res)
    if chL is not None:
        next_res = _copy_into(chAb, chL, start_resseq=next_res)

    _copy_into(chBg, chAg, start_resseq=1)

    output_pdb = Path(output_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(s)
    io.save(str(output_pdb))


def convert_pdb_to_fastas(
    pdb_file: Path,
    fasta_outdir: Path,
    heavy_chain_id="H",
    light_chain_id="L",
    antigen_chain_id="A",
):
    """
    Create:
      - annotation FASTA: separate H and L (if present)
      - ESM FASTA: merged antibody sequence as A, antigen as B
    """
    fasta_outdir.mkdir(parents=True, exist_ok=True)

    HID = _norm_chain(heavy_chain_id)
    LID = _norm_chain(light_chain_id)
    AID = _norm_chain(antigen_chain_id)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)

    root = pdb_file.stem

    seq_H = get_chain_sequence(structure, HID)
    seq_L = get_chain_sequence(structure, LID)
    seq_Ag = get_chain_sequence(structure, AID)

    fasta_annot = fasta_outdir / f"{root}_HL.fasta"
    with open(fasta_annot, "w") as f:
        if seq_H:
            f.write(f">{root}.H\n{seq_H}\n")
        if seq_L:
            f.write(f">{root}.L\n{seq_L}\n")

    fasta_esm = fasta_outdir / f"{root}_merged_A_B.fasta"
    with open(fasta_esm, "w") as f:
        f.write(f">{root}.A\n{seq_H + seq_L}\n")
        f.write(f">{root}.B\n{seq_Ag}\n")

    log.info(f"FASTA files created for {pdb_file.name}")
    return fasta_annot, fasta_esm


def preprocess_input_pdb(
    work_dir: Path,
    pdb_file: Path,
    heavy_chain_id: str,
    light_chain_id: str,
    antigen_chain_id: str,
):
    fasta_out = work_dir / "fastas"
    fasta_annot, fasta_esm = convert_pdb_to_fastas(
        pdb_file,
        fasta_out,
        heavy_chain_id=heavy_chain_id,
        light_chain_id=light_chain_id,
        antigen_chain_id=antigen_chain_id,
    )

    merged_pdb = work_dir / "processed" / f"{pdb_file.stem}.pdb"
    create_merged_pdb(
        pdb_file, heavy_chain_id, light_chain_id, antigen_chain_id, merged_pdb
    )

    return merged_pdb, (fasta_annot, fasta_esm)


# ===============================================================
# ESM EMBEDDINGS
# ===============================================================

def calculate_checksum(path: str, algo="sha256") -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def download_weights(url: str, dest: str) -> str:
    resp = requests.get(url, stream=True, timeout=10)
    resp.raise_for_status()
    with open(dest, "wb") as f:
        for chunk in resp.iter_content(8192):
            f.write(chunk)
    return dest


def fetch_weights() -> str:
    models = [
        (
            "WEIGHT_PATH",
            f"{ESM_MODEL}.pt",
            f"https://dl.fbaipublicfiles.com/fair-esm/models/{ESM_MODEL}.pt",
            EXPECTED_CHECKSUMS[0],
        ),
        (
            "REG_WEIGHT_PATH",
            f"{ESM_MODEL}-contact-regression.pt",
            f"https://dl.fbaipublicfiles.com/fair-esm/regression/{ESM_MODEL}-contact-regression.pt",
            EXPECTED_CHECKSUMS[1],
        ),
    ]

    for env_var, fname, url, checksum in models:
        path = os.getenv(env_var) or fname
        if not os.path.exists(path):
            download_weights(url, path)

        if calculate_checksum(path) != checksum:
            log.warning(f"Checksum mismatch for {path}, re-downloading.")
            download_weights(url, path)

    return f"{ESM_MODEL}.pt"


def get_model_output(toks, model, layers):
    out = model(toks, repr_layers=layers, return_contacts="contacts" in INCLUDE)
    return {layer: t.cpu() for layer, t in out["representations"].items()}


def process_batch(labels, strs, outdir, representations):
    paths = []
    for i, label in enumerate(labels):
        data = {}
        trunc = min(TRUNCATION_SEQ_LENGTH, len(strs[i]))

        if "per_tok" in INCLUDE:
            data["representations"] = {
                l: t[i, 1 : trunc + 1].clone() for l, t in representations.items()
            }

        if "mean" in INCLUDE:
            data["mean_representations"] = {
                l: t[i, 1 : trunc + 1].mean(0).clone()
                for l, t in representations.items()
            }

        outpath = outdir / f"{label}.pt"
        torch.save(data, outpath)
        paths.append(outpath)

    return paths


def get_embedding(fasta_file: Path, output_dir: Path) -> List[Path]:
    log.info("Generating ESM embeddings...")

    esm_path = fetch_weights()
    model, alphabet = pretrained.load_model_and_alphabet(esm_path)
    model.eval()
    if torch.cuda.is_available():
        model.cuda()

    dataset = FastaBatchedDataset.from_file(fasta_file)
    batches = dataset.get_batch_indices(TOKS_PER_BATCH, extra_toks_per_seq=1)
    loader = torch.utils.data.DataLoader(
        dataset,
        collate_fn=alphabet.get_batch_converter(TRUNCATION_SEQ_LENGTH),
        batch_sampler=batches,
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    repr_layers = [
        (i + model.num_layers + 1) % (model.num_layers + 1) for i in REPR_LAYERS
    ]

    emb_paths: List[Path] = []
    with torch.no_grad():
        for labels, strs, toks in loader:
            if torch.cuda.is_available():
                toks = toks.cuda(non_blocking=True)

            reps = get_model_output(toks, model, repr_layers)
            emb_paths.extend(process_batch(labels, strs, output_dir, reps))

    return emb_paths


# ===============================================================
# CLUSTERING AND GRAPH GENERATION
# ===============================================================

def cluster(hdf5_path: str):
    dataset = HDF5DataSet(name="Train", root="./", database=hdf5_path)
    PreCluster(dataset, method="mcl")
    log.info("Clustering completed.")


def gen_graph_cdrs_orientation_contacts_one_by_one(
    pdb_dir: str,
    n_cores: int,
    outfile_name: str,
    use_regions=True,
    region_json=None,
    antigen_chainid=None,
    graph_type="atom",
    use_voro=False,
    embedding_path=None,
    contact_features=True,
):
    base = Path(pdb_dir).name
    tmpdir = tempfile.mkdtemp(prefix=f"{base}_")

    GraphHDF5(
        pdb_path=pdb_dir,
        graph_type=graph_type,
        outfile=outfile_name,
        nproc=n_cores,
        tmpdir=tmpdir,
        use_regions=use_regions,
        region_json=region_json,
        add_orientation=True,
        contact_features=contact_features,
        antigen_chainid=antigen_chainid,
        use_voro=use_voro,
        embedding_path=embedding_path,
    )
    return outfile_name


def correct_json(work_dir: Path) -> Path:
    region_json = Path(work_dir) / "annotations" / "annotations_cdrs.json"
    if not region_json.is_file():
        raise FileNotFoundError(f"Region JSON not found: {region_json}")

    with open(region_json) as f:
        data = json.load(f)

    new = {k.replace(".pdb", ""): v for k, v in data.items()}
    with open(region_json, "w") as f:
        json.dump(new, f, indent=4)

    return region_json


def add_embedding(work_dir: Path, hdf5_path: str):
    """
    Add ESM embedding feature to each molecule in the graph HDF5.
    """
    base = work_dir / "processed" / "embeddings"

    def _resolve_pt(mol_name: str, chain: str) -> Path:
        candidates = [mol_name, mol_name.replace(".pdb", ""), mol_name.split("/")[-1]]
        if candidates[0].endswith("_graph"):
            candidates.append(candidates[0].removesuffix("_graph"))
        for cand in candidates:
            pt = base / f"{cand}.{chain}.pt"
            if pt.is_file():
                return pt
        pt = base / f"{Path(mol_name).stem}.{chain}.pt"
        return pt

    with h5py.File(hdf5_path, "a") as f:
        for mol in list(f.keys()):
            residues = f[mol]["nodes"][()]
            emb = torch.zeros(len(residues), 1)

            for i, res in enumerate(residues):
                chain, idx = res[0].decode(), int(res[1].decode())
                pt = _resolve_pt(mol, chain)

                if not pt.is_file():
                    raise FileNotFoundError(
                        f"Missing embedding: {pt} (mol={mol}, chain={chain})"
                    )

                data = torch.load(pt, map_location="cpu")
                vecs = data["representations"][33]
                if idx < 1 or idx > vecs.shape[0]:
                    raise IndexError(
                        f"Residue index out of bounds for {pt}: idx={idx}, "
                        f"len={vecs.shape[0]}"
                    )
                emb[i] = vecs[idx - 1].mean()

            grp = f[mol].require_group("node_data")
            if "embedding" in grp:
                del grp["embedding"]
            grp.create_dataset("embedding", data=emb.numpy())

    log.info(f"Added embeddings → {hdf5_path}")


# ===============================================================
# MODEL EVALUATION
# ===============================================================

def deeprank_evaluate_model(
    target_name: str,
    hdf5_test: str,
    model_path: str = str(MODEL_PATH),
    save_name: str = "eval_predictions.hdf5",
) -> Path:
    net = NeuralNet(
        database=hdf5_test,
        Net=egnn,
        node_feature=NODE_FEATURES,
        edge_feature=EDGE_FEATURES,
        target=None,
        task=TASK,
        batch_size=BATCH_SIZE,
        num_workers=NUM_WORKERS,
        device_name=DEVICE,
        shuffle=False,
        pretrained_model=model_path,
        cluster_nodes="mcl",
    )

    out = Path(hdf5_test).parent / save_name
    net.predict(database_test=hdf5_test, hdf5=str(out))
    log.info(f"Predictions saved → {out}")

    return out


def hdf5_to_csv(hdf5_path: str) -> str:
    """Convert DeepRank-Ab prediction HDF5 to CSV with mol and dockq values."""
    hdf5_path = Path(hdf5_path)
    out_csv = hdf5_path.with_suffix(".csv")

    with h5py.File(hdf5_path, "r") as f:
        group = f["epoch_0000"]["pred"]
        mol = group["mol"][()]
        dockq = group["outputs"][()]
        mol = [m.decode("utf-8") for m in mol]

    df = pd.DataFrame({"mol": mol, "dockq": dockq})
    df.to_csv(out_csv, index=False)
    return str(out_csv)


def ca_coords(pdb_file: Path, chain_id: str, max_resid=125) -> np.ndarray:
    coords = []
    residue_index = 0

    with open(pdb_file, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()
            chain = line[21].strip()

            if chain != chain_id or atom_name != "CA":
                continue

            if residue_index > max_resid:
                break

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            coords.append([x, y, z])
            residue_index += 1

    coords_arr = np.array(coords, dtype=float)

    if coords_arr.size == 0:
        coords_arr = coords_arr.reshape((0, 3))

    return coords_arr


# ===============================================================
# LOW-CONTACT FLAG FOR ANTIBODIES
# ===============================================================

def flag_low_contacts(
    pdb_file: Path, heavy_chain_id: str, light_chain_id: str, threshold: int = 15
):
    """
    Flag ONLY if:
      - both heavy (H) and light (L) chains are present
      - AND they have low contacts

    If one is missing → do NOT flag (valid case)
    """
    coords_H = ca_coords(pdb_file, heavy_chain_id, max_resid=125)
    coords_L = ca_coords(pdb_file, light_chain_id, max_resid=115)

    if coords_H.shape[0] == 0 or coords_L.shape[0] == 0:
        return False

    dists = np.linalg.norm(coords_H[:, None, :] - coords_L[None, :, :], axis=-1)
    n_contacts = int(np.sum(dists < 8.0))

    if n_contacts < threshold:
        log.warning(
            f"[FLAG] LOW H-L CONTACTS → {pdb_file.name} with {n_contacts} contacts"
        )
        return True

    return False


def parse_output(
    csv_output: str, original_pdb_dir: Path, heavy_chain_id=None, light_chain_id=None
) -> None:
    """
    Add quality flags based on ORIGINAL PDBs (before merging).
    Skips check if either H or L chain is missing (valid nanobody or single-chain antibody),
    and marks them as 'not_applicable'.
    """
    df = pd.read_csv(csv_output)
    df = df.rename(columns={"mol": "pdb_id", "dockq": "predicted_dockq"})

    flags = []

    for pdb_id in df["pdb_id"]:
        matches = list(original_pdb_dir.glob(f"{pdb_id}*.pdb"))
        pdb_file = matches[0]

        if not heavy_chain_id or not light_chain_id:
            flags.append("not_applicable")
            continue

        is_bad = flag_low_contacts(pdb_file, heavy_chain_id, light_chain_id)
        flags.append("low_HL_contacts" if is_bad else "ok")

    df["quality_flag"] = flags
    df = df.sort_values(by="predicted_dockq", ascending=False)
    df.to_csv(csv_output, index=False)

    log.info(f"Output written → {csv_output}")


# ===============================================================
# MAIN PIPELINE
# ===============================================================

def main():
    parser = argparse.ArgumentParser(
        description="DeepRank-Ab inference pipeline. Chain IDs are auto-detected via "
                    "ANARCI when not provided."
    )
    parser.add_argument("pdb_file", help="Input PDB file (single model or ensemble)")
    parser.add_argument(
        "--heavy_chain_id",
        default=None,
        help="Heavy chain ID (auto-detected via ANARCI if omitted)",
    )
    parser.add_argument(
        "--light_chain_id",
        default=None,
        help="Light chain ID (auto-detected via ANARCI if omitted; "
             "pass '-' to explicitly indicate a nanobody/VHH)",
    )
    parser.add_argument(
        "--antigen_chain_id",
        default=None,
        help="Antigen chain ID (auto-detected via ANARCI if omitted)",
    )
    args = parser.parse_args()

    pdb_file = Path(args.pdb_file)

    # ------------------------------------------------------------------
    # Resolve chain IDs: manual override takes priority, ANARCI as fallback
    # ------------------------------------------------------------------
    manual_heavy = _norm_chain(args.heavy_chain_id)
    manual_light = _norm_chain(args.light_chain_id)
    manual_antigen = _norm_chain(args.antigen_chain_id)

    if manual_heavy and manual_antigen:
        # Both essential chains provided manually — use as-is
        HID = manual_heavy
        LID = manual_light  # may be None (nanobody)
        AID = manual_antigen
        log.info(f"Using manually specified chains → H={HID}, L={LID}, Ag={AID}")
    else:
        # Auto-detect with ANARCI (runs on first model only; chain IDs are
        # consistent across all models in a standard ensemble)
        log.info("Chain IDs not fully specified — running ANARCI auto-detection...")
        HID, LID, AID = detect_chains_anarci(pdb_file)
        # Allow a manual override to win over auto-detection for any single chain
        if manual_heavy:
            log.info(f"Overriding auto-detected heavy chain with --heavy_chain_id={manual_heavy}")
            HID = manual_heavy
        if args.light_chain_id is not None:  # explicit '-' means "no light chain"
            log.info(f"Overriding auto-detected light chain with --light_chain_id={manual_light}")
            LID = manual_light
        if manual_antigen:
            log.info(f"Overriding auto-detected antigen chain with --antigen_chain_id={manual_antigen}")
            AID = manual_antigen
        log.info(f"Final chains → H={HID}, L={LID}, Ag={AID}")

    # Workspace name encodes chain assignment for traceability
    lid_tag = LID if LID else "nb"  # 'nb' = nanobody / no light chain
    identificator = f"{pdb_file.stem}-deeprank_ab_pred_{HID}{lid_tag}_{AID}"
    work = setup_workspace(identificator)

    # Copy PDB
    copied = work / pdb_file.name
    shutil.copy(pdb_file, copied)

    split_dir = work / "original_models"
    split_dir.mkdir(exist_ok=True)

    pdb_models = split_input_pdb(copied, split_dir)

    avail = os.cpu_count() or 1
    cores = max(1, min(MAX_cores, avail))

    processed_dir = work / "processed"
    processed_dir.mkdir(exist_ok=True)

    # ============================================================
    # Build merged PDBs + embeddings
    # ============================================================
    fasta_annot = None
    embed_dir = processed_dir / "embeddings"
    embed_dir.mkdir(exist_ok=True)

    for pdb in pdb_models:
        _merged_pdb, (fasta_annot, fasta_esm) = preprocess_input_pdb(
            work, pdb, HID, LID, AID
        )
        get_embedding(fasta_esm, embed_dir)

    if fasta_annot is None:
        raise RuntimeError("No FASTA files were produced (unexpected).")

    # ============================================================
    # Annotation
    # ============================================================
    anno_dir = work / "annotations"
    anno_dir.mkdir(exist_ok=True)

    n_cores_anno = clamp_cores(cores, len(pdb_models))

    annotate_folder_one_by_one_mp(
        processed_dir,
        Path(fasta_annot.parent),
        output_dir=str(anno_dir),
        n_cores=n_cores_anno,
        antigen_chainid="B",
    )

    region_json = correct_json(work)
    log.info("Region JSON corrected.")

    # ============================================================
    # Graph generation
    # ============================================================
    pdbs_for_graph = sorted(processed_dir.glob("*.pdb"))
    n_cores_graph = clamp_cores(cores, len(pdbs_for_graph))
    log.info(
        f"Using {n_cores_graph} cores for graph generation of "
        f"{len(pdbs_for_graph)} models."
    )

    graph_out = work / f"{identificator}_graph.hdf5"

    gen_graph_cdrs_orientation_contacts_one_by_one(
        pdb_dir=str(processed_dir),
        n_cores=n_cores_graph,
        outfile_name=str(graph_out),
        use_regions=True,
        region_json=str(region_json),
        antigen_chainid="B",
        graph_type="atom",
        use_voro=True,
        embedding_path=str(embed_dir),
        contact_features=True,
    )

    # ============================================================
    # Evaluation
    # ============================================================
    cluster(str(graph_out))

    deeprank_evaluate_model(
        "dockq",
        str(graph_out),
        model_path=str(MODEL_PATH),
        save_name=f"{identificator}_predictions.hdf5",
    )

    pred_h5 = work / f"{identificator}_predictions.hdf5"
    csv_out = hdf5_to_csv(str(pred_h5))

    # ============================================================
    # Contact check
    # ============================================================
    parse_output(csv_out, split_dir, heavy_chain_id=HID, light_chain_id=LID)


if __name__ == "__main__":
    main()