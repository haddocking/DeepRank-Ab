import os
import sys
import glob
import h5py
from tqdm import tqdm
import multiprocessing as mp
from functools import partial
import pickle
import warnings
import json
from Bio.PDB.PDBParser import PDBParser
from Bio import BiopythonWarning


warnings.filterwarnings("ignore", category=BiopythonWarning)

import sys, os
sys.path.append(os.path.dirname(__file__))


"""
Build residue-level graphs from many docking-model PDB files (optionally a reference/native),
and store them in a single HDF5. Optionally add per-residue ESM-2 embeddings.
"""


class GraphHDF5(object):
    def __init__(
        self,
        pdb_path,
        ref_path=None,
        graph_type="residue",
        embedding_path=None,   # optional
        select=None,
        outfile="graph.hdf5",
        nproc=1,
        use_tqdm=True,
        tmpdir="./",
        limit=None,
        biopython=False,
        use_regions=True,
        region_json=None,
        add_orientation=False,
        contact_features: bool = True,
        antigen_chainid: str = 'A',
        use_voro: bool = False,
    ):
        """Entry point that orchestrates graph construction and output."""
        # --- gather PDB files ---
        pdbs = [f for f in os.listdir(pdb_path) if f.endswith(".pdb")]
        for i in pdbs:
            if '_xray_tidy' in i:
                pdbs.remove(i)

        if select:
            pdbs = [f for f in pdbs if f.startswith(select)]
        pdbs = [os.path.join(pdb_path, f) for f in pdbs]
        if limit is not None:
            pdbs = pdbs[limit[0]:limit[1]] if isinstance(limit, list) else pdbs[:limit]

        if len(pdbs) == 0:
            raise FileNotFoundError(f"no .pdb files found under {pdb_path}")


        # --- reference structure ---
        ref = os.path.join(ref_path, base + ".pdb") if ref_path else None

        # --- embeddings are optional ---
        self.embedding_path = embedding_path
        if self.embedding_path is None:
            print("no embeddings path; skipping", flush=True)

        self.antigen_chainid = antigen_chainid 
        self.use_regions = use_regions
        self.graph_type = graph_type
        self.use_voro = use_voro

        # --- region mapping (optional) ---
        if self.use_regions:
            if not region_json or not os.path.isfile(region_json):
                raise FileNotFoundError(f"region_json not found: {region_json}")
            with open(region_json) as f:
                raw = json.load(f)  # model_name -> list of [AA, tag] #no .pdb in ilrais's aannottaion 
           
            # map (model, residue_number) -> region tag using chain A numbering

            parser = PDBParser(QUIET=True)
            self.region_map = {}
            for pdbfile in pdbs:
                model_name = os.path.basename(pdbfile).replace('.pdb', '')
                if model_name not in raw:
                    continue

                annot = raw[model_name]  # e.g. [("A", "FR"), ("S", "L1"), …]
                struct = parser.get_structure(model_name, pdbfile)
                chainA = struct[0][next(c.get_id() for c in struct[0] if c.get_id() != antigen_chainid)]
                pdb_nums = [r.get_id()[1] for r in chainA.get_residues()]

                # for resnum, (_, tag) in zip_longest(pdb_nums, annot, fillvalue=('-', 'UNK')):
                #     self.region_map[(model_name, resnum)] = tag
                num_annot = len(annot)

                for i, resnum in enumerate(pdb_nums):
                    if i < num_annot:
                        _, tag = annot[i]   # use the annotation
                    else:
                        tag = "UNK"         # fill missing residues with UNK
                    self.region_map[(model_name, resnum)] = tag
                    #print(f"mapping {model_name} res {resnum} to region {tag}")



        # feature flags
        self.add_orientation = add_orientation
        self.contact_features = contact_features

        # --- build graphs ---
        if nproc == 1:
            self.get_all_graphs(
                pdbs, ref, outfile,
                use_tqdm=use_tqdm,
                biopython=biopython,
                region_map=self.region_map,
                antigen_chainid=self.antigen_chainid
            )
        else:
            if not os.path.isdir(tmpdir):
                os.mkdir(tmpdir)

            chunksize = max(1, len(pdbs) // (nproc * 8))
            ctx = mp.get_context("spawn")
            
            print(f"running in parallel with {nproc} processes")
            
            pool = ctx.Pool(nproc)
            part = partial(
                self._pickle_one_graph,
                ref=ref,
                tmpdir=tmpdir,
                biopython=biopython,
                region_map=self.region_map,
                antigen_chainid=self.antigen_chainid,
                embedding_path=self.embedding_path 
            )
            
            # FIXED: Replace blocking pool.map with streaming imap_unordered + safe cleanup
            try:
                for _ in pool.imap_unordered(part, pdbs, chunksize=chunksize):
                    pass
            finally:
                pool.close()
                pool.join()

            # merge pickled graphs into the HDF5 file
            with h5py.File(outfile, "w") as f5:
                for fn in glob.glob(os.path.join(tmpdir, "*.pkl")):
                    with open(fn, "rb") as f:
                        g = pickle.load(f)
                    try:
                        g.nx2h5(f5)
                    except Exception as e:
                        print(f"error storing graph {fn}: {e}")
                    os.remove(fn)

        # --- clean up temporary files ---
        for ext in ("*.izone", "*.lzone", "*.refpairs"):
            for fn in glob.glob(ext):
                try:
                    os.remove(fn)
                except OSError:
                    pass

        # --- append embeddings only when provided ---
        if self.embedding_path:
            try:
                self._add_embedding(outfile=outfile, pdbs=pdbs, embedding_path=self.embedding_path)
            except Exception as e:
                print(f"warning: embedding step failed; graphs are ready: {e}", flush=True)


    def get_all_graphs(
            self,
            pdbs,
            ref,
            outfile,
            use_tqdm=True,
            biopython=False,
            region_map=None,
            antigen_chainid='A'
        ):
        """Generate all graphs in memory and save to HDF5."""
        graphs = []

        if use_tqdm:
            desc = "{:25s}".format("   create HDF5")
            lst = tqdm(pdbs, desc=desc, file=sys.stdout)
        else:
            lst = pdbs

        for name in lst:
            try:
                g = self._get_one_graph(
                    name,
                    ref,
                    biopython,
                    region_map=self.region_map if self.use_regions else None,
                    antigen_chainid=self.antigen_chainid,
                    embedding_path=self.embedding_path,
                )
                ### FIX: skip failed graphs
                if g is None:
                    print(f"[WARN] Skipped graph for {name}")
                    continue

                graphs.append(g)

            except Exception as e:
                print(f"error computing graph {name}: {e}")

        with h5py.File(outfile, "w") as f5:
            for g in graphs:
                ### FIX: skip None
                if g is None:
                    continue
                try:
                    g.nx2h5(f5)
                except Exception as e:
                    print(f"error writing graph {getattr(g,'pdb','UNKNOWN')}: {e}")


    def _pickle_one_graph(self, name, ref, tmpdir="./", biopython=False, region_map=None, antigen_chainid=None, embedding_path=None):
        """Build a graph and pickle it."""
        chain_id = antigen_chainid if antigen_chainid is not None else self.antigen_chainid

        try:
            if self.graph_type == 'residue':
                from ResidueGraph import ResidueGraph
                g = ResidueGraph(
                    pdb=name,
                    biopython=biopython,
                    region_map=region_map,
                    add_orientation=self.add_orientation,
                    contact_features=self.contact_features,
                    antigen_chainid=chain_id
                )
            elif self.graph_type == 'atom':
                from AtomGraph import AtomGraph
                g = AtomGraph(
                    pdb=name,
                    antigen_chainid=chain_id,
                    region_map=region_map,
                    use_voro=self.use_voro,
                    contact_features=self.contact_features,
                )

        except Exception as e:
            print(f"[WARN] Failed graph {name}: {e}")
            return  ### FIX: don't write bad graphs

        mol_name = os.path.basename(name).rsplit('.', 1)[0]
        fname = os.path.join(tmpdir, mol_name + ".pkl")

        with open(fname, "wb") as f:
            pickle.dump(g, f)


    def _get_one_graph(self, name, ref, biopython=False, region_map=None, antigen_chainid='A', embedding_path=None):
        """Build one graph and return it."""
        try:
            if self.graph_type == 'residue':
                from ResidueGraph import ResidueGraph
                g = ResidueGraph(
                    pdb=name,
                    biopython=biopython,
                    region_map=self.region_map,
                    add_orientation=self.add_orientation,
                    contact_features=self.contact_features,
                    antigen_chainid=antigen_chainid
                )

            elif self.graph_type == 'atom':
                from AtomGraph import AtomGraph
                g = AtomGraph(
                    pdb=name,
                    antigen_chainid=antigen_chainid,
                    region_map=region_map,
                    use_voro=self.use_voro,
                )
            return g

        except Exception as e:
            print(f"[WARN] Failed graph {name}: {e}")
            return None   ### FIX: return None so caller can skip 

    @staticmethod
    def _add_embedding(outfile, pdbs, embedding_path):
        """Append per-residue embeddings to the HDF5 if available."""
        # import torch here to avoid a hard dependency when embeddings are skipped
        try:
            import torch
        except Exception as e:
            raise RuntimeError(f"PyTorch is required to add embeddings: {e}")

        with h5py.File(outfile, 'r+') as f:
            mols = list(f.keys())
            for mol in mols:
                residues = f[mol]['nodes'][()]
                emb_tensor = torch.zeros(len(residues), 1)

                # find a matching pdb path (best-effort)
                try:
                    pdb_file = next(i for i in pdbs if mol in i)
                except StopIteration:
                    pdb_file = None  # not used further

                for i in range(len(residues)):
                    chainID = residues[i][0].decode()
                    resID   = residues[i][1].decode()
                    pt_name = mol + '.' + chainID + '.pt'
                    pt_path = os.path.join(embedding_path, pt_name)
                    res_number = int(resID)
                    try:
                        data = torch.load(pt_path)
                        embedding = data["representations"][33][res_number - 1].mean()
                        emb_tensor[i] = embedding
                    except Exception:
                        emb_tensor[i] = 0

                # write/overwrite dataset
                grp = f[mol].require_group('node_data')
                if 'embedding' in grp:
                    del grp['embedding']
                grp.create_dataset('embedding', data=emb_tensor.numpy())

