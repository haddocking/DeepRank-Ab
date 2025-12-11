import json
import warnings
from pathlib import Path
from Bio.PDB import PDBParser
from Bio import SeqIO, pairwise2
from anarci import anarci
from itertools import chain as chain_
from multiprocessing import Pool

import sys 
BASE = Path(__file__).resolve().parent 
hmmscan_path = BASE / "ANARCI"
print("Using ANARCI path:", hmmscan_path)




# residue code map (includes MSE -> M)
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M"
}

# IMGT CDR ranges
CDR_L = {"L1": range(27, 39), "L2": range(56, 66), "L3": range(105, 118)}
CDR_H = {"H1": range(27, 39), "H2": range(56, 66), "H3": range(105, 118)}


def extract_chain_sequence(chain) -> str:
    """Return one-letter amino acid sequence for a chain."""
    return ''.join(THREE_TO_ONE.get(res.get_resname(), '') for res in chain)


def region_l(position: int) -> str:
    """Label light-chain IMGT position as FR or CDR."""
    for label, rng in CDR_L.items():
        if position in rng:
            return label
    return "CONST"  # fallback

def region_h(position: int) -> str:
    """Label heavy-chain IMGT position as FR or CDR."""
    for label, rng in CDR_H.items():
        if position in rng:
            return label
    return "CONST"

def locate(subseq: str, fullseq: str, tag: str):
    """Local alignment of subseq in fullseq; return (start, end) indices in fullseq."""
    alignment = pairwise2.align.localms(fullseq, subseq, 2, -1, -0.5, -0.1, one_alignment_only=True)
    if not alignment:
        raise RuntimeError(f"{tag}: alignment failed")
    _, _, _, start, end = alignment[0]
    return start, end




def annotate_single_pdb(fasta_file: Path, pdb_file: Path, antigen_chainid: str = 'A'):
    """
    Annotate CDR/FR/CONST for a single PDB using a reference fasta.

    Args:
        fasta_file: Path to fasta file with exactly two sequences ('H' and 'L').
        pdb_file: Path to PDB file to annotate.
        output_dir: Base output directory to store annotations.
    """
    parser = PDBParser(QUIET=True)
    seqs = [record for record in SeqIO.parse(fasta_file, "fasta") if record.seq]
    if len(seqs) != 2:
        raise ValueError(f"{fasta_file}: expected 2 sequences (H and L), found {len(seqs)}")

    ref = {}
    for rec in seqs:
        role = rec.id.strip().upper()[-1]
        if role not in ("H", "L"):
            raise ValueError(f"{fasta_file}: unexpected record id '{rec.id}'")
        ref[role] = str(rec.seq)

    seqH, seqL = ref["H"], ref["L"]

    model_name = pdb_file.stem
    struct = parser.get_structure(model_name, str(pdb_file))
    if 'A' not in struct[0]:
        raise ValueError(f"{model_name}: chain A not found")
    #antigen named chain A; H+L named chain B
    all_chainids = [chain.id for chain in struct[0]]
    chainid_antibody = [chainid for chainid in all_chainids if chainid != antigen_chainid]
    chainA_seq = extract_chain_sequence(struct[0][chainid_antibody[0]])

    # locate H and L within chain A
    try:
        h0, _ = locate(seqH, chainA_seq, f"{model_name}_H")
        l0, _ = locate(seqL, chainA_seq, f"{model_name}_L")
        #print(f"{model_name}: H at {h0}-{h0+len(seqH)}, L at {l0}-{l0+len(seqL)} in chain B")


        # adjust start if N-term residues are missing
        orig_h0 = h0
        while h0 > 0 and chainA_seq[h0] != seqH[0]:
            h0 -= 1
        if h0 != orig_h0:
            warnings.warn(f"{model_name}: heavy missing {orig_h0 - h0} N-term residues, shifting start")

        orig_l0 = l0
        while l0 > 0 and chainA_seq[l0] != seqL[0]:
            l0 -= 1
        if l0 != orig_l0:
            warnings.warn(f"{model_name}: light missing {orig_l0 - l0} N-term residues, shifting start")
    

    except RuntimeError as e:
        warnings.warn(str(e))
        return

    # variable domains (~first 130 aa)
    varH, varL = seqH[:130], seqL[:130]
    records = []
    mapping = []

    if varH:
        records.append((f"{model_name}_H", varH))
        mapping.append({'model': model_name, 'chain': 'H', 'start': h0, f'chain{chainid_antibody[0]}': chainA_seq})
    if varL:
        records.append((f"{model_name}_L", varL))
        mapping.append({'model': model_name, 'chain': 'L', 'start': l0, f'chain{chainid_antibody[0]}': chainA_seq})

    if not records:
        warnings.warn(f"{model_name}: no variable domains found")
        return


    # run ANARCI with IMGT numbering
    numbering, _, _ = anarci(records, scheme="imgt", assign_germline=True, output=False, hmmerpath=hmmscan_path)

    annotations = {}
    for info, num in zip(mapping, numbering):
        model = f"{info['model']}.pdb"
        chain = info['chain']
        start = info['start']
        seq = info[f'chain{chainid_antibody[0]}']

        annotations.setdefault(model, [(aa, 'CONST') for aa in seq])
        ann = annotations[model]
        if not num:
            continue

        aligned = num[0][0]  # list of ((position, ins), aa)
        region_fn = region_h if chain == 'H' else region_l

        idx = start


        for ((pos, _), aa) in aligned:
            if aa == '-':
                continue


            while ann[idx][0] != aa:
                idx += 1
            ann[idx] = (aa, region_fn(pos))
            idx += 1
    

    return annotations


def annotate_folder(folder: Path, output_dir: Path, fasta_file: Path = None, antigen_chainid: str = 'A'):
    """
    Annotate all PDB files in a folder using corresponding fasta files.

    Args:
        folder: Path to folder containing PDB and fasta files.
        output_dir: Base output directory to store annotations.
    """
    pdb_files = list(folder.glob("*.pdb"))
   
    all_annotations = {}
    for pdb_file in pdb_files:
        try:
            annotations = annotate_single_pdb(fasta_file, pdb_file, antigen_chainid)
            if annotations is None:
                continue
            all_annotations.update(annotations)
        except:
            continue

    #dump all annotations
    output_dir = folder / output_dir
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    output_file = output_dir/ "annotations_cdrs.json"
    with open(output_file, 'w') as f:
        json.dump(all_annotations, f, indent=2)
    
    print(f"Annotations written to {output_file}")


def _annotate_single_wrapper(args):
    """Helper to call annotate_single_pdb with try/except in multiprocessing."""
    fasta_file, pdb_file, antigen_chainid = args    

    annotations = annotate_single_pdb(fasta_file, pdb_file, antigen_chainid)
    return annotations

def annotate_folder_one_by_one_mp(folder: Path, fasta_folder: Path, output_dir: Path, n_cores: int = None, antigen_chainid: str = 'A'):
    """ 
    Annotate all PDB files in a folder using corresponding fasta files in parallel.

    Args:
        folder: Path to folder containing PDB files.
        fasta_folder: Path to folder containing fasta files.
        output_dir: Base output directory to store annotations.
        n_processes: Number of parallel processes. Defaults to cpu_count().
    """
    pdb_files = list(folder.glob("*.pdb"))
    tasks = []

    for pdb_file in pdb_files:
        fasta_file = fasta_folder / f"{pdb_file.stem}_HL.fasta"
        if not fasta_file.exists():
            warnings.warn(f"{fasta_file} not found, skipping {pdb_file}")
            continue
        tasks.append((fasta_file, pdb_file, antigen_chainid))

    n_processes = n_cores if n_cores else os.cpu_count()

    all_annotations = {}
    if tasks:
        with Pool(processes=n_processes) as pool:
            results = pool.map(_annotate_single_wrapper, tasks)


        for res in results:
            all_annotations.update(res)

    # dump all annotations
    output_dir = folder / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "annotations_cdrs.json"
    with open(output_file, 'w') as f:
        json.dump(all_annotations, f, indent=2)

    print(f"Annotations written to {output_file}")

# folder_path = '/trinity/login/xxu/02_rankab/external_test_set/arne_e_set/redo/data/8f5n'
# fasta_path = '/trinity/login/xxu/02_rankab/external_test_set/arne_e_set/redo/data/8f5n/fasta'
# annotate_folder_one_by_one(Path(folder_path), Path('annotations'), Path(fasta_path))



def annotate_folder_one_by_one_mp_single_fasta(folder: Path, fasta_file_path: Path, output_dir: Path, n_cores: int = None, antigen_chainid: str = 'A'):
    """ 
    Annotate all PDB files in a folder using corresponding fasta files in parallel.

    Args:
        folder: Path to folder containing PDB files.
        fasta_folder: Path to folder containing fasta files.
        output_dir: Base output directory to store annotations.
        n_processes: Number of parallel processes. Defaults to cpu_count().
    """
    pdb_files = list(folder.glob("*.pdb"))
    tasks = []

    for pdb_file in pdb_files:
        # fasta_file = fasta_folder / f"{pdb_file.stem}.fasta"
        # if not fasta_file.exists():
        #     warnings.warn(f"{fasta_file} not found, skipping {pdb_file}")
        #     continue
        tasks.append((fasta_file_path, pdb_file, antigen_chainid))

    n_processes = n_cores if n_cores else os.cpu_count()

    all_annotations = {}
    if tasks:
        with Pool(processes=n_processes) as pool:
            results = pool.map(_annotate_single_wrapper, tasks)

        for res in results:
            all_annotations.update(res)

    # dump all annotations
    output_dir = folder / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "annotations_cdrs.json"
    with open(output_file, 'w') as f:
        json.dump(all_annotations, f, indent=2)

    print(f"Annotations written to {output_file}")



#test single 


# fasta_f = '/trinity/login/xxu/02_rankab/diff_analysis/fasta_backup/8s6z_antibody_clean.fasta'
# pdbf = '/trinity/login/xxu/02_rankab/diff_analysis/af3_100seed_data/8s6z/rank-111_seed63-2_merged.pdb'
# annotate_single_pdb(Path(fasta_f), Path(pdbf), antigen_chainid='B')
