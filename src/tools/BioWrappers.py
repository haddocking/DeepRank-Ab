import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import ResidueDepth, get_surface, residue_depth
from Bio.PDB.HSExposure import HSExposureCA

import warnings
from Bio import BiopythonWarning

import tempfile

with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)
    from Bio import SearchIO

from time import time


def get_bio_model(pdbfile):
    """Get the model

    Args:
        pdbfile (str): pdbfile

    Returns:
        [type]: Bio object
        
    """

    """
    At the very start, you need to turn the text of a PDB file into an in-memory object you can 
    interrogate. Biopython’s PDBParser does exactly that:

    You give it the path to a .pdb file.

    It reads all the ATOM/HETATM records, builds up a hierarchy of Structure → Model → Chain → 
    Residue → Atom, and returns that structure.

    Since most PDBs have a single “Model 0,” the helper immediately returns structure[0] 
    so that downstream code works purely against a Model object.
    """
    parser = PDBParser()
    structure = parser.get_structure('_tmp', pdbfile)
    return structure[0]

# Result: a Biopython Model that you can index by chain IDs and residue numbers to grab 
# any residue or atom you like.


# TODO: function too complex, refactor
def get_depth_res(model):
    """Get the residue Depth

    Args:
        model (bio model): model of the strucrture

    Returns:
        dict: depth res
    """
    """
    Once you have the Model, you often want to know how deeply buried each residue is inside the 
    protein core. The ResidueDepth class from Biopython wraps an algorithm that:
    Computes the molecular surface of the whole structure,

    For each residue, finds all its constituent atoms,

    Measures the distance from each atom back to the solvent-accessible surface, and

    Defines the residue’s “depth” as the minimum over its atoms.

    """
    rd = ResidueDepth(model)

    data = {}
    t0 = time()
    for k in list(rd.keys()):
        new_key = (k[0], k[1][1])
        data[new_key] = rd[k][0]

    return data
# Result: a plain Python dict mapping e.g. ("A", 123) → 2.8 (Å), telling you residue 123 
# on chain A sits 2.8 Å beneath the surface.



# TODO: function too complex, refactor
def get_depth_contact_res(model, contact_res):
    """Get the residue Depth only for a specific set of residues

    Args:
        model (bio model): model of the strucrture
        contact_res (list): list of contact residues

    Returns:
        dict: depth res
    """

    """
    Sometimes you only care about a small set of residues—say those you found experimentally 
    to touch another molecule. This helper:

Calls get_surface(model) to compute exactly the same molecular surface that ResidueDepth 
would have used.

Iterates over your provided list of contact-residue identifiers (each also a (chain, resnum) tuple).

For each, grabs the corresponding Residue object from the Model, and calls the standalone 
function residue_depth(res, surface) to measure that single residue’s depth against the precomputed surface.

Returns a dict mapping each of your input residues to its depth.
    """

    surface = get_surface(model)
    data = {}
    for r in contact_res:
        chain = model[r[0]]
        res = chain[r[1]]
        data[r] = residue_depth(res, surface)
    return data

# Result: like get_depth_res, but only for the residues you care about—cheaper 
# if you only have a handful.
def get_hse(model):
    """Get the hydrogen surface exposure

    Args:
        model (bio model): model of the strucrture

    Returns:
        dict: hse data
    """
    """
    Hydrogen-shell exposure (HSE) is another way to quantify how exposed a residue’s Cα 
    atom is, by counting how many of its immediate neighbors lie above or below the Cα 
    plane. Biopython’s HSExposureCA class handles this:
    You give it the full Model and it builds an internal map of every residue’s Cα environment.

Iterating its keys (again (chain, Residue)), you get for each a tuple (up_count, down_count, raw_angle).

The helper rekeys each entry to (chain, residue_number), and ensures you never get a None angle by substituting 0.0.
    """

    hse = HSExposureCA(model)
    data = {}
    for k in list(hse.keys()):
        new_key = (k[0], k[1][1])

        x = hse[k]
        if x[2] is None:
            x = list(x)
            x[2] = 0.0
            x = tuple(x)

        data[new_key] = x
    return data
# esult: a dict where each residue maps to a small 3-tuple, e.g. ("B", 45) → (5, 3, 28.7), 
# meaning “5 neighbors above, 3 below, angle 28.7°.”