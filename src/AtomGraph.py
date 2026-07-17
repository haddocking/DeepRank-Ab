import networkx as nx
from pdb2sql import interface
from tools.edge_orientation import compute_edge_orientation
from tools.VoroArea import VoronotaAreas

import numpy as np
import os
from Graph import Graph
from tools import BSA
import torch

class AtomGraph(Graph):
    """
    Build a atom-level interface graph from a PDB and decorate it with features.
    Nodes are residues; edges include interface and intra-chain contacts (both directions).
    Optional extras: PSSM alignment, region labels, edge orientations, contact features.
    """
    # TODO: function too complex, refactor
    def __init__(
        self,
        pdb=None,
        contact_distance=5,
        internal_contact_distance=3,
        region_map=None,
        add_orientation=True,
        contact_features=True,
        antigen_chainid="B",
        use_regions=True,
        use_voro = False,
        embedding_path=None,

    ):
        super().__init__()
        # basic ids
        self.type = "atom"
        self.pdb = pdb
        self.name = os.path.splitext(os.path.basename(pdb))[0]

        # region mapping (disable if empty/None)
        self.region_map = region_map 
        self.antigen_chainid = antigen_chainid

        self.contact_distance = contact_distance
        self.internal_contact_distance = internal_contact_distance
        self.use_regions = use_regions
        self.add_orientation = add_orientation
        self.use_voro = use_voro
        self.embedding_path = embedding_path
        self.contact_features = contact_features
        
        self.edge_features = []

        # contact features (dist, elec, vdw, covalent)
        if self.contact_features:
            self.edge_features += ["dist", "elec", "vdw", "covalent"]

        # orientation feature (optional)
        if self.add_orientation:
            self.edge_features.append("orientation")

        # voronota area (optional)
        if self.use_voro:
            self.edge_features.append("voro_area")
        
        self.residue_charge = {
            "CYS": 2.5,  "HIS": -3.20, "ASN": -3.50, "GLN": -3.50, "SER": -0.80,
            "THR": -0.70, "TYR": -1.30, "TRP": -0.90, "ALA": 1.8,  "PHE": 2.8,
            "GLY": -0.4, "ILE": 4.50, "VAL": 4.20, "MET": 1.9,   "PRO": -1.60,
            "LEU": 3.80, "GLU": -3.50, "ASP": -3.50, "LYS": -3.90, "ARG": -4.50,
        }

        # residue index mapping
        self.residue_names = {
            "CYS": 0,  "HIS": 1,  "ASN": 2,  "GLN": 3,  "SER": 4,
            "THR": 5,  "TYR": 6,  "TRP": 7,  "ALA": 8,  "PHE": 9,
            "GLY": 10, "ILE": 11, "VAL": 12, "MET": 13, "PRO": 14,
            "LEU": 15, "GLU": 16, "ASP": 17, "LYS": 18, "ARG": 19,
        }

        # residue polarity category
        self.residue_polarity = {
            "CYS": "polar", "HIS": "polar", "ASN": "polar", "GLN": "polar",
            "SER": "polar", "THR": "polar", "TYR": "polar", "TRP": "polar",
            "ALA": "apolar", "PHE": "apolar", "GLY": "apolar", "ILE": "apolar",
            "VAL": "apolar", "MET": "apolar", "PRO": "apolar", "LEU": "apolar",
            "GLU": "neg_charged", "ASP": "neg_charged",
            "LYS": "pos_charged", "ARG": "pos_charged",
        }

        self.polarity_encoding = {
            "apolar": 0, "polar": 1, "neg_charged": 2, "pos_charged": 3,
        }

        self.atom_type = {
        # Backbone atoms
        "N": 0, "CA": 1, "C": 2, "O": 3,

        # Side-chain atoms
        "CB": 4, "CG": 5, "CG1": 6, "CG2": 7,
        "CD": 8, "CD1": 9, "CD2": 10, "CE": 11,
        "CE1": 12, "CE2": 13, "CE3": 14, "CZ": 15,
        "CZ2": 16, "CZ3": 17, "CH2": 18, "ND1": 19,
        "ND2": 20, "NE": 21, "NE1": 22, "NE2": 23,
        "NH1": 24, "NH2": 25, "NZ": 26, "OD1": 27,
        "OD2": 28, "OE1": 29, "OE2": 30, "OG": 31,
        "OG1": 32, "OH": 33, "SD": 34, "SG": 35,
        }

        self.region_labels = [
            "FR", "L1", "L2", "L3",
            "H1", "H2", "H3",
            "CONST", "AG",
        ]
        self.region2idx = {lab: i for i, lab in enumerate(self.region_labels)}
        self.contact_features = contact_features

        if self.use_voro:
            self.voro = VoronotaAreas(self.pdb)

        db = interface(self.pdb)

        #generate graph 
        self.get_graph(db)

        #get node features
        self.get_node_features(db)
        
        #get edge features
        self.get_edge_features(db)

        
        # close db
        db._close()

    # TODO: function too complex, refactor
    def get_contact_atoms(
        self,
        db,
        cutoff=8.5,
        chain1='A',
        chain2='B',
        only_backbone_atoms=False,
        excludeH=False):
        import itertools

        chainIDs = [chain1, chain2]

        # Load chain data
        xyz, chainid, resnum, resname, atomnum, atName = {}, {}, {}, {}, {}, {}
        for chain in chainIDs:
            data = np.array(db.get('x,y,z,chainID,resName,resSeq,name,serial', chainID=chain))
            xyz[chain]      = data[:, :3].astype(float)
            chainid[chain]  = data[:, 3]
            resname[chain]  = data[:, 4]
            resnum[chain]   = data[:, 5].astype(int)
            atName[chain]   = data[:, 6]
            atomnum[chain]  = data[:, 7].astype(int)

        contact_pairs = {}

        # Iterate over chain pairs (A,B)
        for c1, c2 in itertools.combinations(chainIDs, 2):
            xyz1, xyz2 = xyz[c1], xyz[c2]
            atName1, atName2 = atName[c1], atName[c2]
            resname1, resname2 = resname[c1], resname[c2]
            resnum1, resnum2 = resnum[c1], resnum[c2]
            chainid1, chainid2 = chainid[c1], chainid[c2]
            atomnum1, atomnum2 = atomnum[c1], atomnum[c2]

            # Loop over atoms in chain1
            for i, x0 in enumerate(xyz1):

                if atName1[i].startswith('H'):
                    continue

                # Compute distances to chain2
                distances = np.sqrt(np.sum((xyz2 - x0) ** 2, axis=1))
                contacts = np.where(distances <= cutoff)[0]
                if len(contacts) == 0:
                    continue

                contact_list = []
                for k in contacts:

                    if atName2[k].startswith('H'):
                        continue


                    contact_list.append((
                        str(chainid2[k]),
                        int(resnum2[k]),
                        str(resname2[k]),
                        int(atomnum2[k]),
                        str(atName2[k])
                    ))

                if contact_list:
                    key_atom = (
                        str(chainid1[i]),
                        int(resnum1[i]),
                        str(resname1[i]),
                        int(atomnum1[i]),
                        str(atName1[i])
                    )
                    contact_pairs[key_atom] = contact_list

        return contact_pairs
    
    # TODO: function too complex, refactor
    def get_graph(self, db):
        """
        Build the directed residue graph:
        - interface edges for cross-chain contacts (both directions)
        - internal edges per chain within an atom-distance cutoff (both directions)
        """
        # directed graph
        self.nx = nx.DiGraph()

        # cross-chain contacts
        self.atom_contact_pairs = self.get_contact_atoms(
            db,
            cutoff=self.contact_distance)
        
        all_nodes = self._get_all_valid_nodes(self.atom_contact_pairs)
        
        # interface edges (u,v) and (v,u)
        for key, neighbors in self.atom_contact_pairs.items():
            if key not in all_nodes:
                continue
            for v in neighbors:
                if v not in all_nodes:
                    continue
                d = self._get_edge_distance(key, v, db)
                self.nx.add_edge(key, v, dist=d, type=bytes("interface", encoding="utf-8"))
                self.nx.add_edge(v, key, dist=d, type=bytes("interface", encoding="utf-8"))
        

        # internal edges (both directions)
        edges, dist_list = self.get_internal_edges(db)
        
        for (u, v), d in zip(edges, dist_list):
            if u in all_nodes and v in all_nodes:
                self.nx.add_edge(u, v, dist=d, type=bytes("internal", encoding="utf-8"))
                self.nx.add_edge(v, u, dist=d, type=bytes("internal", encoding="utf-8"))


    # TODO: function too complex, refactor
    def get_node_features(self, db):
        """
        Compute per-node features: chain id, mean 3D position, type, charge,
        polarity, BSA, optional PSSM/IC, optional region label, and optional BioPython features.
        """
        bsa_calc = BSA.BSA(self.pdb, db)
        bsa_calc.get_structure()
        bsa_calc.get_contact_residue_sasa(cutoff=self.contact_distance)
        bsa_data = bsa_calc.bsa_data


        #(chain1, resSeq1, resName1, atomID1, atomName1)

        for node_key in self.nx.nodes:
            chainID = node_key[0]
            resName = node_key[2]
            atomid = node_key[3]
            atomname = node_key[4]
            resSeq = int(node_key[1])

            self.nx.nodes[node_key]["chain"] = {"A": 0, "B": 1}[chainID]
            self.nx.nodes[node_key]["pos"] = db.get("x,y,z", chainID=chainID, resSeq=resSeq, serial=atomid)[0]


            self.nx.nodes[node_key]["res_type"] = self.onehot(
                self.residue_names[resName], len(self.residue_names)
            )
            self.nx.nodes[node_key]["charge"] = self.residue_charge[resName]
            self.nx.nodes[node_key]["polarity"] = self.onehot(
                self.polarity_encoding[self.residue_polarity[resName]],
                len(self.polarity_encoding),
            )
            node_key_bsa = (chainID, resSeq, resName)
            self.nx.nodes[node_key]["bsa"] = bsa_data[node_key_bsa]
            
            self.nx.nodes[node_key]["atom_type"] = self.onehot(
                self.atom_type[atomname], len(self.atom_type)
            )
            # region feature: CDRs on chain A from mapping, AG for chain B
            if self.use_regions:
                if node_key[0] == self.antigen_chainid:
                    ag_idx = self.region2idx["AG"]
                    self.nx.nodes[node_key]["region"] = self.onehot(
                        ag_idx, len(self.region_labels)
                    )
                else:
                    key = (f'{self.name}', resSeq)
                    if key not in self.region_map:
                        raise KeyError(f"missing region mapping for {key}")
                    region_label = self.region_map[key]
                    try:
                        idx = self.region2idx[region_label]
                    except KeyError:
                        raise KeyError(f"unknown region label '{region_label}' for {key}")
                    self.nx.nodes[node_key]["region"] = self.onehot(
                        idx, len(self.region_labels)
                    )


    # TODO: function too complex, refactor
    def get_edge_features(self, db):
        """
        Add edge-level features:
        - directed edge_index
        - optional local-frame orientation per edge
        - optional residue contact features
        """
        # edge index (directed)
        node_keys = list(self.nx.nodes)
        directed_idx = [[node_keys.index(u), node_keys.index(v)] for u, v in self.nx.edges]
        self.nx.edge_index = directed_idx

        for idx, (u, v) in enumerate(self.nx.edges):
            # self.nx.edges[u, v]["rij"] = rij[idx][0]
            if self.use_voro:
                area = self._get_voro_contact_area(u, v)
                self.nx.edges[u, v]["voro_area"] = area

        
        #orientation (optional)
        if self.add_orientation:
            pos_CA, pos_C, pos_N = [], [], []
            for (chain, resseq, resname, atomid, atomname), _ in self.nx.nodes(data=True):
                sel = dict(chainID=chain, resSeq=int(resseq))
                pos_CA.append(db.get("x,y,z", name="CA", **sel)[0])
                pos_C.append(db.get("x,y,z", name="C", **sel)[0])
                pos_N.append(db.get("x,y,z", name="N", **sel)[0])

            pos_CA = torch.tensor(pos_CA, dtype=torch.float)
            pos_C = torch.tensor(pos_C, dtype=torch.float)
            pos_N = torch.tensor(pos_N, dtype=torch.float)
            ei_t = torch.tensor(directed_idx, dtype=torch.long).t()

            orient = compute_edge_orientation(pos_CA, pos_C, pos_N, ei_t).cpu().numpy()
            for idx, (u, v) in enumerate(self.nx.edges):
                self.nx.edges[u, v]["orientation"] = orient[idx]



        from tools.contacts_dr2 import add_residue_contacts_atomgraph
        add_residue_contacts_atomgraph(self, db)
    

    # TODO: function too complex, refactor
    def _get_voro_contact_area(self, node1, node2):
        """
        Compute Voronota contact area between two atoms (nodes).
        Each node is (chainID, resSeq, resName, atomID, atomName).
        """
        try:
            at1 = node1
            at2 = node2
            at1k = (str(at1[3]), str(at1[0]), str(at1[1]), str(at1[2]), str(at1[4]))
            at2k = (str(at2[3]), str(at2[0]), str(at2[1]), str(at2[2]), str(at2[4]))
            return self.voro.contact_areas.get(at1k, {}).get(at2k, 0.0) or \
                self.voro.contact_areas.get(at2k, {}).get(at1k, 0.0) or 0.0
        except Exception:
            return 0.0


    # TODO: function too complex, refactor
    def _get_all_valid_nodes(self, atom_contact_pairs):
        """
        Get all valid atom-level nodes across both chains.
        """
        valid_res = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
            "TYR", "VAL", "ASX", "SEC", "GLX",
        }

        if not isinstance(atom_contact_pairs, dict):
            raise TypeError("Expected atom_contact_pairs to be a dictionary.")

        # Copy to avoid mutating original
        filtered_pairs = dict(atom_contact_pairs)

        # Filter invalid residue names in keys
        keys_to_remove = [
            atom for atom in filtered_pairs.keys() if atom[2] not in valid_res
        ]

        for atom in keys_to_remove:
            filtered_pairs.pop(atom, None)

        # Collect all valid partner atoms
        valid_nodes = set(filtered_pairs.keys())

        for _, partners in filtered_pairs.items():
            for atom in partners:
                if atom[2] in valid_res:
                    valid_nodes.add(atom)

        # Filter out atoms not in self.atom_type
        valid_nodes = [
            atom for atom in valid_nodes if atom[4] in self.atom_type
        ]

        # Return sorted list (for reproducibility)
        return sorted(valid_nodes, key=lambda x: (x[0], x[1], x[3]))

    

    # TODO: function too complex, refactor
    def _get_edge_distance(self, node1, node2, db):
        """
        Minimal atom-atom distance between two residues.
        """
        #        # (chain1, resSeq1, resName1, atomID1, atomName1)
        xyz1 = np.array(db.get("x,y,z", chainID=node1[0], resSeq=node1[1], serial=node1[-2]))
        xyz2 = np.array(db.get("x,y,z", chainID=node2[0], resSeq=node2[1], serial=node2[-2]))
        d2 = (
            -2 * np.dot(xyz1, xyz2.T)
            + np.sum(xyz1**2, axis=1)[:, None]
            + np.sum(xyz2**2, axis=1)
        )
        return np.sqrt(np.min(d2))
    
    # TODO: function too complex, refactor
    def _get_internal_edges_chain(self, nodes, db, cutoff):
        """
        Find intra-chain edges and minimal atom-atom distance between residue pairs.
        """
        nn = len(nodes)
        edges, dist = [], []
        for i1 in range(nn):
            xyz1 = np.array(db.get("x,y,z", chainID=nodes[i1][0], resSeq=nodes[i1][1], serial=nodes[i1][-2]))
            for i2 in range(i1 + 1, nn):
                xyz2 = np.array(db.get("x,y,z", chainID=nodes[i2][0], resSeq=nodes[i2][1], serial=nodes[i2][-2]))
                d2 = (
                    -2 * np.dot(xyz1, xyz2.T)
                    + np.sum(xyz1**2, axis=1)[:, None]
                    + np.sum(xyz2**2, axis=1)
                )
                if np.any(d2 < cutoff**2):
                    edges.append((nodes[i1], nodes[i2]))
                    dist.append(np.sqrt(np.min(d2)))
        return edges, dist
    
    def onehot(self, idx, size):
        """One-hot encode an integer index."""
        onehot = torch.zeros(size)
        onehot[idx] = 1
        return np.array(onehot)


    # TODO: function too complex, refactor
    def get_internal_edges(self, db):
        """
        Collect intra-chain edges for both chains using the atom-level cutoff.
        Returns edge tuples and their minimal atom-atom distances.
        """
        nodesA, nodesB = [], []
        for n in self.nx.nodes:
            if n[0] == "A":
                nodesA.append(n)
            elif n[0] == "B":
                nodesB.append(n)

        edgesA, distA = self._get_internal_edges_chain(nodesA, db, self.internal_contact_distance)
        edgesB, distB = self._get_internal_edges_chain(nodesB, db, self.internal_contact_distance)
        return edgesA + edgesB, distA + distB
        
    
