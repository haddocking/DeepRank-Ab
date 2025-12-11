# DeepRank-Ab Inference Pipeline

DeepRank-Ab is a scoring function for ranking antibody-antigen models. 

Refer to our publication for more details at: 

https://www.biorxiv.org/content/10.64898/2025.12.03.691974v1


This repository provides full inference pipleline for model described in the paper. 


---

## Features

* **PDB Processing:** Split ensemble structures, extract chain sequences, and merge chains for analysis.
* **FASTA Conversion:** Generate FASTA for annotation and ESM-compatible merged sequences.
* **ESM Embeddings:** Compute embeddings using `esm2_t33_650M_UR50D`.
* **Graph Generation:** Build atom-level graphs with pre-computed node and edge features.
* **Prediction:** Evaluate models with pretrained weights and output predicted DockQ scores.

---

## Installation

1. Clone this repository:

```bash
git clone https://github.com/haddocking/DeepRank-Ab
cd DeepRank-Ab
```

2. Create and activate a mamba environment: 

```bash
mamba env create -f environment-gpu.yml
mamba activate deeprank-ab
```

3. Install ANARCI

```bash
git clone https://github.com/oxpig/ANARCI
cd ANARCI
python setup.py install 
```





## Usage

The main entry point is `DeepRank-Ab/scripts/inference.py`:


```bash
python3 scripts/inference.py <pdb_file> <antibody_heavy_chain_id> <antibody_light_chain_id> <antigen_chain_id>
```

### Example

```bash
python3 scripts/inference.py example/test.pdb H L A 
```

This command will create a workspace, generate embeddings, annotate CDRs, build graphs, cluster nodes, predict DockQ scores, and save outputs in CSV and HDF5 formats.

---

## Input

* **PDB file:** Structure of the antibody-antigen complex. could be either a single model or an ensemble of models.
* **Heavy chain ID:** Chain identifier for the antibody heavy chain (e.g., `H`).
* **Light chain ID:** Chain identifier for the antibody light chain (e.g., `L`).
* **Antigen chain ID:** Chain identifier for the antigen (e.g., `A`).

---

## Output

* **HDF5 predictions:** Neural network predictions for DockQ scores.
* **CSV predictions:** Columns: `pdb_id,predicted_dockq`.
* **Graphs:** HDF5 graph representations of the complex.
* **Embeddings:** Per-chain ESM embeddings used for prediction.

---

## Pipeline Overview

```text
PDB Input
   │
   ├─> Split ensembles
   │
   ├─> Convert chains to FASTA
   │
   ├─> Generate ESM embeddings
   │
   ├─> Annotate CDRs
   │
   ├─> Build atom-level graph with node & edge features
   │
   ├─> Add embeddings to graph
   │
   ├─> Cluster nodes (MCL)
   │
   └─> Predict DockQ scores using pretrained EGNN
       │
       └─> CSV/HDF5 predictions
```


## 