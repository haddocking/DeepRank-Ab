# DeepRank-Ab Inference Pipeline

DeepRank-Ab is a geometric deep learning scoring function for ranking antibody–antigen docking models and predicting **DockQ scores**.

📄 Publication (preprint):  
https://www.biorxiv.org/content/10.64898/2025.12.03.691974v1

This repository provides a **fully automated inference pipeline** from raw PDB → DockQ prediction + quality flags.

---

## 🚀 Key Features

### 🧬 Fully Automated Structure Processing
- Ensemble PDB splitting (MODEL/ENDMDL aware)
- Automatic antibody/antigen chain detection via ANARCI
- Heavy/light chain inference with fallback manual override
- Multi-chain antigen merging into a single chain

### 🧪 Feature Engineering Pipeline
- FASTA generation (CDR + ESM formats)
- ESM-2 embeddings (esm2_t33_650M_UR50D)
- Atom-level graph construction (EGNN-ready)
- CDR annotation via ANARCI
- Region-aware graph features

### 🧠 Deep Learning Inference
- Pretrained EGNN model inference
- DockQ regression output
- Node + edge feature integration:
  - Atom type, polarity, BSA, region, embeddings
  - Voro area, covalent, VdW, orientation

### ⚠️ Structural Quality Filtering
- VdW clash detection (p01–p99 training bounds)
- Heavy–light chain contact validation
- Automatic quality flags in final CSV

---

## 📦 Installation

### 1. Clone repository
```bash
git clone https://github.com/haddocking/DeepRank-Ab
cd DeepRank-Ab
```

### 2. Environment setup
```bash
mamba env create -f environment-gpu.yml
mamba activate deeprank-ab
```

### 3. Install ANARCI
https://github.com/oxpig/ANARCI

Ensure hmmscan is available.

---

## ⚙️ Usage

### Basic command
```bash
python3 scripts/inference.py <pdb_file>
```

### Example
```bash
python3 scripts/inference.py example/test.pdb
```

---

## 🧬 Input Requirements

- PDB file (single model or ensemble supported)
- Optional chain overrides:
  - --heavy_chain_id
  - --light_chain_id
  - --antigen_chain_id

If not provided, chains are auto-detected via ANARCI.

---

## 🔗 Pipeline Overview

1. Workspace creation
2. PDB splitting
3. Chain detection
4. Antigen merging
5. FASTA generation
6. ESM embeddings
7. CDR annotation
8. Graph construction
9. VdW clash filtering
10. Clustering
11. DockQ prediction
12. CSV output

---

## 📊 Outputs

- *_predictions.hdf5
- *.csv

### CSV columns
- pdb_id
- predicted_dockq
- HL_contact_flag
- vdw_clash_flag

---

## ⚠️ Quality Flags

HL contact:
- ok
- low_HL_contacts
- not_applicable

VdW clash:
- ok
- potential_clash

---

## 🧰 Dependencies

- PyTorch
- BioPython
- h5py
- pandas
- numpy
- esm
- ANARCI
- EGNN

---

## 📫 Support

Open a GitHub issue for help.
