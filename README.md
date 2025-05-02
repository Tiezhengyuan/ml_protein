# Protein Design
computational protein design

## 1. Protein Design
### 1.1 Membrane Protein
 - Rosetta: 3UKM
### 1.2 Antibody Design
 - Rosetta: demo Ab
### 1.3 Design Backbone
 - design monomer using RFDiffusion
 - design homotrimer using ProteinMPNN
### 1.4 motif scaffolding
 - design motif scaffolding
## 2. Docking: Protein - Small Molecule
### 2.1 Protein ~ Ligand
 - Rosetta: G-protein ~ ligand
 
## 3. Docking: Protein ~ Protein
### 3.1 Antigen-Antibody
 - Rosetta: 3GBN ~ 3GBM
### 3.2 design binder
 - denovo design for insulin receptor 3W11s using RFDiffusion
### 3.3 scaffolding
 - insulin receptor 3W11
## 4. Datasets

### 4.1
 - protein-ligand complex from PDB

### 4.5
cd data/protein_folding

For `assembly1` complexes, use the PDB's `20240101` AWS snapshot:
```
aws s3 sync s3://pdbsnapshots/20240101/pub/pdb/data/assemblies/mmCIF/divided/ ./pdb_data/unfiltered_assembly_mmcifs
```
Or as a fallback, use rsync:
```
rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/assemblies/mmCIF/divided/ ./pdb_data/unfiltered_assembly_mmcifs/
```
For asymmetric unit complexes, also use the PDB's `20240101` AWS snapshot:
```
aws s3 sync s3://pdbsnapshots/20240101/pub/pdb/data/structures/divided/mmCIF/ ./pdb_data/unfiltered_asym_mmcifs
```
Or as a fallback, use rsync:
```
rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ ./pdb_data/unfiltered_asym_mmcifs/
```
