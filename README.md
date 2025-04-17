# Protein Design
computational protein design

## 1. Protein Design
### 1.1 Membrane Protein
 - Rosetta: 3UKM
### 1.2 Antibody Design
 - Rosetta: demo Ab
## 2. Docking: Protein - Small Molecule
### 2.1 Protein ~ Ligand
 - Rosetta: G-protein ~ ligand
 - RFDiffusion: design binder of insulin receptor
## 3. Docking: Protein ~ Protein
### 3.1 Antigen-Antibody
 - Rosetta: 3GBN ~ 3GBM

## 5. Datasets
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
