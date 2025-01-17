# ml_protein
computational protein



cd data/protein_folding


# For `assembly1` complexes, use the PDB's `20240101` AWS snapshot:
aws s3 sync s3://pdbsnapshots/20240101/pub/pdb/data/assemblies/mmCIF/divided/ ./pdb_data/unfiltered_assembly_mmcifs
# Or as a fallback, use rsync:
rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/assemblies/mmCIF/divided/ ./pdb_data/unfiltered_assembly_mmcifs/

# For asymmetric unit complexes, also use the PDB's `20240101` AWS snapshot:
aws s3 sync s3://pdbsnapshots/20240101/pub/pdb/data/structures/divided/mmCIF/ ./pdb_data/unfiltered_asym_mmcifs
# Or as a fallback, use rsync:
rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ ./pdb_data/unfiltered_asym_mmcifs/