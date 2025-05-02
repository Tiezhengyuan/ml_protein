#!/bin/sh


# your top level rsync directory
MIRRORDIR=$1
# number of parallel download processes
NPROCS=$2
# location of local rsync
RSYNC=rsync

# RCSB PDB server name
SERVER=rsync.wwpdb.org::ftp
# port RCSB PDB server is using
PORT=33444

# get file list, remove first and last 3 lines 
# from output (double-check that these lines are not needed)
${RSYNC} -lpt -v -z --delete --port=$PORT --no-h \
    --list-only ${SERVER}/data/structures/divided/pdb/ \
    | cut -c 44- | head -n -3 | tail -n +3 > dirlist.txt

cat dirlist.txt | xargs -n1 -P${NPROCS} -I% rsync -rlpt -v -z --port=$PORT -P ${SERVER}/data/structures/divided/pdb/% $MIRRORDIR
rm dirlist.txt
