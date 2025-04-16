#!/usr/bin/env python
# original script: $ROSETTA3/tools/protein_tools/scripts/pdb_renumber.py

from optparse import OptionParser
import sys
from Bio.PDB import PDBParser, PDBIO

from util import load_pdb, universal_open, get_table


usage = "%prog input.pdb output.pdb"
parser= OptionParser(usage)
parser.add_option("-n",
    dest="start",
    help="residue number to start with, default is 1",
    default=1
)
parser.add_option("--preserve",
    dest="preserve",
    help="preserve insertion code and heteroflags",
    default=False,
    action="store_true"
)
parser.add_option("--norestart",
    dest="norestart",
    help="don't start renumbering at each chain, default=False",
    default=False,
    action="store_true"
)
parser.add_option("--keep-table",
    dest="table",
    help="Preserve the rosetta score table at the bottom of the pdb",
    action="store_true",
    default=False
)
(options, args) = parser.parse_args()
print(options)

# validate start number
try:
    residue_id = int(options.start)
except ValueError:
    sys.exit("residue number specified with -n must be an integer")

# get structure from source file
struct = load_pdb(args[0])

# process structure
chain_id = ""
for residue in struct.get_residues():
    chain = residue.get_parent()
    if chain_id != chain.get_id() and not options.norestart:
        chain_id = chain.get_id()
        residue_id = int(options.start)
    #print chain.get_id()
    if options.preserve:
        hetero = residue.id[0]
        insert = residue.id[2]
        residue.id = (hetero, residue_id, insert)
    else:
        residue.id = (' ', residue_id, ' ')
    residue_id += 1

# save data in pdb format
outfile = args[1]
f = universal_open(outfile,'w')
print(f"Try to create {outfile}")

# save atoms
io=PDBIO()
io.set_structure(struct)
io.save(f)

# save table
if(options.table):
    raw_table = get_table(args[0])
    f.writelines(raw_table)

f.close()
