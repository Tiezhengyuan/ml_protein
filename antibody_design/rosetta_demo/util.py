'''
utils
'''
import math
import gzip
import sys

from Bio.PDB import PDBParser

def universal_open(path, mode):
    """
	given a path and a mode, return a filehandle.
    allows you to open files regardless of whether they are gzipped
	"""
    extension = path.split(".")[-1]
    if(extension == "gz"):
        try:
            return gzip.open(path, mode)
        except IOError:
            sys.exit("could not open "+path+" in mode "+mode)
    else:
        try:
            return open(path, mode)
        except IOError:
            sys.exit("could not open "+path+" in mode "+mode)

def load_pdb(path):
	"""
	return a biopython structure object given a pdb file path
	"""
	print('###path:', path)
	pdb_file = universal_open(path, 'r')
	parser = PDBParser(PERMISSIVE=1)
	structure = parser.get_structure(path[0:4], pdb_file)
	pdb_file.close()
	return structure

def get_table(path):
    """
	return the score table from the bottom of a PDB as a list of lines
	"""
    raw_table = []
    infile = universal_open(path,'r')
    table = False
    for line in infile:
        line_split = line.split()
        if len(line_split) < 1:
            break
        if line_split[0] == "#BEGIN_POSE_ENERGIES_TABLE":
            table = True
            raw_table.append(line)
        #elif table and line_split[0] == "#END_POSE_ENERGIES_TABLE":
        #    raw_table.append(line)
        #    break
        elif table:
            raw_table.append(line)
    infile.close()
    return raw_table