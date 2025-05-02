#!/usr/bin/env python

'''
Employ the package PyMol to calculate RMSD
'''

from optparse import OptionParser
import sys
import subprocess
from sys import argv
from array import *
import os
import copy
from pprint import pprint
import math
import gzip
import json
import numpy as np



# initialize pymol in headless mode without GUI
import __main__
__main__.pymol_argv = ['pymol', '-qei']

import pymol
# pymol.finish_launching()

from pymol import cmd
from pymol import stored
from pymol import selector
pymol.cmd.feedback('disable', 'all', 'actions')
pymol.cmd.feedback('disable', 'all', 'results')



# for importing as a plugin into PyMol
### Parse the options ###
parser=OptionParser(usage='calculate_rmsd.py [options] <list of pdbs>', description="\
This script calculates RMS to native for each provided model.\
It does so for 1) the whole complex, 2) all HETATMs or specified chains,\
3) '2' + side-chains within specified angstroms\n\
")

from sys import argv
if len(argv) == 1 : # only the script name
	argv.append('--help')

### Add options
parser.add_option("-n", "--native", dest="native", type='string',
		help="path to the native PDB structure")
parser.add_option("-c", "--chains", dest="chains", type='string',
		help="comma separated list of 1-char chains as found in PDB")
parser.add_option("-a", "--angstroms", dest="angstroms", type='string', default="4",
		help="protein residues within this many angstroms will be included in a separate RMS calculation")
parser.add_option("-o", "--output", dest="output", type='string', default="out.txt",
		help="path to output file")
parser.add_option("--include_hydrogens", dest="include_hydrogens", default="False", action="store_true",
		help="provide this flag to include Hs in the calculations")
parser.add_option("--symmetrical_chains", dest="symmetrical_chains", type='string',
		help="comma separated list of protein chains that are super_imposeable. "+ \
		"This script will report the minimum RMS from aligning each pairing of superimposeable chains")
parser.add_option("--no_align", dest="do_align", 
	default="True", action="store_false", #default does align by protein before RMSD
		help="provide this flag to calculate RMSD without doing alignment")

(options, args)=parser.parse_args()



def optAlign( sel1, sel2 ):
	"""
	optAlign performs the Kabsch alignment algorithm upon the alpha-carbons of two selections.
	Example:   optAlign MOL1 and i. 20-40, MOL2 and i. 102-122
	Example 2: optAlign 1GGZ and i. 4-146 and n. CA, 1CLL and i. 4-146 and n. CA
 
	Two RMSDs are returned.  One comes from the Kabsch algorithm and the other from
	PyMol based upon your selections.
 
	By default, this program will optimally align the ALPHA CARBONS of the selections provided.
	To turn off this feature remove the lines between the commented "REMOVE ALPHA CARBONS" below.
 
	@param sel1: First PyMol selection with N-atoms
	@param sel2: Second PyMol selection with N-atoms
	"""
	cmd.reset()
 
	# make the lists for holding coordinates
	# partial lists
	stored.sel1 = []
	stored.sel2 = []
	# full lists
	stored.mol1 = []
	stored.mol2 = []
 
	# Get the selected coordinates.  We
	# align these coords.
	cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
	cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
 
	# get molecule name
	mol1 = cmd.identify(sel1,1)[0][0]
	mol2 = cmd.identify(sel2,1)[0][0]
 
	# Get all molecule coords.  We do this because
	# we have to rotate the whole molcule, not just
	# the aligned selection
	cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
	cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")
 
	# check for consistency
	print(len(stored.sel1))
	print(len(stored.sel2))
	assert len(stored.sel1) == len(stored.sel2)
	L = len(stored.sel1)
	assert L > 0
 
	# must alway center the two proteins to avoid
	# affine transformations.  Center the two proteins
	# to their selections.
	COM1 = np.sum(stored.sel1,axis=0) / float(L)
	COM2 = np.sum(stored.sel2,axis=0) / float(L)
	stored.sel1 -= COM1
	stored.sel2 -= COM2
 
	# Initial residual, see Kabsch.
	E0 = np.sum(
		np.sum(stored.sel1 * stored.sel1,axis=0),axis=0) + np.sum( np.sum(stored.sel2 * stored.sel2,axis=0),
		axis=0
	)
 
	#
	# This beautiful step provides the answer.  V and Wt are the orthonormal
	# bases that when multiplied by each other give us the rotation matrix, U.
	# S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
	V, S, Wt = np.linalg.svd( np.dot( np.transpose(stored.sel2), stored.sel1))
 
	# we already have our solution, in the results from SVD.
	# we just need to check for reflections and then produce
	# the rotation.  V and Wt are orthonormal, so their det's
	# are +/-1.
	reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))
 
	if reflect == -1.0:
		S[-1] = -S[-1]
		V[:,-1] = -V[:,-1]
 
	RMSD = E0 - (2.0 * sum(S))
	RMSD = np.sqrt(abs(RMSD / L))
 
	#U is simply V*Wt
	U = np.dot(V, Wt)
 
	# rotate and translate the molecule
	stored.sel2 = np.dot((stored.mol2 - COM2), U)
	stored.sel2 = stored.sel2.tolist()
	# center the molecule
	stored.sel1 = stored.mol1 - COM1
	stored.sel1 = stored.sel1.tolist()
 
	# let PyMol know about the changes to the coordinates
	cmd.alter_state(1,mol1,"(x,y,z)=stored.sel1.pop(0)")
	cmd.alter_state(1,mol2,"(x,y,z)=stored.sel2.pop(0)")
 
	print("RMSD=%f" % RMSD)
 
	# we already have our solution, in the results from SVD.
	# we just need to check for reflections and then produce
	# the rotation.  V and Wt are orthonormal, so their det's
	# are +/-1.
	reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))
 
	if reflect == -1.0:
		S[-1] = -S[-1]
		V[:,-1] = -V[:,-1]
 
	RMSD = E0 - (2.0 * sum(S))
	RMSD = np.sqrt(abs(RMSD / L))
 
	#U is simply V*Wt
	U = np.dot(V, Wt)
 
	# rotate and translate the molecule
	stored.sel2 = np.dot((stored.mol2 - COM2), U)
	stored.sel2 = stored.sel2.tolist()
	# center the molecule
	stored.sel1 = stored.mol1 - COM1
	stored.sel1 = stored.sel1.tolist()
 
	# let PyMol know about the changes to the coordinates
	cmd.alter_state(1,mol1,"(x,y,z)=stored.sel1.pop(0)")
	cmd.alter_state(1,mol2,"(x,y,z)=stored.sel2.pop(0)")
 
	print("RMSD=%f" % RMSD)
 
cmd.extend("optAlign", optAlign)


######################################
def load_file(file, name):
	print(f"load file {file} with the name = {name}:")
	pymol.cmd.load(file, name)
	# call all waters by the same name
	pymol.cmd.alter(selection='resn TP3', expression="resn='WAT'")

def get_hydro():
	if options.include_hydrogens:
		return 'and not hydro'
	return ''

def make_selections(file):
	'''
	select atoms using pymol
	'''
	# load pdb
	load_file(options.native, 'native')
	load_file(file, 'docked')
	# set hydrogen atom
	hydro=get_hydro()

	# select ligand atoms
	if options.chains==None:
		pymol.cmd.select(name='docked_lig', selection=f'docked and hetatm {hydro}')
		pymol.cmd.select(name='native_lig', selection=f'native and hetatm {hydro}')
	else:
		chains= options.chains.replace(',','+')
		pymol.cmd.select(
			name='docked_lig',
			selection = f'docked and chain {chains} {hydro}'
		)
		pymol.cmd.select(
			name='native_lig',
			selection = f'native and chain {chains} {hydro}'
		)
	pymol.cmd.select(
		name='docked_lig_sc',
		selection = f'(docked within {options.angstroms} of docked_lig) {hydro}'
	)
	pymol.cmd.select(
		name='native_lig_sc',
		selection = f'(native within {options.angstroms} of native_lig) {hydro}'
	)

def align_for_min_rms(file):
	'''
	align atoms by minimizing RMSD
	'''
	if options.symmetrical_chains is None:
		pymol.cmd.align('docked and not hetatm', 'native and not hetatm')
		return None

	best_rms = 100
	best_alignement = (None,None)
	symmetrical_chains = options.symmetrical_chains.split(',')
	for i in range(len(symmetrical_chains)-1):
		for j in range(i+1, len(symmetrical_chains)):
			pymol.cmd.align(
				'docked and not hetatm and chain ' + symmetrical_chains[i],
				'native and not hetatm and chain ' + symmetrical_chains[j]
			)
			# calculate ligan RMSD
			this_rms = pymol.cmd.rms_cur('docked_lig', 'native_lig')
			if this_rms < best_rms:
				best_rms = this_rms
				best_alignment = (i, j)
	pymol.cmd.align(
		'docked and not hetatm and chain '+ symmetrical_chains[best_alignment[0]],
		'native and not hetatm and chain '+ symmetrical_chains[best_alignment[1]]
	)

def calculate_rms(file):
	# select atoms
	make_selections(file)

	if options.do_align:
		align_for_min_rms(file)

	#
	all_rms= pymol.cmd.rms_cur('docked', 'native')
	lig_and_sc_rms= pymol.cmd.rms_cur('docked_lig_sc', 'native_lig_sc')
	cal_rms = {
		'all': all_rms,
		'ligand_sc': lig_and_sc_rms,
	}
	

	if options.chains is None:
		hetatm_rms = pymol.cmd.rms_cur('docked_lig', 'native_lig')
		cal_rms['hetatm'] = round(hetatm_rms, 2)
	else:
		chains= options.chains.replace(',','+')
		hydro=get_hydro()
		for chain in chains.split('+'):
			# rmsd by chain
			ligand_rms= pymol.cmd.rms_cur(
				'docked and chain '+chain+' '+hydro,
				'native and chain '+chain+' '+hydro
			)
			cal_rms[f"chain_{chain}_ligand"] = ligand_rms

			#Calculate centroid distance if centroid mode
			docked_cent = pymol.cmd.get_coords('docked_lig')
			native_cent = pymol.cmd.get_coords('native_lig')
			dist = np.square(docked_cent - native_cent)
			cal_rms[f'chain_{chain}_centroid'] = dist.sum()**0.05

			# provide this flag to calculate ligand RMSD based on conformer differences 
			#only by aligning ligands onto each other
			#Align ligand then do RMSD Superimposed rmsd
			# optAlign('docked_lig', 'native_lig')
			align_ligand_rms= pymol.cmd.rms_cur('docked_lig', 'native_lig')
			cal_rms[f'chain_{chain}_superimposed'] = align_ligand_rms

		if len(chains)>1:
			all_chain_rms = pymol.cmd.rms_cur('docked_lig', 'native_lig')
			cal_rms['all_chain'] = all_chain_rms

	return {k:round(v, 2) for k,v in cal_rms.items()}

def parse_scores(file):
	'''
	retrieve scores from pdb
	'''
	scores = {}
	with open(file, 'r') as f:
		for line in f:
			line = line.strip()
			items = line.split()
			if len(items) == 2:
				key = items[0]
				if key.startswith('ligand') or key.startswith('total'):
					scores[key] = items[1]
	return scores

#########################################
print('###options:', options)
print('###files:', args)

data = {
	'native': {
		'file': options.native,
	},
	'docked': [],
}
for file in args:
	rec = {
		'file': file,
		'rmsd': calculate_rms(file),
		'scores': parse_scores(file),
	}
	data['docked'].append(rec)
	pymol.cmd.delete('all')
pprint(data['docked'][-1])

outfile = f'{options.output}.json'
with open(outfile, 'w') as f:
	json.dump(data, f, indent=4, sort_keys=True)
# Quit PyMOL
cmd.quit()