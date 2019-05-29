#!/usr/bin/env python

residue_codes = {'ARG':'R','HIS':'H','LYS':'K','ASP':'D','GLU':'E',
	'SER':'S','THR':'T','ASN':'N','GLN':'Q','CYS':'C','SEL':'U','GLY':'G','PRO':'P',
	'ALA':'A','ILE':'I','LEU':'L','MET':'M','PHE':'F','TRP':'W','TYR':'Y','VAL':'V'}

# you must use a list instead of a generator if you are using python 3 and np.in1d
protein_residues = list(residue_codes.keys())
