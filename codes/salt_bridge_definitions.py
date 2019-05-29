#!/usr/bin/env python

import numpy as np

class SaltBridge:
	#!!! this is deprecated by the version below. it should be retired
	"""
	Identify salt bridges.
	"""
	"""
	A note on other software definitions.
	VMD seems permissive: 
		https://www.ks.uiuc.edu/Research/vmd/plugins/saltbr/
	> A salt bridge is considered to be formed if the distance between any of 
	the oxygen atoms of acidic residues and the nitrogen atoms of basic residues 
	are within the cut-off distance (default 3.2 Angstroms) in at least one 
	frame. The default distance cut-off can be changed by the user. This plugin 
	does not attempt to identify hydrogen bonds.
	#! see also 
		http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/
	#! see also 
		https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6016756/
	> Salt-bridge and aromatic-aromatic bond calculation
	>Positively charged atoms are Lys NZ1, Arg NH1, Arg NH2 and His NE2. 
	Negatively charged atoms are Asp OD1, Asp OD2, Glu OE1, Glu OE2. A salt 
	bridge is defined if two oppositely charged atoms lie within 4A across the 
	interface. A pi-pi interaction is defined between two aromatic amino acids 
	if the distance calculated from their centroid is less than 7.5A.
	"""
	defn_salt_bridges = {
		'donor':[
			{'resname':'ARG','atoms':['NH1','NH2']},
			{'resname':'LYS','atoms':['NZ']},
			# we allow HIS which is charged at neutral 
			#   but we discard NE2 which is protonated with HE2
			{'resname':'HIS','atoms':['ND1','NE2'][:1]}],
		'acceptor':[
			{'resname':'GLU','atoms':['OE1','OE2']},
			{'resname':'ASP','atoms':['OD1','OD2']},]}
	def __init__(self): pass
	def filter_salt_bridges(self,bonds,rowspec):
		"""Filter a bonds list for salt bridges."""
		# search should be fast because it is composed of a bunch of filters on bonds
		#! however there are many combinations. are there ways to improve?
		subject_target_combos = [['subject','target'],['target','subject']]
		salt_filter = np.any([
			np.all([
				bonds[:,rowspec.index('%s_resname'%this)]==salt_bridge_donor['resname'],
				np.in1d(bonds[:,rowspec.index('%s_atom'%this)],salt_bridge_donor['atoms']),
				bonds[:,rowspec.index('%s_resname'%that)]==salt_bridge_acceptor['resname'],
				np.in1d(bonds[:,rowspec.index('%s_atom'%that)],salt_bridge_acceptor['atoms']),
				],axis=0)
			# another loop over OR for subject target order, hence symmetry
			for this,that in subject_target_combos
			# loop over possible salt bridges
			for salt_bridge_donor in self.defn_salt_bridges['donor']
			for salt_bridge_acceptor in self.defn_salt_bridges['acceptor']
			],axis=0)
		return salt_filter

class SaltBridge:
	"""
	Identify salt bridges.
	"""
	"""
	A note on other software definitions.
	VMD seems permissive: 
		https://www.ks.uiuc.edu/Research/vmd/plugins/saltbr/
	> A salt bridge is considered to be formed if the distance between any of 
	the oxygen atoms of acidic residues and the nitrogen atoms of basic residues 
	are within the cut-off distance (default 3.2 Angstroms) in at least one 
	frame. The default distance cut-off can be changed by the user. This plugin 
	does not attempt to identify hydrogen bonds.
	#! see also 
		http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/
	#! see also 
		https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6016756/
	> Salt-bridge and aromatic-aromatic bond calculation
	> Positively charged atoms are Lys NZ1, Arg NH1, Arg NH2 and His NE2. 
	Negatively charged atoms are Asp OD1, Asp OD2, Glu OE1, Glu OE2. A salt 
	bridge is defined if two oppositely charged atoms lie within 4A across the 
	interface. A pi-pi interaction is defined between two aromatic amino acids 
	if the distance calculated from their centroid is less than 7.5A.
	"""
	defn_salt_bridges_protein = {
		'donor':[
			{'resname':'ARG','atoms':['NH1','NH2']},
			{'resname':'LYS','atoms':['NZ']},
			# we allow HIS which is charged at neutral 
			#   but we discard NE2 which is protonated with HE2
			{'resname':'HIS','atoms':['ND1','NE2'][:1]}],
		'acceptor':[
			{'resname':'GLU','atoms':['OE1','OE2']},
			{'resname':'ASP','atoms':['OD1','OD2']},]}
	def __init__(self): pass
	def filter_salt_bridges_protein(self,bonds,rowspec):
		"""Filter a bonds list for salt bridges."""
		# search should be fast because it is composed of filters on bonds
		#! however there are many combinations. are there ways to improve?
		subject_target_combos = [['subject','target'],['target','subject']]
		salt_filter = np.any([
			np.all([
				bonds[:,rowspec.index('%s_resname'%this)]==salt_bridge_donor['resname'],
				np.in1d(bonds[:,rowspec.index('%s_atom'%this)],salt_bridge_donor['atoms']),
				bonds[:,rowspec.index('%s_resname'%that)]==salt_bridge_acceptor['resname'],
				np.in1d(bonds[:,rowspec.index('%s_atom'%that)],salt_bridge_acceptor['atoms']),
				],axis=0)
			# another loop over OR for subject target order, hence symmetry
			for this,that in subject_target_combos
			# loop over possible salt bridges
			for salt_bridge_donor in self.defn_salt_bridges_protein['donor']
			for salt_bridge_acceptor in self.defn_salt_bridges_protein['acceptor']
			],axis=0)
		return salt_filter
	def filter_salt_bridges_protein_lipid(self,bonds,rowspec,kind='permissive'):
		"""Filter a bonds list for salt bridges between a protein and a membrane."""
		if kind=='permissive':
			lipid_selection = np.in1d(bonds[:,
				rowspec.index('target_atom')].astype('<U1'),['O'])
		else: raise Exception('unclear protein_lipid bond style: %s'%kind)
		# assume target is lipid and subject is protein
		this,that = 'subject','target'
		salt_filter = np.any([
			np.all([
				bonds[:,rowspec.index('%s_resname'%this)]==salt_bridge_donor['resname'],
				np.in1d(bonds[:,rowspec.index('%s_atom'%this)],salt_bridge_donor['atoms']),
				lipid_selection,
				],axis=0)
			# loop over donors only
			for salt_bridge_donor in self.defn_salt_bridges_protein['donor']
			],axis=0)
		return salt_filter