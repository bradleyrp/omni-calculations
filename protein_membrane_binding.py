#!/usr/bin/env python

"""
Analyze hydrogen bonding and contact data to plot protein-membrane interactions.

UPSTREAM:
~~~~~~~~~
residue_codes,picturesave
WorkSpace for
	plotdir for picturesave
	plotload for data object
	simulation names from metadata
constant names in loader
get_lipid_resnames
"""

from omni import WorkSpace
control.routine = None

#from calcs.codes.salt_bridge_definitions import SaltBridge
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
	>Positively charged atoms are Lys NZ1, Arg NH1, Arg NH2 and His NE2. 
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

from ortho import importer,uniform
from calcs.codes.consts import residue_codes
from omni import picturesave,status
from omni import subdivide_trajectory,zoom_figure,PostAccumulator

from omni.base.store import picturesave
from omni.base.utils import subdivide_trajectory
from omni.plotter.utils import zoom_figure
from omni.base.utils import PostAccumulator
from omni import status

import numpy as np
import pandas as pd
from collections import OrderedDict as odict

import matplotlib as mpl
import matplotlib.pyplot as plt

@loader
def load():
	"""Runs once to provide data."""
	work = WorkSpace(analysis=True)
	#! loading more than necessary here. data.set works best
	data = work.plotload(plotname='contacts')
	sns = sns_ordered = work.sns()
	# custom parameters
	from calcs.codes.fetch_coordinates import get_lipid_resnames
	lipid_resnames = get_lipid_resnames(work)
	lipid_colors = {'DOPE':'gray','DOPS':'blue','PI2P':'red'}
	lipid_canon_order = ['PI2P','DOPS','DOPE']
	post = PostAccumulator()

def compute_protein_lipid_bonds(sn,explicit=True,trim=True,kind=None,**kwargs):
	"""
	Take the bonds and observations list and compute protein_lipid bonds.
	"""
	lipids = lipid_resnames
	if not explicit: raise Exception('dev')
	# must run data.set beforehand to select the right bonds
	dat = data.this[sn]
	bonds_this,obs_this = [dat[i] for i in ['bonds','observations']]
	rowspec_this = dat['defn_rowspec']
	counts = {}
	# loop over lipid types
	for lipid in lipids:
		lipid_this = [lipid]
		# filter peptide-lipid bonds
		# note that this operation is symmetric so either subject and target can
		#   be either a lipid or protein and it will detect them properly
		#   because it uses residue_codes and lipid_this, however we typically
		#   set the subject and target to be specifically one or the other
		#!!!!!if kind=='contacts':
		#! previous note: that we could have used bond_strict_order also
		#!   it might be better to code up a symmetry switch if desired?
		residue_lipid_filter = np.any((
			np.all((np.in1d(bonds_this[:,rowspec_this.index('subject_resname')],lipid_this),
				np.in1d(bonds_this[:,rowspec_this.index('target_resname')],
					residue_codes.keys())),axis=0),
			np.all((np.in1d(bonds_this[:,rowspec_this.index('target_resname')],lipid_this),
				np.in1d(bonds_this[:,rowspec_this.index('subject_resname')],
					residue_codes.keys())),axis=0),
			),axis=0)
		if kind=='salt_bridges' and False:
			valid_salt_bridges = kwargs['valid_salt_bridges']
			inds = np.where(
				np.all((
					#! using symmetry here is definitely not necessary. the target must be a lipid
					np.any((
						np.all((np.in1d(bonds_this[:,rowspec_this.index('subject_resname')],lipid_this),
							np.in1d(bonds_this[:,rowspec_this.index('target_resname')],
								residue_codes.keys())),axis=0),
						np.all((np.in1d(bonds_this[:,rowspec_this.index('target_resname')],lipid_this),
							np.in1d(bonds_this[:,rowspec_this.index('subject_resname')],
								residue_codes.keys())),axis=0),
						),axis=0),
					# match a lipid oxygen by the first letter using the astype trick
					# note that this uses the Oxygen naming convention
					np.in1d(bonds_this[:,rowspec_this.index('target_atom')].astype('<U1'),['O']),
					# match any candidate bond by a list of possible bonds
					np.any(np.array([
						# ensure the right subject resname and atom
						np.all((
							bonds_this[:,rowspec_this.index('subject_resname')]==salt_bridge['resname'],
							np.in1d(bonds_this[:,rowspec_this.index('subject_atom')],salt_bridge['atoms']),
							),axis=0)
						for salt_bridge in valid_salt_bridges
						]),axis=0),
					),axis=0))
		if kind=='salt_bridges':
			salt_filter = SaltBridge().filter_salt_bridges_protein_lipid(
				rowspec=rowspec_this,bonds=bonds_this)
			inds = np.where(np.all((residue_lipid_filter,salt_filter),axis=0))
		elif kind=='contacts':
			inds = np.where(residue_lipid_filter)
		#else: raise Exception('invalid kind: %s'%kind)
		counts[lipid] = obs_this.T[inds].sum(axis=0)
		# trim removes zero observation lipids
		if trim and counts[lipid].sum()==0: del counts[lipid]
	return counts

def plot_protein_lipid_timeseries(sns,kind,cutoff,**kwargs):
	"""
	Plot a segmented timeseries of the current data bond type.
	"""
	global lipid_resnames
	n_segments = 20
	global data
	# compute
	counts = post.get(**{'kind':kind,'cutoff':cutoff})
	lipids_this = [m for m in lipid_canon_order if m in 
			set([j for k in [i.keys() for i in counts.values()] for j in k])]
	# +++ TOO CUSTOM figure size
	fig,axes = plt.subplots(nrows=1,ncols=len(sns),figsize=(12,6))
	for snum,sn in enumerate(sns):
		ax = axes[snum]
		zero = np.zeros(n_segments)
		for lnum,lipid in enumerate(lipids_this):
			traj = counts[sn][lipid]
			values = [traj[subdivide_trajectory(i,n_segments,nframes=len(traj))].mean() 
				for i in range(n_segments)]
			kwargs_bar = {}
			if lipid in lipid_colors: kwargs_bar['color'] = lipid_colors[lipid]
			ax.bar(np.arange(n_segments),values,bottom=zero,width=1.0,align='edge',**kwargs_bar)
			ax.set_title(sn)
			zero += values
	# normalize plot maxima (true max, not the subdivision minima max)
	count_max = max([np.array(counts[sn].values()).sum(axis=0).max() for sn in sns_this])
	for snum,sn in enumerate(sns): 
		ax = axes[snum]
		ax.set_ylim((0,count_max))
	namer_default = lambda kind,plot='protein_lipid_binding.timeseries': (
		'fig.%s.%s'%(plot,kind))
	namer = kwargs.get('namer',namer_default)
	#! zoom_figure(fig,3.0)
	picturesave(namer(kind=kwargs.get('kind_label',kind)),work.plotdir,
		backup=False,version=True,meta={},extras=[],dpi=300,form='png')

def plot_protein_lipid_histograms(sns,kind,cutoff,**kwargs):
	"""
	Histograms of protein-llipid contacts
	"""
	global data,post
	# DIMENSION 1: simulations BY ROW
	sns_this = sns
	# unpack counts
	counts = post.get(**{'kind':kind,'cutoff':cutoff})
	if False:
		data.set('contacts',select={'cutoff':3.4,
			'object':{'predefined':'lipid heavy'},
			'subject':{'predefined':'protein heavy'}})
		# compute
		counts = dict([(sn,
			compute_protein_lipid_bonds(
				data.this[sn],lipids=lipid_resnames,explicit=True,
				kind='protein_lipid_salt',valid_salt_bridges=valid_salt_bridges)
			) for sn in sns_this])
	# DIMENSION 3: LIPIDS by COLOR inside each plot
	lipids_this = [m for m in lipid_canon_order if m in 
		set([j for k in [i.keys() for i in counts.values()] for j in k])]
	# DIMENSION 2: valid PROTEIN RESIDUES by column with an extra filter
	residues_remap = {'nwaspbilayernochl0':slice(0,22)}
	resid_resname = dict([(sn,odict(zip(data.this[sn]['subject_residues_resids'],
		data.this[sn]['subject_residues_resnames']))) for sn in sns_this])
	# remap once to apply a first filter
	resid_resname = dict([(sn,
		odict(resid_resname[sn].items()[residues_remap.get(sn,slice(None,None))])) for sn in sns_this])
	# for salt bridges we modify the table
	#valid_resnames = list(set([i['resname'] for i in valid_salt_bridges]))
	#resid_resname = dict([(sn,dict([(i,j) for i,j in resid_resname[sn].items() 
	#	if j in valid_resnames])) for sn in sns_this])
	n_residues = max([len(j) for i,j in resid_resname.items()])
	nframes = dict([(sn,len(i['observations'])) for sn,i in data.this.items()])

	def stylize(ax,label,peak_count,peak_observations):
		ax.set_xlabel(label,fontsize=18)
		#ax.set_ylim(1,peak_count+1)
		#ax.set_ylim(0-0.5,peak_count-0.5)
		#ax.set_yticks(range(peak_count))
		ax.set_xlim(0,peak_observations)
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['left'].set_visible(False)
		ax.set_xticks([])
		ax.set_xticklabels([])
		#ax.set_yticklabels(['%d'%int(i) for i in np.arange(1,1+peak_observations)])

	# multi-pass execution of this function to get normalization factors
	for what in ['raw','reduce','plot'][:]: 
		if what=='raw': counts_raw = {}
		# get the peak of the observations ()
		if what=='reduce':
			peak = max([max(counts_raw[sn][resid].get(lipid,[0])) 
				for sn in sns_this for resid in resid_resname[sn] for lipid in lipids_this])
			peak_count = int(peak)
			bins = range(1,int(peak_count)+1)
		if what=='plot': 
			# counts are already normalized by nframes here
			#peak_observations = np.concatenate([np.array(i.values()) for i in counts_raw[sn].values()]).max()
			# the maximum observation is the largest sum of salt bridges with each lipid type, which
			#   means that the sum can be larger than unity because you can have more than one bond per frame
			#peak_observations = max([np.array(counts_raw[sn][k].values()).sum(axis=0).max() for k in counts_raw[sn]])
			peak_observations = max([max([np.array(counts_raw[sn][k].values()).sum(axis=0).max() 
				for k in counts_raw[sn]]) for sn in sns_this])
			#peak_observations = int(peak)
			counts_total = dict([(sn,dict([(lipid,np.sum([counts_raw[sn][resid][lipid]
				for resid in resid_resname[sn]],axis=0)) for lipid in lipids_this]))
				for sn in sns_this])
			# sort the total to confirm that many of the bars are near the max
			#peak_observations_total = np.concatenate([i.values() for i in counts_raw[sn].values()]).sum(axis=1).max()
			#peak_observations_total = np.array([np.array(counts_raw[sn][i].values()) for i in counts_raw[sn].keys()]).sum(axis=1).max(axis=0).max()
			#peak_observations_total = max([np.array([np.array(counts_raw[sn][i].values()) for i in counts_raw[sn].keys()]).sum(axis=1).max(axis=0).max() for sn in sns_this])
			#peak_observations_total = max([np.array([np.array(counts_raw[sn][i].values()) for i in counts_raw[sn].keys()]).sum(axis=0).max() for sn in sns_this])
			peak_observations_total = max([np.array([np.array(counts_raw[sn][i].values()) for i in counts_raw[sn].keys()]).sum(axis=0).sum(axis=0).max() for sn in sns_this])
			if True:
				# +++ TOO CUSTOM figure size
				fig,axes = plt.subplots(
					nrows=len(sns_this),ncols=n_residues+1,squeeze=False,figsize=(20,8))
			else:
				from omni.legacy.panels import panelplot
				# +++ TOO CUSTOM figure size
				axes,fig = panelplot(figsize=(20,8),
					layout={'out':{'hspace':0.5,'wspace':0.5,
					'grid':[1,n_residues+1]},'ins':[{'hspace':0.5,'wspace':0.5,
					'grid':[len(sns_this),1]} for i in range(n_residues+1)]})
				# panel plot is terrible. backwards. or they changed it
				axes = np.transpose(axes)

			# get resids with nonzero bonds for each type by classifying simulations
			#   and searching for nonzero bonds for each
			proteins_unique = list(set([tuple(np.concatenate(
				resid_resname[sn].items())) for sn in sns_this]))
			valid_resids = [[] for i in proteins_unique]
			sn_to_unique = {}
			for sn in sns_this:
				index_u = proteins_unique.index(tuple(np.concatenate(resid_resname[sn].items())))
				sn_to_unique[sn] = index_u
				valid_resids[index_u].extend(
					dict([(k,v) for k,v in counts_raw[sn].items() 
						if np.concatenate(v.values()).sum()>0]).keys())
			valid_resids = [sorted(set(i)) for i in valid_resids]
			resid_resname = odict([(sn,odict([(i,resid_resname[sn][i]) for i in resid_resname[sn].keys() 
				if i in valid_resids[sn_to_unique[sn]]])) for sn in sns_this])
			# remove absent columns 
			for snum,sn in enumerate(sns_this):
				for rnum in range(n_residues):
					if rnum>=len(resid_resname[sn]): fig.delaxes(axes[snum][rnum+1])
		# loop over simulations (rows)
		for snum,sn in enumerate(sns_this):
			if what=='raw': counts_raw[sn] = {}
			resids = resid_resname[sn].keys()
			resnames = data.this[sn]['subject_residues_resnames']
			for rnum,(resid,resname) in enumerate(resid_resname[sn].items()):
				# remove extra plots if fewer residues
				if what=='raw': 
					counts_raw[sn][resid] = {}
				if what=='plot':
					ax = axes[snum][rnum+1]
				# loop over lipids
				for lipid in lipids_this:
					if what in ['raw','reduce']:
						# calculation for one small tile
						rowspec_this = data.this[sn]['defn_rowspec']
						bonds = data.this[sn]['bonds']
						obs = data.this[sn]['observations']
						subject_is_resid = bonds[:,rowspec_this.index('subject_resid')]==str(resid)
						target_is_lipid = bonds[:,rowspec_this.index('target_resname')]==lipid
						if kind=='salt_bridges':
							salt_filter = SaltBridge().filter_salt_bridges_protein_lipid(
								rowspec=rowspec_this,bonds=bonds)
							rows = np.where(np.all((target_is_lipid,subject_is_resid,salt_filter),axis=0))
						elif kind=='contacts':
							rows = np.where(np.all((target_is_lipid,subject_is_resid),axis=0))
						else: raise Exception('invalid kind: %s'%kind)
						"""
						# selection logic
						subject_is_resid = bonds[:,rowspec_this.index('subject_resid')]==str(resid)
						subject_has_atom = np.in1d(bonds[:,rowspec_this.index('subject_atom')],
							np.concatenate([i['atoms'] 
								for i in valid_salt_bridges if i['resname']==resname]))
						target_is_lipid = bonds[:,rowspec_this.index('target_resname')]==lipid
						target_atom_is_oxygen = bonds[:,
							rowspec_this.index('target_atom')].astype('<U1')=='O'
						rows = np.where(np.all((
							subject_is_resid,
							subject_has_atom,
							target_is_lipid,
							target_atom_is_oxygen,
							),axis=0))
						"""
					if what=='raw': 
						counts_raw[sn][resid][lipid] = obs.T[rows].sum(axis=0)
					elif what=='reduce':
						counts,vals = np.histogram(counts_raw[sn][resid][lipid],bins=bins)
						# rescale by nframes
						#!!! explain
						counts_raw[sn][resid][lipid] = counts/float(nframes[sn])
				if what=='plot':
					#! do not plot df[1:] because we already account for bins
					df = pd.DataFrame(counts_raw[sn][resid],columns=lipids_this)
					kwargs_bars = {}
					if lipid_colors:
						kwargs_bars['color'] = [lipid_colors[lipid] for lipid in lipids_this]
					df[:].plot.barh(ax=ax,stacked=True,legend=False,width=1.0,**kwargs_bars)
					stylize(ax,'%s%d'%(residue_codes[resname],resid),peak_count,peak_observations)
			# summary column on the left
			if what=='plot':
				ax = axes[snum][0]
				df = pd.DataFrame(counts_total[sn],columns=lipids_this)
				kwargs_bars = {}
				if lipid_colors:
					kwargs_bars['color'] = [lipid_colors[lipid] for lipid in lipids_this]
				df[:].plot.barh(ax=ax,stacked=True,legend=False,width=1.0,**kwargs_bars)
				#! the real one sums more
				#! peak_observations_total = np.concatenate(counts_total[sn].values()).max()
				stylize(ax,'All',peak_count,peak_observations_total)
				ax.set_ylabel(sn)
		if what=='plot':
			# hook: namer(kind,plot='protein_self_contacts')
			namer_default = lambda kind,plot='protein_lipid_binding': (
				'fig.%s.%s'%(plot,kind))
			namer = kwargs.get('namer',namer_default)
			zoom_figure(fig,3.0)
			picturesave(namer(kind=kwargs.get('kind_label',kind)),work.plotdir,
				backup=False,version=True,meta={},extras=[],dpi=300,form='png')

if __name__=='__main__': 

	#! exclude the large protein
	sns_this = [i for i in sns if i not in ['nwaspbilayernochl0']]
	sns_this = sns

	# select a sweep
	surveys = []
	for cutoff in [3.4,5.0]: 
		surveys.append({'kind':'contacts','cutoff':cutoff,
			'kind_label':'contacts_%s'%cutoff})
	# +++ assume salt bridge means cutoff 3.4 here
	surveys.append({'kind':'salt_bridges',
		'cutoff':3.4,'kind_label':'salt_bridges'})

	def collect_post(**kwargs):
		"""Perform post-processing."""
		global post
		cutoff = kwargs.pop('cutoff')
		kind = kwargs.pop('kind')
		if kwargs: status('PostAccumulator ignores: %s'%kwargs)
		# +++ unassume: we index post by kind,cutoff without salt bridge cutoff
		this = {'kind':kind,'cutoff':cutoff}
		# no redundant calculation
		if post.done(**this): return
		# kernel of the postprocess
		data.set('contacts',select={'cutoff':cutoff,
			'target':{'predefined':'lipid heavy'},
			'subject':{'predefined':'protein heavy'}})
		# +++ assume explicit
		incoming = dict([(sn,
			compute_protein_lipid_bonds(sn=sn,explicit=True,kind=kind))
			for sn in sns_this])
		post.add(data=incoming,meta=this)
	# cross surveys with post
	for survey in surveys: collect_post(**survey)

	# plot loop
	for survey in surveys[:1]: 
		cutoff = survey['cutoff']
		kind = survey['kind']
		data.set('contacts',select={'cutoff':cutoff,
			'target':{'predefined':'lipid heavy'},
			'subject':{'predefined':'protein heavy'}})
		plot_protein_lipid_histograms(sns=sns_this,**survey)
		#plot_protein_lipid_timeseries(sns=sns_this,**survey)

