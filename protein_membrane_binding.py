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

#! move this back
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

from ortho import importer,uniform,tracebacker
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
import copy
import re

import matplotlib as mpl
mpl.use(ortho.conf.get('mpl_backend','Agg'))
import matplotlib.pyplot as plt

@loader
def load():
	"""Runs once to provide data."""
	work = WorkSpace(analysis=True)
	# we load the full data set and then use data.set to subselect
	# note that we choose the collection of simulations here, i.e. "cursor"
	data = work.plotload(plotname='contacts')
	data_hbonds = work.plotload(plotname='hydrogen_bonding')
	sns = sns_ordered = work.sns()
	# custom parameters
	from calcs.codes.fetch_coordinates import get_lipid_resnames
	lipid_resnames = get_lipid_resnames(work)
	lipid_colors = {'DOPE':'gray','DOPS':'blue','PI2P':'red'}
	lipid_canon_order = ['PI2P','DOPS','DOPE']
	post = PostAccumulator()

def compute_protein_lipid_bonds(sn,scoring_method,trim=True,kind=None,**kwargs):
	"""
	Take the bonds and observations list and compute protein_lipid bonds.
	"""
	lipids = lipid_resnames
	# must run data.set beforehand to select the right bonds
	# salt bridges uses data
	if kind=='salt_bridges': 
		dat = data.this[sn]
		rowspec_this = dat['defn_rowspec']
		#! no reduction step in the calculation so we include it here
		rowspec_reduced_cols = np.array([0,1,3,4])
	# hydrogen bonds have a different data
	elif kind=='hydrogen_bonds': 
		dat = data_hbonds.this[sn]
		#! hardcoded rowspec here because it was not in the original (it is now)
		rowspec_this = ['subject_resname','subject_resid','subject_atom',
			'target_resname','target_resid','target_atom','hydrogen_atom']
		#! no reduction step in the calculation so we include it here
		rowspec_reduced_cols = np.array([0,1,3,4])
	else: raise Exception('invalid kind: %s'%kind)
	bonds_this,obs_this = [dat[i] for i in ['bonds','observations']]
	counts = {}
	# loop over lipid types
	for lipid in lipids:
		lipid_this = [lipid]
		# filter peptide-lipid bonds
		# note that this operation is symmetric so either subject and target can
		#   be either a lipid or protein and it will detect them properly
		#   because it uses residue_codes and lipid_this, however we typically
		#   set the subject and target to be specifically one or the other
		#! note beware python 3 needs list to wrap keys because in1d ignores
		#!   a dict_keys when a list works, which is a possible bug
		residue_lipid_filter = np.any((
			np.all((np.in1d(bonds_this[:,rowspec_this.index('subject_resname')],lipid_this),
				np.in1d(bonds_this[:,rowspec_this.index('target_resname')],
					list(residue_codes.keys()))),axis=0),
			np.all((np.in1d(bonds_this[:,rowspec_this.index('target_resname')],lipid_this),
				np.in1d(bonds_this[:,rowspec_this.index('subject_resname')],
					list(residue_codes.keys()))),axis=0),
			),axis=0)
		if kind=='salt_bridges':
			salt_filter = SaltBridge().filter_salt_bridges_protein_lipid(
				rowspec=rowspec_this,bonds=bonds_this)
			inds = np.where(np.all((residue_lipid_filter,salt_filter),axis=0))
		elif kind=='contacts':
			inds = np.where(residue_lipid_filter)
		else: 
			# identical to contacts
			inds = np.where(residue_lipid_filter)
		counts[lipid] = obs_this.T[inds].sum(axis=0)
		# trim removes zero observation lipids
		if trim and counts[lipid].sum()==0: del counts[lipid]
		elif scoring_method=='reduced':
			# +++ reduce step to map from explicit atom-atom bonds to residues
			bonds_residues,lookups = np.unique(
				bonds_this[:,rowspec_reduced_cols].astype('S'),
				axis=0,return_inverse=True)
			reduced = np.zeros((len(bonds_residues),len(obs_this))) 
			# note that forgetting inds caused reduced counts to be even higher 
			#   than explicit, which allowed me to catch the bug and add inds
			for ii,i in enumerate(lookups[inds]): reduced[i] += obs_this.T[inds][ii] 
			reduced = (reduced>0).astype(int).T
			counts[lipid] = reduced.T.sum(axis=0)
		elif scoring_method=='explicit': pass
		else: raise Exception('invalid scoring method: %s'%scoring_method)
	return counts

def plot_protein_lipid_timeseries(sns,kind,scoring_method,**kwargs):
	"""
	Plot a segmented timeseries of the current data bond type.
	"""
	cutoff = kwargs.get('cutoff',None)
	global lipid_resnames
	n_segments = kwargs.get('n_segments',20)
	global data
	figplace = kwargs.get('figplace',None)
	def get_ax(sn):
		if not figplace: 
			ax = axes[sns.index(sn)]
		else: 
			row_n,col_n = figplace(sn)
			ax = axes[row_n][col_n]
		return ax
	def stylize(ax):
		ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)
		ax.tick_params(axis='x',which='both',top=False,bottom=False,labelbottom=True)
	figure_zoom = kwargs.get('figure_zoom',3.0)
	# compute
	counts = post.get(**{'kind':kind,'cutoff':cutoff,
		'scoring_method':scoring_method,'merged':kwargs.get('merged',False)})
	lipids_this = [m for m in lipid_canon_order if m in 
			set([j for k in [i.keys() for i in counts.values()] for j in k])]
	subplots_kwargs = {}
	subplots_adjust = kwargs.get('subplots_adjust',{})
	figsize = kwargs.get('figsize',None)
	if figsize: subplots_kwargs['figsize'] = figsize
	fig,axes = plt.subplots(
		nrows=kwargs.get('nrows',1),ncols=kwargs.get('ncols',len(sns)),
		**subplots_kwargs)
	if subplots_adjust: plt.subplots_adjust(**subplots_adjust)
	if not figplace: axes = [i for j in axes for i in j]
	for snum,sn in enumerate(sns):
		ax = get_ax(sn)
		duration = data.extras[sn]['end']/1000.-data.extras[sn]['start']/1000.
		zero = np.zeros(n_segments)
		for lnum,lipid in enumerate(lipids_this):
			if lipid not in counts[sn]: continue
			traj = counts[sn][lipid]
			values = [traj[subdivide_trajectory(i,n_segments,nframes=len(traj))].mean() 
				for i in range(n_segments)]
			kwargs_bar = {}
			if lipid in lipid_colors: kwargs_bar['color'] = lipid_colors[lipid]
			bar_width = duration/n_segments
			xvals = np.arange(0,duration,bar_width)
			ax.bar(xvals,values,bottom=zero,width=bar_width,align='edge',**kwargs_bar)
			ax.set_title(work.metadata.meta.get(sn,{}).get('label_compact',sn))
			ax.set_ylabel('bonds')
			ax.set_xlabel('time')
			zero += values
		ax.set_xlim((0,duration))
		stylize(ax)
	# normalize plot maxima (true max, not the subdivision minima max)
	count_max = max([np.array(list(counts[sn].values())).sum(axis=0).max() for sn in sns])
	for snum,sn in enumerate(sns): 
		ax = get_ax(sn)
		ax.set_ylim((0,count_max))
	namer_default = lambda kind,plot='protein_lipid_binding.timeseries': (
		'fig.%s.%s'%(plot,kind))
	namer = kwargs.get('namer',namer_default)
	if figure_zoom: zoom_figure(fig,figure_zoom)
	collection_name = kwargs.get('collection_name',None)
	collection_tag = '' if not collection_name else '_%s'%collection_name
	if scoring_method=='explicit': scoring_method_tag = '.explicit'
	else: scoring_method_tag = ''
	kind_tags = kwargs.get('kind_label',kind)+scoring_method_tag
	picturesave(namer(plot='binding_series'+collection_tag,kind=kind_tags),work.plotdir,
		backup=False,version=True,meta={'sns':sns_this},extras=[],dpi=300,form='png')

def plot_protein_lipid_histograms(sns,kind,scoring_method,**kwargs):
	"""
	Histograms of protein-llipid contacts
	"""
	global counts_raw,df_out # debugging
	global data,post
	scale_counts = kwargs.get('scale_counts',False)
	cutoff = kwargs.get('cutoff',None)
	figsize = kwargs.get('figsize',None)
	figure_zoom = kwargs.get('figure_zoom',3.0)
	# DIMENSION 1: simulations BY ROW
	sns_this = sns
	if kwargs.get('merged',False): sns_raw = work.sns()
	else: sns_raw = sns
	# unpack counts
	#! 	 fetch = {'kind':kind,'scoring_method':scoring_method}
	#!   if cutoff: fetch['cutoff'] = cutoff
	#!   counts = post.get(**fetch)
	counts = post.get(**{'kind':kind,'cutoff':cutoff,'scoring_method':scoring_method})
	# DIMENSION 3: LIPIDS by COLOR inside each plot
	lipids_this = [m for m in lipid_canon_order if m in 
		set([j for k in [i.keys() for i in counts.values()] for j in k])]
	# DIMENSION 2: valid PROTEIN RESIDUES by column with an extra filter
	residues_remap = {'nwaspbilayernochl0':slice(0,22)}
	resid_resname = dict([(sn,odict(zip(data.this[sn]['subject_residues_resids'].astype(int),
		data.this[sn]['subject_residues_resnames'].astype(str)))) for sn in sns_raw])
	# remap once to apply a first filter
	resid_resname = dict([(sn,
		odict(list(resid_resname[sn].items())[residues_remap.get(
			sn,slice(None,None))])) for sn in sns_raw])
	if kwargs.get('merged',False):
		merge_rule = kwargs['merge_rule']
		resid_resname_merged = {} 
		for sn in merge_rule:
			them = []
			for sn_this in merge_rule[sn]:
				them.append(resid_resname[sn_this])
			if not all([i==them[0] for i in them]):
				raise Exception('inconsistency in merge rules')
			resid_resname_merged[sn] = them[0]
		resid_resname = resid_resname_merged
		# prepare a disambiguator for later
		# +++ assume merge rules are for identical simulations
		sn_disambiguator = lambda sn:merge_rule[sn][0]
		# +++ assume the merge rule is a concatenation (i.e. catting the 
		#   trajectories) and not merging them, since they could have different
		#   numbers of frames anyway
		nframes = dict([(sn,sum([len(data.this[i]['observations']) 
			for i in merge_rule[sn]])) for sn in merge_rule])
	else: 
		sn_disambiguator = lambda sn:sn
		nframes = dict([(sn,len(i['observations'])) for sn,i in data.this.items()])

	n_residues = max([len(j) for i,j in resid_resname.items()])

	def stylize(ax,label,peak_count,peak_observations):
		ax.set_xlabel(label)
		#! ax.set_ylim(1,peak_count+1)
		#! ax.set_ylim(0-0.5,peak_count-0.5)
		#! ax.set_yticks(range(peak_count))
		ax.set_xlim(0,peak_observations)
		ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)
		ax.tick_params(axis='x',which='both',top=False,bottom=False,labelbottom=True)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['left'].set_visible(False)
		ax.set_xticks([])
		ax.set_xticklabels([])
		#! ax.set_yticklabels(['%d'%int(i) for i in np.arange(1,1+peak_observations)])

	# multi-pass execution of this function to get normalization factors
	for what in ['raw','reduce','plot']: 
		if what=='raw': counts_raw = {}
		# get the peak of the observations ()
		if what=='reduce':
			peak = max([max(counts_raw[sn][resid].get(lipid,[0])) 
				for sn in sns_this for resid in resid_resname[sn] for lipid in lipids_this])
			peak_count = int(peak)
			bins = range(0,int(peak_count)+1)
		if what=='plot': 
			# counts are already normalized by nframes here
			# note that we exclude the zero bin, which is otherwise included
			peak_observations = max([max([np.array(list(
				counts_raw[sn][k].values()))[:,1:].sum(axis=0).max() 
				for k in counts_raw[sn]]) for sn in sns_this]) 
			counts_total = dict([(sn,dict([(lipid,np.sum([counts_raw[sn][resid][lipid]
				for resid in resid_resname[sn]],axis=0)) for lipid in lipids_this]))
				for sn in sns_this])
			peak_observations_total = max([np.array([np.array(list(counts_raw[sn][i].values()))[:,1:] 
				for i in list(counts_raw[sn].keys())]).sum(axis=0).sum(axis=0).max() for sn in sns_this])
			# the residues with nonzero counts are displayed only
			resids_consensus = np.unique(np.concatenate([[k for k in counts_raw[sn] 
				if np.array(list(counts_raw[sn][k].values()))[:,1:].sum()>0] for sn in sns_this],-1))
			n_residues = len(resids_consensus)
			nrows,ncols =len(sns_this),n_residues+1
			subplots_kwargs = {}
			if figsize: subplots_kwargs['figsize'] = figsize
			fig,axes = plt.subplots(
				nrows=nrows,ncols=ncols,squeeze=False,**subplots_kwargs)
			# get resids with nonzero bonds for each type by classifying simulations
			#   and searching for nonzero bonds for each
			proteins_unique = list(set([tuple(np.concatenate(
				list(resid_resname[sn].items()))) for sn in sns_this]))
			valid_resids = [[] for i in proteins_unique]
			sn_to_unique = {}
			for sn in sns_this:
				index_u = proteins_unique.index(tuple(np.concatenate(list(resid_resname[sn].items()))))
				sn_to_unique[sn] = index_u
				valid_resids[index_u].extend(
					dict([(k,v) for k,v in counts_raw[sn].items() 
						if np.concatenate(list(v.values())).sum()>0]).keys())
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
			# reduce residues to consensus here
			if what=='plot':
				resid_resname_loop = dict([(resid,resname) 
					for resid,resname in resid_resname[sn].items() if resid in resids_consensus])
			else: resid_resname_loop = resid_resname[sn]
			sn_this = sn_disambiguator(sn)
			#! python 3 we need to cast this?
			resnames = data.this[sn_this]['subject_residues_resnames'].astype(str)
			for rnum,(resid,resname) in enumerate(resid_resname_loop.items()):
				# remove extra plots if fewer residues
				if what=='raw': 
					counts_raw[sn][resid] = {}
				if what=='plot':
					ax = axes[snum][rnum+1]
				# loop over lipids
				for lipid in lipids_this:
					if what in ['raw','reduce']:
						# calculation for one small tile
		
						rowspec_this = data.this[sn_this]['defn_rowspec']
						bonds = data.this[sn_this]['bonds']
						obs = data.this[sn_this]['observations']
						subject_is_resid = bonds[:,rowspec_this.index('subject_resid')]==str(resid)
						target_is_lipid = bonds[:,rowspec_this.index('target_resname')]==lipid
						if kind=='salt_bridges':
							salt_filter = SaltBridge().filter_salt_bridges_protein_lipid(
								rowspec=rowspec_this,bonds=bonds)
							rows = np.where(np.all((target_is_lipid,subject_is_resid,salt_filter),axis=0))
						elif kind=='contacts':
							rows = np.where(np.all((target_is_lipid,subject_is_resid),axis=0))
						# same calculation as contacts
						elif kind=='hydrogen_bonds':
							rows = np.where(np.all((target_is_lipid,subject_is_resid),axis=0))
						else: raise Exception('invalid kind: %s'%kind)
					if what=='raw': 
						if scoring_method=='reduced':
							bonds_this = bonds[rows]
							obs_this = obs.T[rows].T
							if len(rows[0])==0:
								counts_raw[sn][resid][lipid] = np.zeros((len(obs)))
							else:
								# +++ reduce step to map from explicit atom-atom bonds to residues
								bonds_residues,lookups = np.unique(
									bonds_this[:,np.array([0,1,3,4])].astype('S'),axis=0,return_inverse=True)
								reduced = np.zeros((len(bonds_residues),len(obs_this))) 
								for ii,i in enumerate(lookups): 
									reduced[i] += obs_this.T[ii] 
								reduced = (reduced>0).astype(int).T
								counts_raw[sn][resid][lipid] = reduced.sum(axis=1)
						elif scoring_method=='explicit':
							counts_raw[sn][resid][lipid] = obs.T[rows].sum(axis=0)
						else: raise Exception('invalid scoring method')
					elif what=='reduce':
						counts,vals = np.histogram(counts_raw[sn][resid][lipid],bins=bins)
						# values run from 0 to one larger than the max
						if scale_counts:
							counts = counts*vals[:-1]
						# rescale by nframes so that we do not have a bias in the histogram
						#   because the total number of frames is different
						counts_raw[sn][resid][lipid] = counts/float(nframes[sn])
				if what=='plot':
					#! do not plot df[1:] because we already account for bins
					df = pd.DataFrame(counts_raw[sn][resid],columns=lipids_this)
					kwargs_bars = {}
					if lipid_colors:
						kwargs_bars['color'] = [lipid_colors[lipid] for lipid in lipids_this]
					df[1:].plot.barh(ax=ax,stacked=True,legend=False,width=1.0,**kwargs_bars)
					df_out = df
					stylize(ax,'%s%d'%(residue_codes[resname],resid),peak_count,peak_observations)
			# summary column on the left
			if what=='plot':
				ax = axes[snum][0]
				df = pd.DataFrame(counts_total[sn],columns=lipids_this)
				kwargs_bars = {}
				if lipid_colors:
					kwargs_bars['color'] = [lipid_colors[lipid] for lipid in lipids_this]
				df[1:].plot.barh(ax=ax,stacked=True,legend=False,width=1.0,**kwargs_bars)
				#! the real one sums more
				#! peak_observations_total = np.concatenate(counts_total[sn].values()).max()
				stylize(ax,'All',peak_count,peak_observations_total)
				ax.set_ylabel(work.metadata.meta.get(sn,{}).get('label_compact',sn))
		if what=='plot':
			# hook: namer(kind,plot='protein_self_contacts')
			namer_default = lambda kind,plot='protein_lipid_binding': (
				'fig.%s.%s'%(plot,kind))
			if figure_zoom: zoom_figure(fig,figure_zoom)
			#! naming sequence should be another function (repeated elsewhere)
			namer = kwargs.get('namer',namer_default)
			# +++ assume reduced otherwise note that the scoring method is explicit
			if scoring_method=='explicit': scoring_method_tag = '.explicit'
			else: scoring_method_tag = ''
			if scale_counts==False: scale_counts_tag = '.unscaled'
			else: scale_counts_tag = ''
			kind_tags = kwargs.get('kind_label',kind)+scoring_method_tag+scale_counts_tag
			collection_name = kwargs.get('collection_name',None)
			collection_tag = '' if not collection_name else '_%s'%collection_name
			# metadata are either in the filename (via tags) or in the meta below
			if not kwargs.get('merged',False): meta = {'sns':sns_this}
			#! deferring a decision on metadata for merged for now
			else: meta = {}
			picturesave(namer(plot='binding_distn'+collection_tag,kind=kind_tags),work.plotdir,
				backup=False,version=True,meta=meta,
				extras=[],dpi=300,form='png')

if __name__=='__main__': 

	"""
	Note two styles for plot scripts: ones with @function decorators are 
	useful for running a specific plot with one parameter and no arguments.
	This script uses main to loop over many.
	"""

	### COLLECT DATA

	# subselect simulation names
	sns_this = sns

	# +++ select subplots with different simulations
	collection_names = ['mdia2','gelsolin','nwasp']
	sns_groups = dict([(c,work.metadata.collections[c]) 
		for c in collection_names])
	# loop over sets of simulations
	loop_sns = list([('all',sns)] if not sns_groups else sns_groups.items())

	# select a sweep
	surveys = []
	# select cutoffs here (note that loading the 5.0 requires more memory)
	for cutoff in [3.4,5.0][:1]: 
		for scoring_method in ['explicit','reduced']:
			surveys.append({
				'kind':'contacts','cutoff':cutoff,
				'kind_label':'contacts_%s'%cutoff,
				'scoring_method':scoring_method,
				'scale_counts':True})

	surveys = []
	for scoring_method in ['explicit','reduced']:
		# +++ assume salt bridge means cutoff 3.4 here
		surveys.append({
			'kind':'salt_bridges',
			'cutoff':3.4,
			'kind_label':'salt_bridges',
			'scoring_method':scoring_method,
			'scale_counts':True})
		surveys.append({
			'kind':'hydrogen_bonds',
			'kind_label':'hydrogen_bonds',
			'scoring_method':scoring_method,
			'scale_counts':True})

	def collect_post(**kwargs):
		"""Perform post-processing."""
		global post
		cutoff = kwargs.pop('cutoff',None)
		kind = kwargs.pop('kind')
		scoring_method = kwargs.pop('scoring_method')
		if kwargs: status('PostAccumulator ignores: %s'%kwargs)
		# +++ unassume: we index post by kind,cutoff without salt bridge cutoff
		this = {'kind':kind,'scoring_method':scoring_method}
		if cutoff:
			this['cutoff'] = cutoff
		# no redundant calculation
		if post.done(**this): return
		# kernel of the postprocess
		if kind in ['contacts','salt_bridges']:
			data.set('contacts',select={'cutoff':cutoff,
				'target':{'predefined':'lipid heavy'},
				'subject':{'predefined':'protein heavy'},})
		incoming = dict([(sn,
			compute_protein_lipid_bonds(
				sn=sn,scoring_method=scoring_method,kind=kind))
			for sn in sns_this])
		post.add(data=incoming,meta=this)

	# cross surveys with post and save the data (only once)
	for survey in surveys: collect_post(**survey)

	### MERGE REPLICATES

	merge_rule = {
		'gelsolin_nochl':['gelbilayer_nochl','gelbilayer_nochl3'],
		'gelsolin_10':['gelbilayer10','gelbilayer10_2'],
		'gelsolin_20':['gelbilayerphys','gelbilayerphys2'],
		'gelsolin_30':['gelbilayer30','gelbilayer30_2'],
		'gelsolin_mut':['gelmutbilayer20','gelmutbilayer20_2'],

		'mdia2_nochl':['mdia2bilayer_nochl2','mdia2bilayer_nochl3'],
		'mdia2_10':['mdia2bilayer10','mdia2bilayer10_2'],
		'mdia2_20':['mdia2bilayerphys','mdia2bilayerphys2'],
		'mdia2_30':['mdia2bilayer30','mdia2bilayer30_2'],

		'nwasp_nochl':['nwaspbilayer_nochl','nwaspbilayer_nochl2'],
		'nwasp_10':['nwaspbilayer10','nwaspbilayer10_2'],
		'nwasp_20':['nwaspbilayerphys','nwaspbilayerphys2'],
		'nwasp_30':['nwaspbilayer30','nwaspbilayer30_2'],}

	merge_labels = odict({
		'gelsolin_nochl':r'0% CHOL'+'\n'+r'20% $\mathrm{{PIP}_2}$',
		'gelsolin_10':r'10% $\mathrm{{PIP}_2}$',
		'gelsolin_20':r'20% $\mathrm{{PIP}_2}$',
		'gelsolin_30':r'30% $\mathrm{{PIP}_2}$',
		'gelsolin_mut':r'mutant'+'\n'+r'20% $\mathrm{{PIP}_2}$',

		'mdia2_nochl':r'0% CHOL'+'\n'+r'20% $\mathrm{{PIP}_2}$',
		'mdia2_10':r'10% $\mathrm{{PIP}_2}$',
		'mdia2_20':r'20% $\mathrm{{PIP}_2}$',
		'mdia2_30':r'30% $\mathrm{{PIP}_2}$',

		'nwasp_nochl':r'0% CHOL'+'\n'+r'20% $\mathrm{{PIP}_2}$',
		'nwasp_10':r'10% $\mathrm{{PIP}_2}$',
		'nwasp_20':r'20% $\mathrm{{PIP}_2}$',
		'nwasp_30':r'30% $\mathrm{{PIP}_2}$',})

	# loop over surveys
	for target in surveys:
		# complete the metadata
		meta_new = copy.deepcopy(post.get_meta(**target))
		# change the metadata
		meta_new['merged'] = True
		#! need to use post.done but fuzzy logic there is poor
		if meta_new not in post.meta:
			this = post.get(**target)
			that = {}
			# merge replicates into a new post
			for sn_sup,sns_sub in merge_rule.items():
				that[sn_sup] = {}
				for sname in sns_sub: 
					for lipid in this[sname]:
						if lipid not in that[sn_sup]: 
							that[sn_sup][lipid] = np.array([])
						# +++ merge rule: concatenate the frame list 
						# we concatenate the frame list to accomodate histograms
						#   while a timeseries could be merged via the mean
						that[sn_sup][lipid] = np.concatenate((
							that[sn_sup][lipid],this[sname][lipid]))
			post.add(data=that,meta=meta_new)

	if False:
		#! previously		
		# one plot, now with merged data
		survey = {'merged':True,'cutoff':3.4,'kind':'salt_bridges','scoring_method':'explicit'}
		extra = {'kind_label':'salt_bridges.merged'}
		extra['merge_rule'] = merge_rule
		#! note that you cannot include extra in survey because some meta do not
		#!   have the merged value, hence they will also match properly. this design
		#!   is worth considering in more detail
		plot_protein_lipid_histograms(sns=post.get(**survey).keys(),
			**dict(survey,**extra))

	### PLOT LOOP

	# select plots
	do_plot_histograms = 0
	do_plot_timeseries = 0
	do_summary_counter = 0
	do_special = 1

	# customizations by plot style by collection
	customize = {
		'histograms':{
			'gelsolin':{
				'figsize':(26,16),'figure_zoom':1.0,},},
		'timeseries':{},}
	for name,sns_this in loop_sns:
		#! hardcoded replicates here for rows in the timeseries
		n_reps = 2
		n_rows = n_reps
		n_cols = int(len(work.metadata.collections[name])/n_reps)
		customize['timeseries'][name] = {
			'nrows':2,'ncols':n_cols,'figsize':(10,6),'figure_zoom':1.0,
			'subplots_adjust':{'hspace':0.5,'wspace':0.5},'n_segments':150,
			'figplace':lambda sn,nrows=2,ncols=n_cols:
				(sns_this.index(sn)%nrows,int(sns_this.index(sn)/nrows))}
	#! see hardcoded dimensions
	#! currently broken by gelsolin mutant
	customize['timeseries']['all'] = {
		'nrows':9,'ncols':3,'figsize':(6,30),'figure_zoom':1.0,
		'subplots_adjust':{'hspace':0.5,'wspace':0.5},'n_segments':150,
		'figplace':lambda sn,nrows=9,ncols=3:
			(sns_this.index(sn)%nrows,int(sns_this.index(sn)/nrows))}

	# adding an "all" plot. this only applies to timeseries (see ignore below)
	#! timeseries all above is currently broken
	if False:
		if 'all' not in list(zip(*loop_sns))[0]: loop_sns += [('all',sns)]

	# plot loop
	for survey in surveys: 
		for sns_name,sns_this in loop_sns:
			cutoff = survey.get('cutoff',None)
			kind = survey['kind']
			if kind in ['contacts','salt_bridges']:
				data.set('contacts',select={'cutoff':cutoff,
					'target':{'predefined':'lipid heavy'},
					'subject':{'predefined':'protein heavy'}})
			# select plots here
			if do_plot_histograms: 
				#! ignore the all plot for histograms
				if sns_name!='all': 
					extra = {}
					if survey.get('merged',False): 
						extra['merge_rule'] = merge_rule
						raise Exception('dev')
					plot_protein_lipid_histograms(sns=sns_this,
						collection_name=sns_name,**survey,
						**customize.get('histograms',{}).get(sns_name,{}),
						**extra)
			if do_plot_timeseries:
				plot_protein_lipid_timeseries(sns=sns_this,
					collection_name=sns_name,**survey,
					**customize.get('timeseries',{}).get(sns_name,{}))

	# basic summary plot for all simulations
	if do_summary_counter:

		#! loop over meta which has the merged data
		for survey in post.meta:
			kind = survey['kind']
			scoring_method = survey['scoring_method']
			cutoff = survey.get('cutoff',None)
			merged = survey.get('merged',False)

			#! hardcoding the keys for merge
			key_regexes = ['mdia2','gel','nwasp']
			key_regex_to_title = {
				'mdia2':'mDia2','gel':'gelsolin','nwasp':'N-WASP',}
			merge_rule_r = dict([(k,i) for i,j in merge_rule.items() for k in j])

			# get counts
			counts = post.get(**{'kind':kind,'cutoff':cutoff,'scoring_method':scoring_method})
			lipids = list(np.unique(np.concatenate([list(j.keys()) 
				for i,j in counts.items() if j.keys()])))
			lipids_missing = [i for i in lipids if i not in lipid_canon_order]
			if any(lipids_missing): 
				raise Exception('missing some lipids in lipids_canon_order: %s'%lipids_missing)
			lipids = [l for l in lipid_canon_order if l in lipids]
			if not merged:
				raw = dict([(k,np.zeros((len(sns),len(lipids)))) for k in ['mean','std']])
				for snum,sn in enumerate(sns):
					for lnum,lipid in enumerate(lipids):
						if sn in counts and lipid in counts[sn]:
							#! error bars?
							raw['mean'][snum,lnum] = counts[sn][lipid].mean()
							raw['std'][snum,lnum] = counts[sn][lipid].std()
				#! is there a better way to use dataframe for both mean and std?
				df_mean = pd.DataFrame(raw['mean'],columns=lipids,index=sns)
				df_std = pd.DataFrame(raw['std'],columns=lipids,index=sns)
			elif merged:
				raw = [[[] for l in lipids] for s in merge_rule]
				for snum,sn in enumerate(merge_rule):
					for lnum,lipid in enumerate(lipids):
						for sn_this in merge_rule[sn]:
							if sn_this in counts and lipid in counts[sn_this]:
								#! error bars?
								raw[snum][lnum].append(counts[sn_this][lipid])
				# since the merged data sums the mean of the counts, we normalize by the number of sns
				counters = np.array([len(j) for i,j in merge_rule.items()])
				raw_proc = {}
				catted = [[np.concatenate(i) if len(i)>0 else [] for i in j] for j in raw]
				raw_proc['mean'] = np.array([[np.mean(i) if len(i)>0 else 0. 
					for i in j] for j in catted])
				raw_proc['std'] = np.array([[np.std(i) if len(i)>0 else 0. 
					for i in j] for j in catted])
				df_mean = pd.DataFrame(raw_proc['mean'],columns=lipids,index=merge_rule.keys())
				df_std = pd.DataFrame(raw_proc['std'],columns=lipids,index=merge_rule.keys())
				raw = raw_proc

			fig,axes = plt.subplots(nrows=1,ncols=len(key_regexes))
			bar_w = 0.8 if merged else 1.0
			for axnum,(ax,key_regex) in enumerate(zip(axes,key_regexes)):
				indices = [i for i in df_mean.index if re.match('^%s'%key_regex,i)]
				ax.set_title(key_regex_to_title[key_regex])
				df_this = df_mean.loc[indices]
				df_this_std = df_std.loc[indices]
				df_this.plot(kind='bar',stacked=True,
					# very neat magic with a dataframe as the error bar here
					yerr=df_this_std,
					color=[lipid_colors[l] for l in lipids],
					width=bar_w,ax=ax,legend=axnum==len(axes)-1)
				if axnum==len(axes)-1:
					labels_legend = [work.metadata.variables['names']['short'][i] for i in lipids]
					legend = ax.legend(labels=labels_legend,loc='upper left',
						bbox_to_anchor=(1.0,0.0,1.,1.))
					frame = legend.get_frame()
					frame.set_edgecolor('white')
					frame.set_facecolor('white')
				#! ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)
				#! ax.tick_params(axis='x',which='both',top=False,bottom=False,labelbottom=True)
				ax.spines['top'].set_visible(False)
				ax.spines['right'].set_visible(False)
				#! ax.spines['bottom'].set_visible(False)
				if axnum>0: ax.spines['left'].set_visible(False)
				if axnum==0:
					ax.set_ylabel({'salt_bridges':'salt bridges',
						'hydrogen_bonds':'hydrogen bonds'}[kind])
				else: ax.set_yticks([])
				#! note that the error bars may be clipped
				ax.set_ylim((0,raw['mean'].sum(axis=1).max()))
				if merged:
					ax.set_xticklabels([merge_labels[i] for i in indices])
				# non-merged plots use merged tick labels for simplicity
				else:
					#! assume that simulations are identical paired
					ax.set_xticks(np.arange(0,len(indices),2)+0.5)
					labels = [merge_labels[merge_rule_r[sn_this]] 
						for sn_this in indices[::2]]
					ax.set_xticklabels(labels)
				# grouping with gaps
				if not merged:
					eb_jitter = 0.1
					h,l = ax.get_legend_handles_labels()
					# move all stacked bars together
					for tnum,this in enumerate(h):
						# error bars are also indexed this way
						error_bars = ax.collections[tnum]
						error_bars.set_linewidth(0.5)
						gap_w = 0.2
						#! hardcoded where
						gap_where = range(2,8+2,2)
						gaps = dict([(i,gap_w) 
							for i in np.arange(len(this)) if i in gap_where])
						x_shift_sum = 0.
						paths = error_bars.get_paths()
						for rnum,rect in enumerate(this.patches):
							x_move = gaps.get(rnum,0.)
							rect.set_x(rect.get_x()+x_shift_sum+x_move)
							# we index the jitter on the stacked bars i.e. h above
							this_jitter = -eb_jitter*len(h)/2+(0.5+tnum)*eb_jitter
							paths[rnum].vertices += np.array([x_shift_sum+x_move+this_jitter,0.])
							x_shift_sum += x_move
					ax.set_xlim((-1*bar_w/2.-gap_w,len(h[0])*bar_w+x_shift_sum-bar_w/2.+gap_w))
					ax.set_xticks(ax.get_xticks()+np.arange(len(ax.get_xticks()))*gap_w)
				# moving error bars on merged plots
				elif merged:
					for tnum,error_bars in enumerate(ax.collections):
						for pnum,paths in enumerate(error_bars.get_paths()):
							this_jitter = -eb_jitter*len(ax.collections)/2+(0.5+tnum)*eb_jitter
							paths.vertices += np.array([this_jitter,0.])

			if False:
				df_all = [df_this]
				# bar grouping happens here
				n_df = 1
				n_col = len(dfall[0].columns) 
				n_ind = len(dfall[0].index)
				h,l = ax.get_legend_handles_labels()
				if False:
					for i in range(0,n_df*n_col,n_col): # len(h) = n_col * n_df
						for j,pa in enumerate(h[i:i+n_col]):
							for rect in pa.patches: # for each index
								rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
								rect.set_hatch(H * int(i / n_col)) #edited part     
								rect.set_width(1 / float(n_df + 1))
					axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
					axe.set_xticklabels(df.index, rotation = 0)
					axe.set_title(title)


			#! naming sequence is perhaps clumsy
			namer = namer_default = lambda kind,plot: ('fig.%s.%s'%(plot,kind))
			# +++ assume reduced otherwise note that the scoring method is explicit
			if scoring_method=='explicit': scoring_method_tag = '.explicit'
			else: scoring_method_tag = ''
			kind_tags = kwargs.get('kind_label',kind)+scoring_method_tag
			collection_name = kwargs.get('collection_name',None)
			collection_tag = '' if not collection_name else '_%s'%collection_name
			merged_tag = '' if not merged else '.merged'
			picturesave(namer(plot='binding_summary',kind=kind_tags+merged_tag),
				work.plotdir,backup=False,version=True,meta={'sns':sns_this},
				extras=[legend],dpi=300,form='png')

	if do_special: 

		#! pick one for now
		survey = {'kind': 'salt_bridges', 'scoring_method': 
		'reduced', 'cutoff': 3.4, 'merged': True}

		kind = survey['kind']
		scoring_method = survey['scoring_method']
		cutoff = survey.get('cutoff',None)
		merged = survey.get('merged',False)
		if not merged: raise Exception('grouping for this plot requires merged data')

		#! hardcoding the keys for merge
		key_regexes = ['mdia2','gel','nwasp']
		key_regex_to_title = {
			'mdia2':'mDia2','gel':'gelsolin','nwasp':'N-WASP',}
		merge_rule_r = dict([(k,i) for i,j in merge_rule.items() for k in j])

		# get counts
		counts = post.get(**{'kind':kind,'cutoff':cutoff,'scoring_method':scoring_method})
		lipids = list(np.unique(np.concatenate([list(j.keys()) 
			for i,j in counts.items() if j.keys()])))
		lipids_missing = [i for i in lipids if i not in lipid_canon_order]
		if any(lipids_missing): 
			raise Exception('missing some lipids in lipids_canon_order: %s'%lipids_missing)
		lipids = [l for l in lipid_canon_order if l in lipids]
		if not merged:
			raw = dict([(k,np.zeros((len(sns),len(lipids)))) for k in ['mean','std']])
			for snum,sn in enumerate(sns):
				for lnum,lipid in enumerate(lipids):
					if sn in counts and lipid in counts[sn]:
						#! error bars?
						raw['mean'][snum,lnum] = counts[sn][lipid].mean()
						raw['std'][snum,lnum] = counts[sn][lipid].std()
			#! is there a better way to use dataframe for both mean and std?
			df_mean = pd.DataFrame(raw['mean'],columns=lipids,index=sns)
			df_std = pd.DataFrame(raw['std'],columns=lipids,index=sns)
		elif merged:
			raw = [[[] for l in lipids] for s in merge_rule]
			for snum,sn in enumerate(merge_rule):
				for lnum,lipid in enumerate(lipids):
					for sn_this in merge_rule[sn]:
						if sn_this in counts and lipid in counts[sn_this]:
							#! error bars?
							raw[snum][lnum].append(counts[sn_this][lipid])
			# since the merged data sums the mean of the counts, we normalize by the number of sns
			counters = np.array([len(j) for i,j in merge_rule.items()])
			raw_proc = {}
			catted = [[np.concatenate(i) if len(i)>0 else [] for i in j] for j in raw]
			raw_proc['mean'] = np.array([[np.mean(i) if len(i)>0 else 0. 
				for i in j] for j in catted])
			raw_proc['std'] = np.array([[np.std(i) if len(i)>0 else 0. 
				for i in j] for j in catted])
			df_mean = pd.DataFrame(raw_proc['mean'],columns=lipids,index=merge_rule.keys())
			df_std = pd.DataFrame(raw_proc['std'],columns=lipids,index=merge_rule.keys())
			raw = raw_proc

		customs = [
			{'figname':'fig.compare.concentration','figsize':(6,6),
				'indices':[
					'gelsolin_10','gelsolin_20','gelsolin_30',
					'mdia2_10','mdia2_20','mdia2_30',
					'nwasp_10','nwasp_20','nwasp_30',
					],'gap_skip':3,'tags':{
						'where':12.0,'what':['gelsolin','mDia2','N-WASP']}},
			{'figname':'fig.compare.cholesterol','figsize':(6,6),
				'indices':[
					'gelsolin_nochl','gelsolin_20',
					'mdia2_nochl','mdia2_20',
					'nwasp_nochl','nwasp_20',
					],'gap_skip':2,'tags':{
						'where':12.0,'what':['gelsolin','mDia2','N-WASP']}},
			{'figname':'fig.compare.gelsolin_mutant','figsize':(3,6),
				'indices':[
					'gelsolin_20','gelsolin_mut',
					],'gap_skip':1,'tags':{
						'where':10.0,'what':['gelsolin','gelsolin R188L']},
				'rotate_labels':True},
			][:]

		for custom in customs:

			# figure with concentration comparison
			fig,axes = plt.subplots(nrows=1,ncols=1,figsize=custom['figsize'])
			extras = []
			axnum = 0
			axes = [axes]
			ax = axes[0]
			bar_w = 0.8 if merged else 1.0

			indices = custom['indices']
			gap_skip = custom['gap_skip']
			gap_where = range(gap_skip,len(indices)+1,gap_skip)
			bar_w = 1.0

			df_this = df_mean.loc[indices]
			df_this_std = df_std.loc[indices]
			df_this.plot(kind='bar',stacked=True,
				# very neat magic with a dataframe as the error bar here
				#! yerr=df_this_std,
				color=[lipid_colors[l] for l in lipids],
				width=bar_w,ax=ax,legend=axnum==len(axes)-1)
			if axnum==len(axes)-1:
				labels_legend = [work.metadata.variables['names']['short'][i] for i in lipids]
				legend = ax.legend(labels=labels_legend,loc='upper left',
					bbox_to_anchor=(1.0,0.0,1.,1.))
				frame = legend.get_frame()
				frame.set_edgecolor('white')
				frame.set_facecolor('white')
			#! ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)
			#! ax.tick_params(axis='x',which='both',top=False,bottom=False,labelbottom=True)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			#! ax.spines['bottom'].set_visible(False)
			if axnum>0: ax.spines['left'].set_visible(False)
			if axnum==0:
				ax.set_ylabel({'salt_bridges':'salt bridges',
					'hydrogen_bonds':'hydrogen bonds'}[kind])
			else: ax.set_yticks([])
			#! note that the error bars may be clipped
			ax.set_ylim((0,raw['mean'].sum(axis=1).max()))
			if merged:
				ax.set_xticklabels([merge_labels[i] for i in indices])
			# non-merged plots use merged tick labels for simplicity
			else:
				#! assume that simulations are identical paired
				ax.set_xticks(np.arange(0,len(indices),gap_skip)+0.5)
				labels = [merge_labels[merge_rule_r[sn_this]] 
					for sn_this in indices[::gap_skip]]
				ax.set_xticklabels(labels)
			# grouping with gaps
			eb_jitter = 0.1
			h,l = ax.get_legend_handles_labels()
			# move all stacked bars together
			for tnum,this in enumerate(h):
				# error bars are also indexed this way
				#! error_bars = ax.collections[tnum]
				#! error_bars.set_linewidth(0.5)
				gap_w = 0.2
				gaps = dict([(i,gap_w) 
					for i in np.arange(len(this)) if i in gap_where])
				x_shift_sum = 0.
				#! paths = error_bars.get_paths()
				for rnum,rect in enumerate(this.patches):
					x_move = gaps.get(rnum,0.)
					rect.set_x(rect.get_x()+x_shift_sum+x_move)
					# we index the jitter on the stacked bars i.e. h above
					#! this_jitter = -eb_jitter*len(h)/2+(0.5+tnum)*eb_jitter
					#! paths[rnum].vertices += np.array([x_shift_sum+x_move+this_jitter,0.])
					x_shift_sum += x_move
			ax.set_xlim((-1*bar_w/2.-gap_w,len(h[0])*bar_w+x_shift_sum-bar_w/2.+gap_w))
			xticks = ax.get_xticks()
			#! did this real fast
			xticks_new = []
			shift = 0
			for ii,i in enumerate(xticks):
				if ii in gaps:
					shift += gaps[ii]
				xticks_new.append(i+shift)
			ax.set_xticks(np.array(xticks_new))

			tags = custom.get('tags',None)
			if tags:
				xticks = ax.get_xticks()
				xinds = range(0,len(xticks)+gap_skip,gap_skip)
				xvals = [xticks[range(xinds[i],xinds[i+1])].mean() for i in range(len(xinds)-1)]
				height = tags['where']
				tagbox = dict(facecolor='w',alpha=1.0,boxstyle="round,pad=0.3")
				for x,label in zip(xvals,tags['what']):
					tb = ax.text(x,height,label,
						bbox=tagbox,ha='center',va='top',color='k',rotation=custom.get('rotate_labels',False)*90)
					extras.append(tb)
				if False:
					text = work.meta[sn]['ptdins_label']
					tb = ax.text(width*(xcenters[snum]),(ymax-ymin)*0.02-label_drop,text,
						bbox=tagbox_ptdins,rotation=-90 if not no_vert_labels and len(text)>15 else 0,ha="center",va="bottom",
						color='k',fontsize=art['fs']['tags'])
					extras.append(tb)
			extras.append(legend)
			picturesave(custom['figname'],
				work.plotdir,backup=False,version=True,meta={'sns':sns_this},
				extras=extras,dpi=300,form='png')
