#!/usr/bin/env python

"""
Analyze hydrogen bonding and contact data to plot protein-membrane interactions.

UPSTREAM
~~~~~~~~
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

from calcs.codes.salt_bridge_definitions import SaltBridge
from ortho import importer,uniform,tracebacker
from calcs.codes.consts import residue_codes,protein_residues
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
	#! retired for now post = PostAccumulator()

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
					protein_residues)),axis=0),
			np.all((np.in1d(bonds_this[:,rowspec_this.index('target_resname')],lipid_this),
				np.in1d(bonds_this[:,rowspec_this.index('subject_resname')],
					protein_residues)),axis=0),
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
	merged = kwargs.get('merged',False)
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
	if not figplace:
		#! causing problems
		try: axes = [i for j in axes for i in j]
		except:	pass
	for snum,sn in enumerate(sns):
		ax = get_ax(sn)
		#! selecting the first simulation
		if not merged: sn_this = sn
		else: sn_this = merge_rule[sn][0]
		duration = data.extras[sn_this]['end']/1000.-data.extras[sn_this]['start']/1000.
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
	#! this is ill-advised!!!!
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
		if False:
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
			#!!! bug: the counts_total should not be the sum of these distributions
			#!!!   particularly because we are summing a number not actually calculating
			#!!!   the distinct lipids that the entire peptide is bound to, ignoring redundancy
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
						if kwargs.get('merged',False):
							"""
							MERGE SIMULATIONS HERE
							Possible pitfall:
							  What happens when you do reduced scoring and there is a bond between residue A and lipid 100 on one simulation and residue B and lipid 100 on the other? How would this appear in the score? We want these to be distinct bonds that are scored separately, but ideally we would just merge the obs and bonds so that this complicated loop is not perturbed. A better example is this: what if we have A-100 and A-100 in both simulations. This would only get scored as a single residue-lipid bond. I am dwelling on this because I would like to fix bonds and obs rather than complicate this loop any further. My initial impulse is to make a merged bond list and then collate the obs into new rows based on that. Now I am thinking it might be best to just cat the bond list and stack the obs so there are quadrants of zeros for each different simulation. This may have the same problem, however. The benefit of the latter method is that you have distinct frames with entirely distinct row numbers. The problem is that the rows get reduced at some point when we discard the atoms.
							#!!!? should we scale the rate here because we have not catted things?
							"""
							#! cannot set this yet, because it is keyed by literal sn later:
							#!   sns_this = merge_rule[sn]
							#! obs = [data.this[i]['observations'] for i in merge_rule[sn]]
							#! obs = data.this[sn_this]['observations']
							#!!! bug: we should not use data here, because data is not set to kind
							#!!!   instead it refers to the salt bridges, and not the hydrogen bonds
							obs = [data.this[i]['observations'] for i in merge_rule[sn]]
							# get the dimensions (frames by unique bonds) for each simulation
							n_rows,n_cols = np.array([i.shape for i in obs]).T
							obs_merge = np.zeros((np.sum(n_rows),np.sum(n_cols)))
							n_rows,n_cols = np.cumsum([np.concatenate([(0,),i]) 
								for i in [n_rows,n_cols]],axis=1)
							for ii,i in enumerate(n_rows[:-1]):
								obs_merge[
									n_rows[ii]:n_rows[ii+1],
									n_cols[ii]:n_cols[ii+1],
									] = obs[ii]
							#! check: np.all(obs_merge[0:1501,0:1175]==obs[0]) 
							obs = obs_merge
							bonds = np.concatenate([data.this[i]['bonds'] for i in merge_rule[sn]])
						# standard behavior
						else: 
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
				"""
				errorlog
					recomputing the all column here
				previously, a block during the what=='plot' part of  the loop, we used this:
					counts_total = dict([(sn,dict([(lipid,np.sum([counts_raw[sn][resid][lipid]
						for resid in resid_resname[sn]],axis=0)) for lipid in lipids_this]))
						for sn in sns_this])
				we have to compute directly from the bonds
				let us try to use the bonds data directly to check on this before continuing
				actually hold off on this, just use obs and obs merge and check the numbers here
				"""
				# see comment above: reduce step to map from explicit atom-atom bonds to residues
				if 0:
					counts_total_alt = {}
					for sn_this in work.sns():
						counts_total_alt[sn_this] = {}
						for lipid in lipids_this:
							print(sn_this)
							#! this whole decision point needs to be better!!!
							if kind=='hydrogen_bonds': dat = data_hbonds.this[sn_this]
							elif kind=='salt_bridges': 
								data.set('contacts',select={'cutoff':cutoff,
									'target':{'predefined':'lipid heavy'},
									'subject':{'predefined':'protein heavy'},})
								dat = data.this[sn_this]
							else: raise Exception('invalid kind: %s'%kind)
							bonds = dat['bonds']
							obs = dat['observations']
							# subselection here
							if kind=='salt_bridges':
								#! this might not be strict enough see the hydrogen_bonding stuff
								inds = np.where(bonds[:,np.array(
									rowspec_this.index('target_resname'))]==lipid)[0]
								bonds_this = bonds[inds]
								obs_this = obs[:,inds]
								if np.any(np.in1d(bonds_this[:,
									rowspec_this.index('target_resname')],protein_residues)):
									raise Exception('salt bridge inconsistency. '
										'we expect this to be one-sided')
							elif kind=='hydrogen_bonds': 
								this = np.any([np.all([np.in1d(bonds[:,
									rowspec_this.index(['target_resname','subject_resname'][i])],
										[lipid_resnames,protein_residues][j]) 
									for i,j in k],axis=0) 
								for k in [[(0,0),(1,1)],[(0,1),(1,0)]]],axis=0)
								inds = np.where(this)[0]
								bonds_this = bonds[inds]
								obs_this = obs[:,inds]
							else: raise Exception('invalid kind: %s'%kind)
							# modified for the entire peptide
							colsel = np.array([rowspec_this.index(i) 
								for i in ['target_resname','target_resid']])
							bonds_residues,lookups = np.unique(
								bonds_this[:,colsel].astype('S'),axis=0,return_inverse=True)
							reduced = np.zeros((len(bonds_residues),len(obs_this))) 
							for ii,i in enumerate(lookups): 
								reduced[i] += obs_this.T[ii] 
							reduced = (reduced>0).astype(int).T
							counts_total_alt[sn_this][lipid] = reduced.sum(axis=1)
					#! troubleshoot / debug
					if 0:
						sn = 'mdia2bilayerphys'
						results = counts_total_alt[sn_this]['PI2P']
						fr = 0
						result = results[fr] # 7
						np.unique(bonds[np.where(obs[fr])[0]][:,
							np.array([3,4])].astype(str),axis=0).shape[0] # 7
						# so far then it looks correct, but what are the unique residue-lipid bonds
						# basically there are 7 unique lipids bound to the peptide though, 
						#   so who cares about residues?
						# so there are 18 unique residue-lipid pairs
						np.unique(bonds[np.where(obs[fr])[0]][:,np.array([0,1,3,4])].astype(str),axis=0)
						them = np.unique(bonds[np.where(obs[fr])[0]][:,
							np.array([0,1,3,4])].astype(str),axis=0).shape # 18 total
						np.unique(them[:,3]).shape # 7 total
						np.unique(them[:,np.array([2,3])],axis=0) 
						# one of them is a PS though, so it should be 6
						# so all the per-lipid counts are equal and this is a BUG
						np.all(counts_total_alt[sn_this]['PI2P']==counts_total_alt[sn_this]['DOPS'])
					# after the fix above we see mean for mdia2 is 4.84 hydrogen bonds or 3.474
					# change quickly to salt to see what's up
					# actually those *were* the salt bridges at 3.4
					# which is weird I must not be calling what I thought I was calling
					#! debug again
					if 0:
						# counts_total_alt['mdia2bilayerphys2']['PI2P'].mean() # this is 22 because it includes all lipids?
						pass
						# need to do another drilldown here!
						respairs = bonds[:,np.array([0,3])].astype(str)
						# skeleton: np.where(np.any([np.all((),axis=0) for i in [1,-1]],axis=0))[0]
						#! np.in1d(bonds[:,rowspec_this.index([
						#!   'target_resname','subject_resname'][i])],
						#!   [lipid_resnames,residue_codes.keys()][j])
						# do not forget to turn dict keys into a list before in1d
						this = np.any([np.all([np.in1d(bonds[:,rowspec_this.index(['target_resname','subject_resname'][i])],[lipid_resnames,protein_residues][j]) 
							for i,j in k],axis=0) for k in [[(0,0),(1,1)],[(0,1),(1,0)]]],axis=0)
						inds = np.where(this)[0]
						# popped this above for the subselection
						# after that fix, we have 3.605 and 4.64 for each
					import ipdb;ipdb.set_trace()
				df = pd.DataFrame(counts_total[sn],columns=lipids_this)
				kwargs_bars = {}
				if lipid_colors:
					kwargs_bars['color'] = [lipid_colors[lipid] for lipid in lipids_this]
				df[1:].plot.barh(ax=ax,stacked=True,legend=False,width=1.0,**kwargs_bars)
				#! import ipdb;ipdb.set_trace()
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
			if scale_counts==True: scale_counts_tag = '.scaled'
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

#! legacy plots are being refactored below
do_legacy_plots = False

if __name__=='__main__' and do_legacy_plots: 

	"""
	Note two styles for plot scripts: ones with @function decorators are 
	useful for running a specific plot with one parameter and no arguments.
	This script uses main to loop over many images.
	"""

	do_contacts_surveys = False
	do_std_plots = True

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

	if do_contacts_surveys:
		# select cutoffs here (note that loading the 5.0 requires more memory)
		for cutoff in [3.4,5.0][:1]: 
			for scoring_method in ['explicit','reduced']:
				surveys.append({
					'kind':'contacts','cutoff':cutoff,
					'kind_label':'contacts_%s'%cutoff,
					'scoring_method':scoring_method,
					'scale_counts':False})
	if do_std_plots:
		for scoring_method in ['explicit','reduced']:
			# +++ assume salt bridge means cutoff 3.4 here
			surveys.append({
				'kind':'salt_bridges',
				'cutoff':3.4,
				'kind_label':'salt_bridges',
				'scoring_method':scoring_method,
				'scale_counts':False})
			surveys.append({
				'kind':'hydrogen_bonds',
				'kind_label':'hydrogen_bonds',
				'scoring_method':scoring_method,
				'scale_counts':False})

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
						### TRY A DIFFERENT MERGE RULE
						if False: # this is for the histograms
							if lipid not in that[sn_sup]: 
								that[sn_sup][lipid] = np.array([])
							# +++ merge rule: concatenate the frame list 
							# we concatenate the frame list to accomodate histograms
							#   while a timeseries could be merged via the mean
							that[sn_sup][lipid] = np.concatenate((
								that[sn_sup][lipid],this[sname][lipid]))
						else:
							if lipid not in that[sn_sup]: 
								that[sn_sup][lipid] = []
							that[sn_sup][lipid].append(this[sname][lipid])
			#! very hard to get this shit working correctly!!!
			for sn_sup in that:
				for lname in that[sn_sup]:
					#!!! HACK
					min_nframes = min([len(i) for i in that[sn_sup][lname]])
					if min_nframes==0: raise Exception('hack oops')
					that[sn_sup][lname] = np.mean([i[:min_nframes] for i in that[sn_sup][lname]],axis=0)
			post.add(data=that,meta=meta_new)
			#!!! want to refresh the above? set post.meta = [] and then run main
			#!!!   actually maybe not? this isn't working

	# debugging merged timeseries SEE THE ALTERNATE MERGE RULE ABOVE
	if False:
		this = post.get(**{'kind': 'hydrogen_bonds', 'scoring_method': 
			'reduced', 'merged': True})['mdia2_20']['PI2P']
		them = [post.get(**{'kind': 'hydrogen_bonds', 'scoring_method': 
			'reduced', 'merged': False})[sn]['PI2P'] 
		for sn in ['mdia2bilayerphys','mdia2bilayerphys2']]
		np.all(np.concatenate(them)==this) # true

	### PLOT LOOP

	# select plots
	do_plot_histograms = 0
	do_plot_timeseries = 0
	do_summary_counter = 0
	do_special = 0

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

	### DEMOS SERIES

	if 1:
		survey = {'merged':True,'cutoff':3.4,'kind':'salt_bridges','scoring_method':'explicit'}
		extra = {'kind_label':'salt_bridges.merged'}
		extra['merge_rule'] = merge_rule
		#! note that you cannot include extra in survey because some meta do not
		#!   have the merged value, hence they will also match properly. this design
		#!   is worth considering in more detail
		plot_protein_lipid_histograms(sns=post.get(**survey).keys(),
			**dict(survey,**extra))
	if 0:
		survey = {'merged':True,'kind':'hydrogen_bonds','scoring_method':'reduced'}
		extra = {'kind_label':'hydrogen_bonds.merged.mdia2'}
		extra['merge_rule'] = merge_rule
		#! note that you cannot include extra in survey because some meta do not
		#!   have the merged value, hence they will also match properly. this design
		#!   is worth considering in more detail
		#! clumsy?
		sns_this = [i for i in post.get(**survey).keys() if 'mdia2' in i]
		plot_protein_lipid_histograms(sns=sns_this,
			**dict(survey,**extra))
	if 0:
		surveys = []
		survey = {'merged':True,'kind':'hydrogen_bonds','scoring_method':'reduced'}
		extra = {'kind_label':'hydrogen_bonds.merged.mdia2'}
		extra['merge_rule'] = merge_rule
		#! note that you cannot include extra in survey because some meta do not
		#!   have the merged value, hence they will also match properly. this design
		#!   is worth considering in more detail
		#! clumsy?
		sns_this = [i for i in post.get(**survey).keys() if 'mdia2' in i]
		plot_protein_lipid_histograms(sns=sns_this,
			**dict(survey,**extra))
	if 0:
		extra = {'kind_label':'hydrogen_bonds.merged.mdia2','merged':True}
		extra['merge_rule'] = merge_rule
		survey = {'kind': 'hydrogen_bonds', 'scoring_method': 'reduced'}
		sns_name,sns_this = loop_sns[0]
		sns_this = ['mdia2_nochl', 'mdia2_10', 'mdia2_20', 'mdia2_30']
		plot_protein_lipid_timeseries(sns=sns_this,
			**dict(survey,**extra))

	# adding an "all" plot. this only applies to timeseries (see ignore below)
	#! timeseries all above is currently broken
	if False:
		if 'all' not in list(zip(*loop_sns))[0]: 
			loop_sns += [('all',sns)]

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

	# reduced plot prepared for BPS
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
				work.plotdir,backup=False,version=True,
				meta={'sns':sns_this,'survey':survey},
				extras=extras,dpi=300,form='png')

def prepare_span(sn):
	"""
	Prepare a combination of hydrogen bond and salt bridge data.
	Do not forget to set the correct data cursor.
	"""
	# prepare the data
	span = {}
	bonds,obs = [data.this[sn][i] for i in ['bonds','observations']]
	bonds_h,obs_h = [data_hbonds.this[sn][i] for i in ['bonds','observations']]
	rowspec_salt = data.this[sn]['defn_rowspec']
	#! hardcoded rowspec here because it was not in the original (it is now. actually no)
	rowspec_hbonds = ['subject_resname','subject_resid','subject_atom',
		'target_resname','target_resid','target_atom','hydrogen_atom']
	# assemble these data
	span['salt_bridges'] = {'bonds':bonds,'obs':obs,'rowspec':rowspec_salt}
	span['hydrogen_bonds'] = {'bonds':bonds_h,'obs':obs_h,'rowspec':rowspec_hbonds}
	return span

def span_to_protein_residues(span):
	"""Get all protein residues from a span of bond data."""
	residues = []
	for name in span:
		rowspec = span[name]['rowspec']
		for key in ['subject','target']:
			cols = [i%key for i in ['%s_resname','%s_resid']]
			residues.extend(span[name]['bonds'][:,np.array([rowspec.index(c) for c in cols])])
	residues = [(i,j) for i,j in residues if i in protein_residues]
	residues = np.unique(np.array(residues).astype(str),axis=0)
	# cheap sort
	residues = np.array(sorted(residues,key=lambda x:x[1]))
	return residues

def compute_histograms(sn,method='count_bound_lipids',protein_residue_target=None):
	"""
	Generate bond distributions.
	"""
	if isinstance(protein_residue_target,type(None)) and method!='filter_protein': 
		raise Exception('incompatible')
	# prepare the data
	span = prepare_span(sn)

	# first reduction: remove the atoms
	col_targets = ['subject_resname','subject_resid','target_resname','target_resid']
	for key in span:
		rowspec = span[key]['rowspec']
		cols = np.array([rowspec.index(i) for i in col_targets])		
		span[key]['bonds'] = span[key]['bonds'][:,cols]

	# finalize the bond type order here
	kind_order = list(span.keys())
	# prepare the stacks
	stack_bonds,stack_obs = [[span[k][j] for k in kind_order] for j in ['bonds','obs']]

	# second reduction: a filter: consider only protein-lipid bonds
	stacks = []
	for bonds_this,obs_this in zip(stack_bonds,stack_obs):
		#! dev: bonds_this,obs_this = list(zip(stack_bonds,stack_obs))[1]
		resname_cols = [col_targets.index(i) for i in ['subject_resname','target_resname']]
		# get the number of protein residues or lipids in each row
		has_protein,has_lipid = [np.sum([np.in1d(bonds_this[:,i],group)*1 
			for i in resname_cols],axis=0) for group in [protein_residues,lipid_resnames]]
		# consistency check: two items in a pair
		if np.any((has_protein+has_lipid)!=2.): raise Exception('inconsistent pairs')
		# filter here for a single protein and lipid in each pair
		idx = np.where(np.all((has_protein==1,has_lipid==1),axis=0))[0]
		stacks.append([bonds_this[idx],obs_this[:,idx]])
	stack_bonds,stack_obs = zip(*stacks)

	# note that there will be some redundancy handled by the upcoming reduction
	bonds_cat = np.concatenate(stack_bonds).astype(str)
	# formulate all unique bonds between all bond types
	bonds_u,lookups = np.unique(bonds_cat,axis=0,return_inverse=True)
	nbonds = len(bonds_u)
	# map the rows before the cat so we can stack the obs
	ind_major = np.concatenate([np.ones(len(i))*ii for ii,i in enumerate(stack_bonds)])
	ind_minor = np.concatenate([np.arange(len(i)) for ii,i in enumerate(stack_bonds)])

	# third reduction: turn bonds over multiple instances (i.e. atoms) into a single bond
	obs_new_stack = []
	for stackno,obs_this in enumerate(stack_obs):
		# use the inverse and the stack number to get the map from multiple reundant bonds
		#   in the original bond stack to unique bonds. use of lookups does the reduction
		reindex = lookups[np.where(ind_major==stackno)[0]]
		obs_new = np.zeros((len(obs_this),nbonds))
		# numpy allows for fast non-buffered additions
		#! previous method does not sum: obs_new[:,reindex] += obs_this
		#! this method is slow: for ii,i in enumerate(reindex): obs_new[:,i] += obs_this[:,ii]
		np.add.at(obs_new.T,reindex,obs_this.T)
		# reduce to binary. note that we could use the "previous method" above to skip to this
		obs_new[obs_new>0] = 1.0
		obs_new_stack.append(obs_new)

	# recall that we ensure that we have 1 protein and 1 lipid in each pair above
	#   see the "two items in a pair" consistency check and the following filter
	resname_cols = [col_targets.index(i) for i in ['subject_resname','target_resname']]
	# when to reverse so the first column is a protein
	idx_sym_r = np.in1d(bonds_u[:,resname_cols[1]],protein_residues)*1
	# flip the bond list so the columns are protein, then resname
	# be very careful when reindexing: you have to swap the rows in separate brackets
	# note that a minor error in the original code was corrected later, but there was a
	#   minor discrepancy with the original results that took effort to confirm. this was
	#   only a problem because so few records needed to be flipped in the first place
	bonds_u[np.where(idx_sym_r)] = bonds_u[np.where(idx_sym_r)][:,np.array([2,3,0,1])]

	if method=='filter_protein':

		"""
		The "count_bound_lipids" method below does a further reduction and counts lipids bound to the 
		entire peptide. in this section we drill down to one protein. note that we are still only counting
		one or more residue-lipid bonds from the same residue as a "1" in our scoring metric
		later we will relax this in the np.add.at stage
		"""

		# subselect bonds for this protein
		subsel = np.where(np.all(bonds_u[:,:2]==protein_residue_target,axis=1))[0]
		bonds_u = bonds_u[subsel]
		# filter the entire stack
		obs_prev_stack = obs_new_stack
		obs_new_stack = []
		for obs_this in obs_prev_stack:
			obs_new_stack.append(obs_this[:,subsel])

	"""
	STATUS
	first we discarded the atoms
	we then computed unique bonds across types (salt bridge, hydrogen bonds) in bonds_u
	then we filtered for protein-lipid bonds
	a reduce step mapped multiple bond entries into a single bond
	after that step ("third reduction") the observation stack matches the bonds_u list
	a new reduce step will discard the protein side and consider only lipid targets
	note that we previously had distinct bond lists by bond type
		but now we only have one bond list. this is a key difference
	the following block will be modular, and can be removed to recover the previous behavior
	in the previous behavior, two protein residues bonding to one lipid would be counted twice
	"""

	# fourth reduce step: discard the protein side
	if do_discard_protein:
		col_nums = [col_targets.index(i) for i in ['target_resname','target_resid']]
		try: bonds_lipids_u,lookups_lipids = np.unique(bonds_u[:,col_nums],axis=0,return_inverse=True)
		#! dangerous?
		except: return None
		ntargets = len(bonds_lipids_u)
		# compared to last time, we only have one bonds list
		bonds_this = bonds_u
		# iterative naming
		obs_prev_stack = obs_new_stack
		obs_new_stack = []
		# note that stack_bonds is just a repeated stack of the bonds_u
		for stackno,obs_this in enumerate(obs_prev_stack):
			#! dev: stackno = 1; obs_this = obs_new_stack[1]
			obs_new = np.zeros((len(obs_this),ntargets))
			# compared to the previous step, we do not have to subselect the lookups
			#   because we only have one bonds list
			reindex = lookups_lipids
			np.add.at(obs_new.T,reindex,obs_this.T)
			# reduce to binary
			obs_new[obs_new>0] = 1.0
			obs_new_stack.append(obs_new)
		# rename bonds_u to continue. note that it only contains the lipid part of the bond now
		bonds_u = bonds_lipids_u
		# update the col_targets to match the new bonds list
		col_targets = ['target_resname','target_resid']

	# get the intersection and disjoint sets
	def reduce_obs_tables_by_bond(left,right,names=None):
		"""
		Decompose bond types into two disjoint sets plus the intersection.
		"""
		order = ['left','right','intersect']
		shape = left.shape
		decomp = np.zeros([len(order)]+list(shape))
		decomp[order.index('left')] = (np.sum((left==1,right==0),axis=0)==2)*1
		decomp[order.index('right')] = (np.sum((left==0,right==1),axis=0)==2)*1
		decomp[order.index('intersect')] = (np.sum((left==1,right==1),axis=0)==2)*1
		result = dict(decomp=decomp,order=order)
		# ensure the order above is respected and pass along the names if they are given
		if names:
			order_named = [None,None,None]
			for ii,i in enumerate(['left','right']):
				order_named[order.index(i)] = names[ii]
			order_named[order.index('intersect')] = 'intersect'
			result['order_named'] = order_named
		return result

	# fix the observation lists if we have some invalid frames
	if len(set([len(i) for i in obs_new_stack]))!=1:
		valid_frames_salt = data.this[sn]['valid_frames']
		valid_frames_hbonds = data_hbonds.this[sn]['valid_frames']
		consensus_frames = np.intersect1d(valid_frames_hbonds,valid_frames_salt)
		#! lazy way to do this
		reindexer = dict([(name,np.array([np.where(j==i)[0][0] for i in consensus_frames]))
			for name,j in zip(['hydrogen_bonds','salt_bridges',],
				[valid_frames_hbonds,valid_frames_salt])])
		obs_new_stack = [obs_new_stack[ii][reindexer[n],:] for ii,n in enumerate(kind_order)]

	result = reduce_obs_tables_by_bond(*obs_new_stack,names=kind_order)
	decomp,decomp_order = result['decomp'],result['order_named']

	# order by lipid type
	lipid_kinds,kind_idx = np.unique(bonds_u[:,col_targets.index('target_resname')],return_inverse=True)
	kind_subsel = [np.where(kind_idx==knum)[0] for knum in range(len(lipid_kinds))]

	# use the following comments to look specifically at hbonds or salt bridges
	# looking at just hydrogen bonds
	# decomp = [decomp[1],decomp[2]]
	# looking at just the salt bridges
	# decomp = [decomp[0],decomp[2]]
	# if no, reduction, it's both

	nframes = [len(i) for i in obs_new_stack]
	if np.unique(nframes).shape[0]!=1: raise Exception('bond types have different frame counts')
	else: nframes = nframes[0]

	def histograms_by_kind(decomp):
		"""Generate histograms from trajectories of bond counts. This must happen at the very end."""
		hists_by_lipid,counts_by_lipid = [],[]
		for knum,kind in enumerate(lipid_kinds):
			counts = np.array([d[:,kind_subsel[knum]].sum(axis=1) for d in decomp]).astype(int)
			counts_by_lipid.append(counts)
		counts_max = max([c.sum(axis=0).max() for c in counts_by_lipid])
		bins = np.arange(0,counts_max+2).astype(int)	
		for knum,kind in enumerate(lipid_kinds):
			counts = counts_by_lipid[knum]
			hist = np.histogram(counts.sum(axis=0),bins=bins)[0]
			hists_by_lipid.append(hist)
		hists_by_lipid = np.array(hists_by_lipid)/nframes
		return hists_by_lipid,bins,counts_by_lipid

	global decomps
	result = {'histograms':{}}
	for name,decomp_this in decomps.items():
		decomp_out = [decomp[decomp_order.index(i)] for i in decomp_this['items']]
		hists_by_lipid,bins,counts_by_lipid = histograms_by_kind(decomp_out)

		# consistency checks for the full set
		# note that this was mostly useful when developing the code but retained for posterity
		if decomp_this['items']==['salt_bridges','hydrogen_bonds','intersect']:
			for kind_this in kind_order:
				for lipid_this in lipid_kinds:
					# consistency check 1: histogram counts sum to total frames we observe bonds
					knum = list(lipid_kinds).index(lipid_this)
					cols = np.array([decomp_order.index(i) for i in [kind_this,'intersect']])
					idx = np.where(bonds_u[:,col_targets.index('target_resname')]==lipid_kinds[knum])[0]
					# note that we ignore the first bin, which is a bond count of zero
					if not (np.histogram(counts_by_lipid[knum][cols].sum(axis=0),bins=bins)[0][1:].sum()==
						np.sum(obs_new_stack[kind_order.index(kind_this)][:,idx].sum(axis=1)>0)):
						raise Exception('failed consistency check (1): %s,%s'%(kind_this,lipid_this))
					# consistency check 2: total bonds for a type and lipid
					knum = list(lipid_kinds).index(lipid_this)
					cols = np.array([decomp_order.index(i) for i in [kind_this,'intersect']])
					bonds = counts_by_lipid[knum][cols].sum(axis=0)
					total_bonds = bonds.sum()
					obs_this = obs_new_stack[kind_order.index(kind_this)]
					idx = np.where(bonds_u[:,col_targets.index('target_resname')]==lipid_kinds[knum])[0]
					total_bonds_alt = obs_this[:,idx].sum(axis=1).sum()
					if not total_bonds==total_bonds_alt: raise Exception('failed consistency check (2)')

		# save the result
		result['histograms'][name] = dict(
			hists_by_lipid=hists_by_lipid,bins=bins,counts_by_lipid=counts_by_lipid,
			lipid_kinds=lipid_kinds)

	return result

if __name__=='__main__': 

	do_discard_protein = True
	do_scale_histograms_by_frame = True

	# globals used for plotting and calculations
	#! move this to the loader
	lipid_colors = {'PI2P':'r','DOPE':'gray','DOPS':'b'}
	lipid_stack_order = ['PI2P','DOPS','DOPE']
	#! do not forget cholesterol
	#! go fix the Landscape in automacs
	lipid_resnames = ['DOPC', 'DOPS', 'POPC', 'DOPE', 'SAPI', 'PI2P']+['CHL1']

	# identify the target data
	cutoff = 3.4
	data.set('contacts',select={'cutoff':cutoff,
		'target':{'predefined':'lipid heavy'},
		'subject':{'predefined':'protein heavy'}})

	# supergroup names must be in the metadata collections
	supergroups = ['mdia2','gelsolin','nwasp']
	#! move this to the metadata
	conditions = {
		'mdia2':['mdia2_nochl','mdia2_10','mdia2_20','mdia2_30'],
		'gelsolin':['gelsolin_nochl','gelsolin_10','gelsolin_20','gelsolin_30'],
		'nwasp':['nwasp_nochl','nwasp_10','nwasp_20','nwasp_30'],}
	#! move this to the metadata
	merge_rule = {
		# gelsolin
		'gelsolin_nochl':['gelbilayer_nochl','gelbilayer_nochl3'],
		'gelsolin_10':['gelbilayer10','gelbilayer10_2'],
		'gelsolin_20':['gelbilayerphys','gelbilayerphys2'],
		'gelsolin_30':['gelbilayer30','gelbilayer30_2'],
		'gelsolin_mut':['gelmutbilayer20','gelmutbilayer20_2'],
		# mDia2
		'mdia2_nochl':['mdia2bilayer_nochl2','mdia2bilayer_nochl3'],
		'mdia2_10':['mdia2bilayer10','mdia2bilayer10_2'],
		'mdia2_20':['mdia2bilayerphys','mdia2bilayerphys2'],
		'mdia2_30':['mdia2bilayer30','mdia2bilayer30_2'],
		# N-WASP
		'nwasp_nochl':['nwaspbilayer_nochl','nwaspbilayer_nochl2'],
		'nwasp_10':['nwaspbilayer10','nwaspbilayer10_2'],
		'nwasp_20':['nwaspbilayerphys','nwaspbilayerphys2'],
		'nwasp_30':['nwaspbilayer30','nwaspbilayer30_2'],}
	decomps_order = ['salt_x','hbonds_x','salt_bridges','hydrogen_bonds','both']
	decomps = {
		'salt_x':{'title':'salt bridges alone','items':['salt_bridges']},
		'hbonds_x':{'title':'hydrogen bonds alone','items':['hydrogen_bonds']},
		'salt_bridges':{'title':'salt bridges','items':['salt_bridges','intersect']},
		'hydrogen_bonds':{'title':'hydrogen bonds','items':['hydrogen_bonds','intersect']},
		'both':{'title':'hydrogen bond\nand salt bridge','items':
			['salt_bridges','hydrogen_bonds','intersect']},}
	decomps_titles = {'salt_x':'salt bridge\n(exclusive)',
		'hbonds_x':'hydrogen bond\n(exclusive)',
		'hydrogen_bonds':'hydrogen bonds',
		'salt_bridges':'salt bridges',
		'both':'salt bridges +\nhydrogen bonds'}

	def prepare_panels(**kwargs):
		"""
		Map simulations onto panels for a plot.
		"""
		sg = kwargs.pop('supergroup',None)
		style = kwargs.pop('style',None)
		residues = kwargs.pop('residues',None)
		decomp_style = kwargs.pop('decomp_style',None)
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)

		# vertical plot with conditions on the rows and replicats on the columns
		if style=='rows_condition_cols_replicate':
			figdim = 3
			n_conditions = len(conditions[sg])
			n_reps = [len(merge_rule[c]) for c in conditions[sg]]
			if len(set(n_reps))!=1: 
				raise Exception('development: uneven number of columns. need axis deletes')
			else: n_reps = n_reps[0]				
			nrows,ncols = n_conditions,n_reps
			fig,axes = plt.subplots(nrows=nrows,ncols=ncols,figsize=(figdim*ncols,figdim*nrows))
			plt.subplots_adjust(wspace=0.5)
			plt.subplots_adjust(hspace=0.5)
			# map the simulations onto axes here
			axes_map = {}
			for cc,c in enumerate(conditions[sg]):
				for ss,s in enumerate(merge_rule[c]):
					axes_map[s] = axes[cc][ss]
			panels = dict(axes=axes,fig=fig,axes_map=axes_map)

		# panel plot with simulations (rows) and decomposition style (columns)
		elif style=='rows_simulation_cols_decomp':
			figdim = 3
			sns = work.metadata.collections[sg]
			nrows,ncols = len(sns),len(decomps_order)
			fig,axes = plt.subplots(nrows=nrows,ncols=ncols,figsize=(figdim*ncols,figdim*nrows))
			plt.subplots_adjust(wspace=0.6)
			plt.subplots_adjust(hspace=0.6)
			# map the simulations onto axes here
			axes_map = {}
			for snum,sn in enumerate(sns):
				for dnum,d in enumerate(decomps_order):
					axes_map[(sn,'decomp',d)] = axes[snum][dnum]
			panels = dict(axes=axes,fig=fig,axes_map=axes_map)

		# panel plot with protein residues on the rows and simulations on the column
		elif style=='rows_residue_cols_simulation':
			figdim = 3
			sns = work.metadata.collections[sg]
			nrows,ncols = len(residues)+1,len(sns)
			fig,axes = plt.subplots(nrows=nrows,ncols=ncols,figsize=(figdim*ncols,figdim*nrows))
			plt.subplots_adjust(wspace=0.6)
			plt.subplots_adjust(hspace=0.6)
			# map the simulations onto axes here
			axes_map = {}
			for snum,sn in enumerate(sns):
				# top row is the "all" column
				# note that mixing of the different data types happens here instead of in the plot
				#   function. the heterogeneous nature of the data justifies our use of the 
				#   "get the data here" section of the plot_histograms function below
				axes_map[(sn,'decomp',decomp_style)] = axes[0][snum]
				for rnum,(i,j) in enumerate(residues):
					axes_map[(sn,'decomp',decomp_style,'res',i,j)] = axes[rnum+1][snum]
			panels = dict(axes=axes,fig=fig,axes_map=axes_map)

		else: raise Exception('invalid panel style: %s'%style)
		return panels

	def plot_histograms(panels,figname,include_zero=True,uniform_max=False):
		fig = panels['fig']
		max_count = 0
		del_axes = []
		do_decomp_label = False
		# plot histograms for each axis
		for sn_base,ax in panels['axes_map'].items():

			# get the data here
			#! the following is a prototype for a score method multiplexer
			global post,post_protein
			# identify the correct set of results
			if isinstance(sn_base,tuple) and len(sn_base)==3 and sn_base[1]=='decomp':
				sn,decomp_style = sn_base[0],sn_base[2]
				do_decomp_label = True
				post_this = post[sn]
			elif isinstance(sn_base,tuple) and len(sn_base)==6:
				if sn_base[1]=='decomp' and sn_base[3]=='res':
					res_id = (sn_base[4],sn_base[5])
					sn,decomp_style = sn_base[0],sn_base[2]
					post_this = post_protein[sn][res_id]
				else: raise Exception('invalid')
			elif isinstance(sn_base,tuple):
				raise Exception('cannot interpret plot target: %s'%str(sn_base))
			# by default we include hydrogen bonds and salt bridges
			else: 
				sn = sn_base
				decomp_style = 'both'
				post_this = post[sn_original]
	
			# retrieve the data now that we selected the post and decomp
			try: hists_by_lipid,bins,counts_by_lipid,lipid_kinds = [
				post_this['histograms'][decomp_style][i]
				for i in 'hists_by_lipid,bins,counts_by_lipid,lipid_kinds'.split(',')]
			#! see "dangerous" above
			except:
				del_axes.append(sn_base)
				continue

			# begin plotting here			
			ax.set_title(work.metadata.meta[sn]['label_compact'])
			ind0 = 0 if include_zero else 1
			bottom = np.zeros(len(bins)-1)
			for lipid_resname in lipid_stack_order:
				if lipid_resname not in lipid_kinds: continue
				heights = hists_by_lipid[list(lipid_kinds).index(lipid_resname)]
				# note that the sum of the heights should be 1.0 after we divide by frames
				# note that order matters: only histogram after you sum particular bond types
				ax.bar(bins[ind0:-1],heights[ind0:],bottom=bottom[ind0:],
					color=lipid_colors[lipid_resname],width=1.0)
				bottom += heights
			ax.set_xticks(range(ind0,bins.max()))
			ax.set_xlim((ind0-0.5,bins.max()-0.5))
			ax.set_xlabel('bound lipids')
			ax.set_ylabel('observations')
			if do_decomp_label:
				ax_r = ax.twinx()
				ax_r.set_yticks([])
				ax_r.yaxis.set_label_position('right')
				ax_r.set_ylabel(decomps_titles[decomp_style],rotation=-90,labelpad=20.)
			ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)
			ax.tick_params(axis='x',which='both',top=False,bottom=False,labelbottom=True)
			if uniform_max:
				max_count = max([max(bins[:-1]),max_count])
		if uniform_max:
			for sn,ax in panels['axes_map'].items():
				ax.set_xlim((0,max_count))
		for sn in del_axes:
			fig.delaxes(panels['axes_map'][sn])
		picturesave(figname,work.plotdir,backup=False,version=True,dpi=300,form='png',meta={})

	if 'post' not in globals():
		post = {}
		sns_all = work.sns()
		for snum,sn in enumerate(sns_all):
			status('computing histograms',i=snum,looplen=len(sns_all))
			post[sn] = compute_histograms(sn)

	# select a supergroup
	sg = 'mdia2'

	# select plots to render
	do_decomposition_plot = 0
	do_default_histograms = 0

	# PLOT SERIES
	if do_default_histograms:
		panels = prepare_panels(style='rows_condition_cols_replicate',
			supergroup='mdia2')
		plot_histograms(panels=panels,figname='fig.histograms.%s'%sg)
	if do_decomposition_plot:
		panels = prepare_panels(style='rows_simulation_cols_decomp',
			supergroup='mdia2')
		plot_histograms(panels=panels,figname='fig.histograms.mdia2.decompose')

	###!!! DEV: protein plot

	# get the residues for this supergroup
	data.set('contacts',select={'cutoff':cutoff,
		'target':{'predefined':'lipid heavy'},
		'subject':{'predefined':'protein heavy'}})
	if 'post_protein' not in globals():
		post_protein = {}
		for sn in work.metadata.collections[sg]:
			span = prepare_span(sn)
			residues = span_to_protein_residues(span)
			post_protein[sn] = {}
			for res in residues:
				status('computing histograms for %s: %s'%(sn,str(tuple(res))))
				result = compute_histograms(sn,method='filter_protein',protein_residue_target=res)
				post_protein[sn][tuple(res)] = result
	panels = prepare_panels(style='rows_residue_cols_simulation',decomp_style='both',
		residues=residues,supergroup='mdia2')
	plot_histograms(panels=panels,include_zero=False,uniform_max=True,
		figname='fig.histograms.mdia2.residues.hydrogen_bonds_salt_bridges')
