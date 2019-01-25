#!/usr/bin/env python

import os
import MDAnalysis
import numpy as np
from omni import WorkSpace
from omni import HookHandler
from omni import basic_compute_loop
from omni import picturesave
from omni.plotter.utils import zoom_figure
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import itertools
import scipy
import scipy.optimize

### FUNCTIONS

def trim_graph_fat(points,edges):
	"""
	Given a list of points and edges that connect them we reindex both to ignore
	the points which do not participate in an edge.
	"""
	# we get unique points indices for the left and right side
	originals_u_l = np.unique(edges[:,0])
	originals_u_r = np.unique(edges[:,1])
	# reduced indices
	reds = np.concatenate([originals_u_l,originals_u_r])
	edges_subsel = [np.in1d(edges[:,0],originals_u_l),np.in1d(edges[:,1],originals_u_r)]
	# confirm that any unique lipid on the left is also bound to a cation on the right
	if not np.all(np.any(edges_subsel,axis=0)==np.all(edges_subsel,axis=0)): raise Exception('error')
	# get the subset of edges using the old indexing
	edges_subsel_original_inds = edges[np.where(np.all(edges_subsel,axis=0))]
	# the remapper will allow us to put the old edges indices into a mapping from old to new
	remapper = np.ones(len(points)).astype(int)*-1
	# the items in the old list which are not in the new list are -1 while the others are indexed by 
	#   the reduced items stored in the reds variable
	remapper[reds] = np.arange(len(reds)).astype(int)
	edges_reduced = remapper[edges]
	points_reduced = points[reds]
	return points_reduced,edges_reduced,reds

def edges_to_clusters(n_pts,edges):
	"""
	Compute clusters from a set of edges.
	"""
	#! clusters = np.ones(len(inds)).astype(int)*-1
	clusters = np.ones(n_pts).astype(int)*-1
	# loop over edges 
	n_clusters = 0
	for i,j in edges:
		if clusters[i]==-1 and clusters[j]==-1:
			clusters[i] = clusters[j] = n_clusters
			n_clusters += 1
		elif clusters[i]==-1 and clusters[j]!=-1:
			clusters[i] = clusters[j]
		elif clusters[i]!=-1 and clusters[j]==-1:
			clusters[j] = clusters[i]
		# merge clusters
		else:
			existing = [clusters[i],clusters[j]]
			for e in existing: 
				clusters[np.where(clusters==e)] = n_clusters
			n_clusters += 1
	# relabel all clusters with zero indicating singletons
	cluster_inds = dict(zip(range(len(clusters)),np.unique(clusters)))
	clusters_original = np.array(clusters)
	for k,v in cluster_inds.items(): 
		clusters[np.where(clusters_original==v)] = k
	return clusters

def unbounded_center(pts,vec):
	"""Get the centroid (not mass weighted) for points on an unbounded (periodic) plane."""
	# via https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
	com = []
	for d in range(2):
		remap = pts[:,d]/vec[d]*2*np.pi
		# contra the wiki use no negative
		com_this = np.arctan2(np.mean(np.sin(remap)),np.mean(np.cos(remap)))/2/np.pi*vec[d]
		# however we have to rewrap the points to get the absolute com 
		#   because the result above can be negative
		com.append(com_this)
	com = np.array(com)
	com += (com<0)*vec-(com>vec)*vec
	return com

def compute_cluster_network(fr,mn):
	"""
	Extract a lipid-cation-lipid bond network and reduce it to a network of lipid-lipid bonds mediated
	by the cations.
	"""
	global dat,dat_abs,post,data,sn,vecs
	global resids_all,imono,resids_to_mesh_ind
	global top_inds,ind_mono_to_mesh,ind_mesh_to_mono
	# get the targets and subjects in their native index
	left,right = np.array([post[sn]['incoming'][fr][k] for k in ['targets','subjects']])
	# get residues for the target lipids
	left_resids = post[sn]['target'].resids[left]
	# filter by lipids in the top monolayer
	left_resids_inds_top = np.where(np.in1d(left_resids,resids_all[imono==mn]))[0]
	left,right = left[left_resids_inds_top],right[left_resids_inds_top]
	# convert target selection index to resid to index in the lipid mesh monolayer
	left_r = resids_to_mesh_ind[post[sn]['target'].resids[left]]
	if -1 in left_r: raise Exception('indexing error')
	# retain the subject index from the original selection because these are cations
	right_ind = right
	# bonds formulated as monolayer mesh index with subject index
	bonds_this = np.transpose([left_r,right_ind])
	# unique bonds
	bonds_u = np.unique(bonds_this,axis=0)

	# collect distinct indices that could be observed on the left and right assuming disjoint points in each
	left_inds = ind_mono_to_mesh[mn]
	n_cations = len(post[sn]['subject'])
	right_inds = np.arange(n_cations)
	points_inds = np.arange(len(left_inds)+len(right_inds)).astype(int)
	lr_offset = len(left_inds)
	"""
	create a mapping from left/right indices to points_inds indices
	note that we assume that the points on the left (column 0) are lipids and those on the right (column 1)
	are cations. we stack the list of points together so that later we can operate on the clusters of
	nodes connected by edges as a generic structure when fixing jumps over periodic boundaries. the
	generic solution to these problems requires a distinction between lipids and cations (left/right) but
	we also wish to operate on the clusters with generic code
	"""
	remapper = {}
	remapper['left'] = dict(zip(left_inds,np.arange(len(left_inds)).astype(int)))
	remapper['right'] = dict(zip(right_inds,(np.arange(len(right_inds))+lr_offset).astype(int)))

	# reindex bonds for the main set of points. at this stage the points are indexed in a single list
	#   with the lipids (distinct "left" points) first followed by cations. later we concatenate them
	#   into the pts_order structure and operate on that
	bonds_reform = np.array([[remapper[j][i] for i in bonds_u.T[jj]] 
		for jj,j in enumerate(['left','right'])]).T

	# arguments to the cluster finder
	#! inds = points_inds
	edges = bonds_reform
	# cluster algorithm requires edges and a set of base points indexed from zero
	if np.any(points_inds!=np.arange(len(points_inds)).astype(int)):
		raise Exception('incoming points indices must be a proper range')
	clusters = edges_to_clusters(n_pts=len(points_inds),edges=bonds_reform)
	if False:
		# cluster algorithm requires edges and a set of base points indexed from zero
		if np.any(points_inds!=np.arange(len(points_inds)).astype(int)):
			raise Exception('incoming points indices must be a proper range')
		clusters = np.ones(len(inds)).astype(int)*-1
		# loop over edges 
		n_clusters = 0
		for i,j in edges:
			if clusters[i]==-1 and clusters[j]==-1:
				clusters[i] = clusters[j] = n_clusters
				n_clusters += 1
			elif clusters[i]==-1 and clusters[j]!=-1:
				clusters[i] = clusters[j]
			elif clusters[i]!=-1 and clusters[j]==-1:
				clusters[j] = clusters[i]
			# merge clusters
			else:
				existing = [clusters[i],clusters[j]]
				for e in existing: 
					clusters[np.where(clusters==e)] = n_clusters
				n_clusters += 1
		# relabel all clusters with zero indicating singletons
		cluster_inds = dict(zip(range(len(clusters)),np.unique(clusters)))
		clusters_original = np.array(clusters)
		for k,v in cluster_inds.items(): 
			clusters[np.where(clusters_original==v)] = k

	# once we have  the clusters we need to map back to the points themselves
	pts_order = np.concatenate([
		#! limit the points when selecting the dat
		dat['%d.%d.points'%(mn,fr)][:len(left_inds)],
		post[sn]['coords_subject'][fr],])

	# move points across periodic boundaries to minimize their second moment
	vec = vecs[fr][:2]


	#! deprecated in favor of a much more elegant method
	if False:

		clusters_u = np.unique(clusters)
		for cc,cnum in enumerate(clusters_u):
			# ignore the zero group which are singletons
			if cnum==0: continue
			# get points in one cluster
			this_cluster = np.where(clusters==cnum)
			pts_this = pts_order[this_cluster][:,:2]
			arg0 = np.zeros(len(pts_this)*2)
			def func(arg):
				x = pts_this + arg.astype(int).reshape((2,-1)).T*vec
				return np.sum((x-x.mean(axis=0))**2)
			ans = scipy.optimize.minimize(fun=func,x0=(arg0,),method='Powell')
			shifts = ans.x.astype(int).reshape((2,-1))
			if (shifts!=0).any():
				#! no change if you do this mistake: pts_order[this_cluster][:,:2] += shifts.T*vec
				pts_order[this_cluster,:2] += shifts.T*vec

	#! new method here
	clusters_u = np.unique(clusters)
	for cc,cnum in enumerate(clusters_u):
		# ignore the zero group which are singletons
		if cnum==0: continue
		# get points in one cluster
		this_cluster = np.where(clusters==cnum)
		pts_this = pts_order[this_cluster][:,:2]
		n_opts = len(pts_this)*2
		cog = unbounded_center(pts_this,vec)
		# minimize distances to the cog for all points
		if np.any(np.abs((pts_this-cog))>vec/2):
			pts_this += ((pts_this-cog)>vec/2)*vec*-1
			pts_this += ((pts_this-cog)<(-vec/2))*vec
			# save the result
			pts_order[this_cluster,:2] = pts_this

	# collapse borders by finding all cations bound to at least two lipids to make a link
	pts_red = pts_order[:len(left_inds)]
	# we reduce i.e. eliminate the points on the right (the cations)
	right_u = np.unique(edges[:,1])
	reduced_bonds = [j for j in [edges[np.where(edges[:,1]==i)[0],0] for i in right_u] if len(j)>1]
	edges_reduced_basic = np.concatenate([
		list(itertools.combinations(v,2))
		for v in reduced_bonds])
	clusters_lipids = edges_to_clusters(n_pts=len(pts_red),edges=edges_reduced_basic)

	# package
	points_reduced,edges_reduced,_ = trim_graph_fat(pts_order,edges)
	points_reduced_lipids,edges_reduced_lipids,reds = trim_graph_fat(pts_red,edges_reduced_basic)
	n_pivot = len(left_inds)
	n_total = len(points_inds)
	outgoing = dict(n_total=n_total,n_pivot=n_pivot,lipids_reindex=reds,
		clusters=clusters,pts=points_reduced,edges=edges_reduced,
		clusters_lipids=clusters_lipids,
		points_reduced_lipids=points_reduced_lipids,edges_reduced_lipids=edges_reduced_lipids)
	return outgoing

def clusters_to_radius_gyration(cluster_object):
	"""
	Compute radius of gyration for clusters in a single cluster object.
	"""
	pts = cluster_object['points_reduced_lipids']
	edges = cluster_object['edges_reduced_lipids']
	clusters = cluster_object['clusters_lipids']
	reindex = cluster_object['lipids_reindex']
	clusters_this = clusters[reindex]
	clusters_inds = np.unique(clusters_this)
	pts_grouped = [pts[np.where(clusters_this==i)][:,:2] for i in clusters_inds]
	rgs = [
		np.sqrt(np.sum(np.linalg.norm(this-this.mean(axis=0),axis=1)**2))
		for this in pts_grouped]
	return np.array(rgs)

def ptdins_bond_graphs_calc(**kwargs):

	"""
	LIPID ABSTRACTOR
	Reduce a bilayer simulation to a set of points.
	"""

	#@loader
	#def load():
	#	work = WorkSpace(analysis=True)
	#	data = work.plotload(plotname='ptdins_percolate')
	#	sns = work.sns()[::-1]

	#! hacking in quickly from a plot script
	global post,sn,vecs
	global dat,dat_abs,post,data,sn
	global resids_all,imono,resids_to_mesh_ind
	global top_inds,ind_mono_to_mesh,ind_mesh_to_mono
	post = {}
	work = kwargs['workspace']
	sn = kwargs['sn']
	kwargs_incoming = kwargs

	#! running this block in the loop causes failures of some kind
	# custom parameters. fetch_coordinates clones automacs
	#! need a way of sync'ing automacs if there are changes?
	from calcs.codes.fetch_coordinates import get_lipid_resnames
	lipid_resnames = get_lipid_resnames(work)

	distance_cutoff = work.metadata.variables['hydration_cutoffs'][
		work.metadata.meta[sn]['cation']]
	#! suffixes need sorted so we hardcode here
	suffixes = ['gro','xtc']
	if False:
		structure,trajectory = [os.path.join(work.postdir,
			'%s.%s'%(data.extras[sn]['slice_path'],j)) for j in suffixes]
	structure,trajectory = [os.path.join(work.postdir,
		'%s'%(kwargs[j])) for j in ['structure','trajectory']]
	uni = MDAnalysis.Universe(structure,trajectory)

	from calcs.codes.fetch_coordinates import fetch_coordinates_subject_target
	fetch_coordinates_subject_target
	lipid_sel = ' or '.join('resname %s'%r for r in lipid_resnames)
	# requires variables/selectors/cations in the metadata
	cation_sel = work.metadata.variables['selectors']['cations']

	#! this section and the one in contacts.py needs to be more modular
	#! faking the kwargs
	kwargs = dict(subject=cation_sel,object=lipid_sel,
		structure=structure,trajectory=trajectory)
	#!!! copied from contacts.py but this is too clumsy
	# get coordinates from the hook
	default_fetch_coordinates = {
		'import_target':'codes/fetch_coordinates.py',
		'function':'fetch_coordinates_subject_target'}
	fetch_coordinates = kwargs.get(
		'fetch_coordinates',default_fetch_coordinates)
	# hooks get the same arguments as compute functions
	coordinates = HookHandler(meta=kwargs,**fetch_coordinates).solve
	nframes = coordinates['nframes']
	vecs = coordinates['vecs']
	subject = coordinates['uni_selections']['subject']
	target = coordinates['uni_selections']['target']
	times = coordinates['times']

	from calcs.contacts import contacts_framewise
	coords_target = coordinates['coordinates_by_selection']['target']
	coords_subject = coordinates['coordinates_by_selection']['subject']
	compute_function = contacts_framewise
	out_args = {'distance_cutoff':distance_cutoff}
	looper = [dict(fr=fr,target=coords_target[fr],subject=coords_subject[fr],
		vec=vecs[fr],**out_args) for fr in range(nframes)]
	# incoming has atomwise subjects/targets lists for those within cutoff
	incoming = basic_compute_loop(compute_function,looper)
	post[sn] = dict(
		incoming=incoming,coords_subject=coords_subject,
		subject=subject,target=target,times=times,
		vecs=vecs,nframes=nframes)

	#! previously this was the postprocess function

	mn_top = 0 #! hardcoded
	if False:
		data.set('lipid_abstractor')
		dat_abs = data.this[sn]
		data.set('lipid_mesh')
		dat = data.this[sn]
	dat = kwargs_incoming['upstream']['lipid_mesh']
	dat_abs = kwargs_incoming['upstream']['lipid_abstractor']
	imono = dat['monolayer_indices']
	resnames = dat['resnames']
	resids_all = dat_abs['resids']
	# alternate formulation is faster, but this is really the numpy version of:
	#   resids_to_mesh_ind = dict([(i,ii) for ii,i in enumerate(resids_all)])
	resids_to_mesh_ind = (np.ones(max(1+np.unique(resids_all)))*-1).astype(int)
	for ii,i in enumerate(resids_all): resids_to_mesh_ind[i] = ii
	# indices for lipids in the top leaflet
	top_inds = np.where(imono==mn_top)[0]
	# mapping from the monolayer index (e.g. 0-299 for lipids in a leaflet of
	#   a bilayer with cholesterol where we ignore cholesterol) to the mesh 
	#   index (e.g. 0-599 when there might be 800 molecules total)
	ind_mono_to_mesh = [np.arange(imono.shape[0])[np.where(imono==mn)] 
		for mn in range(2)]
	# reverse mapping so we can go from mesh index (e.g. 0-599) to index
	#   in the list of monolayer points
	#! you could use a dict lookup e.g. dict(zip(*np.array([
	#!   ind_mono_to_mesh[mn],np.arange(ind_mono_to_mesh[mn].shape[0])])))
	# just pasting in arange into the blocks of 0 and 1 in imono
	ind_mesh_to_mono = np.array(imono)
	for mn in range(2):
		ind_mesh_to_mono[np.where(imono==mn)] = \
			np.arange(len(np.where(imono==mn)[0]))
	mn = mn_top
	looper = [dict(mn=mn_top,fr=fr) for fr in range(post[sn]['nframes'])]
	incoming = basic_compute_loop(compute_cluster_network,
		looper=looper,n_jobs=4,run_parallel=False)

	if False:

		# compute clusters
		for snum,sn in enumerate(sns):
			if 'cluster_results' not in post[sn]:
				print('status computing clusters for %s %d/%d'%(sn,snum+1,len(sns)))
				post[sn]['cluster_results'] = postprocess(sn)
			if 'cluster_radius_gyration' not in post[sn]:
				post[sn]['cluster_radius_gyration'] = basic_compute_loop(
					clusters_to_radius_gyration,looper=looper,run_parallel=False)
	post[sn]['cluster_results'] = incoming
	looper = [{'cluster_object':o} for o in post[sn]['cluster_results']]
	post[sn]['cluster_radius_gyration'] = basic_compute_loop(
		clusters_to_radius_gyration,looper=looper,run_parallel=False)

	# pack
	attrs,result = {},{}
	#! save from abstractor
	attrs['selector'] = kwargs_incoming['upstream']['lipid_abstractor']['selector']
	attrs['nojumps'] = kwargs_incoming['upstream']['lipid_abstractor']['nojumps']
	#! package in the dumb way
	cluster_results_keys = post[sn]['cluster_results'][0].keys()
	for fr in range(nframes):
		for key in cluster_results_keys:
			result['%d__%s'%(fr,key)] = np.array(post[sn]['cluster_results'][fr][key])
	for fr in range(nframes):
		result['%d__rgs'%(fr)] = np.array(post[sn]['cluster_radius_gyration'][fr])
	return result,attrs	
