#!/usr/bin/env python

import re
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
from ortho import status
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
	global dat,dat_abs,post,data
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
	pts_order_bak = np.array(pts_order)
	# move points across periodic boundaries to minimize their second moment
	vec = vecs[fr][:2]
	clusters_u = np.unique(clusters)
	cogs = [] #! debugging
	for cc,cnum in enumerate(clusters_u):
		# ignore the zero group which are singletons
		if cnum==0: continue
		# get points in one cluster
		this_cluster = np.where(clusters==cnum)
		pts_this = pts_order[this_cluster][:,:2]
		n_opts = len(pts_this)*2

		#! unused block follows
		if False:
			# edges in absolute indices
			edges_this = edges[np.where(np.all([np.in1d(edges[:,i],this_cluster[0]) for i in range(2)],axis=0))] 
			# convert edges from absolute to relative indices for pts_this
			ind_abs_to_rel = dict([(j,jj) for jj,j in enumerate(this_cluster[0])])
			# edges in relative indices
			edges_this_rel = np.array([[ind_abs_to_rel[i] for i in j] for j in edges_this])		
		#! moved to a function
		if False:
			### demo a method to get the real centroid
			# via https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
			com = []
			for d in range(2):
				remap = pts_this[:,d]/vec[d]*2*np.pi
				# contra the wiki use no negative
				com_this = np.arctan2(np.mean(np.sin(remap)),np.mean(np.cos(remap)))/2/np.pi*vec[d]
				# however we have to rewrap the points to get the absolute com 
				#   because the result above can be negative
				com.append(com_this)
			com = np.array(com)
			com += (com<0)*vec-(com>vec)*vec
		
		cog = unbounded_center(pts_this,vec)
		cogs.append(cog) #! debugging only
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
	outgoing = dict(n_total=n_total,n_pivot=n_pivot,lipids_reindex=reds,pts_order=pts_order,
		clusters=clusters,pts=points_reduced,edges=edges_reduced,
		clusters_lipids=clusters_lipids,
		points_reduced_lipids=points_reduced_lipids,edges_reduced_lipids=edges_reduced_lipids,
		cogs=cogs,
		) #! temporary
	#if np.any(pts_order!=pts_order_bak):
	#	import ipdb;ipdb.set_trace()
	return outgoing

def postprocess(sn):
	"""
	Compute the clusters for a simulation.
	"""
	global dat,dat_abs,post,data
	global resids_all,imono,resids_to_mesh_ind
	global top_inds,ind_mono_to_mesh,ind_mesh_to_mono
	mn_top = 0 #! hardcoded
	data.set('lipid_abstractor')
	dat_abs = data.this[sn]
	data.set('lipid_mesh')
	dat = data.this[sn]
	imono = dat['monolayer_indices']
	resnames = dat['resnames'].astype(str)
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
	looper = [dict(mn=mn_top,fr=fr) for fr in range(post[sn]['nframes'])[slice(*frameslice)]]
	incoming = basic_compute_loop(compute_cluster_network,
		looper=looper,n_jobs=4,run_parallel=False)
	return incoming

@loader
def load():
	work = WorkSpace(analysis=True)
	data = work.plotload(plotname='ptdins_percolate')
	sns = work.sns()[::-1]

	#! running this block in the loop causes failures of some kind
	# custom parameters. fetch_coordinates clones automacs
	#! need a way of sync'ing automacs if there are changes?
	from calcs.codes.fetch_coordinates import get_lipid_resnames
	lipid_resnames = get_lipid_resnames(work)

	post = {}
	for sn in sns:

		distance_cutoff = work.metadata.variables['hydration_cutoffs'][
			work.metadata.meta[sn]['cation']]

		#! suffixes need sorted so we hardcode here
		suffixes = ['gro','xtc']
		structure,trajectory = [os.path.join(work.postdir,
			'%s.%s'%(data.extras[sn]['slice_path'],j)) for j in suffixes]
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
		incoming = basic_compute_loop(compute_function,looper,run_parallel=False)
		post[sn] = dict(
			incoming=incoming,coords_subject=coords_subject,
			subject=subject,target=target,times=times,
			vecs=vecs,nframes=nframes)

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

def weird_count(sn,fr):
	outgoing = post[sn]['cluster_results'][fr]
	pts,edges = outgoing['points_reduced_lipids'],outgoing['edges_reduced_lipids']
	lipids_these = np.where(post[sn]['cluster_results'][fr]['clusters_lipids']>0)[0]
	reindex = np.zeros(len(pts)).astype(int)
	reindex[np.unique(outgoing['edges_reduced_lipids'])] = np.arange(len(lipids_these))
	# corresponding cluster index for this frame is
	cinds = post[sn]['cluster_results'][fr]['clusters_lipids'][lipids_these]
	bonds = resnames[lipids_these][reindex[edges]]
	return dict(bonds=bonds,cinds=cinds)

def load():
	work = WorkSpace(analysis=True)
	data = work.plotload(plotname='ptdins_percolate')
	sns = work.sns()[::-1] # no need to go backwards
	#! running this block in the loop causes failures of some kind
	# custom parameters. fetch_coordinates clones automacs
	#! need a way of sync'ing automacs if there are changes?
	from calcs.codes.fetch_coordinates import get_lipid_resnames
	lipid_resnames = get_lipid_resnames(work)

if __name__=='__main__': 

	if 'post_drifts' not in globals():

		global frameslice
		frameslice = (770,790)
		if 'data_up' not in globals():
			data_up = work.plotload(plotname='ptdins_bond_graphs_calc')
			post_rgs,nframes_by_sn = {},{}
			for snum,sn in enumerate(sns):
				data.set('lipid_abstractor')
				dat_abs = data.this[sn]
				data_up.set('ptdins_bond_graphs_calc')
				dat = data_up.this[sn]
				nframes = dat_abs['nframes']
				nframes_by_sn[sn] = nframes
				rgs = np.concatenate([dat['%d__rgs'%fr] for fr in range(nframes)])
				post_rgs[sn] = rgs

		mn_top = 0
		cluster_contents ={}
		for sn in sns:
			imono = data.this[sn]['monolayer_indices']
			data.set('lipid_abstractor')
			resnames = data.this[sn]['resnames'].astype(str)[imono==mn_top]
			resnames_int = np.zeros(len(resnames)).astype(int)
			resname_to_int = dict([(j,i) for i,j in enumerate(np.unique(resnames))])
			for resname in resname_to_int:
				resnames_int[np.where(resnames==resname)] = resname_to_int[resname]
			rgs = post_rgs[sn]
			counts_resnames = np.zeros((len(rgs),len(resname_to_int)))
			n_resnames = np.unique(resnames)
			counter = 0
			nframes = nframes_by_sn[sn]
			for fr in range(nframes):
				clusters_this = data_up.this[sn]['%d__clusters_lipids'%fr]
				resnames_this = [resnames_int[np.where(clusters_this==i)[0]] 
					for i in np.unique(clusters_this)[1:]]
				for r in resnames_this:
					this,counts = np.unique(r,return_counts=True)
					counts_resnames[counter][this] = counts
					counter += 1
			cluster_sizes = counts_resnames.sum(axis=1)
			#! i=2;rgs[np.all((rgs<bins[i],rgs>=bins[i-1]),axis=0)]
			cluster_contents[sn] = dict(cluster_sizes=cluster_sizes,counts_resnames=counts_resnames)

		redo = 1
		# compute clusters
		for snum,sn in enumerate(sns):
			if 'cluster_results' not in post[sn] or redo:
				print('status computing clusters for %s %d/%d'%(sn,snum+1,len(sns)))
				post[sn]['cluster_results'] = postprocess(sn)
			if 'cluster_radius_gyration' not in post[sn] or redo:
				looper = [{'cluster_object':o} for o in post[sn]['cluster_results']]
				post[sn]['cluster_radius_gyration'] = basic_compute_loop(
					clusters_to_radius_gyration,looper=looper,run_parallel=False)
			if 'cluster_bond_propensities' not in post[sn] or redo:
				looper = [dict(fr=fr,sn=sn) for fr in range(len(post[sn]['cluster_results']))]
				post[sn]['cluster_bond_propensities'] = basic_compute_loop(
					weird_count,looper=looper,run_parallel=False)

		cluster_ind_to_fr = {}
		for sn in sns:
			rgs = post_rgs[sn]
			nframes = nframes_by_sn[sn]
			n_clusters_per_frame = [np.unique(data_up.this[sn]['%d__clusters_lipids'%fr]).shape[0]-1 
				for fr in range(nframes)]
			max_rg_ind = rgs.argmax()
			cluster_ind_to_fr[sn] = np.concatenate([np.ones(i)*ii 
				for ii,i in enumerate(n_clusters_per_frame)]).astype(int)
			# 775. adjust frameslice above to render it with the code below
			fr = cluster_ind_to_fr[sn][max_rg_ind]

		sns = ['membrane-v531','membrane-v532']
		colors = {'DOPE':'gray','DOPS':'blue','PI2P':'red'}
		#! broken! so we have to do this
		fig = plt.figure(figsize=(8,8))
		gs = gridspec.GridSpec(1,2,width_ratios=[5,1])
		axes = [plt.subplot(g) for g in gs]
		for snum,sn in enumerate(sns):
			data.set('lipid_abstractor')
			dat = data.this[sn]
			nframes = dat['nframes']
			ax = axes[snum]
			rgs = post_rgs[sn]
			lev = -1 # rounding level
			bins = np.arange(0,np.ceil(rgs.max()/10**lev+1)*10**lev,10**lev)
			counts,_ = np.histogram(rgs,bins=bins)
			mids = (bins[1:]+bins[:-1])/2.
			#! ax.plot(mids,counts)
			#! ax.bar(mids,counts)
			# plot cluster lipid counts instead
			cluster_sizes = cluster_contents[sn]['cluster_sizes']
			counts_resnames = cluster_contents[sn]['counts_resnames']
			bins = np.arange(0,cluster_sizes.max()+2)-0.5
			counts,_ = np.histogram(cluster_sizes,bins=bins)
			mids = (bins[1:]+bins[:-1])/2.
			#! you can plot the total count here
			#! ax.bar(mids,counts)
			# breaking down each bar by the lipid type
			# prepare a stacked bar plot
			counts_mean = np.array([counts_resnames[np.where(cluster_sizes==i)] for i in mids])
			#! need resname_to_int
			counts_mean = np.array([j.mean(axis=0) if len(j)>0 else np.zeros(len(resname_to_int)) for j in counts_mean]).T
			#!!!
			weighted = (counts*counts_mean/mids)
			weighted[np.isnan(weighted)] = 0
			#! no! not a cumulative sum! weighted = np.cumsum(weighted,axis=0)/nframes
			weighted = weighted/nframes
			base = np.zeros((len(weighted.T)))
			for rnum,row in enumerate(weighted):
				ax.bar(mids,row,bottom=base,width=1.0,
					color=colors[dict([(j,i) for i,j in resname_to_int.items()])[rnum]])
				base += row
			ax.set_title(work.metadata.meta[sn]['ion_label'])
			#! ax.set_ylim((0,1.05*weighted.max()))
			ax.set_xlabel('lipids in a cluster')
			ax.set_ylabel(r'$\langle N_{clusters} \rangle$')
			ax.set_xlim((ax.set_xlim(2-0.5,mids[2:].ptp()+2.0+0.5)))
			#! checking that the number of clusters is right
			#! weighted.sum(axis=0).sum() # this is the sum of the histogram and should be the average clusters per frame, 45.8 for v532
		axes[1].set_ylim(axes[0].get_ylim())
		picturesave('TMP7',work.plotdir,backup=False,version=True,meta={},extras=[],dpi=300,form='png')

		colors = {'DOPE':'gray','DOPS':'blue','PI2P':'red'}
		fr = 775
		fr = fr-frameslice[0]
		# figure
		ax_pad = 0.2
		zoom_fac = 1.
		figscale = 10.0
		ncols = 2
		figprop = 0.5*ncols
		zorder = {'pbc':0,'lines':1,'lipids':2,'cations':3}
		color_lipids = {'DOPE':'b','DOPS':'m','PI2P':'r','DOPC':'g','POPC':'g'}
		fig = plt.figure(figsize=(figscale*figprop,figscale))
		gs = gridspec.GridSpec(2,2)
		axes = [plt.subplot(g) for g in gs]
		for snum,sn in enumerate(sns):
			outgoing = post[sn]['cluster_results'][fr]
			title = work.metadata.meta[sn]['ion_label']
			ax = axes[0+ncols*snum]
			ax.set_title('bound cations, '+title)
			pts,edges = outgoing['pts'],outgoing['edges']
			# start color calculation
			cation_color = 'k'
			cmap_name = 'jet'
			clusters = outgoing['clusters']
			n_pivot = outgoing['n_pivot']
			inds = np.arange(outgoing['n_total']).astype(int)
			color_inds = np.unique(clusters)
			np.random.shuffle(color_inds)
			color_hash = dict([(i,mpl.cm.__dict__[cmap_name](float(ii)/len(color_inds)))
				for ii,i in enumerate(color_inds)])
			# the conditional makes the cations red because they are sequentially last
			color = [color_hash[clusters[i]] 
				if i<n_pivot else cation_color 
				for i in inds if clusters[i]!=0]
			# end color caclulation
			ax.scatter(*pts[:,:2].T,color=color,s=10)
			#! cogs=outgoing['cogs'];ax.scatter(*np.array(cogs).T,color='r',s=20)
			lines = np.array([pts[edges[:,0]][:,:2],
				pts[edges[:,1]][:,:2]]).transpose((1,0,2))
			lc = mpl.collections.LineCollection(lines,color='k',lw=0.5,zorder=zorder['lines'])
			ax.add_collection(lc)		
			
			#! repetitive with above minus color calculation
			ax = axes[1+ncols*snum]
			ax.set_title('bridged lipids, '+title)
			#! pts_red,edges_reduced = outgoing['p2'],outgoing['e2']
			pts,edges = outgoing['points_reduced_lipids'],outgoing['edges_reduced_lipids']
			#! pts = outgoing['pts_order']
			#! trying to get colors here
			lipids_these = np.where(post[sn]['cluster_results'][fr]['clusters_lipids']>0)[0]
			pts_these = pts[np.unique(outgoing['edges_reduced_lipids'])]
			#! ax.scatter(*pts[:,:2].T,color='k',s=10)
			imono = data.this[sn]['monolayer_indices']
			data.set('lipid_abstractor')
			resnames = data.this[sn]['resnames'].astype(str)[imono==mn_top]
			colors_this = [colors[resnames[i]] for i in lipids_these]
			ax.scatter(*pts_these[:,:2].T,color=colors_this,s=10)
			lines = np.array([pts[edges[:,0]][:,:2],
				pts[edges[:,1]][:,:2]]).transpose((1,0,2))
			lc = mpl.collections.LineCollection(lines,color='k',lw=0.5,zorder=zorder['lines'])
			ax.add_collection(lc)		
		zoom_figure(fig,zoom_fac)
		plt.subplots_adjust(wspace=0.35,hspace=0.35)
		picturesave('TMP6',work.plotdir,backup=False,version=True,meta={},extras=[],dpi=300,form='png')

	def compute_per_lipid(sn,fr):
		"""
		objective:
			for each lipid, identify its type
			lipid to cluster size
			lipid to degree
			lipid to degree with other types
		development history:
			wrote the following with outgoing and post and then moved it to data_up
			wrapped it in a function for fast computation
		"""
		#! sn = 'membrane-v532'
		#! originally testing this code on post: outgoing = post[sn]['cluster_results'][fr]
		# we have some kind of funky mapping where some of the points are not participating in edges so we reindex
		#! pts,edges = outgoing['points_reduced_lipids'],outgoing['edges_reduced_lipids']
		# switching to get the data from data_up instead of outgoing
		#! fr = 500 disabled fr,sn when we wrap it in a function
		pts = data_up.this[sn]['%d__points_reduced_lipids'%fr]
		edges = data_up.this[sn]['%d__edges_reduced_lipids'%fr]
		# e.g. 373 points and 388 edges. this must include ghosts?
		#! previously: clusters_lipids = post[sn]['cluster_results'][fr]['clusters_lipids']
		clusters_lipids = data_up.this[sn]['%d__clusters_lipids'%fr]
		lipids_these = np.where(clusters_lipids>0)[0]
		# e.g. 261 lipids in clusters
		# determine the size
		#! clusters_lipids = post[sn]['cluster_results'][fr]['clusters_lipids']
		#! old school: dict(zip(*np.unique(clusters_lipids,return_counts=True)))
		indices,sizes_per_index = np.unique(clusters_lipids,return_counts=True)
		# size mapper
		sizes = np.zeros(len(indices))
		sizes[indices] = sizes_per_index
		# e.g. sizes is 41 because there are 40 clusters
		# get the cluster size per lipid over 261 lipids in clusters
		sizes_by_lipid = sizes[clusters_lipids[clusters_lipids>0]]
		# resnames for the lipids in clusters
		resnames_by_lipid = resnames[lipids_these]
		# now we have the cluster sizes and residue types for each lipid participant
		# find the degree of each lipid
		# but first, get the bond types per edge because that might help
		reindex = np.zeros(len(pts)).astype(int)
		reindex[np.unique(edges)] = np.arange(len(lipids_these))
		# corresponding cluster index for this frame is
		#! previously cif = post[sn]['cluster_results'][fr]['clusters_lipids'][lipids_these]
		cif = clusters_lipids[lipids_these]
		# bonds corresponds to edges
		bonds = resnames[lipids_these][reindex[edges]]
		# back to the question of degree: edges is indexed on all lipids
		# need to get the edges list in the reduced space, of e.g. 261 lipids
		#! eg_lipid = lipids_these[-1]
		# degree for this lipid
		#! degree = edges[np.where(np.any(edges==eg_lipid,axis=1))]
		#!? is there redundancy in the edges list? it looks like this might be legitimate double-connections
		degrees = np.array([edges[np.where(np.any(edges==l,axis=1))].shape[0] for l in lipids_these])
		# now we are ready to run this calculation on the entire trajectory
		return dict(rbl=resnames_by_lipid,degrees=degrees,sbl=sizes_by_lipid,resids=lipids_these)

	if 'post_drifts' not in globals():
		post_drifts = {}
		for sn in sns:
			data.set('lipid_abstractor')
			dat = data.this[sn]
			nframes = dat['nframes']
			looper = [dict(fr=fr,sn=sn) for fr in range(nframes)]
			incoming = basic_compute_loop(compute_per_lipid,looper=looper,run_parallel=False)
			# formulate the results into trajectories
			resids_all = np.unique(np.concatenate([incoming[i]['resids'] for i in range(nframes)]))
			resids_all_reindex = (np.ones(np.max(resids_all)+1)*-1).astype(int)
			resids_all_reindex[resids_all] = np.arange(len(resids_all))
			resids_reindexer = []
			for fr in range(nframes):
				resids = incoming[fr]['resids']
				resids_reindex = (np.ones(np.max(resids_all)+1)*-1).astype(int)
				resids_reindex[resids_all_reindex[resids]] = np.arange(len(resids))
				resids_reindexer.append(resids_reindex)
			imono = data.this[sn]['monolayer_indices']
			data.set('lipid_abstractor')
			resnames = data.this[sn]['resnames'].astype(str)[imono==mn_top]
			resnames_all = dict([(i,ii) for ii,i in enumerate(np.unique(resnames))])
			sizes = np.zeros((nframes,len(resids_all))).astype(int)
			degrees = np.zeros((nframes,len(resids_all))).astype(int)
			resnames_summary = np.zeros((nframes,len(resids_all))).astype(int)
			for fr in range(nframes):
				resids = incoming[fr]['resids']
				resids_reindex = resids_reindexer[fr]
				sizes[fr][resids_all_reindex[resids]] = incoming[fr]['sbl']
				degrees[fr][resids_all_reindex[resids]] = incoming[fr]['degrees']
				#! this is redundant by frame I think
				resnames_summary[fr][resids_all_reindex[resids]] = np.array([resnames_all[i] for i in incoming[fr]['rbl']])
			post_drifts[sn] = {}
			post_drifts[sn]['sizes'] = sizes
			post_drifts[sn]['degrees'] = degrees
			#! this is probably redundant above
			post_drifts[sn]['resnames_summary'] = resnames_summary.mean(axis=0).astype(int)
			post_drifts[sn]['resnames_toc'] = resnames_all

	if 'once_through' not in globals():
		# this took less than 1 hour, finishing and converting the nugget for the compute_per_lipid function
		colors = {'DOPE':'gray','DOPS':'blue','PI2P':'red'}
		# figure
		ax_pad = 0.2
		zoom_fac = 1.
		figscale = 10.0
		ncols = 2
		figprop = 0.5*ncols
		zorder = {'pbc':0,'lines':1,'lipids':2,'cations':3}
		color_lipids = {'DOPE':'b','DOPS':'m','PI2P':'r','DOPC':'g','POPC':'g'}
		fig = plt.figure(figsize=(figscale*figprop,figscale))
		gs = gridspec.GridSpec(2,ncols)
		axes = [plt.subplot(g) for g in gs]
		for snum,sn in enumerate(sns):
			title = work.metadata.meta[sn]['ion_label']
			# get colors
			resnames_all = post_drifts[sn]['resnames_toc']
			resnames_all_r = dict([(j,i) for i,j in resnames_all.items()])
			colors_this = [colors[resnames_all_r[i]] for i in post_drifts[sn]['resnames_summary']]
			ax = axes[ncols*snum+0]
			ax.set_title('degree by lipid, %s'%title)
			for ii,i in enumerate(post_drifts[sn]['degrees'].T):
				ax.plot(i,color=colors_this[ii],alpha=0.5)
			ax = axes[ncols*snum+1]
			ax.set_title('cluster size by lipid, %s'%title)
			for ii,i in enumerate(post_drifts[sn]['sizes'].T):
				ax.plot(i,color=colors_this[ii],alpha=0.5)

		zoom_figure(fig,zoom_fac)
		plt.subplots_adjust(wspace=0.35,hspace=0.35)
		picturesave('TMP8',work.plotdir,backup=False,version=True,meta={},extras=[],dpi=300,form='png')

	if 'once_through' not in globals() or True:

		# iterate over the final plots
		sns = work.sns()
		# sns = ['membrane-v565'] #! hacking a bit
		for sn in sns:

			"""
			demo for cluster visualization method
				need a way to identify unique clusters
				problem is that the cluster index might change from one frame to the next
					because the indexing itself has no consistency. clusters are indexed algorithmically/stochastically
				a cluster is uniquely identified by its constituents
				the post_drifts creation routine above has uniquely identified lipids
					or rather each lipid has the cluster size
				but this information (from post_drifts) is not enough to identify the cluster
					because it only links to the cluster size
				first step is to get a single cluster and indices inside 
				note that the clsuters_lipids itemizes a lipid index for all lipids
					these indices could be different on the next frame
					so how to match them up?
					make a list of *participants* in each cluster for each frame
					then you can compare them to see if identical
					this first pass will be slow, I can tell
			"""

			# develop a unique listing of each item in a cluster
			# note that cluster_specs is indexed by absolute lipid residue index
			cluster_specs = []
			nframes = nframes_by_sn[sn]
			for fr in range(nframes):
				status('collecting clusters',i=fr+1,looplen=nframes)
				pts = data_up.this[sn]['%d__points_reduced_lipids'%fr]
				edges = data_up.this[sn]['%d__edges_reduced_lipids'%fr]
				clusters_lipids = data_up.this[sn]['%d__clusters_lipids'%fr]
				cluster_inds = np.arange(1,max(clusters_lipids))
				cl = [tuple(np.where(clusters_lipids==i)[0]) for i in cluster_inds]
				cluster_specs.append(cl)
			# map each cluster to a code
			codes = 0
			cluster_toc = {}
			for fr in range(nframes):
				for i in cluster_specs[fr]:
					if i not in cluster_toc:
						cluster_toc[i] = codes
						codes += 1
			# prepare colors
			sizes_toc = dict([(i,len(i)) for i in cluster_toc])
			cluster_order = dict([(i,ii) for ii,i in enumerate(cluster_toc.keys())])
			codes_by_size = []
			for fr in range(nframes):
				codes = [cluster_toc[i] for i in cluster_specs[fr]]
				sizes = [sizes_toc[i] for i in cluster_specs[fr]]
				codes_by_size.append(list(zip(sizes,codes)))
			# get the maximum size for all clusters with a buffer for length
			total_size_max = max([sum(list(zip(*i))[0])+len(i) for i in codes_by_size])
			# generate colors
			range_this = list(range(len(cluster_toc)))
			np.random.shuffle(range_this)
			colors_mapped = [mpl.cm.__dict__['jet'](float(i)/(len(cluster_toc)-1.)) for i in range_this]
			# get residue names
			imono = data.this[sn]['monolayer_indices']
			data.set('lipid_abstractor')
			resnames = data.this[sn]['resnames'].astype(str)[imono==mn_top]
			resnames_all = dict([(i,ii) for ii,i in enumerate(np.unique(resnames))])
			colors_these = [mpl.colors.to_rgba(i) for i in ['gray','blue','red']]

			# make a blank canvas
			canvas_unsort = np.zeros((nframes,total_size_max,4))
			canvas_sort = np.zeros((nframes,total_size_max,4))
			canvas_resname_color = np.zeros((nframes,total_size_max,4))
			"""
			first built the plot with a "natural" ordering that is somewhat stochastic but largely reflects the
			order in which clusters were identified, and hence is somewhat coherent
			then I changed the below inner for loop to work on size and this made things cohere more, since 
			the size carries some ability to uniquely identify the clusters
			the next modification is to change from colors to residue namein a second axis
			"""
			for fr in range(nframes):
				y = 0
				for size,index in codes_by_size[fr]:
					canvas_unsort[fr,y:y+size,:] = colors_mapped[index]
					y += size + 1
				y = 0
				#? stacking by a sorted size makes the plot look more uniform, but not ncessarily better?
				for size,index in sorted(codes_by_size[fr],key=lambda x:x[0]):
					canvas_sort[fr,y:y+size,:] = colors_mapped[index]
					y += size + 1
				y = 0
				# sorting by size and color by residue name
				for size,index in sorted(codes_by_size[fr],key=lambda x:x[0]):
					this_cluster = list(cluster_toc.keys())[index]
					this_colors = np.array([[colors_these[k] for k in j] 
						for j in [sorted([resnames_all[resnames[i]] 
							for i in this_cluster])]])
					canvas_resname_color[fr,y:y+size,:] = this_colors
					y += size + 1

			if 0:

				# figure
				ax_pad = 0.2
				zoom_fac = 1.
				figscale = 10.0
				ncols = 1
				nrows = 4
				figprop = 0.5*ncols
				fig = plt.figure(figsize=(figscale*figprop,figscale))
				gs = gridspec.GridSpec(nrows,ncols)
				axes = [plt.subplot(g) for g in gs]
				
				ax = axes[0]
				ax.imshow(canvas_unsort.transpose((1,0,2)),origin='lower',interpolation='nearest')

				ax = axes[1]
				ax.imshow(canvas_sort.transpose((1,0,2)),origin='lower',interpolation='nearest')

				ax = axes[2]
				ax.imshow(canvas_resname_color.transpose((1,0,2)),origin='lower',interpolation='nearest')

				ax = axes[3]

				# original funnel plot has big error bars that grow and obscure growth in the cluster size
				stats = np.zeros((nframes,2))
				for fr in range(nframes):
					item = codes_by_size[fr]
					raw = list(zip(*item))[0]
					stats[fr] = np.mean(raw),np.std(raw)
				ax.plot(stats[:,0],c='r',lw=3,zorder=2)
				if 1:
					ax.fill_between(np.arange(nframes),stats[:,0]+stats[:,1],stats[:,0]-stats[:,1],color='b',alpha=0.3,zorder=0,lw=0)

				max_counts = 0
				for fr in range(nframes):
					sizes,counts = np.unique(list(zip(*codes_by_size[fr]))[0],return_counts=True)
					max_counts = max(counts) if max(counts)>max_counts else max_counts
				for fr in range(nframes):
					sizes,counts = np.unique(list(zip(*codes_by_size[fr]))[0],return_counts=True)
					colors = np.array([mpl.colors.to_rgba('k') for s in sizes])
					colors[:,-1] = counts/max_counts
					ax.scatter([fr for i in sizes],sizes,s=1,c=colors,zorder=1)

				zoom_figure(fig,zoom_fac)
				plt.subplots_adjust(wspace=0.35,hspace=0.35)
				picturesave('TMP9',work.plotdir,backup=False,version=True,meta={},extras=[],dpi=300,form='png')


			max_sizes = 0
			max_counts = 0
			max_sc = 0
			for fr in range(nframes):
				sizes,counts = np.unique(list(zip(*codes_by_size[fr]))[0],return_counts=True)
				max_sizes = max(sizes) if max(sizes)>max_sizes else max_sizes
				max_counts = max(counts) if max(counts)>max_counts else max_counts
				combo = sizes*counts
				max_sc = max(combo) if max(combo)>max_sc else max_sc

			n_lipids = 3 #! hardcoded for now

			# stacked bar chart of composition of each sized cluster
			bars_comp = []
			for fr in [0,nframes-1][:]:
				sizes,counts = np.unique(list(zip(*codes_by_size[fr]))[0],return_counts=True)
				bars = []
				for size,count in zip(sizes,counts):
					props = np.zeros(n_lipids)
					participants = [i for i in cluster_specs[fr] if len(i)==size]
					if len(participants)==0: continue
					kinds,counts_residue = np.unique(
						# codes for each lipid participating in all clusters of this size
						[resnames_all[resnames[i]] for i in np.concatenate(participants)],
						return_counts=True)
					# save the cluster size, the residue codes, the counts for those codes, 
					#   and the total participants
					bars.append((size,kinds,counts_residue,size*count))
				bars_comp.append(bars)

			pixels_w = 5000.
			# build pixel maps of the lipids participating in each size
			# note that this was retired in favor of stacked horizontal bar graphs
			#   and then it was un-retired when I couldn't get the axes to look right
			if False:
				canvases = []
				for fr in [0,nframes-1][:]:
					sizes,counts = np.unique(list(zip(*codes_by_size[fr]))[0],return_counts=True)
					#! changing to nframes size canvas = np.zeros((max_sizes+1,max_sc,4))
					fac = pixels_w/max_sc
					canvas = np.zeros((max_sizes+1,int(pixels_w),4))
					for ii,i in enumerate(sizes):
						out = sizes[ii]*counts[ii]
						out = int(sizes[ii]*counts[ii]*fac)
						# plotting the equivalent of a stacked bar plot
						canvas[i][:out] = mpl.colors.to_rgba('k')
					canvases.append(canvas)
				canvases[0] = canvases[0][:,::-1,:]
			# unretired version uses bars_comp above
			canvases,counts_max = [],[]
			for side,bars in enumerate(bars_comp):
				# interject a way to figure out the 
				vals = np.zeros((n_lipids,max_sizes+1))
				for lnum in range(n_lipids):
					for size,kinds,counts,total in bars:
						if lnum in kinds:
							vals[lnum][size] = counts[list(kinds).index(lnum)]/sum(counts)*total
				vals = vals[::-1] #! reverse lipid order here and with colors_these
				vals_cumulative = np.cumsum(np.concatenate(([np.zeros(max_sizes+1)],vals)),axis=0)
				canvas = np.zeros((max_sizes+1,int(pixels_w),4))
				fac = pixels_w/vals_cumulative.max()
				# save this for the x-axis limits later
				counts_max.append(vals_cumulative.max())
				for lnum in range(n_lipids):
					for ii,(i,j) in enumerate(zip(
						(vals_cumulative[lnum]*fac).astype(int),(vals_cumulative[lnum+1]*fac).astype(int)
						)):
						canvas[ii,i:j] = colors_these[::-1][lnum]
				canvases.append(canvas)
			#! invert_xaxis instead: canvases[0] = canvases[0][:,::-1]
			
			# build a pixel map of the numbers of clusters of each size by time
			canvas_traj = np.zeros((max_sizes+1,nframes,4))
			canvas_traj_alt = np.zeros((max_sizes+1,nframes))
			for fr in range(nframes):
				sizes,counts = np.unique(list(zip(*codes_by_size[fr]))[0],return_counts=True)
				for s,c in zip(sizes,counts):
					canvas_traj[s,fr] = mpl.cm.jet(float(c)/max_counts)
					canvas_traj_alt[s,fr] = float(c)

			# calculate the total clustered lipids
			totals = np.zeros(nframes)
			totals_by_lipid = np.zeros((n_lipids,nframes,))
			for fr in range(nframes):
				sizes,counts = np.unique(list(zip(*codes_by_size[fr]))[0],return_counts=True)
				totals[fr] = sum(sizes*counts)
				resnames_inds,counts = np.unique(
					resnames_int[np.concatenate(cluster_specs[fr])],return_counts=True)
				totals_by_lipid[resnames_inds,fr] = counts

			if 0:

				layout = {
					'panels':{'start':2,'end':3,'trajectory':1,'counts':0},
					'ncols':2,'nrows':2}

				# figure
				ax_pad = 0.2
				zoom_fac = 1.
				figscale = 10.0
				ncols = layout['ncols']
				nrows = layout['nrows']
				figprop = 0.5*ncols
				fig = plt.figure(figsize=(figscale*figprop,figscale))
				gs = gridspec.GridSpec(nrows,ncols)
				axes = [plt.subplot(g) for g in gs]

				time_start,time_end = data.extras[sn]['start']/1000.,data.extras[sn]['end']/1000.
				max_time = (time_end-time_start) # ns

				# previous 1-plot option
				if False:
					catted = np.concatenate((
						canvases[0].transpose((1,0,2)),
						canvas_traj.transpose((1,0,2)),
						canvases[1].transpose((1,0,2))))
					ax = axes[0]
					ax.imshow(catted.transpose((1,0,2)),origin='lower',interpolation='nearest')
				aspect,aspect_flank = max_time/max_sizes,pixels_w/max_sizes

				#! retired in favor of stacked bar plots
				if 1:
					for i,j,when,flop in [
						(layout['panels']['start'],0,r'$t=%d\,ns$'%time_start,-1),
						(layout['panels']['end'],1,r'$t=%d\,ns$'%time_end,1)]:
						ax = axes[i]
						ax.set_title(when)
						ax.set_ylabel('lipids per cluster')
						ax.set_xlabel(r'lipids per cluster $\times$ clusters')
						aspect_flank = float(counts_max[j])/max_sizes
						ax.imshow(canvases[j],origin='lower',interpolation='nearest',
							aspect=aspect_flank,extent=[0,counts_max[j],0-0.5,max_sizes+1-0.5])
						ax.set_ylim((2-0.5,max_sizes+1-0.5))
						if flop==-1: ax.invert_xaxis()
				#! stacked bar plot method is difficult to align with imshow
				else:
					for axnum,which,flop,when in [(0,0,-1,r'$t=%d\,ns$'%time_start),(2,1,1,r'$t=%d\,ns$'%time_end)]:
						ax = axes[axnum]
						bars = bars_comp[which]
						vals = np.zeros((n_lipids,max_sizes+1))
						for lnum in range(n_lipids):
							for size,kinds,counts,total in bars:
								if lnum in kinds:
									vals[lnum][size] = counts[list(kinds).index(lnum)]/sum(counts)*total
						vals = vals[::-1] #! reverse lipid order here and with colors_these
						vals_cumulative = np.cumsum(np.concatenate(([np.zeros(max_sizes+1)],vals)),axis=0)
						for lnum in range(n_lipids)[::flop]:
							ax.barh(range(max_sizes+1),flop*(vals_cumulative[lnum+1]-vals_cumulative[lnum]),left=flop*vals_cumulative[lnum],
								color=colors_these[::-1][lnum],zorder=2,align='center',height=1.0)
						if False:
							for i in [1,2,3]:
								ax.plot(flop*np.arange((max_sizes+1/i))*i,np.arange((max_sizes+1/i)),c='k',zorder=1)
						#! from barh: ax.set_ylim((1.5,max_sizes-0.5))
						ax.set_title(when)
						ax.set_ylabel('lipids per cluster')
						ax.set_xlabel('lipids (total)')
						if flop==-1:
							ax.set_xticklabels([int(-1*i) for i in ax.get_xticks()])

				ax = axes[layout['panels']['trajectory']]
				ax.set_title('clusters by size')
				ax.set_ylabel('lipids per cluster')
				ax.set_xlabel('time (ns)')
				im = ax.imshow(canvas_traj,origin='lower',interpolation='nearest',aspect=aspect,extent=[time_start,time_end,0-0.5,max_sizes+1-0.5])
				ax.set_ylim((2-0.5,max_sizes+1-0.5))

				axins = inset_axes(ax,width="5%",height="100%",loc=3,
					bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
				#! cbar = plt.colorbar(im,cax=axins,orientation="vertical")
				# discrete color bar via https://stackoverflow.com/a/14779462/3313859 with slight mods
				cmap = mpl.cm.__dict__['jet']
				bounds = np.arange(1,max_counts+2)
				ticks = np.arange(1,max_counts+2)+0.5
				norm = mpl.colors.BoundaryNorm(bounds,cmap.N)
				cbar = mpl.colorbar.ColorbarBase(axins,cmap=cmap,norm=norm,
					spacing='proportional',ticks=ticks,boundaries=bounds,format='%1i')
				#! cbar.set_label('number\nof\nclusters',labelpad=-40,y=1.1,rotation=0)
				axins.set_title('observed\nclusters')
				axins.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
				axins.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
				axins.set_yticks([])

				times = np.linspace(time_start,time_end,nframes)
				ax = axes[layout['panels']['counts']]
				ax.set_title('lipids in clusters')
				ax.set_ylabel('number of lipids')
				ax.set_xlabel('time (ns)')
				base = np.zeros(nframes)
				#! reorder lipids hardcoded
				#! resnames_all_possible = work.metadata.variables.get('selectors',{})['resnames_lipid_chol'] 
				resnames_all_possible = ['DOPE']+['POPC', 'DOPC', 'DOPS', 'CHL1', 'PI2P', 'P35P', 'PIPU', 'PIPP', 'SAPI']
				reorder = np.array([resnames_all[i] for i in resnames_all_possible if i in resnames_all])
				totals_by_lipid = np.cumsum([totals_by_lipid[i] for i in reorder],axis=0)
				for ii,i in enumerate(totals_by_lipid):
					if ii==0: y2 = np.zeros(nframes)
					else: y2 = totals_by_lipid[ii-1]
					y1 = totals_by_lipid[ii]
					ax.fill_between(times,y1=y1,y2=y2,color=colors_these[reorder[ii]],zorder=3-ii)
				ax.set_aspect(np.ptp(ax.get_xlim())/np.ptp(ax.get_ylim()))

				# share y-axes
				#! axes[1].get_shared_y_axes().join(axes[1],axes[0])

				for ax in axes:
					ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
					ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')

				ax = axes[layout['panels']['counts']]
				legendspec = [
					dict(patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=colors_these[lnum]),
						name=resnames_all_r[lnum])
					for lnum in range(n_lipids)]
				patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
				legend = ax.legend(patches,labels,
					#! adjust the bbox to taste
					loc='upper center', bbox_to_anchor=(0.5,-0.15),
					ncol=n_lipids,fontsize=14)
				frame = legend.get_frame()
				frame.set_edgecolor('w')
				frame.set_facecolor('white')
				# via https://stackoverflow.com/a/4701285/3313859 
				#! this trick does not really shrink the axis however
				box = ax.get_position()
				ax.set_position([box.x0, box.y0 + box.height * 0.1,
					box.width, box.height * 0.9])

				zoom_figure(fig,zoom_fac)
				plt.subplots_adjust(wspace=0.45,hspace=0.35)
				picturesave('fig.cluster_crude.%s'%sn,work.plotdir,
					backup=False,version=True,meta={},extras=[legend],dpi=300,form='png')

			# preliminary strippplot option
			if 0:

				#import seaborn as sb
				#import pandas as pd

				# one timepoint is: canvas_traj_alt[:,0]
				whens = np.linspace(0,canvas_traj_alt.shape[1]-1,10).astype(int)
				this = canvas_traj_alt[:,0]
				raw = np.concatenate([np.ones(int(i*j))*i for i,j in zip(np.arange(len(this)),this)])
				raws = [np.concatenate([np.ones(int(i*j))*i for i,j in zip(np.arange(len(this)),this)])
					for this in [canvas_traj_alt[:,i] for i in whens]]
				# send the actual histograms to the df
				raw_prep = pd.DataFrame(dict([(ii,canvas_traj_alt[:,i]) for ii,i in enumerate(whens)]))
				dat_raw = np.concatenate(np.array([[(ii,j) for j in i] for ii,i in zip(range(len(raws)),raws)]))
				df = pd.DataFrame(data=dat_raw,columns=['x','y'])

				fig = plt.figure(figsize=(8,8))
				gs = gridspec.GridSpec(1,1)
				axes = [plt.subplot(g) for g in gs]
				ax = axes[0]
				sb.set(style="whitegrid")
				ax = sb.stripplot(x='x',y='y',data=df,size=3,jitter=0.2)
				#! ax = sb.violinplot(x='x',y='y',data=df)
				picturesave('TMP11-%s'%sn,work.plotdir,
					backup=False,version=True,meta={},extras=[legend],dpi=300,form='png')

			# above has one dot per lipid. below changes to one dot per cluster and uses a sizer
			# preliminary strippplot option
			if 0:

				import seaborn as sb
				import pandas as pd

				# one timepoint is: canvas_traj_alt[:,0]
				whens = np.linspace(0,canvas_traj_alt.shape[1]-1,10).astype(int)
				this = canvas_traj_alt[:,0]
				raw = np.concatenate([np.ones(int(i*j))*i for i,j in zip(np.arange(len(this)),this)])
				raws = [np.concatenate([np.ones(int(i*j))*i for i,j in zip(np.arange(len(this)),this)])
					for this in [canvas_traj_alt[:,i] for i in whens]]
				# send the actual histograms to the df
				raw_prep = pd.DataFrame(dict([(ii,canvas_traj_alt[:,i]) for ii,i in enumerate(whens)]))
				dat_raw = np.concatenate(np.array([[(ii,j) for j in i] for ii,i in zip(range(len(raws)),raws)]))
				# above dat_raw is one lipid per cluster size while raw_prep is the number of each cluster size

				df = pd.DataFrame(data=dat_raw,columns=['x','y'])

				# size mapper
				dot_s_min,dot_s_max = 1.,dat_raw[:,1].max()
				def sizer(i,small=0.5,big=10):
					return float((i-dot_s_min))/(dot_s_max-dot_s_min)*(big-small)+small
				sizes = np.array([sizer(i) for i in dat_raw[:,1]])

				fig = plt.figure(figsize=(8,8))
				gs = gridspec.GridSpec(1,1)
				axes = [plt.subplot(g) for g in gs]
				ax = axes[0]
				sb.set(style="whitegrid")
				# ax = sb.stripplot(x='x',y='y',data=df,size=3,jitter=0.2,s=sizes)
				ax = sb.stripplot(x='x',y='y',data=raw_prep,jitter=0.2,s=sizes)
				#! ax = sb.violinplot(x='x',y='y',data=df)
				picturesave('TMP12-%s'%sn,work.plotdir,
					backup=False,version=True,meta={},extras=[legend],dpi=300,form='png')

			# final cluster blob plot
			if 0: 

				from matplotlib.collections import PatchCollection

				class dotrow:
					def __init__(self,ax,size,center):
						self.ax = ax
						self.center = center
						self.size = size

				import seaborn as sb
				import pandas as pd
				sb.set(style="white")

				# one timepoint is: canvas_traj_alt[:,0]
				whens = np.linspace(0,canvas_traj_alt.shape[1]-1,10).astype(int)
				this = canvas_traj_alt[:,0]
				raw = np.concatenate([np.ones(int(i*j))*i for i,j in zip(np.arange(len(this)),this)])
				raws = [np.concatenate([np.ones(int(i*j))*i for i,j in zip(np.arange(len(this)),this)])
					for this in [canvas_traj_alt[:,i] for i in whens]]
				# send the actual histograms to the df
				raw_prep = pd.DataFrame(dict([(ii,canvas_traj_alt[:,i]) for ii,i in enumerate(whens)]))
				dat_raw = np.concatenate(np.array([[(ii,j) for j in i] for ii,i in zip(range(len(raws)),raws)]))
				# above dat_raw is one lipid per cluster size while raw_prep is the number of each cluster size
				df = pd.DataFrame(data=dat_raw,columns=['x','y'])

				# size mapper
				dot_s_min,dot_s_max = 1.,dat_raw[:,1].max()
				def sizer(i,small=0.5,big=10):
					return float((i-dot_s_min))/(dot_s_max-dot_s_min)*(big-small)+small
				sizes = np.array([sizer(i) for i in dat_raw[:,1]])

				fig = plt.figure(figsize=(8,8))
				gs = gridspec.GridSpec(1,1)
				axes = [plt.subplot(g) for g in gs]
				ax = axes[0]

				# reform thing
				source = np.array(raw_prep)
				width = 5.0
				centers = np.arange(1,1+2*width)*2*width
				width_max = width*2
				obs_max = source.max()
				obs_max_inds = np.transpose(np.where(obs_max==raw_prep))[0]
				scale_cluster_size = obs_max_inds[0]
				scale = obs_max/width_max/4.
				yvals = np.arange(source.shape[0])

				xlims = (0+width-2*scale,centers[-1]+width+2*scale)
				ylims = (0,70)
				stretch = (ylims[1]-ylims[0])/(xlims[1]-xlims[0])

				cmap = mpl.cm.__dict__['jet']
				colors = [cmap(float(ii)/(len(source.T)-1)) for ii in range(len(source.T))]
				def scale_calc(nlipids):
					base_size = (scale/2.) # radius for scale_cluster_size
					apl = base_size**2*np.pi/scale_cluster_size
					"""
					apl = base_size**2*pi/scale_cluster_size
					area = pi*r**2 = nlipids * apl
					radius = np.sqrt((nlipids*apl)/np.pi)
					"""
					return np.sqrt((nlipids*apl)/np.pi)
				# ii is the time index
				for ii,this in list(enumerate(source.T))[:]:
					# j is the number of clusters and jj is the indexcluster size index
					for jj,j in enumerate(this):
						radius = scale_calc(jj) # off by a factor of 2 but whatever
						if (ii==0 and jj==2) or (ii==9 and jj==61):
							print(radius)
						if j==0: continue
						center = centers[ii]
						#! ax.axvline(center,color='k',zorder=1)
						spots = np.arange(j)*radius*2
						spots += center-(spots.max()-spots.min())/2.
						patches = []
						for x in spots:
							detail = (x,yvals[jj])
							patch = mpl.patches.Circle(detail,radius=radius,lw=0)
							patches.append(patch)
						collection = PatchCollection(patches,alpha=1.0,color='k',zorder=10,lw=0)#,color=colors[ii]
						ax.add_collection(collection)
				ax.set_aspect('equal')
				ax.set_xlim(xlims)
				ax.set_ylim(ylims)
				times = np.linspace(time_start,time_end,nframes)
				ax.set_xticks(centers)
				ax.set_xticklabels(times[whens].astype(int)-1)
				ax.set_xlabel('time (ns)')
				ax.set_ylabel('$N_{lipids}$ per cluster')
				title = '%s, %s'%(work.metadata.meta[sn]['ptdins_label'],work.metadata.meta[sn]['ion_label'])
				ax.set_title(title)
				picturesave('fig.cluster_blob.%s'%sn,work.plotdir,
					backup=False,version=True,meta={},extras=[],dpi=300,form='png')

	# mark that we have completed one interation
	once_through = True
