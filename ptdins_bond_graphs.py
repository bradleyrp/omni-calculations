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

# CON-to-new this loader works when you bring in new data
#! @loader 
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

	if 1:
		if 'post' not in globals(): post = {}
		global frameslice
		frameslice = (770,790)
		
		# CON-to-new reformulate the post
		if 1:
			if 'post' not in globals():
				post2 = {}
				for sn in sns:
					data.set('ptdins_bond_graphs_calc')
					keys = [re.match('^0__(.+)$',i).group(1) for i in dat if re.match('^0__',i)]
					dat = data.this[sn]
					data.set('lipid_abstractor')
					dat_abs = data.this[sn]
					nframes = dat_abs['nframes']
					post2[sn] = dict(cluster_results=[{k:dat['%d__%s'%(fr,k)] for k in keys} 
						for fr in range(nframes)])
		# CON-to-new reformulate the post
		if 1:
			if 'data2' not in globals():
				data2 = work.plotload(plotname='ptdins_bond_graphs_calc')
				post_rgs = {}
				for snum,sn in enumerate(sns):
					data.set('lipid_abstractor')
					dat_abs = data.this[sn]
					data2.set('ptdins_bond_graphs_calc')
					dat = data2.this[sn]
					nframes = dat_abs['nframes']
					rgs = np.concatenate([dat['%d__rgs'%fr] for fr in range(nframes)])
					post_rgs[sn] = rgs
		# CON-to-new
		if 1:
			mn_top = 0
			cluster_contents ={}
			for sn in sns:
				imono = data.this[sn]['monolayer_indices']
				data.set('lipid_abstractor')
				resnames = data.this[sn]['resnames'][imono==mn_top]
				resnames_int = np.zeros(len(resnames)).astype(int)
				resname_to_int = dict([(j,i) for i,j in enumerate(np.unique(resnames))])
				for resname in resname_to_int:
					resnames_int[np.where(resnames==resname)] = resname_to_int[resname]
				rgs = post_rgs[sn]
				counts_resnames = np.zeros((len(rgs),len(resname_to_int)))
				n_resnames = np.unique(resnames)
				counter = 0
				for fr in range(nframes):
					clusters_this = data2.this[sn]['%d__clusters_lipids'%fr]
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

		# CON-to-new finding the largest cluster to see if something is weird
		if 1:
			cluster_ind_to_fr = {}
			for sn in sns:
				rgs = post_rgs[sn]
				n_clusters_per_frame = [np.unique(data2.this[sn]['%d__clusters_lipids'%fr]).shape[0]-1 
				for fr in range(nframes)]
				max_rg_ind = rgs.argmax()
				cluster_ind_to_fr[sn] = np.concatenate([np.ones(i)*ii 
					for ii,i in enumerate(n_clusters_per_frame)]).astype(int)
				# 775. adjust frameslice above to render it with the code below
				fr = cluster_ind_to_fr[sn][max_rg_ind]
		# CON-to-new rgs compute
		if 1:
			colors = {b'DOPE':'gray',b'DOPS':'blue',b'PI2P':'red'}
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

	if 1:
		
		colors = {b'DOPE':'gray',b'DOPS':'blue',b'PI2P':'red'}
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
			resnames = data.this[sn]['resnames'][imono==mn_top]
			colors_this = [colors[resnames[i]] for i in lipids_these]
			ax.scatter(*pts_these[:,:2].T,color=colors_this,s=10)
			lines = np.array([pts[edges[:,0]][:,:2],
				pts[edges[:,1]][:,:2]]).transpose((1,0,2))
			lc = mpl.collections.LineCollection(lines,color='k',lw=0.5,zorder=zorder['lines'])
			ax.add_collection(lc)		
		# plot formatting
		if False:
			for pspec,ax in zip(pspecs,axes):
				sn = sn=pspec['sn']
				pspec['func'](ax=ax,sn=sn)
		zoom_figure(fig,zoom_fac)
		plt.subplots_adjust(wspace=0.35,hspace=0.35)
		picturesave('TMP6',work.plotdir,backup=False,version=True,meta={},extras=[],dpi=300,form='png')

	if 1:

		# we have some kind of funky mapping where some of the points are not participating in edges so we reindex
		pts,edges = outgoing['points_reduced_lipids'],outgoing['edges_reduced_lipids']
		lipids_these = np.where(post[sn]['cluster_results'][fr]['clusters_lipids']>0)[0]
		reindex = np.zeros(len(pts)).astype(int)
		reindex[np.unique(outgoing['edges_reduced_lipids'])] = np.arange(len(lipids_these))
		# corresponding cluster index for this frame is
		cif = post[sn]['cluster_results'][fr]['clusters_lipids'][lipids_these]
		bonds = resnames[lipids_these][reindex[edges]]
		