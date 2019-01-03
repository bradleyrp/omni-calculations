#!/usr/bin/env python

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

@loader
def load():
	work = WorkSpace(analysis=True)
	data = work.plotload(plotname='ptdins_percolate')
	sns = work.sns()

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
		incoming = basic_compute_loop(compute_function,looper)
		post[sn] = dict(
			incoming=incoming,coords_subject=coords_subject,
			subject=subject,target=target,times=times,
			vecs=vecs,nframes=nframes)

def trim_graph_fat(points,edges):
	"""
	"""
	### need a fast reindex operation to select only the points in edges
	# currently pts_order and edges includes the entire set of lipids and cations while only some are 
	#   connected by bonds. hence we need a reducer step that reindexes. how to do a fast reindex?
	# fast reindexer method (disjoint case) 
	#! wut is disjoint case? explain

	# arguments are the original list of points and the edges which may not include all of them
	#!! originals = pts_order
	#!! edges = edges
	#!! we mark the points that are necessary for the left and right items in the edges list
	#! retained_l = np.zeros(len(originals)).astype(int)
	#! retained_r = np.zeros(len(originals)).astype(int)
	# we get unique points indices for the left and right side
	originals_u_l = np.unique(edges[:,0])
	originals_u_r = np.unique(edges[:,1])
	#!! find the edges that are in the unique l/r points list
	#! retained_l[originals_u_l[np.where(np.in1d(edges[:,0],originals_u_l))]] = 1
	#! retained_r[originals_u_r[np.where(np.in1d(edges[:,1],originals_u_r))]] = 1
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
	return points_reduced,edges_reduced

def compute_cluster_network(fr):
	"""
	...!!!
	"""
	# get the targets and subjects in their native index
	left,right = np.array([post[sn]['incoming'][fr][k] for k in ['targets','subjects']])
	# get residues for the target lipids
	left_resids = post[sn]['target'].resids[left]
	# filter by lipids in the top monolayer
	left_resids_inds_top = np.where(np.in1d(left_resids,resids_all[imono==mn_top]))[0]
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
	left_inds = ind_mono_to_mesh[mn_top]
	n_cations = len(post[sn]['subject'])
	right_inds = np.arange(n_cations)
	points_inds = np.arange(len(left_inds)+len(right_inds)).astype(int)
	lr_offset = len(left_inds)
	# create a mapping from left/right indices to points_inds indices
	remapper = {}
	remapper['left'] = dict(zip(left_inds,np.arange(len(left_inds)).astype(int)))
	remapper['right'] = dict(zip(right_inds,(np.arange(len(right_inds))+lr_offset).astype(int)))

	# reindex bonds for the main set of points
	bonds_reform = np.array([[remapper[j][i] for i in bonds_u.T[jj]] 
		for jj,j in enumerate(['left','right'])]).T

	# arguments to the cluster finder
	inds = points_inds
	edges = bonds_reform
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
	for k,v in cluster_inds.items():
		clusters[np.where(clusters==v)] = k
	# once we have  the clusters we need to map back to the points themselves
	pts_order = np.concatenate([
		#! limit the points when selecting the dat
		dat['%d.%d.points'%(mn,fr)][:len(left_inds)],
		post[sn]['coords_subject'][fr],])
	pts_order_bak = np.array(pts_order)

	vec = vecs[fr][:2]

	#! this section is slow: print('before opt')
	clusters_u = np.unique(clusters)
	for cc,cnum in enumerate(clusters_u):
		#! print('optimizing %d/%d'%(cc+1,len(clusters_u)))
		#! ignore loners
		if cnum==0: continue
		# get points in one cluster
		this_cluster = np.where(clusters==cnum)
		#! np.linalg.norm(pts_order[this_cluster][:,:2]-pts_order[this_cluster][:,:2].mean(axis=0))
		#! this = pts_order[this_cluster][:,:2]-pts_order[this_cluster][:,:2].mean(axis=0)
		#! np.sqrt(np.sum(this.T[0]**2+this.T[1]**2)) 
		#! or: np.linalg.norm(pts_order[this_cluster][:,:2]-pts_order[this_cluster][:,:2].mean(axis=0))
		#! def centroid(shifts)
		pts_this = pts_order[this_cluster][:,:2]
		arg0 = np.zeros(len(pts_this)*2)
		def func(arg):
			x = pts_this + arg.astype(int).reshape((2,-1)).T*vec
			#! return np.linalg.norm(x-x.mean(axis=0))
			return np.sum((x-x.mean(axis=0))**2)
		ans = scipy.optimize.minimize(fun=func,x0=(arg0,),method='Powell')
		shifts = ans.x.astype(int).reshape((2,-1))
		if (shifts!=0).any():
			print('moving %d'%cc)
			print(shifts)
			#! no change if you do: pts_order[this_cluster][:,:2] += shifts.T*vec
			pts_order[this_cluster,:2] += shifts.T*vec
	#! print('after opt')

	# collapse borders by finding all cations bound to at least two lipids to make a link
	#! hacking again. limit the points when selecting the dat
	#! pts_red = dat['%d.%d.points'%(mn,fr)][:len(left_inds)]
	#! changing to pts_order so things are fixed over PBCs
	pts_red = pts_order[:len(left_inds)]
	#! use the fact  that the previous edge list had cations on the right
	right_u = np.unique(edges[:,1])
	reduced_bonds = [j for j in [edges[np.where(edges[:,1]==i)[0],0] for i in right_u] if len(j)>1]

	edges_reduced_basic = np.concatenate([
		list(itertools.combinations(v,2))
		for v in reduced_bonds])
	#! pts_lipids = pts_red[np.unique(edges_reduced)][:,:2]

	# package
	points_reduced,edges_reduced = trim_graph_fat(pts_order,edges)
	points_reduced_lipids,edges_reduced_lipids = trim_graph_fat(pts_red,edges_reduced_basic)
	outgoing = dict(
		pts=points_reduced,edges=edges_reduced,
		points_reduced_lipids=points_reduced_lipids,edges_reduced_lipids=edges_reduced_lipids,
		p2=pts_red,e2=edges_reduced_basic)
	return outgoing

def proccer(sn):
	mn_top = 0

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
	fr = 200
	outgoing = compute_cluster_network(fr)
	return outgoing

if __name__=='__main__':

	if 'postpost' not in globals():
		postpost = {}
		for sn in ['membrane-v531','membrane-v532']:
			postpost[sn] = proccer(sn)

	# figure
	ax_pad = 0.2
	zoom_fac = 1.
	figscale = 10.0
	figprop = 1.
	zorder = {'pbc':0,'lines':1,'lipids':2,'cations':3}
	color_lipids = {'DOPE':'b','DOPS':'m','PI2P':'r','DOPC':'g','POPC':'g'}
	fig = plt.figure(figsize=(figscale*figprop,figscale))
	gs = gridspec.GridSpec(2,2)
	axes = [plt.subplot(g) for g in gs]

	# basic view with lines across PBCs is suppressed here
	if False:
		ax = axes[0]
		ax.set_title('lipids (black) bound to cations (red)')
		# plot the lipid points
		pts = dat['%d.%d.points'%(mn,fr)]
		pts_lipids = pts[np.unique(ind_mesh_to_mono[bonds_u[:,0]])][:,:2]
		ax.scatter(*pts_lipids.T,s=2,c='k')
		pts_cations = post[sn]['coords_subject'][fr][np.unique(bonds_u[:,1])][:,:2]
		ax.scatter(*pts_cations.T,s=2,c='r')
		# draw lines between all points
		lines = np.array([pts[ind_mesh_to_mono[bonds_u[:,0]]][:,:2],
			post[sn]['coords_subject'][fr][bonds_u[:,1]][:,:2]]).transpose((1,0,2))
		lc = mpl.collections.LineCollection(lines,color='k',lw=0.5,zorder=zorder['lines'])
		ax.add_collection(lc)		

	if False:
		ax = axes[0]
		ax.set_title('clusters after correcting PBCs')
		# randomized order
		color_inds = np.unique(clusters)
		np.random.shuffle(color_inds)
		color_hash = dict([(i,mpl.cm.__dict__['winter'](float(ii)/len(color_inds)))
			for ii,i in enumerate(color_inds)])
		# the conditional makes the cations red because they are sequentially last
		color = [color_hash[clusters[i]] 
			if i<len(left_inds) else 'r' 
			for i in inds if clusters[i]!=0]
		pts_order_filter = np.array([pts_order[ii,:2] for ii,i in enumerate(pts_order) if clusters[ii]!=0])
		ax.scatter(*pts_order_filter.T,color=color,s=10)
		lines = np.array([pts_order[edges[:,0]][:,:2],pts_order[edges[:,1]][:,:2]]).transpose((1,0,2))
		lc = mpl.collections.LineCollection(lines,color='k',lw=0.5,zorder=zorder['lines'])
		ax.add_collection(lc)		

	for snum,sn in enumerate(sns):

		outgoing = postpost[sn]

		ax = axes[0+2*snum]
		ax.set_title('lipid-cation-lipid bonds')
		# randomized order
		#! we have to plot dots for pts_lipids which are the unique cation-involved lipids
		#!   and then use pts_red to be indexed by edges_reduced
		pts_red,edges_reduced = outgoing['pts'],outgoing['edges']
		ax.scatter(*pts_red[:,:2].T,color='k',s=10)
		lines = np.array([pts_red[edges_reduced[:,0]][:,:2],
			pts_red[edges_reduced[:,1]][:,:2]]).transpose((1,0,2))
		lc = mpl.collections.LineCollection(lines,color='k',lw=0.5,zorder=zorder['lines'])
		ax.add_collection(lc)		

		ax = axes[1++2*snum]
		ax.set_title('?')
		# randomized order
		#! we have to plot dots for pts_lipids which are the unique cation-involved lipids
		#!   and then use pts_red to be indexed by edges_reduced
		pts_red,edges_reduced = outgoing['p2'],outgoing['e2']
		pts_red = outgoing['points_reduced_lipids']
		edges_reduced = outgoing['edges_reduced_lipids']
		ax.scatter(*pts_red[:,:2].T,color='k',s=10)
		lines = np.array([pts_red[edges_reduced[:,0]][:,:2],
			pts_red[edges_reduced[:,1]][:,:2]]).transpose((1,0,2))
		lc = mpl.collections.LineCollection(lines,color='k',lw=0.5,zorder=zorder['lines'])
		ax.add_collection(lc)		

	if False:
		for pspec,ax in zip(pspecs,axes):
			sn = sn=pspec['sn']
			pspec['func'](ax=ax,sn=sn)
	zoom_figure(fig,zoom_fac)
	plt.subplots_adjust(wspace=0.35,hspace=0.35)
	picturesave('TMP6',work.plotdir,backup=False,version=True,meta={},extras=[],dpi=300,form='png')

