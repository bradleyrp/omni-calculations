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

#! previously followed the computation of post
if False:

	class PeriodicMesh:
		def __init__(self,pairs,pairs_inds,vec):
			self.pairs = pairs
			self.pairs_inds = pairs_inds
			self.vec = vec
			# left and right points come from the pairs
			self.ptsl,self.ptsr = self.pairs.transpose((1,0,2))
			self.points = {}
			self.lines = {}
			self.inds = {}

			# make lipid-cation links
			self.make_links(name='lipid_cation',left=self.ptsl,right=self.ptsr,
				pairs=self.pairs,inds=self.pairs_inds)

			# perform the reduction
			left,right = self.pairs_inds.T
			# reduce collcetions
			hashed = dict([(l,list(right[np.where(left==l)[0]])) for l in np.unique(left)])
			hashed_coupled = dict([(k,v) for k,v in hashed.items() if len(v)>1])
			# itertools is sorted so we have no danger of repeats
			if not hashed_coupled: 
				self.lines['reduced'] = []
				return
			self.links = np.unique(np.concatenate([list(itertools.combinations(v,2)) 
				for v in hashed_coupled.values()]),axis=0)
			links_lines = dat['points'][fr][self.links][...,:2]
			self.lines['reduced'] = links_lines
			self.points['reduced'] = links_lines.transpose((1,0,2))

			# make lipid-lipid links
			self.make_links(name='lipid_lipid',
				left=self.points['reduced'][:,0],right=self.points['reduced'][:,1],
				pairs=self.lines['reduced'],inds=self.links)

			#! debug
			self.hashed = hashed
			self.hashed_coupled = hashed_coupled

		def make_links(self,name,left,right,pairs,inds):
			ptsl,ptsr,pairs = left,right,pairs
			vec = self.vec
			# convert points into lines
			n_links = len(pairs)
			# reformulate into pairs
			lines = np.array(pairs)
			# dealing with PBCs in a reactive way: find long lines and duplicate points
			deltas = lines.transpose((1,0,2)).ptp(axis=0)
			# identify boundary jumps
			jumps = np.transpose(np.where(deltas>vec[:2]/2.))
			extra_lines = []
			ptsl_extra = []
			ptsr_extra = []
			inds_extra = []
			# for each jump we identify the axis of the jump and then duplicate the point along that axis
			for ind,dim in jumps:
				# dim tells us the dimension that jumps, now we find direction, which
				#   explains the order in which we need to move the first or second 
				#   points in the pair either positive or negative directions
				way = 1*np.diff(lines[ind][:,dim])[0]>0
				shifts = np.zeros((2,2,2))
				shifts[0,0,dim] = 2*way-1
				shifts[1,1,dim] = 1-2*way
				for s in shifts:
					this = lines[ind]+s*vec[:2]
					inds_extra.append(inds[ind])
					extra_lines.append(this)
					ptsl_extra.append(this[0])
					ptsr_extra.append(this[1])
			extra_lines = np.array(extra_lines)
			inds_extra = np.array(inds_extra)
			if len(ptsl_extra)>0:
				ptsl = np.concatenate((ptsl[:,:2],ptsl_extra))
				ptsr = np.concatenate((ptsr[:,:2],ptsr_extra))
			lines_keep = np.ones(lines.shape[0],dtype=bool)
			lines_keep[jumps[:,0]] = False
			lines = lines[lines_keep]
			if len(extra_lines)>0:
				lines = np.concatenate((lines,extra_lines))
				inds = np.concatenate((inds,inds_extra))
			# package
			self.lines[name] = lines
			self.points[name] = np.array([ptsl,ptsr])
			self.inds[name] = inds

	data.set('lipid_abstractor')
	results = {}
	for sn in work.sns():

		incoming = post[sn]['incoming']
		coords_subject = post[sn]['coords_subject']

		dat = data.this[sn]
		nmol = len(dat['resids'])
		resid_to_index = dict(np.transpose([dat['resids'],np.arange(nmol)]))

		fr = 500 #! check out frame 100 it breaks your algo bro!

		# pairs are indexed by the selections
		pairs_raw = np.transpose([incoming[fr][i] for i in ['subjects','targets']])
		pairs_u = np.unique(pairs_raw,axis=0)
		# subject (cation) is indexed by the selection
		left = pairs_u[:,0]
		# target (lipid) is re-indexed by the center of  the lipid
		# the following mapping is from pairs_u (selection index) to target 
		#   (selection) to resids to index in the lipid_abstractor results, which
		#   is the centers of mass
		right = np.array(list(map(lambda x:resid_to_index[x],target[pairs_u[:,1]].resids)))
		# unique pairs
		left,right = np.unique(np.transpose([left,right]),axis=0).T
		# pairs by index for later reduction
		pairs_inds = np.array([left,right]).T
		sub_pts = coords_subject[fr][left]
		obj_pts = dat['points'][fr][right]
		pairs = np.array([sub_pts,obj_pts]).transpose((1,0,2))[...,:2]
		# formulate a mesh object
		vec = dat['vecs'][fr]
		pm = PeriodicMesh(pairs=pairs,pairs_inds=pairs_inds,vec=vec)
		results[sn] = pm

	### plot functions

	def birdseye(ax,sn):
		global results
		pm = results[sn]
		# plot pbc
		vec = pm.vec
		ax.plot(*np.array(
			[[0,0],[vec[0],0],[vec[0],vec[1]],[0,vec[1]]])[np.arange(0,5)%4].T,
			zorder=zorder['pbc'],lw=0.5,c='k')
		ax.set_xlim((vec[0]*(-1*ax_pad),vec[0]*(1.+ax_pad)))
		ax.set_ylim((vec[1]*(-1*ax_pad),vec[1]*(1.+ax_pad)))
		ax.set_aspect('equal')
		ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)
		ax.tick_params(axis='x',which='both',top=False,bottom=False,labelbottom=True)
		ax.set_xlabel('x (nm)')
		ax.set_ylabel('y (nm)')

	def subplot_lipid_cation(ax,sn):
		"""Plot the lipid-cation links."""
		global results,zorder
		pm = results[sn]
		colors = 'rk'
		for ii,pts in enumerate(pm.points['lipid_cation']):
			if False: # colors not working. probably have the outer leaflet in the data
				if ii==1: 
					color = [color_lipids[i] for i in dat['resnames'][pm.inds['lipid_cation'][:,1]].astype(str)]
				else: color = 'r'
			else: color = colors[ii]
			ax.scatter(pts[:,0],pts[:,1],s=2,color=color,zorder=zorder['cations'])
		# links between lipids
		lines = pm.lines['lipid_cation']
		if len(lines)!=None:
			lc = mpl.collections.LineCollection(lines,color='k',lw=0.5,zorder=zorder['lines'])
			ax.add_collection(lc)
		ax.set_title('bound cations, '+work.metadata.meta[sn]['ion_label'])

	def subplot_lipid_lipid(ax,sn):
		"""Plot the lipid-cation links."""
		global results,zorder
		pm = results[sn]
		#! for some reason the lipid_lipid points are only the ones crossing PBCs
		for name in ['reduced','lipid_lipid']:
			for ii,pts in enumerate(pm.points[name]):
				ax.scatter(pts[:,0],pts[:,1],s=2,color='k',zorder=zorder['cations'])
		# links between lipids
		lines = pm.lines['lipid_lipid']
		if len(lines)!=None:
			lc = mpl.collections.LineCollection(lines,color='k',lw=0.5,zorder=zorder['lines'])
			ax.add_collection(lc)
		ax.set_title('lipid-cation-lipid links, '+work.metadata.meta[sn]['ion_label'])

	# figure layout
	pspecs_base = [{'func':subplot_lipid_cation},{'func':subplot_lipid_lipid}]
	pspecs = [{'func':pbase['func'],'sn':sn}
		for pbase in pspecs_base
		for sn in sns]
		
	import matplotlib.gridspec as gridspec
	
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
	for pspec,ax in zip(pspecs,axes):
		sn = sn=pspec['sn']
		pspec['func'](ax=ax,sn=sn)
		birdseye(ax=ax,sn=sn)
	zoom_figure(fig,zoom_fac)
	plt.subplots_adjust(wspace=0.35,hspace=0.35)
	picturesave('TMP4',work.plotdir,backup=False,version=True,meta={},extras=[],dpi=300,form='png')

if __name__=='__main__':

	mn_top = 0
	fr = 200 #! check out frame 100 it breaks your algo bro!

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

	"""
	pseudocode 
		objective: identify clusters and wrap points
		steps
			1. identify clusters
				get a list of all points assigned an index of -1
				for each line
					if both points are -1 give it a new index
					else if either point is -1 give it the index of the other
					else merge the clusters to the lower value
			2. move points
	incoming data structure
		problem is that we want to have a set of non-redundant points so the clusters are well-defined
		the left points are lipids and the right ones are cations
		options:
			A. start with a canonical listing of all possible points and then make a mapping for each one
			B. at each frame collect the distinct points for each side and reindex them
		sticking with option A because that might be easier for later to track changes
		but are the points in the mesh object the same for each frame?
		it might actually be harder to formulate a single structure
	"""

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

	clusters_u = np.unique(clusters)
	for cc,cnum in enumerate(clusters_u):
		print('optimizing %d/%d'%(cc+1,len(clusters_u)))
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

	# collapse borders by finding all cations bound to at least two lipids to make a link
	#! hacking again. limit the points when selecting the dat
	#! pts_red = dat['%d.%d.points'%(mn,fr)][:len(left_inds)]
	#! changing to pts_order so things are fixed over PBCs
	pts_red = pts_order[:len(left_inds)]
	#! use the fact  that the previous edge list had cations on the right
	right_u = np.unique(edges[:,1])
	reduced_bonds = [j for j in [edges[np.where(edges[:,1]==i)[0],0] for i in right_u] if len(j)>1]

	edges_reduced = np.concatenate([
		list(itertools.combinations(v,2))
		for v in reduced_bonds])
	pts_lipids = pts_red[np.unique(edges_reduced)][:,:2]

	if 1:

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

		ax = axes[0]
		ax.set_title('clusters after correcting PBCs')

		# randomized order
		color_inds = np.unique(clusters)
		np.random.shuffle(color_inds)
		def get_color(i,n,name='winter'): 
			return mpl.cm.__dict__[name](float(color_inds[i])/n)
		# the conditional makes the cations red because they are sequentially last
		color = [get_color(clusters[i],n=np.unique(clusters).shape[0]) 
			if i<len(left_inds) else 'r' 
			for i in inds if clusters[i]!=0]
		pts_order_filter = np.array([pts_order[ii,:2] for ii,i in enumerate(pts_order) if clusters[ii]!=0])
		ax.scatter(*pts_order_filter.T,color=color,s=10)
		lines = np.array([pts_order[edges[:,0]][:,:2],pts_order[edges[:,1]][:,:2]]).transpose((1,0,2))
		lc = mpl.collections.LineCollection(lines,color='k',lw=0.5,zorder=zorder['lines'])
		ax.add_collection(lc)		

		ax = axes[1]
		ax.set_title('lipid-cation-lipid bonds')

		# randomized order
		#! we have to plot dots for pts_lipids which are the unique cation-involved lipids
		#!   and then use pts_red to be indexed by edges_reduced
		ax.scatter(*pts_lipids.T,color='k',s=10)
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

