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
import itertools

@loader
def load():
	work = WorkSpace(analysis=True)
	data = work.plotload(plotname='ptdins_percolate')
	sns = work.sns()

if __name__=='__main__':

	#! development
	if 'uni' not in globals():

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
			post[sn] = {'incoming':incoming,'coords_subject':coords_subject}

	class PeriodicMesh:
		def __init__(self,pairs,vec,pairs_inds=None):
			"""
			We retain the distinction between two sets of points for aesthetics.
			"""
			self.left = pairs[:,0]
			self.right = pairs[:,1]
			self.points_pairs = pairs
			self.vec = vec
			self.left_i = pairs_inds[:,0]
			self.right_i = pairs_inds[:,1]

		def pbc_lines(self):
			"""..."""
			sub_pts,obj_pts = self.left,self.right
			vec = self.vec
			# convert points into lines
			n_links = len(sub_pts)
			# reformulate into pairs
			lines = np.array([sub_pts,obj_pts]).transpose((1,0,2))[...,:2]
			# dealing with PBCs in a reactive way: find long lines and duplicate points
			deltas = lines.transpose((1,0,2)).ptp(axis=0)
			# identify boundary jumps
			jumps = np.transpose(np.where(deltas>vec[:2]/2.))
			extra_lines = []
			sub_pts_extra = []
			obj_pts_extra = []
			# for each jump we identify the axis of the jump and then duplicate the point along that axis
			for ind,dim in jumps:
				# dim tells us the dimension that jumps, now we find direction, which
				#   explains the order in which we need to move the first or second 
				#   points in the pair either positive or negative directions
				way = 1*np.diff(lines[ind][:,dim])[0]>0
				shifts = np.zeros((2,2,2))
				#shifts[0,way,dim] = 1-way
				#shifts[1,way,dim] = way
				shifts[0,0,dim] = 2*way-1
				shifts[1,1,dim] = 1-2*way
				for s in shifts:
					this = lines[ind]+s*vec[:2]
					extra_lines.append(this)
					sub_pts_extra.append(this[0])
					obj_pts_extra.append(this[1])
			extra_lines = np.array(extra_lines)
			if len(sub_pts_extra)>0:
				sub_pts = np.concatenate((sub_pts[:,:2],sub_pts_extra))
				obj_pts = np.concatenate((obj_pts[:,:2],obj_pts_extra))
			lines_keep = np.ones(lines.shape[0],dtype=bool)
			lines_keep[jumps[:,0]] = False
			lines = lines[lines_keep]
			if len(extra_lines)>0:
				lines = np.concatenate((lines,extra_lines))
			self.lines = lines
			self.pts_left = sub_pts
			self.pts_right = obj_pts

		def reduce_lines(self):
			"""..."""
			right,left = self.right_i,self.left_i
			# reduce collcetions
			hashed = dict([(l,list(right[np.where(left==l)[0]])) for l in np.unique(left)])
			hashed_coupled = dict([(k,v) for k,v in hashed.items() if len(v)>1])
			# itertools is sorted so we have no danger of repeats
			links = np.unique(np.concatenate([list(itertools.combinations(v,2)) 
				for v in hashed_coupled.values()]),axis=0)
			links_lines = dat['points'][fr][links][...,:2]
			self.lines2 = links_lines

	class PeriodicMesh:
		def __init__(self,pairs,pairs_inds,vec):
			self.pairs = pairs
			self.pairs_inds = pairs_inds
			self.vec = vec
			# left and right points come from the pairs
			self.ptsl,self.ptsr = self.pairs.transpose((1,0,2))
			self.points = {}
			self.lines = {}

			# make lipid-cation links
			self.make_links(name='lipid_cation',left=self.ptsl,right=self.ptsr,pairs=self.pairs)

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
				pairs=self.lines['reduced'])

		def make_links(self,name,left,right,pairs):
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
					extra_lines.append(this)
					ptsl_extra.append(this[0])
					ptsr_extra.append(this[1])
			extra_lines = np.array(extra_lines)
			if len(ptsl_extra)>0:
				ptsl = np.concatenate((ptsl[:,:2],ptsl_extra))
				ptsr = np.concatenate((ptsr[:,:2],ptsr_extra))
			lines_keep = np.ones(lines.shape[0],dtype=bool)
			lines_keep[jumps[:,0]] = False
			lines = lines[lines_keep]
			if len(extra_lines)>0:
				lines = np.concatenate((lines,extra_lines))
			# package
			self.lines[name] = lines
			self.points[name] = np.array([ptsl,ptsr])

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
		self = pm #! debug

	if True:

		ax_pad = 0.2
		zoom_fac = 1.5
		zorder = {'pbc':0,'lines':1,'lipids':2,'cations':3}
		fig = plt.figure()
		for pnum,(sn,result) in enumerate(results.items()):
			#! globals().update(**result)
			pm = result
			ax = fig.add_subplot([121,122][pnum],aspect='equal')

			#! good for lipid_cation
			if False:
				for color,pts in zip(['b','k'],pm.points['lipid_cation']):
					ax.scatter(pts[:,0],pts[:,1],s=2,color=color,zorder=zorder['cations'])

			for name in ['lipid_cation','reduced','lipid_lipid'][-1:]:
				for pts in pm.points['lipid_lipid']:
					ax.scatter(pts[:,0],pts[:,1],s=2,color='k',zorder=zorder['cations'])
				for pts in pm.points['reduced']:
					ax.scatter(pts[:,0],pts[:,1],s=2,color='k',zorder=zorder['cations'])
				# links between lipids
				lines = pm.lines[name]
				if len(lines)!=None:
					lc = mpl.collections.LineCollection(lines,color='b',lw=0.75,zorder=zorder['lines'])
					ax.add_collection(lc)

			# plot pbc
			vec = pm.vec
			ax.plot(*np.array(
				[[0,0],[vec[0],0],[vec[0],vec[1]],[0,vec[1]]])[np.arange(0,5)%4].T,
				zorder=zorder['pbc'],lw=0.5,c='k')

			ax.set_xlim((vec[0]*(-1*ax_pad),vec[0]*(1.+ax_pad)))
			ax.set_ylim((vec[1]*(-1*ax_pad),vec[1]*(1.+ax_pad)))
			ax.set_title(work.metadata.meta[sn]['ion_label'])

		zoom_figure(fig,zoom_fac)
		picturesave('TMP3',work.plotdir,backup=False,version=True,meta={},extras=[],dpi=300,form='png')

		if False:

			# reduce collections
			hashed = dict([(l,list(right[np.where(left==l)[0]])) for l in np.unique(left)])
			hashed_coupled = dict([(k,v) for k,v in hashed.items() if len(v)>1])
			# itertools is sorted so we have no danger of repeats
			links = np.unique(np.concatenate([list(itertools.combinations(v,2)) 
				for v in hashed_coupled.values()]),axis=0)
			links_lines = dat['points'][fr][links][...,:2]

			# package
			results[sn] = dict(sub_pts=sub_pts,obj_pts=obj_pts,lines=lines,links_lines=links_lines,vec=vec)

		if False:

			zorder = {'pbc':0,'lines':1,'lipids':2,'cations':3}
			fig = plt.figure()
			for pnum,(sn,result) in enumerate(results.items()):
				globals().update(**result)
				ax = fig.add_subplot([121,122][pnum],aspect='equal')

				ax.scatter(sub_pts[:,0],sub_pts[:,1],s=2,color='b',zorder=zorder['cations'])
				ax.scatter(obj_pts[:,0],obj_pts[:,1],s=2,color='k',zorder=zorder['lipids'])

				# lines between lipids and cations
				lc = mpl.collections.LineCollection(lines,color='k',lw=0.5,zorder=zorder['lines'])
				ax.add_collection(lc)

				# links between lipids
				lc = mpl.collections.LineCollection(links_lines,color='r',lw=0.75,zorder=zorder['lines'])
				ax.add_collection(lc)

				# plot pbc
				ax.plot(*np.array(
					[[0,0],[vec[0],0],[vec[0],vec[1]],[0,vec[1]]])[np.arange(0,5)%4].T,
					zorder=zorder['pbc'],lw=0.5,c='k')

			picturesave('TMP2',work.plotdir,backup=False,version=True,meta={},extras=[],dpi=300,form='png')

"""
moment of inertia
	radius of gyration
	square of distance
	(x-x0)**2

	moment of inertia

	not overall
		do radius of gyration tensor in 2D ixx iyy ixy
		then diagonolize
		lambda x and y tells you anisotropy

	ryan thinks still no connectivity

"""