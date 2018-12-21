
#!/usr/bin/env python

"""
Contact code derived from hydrogen_bonding and/or salt_bridges.
Note that the hydrogen bonding and salt bridges code explicitly codes for three-party interactions of 
a highly specific nature, whereas this code generalizes this to nearby contacts between two elements
with (1) the option to specify that the subject is e.g. the protein residues (instead of checking all 
contacts in the system, which would be costly) and (2) the same exact data structure as the two progenitor
codes so that later these contacts can be easily summed.
"""

import os,sys,glob,re,time
import numpy as np
import scipy
import scipy.spatial
import MDAnalysis
import numpy as np

from ortho import str_types,status
from omni import get_automacs,basic_compute_loop,uniquify
from omni import vecangle,vecnorm,HookHandler

# define the columns for a row in the master dataset
defn_rowspec = [
	'subject_resname','subject_resid','subject_atom',
	'target_resname','target_resid','target_atom']
defn_bonds_to_reduced = np.array([0,1,3,4])
defn_rowspec_reduced = [defn_rowspec[i] for i in defn_bonds_to_reduced]

def contacts_framewise(fr,target,subject,vec,**kwargs):
	"""
	Compute close contacts using global subject and target coordinates. 
	Called by the contacts function.
	"""
	debug = kwargs.get('debug',False)
	distance_cutoff = kwargs['distance_cutoff']
	lenscale = kwargs.get('lenscale',10.0)
	# ensure that the points are inside the box
	boxstuff = lambda pts,vec : pts-np.floor(pts/vec)*vec
	# convert back to advanced indexing
	aind = lambda x : tuple(x.T)
	# set the foreground and background points
	pts_back_unstuffed = np.array(target)
	pts_fore_unstuffed = np.array(subject)
	pts_back = boxstuff(pts_back_unstuffed,vec)
	pts_fore = boxstuff(pts_fore_unstuffed,vec)
	"""
	historical note: previously the boxsize to cKDTree had to be twice as long 
	as the regular box vectors for an unknown reason. this bug was eventually
	fixed and the cKDTree with a "torus" torus norm works great. the original
	implementation also implemented a foolish method for boxstuff:
	  lambda pts,vec : pts-(pts>vec)*vec+(pts<np.array([0.,0.,0.]))*vec
	which did not account for points that strayed by *two* box vectors!
	however this is now corrected so we dropped the exception
	  except: return {'subjects':np.array([]),'targets':np.array([])}
	note that the exception required a ton of valid-frames tabulation downstream
	of this code, particularly when synchronizing this calculation with others
	note that future checks of boxstuff results with the following line
	  np.any(pts_back>vec) or np.any(pts_back<0)
	showed that some incoming points were still outside the box for some reason
	possible related to precision. hence ERROR REMAINS for some frames!
	"""
	try: tree = scipy.spatial.ckdtree.cKDTree(pts_back,boxsize=vec)
	except: return {'subjects':np.array([]),'targets':np.array([])}
	close,nns = tree.query(pts_fore,k=10,
		distance_upper_bound=distance_cutoff/lenscale)
	# index pairs within the cutoff distance
	close_pairs = np.transpose(np.where(np.all((
		close<distance_cutoff/lenscale,close>0),axis=0)))
	# list of close donors
	close_targets = nns[aind(close_pairs)]
	close_subjects = close_pairs[:,0]
	return {'subjects':close_subjects,'targets':close_targets}

def count_reduced_contact(resid,resname_set,bonds,obs,bonds_red,mode='full'):
	"""
	Tally up the contacts based on a filter over subject resid and lipid residue name.
	We have two modes: explicit retains full residue/atom specificify while reduced only tells you if a 
	bond between two residues exists.
	Note that this is the kernel of the part of the contact code that counts the unique bonds identified 
	in the beginning of the contact code.
	"""
	if mode=='full':
		# filter the observations by the protein residue (subject_resid) and target resname
		#   providing a result
		which = np.where(np.all((bonds[:,
			defn_rowspec.index('subject_resid')].astype(int)==
			resid,np.in1d(bonds[:,defn_rowspec.index('target_resname')],resname_set)),axis=0))
		result = obs.T[which].sum(axis=0)
	elif mode=='reduced':
		# the explicit result has rows over specific bonds (in reduced form, minus atom names)
		#! deprecated: explicit = obs.T[np.where(np.all((bonds_red[:,
		#!   rowspec_red.index('subject_resid')].astype(int)==resid,np.in1d(
		#!     bonds_red[:,rowspec_red.index('target_resname')],resname_set)),axis=0))]
		# we sum over the rows and only test for non-zero, since we do not care which specific bond is 
		#   formed, hence returning whether or not there was a bond between protein residue and lipid
		#! deprecated, incorrect: result = (explicit.sum(axis=0)>1)*1
		which = np.where(np.all((bonds_red[:,
			defn_rowspec_reduced.index('subject_resid')].astype(int)==
			resid,np.in1d(bonds_red[:,defn_rowspec_reduced.index('target_resname')],
				resname_set)),axis=0))
		idx,counts = uniquify(bonds_red[which].astype(str))
		result = obs.T[which][idx].sum(axis=0)
	else: raise Exception('mode is either full or reduced')
	return result

def count_reduced_contact_reduced(**kwargs): 
	"""Decorator for the compute loop."""
	return count_reduced_contact(mode='reduced',**kwargs)
	
def contacts(**kwargs):
	"""
	Identify, catalog, and count contacts in a simulation.
	Note that this code was developed to mimic the data structures used by hydrogen bonding and salt bridging
	codes, and stores the various contacts up to a very high level of specificity (i.e. the specific residue
	and atom names).
	"""
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	debug = kwargs.get('debug',False)
	run_parallel = kwargs.get('run_parallel',True)
	# distance cutoff stays in angstroms until the compute function
	distance_cutoff = calc['specs']['cutoff']

	# get coordinates from the hook
	default_fetch_coordinates = {
		'import_target':'codes/fetch_coordinates.py',
		'function':'fetch_coordinates_via_metadata'}
	fetch_coordinates = kwargs.get(
		'fetch_coordinates',default_fetch_coordinates)
	# hooks get the same arguments as compute functions
	coordinates = HookHandler(meta=kwargs,**fetch_coordinates).solve
	nframes = coordinates['nframes']
	vecs = coordinates['vecs']
	subject = coordinates['uni_selections']['subject']
	targets = coordinates['uni_selections']['target']
	times = coordinates['times']

	# save topology for later
	uni = coordinates['uni']
	_,idx,counts = np.unique(uni.residues.resnames,return_index=True,return_counts=True)
	resnames = uni.residues.resnames[np.sort(idx)]
	resnames_master = np.array(resnames)
	rescounts = counts[np.argsort(idx)]

	# parallel compute
	coords_target = coordinates['coordinates_by_selection']['target']
	coords_subject = coordinates['coordinates_by_selection']['subject']
	compute_function = contacts_framewise
	out_args = {'distance_cutoff':distance_cutoff}
	looper = [dict(fr=fr,target=coords_target[fr],subject=coords_subject[fr],
		vec=vecs[fr],**out_args) for fr in range(nframes)]
	incoming = basic_compute_loop(compute_function,looper,run_parallel=run_parallel)
	# data reduction, get valid frames
	valid_frames = np.where([len(i['subjects'])>0 for i in incoming])[0]
	if all([len(incoming[i]['subjects'])==0 for i in valid_frames]):
		raise Exception('no contacts. maybe your cutoff is too short: %s'%distance_cutoff)
	obs_by_frames = np.array([len(incoming[i]['subjects']) for i in valid_frames]).astype(int)
	# concatenate the donor/acceptor indices across all frames
	subject_cat = np.concatenate([incoming[i]['subjects'] for i in valid_frames]).astype(int)
	target_cat = np.concatenate([incoming[i]['targets'] for i in valid_frames]).astype(int)

	start_time = time.time()
	# tabulate each bond observation
	tabulation = np.transpose((subject.resnames[subject_cat],subject.resids[subject_cat],
		subject.names[subject_cat],targets.resnames[target_cat],targets.resids[target_cat],
		targets.names[target_cat],))
	status('stopwatch: %.1fs'%(time.time()-start_time),tag='compute')

	idx,counts = uniquify(tabulation.astype(str))
	bonds_catalog = tabulation[idx]

	start_time = time.time()
	# preallocate bond counts per frame
	counts_per_frame = np.zeros((len(valid_frames),len(idx)))
	# hash the binds over the indices
	bonds_to_idx = dict([(tuple(b),bb) for bb,b in enumerate(bonds_catalog)])
	frame_lims = np.concatenate(([0],np.cumsum(obs_by_frames)))
	for fr,i in enumerate(frame_lims[:-1]):
		status('counting observations per frame',i=fr,looplen=len(valid_frames),
			tag='compute',start=start_time)
		obs_this = tabulation[frame_lims[fr]:frame_lims[fr+1]]
		counts_per_frame[fr][np.array([bonds_to_idx[tuple(o)] for o in obs_this])] += 1
	status('stopwatch: %.1fs'%(time.time()-start_time),tag='compute')
	status('done heavy lifting',tag='compute')
	# note the size of the outgoing data. we could shrink this by discarding atom names
	status('observation array for cutoff %.1f is %.1fMB'%(
		distance_cutoff,sys.getsizeof(counts_per_frame)/10**6.),tag='note')

	# package the dataset
	result,attrs = {},{}
	# everything is indexed by idx
	result['bonds'] = bonds_catalog
	result['observations'] = counts_per_frame
	result['valid_frames'] = valid_frames
	result['nframes'] = np.array(nframes)
	result['resnames'] = resnames_master
	result['subject_residues_resnames'] = subject.residues.resnames
	result['target_residues_resnames'] = targets.residues.resnames
	result['subject_residues_resids'] = subject.residues.resids
	result['target_residues_resids'] = targets.residues.resids
	result['nmols'] = rescounts
	result['times'] = np.array(times)
	# save the rowspec so we know up to unpack the data correctly
	attrs['defn_rowspec'] = defn_rowspec
	attrs['defn_rowspec_reduced'] = defn_rowspec

	#! this needs refactored
	if False:

		# some basic post-processing common to many of the plots
		bonds,obs = bonds_catalog,counts_per_frame
		# post: generate timewise trajectories for the number of contacts between protein residues and lipids
		# methodology note: in a basic version of this calculation we simply count all of the bonds between 
		#   any lipid-protein residue pair. this means that being more close to a lipid might result in more 
		#   contacts and hence generates a higher score. hence we have two versions of the calculation. one 
		#   counts the total number of contacts, and the other discards atom information and scores contacts 
		#   with a maximum of one per protein residue-lipid pair. this calculation does both
		#! need to check for atom-name resolution otherwise this is moot
		resids = result['subject_residues_resids']
		lipid_resnames = np.unique(bonds[:,defn_rowspec.index('target_resname')])
		resname_combos = [(r,np.array([r])) for r in lipid_resnames]+[('all lipids',np.array(lipid_resnames))]
		# reduce the data for the modified count described above
		bonds_red = bonds[:,defn_bonds_to_reduced]
		looper = [{'resid':resid,'resname_set':resname_set,'bonds':bonds,'obs':obs,'bonds_red':bonds_red} 
			for resid in resids for resname_name,resname_set in resname_combos]
		compute_function = count_reduced_contact
		incoming = basic_compute_loop(compute_function,looper,run_parallel=run_parallel)
		# package this as a list of resid/resname pairs and the counts for them
		result['pairs_resid_resname'] = np.array([(resid,resname_name) 
			for resid in resids for resname_name,resname_set in resname_combos]).astype(str)
		result['counts_resid_resname'] = np.array(incoming)
		compute_function = count_reduced_contact_reduced
		incoming = basic_compute_loop(compute_function,looper,run_parallel=run_parallel)
		result['counts_resid_resname_singleton'] = np.array(incoming)
	
	return result,attrs
