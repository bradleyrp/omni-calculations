#!/usr/bin/env python

import re,time
from ortho import Handler,requires_python,status
from omni import get_automacs

class MultiCoordinates:
	"""
	Simple class that wraps MDAnalysis for getting coordinates
	and returning selections. Used as a hook to different functions so that
	this code is modular and interchangeable.
	"""
	def __init__(self,structure,trajectory,**kwargs):
		# all remaining keyword arguments are selections
		self.structure = structure
		self.trajectory = trajectory
		self.selections = kwargs
	@property
	def result(self):
		self.load()
		return self._result
	@requires_python('MDAnalysis')
	def load(self):
		"""Load a trajectory from MDAnalysis."""
		data = {}
		import MDAnalysis
		uni = MDAnalysis.Universe(self.structure,self.trajectory)
		lenscale = 10.0 # MDAnalysis uses Angstroms but we prefer nm
		nframes = len(uni.trajectory)
		coords = dict([(name,[]) for name,sel in self.selections.items()])
		selections_mda = dict([(name,uni.select_atoms(sel)) for name,sel in self.selections.items()])
		# prepare coordinates for each frame
		st = time.time()
		vecs,times = [],[]
		# purposefully profligate with the memory so this goes quickly
		for fr in range(nframes):
			status('caching coordinates',tag='compute',i=fr,looplen=nframes,start=st)	
			uni.trajectory[fr]
			times.append(uni.trajectory.time)
			vecs.append(uni.dimensions[:3]/lenscale)
			for name in coords: 
				coords[name].append(selections_mda[name].positions/lenscale)
		status('completed caching in %.1f minutes'%((time.time()-st)/60.),tag='status')
		data['uni'] = uni
		data['uni_selections'] = selections_mda
		data['nframes'] = nframes
		data['coordinates_by_selection'] = coords
		data['selections'] = self.selections
		data['vecs'] = vecs
		data['times'] = times
		self._result = data
		return data

def get_lipid_resnames(work):
	"""Use automacs and lipid definitions to return MDAnalysis selection for lipids."""
	mod = get_automacs()
	ff_name = work.metadata.variables.get('force_field',None)
	if not ff_name: raise Exception('we must be very careful with the residue naming. '
		'you must add `force_field` to the `variables` dictionary in your metadata to continue.')
	mod['state'].force_field = 'charmm'
	Landscape = mod['Landscape']
	land = Landscape(ff=ff_name)
	target_resnames = land.objects_by_category('lipid')
	return target_resnames

class PredefinedSelectionMDA:
	def __init__(self,name,work): 
		self.name = name
		self.work = work
	@property
	def result(self):
		return getattr(self,re.sub(' ','_',self.name))()
	def lipid_heavy(self):
		"""Return selection of lipid residue names."""
		target_resnames = get_lipid_resnames(self.work)
		subject_selection = '(%s) and not name H*'%' or '.join(['resname %s'%i for i in target_resnames])
		return subject_selection
	def protein_heavy(self): 
		"""MDAnalysis selection for protein with no hydrogens."""
		return "protein and not name H*"

def fetch_coordinates_via_metadata(**kwargs):
	"""
	Hook for getting coordinates of a subject and object specified 
	in the calc specs.
	"""
	kwargs_compute = kwargs.pop('kwargs',{})
	structure = kwargs_compute.get('structure',None)
	if not structure: structure = kwargs_compute.get('grofile',None)
	trajectory = kwargs_compute.get('trajectory',None)
	if not trajectory: trajectory = kwargs.get('trajectory',None)
	subject = kwargs_compute.get('calc',{}).get('specs',{}).get('subject',None)
	obj = kwargs_compute.get('calc',{}).get('specs',{}).get('target',None)
	if not subject or not obj:
		import ipdb;ipdb.set_trace()
		raise Exception('subject or target is empty: %s, %s'%(subject,obj))
	selections = {'subject':subject,'target':obj}
	# special selections here
	for key,this in selections.items():
		if isinstance(this,dict) and set(this.keys())=={'predefined'}:
			selections[key] = PredefinedSelectionMDA(
				this['predefined'],work=kwargs_compute['workspace']).result
	mc = MultiCoordinates(
		structure=structure,trajectory=trajectory,
		**selections)
	return mc.load()

#! the above works in contacts.py and uses the metadata from kwargs_compute
#! the below works more generally

def fetch_coordinates_subject_target(**kwargs):
	"""
	Hook for getting coordinates of a subject and object specified 
	in the calc specs.
	"""
	kwargs_this = kwargs['kwargs'] #! obviously foolish; needs fixed
	structure = kwargs_this.get('structure',None)
	trajectory = kwargs_this.get('trajectory',None)
	subject = kwargs_this.get('subject',None)
	obj = kwargs_this.get('object',None)
	selections = {'subject':subject,'target':obj}
	# special selections here
	#! special selections is repetitive with fetch_coordinates_subject_target
	for key,this in selections.items():
		if isinstance(this,dict) and set(this.keys())=={'predefined'}:
			selections[key] = PredefinedSelectionMDA(
				this['predefined'],work=kwargs_this['workspace']).result
	mc = MultiCoordinates(
		structure=structure,trajectory=trajectory,
		**selections)
	return mc.load()
