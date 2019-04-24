#!/usr/bin/env python

"""
IMPORT SAMANEH'S DATA!

previously used the omnicalc that includes maps.py and now updated for the one that uses structs.py
"""

import os,sys,glob,re
from ortho import status
import numpy as np

def import_membrane_mesh(**kwargs):
	"""
	Import membrane XYZ data and send it to a calculation that mimics `undulations` for Samaneh's data.
	"""
	sn = kwargs.pop('sn',None)
	calc = kwargs.pop('calc',None)
	work = kwargs.pop('work',None)
	if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
	#---location data can be found in the slices dictionary
	#---! note that the slice name is hard-coded here: "current"
	#!!! the work.slices is now a SliceMeta which is not subscriptable
	if not isinstance(work.slices,dict): work.slices = work.slices.__dict__['raw']
	location = work.slices[sn]['readymade_meso_v1']['current']
	#---get all of the filenames
	fns = sorted(glob.glob(os.path.join(location['path'],location['directory'],
		location['membrane_xyz_glob'])))
	#---! save timestamps here if desired. ensure synchronicity with the nanogel inputs
	#---read each file
	nframes = len(fns)
	points = []
	for fnum,fn in enumerate(fns):
		status('reading %s'%os.path.basename(fn),i=fnum,looplen=nframes,tag='load')
		with open(fn) as fp: text = fp.read()
		lines = text.splitlines()
		#---first line is metadata
		#---! should we drop the first line?
		topline = lines[0]
		#---the frame is xyz in columns plus fourth column
		frame = np.array([[float(j) for j in line.split()] for line in lines[1:]])
		points.append(frame)
	return points

def import_nanogel_positions_standard(**kwargs):
	"""
	Import nanogel data and send it to a calculation that mimics `protein_abstractor` for Samaneh's data.
	"""
	sn = kwargs.pop('sn',None)
	calc = kwargs.pop('calc',None)
	work = kwargs.pop('work',None)
	if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
	#---location data can be found in the slices dictionary
	#---! note that the slice name is hard-coded here: "current"
	#! fixing the SliceMeta object
	if not isinstance(work.slices,dict): work.slices = work.slices.__dict__['raw']
	location = work.slices[sn]['readymade_meso_v1']['current']
	with open(os.path.join(location['path'],location['directory'],location['nanogel_dat'])) as fp:
		text = fp.read()
	#---nanogel is saved with the step number not the frame number
	#! retiring this step_to_frame = lambda x: x/1000000
	step_to_frame = lambda x: x
	regex_frame = r'(\d+)\n(.*?)(?=\n\d+\n|\Z)'
	frames = re.findall(regex_frame,text,flags=re.M+re.DOTALL)
	framenos,points = [],[]
	for fnum,frame in enumerate(frames):
		status('reading nanogel frame',i=fnum,looplen=len(frames),tag='load')
		framenos.append(step_to_frame(int(frame[0])))
		ixyz = np.array([[float(j) for j in i.split()] for i in frame[1].splitlines()])
		if not np.all(ixyz[:,0].astype(int)==np.arange(1,len(ixyz)+1)):
			raise Exception('indexing problem in the nanogel')
		points.append(ixyz[:,1:])
	return {'framenos':framenos,'points':np.array(points)}

def import_nanogel_positions_rigid(**kwargs):
	"""
	Import nanogel data and send it to a calculation that mimics `protein_abstractor` for Samaneh's data.
	"""
	sn = kwargs.pop('sn',None)
	calc = kwargs.pop('calc',None)
	work = kwargs.pop('work',None)
	if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
	#---location data can be found in the slices dictionary
	#---! note that the slice name is hard-coded here: "current"
	#! fixing the SliceMeta object
	if not isinstance(work.slices,dict): work.slices = work.slices.__dict__['raw']
	location = work.slices[sn]['readymade_meso_v1']['current']
	with open(os.path.join(location['path'],location['directory'],location['nanogel_dat'])) as fp:
		text = fp.read()
	#---nanogel is saved with the step number not the frame number
	regex_frame = r'^(\d+)\s+(\d+)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(.*?)(?:\Z|\n)'
	frames = re.findall(regex_frame,text,flags=re.M+re.DOTALL)
	framenos,points = [],[]
	for fnum,frame in enumerate(frames):
		status('reading nanogel frame',i=fnum,looplen=len(frames),tag='load')
		#! record the raw frame here instead of dividing
		#! note that we will have some lack of synchronicity but that cannot be avoided for this dataset
		framenos.append(fnum)
		points.append([np.array(frame[2:5]).astype(float)])
	return {'framenos':framenos,'points':np.array(points)}

def import_nanogel_positions(**kwargs):
	"""
	Import nanogel data and send it to a calculation that mimics `protein_abstractor` for Samaneh's data.
	"""
	sn = kwargs['sn']
	work = kwargs['work']
	alt_call = work.meta.get(sn,{}).get('import_nanogel_alt',None)
	if alt_call: return globals()[alt_call](**kwargs)
	else: return import_nanogel_positions_standard(**kwargs)
