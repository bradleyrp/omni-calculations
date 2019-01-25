#!/usr/bin/env python

import sys,time
import numpy as np
from omni import basic_compute_loop 
from calcs.codes.mesh import makemesh

def lipid_mesh(**kwargs):
	"""
	Compute monolayer mesh objects.
	"""
	# parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	dat = kwargs['upstream']['lipid_abstractor']
	resnames = dat['resnames']
	monolayer_indices = dat['monolayer_indices']
	nframes = dat['nframes']
	debug = kwargs.pop('debug',False)
	kwargs_out = dict(curvilinear=
		calc.get('specs',{}).get('curvilinear',False))
	imono = [np.where(monolayer_indices==mn) for mn in range(2)]
	indices = [(mn,fr) for fr in range(nframes) for mn in range(2)]
	looper = [dict(kwargs_out,**{
		'pts':dat['points'][fr][imono[mn]],'vec':dat['vecs'][fr]})
		for (mn,fr) in indices]
	incoming = basic_compute_loop(compute_function=makemesh,looper=looper)
	# pack
	attrs,result = {},{}
	result['nframes'] = np.array(nframes)
	result['vecs'] = dat['vecs']
	result['resnames'] = resnames
	result['monolayer_indices'] = monolayer_indices
	# pack mesh objects
	# keys should be uniform and include: vertnorms, simplices, nmol,
	#   facenorms, gauss, points, vec, ghost_ids, mean, principals, areas
	keylist = incoming[0].keys()
	for key in keylist:
		# same loop nesting as looper above
		for ii,(mn,fr) in enumerate(indices):
			result['%d.%d.%s'%(mn,fr,key)] = incoming[ii][key]
	return result,attrs	
