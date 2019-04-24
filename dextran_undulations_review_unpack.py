#!/usr/bin/env python

"""
Unpack the curvature-undulation coupling data.

Usage:

cd /home/localshare/factory/calc/dextran
source ../../../env/bin/activate py3
make set meta_filter dextran_20190407_v5_neighbor_100.yaml
make go script=calcs/dextran_undulations_review_unpack.py
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from omni import WorkSpace
import calcs
from calcs.codes.curvature_coupling.curvature_coupling_plots \
	import individual_reviews_plotter
from omni.base.tools import gopher
from omni.base.store import picturesave
import numpy as np
import pandas as pd

# plot details in the spcs file
plotname = 'curvature_undulation_coupling'

# desired calculation details that reflect the specs file
spec_this = ('curvature_undulation_coupling_dex',
 {'design': {'curvature_positions': {'distance_cutoff': 100.0,
                                     'method': 'neighborhood',
                                     'spacer': 40.0},
             'curvature_sum': 'mean',
             'extents': {'extent': 20.0, 'method': 'fixed_isotropic'},
             'style': 'wilderness'},
  'fitting': {'high_cutoff': 0.2, 'initial_kappa': 10.0},
  'upstream': {'import_readymade_meso_v1_membrane': {},
               'import_readymade_meso_v1_nanogel': {}}})

@loader
def load():
	work = WorkSpace(analysis=True)
	plotspecs = work.metadata.plots[plotname].get('specs',{})
	calcname = plotspecs.get('calcname',plotname)
	data = work.plotload(plotname=plotname)
	sns = sns_ordered = work.sns()


if __name__=='__main__':
	data.set('curvature_undulation_coupling_dex')
	# since the analysis depends on a setting in the specs we confirm
	if not data.cursor==spec_this:
		raise Exception('wrong specs file')
	# the kappa and gamma (or sigma) are the first two values of the fit ('x')
	params = dict([(sn,dict(list(zip(('kappa','gamma'),
		data.this[sn]['x'][:2])))) for sn in sns])
	fits = pd.DataFrame(params).T
	print('status Fitted parameters are stored in `fits`:')
	print(fits)
	print(('status Instantaneous curvature fields are '
		'stored in the dict `c0` while average fields are in `c0_avg` '
		'which can be saved from the terminal'))
	c0_avg,c0 = [dict([(sn,
		data.this[sn][key])
		for sn in sns])
		for key in ['cf','cf_first']]

data.this['SemiRigidE1'].keys()
dict_keys(['cf', 'cf_first', 'drop_gaussians_points', 'jac', 'qs', 'ratios', 'x', 'bundle', 'spec'])

>>> data.this['SemiRigidE1']['qs']
array([0.01208305, 0.0241661 , 0.03624915, ..., 0.03820995, 0.02701852,
       0.01708801])
>>> data.this['SemiRigidE1']['x']
array([ 1.68908262e+01,  3.20481392e-04,  6.32755638e+00, -4.65822755e-01,
        9.29523225e-01,  1.31834903e+00,  8.15226137e-01,  1.21394647e-01,
        7.65267589e-01, -4.09788641e-01, -2.36197942e-01, -2.07752748e-01,
        1.07601565e+00,  6.25013665e-01, -1.94292366e-01, -2.03905978e+00,
       -3.77792485e-01,  8.47226593e-01,  9.18135134e-02,  7.88798517e-01,
       -2.93867055e-01,  3.52463297e-01,  4.95388197e-01,  7.47588380e-02,
       -1.30486786e-01,  3.81546891e-01,  2.93538171e-01, -1.04181256e-01])