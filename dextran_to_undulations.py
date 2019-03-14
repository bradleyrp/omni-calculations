#!/usr/bin/env python

"""
Mistakenly started this file when I was on the wrong commits for some old codes

"""

import glob
import re

import numpy as np

base = "/Users/rpb/worker/simulations/samrad-curvature-2019.03.00/ResultsForDex/CL160ENS-1"
positions = 'NG_Position-0_Frame-0.dat'

with open(os.path.join(base,positions)) as fp: pos_raw = fp.read()

regex_entry = r'(\d+)\n(.*?)\n(?:\d+\n|\Z)'
entries = re.findall(r'(\d+)\n(.*?)\n(?:\d+\n|\Z)',pos_raw,flags=re.M+re.DOTALL)
regex_points = r'(\d+)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)'
proc = [(int(i),np.array(
	re.findall(regex_points,j,flags=re.M+re.DOTALL)).astype(float)) 
	for i,j in entries]

fns = glob.glob(os.path.join(base,'memb_xyz-0-0-*.xyz'))
fn = fns[0]

with open(fn) as fp: raw = fp.read()

regex_snap = r'^\s*(\d+)\s+(\d+)\s+(\d+(?:\.\d+)?)\s*\n(.+)'
regex_col4 = r'\s*'+r'(\d+(?:\.\d+)?)(?:\s+|\Z)'*4
this = re.match(regex_snap,raw,flags=re.M+re.DOTALL).groups()
detail,data = [np.array(this[:3]).astype(float)]+[np.array(re.findall(regex_col4,this[3])).astype(float)]

import matplotlib as mpl
from matplotlib import pylab as plt

# plt.scatter(*data[:,:2].T)
# plt.show()

#! . ../../env/bin/activate py3
if 1:
	import mayavi
	from mayavi import mlab
	x,y,z = data[:,:3].T
	mlab.points3d(x,y,z,s=0.5)

