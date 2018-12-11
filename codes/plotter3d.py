#!/usr/bin/python

import os
#os.environ['ETS_TOOLKIT']='qt4'
#try: 
from mayavi import mlab
#except: print 'mayavi missing'
import numpy
from numpy import *

def meshplot(data,vecs=None,translate=[0.0,0.0,0.0],show='both',wirecolor=None,
	tri=None,surfcolor=None,lsize=None,maxmin=None,opacity=1.,colormap=None):
	'''Plots a mesh given points in 3-space.'''
	if type(data) == list:
		data = array(data)
	if lsize == None: lsize=2.
	if colormap == None: colormap = 'blue-red'
	if maxmin == None:
		maxmin = [None,None]
	#---Automatically creates the mesh
	if tri == None:
		if shape(data)[1] != 3:
			data = unzipgrid(data,vecs=vecs)
		X,Y,Z = data[:,0]+translate[0],data[:,1]+translate[1],data[:,2]+translate[2]
		pts = mlab.points3d(X, Y, Z, Z)
		mesh = mlab.pipeline.delaunay2d(pts)
		pts.remove()
		if show == 'both':
			surf1 = mlab.pipeline.surface(mesh,representation='wireframe',line_width=lsize,color=(0,0,0),
				vmin=maxmin[0],vmax=maxmin[1],opacity=opacity,colormap='blue-red')
			surf2 = mlab.pipeline.surface(mesh,representation='surface',color=surfcolor,vmin=maxmin[0],
				vmax=maxmin[1],opacity=opacity,colormap=colormap)
		elif show == 'wire':
			if wirecolor == None:
				surf1 = mlab.pipeline.surface(mesh,representation='wireframe',line_width=lsize,color=(0,0,0),
					vmin=maxmin[0],vmax=maxmin[1],opacity=opacity,colormap='blue-red')
			elif wirecolor == 0:
				surf1 = mlab.pipeline.surface(mesh,representation='wireframe',line_width=lsize,color=None,
					vmin=maxmin[0],vmax=maxmin[1],opacity=opacity,colormap='blue-red')
			else:
				surf1 = mlab.pipeline.surface(mesh,representation='wireframe',
					line_width=lsize,color=wirecolor,vmin=maxmin[0],vmax=maxmin[1],opacity=opacity,
					colormap='blue-red')
		elif show == 'surf':
			surf2 = mlab.pipeline.surface(mesh,representation='surface',color=surfcolor,line_width=lsize,
				vmin=maxmin[0],vmax=maxmin[1],opacity=opacity,colormap='blue-red')
		mlab.axes(y_axis_visibility=False,x_axis_visibility=False,z_axis_visibility=False)
		mlab.xlabel("x")
		mlab.ylabel("y")
		mlab.zlabel("z")
	#---Custom mesh
	else:
		if shape(data)[1] != 3:
			data = unzipgrid(data,vecs=vecs)
		X,Y,Z = data[:,0]+translate[0],data[:,1]+translate[1],data[:,2]+translate[2]
		if show == 'both':
			surf1 = mlab.triangular_mesh(X,Y,Z,tri,representation='wireframe',line_width=lsize,color=(0,0,0))
			surf2 = mlab.triangular_mesh(X,Y,Z,tri,representation='surface')
		elif show == 'wire':
			if wirecolor == None:
				surf1 = mlab.triangular_mesh(X,Y,Z,tri,representation='wireframe',
					line_width=lsize,color=(0,0,0))
			elif wirecolor == 0:
				surf1 = mlab.triangular_mesh(X,Y,Z,tri,representation='wireframe',line_width=lsize,color=None)
			else:
				surf1 = mlab.triangular_mesh(X,Y,Z,tri,representation='wireframe',
					line_width=lsize,color=wirecolor)
		elif show == 'surf':
			surf2 = mlab.triangular_mesh(X,Y,Z,tri,representation='surface')
		mlab.axes(y_axis_visibility=False,x_axis_visibility=False,z_axis_visibility=False)
		mlab.xlabel("x")
		mlab.ylabel("y")
		mlab.zlabel("z")

def pbcwire(vecs):

	"""
	Plot a wireframe around the periodic boundary conditions assuming that they begin at the origin.
	"""

	mv = vecs
	mlab.plot3d(array([0.,mv[0],mv[0],0,0.,0.,mv[0],mv[0],0,0]),
		array([0.,0.,mv[1],mv[1],0,0,0,mv[1],mv[1],0]),
		array([0.,0.,0.,0.,0.,mv[2],mv[2],mv[2],mv[2],mv[2]]))
	mlab.plot3d(array([0,0,mv[0],mv[0],mv[0],mv[0]]),
		array([mv[1],mv[1],mv[1],mv[1],0,0]),array([0,mv[2],mv[2],0,0,mv[2]]))

def meshpoints(data,vecs=None,translate=[0.0,0.0,0.0],scale_mode=None,scale_factor=None,color=(0,0,0),
	opacity=1,resolution=8,**kwargs):
	'''Add points to a mesh plot.'''
	if len(data) == 1 or (type(data) == numpy.ndarray and len(data) == 3):
		data = [data]
	if type(data) == list:
		data = array(data)
	X,Y,Z = data[:,0]+translate[0],data[:,1]+translate[1],data[:,2]+translate[2]
	if scale_factor == None:
		scale_factor = 3.
		scale_mode = 'none'
		s = Z
	elif type(scale_factor) == int or type(scale_factor) == float:
		scale_factor = [scale_factor for i in range(len(data))]
		scale_mode = 'scalar'
		s = scale_factor
		scale_factor = 0.0005
	else:
		scale_mode = 'scalar'
		s = scale_factor
		scale_factor = 0.0005
	if shape(data)[1] != 3:
		data = unzipgrid(data,vecs=vecs)
	vec = kwargs.pop('vec',[])
	if len(vec)>0: pbcwire(vec)
	pts = mlab.points3d(X,Y,Z,s,scale_mode=scale_mode,mode='sphere',
		scale_factor=0.5,color=color,opacity=opacity,resolution=resolution)
	return pts
	
def unzipgrid(surf,vecs=None,grid=None,rounder_vecs=[1.0,1.0],reverse=0):
	'''Turns a 2D array into a set of points in 3-space for meshplot and meshpoints.'''
	#---Nb: this function is also a child of the MembraneSet class.
	#---Nb: duplicated here for plotting purposes only.
	if type(surf) != ndarray:
		surf = array(surf)
	grid = [shape(surf)[i] for i in range(2)]
	if reverse != 0: grid = grid[::-1];
	if vecs != None:
		rounder_vecs = [vecs[i]/(grid[i]-1) for i in range(2)]
	replotin = surf
	surfreplot = []
	for i in range(grid[0]):
			for j in range(grid[1]):
				surfreplot.append([i*rounder_vecs[0],j*rounder_vecs[1],replotin[i,j]])
	surfreplot = array(surfreplot)
	return surfreplot

