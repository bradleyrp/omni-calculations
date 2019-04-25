#!/usr/bin/env python

"""
Refactored curvture-undulation coupling code.

USAGE:
  make set meta_filter=dextran_20190415.yaml
  make go script=calcs/dextran_undulations_review.py

You can develop new code in the __main__ section.
Reexecute your new code with the `main` command.
"""

from omni import WorkSpace
import calcs
from calcs.codes.curvature_coupling.curvature_coupling_plots \
import individual_reviews_plotter
from omni.base.tools import gopher
from ortho import delve
import copy

# keys for the specs file
plotname = 'curvature_undulation_coupling'
calcname = 'curvature_undulation_coupling_dex'

# settings for various plots
plotspec = {
	'coupling_review':{
		'viewnames':['average_height','average_height_pbc',
			'neighborhood_static',
			'neighborhood_dynamic','average_field','example_field',
			'example_field_pbc',
			'spectrum','spectrum_zoom']},
	'coupling_review.simple':{
		'viewnames':['average_height','example_field'],'figsize':(6,6)},
	'coupling_review.center_debug':{
		'viewnames':['neighborhood_static','neighborhood_dynamic',
			'average_height','average_height_pbc','average_field',
			'example_field','example_field_pbc',
			'average_height_center','average_field_center',
			'average_field_center_smear','curvature_field_center'],
			'figsize':(16,16)},
	'coupling_review.simple_centered':{
		'viewnames':['average_height_center','curvature_field_center',
		'coupling_review.simple'],
		'figsize':(8,8),'horizontal':True,'wspace':0.7},}


@loader
def load():
	"""
	Load the data once.
	Develop below in __main__ and reexecute with the `main` command.
	"""
	work = WorkSpace(analysis=True)
	plotspecs = work.metadata.plots[plotname].get('specs',{})
	calcname = plotspecs.get('calcname',plotname)
	# collect data
	data = work.plotload(plotname=plotname)
	sns = sns_ordered = work.sns()
	# separate plotload for the entire undulation coupling loop
	data_ucc = this = work.plotload(calcname)
	# fetch all hypotheses for the calculation
	upstreams,upstreams_stubs = work.calcs.unroll_loops(
		work.metadata.calculations[calcname],return_stubs=True)
	# get tags for the hypotheses
	tags = [delve(i,'specs','design') for i in upstreams_stubs]

def plot_results(tag,plotspec):
	"""
	Plot a few versions of the main figure.
	"""
	global data,upstreams,tags
	plotspec = copy.deepcopy(plotspec)
	# turn some off when developing
	for i in [
		'coupling_review.center_debug','coupling_review.simple_centered'
		]: plotspec.pop(i)
	# lookups require tags in the same order as upstreams
	data_ucc.set('curvature_undulation_coupling_dex',
		select=upstreams[tags.index(tag)]['specs'])
	# merging the data into an object expected by the plotter
	#! this merge is required for the outdated plotter code
	datas = {tag:data_ucc.this}
	postdat = {}
	for sn in sns:
		postdat[sn] = {}
		#! reformulating postdat for the plot scripts here too
		for key in [
			#! order matters here because vecs might get overwritten!
			'import_readymade_meso_v1_nanogel',
			'import_readymade_meso_v1_membrane',]:
			data.set(key)
			for i,j in data.this[sn].items():
				postdat[sn][i] = j
	#! patching in the extra stuff done by curvature_coupling_lodaer_dex
	for sn in sns:
		data.set('import_readymade_meso_v1_nanogel')
		postdat[sn]['points_protein'] = data.this[sn]['points_all']
	#! vecs are now different too
	for sn in sns: postdat[sn]['vecs'] = postdat[sn]['vecs'].mean(axis=0)
	#! yet another item from curvature_coupling_loader_dex
	for sn in sns: 
		postdat[sn]['points_protein_mean'] = \
			postdat[sn]['points_protein'].mean(axis=1).mean(axis=0)
	seep = dict(data=data,datas=datas,work=work,postdat=postdat,
		undulations_name=plotspecs['undulations_name'],
		protein_abstractor_name=plotspecs['protein_abstractor_name'])
	for out_fn,details in plotspec.items(): 
		report_this = individual_reviews_plotter(
			out_fn=out_fn,seep=seep,**details)
		if out_fn=='coupling_review': report = report_this
	# save the report
	return report

def plot_recap(report):
	"""
	Bar plots of curvature and error across hypotheses.
	"""
	import re
	#! many hardcoded features of this plot
	for key,val in report.items(): 
		report[key]['style'] = re.match(
			'^(SemiRigid|Rigid|CL160ENS)',key[0]).group(1)
		if re.match(r'^SemiRigidE\d+_50',key[0]):
			report[key]['radius'] = 50
		elif re.match(r'^RigidNCE\d+',key[0]):
			report[key]['radius'] = 50
		elif re.match(r'^SemiRigid',key[0]):
			report[key]['radius'] = 25
		else: report[key]['radius'] = -1
	groups = [
		('SemiRigidE1','SemiRigidE2','SemiRigidE3'),
		('RigidNCE1','RigidNCE2','RigidNCE2',),
		('SemiRigidE1_50','SemiRigidE2_50','SemiRigidE3_50'),
		('CL160ENS-1','CL160ENS-2','CL160ENS-3',),]
	maxes = {'max_curvature':0.15,'error':0.014}
	titles = {'v3':'neighborhood\n(50,40x20)',
		'v5':'neighborhood\n(100,40x20)','v4':'pixel\n(80x40)'}
	import pandas as pd
	import matplotlib as mpl
	import matplotlib.pyplot as plt
	from omni.plotter.panels import square_tiles
	from omni.base.store import picturesave
	counter = 0
	#! axes,fig = square_tiles(4,(10,10),)
	tags = ['v3','v5','v4']
	fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(12,12))
	color_list = ['krb','krb','krb']
	# seaborn for better colors
	import seaborn as sb
	pal = sb.color_palette()
	color_list = [pal[:len(groups)] for i in range(len(titles))]
	for inum,item in enumerate(['max_curvature','error']):
		for tnum,tag in enumerate(tags):
			ax = axes[inum][tnum]
			result = {}
			for key in ['error','style','radius','max_curvature']:
				result[key] = [[report[(i,tag)][key] 
					for i in g] for g in groups]
			df = pd.DataFrame(result[item])
			ax.set_title('%s %s'%(item,titles[tag]))
			#! df.groupby('styles').errors.value_counts().plot.bar()
			#! counter intuitiuve color ordering
			plot = df.plot.bar(color=color_list,ax=ax)
			ax.get_legend().remove()
			ax.set_xticklabels(['SemiRigid','Rigid',
				'Semi-\nRigid (50)','Dextran'],rotation=0)
			ax.set_ylim((0,maxes[item]))
	picturesave('fig.summary',work.plotdir,backup=False,version=True,meta={})

def make_images_and_report():
	"""Make all images and summary bar plot."""
	global plotspec
	report = {}
	for tag in tags:
		# save features of each hypothesis set while images are rendered
		result = plot_results(tag=tag,plotspec=plotspec)
		for sn,tag in result:
			report.update(result)
	# plot a summary bar plot from the features we saved during plots
	plot_recap(report=report)
	return report

if __name__=='__main__':

	# to remake images delete the report otherwise they are only made once
	if 'report' not in globals(): report = make_images_and_report()

	kwargs = {}

	kB = 1.38E-23
	temp = 300.00
	mass = 5.5E-22
	Plank = 6.3E-34
	Area_mem = 0.25E-12
	my_const = 3.1415


	if len(args)!=2: raise Exception('arguments list must be two items long: qs and energies')
	high_cutoff = kwargs.pop('high_cutoff',1.0)
	kappa,sigma = kwargs.pop('kappa'),kwargs.pop('sigma')
	if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
	#---the wavevectors and corresponding energies are the arguments
	qs,energy = args
	freq, oper = np.zeros(qs.shape,'d'), np.zeros(qs.shape,'d')
	nqs = len(qs)
	for j in range(0,nqs):
		if qs[j] > 0.0:
			freq[j] = np.sqrt(Area_mem*(kappa*qs[j]**4+sigma*qs[j]**2))/np.pi
			oper[j] = freq[j]*Plank/(kB*temp)

	Entropy_term = 0.0
	for j in range (0,nqs):
	if qs[j] > 0.0:
	Entropy_term += (oper[j]/(np.exp(oper[j])-1.0))-np.log(1.0-np.exp(-oper[j]))

	Entropy = Entropy_term*kB*6.022140857E23    # un j/K/mol
