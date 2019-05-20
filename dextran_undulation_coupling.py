#!/usr/bin/env python
# vim: set ts=4 sts=4 sw=4 noet

"""
Refactored undulation-curvature coupling (ucc) code circa 2019.04.15.
Use this code to analyze the undulation-curvature coupling results.

USAGE:
  Run the following once:
  	make set meta_filter=dextran_20190424.yaml
  Run this command to enter the plot environment:
  	make go script=calcs/dextran_undulation_coupling.py
  When developing code at the end (__main__), use the `main` command to execute.
"""

from omni import WorkSpace
import calcs
from calcs.codes.curvature_coupling.curvature_coupling_plots \
import individual_reviews_plotter
from calcs.codes.undulate import calculate_undulations
from omni.base.tools import gopher
from ortho import delve,treeview
import copy
import pandas as pd
import numpy as np
from ortho import status


# no plots by default
control.routine = []

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
	# prepare a readable set of hypotheses
	contents,contents_full = {},{}
	# the tags and specs are in sequence from unroll_loops
	for tag,spec,stub in zip(tags,upstreams,upstreams_stubs):
		contents[tag] = dict([(k,spec['specs'][k]) 
			for k in ['design','fitting']])
		contents_full[tag] = spec

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

@autoplot
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

@function
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

def package_ucc(dat,sns):
	"""
	Organize the results from one hypothesis.
	"""
	# prepare summary
	# the kappa and gamma (or sigma) are the first two values of the fit ('x')
	params = dict([(sn,dict(list(zip(('kappa','gamma'),
		dat[sn]['x'][:2])))) for sn in sns])
	fits = pd.DataFrame(params).T
	print('status Fitted parameters are stored in `fits`:')
	print(fits)
	print(('status Instantaneous curvature fields are '
		'stored in the dict `c0` while average fields are in `c0_avg` '
		'which can be saved from the terminal'))
	c0_avg,c0 = [dict([(sn,
		dat[sn][key])
		for sn in sns])
		for key in ['cf','cf_first']]
	qs =  dict([(sn,dat[sn]['qs']) for sn in sns])
	energy =  dict([(sn,dat[sn]['ratios']) for sn in sns])
	# all results are keyed by simulation name
	return dict(c0_avg=c0_avg,c0=c0,fits=fits,qs=qs,energy=energy)

def calculate_entropy(*args,**kwargs):
	"""
	Compute the entropy with wavevectors and corresponding energies.
	Requires inputs from curvature-undulation coupling.
	"""
	# previous method, predates the one coded in __main__ below
	kB = 1.38E-23
	temp = 300.00
	mass = 5.5E-22
	Plank = 6.3E-34
	Area_mem = 0.25E-12
	my_const = 3.1415

	if len(args)!=2: 
		raise Exception('arguments list must be two items long: qs,energies')
	high_cutoff = kwargs.pop('high_cutoff',1.0)
	kappa,sigma = kwargs.pop('kappa'),kwargs.pop('sigma')
	if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
	# the wavevectors and corresponding energies are the arguments
	qs,energy = args
	freq,oper = np.zeros(qs.shape,'d'), np.zeros(qs.shape,'d')
	nqs = len(qs)
	for j in range(0,nqs):
		if qs[j] > 0.0:
			freq[j] = np.sqrt(Area_mem*(kappa*qs[j]**4+sigma*qs[j]**2))/np.pi
			oper[j] = freq[j]*Plank/(kB*temp)
	Entropy_term = 0.0
	for j in range (0,nqs):
		if qs[j] > 0.0:
			Entropy_term += (oper[j]/(
				np.exp(oper[j])-1.0))-np.log(1.0-np.exp(-oper[j]))
	Entropy = Entropy_term*kB*6.022140857E23 # j/K/mol
	return Entropy

if __name__=='__main__':

	"""
	USAGE
	You can develop code here. Run it with `main` after saving the file.
	Run `make_images_and_report` to remake images and summary bar plots.
	"""

	# ENTROPY METHOD 1 (set True to try this)
	if 1:

		# review the hypotheses
		treeview(dict(hypotheses=contents))

		# select a hypothesis by tag from the hypotheses dict
		tag_this = 'v5'

		# get the data from the undulation-curvature coupling object
		data_ucc.set(select=contents_full[tag_this]['specs'])
		dat = data_ucc.this
		result = package_ucc(dat,sns)

		# select a simulation
		#! note that this uses the original entropy function above
		#!   this was replaced with the entropy method below on 2019.05.12
		#!   and that function was added to dextran_undulation.py on 2019.04.26
		for sn in sns:
			#! result is currently inf
			entropy = calculate_entropy(
				result['qs'][sn],
				result['energy'][sn],
				kappa=result['fits']['kappa'][sn],
				sigma=result['fits']['gamma'][sn],)

		# entropy calculation starts here
		# this mimics the version in dextran_undulation.py from 2019.04.26
		# that script gets uspec from a function called calculate_undulations
		# we populate uspec here, from the curvature coupling data
		#! note result is populated above while results are populated below, to
		#!   match dextran_undulation.py
		results = {}
		for snum,sn in enumerate(sns):
			status('calculating entropy for %s'%sn,i=snum,looplen=len(sns))

			"""
			Notes on implementation:
			At this point we are getting the variable "result" from data_ucc.
			The data_ucc contains the undulation coupling results which contain the fitted kappa and sigma. However, that method is fitting the energy (and not the hqhq). The following construction of uspec is designed to mimic the result of a calculate_undulation function used in the dextran_undulation.py code. It puts the data in the same exact structure used by the two entropy calculations below. After we build it, we are going to replace energy_raw, which is an energy value optimized to match kBT, with the hqs, because the entropy calculation needs the hqs and not the energy.
			"""
			uspec = dict(q_raw=result['qs'][sn],
				energy_raw=result['energy'][sn],
				kappa=result['fits']['kappa'][sn],
				sigma=result['fits']['gamma'][sn])

			### INTERVENTION
			"""
			this code was taken from dextran_undulation.py and provides the hqs directly
			BE CAREFUL THAT YOU ARE ONLY COMPARING THE SAME midplane_method!
			"""
			# wavevector limits
			lims = (0.,1.0)
			# select a midplane method
			midplane_method = [
				'flat','average','average_normal',
				][0]

			if midplane_method=='average_normal':
				data.set('import_readymade_meso_v1_membrane')
				dat = data.this[sn]
				vecs = dat['vecs']
				data.set('undulations_average_normal')
				dat = data.this[sn]
				surf = mesh = dat['average_normal_heights']
				custom_heights = surf
			else: 
				data.set('import_readymade_meso_v1_membrane')
				dat = data.this[sn]
				mesh = dat['mesh']
				custom_heights = None
				vecs = dat['vecs']
				# assume bilayer, even if one mesh, then take the average
				surf = np.mean(mesh,axis=0)

			if 1: surf -= np.tile(surf.reshape(len(surf),-1).mean(axis=1),
				(surf.shape[1],surf.shape[2],1)).transpose((2,0,1))

			# kernel of this plot/calculation: calculate the spectra here
			uspec_no_coupling = calculate_undulations(surf,vecs,chop_last=True,
				custom_heights=custom_heights,
				perfect=True,lims=lims,raw=False,midplane_method=midplane_method,
				fit_style='band,perfect,curvefit',fit_tension=True)
			# note that this value is called "energy_raw" but it is really hqs
			#!!! correct this nomenclature error
			wrong_method = False
			if not wrong_method: uspec['energy_raw'] = uspec_no_coupling['energy_raw']
			### END INTERVENTION

			kB = 1.38E-23
			temp = 300.00
			mass = 5.5E-22
			area_mem = 0.25E-12 
			Plank = 6.62607015*10**-34

			# note that uspec includes the zeroth mode which we omit
			#   this is necessary to match the method from 
			#   dextran_undulation_coupling.py
			kappa = uspec['kappa']
			sigma = uspec['sigma']
			qs = uspec['q_raw']
			hqs = hq2 = uspec['energy_raw']
			"""
			# add back the zero mode:
			The above is deprecated. You added the zero mode because the energy value coming from the coupling code had dropped the zero mode but we need the hqs anyway, not the energy, hence we skip this step.
			"""
			if wrong_method:
				hqs = hq2 = np.concatenate(([0.],hqs))
				qs = np.concatenate(([0.,],qs))
				square_dim = np.sqrt(qs.shape)[0]
				if square_dim.astype(int)!=square_dim:
					raise Exception('dimension error! data might not be square '
						'because we have %s items'%qs.shape)
				else: square_dim = square_dim.astype(int)
				nx,ny = (square_dim,square_dim)
				hq2_square = np.reshape(hq2,(nx,ny))

			#! an older code
			if 0:
				w_q1 = np.sqrt(area_mem*(kappa*qs**4+sigma*qs**2)/mass)
				w_q2 = np.sqrt(kB*temp/(hq2**2*mass)) 

			# CALCULATE ENTROPY
			# code from 2019.04.26

			freq1, oper1 = np.zeros(qs.shape,'d'), np.zeros(qs.shape,'d')
			freq2, oper2 = np.zeros(qs.shape,'d'), np.zeros(qs.shape,'d')
			nqs = len(qs)
			for j in range(0,nqs):
				if qs[j] > 0.0:
					freq1[j] = np.sqrt(area_mem*(
							kappa*kB*temp*qs[j]**4*10**36+
							sigma*kB*temp*qs[j]**2*10**36)/mass)
					oper1[j] = freq1[j]*Plank/(kB*temp)/(2*np.pi)
					freq2[j] = np.sqrt(kB*temp/(hqs[j]*10**-18*mass))
					oper2[j] = freq2[j]*Plank/(kB*temp)/(2*np.pi)

			Entropy_term_1 = 0.0
			Entropy_term_2 = 0.0    

			for j in range (0,nqs):
				if qs[j] > 0.0:
					Entropy_term_1 += (oper1[j]/(
						np.exp(oper1[j])-1.0))-np.log(1.0-np.exp(-oper1[j]))
					Entropy_term_2 += (oper2[j]/(
						np.exp(oper2[j])-1.0))-np.log(1.0-np.exp(-oper2[j]))

			Entropy_1 = Entropy_term_1*kB*6.022140857E23    # un j/K/mol
			Entropy_2 = Entropy_term_2*kB*6.022140857E23    # un j/K/mol        

			# save the results here
			results[sn] = dict(
				kappa=uspec['kappa'],
				sigma=uspec['sigma'],
				Entropy_1=Entropy_1,
				Entropy_2=Entropy_2,
			)

		### NOTE! Ryan is checking to see why the zeroth mode is not zero in dextran_undulation.py!
