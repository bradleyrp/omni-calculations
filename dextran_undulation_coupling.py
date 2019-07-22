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
import seaborn as sb
import numpy as np
from ortho import status
import matplotlib as mpl
import re
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from omni.plotter.panels import square_tiles,panelplot
from omni.base.store import picturesave
import scipy
import json

#! https://stackoverflow.com/questions/31908982
#! requirements for the fancy legend
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle,Circle
from matplotlib.legend_handler import HandlerBase
class HandlerColormap(HandlerBase):
	def __init__(self, custom_colors, **kw):
		HandlerBase.__init__(self, **kw)
		self.custom_colors = custom_colors
		self.num_stripes = len(custom_colors)
	def create_artists(self, legend, orig_handle, 
		xdescent, ydescent, width, height, fontsize, trans):
		stripes = []
		for i in range(len(self.custom_colors)):
			s = Rectangle([xdescent + i * width / self.num_stripes, ydescent], 
				width / self.num_stripes, 
				height, 
				#! fc=self.cmap((2 * i + 1) / (2 * self.num_stripes)), 
				fc = self.custom_colors[i], 
				transform=trans)
			stripes.append(s)
		return stripes

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
		'figsize':(8,8),'horizontal':True,'wspace':0.7},
	'coupling_review.three':{
		'viewnames':['average_height','average_field','spectrum'],
		'panelplot_layout':{'out':{'grid':[1,1]},'ins':[
			{'grid':[3,1],'hratios':[1,1,0.5]}]},'tickoff':False,
		#! 'panelplot_layout':{'out':{'grid':[2,1],
		#!   'hratios':[1,0.5]},'ins':[{'grid':[1,2]},{'grid':[1,1]}]},
		'figsize':(3,10),'horizontal':True,'wspace':0.7},}

@loader
def load():
	"""
	Load the data once.
	Develop below in __main__ and reexecute with the `main` command.
	"""
	work = WorkSpace(analysis=True)
	plotspecs = work.metadata.plots[plotname].get('specs',{})
	calcname = plotspecs.get('calcname',plotname)
	#! this is necessary for dextran_2019711_combine_plot.yaml
	colcomb = work.metadata.director.get('collections_combine',False) 
	plargs = {}
	if colcomb: plargs = dict(collections=['all150','all100'])
	# collect data
	data = work.plotload(plotname=plotname,**plargs)
	sns = sns_ordered = work.sns()
	# custom method for merging for a specific meta filter
	is_custom = ('dextran_20190711_combine_plot.yaml' in 
		[os.path.basename(i) for i in work.metadata.specs_files])
	if is_custom:
		# custom override
		request_calc = {
			'curvature_undulation_coupling_dex_1':
				work.metadata.calculations['curvature_undulation_coupling_dex_1'].get('specs',{}),
			'curvature_undulation_coupling_dex_2':
				work.metadata.calculations['curvature_undulation_coupling_dex_2'].get('specs',{})}
		calcnames_this = [
			'curvature_undulation_coupling_dex_1',
			'curvature_undulation_coupling_dex_2']
		plargs['request_calc_override'] = request_calc
	else: calcnames_this = [calcname]
	# separate plotload for the entire undulation coupling loop
	data_ucc = this = work.plotload(calcname,**plargs)
	contents,contents_full = {},{}
	# loop over calculations contributing to this sweep
	for calcname in calcnames_this:
		upstreams,upstreams_stubs = work.calcs.unroll_loops(
			work.metadata.calculations[calcname],return_stubs=True)
		# get tags for the hypotheses
		tags = [delve(i,'specs','design') for i in upstreams_stubs]
		# the tags and specs are in sequence from unroll_loops
		for tag,spec,stub in zip(tags,upstreams,upstreams_stubs):
			if tag in contents: 
				raise Exception(('contents already has this tag (%s) '
					'hence the merge failed!')%tag)
			contents[tag] = dict([(k,spec['specs'][k]) 
				for k in ['design','fitting']])
			contents_full[tag] = spec
		# more sensible name
		design_sweep = contents_full
	if is_custom:
		#! this merges simulations from the two calculations (_1, _2) 
		#!   so they are available in one "this". note that data_ucc.set 
		#!   should not be needed for data_ucc. be careful with this!
		data_ucc.this = dict([(sn,data_ucc.data[ii][sn]['data']) 
			for ii,i in enumerate(data_ucc.names) for sn in data_ucc.data[ii]]) 
		#! !!! sns = data_ucc.this.keys()
		sns = ['RigidNCE1', 'RigidNCE2', 'RigidNCE3', 'SemiRigidE1_50', 
			'SemiRigidE2_50', 'SemiRigidE3_50', 'CL160ENS-1', 'CL160ENS-2', 'CL160ENS-3']

def color_vary(colors,n,p=0.2):
	"""
	Given a list of colors, make some variations around those colors.
	Originally designed to trace a path through color space, however it makes
	more sense to just vary around a base color by darkness levels.
	"""
	colors_out = []
	for cnum,color in enumerate(colors):
		left,center,right = [np.array(mpl.colors.to_rgb(i)) for i in ['w',color,'k']]
		gaps = [np.linalg.norm(i-center) for i in [left,right]]
		# discretize the progress along two lines
		subdiv = np.linspace(1-p*1.,1+p*1.,n)
		colors_sub = []
		for ii,i in enumerate(subdiv):
			if i<=1: x = left+i*(center-left)
			else: x = center+(i-1.0)*(right-center)
			#! scaling by gaps was wrong here! write it in your log!
			if any([i<0 or i>1 for i in x]):
				raise Exception('error')
			colors_sub.append(x)
		colors_out.append(colors_sub)
	# reorder for pandas
	color_list = [colors_out[j][i] for i in range(n) for j in range(len(colors))]
	return color_list

def plot_results(tag,plotspec):
	"""
	Plot a few versions of the main figure.
	"""
	report = {}
	global data,upstreams,tags
	plotspec = copy.deepcopy(plotspec)
	# turn some off when developing
	for i in [
		'coupling_review.center_debug','coupling_review.simple_centered',
		'coupling_review.simple',
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
		yticks = np.arange(-0.08,0.08+0.02,0.02).round(10)
		report_this = individual_reviews_plotter(
			#! bootstrap the following maximum from the report if the data change
			#! true maximum curvature for v5 is: max([report[(sn,'v5')]['max_curvature'] for sn in data.this]
			#!   which is 0.08865340495632093 or 0.09
			#! custom color bar limits
			custom_cbar=dict(
				title=r'$\mathrm{\langle C_0(x,y) \rangle}$'+'\n'+r'$\mathrm{({nm}^{-1}\times{10}^2})$',
				yticks=yticks,ytick_labels=['%d'%(i*10**2) for i in yticks],
				shift_up=1.2),
			max_curvature_overall=0.09,
			use_instantaneous_grid_on_average=True,
			out_fn=out_fn,seep=seep,**details)
		if out_fn=='coupling_review': report = report_this
	# save the report
	return report

@autoplot
def plot_recap(report,subset=False,bw=False,which='error_curvature'):
	"""
	Bar plots of curvature and error across hypotheses.
	"""
	if which=='error_curvature':
		items = ['max_curvature','error']
		figsize = 8
		fn_base = 'summary'
		legend_args = dict(loc='upper left',
			bbox_to_anchor=(1.0,0.0,1.,1.))
	elif which=='strength':
		items = ['strength']
		figsize = 6
		fn_base = 'summary_strength'
		legend_args = dict(loc='upper right')
	else: raise Exception('inclear plot target: %s'%which)
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
	if not subset:
		groups = [
			('SemiRigidE1','SemiRigidE2','SemiRigidE3'),
			('RigidNCE1','RigidNCE2','RigidNCE3',),
			('SemiRigidE1_50','SemiRigidE2_50','SemiRigidE3_50'),
			('CL160ENS-1','CL160ENS-2','CL160ENS-3',),]
		groups_names = ['Semi-Rigid','Rigid','Semi-Rigid 50','Nanogel']
		tick_label_list = ['SemiRigid','Rigid','Semi-\nRigid (50)','Dextran']
		groups_names = ['Rigid-Tethered','Rigid','Rigid_Tethered 50','Flexible']
		tick_label_list = ['Rigid-Tethered','Rigid','Rigid-\nTethered (50)','Flexible']
	else:
		groups = [
			('RigidNCE1','RigidNCE2','RigidNCE3',),
			('SemiRigidE1_50','SemiRigidE2_50','SemiRigidE3_50'),
			('CL160ENS-1','CL160ENS-2','CL160ENS-3',),]
		tick_label_list = ['Rigid','Semi-\nRigid','Dextran']
		groups_names = ['Rigid','Semi-Rigid','Nanogel']
		tick_label_list = ['Rigid','Rigid-\nTethered','Flexible']
		groups_names = ['Rigid','Rigid-Tethered','Flexible']
	#! before the latest round of dextran calculations, now being plotted on
	#!   2019.06.03, the maximum error was beneath 0.014 so the recent
	#!   calculations from roughly last week might have high error for some 
	#!   unknown reason that needs investigated. 
	#! once I recovered the old data I set if rom 0.09 back to 0.014
	#! max strength for the subset is 
	#!   np.array([[report[(i,'v5')]['strength'] for i in g] for g in groups]).max()
	maxes = {'max_curvature':0.15,'error':0.014,'strength':0.0009}
	titles = {'v3':'neighborhood\n(50,40x20)',
		'v5':'neighborhood\n(100,40x20)','v4':'pixel\n(80x40)'}
	counter = 0
	if not subset:
		tags = ['v3','v5','v4']
		fig,axes = plt.subplots(nrows=len(items),ncols=len(tags),figsize=(12,12))
		extra_tag = ''
	else:
		#! limiting tags here now that we have settled on one
		tags = ['v5']
		axes,fig = square_tiles(len(items),figsize=figsize,
			favor_rows=True,wspace=0.4,hspace=0.1)
		axes = [[axes[i]] for i in range(len(items))]
		extra_tag = '.spec_v5'

	# seaborn for better colors
	import seaborn as sb
	pal = sb.color_palette()
	if subset: colors = ['#e41a1c','#377eb8','#4daf4a']
	else: colors = 'rgbm'
	color_list = color_vary(colors=colors,n=3,p=0.2)
	for inum,item in enumerate(items):
		for tnum,tag in enumerate(tags):
			ax = axes[inum][tnum]
			result = {}
			for key in ['error','style','radius','max_curvature','strength']:
				result[key] = [[report[(i,tag)][key] 
					for i in g] for g in groups]
				#! notes on "replicate-grouped bars": flatten it
				#! result[key] = [report[(g,tag)][key] for group in groups for g in group]
			df = pd.DataFrame(result[item])
			#! alternate formulation: prepare the right tags
			sn_to_group = [(g,report[(g,tag)]['style']) for group in groups for g in group]
			index = [i[1] for i in sn_to_group]
			#! df = pd.DataFrame(result[item],index=index)
			if 0: df = pd.DataFrame(zip(index,result[item]),columns=['index','error'])
			#! ax.set_title('%s %s'%(item,titles[tag]))
			if item=='error':
				ax.set_ylabel('error')
			elif item=='max_curvature':
				ax.set_ylabel(r'$\mathrm{C_{0,max}}\,({nm}^{-1})$')
			elif item=='strength':
				ax.set_ylabel(r'$\mathrm{\frac{\int{}{} C_0 dA}{\int{}{} dA}\,({nm}^{-2})}$')
			#! notes on "replicate-grouped bars": try using groupby
			#! df.groupby('styles').errors.value_counts().plot.bar()
			#! counter intuitiuve color ordering
			#! df = df.unstack()
			plot = df.plot.bar(ax=ax,edgecolor='black',stacked=False)
			#! notes on "replicate-grouped bars": we have to set the bars manually
			#!   because pandas expects each member of a group to be the same color
			bars = [this for this in plot.get_children()
				if this.__class__.__name__=='Rectangle' and this._width<1]
			if len(color_list)<len(bars):
				raise Exception('we have %d bars and only %d colors'%(len(bars),len(color_list)))
			for bnum,bar in enumerate(bars):
				bar.set_facecolor(color_list[bnum])
			ax.get_legend().remove()
			ax.set_xticklabels(tick_label_list,rotation=0)
			ax.set_ylim((0,maxes[item]))
			if item in ['error','strength']:
				plt.ticklabel_format(style='sci',axis='y',scilimits=(0,3))
			ax.tick_params(axis='x',which='both',bottom=False)

	# custom legend for replicates by group
	labels_legend = [i for j in list(zip(*groups)) for i in j]
	n_reps = int(len(labels_legend)/len(colors))
	labels = ['%s'%i for i in groups_names]
	reorder = [j*n_reps+i for i in range(n_reps) 
		for j in range(len(colors))]
	color_list = [color_list[i] for i in reorder]
	# create proxy artists as handles:
	cmap_handles = [Rectangle((0, 0), 1, 1) for _ in labels]
	handler_map = dict(zip(cmap_handles,
		[HandlerColormap(custom_colors=color_list[ii*n_reps:(ii+1)*n_reps]) for ii,i in enumerate(labels)]))
	legend_args = dict(loc='upper right')
	legend = ax.legend(handles=cmap_handles, 
		labels=labels, 
		handler_map=handler_map,
		**legend_args)
	frame = legend.get_frame()
	legend.set_title('replicates')
	frame.set_edgecolor('white')
	frame.set_facecolor('white')

	#! original legend
	if 0:
		# last plot gets a legend
		# note that the bars are ordered by item in the group, then by group so we have to re-reorder
		labels_legend = [i for j in list(zip(*groups)) for i in j]
		labels_legend = [work.metadata.meta[sn].get('name',sn) for sn in labels_legend]
		patches = [mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_list[cnum]) 
			for cnum,__name__ in enumerate(labels_legend)]
		#! assume equal replicates in each group
		n_reps = int(len(labels_legend)/len(colors))
		reorder = [j*n_reps+i for i in range(n_reps) 
			for j in range(len(colors))]
		patches = [patches[i] for i in reorder]
		labels_legend = [labels_legend[i] for i in reorder]
		legend = ax.legend(patches,labels_legend,loc='upper left',
			bbox_to_anchor=(1.0,0.0,1.,1.))
		frame = legend.get_frame()
		frame.set_edgecolor('white')
		frame.set_facecolor('white')
	picturesave('fig.%s%s'%(fn_base,extra_tag),work.plotdir,
		backup=False,version=True,meta={},extras=[legend])

def compute_strengths():
	"""
	Must follow the creation of the report.
	"""
	import scipy.integrate
	import scipy.interpolate
	def integral2d(z,lx,ly,spacer):
		xx,yy = np.meshgrid(np.arange(0,lx,spacer),np.arange(0,ly,spacer))
		fi = scipy.interpolate.interp2d(xx,yy,z,kind='cubic')
		sample_recap = np.array([fi(x,y) for ii,(x,y) 
			in enumerate(zip(xx.reshape(-1),yy.reshape(-1)))]).reshape(xx.shape)
		return scipy.integrate.dblquad(fi,0,lx,lambda x:0,lambda x: ly)
	global report
	for tag in tags:
		for sn in data.this:
			data_ucc.set('curvature_undulation_coupling_dex',
				select=upstreams[tags.index(tag)]['specs'])
			dat = data_ucc.this[sn]['cf_first']
			lx,ly = data.this[sn]['vecs'].mean(axis=0)[:2]
			spacing = data.this[sn]['grid_spacing']
			dA = spacing**2
			# numerical integration with simpson's rule
			#! to improve the resultion you would have to do some interpolation here
			denom = scipy.integrate.simps(scipy.integrate.simps(
				np.ones([i+1 for i in dat.shape])*np.sqrt(dA)))
			# the denominator above should be area, or lx*ly
			ans_numerator = scipy.integrate.simps(scipy.integrate.simps(dat**2*dA))
			#! interpolation is unreliable and does not improve much on Simpson
			if 0:
				print('status interpolating %s (slow)'%sn)
				ans_alt = integral2d(dat**2,lx,ly,spacing)
				"""
				the interpolation method gives very similar answers
				ans_numerator is 26.646672781514038
				ans_alt is 27.066673349948672
				difference is 1.5%
				"""
				ans_numerator = ans_alt[0]
			report[(sn,tag)]['strength'] = ans_numerator/denom

@function
def make_images_and_report():
	"""Make all images and summary bar plot."""
	global plotspec,report
	report = {}
	for tag in tags:
		# save features of each hypothesis set while images are rendered
		result = plot_results(tag=tag,plotspec=plotspec)
		for sn,tag in result:
			report.update(result)
	# plot a summary bar plot from the features we saved during plots
	try: 
		plot_recap(report=report)
		status("maximum error is %.3f"%max([i['error'] for i in report.values()]))
		return report
	#! protect against failure during development
	except: pass	

def package_ucc(dat,sns):
	"""
	Organize the results from one hypothesis.
	"""
	# prepare summary
	# the kappa and gamma (or sigma) are the first two values of the fit ('x')
	params = dict([(sn,dict(list(zip(('kappa','gamma'),
		dat[sn]['x'][:2])))) for sn in sns])
	fits = pd.DataFrame(params).T
	if 0:
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

def calculate_undulations_wrapper(sn,**kwargs):
	"""
	Fit the undulations for both the undulations survey and the entropy calculation.
	"""
	global data,calculate_undulations,data_average_normal
	fit_style = kwargs.pop('fit_style')
	residual_form = kwargs.pop('residual_form')
	midplane_method = kwargs.pop('midplane_method')
	fit_tension = kwargs.pop('fit_tension')
	lims = kwargs.pop('lims')
	bin_size = kwargs.pop('bin_size',None)
	if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
	if midplane_method=='average_normal': 
		data.set('undulations_average_normal')
		data_average_normal = data.this
		#! custom_heights = data_average_normal[sn]['data']['average_normal_heights']
		custom_heights = data_average_normal[sn]['average_normal_heights']
	else: custom_heights = None
	#! updating for datapacK dat = data[sn]['data']
	data.set('import_readymade_meso_v1_membrane')
	dat = data.this[sn]
	mesh = dat['mesh']
	vecs = dat['vecs']
	surf = np.mean(dat['mesh'],axis=0)
	#! when prepping the data structure you cannot forget to remove the zero mode!
	#!   this was an oversight until we noticed a nonzero zeroth mode
	#! the following gives axes don't match array:
	#!   surf -= np.tile(surf.reshape(len(surf),-1).mean(axis=1),
	#!   (surf.shape[0],surf.shape[1])).transpose((2,0,1))
	surf -= np.tile(surf.reshape(len(surf),-1).mean(axis=1),
		(surf.shape[1],surf.shape[2],1)).transpose((2,0,1))
	uspec = calculate_undulations(surf,vecs,
		fit_style=fit_style,custom_heights=custom_heights,lims=lims,
		midplane_method=midplane_method,residual_form=residual_form,fit_tension=fit_tension,
		bin_size=bin_size)
	return uspec

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

if __name__=='__main__' and 0:

	"""
	USAGE
	You can develop code here. Run it with `main` after saving the file.
	Run `make_images_and_report` to remake images and summary bar plots.
	"""

	do_full_plots = False
	do_meta_comparison = False
	do_combined_undulation_plot = False
	do_meta_comparison_bars = False

	# ENTROPY METHOD 1 (set True/1 to try this)
	if do_meta_comparison or 'meta_results' not in globals():

		# wrapping the calculation of results in a loop over meta hypotheses
		meta_results = {}
		# enumerate the meta hypotheses
		hypos = [
			{'fit_style': 'band,perfect,curvefit', 
				'midplane_method': 'average_normal', 'lims': (0.0, 0.1), 
				'residual_form': 'log', 'fit_tension': True},
			{'fit_style': 'band,perfect,curvefit', 
				'midplane_method': 'average_normal', 'lims': (0.0, 1.0), 
				'residual_form': 'log', 'fit_tension': True},
			{'fit_style': 'band,perfect,curvefit', 
				'midplane_method': 'flat', 'lims': (0.0, 1.0), 
				'residual_form': 'log', 'fit_tension': True},
			{'source':'package_ucc'}]
		# once you select a hypothesis for publication, make the full
		#   undulation plot below by choosing underneath where it says
		#   "make the midplane method and limit choices here"
		hypo_names = ['average normal','average normal 0-1','flat 0-1','UCC']
		for hnum,hypo in enumerate(hypos):

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
				#! note that this is ignored because we get the results from the ucc
				lims = (0.,1.0)
				# select a midplane method
				midplane_method = [
					'flat','average','average_normal',
					][0]

				# implement the hypotheses here
				# the following flag specifies whether we are using the coupling results
				#   that come from package_ucc or calculating the no-coupling data on the fly
				hypo_ucc = hypo.get('source',False)=='package_ucc'
				if not hypo_ucc:
					lims = hypo['lims']
					midplane_method = hypo['midplane_method']


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
				if not wrong_method: 
					uspec['energy_raw'] = uspec_no_coupling['energy_raw']

				### END INTERVENTION

				"""
				recall that the method above is called the "wrong" method because if you do not get the energy values
				from the calculate_undulations code, then it comes from the uspec from the result from package_ucc which
				means that the values are really energy values
				adding the following line helps to ensure that the right undulation spectra values are contributing to the entropy
				!!! beware that this might change the coupling results. audit the script to be sure
				"""
				if not hypo_ucc: 
					uspec = uspec_no_coupling

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
					if qs[j] < 0.184:
						qs[j] = qs[j]
					else: qs[j] = 0.0

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

			# archive this set of results
			meta_results[hnum] = copy.deepcopy(results)

		# the following plots require 
		results = meta_results[hypo_names.index('average normal')]

	if 'report' not in globals() or do_full_plots:
		make_images_and_report()
		compute_strengths()

	if do_full_plots:
		plot_recap(report,subset=True)
		plot_recap(report)
		plot_recap(report,subset=True,which='strength')

	if do_combined_undulation_plot:

		"""
		Composite undulation plot.
		Under development for the manuscript. Show the height undulations 
		with various fits, even if they include curvature.
		"""

		subset = True
		#! order must match the colors
		sns_this = ['RigidNCE1', 'RigidNCE2', 'RigidNCE3', 
			'SemiRigidE1_50', 'SemiRigidE2_50', 'SemiRigidE3_50', 
			'CL160ENS-1', 'CL160ENS-2', 'CL160ENS-3']
		#! sns_this = ['CL160ENS-1', 'CL160ENS-2', 'CL160ENS-3']
		groups = [
			('RigidNCE1','RigidNCE2','RigidNCE3',),
			('SemiRigidE1_50','SemiRigidE2_50','SemiRigidE3_50'),
			('CL160ENS-1','CL160ENS-2','CL160ENS-3',),]
		groups_names = ['Rigid','Rigid-Tethered','Flexible']
		group_to_super = {}
		for s,group in zip(['rigid','semi','flex'],groups):
			for g in group: group_to_super[g] = s
		colors = ['#e41a1c','#377eb8','#4daf4a']
		color_list = color_vary(colors=colors,n=3,p=0.2)
		super_colors = {'rigid':'r','semi':''}
		from ortho import sweeper,status
		from omni.plotter.panels import square_tiles
		# make the midplane method and limit choices here (see hypos above)
		midplane_method = ['flat','average','average_normal'][1]
		high_cutoff_undulation = 1.0
		from calcs.codes.undulate_plot import add_undulation_labels,add_axgrid,add_std_legend
		def hqhq(q_raw,kappa,sigma,area,exponent=4.0):
			return 1.0/(area/2.0*(kappa*q_raw**(exponent)+sigma*q_raw**2))
		sns = work.sns()
		sn = sns[0]
		art = {'fs':{'legend':8}}
		lims = (0.0,high_cutoff_undulation)
		plotspec = {
			'fit_style':['band,perfect,curvefit','band,perfect,fit'][0],
			'midplane_method':midplane_method,
			'lims':[(0.0,high_cutoff_undulation),(0.04,high_cutoff_undulation)][0],
			'residual_form':['log','linear'][0],
			'fit_tension':[False,True][1]}
		uspecs = dict([(sn,calculate_undulations_wrapper(sn=sn,**plotspec))
			for sn in sns_this])
		axes,fig = square_tiles(1,figsize=8,favor_rows=True,wspace=0.4,hspace=0.1)
		ax = axes[0]
		for style in ['std','ucc']:
			for snum,(sn,uspec) in enumerate(uspecs.items()):

				if style=='ucc':
					results = meta_results[hypo_names.index('UCC')]
					# get the true fitted kappa and sigma from each method
					uspecs[sn]['kappa'] = results[sn]['kappa'] #! *2 #! WHY?
					#! factor of two!!!
					#! [uspecs[sn]['kappa'] for sn in sns_this]
					uspecs[sn]['gamma'] = uspecs[sn]['sigma'] = results[sn]['sigma']
					uspecs[sn]['area'] = uspecs[sn]['area']
					#! uspecs[sn]['crossover'] = 0.03
					#!!! now it is in the right spot

				q_binned,energy_binned = uspec['q_binned'],uspec['energy_binned']
				# omitting the binned plot here, which might have more regular spacing?
				q_fit,energy_fit = np.transpose(uspec['points'])
				#! plot either the cutoff or complete spectra
				if 0: ax.plot(q_fit,energy_fit,
					'.',lw=0,markersize=10,markeredgewidth=0.5,markeredgecolor='k',
					c=color_list[snum],label=None,zorder=4)
				else: ax.plot(q_binned,energy_binned,
					'.',lw=0,markersize=10,markeredgewidth=0.5,markeredgecolor='k',
					c=color_list[snum],label=None,zorder=4,alpha=1.0)
				#! override q_fit to show the full spectrum for the undulation coupling method
				#!   note that this reflects the fact that the limit was probably too high on that method
				if style=='ucc': q_fit = np.unique(q_binned)
				if 'linear_fit_in_log' in uspec:
					exponent = uspec['linear_fit_in_log']['c0']
					status('alternate exponent %.3f'%exponent,tag='note')
					ax.plot(q_fit,hqhq(q_fit,kappa=uspec['kappa'],sigma=uspec['sigma'],
						exponent=-1.0*exponent,area=uspec['area']),lw=3,zorder=5,c='g')
				elif 'crossover' in uspec:
					ax.plot(q_fit,hqhq(q_fit,kappa=uspec['kappa'],sigma=0.0,
						area=uspec['area']),lw=3,zorder=1,c=color_list[snum])
					# find the first value above the crossover and shift the tension down
					#! the following is not working properly
					ind = np.where(q_fit>uspec['crossover'])[0][0]
					fac = abs((hqhq(q_fit,kappa=0,sigma=uspec['sigma'],
						area=uspec['area'])/hqhq(q_fit,kappa=uspec['kappa'],
						sigma=0,area=uspec['area']))[ind])
					ax.plot(q_fit,fac*hqhq(q_fit,kappa=0,sigma=uspec['sigma'],
						area=uspec['area']),lw=3,zorder=1,c=color_list[snum])
					ax.axvline(uspec['crossover'],c='k',lw=1.0)
				# ignoring the plot styles above
				else:
					# standard exponent plotting method
					ax.plot(q_fit,hqhq(q_fit,kappa=uspec['kappa'],sigma=uspec['sigma'],
						area=uspec['area']),'-' if style=='std' else '--',lw=2,
						zorder=3,c=color_list[snum],alpha=0.5 if style=='std' else 1.0,)
				# omitting the residuals
				add_undulation_labels(ax,art=art)
				add_std_legend(ax,loc='upper right',art=art)
				add_axgrid(ax,art=art)
		ax.set_ylim((10**-5,10**1))
		# custom legend for replicates by group
		labels_legend = [i for j in list(zip(*groups)) for i in j]
		n_reps = int(len(labels_legend)/len(colors))
		labels = ['%s'%i for i in groups_names]
		reorder = [j*n_reps+i for i in range(n_reps) 
			for j in range(len(colors))]
		color_list = [color_list[i] for i in reorder]
		# create proxy artists as handles:
		cmap_handles = [Rectangle((0, 0), 1, 1) for _ in labels]
		# explain the symbols after thecolors
		cmap_handles.append(mpl.lines.Line2D([0],[0],
			marker='.',lw=0,markersize=10,markeredgewidth=0.5,markeredgecolor='k',color='w'))
		cmap_handles.append(mpl.lines.Line2D([0],[0],
			linestyle='-',lw=2,markersize=10,markeredgewidth=1,color='k',alpha=0.5))
		cmap_handles.append(mpl.lines.Line2D([0],[0],
			linestyle='--',lw=2,markersize=10,markeredgewidth=1,color='k'))
		labels += ['observed','S1','S2']
		handler_map = dict(zip(cmap_handles,
			[HandlerColormap(custom_colors=color_list[ii*n_reps:(ii+1)*n_reps]) 
			for ii,i in enumerate(labels[:3])]))
		legend_args = dict(loc='upper right',framealpha=1.0)
		legend = ax.legend(handles=cmap_handles, 
			labels=labels, 
			handler_map=handler_map,
			fontsize=12,
			**legend_args)
		frame = legend.get_frame()
		legend.set_title('replicates')
		frame.set_edgecolor('black')
		frame.set_facecolor('white')
		# select the relevant metadata here by hypothesis name in the manuscript
		#   in which S1 is the no-curvature and S2 is the curvature case
		#   saving the relevant data for each
		# note that we are currently using v5 which came from dextran_undulation_coupling.py
		#   and which probably used the "flat" procedure which does not look very good on the 
		#   undulation spectrum for the no-curvature (S1) case
		#! the above discrepancy (flat vs average, which looks better) should be checked!
		# remove tuple keys from the report
		report_out = dict()
		for tag in tags:
			report_out[tag] = {}
			for sn in sns:
				report_out[tag][sn] = report[(sn,tag)]
		# consistency check: meta_results[0]==results
		meta = {'s1':{'plotspec':plotspec,'results':results},
			's2':{'tag_this':tag_this,'specs':contents_full,'report':report_out}}
		#! small differences in the meta cause increasing versions for no reason
		#! recall that in the first round of plots we looked at three different midplane methods
		#!   however now we have settled on the "correct" versions
		#! resetting the meta here
		meta = {}
		picturesave('fig.undulation_comparison',work.plotdir,
			backup=False,version=True,meta=meta)

	if do_meta_comparison_bars:

		"""
		Accompany the undulation plot with a comparison of the different entropy, kappa, sigma from 
		the hypotheses we are using.
		"""

		# this is mostly a conduit after adding meta_results above
		post_meta = {}
		for hnum,hypo in enumerate(hypos):
			post_meta[hnum] = meta_results[hnum]

		# enumerate the readouts
		tag_this = 'v5'
		readouts = ['kappa','sigma','Entropy_1','Entropy_2']
		axes,fig = panelplot(layout={'out':{'grid':[1,len(readouts)]},
			'ins':[{'grid':[1,1]} for i in readouts]},figsize=(12,4))
		axes = [i for j in axes for i in j]
		hypo_colors = ['#66c2a5','#fc8d62','#8da0cb']
		hypo_colors = ['#e41a1c','#377eb8','#4daf4a']
		hypo_colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3']
		hypo_colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3']
		for rnum,readout in enumerate(readouts):
			ax = axes[rnum]
			for hnum,hypo in enumerate(hypos):
				kwargs = {}
				ax.bar(range(hnum,len(sns)*len(hypos)+hnum,len(hypos)),
					[post_meta[hnum][sn][readout] for sn in sns],
					color=hypo_colors[hnum],width=1.0,
					label=hypo_names[hnum])
			ax.set_xticklabels(sns,rotation=90)
			ax.set_xticks(np.arange(0,len(sns)*len(hypos),len(hypos)))
			ax.set_title(readout)
		legend = axes[-1].legend(bbox_to_anchor=(1.04,1), borderaxespad=0)
		picturesave('fig.meta_comparison',work.plotdir,
			backup=False,version=True,meta={},extras=[legend])
		# package the meta_results in a text file
		meta_results_out = dict([(name,meta_results[hnum]) for hnum,name in enumerate(hypo_names)])
		with open(work.plotdir+'/meta_comparison.json','w') as fp: 
			fp.write(json.dumps(meta_results_out))

if __name__=='__main__':

	"""
	Refactoring NOW.
	"""

	do_survey_one = False
	# the alternateive residual map for do_survey_one uses imshow
	do_alt_residual_map = True
	do_survey_detailed = True
	do_ultimate_comparison = True
	if is_custom:
		work.request_calc = {'curvature_undulation_coupling_dex_1':work.metadata.calculations['curvature_undulation_coupling_dex_1'].get('specs',{}),'curvature_undulation_coupling_dex_2':work.metadata.calculations['curvature_undulation_coupling_dex_2'].get('specs',{})}
		work.calcs.unroll_loops(work.request_calc)

	if do_survey_one:

		"""
		the tags and contents_full variables outline the design space
		select a single simulation and plot the curvature fields and undulations 
		"""
		# assume all is set to a single simulation
		sns = work.metadata.collections['all']
		if len(sns)!=1: 
			print('warning the "all" collection must have one simulation to do the opt survey')
		else:

			def data_ucc_set(specs):
				"""Note that adding flags to the design spec causes the select flag in the set 
				command to fail, hence we have to HACK through it."""
				#! previously: data_ucc.set(select=design_sweep[tag_this]['specs'])
				data_ucc.this = data_ucc.data[np.where([specs==j 
					for i,j in data_ucc.names])[0][0]]
				data_ucc.this = dict([(sn,data_ucc.this[sn]['data']) for sn in data_ucc.this])
	
			sn = sns[0]
			axes,fig = panelplot(layout={'out':{'grid':[2,1]},
				'ins':[{'grid':[1,6]},{'grid':[1,len(tags)]}]},
				figsize=(14,8))
			axes = {
				'error':axes[0][0],'kappa':axes[0][1],'sigma':axes[0][2],
				'undulations':axes[0][4],'spectrum':axes[0][5],'residuals':axes[0][3],
				'fields':axes[1]}
			colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3']
			colors = ['#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0',
				'#f0027f','#bf5b17','#666666']
			colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
				'#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
			colors += colors
			cmap_name = 'RdBu_r'
			scatter_size = 4
			labels = {}
			for tnum,tag_this in enumerate(tags):
				selector = design_sweep[tag_this]['specs']
				data_ucc_set(specs=selector)
				dat = data_ucc.this
				result = package_ucc(dat,sns)
				high_cutoff = design_sweep[tag_this]['specs']['fitting']['high_cutoff']
				qs,energy = result['qs'][sn],result['energy'][sn]
				ax = axes['spectrum']
				opt_method = design_sweep[tag_this]['specs']['design']['optimize_method']
				samples = design_sweep[tag_this]['specs']['design']['samples']
				label = '%s: %s,%s'%(tag_this,
					design_sweep[tag_this]['specs']['design']['optimize_method'],
					{-1:'all'}.get(samples,'%d'%samples))
				#! the labels are too busy so we override it
				label = tag_this
				labels[tag_this] = label
				#!!! hack to rescue the qs because they were not saved properly! for v18. did not help
				#!!! THIS WAS FIXED SO IT CAN BE REMOVED
				if len(qs)!=len(energy):
					from calcs.codes.undulate import perfect_collapser
					perfect_collapser(energy,energy)
					qs = np.unique(qs)
					selector = design_sweep['v02']['specs']
					data_ucc_set(specs=selector)
					dat = data_ucc.this
					result = package_ucc(dat,sns)
					qs2,energy2 = result['qs'][sn],result['energy'][sn]
					qs = qs2
				else: qs_prev = qs
				if selector['design'].get('perfect_collapse_weighted',False):
					from calcs.codes.undulate import perfect_collapser
					qs_col,energy_col,inds = perfect_collapser(qs,energy)
					ax.scatter(qs_col,energy_col,label=label,c=colors[tnum],s=scatter_size,zorder=3)
				else: ax.scatter(qs,energy,label=label,c=colors[tnum],s=scatter_size,zorder=3)
				ax.axvline(high_cutoff,c='k')
				ax.axhline(10**0,c='k')
				ax.set_xscale('log')
				ax.set_yscale('log')
				#! standard residual plot
				if not do_alt_residual_map:
					ax = axes['residuals']
					fitted_inds = np.where(qs<=high_cutoff)
					x,y = qs,energy-10**0
					ax.scatter(x[fitted_inds],y[fitted_inds],label=label,
						c=colors[tnum],s=scatter_size,zorder=3)
					ax.axhline(0,c='k')
					#! ax.set_xscale('log')
			# alternate residual plot
			if do_alt_residual_map:
				# collect the number of fitted points for each hypothesis
				fitted_inds_counts = np.zeros(len(tags))
				for tnum,tag_this in enumerate(tags):
					selector = design_sweep[tag_this]['specs']
					data_ucc_set(specs=selector)
					dat = data_ucc.this
					result = package_ucc(dat,sns)
					high_cutoff = design_sweep[tag_this]['specs']['fitting']['high_cutoff']
					qs,energy = result['qs'][sn],result['energy'][sn]
					fitted_inds = np.where(qs<=high_cutoff)[0]
					fitted_inds_counts[tnum] = len(fitted_inds)
				#! uses the fitted_inds and qs above
				residual_map = np.zeros((len(tags),fitted_inds_counts.max().astype(int)))
				for tnum,tag_this in enumerate(tags):
					data_ucc_set(specs=design_sweep[tag_this]['specs'])
					dat = data_ucc.this
					result = package_ucc(dat,sns)
					qs,energy = result['qs'][sn],result['energy'][sn]
					fitted_inds = np.where(qs<=high_cutoff)[0]
					x,y = qs,np.log10(energy)
					residual_map[tnum][:len(fitted_inds)] = y[fitted_inds]
				ax = axes['residuals']
				max_residual_mag = np.abs(residual_map).max()
				ax.imshow(residual_map,interpolation='nearest',aspect='auto',
					cmap=mpl.cm.__dict__[cmap_name],
					vmax=max_residual_mag,vmin=max_residual_mag*-1)
			# collect errors
			errors,kappas,sigmas = [],[],[]
			for tnum,tag_this in enumerate(tags):
				data_ucc_set(specs=design_sweep[tag_this]['specs'])
				dat = data_ucc.this
				error = data_ucc.this[sn]['bundle'][sn]['fun']
				errors.append(error)
				kappas.append(data_ucc.this[sn]['x'][0])
				do_tension = design_sweep[tag_this]['specs']['design'].get('tension',True)
				sigmas.append(data_ucc.this[sn]['x'][1] if do_tension else 0.0)
			# transform errors to errors per mode
			#!!! check this
			errors = [np.sqrt(i)/len(dat[sn]['qs']) for i in errors]
			ax = axes['error']
			ax.bar(range(len(tags)),errors,color=colors)
			ax.set_ylabel('error per mode')
			ax = axes['kappa']
			ax.bar(range(len(tags)),kappas,color=colors)
			ax = axes['sigma']
			ax.bar(range(len(tags)),sigmas,color=colors)
			legend = axes['spectrum'].legend(bbox_to_anchor=(1.04,1), borderaxespad=0)
			curvatures = {}
			for tnum,tag_this in enumerate(tags):
				data_ucc_set(specs=design_sweep[tag_this]['specs'])
				dat = data_ucc.this
				cf = data_ucc.this[sn]['cf']
				curvatures[tag_this] = cf
			c0_max = max([np.abs(i).max() for i in curvatures.values()])
			for tnum,tag_this in enumerate(tags):
				cf = curvatures[tag_this]
				ax = axes['fields'][tnum]
				# if you do not zero-center this, the background color indicates the net curvature
				#   direction, for example, if it is salmon (red) then the net curvature is more red
				ax.imshow(cf.T,interpolation='nearest',cmap=mpl.cm.__dict__[cmap_name],
					vmax=c0_max,vmin=-1*c0_max)
				ax.set_title(labels.get(tag_this,tag_this))

			sns_this = [sn]
			low_cutoff = 0.0
			color_list = colors
			high_cutoff = design_sweep[tag_this]['specs']['fitting']['high_cutoff']
			# plot undulations (repetitive with code from the detailed survey)
			ax = axes['undulations']
			fit_style = [
				'band,perfect,curvefit',
				'band,perfect,fit',
				'band,blurry,curvefit',
				][-1]
			spec_undulations = dict(fit_style=fit_style,residual_form='log',
				midplane_method='average_normal',fit_tension=True,lims=(low_cutoff,high_cutoff))
			spec_undulations['bin_size'] = 0.01
			uspecs = {}
			for snum,sn in enumerate(sns_this):
				uspec = calculate_undulations_wrapper(sn=sn,**spec_undulations)
				uspecs[sn] = uspec
				color = color_list[snum]
				# remove zero
				q_binned = uspec['q_binned']
				energy_binned = uspec['energy_binned']
				inds = q_binned>0
				x,y = q_binned[inds],energy_binned[inds]
				uspec['points'] = x,y
				ax.scatter(x,y,color='k',s=1,zorder=1)
			ax.set_xscale('log')
			ax.set_yscale('log')
			# plot curvature fits on the undulation plot
			def hqhq(q_raw,kappa,sigma,area,exponent=4.0):
				return 1.0/(area/2.0*(kappa*q_raw**(exponent)+sigma*q_raw**2))
			for tnum,tag_this in enumerate(tags):
				data_ucc_set(specs=design_sweep[tag_this]['specs'])
				uspec = uspecs[sn]
				x,y = uspec['points']
				area = uspec['area']
				do_tension = design_sweep[tag_this]['specs']['design'].get('tension',True)
				y = hqhq(q_raw=x,kappa=data_ucc.this[sn]['x'][0],
					sigma=data_ucc.this[sn]['x'][1] if do_tension else 0.0,
					area=area)
				x,y = x[y>0],y[y>0]
				ax.plot(x,y,linestyle='--',color=color_list[tnum],zorder=tnum,lw=1)
			if 0: ax.plot(x,hqhq(q_raw=x,kappa=uspec['kappa'],sigma=uspec['sigma'],
				area=area),linestyle='-',color=color_list[snum])
			ax.axvline(high_cutoff,c='k')
			# plot undulation spectra with the no-curvature fit
			def hqhq(q_raw,kappa,sigma,area,exponent=4.0):
				return 1.0/(area/2.0*(kappa*q_raw**(exponent)+sigma*q_raw**2))
			# for this survey, there should only be one simulation anyway
			for snum,sn in enumerate(sns_this):
				uspec = uspecs[sn]
				x,y = uspec['points']
				area = uspec['area']
				ax.plot(x,hqhq(q_raw=x,kappa=uspec['kappa'],sigma=uspec['sigma'],
					area=area),linestyle=':',color='k',alpha=0.35)
			picturesave('fig.DEV01',work.plotdir,
				backup=False,version=True,meta={},extras=[legend])

	if do_survey_detailed:

		if do_ultimate_comparison:

			"""
			transplanting the objective function here to get the error
			"""

			from calcs.codes.curvature_coupling.curvature_coupling import InvestigateCurvature

			do_vibe = False
			do_tension = True
			do_collapse = True

			sn = 'CL160ENS-3'
			calc = work.metadata.calculations['curvature_undulation_coupling_dex_2']['specs']
			specs = calc
			data.set('import_readymade_meso_v1_membrane')
			undulations_data = data.this
			data.set('import_readymade_meso_v1_nanogel')
			protein_data = data.this
			design = specs['design']['loop']['v24']
			fitting = specs['fitting']
			# backfill the right keys because InvestigateCurvature
			work.meta = work.metadata.meta

			lowcut,hicut = 0.0,0.2
			do_blurry = False
			do_collapse = True
			do_collapse_weighted = False
			do_positive_tension = False

			signterm = -1.0
			machine_eps = eps = np.finfo(float).eps
			ndrops = 0
			Nfeval = 0

			calcs_dn = os.path.join(os.getcwd(),'calcs')
			if calcs_dn not in sys.path: sys.path.insert(0,calcs_dn)
			ic = InvestigateCurvature(sn=sn,work=work,
				design=design,fitting=fitting,
				protein_abstractor=protein_data[sn],
				undulations=undulations_data[sn],
				do_calculation=False)
			hqs = ic.memory[(sn,'hqs')]
			vecs = ic.memory[(sn,'vecs')]
			# formulate the wavevectors
			lenscale = 1.0
			m,n = mn = np.shape(hqs)[1:]
			Lx,Ly = np.mean(vecs,axis=0)[:2]
			q2d = lenscale*np.array([[np.sqrt(
				((i-m*(i>m/2))/((Lx)/1.)*2*np.pi)**2+
				((j-n*(j>n/2))/((Ly)/1.)*2*np.pi)**2)
				for j in range(0,n)] for i in range(0,m)])
			q_raw = q_raw_expl = np.reshape(q2d,-1)[1:]
			area = (Lx*Ly/lenscale**2)
			if do_collapse:
				from codes.undulate import perfect_collapser
				x = q_raw
				x_red,_,reduce_inds = perfect_collapser(x,x)
				reduce_inds_lazy = [np.where(reduce_inds==i) for i in np.unique(reduce_inds)]
				# override the band forthe new points
				# choosing greater than lowcut to exclude the zero mode
				band = np.where(np.all((x_red>lowcut,x_red<hicut),axis=0))

			dotplace = lambda n : re.compile(r'(\d)0+$').sub(r'\1',"%3.5f"%float(n)).rjust(8)
			def residual(values): 
				return np.sum(np.log10(values.clip(min=machine_eps))**2)/float(len(values))
			residual_local = residual

			def callback(args):
				"""Watch the optimization."""
				global Nfeval
				name_groups = ['kappa',]+(['gamma'] if do_tension else [])+(
					['vibe'] if do_vibe else [])+['curve(%d)'%i for i in range(ndrops)]
				text = ' step = %d '%Nfeval+' '.join([name+' = '+dotplace(val)
					for name,val in list(zip(name_groups,args))+[('error',objective(args))]])
				status('searching! '+text,tag='optimize')
				Nfeval += 1

			def multipliers(x,y): 
				"""Multiplying complex matrices in the list of terms that contribute to the energy."""
				return x*np.conjugate(y)

			def objective(args,mode='residual'):
				"""
				Fit parameters are defined in sequence for the optimizer.
				They are: kappa,gamma,vibe,*curvatures-per-dimple.
				"""
				if do_vibe and do_tension:
					(kappa,gamma,vibe),curvatures = args[:3],args[3:]
				elif do_vibe and not do_tension:
					(kappa,vibe),curvatures = args[:2],args[2:]
				elif not do_vibe and do_tension:
					(kappa,gamma),curvatures = args[:2],args[2:]
				elif not do_vibe and not do_tension:
					(kappa,),curvatures = args[:1],args[1:]
				else: raise Exception('unclear request')
				if not do_tension: gamma = 0.0
				#! composite = self.curvature_sum(cfs,curvatures,method=curvature_sum_method)
				#! cqs = cctools.fft_field(composite)
				termlist = [multipliers(x,y) for x,y in [(hqs,hqs),]]
				termlist = [np.reshape(np.mean(k,axis=0),-1)[1:] for k in termlist]
				# skipping assertion and dropping imaginary
				termlist = [np.real(k) for k in termlist]
				# reduce the termlist if we are doing any binning
				if do_blurry or (do_collapse and not do_collapse_weighted):
					termlist = [np.array([termlist[i][r].mean() 
						for r in reduce_inds_lazy]) for i in range(len(termlist))]
					q_raw = x_red
				# if you refer to q_raw in a conditional below and it exists in the parent
				#   scope as well, then whenever you skip the conditional, it will not be defined
				#   because python knows it is a local and hence ignores the parent scope. it is
				#   best to be very careful with scope
				else: q_raw = q_raw_expl
				if do_positive_tension: gamma = np.abs(gamma)
				if 0: hel = (kappa/2.0*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
					+signterm*termlist[2]*q_raw**2+termlist[3])
					+gamma*area*(termlist[0]*q_raw**2))
				hel = (kappa/2.0*area*(termlist[0]*q_raw**4)+gamma*area*(termlist[0]*q_raw**2))
				if do_vibe:
					ratio = hel/((vibe*q_raw+machine_eps)/(np.exp(vibe*q_raw)-1)+machine_eps)
				else: ratio = hel
				if mode=='residual': return residual_local(ratio[band])
				elif mode=='ratio': return ratio
				else: raise Exception('invalid mode %s'%mode)

			args = data_ucc.this[sn]['x'][0],data_ucc.this[sn]['x'][1]
			error_alt = objective(args)
			error_real = data_ucc.this[sn]['bundle'][sn]['fun']
			method_opt = ic.design.get('optimize_method','BFGS')
			# select an optimization method
			methods_opt = dict([
				('SLSQP',{'maxiter':1000}),
				('BFGS',{})])
			method_args = methods_opt[method_opt]
			fit = scipy.optimize.minimize(objective,
				x0=tuple(args),method=method_opt,
				callback=callback,options=method_args)

		groups = {
			'Rigid':('RigidNCE1','RigidNCE2','RigidNCE3'),
			'Rigid-Tethered':('SemiRigidE1_50','SemiRigidE2_50','SemiRigidE3_50'),
			'Flexible':('CL160ENS-1','CL160ENS-2','CL160ENS-3',),}
		groups_order = ['Rigid','Rigid-Tethered','Flexible']
		colors = ['#e41a1c','#377eb8','#4daf4a']
		color_list = color_vary(colors=colors,n=3,p=0.2)

		def dextran_bar(ax,data,color_list,
			groups_order,sci_y=False,negatives=False,do_legend=True,
			legend_loc='upper right',bbox_legend=False):
			"""
			Stylized plot for the dextran manuscript.
			"""
			df = pd.DataFrame(data)
			plot = df.plot.bar(ax=ax,edgecolor='black',stacked=False)
			if negatives:
				df_below = pd.DataFrame(negatives)
				df.plot.bar(ax=ax)
			#! notes on "replicate-grouped bars": we have to set the bars manually
			#!   because pandas expects each member of a group to be the same color
			#!   for example repeats of sets of three colors for each of three different
			#!   slices of the original dataframe. this is not useful for replicates
			bars = [this for this in plot.get_children()
				if this.__class__.__name__=='Rectangle' and this._width<1]
			if len(color_list)<len(bars):
				raise Exception('we have %d bars and only %d colors'%(len(bars),len(color_list)))
			for bnum,bar in enumerate(bars):
				bar.set_facecolor(color_list[bnum])
			ax.get_legend().remove()
			ax.set_xticklabels(groups_order,rotation=0)
			if sci_y:
				plt.ticklabel_format(style='sci',axis='y',scilimits=(0,3))
			ax.tick_params(axis='x',which='both',bottom=False)
			# custom legend for replicates by group
			if do_legend:
				n_reps = list(set([len(i) for i in data]))
				if len(n_reps)!=1: raise Exception('non-uniform number of replicates')
				else: n_reps = n_reps[0]
				labels = ['%s'%i for i in groups_order]
				reorder = [j*n_reps+i for i in range(n_reps) 
					for j in range(len(colors))]
				color_list = [color_list[i] for i in reorder]
				# create proxy artists as handles:
				cmap_handles = [Rectangle((0, 0), 1, 1) for _ in labels]
				handler_map = dict(zip(cmap_handles,
					[HandlerColormap(custom_colors=color_list[ii*n_reps:(ii+1)*n_reps]) 
					for ii,i in enumerate(labels)]))
				legend_args = dict(loc=legend_loc)
				if bbox_legend: legend_args['bbox_to_anchor'] = bbox_legend
				legend = ax.legend(handles=cmap_handles, 
					labels=labels, 
					handler_map=handler_map,
					**legend_args)
				frame = legend.get_frame()
				legend.set_title('replicates')
				frame.set_edgecolor('white')
				frame.set_facecolor('white')
				return legend

		cmap_name = 'RdBu_r'

		# select the design
		"""
		note on loops:
			...!!!
		"""
		if 0:
			#!!!! did we ever test v07 !!!??? or is this the first time?
			tag_this = 'v02'
			data_ucc.set(select=design_sweep[tag_this]['specs'])
			high_cutoff = design_sweep[tag_this]['specs']['fitting']['high_cutoff']
			low_cutoff = 0.0

			tag_this = 'v21'
			data_ucc.set(select=design_sweep[tag_this]['specs'])
			high_cutoff = design_sweep[tag_this]['specs']['fitting']['high_cutoff']
			low_cutoff = 0.0
			data_ucc_fold = this = work.plotload(calcname+'_fold')

		sns_this = [i for j in groups.values() for i in j]
		nrows_grid = ncols_grid = int(np.ceil(np.sqrt(len(sns_this))))
		axes,fig = panelplot(layout={'out':{'grid':[2,3],
			'hspace':0.2,'wspace':0.2},
			'ins':[
				{'grid':[nrows_grid,ncols_grid],'hspace':0.5,'wspace':0.5},
				{'grid':[2,2],'hspace':0.5,'wspace':0.5},{'grid':[0,0]},
				{'grid':[nrows_grid,ncols_grid],'hspace':0.5,'wspace':0.5},
				{'grid':[1,1]},{'grid':[1,1]},]},
			figsize=(20,14))
		axes = {
			'kappa':axes[1][0],'gamma':axes[1][1],
			'error':axes[1][2],'strengths':axes[1][3],
			'fields':axes[0],'z':axes[3],'energy':axes[4][0],
			'undulations':axes[5][0]}
		for a in range(nrows_grid*ncols_grid-1,len(sns)-1,-1):
			fig.delaxes(axes['fields'][a])

		def get_kappa(sn):
			return data_ucc.this[sn]['x'][0]
		def get_gamma(sn):
			return data_ucc.this[sn]['x'][1]
		def get_error(sn):
			return data_ucc.this[sn]['bundle'][sn]['fun']
		def get_c0_max(sn):
			return data_ucc.this[sn]['cf'].max()
		def get_c0_min(sn):
			return data_ucc.this[sn]['cf'].min()

		# bar plots of relevant readouts
		extras = []
		for key in ['kappa','gamma','error']:
			ax = axes[key]
			ax.set_title({'kappa':'bending rigidity','gamma':'tension','error':'error'}[key])
			ax.set_ylabel({
				'kappa':r'$\mathrm{\kappa\,(k_B T)}$',
				'gamma':r'$\mathrm{\sigma\,(k_B T {nm}^{-2})}$','error':'MSE'}[key])
			this = [[globals()['get_%s'%key](sn) 
				for sn in groups[g]] for g in groups_order]
			if key=='gamma':
				kwargs['bbox_legend'] = (1.0,0.0,1.,1.)
				kwargs['legend_loc'] = 'upper left'
				kwargs['do_legend'] = True
			else: kwargs = {'do_legend':False}
			out = dextran_bar(ax=ax,data=this,
				color_list=color_list,groups_order=groups_order,**kwargs)
			if out: extras.append(out)

		# plot C_0 histograms
		sb.violinplot(data=[data_ucc.this[sn]['cf'].reshape(-1) 
			for sn in sns_this],ax=axes['strengths'])
		ax = axes['strengths']
		ax.set_ylabel(r'$C_0\,(nm)$')
		ax.set_title('curvature')

		# plot curvature fields
		curvatures = {}
		for snum,sn in enumerate(sns_this):
			cf = data_ucc.this[sn]['cf']
			curvatures[sn] = cf
		c0_max = max([np.abs(i).max() for i in curvatures.values()])
		for snum,sn in enumerate(sns_this):
			cf = curvatures[sn]
			ax = axes['fields'][snum]
			ax.imshow(cf.T,interpolation='nearest',cmap=mpl.cm.__dict__[cmap_name],
				vmax=c0_max,vmin=-1*c0_max)
			ax.set_title(sn)
			ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=True)
			ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)

		# plot mean surfaces
		z_surfaces = {}
		for snum,sn in enumerate(sns_this):
			data.set('import_readymade_meso_v1_membrane')
			mesh = data.this[sn]['mesh']
			surf = np.mean(mesh,axis=0).mean(axis=0)
			surf -= surf.mean()
			z_surfaces[sn] = surf
		z_max = max([np.abs(i).max() for i in z_surfaces.values()])
		for snum,sn in enumerate(sns_this):
			ax = axes['z'][snum]
			surf = z_surfaces[sn]
			ax.imshow(surf.T,interpolation='nearest',cmap=mpl.cm.__dict__[cmap_name],
				vmax=z_max,vmin=-1*z_max)
			ax.set_title(sn)
			ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=True)
			ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)

		# plot energy spectra
		ax = axes['energy']
		ax.set_title('energy spectra')
		#ax.set_ylabel(r'$\mathrm{\mathcal{H_{el}}}$')
		#ax.set_xlabel(r'$\mathrm{\vec{q}}$')
		for snum,sn in enumerate(sns_this):
			qs = data_ucc.this[sn]['qs']
			energy = data_ucc.this[sn]['ratios']
			ax.scatter(qs,energy,color=color_list[snum])
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.axhline(10**0,c='k')
		#!!! hardcoded because we have mangled the data to fit it in this script
		if is_custom: 
			high_cutoff,low_cutoff = 0.2,0.0
		else: #! added this for pixel ...!!!
			high_cutoff,low_cutoff = 0.2,0.0
		ax.axvline(high_cutoff,c='k')

		#! sns_this = sns_this[:1]
		# plot undulations
		zero_limit = 10**-5 # enforce no tension-drop below zero
		ax = axes['undulations']
		ax.set_xlabel(r'$\mathrm{\vec{q}}$')
		#ax.set_ylabel(r'$\mathrm{ \langle h_{\vec{q}} h_{\vec{q^{\prime}}} \rangle}$')
		ax.set_title('undulations')
		fit_style = [
			'band,perfect,curvefit',
			'band,perfect,fit',
			'band,blurry,curvefit',
			][-1]
		spec_undulations = dict(fit_style=fit_style,residual_form='log',
			midplane_method='average_normal',fit_tension=True,lims=(low_cutoff,high_cutoff))
		spec_undulations['bin_size'] = 0.01
		uspecs = {}
		for snum,sn in enumerate(sns_this):
			uspec = calculate_undulations_wrapper(sn=sn,**spec_undulations)
			uspecs[sn] = uspec
			color = color_list[snum]
			# remove zero
			q_binned = uspec['q_binned']
			energy_binned = uspec['energy_binned']
			inds = q_binned>0
			x,y = q_binned[inds],energy_binned[inds]
			uspec['points'] = x,y
			ax.scatter(x,y,color=color)
		ax.set_xscale('log')
		ax.set_yscale('log')
		# plot fits
		def hqhq(q_raw,kappa,sigma,area,exponent=4.0):
			return 1.0/(area/2.0*(kappa*q_raw**(exponent)+sigma*q_raw**2))
		for snum,sn in enumerate(sns_this):
			uspec = uspecs[sn]
			x,y = uspec['points']
			area = uspec['area']
			y_this = hqhq(q_raw=x,kappa=get_kappa(sn),sigma=get_gamma(sn),
				area=area)
			inds = y_this>zero_limit
			ax.plot(x[inds],y_this[inds],linestyle='--',color=color_list[snum])
			ax.plot(x,hqhq(q_raw=x,kappa=uspec['kappa'],sigma=uspec['sigma'],
				area=area),linestyle='-',color=color_list[snum])
		ax.axvline(high_cutoff,c='k')

		picturesave('fig.DEV02',work.plotdir,
			backup=False,version=True,meta={},extras=extras)
