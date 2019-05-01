#!/usr/bin/env python

from omni import WorkSpace
import calcs
from calcs.codes import undulate
from calcs.codes.undulate import calculate_undulations
from calcs.codes.undulate_plot import undulation_panel
import numpy as np
from ortho import sweeper,status
from omni.legacy.panels import panelplot
from omni.base.store import picturesave
import re

# keys for the specs file
plotname = 'undulations'

@loader
def load():
	work = WorkSpace(analysis=True)
	data = work.plotload(plotname=plotname)
	sns = work.sns()	

if __name__=='__main__':

	if False:

		# wavevector limits
		lims = (0.,1.0)

		# select a simulation from sns
		sn = 'CL160ENS-1'
		# select a midplane method
		midplane_method = [
			'flat','average','average_normal',
			][0]

		if midplane_method=='average_normal':
			data.set('undulations_average_normal')
			dat = data.this[sn]
			mesh = dat['mesh']
			custom_heights = None
		else: 
			data.set('import_readymade_meso_v1_membrane')
			dat = data.this[sn]
			mesh = dat['mesh']
			custom_heights = None

		vecs = dat['vecs']
		# assume bilayer, even if one mesh, then take the average
		surf = np.mean(mesh,axis=0)
		# kernel of this plot/calculation: calculate the spectra here
		uspec = calculate_undulations(surf,vecs,chop_last=True,custom_heights=custom_heights,
			perfect=True,lims=lims,raw=False,midplane_method=midplane_method)

		layout = {'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]}
		axes,fig = panelplot(layout,figsize=(8,8))
		ax = axes[0]

		kB = 1.38E-23
		temp = 300.00
		mass = 1.44E-22
		area_mem = 0.25E-12 

		kappa = uspec['kappa']
		sigma = uspec['sigma']
		qs = uspec['q_raw']
		hqs = uspec['energy_raw']

		w_q1 = np.sqrt(area_mem*(kappa*qs**4+sigma*qs**2)/mass)
		w_q2 = np.sqrt(kB*temp/(hqs**2*mass)) 

	#! from: edit ../../../retired/factory/calc/dextran/calcs/plot-dextran_analyze_undulations.py
	if 1:

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

			uspec = calculate_undulations(surf,vecs,
				fit_style=fit_style,custom_heights=custom_heights,lims=lims,
				midplane_method=midplane_method,residual_form=residual_form,fit_tension=fit_tension)
			return uspec


		#! more imports
		from omni.plotter.panels import square_tiles

		#---high cutoff for undulations since it lacks this setting
		high_cutoff_undulation = 0.1
		#---subjects of the analysis: different ways to fit the undulations
		manyspectra = sweeper(**{
			'fit_style':['band,perfect,curvefit','band,perfect,fit'][:1],
			'midplane_method':['flat','average','average_normal'],
			'lims':[(0.0,high_cutoff_undulation),(0.04,high_cutoff_undulation)],
			'residual_form':['log','linear'][:1],
			'fit_tension':[False,True]})

		#! from the spectra_comparison section

		#! from codes.undulate_plot import add_undulation_labels,add_axgrid,add_std_legend
		from calcs.codes.undulate_plot import add_undulation_labels,add_axgrid,add_std_legend

		def hqhq(q_raw,kappa,sigma,area,exponent=4.0):
			return 1.0/(area/2.0*(kappa*q_raw**(exponent)+sigma*q_raw**2))

		sns = work.sns()
		sn = sns[0]
		art = {'fs':{'legend':8}}
		lims = (0.0,high_cutoff_undulation)
		plotspecs = manyspectra
		axes,fig = square_tiles(len(plotspecs),figsize=18,favor_rows=True,wspace=0.4,hspace=0.1)
		for pnum,plotspec in enumerate(plotspecs):
			ax = axes[pnum]

			#! previously sn was not required
			uspec = calculate_undulations_wrapper(sn=sn,**plotspec)

			label = 'structure: %s'%re.sub('_',' ',plotspec['midplane_method'])
			label += '\n residuals: %s'%plotspec['residual_form']
			label += '\n method: \n%s'%plotspec['fit_style']
			label += '\n'+r'$\mathrm{\kappa='+('%.1f'%uspec['kappa'])+'\:k_BT}$'
			if uspec['sigma']!=0.0:
				label += '\n'+r'$\mathrm{\sigma='+('%.3f'%uspec['sigma'])+'\:{k}_{B} T {nm}^{-2}}$'
			colors = ['b','k','r']

			q_binned,energy_binned = uspec['q_binned'],uspec['energy_binned']
			ax.plot(q_binned,energy_binned,'.',lw=0,markersize=10,markeredgewidth=0,
				c=colors[0],label=None,alpha=0.2)
			q_fit,energy_fit = np.transpose(uspec['points'])
			ax.plot(q_fit,energy_fit,'.',lw=0,markersize=4,markeredgewidth=0,
				c=colors[1],label=label,alpha=1.,zorder=4)
			#---alternate exponent
			if 'linear_fit_in_log' in uspec:
				exponent = uspec['linear_fit_in_log']['c0']
				status('alternate exponent %.3f'%exponent,tag='note')
				ax.plot(q_fit,hqhq(q_fit,kappa=uspec['kappa'],sigma=uspec['sigma'],
					#---distinctive green color for alternate exponents
					exponent=-1.0*exponent,area=uspec['area']),lw=3,zorder=5,c='g')
			elif 'crossover' in uspec:
				ax.plot(q_fit,hqhq(q_fit,kappa=uspec['kappa'],sigma=0.0,
					area=uspec['area']),lw=3,zorder=5,c='g')
				ax.plot(q_fit,hqhq(q_fit,kappa=0,sigma=uspec['sigma'],
					area=uspec['area']),lw=3,zorder=5,c='g')
				ax.axvline(uspec['crossover'],c='k',lw=1.0)
			#---standard exponent
			else:
				ax.plot(q_fit,hqhq(q_fit,kappa=uspec['kappa'],sigma=uspec['sigma'],
					area=uspec['area']),lw=1,zorder=5,c=colors[2])
				#---inset axis shows the residuals
				#---! note that there is a big whereby the log tick marks still show up sometimes
				from mpl_toolkits.axes_grid1.inset_locator import inset_axes
				axins = inset_axes(ax,width="30%",height="30%",loc=3)
				diffs = (energy_fit-
					hqhq(q_fit,kappa=uspec['kappa'],sigma=uspec['sigma'],area=uspec['area']))
				axins.set_title('residuals',fontsize=6)
				if plotspec['residual_form']=='log': 
					axins.plot(10**diffs,'.-',lw=0.5,ms=1,color='k')
					axins.set_xscale('log')
					axins.set_yscale('log')
					axins.axhline(1.0,c='k',lw=0.5,alpha=0.5)
				else: 
					axins.axhline(0.0,c='k',lw=0.5,alpha=0.5)
					axins.plot(diffs,'.-',lw=0.5,ms=1,color='k')
				axins.set_xticks([])
				axins.set_yticks([])
				axins.set_xticklabels([])
				axins.set_yticklabels([])
				axins.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
				axins.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			add_undulation_labels(ax,art=art)
			add_std_legend(ax,loc='upper right',art=art)
			add_axgrid(ax,art=art)

		picturesave('fig.undulation_survey',work.plotdir,backup=False,version=True,meta={})
