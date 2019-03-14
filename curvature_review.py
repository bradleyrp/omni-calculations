#!/usr/bin/env python

"""
INSPECT THE CURVATURE COUPLING RESULTS
"""

#! macos is fucked
import matplotlib as mpl
mpl.use("TkAgg")


#! somewhat counterintuitive to do this directly. more elegant solution
from calcs.codes import curvature_coupling # .InvestigateCurvature
from calcs.codes.curvature_coupling import InvestigateCurvature

#! this was part of some thing where we hooked the plotload? not sure
#! codes.curvature_coupling.InvestigateCurvature.plotload = plotload
#! codes.curvature_coupling.InvestigateCurvature.plotname = plotname
#! codes.curvature_coupling.InvestigateCurvature.work = work

# alternate loader
from calcs.dextran_importer import plotloader_for_dextran

if __name__=='__main__':
	ic = codes.curvature_coupling.InvestigateCurvature.\
		InvestigateCurvature(plotloader=plotloader_for_dextran)

	if False:
		art = {'fs':{'legend':10}}
		labels = dict([(i,i) for i in work.sns()])
		data,calc = plotloader_for_dextran('curvature_dextran')
		from codes.undulate_plot import undulation_panel
		axes,fig = panelplot(layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]},figsize=(10,10))
		ax = axes[0]
		uspecs = undulation_panel(ax,data['undulations'],art=art,labels=labels)
		picturesave('fig.undulation_dev',work.plotdir,backup=False,version=True,meta={})
