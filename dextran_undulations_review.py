#!/usr/bin/env python

"""
Refactored curvture-undulation coupling code.
"""

from omni import WorkSpace
import calcs
from calcs.codes.curvature_coupling.curvature_coupling_plots \
	import individual_reviews_plotter
#! is the gopher a proto-hook?
from omni.base.tools import gopher

plotname = 'curvature_undulation_coupling'

@loader
def load():
	work = WorkSpace(analysis=True)
	#! we always get plotname? where is this documented
	plotspecs = work.metadata.plots[plotname].get('specs',{})
	calcname = plotspecs.get('calcname',plotname)

if __name__=='__main__':
	data = work.plotload(plotname=plotname)
	sns = sns_ordered = work.sns()
	#! put this in loader?
	#alt_loader = plotspecs.get('loader',None)
	#if alt_loader:
	#	postdat = gopher(alt_loader)(data=data)

	# make the plots
	if 0:

		#! previous function
		"""
		Plot a few versions of the main figure.
		"""
		plotspec = {
			'coupling_review':{
				'viewnames':['average_height','average_height_pbc','neighborhood_static',
					'neighborhood_dynamic','average_field','example_field','example_field_pbc',
					'spectrum','spectrum_zoom']},
			'coupling_review.simple':{
				'viewnames':['average_height','example_field'],'figsize':(6,6)},
			'coupling_review.center_debug':{
				'viewnames':['neighborhood_static','neighborhood_dynamic',
					'average_height','average_height_pbc','average_field','example_field','example_field_pbc',
					'average_height_center','average_field_center','average_field_center_smear','curvature_field_center'],
					'figsize':(16,16)},
			'coupling_review.simple_centered':{
				'viewnames':['average_height_center','curvature_field_center',
				'coupling_review.simple'],
				'figsize':(8,8),'horizontal':True,'wspace':0.7},}
		# turn some off when developing
		for i in [
			'coupling_review.center_debug','coupling_review.simple_centered'
			]: plotspec.pop(i)
		#global seepspace
		#! seepspace = 'data,datas,postdat,undulations_name,calcs,protein_abstractor_name'.split(',')
		#! seep = dict([(key,globals()[key]) for key in seepspace])
		#! custom name?
		data.set(plotspecs['calcname'])
		#! reformulate datas here which is a loop over items in a sweep
		#! cannot remember easier way how to get the real key out of things so change it manually
		tag = work.metadata.plots['curvature_undulation_coupling'][
			'calculation']['curvature_undulation_coupling_dex']['design']
		datas = {tag:data.this}
		#! merging the data
		#! make the merge into a fundamental operation for backwards compatibility?
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
		for sn in sns: postdat[sn]['points_protein_mean'] = postdat[sn]['points_protein'].mean(axis=1).mean(axis=0)
		seep = dict(data=data,datas=datas,work=work,postdat=postdat,
			undulations_name=plotspecs['undulations_name'],
			protein_abstractor_name=plotspecs['protein_abstractor_name'])
		for out_fn,details in plotspec.items(): 
			report_this = individual_reviews_plotter(out_fn=out_fn,seep=seep,**details)
			#! save the report from the full plot not sure if the data is there otherwise
			if out_fn=='coupling_review': report = report_this

	# summary bar plot
	if 1:

		rv3 = {('SemiRigidE1', 'v3'): {'max_curvature': 0.14138229165211977, 'error': 0.009939512501818144}, ('SemiRigidE2', 'v3'): {'max_curvature': 0.09590345687591227, 'error': 0.012985067764987571}, ('SemiRigidE3', 'v3'): {'max_curvature': 0.14138229165211977, 'error': 0.009939512501818144}, ('RigidNCE1', 'v3'): {'max_curvature': 0.08306147103292033, 'error': 0.009246283692621288}, ('RigidNCE2', 'v3'): {'max_curvature': 0.09303786470255475, 'error': 0.0068183682878701305}, ('RigidNCE3', 'v3'): {'max_curvature': 0.07973297644861192, 'error': 0.0060031490273046465}, ('SemiRigidE1_50', 'v3'): {'max_curvature': 0.08554759975457213, 'error': 0.007391262093393689}, ('SemiRigidE2_50', 'v3'): {'max_curvature': 0.06798073221576455, 'error': 0.009048710685892464}, ('SemiRigidE3_50', 'v3'): {'max_curvature': 0.07182848216483961, 'error': 0.01109768205995741}}
		rv3 = {('SemiRigidE1', 'v3'): {'max_curvature': 0.14138229165211977, 'error': 0.009939512501818144}, ('SemiRigidE2', 'v3'): {'max_curvature': 0.09590345687591227, 'error': 0.012985067764987571}, ('SemiRigidE3', 'v3'): {'max_curvature': 0.09899645065259409, 'error': 0.012048776986263728}, ('RigidNCE1', 'v3'): {'max_curvature': 0.08306147103292033, 'error': 0.009246283692621288}, ('RigidNCE2', 'v3'): {'max_curvature': 0.09303786470255475, 'error': 0.0068183682878701305}, ('RigidNCE3', 'v3'): {'max_curvature': 0.07973297644861192, 'error': 0.0060031490273046465}, ('SemiRigidE1_50', 'v3'): {'max_curvature': 0.08554759975457213, 'error': 0.007391262093393689}, ('SemiRigidE2_50', 'v3'): {'max_curvature': 0.06798073221576455, 'error': 0.009048710685892464}, ('SemiRigidE3_50', 'v3'): {'max_curvature': 0.07182848216483961, 'error': 0.01109768205995741}}
		rv4 = {('SemiRigidE1', 'v4'): {'max_curvature': 0.05971304047033326, 'error': 0.008322770673742149}, ('SemiRigidE2', 'v4'): {'max_curvature': 0.054651103152930515, 'error': 0.013003120778768753}, ('SemiRigidE3', 'v4'): {'max_curvature': 0.05971304047033326, 'error': 0.008322770673742149}, ('RigidNCE1', 'v4'): {'max_curvature': 0.01353342905131086, 'error': 0.005733248904276125}, ('RigidNCE2', 'v4'): {'max_curvature': 0.013438574649747445, 'error': 0.003816312464565079}, ('RigidNCE3', 'v4'): {'max_curvature': 0.01955441827998415, 'error': 0.0037397354822525345}, ('SemiRigidE1_50', 'v4'): {'max_curvature': 0.02162732348804522, 'error': 0.008484450277008799}, ('SemiRigidE2_50', 'v4'): {'max_curvature': 0.028833551566472864, 'error': 0.00995067324451432}, ('SemiRigidE3_50', 'v4'): {'max_curvature': 0.02262899299946414, 'error': 0.01284419400338599}}
		rv4 = {('SemiRigidE1', 'v4'): {'max_curvature': 0.05971304047033326, 'error': 0.008322770673742149}, ('SemiRigidE2', 'v4'): {'max_curvature': 0.054651103152930515, 'error': 0.013003120778768753}, ('SemiRigidE3', 'v4'): {'max_curvature': 0.03697780564422796, 'error': 0.012660110100315187}, ('RigidNCE1', 'v4'): {'max_curvature': 0.01353342905131086, 'error': 0.005733248904276125}, ('RigidNCE2', 'v4'): {'max_curvature': 0.013438574649747445, 'error': 0.003816312464565079}, ('RigidNCE3', 'v4'): {'max_curvature': 0.01955441827998415, 'error': 0.0037397354822525345}, ('SemiRigidE1_50', 'v4'): {'max_curvature': 0.02162732348804522, 'error': 0.008484450277008799}, ('SemiRigidE2_50', 'v4'): {'max_curvature': 0.028833551566472864, 'error': 0.00995067324451432}, ('SemiRigidE3_50', 'v4'): {'max_curvature': 0.02262899299946414, 'error': 0.01284419400338599}}
		rv5 = {('SemiRigidE1', 'v5'): {'max_curvature': 0.08523669890680102, 'error': 0.006924914329886742}, ('SemiRigidE2', 'v5'): {'max_curvature': 0.06947206800599949, 'error': 0.008977696057825358}, ('SemiRigidE3', 'v5'): {'max_curvature': 0.08523669890680102, 'error': 0.006924914329886742}, ('RigidNCE1', 'v5'): {'max_curvature': 0.049633332753975454, 'error': 0.007487852076472254}, ('RigidNCE2', 'v5'): {'max_curvature': 0.030629816745371068, 'error': 0.004490451284959868}, ('RigidNCE3', 'v5'): {'max_curvature': 0.062337499205996634, 'error': 0.003970682237691373}, ('SemiRigidE1_50', 'v5'): {'max_curvature': 0.07504472543247931, 'error': 0.005085518807699702}, ('SemiRigidE2_50', 'v5'): {'max_curvature': 0.05366675331490715, 'error': 0.005376297046617336}, ('SemiRigidE3_50', 'v5'): {'max_curvature': 0.08480847096132728, 'error': 0.008098087497265188}}
		rv5 = {('SemiRigidE1', 'v5'): {'max_curvature': 0.08523669890680102, 'error': 0.006924914329886742}, ('SemiRigidE2', 'v5'): {'max_curvature': 0.06947206800599949, 'error': 0.008977696057825358}, ('SemiRigidE3', 'v5'): {'max_curvature': 0.08865340495632093, 'error': 0.006582880323167377}, ('RigidNCE1', 'v5'): {'max_curvature': 0.049633332753975454, 'error': 0.007487852076472254}, ('RigidNCE2', 'v5'): {'max_curvature': 0.030629816745371068, 'error': 0.004490451284959868}, ('RigidNCE3', 'v5'): {'max_curvature': 0.062337499205996634, 'error': 0.003970682237691373}, ('SemiRigidE1_50', 'v5'): {'max_curvature': 0.07504472543247931, 'error': 0.005085518807699702}, ('SemiRigidE2_50', 'v5'): {'max_curvature': 0.05366675331490715, 'error': 0.005376297046617336}, ('SemiRigidE3_50', 'v5'): {'max_curvature': 0.08480847096132728, 'error': 0.008098087497265188}}

		rv3 = {('SemiRigidE1', 'v3'): {'max_curvature': 0.14138229165211977, 'error': 0.009939512501818144}, ('SemiRigidE2', 'v3'): {'max_curvature': 0.09590345687591227, 'error': 0.012985067764987571}, ('SemiRigidE3', 'v3'): {'max_curvature': 0.09899645065259409, 'error': 0.012048776986263728}, ('RigidNCE1', 'v3'): {'max_curvature': 0.08306147103292033, 'error': 0.009246283692621288}, ('RigidNCE2', 'v3'): {'max_curvature': 0.09303786470255475, 'error': 0.0068183682878701305}, ('RigidNCE3', 'v3'): {'max_curvature': 0.07973297644861192, 'error': 0.0060031490273046465}, ('SemiRigidE1_50', 'v3'): {'max_curvature': 0.08554759975457213, 'error': 0.007391262093393689}, ('SemiRigidE2_50', 'v3'): {'max_curvature': 0.06798073221576455, 'error': 0.009048710685892464}, ('SemiRigidE3_50', 'v3'): {'max_curvature': 0.07182848216483961, 'error': 0.01109768205995741}, ('CL160ENS-1', 'v3'): {'max_curvature': 0.055212865915820225, 'error': 0.0049689136139662105}, ('CL160ENS-2', 'v3'): {'max_curvature': 0.06258987318989456, 'error': 0.0029506162443664495}, ('CL160ENS-3', 'v3'): {'max_curvature': 0.04762465625474257, 'error': 0.006305119921132189}}
		rv4 = {('SemiRigidE1', 'v4'): {'max_curvature': 0.05971304047033326, 'error': 0.008322770673742149}, ('SemiRigidE2', 'v4'): {'max_curvature': 0.054651103152930515, 'error': 0.013003120778768753}, ('SemiRigidE3', 'v4'): {'max_curvature': 0.03697780564422796, 'error': 0.012660110100315187}, ('RigidNCE1', 'v4'): {'max_curvature': 0.01353342905131086, 'error': 0.005733248904276125}, ('RigidNCE2', 'v4'): {'max_curvature': 0.013438574649747445, 'error': 0.003816312464565079}, ('RigidNCE3', 'v4'): {'max_curvature': 0.01955441827998415, 'error': 0.0037397354822525345}, ('SemiRigidE1_50', 'v4'): {'max_curvature': 0.02162732348804522, 'error': 0.008484450277008799}, ('SemiRigidE2_50', 'v4'): {'max_curvature': 0.028833551566472864, 'error': 0.00995067324451432}, ('SemiRigidE3_50', 'v4'): {'max_curvature': 0.02262899299946414, 'error': 0.01284419400338599}, ('CL160ENS-1', 'v4'): {'max_curvature': 0.012734906303023503, 'error': 0.003922456478143239}, ('CL160ENS-2', 'v4'): {'max_curvature': 0.011181604877787046, 'error': 0.002656720708137164}, ('CL160ENS-3', 'v4'): {'max_curvature': 0.014730482053864201, 'error': 0.007221131812825215}}
		rv5 = {('SemiRigidE1', 'v5'): {'max_curvature': 0.08523669890680102, 'error': 0.006924914329886742}, ('SemiRigidE2', 'v5'): {'max_curvature': 0.06947206800599949, 'error': 0.008977696057825358}, ('SemiRigidE3', 'v5'): {'max_curvature': 0.08865340495632093, 'error': 0.006582880323167377}, ('RigidNCE1', 'v5'): {'max_curvature': 0.049633332753975454, 'error': 0.007487852076472254}, ('RigidNCE2', 'v5'): {'max_curvature': 0.030629816745371068, 'error': 0.004490451284959868}, ('RigidNCE3', 'v5'): {'max_curvature': 0.062337499205996634, 'error': 0.003970682237691373}, ('SemiRigidE1_50', 'v5'): {'max_curvature': 0.07504472543247931, 'error': 0.005085518807699702}, ('SemiRigidE2_50', 'v5'): {'max_curvature': 0.05366675331490715, 'error': 0.005376297046617336}, ('SemiRigidE3_50', 'v5'): {'max_curvature': 0.08480847096132728, 'error': 0.008098087497265188}, ('CL160ENS-1', 'v5'): {'max_curvature': 0.03458151820333353, 'error': 0.003736936270495518}, ('CL160ENS-2', 'v5'): {'max_curvature': 0.03041611269641718, 'error': 0.002623271833668978}, ('CL160ENS-3', 'v5'): {'max_curvature': 0.03379066556325716, 'error': 0.005217957665003501}}
		report = rv3
		report.update(rv4)
		report.update(rv5)

		import re
		for key,val in report.items(): 
			report[key]['style'] = re.match('^(SemiRigid|Rigid|CL160ENS)',key[0]).group(1)
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
		titles = {'v3':'neighborhood\n(50,40x20)','v5':'neighborhood\n(100,40x20)','v4':'pixel\n(80x40)'}
		import pandas as pd
		import matplotlib as mpl
		import matplotlib.pyplot as plt
		from omni.plotter.panels import square_tiles
		from omni.base.store import picturesave
		#!!!!! WHY IS THIS SO HARD! to get the bar plots
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
					result[key] = [[report[(i,tag)][key] for i in g] for g in groups]
				df = pd.DataFrame(result[item])
				ax.set_title('%s %s'%(item,titles[tag]))
				#! df.groupby('styles').errors.value_counts().plot.bar()
				#! counter intuitiuve color ordering
				plot = df.plot.bar(color=color_list,ax=ax)
				ax.get_legend().remove()
				ax.set_xticklabels(['SemiRigid','Rigid','Semi-\nRigid (50)','Dextran'],rotation=0)
				ax.set_ylim((0,maxes[item]))
		picturesave('fig.summary',work.plotdir,backup=False,version=True,meta={})

