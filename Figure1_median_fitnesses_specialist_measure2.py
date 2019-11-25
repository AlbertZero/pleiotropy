import matplotlib.pylab as pt
import numpy
import matplotlib.gridspec as gridspec
import matplotlib
import scipy.stats
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm
import import_utilities

fit_dict, gen_dict, std_err_dict, random_reps = import_utilities.import_fitness_data('data/pleiotropy_fitness_data_7_4_2019.txt')

####Get fitnesses across set of diagnostic measurement environments
evolution_environments = ['Evolved in SC','Evolved in .07% glu', 'Evolved in gal','Evolved in 37 C' ,'Evolved in 21 C','Evolved in .8 M','Evolved in pH 3','Evolved in pH 7.3']

evolution_env_labels = ['SC','Low\nglu','Gal','High\ntemp','Low\ntemp','High\nsalt','pH 3','pH 7.3']
measurement_environments = ['fitness in SC','fitness in .07% glu','fitness in gal','fitness in 37 C','fitness in 21 C','fitness in .8 M','fitness in pH 3','fitness in pH 7.3']
measurement_env_labels = ['SC','Low glu','Gal','High temp','Low temp','High salt','pH 3','pH 7.3']

##Get a table of the median fitness in each measurement and evolution condition

nbootstrap = 500

fitness_table = numpy.zeros((8,8),dtype='float')
fitness_table_vplus = numpy.zeros((8,8),dtype='float')
row_stat = []
row_stat_err = []
#print(fit_dict)
for i in range(8):
	fit_arr_temp = []
	err_arr_temp = []
	specialist_stat_bs = []
	for j in range(8):
		
		evolution_environment = evolution_environments[i]
		measurement_environment = measurement_environments[j]
		
		gen_scalefactor = gen_dict[evolution_environment]/700.
		
		###Evolution env on rows, measurement env on columns; in general fitness_table[i,j] != fitness_table[j,i]
		
		fitlist = []
		errlist = []
		fitlist_vplus = []
		for clone in fit_dict[evolution_environment]:
			fitlist.append(fit_dict[evolution_environment][clone][measurement_environment]*100./gen_scalefactor) ###fitness gained/lost per 700 gens of evolution, in percent
			errlist.append(std_err_dict[evolution_environment][clone][measurement_environment]*100./gen_scalefactor)
			if fit_dict[evolution_environment][clone]["virus"] == "Virus retained":
				fitlist_vplus.append(fit_dict[evolution_environment][clone][measurement_environment]*100./gen_scalefactor)
		fit_arr_temp.append(fitlist)
		err_arr_temp.append(errlist)
				
		fitness_table[i,j] = numpy.median(fitlist) ###Median fitness across the clones from evolution_environment i, in measurement environment j
		fitness_table_vplus[i,j] = numpy.median(fitlist_vplus)
		
	fit_arr_temp = numpy.array(fit_arr_temp).T
	err_arr_temp = numpy.array(err_arr_temp).T
	nclones,nenv = fit_arr_temp.shape
	if i==0:
		print(fit_arr_temp)
	#print(nclones,nenv)
	row_stat1 = numpy.mean(numpy.sum(fit_arr_temp < -2*err_arr_temp,axis=1)/nenv)
	for n in range(nbootstrap):
		clone_sel = numpy.random.choice(numpy.arange(nclones),size=nclones)
		fits_bs = fit_arr_temp[clone_sel,:]
		errs_bs = err_arr_temp[clone_sel,:]
		#if n == 1:
		#	print(numpy.sum(fit_arr_temp > 0,axis=1)/nenv)
		stat_bs = numpy.mean(numpy.sum(fits_bs < -2*errs_bs,axis=1)/nenv) ##Avg number of envts where fitness was lost
		specialist_stat_bs.append(stat_bs)
	#print(specialist_stat_bs)
	row_stat.append(row_stat1)
	#row_stat.append(numpy.percentile(specialist_stat_bs,50))
	row_stat_err.append([row_stat1-numpy.percentile(specialist_stat_bs,5),numpy.percentile(specialist_stat_bs,95)-row_stat1])

#####
col_stat = []
col_stat_err = []
for i in range(8):
	fit_arr_temp = []
	competition_stat_bs = []
	for j in range(8):
		
		evolution_environment = evolution_environments[j]
		measurement_environment = measurement_environments[i]
		
		gen_scalefactor = gen_dict[evolution_environment]/700.
		
		###Evolution env on rows, measurement env on columns; in general fitness_table[i,j] != fitness_table[j,i]
		
		fitlist = []
		fitlist_vplus = []
		for clone in fit_dict[evolution_environment]:
			fitlist.append(fit_dict[evolution_environment][clone][measurement_environment]*100./gen_scalefactor) ###fitness gained/lost per 700 gens of evolution, in percent
			if fit_dict[evolution_environment][clone]["virus"] == "Virus retained":
				fitlist_vplus.append(fit_dict[evolution_environment][clone][measurement_environment]*100./gen_scalefactor)
		fit_arr_temp.append(fitlist)
	
	other_envs = []
	evol_env = []
	
	for j in range(8):
		nclones = len(fit_arr_temp[j])
		clone_sel = fit_arr_temp[j]
		if i != j:
			other_envs.extend(clone_sel)
		else:
			evol_env.extend(clone_sel)
	stat = 0
	other_envs = numpy.array(other_envs)
	evol_env = numpy.array(evol_env)
	
	for clone_fit in evol_env:
		stat += numpy.sum(other_envs < clone_fit)/len(other_envs)
	stat = stat/len(evol_env)
					
	for n in range(nbootstrap):
		evol_env_bs = []
		other_envs_bs = []
		for j in range(8):
			nclones = len(fit_arr_temp[j])
			clone_sel = numpy.random.choice(fit_arr_temp[j],size=nclones)
			if i != j:
				other_envs_bs.extend(clone_sel)
			else:
				evol_env_bs.extend(clone_sel)
		evol_env_bs = numpy.array(evol_env_bs)
		other_envs_bs = numpy.array(other_envs_bs)
		bs_stat = 0
		for clone_fit in evol_env_bs:
			bs_stat += numpy.sum(other_envs_bs < clone_fit)/len(other_envs_bs)
		bs_stat = bs_stat/len(evol_env_bs)
		
		
		competition_stat_bs.append(bs_stat)
	col_stat.append(stat)
	#col_stat.append(numpy.percentile(competition_stat_bs,50))
	col_stat_err.append([stat-numpy.percentile(competition_stat_bs,5),numpy.percentile(competition_stat_bs,95)-stat])

###Construct axes to plot 

fig = pt.figure(figsize=(6,5.8))

ax = pt.axes([.16,.16,.74,.74])
cb_ax = pt.axes([.94,.16,.04,.74])
ax_v = pt.axes([0,.16,.14,.74])
ax_h = pt.axes([.16,0.,.74,.14])

row_stat_err = numpy.array(row_stat_err).T
col_stat_err = numpy.array(col_stat_err).T

ax_v.barh(numpy.arange(8),numpy.array(row_stat)[::-1],.7,xerr=row_stat_err[::-1],color='grey')
ax_h.bar(numpy.arange(8),col_stat,.7,yerr=col_stat_err,color='grey')

ax_v.set_ylim(-.5,7.5)
ax_h.set_xlim(-.5,7.5)



ax_v.set_xlim(0,1)
ax_h.set_ylim(0,1)
ax_v.xaxis.tick_top()
ax_v.set_xticks([0,.5,1])
ax_v.set_xlabel('Degree of\nSpecialization',fontsize=8)
ax_v.xaxis.set_label_position('top') 
ax_h.set_ylabel('Competitiveness\nin Home Env.',fontsize=8)

ax_v.set_yticks(numpy.arange(8))
ax_h.set_xticks(numpy.arange(8))

ax_h.set_xticklabels(evolution_env_labels,fontsize=10)
ax_v.set_yticklabels(measurement_env_labels[::-1],fontsize=10)

ax_h.set_xlabel('Measurement environment',fontsize=14)
ax_v.set_ylabel('Evolution environment',fontsize=14)

##Plot the table with a diverging colormap so we can see fitness gains/losses

edgecolor_constructor = []
linewidth_constructor = []
counter = 0
for i in range(8):
	for j in range(8):
		if i == 7-j or i == 8-j:
			edgecolor_constructor.append('w')
			linewidth_constructor.append(2)
		else:
			edgecolor_constructor.append('w')
			linewidth_constructor.append(.5)
			

im = ax.pcolor( fitness_table[::-1], cmap = 'RdBu', vmin=-8,vmax=8,edgecolors='w', linewidth=.5)
pt.colorbar(im, label='Median fitness (%)',cax=cb_ax)
#ax = pt.gca()
#ax.set_xticks(numpy.arange(.5,8.7,1))
#ax.set_yticks(numpy.arange(.5,8.7,1))
ax.set_xticklabels([])
ax.set_yticklabels([])

ax.set_xlim(0,8)
ax.set_ylim(0,8)
ax.tick_params(axis=u'both', which=u'both',length=0)

ax.plot([0,8],[8,0],'--',color='Gray',linewidth=.5)
for i in range(8):
	for j in range(8):
		if fitness_table[i,7-j] < 0 and fitness_table[i,7-j] > -.1: ###This will get printed as -0.0, so print 0.0 instead
			ax.text(7-j+.5,7-i+.5,abs(numpy.round(fitness_table[i,7-j],1)),color='Gray',fontsize=9,horizontalalignment='center',verticalalignment='center')
		elif i==7-j:
			ax.text(7-j+.5,7-i+.5,numpy.round(fitness_table[i,7-j],1),color='k',fontsize=9,horizontalalignment='center',verticalalignment='center')

		else:
			ax.text(7-j+.5,7-i+.5,numpy.round(fitness_table[i,7-j],1),color='Gray',fontsize=9,horizontalalignment='center',verticalalignment='center')
		
pt.savefig('figures/median_fitness_heatmap_with_specialist_stats_updated.pdf',bbox_inches='tight')

##Plot the table with clones that lost the virus removed

# pt.figure()
# 
# pt.pcolor( fitness_table_vplus[::-1], cmap = 'RdBu', vmin=-5,vmax=5,edgecolors='w' )
# pt.colorbar(label='Median fitness (%)')
# ax = pt.gca()
# ax.set_xticks(numpy.arange(.5,8.7,1))
# ax.set_yticks(numpy.arange(.5,8.7,1))
# ax.set_xticklabels(evolution_env_labels,fontsize=8)
# ax.set_yticklabels(measurement_env_labels[::-1],fontsize=8)
# ax.set_xlim(0,8)
# ax.set_ylim(0,8)
# ax.tick_params(axis=u'both', which=u'both',length=0)
# ax.set_xlabel('Measurement environment',fontsize=14)
# ax.set_ylabel('Evolution environment',fontsize=14)
# pt.plot([0,8],[8,0],'--',color='Gray',linewidth=.5)
# for i in range(8):
# 	for j in range(8):
# 		if fitness_table_vplus[i,7-j] < 0 and fitness_table_vplus[i,7-j] > -.1: ###This will get printed as -0.0, so print 0.0 instead
# 			ax.text(7-j+.5,7-i+.5,abs(numpy.round(fitness_table_vplus[i,7-j],1)),color='Gray',fontsize=7,horizontalalignment='center',verticalalignment='center')
# 		elif i==7-j:
# 			ax.text(7-j+.5,7-i+.5,numpy.round(fitness_table_vplus[i,7-j],1),color='k',fontsize=7,horizontalalignment='center',verticalalignment='center')
# 
# 		else:
# 			ax.text(7-j+.5,7-i+.5,numpy.round(fitness_table_vplus[i,7-j],1),color='Gray',fontsize=7,horizontalalignment='center',verticalalignment='center')
# pt.savefig('figures/median_fitness_heatmap_vplusonly_check.pdf',bbox_inches='tight')