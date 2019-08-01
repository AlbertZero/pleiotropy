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

fitness_table = numpy.zeros((8,8),dtype='float')
fitness_table_vplus = numpy.zeros((8,8),dtype='float')
#print(fit_dict)
for i in range(8):
	for j in range(8):
		
		evolution_environment = evolution_environments[i]
		measurement_environment = measurement_environments[j]
		
		gen_scalefactor = gen_dict[evolution_environment]/700.
		
		###Evolution env on rows, measurement env on columns; in general fitness_table[i,j] != fitness_table[j,i]
		
		fitlist = []
		fitlist_vplus = []
		for clone in fit_dict[evolution_environment]:
			fitlist.append(fit_dict[evolution_environment][clone][measurement_environment]*100./gen_scalefactor) ###fitness gained/lost per 700 gens of evolution, in percent
			if fit_dict[evolution_environment][clone]["virus"] == "Virus retained":
				fitlist_vplus.append(fit_dict[evolution_environment][clone][measurement_environment]*100./gen_scalefactor)
				
		fitness_table[i,j] = numpy.median(fitlist) ###Median fitness across the clones from evolution_environment i, in measurement environment j
		fitness_table_vplus[i,j] = numpy.median(fitlist_vplus)

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
			
pt.figure()

pt.pcolor( fitness_table[::-1], cmap = 'RdBu', vmin=-8,vmax=8,edgecolors='w', linewidth=.5)
pt.colorbar(label='Median fitness (%)')
ax = pt.gca()
ax.set_xticks(numpy.arange(.5,8.7,1))
ax.set_yticks(numpy.arange(.5,8.7,1))
ax.set_xticklabels(evolution_env_labels,fontsize=8)
ax.set_yticklabels(measurement_env_labels[::-1],fontsize=8)
ax.set_xlim(0,8)
ax.set_ylim(0,8)
ax.tick_params(axis=u'both', which=u'both',length=0)
ax.set_xlabel('Measurement environment',fontsize=14)
ax.set_ylabel('Evolution environment',fontsize=14)
pt.plot([0,8],[8,0],'--',color='Gray',linewidth=.5)
for i in range(8):
	for j in range(8):
		if fitness_table[i,7-j] < 0 and fitness_table[i,7-j] > -.1: ###This will get printed as -0.0, so print 0.0 instead
			ax.text(7-j+.5,7-i+.5,abs(numpy.round(fitness_table[i,7-j],1)),color='Gray',fontsize=7,horizontalalignment='center',verticalalignment='center')
		elif i==7-j:
			ax.text(7-j+.5,7-i+.5,numpy.round(fitness_table[i,7-j],1),color='k',fontsize=7,horizontalalignment='center',verticalalignment='center')

		else:
			ax.text(7-j+.5,7-i+.5,numpy.round(fitness_table[i,7-j],1),color='Gray',fontsize=7,horizontalalignment='center',verticalalignment='center')
		
pt.savefig('figures/median_fitness_heatmap_check.pdf',bbox_inches='tight')

##Plot the table with clones that lost the virus removed

pt.figure()

pt.pcolor( fitness_table_vplus[::-1], cmap = 'RdBu', vmin=-5,vmax=5,edgecolors='w' )
pt.colorbar(label='Median fitness (%)')
ax = pt.gca()
ax.set_xticks(numpy.arange(.5,8.7,1))
ax.set_yticks(numpy.arange(.5,8.7,1))
ax.set_xticklabels(evolution_env_labels,fontsize=8)
ax.set_yticklabels(measurement_env_labels[::-1],fontsize=8)
ax.set_xlim(0,8)
ax.set_ylim(0,8)
ax.tick_params(axis=u'both', which=u'both',length=0)
ax.set_xlabel('Measurement environment',fontsize=14)
ax.set_ylabel('Evolution environment',fontsize=14)
pt.plot([0,8],[8,0],'--',color='Gray',linewidth=.5)
for i in range(8):
	for j in range(8):
		if fitness_table_vplus[i,7-j] < 0 and fitness_table_vplus[i,7-j] > -.1: ###This will get printed as -0.0, so print 0.0 instead
			ax.text(7-j+.5,7-i+.5,abs(numpy.round(fitness_table_vplus[i,7-j],1)),color='Gray',fontsize=7,horizontalalignment='center',verticalalignment='center')
		elif i==7-j:
			ax.text(7-j+.5,7-i+.5,numpy.round(fitness_table_vplus[i,7-j],1),color='k',fontsize=7,horizontalalignment='center',verticalalignment='center')

		else:
			ax.text(7-j+.5,7-i+.5,numpy.round(fitness_table_vplus[i,7-j],1),color='Gray',fontsize=7,horizontalalignment='center',verticalalignment='center')
pt.savefig('figures/median_fitness_heatmap_vplusonly_check.pdf',bbox_inches='tight')