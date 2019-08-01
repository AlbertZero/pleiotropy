import numpy
import matplotlib.pylab as pt
import import_utilities

fit_dict_cured_ref_fa = import_utilities.import_cured_ref_FA('data/Cured_ref_21C_fitness_table.txt')
fit_dict, gen_dict, std_err_dict, random_reps = import_utilities.import_fitness_data('data/pleiotropy_fitness_data_7_4_2019.txt')


extr_envs = ['Evolved in SC','Evolved in .8 M','Evolved in 37 C', 'Evolved in .07% glu', 'Evolved in gal','Evolved in 21 C','Evolved in pH 3','Evolved in pH 7.3']


pt.figure()
for evol_env in fit_dict_cured_ref_fa:
	if evol_env in extr_envs:
		for clone in fit_dict_cured_ref_fa[evol_env]:
		
			if fit_dict[evol_env][clone]["virus"] == "Virus lost":
			
				handle1,=pt.plot( [.1*numpy.random.random(), .5 + .1*numpy.random.random()], [(fit_dict_cured_ref_fa[evol_env][clone]['virusplusref'])*100,(fit_dict_cured_ref_fa[evol_env][clone]['virusminusref'])*100], linestyle='-', marker = 'o', markersize=4, markeredgewidth = 0, color = 'Brown', alpha = .6, linewidth=.7)
			
			else:
				handle2,=pt.plot( [.1*numpy.random.random(), .5 + .1*numpy.random.random()], [(fit_dict_cured_ref_fa[evol_env][clone]['virusplusref'])*100,(fit_dict_cured_ref_fa[evol_env][clone]['virusminusref'])*100], linestyle='-', marker = 'o', markersize=4, markeredgewidth = 0, color = 'MidnightBlue', alpha = .6, linewidth=.7)

pt.legend([handle2,handle1],['V+','V-'])				 
pt.xlim([-.05,.65])
ax = pt.gca()
ax.set_xticks([.05,.55])
ax.set_xticklabels(['Ancestor with virus', 'Cured ancestor'],fontsize=14)
ax.set_ylabel('Fitness at 21 C (%)',fontsize=14)
pt.savefig('figures/Fitness_v_cured_ref_21C_check.pdf',bbox_inches='tight')