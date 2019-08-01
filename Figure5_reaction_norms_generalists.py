import matplotlib.pylab as pt
import numpy
import matplotlib.gridspec as gridspec
import matplotlib
import scipy.stats
from matplotlib.colors import LinearSegmentedColormap
import import_utilities

fit_dict, gen_dict, std_err_dict, random_reps = import_utilities.import_fitness_data('data/pleiotropy_fitness_data_7_4_2019.txt')

fig = pt.figure(figsize=(10,10))
gs = gridspec.GridSpec(2,3,bottom=.37)
gs.update(wspace=.2,hspace=.2)

axa = pt.subplot(gs[0,0])
axb = pt.subplot(gs[0,1])
axc = pt.subplot(gs[0,2])

ax1 = pt.subplot(gs[1,0])
ax2 = pt.subplot(gs[1,1])
ax3 = pt.subplot(gs[1,2])



#pt.figure(figsize=(12,4))
gs2 = gridspec.GridSpec(1,3,top=.3)
gs2.update(wspace=.2)

ax4 = pt.subplot(gs2[0,0])
ax5 = pt.subplot(gs2[0,1])
ax6 = pt.subplot(gs2[0,2])



pad = .05 ##Fraction of x range to put as white space at the start and end of each x axis

###Salt gradient

coolwarm = pt.get_cmap('coolwarm')

measurement_envs = ['fitness in SC','fitness in .2 M','fitness in .4 M','fitness in .8 M']
evolution_envs = ['Evolved in SC','Evolved in .2 M','Evolved in .8 M']
color_list = ['peru','saddlebrown','brown']#,'MediumSlateBlue','Tomato','black','brown','gray','seagreen','darkmagenta']
axis_label_list = ['SC','Med salt','Med high salt','High salt']
legend_label_list = ['Evolved in SC (0 M)','Evolved in .2 M salt','Evolved in .8 M salt']

x_scale_list = numpy.array([0,.2,.4,.8])


nlevels = 4
mean_fit_vec = numpy.zeros( (4,), dtype='float' )

legend_handles_list = []
legend_handles_list2 = []

nacl_fit_dict = {}
nacl_err_dict = {}

env_ind = 0
for evol_env in evolution_envs:
	
	gen_time_norm_factor = gen_dict[evol_env]/700.
	clone_count = 0
	
	mean_fit_vec[:] = 0
	nacl_fit_dict[evol_env] = []
	nacl_err_dict[evol_env] = []
	for clone in fit_dict[evol_env]:
		
		if (fit_dict[evol_env][clone]["virus"] == "Virus retained" and fit_dict[evol_env][clone]["hapdip"] != "Diploid"): #Including only v+ clones and presumptive haploids (note that we are including clones where the staining was indeterminate as haploids)
			
			fit_vec = numpy.array([ fit_dict[evol_env][clone][ measurement_envs[0] ]/gen_time_norm_factor, fit_dict[evol_env][clone][measurement_envs[1]]/gen_time_norm_factor, fit_dict[evol_env][clone][measurement_envs[2]]/gen_time_norm_factor, fit_dict[evol_env][clone][measurement_envs[3]]/gen_time_norm_factor ])*100
			std_err_vec = numpy.array([ std_err_dict[evol_env][clone][ measurement_envs[0] ]/gen_time_norm_factor, std_err_dict[evol_env][clone][measurement_envs[1]]/gen_time_norm_factor, std_err_dict[evol_env][clone][measurement_envs[2]]/gen_time_norm_factor, std_err_dict[evol_env][clone][measurement_envs[3]]/gen_time_norm_factor ])*100
			
			mean_fit_vec += fit_vec
			clone_count += 1.
			
			
			nacl_fit_dict[evol_env].append(fit_vec)
			nacl_err_dict[evol_env].append(std_err_vec)
			
	
	nacl_fit_dict[evol_env] = numpy.array(nacl_fit_dict[evol_env])
	nacl_err_dict[evol_env] = numpy.array(nacl_err_dict[evol_env])
	n = float(len(nacl_fit_dict[evol_env][:,0]))
	
	frac_beneficial = numpy.array([ numpy.sum( numpy.greater(nacl_fit_dict[evol_env][:,i], 0) )/n for i in range(nlevels) ])
	
	err = numpy.sqrt(n*frac_beneficial*(1-frac_beneficial))/n
	
	legend_handles2 = ax1.errorbar( x_scale_list+.005*env_ind, frac_beneficial, yerr=err, color = color_list[env_ind], capsize = 0, marker='o', markeredgewidth = 0, linewidth = 3, elinewidth = 2)
	legend_handles_list2.append(legend_handles2)
	
	#if env_ind < 1.5:
	#ax1.plot(x_scale_list[env_ind]+.005*env_ind, frac_beneficial[env_ind], color = color_list[env_ind], marker='o', markersize=10, markeredgecolor='k')
	#else:
	#	ax1.plot(x_scale_list[env_ind+1]+.005*env_ind, frac_beneficial[env_ind+1], color = color_list[env_ind], marker='o', markersize=10, markeredgecolor='k')
	#####
	
	generalists = numpy.sum( [numpy.all(nacl_fit_dict[evol_env][i,:] > 0) for i in range(int(n))] )/n
	print(generalists)
	gen_err = numpy.sqrt(n*generalists*(1-generalists))/n
	
	legend_handle = ax4.bar( env_ind,  generalists, yerr=gen_err, width=.48,color='C0' )
	legend_handles_list.append(legend_handle)
	
	legend_handle3 = axa.errorbar(x_scale_list+.005*env_ind, numpy.mean(nacl_fit_dict[evol_env],axis=0), yerr=numpy.std(nacl_fit_dict[evol_env],axis=0)/numpy.sqrt(clone_count - 1), marker='o', color = color_list[env_ind], linewidth = 3, elinewidth = 2)
	
	env_ind += 1
	
# env_ind = 0
# for evol_env in evolution_envs:
# 	
# 	
# 	if env_ind < 1.5:
# 		axa.plot(x_scale_list[env_ind]+.005*env_ind, numpy.mean(nacl_fit_dict[evol_env],axis=0)[env_ind], color = color_list[env_ind], marker='X', markersize=14, markeredgecolor='k')
# 	else:
# 		axa.plot(x_scale_list[env_ind+1]+.005*env_ind, numpy.mean(nacl_fit_dict[evol_env],axis=0)[env_ind+1], color = color_list[env_ind], marker='X', markersize=14, markeredgecolor='k')
# 	env_ind += 1
	
ax4.set_xticks(numpy.arange(3))
ax4.set_xticklabels( ['0','.2', '.8'], fontsize=12 )

ax1.set_xlim([-.8*pad, .8 + .8*pad])
ax1.set_xticks( x_scale_list)
ax1.set_xticklabels( ['0','.2', '.4', '.8'], fontsize=12 )


ax4.set_ylabel( 'Fraction of generalists' )
ax4.set_xlabel( 'Salt concentration, home env.' )
ax1.set_xlabel( 'Salt concentration (M)' )
ax1.set_ylabel('Fraction beneficial')
ax4.set_ylim([0,1])
#ax4.legend(legend_handles_list,legend_label_list,loc='upper left',fontsize=8)
	
axa.legend( legend_handles_list2, legend_label_list, numpoints = 1, loc='upper left', fontsize=8)

axa.set_xlim([-.8*pad, .8 + .8*pad])

axa.set_xticks( x_scale_list)
axa.set_xticklabels([])

axa.set_ylabel( 'Mean relative fitness (%)' )

####pH (low to high)

measurement_envs = ['fitness in pH 3','fitness in pH 3.8','fitness in SC','fitness in pH 6','fitness in pH 7.3']
evolution_envs = ['Evolved in pH 3','Evolved in pH 3.8','Evolved in SC','Evolved in pH 6','Evolved in pH 7.3']
color_list = ['dodgerblue','mediumblue','midnightblue','Tomato','brown']#,'blueviolet','seagreen','gray','Tomato','brown']
legend_label_list = ['Evolved at pH 3','Evolved at pH 3.8','Evolved in SC (4.5)','Evolved at pH 6','Evolved at pH 7.3']
x_scale_list = numpy.array([3, 3.8, 4.5, 6, 7.3])

mean_fit_vec = numpy.zeros( (5,), dtype='float' )
ph_err_dict = {}
ph_fit_dict = {}
legend_handles_list = []
legend_handles_list2 = []
nlevels = 5
env_ind = 0
for evol_env in evolution_envs:
	
	gen_time_norm_factor = gen_dict[evol_env]/700.
	clone_count = 0
	
	mean_fit_vec[:] = 0
	ph_fit_dict[evol_env] = []
	ph_err_dict[evol_env] = []
	for clone in fit_dict[evol_env]:
		
		if (fit_dict[evol_env][clone]["virus"] == "Virus retained" and fit_dict[evol_env][clone]["hapdip"] != "Diploid"): #Including only v+ clones
			
			fit_vec = numpy.array([ fit_dict[evol_env][clone][ measurement_envs[0] ]/gen_time_norm_factor, fit_dict[evol_env][clone][measurement_envs[1]]/gen_time_norm_factor, fit_dict[evol_env][clone][measurement_envs[2]]/gen_time_norm_factor, fit_dict[evol_env][clone][measurement_envs[3]]/gen_time_norm_factor, fit_dict[evol_env][clone][measurement_envs[4]]/gen_time_norm_factor ])*100
			std_err_vec = numpy.array([ std_err_dict[evol_env][clone][ measurement_envs[0] ]/gen_time_norm_factor, std_err_dict[evol_env][clone][measurement_envs[1]]/gen_time_norm_factor, std_err_dict[evol_env][clone][measurement_envs[2]]/gen_time_norm_factor, std_err_dict[evol_env][clone][measurement_envs[3]]/gen_time_norm_factor, std_err_dict[evol_env][clone][measurement_envs[4]]/gen_time_norm_factor ])*100
			
			mean_fit_vec += fit_vec
			clone_count += 1.
			ph_fit_dict[evol_env].append(fit_vec)
			ph_err_dict[evol_env].append(std_err_vec)
			
			
	
	ph_fit_dict[evol_env] = numpy.array(ph_fit_dict[evol_env])
	ph_err_dict[evol_env] = numpy.array(ph_err_dict[evol_env])
	n = float(len(ph_fit_dict[evol_env][:,0]))
	frac_beneficial = numpy.array([ numpy.sum( numpy.greater(ph_fit_dict[evol_env][:,i], 0) )/n for i in range(nlevels) ])
	
	err = numpy.sqrt(n*frac_beneficial*(1-frac_beneficial))/n
	
	legend_handles2 = ax2.errorbar( x_scale_list+.02*env_ind, frac_beneficial,  color = color_list[env_ind], yerr=err, capsize = 0, marker='o', markeredgewidth = 0, linewidth = 3, elinewidth = 2)
	legend_handles_list2.append(legend_handles2)
	
	#ax2.plot(x_scale_list[env_ind]+.02*env_ind, frac_beneficial[env_ind], color = color_list[env_ind], marker='o', markersize=8)

	#####
	generalists = numpy.sum( [numpy.all(ph_fit_dict[evol_env][i,:] > 0) for i in range(int(n))] )/n
	
	gen_err = numpy.sqrt(n*generalists*(1-generalists))/n
	
	legend_handle = ax5.bar( env_ind,  generalists, yerr=gen_err, width=.8,color='C0' )
	legend_handles_list.append(legend_handle)
	legend_handle = axb.errorbar(x_scale_list+.02*env_ind, numpy.mean(ph_fit_dict[evol_env],axis=0), yerr=numpy.std(ph_fit_dict[evol_env],axis=0)/numpy.sqrt(clone_count - 1), marker='o', color = color_list[env_ind], linewidth = 3, elinewidth = 2)
	#axb.plot(x_scale_list[env_ind]+.02*env_ind, numpy.mean(ph_fit_dict[evol_env],axis=0)[env_ind], color = color_list[env_ind], marker='o', markersize=8)

	env_ind += 1


ax5.set_xticks(numpy.arange(5))
ax5.set_xticklabels( ['3','3.8', '4.5', '6', '7.3'], fontsize=12 )
ax5.set_ylim([0,1])
#ax5.legend(legend_handles_list,legend_label_list,loc='upper left',fontsize=8)
ax2.set_xlim([3 - 4.3*pad, 7.3 + 4.3*pad])
ax2.set_xticks( 	x_scale_list)
ax2.set_xticklabels( ['3','3.8', '4.5', '6', '7.3'], fontsize=12 )


ax2.set_xlabel( 'pH' )
ax5.set_xlabel( 'pH, home env.' )


axb.legend( legend_handles_list2, legend_label_list, numpoints = 1, loc='upper center', fontsize=8)
axb.set_ylim([-6,8])
axb.set_xlim([3 - 4.3*pad, 7.3 + 4.3*pad])
axb.set_xticks(x_scale_list)
axb.set_xticklabels([])
###Temperature (low to high)

measurement_envs = ['fitness in 21 C','fitness in SC','fitness in 34 C','fitness in 37 C']
evolution_envs = ['Evolved in 21 C','Evolved in SC','Evolved in 37 C']
legend_label_list = ['Evolved at 21 C','Evolved in SC (30 C)','Evolved at 37 C']
color_list = ['midnightblue','peru','Brown']

x_scale_list = numpy.array([21,30,34,37])

nlevels = 4

mean_fit_vec = numpy.zeros( (4,), dtype='float' )
temp_fit_dict = {}
temp_err_dict = {}
legend_handles_list = []
legend_handles_list2 = []
env_ind = 0
for evol_env in evolution_envs:
	
	gen_time_norm_factor = gen_dict[evol_env]/700.
	clone_count = 0
	
	mean_fit_vec[:] = 0
	temp_fit_dict[evol_env] = []
	temp_err_dict[evol_env] = []
	for clone in fit_dict[evol_env]:
		
		if (fit_dict[evol_env][clone]["virus"] == "Virus retained" and fit_dict[evol_env][clone]["hapdip"] != "Diploid"): #Including only v+ clones
			
			fit_vec = numpy.array([ fit_dict[evol_env][clone][ measurement_envs[0] ]/gen_time_norm_factor, fit_dict[evol_env][clone][measurement_envs[1]]/gen_time_norm_factor, fit_dict[evol_env][clone][measurement_envs[2]]/gen_time_norm_factor, fit_dict[evol_env][clone][measurement_envs[3]]/gen_time_norm_factor ])*100
			std_err_vec = numpy.array([ std_err_dict[evol_env][clone][ measurement_envs[0] ]/gen_time_norm_factor, std_err_dict[evol_env][clone][measurement_envs[1]]/gen_time_norm_factor, std_err_dict[evol_env][clone][measurement_envs[2]]/gen_time_norm_factor, std_err_dict[evol_env][clone][measurement_envs[3]]/gen_time_norm_factor ])*100
			temp_fit_dict[evol_env].append(fit_vec)
			mean_fit_vec += fit_vec
			clone_count += 1.
			temp_err_dict[evol_env].append(std_err_vec)
			
	
	temp_fit_dict[evol_env] = numpy.array(temp_fit_dict[evol_env])
	temp_err_dict[evol_env] = numpy.array(temp_err_dict[evol_env])
	n = float(len(temp_fit_dict[evol_env][:,0]))
	
	frac_beneficial = numpy.array([ numpy.sum(numpy.greater( temp_fit_dict[evol_env][:,i], 0) )/n for i in range(nlevels) ])
	if evol_env == 'Evolved in SC':
		print(temp_fit_dict[evol_env])
	err = numpy.sqrt(n*frac_beneficial*(1-frac_beneficial))/n
	
	legend_handles2 = ax3.errorbar( x_scale_list+.1*env_ind, frac_beneficial, yerr=err, color = color_list[env_ind], capsize = 0, markeredgewidth = 0, marker='o', linewidth = 3, elinewidth = 2)
	legend_handles_list2.append(legend_handles2)
	#if env_ind < 1.5:
	#	ax3.plot(x_scale_list[env_ind]+.1*env_ind, frac_beneficial[env_ind], color = color_list[env_ind], marker='o', markersize=8)
	#else:
	#	ax3.plot(x_scale_list[env_ind+1]+.1*env_ind, frac_beneficial[env_ind+1], color = color_list[env_ind], marker='o', markersize=8)

	#####
	generalists = numpy.sum( [numpy.all(temp_fit_dict[evol_env][i,:] > 0) for i in range(int(n))] )/n
	
	gen_err = numpy.sqrt(n*generalists*(1-generalists))/n
	
	legend_handle = ax6.bar( env_ind,  generalists, yerr=gen_err, width=.48,color='C0' )
	legend_handles_list.append(legend_handle)
	legend_handle = axc.errorbar(x_scale_list+.1*env_ind, numpy.mean(temp_fit_dict[evol_env],axis=0), yerr=numpy.std(temp_fit_dict[evol_env],axis=0)/numpy.sqrt(clone_count - 1), marker='o', color = color_list[env_ind], linewidth = 3, elinewidth = 2)
	#if env_ind < 1.5:
	#	axc.plot(x_scale_list[env_ind]+.1*env_ind, numpy.mean(temp_fit_dict[evol_env],axis=0)[env_ind], color = color_list[env_ind], marker='o', markersize=8)
	#else:
	#	axc.plot(x_scale_list[env_ind+1]+.1*env_ind, numpy.mean(temp_fit_dict[evol_env],axis=0)[env_ind+1], color = color_list[env_ind], marker='o', markersize=8)

	env_ind += 1


ax1.set_ylim(-.05,1.05)
ax2.set_ylim(-.05,1.05)
ax3.set_ylim(-.05,1.05)

ax3.set_xlim([ 21 - 16*pad, 37 + 16*pad])
ax3.set_xticks( x_scale_list)
ax3.set_xticklabels( ['21','30', '34', '37'], fontsize=12 )
ax6.set_xticks(numpy.arange(3))
ax6.set_xticklabels( ['21','30', '37'], fontsize=12 )
#ax6.legend(legend_handles_list,legend_label_list,loc='upper left',fontsize=8)

ax6.set_xlabel( 'Temperature, home env.' )
ax3.set_xlabel( 'Temperature (degrees C)' )

axc.legend( legend_handles_list2, legend_label_list, numpoints = 1, loc='upper left', fontsize=8)

axc.set_xlim([ 21 - 16*pad, 37 + 16*pad])
axc.set_xticks( x_scale_list)
axc.set_xticklabels([])

ax_labels = ['A','B','C','D','E','F','G','H','I']
ax_handles = [axa,axb,axc,ax1,ax2,ax3,ax4,ax5,ax6]
for i in range(len(ax_handles)):
	ax = ax_handles[i]
	ax.text(-.1,1.1,ax_labels[i],horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, fontsize=14, fontname='Arial')

pt.savefig('figures/reactionnorms_fracbeneficial_generalists_check.pdf',bbox_inches='tight')
pt.close()
