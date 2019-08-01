import matplotlib.pylab as pt
import numpy
import matplotlib.gridspec as gridspec
import matplotlib
import scipy.stats
from matplotlib.colors import LinearSegmentedColormap
import import_utilities

##Import fitness data

fit_dict, gen_dict, std_err_dict, random_reps = import_utilities.import_fitness_data('data/pleiotropy_fitness_data_7_4_2019.txt')

##Set up giant matrix figure

measurement_envs = ['fitness in SC','fitness in 21 C','fitness in .07% glu','fitness in 37 C','fitness in gal','fitness in pH 3','fitness in pH 7.3','fitness in .8 M']

evolution_envs = ['Evolved in SC','Evolved in 37 C', 'Evolved in .8 M','Evolved in .07% glu','Evolved in 21 C', 'Evolved in gal', 'Evolved in pH 3','Evolved in pH 7.3','Evolved in high eth']
color_list = ['Orange','MediumSlateBlue','midnightblue','Tomato','black','brown','gray','seagreen','darkmagenta']

##Set up background color map

cdict1 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.1),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 1.0),
                   (0.5, 0.1, 0.0),
                   (1.0, 0.0, 0.0))
        }
cdict1['alpha'] = ((0, .5, .5), (.5, .5, .5), (1, .5, .5))
cmap1 = LinearSegmentedColormap('cmap1',cdict1)


coolwarm = pt.get_cmap('coolwarm')

measurement_envs = ['fitness in SC','fitness in .2 M','fitness in .4 M','fitness in .8 M']
evolution_envs = ['Evolved in SC','Evolved in .2 M','Evolved in .8 M']
color_list = ['peru','saddlebrown','brown']#,'MediumSlateBlue','Tomato','black','brown','gray','seagreen','darkmagenta']
axis_label_list = ['SC','.2 M salt','.4 M salt','.8 M salt']
legend_label_list = ['Evolved in SC (0 M)','Evolved in .2 M salt','Evolved in .8 M salt']

dummy_figure = pt.figure()
ax = pt.gca()
sc_handle = ax.errorbar(0, 0, xerr = .1, yerr =.1, fmt='o', color= color_list[0], markeredgewidth=0, capsize=0)
med_salt_handle = ax.errorbar(0, 0, xerr = .1, yerr =.1, fmt='o', color= color_list[1], markeredgewidth=0, capsize=0)
high_salt_handle = ax.errorbar(0, 0, xerr = .1, yerr =.1, fmt='o', color= color_list[2], markeredgewidth=0, capsize=0)

fig = pt.figure(figsize=(6,6))
gs_object = gridspec.GridSpec(4,4)
#gs_object.update(wspace=.15,hspace=.15)
gs_object.update(wspace=.05,hspace=.05)

for i in numpy.arange(4):
	for	j in range(i+1):
		
		meas_env1 = measurement_envs[3-j]
		meas_env2 = measurement_envs[3-i]
		
		ax = pt.subplot(gs_object[i,j])
		ax2 = pt.subplot(gs_object[j,i])
		if j < i:
			xvals = []
			yvals = []
			for k in range(3):
				
				evol_env = evolution_envs[k]
				gen_time_norm_factor = gen_dict[evol_env]/700.
				
				for clone in fit_dict[evol_env]:
					if fit_dict[evol_env][clone]["virus"] == "Virus retained": #Including only v+ clones
						
						ax.errorbar(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor, fit_dict[evol_env][clone][meas_env2]/gen_time_norm_factor, xerr = std_err_dict[evol_env][clone][meas_env1], yerr = std_err_dict[evol_env][clone][meas_env2], fmt='o', color= color_list[k], markeredgewidth=0, capsize=0, markersize=4.8)
						xvals.append(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor)
						yvals.append(fit_dict[evol_env][clone][meas_env2]/gen_time_norm_factor)
			
			corr_coeff = (scipy.stats.pearsonr(numpy.array(xvals),numpy.array(yvals))[0])
			color_scale = (corr_coeff + 1)/2.
			ax2.patch.set_facecolor(color=cmap1(color_scale))
			ax2.text(.5,.5,str(numpy.round(corr_coeff,1)),horizontalalignment='center',verticalalignment='center',fontsize=14)
			ax2.set_xticks([])
			ax2.set_yticks([])
			
			y_lb = numpy.percentile(yvals, 0) - .02
			y_ub = numpy.percentile(yvals, 100) + .02
			x_lb = numpy.percentile(xvals, 0) - .02
			x_ub = numpy.percentile(xvals, 100) + .02
			ax.set_xlim([x_lb,x_ub])
			ax.set_ylim([y_lb,y_ub])
			if j < .5:
				ax.set_yticks([.8*y_lb,.8*y_ub])
				ax.set_yticklabels([round(100*.8*y_lb,1),round(100*.8*y_ub,1)],fontsize=9)
				#ax.set_yticks([])
				ax.set_ylabel(axis_label_list[3-i],fontsize=10)
				if i == 2:
					ax.text(-.3,.23,'Relative Fitness (%)',rotation=90,fontsize=14)
			else:
				ax.set_yticks([])
			if i > 2.5:
				ax.set_xticks([.8*x_lb,.8*x_ub])
				ax.set_xticklabels([round(100*.8*x_lb,1),round(100*.8*x_ub,1)],fontsize=9)
				#ax.set_xticks([])
				ax.set_xlabel(axis_label_list[3-j],fontsize=10)
				if j == 1:
					ax.text(-.105,-.22,'Relative Fitness (%)',fontsize=14)
			else:
				ax.set_xticks([])
			ax.axhline(0,0,1,color='k')
			ax.axvline(0,0,1,color='k')
		else:
			xvals = []
			hist_list = []
			for k in range(3):
				evol_env = evolution_envs[k]
				gen_time_norm_factor = gen_dict[evol_env]/700.
				hist_list.append([])
				for clone in fit_dict[evol_env]:
					if fit_dict[evol_env][clone]["virus"] == "Virus retained": #Including only v+ clones
					
						ax.plot(random_reps[evol_env][clone][meas_env1][0]/gen_time_norm_factor, random_reps[evol_env][clone][meas_env2][1]/gen_time_norm_factor, 'o', color= color_list[k], markeredgewidth=0, markersize=4.8)
						xvals.append(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor)
						hist_list[k].append(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor)
			lb = numpy.percentile(xvals, 0) - .02
			ub = numpy.percentile(xvals, 100) + .02
			ax.set_xlim([lb,ub])
			ax.set_ylim([lb,ub])
		
			if j < .5:
				ax.set_yticks([.8*lb,.8*ub])
				ax.set_yticklabels([round(100*.8*lb,1),round(100*.8*ub,1)],fontsize=9)
				#ax.set_yticks([])
				ax.set_ylabel(axis_label_list[3-i],fontsize=10)
			else:
				ax.set_yticks([])
			if i > 2.5:
				ax.set_xticks([.8*lb,.8*ub])
				ax.set_xticklabels([round(100*.8*lb,1),round(100*.8*ub,1)],fontsize=9)
				#ax.set_xticks([])
				ax.set_xlabel(axis_label_list[3-j],fontsize=10)
			else:
				ax.set_xticks([])
			ax.axhline(0,0,1,color='k')
			ax.axvline(0,0,1,color='k')
			
			
pt.legend([sc_handle,med_salt_handle,high_salt_handle], legend_label_list,  bbox_to_anchor=(-1.78,4.7),numpoints=1, fontsize=8)

pt.savefig('figures/Salt_gradient_figure_check.pdf',bbox_inches='tight')

measurement_envs = ['fitness in pH 3','fitness in pH 3.8','fitness in SC','fitness in pH 6','fitness in pH 7.3']
evolution_envs = ['Evolved in pH 3','Evolved in pH 3.8','Evolved in SC','Evolved in pH 6','Evolved in pH 7.3']
color_list = ['Orange','MediumSlateBlue','midnightblue','black','darkgray']#,'blueviolet','seagreen','gray','Tomato','brown']


pt.figure()
ax = pt.gca()
color_list = ['dodgerblue','mediumblue','midnightblue','Tomato','brown']#,'blueviolet','seagreen','gray','Tomato','brown']
sc_handle = ax.errorbar(0, 0, xerr = .1, yerr =.1, fmt='o', color= color_list[2], markeredgewidth=0, capsize=0)
pH3_handle = ax.errorbar(0, 0, xerr = .1, yerr =.1, fmt='o', color= color_list[0], markeredgewidth=0, capsize=0)
pH38_handle = ax.errorbar(0, 0, xerr = .1, yerr =.1, fmt='o', color= color_list[1], markeredgewidth=0, capsize=0)
pH6_handle = ax.errorbar(0, 0, xerr = .1, yerr =.1, fmt='o', color= color_list[3], markeredgewidth=0, capsize=0)
pH7_handle = ax.errorbar(0, 0, xerr = .1, yerr =.1, fmt='o', color= color_list[4], markeredgewidth=0, capsize=0)
axis_label_list = ['pH 3','pH 3.8','SC (pH 4.5)','pH 6','pH 7.3']
legend_label_list = ['Evolved at pH 3','Evolved at pH 3.8','Evolved in SC','Evolved at pH 6','Evolved at pH 7.3']
fig = pt.figure(figsize=(7,7))
gs_object = gridspec.GridSpec(5,5)
gs_object.update(wspace=.05,hspace=.05)
for i in range(5):
	for j in range(i+1):
		
		meas_env1 = measurement_envs[j]
		meas_env2 = measurement_envs[i]
		
		ax = pt.subplot(gs_object[i,j])
		ax2 = pt.subplot(gs_object[j,i])
		if j < i:
			xvals = []
			yvals = []
			for k in range(5):
				evol_env = evolution_envs[k]
				gen_time_norm_factor = gen_dict[evol_env]/700.
				
				for clone in fit_dict[evol_env]:
					if fit_dict[evol_env][clone]["virus"] == "Virus retained": #Including only v+ clones
					#print meas_env1, meas_env2, evol_env
						ax.errorbar(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor, fit_dict[evol_env][clone][meas_env2]/gen_time_norm_factor, xerr = std_err_dict[evol_env][clone][meas_env1], yerr = std_err_dict[evol_env][clone][meas_env2], fmt='o', color= color_list[k], markeredgewidth=0, capsize=0, markersize=4.8)
						xvals.append(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor)
						yvals.append(fit_dict[evol_env][clone][meas_env2]/gen_time_norm_factor)
						
			y_lb = numpy.percentile(yvals, 0) - .02
			y_ub = numpy.percentile(yvals, 100) + .02
			x_lb = numpy.percentile(xvals, 0) - .02
			x_ub = numpy.percentile(xvals, 100) + .02
			
			ax.set_xlim([x_lb,x_ub])
			ax.set_ylim([y_lb,y_ub])
			
			if j < .5:
				ax.set_yticks([.8*y_lb,.8*y_ub])
				ax.set_yticklabels([round(100*.8*y_lb,1),round(100*.8*y_ub,1)],fontsize=9)
				#ax.set_yticks([])
				ax.set_ylabel(axis_label_list[i],fontsize=10)
				if i == 2:
					ax.text(-.28,.1,'Relative Fitness (%)',rotation=90,fontsize=14)
			else:
				ax.set_yticks([])
			if i > 3.5:
				ax.set_xticks([.8*x_lb,.8*x_ub])
				ax.set_xticklabels([round(100*.8*x_lb,1),round(100*.8*x_ub,1)],fontsize=9)
				#ax.set_xticks([])
				ax.set_xlabel(axis_label_list[j],fontsize=10)
				#ax.set_xticklabels([])
				if j == 1:
					ax.text(-.08,-.24,'Relative Fitness (%)',fontsize=14)
			else:
				ax.set_xticks([])
			ax.axhline(0,0,1,color='k')
			ax.axvline(0,0,1,color='k')
			corr_coeff = scipy.stats.pearsonr(numpy.array(xvals),numpy.array(yvals))[0]
			color_scale = (corr_coeff + 1)/2.
			ax2.patch.set_facecolor(color=cmap1(color_scale))
			ax2.text(.5,.5,str(numpy.round(corr_coeff,1)),horizontalalignment='center',verticalalignment='center',fontsize=14)
			ax2.set_xticks([])
			ax2.set_yticks([])
		else:
			xvals = []
			hist_list = []
			for k in range(5):
				evol_env = evolution_envs[k]
				gen_time_norm_factor = gen_dict[evol_env]/700.
				hist_list.append([])
				for clone in fit_dict[evol_env]:
					if fit_dict[evol_env][clone]["virus"] == "Virus retained": #Including only v+ clones
					
						ax.plot(random_reps[evol_env][clone][meas_env1][0]/gen_time_norm_factor, random_reps[evol_env][clone][meas_env2][1]/gen_time_norm_factor, 'o', color= color_list[k], markeredgewidth=0, markersize=4.8)
						xvals.append(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor)
						hist_list[k].append(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor)
			lb = numpy.percentile(xvals, 0) - .02
			ub = numpy.percentile(xvals, 100) + .02
			ax.set_xlim([lb,ub])
			ax.set_ylim([lb,ub])
			
			if j < .5:
				ax.set_yticks([.8*lb,.8*ub])
				#ax.set_yticks([])
				ax.set_yticklabels([round(100*.8*lb,1),round(100*.8*ub,1)],fontsize=9)
				#ax.set_yticklabels([])
				ax.set_ylabel(axis_label_list[i],fontsize=10)
			else:
				ax.set_yticks([])
			if i > 3.5:
				ax.set_xticks([.8*lb,.8*ub])
				#ax.set_xticks([])
				ax.set_xlabel(axis_label_list[j],fontsize=10)
				ax.set_xticklabels([round(100*.8*lb,1),round(100*.8*ub,1)],fontsize=9)
				
			else:
				ax.set_xticks([])
			
			
				
			ax.axhline(0,0,1,color='k')
			ax.axvline(0,0,1,color='k')
			
pt.legend([pH3_handle, pH38_handle, sc_handle,pH6_handle,pH7_handle], legend_label_list,  bbox_to_anchor=(-4,5.2),numpoints=1,ncol=2,fontsize=10,loc = 3)

pt.savefig('figures/pH_gradient_check.pdf',bbox_inches='tight')

measurement_envs = ['fitness in 21 C','fitness in SC','fitness in 34 C','fitness in 37 C']
evolution_envs = ['Evolved in 21 C','Evolved in SC','Evolved in 37 C']
color_list = ['midnightblue','peru','Brown']
axis_label_list = ['21 C','30 C','34 C','37 C']
legend_label_list = ['Evolved at 21 C','Evolved at 30 C','Evolved at 37 C']

pt.figure()
ax = pt.gca()
highT_handle = ax.errorbar(0, 0, xerr = .1, yerr =.1, fmt='o', color= color_list[2], markeredgewidth=0, capsize=0)
lowT_handle = ax.errorbar(0, 0, xerr = .1, yerr =.1, fmt='o', color= color_list[0], markeredgewidth=0, capsize=0)
SC_handle = ax.errorbar(0, 0, xerr = .1, yerr =.1, fmt='o', color= color_list[1], markeredgewidth=0, capsize=0)




fig = pt.figure(figsize=(7,7))
gs_object = gridspec.GridSpec(5,5)
gs_object.update(wspace=.05,hspace=.05)
for i in range(4):
	for j in range(i+1):
		
		meas_env1 = measurement_envs[j]
		meas_env2 = measurement_envs[i]
		
		ax = pt.subplot(gs_object[i,j])
		ax2 = pt.subplot(gs_object[j,i])
		if j < i:
			xvals = []
			yvals = []
			for k in range(3):
				evol_env = evolution_envs[k]
				gen_time_norm_factor = gen_dict[evol_env]/700.
				
				for clone in fit_dict[evol_env]:
					if fit_dict[evol_env][clone]["virus"] == "Virus retained": #Including only v+ clones
					#print meas_env1, meas_env2, evol_env
						ax.errorbar(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor, fit_dict[evol_env][clone][meas_env2]/gen_time_norm_factor, xerr = std_err_dict[evol_env][clone][meas_env1], yerr = std_err_dict[evol_env][clone][meas_env2], fmt='o', color= color_list[k], markeredgewidth=0, capsize=0, markersize=4.8)
						xvals.append(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor)
						yvals.append(fit_dict[evol_env][clone][meas_env2]/gen_time_norm_factor)
			y_lb = numpy.percentile(yvals, 0) - .02
			y_ub = numpy.percentile(yvals, 100) + .02
			x_lb = numpy.percentile(xvals, 0) - .02
			x_ub = numpy.percentile(xvals, 100) + .02
			
			ax.set_xlim([x_lb,x_ub])
			ax.set_ylim([y_lb,y_ub])
			
			if j < .5:
				ax.set_yticks([.8*y_lb,.8*y_ub])
				ax.set_yticklabels([round(100*.8*y_lb,1),round(100*.8*y_ub,1)],fontsize=9)
				#ax.set_yticks([])
				ax.set_ylabel(axis_label_list[i],fontsize=10)
				if i == 2:
					ax.text(-.28,.135,'Relative Fitness (%)',rotation=90,fontsize=14)
			else:
				ax.set_yticks([])
			if i > 2.5:
				ax.set_xticks([.8*x_lb,.8*x_ub])
				ax.set_xticklabels([round(100*.8*x_lb,1),round(100*.8*x_ub,1)],fontsize=9)
				#ax.set_xticks([])
				ax.set_xlabel(axis_label_list[j],fontsize=10)
				if j == 1:
					ax.text(-.08,-.23,'Relative Fitness (%)',fontsize=14)
			else:
				ax.set_xticks([])
			ax.axhline(0,0,1,color='k')
			ax.axvline(0,0,1,color='k')
			corr_coeff = scipy.stats.pearsonr(numpy.array(xvals),numpy.array(yvals))[0]
			color_scale = (corr_coeff + 1)/2
			ax2.patch.set_facecolor(color=cmap1(color_scale))
			ax2.text(.5,.5,str(numpy.round(corr_coeff,1)),horizontalalignment='center',verticalalignment='center',fontsize=14)
			ax2.set_xticks([])
			ax2.set_yticks([])
		else:
			xvals = []
			hist_list = []
			for k in range(3):
				evol_env = evolution_envs[k]
				gen_time_norm_factor = gen_dict[evol_env]/700.
				hist_list.append([])
				for clone in fit_dict[evol_env]:
					if fit_dict[evol_env][clone]["virus"] == "Virus retained": #Including only v+ clones
					
						ax.plot(random_reps[evol_env][clone][meas_env1][0]/gen_time_norm_factor, random_reps[evol_env][clone][meas_env2][1]/gen_time_norm_factor, 'o', color= color_list[k], markeredgewidth=0, markersize=4.8)
						xvals.append(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor)
						hist_list[k].append(fit_dict[evol_env][clone][meas_env1]/gen_time_norm_factor)
			lb = numpy.percentile(xvals, 0) - .02
			ub = numpy.percentile(xvals, 100) + .02
			ax.set_xlim([lb,ub])
			ax.set_ylim([lb,ub])
			
			if j < .5:
				ax.set_yticks([.8*lb,.8*ub])
				ax.set_yticklabels([round(100*lb,1),round(100*ub,1)],fontsize=9)
				#ax.set_yticks([])
				ax.set_ylabel(axis_label_list[i],fontsize=10)
			else:
				ax.set_yticks([])
			if i > 2.5:
				ax.set_xticks([.8*lb,.8*ub])
				ax.set_xticklabels([round(100*.8*lb,1),round(100*.8*ub,1)],fontsize=9)
				#ax.set_xticks([])
				ax.set_xlabel(axis_label_list[j],fontsize=10)
			else:
				ax.set_xticks([])
			ax.axhline(0,0,1,color='k')
			ax.axvline(0,0,1,color='k')
			
			
#gs_object.tight_layout(fig)
pt.legend([lowT_handle, SC_handle, highT_handle], legend_label_list,  bbox_to_anchor=(-3.22,4.2),numpoints=1,ncol=1,fontsize=10,loc = 3)

pt.savefig('figures/Temp_gradient_check.pdf',bbox_inches='tight')