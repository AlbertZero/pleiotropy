import matplotlib.pylab as pt
import numpy
import matplotlib.gridspec as gridspec
import matplotlib
import scipy.stats
from matplotlib.colors import LinearSegmentedColormap

import scipy.cluster.hierarchy as cluster
import scipy.spatial.distance as distance
import matplotlib.cm as cm

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import random

import import_utilities

def assign_color(xy_loc):
	###This function assigns a color to an (x,y) coordinate location
	norm_instance = matplotlib.colors.Normalize(vmin=0,vmax=1)
	colormap_obj1 = cm.ScalarMappable(cmap="viridis",norm=norm_instance)
	c1 = colormap_obj1.to_rgba(xy_loc[0])[0:3]
	grayscale_factor = xy_loc[1]
	
	c2 = grayscale_factor*numpy.array(c1) + (1. - grayscale_factor)*numpy.array([.5,.5,.5])
	
	color=c2
	return color
	
##Import fitness data

fit_dict, gen_dict, std_err_dict, random_reps = import_utilities.import_fitness_data('data/pleiotropy_fitness_data_7_4_2019.txt')


####Get phenotypic clusters across set of diagnostic measurement environments

evolution_environments = ['Evolved in SC','Evolved in .8 M','Evolved in 37 C', 'Evolved in .07% glu', 'Evolved in gal','Evolved in 21 C','Evolved in pH 3','Evolved in pH 7.3']
envt_labels = ['SC','High salt','High temp','Low glu','Gal','Low temp','pH 3','pH 7.3']
measurement_environments = ['fitness in SC','fitness in gal','fitness in 21 C','fitness in 37 C','fitness in .07% glu','fitness in pH 3','fitness in pH 7.3','fitness in .8 M']

##Get all fitness vectors over the measurement environments

fit_dict_extr = {}
fit_mat = []
v_list = []
hd_list = []
env_list = []
virus_dict = {}
hapdip_dict = {}
err_dict_extr = {}

fit_mat_vplus = []
hd_list_vplus = []
env_list_vplus = []

for evol_env in evolution_environments:

	gen_norm_factor = gen_dict[evol_env]/700.
	fit_dict_extr[evol_env] = []
	err_dict_extr[evol_env] = []
	virus_dict[evol_env] = []
	hapdip_dict[evol_env] = []
	
	for clone in fit_dict[evol_env]:
		
		fit_vec = [ fit_dict[evol_env][clone][m]/gen_norm_factor for m in measurement_environments ]
		err_vec = [ std_err_dict[evol_env][clone][m]/gen_norm_factor for m in measurement_environments ]
		
		fit_dict_extr[evol_env].append(fit_vec)
		err_dict_extr[evol_env].append(err_vec)
		
		fit_mat.append(fit_vec)
		v_list.append(fit_dict[evol_env][clone]["virus"])
		
		
		hd_list.append(fit_dict[evol_env][clone]["hapdip"])
		env_list.append(evol_env)
		virus_dict[evol_env].append(fit_dict[evol_env][clone]["virus"])
		hapdip_dict[evol_env].append(fit_dict[evol_env][clone]["hapdip"])
		
		if fit_dict[evol_env][clone]["virus"] == "Virus retained":
			fit_mat_vplus.append(fit_vec)
			hd_list_vplus.append(fit_dict[evol_env][clone]["hapdip"])
			env_list_vplus.append(evol_env)
			
	fit_dict_extr[evol_env] = numpy.array( fit_dict_extr[evol_env] )
	err_dict_extr[evol_env] = numpy.array(err_dict_extr[evol_env] )
	
fit_mat = numpy.array(fit_mat) ##All the clone fitness vectors concatenated
fit_mat_vplus = numpy.array(fit_mat_vplus)


	
####Run a tsne

tsne_obj = TSNE(n_components = 2, random_state=11) ###random state determined via minimizing over kl divergence in 20 random seeds
fit_mat_tsne = tsne_obj.fit_transform(fit_mat)

kl = tsne_obj.kl_divergence_

# for rs in numpy.arange(19)+1:
# 	tsne_obj_new = TSNE(n_components = 2, random_state=rs)
# 	fit_mat_tr = tsne_obj_new.fit_transform(fit_mat)
# 	kl_new = tsne_obj_new.kl_divergence_
# 	if kl_new < kl:
# 		fit_mat_tsne = fit_mat_tr
# 		kl = kl_new
# 		print(rs)

fit_mat_tr_scaled = fit_mat_tsne.copy()
fit_mat_tr_scaled[:,0] = fit_mat_tr_scaled[:,0] - numpy.min(fit_mat_tsne[:,0])
fit_mat_tr_scaled[:,1] = fit_mat_tr_scaled[:,1] - numpy.min(fit_mat_tsne[:,1])
fit_mat_tr_scaled[:,0] = fit_mat_tr_scaled[:,0]/numpy.max(fit_mat_tr_scaled[:,0])
fit_mat_tr_scaled[:,1] = fit_mat_tr_scaled[:,1]/numpy.max(fit_mat_tr_scaled[:,1])

colors = ['MidnightBlue','silver','plum']
colormap_obj = cm.ScalarMappable(cmap="Paired")
envt_colors = dict(zip(evolution_environments,colormap_obj.to_rgba(numpy.arange(8))))
#print(envt_colors)

color_mat = {}
color_mat_envs = {}
tsne_divided = {}

c_list_env = []

for i in range(len(v_list)):
	tsne_loc = fit_mat_tr_scaled[i,:]
	env = env_list[i]
	env_color = envt_colors[env]
	if v_list[i] == "Virus lost" and hd_list[i] != "Diploid":
		
		if "vminus" not in color_mat:
			color_mat["vminus"] = []
			color_mat_envs["vminus"] = []
			tsne_divided["vminus"] = []
		color_mat["vminus"].append( assign_color(tsne_loc) )
		color_mat_envs["vminus"].append(env_color)
		tsne_divided["vminus"].append( fit_mat_tr_scaled[i,:] )
	
	elif v_list[i] == "Virus lost":
		if "vminus-dip" not in color_mat:
			color_mat["vminus-dip"] = []
			color_mat_envs["vminus-dip"] = []
			tsne_divided["vminus-dip"] = []
		color_mat["vminus-dip"].append( assign_color(tsne_loc) )
		color_mat_envs["vminus-dip"].append(env_color)
		tsne_divided["vminus-dip"].append( fit_mat_tr_scaled[i,:] )
		
	elif hd_list[i] == "Diploid":
		
		if "vplus-dip" not in color_mat:
			color_mat["vplus-dip"] = []
			color_mat_envs["vplus-dip"] = []
			tsne_divided["vplus-dip"] = []
		color_mat["vplus-dip"].append( assign_color(tsne_loc) )
		color_mat_envs["vplus-dip"].append(env_color)
		tsne_divided["vplus-dip"].append( fit_mat_tr_scaled[i,:] )
		
	else:
		if "vplus" not in color_mat:
			color_mat["vplus"] = []
			color_mat_envs["vplus"] = []
			tsne_divided["vplus"] = []
		color_mat["vplus"].append( assign_color(tsne_loc) )
		color_mat_envs["vplus"].append(env_color)
		tsne_divided["vplus"].append( fit_mat_tr_scaled[i,:] )
		
		
	c_list_env.append( envt_colors[env_list[i]][0:3])
	
dummy_fig = pt.figure()
ax = pt.gca()

key_list = ["vplus","vplus-dip","vminus","vminus-dip"]
grp_labels = ["V+/Haploid","V+/Diploid","V-/Haploid","V-/Diploid"]
marker_list = ["o","v","s","^"]

legend_handles2 = []
for i in range(len(evolution_environments)):
	env = evolution_environments[i]
	color = envt_colors[env]
	leg, = ax.plot(random.random(),random.random(), 'o', c=color)
	legend_handles2.append(leg)

legend_handles1 = []
for i in range(len(key_list)):
	grp = grp_labels[i]
	grp_marker = marker_list[i]
	leg, = ax.plot(random.random(),random.random(), marker=grp_marker, c='k', linestyle='None')
	legend_handles1.append(leg)

pt.savefig('figures/dummy_legend_fig.pdf',bbox_inches='tight')

fig1 = pt.figure(figsize=(14,13))

gs1 = gridspec.GridSpec(1, 2, bottom=.55, wspace = .2)

ax1 = pt.subplot(gs1[0,0])
ax2 = pt.subplot(gs1[0,1])

#fig, (ax1,ax2) = pt.subplots(1,2,figsize=(14,6))

for i in range(len(key_list)):
	key = key_list[i]
	tsne_divided[key] = numpy.array(tsne_divided[key])
	grp_marker = grp_labels[i]
	npts = tsne_divided[key].shape[0]
	ax1.scatter(tsne_divided[key][:,0],tsne_divided[key][:,1], c=color_mat_envs[key], marker=marker_list[i], alpha=.9)

ax1.legend(legend_handles2,envt_labels)
ax1.set_xlabel('tsne1')
ax1.set_ylabel('tsne2')

for i in range(len(key_list)):
	key = key_list[i]
	tsne_divided[key] = numpy.array(tsne_divided[key])
	grp_marker = grp_labels[i]
	npts = tsne_divided[key].shape[0]
	ax2.scatter(tsne_divided[key][:,0],tsne_divided[key][:,1], c=color_mat[key], marker=marker_list[i], alpha=.9)

	#for pt in range(npts):
		#alpha_ts = tsne_divided[key][pt,1]
		#ax2.scatter(tsne_divided[key][pt,0],tsne_divided[key][pt,1], c=color_mat[key][pt], marker=marker_list[i], alpha=alpha_ts)
	
legend1 = pt.legend(legend_handles1[2:],grp_labels[2:],loc=(.7,.5))
legend2 = pt.legend(legend_handles1[0:2],grp_labels[:2],loc=(.5,.88))
ax2.add_artist(legend1)
ax2.add_artist(legend2)
#ax2.legend(legend_handles1[3:],grp_labels[3:],loc=(.8,.6))

ax2.set_xlabel('tsne1')
ax2.set_ylabel('tsne2')

ax1.text(-.12,1.06,'A',fontname='Arial',fontsize=18)
ax2.text(-.12,1.06,'B',fontname='Arial',fontsize=18)
####Scale the fit_mat_tr so the range is from 0 to 1 in each axis

tsne_dict_extr = {}
for evol_env in evolution_environments:
	tsne_dict_extr[evol_env] = fit_mat_tr_scaled[ numpy.where(numpy.array(env_list) == evol_env)[0] ]

evolution_envs_name_dict = dict(zip(evolution_environments, ['Evolved in SC','Evolved in High Salt','Evolved in High Temp','Evolved in Low Glu','Evolved in Gal','Evolved in Low Temp','Evolved in pH 3','Evolved in pH 7.3','Evolved in High Eth']))

####Plot in each envt with colors following post-facto V, D, M labeling

gs2 = gridspec.GridSpec(2, 4, wspace=0,hspace=0,top=.5)

ax1 = pt.subplot(gs2[0,0])
ax2 = pt.subplot(gs2[0,1])
ax3 = pt.subplot(gs2[0,2])
ax4 = pt.subplot(gs2[0,3])


ax5 = pt.subplot(gs2[1,0])
ax6 = pt.subplot(gs2[1,1])
ax7 = pt.subplot(gs2[1,2])
ax8 = pt.subplot(gs2[1,3])

#fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = pt.subplots(3,3, figsize = (14,12))

ax_handles = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]

colors = ['MidnightBlue','silver','plum']
env_counter = 0

###This plotting loop is broken up into two pieces so that the V+ clones are plotted in order of their tsne coordinate, so that the layering of the lines and the colors accord

for evol_env in evolution_environments:
	color_order = []
	
	for i in range(len(fit_dict_extr[evol_env][:,0])):
		vir_status = virus_dict[evol_env][i]
		hapdip_status = hapdip_dict[evol_env][i]
		tsne_loc = tsne_dict_extr[evol_env][i]
		
		color_tsne = assign_color(tsne_loc)
		color_order.append(color_tsne[0])
		if vir_status == "Virus lost" and hapdip_status != "Diploid":
			marker_grp = 'o'
			ax_handles[env_counter].errorbar( numpy.arange(8), fit_dict_extr[evol_env][i,:]*100, yerr = err_dict_extr[evol_env][i,:]*100, color = color_tsne, marker = marker_grp, markeredgewidth=0, capsize=0, alpha = 1 )
		elif vir_status == "Virus lost":
			marker_grp = '^'
			ax_handles[env_counter].errorbar( numpy.arange(8), fit_dict_extr[evol_env][i,:]*100, yerr = err_dict_extr[evol_env][i,:]*100, color = color_tsne, marker = marker_grp, markeredgewidth=0, capsize=0, alpha = 1 )
		elif hapdip_status == "Diploid":
			marker_grp = 'v'
		else:
			marker_grp = 's'
	
	color_ordering = numpy.argsort(color_order)
	
	for i in color_ordering:
		vir_status = virus_dict[evol_env][i]
		hapdip_status = hapdip_dict[evol_env][i]
		tsne_loc = tsne_dict_extr[evol_env][i]
		
		color_tsne = assign_color(tsne_loc)
		color_order.append(color_tsne[0])
		if vir_status == "Virus lost" and hapdip_status != "Diploid":
			marker_grp = 'v'
		elif vir_status == "Virus lost":
			marker_grp = '^'
		
		elif hapdip_status == "Diploid":
			marker_grp = 's'
			ax_handles[env_counter].errorbar( numpy.arange(8), fit_dict_extr[evol_env][i,:]*100, yerr = err_dict_extr[evol_env][i,:]*100, color = color_tsne, marker = marker_grp, markeredgewidth=0, capsize=0, alpha = 1 )

		else:
			marker_grp = 'o'
			ax_handles[env_counter].errorbar( numpy.arange(8), fit_dict_extr[evol_env][i,:]*100, yerr = err_dict_extr[evol_env][i,:]*100, color = color_tsne, marker = marker_grp, markeredgewidth=0, capsize=0, alpha = 1 )

		
				
		
	ax_handles[env_counter].set_xlim([-.5,7.5])
	ax_handles[env_counter].set_xticks(numpy.arange(8))
	ax_handles[env_counter].grid(color='grey',alpha=.3)
	if env_counter > 3.5:
		ax_handles[env_counter].set_xticklabels( ['SC','Gal','Low\nTemp','High\nTemp','Low\nGlu','pH\n3','pH\n7.3','High\nSalt'],fontsize=7)
	else:
		ax_handles[env_counter].set_xticklabels([])
	ax_handles[env_counter].set_ylim([-28,20])
	ax_handles[env_counter].set_yticks([-20,-10,0,10])
	if env_counter == 0 or env_counter == 4:
		ax_handles[env_counter].set_ylabel('Relative fitness (%)')
	else:
		ax_handles[env_counter].set_yticklabels([])
		
	label = evolution_envs_name_dict[evol_env]
	
	ax_handles[env_counter].text( 3.5,15, label, horizontalalignment='center',verticalalignment='center' )
	
	env_counter += 1
ax1.text(-1.2,22,'C',fontname='Arial',fontsize=18)	
pt.savefig('figures/tsne_VDM_squiggle_check.pdf', bbox_inches='tight')