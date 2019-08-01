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

fit_dict, gen_dict, std_err_dict, random_reps = import_utilities.import_fitness_data('data/pleiotropy_fitness_data_7_4_2019.txt')


####Get phenotypic clusters across set of diagnostic measurement environments

evolution_environments = ['Evolved in SC','Evolved in .8 M','Evolved in 37 C', 'Evolved in .07% glu', 'Evolved in gal','Evolved in 21 C','Evolved in pH 3','Evolved in pH 7.3','Evolved in pH 6','Evolved in pH 3.8','Evolved in .2 M']
envt_labels = ['SC','.8 M','37 C', '.07% glu', 'gal','21 C','pH 3','pH 7.3']
measurement_environments = ['fitness in SC','fitness in gal','fitness in 21 C','fitness in 37 C','fitness in .07% glu','fitness in pH 3','fitness in pH 7.3','fitness in .8 M']

####Count the number of V+/H,V+/U,V+/D,V-/H,V-/U,V-/D from each env

vdm_count_mat = []
clones_per_env = []
env_list = []

for env in evolution_environments:
	env_list.append(env)
	clones_per_env.append(0)
	vdm_count_mat.append([0,0,0,0,0,0])
	
	for clone in fit_dict[env]:
		clones_per_env[-1] += 1
		v_stat = fit_dict[env][clone]["virus"]
		dip_stat = fit_dict[env][clone]["hapdip"]
		
		if v_stat == "Virus retained" and dip_stat == "Haploid":
			vdm_count_mat[-1][0] += 1
		elif v_stat == "Virus retained" and dip_stat == "Undetermined":
			vdm_count_mat[-1][1] += 1
		elif v_stat == "Virus retained" and dip_stat == "Diploid":
			vdm_count_mat[-1][2] += 1
		elif v_stat == "Virus lost" and dip_stat == "Haploid":
			vdm_count_mat[-1][3] += 1
		elif v_stat == "Virus lost" and dip_stat == "Undetermined":
			vdm_count_mat[-1][4] += 1
		elif v_stat == "Virus lost" and dip_stat == "Diploid":
			vdm_count_mat[-1][5] += 1
		else:
			print(env, clone, v_stat, dip_stat)

vdm_count_mat = numpy.array(vdm_count_mat)
print(numpy.sum(clones_per_env))
print(numpy.sum(vdm_count_mat, axis=0))
print(numpy.sum(vdm_count_mat, axis=0)/numpy.sum(clones_per_env))
print(vdm_count_mat)
print(env_list)
chi2, p, dof, expected = scipy.stats.chi2_contingency(vdm_count_mat)
print(chi2, p, dof)

v_count_mat = numpy.zeros((vdm_count_mat.shape[0], 2))
v_count_mat[:,0] = numpy.sum(vdm_count_mat[:,0:3], axis=1)
v_count_mat[:,1] = numpy.sum(vdm_count_mat[:,3:], axis=1)

chi2vm, pvm, dofvm, expected = scipy.stats.chi2_contingency(v_count_mat)
print(chi2vm, pvm, dofvm)

hd_count_mat = numpy.zeros((vdm_count_mat.shape[0], 2))
hd_count_mat[:,0] = vdm_count_mat[:,0] + vdm_count_mat[:,3]
hd_count_mat[:,1] = vdm_count_mat[:,2] + vdm_count_mat[:,5]

chi2vm, pvm, dofvm, expected = scipy.stats.chi2_contingency(hd_count_mat)
print(chi2vm, pvm, dofvm)

####


