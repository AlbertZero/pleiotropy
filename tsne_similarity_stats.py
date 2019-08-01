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
from sklearn.neighbors import KNeighborsClassifier

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

def find_neighborhood_sharing(neighbors_array, env_list):
	
	nclones = len(env_list)
	env_list = numpy.array(env_list)
	sharing = []
	for i in range(nclones):
		env = env_list[i]
		neighbors = neighbors_array[i,:]
		neighbor_envs = env_list[neighbors]
		nshared = numpy.sum(neighbor_envs == env)
		sharing.append(nshared)
		
	return sharing

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

fit_mat_vminus = []
hd_list_vminus = []
env_list_vminus = []

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
		 
				
	fit_dict_extr[evol_env] = numpy.array( fit_dict_extr[evol_env] )
	err_dict_extr[evol_env] = numpy.array(err_dict_extr[evol_env] )
	
fit_mat = numpy.array(fit_mat) ##All the clone fitness vectors concatenated

print(len(v_list))
v_list = numpy.array(v_list)
hd_list = numpy.array(hd_list)

print( numpy.sum(numpy.logical_and(v_list == "Virus retained", hd_list == "Diploid")) )
print( numpy.sum(numpy.logical_and(v_list == "Virus lost", hd_list == "Diploid")))
print( numpy.sum(numpy.logical_and(v_list == "Virus retained", hd_list == "Haploid")) )
print( numpy.sum(numpy.logical_and(v_list == "Virus lost", hd_list == "Haploid")))
print( numpy.sum(numpy.logical_and(v_list == "Virus retained", hd_list == "Undetermined")) )
print( numpy.sum(numpy.logical_and(v_list == "Virus lost", hd_list == "Undetermined")))

	
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

####Look at the k nearest neighbors of each point, and ask what proportion of them are from the same environment, compared with random permutations of labels

neigh = KNeighborsClassifier(n_neighbors=5)
neigh.fit(fit_mat_tsne,numpy.ones((fit_mat_tsne.shape[0],)))
neighbors = neigh.kneighbors(return_distance=False)
print(neighbors)
shared_env_neighbors = find_neighborhood_sharing(neighbors, env_list)
print(numpy.mean(shared_env_neighbors))

npermute = 1000
nshared = []
for i in range(npermute):
	
	permuted_env_list = numpy.random.permutation(env_list)
	shared_env_neighbors = find_neighborhood_sharing(neighbors, permuted_env_list)
	nshared.append(numpy.mean(shared_env_neighbors))

pt.hist(nshared)
print(numpy.percentile(nshared,50),numpy.percentile(nshared,5),numpy.percentile(nshared,95))
pt.show()
