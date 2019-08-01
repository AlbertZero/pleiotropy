import matplotlib.pylab as pt
import numpy
import matplotlib.gridspec as gridspec
import matplotlib
import scipy.stats
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import pearsonr
import seaborn as sns
from numpy.random import multinomial
import import_utilities

####This script imports the lengths of the yeast ORFS and redistributes nonsynonymous genic de novo mutations amongst ORFS to test for parallelism.

def import_w303_gene_lengths(gff = 'data/w303_ref.gff'):

	gene_length_dict = {}
	
	file = open(gff,'r')
	for line in file:
		if line.startswith('#'):
			continue
		linelist = line.strip().split('\t')
		
		entry_type = linelist[2]
		if entry_type == 'CDS':
			
			cds_start = int(linelist[3])
			cds_end = int(linelist[4])
			annot_dict = {}
			
			for item in linelist[-1].split(';'):
				item_list = item.split('=')
				if len(item_list) > 1.5 and item_list[1] != '':
				
					annot_dict[item_list[0]] = item_list[1]
			
			if 'gene' in annot_dict:
				gene_name = annot_dict['gene']
			else:
				gene_name = annot_dict['Parent']
				
			if gene_name not in gene_length_dict:
				gene_length_dict[gene_name] = cds_end - cds_start
			
			else:
				gene_length_dict[gene_name] += cds_end - cds_start
		
	return gene_length_dict

gene_length_dict = import_w303_gene_lengths(gff = 'data/w303_ref.gff')	
###Import mutation data
gene_dict, clone_mut_num_dict, clones_per_env_dict, mut_dict = import_utilities.import_mutation_data(file_in='data/mutation_table_7_3_2019.txt')

num_nonsyn_genes = 0
for gene in gene_dict:
	for env in gene_dict[gene]:
		num_nonsyn_genes += numpy.sum(gene_dict[gene][env])
		print(gene_dict[gene][env])
print(num_nonsyn_genes)

###Redistribute the observed nonsynonymous mutations amongst the genes

gene_length_list = numpy.array([val for key,val in gene_length_dict.items()])

fractional_lengths = gene_length_list/numpy.sum(gene_length_list)

niter = 1000
nbins = 8
####
muts_per_env = []
env_list = []
for env in clone_mut_num_dict:
	nmuts = 0
	env_list.append(env)
	for clone in clone_mut_num_dict[env]:
		nmuts += clone_mut_num_dict[env][clone]
	muts_per_env.append(nmuts)
	
within_env_dists = []
overall_dists = []

for i in range(niter):
	redist_envs = []
	for j in range(len(env_list)):
		redist_env = multinomial(muts_per_env[j], fractional_lengths, 1)
		redist_envs.append(redist_env[0,:])
	
	redist_envs = numpy.array(redist_envs)
	
	####Construct the expected number of certain numbers of hits, within each environment and across the whole experiment
	
	nwithin = []
	noverall = []
	for nhits in range(nbins):
		
		ngenes_within = numpy.mean( numpy.sum(redist_envs == nhits, axis=1) )
		ngenes_overall = numpy.sum( numpy.sum(redist_envs, axis=0) == nhits )
		
		nwithin.append(ngenes_within)
		noverall.append(ngenes_overall)
	
	within_env_dists.append(nwithin)
	overall_dists.append(noverall)

within_env_dists = numpy.array(within_env_dists)
overall_dists = numpy.array(overall_dists)

####
mut_mat = numpy.zeros((len(env_list), len(gene_dict)),dtype='int')
gene_list = list(gene_dict.keys())
for i in range(len(gene_list)):
	gene = gene_list[i]
	for j in range(len(env_list)):
		env = env_list[j]
		if env in gene_dict[gene]:
			nmuts = numpy.sum(gene_dict[gene][env])
		
			mut_mat[j,i] = nmuts
####

within_env_dist_obs = []
overall_dist_obs = []

for nhits in range(nbins):
		
	ngenes_within = numpy.mean( numpy.sum(mut_mat == nhits, axis=1) )
	ngenes_overall = numpy.sum( numpy.sum(mut_mat, axis=0) == nhits )
		
	within_env_dist_obs.append(ngenes_within)
	overall_dist_obs.append(ngenes_overall)


fig, (ax1,ax2) = pt.subplots(1,2,figsize=(8,4))

shift = .2

ax1.bar(numpy.arange(nbins-1)-shift+1,numpy.mean(overall_dists,axis=0)[1:],width=.3)
ax1.bar(numpy.arange(nbins-1)+shift+1,overall_dist_obs[1:],width=.3)
ax1.set_xticks(numpy.arange(nbins-1)+1)
ax1.set_ylim(0,20)
ax1.legend(['Null','Observed'])
ax1.set_xlabel('Hits per gene')
ax1.set_ylabel('Number of genes')

ax1.text(-.05,1.05,'A',fontname='Arial',transform=ax1.transAxes,fontsize=14)

ax2.bar(numpy.arange(nbins-1)-shift+1,numpy.mean(within_env_dists,axis=0)[1:],width=.3)
ax2.bar(numpy.arange(nbins-1)+shift+1,within_env_dist_obs[1:],width=.3)
ax2.set_xticks(numpy.arange(nbins-1)+1)
ax2.set_ylim(0,20)
ax2.legend(['Null','Observed'])
ax2.set_xlabel('Hits per gene')
ax2.set_ylabel('Average number of genes per environment')
ax2.text(-.05,1.05,'B',fontname='Arial',transform=ax2.transAxes,fontsize=14)
pt.tight_layout()
pt.savefig('figures/gene_hit_enrichment_check.pdf',bbox_inches='tight')