import matplotlib.pylab as pt
import numpy
import matplotlib.gridspec as gridspec
import matplotlib
import scipy.stats
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import pearsonr
import seaborn as sns
from numpy.random import permutation
import import_utilities

def permute_constraining_sums(table):

	####Construct list of all (env, gene) pairs
	nenvs, ngenes = table.shape
	env_label_list = []
	gene_label_list = []
	
	for i in range(nenvs):
		for j in range(ngenes):
			nhits = table[i,j]
			for k in range(nhits):
				env_label_list.append( i )
				gene_label_list.append( j )
	
	###Permute genes
	permuted_gene_labels = permutation( gene_label_list )
	
	###Figure out which envs they are now assigned to
	
	new_table = numpy.zeros_like(table)
	
	for k in range(len(gene_label_list)):
		i = env_label_list[k]
		j = permuted_gene_labels[k]
		new_table[i,j] += 1
	
	###Check on row and column sums 
	
	new_table_sum1 = numpy.sum(new_table, axis=0)
	new_table_sum2 = numpy.sum(new_table, axis=1)
	table_sum1 = numpy.sum(table, axis=0)
	table_sum2 = numpy.sum(table, axis=1)
	
	if (new_table_sum1 != table_sum1).any():
		print('Something is wrong on axis0 sum, ', table_sum1, new_table_sum1)
	if (new_table_sum2 != table_sum2).any():
		print('Something is wrong on axis1 sum, ', table_sum2, new_table_sum2)
	
	return new_table

def calculate_MI(table, clones_per_env):
	
	nenvs, ngenes = table.shape
	muts_per_gene = numpy.sum(table, axis=0)
	tot_num_clones = numpy.sum(clones_per_env)
	
	
	MI = 0
	
	for j in range(ngenes):
		p_gene = muts_per_gene[j]/tot_num_clones
		
		for i in range(nenvs):
			p_gene_given_env = table[i,j]/clones_per_env[i]
			
			p_env = clones_per_env[i]/tot_num_clones
			#p_joint = table[i,j]/clones_per_env[i]*p_env
			if p_gene_given_env > 0:
				MI_increment = p_env*(p_gene_given_env*numpy.log2(p_gene_given_env/p_gene) + (1 - p_gene_given_env)*numpy.log2((1 - p_gene_given_env)/(1-p_gene)))
				#MI_increment = p_joint*numpy.log2(p_gene_given_env/p_gene) + (1 - p_joint)*numpy.log2((1 - p_gene_given_env)/(1-p_gene))
				MI += MI_increment
	
	return MI

###Data import

gene_dict, clone_mut_num_dict, clones_per_env_dict, mut_dict = import_utilities.import_mutation_data(file_in='data/mutation_table_7_3_2019.txt')
family_list_go_ordered, go_slim_data = import_utilities.import_GO_slim(file_in = 'data/GO_Slim_output_7_4_2019.txt')


mut_counts = {}
mut_list = []
for env in mut_dict:
	
	for clone in mut_dict[env]:
		
		for mut in mut_dict[env][clone]:
			
			if mut in mut_list:
				
				
				if tuple(mut) not in mut_counts:
					mut_counts[tuple(mut)] = 2
				else:
					mut_counts[tuple(mut)] += 1
			else:
				mut_list.append(mut)
####
n_shared = 0
for mut in mut_counts:
	n_shared += mut_counts[mut]
	
print('Number of unique mutations shared amongst multiple clones',len(mut_counts))
print('Number of shared mutations', n_shared)
print('Number of unique mutations', len(mut_list) - len(mut_counts))
print('Sanity check of the total number of mutations', n_shared + len(mut_list) - len(mut_counts))
print('Fraction of shared mutations', n_shared/(len(mut_list) - len(mut_counts)))




###List genes mutated 4 or more times across the experiment or 2 or more times within an environment

gene_list = []
gene_num_hits = []
twohit_within_list = []
counter = 0
for gene in gene_dict:
	gene_list.append(gene)
	gene_num_hits.append(0)
	for env in gene_dict[gene]:
		
		gene_num_hits[counter] += numpy.sum( gene_dict[gene][env] )
		
		if numpy.sum(gene_dict[gene][env]) > 1.5 and gene not in twohit_within_list: ###2 or more hits within envt
			
			twohit_within_list.append(gene)
			
	counter += 1
	
gene_order = numpy.argsort(gene_num_hits)[::-1]

gene_num_hits = numpy.array(gene_num_hits)

gene_num_hits_ordered = gene_num_hits[ gene_order ]

three_hits_plus_inds = gene_order[ numpy.greater( gene_num_hits_ordered, 3.5 ) ]


gene_list_threeplus = numpy.array(gene_list)[ three_hits_plus_inds ]
print(len(gene_list_threeplus))
print(gene_list_threeplus)
gene_num_hits_threeplus = gene_num_hits[ three_hits_plus_inds ]
print(sum(gene_num_hits_threeplus))
two_hit_within_unique = list(set(twohit_within_list)-set(gene_list_threeplus))
print(two_hit_within_unique)
print('Number of significant genes (4 hit overall or 2 hit within envt), ',len(two_hit_within_unique) + len(gene_list_threeplus))
print('Number of clones with mutations in these genes, ', sum(gene_num_hits_threeplus) + 2*len(two_hit_within_unique))

nsiggenes = len( list(set(twohit_within_list)-set(gene_list_threeplus) )) + len(gene_list_threeplus)

sig_gene_list = list(set(gene_list_threeplus).union(set(twohit_within_list)))

outfile = open('data/significant_genes_check.txt','w')
outfile.write((' ').join(sig_gene_list))
outfile.close()


####Initialize with the Ras family, which the go slim analysis does not pick up, so must be added by hand
gene_families = {'Ras':['IRA1','IRA2','GPB1','GPB2','PDE2','CDC25']}
gene_list_by_families = ['IRA1','IRA2','GPB1','GPB2','PDE2','CDC25']
families = ['Ras']
family_borders = [0,6]

#gene_families = {'Ras':['IRA1','IRA2','GPB1','GPB2','PDE2','CDC25'],'Chromatin org.':['SIR1','SIR2','SIR3','SIR4','ARP5','ORC5'],'Transmem. transport':['HNM1','ENA2','SMF1','SMF2','VHT1','MCX1'],'Cell wall':['ROT2','CCW12','KRE6','ANP1','FKS1','OSW2']}
#gene_list_by_families = ['IRA1','IRA2','GPB1','GPB2','PDE2','CDC25','SIR1','SIR2','SIR3','SIR4','ARP5','ORC5','HNM1','ENA2','SMF1','SMF2','VHT1','MCX1','ROT2','CCW12','KRE6','ANP1','FKS1','OSW2']
#families = ['Ras','Chromatin org.','Transmem. transport','Cell wall']
#family_borders = [0,6,12,18,24]

ngenes = family_borders[-1]

for family in family_list_go_ordered:
	genes = go_slim_data[family]['genes']
	nunique = len( list( set(genes) - set(gene_list_by_families) ) )
	
	if nunique > 1.5:
		unique_fam_genes = list( set(genes) - set(gene_list_by_families) )
	
		gene_families[family] = unique_fam_genes
		gene_list_by_families.extend(unique_fam_genes)
		family_borders.append(family_borders[-1] + len(unique_fam_genes))
		families.append(family)
		ngenes += len(unique_fam_genes)

unassigned_genes = list(set(sig_gene_list) - set(gene_list_by_families))


gene_family_lookup = {}

for gene in gene_list_by_families:
	for family in gene_families:
		if gene in gene_families[family]:
			gene_family_lookup[gene] = family

for gene in unassigned_genes:
	gene_family_lookup[gene] = 'other'		

families.append('other')
family_borders.append( family_borders[-1] + len(unassigned_genes) )
gene_list_by_families.extend(unassigned_genes)

family_labels = []
for family in families:
	family_strs = family.split(' ')
	new_family_strs = ('\n').join(family_strs)
	family_labels.append(new_family_strs)

family_inds = dict(zip(families,range(len(families))))

####


evolution_environments = ['Evolved in pH 3', 'Evolved in pH 3.8', 'Evolved in 37 C', 'Evolved in .8 M','Evolved in pH 6', 'Evolved in pH 7.3', 'Evolved in SC', 'Evolved in .07% glu','Evolved in gal',  'Evolved in 21 C', 'Evolved in .2 M']
envs_short = ['pH 3','pH 3.8','High\ntemp','High\nsalt','pH 6','pH 7.3','SC','Low glu','Gal','Low\ntemp','Med\nsalt']

#evolution_environments = ['Evolved in pH 3', 'Evolved in pH 3.8', 'Evolved in pH 6', 'Evolved in pH 7.3','Evolved in SC', 'Evolved in .07% glu','Evolved in gal', 'Evolved in 37 C', 'Evolved in 21 C','Evolved in .8 M', 'Evolved in .2 M']
#envs_short = ['pH 3','pH 3.8','pH 6','pH 7.3','SC','Low glu','Gal','High\ntemp','Low\ntemp','High\nsalt','Med\nsalt']
gene_family_matrix = numpy.zeros( (len(gene_list_by_families), len(evolution_environments)),dtype=float)
gene_family_matrix_counts = numpy.zeros( (len(gene_list_by_families), len(evolution_environments)),dtype=int)
gene_family_matrix_summed = numpy.zeros( (len(gene_list_by_families), len(evolution_environments)),dtype=float)
#gene_family_matrix_nonsense = numpy.zeros( (len(families), len(evolution_environments)),dtype=int)

####Make a table that lists the fraction of clones from a particular environment that got a nonsynonymous mutation in a particular gene

num_clones_per_env = []
for i in range(len(evolution_environments)):
	env = evolution_environments[i]
	num_clones_per_env.append(clones_per_env_dict[env])
	
for gene in gene_list_by_families:

	family = gene_family_lookup[gene]
	family_ind = family_inds[family]
	gene_ind = gene_list_by_families.index(gene)
	
	for i in range(len(evolution_environments)):
		env = evolution_environments[i]
		if env in gene_dict[gene]:
			
			gene_family_matrix[gene_ind, i] += (gene_dict[gene][env][0] + gene_dict[gene][env][1])/clones_per_env_dict[env]###All nonsyn muts
			gene_family_matrix_counts[gene_ind, i] += (gene_dict[gene][env][0] + gene_dict[gene][env][1])
			#gene_family_matrix_nonsense[family_ind, i] += gene_dict[gene][env][1]

MI_base = calculate_MI( gene_family_matrix_counts.T, num_clones_per_env )
MI_list = []
niter = 10000

for n in range(niter):
	
	new_table = permute_constraining_sums( gene_family_matrix_counts.T )
	
	MI = calculate_MI(new_table, num_clones_per_env)		
	MI_list.append(MI)
	
p = numpy.sum( numpy.array(MI_list) > MI_base )/niter
print(p)
print(MI_base)
print(numpy.mean(MI_list))
pt.figure()
pt.hist(MI_list)

pt.show()