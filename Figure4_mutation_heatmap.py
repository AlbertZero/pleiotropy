import matplotlib.pylab as pt
import numpy
import matplotlib.gridspec as gridspec
import matplotlib
import scipy.stats
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import pearsonr
import seaborn as sns
import import_utilities

gene_dict, clone_mut_num_dict, clones_per_env_dict, mut_dict = import_utilities.import_mutation_data(file_in='data/mutation_table_7_3_2019.txt')
family_list_go_ordered, go_slim_data = import_utilities.import_GO_slim(file_in = 'data/GO_Slim_output_7_4_2019.txt')

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


####Initialize with the Ras family, which the go slim analysis does not pick up, so must be added by hand
gene_families = {'Ras':['IRA1','IRA2','GPB1','GPB2','PDE2']}
gene_list_by_families = ['IRA1','IRA2','GPB1','GPB2','PDE2']
families = ['Ras']
family_borders = [0,5]


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


evolution_environments = ['Evolved in pH 3', 'Evolved in pH 3.8', 'Evolved in pH 6', 'Evolved in pH 7.3', 'Evolved in 37 C', 'Evolved in .8 M', 'Evolved in SC', 'Evolved in .07% glu','Evolved in gal',  'Evolved in 21 C', 'Evolved in .2 M']
envs_short = ['pH 3','pH 3.8','pH 6','pH 7.3','High\ntemp','High\nsalt','SC','Low glu','Gal','Low\ntemp','Med\nsalt']

gene_family_matrix = numpy.zeros( (len(gene_list_by_families), len(evolution_environments)),dtype=float)
gene_family_matrix_summed = numpy.zeros( (len(gene_list_by_families), len(evolution_environments)),dtype=float)

####Make a table that lists the fraction of clones from a particular environment that got a nonsynonymous mutation in a particular gene

for gene in gene_list_by_families:

	family = gene_family_lookup[gene]
	family_ind = family_inds[family]
	gene_ind = gene_list_by_families.index(gene)
	
	for i in range(len(evolution_environments)):
		env = evolution_environments[i]
		if env in gene_dict[gene]:
			
			gene_family_matrix[gene_ind, i] += (gene_dict[gene][env][0] + gene_dict[gene][env][1])/clones_per_env_dict[env]###All nonsyn muts
			
			
family_labels = {'Ras':'Ras','invasive growth in response to glucose limitation':'Invasive\ngrowth\nin gluc lim','conjugation':'Conjugation','response to osmotic stress':'Osmotic\nstress','cellular ion homeostasis':'Ion\nhomeostasis','response to chemical':'Resp. to\nchemical','protein modification by small protein conjugation or removal':'Protein\nmodification','chromatin organization':'Chromatin\norganization','cell wall organization or biogenesis':'Cell wall','lipid metabolic process':'Lipid met.\nprocess','other':'Other'}#'Transcription\nby RNA Pol II',

pt.figure(figsize=(8,6))
print(gene_family_matrix)
zz = numpy.ma.masked_where(gene_family_matrix < 0.001, gene_family_matrix)

my_plot=pt.pcolormesh(numpy.flip(zz,axis=0),vmin=0, vmax=numpy.max(gene_family_matrix), cmap='Reds',edgecolors='w',linewidth=.5)
my_plot.cmap.set_bad("white")
ax = pt.gca()

ax.set_yticks(numpy.arange(gene_family_matrix.shape[0])+.5)
ax.set_xticks(numpy.arange(gene_family_matrix.shape[1])+.5)
ax.set_xticklabels(envs_short,fontsize=8)
ax.set_yticklabels(gene_list_by_families[::-1], fontsize=5)
ax.tick_params(axis='both',length=0)

ny_grid = gene_family_matrix.shape[0]

for i in range(len(families)):
	lind = family_borders[i]
	uind = family_borders[i+1]
	family = family_labels[families[i]] ###Translation of the GO slim category to a slightly shorter name for plotting
	
	ax.axhline(ny_grid - uind,-.24,1,color='k',linewidth=.5,clip_on=False)
	
	ax.text(-1.6,ny_grid - (uind+lind)/2,family,verticalalignment='center',horizontalalignment='center',fontsize=7,fontvariant='small-caps')

ax.axhline(ny_grid,-.24,1,color='k',linewidth=.5,clip_on=False)

for i in range(len(envs_short)-1):
	ax.axvline(i+1,0,1,color='k',linewidth=.5)

cbar = pt.colorbar()
cbar.set_label('Fraction of clones')

pt.savefig('figures/mutation_heatmap_check.pdf',bbox_inches='tight')