import matplotlib.pylab as pt
import numpy
from sklearn.linear_model import LinearRegression
import import_utilities
import matplotlib

# Set the default color cycle
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]) 

###This script explores a series of linear models to predict the fitness increment of a clone after adaptation.
###In particular, the fitness of clone j in environment i, Yij, will be predicted as follows:

###Y_ij = a_i + sum_k alpha_ki E_j + zeta_i V_j + beta_i H_j + sum_k gamma_ki G_j + eps_ij
###Where E_j are indicator variables on the environment where the clone evolved; V_j are indicator variables on virus status; H_j are indicator variables on haploid/diploid status; and G_kj are indicator variables on having a mutation in one of the gene categories.
###We will explore the fraction of variance predicted by environment alone; virus status alone; and the full model for the V+ clones (because the V- ones have frequency-dependent fitnesses), in each measurement environment

output_stats_file = 'data/linear_model_output.txt'

##Import fitness data

fit_dict, gen_dict, std_err_dict, random_reps = import_utilities.import_fitness_data('data/pleiotropy_fitness_data_7_4_2019.txt')

##Make linear predictor plot with fitness data
evolution_environments = ['Evolved in SC','Evolved in gal','Evolved in 21 C','Evolved in 37 C', 'Evolved in .07% glu', 'Evolved in pH 3','Evolved in pH 7.3','Evolved in .8 M']
envt_labels = ['SC','Gal','Low\ntemp','High\ntemp','Low\nglu','pH 3','pH 7.3','High\nsalt']
envt_labels_oneline = ['SC','Gal','Low temp','High temp','Low glu','pH 3','pH 7.3','High salt']
measurement_environments = ['fitness in SC','fitness in gal','fitness in 21 C','fitness in 37 C', 'fitness in .07% glu', 'fitness in pH 3','fitness in pH 7.3','fitness in .8 M']

##

###Predict fitnesses for all non-home clones based on -their evolution environment, -their virus status, All remaining variance not due to measurement error is evolutionary stochasticity.

##

nbootstrap = 500

vV_allenv = []
vE_allenv = []
vVE_allenv = []
vVED_allenv = []

vV_allenv_err = []
vE_allenv_err = []
vVE_allenv_err = []
vVED_allenv_err = []

err_var_by_env = []

for i in range(len(evolution_environments)):
	
	meas_env = measurement_environments[i]
	predictor_mat_all = []
	
	outcome_list = []
	meas_err_list = []
	for j in range(len(evolution_environments)):
		if j != i:
			evol_env = evolution_environments[j]
			gen_norm_factor = gen_dict[evol_env]
			for clone in fit_dict[evol_env]:
			
				clone_fit = fit_dict[evol_env][clone][meas_env]*700/gen_norm_factor
				clone_err = std_err_dict[evol_env][clone][meas_env]*700/gen_norm_factor
			
				outcome_list.append(clone_fit)
				meas_err_list.append(clone_err)
			
				
				predictor_list = [0.,0.,0.,0.,0.,0.,0.,0.,0.]
				
				predictor_list[j] = 1 ###Indicator variable on the evolution environment
			
				if fit_dict[evol_env][clone]['virus'] == 'Virus retained':
					v_indicator = 1.
					#print('found V')
				else:
					v_indicator = 0.
				#print(fit_dict[evol_env][clone]['hapdip'])
				#if fit_dict[evol_env][clone]['hapdip'] == 'Diploid':
				#	d_indicator = 1.
					#print('found D')
				#else:
				#	d_indicator = 0.
			
				predictor_list[8] = v_indicator ##Indicator variable on virus status
				#predictor_list[9] = d_indicator ##Indicator variable on diploid status
				
				predictor_mat_all.append(predictor_list)
	
	predictor_mat_all = numpy.array(predictor_mat_all)
	
	outcome_list = numpy.array(outcome_list)
	
	meas_err_list = numpy.array(meas_err_list)
	meas_err_var = numpy.mean(meas_err_list**2)/numpy.var(outcome_list)
	err_var_by_env.append(meas_err_var)
	
	nclones = predictor_mat_all.shape[0]
	vV_bs = []
	vE_bs = []
	vVE_bs = []
	vVED_bs = []
	
	vreg = LinearRegression().fit(predictor_mat_all[:,8].reshape(-1,1), outcome_list[:])
	ereg = LinearRegression().fit(predictor_mat_all[:,:8], outcome_list[:])
	vereg = LinearRegression().fit(predictor_mat_all[:,:], outcome_list[:])
	
	vV = vreg.score(predictor_mat_all[:,8].reshape(-1,1), outcome_list[:])
	vE = ereg.score(predictor_mat_all[:,:8], outcome_list[:])
	vVE = vereg.score(predictor_mat_all[:,:], outcome_list[:])
		
	# for n in range(nbootstrap):
# 		bs_inds = numpy.random.choice(numpy.arange(nclones),size=nclones)
# 		vreg = LinearRegression().fit(predictor_mat_all[bs_inds,8].reshape(-1,1), outcome_list[bs_inds])
# 		ereg = LinearRegression().fit(predictor_mat_all[bs_inds,:8], outcome_list[bs_inds])
# 		vereg = LinearRegression().fit(predictor_mat_all[bs_inds,:9], outcome_list[bs_inds])
# 		vedreg = LinearRegression().fit(predictor_mat_all[bs_inds,:], outcome_list[bs_inds])
# 		
# 		vV = vreg.score(predictor_mat_all[bs_inds,8].reshape(-1,1), outcome_list[bs_inds])
# 		vE = ereg.score(predictor_mat_all[bs_inds,:8], outcome_list[bs_inds])
# 		vVE = vereg.score(predictor_mat_all[bs_inds,:9], outcome_list[bs_inds])
# 		vVED = vedreg.score(predictor_mat_all[bs_inds,:], outcome_list[bs_inds])
# 		
# 		vV_bs.append(vV)
# 		vE_bs.append(vE)
# 		vVE_bs.append(vVE)
# 		vVED_bs.append(vVED)
# 		
# 		vVm = numpy.percentile(vV_bs,50)
# 		vV_err = [vVm-numpy.percentile(vV_bs,25),numpy.percentile(vV_bs,75)-vVm]
# 		vEm = numpy.percentile(vE_bs,50)
# 		vE_err = [vEm-numpy.percentile(vE_bs,25),numpy.percentile(vE_bs,75)-vEm]
# 		vVEm = numpy.percentile(vVE_bs,50)
# 		vVE_err = [vVEm-numpy.percentile(vVE_bs,25),numpy.percentile(vVE_bs,75)-vVEm]
# 		vVEDm = numpy.percentile(vVED_bs,50)
# 		vVED_err = [vVEDm-numpy.percentile(vVED_bs,25),numpy.percentile(vVED_bs,75)-vVEDm]
# 		
	# vV_allenv.append(vVm)
# 	vV_allenv_err.append(vV_err)
# 	vE_allenv.append(vEm)
# 	vE_allenv_err.append(vE_err)
# 	vVE_allenv.append(vVEm)
# 	vVE_allenv_err.append(vVE_err)
# 	vVED_allenv.append(vVEDm)
# 	vVED_allenv_err.append(vVED_err)
	
	vV_allenv.append(vV)
	vE_allenv.append(vE)
	vVE_allenv.append(vVE)

vV_allenv = numpy.array(vV_allenv)
vE_allenv = numpy.array(vE_allenv)
vVE_allenv = numpy.array(vVE_allenv)
#vVED_allenv = numpy.array(vVED_allenv)

# vV_allenv_err = numpy.array(vV_allenv_err)
# vE_allenv_err = numpy.array(vE_allenv_err)
# vVE_allenv_err = numpy.array(vVE_allenv_err)
# vVED_allenv_err = numpy.array(vVED_allenv_err)
err_var_by_env = numpy.array(err_var_by_env)

#print(vV_allenv_err.shape)

pt.figure(figsize=(6,5))

width = .8
vbar = pt.bar(numpy.arange(8),vV_allenv,width,color='C2')
ebar = pt.bar(numpy.arange(8),vVE_allenv-vV_allenv,width,bottom=vV_allenv,color='C0')
#dbar = pt.bar(numpy.arange(8),vVED_allenv-vVE_allenv,width,bottom=vVE_allenv,color='C3')
stochbar = pt.bar(numpy.arange(8),1-err_var_by_env-vVE_allenv,width,bottom=vVE_allenv,color='C4')
errbar = pt.bar(numpy.arange(8),err_var_by_env,width,bottom=1-err_var_by_env,color='C5')
ax = pt.gca()
ax.set_xticks(numpy.arange(8))
ax.set_xticklabels(envt_labels)
ax.set_ylabel('Proportion of variance',fontsize=12)
ax.set_xlabel('Measurement environment',fontsize=12)
##ax.set_xlim(-1,10)
pt.legend([vbar[0],ebar[0],stochbar[0],errbar[0]],['V+/V-','Home Env','Idiosyncratic\npleiotropy','Meas Err'],loc=[.65,.6],fontsize=9)

pt.savefig('figures/fitness_predictors_nohome_vfirst_11_21_2019.pdf',bbox_inches='tight')

pt.figure(figsize=(6,5))

width = .8
ebar = pt.bar(numpy.arange(8),vE_allenv,width,color='C0')
vbar = pt.bar(numpy.arange(8),vVE_allenv-vE_allenv,width,bottom=vE_allenv,color='C2')
#dbar = pt.bar(numpy.arange(8),vVED_allenv-vVE_allenv,width,bottom=vVE_allenv,color='C3')
stochbar = pt.bar(numpy.arange(8),1-err_var_by_env-vVE_allenv,width,bottom=vVE_allenv,color='C4')
errbar = pt.bar(numpy.arange(8),err_var_by_env,width,bottom=1-err_var_by_env,color='C5')
ax = pt.gca()
ax.set_xticks(numpy.arange(8))
ax.set_xticklabels(envt_labels)
ax.set_ylabel('Proportion of variance',fontsize=12)
ax.set_xlabel('Measurement environment',fontsize=12)
#ax.set_xlim(-1,10)
pt.legend([ebar[0],vbar[0],stochbar[0],errbar[0]],['Home Env','V+/V-','Other\npleiotropy','Meas Err'],loc=[.65,.6],fontsize=9)

pt.savefig('figures/fitness_predictors_nohome_efirst_11_21_2019.pdf',bbox_inches='tight')

file = open(output_stats_file,'w')
file.write('Model' + '\t' + ('\t').join(envt_labels_oneline) + '\n')
file.write('Evolution env.' + '\t' + ('\t').join([str(p) for p in vE_allenv]) + '\n')
file.write('Evolution env. + Virus' + '\t' + ('\t').join([str(p) for p in vVE_allenv]) + '\n')
file.write('Virus (additional)' + '\t' + ('\t').join([str(p) for p in vVE_allenv-vE_allenv]) + '\n')
file.write('All stochasticity' + '\t' + ('\t').join([str(p) for p in 1 - vE_allenv-err_var_by_env]) + '\n')
file.write('Other pleiotropy' + '\t' + ('\t').join([str(p) for p in 1 - vVE_allenv-err_var_by_env]) + '\n')
file.write('Measurement err.' + '\t' + ('\t').join([str(p) for p in err_var_by_env]) + '\n')

print('Env variance ', vE_allenv)
print(numpy.mean(vE_allenv))
print('Stoch variance ', 1 - vE_allenv - err_var_by_env)
print(numpy.mean(1 - vE_allenv - err_var_by_env))
###Now include mutation data. We will have an indicator variable for the presence/absence of a nonsynonymous mutation in one of the categories we identified, for each clone.
##Import mutation data

gene_dict, clone_mut_num_dict, clones_per_env_dict, mut_dict = import_utilities.import_mutation_data(file_in='data/mutation_table_7_3_2019.txt')
family_list_go_ordered, go_slim_data = import_utilities.import_GO_slim(file_in = 'data/GO_Slim_output_7_4_2019.txt')

#print(gene_dict)
print(mut_dict.keys())
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
#print(len(gene_list_threeplus))
#print(gene_list_threeplus)
gene_num_hits_threeplus = gene_num_hits[ three_hits_plus_inds ]
#print(sum(gene_num_hits_threeplus))
two_hit_within_unique = list(set(twohit_within_list)-set(gene_list_threeplus))
#print(two_hit_within_unique)
#print('Number of significant genes (4 hit overall or 2 hit within envt), ',len(two_hit_within_unique) + len(gene_list_threeplus))
#print('Number of clones with mutations in these genes, ', sum(gene_num_hits_threeplus) + 2*len(two_hit_within_unique))

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

#####

vG_allenv = []
vG_allenv_err = []
vGD_allenv = []
vGD_allenv_err = []
vGDV_allenv = []
vGDV_allenv_err = []
vGDVE_allenv = []
vGDVE_allenv_err = []

err_var_by_env = []

for i in range(len(evolution_environments)):
	
	meas_env = measurement_environments[i]
	predictor_mat_all = []
	
	outcome_list = []
	meas_err_list = []
	for j in range(len(evolution_environments)):
		if j != i:
			evol_env = evolution_environments[j]
			gen_norm_factor = gen_dict[evol_env]
			for clone in fit_dict[evol_env]:
				if clone in mut_dict[evol_env]: #fit_dict[evol_env][clone]['virus'] == 'Virus retained' and
					
					muts = mut_dict[evol_env][clone]
					
					clone_fit = fit_dict[evol_env][clone][meas_env]*700/gen_norm_factor
					clone_err = std_err_dict[evol_env][clone][meas_env]*700/gen_norm_factor
			
					outcome_list.append(clone_fit)
					meas_err_list.append(clone_err)
			
					predictor_list = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
					predictor_list[j] = 1. ###Indicator variable on the evolution environment
			
				
					#print(fit_dict[evol_env][clone]['hapdip'])
					
					
					if fit_dict[evol_env][clone]['virus'] == 'Virus retained':
						v_indicator = 1.
						#print('found V')
					else:
						v_indicator = 0.
					#print(fit_dict[evol_env][clone]['hapdip'])
					if fit_dict[evol_env][clone]['hapdip'] == 'Diploid':
						d_indicator = 1.
						#print('found D')
					else:
						d_indicator = 0.
			
					predictor_list[8] = v_indicator ##Indicator variable on virus status
					predictor_list[9] = d_indicator ##Indicator variable on diploid status
				
					
					for m in muts:
						if m[4] in gene_family_lookup:
							cat = gene_family_lookup[m[4]]
							if cat != 'other':
								cat_index = family_inds[cat]
								predictor_list[10+cat_index] = 1
								
					predictor_mat_all.append(predictor_list)
	
	predictor_mat_all = numpy.array(predictor_mat_all)
	outcome_list = numpy.array(outcome_list)
	meas_err_list = numpy.array(meas_err_list)
	meas_err_var = numpy.mean(meas_err_list**2)/numpy.var(outcome_list)
	err_var_by_env.append(meas_err_var)
	
	nclones = predictor_mat_all.shape[0]
	vG_bs = []
	vGD_bs = []
	vGDV_bs = []
	vGDVE_bs = []
	
	vGreg = LinearRegression().fit(predictor_mat_all[:,10:], outcome_list[:])
	vGDreg = LinearRegression().fit(predictor_mat_all[:,9:], outcome_list[:])
	vGDVreg = LinearRegression().fit(predictor_mat_all[:,8:], outcome_list[:])
	vGDVEreg = LinearRegression().fit(predictor_mat_all[:,:], outcome_list[:])
	
	
	vG = vGreg.score(predictor_mat_all[:,10:], outcome_list[:])
	vGD = vGDreg.score(predictor_mat_all[:,9:], outcome_list[:])
	vGDV = vGDVreg.score(predictor_mat_all[:,8:], outcome_list[:])
	vGDVE = vGDVEreg.score(predictor_mat_all[:,:], outcome_list[:])
	
	# for n in range(nbootstrap):
# 		bs_inds = numpy.random.choice(numpy.arange(nclones),size=nclones)
# 		
# 		vGreg = LinearRegression().fit(predictor_mat_all[bs_inds,10:], outcome_list[bs_inds])
# 		vGDreg = LinearRegression().fit(predictor_mat_all[bs_inds,9:], outcome_list[bs_inds])
# 		vGDVreg = LinearRegression().fit(predictor_mat_all[bs_inds,8:], outcome_list[bs_inds])
# 		vGDVEreg = LinearRegression().fit(predictor_mat_all[bs_inds,:], outcome_list[bs_inds])
# 		
# 		
# 		vG = vGreg.score(predictor_mat_all[bs_inds,10:], outcome_list[bs_inds])
# 		vGD = vGDreg.score(predictor_mat_all[bs_inds,9:], outcome_list[bs_inds])
# 		vGDV = vGDVreg.score(predictor_mat_all[bs_inds,8:], outcome_list[bs_inds])
# 		vGDVE = vGDVEreg.score(predictor_mat_all[bs_inds,:], outcome_list[bs_inds])
# 		
# 		vG_bs.append(vG)
# 		vGD_bs.append(vGD)
# 		vGDV_bs.append(vGDV)
# 		vGDVE_bs.append(vGDVE)
# 		
# 		vGm = numpy.percentile(vG_bs,50)
# 		vG_err = [vGm-numpy.percentile(vG_bs,25),numpy.percentile(vG_bs,75)-vGm]
# 		vGDm = numpy.percentile(vGD_bs,50)
# 		vGD_err = [vGDm-numpy.percentile(vGD_bs,25),numpy.percentile(vGD_bs,75)-vGDm]
# 		vGDVm = numpy.percentile(vGDV_bs,50)
# 		vGDV_err = [vGDVm-numpy.percentile(vGDV_bs,25),numpy.percentile(vGDV_bs,75)-vGDVm]
# 		vGDVEm = numpy.percentile(vGDVE_bs,50)
# 		vGDVE_err = [vGDVEm-numpy.percentile(vGDVE_bs,25),numpy.percentile(vGDVE_bs,75)-vGDVEm]
# 		
	# vG_allenv.append(vGm)
# 	vG_allenv_err.append(vG_err)
# 	vGD_allenv.append(vGDm)
# 	vGD_allenv_err.append(vGD_err)
# 	vGDV_allenv.append(vGDVm)
# 	vGDV_allenv_err.append(vGDV_err)
# 	vGDVE_allenv.append(vGDVEm)
# 	vGDVE_allenv_err.append(vGDVE_err)
	
	vG_allenv.append(vG)
	vGD_allenv.append(vGD)
	vGDV_allenv.append(vGDV)
	vGDVE_allenv.append(vGDVE)

vG_allenv = numpy.array(vG_allenv)
#vG_allenv_err = numpy.array(vG_allenv_err)
vGD_allenv = numpy.array(vGD_allenv)
#vGD_allenv_err = numpy.array(vGD_allenv_err)
vGDV_allenv = numpy.array(vGDV_allenv)
#vGDV_allenv_err = numpy.array(vGDV_allenv_err)
vGDVE_allenv = numpy.array(vGDVE_allenv)
#vGDVE_allenv_err = numpy.array(vGDVE_allenv_err)
err_var_by_env = numpy.array(err_var_by_env)

pt.figure(figsize=(6,6))

width = .8
gbar = pt.bar(numpy.arange(8),vG_allenv,width,color='C1')
dbar = pt.bar(numpy.arange(8),vGD_allenv-vG_allenv,width,bottom=vG_allenv,color='C3')
vbar = pt.bar(numpy.arange(8),vGDV_allenv-vGD_allenv,width,bottom=vGD_allenv,color='C2')
ebar = pt.bar(numpy.arange(8),vGDVE_allenv-vGDV_allenv,width,bottom=vGDV_allenv,color='C0')


stochbar = pt.bar(numpy.arange(8),1-err_var_by_env-vGDVE_allenv,width,bottom=vGDVE_allenv,color='C4')
errbar = pt.bar(numpy.arange(8),err_var_by_env,width,bottom=1-err_var_by_env,color='C5')
ax = pt.gca()
ax.set_xticks(numpy.arange(8))
ax.set_xticklabels(envt_labels)
#ax.set_xlim(-1,10)
pt.legend([gbar[0],dbar[0],vbar[0],ebar[0],stochbar[0],errbar[0]],['Multihit genes','Hap/Dip','Virus','Home Env','Other\npleiotropy','Meas Err'],loc=[.65,.6],fontsize=9)
ax.set_ylabel('Proportion of variance',fontsize=12)
ax.set_xlabel('Measurement environment',fontsize=12)
pt.savefig('figures/fitness_predictors_nohome_genes_11_21_2019.pdf',bbox_inches='tight')

print('Multihit Gene variance ',vG_allenv)
print(numpy.mean(vG_allenv))
print('All genetic varaince ',vGDV_allenv)
print(numpy.mean(vGDV_allenv))
print('Additional env ',vGDVE_allenv-vGDV_allenv)
print(numpy.mean(vGDVE_allenv-vGDV_allenv))
print('Other stoch ',1-err_var_by_env-vGDVE_allenv)
print(numpy.mean(1-err_var_by_env-vGDVE_allenv))

file.write('Genes' + '\t' + ('\t').join([str(p) for p in vG_allenv]) + '\n')
file.write('Genes + Hap/Dip' + '\t' + ('\t').join([str(p) for p in vGD_allenv]) + '\n')
file.write('Genes + Hap/Dip + Virus' + '\t' + ('\t').join([str(p) for p in vGDV_allenv]) + '\n')
file.write('Genes + Hap/Dip + Virus + Evolution env.' + '\t' + ('\t').join([str(p) for p in vGDVE_allenv]) + '\n')
file.write('Evolution env. above genetics' + '\t' + ('\t').join([str(p) for p in vGDVE_allenv - vGDV_allenv]) + '\n')
file.write('Other pleiotropy' + '\t' + ('\t').join([str(p) for p in 1 - vGDVE_allenv - err_var_by_env]) + '\n')
file.write('Measurement err.' + '\t' + ('\t').join([str(p) for p in err_var_by_env]))

file.close()

####Same, restricted to V+ clones
#####

# vG_allenv = []
# vG_allenv_err = []
# vGD_allenv = []
# vGD_allenv_err = []
# vGDV_allenv = []
# vGDV_allenv_err = []
# vGDE_allenv = []
# vGDE_allenv_err = []
# 
# err_var_by_env = []
# 
# for i in range(len(evolution_environments)):
# 	
# 	meas_env = measurement_environments[i]
# 	predictor_mat_all = []
# 	
# 	outcome_list = []
# 	meas_err_list = []
# 	for j in range(len(evolution_environments)):
# 		if j != i:
# 			evol_env = evolution_environments[j]
# 			gen_norm_factor = gen_dict[evol_env]
# 			for clone in fit_dict[evol_env]:
# 				if fit_dict[evol_env][clone]['virus'] == 'Virus retained' and clone in mut_dict[evol_env]: #
# 					
# 					muts = mut_dict[evol_env][clone]
# 					
# 					clone_fit = fit_dict[evol_env][clone][meas_env]*700/gen_norm_factor
# 					clone_err = std_err_dict[evol_env][clone][meas_env]*700/gen_norm_factor
# 			
# 					outcome_list.append(clone_fit)
# 					meas_err_list.append(clone_err)
# 			
# 					predictor_list = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
# 					predictor_list[j] = 1. ###Indicator variable on the evolution environment
# 			
# 				
# 					#print(fit_dict[evol_env][clone]['hapdip'])
# 					
# 					
# 					#if fit_dict[evol_env][clone]['virus'] == 'Virus retained':
# 					#	v_indicator = 1.
# 						#print('found V')
# 					#else:
# 					#	v_indicator = 0.
# 					#print(fit_dict[evol_env][clone]['hapdip'])
# 					if fit_dict[evol_env][clone]['hapdip'] == 'Diploid':
# 						d_indicator = 1.
# 						#print('found D')
# 					else:
# 						d_indicator = 0.
# 			
# 					#predictor_list[8] = v_indicator ##Indicator variable on virus status
# 					predictor_list[8] = d_indicator ##Indicator variable on diploid status
# 				
# 					
# 					for m in muts:
# 						if m[4] in gene_family_lookup:
# 							cat = gene_family_lookup[m[4]]
# 							if cat != 'other':
# 								cat_index = family_inds[cat]
# 								predictor_list[9+cat_index] = 1
# 								
# 					predictor_mat_all.append(predictor_list)
# 	
# 	predictor_mat_all = numpy.array(predictor_mat_all)
# 	outcome_list = numpy.array(outcome_list)
# 	meas_err_list = numpy.array(meas_err_list)
# 	meas_err_var = numpy.mean(meas_err_list**2)/numpy.var(outcome_list)
# 	err_var_by_env.append(meas_err_var)
# 	
# 	nclones = predictor_mat_all.shape[0]
# 	vG_bs = []
# 	vGD_bs = []
# 	vGDV_bs = []
# 	vGDE_bs = []
# 	
# 	for n in range(nbootstrap):
# 		bs_inds = numpy.random.choice(numpy.arange(nclones),size=nclones)
# 		vGreg = LinearRegression().fit(predictor_mat_all[bs_inds,9:], outcome_list[bs_inds])
# 		vGDreg = LinearRegression().fit(predictor_mat_all[bs_inds,8:], outcome_list[bs_inds])
# 		#vGDVreg = LinearRegression().fit(predictor_mat_all[bs_inds,8:], outcome_list[bs_inds])
# 		vGDEreg = LinearRegression().fit(predictor_mat_all[bs_inds,:], outcome_list[bs_inds])
# 		
# 		
# 		vG = vGreg.score(predictor_mat_all[bs_inds,9:], outcome_list[bs_inds])
# 		vGD = vGDreg.score(predictor_mat_all[bs_inds,8:], outcome_list[bs_inds])
# 		#vGDV = vGDVreg.score(predictor_mat_all[bs_inds,8:], outcome_list[bs_inds])
# 		vGDE = vGDEreg.score(predictor_mat_all[bs_inds,:], outcome_list[bs_inds])
# 		
# 		vG_bs.append(vG)
# 		vGD_bs.append(vGD)
# 		#vGDV_bs.append(vGDV)
# 		vGDE_bs.append(vGDE)
# 		
# 		vGm = numpy.percentile(vG_bs,50)
# 		vG_err = [vGm-numpy.percentile(vG_bs,25),numpy.percentile(vG_bs,75)-vGm]
# 		vGDm = numpy.percentile(vGD_bs,50)
# 		vGD_err = [vGDm-numpy.percentile(vGD_bs,25),numpy.percentile(vGD_bs,75)-vGDm]
# 		#vGDVm = numpy.percentile(vGDV_bs,50)
# 		#vGDV_err = [vGDVm-numpy.percentile(vGDV_bs,25),numpy.percentile(vGDV_bs,75)-vGDVm]
# 		vGDEm = numpy.percentile(vGDE_bs,50)
# 		vGDE_err = [vGDEm-numpy.percentile(vGDE_bs,25),numpy.percentile(vGDE_bs,75)-vGDEm]
# 		
# 	vG_allenv.append(vGm)
# 	vG_allenv_err.append(vG_err)
# 	vGD_allenv.append(vGDm)
# 	vGD_allenv_err.append(vGD_err)
# 	#vGDV_allenv.append(vGDVm)
# 	#vGDV_allenv_err.append(vGDV_err)
# 	vGDE_allenv.append(vGDEm)
# 	vGDE_allenv_err.append(vGDE_err)
# 	
# 
# vG_allenv = numpy.array(vG_allenv)
# vG_allenv_err = numpy.array(vG_allenv_err)
# vGD_allenv = numpy.array(vGD_allenv)
# vGD_allenv_err = numpy.array(vGD_allenv_err)
# #vGDV_allenv = numpy.array(vGDV_allenv)
# #vGDV_allenv_err = numpy.array(vGDV_allenv_err)
# vGDE_allenv = numpy.array(vGDE_allenv)
# vGDE_allenv_err = numpy.array(vGDE_allenv_err)
# err_var_by_env = numpy.array(err_var_by_env)
# 
# pt.figure(figsize=(6,5))
# 
# width = .8
# gbar = pt.bar(numpy.arange(8),vG_allenv,width,color='C1')
# dbar = pt.bar(numpy.arange(8),vGD_allenv-vG_allenv,width,bottom=vG_allenv,color='C3')
# #vbar = pt.bar(numpy.arange(8),vGDV_allenv-vGD_allenv,width,bottom=vGD_allenv,yerr=vGDV_allenv_err.T)
# ebar = pt.bar(numpy.arange(8),vGDE_allenv-vGD_allenv,width,bottom=vGD_allenv,color='C0')
# 
# 
# stochbar = pt.bar(numpy.arange(8),1-err_var_by_env-vGDE_allenv,width,bottom=vGDE_allenv,color='C4')
# errbar = pt.bar(numpy.arange(8),err_var_by_env,width,bottom=1-err_var_by_env,color='C5')
# ax = pt.gca()
# ax.set_xticks(numpy.arange(8))
# ax.set_xticklabels(envt_labels)
# #ax.set_xlim(-1,10)
# pt.legend([gbar[0],dbar[0],ebar,stochbar[0],errbar[0]],['Multihit genes','Hap/Dip','Home Env','Stochasticity','Meas Err'],loc=[.65,.6],fontsize=9)
# ax.set_ylabel('Proportion of variance',fontsize=12)
# ax.set_xlabel('Measurement environment',fontsize=12)
# pt.savefig('figures/fitness_predictors_nohome_genes_vplusonly.pdf',bbox_inches='tight')
