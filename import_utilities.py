import numpy

def import_fitness_data(file_in = 'data/pleiotropy_fitness_data_7_4_2019.txt'):
	
	##Import fitness data

	file = open(file_in,'r')

	fit_dict = {}
	err_dict = {}
	gen_dict = {}

	fitness_measurement_columns = numpy.arange(2,32,2)
	replicate_measurement_columns = numpy.arange(3,33,2)

	firstline = True
	for line in file:
		if firstline:
			header = line.strip().split('\t')
			firstline = False
		else:
			line_list = line.strip().split('\t')
		
			evol_env = line_list[0]
			clone = line_list[1]
			virus_retained = line_list[-3]
			num_gens = float(line_list[-2])
			hap_dip = line_list[-1]
		
			if evol_env not in fit_dict:
				fit_dict[evol_env] = {}
				err_dict[evol_env] = {}
				gen_dict[evol_env] = num_gens
			
			if clone not in fit_dict[evol_env]:
				fit_dict[evol_env][clone] = {}
				err_dict[evol_env][clone] = {}
				fit_dict[evol_env][clone]["virus"] = virus_retained
				fit_dict[evol_env][clone]["hapdip"] = hap_dip
		
			for column in fitness_measurement_columns:
				if line_list[column] != "NA":
					fit_dict[evol_env][clone][header[column]] = float(line_list[column])
			
			for column in replicate_measurement_columns:
				if line_list[column] != "NA":
					err_dict[evol_env][clone][header[column-1]] = [float(m) for m in line_list[column].split(',')]

	###Find standard errors for each measured fitness; resample fitnesses

	std_err_dict = {}
	random_reps = {}
	for evol_env in err_dict:
		std_err_dict[evol_env] = {}
		random_reps[evol_env] = {}
		for clone in err_dict[evol_env]:
			std_err_dict[evol_env][clone] = {}
			random_reps[evol_env][clone] = {}
			for meas_env in err_dict[evol_env][clone]:
				deviations = numpy.array(err_dict[evol_env][clone][meas_env])
			
				std_err = numpy.sqrt( numpy.sum(deviations**2)/(float(len(deviations - 1))) )/numpy.sqrt(len(deviations)) #Unbiased variance estimate divided by sqrt of number of replicates
				std_err_dict[evol_env][clone][meas_env] = std_err
			
				random_reps[evol_env][clone][meas_env] = [err_dict[evol_env][clone][meas_env][0] + fit_dict[evol_env][clone][meas_env], err_dict[evol_env][clone][meas_env][1] + fit_dict[evol_env][clone][meas_env] ] ###Replicate 1 vs. replicate 2
				
	return fit_dict, gen_dict, std_err_dict, random_reps

def import_mutation_data(file_in='data/mutation_table_7_3_2019.txt'):

	###Import mutation data

	file = open(file_in, 'r')

	mut_dict = {}

	firstline = True

	for line in file:

		if firstline:
			firstline = False
		else:
			line_list = line.strip().split('\t')
			env = line_list[0]
			clone = line_list[1]
		
			if env not in mut_dict:
				mut_dict[env] = {}
				mut_dict[env][clone] = [ line_list[2:] ]
			elif clone not in mut_dict[env]:
				mut_dict[env][clone] = [ line_list[2:] ]
			else:
				mut_dict[env][clone].append( line_list[2:] )

	file.close()



	###Determine list of potential SGV mutations, which are identical at the nucleotide level

	mut_list = []
	all_clones_list = []
	sgv_list1 = []
	sgv_list2 = []
	sgv_list_envs_clones = []

	for env in mut_dict:
	
		for clone in mut_dict[env]:
		
			for mut in mut_dict[env][clone]:
			
				if mut in mut_list:
				
					ind = mut_list.index(mut)
				
					prev_clone = all_clones_list[ind]
					sgv_list1.append(mut)
					sgv_list_envs_clones.append(prev_clone)
					sgv_list1.append(mut)
					sgv_list_envs_clones.append((env,clone))
				
				else:
				
					mut_list.append(mut)
					all_clones_list.append((env,clone))

	sgv_list_non = []
	env_clone_list_non = []
	for i in range(len(sgv_list1)):
		entry = sgv_list1[i]
		id = sgv_list_envs_clones[i]

		if entry[6] == 'Non':
			sgv_list_non.append(entry[4])
			env_clone_list_non.append(id)

	###Count nonsynonymous mutations in genes by environment. Note that we are only counting each gene once for each clone, as most instances where there are apparently 2 nearby mutations are likely to be a compound event.

	gene_dict = {}
	clone_mut_num_dict = {}
	clones_per_env_dict = {}
	for env in mut_dict:
		clones_per_env_dict[env] = 0
		clone_mut_num_dict[env] = {}
		for clone in mut_dict[env]:

			local_gene_list = []
			clone_mut_num_dict[env][clone] = 0
			clones_per_env_dict[env] += 1

			for mut in mut_dict[env][clone]:
	
				if mut not in sgv_list1:
		
					annot = mut[6]
					gene = mut[4]
		
					if (annot == 'Non' and gene not in local_gene_list):
			
						local_gene_list.append(gene)
			
						if ('*' in mut or 'FS' in mut):
				
							tally = numpy.array([0,1]) ##nonsense mutation
						else:
							tally = numpy.array([1,0]) ##missense mutation
				
			
						if gene not in gene_dict:
				
							gene_dict[gene] = {}
							gene_dict[gene][env] = tally
			
						elif env not in gene_dict[gene]:
				
							gene_dict[gene][env] = tally
			
						else:
				
							gene_dict[gene][env] += tally
				
						clone_mut_num_dict[env][clone] += 1

	return gene_dict, clone_mut_num_dict, clones_per_env_dict, mut_dict

def import_GO_slim(file_in = 'data/GO_Slim_output_7_4_2019.txt'):

	#####Read in the GO slim analysis and assign genes to families

	go_file = open(file_in,'r')

	go_slim_data = {}

	family_list_go = []
	enrichment_list = []

	firstline = True
	for line in go_file:
		if firstline:
			firstline = False
			continue
		entry = line.strip().split('\t')
		family = entry[1]
		genes = entry[8]
		percent_in_list = float(entry[2])/float(entry[3])
		percent_in_genome = float(entry[5])/float(entry[6])
		go_slim_data[family] = {}
		go_slim_data[family]['genes'] = genes.split(', ')
		go_slim_data[family]['enrichment'] = percent_in_list/percent_in_genome
	
		family_list_go.append(family)
		enrichment_list.append(percent_in_list/percent_in_genome)

	family_list_go_ordered = numpy.array(family_list_go)[numpy.argsort(enrichment_list)[::-1]]

	go_file.close()
	
	return family_list_go_ordered, go_slim_data

def import_cured_ref_FA(file_in = 'data/Cured_ref_21C_fitness_table.txt'):
	
	
	file = open(file_in, 'r')

	fit_dict_cured_ref_fa = {}

	firstline = True
	for line in file:
		if firstline:
			firstline = False
		else:
			line_list = line.strip().split('\t')
			evol_env = line_list[0]
			clone = line_list[1]
		
			if evol_env not in fit_dict_cured_ref_fa:
				fit_dict_cured_ref_fa[evol_env] = {}
				fit_dict_cured_ref_fa[evol_env] = {}
				fit_dict_cured_ref_fa[evol_env][clone] = {}
				fit_dict_cured_ref_fa[evol_env][clone]['virusplusref'] = float(line_list[2])
				fit_dict_cured_ref_fa[evol_env][clone]['virusminusref'] = float(line_list[3])
			else:
				fit_dict_cured_ref_fa[evol_env][clone] = {}
				fit_dict_cured_ref_fa[evol_env][clone]['virusplusref'] = float(line_list[2])
				fit_dict_cured_ref_fa[evol_env][clone]['virusminusref'] = float(line_list[3])

	file.close()
	
	return fit_dict_cured_ref_fa