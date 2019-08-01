import numpy
import matplotlib.pylab as pt
import import_utilities

####Import clone fitness data

fit_dict, gen_dict, std_err_dict, random_reps = import_utilities.import_fitness_data('data/pleiotropy_fitness_data_7_4_2019.txt')

####Import recon fitness data

gene_list = ['PEP1','SIR2','ASG1','HAL5','STP3','CUE4']

file_in = open('data/reconstructions_fitness_table_Salt17.txt', 'r')

recon_fit_mat = {}

firstline = True

for line in file_in:
	line_list = line.strip().split('\t')
	
	if firstline:
		header = line_list
		firstline = False
	else:
		gene = line_list[0]
		
		recon_type = line_list[1]
		if gene in gene_list:
			if gene not in recon_fit_mat:
				recon_fit_mat[gene] = {}
			recon_fit_mat[gene][recon_type] = {}
			
			for i in range(4):
				
				column = 2*(i+1)
				key = header[column]
				
				recon_fit_mat[gene][recon_type][key] = float(line_list[column])

file_in.close()

####

measurement_envs = ["fitness in SC","fitness in 21 C","fitness in 37 C","fitness in .8 M"]
measurement_envs_r = ["Fitness in SC","Fitness in 21 C","Fitness in 37 C","Fitness in .8 M"]
clone_fit_vec = []
clone_err_vec = []

for m in measurement_envs:
	
	evol_env = "Evolved in .8 M"
	clone = "clone17"
	
	clone_fit_vec.append( fit_dict[evol_env][clone][m] )
	clone_err_vec.append( std_err_dict[evol_env][clone][m] )
	
#####

fig = pt.figure()
ax = pt.gca()

pt.figure()
aux = pt.gca()

ax.errorbar(numpy.arange(4), numpy.array(clone_fit_vec)*100, yerr=clone_err_vec, marker='o', color = 'k')
legend_handles = []
legend_handles.append(aux.scatter(numpy.arange(4), numpy.array(clone_fit_vec)*100, color = 'k'))

symbol_key = ['-','--']
color_key = dict(zip(gene_list, ['MediumSlateBlue','darkblue','Tomato','Brown','darkgreen','indigo']))
genes_already_in_legend = []
for gene in gene_list:
	for recon_type in recon_fit_mat[gene]:
		if gene == 'CUE4' and '1' in recon_type:
			continue
		elif 'KO' in recon_type:
			continue
			symbol = symbol_key[1]
		else:
			symbol = symbol_key[0]
		fit_list = []
		for m in measurement_envs_r:
			fit_list.append( recon_fit_mat[gene][recon_type][m] )
		
		
		if gene not in genes_already_in_legend:
			legend_handles.append(aux.scatter(numpy.arange(4), numpy.array(fit_list)*100, color = color_key[gene]))
			genes_already_in_legend.append(gene)
			ax.plot(numpy.arange(4), numpy.array(fit_list)*100, marker='o', linestyle=symbol, color = color_key[gene])
			
ax.set_xticks(numpy.arange(4))
ax.set_xticklabels(['SC','Low Temp','High Temp','High Salt'])
pt.figure(fig.number)
pt.ylabel('Relative fitness (%)')
pt.xlabel('Measurement environment')
pt.legend(legend_handles, ['High Salt Clone 17','PEP1','SIR2','ASG1','HAL5','STP3','CUE4'])
pt.savefig('figures/salt17_recons_check.pdf',bbox_inches='tight')

		
		