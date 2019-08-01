import matplotlib.pylab as pt
import numpy

def get_column( header, field_name ):
	
	header = numpy.array(header)
	
	col = numpy.where( header == field_name )[0][0]
	
	return col
	
file_in = 'data/Backcross_fitnesses.txt'

file = open(file_in, 'r')

firstline = True

mating_dict = {}

for line in file:
	
	if firstline:
	
		header = line.strip().split('\t')
		
		mat_label_col = get_column( header, 'Cross Number')
		
		parent1_col = get_column( header, 'Parent 1 fitness, 21 C')
		print(parent1_col)
		parent1_rep_col = get_column( header, 'Parent 1 replicate deviation')
		parent2_col = get_column( header, 'Parent 2 fitness, 21 C')
		parent2_rep_col = get_column( header, 'Parent 2 replicate deviation')
		
		spore_fit_col = get_column( header, 'Spore fitness in 21 C')
		spore_rep_col = get_column( header, 'replicate deviation')
		
		firstline = False
	
	else:
		
		line_list = line.strip().split('\t')
		
		mating = line_list[mat_label_col]
		
		if mating not in mating_dict:
			
			mating_dict[mating] = {}
			mating_dict[mating]['spores'] = {}
			mating_dict[mating]['parents'] = {}
			mating_dict[mating]['spores']['fitnesses'] = [ float(line_list[spore_fit_col]) ]
			mating_dict[mating]['spores']['deviations'] = [ float(line_list[spore_rep_col]) ]
			mating_dict[mating]['parents']['fitnesses'] = [ float(line_list[parent1_col]), float(line_list[parent2_col]) ]
			mating_dict[mating]['parents']['deviations'] = [ float(line_list[parent1_rep_col]), float(line_list[parent2_rep_col]) ]
		
		else:
			
			mating_dict[mating]['spores']['fitnesses'].append( float(line_list[spore_fit_col]) )
			mating_dict[mating]['spores']['deviations'].append( float(line_list[spore_rep_col]) )

file.close()

####

fig, (ax1, ax2) = pt.subplots(1,2, figsize=(4,4))

mating = 'mat5'

for i in range(len(mating_dict[mating]['spores']['fitnesses'])):
	
	fitness = mating_dict[mating]['spores']['fitnesses'][i]
	err = mating_dict[mating]['spores']['deviations'][i]

	r = numpy.random.random()
	
	ax2.errorbar( r, fitness, yerr = err, marker = 'o', capsize=0, color = 'C0', label='Spore', alpha = .8)#, markeredgewidth=0, markersize=4, linewidth=.5)

ax2.errorbar( [.5, .5], mating_dict[mating]['parents']['fitnesses'], yerr = mating_dict[mating]['parents']['deviations'], fmt = 'o', capsize=0, color = 'Black', label='Parent')#, markeredgewidth=0, markersize=4, linewidth=.5)

ax2.set_ylim(-.5,.05)

mating = 'mat1'


for mating in ['mat1','mat8','mat11']:
	
	for i in range(len(mating_dict[mating]['spores']['fitnesses'])):
	
		fitness = mating_dict[mating]['spores']['fitnesses'][i]
		err = mating_dict[mating]['spores']['deviations'][i]
		r = numpy.random.random()
		ax1.errorbar( r, fitness, yerr = err, marker = 'o', capsize=0, label='Spore',color = 'C0', alpha = .8)#, markeredgewidth=0, markersize=4, linewidth=.5)
ax1.errorbar( [.5, .5], mating_dict[mating]['parents']['fitnesses'], yerr = mating_dict[mating]['parents']['deviations'], fmt = 'o', capsize=0, color = 'Black', label='Parent')#, markeredgewidth=0, markersize=4, linewidth=.5)

ax1.set_ylim(-.5,.05)
ax1.set_ylabel('Fitness')
ax1.set_title('Backcross')
ax2.set_title('Control cross')
ax1.set_xticks([])
ax2.set_xticks([])
ax2.legend()
pt.tight_layout()
pt.savefig('figures/backcross_figure.pdf',bbox_inches='tight')
	
	
	
	