import matplotlib.pylab as pt
import numpy

file = open('data/frequency_dependence_vplus_vs_vminus_21C.txt','r')
firstline = True

fit_list = []
freq_list = []
err_list = []
for line in file:
	
	if firstline:
		firstline = False
	else:
		linelist = line.strip().split('\t')
		if linelist[0] == 'Measured in 21 C' and linelist[1] == 'cit2':
			freq_list.append(float(linelist[2]))
			fit_list.append(float(linelist[3]))
			err_list.append(float(linelist[4])) ###This is the difference in fitness between the two replicate measurements

freq_list = numpy.array(freq_list)
fit_list = numpy.array(fit_list)
err_list = numpy.array(err_list)

pt.figure()
pt.errorbar(freq_list,-1*fit_list,yerr=err_list/2., marker='o') ###Fitness of V- relative to V+ is -1*fitness as listed in file
pt.xlabel('Frequency of V- clone')
pt.ylabel('Fitness of V- clone')
pt.savefig('figures/frequency_dependence_vplus_vminus_21C.pdf',bbox_inches='tight')