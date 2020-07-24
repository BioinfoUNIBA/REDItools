#!/usr/bin/env python
#################################### REDI OUT TABLE ########################################################
#Region		Position	Reference	Strand	Coverage-q30	MeanQ	BaseCount[A,C,G,T]	   #
#AllSubs	Frequency	gCoverage-q30	gMeanQ	gBaseCount[A,C,G,T]	gAllSubs	gFrequency #
############################################################################################################

###################################GET_DE_events_table###########################################################
#chromosome	position	type_editing	SRR3306830_CTRL		SRR3306831_CTRL		SRR3306832_CTRL #	
#SRR3306833_CTRL	SRR3306834_CTRL	SRR3306835_CTRL   SRR3306836_CTRL	SRR3306823_DIS	SRR3306824_DIS	#
#SRR3306825_DIS	SRR3306826_DIS	SRR3306827_DIS	SRR3306828_DIS	SRR3306829_DIS	[num_controls/num_disease]	#
#delta_diff	pvalue (Mannwhitney)       									#
#################################################################################################################

import os, sys, argparse
from scipy import stats
from scipy.stats import wilcoxon, mannwhitneyu, fisher_exact
import numpy as np
import pandas as pd
import math

parser = argparse.ArgumentParser()
parser.add_argument("-c", action = 'store', dest = 'min_coverage', 
		type = int, default=10,  help='Coverage-q30')
parser.add_argument("-cpval", action = 'store', dest = 'pvalue_correction',
                type = int, default = 0, help = '1 --> Bonferroni correction / 2 --> Benjamini hochberg')
parser.add_argument("-input_file", action = 'store', dest = 'samples_informations_file',
		type = str, default= 'empty', help = 'Comma separated file  e.g: SRR3306830,Control \
		SRR3306829,Healthy...etc')
parser.add_argument("-f", action = 'store', dest = 'min_edit_frequency',
		type = float, default=0.1, help='Editing Frequency')
parser.add_argument("-mts", action = 'store', dest = 'min_sample_testing',
		type = float, default=50.0, help="min percentage of each sample category")
parser.add_argument("-sig", action = 'store', dest = 'only_significant',
		type = str, default = 'no', help = 'Return only significant editing events')
parser.add_argument("-linear", action = 'store_true', help = 'Enable linear model')

args = parser.parse_args()
min_coverage = args.min_coverage
min_edit_frequency = args.min_edit_frequency
min_sample_testing = args.min_sample_testing
only_significants = args.only_significant
pvalue_correction = args.pvalue_correction
samples_informations_file = args.samples_informations_file
enable_linear_model = args.linear

if args.samples_informations_file == 'empty':
	parser.error('sample_informations_file is MISSING!' + '\n' + \
	'Please type "python get_DE_events.py -h" for more details on usage of this script')


def call_differential_editing_sites(config_file):
	stability_value = 0.03 #value below which you may use a lower coverage for adding more samples to increase power
	min_disease_people = 5 #min number people supporting higher coverage for whch you may base stability off measurements off of
	min_control_people = 5  #min number control poeple supporting higher coverage for which you may base stability off of
	min_disease_people_5_cov = 10 #min disease number of people of 5 coverage you must have if needing to use unstable 5x coverage
	min_control_people_5_cov = 10 #min control number of people of 5 coverage you must have if needing to use unstable 5x coverage
	editing_file= './temp.csv'
	output_file = './editing_sites.with_stats_converted_disease.csv'
	#read in files
	editing_table = pd.read_csv(editing_file,sep='\t')
	#config_table = pd.read_csv(config_file,sep=',',header=None)
	config_table = pd.read_csv(config_file,sep=',',skiprows=1,header=None)
	all_people = config_table[0]
	disease_people = config_table[0][config_table[1] == "DIS"].reset_index(drop = True) #TODO Change do disease!!!
	control_people = config_table[0][config_table[1] == "CTRL"].reset_index(drop = True) #TODO Change to control!!!

	#now get just an editing table and coverage table
	edit_level_table = editing_table[all_people]
	#edit_level_table = editing_table[np.r_[all_people]]

	def get_editing_levels_for_cov_table(i):
	  info = i.astype(str).str.split(pat="\\^")
	  editing_levels = info.apply(lambda x: float('nan') if x[0] == "nan" else x[2])
	  return editing_levels
	cov_table = edit_level_table.apply(get_editing_levels_for_cov_table)
	cov_table = cov_table.apply(lambda x: pd.to_numeric(x)) #TODO check if as.numeric and pandas to_numeric do the same.

	def get_editing_levels(i):
	  info = i.astype(str).str.split(pat="\\^")
	  editing_levels = info.apply(lambda x: float('nan') if x[0] == "nan" else x[0])
	  return editing_levels
	edit_level_table = edit_level_table.apply(get_editing_levels)
	edit_level_table = edit_level_table.apply(lambda x: pd.to_numeric(x)) #TODO check precision on R and python

	#go down line by line and get the prevalence info and mean editing levels based off of stable coverages
	#WARNING I'm using float here, not integer allowing NaN values. Is ok?
	coverage_threshold_used = np.repeat(0.,edit_level_table.shape[0]) #will hold the coverage threshold required for this editing site
	stability_based_on = np.repeat(0.,edit_level_table.shape[0]) #will hold what coverage stability requirements were determined
	stable_mean_disease_editing_level = np.repeat(0.,edit_level_table.shape[0]) #mean autistic editing level using individuals passing coverage threshold
	stable_std_dev_disease_editing_level = np.repeat(0.,edit_level_table.shape[0]) #standard deviation of autistic editing level using individuals passing coverage threshold
	stable_mean_control_editing_level = np.repeat(0.,edit_level_table.shape[0]) #mean control editing level using individuals passing coverage threshold
	stable_std_dev_control_editing_level = np.repeat(0.,edit_level_table.shape[0]) #standard deviation of control editing level using individuals passing coverage threshold
	stable_number_disease_with_at_least_min_coverage = np.repeat(0.,edit_level_table.shape[0]) #number of autistic individuals passing the coverage threshold
	stable_number_disease_nonzero_editing_and_min_coverage = np.repeat(0.,edit_level_table.shape[0]) #number of autistic individuals without non zero editing level and passing coverage threshold
	stable_disease_prevalence = np.repeat(0.,edit_level_table.shape[0]) #proportion autistic individuals with nonzero editing
	stable_number_control_with_at_least_min_coverage = np.repeat(0.,edit_level_table.shape[0]) #same as disease but for control subjects
	stable_number_control_nonzero_editing_and_min_coverage = np.repeat(0.,edit_level_table.shape[0])
	stable_control_prevalence = np.repeat(0.,edit_level_table.shape[0])
	stable_total_number_individuals_nonzero_editing_and_min_coverage = np.repeat(0.,edit_level_table.shape[0]) #total number of disease and control subjects passing the coverage threshold and having nonzero editing level
	stable_mann_whitney_p_value = np.repeat(0.,edit_level_table.shape[0]) #wilcoxon rank sum test p value using individuals passing the coverage threshold
	stable_editing_level_effect_size = np.repeat(0.,edit_level_table.shape[0]) #difference between mean disease and mean control
	stable_frequency_fishers_p_value = np.repeat(0.,edit_level_table.shape[0]) #prevalence p value determined using two-tailed fisher's exact test
	stable_frequency_OR = np.repeat(0.,edit_level_table.shape[0]) #odds ratio of the fisher's exact teest
	stable_prevalence_effect_size = np.repeat(0.,edit_level_table.shape[0]) #difference in editing level prevalences between disease and control subjects
	#WARNING those are np arrays.

	for i in range(0,edit_level_table.shape[0]):
	  print i  #keep track of progress
	  disease_edit_row = edit_level_table.loc[i, disease_people]
	  control_edit_row = edit_level_table.loc[i, control_people]
	  disease_cov_row = cov_table.loc[i, disease_people]
	  control_cov_row = cov_table.loc[i, control_people]
	  #find what coverage we can base stability off of
	  number_disease_20_cov = disease_cov_row[disease_cov_row >= 20].count()
	  number_control_20_cov = control_cov_row[control_cov_row >=20].count()
	  number_disease_15_cov = disease_cov_row[disease_cov_row >= 15].count()
	  number_control_15_cov = control_cov_row[control_cov_row >= 15].count()
	  number_disease_10_cov = disease_cov_row[disease_cov_row >= 10].count()
	  number_control_10_cov = control_cov_row[control_cov_row >= 10].count()
	  number_disease_5_cov = disease_cov_row[disease_cov_row >= 5].count()
	  number_control_5_cov = control_cov_row[control_cov_row >= 5].count()
	  if number_disease_20_cov >= min_disease_people and number_control_20_cov >= min_control_people:
		stability_based_on[i] = 20
	  elif number_disease_15_cov >= min_disease_people and number_control_15_cov >= min_control_people:
		stability_based_on[i] = 15
	  elif number_disease_10_cov >= min_disease_people and number_control_10_cov >= min_control_people:
		stability_based_on[i] = 10
	  elif number_disease_5_cov >= min_disease_people_5_cov and number_control_5_cov >= min_control_people_5_cov:
		stability_based_on[i] = 5
	  else:
		#stability_based_on[i] = -99999 # there's no np.nan integer representation, only float. We use an invalid value.
		stability_based_on[i] = float('nan')

	  #need to deal with cases where there just are not enough disease individuals or control individuals to calculate mean
	  if np.isnan(stability_based_on[i]):

		coverage_threshold_used[i] = 5 #I warn users not to use editing sites that don't have any stability_based_on measurement. We include min coverage of 5 just to get statistical information anyways
		#stable_min_cov=5
		#otherwise we can now try to find the stable_min_cov that'll be used for calculation of all statistics'

	  else:
		current_stability_cov =  stability_based_on[i]
		stability_disease_mean = disease_edit_row[disease_cov_row >= current_stability_cov].mean()
		stability_control_mean = control_edit_row[control_cov_row >= current_stability_cov].mean()
		#print np.arange(5,stability_based_on[i]+1e-4,5)
		for j in np.arange(5,stability_based_on[i]+1e-4,5): #WARNING using 1e-4 allowing to include stop
		  disease_mean = disease_edit_row[disease_cov_row >= j].mean()
		  control_mean = control_edit_row[control_cov_row >= j].mean()
		  if np.absolute(disease_mean-stability_disease_mean) <=stability_value and np.absolute(control_mean-stability_control_mean) <=stability_value :
		    coverage_threshold_used[i] = j
		    break
	  #now let's calculate all our statics based on the stable coverage threshold
	  stable_min_cov = coverage_threshold_used[i]
	  disease_adju_edit_row = disease_edit_row[np.logical_and(np.logical_and((~np.isnan(disease_edit_row)), (~np.isnan(disease_cov_row))), (disease_cov_row >= stable_min_cov))]
	  disease_adju_cov_row = disease_cov_row[np.logical_and((~np.isnan(disease_cov_row)), (disease_cov_row >= stable_min_cov))]
	  control_adju_edit_row = control_edit_row[ np.logical_and(np.logical_and((~np.isnan(control_edit_row)), (~np.isnan(control_cov_row))), (control_cov_row >= stable_min_cov))]
	  control_adju_cov_row = control_cov_row[np.logical_and((~np.isnan(control_cov_row)), (control_cov_row >= stable_min_cov))]
	  stable_mean_disease_editing_level[i] = disease_adju_edit_row.mean()
	  stable_std_dev_disease_editing_level[i] = disease_adju_edit_row.std()
	  stable_mean_control_editing_level[i] = control_adju_edit_row.mean()
	  stable_std_dev_control_editing_level[i] = control_adju_edit_row.std()
	  stable_number_disease_with_at_least_min_coverage[i] = disease_adju_cov_row[disease_adju_cov_row >=stable_min_cov].count()
	  stable_number_disease_nonzero_editing_and_min_coverage[i] = disease_adju_cov_row[ (~np.isnan(disease_adju_cov_row)) & (disease_adju_cov_row >= stable_min_cov) & (disease_adju_edit_row > 0) ].count()
	  stable_disease_prevalence[i] = stable_number_disease_nonzero_editing_and_min_coverage[i]/stable_number_disease_with_at_least_min_coverage[i]
	  stable_number_control_with_at_least_min_coverage[i] = control_adju_cov_row[control_adju_cov_row >=stable_min_cov].count()
	  stable_number_control_nonzero_editing_and_min_coverage[i] = control_adju_cov_row[(~np.isnan(control_adju_cov_row)) & (control_adju_cov_row >= stable_min_cov) & (control_adju_edit_row > 0)].count()
	  stable_control_prevalence[i] = stable_number_control_nonzero_editing_and_min_coverage[i]/stable_number_control_with_at_least_min_coverage[i]
	  stable_total_number_individuals_nonzero_editing_and_min_coverage[i] = (stable_number_disease_nonzero_editing_and_min_coverage[i] + stable_number_control_nonzero_editing_and_min_coverage[i]).sum()
	  if (len(disease_adju_edit_row) >=1) & (len(control_adju_edit_row) >=1):
		if (np.all(disease_adju_edit_row.values == control_adju_edit_row.values)):
		  stable_mann_whitney_p_value[i] = float('nan')
		else:
		  temp, stable_mann_whitney_p_value[i] = mannwhitneyu(disease_adju_edit_row,control_adju_edit_row, alternative='two-sided')
	  else:
		stable_mann_whitney_p_value[i] = float('nan')
	  stable_editing_level_effect_size[i] =  np.absolute(stable_mean_disease_editing_level[i] - stable_mean_control_editing_level[i])
	  fisher_matrix = np.matrix([[stable_number_disease_nonzero_editing_and_min_coverage[i], stable_number_disease_with_at_least_min_coverage[i]-stable_number_disease_nonzero_editing_and_min_coverage[i]], [stable_number_control_nonzero_editing_and_min_coverage[i], stable_number_control_with_at_least_min_coverage[i]-stable_number_control_nonzero_editing_and_min_coverage[i]]])
	  stable_frequency_OR[i], stable_frequency_fishers_p_value[i] = fisher_exact(fisher_matrix)  
	  #print stable_frequency_OR[i]
	  #print stable_frequency_fishers_p_value[i]
	  stable_prevalence_effect_size[i] = np.absolute(stable_disease_prevalence[i] - stable_control_prevalence[i])

	#now put everything back together as a table
	header_info = editing_table[['chromosome','position','type_editing']]
	stats_table = pd.DataFrame(coverage_threshold_used)
	stats_table = stats_table.rename(columns={stats_table.columns[0]: 'coverage_threshold_used'})
	stats_table['stability_based_on'] = pd.DataFrame(stability_based_on)
	stats_table['stable_mean_disease_editing_level'] = pd.DataFrame(stable_mean_disease_editing_level)
	stats_table['stable_std_dev_disease_editing_level'] = pd.DataFrame(stable_std_dev_disease_editing_level)
	stats_table['stable_mean_control_editing_level'] = pd.DataFrame(stable_mean_control_editing_level)
	stats_table['stable_std_dev_control_editing_level'] = pd.DataFrame(stable_std_dev_control_editing_level)
	stats_table['stable_number_disease_with_at_least_min_coverage'] = pd.DataFrame(stable_number_disease_with_at_least_min_coverage)
	stats_table['stable_number_disease_nonzero_editing_and_min_coverage'] = pd.DataFrame(stable_number_disease_nonzero_editing_and_min_coverage)
	stats_table['stable_disease_prevalence'] = pd.DataFrame(stable_disease_prevalence)
	stats_table['stable_number_control_with_at_least_min_coverage'] = pd.DataFrame(stable_number_control_with_at_least_min_coverage)
	stats_table['stable_number_control_nonzero_editing_and_min_coverage'] = pd.DataFrame(stable_number_control_nonzero_editing_and_min_coverage)
	stats_table['stable_control_prevalence'] = pd.DataFrame(stable_control_prevalence)
	stats_table['stable_total_number_individuals_nonzero_editing_and_min_coverage'] = pd.DataFrame(stable_total_number_individuals_nonzero_editing_and_min_coverage)
	stats_table['stable_mann_whitney_p_value'] = pd.DataFrame(stable_mann_whitney_p_value)
	stats_table['stable_editing_level_effect_size'] = pd.DataFrame(stable_editing_level_effect_size)
	stats_table['stable_frequency_fishers_p_value'] = pd.DataFrame(stable_frequency_fishers_p_value)
	stats_table['stable_frequency_OR'] = pd.DataFrame(stable_frequency_OR)
	stats_table['stable_prevalence_effect_size'] = pd.DataFrame(stable_prevalence_effect_size)

	full_table = pd.concat([header_info, stats_table, editing_table[all_people]], axis=1)

	#write the full_table to output
	full_table.to_csv(output_file, sep='\t', index=False)

	print "job completed\n"




def Set_Chr_Nr(Chr):
    """ Sort by chromosome """
    if Chr: 
	New = Chr.lstrip('chr').split('_')[0]
        if New == 'X': New = 23
        elif New == 'Y': New = 24
        elif New == 'M': New = 25
        else: New = int(New)
    else:
        New = 0
    return New

def Sample_percentage(row):
	"""Percentage of samples from each type"""
	percentage = (len(filter(lambda x: x!= '-', row))/float(len(row)))*100
	return round(percentage)

def Sample_count(row):
	"""Number of samples from each type"""
	count = len(filter(lambda x: x!= '-', row))
	return count 

def get_bh(pvalue,siglevel):
	"""B-H correction """
	pvalue.sort()
	x=1
	y=0
	p=0
	for i in pvalue:
		nf=i[0]*len(pvalue)
		fdr=nf/x
		if fdr<=siglevel:
			i[1].append('True')
			p=i[0]
			y+=1
		else: i[1].append('False')
		x+=1
	return pvalue,y,p
	
def get_b(pvalue,siglevel): 
	"""Bonferroni correction"""
	pvalue.sort()
	y=0
	pp=1.0
	for i in pvalue:
		p=i[0]*len(pvalue)
		if p<=siglevel:
			i[1].append('True')
			y+=1
			if p<pp: pp=p
		else: i[1].append('False')
	return pvalue,y,pp

def only_sig(row_a,row):
	"""Returns only significant events"""
	if(row_a[-1] != '-' and row_a[-1] != 0.0 and row_a[-1] <= 0.05):
		row =  row[0].split('_') + row[2:]
		row.insert(2, 'A.to.G')
		print '\t'.join(map(str,row))

def tuple_replace(i):
	if type(i) == tuple:
		return i[0]
	else:
		return i

def tuple_replace_bis(k):
        if type(k) == tuple:
                return k[1]
        else:
             	return k

def remove_underscore(lis):
	lis = lis[:lis.index('_')]
	return lis


sample_informations = {}
with open(samples_informations_file, 'r') as f:
    for line in f:
        if line.startswith('SRR'):
            line = map(str.strip, line.split(','))
            sample_informations.setdefault(line[0], line[1])


cwd = filter(os.path.isdir, os.listdir(os.getcwd()))
all_available_sites = []
sample_edited_sites = {}
for directory in cwd:
    if directory.startswith('SRR'):
        path = list(os.walk(directory + '/editing/'))
        table = path[1][0] + '/' + path[1][-1][-1] 
        with open(table,'r') as a:
            for line in a:
                if line.startswith('chr'):
                    s = map(str.strip, line.split("\t"))
		    if s[7] == 'AG':
	                site, freq, coverage = s[0] + "_" + s[1], s[8], s[4]
			freq_gnum_cov = '%s^%s^%s' %(s[8],eval(s[6])[2],s[4]) 
			if site not in all_available_sites: all_available_sites.append(site)
			if (int(coverage) >= min_coverage) and (float(freq) >= min_edit_frequency):
                		sample_edited_sites.setdefault((directory, site), []).append((freq, freq_gnum_cov))

table_columns = map(lambda x: x + '_' + sample_informations[x], sorted(sample_informations.keys()))

disease = [i for i in table_columns if i.upper().find('DIS') != -1]
controls = [i for i in table_columns if i.upper().find('CTRL') != -1]

if enable_linear_model:
	outtable=''
        header = ['chromosome', 'position', 'type_editing'] + map(remove_underscore, controls) + map(remove_underscore, disease)
        outtable += '\t'.join(header)
	outtable += '\n'
        #print '\t'.join(header)
        for chrom in sorted(all_available_sites, key = lambda x: Set_Chr_Nr(x)):
                row = [chrom]
                for col in header[2:]:#header.index('[num_controls/num_disease]')]:
                        row.append(sample_edited_sites.get((col.split('_')[0],chrom), ['-'])[0])
                ctrls = zip(*(zip(controls,row[1:])))[1]
                dss = zip(*(zip(disease,row[len(ctrls)+1:])))[1]
                ctrls_freq = map(tuple_replace, ctrls)
                dss_freq = map(tuple_replace, dss)
                row.append(str([Sample_count(ctrls), Sample_count(dss)]))

                row_b = map(tuple_replace_bis, row)
                row_b = row_b[0].split('_') + row_b[2:]
                row_b.insert(2, 'A.to.G')
		final_list = row_b[:-1]
                #print '\t'.join(map(str,final_list))
		outtable += '\t'.join(map(str,final_list)).replace('-','NA')
		outtable += '\n'

	with open('temp.csv','w') as t:
		t.write(outtable)
		t.close()

	# call linear model script
	call_differential_editing_sites(samples_informations_file) 
	

else:
	header = ['chromosome', 'position', 'type_editing'] + controls + disease + ['[num_controls/num_disease]'] + ['delta_diff'] + ['pvalue (Mannwhitney)']

	if pvalue_correction == 1:
		header += ['pvalue Bonferroni corrected']
	if pvalue_correction == 2:
		header += ['pvalue BH corrected']
		
	print '\t'.join(header)
	
	for chrom in sorted(all_available_sites, key = lambda x: Set_Chr_Nr(x)):
		row = [chrom]
		for col in header[3:header.index('[num_controls/num_disease]')]:
			row.append(sample_edited_sites.get((col.split('_')[0],chrom), ['-'])[0])
		ctrls = zip(*(zip(controls,row[1:])))[1]
		dss = zip(*(zip(disease,row[len(ctrls)+1:])))[1] 
		ctrls_freq = map(tuple_replace, ctrls)
		dss_freq = map(tuple_replace, dss)
		row.append(str([Sample_count(ctrls), Sample_count(dss)]))
		if (Sample_percentage(ctrls) >= min_sample_testing) and (Sample_percentage(dss) >= min_sample_testing):
			ctrls_mean = sum(map(float, filter(lambda x: x!= '-', ctrls_freq)))/len(filter(lambda x: x!= '-', ctrls_freq))
	                dss_mean = sum(map(float, filter(lambda x: x!= '-', dss_freq)))/len(filter(lambda x : x!= '-', dss_freq))
			delta_diff =  abs(ctrls_mean - dss_mean)
			pvalue=stats.mannwhitneyu(ctrls_freq, dss_freq, alternative='two-sided')
			row.append(round(delta_diff, 3))
			row.append(str(round(pvalue[1], 3)))
			correction_argmnt = [(pvalue[1], ctrls_freq+dss_freq)]
		
			if pvalue_correction == 1:
				row.append(round(get_b(correction_argmnt, 0.05)[-1], 6))
			elif pvalue_correction == 2:
				row.append(round(get_bh(correction_argmnt, 0.05)[-1], 6))
		else:
			if pvalue_correction == 0:
				row += ['-', '-']
			else:
				row += ['-', '-', '-']
		row_a = map(tuple_replace, row)
		row_b = map(tuple_replace_bis, row)
		if pvalue_correction != 0 and only_significants == 'yes':
			only_sig(row_a,row_b)
		else:
			row_b =  row_b[0].split('_') + row_b[2:]
	                row_b.insert(2, 'A.to.G')
			print '\t'.join(map(str,row_b))

