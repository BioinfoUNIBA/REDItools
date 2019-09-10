#!/usr/bin/env python
"""
"""

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu, fisher_exact
import math

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-input_file", action = 'store', dest = 'config_file',
                type = str, default= 'empty', help = 'Comma separated file  e.g: SRR3306830,Control \
                SRR3306829,Healthy...etc')
args = parser.parse_args()
config_file = args.config_file


stability_value = 0.03 #value below which you may use a lower coverage for adding more samples to increase power
min_disease_people = 5 #min number people supporting higher coverage for whch you may base stability off measurements off of
min_control_people = 5  #min number control poeple supporting higher coverage for which you may base stability off of
min_disease_people_5_cov = 10 #min disease number of people of 5 coverage you must have if needing to use unstable 5x coverage
min_control_people_5_cov = 10 #min control number of people of 5 coverage you must have if needing to use unstable 5x coverage
editing_file= './temp.csv'
output_file = './editing_sites.with_stats_converted_disease.csv'


#######DONE CANGING PARAMETERS###################

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
full_table.to_csv(output_file, sep='\t')

print "job completed\n"

