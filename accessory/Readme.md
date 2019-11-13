<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />  
  </head>
  <body>
<h1>get_DE_events.py</h1>
<h5>This scripts and its related files are part of the supplemental material for the paper<br>
  "Investigating RNA editing in deep transcriptome datasets with REDItools and REDIportal"</h5>
<p class-text="justify">
For control case studies by launching the get_DE_events.py script the user can filter REDItoolDnaRna.py outputs according to the following criteria:
<ul>
<li>RNAseq coverage per position (default <b>10 reads</b>)</li>
<li>Minimum editing frequency per position (default <b>10%</b>)</li>
For each editing candidate, the script applies the Mann–Whitney test to check the significance between the two conditions, 
control and HD. By default the test is carried out only if the number of editing events per position is equal to 50% of the samples per group. 
Optionally, p-values can be corrected using Benjamini–Hochberg or Bonferroni tests. 
</ul>
<p>Usage:</p> 
<pre>
get_DE_events.py [-h] [-c MIN_COVERAGE] [-cpval PVALUE_CORRECTION]
                        [-input_file SAMPLES_INFORMATIONS_FILE]
                        [-f MIN_EDIT_FREQUENCY] [-mts MIN_SAMPLE_TESTING]
                        [-sig ONLY_SIGNIFICANT] [-linear]
  
optional arguments:
  -h, --help                             show this help message and exit
  -c MIN_COVERAGE                        Coverage-q30
  -cpval PVALUE_CORRECTION 1 -->         Bonferroni correction / 2 --> Benjamini hochberg
  -input_file SAMPLES_INFORMATIONS_FILE  Comma separated file e.g: <b>Sample,Status</b>
  -f MIN_EDIT_FREQUENCY                  Editing Frequency
  -mts MIN_SAMPLE_TESTING                min percentage of each sample category
  -sig ONLY_SIGNIFICANT                  Return only significant editing events 
                                         (if -cpval flag is activated)
  -linear                                Calculate differential RNA editing using 
                                                                                       the methodology by Tran et al. (2019)
 
<b>e.g.</b> python ../REDItools/accessory/get_DE_events.py -cpval 2 -input_file  sample_information.csv -sig yes
<p class-text="justify">The script will filter REDItoolDnaRna.py outputs for each sample contained in the 
SAMPLES_INFORMATIONS_FILE returning only significant editing events (pval <= 0.05)
in accordance with Benjamini hochberg correction.</p>

</pre>
</body>
</html> 
