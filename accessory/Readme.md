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
For each editing candidate, the script applies the Mann–Whitne<b>y test to check the significance between the two conditions, 
control and HD. By default the test is carried out only if the number of editing events per position is equal to 50% of the samples per group. 
Optionally, p-values can be corrected using Benjamini–Hochberg or Bonferroni tests. 
<ul>
<p>Usage:</p>
<pre><pre>
</body>
</html>
