<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />  
  </head>
  <body>
  <h3><a>Table Of Contents</a></h3>
  <ul>
<li><a class="reference internal" href="#">REDItools: python scripts for RNA editing detection by RNA-Seq data</a></li>
<li><a class="reference internal" href="#introduction">Introduction</a></li>
<li><a class="reference internal" href="#publications">Publications</a></li>
<li><a class="reference internal" href="#download">Download</a><ul>
<li><a class="reference internal" href="#download-reditools">Download REDItools</a></li>
<li><a class="reference internal" href="#download-testing-dataset">Download testing dataset</a></li>
</ul>
</li>
<li><a class="reference internal" href="#prerequisites">Prerequisites</a></li>
<li><a class="reference internal" href="#installation">Installation</a></li>
<li><a class="reference internal" href="#input-files">Input files</a></li>
<li><a class="reference internal" href="#preparing-input-files">Preparing input files</a><ul>
<li><a class="reference internal" href="#fasta">FASTA</a></li>
<li><a class="reference internal" href="#bam">BAM</a></li>
<li><a class="reference internal" href="#gtf">GTF</a></li>
<li><a class="reference internal" href="#tab">TAB</a></li>
<li><a class="reference internal" href="#splicesites">SpliceSites</a></li>
</ul>
</li>
<li><a class="reference internal" href="#output-files">Output files</a></li>
<li><a class="reference internal" href="#usage-information">Usage Information</a><ul>
<li><a class="reference internal" href="#reditooldnarna-py">REDItoolDnaRna.py</a></li>
<li><a class="reference internal" href="#reditoolknown-py">REDItoolKnown.py</a></li>
<li><a class="reference internal" href="#reditooldenovo-py">REDItoolDenovo.py</a></li>
</ul>
</li>
<li><a class="reference internal" href="#additional-scripts">Additional scripts:</a><ul>
<li><a class="reference internal" href="#reditoolblatcorrection-py">REDItoolBlatCorrection.py</a></li>
<li><a class="reference internal" href="#filtertable-py">FilterTable.py</a></li>
<li><a class="reference internal" href="#annotatetable-py">AnnotateTable.py</a></li>
<li><a class="reference internal" href="#searchintable-py">SearchInTable.py</a></li>
<li><a class="reference internal" href="#selectpositions-py">selectPositions.py</a></li>
<li><a class="reference internal" href="#sorttable-py-new-in-version-1-0-3">SortTable.py (new in version 1.0.3)</a></li>
<li><a class="reference internal" href="#tabletotabix-py-new-in-version-1-0-3">tableToTabix.py (new in version 1.0.3)</a></li>
<li><a class="reference internal" href="#sortgff-py-new-in-version-1-0-3">SortGFF.py (new in version 1.0.3)</a></li>
<li><a class="reference internal" href="#gfftotabix-py-new-in-version-1-0-3">GFFtoTabix.py (new in version 1.0.3)</a></li>
<li><a class="reference internal" href="#tabletogff-py-new-in-version-1-0-3">TableToGFF.py (new in version 1.0.3)</a></li>
</ul>
</li>
<li><a class="reference internal" href="#usage-example">Usage Example</a><ul>
<li><a class="reference internal" href="#test-dataset">Test dataset</a></li>
</ul>
</li>
<li><a class="reference internal" href="#contact">Contact</a></li>
</ul>


            
  <div class="toctree-wrapper compound">
<ul class="simple">
</ul>
</div>
<div class="section" id="reditools-python-scripts-for-rna-editing-detection-by-rna-seq-data">
<h1>REDItools: python scripts for RNA editing detection by RNA-Seq data<a class="headerlink" href="#reditools-python-scripts-for-rna-editing-detection-by-rna-seq-data" title="Permalink to this headline"></a></h1>
</div>
<div class="section" id="introduction">
<h1>Introduction<a class="headerlink" href="#introduction" title="Permalink to this headline"></a></h1>
<p>REDItools are python scripts developed with the aim to study RNA editing at genomic scale
by next generation sequencing data. RNA editing is a post-transcriptional phenomenon
involving the insertion/deletion or substitution of specific bases in precise RNA localizations.
In human, RNA editing occurs by deamination of cytosine to uridine (C-to-U) or mostly by the
adenosine to inosine (A-to-I) conversion through ADAR enzymes. A-to-I substitutions may have
profound functional consequences and have been linked to a variety of human diseases including
neurological and neurodegenerative disorders or cancer. Next generation sequencing technologies
offer the unique opportunity to investigate in depth RNA editing even though no dedicated
software has been released up to now.</p>
<p>REDItools are simple python scripts conceived to facilitate the investigation of RNA editing
at large-scale and devoted to research groups that would to explore such phenomenon in own
data but don’t have sufficient bioinformatics skills.
They work on main operating systems (although unix/linux-based OS are preferred), can handle reads from whatever
platform in the standard BAM format and implement a variety of filters.</p>
</div>
<div class="section" id="publications">
<h1>Publications<a class="headerlink" href="#publications" title="Permalink to this headline"></a></h1>
<ul>
<li><p class="first"><em>Picardi E, D&#8217;Erchia AM, Gallo A and Pesole G</em></p>
<p><strong>Detection of post-transcriptional RNA editing events.</strong></p>
<p><em>Methods Mol Biol. 2015;1269:189-205</em></p>
<p>PubMed PMID: <a class="reference external" href="http://www.ncbi.nlm.nih.gov/pubmed/25577380">25577380</a></p>
</li>
<li><p class="first"><em>Picardi E, D&#8217;Erchia AM, Gallo A, Montalvo A and Pesole G.</em></p>
<p><strong>Uncovering RNA Editing Sites in Long Non-Coding RNAs.</strong></p>
<p><em>Front Bioeng Biotechnol. 2014 Dec 5;2:64</em></p>
<p>PubMed PMID: <a class="reference external" href="http://www.ncbi.nlm.nih.gov/pubmed/25538940">25538940</a></p>
</li>
<li><p class="first"><em>Picardi E and Pesole G.</em></p>
<p><strong>REDItools: high-throughput RNA editing detection made easy.</strong></p>
<p><em>Bioinformatics. 2013 Jul 15;29(14):1813-4</em></p>
<p>PubMed PMID: <a class="reference external" href="http://www.ncbi.nlm.nih.gov/pubmed/23742983">23742983</a></p>
</li>
</ul>
</div>
<div class="section" id="download">
<h1>Download<a class="headerlink" href="#download" title="Permalink to this headline"></a></h1>
<div class="section" id="download-reditools">
<h2>Download old REDItools packages<a class="headerlink" href="#download-reditools" title="Permalink to this headline"></a></h2>
<ul class="simple">
<li><a class="reference external" href="https://sourceforge.net/projects/reditools/files/REDItools-1.2.1.zip/download">REDItools v1.2.1</a></li>
<li><a class="reference external" href="https://sourceforge.net/projects/reditools/files/REDItools-1.2.tar.gz/download">REDItools v1.2.0</a></li>
<li><a class="reference external" href="http://sourceforge.net/projects/reditools/files/REDItools-1.0.4.tar.gz/download">REDItools v1.0.4</a></li>
<li><a class="reference external" href="http://sourceforge.net/projects/reditools/files/REDItools-1.0.3.tar.gz/download">REDItools v1.0.3</a></li>
<li><a class="reference external" href="http://sourceforge.net/projects/reditools/files/REDItools-1.0.2.tar.gz/download">REDItools v1.0.2</a></li>
</ul>
</div>
<div class="section" id="download-testing-dataset">
<h2>Download testing dataset<a class="headerlink" href="#download-testing-dataset" title="Permalink to this headline"></a></h2>
<ul class="simple">
<li><a class="reference external" href="http://srv00.recas.ba.infn.it/reditools/data/testREDItools.tar.gz">testREDItools.tar.gz</a></li>
</ul>
</div>
</div>
<div class="section" id="prerequisites">
<h1>Prerequisites<a class="headerlink" href="#prerequisites" title="Permalink to this headline"></a></h1>
<p>REDItools require python 2.7 (available at the official python <a class="reference external" href="http://www.python.org">web-site</a>) while python 3 is not yet supported.
In addition, REDItools need two external modules:</p>
<ul class="simple">
<li>pysam (mandatory) version >= 0.15 available <a class="reference external" href="https://pysam.readthedocs.io/en/latest/installation.html">here</a></li>
<li>fisher (optional) version 0.1.4 (optional) available at python <a class="reference external" href="http://pypi.python.org/pypi/fisher/">web site</a></li>
</ul>
<p>To perform Blat correction and format alignment exchanges (SAM to BAM and vice versa) the
following packages should be installed or already present in your path:</p>
<ul class="simple">
<li><a class="reference external" href="http://hgdownload.cse.ucsc.edu/admin/exe/">Blat</a> package including gfServer and gfClient executables or <a class="reference external" href="http://icebert.github.io/pblat/">pblat</a></li>  
<li><a class="reference external" href="http://www.htslib.org/">Samtools</a> and <a class="reference external" href="http://www.htslib.org/">tabix</a></li>
</ul>
</div>
<div class="section" id="installation">
<h1>Installation<a class="headerlink" href="#installation" title="Permalink to this headline"></a></h1>
<p>REDItools require pysam module. To check its presence in your system try the
command:</p>
<div class="highlight-python"><pre>python -c 'import pysam'</pre>
</div>
<p>If no errors are raised, the module is already installed. Alternatively, download a copy
from the above link and follow instructions inside.</p>
<p>REDItools are stand-alone scripts and each one can be simply launched by:</p>
<div class="highlight-python"><pre>python REDItoolscript.py</pre>
</div>
<p>Alternatively you can easily install REDItools by the following commands:</p>
<div class="highlight-python">
<pre>
git clone https://github.com/BioinfoUNIBA/REDItools
cd REDItools
python setup.py install
</pre>
<p>For old packages use the following lines:</p>
<pre>gunzip reditools-1.0.2.tar.gz
tar –xvf reditools-1.0.2.tar
cd reditools-1.0.2
python setup.py install
</pre>
</div>
<p>To install scripts in a specific location the &#8211;prefix option can be used:</p>
<div class="highlight-python"><pre>python setup.py install --prefix=/my/path</pre>
</div>
<p>All scripts will be located in /my/path/bin.</p>
</div>
<div class="section" id="input-files">
<h1>Input files<a class="headerlink" href="#input-files" title="Permalink to this headline"></a></h1>
<p>REDItools have been designed to handle mapped reads from next generation sequencing
technologies in the standard BAM format.
The following input files are accepted:</p>
<ul class="simple">
<li>BAM files containing RNA-Seq or DNA-Seq reads aligned onto a reference genome (mandatory);</li>
<li>FASTA file of reference genome (mandatory);</li>
<li>GTF file containing specific genomic positions (optional);</li>
<li>TAB file containing individual genomic positions (optional);</li>
<li>SpliceSites file containing splice sites information (optional).</li>
</ul>
</div>
<div class="section" id="preparing-input-files">
<h1>Preparing input files<a class="headerlink" href="#preparing-input-files" title="Permalink to this headline"></a></h1>
<div class="section" id="fasta">
<h2>FASTA<a class="headerlink" href="#fasta" title="Permalink to this headline"></a></h2>
<p>Reference genome in fasta format needs to be indexed using samtools:</p>
<div class="highlight-python"><pre>samtools faidx reference.fa</pre>
</div>
</div>
<div class="section" id="bam">
<h2>BAM<a class="headerlink" href="#bam" title="Permalink to this headline"></a></h2>
<p>Starting from alignments in SAM format, you first need to convert SAM to BAM:</p>
<div class="highlight-python"><pre>samtools view –bS –T reference.fa myfile.sam &gt; myfile.bam</pre>
</div>
<p>then sort your BAM:</p>
<div class="highlight-python"><pre>samtools sort myfile.bam myfile.sorted</pre>
</div>
<p>and finally index your sorted BAM:</p>
<div class="highlight-python"><pre>samtools index myfile.sorted.bam</pre>
</div>
</div>
<div class="section" id="gtf">
<h2>GTF<a class="headerlink" href="#gtf" title="Permalink to this headline"></a></h2>
<p>GTF files for specific annotations and genomic regions can be downloaded from <a class="reference external" href="http://genome.ucsc.edu/">UCSC</a>.
Alternatively, UCSC tables in psl format may be used and converted in GTF by the following awk/gawk commands:</p>
<ul>
<li><p class="first">For RepeakMask table:</p>
<div class="highlight-python"><pre>$ gunzip rmsk.txt.gz
$ gawk 'OFS="\t"{print $6,"rmsk_hg19",$12,$7+1,$8,".",$10,".","gene_id \""$11"\"; transcript_id \""$13"\";"}' rmsk.txt &gt; rmsk.gtf</pre>
</div>
</li>
<li><p class="first">For simpleRepeat table:</p>
<div class="highlight-python"><pre>$ gunzip simpleRepeat.txt.gz
$ gawk 'OFS="\t"{print $2,"trf_hg19","trf",$3+1,$4,".","+",".","gene_id \""$5"\"; transcript_id \""$17"\";"}' simpleRepeat.txt &gt; simpleRepeat.gtf</pre>
</div>
</li>
<li><p class="first">For RefSeq annotations you can download from UCSC and use the <a class="reference external" href="http://hgdownload.soe.ucsc.edu/admin/exe/">genePredToGtf</a> utility:</p>
<div class="highlight-python"><pre>$ gunzip refGene.txt.gz
$ cut -f 2- refGene.txt | genePredToGtf -utr -source=hg19_refseq file stdin refGene.gtf</pre>
</div>
</li>
</ul>
<p>Other annotations can be generated accordingly. Please, consider that the GTF format includes 9 Tab-delimited fields:</p>
<ul class="simple">
<li>chromosome/region name</li>
<li>source (for example hg19_refseq)</li>
<li>feature (for example CDS, exon, snp ...)</li>
<li>start coordinate (1-based)</li>
<li>end coordinate (1-based)</li>
<li>score (use a dot if unknown)</li>
<li>strand (+ or - use + if unknown)</li>
<li>frame for CDS only or a dot (a dot also if unknown)</li>
<li>attributes. For REDItools there are two mandatory attributes: gene_id and transcript_id</li>
</ul>
<p>Following lines show a GTF example for RefSeq annotations:</p>
<div class="highlight-python"><pre>$ head -n6 refGene.gtf
chr7    hg19_refseq     5UTR    139025878       139026130       .       +       .       gene_id "C7orf55"; transcript_id "NM_197964";
chr7    hg19_refseq     CDS     139026131       139026268       .       +       0       gene_id "C7orf55"; transcript_id "NM_197964";
chr7    hg19_refseq     CDS     139030247       139030447       .       +       0       gene_id "C7orf55"; transcript_id "NM_197964";
chr7    hg19_refseq     3UTR    139030451       139031065       .       +       .       gene_id "C7orf55"; transcript_id "NM_197964";
chr7    hg19_refseq     start_codon     139026131       139026133       .       +       0       gene_id "C7orf55"; transcript_id "NM_197964";
chr7    hg19_refseq     stop_codon      139030448       139030450       .       +       0       gene_id "C7orf55"; transcript_id "NM_197964";</pre>
</div>
<p>GTF files can be sorted by the command:</p>
<div class="highlight-python"><pre>$ sort -k1,1 -k4,4n mygft &gt; mygtf_sorted</pre>
</div>
</div>
<div class="section" id="tab">
<h2>TAB<a class="headerlink" href="#tab" title="Permalink to this headline"></a></h2>
<p>TAB files are simple textual files with at least three tabulated columns including:</p>
<ul class="simple">
<li>genomic region (generally the chromosome name according to the reference genome)</li>
<li>coordinate of the position (1-based)</li>
<li>strand (+ or -). You can also indicate strand by 0 (strand -), 1 (strand +) or 2 (+ and - or unknown)</li>
</ul>
<table border="0" class="docutils">
<thead valign="bottom">
<tr class="row-odd"><th class="head">genomic region</th>
<th class="head">coordinate</th>
<th class="head">strand</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>chr21</td>
<td>10205589</td>
<td>-</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>10205629</td>
<td>-</td>
</tr>
<tr class="row-even"><td>chr21</td>
<td>15411496</td>
<td>+</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15412990</td>
<td>+</td>
</tr>
<tr class="row-even"><td>chr21</td>
<td>15414553</td>
<td>+</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15415901</td>
<td>+</td>
</tr>
<tr class="row-even"><td>chr21</td>
<td>15417667</td>
<td>+</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15423330</td>
<td>+</td>
</tr>
</tbody>
</table>
<p>TAB files must be coordinate sorted. In unix/linux environment they can be sorted by the sort command:</p>
<div class="highlight-python"><pre>sort -k1,1 -k2,2n mytable.txt &gt; mytable.sorted.txt</pre>
</div>
</div>
<div class="section" id="splicesites">
<h2>SpliceSites<a class="headerlink" href="#splicesites" title="Permalink to this headline"></a></h2>
<p>SpliceSites files are simple textual file including 5 columns separated by space:</p>
<ul class="simple">
<li>genomic region (chromosome/region name according to reference and input BAMs)</li>
<li>start coordinate of splice site (1-based)</li>
<li>end coordinate of splice site (1-based)</li>
<li>splice site type, A for acceptor and D for donor</li>
<li>strand (+ or -)</li>
</ul>
<table border="0" class="docutils">
<thead valign="bottom">
<tr class="row-odd"><th class="head">genomic region</th>
<th class="head">start</th>
<th class="head">end</th>
<th class="head">type</th>
<th class="head">strand</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>chr1</td>
<td>13224</td>
<td>13225</td>
<td>A</td>
<td>+</td>
</tr>
<tr class="row-odd"><td>chr1</td>
<td>12227</td>
<td>12228</td>
<td>D</td>
<td>+</td>
</tr>
<tr class="row-even"><td>chr1</td>
<td>12594</td>
<td>12595</td>
<td>A</td>
<td>+</td>
</tr>
<tr class="row-odd"><td>chr1</td>
<td>12721</td>
<td>12722</td>
<td>D</td>
<td>+</td>
</tr>
<tr class="row-even"><td>chr1</td>
<td>13402</td>
<td>13403</td>
<td>A</td>
<td>+</td>
</tr>
<tr class="row-odd"><td>chr1</td>
<td>13655</td>
<td>13656</td>
<td>D</td>
<td>+</td>
</tr>
<tr class="row-even"><td>chr1</td>
<td>29320</td>
<td>29321</td>
<td>D</td>
<td>-</td>
</tr>
<tr class="row-odd"><td>chr1</td>
<td>24891</td>
<td>24892</td>
<td>A</td>
<td>-</td>
</tr>
<tr class="row-even"><td>chr1</td>
<td>24737</td>
<td>24738</td>
<td>D</td>
<td>-</td>
</tr>
<tr class="row-odd"><td>chr1</td>
<td>18379</td>
<td>18380</td>
<td>A</td>
<td>-</td>
</tr>
</tbody>
</table>
<p>SpliceSites files can be obtained starting from UCSC annotations in psl format and using
the perl script psl_splicesites included in the <a class="reference external" href="http://research-pub.gene.com/gmap/">GMAP/GSNAP</a> package:</p>
<div class="highlight-python"><pre>gunzip -c refGene.txt.gz | psl_splicesites -s 1 &gt; mysplicesites</pre>
</div>
<p>The format of mysplicesites is:</p>
<div class="highlight-python"><pre>&gt;NM_005236.exon11/11 chr16:14041470..14041471 acceptor 2778
&gt;NM_005235.exon1/28 chr2:213403173..213403172 donor 413544
&gt;NM_005235.exon2/28 chr2:212989629..212989628 acceptor 413544
&gt;NM_005235.exon2/28 chr2:212989477..212989476 donor 177135
&gt;NM_005235.exon3/28 chr2:212812342..212812341 acceptor 177135
&gt;NM_005235.exon3/28 chr2:212812155..212812154 donor 159270
&gt;NM_005235.exon4/28 chr2:212652885..212652884 acceptor 159270
&gt;NM_005235.exon4/28 chr2:212652750..212652749 donor 37320
&gt;NM_005235.exon5/28 chr2:212615430..212615429 acceptor 37320</pre>
</div>
<p>Then, the following awk/gawk command line can be used to get the final SpliceSite file:</p>
<div class="highlight-python"><pre>$ gawk -F" " '{split($2,a,":"); split(a[2],b,"."); if (b[1]&gt;b[3]) print a[1],b[3],b[1],toupper(substr($3,1,1)),"-"; else print a[1],b[1],b[3],toupper(substr($3,1,1)),"+"}' mysplicesites &gt; mysplicesites.ss
$ more mysplicesites.ss
chr16 14041470 14041471 A +
chr2 213403172 213403173 D -
chr2 212989628 212989629 A -
chr2 212989476 212989477 D -
chr2 212812341 212812342 A -
chr2 212812154 212812155 D -
chr2 212652884 212652885 A -
chr2 212652749 212652750 D -
chr2 212615429 212615430 A -</pre>
</div>
</div>
</div>
<div class="section" id="output-files">
<h1>Output files<a class="headerlink" href="#output-files" title="Permalink to this headline"></a></h1>
<p>REDItools print out results in simple textual tables:</p>
<table border="0" class="docutils">
<thead valign="bottom">
<tr class="row-odd"><th class="head">Region</th>
<th class="head">Position</th>
<th class="head">Reference</th>
<th class="head">Strand</th>
<th class="head">Coverage-q25</th>
<th class="head">MeanQ</th>
<th class="head">BaseCount[A,C,G,T]</th>
<th class="head">AllSubs</th>
<th class="head">Frequency</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>chr21</td>
<td>15412990</td>
<td>A</td>
<td>1</td>
<td>18</td>
<td>37.22</td>
<td>[3, 0, 15, 0]</td>
<td>AG</td>
<td>0.83</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15415901</td>
<td>A</td>
<td>1</td>
<td>13</td>
<td>37.15</td>
<td>[2, 0, 11, 0]</td>
<td>AG</td>
<td>0.85</td>
</tr>
<tr class="row-even"><td>chr21</td>
<td>15423330</td>
<td>A</td>
<td>1</td>
<td>11</td>
<td>38.27</td>
<td>[4, 0, 7, 0]</td>
<td>AG</td>
<td>0.64</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15425640</td>
<td>A</td>
<td>1</td>
<td>8</td>
<td>36.12</td>
<td>[0, 0, 8, 0]</td>
<td>AG</td>
<td>1.00</td>
</tr>
<tr class="row-even"><td>chr21</td>
<td>15456434</td>
<td>T</td>
<td>1</td>
<td>90</td>
<td>34.96</td>
<td>[0, 6, 1, 83]</td>
<td>TC TG</td>
<td>0.07</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15461406</td>
<td>A</td>
<td>1</td>
<td>83</td>
<td>37.27</td>
<td>[73, 0, 10, 0]</td>
<td>AG</td>
<td>0.12</td>
</tr>
<tr class="row-even"><td>chr21</td>
<td>15461417</td>
<td>A</td>
<td>1</td>
<td>90</td>
<td>36.26</td>
<td>[72, 0, 18, 0]</td>
<td>AG</td>
<td>0.20</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15461444</td>
<td>A</td>
<td>1</td>
<td>64</td>
<td>37.22</td>
<td>[26, 0, 38, 0]</td>
<td>AG</td>
<td>0.59</td>
</tr>
<tr class="row-even"><td>chr21</td>
<td>15461479</td>
<td>A</td>
<td>1</td>
<td>70</td>
<td>36.96</td>
<td>[66, 0, 4, 0]</td>
<td>AG</td>
<td>0.06</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15461486</td>
<td>A</td>
<td>1</td>
<td>68</td>
<td>37.06</td>
<td>[61, 0, 7, 0]</td>
<td>AG</td>
<td>0.10</td>
</tr>
<tr class="row-even"><td>chr21</td>
<td>15461503</td>
<td>A</td>
<td>1</td>
<td>76</td>
<td>37.26</td>
<td>[69, 0, 7, 0]</td>
<td>AG</td>
<td>0.09</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15461511</td>
<td>A</td>
<td>1</td>
<td>81</td>
<td>37.68</td>
<td>[55, 0, 26, 0]</td>
<td>AG</td>
<td>0.32</td>
</tr>
</tbody>
</table>
<dl class="docutils">
<dt>where:</dt>
<dd><ul class="first last simple">
<li>Region: is the genomic region according to reference</li>
<li>Position: is the exact genomic coordinate (1-based)</li>
<li>Reference: is the nucleotide base in reference genome</li>
<li>Strand: is strand information with notation 1 for + strand, 0 for - strand and 2 for unknown or not defined strand</li>
<li>Coverage-qxx: is the depth per site at a given xx quality score (min. value)</li>
<li>MeanQ: is the mean quality score per site</li>
<li>BaseCount[A,C,G,T]: is the base distribution per site in the order A, C, G and T</li>
<li>AllSubs: is the list of observed substitutions at a given site, separated by a space. A character &#8220;-&#8221; is included in case of invariant sites.</li>
<li>Frequency: is the observed frequency of substitution. In case of multiple substitutions, it refers to the first in the AllSubs field.
In the table above, for example, at line 5 the frequency 0.07 is linked to TC substitution and not to TG that will be always lower than this value.</li>
</ul>
</dd>
<dt>REDItoolDnaRna.py includes five additional columns to take into account information from DNA-Seq reads. Such columns are indicated as:</dt>
<dd><ul class="first last simple">
<li>gCoverage-q25: is the depth per site at a given xx quality score (min. value) in DNA-Seq</li>
<li>gMeanQ: is the mean quality score per site in DNA-Seq</li>
<li>gBaseCount[A,C,G,T]: is the base distribution per site in the order A, C, G and T</li>
<li>gAllSubs: is the list of observed substitutions at a given site, separated by a space. A character &#8220;-&#8221; is included in case of invariant sites.</li>
<li>gFrequency: is the observed frequency of substitution. In case of multiple substitutions, it refers to the first in the gAllSubs field.</li>
</ul>
</dd>
</dl>
<p>In case of positions not supported by DNA-Seq reads, a character &#8220;-&#8221; will be added to each extra column.</p>
<dl class="docutils">
<dt>REDItoolDenovo.py includes one additional column concerning the reliability of editing prediction.</dt>
<dd><ul class="first last simple">
<li>Pvalue: is the pvalue per site calculated according to Fisher exact test. It indicates how much the observed base distribution for a change is different from the expected,
calculated by the empirical base substitution for the entire RNA-Seq experiment.</li>
</ul>
</dd>
</dl>
<p><strong>Note</strong> that BaseCount and Reference column are modified according to Strand column. In case of Strand value of 2, the genomic Reference is used.</p>
</div>
<div class="section" id="usage-information">
<h1>Usage Information<a class="headerlink" href="#usage-information" title="Permalink to this headline"></a></h1>
<div class="section" id="reditooldnarna-py">
<h2>REDItoolDnaRna.py<a class="headerlink" href="#reditooldnarna-py" title="Permalink to this headline"></a></h2>
<p>REDItoolDnaRna.py is the main script devoted to the identification of RNA editing events
taking into account the combined information from RNA-Seq and DNA-Seq data in BAM format.
To look at potential RNA editing candidates, RNA-Seq data alone can be used.</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>RNA-Seq BAM file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-j</span></kbd></td>
<td>DNA-Seq BAM files separated by comma or folder
containing BAM files. <strong>Note</strong> that each chromosome/region
must be present in a single BAM file only.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-I</span></kbd></td>
<td>Sort input RNA-Seq BAM file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-J</span></kbd></td>
<td>Sort input DNA-Seq BAM file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-f</span></kbd></td>
<td>Reference file in fasta format. <strong>Note</strong> that chromosome/region names
in the reference must match chromosome/region names in BAMs files.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-C</span></kbd></td>
<td>Base interval to explore [100000]. It indicates how many bases have to be loaded
during the run.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-k</span></kbd></td>
<td>List of chromosomes to skip separated by comma or file (each line must contain a chromosome/region name).</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-t</span></kbd></td>
<td>Number of threads [1]. It indicates how many processes should be launched. Each process will work on an individual chromosome/region.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Output folder [rediFolder_XXXX] in which all results will be stored. XXXX is a random number generated at each run.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-F</span></kbd></td>
<td>Internal folder name [null] is the main folder containing output tables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-M</span></kbd></td>
<td>Save a list of columns with quality scores. It produces at most two files in the pileup-like format.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-c</span></kbd></td>
<td>Minimum read coverage (dna,rna) [10,10]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-Q</span></kbd></td>
<td>Fastq offset value (dna,rna) [33,33]. For Illumina fastq 1.3+ 64 should be used.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-q</span></kbd></td>
<td>Minimum quality score (dna,rna) [25,25]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-m</span></kbd></td>
<td>Minimum mapping quality score (dna,rna) [25,25]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-O</span></kbd></td>
<td>Minimum homoplymeric length (dna,rna) [5,5]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Infer strand (for strand oriented reads) [1]. It indicates which read is in line with RNA. Available values are: 1:read1 as RNA,read2 not as RNA;
2:read1 not as RNA,read2 as RNA; 12:read1 as RNA,read2 as RNA; 0:read1 not as RNA,read2 not as RNA.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-g</span></kbd></td>
<td>Strand inference type 1:maxValue 2:useConfidence [1]; maxValue: the most prominent strand count will be used; useConfidence: strand is assigned if over a prefixed frequency confidence (-x option)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-x</span></kbd></td>
<td>Strand confidence [0.70]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-S</span></kbd></td>
<td>Strand correction. Once the strand has been inferred, only bases according to this strand will be selected.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-G</span></kbd></td>
<td>Infer strand by GFF annotation (must be GFF and sorted, otherwise use -X). Sorting requires grep and sort unix executables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-K</span></kbd></td>
<td>GFF File with positions to exclude (must be GFF and sorted, otherwise use -X). Sorting requires grep and sort unix executables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-T</span></kbd></td>
<td>Work only on given GFF positions (must be GFF and sorted, otherwise use -X). Sorting requires grep and sort unix executables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-X</span></kbd></td>
<td>Sort annotation files. It requires grep and sort unix executables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-e</span></kbd></td>
<td>Exclude multi hits in RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-E</span></kbd></td>
<td>Exclude multi hits in DNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-d</span></kbd></td>
<td>Exclude duplicates in RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-D</span></kbd></td>
<td>Exclude duplicates in DNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-p</span></kbd></td>
<td>Use paired concardant reads only in RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-P</span></kbd></td>
<td>Use paired concardant reads only in DNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-u</span></kbd></td>
<td>Consider mapping quality in RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-U</span></kbd></td>
<td>Consider mapping quality in DNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-a</span></kbd></td>
<td>Trim x bases up and y bases down per read [0-0] in RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-A</span></kbd></td>
<td>Trim x bases up and y bases down per read [0-0] in DNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-b</span></kbd></td>
<td>Blat folder for correction in RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-B</span></kbd></td>
<td>Blat folder for correction in DNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-l</span></kbd></td>
<td>Remove substitutions in homopolymeric regions in RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-L</span></kbd></td>
<td>Remove substitutions in homopolymeric regions in DNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-v</span></kbd></td>
<td>Minimum number of reads supporting the variation [3] for RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-n</span></kbd></td>
<td>Minimum editing frequency [0.1] for RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-N</span></kbd></td>
<td>Minimum variation frequency [0.1] for DNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-z</span></kbd></td>
<td>Exclude positions with multiple changes in RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-Z</span></kbd></td>
<td>Exclude positions with multiple changes in DNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-W</span></kbd></td>
<td>Select RNA-Seq positions with defined changes (separated by comma ex: AG,TC) [default all]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-R</span></kbd></td>
<td>Exclude invariant RNA-Seq positions</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-V</span></kbd></td>
<td>Exclude sites not supported by DNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-w</span></kbd></td>
<td>File containing splice sites annotations (SpliceSite file format see above for details)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-r</span></kbd></td>
<td>Num. of bases near splice sites to explore [4]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">--gzip</span></kbd></td>
<td>Gzip output files</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span>, <span class="option">--help</span></kbd></td>
<td>Print the help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example:</p>
<div class="highlight-python"><pre>REDItoolDnaRna.py -i rnaseq.bam -j dnaseq.bam -f myreference.fa -o myoutputfolder</pre>
</div>
</div>
<div class="section" id="reditoolknown-py">
<h2>REDItoolKnown.py<a class="headerlink" href="#reditoolknown-py" title="Permalink to this headline"></a></h2>
<p>REDItoolKnown.py has been developed to explore the RNA editing potential of RNA-Seq data
sets using known editing events. Such events can be downloaded from DARNED database or
generated from supplementary materials of a variety of publications. Known RNA editing events
have to be stored in TAB files (see above for details).</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>BAM file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-I</span></kbd></td>
<td>Sort input BAM file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-f</span></kbd></td>
<td>Reference in fasta file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-l</span></kbd></td>
<td>List of known RNA editing events</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-C</span></kbd></td>
<td>Base interval to explore [100000]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-k</span></kbd></td>
<td>List of chromosomes to skip separated by comma or file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-t</span></kbd></td>
<td>Number of threads [1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Output folder [rediFolder_XXXX] in which all results will be stored. XXXX is a random number generated at each run.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-F</span></kbd></td>
<td>Internal folder name [null] is the main folder containing output tables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-c</span></kbd></td>
<td>Min. read coverage [10]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-Q</span></kbd></td>
<td>Fastq offset value [33]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-q</span></kbd></td>
<td>Minimum quality score [25]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-m</span></kbd></td>
<td>Minimum mapping quality score [25]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-O</span></kbd></td>
<td>Minimum homoplymeric length [5]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Infer strand (for strand oriented reads) [1]. It indicates which read is in line with RNA. Available values are: 1:read1 as RNA,read2 not as RNA;
2:read1 not as RNA,read2 as RNA; 12:read1 as RNA,read2 as RNA; 0:read1 not as RNA,read2 not as RNA.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-g</span></kbd></td>
<td>Strand inference type 1:maxValue 2:useConfidence [1]; maxValue: the most prominent strand count will be used; useConfidence: strand is assigned if over a prefixed frequency confidence (-x option)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-x</span></kbd></td>
<td>Strand confidence [0.70]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-S</span></kbd></td>
<td>Strand correction. Once the strand has been inferred, only bases according to this strand will be selected.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-G</span></kbd></td>
<td>Infer strand by GFF annotation (must be sorted, otherwise use -X). Sorting requires grep and sort unix executables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-X</span></kbd></td>
<td>Sort annotation files. It requires grep and sort unix executables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-K</span></kbd></td>
<td>File with positions to exclude (chromosome_name coordinate)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-e</span></kbd></td>
<td>Exclude multi hits</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-d</span></kbd></td>
<td>Exclude duplicates</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-p</span></kbd></td>
<td>Use paired concardant reads only</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-u</span></kbd></td>
<td>Consider mapping quality</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-T</span></kbd></td>
<td>Trim x bases up and y bases down per read [0-0]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-B</span></kbd></td>
<td>Blat folder for correction</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-U</span></kbd></td>
<td>Remove substitutions in homopolymeric regions</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-v</span></kbd></td>
<td>Minimum number of reads supporting the variation [3]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-n</span></kbd></td>
<td>Minimum editing frequency [0.1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-E</span></kbd></td>
<td>Exclude positions with multiple changes</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-P</span></kbd></td>
<td>File containing splice sites annotations (SpliceSite file format see above for details)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-r</span></kbd></td>
<td>Num. of bases near splice sites to explore [4]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print the help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example:</p>
<div class="highlight-python"><pre>REDItoolKnown.py -i rnaseq.bam -f reference.fa -l knownEditingSites.tab</pre>
</div>
</div>
<div class="section" id="reditooldenovo-py">
<h2>REDItoolDenovo.py<a class="headerlink" href="#reditooldenovo-py" title="Permalink to this headline"></a></h2>
<p>REDItoolDenovo.py has been conceived to predict potential RNA editing events using RNA-Seq
data alone and without any a priori knowledge about genome information and biological properties
of RNA editing phenomenon.</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>BAM file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-I</span></kbd></td>
<td>Sort input BAM file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-f</span></kbd></td>
<td>Reference in fasta file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-k</span></kbd></td>
<td>List of chromosomes to skip separated by comma or file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-t</span></kbd></td>
<td>Number of threads [1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Output folder [rediFolder_XXXX] in which all results will be stored. XXXX is a random number generated at each run.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-F</span></kbd></td>
<td>Internal folder name [null] is the main folder containing output tables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-b</span></kbd></td>
<td>Use input distribution file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-a</span></kbd></td>
<td>Fisher Tail [l, r, t] [default l] [l=left_tail, r=right_tail, t=two_tail]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-c</span></kbd></td>
<td>Min. read coverage [10]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-Q</span></kbd></td>
<td>Fastq offset value [33]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-q</span></kbd></td>
<td>Min. quality score [25]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-m</span></kbd></td>
<td>Min. mapping quality score [25]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-O</span></kbd></td>
<td>Min. homoplymeric length [5]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Infer strand (for strand oriented reads) [1]. It indicates which read is in line with RNA. Available values are: 1:read1 as RNA,read2 not as RNA;
2:read1 not as RNA,read2 as RNA; 12:read1 as RNA,read2 as RNA; 0:read1 not as RNA,read2 not as RNA.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-g</span></kbd></td>
<td>Strand inference type 1:maxValue 2:useConfidence [1]; maxValue: the most prominent strand count will be used; useConfidence: strand is assigned if over a prefixed frequency confidence (-x option)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-x</span></kbd></td>
<td>Strand confidence [0.70]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-S</span></kbd></td>
<td>Strand correction. Once the strand has been inferred, only bases according to this strand will be selected.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-G</span></kbd></td>
<td>Infer strand by gff annotation (must be sorted, otherwise use -X)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-X</span></kbd></td>
<td>Sort annotation files</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-K</span></kbd></td>
<td>File with positions to exclude</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-e</span></kbd></td>
<td>Exclude multi hits</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-d</span></kbd></td>
<td>Exclude duplicates</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-l</span></kbd></td>
<td>Select significant sites</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-V</span></kbd></td>
<td>Significant value [0.05]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-w</span></kbd></td>
<td>Statistical test [BH, BO, NO] [default BH] [BH=Benjamini, BO=Bonferroni, NO=No Correction]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-U</span></kbd></td>
<td>Use specific substitutions separated by comma [example: AG,TC]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-p</span></kbd></td>
<td>Use paired concardant reads only</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-u</span></kbd></td>
<td>Consider mapping quality</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-T</span></kbd></td>
<td>Trim x bases up and y bases down per read [0-0]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-B</span></kbd></td>
<td>Blat folder for correction</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-W</span></kbd></td>
<td>Remove substitutions in homopolymeric regions</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-v</span></kbd></td>
<td>Minimum number of reads supporting the variation [3]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-n</span></kbd></td>
<td>Minimum editing frequency [0.1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-E</span></kbd></td>
<td>Exclude positions with multiple changes</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-P</span></kbd></td>
<td>File containing splice sites annotations (SpliceSite file format see above for details)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-r</span></kbd></td>
<td>Number of bases near splice sites to explore [4]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print the help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example:</p>
<div class="highlight-python"><pre>REDItoolDenovo.py -i rnaseq.bam -f reference.fa</pre>
</div>
</div>
</div>
<div class="section" id="additional-scripts">
<h1>Additional scripts:<a class="headerlink" href="#additional-scripts" title="Permalink to this headline"></a></h1>
<div class="section" id="reditoolblatcorrection-py">
<h2>REDItoolBlatCorrection.py<a class="headerlink" href="#reditoolblatcorrection-py" title="Permalink to this headline"></a></h2>
<p>REDItoolBlatCorrection.py requires gfServer and gfClient programs from <a class="reference external" href="http://hgdownload.cse.ucsc.edu/admin/exe/">Blat package</a> .
Reference fasta file can be converted in .2bit format using the utility <a class="reference external" href="http://hgdownload.cse.ucsc.edu/admin/exe/">faToTwoBit</a></p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>BAM file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-I</span></kbd></td>
<td>Sort input BAM file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-f</span></kbd></td>
<td>Genomic fasta file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-F</span></kbd></td>
<td>Genomic fasta file in 2bit format for gfServer</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-t</span></kbd></td>
<td>Num. of working threads [1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Output Folder [BlatCorrection_XXXX]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-k</span></kbd></td>
<td>List of chromosomes to skip separated by comma or file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-r</span></kbd></td>
<td>Regions in GFF in which Blat will be launched</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Sort GFF (if unsorted). It requires grep and sort unix executables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-q</span></kbd></td>
<td>Minimum quality score [25]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-Q</span></kbd></td>
<td>Fastq offset value [33]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-V</span></kbd></td>
<td>Verify if gfServer is alive</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-T</span></kbd></td>
<td>Stop gfServer at script end</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print the help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example:</p>
<div class="highlight-python"><pre>REDItoolBlatCorrection.py -i rnaseq.bam -f reference.fa -F reference.2bit -o BlatCorrection -V -T</pre>
</div>
<p>At the end of the run, in the BlatCorrection folder, there will be different .bad files, one per chromosome/region.
Each .bad file includes two space separated columns:</p>
<ul class="simple">
<li>read name</li>
<li>a number from 0 to 2 indicating the read mate; 1: if read is first in pair, 2: if read is second in pair, 0: if read is a single end</li>
</ul>
<p>And the content looks like:</p>
<div class="highlight-python"><pre>TUPAC_0006:2:41:10999:19783#0 1
TUPAC_0006:2:41:10999:19783#0 2
TUPAC_0006:3:54:8655:18923#0 2
PRESLEY_0005:1:73:13079:17565#0 1
TUPAC_0006:2:9:12695:14552#0 1
PRESLEY_0005:1:73:13079:17565#0 1
TUPAC_0006:2:9:12695:14552#0 1
PRESLEY_0005:1:12:1152:17918#0 2
SINATRA_0006:7:15:8730:10887#0 1
SINATRA_0006:7:67:12736:11713#0 1
SINATRA_0006:7:50:9125:3151#0 1
SINATRA_0006:8:49:4592:20505#0 1
SINATRA_0006:7:14:10225:3766#0 1
PRESLEY_0005:1:118:4480:20145#0 1
SINATRA_0006:7:6:19272:9901#0 1</pre>
</div>
</div>
<div class="section" id="filtertable-py">
<h2>FilterTable.py<a class="headerlink" href="#filtertable-py" title="Permalink to this headline"></a></h2>
<p>FilterTable.py filters positions of a input table according to specific annotations indexed by tabix tool.
Filtered out positions will be marked with &#8220;#&#8221; at the beginning of each line. To exclude such lines the
option -E should be used. Features are the same as indicated in the third field of GTF annotation file.</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>Table file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-f</span></kbd></td>
<td>Sorted file with positions to filter in (GTF sorted file as above)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Sorted file with positions to filter out (GTF sorted file as above)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-F</span></kbd></td>
<td>Features to filter in (separated by comma)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-S</span></kbd></td>
<td>Features to filter out (separated by comma)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-E</span></kbd></td>
<td>Exclude positions filtered out</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Save filtered lines to a file [stdout]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-p</span></kbd></td>
<td>Print simple statistics</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print the help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example:</p>
<div class="highlight-python"><pre>FilterTable.py -i mytable -s dbsnp137.gtf.gz -S snp -o mytable_filtered -E -p</pre>
</div>
</div>
<div class="section" id="annotatetable-py">
<h2>AnnotateTable.py<a class="headerlink" href="#annotatetable-py" title="Permalink to this headline"></a></h2>
<p>AnnotateTable.py annotates individual positions of a table file according to annotations indexed by tabix tool.
It adds from 1 to 3 additional columns according to -c option.</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-a</span></kbd></td>
<td>Sorted Annotation file in GTF format</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>Annotate a file of positions [column1=region, column2=coordinate (1 based)]
or a single position [region:coordinate (1 based)]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Strand column in annotation file [4]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-u</span></kbd></td>
<td>Not use table strand info (fix it to 2)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-c</span></kbd></td>
<td>Add columns separated by comma (feature:1, gene_id:2, transcript_id:3) [1,2]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-n</span></kbd></td>
<td>Column name [Col]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-S</span></kbd></td>
<td>Correct strand by annotation</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-C</span></kbd></td>
<td>Columns with base distribution [7,12] (in combination with -S)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Save lines to a file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print the help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example:</p>
<div class="highlight-python"><pre>AnnotateTable.py -i mytable -a rmsk.gtf.gz -u -c1,2,3 -n rmsk -o mytable_rmsk</pre>
</div>
</div>
<div class="section" id="searchintable-py">
<h2>SearchInTable.py<a class="headerlink" href="#searchintable-py" title="Permalink to this headline"></a></h2>
<p>SearchInTable.py looks for individual positions in a list or table of positions.</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>Sorted table file (first col=reference; second col=coordinate 1 based)
or tabix indexed table (ending with .gz) (TAB format is required)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-j</span></kbd></td>
<td>Query (file or single positions as chr21:123456)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-C</span></kbd></td>
<td>Sequence name column [1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-S</span></kbd></td>
<td>Start column [2]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-E</span></kbd></td>
<td>End column; can be identical to &#8216;-S&#8217; [2]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-p</span></kbd></td>
<td>Print position header (like a fasta header &gt;chr21:123456)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-n</span></kbd></td>
<td>Print &#8220;Not found&#8221;</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Print simple statistics on standard error</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Save found/not found positions on file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print this help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example:</p>
<div class="highlight-python"><pre>SearchInTable.py -i mytable.gz -j mypositions -s -o results</pre>
</div>
</div>
<div class="section" id="selectpositions-py">
<h2>selectPositions.py<a class="headerlink" href="#selectpositions-py" title="Permalink to this headline"></a></h2>
<p>selectPositions.py can filter an output REDItool table according to different criteria.</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>Table file from REDItools</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-d</span></kbd></td>
<td>Base distribution column for DNA-Seq (-1: no DNA-Seq) [-1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-c</span></kbd></td>
<td>Coverage RNA-Seq [5]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-C</span></kbd></td>
<td>Coverage DNA-Seq [5]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-v</span></kbd></td>
<td>Bases supporting RNA-Seq variation [1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-V</span></kbd></td>
<td>Bases supporting DNA-Seq variation [0]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Substitutions to select in RNA-Seq (separated by comma AG,CT) [all]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-f</span></kbd></td>
<td>Frequency of variation in RNA-Seq [0.1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-F</span></kbd></td>
<td>Frequency of non-variation in DNA-Seq [0.95]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-e</span></kbd></td>
<td>Exclude multiple substitutions in RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-r</span></kbd></td>
<td>Exclude invariant sites in RNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-R</span></kbd></td>
<td>Exclude variant sites in DNA-Seq #</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-u</span></kbd></td>
<td>Use only positions supported by DNA-Seq</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Save selected positions on outTable_533864766</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print this help</td></tr>
</tbody>
</table>
</dd>
</dl>
</div>
<div class="section" id="sorttable-py-new-in-version-1-0-3">
<h2>SortTable.py (new in version 1.0.3)<a class="headerlink" href="#sorttable-py-new-in-version-1-0-3" title="Permalink to this headline"></a></h2>
<p>SortTable.py is a facility to sort an input table according to genomic region and coordinates.
Input table may be TAB-delimited or use whatever delimiter. This script has been introduced for
users working on operating systems in which the sort program is not present. On very large input files
sorting time may increase respect to the unix sort program. On common GFFs including gene annotations
sorting time is lower or equal than the sort program.
Optionally, an input table can be outputted as TAB-delimited. Please tune the -b option according to
your memory capability. Default value should work well for many computers.</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>Table file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-d</span></kbd></td>
<td>Delimiter character [t] (default TAB)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Sequence name column [1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-c</span></kbd></td>
<td>Start column [4]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-e</span></kbd></td>
<td>End column (can be identical to -c) [5]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-m</span></kbd></td>
<td>Skip lines starting with [#]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Sorted output file [sortedTable_%s]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-O</span></kbd></td>
<td>Output as TAB-delimited</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-b</span></kbd></td>
<td>Buffer size (as number of lines) [32000]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-t</span></kbd></td>
<td>Temporary directory to use (multiple -t may be used)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print this help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example for input table with space-based columns:</p>
<div class="highlight-python"><pre>SortTable.py -i mytable.txt -d ' ' -s 1 -c 2 -e 2 -O</pre>
</div>
</div>
<div class="section" id="tabletotabix-py-new-in-version-1-0-3">
<h2>tableToTabix.py (new in version 1.0.3)<a class="headerlink" href="#tabletotabix-py-new-in-version-1-0-3" title="Permalink to this headline"></a></h2>
<p>tableToTabix.py create a tabix table and by default sort the input table. Useful to generate tabix tables
from GFF file to prepare in advance annotation tables for REDItools. As a specific alternative, GFFtoTabix.py may be used. Option -b should be tuned according to
memory capabilities in case of sorting. Tabix by default compresses the input table and then creates indices. If
a copy of the sorted and uncompressed table is needed, -u option should be used.</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>TAB-delimited file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Sequence name column [1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-c</span></kbd></td>
<td>Start column [4]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-e</span></kbd></td>
<td>End column (can be identical to -c) [5]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-m</span></kbd></td>
<td>Skip lines starting with [#]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-0</span></kbd></td>
<td>Zero based coordinates</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-S</span></kbd></td>
<td>Do not sort input file (sort by default)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-b</span></kbd></td>
<td>Buffer size (as number of lines) [32000]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-t</span></kbd></td>
<td>Temporary directory to use (multiple -t may be used)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-u</span></kbd></td>
<td>Save an uncompressed GFF copy (add _copy suffix)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print this help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example for GFF file:</p>
<div class="highlight-python"><pre>tableToTabix.py -i mytable.gff -u</pre>
</div>
</div>
<div class="section" id="sortgff-py-new-in-version-1-0-3">
<h2>SortGFF.py (new in version 1.0.3)<a class="headerlink" href="#sortgff-py-new-in-version-1-0-3" title="Permalink to this headline"></a></h2>
<p>SortGFF.py is a facility to sort only GFF files and works as SortTable.py. Option -b should be tuned
according to memory capabilities. Useful for users working on machines in which the sort program is not present.</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>GFF file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Sorted output file [GFF_sorted_%s]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-b</span></kbd></td>
<td>Buffer size (as number of lines) [32000]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-t</span></kbd></td>
<td>Temporary directory to use (multiple -t may be used)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print this help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example:</p>
<div class="highlight-python"><pre>SortGFF.py -i mytable.gff -u</pre>
</div>
</div>
<div class="section" id="gfftotabix-py-new-in-version-1-0-3">
<h2>GFFtoTabix.py (new in version 1.0.3)<a class="headerlink" href="#gfftotabix-py-new-in-version-1-0-3" title="Permalink to this headline"></a></h2>
<p>GFFtoTabix.py is a script to convert GFF files to tabix files. Useful to generate GFF annotations for REDItools.
It works as tableToTabix.py but specific for GFF files. Use -u to store an uncompressed version of input GFF.
Sorting is activated by default.</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>GFF file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-S</span></kbd></td>
<td>Do not sort GFF (sort by default)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-b</span></kbd></td>
<td>Buffer size (as number of lines) [32000]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-t</span></kbd></td>
<td>Temporary directory to use (multiple -t may be used)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-u</span></kbd></td>
<td>Save an uncompressed GFF copy (add _copy suffix)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print this help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example:</p>
<div class="highlight-python"><pre>GFFtoTabix.py -i mytable.gff -u</pre>
</div>
</div>
<div class="section" id="tabletogff-py-new-in-version-1-0-3">
<h2>TableToGFF.py (new in version 1.0.3)<a class="headerlink" href="#tabletogff-py-new-in-version-1-0-3" title="Permalink to this headline"></a></h2>
<p>TableToGFF.py is a facility to convert an output REDItool table to GFF in order to be used as input in REDItoolDnaRna.py.
Optionally, it can sort and tabix the output GFF. Options -b and -T work if -s is in effect.</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>Table file from REDItools</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Sort output GFF</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-t</span></kbd></td>
<td>Tabix output GFF (requires Pysam module)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-b</span></kbd></td>
<td>Buffer size (as number of lines) [32000] (requires -s)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-T</span></kbd></td>
<td>Temporary directory (requires -s)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Outfile [outTable_%s.gff]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print this help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example:</p>
<div class="highlight-python"><pre>TableToGFF.py -i mytable -s -t -o myoutputTable.gff</pre>
</div>
</div>
</div>
<div class="section" id="usage-example">
<h1>Usage Example<a class="headerlink" href="#usage-example" title="Permalink to this headline"></a></h1>
<div class="section" id="test-dataset">
<h2>Test dataset<a class="headerlink" href="#test-dataset" title="Permalink to this headline"></a></h2>
<p>The following example shows how to use REDItools to detect RNA editing sites in Alu regions. To this aim, please
download and extract testREDItools.tar.gz:</p>
<div class="highlight-python"><pre>gunzip testREDItools.tar.gz
tar -xvf testREDItools.tar
cd testREDItools
ls
dna.bam      reditool-output-sample.tar.gz  reference.fa.fai       refGene.sorted.gtf.gz.tbi  rmsk.gtf.gz.tbi  rna.bam.bai
dna.bam.bai  reference.fa            refGene.sorted.gtf.gz  rmsk.gtf.gz                rna.bam</pre>
</div>
<dl class="docutils">
<dt>The testREDItools folder contains:</dt>
<dd><ul class="first last simple">
<li>rna.bam: aligned RNA-Seq reads in BAM format (coordinate sorted by samtools sort)</li>
<li>rna.bam.bai: index file for rna.bam (indexed by samtools index)</li>
<li>dna.bam: aligned DNA-Seq reads (coordinate sorted by samtools sort)</li>
<li>dna.bam.bai: index file for dna.bam (indexed by samtools index)</li>
<li>reference.fa: reference genome in fasta format</li>
<li>reference.fa.fai: index file for reference.fa (indexed by samtools faidx)</li>
<li>rmsk.gtf.gz: gzipped gtf file of repeated elements (sorted by sort -k1,1 -k2,4n)</li>
<li>rmsk.gtf.gz.tbi: tabix index file for rmsk.gtf.gz</li>
<li>refGene.sorted.gtf.gz: gzipped gtf file of RefSeq genes (sorted by sort -k1,1 -k2,4n)</li>
<li>refGene.sorted.gtf.gz.tbi: tabix index file for refGene.sorted.gtf.gz</li>
<li>reditool-output-sample.tar.gz: folder including output examples</li>
</ul>
</dd>
</dl>
<p>RNA and DNA reads were extracted from the region chr21:47721047-47743785.
Reads were obtained according to Ramaswami et al. paper. RNA reads mapped by BWA where kindly provided by
Gokul Ramaswami.
If REDItools have been correctly installed in your path, you can call all REDItoolDnaRNA.py options:</p>
<div class="highlight-python"><div class="highlight"><pre><span class="n">REDItoolDnaRna</span><span class="o">.</span><span class="n">py</span> <span class="o">-</span><span class="n">h</span>
</pre></div>
</div>
<p>If everything is ok, REDItoolDnaRna.py can be launched by:</p>
<div class="highlight-python"><pre>[epicardi@srv00 sample]$ REDItoolDnaRna.py -i rna.bam -j dna.bam -f reference.fa -o reditool-test -c 10,1 -q 25,25 -m 20,20 -s 2 -g 1 -u -a 6-0 -v 2 -n0.0 -N0.0 -V
Script time --&gt; START: 06/02/2013 20:44:48
Analysis ID: 891206177
Analysis on 1 regions.
Started analysis on region: chr21
Job completed for region: chr21
Merging Tables.
Results saved on reditool-test/DnaRna_891206177/outTable_891206177
Script time --&gt; END: 06/02/2013 20:45:13</pre>
</div>
<p>In this case we are requiring to extract RNA and DNA positions with a minimal coverage of 10 for DNA and 1 for RNA, a minimum quality score of 25 for both,
a minimum mapping quality of 20 for both. Please consider that quality scores in DNA are in Sanger format, so 33 is used as offset to calculate the phred score
per site genomic. On the contrary, quality scores for RNA reads are in the Illumina 1.3+ format and, thus, the value 64 has to be used through the option -Q.
Since RNA reads are strand oriented and the second read of the pair maintains the RNA orientation, we ask REDItoolDnaRna.py to infer the strand by -s 2 (second in pair good for orientation),
and assign to each position the overrepresented strand (-g 1 option). In addition, we remove first 6 nucleotides from each read (-a 6-0) and require sites supported by at least 2
variant bases (-v 2) without taking into account the frequency of variation in both RNA (-n 0.0) and DNA (-N 0.0). Finally, we exclude positions not supported by DNA-Seq reads (-V).
At the end of the run, REDItoolDnaRna.py will create the folder reditool-test/DnaRna_891206177/ and the output table outTable_891206177. The folder will contain also the parameters.txt
file summarizing all parameter values used for the run.
Now enter into the output folder:</p>
<div class="highlight-python"><pre>cd reditool-test/DnaRna_891206177/
ls
outTable_891206177  parameters.txt</pre>
</div>
<p>and run the accessory script selectPositions.py to filter out other positions:</p>
<div class="highlight-python"><pre>selectPositions.py -i outTable_891206177 -d 12 -c 2 -C 10 -v 2 -V 0 -f 0.1 -F 1.0 -e -u -o candidates.txt
Script time --&gt; START: 06/02/2013 20:46:20
Reading table...
Total lines: 17879
Filtered in lines: 10
Selected lines saved on candidates.txt
Script time --&gt; END: 06/02/2013 20:46:21</pre>
</div>
<p>In this case we are excluding positions showing an editing frequency lower that 0.1 and supported by heterozigous DNA sites.
Option -u is used to collect only positions supported by DNA-Seq. Option -e is used to remove RNA positions showing multiple substitutions.
According to these criteria only 10 positions will be retained:</p>
<div class="highlight-python"><pre>more candidates.txt
Region  Position        Reference       Strand  Coverage-q25    MeanQ   BaseCount[A,C,G,T]      AllSubs Frequency       gCoverage-q25   gMeanQ  gBaseCount[A,C,G,T]  gAllSubs        gFrequency
chr21   47739131        A       0       20      34.60   [18, 0, 2, 0]   AG      0.10    30      30.40   [30, 0, 0, 0]   -       0.00
chr21   47739578        A       0       14      38.50   [11, 0, 3, 0]   AG      0.21    26      30.27   [26, 0, 0, 0]   -       0.00
chr21   47739644        A       0       18      36.61   [13, 0, 5, 0]   AG      0.28    23      30.22   [23, 0, 0, 0]   -       0.00
chr21   47739647        A       0       16      36.12   [7, 0, 9, 0]    AG      0.56    18      30.67   [18, 0, 0, 0]   -       0.00
chr21   47739724        A       0       8       37.00   [6, 0, 2, 0]    AG      0.25    30      30.10   [30, 0, 0, 0]   -       0.00
chr21   47739725        A       0       8       37.75   [5, 0, 3, 0]    AG      0.38    30      29.93   [30, 0, 0, 0]   -       0.00
chr21   47739764        A       0       6       37.33   [4, 0, 2, 0]    AG      0.33    19      30.37   [19, 0, 0, 0]   -       0.00
chr21   47740295        A       0       10      34.00   [6, 0, 4, 0]    AG      0.40    30      29.77   [30, 0, 0, 0]   -       0.00
chr21   47741150        A       0       28      36.57   [25, 0, 3, 0]   AG      0.11    24      31.54   [24, 0, 0, 0]   -       0.00
chr21   47741221        A       0       49      36.33   [44, 0, 5, 0]   AG      0.10    40      29.50   [40, 0, 0, 0]   -       0.00</pre>
</div>
<p>Now we can annotate these positions using the information in the Repeat Mask database to look at Alu sites:</p>
<div class="highlight-python"><pre>AnnotateTable.py -a ../../rmsk.gtf.gz -i candidates.txt -u -c 1,2,3 -n RepMask -o candidates.rmsk.txt
Script time --&gt; START: 06/02/2013 20:48:29
Table saved on candidates.rmsk.txt
Script time --&gt; END: 06/02/2013 20:48:29

more candidates.rmsk.txt
Region  Position        Reference       Strand  Coverage-q25    MeanQ   BaseCount[A,C,G,T]      AllSubs Frequency       gCoverage-q25   gMeanQ  gBaseCount[A,C,G,T]  gAllSubs        gFrequency      RepMask_feat    RepMask_gid     RepMask_tid
chr21   47739131        A       0       20      34.60   [18, 0, 2, 0]   AG      0.10    30      30.40   [30, 0, 0, 0]   -       0.00    -       -       -
chr21   47739578        A       0       14      38.50   [11, 0, 3, 0]   AG      0.21    26      30.27   [26, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
chr21   47739644        A       0       18      36.61   [13, 0, 5, 0]   AG      0.28    23      30.22   [23, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
chr21   47739647        A       0       16      36.12   [7, 0, 9, 0]    AG      0.56    18      30.67   [18, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
chr21   47739724        A       0       8       37.00   [6, 0, 2, 0]    AG      0.25    30      30.10   [30, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
chr21   47739725        A       0       8       37.75   [5, 0, 3, 0]    AG      0.38    30      29.93   [30, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
chr21   47739764        A       0       6       37.33   [4, 0, 2, 0]    AG      0.33    19      30.37   [19, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE
chr21   47740295        A       0       10      34.00   [6, 0, 4, 0]    AG      0.40    30      29.77   [30, 0, 0, 0]   -       0.00    SINE    AluSz   Alu-SINE
chr21   47741150        A       0       28      36.57   [25, 0, 3, 0]   AG      0.11    24      31.54   [24, 0, 0, 0]   -       0.00    SINE    AluSx4  Alu-SINE
chr21   47741221        A       0       49      36.33   [44, 0, 5, 0]   AG      0.10    40      29.50   [40, 0, 0, 0]   -       0.00    SINE    AluSx4  Alu-SINE</pre>
</div>
<p>To remove positions not annotated in SINE regions we use the FilterTable.py script:</p>
<div class="highlight-python"><pre>FilterTable.py -i candidates.rmsk.txt -f ../../rmsk.gtf.gz -F SINE -E -o candidates.rmsk.alu.txt -p
Script time --&gt; START: 06/02/2013 20:50:25
Reading Table file...
All positions: 10
Positions filtered in: 9
Positions filtered out: 1
Script time --&gt; END: 06/02/2013 20:50:25</pre>
</div>
<p>Finally we can add gene annotations using RefSeq database:</p>
<div class="highlight-python"><pre>AnnotateTable.py -a /home/epicardi/annotation/hg19/annotation/tabixAnn/refseq/refGene.sorted.gtf.gz -i candidates.rmsk.alu.txt -u -c 1,2 -n RefSeq -o candidates.rmsk.alu.ann.txt
Script time --&gt; START: 06/02/2013 20:52:48
Table saved on candidates.rmsk.alu.ann.txt
Script time --&gt; END: 06/02/2013 20:52:48

more candidates.rmsk.alu.ann.txt
Region  Position        Reference       Strand  Coverage-q25    MeanQ   BaseCount[A,C,G,T]      AllSubs Frequency       gCoverage-q25   gMeanQ  gBaseCount[A,C,G,T]  gAllSubs        gFrequency      RepMask_feat    RepMask_gid     RepMask_tid     RefSeq_feat     RefSeq_gid
chr21   47739578        A       0       14      38.50   [11, 0, 3, 0]   AG      0.21    26      30.27   [26, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
chr21   47739644        A       0       18      36.61   [13, 0, 5, 0]   AG      0.28    23      30.22   [23, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
chr21   47739647        A       0       16      36.12   [7, 0, 9, 0]    AG      0.56    18      30.67   [18, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
chr21   47739724        A       0       8       37.00   [6, 0, 2, 0]    AG      0.25    30      30.10   [30, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
chr21   47739725        A       0       8       37.75   [5, 0, 3, 0]    AG      0.38    30      29.93   [30, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
chr21   47739764        A       0       6       37.33   [4, 0, 2, 0]    AG      0.33    19      30.37   [19, 0, 0, 0]   -       0.00    SINE    AluSq4  Alu-SINE     intron  C21orf58
chr21   47740295        A       0       10      34.00   [6, 0, 4, 0]    AG      0.40    30      29.77   [30, 0, 0, 0]   -       0.00    SINE    AluSz   Alu-SINE     intron  C21orf58
chr21   47741150        A       0       28      36.57   [25, 0, 3, 0]   AG      0.11    24      31.54   [24, 0, 0, 0]   -       0.00    SINE    AluSx4  Alu-SINE     intron  C21orf58
chr21   47741221        A       0       49      36.33   [44, 0, 5, 0]   AG      0.10    40      29.50   [40, 0, 0, 0]   -       0.00    SINE    AluSx4  Alu-SINE     intron  C21orf58</pre>
</div>
</div>
</div>
<div class="section" id="contact">
<h1>Contact<a class="headerlink" href="#contact" title="Permalink to this headline"></a></h1>
<ul class="simple">
<li><strong>Ernesto Picardi</strong>: ernesto.picardi@gmail.com</li>
</ul>
</div>
<a href="#">Go to top of page</a>
  </body>
</html>
