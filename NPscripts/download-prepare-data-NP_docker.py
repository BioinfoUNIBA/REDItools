#!/usr/bin/env python 
import sys, os, time
import commands
import distutils.spawn

try:
	wdir=sys.argv[1] # working directory
	redipath=sys.argv[2] # path to REDItools folder
	usepath=sys.argv[3]
except:
	sys.exit('<working directory> <path to REDItools folder> <use file with paths 0/1>')

def getData(cmd):
	tr=0
	while 1:
		st,out=commands.getstatusoutput(cmd)
		if st==0:
			return 0
		tr+=1
		if tr==10: break
	if tr>0: return 1

def is_tool(name):
	wn=distutils.spawn.find_executable(name)
	if wn==None: return 1
	else: return wn

def get_time(tstart,tend):
	telapsed=tend - tstart
	t_taken=time.strftime("%H:%M:%S", time.gmtime(telapsed))
	return t_taken

if usepath!='1':
	exe=['bwa','STAR','awk','bgzip','tabix','sort','gtf_splicesites','wget -c --retry-connrefused --tries=0 --timeout=5 ','python','gunzip']
	nt=[]
	prg={}
	for i in exe:
		p=is_tool(i)
		if p==1: nt.append(i)
		prg[i]=p
	if len(nt)>0:
		for i in nt:
			sys.stderr.write('Program %s NOT found\n' %(i))
		sys.exit('Install required software first.')
else:
	if not os.path.exists('mypaths'): sys.exit('File mypaths does not exists.')
	nt=[]
	f=open('mypaths')
	prg={}
	for i in f:
		l=(i.strip()).replace(' ','')
		l=l.split('=')
		prg[l[0]]=l[1]	

if not os.path.exists(redipath): sys.exit('REDItools path does not exist.')
redirec=os.path.join(redipath,'accessory','rediportal2recoding.py')
if not os.path.exists(redipath): sys.exit('rediportal2recoding.py script not found.')
prg['redirec']=redirec

cdir=os.getcwd()
sys.stderr.write('Current directory: %s\n' %(cdir))
folder=os.path.join(cdir,wdir)
if not os.path.exists(folder):
	os.mkdir(folder)
	sys.stderr.write('Directory %s created.\n' %(wdir))
else:
	sys.stderr.write('Found working directory.\n')
sys.stderr.write('Entering %s\n' %(wdir))
os.chdir(folder)

#human genome
sys.stderr.write('Getting human genome\n')
tstart = time.time()
os.mkdir('genome_hg19')
os.chdir('genome_hg19')
wcmd='%s ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz' %(prg['wget'])
ot=getData(wcmd)
if ot==1: sys.stderr.write('I cannot download the human genome.\n')
else: sys.stderr.write('Human genome complete.\n')
tend = time.time()
sys.stderr.write('Human genome - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')
#Gencode
sys.stderr.write('Getting GENCODE genes\n')
tstart = time.time()
os.mkdir('Gencode_annotation')
os.chdir('Gencode_annotation')
gcmd='%s ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/gencode.v30lift37.annotation.gtf.gz' %(prg['wget'])
ot=getData(gcmd)
if ot==1: sys.stderr.write('I cannot download GENCODE annotations.\n')
else: sys.stderr.write('GENCODE annotations ready.\n')
tend = time.time()
sys.stderr.write('GENCODE annotations - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')
#RefSeq
sys.stderr.write('Getting RefSeq hg19\n')
tstart = time.time()
os.mkdir('Strand_detection')
os.chdir('Strand_detection')
gcmd='%s --no-check-certificate https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg19_RefSeq.bed.gz' %(prg['wget'])
ot=getData(gcmd)
if ot==1: sys.stderr.write('I cannot download REFSEQ annotations.\n')
else: sys.stderr.write('REFSEQ annotations ready.\n')
tend = time.time()
sys.stderr.write('REFSEQ annotations - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')
#RepeatMasker
sys.stderr.write('Getting RepeatMasker\n')
tstart = time.time()
os.mkdir('rmsk')
os.chdir('rmsk')
gcmd='%s http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz' %(prg['wget'])
ot=getData(gcmd)
if ot==1: sys.stderr.write('I cannot download RepeatMasker annotations.\n')
else: sys.stderr.write('RepeatMasker annotations ready.\n')
tend = time.time()
sys.stderr.write('RepeatMasker annotations - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')
#dbSNP
sys.stderr.write('Getting dbSNP\n')
tstart = time.time()
os.mkdir('snp151')
os.chdir('snp151')
gcmd='%s http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151.txt.gz' %(prg['wget'])
ot=getData(gcmd)
if ot==1: sys.stderr.write('I cannot download dbSNP annotations.\n')
else: sys.stderr.write('dbSNP annotations ready.\n')
tend = time.time()
sys.stderr.write('dbSNP annotations - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')
#REDIportal
sys.stderr.write('Getting REDIportal\n')
tstart = time.time()
os.mkdir('rediportal')
os.chdir('rediportal')
gcmd='%s http://srv00.recas.ba.infn.it/webshare/rediportalDownload/table1_full.txt.gz' %(prg['wget'])
ot=getData(gcmd)
if ot==1: sys.stderr.write('I cannot download REDIportal annotations.\n')
else: sys.stderr.write('REDIportal annotations ready.\n')
tend = time.time()
sys.stderr.write('REDIportal annotations - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')
#NA12878 - WGS
sys.stderr.write('Getting NA12878 data - WGS\n')
tstart = time.time()
os.mkdir('WGS_ERR262997')
os.chdir('WGS_ERR262997')
fq1cmd='%s ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR262/ERR262997/ERR262997_1.fastq.gz' %(prg['wget'])
fq2cmd='%s ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR262/ERR262997/ERR262997_2.fastq.gz' %(prg['wget'])
f1=getData(fq1cmd)
f2=getData(fq2cmd)
if f1==1: sys.stderr.write('I cannot download READ1.\n')
else:
	gu=getData('%s ERR262997_1.fastq.gz' %(prg['gunzip']))
	sys.stderr.write('READ1 ready.\n')
if f2==1: sys.stderr.write('I cannot download READ2.\n')
else:
	gu=getData('%s ERR262997_2.fastq.gz' %(prg['gunzip']))
	sys.stderr.write('READ2 ready.\n')
tend = time.time()
sys.stderr.write('NA12878 data - WGS - time taken: %s\n' %(get_time(tstart,tend)))	
os.chdir('..')

#NA12878 - RNAseq
sys.stderr.write('Getting NA12878 data - RNAseq\n')
tstart = time.time()
os.mkdir('RNASeq_SRR1258218')
os.chdir('RNASeq_SRR1258218')
fq1cmd='%s ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/008/SRR1258218/SRR1258218_1.fastq.gz' %(prg['wget'])
fq2cmd='%s ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/008/SRR1258218/SRR1258218_2.fastq.gz' %(prg['wget'])
f1=getData(fq1cmd)
f2=getData(fq2cmd)
if f1==1: sys.stderr.write('I cannot download READ1.\n')
else:
	#gu=getData('%s SRR1258218_1.fastq.gz' %(prg['gunzip']))
	sys.stderr.write('READ1 ready.\n')
if f2==1: sys.stderr.write('I cannot download READ2.\n')
else:
	#gu=getData('%s SRR1258218_2.fastq.gz' %(prg['gunzip']))
	sys.stderr.write('READ2 ready.\n')
tend = time.time()
sys.stderr.write('NA12878 data - RNAseq - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')

#PRJNA316625
sys.stderr.write('Getting PRJNA316625 data\n')
tstart = time.time()
os.mkdir('PRJNA_316625')
os.chdir('PRJNA_316625')
fqlist=['ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/003/SRR3306823/SRR3306823_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/003/SRR3306823/SRR3306823_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/004/SRR3306824/SRR3306824_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/004/SRR3306824/SRR3306824_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/005/SRR3306825/SRR3306825_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/005/SRR3306825/SRR3306825_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/006/SRR3306826/SRR3306826_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/006/SRR3306826/SRR3306826_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/007/SRR3306827/SRR3306827_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/007/SRR3306827/SRR3306827_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/008/SRR3306828/SRR3306828_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/008/SRR3306828/SRR3306828_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/009/SRR3306829/SRR3306829_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/009/SRR3306829/SRR3306829_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/000/SRR3306830/SRR3306830_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/000/SRR3306830/SRR3306830_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/001/SRR3306831/SRR3306831_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/001/SRR3306831/SRR3306831_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/002/SRR3306832/SRR3306832_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/002/SRR3306832/SRR3306832_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/003/SRR3306833/SRR3306833_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/003/SRR3306833/SRR3306833_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/004/SRR3306834/SRR3306834_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/004/SRR3306834/SRR3306834_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/005/SRR3306835/SRR3306835_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/005/SRR3306835/SRR3306835_2.fastq.gz','ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/006/SRR3306836/SRR3306836_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR330/006/SRR3306836/SRR3306836_2.fastq.gz']
for i in fqlist:
	fq1,fq2=i.split(',')
	base=(os.path.basename(fq1)).split("_")[0]
	os.mkdir(base)
	os.chdir(base)
	fq1cmd='%s %s' %(prg['wget'],fq1)
	fq2cmd='%s %s' %(prg['wget'],fq2)
	fq1_=getData(fq1cmd)
	fq2_=getData(fq2cmd)
	if fq1_+fq2_>0:
		sys.stderr.write('I cannot download all files in %s.\n' %(base))
		os.chdir('..')
	else:
		#gu1=getData('%s %s' %(prg['gunzip'],os.path.basename(fq1))) 
		#gu2=getData('%s %s' %(prg['gunzip'],os.path.basename(fq2)))
		sys.stderr.write('Files in %s ready.\n' %(base))
		os.chdir('..')
tend = time.time()
sys.stderr.write('PRJNA316625 - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')

sys.stderr.write('Preparing data ...\n')
sys.stderr.write('BWA indexing...\n')
tstart = time.time()
os.chdir('genome_hg19')
cmd='%s GRCh37.primary_assembly.genome.fa.gz' %(prg['gunzip'])
cmd1='%s index GRCh37.primary_assembly.genome.fa' %(prg['bwa'])
ot=getData(cmd)
ot1=getData(cmd1)
if ot+ot1>0: sys.stderr.write('BWA indexing error.\n')
else: sys.stderr.write('BWA indices ready.\n')
tend = time.time()
sys.stderr.write('BWA indexing - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')

sys.stderr.write('STAR indexing...\n')
cmd='%s Gencode_annotation/gencode.v30lift37.annotation.gtf.gz' %(prg['gunzip'])
ot=getData(cmd)
if ot==1: sys.stderr.write('Gunzipping gencode error.\n')
else: sys.stderr.write('Gunzipping gencode ready.\n')
tstart = time.time()
if not os.path.exists('STAR'): os.mkdir('STAR')
os.chdir('STAR')
os.mkdir('STAR_genome_index_ucsc')
cmd='%s --runMode genomeGenerate --genomeDir STAR_genome_index_ucsc --genomeFastaFiles ../genome_hg19/GRCh37.primary_assembly.genome.fa --sjdbGTFfile ../Gencode_annotation/gencode.v30lift37.annotation.gtf --sjdbOverhang 75' %(prg['STAR'])
ot=getData(cmd)
if ot==1: sys.stderr.write('STAR indexing error.\n')
else: sys.stderr.write('STAR indices ready.\n')
tend = time.time()
sys.stderr.write('STAR indexing - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')

sys.stderr.write('Prepare RepeatMasker annotations ...\n')
tstart = time.time()
os.chdir('rmsk')
cmd4='%s rmsk.txt.gz' %(prg['gunzip'])
cmd='%s \'OFS="\t"{print $6,"rmsk_hg19",$12,$7+1,$8,".",$10,".","gene_id \""$11"\"; transcript_id \""$13"\";"}\' rmsk.txt > rmsk.gtf' %(prg['awk'])
cmd1='%s -k1,1 -k4,4n rmsk.gtf > rmsk.sorted.gtf' %(prg['sort'])
cmd2='%s rmsk.sorted.gtf' %(prg['bgzip'])
cmd3='%s -p gff rmsk.sorted.gtf.gz' %(prg['tabix'])
ot4=getData(cmd4)
ot=getData(cmd)
ot1=getData(cmd1)
ot2=getData(cmd2)
ot3=getData(cmd3)
if ot==4: sys.stderr.write('RepeatMasker gunzip error.\n')
if ot==1: sys.stderr.write('RepeatMasker awk error.\n')
if ot1==1: sys.stderr.write('RepeatMasker sort error.\n')
if ot2==1: sys.stderr.write('RepeatMasker bgzip error.\n')
if ot3==1: sys.stderr.write('RepeatMasker tabix error.\n')
if ot+ot1+ot2+ot3+ot4==0: sys.stderr.write('RepeatMasker ready.\n')
tend = time.time()
sys.stderr.write('Prepare RepeatMasker annotations - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')

sys.stderr.write('Prepare dbSNP annotations ...\n')
tstart = time.time()
os.chdir('snp151')
cmd4='%s snp151.txt.gz' %(prg['gunzip'])
cmd='%s \'OFS="\t"{if ($11=="genomic" && $12=="single") print $2,"ucsc_snp151_hg19","snp",$4,$4,".",$7,".","gene_id \""$5"\"; transcript_id \""$5"\";"}\' snp151.txt > snp151.gtf' %(prg['awk'])
cmd1='%s -k1,1 -k4,4n snp151.gtf > snp151.sorted.gtf' %(prg['sort'])
cmd2='%s snp151.sorted.gtf' %(prg['bgzip'])
cmd3='%s -p gff snp151.sorted.gtf.gz' %(prg['tabix'])
ot4=getData(cmd4)
ot=getData(cmd)
ot1=getData(cmd1)
ot2=getData(cmd2)
ot3=getData(cmd3)
if ot==4: sys.stderr.write('dbSNP gunzip error.\n')
if ot==1: sys.stderr.write('dbSNP awk error.\n')
if ot1==1: sys.stderr.write('dbSNP sort error.\n')
if ot2==1: sys.stderr.write('dbSNP bgzip error.\n')
if ot3==1: sys.stderr.write('dbSNP tabix error.\n')
if ot+ot1+ot2+ot3+ot4==0: sys.stderr.write('dbSNP ready.\n')
tend = time.time()
sys.stderr.write('Prepare dbSNP annotations - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')

sys.stderr.write('Prepare splice sites annotations ...\n')
tstart = time.time()
os.chdir('Gencode_annotation')
cmd='%s gencode.v30lift37.annotation.gtf > splicesites' %(prg['gtf_splicesites'])
cmd1='%s -F" " \'{split($2,a,":"); split(a[2],b,"."); if (b[1]>b[3]) print a[1],b[3],b[1],toupper(substr($3,1,1)),"-"; else print a[1],b[1],b[3],toupper(substr($3,1,1)),"+"}\' splicesites > gencode.v30lift37.splicesites.txt' %(prg['awk'])
ot=getData(cmd)
ot1=getData(cmd1)
if ot==1: sys.stderr.write('Splice sites gtf_splicesites error.\n')
if ot1==1: sys.stderr.write('Splice sites sort error.\n')
if ot+ot1==0: sys.stderr.write('Splice sites ready.\n')
tend = time.time()
sys.stderr.write('Prepare splice sites annotations - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')

sys.stderr.write('Prepare REDIportal annotations ...\n')
tstart = time.time()
os.chdir('rediportal')
cmd7='%s table1_full.txt.gz' %(prg['gunzip'])
cmd='%s \'OFS="\t"{sum+=1; print $1,"rediportal","ed",$2,$2,".",$5,".","gene_id \""sum"\"; transcript_id \""sum"\";"}\' table1_full.txt > atlas.gtf' %(prg['awk'])
cmd1='%s atlas.gtf' %(prg['bgzip'])
cmd2='%s -p gff atlas.gtf.gz' %(prg['tabix'])
cmd3='%s %s table1_full.txt > atlas_recoding.gff' %(prg['python'],redirec)
cmd4='%s -V -k1,1 -k4,4n atlas_recoding.gff > srtd_atlas_recoding.gff' %(prg['sort'])
cmd5='%s srtd_atlas_recoding.gff' %(prg['bgzip'])
cmd6='%s -p gff srtd_atlas_recoding.gff.gz' %(prg['tabix'])
ot7=getData(cmd7)
ot=getData(cmd)
ot1=getData(cmd1)
ot2=getData(cmd2)
ot3=getData(cmd3)
ot4=getData(cmd4)
ot5=getData(cmd5)
ot6=getData(cmd6)
if ot==7: sys.stderr.write('REDIportal gunzip error.\n')
if ot==1: sys.stderr.write('REDIportal awk error.\n')
if ot1==1: sys.stderr.write('REDIportal bgzip error.\n')
if ot2==1: sys.stderr.write('REDIportal tabix error.\n')
if ot3==1: sys.stderr.write('REDIportal python error.\n')
if ot4==1: sys.stderr.write('REDIportal sort error.\n')
if ot5==1: sys.stderr.write('REDIportal bgzip error.\n')
if ot6==1: sys.stderr.write('REDIportal tabix error.\n')
if ot+ot1+ot2+ot3+ot4+ot5+ot6+ot7==0: sys.stderr.write('REDIportal ready.\n')
tend = time.time()
sys.stderr.write('Prepare REDIportal annotations - time taken: %s\n' %(get_time(tstart,tend)))
os.chdir('..')

sys.stderr.write('ALL DONE. ENJOY REDItools.\n')
