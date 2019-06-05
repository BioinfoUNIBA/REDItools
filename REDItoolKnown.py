#!/home/epicardi/bin/python27/bin/python
# Copyright (c) 2013-2014 Ernesto Picardi <ernesto.picardi@uniba.it>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import sys, os, time, math, random, getopt, operator, string, errno
try: import pysam
except: sys.exit('Pysam module not found.')
from multiprocessing import Process, Queue
from Queue import Empty

pysamVersion=pysam.__version__

sys.stderr.write('Pysam version used: %s\n' %(pysamVersion))

version='1.1'

pid=str(os.getpid()+random.randint(0,999999999))

def usage():
	print """
USAGE: python REDItoolKnown.py [options]
Options:
-i		BAM file
-I		Sort input BAM file
-f		Reference in fasta file
-l		List of known RNA editing events
-C		Base interval to explore [100000]
-k		List of chromosomes to skip separated by comma or file
-t		Number of threads [1]
-o		Output folder [rediFolder_%s]
-F		Internal folder name [null]
-c		Min. read coverage [10]
-q		Min. quality score [30]
-m		Min. mapping quality score [30]*
-O		Min. homoplymeric length [5]
-s		Infer strand (for strand oriented reads) [1]
-g		Strand inference type 1:maxValue 2:useConfidence [1]
-x		Strand confidence [0.70]
-S		Strand correction
-G		Infer strand by gff annotation (must be sorted, otherwise use -X)
-X		Sort annotation files
-K		File with positions to exclude
-e		Exclude multi hits
-d		Exclude duplicates
-p		Use paired concardant reads only
-u		Consider mapping quality
-T		Trim x bases up and y bases down per read [0-0]
-B		Blat folder for correction
-U		Remove substitutions in homopolymeric regions
-v		Min. num. of reads supporting the variation [3]
-n		Min. editing frequency [0.1]
-E		Exclude positions with multiple changes
-P		File containing splice sites annotations
-r		Num. of bases near splice sites to explore [4]
-H		No Table Header
-h		Print this help

*This value may change according to the aligner:
	- For Bowtie use 255
	- For Bowtie2 use 40
	- For BWA use 30
	- For RNA-STAR use 255
	- For HiSAT2 use 60
	- For Tophat1 use 255
	- For Tophat2 use 50
	- For GSNAP use 30

"""%(pid)

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:f:k:t:o:c:q:m:O:s:edpuT:B:Sv:n:EP:r:hHIXG:K:l:C:F:x:g:U")
except getopt.GetoptError as err:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

MAX_DEPTH=100000
corrstr=0
strconf=0.70 #confidenza strand
useconf=0
bamfile=''
fastafile=''
sortbam=0
kfile=''
nochrs=[]
NCPU=1
infolder=''
outfolder_='rediFolder_%s' %(pid)
MINCOV=10
QVAL=33 #NOT USED
MQUAL=30
MAPQ=30
homo=5
rmpv = '0-0'
rmp = [int(x) for x in rmpv.split('-')]
getstrand=0 # considera la strand
exh=0 # escludi multi hits
exd=0 # escludi duplicati
conc=0 # se presenti paired-end, usa solo quelle concordanti
mq=0 # considera il map quality
rmnuc=0 # rimuovi nucleotide a monte ed a valle delle read; connesso a rmp e rmpv
blatr=0 # applica la correzione blat
blatfolder=''
rmsh=0 # rimuovi sostituzioni in omopolimeri di lunghezza maggiore o uguale a homo
vnuc=3 # numero minimo di basi che supportano la variazione
mmf=0.1 # frequenza minima della variazione
exms=0 # escludi sostituzioni multiple
exss=0 # escludi posizioni introniche nei pressi dei siti di splicing a nss nucleotidi 
nss=4 # basi introniche da esplorare per ogni sito si splicing
splicefile='' #'splicesites.hg18.sorted.txt'
#custsub=0  # use custom distribution 
#custfile='' # custom distribution file
#sigsites=0 # select significant sites
#test = 'bh' # select statistical test
usubs=[x+y for x in 'ACGT' for y in 'ACGT' if x!=y] # use these substitutions [default all]
#sval=0.05 # significant value
annfile='' # use annotation file for strand correction and features
sortann=0 # sort annotation file
uann=0 # use annotation
exfile='' # use annotations to exclude positions
expos=0 #
chunckval=100000
unchange1=1
unchange2=0
noheader=0

for o, a in opts:
	if o in ("-h","--help"):
		usage()
		sys.exit()
	elif o == "-H": noheader=1
	elif o == "-i": bamfile=a
	elif o == "-f": fastafile=a
	elif o == "-l": kfile=a
	elif o == "-k":
		if os.path.exists(a):
			f=open(a)
			nochrs=[x.strip() for x in f if x.strip()!='']
			f.close()
		else: nochrs=[x for x in a.split(',') if x.strip()!='']
	elif o == "-t": NCPU=int(a)
	elif o == "-F": infolder=a	
	elif o == "-o": outfolder_=a
	elif o == "-c": MINCOV=int(a)
	elif o == "-Q": QVAL=int(a)
	elif o == "-q": MQUAL=int(a)
	elif o == "-m": MAPQ=int(a)	
	elif o == "-O": homo=int(a)
	elif o == "-x": strconf=float(a)
	elif o == "-g":
		if a=='2': useconf=1
	elif o == "-s":
		getstrand=1
		if int(a)==1: unchange1,unchange2=1,0
		elif int(a)==0: unchange1,unchange2=0,0
		elif int(a)==2: unchange1,unchange2=0,1
		elif int(a)==12: unchange1,unchange2=1,1
	elif o == "-U": usubs=[x.upper() for x in a.split(',') if a.strip()!='']
	elif o == "-e": exh=1
	elif o == "-d": exd=1
	elif o == "-p": conc=1
	elif o == "-I": sortbam=1
	elif o == "-X": sortann=1
	elif o == "-C": chunckval=int(a)
	elif o == "-u": mq=1
	elif o == "-T":
		rmpv = a
		try:
			rmp = [int(x) for x in rmpv.split('-')]
			rmnuc=1
		except: rmnuc=0
	elif o == "-B":
		blatfolder=a
		if os.path.exists(blatfolder): blatr=1
	elif o == "-S": corrstr=1
	elif o == "-U": rmsh=1
	elif o == "-v": vnuc=int(a)
	elif o == "-n": mmf=float(a)
	elif o == "-E": exms=1
	elif o == "-P":
		splicefile=a
		if os.path.exists(splicefile): exss=1
	elif o == "-K":
		exfile=a
		if os.path.exists(exfile): expos=1	
	elif o == "-r": nss=int(a)
	elif o == "-G":
		annfile=a
		uann=1
	else:
		assert False, "Unhandled Option"

#######
commandLine=' '.join(sys.argv[1:])
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
params=[]
#Input parameters
params.append('REDItoolKnown version %s\n' %(version))
params.append('User command line: %s\n' %(commandLine))
params.append('Analysis ID: %s\n' %(pid))
params.append('Analysis time: %s\n' %(script_time))
params.append('-i --> BAM file: %s\n' %(bamfile))
params.append('-f --> Reference file: %s\n' %(fastafile))
params.append('-I --> Sort input BAM file: %i\n' %(sortbam))
params.append('-l --> File with known RNA editing positions: %s\n' %(kfile))
params.append('-X --> Sort annotation files: %i\n' %(sortann))
params.append('-k --> Regions to exclude: %s\n' %(','.join(nochrs)))
params.append('-t --> Number of working threads: %i\n' %(NCPU))
params.append('-C --> Base interval to explore: %i\n' %(chunckval))
params.append('-o --> Output folder: %s\n' %(outfolder_))
params.append('-F --> Infolder folder: %s\n' %(infolder))
params.append('-c --> Min. per base coverage: %i\n' %(MINCOV))
#params.append('-Q --> FastQ offset value: %i\n' %(QVAL))
params.append('-q --> Min. per base quality: %i\n' %(MQUAL))
params.append('-m --> Min. mapping quality: %i\n' %(MAPQ))
params.append('-O --> Min. homoplymeric length: %i\n' %(homo))
params.append('-s --> Infer strand: %i - %i-%i\n' %(getstrand,unchange1,unchange2))
params.append('-g --> Use confidence: %i\n' %(useconf))
params.append('-x --> Strand confidence: %.2f\n' %(strconf))
params.append('-S --> Strand correction : %i\n' %(corrstr))
params.append('-G --> GFF annotation to infer strand: %s\n' %(annfile))
params.append('-e --> Exclude multi hits: %i\n' %(exh))
params.append('-d --> Exclude duplicates: %i\n' %(exd))
params.append('-p --> Use paired concardant reads only: %i\n' %(conc))
params.append('-u --> Consider mapping quality: %i\n' %(mq))
params.append('-T --> Trim x bases up and y bases down per read: %i - %i-%i\n' %(rmnuc,rmp[0],rmp[1]))
params.append('-B --> Blat folder for correction: %s\n' %(blatfolder))
params.append('-S --> Remove substitutions in homopolymeric regions: %i\n' %(rmsh))
params.append('-v --> Min. num. of reads supporting the variation: %i\n' %(vnuc))
params.append('-n --> Min. editing frequency: %.2f\n' %(mmf))
params.append('-E --> Exclude positions with multiple changes: %i\n' %(exms))
params.append('-P --> File containing splice sites annotations: %s\n' %(splicefile))
params.append('-r --> Num. of bases near splice sites to explore: %i\n' %(nss))
params.append('-K --> File with positions to exclude: %s\n' %(exfile))
#######

def pid_exists(pid):
    """Check whether pid exists in the current process table."""
    if pid < 0:
        return False
    try:
        os.kill(pid, 0)
    except OSError, e:
        return e.errno == errno.EPERM
    else:
        return True
        
def get_no(pvalue,siglevel,ngenes): # No Correction
	lista=[]
	pp=siglevel
	y=0
	for i in pvalue:
		p=i[0]
		if p<=siglevel:
			lista.append(i)
			y+=1
	return lista,y,pp

def get_b(pvalue,siglevel,ngenes): # Bonferroni
	pvalue.sort()
	lista=[]
	y=0
	#bcorr=siglevel/ngenes
	pp=1.0
	for i in pvalue:
		p=i[0]*ngenes
		if p<=siglevel:
			lista.append(i)
			#lista[i[1]]=i[0]
			y+=1
			if p<pp: pp=p
	#print "Passed:",y,pp
	return lista,y,pp

def get_bh(pvalue,siglevel,ngenes): # B-H
	pvalue.sort()
	#print ngenes
	lista=[]
	x=1
	y=0
	p=0
	for i in pvalue:
		nf=i[0]*ngenes
		fdr=nf/x
		if fdr<=siglevel:
			#dic[i[1]]=i[0]
			lista.append(i)
			p=i[0]
			y+=1
		x+=1
	#print "Passed:",y,p
	return lista,y,p

def getTail(pp):
	if ftail=='l': return pp.left_tail
	elif ftail=='r': return pp.right_tail
	elif ftail=='t': return pp.two_tail
	
def getDicSS(dicp): # dicp = dizionario con le frequenze di sostituzione
	dicpp={}
	for i in dicp:
		if i[0]!=i[1]:
			dicpp[i]=1-dicp[i]
	return dicpp

def getFreads(bases):
	fread={'A':0,'C':0,'G':0,'T':0}
	for i in range(4):
		if i==0: fread['A']=bases[i]
		elif i==1: fread['C']=bases[i]
		elif i==2: fread['G']=bases[i]
		elif i==3: fread['T']=bases[i]
	return fread

def getSub(ref,fread,dics):
	#fread={A,C,G,T}
	nref=fread[ref.upper()]
	sub=[(ref.upper()+i,nref,fread[i]) for i in fread if i!=ref.upper() and fread[i]!=0]
	allsub=' '.join([x[0] for x in sub])
	# lista del tipo [('AT', 50, 10), ('AG', 50, 2)]
	res=[] #[(int(dics[i[0]]*(i[1]+i[2])),((i[1]+i[2])-exp1),pvalue(i[1],i[2],int(dics[i[0]]*(i[1]+i[2])),((i[1]+i[2])-exp1))) for i in sub]
	for i in sub:
		#if binomial:
		#	pval=bdtrc(i[2],i[1]+i[2],(1.0-dics[i[0]]))
		#	#pval=Bprob(i[2],i[1]+i[2],(1.0-dics[i[0]]))
		#	#print i[2],i[1]+i[2],(1.0-dics[i[0]]),pval
		#	obs1,obs2,exp1,exp2=0,0,0,0
		obs1=i[1]
		obs2=i[2]
		exp1=int(dics[i[0]]*(i[1]+i[2]))
		exp2=((i[1]+i[2]) - exp1)
		pval=pvalue(obs1,obs2,exp1,exp2)
		pval=getTail(pval)
		res.append((i[0],obs1,obs2,exp1,exp2,str(pval)))
	if len(res)==1: return res[0][5] #,allsub,fread
	elif len(res) > 1:
		rr=[float(x[-1]) for x in res]
		idx=rr.index(min(rr))
		return res[idx][5] #,allsub,fread
	else: return '1.0' #,0,0
	
def BaseCount(seq,ref):
	b={'A':0,'C':0,'G':0,'T':0}
	subs=[]
	subv=[]
	for i in seq.upper():
		if b.has_key(i): b[i]+=1
	for i in b:
		if not b.has_key(ref): continue
		if b[i]!=0 and i!=ref:
			vv=float(b[i])/(b[i]+b[ref])
			subv.append((b[i],vv,ref+i))
	subv.sort()
	subv.reverse()
	for i in subv:
		if i[0]>=vnuc and i[1]>=mmf: subs.append(i[2])
	freq=0.0
	if len(subs)==0: subs.append('-')
	else: freq=subv[0][1]	
	return sum(b.values()),[b['A'],b['C'],b['G'],b['T']],' '.join(subs),'%.2f'%(freq)

def meanq(v,n):
	try:m=float(v)/n
	except: m=0.0
	return '%.2f'%(m)
	
def rmHomo(sequp,seqdw,gh,ref):
	if len(sequp)==0 and len(seqdw)==0: return 0
	up,dw=0,0
	for i in seqdw:
		if i==ref:dw+=1
		else:break
	for i in sequp[::-1]:
		if i==ref:up+=1
		else:break
	if up+dw+1 >= gh : return 1
	return 0

def prop(tot,va):
	try: av=float(va)/tot
	except: av=0.0
	return av

def vstand(strand):
	vv=[(strand.count('+'),'+'),(strand.count('-'),'-'),(strand.count('*'),'*')]
	if vv[0][0]==0 and vv[1][0]==0: return '*'
	if useconf:
		totvv=sum([x[0] for x in vv[:2]])
		if prop(totvv,vv[0][0])>=strconf: return '+'
		if prop(totvv,vv[1][0])>=strconf: return '-'
		return '*'
	else:
		if vv[0][0]==vv[1][0] and vv[2][0]==0: return '+'
		return max(vv)[1]

def comp(s):
	a={'A':'T','T':'A','C':'G','G':'C'}
	ss=''
	for i in s.upper():
		if a.has_key(i): ss+=a[i]
		else: ss+='N'
	return ss	

def whereis(program):
	for path in os.environ.get('PATH', '').split(':'):
		if os.path.exists(os.path.join(path, program)) and not os.path.isdir(os.path.join(path, program)): return 1
	return 0

def vstrand(lista):
	if len(lista)==0: return '2'
	p=lista.count('+')
	m=lista.count('-')
	if p==len(lista): return '1'
	elif m==len(lista): return '0'
	else: return '2'

def getd(lines): #fixed error in reading strand 6/3/2014
	d={}
	for i in lines:
		l=(i.strip('\n\r')).split('\t')
		if len(l)>=3:
			if l[2]=='+': strand='1'
			elif l[2]=='-': strand='0'
			elif l[2] in '012': strand=l[2]
			else: strand='2'
		else: strand='2'
		d[int(l[1])]=strand
	return d

def normByStrand(seq_,strand_,squal_,mystrand_):
	st='+'
	if mystrand_=='0': st='-'
	seq,qual,squal='',0,[]
	for i in range(len(seq_)):
		if strand_[i]==st:
			seq+=seq_[i]
			if type(i)==int: qual+=ord(squal_[str(i)])-QVAL
			else: qual+=squal_[i]-QVAL
			squal.append(squal_[i])
	return seq,qual,squal

def normByBlat(seq_,strand_,squal_,blatc_):
	seq,qual,squal,strand='',0,[],''
	for i in range(len(seq_)):
		if blatc_[i]=='1':
			seq+=seq_[i]
			if type(i)==int: qual+=ord(squal_[str(i)])-QVAL
			else: qual+=squal_[i]-QVALù
			squal.append(squal_[i])
			strand+=strand_[i]
	return seq,qual,squal,strand

def testBlat(blc):
	if blc.count('1') > blc.count('0'): return 1
	return 0
		
#######
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> START: %s\n"%(script_time))
sys.stderr.write("Analysis ID: %s\n"%(pid))

if not os.path.exists(bamfile):
	usage()
	sys.exit('BAM file %s not found.' %(bamfile))
if sortbam:
	sys.stderr.write('Sorting BAM file.\n')
	pysam.sort(bamfile,'sorted_%s'%(pid))
	os.rename(bamfile,bamfile+'_old')
	os.rename('sorted_%s.bam'%(pid),bamfile)
	sys.stderr.write('Indexing BAM file.\n')
	pysam.index(bamfile)
if not os.path.exists(bamfile+'.bai') and not sortbam:
	sys.stderr.write('Indexing BAM file.\n')
	pysam.index(bamfile)	
if not os.path.exists(fastafile):
	usage()
	sys.exit('Fasta file %s not found.' %(fastafile))
if not os.path.exists(fastafile+'.fai'):
	sys.stderr.write('Indexing Fasta file.\n')
	pysam.faidx(fastafile)
if not os.path.exists(kfile): sys.exit('File containing RNA editing positions not found.')
if sortann:
		if not whereis('grep'): sys.exit('grep command not found.')
		if not whereis('sort'): sys.exit('sort command not found.')
		sys.stderr.write('Sorting file with known editing positions.\n')
		scmd='grep -v ^"chrom" %s | grep -v "^[[:space:]]*$" | sort -k1,1 -k2,2n > %s' %(kfile,'positions_%s'%(pid))
		os.system(scmd)
		os.rename(kfile,kfile+'_old')
		os.rename('positions_%s'%(pid),kfile)
if not os.path.exists(kfile+'.tbi'):
		sys.stderr.write('Indexing file with known positions.\n')
		kfile=pysam.tabix_index(kfile, seq_col=0, start_col=1, end_col=1)
		# Format for tabfile with positions:
		# chr start strand
##################################
# check reference names
rrefs={}
ridxinfo=pysam.idxstats(bamfile)
for j in ridxinfo.split('\n'): #MOD
	l=(j.strip()).split('\t')
	if l[0] in ['*','']: continue #MOD
	if int(l[2])+int(l[3]) > 0: rrefs[l[0]]=int(l[1])
frefs=[]
fidxinfo=open(fastafile+'.fai')
for j in fidxinfo:
	l=(j.strip()).split('\t')
	if l[0]=='': continue
	frefs.append(l[0])
fidxinfo.close()
# in rna-seq
rnof=[]
for i in rrefs.keys():
	if i not in frefs: sys.stderr.write('WARNING: Region %s in RNA-Seq not found in reference file.\n' %(i))
##################################

if uann:
	getstrand=0
	if not os.path.exists(annfile):
		usage()
		sys.exit('Annotation file %s not found.' %(annfile))
	if sortann:
		if not whereis('grep'): sys.exit('grep command not found.')
		if not whereis('sort'): sys.exit('sort command not found.')
		sys.stderr.write('Sorting annotation file.\n')
		scmd='grep ^"#" %s; grep -v ^"#" %s | sort -k1,1 -k4,4n > %s' %(annfile,annfile,'annotation_%s'%(pid))
		os.system(scmd)
		os.rename(annfile,annfile+'_old')
		os.rename('annotation_%s'%(pid),annfile)
	if not os.path.exists(annfile+'.tbi'):
		sys.stderr.write('Indexing annotation file.\n')
		annfile=pysam.tabix_index(annfile, preset='gff')
###########################################################
# Annotation file to exclude positions
if expos:
	if not os.path.exists(exfile):
		usage()
		sys.exit('File %s not found.' %(exfile))
	if sortann:
		if not whereis('grep'): sys.exit('grep command not found.')
		if not whereis('sort'): sys.exit('sort command not found.')
		sys.stderr.write('Sorting file.\n')
		scmd='grep ^"#" %s; grep -v ^"#" %s | sort -k1,1 -k4,4n > %s' %(exfile,exfile,'exfile_%s'%(pid))
		os.system(scmd)
		os.rename(exfile,exfile+'_old')
		os.rename('exfile_%s'%(pid),exfile)
	if not os.path.exists(exfile+'.tbi'):
		sys.stderr.write('Indexing %s file.\n' %(exfile))
		exfile=pysam.tabix_index(exfile, preset='gff')	
###########################################################
#mainbam=pysam.Samfile(bamfile,"rb")
#regions=mainbam.references
#regionslens=mainbam.lengths
#mainbam.close()
#dicregions=dict([(regions[x],regionslens[x]) for x in range(len(regions))])
#chrs=[x for x in regions if x not in nochrs]
dicregions=dict(rrefs.items())
chrs=[x for x in dicregions.keys() if x not in nochrs]
sys.stderr.write('Analysis on %i regions.\n' %(len(chrs)))

if infolder!='': outfolder=os.path.join(outfolder_,'known_%s_%s' %(infolder,pid))
else: outfolder=os.path.join(outfolder_,'known_%s' %(pid))
if not os.path.exists(outfolder):
	splitfolder=os.path.split(outfolder)
	if not os.path.exists(splitfolder[0]): os.mkdir(splitfolder[0])
	os.mkdir(outfolder)	
outtable=os.path.join(outfolder,'outTable_%s' %(pid))
#write command line and input parameters
f=open(os.path.join(outfolder,'parameters.txt'),'w')
f.writelines(params)
f.close()

def exploreBAM(myinput):
	inputs=myinput.split('$')
	chr,bamfile=inputs[0],inputs[1]
	outfile=os.path.join(outfolder,'table_%s_%s'%(chr,pid))
	#outfile2=os.path.join(outfolder,'subs_%s_%s'%(chr,pid))
	d,di={},{}
	bam=pysam.Samfile(bamfile,"rb")
	fasta=pysam.Fastafile(fastafile)
	ktabix=pysam.Tabixfile(kfile)
	lenregion=dicregions[chr]
	if uann: tabix=pysam.Tabixfile(annfile)
	if expos: extabix=pysam.Tabixfile(exfile)
	out=open(outfile,'w')
	#if not custsub:
	#	dsubs=dict([(x+y, 0) for x in 'ACGT' for y in 'ACGT'])
	#	out2=open(outfile2,'w')
	#header='Region\tPosition\tReference\tCoverage\tMeanQuality\tBaseCount\tSubs\tFrequency\n'
	#out.write(header)
	sys.stderr.write('Started analysis on region: %s\n'%(chr))
	if blatr:
		badblat=os.path.join(blatfolder,'blatseqs_%s.bad'%(chr))
		if os.path.exists(badblat):
			sys.stderr.write('Using Blat mapping for region %s\n'%(chr))
			f=open(badblat)
			for i in f:
				l=(i.strip()).split()
				d[l[0]+'_'+l[1]]=int(l[1])
			f.close()
			sys.stderr.write('Found %i reads for region %s\n'%(len(d),chr))
	if exss:
		if os.path.exists(splicefile):
			sys.stderr.write('Loading known splice sites for region %s\n'%(chr))
			f=open(splicefile)
			for i in f:
				l=(i.strip()).split()
				if l[0]!=chr: continue
				st,tp,cc=l[4],l[3],int(l[1])
				if st=='+' and tp=='D':
					for j in range(nss): di[cc+(j+1)]=0
				if st=='+' and tp=='A':
					for j in range(nss): di[cc-(j+1)]=0
				if st=='-' and tp=='D': 	
					for j in range(nss): di[cc-(j+1)]=0
				if st=='-' and tp=='A':
					for j in range(nss): di[cc+(j+1)]=0	
			f.close()
			sys.stderr.write('Loaded %i positions for %s\n'%(len(di),chr))
	if chr in ktabix.contigs:
		for kpos in range(0,lenregion,chunckval):
			startk,endk=kpos,(kpos+chunckval)-1
			kres=[kk for kk in ktabix.fetch(reference=chr,start=startk,end=endk)]
			if len(kres)==0: continue
			kdic=getd(kres)
			#print kdic
			# else explore bam to find exact positions
			for pileupcolumn in bam.pileup(chr,startk,endk,stepper='nofilter', max_depth=MAX_DEPTH):
				if not startk<=pileupcolumn.reference_pos<=endk: continue
				if not kdic.has_key(pileupcolumn.reference_pos+1): continue
				ref=fasta.fetch(chr,pileupcolumn.reference_pos,pileupcolumn.reference_pos+1).upper()
				seq,qual,strand,squal,blatc='',0,'',[],''
				if rmsh:
					if ((pileupcolumn.reference_pos+1)-homo)-1 < 0: sequp=''
					else: sequp=(fasta.fetch(chr,((pileupcolumn.reference_pos+1)-homo)-1,(pileupcolumn.reference_pos+1)-1)).upper()
					seqdw=(fasta.fetch(chr,pileupcolumn.reference_pos+1,(pileupcolumn.reference_pos+1)+homo)).upper()
				for pileupread in pileupcolumn.pileups: # per ogni base dell'allineamento multiplo
					#s,q,t,qq=pileupread.alignment.seq[pileupread.qpos].upper(),ord(pileupread.alignment.qual[pileupread.qpos])-QVAL,'*',pileupread.alignment.qual[pileupread.qpos]
					# escludi posizioni introniche nei pressi di splice sites
					if exss and di.has_key(pileupcolumn.reference_pos+1): continue
					# multiple hit
					if exh:
						if pileupread.alignment.is_secondary: continue
						if pileupread.alignment.has_tag('NH'):
							if pileupread.alignment.get_tag('NH') > 1: continue
					# duplicates
					if exd and pileupread.alignment.is_duplicate: continue
					# se paired end
					if conc: # se devi usare solo le paired
						# se non sono paired
						if not pileupread.alignment.is_paired: continue
						# se non sono concordanti
						if not pileupread.alignment.is_proper_pair: continue
						# se concordanti ma nello stesso orientamento
						flag=pileupread.alignment.flag
						if pileupread.alignment.is_duplicate: flag=flag-1024
						if pileupread.alignment.is_secondary: flag=flag-256
						if flag in [67,131,115,179]: continue
					# mapping quality
					if mq and pileupread.alignment.mapping_quality < MAPQ: continue
					#se la qualita' >= alla qualita' minima
					s,q,t,qq=pileupread.alignment.query_sequence[pileupread.query_position].upper(),pileupread.alignment.query_qualities[pileupread.query_position],'*',pileupread.alignment.query_qualities[pileupread.query_position]
					if q >= MQUAL and pileupcolumn.reference_pos in pileupread.alignment.get_reference_positions():
						#tags=dict(pileupread.alignment.tags)
						#deduci la strand per ogni posizione
						if getstrand:
							#usa le info del mapping se strand oriented 
							if pileupread.alignment.is_read1:
								if unchange1:
									if pileupread.alignment.is_reverse: t='-'
									else: t='+'
								else:
									if pileupread.alignment.is_reverse: t='+'
									else: t='-'
							elif pileupread.alignment.is_read2:
								if unchange2:
									if pileupread.alignment.is_reverse: t='-'
									else: t='+'
								else:
									if pileupread.alignment.is_reverse: t='+'
									else: t='-'
							else: # for single ends
								if unchange1:
									if pileupread.alignment.is_reverse: t='-'
									else: t='+'
								else:
									if pileupread.alignment.is_reverse: t='+'
									else: t='-'	
						if rmnuc:
							#rlen=pileupread.alignment.rlen #pileupread.alignment.qlen #lunghezza della specifica read
							#print rlen,pileupread.qpos,pileupread.alignment.qstart,pileupread.alignment.qend
							# verifica se il nuc deve essere rimosso alle estremita' nel range x-y
							# testare il forward
							#qp=pileupread.qpos #pileupread.qpos-pileupread.alignment.qstart
							#print pileupread.qpos,pileupread.alignment.rlen,len(pileupread.alignment.seq)
							#if pileupread.alignment.is_reverse:
							#	if (rlen-qp)-1 < rmp[0]:continue
							#	if (rlen-qp)-1 > ((rlen)-rmp[1])-1: continue
							#else:
							#	if qp<rmp[0]:continue
							#	if qp>(rlen-rmp[1])-1: continue
							rlen=pileupread.alignment.query_length #pileupread.alignment.qlen #lunghezza della specifica read
							qp=pileupread.query_position #pileupread.qpos-pileupread.alignment.qstart
							if pileupread.alignment.is_reverse:
								if qp>(rlen-rmp[0])-1: continue
								if qp<rmp[1]:continue
							else:
								if qp<rmp[0]:continue
								if qp>(rlen-rmp[1])-1: continue				
						# se la read di appartenenza non mappa in modo univoco con Blat
						if blatr:
							rt=0
							if pileupread.alignment.is_read1: rt=1
							elif pileupread.alignment.is_read2: rt=2
							rname=pileupread.alignment.query_name+'_%i'%(rt)
							if d.has_key(rname): blatc+='0' #continue
							else: blatc+='1'
						# se la base e' diversa dal reference
						# se in regione omopolimerica scarta
						if rmsh and rmHomo(sequp,seqdw,homo,ref): continue
						seq+=s
						qual+=q
						strand+=t
						squal.append(qq)
				if seq.strip()!='':
					if blatr:
						if testBlat(blatc): seq,qual,squal,strand=normByBlat(seq,strand,squal,blatc)
						else: continue
					#print pileupcolumn.reference_pos+1,seq,squal
					#mystrand=kdic[pileupcolumn.reference_pos+1]
					#print mystrand
					try: mystrand=kdic[pileupcolumn.reference_pos+1]
					except: mystrand='2'
					#print chr,pileupcolumn.reference_pos+1,seq,strand, mystrand
					if uann and not getstrand:
						if chr in tabix.contigs:
							sres=[kk.strand for kk in tabix.fetch(reference=chr,start=(pileupcolumn.reference_pos),end=(pileupcolumn.reference_pos+1),parser=pysam.asGTF())]
							mystrand=vstrand(sres)
					if getstrand and not uann:
						mystr=vstand(strand)
						if mystr=='-': mystrand='0'
						elif mystr=='+': mystrand='1'
						else: mystrand='2'
					if mystrand=='0':
						seq=comp(seq)
						ref=comp(ref)
					#if getstrand and mystrand in ['1','0'] and not useconf: seq,qual,squal=normByStrand(seq,strand,squal,mystrand)
					if getstrand and mystrand in ['1','0'] and corrstr: seq,qual,squal=normByStrand(seq,strand,squal,mystrand)
					if uann and mystrand in ['1','0'] and corrstr: seq,qual,squal=normByStrand(seq,strand,squal,mystrand)
					#if not getstrand and not uann and mystrand in ['1','0']: seq,qual,squal=normByStrand(seq,strand,squal,mystrand)
					#print chr,pileupcolumn.reference_pos+1,seq,strand,mystrand
					cov,bcomp,subs,freq=BaseCount(seq,ref)
					if cov < MINCOV: continue
					if exms and subs.count(' ')>0: continue
					mqua=meanq(qual,len(seq))
					if expos:
						if chr in extabix.contigs:
							exres=[kk for kk in extabix.fetch(reference=chr,start=(pileupcolumn.reference_pos),end=(pileupcolumn.reference_pos+1))]
							if len(exres)>0: continue			
					line='\t'.join([chr,str(pileupcolumn.reference_pos+1),ref,mystrand,str(cov),(mqua),str(bcomp),subs,freq])+'\n'
					out.write(line)
	bam.close()
	fasta.close()
	ktabix.close()
	out.close()
	if uann: tabix.close()
	if expos: extabix.close()
	sys.stderr.write('Job completed for region: %s\n'%(chr))

def addPvalue(myinput2): # not used here 
	inputs=myinput2.split('$')
	f=open(inputs[0])
	subs=eval((f.readline()).strip())
	f.close()
	dsubs={}
	for i in subs: dsubs[i]=float(subs[i])/sum(subs.values())
	dsubss=getDicSS(dsubs)
	#print dsubss
	o=open(inputs[2],'w')
	f=open(inputs[1])
	for i in f:
		l=(i.strip()).split('\t')
		if i.strip()=='': continue
		#if i.startswith('Region'):
		#	l.append('Pvalue')
		#	o.write('\t'.join(l)+'\n')
		#	continue
		if l[6]!='-': pval=getSub(l[2],getFreads(eval(l[6])),dsubss)
		else: pval='1.0'
		l.append(pval)
		o.write('\t'.join(l)+'\n')
	o.close()
	
def do_work(q):
	while True:
		try:
			x=q.get(block=False)
			exploreBAM(x)
		except Empty:
			break

work_queue = Queue()
for i in chrs:
	strinput=i+'$'+bamfile
	work_queue.put(strinput)
processes=[Process(target=do_work, args=(work_queue,)) for i in range(NCPU)]
for t in processes:
	t.start()
for t in processes:
	t.join()
time.sleep(0.5)
#
head='Region\tPosition\tReference\tStrand\tCoverage-q%i\tMeanQ\tBaseCount[A,C,G,T]\tAllSubs\tFrequency\n' %(MQUAL)
sys.stderr.write('Merging Tables.\n')
o=open(outtable,'w')
if noheader==0: o.write(head)
for i in chrs:
	#tabfile=os.path.join(outfolder,'outTable_%s_%s' %(i,pid))
	tabfile=os.path.join(outfolder,'table_%s_%s' %(i,pid))
	if os.path.exists(tabfile):
		f=open(tabfile)
		for j in f: o.write(j)
		f.close()
		os.remove(tabfile)
		#os.remove(intabfile)
o.close()

#if sigsites:
#	sys.stderr.write('Selecting significant sites.\n')
#	outsig=os.path.join(outfolder,'outTableSig_%s' %(pid))
#	f=open(outtable)
#	o=open(outsig,'w')
#	o.write(head)
#	allv=[]
#	for i in f:
#		if i.startswith('Region'): continue
#		if i.strip()=='': continue
#		l=(i.strip()).split('\t')
#		if l[7]=='-': continue
#		if l[7] not in usubs: continue
#		pp=float(l[9])
#		allv.append((pp,i))
#	if test=='bh': rr=get_bh(allv,sval,len(allv))
#	elif test=='bo': rr=get_b(allv,sval,len(allv))
#	else: rr=get_no(allv,sval,len(allv))
#	for i in rr[0]: o.write(i[1])		
#	f.close()
#	o.close()

sys.stderr.write('Results saved on %s\n'%(outtable))
#if sigsites: sys.stderr.write('Significant sites saved on %s\n'%(outsig))

script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> END: %s\n"%(script_time))

