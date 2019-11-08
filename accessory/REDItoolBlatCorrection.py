#!/usr/bin/env/python
# UPDATE 07-11-2019
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

import sys, os, time, math, random, getopt
try: import pysam
except: sys.exit('pysam module not found.\nPlease install it before.')
from multiprocessing import Process, Queue
from Queue import Empty

pid=str(os.getpid()+random.randint(0,999999999))

def usage():
	print """
USAGE: python REDItoolBlatCorrection.py [options]
Options:
-i		BAM file
-I		Sort input BAM file
-f		Genomic fasta file
-F		Genomic fasta file in 2bit format for gfServer
-t		Num. of working threads [1]
-o		Output Folder [BlatCorrection_%s]
-k		List of chromosomes to skip separated by comma or file
-r		Regions in GFF
-s		Sort GFF
-q		Min. quality score [25]
-Q		Fastq offset value [33]
-V		Verify if gfServer is alive
-T		Stop gfServer at script end
-h		Print this help

"""%(pid)

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:f:t:o:k:q:Q:hIF:VTr:s")
except getopt.GetoptError as err:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

sortann=0
regfile='' # file with regions in GFF
usereg=0
bamfile=''
fastafile=''
fastafile2bit=''
nochrs=[] #(sys.argv[4]).split(',')
NCPU=1
BLATdb='/'
QVAL=33
MQUAL=25
outbase='BlatCorrection_%s' %(pid)
sortbam=0
gfa,isalive,sgf=0,0,0
for o, a in opts:
	if o == "-h":
		usage()
		sys.exit()
	elif o == "-i": bamfile=a
	elif o == "-I": sortbam=1
	elif o == "-f": fastafile=a
	elif o == "-r": 
		regfile=a
		usereg=1
	elif o == "-F": fastafile2bit=a
	elif o == "-t": NCPU=int(a)
	elif o == "-o": outbase=a
	elif o == "-V": gfa=1
	elif o == "-T": sgf=1
	elif o == "-s": sortann=1
	elif o == "-k":
		if os.path.exists(a):
			f=open(a)
			nochrs=[x.strip() for x in f if x.strip()!='']
			f.close()
		else: nochrs=[x for x in a.split(',') if x.strip()!='']
	elif o == "-q": MQUAL=int(a)
	elif o == "-Q": QVAL=int(a)
	else:
		assert False, "Unhandled Option"

def whereis(program):
	for path in os.environ.get('PATH', '').split(':'):
		if os.path.exists(os.path.join(path, program)) and not os.path.isdir(os.path.join(path, program)): return 1
	return 0
#######
def BaseCount(seq,ref):
	b={'A':0,'C':0,'G':0,'T':0}
	subs=[]
	for i in seq.upper():
		if b.has_key(i): b[i]+=1
	for i in b:
		if b[i]!=0 and i!=ref: subs.append(ref+i)
	if len(subs)==0: subs.append('-')
	return str([b['A'],b['C'],b['G'],b['T']]),' '.join(subs)

#### for blat
def getPS(line):
	pid = (100.0 - (pslCalcMilliBad(line) * 0.1))
	score = pslScore(line)
	#print "The percentage:",pid
	#print "Score:",score
	return pid,score

def pslScore(cols):
	sizeMul =  1
	return sizeMul * (int(cols[0]) + (int(cols[2]))) - sizeMul * int(cols[1]) - int(cols[4]) - int(cols[6])

def round(number):
	return int(number + .5);

def pslCalcMilliBad(cols):
	sizeMul = 1
	# cols[0]  matches
	# cols[1]  misMatches
	# cols[2]  repMaches
	# cols[4]  qNumInsert
	# cols[6]  tNumInsert
	# cols[11] qStart
	# cols[12] qEnd
	# cols[15] tStart
	# cols[16] tEnd
	qAliSize = sizeMul * (int(cols[12]) - int(cols[11]))
	tAliSize = int(cols[16]) - int(cols[15])
	# I want the minimum of qAliSize and tAliSize
	if qAliSize < tAliSize: aliSize = qAliSize #? $aliSize = $qAliSize : $aliSize =  $tAliSize;
	else: aliSize = tAliSize
	# return 0 is AliSize == 0
	if aliSize <= 0: return 0
	# size diff
	sizeDiff = qAliSize - tAliSize
	if sizeDiff < 0: sizeDiff = 0  
	# insert Factor
	insertFactor = int(cols[4]) 
	# $insertFactor += $cols[6];
	milliBad = (1000 * (int(cols[1])*sizeMul + insertFactor + round(3*math.log( 1 + sizeDiff)))) / (sizeMul * (int(cols[0]) + int(cols[2])+ int(cols[1])))
	return milliBad

def com(num,list):
	for i in list:
		if i[0]<=num<=i[1]: return 1
	return 0

def min95(val,score):
	#if val < (score*95.0)/100: return 1
	if val < (score*0.95): return 1
	return 0

def readLines(lines):
	res=[]
	for line in lines:
		pidd,score=getPS(line)
		#print pidd,score
		sp=[int(x) for x in (line[18].strip(',')).split(',')]
		tstarts=[int(x) for x in (line[20].strip(',')).split(',')]
		ex=[(tstarts[x]+1,tstarts[x]+sp[x]) for x in range(len(sp))]
		nl=[line[9],score,str(int(line[11])+1),line[12],str(line[10]),pidd,line[13],line[8],int(line[15])+1,int(line[16]),ex,int(line[0])]
		res.append((int(line[0]),nl)) #score
		#if d.has_key(line[9]): d[line[9]].append((score,nl))
		#else: d[line[9]]=[(score,nl)]
	return res

def comp(ri,hits):
	g,ng=0,0
	hits.sort()
	hits.reverse()
	if len(hits)==1: #unique hit with editing candidate position included
		if hits[0][1][6]==ri[2] and com(ri[1],hits[0][1][10]): g+=1 #float(hits[0][1][5])>=90.0
		else: ng+=1
	elif len(hits)>1: #multiple hits
		if hits[0][1][6]==ri[2] and min95(hits[1][0],hits[0][0]): # if second best score less than 95% of first best score
			if com(ri[1],hits[0][1][10]): g+=1 # if first best hit include editing position
			else: ng+=1
		else: ng+=1
	if g>ng: return 1
	return 0

def readPSL(infile,dread,outfile):
	f=open(infile)
	o=open(outfile,'w')
	name,lines,xx='',[],0
	while 1:
		line=f.readline()
		if not line:
			if name=='': break
			oread=dread[name]
			bread=readLines(lines)
			badr=''
			if len(bread)==0: badr=name
			else: 
				if not comp(oread,bread): badr=name
			if badr!='':
				o.write(name[:-2]+' '+name[-1]+'\n')
				xx+=1
			break
		if line.strip()=='': continue	
		if line.startswith('psL'): continue
		if (line.strip()).startswith('match'): continue
		if line.startswith('-'): continue
		l=(line.strip()).split('\t')
		if l[9]!=name:
			if len(lines)!=0:
				oread=dread[name]
				bread=readLines(lines)
				badr=''
				if len(bread)==0: badr=name
				else: 
					if not comp(oread,bread): badr=name
				if badr!='':
					o.write(name[:-2]+' '+name[-1]+'\n')
					xx+=1
			lines=[l]
			name=l[9]
		else: lines.append(l)
	f.close()
	o.close()
	return xx

def readgf(infile):
	f=open(infile)
	for i in f:
		if 'Server ready for queries!' in i:
			f.close()
			return 1
	f.close()
	return 0

def parse(line):
	l=(line.strip()).split('\t')
	cc=(int(l[3]),int(l[4]))
	return cc
	
####

#######
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> START: %s\n"%(script_time))
if not whereis('samtools'): sys.exit('Samtools not found in your path.')
if not whereis('gfClient'): sys.exit('gfClient not found in your path.')
if not whereis('gfServer'): sys.exit('gfServer not found in your path.')
if not os.path.exists(bamfile):
	usage()
	sys.exit('BAM file %s not found.' %(bamfile))
if sortbam:
	sys.stderr.write('Sorting input BAM file.\n')
	pysam.sort(bamfile,'sorted_%s'%(pid))
	os.rename(bamfile,bamfile+'_old')
	os.rename('sorted_%s.bam'%(pid),bamfile)
	sys.stderr.write('Indexing input BAM file.\n')
	pysam.index(bamfile)
if not os.path.exists(bamfile+'.bai') and not sortbam:
	sys.stderr.write('Indexing input BAM file.\n')
	pysam.index(bamfile)

if not os.path.exists(fastafile):
	usage()
	sys.exit('Fasta file %s not found.' %(fastafile))
if not os.path.exists(fastafile+'.fai'):
	sys.stderr.write('Indexing Fasta file.\n')
	pysam.faidx(fastafile)

if not os.path.exists(fastafile2bit):
	usage()
	sys.exit('Fasta file %s in 2bit format not found.' %(fastafile2bit))
	
if usereg:
	if not os.path.exists(regfile):
		usage()
		sys.exit('GFF file %s not found.' %(regfile))
	if sortann:
		if not whereis('grep'): sys.exit('grep command not found.')
		if not whereis('sort'): sys.exit('sort command not found.')
		sys.stderr.write('Sorting GFF file.\n')
		scmd='grep ^"#" %s; grep -v ^"#" %s | sort -k1,1 -k4,4n > %s' %(regfile,regfile,'workf_%s'%(pid))
		os.system(scmd)
		os.rename(regfile,regfile+'_old')
		os.rename('workf_%s'%(pid),regfile)
	if not os.path.exists(regfile+'.tbi'):
		sys.stderr.write('Indexing GFF file.\n')
		regfile=pysam.tabix_index(regfile, preset='gff')	

###
#simple test to view if gfServer is alive
#important info for gfServer and gfClient added 12-4-2014
if not os.path.isabs(fastafile2bit):
	fastafile2bit=os.path.abspath(fastafile2bit)
BLATdb=os.path.dirname(fastafile2bit)
BLATdbfile=os.path.basename(fastafile2bit)

if gfa:
	#testseq=''
	#for i in range(50): testseq+=random.choice('ACGT')
	#f=open('a','w')
	#f.write('>testseq\n%s\n'%(testseq))
	#f.close()
	#cmd='gfClient localhost 3333 %s -minScore=20 -minIdentity=0 -nohead a b > c 2> d' %(BLATdb)
	cmd='gfServer status localhost 3333 > c 2> d'
	sys.stderr.write('gfServer test...\n')
	os.system(cmd)	
	f=open('d')
	if not (f.readline()).startswith('Error'):
		sys.stderr.write('gfServer is alive.\n')
		isalive=1
	f.close()
	#os.remove('a')
	#os.remove('b')
	os.remove('c')
	os.remove('d')
###
if not isalive:
	sys.stderr.write('Starting gfServer\n')
	cpath=os.getcwd()
	os.chdir(BLATdb)
	os.system('gfServer -stepSize=5 -repMatch=2253 -canStop -log=%s start localhost 3333 %s > %s &' %(os.path.join(cpath,'gfserver_'+pid),BLATdbfile,os.path.join(cpath,'gfout_'+pid)))
	os.chdir(cpath)
	time.sleep(5) # aspetta 5 sec 
	gfready=0
	while not gfready:
		gfready=readgf('%s' %(os.path.join(cpath,'gfserver_'+pid)))
		time.sleep(5)
	sys.stderr.write('gfServer ready\n')

mainbam=pysam.Samfile(bamfile,"rb")
regions=mainbam.references
contigs=[]
if usereg:
	rfile=pysam.Tabixfile(regfile)
	contigs=rfile.contigs
	rfile.close()
chrs={}
for i in regions:
	if i in nochrs: continue
	if usereg and not i in contigs: continue
	conta=mainbam.count(i)
	inblat=os.path.join(outbase,'blatseqs_%s.fa'%(i))
	outblat=os.path.join(outbase,'blatseqs_%s.psl'%(i))
	infoblat=os.path.join(outbase,'blatseqs_%s.info'%(i))
	badblat=os.path.join(outbase,'blatseqs_%s.bad'%(i))
	chrbam=bamfile #os.path.join(outbase,'%s.bam'%(i))
	if conta > 0:
		#if chrs.has_key(i): genes[i].append((conta,bname))
		chrs[i]=[conta,inblat,outblat,badblat,infoblat,chrbam,regfile]
nreads=sum([x[0] for x in chrs.values()])
sys.stderr.write('Found %i reads in %i chromosomes.\n' %(nreads,len(chrs)))
sys.stderr.write('Excluded %i chromosomes.\n' %(len(nochrs)))
mainbam.close()
if not os.path.exists(outbase): os.mkdir(outbase)
# main function

def createFasta(chr,inblat,infoblat,chrbam,fastafile,rgfile):
	d={}
	fasta=pysam.Fastafile(fastafile)
	bam=pysam.Samfile(chrbam,"rb")
	if usereg: tabix=pysam.Tabixfile(rgfile)
	sys.stderr.write('Started analysis on region: %s\n'%(chr))
	blatfasta=open(inblat,'w')
	inb=open(infoblat,'w')
	if usereg:	
		for reg in tabix.fetch(reference=chr):
			cc=parse(reg)
			xx=1
			for pileupcolumn in bam.pileup(chr,cc[0]-1,cc[1]):
				if not cc[0]<=pileupcolumn.pos+1<=cc[1]: continue
				#print cc
				#print 'Found',pileupcolumn.pos+1
				ref=fasta.fetch(chr,pileupcolumn.pos,pileupcolumn.pos+1).upper()
				#print pileupcolumn.pos+1,ref
				for pileupread in pileupcolumn.pileups:
					if ord(pileupread.alignment.qual[pileupread.query_position_or_next])-QVAL >= MQUAL and pileupcolumn.pos in pileupread.alignment.positions:
						seq=pileupread.alignment.seq[pileupread.query_position_or_next].upper()
						if ref!=seq: #substitution
							rt=0
							if pileupread.alignment.is_read1: rt=1
							elif pileupread.alignment.is_read2: rt=2
							rname=pileupread.alignment.qname+'_%i'%(rt)
							#rname=pileupread.alignment.qname+"$$"+chr+"-"+str(pileupcolumn.pos+1)+'-'+str(xx)+'_%i'%(rt)
							if not d.has_key(rname):
								#iline=(rname,pileupread.alignment.positions[0]+1,chr)
								iline=(rname,pileupcolumn.pos+1,chr)
								inb.write(str(iline)+'\n')
								fline='>'+rname+'\n'+pileupread.alignment.seq+'\n'				
								blatfasta.write(fline)
								xx+=1		
	else:
		for pileupcolumn in bam.pileup(chr):
			ref=fasta.fetch(chr,pileupcolumn.pos,pileupcolumn.pos+1).upper()
			for pileupread in pileupcolumn.pileups:
				if ord(pileupread.alignment.qual[pileupread.query_position_or_next])-QVAL >= MQUAL and pileupcolumn.pos in pileupread.alignment.positions:
					seq=pileupread.alignment.seq[pileupread.query_position_or_next].upper()
					if ref!=seq: #substitution
						rt=0
						if pileupread.alignment.is_read1: rt=1
						elif pileupread.alignment.is_read2: rt=2
						rname=pileupread.alignment.qname+'_%i'%(rt)
						if not d.has_key(rname):
							#iline=(rname,pileupread.alignment.positions[0]+1,chr)
							iline=(rname,pileupcolumn.pos+1,chr)
							inb.write(str(iline)+'\n')
							fline='>'+rname+'\n'+pileupread.alignment.seq+'\n'				
							blatfasta.write(fline)		
	blatfasta.close()
	inb.close()
	bam.close()
	fasta.close()
	if usereg: tabix.close()

def runSamtools(chrbam,rgfile):
	sys.stderr.write('Creating %s\n' %(chrbam))
	cmd1='samtools view -b %s %s > %s' %(bamfile,os.path.split(chrbam[:-4])[-1],chrbam)
	cmd2='samtools index %s' %(chrbam)
	os.system(cmd1)
	os.system(cmd2)

def createBad(strline):
	info=strline.split('$')
	#sys.stderr.write('Splitting BAM...\n')
	chr,inblat,outblat,badblat,infoblat,chrbam,rgfile=info[0],info[1],info[2],info[3],info[4],info[5],info[6]
	#runSamtools(chrbam,rgfile)
	#cmd='python createFasta.py %s %s %s %s %s %i %i' %(chr,inblat,infoblat,chrbam,fastafile,QVAL,MQUAL)
	#os.system(cmd)
	createFasta(chr,inblat,infoblat,chrbam,fastafile,rgfile)
	#blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -noHead hg19_softmasked_ALLCHR.fa $fafile $pslfile
	cmd='gfClient localhost 3333 %s -minScore=20 -minIdentity=0 -nohead %s %s > /dev/null' %(BLATdb,inblat,outblat)
	sys.stderr.write('Running Blat on %s.\n'%(chr))
	os.system(cmd)
	dread={}
	f=open(infoblat)
	for i in f:
		l=eval(i.strip())
		dread[l[0]]=l
	f.close()
	sys.stderr.write('Parsing Blat on %s.\n'%(chr))
	numbadr=readPSL(outblat,dread,badblat)
	os.remove(outblat)
	os.remove(inblat)
	os.remove(infoblat)
	##os.remove(chrbam)
	##os.remove(chrbam+'.bai')
	sys.stderr.write('Analysis on %s finished. Bad reads: %i.\n'%(chr,numbadr))

def do_work(q):
	while True:
		try:
			x=q.get(block=False)
			createBad(x)
		except Empty:
			break

work_queue = Queue()
for i in chrs:
	files=chrs[i]
	strline=i+'$'+'$'.join(files[1:])
	work_queue.put(strline)
processes=[Process(target=do_work, args=(work_queue,)) for i in range(NCPU)]
for t in processes:
	t.start()
for t in processes:
	t.join()
	time.sleep(2)
time.sleep(2)

if sgf:
	sys.stderr.write('Stopping gfServer\n')
	os.system('gfServer stop localhost 3333 > gfstop_%s &' %(pid))
	time.sleep(2)
	for i in os.listdir('.'):
		if i.startswith('gfstop_'): os.remove('gfstop_%s' %(pid))
		if i.startswith('gfserver_'): os.remove('gfserver_%s' %(pid))
		if i.startswith('gfout_'): os.remove('gfout_%s' %(pid))
else:
	sys.stderr.write('gfServer not stopped!\nUse the command "gfServer stop localhost 3333" to stop the process.\n')
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> END: %s\n"%(script_time))
