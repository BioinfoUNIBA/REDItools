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

import sys, os, getopt, random, time
try: import pysam
except: sys.exit('Pysam module not found.')
pid=str(os.getpid()+random.randint(0,999999999))

pysamVersion=pysam.__version__
sys.stderr.write('Pysam version used: %s\n' %(pysamVersion))

def usage():
	print """
USAGE: python AnnotateTable.py [options]
Options:
-a		Sorted Annotation file
-i		Annotate a file of positions [column1=region, column2=coordinate (1 based)]
		or a single position [region:coordinate (1 based)]
-k		skip lines starting with: #
-r		Add a prefix to chromosome name [] (chr when the name is a number)
-s		Strand column in annotation file [4]
-u		Not use table strand info (fix it to 2)
-c		Add columns separated by comma (feature:1, gene_id:2, transcript_id:3) [1,2]
-n		Column name [Col]
-S		Correct strand by annotation
-C		Columns with base distribution [7,12] (in combination with -S)
-o		Save lines to a file
-h		Print this help
"""

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:a:o:hs:c:n:SC:uk:r:',["help"])
except getopt.GetoptError, err:
	print str(err) 
	usage()
	sys.exit()

if len(opts)==0:
	usage()
	sys.exit()
tablefile,outfile,annfile='','',''
save,ap,af,addc,cs,nos=0,0,0,[0,1],0,0
csc=[6,11]
strcol=3
colname='Col'
skip='Region'
addchr=''
for o,a in opts:
	if o in ("-h","--help"):
		usage()
		sys.exit()
	elif o == "-n": colname = a
	elif o == "-k": skip = a
	elif o == "-r": addchr = a	
	elif o == "-i":
		tablefile = a
		if not os.path.exists(tablefile): ap,af=1,0
		else: ap,af=0,1
	elif o == "-o":
		outfile = a
		save=1
	elif o == "-s": strcol = int(a)-1
	elif o == "-S": cs = 1
	elif o == "-u": nos = 1
	elif o == "-C": csc=[int(x)-1 for x in a.split(',')]
	elif o == "-c":
		addc = [int(x)-1 for x in a.split(',') if x in ['1','2','3']]
		addc.sort()
	elif o == "-a":
		annfile = a
		if annfile=='':
			usage()
			sys.exit('Sorted annotation file not found.')		
	else:
		assert False, "unhandled option"

##############
def gstr(v):
	if v=='-': return '0'
	else: return '1'

def comp(s):
	a={'A':'T','T':'A','C':'G','G':'C'}
	ss=''
	for i in s.upper():
		if a.has_key(i): ss+=a[i]
		elif i==' ': ss+=' '
		elif i=='-': ss+='-'
		else: ss+='N'
	return ss

def bcomp(b):
	bb=eval(b)
	return str([bb[3],bb[2],bb[1],bb[0]])

def checkstr(stringa):
	strand='+-'
	if stringa=='0': strand='-'
	elif stringa=='1': strand='+'
	elif stringa=='2': strand='+-'
	elif stringa=='-': strand='-'
	elif stringa=='+': strand='+'
	return strand

def parse(res):
	d={'+':{},'-':{}}
	anns='+'
	for i in res:
		if i[3]=='+':
			if d['+'].has_key(i[1]):
				if i[0] not in d['+'][i[1]][0]: d['+'][i[1]][0]=d['+'][i[1]][0]+','+i[0]
				if i[2]+'-'+i[0] not in d['+'][i[1]][1]: d['+'][i[1]][1]=d['+'][i[1]][1]+','+i[2]+'-'+i[0]
			else:
				d['+'][i[1]]=[i[0],i[2]+'-'+i[0]]
		elif i[3]=='-':
			if d['-'].has_key(i[1]):
				if i[0] not in d['-'][i[1]][0]: d['-'][i[1]][0]=d['-'][i[1]][0]+','+i[0]
				if i[2]+'-'+i[0] not in d['-'][i[1]][1]: d['-'][i[1]][1]=d['-'][i[1]][1]+','+i[2]+'-'+i[0]
			else:
				d['-'][i[1]]=[i[0],i[2]+'-'+i[0]]
	gip='$'.join(d['+'].keys())
	featp='$'.join([d['+'][x][0] for x in d['+'].keys()])
	tip='$'.join([d['+'][x][1] for x in d['+'].keys()])
	gim='$'.join(d['-'].keys())
	featm='$'.join([d['-'][x][0] for x in d['-'].keys()])
	tim='$'.join([d['-'][x][1] for x in d['-'].keys()])
	p=[featp,gip,tip]
	m=[featm,gim,tim]
	pm=[(featp+'&'+featm).strip('&'),(gip+'&'+gim).strip('&'),(tip+'&'+tim).strip('&')]
	if len(d['+'])==0 and len(d['-'])!=0: anns='-'
	if len(d['+'])==0: p=['-','-','-']
	if len(d['-'])==0: m=['-','-','-']
	if len(d['+'])==0 and len(d['-'])==0:
		pm=['-','-','-']
		anns='+-'
	if len(d['+'])!=0 and len(d['-'])!=0: anns='+-'
	return (p,m,pm,anns)

#chr17:7590770
	
###############
if ap and af:
	usage()
	sys.exit('You can annotate a file of positions or a single positions but not both in one run.')
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> START: %s\n" %(script_time))


if not os.path.exists(annfile+'.tbi'):
	sys.stderr.write('Indexing %s file.\n' %(annfile))
	annfile=pysam.tabix_index(annfile, preset='gff')

tabix=pysam.Tabixfile(annfile)
contig=tabix.contigs

if ap:
	prinfo=['Feature --> ','Gid --> ','Tid --> ']
	try:
		query=tablefile.split(':')
		chr,pos=addchr+query[0],int(query[1])-1
		try: strand=checkstr(query[2])
		except: strand=checkstr('')
		if nos: strand='+-'
		sres=[]
		if chr in contig:
			sres=[(kk.feature,kk.gene_id,kk.transcript_id,kk.strand) for kk in tabix.fetch(reference=chr,start=pos,end=pos+1,parser=pysam.asGTF())]	
		ann=parse(sres)
		if strand=='+': res=ann[0]
		elif strand=='-': res=ann[1]
		else: res=ann[2]		
		for i in addc:
			print prinfo[i]+ res[i]
	except: sys.exit('Error: not correct position.')
	
if af:
	if save: o=open(outfile,'w')
	f=open(tablefile)
	hinfo=['%s_feat' %(colname),'%s_gid' %(colname),'%s_tid' %(colname)]
	for i in f:
		if i.strip()=='': continue
		if i.startswith('Region'):
			h=[i.strip()]
			for k in addc: h.append(hinfo[k])
			if save: o.write('\t'.join(h)+'\n')
			else: print '\t'.join(h)
			continue
		if i.startswith(skip): continue
		l=(i.strip()).split('\t')
		chr,pos=addchr+l[0],int(l[1])-1
		try: strand=checkstr(l[strcol])
		except: strand='+-'
		if nos: strand='+-'
		sres=[]
		#print chr,pos,pos+1
		if chr in contig:
			sres=[(kk.feature,kk.gene_id,kk.transcript_id,kk.strand) for kk in tabix.fetch(reference=chr,start=pos,end=pos+1,parser=pysam.asGTF())]	
		ann=parse(sres) #(p,m,pm,anns)
		if cs:
			if ann[3]=='+-': pass
			elif ann[3]==strand: pass
			elif ann[3]!=strand:
				l[2]=comp(l[2])
				l[strcol]=gstr(ann[3])
				strand=l[strcol]
				for j in csc:
					try:
						l[j]=bcomp(l[j])
						l[j+1]=comp(l[j+1])
					except: pass
		if strand=='+': res=ann[0]
		elif strand=='-': res=ann[1]
		else: res=ann[2]
		for j in addc: l.append(res[j])
		if save: o.write('\t'.join(l)+'\n')
		else: print '\t'.join(l)
tabix.close()
if save:
	o.close()
	sys.stderr.write("Table saved on %s\n" %(outfile))
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> END: %s\n" %(script_time))
