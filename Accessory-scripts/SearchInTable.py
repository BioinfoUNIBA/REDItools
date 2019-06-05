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

import sys, os, getopt, time
try: import pysam
except: sys.exit('Pysam module not found.')
#pid=str(os.getpid()+random.randint(0,999999999))

def usage():
	print """
USAGE: python SearchInTable.py [options]
Options:
-i		Sorted table file (first col=reference; second col=coordinate 1 based)
		or tabix indexed table (ending with .gz)
-q		Query (file or single positions: chr21:123456)
-C		Sequence name column [1]
-S		Start column [2]
-E		End column; can be identical to '-S' [2]
-P		Print to stdout found lines
-p		Print position header (like a fasta header >chr21:123456)
-n		Print "Not found"
-s		Print simple statistics on standard error
-k		Skip lines starting with in query file
-o		Save found/not found positions on file
-h		Print this help

"""
#-k		skip first INT lines [0]

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:q:k:pso:hnC:S:E:O:P",["help"])
	if len(opts)==0:
		usage()
		sys.exit(2)
except getopt.GetoptError as err:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

tablefile=''
query=''
outfile=''
outfile2=''
pr,prn,prf=0,0,0
ps=0
sv,sv2=0,0
sk=0
ski=''
skil=0

scol,stcol,ecol=0,1,1
for o, a in opts:
	if o in ("-h","--help"):
		usage()
		sys.exit()
	elif o == "-i":
		tablefile=a
		if not os.path.exists(tablefile):
			usage()
			sys.exit('Table file not found')
	elif o == "-q":
		query=a
		if query=='':
			usage()
			sys.exit('Query not found.')
	elif o == "-p": pr=1
	elif o == "-C": scol=int(a)-1
	elif o == "-S": stcol=int(a)-1
	elif o == "-E": ecol=int(a)-1	
	elif o == "-n": prn=1
	elif o == "-P": prf=1
	elif o == "-k":
		ski=a
		skil=1
	elif o == "-s": ps=1
	elif o == "-o":
		outfile=a
		sv=1
	elif o == "-O":
		outfile2=a
		sv2=1
	else:
		assert False, "Unhandled Option"


script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stdout.write("Script time --> START: %s\n"%(script_time))
if not os.path.exists(tablefile):
	#sys.stderr.write('Compressing table file.\n')
	#pysam.tabix_index(tablefile, tablefile+'.gz')
	sys.stderr.write('Indexing table file.\n')
	tablefile=pysam.tabix_index(tablefile, seq_col=scol, start_col=stcol, end_col=ecol)
#if tablefile.endswith('.gz') and not tablefile.endswith('.tbi'):
#	tablefile=pysam.tabix_index(tablefile, seq_col=scol, start_col=stcol, end_col=ecol)

tabix=pysam.Tabixfile(tablefile)
allref=tabix.contigs
positions=[]
if os.path.exists(query):
	f=open(query)
	for i in f:
		if i.strip()=='': continue
		if i.startswith('#'): continue
		if i.startswith('Region'): continue
		if skil:
			if i.startswith(ski): continue
		l=(i.strip()).split()
		positions.append((l[0],int(l[1])-1))
	f.close()
elif query.count(":")==1:
	l=(query.strip()).split(':')
	positions.append((l[0],int(l[1])-1))
else: sys.exit('I cannot read the query.')

if sv:
	outf=open(outfile+'_found','w')
	outnf=open(outfile+'_notfound','w')
if sv2:
	outf2=open(outfile2+'_foundInSortedTable','w')
xx=0
for pos in positions:
	res=[]
	if pos[0] in allref:
		res=[kk for kk in tabix.fetch(reference=pos[0],start=pos[1],end=pos[1]+1)]
	if pr: sys.stdout.write('>%s:%i\n' %(pos[0],pos[1]+1))
	if len(res)==0:
		if prn: sys.stdout.write('Not Found\n')
		if sv: outnf.write('%s\t%i\n' %(pos[0],pos[1]+1))
	else:
		#if sv: outf.write(res[0]+'\n')
		if sv: outf.write(res[0]+'\n')
		if prf: sys.stdout.write(res[0]+'\n')
		xx+=1	
tabix.close()
if sv:
	outf.close()
	outnf.close()
if ps:
	sys.stdout.write('Positions in query: %i\n' %(len(positions)))
	sys.stdout.write('Positions found: %i\n' %(xx))
	sys.stdout.write('Positions not found: %i\n' %(len(positions)-xx))
if sv:
	sys.stdout.write('Found line(s) saved on: %s\n' %(outfile+'_found'))
	sys.stdout.write('Not found line(s) saved on: %s\n' %(outfile+'_notfound'))
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stdout.write("Script time --> END: %s\n"%(script_time))