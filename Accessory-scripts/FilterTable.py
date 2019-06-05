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

"""
To do: filtering according to strand of positions in table file
"""

import sys, time, getopt, string, os, random
try: import pysam
except: sys.exit('Pysam module not found.')

pid=str(os.getpid()+random.randint(0,999999999))

def usage():
	print """
USAGE: python FilterTable.py [options]
Options:
-i		Table file
-f		Sorted file with positions to filter in
-s		Sorted file with positions to filter out
-F		Features to filter in (separated by comma)
-S		Features to filter out (separated by comma)
-E		Exclude positions filtered out
-o		Save filtered lines to a file [stdout]
-p		Print simple statistics
-h		Print this help

"""

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:f:hs:F:S:Ep',["help"])
except getopt.GetoptError, err:
	print str(err) 
	usage()
	sys.exit()

if len(opts)==0:
	usage()
	sys.exit()
tablefile,outfile='',''
ffile,ofile='',''
save,ff,fo,exp,ps=0,0,0,0,0
infeat,outfeat=[],[]
for o,a in opts:
	if o in ("-h","--help"):
		usage()
		sys.exit()
	elif o == "-i":
		tablefile = a
		if not os.path.exists(tablefile):
			usage()
			sys.exit('Table file not found')
	elif o == "-o":
		outfile = a
		save=1
	elif o == "-s":
		ofile = a
		fo=1
		if ofile=='':
			usage()
			sys.exit('Sorted file with positions to filter out not found.')		
	elif o == "-f":
		ffile = a
		ff=1
		if ffile=='':
			usage()
			sys.exit('Sorted file with positions to filter in not found.')		
	elif o == "-F":
		infeat=[x.lower() for x in a.split(',')]
	elif o == "-S":
		outfeat=[x.lower() for x in a.split(',')]
	elif o == "-E": exp=1
	elif o == "-p": ps=1	
	else:
		assert False, "unhandled option"

# Funzioni
def filterIn(chr,exfeat,pos):
	if len(exfeat)==0: return 1
	if ff and not chr in contigf: return 0 
	elif not ff: return 1
	res=[(kk.feature).lower() for kk in tabixf.fetch(reference=chr,start=pos,end=pos+1,parser=pysam.asGTF())]
	for i in exfeat:
		if i in res: return 1
	return 0

def filterOut(chr,exfeat,pos):
	if len(exfeat)==0: return 0
	if fo and not chr in contigo: return 0
	elif not fo: return 0 
	res=[(kk.feature).lower() for kk in tabixo.fetch(reference=chr,start=pos,end=pos+1,parser=pysam.asGTF())]	
	for i in exfeat:
		if i in res: return 1
	return 0
	
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> START: %s\n" %(script_time))

if fo:
	if not os.path.exists(ofile+'.tbi'):
		sys.stderr.write('Indexing %s file.\n' %(ofile))
		ofile=pysam.tabix_index(ofile, preset='gff')
if ff:
	if not os.path.exists(ffile+'.tbi'):
		sys.stderr.write('Indexing %s file.\n' %(ffile))
		ffile=pysam.tabix_index(ffile, preset='gff')

if fo:
	tabixo=pysam.Tabixfile(ofile)
	contigo=tabixo.contigs
if ff:
	tabixf=pysam.Tabixfile(ffile)
	contigf=tabixf.contigs
	
sys.stderr.write('Reading Table file...\n')
if save: o=open(outfile,'w')
f=open(tablefile)
y,x,xx=0,0,0
for i in f:
	if i.strip()=='': continue
	if i.startswith('#'): continue
	if i.startswith('Region'):
		if save: o.write(i.strip()+'\n')
		else: sys.stdout.write(i)
		continue
	l=(i.strip('\n')).split('\t')
	xx+=1
	reg,pos = l[0],int(l[1]) # sottrarre -1 per la ricerca nella tabella
	fin=filterIn(reg,infeat,pos-1)
	fout=filterOut(reg,outfeat,pos-1)
	if fin:
		if fout:
			x+=1
			if exp: continue
			if save: o.write('#'+i)
			else: sys.stdout.write('#'+i)
		else:
			y+=1
			if save: o.write(i)
			else: sys.stdout.write(i)				
	else:
		x+=1
		if exp: continue
		if save: o.write('#'+i)
		else: sys.stdout.write('#'+i)
	
f.close()
if save: o.close()
if ff: tabixf.close()
if fo: tabixo.close()
if ps:
	sys.stdout.write("All positions: %i\n" %(xx))
	sys.stdout.write("Positions filtered in: %i\n" %(y))
	sys.stdout.write("Positions filtered out: %i\n" %(x))	
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> END: %s\n" %(script_time))

