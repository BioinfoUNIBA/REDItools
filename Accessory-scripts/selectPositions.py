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

import sys, getopt, os, time, random, gzip

version='1.0'
pid=str(os.getpid()+random.randint(0,999999999))

def usage():
	print """
USAGE: python selectPositions.py [options]
Options:
-i		Table file from REDItools
-d		Base distribution column for DNA-Seq (-1: no DNA-Seq) [-1]
-c		Coverage RNA-Seq [5]
-C		Coverage DNA-Seq [5]
-v		Bases supporting RNA-Seq variation [1]
-V		Bases supporting DNA-Seq variation [0]
-s		Substitutions to select in RNA-Seq (separated by comma AG,CT) [all]
-f		Frequency of variation in RNA-Seq [0.1]
-F		Frequency of non-variation in DNA-Seq [0.95]
-e		Exclude multiple substitutions in RNA-Seq
-r		Exclude invariant sites in RNA-Seq
-R		Exclude variant sites in DNA-Seq #
-u		Use only positions supported by DNA-Seq
-o		Save selected positions on outTable_%s
-h		Print this help

"""%(pid)

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:c:C:v:s:f:F:euo:hrd:RV:",["help"])
	if len(opts)==0:
		usage()
		sys.exit(2)
except getopt.GetoptError as err:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

tablefile=''
outfile='outTable_%s' %(pid)
#rna-seq
cov=5
bvar=1
sfreq=0.1
expos=0
upos=0
exinv=0
subs=[x+y for x in 'ACGT' for y in 'ACGT' if x!=y]
#dna-seq
dnacol=11
dnacols=[x for x in range(dnacol-2,dnacol+3,1)]
isdna=0
gcov=5
gsfreq=0.95
gexvar=0
gbvar=0

for o, a in opts:
	if o in ("-h","--help"):
		usage()
		sys.exit()
	elif o == "-i":
		tablefile=a
		if not os.path.exists(tablefile):
			usage()
			sys.exit('Table file not found')
	elif o == "-c": cov=int(a)
	elif o == "-C": gcov=int(a)
	elif o == "-v": bvar=int(a)
	elif o == "-V": gbvar=int(a)
	elif o == "-s": subs=[x.upper() for x in a.split(',') if x.strip()!='']
	elif o == "-f": sfreq=float(a)
	elif o == "-F": gsfreq=float(a)	
	elif o == "-e": expos=1
	elif o == "-u": upos=1
	elif o == "-r": exinv=1
	elif o == "-R": gexvar=1	
	elif o == "-d":
		dnacol=int(a)
		if dnacol>3:
			isdna=1
			dnacols=[x-1 for x in range(dnacol-2,dnacol+3,1)]
	elif o == "-o": outfile=a
	else:
		assert False, "Unhandled Option"

def isnvar(nuc,idx,val):
	n=eval(nuc)
	x=0
	for j in range(4):
		if j!=idx and n[j]>=val:
			x+=1
	if x>0: return 1
	return 0

def isnvar2(nuc,idx,val):
	n=eval(nuc)
	x=0
	for j in range(4):
		if j!=idx: x+=n[j]
	if x<=val: return 1
	return 0

def issub(osubs,esubs):
	x=0
	for i in osubs:
		if i in esubs: x+=1
	if x>0: return 1
	return 0

def vinv(nuc,idx,val):
	n=eval(nuc)
	try: v=float(n[idx])/sum(n)
	except: v=0.0
	if v>=val: return 1
	return 0
	
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stdout.write("Script time --> START: %s\n"%(script_time))
sys.stdout.write("Reading table...\n")

if tablefile.endswith('.gz'): f=gzip.open(tablefile,'rb')
else: f=open(tablefile)
db={'A':0,'C':1,'G':2,'T':3}
o=open(outfile,'w')
xx,yy=0,0
for i in f:
	if i.startswith('Region'):
		o.write(i)
		continue
	if i.strip()=='': continue
	l=(i.strip()).split('\t')
	xx+=1
	if l[2] not in 'ACGTacgt': continue
	if exinv and l[7]=='-': continue
	if int(l[4])<cov: continue
	if not isnvar(l[6],db[l[2]],bvar): continue
	if l[7]!='-':
		osubs=[x.upper() for x in l[7].split()]
		if expos and len(osubs)>1: continue
		if not issub(osubs,subs): continue
	if float(l[8])<sfreq: continue
	#DNA-Seq
	if upos and not isdna: continue
	if upos and isdna and l[dnacols[0]]=='-': continue
	if isdna and l[dnacols[0]]!='-':
		if int(l[dnacols[0]])<gcov: continue
		if not isnvar2(l[dnacols[2]],db[l[2]],gbvar): continue
		if not vinv(l[dnacols[2]],db[l[2]],gsfreq): continue
		else: l[dnacols[3]]='$'
		if gexvar and l[dnacols[3]]!='-': continue #rivedere questa opzione
	o.write(i)
	yy+=1	
f.close()
o.close()
sys.stdout.write("Total lines: %i\n"%(xx))
sys.stdout.write("Filtered in lines: %i\n"%(yy))
sys.stdout.write("Selected lines saved on %s\n"%(outfile))
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stdout.write("Script time --> END: %s\n"%(script_time))

