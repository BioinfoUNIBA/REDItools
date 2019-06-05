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

import sys, os, getopt, time, random, heapq, gzip
from tempfile import gettempdir
from itertools import islice, cycle
from collections import namedtuple
from operator import itemgetter

version='1.0'
pid=str(os.getpid()+random.randint(0,999999999))

def usage():
	print """
USAGE: python TableToGFF.py [options]
Options:
-i		Table file from REDItools
-s		Sort output GFF
-t		Tabix output GFF (requires Pysam module)
-b		Buffer size (as number of lines) [32000] (requires -s)
-T		Temporary directory (requires -s)
-o		Outfile [outTable_%s.gff]
-h		Print this help

"""%(pid)

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:o:sthT:b:",["help"])
	if len(opts)==0:
		usage()
		sys.exit(2)
except getopt.GetoptError as err:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

tablefile=''
outfile='outTable_%s.gff' %(pid)
sort=0
tabix=0
buffer_size=32000
tempdirs=[]
for o, a in opts:
	if o in ("-h","--help"):
		usage()
		sys.exit()
	elif o == "-i":
		tablefile=a
		if not os.path.exists(tablefile):
			usage()
			sys.exit('Table file not found')
	elif o == "-o": outfile=a
	elif o == "-s": sort=1
	elif o == "-t": tabix=1
	elif o == "-b": buffer_size=int(a)
	elif o == "-T": tempdirs.append(a)
	else:
		assert False, "Unhandled Option"

#Sorting code from SortGFF.py

Keyed = namedtuple("Keyed", ["key", "obj"])
key_=eval('lambda line : (%s)' %('line[:]'))

def gk(key,obj):
	ik=itemgetter(0,3,4)(obj.split('\t'))
	return key((ik[0],int(ik[1]),int(ik[2])))

def merge(key=None, *iterables):
	# based on code posted by Scott David Daniels in c.l.p.
	# http://groups.google.com/group/comp.lang.python/msg/484f01f1ea3c832d
	#print iterables
	if key is None:
		keyed_iterables = iterables
	else:
		keyed_iterables = [(Keyed(gk(key,obj), obj) for obj in iterable) for iterable in iterables]
		#print keyed_iterables
	for element in heapq.merge(*keyed_iterables):
		yield element.obj

def batch_sort(input, output, key=None, buffer_size=32000, tempdirs=None):
	if tempdirs is None:
		tempdirs = []
	if not tempdirs:
		tempdirs.append(gettempdir())
	chunks = []
	xx=0
	try:
		with open(input,'rb',64*1024) as input_file:
			input_iterator = iter(input_file)
			for tempdir in cycle(tempdirs):
				current_chunk2=[]
				for j in islice(input_iterator,buffer_size):
					l=(j.strip()).split('\t')
					l[3]=int(l[3])
					l[4]=int(l[4])
					current_chunk2.append(l)
				current_chunk3=[]
				for j in sorted(current_chunk2, key=itemgetter(0,3,4)):
					j[3]=str(j[3])
					j[4]=str(j[4])
					current_chunk3.append('\t'.join(j)+'\n')
				xx+=len(current_chunk3)
				if not current_chunk3: break
				sys.stdout.write("Loaded and sorted %i lines.\n"%(xx))
				output_chunk = open(os.path.join(tempdir,'%06i_%s'%(len(chunks),pid)),'w+b',64*1024)
				chunks.append(output_chunk)
				output_chunk.writelines(current_chunk3)
				output_chunk.flush()
				output_chunk.seek(0)
		sys.stdout.write("Merging from %i files.\n"%(len(chunks)))
		with open(output,'wb',64*1024) as output_file:
			output_file.writelines(merge(key, *chunks))
	finally:
		for chunk in chunks:
			try:
				chunk.close()
				os.remove(chunk.name)
			except Exception:
				pass
#END sorting code

script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stdout.write("Script time --> START: %s\n"%(script_time))
sys.stdout.write("Reading table...\n")
if tablefile.endswith('.gz'): f=gzip.open(tablefile,'rb')
else: f=open(tablefile)
o=open(outfile,'w')
xx=0
#chr21	10205589	C	0	12	34.75	[0, 3, 0, 9]	CT	0.75	16	28.56	[0, 16, 0, 0]	-	0.00	-
for i in f:
	if i.startswith('Region'): continue
	if i.strip()=='': continue
	l=(i.strip()).split('\t')
	strand='+'
	if l[3]=='0': strand='-'
	gffLine=[l[0],'reditoolTable','pos',l[1],l[1],'.',strand,'.',l[0]+'-'+l[1]]
	o.write('\t'.join(gffLine)+'\n')
	xx+=1
f.close()
o.close()
sys.stdout.write("Converted %i lines.\n"%(xx))
sys.stdout.write("GFF saved on %s\n"%(outfile))
if sort:
	sys.stdout.write("Sorting GFF file...\n")
	outfileS='.'.join(outfile.split('.')[:-1])+'.sorted.gff'
	batch_sort(outfile,outfileS,key_,buffer_size,tempdirs)
	outfile=outfileS
	sys.stdout.write("Sorted GFF saved on %s\n"%(outfileS))
if tabix:
	try:
		import pysam
		sys.stdout.write("Indexing GFF file...\n")
		outfileS=pysam.tabix_index(outfile, preset='gff')
		sys.stdout.write("Tabix file saved on %s.\n" %(outfileS))
		sys.stdout.write("Indices saved on %s.tbi.\n" %(outfileS))
	except: sys.exit('Pysam module not found.\nTabix indexing not available.')

script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stdout.write("Script time --> END: %s\n"%(script_time))

