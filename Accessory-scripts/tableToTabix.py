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

import sys, os, getopt, time, random, heapq, shutil
from tempfile import gettempdir
from itertools import islice, cycle
from collections import namedtuple
from operator import itemgetter
try: import pysam
except: sys.exit('Pysam module not found.')

version='1.0'
pid=str(os.getpid()+random.randint(0,999999999))

def usage():
	print """
USAGE: python tableToTabix.py [options]
Options:
-i		TAB-delimited file
-s		Sequence name column [1]
-c		Start column [4]
-e		End column (can be identical to -c) [5]
-m		Skip lines starting with [#]
-0		Zero based coordinates
-S		Do not sort input file (sort by default)
-b		Buffer size (as number of lines) [32000]
-t		Temporary directory to use (multiple -t may be used)
-u		Save an uncompressed GFF copy (add _copy suffix)
-h		Print this help

"""

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:Sb:t:hus:c:e:m:0",["help"])
	if len(opts)==0:
		usage()
		sys.exit(2)
except getopt.GetoptError as err:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

GFFfile=''
buffer_size=32000
tempdirs=[]
sort=1
mc=0 # save an uncompressed GFF copy, default no
scol=0 # sequence column name
bcol=3 # start column
ecol=4 # end column
schar='#' # skip lines starting with this character
zcoord=False # zero based coordinated 
for o, a in opts:
	if o in ("-h","--help"):
		usage()
		sys.exit()
	elif o == "-i":
		GFFfile=a
		outfile='.'.join(GFFfile.split('.')[:-1])+'.sorted.gff'
		if not os.path.exists(GFFfile):
			usage()
			sys.exit('GFF file not found')
	elif o == "-b": buffer_size=int(a)
	elif o == "-t": tempdirs.append(a)
	elif o == "-S": sort=0
	elif o == "-u": mc=1
	elif o == "-m": schar=a
	elif o == "-s": scol=int(a)-1
	elif o == "-c": bcol=int(a)-1
	elif o == "-e": ecol=int(a)-1
	elif o == "-0": zcoord=True
	else:
		assert False, "Unhandled Option"

Keyed = namedtuple("Keyed", ["key", "obj"])
key_=eval('lambda line : (%s)' %('line[:]'))

def gk(key,obj):
	ik=itemgetter(scol,bcol,ecol)(obj.split('\t'))
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
					if j.startswith('Region'): continue
					if j.startswith(schar): continue
					l=(j.strip()).split('\t')
					l[bcol]=int(l[bcol])
					l[ecol]=int(l[ecol])
					current_chunk2.append(l)
				current_chunk3=[]
				for j in sorted(current_chunk2, key=itemgetter(scol,bcol,ecol)):
					j[bcol]=str(j[bcol])
					j[ecol]=str(j[ecol])
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

script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stdout.write("Script time --> START: %s\n"%(script_time))
if sort:
	sys.stdout.write("Sorting GFF file...\n")
	batch_sort(GFFfile,outfile,key_,buffer_size,tempdirs)
	GFFfile=outfile
if mc:
	copyfile=GFFfile+'_copy'
	shutil.copyfile(GFFfile,copyfile)
	sys.stdout.write("A copy of uncompressed GFF file has been saved on %s.\n" %(copyfile))
sys.stdout.write("Indexing GFF file...\n")
GFFfile=pysam.tabix_index(GFFfile,seq_col=scol, start_col=bcol, end_col=ecol, zerobased=zcoord)
sys.stdout.write("Tabix file saved on %s.\n" %(GFFfile))
sys.stdout.write("Indices saved on %s.tbi.\n" %(GFFfile))
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stdout.write("Script time --> END: %s\n"%(script_time))