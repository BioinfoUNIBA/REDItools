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

## {{{ http://code.activestate.com/recipes/576755/ (r3)
# based on Recipe 466302: Sorting big files the Python 2.4 way
# by Nicolas Lehuen

# Works on python 2.7+ no 3.x

import sys, os, getopt, heapq, time, random
from tempfile import gettempdir
from itertools import islice, cycle
from collections import namedtuple
from operator import itemgetter

version='1.0'
pid=str(os.getpid()+random.randint(0,999999999))

def usage():
	print """
USAGE: python SortGFF.py [options]
Options:
-i		GFF file
-o		Sorted output file [GFF_sorted_%s]
-b		Buffer size (as number of lines) [32000]
-t		Temporary directory to use (multiple -t may be used)
-h		Print this help

"""%(pid)

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:o:b:t:h",["help"])
	if len(opts)==0:
		usage()
		sys.exit(2)
except getopt.GetoptError as err:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

GFFfile=''
outfile='GFF_sorted_%s' %(pid)
buffer_size=32000
tempdirs=[]
for o, a in opts:
	if o in ("-h","--help"):
		usage()
		sys.exit()
	elif o == "-i":
		GFFfile=a
		if not os.path.exists(GFFfile):
			usage()
			sys.exit('GFF file not found')
	elif o == "-o": outfile=a
	elif o == "-b": buffer_size=int(a)
	elif o == "-t": tempdirs.append(a)	
	else:
		assert False, "Unhandled Option"

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

script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stdout.write("Script time --> START: %s\n"%(script_time))
batch_sort(GFFfile,outfile,key_,buffer_size,tempdirs)
sys.stdout.write("Sorted GFF saved on %s\n"%(outfile))
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stdout.write("Script time --> END: %s\n"%(script_time))