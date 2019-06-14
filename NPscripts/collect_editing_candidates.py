
import sys, os
import glob

atab=glob.glob('firstalu/DnaRna_*/outTable_*')[0] #alu refined
ftab=glob.glob('second/DnaRna_*/outTable_*')[0]
if not os.path.exists('knownEditing'): sys.exit('knownEditing file not found.')
if not os.path.exists('pos.txt'): sys.exit('pos.txt file not found.')
if not os.path.exists('posalu.txt'): sys.exit('posalu.txt file not found.')

o=open('editing.txt','w')
f=open('knownEditing')
for i in f: o.write(i)
f.close()
if os.path.exists(ftab):
	f=open(ftab)
	d={}
	for i in f:
		if i.startswith('Region'): continue
		l=(i.strip()).split('\t')
		d[(l[0],l[1])]=0
	f.close()
	f=open('pos.txt')
	for i in f:
		if i.startswith('Region'): continue
		l=(i.strip()).split('\t')
		if d.has_key((l[0],l[1])): o.write(i)
	f.close()
f=open(atab)
d={}
for i in f:
	if i.startswith('Region'): continue
	l=(i.strip()).split('\t')
	d[(l[0],l[1])]=0
f.close()
f=open('posalu.txt')
for i in f:
	if i.startswith('Region'): continue
	l=(i.strip()).split('\t')
	if d.has_key((l[0],l[1])): o.write(i)
f.close()
