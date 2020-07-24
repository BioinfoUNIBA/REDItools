#!/usr/bin/env python
import sys

try:
	infile=sys.argv[1]
except:
	sys.exit('USAGE: <REDItool table>')

s={}
for i in 'ACGT':
	for j in 'ACGT':
		if i!=j: s[i+j]=0
n={}
x=0
for i in 'ACGT':
	n[i]=x
	x+=1
all=0
f=open(infile)
for i in f:
	if i.startswith('Reg'): continue
	l=(i.strip()).split('\t')
	if l[7]=='-': continue
	sub=l[7].split()[0]
	s[sub]+=1
	all+=1
f.close()

for i in s:
	try: v=(s[i]/float(all))*100
	except: v=0.0
	print i,s[i],all,v

