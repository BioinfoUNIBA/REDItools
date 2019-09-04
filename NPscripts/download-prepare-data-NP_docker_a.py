#!/usr/bin/env python 
import sys, os, time
import commands
import distutils.spawn

wdir='rna_editing_protocol'
redipath='./REDItools/'

def is_tool(name):
	wn=distutils.spawn.find_executable(name)
	if wn==None: return 1
	else: return wn

exe=['bwa','STAR','awk','bgzip','tabix','sort','gtf_splicesites','wget','python','gunzip']
nt=[]
prg={}
for i in exe:
  p=is_tool(i)
  if p==1: nt.append(i)
  prg[i]=p
if len(nt)>0:
  for i in nt:
    sys.stderr.write('Program %s NOT found\n' %(i))
  sys.exit('Install required software first.')

redirec=os.path.join(redipath,'accessory','rediportal2recoding.py')
if not os.path.exists(redipath): sys.exit('rediportal2recoding.py script not found.')
prg['redirec']='../../' + redirec.lstrip('./')

cdir=os.getcwd()
sys.stderr.write('Current directory: %s\n' %(cdir))
folder=os.path.join(cdir,wdir)
os.mkdir(folder)
sys.stderr.write('Directory %s created.\n' %(wdir))
sys.stderr.write('Entering %s\n' %(wdir))
os.chdir(folder)
