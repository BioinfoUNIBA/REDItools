#!/usr/bin/env python
import sys, os
import math

try:
	pslfile=sys.argv[1]
	outfile=sys.argv[2]
except:
	sys.exit('USAGE: <psl file> <output file>')
	
#### for blat
def getPS(line):
	pid = (100.0 - (pslCalcMilliBad(line) * 0.1))
	score = pslScore(line)
	#print "The percentage:",pid
	#print "Score:",score
	return pid,score

def pslScore(cols):
	sizeMul =  1
	return sizeMul * (int(cols[0]) + (int(cols[2]))) - sizeMul * int(cols[1]) - int(cols[4]) - int(cols[6])

def round(number):
	return int(number + .5);

def pslCalcMilliBad(cols):
	sizeMul = 1
	# cols[0]  matches
	# cols[1]  misMatches
	# cols[2]  repMaches
	# cols[4]  qNumInsert
	# cols[6]  tNumInsert
	# cols[11] qStart
	# cols[12] qEnd
	# cols[15] tStart
	# cols[16] tEnd
	qAliSize = sizeMul * (int(cols[12]) - int(cols[11]))
	tAliSize = int(cols[16]) - int(cols[15])
	# I want the minimum of qAliSize and tAliSize
	if qAliSize < tAliSize: aliSize = qAliSize #? $aliSize = $qAliSize : $aliSize =  $tAliSize;
	else: aliSize = tAliSize
	# return 0 is AliSize == 0
	if aliSize <= 0: return 0
	# size diff
	sizeDiff = qAliSize - tAliSize
	if sizeDiff < 0: sizeDiff = 0  
	# insert Factor
	insertFactor = int(cols[4]) 
	# $insertFactor += $cols[6];
	milliBad = (1000 * (int(cols[1])*sizeMul + insertFactor + round(3*math.log( 1 + sizeDiff)))) / (sizeMul * (int(cols[0]) + int(cols[2])+ int(cols[1])))
	return milliBad

def com(num,list):
	for i in list:
		if i[0]<=num<=i[1]: return 1
	return 0

def min95(val,score):
	#if val < (score*95.0)/100: return 1
	if val < (score*0.95): return 1
	return 0

def readLines(lines):
	res=[]
	for line in lines:
		pidd,score=getPS(line)
		#print pidd,score
		sp=[int(x) for x in (line[18].strip(',')).split(',')]
		tstarts=[int(x) for x in (line[20].strip(',')).split(',')]
		ex=[(tstarts[x]+1,tstarts[x]+sp[x]) for x in range(len(sp))]
		nl=[line[9],score,str(int(line[11])+1),line[12],str(line[10]),pidd,line[13],line[8],int(line[15])+1,int(line[16]),ex,int(line[0])]
		res.append((int(line[0]),nl)) #score
		#if d.has_key(line[9]): d[line[9]].append((score,nl))
		#else: d[line[9]]=[(score,nl)]
	return res

def comp(ri,hits):
	g,ng=0,0
	hits.sort()
	hits.reverse()
	if len(hits)==1: #unique hit with editing candidate position included
		if hits[0][1][6]==ri[2] and com(ri[1],hits[0][1][10]): g+=1 #float(hits[0][1][5])>=90.0
		else: ng+=1
	elif len(hits)>1: #multiple hits
		if hits[0][1][6]==ri[2] and min95(hits[1][0],hits[0][0]): # if second best score less than 95% of first best score
			if com(ri[1],hits[0][1][10]): g+=1 # if first best hit include editing position
			else: ng+=1
		else: ng+=1
	if g>ng: return 1
	return 0

def readPSL(infile,outfile):
	f=open(infile)
	o=open(outfile,'w')
	name,lines,xx='',[],0
	while 1:
		line=f.readline()
		if not line:
			if name=='': break
			nn=name.split('$')
			oread=(name,int(nn[2]),nn[1])
			bread=readLines(lines)
			badr=''
			if len(bread)==0: badr=name
			else: 
				if not comp(oread,bread): badr=name
			if badr!='':
				#o.write(name[:-2]+' '+name[-1]+'\n')
				o.write(name.split('_')[0]+' '+name.split('$')[0][-1]+'\n')
				xx+=1
			break
		if line.strip()=='': continue	
		if line.startswith('psL'): continue
		if (line.strip()).startswith('match'): continue
		if line.startswith('-'): continue
		l=(line.strip()).split('\t')
		if l[9]!=name:
			if len(lines)!=0:
				nn=name.split('$')
				#(rname,pileupcolumn.pos+1,chr)
				oread=(name,int(nn[2]),nn[1]) #dread[name]
				bread=readLines(lines)
				badr=''
				if len(bread)==0: badr=name
				else: 
					if not comp(oread,bread): badr=name
				if badr!='':
					#o.write(name[:-2]+' '+name[-1]+'\n')
					o.write(name.split('_')[0]+' '+name.split('$')[0][-1]+'\n')
					xx+=1
			lines=[l]
			name=l[9]
		else: lines.append(l)
	f.close()
	o.close()
	return xx

def readgf(infile):
	f=open(infile)
	for i in f:
		if 'Server ready for queries!' in i:
			f.close()
			return 1
	f.close()
	return 0

def parse(line):
	l=(line.strip()).split('\t')
	cc=(int(l[3]),int(l[4]))
	return cc

readPSL(pslfile,outfile)

