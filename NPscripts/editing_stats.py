import os 
import glob 
import commands

def getDistro(lines):
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
	for i in lines:
		sub=i[7].split()[0]
		nuc=eval(i[6])
		nv= nuc[n[sub[1]]]
		s[sub]+=nv
		all+=nv
	d={}
	for i in s:
		try: v=(s[i]/float(all))*100
		except: v=0.0
		d[i]=(s[i],all,v)	
	return d




atab=glob.glob('firstalu/DnaRna_*/outTable_*')[0]
ftab=glob.glob('second/DnaRna_*/outTable_*')[0]
os.system('cp %s refinedEditing' %(ftab))

o=open('editing.txt','w')
f=open('knownEditing')
for i in f: o.write(i)
f.close()

if os.path.exists('refinedEditing'):
	f=open('refinedEditing')
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
o.close()

if os.path.exists('refinedEditing'):
	st,out=commands.getstatusoutput('wc -l refinedEditing')
	fpos=int(out.split()[0])
	if fpos==0:
		sys.exit('NO FURTHER EDITING')
	st,out=commands.getstatusoutput('python subCount.py refinedEditing > sub1')
	if st!=0: sys.exit('ERROR in sub1.\n')
	st,out=commands.getstatusoutput('python subCount2.py refinedEditing > sub2')
	if st!=0: sys.exit('ERROR in sub2.\n')

alu,nonalu,nonrep,kn=[],[],[],0
f=open('editing.txt')
for i in f:
	if i.startswith('Reg'): continue
	l=(i.strip()).split('\t')
	if l[18]=='ed': kn+=1
	if l[14]=='SINE' and l[15][:3]=='Alu': alu.append(l)
	elif l[14]!='-' and l[15][:3]!='Alu': nonalu.append(l)
	elif l[14]=='-' and l[15]=='-': nonrep.append(l)
f.close()

alust=getDistro(alu)
nonalust=getDistro(nonalu)
nonrepst=getDistro(nonrep)
all=getDistro(alu+nonalu+nonrep)

f=open('editingStats.txt','w')
h=['Sub','#Alu','Alu','%Alu','#NonAlu','NonAlu','%NonAlu','#NonRep','NonRep','%NonRep','#All','All','%All']
f.write('\t'.join(h)+'\n')
for i in alust:
	r=[i,alust[i][0],alust[i][1],alust[i][2],nonalust[i][0],nonalust[i][1],nonalust[i][2],nonrepst[i][0],nonrepst[i][1],nonrepst[i][2],all[i][0],all[i][1],all[i][2]]
	r=[str(x) for x in r]
	f.write('\t'.join(r)+'\n')
f.close()
