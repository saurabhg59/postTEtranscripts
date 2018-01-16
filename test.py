#! /home/saurabhg59/anaconda2/bin/python

import numpy as np
import sys

files={}
geneValues={}
teValues={}
l1Values={}
h2bValues={}
sums={}
mean=0

H2Bindex = {"NM_001002916" : 0, "NM_170610" : 1,
"NM_021062" : 2, "NM_003526" : 3,
"NM_021063" : 4, "NM_138720" : 5,
"NM_003523" : 6, "NM_003522" : 7,
"NM_003518" : 8, "NM_003524" : 9,
"NM_003525" : 10, "NM_021058" : 11,
"NM_001312653" : 12, "NM_080593" : 13,
"NM_003519" : 14, "NM_003521" : 15,
"NM_003520" : 16, "NM_003527" : 17,
"NM_003528" : 18}

for filename in sys.argv[1:]:
	lineCount=0
	with open(filename,'r') as INPUT:
		for line in INPUT:
			line=line.strip()
			if (lineCount==0):
				a=line.split("\t")
				a[1]=a[1][:-2]
				a[2]=a[2][:-2]
				files[a[1]]=a[2]
				geneValues[a[1]]=[]
				teValues[a[1]]=[]
				h2bValues[a[1]]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
				geneValues[a[2]]=[]
				teValues[a[2]]=[]
				h2bValues[a[2]]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
				sums[a[1]]=0
				sums[a[2]]=0
				lineCount+=1
			else:
				line=line.strip()
				b=line.split("\t")
				if (":" in b[0]):
					teValues[a[1]].append(int(b[1]))
					teValues[a[2]].append(int(b[2]))
				else:
					geneValues[a[1]].append(int(b[1]))
					geneValues[a[2]].append(int(b[2]))
					for i in H2Bindex:
						if i in b[0]:
							h2bValues[a[1]][H2Bindex[i]]=int(b[1])
							h2bValues[a[2]][H2Bindex[i]]=int(b[2])


#delete all 0 value lines here



###################################################################################################################
### Normalization part
###################################################################################################################

for i in files:
	for j in geneValues[i]:
		sums[i]+=j
	for k in teValues[i]:
		sums[i]+=k
	for l in geneValues[files[i]]:
		sums[files[i]]+=l
	for m in teValues[files[i]]:
		sums[files[i]]+=m

for i in sums:
	# print i,"--->",sums[i]
	# for j in h2bValues[i]:
	# 	print j
	mean+=sums[i]

mean/=float(len(sums))

factors=sums

for i in factors:
	factors[i]=mean/factors[i]
	for j in geneValues[i]:
		j=j*factors[i]
	for k in teValues[i]:
		k=k*factors[i]

###################################################################################################################
### print to file and call R script
###################################################################################################################
