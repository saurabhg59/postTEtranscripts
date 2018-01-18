#! /usr/bin/python

import subprocess
import numpy as np
import sys
from subprocess import call

files={}
geneValues={}
teValues={}
l1Values={}
h2bValues={}
sums={}
NEWsums={}
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

#check for rmsk or l1 file is missing right now

for filename in sys.argv[1:]:
	lineCount=0
	with open(filename,'r') as INPUT:
		for line in INPUT:
			line=line.strip()
			if (lineCount==0):
				a=line.split("\t")
				a[1]=a[1][:-2]
				a[2]=a[2][:-2]
				tempFilename=a[1]
				files[a[1]]=a[2] #treatment ---> control
				geneValues[a[1]]=[]
				teValues[a[1]]=[]
				h2bValues[a[1]]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
				geneValues[a[2]]=[]
				teValues[a[2]]=[]
				h2bValues[a[2]]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
				sums[a[1]]=0
				sums[a[2]]=0
				NEWsums[a[1]]=0
				NEWsums[a[2]]=0
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

###################################################################################################################
### Deleting lines with all 0s
###################################################################################################################

geneLength=len(geneValues[tempFilename])

emptyLines=[]

for i in range(geneLength):
	geneSum=0
	for j in geneValues:
		geneSum+=geneValues[j][i]
	if(geneSum==0):
		emptyLines.append(i)

for i in geneValues:
	geneValues[i]=(np.delete(geneValues[i],emptyLines)).tolist()

teLength=len(teValues[tempFilename])

emptyLinesTE=[]

for i in range(teLength):
	teSum=0	
	for j in teValues:
		teSum+=teValues[j][i]
	if(teSum==0):
		emptyLinesTE.append(i)

for i in teValues:
	teValues[i]=(np.delete(teValues[i],emptyLinesTE)).tolist()

# l1Length=len(l1Values[tempFilename])

# emptyLinesL1=[]

# for i in range(l1Length):
# 	l1Sum=0	
# 	for j in l1Values:
# 		l1Sum+=l1Values[j][i]
# 	if(l1Sum==0):
# 		emptyLines.append(i)

# for i in l1Values:
# 	l1Values[i]=(np.delete(l1Values[i],emptyLines)).tolist()

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
	mean+=sums[i]

mean/=float(len(sums))

factors=dict(sums)

geneLengths=len(geneValues[tempFilename])
teLengths=len(teValues[tempFilename])

for i in factors:
	factors[i]=mean/factors[i]
	for j in range(geneLengths):
		geneValues[i][j]*=factors[i]
	for k in range(teLengths):
		teValues[i][k]*=factors[i]

###################################################################################################################
### determine H2B barplot colors somehow
###################################################################################################################

# mutants=dict(files)

# with open("mutantsFile.txt",'r') as MUTANT:
# 	for line in MUTANT:
# 		line=line.strip()
# 		temp=line.split("\t")
# 		mutants[temp[0]]=[]
# 		for i in temp[1:]:
# 			mutants[temp[0]].append(i)

# for i in mutants:
# 	print i, mutants[i]

###################################################################################################################
### print to file and call R script
###################################################################################################################

for i in files:
	with open("tempViolin.csv","w+") as VIOLIN:
		VIOLIN.write("value,id\n")
		for j in teValues[i]:
			VIOLIN.write(str(j)+",2\n")
		for k in teValues[files[i]]:
			VIOLIN.write(str(k)+",1\n")

	with open("tempBar1.csv","w+") as BARPLOT1:
		BARPLOT1.write()
		for j in h2bValues[i]:


	with open("tempBar2.csv","w+") as BARPLOT2:
		BARPLOT2.write()
		for k in h2bValues[i]:


	arguements=i+"VS"+files[i]+".pdf"
	subprocess.check_call(["./test.R",arguements])
	subprocess.check_call(["rm","tempViolin.csv"])	

# add part to print the required values to a csv file which will be read by the R script

# arguements = violin plot filename <SPACE> barplot1 filename <SPACE> barplot2 filename

# call("./test.R " + arguements ) # this should call the R script and give the output filename as an arguement(ideally, needs to be tested)

# delete the temp csv file created above and loop this entire process

# need 2nd R script to make barplots here