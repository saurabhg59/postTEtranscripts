#!/usr/bin/env python

from subprocess import check_call
import sys
import numpy as np
import pandas as pd

files={}
geneValues={}
oldValues={}
youngValues={}
middleValues={}
sums={}
l1baseValues={}

###################################################################################################################
### Reading in all the .cntTable files and storing the values into hashes
###################################################################################################################

for filename in sys.argv[1:]:
	if("l1Age" in filename):	
		with open(filename,'r') as INPUT:
			data=INPUT.read().strip().split('\n')
		temp1=data[0].split('\t')
		# print(filename)
		# print(data[0])
		# print(temp1)
		temp1[1]=temp1[1][:-2]
		temp1[2]=temp1[2][:-2]
		files[temp1[1]]=temp1[2]
		geneValues[temp1[1]]=[]
		geneValues[temp1[2]]=[]
		oldValues[temp1[1]]=[]
		oldValues[temp1[2]]=[]
		middleValues[temp1[1]]=[]
		middleValues[temp1[2]]=[]
		youngValues[temp1[1]]=[]
		youngValues[temp1[2]]=[]
		sums[temp1[1]]=0
		sums[temp1[2]]=0

		for line in data[1:]:
			temp2=line.split('\t')
			if(':' in temp2[0]):
				if('old' in temp2[0]):
					oldValues[temp1[1]].append(int(temp2[1]))
					oldValues[temp1[2]].append(int(temp2[2]))
					sums[temp1[1]]+=int(temp2[1])
					sums[temp1[2]]+=int(temp2[2])
				elif('middle' in temp2[0]):
					middleValues[temp1[1]].append(int(temp2[1]))
					middleValues[temp1[2]].append(int(temp2[2]))
					sums[temp1[1]]+=int(temp2[1])
					sums[temp1[2]]+=int(temp2[2])
				elif('young' in temp2[0]):
					youngValues[temp1[1]].append(int(temp2[1]))
					youngValues[temp1[2]].append(int(temp2[2]))
					sums[temp1[1]]+=int(temp2[1])
					sums[temp1[2]]+=int(temp2[2])
			else:
				geneValues[temp1[1]].append(int(temp2[1]))
				geneValues[temp1[2]].append(int(temp2[2]))
				sums[temp1[1]]+=int(temp2[1])
				sums[temp1[2]]+=int(temp2[2])

	elif("l1_" in filename):
		with open(filename,'r') as INPUT:
			data=INPUT.read().strip().split('\n')
		temp1=data[0].split('\t')
		# print(filename)
		# print(data[0])
		# print(temp1)
		temp1[1]=temp1[1][:-2]
		temp1[2]=temp1[2][:-2]
		l1baseValues[temp1[1]]=[]
		l1baseValues[temp1[2]]=[]
		for line in data[1:]:
			temp2=line.split('\t')
			if(':' in temp2[0]):
				l1baseValues[temp1[1]].append(int(temp2[1]))
				l1baseValues[temp1[2]].append(int(temp2[2]))
				sums[temp1[1]]+=int(temp2[1])
				sums[temp1[2]]+=int(temp2[2])
# print(len(sums))
###################################################################################################################
### Deleting lines with all 0s
###################################################################################################################

geneLength=len(geneValues["TCGA-VP-A87D-01.bam"])
oldLength=len(oldValues["TCGA-VP-A87D-01.bam"])
middleLength=len(middleValues["TCGA-VP-A87D-01.bam"])
youngLength=len(youngValues["TCGA-VP-A87D-01.bam"])
l1baseLength=len(l1baseValues["TCGA-VP-A87D-01.bam"])

emptyLinesGenes=[]
emptyLinesOld=[]
emptyLinesYoung=[]
emptyLinesMiddle=[]
emptyLinesL1base=[]

for i in range(geneLength):
	geneSum=0
	for j in geneValues:
		geneSum+=geneValues[j][i]
	if(geneSum==0):
		emptyLinesGenes.append(i)

for i in range(oldLength):
	oldSum=0
	for j in oldValues:
		oldSum+=oldValues[j][i]
	if(oldSum==0):
		emptyLinesOld.append(i)

for i in range(middleLength):
	middleSum=0
	for j in middleValues:
		middleSum+=middleValues[j][i]
	if(middleSum==0):
		emptyLinesMiddle.append(i)

for i in range(youngLength):
	youngSum=0
	for j in youngValues:
		youngSum+=youngValues[j][i]
	if(youngSum==0):
		emptyLinesYoung.append(i)

for i in range(l1baseLength):
	l1baseSum=0
	for j in l1baseValues:
		l1baseSum+=l1baseValues[j][i]
	if(l1baseSum==0):
		emptyLinesL1base.append(i)

for i in geneValues:
	geneValues[i]=(np.delete(geneValues[i],emptyLinesGenes)).tolist()

for i in oldValues:
	oldValues[i]=(np.delete(oldValues[i],emptyLinesOld)).tolist()

for i in middleValues:
	middleValues[i]=(np.delete(middleValues[i],emptyLinesMiddle)).tolist()

for i in youngValues:
	youngValues[i]=(np.delete(youngValues[i],emptyLinesYoung)).tolist()

for i in l1baseValues:
	l1baseValues[i]=(np.delete(l1baseValues[i],emptyLinesL1base)).tolist()

###################################################################################################################
### Normalization part
###################################################################################################################
mean=0
for i in sums:
	mean+=sums[i]

mean/=float(len(sums))

factors=dict(sums)

geneLength=len(geneValues["TCGA-VP-A87D-01.bam"])
oldLength=len(oldValues["TCGA-VP-A87D-01.bam"])
middleLength=len(middleValues["TCGA-VP-A87D-01.bam"])
youngLength=len(youngValues["TCGA-VP-A87D-01.bam"])
l1baseLength=len(l1baseValues["TCGA-VP-A87D-01.bam"])

for i in factors:
	factors[i]=mean/factors[i]
	for j in range(geneLength):
		geneValues[i][j]*=factors[i]
	for j in range(oldLength):
		oldValues[i][j]*=factors[i]
	for j in range(middleLength):
		middleValues[i][j]*=factors[i]
	for j in range(youngLength):
		youngValues[i][j]*=factors[i]
	for j in range(l1baseLength):
		l1baseValues[i][j]*=factors[i]

###################################################################################################################
### This part will write the values for each pair to files and create boxplots (TBD)
###################################################################################################################

for i in files:
	with open(i+"_VS_"+files[i]+"_L1_AGE_PLOT_GG.csv","w+") as OUTPUT:
		OUTPUT.write("value,id\n")
		for j in l1baseValues[i]:
			OUTPUT.write(str(j)+",2\n")
		for k in l1baseValues[files[i]]:
			OUTPUT.write(str(k)+",1\n")
		for l in youngValues[i]:
			OUTPUT.write(str(l)+",4\n")
		for m in youngValues[files[i]]:
			OUTPUT.write(str(m)+",3\n")
		for l in middleValues[i]:
			OUTPUT.write(str(l)+",6\n")
		for m in middleValues[files[i]]:
			OUTPUT.write(str(m)+",5\n")
		for l in oldValues[i]:
			OUTPUT.write(str(l)+",8\n")
		for m in oldValues[files[i]]:
			OUTPUT.write(str(m)+",7\n")

	arguements=["./agePlot.R"]
	arguements.append(i+"_VS_"+files[i]+"_L1_AGE_PLOT.pdf")
	arguements.append(i+"_VS_"+files[i]+"_L1_AGE_PLOT_GG.csv")
	print(i+"_VS_"+files[i]+"_L1_AGE_PLOT.pdf")
	check_call(arguements)