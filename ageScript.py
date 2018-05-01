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

###################################################################################################################
### Reading in all the .cntTable files and storing the values into hashes
###################################################################################################################

for filename in sys.argv[1:]:
	lineCount=0
	with open(filename,'r') as INPUT:
		data=INPUT.read().strip().split('\n')
	temp1=data[0].split('\t')
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
	sums[temp2[2]]=0

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

###################################################################################################################
### Deleting lines with all 0s
###################################################################################################################

geneLength=len(geneValues[geneValues.keys()[0]])
oldLength=len(oldValues[oldValues.keys()[0]])
middleLength=len(middleValues[middleValues.keys()[0]])
youngLength=len(youngValues[youngValues.keys()[0]])

emptyLinesGenes=[]
emptyLinesOld=[]
emptyLinesYoung=[]
emptyLinesMiddle=[]

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

for i in geneValues:
	geneValues[i]=(np.delete(geneValues[i],emptyLinesGenes)).tolist()

for i in oldValues:
	oldValues[i]=(np.delete(oldValues[i],emptyLinesOld)).tolist()

for i in middleValues:
	middleValues[i]=(np.delete(middleValues[i],emptyLinesMiddle)).tolist()

for i in youngValues:
	youngValues[i]=(np.delete(youngValues[i],emptyLinesYoung)).tolist()

###################################################################################################################
### Normalization part
###################################################################################################################

for i in sums:
	mean+=sums[i]

mean/=float(len(sums))

factors=dict(sums)

geneLength=len(geneValues[geneValues.keys()[0]])
oldLength=len(oldValues[oldValues.keys()[0]])
middleLength=len(middleValues[middleValues.keys()[0]])
youngLength=len(youngValues[youngValues.keys()[0]])

for i in factors:
	factors[i]=mean/factors[i]
	for j in range(geneLengths):
		geneValues[i][j]*=factors[i]
	for j in range(oldLength):
		oldValues[i][j]*=factors[i]
	for j in range(middleLength):
		middleValues[i][j]*=factors[i]
	for j in range(youngLength):
		youngValues[i][j]*=factors[i]

###################################################################################################################
### This part will write the values for each pair to files and create boxplots (TBD)
###################################################################################################################