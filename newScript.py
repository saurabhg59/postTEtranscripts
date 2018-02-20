#!/usr/bin/python3

### This scripts takes in as input list of .cntTable files, the rmsk files have to be modified to keep only the TE lines right now, it then normalizes them based on sum of all gene+TEs+L1s
### and then prints the output to 2 csv files 1 for gene+l1 and 1 for RMSK annotations, right now these files need to be merged manually.

import subprocess
import numpy as np
import sys
import pandas as pd

l1Filenumber=0
rmskFilenumber=0

for filename in sys.argv[1:]:
	if 'l1' in filename:	
		if l1Filenumber==0:
			df1=pd.read_table(filename,sep='\t')
			df1.set_index('gene/TE',inplace=True)
			l1Filenumber+=1
		else:
			temp=pd.read_table(filename,sep='\t')
			temp.set_index('gene/TE',inplace=True)
			df1=pd.concat([df1,temp],axis=1,join_axes=[df1.index])
	elif 'rmsk' in filename:
		if rmskFilenumber==0:
			df2=pd.read_table(filename,sep='\t')
			df2.set_index('gene/TE',inplace=True)
			rmskFilenumber+=1
		else:
			temp=pd.read_table(filename,sep='\t')
			temp.set_index('gene/TE',inplace=True)
			df2=pd.concat([df2,temp],axis=1,join_axes=[df2.index])

sums=[]
for i in range(len(list(df1))):
	sums.append(np.sum(df1[list(df1)[i]])) ### calculates the sum for column i of the dataframe1
	sums[i]+=np.sum(df2[list(df2)[i]])
	
	### this is giving the same sums as the the old script

meanSum=np.mean(sums)

for i in range(len(sums)):
	sums[i]=meanSum/sums[i]
	df1[list(df1)[i]]=df1[list(df1)[i]].apply(lambda x: x*sums[i])
	df2[list(df2)[i]]=df2[list(df2)[i]].apply(lambda x: x*sums[i])

df1.to_csv('gene&L1.csv',sep=',')
df2.to_csv('TEs.csv',sep=',')