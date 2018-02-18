#!/usr/bin/python3

import subprocess
import numpy as np
import sys
import pandas as pd

filenumber=0

for filename in sys.argv[1:]:
	if filenumber==0:
		df=pd.read_table(filename,sep='\t')
		df.set_index('gene/TE',inplace=True)
		filenumber+=1
	else:
		temp=pd.read_table(filename,sep='\t')
		temp.set_index('gene/TE',inplace=True)
		df=pd.concat([df,temp],axis=1,join_axes=[df.index])

df.to_csv('test.csv',sep=',')


# sums=[]
# for i in range(len(list(DATAFRAME))):
# 	sums.append(np.sum(DATAFRAME[list(DATAFRAME)[i]])) ### calculates the sum for column i of the dataframe