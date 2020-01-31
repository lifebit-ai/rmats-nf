#!/usr/bin/env python
"""
License: Proprietary
Copyright: The Jackson Laboratory (2019)
Authors: Mohan Bolisetty, Anne Deslattes Mays

Description: This function, which ultimately could be part of a large class, recursively for each file,
  appends the samples together.   A large assumption is that each file has the same number or rows.

"""

import pandas as pd
import sys

filelist=sys.argv[1]
outputrmatsmatrix=sys.argv[2]
separator=sys.argv[3]

def create_matrix(filelist):

    count=0
    data_df1 = pd.DataFrame()
    fhandle = open(filelist, 'r')
    for filename in fhandle:
       count += 1
       filename = filename.strip()
       print ( count,filename)
       if data_df1.empty:
           data_df1 = pd.read_csv(filename, header=0, index_col=0, sep=separator)
           continue
       data_df2 = pd.read_csv(filename, header=0, index_col=0, sep=separator)
       print(data_df2)
       data_df1 = pd.concat([data_df1, data_df2], axis=1)
       print("---"*5)
       print(data_df1)
       print(data_df2)
    return data_df1

final_df = create_matrix(filelist)
print (final_df)
print (final_df.to_csv(outputrmatsmatrix, sep=','))
