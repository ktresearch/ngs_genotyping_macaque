# -*- coding: utf-8 -*-
"""
Created on Sat May 23 2020

@author: TANIMOTO Kousuke
"""

import re
import numpy as np
from sklearn.cluster import KMeans
import sys
import datetime

input = sys.argv
process = input[1]
outputfile = 'log.txt'
with open(outputfile,'a') as file_out:
    dt = datetime.datetime.now()
    print(process,dt,sep='\t',file=file_out)
