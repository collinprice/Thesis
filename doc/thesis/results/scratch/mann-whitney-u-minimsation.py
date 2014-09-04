import scipy
import scipy.stats
import numpy
import os.path
import time
import shutil
import copy
import glob

import sys

from os import listdir
from os.path import isfile, join

from itertools import product, groupby

import numpy as np
from numpy import array
from numpy import array, asarray, dot, ma, zeros, sum

import scipy
from scipy.stats import chi2, t
from scipy.stats.mstats import rankdata
from scipy.stats import tiecorrect
from scipy.stats import distributions

def mann_whitney_u(x, y):
    x = asarray(x)
    y = asarray(y)
    n1 = len(x)
    n2 = len(y)
    ranked = rankdata(np.concatenate((x,y)))
    rankx = ranked[0:n1]       # get the x-ranks
    # ranky = ranked[n1:]        # the rest are y-ranks
    u1 = n1*n2 + (n1*(n1+1))/2.0 - np.sum(rankx,axis=0)  # calc U for x
    u2 = n1*n2 - u1                            # remainder is U for y
    #bigu = max(u1,u2)
    smallu = min(u1,u2)
    # T = np.sqrt(tiecorrect(ranked))  # correction factor for tied scores
    T = tiecorrect(ranked)
    #print T
    if T == 0:
        #raise ValueError('All numbers are identical in amannwhitneyu')
        z = 0
    else:
        sd = np.sqrt(T*n1*n2*(n1+n2+1)/12.0)
        z = (smallu-n1*n2/2.0) / sd  # normal approximation for prob calc
    
    return u1, u2, z, distributions.norm.sf(abs(z))  # (1.0 - zprob(z))

def getValues(file):
    return [float(f) for f in open(os.path.join(path, file)).readlines()[-1].split(" ")[1:]]

def read_floats(file):
    with open(os.path.join(path, file)) as f:
        return [float(x) for x in f]

path = sys.argv[1]
files_list = [f for f in listdir(path) if isfile(join(path,f)) ]

# dim = sys.argv[2]
# files_list = [f for f in listdir(path) if isfile(join(path,f)) and str(f).endswith(dim + '.txt') ]

print 'algorithm\twins\tlosses'
lines = {}
for a1 in files_list:
    wins = 0
    losses = 0

    n1 = str(a1)[:-4]
    a1_values = read_floats(a1)

    
    # print n1, ', ',
    for a2 in files_list:
        n2 = str(a2)[:-4]
        if n1 == n2: continue

        a2_values = read_floats(a2)
 
        # print a1_values
        # print a2_values
        # print n1 + " vs. " + n2
        mw = mann_whitney_u(a1_values, a2_values)
        print mw
        alpha = 0.05
        if float(mw[3]) < float(alpha):
            if mw[0] > mw[1]:
                wins += 1
            elif mw[0] < mw[1]:
                losses += 1

    #print wins, ', ', losses
    lines[n1] = [wins, losses]

#print lines
srtd = sorted(lines.items(), key=lambda x: x[0])

for l in srtd:
    print l[0], '\t', l[1][0], '\t', l[1][1]
