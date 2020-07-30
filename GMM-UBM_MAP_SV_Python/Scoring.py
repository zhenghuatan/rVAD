from __future__ import division
import numpy as numpy
import pickle
import os
import sys
import code
import math
from multiprocessing import Pool, cpu_count
from accStat import Collect_Stats, MAPaDapt, Loglikelihood, htkread, multi_thread



#=====================================
# 3. GMM-UBM: scoring               ||
#=====================================

print "scoring ..."

nmix=4
ubmDir= 'GMM' + str(nmix)                     #directory of "GMM-UBM" model
Tardest='MAP3_Tau10.0'                        #directory of target model
Scorefile='score.txt'                         #output file : scores
CORES=2                                       #[number of threads to be used]



Tstndx=numpy.loadtxt('Trial.lst', dtype=str)  #loading the test trial list


with open(ubmDir + '/' + 'ubm') as f:       #read ubm (once)
      print "Load ubm .. %s" %(f)
      ubm_mu, ubm_cov, ubm_w = pickle.load(f)

multi_thread(CORES, Tstndx, Tardest, Scorefile, ubm_mu, ubm_cov, ubm_w)     
  
print("score --> %s\n" %(Scorefile))
