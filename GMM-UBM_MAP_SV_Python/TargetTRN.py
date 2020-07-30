from __future__ import division
import numpy as numpy
import pickle
import os
import sys
import code
import math
from accStat import Collect_Stats, MAPaDapt, Loglikelihood, htkread, multi_thread



#======================================
# 3. MAP - target model formation    ||
#======================================
nmix=4
ubmDir= 'GMM' + str(nmix)
MapItr, Tau =3, 10.0                              #[no of MAP iteration], [value of relevance factor]
Tardest='MAP' + str(MapItr) + '_Tau' + str(Tau)   #output directory of target model


if not os.path.exists(Tardest):
   os.makedirs(Tardest)



#load UBM (previous trained)
with open(ubmDir + '/' + 'ubm') as f:             
      print "load ubm .. %s" %(ubmDir + '/' + 'ubm')
      ubm_mu, ubm_cov, ubm_w = pickle.load(f)
 

print ubm_mu.shape, ubm_cov.shape, ubm_w.shape



#load the target training list file
tarList=[]; scpFile=[]

with open('target.ndx','r') as fin:
      for x in fin:
           x=x.split(',');  mFile=x[0]; fFile=x[1].strip('\n')
           tarList.append(mFile)
           scpFile.append(fFile)

#for target wise model
unQtarget=set(tarList)
for x in unQtarget:
    print x
    index = [i for i in range(len(tarList)) if x in tarList[i]]
    trnFile=[]
    for i,j in enumerate(index): 
        trnFile.append(scpFile[j])

    #call map subroutine
    m_, v_, w_=MAPaDapt(trnFile,  ubm_mu, ubm_cov, ubm_w, Tau, MapItr, 'Mean')

    #save target model
    with open(Tardest + '/' + x, 'w') as f:
          pickle.dump([m_, v_, w_], f)
