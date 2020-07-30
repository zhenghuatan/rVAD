from __future__ import division
import numpy as numpy
import pickle
import os
import sys
import code
import math
from multiprocessing import Pool, cpu_count
from em_ubm import GMM_train
from accStat import Collect_Stats, MAPaDapt, Loglikelihood, htkread, multi_thread



#=================================
#2. GMM-UBM training            ||  
#=================================
nmix, dsfactor, rmd, emIter =4, 10, 0, 5 #nmix=[mixture,power of 2],dfactor= decimination of frames during itermediate UBM training/file (speed up),[EM iter]
                                           # rmd =1 ;  1) randomize frames  --> 2) decimination [llh may not increasing in EM for interm. model]         
ubmDir= 'GMM' + str(nmix)  #directory where GMM has to be stored



#reading and checking the UBM training files in list
tmpScp=[]
with open('UBM.lst', 'r') as fin:  
     for x in fin:
         try:  
              if len(htkread(x.strip())) >1:
                 tmpScp.append(x.strip())
              else:
                 print("..File discarded GMM trn due to 1-frame/empty. ...%s\n" %(x)) 
         except:
              print("file does not exist.. %s\n" %(x))




#checking no of files survive for the GMM training
if len(tmpScp) ==0:
   print "no file available for GMM training"
   exit()
else:
   print "No. of files for UBM training ...%d" %(len(tmpScp))     





#call GMM training subroutine
ubm_mu, ubm_cov, ubm_w= GMM_train(tmpScp, nmix, dsfactor, rmd, emIter)



#save GMM
if not os.path.exists(ubmDir):
   os.makedirs(ubmDir)

with open(ubmDir + '/' + 'ubm', 'w') as f:
      pickle.dump([ubm_mu, ubm_cov, ubm_w], f)


