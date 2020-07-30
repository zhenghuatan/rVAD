from __future__ import division
import numpy as numpy
import pickle
import os
import sys
import code
import math
from mfcc import  mfcc, vad_thr, cmvn, writehtk




#===============================================
# 1. Feature extraction, Vad, cmvn            ||
#===============================================

winlen, ovrlen, pre_coef, nfilter, nftt = 0.025, 0.01, 0.97, 20, 512  #[window size (sec)], [frame shift(sec)], [pre-emp coeff],
                                                                      #[no. of filter in MFCC], [N-point FFT]

opts=1  #for e-VAD

with open('feat.lst', 'r') as fin: #[[load list of speech file for feature extraction]]

    for x in fin:

          x=x.split(',');  wFile=x[0]; fFile=x[1].strip('\n')
    
          try:

               #call MFCC feature extraction subroutine
               f, E, fs=mfcc(wFile,winlen, ovrlen, pre_coef, nfilter, nftt)

               # VAD part 
               if opts == 1: 

                   f=vad_thr(f,E)       #Energy threshold based VAD [comment this  line if you would like to plugin the rVAD labels]

               elif opts == 0:

                   l=numpy.loadtxt('..corresponding vad label file');     #[Pluggin the VAD label generated by rVAD matlab]
        
                   if (len(f) - len(l)) ==1: #1-[end-frame] correction [matlab/python]
                       l= numpy.append(l,l[-1:,])
                   elif (len(f) -len(l)) == -1:
                       l=numpy.delete(l,-1)
        
                   if (len(l) == len(f)) and (len(numpy.argwhere(l==1)) !=0):
                       idx=numpy.where(l==1)
                       f=f[idx]
                   else: 
                       print "mismatch frames between: label and feature files or no voice-frame in VAD"
                       exit()



               # Zero mean unit variance  normalize after VAD  
               f=cmvn(f)

               #write the VAD+normalized features  in file 
               if not os.path.exists(os.path.dirname(fFile)): # create director for the feature file
                   os.makedirs(os.path.dirname(fFile))

               print("%s --> %s\n" %(wFile,fFile))

               writehtk(fFile, f , 0.01) 
  
          except:
               print("Fail .. %s ---> %s\n" %(wFile, fFile))        


