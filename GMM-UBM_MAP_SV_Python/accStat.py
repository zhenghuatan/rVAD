from __future__ import division
import numpy as numpy
import sys 
import math
import struct
import code
from multiprocessing import Process
from multiprocessing.queues import Queue
import pickle
import os, sys
import shutil
# tested on python 2.7




def Collect_Stats(data,  m, v, w):
    #estimate zero, first, second order statistics 
    
    try:
       n_mixtures=numpy.size(w)  #512x1
    except:
       n_mixtures=len(w)
    
    try:
       dim=numpy.size(m, axis=0) #57x512
    except:
       dim=len(m)
    
    if dim != numpy.size(data, axis=0):
       print "Feature dimension mismatch with model!!!"
       exit() 


    # removing "nan or inf" frames from data
    idx=numpy.argwhere(numpy.isnan(data))
    data = numpy.delete(data, numpy.unique(idx[:,0]), axis=0) #axis=0 i.e. indicate rows, 1 - column
   
    idx=numpy.argwhere(numpy.isinf(data)) #inf
    data = numpy.delete(data, numpy.unique(idx[:,0]), axis=0)
  
    n_frames=numpy.size(data, axis=1)  # 57 X frames

    if n_frames == 0:
       print "avoiding empty file after remove Nan/Inf"
       n, f, s, llk = [], [], [], []

    else:   

       #compute the GMM posteriors for the given data
       gammas, llk =GaussPosterior(data, m, v, w, n_mixtures, dim, n_frames)

       #
       n=sum(gammas.T)    # over Mixtures, zero order stat
       f=numpy.matrix(data) * numpy.matrix(gammas.T) # first order stat : [feat * GAMMAS']
       s=numpy.matrix(numpy.multiply(data, data)) * numpy.matrix(gammas.T)  #%2nd order stat (feat .* feat) * GAMMAS';

    return n, f, s, llk, n_frames




def GaussPosterior(data, m, v, w, n_mixtures, dim, n_frames):
    # Calculate posterior for a given data
    

    # denominator part of gaussian pdf
    B=(2*math.pi)
    B=math.pow(B,dim/2)
    
    A=B*(numpy.sqrt(numpy.prod(v.T, axis=1)))
    a=w/numpy.matrix(A.T) 


    #---
    B=(n_mixtures, n_frames)
    gammas = numpy.zeros(B, dtype=numpy.float64)

    

    for ii in range (0, n_mixtures):
            gammas[ii,:] = GaussFun(data, a[:,ii], m[:,ii], v[:,ii],  n_frames, dim, n_mixtures)
    
    gamasum = sum(gammas)
    llk=numpy.log(gamasum)

    #normalize
    gammas=gammas/gamasum
    
    return gammas, llk



    
def GaussFun(data, aa, b, c, n_frames, dim, n_mixtures):
    #fit data to pdf
    auxC = -0.5/c
    aux = numpy.subtract(data.T, b.T)
    aux = numpy.power(aux, 2)         
    aux=aux.T
    Y = numpy.matrix(auxC.T) * aux    

     
    Y = numpy.exp(Y) * numpy.asscalar(aa)

    return Y



def MAPaDapt(ScpFileList,  m, v, w, Tau, MapItr, parameter): # only mean adaptation implemented
    #Douglas A. Reynolds, Thomas F. Quatieri, and Robert B. Dunn, "Speaker Verification Using Adapted Gaussian Mixture Models",Digital Signal Processing 10, 19-41 (2000)
    
    for iter in range(0, MapItr):

          #initialize
          ixx=(numpy.size(w))
          N_tot=numpy.zeros(ixx)
 
          ixx=(numpy.size(m, axis=0), numpy.size(w))
          F_tot, S_tot = numpy.zeros(ixx), numpy.zeros(ixx)
          frm_tot, llk_tot = 0, 0.0 #for average llk
          

          
          for file_ in  ScpFileList: # number of files
                
                data=htkread(file_)
                data=data.T  # dim x frames

                N, F, S, llk, n_frames =Collect_Stats(data,  m, v, w) #  for statistics
                print "MAP %d ...%s, frame %d  llh %f" %(iter, file_, n_frames, sum(llk)/len(llk))
                
                if (numpy.size(N, axis=0) == numpy.size(w) ) and (numpy.size(F, axis=0) == numpy.size(m, axis=0)): # avoid empty file     

                     N_tot, F_tot, S_tot = N_tot +  N, F_tot +  F, S_tot +  S
                     llk_tot, frm_tot = llk_tot + (sum(llk)/len(llk)), frm_tot + n_frames
                     

                else:
                     continue

           #all stat estimated above
          

          if parameter == 'Mean':
                alpha, mean_ml = N_tot/(N_tot + Tau), F_tot/N_tot # ml estimate 'mean'
                mu1, mu2=numpy.multiply(m, (1-alpha)), numpy.multiply(mean_ml, alpha)   
                m_new=mu1 + mu2

          else:
             print "code only support - mean adaptation!"
             print "config parameter has to be 'm' -program terminated!"
             exit()

          # replace the old UBM parameters with re-estimated one
          m=m_new
 
    return m, v, w




def Loglikelihood(data,  m, v, w):
 
    N, F, S, llk, n_frames =Collect_Stats(data,  m, v, w)

    return llk       


def htkread(Filename):  
    
    fid=open(Filename,'rb')
    header = fid.read(12)
    (htk_size, htk_period, vec_size, htk_kind) = struct.unpack('>iihh', header) #big endean data format        
    data = numpy.fromfile(fid, dtype='f')
    param = data.reshape((htk_size, int(vec_size / 4))).byteswap()
 
    return param 


def multi_thread(CORES, Tstndx, Tardest, Scorefile, ubm_mu, ubm_cov, ubm_w):
    #for scoring split the job among the different threads
    

    if CORES <=0:
       print "No of cores can't be <=0"
       exit()

    if Tstndx.size == 1: #for single trial only
       y=(numpy.array_str(Tstndx)).split(',')
       data=htkread(y[1])
       with open( Tardest + '/' + y[0]) as f: #read target model
                      m_, v_, w_ = pickle.load(f)


       llk_ubm, llk_tr=Loglikelihood(data.T,  ubm_mu, ubm_cov, ubm_w), Loglikelihood(data.T,  m_, v_, w_)
       llr=(sum(llk_tr)  - sum(llk_ubm))/numpy.size(llk_tr, axis=0)
       txt= y[0] + ' ' + y[1]  +  ' ' + str(llr)
       fp=open(Scorefile,'w')
       fp.write(txt)
       return 


    tmp='tmp'

    if not os.path.exists(tmp):
       os.makedirs(tmp)

    queues = [Queue() for i in range(CORES)]
    part=[]
    if (Tstndx.size >= CORES) and (CORES > 1):
       part=range(0,Tstndx.size, numpy.int(Tstndx.size/CORES))
       part[0]=-1
  
       if len(part) < (CORES+1):
          part=numpy.append(part, Tstndx.size-1)
       elif part[-1] != Tstndx.size-1:
          part[-1]=Tstndx.size-1  

       args = [(part[i]+1,part[i+1], Tstndx, Tardest, i, tmp, ubm_mu, ubm_cov, ubm_w, True, queues[i]) for i in range(CORES)]
    else:
       args = [(0,Tstndx.size-1,Tstndx, Tardest, i, tmp, ubm_mu, ubm_cov, ubm_w, True, queues[i]) for i in range(CORES)]


    jobs = [Process(target=Scoring, args=(a)) for a in args]
    for j in jobs: j.start()
    for j in jobs: j.join()
    
    nt=0 
    with open(Scorefile, 'w') as outfile:
        for j in range(len(jobs)): 
            tx=tmp+ '/part' + str(j) + '.score'
            print tx
            with open(tx) as infile:
                 x=infile.read()
                 outfile.write(x)
                 x=x.count("\n")               
                 nt=nt+x
            os.remove(tx)
    outfile.close()     
    if os.path.isdir(tmp):
       shutil.rmtree(tmp)         
    print("No. of trial %d  --no. of score %d\n" %(Tstndx.size, nt))


def Scoring(iStart, iEnd, Tstndx, Tardest, j, tmp, ubm_mu, ubm_cov, ubm_w, multi=False, queue=0):


    with open(tmp + '/part' + str(j) + '.score', 'w') as f2:
         for i in range(iStart,iEnd+1): 

                x_=(Tstndx[i]).split(',')
                with open( Tardest + '/' + x_[0]) as f: #read target model
                      m_, v_, w_ = pickle.load(f)
        

                try:
  
                      data=htkread(x_[1])
                      llk_ubm, llk_tr=Loglikelihood(data.T,  ubm_mu, ubm_cov, ubm_w), Loglikelihood(data.T,  m_, v_, w_)
 
                      if  numpy.size(llk_ubm, axis=0) != numpy.size(llk_tr, axis=0):
                           print "number of scores from target model differs from ubm!!!"
                           print x_[1]
                           exit()

                      llr=(sum(llk_tr)  - sum(llk_ubm))/numpy.size(llk_tr, axis=0) 

                      txt= x_[0] + ' ' + x_[1]  +  ' ' + str(llr)
                      f2.write("%s\n" %(txt))
                except:
                      print("File avoided due to error ..%s\n" %(x_[1]))
    f2.close()
     
       

