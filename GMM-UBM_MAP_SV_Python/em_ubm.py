from __future__ import division
import numpy as numpy
import pickle
import os
import sys
import numbers
import math
import code
from accStat import Collect_Stats, MAPaDapt, Loglikelihood, htkread 



def GMM_train(trnFile, Nmix, dfactor, rmd, emIter):
   #Dempster, A., Laird, N., and Rubin, D., Maximum likelihood from incomplete data via the EM algorithm, J. Roy. Stat. Soc. 39 (1977), -38.

   if isinstance(math.log(Nmix,2),  numbers.Integral) == 'False': 
       print "No. of mixture has to be power of 2!!"


    
   #Initial UBM -
   dim = numpy.size(htkread(trnFile[0]), axis=1)
   mn_x, cov_xx, w_x, n_frame = numpy.zeros(dim), numpy.zeros(dim), 1.0, 0.0
   
 
   
   for i in range(0, len(trnFile)): 
       data=htkread(trnFile[i]) 
       mn_x = mn_x + data.sum(axis=0)  #for mean
       n_frame= n_frame + numpy.size(data, axis=0)
   mn_x = numpy.matrix(mn_x / n_frame)


   for i in range(0, len(trnFile)):
       data=htkread(trnFile[i])
       cov_xx = cov_xx + sum(numpy.power((data - mn_x),2))  # for diagonal cov
   cov_xx = numpy.matrix(cov_xx/(n_frame-1))
   vFloor = numpy.array(cov_xx*0.1) # 10 % of global variance

   mn_x, cov_xx, w_x, mix = mn_x.T, cov_xx.T, numpy.matrix(w_x), 1
   #------------------------------------------------   
   print "Initial UBM"
   print ""
   print mn_x.T
   print cov_xx.T  
   
   while (1):
 
          print ""
          print "Estimating GMM parameters for %d mixtures" %(mix)
          print "================================================="
          print ""

         
          if Nmix == mix:
             dfactor =1 #[last case consider all data]

          for iter in range(0,emIter):
               ixx=(dim, mix)
               N_tot, F_tot, S_tot, frm_tot, llk_tot = numpy.zeros(mix), numpy.zeros(ixx), numpy.zeros(ixx), 0, 0.0

               for file_ in range(0, len(trnFile)):
                   data=htkread(trnFile[file_])

                   #[for randomize part]
                   if rmd == 1:
                      i_=numpy.random.permutation(len(data))
                      data=data[i_]
                   #--------------------
 
                   x_dec = [i_ for i_ in range(0, numpy.size(data, axis=0), dfactor)]
                   data=data[x_dec,:]

                   #statistics using the current model parameters i.e. expectation
                   N, F, S, llk, n_frames =Collect_Stats(data.T,  mn_x, cov_xx, w_x) #  for statistics

                   #accumulation
                   if (numpy.size(N, axis=0) == numpy.size(w_x, axis=1) ) and (numpy.size(F, axis=0) == numpy.size(mn_x, axis=0)): # not empty file     
                        N_tot, F_tot, S_tot = N_tot +  N, F_tot +  F, S_tot +  S
                        llk_tot, frm_tot = llk_tot + sum(llk), frm_tot + n_frames

                   else:
                        continue

               print "em iter: %d llk = %f" %(iter,llk_tot/frm_tot)
                   
               # update new paratmeter  
               w_x, mn_x  = N_tot / sum(N_tot), F_tot/N_tot # ml estimate 'mean'
               cov_xx = S_tot/N_tot - numpy.power(mn_x,2)   #new estimated "cov"
               cov_xx, w_x, mn_x = numpy.matrix(cov_xx), numpy.matrix(w_x), numpy.matrix(mn_x)
               cov_xx = numpy.matrix(constrain_varfloor(cov_xx, vFloor))
               

                  
          if mix >= Nmix:
             break
          else:
             mn_x, cov_xx, w_x = split_gmm(mn_x, cov_xx, w_x)
    
          mix = int(mix* 2)

   return mn_x, cov_xx, w_x          



def constrain_varfloor(sigma, vFloor):

    B=(numpy.size(sigma, axis=0), numpy.size(sigma, axis=1))
    kxx=numpy.zeros(B)                    
    for j in range(0, numpy.size(sigma, axis=0)):
        thr, y=vFloor.T[j], sigma[j]
        idx=numpy.where(y < thr)
        y[idx]=thr
        kxx[j,:]=y
    sigma=kxx

    return sigma
    


def split_gmm(gm, gv, gw):

    ndim, nmix = numpy.size(gv, axis=0), numpy.size(gv, axis=1)
    arg_max, sig_max= numpy.argmax(gv, axis=0), numpy.max(gv, axis=0)
 
    B=(ndim, nmix)
    eps= numpy.matrix(numpy.zeros(B))
    
    
    ##-matrix and array define can create problem to get access the "var' element and assign properly
    for idxx_ in range(0,numpy.size(arg_max,axis=1)): #fdim
        #print idxx_, arg_max[0,idxx_]
        eps[arg_max[0,idxx_],idxx_]=numpy.sqrt(sig_max[0,idxx_])


    mu=[]
    x1, x2 = gm - eps, gm + eps
    mu, sigma, w = numpy.column_stack((x1,x2)), numpy.column_stack(( gv, gv )) , numpy.column_stack((gw, gw)) * 0.5
    
    return mu, sigma, w

     


     
