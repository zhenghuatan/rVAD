from __future__ import division
import numpy
import sys
import os
import math
import struct
import scipy.io.wavfile as wav
from scipy.fftpack import dct
from scipy.signal import lfilter
import code


def mfcc(file_wav, winlen, ovrlen, pre_coef, nfilter, nftt):
    #	S. Davis ; P. Mermelstein, "Comparison of parametric representations for monosyllabic word recognition in continuously spoken sentences", IEEE Transactions on Acoustics, Speech, and Signal Processing ( Volume: 28, Issue: 4, Aug 1980 )
    
     fs, speech = speech_wave(file_wav) 

     eps = numpy.finfo(float).eps

     #for simple  - enegry threshold vad 
     Espeech= enframe(speech, fs, winlen, ovrlen) #framing (before pre-emphasis)
     Espeech= 20*numpy.log10(numpy.std(Espeech, axis=1, ddof=1)  + eps)
     #====
 
     speech = numpy.append(speech[0],speech[1:]-pre_coef*speech[:-1]) #pre-emphasis
     speech = enframe(speech, fs, winlen, ovrlen) #framing

     if numpy.size(speech, axis=0) != numpy.size(Espeech, axis=0):
        print "Mismatch frame numbers for  enegry and feature vectors"


     w = numpy.matrix(numpy.hamming(int(fs*winlen)) )
     w = numpy.tile(w,(numpy.size(speech, axis=0), 1))

     speech = numpy.multiply (speech, w) #apply window


     ff=(fs/2)* (numpy.linspace(0, 1 , int(nftt/2 +1) ))
     fmel=2595*numpy.log10(1+ ff/700) #mel-scale
     fmelmax, fmelmin = numpy.max(fmel), numpy.min(fmel)
     
     filtbankMel= numpy.linspace(fmelmin,fmelmax, nfilter+2) #define filter in mel domain
     filbankF=700*( numpy.power(10, (filtbankMel/2595)) -1)

     #fft
     ffy=numpy.abs(numpy.fft.fft(speech,nftt))     
     ffy, idx =numpy.power(ffy, 2), range(1,int(nftt/2) +1)
     ffy=ffy[:,idx]
 
     BB=(len(ff), nfilter)
     fbank=numpy.zeros(BB)

     for nf in range(0, nfilter):
          fbank[:,nf] = trimf(ff, filbankF[nf], filbankF[nf+1], filbankF[nf+2])
     
     #discard "filter bank energy" 
     fbank = fbank[1:] 
     fbnkSum = numpy.matrix(ffy) * numpy.matrix(fbank)  
     
     #dct
     fbankSum_eps = numpy.log10(fbnkSum.T + eps)   
     t=(dct(fbankSum_eps.T, norm = 'ortho')).T
     t= t[1:]  #dicard "c0"
     
     #rasta filtering
     t=rastaFilter(t).T

      
     # d, dd
     d=delta(t.T, 3).T         
     dd=delta(d.T, 3).T
    
 
     return numpy.column_stack((t,d,dd)), Espeech, fs
    


def speech_wave(fileName_):

     (fs,sig) = wav.read(fileName_)
     sig=sig/numpy.amax(numpy.abs(sig)) #normalize the signal (for the feature extraction)
     return fs, sig
 
def enframe(speech, fs, winlen, ovrlen):
     #split the speech data into frames  
     N, flth, foVr = len(speech), int(numpy.fix(fs*winlen)),  int(numpy.fix(fs*ovrlen))
     
     if len(speech) < flth:
        print "speech file length shorter than 1-frame"
        exit()
     

     frames = int(numpy.ceil( (N - flth + foVr)/foVr))
     slen = (frames-1)*foVr + flth


     if len(speech) < slen:
        signal = numpy.concatenate((speech, numpy.zeros((slen - N))))

     else:
        signal = deepcopy(speech)
  

     idx = numpy.tile(numpy.arange(0,flth),(frames,1)) + numpy.tile(numpy.arange(0,(frames)*foVr,foVr),(flth,1)).T
     idx = numpy.array(idx,dtype=numpy.int64)
    
 
     return signal[idx]



def trimf(Xx, aA_, bB_, cC_):
    
    if aA_ > bB_:
          print "Parameter: a > b"
          exit()
    elif  bB_ > cC_:
          print "Parameter: b > c"
          exit()
    elif  aA_ > cC_:
          print "Parameter: a > c"
          exit()
    
    BB=len(Xx) ## ff
    ky = numpy.zeros(BB)

    index=numpy.where( (Xx <= aA_) | (cC_<= Xx))     
    ky[index] = numpy.zeros(len(index[0]))            

    # slope 1
    if aA_ != bB_:
       index = numpy.where((aA_ < Xx) & (Xx < bB_))
       ky[index] = (Xx[index]- aA_)/(bB_ - aA_)

    # slope 2
    if bB_ != cC_:
       index = numpy.where( (bB_ < Xx) & (Xx < cC_))
       ky[index] = (cC_ - Xx[index])/(cC_ - bB_)


    #Center 
    index = numpy.where(Xx == bB_)
    ky[index] = numpy.ones(len(index[0]))  
    
    return ky  
 

def rastaFilter(Yx):
    #H. Hermansky and N. Morgan, "RASTA processing of speech", IEEE Trans. on Speech and Audio Proc., vol. 2, no. 4, pp. 578-589, Oct. 1994. 
     
    numer = numpy.array([i_ for i_ in range(-2,2+1, 1) ])  
    numer, denom = (-1.0*numer)/ sum(numpy.power(numer,2)), numpy.array([1, -0.94])
    

    BB=(numpy.size(Yx, axis=0),4)
    z, y =numpy.zeros(BB), numpy.zeros(BB)
  
    for fdim in range(0, numpy.size(Yx, axis=0)):
       y[fdim,:], z[fdim,:] = lfilter(numer, 1, Yx[fdim,0:4].T, axis=-1, zi=[0, 0, 0, 0]) 

    BB=(numpy.size(Yx, axis=0), numpy.size(Yx, axis=1) -4)
    txx=numpy.zeros(BB)

    for fdim in range(0, numpy.size(Yx, axis=0)):
        tmp=lfilter(numer, denom, Yx[fdim,4:numpy.size(Yx, axis=1)].T, axis=-1, zi=z[fdim,:])  
        txx[fdim,:]=tmp[0].T
        
    
    return numpy.column_stack((y*0,txx))

def delta(xY_, wx_):

    nr, nc = xY_.shape

    hlen = int(math.floor(wx_/2))
    w, win, nc = int(2*hlen) + 1,  range(hlen,-hlen-1, -1), numpy.size(xY_, axis=1)

    y1, y2 = numpy.tile(numpy.matrix(xY_[:,0]).T, (1,hlen)), numpy.tile( numpy.matrix(xY_[:,-1]).T, (1,hlen))
    yy=numpy.column_stack((y1, xY_, y2) )
    
    d = lfilter(win, 1, yy, axis=-1, zi=None)
    d = d[:,2*hlen + numpy.array(range(0,nc)) ]

    return d


def cmvn(px):

    _mean, _std = numpy.mean(px, axis=0, dtype='float64'),numpy.std(px, axis=0, ddof=1, dtype='float64') 
    
    return (px -_mean)/_std


def writehtk(filename, data, fp):

    htk_size, fperiod, fdim, paramKind = numpy.size(data, axis=0), numpy.round(fp*1.E7), (numpy.size(data, axis=1) *4), 6
    
    fid=open(filename,'wb')
    fid.write(struct.pack(">iihh", htk_size, fperiod, fdim, paramKind)) # ">" big endian
    numpy.array(data, dtype="f").byteswap().tofile(fid)
    fid.close()
          

def vad_thr(feat, E):

    idx=numpy.logical_and(E > numpy.max(E)-30,  E> -55 )
    feat=feat[idx,:]    

    return feat




