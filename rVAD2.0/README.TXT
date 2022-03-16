Noise-robust voice activity detection algorithm (rVAD).
Version 2.0

28 Nov 2017, Zheng-Hua Tan

Usage: vad(finwav, fvad) 
       vad(finwav, fvad, opts) 
       vad(finwav, fvad, opts, vadThres). 

where finwav is the input WAVE file path and name, fvad is the output VAD file path and name, opts can be 0 for using pitch (default option) or 1 for using flatness (significantly faster at the cost of slightly reduced accuracy), and finally vadThres is the threshold for VAD. Refer to vad.m for more detailed explanation. 

The code has been tested on Matlab R2016a.

Refs:
[1] Z.-H. Tan, A.k. Sarkara and N. Dehak, "rVAD: an unsupervised segment-based robust voice activity detection method," Computer Speech and Language, vol. 59, pp. 1-21, 2020. 
[2] Z.-H. Tan and B. Lindberg, "Low-complexity variable frame rate analysis for speech recognition and voice activity detection,” IEEE Journal of Selected Topics in Signal Processing, vol. 4, no. 5, pp. 798-807, 2010.

Contact:
  Prof Zheng-Hua Tan
  Aalborg University, Denmark
  zt@es.aau.dk
  https://vbn.aau.dk/en/persons/107665  
