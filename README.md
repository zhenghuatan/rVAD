# rVAD

## Description
This is the Matlab and Python library for an unsupervised method for robust voice activity detection (rVAD) or speech activity detection (SAD), as presented in [rVAD: An Unsupervised Segment-Based Robust Voice Activity Detection Method](https://www.sciencedirect.com/science/article/pii/S0885230819300920). 

The rVAD method consists of two passes of denoising followed by a VAD stage. It has been applied as a preprocessor for speech recognition, speaker identification, language identification, age and gender identification, human-robot interaction, audio archive segmentation, and so on. More info on [the rVAD webpage](http://kom.aau.dk/~zt/online/rVAD/). 

## Source code for rVAD: 
Source code in Matlab for rVAD (including rVAD-fast) is available under the [rVAD2.0](rVAD2.0/) folder. It is straightforward to use: Simply call the function vad.m. Some Matlab functions and their modified versions from the publicly available VoiceBox are included with kind permission of Mike Brookes.  

Source code in Python for rVAD-fast is available under the [rVADfast_py_2.0](rVADfast_py_2.0/) folder. 

## Python code for GMM-UBM based speaker verification
Source code in Python for training and testing of GMM-UBM and maximum a posterirori (MAP) adapation based speaker verification is available under the [GMM-UBM_MAP_SV_Python](GMM-UBM_MAP_SV_Python/) folder. 

## Reference VAD for Aurora 2 database:
The frame-by-frame reference VAD was generated from the clean set of Aurora 2 using forced-alignment speech recognition and has been used as a 'ground truth' for evaluating VAD algorithms. Our study shows that forced-alignment ASR performs as well as a human expert labeler for generating VAD references, as detailed in [Comparison of Forced-Alignment Speech Recognition and Humans for Generating Reference VAD](https://www.isca-speech.org/archive/interspeech_2015/papers/i15_2937.pdf). Here are the generated [reference VAD for the training set](Aurora2TrainSet-ReferenceVAD.zip) and the [reference VAD for the test set](Aurora2TestSet-ReferenceVAD.zip). 
