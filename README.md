# rVAD

## Description
Matlab and Python libraries for an unsupervised method for robust voice activity detection (rVAD) or speech activity detection (SAD), as presented in [rVAD: An Unsupervised Segment-Based Robust Voice Activity Detection Method](https://www.sciencedirect.com/science/article/pii/S0885230819300920). 

The rVAD method consists of two passes of denoising followed by a VAD stage. It has been applied as a preprocessor for a wide range of applications, such as speech recognition, speaker identification, language identification, age and gender identification, self-supervised learning, human-robot interaction, audio archive segmentation, and so on.  

## Source code for rVAD: 
Source code in Matlab for rVAD (including both rVAD and rVAD-fast) is available under the [rVAD2.0](rVAD2.0/) folder. It is straightforward to use: Simply call the function vad.m. Some Matlab functions and their modified versions from the publicly available VoiceBox are included with kind permission of Mike Brookes.  

Source code in Python for rVAD-fast is available under the [rVADfast_py_2.0](rVADfast_py_2.0/) folder. 

rVAD-fast is 10+ times faster than rVAD while rVAD has superior performance. 

## Reference VAD for Aurora 2 database:
The frame-by-frame reference VAD was generated from the clean set of Aurora 2 using forced-alignment speech recognition and has been used as a 'ground truth' for evaluating VAD algorithms. Our study shows that forced-alignment ASR performs as well as a human expert labeler for generating VAD references, as detailed in [Comparison of Forced-Alignment Speech Recognition and Humans for Generating Reference VAD](https://www.isca-speech.org/archive/interspeech_2015/papers/i15_2937.pdf). Here are the generated [reference VAD for the training set](Aurora2TrainSet-ReferenceVAD.zip) and the [reference VAD for the test set](Aurora2TestSet-ReferenceVAD.zip). 

