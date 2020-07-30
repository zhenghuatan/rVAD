function []=vad(finwav, fvad, opts, vadThres)

% Usage: vad(finwav, fvad) 
%        vad(finwav, fvad, opts) 
%        vad(finwav, fvad, opts, vadThres). 
%
% finwav: The input WAVE file path and name.
%
% fvad: The output VAD file path and name [optional]. If the output is in 0-1 format, each line in the file is the label for that frame (0 for non-speech and 1 for speech), while if the output is in the segment format, each line contains the start frame number and the end frame number for a speech segment. The default is 0-1 format, and one can switch to the segment format by choosing another line of fprintf in the end of this code. The frame shift is 10ms.
%
% opts: 0 for using pitch (default option), and 1 for using flatness (significantly faster at the cost of slightly reduced accuracy). 
%
% vadThres: The threshold for VAD. The default value is 0.4. Increasing vadThres (e.g. to 0.5) makes the VAD more aggressive, i.e. the number of frames to be detected as speech will be reduced.
%
% Refs:
%  [1] Z.-H. Tan, A.k. Sarkara and N. Dehak, "rVAD: an unsupervised segment-based robust voice activity detection method," Computer Speech and Language, 2019. 
%  [2] Z.-H. Tan and B. Lindberg, "Low-complexity variable frame rate analysis for speech recognition and voice activity detection,” IEEE Journal of Selected Topics in Signal Processing, vol. 4, no. 5, pp. 798-807, 2010.
%
% 2017-11-28, Zheng-Hua Tan

if nargin < 2; error('Usage: vad(finwav, fvad)'); end
if nargin == 2
  opts = 0; vadThres = 0.4; 
elseif nargin == 3
  vadThres = 0.4; 
end

[data,fs]= audioread(finwav);
% [data,fs]=wavread(finwav);
% [data, fs]=aurora2read(finwav);

% Parameter setting
ENERGYFLOOR = exp(-50);
flen=floor(fs/40); % 25ms frame length 
fsh10=fs/100; % 10ms frame shift
nfr10=floor((length(data)-(flen-fsh10))/fsh10);

b=[0.9770   -0.9770]; a=[ 1.0000   -0.9540];
fdata=filter(b,a,data);

if opts == 0
  [pv01, pitch]=pitchestm(data, fs, nfr10);
else              % using flatness 
  ftThres = 0.5;  % Default threshold. It can range from 0 to 1. Increasing ftThres increases the number of frames being detected as speech.
  [ft]= sflux(data,flen,fsh10);
  pv01 = (ft <= ftThres);  % <= threshold would give  1( meaning a speech frame)
  pitch=ft;
end

pvblk=pitchblockdetect(pv01, nfr10, pitch, opts);

[noise_samp, n_noise_samp, noise_seg]=snre_highenergy(fdata, nfr10, flen, fsh10, ENERGYFLOOR, pv01);

%% Set high energy segments to zero 
for i=1:n_noise_samp
    fdata(noise_samp(i,1):noise_samp(i,2)) = 0;
end

[dfdatarm]=specsub(fdata,fs);
% [dfdatarm]=specsub(fdata,fs,noise_seg,pv01);

[vad_seg]=snre_vad(dfdatarm, nfr10, flen, fsh10, ENERGYFLOOR, pv01, pvblk, vadThres);

%% Output VAD results in 0-1 format (1 for speech frames and 0 for non-speech ones) 
if isempty(vad_seg) ==1
   z=zeros(nfr10,1);
else
   y=[];
   for i=1:size(vad_seg,1)
       y=[ y ; [ vad_seg(i,1):vad_seg(i,2)]' ];
   end
   z=zeros(nfr10,1);
   z([y],1)=1;

   if sum(z) ~= size(y,1) % checking
      error('The number of labeled speech frames does not matched the results of detected speech segments!');
   end
end

fid=fopen(fvad,'w');
fprintf(fid, '%d\n',z'); % 0-1 VAD output
% fprintf(fid, '%d\n',vad_seg); % segment-label VAD output
fclose(fid);

