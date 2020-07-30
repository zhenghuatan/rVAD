function [noise_samp, n_noise_samp, noise_seg, D, Dsmth, snre_vad, e]=snre_highenergy(dfdata, nfr10, flen, fsh10, ENERGYFLOOR, pv01) 

% Ref:
%   Z-H Tan and B Lindberg
%   Low-Complexity Variable Frame Rate Analysis for Speech Recognition and Voice Activity Detection
%   IEEE Journal of Selected Topics in Signal Processing, 4(5), Oct. 2010.

% square root of a posteriori SNR weighted engergy difference, block based 
%
% Modified 01 Mar 2013

Dexpl=18;
Dexpr=18;
segThres = 0.25; 
e=zeros(nfr10,1);
for i=1:nfr10
    for j=1:flen 
        e(i)=e(i)+dfdata((i-1)*fsh10+j)*dfdata((i-1)*fsh10+j);  
    end
    if e(i) <= ENERGYFLOOR
        e(i)=ENERGYFLOOR;
    end
end



emin=ones(nfr10,1);
NESEG = 200;
if nfr10 < NESEG; NESEG=nfr10; end
for i=1:floor(nfr10/NESEG)
	[eY,eI]=sort(e((i-1)*NESEG+1:i*NESEG));
	emin((i-1)*NESEG+1:i*NESEG)=eY(floor(NESEG*0.1));
    if i~=1
        emin((i-1)*NESEG+1:i*NESEG)=0.9*emin((i-1)*NESEG)+0.1*emin((i-1)*NESEG+1); 
    end
end
if i*NESEG~=nfr10
    [eY,eI]=sort(e((i-1)*NESEG+1:nfr10));
    emin(i*NESEG+1:nfr10)=eY(floor((nfr10-(i-1)*NESEG)*0.1));
    emin(i*NESEG+1:nfr10)=0.9*emin(i*NESEG)+0.1*emin(i*NESEG+1);
end




D=zeros(nfr10,1);   
postsnr=zeros(nfr10,1);
for i=2:nfr10
    postsnr(i) =log10(e(i))-log10(emin(i));
    if postsnr(i)<0
        postsnr(i)=0; 
    end 
    D(i)=sqrt(abs(e(i)-e(i-1))*postsnr(i));
end
D(1)=D(2);

Dexp = vertcat(ones(Dexpl,1)*D(1), D, ones(Dexpr,1)*D(nfr10));
Dsmth = zeros(nfr10,1);
for i=1:nfr10
    Dsmth(i)=sum(Dexp(i:i+Dexpl+Dexpr));
end

for i=1:floor(nfr10/NESEG)
    Dsmth_max((i-1)*NESEG+1:i*NESEG)=max(e((i-1)*NESEG+1:i*NESEG));
end
if i*NESEG~=nfr10
    Dsmth_max(i*NESEG+1:nfr10)=max(e((i-1)*NESEG+1:nfr10)); 
end

snre_vad = zeros(nfr10,1);
for i=1:nfr10
   if Dsmth(i)>Dsmth_max(i)*segThres; snre_vad(i)=1; end
end



% block based processing to remove noise part by using snre_vad1.
sign_vad = 0;
noise_seg=zeros(floor(nfr10/1.6),1);
noise_samp=zeros(nfr10,2);
n_noise_samp=0;
for i=1:nfr10
    if snre_vad(i) == 1 && sign_vad == 0 % start of a segment
        sign_vad = 1;
        nstart=i;
    elseif (snre_vad(i) ==0 || i==nfr10) && sign_vad == 1 % end of a segment
        sign_vad = 0;
        nstop=i-1;
        if sum(pv01(nstart:nstop))==0
            noise_seg(round(nstart/1.6):floor(nstop/1.6)) = 1;
            n_noise_samp=n_noise_samp+1;
            noise_samp(n_noise_samp,:)=[(nstart-1)*fsh10+1 nstop*fsh10];
        end
    end
end
noise_samp(n_noise_samp+1:nfr10,:)=[];

