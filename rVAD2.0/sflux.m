function [ft, d, sVar]= sflux(data,flen,fsh10);

%% output -
% d - spectral flux
% ft - spectral flatness
% sVar - spectral variance

nftt=pow2(nextpow2(flen)); %% FFT point 

%% sf-> spectral flux, ft-> spectral flatness, sVar-> spectral variance
x=enframe(data,flen,fsh10);
w=hamming(flen);
x=x.*repmat(w',size(x,1),1); 
 
ak=abs(fft(x',nftt)); % spectrum
ak=ak'; 
ak=ak(:, 1:fix(nftt/2)+1); 

ak_1=ak(2:end,:); % ak-1
ak_1=[ ak_1 ; ak_1(end,:)]; % ak(t-1)
 
d= sum((ak - ak_1).^2, 2); % sum_k [ak(t) -ak(t-1)]
denA= sqrt(sum(ak.^2, 2)) .* sqrt( sum(ak_1.^2, 2) );
d=(d+eps)./(denA+eps);  
          
%% flatness
win=size(ak,2); % number of bands in spectra
num= exp( (1/win) * sum( log(ak),2) );
den= (1/win) * sum( ak,2);
ft= (num+eps)./(den+eps);
        
%% Spectral Variation is the normalized by the correlation of spectrum between consecutive frames
num= (sum(ak.*ak_1,2) +eps)./(denA+eps);
sVar= 1- num;
        
        

        
