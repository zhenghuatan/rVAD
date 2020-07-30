function [vad_seg, D, Dsmth, snre_vad, pv_vad, e, segsnr]=snre_vad(dfdata, nfr10, flen, fsh10, ENERGYFLOOR, pv01, pvblk, vadThres)

% Ref:
%   Zheng-Hua Tan and Børge Lindberg
%   Low-Complexity Variable Frame Rate Analysis for Speech Recognition and Voice Activity Detection
%   IEEE Journal of Selected Topics in Signal Processing, 4(5), Oct. 2010.

% a posteriori SNR weighted engergy difference
Dexpl=18;
Dexpr=18;
% vadThres = 0.1; %0.1 best  %0.125 %18, 18, 0.1 for maximal, 0.16 for mean, 0.12 for ss-seg
Dsmth=zeros(nfr10,1); % smoothed energy distance 

% energy estimation
e=zeros(nfr10,1);
for i=1:nfr10
    for j=1:flen
        e(i)=e(i)+dfdata((i-1)*fsh10+j)*dfdata((i-1)*fsh10+j);
    end
    if e(i) <= ENERGYFLOOR
        e(i)=ENERGYFLOOR;
    end
end

segsnr=zeros(nfr10,1);
segsnrsmth=1; sign_segsnr=0;
D=zeros(nfr10,1);
postsnr=D;
snre_vad=zeros(nfr10,1);
sign_pv=0;
for i=1:nfr10
    if pvblk(i)==1 && sign_pv==0
        nstart=i;
        sign_pv=1; % a pitch segment starts 
    elseif (pvblk(i)==0 || i==nfr10) && sign_pv==1 
        nstop=i-1;    % a pitch segment ends
        if i==nfr10; nstop=i; end
        sign_pv=0;
        
        %if nstart>1 && nstop<nfr10
        %datai=dfdata((nstart-1)*fsh10+1:nstop*fsh10+flen-fsh10);
        datai=dfdata((nstart-1)*fsh10+1:(nstop-1)*fsh10+flen-fsh10);
        
        %[datai,snr]=specsub_segment(datai,8000);
        %segsnr(nstart:nstop)=10*log10(snr);
        %if isinf(segsnr(nstart+1))
        %    segsnr(nstart:nstop)=200;
        %end
        for j=nstart:nstop-1 % previously it was for j=nstart:nstop-1
            for h=1:flen
                e(j)=e(j)+datai((j-nstart)*fsh10+h)*datai((j-nstart)*fsh10+h);
            end
            if e(j) <= ENERGYFLOOR; e(j)=ENERGYFLOOR; end
        end
        e(nstop)=e(nstop-1);
        %end
        
        [eY,eI]=sort(e(nstart:nstop));
        emin=eY(floor((nstop-nstart+1)*0.1));
        for j=nstart+1:nstop
            postsnr(j) =log10(e(j))-log10(emin); % calculate a posteriori SNR within a pitch segment
            if postsnr(j)<0
                postsnr(j)=0;
            end
            D(j)=sqrt(abs(e(j)-e(j-1))*postsnr(j)); % weighted energy distance 
        end
        D(nstart)=D(nstart+1);
       
        Dexp = vertcat(ones(Dexpl,1)*D(nstart), D(nstart:nstop), ones(Dexpr,1)*D(nstop));
        for j=0:nstop-nstart
            Dsmth(nstart+j)=sum(Dexp(j+1:j+Dexpl+Dexpr));
        end
        
        % Dsmth_thres = max(Dsmth(nstart:nstop));
        Dsmth_thres=sum(Dsmth(nstart:nstop).*pv01(nstart:nstop))/sum(pv01(nstart:nstop));
        
        for j=nstart:nstop
            if Dsmth(j)>Dsmth_thres*vadThres
                snre_vad(j)=1;
            end
        end
    end
end
pv_vad=snre_vad;

nexpl=33;
nexpr=47; % 29 and 39, estimated statistically, 95% ; 33, 47 %98 for voicebox pitch
sign_vad=0;
for i=1:nfr10
    if snre_vad(i)==1 && sign_vad==0
        nstart=i;
        sign_vad=1;
    elseif (snre_vad(i)==0 || i==nfr10) && sign_vad==1
        nstop=i-1;
        if i==nfr10; nstop=i; end
        sign_vad=0;
        for j=nstart:nstop
            if pv01(j)==1
                break;
            end
        end
        pv_vad(nstart:max(j-nexpl-1,1))=0; % beyond 33 frames to the left, non speech 
        for j=0:(nstop-nstart)
            if pv01(nstop-j)==1
                break;
            end
        end
        pv_vad(nstop-j+1+nexpr:nstop)=0; % beyond 47 frames to the right, non speech
    end
end

nexpl =5; nexpr=12; % 9 and 13, estimated statistically 5%; 5, 12 %2 for voicebox pitch
sign_vad=0;
for i=1:nfr10
    if snre_vad(i)==1 && sign_vad==0
        nstart=i;
        sign_vad=1;
    elseif (snre_vad(i)==0 || i==nfr10) && sign_vad==1
        nstop=i-1;
        if i==nfr10; nstop=i; end
        sign_vad=0;
        if  sum(pv01(nstart:nstop)) > 4
            for j=nstart:nstop
                if pv01(j)==1
                    break;
                end
            end
            pv_vad(max(j-nexpl,1):j-1)=1;
            for j=0:(nstop-nstart)
                if pv01(nstop-j)==1
                    break;
                end
            end
            pv_vad(nstop-j+1:min(nstop-j+nexpr,nfr10))=1;
        end
        esegment=sum(e(nstart:nstop))/(nstop-nstart+1);
        if esegment < 0.001
            pv_vad(nstart:nstop)=0;
        end
        if sum(pv01(nstart:nstop)) <= 2
            pv_vad(nstart:nstop) = 0;
        end
    end
end

sign_vad=0;
esum=0;
for i=1:nfr10
    if pv_vad(i)==1 && sign_vad==0
        nstart=i;
        sign_vad=1;
    elseif (pv_vad(i)==0 || i==nfr10) && sign_vad==1
        nstop=i-1;
        if i==nfr10; nstop=i; end
        sign_vad=0;
        esum=esum+sum(e(nstart:nstop)); 
    end
end
eave=esum/(sum(pv_vad)+eps); % average pitch segment energy over the utterance 
sign_vad=0;
for i=1:nfr10
    if pv_vad(i)==1 && sign_vad==0
        nstart=i;
        sign_vad=1;
    elseif (pv_vad(i)==0 || i==nfr10) && sign_vad==1
        nstop=i-1;
        if i==nfr10; nstop=i; end
        sign_vad=0;
        % if sum(e(nstart:nstop))/(nstop-nstart+1)<eave*0.05
        %     pv_vad(nstart:nstop) = 0;  % detected speech segment has an energy smaller than 5% of average pitch segment energy, classify as non-speech 
        % end
        % This has an impact on long-duration recordings only. 
    end
end



sign_vad=0;
vad_seg=zeros(nfr10,2);
n_vad_seg=0;
for i=1:nfr10
    if pv_vad(i)==1 && sign_vad==0
        nstart=i;
        sign_vad=1;
    elseif (pv_vad(i)==0 || i==nfr10) && sign_vad==1
        nstop=i-1;
        sign_vad=0;
        n_vad_seg=n_vad_seg+1;
        vad_seg(n_vad_seg,:)=[nstart nstop];
    end
end
vad_seg(n_vad_seg+1:nfr10,:)=[];

