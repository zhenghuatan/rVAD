function [pv01, fx]=pitchestm(data, fs, nfr10, pv01)
% function [pv01, pvblk, pvblkb, pv]=pitchestm(data, fs, nfr10, pv01)

[fx,tt,pv,fv]=fxpefac(data, fs); % should plus 3
npv=length(pv);
pv01=zeros(nfr10,1); sign_pv=0;
for i=1:npv
    if pv(i)>0.25  
        pv01(i+3) =1;    
        if sign_pv==0
            sign_pv=1;
            nstart=i;
        end
    else
        if sign_pv==1
            sign_pv=0;
            nstop=i-1;
            if nstop-nstart<3
                pv01(nstart+3:nstop+3)=0;  % Remove 2 frames only pitch
            end
        end
    end
end
pv01(1:3)=pv01(4);
fxtmp(1:3)=fx(4); fxtmp(4:npv+3)=fx(1:npv);
if (npv+3) < nfr10
    pv01(npv+4:nfr10)=pv01(npv+3);
    fxtmp(npv+4:nfr10)=fx(npv);
else
    pv01=pv01(1:nfr10);
    fxtmp=fxtmp(1:nfr10);
end
fx=fxtmp;

