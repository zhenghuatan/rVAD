function [pvblk]=pitchblockdetect(pv01, nfr10, pitch, opts)

%pitch block detection

if nfr10==length(pv01)+1
    pv01(nfr10)=pv01(nfr10-1); 
end 

if opts == 0

sign_pv=0;
for i=1:nfr10
    if pv01(i)==1 && sign_pv==0
        nstart=i;
        sign_pv=1;
    elseif (pv01(i)==0 || i==nfr10) && sign_pv==1
        nstop=i;
        if i==nfr10; nstop=i+1; end
        sign_pv=0;
        pitchseg=zeros(nstop-nstart,1);
        for j=nstart:nstop-1

            if isstring(pitch(j))
                pitchseg(j-nstart+1)=str2double(pitch(j));
            else
                 pitchseg(j-nstart+1)=pitch(j);
            end

        end
        if sum(abs(round(pitchseg-mean(pitchseg))))==0 && nstop-nstart+1>=10
            pv01(nstart:nstop-1)=0;
        end
    end
end

end %opts


sign_pv=0;
pvblk=pv01;
for i=1:nfr10
    if pv01(i)==1 && sign_pv==0
        nstart=i;
        sign_pv=1;
        pvblk(max(nstart-60,1):nstart)=1;
    elseif (pv01(i)==0 || i==nfr10) && sign_pv==1
        nstop=i;
        sign_pv=0;
        pvblk(nstop:min(nstop+60,nfr10))=1;
    end
end

