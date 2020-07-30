function []=vadbatch_1folder_diffpathes(wavfold,nfiles1,nfiles2)

dinwav=strcat('/data/scratch/najim/RATS/Correct_data/',wavfold,'/');
dpitch=strcat('/data/scratch2/najim/RATS/PEF/',wavfold,'/');

doutwav1=strcat('/data/scratch2/zhenghua/RATS/wav1/',wavfold);
doutwav2=strcat('/data/scratch2/zhenghua/RATS/wav2/',wavfold);
dvad=strcat('/data/scratch2/zhenghua/RATS/vad/',wavfold);

d1=dir(dinwav); 
n1=length(d1);
if nargin==1
    nfiles1=1; nfiles2=n1-2;
elseif nargin==2
    nfiles2=n1-2;
    doutwav1=strcat(doutwav1,'_',num2str(nfiles1));
    doutwav2=strcat(doutwav2,'_',num2str(nfiles1));
    dvad=strcat(dvad,'_',num2str(nfiles1));
elseif nargin==3
    doutwav1=strcat(doutwav1,'_',num2str(nfiles1));
    doutwav2=strcat(doutwav2,'_',num2str(nfiles1));
    dvad=strcat(dvad,'_',num2str(nfiles1));
end
doutwav1=strcat(doutwav1,'/')
doutwav2=strcat(doutwav2,'/')
dvad=strcat(dvad,'/')

if ~isdir(doutwav1); mkdir(doutwav1); end
if ~isdir(doutwav2); mkdir(doutwav2); end
if ~isdir(dvad); mkdir(dvad); end
if nfiles2>n1-2; nfiles2=n1-2; end
for i1=2+nfiles1:2+nfiles2
    [str1, str2]=strread(d1(i1).name,'%s%s','delimiter','.');
    finwav=strcat(dinwav,d1(i1).name) 
    fpitch=strcat(dpitch,str1{1},'.PEF');
    foutwav1=strcat(doutwav1,d1(i1).name);
    foutwav2=strcat(doutwav2,d1(i1).name);
    fvad=strcat(dvad,str1{1},'.vad');
    vad(finwav,fpitch,foutwav1,foutwav2,fvad);
end

clear all;



