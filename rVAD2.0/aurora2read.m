function [data_float, fs]=aurora2read(fname)

% Read data from the Aurora2 database 

fid=fopen(fname,'r','b');
data=fread(fid,'int16');
fclose(fid);
fs=8000;

%str1=strread(fname,'%s','delimiter','.'); 

data_float = double(data)/2^15; %% Normalize int16(y) by 2^15

% wavwrite(data_float,fs,strcat(str1{1},'.wav'));

