function [posSum,negSum,tim,dp,mob,posDat,posErr,negDat,negErr]=hm_load_ais_nds(fullPath2file)
%Load ais data that is saved in the nds-format
%

%1st-3rd row are info header 
%4th row is legend
%5-end rows are data

dp=load('dpais.dat')*1e-9;

fid=fopen(fullPath2file,'r');
%read first header

for i = 1:4
    hdr{i} = fgets(fid);
end

IdatHdr=regexp(hdr{4},'sp');
IerrHdr=regexp(hdr{4},'sperr');
for i=1:length(IdatHdr)-1
    if any(IerrHdr==IdatHdr(i)),
        Isrt=IdatHdr(i)+6;
    else
        Isrt=IdatHdr(i)+3;
    end
    mobAll(i)=str2num(hdr{4}(Isrt:IdatHdr(i+1)-1)); %(cm^2/V/m)
end
mobAll(i+1)=str2num(hdr{4}(IdatHdr(i+1)+6:end)); %(cm^2/V/m)

mob=mobAll(1:27);

%find if opmode column is used or not
if ~isempty(strmatch(hdr{4},'opmode'))
    dataStartCol=3;
    strformat='%s%s';
else
    dataStartCol=4;
    strformat='%s%s%*s';
end

fseek(fid, 0, -1);
c= textscan(fid,[strformat,repmat('%f',1,length(mobAll))],'delimiter',',','headerLines',4);
fclose(fid);


%process start time
cated=cat(1,c{1}{1:end});
yyyy=str2num(cated(:,1:4));
mm=str2num(cated(:,6:7));
dd=str2num(cated(:,9:10));
hh=str2num(cated(:,12:13));
mi=str2num(cated(:,15:16));
ss=str2num(cated(:,18:26));

startTim=datenum(yyyy,mm,dd,hh,mi,ss);

%process stop time
cated=cat(1,c{2}{1:end});
yyyy=str2num(cated(:,1:4));
mm=str2num(cated(:,6:7));
dd=str2num(cated(:,9:10));
hh=str2num(cated(:,12:13));
mi=str2num(cated(:,15:16));
ss=str2num(cated(:,18:26));

stopTim=datenum(yyyy,mm,dd,hh,mi,ss);

% if third column is operation mode
%     data starts from 4th cloumn
%     other wise from 3rd
datRaw=cat(2,c{3:end});
% datSums=[stopTim,datRaw];
tim=stopTim; %winter time? UTC time?
posDat=datRaw(:,1:27);
posErr=datRaw(:,28:27*2);
negDat=datRaw(:,55:27+55-1);
negErr=datRaw(:,82:27+82-1);

posSum=[dp;posDat];
negSum=[dp;negDat];

posSum=[[0;tim],[0;sum(posSum(2:end,:),2)],posSum];
negSum=[[0;tim],[0;sum(negSum(2:end,:),2)],negSum];

