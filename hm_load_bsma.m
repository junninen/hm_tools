function strOut=hm_load_bsma(fullPath2file,day)
%
% read bsma 2E-files
%
% used by hm_load
%

%Heikki Junninen
%11.06.2007

fid=fopen(fullPath2file,'r');
hdr = fgets(fid);
Itab=regexp(hdr,'\t');
Itab=[0,Itab,length(hdr)+1];
for i=1:length(Itab)-1
    hdrs{i}=hdr(Itab(i)+1:Itab(i+1)-1);
end
fseek(fid, 0, -1);
c= textscan(fid,[repmat('%f',1,length(hdrs))],'headerLines',1);
fclose(fid);

dat=cell2mat(c);

y=floor(dat(:,1)/10000);
m=floor((dat(:,1)-y*10000)/100);
d=dat(:,1)-y*10000-m*100;
y=y+2000;

hh=floor(dat(:,2)/100);
mn=dat(:,2)-hh*100;

tim=datenum(y,m,d,hh,mn,0);
clear y m d hh mn

Itim=find(floor(tim)==day);

IdimPos=(8:17);
IdimNeg=(18:27);
ImobPos=(28:43);
ImobNeg=(44:59);

%extract header
Iviiva=regexp(hdrs(IdimPos),'-');
hdrDimPos=hdrs(IdimPos);
for i=1:length(hdrDimPos)
    dpPos1(i)=str2num(hdrDimPos{i}(3:Iviiva{i}-1))*1e-9;
    dpPos2(i)=str2num(hdrDimPos{i}(Iviiva{i}+1:end))*1e-9;
end
avrPdPos=exp(mean(log([dpPos1;dpPos2])));
dDPpos=log10(dpPos2)-log10(dpPos1);

Iviiva=regexp(hdrs(IdimNeg),'-');
hdrDimNeg=hdrs(IdimNeg);
for i=1:length(hdrDimNeg)
    dpNeg1(i)=str2num(hdrDimNeg{i}(3:Iviiva{i}(end)-1))*1e-9;
    dpNeg2(i)=str2num(hdrDimNeg{i}(Iviiva{i}(end)+1:end))*1e-9;
end
avrPdNeg=exp(mean(log([dpNeg1;dpNeg2])));
dDPneg=log10(dpNeg2)-log10(dpNeg1);

Iviiva=regexp(hdrs(ImobPos),'-');
hdrMobPos=hdrs(ImobPos);
for i=1:length(hdrMobPos)
    mobPos1(i)=str2num(hdrMobPos{i}(3:Iviiva{i}(end)-1));
    mobPos2(i)=str2num(hdrMobPos{i}(Iviiva{i}(end)+1:end));
end
avrMobPos=mean([mobPos1;mobPos2]);
dMobPos=mobPos2-mobPos1;

Iviiva=regexp(hdrs(ImobNeg),'-');
hdrMobNeg=hdrs(ImobNeg);
for i=1:length(hdrMobNeg)
    mobNeg1(i)=str2num(hdrMobNeg{i}(3:Iviiva{i}(end)-1));
    mobNeg2(i)=str2num(hdrMobNeg{i}(Iviiva{i}(end)+1:end));
end
avrMobNeg=mean([mobNeg1;mobNeg2]);
dMobNeg=mobNeg2-mobNeg1;
%====

%convert absolute concentrations to dN/dlogDP and dN/dlogZ
dimPos=[[0,0,avrPdPos];[tim(Itim),NaN(length(Itim),1),dat(Itim,IdimPos)./repmat(dDPpos,length(Itim),1)]];
dimNeg=[[0,0,avrPdNeg];tim(Itim),NaN(length(Itim),1),dat(Itim,IdimNeg)./repmat(dDPneg,length(Itim),1)];
mobPos=[[0,0,avrMobPos];tim(Itim),NaN(length(Itim),1),dat(Itim,ImobPos)./repmat(dMobPos,length(Itim),1)];
mobNeg=[[0,0,avrMobNeg];tim(Itim),NaN(length(Itim),1),dat(Itim,ImobNeg)./repmat(dMobNeg,length(Itim),1)];

mobPosAbs=[[0,0,avrMobPos];tim(Itim),NaN(length(Itim),1),dat(Itim,ImobPos)];
mobNegAbs=[[0,0,avrMobNeg];tim(Itim),NaN(length(Itim),1),dat(Itim,ImobNeg)];

% %reduced mobility
% 
% mobN=repmat(avrMobNeg,size(Itim,1),1);
% mobP=repmat(avrMobPos,size(Itim,1),1);
% 
% P0=1013.25; %mbar
% T0=273.15; %K

strOut.datPos=dimPos;
strOut.datNeg=dimNeg;

strOut.mobPos=mobPos;
strOut.mobNeg=mobNeg;

strOut.mobPosAbs=mobPosAbs;
strOut.mobNegAbs=mobNegAbs;

strOut.dpPos_low=dpPos1;
strOut.dpPos_hgh=dpPos2;
strOut.dpNeg_low=dpNeg1;
strOut.dpNeg_hgh=dpNeg2;
strOut.mobPos_low=mobPos1;
strOut.mobPos_hgh=mobPos2;
strOut.mobNeg_low=mobNeg1;
strOut.mobNeg_hgh=mobNeg2;
