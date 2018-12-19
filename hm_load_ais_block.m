function [neg,pos,dp,mob,blockStartTim,blockStopTim]=hm_load_ais_block(fullPath2file)
%
% Load AIS block filse saved with spectops
%
%

%Heikki Junninen
%Dec 2014

nais_i=importdata(fullPath2file);
%         nais_p=importdata([dPath,dayStr,'-block-particles.spectra']);

%remove notification about the difference from UTC
%time for ions
%find the first data row
dt=cellfun(@(x) strcmp(x(1),'#'),nais_i.textdata(:,1),'UniformOutput',false);
Idat=find(cell2mat(dt),1,'last')+2;
dt=cellfun(@(x) x(1:26),nais_i.textdata(Idat:end,1),'UniformOutput',false);
blockStartTim=datenum(dt);

dt=cellfun(@(x) x(1:26),nais_i.textdata(Idat:end,2),'UniformOutput',false);
blockStopTim=datenum(dt);

mob=[3.16, 2.37, 1.78, 1.33, 1, 0.75, 0.562, 0.422, 0.316, 0.237, 0.178, 0.133, 0.1, 0.075, 0.0562, 0.0422, 0.0316, 0.0237, 0.0178, 0.0133, 0.01, 0.0075, 0.00562, 0.00422, 0.00316, 0.00237, 0.00178, 0.00133];
dp=[0.46 0.62 0.81 1.03 1.27 1.49 1.72 2.00 2.35 2.78 3.30 3.90 4.60 5.41 6.35 7.43 8.69 10.2  11.9 13.8 16.1 18.7 21.8 25.4 29.6 34.6 40.1 45]*1e-9;

neg=nais_i.data(:,1:28);
pos=nais_i.data(:,28+28+1:28+28+28);


