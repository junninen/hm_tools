function [hmD] = hm_regrid(hmD,newDp,newTim,srs)
%
% interpolate data to new grid
%
% hmD = hm_regrid(hmD,newDp,newTim,[srs])
%
%INPUTS ([]'s are optional)::
%     hmD (struct) data structure loaded with hm_load
%   newDp (double) new particle size vector
% newTime (double) new time vector
%   [srs] (string) name of source, default is all instruments
%                 different instrument names see hm_load help.
%
% interpolates data to new grid (0.5h, 20bins)
%
% See also: hm_load

% Heikki Junninen
% Oct 2016

makeTim=0;
makeDp=0;
if nargin<=3
    srs='dmps';
end
if nargin<=2
    newTim=[];
    makeTim=1;
end
if nargin==1
    newDp=[];
    makeDp=1;
end

%interpolate data to new grid
eval(['datC=hmD.',srs,';'])

nrDays=length(datC);
for d=1:nrDays
    
    dP=datC{d}(1,3:end);
    if makeDp
        nrDp=30;
        %     newDp=10.^(log10(3.1e-9):(log10(1000e-9)-log10(3e-9))/(nrDp+1):log10(500e-9));
        newDp=logspace(log10(3.1e-9),log10(700e-9),nrDp);
    end
    
    
    tim=datC{d}(2:end,1);
    if makeTim
        newTim=[(10/60)/24:10/60/(24):ceil(1)]';
        newTim(end)=[];
        newTim=newTim+floor(tim(1));
    end
    dmpsT=repmat(tim,1,length(dP));
    ltm=length(tim);
    dum=datC{d}(2:end,3:end);
    %     dum=medfilt1(datC{1}(2:end,3:end),13);
    dmpsi = griddata(dmpsT,repmat(log10(dP),ltm,1),dum,...
        repmat(newTim,1,length(newDp)),repmat(log10(newDp),length(newTim),1),'linear');
    
    %     save regrided data to structure
    dmps=[newTim,zeros(length(newTim),1),dmpsi];
    dmps=[[0 0 newDp];dmps];
    eval(['hmD.',srs,'{d}=dmps;']);
    eval(['hmD.meta.',srs,'.dp{d}=newDp;']);
    eval(['dv=datevec(hmD.meta.',srs,'.tim{d}(1));']);
    eval(['hmD.meta.',srs,'.tim{d}=newTim+datenum(dv(1),0,0);']);
end