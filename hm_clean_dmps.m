function hmD=hm_clean_dmps(hmD)
%
% Search for repeating numbers and replace the most frequent one with NAN
% if relative frequency is higher than 15%
%
% hmD=hm_clean_dmps(hmD)
%

%removed
% Clean DMPS spectrum for unrealistic big values
% hmD=hm_clean_dmps(hmD,lim)
% if lim not given the default is 1e7
%
%


%Heikki Junninen
%Jan 2016

% if nargin==1
%     lim=1e7;
% end
nrDays=length(hmD.dmps);

% for d=1:nrDays
%     dat=hmD.dmps{d}(2:end,3:end);
%     Ibad=dat>lim;
%     dat(Ibad)=NaN;
%     hmD.dmps{d}(2:end,3:end)=dat;
% end

% find most frequent numbers

for d=1:nrDays
    totSum=0;
    dat=hmD.dmps{d}(2:end,3:end);
    nrD=size(dat,2);
    %column wise
    [m f]=mode(dat);
    fRel=f/nrD;
    for i=1:length(f)
        if fRel(i)>0.15 & m(i)~=0
            Ibad=dat(:,i)==m(i);
            dat(Ibad,i)=NaN;
            totSum=totSum+sum(Ibad);
        end
    end
    %row wise
    fRel=1;
    while any(fRel>0.15 & m~=0 & ~isnan(m))
        [m f]=mode(dat,2);
        nrD=size(dat,1);
        fRel=f/nrD;
        for i=1:length(f)
            if fRel(i)>0.15 & m(i)~=0
                Ibad=dat(i,:)==m(i);
                dat(i,Ibad)=NaN;
                totSum=totSum+sum(Ibad);
            end
        end
    end
    hmD.dmps{d}(2:end,3:end)=dat;
    if totSum>0
        disp(['hm_clean_dmps: ' ,num2str(totSum), ' repeating values replaced with NaN'])
    end
end