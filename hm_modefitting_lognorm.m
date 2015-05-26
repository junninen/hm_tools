function hmD=hm_modefitting_lognorm(hmD,inst,mxNrM,dpLim)
%
% Lognormal mode fitting
%
%hmD=hm_modefitting_lognorm(hmD,inst,mxNrM,dpRange)
%
% Uses: MF_lognorm

doSmooth=1;

%check if instrument is loaded
if ~isfield(hmD,inst)
    disp(['hm_modefitting_lognorm: ',inst,' not found!'])
    return
end

if nargin<=3
    dpLim=[0.4e-9,1000e-9];
end
%smooting
if doSmooth
    eval(['dat=H_2dmedfilt(hmD.',inst,'(2:end,3:end),[3,1]);']);
end
% dat=hmD.dmps(2:end,3:end);
eval(['dp=hmD.',inst,'(1,3:end);']);
eval(['tim=hmD.meta.',inst,'.tim{1};']);

%select dp-range
Idp=find(dp>=dpLim(1) & dp<=dpLim(2));
dp=dp(Idp);
dat=dat(:,Idp);

% %if fitting is already done but only for some time period,
% %append fittings
% if isfield(hmD,'fits')
%
%     if isfield(hmD.fits,inst)
%         eval(['[n,m]=size(hmD.fits.',inst,'.zs);']);
%
%         eval(['zs=hmD.fits.',inst,'.zs;']);
%         eval(['ws=hmD.fits.',inst,'.ws;']);
%         eval(['Ns=hmD.fits.',inst,'.Ns;']);
%     else
%         zs=NaN(size(tim,1),mxNrM);
%         ws=NaN(size(tim,1),mxNrM);
%         Ns=NaN(size(tim,1),mxNrM);
%     end
% else
zs=NaN(size(tim,1),mxNrM);
ws=NaN(size(tim,1),mxNrM);
Ns=NaN(size(tim,1),mxNrM);
% end

%if only part of data is fitted
% [mn,I1]=min(abs(tim-x(1)));
% [mn,I2]=min(abs(tim-x(2)));

I1=1;
I2=length(tim);
for i=I1:I2
    [param,fval,yhat,N]=MF_lognorm(dat(i,:),dp,0,mxNrM);
    %     drawnow
    nrM=length(N);
    zs(i,1:nrM)=param(nrM+1:end);
    ws(i,1:nrM)=param(1:nrM);
    Ns(i,1:nrM)=N;
end

Idel=(zs<=dpLim(1) | zs>=dpLim(2));
if any(Idel)
    zs(Idel)=NaN;
    ws(Idel)=NaN;
    Ns(Idel)=NaN;
end

if isfield(hmD,'fits')
    if isfield(hmD.fits,inst)
        eval(['hmD.fits.',inst,'.tm=tim;']);
        eval(['hmD.fits.',inst,'.zs=[hmD.fits.',inst,'.zs,zs];']);
        eval(['hmD.fits.',inst,'.ws=[hmD.fits.',inst,'.ws,ws];']);
        eval(['hmD.fits.',inst,'.Ns=[hmD.fits.',inst,'.Ns,Ns];']);
    else
        eval(['hmD.fits.',inst,'.tm=tim;']);
        eval(['hmD.fits.',inst,'.zs=zs;']);
        eval(['hmD.fits.',inst,'.ws=ws;']);
        eval(['hmD.fits.',inst,'.Ns=Ns;']);
    end
else
    eval(['hmD.fits.',inst,'.tm=tim;']);
    eval(['hmD.fits.',inst,'.zs=zs;']);
    eval(['hmD.fits.',inst,'.ws=ws;']);
    eval(['hmD.fits.',inst,'.Ns=Ns;']);
end