function hm_mf_write_ascii(hmD,srs,savePath,hdr,fnm)
%save fitting results as ascii file
%
% hm_mf_write_ascii(hmD,srs,savePath,[hdr],[fnm])
%
% [] - inputs are optional
% hmD - hm-structure with mode fitting results
% srs - name of instrument, e.g. dmps, aisn
% savePath - full path where to save
% hdr - 0 - no header written
%       1 - write variable names (default)
%       2 - convert variable names to numbers: double(variableNameRow)
%           when reloaded to MatLab this can be converted to
%           variable names by char(firstRow)
% fnm - file name model (optional),
%       if not given file name will be hmMF_yyyymmdd.txt
%       for the name agove fnm='sprintf(''hmMF_%04d%02d%02d.txt'',yyyy,mm,dd)'
% file structure will be:
% y m d h m s z1 z2 zn w1 w2 wn H1 H2 Hn N1 N2 Nn
%
% n-varies from day to day
%
%


% fls=dir('C:\1h\da\mf\da2\mat');
%
% for f=1:length(fls)


%     datestr(d,30)
% end

%% check variables
if nargin<3
    disp('hm_mf_write_ascii: Not enough input parameters. Nothing done!')
    return
elseif nargin==3
    hdr=1;
    fnm='sprintf(''hmMF%04d%02d%02d.txt'',yyyy,mm,dd)';
elseif nargin==4
    if ischar(hdr)
        fnm=hdr;
        hdr=1;
    else
        fnm='sprintf(''hmMF_%04d%02d%02d.txt'',yyyy,mm,dd)';
    end
end

if hdr>2
    disp('hm_mf_write_ascii: wrong hdr, it must be 0, 1 or 2. Using hdr=1!')
end

if ~isfield(hmD,'fit')
    disp('hm_mf_write_ascii: field ''fit'' is not found. Nothing done!')
    return
end



%% reminder
disp('hm_mf_write_ascii: REMINDER only one day structures implemented')

%% action

dvec=datevec(hmD.meta.startTime);
yyyy=dvec(1);
mm=dvec(2);
dd=dvec(3);

eval(['fnam=',fnm]);
[fid,msg]=fopen([savePath,fnam],'w');

if fid==2
    disp(msg)
end

eval(['fits=hmD.fit.',srs,';']);

eval(['tim=hmD.meta.',srs,'.tim{1};']);
dvec=datevec(tim);
zs=fits.zs{1};
ws=fits.ws{1};
Hs=fits.Hs{1};
Ns=fits.Ns{1};
nrM=size(zs,2);

if hdr==1
    timH=['year\tmonth\tday\thour\tmin\tsec\t'];
    zsH='GMD1\t';
    wsH='sigma1\t';
    HsH='area1\t';
    NsH='numberConc1\t';
    for m=2:nrM-1
        zsH=[zsH,'GMD',num2str(m),'\t'];
        wsH=[wsH,'sigma',num2str(m),'\t'];
        HsH=[HsH,'area',num2str(m),'\t'];
        NsH=[NsH,'numberConc',num2str(m),'\t'];
    end
    zsH=[zsH,'GMD',num2str(m+1),'\t'];
    wsH=[wsH,'sigma',num2str(m+1),'\t'];
    HsH=[HsH,'area',num2str(m+1),'\t'];
    NsH=[NsH,'numberConc',num2str(m+1)];

    fprintf(fid,[timH,zsH,wsH,HsH,NsH,'\r\n'])
elseif hdr==2
    disp('hm_mf_write_ascii: REMINDER. hdr=2 not implemented, Using hdr=0')
end

fprintf(fid,['%4d\t%02d\t%02d\t%02d\t%02d\t%2.2f\t',repmat('%1.8E\t',1,nrM),repmat('%1.3f\t',1,nrM),repmat('%5.4f\t',1,nrM),repmat('%5.4f\t',1,nrM-1),'%5.4f\r\n'],[dvec,zs,ws,Hs,Ns]');

fclose(fid);
