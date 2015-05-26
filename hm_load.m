function AEdat=hm_load(dates,varargin)
%
% hmD=hm_load(srtTim,[[argID,] value, ...])
% hmD=hm_load(srtEndTim,[[argID,] value, ...])
%
%INPUTS ([]'s are optional):
% srtTim    (string) 'yyyy-mm-dd HH:MM:SS' or 'yyyy-mm-dd' or 'yyyymmdd' or
%           (double) datenum
% srtEndTim    (string) 'yyyy-mm-dd HH:MM:SS' or 'yyyy-mm-dd' or 'yyyymmdd' or
%           (double) datenum
%                    if endTim empty the current day will be loaded
% [argID,   (string) See below. The values which are unambiguous can
%   value]  (varies) be given without the preceeding argID.
%
% Argument IDs and corresponding values. The values
% which are unambiguous (marked with '*') can be given without the
% preceeding argID.
%
% parPath   (string)  full path to parameters file
%                      parameters file contains all pathes for files
%                      and file-name-model
% dmps       *(string)  load dmps data
% ais        *(string)  load ais negative and positive data
% aisn       *(string)  load ais negative data
% aisp       *(string)  load ais positive data
% nais       *(string)  load nais negative and positive data
% naisn      *(string)  load nais negative data
% naisp      *(string)  load nais positive data
% ais_block  *(string)  load ion data from block-files
% nais_block *(string)  load particles data from block-files
% bsma       *(string)  load bsma negative and positive data
% bsman      *(string)  load bsma negative data
% bsmap      *(string)  load bsma positive data
% vol        *(string)  load volatility hot and cold data
% volh       *(string)  load volatility hot data
% volc       *(string)  load volatility cold data
%
%OUTPUT:
%
% hmD       (struct) .type %structure type
%                    .StartTimeString
%                    .EndTimeString
%                    .meta  %structure for meta information
%                         .startTime %time in datenum format
%                         .endTime   %time in datenum format
%                         .history %history of data manipulation functions
%                         .insts %list of instruments, hard coded
%                         .insts_loaded %loaded instruments = 1
%                         .insts % list of all instruments, used by other
%                                  HMTools functions
%                               .paths %data path of an instrument
%                               .nfo   %info
%                               .dp    %Millikan diameter, for DMPS, AIS,
%                                      %NAIS
%                               .tim   %time in datenum-format
%                               .raw   %raw, not modified data 
%                               .smoothing %smoothing done = 1 (if more than one time, it will be a vector) 
%                               .win   %smoothing windowses
%                               .smoothMeth % %smoothing methods
%                    .dmps  % the dmps raw data (= inverted data)
%                    .aisn  % the ais negative raw data
%                    .aisp  % the ais positive raw data
%                    .cpc   % the cpc raw data
%                    .naisn % the nais negative raw data
%                    .naisp % the nais positive raw data
%                    .bsman % the BSMA neg data
%                    .bsmap % the BSMA pos dta
%                    .volh  % volatility data hot
%                    .volc  % volatility data cold
%                    .htdma % HTDMA data (NOT IMPLEMENTED)
%                    .metData  %structure of selected meteorological data
%
%EXAMPLES:
%
%Heikki Junninen 20.02.07
%update: 18.03.2007

%logBook
%29.05.2007
%  AIS and NAIS meta.dp is millikan
%  hm_load will not load the external file, it will be loaded here
%  and written to meta.dp field
%15.05.2007
%  add index of loaded instruments to meta
%29.05.2007
%  add BSMA dp calculations
%  add dp to meta data
%01.06.2007
%  convert all datas time to datenum and add it to meta
%5.12.2007
%  load all datas automaticallly to cell, not like before, to matric if one
%  day and to cell if many days
%
%Acknowledgement
%Lots of coding tricks are copied from the excelent SOM toolbox.
% http://www.cis.hut.fi/projects/somtoolbox/

% versn=0.01;

%defaults
%if remains all zeros after input evaluation, will be replaced to ones.
%1=loads instrument data 0=does not load
instruments={...
    'dmps',...
    'aisn',...
    'aisp',...
    'cpc',...
    'naisn',...
    'naisp',...
    'bsman',...
    'bsmap',...
    'volh',...
    'volc',...
    'aisn_block',...
    'aisp_block',...
    'naisn_block',...
    'naisp_block'...
    };
nrI=length(instruments);
loadInstrument=zeros(nrI,1); %looks for and loads all instrumnets

global PARAPATH;
%  global PARPATH; PARAPATH='C:\1h\dat\parameters\';
parPath=[PARAPATH,'hm_param.txt'];

%extract inputs

if length(dates)==1
    srtTim=dates;
    endTim=dates;
elseif length(dates)==2
    srtTim=dates(1);
    endTim=dates(2);
end

if ~iscell(varargin)
    endTim=varargin;
else
    i=1;
    while i<=length(varargin),
        argok = 1;
        if ischar(varargin{i}),
            switch varargin{i},
                % argument IDs
                case 'parPath',       i=i+1; parPath = varargin{i};

                    % unambiguous values
                case {'dmps'}, loadInstrument(1)=1;
                case {'aisn'}, loadInstrument(2)=1;
                case {'aisp'}, loadInstrument(3)=1;
                case {'ais'}, loadInstrument(2)=1;loadInstrument(3)=1;
                case {'cpc'},  loadInstrument(4)=1;
                case {'nais'}, loadInstrument(5)=1;loadInstrument(6)=1;
                case {'naisn'}, loadInstrument(5)=1;
                case {'naisp'}, loadInstrument(6)=1;
                case {'bsma'}, loadInstrument(7)=1;loadInstrument(8)=1;
                case {'bsman'}, loadInstrument(7)=1;
                case {'bsmap'}, loadInstrument(8)=1;
                case {'vol'}, loadInstrument(9)=1;loadInstrument(10)=1;
                case {'volh'}, loadInstrument(9)=1;
                case {'volc'}, loadInstrument(10)=1;
                case {'aisn_block'}, loadInstrument(11)=1;
                case {'aisp_block'}, loadInstrument(12)=1;
                case {'naisn_block'}, loadInstrument(13)=1;
                case {'naisp_block'}, loadInstrument(14)=1;
                case {'ais_block'}, loadInstrument(11)=1;loadInstrument(12)=1;
                case {'nais_block'}, loadInstrument(13)=1;loadInstrument(14)=1;
                    
                    %case {'htdma'}, loadInstrument(11)=1;
                otherwise, argok=0;
            end
        else
            argok = 0;
        end
        if ~argok,
            disp(['hm_load: Ignoring invalid argument #' num2str(i+1)]);
        end
        i = i+1;
    end
    if ~any(loadInstrument)
        loadInstrument=loadInstrument+1;
    end
end

%extract date
if ischar(srtTim)
    L=length(srtTim);
    if L==19
        tim{1}=datenum(srtTim, 'yyyy-mm-dd HH:MM:SS');
    elseif L==10
        tim{1}=round(datenum(srtTim, 'yyyy-mm-dd'));
    elseif L==8 && isempty(strfind('srtTim','-'))
        tim{1}=round(datenum(srtTim, 'yyyymmdd'));
    else
        error('hm_load: Strange date form. Should be ''yyyy-mm-dd HH:MM:SS'' or ''yyyy-mm-dd'' or ''yyyymmdd'' or datenum')
    end
else
    tim{1}=srtTim;
end

if ischar(endTim)
    L=length(endTim);
    if L==19
        tim{2}=datenum(endTim, 'yyyy-mm-dd HH:MM:SS');
    elseif L==10
        tim{2}=round(datenum(endTim, 'yyyy-mm-dd'));
    elseif L==8 && isempty(strfind('endTim','-'))
        tim{2}=round(datenum(endTim, 'yyyymmdd'));
    else
        error('hm_load: Strange date form! Should be ''yyyy-mm-dd HH:MM:SS'' or ''yyyy-mm-dd'' or ''yyyymmdd'' or datenum')
    end
else
    tim{2}=endTim;
end

if isempty(tim{2}),tim{2}=tim{1};end

%find out how many files needed to load
tims=tim{1}:tim{2};
nrDays=length(tims);
yyyyD=NaN(nrDays,1);
yyD=NaN(nrDays,1);
mmD=NaN(nrDays,1);
ddD=NaN(nrDays,1);
for i=1:nrDays,
    ti=tims(i);
    dv=datevec(ti);
    yyyyD(i)=dv(1);
    mmD(i)=dv(2);
    ddD(i)=dv(3);
    
    if yyyyD(i)>=2000
        yyD(i)=yyyyD(i)-2000;
    else
        yyD(i)=yyyyD(i)-1900;
    end
end

%load parameters from parameters file
params=hm_readInstrumentParam(parPath);

%collect parameters
[names{1:length(params)}]=deal(params.name);
names=lower(names);
inOk=ones(size(loadInstrument));
for i=1:length(loadInstrument)
    if loadInstrument(i)
        Iin=strncmp(instruments{i},names,length(instruments{i}));
        if sum(Iin)==1
            inPath{i}=params(Iin).path;
            fnm{i}=params(Iin).fnm;
            if isfield(params(Iin),'alt_fnm')
                alt_fnm{i}=params(Iin).alt_fnm;
            else
                alt_fnm{i}=[];
            end
            inOk(i)=1;
        elseif sum(Iin)==2
            error(['hm_load: ',instruments{i},' defined twice in parameter file'] )
        else
            inPath{i}=[];
            fnm{i}=[];
            alt_fnm{i}=[];
            inOk(i)=0;
        end
    end
end

% %check if instruments that have been requested to load have parameters
% if any(inOk==0),
%     disp('hm_load: Some of the requested instruments dont have parameters in')
%         disp(parPath)
%     I=find(~inOk);
%     for i=1:length(I)
%         disp(['Not loaded: ', instruments{I(i)}]);
%     end
% end
%
%

%load instruments

nrFailed=zeros(length(loadInstrument),1);
for i=1:length(loadInstrument),
    if inOk(i)&&loadInstrument(i)
        for d=1:nrDays
            yyyy=yyyyD(d);yy=yyD(d);mm=mmD(d);dd=ddD(d);
            
            inFile=eval(fnm{i});
            inFile2=' ';
            if ~isempty(alt_fnm{i}),
                inFile2=eval(alt_fnm{i});
            end
            if ~exist([inPath{i},inFile],'file') && exist([inPath{i},inFile2],'file')==2
                inFile=inFile2; %use alternative fnm
            end
            
            day=datenum(yyyy,mm,dd);
            try                
                %for ais after July 2006 need special load function
                %only if read nd files, the sum files are normal ones
                if strcmp('ais',instruments{i}) & day>datenum(2006,07,31) & ~strcmp('sum',inFile(end-2:end))
%                     if ~exist('posSum','var') %check if already loaded
                        [posSum,negSum,tim_ais,dp_ais,mob_ai,posDat,posErr,negDat,negErr]=hm_load_ais_nd([inPath{i},inFile]);
                    intim{i}{d}=tim_ais;
%                     end
                    if strmatch('aisn',instruments{i})
                        datSum=negSum;
                        intim{i}{d}=tim_ais;
                        dp{i}{d}=dp_ais;
                        mob_ais{i}{d}=mob_ai;
%                         clear negSum
                    else
                        datSum=posSum;
                        intim{i}{d}=tim_ais;
                        dp{i}{d}=dp_ais;
                        mob_ais{i}{d}=mob_ai;
%                         clear posSum
                    end
                elseif strmatch('bsma',instruments{i})
                    [bsmaDat]=hm_load_bsma([inPath{i},inFile],day);
                    if strmatch('bsman',instruments{i})
                        datSum=bsmaDat.datNeg;
                        datMobNeg{i}{d}=bsmaDat.mobNeg;
                        datMobNegAbs{i}{d}=bsmaDat.mobNegAbs;
                        
                    else
                        datSum=bsmaDat.datPos;
                        datMobPos{i}{d}=bsmaDat.mobPos;
                        datMobPosAbs{i}{d}=bsmaDat.mobPosAbs;
                    end
                    dp{i}=datSum(1,3:end);
                    intim{i}{d}=datSum(2:end,1);
                elseif strfind(instruments{i},'block')
                    [neg,pos,dps,mob,blockStartTim,blockStopTim]=hm_load_ais_block([inPath{i},inFile]);
                    
                     if strcmp('aisn_block',instruments{i})
                        
                         neg=[blockStartTim,zeros(size(neg,1),1),neg];
                         neg=[0,0,dps;neg];
                         
                         datSum=neg;
                        intim{i}{d}=blockStartTim;
                        dp{i}{d}=dps;
                        mob_ais{i}{d}=mob;
%                         clear negSum
                     elseif strcmp('aisp_block',instruments{i})
                         pos=[blockStartTim,zeros(size(pos,1),1),pos];
                         pos=[0,0,dps;pos];
                        datSum=pos;
                        intim{i}{d}=blockStartTim;
                        dp{i}{d}=dps;
                        mob_ais{i}{d}=mob;
%                         clear posSum
                     elseif strcmp('naisp_block',instruments{i})
                         pos=[blockStartTim,zeros(size(pos,1),1),pos];
                         pos=[0,0,dps;pos];
                        datSum=pos;
                        intim{i}{d}=blockStartTim;
                        dp{i}{d}=dps;
                        mob_ais{i}{d}=mob;
%                         clear posSum
                     elseif strcmp('naisn_block',instruments{i})
                         neg=[blockStartTim,zeros(size(neg,1),1),neg];
                         neg=[0,0,dps;neg];
                        datSum=neg;
                        intim{i}{d}=blockStartTim;
                        dp{i}{d}=dps;
                        mob_ais{i}{d}=mob;
%                         clear negSum
                     end
                    
                else %normal sum files of other instruments
                    datSum=load([inPath{i},inFile]);
                    intim{i}{d}=datSum(2:end,1);
                    dp{i}{d}=datSum(1,3:end);
                    if all(dp{i}{d}>100e-9)
                        %convert to meters
                        dp{i}{d}=dp{i}{d}*1e-9;
                    end
                    %if time is in julian day convert it to datenum
                    if intim{i}{d}(1)<400
                        intim{i}{d}=intim{i}{d}+datenum(yyyy,1,1)-1;
                    elseif intim{i}{d}(1)>1e10
                        yy=floor(intim{i}{d}/10000000000);
                        mm=floor(intim{i}{d}/100000000)-yy*100;
                        dd=floor(intim{i}{d}/1000000)-(yy*10000+mm*100);
                        HH=floor(intim{i}{d}/10000)-(yy*1000000+mm*10000+dd*100);
                        MM=floor(intim{i}{d}/100)-(yy*100000000+mm*1000000+dd*10000+HH*100);
                        SS=floor(intim{i}{d}/1)-(yy*10000000000+mm*100000000+dd*1000000+HH*10000+MM*100);
                        intim{i}{d}=datenum(yy,mm,dd,HH,MM,SS);
                    end
                end

                %load external file with diameter conversions
                    %replace dp and write this dp to meta-field if file
                    %found
                    %this applies only for ais and bsma

                if strcmp('aisn',instruments{i}) || strcmp('aisp',instruments{i})                 
                    convFile='diameters_AIS_mobility_tammet_milligan.dat';
                    if exist(convFile,'file')
                        %     first col mobility
                        %     second Tammet
                        %     third Millikan diameter
                        convD=load(convFile);
                        if exist('dp_ais','var')
                            dp{i}{d}=convD(1:length(dp_ais),3)';%Millikan
                        
                        end
                    end
                elseif strcmp('naisn',instruments{i}) || strcmp('naisp',instruments{i})                 
                    convFile='diameters_NAIS_mobility_tammet_milligan.dat';
                    if exist(convFile,'file')
                        %     first col mobility
                        %     second Tammet
                        %     third Millikan diameter

                        convD=load(convFile);
                        if exist('dp_ais','var')
                            dp{i}{d}=convD(1:length(dp_ais),3)';%Millikan
                        end
                    end
%                 elseif strcmp('bsman',instruments{i}) || strcmp('bsmap',instruments{i})                 
%                     convFile='diameters_BSMA_mobility_tammet_milligan.dat';
%                     if exist(convFile,'file')
%                         %     first col mobility
%                         %     second Tammet
%                         %     third Millikan diameter
%                         convD=load(convFile);
%                         dp{i}{d}=convD(:,3)';%Millikan
%                     end
                end

            catch ER
                disp(['hm_load: ',ER.message]);
                disp(['hm_load: Error when reading ',[inPath{i},inFile],' ']);
                datSum=[];
                nrFailed(i)=nrFailed(i)+1;
                intim{i}{d}=[];
            end
            datSums{i,d}=datSum;
            
        end
    end
end

inOk(find(nrFailed==nrDays))=0;


meta=struct('startTime',tim{1},'endTime',tim{2},'history',[]);
histry={datestr(now),mfilename};
meta.history=histry;
AEdat=struct('Type','hm_dat','startTime',datestr(tim{1}),'endTime',datestr(tim{2}),'meta',meta);

eval(['AEdat.meta.insts=instruments;']);
eval(['AEdat.meta.insts_loaded=logical(inOk);']);
for i=1:length(loadInstrument)
    if inOk(i)&&loadInstrument(i)
%         if nrDays==1
%             eval(['AEdat.',instruments{i},'=datSums{i};']);
%         else
            eval(['AEdat.',instruments{i},'=datSums(i,:);']);
%         end
        eval(['AEdat.meta.',instruments{i},'.paths=inPath{i};']);
        eval(['AEdat.meta.',instruments{i},'.nfo=[];']);
        eval(['AEdat.meta.',instruments{i},'.dp=dp{i};']);   
        eval(['AEdat.meta.',instruments{i},'.tim=intim{i};']);
        eval(['AEdat.meta.',instruments{i},'.raw=[];']);
        eval(['AEdat.meta.',instruments{i},'.smoothing=0;']);
        if strmatch('bsma',instruments{i})
            %if bsma, add also the mobility matricies
            if strmatch('bsman',instruments{i})
                eval(['AEdat.meta.',instruments{i},'.datMobNeg=datMobNeg{i};']);
                eval(['AEdat.meta.',instruments{i},'.datMobNegAbs=datMobNegAbs{i};']);
                eval(['AEdat.meta.',instruments{i},'.mob_limits=[bsmaDat.mobNeg_low;bsmaDat.mobNeg_hgh];']);
            else
                eval(['AEdat.meta.',instruments{i},'.datMobPos=datMobPos{i};']);
                eval(['AEdat.meta.',instruments{i},'.datMobPosAbs=datMobPosAbs{i};']);
                eval(['AEdat.meta.',instruments{i},'.mob_limits=[bsmaDat.mobPos_low;bsmaDat.mobPos_hgh];']);
            end
        end
        if strmatch('aisn',instruments{i},'exact') & exist('negErr')
            %if ais, and time is after july 2006, save also errors
            %and mobilities and size bins         
            eval(['AEdat.meta.',instruments{i},'.nfo=''mob unit (cm^2/V/s)'';']);
%             eval(['AEdat.meta.',instruments{i},'.dp=dp_ais;']);    
            eval(['AEdat.meta.',instruments{i},'.mob=mob_ais{i};']);
%             eval(['AEdat.meta.',instruments{i},'.tim=tim_ais;']);
            eval(['AEdat.meta.',instruments{i},'.err=negErr;']);          
            
        end
        if strmatch('aisp',instruments{i},'exact') & exist('posErr')
            %if ais, and time is after july 2006, save also errors
            %and mobilities and size bins
            eval(['AEdat.meta.',instruments{i},'.nfo=''mob unit (cm^2/V/s)'';']);
%             eval(['AEdat.meta.',instruments{i},'.dp=dp_ais;']);    
            eval(['AEdat.meta.',instruments{i},'.mob=mob_ais{i};']); 
%             eval(['AEdat.meta.',instruments{i},'.tim=tim_ais;']);
            eval(['AEdat.meta.',instruments{i},'.err=posErr;']);          
        end
    end
end

function dp=calc_bsma_dp
%calculate dp for bsma
%values from BSMA2PLT.M
%

dpY=[0.5,0.7,1.0,1.5,2.5,4.0,6.5]*1e-9;
dpX=[1 2.18 3.43 4.85 6.65 8.3 10];
bin=1:10;
dp=interp1(dpX,dpY,bin);

% if doComb
%     [AEdat,chProb]=AE_interp(AEdat,0);
% end