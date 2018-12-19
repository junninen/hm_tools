function AEdat=hm_comb_days(AEdat)
%
% combine days to one matrix and interpolate to same grid if .
%
%hmD=hm_comb_days(hmD)
%
%hmD - data struct, made by hm_laod (type: hm_dat)
%
%

%heikki Junninen
%12.03.2007
% Update 14.03.2007
%log
%29.May.2008
% data saved in cell, for compatability with other functions
% splits big data sets to many blocks, in order to fit in memory

%check inputs
if ~isstruct(AEdat)
    error('hm_comb_days: input must be structure!')
end

if ~strcmp(AEdat.Type,'hm_dat')
    error('hm_comb_days: input must be structure made by hm_load!')
end



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
nrInstr=length(instruments);
processInstrument=zeros(nrInstr,1); %looks for and loads all instrumnets
%if remains all zeros after input evaluation, will be replaced to ones.
%1=loads instrument data 0=does not load

for i=1:nrInstr
    if isfield(AEdat,instruments{i})
        processInstrument(i)=1;
    end
end

%particle size bins for new grid
% startDp=log10(3e-9)
% endDp=log10(990e-9)
% nrBins=40
% int=(endDp-startDp)/(nrBins-1);
%
% newDp=startDp:int:endDp;
%
% timeVector=AEdat.meta.startTime:AEdat.meta.endTime;
% timeVector=AEdat.meta.dmps.tim{1}(1):AEdat.meta.dmps.tim{end}(end);
for i=1:nrInstr
    if processInstrument(i)
        %         dpBig=[];
        %         timBig=[];
        %         tiBig=[];
        %         datBig=[];
        
        %check how many days are loaded
        if ~iscell(eval(['AEdat.',instruments{i}]))
            disp('hm_comb_days: Days are already combined or there is only one day!')
            disp('no changes made');
            return
        end
        
        datStr=['AEdat.',instruments{i}];
        timStr=['AEdat.meta.',instruments{i}];
        nrDays=length(eval(datStr));
        eval(['Iempty=cellfun(@(x) isempty(x),',datStr,');'])
        if any(Iempty)
            disp('hm_comb_days: Some days are missing. Can''t combine!')
            return
        end
        sumf=cell(nrDays,1);
        dp=sumf;        
        temp=sumf;
        for id=1:nrDays
            sumf{id}=eval([datStr,'{',num2str(id),'}']);
            %             [n,m]=size(sumf{id});
%             dp{id}=sumf{id}(1,3:end);
            dp{id}=eval([timStr,'.dp','{',num2str(id),'}']);
            %             nrBins(id)=length(dp{id});
            %             ti{id}=sumf{id}(2:n,1); %julian date
            %             tim{id}=timeVector(id)+(ti{id}-floor(ti{id})); %matlab date
        end
        
        %if data has same number of bins
        %collect one big matrix for dp,for tim, and for data
        %make one interpolation instead of one for each day
        %         if length(unique(nrBins))==1
        %             for id=1:nrDays
        %                 [n,m]=size(sumf{id});
        %                 dpBig=[dpBig;repmat(dp{id},n-1,1)];
        %                 timBig=[timBig;repmat(tim{id},1,m-2)];
        %                 tiBig=[tiBig;repmat(ti{id},1,m-2)];
        %                 datBig=[datBig;sumf{id}(2:n,3:m)];
        %             end
        %         end
        %         disp('interpolating for new grid ...')
        
        %time interval will be median of current intervals
        %for matlab time
        %         medInterval=median(timBig(2:end,1)-timBig(1:end-1,1));
        %
        %         newTim=[timBig(1):medInterval:timBig(end,1)]';
        %         dati = griddata(timBig,dpBig,datBig,...
        %             repmat(newTim,1,length(newDp)),repmat(newDp,size(newTim,1),1),'linear');
        
        %         for julian time
        %         disp('think')
        %         medInterval=median(tiBig(2:end,1)-tiBig(1:end-1,1));
        
        %         newTi=[tiBig(1):medInterval:tiBig(end,1)]';
        %         newTim=[timBig(1):medInterval:timBig(end,1)]';
        
        %         %split to blocks if too big matrix
        %         [n,m]=size(newTi);
        %         if n>3000
        %             %split
        %             lim=0:3000:n;
        %             if lim(end)~=n,lim=[lim,n];end
        %             for si=1:length(lim)-1
        %                 dati_tmp{si} = griddata(tiBig,log10(dpBig),datBig,...
        %                     repmat(newTi(lim(si)+1:lim(si+1),:),1,length(newDp)),repmat(newDp,size(newTi(lim(si)+1:lim(si+1),:),1),1),'linear');
        %             end
        %             dati=cat(1,dati_tmp{:});
        %         else
        %             %make all in one go
        %             dati = griddata(tiBig,log10(dpBig),datBig,...
        %                 repmat(newTi,1,length(newDp)),repmat(newDp,size(newTi,1),1),'linear');
        
        %pick the longest as common
        L=cellfun(@(x) length(x),dp);
        [~,Imx]=max(L);
        commonDp=dp{Imx};
        
        %         temp{1}=sumf{1}(2:end,3:end);
        for id=1:nrDays
            rawDat=sumf{id}(2:end,3:end);
            rawDp=dp{id};
            temp{id}=interp1(rawDp,rawDat', commonDp)';
        end
        dati=cat(1,temp{:});
        clear temp
        eval(['newTi=cat(1,AEdat.meta.',instruments{i},'.tim{:});']);
        %         end
        
        %put nans back if gap longer than hour
        for id=1:nrDays-1
            t1=sumf{id}(end,1);
            t2=sumf{id+1}(2,1);
            
            if (t2-t1)*24>1
                In=newTi>t1 & newTi<t2;
                dati(In,:)=NaN;
            end
        end
        
        %replace original data with interpolated one
        %         eval([datStr,'{',num2str(id),'}=dati;']);
        dati=[newTi,NaN(size(newTi)),dati];
        dati=[[NaN,NaN,commonDp];dati];
        eval([datStr,'={dati};']);
        eval([timStr,'.tim={newTi};']);
        %         eval([timStr,'.dp={10.^newDp};']);
        eval([timStr,'.dp={commonDp};']);
        
    end
end

AEdat=hm_history(AEdat,mfilename);