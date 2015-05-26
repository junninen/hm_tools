function hmD=hm_gr_auto(hmD,sizeLim,concLim,timLim)
%
% Read in hmD structure with mode fitting done
% define growth rates
%
% hmD=hm_gr_auto(hmD,[sizeLim],[concLim])

%Heikki Junninen
if nargin==1
    sizeLim=20*1e-9;
    concLim=1000;
    timLim=[0 24];
elseif nargin==2 | nargin==3
    error('hm_gr_auto: wrong number of inputs. Can be 1 or 4')
end

disp('currently only for DMSP')
if ~isfield(hmD,'fit')
    %make mode fitting if not done yet
    disp('hm_gr_auto: makes mode fitting to DMPS data')
    if hmD.meta.dmps.smoothing
        hmD=hm_smoothing(hmD);
    end
    hmD=hm_mf(hmD,'dmps','lognorm',3,0);
end
tm=hmD.meta.dmps.tim{1};
dp=hmD.meta.dmps.dp{1};

dv=datevec(tm);
nrM=size(hmD.fit.dmps.zs{1},2);
hh=repmat(dv(:,4),1,nrM);
zs=hmD.fit.dmps.zs{1};
zs(zs<1e-9)=NaN;
I=(zs<sizeLim & hmD.fit.dmps.Ns{1}>concLim & hh>timLim(1) & hh<timLim(2));
Inucl=find(any(I,2));


if ~isempty(Inucl)
    [mx,Imx]=max(hmD.fit.dmps.Ns{1}.*I,[],2);

    I=hmD.fit.dmps.Ns{1}==repmat(mx,1,nrM);
    zok=hmD.fit.dmps.zs{1}.*I;
    zok(isnan(zok))=0;
    Inucl_first=Inucl(1);
    Inucl_last=Inucl(end);

    sizeFirst=hmD.fit.dmps.zs{1}(Inucl_first,:);
    sizeLast=hmD.fit.dmps.zs{1}(Inucl_last,:);

    [mnSizeFirst]=min(sizeFirst);

    okLast=sizeLast>mnSizeFirst;

    if any(okLast)

        %     [mx,Imx]=max(hmD.fit.dmps.Ns{1}(Inucl_last,:).*I(Inucl_last,:));

        %     zs_sel1=hmD.fit.dmps.zs{1}(Inucl_first,:);
        zs_sel1=mnSizeFirst;
        
        zs_sel2=min(hmD.fit.dmps.zs{1}(Inucl_last,okLast));

%         zs_sel1=(hmD.fit.dmps.zs{1}(Inucl_first,:).*I(Inucl_first,:));
%         zs_sel2=(hmD.fit.dmps.zs{1}(Inucl_last,:).*I(Inucl_last,:));
% 
%         zs_sel1(zs_sel1==0)=1;
%         zs_sel2(zs_sel2==0)=1;
        %
        a=[tm([Inucl_first,Inucl_last]),[min(zs_sel1);min(zs_sel2)]];

        [hmD,gr_str]=hm_gr(hmD,'dmps',a);

        figure(111)
        hm_plot(hmD)
        hold on
        plot(hmD.meta.dmps.tim{1},hmD.fit.dmps.zs{1},'.k')
        plot(gr_str.tm,gr_str.z,'w.')
        x=(gr_str.tm-floor(gr_str.tm(1)))*24;
        plot(gr_str.tm,x*gr_str.gr_avr(1)+gr_str.int_avr(1),'w')
        % find datapoints falling in between 1-5 nm/h
        hold off
        hmD.fit.dmps.gr_str=gr_str;
    end
end

