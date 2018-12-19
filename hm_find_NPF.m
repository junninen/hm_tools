function [res,out]=hm_find_NPF(hmD,srs,doPlot)
%
% Find regional new particle formation events
%
%

%check inputs
if ~isstruct(hmD)
    error('hm_find_NPF: input must be structure!')
end

if ~strcmp(hmD.Type,'hm_dat')
    error('hm_find_NPF: input must be structure made by hm_load!')
end

if nargin==1
    srs='dmps';
    doPlot=1;
end

if nargin==2
    doPlot=1;
end

% size bins for EF
dps=10.^(log10(3e-9):((log10(25e-9)-log10(3e-9))/(9-1)):log10(25e-9));
%starting from 7 nm
%dps=10.^(log10(7e-9):((log10(25e-9)-log10(7e-9))/(9-1)):log10(25e-9));


%sun rising hours
tims=eval(['hmD.meta.',srs,'.tim{1}']);
tim=tims(1);
[y,m,d,h]=datevec(tim);
hs=(tims-floor(tims))*24; %decimal hours
lat=51.53;
lon=12.933;
[~,tr,sn,ts]=aurinko(y,m,d,h,lat,lon,round(lon));
tr=6;
ts=18;
datStr=['hmD.',srs];
timStr=['hmD.meta.',srs];
nrDays=length(eval(datStr));

c=hm_conc(hmD,srs,[dps(1),dps(2)]);

c10=NaN(size(c,1),length(dps)-1);
c10(:,1)=c(:,2);
c10_tm=c(:,1);
for si=2:length(dps)-1
    c=hm_conc(hmD,srs,[dps(si),dps(si+1)]);
    c10(:,si)=c(:,2);
end
difSum=sum(c10-repmat(c10(:,end),1,size(c10,2)),2);
nanoSum=sum(c10(:,1:end-2),2);
nanoSum2=sum(c10(:,1:3),2);

Ioutlyer=nanoSum>mean(nanoSum)+std(nanoSum)*4;
nanoSum(Ioutlyer)=NaN;

Idark=hs<=tr | hs>=ts;
Ilight=find(hs>tr & hs<ts);
c_dark=nanmedian(nanoSum(Idark));
c_light=nanmedian(nanoSum(Ilight));

cstd_all=nanstd(nanoSum);

cstd_dark=nanstd(nanoSum(Idark));
cstd_light=nanstd(nanoSum(Ilight));

cmax_dark=nanmax(nanoSum(Idark));
cmax_light=nanmax(nanoSum(Ilight));

c90_dark=prctile(nanoSum(Idark),90);
c90_light=prctile(nanoSum(Ilight),90);

% [hyp,p1]=ttest2(nanoSum(Ilight),(nanoSum(Idark)));
[p1, hyp]=ranksum(nanoSum(Ilight),(nanoSum(Idark)));
% [hyp,p1]=ztest(nanoSum(Ilight),mean(nanoSum(Idark)),std(nanoSum(Idark)))
hyp=hyp+1;
hypStr={'Not different','Different'};

% [hyp2,p2]=ztest(nanoSum(Idark),mean(nanoSum(Ilight)),std(nanoSum(Ilight)))
% hyp2=hyp2+1;
% hypStr={'Not different','Different'};

[mxLight ImxLight]=max(difSum(Ilight));
[mnLight ImnLight]=min(difSum(Ilight));
tm_mxLight=tims(Ilight(ImxLight));
tm_mnLight=tims(Ilight(ImnLight));
mxH=(tm_mxLight-floor(tm_mxLight))*24;
mnH=(tm_mnLight-floor(tm_mnLight))*24;


[apTimes,dpTimes]=H_apperanceTime(c10_tm,nanoSum,2,max(nanoSum)/4);
[apTimes2,dpTimes2]=H_apperanceTime(c10_tm,nanoSum2,2,max(nanoSum2)/4);
apH=(apTimes-floor(apTimes))*24; %decimal hours
apH2=(apTimes2-floor(apTimes2))*24; %decimal hours

dat=hmD.dmps{1}(2:end,3:end);

diffSizeDistribution=[median(dat(Ilight,:))-median(dat(Idark,:))];
I25=hmD.meta.dmps.dp{1}<25e-9;
I10=hmD.meta.dmps.dp{1}<10e-9;
diffSz25=diffSizeDistribution(I25);
%has maximum below 25nm
Imx25=diffSz25(1:end-4)<diffSz25(2:end-3) & diffSz25(2:end-3)<diffSz25(3:end-2) & diffSz25(3:end-2)>diffSz25(4:end-1) & diffSz25(4:end-1)>diffSz25(5:end);
if any(Imx25)
    mx25_dp=hmD.meta.dmps.dp{1}((Imx25));
    mx25_c=diffSz25(find(Imx25));
    mx10_c=max(diffSz25(I10));
else
    mx25_dp=NaN;
    mx25_c=NaN;
    mx10_c=NaN;
end

out=[mxH,mnH,mxLight,mnLight, c_light,c_dark,cstd_light,cstd_dark,cmax_light,cmax_dark,c90_light,c90_dark,cstd_all,mx25_dp,mx25_c,mx10_c];
%% plotting
if doPlot
    figure,
    subplot(4,1,1)
    hm_plot(hmD)
    xl=xlim;
    hold on
    plot(xl,dps(end-1)*[1 1],'k')
    plot(xl,dps(end)*[1 1],'k')
    
    
    subplot(4,1,2)
    plot(c10_tm(:,1),difSum)
    hold all
    yl=ylim;
    xlim(xl);
    plot((floor(tim)+tr/24)*[1 1],yl)
    plot((floor(tim)+ts/24)*[1 1],yl)
    plot(xl,[0 0],'k')
    
    plot(tm_mxLight,difSum(Ilight(ImxLight)),'or')
    plot(tm_mnLight,difSum(Ilight(ImnLight)),'ob')
    
    plot(apTimes,difSum(Ilight(ImnLight)),'ok')
    
    
    H_xdatetick_int(2);
    
    
    subplot(4,1,3)
    plot(c10_tm(:,1),nanoSum)
    H_xdatetick_int(2);
    xlim(xl);
    title(['C_{light}=',num2str(c_light),'C_{dark}=',num2str(c_dark),', ',hypStr{hyp}])
    
    subplot(4,1,4)
    plot(c10_tm(:,1),cumsum(nanoSum))
    H_xdatetick_int(2);
    ylim([-inf inf])
    xlim(xl);
    
end
%% classification

% if mxH < mnH && c_light > c_dark+cstd_dark*2 && (sum(nanoSum>300)>length(nanoSum)*0.2)
if length(unique(nanoSum))<length(nanoSum)*0.85
    resStr=('Bad');
    res=-1;
elseif mxH < mnH && c_light > c_dark+cstd_dark*1 && apH> tr && apH<ts && apH<mxH && (hyp==2 | (hyp==1 && p1>0.05)) && ~(sum(nanoSum<500)>length(nanoSum)*0.95)
    resStr=('nucleation & growth');
    res=1;
elseif mxH > mnH && c_light > c_dark+cstd_dark*1 && apH> tr && apH<ts
    resStr=('nucleation');
    res=0.5;
    % elseif sum(nanoSum<200)>length(nanoSum)*0.95 && isnan(apH) && c_light < c_dark+cstd_dark*1
elseif (hyp==1 || (hyp==2 && p1>0.01)) && isnan(apH) && sum(nanoSum<500)>length(nanoSum)*0.95
    resStr=('NE low conc');
    res=0;
elseif (hyp==1 || (hyp==2 && p1>0.01)) && isnan(apH)
    resStr=('NE high conc');
    res=0.2;
    
elseif hyp==2 && isnan(apH)
    resStr=('UD day different');
    res=0.4;
    
else
    resStr=('UD');
    res=0.3;
end

if doPlot
    subplot(4,1,2)
    title(resStr)
    
    figure,
    subplot(2,2,1)
    semilogx(hmD.meta.dmps.dp{1},[median(dat(Ilight,:));median(dat(Idark,:))]')
    hold all
    semilogx(hmD.meta.dmps.dp{1},diffSizeDistribution')
    yl=ylim;
    plot(10e-9*[1 1],yl,'k')
    plot(25e-9*[1 1],yl,'k')
    legend('Light','Dark','diff','Location','NorthWest')
    xlim([1e-9 1e-6])
    title(hmD.startTime)
    
    subplot(2,2,2)
    I10=hmD.meta.dmps.dp{1}<10e-9;
    I25=hmD.meta.dmps.dp{1}>=10e-9 & hmD.meta.dmps.dp{1}<25e-9;
    
    plot(hmD.meta.dmps.tim{1}(1:end-1),sum(max(diff(dat(:,I10),1),0),2))
    hold on
    plot(hmD.meta.dmps.tim{1}(1:end-1),sum(max(diff(dat(:,I25),1),0),2))
    yl=ylim;
    plot(floor(hmD.meta.dmps.tim{1}(1))+tr/24*[1 1],yl,'k')
    plot(floor(hmD.meta.dmps.tim{1}(1))+ts/24*[1 1],yl,'k')
    datetick('x','HH:MM','keeplimits')
    legend('<10nm','10-25nm')
    title('Cumulative positive difference for each hour')
    xlabel('Time')
    
    subplot(2,2,3)
    semilogx(hmD.meta.dmps.dp{1},sum(abs(max(diff(dat,1),0)),1))
    title('Cumulative positive difference for each size')
    xlabel('Diameter (m)')
    xlim([1e-9 1e-6])
    yl=ylim;
    hold on
    plot(10e-9*[1 1],yl,'k')
    plot(25e-9*[1 1],yl,'k')
    
    subplot(2,2,4)
    semilogx(hmD.meta.dmps.dp{1},sum(abs(max(diff(dat(Ilight,:),1),0)),1))
    title('Cumulative positive difference during light hours for each size')
    xlabel('Diameter (m)')
    xlim([1e-9 1e-6])
    yl=ylim;
    hold on
    plot(10e-9*[1 1],yl,'k')
    plot(25e-9*[1 1],yl,'k')
end

