function hm_plot(AEstr,varargin)
%
% hm_plot(hmD,plotType)
% Plot pcolor plots
%
% if all 3 data present 2,2 subplots
% if only ais data 2,1 subplot or input 'ais'
% if only dmps 1 plot ot input 'dmps'
%
%plotTypes:
%     'all'
%     'ais',
%     'dmps',
%     'aisn',
%     'aisp',
%     'naisn',
%     'naisp',
%     'nais4n',
%     'nais4p',
%     'nais',
%     'ions',
%     'comb',
%     'diff',
%     'vol',
%     'volc',
%     'volh',
%     'meas',
%     'fits'


caxisMax=4;

% if AEstr.aisn==0, plotaisn=0;end
% if AEstr.aisp==0, plotaisp=0;end
% if AEstr.dmps==0, plotdmps=0;end

if nargin==1
    plotType='all_autom';
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
    
    h=0;
    plotInstr=[];
    for i=1:length(instruments),
        if isfield(AEstr,instruments{i})
            h=h+1;
            plotInstr{h}=instruments{i};
        end
    end
    if isempty(plotInstr)
        disp('hm_plot: no data to plot')
        return
    end
else
    plotType=varargin{1};
end
plotFits=0;
if length(varargin)==2
    if strcmp(varargin{2},'fits')
        plotFits=1;
    end
end
switch plotType
    case 'all_autom'
        lowDP=0.3e-9;
        hghDP=1100e-9;
        nrPlots=length(plotInstr);
        nsub=ceil(sqrt(nrPlots));
        msub=round(nrPlots/nsub);
        for i=1:nrPlots
            if nrPlots~=1
                subplot(nsub,msub,i)
            end
            if nrPlots==1 & strcmp(plotInstr{i},'dmps')
                lowDP=3e-9;
            end
            eval(['plot_',plotInstr{i},'(AEstr,caxisMax,lowDP,hghDP);']);
        end
    case 'all'
        lowDP=0.3e-9;
        hghDP=1100e-9;
        
        subplot(3,2,1)
        plot_dmps(AEstr,caxisMax,lowDP,hghDP,plotFits);
        
        subplot(3,2,2)
        plot_comb(AEstr,caxisMax,lowDP,hghDP);
        
        subplot(3,2,3)
        plot_aisp(AEstr,caxisMax,lowDP,hghDP);
        
        subplot(3,2,4)
        plot_aisn(AEstr,caxisMax,lowDP,hghDP);
        
        subplot(3,2,5)
        plot_naisp(AEstr,caxisMax,lowDP,hghDP);
        
        subplot(3,2,6)
        plot_naisn(AEstr,caxisMax,lowDP,hghDP);
        
    case 'ais',
        lowDP=0.3e-9;
        hghDP=50e-9;
        
        subplot(2,1,1)
        plot_aisn(AEstr,caxisMax,lowDP,hghDP);
        
        subplot(2,1,2)
        plot_aisp(AEstr,caxisMax,lowDP,hghDP);
        
        
        
    case 'dmps',
        lowDP=2e-9;
        hghDP=1100e-9;
        
        plot_dmps(AEstr,caxisMax,lowDP,hghDP,plotFits)
        
    case 'aisn',
        lowDP=0.3e-9;
        hghDP=50e-9;
        
        plot_aisn(AEstr,caxisMax,lowDP,hghDP)
        
    case 'aisp',
        lowDP=0.3e-9;
        hghDP=50e-9;
        
        plot_aisp(AEstr,caxisMax,lowDP,hghDP)
        
    case 'naisn',
        lowDP=0.4e-9;
        hghDP=50e-9;
        
        plot_naisn(AEstr,caxisMax,lowDP,hghDP)
    case 'naisp',
        lowDP=0.4e-9;
        hghDP=50e-9;
        
        plot_naisp(AEstr,caxisMax,lowDP,hghDP)
    case 'nais',
        lowDP=0.4e-9;
        hghDP=50e-9;
        subplot(2,1,1)
        plot_naisn(AEstr,caxisMax,lowDP,hghDP)
        subplot(2,1,2)
        plot_naisp(AEstr,caxisMax,lowDP,hghDP)
    case 'nais4n',
        lowDP=0.3e-9;
        hghDP=50e-9;
        
        plot_nais4n(AEstr,caxisMax,lowDP,hghDP)
    case 'nais4p',
        lowDP=0.3e-9;
        hghDP=50e-9;
        
        plot_nais4p(AEstr,caxisMax,lowDP,hghDP)
    case 'nais4nA',
        lowDP=0.3e-9;
        hghDP=50e-9;
        
        plot_nais4nA(AEstr,caxisMax,lowDP,hghDP)
    case 'nais4pA',
        lowDP=0.3e-9;
        hghDP=50e-9;
        
        plot_nais4pA(AEstr,caxisMax,lowDP,hghDP)
        
    case 'nais4',
        lowDP=0.3e-9;
        hghDP=50e-9;
        subplot(2,2,1)
        plot_nais4n(AEstr,caxisMax,lowDP,hghDP)
        subplot(2,2,2)
        plot_nais4p(AEstr,caxisMax,lowDP,hghDP)
        subplot(2,2,3)
        plot_nais4nA(AEstr,caxisMax,lowDP,hghDP)
        subplot(2,2,4)
        plot_nais4pA(AEstr,caxisMax,lowDP,hghDP)
    case 'ions',
        lowDP=0.3e-9;
        hghDP=50e-9;
        
        subplot(2,2,1)
        plot_aisp(AEstr,caxisMax,lowDP,hghDP);
        
        subplot(2,2,2)
        plot_aisn(AEstr,caxisMax,lowDP,hghDP);
        
        subplot(2,2,3)
        plot_naisp(AEstr,caxisMax,lowDP,hghDP);
        
        subplot(2,2,4)
        plot_naisn(AEstr,caxisMax,lowDP,hghDP);
        
    case 'comb',
        lowDP=0.4e-9;
        hghDP=1100e-9;
        
        plot_comb(AEstr,caxisMax,lowDP,hghDP)
        
    case 'diff',
        lowDP=0.3e-9;
        hghDP=50e-9;
        
        tima=AEstr.aisp(2:end,1);
        dPa=AEstr.aisp(1,3:end);
        
        subplot(2,2,1)
        plot_aisp(AEstr,caxisMax,lowDP,hghDP)
        subplot(2,2,2)
        plot_aisn(AEstr,caxisMax,lowDP,hghDP)
        subplot(2,2,3)
        plot_dif(AEstr,caxisMax,lowDP,hghDP)
        
    case 'vol',
        lowDP=2e-9;
        hghDP=500e-9;
        
        subplot(2,1,1)
        plot_volC(AEstr,caxisMax,lowDP,hghDP)
        
        subplot(2,1,2)
        plot_volH(AEstr,caxisMax,lowDP,hghDP)
        
    case 'volc',
        lowDP=2e-9;
        hghDP=500e-9;
        
        plot_volC(AEstr,caxisMax,lowDP,hghDP)
        
    case 'volh',
        lowDP=2e-9;
        hghDP=500e-9;
        
        plot_volH(AEstr,caxisMax,lowDP,hghDP)
        
    case 'meas',
        lowDP=0.3e-9;
        hghDP=1100e-9;
        
        subplot(3,2,1)
        plot_dmps(AEstr,caxisMax,lowDP,hghDP,plotFits);
        
        subplot(3,2,3)
        plot_aisp(AEstr,caxisMax,lowDP,hghDP);
        
        subplot(3,2,4)
        plot_aisn(AEstr,caxisMax,lowDP,hghDP);
        
        subplot(3,2,5)
        plot_naisp(AEstr,caxisMax,lowDP,hghDP);
        
        subplot(3,2,6)
        plot_naisn(AEstr,caxisMax,lowDP,hghDP);
        
    case 'bsman',
        lowDP=0.3e-9;
        hghDP=50e-9;
        
        plot_bsman(AEstr,caxisMax,lowDP,hghDP)
    case 'bsmap',
        lowDP=0.3e-9;
        hghDP=50e-9;
        
        plot_bsmap(AEstr,caxisMax,lowDP,hghDP)
    case 'bsma',
        lowDP=0.3e-9;
        hghDP=50e-9;
        subplot(2,1,1)
        plot_bsman(AEstr,caxisMax,lowDP,hghDP)
        subplot(2,1,2)
        plot_bsmap(AEstr,caxisMax,lowDP,hghDP)
    case 'aisn_block'
        lowDP=0.3e-9;
        hghDP=50e-9;
        plot_aisn_block(AEstr,caxisMax,lowDP,hghDP)
    case 'aisp_block'
        lowDP=0.3e-9;
        hghDP=50e-9;
        plot_aisp_block(AEstr,caxisMax,lowDP,hghDP)
    case 'ais_block'
        lowDP=0.3e-9;
        hghDP=50e-9;
        subplot(2,1,1)
        plot_aisn_block(AEstr,caxisMax,lowDP,hghDP)
        subplot(2,1,2)
        plot_aisp_block(AEstr,caxisMax,lowDP,hghDP)
    case 'naisn_block'
        lowDP=1.3e-9;
        hghDP=50e-9;
        plot_naisn_block(AEstr,caxisMax,lowDP,hghDP)
    case 'naisp_block'
        lowDP=1.3e-9;
        hghDP=50e-9;
        plot_naisp_block(AEstr,caxisMax,lowDP,hghDP)
    case 'nais_block'
        lowDP=1.3e-9;
        hghDP=50e-9;
        subplot(2,1,1)
        plot_naisn_block(AEstr,caxisMax,lowDP,hghDP)
        subplot(2,1,2)
        plot_naisp_block(AEstr,caxisMax,lowDP,hghDP)
        
    otherwise
        disp('AE_plot: unknown plotType')
end


%SUBFUNCTIONS=======================================================
%% plot_dmps
function h=plot_dmps(AEstr,caxisMax,lowDP,hghDP,plotFits)
if isempty(AEstr.dmps)
    disp('AE_plot: no data');
    return
end

if nargin==4
    plotFits=0;
end
Iday=[];
if iscell(AEstr.dmps)
    if isempty(Iday)
        for i=1:length(AEstr.dmps)
            if ~isempty(AEstr.dmps{i})
                %             timd=AEstr.dmps{i}(2:end,1);
                timd=AEstr.meta.dmps.tim{i};
                
                %take 0.6 from size in order to match with AIS
                %         dPd=AEstr.dmps{i}(1,3:end)-0.6e-9;
                dPd=AEstr.dmps{i}(1,3:end);
                if all(dPd>100e-9)
                    %convert to meters
                    dPd=dPd*1e-9;
                end
                
                plotD=log10(max(AEstr.dmps{i}(2:end,3:end)',1e-10));
                plotD(isnan(AEstr.dmps{i}(2:end,3:end)'))=NaN;
                if all(plotD(:)==-10)
                    plotD=log10(AEstr.dmps{i}(2:end,3:end)');
                    %                 caxisMax=10*1e-18
                end
                pcolor(timd,dPd,plotD),shading flat
                colormap('jet')
                hold on
                if plotFits
                    try
                        if isfield(AEstr,'fit')
                            plot(AEstr.fit.dmps.tm{i},AEstr.fit.dmps.zs{i},'.k')
                        else
                            [Ps]=hm_pf(AEstr,'dmps',0,0);
                            plot(Ps(:,1),10.^Ps(:,2),'.k')
                        end
                    catch ME
                        disp(ME.message)
                    end
                end
            end
        end
    else
        i=Iday;
        if ~isempty(AEstr.dmps{i})
            %             timd=AEstr.dmps{i}(2:end,1);
            timd=AEstr.meta.dmps.tim{i};
            
            %take 0.6 from size in order to match with AIS
            %         dPd=AEstr.dmps{i}(1,3:end)-0.6e-9;
            dPd=AEstr.dmps{i}(1,3:end);
            if all(dPd>100e-9)
                %convert to meters
                dPd=dPd*1e-9;
            end
            
            plotD=log10(max(AEstr.dmps{i}(2:end,3:end)',1e-10));
            if all(plotD(:)==-10)
                plotD=log10(AEstr.dmps{i}(2:end,3:end)');
                %                 caxisMax=10*1e-18
            end
            pcolor(timd,dPd,plotD),shading flat
            colormap('jet')
            
            hold on
            if plotFits
                try
                    plot(hmD.fit.dmps.tm{i},hmD.fit.dmps.z{i},'.k')
                end
            end
        end
    end
    hold off
    title(['DMPS: ',AEstr.startTime,' - ',AEstr.endTime])
    shading flat; grid on
    set(gca,'YScale','log')
    %     axis([floor(AEstr.dmps{1}(2,1))-0.03 AEstr.dmps{i}(end,1)+0.03 lowDP hghDP])
    
    %     axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])
    
    caxis([1 caxisMax])
    %     timPer=(AEstr.meta.startTime)-(AEstr.meta.endTime);
else
    %     timd=AEstr.dmps(2:end,1);
    timd=AEstr.meta.dmps.tim{1};
    %take 0.6 from size in order to match with AIS
    %     dPd=AEstr.dmps(1,3:end)-0.6e-9;
    
    dPd=AEstr.dmps(1,3:end);
    if all(dPd>100e-9)
        %convert to meters
        dPd=dPd*1e-9;
    end
    
    plotD=log10(max(AEstr.dmps(2:end,3:end)',1e-10));
    pcolor(timd,dPd,plotD),shading flat
    colormap('jet')
    
    if ~isempty(AEstr.startTime)
        title(['DMPS: ',AEstr.startTime])
    else
        title(['DMPS: ',AEstr.srtWhen,' - ',AEstr.endWhen])
    end
    shading flat; grid on
    set(gca,'YScale','log')
    %     axis([AEstr.dmps(2,1)-0.03 AEstr.dmps(end,1)+0.03 0.4e-9 1e-6])
    %     axis([floor(AEstr.dmps(2,1))-0.03 AEstr.dmps(end,1)+0.03 lowDP hghDP])
    %     axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+0.03 lowDP hghDP])
    
    caxis([1 caxisMax])
    % timPer=AEstr.dmpsT(end,1)-AEstr.dmpsT(1,1);
    
end

timPer=(AEstr.meta.endTime)-(AEstr.meta.startTime);
axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])

if timd(1)>370 & timPer<2.5
    %     datetick('x',15)
    xdatetick(6)
elseif timd(1)>370 & timPer>=2.5
    %     datetick('x',19,'keepticks')
    %     datetick('x',19)
    xdatetick(24)
    datetick('x',7,'keepticks')
end
% 	caxis([1 5])
axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])
ylabel('Millikan diameter (m)')

%====================
%% plot_aisp
function plot_aisp(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.aisp)
    disp('AE_plot: no data');
    return
end

% %load external file with diameter conversions
%  dPa=[];
% convFile='diameters_AIS_mobility_tammet_milligan.dat';
% if exist(convFile)
%     %     first col mobility
%     %     second Tammet
%     %     third Millikan diameter
%     convD=load(convFile);
%     dPa=convD(:,3)';%Millikan
% end
%
for i=1:length(AEstr.aisp)
    if ~isempty(AEstr.aisp{i})
        %         tima=AEstr.aisp{i}(2:end,1);
        tima=AEstr.meta.aisp.tim{i};
        dPa=AEstr.meta.aisp.dp{i}(1,:);
        
        plotD=log10(max(AEstr.aisp{i}(2:end,3:end)',1e-10));
        %             plotD=log10(AEstr.aisp{i}(2:end,3:end)');
        pcolor(tima,dPa,plotD)
        colormap('jet')
        
        hold on
    end
end
hold off
title(['AIS (+): ',AEstr.startTime,' - ',AEstr.endTime])
shading flat; grid on
set(gca,'YScale','log')
%     axis([floor(AEstr.aisp{1}(2,1))-0.03 AEstr.aisp{i}(end,1)+0.03 lowDP hghDP])
%          caxis([1 caxisMax]);
%          axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])

%     xdatetick(6)

timPer=(AEstr.meta.endTime)-(AEstr.meta.startTime);

%     xdatetick(6)
if tima(1)>370 & timPer<2.5
    datetick('x',15)
elseif tima(1)>370 & timPer>=2.5
    %     datetick('x',19,'keepticks')
    %     datetick('x',19)
    xdatetick(24)
    datetick('x',7,'keepticks')
end
axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])
caxis([1 caxisMax]);
xdatetick(6)
ylabel('Millikan diameter (m)')

%====================
%% plot_aisn
function plot_aisn(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.aisn)
    disp('AE_plot: no data');
    return
end

% %load external file with diameter conversions
% dPa=[];
% convFile='diameters_AIS_mobility_tammet_milligan.dat';
% if exist(convFile)
%     %     first col mobility
%     %     second Tammet
%     %     third Milligan diameter
%     convD=load(convFile);
%     dPa=convD(:,3)';%Milligan
% end
%
for i=1:length(AEstr.aisn)
    if ~isempty(AEstr.aisn{i})
        %         tima=AEstr.aisn{i}(2:end,1);
        tima=AEstr.meta.aisn.tim{i};
        dPa=AEstr.meta.aisn.dp{i}(1,:);
        plotD=log10(max(AEstr.aisn{i}(2:end,3:end)',1e-10));
        pcolor(tima,dPa,plotD)
        colormap('jet')
        
        hold on
    end
end
hold off
title(['AIS (-): ',AEstr.startTime,' - ',AEstr.endTime])
shading flat; grid on
set(gca,'YScale','log')


timPer=(AEstr.meta.endTime)-(AEstr.meta.startTime);

if tima(1)>370 & timPer<2.5
    datetick('x',15)
elseif tima(1)>370 & timPer>=2.5
    %     datetick('x',19,'keepticks')
    %     datetick('x',19)
    xdatetick(24)
    datetick('x',7,'keepticks')
end
axis([floor(AEstr.meta.startTime)-0.04 floor(AEstr.meta.endTime)+1.04 lowDP hghDP])
caxis([1 caxisMax]);
xdatetick(6)
ylabel('Millikan diameter (m)')
%====================
%% plot_ais
function plot_ais(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.ais)
    disp('AE_plot: no data');
    return
end

tima=AEstr.ais(2:end,1);
dPa=AEstr.ais(1,3:end);

plotD=log10(max(AEstr.ais(2:end,3:end)',1e-10));
pcolor(tima,dPa,plotD),shading flat
colormap('jet')

title(['AIS total: ',AEstr.startTime])
shading flat; grid on
set(gca,'YScale','log')
axis([AEstr.aisn(2,1)-0.03 AEstr.aisn(end,1)+0.03 lowDP hghDP])
caxis([1 caxisMax])
%colorbar

%====================
%% plot_dif
function plot_dif(AEstr,caxisMax,lowDP,hghDP)
tima=AEstr.aisn(2:end,1);
dPa=AEstr.aisn(1,3:end);

% 	plotD=log10(max(AEstr.aisn(2:end,3:end)',1e-10))-log10(max(AEstr.aisp(2:end,3:end)',1e-10));
plotD=max(AEstr.aisn(2:end,3:end)',1e-10)-max(AEstr.aisp(2:end,3:end)',1e-10);
pcolor(tima,dPa,plotD),shading flat
colormap('jet')

title('AIS difference; neg-pos')
shading flat; grid on
set(gca,'YScale','log')
axis([AEstr.aisn(2,1)-0.03 AEstr.aisn(end,1)+0.03 lowDP hghDP])
caxis([-std(plotD(:))*2 std(plotD(:))*2])

%====================
%% plot_naisp
function plot_naisp(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.naisp)
    disp('AE_plot: no data');
    return
end
% %load external file with diameter conversions
% dPa=[];
% convFile='diameters_NAIS_mobility_tammet_milligan.dat';
% if exist(convFile)
%     %     first col mobility
%     %     second Tammet
%     %     third Milligan diameter
%     convD=load(convFile);
%     dPa=convD(:,3)';%Milligan
% end
%
for i=1:length(AEstr.naisp)
    if ~isempty(AEstr.naisp{i})
        tima=AEstr.meta.naisp.tim{i};
        dPa=AEstr.meta.naisp.dp{i}(1,:);
        dat=AEstr.naisp{i}(2:end,3:end);
        
        plotD=log10(max(dat',1e-10));
        pcolor(tima,dPa,plotD),shading flat
        colormap('jet')
        
        hold on
    end
end
hold off
if i>1
    title(['NAIS: (+) ',AEstr.startTime,' - ',AEstr.endTime])
else
    title(['NAIS: (+) ',AEstr.startTime])
end
shading flat; grid on
set(gca,'YScale','log')
axis([floor(AEstr.meta.naisp.tim{1}(1))-0.03 AEstr.meta.naisp.tim{end}(end)+0.03  lowDP hghDP])
caxis([1 caxisMax])


timPer=(AEstr.meta.endTime)-(AEstr.meta.startTime);

if tima(1)>370 & timPer<2.5
    %     datetick('x',15)
    xdatetick(6)
elseif tima(1)>370 & timPer>=2.5
    %     datetick('x',19,'keepticks')
    %     datetick('x',19)
    xdatetick(24)
    datetick('x',7,'keepticks')
end
ylabel('Millikan diameter (m)')


%colorbar

%====================
%% plot_naisn
function plot_naisn(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.naisn)
    disp('AE_plot: no data');
    return
end

%load external file with diameter conversions
% dPa=[];
% convFile='diameters_NAIS_mobility_tammet_milligan.dat';
% if exist(convFile)
%     %     first col mobility
%     %     second Tammet
%     %     third Milligan diameter
%     convD=load(convFile);
%     dPa=convD(:,3)';%Milligan
% end

for i=1:length(AEstr.naisn)
    if ~isempty(AEstr.naisn{i})
        tima=AEstr.meta.naisn.tim{i};
        
        dPa=AEstr.meta.naisn.dp{i}(1,:);
        dat=AEstr.naisn{i}(2:end,3:end);
        
        plotD=log10(max(dat',1e-10));
        % plotD=(max(AEstr.naisn(2:end,3:end)',1e-10));
        pcolor(tima,dPa,plotD),shading flat
        colormap('jet')
        
        hold on
    end
end
hold off

if i>1
    title(['NAIS: (-) ',AEstr.startTime,' - ',AEstr.endTime])
else
    title(['NAIS: (-) ',AEstr.startTime])
end

shading flat; grid on
set(gca,'YScale','log')
axis([floor(AEstr.meta.naisn.tim{1}(1))-0.03 AEstr.meta.naisn.tim{end}(end)+0.03  lowDP hghDP])
caxis([1 caxisMax])
%colorbar

timPer=(AEstr.meta.endTime)-(AEstr.meta.startTime);

if tima(1)>370 & timPer<2.5
    %     datetick('x',15)
    xdatetick(6)
elseif tima(1)>370 & timPer>=2.5
    %     datetick('x',19,'keepticks')
    %     datetick('x',19)
    xdatetick(24)
    datetick('x',7,'keepticks')
end
ylabel('Millikan diameter (m)')

%====================
%% plot_nais4p
function plot_nais4p(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.nais4p)
    disp('AE_plot: no data');
    return
end

tima=AEstr.nais4p(2:end,1);
dPa=AEstr.nais4p(1,3:end);

plotD=log10(max(AEstr.nais4p(2:end,3:end)',1e-10));
pcolor(tima,dPa,plotD),shading flat
colormap('jet')

title(['NAIS ion: (+) ',AEstr.startTime])
shading flat; grid on
set(gca,'YScale','log')
axis([AEstr.nais4p(2,1)-0.03 AEstr.nais4p(end,1)+0.03  lowDP hghDP])
caxis([1 caxisMax])
%colorbar
%====================
%% plot_nais4n
function plot_nais4n(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.nais4n)
    disp('AE_plot: no data');
    return
end

tima=AEstr.nais4n(2:end,1);
dPa=AEstr.nais4n(1,3:end);

plotD=log10(max(AEstr.nais4n(2:end,3:end)',1e-10));
pcolor(tima,dPa,plotD),shading flat
colormap('jet')

title(['NAIS ion: (-) ',AEstr.startTime])
shading flat; grid on
set(gca,'YScale','log')
axis([AEstr.nais4n(2,1)-0.03 AEstr.nais4n(end,1)+0.03  lowDP hghDP])
caxis([1 caxisMax])
%colorbar
%====================
%% plot_nais4pA
function plot_nais4pA(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.nais4pA)
    disp('AE_plot: no data');
    return
end

tima=AEstr.nais4pA(2:end,1);
dPa=AEstr.nais4pA(1,3:end);

plotD=log10(max(AEstr.nais4pA(2:end,3:end)',1e-10));
pcolor(tima,dPa,plotD),shading flat
colormap('jet')

title(['NAIS particles: (+) ',AEstr.startTime])
shading flat; grid on
set(gca,'YScale','log')
axis([AEstr.nais4pA(2,1)-0.03 AEstr.nais4pA(end,1)+0.03  lowDP hghDP])
caxis([1 caxisMax])
%colorbar
%====================
%% plot_nais4nA
function plot_nais4nA(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.nais4nA)
    disp('AE_plot: no data');
    return
end

tima=AEstr.nais4nA(2:end,1);
dPa=AEstr.nais4nA(1,3:end);

plotD=log10(max(AEstr.nais4nA(2:end,3:end)',1e-10));
pcolor(tima,dPa,plotD),shading flat
colormap('jet')

title(['NAIS particles: (-) ',AEstr.startTime])
shading flat; grid on
set(gca,'YScale','log')
axis([AEstr.nais4nA(2,1)-0.03 AEstr.nais4nA(end,1)+0.03  lowDP hghDP])
caxis([1 caxisMax])
%colorbar

%====================
%% plot_volC
function plot_volC(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.volC)
    disp('AE_plot: no data');
    return
end
timd=AEstr.volC(2:end,1);
dPd=AEstr.volC(1,3:end);
caxisMax,lowDP,hghDP=4;

plotD=log10(max(AEstr.volC(2:end,3:end)',1e-10));
pcolor(timd,dPd,plotD),shading flat
colormap('jet')

title(['Volatility, cold: ',AEstr.startTime])
shading flat; grid on
set(gca,'YScale','log')
axis([AEstr.dmps(2,1)-0.03 AEstr.dmps(end,1)+0.03  lowDP hghDP])
% axis([AEstr.volC(2,1)-0.03 AEstr.volC(end,1)+0.03 2e-9 50e-9])
caxis([1 caxisMax])

%====================
%% plot_volH
function plot_volH(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.volH)
    disp('AE_plot: no data');
    return
end

timd=AEstr.volH(2:end,1);
dPd=AEstr.volH(1,3:end)-0.6e-9;

plotD=log10(max(AEstr.volH(2:end,3:end)',1e-10));
pcolor(timd,dPd,plotD),shading flat
colormap('jet')

title(['Volatility, hot: ',AEstr.startTime])
shading flat; grid on
set(gca,'YScale','log')
axis([AEstr.dmps(2,1)-0.03 AEstr.dmps(end,1)+0.03  lowDP hghDP])
% axis([AEstr.volH(2,1)-0.03 AEstr.volH(end,1)+0.03 2e-9 50e-9])
caxis([1 caxisMax])

%====================
%% plot_comb
function plot_comb(AEstr,caxisMax,lowDP,hghDP)
if isempty(AEstr.comb)
    AEstr=AE_interp(AEstr);
end

plotDc=log10(max(AEstr.comb(2:end,3:end),1));
timc=repmat(AEstr.comb(2:end,1),1,length(AEstr.comb(1,3:end)));
dPc=repmat(AEstr.comb(1,3:end),length(AEstr.comb(2:end,1)),1);
pcolor(timc,dPc,plotDc)
colormap('jet')

title('Combined; dmps+ais')
shading flat; grid on
set(gca,'YScale','log')
axis([AEstr.comb(2,1)-0.03 AEstr.comb(end,1)+0.03 lowDP hghDP])
caxis([1 caxisMax])

%====================
%% plot_bsman
function plot_bsman(AEstr,caxisMax,lowDP,hghDP)

inst='bsman';
eval(['instDat=AEstr.',inst,';'])
if ~iscell(instDat)
    plotDc=log10(max(instDat(2:end,3:end),1));
    
    eval(['dp=AEstr.meta.',inst,'.dp;'])
    eval(['tim=AEstr.meta.',inst,'.tim{1};'])
    
    timc=repmat(tim,1,length(dp));
    dPc=repmat(dp,length(tim),1);
    
    pcolor(timc,dPc,plotDc)
    colormap('jet')
    
    shading flat; grid on
    set(gca,'YScale','log')
    axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])
    xdatetick(6)
    
    caxis([1 caxisMax])
    switch inst
        case 'bsman'
            title([upper(inst(1:4)),' (-) : ',AEstr.startTime])
        case 'bsmap'
            title([upper(inst(1:4)),' (+) : ',AEstr.startTime])
    end
    
else
    for d=1:length(instDat)
        eval(['p = size (AEstr.',inst,');'])
        plotDc=log10(max(instDat{d}(2:end,3:end),1));
        
        eval(['dp=AEstr.meta.',inst,'.dp;'])
        eval(['tim=AEstr.meta.',inst,'.tim{d};'])
        
        timc=repmat(tim,1,length(dp));
        dPc=repmat(dp,length(tim),1);
        
        pcolor(timc,dPc,plotDc)
        colormap('jet')
        
        hold on
    end
    hold off
    
    shading flat; grid on
    set(gca,'YScale','log')
    axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])
    xdatetick(6)
    
    caxis([1 caxisMax])
    switch inst
        case 'bsman'
            title([upper(inst(1:4)),' (-) : ',AEstr.startTime,' - ',AEstr.endTime])
        case 'bsmap'
            title([upper(inst(1:4)),' (+) : ',AEstr.startTime,' - ',AEstr.endTime])
    end
end


%====================
%% plot_bsmap
function plot_bsmap(AEstr,caxisMax,lowDP,hghDP)

inst='bsmap';
eval(['instDat=AEstr.',inst,';'])
if ~iscell(instDat)
    plotDc=log10(max(instDat(2:end,3:end),1));
    
    eval(['dp=AEstr.meta.',inst,'.dp;'])
    eval(['tim=AEstr.meta.',inst,'.tim{1};'])
    
    timc=repmat(tim,1,length(dp));
    dPc=repmat(dp,length(tim),1);
    
    pcolor(timc,dPc,plotDc)
    colormap('jet')
    
    shading flat; grid on
    set(gca,'YScale','log')
    axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])
    xdatetick(6)
    
    caxis([1 caxisMax])
    switch inst
        case 'bsman'
            title([upper(inst(1:4)),' (-) : ',AEstr.startTime])
        case 'bsmap'
            title([upper(inst(1:4)),' (+) : ',AEstr.startTime])
    end
    
else
    for d=1:length(instDat)
        plotDc=log10(max(instDat{d}(2:end,3:end),1));
        
        eval(['dp=AEstr.meta.',inst,'.dp;'])
        eval(['tim=AEstr.meta.',inst,'.tim{d};'])
        
        timc=repmat(tim,1,length(dp));
        dPc=repmat(dp,length(tim),1);
        
        pcolor(timc,dPc,plotDc)
        colormap('jet')
        
        hold on
    end
    hold off
    
    shading flat; grid on
    set(gca,'YScale','log')
    axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])
    xdatetick(6)
    
    caxis([1 caxisMax])
    switch inst
        case 'bsman'
            title([upper(inst(1:4)),' (-) : ',AEstr.startTime,' - ',AEstr.endTime])
        case 'bsmap'
            title([upper(inst(1:4)),' (+) : ',AEstr.startTime,' - ',AEstr.endTime])
    end
    
end

function plot_aisp_block(AEstr,caxisMax,lowDP,hghDP)
for i=1:length(AEstr.aisp_block)
    if ~isempty(AEstr.aisp_block{i})
        %         tima=AEstr.aisp{i}(2:end,1);
        tima=AEstr.meta.aisp_block.tim{i};
        dPa=AEstr.meta.aisp_block.dp{i}(1,:);
        
        plotD=log10(max(AEstr.aisp_block{i}(2:end,3:end)',1e-10));
        %             plotD=log10(AEstr.aisp{i}(2:end,3:end)');
        pcolor(tima,dPa,plotD)
        colormap('jet')
        
        hold on
    end
end
hold off
title(['AIS (+): ',AEstr.startTime,' - ',AEstr.endTime])
shading flat; grid on
set(gca,'YScale','log')
%     axis([floor(AEstr.aisp{1}(2,1))-0.03 AEstr.aisp{i}(end,1)+0.03 lowDP hghDP])
%          caxis([1 caxisMax]);
%          axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])

%     xdatetick(6)

timPer=(AEstr.meta.endTime)-(AEstr.meta.startTime);

%     xdatetick(6)
if tima(1)>370 & timPer<2.5
    datetick('x',15)
elseif tima(1)>370 & timPer>=2.5
    %     datetick('x',19,'keepticks')
    %     datetick('x',19)
    xdatetick(24)
    datetick('x',7,'keepticks')
end
axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])
caxis([1 caxisMax]);
xdatetick(6)
ylabel('Millikan diameter (m)')

function plot_aisn_block(AEstr,caxisMax,lowDP,hghDP)
for i=1:length(AEstr.aisn_block)
    if ~isempty(AEstr.aisn_block{i})
        %         tima=AEstr.aisp{i}(2:end,1);
        tima=AEstr.meta.aisn_block.tim{i};
        dPa=AEstr.meta.aisn_block.dp{i}(1,:);
        
        plotD=log10(max(AEstr.aisn_block{i}(2:end,3:end)',1e-10));
        %             plotD=log10(AEstr.aisp{i}(2:end,3:end)');
        pcolor(tima,dPa,plotD)
        colormap('jet')
        
        hold on
    end
end
hold off
title(['AIS (-): ',AEstr.startTime,' - ',AEstr.endTime])
shading flat; grid on
set(gca,'YScale','log')
%     axis([floor(AEstr.aisp{1}(2,1))-0.03 AEstr.aisp{i}(end,1)+0.03 lowDP hghDP])
%          caxis([1 caxisMax]);
%          axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])

%     xdatetick(6)

timPer=(AEstr.meta.endTime)-(AEstr.meta.startTime);

%     xdatetick(6)
if tima(1)>370 & timPer<2.5
    datetick('x',15)
elseif tima(1)>370 & timPer>=2.5
    %     datetick('x',19,'keepticks')
    %     datetick('x',19)
    xdatetick(24)
    datetick('x',7,'keepticks')
end
axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])
caxis([1 caxisMax]);
xdatetick(6)
ylabel('Millikan diameter (m)')


function plot_naisp_block(AEstr,caxisMax,lowDP,hghDP)
for i=1:length(AEstr.naisp_block)
    if ~isempty(AEstr.naisp_block{i})
        %         tima=AEstr.aisp{i}(2:end,1);
        tima=AEstr.meta.naisp_block.tim{i};
        dPa=AEstr.meta.naisp_block.dp{i}(1,:);
        
        plotD=log10(max(AEstr.naisp_block{i}(2:end,3:end)',1e-10));
        %             plotD=log10(AEstr.aisp{i}(2:end,3:end)');
        pcolor(tima,dPa,plotD)
        colormap('jet')
        
        hold on
    end
end
hold off
title(['AIS (+): ',AEstr.startTime,' - ',AEstr.endTime])
shading flat; grid on
set(gca,'YScale','log')
%     axis([floor(AEstr.aisp{1}(2,1))-0.03 AEstr.aisp{i}(end,1)+0.03 lowDP hghDP])
%          caxis([1 caxisMax]);
%          axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])

%     xdatetick(6)

timPer=(AEstr.meta.endTime)-(AEstr.meta.startTime);

%     xdatetick(6)
if tima(1)>370 & timPer<2.5
    datetick('x',15)
elseif tima(1)>370 & timPer>=2.5
    %     datetick('x',19,'keepticks')
    %     datetick('x',19)
    xdatetick(24)
    datetick('x',7,'keepticks')
end
axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])
caxis([1 caxisMax]);
xdatetick(6)
ylabel('Millikan diameter (m)')

function plot_naisn_block(AEstr,caxisMax,lowDP,hghDP)
for i=1:length(AEstr.naisn_block)
    if ~isempty(AEstr.naisn_block{i})
        %         tima=AEstr.aisp{i}(2:end,1);
        tima=AEstr.meta.naisn_block.tim{i};
        dPa=AEstr.meta.naisn_block.dp{i}(1,:);
        
        plotD=log10(max(AEstr.naisn_block{i}(2:end,3:end)',1e-10));
        %             plotD=log10(AEstr.aisp{i}(2:end,3:end)');
        pcolor(tima,dPa,plotD)
        colormap('jet')
        
        hold on
    end
end
hold off
title(['NAIS (-): ',AEstr.startTime,' - ',AEstr.endTime])
shading flat; grid on
set(gca,'YScale','log')
%     axis([floor(AEstr.aisp{1}(2,1))-0.03 AEstr.aisp{i}(end,1)+0.03 lowDP hghDP])
%          caxis([1 caxisMax]);
%          axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])

%     xdatetick(6)

timPer=(AEstr.meta.endTime)-(AEstr.meta.startTime);

%     xdatetick(6)
if tima(1)>370 & timPer<2.5
    datetick('x',15)
elseif tima(1)>370 & timPer>=2.5
    %     datetick('x',19,'keepticks')
    %     datetick('x',19)
    xdatetick(24)
    datetick('x',7,'keepticks')
end
axis([floor(AEstr.meta.startTime)-0.03 floor(AEstr.meta.endTime)+1.03 lowDP hghDP])
caxis([1 caxisMax]);
xdatetick(6)
ylabel('Millikan diameter (m)')

%=======================
%% xdatetick
function xdatetick(int)
%set x tick to sertain interval (hours)
%

xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
set(gca,'xtick',floor(xlim(1)):int/24:ceil(xlim(2)))
datetick('x',15,'keepticks')
set(gca,'xlim',xlim,'ylim',ylim)
grid on