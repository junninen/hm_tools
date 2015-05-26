function [param,fval,y1_F,h_out,N]=MF_lognorm(dat,dp,doPlot,mxNrM,init)
%
% fit lognormal distribution to size distribution
%
% [param,fval,y1_F,h_out,N]=MF_lognorm(dat,dp,[doPlot],[mxNrM],[initPar])
%
% param (1x4 vector) parameters, [width1,width2,gmd1,gmd2]
% fval  (vector)     index of fit
% yhat  (multidim matrix) extracted modes
%
% dat   (vector)    data
% dp    (vector)    same size as dat, size bins
% doPlot (double)   1=plot result, 0 = dont plot anything
% mxNrM  (double)   number of maximum modes extracted, default=3
% init (vector)     initialization

%Heikki Junninen 10.04.2007

%defaults and criterias
maxIter=max(mxNrM-1,1); %5

global scaleNegErr
global smlLimR
global smlLimA
global errLim
global rmseLim
global armseLim
global okRangeCoef
% global smallAeroWeight % below 20nm fit is multiplied by constant
global maxIt
global sameCntLim
global minErr
global rang
global doIterPlot
global doIterDisp

% scaleNegErr=1.6; %4 %scale negative residuals by factor
scaleNegErr=2.6; %4 %scale negative residuals by factor
% scaleNegErr=1; %4 %scale negative residuals by factor

smlLimR=0.001; %relative peak area approximation limit
%0.05 smaller than this relative conc mode will be removed
smlLimA=200; %200 smaller than this absolute conc will be removed
errLim=0.4; %0.4 modes that have local error bigger than that will be removed

rmseLim=50; %50
armseLim=0.001; %0.01

okRangeCoef=1.0; %1.3,1.8
% smallAeroWeight=2;

%inside mlr_iteration parameters
maxIt=60; %number of mlr_iterations
sameCntLim=2; %max number of no improvement iterations
minErr=0.04; %error limit
rang=0.3; %range parameter, how big is the window to search
doIterPlot=0; %plot each iteration step
doIterDisp=0; %display stoping criterias

%% check inputs
doinit=0;
if nargin<=3
    mxNrM=3;
    doinit=1;
end

if nargin<=4
    doinit=1;
end

dat=dat(:)';

%if not provided do initialisation
%make initial chromosome
%equally distributed z values and few sigmas
%starting with 1 mode

if doinit
    nrPeaks=1;
    [param1]=makeGood1modeInit(dat,dp,3);
end

if ~doinit
    nrPeaks=length(init)/2;
    param1=init;
end

%check goodness and evaluate the first step
[param1,fval,yhat,yshat,nrPeaks]=mlr_iteration(dat,dp,nrPeaks,param1,10);

% if positive residual is more than 10% of original data fit new mode
%below 1nm dont add new mode
param_new=param1;
fval_new=fval;
% Iabove1nm=find(dp>1e-9); %index for values above 1nm
% notTooClose=1; %new mode have to be further than 1 sigma

% resDat=dat-yhat;
[rmse,armse,areaAprox,resDat]=fitness(dat,yhat,dp);

nrIter=0;

if (rmse>rmseLim | armse>armseLim) & any(areaAprox>smlLimA)
    do=1;
    %if max number of modes already from init, but errors are big
    %delete smallest one, in order to make new fit
    if nrPeaks==mxNrM
        [mn,Isml]=min(sum(yshat'));
        param_new([Isml,Isml+nrPeaks])=[];
        nrPeaks=nrPeaks-1;
    end
else
    do=0;
end

while nrPeaks<mxNrM & do
    nrIter=nrIter+1;
    nrP=nrPeaks;
    %     errs(nrIter,:)=[fval_new,rmse,armse];
    if nrP>1
        [param_rem,nrPeaks_rem]=rem_mode(dat,dp,nrPeaks,param_new);
        if nrPeaks_rem~=nrPeaks
            %new iter with bad mode removed
            [param_rem,fval_rem,yhat_rem,yshat,nrPeaks_rem]=mlr_iteration(dat,dp,nrPeaks_rem,param_rem,scaleNegErr);
            [rmse_rem,armse_rem,areaAprox_rem,resDat]=fitness(dat,yhat_rem,dp);
            [rmse,armse,areaAprox,resDat]=fitness(dat,yhat,dp);

            if rmse_rem<rmse
                param_new=param_rem;
                yhat=yhat_rem;
                nrPeaks=nrPeaks_rem;
                fval_new=fval_rem;
            end
        end
    end

    %add new mode
    if (rmse>rmseLim | armse>armseLim) | any(areaAprox>smlLimA)
        [param_add,nrPeaks_add]=add_mode(dat,resDat,dp,nrPeaks,param_new);
        %if add_mode-function did not add a mode, means there is no place to
        %add
        if nrPeaks_add~=nrPeaks
            [param_add,fval_add,yhat_add,yshat,nrPeaks_add]=mlr_iteration(dat,dp,nrPeaks_add,param_add,scaleNegErr);
            %             hold all, plot(yhat),plot(dat,'k')
            %             %compare fitness
            [rmse,armse,areaAprox,resDat]=fitness(dat,yhat,dp);
            [rmse_add,armse_add,areaAprox_add,resDat]=fitness(dat,yhat_add,dp);

            if rmse_add<rmse
                param_new=param_add;
                yhat=yhat_add;
                nrPeaks=nrPeaks_add;
                fval_new=fval_add;
            end
        else
            %if cant add mode terminate the loop
            do=0;
        end

    end

    [rmse,armse,areaAprox,resDat]=fitness(dat,yhat,dp);

    %stop criterions
    %nrP==nrPeaks |
    if  ((rmse<rmseLim | armse<armseLim) & all(areaAprox<smlLimA))
        do=0;
        %         disp(['Error:',num2str(rmse),'(',num2str(rmseLim),')'])
        %         disp(['Absolute error:',num2str(armse),'(',num2str(armseLim),')'])
        %         disp(['Area:',num2str(areaAprox),'(',num2str(smlLimA),')'])
    end
    if nrIter>maxIter
        do=0;
        %         disp(['iterations:', num2str(maxIter)])
    end
end

rmse=[];
if fval_new(end)<fval
    param=param_new;

else
    param=param1;
    %     nrPeaks=1;
    %     nrPeaks=length(param)/2;
end
nrPeaks=length(param)/2;

%check if there are bad modes
[param_old,fval_old,yhat_old,yshat,nrPeaks_old]=mlr_iteration(dat,dp,nrPeaks,param,scaleNegErr);
[param_rem,nrPeaks_rem]=rem_mode(dat,dp,nrPeaks,param);
if nrPeaks_rem~=nrPeaks_old
    %new iter with bad mode removed and then stop the while loop
    fval_old=fval;
    [param_rem,fval_rem,yhat_rem,yshat,nrPeaks]=mlr_iteration(dat,dp,nrPeaks_rem,param_rem,scaleNegErr);
    [rmse_new,armse,areaAprox,resDat]=fitness(dat,yhat_rem,dp);
    if fval_old>fval_rem
        param=param_rem;
    else
        param=param_old;
        fval=fval_old;
        yhat=yhat_old;
        nrPeaks=nrPeaks_old;
    end
end

nrPeaks=length(param)/2;
[fvalF,y1_F,yF,h_chF]=H_lognorm2pAE4(dat,dp,nrPeaks,param,scaleNegErr);
[rmse,armse,areaAprox,resDat]=fitness(dat,yF,dp);

 [par1]=nlinOptim(dat,dp,param);
% par1=param;
param=par1;
[fvalF_nl,y1_F_nl,yF_nl,h_chF_nl]=H_lognorm2pAE4(dat,dp,nrPeaks,par1,scaleNegErr);

% disp(['Weighted error: ',num2str(fval)])
% disp(['         Error: ',num2str(rmse),'(',num2str(rmseLim),')'])
% disp(['Absolute error: ',num2str(armse),'(',num2str(armseLim),')'])
% disp(['          Area: ',num2str(areaAprox),'(',num2str(smlLimA),')'])

h_out=h_chF;

calcNumConc=1;
N=h_out;
if calcNumConc
    dmin=log10(dp(1));
    dmax=log10(dp(end));
    logDp=log10(dp);

    dpi=[dmin:0.001:dmax]';
    for i=1:nrPeaks
        N(i)=sum(interp1(logDp',squeeze(y1_F(:,i,:)),dpi)*0.001);
    end
end

%% doPlot
if doPlot
    plot(dp,squeeze(y1_F(:,1,:))','r')
    hold on
    for i=2:nrPeaks
        plot(dp,squeeze(y1_F(:,i,:))','r')
    end
    plot(dp,yF,'k','linewidth',1.5)

    plot(dp,[dat],'.')
    hold off
    set(gca,'xscale','log')
    drawnow
end

%==========================================================================
%% SUBFUNCTIONS

%% makeGood1modeInit
function [param1]=makeGood1modeInit(dat,dp,scaleNegErr)
%
%make population for one mode that covers whole space
%and select the best one
%

scaleNegErr=10; %alwais the same in first initialzation

nrPeaks=1;


% dp1=log(1); %in nm
% dp2=log(400); %in nm

dp1=log((dp(1)*1e9/2)); %in nm
dp2=log((dp(end)*1e9)); %in nm

step=(dp2-dp1)/70;

zs=exp([dp1:step:dp2])';

ch=[ones(size(zs,1),nrPeaks)*1.13,zs*1e-9];
ch=[ones(size(zs,1),nrPeaks)*1.2,zs*1e-9];
ch=[ch;ones(size(zs,1),nrPeaks)*1.42,zs*1e-9];
ch=[ch;ones(size(zs,1),nrPeaks)*1.64,zs*1e-9];
ch=[ch;ones(size(zs,1),nrPeaks)*1.9,zs*1e-9];
ch=[ch;ones(size(zs,1),nrPeaks)*2.2,zs*1e-9];
%   ch=[ch;ones(size(zs,1),nrPeaks)*1.01,zs*1e-9];
% ch=[ch;ones(size(zs,1),nrPeaks)*1.03,zs*1e-9];


[fval,y1,y,h_ch]=H_lognorm2pAE4(dat,dp,nrPeaks,ch,scaleNegErr);
[mn,Imn]=min(fval);
yhat=y(Imn,:);
fval=fval(Imn);
param1=ch(Imn,:);

% h=0;
% for i=1:0.2:3
% h=h+1;
%     [fval,y1,y,h_ch]=H_lognorm2pAE4(dat,dp,nrPeaks,ch,i);
%     [mn(h),Imn(h)]=min(fval);
%     yhat=y(Imn(h),:);
%     [rmse(h),armse,areaAprox,resDat]=fitness(dat,yhat,dp);
%     figure,plot([dat;yhat]')
% end
%% add_mode
function [param_new,nrPeaks]=add_mode(dat,resDat,dp,nrPeaks,param)
%
%
%


%fit just positive residuals
global okRangeCoef
% global smlLimR
global smlLimA
% global errLim
% global rmseLim
% global armseLim


if 1
    %just find where is smoothed maximum in residual
    resDat=medfilt1(max(resDat,0),5); %only positive residuals
    resDats=resDat;
    %     resDats=(resDat(1:end-2)+resDat(2:end-1)+resDat(3:end))/3;
    %     resDats=[resDat(1),resDats,resDat(end)];

    %     %weigh residuals bellow 20nm
    %     Isml=dp<20e-9 & dp>3e-9;
    %     resDat(Isml)=resDat(Isml)*3;

    %     plot(resDat),hold all,plot(resDats),
    %find the solutions only one sigma away from any existing modes

    lims=zeros(nrPeaks,2);
    for i=1:nrPeaks
        %         ws=param(i);
        ws=okRangeCoef*param(i);
        zs=param(nrPeaks+i);
        % plot(10.^([log10(param1(2))-log10(param1(1)),log10(param1(2))+log
        lims(i,:)=10.^([log10(zs)-log10(ws),log10(zs)+log10(ws)]);
        Idel=dp>lims(i,1) & dp<lims(i,2);
        resDats(Idel)=0;
        %         plot(resDats),
    end

    %calculate sums of all adjacent positive residuals and select the one
    %that has higer sum
    [areaAprox,Isrt,Istp]=calc_residual_areas(resDats);
    [mxa,Imxa]=max(areaAprox);

    mask=zeros(size(resDats));
    mask(Isrt(Imxa):Istp(Imxa))=1;
    [mx,Imx]=max(mask.*resDats); %select max point

    if mxa>smlLimA
        %set all widths to 1.4
        %         param_new=[ones(1,nrPeaks)*1.4,1.4,param(nrPeaks+1:end),dp(Imx)];
        %use already fitted values
        param_new=[param(1:nrPeaks),1.4,param(nrPeaks+1:end),dp(Imx)];
        nrPeaks=nrPeaks+1;
    else
        param_new=param;
    end
end

%% calc_residual_areas
function [areaAprox,Isrt,Istp]=calc_residual_areas(resDats)
%
% calculate areas of adjacent positive regions
%
areaAprox=[];
Isrt=[];
Istp=[];

Ipos=resDats>0;

Iedge=[Ipos(2:end)-Ipos(1:end-1)];

Inz=Iedge(Iedge~=0);

if isempty(Inz)
    Iedge(1)=1;
    Iedge(end)=-1;
else
    if Inz(1)==-1
        Iedge(1)=1;
    end

    if Inz(end)==1
        Iedge(end)=-1;
    end
end
Ie=find(Iedge~=0);
Idel=diff(Iedge(Ie))==0;
Iedge(Ie(Idel))=0;

Isrt=find(Iedge==1);
Istp=find(Iedge==-1);

if length(Isrt)==1 & length(Istp)==0
    Istp=Isrt+1;
end

for i=1:length(Isrt)
    areaAprox(i)=sum(resDats(Isrt(i):Istp(i)));
end

%% fitness
function [rmse,armse,areaAprox,resDat]=fitness(dat,yhat,dp)
%
%Calculate absolute rmse, relative rmse and separate sums of adjasent positive
%residual regions
%

resDat=dat-yhat;
% Isml=dp<20e-9;
Iz=resDat<0;
resDats=medfilt1(max(resDat,0),5);
resDats(Iz)=0;
[areaAprox,Isrt,Istp]=calc_residual_areas(resDats);
Isml=dp(Istp)<20e-9;
areaAprox(Isml)=areaAprox(Isml)*3;
rmse=sqrt(mean(resDat.^2));
armse=abs(1-sum(yhat)/sum(dat));


%% dummy_iteration
function [param_new,fval,yhat]=dummy_iteration(dat,dp,nrPeaks,param)
%make some random variation, check the best
%iterate few times
%variation window decreasing after each iteration
grad=100;
h=0;
i=0.02;
while (grad>1e-12 & h<20)
    interc=1-i/2;
    slp=i;
    varMat=(rand(99,nrPeaks*2)*slp+interc);
    varMat=[ones(1,nrPeaks*2);varMat];
    %     param(:,3:4)=log(param(:,3:4));
    ch=varMat.*repmat(param,100,1);
    %     ch(:,3:4)=exp(ch(:,3:4));
    [fval3,y1_3,y3,h_ch3]=H_lognorm2pAE4(dat,dp,nrPeaks,ch,3);
    [mn,Imn3]=min(fval3);
    param=ch(Imn3,:);
    h=h+1;
    fval(h)=fval3(Imn3);
    if h>=4
        grad=mean(fval(end-3:end-1))-fval(end);
    end
    g(h)=grad;
    i=i-0.005;
end
yhat=y3(Imn3,:);
param_new=param;

% function checkDistance()

%% mlr_iteration
function [param_new,fval,yhat,yshat,nrPeaks]=mlr_iteration(dat,dp,nrPeaks,param,scaleNegErr)

%calculate

%criterias

% maxIt=60;
% sameCntLim=3;
% minErr=0.04;
% rang=0.1; %range parameter
% % scaleNegErr=2;
% doPlot=1;
% doDisp=0;

global maxIt
global sameCntLim
global minErr
global rang
global doIterPlot
global doIterDisp

bet=zeros(1,nrPeaks*2);
iter=1;
nrIt=1;
d=zeros(nrPeaks*2,1);
delta_new=1;

if nrPeaks==0
    a=1
end

%make initial population
%cube for all modes
[ch_all]=makeCube(param,rang,nrPeaks);
%cube only with new mode (last one)
[ch_new]=makeCube(param([nrPeaks,end]),rang,1);
ch_n=[repmat(param(1:nrPeaks-1),4,1),ch_new(:,1),repmat(param(nrPeaks+1:end-1),4,1),ch_new(:,2)];
%combine all
ch=[ch_all;ch_n;param];

while iter
    nrIt=nrIt+1;

    %     ch=makePopulation(param,nrPeaks,bet,49,rang);

    win=rang;
    if nrIt>2 %first iteration use initial population
        %cube for all modes
        [ch]=makeCube(param,win,nrPeaks);
        %combine all
        ch=[ch;param];
    end

    %     %set lower limit for z, experimental!!!!!!!!
    %     disp('limit')
    %     ch(:,(nrPeaks+1):end)=max(ch(:,(nrPeaks+1):end),1.7e-9);

    [lims,Idels]=checkParams(ch, nrPeaks);


    %first look if any of the modes have all values inside limits
    %if so delete that mode
    IallBad=all(Idels,1);
    %if all bad remove one completely
    %(normally mode deletion is out side the mlr-iteration, but this is exception)
    if all(IallBad)
        %         disp('all bad')
        [param,nrPeaks]=rem_mode(dat,dp,nrPeaks,param);
        [ch]=makeCube(param,win,nrPeaks);
        ch=[ch;param];
        [lims,Idels]=checkParams(ch, nrPeaks);
        IallBad=all(Idels,1);
    end

    if any(IallBad)
        ch(:,[IallBad,IallBad])=[];
        %         nrPeaks=nrPeaks-1;
        nrPeaks=sum(~IallBad);
        Idels(:,IallBad)=[];
    end

    %now delete all genes that have atleast one bed entry
    IbadGene=any(Idels,2);
    if any(IbadGene) & sum(IbadGene)~=length(IbadGene)
        ch(IbadGene,:)=[];
    end

    %check if any genes left

    if isempty(ch)
        a=1
    end

    [fval3,y1_3,y3,h_ch3]=H_lognorm2pAE4(dat,dp,nrPeaks,ch,scaleNegErr);
    %     [fval,y1,h_out]=H_lognorm2pAE2(gaD,newCh);
    %remove solutions that are negative
    Ineg=any(any(y1_3<0,3),2);

    %if all are negative start again with one mode
    if all(Ineg)
        nrPeaks=1;
        [mx,Imx]=max(dat);
        ch=[1.4,dp(Imx)];
        [fval3,y1_3,y3,h_ch3]=H_lognorm2pAE4(dat,dp,nrPeaks,ch,scaleNegErr);
        %     [fval,y1,h_out]=H_lognorm2pAE2(gaD,newCh);
        %remove solutions that are negative
        Ineg=any(any(y1_3<0,3),2);
    end
    fval3(Ineg)=[];
    y1_3(Ineg,:,:)=[];
    y3(Ineg,:)=[];
    h_ch3(Ineg,:,:)=[];
    ch(Ineg,:)=[];

    [srt Isrt]=sort(fval3);
    %     best10=ch(Isrt(1:10),:);
    fval=fval3;
    %     fvals(nrIt)=fval;
    %calculate mlr and use the terms for step direction evaluation
    %only sign counts

    %     stdNewCh=std(ch);
    %     len=size(ch,1);
    %
    %     Iok=find(stdNewCh>1e-30); %find small std
    %     % calc mlr only for variables that have std
    %     newChok=ch(:,Iok);
    %     %autoscale
    %     x=(newChok-repmat(mean(newChok),len,1))./repmat(stdNewCh(:,Iok),len,1);
    %     y_fval=(fval-min(fval))/max(fval-min(fval));
    %     %     beta = inv(x'*x)*x'*y_fval;
    %     beta = (x'*x)\x'*y_fval;
    %     betok=(beta<0)-(beta>0); %look only direction
    %     bet=ones(1,nrPeaks*2);
    %     bet(Iok)=betok;

    %     z1L=best10(1,4)-step(1)+step(1)*bet(4);
    %     z1H=best10(1,4)+step(1)+step(1)*bet(4);
    %     z2L=best10(1,5)-step(2)+step(2)*bet(5);
    %     z2H=best10(1,5)+step(2)+step(2)*bet(5);
    %     z3L=best10(1,6)-step(3)+step(3)*bet(6);
    %     z3H=best10(1,6)+step(3)+step(3)*bet(6);
    %

    [mn,Imn]=min(fval);

    param=ch(Imn,:);
    mnF(nrIt)=round(mn*1000)/1000;
    % 	(newCh(Imn,5)-newCh(Imn,4))*1e9
    %      	disp(mnF(nrIt));
    %count how many times the fval have been the same till two desimal
    if mnF(nrIt-1)==mnF(nrIt),
        sameCnt=sameCnt+1;
        rang=max(rang-0.03,0.01); %if no improvement decreas the win
    else
        sameCnt=0;
    end
    if nrPeaks==0
        a=1
    end
    %============
    %% plot (inside loop)
    %============
    if doIterPlot
        [mn,Imn]=min(fval);
        bst=ch(Imn,:);
        bestFval=fval(Imn);
        [fval,y1,y,h_ch3]=H_lognorm2pAE4(dat,dp,nrPeaks,bst,scaleNegErr);
        % 	nuc2Atk(Imn)*1e9
        % 	log10(newCh(Imn,5))-log10(newCh(Imn,4))
        lenX=length(dp);
        i=1;
        plot(dp,reshape(y1(i,:,:),nrPeaks,lenX)'),
        set(gca,'xscale','log'),
        hold on,
        plot(dp,dat,'k.'),
        plot(dp,sum(reshape(y1(i,:,:),nrPeaks,lenX),1),'k'),
        %         plot(lims,[1000,1000],'*'),
        t=[max(squeeze(y1(i,:,:))')/2]';

        [lim,Ibad]=checkParams(bst, nrPeaks);
        lims=reshape(lim,nrPeaks,2)';
        ts=repmat(t,1,2)+repmat([1:nrPeaks]',1,2);

        plot(lims,ts','*-'),

        %         plot(lims(:,1),[1000],'*'),
        %         plot(lims(:,2),[1000],'*'),
        % 		plot(initVal,ones(1,3),'r.')
        hold off
        drawnow
    end

    %stoping criteria
    %     if sameCnt==round(sameCntLim/3)
    %         rang=0.5; %increas bracheting window when some times no improvement
    %         disp('Stroke!')
    %     end

    if nrIt==maxIt | sameCnt>=sameCntLim | mn<minErr
        % 		if nrIt==maxIt  | mn<minErr
        iter=0;
        % 	nrIt
        if doIterDisp
            if nrIt==maxIt,
                disp('maxIt')
            elseif sameCnt>=sameCntLim
                disp(['noChange in, ',num2str(sameCnt),' consequent iterations; ',num2str(nrIt),' iters'])
            elseif mn<minErr,
                disp('minErr')
            end
        end
    end
end

param_new=ch(Imn,:);
fval=fval3(Imn);
yhat=y3(Imn,:);
yshat=squeeze(y1_3(Imn,:,:));

%% check parameters
function [lims,Ibad]=checkParams(param, nrPeaks)
%check if some peaks are inside the limits of other peaks
global okRangeCoef
% okRangeCoef=1.4
Idels=zeros(size(param,1),nrPeaks);
for i=1:nrPeaks
    ws=okRangeCoef*param(:,i);
    zs=param(:,nrPeaks+i);
    % plot(10.^([log10(param1(2))-log10(param1(1)),log10(param1(2))+log10(param1(1))]),[3000,3000],'*')
    limsL(:,i)=10.^(log10(zs)-log10(ws));
    limsU(:,i)=10.^(log10(zs)+log10(ws));

    Idel=param(:,nrPeaks+1:end)>repmat(limsL(:,i),1,nrPeaks) & param(:,nrPeaks+1:end)<repmat(limsU(:,i),1,nrPeaks);
    Idel(:,i)=0; %ignore the mode in question
    Idels=Idels+Idel;

    %     Idel2=param(:,nrPeaks+1:end)>lims(1) & param(:,nrPeaks+1:end)<lims(2);
    %     Idel(:,i)=0; %ignore the mode in question
end
Idels=Idels>0;

%find if peaks limits are inside other peaks limits
Iinside=zeros(size(param,1),nrPeaks);
for i=1:nrPeaks
    I=1:nrPeaks;
    Iinside(:,I~=i)=repmat(limsL(:,i),1,nrPeaks-1)<limsL(:,I~=i) & repmat(limsU(:,i),1,nrPeaks-1)>limsU(:,I~=i);
end

% Iinside=Iinside>0;

Ibad=Idels|Iinside;
lims=[limsL,limsU];

%% make inital population
function ch=makePopulation(param,nrPeaks,bet,nrMemb,range)
%
% make population by varing randomly around given parameters
%


i=ones(1,nrPeaks*2)*range; %0.2=about 10% variability
interc=1-i/2;
slp=i;
varMat=[(rand(nrMemb,nrPeaks*2).*repmat(slp,nrMemb,1)+repmat(interc,nrMemb,1));ones(1,nrPeaks*2)];
varMat(:,bet==-1)=max(varMat(:,bet==-1),1); %Only positive window
varMat(:,bet==1)=min(varMat(:,bet==1),1); %Only negative window
ch=repmat(param,nrMemb+1,1).*varMat;

%% rem_mode
function [param,nrPeaks]=rem_mode(dat,dp,nrPeaks,param)
%remove modes that are bad
%

global smlLimR
global smlLimA
global errLim

%criterias
% smlLimR=0.05;
% smlLimA=200;
% errLim=0.3;


%out of limit

%     Check if final best solution is acceptable
[lims,Ibad]=checkParams(param, nrPeaks);
%if not remove the bad mode
[fval,y1,y,h_ch3]=H_lognorm2pAE4(dat,dp,nrPeaks,param,3);


%too small
%sum of mode is less than smlLimR of sum of data
squeezedY=squeeze(y1(1,:,:));
[n,m]=size(squeezedY);
if n>m
    squeezedY=squeezedY';
end
msum=sum(squeezedY,2)';
dsum=sum(dat)';

%relative peak area approximation
IsmlR=msum<dsum*smlLimR;

%absolute
IsmlA=msum<smlLimA;



%error of the mode to big

res=abs(dat-y);
err=zeros(1,nrPeaks);
for i=1:nrPeaks
    Ihgh=dp>lims(1,i) & dp<lims(1,i+nrPeaks);
    Iz=dat(Ihgh)==0;
    res_tmp=res(Ihgh);
    dat_tmp=dat(Ihgh);
    res_tmp(Iz)=[];
    dat_tmp(Iz)=[];
    if isempty(res_tmp)
        err(i)=0;
    else
        err(i)=mean(res_tmp./dat_tmp);
    end
end

Ierr=err>errLim;

% Idel=Ibad|(IsmlR&IsmlA)|Ierr;
Idel=Ibad|(IsmlR)|Ierr;


%if all are bad keep the bigger mode
if all(Idel)
    [mx Imx]=max(max(squeeze(y1)'));
    Idel(Imx)=0;
end
if any(Idel)
    param(:,[Idel,Idel])=[];
    %     bet(:,[Idel,Idel])=[];
    nrPeaks=sum(~Idel);
    lims(:,[Idel,Idel])=[];
end
% [fval,y1,y,h_ch3]=H_lognorm2pAE4(dat,dp,nrPeaks,param,3);

%% make cube
function [ch]=makeCube(param,win,nrPeaks)
%
% make window around parameters
%win-relative window,  same size as param


eval(['global HLmat',num2str(nrPeaks)])
eval(['HLmat=','HLmat',num2str(nrPeaks),';'])


step=param.*win;
lowLims=param-step;
hghLims=param+step;

lims=[lowLims;hghLims];
if isempty(HLmat)
    %make matrix with all permutations of high an low boundaries
    startFor=['h=0;'];
    statement=['h=h+1;'];
    endFor=[];

    for i=1:nrPeaks*2
        startFor=[startFor,'for i',num2str(i),'=1:2,'];
        statement=[statement,'HLmat(h,',num2str(i),')=i',num2str(i),';'];
        endFor=[endFor,'end,'];
    end

    eval([startFor,statement,endFor])

    eval(['HLmat',num2str(nrPeaks),'=HLmat;'])
end

%make final matrix
for i=1:nrPeaks*2
    ch(:,i)=lims(HLmat(:,i),i);
end



%% nlinOptim
function [par1]=nlinOptim(dat,dp,param)

% NOT FINISHED

%     limLow=[ones(1,nrPeaks)*1.1,param(nrPeaks+1:end)*1.01];
%     limHgh=[ones(1,nrPeaks)*1.5,param(nrPeaks+1:end)*0.99];


xy=[dp',dat'];

Isel=all(~isnan(xy),2);
xy=xy(Isel,:);
par1=param;

%optimize
opt = optimset('GradObj','off','Display','off','Diagnostics','off','TolFun',1,'LargeScale','off');
%     [par1]=fmincon(@(par1) H_lognorm2p(par1,xy),init,[],[],[],[],limLow,limHgh,@H_lognormCon,opt);
% [par1]=nlinfit(xy,xy(:,2),@H_lognorm2p,param);
[par1]=fminsearch(@(par1) H_lognorm2p(par1,xy),param);