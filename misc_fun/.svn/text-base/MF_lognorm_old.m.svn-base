function [param,fval,y1_F,h_out,N]=MF_lognorm(dat,dp,doPlot,mxNrM,init)
%
% fit lognormal distribution to size distribution
%
%[param,fval,yhat]=MF_lognorm(dat,dp,[doPlot],[initPar])
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

if ~doinit
    nrPeaks=length(init)/2;
    ch=init;
    [ch,fval_new,yhat,yshat,nrPeaks]=mlr_iteration(dat,dp,nrPeaks,init);
    Ineg=find(any(yshat<0,2));
    if ~isempty(Ineg)
        %         init([Ineg,Ineg+nrPeaks])=[];
        %         Ineg=[];
        %         ch=init;
        %         nrPeaks=nrPeaks-1;
        doinit=1;
    end
end

if doinit
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

end

[fval,y1,y,h_ch]=H_lognorm2pAE4(dat,dp,nrPeaks,ch,3);
[mn,Imn]=min(fval);
yhat=y(Imn,:);
fval=fval(Imn);
param1=ch(Imn,:);

resDat=dat-yhat;

% if positive residual is more than 10% of original data fit new mode
%below 1nm dont add new mode
param_new=param1;
fval_new=fval;
Iabove1nm=find(dp>1e-9); %index for values above 1nm
% notTooClose=1; %new mode have to be further than 1 sigma
% plot(10.^([log10(param1(2))-log10(param1(1)),log10(param1(2))+log10(param1(1))]),[3000,3000],'*')


% while sum(max(resDat(Iabove1nm),0))>sum(dat(Iabove1nm))*0.08 & nrPeaks<mxNrM
do=1;
w=0;
rmse=sqrt(mean(resDat.^2));
armse=1-sum(yhat)/sum(dat);

while (rmse>200 | armse>0.1) & nrPeaks<mxNrM & do
    w=w+1;
    nrP=nrPeaks;
    [param_new,nrPeaks]=add_mode(dat,resDat,dp,nrPeaks,ch,param_new);
    %     [param_new,fval_new,yhat]=dummy_iteration(dat,dp,nrPeaks,param_new);
    %if add_mode-function did not add a mode, means there is no place to
    %add

    [param_new,fval_new,yhat,yshat,nrPeaks]=mlr_iteration(dat,dp,nrPeaks,param_new);
    resDat=dat-yhat;
    rmse(w)=sqrt(mean(resDat.^2));

    %stop criterions

    rmse=sqrt(mean(resDat.^2));
    armse=1-sum(yhat)/sum(dat);

    if nrP==nrPeaks
        do=0;
    end
    if w>5
        do=0;
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
[fvalF,y1_F,yF,h_chF]=H_lognorm2pAE4(dat,dp,length(param)/2,param,4);

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
%SUBFUNCTIONS

function [param_new,nrPeaks]=add_mode(dat,resDat,dp,nrPeaks,ch,param)

%fit just positive residuals
if 0
    %%find solution from provided ch
    %find the solutions only one sigma away from any existing modes
    for i=1:nrPeaks
        ws=param(i);
        zs=param(nrPeaks+i);
        % plot(10.^([log10(param1(2))-log10(param1(1)),log10(param1(2))+log10(param1(1))]),[3000,3000],'*')
        lims(i,:)=10.^([log10(zs)-log10(ws),log10(zs)+log10(ws)]);
        Idel=ch(:,2)>lims(1) & ch(:,2)<lims(2);
        ch(Idel,:)=[];
    end
    [fval2,y1_2,y2,h_ch2]=H_lognorm2pAE4(max(resDat,0),dp,1,ch,10);
    [mn,Imn2]=min(fval2);

    param2=ch(Imn2,:);
    yhat=y2(Imn2,:);
    y1hat=y2(Imn2,:,:);
    fval=fval2(Imn2);

    param_new=[param(1:nrPeaks),param2(1),param(nrPeaks+1:end),param2(2)];
    nrPeaks=nrPeaks+1;
end

if 1
    %just find where is smoothed maximum in residual
    resDat=medfilt1(max(resDat,0),5); %only positive residuals
    resDats=resDat(1:end-2)+resDat(2:end-1)+resDat(3:end);
    resDats=[resDat(1),resDats,resDat(end)];

    %     plot(resDat),hold all,plot(resDats),
    %find the solutions only one sigma away from any existing modes
    for i=1:nrPeaks
        ws=param(i);
        zs=param(nrPeaks+i);
        % plot(10.^([log10(param1(2))-log10(param1(1)),log10(param1(2))+log10(param1(1))]),[3000,3000],'*')
        lims(i,:)=10.^([log10(zs)-log10(ws),log10(zs)+log10(ws)]);
        Idel=dp>lims(1) & dp<lims(2);
        resDats(Idel)=0;
        %         plot(resDats),
    end

    [mx,Imx]=max(resDats);
    if mx>100
        param_new=[ones(1,nrPeaks)*1.4,1.4,param(nrPeaks+1:end),dp(Imx)];
        nrPeaks=nrPeaks+1;
    else
        param_new=param;
    end
end

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
function [param_new,fval,yhat,yshat,nrPeaks]=mlr_iteration(dat,dp,nrPeaks,param)

%calculate

bet=zeros(1,nrPeaks*2);
iter=1;
nrIt=1;
while iter
    nrIt=nrIt+1;
    i=ones(1,nrPeaks*2)*0.5; %0.2=about 10% variability
    interc=1-i/2;
    slp=i;
    varMat=[(rand(49,nrPeaks*2).*repmat(slp,49,1)+repmat(interc,49,1));ones(1,nrPeaks*2)];
    varMat(:,bet==-1)=max(varMat(:,bet==-1),1); %Only positive window
    varMat(:,bet==1)=min(varMat(:,bet==1),1); %Only negative window
    ch=repmat(param,50,1).*varMat;

    [lims,Idels]=checkParams(ch, nrPeaks);


    %first look if any of the modes have all values inside limits
    %if so delete that mode
    IallBad=all(Idels,1);
    if any(IallBad)
        ch(:,[IallBad,IallBad])=[];
        %         nrPeaks=nrPeaks-1;
        nrPeaks=sum(~Ibad);
        Idels(:,IallBad)=[];
    end

    %now delete all genes that have atleast one bed entry
    IbadGene=any(Idels,2);
    if any(IbadGene)
        ch(IbadGene,:)=[];
    end

    %check if any genes left


    [fval3,y1_3,y3,h_ch3]=H_lognorm2pAE4(dat,dp,nrPeaks,ch,3);
    %     [fval,y1,h_out]=H_lognorm2pAE2(gaD,newCh);
    %remove solutions that are negative
    Ineg=any(any(y1_3<0,3),2);
    fval3(Ineg)=[];
    y1_3(Ineg,:,:)=[];
    y3(Ineg,:)=[];
    h_ch3(Ineg,:,:)=[];
    ch(Ineg,:)=[];

    [srt Isrt]=sort(fval3);
    %     best10=ch(Isrt(1:10),:);
    fval=fval3;

    %calculate mlr and use the terms for step direction evaluation
    %only sign counts

    stdNewCh=std(ch);
    len=size(ch,1);

    Iok=find(stdNewCh>1e-30); %find small std
    % calc mlr only for variables that have std
    newChok=ch(:,Iok);
    x=(newChok-repmat(mean(newChok),len,1))./repmat(stdNewCh(:,Iok),len,1);
    y_fval=(fval-min(fval))/max(fval-min(fval));
    beta = inv(x'*x)*x'*y_fval;
    betok=(beta<0)-(beta>0); %look only direction
    bet=ones(1,nrPeaks*2);
    bet(Iok)=betok;


    [mn,Imn]=min(fval);
    param=ch(Imn,:);
    % 	(newCh(Imn,5)-newCh(Imn,4))*1e9
    mnF(nrIt)=round(mn*1000)/1000;
    % 	disp(mnF(nrIt));
    %count how many times the fval have been the same till two desimal
    if mnF(nrIt-1)==mnF(nrIt),
        sameCnt=sameCnt+1;
    else
        sameCnt=0;
    end

    %     Check if final best solution is acceptable
    [lims,Ibad]=checkParams(param, nrPeaks);
    %if not remove the bad mode
    if all(Ibad)
        [fval,y1,y,h_ch3]=H_lognorm2pAE4(dat,dp,nrPeaks,param,3);
        [mx Imx]=max(max(squeeze(y1)'));
        %keep the bigger mode
        Ibad(Imx)=0;
    end
    if any(Ibad)
        param(:,[Ibad,Ibad])=[];
        bet(:,[Ibad,Ibad])=[];
        nrPeaks=sum(~Ibad);
        lims(Ibad,:)=[];
    end
    [fval,y1,y,h_ch3]=H_lognorm2pAE4(dat,dp,nrPeaks,param,3);


    %============
    %plot (inside loop)
    %============
    if 0
        %         [mn,Imn]=min(fval);
        %         bst=ch(Imn,:);
        %         bestFval=fval(Imn);
        %         [fval,y1,y,h_ch3]=H_lognorm2pAE4(dat,dp,nrPeaks,bst,3);
        % 	nuc2Atk(Imn)*1e9
        % 	log10(newCh(Imn,5))-log10(newCh(Imn,4))
        lenX=length(dp);
        i=1;
        plot(dp,reshape(y1(i,:,:),nrPeaks,lenX)'),set(gca,'xscale','log'),
        hold on,
        plot(dp,dat,'k.'),
        plot(dp,sum(reshape(y1(i,:,:),nrPeaks,lenX),1),'k'),
        %         plot(lims,[1000,1000],'*'),
        t=[max(squeeze(y1(i,:,:))')/2]';

        ts=repmat(t,1,size(lims,2))+repmat([1:size(lims,1)]',1,2);
        plot(lims',ts','*-'),
        %         plot(lims(:,1),[1000],'*'),
        %         plot(lims(:,2),[1000],'*'),
        % 		plot(initVal,ones(1,3),'r.')
        hold off
        drawnow
    end

    %stoping criteria
    maxIt=45;
    sameCntLim=10;
    minErr=0.04;
    %     if sameCnt==round(sameCntLim/3)
    %         win=win+16; %increas bracheting window when some times no improvement
    %         % 		disp('Stroke!')
    %     end

    doDisp=0;
    if nrIt==maxIt | sameCnt>=sameCntLim | mn<minErr
        % 		if nrIt==maxIt  | mn<minErr
        iter=0;
        % 	nrIt
        if doDisp
            if nrIt==maxIt,
                disp('maxIt')
            elseif sameCnt>=sameCntLim
                disp('noChange')
            elseif mn<minErr,
                disp('minErr')
            end
        end
    end
end

% param_new=ch(Imn,:);
% fval=fval3(Imn);
% yhat=y3(Imn,:);
% yshat=squeeze(y1_3(Imn,:,:));

%     Check if final best solution is acceptable
[lims,Ibad]=checkParams(param, nrPeaks);
%if not remove the bad mode
if any(Ibad)
    param(:,[Ibad,Ibad])=[];
    bet(:,[Ibad,Ibad])=[];
    nrPeaks=sum(~Ibad);
end
[fval,y1,y,h_ch3]=H_lognorm2pAE4(dat,dp,nrPeaks,param,3);

param_new=param;
yhat=y;
yshat=squeeze(y1(1,:,:));

%% check parameters
function [lims,Idels]=checkParams(param, nrPeaks)
Idels=zeros(size(param,1),nrPeaks);
for i=1:nrPeaks
    ws=1.2*param(i);
    zs=param(nrPeaks+i);
    % plot(10.^([log10(param1(2))-log10(param1(1)),log10(param1(2))+log10(param1(1))]),[3000,3000],'*')
    lims(i,:)=10.^([log10(zs)-log10(ws),log10(zs)+log10(ws)]);

    Idel=param(:,nrPeaks+1:end)>lims(i,1) & param(:,nrPeaks+1:end)<lims(i,2);
    Idel(:,i)=0; %ignore the mode in question
    Idels=Idels+Idel;

    %     Idel2=param(:,nrPeaks+1:end)>lims(1) & param(:,nrPeaks+1:end)<lims(2);
    %     Idel(:,i)=0; %ignore the mode in question
end
Idels=Idels>0;