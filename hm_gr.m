function [hmD,gr_str]=hm_gr(hmD,srs,a)
%
%Calculate growth rate for fitted modes along the
%virtual selected line
%
%[hmD,gr_str]=hm_gr(hmD,srs,[a])
%
% hmD - hm_tools structure
% srs - source, 'dmps'
% a   - endpoints of the line, if not given have to selec manually


if ~isfield(hmD,'fit')
    disp(['hm_gr: Fitting have to be made first'])
    gr_str=[];
    return
elseif ~isfield(hmD.fit,srs)
    error(['hm_gr: Field ',srs,' does not exist!'])
end

nrDays=eval(['length(hmD.',srs,')']);

for d=1:nrDays
    %make plot, one day at the time
    if nargin==2
    a=ginput(2);
    end

%     tmBig=eval(['hmD.fit.',srs,'.tm']);
    tmBig=eval(['hmD.meta.',srs,'.tim{d}']);
    %convert time to hours starting
    mnH=min(floor(tmBig*24));
    tmBig=(tmBig*24)-mnH;
    a(:,1)=a(:,1)*24-mnH;
    zBig=eval(['hmD.fit.',srs,'.zs{d}']);
    wBig=eval(['hmD.fit.',srs,'.ws{d}']);
    nBig=eval(['hmD.fit.',srs,'.Ns{d}']);
    % dis=eval(['hmD.fit.',srs,'.yhat']);

    nrM=size(zBig,2);

    [p,s]=polyfit(a(:,1),a(:,2),1);
    ItmWin=tmBig>=a(1,1)&tmBig<=a(2,1);
    tm=tmBig(ItmWin);
    z=zBig(ItmWin,:);
    w=wBig(ItmWin,:);
    n=nBig(ItmWin,:);
    % dists=dis(ItmWin)';
    tmHuge=repmat(tm,1,nrM);


    nrTm=length(tm); %how many time points used for calc
    zhat=p(1).*tm+p(2);

    % hold on
    % plot(tm,zhat,'w')
    % hold off

    er=abs(log10(repmat(zhat,1,nrM))-log10(z));

    %select smallest by min-fun indexes

    [mn,Imn]=min(er');
    Itm=[1:nrTm]';
    Iok=mn<0.2;
    Itm=Itm(Iok);
    Imn=Imn(Iok);
    Isel=size(tm,1)*(Imn-1)+Itm';

    tmSel=tm(Iok);

    %growth rate
    zSel=z(Isel);

    %formation rate
    nSel=n(Isel);


    % hold on
    % plot(tmSel,zSel,'w*')
    % hold off

    %final fit gr
    % x=tmSel*24;
    x=tmSel; %conversion to hours made already
    y=zSel;
    x=x(:);
    y=y(:);

    %bootstrap
    lx=length(x);
    indx=1:lx;
    nrbs=round(lx*5/6);
    for i=1:100
        randIndx=randperm(lx);
        Ibs=randIndx(1:nrbs);
        [p_bs(i,:),s_bs{i}]=polyfit(x(Ibs),y(Ibs),1);
    end
    gr_prc25=prctile(p_bs(:,1),25);
    gr_prc50=prctile(p_bs(:,1),50);
    gr_prc75=prctile(p_bs(:,1),75);

    gr_avr=mean(p_bs(:,1));
    gr_std=std(p_bs(:,1));

    int_prc25=prctile(p_bs(:,2),25);
    int_prc50=prctile(p_bs(:,2),50);
    int_prc75=prctile(p_bs(:,2),75);

    int_avr=mean(p_bs(:,2));
    int_std=std(p_bs(:,2));


    yhat=x.*mean(p_bs(:,1))+mean(p_bs(:,2));
    % hold on
    % plot(tmSel,yhat,'b','linewidth',2)
    % hold off
    % title([srs,' ',hmD.FNam,' gr\_avr=',num2str(round(gr_avr*1e9*100)/100),' [nm/h]'])

    % %bootstrap formation rate
    % x=tmSel*24;
    x=tmSel;
    y=nSel;
    x=x(:);
    y=y(:);

    for i=1:100
        randIndx=randperm(lx);
        Ibs=randIndx(1:nrbs);
        [p_bs(i,:),s_bs{i}]=polyfit(x(Ibs),y(Ibs),1);
    end
    fr_prc25=prctile(p_bs(:,1),25);
    fr_prc50=prctile(p_bs(:,1),50);
    fr_prc75=prctile(p_bs(:,1),75);

    fr_avr=mean(p_bs(:,1));
    fr_std=std(p_bs(:,1));

    %save results
    gr_str(d).dp_range=[min(zSel),max(zSel)];
    gr_str(d).tm_range=[(x(1)+mnH)/24, (x(end)+mnH)/24];
    gr_str(d).gr_avr=[gr_avr,gr_std];
    gr_str(d).gr_med=[gr_prc50,gr_prc25,gr_prc75];
    gr_str(d).int_avr=[int_avr,int_std];
    gr_str(d).int_med=[int_prc50,int_prc25,int_prc75];
    % gr_str.totConc=totConc;
    gr_str(d).fr_avr=[fr_avr,fr_std];
    gr_str(d).fr_med=[fr_prc50,fr_prc25,fr_prc75];
    gr_str(d).tm=(tmSel(:)+mnH)/24;
    gr_str(d).z=zSel(:);
    gr_str(d).N=nSel(:);
    gr_str(d).zhat=yhat;
end

% tmNow=datestr(now,30);
% eval(['hmD.fit.',srs,'.gr','_',tmNow])=gr_str;
%
% if isfield(eval(['hmD.fit.',srs,'.gr']))
%     fnams=fieldnames(eval(['hmD.fit.',srs,'.gr']));
%     nrStr=length(fnams);
% else
%     eval(['hmD.fit.',srs,'.gr.gr_str0']);
% end

