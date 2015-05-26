function [hmD,gr_str]=hm_gr_max_auto(hmD,sizeLim, concLim,timeLim)

% timeLim=[9 17];
% sizeLim=[3 30]*1e-9;
% concLim=1000;
tm=hmD.meta.dmps.tim{1};
dp=hmD.meta.dmps.dp{1};
dat=hmD.dmps{1}(2:end,3:end);
dv=datevec(tm);
Itm=dv(:,4)>timeLim(1) & dv(:,4)<timeLim(2);
Idp=dp>sizeLim(1) & dp<sizeLim(2);
datF=NaN(size(dat));
datF(Itm,Idp)=dat(Itm,Idp);
datF=max(datF,concLim);
datF(datF==concLim)=NaN;
figure(111)
% pcolor(tm,dp,datF'),set(gca,'yscale','log')
pcolor(tm,dp,log10(max(dat,1))'),
set(gca,'yscale','log')
caxis([1 4])

[mx2,Imx2]=max(datF,[],2);

dpSel=dp(Imx2)';
Idel=Imx2==1 | dpSel==sizeLim(2);
dpSel(Idel)=[];
tmSel=tm;
tmSel(Idel)=[];
hold on
plot(tmSel,dpSel,'ow')
%bootstrap
x=(tmSel-floor(tmSel))*24; %conversion to hours
y=dpSel;
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

int_avr=mean(p_bs(:,2));
int_std=std(p_bs(:,1));

plot(tmSel,x*gr_avr+int_avr,'-w')
hold off
gr_str.gr_avr=gr_avr;
gr_str.int_avr=int_avr;

hmD.fit.dmps.gr_str=gr_str;

