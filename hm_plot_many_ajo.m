
d=datenum(2004,4,15);
d2=datenum(2004,4,15);

parPath1='C:\1h\dat\parameters\hm_param.txt';
parPath2='C:\1h\dat\parameters\hm_param_varrio.txt';


hmD1=hm_load([d,d2],'parPath',parPath1,'dmps');
hmD2=hm_load([d,d2],'parPath',parPath2,'dmps');

H_subplot(2,1,1)
hm_plot(hmD1,'dmps')
set(gca,'xticklabel',' ')
H_subplot(2,1,2)
hm_plot(hmD2,'dmps')
title(' ')

axes('position',[0 0 1 1],'visible','off')
text(0.15,0.91,hmD1.meta.dmps.paths(11:end-6),'fontsize',12,'fontweight','bold')
text(0.15,0.48,hmD2.meta.dmps.paths(11:end-6),'fontsize',12,'fontweight','bold')