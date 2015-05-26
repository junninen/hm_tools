
function out=H_rm_spikes(dat,step)
%
% Clean vector (column if matrix) from spiks, by calculating moving median over steps sised window 
%
% out=H_rm_spikes(dat,step)
%


%Heikki Junninen
% Sept 2010

[n m]=size(dat);
tr=0;
if n==1
    dat=dat';
    [n m]=size(dat);
    tr=1;
end
out=dat;

if nargin==1
    step=3;
end
for i=1:m
    y=dat(:,i);
%     y_mid=y(2:n-1);
%     y_neig_avr=(y(1:n-2)+y(3:n))/2;
%     y_neig_avr=median([y(1:n-2),y(3:n)],2);
%     Ibad=abs(y_mid-y_neig_avr)./abs(y_neig_avr)>lim;
%     if any(Ibad)
%         y_mid(Ibad)=y_neig_avr(Ibad);
%     end
%     out(2:n-1,i)=median([y(1:n-2),y(2:n-1),y(3:n)],2);
    
    emdata=H_embending(y,step);
    tmp=median(emdata,2);
     [nt mt]=size(tmp);
    halfStep=round(step/2);
    out(halfStep:nt+halfStep-1,i)=tmp;
    out(1:halfStep-1,i)=tmp(1);
    out(nt+halfStep:n,i)=tmp(end);
end
if tr
    out=out';
end