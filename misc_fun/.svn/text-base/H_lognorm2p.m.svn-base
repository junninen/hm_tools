function [fval,ytot,y1,h_out]=H_lognorm2p(param,xy)
%
%Log-normal distribution with 2 param
%height of the distribution is calculated analytically from y
%Only one set of parameters given at the time
%Can be used with nlinfit-routine
%
% [fval,ytot,y1,h_out]=H_lognorm2p(param,xy)
%
% xy - x values in 1. column and y val in 2.
%
% w - width
% z - position
%

nr_peaks=length(param)/2;

param=real(param);

 x  = repmat(xy(:,1)',nr_peaks,1);
y  = repmat(xy(:,2)',nr_peaks,1);
m=length(x);
w  = param(1:nr_peaks);              w=repmat(w',1,m);
z  = param(nr_peaks+1:nr_peaks*2);   z=repmat(z',1,m);
%z=real(z);

% check if param values in side limits
% Iwout=find(w>2.2| w<1);
% Izout=find(z>x(end,end) | z<x(1));
% if any(Iwout) | any(Izout)
%    ytot=ones(size(xy(:,2)))*1e10;
%    return
% end
% 
%calculate h
% find closest known concentrations at z value
[dummy ind]=min(abs(x-z)'); %index of closest z
yz=y(1,ind);
Idif=find(dummy~=0);
% maybe better to interpolate the z value if not exact

for i=1:nr_peaks
    Is=find(x(i,:)<z(i,:)); %index of elements of smaller x tnan z
    Ib=find(x(i,:)>z(i,:)); %index of elements of  bigger x tnan z
    if ~isempty(Is) & ~isempty(Ib)
%         Isml(i)=max(Is);
%         Ibgr(i)=min(find(x(i,:)>z(i,:)));
%         y_vana=[y(i,Isml(i)),y(i,Ibgr(i))];
%         x_vana=[x(i,Isml(i)),x(i,Ibgr(i))];
%         x_uus=[x(i,Isml(i)),z(i,1),x(i,Ibgr(i))];
%         y_uus=interp1(x_vana,y_vana,x_uus);
%         y_interp(i)=(y_uus(2));

        
        Isml=max(Is);
        Ibgr=min(find(x(i,:)>z(i,:)));
        k=(y(i,Ibgr)-y(i,Isml))/(x(i,Ibgr)-x(i,Isml));
        b=y(i,Isml)-k*x(i,Isml);
        y_interp(i)=k*z(i,1)+b;
    elseif isempty(Is)
        y_interp(i)=y(1);
        %check what should be done when z smaller than x
    elseif isempty(Ib)
        y_interp(i)=y(end);

    end
end

%lognorm distribution
yz(Idif)=y_interp(Idif);
lgW=log(w);
c1=1./(lgW.*2.5066);
c2=log(x./z).^2;
c3=2.*(lgW.^2);

b=c1.*exp(-c2./c3);
h_out=max(yz*pinv(b(:,ind)'),0)';

h  = repmat(h_out,1,m);


%calculate
c1=h./(lgW.*2.5066);
% c2=log(x./z);

y1=c1.*exp(-(c2./c3));
ytot=sum(y1,1);

orig=xy(:,2);
ytot=ytot';

L=length(ytot);
%RMSE
dif=orig-ytot;
if 0
    fval=sqrt(sum(dif.^2)./L);
else
    %weighter RMSE
    wgt=ones(size(orig)); %weight with negative error    
    In=find(dif<0);
    wgt(In)=wgt(In)*2;
    fval=sqrt(sum((dif.*wgt).^2)./L);
end


