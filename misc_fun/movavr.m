function y=movavr(x,m)
%y=movavr(x,m)
%moving average
%  y = smoothed
%  x = orig. data
%  m = window
%
% missing data (NaNs) are linearly interpolated during the calculation
% and NaN are put back afterwards.


x=x(:);

%find all nans
Inan=isnan(x);

%if first or last value is NaN replase it with first available value
%this is only for interpolation and will be replaced by NaN eventually
if isnan(x(1))
    Inn=find(~isnan(x));
    x(1)=x(Inn(1));
end

if isnan(x(end))
    Inn=find(~isnan(x));
    x(end)=x(Inn(end));
end

%replace NaNs with linear interpolation
if any(isnan(x))
    x=MDrepl1(x);
end

x=x(:)';
z = [0 cumsum(x)];
n=length(x);


%emb=som_normalize(std(H_embending(z',m)')','range');

yt = ( z(m+1:n+1) - z(1:n-m+1) ) / m;
y=zeros(1,n);
y(:)=nan;
y(1)=mean(x(1:m));
y(n)=mean(x(n-m:n));

me=round((length(x)-length(yt))/2);
ne=length(yt);
y(me+1:me+ne)=yt;

%place in the NaNs
y(Inan)=NaN;
