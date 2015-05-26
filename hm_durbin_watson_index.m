function dw=hm_durbin_watson_index(dat,dats)
%
%Calculate Durbin-Watson index for automatic smoothing window evaluation
%
%  dw=hm_durbin_watson_index(dat,dats)
%
%    dw - Durbin-Watson index
%   dat - raw data
%  dats - smoothed data
%
%if input is matrix dw will be average of row and col-wise index
%
%Ref: J. Durbin and G.S. Watson, Biometrika 37 (1950) 409
%G. Vivo-Truyols et al. (2005) Automatic program for peak detection and
%deconvolution of multi-overlapped chromatographic signals. Part I: Peak
%detection. Journal of Chromatography A, 1096, 133-145

%Heikki Junninen
%Update 20.Nov.2007
%  - min value for sum(res) set to 1e-18, to awoid warning

if ~isequal(size(dat),size(dats))
    error('hm_durbin_watson_index: data and smoothed data must be same size')
end

[n,m]=size(dat);

%replace NaNs by 0
Idn=isnan(dat);
Isn=isnan(dats);
dat(Idn)=0;
dat(Isn)=0;
dats(Idn)=0;
dats(Isn)=0;

if n==1 | m==1
    n=length(dat);
    res=dat-dats;
    dw1=sum(((res(2:end))-(res(1:end-1))).^2);
    dw2=sum(res.^2);

    dw3=(n/(n-1));
    dw=dw1./dw2.*dw3;
else
    res=dat-dats;
    dw1=sum(((res(2:end,:))-(res(1:end-1,:))).^2,1);
    dw2=sum(res.^2,1);

    dw3=(n/(n-1));
    dwR=mean(dw1./dw2.*dw3);

    dw1=sum(((res(:,2:end))-(res(:,1:end-1,:))).^2,2);
    dw2=max(sum(res.^2,2),1e-18);
    dw3=(m/(m-1));
    dwC=mean(dw1./dw2.*dw3);
    dw=mean([dwR,dwC]);
end