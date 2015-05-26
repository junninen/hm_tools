function [fval,y1,y,h_ch]=H_lognorm2pAE4(dat,dp,nrPeaks,param,negWeight)
%
%Modified H_lognorm2pAE3 to run without structures, just raw imputs
%Log-normal distribution with 2 param
%height of the distribution is calculated analytically
%optimized for aerosol data, 
%
%calculates all chromosomes at the same time not in for loop
%
%Initialization is made with H_lognorm2pAE1_init.m
%
% h - height
% w - width
% z - position
%

% xx=zeros(x.xlng,nr_peaks);

if nargin<=4
    negWeight=1;
end
mpl=ones(length(dat),1);
x=dp;
lengthX=length(x);
nrOfChroms=size(param,1);
% fval=zeros(nrOfChroms,1);

x1=repmat(x,nrPeaks,1)';
x2=repmat(x1,1,nrOfChroms);
xx=reshape(x2',nrOfChroms,nrPeaks,lengthX);


%Set parameters
w  = param(:,1:nrPeaks,ones(1,lengthX));%               w  = w(mpl,:);
z  = param(:,nrPeaks+1:nrPeaks*2,ones(1,lengthX));%    z  = z(mpl,:);

%calculate h

% sqrt(2*pi)=2.5066
% (3*pi)^(1/3)=2.1123
logW=log(w);

c1=1./(logW.*2.5066);
% c2=log(xx./z).^2;
c2=log((xx./z)).^2;
c3=2.*(logW.^2);

b=c1.*exp(-c2./c3);
%matrix algebra can be done only for 2-D matrix :(
h=ones(nrOfChroms,nrPeaks,lengthX);

%calculate at each value of size distribution
h_ch=zeros(nrOfChroms,nrPeaks);
for i=1:nrOfChroms,
    b_tmp=b(i,:,:);
    b_tmp(isnan(b_tmp))=0;
    b_ch=reshape(b_tmp,nrPeaks,lengthX)'; %B of one chromosome
        h_out=max(dat*pinv(b_ch'),0); %original
%     h_out=max(dat(i,:)*pinv(b_ch'),0); %original
%     h_out=max(str.dat(i,:)*pinv(b_ch'),0); %original
    h_ch(i,:)= h_out;
      h(i,:,:)  = h_out(mpl,:)';
end

%calculate
%all chromosomes at the same time

% logW=log(w);
c1=h./(logW.*2.5066);
% c2=log(xx./z); %already calculated
% c3=2.*(logW.^2); %already calculated

% y1=[c1.*exp(-0.5*(c2./logW))];

y1=c1.*exp(-c2./c3);

% ytemp=y1(:,1,:)+y1(:,2,:)+y1(:,3,:)+y1(:,4,:);
In=isnan(y1);
y1(In)=0;
ytemp=sum(y1,2);

%ytemp=nansum(y1);
y=reshape(ytemp,nrOfChroms,lengthX);

% 				case 'rmse'
L=repmat(lengthX,nrOfChroms,1);
orig=repmat(dat,nrOfChroms,1);
% m=max(orig(1,:));

%  wgt=str.nuclx+1; %will weight nuclaetin area (3-20) nm
%  wgt=repmat(wgt,nrOfChroms,1);
   wgt=ones(size(orig)); %weight with negative error
      In=find(y-orig>0);
      wgt(In)=wgt(In)*negWeight;
%  fval=sqrt(sum(((orig/m-y/m).*wgt).^2,2)./L);
% fval=sqrt(sum(((orig/m-y/m)).^2,2)./L);

%  fval=(sum(((orig/m-y/m).*wgt).^2,2)./L).^(1/4);

%scaling each residual with original concentration
%   fval=sqrt(sum((((orig-y)./orig).*wgt).^2,2)./L);

% fval=sqrt(sum(((orig-y)).^2,2)./L);
div=orig-y;
  fval=sqrt(sum((((div)).*wgt).^2,2)./L);