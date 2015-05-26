function mmat=H_2dmedfilt(mat,win)
%2d median filter, win must be odd
%
% H_2dmedfilt(mat,win)
%   mat-data matrix
%   win-vector with two elements window in cloumn wise and window in row
%       wise. Both must be odd!
%

%
%

[n,m]=size(mat);

bigTemp=NaN(n,m,prod(win));

wx1=-(win(1)-1)/2;
wx2=(win(1)-1)/2;

wy1=-(win(2)-1)/2;
wy2=(win(2)-1)/2;

ix=1:m;
iy=1:n;

h=0;
for wx=wx1:wx2
    for wy=wy1:wy2
        h=h+1;
        Ix=(min(max(iy+wy,1),n));
        Iy=(min(max(ix+wx,1),m));
        bigTemp(:,:,h)= mat(Ix,Iy);
    end
end

mmat=median(bigTemp,3);