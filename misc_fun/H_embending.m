function emdata=H_embending(data,step)
%
%
%
%emdata=H_embending(data,step)
%

[n m]=size(data);
ste=step-1;
if m~=1 && n~=1
    error('Embending can only be done to single variable');
end
data=data(:);
[n m]=size(data);
emb_dim=step;

%Find the number of rows to be transformed
temp_n=emb_dim*round(n/emb_dim);
if temp_n>=n,
   temp_n=emb_dim*(round(n/emb_dim)-1);
end

% h=0;
% for j=1:emb_dim,
%    h=h+1;
%    
%    emb_data(:,j)=data(j:(temp_n-1+h),1);  %transform the time series to matrix [1 2 3;4 5 6; etc.]
% end


%c=a(ones(20,1),:)+ad(:,ones(1,8));
a=1:step;
ad=(1:n)';

inx=a(ones(n,1),:)+ad(:,ones(1,step))-1;
ir=find(inx(:,end)==n);
inx(ir:end,:)=[];
emb_data=data(inx);
emdata=emb_data;
