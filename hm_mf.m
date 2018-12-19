function hmD=hm_mf(hmD,srs,meth,maxNrM,doPlot)
%
%mode fitting with MF_lognorm
%
%hmD=hm_mf(hmD,srs,meth,maxNrM,doPlot)
%
% hmD    - hm-structure
% meth   - method to use, options:
%          lognorm - multible lognormal fits, maximum number of modes is 4
%          peaks - locate peaks and evaluate fitting using peak tops only
% srs    - data source, if not provided fits all data:
%          eg. enter 'dmps' to fit only dmps data
% maxNrM - maximum number of modes
% doPlot - plot during the fitting, slower
%
%result will be in fit-fieled in hmD structure
%each day in separate cell
%     zs:  -  geometric meand diameter
%     ws:  -  width of the mode
%     Hs:  -  area of the mode (dN/dlogDp)
%     Ns:  -  number concentration of mode (dN)
%
% Example:
%
% hmD=hm_load(datenum(2007,4,4),'dmps');
% hmD=hm_mf(hmD,'dmps','lognorm',3,1);
% hm_plot(hmD,'dmps')
% hold on
% for i=1:3,plot(hmD.meta.dmps.tim{1},hmD.fit.dmps.zs{1}(:,i),'k.');end
% hold off

% if nargin<=3
%     startT=eval(['hmD.',srs,'(2,1)']);
%     endT=eval(['hmD.',srs,'(end,1)']);
% end

eval(['nrDays=length(hmD.',srs,');'])
if nargin<4
    maxNrM=3;
end
if nargin<5
    doPlot=0;
end
for d=1:nrDays
    tim=eval(['hmD.',srs,'{d}(2:end,1)']);
    dat=eval(['hmD.',srs,'{d}(2:end,3:end)']);
%     dp=eval(['hmD.',srs,'{d}(1,3:end)']); %diameter in original file (Mobility or Tammet or aerodynamic)
    dp=eval(['hmD.meta.',srs,'.dp{d}']); %diameter in Millikan diameter 
    tot=eval(['hmD.',srs,'{d}(2:end,2)']);

    % dlogDp=10.^([log10(dp(2:end))-log10(dp(1:end-1)),log10(dp(end))-log10(dp(end-1))]);
    % dat=dat.*repmat(dlogDp,size(dat,1),1);
    % Inan=isnan(dat);
    % dat(Inan)=0;
    %
    % dat=dat./repmat(sum(dat,2),1,size(dat,2));
    % dat=dat.*repmat(tot,1,size(dat,2));
    %
    % dat(Inan)=NaN;
    L=length(tim);
    h_out=cell(L,1);
    N=cell(L,1);
    param=cell(L,1);
    leg=zeros(L,1);
    switch meth

        case 'lognorm'
            for i=1:length(tim),
                try
                    Iok=~isnan(dat(i,:));
                    if i==1 & d==1
                        [param{i},fval,yhat,h_out{i},N{i}]=MF_lognorm(dat(i,Iok),dp(Iok),doPlot,maxNrM);
                    elseif i==1 & d>1
                        %take init from previous day
                        eval(['zs=hmD.fit.',srs,'.zs{d-1}(end,:);']);
                        eval(['ws=hmD.fit.',srs,'.ws{d-1}(end,:);']);
                        init=[ws,zs];
                        [param{i},fval,yhat,h_out{i},N{i}]=MF_lognorm(dat(i,Iok),dp(Iok),doPlot,maxNrM,init);
                    else
                        %                 [param{i},fval,yhat,h_out{i},N{i}]=MF_lognorm(dat(i,Iok),dp(Iok),0,4);
                        [param{i},fval,yhat,h_out{i},N{i}]=MF_lognorm(...
                            dat(i,Iok),dp(Iok),doPlot,maxNrM,param{i-1});
%                         [param{i},fval,yhat,h_out{i},N{i}]=MF_lognorm(dat(i,Iok),dp(Iok),doPlot,maxNrM);
                        %                 i
                        %                 pause
                    end
                catch
                    disp(lasterr)
                    param{i}=NaN(1,maxNrM*2);
                end
                %             drawnow
                leg(i)=length(param{i});
            end
            zs=NaN(length(param),max(leg)/2);
            ws=zs;
            Ns=zs;
            Hs=zs;

            for i=1:length(param),
                zs(i,1:leg(i)/2)=param{i}((leg(i)/2)+1:leg(i));
                ws(i,1:leg(i)/2)=param{i}(1:leg(i)/2);
                if ~isempty(N{i})
                    Ns(i,1:leg(i)/2)=N{i}(1:leg(i)/2);
                    Hs(i,1:leg(i)/2)=h_out{i}(1:leg(i)/2);
                end
            end
            [Isort,zs,Hs,ws,Ns]=sortModes(zs,Hs,ws,Ns);
        case 'peaks'
            doPlot=0;
            [Ps,Pts]=hm_pf(hmD,srs,doPlot);
            eval(['hmD.fit.',srs,'.tm{d}=Ps(:,1);']);
            zs=Ps(:,2);
            ws=Ps(:,4);
            Hs=Ps(:,3);
            Ns=Ps(:,5);
    end %switch

    eval(['hmD.fit.',srs,'.zs{d}=zs;']);
    eval(['hmD.fit.',srs,'.ws{d}=ws;']);
    eval(['hmD.fit.',srs,'.Hs{d}=Hs;']);
    eval(['hmD.fit.',srs,'.Ns{d}=Ns;']);
    eval(['hmD.fit.',srs,'.Algorithm=''',meth,''';']);
end %days

%% sortModes
function [Isort,zs,Hs,ws,Ns]=sortModes(zs,Hs,ws,Ns)
%
%sort modes accoring to mean diameter and area
%

zsn=(log10(zs)-log10(1e-9))./log10(999e-9); %scale between 1nm = 0 and 1000nm=1
Hsn=Hs./100000; %scale between 0 = 0 1=10000;
[n,m]=size(zs);
Isort=NaN(n,m);
zsOld=zs;
%arrange only by mode mean using matlab contest
if 0
   for i=2:n 
       
       z1=zs(i-1,:);z1(isnan(z1))=0;
       z2=zs(i,:);z2(isnan(z2))=0;
    
       moves = H_solver_ano(z2, z1, 10);
       movs=unique(moves,'rows');
       nrMov=size(movs,1);
       z2s= doMoves(z2,movs,nrMov);
       z2s(z2s==1)=NaN;
       
       zs(i,:)= doMoves(zs(i,:),movs,nrMov);
       Hs(i,:)= doMoves(Hs(i,:),movs,nrMov);
       ws(i,:)= doMoves(ws(i,:),movs,nrMov);
       Ns(i,:)= doMoves(Ns(i,:),movs,nrMov);

       %        [mn,I]=min(abs(repmat(z1,m,1)-repmat(z2',1,m)),[],1);
%              
%        Isort(i,:)=I;Isort(i,In)=NaN;
%        zs(i,:)=zs(i,I);zs(i,In)=NaN;
%        Hs(i,:)=Hs(i,I);Hs(i,In)=NaN;
%        ws(i,:)=ws(i,I);ws(i,In)=NaN;
%        Ns(i,:)=Ns(i,I);Ns(i,In)=NaN;
   end
end

%for loop through all modes starting from bigest
if 1
    NsTemp=Ns;
    NsTemp(isnan(Ns))=0;
    for i=2:n

        z1=zs(i-1,:);z1(isnan(z1))=0;
        z2=zs(i,:);z2(isnan(z2))=0;
        [s,Is]=sort(NsTemp(i,:),'descend');
        if any(z2==0)
            a=1;
        end
        z2in=zeros(1,m);
        Hsin=z2in;wsin=z2in;Nsin=z2in;
        for im=1:m
            [mn,Imn]=min(abs(z1-z2(Is(im))));
            z1(Imn)=1;
            z2in(1,Imn)=zs(i,Is(im));
            Hsin(1,Imn)=Hs(i,Is(im));
            wsin(1,Imn)=ws(i,Is(im));
            Nsin(1,Imn)=Ns(i,Is(im));
        end
        zs(i,:)=z2in;
        Hs(i,:)=Hsin;
        ws(i,:)=wsin;
        Ns(i,:)=Nsin;
    end
end

%arange by matching vectors
if 0
    for i=2:n
    a1=[zsn(i-1,:)',Hsn(i-1,:)'];
    a2=[zsn(i,:)',Hsn(i,:)'];

    
    [I,srtM]=H_matchVect(a1,a2);
    Isort(i,:)=I';
    zs(i,:)=zs(i,I);
    Hs(i,:)=Hs(i,I);
    ws(i,:)=ws(i,I);
    Ns(i,:)=Ns(i,I);
    % [orderInd,sortedMat,dis1,dis2,disInd]=H_matchVect(mat1,mat2)
end
end