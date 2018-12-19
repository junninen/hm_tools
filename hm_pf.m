function [Ps,Pts]=hm_pf(hmD,srs,doPlot,doArea)
%
%mode fitting with derivatives
%

eval(['nrDays=length(hmD.',srs,');'])
if nargin<3
    doPlot=0;
    doArea=1;
end
if nargin<4
    doArea=1;
end
for d=1:nrDays
    %     tim=eval(['hmD.',srs,'{d}(2:end,1)']);
    tim=eval(['hmD.meta.',srs,'.tim{d}']);
    dat=eval(['hmD.',srs,'{d}(2:end,3:end)']);
    %     dp=eval(['hmD.',srs,'{d}(1,3:end)']); %diameter in original file (Mobility or Tammet or aerodynamic)
    dp=eval(['hmD.meta.',srs,'.dp{d}']); %diameter in Millikan diameter
    tot=eval(['hmD.',srs,'{d}(2:end,2)']);
    Ps=[];
    for i=1:length(tim),
        try
            Iok=~isnan(dat(i,:));
            dati=dat(i,Iok);
            dpi=log10(dp(Iok));
            
            %             varargout=tof_locate_peaks(mz,spec,param)
            SlopeThreshold=.01;
            AmpThreshold=10;
            smoothwidth=3;
            peakgroup=4;
            
            [P,d]=findpeaks(dpi,dati,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup);
            P(:,1)=tim(i);
            
            N=NaN(size(P,1),1);

            % fit width on linear scale
            if doArea
                global z_in
                z_in=10.^(P(:,2))';
                w_in=ones(size(z_in))*1.1;
                
                xy=[dp(Iok)',dati'];
                maxIt=1000;
                opt = optimset('GradObj','off','Display','off','maxiter',maxIt);
                [par1]=fminsearch(@(par1) H_lognorm1p_globalZ(par1,xy),w_in,opt);
                [fval,ytot,y1,h_out]=H_lognorm1p_globalZ(par1,xy);
                
                if doPlot
                    plot(dp(Iok),dati,'.k')
                    hold all
                    plot(dp,squeeze(y1))
                    hold off
                    set(gca,'xscale','log')
                    drawnow
                    %             pause(0.3)
                end
                % mode area
                for p=1:size(P,1)
                    %         param=[h,w,z];
                    w=par1(p);
                    z=z_in(p);
                    param=[w,z];
                    N(p)=mode_area(dati,dp(Iok),param);
                    
                end
                P(:,4)=par1';
            P(:,2)=z_in;
            end
            Ps=[Ps;[P,N]];
            
            %                         figure(1221)
            %             clf
            %             plot(dpi,dati)
            %             line(P(:,2),P(:,3),'marker','.','color','r','LineStyle','-')
            % %             pause(0.05)
        end
    end
    
    Pts=[];
    for t=1:length(dp)
        if dp(t)<25e-9
            dati=dat(:,t);
            
            SlopeThreshold=.01;
            AmpThreshold=10;
            smoothwidth=3;
            peakgroup=5;
            
            [Pt,d]=findpeaks(tim,dati,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup);
            Pt(:,1)=dp(t);
            
            %         plot(tim,dati,'.-'),hold on,plot(Pt(:,2),Pt(:,3),'or'),hold off
            
            Pts=[Pts;Pt];
        end
    end
end

%% calculate mode areas

function N=mode_area(dat,dp,param)

%fit with final parameters, and original data
scaleNegErr=1; %has effect only on error calculation, not used here
nrPeaks=1;
[fvalF,y1_F,yF,h_chF]=H_lognorm2pAE4(dat,dp,nrPeaks,param,scaleNegErr);

calcNumConc=1;
N=[];
if calcNumConc
    dmin=log10(dp(1));
    dmax=log10(dp(end));
    logDp=log10(dp);
    
    dpi=[dmin:0.001:dmax]';
    for i=1:nrPeaks
        N(i)=sum(interp1(logDp',squeeze(y1_F(:,i,:)),dpi)*0.001);
    end
end
