function [AEstr,pmfR]=hm_modefitting(AEstr,srs)
% Mode fitting using PMF
%based on MF_pmf_v6
% pmfR=MF_pmf_v6(AEstr,srs)
%
%

%symetry check added to mode combine part
sig1=1.18; %1.18
deltaZ=0.08; %0.08

if nargin==1
%     srs='all_autom';
    instruments={...
        'dmps',...
        'aisn',...
        'aisp',...
        'cpc',...
        'naisn',...
        'naisp',...
        'bsman',...
        'bsmap',...
        'volh',...
        'volc'...
        };

    h=0;
    prcInstr=[];
    for i=1:length(instruments),
        if isfield(AEstr,instruments{i})
            h=h+1;
            prcInstr{h}=instruments{i};
        end
    end
    if isempty(prcInstr)
        disp('hm_modefitting: no data to plot')
        return
    end
else
    prcInstr{1}=srs;
end


for in=1:length(prcInstr)
    sr=prcInstr{in};
    switch sr
        case 'aisn'
            dat=AEstr.aisn(2:end,3:end);
            %         dat=H_2dmedfilt(dat,[1,5]);
            dP=AEstr.aisn(1,3:end);
            tm=AEstr.aisn(2:end,1);
            tm=AEstr.meta.aisn.tim{1};
        case 'aisp'
            dat=AEstr.aisp(2:end,3:end);
            %         dat=H_2dmedfilt(dat,[1,5]);
            dP=AEstr.aisp(1,3:end);
            tm=AEstr.aisp(2:end,1);

        case 'naisn'
            dat=AEstr.naisn(2:end,3:end);
            dP=AEstr.naisn(1,3:end);
            tm=AEstr.naisn(2:end,1);

        case 'naisp'
            dat=AEstr.naisp(2:end,3:end);
            dP=AEstr.naisp(1,3:end);
            tm=AEstr.naisp(2:end,1);
        case 'nais4n'
            dat=AEstr.nais4n(2:end,3:end);
            dP=AEstr.nais4n(1,3:end);
            tm=AEstr.nais4n(2:end,1);

        case 'nais4p'
            dat=AEstr.nais4p(2:end,3:end);
            dP=AEstr.nais4p(1,3:end);
            tm=AEstr.nais4p(2:end,1);
        case 'nais4nA'
            dat=AEstr.nais4nA(2:end,3:end);
            dP=AEstr.nais4nA(1,3:end);
            tm=AEstr.nais4nA(2:end,1);

        case 'nais4pA'
            dat=AEstr.nais4pA(2:end,3:end);
            dP=AEstr.nais4pA(1,3:end);
            tm=AEstr.nais4pA(2:end,1);

        case 'ais'
            dat=AEstr.ais(2:end,3:end);
            %         dat=H_2dmedfilt(dat,[1,5]);
            dP=AEstr.ais(1,3:end);
            tm=AEstr.ais(2:end,1);

        case 'dmps'
            dat=AEstr.dmps{1}(2:end,3:end);
            %         dat=H_2dmedfilt(dat,[1,3]);
            dP=AEstr.dmps{1}(1,3:end);
            tm=AEstr.dmps{1}(2:end,1);
            sig1=1.26;
            deltaZ=0.088;
        case 'comb'
            if isempty(AEstr.comb)
                AEstr=AE_interp(AEstr);
            end
            sig1=1.18;
            deltaZ=0.088;
            if isempty(AEstr.comb)
                pmfR=[];
                return
            end
            dat=AEstr.comb(2:end,3:end);
            dP=AEstr.comb(1,3:end);
            tm=AEstr.comb(2:end,1);

        case 'volC'
            dat=AEstr.volC(2:end,3:end);
            dP=AEstr.volC(1,3:end);
            tm=AEstr.volC(2:end,1);

        case 'volH'
            dat=AEstr.volH(2:end,3:end);
            dP=AEstr.volH(1,3:end);
            tm=AEstr.volH(2:end,1);

        case 'test'
            dat=AEstr.dat;
            dP=AEstr.dP;
            tm=AEstr.tim;
        otherwise
            error('unknown source')
    end

    %find rows with all values negative
    %and replace them by NaNs
    dat(all(dat<0,2),:)=NaN;
    dat(all(dat>1e5,2),:)=NaN;

    doPlot=0;
    doCalc=1;
    if doCalc
        %     x=exp(-20.5:.05:-13);
        %     deltaZ=0.1;
        %     z=exp(-20.2:deltaZ:-13.2);

        x=exp(-21.7:.05:-13);
        %     deltaZ=(log(dP(end))-log(dP(1)))/(length(dP)*1.3);
        %     deltaZ=(log(dP(end))-log(dP(1)))/(length(dP)*3);
        %      deltaZ=0.08;
        %      deltaZ=(log(dP(end))-log(dP(1)))/(length(dP)*1.5);
        %     z=exp(-21.9:deltaZ:-13.2);
        %     newDp=exp(log(dP(1)):(log(dP(end))-log(dP(1)))/69:log(dP(end)));
        newDp=exp(log(dP(1)):deltaZ:log(dP(end))*1.0001);

        %     newDp=exp(log(dP(1)):0.0838:log(dP(end)));
        % x=exp(log(dP(1))*0.9:.05:log(dP(end))*1.1);
        %     deltaZ=0.1;
        %     z=exp(log(dP(1))*0.95:deltaZ:log(dP(end))*1.05);
        z=exp(log(dP(1))*0.99:deltaZ:log(dP(end))*1.01);


        p=length(z);

        %     sig1=1.4;
        %     sig2=1.4;
        % Best berformance when ratio deltaZ=0.2 and sig1=1.11
        %     sig1=1.12; %def 1.25 %assigned above
        %     sig1=1.15;
        sig2=sig1;

        %     f=[ones(1,p)*1.5,z];
        f1=[ones(1,p)*sig1,z];
        f2=[ones(1,p)*sig2,z];


        itStr=struct('dat',ones(1,length(x)),'x',x,'xx',[],'mpl',ones(length(x),1));
        peaks=struct('nr_peaks',p);
        itStr.peaks=peaks;


        %[fval,y1,h_out]=H_lognorm2pAE2(itStr,f);y=squeeze(y1)';
        [fval,y1,h_out]=H_lognorm2pAE2(itStr,f1);ys1=squeeze(y1)';
        [fval,y1,h_out]=H_lognorm2pAE2(itStr,f2);ys2=squeeze(y1)';

        ys1=som_normalize(ys1,'range');
        ys2=som_normalize(ys2,'range');

        %     %find closest bins in real data
        %      [mn Imn]=min(abs(repmat(x',1,length(dP))-repmat(dP,length(x),1)));
        %
        %     %find closest bins in 3x higher resolution than real data
        %  [mn Imn]=min(abs(repmat(x',1,length(x(1:2:end)))-repmat(x(1:2:end),length(x),1)));
        [mn Imn]=min(abs(repmat(x',1,length(newDp))-repmat(newDp,length(x),1)));
        inAAhgh=ys2(Imn,:)*1;
        inAAlow=ys1(Imn,:)*1;

        Isml=find(sum(inAAhgh)<1);

        inAAhgh(:,Isml)=[];
        inAAlow(:,Isml)=[];
        z(Isml)=[];
        p=size(inAAhgh,2);
        inAA=inAAhgh;
        
        len=size(dat,2);
        minD=repmat(min(dat')',1,len);
        %     datn=(dat-minD)./repmat(max(dat')',1,len);

        %outlayers replaced by 0 and are treated like BDL
        %     dat(find((mean(dat(:))+4*std(dat(:)))<dat))=0;

        dati=AE_interp_test(dat,dP,newDp);
        dati=max(dati,1e-3);
        %data not normalized
        
        pmfD=Hpmf_dat(dati,p,'inAAlow',inAAlow(:,1:p),'inAAhgh',inAAhgh(:,1:p),...
            'inAA',inAA,'normA',2,'maxit',800,'DL',10,'smth',1);
        pmfR=Hpmf_opt(pmfD);
        %          [newZs,sym,N,Nrel]=MF_pmf_combPeaks(pmfR,z,dP);
        [newZs,sym,N,Nrel,Icomb]=MF_pmf_combPeaks(pmfR,newDp,dP,z,deltaZ);
        %     pmfR.when=AEstr.FNam;
        %     pmfR.where=AEstr.Path;
        pmfR.source=sr;
        pmfR.nrM=size(newZs,2);
        pmfR.dP=newDp;
        %     pmfR.srs=srs;

        fit.z=[newZs];
        fit.sym=[sym];
        fit.N=[N];
        fit.Nrel=[Nrel];
        fit.tm=repmat(tm,1,size(newZs,2));
        fit.Icomb=Icomb;
        pmfR.fit=fit;
        eval(['AEstr.fit.',sr,'=fit;']);
    end

    %plot
    if doPlot
        AE_plot(AEstr,srs),
        hold on,
        Isel=N>10;
        %              Isel=N>0;
        t=repmat(tm,1,size(newZs,2));
        plot(t(Isel),(newZs(Isel)),'k.'),
        %     plot(t(Isel),(newZs(Isel)./sym(Isel)),'w*'),
        hold off
    end
end %for instruments

%SUBFUNCTION============================================================
function [newZs,syms,Ns,Nrels,Icomb]=MF_pmf_combPeaks(pmfR,newDp,dP,z,deltaZ)
%
% Combine peaks from pmf fitting
%
%

newZs=zeros(size(pmfR.data,1),10);
syms=zeros(size(pmfR.data,1),10);

IallDist=1:size(pmfR.G,2);
doPlot=0;

% sData=sgolayfilt(pmfR.data,2,21,[],1);
% sData=sgolayfilt(sData,2,5,[],2);
%
% pmfR.data=sData;

for tm=1:1:size(pmfR.data,1);
    %       for tm=144
    if ~all(isnan(pmfR.data(tm,:)))

        nrSmlDis=size(pmfR.F,2);
        %derivate
        dat=pmfR.G(tm,:)*pmfR.F';
        %         dat=sData(tm,:);
        df=diff(dat);

        %local maximum
        %using derivate
        %         Ipeak=find(df(1:end-1)>0 & df(2:end)<0);

        %using maalaisjärki
        Ipeak=find(dat(1:end-2)<dat(2:end-1) & dat(2:end-1)>dat(3:end))+1;
        if dat(1)>dat(2),Ipeak=[1,Ipeak];end
        if dat(end)>dat(end-1),Ipeak=[Ipeak,size(dat,2)];end

        if 1
            %look for a local minimas in G
            tempG=pmfR.G(tm,:);

            %         tempGs=sgolayfilt(tempG,2,7,[],2);
            %         tempG=[0.4*tempGs+0.6*tempG];
            %             difFact=0.8;
            difFact=1;
            mins=(tempG(2:end-1)<tempG(3:end)*difFact & tempG(1:end-2)*difFact>tempG(2:end-1));
            %            if tm==222
            %                a=1;
            %            end
            if mins(end)==length(tempG)-2
                mins=[1,mins,0]; %mins(2)=0;mins(end-1)=0;
            else
                mins=[1,mins,1]; %mins(2)=0;mins(end-1)=0;
            end
            %             Imins=find(mins);
            Imins=find(mins | pmfR.G(tm,:)==min(pmfR.G(tm,:)));

            Imins(Imins(1:end-1)-Imins(2:end)==-1)=[]; %remove adjasent mins
            nrP=length(Imins)-1;
            newZ=[];
            newZM=[];
            N=[];
            Nrel=[];
            %             Icomb=[];
            i=0;
            while i<nrP
                %             for i=1:nrP,
                %             plot(pmfR.G(tm,[Imins(i):Imins(i+1)])*pmfR.F(:,[Imins(i):Imins(i+1)])');
                %             hold on
                i=i+1;
                if Imins(i+1)==Imins(end)
                    Icomb{tm}{i}=Imins(i):Imins(i+1);
                else
                    Icomb{tm}{i}=Imins(i):Imins(i+1)-1;
                end

                %calc symmetry parameter
                Idist=Icomb{tm}{i}(find(pmfR.G(tm,Icomb{tm}{i})>0));
                %             newZM(pks)=exp(sum(log(z(Idist)).*pmfR.G(tm,Idist))./sum(pmfR.G(tm,Idist)));
                w=pmfR.G(tm,Idist);
                newZM(i)=sum(z(Idist).*w)/sum(w);
                dat=pmfR.G(tm,Icomb{tm}{i})*pmfR.F(:,Icomb{tm}{i})';

                [sm,nZ]=calcSym(dat,newZM(i),dP,deltaZ);

                %                  [mx,Imx]=max(dat);
                %                  sm=sum(dat(1:Imx))/sum(dat(Imx:end));

                %                 newZ(i)=nZ;
                sym(i)=sm;
                %                 plot(dat)
                %shoulder detection if mode very skewed
                if tm==120
                    a=1;
                end
                if (sm>1.20 | sm<0.80) & length(Icomb{tm}{i})>2
                    doSoulderDetection=1;
                else
                    doSoulderDetection=0;
                    Ishou{i}=[];
                end
                if doSoulderDetection

                    %NEW; iterate till both peaks are symmetric

                    %fronting, start iteration from front
                    if sm<1
                        count=0;
                        %                         [mx,Imx]=max(dat);
                        %                         mxdP=newdP(Imx);
                        for sh=2:length(Icomb{tm}{i}),
                            count=count+1;
                            dat=pmfR.G(tm,Icomb{tm}{i}(1:sh))*pmfR.F(:,Icomb{tm}{i}(1:sh))';
                            %                             [sh_sm(count),sh_nZ(count)]=calcSym(dat,mxdP,dP,deltaZ);
                            [mx,Imx]=max(dat);
                            sh_sm2(count)=sum(dat(1:Imx))/sum(dat(Imx:end));
                            %                             subplot(5,5,sh)
                            %                             plot(dat),
                        end
                        Iok=find((sh_sm2(1:end-1)>1 & sh_sm2(2:end)<1) | (sh_sm2(1:end-1)<1 & sh_sm2(2:end)>1));

                        %solution closest to one
                        if ~isempty(Iok)
                            Imn=Iok(end);
                            dummy=Icomb{tm}{i};
                            Icomb{tm}{i}=dummy(1:Imn);
                            Ishou{i}=dummy(Imn+1:end);
                        else
                            Ishou{i}=[];
                        end
                        %tailing, start iteration from back
                    elseif sm>1
                        count=0;
                        for sh=length(Icomb{tm}{i}):-1:2,
                            count=count+1;
                            dat=pmfR.G(tm,Icomb{tm}{i}(sh:end))*pmfR.F(:,Icomb{tm}{i}(sh:end))';
                            %                             [sh_sm(count),sh_nZ(count)]=calcSym(dat,newZM(i),dP,deltaZ);
                            %                             subplot(5,5,sh)
                            %                             plot(dat),
                            [mx,Imx]=max(dat);
                            sh_sm2(count)=sum(dat(1:Imx))/sum(dat(Imx:end));
                        end
                        Iok=find((sh_sm2(1:end-1)>1 & sh_sm2(2:end)<1) | (sh_sm2(1:end-1)<1 & sh_sm2(2:end)>1));

                        %solution closest to one
                        if ~isempty(Iok)
                            Imn=Iok(1);
                            dummy=Icomb{tm}{i};
                            Ishou{i}=dummy(1:Imn-1);
                            Icomb{tm}{i}=dummy(Imn:end);
                        else
                            Ishou{i}=[];
                        end
                    end
                    clear sh_sm Imn Imx Iok sh_sm2
                end %shoulder
            end %while
        end %1

        %         %if last one is not selected
        %         if Icomb{tm}{end}(end)~=nrSmlDis
        %             nrP=nrP+1;
        %             Icomb{tm}{nrP}=IfrstMemb:nrSmlDis;
        %         end

        %calculate new Z
        %include for new Z calculations only the small distributions that
        %have G bigger than 0.2
        %consist of 90% of total area
        %
        nrP1=nrP;
        for pks=1:nrP1
            if ~isempty(Ishou{pks})
                Icomb{tm}=[Icomb{tm},Ishou(pks)];
                nrP=nrP+1;
            end
        end

        for pks=1:nrP,
            %             Idist=Icomb{pks}(find(Icomb{pks}>0.3));
            Idist=Icomb{tm}{pks}(find(pmfR.G(tm,Icomb{tm}{pks})>0));
            %             newZM(pks)=exp(sum(log(z(Idist)).*pmfR.G(tm,Idist))./sum(pmfR.G(tm,Idist)));
            w=pmfR.G(tm,Idist);
            newZM(pks)=sum(z(Idist).*w)/sum(w);

            dat=pmfR.G(tm,Icomb{tm}{pks})*pmfR.F(:,Icomb{tm}{pks})';
            datx=1:length(dat);
            [s Is]=sort(dat,'descend');

            datF=dat(sort(Is(1:6)));
            xF=datx(sort(Is(1:6)));
            x=xF(1):0.01:xF(end);

            %             p=polyfit(xF,datF,3);
            %             peak=p(1)*x.^3+p(2)*x.^2+p(3)*x+p(4);

            p=polyfit(xF,datF,2);
            %             peak=p(1)*x.^2+p(2)*x+p(3);

            %                         nrM=length(par2)/2;
            %             limLow=[ones(1,nrM)*1.1,par2(nrM+1:end)-0.15];
            %             limHgh=[ones(1,nrM)*1.5,par2(nrM+1:end)+0.15];
            %
            %             [par2]=fmincon(@(par2) H_lognorm2p(par2,xy),par2,[],[],[],[],limLow,limHgh,@H_lognormCon,opt);
            %             [fval2,ytot,y1,h_out]=H_lognorm2p(par2,xy);


            %             [mx Imx]=max(peak);
            %             newZ(pks)=exp((x(Imx)-1)*deltaZ+log(dP(1)));

            %works only for 2order polynome
            newZ(pks)=exp((mean(roots(p))-1)*deltaZ+log(dP(1)));

            %if symmetry very skewed, devide the peak


            %                         plot(dat,'.'),hold on
            %                         plot(xF,datF,'.k'),
            %                         plot(x,peak,'r'),hold off
            %                         drawnow
            %                         pause(0.15)
            %
            %calculate the area, partcls/cc/log10(dP)
            N(pks)=trapz(log10(newDp),pmfR.G(tm,Icomb{tm}{pks})*pmfR.F(:,Icomb{tm}{pks})');
            %relative concentration on this mode
            Nrel(pks)=N(pks)/trapz(log10(newDp),pmfR.G(tm,:)*pmfR.F(:,:)');
        end

        sym=newZ./newZM;

        if doPlot
            %             plot(z,pmfR.G(tm,:)*pmfR.F(:,:)','-k')
            h1=plot(newDp,pmfR.G(tm,:)*pmfR.F(:,:)','-k');
            set(h1,'linewidth',2)
            hold on
            h2=plot(newDp,pmfR.data(tm,:),'k.');
            set(h2,'markersize',10)
            for pks=1:nrSmlDis,
                h3=plot(newDp,pmfR.G(tm,pks)*pmfR.F(:,pks)','color',[0.5 0.5 0.5]);
                set(h3,'linewidth',2)
            end
            for pks=1:nrP,
                h4=plot(newDp,pmfR.G(tm,Icomb{tm}{pks})*pmfR.F(:,Icomb{tm}{pks})','r');
                set(h4,'linewidth',2)
            end
            % plot(mfD.dP,pmfR.G(tm,Icomb{2})*pmfR.F(:,Icomb{2})','b')
            plot(newZ,max(pmfR.data(tm,:)),'*r')
            plot(newZM,max(pmfR.data(tm,:)),'*b')
            axis([min(newDp) max(newDp) 0 max(pmfR.data(:))])
            set(gca,'xscale','log')
            hold off
            %             title(['data row ',num2str(tm)])
            drawnow
            pause(0.1)
        end

        Iok=newZ>=dP(1) & newZ<=dP(end);

        newZs(tm,1:nrP)=newZ.*Iok;
        syms(tm,1:nrP)=sym.*Iok;
        Ns(tm,1:nrP)=N.*Iok;
        Nrels(tm,1:nrP)=Nrel.*Iok;

        clear Ipeak newZ newZM sym N Nrel Ishou sumAbsDif IsmlPeak overLapParent
    end
end
newZs(:,all(newZs==0,1))=[];
syms(:,all(syms==0,1))=[];
Ns(:,all(Ns==0,1))=[];
Nrels(:,all(Nrels==0,1))=[];

%==========================================================================
%SUBFUNCTIONS

function [sym,newZ]=calcSym(dat,newZM,dP,deltaZ)
%
%Calculate symmetry parameter
%

datx=1:length(dat);
[s Is]=sort(dat,'descend');

datF=dat(sort(Is(1:6)));
xF=datx(sort(Is(1:6)));
x=xF(1):0.01:xF(end);

p=polyfit(xF,datF,2);

%works only for 2order polynome
newZ=exp((mean(roots(p))-1)*deltaZ+log(dP(1)));
sym=newZ./newZM;
