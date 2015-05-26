function hmD=hm_fitFit(zs,Hs)
%
%find best fit for the parts
%
%
%


% AE_plot(AE_load(mfD.when,mfD.where),mfD.source)
% hold on
fixFit=mfD.fit;
col=['b','g','r','m','y'];

%split data and do fitting to the sections
if 1
    for m=mfD.nrM+1:size(mfD.fit,2)
        datY=log10(mfD.fit(:,m));

        datX=mfD.tm;
        [n,dummy]=size(datY);

        %divide in to hour parts
        interval=5;
        Ih=1:interval:n;
        if Ih(end)~=n,
            Ih=[Ih,n];
        end
        nh=length(Ih); %how many sections
        %     plot(datX,10.^datY,[col(m-3),'.'])
        plot(datX,10.^datY,'ko')
        for i=1:nh-1
            tempIh=Ih(i):Ih(i+1);
            x=datX(tempIh);
            y=datY(tempIh);
            %         x=datX(Ih(i):Ih(i+1));
            %         y=datY(Ih(i):Ih(i+1));
            %
            nInt=length(tempIh);
            [p,s]=polyfit(x,y,1);
            yhat{i}=10.^(x*p(1)+p(2));


            resid=10.^y-yhat{i};
            Iok=find(abs(resid)<(mean(abs(resid))+1*std(abs(resid))));
            if length(Iok)~=nInt;
                [p,s]=polyfit(x(Iok),y(Iok),1);
                yhat{i}=10.^(x(Iok)*p(1)+p(2));
            end
            if s.normr<0.2
                plot(x(Iok),10.^y((Iok)),[col(m-mfD.nrM),'.'])
                plot(x(Iok),yhat{i},col(m-mfD.nrM))
                fixFit(tempIh(Iok),m)=yhat{i};
                %         else
                %             plot(x(Iok),10.^y((Iok)),'.w')
                %             plot(x(Iok),yhat{i},'w')

            end
        end
    end
    hold off
    mfD.fixFit=fixFit;
end