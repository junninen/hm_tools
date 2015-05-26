function hmD=hm_mf_number_conc(hmD)
%
%Calculate number concentration for fitted modes
% hmD=hm_mf_number_conc(hmD)
%
% Replaces old Ns (dN/dlogDp) with new values (dN/mode, cm-3)!
%

%Heikki Junninen
% 22.May.2008


for h=1:length(hmD.dmps(2:end,1))
    dat=hmD.dmps(h+1,3:end);
    dp=hmD.dmps(1,3:end);
    param=[hmD.fits.dmps.ws(h,:),hmD.fits.dmps.zs(h,:)];
    nrPeaks=size(hmD.fits.dmps.ws,2);

    [fvalF,y1_F,yF,h_chF]=H_lognorm2pAE4(dat,dp,nrPeaks,param,1);

    h_out=h_chF;

    calcNumConc=1;
    N=h_out;
    if calcNumConc
        dmin=log10(dp(1));
        dmax=log10(dp(end));
        logDp=log10(dp);

        dpi=[dmin:0.001:dmax]';
        for i=1:nrPeaks
            N(i)=sum(interp1(logDp',squeeze(y1_F(:,i,:)),dpi)*0.001);
        end
    end
    hmD.fits.dmps.Nns(h,:)=N;
end
