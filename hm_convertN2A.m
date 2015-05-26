function hmD=hm_convertN2A(hmD)
%convert number concentration to volume concentration
%
%
%

for i=1:length(hmD.dmps)
    dP=hmD.dmps{i}(1,3:end);
    dN=hmD.dmps{i}(2:end,3:end);
    dA=dN.*pi.*repmat(dP.^2,size(dN,1),1);
    hmD.dmps{i}(2:end,3:end)=dA;

    %calculate the total
    dmin=log10(dP(1));
    dmax=log10(dP(end));
    
    % disp(['Minimum diameter is: ' num2str(1e9*10 .^dmin,3) ' nm'])
    % disp(['Maximum diameter is: ' num2str(1e9*10 .^dmax,3) ' nm'])

    dpi=[dmin:0.001:dmax]';
    conci=sum(interp1(log10(dP),dA',dpi)*0.001);
    hmD.dmps{i}(2:end,2)=conci;
    
    %     dS=dN.*pi.*repmat(dP.^2,size(dN,1),1);
    %     hmDs.dmps{i}(2:end,3:end)=dS;

end


