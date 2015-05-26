function [out]=AE_conc(hmD,srs,dP_lim)
%
% calculate integral-concentration of given range
%
% [c]=hm_conc(hmD,srs,dP_lim)
%
% c - [time,data] - total number concentration between limits
% hmD - data structure loaded by hm_load
% srs - instrument name, string
% dP_lim - 2-number vector of aerosol size limits in meters
%
%EXAMPLE:
%calculate number concentration between 3 and 20 nanometers
% [c]=hm_conc(hmD,'dmps',[3,20]*1e-9)

% 14.03.2007
% Heikki Junninen

if ~strcmp(hmD.Type,'hm_dat')
    error('hm_conc: wrong data structure, use hm_load for loading!')
end

data=eval(['hmD.',srs]);
timeC=eval(['hmD.meta.',srs,'.tim']);
out=[];
if iscell(data)
    for d=1:length(data)
        if ~isempty(data{d})
            %         tm=data{d}(2:end,1);
            tm=timeC{d};
%             dp=log10(data{d}(1,3:end))';
            dp=eval(['log10(hmD.meta.',srs,'.dp{d});']);
            % dat=medfilt1(data(2:end,3:end),3);
            dat=data{d}(2:end,3:end);

            dmin=max(log10(dP_lim(1)),dp(1));
            dmax=min(log10(dP_lim(2)),dp(end));
            % disp(['Minimum diameter is: ' num2str(1e9*10 .^dmin,3) ' nm'])
            % disp(['Maximum diameter is: ' num2str(1e9*10 .^dmax,3) ' nm'])

            dpi=[dmin:0.001:dmax]';
            dat(isnan(dat))=0;
            try
                conci=sum(interp1(dp,dat',dpi)*0.001);
                %plot(time,conci,time,v(2:m,2))
                %grid
                %xlabel('time')
                %ylabel('concentration [#/cm^3]')
                %legend('fraction','total',0)
                if length(conci)>1
                    ou=[tm, conci'];
                    out=[out;ou];
                end
            end
        end
    end

else
    tm=data(2:end,1);
    dp=log10(data(1,3:end))';
    % dat=medfilt1(data(2:end,3:end),3);
    dat=data(2:end,3:end);

    dmin=max(log10(dP_lim(1)),dp(1));
    dmax=min(log10(dP_lim(2)),dp(end));
    % disp(['Minimum diameter is: ' num2str(1e9*10 .^dmin,3) ' nm'])
    % disp(['Maximum diameter is: ' num2str(1e9*10 .^dmax,3) ' nm'])

    dpi=[dmin:0.001:dmax]';
    conci=sum(interp1(dp,dat',dpi)*0.001);

    %plot(time,conci,time,v(2:m,2))
    %grid
    %xlabel('time')
    %ylabel('concentration [#/cm^3]')
    %legend('fraction','total',0)

    out=[tm, conci'];
end