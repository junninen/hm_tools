function hmD=hm_normalize(hmD,varargin)
%
% Normalize data so that the total concentration is 1
%
%

srs='dmps';
%varargin
i=1;
while i<=length(varargin),
    argok = 1;
    if ischar(varargin{i}),
        switch varargin{i},
            % argument IDs
            case 'srs',                 i=i+1; srs = varargin{i};
            otherwise, argok=0;
        end
    else
        argok = 0;
    end
    if ~argok,
        disp(['Ignoring invalid argument #' num2str(i+1)]);
    end
    i = i+1;
end

eval(['nrD=length(hmD.',srs,');'])
for i=1:nrD
    eval(['d=hmD.',srs,'{i};'])
    c=hm_conc(hmD,srs,[1e-9,1000e-9]);
    totC=nansum(c(:,2));
    d(2:end,3:end)=d(2:end,3:end)/totC;
    eval(['hmD.',srs,'{i}=d;'])
end

hmD=hm_history(hmD,mfilename);