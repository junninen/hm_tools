function hmD=hm_reset(hmD)
%
%Return original raw data
%
%  hmD=hm_reset(hmD,srs)
%

%Heikki Junninen 29.03.2007
%

disp('hm_reset: NOT READY')
if nargin==1
%check which instruments are available
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
    procInstr=[];
    for i=1:length(instruments),
        if isfield(AEstr,instruments{i})
            h=h+1;
            procInstr{h}=instruments{i};
        end
    end
else
    
    
    
end
        
        