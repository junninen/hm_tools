function hmD=hm_smooting(hmD,meth,wins,srs)
%
% smoothing aerosol data
% hmD=hm_smoothing(hmD,meth,wins,srs)
%
%   hmD  (struct) data structure loaded with hm_load
%   srs  (string) name of source, default is all instruments
%                 different instrument names see hm_load help.
%   meth (string) name of smoothing methods,
%                 '2dmedfilt' or '2dmeanfilt'
%   wins (double) window for smoothing, size depends on method used
%                 for 2dmedfilt - [w1, w2], w2=time window, steps, w1=particle
%                 size window, bins, see help H_2dmedfilt for help
%   new data will replace the old one.
%   Rawdata is back-uped to meta.raw
%
%
% if wins are not given evaluate them using durbin_watson index
%
% see also hm_load, H_2dmedfilt

%Heikki Junninen
%18.03.2007
%update 15.10.2007
%update 20.11.2007
%   - added multible smoothing pile up
%   - copy raw data only during first time

if ~strcmp(hmD.Type,'hm_dat')
    error('hm_smoothing: wrong data structure, use hm_load for loading!')
end

if nargin==3
    srs=[];
end

if nargin==2
    wins=[];
    srs=[];
end

if nargin==1
    meth='2dmedfilt';
    wins=[];
    srs=[];
end

instruments=hmD.meta.insts;
nrIns=length(instruments);

%smoothInstr=zeros(nrIns,1);

if isempty(srs)
    h=0;
    for i=1:nrIns,
        if isfield(hmD,instruments{i})
            h=h+1;
            smoothInstr{h}=instruments{i};
        end
    end
else
    Iinstr=strcmp(instruments,srs);
    smoothInstr=instruments(Iinstr);

    if isempty(smoothInstr)
        error('hm_smoothing: could not find instrument. Check typing!');
    end
end

for i=1:length(smoothInstr)
    eval(['doneAlready=hmD.meta.',smoothInstr{i},'.smoothing;']);

    data=eval(['hmD.',smoothInstr{i}]);
    datas=data;

    if ~iscell(data)
        %if wins are not given evaluate them using durbin_watson index
        if isempty(wins)
            w=3:2:9;
            for h1=1:length(w)
                for h2=1:length(w)
                    dws(h1,h2)=hm_durbin_watson_index(data(2:end,3:end),H_2dmedfilt(data(2:end,3:end),[w(h1),w(h2)]));
                    w1(h1,h2)=w(h1);
                    w2(h1,h2)=w(h2);
                end
            end
            [mn,Imn1]=min(abs(dws(:)-2));
            wins=[w1(Imn1),w2(Imn1)];
        end
        datas(2:end,3:end)=H_2dmedfilt(data(2:end,3:end),wins);

    else
        for di=1:length(data)
            if isempty(wins)
                w=3:2:9;
                for h1=1:length(w)
                    for h2=1:length(w)
                        dws(h1,h2)=hm_durbin_watson_index(data{di}(2:end,3:end),H_2dmedfilt(data{di}(2:end,3:end),[w(h1),w(h2)]));
                        w1(h1,h2)=w(h1);
                        w2(h1,h2)=w(h2);
                    end
                end
                [mn,Imn1]=min(abs(dws(:)-2));
                wins=[w1(Imn1),w2(Imn1)];
            end

            switch meth
                case '2dmedfilt'
                    datas{di}(2:end,3:end)=H_2dmedfilt(data{di}(2:end,3:end),wins);
                case '2dmeanfilt'
                    datas{di}(2:end,3:end)=H_2dmeanfilt(data{di}(2:end,3:end),wins);
                otherwise
                    disp('hm_smoothing: Unknown method. Nothing changed!')
                    return
            end
        end
    end

    if doneAlready
        %if done already dont copy raw data to meta
        %and expand the smoothing and win matrix
        eval(['hmD.meta.',smoothInstr{i},'.smoothing=[hmD.meta.',smoothInstr{i},'.smoothing;1];']);
        eval(['hmD.meta.',smoothInstr{i},'.win=[hmD.meta.',smoothInstr{i},'.win;wins];']);
        eval(['hmD.meta.',smoothInstr{i},'.smoothMeth=[hmD.meta.',smoothInstr{i},'.smoothMeth;{meth}];']);
            eval(['hmD.',smoothInstr{i},'=datas;']);
else
        eval(['hmD.meta.',smoothInstr{i},'.raw=data;']);
        eval(['hmD.meta.',smoothInstr{i},'.smoothing=1;']);
        eval(['hmD.meta.',smoothInstr{i},'.win=wins;']);
        eval(['hmD.meta.',smoothInstr{i},'.smoothMeth={''',meth,'''};']);
        eval(['hmD.',smoothInstr{i},'=datas;']);
    end %if doneAlready

end

hmD=hm_history(hmD,mfilename);