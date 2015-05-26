function xtick=H_xdatetick_int(int,labInt,mode,ax)
%set x tick to sertain interval (hours)
%H_xdatetick_int(int,[labInt],[ax])
%int - interval of ticks and grids
%labInt - interval of xticklabels
%
% Alternative usage:
% if int given as charagter sets tick at beginning of day, month or year
% for days
% H_xdatetick_int('d')
%
% every second day
% H_xdatetick_int('d',2)
%
% every thrd month
% H_xdatetick_int('m',2)
%


%see also: H_xdatetick
xtick=[];

if nargin==1
    labInt=1;
    mode='DD/mm';
    ax=gca;
end
if nargin==2
    mode='DD';
    ax=gca;
end
if nargin==3
    ax=gca;
end
if ~ischar(mode)
    mode='DD';
end

if ischar(int)
    xlim=get(ax,'xlim');
    tm=floor(xlim(1)):24/24:ceil(xlim(2));
    
    if regexp(lower(int),'d')
        %set tick and label at every midnight
        xtk=tm;
    elseif regexp(lower(int),'m')
        [~, ~, d]=datevec(tm);
        xtk=tm(d==1);
        
    elseif regexp(lower(int),'y')
        [~, m, d]=datevec(tm);
        xtk=tm(m==1 & d==1);
        
    else
        error('H_xdatetick_int: unknown input')
    end
    
    xtkLab=datestr(xtk',mode);
    
    for i=2:labInt:length(xtkLab)
        xtkLab(i:(i+labInt-2),:)=' ';
    end
    
    set(gca,'xtick',xtk,'xticklabel',xtkLab)
else
    xlim=get(ax,'xlim');
    ylim=get(ax,'ylim');
    xtick=floor(xlim(1)):int/24:ceil(xlim(2));
    set(ax,'xtick',xtick)
    set(ax,'xlim',xlim,'ylim',ylim);
    grid on
    
    
    % datetick('x',19,'keepticks','keeplimits');
    
    if int<24
        %     ticks=get(gca,'Xtick');
        %     lab=get(gca,'XTickLabel')
        datetick(ax,'x',15,'keepticks','keeplimits');
%         mode='HH:MM';
    else
        datetick(ax,'x',19,'keepticks','keeplimits');
%         mode='dd/mm';
    end
    
    
    if nargin==2
        labInt=round(labInt/int);
        xlab=get(gca,'xticklabel');
        xtick=get(gca,'xtick')';
        for i=1:labInt:length(xlab)
            xlab(i:(i+labInt-2),:)=' ';
            xtick(i:(i+labInt-2),:)=NaN;
        end
        set(gca,'xticklabel',xlab);
        xtick(isnan(xtick))=[];
        if isinf(ylim(1))
            ytick=get(gca,'ytick');
            ytmin=min(ytick);
            ytmax=max(ytick);
        else
            ytmin=ylim(1);
            ytmax=ylim(2);
        end
        
        hl=line([xtick,xtick]',[repmat(ytmin,length(xtick),1),repmat(ytmax,length(xtick),1)]','color',[0.6 0.6 0.6]);
        uistack(hl,'bottom')
    end
end