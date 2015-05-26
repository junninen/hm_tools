function varargout = eventClassifier(varargin)
% EVENTCLASSIFIER M-file for eventClassifier.fig
%      EVENTCLASSIFIER, by itself, creates a new EVENTCLASSIFIER or raises the existing
%      singleton*.
%
%      H = EVENTCLASSIFIER returns the handle to a new EVENTCLASSIFIER or the handle to
%      the existing singleton*.
%
%      EVENTCLASSIFIER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVENTCLASSIFIER.M with the given input arguments.
%
%      EVENTCLASSIFIER('Property','Value',...) creates a new EVENTCLASSIFIER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before eventClassifier_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eventClassifier_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help eventClassifier

% Last Modified by GUIDE v2.5 04-Jun-2007 16:24:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @eventClassifier_OpeningFcn, ...
    'gui_OutputFcn',  @eventClassifier_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before eventClassifier is made visible.
function eventClassifier_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eventClassifier (see VARARGIN)

% Choose default command line output for eventClassifier
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using eventClassifier.
if strcmp(get(hObject,'Visible'),'off')
    plot(1);
    set(gca,'fontSize',12)
end

%make empty data variable
global hmD
hmD=[];

% UIWAIT makes eventClassifier wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%==========================================================================
% --- Outputs from this function are returned to the command line.
function varargout = eventClassifier_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
%% CloseMenuItem_Callback
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global currentDate
global hmD
global parPath

clear currentDate
clear hmD

if exist([fileparts(parPath),'\hmTemp'],'dir')
    disp(['remove temporary directory: ',[fileparts(parPath),'\hmTemp']])
    rmdir([fileparts(parPath),'\hmTemp'],'s')
end

selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

%==========================================================================
% --- Executes on selection change in popupmenu1.
%% popupmenu1_Callback
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

global currentDate;

popup_sel_index = get(handles.popupmenu1, 'Value');
yy=get(handles.popupmenu1, 'String');
yyy=str2num(yy{popup_sel_index});
currentDate=datenum(yyy,1,1);
loadData(handles);
makePlot(handles);

%slider in julian day
%for some reason does not work in datenum

mn=1;
mx=365;

slider_step(1) = 1/(mx-mn);
slider_step(2) = 30/(mx-mn);
set(handles.slider1,'sliderstep',slider_step,...
    'max',mx,'min',mn,'Value',mn)


%if file-name-model is given make slider
% dPath=get(handles.dataPath,'String');
% %check if dPath end to \ or /
% if ispc
%     if double(dPath(end))~=92
%         dPath=[dPath,'\'];
%     end
% else
%     if double(dPath(end))~=47
%         dPath=[dPath,'/'];
%     end
% end
%
% dats=dir([dPath,'sum\',yy]);

%==========================================================================
% --- Executes during object creation, after setting all properties.
%% popupmenu1_CreateFcn
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% function popupmenu1_populate(hObject, eventdata, handles)
% global parPath
%
% dPath=get(handles.dataPath,'String');
% fls=dir(dPath);
%
% h=1;
% folderNames{h}='Select year';
% for i=1:length(fls)
%     if fls(i).isdir & isempty(strmatch('.',fls(i).name))
%         h=h+1;
%         folderNames{h}=fls(i).name;
%     end
% end
% set(handles.popupmenu1, 'String', folderNames);
%

%==========================================================================
% --- Executes on button press in ev1a.
%% ev1a_Callback
function ev1a_Callback(hObject, eventdata, handles)
% hObject    handle to ev1a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ev1a
global hmD;
global currentDate;

handles=setRadioPutton(handles,'ev1a');
saveEvTab(handles,'ev1a')

set_sumEv(handles)

%==========================================================================
% --- Executes on button press in nonEv.
%% nonEv_Callback
function nonEv_Callback(hObject, eventdata, handles)
% hObject    handle to nonEv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nonEv
global hmD;
global currentDate;

handles=setRadioPutton(handles,'nonEv');
saveEvTab(handles,'nonEv')

set_sumEv(handles)


%==========================================================================
% --- Executes on button press in ev2.
%% ev2_Callback
function ev2_Callback(hObject, eventdata, handles)
% hObject    handle to ev2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ev2
global hmD;
global currentDate;

handles=setRadioPutton(handles,'ev2');
saveEvTab(handles,'ev2')

set_sumEv(handles)


%==========================================================================
% --- Executes on button press in ev1b.
%% ev1b_Callback
function ev1b_Callback(hObject, eventdata, handles)
% hObject    handle to ev1b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ev1b
global hmD;
global currentDate;

handles=setRadioPutton(handles,'ev1b');
saveEvTab(handles,'ev1b')

set_sumEv(handles)

%==========================================================================
%% evApple_Callback
function evApple_Callback(hObject, eventdata, handles)
% hObject    handle to undef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hmD;
global currentDate;

handles=setRadioPutton(handles,'evApple');
saveEvTab(handles,'evApple')

set_sumEv(handles)

%==========================================================================
%% evBump_Callback
function evBump_Callback(hObject, eventdata, handles)
% hObject    handle to undef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hmD;
global currentDate;

handles=setRadioPutton(handles,'evBump');
saveEvTab(handles,'evBump')

set_sumEv(handles)

%==========================================================================
%% evRain_Callback
function evRain_Callback(hObject, eventdata, handles)
% hObject    handle to undef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hmD;
global currentDate;

handles=setRadioPutton(handles,'evRain');
saveEvTab(handles,'evRain')

set_sumEv(handles)

%==========================================================================
%% evFeatureless_Callback
function evFeatureless_Callback(hObject, eventdata, handles)
% hObject    handle to undef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hmD;
global currentDate;

handles=setRadioPutton(handles,'evFeatureless');
saveEvTab(handles,'evFeatureless')

set_sumEv(handles)

%==========================================================================
% --- Executes on button press in undef.
%% undef_Callback
function undef_Callback(hObject, eventdata, handles)
% hObject    handle to undef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of undef
global hmD;
global currentDate;

handles=setRadioPutton(handles,'undef');
saveEvTab(handles,'undef')

set_sumEv(handles)

%==========================================================================
%% setRadioPutton
function handles=setRadioPutton(handles,radPutton)
%set all radioputtons to zero, but the one indicated
%

set(handles.ev1a,'Value',0)
set(handles.ev1b,'Value',0)
set(handles.ev2,'Value',0)
set(handles.evApple,'Value',0)
set(handles.evBump,'Value',0)
set(handles.evRain,'Value',0)
set(handles.evFeatureless,'Value',0)
set(handles.nonEv,'Value',0)
set(handles.undef,'Value',0)
set(handles.badData,'Value',0)

eval(['set(handles.',radPutton,',''Value'',1)'])

%==========================================================================
%% saveEvTab
function saveEvTab(handles,radPutton)
%save evTab to hmD-structure
%

global hmD

cols={'date',...
    'ev1a','ev1b','ev2','evApple','evBump','evRain','evFeatureless',...
    'nonEv','undef','badData','partlyBad','checkSum'};
evTab=zeros(1,12);
Icol=strmatch(radPutton,cols);

evTab(Icol)=1;
evTab(1)=hmD.meta.startTime;
partlyBad=get(handles.partlyBad,'Value');

evTab(strmatch('partlyBad',cols))=partlyBad;
evTab(strmatch('checkSum',cols))=sum(evTab(2:strmatch('badData',cols)));

instr_list=get(handles.popupmenu4_instruments,'String');
sel=get(handles.popupmenu4_instruments,'Value');
instr=instr_list{sel};

eval(['hmD.eventClass.',instr,'.evTabHdr=cols;']);
eval(['hmD.eventClass.',instr,'.evTab=evTab;']);


% instr_list=get(handles.popupmenu4_instruments,'String');
% sel=get(handles.popupmenu4_instruments,'Value');
% instr=instr_list{sel};
% evTab=zeros(1,12);
% partlyBad=get(handles.partlyBad,'Value');
% eval(['jDay=floor(hmD.',instr,'(2,1));']);
% eval(['hmD.eventClass.',instr,'.evTab=[hmD.meta.startTime,str2num(datestr(currentDate,10)),jDay,1,0,0,0,0,0,0,partlyBad,0];']);


%==========================================================================
% --- Executes on button press in partlyBad.
%% partlyBad_Callback
function partlyBad_Callback(hObject, eventdata, handles)
% hObject    handle to partlyBad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of partlyBad
global hmD;
global currentDate;

instr_list=get(handles.popupmenu4_instruments,'String');
sel=get(handles.popupmenu4_instruments,'Value');
instr=instr_list{sel};

partlyBad=1;
eval(['jDay=floor(hmD.',instr,'(2,1));']);
eval(['hmD.eventClass.',instr,'.evTab=[hmD.meta.startTime,str2num(datestr(currentDate,10)),jDay,0,0,0,0,0,0,0,partlyBad,0];']);

set_sumEv(handles)

%==========================================================================
% --- Executes on button press in badData.
%% badData_Callback
function badData_Callback(hObject, eventdata, handles)
% hObject    handle to badData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of badData
global hmD;
global currentDate;

handles=setRadioPutton(handles,'badData');
saveEvTab(handles,'badData')

set_sumEv(handles)
%%=============================================================
% --- Executes on button press in prev.
%% prev_Callback
function prev_Callback(hObject, eventdata, handles)
% hObject    handle to prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global currentDate

saveTempData(handles);
currentDate=currentDate-1;
slider1_update(handles);
makePlot(handles);


%==========================================================================
% --- Executes on button press in next.
%% next_Callback
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global currentDate

saveTempData(handles);
currentDate=currentDate+1;
slider1_update(handles);
makePlot(handles);

%==========================================================================
% --- Executes on slider movement.
%% slider1_Callback
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global currentDate

jday=floor(get(hObject,'Value'));
currentDate=jday+datenum(str2num(datestr(currentDate,10)),1,1)-1;
makePlot(handles);

% --- Executes on slider movement.
function slider1_update(handles)

global currentDate
jday=currentDate-datenum(str2num(datestr(currentDate,10)),1,1)+1;
set(handles.slider1,'Value',jday);

%%===============================================================
% --- Executes on button press in quit.
%% quit_Callback
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global currentDate
global hmD
global parPath

clear currentDate
clear hmD

disp(['remove temporary directory: ',[fileparts(parPath),'\hmTemp']])
rmdir([fileparts(parPath),'\hmTemp'],'s')

CloseMenuItem_Callback(hObject, eventdata, handles)


%==========================================================================
% --- Executes during object creation, after setting all properties.
%% slider1_CreateFcn
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%==========================================================================
% --- Executes on selection change in popupmenu4_instruments.
%% popupmenu4_instruments_Callback
function popupmenu4_instruments_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4_instruments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu4_instruments contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4_instruments

global hmD


instr_list=get(handles.popupmenu4_instruments,'String');
sel=get(handles.popupmenu4_instruments,'Value');
instr=instr_list{sel};

classif_populate(handles,instr);

makePlot(handles);

%==========================================================================
%% classif_populate
function classif_populate(handles,instr)

global hmD
%

if isfield(hmD,'eventClass')
    if isfield(hmD.eventClass,instr)
        try
            %try first load headers from structure
            eval(['cols=hmD.eventClass.',instr,'.evTabHdr;'])
        catch
            cols={'date',...
                'ev1a','ev1b','ev2','evApple','evBump','evRain','evFeatureless',...
                'nonEv','undef','badData','partlyBad','checkSum'};
        end
        for i=2:length(cols)-2
            val=eval(['hmD.eventClass.',instr,'.evTab(',num2str(strmatch(cols{i},cols)),')']);
            eval(['set(handles.',num2str(cols{i}),',''Value'',val)'])
        end

        %         set(handles.ev1a,'Value',eval(['hmD.eventClass.',instr,'.evTab(4)']))
    end
else
    %set all to zero
    cols={'date',...
        'ev1a','ev1b','ev2','evApple','evBump','evRain','evFeatureless',...
        'nonEv','undef','badData','partlyBad','checkSum'};

    for i=2:length(cols)-2
        val=0;
        eval(['set(handles.',num2str(cols{i}),',''Value'',val)'])
    end

    %     set(handles.ev1a,'Value',0)

end


% --- Executes during object creation, after setting all properties.
%% popupmenu4_instruments_CreateFcn
function popupmenu4_instruments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4_instruments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%========================================================================
%load, save and plot data
%% loadData
function loadData(handles)
global currentDate;
global hmD;
global parPath;
global savePath;

if ispc
    delim='\';
else
    delim='/';
end

tempFile=[fileparts(parPath),delim,'hmTemp',delim,'hm_',num2str(currentDate),'.mat'];
savedTempFile=[savePath,delim,'hm_',H_datestr(currentDate,34),'.mat'];

%if temp file exist load this
if exist(tempFile,'file')
    tmp=load(tempFile);
    hmD=tmp.hmD;
    clear tmp
elseif exist(savedTempFile,'file')
    tmp=load(savedTempFile);
    hmD=tmp.hmD;
    clear tmp
else
    hmD=hm_load(currentDate,'parPath',parPath);
end

Iinst=find(hmD.meta.insts_loaded);

for i=1:length(hmD.meta.insts_loaded)
    if hmD.meta.insts_loaded(i)
        classif_populate(handles,hmD.meta.insts{i})
    end
end

if any(hmD.meta.insts_loaded)
    insts=hmD.meta.insts(hmD.meta.insts_loaded);

    val=get(handles.popupmenu4_instruments, 'Value');
    if val>length(insts)
        set(handles.popupmenu4_instruments, 'Value',1)
    end
    set(handles.popupmenu4_instruments, 'String', insts);
    set(handles.inst_available,'string',['Instruments available ',num2str(sum(hmD.meta.insts_loaded))])

else
    set(handles.popupmenu4_instruments, 'String', {'missing'});
    set(handles.inst_available,'string',['Instruments available ',num2str(sum(hmD.meta.insts_loaded))])
end

%========================================================================
%% saveTempData
function saveTempData(handles)
%save data to temp folder that is created where parameter file is
%parPath\hmTemp

global parPath
global currentDate
global hmD

tempPath=[fileparts(parPath),'\hmTemp\'];
if ~exist(tempPath,'dir')
    mkdir(tempPath);
end

% if ~exist(['save ',resPath,'\',yy,'\hmRes\','hm',yy,mm,dd,'.mat'],'file')
eval(['save(''',tempPath,'hm_',num2str(currentDate),''',''hmD'')'])
% end


%========================================================================
%% makePlot
function makePlot(handles)

global hmD
global currentDate;


if currentDate~=hmD.meta.startTime
    loadData(handles);
end
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu4_instruments, 'Value');
instList=get(handles.popupmenu4_instruments, 'String');
inst=instList{popup_sel_index};

if isfield(hmD,inst)
    set(handles.smoothing,'Value',eval(['hmD.meta.',inst,'.smoothing']));

    hm_plot(hmD,inst);

    set(gca,'fontSize',12)
    xlim=get(gca,'xlim');
    % set(gca,'xlim',xlim,'ylim',[1e-9 1e-6])
    hold on
    plot(xlim,[25e-9,25e-9],'k')
    hold off

    title(' ')
    set(handles.sum_inst,'String',[upper(inst),': ',datestr(currentDate,1)])

    if isfield(hmD,'fit')
        if isfield(hmD.fit,inst)
            hold on
            eval(['plot(hmD.fit.',inst,'.tm,hmD.fit.',inst,'.zs{1},''.'')'])
            hline=findobj(gca,'Type','line');
            for i=1:length(hline)
                set(hline(i),'color',[0.5,0.5,0.5])
            end
            hold off
        end
    end

    set_sumEv(handles)
else
    plot([0 NaN 1],[0 NaN 1])
    set(gca,'xticklabel',[' '],'yticklabel',[' '])
    text(0.3,0.5,'No data to plot','fontsize',20)
    set(handles.sum_inst,'String',datestr(currentDate,1))

end

%remove previous day or instrument handles
set(handles.figure1,'UserData',[])
set(handles.sum_modes,'Value',1)
%=========================================================================
%% set_sumEv
function set_sumEv(handles)
%Set event calssification to summary string
%
%

clsStr={'Event 1a','Event 1b','Event 2',...
    'Apple','Bump','Rain','Featureless',...
    'Non event','Undefined','Bad data'};
cls(1)=get(handles.ev1a,'Value');
cls(2)=get(handles.ev1b,'Value');
cls(3)=get(handles.ev2,'Value');
cls(4)=get(handles.evApple,'Value');
cls(5)=get(handles.evBump,'Value');
cls(6)=get(handles.evRain,'Value');
cls(7)=get(handles.evFeatureless,'Value');
cls(8)=get(handles.nonEv,'Value');
cls(9)=get(handles.undef,'Value');
cls(10)=get(handles.badData,'Value');

if any(cls)
    for i=1:10
        if cls(i)
            set(handles.sum_ev,'String',clsStr{i})
        end
    end
else
    set(handles.sum_ev,'String',' ')
end
set_sumModes(handles)

%=======================================================================
%% set_sumModes
function set_sumModes(handles,newFitMade)
%
% set modes summary to summary panel
%
global currentTime
global hmD

if nargin ==1
    newFitMade=0;
end

popup_sel_index = get(handles.popupmenu4_instruments, 'Value');
instList=get(handles.popupmenu4_instruments, 'String');
inst=instList{popup_sel_index};

if isfield(hmD,'fit')
    if isfield(eval(['hmD.fit']),inst)
        if isfield(eval(['hmD.fit.',inst]),'gr')
            gr=eval(['hmD.fit.',inst,'.gr']);
            if isempty(gr)
                set(handles.sum_modes,'String',{' '})
                set(handles.sum_modes,'Value',1)
            else
                for i=1:length(gr)
                    dp=round(gr(i).dp_range*1e9*10)/10;
                    grAvr=num2str(round(gr(i).gr_avr(1)*1e9*10)/10);
                    grStd=num2str(round(gr(i).gr_avr(2)*1e9*100)/100);
                    mods{i}=['GR: ', grAvr,  '+-' ,grStd,'[nm/h] range: ',num2str(dp(1)),'-',num2str(dp(2)),' nm'];
                end

                set(handles.sum_modes,'String',mods)
                if newFitMade
                    set(handles.sum_modes,'Value',length(mods))
                end
            end
        else
            set(handles.sum_modes,'String',{' '})
            set(handles.sum_modes,'Value',1)
        end
    else
        set(handles.sum_modes,'String',{' '})
        set(handles.sum_modes,'Value',1)
    end
else
    set(handles.sum_modes,'String',{' '})
    set(handles.sum_modes,'Value',1)
end


% --- Executes on button press in setParameters.
%% setParameters_Callback
function setParameters_Callback(hObject, eventdata, handles)
% hObject    handle to setParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global parPath

h=hm_setParamGUI;
waitfor(h);
parseParams(hObject, eventdata, handles);

%% parseParams
function parseParams(hObject, eventdata, handles)
%parse parameters load or set from hm_setParamsGUI
global parPath
params=hm_readInstrumentParam(parPath);

%find years available
%only check if subfolder under path are years
%otherwise use just years 1996-present
h=0;
for i=1:length(params)
    if strcmp(params(i).fnm(10:12),'%4d')
        fls=dir(params(i).path);
        for ii=1:length(fls)
            if fls(ii).isdir & isempty(strfind(fls(ii).name,'.'))
                y=str2num(fls(ii).name);
                if ~isempty(y)
                    h=h+1;
                    y_list(h)=y;
                end
            end
        end
    end
    inst_list{i}=params(i).name;
end

y_list=unique(y_list);
y_cell{1}='Select year';
for i=1:length(y_list)
    y_cell{i+1}=num2str(y_list(i));
end

set(handles.popupmenu1, 'String', y_cell);
set(handles.popupmenu4_instruments, 'String', inst_list);

%=======================================================================
% --- Executes on button press in makeFitting.
%% makeFitting_Callback
function makeFitting_Callback(hObject, eventdata, handles)
% hObject    handle to makeFitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hmD

popup_sel_index = get(handles.popupmenu4_instruments, 'Value');
instList=get(handles.popupmenu4_instruments, 'String');
inst=instList{popup_sel_index};

set_info_txt(handles,'Select starting and ending time by clicking surface plot.')

axes(handles.axes1)
[x,y]=ginput(2);

axes(handles.axes2)
% eval(['dat=H_2dmedfilt(hmD.',inst,'{1}(2:end,3:end),[3,1]);']);
eval(['dat=hmD.',inst,'{1}(2:end,3:end);']);
% dat=hmD.dmps(2:end,3:end);
% eval(['dp=hmD.',inst,'(1,3:end);']);
% eval(['tim=hmD.',inst,'(2:end,1);']);
eval(['dp=hmD.meta.',inst,'.dp{1};']);
eval(['tim=hmD.meta.',inst,'.tim{1};']);
mxNrM=4;

%if fitting is already done append fittings
if isfield(hmD,'fit')
    if isfield(hmD.fit,inst)
        eval(['zs=hmD.fit.',inst,'.zs{1};']);
        eval(['ws=hmD.fit.',inst,'.ws{1};']);
        eval(['Hs=hmD.fit.',inst,'.Hs{1};']);
        eval(['Ns=hmD.fit.',inst,'.Ns{1};']);
    else
        zs=NaN(size(tim,1),mxNrM);
        ws=NaN(size(tim,1),mxNrM);
        Hs=NaN(size(tim,1),mxNrM);
        Ns=NaN(size(tim,1),mxNrM);
    end
else
    zs=NaN(size(tim,1),mxNrM);
    ws=NaN(size(tim,1),mxNrM);
    Hs=NaN(size(tim,1),mxNrM);
    Ns=NaN(size(tim,1),mxNrM);
end

[mn,I1]=min(abs(tim-x(1)));
[mn,I2]=min(abs(tim-x(2)));
for i=I1:I2
    Iok=~isnan(dat(i,:));
    if i==I1
        [param{i},fval,yhat,h_out,N]=MF_lognorm(dat(i,Iok),dp(Iok),1,4);
    else
        [param{i},fval,yhat,h_out,N]=MF_lognorm(dat(i,Iok),dp(Iok),1,4,param{i-1});
    end

%     [param,fval,yhat,H,N]=MF_lognorm(dat(i,:),dp,0,mxNrM);
    nrM=length(N);
    zs(i,1:nrM)=param{i}(nrM+1:end);
    ws(i,1:nrM)=param{i}(1:nrM);
    Ns(i,1:nrM)=N;
    Hs(i,1:nrM)=h_out;

%     axis([0.5e-9 1000e-9,-inf inf])
%     set(gca,'xtick',[1e-9,1e-8,1e-7,1e-6],'xgrid','on','xminorgrid','off')
%     set(gca,'yticklabel',' ')

%     drawnow
    axes(handles.axes1)
    hold on
    
    for m=1:nrM
%         log10(Ns(i,m))
%         markerSize=min(round(log10(Ns(i,m))+2),9)
        markerSize=round(min(max(((Ns(i,m)-10)/6000)*16,1),10));
        %     plot(handles.axes1,repmat(tim(i),1,length(param{i})/2),param{i}(length(param{i})/2+1:end),...
        %         'k.','markersize',markerSize)
        plot(handles.axes1,tim(i),zs(i,m),'ko','markersize',markerSize,'markerfacecolor','k')
    end
    hold off
    axes(handles.axes2)
end
eval(['hmD.fit.',inst,'.tm{1}=tim;']);
eval(['hmD.fit.',inst,'.zs{1}=zs;']);
eval(['hmD.fit.',inst,'.ws{1}=ws;']);
eval(['hmD.fit.',inst,'.Ns{1}=Ns;']);
eval(['hmD.fit.',inst,'.Hs{1}=Hs;']);


set_info_txt(handles,' ')

%==========================================================================
% --- Executes on button press in smoothing.
%% smoothing_Callback
function smoothing_Callback(hObject, eventdata, handles)
% hObject    handle to smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of smoothing
global hmD

popup_sel_index = get(handles.popupmenu4_instruments, 'Value');
instList=get(handles.popupmenu4_instruments, 'String');
inst=instList{popup_sel_index};

if get(hObject,'Value')
    try
        eval(['hmD=hm_smoothing(hmD,''2dmeanfilt'',[3,5],''',inst,''');']);
    catch
        disp(lasterr)
    end
else
    if eval(['hmD.meta.',inst,'.smoothing'])
        eval(['hmD.',inst,'=hmD.meta.',inst,'.raw;'])
        eval(['hmD.meta.',inst,'.smoothing=0;'])
    end
end

makePlot(handles);


%==========================================================================
% --- Executes on button press in pushbutton13_gr.
%% pushbutton13_gr_Callback
function pushbutton13_gr_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13_gr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hmD

popup_sel_index = get(handles.popupmenu4_instruments, 'Value');
instList=get(handles.popupmenu4_instruments, 'String');
inst=instList{popup_sel_index};
axes(handles.axes1)

[hmD,gr_str]=hm_gr(hmD,inst);

% axis(handles.axes2)
% plot(gr_str.)

if isempty(gr_str)
    set_info_txt(handles,'Mode fitting have to be made first!')
else
    instStr=eval(['hmD.fit.',inst]);
    if isfield(instStr,'gr')
        nrOldM=length(instStr.gr);
        instStr.gr(nrOldM+1)=gr_str;
    else
        instStr.gr=gr_str;
    end
    eval(['hmD.fit.',inst,'=instStr;']);
end

set_sumModes(handles,1)
sum_modes_Callback(handles.sum_modes, eventdata, handles)

%% set_info_txt
function set_info_txt(handles,txt)
%set info text
%

set(handles.info_text,'String',txt)

% --------------------------------------------------------------------
%% loadParam_Callback
function loadParam_Callback(hObject, eventdata, handles)
% % hObject    handle to loadParam (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)

[fileName, pathName]=uigetfile('*.txt','Pick a parameter file');

if fileName
    global parPath

    parPath=[pathName,fileName];
    parseParams(hObject, eventdata, handles)
end

% % --------------------------------------------------------------------
% function setParam_Callback(hObject, eventdata, handles)
% % hObject    handle to setParam (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
%% save_Callback
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%copy temp folder to given path,
%if not temp folder rename the files
global parPath
global savePath

savePath=uigetdir(fileparts(parPath),'Select the folder to save work to');

if ispc
    delim='\';
else
    delim='/';
end
nrFls=0;
fls=dir([fileparts(parPath),delim,'hmTemp']);
for f=1:length(fls)
    if ~fls(f).isdir & strmatch('hm_',fls(f).name)
        try
            nrFls=nrFls+1;
            dn=[H_datestr(str2num(fls(f).name(4:end-3)),34),'.mat'];
            srs=[fileparts(parPath),delim,'hmTemp',delim,fls(f).name];
            dest=[savePath,delim,'hm_',dn];
            [success(nrFls),message,messageID] = copyfile(srs,dest);
            if ~success(nrFls)
                disp(['eventClassifier: save_Callback: ' message])
            end
        catch
        end
    end
end

%if all ok delete temp folder
if all(success)
    %copy also the param file
    [pth,nm,ext]=fileparts(parPath);
    [success,message,messageID] = copyfile(parPath,[savePath,delim,nm,ext]);
    if ~strcmp(parPath,[savePath,delim,nm,ext])
        disp(['eventClassifier: parameters file is copied to ',savePath])
        disp(['new param file (',[savePath,delim,nm,ext],') will be used'])
    end
    rmdir([fileparts(parPath),'\hmTemp'],'s')
    parPath=[savePath,delim,nm,ext];
end

%==========================================================================
% --- Executes on selection change in sum_modes.
%% sum_modes_Callback
function sum_modes_Callback(hObject, eventdata, handles)
% hObject    handle to sum_modes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns sum_modes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sum_modes

global hmD
set_sumModes(handles)

popup_sel_index = get(handles.popupmenu4_instruments, 'Value');
instList=get(handles.popupmenu4_instruments, 'String');
inst=instList{popup_sel_index};

popup_sel_index = get(handles.sum_modes, 'Value');
modesList=get(handles.sum_modes, 'String');
modes=modesList{popup_sel_index};

eval(['tm=hmD.fit.',inst,'.gr(',num2str(popup_sel_index),').tm;']);
eval(['z=hmD.fit.',inst,'.gr(',num2str(popup_sel_index),').z;']);
eval(['zhat=hmD.fit.',inst,'.gr(',num2str(popup_sel_index),').zhat;']);
eval(['N=hmD.fit.',inst,'.gr(',num2str(popup_sel_index),').N;']);

plotGR(handles,tm,z,zhat,N)

%==========================================================================
%% plotGR
function plotGR(handles,tm,z,zhat,N)
h=get(handles.figure1,'UserData');
if ~isempty(h)
    for i=1:length(h)
        delete(h(i))
    end
    set(handles.figure1,'UserData',[])
end

axes(handles.axes1)
hold on
htemp=plot(tm,z,'w*',tm,zhat,'b');
% plot(tm,zhat,'b','linewidth',1.5)
hold off
set(handles.figure1,'UserData',htemp);

axes(handles.axes2)
plot(tm,N,'linewidth',1.5)
axis([-inf inf -inf inf])
set_sumModes(handles)



%==========================================================================
% --- Executes during object creation, after setting all properties.
%% sum_modes_CreateFcn
function sum_modes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sum_modes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
%% DelMode_Callback
function DelMode_Callback(hObject, eventdata, handles)
% hObject    handle to DelMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hmD
set_sumModes(handles)

popup_sel_index = get(handles.popupmenu4_instruments, 'Value');
instList=get(handles.popupmenu4_instruments, 'String');
inst=instList{popup_sel_index};

popup_sel_index = get(handles.sum_modes, 'Value');
modesList=get(handles.sum_modes, 'String');
modes=modesList{popup_sel_index};
set(handles.sum_modes, 'Value',1);
eval(['hmD.fit.',inst,'.gr(',num2str(popup_sel_index),')=[];']);
set_sumModes(handles)

h=get(handles.figure1,'UserData');
if ~isempty(h)
    for i=1:length(h)
        delete(h(i))
    end
    set(handles.figure1,'UserData',[])
end


%==========================================================================
% --- Executes on button press in delFitting.
%% delFitting_Callback
function delFitting_Callback(hObject, eventdata, handles)
% hObject    handle to delFitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hmD

popup_sel_index = get(handles.popupmenu4_instruments, 'Value');
instList=get(handles.popupmenu4_instruments, 'String');
inst=instList{popup_sel_index};

if isfield(hmD,'fit')
    if isfield(hmD.fit,inst)
        fts=rmfield(hmD.fit,inst);
        hmD.fit=fts;
    end
end

saveTempData(handles);
makePlot(handles);

%========================================================================
%% showDate
function showDate(hObject, eventdata, handles)
global currentDate

jday=floor(get(hObject,'Value'));
currentDate=jday+datenum(str2num(datestr(currentDate,10)),1,1)-1;
disp(datestr(currentDate))


% --------------------------------------------------------------------
%% loadWork_Callback
function loadWork_Callback(hObject, eventdata, handles)
% hObject    handle to loadWork (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global parPath
[fileName, pathName]=uigetfile('*.txt','Pick a parameter file of saved work');

if fileName
    global parPath

    parPath=[pathName,fileName];
    parseParams(hObject, eventdata, handles)


    fls=dir(fileparts(parPath));
    nrD=0;
    for d=1:length(fls)
        if strmatch('hm_',fls(d).name) & ~isempty(findstr('.mat',fls(d).name))
            nrD=nrD+1;
        end
    end
    set_info_txt(handles,['Number of days processed in saved folder: ',num2str(nrD)])
end