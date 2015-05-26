function varargout = hm_setParamGUI(varargin)
% HM_SETPARAMGUI M-file for hm_setParamGUI.fig
%      HM_SETPARAMGUI, by itself, creates a new HM_SETPARAMGUI or raises the existing
%      singleton*.
%
%      H = HM_SETPARAMGUI returns the handle to a new HM_SETPARAMGUI or the handle to
%      the existing singleton*.
%
%      HM_SETPARAMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HM_SETPARAMGUI.M with the given input arguments.
%
%      HM_SETPARAMGUI('Property','Value',...) creates a new HM_SETPARAMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before hm_setParamGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to hm_setParamGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help hm_setParamGUI

% Last Modified by GUIDE v2.5 16-May-2007 02:25:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @hm_setParamGUI_OpeningFcn, ...
    'gui_OutputFcn',  @hm_setParamGUI_OutputFcn, ...
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


% --- Executes just before hm_setParamGUI is made visible.
function hm_setParamGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hm_setParamGUI (see VARARGIN)

% Choose default command line output for hm_setParamGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes hm_setParamGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = hm_setParamGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function path1_Callback(hObject, eventdata, handles)
% hObject    handle to path1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path1 as text
%        str2double(get(hObject,'String')) returns contents of path1 as a double

pathName=uigetdir;
set(handles.path1,'String',pathName)

% --- Executes during object creation, after setting all properties.
function path1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
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



function fnm1_Callback(hObject, eventdata, handles)
% hObject    handle to fnm1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fnm1 as text
%        str2double(get(hObject,'String')) returns contents of fnm1 as a double


% --- Executes during object creation, after setting all properties.
function fnm1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fnm1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function path2_Callback(hObject, eventdata, handles)
% hObject    handle to path2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path2 as text
%        str2double(get(hObject,'String')) returns contents of path2 as a double

pathName=uigetdir;
set(hObject,'String',pathName)


% --- Executes during object creation, after setting all properties.
function path2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function fnm2_Callback(hObject, eventdata, handles)
% hObject    handle to fnm2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fnm2 as text
%        str2double(get(hObject,'String')) returns contents of fnm2 as a double


% --- Executes during object creation, after setting all properties.
function fnm2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fnm2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function path3_Callback(hObject, eventdata, handles)
% hObject    handle to path3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path3 as text
%        str2double(get(hObject,'String')) returns contents of path3 as a double

pathName=uigetdir;
set(hObject,'String',pathName)

% --- Executes during object creation, after setting all properties.
function path3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function fnm3_Callback(hObject, eventdata, handles)
% hObject    handle to fnm3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fnm3 as text
%        str2double(get(hObject,'String')) returns contents of fnm3 as a double


% --- Executes during object creation, after setting all properties.
function fnm3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fnm3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function path4_Callback(hObject, eventdata, handles)
% hObject    handle to path4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path4 as text
%        str2double(get(hObject,'String')) returns contents of path4 as a double

pathName=uigetdir;
set(hObject,'String',pathName)

% --- Executes during object creation, after setting all properties.
function path4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function fnm4_Callback(hObject, eventdata, handles)
% hObject    handle to fnm4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fnm4 as text
%        str2double(get(hObject,'String')) returns contents of fnm4 as a double


% --- Executes during object creation, after setting all properties.
function fnm4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fnm4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in setParam.
function setParam_Callback(hObject, eventdata, handles)
% hObject    handle to setParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% save the parameter file to hard disk and copy the path of it to global
%for usage later

global parPath

[instrs,paths,fnms]=extractContent(handles);
if isempty(fnms)
    return
end

[fileName,pathName]=uiputfile('*.txt','Save parameter file', 'hm_param.txt');

parPath=[pathName,fileName];
fid=fopen([pathName,fileName],'w');
if ispc,kono='\';else,kono='/';end

for i=1:length(instrs)
    fprintf(fid,'<instrument>\n');
    fprintf(fid,['  <name>',instrs{i},'</name>\n']);

    %replace all / with //
    Ikono=findstr(kono,paths{i});
    Ikono=[1,Ikono];
    newPath=[];
    for p=1:length(Ikono)-1
        newPath=[newPath,paths{i}(Ikono(p):Ikono(p+1))];
    end
    newPath=[newPath,paths{i}(Ikono(p+1):end),kono,kono];
    fprintf(fid,['  <path>',newPath,'</path>\n']);

    %extract fnm
    Ikono=findstr(kono,fnms{i});
    Ikono=[0,Ikono];
    if Ikono~=0
    newfnm=[];
    for p=1:length(Ikono)-1
        part{p}=fnms{i}(Ikono(p)+1:Ikono(p+1)-1);
        newfnm=[newfnm,part{p},kono,kono,kono,kono];
    end
    newfnm=[newfnm,fnms{i}(Ikono(p+1)+1:end)];
    else
        newfnm=fnms{i};
    end  
    subsT={'dd','yyyy','yy','mm'};
    subsW={'%%02d','%%4d','%%02d','%%02d'};
    subsX={'xx','xxxx','xx','xx'};
    vars=[];
    %     newfnm=fnms{i};

    %find positions
    dummyfnm=newfnm;
    possOrdered=[];
    h=0;
    for s=1:length(subsT),
        [dummyfnm,poss{s}]=hm_replace_char(dummyfnm,subsT{s},subsX{s});
        if ~isempty(poss{s})
            for si=1:length(poss{s})
                h=h+1;
                subsTordered{h}=subsT{s};
                possOrdered(h)=poss{s}(si);
            end
        end
    end

    %replace text
    for s=1:length(subsT),
        [newfnm,dummy]=hm_replace_char(newfnm,subsT{s},subsW{s});
    end

    
    [s,Is]=sort(possOrdered);
    vars=[];
    for v=1:length(Is)
        vars=[vars,',',subsTordered{Is(v)}];
    end

    fprintf(fid,['  <fnm>sprintf(''',newfnm,'''',vars,')</fnm>\n']);
    fprintf(fid,'</instrument>\n');
end
fclose(fid);

delete(handles.figure1)

% % --- Executes on button press in loadParam.
% function loadParam_Callback(hObject, eventdata, handles)
% % hObject    handle to loadParam (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% [fileName, pathName]=uigetfile('*.txt','Pick a parameter file');
% global parPath
% 
% parPath=[pathName,fileName];
% % delete(handles.figure1)
% 
% --- Executes on button press in addInstrument1.
function addInstrument1_Callback(hObject, eventdata, handles)
% hObject    handle to addInstrument1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.path2,'Visible','on');
set(handles.fnm2,'Visible','on');
set(handles.text5,'Visible','on');
set(handles.popupmenu4,'Visible','on');
set(handles.addInstrument1,'Visible','off');
set(handles.addInstrument2,'Visible','on');

% --- Executes on button press in addInstrument2.
function addInstrument2_Callback(hObject, eventdata, handles)
% hObject    handle to addInstrument2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.path3,'Visible','on');
set(handles.fnm3,'Visible','on');
set(handles.text6,'Visible','on');
set(handles.popupmenu5,'Visible','on');
set(handles.addInstrument2,'Visible','off');
set(handles.addInstrument3,'Visible','on');



% --- Executes on button press in addInstrument3.
function addInstrument3_Callback(hObject, eventdata, handles)
% hObject    handle to addInstrument3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.path4,'Visible','on');
set(handles.fnm4,'Visible','on');
set(handles.text7,'Visible','on');
set(handles.popupmenu6,'Visible','on');
set(handles.addInstrument3,'Visible','off');

function [instrs,paths,fnms]=extractContent(handles)
%Extract content from the figure
%
%


noInst=1;
h=0;
%1st instrument
if get(handles.popupmenu1,'Visible') & ~strcmp(get(handles.path1,'String'),'Part of the path that stays constant') & ~strcmp(get(handles.fnm1,'String'),'part of the path and file name that can change')
    h=h+1;
    ind=get(handles.popupmenu1,'Value');
    instList=get(handles.popupmenu1,'String');
    instrs{h}=instList{ind};

    paths{h}=get(handles.path1,'String');
    fnms{h}=get(handles.fnm1,'String');
    noInst=0;
end

%2st instrument
if get(handles.popupmenu4,'Visible') & ~strcmp(get(handles.path2,'String'),'Part of the path that stays constant') & ~strcmp(get(handles.fnm2,'String'),'part of the path and file name that can change')
    h=h+1;
    ind=get(handles.popupmenu4,'Value');
    instList=get(handles.popupmenu4,'String');
    instrs{h}=instList{ind};

    paths{h}=get(handles.path2,'String');
    fnms{h}=get(handles.fnm2,'String');
    noInst=0;
end

%3st instrument
if get(handles.popupmenu5,'Visible') & ~strcmp(get(handles.path3,'String'),'Part of the path that stays constant') & ~strcmp(get(handles.fnm3,'String'),'part of the path and file name that can change')
    h=h+1;
    ind=get(handles.popupmenu5,'Value');
    instList=get(handles.popupmenu5,'String');
    instrs{h}=instList{ind};

    paths{h}=get(handles.path3,'String');
    fnms{h}=get(handles.fnm3,'String');
    noInst=0;
end

%4st instrument
if get(handles.popupmenu6,'Visible') & ~strcmp(get(handles.path4,'String'),'Part of the path that stays constant') & ~strcmp(get(handles.fnm4,'String'),'part of the path and file name that can change')
    h=h+1;
    ind=get(handles.popupmenu6,'Value');
    instList=get(handles.popupmenu6,'String');
    instrs{h}=instList{ind};

    paths{h}=get(handles.path4,'String');
    fnms{h}=get(handles.fnm4,'String');
    noInst=0;
end

%find if / and \ are correct
for i=1:length(paths)
    if ispc,noKono='/';else,noKono='\';end
    if ~isempty(strfind(fnms{i},noKono))
        eval(['set(handles.fnm',num2str(i),',''BackgroundColor'',[1 0 0])'])
        msgbox(['No good! You are using wrong character ',noKono])
        fnms=[];
        uiwait
    else
        eval(['set(handles.fnm',num2str(i),',''BackgroundColor'',[1 1 1])'])
    end
end
if noInst
    msgbox('No good! You should choose at least one instrument...')
end
