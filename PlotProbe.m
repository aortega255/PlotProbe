function varargout = PlotProbe(varargin)
% PLOTPROBE MATLAB code for PlotProbe.fig
%      PLOTPROBE, by itself, creates a new PLOTPROBE or raises the existing
%      singleton*.
%
%      H = PLOTPROBE returns the handle to a new PLOTPROBE or the handle to
%      the existing singleton*.
%
%      PLOTPROBE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTPROBE.M with the given input arguments.
%
%      PLOTPROBE('Property','Value',...) creates a new PLOTPROBE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlotProbe_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlotProbe_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlotProbe

% Last Modified by GUIDE v2.5 11-Oct-2021 15:39:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlotProbe_OpeningFcn, ...
                   'gui_OutputFcn',  @PlotProbe_OutputFcn, ...
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


% --- Executes just before PlotProbe is made visible.
function PlotProbe_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlotProbe (see VARARGIN)
debugFlag=1;

handles.debugFlag=debugFlag;
if debugFlag
%     Homer3Path='C:\Users\oanto\Documents\MATLAB\work\BU\Homer3tucker';
%     run([Homer3Path,filesep,'setpaths'])
    load('out_hrf.mat','yavg','data','probe')
    datos.timeSeries=data;
    datos.averages=yavg;
    datos.probe=probe;
else
    if exist("varargin","var")
        datos=varargin{1};
    end
end

handles.data=datos;
% Choose default command line output for PlotProbe
handles.output = hObject;

%if data is present, send to axis to plot
if ~isempty(datos)
    plotData(hObject,handles);
end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PlotProbe wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PlotProbe_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse press over axes background.
function probeAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to probeAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.currentPoint = get(hObject,'CurrentPoint');


% --------------------------------------------------------------------
function fileTag_Callback(hObject, eventdata, handles)
% hObject    handle to fileTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function loadTag_Callback(hObject, eventdata, handles)
% hObject    handle to loadTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% function called when loading the data
fName=uigetfile;
load(fNAme,'yavg','data','probe')
datos.timeSeries=data;
datos.averages=yavg;
datos.probe=probe;
handles.data=datos;
%if data is present, send to axis to plot
if ~isempty(datos)
    plotData(hObject,handles);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function plotData(hObject,handles)

probe= handles.data.probe;
yavg=handles.data.averages;
rhoRange=[str2double(handles.loLimrho.String) str2double(handles.hiLimrho.String)];
fig =handles.probeAxes;

if isempty(rhoRange)
    rhoMin = 0;
    rhoMax = 100;
else
    rhoMin = rhoRange(1);
    if length(rhoRange>=2)
        rhoMax = rhoRange(2);
    else
        rhoMax = 100;
    end
end


% Get output from BasicProcSteam.m

measurementList = yavg.measurementList;
dAvg = yavg.dataTimeSeries;
tHRF = yavg.time;

% create the ml matrix and get source detector separation (SDS) for each channel 

nMeas = length(measurementList);

mlHRF = zeros(nMeas,4);
rhoHRF = zeros(nMeas,1);
pCH = zeros(nMeas,2);

for iM = 1:nMeas
    mlHRF(iM,1) = measurementList(iM).sourceIndex;
    mlHRF(iM,2) = measurementList(iM).detectorIndex;
    mlHRF(iM,4) = measurementList(iM).wavelengthIndex;

    switch measurementList(iM).dataTypeLabel
        
        case 'HRF HbO'
            mlHRF(iM,3) = 1;
        case 'HRF HbR'
            mlHRF(iM,3) = 2;
        case 'HRF HbT'
            mlHRF(iM,3) = 3;
            
        case 'HRF dOD'
            mlHRF(iM,3) = 1;
        case 'HRF dMTF'
            mlHRF(iM,3) = 2;
        case 'HRF dVar'
            mlHRF(iM,3) = 3;
    end
    
    ps = probe.sourcePos3D( mlHRF(iM,1), : );
    pd = probe.detectorPos3D( mlHRF(iM,2), : );
    rhoHRF(iM) = norm( ps-pd);
    
    pCH(iM,:) = (probe.sourcePos2D( mlHRF(iM,1), 1:2 ) + probe.detectorPos2D( mlHRF(iM,2), 1:2 )) / 2;
end 

% plot HRFs in probe arrangement

lstKeep1 = find(rhoHRF>=rhoMin & rhoHRF<=rhoMax & mlHRF(:,3)==1 );
lstKeep2 = find(rhoHRF>=rhoMin & rhoHRF<=rhoMax & mlHRF(:,3)==2 );
lstKeep3 = find(rhoHRF>=rhoMin & rhoHRF<=rhoMax & mlHRF(:,3)==3 );

lstKeep1a = find(rhoHRF>=rhoMin & rhoHRF<=rhoMax & mlHRF(:,3)==1 & mlHRF(:,4)==1 );
lstKeep2a = find(rhoHRF>=rhoMin & rhoHRF<=rhoMax & mlHRF(:,3)==2 & mlHRF(:,4)==1 );
lstKeep3a = find(rhoHRF>=rhoMin & rhoHRF<=rhoMax & mlHRF(:,3)==3 & mlHRF(:,4)==1 );

lstKeep1b = find(rhoHRF>=rhoMin & rhoHRF<=rhoMax & mlHRF(:,3)==1 & mlHRF(:,4)==2 );
lstKeep2b = find(rhoHRF>=rhoMin & rhoHRF<=rhoMax & mlHRF(:,3)==2 & mlHRF(:,4)==2 );
lstKeep3b = find(rhoHRF>=rhoMin & rhoHRF<=rhoMax & mlHRF(:,3)==3 & mlHRF(:,4)==2 );

% d1 = dAvg(:,lstKeep1);
% d2 = dAvg(:,lstKeep2);
% d3 = dAvg(:,lstKeep3);

d1a = dAvg(:,lstKeep1a);
d2a = dAvg(:,lstKeep2a);
d3a = dAvg(:,lstKeep3a);

d1b = dAvg(:,lstKeep1b);
d2b = dAvg(:,lstKeep2b);
d3b = dAvg(:,lstKeep3b);

axes(fig)
hp=plot( probe.sourcePos2D(mlHRF(lstKeep1,1),1), probe.sourcePos2D(mlHRF(lstKeep1,1),2), 'r.' );
set(hp,'color',[1 0.5 0.5])
set(hp,'markersize',8)
hold on
hp=plot( probe.detectorPos2D(mlHRF(lstKeep1,2),1), probe.detectorPos2D(mlHRF(lstKeep1,2),2), 'b.' );
set(hp,'color',[0.5 0.5 1])
set(hp,'markersize',8)
%plot( pCH(lstKeep1,1), pCH(lstKeep1,2), 'go' )

xx = ((tHRF - min(tHRF)) / (max(tHRF) - min(tHRF)) - 0.5) * 0.05;
yy = dAvg(:,lstKeep1);
yy(:,:,2) = dAvg(:,lstKeep2);
yy(:,:,3) = dAvg(:,lstKeep3);

if 0
    minyy = min(yy(:));
    maxyy = max(yy(:));
    yy1 = ((yy(:,:,1)) / (maxyy-minyy)) * 0.05;
    yy2 = ((yy(:,:,2)) / (maxyy-minyy)) * 0.05;
    yy3 = ((yy(:,:,3)) / (maxyy-minyy)) * 0.05;
    
    for iM = 1:length(lstKeep1)
        plot( xx+pCH(lstKeep1(iM),1), yy1(:,iM)+pCH(lstKeep1(iM),2), 'r-')
    end
    for iM = 1:length(lstKeep2)
        plot( xx+pCH(lstKeep2(iM),1), yy2(:,iM)+pCH(lstKeep2(iM),2), 'b-')
    end
    for iM = 1:length(lstKeep3)
        plot( xx+pCH(lstKeep3(iM),1), yy3(:,iM)+pCH(lstKeep3(iM),2), 'c-')
    end
    
else
    minyy = min(min(d1a));
    maxyy = max(max(d1a));
    yy1 = ((d1a) / (maxyy-minyy)) * 0.05;
    for iM = 1:length(lstKeep1a)
        plot( xx+pCH(lstKeep1a(iM),1), yy1(:,iM)+pCH(lstKeep1a(iM),2), 'r-.')
    end

    minyy = min(min(d1b));
    maxyy = max(max(d1b));
    yy1 = ((d1b) / (maxyy-minyy)) * 0.05;
    for iM = 1:length(lstKeep1b)
        plot( xx+pCH(lstKeep1b(iM),1), yy1(:,iM)+pCH(lstKeep1b(iM),2), 'r-')
    end
    
    minyy = min(min(d2a));
    maxyy = max(max(d2a));
    yy1 = ((d2a) / (maxyy-minyy)) * 0.05;
    for iM = 1:length(lstKeep2a)
        plot( xx+pCH(lstKeep2a(iM),1), yy1(:,iM)+pCH(lstKeep2a(iM),2), 'b-.')
    end

    minyy = min(min(d2b));
    maxyy = max(max(d2b));
    yy1 = ((d2b) / (maxyy-minyy)) * 0.05;
    for iM = 1:length(lstKeep2b)
        plot( xx+pCH(lstKeep2b(iM),1), yy1(:,iM)+pCH(lstKeep2b(iM),2), 'b-')
    end
    
    minyy = min(min(d3a));
    maxyy = max(max(d3a));
    yy1 = ((d3a) / (maxyy-minyy)) * 0.05;
    for iM = 1:length(lstKeep3a)
        plot( xx+pCH(lstKeep3a(iM),1), yy1(:,iM)+pCH(lstKeep3a(iM),2), 'c-.')
    end

    minyy = min(min(d3b));
    maxyy = max(max(d3b));
    yy1 = ((d3b) / (maxyy-minyy)) * 0.05;
    for iM = 1:length(lstKeep3b)
        plot( xx+pCH(lstKeep3b(iM),1), yy1(:,iM)+pCH(lstKeep3b(iM),2), 'c-')
    end
end

guidata(hObject, handles);



function hiLimrho_Callback(hObject, eventdata, handles)
% hObject    handle to hiLimrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hiLimrho as text
%        str2double(get(hObject,'String')) returns contents of hiLimrho as a double


% --- Executes during object creation, after setting all properties.
function hiLimrho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hiLimrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loLimrho_Callback(hObject, eventdata, handles)
% hObject    handle to loLimrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loLimrho as text
%        str2double(get(hObject,'String')) returns contents of loLimrho as a double


% --- Executes during object creation, after setting all properties.
function loLimrho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loLimrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filterDistances.
function filterDistances_Callback(hObject, eventdata, handles)
% hObject    handle to filterDistances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filterDistances


% --- Executes on button press in debug.
function debug_Callback(hObject, eventdata, handles)
% hObject    handle to debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('debug!')
