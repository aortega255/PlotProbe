function varargout = plotProbeV2(varargin)
% PLOTPROBEV2 MATLAB code for plotProbeV2.fig
%      PLOTPROBEV2, by itself, creates a new PLOTPROBEV2 or raises the existing
%      singleton*.
%
%      H = PLOTPROBEV2 returns the handle to a new PLOTPROBEV2 or the handle to
%      the existing singleton*.
%
%      PLOTPROBEV2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTPROBEV2.M with the given input arguments.
%
%      PLOTPROBEV2('Property','Value',...) creates a new PLOTPROBEV2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plotProbeV2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plotProbeV2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plotProbeV2

% Last Modified by GUIDE v2.5 07-Dec-2021 13:26:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @plotProbeV2_OpeningFcn, ...
    'gui_OutputFcn',  @plotProbeV2_OutputFcn, ...
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


% --- Executes just before plotProbeV2 is made visible.
function plotProbeV2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plotProbeV2 (see VARARGIN)
debugFlag=0;

handles.debugFlag=debugFlag;
if debugFlag
    %     Homer3Path='C:\Users\oanto\Documents\MATLAB\work\BU\Homer3tucker';
    %     run([Homer3Path,filesep,'setpaths'])
    datos=conditionData('out_hrf.mat');
else
    if ~isempty(varargin)
        datos=conditionData(varargin);
    end
end
if exist('datos','var')
    handles.data=datos;
    %if data is present, send to axis to plot
    if ~isempty(datos)
        plotData(hObject,handles);
    end
end
% Choose default command line output for plotProbeV2
handles.output = hObject;




% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plotProbeV2 wait for user response (see UIRESUME)
% uiwait(handles.mainFigure);


% --- Outputs from this function are returned to the command line.
function varargout = plotProbeV2_OutputFcn(hObject, eventdata, handles)
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
%handles.currentPoint = get(hObject,'CurrentPoint');
stimIndex=1;  %<-change this to be a control reading from the available types of stims
%Identify curve parameters
HRFNum=str2num(get(eventdata.Source,'Tag')); %#ok<ST2NM>
wavelengthIndex=handles.data.averages.measurementList(HRFNum).wavelengthIndex;
dataTypeLabel=handles.data.averages.measurementList(HRFNum).dataTypeLabel;
dataType=handles.data.averages.measurementList(HRFNum).dataType;
dataTypeIndex=handles.data.averages.measurementList(HRFNum).dataTypeIndex;
sourceIndex=handles.data.averages.measurementList(HRFNum).sourceIndex;
detectorIndex=handles.data.averages.measurementList(HRFNum).detectorIndex;
rho=handles.data.rhoHRF(HRFNum);

%find the other related curves
ml=handles.data.measurementList;
relatedChannels=find(ml(:,1)==sourceIndex&ml(:,2)==detectorIndex); %thes are all same source and detector and thus same rho; indexed to the time series


%figure out which time series I'm plotting
activePlots=[handles.al1.Value||handles.al2.Value,...
    handles.ml1.Value||handles.ml2.Value,...
    handles.vl1.Value||handles.vl2.Value];
indiPlots=find(activePlots==1);
nsubplots=sum(activePlots);

%define time series to plot
t=handles.data.timeSeries.time;
plotBuffer=handles.data.timeSeries.dataTimeSeries(:,relatedChannels);


%plot
figure
switch handles.seriesSelector.SelectedObject.String
    case 'Delta moment'
        if nsubplots
            for ki=1:nsubplots
                subplot(nsubplots,1,ki)
                switch indiPlots(ki)
                    case 1
                        scaling=1;
                        offset=0;
                        if handles.al1.Value
                            yyaxis left
                            plot(t,plotBuffer(:,1))
                            scaling=max(plotBuffer(:,1));
                            offset=mean(plotBuffer(:,1));
                            ylabel('Intensity \lambda_1')
                            hold on
                        end
                        if handles.dispStims.Value
                            plot(t,handles.data.stims{stimIndex}.onsets*(scaling-offset)+offset,'-k')
                            Ys=[zeros(2,size(handles.data.stims{stimIndex}.s1,2))+offset;...
                                scaling*ones(2,size(handles.data.stims{stimIndex}.s1,2))];
                            patch(handles.data.stims{stimIndex}.s1,Ys,[1,.6,.6])
                        end
                        if handles.al2.Value
                            yyaxis right
                            plot(t,plotBuffer(:,4))
                            ylabel('Intensity \lambda_2')
                        end
                    case 2
                        scaling=1;
                        offset=0;
                        if handles.ml1.Value
                            yyaxis left
                            plot(t,plotBuffer(:,2))
                            scaling=max(plotBuffer(:,2));
                            offset=mean(plotBuffer(:,2));
                            ylabel('First moment \lambda_1')
                            hold on
                        end
                        if handles.dispStims.Value
                            plot(t,handles.data.stims{stimIndex}.onsets*(scaling-offset)+offset,'-k')
                            Ys=[zeros(2,size(handles.data.stims{stimIndex}.s1,2))+offset;...
                                scaling*ones(2,size(handles.data.stims{stimIndex}.s1,2))];
                            patch(handles.data.stims{stimIndex}.s1,Ys,[1,.6,.6])
                        end
                        hold off
                        if handles.ml2.Value
                            yyaxis right
                            plot(t,plotBuffer(:,5))
                            ylabel('First moment \lambda_2')
                        end
                    case 3
                        scaling=1;
                        offset=0;
                        if handles.vl1.Value
                            yyaxis left
                            plot(t,plotBuffer(:,3))
                            scaling=max(plotBuffer(:,3));
                            offset=mean(plotBuffer(:,3));
                            ylabel('Variance \lambda_1')
                            hold on
                        end
                        if handles.dispStims.Value
                            plot(t,handles.data.stims{stimIndex}.onsets*(scaling-offset)+offset,'-k')
                            Ys=[zeros(2,size(handles.data.stims{stimIndex}.s1,2))+offset;...
                                scaling*ones(2,size(handles.data.stims{stimIndex}.s1,2))];
                            patch(handles.data.stims{stimIndex}.s1,Ys,[1,.6,.6])
                            
                        end
                        hold off
                        if handles.vl2.Value
                            yyaxis right
                            plot(t,plotBuffer(:,6))
                            ylabel('Variance \lambda_2')
                        end
                end
            end
            sgtitle(['Source ', num2str(sourceIndex),...
                '. Detector ',num2str(detectorIndex),...
                '. (',' \rho=',num2str(rho),')']);
            xlabel('Time [s]')
        else
            plot(t,handles.data.timeSeries.dataTimeSeries(:,HRFNum))
            xlabel('Time [s]')
            ytags={'Intensity','m1','V'};
            ylabel([ytags{handles.data.probe.momentOrders(dataTypeIndex)+1},...
                '\lambda_',num2str(wavelengthIndex),...
                ' \rho=',num2str(rho)])
            title(['Source ', num2str(sourceIndex),...
                '. Detector ',num2str(detectorIndex),...
                '.']);
        end
    case 'Hb time series'
        
    case 'HRF'
        
end





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
[fName,path]=uigetfile;
if ~isempty(fName)
    handles.mainFigure.Name = ['plotProbeV2 - ',fName];
    datos=conditionData({[path,fName]});
    handles.data=datos;
    %if data is present, send to axis to plot
    if ~isempty(datos)
        plotData(hObject,handles);
    end
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function mainFigure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function [outputData,params]=conditionData(inputData)
%this function will be used either to process the data received either from
%the load function or as the passed argument. It will either accept a
%filename to load, or a data structure containing the desired variables. It
%will be assumed that a string input is a filename, otherwise it will be a
%structure
%This should be able to read snirf files and not only mat files
if length(inputData)==1 %assumes it's a file name
    load(inputData{1},'yavg','data','probe','stims')
    datos.timeSeries=data;
    datos.averages=yavg;
    datos.probe=probe;
    datos.stims=stims;
elseif length(inputData)>1 %assumes it's a varargin argument
    %however, said varargin argument should have up to 7 entries; non
    %provided arguments should be []
    %[probe,stim,HRFdod,HRFconc,dod,conc,Intensity]
    %dod is implicitly permitted to be delta of a moment of higher order;
    %Intensity is implicitly permitted to be a moment of higher order
    datos.probe=inputData{1};
    datos.stims=inputData{2};
    flagConc=0;
    flagDod=0;
    try
        datos.averages=inputData{3};                
        flagDod=1; %flag indicating what type of data was loaded
    catch
        datos.averages=inputData{4};
        flagConc=1;
    end
    %for now ignore HRF in conc space if it's present in dod space; at some
    %point it should do different things if we provide both
    flagTimeSeries=0;  %if this is zero, disable clicking on curves
    if flagDod
        try
            datos.timeSeries=inputData{5};
            flagTimeSeries=1;
        end
    elseif flagConc%for now, assume only one or the other is active
        try
            datos.timeSeries=inputData{6};
            flagTimeSeries=1;
        end
    else
        %if neither is present,try to use intensity data
        try
            datos.timeSeries=inputData{7};
            flagTimeSeries=1;
        end
    end            
end


%do some data and variable validation. I basically need to find if the data
%is TD fNIRS or regular. Also, find what fields are included so options can
%be enabled in the GUI

%preprocess stuff
stimClassNumber=length(datos.stims);
stims=cell(1,stimClassNumber);
stimNames=cell(1,stimClassNumber);
tt=datos.timeSeries.time;
%make stimulus time series
for kii=1:stimClassNumber
    onsets=zeros(size(tt));  %stimulus onsets
    s1=zeros(4,size(datos.stims(kii).data,1)); %used for creating patches
    for ki=1:size(datos.stims(kii).data,1)
        sIdx=find(tt>=datos.stims.data(ki,1)&tt<=datos.stims.data(ki,1)+datos.stims.data(ki,2));
        onsets(sIdx(1))=1;
        s1(:,ki)=[tt(sIdx(1));tt(sIdx(end));tt(sIdx(end));tt(sIdx(1))];
    end
    stims{kii}.onsets=onsets;
    stims{kii}.s1=s1;
    stims{kii}.name=datos.stims(kii).name;
    stimNames{kii}=datos.stims(kii).name;
end
datos.stims=stims;
[ml,rho] = getMeasurementListArray(datos.timeSeries,datos.probe);
datos.measurementList=ml;
[ml,rho] = getMeasurementListArray(datos.averages,datos.probe);
datos.measurementListHRF=ml;
datos.rhoHRF=rho;

%output a structure indicating what things will be enabled or
%disabled based on the data available. Not sure if it's better to send the
%handles structure to this function directly or make the function caller do
%the GUI adaptation
params.stims.number=stimClassNumber;
params.stims.names=stimNames;

outputData=datos;


function plotData(hObject,handles)
%data types at https://github.com/fNIRS/snirf/blob/master/snirf_specification.md
probe= handles.data.probe;
yavg=handles.data.averages;
rhoRange=[str2double(handles.loLimrho.String) str2double(handles.hiLimrho.String)];
fig =handles.probeAxes;
buttonPressFuncHandle=@(hObject,eventdata)plotProbeV2('probeAxes_ButtonDownFcn',hObject,eventdata,guidata(hObject));
cla(fig)

if isempty(rhoRange)||~handles.filterDistances.Value
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
ax=fig;
hp=plot(ax, probe.sourcePos2D(mlHRF(lstKeep1,1),1), probe.sourcePos2D(mlHRF(lstKeep1,1),2), 'r.' );
set(hp,'color',[1 0.5 0.5])
set(hp,'markersize',8)
set(hp,'ButtonDownFcn',buttonPressFuncHandle);
hold on
hp=plot(ax, probe.detectorPos2D(mlHRF(lstKeep1,2),1), probe.detectorPos2D(mlHRF(lstKeep1,2),2), 'b.' );
set(hp,'color',[0.5 0.5 1])
set(hp,'markersize',8)
set(hp,'ButtonDownFcn',buttonPressFuncHandle);
%plot( pCH(lstKeep1,1), pCH(lstKeep1,2), 'go' )

xx = ((tHRF - min(tHRF)) / (max(tHRF) - min(tHRF)) - 0.5) * 0.05;
yy = dAvg(:,lstKeep1);
yy(:,:,2) = dAvg(:,lstKeep2);
yy(:,:,3) = dAvg(:,lstKeep3);

if 0
    minyy = min(yy(:));
    maxyy = max(yy(:));
    yy1 = ((yy(:,:,1)) / (maxyy-minyy)) * 0.05 * str2double(handles.edit3.String);
    yy2 = ((yy(:,:,2)) / (maxyy-minyy)) * 0.05 * str2double(handles.edit3.String);
    yy3 = ((yy(:,:,3)) / (maxyy-minyy)) * 0.05 * str2double(handles.edit3.String);
    for iM = 1:length(lstKeep1)
        h=plot(ax, xx+pCH(lstKeep1(iM),1), yy1(:,iM)+pCH(lstKeep1(iM),2), 'r-');
        set(h,'ButtonDownFcn',buttonPressFuncHandle);
        set(h,'Tag',num2str(lstKeep1(iM)));
    end
    for iM = 1:length(lstKeep2)
        h=plot(ax, xx+pCH(lstKeep2(iM),1), yy2(:,iM)+pCH(lstKeep2(iM),2), 'b-');
        set(h,'ButtonDownFcn',buttonPressFuncHandle);
        set(h,'Tag',num2str(lstKeep2(iM)));
    end
    for iM = 1:length(lstKeep3)
        h=plot(ax, xx+pCH(lstKeep3(iM),1), yy3(:,iM)+pCH(lstKeep3(iM),2), 'c-');
        set(h,'ButtonDownFcn',buttonPressFuncHandle);
        set(h,'Tag',num2str(lstKeep3(iM)));
    end
    
else
    minyy = min(min(d1a));
    maxyy = max(max(d1a));
    yy1 = ((d1a) / (maxyy-minyy)) * 0.05 * str2double(handles.edit3.String);
    for iM = 1:length(lstKeep1a)
        h=plot(ax, xx+pCH(lstKeep1a(iM),1), yy1(:,iM)+pCH(lstKeep1a(iM),2), 'r-.');
        set(h,'ButtonDownFcn',buttonPressFuncHandle);
        set(h,'Tag',num2str(lstKeep1a(iM)));
    end
    
    minyy = min(min(d1b));
    maxyy = max(max(d1b));
    yy1 = ((d1b) / (maxyy-minyy)) * 0.05 * str2double(handles.edit3.String);
    for iM = 1:length(lstKeep1b)
        h=plot(ax, xx+pCH(lstKeep1b(iM),1), yy1(:,iM)+pCH(lstKeep1b(iM),2), 'r-');
        set(h,'ButtonDownFcn',buttonPressFuncHandle);
        set(h,'Tag',num2str(lstKeep1b(iM)));
    end
    
    minyy = min(min(d2a));
    maxyy = max(max(d2a));
    yy1 = ((d2a) / (maxyy-minyy)) * 0.05 * str2double(handles.edit3.String);
    for iM = 1:length(lstKeep2a)
        h=plot(ax, xx+pCH(lstKeep2a(iM),1), yy1(:,iM)+pCH(lstKeep2a(iM),2), 'b-.');
        set(h,'ButtonDownFcn',buttonPressFuncHandle);
        set(h,'Tag',num2str(lstKeep2a(iM)));
    end
    
    minyy = min(min(d2b));
    maxyy = max(max(d2b));
    yy1 = ((d2b) / (maxyy-minyy)) * 0.05 * str2double(handles.edit3.String);
    for iM = 1:length(lstKeep2b)
        h=plot(ax, xx+pCH(lstKeep2b(iM),1), yy1(:,iM)+pCH(lstKeep2b(iM),2), 'b-');
        set(h,'ButtonDownFcn',buttonPressFuncHandle);
        set(h,'Tag',num2str(lstKeep2b(iM)));
    end
    
    minyy = min(min(d3a));
    maxyy = max(max(d3a));
    yy1 = ((d3a) / (maxyy-minyy)) * 0.05 * str2double(handles.edit3.String);
    for iM = 1:length(lstKeep3a)
        h=plot(ax, xx+pCH(lstKeep3a(iM),1), yy1(:,iM)+pCH(lstKeep3a(iM),2), 'c-.');
        set(h,'ButtonDownFcn',buttonPressFuncHandle);
        set(h,'Tag',num2str(lstKeep3a(iM)));
    end
    
    minyy = min(min(d3b));
    maxyy = max(max(d3b));
    yy1 = ((d3b) / (maxyy-minyy)) * 0.05 * str2double(handles.edit3.String);
    for iM = 1:length(lstKeep3b)
        h=plot(ax, xx+pCH(lstKeep3b(iM),1), yy1(:,iM)+pCH(lstKeep3b(iM),2), 'c-');
        set(h,'ButtonDownFcn',buttonPressFuncHandle);
        set(h,'Tag',num2str(lstKeep3b(iM)));
    end
end
axis(ax,'image')

% Update handles structure
guidata(hObject, handles);



function hiLimrho_Callback(hObject, eventdata, handles)
% hObject    handle to hiLimrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hiLimrho as text
%        str2double(get(hObject,'String')) returns contents of hiLimrho as a double
plotData(hObject,handles)

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
plotData(hObject,handles)


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
plotData(hObject,handles)

% --- Executes on button press in debug.
function debug_Callback(hObject, eventdata, handles)
% hObject    handle to debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('debug!')


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function probeAxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to probeAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate probeAxes


% --- Executes during object deletion, before destroying properties.
function probeAxes_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to probeAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in al1.
function al1_Callback(hObject, eventdata, handles)
% hObject    handle to al1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of al1


% --- Executes on button press in ml1.
function ml1_Callback(hObject, eventdata, handles)
% hObject    handle to ml1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ml1


% --- Executes on button press in vl1.
function vl1_Callback(hObject, eventdata, handles)
% hObject    handle to vl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vl1


% --- Executes on button press in al2.
function al2_Callback(hObject, eventdata, handles)
% hObject    handle to al2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of al2


% --- Executes on button press in ml2.
function ml2_Callback(hObject, eventdata, handles)
% hObject    handle to ml2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ml2


% --- Executes on button press in vl2.
function vl2_Callback(hObject, eventdata, handles)
% hObject    handle to vl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vl2


% --- Executes on button press in dispStims.
function dispStims_Callback(hObject, eventdata, handles)
% hObject    handle to dispStims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dispStims


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

%plot data with the new scale
plotData(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
