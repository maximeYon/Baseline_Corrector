function varargout = Baseline_Corrector(varargin)
% BASELINE_CORRECTOR MATLAB code for Baseline_Corrector.fig
%      BASELINE_CORRECTOR, by itself, creates a new BASELINE_CORRECTOR or raises the existing
%      singleton*.
%
%      H = BASELINE_CORRECTOR returns the handle to a new BASELINE_CORRECTOR or the handle to
%      the existing singleton*.
%
%      BASELINE_CORRECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BASELINE_CORRECTOR.M with the given input arguments.
%
%      BASELINE_CORRECTOR('Property','Value',...) creates a new BASELINE_CORRECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Baseline_Corrector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Baseline_Corrector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Baseline_Corrector

% Last Modified by GUIDE v2.5 04-Nov-2020 16:55:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Baseline_Corrector_OpeningFcn, ...
    'gui_OutputFcn',  @Baseline_Corrector_OutputFcn, ...
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


% --- Executes just before Baseline_Corrector is made visible.
function Baseline_Corrector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Baseline_Corrector (see VARARGIN)


% Choose default command line output for NUSprocess
handles.output = hObject;
result = struct;
result.handles=handles;

% Update handles structure
guidata(hObject, handles);

result.argument = 0;
if size(varargin,1) ~=0
    result.location = varargin{1,1};
    if size(varargin,2)>1
        for mot = 2:size(varargin,2)
            result.location = [result.location ' ' varargin{1,mot}];
        end
    end
    result.argument = 1;
end

set(result.handles.AdvanceFilter,'visible','off');

% %Create a zoom handle and define the post-callback process
% handles.zhndl = zoom;
% handles.zhndl.ActionPostCallback = {@zoomMapAspect,handles};

% Store all data in the Gui
setMyData(result);

if result.argument == 1
    argument()
end



% --- Outputs from this function are returned to the command line.
function varargout = Baseline_Corrector_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function setMyData(result)
% Store data struct in figure
setappdata(gcf,'result',result);

function result=getMyData()
% Get data struct stored in figure
result=getappdata(gcf,'result');

function argument()
result=getMyData(); % Get data struct stored in figure

cla(result.handles.TemporalAxe);
cla(result.handles.AffichageSpectre);
cla(result.handles.Affichage_res);
axes(result.handles.AffichageSpectre);
setMyData(result);
result=getMyData();


%% Ouverture du fichier

result.dfile = result.location;

fid=fopen(result.dfile,'r','l');
[data,~]=fread(fid,'int32');
fclose(fid);

result.dfileImag = [result.dfile(1:end-2) '1i'];
fid=fopen(result.dfileImag,'r','l');
[dataImag,~]=fread(fid,'int32');
fclose(fid);
result.spectreImag = dataImag;

%%Ouverture acqus
dirname = result.dfile(1:end-2);
fname = result.dfile(end-2:end);
index = find(isletter(dirname));
sizeindex = size(isletter(dirname),2);
numberProc = sizeindex -index(1,end)-2;
dirnameEXPNO = dirname(1:end-(numberProc+7));
if exist([dirnameEXPNO,'acqus'],'file')==2
    head0=textread([dirnameEXPNO,'acqus'],'%s');
    n = 2:33;
    position = strmatch('$$',head0');
    acqu.file       = char(head0(position(2)+1));
    
    result.DEvalue  = str2num(char(head0(strmatch('##$DE=',head0')+1)));
    result.SWH  = str2num(char(head0(strmatch('##$SW_h=',head0')+1)));
    result.DWvalue = (1/(result.SWH*2))*10^6;
    
    set(result.handles.DW,'String',result.DWvalue);
    set(result.handles.DE,'String',result.DEvalue);
    result.TruncOnly =1;
    set(result.handles.Trunc,'Value',1);
    result.frequSinus = 3.1416*result.DEvalue/(result.DWvalue*(size(data,1)));
    set(result.handles.frequSinus,'String',result.frequSinus);
end

%%Ouverture acqs (old names)
if exist([dirname,fname(1:end-3),'.aqs'],'file')==2
    head0=textread([dirname,fname(1:end-3),'.aqs'],'%s');
    n = 2:33;
    position = strmatch('$$',head0');
    acqu.file       = char(head0(position(2)+1));
    
    result.DEvalue  = str2num(char(head0(strmatch('##$DE=',head0')+1)));
    result.SWH  = str2num(char(head0(strmatch('##$SW_h=',head0')+1)));
    result.DWvalue = (1/(result.SWH*2))*10^6;
    
    set(result.handles.DW,'String',result.DWvalue);
    set(result.handles.DE,'String',result.DEvalue);
    result.TruncOnly =1;
    set(result.handles.Trunc,'Value',1);
    result.frequSinus = 3.1416*result.DEvalue/(result.DWvalue*(size(data,1)));
    set(result.handles.frequSinus,'String',result.frequSinus);
end

result.autoDEDW = 0;
if isfield(result, 'DEvalue')==1 && isfield(result, 'DWvalue')==1;
    result.autoDEDW = 1;
end

% definition de spectre et frq
result.spectre = data;
frq = (1:size(result.spectre,1))';
meanfrq = mean(frq);
result.frq = frq-meanfrq;

%% Affichage
result.spectreF1 = result.spectre;
result.spectreFinal = result.spectreF1;
result.FitFinal = zeros(size(result.spectreFinal,1),1);
result.valfitSinus= zeros(size(result.spectreFinal,1),1);

axes(result.handles.AffichageSpectre);
cla(result.handles.AffichageSpectre);
hold on
plot(result.frq,result.spectreFinal,'Color',[0 0 1],'Parent',result.handles.AffichageSpectre);
plot(result.frq,result.FitFinal,'Color',[1 0 0],'Parent',result.handles.AffichageSpectre);
hold off
set(result.handles.AffichageSpectre,'xlim',[min(result.frq) max(result.frq)]);
set(result.handles.AffichageSpectre,'ylim',[min(result.spectreFinal)-max(result.spectreFinal)/15 max(result.spectreFinal)]);
set(result.handles.AffichageSpectre,'XtickLabel',[],'YtickLabel',[]);

cla(result.handles.Affichage_res);
set(result.handles.Affichage_res,'visible','off');

setMyData(result);
DE_Callback
DW_Callback
result=getMyData(); % Get data struct stored in figure

%%%%%%%%%%% truncation only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PointCut = ceil(result.DEvalue/result.DWvalue);
HalfPointCut = round(PointCut/2);

FIDSpectre = ifftshift(ifft(fftshift(result.spectre)));
FIDSpectreSum = abs(real(FIDSpectre))+abs(imag(FIDSpectre));

ZoneInf = FIDSpectreSum((size(FIDSpectre,1)/2+1-HalfPointCut)-6:(size(FIDSpectre,1)/2+1),1);
ZoneSup = FIDSpectreSum((size(FIDSpectre,1)/2+1):(size(FIDSpectre,1)/2+1+(HalfPointCut)+6),1);
ZoneInf = flip(ZoneInf);

DiffZoneInf = diff(ZoneInf);
if result.autoDEDW == 1
    LOCinfneg=find(DiffZoneInf(HalfPointCut:end, 1)<0);
    if isempty(LOCinfneg)
        LOCinfneg =1;
    end
    LOCinfneg = LOCinfneg+HalfPointCut;
else
    LOCinfneg=find(DiffZoneInf<0);
    if isempty(LOCinfneg)
        LOCinfneg =1;
    end
end
ZoneInfGradPos = DiffZoneInf(1:LOCinfneg(1,1)-1,1);
if size(ZoneInfGradPos,1)==0
    ZoneInfGradPos =DiffZoneInf(1,1);
end
[Valmax,PosInf] = max(ZoneInfGradPos);
if size(ZoneInfGradPos,1)> PosInf && ZoneInfGradPos(PosInf+1,1)>Valmax/2;
    PosInf = PosInf+1;
end

DiffZoneSup = diff(ZoneSup);
if result.autoDEDW == 1
    LOCSupneg=find(DiffZoneSup(HalfPointCut:end, 1)<0);
    if isempty(LOCSupneg)
        LOCSupneg =1;
    end
    LOCSupneg = LOCSupneg+HalfPointCut;
else
    LOCSupneg=find(DiffZoneSup<0);
    if isempty(LOCSupneg)
        LOCSupneg =1;
    end
end
ZoneSupGradPos = DiffZoneSup(1:LOCSupneg(1,1)-1,1);
if size(ZoneSupGradPos,1)==0
    ZoneSupGradPos =DiffZoneSup(1,1);
end
[ValmaxSup,PosSup] = max(ZoneSupGradPos);
if size(ZoneSupGradPos,1)> PosSup && ZoneSupGradPos(PosSup+1,1)>ValmaxSup/2;
    PosSup = PosSup+1;
end

HalfPointCut = max([PosInf PosSup]);
PointCut = 2*HalfPointCut;

FIDtrunc = (complex(ones(size(FIDSpectre,1),1),ones(size(FIDSpectre,1),1)));
FIDtrunc((size(FIDSpectre,1)/2+1-HalfPointCut+1):(size(FIDSpectre,1)/2+1+HalfPointCut-1),1) = complex(zeros(PointCut-1,1),zeros(PointCut-1,1));
FIDtrunc = FIDtrunc.*max(real([FIDSpectreSum((size(FIDSpectre,1)/2+1-HalfPointCut),1) FIDSpectreSum((size(FIDSpectre,1)/2+1+HalfPointCut),1)]));

set(result.handles.TemporalAxe,'visible','on');
cla(result.handles.TemporalAxe);
axes(result.handles.TemporalAxe);
hold on
plot(real(FIDSpectreSum),'b','Parent',result.handles.TemporalAxe);
plot(real(FIDtrunc),'r','Parent',result.handles.TemporalAxe);
set(result.handles.TemporalAxe,'xlim',[(size(FIDSpectre,1)/2+1-HalfPointCut)-10 (size(FIDSpectre,1)/2+1+HalfPointCut+10)]);
set(result.handles.TemporalAxe,'XtickLabel',[],'YtickLabel',[]);
hold off
axes(result.handles.AffichageSpectre);
result.HalfPointCut = HalfPointCut;
result.PointCut = PointCut;
set(result.handles.NumTruncPoint,'String',result.PointCut);

setMyData(result);
frequSinus_Callback
ModFit2_Callback
numbpoints_Callback
PowerStop_Callback
AdvanceFilter_Callback
Trunc_Callback

% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
axes(result.handles.AffichageSpectre);

%% ecriture du fichier 1r
% mkdir(result.dfile)
NewRfile = result.dfile;
NewR=fopen(NewRfile,'w+','l');
fwrite(NewR,result.spectreFinal,'int32')
fclose(NewR);

%% ecriture du fichier 1i
%Calcul de la partie imaginaire par Hilbert transform
if isfield(result, 'spectreImag')==1
    n = length(result.spectreImag);
    f = fft(result.spectreImag);
    h = [0 1i*ones(1,n/2-1) 0 -1i*ones(1,n/2-1)]';
    hf = h.*f;
    spectreImagShift = real(ifft(hf));
    spectreImagShiftcorr = spectreImagShift - result.FitFinal;
    f2 = fft(spectreImagShiftcorr);
    h = [0 1i*ones(1,n/2-1) 0 -1i*ones(1,n/2-1)]';
    hf2 = h.*f2;
    result.spectreImagFinal = -real(ifft(hf2));
    
    NewIfile = [result.dfile(1:end-2) '1i'];
    NewI=fopen(NewIfile,'w+','l');
    fwrite(NewI,result.spectreImagFinal,'int32')
    fclose(NewI);
end

setMyData(result)



function ModFit2_Callback(hObject, eventdata, handles)
% hObject    handle to ModFit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if isfield(result, 'ModFit2')==0
    result.ModFit2 = ceil(log10(size(result.spectre,1)*size(result.spectre,1)));
    set(result.handles.ModFit2,'String',result.ModFit2);
end
result.ModFit2 = str2double(get(result.handles.ModFit2,'String'));
setMyData(result)
% Hints: get(hObject,'String') returns contents of ModFit2 as text
%        str2double(get(hObject,'String')) returns contents of ModFit2 as a double


% --- Executes during object creation, after setting all properties.
function ModFit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ModFit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RetryFit2.
function RetryFit2_Callback(hObject, eventdata, handles)
% hObject    handle to RetryFit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get data struct stored in figure
axes(result.handles.AffichageSpectre);
clear result.spectreFinal;
setMyData(result);
second_fit()



function frequSinus_Callback(hObject, eventdata, handles)
% hObject    handle to frequSinus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if isfield(result, 'frequSinus')==0
    result.frequSinus = 1/size(result.spectre,1)*100;
    set(result.handles.frequSinus, 'ForegroundColor',[1 0 0]);
    set(result.handles.frequSinus,'String',result.frequSinus);
else
    set(result.handles.frequSinus, 'ForegroundColor',[0 0 0]);
end
result.frequSinus = str2double(get(result.handles.frequSinus,'String'));
setMyData(result)
% Hints: get(hObject,'String') returns contents of frequSinus as text
%        str2double(get(hObject,'String')) returns contents of frequSinus as a double


% --- Executes during object creation, after setting all properties.
function frequSinus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frequSinus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in retrySinusFit.
function retrySinusFit_Callback(hObject, eventdata, handles)
% hObject    handle to retrySinusFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get data struct stored in figure
clear result.valfitSinus;
clear result.spectreF1;
setMyData(result);
sinus_fit()


% --- Executes on button press in OpenFile.
function OpenFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

result=getMyData(); % Get data struct stored in figure

cla(result.handles.TemporalAxe);
cla(result.handles.AffichageSpectre);
cla(result.handles.Affichage_res);
axes(result.handles.AffichageSpectre);
setMyData(result);
result=getMyData();

clear result;

handles.output = hObject;
result = struct;
result.handles=handles;
% Update handles structure
guidata(hObject, handles);
% Store all data in the Gui
setMyData(result);

result=getMyData(); % Get data struct stored in figure

%% Ouverture du fichier
[fname,dirname] = uigetfile('*1r', 'Open Bruker 1r file');
result.dfile = [dirname fname];
fid=fopen(result.dfile,'r','l');
[data,~]=fread(fid,'int32');
fclose(fid);

result.dfileImag = [result.dfile(1:end-2) '1i'];
fid=fopen(result.dfileImag,'r','l');
[dataImag,~]=fread(fid,'int32');
fclose(fid);
result.spectreImag = dataImag;

%%Ouverture acqus
index = find(isletter(dirname));
sizeindex = size(isletter(dirname),2);
numberProc = sizeindex -index(1,end)-2;
dirnameEXPNO = dirname(1:end-(numberProc+7));
if exist([dirnameEXPNO,'acqus'],'file')==2
    %     if exist([dirname,'acqus'],'file')==2
    head0=textread([dirnameEXPNO,'acqus'],'%s');
    n = 2:33;
    position = strmatch('$$',head0');
    acqu.file       = char(head0(position(2)+1));
    
    result.DEvalue  = str2num(char(head0(strmatch('##$DE=',head0')+1)));
    result.SWH  = str2num(char(head0(strmatch('##$SW_h=',head0')+1)));
    result.DWvalue = (1/(result.SWH*2))*10^6;
    
    set(result.handles.DW,'String',result.DWvalue);
    set(result.handles.DE,'String',result.DEvalue);
    result.TruncOnly =1;
    set(result.handles.Trunc,'Value',1);
    result.frequSinus = 3.1416*result.DEvalue/(result.DWvalue*(size(data,1)));
    set(result.handles.frequSinus,'String',result.frequSinus);
end

%%Ouverture acqs (old names)
if exist([dirname,fname(1:end-3),'.aqs'],'file')==2
    head0=textread([dirname,fname(1:end-3),'.aqs'],'%s');
    n = 2:33;
    position = strmatch('$$',head0');
    acqu.file       = char(head0(position(2)+1));
    
    result.DEvalue  = str2num(char(head0(strmatch('##$DE=',head0')+1)));
    result.SWH  = str2num(char(head0(strmatch('##$SW_h=',head0')+1)));
    result.DWvalue = (1/(result.SWH*2))*10^6;
    
    set(result.handles.DW,'String',result.DWvalue);
    set(result.handles.DE,'String',result.DEvalue);
    result.TruncOnly =1;
    set(result.handles.Trunc,'Value',1);
    result.frequSinus = 3.1416*result.DEvalue/(result.DWvalue*(size(data,1)));
    set(result.handles.frequSinus,'String',result.frequSinus);
end

result.autoDEDW = 0;
if isfield(result, 'DEvalue')==1 && isfield(result, 'DWvalue')==1
    result.autoDEDW = 1;
end

% definition de spectre et frq
result.spectre = data;
frq = (1:size(result.spectre,1))';
meanfrq = mean(frq);
result.frq = frq-meanfrq;

%% Affichage
result.spectreF1 = result.spectre;
result.spectreFinal = result.spectreF1;
result.FitFinal = zeros(size(result.spectreFinal,1),1);
result.valfitSinus= zeros(size(result.spectreFinal,1),1);

axes(result.handles.AffichageSpectre);
cla(result.handles.AffichageSpectre);
hold on
plot(result.frq,result.spectreFinal,'Color',[0 0 1],'Parent',result.handles.AffichageSpectre);
plot(result.frq,result.FitFinal,'Color',[1 0 0],'Parent',result.handles.AffichageSpectre);
hold off
set(result.handles.AffichageSpectre,'xlim',[min(result.frq) max(result.frq)]);
set(result.handles.AffichageSpectre,'ylim',[min(result.spectreFinal)-max(result.spectreFinal)/15 max(result.spectreFinal)]);
set(result.handles.AffichageSpectre,'XtickLabel',[],'YtickLabel',[]);

cla(result.handles.Affichage_res);
set(result.handles.Affichage_res,'visible','off');

setMyData(result);
DE_Callback
DW_Callback
result=getMyData(); % Get data struct stored in figure


%%%%%%%%%%% truncation only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PointCut = ceil(result.DEvalue/result.DWvalue);
HalfPointCut = round(PointCut/2);

FIDSpectre = ifftshift(ifft(fftshift(result.spectre)));
FIDSpectreSum = abs(real(FIDSpectre))+abs(imag(FIDSpectre));

ZoneInf = FIDSpectreSum((size(FIDSpectre,1)/2+1-HalfPointCut)-6:(size(FIDSpectre,1)/2+1),1);
ZoneSup = FIDSpectreSum((size(FIDSpectre,1)/2+1):(size(FIDSpectre,1)/2+1+(HalfPointCut)+6),1);
ZoneInf = flip(ZoneInf);

DiffZoneInf = diff(ZoneInf);
if result.autoDEDW == 1
    LOCinfneg=find(DiffZoneInf(HalfPointCut:end, 1)<0);
    if isempty(LOCinfneg)
        LOCinfneg =1;
    end
    LOCinfneg = LOCinfneg+HalfPointCut;
else
    LOCinfneg=find(DiffZoneInf<0);
    if isempty(LOCinfneg)
        LOCinfneg =1;
    end
end
ZoneInfGradPos = DiffZoneInf(1:LOCinfneg(1,1)-1,1);
if size(ZoneInfGradPos,1)==0
    ZoneInfGradPos =DiffZoneInf(1,1);
end
[Valmax,PosInf] = max(ZoneInfGradPos);
if size(ZoneInfGradPos,1)> PosInf && ZoneInfGradPos(PosInf+1,1)>Valmax/2;
    PosInf = PosInf+1;
end

DiffZoneSup = diff(ZoneSup);
if result.autoDEDW == 1
    LOCSupneg=find(DiffZoneSup(HalfPointCut:end, 1)<0);
    if isempty(LOCSupneg)
        LOCSupneg =1;
    end
    LOCSupneg = LOCSupneg+HalfPointCut;
else
    LOCSupneg=find(DiffZoneSup<0);
    if isempty(LOCSupneg)
        LOCSupneg =1;
    end
end
ZoneSupGradPos = DiffZoneSup(1:LOCSupneg(1,1)-1,1);
if size(ZoneSupGradPos,1)==0
    ZoneSupGradPos =DiffZoneSup(1,1);
end
[ValmaxSup,PosSup] = max(ZoneSupGradPos);
if size(ZoneSupGradPos,1)> PosSup && ZoneSupGradPos(PosSup+1,1)>ValmaxSup/2;
    PosSup = PosSup+1;
end

HalfPointCut = max([PosInf PosSup]);
PointCut = 2*HalfPointCut;

FIDtrunc = (complex(ones(size(FIDSpectre,1),1),ones(size(FIDSpectre,1),1)));
FIDtrunc((size(FIDSpectre,1)/2+1-HalfPointCut+1):(size(FIDSpectre,1)/2+1+HalfPointCut-1),1) = complex(zeros(PointCut-1,1),zeros(PointCut-1,1));
FIDtrunc = FIDtrunc.*max(real([FIDSpectreSum((size(FIDSpectre,1)/2+1-HalfPointCut),1) FIDSpectreSum((size(FIDSpectre,1)/2+1+HalfPointCut),1)]));

set(result.handles.TemporalAxe,'visible','on');
cla(result.handles.TemporalAxe);
axes(result.handles.TemporalAxe);
hold on
plot(real(FIDSpectreSum),'b','Parent',result.handles.TemporalAxe);
plot(real(FIDtrunc),'r','Parent',result.handles.TemporalAxe);
set(result.handles.TemporalAxe,'xlim',[(size(FIDSpectre,1)/2+1-HalfPointCut)-10 (size(FIDSpectre,1)/2+1+HalfPointCut+10)]);
set(result.handles.TemporalAxe,'XtickLabel',[],'YtickLabel',[]);
hold off
axes(result.handles.AffichageSpectre);
result.HalfPointCut = HalfPointCut;
result.PointCut = PointCut;
set(result.handles.NumTruncPoint,'String',result.PointCut);

setMyData(result);
frequSinus_Callback
ModFit2_Callback
numbpoints_Callback
PowerStop_Callback
AdvanceFilter_Callback
Trunc_Callback



function sinus_fit()
result=getMyData(); % Get data struct stored in figure
axes(result.handles.AffichageSpectre);
spectre = result.spectre;
frq = result.frq;
cla(result.handles.TemporalAxe);
set(result.handles.TemporalAxe,'visible','off');

h = msgbox({'Please wait' 'Fitting in progress'});

MINTAB = [frq  spectre];
% PeakFinder_Callback;
% result.PeakFinder = get(result.handles.PeakFinder, 'Value');
% if result.PeakFinder == 1
%     % find peaks
%     [pks,locs,widths,proms] = findpeaks(spectre);
%     [sortproms, indproms] = sort(proms,'descend');
%     maxssb = 100;
%     if length(proms)>maxssb
%         locs = locs(indproms(1:maxssb));
%         widths = widths(indproms(1:maxssb));
%     end
%     locs_whole = [];
%     fw = 1.0;
%     for i = 1:length(locs)
%         locs_whole = [locs_whole , locs(i)-fw*widths(i):locs(i)+fw*widths(i)];
%     end
%     locs_whole = unique(round(locs_whole));
%     locs_whole(locs_whole<=0)=[];
%     locs_whole(locs_whole>length(spectre))=[];
%     MINTAB(locs_whole,:) = [];
% end

% filtre des pics intenses
bins = size(MINTAB,1)/10;
[hi xhi] = hist(MINTAB(:,2),bins);
hi = smooth(hi,8)';
[maxH,coordH] = max(hi);
NoiseH = hi(hi>= maxH/2);
LimValpic = xhi(min([coordH + round(size(NoiseH,2)*1.0) length(xhi)]));
LimValpic = LimValpic *result.AdvanceFilterValue/100;
Mintab1 = MINTAB(:,1);
Mintab2 = MINTAB(:,2);
frqF = Mintab1(Mintab2< LimValpic);
spectreF = Mintab2(Mintab2< LimValpic);


%% premier fit sinus card -------------------------------------------
% diminution taille
if size(spectreF,1)> 10000;
    sizefactor = size(spectreF,1)/10000;
    spectreFS = (interp1(1:size(spectreF,1),spectreF,1:sizefactor:size(spectreF,1)))';
    frqFS = (interp1(1:size(spectreF,1),frqF,1:sizefactor:size(spectreF,1)))';
else
    spectreFS = spectreF;
    frqFS = frqF;
end
%Paramètres
per = result.frequSinus; % facteur multiplicatif period premier fit

amplfit = max(spectre(:))/2;
halfspectre = find(spectre==max(spectre))-size(spectre,1)/2;
%halfspectre = 0;
offset = size(spectre,1)*0.10;

myfittype = fittype('(-a*(sin((x-b)*c))/((x-b)*c))','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c'});
myfit = fit(frqFS,spectreFS,myfittype ,'StartPoint', [-amplfit, halfspectre, per],'Lower', [-amplfit*10, -offset, 0.00001],'Upper', [amplfit*10, offset, 1],'Robust','LAR','MaxFunEvals',1000,'TolFun',10^-8,'TolX',10^-10);
result.valfitSinus = feval(myfit,frq);

%%%%% test
set(result.handles.frequSinus,'String',myfit.c);

%%%%%%%%%%% truncation only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HalfPointCut = result.HalfPointCut;
PointCut = result.PointCut;

FidBasline = ifftshift(ifft(fftshift(result.valfitSinus)));
FIDBasTruncOnly = complex(zeros(size(FidBasline,1),1),zeros(size(FidBasline,1),1));
%FIDtrunc((size(FIDSpectre,1)/2+1-HalfPointCut+1):(size(FIDSpectre,1)/2+1+HalfPointCut-1),1) = complex(zeros(PointCut-1,1),zeros(PointCut-1,1));
FIDBasTruncOnly((size(FidBasline,1)/2+1-HalfPointCut+1:(size(FidBasline,1)/2+1+HalfPointCut-1)),1) = FidBasline((size(FidBasline,1)/2+1-HalfPointCut)+1:(size(FidBasline,1)/2+1+HalfPointCut-1),1);
FitSinusTruncOnly = fftshift(fft(fftshift(FIDBasTruncOnly)));

if result.TruncOnly ==1;
    result.valfitSinus = real(FitSinusTruncOnly);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(h)

cla(result.handles.AffichageSpectre);
plot(frq,spectre,'Color',[0 0 1],'Parent',result.handles.AffichageSpectre);
hold on
plot(frq,result.valfitSinus,'Color',[1 0 0],'Parent',result.handles.AffichageSpectre);
hold off
set(result.handles.AffichageSpectre,'xlim',[min(frq) max(frq)]);
set(result.handles.AffichageSpectre,'XtickLabel',[],'YtickLabel',[]);
result.spectreF1 = spectre-result.valfitSinus;
result.spectreFinal = result.spectreF1;
result.FitFinal = result.valfitSinus;
% handles    structure with handles and user data (see GUIDATA)

cla(result.handles.Affichage_res);
set(result.handles.Affichage_res,'visible','on');
plot(frq,result.spectreFinal,'Color',[0 0 0],'Parent',result.handles.Affichage_res);
set(result.handles.Affichage_res,'xlim',[min(frq) max(frq)]);
set(result.handles.Affichage_res,'ylim',[min(result.spectreFinal)-max(result.spectreFinal)/15 max(result.spectreFinal)/5]);
set(result.handles.Affichage_res,'XtickLabel',[],'YtickLabel',[]);

linkaxes([result.handles.AffichageSpectre,result.handles.Affichage_res],'x')

setMyData(result)


function second_fit()

result=getMyData(); % Get data struct stored in figure
axes(result.handles.AffichageSpectre);
valfit = result.valfitSinus;
spectreF1 = result.spectreF1;
spectre = result.spectre;
frq = result.frq;

spectrecorr = spectreF1;
critere = 1;
valfittot = zeros(size(spectre,1),2);
valfittot(:,1) = valfit;
FitParam = result.ModFit2;

global true;
true = 1;


%%%%%%%%%%% truncation only %%%%%%%%%%%%%%%%%%%%
HalfPointCut = result.HalfPointCut;
PointCut = result.PointCut;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iter = 1;
result.iter = iter;
set(result.handles.iterations,'String',num2str(iter));
LimValOld = 0;
frqOld = frq;
PeakFinder_Callback;
result.PeakFinder = get(result.handles.PeakFinder, 'Value');
result.SymThresh = get(result.handles.SymThresh, 'Value');
result.SW = get(result.handles.PeakFinder, 'Value');
checkbox6_Callback;
result.ShowHist = get(result.handles.checkbox6, 'Value');
while critere > (1*10^-result.PowerStopValue) && true == 1;
    %     iter
    result.iter=iter;
    set(result.handles.iterations,'String',num2str(iter));    iter = iter + 1;
    FitParam = result.ModFit2;
    
    MINTAB = [frq  spectrecorr];
    if result.PeakFinder == 1
        % find peaks
        [pks,locs,widths,proms] = findpeaks(spectrecorr);
        [sortproms, indproms] = sort(proms,'descend');
        maxssb = 100;
        if length(proms)>maxssb
            locs = locs(indproms(1:maxssb));
            widths = widths(indproms(1:maxssb));
        end
        locs(widths>length(spectrecorr)/100)=[];
        widths (widths>length(spectrecorr)/100)=[];
        %widths(widths>length(spectrecorr)/100) = length(spectrecorr)/100;
        locs_whole = [];
        fw = 2.0;
        for i = 1:length(locs)
            locs_whole = [locs_whole , locs(i)-fw*widths(i):locs(i)+fw*widths(i)];
        end
        locs_whole = unique(round(locs_whole));
        locs_whole(locs_whole<=0)=[];
        locs_whole(locs_whole>length(spectrecorr))=[];
        MINTAB(locs_whole,:) = [];
    end
    
    
    bins = size(MINTAB,1)/10;
    [hi xhi] = hist(MINTAB(:,2),bins);
    %    [hi xhi] = hist(spectrecorr,bins);
    hi = smooth(hi,8)';
    [maxH,coordH] = max(hi);
    NoiseH = hi(hi>= maxH/2);
    LimVal = xhi(min([coordH + round(size(NoiseH,2)*1.0); length(hi)]));
    LimVal = LimVal *result.AdvanceFilterValue/100;
    if result.SymThresh
        LimValDown = xhi(max([coordH - size(NoiseH,2)*1; 1]));
        LimValDown = LimValDown *result.AdvanceFilterValue/100;
    end
    if (abs(LimVal-LimValOld)/LimVal)<0.05 && critere<1e-4% && iter>100% to avoid unstability
        LimVal = LimValOld;
        Mintab1 = MINTAB(:,1);
        Mintab1 = Mintab1(ismember(Mintab1,frqOld));
        Mintab2 = spectrecorr(ismember(frq,Mintab1));
        clear MINTABf;
        MINTABf(:,1) = Mintab1;
        MINTABf(:,2) = Mintab2;
        frqOld = MINTABf(:,1); %needed for next iter in this case
    else
        LimValOld = LimVal;
        Mintab1 = MINTAB(:,1);
        Mintab2 = MINTAB(:,2);
        clear MINTABf;
        if result.SymThresh
            MINTABf(:,1) = Mintab1(Mintab2< LimVal & Mintab2> LimValDown);
            MINTABf(:,2) = Mintab2(Mintab2< LimVal & Mintab2> LimValDown);
        else
            MINTABf(:,1) = Mintab1(Mintab2< LimVal);
            MINTABf(:,2) = Mintab2(Mintab2< LimVal);
        end
        frqOld = MINTABf(:,1); %needed for next iter in the other case
    end
    
    npt = 2000;
    if size(MINTABf,1)>npt
        num0 = pi*1e-6; %use NaN too long for smooth; unique number instead
        spectre0 = ones(length(spectrecorr),1)*num0;
        spectre0(ismember(frq,MINTABf(:,1)))=MINTABf(:,2);
        spectreS = smooth(spectre0,round(size(MINTABf,1)/npt));
        spectreS(spectre0 == num0) = num0;
        MINTABf(:,2) = spectreS(spectreS~=num0);
        MINTABf = MINTABf(1:round(size(MINTABf,1)/npt):end,:);
    end
    
    
    cla(result.handles.TemporalAxe);
    if result.ShowHist ==1;
        set(result.handles.TemporalAxe,'visible','on');
        limH = zeros(size(hi));
        if result.SymThresh
            limH(xhi>=LimVal | xhi<=LimValDown) = hi(xhi>=LimVal | xhi<=LimValDown);
        else
            limH(xhi>=LimVal) = hi(xhi>=LimVal);
        end
        hifilt = hi;
        hifilt(limH~=0) = 0;
        Hbar = [hifilt;limH]';
        bar(result.handles.TemporalAxe,(xhi-xhi(1))/(xhi(end)-xhi(1)),Hbar,'stacked')
        set(result.handles.TemporalAxe,'xlim',([0,min([(xhi(coordH+1)-xhi(1))*5,xhi(end)])/xhi(end)]))
        set(result.handles.TemporalAxe,'Yscale','log')
        set(result.handles.TemporalAxe,'Xscale','linear')
        set(result.handles.TemporalAxe,'XtickLabel',[],'YtickLabel',[]);
    else
        set(result.handles.TemporalAxe,'visible','off');
    end
    
    
    %%%------------ weighting for negative points------------------------
    
    
    weight = MINTABf(:,2);
    if result.SW == 1
        
        %Gaussian broadening
        gauss_creation = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
        gauss_b = (gauss_creation(1:size(spectrecorr,1),(size(spectrecorr,1)/2)+0.5,25,1,0))';
        fidgauss =fftshift(ifft(spectrecorr)).*gauss_b;
        %         fidgauss =fftshift(ifft(spectrecorr)).*gausswin(size(spectrecorr,1),size(spectrecorr,1)/50);
        spectregauss = fft(fftshift(fidgauss));
        weight0 = abs(spectregauss);
        weight0 = weight0/max(weight0)*9 + 1;
        weight=weight0(ismember(frq,MINTABf(:,1)));
        
        %         weight(weight>=0) = 0;
        %         weight = abs(weight)./max(abs(weight(:)));
        %         weight = weight./mean(weight(:));
        %         %weight = 1;
        %         %-------------- gausian weighting of the weight---------------------
        %         x = MINTABf(:,1);
        %         mu = 0;
        %         mu = sum(frq.*spectrecorr)./sum(spectrecorr);
        %         if (mu>frq(end) || mu<frq(1))
        %             mu = 0  ;
        %         end
        %         Imax = max(spectrecorr);
        %         IUpHalf = find(spectrecorr>= Imax/5);
        %
        %         s = (IUpHalf(end)-IUpHalf(1))*10;
        %         %s = size(spectrecorr,1)/12;
        %
        %         p1 = -.5 * ((x - mu)/s) .^ 2;
        %         p2 = (s * sqrt(2*pi));
        %         fgauss = exp(p1) ./ p2;
        %         fgauss = fgauss./max(fgauss(:));
        %         weight = weight.*fgauss;
        %         %-------------- gausian weighting of the weight---------------------
        %         weight = weight+1;
        %         weight(1,1) = 1;weight(end,1) = 1;
    else
        weight(:,1) = 1;
    end
    %
    weightdis = (weight-1)/(max(weight(:)))*(max(spectrecorr)-min(spectrecorr))/10 + min(spectrecorr);
    
    %%%------------ weighting for negative points------------------------
    
    %     [myfit,~] = fit(MINTABf(:,1),MINTABf(:,2),'smoothingspline','SmoothingParam',1/(10^FitParam));
    
    %     %%%------------ divide spectrum in domains if size >2000 pts --------
    %
    %     if size(MINTABf,1) > 1000
    %         clear frqD valfitD;
    %         nDomains = fix(size(MINTABf,1)/500);
    %         for i = 1:nDomains-1
    %             MINTABfD = MINTABf((i-1)*500+1:(i-1)*500+1000,:);
    %             weightD = weight((i-1)*500+1:(i-1)*500+1000);
    %             [myfit,~] = fit(MINTABfD(:,1),MINTABfD(:,2),'smoothingspline','SmoothingParam',1/(10^FitParam),'Weight',weightD);
    %             if i==1
    %                 frqD{1,1} =  frq(1:find(frq==MINTABfD(end,1)));
    %             else
    %                 frqD{i,1} = frq(find(frq==MINTABfD(1,1)):find(frq==MINTABfD(end,1)));
    %             end
    %             valfitD{i,1} = feval(myfit,frqD{i,1});
    %         end
    %         MINTABfD = MINTABf((nDomains-1)*500+1:end,:);
    %         weightD = weight((nDomains-1)*500+1:end);
    %         [myfit,~] = fit(MINTABfD(:,1),MINTABfD(:,2),'smoothingspline','SmoothingParam',1/(10^FitParam),'Weight',weightD);
    %         frqD{nDomains,1} = frq(find(frq==MINTABfD(1,1)):end);
    %         valfitD{nDomains,1} = feval(myfit,frqD{nDomains,1});
    %
    %         frqC = [];
    %         valfitC = [];
    %         for i = 1:nDomains
    %             frqC = [frqC ; frqD{i,1}];
    %             valfitC = [valfitC ; valfitD{i,1}];
    %         end
    %         [frqUni, ia, ic] = unique(frqC);
    %         valfit = valfitC(ia);
    %     else
    [myfit,~] = fit(MINTABf(:,1),MINTABf(:,2),'smoothingspline','SmoothingParam',1/(10^FitParam),'Weight',weight);
    valfit = feval(myfit,frq);
    %     end
    
    
    %     [myfit,~] = fit(MINTABf(:,1),MINTABf(:,2),'smoothingspline','SmoothingParam',1/(10^FitParam),'Weight',weight);
    %     valfit = feval(myfit,frq);
    
    %%%%%%%%%%% truncation only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FidBasline2 = ifftshift(ifft(fftshift(valfit)));
    FIDTruncOnly = complex(zeros(size(FidBasline2,1),1),zeros(size(FidBasline2,1),1));
    FIDTruncOnly((size(FidBasline2,1)/2+1-HalfPointCut)+1:(size(FidBasline2,1)/2+1+HalfPointCut-1),1) = FidBasline2((size(FidBasline2,1)/2+1-HalfPointCut)+1:(size(FidBasline2,1)/2+1+HalfPointCut-1),1);
    FitTruncOnly = fftshift(fft(fftshift(FIDTruncOnly)));
    
    if result.TruncOnly ==1;
        valfit = real(FitTruncOnly);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cla(result.handles.AffichageSpectre);
    hold on
    plot(frq,spectrecorr,'Color',[0 0 1],'Parent',result.handles.AffichageSpectre);
    plot(frq,valfit,'Color',[1 0 0],'Parent',result.handles.AffichageSpectre);
    plot(MINTABf(:,1),MINTABf(:,2),'Color',[1 0 0],'Marker','.','linestyle','none','Parent',result.handles.AffichageSpectre)
    %     plot(MINTABf(:,1),weightdis,'Color',[0 1 0],'Marker','.','linestyle','none','Parent',result.handles.AffichageSpectre)
    hold off
    set(result.handles.AffichageSpectre,'ylim',[min(spectrecorr) max(spectrecorr)/5]);
    set(result.handles.AffichageSpectre,'xlim',[min(frq) max(frq)]);
    set(result.handles.AffichageSpectre,'XtickLabel',[],'YtickLabel',[]);
    drawnow
    
    spectrecorr = spectrecorr-valfit;
    valfittot(:,2) = valfittot(:,2) + valfit;
    
    critere = (sum(abs(valfit),1))./(size(spectre,1).*max(spectre(:)));
    result.critere = critere;
    crit(iter-1) = critere;
    setMyData(result);
    result.Fitstop = get(result.handles.Stopfit, 'Value');
    result.PeakFinder = get(result.handles.PeakFinder, 'Value');
    result.SW = get(result.handles.SmartW, 'Value');
    critere_Callback;
    result.AdvanceFilterValue = str2double(get(result.handles.AdvanceFilter,'String'));
    AdvanceFilter_Callback;
    result.SymThresh = get(result.handles.SymThresh, 'Value');
    result.PowerStopValue = str2double(get(result.handles.PowerStop,'String'));
    result.ModFit2 = str2double(get(result.handles.ModFit2,'String'));
end
%result.spectreF1 = spectrecorr;
result.iter = iter;

%% plot final
% disp(['Stopping baseline correction after ' num2str(iter) ' iterations'])
valfittotsum = sum(valfittot,2);
result.spectreFinal = spectrecorr;
result.FitFinal = valfittotsum;
minspectrefit = min(min(spectre),min(valfittotsum));
cla(result.handles.AffichageSpectre);
hold on
plot(frq,spectre,'Color',[0 0 1],'Parent',result.handles.AffichageSpectre);
plot(frq,valfittotsum,'Color',[1 0 0],'Parent',result.handles.AffichageSpectre);
hold off
set(result.handles.AffichageSpectre,'xlim',[min(frq) max(frq)]);
set(result.handles.AffichageSpectre,'ylim',[min(minspectrefit) max(spectre)/5]);
set(result.handles.AffichageSpectre,'XtickLabel',[],'YtickLabel',[]);

cla(result.handles.TemporalAxe);
set(result.handles.TemporalAxe,'visible','off');
cla(result.handles.Affichage_res);
set(result.handles.Affichage_res,'visible','on');
plot(frq,result.spectreFinal,'Color',[0 0 0],'Parent',result.handles.Affichage_res);
set(result.handles.Affichage_res,'xlim',[min(frq) max(frq)]);
set(result.handles.Affichage_res,'ylim',[min(result.spectreFinal)-max(result.spectreFinal)/15 max(result.spectreFinal)/1]);
set(result.handles.Affichage_res,'XtickLabel',[],'YtickLabel',[]);

% figure(100)
% semilogy(crit)
% title('Criterion convergence')
% xlabel('Iteration')
linkaxes([result.handles.AffichageSpectre,result.handles.Affichage_res],'x')
setMyData(result);
logfile;



function critere_Callback(hObject, eventdata, handles)
% hObject    handle to critere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
set(result.handles.critere,'String',result.critere);
setMyData(result)
% Hints: get(hObject,'String') returns contents of critere as text
%        str2double(get(hObject,'String')) returns contents of critere as a double


% --- Executes during object creation, after setting all properties.
function critere_CreateFcn(hObject, eventdata, handles)
% hObject    handle to critere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Stopfit.
function Stopfit_Callback(hObject, eventdata, handles)
% hObject    handle to Stopfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global true;
true = 0;


% --- Executes on button press in PointMod.
function PointMod_Callback(hObject, eventdata, handles)
% hObject    handle to PointMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
axes(result.handles.AffichageSpectre);
FitFinal = result.FitFinal;
spectre = result.spectre;
frq = result.frq;
spectreFinal = spectre-FitFinal;
nbpoint = result.numbpoints;
minspectrefit = min(min(spectre),min(FitFinal));
cla(result.handles.TemporalAxe);
set(result.handles.TemporalAxe,'visible','off');


% affichage
result.limX = get(result.handles.AffichageSpectre,'xlim');
result.limY = get(result.handles.AffichageSpectre,'ylim');

cla(result.handles.AffichageSpectre);
hold on
plot(frq,spectre,'Color',[0 0 1],'Parent',result.handles.AffichageSpectre);
plot(frq,FitFinal,'Color',[0 0 0],'Parent',result.handles.AffichageSpectre);
hold off
set(result.handles.AffichageSpectre,'xlim',result.limX);
set(result.handles.AffichageSpectre,'ylim',result.limY);
set(result.handles.AffichageSpectre,'XtickLabel',[],'YtickLabel',[]);

% definition des points et affichage

pas = floor(size(spectre,1)/nbpoint);
Indicepos = 1:pas:size(spectre,1);
Indicepos = Indicepos(1,1:nbpoint);

names = cell(1,nbpoint+1);
for ind = 1:nbpoint
    names{1,ind} = ['p' num2str(ind)];
    result.points.(names{ind}) = impoint(result.handles.AffichageSpectre,frq(Indicepos(1,ind),1),FitFinal(Indicepos(1,ind),1));
    setColor(result.points.(names{ind}),[1 0.55 0])
end
names{1,nbpoint+1} = ['p' num2str(nbpoint+1)];
result.points.(names{nbpoint+1}) = impoint(result.handles.AffichageSpectre,frq(end,1),FitFinal(end,1));
setColor(result.points.(names{nbpoint+1}),[1 0.55 0])

pos = zeros(nbpoint+1,2);
for ind = 1:nbpoint+1
    pos(ind,:) = getPosition(result.points.(names{1,ind}));
end

result.pos = pos;
result.nbpoint = nbpoint+1;
result.names = names;
[myfitPoint,~] = fit(pos(:,1),pos(:,2),'smoothingspline','SmoothingParam',1);
valfitPoint = feval(myfitPoint,frq);
result.FitFinal = valfitPoint;

hold on
result.fi = plot(frq,valfitPoint,'Color',[1 0 0],'Parent',result.handles.AffichageSpectre);
uistack(result.fi,'bottom')
hold off

if isfield(result.handles.Affichage_res,'ylim') ==1;
    result.limYres = get(result.handles.Affichage_res,'ylim');
else
    result.limYres = [min(result.spectreFinal)-max(result.spectreFinal)/15 max(result.spectreFinal)/5];
end

if isfield(result.handles.Affichage_res,'xlim') ==1;
    result.limXres = get(result.handles.Affichage_res,'xlim');
else
    result.limXres = result.limX;
end

cla(result.handles.Affichage_res);
set(result.handles.Affichage_res,'visible','on');
plot(frq,spectreFinal,'Color',[0 0 0],'Parent',result.handles.Affichage_res);
set(result.handles.Affichage_res,'xlim',result.limXres);
set(result.handles.Affichage_res,'ylim',result.limYres);
set(result.handles.Affichage_res,'XtickLabel',[],'YtickLabel',[]);

linkaxes([result.handles.AffichageSpectre,result.handles.Affichage_res],'x')

setMyData(result)
result.t = timer;

result.t.TimerFcn = @(~,~) refreshFit();
result.t.Period = 2;
result.t.ExecutionMode = 'fixedDelay';
start(result.t)
setMyData(result)


% --- Executes on button press in ApplyPointmod.
function ApplyPointmod_Callback(hObject, eventdata, handles)
% hObject    handle to ApplyPointmod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
axes(result.handles.AffichageSpectre);
delete(result.t)
frq = result.frq;
spectre = result.spectre;
posActual = zeros(result.nbpoint,2);
for ind = 1:result.nbpoint
    posActual(ind,:) = getPosition(result.points.(result.names{1,ind}));
end

pos = posActual;
[myfitPoint,~] = fit(pos(:,1),pos(:,2),'smoothingspline','SmoothingParam',1);
valfitPoint = feval(myfitPoint,frq);

%%%%%%%%%%% truncation only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HalfPointCut = result.HalfPointCut;
PointCut = result.PointCut;

FidfitPoint = ifftshift(ifft(fftshift(valfitPoint)));
FIDPointTruncOnly = complex(zeros(size(FidfitPoint,1),1),zeros(size(FidfitPoint,1),1));
FIDPointTruncOnly((size(FidfitPoint,1)/2+1-HalfPointCut)+1:(size(FidfitPoint,1)/2+1+HalfPointCut-1),1) = FidfitPoint((size(FidfitPoint,1)/2+1-HalfPointCut)+1:(size(FidfitPoint,1)/2+1+HalfPointCut-1),1);
FitPointTruncOnly = fftshift(fft(fftshift(FIDPointTruncOnly)));

if result.TruncOnly ==1
    valfitPoint = real(FitPointTruncOnly);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result.FitFinal = valfitPoint;
result.spectreFinal = spectre-result.FitFinal;
minspectrefit = min(min(spectre),min(result.FitFinal));
% affichage
cla(result.handles.AffichageSpectre);
hold on
plot(frq,spectre,'Color',[0 0 1],'Parent',result.handles.AffichageSpectre);
plot(frq,result.FitFinal,'Color',[1 0 0],'Parent',result.handles.AffichageSpectre);
hold off
set(result.handles.AffichageSpectre,'xlim',[min(frq) max(frq)]);
set(result.handles.AffichageSpectre,'ylim',[min(minspectrefit) max(spectre)/5]);
set(result.handles.AffichageSpectre,'XtickLabel',[],'YtickLabel',[]);

cla(result.handles.Affichage_res);
result.limX = get(result.handles.Affichage_res,'xlim');
result.limY = get(result.handles.Affichage_res,'ylim');
plot(frq,result.spectreFinal,'Color',[0 0 0],'Parent',result.handles.Affichage_res);
set(result.handles.Affichage_res,'xlim',result.limX);
set(result.handles.Affichage_res,'ylim',result.limY);
set(result.handles.Affichage_res,'XtickLabel',[],'YtickLabel',[]);

linkaxes([result.handles.AffichageSpectre,result.handles.Affichage_res],'x')

setMyData(result)

function refreshFit()
result=getMyData(); % Get back the result structure
frq = result.frq;
spectre = result.spectre;

%%%%%%%%%%% truncation only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HalfPointCut = result.HalfPointCut;
PointCut = result.PointCut;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(result,'pos') ==1
    posActual = zeros(result.nbpoint,2);
    for ind = 1:result.nbpoint
        posActual(ind,:) = getPosition(result.points.(result.names{1,ind}));
    end
    if isequal(result.pos,posActual) ==0
        pos = posActual;
        [myfitPoint,~] = fit(pos(:,1),pos(:,2),'smoothingspline','SmoothingParam',1);
        valfitPoint = feval(myfitPoint,frq);
        
        %%%%%%%%%%% truncation only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        FidfitPoint = ifftshift(ifft(fftshift(valfitPoint)));
        FIDPointTruncOnly = complex(zeros(size(FidfitPoint,1),1),zeros(size(FidfitPoint,1),1));
        FIDPointTruncOnly((size(FidfitPoint,1)/2+1-HalfPointCut)+1:(size(FidfitPoint,1)/2+1+HalfPointCut-1),1) = FidfitPoint((size(FidfitPoint,1)/2+1-HalfPointCut)+1:(size(FidfitPoint,1)/2+1+HalfPointCut-1),1);
        FitPointTruncOnly = fftshift(fft(fftshift(FIDPointTruncOnly)));
        
        if result.TruncOnly ==1
            valfitPoint = real(FitPointTruncOnly);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        result.FitFinal = valfitPoint;
        % affichage
        delete(result.fi)
        hold on
        result.fi = plot(frq,valfitPoint,'Color',[1 0 0],'Parent',result.handles.AffichageSpectre);
        hold off
        uistack(result.fi,'bottom')
        result.pos = pos;
        setMyData(result);
        spectreFinal = result.spectre-valfitPoint;
        result.limX = get(result.handles.Affichage_res,'xlim');
        result.limY = get(result.handles.Affichage_res,'ylim');
        cla(result.handles.Affichage_res);
        plot(frq,spectreFinal,'Color',[0 0 0],'Parent',result.handles.Affichage_res);
        set(result.handles.Affichage_res,'xlim',result.limX);
        set(result.handles.Affichage_res,'ylim',result.limY);
        set(result.handles.Affichage_res,'XtickLabel',[],'YtickLabel',[]);
        if isfield(result,'agreement') ==1
            spectreFinaltemp = spectre-result.FitFinal;
            HalfPointCut = result.HalfPointCut;
            PointCut = result.PointCut;
            Fidfull = ifftshift(ifft(fftshift(spectreFinaltemp)));
            FIDtronc = Fidfull;
            Fidfullref = ifftshift(ifft(fftshift(spectre)));
            %             FIDtronc((size(Fidfull,1)/2-HalfPointCut)+2:(size(Fidfull,1)/2+(PointCut-HalfPointCut-1))+1,1) = complex(zeros(PointCut-1,1),zeros(PointCut-1,1));
            FIDtronc((size(Fidfull,1)/2+1-HalfPointCut)+1:(size(Fidfull,1)/2+1+HalfPointCut-1),1) = Fidfullref((size(Fidfull,1)/2+1-HalfPointCut)+1:(size(Fidfull,1)/2+1+HalfPointCut-1),1);
            spectreCorrcalc = fftshift(fft(fftshift(FIDtronc)));
            diff = spectre-real(spectreCorrcalc);
            result.agreement = ((sum(abs(spectre))-sum(abs(diff)))/sum(abs(spectre)))*100;
            set(result.handles.AgreementValue,'String',result.agreement);
        end
    end
end

setMyData(result);


function numbpoints_Callback(hObject, eventdata, handles)
% hObject    handle to numbpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if isfield(result, 'numbpoints')==0;
    result.numbpoints = 100;
    set(result.handles.numbpoints,'String',result.numbpoints);
end
result.numbpoints = str2double(get(result.handles.numbpoints,'String'));
setMyData(result)
% Hints: get(hObject,'String') returns contents of numbpoints as text
%        str2double(get(hObject,'String')) returns contents of numbpoints as a double


% --- Executes during object creation, after setting all properties.
function numbpoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numbpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OpenASCII.
function OpenASCII_Callback(hObject, eventdata, handles)
% hObject    handle to OpenASCII (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get data struct stored in figure

clear result;

handles.output = hObject;
result = struct;
result.handles=handles;
% Update handles structure
guidata(hObject, handles);
% Store all data in the Gui
setMyData(result);

result=getMyData(); % Get data struct stored in figure
%% Ouverture du fichier
[fname,dirname] = uigetfile('*.txt', 'Open ASCII file');
result.dfile = [dirname fname];
result.fname = fname;
result.dirname = dirname;

result.dataASCII = load(result.dfile,'-ascii');
result.spectre = result.dataASCII(:,2);
frq = (1:size(result.spectre,1))';
meanfrq = mean(frq);
result.frq = frq-meanfrq;

result.TruncOnly =1;
set(result.handles.Trunc,'Value',1);
fid = real(ifft(result.spectre));
[pks,locs,w,p] =findpeaks(real(fid));
pmax = max(p);
pthresh = pmax/50;
ind = locs(find(p>pthresh));
result.NumTruncPointValue = 2*(ind(1)-1);
set(result.handles.NumTruncPoint,'String',result.NumTruncPointValue);



%% Affichage
result.spectreF1 = result.spectre;
result.spectreFinal = result.spectreF1;
result.FitFinal = zeros(size(result.spectreFinal,1),1);
result.valfitSinus= zeros(size(result.spectreFinal,1),1);

cla(result.handles.AffichageSpectre);
hold on
plot(frq,result.spectreFinal,'Color',[0 0 1],'Parent',result.handles.AffichageSpectre);
hold off
set(result.handles.AffichageSpectre,'xlim',[min(frq) max(frq)]);
set(result.handles.AffichageSpectre,'ylim',[min(result.spectreFinal) max(result.spectreFinal)/5]);
set(result.handles.AffichageSpectre,'XtickLabel',[],'YtickLabel',[]);

cla(result.handles.Affichage_res);
set(result.handles.Affichage_res,'visible','off');


setMyData(result);
frequSinus_Callback
ModFit2_Callback
numbpoints_Callback
PowerStop_Callback
AdvanceFilter_Callback
NumTruncPoint_Callback
Trunc_Callback



% --- Executes on button press in SaveASCII.
function SaveASCII_Callback(hObject, eventdata, handles)
% hObject    handle to SaveASCII (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
axes(result.handles.AffichageSpectre);

%% Calcul de la partie imaginaire par Hilbert transform
if isfield(result, 'spectreImag')==1
    n = length(result.spectreImag);
    f = fft(result.spectreImag);
    h = [0 1i*ones(1,n/2-1) 0 -1i*ones(1,n/2-1)]';
    hf = h.*f;
    spectreImagShift = real(ifft(hf));
    spectreImagShiftcorr = spectreImagShift - result.FitFinal;
    f2 = fft(spectreImagShiftcorr);
    h = [0 1i*ones(1,n/2-1) 0 -1i*ones(1,n/2-1)]';
    hf2 = h.*f2;
    result.spectreFinalImag = -real(ifft(hf2));
else
    n = length(result.spectreFinal);
    f = fft(result.spectreFinal);
    h = [0 1i*ones(1,n/2-1) 0 -1i*ones(1,n/2-1)]';
    hf = h.*f;
    result.spectreFinalImag = -real(ifft(hf));
end

%% ecriture du fichier
% mkdir(result.dfile)
[dirname,name,ext] = fileparts(result.dfile);
NewASCIIfile = [dirname filesep name  '_corr.txt'];

%VSK export data
dataASCII = [result.dataASCII(:,1) , result.spectreFinal];
fid = fopen(NewASCIIfile,'w');
fprintf(fid,'%d\t%d\n',dataASCII');
fclose(fid);
%VSK export data end
setMyData(result)



% --- Executes on button press in Reload.
function Reload_Callback(hObject, eventdata, handles)
% hObject    handle to Reload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get data struct stored in figure
axes(result.handles.AffichageSpectre);

clear result.spectre result.spectreImag result.spectreF1 result.spectreFinal;

%% Reload
fid=fopen(result.dfile,'r','l');
[data,~]=fread(fid,'int32');
fclose(fid);
result.spectre = data;

result.dfileImag = [result.dfile(1:end-2) '1i'];
fid=fopen(result.dfileImag,'r','l');
[dataImag,~]=fread(fid,'int32');
fclose(fid);
result.spectreImag = dataImag;



%% Affichage
result.spectreF1 = result.spectre;
result.spectreFinal = result.spectreF1;
result.FitFinal = zeros(size(result.spectreFinal,1),1);
result.valfitSinus= zeros(size(result.spectreFinal,1),1);

cla(result.handles.AffichageSpectre);
hold on
plot(result.frq,result.spectreFinal,'Color',[0 0 1],'Parent',result.handles.AffichageSpectre);
plot(result.frq,result.FitFinal,'Color',[1 0 0],'Parent',result.handles.AffichageSpectre);
hold off
set(result.handles.AffichageSpectre,'xlim',[min(result.frq) max(result.frq)]);
set(result.handles.AffichageSpectre,'ylim',[min(result.spectreFinal) max(result.spectreFinal)/5]);
set(result.handles.AffichageSpectre,'XtickLabel',[],'YtickLabel',[]);

cla(result.handles.Affichage_res);
set(result.handles.Affichage_res,'visible','off');
set(result.handles.TemporalAxe,'visible','on');
axes(result.handles.AffichageSpectre);
setMyData(result);
NumTruncPoint_Callback()


function DE_Callback(hObject, eventdata, handles)
% hObject    handle to DE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DE as text
%        str2double(get(hObject,'String')) returns contents of DE as a double
result=getMyData(); % Get data struct stored in figure
if isfield(result, 'DEvalue')==0;
    result.DEvalue = 10;
    set(result.handles.DE, 'ForegroundColor',[1 0 0]);
    set(result.handles.DE,'String',result.DEvalue);
else
    set(result.handles.DE, 'ForegroundColor',[0 0 0]);
    result.frequSinus = 3.1416*result.DEvalue/(result.DWvalue*(size(result.spectre,1)));
    set(result.handles.frequSinus,'String',result.frequSinus);
end
result.DEvalue = str2double(get(result.handles.DE,'String'));

setMyData(result);


% --- Executes during object creation, after setting all properties.
function DE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DW_Callback(hObject, eventdata, handles)
% hObject    handle to DW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DW as text
%        str2double(get(hObject,'String')) returns contents of DW as a double
result=getMyData(); % Get data struct stored in figure
if isfield(result, 'DWvalue')==0;
    result.DWvalue = 0.25;
    set(result.handles.DW, 'ForegroundColor',[1 0 0]);
    set(result.handles.DW,'String',result.DWvalue);
else
    set(result.handles.DW, 'ForegroundColor',[0 0 0]);
    result.frequSinus = 3.1416*result.DEvalue/(result.DWvalue*(size(result.spectre,1)));
    set(result.handles.frequSinus,'String',result.frequSinus);
end
result.DWvalue = str2double(get(result.handles.DW,'String'));

setMyData(result);


% --- Executes during object creation, after setting all properties.
function DW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PowerStop_Callback(hObject, eventdata, handles)
% hObject    handle to PowerStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if isfield(result, 'PowerStopValue')==0;
    result.PowerStopValue = 6;
    set(result.handles.PowerStop,'String',result.PowerStopValue);
end
result.PowerStopValue = str2double(get(result.handles.PowerStop,'String'));
setMyData(result)
% Hints: get(hObject,'String') returns contents of PowerStop as text
%        str2double(get(hObject,'String')) returns contents of PowerStop as a double


% --- Executes during object creation, after setting all properties.
function PowerStop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PowerStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Trunc.
function Trunc_Callback(hObject, eventdata, handles)
% hObject    handle to Trunc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if isfield(result, 'TruncOnly')==0;
    result.TruncOnly = 0;
    set(result.handles.Trunc,'Value',result.TruncOnly);
end
result.TruncOnly = get(result.handles.Trunc,'Value');
setMyData(result)
% Hint: get(hObject,'Value') returns toggle state of Trunc



function NumTruncPoint_Callback(hObject, eventdata, handles)
% hObject    handle to NumTruncPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
result.NumTruncPointValue = str2double(get(result.handles.NumTruncPoint,'String'));
if result.NumTruncPointValue/2 ~= round(result.NumTruncPointValue/2);
    result.NumTruncPointValue = (ceil(result.NumTruncPointValue/2))*2;
end
set(result.handles.NumTruncPoint,'String',result.NumTruncPointValue);

PointCut = result.NumTruncPointValue;
HalfPointCut = PointCut/2;

FIDSpectre = ifftshift(ifft(fftshift(result.spectre)));
FIDSpectreSum = abs(real(FIDSpectre))+abs(imag(FIDSpectre));
ZoneInf = FIDSpectreSum((size(FIDSpectre,1)/2+1-HalfPointCut)-10:(size(FIDSpectre,1)/2+1),1);
ZoneSup = FIDSpectreSum((size(FIDSpectre,1)/2+1):(size(FIDSpectre,1)/2+1+(HalfPointCut)+10),1);
FIDtrunc = (complex(ones(size(FIDSpectre,1),1),ones(size(FIDSpectre,1),1)));
FIDtrunc((size(FIDSpectre,1)/2+1-HalfPointCut+1):(size(FIDSpectre,1)/2+1+HalfPointCut-1),1) = complex(zeros(PointCut-1,1),zeros(PointCut-1,1));
FIDtrunc = FIDtrunc.*max(real([FIDSpectreSum((size(FIDSpectre,1)/2+1-HalfPointCut),1) FIDSpectreSum((size(FIDSpectre,1)/2+1+HalfPointCut),1)]));

cla(result.handles.Affichage_res);
set(result.handles.Affichage_res,'visible','off');
cla(result.handles.TemporalAxe);
set(result.handles.TemporalAxe,'visible','on');
axes(result.handles.TemporalAxe);
hold on
plot(real(FIDSpectreSum),'b','Parent',result.handles.TemporalAxe);
plot(real(FIDtrunc),'r','Parent',result.handles.TemporalAxe);
set(result.handles.TemporalAxe,'xlim',[(size(FIDSpectre,1)/2+1-HalfPointCut)-10 (size(FIDSpectre,1)/2+1+HalfPointCut+10)]);
set(result.handles.TemporalAxe,'XtickLabel',[],'YtickLabel',[]);
hold off
axes(result.handles.AffichageSpectre);
result.HalfPointCut = HalfPointCut;
result.PointCut = PointCut;

setMyData(result)
% Hints: get(hObject,'String') returns contents of NumTruncPoint as text
%        str2double(get(hObject,'String')) returns contents of NumTruncPoint as a double


% --- Executes during object creation, after setting all properties.
function NumTruncPoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumTruncPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- log file.
function logfile()
result=getMyData();
[dirname,name,ext] = fileparts(result.dfile);
% option 1r bruker format -> put log file in title
if strcmp(name,'1r')
    [dirname,name,ext] = fileparts(result.dfile);
    %in title
    fname = [dirname filesep 'title'];
    [fid, errmsg] =fopen(fname,'a');
    line0 = '-----------------------';
    line1 = ['Baseline-corrected after ' num2str(result.iter) ' iterations'];
    line2 = ['Convergence criterion value : ', num2str(result.critere)];
    if result.TruncOnly == 1
        line3 = 'Truncation only option : yes';
        line4 = ['Truncation filter : ' num2str(result.PointCut) ' points'];
    else
        line3 = 'Truncation only option : no';
        line4 = '';
    end
    line5 = ['Threshold level : ' num2str(result.AdvanceFilterValue) ' %'];
    fprintf(fid,'\n%s\n%s\n%s\n%s\n%s\n%s\n%s',line0,line1,line2,line3,line4,line5,date);
    fclose(fid);
    
else
    fname = [dirname filesep name '.log'];
    [fid, errmsg] =fopen(fname,'w');
    line0 = '-----------------------';
    line1 = ['Baseline-corrected after ' num2str(result.iter) ' iterations'];
    line2 = ['Convergence criterion value : ', num2str(result.critere)];
    if result.TruncOnly == 1
        line3 = 'Truncation only option : yes';
        line4 = ['Truncation filter : ' num2str(result.PointCut) ' points'];
    else
        line3 = 'Truncation only option : no';
        line4 = '';
    end
    line5 = ['Threshold level : ' num2str(result.AdvanceFilterValue) ' %'];
    fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s',line1,line2,line3,line4,line5,date);
    fclose(fid);
end



%% Expert mode options

% --- Executes on button press in AdvanceMode.
function AdvanceMode_Callback(hObject, eventdata, handles)
% hObject    handle to AdvanceMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if  get(result.handles.AdvanceMode,'value')==1
    set(result.handles.AdvanceFilter,'visible','on');
    set(result.handles.textAdvanceValue,'visible','on');
    set(result.handles.PeakFinder,'visible','on');
    set(result.handles.SmartW,'visible','on');
    set(result.handles.SymThresh,'visible','on');
    set(result.handles.checkbox6,'visible','on');
else
    set(result.handles.AdvanceFilter,'visible','off');
    set(result.handles.textAdvanceValue,'visible','off');
    set(result.handles.PeakFinder,'visible','off');
    set(result.handles.SmartW,'visible','off');
    set(result.handles.SymThresh,'visible','off');
    set(result.handles.checkbox6,'visible','off');
end
% Store all data in the Gui
setMyData(result);

function AdvanceFilter_Callback(hObject, eventdata, handles)
% hObject    handle to AdvanceFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if isfield(result, 'AdvanceFilterValue')==0
    result.AdvanceFilterValue = 100;
    set(result.handles.AdvanceFilter,'String',result.AdvanceFilterValue);
end
result.AdvanceFilterValue = str2double(get(result.handles.AdvanceFilter,'String'));
setMyData(result)
% Hints: get(hObject,'String') returns contents of AdvanceFilter as text
%        str2double(get(hObject,'String')) returns contents of AdvanceFilter as a double

% --- Executes during object creation, after setting all properties.
function AdvanceFilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AdvanceFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in PeakFinder.
function PeakFinder_Callback(hObject, eventdata, handles)
% hObject    handle to PeakFinder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if isfield(result, 'PeakFinder')==0
    result.PeakFinder = 0;
end
setMyData(result)
% Hint: get(hObject,'Value') returns toggle state of PeakFinder


% --- Executes on button press in SmartW.
function SmartW_Callback(hObject, eventdata, handles)
% hObject    handle to SmartW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if isfield(result, 'SW')==0
    result.SW = 1;
end
setMyData(result)
% Hint: get(hObject,'Value') returns toggle state of SmartW



function iterations_Callback(hObject, eventdata, handles)
% hObject    handle to iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iterations as text
%        str2double(get(hObject,'String')) returns contents of iterations as a double
result=getMyData(); % Get back the result structure
if isfield(result, 'iter')==0
    result.iter = 0;
end
setMyData(result)


% --- Executes on button press in SymThresh.
function SymThresh_Callback(hObject, eventdata, handles)
% hObject    handle to SymThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if isfield(result, 'SymThresh')==0
    result.SymThresh = 0;
    set(result.handles.SymThresh,'Value',result.SymThresh);
end
setMyData(result)
% Hint: get(hObject,'Value') returns toggle state of SymThresh



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if isfield(result, 't')==1;
    delete(result.t)
end
setMyData(result)
% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getMyData(); % Get back the result structure
if isfield(result, 'ShowHist')==0
    result.ShowHist = 0;
end
setMyData(result)
% Hint: get(hObject,'Value') returns toggle state of checkbox6
