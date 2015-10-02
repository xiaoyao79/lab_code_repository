function varargout = ThreePI(varargin)
% THREEPI M-file for ThreePI.fig
%      THREEPI, by itself, creates a new THREEPI or raises the existing
%      singleton*.
%
%      H = THREEPI returns the handle to a new THREEPI or the handle to
%      the existing singleton*.
%
%      THREEPI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THREEPI.M with the given input arguments.
%
%      THREEPI('Property','Value',...) creates a new THREEPI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ThreePI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ThreePI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ThreePI

% Last Modified by GUIDE v2.5 07-May-2009 09:11:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ThreePI_OpeningFcn, ...
                   'gui_OutputFcn',  @ThreePI_OutputFcn, ...
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


% --- Executes just before ThreePI is made visible.
function ThreePI_OpeningFcn(hObject, eventdata, handles, varargin)%#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ThreePI (see VARARGIN)

% Choose default command line output for ThreePI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ThreePI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ThreePI_OutputFcn(hObject, eventdata, handles)%#ok 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Load_Button.
function Load_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[uijobfile, uijobpath] =uigetfile('Job File');
if uijobfile==0
    return
end
load([uijobpath uijobfile]);

set(handles.ImDir,'string',Data.plt_info.pltdir);
set(handles.BaseName,'string',Data.plt_info.pltname);
set(handles.NumZeros,'string',Data.plt_info.numzeros);
set(handles.FrmSkip,'string',Data.plt_info.skipframes);
set(handles.StrFrm,'string',Data.plt_info.startframe);
set(handles.EndFrm,'string',Data.plt_info.endframe);

set(handles.Unit_drop,'Value',Data.phys.unitvar);
set(handles.NonDim_Check,'value',Data.phys.nondim);
set(handles.Char_Vel,'string',Data.phys.charvel);
set(handles.Char_Length,'string',Data.phys.charlen);
set(handles.Freq_in,'string',Data.phys.freq);
set(handles.Pix_size,'string',Data.phys.micrperpix);
set(handles.Pulse_fillin,'string',Data.phys.pulsesep);
set(handles.Visc_val,'string',Data.phys.kinvisc);

set(handles.Steady_check,'value',Data.calc.steady);
set(handles.Unsteady_check,'value',Data.calc.unsteady);
set(handles.Vel_Mag_Check,'value',Data.calc.velmagvar);
set(handles.Vort_Check,'value',Data.calc.vorticityvar);
set(handles.Strain_Check,'value',Data.calc.strainratevar);
set(handles.TurbKE_Check,'value',Data.calc.turbKEvar);
set(handles.TotKE_Check,'value',Data.calc.totalKEvar);
set(handles.ReStress_Check,'value',Data.calc.reynoldsstressvar);
set(handles.TurbStat_Check,'value',Data.calc.turbstats);
set(handles.Diss_1_Check,'value',Data.calc.diss1term);
set(handles.Diss_1_drop,'value',Data.calc.diss1var);
set(handles.Diss_5_Check,'value',Data.calc.diss5term);
set(handles.Diss_LES_Check,'value',Data.calc.dissLES);
set(handles.Diss_LES_drop,'value',Data.calc.dissLESmeth);
set(handles.Diss_LES_Coef,'string',Data.calc.dissLEScoef);
set(handles.Diss_LES_Filt,'string',Data.calc.dissLESfilt);
set(handles.Diss_Comp_Check,'value',Data.calc.dissCOMP);
set(handles.Grad_Stat,'value',Data.calc.gradientsvar);
set(handles.Grad_type_Check,'value',Data.calc.gradienttype);
set(handles.Grad_Drop,'value',Data.calc.gradientmethod);
set(handles.VortID,'value',Data.calc.vortexID);
set(handles.VortID_Method,'value',Data.calc.VortexIDmethod);

set(handles.Job_name,'string',Data.output.job);
set(handles.Out_Name,'string',Data.output.savename);
set(handles.Out_Type_Drop,'value',Data.output.fileformat);
set(handles.SaveLoc_fillin,'string',Data.output.saveloc);

fprintf('Loaded\n')



% --- Executes on button press in Save_Button.
function Save_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Data.plt_info.pltdir          = get(handles.ImDir,'string');
Data.plt_info.pltname         = get(handles.BaseName,'string');
Data.plt_info.numzeros        = get(handles.NumZeros,'string');
Data.plt_info.skipframes      = get(handles.FrmSkip,'string');
Data.plt_info.startframe      = get(handles.StrFrm,'string');
Data.plt_info.endframe        = get(handles.EndFrm,'string');

Data.phys.unitvar             = get(handles.Unit_drop,'Value');
Data.phys.nondim              = get(handles.NonDim_Check,'value');%%%%%%%
Data.phys.charvel             = get(handles.Char_Vel,'string');
Data.phys.charlen             = get(handles.Char_Length,'string');
Data.phys.freq                = get(handles.Freq_in,'string');
Data.phys.micrperpix          = get(handles.Pix_size,'string');
Data.phys.pulsesep            = get(handles.Pulse_fillin,'string');
Data.phys.kinvisc             = get(handles.Visc_val,'string');

Data.phys.xoffset             = 0;
Data.phys.yoffset             = 0;
Data.phys.zoffset             = 0;

Data.calc.steady              = get(handles.Steady_check,'value');
Data.calc.unsteady            = get(handles.Unsteady_check,'value');
Data.calc.velmagvar           = get(handles.Vel_Mag_Check,'value');
Data.calc.vorticityvar        = get(handles.Vort_Check,'value');
Data.calc.strainratevar       = get(handles.Strain_Check,'value');
Data.calc.turbKEvar           = get(handles.TurbKE_Check,'value');
Data.calc.totalKEvar          = get(handles.TotKE_Check,'value');
Data.calc.reynoldsstressvar   = get(handles.ReStress_Check,'value');
Data.calc.turbstats           = get(handles.TurbStat_Check,'value');
Data.calc.diss1term           = get(handles.Diss_1_Check,'value');
Data.calc.diss1var            = get(handles.Diss_1_drop,'value');%%%%%%%%%%%
Data.calc.diss5term           = get(handles.Diss_5_Check,'value');
Data.calc.dissLES             = get(handles.Diss_LES_Check,'value');
Data.calc.dissLESmeth         = get(handles.Diss_LES_drop,'value');%%%%%%%%%
Data.calc.dissLEScoef         = get(handles.Diss_LES_Coef,'string');%%%%%%%%%
Data.calc.dissLESfilt         = get(handles.Diss_LES_Filt,'string');%%%%%%%%%
Data.calc.dissCOMP            = get(handles.Diss_Comp_Check,'value');
Data.calc.gradientsvar        = get(handles.Grad_Stat,'value');
Data.calc.gradienttype        = get(handles.Grad_type_Check,'value');
Data.calc.gradientmethod      = get(handles.Grad_Drop,'value');
Data.calc.vortexID            = get(handles.VortID,'value');
Data.calc.VortexIDmethod      = get(handles.VortID_Method,'value');

Data.output.job               = get(handles.Job_name,'string');
Data.output.savename          = get(handles.Out_Name,'string');
Data.output.fileformat        = get(handles.Out_Type_Drop,'value');
Data.output.saveloc           = get(handles.SaveLoc_fillin,'string');

if ispc
    if ~strcmpi(Data.plt_info.pltdir(end),'\')
        Data.plt_info.pltdir = sprintf('%s%s',Data.plt_info.pltdir,'\');
    end
    if ~strcmpi(Data.output.saveloc(end),'\')
        Data.output.saveloc = sprintf('%s%s',Data.output.saveloc,'\');
    end
else
    if ~strcmpi(Data.plt_info.pltdir(end),'/')
        Data.plt_info.pltdir = sprintf('%s%s',Data.plt_info.pltdir,'/');
    end
    if ~strcmpi(Data.output.saveloc(end),'/')
        Data.output.saveloc = sprintf('%s%s',Data.output.saveloc,'/');
    end
end

[saveName,saveLoc] = uiputfile(sprintf('%s',Data.output.job));
save(sprintf('%s%s',saveLoc,saveName),'Data')
fprintf('%s%s\n',saveLoc,saveName)
fprintf('Saved\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Button
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Run_Button.
function Run_Button_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Run_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Data.plt_info.pltdir          = get(handles.ImDir,'string');
Data.plt_info.pltname         = get(handles.BaseName,'string');
Data.plt_info.numzeros        = get(handles.NumZeros,'string');
Data.plt_info.skipframes      = get(handles.FrmSkip,'string');
Data.plt_info.startframe      = get(handles.StrFrm,'string');
Data.plt_info.endframe        = get(handles.EndFrm,'string');

Data.phys.unitvar             = get(handles.Unit_drop,'Value');
Data.phys.nondim              = get(handles.NonDim_Check,'value');%%%%%%%
Data.phys.charvel             = get(handles.Char_Vel,'string');
Data.phys.charlen             = get(handles.Char_Length,'string');
Data.phys.freq                = get(handles.Freq_in,'string');
Data.phys.micrperpix          = get(handles.Pix_size,'string');
Data.phys.pulsesep            = get(handles.Pulse_fillin,'string');
Data.phys.kinvisc             = get(handles.Visc_val,'string');

Data.phys.xoffset             = 0;
Data.phys.yoffset             = 0;
Data.phys.zoffset             = 0;

Data.calc.steady              = get(handles.Steady_check,'value');
Data.calc.unsteady            = get(handles.Unsteady_check,'value');
Data.calc.velmagvar           = get(handles.Vel_Mag_Check,'value');
Data.calc.vorticityvar        = get(handles.Vort_Check,'value');
Data.calc.strainratevar       = get(handles.Strain_Check,'value');
Data.calc.turbKEvar           = get(handles.TurbKE_Check,'value');
Data.calc.totalKEvar          = get(handles.TotKE_Check,'value');
Data.calc.reynoldsstressvar   = get(handles.ReStress_Check,'value');
Data.calc.turbstats           = get(handles.TurbStat_Check,'value');
Data.calc.diss1term           = get(handles.Diss_1_Check,'value');
Data.calc.diss1var            = get(handles.Diss_1_drop,'value');%%%%%%%%%%%
Data.calc.diss5term           = get(handles.Diss_5_Check,'value');
Data.calc.dissLES             = get(handles.Diss_LES_Check,'value');
Data.calc.dissLESmeth         = get(handles.Diss_LES_drop,'value');%%%%%%%%%
Data.calc.dissLEScoef         = get(handles.Diss_LES_Coef,'string');%%%%%%%%%
Data.calc.dissLESfilt         = get(handles.Diss_LES_Filt,'string');%%%%%%%%%
Data.calc.dissCOMP            = get(handles.Diss_Comp_Check,'value');
Data.calc.gradientsvar        = get(handles.Grad_Stat,'value');
Data.calc.gradienttype        = get(handles.Grad_type_Check,'value');
Data.calc.gradientmethod      = get(handles.Grad_Drop,'value');
Data.calc.vortexID            = get(handles.VortID,'value');
Data.calc.VortexIDmethod      = get(handles.VortID_Method,'value');

Data.output.job               = get(handles.Job_name,'string');
Data.output.savename          = get(handles.Out_Name,'string');
Data.output.fileformat        = get(handles.Out_Type_Drop,'value');
Data.output.saveloc           = get(handles.SaveLoc_fillin,'string');

if ispc
    if ~strcmpi(Data.plt_info.pltdir(end),'\')
        Data.plt_info.pltdir = sprintf('%s%s',Data.plt_info.pltdir,'\');
    end
    if ~strcmpi(Data.output.saveloc(end),'\')
        Data.output.saveloc = sprintf('%s%s',Data.output.saveloc,'\');
    end
else
    if ~strcmpi(Data.plt_info.pltdir(end),'/')
        Data.plt_info.pltdir = sprintf('%s%s',Data.plt_info.pltdir,'/');
    end
    if ~strcmpi(Data.output.saveloc(end),'/')
        Data.output.saveloc = sprintf('%s%s',Data.output.saveloc,'/');
    end
end

ThreePIcode(Data);



























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Support Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ImDir_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to ImDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImDir as text
%        str2double(get(hObject,'String')) returns contents of ImDir as a double


% --- Executes during object creation, after setting all properties.
function ImDir_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to ImDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Browse_PLT.
function Browse_PLT_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Browse_PLT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path=uigetdir('PLT Directory');
set(handles.ImDir,'String',path);
%If the user hits "cancel":
if path==0
    set(handles.ImDir,'String',' ');
    return
end



function NumZeros_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to NumZeros (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumZeros as text
%        str2double(get(hObject,'String')) returns contents of NumZeros as a double


% --- Executes during object creation, after setting all properties.
function NumZeros_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to NumZeros (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FrmSkip_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to FrmSkip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrmSkip as text
%        str2double(get(hObject,'String')) returns contents of FrmSkip as a double


% --- Executes during object creation, after setting all properties.
function FrmSkip_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to FrmSkip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StrFrm_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to StrFrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StrFrm as text
%        str2double(get(hObject,'String')) returns contents of StrFrm as a double


% --- Executes during object creation, after setting all properties.
function StrFrm_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to StrFrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EndFrm_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to EndFrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EndFrm as text
%        str2double(get(hObject,'String')) returns contents of EndFrm as a double


% --- Executes during object creation, after setting all properties.
function EndFrm_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to EndFrm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BaseName_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to BaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BaseName as text
%        str2double(get(hObject,'String')) returns contents of BaseName as a double


% --- Executes during object creation, after setting all properties.
function BaseName_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to BaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Char_Length_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Char_Length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Char_Length as text
%        str2double(get(hObject,'String')) returns contents of Char_Length as a double


% --- Executes during object creation, after setting all properties.
function Char_Length_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Char_Length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Unit_drop.
function Unit_drop_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Unit_drop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Unit_drop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Unit_drop


% --- Executes during object creation, after setting all properties.
function Unit_drop_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Unit_drop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Char_Vel_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Char_Vel_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Char_Vel_text as text
%        str2double(get(hObject,'String')) returns contents of Char_Vel_text as a double


% --- Executes during object creation, after setting all properties.
function Char_Vel_text_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Char_Vel_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Out_Name_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Out_Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Out_Name as text
%        str2double(get(hObject,'String')) returns contents of Out_Name as a double


% --- Executes during object creation, after setting all properties.
function Out_Name_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Out_Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Steady_check.
function Steady_check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Steady_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Steady_check


% --- Executes on button press in Unsteady_check.
function Unsteady_check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Unsteady_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Unsteady_check


% --- Executes on selection change in Out_Type_Drop.
function Out_Type_Drop_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Out_Type_Drop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Out_Type_Drop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Out_Type_Drop


% --- Executes during object creation, after setting all properties.
function Out_Type_Drop_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Out_Type_Drop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Job_name_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Job_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Job_name as text
%        str2double(get(hObject,'String')) returns contents of Job_name as a double


% --- Executes during object creation, after setting all properties.
function Job_name_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Job_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Vel_Mag_Check.
function Vel_Mag_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Vel_Mag_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Vel_Mag_Check


% --- Executes on button press in Vort_Check.
function Vort_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Vort_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Vort_Check


% --- Executes on button press in Strain_Check.
function Strain_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Strain_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Strain_Check


% --- Executes on button press in TurbKE_Check.
function TurbKE_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to TurbKE_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TurbKE_Check


% --- Executes on button press in TotKE_Check.
function TotKE_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to TotKE_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TotKE_Check


% --- Executes on button press in ReStress_Check.
function ReStress_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to ReStress_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ReStress_Check


% --- Executes on button press in TurbStat_Check.
function TurbStat_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to TurbStat_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TurbStat_Check


% --- Executes on button press in Diss_1_Check.
function Diss_1_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_1_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Diss_1_Check


% --- Executes on button press in Diss_5_Check.
function Diss_5_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_5_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Diss_5_Check


% --- Executes on button press in Diss_LES_Check.
function Diss_LES_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_LES_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Diss_LES_Check


% --- Executes on selection change in Grad_Drop.
function Grad_Drop_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Grad_Drop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Grad_Drop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Grad_Drop


% --- Executes during object creation, after setting all properties.
function Grad_Drop_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Grad_Drop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Freq_in_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Freq_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Freq_in as text
%        str2double(get(hObject,'String')) returns contents of Freq_in as a double


% --- Executes during object creation, after setting all properties.
function Freq_in_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Freq_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pulse_fillin_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Pulse_fillin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pulse_fillin as text
%        str2double(get(hObject,'String')) returns contents of Pulse_fillin as a double


% --- Executes during object creation, after setting all properties.
function Pulse_fillin_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Pulse_fillin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pix_size_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Pix_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pix_size as text
%        str2double(get(hObject,'String')) returns contents of Pix_size as a double


% --- Executes during object creation, after setting all properties.
function Pix_size_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Pix_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Visc_val_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Visc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Visc_val as text
%        str2double(get(hObject,'String')) returns contents of Visc_val as a double


% --- Executes during object creation, after setting all properties.
function Visc_val_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Visc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Grad_Stat.
function Grad_Stat_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Grad_Stat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Grad_Stat


% --- Executes during object creation, after setting all properties.
function Char_Vel_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Char_Vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Diss_Comp_Check.
function Diss_Comp_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_Comp_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Diss_Comp_Check


% --- Executes on button press in Grad_type_Check.
function Grad_type_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Grad_type_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Grad_type_Check


% --- Executes on selection change in Diss_1_drop.
function Diss_1_drop_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_1_drop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Diss_1_drop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Diss_1_drop


% --- Executes during object creation, after setting all properties.
function Diss_1_drop_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_1_drop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Diss_LES_drop.
function Diss_LES_drop_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_LES_drop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Diss_LES_drop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Diss_LES_drop


% --- Executes during object creation, after setting all properties.
function Diss_LES_drop_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_LES_drop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Diss_LES_Coef_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_LES_Coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Diss_LES_Coef as text
%        str2double(get(hObject,'String')) returns contents of Diss_LES_Coef as a double


% --- Executes during object creation, after setting all properties.
function Diss_LES_Coef_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_LES_Coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Diss_LES_Filt_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_LES_Filt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Diss_LES_Filt as text
%        str2double(get(hObject,'String')) returns contents of Diss_LES_Filt as a double


% --- Executes during object creation, after setting all properties.
function Diss_LES_Filt_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to Diss_LES_Filt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NonDim_Check.
function NonDim_Check_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to NonDim_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NonDim_Check


% --- Executes on button press in VortID.
function VortID_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to VortID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of VortID


% --- Executes on selection change in VortID_Method.
function VortID_Method_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to VortID_Method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns VortID_Method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from VortID_Method


% --- Executes during object creation, after setting all properties.
function VortID_Method_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to VortID_Method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SaveLoc_fillin_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to SaveLoc_fillin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SaveLoc_fillin as text
%        str2double(get(hObject,'String')) returns contents of SaveLoc_fillin as a double


% --- Executes during object creation, after setting all properties.
function SaveLoc_fillin_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to SaveLoc_fillin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Saveloc_button.
function Saveloc_button_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to Saveloc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path=uigetdir('PLT Directory');
set(handles.SaveLoc_fillin,'String',path);
%If the user hits "cancel":
if path==0
    set(handles.SaveLoc_fillin,'String',' ');
    return
end


