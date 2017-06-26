function varargout = minusConver(varargin)
% MINUSCONVER MATLAB code for minusConver.fig
%      MINUSCONVER, by itself, creates a new MINUSCONVER or raises the existing
%      singleton*.
%
%      H = MINUSCONVER returns the handle to a new MINUSCONVER or the handle to
%      the existing singleton*.
%
%      MINUSCONVER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MINUSCONVER.M with the given input arguments.
%
%      MINUSCONVER('Property','Value',...) creates a new MINUSCONVER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before minusConver_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to minusConver_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help minusConver

% Last Modified by GUIDE v2.5 19-Oct-2016 09:41:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @minusConver_OpeningFcn, ...
                   'gui_OutputFcn',  @minusConver_OutputFcn, ...
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


% --- Executes just before minusConver is made visible.
function minusConver_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to minusConver (see VARARGIN)

% Choose default command line output for minusConver
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes minusConver wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = minusConver_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filename;
global pathname;
global filedata;
if ~ischar(pathname)
    pathname=cd();
end
[filename,pathname]=uigetfile(strcat(pathname,'\*.txt'));
open_pathname=pathname;
if isequal(filename,0) || isequal(open_pathname,0)
else
    fidin = fopen(strcat(pathname,filename));
    while ~feof(fidin)
        tline=fgetl(fidin);
        if double(tline(1)) >= 48 && double(tline(1)) <= 57
            filedata=[filedata;str2num(tline)];
            continue;
        end
    end
    fclose(fidin);
    x = filedata(:,1);
    y = filedata(:,2);
    
    axes(handles.axes);
    cla reset;
    plot(x,y,'r.');
end
    
    


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filedata;
global pathname;
global filename;

if ~isequal(filedata,0)
    convData = filedata;
    convData(:,2) = convData(:,2)*(-1);
    filename = strcat('Conv_',filename);
    save(strcat(pathname,filename),'convData','-ascii');
    x = convData(:,1);
    y = convData(:,2);
    
    axes(handles.axes);
    cla reset;
    plot(x,y,'r.');
    filedata=[];
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir;
path=strcat(path,'\');
filelist = dir(fullfile(path,'*.txt'));
n=length(filelist);
convData=[];

if n==0
    msgbox('can not find the file!');
else
    resultPath = uigetdir;
    if ~isequal(resultPath,0)
        h=waitbar(0,'please wait....');
        for i=1:n
            fidin = fopen(strcat(path,filelist(i).name));
            while ~feof(fidin)
                tline=fgetl(fidin);
                if (double(tline(1)) >= 48 && double(tline(1)) <= 57) || '-'==tline(1)
                    convData=[convData;str2num(tline)];
                    continue;
                end
            end
            fclose(fidin);
            convData(:,2) = convData(:,2)*(-1);
            resultFilename = strcat(resultPath,'\');
            resultFilename = strcat(resultFilename,filelist(i).name);
             fid = fopen(resultFilename,'w+');
             for j=1:length(convData(:,1))
                 fprintf(fid,'%f\t%E\n',convData(j,1),convData(j,2));
             end
 %           save(resultFilename,'convData','-ascii');
             fclose(fid);
             convData=[];
            waitbar(i/n,h);
        end
        close(h);
    end
end
    

            
