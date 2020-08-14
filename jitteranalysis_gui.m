
% jitteranalysis_gui.m is a user interface to display and save jitter analysis between bursts 
%   [spike_jitter,jitter_ave,jitter_overall] = jitteranalysis_gui(spike_times)   
%           spike_jitter - returns the jitter between different bursts in the waveform; 
%                   spike_jitter(:,1,2) = the jitter between burst 1 and burst 2, the first dimension contains the jitter between each spike_interval
%                   Note: Since each burst comparison uses a different number of spike times to calculate (different bursts
%                         have different number of spikes recorded), many of the ending rows will store 0 (matrices must have same number of columns,rows,etc)
%                   Note: spike_jitter(:,2,1) = an array of zeros since the jitter is already stored between burst 1 and 2 in spike_jitter(:,1,2)
%                   
%           jitter_ave - returns the average jitter between two bursts, jitter_ave(1,2) = the average jitter between burst 1 and 2
%                   Note: jitter(2,1) = 0 since the jitter average is already stored between burst 1 and 2 in jitter(1,2)
%           jitter_overall - returns the overall jitter between the two bursts
%
%           spike_time - 2d matrix of time points that spikes occurred (each column represents spike times for a different burst)
%           align_time - vector of alignment spikes for each burst (used to align two bursts to each other before calculating jitter)


function varargout = jitteranalysis_gui(varargin)
% JITTERANALYSIS_GUI MATLAB code for jitteranalysis_gui.fig
%      JITTERANALYSIS_GUI, by itself, creates a new JITTERANALYSIS_GUI or raises the existing
%      singleton*.
%
%      H = JITTERANALYSIS_GUI returns the handle to a new JITTERANALYSIS_GUI or the handle to
%      the existing singleton*.
%
%      JITTERANALYSIS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JITTERANALYSIS_GUI.M with the given input arguments.
%
%      JITTERANALYSIS_GUI('Property','Value',...) creates a new JITTERANALYSIS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before jitteranalysis_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to jitteranalysis_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help jitteranalysis_gui

% Last Modified by GUIDE v2.5 22-Jul-2014 15:40:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @jitteranalysis_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @jitteranalysis_gui_OutputFcn, ...
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

% --- Executes just before jitteranalysis_gui is made visible.
function jitteranalysis_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to jitteranalysis_gui (see VARARGIN)

% Choose default command line output for jitteranalysis_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Spike times must be provided
spike_times = varargin{1};
align_spike = varargin{2};
        
%Initialize userdata in figure as empty
set(handles.maingui_fig,'UserData',[]);
%Store the spike times parameters
userData.spike_times = spike_times;
userData.align_spike = align_spike;

%Initialize spike jitter variables
userData.spike_jitter = [];
userData.jitter_ave = [];
userData.jitter_overall = [];
set(handles.maingui_fig,'UserData',userData);

%Display the number of bursts in the edit text
set(handles.numbursts_text,'String',size(align_spike,2))

% UIWAIT makes jitteranalysis_gui wait for user response (see UIRESUME)
uiwait(handles.maingui_fig);


% --- Outputs from this function are returned to the command line.
function varargout = jitteranalysis_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
userData = get(handles.maingui_fig,'UserData');
varargout{1} = userData.spike_jitter;
varargout{2} = userData.jitter_ave;
varargout{3} = userData.jitter_overall;
delete(handles.maingui_fig)


% --- Executes on button press in calculatejitter.
%Calculates the jitter between the provided spike times
function calculatejitter_Callback(hObject, eventdata, handles)
% hObject    handle to addspike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.maingui_fig,'UserData');
spike_times = userData.spike_times;
align_spike = userData.align_spike;

if isempty(spike_times) || isempty(align_spike)
    msgbox('Please select spike time data before calculating jitter.')
elseif size(spike_times,2) == 1 || size(align_spike,2) == 1
    msgbox('Please select spike time data that has more than 1 burst to perform jitter analysis.')
else
    %Pre-allocate space for spike_jitter and jitter_ave
    num_bursts = size(align_spike,2);
    spike_jitter = zeros(size(spike_times,1),num_bursts,num_bursts);
    jitter_ave = zeros(num_bursts,num_bursts);

    %Calculate the interval between spikes
    spike_interval(1:length(spike_times)-1,:) = spike_times(2:end,:)-spike_times(1:end-1,:);

    %Calculate the jitter for each of the bursts
    num_jitter=1;
    for burst1 = 1:num_bursts
        for burst2 = 1:num_bursts
            if burst1<burst2
                %Take out the non-zero spike_times from the matrix
                spike_time1 = spike_times(:,burst1);
                spike_time1 = spike_time1(spike_time1~=0);
                spike_time2 = spike_times(:,burst2);
                spike_time2 = spike_time2(spike_time2~=0);
                %Number of spikes before alignment spike
                [num_spikes_before, less_spikes] = min([length(spike_time1(spike_time1<align_spike(burst1))) length(spike_time2(spike_time2<align_spike(burst2)))]);
                %Number of spikes omitted from jitter analysis from the burst with more spikes before alignment spike
                spikes_omitted = max([length(spike_time1(spike_time1<align_spike(burst1))) length(spike_time2(spike_time2<align_spike(burst2)))]) - num_spikes_before;
                %Number of spikes after alignment spike
                num_spikes_after = min([length(spike_time1(spike_time1>align_spike(burst1))) length(spike_time2(spike_time2>align_spike(burst2)))]);
                num_spikes = num_spikes_before+num_spikes_after;

                if less_spikes == 1 %Seizure 1 has fewer spikes than seizure 2.
                    spike_jitter(1:num_spikes-1,burst1,burst2) = abs(spike_interval(1:num_spikes-1,burst1)-spike_interval(spikes_omitted+1:num_spikes+spikes_omitted-1,burst2))*200./...
                        (spike_interval(1:num_spikes-1,burst1)+spike_interval(spikes_omitted+1:num_spikes+spikes_omitted-1,burst2));
                else
                    spike_jitter(1:num_spikes-1,burst1,burst2) = abs(spike_interval(spikes_omitted+1:num_spikes+spikes_omitted-1,burst1)-spike_interval(1:num_spikes-1,burst2))*200./...
                        (spike_interval(spikes_omitted+1:num_spikes+spikes_omitted-1,burst1)+spike_interval(1:num_spikes-1,burst2));
                end
                jitter_ave(burst1,burst2) = mean(spike_jitter(1:num_spikes-1,burst1,burst2));
                %Save the data into table format
                spike_jitter_table(num_jitter,:) = [burst1,burst2,jitter_ave(burst1,burst2)];
                num_jitter = num_jitter+1;
            end
        end
    end
    %Add the spike jitter ave into the GUI table
    set(handles.spikejitter_table,'Data',spike_jitter_table)

    jitter_ave_temp = reshape(jitter_ave',size(jitter_ave,1)*size(jitter_ave,2),1);
    jitter_overall = mean(jitter_ave_temp(jitter_ave_temp~=0));
    %Display the overall  jitter on the GUI 
    set(handles.overalljitter_text,'String',jitter_overall)

    %Save jitter into userData
    userData.spike_jitter = spike_jitter;
    userData.jitter_ave = jitter_ave;
    userData.jitter_overall = jitter_overall;
    set(handles.maingui_fig,'UserData',userData);
end



% --- Executes on button press in saveanalysis.
% Saves the jitter analysis (jitter_ave and jitter_overall) to an excel spreadsheet
function saveanalysis_Callback(hObject, eventdata, handles)
% hObject    handle to removespike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.maingui_fig,'UserData');

filename = get(handles.filename,'String');
if isempty(filename)
    msgbox('Please enter a filename for the jitter analysis.')
else
    answer = 'Yes'; %Default to logical 1 to create/overwrite a file.
    if isfile([filename '.xlsx'])
        answer = questdlg(['A file with the name "' filename '" already exists. Would you like to overwrite it?'],'Confirm Save','Yes','No','No');
    end
    if strcmp(answer,'Yes')
        if isfile([filename '.xlsx'])
            delete([filename '.xlsx'])
        end
        jitter_ave = userData.jitter_ave;
        
        %Get the version number. If before R2019a, use xlswrite instead of writecell
        mat_version = version;
        mat_version_year = str2num(mat_version(find(mat_version=='R')+1:find(mat_version=='R')+4));
        
        if mat_version_year >= 2019
            writecell({'Burst Number','Burst Number','Jitter %'},[filename '.xlsx'])
            %Store the spike jitter in an excel spreadsheet
            jitter_count=0;
            for burst1 = 1:size(jitter_ave,1)
                for burst2 = 1:size(jitter_ave,2)
                    if burst1<burst2
                        row = num2str(jitter_count+2);
                        writematrix([burst1;burst2;jitter_ave(burst1,burst2)]',[filename '.xlsx'],'Range',['A' row ':C' row])
                        jitter_count = jitter_count+1;
                    end
                end
            end
            %Store the overall jitter in an excel spreadsheet
            row = num2str(jitter_count+3);
            writecell({'Overall Jitter',userData.jitter_overall},[filename '.xlsx'],'Range',['A' row ':B' row])
        else
            xlswrite([filename '.xlsx'],{'Burst Number','Burst Number','Jitter %'})
            %Store the spike jitter in an excel spreadsheet
            jitter_count=0;
            for burst1 = 1:size(jitter_ave,1)
                for burst2 = 1:size(jitter_ave,2)
                    if burst1<burst2
                        row = num2str(jitter_count+2);
                        xlswrite([filename '.xlsx'],[burst1;burst2;jitter_ave(burst1,burst2)]',1,['A' row])
                        jitter_count = jitter_count+1;
                    end
                end
            end
            %Store the overall jitter in an excel spreadsheet
            row = num2str(jitter_count+3);
            xlswrite([filename '.xlsx'],{'Overall Jitter',userData.jitter_overall},1,['A' row])
        end
    end
end


% --- Executes on button press in loadspiketimes.
% Load spike times from an excel sheet in order to perform the jitter analysis
function loadspiketimes_Callback(hObject, eventdata, handles)
% hObject    handle to removespike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.maingui_fig,'UserData');

%Prompt user for an excel that contains spike times
[filename,path] = uigetfile('*.xlsx');

%Get the sheet names of the file
sheet_names = sheetnames([path filename]);
clear spike_time align_spike
for i = 1:length(sheet_names)
    num_spikes = length(readmatrix([path filename],'Sheet',sheet_names(i),'Range','A2:A1000'));
    spike_times(1:num_spikes,i) = readmatrix([path filename],'Sheet',sheet_names(i),'Range','A2:A1000');
    temp_align = readmatrix([path filename],'Sheet',sheet_names(i),'Range','C2:C2');
    align_spike(1,i)=temp_align(1);
end

%Save the spikes and align spike into user data
userData.spike_times = spike_times;
userData.align_spike = align_spike;
set(handles.maingui_fig,'UserData',userData);

%Display the number of bursts for these spike times
set(handles.numbursts_text,'String',size(align_spike,2));

%Clear the data in the table
set(handles.spikejitter_table,'Data',{})

%Clear the data in overall jitter
set(handles.overalljitter_text,'String','')



% --- Executes on button press in done.
% Exit the gui when the user presses done
function done_Callback(hObject, eventdata, handles)
% hObject    handle to find_spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Close out of the GUI
maingui_fig_CloseRequestFcn(handles.maingui_fig, eventdata, handles);



% --- Executes when user attempts to close maingui_fig.
function maingui_fig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to maingui_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


