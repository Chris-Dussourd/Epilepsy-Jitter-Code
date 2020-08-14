
% selectspikes_gui.m is a user interface to find spikes in bursts 
%   [spike_time, peak_value] = selectspikes_gui(data,scanFreq,upper_threshold,lower_threshold,min_spikeperiod,max_spikeperiod,min_ISI)   
%           spike_time - returns the time points that spikes occurred (each column represents spike times for a different burst/time period)
%           peak_value - returns the peak values of the spikes (each column represents peak values for a different burst/time period)
%           data - vector of voltage data points 
%           upper_threshold - used to locate possible spike times (points above this threshold might be spikes)
%           lower_threshold - used to determine the spike period of the spike
%           min_spikeperiod - minimum allowed spike period to be considered a spike
%           max_spikeperiod - maximum allowed spike period to be considered a spike
%           min_ISI - minimum inter-spike interval allowed between spikes (larger spike is taken for close spikes)

%   [spike_time, peak_value] = selectspikes_gui(data,scanFreq,current_spikes,align_spike,varargin)
%           current_spikes - spike times provided by the user to be edited in the gui
%                           - if current_spikes is passed in, it must always be the 3rd variable %%%%important
%           align_spike - the alignment spike provided by the user (to align seizures for plotting)
%                           - alignment spike must always be the 4th variable if current_spikes is passed in %%%%important
%                           - current spikes and align spike are in units of seconds

function varargout = selectspikes_gui(varargin)
% SELECTSPIKES_GUI MATLAB code for selectspikes_gui.fig
%      SELECTSPIKES_GUI, by itself, creates a new SELECTSPIKES_GUI or raises the existing
%      singleton*.
%
%      H = SELECTSPIKES_GUI returns the handle to a new SELECTSPIKES_GUI or the handle to
%      the existing singleton*.
%
%      SELECTSPIKES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTSPIKES_GUI.M with the given input arguments.
%
%      SELECTSPIKES_GUI('Property','Value',...) creates a new SELECTSPIKES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before selectspikes_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to selectspikes_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help selectspikes_gui

% Last Modified by GUIDE v2.5 09-Aug-2020 16:51:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @selectspikes_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @selectspikes_gui_OutputFcn, ...
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

% --- Executes just before selectspikes_gui is made visible.
function selectspikes_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to selectspikes_gui (see VARARGIN)

% Choose default command line output for selectspikes_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Data and sampling frequency must be provided
data = varargin{1};
scanFreq = varargin{2};
upperthresh = 0.1;  %Default upper threshold to 100uV
lowerthresh = 0.05; %Default lower threshold to 50uV
min_spikeperiod = 5; %Default min spike period to 5ms
max_spikeperiod = 150; %Default max spike period to 150ms
min_ISI = 100; %Default min inter-spiker interval to 100ms
sign = 1; %Default sign to positive (sign = -1  means negative going spikes)

%If any spike detection parameters are passed in set
%Note: nargin is 3 more than number of variables passed in for MATLAB guis
if nargin >= 6
    if length(varargin{3}) == 1
        upperthresh = varargin{3};
        if nargin >= 7
            lowerthresh = varargin{4};
        end
        if nargin >= 8
            min_spikeperiod = varargin{5};
        end
        if nargin >= 9
            max_spikeperiod = varargin{6};
        end
        if nargin == 10
            min_ISI = varargin{7};
        end
        current_spikes = []; %User didn't provide spikes to edit
    else
        %If length of 3rd variable passed in is a vector - spike times are 
        %   provided by the user to be edited in the gui
        current_spikes = varargin{3};
        if nargin >= 7
            if length(varargin{4}) == 1
                upperthresh = varargin{4};
                if nargin >= 8
                    lowerthresh = varargin{5};
                end
                if nargin >= 9
                    min_spikeperiod = varargin{6};
                end
                if nargin >= 10
                    max_spikeperiod = varargin{7};
                end
                if nargin == 11
                    min_ISI = varargin{8};
                end
                align_spike = []; %User didn't provide alignment spikes
            else
                %If length of 4th variable passed in is a vector - align _spikes 
                %  are provided by the user to be edited by gui
                align_spike = varargin{4};
                if nargin >= 8
                    upperthresh = varargin{5};
                end
                if nargin >= 9
                    lowerthresh = varargin{6};
                end
                if nargin >= 10
                    min_spikeperiod = varargin{7};
                end
                if nargin >= 11
                    max_spikeperiod = varargin{8};
                end
                if nargin == 12
                    min_ISI = varargin{9};
                end
            end
        end
    end
else
    current_spikes = [];  %User didn't provide spikes to edit
end
        
%Initialize userdata in figure as empty
set(handles.maingui_fig,'UserData',[]);
%Store the data parameters
userData.data = data;
userData.scanFreq = scanFreq;

%Save the parameters for automatic spike detection
userData.upperthresh = upperthresh;
userData.lowerthresh = lowerthresh;
userData.min_spikeperiod = min_spikeperiod;
userData.max_spikeperiod = max_spikeperiod;
userData.min_ISI = min_ISI;
userData.sign = sign;
%Set the parameters for spike detection in edit texts of gui
set(handles.upperthresh,'String',upperthresh)
set(handles.lowerthresh,'String',lowerthresh)
set(handles.min_spikeperiod,'String',min_spikeperiod)
set(handles.max_spikeperiod,'String',max_spikeperiod)
set(handles.min_ISI,'String',min_ISI)
if sign == 1
    set(handles.positive_sign,'Value',get(handles.positive_sign,'Max'));
else
    set(handles.negative_sign,'Value',get(handles.negative_sign,'Max'));
end

%Initialize saved variable (used to determine if analysis was saved
userData.saved = 0;

%Initialize spike time and peak value
userData.burst = 1;
if length(current_spikes)> 1 %Spike times have been provided by the user
    userData.spike_time = current_spikes;
    for i = 1:size(current_spikes,2)
        temp_spikes = current_spikes(:,i);
        temp_spikes = temp_spikes(temp_spikes~=0);
        %Find and store the peak values
        userData.peak_value(1:length(temp_spikes),i) = data(round(temp_spikes*scanFreq));
        if length(align_spike) < i
            userData.align_spike(:,i) = current_spikes(1,i); %Default the alignment spike to the first spike
        else
            userData.align_spike(:,i) = align_spike(:,i);
        end
        %Add burst to listbox
        listbox_string = get(handles.burst_list,'String');
        listbox_string = [listbox_string; {['Burst:  ' num2str(floor(userData.spike_time(1,userData.burst))) ' seconds']}];
        set(handles.burst_list,'String',listbox_string,'Value',[]);
        
        userData.burst = userData.burst+1;
    end
    userData.burst = size(current_spikes,2)+1;
else
    userData.spike_time = []; 
    userData.peak_value = [];
    userData.align_spike = [];
end
%Set the data stored into the UserData of the figure
set(handles.maingui_fig,'UserData',userData);
%Plot the data in the gui axes (axes are named "burst_recording")
axes(handles.burst_recording)
plot((1:length(data))/(scanFreq),data)
xlabel('Time (sec)')
ylabel('Recording Voltage (mV)')
zoom reset; %Allows zoom in/out of current axes

% UIWAIT makes selectspikes_gui wait for user response (see UIRESUME)
uiwait(handles.maingui_fig);


% --- Outputs from this function are returned to the command line.
function varargout = selectspikes_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
userData = get(handles.maingui_fig,'UserData');
varargout{1} = userData.spike_time;
varargout{2} = userData.peak_value;
varargout{3} = userData.align_spike;
delete(handles.maingui_fig)

% --- Executes on button press in addspike.
%Adds the spike selected using the data cursor to spike times
function addspike_Callback(hObject, eventdata, handles)
% hObject    handle to addspike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.burst_recording)

userData = get(handles.maingui_fig,'UserData');
burst = userData.burst;
%Get the data point from the cursor
[spike_time_temp,peak_value_temp] = ginput(1);

%Add the data point to the spike time and peak value array (and sort)
spike_time = userData.spike_time(:,burst);
peak_value = userData.peak_value(:,burst);
peak_value = [peak_value(spike_time~=0); peak_value_temp];
[spike_time, sort_index] = unique([spike_time(spike_time ~= 0); spike_time_temp]);
peak_value = peak_value(sort_index);

%Add an annotation for the new spike time added
plot_peakpoints(spike_time_temp,peak_value_temp,handles.burst_recording)

%Store the spike times and peak values
userData.spike_time(1:length(spike_time),burst) = spike_time;
userData.peak_value(1:length(peak_value),burst) = peak_value;
set(handles.maingui_fig,'UserData',userData);


% --- Executes on button press in removespike.
% Removes the nearest spike time to current location of data cursor
function removespike_Callback(hObject, eventdata, handles)
% hObject    handle to removespike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.burst_recording)

userData = get(handles.maingui_fig,'UserData');
burst = userData.burst;
%Get the data point from a user click
[spike_time_temp,~] = ginput(1);

%Find data point to remove (nearest spike time to selected time)
spike_time = userData.spike_time(:,burst);
peak_value = userData.peak_value(:,burst);
[~, remove_index] = min(abs(spike_time-spike_time_temp));

%Remove the annotation for the spike
remove_peakpoints(spike_time(remove_index),peak_value(remove_index),handles.burst_recording)

%Remove the data point with nearest spike time 
spike_time(remove_index) = [];
peak_value(remove_index) = [];
spike_time = [spike_time; 0];
peak_value = [peak_value; 0];

%Store the spike times and peak values
userData.spike_time(1:length(spike_time),burst) = spike_time;
userData.peak_value(1:length(peak_value),burst) = peak_value;
set(handles.maingui_fig,'UserData',userData);


% --- Executes on button press in find_spikes.
% Use automatic spike detection algorithm to find spikes (clears any previous spikes not saved)
function find_spikes_Callback(hObject, eventdata, handles)
% hObject    handle to find_spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(findall(gcf,'Tag','Peak Arrows'))

%Get the data points and xlimits of the burst in the axes
userData = get(handles.maingui_fig,'UserData');
xlimit = get(handles.burst_recording,'xlim');
scanFreq = userData.scanFreq;
data = userData.data;
data_x = (1:length(data))/(scanFreq);
ydata = data(data_x > xlimit(1));
xdata = data_x(data_x > xlimit(1));
ydata = ydata(xdata < xlimit(2));

%Use the algorithm to find the spike times in burst interval
[spike_time, peak_value] = findspikes(ydata*userData.sign,scanFreq,userData.upperthresh,userData.lowerthresh,...
    userData.min_spikeperiod,userData.max_spikeperiod,userData.min_ISI);
spike_time = (spike_time-1)'/(scanFreq)+xdata(1);
peak_value = peak_value'*userData.sign;
%Save the spike times and peak values in userdata 
userData.spike_time(1:max([length(userData.spike_time) length(spike_time)]),userData.burst) = [spike_time; zeros(max([length(userData.spike_time)-length(spike_time) 0]),1)];
userData.peak_value(1:max([length(userData.spike_time) length(spike_time)]),userData.burst) = [peak_value; zeros(max([length(userData.spike_time)-length(peak_value) 0]),1)];
set(handles.maingui_fig,'UserData',userData);
%Plot the spike times on the figure
plot_peakpoints(spike_time,peak_value,handles.burst_recording)


% --- Executes on button press in add_spikes.
% Use automatic spike detection algorithm to add spikes
function add_spikes_Callback(hObject, eventdata, handles)
% hObject    handle to add_spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(findall(gcf,'Tag','Peak Arrows'))

%Get the data points and xlimits of the burst in the axes
userData = get(handles.maingui_fig,'UserData');
xlimit = get(handles.burst_recording,'xlim');
scanFreq = userData.scanFreq;
data = userData.data;
data_x = (1:length(data))/(scanFreq);
ydata = data(data_x > xlimit(1));
xdata = data_x(data_x > xlimit(1));
ydata = ydata(xdata < xlimit(2));

%Use the algorithm to find the spike times in burst interval
[new_spike_time, new_peak_value] = findspikes(ydata*userData.sign,scanFreq,userData.upperthresh,userData.lowerthresh,...
    userData.min_spikeperiod,userData.max_spikeperiod,userData.min_ISI);
new_spike_time = (new_spike_time-1)'/(scanFreq)+xdata(1);
new_peak_value = new_peak_value'*userData.sign;
%Save the spike times and peak values in userdata 
%Add the data point to the spike time and peak value array (and sort)
if size(userData.spike_time)>=userData.burst
    spike_time = userData.spike_time(:,userData.burst);
    peak_value = userData.peak_value(:,userData.burst);
    peak_value = [peak_value(spike_time~=0); new_peak_value];
    %Sort and remove duplicates
    [spike_time, sort_index] = unique([spike_time(spike_time ~= 0); new_spike_time]);
    peak_value = peak_value(sort_index);
else
    spike_time = new_spike_time;
    peak_value = new_peak_value;
end

%Remove spike times that are less than the min inter-spike interval
i = 2;
while i<length(spike_time)
    if spike_time(i)-spike_time(i-1) < userData.min_ISI/1000
        %Take the largest peak value
        if abs(peak_value(i)) >= abs(peak_value(i-1))
            spike_time = [spike_time(1:i-2); spike_time(i:end)];
            peak_value = [peak_value(1:i-2); peak_value(i:end)];
            i=i-1;
        else
            spike_time = [spike_time(1:i-1); spike_time(i+1:end)];
            peak_value = [peak_value(1:i-1); peak_value(i+1:end)];
            i=i-1;
        end
    end
    i=i+1;
end
userData.spike_time(1:max([length(userData.spike_time) length(spike_time)]),userData.burst) = [spike_time; zeros(max([length(userData.spike_time)-length(spike_time) 0]),1)];
userData.peak_value(1:max([length(userData.spike_time) length(spike_time)]),userData.burst) = [peak_value; zeros(max([length(userData.peak_value)-length(peak_value) 0]),1)];
set(handles.maingui_fig,'UserData',userData);
%Plot the spike times on the figure
plot_peakpoints(spike_time(spike_time>xlimit(1) & spike_time<xlimit(2)),...
    peak_value(spike_time>xlimit(1) & spike_time<xlimit(2)),handles.burst_recording)


% --- Executes on button press in plot_peaklabels.
function plot_peaklabels_Callback(hObject, eventdata, handles)
% hObject    handle to plot_peaklabels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Remove Peak Arrow Annotations from previous burst
delete(findall(gcf,'Tag','Peak Arrows'))

%Get the data points and xlimits of the burst in the axes
userData = get(handles.maingui_fig,'UserData');
spike_time = userData.spike_time(:,userData.burst);
peak_value = userData.peak_value(:,userData.burst);
xlimit = get(handles.burst_recording,'xlim');

%Plot the spike times on the figure
plot_peakpoints(spike_time(spike_time>xlimit(1) & spike_time<xlimit(2)),...
    peak_value(spike_time>xlimit(1) & spike_time<xlimit(2)),handles.burst_recording)

if size(userData.align_spike,2)>=userData.burst
    align_spike=userData.align_spike(1,userData.burst);
    align_index=spike_time==align_spike;
    %Plot the alignment spike in red
    if align_spike>xlimit(1) && align_spike<xlimit(2)
        plot_peakpoints(align_spike,peak_value(align_index),handles.burst_recording,'r');
    end
end


% --- Executes on button press in saveburst.
% Saves the burst and displays it in the listbox
function saveburst_Callback(hObject, eventdata, handles)
% hObject    handle to saveburst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Remove Peak Arrow Annotations from previous burst
delete(findall(gcf,'Tag','Peak Arrows'))
%Increase burst number by one
userData = get(handles.maingui_fig,'UserData');
if length(userData.align_spike) == userData.burst
    spike_time = userData.spike_time(:,userData.burst);
    peak_value = userData.peak_value(:,userData.burst);
    peak_value = peak_value(spike_time~=0);
    %Remove any duplicates in spike times and peak values
    [spike_time, sort_index] = unique(spike_time(spike_time~=0));
    peak_value = peak_value(sort_index);
    userData.burst = userData.burst+1;
    %Store the spike times and peak values
    userData.spike_time(1:length(spike_time),userData.burst-1) = spike_time;
    userData.peak_value(1:length(peak_value),userData.burst-1) = peak_value;
    set(handles.maingui_fig,'UserData',userData);
    
    %Add burst to listbox
    listbox_string = get(handles.burst_list,'String');
    listbox_string = [listbox_string; {['Burst:  ' num2str(floor(userData.spike_time(1,userData.burst-1))) ' seconds']}];
    set(handles.burst_list,'String',listbox_string,'Value',[]);

else
    msgbox('Please set the alignment spike to align with other bursts before saving burst.')
end


%Plot bursts in order to find the best regions to compare
% --- Executes on button press in comparebursts.
function comparebursts_Callback(hObject, eventdata, handles)
% hObject    handle to comparebursts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Save the bursts selected in the listbox
burst = get(handles.burst_list,'value');
if isempty(burst)
    %Defaults to first burst if none are selected
    burst = 1;
end

userData = get(handles.maingui_fig,'UserData');
if isempty(userData.align_spike)
    align_spike = userData.spike_time(1,:);
else
    align_spike = userData.align_spike;
end
%Uses the same plot length for all bursts
plot_length_1 = max(align_spike(:,burst)-userData.spike_time(1,burst)); %Before alignment spike
plot_length_2 = max(max(userData.spike_time(:,burst))-align_spike(:,burst));     %After alignment spike
data = userData.data;
scanFreq = userData.scanFreq;
data_x = (1:length(data))/(scanFreq);
figure('Units','Normalized','Position',[0 0.1 1 0.9]);

%Make plots of all bursts selected
for i = 1:length(burst)
    %Take out spike times and peak values from stored matrices
    spike_time = userData.spike_time(:,burst(i));
    peak_value = userData.peak_value(:,burst(i));
    peak_value = peak_value(spike_time~=0);
    spike_time = spike_time(spike_time~=0);
    %The limits extend slightly before and after spike times
    xlimit = [align_spike(:,burst(i))-plot_length_1-5 align_spike(:,burst(i))+plot_length_2+5];
    
    %Plot the burst in a subplot with peak points labeled
    h_1 = subplot(max([length(burst) 2]),1,i);
    ydata = data(data_x > xlimit(1));
    xdata = data_x(data_x > xlimit(1));
    ydata = ydata(xdata < xlimit(2));
    xdata = xdata(xdata < xlimit(2));
    plot(xdata,ydata)
    xlim([xlimit(1) xlimit(2)])
    %Label the peak points
    plot_peakpoints(spike_time,peak_value,h_1)
    %Plot the alignment spike in red
    align_index = find(spike_time == align_spike(:,burst(i)));
    plot_peakpoints(spike_time(align_index),peak_value(align_index),h_1,'r');
end
%Compare burst with current axes in figure if only one burst is selected
if length(burst) == 1
    xlimit2 = get(handles.burst_recording,'xlim');
    %Make plots the same length
    xlimit2(2) = xlimit2(1)+xlimit(2)-xlimit(1);

    %Plot burst two in second subplot
    subplot(2,1,2);
    ydata = data(data_x > xlimit2(1));
    xdata = data_x(data_x > xlimit2(1));
    ydata = ydata(xdata < xlimit2(2));
    xdata = xdata(xdata < xlimit2(2));
    plot(xdata,ydata)
    xlim([xlimit2(1) xlimit2(2)])
end


% --- Executes during object creation, after setting all properties.
function burst_recording_CreateFcn(hObject, eventdata, handles)
% hObject    handle to burst_recording (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate burst_recording


%Stores the upper threshold if user changes the edit box
function upperthresh_Callback(hObject, eventdata, handles)
% hObject    handle to upperthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upperthresh as text
%        str2double(get(hObject,'String')) returns contents of upperthresh as a double
userData = get(handles.maingui_fig,'UserData');
userData.upperthresh = str2double(get(hObject,'String'));
set(handles.maingui_fig,'UserData',userData);


% --- Executes during object creation, after setting all properties.
function upperthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upperthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Stores the lower threshold if user changes the edit box
function lowerthresh_Callback(hObject, eventdata, handles)
% hObject    handle to lowerthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowerthresh as text
%        str2double(get(hObject,'String')) returns contents of lowerthresh as a double
userData = get(handles.maingui_fig,'UserData');
userData.lowerthresh = str2double(get(hObject,'String'));
set(handles.maingui_fig,'UserData',userData);


% --- Executes during object creation, after setting all properties.
function lowerthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowerthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Stores the min spike period if user changes the edit box
function min_spikeperiod_Callback(hObject, eventdata, handles)
% hObject    handle to min_spikeperiod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_spikeperiod as text
%        str2double(get(hObject,'String')) returns contents of min_spikeperiod as a double
userData = get(handles.maingui_fig,'UserData');
userData.min_spikeperiod = str2double(get(hObject,'String'));
set(handles.maingui_fig,'UserData',userData);

% --- Executes during object creation, after setting all properties.
function min_spikeperiod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_spikeperiod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Stores the max spike period if user changes the edit box
function max_spikeperiod_Callback(hObject, eventdata, handles)
% hObject    handle to max_spikeperiod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_spikeperiod as text
%        str2double(get(hObject,'String')) returns contents of max_spikeperiod as a double
userData = get(handles.maingui_fig,'UserData');
userData.max_spikeperiod = str2double(get(hObject,'String'));
set(handles.maingui_fig,'UserData',userData);

% --- Executes during object creation, after setting all properties.
function max_spikeperiod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_spikeperiod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Stores the min inter-spike interval if user changes the edit box
function min_ISI_Callback(hObject, eventdata, handles)
% hObject    handle to min_ISI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_ISI as text
%        str2double(get(hObject,'String')) returns contents of min_ISI as a double
userData = get(handles.maingui_fig,'UserData');
userData.min_ISI = str2double(get(hObject,'String'));
set(handles.maingui_fig,'UserData',userData);

% --- Executes during object creation, after setting all properties.
function min_ISI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_ISI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%Stores the filename if the user changes the edit box
function filename_Callback(hObject, eventdata, handles)
% hObject    handle to min_ISI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_ISI as text
%        str2double(get(hObject,'String')) returns contents of min_ISI as a double
userData = get(handles.maingui_fig,'UserData');
userData.filename = str2double(get(hObject,'String'));
set(handles.maingui_fig,'UserData',userData);


% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in burst_list.
function burst_list_Callback(hObject, eventdata, handles)
% hObject    handle to burst_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns burst_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from burst_list

%If item on listbox was double_clicked
if strcmp(get(handles.maingui_fig,'SelectionType'),'open')
    %Remove Peak Arrow Annotations from previous burst
    delete(findall(gcf,'Tag','Peak Arrows'))
    
    %Set the burst selected to the last burst in list
    burst = get(handles.burst_list,'value');
    listbox_string = get(handles.burst_list,'String');
    listbox_string(burst) = [];
    set(handles.burst_list,'String',listbox_string,'Value',length(listbox_string));
    
    %Remove data from a burst previously worked on before editing this burst
    userData = get(handles.maingui_fig,'UserData');
    if size(userData.spike_time,2)==userData.burst
        userData.spike_time(:,userData.burst) = [];
        userData.peak_value(:,userData.burst) = [];
    end
    
    %Set burst selected to last burst in list
    num_bursts = size(userData.spike_time,2);
    userData.spike_time = [userData.spike_time(:,1:burst-1) userData.spike_time(:,burst+1:end) userData.spike_time(:,burst)];
    userData.peak_value = [userData.peak_value(:,1:burst-1) userData.peak_value(:,burst+1:end) userData.peak_value(:,burst)];
    userData.align_spike = [userData.align_spike(:,1:burst-1) userData.align_spike(:,burst+1:end) userData.align_spike(:,burst)];
    userData.burst = userData.burst-1;
    set(handles.maingui_fig,'UserData',userData);
    
    %Display Burst in axes with peak points plotted
    axes(handles.burst_recording)
    spike_time = userData.spike_time(:,num_bursts);
    peak_value = userData.peak_value(:,num_bursts);
    peak_value = peak_value(spike_time~=0);
    spike_time = spike_time(spike_time~=0);
    axis([spike_time(1)-0.5 spike_time(end) -max(abs(peak_value))-min(abs(peak_value)) max(abs(peak_value))+min(abs(peak_value))])
    %Label the peak points
    plot_peakpoints(spike_time,peak_value,handles.burst_recording)
    %Label the alignment spike in red
    align_index = find(spike_time == userData.align_spike(:,num_bursts));
    plot_peakpoints(spike_time(align_index),peak_value(align_index),handles.burst_recording,'r');
end

% --- Executes during object creation, after setting all properties.
function burst_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to burst_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in deleteburst.
function deleteburst_Callback(hObject, eventdata, handles)
% hObject    handle to deleteburst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Delete burst selected from the listbox
burst = get(handles.burst_list,'value');
listbox_string = get(handles.burst_list,'String');
listbox_string(burst(1)) = [];
set(handles.burst_list,'String',listbox_string,'Value',[]);

%Delete burst from userData
userData = get(handles.maingui_fig,'UserData');
userData.spike_time(:,burst(1)) = [];
userData.peak_value(:,burst(1)) = [];
userData.align_spike(:,burst(1)) = [];
userData.burst = userData.burst-1;
set(handles.maingui_fig,'UserData',userData);


% --- Executes when selected object is changed in sign_uipanel.
%Stores the sign if user changes the selection of the radio buttons
function sign_uipanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in sign_uipanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.maingui_fig,'UserData');
if get(handles.positive_sign,'Value') == get(handles.positive_sign,'Max')
    userData.sign = 1;
else
    userData.sign = -1;
end
set(handles.maingui_fig,'UserData',userData);


% --- Executes on button press in alignment_spike.
% Set the spike which aligns the seizures
function alignment_spike_Callback(hObject, eventdata, handles)
% hObject    handle to alignment_spike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.maingui_fig,'UserData');
burst = userData.burst;
%Get the data point from the data cursor
[spike_time_temp,~] = ginput(1);
[~,align_index] = min(abs(spike_time_temp-userData.spike_time(:,burst)));
userData.align_spike(1,burst) = userData.spike_time(align_index,burst);
%Plot the alignment spike in red
h_axes = handles.burst_recording;
plot_peakpoints(userData.spike_time(align_index,burst),userData.peak_value(align_index,burst),h_axes,'r');
set(handles.maingui_fig,'UserData',userData);


% --- Executes on button press in saveanalysis.
% Saves the spike time data to an excel spreadsheet
function saveanalysis_Callback(hObject, eventdata, handles)
% hObject    handle to removespike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.maingui_fig,'UserData');

filename = get(handles.filename,'String');
if isempty(filename) || strcmp(filename,'<Enter File Name Here>')
    msgbox('Please enter a filename to save the spike time data.')
elseif isempty(userData.spike_time)
    msgbox('Please save spike time data for at least one burst before saving to excel.')
else
    answer = 'Yes'; %Default to 'Yes' to create/overwrite a file.
    if isfile([filename '.xlsx'])
        answer = questdlg(['A file with the name "' filename '" already exists. Would you like to overwrite it?'],'Confirm Save','Yes','No','No');
    end
    if strcmp(answer,'Yes')
        if isfile([filename '.xlsx'])
            delete([filename '.xlsx'])
        end
        %If alignment spike is empty, default alignment spike to first spike
        if isempty(userData.align_spike)
            userData.align_spike = userData.spike_time(1,:);
        end
        
        %Sort the spike times before saving it
        [~,sort_index] = sort(userData.spike_time(1,:));
        userData.spike_time = userData.spike_time(:,sort_index);
        userData.align_spike = userData.align_spike(:,sort_index);
        
        index=0;
        
        %Get the version number. If before R2019a, use xlswrite instead of writecell
        mat_version = version;
        mat_version_year = str2num(mat_version(find(mat_version=='R')+1:find(mat_version=='R')+4));
        
        if mat_version_year >= 2019
            %For each burst, save the spike times in a new sheet of an excel file
            for burst = 1:size(userData.spike_time,2)
                %Take out the non-zero spike_times from the matrix
                spike_time1 = userData.spike_time(:,burst);
                spike_time1 = spike_time1(spike_time1~=0);
                index = index+1;
                %Save the seizure spike time data to an excel file
                writecell(['Spike Timings (s)';num2cell(spike_time1)],[filename '.xlsx'],'Sheet',index,'Range','A:A')
                writecell(['Alignment Spike';num2cell(userData.align_spike(:,burst))],[filename '.xlsx'],'Sheet',index,'Range','C:C')
            end
        else
            %For each burst, save the spike times in a new sheet of an excel file
            for burst = 1:size(userData.spike_time,2)
                %Take out the non-zero spike_times from the matrix
                spike_time1 = userData.spike_time(:,burst);
                spike_time1 = spike_time1(spike_time1~=0);
                index = index+1;
                %Save the seizure spike time data to an excel file
                xlswrite([filename '.xlsx'],['Spike Timings (s)';num2cell(spike_time1)],index, 'A1')
                xlswrite([filename '.xlsx'],['Alignment Spike';num2cell(userData.align_spike(:,burst))],index,'C1')
            end
        end
        userData.saved = 1; %We saved data
        set(handles.maingui_fig,'UserData',userData);
    end
end


% --- Executes on button press in done.
% Exit the gui when the user presses done.
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

if isempty(handles)
    answer = 'Yes';
else
    answer = 'Yes'; %Default to yes to exit.
    userData = get(handles.maingui_fig,'UserData');
    if userData.saved == 0 
        answer = questdlg('You have not saved any spike time data. Are you sure you want to exit?','Confirm Close','Yes','No','No');
    end
end
if strcmp(answer,'Yes')
    if isequal(get(hObject, 'waitstatus'), 'waiting')
        % The GUI is still in UIWAIT, us UIRESUME
        uiresume(hObject);
    else
        % The GUI is no longer waiting, just close it
        delete(hObject);
    end
end
