%{
Spike Detection and Jitter Analysis
1. Prompts a user for .mat file that contains two variables:
    - channel: an array of the electrical recording data
    - scanFreq: the frequency that the data was sampled
2. It filters the data using a 4th order bandpass filter and a low pass filter
3. Prompts user for .xlsx that contains spike time data. Optional.
3. The spike detection GUI is run and displayed. Users can zoom into seizures and select seizure spike times and an alignment spike.
4. After the user selects done with spike detection, the user can then perform jitter analysis on the data from the spike detection GUI. 
   They can also browse to new spike time data (in the form of excel files).
%}
clc;

low_cut = 3;    %Low Pass Filter
    
%Prompt user for a .mat file that contains the electrical recording and scan frequency
[filename,path] = uigetfile('*.mat');
load([path filename])

%Filter the data
channel_filtered = filter_data(channel,scanFreq,low_cut);

%Default spike_time and align_spike data to empty
spike_times=[];
align_spike=[];

%Ask the user if they want to provide an excel file of spike time data
answer = questdlg('Would you like to provide spike time data for this recording? Answer "No" if you do not have spiket time data for this recording yet.','Provide Spike Time Data','Yes','No','No');

%Retrieve spike time data from the excel file
if strcmp(answer,'Yes')
    %Prompt user for a .xlsx file that contains the spike time data
    [filename,path] = uigetfile('*.xlsx');
    
    %Get the version number of MATLAB. If before R2019a, use xlsread instead of readmatrix
    mat_version = version;
    mat_version_year = str2num(mat_version(find(mat_version=='R')+1:find(mat_version=='R')+4));

    if mat_version_year >= 2019
        %Get the sheet names of the file
        sheet_names = sheetnames([path filename]);
        %Loop through each sheet (burst) and get the spike time and align spike
        for i = 1:length(sheet_names)
            num_spikes = length(readmatrix([path filename],'Sheet',sheet_names(i),'Range','A2:A1000'));
            spike_times(1:num_spikes,i) = readmatrix([path filename],'Sheet',sheet_names(i),'Range','A2:A1000');
            temp_align = readmatrix([path filename],'Sheet',sheet_names(i),'Range','C2:C2');
            align_spike(1,i)=temp_align(1);
        end
    else
        %Get the sheet names of the file
        [~,sheets] = xlsfinfo([path filename]);
        num_sheets = numel(sheets);
        j=1;
        %Loop through each sheet (burst) and get the spike time and align spike
        for i = 1:num_sheets
            num_spikes = length(xlsread([path filename],i,'A2:A1000'));
            spike_times(1:num_spikes,j) = xlsread([path filename],i,'A2:A1000');
            align_spike(1,j) = xlsread([path filename],i,'C2');
            j = j+1;
        end
    end
end


%Call the spike detection GUI.
[spike_times, peak_value,align_spike] = selectspikes_gui(channel_filtered,scanFreq,spike_times,align_spike);

%Call the jitter analysis GUI.
jitteranalysis_gui(spike_times,align_spike);

