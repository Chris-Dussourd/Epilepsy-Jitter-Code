%{
Spike Detection and Jitter Analysis
1. Prompts a user for .mat file that contains four variables:
    - channel: an array of the electrical recording data
    - scanFreq: the frequency that the data was sampled
2. It filters the data using a 4th order bandpass filter and a low pass filter
3. The spike detection GUI is displayed. Users can zoom into seizures and select seizure spike times and an alignment spike.
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


%Call the spike detection GUI.
[spike_time, peak_value,align_spike] = selectspikes_gui(channel_filtered,scanFreq);

%Call the jitter analysis GUI.
jitteranalysis_gui(spike_time,align_spike);

