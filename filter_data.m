function [data_filtered] = filter_data(data,sample_frequency,low_cut)
    %processing parameters for band-stop filtering
    band_stop1 = 58; % stop band starts at 55 Hz
    band_stop2 = 62; % stop band stops at 65 Hz
    band_stop = [band_stop1*(1/(0.5*sample_frequency)) band_stop2*(1/(0.5*sample_frequency))];
    [z, p, k] = butter(4, band_stop, 'stop');  % create 4th-order Butterworth filter
    [sos,g] = zp2sos(z,p,k);
%     Hd = dfilt.df2tsos(sos,g);
    band_filtered_channel = filtfilt(sos,g,data);
    
    %processing parameters for low-cut filtering
    [z,p,k] = butter(4, low_cut/(0.5*sample_frequency), 'high');  % create 4th-order Butterworth filter
    [sos,g] = zp2sos(z,p,k);
%     Hd = dfilt.df2tsos(sos,g);
    high_filtered_channel = filtfilt(sos,g,band_filtered_channel);              % filter channel 2
    data_filtered = high_filtered_channel;           
end

