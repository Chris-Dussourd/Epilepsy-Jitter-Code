
% findspikes.m finds the peak points of spikes using the threshold passed in
%   [spike_time, peak_value] = findspikes(data_y,scanFreq,upper_threshold,lower_threshold,min_spikeperiod,max_spikeperiod,min_ISI)   
%           spike_time - function returns the time points that spikes occurred (each column represents spike times for a different burst/time period)
%           peak_value - function returns the peak values of the spikes (each column represents peak values for a different burst/time period)
%           data_y - vector of data points 
%           upper_threshold - used to locate possible spike times (points above this threshold might be spikes)
%           lower_threshold - used to determine the spike period of the spike
%           min_spikeperiod - minimum allowed spike period to be considered a spike
%           max_spikeperiod - maximum allowed spike period to be considered a spike
%           min_ISI - minimum inter-spike interval allowed between spikes (larger spike is taken for close spikes)

function [spike_time, peak_value] = findspikes(data_y,scanFreq,upper_threshold,varargin)

spike_time = [];
peak_value = [];

lower_threshold = upper_threshold/2; %Default lower threshold to half of upper
min_spikeperiod = 5; %Default min spike period to 5ms
max_spikeperiod = 150; %Default max spike period to 150ms
min_ISI = 100; %Default min inter-spiker interval to 100ms

switch nargin
    case 3
    case 4
        lower_threshold = varargin{1};
    case 5
        lower_threshold = varargin{1};
        min_spikeperiod = varargin{2};
    case 6
        lower_threshold = varargin{1};
        min_spikeperiod = varargin{2};
        max_spikeperiod = varargin{3};
    case 7
        lower_threshold = varargin{1};
        min_spikeperiod = varargin{2};
        max_spikeperiod = varargin{3};
        min_ISI = varargin{4};
    otherwise
        error(message('Too many input arguments'))
end

data_x = 1:length(data_y);
%Find the locations where the data crosses the threshold
temp1 = data_y(1:(end-1))-upper_threshold;
temp2 = data_y(2:end)-upper_threshold;
%Since we subtracted the threshold, the zero crossing is where the bursts pass the threshold
threshold_reached = data_x(temp1.*temp2 < 0);
%Now find cossing points for the lower_threshold
temp1 = data_y(1:(end-1))-lower_threshold;
temp2 = data_y(2:end)-lower_threshold;
lowthreshold_reached = data_x(temp1.*temp2<0);

if ~isempty(threshold_reached)
    %Indices
    i = 1;  %i is the index for time_point arrays
    j = 1;  %j is the index for spike times
    nearest_crossing = 1;
    nearest_crossing_prev = 1;
    while i<length(threshold_reached)-1
        %time_points contain the two threshold crossing points:  the max value between the two points represents the peak
        time_point = [threshold_reached(i) threshold_reached(i+1)];
        %If the time_points are not too far or too close to be considered a spike (<max_spikeperiod & >1ms)
        %and the time points enclose a local maximum (not a local minimum), store as a spike
        if time_point(2)-time_point(1)<scanFreq*max_spikeperiod/1000 && time_point(2)-time_point(1)>scanFreq/1000 &&...
                mean(data_y(time_point(1)+1:time_point(2)))>upper_threshold && mean(data_y(time_point(1)+1:time_point(2)))>upper_threshold
            %Find where the peak value occurs - label as spike time
            [peak_value(j),spike_time(j)] = max(data_y(time_point(1):time_point(2)));
            spike_time(j) = spike_time(j)+time_point(1)-1;
            
            %Find the closest lower threshold point before time point 1
            nearest_crossing_temp = lowthreshold_reached(lowthreshold_reached<=time_point(1));
            if ~isempty(nearest_crossing_temp)
                nearest_crossing = nearest_crossing_temp(end);
            end
            %Find the closest lower threshold point after time point 2
            nearest_crossing_end_temp = lowthreshold_reached(lowthreshold_reached>=time_point(2));
            if ~isempty(nearest_crossing_end_temp)
                nearest_crossing_end = nearest_crossing_end_temp(1);
            else
                nearest_crosssing_end = length(data_y);
            end
            
            %If the spike times are too close to be considered one spike (<min_ISI), take the max spike value
            if j>1 && spike_time(j)-spike_time(j-1)<scanFreq*min_ISI/1000
                if peak_value(j) > peak_value(j-1)
                    %Replace previous spike with current spike if current spike has higher amplitude
                    spike_time(j-1) = spike_time(j);
                    peak_value(j-1) = peak_value(j);
                    spike_time(j) = [];
                    peak_value(j) = [];
                    j = j-1;
                else
                    %Remove current spike if previous spike has higher amplitude
                    spike_time(j) = [];
                    peak_value(j) = [];
                    j = j-1;
                end
            %Checks if spike is part of same burst as previous spike (reaches zero the lower threshold at the same point)
            elseif j>1 && nearest_crossing == nearest_crossing_prev
                %Remove one of the spikes since they are part of the same burst
                if peak_value(j) > peak_value(j-1)
                    spike_time(j-1) = spike_time(j);
                    peak_value(j-1) = peak_value(j);
                    spike_time(j) = [];
                    peak_value(j) = [];
                    j = j-1;
                else
                    spike_time(j) = [];
                    peak_value(j) = [];
                    j = j-1;
                end
            %If time period of entire spike is too long (>max_spikeperiod) or too short (<min_spikeperiod), remove spike
            elseif nearest_crossing_end-nearest_crossing>scanFreq*max_spikeperiod/1000 || nearest_crossing_end-nearest_crossing<scanFreq*min_spikeperiod/1000
                spike_time(j) = [];
                peak_value(j) = [];
                j = j-1;
            %If there is a larger peak within the burst, then find the largest peak spike time
            elseif peak_value(j)<max(data_y(nearest_crossing:nearest_crossing_end))
                [peak_value(j) spike_time(j)] = max(data_y(nearest_crossing:nearest_crossing_end));
                spike_time(j) = spike_time(j)+nearest_crossing-1;
                i = find(threshold_reached<nearest_crossing_end,1,'last')-1;
            else
                nearest_crossing_prev = nearest_crossing;
            end
            i = i+2;
            j = j+1;
        elseif time_point(2)-time_point(1)<scanFreq/1000 %Time points are less than 1ms apart
            nearest_crossing_end_temp = lowthreshold_reached(lowthreshold_reached>time_point(1));
            if ~isempty(nearest_crossing_end_temp)
                nearest_crossing_end = nearest_crossing_end_temp(1);
                %Due to high frequency components in waveform (bursty activity), the threshold might be crossed several times in one spike
                %Remove threshold points that separate beginning and ending of spike (waveform must not return to lower threshold before spike ends)
                if nearest_crossing_end > threshold_reached(i+2)
                    threshold_reached(i+1) = threshold_reached(i);
                end
            else
                nearest_crossing_end = length(data_y);
            end
            i = i+1;
        else
            i = i+1;
        end
        
    end
end

end

