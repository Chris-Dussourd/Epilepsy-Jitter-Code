%Plots arrows at spike time points - visual inspection of spike detection algorithms
%   data - entire neural recording or waveform to be plotted
%   spike_times - spike positions in the data set
%   peak_values - the peak voltage values at the spike times
%   scanFreq - sampling frequency when recording data
%   varargin{1} - axes of figure to plot peak points in if provided
function plot_peakpoints(spike_times,peak_values,varargin)

if nargin == 3
    %Figure is already provided
    axes(varargin{1})
    color = 'k'; %Default to black arrows
elseif nargin == 4
    if length(varargin{1})>1
        %Data and sampling frequency are provided
        data = varargin{1};
        scanFreq = varargin{2};
        figure;
        plot((1:length(data))/scanFreq,data);
        color = 'k';
    else
        axes(varargin{1})
        color = varargin{2};
    end
elseif nargin == 5
    data = varargin{1};
    scanFreq = varargin{2};
    figure;
    plot((1:length(data))/scanFreq,data);
    color = varargin{3};
end
axis_position = get(gca,'Position');
xlimit = get(gca,'xlim');
ylimit = get(gca,'ylim');

%The sign determines direction the arrow points
sign = ones(length(spike_times),1);
sign(peak_values<0) = -1;

for i = 1:length(spike_times)
    %Plot an arrow at each of the spike times
    xAnnotation = axis_position(1)+[(spike_times(i)-xlimit(1))/(xlimit(2)-xlimit(1)) (spike_times(i)-xlimit(1))/(xlimit(2)-xlimit(1))]*axis_position(3);
    yAnnotation = axis_position(2)+[(peak_values(i)-ylimit(1))/(ylimit(2)-ylimit(1))+sign(i)/20 (peak_values(i)-ylimit(1))/(ylimit(2)-ylimit(1))]*axis_position(4);
    if xAnnotation(1) >= 0 && xAnnotation(1) <= 1 && yAnnotation(2) >= 0 && yAnnotation(1) <= 1
        hAnnotation = handle(annotation('arrow',xAnnotation,yAnnotation,'Tag','Peak Arrows','Color',color));
    end
%     hAnnotation.pinAtAffordance(1);
%     hAnnotation.pinAtAffordance(2);

end