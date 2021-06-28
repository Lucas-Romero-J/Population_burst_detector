%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Software developed by Javier Lucas-Romero for the paper: (Título)
% July 2021
% This script is designed to analyse sorted spike trains exported from spike2 in order to detect
% synchronic population activity
% The input file architecture is critical to avoid an interrupted execution.
% Example files are provided in the repositoy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
close all
[filename,folderpath1]=uigetfile('*.mat','Select a Matlab file exported from Spike 2');
spike_trains=struct2cell(load([folderpath1 filename])); % File containing the spike trains


number_units=size(spike_trains,1)-1; % amount of neurones in the experiment


tb=0.2; % bin size for the mean frequency plot

%Sometimes the first element in the cell does not contain time values and has to be deleted

try
    spike_trains{number_units+1,1}.times;
    cleaning=1;
catch
    cleaning=2;
end
if cleaning==1
    spike_trains(1)=[];
end

% In this section the spike trains are collapsed into a single vector and sorted
accum=[];
for n=1:number_units
    accum=[accum;spike_trains{n,1}.times];
end
accum=sort(accum);

% Some events can appear at the exact same time point
% To avoid mistakes they are moved 0.0001 s
while length(accum)~=length(unique(accum))
    counter_aux=0.0001;
    for n=1:length(accum)-1
        for i=n+1:length(accum)
            if accum(n)==accum(i)
                accum(i)=accum(i)+counter_aux;
                counter_aux=counter_aux+0.0001;
            else
                break
            end
        end
        counter_aux=0.0001;
    end
end

% The file lenght is stimated as the last spike recorded
T_max=max(accum);


intervals=0:tb:(round(T_max)+2); % The +2 is added to ensure every time point is included
accum_intervals=cell(length(intervals)-1,1);
% In this cell the events will be grouped according to the interval they belong

% This will be used to define a time axis for the mean frequency plot
accum_time_axis=accum;

for m=1:length(intervals)-1
    cond=accum>=intervals(m)&accum<intervals(m+1); % cond and variants will be used subsequently to define logical vectors used to extract values
    if sum(cond)==0 % If there are not events in the interval, at least the time value has to be stored to associate that interval a frequency=0
        accum_time_axis=[accum_time_axis;intervals(m)];
    end
    interval_events=accum(cond);
    accum_intervals{m,1}=interval_events;
end
accum_time_axis=sort(accum_time_axis);

% For each event or empty interval, the associated  mean frequency
% This section transcribes the Spike2 mean frequency formula
Freq_values=[];
for j=1:size(accum_intervals)
    current_interval=accum_intervals{j,1};
    if isempty(current_interval) %If there are not events the frequency value=0
        Freq_values=[Freq_values,0];
        continue
    end
    tl=current_interval(1); %tl=time of the first event in the time range
    n=length(current_interval); % n=events in the time range
    for i=1:length(current_interval)
        %te=time of the current event
        te=current_interval(i);
        if (te-tl)>tb/2
            freq=(n-1)/(te-tl);
        else
            freq=n/tb;
        end
        Freq_values=[Freq_values,freq];
    end
end

% The time axis must start at 0 with an associated frequency value of 0
if accum_time_axis(1)~=0
    accum_time_axis=[0;accum_time_axis];
    Freq_values=[0,Freq_values];
end

% Upsampling at 20000 Hz by interpolation. Will be used in the future for DC correlation
interp_x=0:(1/20000):T_max;
interp_y=interp1(accum_time_axis,Freq_values,interp_x);
interp_x(end)=[];
interp_y(end)=[];
interp_y=interp_y';

% Plot representing the mean frequency and individual units firing
hold on
plot(interp_x,interp_y)
for n=1:number_units
    current_unit=spike_trains{n,1}.times;
    A=zeros(1,length(current_unit));
    A(1,:)=max(interp_y)+2*n;
    plot(spike_trains{n,1}.times,A,'*')
end

hold off % Graph visualitation is needed to determine the manual threshold
% The suggested value for the example file is 75 Hz
threshold_accum=input('Select detection threshold (Hz): ');
threshold_cond_acum=interp_y>threshold_accum;
times_above_threshold_accum=interp_x(threshold_cond_acum); % Time values in the graph with an associated frequency above the threshold
times_above_threshold_accum2=times_above_threshold_accum;% This will be used to extract the maximum firing concentration for each event


% The next section allows the detection of the rising phase of the events crossing the threshold that are separated at least 500 ms
for n=1:length(times_above_threshold_accum)
    try
        times_above_threshold_accum=[times_above_threshold_accum(1:n),times_above_threshold_accum(times_above_threshold_accum>times_above_threshold_accum(n)+0.5)];
    catch
        break
    end
end

% In this section the maximum firing concentration is detected according to mean frequency plot or firing histograms
% Several logical vectors named as condition(n) will be used
amplitudes_above_threshold=interp_y(threshold_cond_acum)';
maximums_mean_freq=[];
for n=1:length(times_above_threshold_accum)
    % condition1: times associated with freq values above threshold in a 0.5 s interval from the detection point
    condition1=times_above_threshold_accum2>times_above_threshold_accum(n)&times_above_threshold_accum2<=times_above_threshold_accum(n)+0.5;
    if sum(condition1)==0 % If the threshold only cuts the event at the peak there would be a single point
        continue
    end
    % Freq values above threshold in the 0.5 s interval from the event detection point
    amplitudes_interval=amplitudes_above_threshold(condition1);
    maximum=max(amplitudes_interval); % local maximum in the interval
    % condition2: Freq values equal to the local maximum
    condition2=amplitudes_above_threshold==maximum;
    % condition3: Freq values in the current interval equal to the local maximum
    condition3=condition1.*condition2;
    condition3=condition3==1; % Conversion to logic vector
    cleaning=find(condition3); % If there are more than 1 maximum, only the first one is considered
    condition3(cleaning(2:length(cleaning)))=0;
    if sum(condition3)>1
        display('error')
        return
    end
    % Time value associated to the detected maximum frequency
    maximums_mean_freq=[maximums_mean_freq, times_above_threshold_accum2(condition3)];
end

% Plot for the accumulate firing histogram
edges2=(0:0.02:T_max);
hist_temp=hist(accum,edges2);
figure
hist(accum,edges2);

% Definition of the intervals around the maximum values of mean frequency
hist_intervals=[];
for n=1:length(maximums_mean_freq)
    hist_intervals(n,1)=maximums_mean_freq(n)-0.25;
    hist_intervals(n,2)=maximums_mean_freq(n)+0.25;
end

% Maximums detection using histogram plot
counter_aux_2=0;
maximums_histogram=[];
for n=1:length(hist_intervals)
    % condition4: time points included in the previously defined intervals
    condition4=edges2>=hist_intervals(n,1)&edges2<=hist_intervals(n,2);
    max_local_hist=max(hist_temp(condition4)); % maximum firing concentration point for the current event
    % condition5: values in the general histogram that match the detected maximum
    condition5=hist_temp==max_local_hist;
    % condition6: position of the detected maximum in the current interval
    condition6=condition4.*condition5;
    condition6=condition6==1; % transformation into a logtic vector
    if sum(condition6)>1
        cleaning=find(condition6); % If there are more than 1 maximum in the interval, only the first one is considered
        condition6(cleaning(2:length(cleaning)))=0;
        counter_aux_2=counter_aux_2+1;
    end
    maximums_histogram=[maximums_histogram, edges2(condition6)];
end


% This line sets the reference point for the following processing
% In this work the reference point is the maximum values detected in the histogram plot
times_above_threshold_accum=maximums_histogram;

% Definition of interval limits for the following processing
left_limit=0.05;
right_limit=0.150;

population_burst_intervals=[];
for n=1:length(times_above_threshold_accum)
    population_burst_intervals(n,1)=times_above_threshold_accum(n)-left_limit;
    population_burst_intervals(n,2)=times_above_threshold_accum(n)+right_limit;
end

% Security step
if size(population_burst_intervals,1)==0
    return
end

% Percentage of spikes falling within synchronic events
percentages_units=zeros(number_units,1);
% Spike distribution around population bursts for each unit
population_spikes_per_unit=cell(number_units,1);

for n=1:number_units
    current_unit=spike_trains{n,1}.times;
    events_accum=[];
    spike_counter=0;
    for i=1:length(current_unit)
        if sum(current_unit(i)>=population_burst_intervals(:,1)&current_unit(i)<=population_burst_intervals(:,2))>=1
            % Chech if each neurone spike is within the defined intervals.
            spike_counter=spike_counter+1;
        end
    end
    percentages_units(n,1)=spike_counter/length(current_unit);
    
    % Generation of spike distributions around population burst for each neurone
    for j=1:size(population_burst_intervals,1)
        firings_current_interval=current_unit(current_unit>=(population_burst_intervals(j,1))&current_unit<(population_burst_intervals(j,2)));
        if j>=2
            % The next line avoids the duplication of spikes if the intervals are overlapping
            firings_current_interval=firings_current_interval(firings_current_interval>population_burst_intervals(j-1,2));
        end
        % The left value of the interval and 0.05 is substracted to have the 0 as the reference point
        firings_current_interval=firings_current_interval-population_burst_intervals(j,1)-left_limit;
        events_accum=[events_accum;firings_current_interval];
    end
    population_spikes_per_unit{n,1}=events_accum';
end


% Output file with the raw times for the distribution
% Each row represents the times for each unit
% If several experiments are analysed sequentially the new units are added in the successive rows
output_file=fopen('output_file.csv','a');
for n=1:number_units
    fprintf(output_file,'%f, ',population_spikes_per_unit{n,1});
    fprintf(output_file,'\n');
end
fclose(output_file);

% Histograms for previsualization of spike distribution per neurone
edges=(-left_limit:0.005:right_limit);
H_neurones=[];
for n=1:number_units
    H_neurones=[H_neurones;hist(population_spikes_per_unit{n,1},edges)];
end

% Data for optional output containing an experiment summary per neurone
num_events=[];
experiment_duration=[];
num_events(1:size(H_neurones,1),1)=length(times_above_threshold_accum); % Total number of detected population bursts
experiment_duration(1:size(H_neurones,1),1)=T_max; % Approximate experiment duration
H_neurones_per_burst=H_neurones./num_events; % For each histogram interval, spikes per population burst
Output_per_neurone=[num_events,experiment_duration,H_neurones_per_burst];

% Data for optional output containing an experiment summary
H_total_exp=mean(H_neurones_per_burst); % Mean distribution histogram for the experiment
event_frequenzy=num_events(1,1)/T_max; % Event frequenzy in Hz
mean_burst_duration=0.005*sum(H_total_exp>(max(H_total_exp)/2)); % Population burst duration at mid amplitude (approximation limited by bin size)
total_spikes_in_bursts=sum(sum(H_neurones(:,:))); % Amount of spikes falling within population bursts
spikes_per_bursts=total_spikes_in_bursts/num_events(1,1); % Mean amount of spikes falling within population bursts
spikes_per_bursts_and_neurone=spikes_per_bursts/number_units; % Mean amount of spikes falling within population bursts per unit
Experiment_summary=[num_events(1,1),event_frequenzy,mean_burst_duration,number_units,spikes_per_bursts,spikes_per_bursts_and_neurone,H_total_exp];

