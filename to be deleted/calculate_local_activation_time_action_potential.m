function [lat,cl] = calculate_local_activation_time_action_potential(t_start,action_potential,cycle_length_percentage)

% find out cycle length
threshold = 0.13;
markers = zeros(size(action_potential,2),1);
cycle_lengths = zeros(size(action_potential,2),1);
for n = 1:size(action_potential,2)
    temp = find_action_potential_marker(action_potential(:,n),threshold);
    cl = diff(temp);
    
    markers(n,1:length(temp)) = temp;
    cycle_lengths(n,1:length(cl)) = cl;
end
markers(markers==0) = NaN;
cycle_lengths(cycle_lengths==0) = NaN;

debug_plot = 0;
if debug_plot == 1
    cl = cycle_lengths(:);
    figure;
    histogram(cl,10);
    xlabel('cycle length, ms');
    ylabel('counts');
    title('cycle length histogram');
end

cl = cycle_lengths(:);
cl(isnan(cl)) = [];
CL = mean(cl);

if isnan(CL) % this means action_potential only contains 1 cycle of time
    % give an generic value
    CL = 180;
end

% window of interest
woi = [t_start, t_start + round(CL*cycle_length_percentage)]; 
lat = nan(size(action_potential,2),1);
for n = 1:size(action_potential,2)
    marker = markers(n,:);
    marker(isnan(marker)) = [];
    m = marker(marker>=woi(1) & marker<=woi(2));
    if ~isempty(m)
        lat(n) = m(1);
    end
end

% shift values so that it starts at 1 ms
lat = lat - min(lat) + 1;

end