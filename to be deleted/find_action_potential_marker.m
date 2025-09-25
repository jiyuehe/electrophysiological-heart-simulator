function m = find_action_potential_marker(s,threshold)

% m is the activation marker. diff(m) will be the cycle length

above_threshold_id = find(s>threshold);

debug_plot = 0;
if debug_plot == 1
    figure;
    plot(s,'k');
    hold on;
    plot(s,'.b');
    plot(above_threshold_id,s(above_threshold_id),'.r');
    hold off;
end

if sum(above_threshold_id) ~= 0 % there is at least one activation
    b = find(diff(above_threshold_id)>1); % beat separator 
    
    debug_plot = 0;
    if debug_plot == 1
        figure;
        plot(s,'k');
        hold on;
        plot(s,'.b');
        plot(above_threshold_id,s(above_threshold_id),'.r');
        scatter(above_threshold_id(b),s(above_threshold_id(b)),'g');
        hold off;
    end
    
    segment_marker = [above_threshold_id(1); above_threshold_id(b); above_threshold_id(b+1); above_threshold_id(end)];
    segment_marker = sort(segment_marker,'ascend');

    segment_marker = unique(segment_marker);
    if rem(length(segment_marker),2) == 1
        segment_marker(end) = [];
    end
    
    debug_plot = 0;
    if debug_plot == 1
        figure;
        plot(s,'k');
        hold on;
        plot(s,'.b');
        plot(above_threshold_id,s(above_threshold_id),'.r');
        scatter(segment_marker,s(segment_marker),'g');
        hold off;
    end
    
    m = zeros(1,length(segment_marker)/2);
    value = zeros(size(m));
    for i = 1:length(m)
        [value(i),id] = max( diff(s(segment_marker((i-1)*2+1):segment_marker(i*2))) ); % max ds/dt
        m(i) = id + segment_marker((i-1)*2+1) - 1;
    end
    
    id_to_delete = value <= 0;
    m(id_to_delete) = []; % sometimes s starts with a portion of the action potential, need to delete marker on such case, only mark at positive ds/dt 

    debug_plot = 0;
    if debug_plot == 1
        figure;
        plot(s,'k');
        hold on;
        plot(s,'.b');
        scatter(m,s(m),'r');
        hold off;
    end
elseif sum(above_threshold_id) == 0 % there is no activation
    m = NaN;
end

end
