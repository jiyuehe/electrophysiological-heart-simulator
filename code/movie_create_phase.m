function [action_potential_phase,activation_phase] = movie_create_phase(action_potential,v_gate)

% action potential phase: 0-2pi of the action potential duration
% ----------------------------------------------------------------------------------------------------
debug_plot = 0;
if debug_plot == 1
    figure;
    plot(action_potential,'b');
end

a = find( action_potential > v_gate ); % time index

debug_plot = 0;
if debug_plot == 1
    figure;
    plot(action_potential,'b');
    hold on;
    plot(a,action_potential(a),'r.');
    hold off;
end

if ~isempty(a)
    b = find( abs(diff(a)) > 1);
    p = [a(1); a(b); a(b+1); a(end)];
    p = sort(p,'ascend');

    debug_plot = 0;
    if debug_plot == 1
        figure;
        plot(action_potential,'b');
        hold on;
        scatter(p,action_potential(p),100,'.r');
        hold off;
    end
    
    % phase interval
    phase_interval = zeros(length(p)/2,2);
    for n = 1:length(p)/2
        phase_interval(n,:) = p([(n-1)*2+1 (n-1)*2+2]);
    end

    % phase
    L_median = ceil( median(phase_interval(:,2)-phase_interval(:,1)) );
    action_potential_phase = zeros(size(action_potential));
    for n = 1:size(phase_interval,1)
        m = phase_interval(n,:);
        L = length(m(1):m(2)-1);

        if m(1) == 1
            action_potential_phase(m(1):m(2)-1) = linspace(m(1)/L_median,1,L); % linspace(X1, X2, N) generates N points between X1 and X2.
        end

        if m(1) ~= 1 && m(2) ~= length(action_potential_phase)        
            action_potential_phase(m(1):m(2)-1) = linspace(0,1,L);
        end

        if m(2) == length(action_potential_phase)
            action_potential_phase(m(1):m(2)-1) = linspace(0,min(1,L/L_median),L);
        end
    end
elseif isempty(a)
    action_potential_phase = zeros(size(action_potential));
end

debug_plot = 0;
if debug_plot == 1
    figure;
    plot(action_potential,'b');
    hold on;
    plot(action_potential_phase,'r');
    hold off;
    axis tight;
end

% activation phase: 0-2pi in between 2 activation peaks
% ----------------------------------------------------------------------------------------------------
[~,p] = findpeaks(action_potential,'MinPeakHeight',mean(action_potential)); % time index

debug_plot = 0;
if debug_plot == 1
    figure;
    plot(action_potential,'b');
    hold on;
    scatter(p,action_potential(p),'r');
    hold off;
end

if ~isempty(p)
    % include t start and t final
    if p(1) ~= 1
        p = [1; p];
    end
    if p(end) ~= length(action_potential)
        p = [p; length(action_potential)];
    end
    
    % phase interval
    phase_interval = zeros(length(p)-1,2);
    for n = 1:length(p)-1
        phase_interval(n,:) = p([n n+1]);
    end

    % phase
    L_median = ceil( median(phase_interval(:,2)-phase_interval(:,1)) );
    activation_phase = zeros(size(action_potential));
    for n = 1:size(phase_interval,1)
        m = phase_interval(n,:);
        L = length(m(1):m(2)-1);

        if m(1) == 1
            activation_phase(m(1):m(2)-1) = linspace(L/L_median,1,L); % linspace(X1, X2, N) generates N points between X1 and X2.
            % NOTE: this is different that creating action_potential_phase
        end

        if m(1) ~= 1 && m(2) ~= length(activation_phase)        
            activation_phase(m(1):m(2)-1) = linspace(0,1,L);
        end

        if m(2) == length(activation_phase)
            activation_phase(m(1):m(2)) = linspace(0,min(1,L/L_median),L+1);
        end
    end
elseif isempty(p)
    activation_phase = zeros(size(action_potential));
end

debug_plot = 0;
if debug_plot == 1
    figure;
    plot(action_potential,'b');
    hold on;
    plot(activation_phase,'r');
    hold off;
    axis tight;
end

end