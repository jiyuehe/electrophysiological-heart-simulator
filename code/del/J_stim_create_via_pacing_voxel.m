function signal = J_stim_create_via_pacing_voxel(pacing_voxel_id,t_final,dt,pacing_start_time,pacing_cycle_length)

% pacing duration
n = 10; % ms
pacing_duration = n/dt; % make sure it is n ms no matter what dt is

% pacing starting times
N = length(pacing_voxel_id);

pacing_starts = cell(N,1);
pacing_ends = cell(N,1);
for n = 1:N
    pacing_starts{n} = pacing_start_time(n)/dt:pacing_cycle_length(n)/dt:t_final/dt-pacing_duration;
    for i = 1:length(pacing_starts{n})
        pacing_ends{n}(i) = pacing_starts{n}(i) + pacing_duration - 1;
    end
end

% pacing strength
J_stim_value = 20; % 20 is good. 10 is not large enough if space step is 0.1 mm and time step is 0.001

% pacing signal
n_step = length(dt:dt:t_final);
signal = zeros(N,n_step);
for n = 1:N
    for i = 1:length(pacing_starts{n})
        signal(n,pacing_starts{n}(i):pacing_ends{n}(i)) = J_stim_value;
    end
end

debug_plot = 0;
if debug_plot == 1
    figure;
    hold on;
    plot(signal(1,:),'r');
    plot(signal(15,:),'b');
    hold off;
end

end