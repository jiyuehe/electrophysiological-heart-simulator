function stimulus = j_stim_compress(pacing_voxel_id,signal)

N = length(pacing_voxel_id);
T = size(signal,2);

voxel_id = [];
time_id = [];
non_zero_value = [];

k = 1;
for n = 1:N
    v_id = pacing_voxel_id(n);

    for t_id = 1:T
        if signal(n,t_id) ~= 0
            voxel_id(k) = v_id;
            time_id(k) = t_id;
            non_zero_value(k) = signal(n,t_id);
            
            k = k + 1;
        end
    end
end

stimulus = [];
stimulus.N = N;
stimulus.T = T;
stimulus.pacing_voxel_id = pacing_voxel_id;
stimulus.voxel_id = voxel_id;
stimulus.time_id = time_id;
stimulus.non_zero_value = non_zero_value;

end