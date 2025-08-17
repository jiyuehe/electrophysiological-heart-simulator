function simulation_input = assign_simulation_input(n_voxel,dt,t_final,fiber_flag,fiber_orientation,pacing_voxel_id,c_voxel,pacing_start_time,pacing_cycle_length,rotor_method,r)

% create J_stim
if rotor_method == 0 || rotor_method == 2  || rotor_method == 3  || rotor_method == 4
    signal = J_stim_create_via_pacing_voxel(pacing_voxel_id,t_final,dt,pacing_start_time,pacing_cycle_length);
    stimulus = j_stim_compress(pacing_voxel_id,signal); % compress J_stim: only save the non-zero values
elseif rotor_method == 1 || rotor_method == 5
    % put something to be compatible with focal code
    stimulus.N = 1; 
    stimulus.voxel_id = 1;
    stimulus.time_id = 1;
    stimulus.non_zero_value = 1;
    stimulus.pacing_voxel_id = 1;
end

debug_plot = 0;
if debug_plot == 1
    [pacing_voxel_id,signal] = j_stim_decompress(stimulus);
    
    figure;
    hold on;
    plot(signal(1,:),'r');
    plot(signal(15,:),'b');
    hold off;
end

tau_in_voxel = ones(n_voxel,1) * 0.3;
tau_out_voxel = ones(n_voxel,1) * 6;
tau_open_voxel = ones(n_voxel,1) * 120;
tau_close_voxel = ones(n_voxel,1) * 80;
v_gate_voxel = ones(n_voxel,1) * 0.13;

% fiber
D0 = cell(n_voxel,1);
for n = 1:n_voxel
    if fiber_flag == 1
        e1 = fiber_orientation(n,:)';
        D0{n} = r*eye(3) + (1-r)*(e1*e1');
    elseif fiber_flag == 0
        % here r = 1
        D0{n} = eye(3);
    end
end

% simulation input
simulation_input.n_voxel = n_voxel;
simulation_input.dt = dt;
simulation_input.t_final = t_final;
simulation_input.stimulus = stimulus;
simulation_input.D0 = D0;
simulation_input.tau_in_voxel = tau_in_voxel;
simulation_input.tau_out_voxel = tau_out_voxel;
simulation_input.tau_open_voxel = tau_open_voxel;
simulation_input.tau_close_voxel = tau_close_voxel;
simulation_input.v_gate_voxel = v_gate_voxel;
simulation_input.c_voxel = c_voxel;
simulation_input.rotor_method = rotor_method;

end
