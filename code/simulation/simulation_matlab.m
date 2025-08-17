function [sim_u_voxel,sim_h_voxel] = simulation_matlab(simulation_input,geometry)

n_voxel = simulation_input.n_voxel; % total number of nodes

neighbor_id_2d = geometry.volume.voxel_based_voxels;
P_2d = equation_parts(n_voxel,neighbor_id_2d,simulation_input);

neighbor_id_2d_2 = neighbor_id_2d; % Matlab cannot index 0, so 0 is replaced with 1, 
neighbor_id_2d_2(neighbor_id_2d_2==0) = 1; % and this does not matter because the associated delta_2d will be 0

% resonstruct J_stim signal
stimulus = simulation_input.stimulus;
[pacing_voxel_id,signal] = j_stim_decompress(stimulus);

% compute simulation
u_current = zeros(n_voxel,1); % initial value 0, set all voxel at rest
h_current = ones(n_voxel,1); % initial value 1, set all voxel at rest
u_next = zeros(n_voxel,1);
h_next = zeros(n_voxel,1);
t_final = simulation_input.t_final; % unit: ms
dt = simulation_input.dt;
T = length(dt:dt:t_final); % number of simulation time steps
Delta = geometry.volume.delta; % voxel spatial distance
id_save = 1;
sim_u_voxel = zeros(n_voxel,t_final); % sampling frequency at 1 kHz
sim_h_voxel = zeros(n_voxel,t_final); % sampling frequency at 1 kHz
for t = 1:T
    do_flag = 1;    
    if do_flag == 1 && mod(t,T/5) == 1
        disp(['simulation ',num2str((t-1)/T*100),'%']);
    end

    % compute diffusion term
    diffusion_term = P_2d(:,21) / (4*Delta^2) .* ...
        (P_2d(:,1) .* (u_current(neighbor_id_2d_2(:,1))-u_current) + P_2d(:,2) .* (u_current(neighbor_id_2d_2(:,2))-u_current) + ...
        P_2d(:,3) .* (u_current(neighbor_id_2d_2(:,3))-u_current) + P_2d(:,4) .* (u_current(neighbor_id_2d_2(:,4))-u_current) + ...
        P_2d(:,5) .* (u_current(neighbor_id_2d_2(:,5))-u_current) + P_2d(:,6) .* (u_current(neighbor_id_2d_2(:,6))-u_current) + ...
        P_2d(:,7) .* (u_current(neighbor_id_2d_2(:,1))-u_current(neighbor_id_2d_2(:,2))) + ...
        P_2d(:,8) .* (u_current(neighbor_id_2d_2(:,3))-u_current(neighbor_id_2d_2(:,4))) + ...
        P_2d(:,9) .* (u_current(neighbor_id_2d_2(:,5))-u_current(neighbor_id_2d_2(:,6))) + ...
        P_2d(:,10) .* (u_current(neighbor_id_2d_2(:,7))-u_current(neighbor_id_2d_2(:,9))) + ...
        P_2d(:,11) .* (u_current(neighbor_id_2d_2(:,10))-u_current(neighbor_id_2d_2(:,8))) + ...
        P_2d(:,12) .* (u_current(neighbor_id_2d_2(:,15))-u_current(neighbor_id_2d_2(:,17))) + ...
        P_2d(:,13) .* (u_current(neighbor_id_2d_2(:,18))-u_current(neighbor_id_2d_2(:,16))) + ...
        P_2d(:,14) .* (u_current(neighbor_id_2d_2(:,11))-u_current(neighbor_id_2d_2(:,13))) + ...
        P_2d(:,15) .* (u_current(neighbor_id_2d_2(:,14))-u_current(neighbor_id_2d_2(:,12))) );
    
    % compute the next time step value of u
    J_stim = zeros(n_voxel,1);
    J_stim(pacing_voxel_id) = signal(:,t);
    u_next = ((h_current .* u_current.^2 .* (1-u_current) ./ P_2d(:,18)) + ...
        (-u_current ./ P_2d(:,19)) + J_stim + diffusion_term) * dt + u_current;
    
    % compute the next time step value of h
    h_next_1 = ((1-h_current) ./ P_2d(:,16)) * dt + h_current;
    h_next_2 = (-h_current ./ P_2d(:,17)) * dt + h_current;
    id_1 = u_current < P_2d(:,20);
    id_2 = u_current >= P_2d(:,20);
    h_next(id_1) = h_next_1(id_1);
    h_next(id_2) = h_next_2(id_2);

    % update value
    u_current = u_next;
    h_current = h_next;

    % save value at 1 kHz
    if mod(t,1/dt) == 1
        sim_u_voxel(:,id_save) = u_current;
        sim_h_voxel(:,id_save) = h_current;
        id_save = id_save + 1;
    end
end

sim_u_voxel = sim_u_voxel';
sim_h_voxel = sim_h_voxel';

end
