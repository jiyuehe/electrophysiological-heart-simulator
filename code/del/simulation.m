function [sim_u_voxel,sim_h_voxel] = simulation(simulation_input,geometry,voxel_flag)

n_voxel = simulation_input.n_voxel;
rotor_method = simulation_input.rotor_method;

if rotor_method == 0 % focal
    u_current = zeros(n_voxel,1); % initial value 0, set all voxel at rest
    h_current = ones(n_voxel,1); % initial value 1, set all voxel at rest
elseif rotor_method == 1 % rotor, 2-patch method
    voxel_id_21 = find(voxel_flag==21);
    voxel_id_22 = find(voxel_flag==22);
    voxel_id_23 = find(voxel_flag==23);

    debug_plot = 0;
    if debug_plot == 1
        voxel = geometry.volume.voxel;
        
        figure;
        hold on;
        scatter3(voxel(:,1),voxel(:,2),voxel(:,3),2,'k','filled');
        scatter3(voxel(voxel_id_21,1),voxel(voxel_id_21,2),voxel(voxel_id_21,3),2,'r','filled');
        scatter3(voxel(voxel_id_22,1),voxel(voxel_id_22,2),voxel(voxel_id_22,3),2,'b','filled');
        scatter3(voxel(voxel_id_23,1),voxel(voxel_id_23,2),voxel(voxel_id_23,3),2,'g','filled');
        hold off;
        view(3);
        axis equal tight vis3d;
        set(gcf,'color','w');
        rotate3d on;
        xlabel('x'); ylabel('y'); zlabel('z');
    end

    u_current = zeros(n_voxel,1);
    h_current = zeros(n_voxel,1); % set to zero, not set to at rest
    u_current(voxel_id_21) = 1;
    h_current(voxel_id_22) = 1;
    u_current(voxel_id_23) = 1;
    h_current(voxel_id_23) = 1;
elseif rotor_method == 2 % rotor, S1-S2 pacing method
    u_current = zeros(n_voxel,1); % initial value 0, set all voxel at rest
    h_current = ones(n_voxel,1); % initial value 1, set all voxel at rest
elseif rotor_method == 3 % rotor, reset region method
    u_current = zeros(n_voxel,1); % initial value 0, set all voxel at rest
    h_current = ones(n_voxel,1); % initial value 1, set all voxel at rest
    
    reset_time_start = 30 / simulation_input.dt; % unit: millisecond/dt
    reset_time_end = 180 / simulation_input.dt; % unit: millisecond/dt
elseif rotor_method == 4 % rotor, reset half mesh method
    u_current = zeros(n_voxel,1); % initial value 0, set all voxel at rest
    h_current = ones(n_voxel,1); % initial value 1, set all voxel at rest
    
    reset_time_start = 150 / simulation_input.dt; % unit: millisecond/dt
    reset_time_end = 200 / simulation_input.dt; % unit: millisecond/dt
    
    voxel = geometry.volume.voxel;
    id = find(voxel_flag==11);
    id = id(1);
    x_value = voxel(id,1);
    reset_voxel_id = find(voxel(:,1) < x_value/2);
    
    debug_plot = 0;
    if debug_plot == 1
        figure;
        hold on;
        scatter3(voxel(:,1),voxel(:,2),voxel(:,3),2,'k','filled');
        scatter3(voxel(reset_voxel_id,1),voxel(reset_voxel_id,2),voxel(reset_voxel_id,3),2,'r','filled');
        hold off;
        view(3);
        axis equal tight vis3d;
        set(gcf,'color','w');
        rotate3d on;
        xlabel('x'); ylabel('y'); zlabel('z');
    end
elseif rotor_method == 5 % rotor, given initial rotating u and h values
    u_current = simulation_input.u_initial;
    h_current = simulation_input.h_initial;
end

neighbor_id_2d = geometry.volume.voxel_based_voxels;
P_2d = equation_parts(n_voxel,neighbor_id_2d,simulation_input);

% convert 2D variable into 1D for CUDA GPU code
P = zeros(n_voxel*21,1);
for v_id = 1:n_voxel
    P((v_id-1) * 21 + 1) = P_2d(v_id,1);
    P((v_id-1) * 21 + 2) = P_2d(v_id,2);
    P((v_id-1) * 21 + 3) = P_2d(v_id,3);
    P((v_id-1) * 21 + 4) = P_2d(v_id,4);
    P((v_id-1) * 21 + 5) = P_2d(v_id,5);
    P((v_id-1) * 21 + 6) = P_2d(v_id,6);
    P((v_id-1) * 21 + 7) = P_2d(v_id,7);
    P((v_id-1) * 21 + 8) = P_2d(v_id,8);
    P((v_id-1) * 21 + 9) = P_2d(v_id,9);
    P((v_id-1) * 21 + 10) = P_2d(v_id,10);
    P((v_id-1) * 21 + 11) = P_2d(v_id,11);
    P((v_id-1) * 21 + 12) = P_2d(v_id,12);
    P((v_id-1) * 21 + 13) = P_2d(v_id,13);
    P((v_id-1) * 21 + 14) = P_2d(v_id,14);
    P((v_id-1) * 21 + 15) = P_2d(v_id,15);
    P((v_id-1) * 21 + 16) = P_2d(v_id,16);
    P((v_id-1) * 21 + 17) = P_2d(v_id,17);
    P((v_id-1) * 21 + 18) = P_2d(v_id,18);
    P((v_id-1) * 21 + 19) = P_2d(v_id,19);
    P((v_id-1) * 21 + 20) = P_2d(v_id,20);
    P((v_id-1) * 21 + 21) = P_2d(v_id,21);
end

neighbor_id_2d = neighbor_id_2d - 1; % subtract 1: convert Matlab index to CUDA index (CUDA code index starts with 0)
neighbor_id_2d(neighbor_id_2d == -1) = 0; %  This will not be a problem, because these term will be multiply by 0 (the associated P value will be 0)
neighbor_id = zeros(n_voxel*18,1); 
for v_id = 1:n_voxel
    neighbor_id((v_id-1) * 18 + 1) = neighbor_id_2d(v_id,1);
    neighbor_id((v_id-1) * 18 + 2) = neighbor_id_2d(v_id,2);
    neighbor_id((v_id-1) * 18 + 3) = neighbor_id_2d(v_id,3);
    neighbor_id((v_id-1) * 18 + 4) = neighbor_id_2d(v_id,4);
    neighbor_id((v_id-1) * 18 + 5) = neighbor_id_2d(v_id,5);
    neighbor_id((v_id-1) * 18 + 6) = neighbor_id_2d(v_id,6);
    neighbor_id((v_id-1) * 18 + 7) = neighbor_id_2d(v_id,7);
    neighbor_id((v_id-1) * 18 + 8) = neighbor_id_2d(v_id,8);
    neighbor_id((v_id-1) * 18 + 9) = neighbor_id_2d(v_id,9);
    neighbor_id((v_id-1) * 18 + 10) = neighbor_id_2d(v_id,10);
    neighbor_id((v_id-1) * 18 + 11) = neighbor_id_2d(v_id,11);
    neighbor_id((v_id-1) * 18 + 12) = neighbor_id_2d(v_id,12);
    neighbor_id((v_id-1) * 18 + 13) = neighbor_id_2d(v_id,13);
    neighbor_id((v_id-1) * 18 + 14) = neighbor_id_2d(v_id,14);
    neighbor_id((v_id-1) * 18 + 15) = neighbor_id_2d(v_id,15);
    neighbor_id((v_id-1) * 18 + 16) = neighbor_id_2d(v_id,16);
    neighbor_id((v_id-1) * 18 + 17) = neighbor_id_2d(v_id,17);
    neighbor_id((v_id-1) * 18 + 18) = neighbor_id_2d(v_id,18);
end

u_next = zeros(n_voxel,1);
h_next = zeros(n_voxel,1);
dt = simulation_input.dt;
Delta = geometry.volume.delta; % voxel spatial distance
N = simulation_input.stimulus.N;
L = length(simulation_input.stimulus.voxel_id);
voxel_id = simulation_input.stimulus.voxel_id - 1;
time_id = simulation_input.stimulus.time_id - 1;
non_zero_value = simulation_input.stimulus.non_zero_value;

% load kernel
threads_per_block = 32; % threads_per_block is an integer range from 1 to 1024, it is found that set it to 32 runs quite fast

cudaFilename = 'kernel_heart_simulation_has_stim.cu';
ptxFilename = 'kernel_heart_simulation_has_stim.ptx'; % this is the compiled cuda file
kernel_has_stim = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename );
kernel_has_stim.ThreadBlockSize = [threads_per_block, 1, 1]; % threads per block
kernel_has_stim.GridSize = [ceil(n_voxel/threads_per_block), 1]; % blocks per grid

cudaFilename = 'kernel_heart_simulation_no_stim.cu';
ptxFilename = 'kernel_heart_simulation_no_stim.ptx'; % this is the compiled cuda file
kernel_no_stim = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename );
kernel_no_stim.ThreadBlockSize = [threads_per_block, 1, 1]; % threads per block
kernel_no_stim.GridSize = [ceil(n_voxel/threads_per_block), 1]; % blocks per grid

% put variables into gpu
u_current = gpuArray(double(u_current));
u_next = gpuArray(double(u_next));
h_current = gpuArray(double(h_current));
h_next = gpuArray(double(h_next));
dt = double(dt);
Delta = double(Delta);
neighbor_id = gpuArray(int32(neighbor_id));
P = gpuArray(double(P));
N = int32(N);
L = int32(L);
voxel_id = gpuArray(int32(voxel_id));
time_id = gpuArray(int32(time_id));
non_zero_value = gpuArray(double(non_zero_value));
n_voxel = int32(n_voxel);

% simulation
t_final = simulation_input.t_final; % unit: ms
T = length(dt:dt:t_final);
id_save = 1;
sim_u_voxel = zeros(n_voxel,t_final); % sampling frequency at 1 kHz
sim_h_voxel = zeros(n_voxel,t_final);
for t = 0:T-1 % for CUDA GPU computing, index starts at 0
    do_flag = 1;
    if do_flag == 1 && mod(t,round(T/5)) == 1
        disp(['simulation ',num2str((t-1)/T*100),'%']);
    end
    
    % rotor related
    if rotor_method == 3 && t >= reset_time_start && t <= reset_time_end
        % set region 4 to at rest
        reset_id = (voxel_flag == 4);
        u_current(reset_id) = 0;
        h_current(reset_id) = 1;
    end
    if rotor_method == 4 && t >= reset_time_start && t <= reset_time_end
        % set half mesh to at rest
        u_current(reset_voxel_id) = 0;
        h_current(reset_voxel_id) = 1;
    end
    
%     % assign non-conductance via setting u = 0
%     if ~isempty(line_of_block_id) && t <= 20/dt
%         u_current(line_of_block_id) = 0;
%     end
    
    % solve differential equations
    if ~isempty(simulation_input.stimulus.pacing_voxel_id)
        id = time_id==t;
        if sum(id)
            voxel_id2 = voxel_id(id);
            L = length(voxel_id2);
            [u_next, h_next] = feval(kernel_has_stim, u_next, h_next, u_current, h_current, t, dt, Delta, neighbor_id, P, ...
                N, L, voxel_id2, time_id, non_zero_value, n_voxel);
        elseif ~sum(id)
            [u_next, h_next] = feval(kernel_no_stim, u_next, h_next, u_current, h_current, t, dt, Delta, neighbor_id, P, N, L, n_voxel);
        end
    elseif isempty(simulation_input.stimulus.pacing_voxel_id)
        [u_next, h_next] = feval(kernel_no_stim, u_next, h_next, u_current, h_current, t, dt, Delta, neighbor_id, P, N, L, n_voxel);
    end

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

debug_plot = 0;
if debug_plot == 1
    v_id = 5;
    figure;
    plot(sim_u_voxel(v_id,:));
end

sim_u_voxel = sim_u_voxel';
sim_h_voxel = sim_h_voxel';

end
