function electrogram_unipolar = compute_uni_egm_gpu_cascade_electrode(electrode_xyz,voxel_for_each_vertex,voxel,fiber_flag,...
    fiber_orientation,r,use_all_voxel_flag,c_voxel,action_potential,delta,voxel_based_voxels)

% equation parts
if use_all_voxel_flag == 1
    node_xyz = voxel;
elseif use_all_voxel_flag == 0
    node_xyz = voxel(voxel_for_each_vertex,:);
end
n_node = size(node_xyz,1);

n_voxel = size(voxel,1);
D0_all = cell(n_voxel,1);
for n = 1:n_voxel
    if fiber_flag == 1
        e1 = fiber_orientation(n,:)';
        D0_all{n} = r*eye(3) + (1-r)*(e1*e1');
    elseif fiber_flag == 0
        D0_all{n} = eye(3); % here r = 1
    end
end
if use_all_voxel_flag == 1
    d = c_voxel;
    D0 = D0_all;
elseif use_all_voxel_flag == 0
    d = c_voxel(voxel_for_each_vertex);
    D0 = D0_all(voxel_for_each_vertex);
end
D11 = zeros(n_node,1);
D12 = zeros(n_node,1);
D13 = zeros(n_node,1);
D21 = zeros(n_node,1);
D22 = zeros(n_node,1);
D23 = zeros(n_node,1);
D31 = zeros(n_node,1);
D32 = zeros(n_node,1);
D33 = zeros(n_node,1);
for n = 1:n_node
    D11(n) = D0{n}(1,1);
    D12(n) = D0{n}(1,2);
    D13(n) = D0{n}(1,3);
    D21(n) = D0{n}(2,1);
    D22(n) = D0{n}(2,2);
    D23(n) = D0{n}(2,3);
    D31(n) = D0{n}(3,1);
    D32(n) = D0{n}(3,2);
    D33(n) = D0{n}(3,3);
end

% compute action potential gradient of each time frame
px = voxel_based_voxels(:,1); px_0_id = px==0; px(px_0_id) = 1;
mx = voxel_based_voxels(:,2); mx_0_id = mx==0; mx(mx_0_id) = 1;
py = voxel_based_voxels(:,3); py_0_id = py==0; py(py_0_id) = 1;
my = voxel_based_voxels(:,4); my_0_id = my==0; my(my_0_id) = 1;
pz = voxel_based_voxels(:,5); pz_0_id = pz==0; pz(pz_0_id) = 1;
mz = voxel_based_voxels(:,6); mz_0_id = mz==0; mz(mz_0_id) = 1;
dvdx = (action_potential(:,px)-action_potential(:,mx)) / (2*delta);
dvdx(:,px_0_id) = 0;
dvdx(:,mx_0_id) = 0;
dvdy = (action_potential(:,py)-action_potential(:,my)) / (2*delta);
dvdy(:,py_0_id) = 0;
dvdy(:,my_0_id) = 0;
dvdz = (action_potential(:,pz)-action_potential(:,mz)) / (2*delta);
dvdz(:,pz_0_id) = 0;
dvdz(:,mz_0_id) = 0;

if use_all_voxel_flag == 1
    dvdx_b = dvdx';
    dvdy_b = dvdy';
    dvdz_b = dvdz';
elseif use_all_voxel_flag == 0
    dvdx_b = dvdx(:,voxel_for_each_vertex)';
    dvdy_b = dvdy(:,voxel_for_each_vertex)';
    dvdz_b = dvdz(:,voxel_for_each_vertex)';
end

% distance from electrode to nodes
n_electrode = size(electrode_xyz,1);
l = zeros(n_node,n_electrode);
l_x = zeros(n_node,n_electrode);
l_y = zeros(n_node,n_electrode);
l_z = zeros(n_node,n_electrode);
for e_id = 1:n_electrode
    l_x(:,e_id) = node_xyz(:,1) - electrode_xyz(e_id,1);
    l_y(:,e_id) = node_xyz(:,2) - electrode_xyz(e_id,2);
    l_z(:,e_id) = node_xyz(:,3) - electrode_xyz(e_id,3);
    l(:,e_id) = sqrt(l_x(:,e_id).^2 + l_y(:,e_id).^2 + l_z(:,e_id).^2);
end
l(l<1) = 1; % electrode is at least 1 mm away from tissue

% cascade all electrodes into 1D array
n_parts = 14;
n_node2 = n_node*n_electrode;
egm_part2 = zeros(n_node2,1);
P2 = zeros(n_node2*n_parts,1);
for e_id = 1:n_electrode
    node_id = (e_id-1)*n_node+1 : (e_id-1)*n_node+n_node;
    
    P2(((node_id)-1) * n_parts + 1) = d;
    P2(((node_id)-1) * n_parts + 2) = D11;
    P2(((node_id)-1) * n_parts + 3) = D12;
    P2(((node_id)-1) * n_parts + 4) = D13;
    P2(((node_id)-1) * n_parts + 5) = D21;
    P2(((node_id)-1) * n_parts + 6) = D22;
    P2(((node_id)-1) * n_parts + 7) = D23;
    P2(((node_id)-1) * n_parts + 8) = D31;
    P2(((node_id)-1) * n_parts + 9) = D32;
    P2(((node_id)-1) * n_parts + 10) = D33;
    P2(((node_id)-1) * n_parts + 11) = l(:,e_id);
    P2(((node_id)-1) * n_parts + 12) = l_x(:,e_id);
    P2(((node_id)-1) * n_parts + 13) = l_y(:,e_id);
    P2(((node_id)-1) * n_parts + 14) = l_z(:,e_id);
end

% load kernel
threads_per_block = 32; % threads_per_block is an integer range from 1 to 1024, it is found that set it to 32 runs quite fast
cudaFilename = 'kernel_compute_unipolar_electrogram.cu';
ptxFilename = 'kernel_compute_unipolar_electrogram.ptx'; % this is the compiled cuda file
kernel_egm = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename );
kernel_egm.ThreadBlockSize = [threads_per_block, 1, 1]; % threads per block
kernel_egm.GridSize = [ceil(n_node2/threads_per_block), 1]; % blocks per grid

% put variables into gpu
egm_part2 = gpuArray(double(egm_part2));
P2 = gpuArray(double(P2));

% compute electrogram
T = size(action_potential,1);
n_electrode = size(electrode_xyz,1);
electrogram_unipolar = zeros(n_electrode,T);
egm_part2_all_time = zeros(n_node2,T);
n_node2 = gpuArray(double(n_node2));
n_parts = gpuArray(double(n_parts));
for t_id = 1:T
    if mod(t_id,T/5) == 0
        disp(['compute electrogram ',num2str(t_id/T*100),'%']);
    end
    
    dvdx = repmat(dvdx_b(:,t_id),n_electrode,1);
    dvdy = repmat(dvdy_b(:,t_id),n_electrode,1);
    dvdz = repmat(dvdz_b(:,t_id),n_electrode,1);
    dvdx = gpuArray(double(dvdx));
    dvdy = gpuArray(double(dvdy));
    dvdz = gpuArray(double(dvdz));
    egm_part2_all_time(:,t_id) = feval(kernel_egm, egm_part2, P2, n_node2, n_parts, dvdx, dvdy, dvdz);
end

for e_id = 1:n_electrode
    node_id = (e_id-1)*n_node+1 : (e_id-1)*n_node+n_node;
    electrogram_unipolar(e_id,:) = sum(egm_part2_all_time(node_id,:));
end

%{
% NOTE: this code may be used to compute only a segment of electrodes' electrograms
%       it is better to scale the egm after all electrodes' egm are computed
% adjust the electrogram magnitude
typical_magnitude = 1; % unit: volt

voltage = zeros(n_electrode,1);
for n = 1:n_electrode
    voltage(n) = range(electrogram_unipolar(n,:));
end
voltage_mean = mean(voltage);

scale = typical_magnitude / voltage_mean;
electrogram_unipolar = electrogram_unipolar * scale;
%}

debug_plot = 0;
if debug_plot == 1
    t = 1:size(electrogram_unipolar,2);
    e_id = [1 2 3];
    figure;
    plot(t,electrogram_unipolar(e_id,:));
end

end
