function electrogram_unipolar = compute_uni_egm_vectorized(electrode_xyz,voxel_for_each_vertex,voxel,fiber_flag,...
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

D11_b = repmat(D11,1,n_electrode);
D12_b = repmat(D12,1,n_electrode);
D13_b = repmat(D13,1,n_electrode);
D21_b = repmat(D21,1,n_electrode);
D22_b = repmat(D22,1,n_electrode);
D23_b = repmat(D23,1,n_electrode);
D31_b = repmat(D31,1,n_electrode);
D32_b = repmat(D32,1,n_electrode);
D33_b = repmat(D33,1,n_electrode);

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

% compute electrogram
T = size(action_potential,1);
electrogram_unipolar = zeros(n_electrode,T);
parfor t_id = 1:T
    if mod(t_id,T/5) == 0
        disp(['compute electrogram ',num2str(t_id/T*100),'%']);
    end

    if use_all_voxel_flag == 1
        dvdx_b = repmat(dvdx(t_id,:)',1,n_electrode);
        dvdy_b = repmat(dvdy(t_id,:)',1,n_electrode);
        dvdz_b = repmat(dvdz(t_id,:)',1,n_electrode);
    elseif use_all_voxel_flag == 0
        dvdx_b = repmat(dvdx(t_id,voxel_for_each_vertex)',1,n_electrode);
        dvdy_b = repmat(dvdy(t_id,voxel_for_each_vertex)',1,n_electrode);
        dvdz_b = repmat(dvdz(t_id,voxel_for_each_vertex)',1,n_electrode);
    end

    electrogram_unipolar(:,t_id) = sum( d ./ (l.^3) .* ... 
        ( (D11_b.*dvdx_b + D12_b.*dvdy_b + D13_b.*dvdz_b).*l_x + ...
        (D21_b.*dvdx_b + D22_b.*dvdy_b + D23_b.*dvdz_b).*l_y + ...
        (D31_b.*dvdx_b + D32_b.*dvdy_b + D33_b.*dvdz_b).*l_z ) );
end

% adjust the electrogram magnitude
typical_magnitude = 1; % unit: volt

voltage = zeros(n_electrode,1);
for n = 1:n_electrode
    voltage(n) = range(electrogram_unipolar(n,:));
end
voltage_mean = mean(voltage);

scale = typical_magnitude / voltage_mean;
electrogram_unipolar = electrogram_unipolar * scale;

debug_plot = 0;
if debug_plot == 1
    t = 1:size(electrogram_unipolar,2);
    figure;
    plot(t,electrogram_unipolar);
end

end
