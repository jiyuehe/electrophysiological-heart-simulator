function voxel_flag = get_voxel_flag_multi_pacing(directory,file_name_vertex_flag,geometry,az,el)

% transfer vertex_flag to voxel_flag
cd(directory.data_dir); % it is necessary to re-load the variable, or Matlab could still using the old one before vertex_flag edited
load(file_name_vertex_flag); % the variable loaded is named vertex_flag
cd(directory.home_dir);

% find out the different pacing sites
pacing_sites = find(vertex_flag==2);
voxel_flag = cell(1,length(pacing_sites));
n_voxel = size(geometry.volume.voxel,1);
for pacing_id = 1:length(pacing_sites)
    pacing_vertex_id = pacing_sites(pacing_id);
    pacing_voxel_id = geometry.voxel_for_each_vertex(pacing_vertex_id);

    voxel_flag{pacing_id} = zeros(n_voxel,1);
    voxel_flag{pacing_id}(pacing_voxel_id) = 2;
end

debug_plot = 0;
if debug_plot == 1
    pacing_id = 1;
    
    % check voxel flags
    id = find(voxel_flag{pacing_id}==2);
    neighbor_id = geometry.volume.voxel_based_voxels(id,:);
    neighbor_id(neighbor_id==0) = [];
    pacing_voxel_id = [id(:); neighbor_id(:)];
    pacing_voxel_id = unique(pacing_voxel_id);

    figure;
    voxel = geometry.volume.voxel;
    scatter3(voxel(:,1),voxel(:,2),voxel(:,3),3,'k','filled');
    hold on;
    scatter3(voxel(pacing_voxel_id,1),voxel(pacing_voxel_id,2),voxel(pacing_voxel_id,3),10,'r','filled');
    hold off;
    axis tight equal vis3d;
    view(az,el); 
    rotate3d on;
    xlabel('x');
    ylabel('y');
end

end