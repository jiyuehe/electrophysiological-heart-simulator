% NOTE: 
% To run this code, go to the Editor tab at the top of the MATLAB UI and click Run. 
% If a pop-up window appears, click Change Folder.

clear; close all; clc;
directory.home_dir = pwd; % pwd returns the path to the current folder
directory.data_dir = fullfile(directory.home_dir,'data'); % generate the file path based on your operating system: use '\' for Windows, '/' for Linux, and so on.
directory.result_dir = fullfile(directory.home_dir,'result');
addpath(genpath('code')); % add the folder to the MATLAB path so that it can recognize and access the code files inside

file_name = 'heart_example';
cd(directory.data_dir);
load([file_name,'.mat']);
cd(directory.home_dir);

az = -170;
el = -90;

arrhythmia_type = 'focal'; % 'focal', 'rotor'
if strcmp(arrhythmia_type,'focal')
    rotor_method = 0; % 0 is focal source
    file_name_vertex_flag = [file_name,'_',arrhythmia_type,'_flag'];
elseif strcmp(arrhythmia_type,'rotor')
    rotor_method = 2; % 2 is rotor via S1-S2 pacing method
    file_name_vertex_flag = [file_name,'_',arrhythmia_type,'_',num2str(rotor_method),'_flag'];
end

% simulation rhythm setup, assign vertex flag
do_flag = 0;
if do_flag == 1
    cd(directory.data_dir);
    if ~isfile([file_name_vertex_flag,'.mat'])
        vertex_flag = zeros(size(data.geometry.edited.vertex,1),1);
        save(file_name_vertex_flag,'vertex_flag');
    end
    cd(directory.home_dir);

    assign_vertex_flag(directory,file_name_vertex_flag,data.geometry.edited,az,el);
end

dt = 0.05; % unit: millisecond. if dt is not small enough, simulation will give NaN. Generally, if c <= 1.0, can use dt = 0.05
t_final = 500;
c = 0.5; % diffusion coefficient. since the 2D plane is not big, c = 1 will be a little too large
n_voxel = size(data.geometry.edited.volume.voxel,1);
fiber_flag = 0; % no fiber
r = []; % no fiber
fiber_orientation = []; % no fiber

voxel_flag_multi_pacing = get_voxel_flag_multi_pacing(directory,file_name_vertex_flag,data.geometry.edited,az,el);

c_voxel = c * ones(n_voxel,1);

pacing_id = 1;

if strcmp(arrhythmia_type,'focal')
    voxel_flag = voxel_flag_multi_pacing{pacing_id};

    id = find(voxel_flag==2);
    neighbor_id = data.geometry.edited.volume.voxel_based_voxels(id,:);
    neighbor_id(neighbor_id==0) = [];
    pacing_voxel_id = [id(:); neighbor_id(:)];
    pacing_voxel_id = unique(pacing_voxel_id);
    
    pacing_start_time = zeros(length(pacing_voxel_id),1);
    pacing_start_time(:) = 1; % unit: millisecond
    pacing_cycle_length = zeros(length(pacing_voxel_id),1);
    pacing_cycle_length(:) = 180; % unit: millisecond
elseif strcmp(arrhythmia_type,'rotor')
    if rotor_method == 2 % s1-s2 pacing method. 
        t_final = 2500;
        cl_1 = 2*t_final;
        cl_2 = 2*t_final;
        s2_time = 230; % unit: millisecond
        [pacing_voxel_id,pacing_start_time,pacing_cycle_length] = s1s2_pacing_setting(geometry,voxel_flag,s2_time,cl_1,cl_2);
    end
end

do_flag = 1;
if do_flag == 1
    simulation_input = assign_simulation_input(n_voxel,dt,t_final,fiber_flag,fiber_orientation,pacing_voxel_id,c_voxel,pacing_start_time,...
        pacing_cycle_length,rotor_method,r);
    [action_potential,h_voxel] = simulation_matlab(simulation_input,data.geometry.edited);

    % % save
    % cd(directory.result_dir);
    % save(['action potential, ',arrhythmia_type,', ',num2str(pacing_id)],'action_potential');
    % save(['h, ',arrhythmia_type,', ',num2str(pacing_id)],'h_voxel');
    % cd(directory.home_dir);
elseif do_flag == 0
    cd(directory.result_dir);
    load(['action potential, ',arrhythmia_type,', ',num2str(pacing_id)]); % the variable loaded is named action_potential
    cd(directory.home_dir);
end

debug_plot = 0;
if debug_plot == 1
    % action potential of some voxel
    id = 1;
    figure;
    hold on;
    for i = 1:length(id)
        plot(action_potential(:,id(i)));
    end
    hold off;
end

debug_plot = 0;
if debug_plot == 1
    % movie
    action_potential2 = action_potential(:,data.geometry.edited.voxel_for_each_vertex);
    v_gate = 0.13;
    action_potential_phase = zeros(size(action_potential2));
    for id = 1:size(action_potential2,2) % convert action potential to phase
        [action_potential_phase(:,id),~] = movie_create_phase(action_potential2(:,id),v_gate);
    end
    movie_data = [];
    movie_data.vertex = data.geometry.edited.vertex;
    movie_data.face = data.geometry.edited.face;
    movie_data.signal = action_potential_phase;
    movie_data.az = az;
    movie_data.el = el;
    file_name = [arrhythmia_type, ' ', num2str(pacing_id)];
    map_data_min = 0;
    map_data_max = 1;
    frame_id = 1:10:t_final;
    movie_type = 'triangular_mesh';
    color_type = 'vertex';
    movie_play(directory,file_name,movie_data,map_data_min,map_data_max,frame_id,movie_type,color_type);
end

debug_plot = 1;
if debug_plot == 1
    % local activation time map
    cd(directory.data_dir);
    load(file_name_vertex_flag); % the variable loaded is named vertex_flag
    cd(directory.home_dir);
    
    if strcmp(arrhythmia_type,'focal')
        t_start = 1;%min(pacing_start_time);
        cycle_length_percentage = 0.99;
    elseif strcmp(arrhythmia_type,'rotor')
        t_start = 500;
        cycle_length_percentage = 0.8;
    end
    action_potential_vertex = action_potential(:,data.geometry.edited.voxel_for_each_vertex);
    lat = calculate_local_activation_time_action_potential(t_start,action_potential_vertex,cycle_length_percentage);

    vertex = data.geometry.edited.vertex;
    face = data.geometry.edited.face;
    map_data = lat;
    map_data_min = min(map_data);
    map_data_max = max(map_data);
    data_threshold = min(map_data)-0.01;
    color_wheel_type = 'red_blue';
    black_back_flag = 0;
    figure;
    set(gcf,'Position',get(0,'Screensize'));
    plot_color_map(face,vertex,map_data,map_data_min,map_data_max,data_threshold,color_wheel_type,black_back_flag,0,az,el);
    
    cd(directory.result_dir);
    exportgraphics(gca,[arrhythmia_type,', t_start = ',num2str(t_start),', ',num2str(pacing_id),'.png'],'Resolution',300);
    cd(directory.home_dir);
    close;
end
