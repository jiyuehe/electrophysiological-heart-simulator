% local activation time map

t_start = 1; % min(pacing_start_time);
cycle_length_percentage = 0.99;
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
