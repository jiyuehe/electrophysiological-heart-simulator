function movie_play(directory,file_name,data,data_min,data_max,frame_id,movie_type,color_type)

save_flag = 1;

signal = data.signal;
az = data.az;
el = data.el;

if length(frame_id) > 1 && save_flag == 1
    cd(directory.result_dir);
    writerObj = VideoWriter(file_name);
    writerObj.FrameRate = 10; % frames per second
    open(writerObj);
    figure('units','pixels','position',[100 100 1280 720]);
end

if length(frame_id) == 1
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
end

movegui('center');
set(gcf,'color','w');

if length(frame_id) > 1  && save_flag == 1
    pause(2); % pause a little bit to wait for the figure to fully expanded
end
set(gcf,'Renderer','zbuffer');

if strcmp(movie_type,'point_cloud')
    voxel = data.voxel;
    h = scatter3(voxel(:,1),voxel(:,2),voxel(:,3),15,'b','filled'); 
elseif strcmp(movie_type,'triangular_mesh')
    vertex = data.vertex;
    face = data.face;
    if strcmp(color_type,'face')
        h = patch('Faces',face,'Vertices',vertex,'FaceColor','flat','EdgeAlpha',0.1);
    elseif strcmp(color_type,'vertex')
        h = patch('Faces',face,'Vertices',vertex,'FaceColor','interp','EdgeAlpha',0.1);
    end
elseif strcmp(movie_type,'triangular_mesh_rotor') % need to assign activation color to vertices
    vertex = data.vertex;
    face = data.face;
    face_type = data.face_type;    
    face_center = data.face_center;
    h = patch('Faces',face,'Vertices',vertex,'FaceColor','interp','EdgeAlpha',0.1);
end

view(az,el);
axis tight equal vis3d;
axis off;
rotate3d on;

for n = frame_id
    [~,movie_color] = convert_data_to_color(data_max,data_min,data_min,signal(n,:),'red_blue');
    
    if strcmp(movie_type,'point_cloud')
        set(h,'CData',movie_color);
    elseif strcmp(movie_type,'triangular_mesh')
        % set back face black
        if strcmp(color_type,'face')
            [back_vertex_id,back_face_id] = find_vertex_face_id_on_the_back(face,vertex);
            movie_color(back_face_id,:) = 0.5;
        elseif strcmp(color_type,'vertex')
            [back_vertex_id,back_face_id] = find_vertex_face_id_on_the_back(face,vertex);
            movie_color(back_vertex_id,:) = 0.5;
        end
        
        set(h,'FaceVertexCData',movie_color);
    elseif strcmp(movie_type,'triangular_mesh_rotor')
        cla;
        h = patch('Faces',face,'Vertices',vertex,'FaceColor','interp','EdgeAlpha',0.1);
        set(h,'FaceVertexCData',movie_color);
        
        face_id_cw = face_type{n}.cw(:);
        face_id_ccw = face_type{n}.ccw(:);
        
        hold on;
        scatter3(face_center(face_id_cw,1),face_center(face_id_cw,2),face_center(face_id_cw,3),50,'r','filled');
        scatter3(face_center(face_id_ccw,1),face_center(face_id_ccw,2),face_center(face_id_ccw,3),50,'b','filled');
        hold off;
    end    
    
    drawnow;
    
    if length(frame_id) > 1 && save_flag == 1
        title(['time: ',num2str(n-1),' ms']);
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
end

if length(frame_id) > 1 && save_flag == 1
    close(writerObj);
    close;
    cd(directory.home_dir);
end
    
end
    