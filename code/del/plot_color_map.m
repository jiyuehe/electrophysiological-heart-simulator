function plot_color_map(face,vertex,data,data_min,data_max,data_threshold,color_wheel_type,black_back_flag,color_on_face_flag,az,el)

[~,map_color] = convert_data_to_color(data_max,data_min,data_threshold,data,color_wheel_type);

if color_on_face_flag == 0
    h = patch('Faces',face,'Vertices',vertex,'FaceColor','interp','EdgeAlpha',0.1,'FaceVertexCData',map_color);
elseif color_on_face_flag == 1
    h = patch('Faces',face,'Vertices',vertex,'FaceColor','flat','EdgeAlpha',0.1,'FaceVertexCData',map_color);
end

view(az,el);
axis off tight equal vis3d;
rotate3d on;
set(gcf,'color','w');

if black_back_flag == 1
    [back_vertex_id,back_face_id] = find_vertex_face_id_on_the_back(face,vertex);
    
    if color_on_face_flag == 0
        map_color(back_vertex_id,:) = 0.5;
    elseif color_on_face_flag == 1
        map_color(back_face_id,:) = 0.5;
    end
    
    set(h,'FaceVertexCData',map_color);
end

end
