function handles = assign_region_flags(handles)

mesh.vertex = handles.vertex;
mesh.face = handles.face;
mesh.edge = handles.edge;
mesh.vertex_based_vertices = handles.vertex_based_vertices;

% obtain selected vertices
c_info = getCursorInfo(handles.dcm_obj);
scatter_vertex_id = [];
for i = 1:length(c_info)
    pos = c_info(i).Position;
    
    for j = 1:size(mesh.vertex,1)
        if isequal(pos,mesh.vertex(j,:))
            scatter_vertex_id(i) = j;
            break;
        end
    end
end

% obtain center vertex, the last vertex clicked by user
center_vertex = scatter_vertex_id(1); % the last one clicked is stored at the first place

% obtain the scattered segments from the selected vertices
scatter_vertex_id(1) = []; % delete the center vertex, so that remains boundry vertices
scattered_segment = [];
for i = 1:length(scatter_vertex_id)-1
    scattered_segment(i,:) = [scatter_vertex_id(i) scatter_vertex_id(i+1)];
end
scattered_segment(end+1,:) = [scatter_vertex_id(end) scatter_vertex_id(1)];

% obtain every vertices in between each scattered segments
step = 0;
steps = size(scattered_segment,1);
h = waitbar(step,'Please wait. Computing...','WindowStyle', 'modal');
vertices_XYZ = mesh.vertex;
boundry_vertex_id = [];
cost_matrix = handles.geometry.cost_matrix;
for i = 1:size(scattered_segment,1)    
    step = step +1;
    
    [~,trace{i}] = dijkstra(vertices_XYZ,mesh.edge,scattered_segment(i,1),scattered_segment(i,2),cost_matrix);
    vertices2add = trace{i};
    vertices2add(1) = [];
    boundry_vertex_id(end+1:end+length(vertices2add)) = vertices2add;
    
    % assign coordinate infinity, so that it won't be used again
    boundryVertex2 = boundry_vertex_id(1:end-1);
    vertices_XYZ(boundryVertex2,:) = repmat([Inf Inf Inf],length(boundryVertex2),1);
    
    waitbar(step/steps);
end
close(h);
boundryVertexXYZ = mesh.vertex(boundry_vertex_id,:);

handles.center_vertex = center_vertex;
handles.boundry_vertex_id = boundry_vertex_id;
handles.boundryVertexXYZ = boundryVertexXYZ;

% obtain vertices inside the boundry, that is vertices within the selected region
vertex_ids_selected = handles.boundry_vertex_id;
vertex_ids_selected(end+1) = handles.center_vertex;
expandList1 = handles.center_vertex;
while ~isempty(expandList1)
    expandList2 = [];    
    for i = 1:length(expandList1)
        % obtain vertices around the vertex
        vertices = mesh.vertex_based_vertices(expandList1(i),:);
        vertices(isnan(vertices)) = [];
        
        % delete vertices already included
        for j = 1:length(vertices)
            if any(vertices(j)==vertex_ids_selected)
                vertices(j) = NaN;
            end
        end
        vertices(isnan(vertices)) = [];
        
        expandList2(end+1:end+length(vertices)) = vertices;
        vertex_ids_selected(end+1:end+length(vertices)) = vertices;
    end    
    expandList1 = expandList2;
end

handles.vertex_flag(vertex_ids_selected) = handles.vertex_flag_number;

end