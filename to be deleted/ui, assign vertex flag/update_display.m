function handles = update_display(handles)

vertex = handles.vertex;
face = handles.face;
vertex_flag = handles.vertex_flag;

c = [];
n_color = max(vertex_flag);
if n_color == 1
    c = 0.8 * [1 1 1];
elseif n_color == 2
    c(1,:) = 0.8 * [1 1 1];
    c(2,:) = [1 0 0];
elseif n_color == 3
    c(1,:) = 0.8 * [1 1 1];
    c(2,:) = [1 0 0];
    c(3,:) = [0 0 0];
elseif n_color > 3
    c(1,:) = 0.8 * [1 1 1];
    c(2,:) = [1 0 0];
    c(3,:) = [0 0 0];
    c(3+1:n_color,:) = jet(n_color-3);    
end

vertex_color = 0.9*ones(length(vertex_flag),3);
for i = 1:n_color
    vertex_color(vertex_flag==i,:) = repmat(c(i,:),sum(vertex_flag==i),1);
end

vertex_color(vertex_flag==21,:) = repmat(0.3*[1 1 1],sum(vertex_flag==21),1); % for initialize rotor
vertex_color(vertex_flag==22,:) = repmat(0.7*[1 1 1],sum(vertex_flag==22),1); % for initialize rotor
vertex_color(vertex_flag==23,:) = repmat([0 1 0],sum(vertex_flag==23),1); % for initialize rotor

% back_vertex_id = find_vertices_id_on_the_back(face,vertex);
% vertex_color(back_vertex_id,:) = 0;

axes(handles.axes1);
cla;
h = patch('Faces',face,'Vertices',vertex,'FaceColor','interp','EdgeAlpha',0.2);
set(h,'FaceVertexCData',vertex_color);
hold on;
pacing_vertex_id = find(vertex_flag==2);
scatter3(vertex(pacing_vertex_id,1),vertex(pacing_vertex_id,2),vertex(pacing_vertex_id,3),60,'r','filled'); % pacing vertex
hold off;
axis vis3d equal;
rotate3d on;
xlabel('x');
ylabel('y');

end
