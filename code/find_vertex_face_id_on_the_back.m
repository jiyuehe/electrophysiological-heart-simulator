function [back_vertex_id,back_face_id] = find_vertex_face_id_on_the_back(face,vertex)

TR = triangulation(face,vertex);
F = faceNormal(TR);
back_face_id = find(F*campos' < 0);

back_vertex_id = [];
for i = 1:length(back_face_id)
    v_ids = face(back_face_id(i),:);
    back_vertex_id(end+1:end+3) = v_ids;
end
back_vertex_id = unique(back_vertex_id);

end