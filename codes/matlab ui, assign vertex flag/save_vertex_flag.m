function save_vertex_flag(handles)

cd(handles.directory.data_dir);
vertex_flag = handles.vertex_flag;
save(handles.file_name,'vertex_flag');
cd(handles.directory.home_dir);

end