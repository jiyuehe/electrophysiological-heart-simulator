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

% assign vertex flag
cd(directory.data_dir);
if ~isfile([file_name_vertex_flag,'.mat'])
    vertex_flag = zeros(size(data.geometry.edited.vertex,1),1);
    save(file_name_vertex_flag,'vertex_flag');
end
cd(directory.home_dir);

assign_vertex_flag(directory,file_name_vertex_flag,data.geometry.edited,az,el);
