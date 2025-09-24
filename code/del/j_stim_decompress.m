function [pacing_voxel_id,signal] = j_stim_decompress(stimulus)

pacing_voxel_id = stimulus.pacing_voxel_id;

T = stimulus.T;
voxel_id = stimulus.voxel_id;
time_id = stimulus.time_id;
non_zero_value = stimulus.non_zero_value;

signal_temp = zeros(max(voxel_id),T);
for i = 1:length(voxel_id)
    signal_temp(voxel_id(i),time_id(i)) = non_zero_value(i);
end
signal = signal_temp(pacing_voxel_id,:);

end