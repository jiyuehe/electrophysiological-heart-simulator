# %%
import os
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import fn_set_axes_equal

# %%
script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
os.chdir(script_dir) # change the working directory

data_path = "data/"

# load the .mat file
mat_data = scipy.io.loadmat(data_path + "heart_example.mat")
# %%
voxel = mat_data['data']['geometry'][0,0]['edited'][0,0]['volume'][0,0]['voxel'][0,0]

debug_plot = 0
if debug_plot == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(voxel[:, 0], voxel[:, 1], voxel[:, 2], color='blue', s=1, marker='.')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1,1,1])
    fn_set_axes_equal.execute(ax)
    plt.show()

# %%
dt = 0.05  # unit: millisecond. if dt is not small enough, simulation will give NaN. Generally, if c <= 1.0, can use dt = 0.05
t_final = 100
c = 1  # diffusion coefficient. c = 1 is good for atrium
n_voxel = voxel.shape[0] 
fiber_flag = 0  # no fiber
r = []  # no fiber
fiber_orientation = []  # no fiber

# pacing
pacing_voxel_id = np.array([36184, 36190, 36191, 36198, 36693, 36694, 36695, 36699, 36700, 36701, 36705, 36706, 36707, 37187, 37192, 37193, 37194, 37199])

# duration
pacing_start_time = 1 # ms
pacing_cycle_length = 180 # ms

n = 10 # ms
pacing_duration = n/dt # make sure it is n ms no matter what dt is

pacing_starts = np.arange(pacing_start_time/dt, t_final/dt - pacing_duration + 1, pacing_cycle_length/dt)
pacing_ends = pacing_starts + pacing_duration - 1

