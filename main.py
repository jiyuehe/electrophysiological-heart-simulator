# %%
import os
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import fn_set_axes_equal

# %% load data
# --------------------------------------------------
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

# %% settings
# --------------------------------------------------
dt = 0.05 # ms. if dt is not small enough, simulation will give NaN. Generally, if c <= 1.0, can use dt = 0.05
t_final = 300 # ms
pacing_voxel_id = np.array([36184, 36190, 36191, 36198, 36693, 36694, 36695, 36699, 36700, 36701, 36705, 36706, 36707, 37187, 37192, 37193, 37194, 37199])
pacing_start_time = 1 # ms
pacing_cycle_length = 180 # ms
pacing_duration = 10 # ms
c = 1 # diffusion coefficient. c = 1 is good for atrium

# create the pacing signal
# --------------------------------------------------
pacing_duration_time_steps = pacing_duration/dt # make sure it is n ms no matter what dt is
J_stim_value = 20 # pacing strength. 20 is good. 10 is not large enough if space step is 0.1 mm and time step is 0.001
pacing_starts = np.arange(pacing_start_time/dt, t_final/dt - pacing_duration_time_steps + 1, pacing_cycle_length/dt)
pacing_ends = pacing_starts + pacing_duration_time_steps - 1

t_step = len(np.arange(dt, t_final + dt, dt)) # time steps
pacing_signal = np.zeros(t_step)
for n in range(len(pacing_starts)):
    pacing_signal[int(pacing_starts[n]):int(pacing_ends[n])+1] = J_stim_value

debug_plot = 0
if debug_plot == 1:
    t = np.arange(dt, t_final + dt, dt)
    plt.figure()
    plt.plot(t,pacing_signal, 'b')
    plt.xlabel('Time (ms)')
    plt.title('Pacing signal')
    plt.show()

# parameters
n_voxel = voxel.shape[0] 
tau_in_voxel = np.ones((n_voxel, 1)) * 0.3
tau_out_voxel = np.ones((n_voxel, 1)) * 6
tau_open_voxel = np.ones((n_voxel, 1)) * 120
tau_close_voxel = np.ones((n_voxel, 1)) * 80
v_gate_voxel = np.ones((n_voxel, 1)) * 0.13

# fiber orientation
fiber_flag = 0 # 0: no fiber, 1: fiber
r = [] # no fiber
fiber_orientation = [] # no fiber
D0 = [None] * n_voxel  # Create list of None values (equivalent to cell array)
for n in range(n_voxel):  # 0-based indexing in Python
    if fiber_flag == 1:
        e1 = fiber_orientation[n, :].reshape(-1, 1)  # Make column vector
        D0[n] = r * np.eye(3) + (1-r) * (e1 @ e1.T)  # @ is matrix multiplication
    elif fiber_flag == 0:
        # here r = 1
        D0[n] = np.eye(3)

# equation parts

