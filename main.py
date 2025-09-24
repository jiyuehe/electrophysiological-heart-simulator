# %%
# load libraries
# --------------------------------------------------
import os
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import utils.fn_set_axes_equal as fn_set_axes_equal
import utils.fn_create_pacing_signal as fn_create_pacing_signal
import utils.fn_equation_parts as fn_equation_parts
import utils.fn_compute_simulation as fn_compute_simulation
import utils.fn_create_phase as fn_create_phase
import utils.fn_convert_data_to_color as fn_convert_data_to_color
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter

# %% 
# load the .mat data file
# --------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
os.chdir(script_dir) # change the working directory

data_path = "data/"
mat_data = scipy.io.loadmat(data_path + "heart_example.mat")

voxel = mat_data['data']['geometry'][0,0]['edited'][0,0]['volume'][0,0]['voxel'][0,0]
neighbor_id_2d = mat_data['data']['geometry'][0,0]['edited'][0,0]['volume'][0,0]['voxel_based_voxels'][0,0].astype(np.int32) -1 # -1 is to convert Matlab 1-based index to Python 0-based index
Delta = mat_data['data']['geometry'][0,0]['edited'][0,0]['volume'][0,0]['delta'][0,0]

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
# settings
# --------------------------------------------------
# simulation parameters
dt = 0.05 # ms. if dt is not small enough, simulation will give NaN. Generally, if c <= 1.0, can use dt = 0.05
t_final = 200 # ms
pacing_start_time = 1 # ms
pacing_cycle_length = 250 # ms
c = 1 # diffusion coefficient. c = 1 is good for atrium
pacing_voxel_id = 100

# parameters of the heart model
n_voxel = voxel.shape[0] 

model_flag = 1
if model_flag == 1: # Mitchell-Schaeffer
    tau_in_voxel = np.ones(n_voxel) * 0.3
    tau_out_voxel = np.ones(n_voxel) * 6
    tau_open_voxel = np.ones(n_voxel) * 120
    tau_close_voxel = np.ones(n_voxel) * 80
    v_gate = 0.13
    v_gate_voxel = np.ones(n_voxel) * v_gate
elif model_flag == 2: # Aliev–Panfilov
    a = 1 # just put something here for now

# %% 
# compute simulation
# --------------------------------------------------
do_flag = 0 # 1: compute simulation, 0: load existing result
if do_flag == 1:
    # fiber orientations
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

    # compute heart model equation parts
    P_2d = fn_equation_parts.execute(n_voxel, D0, neighbor_id_2d, tau_open_voxel, tau_close_voxel, tau_in_voxel, tau_out_voxel, v_gate_voxel, c)

    # create the pacing signal
    pacing_signal = fn_create_pacing_signal.execute(dt, t_final, pacing_start_time, pacing_cycle_length)

    debug_plot = 0
    if debug_plot == 1: # plot pacing signal
        t = np.arange(dt, t_final + dt, dt)
        plt.figure()
        plt.plot(t,pacing_signal, 'b')
        plt.xlabel('Time (ms)')
        plt.title('Pacing signal')
        plt.show()

    # compute simulation
    action_potential, h = fn_compute_simulation.execute(neighbor_id_2d, pacing_voxel_id, n_voxel, dt, t_final, pacing_signal, P_2d, Delta)
    np.save('result/action_potential.npy', action_potential)

    # create phase from action potential
    action_potential_phase = np.zeros_like(action_potential)
    activation_phase = np.zeros_like(action_potential)
    for id in range(action_potential.shape[0]):
        action_potential_phase[id,:], activation_phase[id,:] = fn_create_phase.execute(action_potential[id,:], v_gate)
    np.save('result/action_potential_phase.npy', action_potential_phase)
elif do_flag == 0:
    action_potential = np.load('result/action_potential.npy')
    action_potential_phase = np.load('result/action_potential_phase.npy')

debug_plot = 0
if debug_plot == 1:
    # action potential of some voxel
    voxel_id = 1000
    plt.figure()
    plt.plot(action_potential[voxel_id, :],'b')
    plt.show()

# %% 
# display result
# --------------------------------------------------
# display simulation phase movie
movie_data = action_potential_phase
xyz = voxel
t = np.arange(1, t_final + 1, 1)
num_particles, num_time_steps = movie_data.shape
v_min = np.min(movie_data)
v_max = np.max(movie_data)

d_buffer = 5 
x_min = np.min(xyz[:,0]) - d_buffer
y_min = np.min(xyz[:,1]) - d_buffer
z_min = np.min(xyz[:,2]) - d_buffer
x_max = np.max(xyz[:,0]) + d_buffer
y_max = np.max(xyz[:,1]) + d_buffer
z_max = np.max(xyz[:,2]) + d_buffer

do_flag = 1
if do_flag == 1:
    print("display movie")

    # dictionary to store view angles for each frame
    view_angles = {}

    interval = 0.01
    plt.figure(figsize=(10, 8))
    ax = plt.axes(projection='3d')
    ax.view_init(elev = -50, azim = 100)

    for n in range(num_time_steps):
        print(n/num_time_steps)

        ax.clear()

        # assign color based on phase to each voxel
        data = action_potential_phase[:, n]
        data_min = v_min
        data_max = v_max
        data_threshold = v_min
        map_color = fn_convert_data_to_color.execute(data, data_min, data_max, data_threshold)

        ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=map_color, s=2, alpha=1)
        
        # set title with current time step
        ax.set_title(f'Time: {t[n]}/{t[-1]}')
        
        # reset axis properties
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        fn_set_axes_equal.execute(ax)
        ax.set_xlim([x_min, x_max])
        ax.set_ylim([y_min, y_max])
        ax.set_zlim([z_min, z_max])

        # capture current view angles
        elev = ax.elev  # elevation angle
        azim = ax.azim  # azimuth angle
        view_angles[n] = {'elev': elev, 'azim': azim}

        plt.pause(interval)

    # save simulation movie as mp4
    do_flag = 0
    if do_flag == 1:
        print("saving movie as mp4")

        fig = plt.figure(figsize=(10, 8))
        ax = plt.axes(projection='3d')

        def animate(n):
            print(n/num_time_steps)

            ax.clear()
            
            # assign color based on phase to each voxel
            data = action_potential_phase[:, n]
            data_min = v_min
            data_max = v_max
            data_threshold = v_min
            map_color = fn_convert_data_to_color.execute(data, data_min, data_max, data_threshold)

            ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=map_color, s=2, alpha=1)
            
            # set title with current time step
            ax.set_title(f'Time: {t[n]}/{t[-1]}')
            
            # reset axis properties
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            fn_set_axes_equal.execute(ax)
            ax.set_xlim([x_min, x_max])
            ax.set_ylim([y_min, y_max])
            ax.set_zlim([z_min, z_max])

            # restore view angle to maintain user's rotation
            ax.view_init(elev=view_angles[n]['elev'], azim=view_angles[n]['azim'])

        anim = animation.FuncAnimation(fig, animate, frames=num_time_steps, interval=10, blit=False, repeat=False)
        # the interval parameter specifies the delay between frames in milliseconds

        # save
        writer = FFMpegWriter(fps=10, bitrate=1800)
        anim.save('result/simulation movie.mp4', writer=writer)

        print("movie saved as mp4")
