# %%
# load libraries
# --------------------------------------------------
import codes
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import plotly.graph_objects as go # pip install plotly. For 3D interactive plot: triangular mesh, and activation movie
import plotly.io as pio
pio.renderers.default = "browser" # simulation result mesh display in internet browser

script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
os.chdir(script_dir) # change the working directory

# %% 
# load the .mat data file
# --------------------------------------------------
data_path = script_dir + "/data/"
voxel, neighbor_id_2d, Delta, voxel_for_each_vertex, vertex, face = codes.processing.prepare_geometry.execute(data_path)

# %% 
# simulation parameters
# --------------------------------------------------
dt = 0.05 # ms. if dt is not small enough, simulation will result nan. Generally, if c <= 1.0, can use dt = 0.05
t_final = 300 # ms. NOTE: need to be at least long enough to have two pacings, or cannot compute phase from action potential
pacing_voxel_id = 100
pacing_start_time = 1 # ms
pacing_cycle_length = 250 # ms

# parameters of the heart model
n_voxel = voxel.shape[0] 

model_flag = 1 # 1: Mitchell-Schaeffer, 2: Aliev–Panfilov
if model_flag == 1: # Mitchell-Schaeffer
    parameter = {}
    parameter['tau_in_voxel'] = np.ones(n_voxel) * 0.3
    parameter['tau_out_voxel'] = np.ones(n_voxel) * 6
    parameter['tau_open_voxel'] = np.ones(n_voxel) * 120
    parameter['tau_close_voxel'] = np.ones(n_voxel) * 80
    c = 1 # diffusion coefficient. c = 1 is good for atrium
    parameter['c_voxel'] = c * np.ones(n_voxel)
    v_gate = 0.13
    parameter['v_gate_voxel'] = np.ones(n_voxel) * v_gate
elif model_flag == 2: # Aliev–Panfilov
    parameter = {}
    c = 0.1 # diffusion coefficient
    parameter['c_voxel'] = c * np.ones(n_voxel)
    v_gate = 0.13
    parameter['v_gate_voxel'] = np.ones(n_voxel) * v_gate

# %% 
# compute simulation
# --------------------------------------------------
do_flag = 1 # 1: compute, 0: load existing result
if do_flag == 1:
    # fiber orientations
    D0 = codes.simulation.fibers.execute(n_voxel)

    # compute heart model equation parts
    P_2d = codes.compute_equation_parts.execute(n_voxel, D0, neighbor_id_2d, parameter, model_flag)

    # create the pacing signal
    pacing_signal = codes.create_pacing_signal.execute(dt, t_final, pacing_start_time, pacing_cycle_length, model_flag)

    # compute simulation
    action_potential, h = codes.compute_simulation.execute_CPU_parallel(neighbor_id_2d, pacing_voxel_id, n_voxel, dt, t_final, pacing_signal, P_2d, Delta, model_flag)
    np.save('result/action_potential.npy', action_potential)

    # compute unipolar electrogram
    electrode_id = [0, 5000, 10000]
    electrode_xyz = voxel[electrode_id, :]
    electrogram_unipolar = codes.compute_unipolar_electrogram.execute_CPU_parallel(electrode_xyz, voxel, D0, parameter['c_voxel'], action_potential, Delta, neighbor_id_2d)
    np.save('result/electrogram_unipolar.npy', electrogram_unipolar)
elif do_flag == 0:
    action_potential = np.load('result/action_potential.npy')
    electrogram_unipolar = np.load('result/electrogram_unipolar.npy')

debug_plot = 0
if debug_plot == 1:
    # action potential
    plt.figure()
    plt.plot(action_potential[electrode_id, :].T)
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (scaled)')
    plt.title('examples of simulated action potential')
    plt.savefig('result/action_potential.png')

    # unipolar electrogram
    plt.figure()
    # e_id = 0
    # plt.plot(electrogram_unipolar[e_id,:], 'b')
    plt.plot(electrogram_unipolar.T)
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (scaled)')
    plt.title('examples of simulated unipolar electrogram')
    plt.savefig('result/unipolar_electrogram.png')

    plt.show()

#%%
# display result
# --------------------------------------------------
# create phase from action potential
do_flag = 1 # 1: compute, 0: load existing result
if do_flag == 1:
    
    action_potential_phase = np.zeros_like(action_potential)
    activation_phase = np.zeros_like(action_potential)
    for id in range(action_potential.shape[0]):
        action_potential_phase[id,:], activation_phase[id,:] = codes.create_phase.execute(action_potential[id,:], v_gate)
    np.save('result/action_potential_phase.npy', action_potential_phase)
elif do_flag == 0:
    action_potential_phase = np.load('result/action_potential_phase.npy')

debug_plot = 0
if debug_plot == 1:
    voxel_id = 2000
    ap = action_potential[voxel_id,:]
    ap_phase = action_potential_phase[voxel_id,:]
    plt.figure()
    plt.plot(ap, 'b')
    plt.plot(ap_phase, 'g')
    plt.show()

# %%
# activation activation movie display using plotly
do_flag = 1
if do_flag == 1:
    # activation activation movie display on mesh using plotly
    movie_data = action_potential_phase[voxel_for_each_vertex,:]
    data_min = np.min(movie_data)
    data_max = np.max(movie_data)
    data_threshold = data_min
    map_color = {}
    n_time = movie_data.shape[1]
    for n in range(n_time):
        if (n % (n_time//5)) == 0:
            print(f'compute color map {n/n_time*100:.1f}%')
        data = movie_data[:, n]
        color = codes.convert_data_to_color.execute(data, data_min, data_max, data_threshold)
        map_color[n] = color
    codes.display_activation_movie.execute_on_mesh(vertex, face, map_color)

    # activation movie display on volume using plotly
    do_flag = 0
    if do_flag == 1:
        movie_data = action_potential_phase[:,0:50] # if too large of data, will not display
        data_min = np.min(movie_data)
        data_max = np.max(movie_data)
        data_threshold = data_min
        map_color = {}
        n_time = movie_data.shape[1]
        for n in range(n_time):
            if (n % (n_time//5)) == 0:
                print(f'compute color map {n/n_time*100:.1f}%')
            data = movie_data[:, n]
            color = codes.convert_data_to_color.execute(data, data_min, data_max, data_threshold)
            map_color[n] = color
        codes.display_activation_movie.execute_on_volume(voxel, map_color)

# activation phase movie display on volume using matplotlib, can save as mp4
do_flag = 0
if do_flag == 1: 
    movie_data = action_potential_phase
    data_min = np.min(movie_data)
    data_max = np.max(movie_data)
    data_threshold = data_min
    map_color = {}
    n_time = movie_data.shape[1]
    for n in range(n_time):
        if (n % (n_time//5)) == 0:
            print(f'compute color map {n/n_time*100:.1f}%')
        data = movie_data[:, n]
        color = codes.convert_data_to_color.execute(data, data_min, data_max, data_threshold)
        map_color[n] = color

    d_buffer = 5 
    x_min = np.min(voxel[:,0]) - d_buffer
    y_min = np.min(voxel[:,1]) - d_buffer
    z_min = np.min(voxel[:,2]) - d_buffer
    x_max = np.max(voxel[:,0]) + d_buffer
    y_max = np.max(voxel[:,1]) + d_buffer
    z_max = np.max(voxel[:,2]) + d_buffer

    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection='3d')
    ax.view_init(elev = -50, azim = 100)
    plot_handle = ax.scatter(voxel[:, 0], voxel[:, 1], voxel[:, 2], c=map_color[0], s=2, alpha=1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_zlim([z_min, z_max])
    codes.set_axes_equal.execute(ax)

    pause_interval = 0.001
    view_angles = {} # dictionary to store view angles for each frame
    n_time = movie_data.shape[1]
    for n in range(n_time):
        if (n % (n_time//10)) == 0:
            print(f'playing movie {n/n_time*100:.1f}%')

        plot_handle.set_facecolor(map_color[n]) # set color based on phase to each voxel
        ax.set_title(f'Time: {n}/{n_time} ms') # set title with current time step

        # capture current view angles
        elev = ax.elev  # elevation angle
        azim = ax.azim  # azimuth angle
        view_angles[n] = {'elev': elev, 'azim': azim}

        plt.pause(pause_interval)

    # save simulation movie as mp4
    do_flag = 1
    if do_flag == 1:
        print("saving movie as mp4")

        def animate(n):
            if (n % (n_time//10)) == 0:
                print(f'saving movie {n/n_time*100:.1f}%')

            plot_handle.set_facecolor(map_color[n]) # set color based on phase to each voxel
            ax.set_title(f'Time: {n}/{n_time} ms') # set title with current time step

            ax.view_init(elev=view_angles[n]['elev'], azim=view_angles[n]['azim']) # restore view angle

        anim = animation.FuncAnimation(fig, animate, frames=n_time, interval=10, blit=False, repeat=False)
        # the interval parameter specifies the delay between frames in milliseconds

        # save
        writer = FFMpegWriter(fps=10, bitrate=1800)
        anim.save('result/simulation movie.mp4', writer=writer)

        print("movie saved as mp4")
