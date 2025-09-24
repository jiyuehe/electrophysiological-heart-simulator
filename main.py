# %%
import os
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import utils.fn_set_axes_equal as fn_set_axes_equal
import utils.fn_create_pacing_signal as fn_create_pacing_signal
import utils.fn_equation_parts as fn_equation_parts
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter

# %% load the .mat data file
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

# %% settings
# --------------------------------------------------
# simulation parameters
dt = 0.05 # ms. if dt is not small enough, simulation will give NaN. Generally, if c <= 1.0, can use dt = 0.05
t_final = 400 # ms
pacing_start_time = 1 # ms
pacing_cycle_length = 190 # ms
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
    v_gate_voxel = np.ones(n_voxel) * 0.13
elif model_flag == 2: # Aliev–Panfilov
    a = 1 # just put something here for now

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

# %% heart model equation parts
# --------------------------------------------------
neighbor_id_2d_2 = neighbor_id_2d # neighbor indices
neighbor_id_2d_2[neighbor_id_2d_2 == -1] = 0 # change -1 to 0, so that it can be used as index in Python
P_2d = fn_equation_parts.execute(n_voxel, D0, neighbor_id_2d_2, tau_open_voxel, tau_close_voxel, tau_in_voxel, tau_out_voxel, v_gate_voxel, c)

# %% create the pacing signal
# --------------------------------------------------
pacing_signal = fn_create_pacing_signal.execute(dt, t_final, pacing_start_time, pacing_cycle_length)

debug_plot = 0
if debug_plot == 1:
    t = np.arange(dt, t_final + dt, dt)
    plt.figure()
    plt.plot(t,pacing_signal, 'b')
    plt.xlabel('Time (ms)')
    plt.title('Pacing signal')
    plt.show()

# %% compute simulation
# --------------------------------------------------
neighbor_id = neighbor_id_2d[pacing_voxel_id, :] # add all the neighbors of the pacing voxel to be paced
pacing_voxel_id = neighbor_id[neighbor_id != -1] # remove the -1s

u_current = np.zeros(n_voxel) # initial value 0, set all voxel at rest
h_current = np.ones(n_voxel)  # initial value 1, set all voxel at rest
u_next = np.zeros(n_voxel)
h_next = np.zeros(n_voxel)

T = int(t_final/dt) # number of simulation time steps
id_save = 0

sim_u_voxel = np.zeros((n_voxel, t_final)) # sampling frequency at 1 kHz
sim_h_voxel = np.zeros((n_voxel, t_final)) # sampling frequency at 1 kHz

for t in range(T): 
    do_flag = 1
    if do_flag == 1 and (t % (T//10)) == 0:
        print(f'simulation {t/T*100:.1f}%')
    
    # compute diffusion term
    diffusion_term = P_2d[:, 20].flatten() / (4*Delta**2) * \
    ( \
        P_2d[:, 0].flatten() * (u_current[neighbor_id_2d_2[:, 0]] - u_current).flatten() + \
        P_2d[:, 1].flatten() * (u_current[neighbor_id_2d_2[:, 1]] - u_current).flatten() + \
        P_2d[:, 2].flatten() * (u_current[neighbor_id_2d_2[:, 2]] - u_current).flatten() + \
        P_2d[:, 3].flatten() * (u_current[neighbor_id_2d_2[:, 3]] - u_current).flatten() + \
        P_2d[:, 4].flatten() * (u_current[neighbor_id_2d_2[:, 4]] - u_current).flatten() + \
        P_2d[:, 5].flatten() * (u_current[neighbor_id_2d_2[:, 5]] - u_current).flatten() + \
        P_2d[:, 6].flatten() * (u_current[neighbor_id_2d_2[:, 0]] - u_current[neighbor_id_2d_2[:, 1]]).flatten() + \
        P_2d[:, 7].flatten() * (u_current[neighbor_id_2d_2[:, 2]] - u_current[neighbor_id_2d_2[:, 3]]).flatten() + \
        P_2d[:, 8].flatten() * (u_current[neighbor_id_2d_2[:, 4]] - u_current[neighbor_id_2d_2[:, 5]]).flatten() + \
        P_2d[:, 9].flatten() * (u_current[neighbor_id_2d_2[:, 6]] - u_current[neighbor_id_2d_2[:, 8]]).flatten() + \
        P_2d[:, 10].flatten() * (u_current[neighbor_id_2d_2[:, 9]] - u_current[neighbor_id_2d_2[:, 7]]).flatten() + \
        P_2d[:, 11].flatten() * (u_current[neighbor_id_2d_2[:, 14]] - u_current[neighbor_id_2d_2[:, 16]]).flatten() + \
        P_2d[:, 12].flatten() * (u_current[neighbor_id_2d_2[:, 17]] - u_current[neighbor_id_2d_2[:, 15]]).flatten() + \
        P_2d[:, 13].flatten() * (u_current[neighbor_id_2d_2[:, 10]] - u_current[neighbor_id_2d_2[:, 12]]).flatten() + \
        P_2d[:, 14].flatten() * (u_current[neighbor_id_2d_2[:, 13]] - u_current[neighbor_id_2d_2[:, 11]]).flatten() \
    ) \
    
    # compute the next time step value of u
    J_stim = np.zeros(n_voxel)
    J_stim[pacing_voxel_id] = pacing_signal[t]
    
    u_next = \
    ( \
        (h_current * (u_current**2) * (1 - u_current)).flatten() / P_2d[:, 17].flatten() - \
        u_current.flatten() / P_2d[:, 18].flatten() + \
        J_stim + \
        diffusion_term \
    ) \
    * dt + u_current \
    
    # compute the next time step value of h
    h_next_1 = ((1 - h_current) / P_2d[:, 15]) * dt + h_current
    h_next_2 = (-h_current / P_2d[:, 16]) * dt + h_current
    
    id_1 = u_current < P_2d[:, 19]
    id_2 = u_current >= P_2d[:, 19]
    
    h_next[id_1] = h_next_1[id_1]
    h_next[id_2] = h_next_2[id_2]
    
    # update value
    u_current = u_next.flatten()
    h_current = h_next.flatten()
    
    # save value at 1 kHz
    if (t % int(1/dt)) == 0:
        sim_u_voxel[:, id_save] = u_current
        sim_h_voxel[:, id_save] = h_current
        id_save += 1

# %% display result
debug_plot = 0
if debug_plot == 1:
    # action potential of some voxel
    voxel_id = 100
    plt.figure()
    plt.plot(sim_u_voxel[voxel_id, :])
    plt.show()

# display simulation movie
voltage = sim_u_voxel
xyz = voxel
t = np.arange(1, t_final + 1, 1)

num_particles, num_time_steps = voltage.shape
v_min = 0.13
v_max = 1.0

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
    ax.view_init(elev=-90, azim=-170)

    for n in range(num_time_steps):
        print(n/num_time_steps)

        ax.clear()

        ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], 
                   c=voltage[:, n], s=2, marker='.', alpha=1, cmap='coolwarm', vmin=v_min, vmax=v_max)
            
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
            
            ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], 
                    c=voltage[:, n], s=2, marker='.', alpha=1, cmap='coolwarm', vmin=v_min, vmax=v_max)
                
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
