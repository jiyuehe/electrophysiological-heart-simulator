import numpy as np
from numba import njit, prange # pip install numba

# vectorized computation
# --------------------------------------------------
def execute_vectorized(neighbor_id_2d, pacing_voxel_id, n_voxel, dt, t_final, pacing_signal, P_2d, Delta):
    neighbor_id = neighbor_id_2d[pacing_voxel_id, :] # add all the neighbors of the pacing voxel to be paced
    pacing_voxel_id = neighbor_id[neighbor_id != -1] # remove the -1s, which means no neighbors

    u_current = np.zeros(n_voxel) # initial value 0, set all voxel at rest
    h_current = np.ones(n_voxel)  # initial value 1, set all voxel at rest
    u_next = np.empty_like(u_current)
    h_next = np.empty_like(h_current)
    J_stim = np.zeros(n_voxel)
    diffusion_term = np.zeros(n_voxel)

    T = int(t_final/dt) # number of simulation time steps
    id_save = 0

    sim_u_voxel = np.zeros((n_voxel, t_final)) # sampling frequency at 1 kHz
    sim_h_voxel = np.zeros((n_voxel, t_final)) # sampling frequency at 1 kHz

    neighbor_id_2d_2 = neighbor_id_2d.copy() # NOTE: without .copy(), changes of neighbor_id_2d_2 will also change neighbor_id_2d
    neighbor_id_2d_2[neighbor_id_2d_2 == -1] = 0 # change -1 to 0, so that it can be used as index
        # this does not matter because the corresponding P_2d will be 0, 
        # so these terms will be eliminated anyway

    for t in range(T): 
        do_flag = 1
        if do_flag == 1 and (t % (T//10)) == 0:
            print(f'simulating {t/T*100:.1f}%')
        
        # compute diffusion term
        diffusion_term = P_2d[:, 20] / (4*Delta**2) * \
        ( \
            P_2d[:, 0] * (u_current[neighbor_id_2d_2[:, 0]] - u_current) + \
            P_2d[:, 1] * (u_current[neighbor_id_2d_2[:, 1]] - u_current) + \
            P_2d[:, 2] * (u_current[neighbor_id_2d_2[:, 2]] - u_current) + \
            P_2d[:, 3] * (u_current[neighbor_id_2d_2[:, 3]] - u_current) + \
            P_2d[:, 4] * (u_current[neighbor_id_2d_2[:, 4]] - u_current) + \
            P_2d[:, 5] * (u_current[neighbor_id_2d_2[:, 5]] - u_current) + \
            P_2d[:, 6] * (u_current[neighbor_id_2d_2[:, 0]] - u_current[neighbor_id_2d_2[:, 1]]) + \
            P_2d[:, 7] * (u_current[neighbor_id_2d_2[:, 2]] - u_current[neighbor_id_2d_2[:, 3]]) + \
            P_2d[:, 8] * (u_current[neighbor_id_2d_2[:, 4]] - u_current[neighbor_id_2d_2[:, 5]]) + \
            P_2d[:, 9] * (u_current[neighbor_id_2d_2[:, 6]] - u_current[neighbor_id_2d_2[:, 8]]) + \
            P_2d[:, 10] * (u_current[neighbor_id_2d_2[:, 9]] - u_current[neighbor_id_2d_2[:, 7]]) + \
            P_2d[:, 11] * (u_current[neighbor_id_2d_2[:, 14]] - u_current[neighbor_id_2d_2[:, 16]]) + \
            P_2d[:, 12] * (u_current[neighbor_id_2d_2[:, 17]] - u_current[neighbor_id_2d_2[:, 15]]) + \
            P_2d[:, 13] * (u_current[neighbor_id_2d_2[:, 10]] - u_current[neighbor_id_2d_2[:, 12]]) + \
            P_2d[:, 14] * (u_current[neighbor_id_2d_2[:, 13]] - u_current[neighbor_id_2d_2[:, 11]]) \
        )
        
        diffusion_term = diffusion_term.flatten()
        
        # compute the next time step value of u
        J_stim.fill(0.0) # reset values to 0s
        J_stim[pacing_voxel_id] = pacing_signal[t]
        
        u_next = u_current + dt * \
        ( \
            (h_current * (u_current**2) * (1 - u_current)) / P_2d[:, 17] - \
            u_current / P_2d[:, 18] + \
            J_stim + \
            diffusion_term \
        )
        
        # compute the next time step value of h
        h_next_1 = ((1 - h_current) / P_2d[:, 15]) * dt + h_current
        h_next_2 = (-h_current / P_2d[:, 16]) * dt + h_current
        
        id_1 = u_current < P_2d[:, 19]
        id_2 = u_current >= P_2d[:, 19]
        
        h_next[id_1] = h_next_1[id_1]
        h_next[id_2] = h_next_2[id_2]
        
        # update value
        u_current = u_next
        h_current = h_next
        
        # save value at 1 kHz
        if (t % int(1/dt)) == 0:
            sim_u_voxel[:, id_save] = u_current
            sim_h_voxel[:, id_save] = h_current
            id_save = id_save + 1

    return sim_u_voxel, sim_h_voxel

# CPU paralleled computation
# --------------------------------------------------
@njit(parallel=True)
def compute_voxel(u_current, h_current, P_2d, neighbor_id_2d_2, J_stim, dt, Delta):
    n_voxel = u_current.shape[0]
    u_next = np.empty_like(u_current)
    h_next = np.empty_like(h_current)
    diffusion_term = np.zeros(n_voxel)
    
    for n in prange(n_voxel): # parallel loop
        # compute diffusion term
        diffusion_term[n] = P_2d[n, 20] / (4*Delta**2) * \
        ( \
            P_2d[n, 0] * (u_current[neighbor_id_2d_2[n, 0]] - u_current[n]) + \
            P_2d[n, 1] * (u_current[neighbor_id_2d_2[n, 1]] - u_current[n]) + \
            P_2d[n, 2] * (u_current[neighbor_id_2d_2[n, 2]] - u_current[n]) + \
            P_2d[n, 3] * (u_current[neighbor_id_2d_2[n, 3]] - u_current[n]) + \
            P_2d[n, 4] * (u_current[neighbor_id_2d_2[n, 4]] - u_current[n]) + \
            P_2d[n, 5] * (u_current[neighbor_id_2d_2[n, 5]] - u_current[n]) + \
            P_2d[n, 6] * (u_current[neighbor_id_2d_2[n, 0]] - u_current[neighbor_id_2d_2[n, 1]]) + \
            P_2d[n, 7] * (u_current[neighbor_id_2d_2[n, 2]] - u_current[neighbor_id_2d_2[n, 3]]) + \
            P_2d[n, 8] * (u_current[neighbor_id_2d_2[n, 4]] - u_current[neighbor_id_2d_2[n, 5]]) + \
            P_2d[n, 9] * (u_current[neighbor_id_2d_2[n, 6]] - u_current[neighbor_id_2d_2[n, 8]]) + \
            P_2d[n, 10] * (u_current[neighbor_id_2d_2[n, 9]] - u_current[neighbor_id_2d_2[n, 7]]) + \
            P_2d[n, 11] * (u_current[neighbor_id_2d_2[n, 14]] - u_current[neighbor_id_2d_2[n, 16]]) + \
            P_2d[n, 12] * (u_current[neighbor_id_2d_2[n, 17]] - u_current[neighbor_id_2d_2[n, 15]]) + \
            P_2d[n, 13] * (u_current[neighbor_id_2d_2[n, 10]] - u_current[neighbor_id_2d_2[n, 12]]) + \
            P_2d[n, 14] * (u_current[neighbor_id_2d_2[n, 13]] - u_current[neighbor_id_2d_2[n, 11]]) \
        )
        
        # compute the next time step value of u
        u_next[n] = u_current[n] + dt * \
        ( \
            (h_current[n] * (u_current[n]**2) * (1 - u_current[n])) / P_2d[n, 17] - \
            u_current[n] / P_2d[n, 18] + \
            J_stim[n] + \
            diffusion_term[n] \
        )

        # compute the next time step value of h
        h_next_1 = ((1 - h_current[n]) / P_2d[n, 15]) * dt + h_current[n]
        h_next_2 = (-h_current[n] / P_2d[n, 16]) * dt + h_current[n]
        
        if u_current[n] < P_2d[n, 19]:
            h_next[n] = h_next_1
        elif u_current[n] >= P_2d[n, 19]:
            h_next[n] = h_next_2

    return u_next, h_next

def execute_CPU_parallel(neighbor_id_2d, pacing_voxel_id, n_voxel, dt, t_final, pacing_signal, P_2d, Delta):
    neighbor_id = neighbor_id_2d[pacing_voxel_id, :] # add all the neighbors of the pacing voxel to be paced
    pacing_voxel_id = neighbor_id[neighbor_id != -1] # remove the -1s, which means no neighbors

    u_current = np.zeros(n_voxel) # initial value 0, set all voxel at rest
    h_current = np.ones(n_voxel)  # initial value 1, set all voxel at rest
    u_next = np.empty_like(u_current)
    h_next = np.empty_like(h_current)
    J_stim = np.zeros(n_voxel)
    diffusion_term = np.zeros(n_voxel)

    T = int(t_final/dt) # number of simulation time steps
    id_save = 0

    sim_u_voxel = np.zeros((n_voxel, t_final)) # sampling frequency at 1 kHz
    sim_h_voxel = np.zeros((n_voxel, t_final)) # sampling frequency at 1 kHz

    neighbor_id_2d_2 = neighbor_id_2d.copy() # NOTE: without .copy(), changes of neighbor_id_2d_2 will also change neighbor_id_2d
    neighbor_id_2d_2[neighbor_id_2d_2 == -1] = 0 # change -1 to 0, so that it can be used as index
        # this does not matter because the corresponding P_2d will be 0, 
        # so these terms will be eliminated anyway

    for t in range(T): 
        do_flag = 1
        if do_flag == 1 and (t % (T//5)) == 0:
            print(f'simulating {t/T*100:.1f}%')
        
        J_stim.fill(0.0) # reset values to 0s
        J_stim[pacing_voxel_id] = pacing_signal[t]

        u_next, h_next = compute_voxel(u_current, h_current, P_2d, neighbor_id_2d_2, J_stim, dt, Delta)
        u_current = u_next
        h_current = h_next
        
        # save value at 1 kHz
        if (t % int(1/dt)) == 0:
            sim_u_voxel[:, id_save] = u_current
            sim_h_voxel[:, id_save] = h_current
            id_save = id_save + 1

    return sim_u_voxel, sim_h_voxel
