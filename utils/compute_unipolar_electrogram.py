import numpy as np
import matplotlib.pyplot as plt

def execute_vectorized(electrode_xyz, voxel, D0, c_voxel, action_potential, delta, voxel_based_voxels):    
    # equation parts
    node_xyz = voxel
    n_node = node_xyz.shape[0]
    n_voxel = voxel.shape[0]
    
    # Extract diffusion tensor components
    D11 = np.zeros(n_node)
    D12 = np.zeros(n_node)
    D13 = np.zeros(n_node)
    D21 = np.zeros(n_node)
    D22 = np.zeros(n_node)
    D23 = np.zeros(n_node)
    D31 = np.zeros(n_node)
    D32 = np.zeros(n_node)
    D33 = np.zeros(n_node)
    
    for n in range(n_node):
        D11[n] = D0[n][0, 0]
        D12[n] = D0[n][0, 1]
        D13[n] = D0[n][0, 2]
        D21[n] = D0[n][1, 0]
        D22[n] = D0[n][1, 1]
        D23[n] = D0[n][1, 2]
        D31[n] = D0[n][2, 0]
        D32[n] = D0[n][2, 1]
        D33[n] = D0[n][2, 2]
    
    # distance from electrode to nodes
    n_electrode = electrode_xyz.shape[0]
    l = np.zeros((n_node, n_electrode))
    l_x = np.zeros((n_node, n_electrode))
    l_y = np.zeros((n_node, n_electrode))
    l_z = np.zeros((n_node, n_electrode))
    
    for e_id in range(n_electrode):
        l_x[:, e_id] = node_xyz[:, 0] - electrode_xyz[e_id, 0]
        l_y[:, e_id] = node_xyz[:, 1] - electrode_xyz[e_id, 1]
        l_z[:, e_id] = node_xyz[:, 2] - electrode_xyz[e_id, 2]
        l[:, e_id] = np.sqrt(l_x[:, e_id]**2 + l_y[:, e_id]**2 + l_z[:, e_id]**2)
    
    l[l < 1] = 1  # electrode is at least 1 mm away from tissue
    
    # Broadcast diffusion tensor components
    D11_b = np.tile(D11[:, np.newaxis], (1, n_electrode))
    D12_b = np.tile(D12[:, np.newaxis], (1, n_electrode))
    D13_b = np.tile(D13[:, np.newaxis], (1, n_electrode))
    D21_b = np.tile(D21[:, np.newaxis], (1, n_electrode))
    D22_b = np.tile(D22[:, np.newaxis], (1, n_electrode))
    D23_b = np.tile(D23[:, np.newaxis], (1, n_electrode))
    D31_b = np.tile(D31[:, np.newaxis], (1, n_electrode))
    D32_b = np.tile(D32[:, np.newaxis], (1, n_electrode))
    D33_b = np.tile(D33[:, np.newaxis], (1, n_electrode))
    
    # Handle voxel-based neighbors for gradient computation
    px = voxel_based_voxels[:, 0].copy()  # x +
    px_0_id = px == 0
    px[px_0_id] = 1
    
    mx = voxel_based_voxels[:, 1].copy()  # x -
    mx_0_id = mx == 0
    mx[mx_0_id] = 1
    
    py = voxel_based_voxels[:, 2].copy()  # y +
    py_0_id = py == 0
    py[py_0_id] = 1
    
    my = voxel_based_voxels[:, 3].copy()  # y -
    my_0_id = my == 0
    my[my_0_id] = 1
    
    pz = voxel_based_voxels[:, 4].copy()  # z +
    pz_0_id = pz == 0
    pz[pz_0_id] = 1
    
    mz = voxel_based_voxels[:, 5].copy()  # z -
    mz_0_id = mz == 0
    mz[mz_0_id] = 1
    
    # Convert to 0-based indexing for Python
    px = px - 1
    mx = mx - 1
    py = py - 1
    my = my - 1
    pz = pz - 1
    mz = mz - 1
    
    # Compute gradients
    dvdx = (action_potential[:, px] - action_potential[:, mx]) / (2 * delta)
    dvdx[:, px_0_id] = 0
    dvdx[:, mx_0_id] = 0
    
    dvdy = (action_potential[:, py] - action_potential[:, my]) / (2 * delta)
    dvdy[:, py_0_id] = 0
    dvdy[:, my_0_id] = 0
    
    dvdz = (action_potential[:, pz] - action_potential[:, mz]) / (2 * delta)
    dvdz[:, pz_0_id] = 0
    dvdz[:, mz_0_id] = 0
    
    # compute electrogram
    T = action_potential.shape[0]
    electrogram_unipolar = np.zeros((n_electrode, T))
    
    for t_id in range(T):
        if (t_id + 1) % (T // 5) == 0:
            print(f'compute electrogram {(t_id + 1) / T * 100:.1f}%')
        
        dvdx_b = np.tile(dvdx[t_id, :].reshape(-1, 1), (1, n_electrode))
        dvdy_b = np.tile(dvdy[t_id, :].reshape(-1, 1), (1, n_electrode))
        dvdz_b = np.tile(dvdz[t_id, :].reshape(-1, 1), (1, n_electrode))
        
        electrogram_unipolar[:, t_id] = np.sum(
            c_voxel / (l**3) * (
                (D11_b * dvdx_b + D12_b * dvdy_b + D13_b * dvdz_b) * l_x +
                (D21_b * dvdx_b + D22_b * dvdy_b + D23_b * dvdz_b) * l_y +
                (D31_b * dvdx_b + D32_b * dvdy_b + D33_b * dvdz_b) * l_z
            ), axis=0
        )
    
    # adjust the electrogram magnitude
    typical_magnitude = 1  # unit: volt
    
    voltage = np.zeros(n_electrode)
    for n in range(n_electrode):
        voltage[n] = np.ptp(electrogram_unipolar[n, :])  # np.ptp is equivalent to range()
    
    voltage_mean = np.mean(voltage)
    
    scale = typical_magnitude / voltage_mean
    electrogram_unipolar = electrogram_unipolar * scale
    
    debug_plot = False
    if debug_plot:
        t = np.arange(electrogram_unipolar.shape[1])
        plt.figure()
        plt.plot(t, electrogram_unipolar.T)
        plt.show()
    
    return electrogram_unipolar