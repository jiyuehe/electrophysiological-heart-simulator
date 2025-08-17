__global__ void simulation_kernel(double* u_next, double* h_next, const double* u_current, const double* h_current, const int t_id, const double dt, const double delta,
	const int* indices, const double* parts, const int N, const int L, const int* voxel_id, const int* time_id, const double* non_zero_value, const int n_voxel) {
    // note: the first 2 arguments are input and output arguments, the rest arguments with "const" are input arguments only
    
	int v_id = threadIdx.x + (blockIdx.x * blockDim.x); // blockDim.x is the built-in variable for threads per block

    if (v_id > n_voxel - 1) {
		return;
	}

    // NOTE: a^3 should write as a*a*a, or will have compile error
	double diffusion_term = 1.0 / (4.0 * delta * delta) *
		(parts[v_id * 21 + 0] * (u_current[indices[v_id * 18 + 0]] - u_current[v_id]) + parts[v_id * 21 + 1] * (u_current[indices[v_id * 18 + 1]] - u_current[v_id]) +
			parts[v_id * 21 + 2] * (u_current[indices[v_id * 18 + 2]] - u_current[v_id]) + parts[v_id * 21 + 3] * (u_current[indices[v_id * 18 + 3]] - u_current[v_id]) +
			parts[v_id * 21 + 4] * (u_current[indices[v_id * 18 + 4]] - u_current[v_id]) + parts[v_id * 21 + 5] * (u_current[indices[v_id * 18 + 5]] - u_current[v_id]) +
			parts[v_id * 21 + 6] * (u_current[indices[v_id * 18 + 0]] - u_current[indices[v_id * 18 + 1]]) +
			parts[v_id * 21 + 7] * (u_current[indices[v_id * 18 + 2]] - u_current[indices[v_id * 18 + 3]]) +
			parts[v_id * 21 + 8] * (u_current[indices[v_id * 18 + 4]] - u_current[indices[v_id * 18 + 5]]) +
			parts[v_id * 21 + 9] * (u_current[indices[v_id * 18 + 6]] - u_current[indices[v_id * 18 + 8]]) + parts[v_id * 21 + 10] * (u_current[indices[v_id * 18 + 9]] - u_current[indices[v_id * 18 + 7]]) +
			parts[v_id * 21 + 11] * (u_current[indices[v_id * 18 + 14]] - u_current[indices[v_id * 18 + 16]]) + parts[v_id * 21 + 12] * (u_current[indices[v_id * 18 + 17]] - u_current[indices[v_id * 18 + 15]]) +
			parts[v_id * 21 + 13] * (u_current[indices[v_id * 18 + 10]] - u_current[indices[v_id * 18 + 12]]) + parts[v_id * 21 + 14] * (u_current[indices[v_id * 18 + 13]] - u_current[indices[v_id * 18 + 11]]));
	diffusion_term = parts[v_id * 21 + 20] * diffusion_term;

	// J_stim
	int id = -1;
	for (int i = 0; i < L; i++) {
		if (v_id == voxel_id[i]) { // is pacing voxel
			id = i;
		}
	}
	double J_stim = 0.0;
	if (id != -1) {
		J_stim = non_zero_value[id];
	}

	// update u
	u_next[v_id] = ((h_current[v_id] * u_current[v_id] * u_current[v_id] * (1 - u_current[v_id]) / parts[v_id * 21 + 17]) + (-u_current[v_id] / parts[v_id * 21 + 18]) + J_stim + diffusion_term) * dt + u_current[v_id];

	// update h
	if (u_current[v_id] < parts[v_id * 21 + 19]) {
		h_next[v_id] = ((1 - h_current[v_id]) / parts[v_id * 21 + 15]) * dt + h_current[v_id];
	}
	if (u_current[v_id] >= parts[v_id * 21 + 19]) {
		h_next[v_id] = (-h_current[v_id] / parts[v_id * 21 + 16]) * dt + h_current[v_id];
	}

	return;
}
