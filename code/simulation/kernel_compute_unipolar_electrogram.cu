__global__ void egm_kernel(double* egm_part, const double* P, const int n_node, const int n_parts, const double* dvdx, const double* dvdy, const double* dvdz) {
    // note: the first argument is input and output argument, the argument with "const" is input argument only
    
	int n_id = threadIdx.x + (blockIdx.x * blockDim.x); // blockDim.x is the built-in variable for threads per block

    if (n_id > n_node - 1) {
		return;
	}

    // NOTE: a^3 should write as a*a*a, or will have compile error
    // matlab index starts at 1, c index starts at 0
    egm_part[n_id] = P[n_id*n_parts+1-1] / (P[n_id*n_parts+11-1]*P[n_id*n_parts+11-1]*P[n_id*n_parts+11-1]) *
        ( (P[n_id*n_parts+2-1] * dvdx[n_id] + P[n_id*n_parts+3-1] * dvdy[n_id] + P[n_id*n_parts+4-1] * dvdz[n_id]) * P[n_id*n_parts+12-1] +
        (P[n_id*n_parts+5-1] * dvdx[n_id] + P[n_id*n_parts+6-1] * dvdy[n_id] + P[n_id*n_parts+7-1] * dvdz[n_id]) * P[n_id*n_parts+13-1] +
        (P[n_id*n_parts+8-1] * dvdx[n_id] + P[n_id*n_parts+9-1] * dvdy[n_id] + P[n_id*n_parts+10-1] * dvdz[n_id]) * P[n_id*n_parts+14-1] );

	return;
}
