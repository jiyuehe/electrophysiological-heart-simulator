import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def execute(action_potential, v_gate):
    # Action potential phase: 0-1 of the action potential duration
    # ----------------------------------------------------------------------------------------------------
    debug_plot = False
    if debug_plot:
        plt.figure()
        plt.plot(action_potential, 'b')
        plt.show()

    a = np.where(action_potential > v_gate)[0]  # time index

    if debug_plot:
        plt.figure()
        plt.plot(action_potential, 'b')
        plt.plot(a, action_potential[a], 'r.')
        plt.show()

    if len(a) > 0:
        # Find discontinuities in the indices
        diff_a = np.diff(a)
        b = np.where(np.abs(diff_a) > 1)[0]
        
        # Create p array with start and end points of continuous segments
        if len(b) > 0:
            p = np.concatenate([[a[0]], a[b], a[b+1], [a[-1]]])
        else:
            p = np.array([a[0], a[-1]])
        
        p = np.sort(p)

        if debug_plot:
            plt.figure()
            plt.plot(action_potential, 'b')
            plt.scatter(p, action_potential[p], s=100, c='r', marker='.')
            plt.show()
        
        # Phase interval
        phase_interval = np.zeros((len(p)//2, 2), dtype=int)
        for n in range(len(p)//2):
            phase_interval[n, :] = p[n*2:n*2+2]

        # Phase
        L_median = int(np.ceil(np.median(phase_interval[:, 1] - phase_interval[:, 0])))
        action_potential_phase = np.zeros_like(action_potential, dtype=float)
        
        for n in range(phase_interval.shape[0]):
            m = phase_interval[n, :]
            L = len(range(m[0], m[1]))

            if m[0] == 0:  # Python uses 0-based indexing
                action_potential_phase[m[0]:m[1]] = np.linspace(m[0]/L_median, 1, L)
            elif m[0] != 0 and m[1] != len(action_potential_phase):
                action_potential_phase[m[0]:m[1]] = np.linspace(0, 1, L)
            elif m[1] == len(action_potential_phase):
                action_potential_phase[m[0]:m[1]] = np.linspace(0, min(1, L/L_median), L)
                
    else:  # if a is empty
        action_potential_phase = np.zeros_like(action_potential, dtype=float)

    if debug_plot:
        plt.figure()
        plt.plot(action_potential, 'b')
        plt.plot(action_potential_phase, 'r')
        plt.tight_layout()
        plt.show()

    # Activation phase: 0-1 in between 2 activation peaks
    # ----------------------------------------------------------------------------------------------------
    peaks, _ = find_peaks(action_potential, height=np.mean(action_potential))  # time index

    if debug_plot:
        plt.figure()
        plt.plot(action_potential, 'b')
        plt.scatter(peaks, action_potential[peaks], c='r')
        plt.show()

    if len(peaks) > 0:
        # Include t start and t final
        if peaks[0] != 0:  # Python uses 0-based indexing
            peaks = np.concatenate([[0], peaks])
        if peaks[-1] != len(action_potential) - 1:
            peaks = np.concatenate([peaks, [len(action_potential) - 1]])
        
        # Phase interval
        phase_interval = np.zeros((len(peaks)-1, 2), dtype=int)
        for n in range(len(peaks)-1):
            phase_interval[n, :] = [peaks[n], peaks[n+1]]

        # Phase
        L_median = int(np.ceil(np.median(phase_interval[:, 1] - phase_interval[:, 0])))
        activation_phase = np.zeros_like(action_potential, dtype=float)
        
        for n in range(phase_interval.shape[0]):
            m = phase_interval[n, :]
            L = len(range(m[0], m[1]))

            if m[0] == 0:  # Python uses 0-based indexing
                activation_phase[m[0]:m[1]] = np.linspace(L/L_median, 1, L)
                # NOTE: this is different than creating action_potential_phase
            elif m[0] != 0 and m[1] != len(activation_phase):
                activation_phase[m[0]:m[1]] = np.linspace(0, 1, L)
            elif m[1] == len(activation_phase):
                activation_phase[m[0]:m[1]+1] = np.linspace(0, min(1, L/L_median), L+1)
                
    else:  # if peaks is empty
        activation_phase = np.zeros_like(action_potential, dtype=float)

    if debug_plot:
        plt.figure()
        plt.plot(action_potential, 'b')
        plt.plot(activation_phase, 'r')
        plt.tight_layout()
        plt.show()

    return action_potential_phase, activation_phase