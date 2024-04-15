# from ..basic import *
#from .plot_utils_style import *
import matplotlib
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, LogLocator
import scipy.stats.qmc as qmc
import numpy as np

def svd_mode_selector(data, tolerance=1e-3, modes=False):
    """
    Selects the number of singular value decomposition (SVD) modes based on a tolerance.
    
    Parameters:
    - data: The input data for SVD.
    - tolerance: The threshold for cumulative energy content in the SVD spectrum.
    - modes: If True, prints the number of selected modes.
    
    Returns:
    - The number of selected modes and the matrix of SVD left singular vectors.
    """

    # Convert input data to a NumPy array and compute SVD
    data_array = np.asarray(data)
    U, singular_values, _ = np.linalg.svd(data_array.T, full_matrices=False)
    singular_values_cumsum = np.cumsum(singular_values) / np.sum(singular_values)
    singular_values_cumsum_tol = 1.0-singular_values_cumsum
    singular_values_cumsum_tol[singular_values_cumsum_tol<np.finfo(float).eps] = np.finfo(float).eps

    # Determine the number of modes where the cumulative sum meets the tolerance
    selected_indices = np.where(singular_values_cumsum_tol < tolerance)[0]
    num_selected_modes = selected_indices[0] + 1 if selected_indices.size > 0 else 1
    
    # Plot the cumulative sum of singular values
    max_ticker = 7
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.semilogy(np.arange(1,len(singular_values)+1), singular_values_cumsum_tol, 's-', color='orange')
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=4))

    ax.axhline(y=tolerance, color="black", linestyle="--")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=max_ticker))
    ax.autoscale(tight=True)
    ax.margins(x=0.02,y=0.02)
    #plt.autoscale(tight=rue)
    
    # Print the number of modes if requested
    if modes:
        print(f"Number of modes selected: {num_selected_modes}")
    
    return num_selected_modes, U


def train_test_split(N_snap, N_sel=None,train_percentage = 0.8, seed = 42):

    np.random.seed(seed)

    # Generate a random permutation of indices from 0 to data_size - 1
    indices = np.random.permutation(N_snap)
    
    if N_sel is not None:
        indices = np.random.choice(indices, size=min(N_sel, N_snap), replace=False)


    # Calculate the number of samples in the training set
    train_set_size = int(len(indices) * train_percentage)

    # Initialize boolean masks
    train_mask = np.zeros(N_snap, dtype=bool)
    test_mask = np.zeros(N_snap, dtype=bool)

    # Set the first train_set_size indices to True for the training mask
    train_mask[indices[:train_set_size]] = True

    # Set the remaining indices to True for the testing mask
    test_mask[indices[train_set_size:]] = True
    
    return train_mask, test_mask


def train_test_split_sobol(N_snap, N_sel=None, train_percentage=0.8):
    # Sobol sequence generator
    sobol = qmc.Sobol(d=1)  # 1D Sobol sequence

    # Generate Sobol sequence points
    sobol_points = sobol.random_base2(m=int(np.ceil(np.log2(N_snap))))
    
    # Scale and convert to indices
    indices = (sobol_points.flatten() * N_snap).astype(int)
    
    # Ensure uniqueness and sufficiency of indices
    indices = np.unique(indices)[:N_snap]
    
    while len(indices) < N_snap:
        extra_points = sobol.random(n=N_snap - len(indices))
        extra_indices = (extra_points.flatten() * N_snap).astype(int)
        indices = np.unique(np.concatenate([indices, extra_indices]))[:N_snap]

    if N_sel is not None:
        # Select a subset of indices if N_sel is specified
        sobol_subset = qmc.Sobol(d=1)
        subset_points = sobol_subset.random_base2(m=int(np.ceil(np.log2(N_sel))))
        subset_indices = (subset_points.flatten() * len(indices)).astype(int)
        indices = indices[np.unique(subset_indices)[:N_sel]]

    # Calculate the number of samples in the training set
    train_set_size = int(N_snap * train_percentage)

    # Initialize boolean masks
    train_mask = np.zeros(N_snap, dtype=bool)
    test_mask = np.zeros(N_snap, dtype=bool)

    # Set the appropriate indices to True for the training and testing masks
    train_mask[indices[:train_set_size]] = True
    test_mask[indices[train_set_size:]] = True

    return train_mask, test_mask