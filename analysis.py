def animate(s):
    
    print('Processing frame', s, 'of', np.max(frames))
        
    #remove previous contourplot by removing all collections in each ct object
    #Using ax.clear() also works, but this removes also axes labels and titles
    global ct1, ct2, ct3
    
    for coll in ct1.collections:
        coll.remove()
        
    for coll in ct2.collections:
        coll.remove()

    for coll in ct3.collections:
        coll.remove()
    
    ct1 = ax1.contourf(x, y, vort[s], 100)
    ct2 = ax2.contourf(x, y, r[s], 100)
    ct3 = ax3.contourf(x, y, covar[s], 100) 
    
def auto_correlation_function(X, max_lag):
    """
    Compute the autocorrelation of X over max_lag time steps
    
    Parameters:
        - X (array, size (N,)): the samples from which to compute the ACF 
        - max_lag (int): the max number of time steps, determines max 
          lead time
          
    Returns:
        - R (array): array of ACF values
    """
    lags = np.arange(1, max_lag)
    R = np.zeros(max_lag)
    R[0] = 1.0
    
    idx = 0
    
    print('Computing auto-correlation function')
    
    #for every lag, compute autocorrelation:
    # R = E[(X_t - mu_t)*(X_s - mu_s)]/(std_t*std_s)
    for lag in lags:
    
        X_t = X[0:-lag]
        X_s = X[lag:]
    
        mu_t = np.mean(X_t)
        std_t = np.std(X_t)
        mu_s = np.mean(X_s)
        std_s = np.std(X_s)
    
        R[idx] = np.mean((X_t - mu_t)*(X_s - mu_s))/(std_t*std_s)
        idx += 1
        
    print('done')

    return R
    
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
import h5py
from matplotlib.animation import FuncAnimation, writers

plt.close('all')
plt.rcParams['image.cmap'] = 'seismic'

#open file dialog, select .hdf5 data
root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename(filetypes=[("HDF5 files", "*.hdf5")])

#load data
h5f = h5py.File(file_path, 'r')
print('Loaded', h5f.keys())

#eddy forcing
r = h5f['r_hat_nm1'][()]
#resolved vorticity
vort = h5f['w_hat_n_LF'][()]
#Zanna covariate (see solver.py)
covar = h5f['covar_zanna_hat_n']

#time scale
Omega = 7.292*10**-5
day = 24*60**2*Omega
#time step model
dt = 0.01
#time step data
d_tau = 0.25
#store 1 sample every 'store_frame_rate' steps of the model
store_frame_rate = np.floor(d_tau*day/dt).astype('int')

#number of samples
n_samples = vort.shape[0]

#generate 2D grid
I = 7
N = 2**I
h = 2*np.pi/N
axis = h*np.arange(1, N+1)
axis = np.linspace(0, 2.0*np.pi, N)
[x , y] = np.meshgrid(axis , axis)

#flags
make_movie = False
compute_acf = True

if make_movie == True:

    #make movie of vorticity and jacobian fields
    fig = plt.figure(figsize=[12,4])
    ax1 = fig.add_subplot(131, xlabel='x', ylabel='y', title='vorticity')
    ax2 = fig.add_subplot(132, xlabel='x', ylabel='y', title='SGS term')
    ax3 = fig.add_subplot(133, xlabel='x', ylabel='y', title='Covariate')

    plt.tight_layout()
    
    vort_contour_levels = np.linspace(-0.7, 0.7, 100)
    
    ct1 = ax1.contourf(x, y, vort[0], 100)
    ct2 = ax2.contourf(x, y, r[0], 100)
    ct3 = ax2.contourf(x, y, covar[0], 100)
    
    #which frames to include in movie
    frames = np.arange(0, 0.1*n_samples, 2).astype('int')
    
    # Set up formatting for the movie files
    anim = FuncAnimation(fig, animate, frames=frames)    
    Writer = writers['ffmpeg']
    writer = Writer(fps=10)#, bitrate=1800)
    anim.save('flow.mp4', writer = writer, dpi=100)
    
if compute_acf:
    
    #choose spatial index
    I = 10; J = 10
    
    fig2 = plt.figure(figsize = [8, 4])
    ax1_acf = fig2.add_subplot(121, title=r'ACF vorticity', xlabel=r'time [days]')
    ax2_acf = fig2.add_subplot(122, title=r'ACF eddy forcing', xlabel=r'time [days]')
    
    lags = 500
    acf_r = auto_correlation_function(r[:, I, J], lags)
    acf_vort = auto_correlation_function(vort[:, I, J], lags)
    
    #time
    tau = np.arange(lags)*d_tau
    
    ax1_acf.plot(tau, acf_vort)
    ax2_acf.plot(tau, acf_r)
    
    plt.tight_layout()
    
plt.show()