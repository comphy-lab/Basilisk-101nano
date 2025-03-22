import numpy as np
import os
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import multiprocessing as mp
from functools import partial
import argparse  # Add at top with other imports
import sys

import matplotlib.colors as mcolors
custom_colors = ["white", "#DA8A67", "#A0522D", "#400000"]
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_hot", custom_colors)

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

def gettingfield(filename, zmin, zmax, rmin, rmax, nr):
    exe = ["./getData-2D_rot_cylinder_ink", filename, str(zmin), str(rmin), str(zmax), str(rmax), str(nr)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    # print(temp2) #debugging
    Rtemp, Ztemp, Ttemp, cstemp  = [],[],[],[]

    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            Ztemp.append(float(temp3[0]))
            Rtemp.append(float(temp3[1]))
            Ttemp.append(float(temp3[2]))
            cstemp.append(float(temp3[3]))

    R = np.asarray(Rtemp)
    Z = np.asarray(Ztemp)
    T = np.asarray(Ttemp)
    cs = np.asarray(cstemp)
    nz = int(len(Z)/nr)

    # print("nr is %d %d" % (nr, len(R))) # debugging
    print("nz is %d" % nz)

    R.resize((nz, nr))
    Z.resize((nz, nr))
    T.resize((nz, nr))
    cs.resize((nz, nr))

    # rotate by 270 degrees
    R = np.rot90(R, k=1)
    Z = np.rot90(Z, k=1)
    T = np.rot90(T, k=1)
    cs = np.rot90(cs, k=1)
    # flip the array
    R = np.flip(R, axis=0)
    Z = np.flip(Z, axis=0)
    T = np.flip(T, axis=0)
    cs = np.flip(cs, axis=0)

    return R, Z, T, cs, nz
# ----------------------------------------------------------------------------------------------------------------------

def process_timestep(ti, caseToProcess, folder, tsnap, GridsPerR, rmin, rmax, zmin, zmax, lw):
    t = tsnap * ti
    place = f"{caseToProcess}/intermediate/snapshot-{t:.4f}"
    name = f"{folder}/{int(t*1000):08d}.png"

    if not os.path.exists(folder):
        os.makedirs(folder)

    if not os.path.exists(place):
        print(f"{place} File not found!")
        return
    
    nr = int(GridsPerR * rmax)
    R, Z, T, vel, nz = gettingfield(place, zmin, zmax, rmin, rmax, nr)
    zminp, zmaxp, rminp, rmaxp = Z.min(), Z.max(), R.min(), R.max()

    # Plotting
    AxesLabel, TickLabel = 50, 20
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(19.20, 10.80))

    # Temperature plot (ax1)
    cntrl1 = ax1.imshow(T, cmap="coolwarm", interpolation='Bilinear', origin='lower', 
                       extent=[rminp, rmaxp, zminp, zmaxp], vmax=1.0, vmin=0.0)

    # Draw two concentric circles centered at (0,0)
    inner_radius = 1.0
    outer_radius = 1.5
    
    # Clip the temperature field to the outer circle region
    outer_circle1 = plt.Circle((0, 0), outer_radius, transform=ax1.transData)
    cntrl1.set_clip_path(outer_circle1)
    
    # Create a white circle for the inner region
    inner_circle1 = plt.Circle((0, 0), inner_radius, fill=True, color='white', linestyle='-', linewidth=lw)
    ax1.add_patch(inner_circle1)

    ax1.set_aspect('equal')
    ax1.set_xlim(rmin, rmax)
    ax1.set_ylim(zmin, zmax)
    ax1.set_title("Dye", fontsize=AxesLabel)
    ax1.axis('off')

    # Velocity plot (ax2)
    cntrl2 = ax2.imshow(vel, cmap="viridis", interpolation='Bilinear', origin='lower', 
                      extent=[rminp, rmaxp, zminp, zmaxp], vmin=0.0,vmax=1.0)

    # Clip the velocity field to the outer circle region
    outer_circle2 = plt.Circle((0, 0), outer_radius, transform=ax2.transData)
    cntrl2.set_clip_path(outer_circle2)
    
    # Create a white circle for the inner region
    inner_circle2 = plt.Circle((0, 0), inner_radius, fill=True, color='white', linestyle='-', linewidth=lw)
    ax2.add_patch(inner_circle2)

    ax2.set_aspect('equal')
    ax2.set_xlim(rmin, rmax)
    ax2.set_ylim(zmin, zmax)
    ax2.set_title("Velocity", fontsize=AxesLabel)
    ax2.axis('off')

    # Adjust spacing between subplots
    plt.tight_layout()
    
    # Add colorbars
    l1, b1, w1, h1 = ax1.get_position().bounds
    l2, b2, w2, h2 = ax2.get_position().bounds
    
    # Temperature colorbar
    cb1_ax = fig.add_axes([l1-0.04, b1, 0.01, h1])
    c1 = plt.colorbar(cntrl1, cax=cb1_ax, orientation='vertical')
    c1.set_label(r'$T$', fontsize=TickLabel, labelpad=5)
    c1.ax.tick_params(labelsize=TickLabel)
    c1.ax.yaxis.set_ticks_position('left')
    c1.ax.yaxis.set_label_position('left')
    c1.ax.yaxis.set_major_formatter(plt.matplotlib.ticker.StrMethodFormatter('{x:,.1f}'))
    
    # Velocity colorbar
    cb2_ax = fig.add_axes([l2+w2+0.01, b2, 0.01, h2])
    c2 = plt.colorbar(cntrl2, cax=cb2_ax, orientation='vertical')
    c2.set_label(r'$\|u_i\|$', fontsize=TickLabel, labelpad=5)
    c2.ax.tick_params(labelsize=TickLabel)
    c2.ax.yaxis.set_major_formatter(plt.matplotlib.ticker.StrMethodFormatter('{x:,.1f}'))
    
    # Set figure background to white for saving
    fig.patch.set_facecolor('white')
    plt.savefig(name, bbox_inches="tight", facecolor='white')
    plt.close()



def main():
    # Get number of CPUs from command line argument, or use all available
    parser = argparse.ArgumentParser()
    parser.add_argument('--CPUs', type=int, default=mp.cpu_count(), help='Number of CPUs to use')
    parser.add_argument('--nGFS', type=int, default=100, help='Number of restart files to process')
    parser.add_argument('--GridsPerR', type=int, default=512, help='Number of grids per R')
    parser.add_argument('--ZMAX', type=float, default=1.5, help='Maximum Z value')
    parser.add_argument('--RMAX', type=float, default=1.5, help='Maximum R value')
    parser.add_argument('--ZMIN', type=float, default=-1.5, help='Minimum Z value')
    parser.add_argument('--RMIN', type=float, default=-1.5, help='Minimum R value')
    parser.add_argument('--tsnap', type=float, default=0.1, help='Time snap')
    parser.add_argument('--caseToProcess', type=str, default='../testCases/2-2D_rot_cylinder_ink', help='Case to process')  
    parser.add_argument('--folderToSave', type=str, default='2-2D_rot_cylinder_ink', help='Folder to save')
    args = parser.parse_args()

    num_processes = args.CPUs
    nGFS = args.nGFS
    tsnap = args.tsnap
    ZMAX = args.ZMAX
    RMAX = args.RMAX
    ZMIN = args.ZMIN
    RMIN = args.RMIN
    rmin, rmax, zmin, zmax = [RMIN, RMAX, ZMIN, ZMAX]
    GridsPerR = args.GridsPerR

    lw = 2
    folder = args.folderToSave
    caseToProcess = args.caseToProcess

    if not os.path.isdir(folder):
        os.makedirs(folder)

    # Create a pool of worker processes
    with mp.Pool(processes=num_processes) as pool:
        # Create partial function with fixed arguments
        process_func = partial(process_timestep, caseToProcess=caseToProcess,
                             folder=folder, tsnap=tsnap,
                             GridsPerR=GridsPerR, rmin=rmin, rmax=rmax, 
                             zmin=zmin, zmax=zmax, lw=lw)
        # Map the process_func to all timesteps
        pool.map(process_func, range(nGFS))

if __name__ == "__main__":
    main()