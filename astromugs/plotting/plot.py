#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
file name: plot
@author: Sacha Gavino
last update: Jan 22
language: PYTHON 3
__________________________________________________________________________________________
short description:  plotting of the disk thermal model
__________________________________________________________________________________________
"""
import glob, sys

import os
import numpy as np
import pandas as pd


import re
from matplotlib.collections import PolyCollection


import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astromugs.constants.constants import autocm


def density2D_grid(path='thermal/', vmin=1e-30, vmax=1e-15, cmap='gnuplot2', dens_type='mass',
                    xlim=None, ylim=None, dust=None, figsize=(10, 14)):
    """Plot all dust species on a single figure with subplots, plus total density."""
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[grid.columns[0]].values

    dens = pd.read_table(path + 'dust_density.inp', engine='python', header=None, skiprows=3)
    dens = dens[0].values
    nspecies = int(len(dens) / (nr * nt))
    dens = np.reshape(dens, (nspecies, nt, nr))

    grid = np.array(grid, copy=True)

    r_edge = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    rr_edge, tt_edge = np.meshgrid(r_edge, theta_edge)
    R = rr_edge * np.sin(tt_edge)
    Z = rr_edge * np.cos(tt_edge)
    dens = np.array(dens,copy=True)
    dens[dens <= 1e-100] = 1e-100

    # Try to read grain sizes for subplot labels
    import os
    sizes_file = path + 'dust_sizes.inp'
    if os.path.isfile(sizes_file):
        sizes = np.loadtxt(sizes_file)
        sizes = np.atleast_1d(sizes)
    elif dust != None:
        rho_m = dust.rho_m #g.cm3
        sizes = dust.sizes()[0] # microns
        grain_mass = dust.grainmass() # in gram

    else:
        sizes = None

    # Layout: enough panels for all species + 1 total
    npanels = nspecies + 1
    ncols = min(nspecies, 4)
    nrows = int(np.ceil(npanels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])

    for idx in range(nrows * ncols):
        ax = axes.flat[idx]
        if idx < nspecies:
            if dens_type == 'number':
                im = ax.pcolormesh(R, Z, dens[idx]/grain_mass[idx], cmap=cmap, shading='auto',
                                   norm=LogNorm(vmin=vmin, vmax=vmax))
                fig.colorbar(im, cax=cbar_ax, label=r'n$_\mathrm{d}$ [cm$^{-3}$]')
            if dens_type == 'surface':
                im = ax.pcolormesh(R, Z, 4*np.pi*sizes[idx]*1e-4*dens[idx]/grain_mass[idx], cmap=cmap, shading='auto',
                                   norm=LogNorm(vmin=vmin, vmax=vmax))
                fig.colorbar(im, cax=cbar_ax, label=r'surfaces [cm$^{-1}$]')
            elif dens_type == 'mass':
                im = ax.pcolormesh(R, Z, dens[idx], cmap=cmap, shading='auto',
                                   norm=LogNorm(vmin=vmin, vmax=vmax))
                fig.colorbar(im, cax=cbar_ax, label=r'$\rho_\mathrm{d}$ [g cm$^{-3}$]')
            ax.set_title(f'bin {idx+1}', fontsize=12)
            if sizes is not None and idx < len(sizes):
                s = sizes[idx]
                if s >= 1e3:
                    size_label = f'{s/1e3:.1f} mm'
                else:
                    size_label = f'{s:.2f} ' + r'$\mu$m'
                ax.text(0.05, 0.95, size_label, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top',
                        horizontalalignment='left', bbox=props)
        elif idx == nspecies:
            total = dens.sum(axis=0)
            im = ax.pcolormesh(R, Z, total, cmap=cmap, shading='auto',
                               norm=LogNorm(vmin=vmin, vmax=vmax))
            ax.set_title('total', fontsize=12)
        else:
            ax.set_visible(False)
            continue

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    # Shared labels
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('z [au]', fontsize=14)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    
    #fig.colorbar(im, cax=cbar_ax, label=r'$\rho_\mathrm{d}$ [g cm$^{-3}$]')

    plt.show()


def density1D_midplane(path='thermal/', vmin=1e-30, vmax=1e-15, dens_type='mass',
                        xlim=None, dust=None, figsize=(12, 8)):
    """
    Plots the 1D dust density profile in the midplane (z=0 / theta=pi/2) 
    as a function of radius for each dust species.

    Parameters:
    -----------
    path : str
        Path to the directory containing RADMC-3D files.
    vmin, vmax : float
        Limits for the Y-axis (density).
    dens_type : str
        Type of density to plot: 'mass' (g/cm^3), 'number' (cm^-3), or 'surface' (cm^-1).
    xlim : tuple/list, optional
        Limits for the X-axis (radius in au).
    dust : object, optional
        An external dust object containing grain sizes and masses if files are missing.
    figsize : tuple
        Size of the output matplotlib figure.
    """
    
    # 1. Read grid structure and dust density data
    # Read the AMR grid file to extract dimensions (nr = radial bins, nt = theta bins)
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = np.array(grid[grid.columns[0]].values, copy=True)
    
    # Read the raw dust density file (flat 1D array of values)
    dens = pd.read_table(path + 'dust_density.inp', engine='python', header=None, skiprows=3)
    dens = dens[0].values
    
    # Deduce the number of dust species and reshape into a 3D array: (species, theta, radius)
    nspecies = int(len(dens) / (nr * nt))
    dens = np.reshape(dens, (nspecies, nt, nr))

    # 2. Extract radial coordinates at cell centers (convert from cm to au)
    # autocm is assumed to be a globally defined constant (1 au = 1.496e13 cm)
    r_edge = grid[:nr+1] / autocm
    r_center = 0.5 * (r_edge[:-1] + r_edge[1:])

    # 3. Identify the midplane index (theta = pi/2)
    # In RADMC-3D spherical coordinates, the equator is exactly at the midpoint of the theta axis
    idx_midplane = nt // 2 

    # 4. Read grain sizes and masses for plotting labels and conversions
    sizes_file = path + 'dust_sizes.inp'
    if os.path.isfile(sizes_file):
        sizes = np.loadtxt(sizes_file)
        sizes = np.atleast_1d(sizes)
    elif dust != None:
        rho_m = dust.rho_m #g.cm3
        sizes = dust.sizes()[0] # microns
        grain_mass = dust.grainmass() # in gram

    else:
        sizes = None

    # 5. Configure the figure layout (Grid of subplots)
    npanels = nspecies + 1  # Number of species + 1 extra panel for the total sum
    ncols = min(npanels, 3) # Maximum of 3 columns
    nrows = int(np.ceil(npanels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes) # Ensure axes is always a 2D array even for a single row

    # Determine Y-axis label depending on the requested density type
    if dens_type == 'number':
        ylabel = r'$n_\mathrm{d}$ [cm$^{-3}$]'
    elif dens_type == 'surface':
        ylabel = r'Surfaces [cm$^{-1}$]'
    else:
        ylabel = r'$\rho_\mathrm{d}$ [g cm$^{-3}$]'

    # 6. Plotting loop over all available subplot slots
    for idx in range(nrows * ncols):
        ax = axes.flat[idx]
        
        if idx < nspecies:
            # Extract 1D radial profile at the midplane for the current dust species
            profile = dens[idx, idx_midplane, :]
            
            # Apply conversion factors based on the selected density type
            if dens_type == 'number':
                y_data = profile / grain_mass[idx] # Mass density / mass of one grain
            elif dens_type == 'surface':
                # Cross-sectional area calculation (converting size from micron to cm)
                y_data = 4 * np.pi * (sizes[idx] * 1e-4) * profile / grain_mass[idx]
            elif dens_type == 'mass':
                y_data = profile # Default is raw mass density
                
            ax.plot(r_center, y_data, color='darkblue', lw=2)
            ax.set_title(f'Bin {idx+1}', fontsize=12)
            
            # Add text box indicating the grain size for this specific bin
            if sizes is not None and idx < len(sizes):
                s = sizes[idx]
                size_label = f'{s/1e3:.1f} mm' if s >= 1e3 else f'{s:.2f} ' + r'$\mu$m'
                ax.text(0.05, 0.95, size_label, transform=ax.transAxes,
                        fontsize=12, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        elif idx == nspecies:
            # Plot total cumulative density (only relevant/calculated for mass density)
            if dens_type == 'mass':
                total_profile = dens[:, idx_midplane, :].sum(axis=0)
                ax.plot(r_center, total_profile, color='black', lw=2.5, linestyle='--')
                ax.set_title('Total Mass', fontsize=12)
            else:
                ax.axis('off') # Hide total panel if it's not mass density
        else:
            ax.axis('off') # Hide any remaining empty subplots in the grid

        # Configure axes scales and limits
        ax.set_yscale('log')
        ax.set_ylim(vmin, vmax)
        if xlim:
            ax.set_xlim(xlim)
        else:
            ax.set_xscale('log') # Logarithmic scale is standard for protoplanetary disks

    # Add global outer axis labels (only on edge plots thanks to sharex/sharey)
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=12)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel, fontsize=12)

    fig.tight_layout()
    plt.show()


def density2D_grid_interactive(path='thermal/', vmin=1e-30, vmax=1e-15, cmap='gnuplot2', dens_type='mass',
                                xlim=None, ylim=None, dust=None, figsize=(10, 14)):
    """Interactive version of density2D_grid with sliders for vmin/vmax.
    Requires %matplotlib widget in the notebook."""
    import ipywidgets as widgets
    from IPython.display import display
    import os

    # --- Load data once ---
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[grid.columns[0]].values

    dens = pd.read_table(path + 'dust_density.inp', engine='python', header=None, skiprows=3)
    dens = dens[0].values
    nspecies = int(len(dens) / (nr * nt))
    dens = np.reshape(dens, (nspecies, nt, nr))

    r_edge = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    rr_edge, tt_edge = np.meshgrid(r_edge, theta_edge)
    R = rr_edge * np.sin(tt_edge)
    Z = rr_edge * np.cos(tt_edge)

    dens[dens <= 1e-100] = 1e-100

    sizes_file = path + 'dust_sizes.inp'
    if os.path.isfile(sizes_file):
        sizes = np.loadtxt(sizes_file)
        sizes = np.atleast_1d(sizes)
    elif dust is not None:
        sizes = dust.sizes()[0]
        grain_mass = dust.grainmass()
    else:
        sizes = None

    # Precompute plot data for each panel
    plot_data = []
    for idx in range(nspecies):
        if dens_type == 'number' and dust is not None:
            plot_data.append(4 * np.pi * sizes[idx] * 1e-4 * dens[idx] / grain_mass[idx])
        else:
            plot_data.append(dens[idx])
    plot_data.append(dens.sum(axis=0))  # total

    npanels = nspecies + 1
    ncols = min(nspecies, 4)
    nrows = int(np.ceil(npanels / ncols))

    # --- Build figure ---
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    meshes = []
    for idx in range(nrows * ncols):
        ax = axes.flat[idx]
        if idx < nspecies:
            im = ax.pcolormesh(R, Z, plot_data[idx], cmap=cmap, shading='auto',
                               norm=LogNorm(vmin=vmin, vmax=vmax))
            meshes.append(im)
            ax.set_title(f'bin {idx+1}', fontsize=12)
            if sizes is not None and idx < len(sizes):
                s = sizes[idx]
                if s >= 1e3:
                    size_label = f'{s/1e3:.1f} mm'
                else:
                    size_label = f'{s:.2f} ' + r'$\mu$m'
                ax.text(0.05, 0.95, size_label, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top',
                        horizontalalignment='left', bbox=props)
        elif idx == nspecies:
            im = ax.pcolormesh(R, Z, plot_data[nspecies], cmap=cmap, shading='auto',
                               norm=LogNorm(vmin=vmin, vmax=vmax))
            meshes.append(im)
            ax.set_title('total', fontsize=12)
        else:
            ax.set_visible(False)
            continue
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('z [au]', fontsize=14)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(meshes[-1], cax=cbar_ax, label=r'$\rho_\mathrm{d}$ [g cm$^{-3}$]')

    # --- Sliders ---
    log_vmin = np.log10(vmin)
    log_vmax = np.log10(vmax)

    vmin_slider = widgets.FloatSlider(
        value=log_vmin, min=-50, max=0, step=0.5,
        description='log(vmin)', continuous_update=False,
        style={'description_width': 'initial'}, layout=widgets.Layout(width='400px'))
    vmax_slider = widgets.FloatSlider(
        value=log_vmax, min=-50, max=0, step=0.5,
        description='log(vmax)', continuous_update=False,
        style={'description_width': 'initial'}, layout=widgets.Layout(width='400px'))

    def update_clim(change):
        new_vmin = 10**vmin_slider.value
        new_vmax = 10**vmax_slider.value
        if new_vmin >= new_vmax:
            return
        new_norm = LogNorm(vmin=new_vmin, vmax=new_vmax)
        for m in meshes:
            m.set_norm(new_norm)
        cbar.update_normal(meshes[-1])
        fig.canvas.draw_idle()

    vmin_slider.observe(update_clim, names='value')
    vmax_slider.observe(update_clim, names='value')

    display(widgets.HBox([vmin_slider, vmax_slider]))
    plt.show()




def temperature2D_grid(path='thermal/', vmin=1e0, vmax=1e3, cmap='gnuplot2',
                    xlim=None, ylim=None, figsize=(10, 14)):
    """Plot all dust species on a single figure with subplots, plus total density."""
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[grid.columns[0]].values

    temp = pd.read_table(path + 'dust_temperature.dat', engine='python', header=None, skiprows=3)
    temp = temp[0].values
    nspecies = int(len(temp) / (nr * nt))
    temp = np.reshape(temp, (nspecies, nt, nr))
    grid = np.array(grid,copy=True)
    r_edge = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    rr_edge, tt_edge = np.meshgrid(r_edge, theta_edge)
    R = rr_edge * np.sin(tt_edge)
    Z = rr_edge * np.cos(tt_edge)

    # Convert edge grid to cell-center grid for contour plotting
    R_center = 0.25 * (R[:-1, :-1] + R[1:, :-1] + R[:-1, 1:] + R[1:, 1:])
    Z_center = 0.25 * (Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])


        # Try to read grain sizes for subplot labels
    import os
    sizes_file = path + 'dust_sizes.inp'
    if os.path.isfile(sizes_file):
        sizes = np.loadtxt(sizes_file)
        sizes = np.atleast_1d(sizes)
    else:
        sizes = None

    # Layout: enough panels for all species + 1 total
    npanels = nspecies + 1
    ncols = min(nspecies, 4)
    nrows = int(np.ceil(npanels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    levels = np.arange(vmin, vmax + 10, 10)

    for idx in range(nrows * ncols):
        ax = axes.flat[idx]
        if idx < nspecies:
            im = ax.pcolormesh(R, Z, temp[idx], cmap=cmap, shading='auto',
                              norm=LogNorm(vmin=vmin, vmax=vmax))
            #im = ax.contourf(R_center, Z_center, temp[idx], levels=levels, cmap=cmap)
            cs = ax.contour(
                R_center,
                Z_center,
                temp[idx],
                levels=[20],
                colors='black',
                linewidths=2
            )
            #ax.clabel(cs, fmt="20 K", fontsize=10)
            ax.set_title(f'bin {idx+1}', fontsize=12)
            if sizes is not None and idx < len(sizes):
                s = sizes[idx]
                if s >= 1e3:
                    size_label = f'{s/1e3:.1f} mm'
                else:
                    size_label = f'{s:.2f} ' + r'$\mu$m'
                ax.text(0.05, 0.95, size_label, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top',
                        horizontalalignment='left', bbox=props)
        elif idx == nspecies:
            total = temp.sum(axis=0)
            im = ax.pcolormesh(R, Z, total, cmap=cmap, shading='auto',
                               norm=LogNorm(vmin=vmin, vmax=vmax))
            ax.set_title('total', fontsize=12)
        else:
            ax.set_visible(False)
            continue

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    # Shared labels
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('z [au]', fontsize=14)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label=r'T [K]')

    plt.show()


def midplane_temp(path='thermal/', xlim=None, ylim=None):
    grid = pd.read_table(path+'amr_grid.inp', engine='python', skiprows=5)
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    try:
        temp = pd.read_table(path+'dust_temperature.dat', engine='python', header=None, skiprows=3)
    except IOError:
        print('plot.midplane_temp: the file dust_temperature.dat is not present. Run a dust thermal simulation first.')
        sys.exit(1)
    temp = temp[0].values
    nbspecies = int(len(temp)/(nr*nt))
    temp = np.reshape(temp, (nbspecies, nt, nr))
    grid = np.array(grid,copy=True)
    dist = grid[:nr+1]/autocm
    theta = grid[nr+1:nr+1+nt+1]
    theta[-1] = np.pi
    dist, tt = np.meshgrid(dist, theta)
    rr = dist*np.sin(tt)
    zz = dist*np.cos(tt)
    midtemp = temp[:, 90, :]
    radii = 0.5*(rr[90][0:rr[90].size-1] + rr[90][1:rr[90].size])

    #--PLOT FIGURE--
    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111)
    #-----profiles
    midtemp = pd.DataFrame(data=midtemp.transpose())
    for ispec in range(0, nbspecies):
        ax.plot(radii, midtemp[ispec].rolling(window=6, center=True).mean(), linewidth=2, linestyle='-', label='bin: {}'.format(ispec+1))
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
    ax.set_xlabel(r'r [au]', fontsize = 20)
    ax.set_ylabel(r'T [K]', fontsize = 20)
    ax.legend(fontsize=15)
    ax.tick_params(labelsize=18)
    plt.show()

def vertical_temp(thermpath='thermal/', chempath='chemistry/', r=100):
    grid = pd.read_table(thermpath+'amr_grid.inp', engine='python', skiprows=5)
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    temp = pd.read_table(thermpath+'dust_temperature.dat', engine='python', header=None, skiprows=3)
    temp = temp[0].values
    nbspecies = int(len(temp)/(nr*nt))

    if nbspecies == 1:
        try:
            temp = pd.read_table(chempath+str(r)+'AU/1D_static.dat', sep=r"\s+", engine='python', header=None, comment='!')
        except IOError:
            print('plot.vertical_temp: radius {} does not exit in the model or path is not correct.'.format(r))
            sys.exit(1)
        #--PLOT FIGURE--
        fig = plt.figure(figsize=(9.6, 8.2))
        ax = fig.add_subplot(111)
        #-----profiles
        ax.plot(temp[5], temp[0], linewidth=2, linestyle='-', label='{} AU'.format(r))
        # ax.set_ylim(0,60)
        # ax.set_xlim(1,350)
        ax.set_ylabel(r'z [au]', fontsize = 20)
        ax.set_xlabel(r'T$_\mathrm{d}$ [K]', fontsize = 20)
        ax.legend(fontsize=15)
        ax.tick_params(labelsize=18)
        plt.show()
    elif nbspecies > 1:
        try:
            static = pd.read_table(chempath+str(r)+'AU/1D_static.dat', sep=r"\s+", engine='python', header=None, comment='!')
            temp = pd.read_table(chempath+str(r)+'AU/temperatures.dat', sep=r"\s+", engine='python', header=None)
        except IOError:
            print('plot.vertical_temp: radius = {} au does not exit in the model or path is not correct.'.format(r))
            sys.exit(1)
        #--PLOT FIGURE--
        fig = plt.figure(figsize=(9.6, 8.2))
        ax = fig.add_subplot(111)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.91, 0.05, '{} AU'.format(r), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=16, bbox=props)
        #-----profiles
        for ai in range(nbspecies):
            ax.plot(temp[ai], static[0], linewidth=2, linestyle='-', label='bin: {}'.format(ai+1))
        # ax.set_ylim(0,60)
        # ax.set_xlim(1,350)
        ax.set_ylabel(r'z [au]', fontsize = 20)
        ax.set_xlabel(r'T$_\mathrm{d}$ [K]', fontsize = 20)
        ax.legend(fontsize=15)
        ax.tick_params(labelsize=18)
        plt.show()

def avz(chempath='thermal/', r=100):
    static = pd.read_table(chempath+str(r)+'AU/1D_static.dat', sep=r"\s+", engine='python', header=None, comment='!', skiprows=1)
    #--PLOT FIGURE--
    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.91, 0.05, '{} AU'.format(r), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=16, bbox=props)
    #-----profiles
    ax.plot(static[3], static[0], linewidth=2, linestyle='-', label='vertical Av')
    # ax.set_ylim(0,60)
    # ax.set_xlim(1,350)
    ax.set_xlabel(r'z [au]', fontsize = 20)
    ax.set_ylabel(r'A$_\mathrm{\nu}$ [mag]', fontsize = 20)
    ax.legend(fontsize=15)
    ax.tick_params(labelsize=18)
    plt.show()

def opacity(path='thermal/'):
    opaclist = sorted(glob.glob(path+'dustkap*'))

    #---absorption
    fig = plt.figure(figsize=(9.6, 8.2)) #fig = plt.figure(figsize=(9.6, 7.2))
    ax = fig.add_subplot(111) 
    ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize=18)
    ax.set_ylabel(r'$\kappa_\mathrm{abs}$ [cm$^2$/g]', fontsize=18)
    ax.set_xlim(1e-1,1e4)
    ax.set_ylim(1e-2,1e5)
    for opac in opaclist:
        name = opac.split("_")[1].split(".")[0]
        kappa = pd.read_table(opac, sep=r"\s+", comment='#', header=None, skiprows=10)
        ax.loglog(kappa[0], kappa[1], linewidth=2, label=name)
    ax.tick_params(labelsize=22)
    ax.legend(fontsize=15)
    plt.show()

    #---scattering
    fig = plt.figure(figsize=(9.6, 8.2)) #fig = plt.figure(figsize=(9.6, 7.2))
    ax = fig.add_subplot(111) 
    ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize=18)
    ax.set_ylabel(r'$\kappa_\mathrm{scat}$ [cm$^2$/g]', fontsize=18)
    ax.set_xlim(1e-1,1e4)
    ax.set_ylim(1e-2,1e5)
    for opac in opaclist:
        name = opac.split("_")[1].split(".")[0]
        kappa = pd.read_table(opac, sep=r"\s+", comment='#', header=None, skiprows=10)
        ax.loglog(kappa[0], kappa[2], linewidth=2, label=name)
    ax.tick_params(labelsize=22)
    ax.legend(fontsize=15)
    plt.show()

    #---angles
    fig = plt.figure(figsize=(9.6, 8.2)) #fig = plt.figure(figsize=(9.6, 7.2))
    ax = fig.add_subplot(111) 
    ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize=18)
    ax.set_ylabel(r'<cos($\theta$)>', fontsize=18)
    ax.set_xlim(1e-1,1e4)
    #ax.set_ylim(0,1)
    for opac in opaclist:
        name = opac.split("_")[1].split(".")[0]
        kappa = pd.read_table(opac, sep=r"\s+", comment='#', header=None, skiprows=10)
        ax.loglog(kappa[0], kappa[3], linewidth=2, label=name)
    ax.tick_params(labelsize=22)
    ax.legend(fontsize=15)
    plt.show()


def localflux(path='thermal/'):
    #---1/ Get grid shape and reshape the local flux array accordingly
    flux = pd.read_table(path+'mean_intensity.out', sep=r"\s+", comment='#', header=None, skiprows=4)
    grid = pd.read_table(path+'amr_grid.inp', engine='python', skiprows=5)
    lam = pd.read_table(path+'mcmono_wavelength_micron.inp', engine='python', header=None, skiprows=1)
    lam = lam[0].values
    flux = flux[0].values
    
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    dist = grid[:nr+1]/autocm
    theta = grid[nr+1:nr+1+nt+1]
    theta[-1] = np.pi
    dist, tt = np.meshgrid(dist, theta)
    rr = dist*np.sin(tt)
    radii = 0.5*(rr[90][0:rr[90].size-1] + rr[90][1:rr[90].size])
    zz = dist*np.cos(tt)
    nlam = int(len(flux)/(nr*nt))
    flux = np.reshape(flux, (nlam, nt, nr))
    midflux = flux[:, 90, :]

    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111)
    #-----profiles
    midflux_df = pd.DataFrame(data=midflux.transpose())
    for ilam in range(0, nlam, 2):
        ax.semilogy(radii, midflux_df[ilam].rolling(window=5, center=True).mean(), linewidth=1, linestyle='-')
    ax.set_xlim(1,200)
    ax.set_ylim(1e-30,1e-10)
    ax.set_xlabel(r'r [au]', fontsize = 20)
    ax.set_ylabel(r'Flux', fontsize = 20)
    ax.grid()
    ax.tick_params(labelsize=22)
    plt.show()


def image(pathfile='thermal/', distance=100, vmin=1e-10, vmax=1e3, cmap='gnuplot2',
          labels=None):
    """Plot a multi-wavelength RADMC3D continuum image (Stokes I).

    Parameters
    ----------
    pathfile : str
        Path to the RADMC3D image.out file.
    distance : float, optional
        Source distance in parsec. Used to convert specific intensity
        [erg/s/cm²/Hz/sr] to flux density [Jy/pixel]. Default is 100 pc.
    vmin, vmax : float, optional
        Color scale limits in Jy/pixel.
    cmap : str, optional
        Colormap name.
    labels : list of str, optional
        Wavelength labels for each panel. If None, labels are generated
        automatically from the wavelengths read in the image header.
    """
    # --- Read RADMC3D image header ---
    with open(pathfile, 'r') as f:
        iformat = int(f.readline())                                 # 1 = I only, 3 = full Stokes
        npix_x, npix_y = [int(x) for x in f.readline().split()]   # image size [pixels]
        nlam = int(f.readline())                                    # number of wavelengths
        pix_cm, _ = [float(x) for x in f.readline().split()]       # pixel size [cm]
        wavelengths = [float(f.readline()) for _ in range(nlam)]   # wavelengths [microns]

    # iformat 3 → full Stokes (I Q U V); anything else → intensity only
    nstokes = 4 if iformat == 3 else 1

    pix_au   = pix_cm / autocm
    box_au   = npix_x * pix_au
    half_box = box_au / 2.0

    # --- Pixel solid angle and Jy/pixel conversion factor ---
    distance_cm = distance * 3.086e18               # pc → cm
    omega_pix   = (pix_cm / distance_cm) ** 2       # sr/pixel
    to_jy       = 1e23 * omega_pix                  # erg/s/cm²/Hz/sr → Jy/pixel

    data = np.loadtxt(pathfile, skiprows=4 + nlam + 1)
    data = np.reshape(data, (nlam, npix_y, npix_x, nstokes))

    extent = [-half_box, half_box, -half_box, half_box]

    # --- Stokes I images in Jy/pixel ---
    imgs = [data[i, :, :, 0] * to_jy for i in range(nlam)]

    # --- Wavelength labels ---
    if labels is None:
        def _wave_label(lam):
            if lam >= 1000:
                return f'Stokes I - {lam/1000:.2g} mm'
            elif lam >= 1:
                return f'Stokes I - {lam:.4g} μm'
            else:
                return f'Stokes I - {lam*1000:.4g} nm'
        labels = [_wave_label(w) for w in wavelengths]

    fig, axes = plt.subplots(1, nlam, figsize=(5*nlam, 5), sharex=True, sharey=True)
    if nlam == 1:
        axes = [axes]

    for i, ax in enumerate(axes):

        im = ax.imshow(
            imgs[i],
            origin='lower',
            extent=extent,
            cmap=cmap,
            norm=LogNorm(vmin=vmin, vmax=vmax),
            interpolation='nearest'
        )

        ax.tick_params(labelsize=17)
        ax.set_xlabel(r'x [au]', fontsize=17)

        if i == 0:
            ax.set_ylabel(r'y [au]', fontsize=17)

        ax.text(
            0.05, 0.95,
            labels[i],
            horizontalalignment='left',
            verticalalignment='top',
            color='red',
            transform=ax.transAxes,
            fontsize=16,
            fontweight='bold'
        )

    # --- Proper colorbar placement ---
    cbar = fig.colorbar(
        im,
        ax=axes,
        location='right',
        fraction=0.025,
        pad=0.02
    )

    cbar.set_label(r'$I_\nu$ [Jy pixel$^{-1}$]', fontsize=17)
    cbar.ax.tick_params(labelsize=14)

    plt.show()


def image_vertical_cut(pathfile='thermal/', distance=100, xlim=None, ylim=None,
                       labels=None, figsize=(9.6, 8.2)):
    """Plot vertical cuts (along y at x=0) of Stokes I for each wavelength.

    Parameters
    ----------
    pathfile : str
        Path to the RADMC3D image.out file.
    distance : float, optional
        Source distance in parsec. Used to convert specific intensity
        [erg/s/cm²/Hz/sr] to flux density [Jy/pixel]. Default is 100 pc.
    xlim, ylim : tuple, optional
        Axis limits.
    labels : list of str, optional
        Wavelength labels. If None, generated automatically from the header.
    figsize : tuple, optional
        Figure size.
    """
    # --- Read RADMC3D image header ---
    with open(pathfile, 'r') as f:
        iformat = int(f.readline())                                 # 1 = I only, 3 = full Stokes
        npix_x, npix_y = [int(x) for x in f.readline().split()]   # image size [pixels]
        nlam = int(f.readline())                                    # number of wavelengths
        pix_cm, _ = [float(x) for x in f.readline().split()]       # pixel size [cm]
        wavelengths = [float(f.readline()) for _ in range(nlam)]   # wavelengths [microns]

    nstokes = 4 if iformat == 3 else 1

    pix_au   = pix_cm / autocm
    box_au   = npix_y * pix_au
    half_box = box_au / 2.0

    # --- Pixel solid angle and Jy/pixel conversion factor ---
    distance_cm = distance * 3.086e18
    omega_pix   = (pix_cm / distance_cm) ** 2
    to_jy       = 1e23 * omega_pix

    data = np.loadtxt(pathfile, skiprows=4 + nlam + 1)
    data = np.reshape(data, (nlam, npix_y, npix_x, nstokes))

    y_au = np.linspace(-half_box, half_box, npix_y)
    ix0  = npix_x // 2  # column at x=0

    if labels is None:
        def _wave_label(lam):
            if lam >= 1000:
                return f'{lam/1000:.2g} mm'
            elif lam >= 1:
                return f'{lam:.4g} μm'
            else:
                return f'{lam*1000:.4g} nm'
        labels = [_wave_label(w) for w in wavelengths]

    fig, ax = plt.subplots(figsize=figsize)

    for i in range(nlam):
        flux_cut = data[i, :, ix0, 0] * to_jy
        ax.semilogy(y_au, flux_cut, linewidth=2, label=labels[i])

    ax.set_xlabel(r'y [au]', fontsize=20)
    ax.set_ylabel(r'$I_\nu$ [Jy pixel$^{-1}$]', fontsize=20)
    ax.legend(fontsize=15)
    ax.tick_params(labelsize=18)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    plt.show()

def numberdens(species='CO', path='thermal/', vmin=1e0, vmax=1e8, cmap='gnuplot2',
               ncols=3, xlim=None, ylim=None, figsize=None,
               save=False, savename='numberdens.pdf'):
    """Plot 2D maps of molecular number densities from ``numberdens_XXX.inp`` files.

    Accepts a single species name or a list of names. Multiple species are
    displayed as a mosaic of subplots sharing the same colour scale.

    Parameters
    ----------
    species : str or list of str, optional
        Species name(s) matching ``numberdens_<species>.inp`` files.
        Default is ``'CO'``.
    path : str, optional
        Path to the directory containing the RADMC-3D files. Default is
        ``'thermal/'``.
    vmin, vmax : float, optional
        Shared colour scale limits [cm :sup:`-3`]. Default is ``1e0``
        and ``1e8``.
    cmap : str, optional
        Colormap name. Default is ``'gnuplot2'``.
    ncols : int, optional
        Maximum number of columns in the mosaic. Default is ``3``.
    xlim, ylim : tuple, optional
        Axis limits (R, Z) in AU applied to every panel.
    figsize : tuple or None, optional
        Figure size. If None, computed automatically from ``ncols`` and
        the number of rows.
    save : bool, optional
        Save the figure to ``savename``. Default is False.
    savename : str, optional
        Output filename when ``save=True``. Default is ``'numberdens.pdf'``.
    """
    species_list = [species] if isinstance(species, str) else list(species)
    nspecies = len(species_list)

    # --- Grid: cell centres (shape nt × nr) — required for shading='gouraud' ---
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = np.array(grid[grid.columns[0]].values, copy=True)

    r_edge     = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    r_cen     = 0.5 * (r_edge[:-1]     + r_edge[1:])
    theta_cen = 0.5 * (theta_edge[:-1] + theta_edge[1:])
    rr, tt = np.meshgrid(r_cen, theta_cen)          # (nt, nr)
    R = rr * np.sin(tt)
    Z = rr * np.cos(tt)

    # --- Layout ---
    ncols = min(ncols, nspecies)
    nrows = int(np.ceil(nspecies / ncols))
    if figsize is None:
        figsize = (5 * ncols + 1, 4 * nrows)

    norm = LogNorm(vmin=vmin, vmax=vmax)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()

    im = None
    for idx, sp in enumerate(species_list):
        ax = axes[idx]
        nd = pd.read_table(path + f'numberdens_{sp}.inp',
                           engine='python', header=None, skiprows=2)
        nd = nd[0].values.reshape(nt, nr)
        nd = np.where(nd <= 0, 1e-100, nd)

        im = ax.pcolormesh(R, Z, nd, cmap=cmap, shading='gouraud',
                           norm=norm, rasterized=True)
        ax.set_title(sp, fontsize=13)
        ax.tick_params(labelsize=12)
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    # Hide unused panels
    for idx in range(nspecies, len(axes)):
        axes[idx].set_visible(False)

    # Shared axis labels
    for ax in axes[(nrows - 1) * ncols:]:
        ax.set_xlabel('R [au]', fontsize=13)
    for i in range(nrows):
        axes[i * ncols].set_ylabel('Z [au]', fontsize=13)

    # Single shared colorbar on the right
    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label=r'$n$ [cm$^{-3}$]')

    if save:
        fig.savefig(savename, bbox_inches='tight')

    plt.show()


def static(chempath='chemistry/', column='nH', vmin=1, vmax=50, iso=None, cmap='gnuplot2',
           xlim=None, ylim=None, figsize=(6, 6), save=False, savename='filename.pdf'):
    """Plot a 2D map of a column from the 1D_static.dat files.

    Scans all XXAU/ folders in chempath, reads each 1D_static.dat,
    and builds a 2D (r, z) map of the chosen column.

    Parameters
    ----------
    chempath : str
        Path to the chemistry directory containing the XXAU/ folders.
    column : str
        Column to plot. One of: 'z', 'nH', 'Tg', 'Av', 'diff', 'Td',
        'inv_ab', 'conv_factor', 'a', 'uv'.
    vmin, vmax : float
        Color scale limits.
    iso : float or list, optional
        Draw contour lines of Td at these levels (e.g. 20 or [20, 50]).
    cmap : str
        Colormap name.
    xlim, ylim : tuple, optional
        Axis limits (r, z) in AU.
    figsize : tuple
        Figure size.
    save : bool
        Save the figure to savename.
    savename : str
        Output filename if save is True.
    """
    import os, re

    columns = ['z', 'nH', 'Tg', 'Av', 'diff', 'Td', 'inv_ab', 'conv_factor', 'a', 'uv']

    # Discover XXAU/ folders and extract radii
    folders = [d for d in os.listdir(chempath)
               if os.path.isdir(os.path.join(chempath, d)) and re.match(r'^\d+AU$', d)]
    rchem = sorted([int(d.replace('AU', '')) for d in folders])

    # Read all files at once (nbz may differ per radius after surface truncation)
    all_data = []
    for r in rchem:
        filepath = os.path.join(chempath, f'{r}AU', '1D_static.dat')
        df = pd.read_table(filepath, sep=r"\s+", comment='!', header=None, engine='python')
        df.columns = columns
        all_data.append(df)

    nbz_max = max(len(d) for d in all_data)

    # Build 2D arrays; NaN for cells above the truncation height of each radius
    static_map = np.full((nbz_max, len(rchem)), np.nan)
    temp_map   = np.full((nbz_max, len(rchem)), np.nan)
    zz         = np.zeros((nbz_max, len(rchem)))

    for idx, data in enumerate(all_data):
        nbz_r = len(data)
        start = nbz_max - nbz_r       # top rows belong to the truncated surface
        static_map[start:, idx] = data[column].values
        temp_map[start:, idx]   = data['Td'].values
        z_col = data['z'].values
        zz[start:, idx] = z_col
        # Extrapolate z above the truncation so pcolormesh has valid coordinates
        # (those cells are NaN in data so they will appear transparent)
        if start > 0:
            dz = (z_col[0] - z_col[1]) if nbz_r > 1 else z_col[0] * 0.1
            zz[:start, idx] = z_col[0] + np.arange(start, 0, -1) * dz

    rr, _ = np.meshgrid(rchem, np.arange(nbz_max))

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal', adjustable='box')
    t = ax.pcolormesh(rr, zz, static_map, cmap=cmap, shading='auto',
                      norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
    clr = fig.colorbar(t, pad=0.01)
    clr.set_label(column, fontsize=16)
    clr.ax.tick_params(labelsize=14)

    if iso is not None:
        if not isinstance(iso, (list, tuple)):
            iso = [iso]
        ax.contour(rr, zz, temp_map, iso, colors='black', linewidths=2.5)

    ax.set_xlabel(r'r [au]', fontsize=20)
    ax.set_ylabel(r'z [au]', fontsize=20)
    ax.tick_params(labelsize=15)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    if save:
        fig.savefig(savename, bbox_inches='tight')

    plt.show()


def nmgc_grainsizes(chempath='chemistry/', quantity='Td', vmin=None, vmax=None, cmap='gnuplot2',
                xlim=None, ylim=None, figsize=(14, 10), save=False, savename='grain_sizes.pdf'):
    """Plot 2D (r, z) maps per grain size from the 1D_grain_sizes.in files.

    Each subplot corresponds to one grain size bin.

    Parameters
    ----------
    chempath : str
        Path to the chemistry directory containing the XXAU/ folders.
    quantity : str
        What to plot: 'Td' for dust temperature, 'nd' for dust number density (nH / inv_ab).
    vmin, vmax : float, optional
        Color scale limits. Auto-determined if None.
    cmap : str
        Colormap name.
    xlim, ylim : tuple, optional
        Axis limits (r, z) in AU.
    figsize : tuple
        Figure size.
    save : bool
        Save the figure to savename.
    savename : str
        Output filename if save is True.
    """
    import os, re

    static_columns = ['z', 'nH', 'Tg', 'Av', 'diff', 'Td', 'inv_ab', 'conv_factor', 'a', 'uv']

    # Discover XXAU/ folders and extract radii
    folders = [d for d in os.listdir(chempath)
               if os.path.isdir(os.path.join(chempath, d)) and re.match(r'^\d+AU$', d)]
    rchem = sorted([int(d.replace('AU', '')) for d in folders])

    # Read the first grain_sizes file to determine N (number of grain sizes)
    # and grain radii from the first data line
    first_gs = os.path.join(chempath, f'{rchem[0]}AU', '1D_grain_sizes.in')
    with open(first_gs, 'r') as f:
        for line in f:
            stripped = line.split('!')[0].strip()
            if not stripped:
                continue
            vals = stripped.split()
            ncols = len(vals)
            grain_radii_cm = np.array([float(v) for v in vals[:ncols // 4]])
            break
    # ncols = 4*N (sizes, inv_ab, Td, CR-peak)
    ngrains = ncols // 4
    grain_radii_um = grain_radii_cm * 1e4  # cm to microns

    # Read all static files first (nbz may differ per radius after surface truncation)
    all_static = []
    for r in rchem:
        static_file = os.path.join(chempath, f'{r}AU', '1D_static.dat')
        sd = pd.read_table(static_file, sep=r"\s+", comment='!', header=None, engine='python')
        sd.columns = static_columns
        all_static.append(sd)

    nbz_max = max(len(sd) for sd in all_static)

    # Build arrays; NaN for cells above the truncation height of each radius
    data_map = np.full((ngrains, nbz_max, len(rchem)), np.nan)
    zz = np.zeros((nbz_max, len(rchem)))

    for idx, r in enumerate(rchem):
        static_data = all_static[idx]
        nbz_r = len(static_data)
        start = nbz_max - nbz_r       # top rows belong to the truncated surface
        nH    = static_data['nH'].values
        z_col = static_data['z'].values
        zz[start:, idx] = z_col
        # Extrapolate z above the truncation so pcolormesh has valid coordinates
        if start > 0:
            dz = (z_col[0] - z_col[1]) if nbz_r > 1 else z_col[0] * 0.1
            zz[:start, idx] = z_col[0] + np.arange(start, 0, -1) * dz

        # Read grain_sizes
        gs_file = os.path.join(chempath, f'{r}AU', '1D_grain_sizes.in')
        gs_lines = []
        with open(gs_file, 'r') as f:
            for line in f:
                stripped = line.split('!')[0].strip()
                if not stripped:
                    continue
                gs_lines.append([float(v) for v in stripped.split()])
        gs_array = np.array(gs_lines)  # (nbz_r, 4*ngrains)

        inv_ab = gs_array[:, ngrains:2*ngrains]       # (nbz_r, ngrains)
        Td     = gs_array[:, 2*ngrains:3*ngrains]     # (nbz_r, ngrains)

        if quantity == 'Td':
            for ig in range(ngrains):
                data_map[ig, start:, idx] = Td[:, ig]
        elif quantity == 'nd':
            for ig in range(ngrains):
                data_map[ig, start:, idx] = nH / inv_ab[:, ig]

    rr, _ = np.meshgrid(rchem, np.arange(nbz_max))

    # Layout
    ncols_plot = min(ngrains, 4)
    nrows_plot = int(np.ceil(ngrains / ncols_plot))

    fig, axes = plt.subplots(nrows_plot, ncols_plot, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    for ig in range(nrows_plot * ncols_plot):
        ax = axes.flat[ig]
        if ig < ngrains:
            ax.set_aspect('equal', adjustable='box')
            im = ax.pcolormesh(rr, zz, data_map[ig], cmap=cmap, shading='gouraud',
                               norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
            # Size label
            s = grain_radii_um[ig]
            if s >= 1e3:
                size_label = f'{s/1e3:.1f} mm'
            else:
                size_label = f'{s:.2f} ' + r'$\mu$m'
            ax.text(0.05, 0.95, size_label, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top',
                    horizontalalignment='left', bbox=props)
            if xlim:
                ax.set_xlim(xlim)
            if ylim:
                ax.set_ylim(ylim)
        else:
            ax.set_visible(False)

    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('z [au]', fontsize=14)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    if quantity == 'Td':
        label = r'T$_\mathrm{d}$ [K]'
    elif quantity == 'nd':
        label = r'n$_\mathrm{d}$ [cm$^{-3}$]'
    else:
        label = quantity
    fig.colorbar(im, cax=cbar_ax, label=label)

    if save:
        fig.savefig(savename, bbox_inches='tight')

    plt.show()



def plot_outputs_nautilus(chempath,
                          main_output_dict,
                          itime=-1,
                          MODE='chemistry',
                          key_list=['CO'],
                          fracab=True,
                          verbose=True,
                          xlim=None,
                          ylim=None,
                          colormap="plasma",
                          vmin=None,
                          vmax=None,
                          common_scale=True):
    """
    Plots a 2D vertical cross-section (poloidal cut) of Nautilus simulation outputs.

    Supports single or multiple keys (physical variables or chemical species) automatically.
    Displays a single layout or a grid of subplots based on input, with options for independent
    or globally shared colorbar scaling.

    Args:
        chempath (str/Path): Path to parent directory containing radius folders (e.g., "5AU/").
        main_output_dict (dict): Nested dictionary where keys are radii and values contain simulation data arrays.
        itime (int): Simulation timestep index to visualize. Defaults to -1 (final timestep).
        MODE (str): Type of variables to plot ('chemistry' or 'physical'). Defaults to 'chemistry'.
        key_list (str/list): Single string or list of keys (species formulas or physical variables) to plot.
        fracab (bool): If True, plots fractional abundances. If False, plots absolute number densities (cm^-3).
        verbose (bool): If True, prints missing file or size mismatch diagnostics.
        xlim (tuple): Custom (min, max) boundaries for the Radius axis.
        ylim (tuple): Custom (min, max) boundaries for the Altitude axis.
        colormap (str): Matplotlib colormap string. Defaults to "gnuplot".
        vmin (float): Forced lower bound for colorbars.
        vmax (float): Forced upper bound for colorbars.
        common_scale (bool): If True, shares identical colorbar scaling bounds across all subplots.

    Returns:
        None: Renders a Matplotlib figure window.
    """
    
    chempath = Path(chempath)

    # Convert a standalone string into a single-element list to ensure consistent iteration
    if isinstance(key_list, str):
        key_list = [key_list]

    # --- INTERNAL HELPERS ---
    def title_mol(mol_name, frac, path, verbose):
        """Formats and returns a LaTeX-compatible string for chemical species titles, including grain environments."""
        m = re.match(r"^([JK])(\d+)", mol_name)
        env = (
            f"{'surface' if m.group(1) == 'J' else 'mantle'} at grain size = {get_grain_size_in_um(Path(path)/'1D_grain_sizes.in', m.group(2), verbose=verbose)} µm"
            if m
            else "none"
        )
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        if env != "none": 
            if frac: return f"${f}$ [$n_{{{f}}}/n_H$]\n({env})"
            else: return f"${f}$ [$n_{{{f}}}$] [cm$^{{-3}}$]\n({env})"
        else:
            if frac: return f"${f}$ [$n_{{{f}}}/n_H$]"
            else: return f"${f}$ [$n_{{{f}}}$] [cm$^{{-3}}$]"

    def title_phys(variable):
        """Formats physical variable names and appends the appropriate scientific units based on string keywords."""
        name = variable.replace("_", " ").title()
        if "temperature" in name.lower(): 
            name = name + " [K]"
        elif "extinction" in name.lower(): 
            name = name + " [mag]"
        elif "density" in name.lower(): 
            name = name + " [$cm^{-3}$]"
        return name

    def get_grain_size_in_um(file_path, bin_index, verbose=False):
        """Parses 1D_grain_sizes.in to retrieve the grain bin radius mapped to micrometers."""
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('!'):
                        continue
                    if '!' in line:
                        line = line.split('!')[0].strip()
                    values = [float(val) for val in line.split()]
                    if not values:
                        continue
                    num_grains = len(values) // 4
                    radii_cm = values[:num_grains]
                    index = int(bin_index) - 1
                    if 0 <= index < num_grains:
                        return radii_cm[index] * 10000.0
            return None
        except FileNotFoundError:
            return None
    
    # --- EXTRACT DATA BY COLUMN (RADIUS) FOR ALL REQUESTED KEYS ---
    species_data = {key: [] for key in key_list}

    for r_value in main_output_dict.keys():
        folder_name = f"{r_value}AU"
        file_path = os.path.join(chempath, folder_name, "1D_static.dat")
        
        if os.path.exists(file_path):
            try:
                z_points = np.loadtxt(file_path, comments='!', usecols=0)
                sub_dict = main_output_dict[r_value]
                
                for key in key_list:
                    if MODE == 'physical':
                        full_array = sub_dict[key]
                        v_points = full_array[itime, :].copy()
                    elif MODE == 'chemistry':
                        abundance_array = sub_dict['abundances']
                        v_points = abundance_array.isel(time=itime).sel(species=key).values.copy()
                        if not fracab:
                            nH = sub_dict["H_number_density"][itime, :]
                            v_points = v_points * nH  
                    
                    # Validate column length alignment before appending
                    if len(z_points) == len(v_points):
                        species_data[key].append({
                            'R': float(r_value),  
                            'z': np.array(z_points),
                            'v': np.array(v_points)
                        })
            except Exception as e:
                if verbose: print(f"Error processing {key} for R={r_value}: {e}")
        else:
            if verbose: print(f"File not found: {file_path}")
    
    # --- GENERATE THE MESH AND VALUES FOR EACH KEY ---
    plot_structures = {}
    all_global_values = [] 

    for key in key_list:
        columns_data = sorted(species_data[key], key=lambda x: x['R'])
        if not columns_data:
            continue
            
        polygons = []
        values = []
        radii = [col['R'] for col in columns_data]
        
        # Calculate horizontal mesh grid edges
        if len(radii) > 1:
            r_midshifts = 0.5 * np.diff(radii)
            r_edges = [radii[0] - r_midshifts[0]] + [radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [radii[-1] + r_midshifts[-1]]
        else:
            r_edges = [radii[0] - 0.5, radii[0] + 0.5]
            
        for i, col in enumerate(columns_data):
            r_left, r_right = r_edges[i], r_edges[i+1]
            z_pts, v_pts = col['z'], col['v']
            
            # Calculate vertical mesh grid edges (assuming z decreases down to midplane)
            z_midshifts = 0.5 * np.diff(z_pts)
            z_edges = [z_pts[0] - z_midshifts[0]] + [z_pts[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_pts[-1] + z_midshifts[-1])]
            
            for j in range(len(v_pts)):
                poly = [(r_left, z_edges[j]), (r_right, z_edges[j]), (r_right, z_edges[j+1]), (r_left, z_edges[j+1])]
                polygons.append(poly)
                values.append(v_pts[j])
                
        vals_array = np.array(values)
        plot_structures[key] = {'polygons': polygons, 'values': vals_array, 'radii': radii, 'all_z': np.concatenate([c['z'] for c in columns_data])}
        
        if common_scale:
            all_global_values.extend(values)

    if not plot_structures:
        if verbose: print("No plottable data found.")
        return

    # --- MANAGEMENT OF SUBPLOT GRID GEOMETRY ---
    num_plots = len(plot_structures)
    
    # Configure canvas architecture dynamically depending on single vs multiple plot targets
    if num_plots == 1:
        fig, ax = plt.subplots(figsize=(10, 6))
        axes = [ax]
    else:
        cols = min(3, num_plots)  
        rows = (num_plots + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4.5 * rows), squeeze=False)
        axes = axes.flatten() 

    # --- CALCULATION OF THE COMMON SCALE (IF ACTIVE) ---
    global_vmin, global_vmax = None, None
    if common_scale and all_global_values and num_plots > 1:
        all_global_values = np.array(all_global_values)
        has_density = any("density" in k.lower() for k in key_list)
        if MODE == 'chemistry' or has_density:
            pos_vals = all_global_values[all_global_values > 0]
            global_vmin = vmin if vmin is not None else (max(1e-15, pos_vals.min()) if len(pos_vals) > 0 else 1e-15)
        else:
            global_vmin = vmin if vmin is not None else all_global_values.min()
        global_vmax = vmax if vmax is not None else all_global_values.max()

    # --- SUBPLOTS PLOTTING ---
    for idx, key in enumerate(key_list):
        if key not in plot_structures:
            continue
            
        ax = axes[idx]
        struct = plot_structures[key]
        vals = struct['values']
        
        # Calculate localized norm scaling limits if global synchronization is inactive
        if common_scale and num_plots > 1:
            actual_vmin, actual_vmax = global_vmin, global_vmax
        else:
            if MODE == 'chemistry' or "density" in key.lower():
                pos_vals = vals[vals > 0]
                actual_vmin = vmin if vmin is not None else (max(1e-15, pos_vals.min()) if len(pos_vals) > 0 else 1e-15)
            else:
                actual_vmin = vmin if vmin is not None else vals.min()
            actual_vmax = vmax if vmax is not None else vals.max()

        # Enforce LogNorm for chemical abundances or densities, otherwise standard Normalize
        if MODE == 'chemistry' or "density" in key.lower() or "extinction" in key.lower():
            color_norm = plt.cm.colors.LogNorm(vmin=actual_vmin, vmax=actual_vmax)
        else:
            color_norm = plt.cm.colors.Normalize(vmin=actual_vmin, vmax=actual_vmax)

        # Draw structured data via non-uniform discrete mesh panels
        coll = PolyCollection(struct['polygons'], array=vals, cmap=colormap, norm=color_norm, edgecolors='none')
        ax.add_collection(coll)
        
        # Format label annotations based on data mode
        if MODE == 'physical':
            lab = title_phys(key)
        elif MODE == 'chemistry':
            first_r = struct['radii'][0] if struct['radii'] else 5
            lab = title_mol(key, fracab, chempath / f"{int(first_r)}AU", verbose=verbose)
            
        sm = plt.cm.ScalarMappable(cmap=colormap, norm=color_norm)
        sm.set_array(vals)
        fig.colorbar(sm, ax=ax, label=lab)
        
        ax.set_xlabel('Radius R [AU]')
        ax.set_ylabel('Altitude z [AU]')
        ax.set_title(lab)
        
        ax.set_xlim(xlim if xlim is not None else (0, max(struct['radii']) * 1.02))
        ax.set_ylim(ylim if ylim is not None else (0, max(struct['all_z']) * 1.07))

    # Remove unneeded empty subplots inside non-perfect grids
    if num_plots > 1:
        for idx in range(num_plots, len(axes)):
            fig.delaxes(axes[idx])

    # Extract simulation timestamp metadata to configure localized headers or suptitle
    try:
        any_key = list(plot_structures.keys())[0]
        first_r = plot_structures[any_key]['radii'][0]
        time_seconds = main_output_dict[first_r]['abundances'].coords['time'].values[itime]
        
        if num_plots == 1:
            axes[0].set_title(f"{axes[0].get_title()} \n $t = {time_seconds/3.156e7:.2f}$ years")
        else:
            fig.suptitle(f'Simulation Output — $t = {time_seconds/3.156e7:.2f}$ years', fontsize=14, y=0.98)
    except:
        pass

    plt.tight_layout()
    plt.show()


def plot_midplane_nautilus_multi(main_output_dict,
                                 itime=-1,
                                 MODE='chemistry',
                                 KEY_NAMES=['CO'],
                                 VARIABLE_LABEL="Fractional Abundance [$n_{X}/n_{H}$]",
                                 fracab=True,
                                 verbose=True,
                                 xlim=None,
                                 ylim=None,
                                 cmap_name="turbo",
                                 vmin=None,
                                 vmax=None):
    """
    Plots 1D radial profiles of multiple variables/species strictly at the disk midplane (z = 0).

    Args:
        main_output_dict (dict): Master dictionary storing the simulation outputs,
            where keys are radii (int/float) and values are sub-dictionaries.
        itime (int, optional): Index of the simulation timestep to visualize. 
            Defaults to -1 (the final timestep).
        MODE (str, optional): Type of variable to plot. Options are 'chemistry' 
            or 'physical'. Defaults to 'chemistry'.
        KEY_NAMES (list of str, optional): List of dict keys for 'physical' variables, 
            or chemical species formula strings for 'chemistry' mode. Defaults to ['CO'].
        VARIABLE_LABEL (str, optional): Label displayed along the vertical axis. 
            Defaults to "Fractional Abundance [$n_{X}/n_{H}$]".
        fracab (bool, optional): If True, plots raw fractional abundances. If False, 
            multiplies abundances by the total hydrogen number density (nH). Defaults to True.
        verbose (bool, optional): If True, prints diagnostic processing errors. Defaults to True.
        xlim (tuple of float, optional): Custom (min, max) boundaries for the Radius axis.
        ylim (tuple of float, optional): Custom (min, max) boundaries for the vertical axis.
        cmap_name (str, optional): Matplotlib colormap name used to color the different lines.
            Defaults to "turbo".
        vmin (float, optional): Forced lower bound for the vertical scale.
        vmax (float, optional): Forced upper bound for the vertical scale.

    Returns:
        None: Displays a Matplotlib pyplot 1D line figure with multiple curves.
    """
    # Force KEY_NAMES to be a list if the user accidentally passes a single string
    if isinstance(KEY_NAMES, str):
        KEY_NAMES = [KEY_NAMES]

    # --- PLOT INITIALIZATION ---
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Generate colors from the specified colormap
    # If there is only one species, pick the center of the colormap; otherwise, distribute evenly
    if len(KEY_NAMES) == 1:
        colors = [plt.get_cmap(cmap_name)(0.5)]
    else:
        colors = plt.get_cmap(cmap_name)(np.linspace(0, 1, len(KEY_NAMES)))

    # Track valid min and max values across all species to handle log-scale limits correctly
    all_valid_mins = []
    all_valid_maxs = []

    # --- LOOP OVER EACH SPECIES / VARIABLE ---
    for idx, key in enumerate(KEY_NAMES):
        radii_list = []
        values_list = []

        for r_value in main_output_dict.keys():
            try:
                sub_dict = main_output_dict[r_value]
                
                # CRITICAL: In Nautilus 1D grid arrays, index -1 represents the midplane (z = 0)
                MIDPLANE_INDEX = -1 
                
                # Extract the single value at the midplane
                if MODE == 'physical':
                    full_array = sub_dict[key]
                    v_midplane = full_array[itime, MIDPLANE_INDEX]
                elif MODE == 'chemistry':
                    abundance_array = sub_dict['abundances']
                    v_midplane = float(abundance_array.isel(time=itime).sel(species=key).values[MIDPLANE_INDEX])
                    
                    if not fracab:
                        nH_midplane = sub_dict["H_number_density"][itime, MIDPLANE_INDEX]
                        v_midplane = v_midplane * nH_midplane
                
                radii_list.append(float(r_value))
                values_list.append(v_midplane)
                
            except Exception as e:
                if verbose:
                    print(f"Error processing midplane data for R={r_value}, Key={key}: {e}")

        if not radii_list:
            if verbose:
                print(f"No data collected for Key: {key}. Skipping.")
            continue

        # Sort arrays by radius to ensure the plotted line connects points sequentially
        sort_indices = np.argsort(radii_list)
        radii_arr = np.array(radii_list)[sort_indices]
        values_arr = np.array(values_list)[sort_indices]

        # Store min/max values for automatic Y-axis scaling
        pos_values = values_arr[values_arr > 0]
        if len(pos_values) > 0:
            all_valid_mins.append(pos_values.min())
        all_valid_maxs.append(values_arr.max())

        # Plot line and marker dots for the current species
        ax.plot(radii_arr, values_arr, color=colors[idx], linestyle='-', marker='o', 
                markersize=5, linewidth=1.5, label=key)

    # --- SCALE AND AXIS LIMIT MANAGEMENT ---
    # Determine whether to use a logarithmic vertical axis
    is_log = MODE == 'chemistry' or any("density" in k.lower() or "extinction" in k.lower() for k in KEY_NAMES)
    
    if is_log:
        ax.set_yscale('log')
        # Fallback value checking to prevent log(0) crash loops across all curves
        global_min = min(all_valid_mins) if all_valid_mins else 1e-15
        global_max = max(all_valid_maxs) if all_valid_maxs else 1.0
        
        actual_vmin = vmin if vmin is not None else max(1e-15, global_min)
        actual_vmax = vmax if vmax is not None else global_max
        ax.set_ylim(actual_vmin, actual_vmax)
    else:
        if vmin is not None or vmax is not None:
            ax.set_ylim(vmin, vmax)

    # --- LABELS AND TITLE CONFIGURATION ---
    ax.set_xlabel('Radius R [AU]')
    ax.set_ylabel(VARIABLE_LABEL)
    
    # Attempt to read and parse the time coordinate in years
    try:
        first_r = list(main_output_dict.keys())[0]
        time_seconds = main_output_dict[first_r]['abundances'].coords['time'].values[itime]
        ax.set_title(f'Midplane ($z = 0$) Profile at $t = {time_seconds/3.156e7:.2e}$ years')
    except:
        ax.set_title('Midplane ($z = 0$) Profile')

    # Apply manual axis overrides if supplied
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)

    # Add legend to identify the plotted species
    ax.legend(loc='best', frameon=True, shadow=False)
    
    ax.grid(True, linestyle=':', alpha=0.5)
    plt.tight_layout()
    plt.show()


def plot_grain_surface_midplane(chempath,
                                main_output_dict,
                                itime=-1,
                                verbose=True,
                                xlim=None,
                                ylim=None,
                                color="darkgreen"):
    """
    Plots the total available grain surface area (sum over all bins of: 4*pi*a^2 * n_grain)
    strictly at the disk midplane (z = 0) as a function of Radius.
    """
    radii_list = []
    surface_list = []

    # Loop over the radius keys (5, 10, 15, etc.)
    for r_value in main_output_dict.keys():
        folder_name = f"{r_value}AU"
        file_path = os.path.join(chempath, folder_name, "1D_grain_sizes.in")
        
        if os.path.exists(file_path):
            try:
                # 1. Load the data from 1D_grain_sizes.in (ignoring comments starting with '!')
                # Every line is a spatial point from atmosphere down to midplane
                grain_data = np.loadtxt(file_path, comments='!')
                
                # Extract only the midplane layer (the last row, index -1)
                midplane_row = grain_data[-1, :]
                
                # 2. Determine the number of grain bins (N)
                total_columns = len(midplane_row)
                N = int(total_columns / 4)
                
                # 3. Extract group 1 (grain radii 'a') and group 2 ('GTODN' values)
                a_array = midplane_row[0:N]          # First N columns [cm]
                gtodn_array = midplane_row[N:2*N]    # Next N columns [dimensionless]
                
                # 4. Get the Hydrogen number density (nH) at the midplane from main_output_dict
                sub_dict = main_output_dict[r_value]
                nH_midplane = sub_dict["H_number_density"][itime, -1]
                
                # 5. Compute the physical formula: Sum over all bins of (4 * pi * a^2 * nH / GTODN)
                # We use numpy vectorized operations for the sum
                total_surface = 4 * np.pi * nH_midplane * np.sum((a_array**2) / gtodn_array)
                
                radii_list.append(float(r_value))
                surface_list.append(total_surface)
                
            except Exception as e:
                if verbose:
                    print(f"Error processing grain data for R={r_value}: {e}")
        else:
            if verbose:
                print(f"File not found: {file_path}")

    if not radii_list:
        if verbose:
            print("No data collected. Check paths or files.")
        return

    # Sort arrays by radius for a clean line plot
    sort_indices = np.argsort(radii_list)
    radii_arr = np.array(radii_list)[sort_indices]
    surface_arr = np.array(surface_list)[sort_indices]

    # --- PLOTTING ---
    fig, ax = plt.subplots(figsize=(9, 5))
    
    # Grain surfaces usually vary over orders of magnitude, a log scale is ideal
    ax.set_yscale('log')
    
    # Plot line + points
    ax.plot(radii_arr, surface_arr, color=color, linestyle='-', marker='s', markersize=5, linewidth=1.5)

    # Labels and formatting
    ax.set_xlabel('Radius R [AU]')
    ax.set_ylabel(r'Total Grain Surface Area $\sim\text{cm}^{2}/\text{cm}^{3}$')
    
    try:
        time_seconds = main_output_dict[radii_list[0]]['abundances'].coords['time'].values[itime]
        ax.set_title(f'Total Grain Surface Area at Midplane ($z=0$) - $t = {time_seconds/3.156e7:.2e}$ yr')
    except:
        ax.set_title('Total Grain Surface Area at Midplane ($z=0$)')

    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)
    
    if ylim is None:
        # Avoid log(0) display errors
        actual_vmin = max(1e-25, surface_arr[surface_arr > 0].min() if any(surface_arr > 0) else 1e-25)
        ax.set_ylim(actual_vmin, surface_arr.max() * 2)

    ax.grid(True, linestyle=':', alpha=0.5)
    plt.tight_layout()
    plt.show()

def plot_vertical_cut_nautilus(chempath,
                              main_output_dict,
                              R,
                              species='CO',
                              itime=-1,
                              fracab=True,
                              col='royalblue',
                              xlim=None,
                              ylim=None,
                              xscale="linear",
                              yscale="linear"):
    """
    Plots the vertical profile (abundance vs. height z) for a given species 
    at a specific disk radius (R) using NAUTILUS simulation outputs.

    Parameters:
    -----------
    chempath : str
        Path to NMGC chemistry outputs
    main_output_dict : dict
        Dictionary containing the NAUTILUS simulation outputs (e.g. pipe.chemistry)
    R : int or float
        The specific radius (in AU) to extract data for.
    species : str, optional
        The chemical species to plot (default is 'CO').
    itime : int, optional
        The time index to plot (default is -1, which corresponds to the last timestep).
    fracab : bool, optional
        If True, plots fractional abundance relative to H (n_sp/nH). 
        If False, plots absolute volume density (cm^-3). Default is True.
    col : str, optional
        Color of the plot line and markers (default is 'royalblue').
    xlim, ylim : tuple, optional
        Limits for the x and y axes (e.g., (min, max)).
    xscale, yscale : str, optional
        Scale of the axes (default is "linear").
    """

    # Check if the requested radius R exists in the provided dictionary keys
    if R not in list(main_output_dict.keys()):
        raise ValueError("R does not exist. Please use an existing R that you can find in list(main_output_dict.keys())")
        
    # Extract the abundance DataArray for the specified radius R
    ab = main_output_dict[R]['abundances']    # DataArray shape: (nb_timesteps, nb_species, nz)
    
    # Select the data for the chosen timestep (itime) and chemical species
    sp_last = ab.isel(time=itime).sel(species=species)   
    sp_arr = sp_last.values                 # Convert to numpy array, shape: (nz,)
    
    # Determine whether to plot absolute density or fractional abundance
    if fracab == False:
        # Get the hydrogen number density for the last physical condition to convert abundance to density
        nH     = pipe.chemistry[30]['H_number_density'][-1]   # shape: (nz,)
        n_sp   = nH * sp_arr
        xlabel = f'n({species}) [cm$^{{-3}}$]'
    else :
        # Use fractional abundance directly
        n_sp = sp_arr
        xlabel = f'n({species})/n$_H$'

    # Load the vertical grid (z) from the static 1D structure file for this radius
    static = pd.read_table(
        f'{chempath}/{R}AU/1D_static.dat',
        sep=r'\s+', comment='!', header=None, engine='python'
    )
    z = static[0].values   # z coordinate in [AU], ranging from surface to midplane
    
    # Retrieve the physical time in seconds for the selected timestep
    time_seconds = main_output_dict[R]['abundances'].coords['time'].values[itime]
    
    # Initialize the matplotlib figure
    fig = plt.figure(figsize=(7, 5))
    
    # Plot the data using both scatter markers and a continuous line
    plt.scatter(n_sp, z, color=col)
    plt.plot(n_sp, z, color=col)
    
    # Apply labels, grid, title, limits, and scales
    plt.xlabel(xlabel)
    plt.ylabel("z [AU]")
    plt.grid(True, linestyle=':', alpha=0.5)
    # Convert time from seconds to years in the title display
    plt.title(f'$R = {R} $ AU - $t = {time_seconds/3.156e7:.2e}$ yr')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xscale(xscale)
    plt.yscale(yscale)
    
    # Adjust layout and display the plot
    plt.tight_layout()
    plt.show()


def plot_atom_ratio_nautilus(chempath,
                             main_output_dict,
                             s1="C",
                             s2="O",
                             itime=-1,
                             verbose=True,
                             xlim=None,
                             ylim=None,
                             colormap="gnuplot",
                             vmin=None,
                             vmax=None):
    """Plots a 2D vertical cross-section of the elemental abundance ratio of two atoms.

    This function builds a non-uniform structured grid (poloidal cut) from discrete 
    radial column blocks derived from Nautilus physical simulation outputs. It 
    dynamically parses molecular formulas for all gas-phase species found at each 
    radius, aggregates the total elemental counts of `s1` and `s2` per grid cell 
    using an optimized matrix multiplication cache, and renders the calculated 
    spatial distribution ratio (s1/s2) via a Matplotlib Polygon Collection.

    Args:
        chempath (str): Path to the parent directory containing the radial simulation 
            folders (e.g., "5AU", "10AU").
        main_output_dict (dict): Master simulation dictionary where keys are radii 
            (int/float) and values are sub-dictionaries containing spatial grid and 
            chemical data properties (e.g., pipe.chemistry).
        s1 (str, optional): Chemical symbol of the numerator element. Defaults to "C".
        s2 (str, optional): Chemical symbol of the denominator element. Defaults to "O".
        itime (int, optional): Index of the simulation timestep to visualize. 
            Defaults to -1 (the final simulation timestep).
        verbose (bool, optional): If True, prints diagnostic mismatched configurations, 
            missing file warnings, or alerts regarding grid cells containing zero `s2` 
            atoms to the console. Defaults to True.
        xlim (tuple of float, optional): Custom (min, max) boundaries for the horizontal 
            Radius axis. Defaults to None (automatically bound to grid geometry).
        ylim (tuple of float, optional): Custom (min, max) boundaries for the vertical 
            Altitude axis. Defaults to None (automatically bound to grid geometry).
        colormap (str, optional): Matplotlib colormap string used to style the discrete 
            mesh and colorbar scale. Defaults to "gnuplot".
        vmin (float, optional): Forced lower bound for the colorbar scale. 
            Defaults to None (computed automatically from the data minimum).
        vmax (float, optional): Forced upper bound for the colorbar scale. 
            Defaults to None (computed automatically from the data maximum).

    Raises:
        ValueError: If either `s1` or `s2` is not present in the allowed elemental 
            network base array ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 
            'Mg', 'Cl', 'P', 'F'].

    Returns:
        None: Displays a Matplotlib pyplot figure window.
    """

    
    elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    if s1 not in elements or s2 not in elements:
        raise ValueError("Please check your atoms, one of them is not existing in the model")
        
    variable_label = f"Atomic Ratio [{s1}/{s2}]"
    
    # --- INTERNAL HELPER ---
    def count_species_elements(species_name, element1, element2):
        """Return a dictionary with the counts of element1 and element2 in the given chemical species."""
        if species_name == 'e-': return {element1: 0, element2: 0}
        formula = species_name.replace('c-', '').replace('l-', '')
        if formula.endswith('+') or formula.endswith('-'):
            formula = formula[:-1]
        pattern = re.compile(r'([A-Z][a-z]?)(-?\d*)')
        composition = {}
        for atom, n in pattern.findall(formula):
            count = int(n) if n else 1
            composition[atom] = composition.get(atom, 0) + count
        return {
            element1: composition.get(element1, 0),
            element2: composition.get(element2, 0)}

    def keep_gas_species_only(species):
        motif = re.compile(r'^[JK]\d{2}|^GRAIN')
        return [e for e in species if not motif.match(e)]

    # --- GLOBAL CACHE INITIALIZATION ---
    atom_cache = {}

    # --- EXTRACT DATA BY COLUMN (RADIUS) ---
    columns_data = []

    for r_value in main_output_dict.keys():
        folder_name = f"{r_value}AU"
        file_path = os.path.join(chempath, folder_name, "1D_static.dat")
        
        if os.path.exists(file_path):
            try:
                z_points = np.loadtxt(file_path, comments='!', usecols=0)
                sub_dict = main_output_dict[r_value]
                abundance_array = sub_dict['abundances']
                
                local_species_list = keep_gas_species_only(list(abundance_array.coords['species'].values))
                
                s1_coeffs = []
                s2_coeffs = []
                for species in local_species_list:
                    if species not in atom_cache:
                        counts = count_species_elements(species, s1, s2)
                        atom_cache[species] = (counts[s1], counts[s2])
                    
                    c1, c2 = atom_cache[species]
                    s1_coeffs.append(c1)
                    s2_coeffs.append(c2)
                
                s1_coeffs = np.array(s1_coeffs)[:, np.newaxis]
                s2_coeffs = np.array(s2_coeffs)[:, np.newaxis]
                
                sliced_abundances = abundance_array.isel(time=itime).sel(species=local_species_list).values
                
                total_s1 = np.sum(sliced_abundances * s1_coeffs, axis=0)
                total_s2 = np.sum(sliced_abundances * s2_coeffs, axis=0)

                if verbose:
                    zero_indices = np.where(total_s2 == 0)[0]
                    if len(zero_indices) > 0:
                        # Retrieve the corresponding altitudes to make the message useful
                        altitudes_zero = z_points[zero_indices]
                        print(f"[R={r_value} AU] No {s2} atoms detected in {len(zero_indices)} vertical grid cell(s).")
                        print(f"   -> Affected altitudes [AU]: {altitudes_zero}")
                
                with np.errstate(divide='ignore', invalid='ignore'):
                    v_points = np.where(total_s2 > 0, total_s1 / total_s2, 0.0)
                
                if len(z_points) == len(v_points):
                    columns_data.append({
                        'R': float(r_value),
                        'z': np.array(z_points),
                        'v': np.array(v_points)
                    })
                elif verbose: 
                    print(f"Size mismatch for R={r_value}: file has {len(z_points)} points, ratio has {len(v_points)}.")
                        
            except Exception as e:
                if verbose: 
                    print(f"Error processing data for R={r_value}: {e}")
        elif verbose: 
            print(f"File not found: {file_path}")
    
    columns_data = sorted(columns_data, key=lambda x: x['R'])
    
    # --- GENERATE THE POLYGON MESH ---
    polygons = []
    values = []
    
    radii = [col['R'] for col in columns_data]
    r_edges = []
    if len(radii) > 1:
        r_midshifts = 0.5 * np.diff(radii)
        r_edges.append(radii[0] - r_midshifts[0])
        for i in range(len(r_midshifts)):
            r_edges.append(radii[i] + r_midshifts[i])
        r_edges.append(radii[-1] + r_midshifts[-1])
    elif len(radii) == 1:
        r_edges = [radii[0] - 0.5, radii[0] + 0.5]
    
    for i, col in enumerate(columns_data):
        r_left = r_edges[i]
        r_right = r_edges[i+1]
        
        z_pts = col['z']
        v_pts = col['v']
        
        z_edges = []
        z_midshifts = 0.5 * np.diff(z_pts)
        z_edges.append(z_pts[0] - z_midshifts[0]) 
        for j in range(len(z_midshifts)):
            z_edges.append(z_pts[j] + z_midshifts[j])
        z_edges.append(max(0.0, z_pts[-1] + z_midshifts[-1])) 
    
        for j in range(len(v_pts)):
            z_top = z_edges[j]
            z_bottom = z_edges[j+1]
            
            poly = [
                (r_left, z_top),
                (r_right, z_top),
                (r_right, z_bottom),
                (r_left, z_bottom)
            ]
            polygons.append(poly)
            values.append(v_pts[j])
    
    # --- PLOTTING ---
    if not polygons:
        if verbose: 
            print("No polygons generated. Check your configuration parameters.")
    else:
        fig, ax = plt.subplots(figsize=(10, 6))
        values = np.array(values)
        
        # 1. Determine data bounds (filtering out <= 0 values for logarithmic checks)
        positive_values = values[values > 0]
        
        actual_vmin = vmin if vmin is not None else values.min()
        actual_vmax = vmax if vmax is not None else values.max()

        # 2. Dynamic normalization selection
        # Check if positive data exists and if the span exceeds one order of magnitude (factor of 10)
        if len(positive_values) > 0 and (actual_vmax / positive_values.min()) > 10.0:
            # Safety check: LogNorm vmin must be strictly greater than 0
            log_vmin = vmin if vmin is not None else positive_values.min()
            color_norm = plt.cm.colors.LogNorm(vmin=log_vmin, vmax=actual_vmax)
            if verbose:
                print(f"[Plot Scale] Dynamic switch to LOGARITHMIC scale (span > 1 order of magnitude).")
        else:
            color_norm = plt.cm.colors.Normalize(vmin=actual_vmin, vmax=actual_vmax)
            if verbose:
                print(f"[Plot Scale] Dynamic switch to LINEAR scale (span <= 1 order of magnitude).")
    
        # 3. Apply the polygon mesh grid and colorbar
        coll = PolyCollection(polygons, array=values, cmap=colormap, norm=color_norm, edgecolors='none')
        ax.add_collection(coll)
    
        sm = plt.cm.ScalarMappable(cmap=colormap, norm=color_norm)
        sm.set_array(values)
        fig.colorbar(sm, ax=ax, label=variable_label)
    
        ax.set_xlabel('Radius R [AU]')
        ax.set_ylabel('Altitude z [AU]')
        
        try:
            first_r = radii[0]
            time_seconds = main_output_dict[first_r]['abundances'].coords['time'].values[itime]
            ax.set_title(f'$t = {time_seconds/3.156e7:.2e}$ years (Ratio {s1}/{s2})')
        except:
            ax.set_title(f'Atomic Ratio {s1}/{s2}')
    
        all_z = np.concatenate([col['z'] for col in columns_data])
        ax.set_xlim(xlim if xlim is not None else (0, max(radii) * 1.02))
        ax.set_ylim(ylim if ylim is not None else (0, max(all_z) * 1.07))
        
        plt.show()


def plot_top_contributing_species(chempath,
                                  main_output_dict,
                                  target_atom="C",
                                  itime=-1,
                                  verbose=True,
                                  spnumber=5,
                                  color="darkred",
                                  phase="gas",
                                  grain_bin=None):
    """
    Computes and plots the top chemical species contributing to the global, 
    volume-integrated total budget of a target chemical element within a protoplanetary disk.

    Parameters:
    -----------
    chempath : str
        Path to the directory containing spatial grid subfolders (e.g., '10AU/', '100AU/').
    main_output_dict : dict
        Nested dictionary mapping radial keys to data sub-structures containing 'abundances' 
        (xarray.DataArray) and 'H_number_density' spatial profiles.
    target_atom : str, default "C"
        The chemical element symbol whose structural reservoir distribution is evaluated.
    itime : int, default -1
        Time index to slice from the multi-epoch data arrays.
    verbose : bool, default True
        If True, outputs warning and missing file notifications to the console.
    spnumber : int, default 5
        Number of top contributing chemical species to display on the bar chart.
    color : str, default "darkred"
        Visual fill color of the generated bar elements.
    phase : str, default "gas"
        The chemical phase environment to isolate ("gas", "surface", "mantle", "grain", or "all").
    grain_bin : int or str, optional
        Specific grain size category to filter when analyzing surface or mantle ice matrices.
    """
    
    allowed_elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    if target_atom not in allowed_elements:
        raise ValueError(f"Target atom '{target_atom}' is not recognized in the chemical network.")
    
    # Updated phase validation list to hold 5 phases
    valid_phases = ["gas", "surface", "mantle", "grain", "all"]
    if phase not in valid_phases:
        raise ValueError(f"Phase '{phase}' unrecognized. Choose among {valid_phases}")
    if grain_bin is not None and phase in ["gas","all"] : raise ValueError("grain_bin and gas phase can not be defined simultaneously")

    # --- INTERNAL HELPERS ---
    def parse_species(species_name):
        if species_name == 'e-':
            return "gas", None, "e-"
            
        grain_match = re.match(r'^([JK])(\d+)(.+)', species_name)
        
        if grain_match:
            p_code, g_bin, raw_formula = grain_match.groups()
            sp_phase = "surface" if p_code == 'J' else "mantle"
            sp_bin = g_bin
        else:
            sp_phase = "gas"
            sp_bin = None
            raw_formula = species_name
        
        # Remove structural isomer descriptors while keeping trailing charge states intact
        clean_formula = raw_formula.replace('c-', '').replace('l-', '')
        return sp_phase, sp_bin, clean_formula

    def count_target_atom(clean_formula, target):
        if clean_formula == 'e-': return 0
        
        calc_formula = clean_formula
        if calc_formula.endswith('+') or calc_formula.endswith('-'):
            calc_formula = calc_formula[:-1]
            
        # Strip out any remaining internal dashes to avoid messing up structural regex lookups
        calc_formula = calc_formula.replace('-', '')
            
        pattern = re.compile(r'([A-Z][a-z]?)(\d*)') 
        composition = {}
        for atom, n in pattern.findall(calc_formula):
            # Fallback mechanism to safely capture malformed/unconventional formatting rules
            try:
                count = int(n) if n else 1
            except ValueError:
                count = 1
            composition[atom] = composition.get(atom, 0) + count
        return composition.get(target, 0)

    def filter_and_cache_species(species_list):
        valid_list = []
        coeffs = []
        
        for sp in species_list:
            if "GRAIN" in sp:
                continue
                
            sp_phase, sp_bin, clean_formula = parse_species(sp)
            
            # --- PHASE FILTERING LOGIC ---
            # "all" passes gas, surface, and mantle environments smoothly
            if phase == "all":
                pass
            # "grain" isolates ice layers (both surface and mantle components)
            elif phase == "grain":
                if sp_phase not in ["surface", "mantle"]:
                    continue
            # Exact key matching for "gas", "surface", or "mantle"
            elif sp_phase != phase:
                continue
                
            if grain_bin is not None and sp_phase in ["surface", "mantle"]:
                if sp_bin != str(grain_bin):
                    continue
            
            coef = count_target_atom(clean_formula, target_atom)
            if coef > 0:
                valid_list.append(sp)
                coeffs.append(coef)
                
        return valid_list, np.array(coeffs)

    global_species_contributions = {}
    AU_to_cm = 1.496e13  

    # --- 1. RECONSTRUCT RADIAL EDGES FOR VOLUME ---
    # Maintains the structural link between the extracted integer and the original dictionary keys
    radii_map = {}
    for original_key in main_output_dict.keys():
        digits = re.findall(r'\d+', str(original_key))
        if digits:
            try:
                rad_int = int(digits[0])
                radii_map[rad_int] = original_key
            except ValueError:
                continue
                
    radii = sorted(list(radii_map.keys()))
    
    if len(radii) == 0:
        if verbose: print("Main output dictionary is empty or contains no valid numerical keys.")
        return
        
    if len(radii) > 1:
        r_midshifts = 0.5 * np.diff(radii)
        r_edges = [radii[0] - r_midshifts[0]] + [radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [radii[-1] + r_midshifts[-1]]
    else:
        r_edges = [radii[0] * 0.9, radii[0] * 1.1]

    # --- 2. DATA ACCUMULATION WITH PHYSICAL VOLUMES ---
    for i, r_value in enumerate(radii):
        folder_name = f"{r_value}AU"
        file_path = os.path.join(chempath, folder_name, "1D_static.dat")
        
        if os.path.exists(file_path):
            try:
                z_points = np.loadtxt(file_path, comments='!', usecols=0)
                
                # Retrieve the exact original key mapping (e.g. string vs int representations)
                orig_key = radii_map[r_value]
                sub_dict = main_output_dict[orig_key]
                
                abundance_array = sub_dict['abundances']
                nH_profile = sub_dict["H_number_density"][itime,:]
                
                if len(z_points) > 1:
                    z_midshifts = 0.5 * np.diff(z_points)
                    z_edges = [z_points[0] - z_midshifts[0]] + [z_points[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_points[-1] + z_midshifts[-1])]
                    dz = np.abs(np.diff(z_edges))
                else:
                    dz = np.array([z_points[0] if z_points[0] > 0 else 1.0])
                
                r_left = r_edges[i] * AU_to_cm
                r_right = r_edges[i+1] * AU_to_cm
                dR = r_right - r_left
                R_center = float(r_value) * AU_to_cm
                
                # Evaluate cylindrical domain shell spatial volume elements: 2 * pi * R * dR * dZ
                cell_volumes = 2 * np.pi * R_center * dR * (dz * AU_to_cm) 
                
                raw_species = list(abundance_array.coords['species'].values)
                local_species_list, target_coeffs = filter_and_cache_species(raw_species)
                
                if not local_species_list:
                    continue 
                
                y_abundances = abundance_array.isel(time=itime).sel(species=local_species_list).values
                
                nH_2d = nH_profile[np.newaxis, :]     
                volumes_2d = cell_volumes[np.newaxis, :] 
                
                # Scale relative molecular fractional abundance values into raw physical particle pools
                physical_atoms = y_abundances * nH_2d * volumes_2d
                total_column_atoms = np.sum(physical_atoms, axis=1)
                total_column_atoms = np.sum(physical_atoms, axis=1)
                contributions_per_species = total_column_atoms * target_coeffs
                
                for idx, species in enumerate(local_species_list):
                    if contributions_per_species[idx] > 0:
                        _, _, clean_formula = parse_species(species)
                        global_species_contributions[clean_formula] = global_species_contributions.get(clean_formula, 0.0) + contributions_per_species[idx]
                        
            except Exception as e:
                if verbose: print(f"Error processing data for R={r_value}: {e}")
        elif verbose: 
            print(f"File not found: {file_path}")

    # --- 3. STATISTICAL PROCESSING & PLOTTING ---
    total_atoms_sum = sum(global_species_contributions.values())
    if total_atoms_sum == 0:
        bin_str = f" (Bin: {grain_bin})" if grain_bin else ""
        print(f"No {target_atom} atoms detected within the chosen criteria [Phase: {phase}{bin_str}].")
        return

    species_percentages = {sp: (val / total_atoms_sum) * 100 for sp, val in global_species_contributions.items()}
    sorted_species = sorted(species_percentages.items(), key=lambda x: x[1], reverse=True)
    
    top = sorted_species[:spnumber]
    others_sum = 100.0 - sum(val for sp, val in top)

    labels = [item[0] for item in top]
    percentages = [item[1] for item in top]

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(labels, percentages, color=color, edgecolor='grey', alpha=0.85)
    
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.4f}%',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 4),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Handle conditional subtitle parsing across the 5 dynamic phases
    is_ice_phase = phase in ["surface", "mantle", "grain"]
    bin_title = f" (Bin {grain_bin})" if (grain_bin and is_ice_phase) else (" (All Grains)" if is_ice_phase else "")
    phase_title = f"Phase: {phase.upper()}{bin_title}"

    ax.set_ylabel(f'Global {phase_title} Budget Contribution (%)', fontsize=11)
    ax.set_xlabel(f'Chemical Species (Top {spnumber})', fontsize=12)
    
    try:
        first_r = radii[0]
        time_seconds = main_output_dict[radii_map[first_r]]['abundances'].coords['time'].values[itime]
        ax.set_title(f'Top {spnumber} Reservoirs for Element [{target_atom}] — {phase_title}\n$t = {time_seconds/3.156e7:.2e}$ years (Volume Integrated)', fontsize=13, pad=15)
    except:
        ax.set_title(f'Top {spnumber} Reservoirs for Element [{target_atom}] — {phase_title}', fontsize=13, pad=15)
        
    ax.set_ylim(0, max(percentages) * 1.15 if percentages else 100)
    ax.grid(axis='y', linestyle='--', alpha=0.4)

    if others_sum > 0.0001:
        plt.figtext(0.15, -0.01, f"* Omitted minor species account for {others_sum:.4f}% of the filtered physical {target_atom} budget.", 
                    fontsize=10, style='italic')

    plt.tight_layout()
    plt.show()

def plot_top_species_per_radius(chempath,
                                main_output_dict,
                                target_atom="C",
                                itime=-1,
                                verbose=True,
                                spnumber=5,
                                phase="gas",
                                grain_bin=None,
                                cmap_name="tab10",
                                rmin=None,
                                rmax=None):
    """
    Computes and plots the top N contributing chemical species for a target element
    as a function of disk radius within an optional [rmin, rmax] radial range. 
    Displays a horizontal bar chart grouped by radius on the Y-axis, with local budget 
    percentages on the X-axis. Each unique chemical species is mapped to a distinct 
    color from the chosen colormap.

    Parameters:
    -----------
    chempath : str
        Path to the directory containing spatial grid subfolders (e.g., '10AU/', '100AU/').
    main_output_dict : dict
        Nested dictionary mapping radial keys to data sub-structures containing 'abundances' 
        and 'H_number_density' spatial profiles.
    target_atom : str, default "C"
        The chemical element symbol whose structural reservoir distribution is evaluated.
    itime : int, default -1
        Time index to slice from the multi-epoch data arrays.
    verbose : bool, default True
        If True, outputs warning and missing file notifications to the console.
    spnumber : int, default 5
        Number of top contributing chemical species to display for each radius.
    phase : str, default "gas"
        The chemical phase environment to isolate ("gas", "surface", "mantle", "grain", or "all").
    grain_bin : int or str, optional
        Specific grain size category to filter when analyzing surface or mantle ice matrices.
    cmap_name : str, default "tab10"
        Name of the Matplotlib colormap used to dynamically assign unique colors to species.
    rmin : float or int, optional
        Minimum inclusive radius boundary (in AU) to process. If None, defaults to the disk innermost grid point.
    rmax : float or int, optional
        Maximum inclusive radius boundary (in AU) to process. If None, defaults to the disk outermost grid point.
    """
    
    # Allowed elements check
    allowed_elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    if target_atom not in allowed_elements:
        raise ValueError(f"Target atom '{target_atom}' is not recognized in the chemical network.")
    
    # Phase validation
    valid_phases = ["gas", "surface", "mantle", "grain", "all"]
    if phase not in valid_phases:
        raise ValueError(f"Phase '{phase}' unrecognized. Choose among {valid_phases}")
    if grain_bin is not None and phase in ["gas","all"] : raise ValueError("grain_bin and gas phase can not be defined simultaneously")

    # --- INTERNAL HELPERS ---
    def parse_species(species_name):
        if species_name == 'e-':
            return "gas", None, "e-"
            
        grain_match = re.match(r'^([JK])(\d+)(.+)', species_name)
        
        if grain_match:
            p_code, g_bin, raw_formula = grain_match.groups()
            sp_phase = "surface" if p_code == 'J' else "mantle"
            sp_bin = g_bin
        else:
            sp_phase = "gas"
            sp_bin = None
            raw_formula = species_name
        
        clean_formula = raw_formula.replace('c-', '').replace('l-', '')
        return sp_phase, sp_bin, clean_formula

    def count_target_atom(clean_formula, target):
        if clean_formula == 'e-': return 0
        
        calc_formula = clean_formula
        if calc_formula.endswith('+') or calc_formula.endswith('-'):
            calc_formula = calc_formula[:-1]
            
        calc_formula = calc_formula.replace('-', '')
            
        pattern = re.compile(r'([A-Z][a-z]?)(\d*)') 
        composition = {}
        for atom, n in pattern.findall(calc_formula):
            try:
                count = int(n) if n else 1
            except ValueError:
                count = 1
            composition[atom] = composition.get(atom, 0) + count
        return composition.get(target, 0)

    def filter_and_cache_species(species_list):
        valid_list = []
        coeffs = []
        for sp in species_list:
            if "GRAIN" in sp:
                continue
            
            sp_phase, sp_bin, clean_formula = parse_species(sp)
            
            # Phase filtering
            if phase == "all":
                pass
            elif phase == "grain":
                if sp_phase not in ["surface", "mantle"]:
                    continue
            elif sp_phase != phase:
                continue
                
            # Size bin filtering
            if grain_bin is not None and sp_phase in ["surface", "mantle"]:
                if sp_bin != str(grain_bin):
                    continue
                    
            coef = count_target_atom(clean_formula, target_atom)
            if coef > 0:
                valid_list.append(sp)
                coeffs.append(coef)
        return valid_list, np.array(coeffs)

    AU_to_cm = 1.496e13  

    # --- 1. RECONSTRUCT RADIAL EDGES FOR VOLUME ---
    radii_map = {}
    for original_key in main_output_dict.keys():
        digits = re.findall(r'\d+', str(original_key))
        if digits:
            try:
                rad_int = int(digits[0])
                radii_map[rad_int] = original_key
            except ValueError:
                continue
                
    extracted_radii = sorted(list(radii_map.keys()))
    
    # --- RADIAL LIMITS FILTERING ---
    radii = []
    for r in extracted_radii:
        if rmin is not None and r < rmin:
            continue
        if rmax is not None and r > rmax:
            continue
        radii.append(r)
    
    if len(radii) == 0:
        if verbose: print(f"No grid radii found within the specified boundaries [rmin={rmin}, rmax={rmax}].")
        return
        
    if len(radii) > 1:
        r_midshifts = 0.5 * np.diff(radii)
        r_edges = [radii[0] - r_midshifts[0]] + [radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [radii[-1] + r_midshifts[-1]]
    else:
        r_edges = [radii[0] * 0.9, radii[0] * 1.1]

    # Structure to hold plot data per radius
    radial_plot_data = {}
    all_encountered_species = set()

    # --- 2. DATA ACCUMULATION PER CELL SHELL ---
    for i, r_value in enumerate(radii):
        folder_name = f"{r_value}AU"
        file_path = os.path.join(chempath, folder_name, "1D_static.dat")
        
        if os.path.exists(file_path):
            try:
                z_points = np.loadtxt(file_path, comments='!', usecols=0)
                orig_key = radii_map[r_value]
                sub_dict = main_output_dict[orig_key]
                
                abundance_array = sub_dict['abundances']
                nH_profile = sub_dict["H_number_density"][itime,:]
                
                if len(z_points) > 1:
                    z_midshifts = 0.5 * np.diff(z_points)
                    z_edges = [z_points[0] - z_midshifts[0]] + [z_points[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_points[-1] + z_midshifts[-1])]
                    dz = np.abs(np.diff(z_edges))
                else:
                    dz = np.array([z_points[0] if z_points[0] > 0 else 1.0])
                
                r_left = r_edges[i] * AU_to_cm
                r_right = r_edges[i+1] * AU_to_cm
                dR = r_right - r_left
                R_center = float(r_value) * AU_to_cm
                
                cell_volumes = 2 * np.pi * R_center * dR * (dz * AU_to_cm) 
                
                raw_species = list(abundance_array.coords['species'].values)
                local_species_list, target_coeffs = filter_and_cache_species(raw_species)
                
                if not local_species_list:
                    continue 
                
                y_abundances = abundance_array.isel(time=itime).sel(species=local_species_list).values
                
                nH_2d = nH_profile[np.newaxis, :]     
                volumes_2d = cell_volumes[np.newaxis, :] 
                
                physical_atoms = y_abundances * nH_2d * volumes_2d
                total_column_atoms = np.sum(physical_atoms, axis=1)
                contributions_per_species = total_column_atoms * target_coeffs
                
                local_radius_contributions = {}
                for idx, species in enumerate(local_species_list):
                    if contributions_per_species[idx] > 0:
                        _, _, clean_formula = parse_species(species)
                        local_radius_contributions[clean_formula] = local_radius_contributions.get(clean_formula, 0.0) + contributions_per_species[idx]
                
                radius_total_sum = sum(local_radius_contributions.values())
                if radius_total_sum > 0:
                    species_percentages = {sp: (val / radius_total_sum) * 100 for sp, val in local_radius_contributions.items()}
                    sorted_species = sorted(species_percentages.items(), key=lambda x: x[1], reverse=True)
                    top_entries = sorted_species[:spnumber]
                    radial_plot_data[r_value] = top_entries
                    
                    for species_name, _ in top_entries:
                        all_encountered_species.add(species_name)
                        
            except Exception as e:
                if verbose: print(f"Error processing data for R={r_value}: {e}")
        elif verbose: 
            print(f"File not found: {file_path}")

    if not radial_plot_data:
        print(f"No data available to plot for target atom {target_atom} within the specified radius range.")
        return

    # --- 3. DYNAMIC COLOR ASSIGNMENT VIA COLORMAP ---
    unique_species_list = sorted(list(all_encountered_species))
    num_unique_species = len(unique_species_list)
    
    try:
        cmap = plt.get_cmap(cmap_name)
    except ValueError:
        if verbose: print(f"Colormap '{cmap_name}' not found. Falling back to default 'tab10'.")
        cmap = plt.get_cmap("tab10")
        
    if hasattr(cmap, 'colors') and len(cmap.colors) >= num_unique_species:
        species_color_mapping = {sp: cmap(idx) for idx, sp in enumerate(unique_species_list)}
    else:
        color_indices = np.linspace(0, 1, num_unique_species) if num_unique_species > 1 else [0.0]
        species_color_mapping = {sp: cmap(color_indices[idx]) for idx, sp in enumerate(unique_species_list)}

    # --- 4. DYNAMIC GRAPHICAL CONSTRUCTION ---
    y_labels = []
    x_percentages = []
    bar_colors = []
    
    reversed_radii = sorted(list(radial_plot_data.keys()), reverse=True)
    group_edges = []
    current_index = 0

    for r_val in reversed_radii:
        top_entries = radial_plot_data[r_val]
        for species_name, pct_val in reversed(top_entries):
            y_labels.append(f"{r_val} AU — {species_name}")
            x_percentages.append(pct_val)
            bar_colors.append(species_color_mapping[species_name])
            current_index += 1
        group_edges.append(current_index)

    fig_height = max(4, len(y_labels) * 0.45)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    
    bars = ax.barh(y_labels, x_percentages, color=bar_colors, edgecolor='grey', alpha=0.85, height=0.7)
    
    for bar in bars:
        width = bar.get_width()
        ax.annotate(f' {width:.2f}%',
                    xy=(width, bar.get_y() + bar.get_height() / 2),
                    xytext=(4, 0),
                    textcoords="offset points",
                    ha='left', va='center', fontsize=9, fontweight='bold')

    for edge in group_edges[:-1]:
        ax.axhline(y=edge - 0.5, color='black', linestyle='-', alpha=0.3, linewidth=1.2)

    is_ice_phase = phase in ["surface", "mantle", "grain"]
    bin_title = f" (Bin {grain_bin})" if (grain_bin and is_ice_phase) else (" (All Grains)" if is_ice_phase else "")
    phase_title = f"Phase: {phase.upper()}{bin_title}"

    ax.set_xlabel('Local Radial Budget Contribution (%)', fontsize=12)
    ax.set_ylabel('Disk Radius & Associated Chemical Carrier', fontsize=12)
    
    try:
        first_r = radii[0]
        time_seconds = main_output_dict[radii_map[first_r]]['abundances'].coords['time'].values[itime]
        ax.set_title(f'Top {spnumber} Radial Reservoirs for Element [{target_atom}] — {phase_title}\n$t = {time_seconds/3.156e7:.2e}$ years (Locally Normalized per Shell)', fontsize=13, pad=15)
    except:
        ax.set_title(f'Top {spnumber} Radial Reservoirs for Element [{target_atom}] — {phase_title}', fontsize=13, pad=15)
        
    # Bounds configuration
    max_pct = max(x_percentages) if x_percentages else 100.0
    ax.set_xlim(0, min(max_pct * 1.05,100) + 8.0)  
    ax.set_ylim(-0.6, len(y_labels) - 0.4)
    
    ax.grid(axis='x', linestyle='--', alpha=0.4)

    plt.tight_layout()
    plt.show()


import os
import re
import numpy as np
import matplotlib.pyplot as plt

def plot_disk_atomic_composition(chempath,
                                 main_output_dict,
                                 itime=-1,
                                 verbose=True,
                                 grain_bin=None,
                                 cmap_name="Set1"):
    """
    Computes and plots the total volume-integrated atomic composition of the protoplanetary 
    disk relative to Hydrogen within each specific phase (X_phase / H_phase). 
    Elements are sorted by atomic mass on the X-axis (excluding 'X').
    The Y-axis displays the local abundance ratio in log scale across 5 chemical phases.

    Parameters:
    -----------
    chempath : str
        Path to the directory containing spatial grid subfolders (e.g., '10AU/', '100AU/').
    main_output_dict : dict
        Nested dictionary mapping radial keys to data sub-structures containing 'abundances' 
        and 'H_number_density' spatial profiles.
    itime : int, default -1
        Time index to slice from the multi-epoch data arrays.
    verbose : bool, default True
        If True, outputs warning and missing file notifications to the console.
    grain_bin : int or str, optional
        Specific grain size category to filter when analyzing ice layers. If specified, 
        gas phase contributions are completely ignored for safety.
    cmap_name : str, default "Set1"
        Name of the Matplotlib colormap used to dynamically assign unique colors to the 5 phases.
    """

    # Dictionary of standard atomic masses used to strictly sort the X-axis
    atomic_masses = {
        'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81, 'C': 12.011, 
        'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305, 
        'Al': 26.982, 'Si': 28.085, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 
        'K': 39.098, 'Ca': 40.078, 'Fe': 55.845, 'Ni': 58.693
    }

    # --- INTERNAL HELPERS ---
    def parse_species(species_name):
        if species_name == 'e-':
            return "gas", None, "e-"
            
        grain_match = re.match(r'^([JK])(\d+)(.+)', species_name)
        
        if grain_match:
            p_code, g_bin, raw_formula = grain_match.groups()
            sp_phase = "surface" if p_code == 'J' else "mantle"
            sp_bin = g_bin
        else:
            sp_phase = "gas"
            sp_bin = None
            raw_formula = species_name
        
        clean_formula = raw_formula.replace('c-', '').replace('l-', '')
        return sp_phase, sp_bin, clean_formula

    def get_chemical_composition(clean_formula):
        if clean_formula == 'e-': 
            return {}
        
        calc_formula = clean_formula
        if calc_formula.endswith('+') or calc_formula.endswith('-'):
            calc_formula = calc_formula[:-1]
            
        calc_formula = calc_formula.replace('-', '')
            
        pattern = re.compile(r'([A-Z][a-z]?)(\d*)') 
        composition = {}
        for atom, n in pattern.findall(calc_formula):
            if atom == 'X':
                continue
            try:
                count = int(n) if n else 1
            except ValueError:
                count = 1
            composition[atom] = composition.get(atom, 0) + count
        return composition

    AU_to_cm = 1.496e13  
    phases_pool = ["all", "grain", "mantle", "surface", "gas"]

    # --- 1. RECONSTRUCT RADIAL GRID ---
    radii_map = {}
    for original_key in main_output_dict.keys():
        digits = re.findall(r'\d+', str(original_key))
        if digits:
            try:
                rad_int = int(digits[0])
                radii_map[rad_int] = original_key
            except ValueError:
                continue
                
    radii = sorted(list(radii_map.keys()))
    
    if len(radii) == 0:
        if verbose: print("Main output dictionary is empty or contains no valid numerical keys.")
        return
        
    if len(radii) > 1:
        r_midshifts = 0.5 * np.diff(radii)
        r_edges = [radii[0] - r_midshifts[0]] + [radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [radii[-1] + r_midshifts[-1]]
    else:
        r_edges = [radii[0] * 0.9, radii[0] * 1.1]

    # --- 2. PRE-PARSING NETWORK SPECIES TO BUILD STATIC MAPS ---
    first_r_key = radii_map[radii[0]]
    all_network_species = list(main_output_dict[first_r_key]['abundances'].coords['species'].values)
    
    species_metadata = []
    unique_elements = set()
    
    for idx, sp in enumerate(all_network_species):
        if "GRAIN" in sp:
            species_metadata.append((None, None, {}))
            continue
        sp_phase, sp_bin, clean_formula = parse_species(sp)
        elem_map = get_chemical_composition(clean_formula)
        species_metadata.append((sp_phase, sp_bin, elem_map))
        for elem in elem_map.keys():
            unique_elements.add(elem)
            
    # Sort elements strictly by atomic mass weights
    unique_elements = sorted(list(unique_elements), key=lambda e: atomic_masses.get(e, 999.0))
    num_elements = len(unique_elements)

    # Dictionary to store accumulated disk-wide total absolute atom sums
    disk_total_budgets = {p: {elem: 0.0 for elem in unique_elements} for p in phases_pool}

    # --- 3. DATA ACCUMULATION PER CELL SHELL (VECTORIZED) ---
    for i, r_value in enumerate(radii):
        folder_name = f"{r_value}AU"
        file_path = os.path.join(chempath, folder_name, "1D_static.dat")
        
        if os.path.exists(file_path):
            try:
                z_points = np.loadtxt(file_path, comments='!', usecols=0)
                orig_key = radii_map[r_value]
                sub_dict = main_output_dict[orig_key]
                
                abundance_array = sub_dict['abundances']
                nH_profile = sub_dict["H_number_density"][itime,:]
                
                if len(z_points) > 1:
                    z_midshifts = 0.5 * np.diff(z_points)
                    z_edges = [z_points[0] - z_midshifts[0]] + [z_points[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_points[-1] + z_midshifts[-1])]
                    dz = np.abs(np.diff(z_edges))
                else:
                    dz = np.array([z_points[0] if z_points[0] > 0 else 1.0])
                
                r_left = r_edges[i] * AU_to_cm
                r_right = r_edges[i+1] * AU_to_cm
                dR = r_right - r_left
                R_center = float(r_value) * AU_to_cm
                
                cell_volumes = 2 * np.pi * R_center * dR * (dz * AU_to_cm) 
                
                y_abundances_2d = abundance_array.isel(time=itime).values
                integrated_molecules_per_species = np.sum(y_abundances_2d * nH_profile * cell_volumes, axis=1)
                
                coef_matrices = {p: np.zeros((num_elements, len(all_network_species))) for p in phases_pool}
                
                for sp_idx, (sp_phase, sp_bin, elem_map) in enumerate(species_metadata):
                    if sp_phase is None: 
                        continue
                        
                    if grain_bin is not None:
                        if sp_phase == "gas":
                            continue
                        if sp_bin != str(grain_bin):
                            continue
                    
                    for elem_idx, elem in enumerate(unique_elements):
                        coef = elem_map.get(elem, 0)
                        if coef > 0:
                            coef_matrices[sp_phase][elem_idx, sp_idx] = coef
                            if sp_phase in ["surface", "mantle"]:
                                coef_matrices["grain"][elem_idx, sp_idx] = coef
                            coef_matrices["all"][elem_idx, sp_idx] = coef

                for p in phases_pool:
                    integrated_atoms_per_element = np.dot(coef_matrices[p], integrated_molecules_per_species)
                    for elem_idx, elem in enumerate(unique_elements):
                        disk_total_budgets[p][elem] += integrated_atoms_per_element[elem_idx]

            except Exception as e:
                if verbose: print(f"Error processing cell matrix for R={r_value}: {e}")
        elif verbose: 
            print(f"File not found: {file_path}")

    # --- 4. COLOR ASSIGNMENT VIA COLORMAP ---
    try:
        cmap = plt.get_cmap(cmap_name)
    except ValueError:
        if verbose: print(f"Colormap '{cmap_name}' not found. Falling back to default 'Set1'.")
        cmap = plt.get_cmap("Set1")
        
    phase_colors = {phase_key: cmap(idx) for idx, phase_key in enumerate(phases_pool)}

# --- 5. GRAPH CONSTRUCTION WITH LOCAL PHASE HYDROGEN NORMALIZATION ---
    fig, ax = plt.subplots(figsize=(11, 6))
    x_indexes = np.arange(num_elements)

    # Ajustement dynamique des phases à tracer selon la présence d'un grain_bin
    phases_to_plot = ["grain", "mantle", "surface"] if grain_bin is not None else phases_pool

    for phase_key in phases_to_plot:
        abundance_ratios_to_local_H = []
        
        # Pull out the local Hydrogen reference count for THIS specific phase pool
        local_hydrogen_reference = disk_total_budgets[phase_key].get("H", 0.0)
        
        for elem in unique_elements:
            absolute_phase_val = disk_total_budgets[phase_key][elem]
            
            # Divide the element phase budget by the local phase Hydrogen count
            if local_hydrogen_reference > 0:
                abundance_ratios_to_local_H.append(absolute_phase_val / local_hydrogen_reference)
            else:
                abundance_ratios_to_local_H.append(0.0)
                
        abundance_ratios_to_local_H = np.array(abundance_ratios_to_local_H)
        
        # Mask 0 values to protect logarithmic scale bounds rendering
        valid_mask = abundance_ratios_to_local_H > 0
        
        if np.any(valid_mask):
            ax.plot(x_indexes[valid_mask], abundance_ratios_to_local_H[valid_mask], 
                    label=phase_key.upper(), 
                    color=phase_colors[phase_key], 
                    linewidth=2.5, 
                    marker='o', 
                    markersize=6, 
                    alpha=0.9)

    # Format graph grids and axis layouts
    ax.set_yscale('log')
    ax.set_ylabel('Abundance Ratio (X_phase / H_phase)', fontsize=12)
    ax.set_xlabel('Chemical Elements (Sorted by Atomic Mass Weights)', fontsize=12)
    
    # Tie the sorted elements text tracking sequence to X ticks
    ax.set_xticks(x_indexes)
    ax.set_xticklabels(unique_elements, fontsize=11, fontweight='bold')
    
    bin_suffix = f" (Grain Size Bin {grain_bin})" if grain_bin else ""
    ax.set_title(f"Protoplanetary Disk phase-dependant Elemental Abundances (X/H per phase){bin_suffix}", fontsize=13, fontweight='bold', pad=12)
    
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend(title="Chemical Phases", title_fontsize='11', loc='upper right', frameon=True)

    try:
        first_r = radii[0]
        time_seconds = main_output_dict[radii_map[first_r]]['abundances'].coords['time'].values[itime]
        plt.figtext(0.15, -0.01, f"Data integrated over disk total volume at epoch: $t = {time_seconds/3.156e7:.2e}$ years", 
                    fontsize=10, style='italic')
    except:
        pass

    plt.tight_layout()
    plt.show()

import os
import re
import numpy as np
import matplotlib.pyplot as plt

def plot_species_evolution_with_grain_size(chempath,
                                           main_output_dict,
                                           target_radius,
                                           itime=-1,
                                           verbose=True,
                                           spnumber=5,
                                           cmap_name="tab10"):
    """
    Computes and plots the evolution of the globally dominant ice chemical species 
    (Surface + Mantle combined) as a function of grain size at a specific disk radius.

    Parameters:
    -----------
    chempath : str
        Path to the directory containing spatial grid subfolders (e.g., '10AU/', '100AU/').
    main_output_dict : dict
        Nested dictionary mapping radial keys to data sub-structures.
    target_radius : int or float
        The specific radius (in AU) at which the grain size analysis is performed.
    itime : int, default -1
        Time index to slice from the multi-epoch data arrays.
    verbose : bool, default True
        If True, outputs warning and missing file notifications to the console.
    spnumber : int, default 5
        Number of top global chemical species to display across all bins.
    cmap_name : str, default "tab10"
        Name of the Matplotlib colormap used to dynamically assign unique colors to species.
    """

    # --- INTERNAL HELPERS ---
    def parse_species(species_name):
        if species_name == 'e-':
            return "gas", None, "e-"
        grain_match = re.match(r'^([JK])(\d+)(.+)', species_name)
        if grain_match:
            p_code, g_bin, raw_formula = grain_match.groups()
            sp_phase = "surface" if p_code == 'J' else "mantle"
            # Keep the raw bin identifier name string (e.g. '01', '1', 'fine'...)
            sp_bin = g_bin 
        else:
            sp_phase = "gas"
            sp_bin = None
            raw_formula = species_name
        clean_formula = raw_formula.replace('c-', '').replace('l-', '')
        return sp_phase, sp_bin, clean_formula

    AU_to_cm = 1.496e13  

    # --- 1. FIND THE TARGET RADIUS DATA ---
    radii_map = {}
    for original_key in main_output_dict.keys():
        digits = re.findall(r'\d+', str(original_key))
        if digits:
            try:
                radii_map[int(digits[0])] = original_key
            except ValueError:
                continue
                
    if int(target_radius) not in radii_map:
        if verbose: print(f"Radius {target_radius} AU not found in main_output_dict.")
        return
        
    orig_key = radii_map[int(target_radius)]
    sub_dict = main_output_dict[orig_key]
    abundance_array = sub_dict['abundances']
    nH_profile = sub_dict["H_number_density"][itime,:]
    raw_species_list = list(abundance_array.coords['species'].values)

    # --- 2. LOAD PHYSICAL SPATIAL VOLUMES FOR THIS RADIUS ---
    sorted_all_radii = sorted(list(radii_map.keys()))
    r_idx = sorted_all_radii.index(int(target_radius))
    if len(sorted_all_radii) > 1:
        r_midshifts = 0.5 * np.diff(sorted_all_radii)
        r_edges = [sorted_all_radii[0] - r_midshifts[0]] + [sorted_all_radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [sorted_all_radii[-1] + r_midshifts[-1]]
    else:
        r_edges = [sorted_all_radii[0] * 0.9, sorted_all_radii[0] * 1.1]

    file_path = os.path.join(chempath, f"{int(target_radius)}AU", "1D_static.dat")
    if not os.path.exists(file_path):
        if verbose: print(f"File not found: {file_path}")
        return
        
    z_points = np.loadtxt(file_path, comments='!', usecols=0)
    if len(z_points) > 1:
        z_midshifts = 0.5 * np.diff(z_points)
        z_edges = [z_points[0] - z_midshifts[0]] + [z_points[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_points[-1] + z_midshifts[-1])]
        dz = np.abs(np.diff(z_edges))
    else:
        dz = np.array([z_points[0] if z_points[0] > 0 else 1.0])

    r_left = r_edges[r_idx] * AU_to_cm
    r_right = r_edges[r_idx+1] * AU_to_cm
    cell_volumes = 2 * np.pi * (float(target_radius) * AU_to_cm) * (r_right - r_left) * (dz * AU_to_cm)
    
    nH_2d = nH_profile[np.newaxis, :]
    volumes_2d = cell_volumes[np.newaxis, :]

    # --- 3. HARVEST AND CACHE ALL VALID DUST BINS AND METALISTS ---
    available_bins_set = set()
    species_metadata = {} 

    for sp in raw_species_list:
        if "GRAIN" in sp:
            continue
        sp_phase, sp_bin, clean_formula = parse_species(sp)
        # Isolate only ice matrices (surface + mantle components)
        if sp_phase in ["surface", "mantle"] and sp_bin is not None:
            available_bins_set.add(sp_bin)
            species_metadata[sp] = (sp_phase, sp_bin, clean_formula)

    # Convert to sorted list to keep stable X tracking indexes. 
    # Tries integer sort fallback if labels are numbers as strings.
    try:
        sorted_bins = sorted(list(available_bins_set), key=lambda x: int(x))
    except ValueError:
        sorted_bins = sorted(list(available_bins_set))

    if not sorted_bins:
        print(f"No ice-coated dust grain bins detected at R={target_radius} AU.")
        return

    # Intermediate storage dictionaries
    bin_raw_data = {b: {} for b in sorted_bins}
    global_species_scores = {} 

    # --- 4. DATA ACCUMULATION PER DUST SIZE CLASS ---
    for sp in raw_species_list:
        if sp not in species_metadata:
            continue
        sp_phase, sp_bin, clean_formula = species_metadata[sp]

        # Extract 1D altitude profile array, compute cumulative volume physical mass particles
        y_values = abundance_array.isel(time=itime).sel(species=sp).values
        absolute_particles = np.sum(y_values * nH_2d * volumes_2d)
        
        if absolute_particles > 0:
            # Combine surface + mantle counts for this clean molecule entry inside its respective size bin
            bin_raw_data[sp_bin][clean_formula] = bin_raw_data[sp_bin].get(clean_formula, 0.0) + absolute_particles
            # Track absolute overall molecule abundance across all sizes to determine global dominance ranking
            global_species_scores[clean_formula] = global_species_scores.get(clean_formula, 0.0) + absolute_particles

    if not global_species_scores:
        print(f"No ice molecules detected on dust grains at R={target_radius} AU.")
        return

    # Extract Top N globally significant molecules tracking
    sorted_global_species = sorted(global_species_scores.items(), key=lambda x: x[1], reverse=True)
    top_global_species = [item[0] for item in sorted_global_species[:spnumber]]

    # --- 5. GRAPH CONSTRUCTION ---
    fig, ax = plt.subplots(figsize=(11, 6))
    
    try:
        cmap = plt.get_cmap(cmap_name)
    except ValueError:
        if verbose: print(f"Colormap '{cmap_name}' not found. Falling back to default 'tab10'.")
        cmap = plt.get_cmap("tab10")
        
    species_colors = {sp: cmap(idx / max(1, spnumber - 1)) for idx, sp in enumerate(top_global_species)}

    plot_percentages = {sp: [] for sp in top_global_species}
    x_positions = np.arange(len(sorted_bins))
    # Read the explicit folder/bin string names directly for the X axis labels
    x_ticks_labels = [f"Bin {b}" for b in sorted_bins]

    for b in sorted_bins:
        total_bin_ice_budget = sum(bin_raw_data[b].values())
        for sp in top_global_species:
            if total_bin_ice_budget > 0:
                # Fraction contribution normalized per specific dust bin size shell
                pct = (bin_raw_data[b].get(sp, 0.0) / total_bin_ice_budget) * 100
                plot_percentages[sp].append(pct)
            else:
                plot_percentages[sp].append(0.0)

    # Plot lines tracking profiles
    for sp in top_global_species:
        ax.plot(x_positions, plot_percentages[sp], 
                label=sp, 
                color=species_colors[sp], 
                linewidth=2.5, 
                marker='o', 
                markersize=6)

    # Format frame axis and labels
    ax.set_xlabel('Grain Bins', fontsize=12)
    ax.set_ylabel('Ice Mantle + Surface Molecule Budget Contribution (%)', fontsize=12)
    
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_ticks_labels, fontsize=10, fontweight='bold')
    
    ax.set_ylim(-2, 105)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(title="Top Global Ice Species", title_fontsize='11', loc='upper right',ncol=2)
    
    try:
        time_seconds = abundance_array.coords['time'].values[itime]
        ax.set_title(f'Top {spnumber} Ice Species (Surface + Mantle) vs Grain Size\n$R = {target_radius}$ AU — $t = {time_seconds/3.156e7:.2e}$ years (Locally Normalized)', fontsize=12, fontweight='bold', pad=12)
    except:
        ax.set_title(f'Top {spnumber} Ice Species (Surface + Mantle) vs Grain Size\n$R = {target_radius}$ AU', fontsize=12, fontweight='bold', pad=12)

    plt.tight_layout()
    plt.show()

import os
import re
import numpy as np
import matplotlib.pyplot as plt

import os
import re
import numpy as np
import matplotlib.pyplot as plt

def plot_ratio_midplane_gas_vs_grain(chempath,
                                    main_output_dict,
                                    s1="C",
                                    s2="O",
                                    itime=-1,
                                    starratio=None,
                                    verbose=True,
                                    xlim=None,
                                    ylim=None):
    """Plots the midplane (z=0) atomic abundance ratio of two elements for gas vs. grain phases.

    This function extracts multi-species chemical abundances at the disk midplane 
    across all simulated radii. It classifies species into either gas-phase or grain-phase 
    (ice surface + mantle reservoirs), computes the aggregated elemental ratio (s1/s2) 
    for each phase, and renders a comparative 1D radial line plot.

    Args:
        chempath (str): Path to the directory containing radial subfolders.
        main_output_dict (dict): Nested dictionary containing simulation outputs.
        s1 (str, optional): Atomic symbol for the numerator element. Defaults to "C".
        s2 (str, optional): Atomic symbol for the denominator element. Defaults to "O".
        itime (int, optional): Time index to slice from the abundance arrays. Defaults to -1.
        starratio (float, optional): Ratio s1/s2 of the star. Defaults to None.
        verbose (bool, optional): If True, prints status and error messages. Defaults to True.
        xlim (tuple of float, optional): Manual limits for the horizontal Radius axis.
        ylim (tuple of float, optional): Manual limits for the vertical Ratio axis.

    Raises:
        ValueError: If either s1 or s2 is not included in the allowed chemical network.
    """
    
    # Define valid chemical network elements
    elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    if s1 not in elements or s2 not in elements:
        raise ValueError("One of the specified elements does not exist in the chemical network.")

    # --- INTERNAL SPECIES PARSER AND ELEMENT COUNTER ---
    def parse_and_count(species_name, element1, element2):
        """Identifies the chemical phase and counts target atoms in a species formula."""
        # Ignore electrons and generic structural grain notations
        if species_name == 'e-' or 'GRAIN' in species_name: 
            return "ignore", 0, 0
            
        # Match surface (J) and mantle (K) ice species strings
        grain_match = re.match(r'^([JK])\d+(.+)', species_name)
        if grain_match:
            sp_phase = "grain"  
            raw_formula = grain_match.group(2)
        else:
            sp_phase = "gas"
            raw_formula = species_name
            
        # Clean up structural isomer markers, trailing charge states, and internal hyphens
        clean_formula = raw_formula.replace('c-', '').replace('l-', '')
        if clean_formula.endswith('+') or clean_formula.endswith('-'):
            clean_formula = clean_formula[:-1]
        clean_formula = clean_formula.replace('-', '')
        
        # Regex parsing to extract atomic counts
        pattern = re.compile(r'([A-Z][a-z]?)(\d*)')
        composition = {}
        for atom, n in pattern.findall(clean_formula):
            try:
                count = int(n) if n else 1
            except ValueError:
                count = 1
            composition[atom] = composition.get(atom, 0) + count
            
        return sp_phase, composition.get(element1, 0), composition.get(element2, 0)

    # --- DATA ACCUMULATION ---
    radii_list = []
    ratio_gas_list = []
    ratio_grain_list = []

    # Map string keys to numerical radius values for proper spatial sorting
    radii_map = {}
    for original_key in main_output_dict.keys():
        digits = re.findall(r'\d+', str(original_key))
        if digits:
            radii_map[int(digits[0])] = original_key
            
    sorted_radii_ints = sorted(list(radii_map.keys()))

    # Loop over sorted radial bins to retrieve midplane abundances
    for r_int in sorted_radii_ints:
        orig_key = radii_map[r_int]
        sub_dict = main_output_dict[orig_key]
        abundance_array = sub_dict['abundances']
        
        # Extract abundances at the disk midplane (deepest vertical grid cell -> index -1)
        # Expected array slice format: (species,)
        midplane_abundances = abundance_array.isel(time=itime, spatial=-1).values
        species_list = list(abundance_array.coords['species'].values)
        
        # Reset atomic accumulation pools for this radius
        total_s1_gas, total_s2_gas = 0.0, 0.0
        total_s1_grain, total_s2_grain = 0.0, 0.0
        
        # Aggregate physical atomic pools across all network species
        for idx, species in enumerate(species_list):
            phase, c1, c2 = parse_and_count(species, s1, s2)
            abundance = midplane_abundances[idx]
            
            if phase == "gas":
                total_s1_gas += abundance * c1
                total_s2_gas += abundance * c2
            elif phase == "grain":
                total_s1_grain += abundance * c1
                total_s2_grain += abundance * c2

        # Compute ratios with protective fallback division checks against empty reservoirs
        gas_ratio = total_s1_gas / total_s2_gas if total_s2_gas > 0 else 0.0
        grain_ratio = total_s1_grain / total_s2_grain if total_s2_grain > 0 else 0.0
        
        radii_list.append(float(r_int))
        ratio_gas_list.append(gas_ratio)
        ratio_grain_list.append(grain_ratio)

    if not radii_list:
        if verbose: print("No matching physical grid data could be collected.")
        return

    # Convert python collections to structured numpy arrays
    radii_arr = np.array(radii_list)
    gas_arr = np.array(ratio_gas_list)
    grain_arr = np.array(ratio_grain_list)

    # --- PLOTTING ---
    fig, ax = plt.subplots(figsize=(9, 5))
    
    # Apply a logarithmic scale to account for sharp chemical variations
    ax.set_yscale('log')
    
    # Render line tracks + markers for each discrete phase
    ax.plot(radii_arr, gas_arr, color="teal", linestyle='-', marker='o', label='Gas', linewidth=1.8)
    ax.plot(radii_arr, grain_arr, color="darkred", linestyle='--', marker='s', label='Grains (Ice)', linewidth=1.8)

    # Label definitions
    ax.set_xlabel('Radius R [AU]', fontsize=11)
    ax.set_ylabel(f'Midplane Atomic Ratio [{s1}/{s2}]', fontsize=11)
    
    # Safely parse coordinate timestamps for dynamic title labeling
    try:
        first_key = radii_map[sorted_radii_ints[0]]
        time_seconds = main_output_dict[first_key]['abundances'].coords['time'].values[itime]
        ax.set_title(f'Atomic Ratio {s1}/{s2} at Midplane ($z=0$) — $t = {time_seconds/3.156e7:.2e}$ yr', fontsize=12, pad=12)
    except:
        ax.set_title(f'Atomic Ratio {s1}/{s2} at Midplane ($z=0$)', fontsize=12, pad=12)

    if starratio is not None:
        ax.axhline(y=starratio,color='black',linestyle='--',label=f'Star {s1}/{s2} ratio')

    # Apply manual axis boundaries if supplied
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)
    
    # Auto-scale vertical bounds while avoiding log(0) clipping exceptions
    if ylim is None:
        all_vals = np.concatenate([gas_arr, grain_arr])
        positive_vals = all_vals[all_vals > 0]
        if len(positive_vals) > 0:
            ax.set_ylim(positive_vals.min() * 0.5, positive_vals.max() * 2)
        else:
            ax.set_ylim(1e-4, 1e2)

    # Render grid lines and legendary context block
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(frameon=True, facecolor='white', edgecolor='gainsboro', loc='best')
    
    plt.tight_layout()
    plt.show()
