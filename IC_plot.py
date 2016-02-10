#!/usr/local/bin/python
# Time-stamp: <2016-02-10 10:00:13 marine>
# Project : IC Dynamics
# Subproject : plot data
# Author : Marine Lasbleis

import numpy as np
import matplotlib.pyplot as plt


def plot_evolution(folder_name="./OUT/", evolution="evolution.txt"):
    """ Read the file evolution from Parody_IC and plot some variables as function of time.

    default file is "./OUT/evolution.txt" (folder_name+evolution)
    plot as function of time: Ra, P, melting and energVP
    """

    data_evolution = np.genfromtxt(folder_name+evolution)
    time = data_evolution[:,0]
    Ra = data_evolution[:,1]
    P = data_evolution[:,2]
    Vpr_mean = data_evolution[:,3]
    Vpr_imag = data_evolution[:,4]
    Vpr_real = data_evolution[:,5]
    melting = data_evolution[:,6]
    energVP = data_evolution[:,7]
    energVP_r = data_evolution[:,8]
    Theta_mean = data_evolution[:,9]
    Theta_rms = data_evolution[:,10]

    fig, ax = plt.subplots(2,2)

    ax[0,0].plot(time, Ra)
    ax[0,0].set_title('Ra')
    ax[0,1].plot(time, P)
    ax[0,1].set_title('P')
    ax[1,0].plot(time, melting)
    ax[1,0].set_title('melting')
    ax[1,1].plot(time, energVP)
    ax[1,1].set_title('energVP')

    plt.show()



def plot_radius(folder_name="./OUT/", radius="rayon"):
    """ plot the repartition of points in radius from the simulation (constant with time)"""
    
    data = np.genfromtxt(folder_name+radius)
    plt.plot(data[:,1], data[:,2], '.')
    plt.title('dr as function of radius')
    plt.show()



if __name__ == '__main__':

    plot_radius()
    plot_evolution()
