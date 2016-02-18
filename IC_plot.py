#!/usr/local/bin/python
# Time-stamp: <2016-02-18 11:07:07 marine>
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


def NS_cross_section(quantity, t, r, fig_info, label="test"):
    """ meridional cross section. """

    nb_contour = 30
    ntheta, nradius = quantity.shape

    RADIUS = np.tile(r, (ntheta, 1))
    THETA = np.tile(np.array([t]).T, (1, nradius))

    X = RADIUS * np.sin(THETA)
    Y = RADIUS * np.cos(THETA)


    CS = fig_info[1].contourf(X, Y, quantity, nb_contour, cmap=plt.cm.hot)
    fig_info[1].axis("equal")
    fig_info[1].set_title(label+" NS section")
    plt.colorbar(CS)
    fig_info[1].get_xaxis().set_visible(False)
    fig_info[1].get_yaxis().set_visible(False)
    
    return

def NS_quiver_plot(Vr, Vt, t, r, fig_info, label='test', sampling=5):

    
    ntheta, nradius = Vr.shape
    RADIUS = np.tile(r, (ntheta, 1))
    THETA = np.tile(np.array([t]).T, (1, nradius))

    X = RADIUS * np.sin(THETA)
    Y = RADIUS * np.cos(THETA)

    Vy = np.sin(THETA)*Vr + np.cos(THETA)*Vt
    Vz = np.cos(THETA)*Vr - np.sin(THETA)*Vt #TODO : verify if we need to change the theta for half of the profile

    fig_info[1].quiver(X[:,::sampling], Y[:,::sampling], Vy[:,::sampling], Vz[:,::sampling], units='x')
    fig_info[1].axis("equal")
    fig_info[1].set_title(label+" NS section")
    fig_info[1].get_xaxis().set_visible(False)
    fig_info[1].get_yaxis().set_visible(False)



    

     
    
    
    

    


if __name__ == '__main__':

    plot_radius()
    plot_evolution()
