#!/usr/local/bin/python
# Time-stamp: <2016-02-10 09:59:24 marine>
# Project : IC Dynamics
# Subproject : read data and plot data
# Author : Marine Lasbleis

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile


def import_data_G(name="G_0.01965", folder_name="./OUT/"):

    f = FortranFile(folder_name+name, 'r')

    version = f.read_reals(dtype='f4')
    time, Ra, Ra_c, P, Ha, Di, Pr, Le = f.read_reals(dtype='f4')
    nradius, ntheta, nphi, azsym = f.read_reals(dtype='f4') # be careful, all of them are reals
    radius = f.read_reals(dtype='f4')
    theta = f.read_reals(dtype='f4') #colatitude

    phi = np.arange(1, int(nphi)+1)/nphi*2.*np.pi/azsym #longitude (not read from file!)

    Vr = np.empty([nphi, ntheta, nradius])
    Vt = np.empty_like(Vr)
    Vp = np.empty_like(Vr)
    Temperature = np.empty_like(Vr)
    Composition = np.empty_like(Vr)

    for ir in np.arange(nradius):
        for it in np.arange(ntheta):
            Vr[:,it,ir]=f.read_reals(dtype='f4')
            Vt[:,it,ir]=f.read_reals(dtype='f4')
            Vp[:,it,ir]=f.read_reals(dtype='f4')
            Composition[:,it,ir]=f.read_reals(dtype='f4')
            Temperature[:,it,ir]=f.read_reals(dtype='f4')
    
    return time, Ra, Ra_c, P, Ha, Di, Pr, Le, nradius, ntheta, nphi, azsym, radius, theta, phi, Vr, Vt, Vp, Temperature, Composition


def average_velocity(Vr, Vt, Vp):
    """ average velocity """
    return np.sqrt(Vr**2.+Vt**2.+Vp**2.)


def cartesian_velocities(Vr, Vt, Vp, theta, phi):
    ## TODO: here, Vp is considered zero, and Vx is zero.

    Vx = np.empty_like(Vr)
    Vy = np.empty_like(Vr)
    Vz = np.empty_like(Vr)
    
    for ip, phi_ in enumerate(phi):
        for it, theta_ in enumerate(theta):
            Vy[ip, it, :] = (Vr[ip, it, :] *np.sin(theta_) + Vt[ip, it, :] *np.cos(theta_))*np.cos(phi_)
            Vz[ip, it, :] = Vr[ip, it, :] *np.cos(theta_) - Vt[ip, it, :] *np.sin(theta_)

    return Vx, Vy, Vz


def average_radius(quantity, theta):
    """ average value of the quantity, for each radius

    INPUT:
    quantity: a 3D numpy-array (longitude, colatitude, radius)
    theta: a 1D array with values of colatitude
    OUTPUT:
    """
    
    nradius, ntheta, nphi = quantity.shape

    THETA = np.tile(np.array([theta]).T, (1, nphi)) #create copies of theta along the dimension of the radius
    THETA = np.tile(THETA, (nradius, 1, 1)) # create copies of theta along the dimension of the longitude
    
    quantity = quantity * np.sin(THETA)

    return np.sum(np.sum(quantity, axis=0), axis=0)/ np.sum(np.sum(np.sin(THETA), axis=0), axis=0)


def average_global(quantity, radius, theta):
    """ mean value over the whole volume """

    THETA = np.tile(np.array([theta]).T, (1, nphi)) #create copies of theta along the dimension of the radius
    THETA = np.tile(THETA, (nradius, 1, 1)) # create copies of theta along the dimension of the longitude
    
    quantity_radius = average_radius(quantity, theta)
    dr = np.diff(radius)
    average_global = np.sum(dr*0.5*(quantity_radius[:-1]*radius[:-1]**2+quantity_radius[1:]*radius[1:]**2))
    volume = np.sum(dr*0.5*(radius[:-1]**2+radius[1:]**2))


    return average_global/volume


if __name__ == '__main__':

    
    time, Ra, Ra_c, P, Ha, Di, Pr, Le, nradius, ntheta, nphi, azsym, radius, theta, phi, Vr, Vt, Vp, Temperature, Composition = import_data_G(name="G_0.02456")

    average_T = average_radius(Temperature, theta)
    
    print 'T0: ', average_global(Temperature, radius, theta)
    print 'Vr rms: ', np.sqrt(average_global(Vr**2., radius, theta))
    print 'Vh rms: ', np.sqrt(average_global(Vt**2.+Vp**2., radius, theta))
    print 'V rms: ', np.sqrt(average_global(average_velocity(Vr, Vt, Vp)**2., radius, theta))
