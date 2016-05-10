#!/usr/local/bin/python
# Time-stamp: <2016-02-18 11:06:34 marine>
# Project : IC Dynamics
# Subproject : read data and plot data
# Author : Marine Lasbleis

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

# personal routines
import IC_plot



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


def average_global(quantity, radius, theta, nphi):
    """ mean value over the whole volume """
 

    nradius = len(radius)
    THETA = np.tile(np.array([theta]).T, (1, nphi)) #create copies of theta along the dimension of the radius
    THETA = np.tile(THETA, (nradius, 1, 1)) # create copies of theta along the dimension of the longitude
    
    quantity_radius = average_radius(quantity, theta)
    dr = np.diff(radius)
    average = np.sum(dr*0.5*(quantity_radius[:-1]*radius[:-1]**2+quantity_radius[1:]*radius[1:]**2))
    volume = np.sum(dr*0.5*(radius[:-1]**2+radius[1:]**2))


    return average/volume


def crossection_data(data, coordinate, choice, sign=1., i=1):
    """ Reshape the data to have field of a quantity `data` to plot over either meridional or equatorial cross section.

    INPUT:
    - data: 3D array (phi, theta, r) (theta is colatitude)
    - coordinate: 1D array (phi or theta).
    - sign needs to be positive (1) if Vr, and negative (-1) if Vt
    OUTPUT:
    - data_final: 2D array (phi/theta, r)
    - coordinate: 1D array (phi or theta, but reshape to fullfilled the whole disc)

    If meridional (and theta) is used, data is reshape to have theta from 0 to 2 pi (input should be theta from 0 to pi).
    If equatorial (and phi) is used, the value for phi[-1] is concatenated to the beginning of the array, to fullfill the disc also. 
     """
    
    iphi, itheta, iradius = data.shape
    if choice == "meridional":
        # in this case, coordinate is theta, and i is the i_phi. Output are functions of (theta, r)
        data_sq1 = np.squeeze(data[i,:,:])
        data_sq2 = sign*np.squeeze(data[iphi/2,:,:])
        
        data_final = np.concatenate((np.array([0.5*(data_sq1[0,:]+data_sq2[0,:])]),
                               data_sq1,
                               np.array([0.5*(data_sq1[-1,:]+data_sq2[-1,:])]),
                               np.flipud(data_sq2),
                               np.array([0.5*(data_sq1[0,:]+data_sq2[0,:])])))
        coordinate = np.concatenate((np.array([0]), coordinate,np.array([np.pi]), np.pi+ coordinate, np.array([2.*np.pi])))

        
    if choice == "equatorial": # in this case, coordinate is phi, and i is not used. Outputs are functions of (phi, r)
        data_sq = np.squeeze(data[:,itheta/2,:])
        data_final = np.concatenate((np.array([data_sq[-1, :]]), data_sq))
        coordinate = np.concatenate((np.array([0]), coordinate))
        

    return data_final, coordinate


def vorticity_phi(Vradius, Vtheta, radius, theta):
    """ compute the vorticity field in the meridional cross section
    
    with the assumptions of cylindrical symmetry (V_\phi = 0 and all \partial/\partial_\phi =0),
    this gives directly the vorticity field.

    INPUT:
    - Vradius, Vtheta:  2D arrays produced by `crossection_data()`. Same size (n_t, n_r)
    - radius: 1D array with radius (size n_r)
    - theta: 1D array, as produced by `crossection_data()`. Size n_t
    OUPUT:
    - vort: vorticity field. 2D array of size (n_t, n_r)
    """

    dr = np.gradient(radius)
    dtheta = np.gradient(theta)
    drVtdr = np.gradient(radius*Vtheta, np.array([dr]))[1]
    dVrdt = np.gradient(Vradius, np.array([dtheta]).T)[0]
    vort = np.empty_like(Vradius)
    vort[:, 1:] = -1./radius[1:] * ( drVtdr[:, 1:] - dVrdt[:, 1:] ) #radius[0]=0, so vort[0]=0
    return vort

def vorticity_theta(Vradius, Vphi, radius, phi, theta=np.pi/2.):
    """ compute the vorticity field in the equatorial cross section
    
    with the assumptions of cylindrical symmetry (V_\theta = 0 and all \partial/\partial_\theta =0),
    this gives directly the vorticity field.
    By default, theta = pi/2 (equatorial plane) and sin(theta)=1.

    INPUT:
    - Vradius, Vphi:  2D arrays produced by `crossection_data()`. Same size (n_p, n_r)
    - radius: 1D array with radius (size n_r)
    - phi: 1D array, as produced by `crossection_data()`. Size n_p
    OUPUT:
    - vort: vorticity field (2D array of size (n_p, n_r)
    """

    dr = np.gradient(radius)
    dphi = np.gradient(phi)
    drVpdr = np.gradient(radius*Vphi, dr)[1]
    dVrdp = np.gradient(Vradius, np.array([dphi]).T)[0]
    vort = np.empty_like(Vradius)
    vort[:,1:] = -1./radius[1:] * (dVrdp[:, 1:]/np.sin(theta) - drVpdr[:, 1:])#radius[0]=0, so vort[0]=0

    return vort



def tests(filename):

    
    time, Ra, Ra_c, P, Ha, Di, Pr, Le, nradius, ntheta, nphi, azsym, radius, theta, phi, Vr, Vt, Vp, Temperature, Composition = import_data_G(name=filename)

    average_T = average_radius(Temperature, theta)

    print Vt.shape
    print nphi
    print nradius, len(radius)

    print 'T0: ', average_global(Temperature, radius, theta, nphi)
    print 'Vr rms: ', np.sqrt(average_global(Vr**2., radius, theta, nphi))
    print 'Vh rms: ', np.sqrt(average_global(Vt**2.+Vp**2., radius, theta, nphi))
    print 'V rms: ', np.sqrt(average_global(average_velocity(Vr, Vt, Vp)**2., radius, theta, nphi))
    
    Vr_me, _ = crossection_data(Vr, theta, choice='meridional', sign=1, i=0)
    Vt_me, theta_total = crossection_data(Vt, theta, choice='meridional', sign=-1,  i=0)

    vorticity_me = vorticity_phi(Vr_me, Vt_me, radius, theta_total)


    Vr_eq, _ = crossection_data(Vr, phi, choice='equatorial', sign=1, i=0)
    Vp_eq, phi_total = crossection_data(Vp, phi, choice='equatorial', sign=-1,  i=0)

    vorticity_eq = vorticity_theta(Vr_eq, Vp_eq, radius, phi_total)


    print vorticity_me.shape, radius.shape, theta_total.shape, phi_total.shape


    #TODO : set the two on same figure
    Temperature_me, theta_total = crossection_data(Temperature, theta, choice='meridional', sign=1, i=0)
    Composition_me, theta_total = crossection_data(Composition, theta, choice='meridional', sign=1, i=0)

    fig, ax = plt.subplots(1)
    IC_plot.NS_cross_section(Temperature_me, theta_total, radius, label="Temperature", fig_info=(fig, ax))
    IC_plot.NS_quiver_plot(Vr_me, Vt_me, theta_total, radius, label="Temperature and velocity", fig_info=(fig, ax))
    
    fig2, ax2 = plt.subplots(1)
    IC_plot.NS_cross_section(Composition_me, theta_total, radius, label="Composition", fig_info=(fig2, ax2))
    IC_plot.NS_quiver_plot(Vr_me, Vt_me, theta_total, radius, label="Composition and velocity", fig_info=(fig2, ax2))
   
    plt.show()



if __name__ == '__main__':


    tests("G_0.25527")

