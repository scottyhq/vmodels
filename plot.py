# -*- coding: utf-8 -*-
"""
Plotting files to accompany vmodels (mogi, yang, okada)
Created on Tue Jun 21 14:59:15 2016

NOTE: Currently set up for plotting Okada

@author: scott
"""
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd

# Add some generic plotting functions for re-use?


def plot_fault(axes, strike=None, dip=None, length=None, width=None, xcen=None, ycen=None, **kwargs):
    ''' Project fault plane onto surface'''
    # Project fault coordinates onto surface
    sins = np.sin(np.deg2rad(strike))
    coss = np.cos(np.deg2rad(strike))
    Lx = 0.5 * length * sins
    Ly = 0.5 * length * coss
    W = width * np.cos(np.deg2rad(dip))

    # Concatenate coordinates
    xb = np.array([-Lx + W * coss, -Lx, Lx, Lx + W * coss, -Lx + W * coss]) + xcen
    yb = np.array([-Ly - W * sins, -Ly, Ly, Ly - W * sins, -Ly - W * sins]) + ycen

    # scale for plotting (km)
    xb = xb * 1e-3
    yb = yb * 1e-3

    end1x = (sins * length/2 * 1e-3)
    end1y = (coss * length/2 * 1e-3)

    # put it on the plots!
    for ax in axes:
        ax.plot(xb, yb, 'k-', lw=2, scalex=False, scaley=False)
        ax.plot(end1x,end1y,'ko',mfc='none', scalex=False, scaley=False) #stirke direction indicator
        #ax.plot([xcen, ], [ycen, ])# add tick at midpoint indicating dip
        # multiple distance times direction perpendicular to stike

def plot_fault3d(axes, strike=None, dip=None, length=None, width=None, xcen=None, ycen=None, **kwargs):
    ''' Rotatable fault in 3d'''
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    print('todo')


def plot_los_indicator(ax, ald):
    ''' Add LOS arrow indicator in axes coordinates
    Inputs:
        ax     axes to add to
        ald    azimuth look direction (second array in geo_incidence
    '''
    L = 0.1
    x0, y0 = (0.8, 0.8)
    dx = L * np.cos(np.pi / 2 - np.deg2rad(ald))
    dy = L * np.sin(np.pi / 2 - np.deg2rad(ald))
    ax.arrow(x0, y0, dx, dy, transform=ax.transAxes, color='k')  # add text too:
    ax.text(0.9, 0.9, 'LOS', ha='right', va='top',
            transform=ax.transAxes, fontweight='bold')
    # ax1.annotate('LOS', (x0,y0), xytext=(x0+dx,x0+dy), xycoords='axes fraction', textcoords='axes fraction',
    # arrowprops=dict(width=1,frac=0.3,headwidth=5,facecolor='black'),
    # fontweight='bold') # NOTE: imshow origin has effect on this

def plot_data_model_residual(data, model, extent, row, northing,
                             vmin=None, vmax=None, unit='cm', point=None,
                             xlim=None, ylim=None, scale_residual=False):
    '''
    Data. Model, Residual with EW profile
    '''
    cmap = 'bwr'
    extent = np.array(extent)/1e3
    northing = northing/1e3

    fig, axes = plt.subplots(2, 2, figsize=(10,10))
    ax = axes.flat[0]
    im = ax.imshow(data, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.axhline(northing, color='k')
    ax.set_title('Data [{}]'.format(unit))
    cb = plt.colorbar(im, ax=ax, orientation='vertical')
    plt.grid(True)

    ax = axes.flat[1]
    im = ax.imshow(model, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.axhline(northing, color='m')
    ax.set_title('Model [{}]'.format(unit))
    cb = plt.colorbar(im, ax=ax, orientation='vertical')
    plt.grid(True)

    ax = axes.flat[2]
    if scale_residual:
        im = ax.imshow(data-model, extent=extent, cmap=cmap) #Autoscale Residual
    else:
        im = ax.imshow(data-model, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_title('Residual [{}]'.format(unit))
    ax.set_xlabel('Easting [km]')
    ax.set_ylabel('Northing [km]')
    cb = plt.colorbar(im, ax=ax, orientation='vertical')
    plt.grid(True)

    # Label Summit
    for ax in axes.flat[:3]:
        if point:
            ax.plot(point[0]/1e3,point[1]/1e3,'k^',scalex=False,scaley=False)
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    ax = axes.flat[3]
    eastings = np.linspace(extent[0], extent[1], data.shape[1])
    ax.plot(eastings, data[row],'ko', mfc='none', label='data')
    ax.plot(eastings, model[row],'m-', lw=2, label='model')
    ax.set_xlabel('Easting [m]')
    ax.axhline(0, color='k')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.savefig('data_model_residual.pdf', bbox_inches='tight')
    plt.show()

def get_clim(data):
    vmax = np.nanmax(np.abs(data))
    vmin = -vmax
    return vmin, vmax

def plot_components(x, y, ux, uy, uz, params=None,vmin=None, vmax=None):
    '''
    show components of deformation, along with projection into LOS
    NOTE: also would be cool to plot 3D surface!
    # just pass all parameters
    '''
    cmap = 'bwr'
    # Convert to km and cm for ploting
    x, y = np.array([x, y]) * 1e-3
    ux, uy, uz = np.array([ux, uy, uz]) * 1e2

    # step size for quiver plot resampling
    nx = 20
    ny = 20

    fig, (ax, ax1, ax2) = plt.subplots(1, 3,
                                subplot_kw=dict(aspect=1.0, adjustable='box-forced'),
                                sharex=True, sharey=True, figsize=(17, 6))  # fig_kw

    extent = [x.min(), x.max(), y.min(), y.max()] #left, right, bottom, top

    #normalize colorscale
    vmin, vmax = get_clim(uz)
    im = ax.imshow(uz, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title('Vertical Displacement, Uz')
    ax.set_xlabel('EW Distance [km]')
    ax.set_ylabel('NS Distance [km]')
    cb = plt.colorbar(im, ax=ax, orientation='horizontal')
    cb.set_label('cm')

    vmin, vmax = get_clim(ux)
    im1 = ax1.imshow(ux, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    #ax1.quiver(x[::nx, ::ny], y[::nx, ::ny], ux[::nx, ::ny], np.zeros_like(uy)[::nx, ::ny])
    ax1.set_title('EW Displacement, Ux')
    cb1 = plt.colorbar(im1, ax=ax1, orientation='horizontal')  # , pad=0.1)
    cb1.set_label('cm')

    vmin, vmax = get_clim(uy)
    im2 = ax2.imshow(uy, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    #ax2.quiver(x[::nx, ::ny], y[::nx, ::ny], np.zeros_like(ux)[::nx, ::ny], uy[::nx, ::ny])
    cb2 = plt.colorbar(im2, ax=ax2, orientation='horizontal')  # , pad=0.1)
    ax2.set_title('NS Displacement, Uy')
    cb2.set_label('cm')

    # Add crosshair grid:
    for a in (ax, ax1, ax2):
        a.axhline(linestyle=':', color='k')
        a.axvline(linestyle=':', color='k')

    if params:
        plot_fault((ax, ax1, ax2), **params)

    plt.suptitle('Components of Deformation', fontsize=16, fontweight='bold')


def plot_los(x, y, ux, uy, uz, los, params,cmap='jet'):
    ''' Separate figure showing displacement and Wrapped Phase in Radar     '''
    # Convert to km and cm for ploting
    x,y = np.array([x,y]) * 1e-3
    ux,uy,uz,los = np.array([ux,uy,uz,los]) * 1e2

    # extract a few varialbes from params dictionary
    inc = params['inc']
    ald = params['ald']
    wavelength = params['wavelength']

    # Set view & reduce number of quiver arrows
    nx=10
    ny=10
    extent = [x.min(), x.max(), y.min(), y.max()]

    # --------------
    fig, (ax,ax1,ax2) = plt.subplots(1,3,subplot_kw=dict(aspect=1.0, adjustable='box-forced'),sharex=True, sharey=True, figsize=(17,6))
    # vertical displacement w/ horizontal vectors # NOTE: add quiver scale arrow!
    im = ax.imshow(uz, extent=extent, cmap=cmap)
    #ax1.quiver(x[::nx, ::ny], y[::nx, ::ny], ux[::nx, ::ny], np.zeros_like(uy)[::nx, ::ny]) #From above
    ax.quiver(x[::nx,::ny], y[::nx,::ny], ux[::nx,::ny], uy[::nx,::ny])
    ax.set_title('Model Displacement Field')
    ax.set_xlabel('EW Distance [km]')
    ax.set_ylabel('NS Distance [km]')
    cb = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.1) #ticks=MaxNLocator(nbins=5) #5 ticks only)
    cb.set_label('cm')

    im = ax1.imshow(los, extent=extent, cmap=cmap)
    # ax1.quiver(x[::nx,::ny], y[::nx,::ny], ux[::nx,::ny], uy[::nx,::ny])
    plot_los_indicator(ax1,ald)
    ax1.set_title('LOS (inc,ald)=(%s,%s)' % (inc,ald) )
    cb = plt.colorbar(im, ax=ax1, orientation='horizontal', pad=0.1) #ticks=MaxNLocator(nbins=5) #5 ticks only)
    cb.set_label('cm')

    # wrapped LOS - think about why this works...
    los_wrapped = np.remainder(los - los.min(), wavelength/2.0) / (wavelength/2.0)
    im = ax2.imshow(los_wrapped, extent=extent, cmap='jet', vmin=0, vmax=1)
    plot_los_indicator(ax2,ald)
    ax2.set_title('LOS Wrapped')
    cb = plt.colorbar(im, ax=ax2, orientation='horizontal', ticks=[], pad=0.1) #ticks=MaxNLocator(nbins=5) #5 ticks only)
    cb.set_label('{0} cm'.format(wavelength/2)) #1 color cycle = lambda/2

    axes = (ax,ax1,ax2)
    plot_fault(axes, **params)
    # Add crosshair grid:
    for a in axes:
        a.axhline(linestyle=':', color='k')
        a.axvline(linestyle=':', color='k')

    plt.suptitle('Model in Radar Line of Sight', fontsize=16, fontweight='bold')


def plot_profile(x,y,ind,ux,uy,uz,los,axis=0,extent=None):
    ''' straight line profile through specified axis of ux,uy,uz,los'''
    # Show profile line in separate LOS plot
    # Convert to km and cm for ploting
    x,y = np.array([x,y]) * 1e-3
    ux,uy,uz,los = np.array([ux,uy,uz,los]) * 1e2

    extent = [x.min(), x.max(), y.min(), y.max()]
    plt.figure()
    plt.imshow(los, cmap='bwr',extent=extent)
    if axis == 0:
        plt.axhline(y[ind],color='k', lw=2)
        # plt.annotate('A',(ind,0))
    else:
        print(x[ind])
        plt.axvline(x[ind],color='k', lw=2)
    plt.title('Profile line')


    # Convert to km and cm for ploting
    # x,y = np.array([x,y]) * 1e-3
    #ux,uy,uz,los = np.array([ux,uy,uz,los]) * 1e2
    dist = x

    if axis==1:
        dist = y
        ux,uy,uz,los = [x.T for x in [ux,uy,uz,los]]

    # extract profiles
    ux_p = ux[ind]
    uy_p = uy[ind]
    uz_p = uz[ind]
    los_p = los[ind]

    fig, axes = plt.subplots(2,2,sharex=True, sharey=True, figsize=(17,6))
    ax = axes.flat[0]
    ax.plot(dist, ux_p,'k.-',lw=2)
    ax.set_title('Ux')
    ax.set_ylabel('Displacement [cm]')
    ax.set_xlabel('Distance [km]')
    #ax.text(0,0.9,'A',fontweight='bold',ma='center',transform=ax.transAxes)
    #ax.text(0.9,0.9,'B',fontweight='bold',ma='center',transform=ax.transAxes)

    ax = axes.flat[1]
    ax.plot(dist, uy_p,'k.-',lw=2)
    ax.set_title('Uy')
    ax = axes.flat[2]
    ax.plot(dist, uz_p,'k.-',lw=2)
    ax.set_title('Uz')
    ax = axes.flat[3]
    ax.plot(dist, los_p,'k.-',lw=2)
    ax.set_title('LOS')

    for ax in axes.flat:
        ax.grid(True)

    # ratio of uz to ur
    # NOTE: just print max ratio
    print('Ur_max/Uz_max (along profile)= ', ux_p.max() / uz_p.max())
    '''
    plt.figure()
    plt.plot(uz_p/ux_p,'k.-',lw=2)
    plt.xlabel('Distance [km]')
    plt.ylabel('Uz/Ur Ratio')
    plt.title('Vertical vs. Horizontal Displacements')
    plt.show()
    '''

#From Yang
'''
def plot_los(x,y,los):
    plt.figure()
    plt.imshow(los.reshape(x.shape), origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
    plt.colorbar()
    plt.title('Yang LOS [cm]')
    plt.xlabel('EW Distance [km]')
    plt.ylabel('NS Distance [km')
    plt.show()
'''

def plot_yang(x,y,ux,uy,uz):
    # Reshape output and subsample for plotting quiver
    '''
    dhx=np.round(ndat / 10)
    dhy=np.round(mdat / 10)
    xsub=x[:-1:dhx,:-1:dhy]
    ysub=y[:-1:dhx,:-1:dhy]
    uxsub = ux[0:-1:dhx,0:-1:dhy]
    uysub = uy[0:-1:dhx,0:-1:dhy]
    '''

    # Simpler: Just take every ten points
    #xsub = xsub[:,:,10]
    nx = 10
    ny = 10

    # Plot model
    #plt.figure()
    #plt.pcolor(yz)
    plt.imshow(uz, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
    #plt.imshow(uz, origin='lower')
    cb = plt.colorbar()
    cb.set_label('cm')
    plt.quiver(x[::ny,::ny], y[::nx,::ny], ux[::nx,::ny], uy[::nx,::ny])
    plt.title('Yang Uz & Ur Surface Deformation')
    plt.xlabel('EW Distance [km]')
    plt.ylabel('NS Distance [km')
