# -*- coding: utf-8 -*-
"""
Plotting files to accompany vmodels (mogi, yang, okada)
Created on Tue Jun 21 14:59:15 2016

@author: scott
"""
import matplotlib as plt
import numpy as np



def get_cart2los(inc, ald, x):
    '''
    NOTE: possible sign convention issues with this function
    '''
    # x is data array
    # converted to LOS
    # los = data.dot(cart2los) * 1e2 # maybe use numpy.tensordot? not sure...
    # For now fake it:
    look = np.deg2rad(inc) * np.ones_like(x)  # incidence
    # heading (degreees clockwise from north)
    head = np.deg2rad(ald) * np.ones_like(x)
    # NOTE: matlab code default is -167 'clockwise from north' (same as hannsen text fig 5.1)
    # This is for descending envisat beam 2, asizmuth look direction (ALD) is
    # perpendicular to heading (-77)

    # however, make_los.pl generates unw file with [Incidence, ALD], ALD for ascending data is 77
    # make_los.pl defines "(alpha) azimuth pointing of s/c"
    EW2los = np.sin(head) * np.sin(look)
    NS2los = np.cos(head) * np.sin(look)
    Z2los = -np.cos(look)
    # NOTE: negative here implies uplift=positive in LOS
    cart2los = -np.dstack([EW2los, NS2los, Z2los])

    return cart2los



def plot_fault(fig, strike=None, delta=None, length=None, width=None, xcen=None, ycen=None, **kwargs):
    ''' matlab way to project fault plane onto surface'''
    # XB = [] #lists for multiple faults in same domain
    # YB = []

    # Project fault coordinates onto surface
    sins = np.sin(np.deg2rad(strike))
    coss = np.cos(np.deg2rad(strike))
    Lx = 0.5 * length * sins
    Ly = 0.5 * length * coss
    W = width * np.cos(np.deg2rad(delta))

    # Concatenate coordinates
    xb = np.array([-Lx + W * coss, -Lx, Lx, Lx + W * coss, -Lx + W * coss]) + xcen
    yb = np.array([-Ly - W * sins, -Ly, Ly, Ly - W * sins, -Ly - W * sins]) + ycen
    # XB.append(xb)
    # YB.append(yb)

    # scale for plotting
    xb = xb * 1e-3
    yb = yb * 1e-3

    # put it on the plots!
    for ax in fig.get_axes():
        ax.plot(xb, yb, 'w-', lw=2)


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


def plot_components(x, y, ux, uy, uz, los, params, profile=False):
    '''
    show components of deformation, along with projection into LOS
    NOTE: also would be cool to plot 3D surface!
    # just pass all parameters
    '''
    cmap = 'bwr'
    # Convert to km and cm for ploting
    x, y = np.array([x, y]) * 1e-3
    ux, uy, uz, los = np.array([ux, uy, uz, los]) * 1e2

    # step size for quiver plot resampling
    nx = 20
    ny = 20

    fig, (ax, ax1, ax2) = plt.subplots(1, 3,
                                subplot_kw=dict(aspect=1.0, adjustable='box-forced'),
                                sharex=True, sharey=True, figsize=(17, 6))  # fig_kw

    extent = [x.min(), x.max(), y.min(), y.max()]
    # plt.locator_params(tight=True, nbins=4) #control number of easting/northing ticks
    # sc = ax.scatter(x_km,y_km,c=data,cmap=plt.cm.bwr) #colormap not centered on zero
    # norm = MidpointNormalize(midpoint=0)
    # sc = ax.scatter(x_km,y_km,c=data,cmap=plt.cm.bwr,norm=norm)
    # im = ax.imshow(uz)
    # im = ax.pcolor(x,y,uz)
    im = ax.imshow(uz, extent=extent, cmap=cmap)
    # ax.quiver(x[::nx,::ny], y[::nx,::ny], ux[::nx,::ny], uy[::nx,::ny])
    # #vector sum - show in second figure
    ax.set_title('Vertical Displacement, Uz')
    ax.set_xlabel('EW Distance [km]')
    ax.set_ylabel('NS Distance [km]')
    # , pad=0.1) #ticks=MaxNLocator(nbins=5) #5 ticks only)
    cb = plt.colorbar(im, ax=ax, orientation='horizontal')
    cb.set_label('cm')

    # NOTE: look into code to see how error and wgt are determined..
    # sc1 = ax1.scatter(x_km,y_km,c=err,cmap=plt.cm.Reds)
    # im1 = ax1.imshow(ux)
    # im1 = ax1.pcolor(x,y,ux)
    im1 = ax1.imshow(ux, extent=extent, cmap=cmap)
    ax1.quiver(x[::nx, ::ny], y[::nx, ::ny], ux[::nx, ::ny], np.zeros_like(uy)[::nx, ::ny])
    ax1.set_title('EW Displacement, Ux')
    cb1 = plt.colorbar(im1, ax=ax1, orientation='horizontal')  # , pad=0.1)
    cb1.set_label('cm')

    # sc2 = ax2.scatter(x_km,y_km,c=wgt,cmap=plt.cm.Blues, norm=LogNorm())
    # cb2 = plt.colorbar(sc2, ax=ax2, orientation='horizontal', pad=0.1)
    # sc2 = ax2.scatter(x_km,y_km,c=wgt,cmap=plt.cm.Blues)
    # im2 = ax2.imshow(uy)
    # im2 = ax2.pcolor(x,y,uy)
    im2 = ax2.imshow(uy, extent=extent, cmap=cmap)
    ax2.quiver(x[::nx, ::ny], y[::nx, ::ny],
               np.zeros_like(ux)[::nx, ::ny], uy[::nx, ::ny])
    cb2 = plt.colorbar(im2, ax=ax2, orientation='horizontal')  # , pad=0.1)
    ax2.set_title('NS Displacement, Uy')
    cb2.set_label('cm')

    # Add crosshair grid:
    for a in (ax, ax1, ax2):
        a.axhline(linestyle=':', color='w')
        a.axvline(linestyle=':', color='w')

    plot_fault(fig, **params)

    plt.suptitle('Components of Deformation', fontsize=16, fontweight='bold')


def plot_los(x, y, ux, uy, uz, los, params, profile=False):
    ''' Separate figure showing displacement and Wrapped Phase in Radar     '''
    # Convert to km and cm for ploting
    cmap='bwr'
    x,y = np.array([x,y]) * 1e-3
    ux,uy,uz,los = np.array([ux,uy,uz,los]) * 1e2

    # extract a few varialbes from params dictionary
    inc = params['inc']
    ald = params['ald']
    wavelength = params['wavelength']

    # Set view
    nx=20
    ny=20
    extent = [x.min(), x.max(), y.min(), y.max()]

    # --------------
    fig, (ax,ax1,ax2) = plt.subplots(1,3,subplot_kw=dict(aspect=1.0, adjustable='box-forced'),sharex=True, sharey=True, figsize=(17,6))
    # vertical displacement w/ horizontal vectors # NOTE: add quiver scale arrow!
    im = ax.imshow(uz, extent=extent, cmap=cmap)
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
    im = ax2.imshow(los_wrapped, extent=extent, cmap='viridis', vmin=0, vmax=1)
    plot_los_indicator(ax2,ald)
    ax2.set_title('LOS Wrapped')
    cb = plt.colorbar(im, ax=ax2, orientation='horizontal', ticks=[], pad=0.1) #ticks=MaxNLocator(nbins=5) #5 ticks only)
    cb.set_label('{0} cm'.format(wavelength/2)) #1 color cycle = lambda/2

    plot_fault(fig,**params)

    # Add crosshair grid:
    for a in (ax,ax1,ax2):
        a.axhline(linestyle=':', color='w')
        a.axvline(linestyle=':', color='w')

    plt.suptitle('Model in Radar Line of Sight', fontsize=16, fontweight='bold')


def plot_profile(ind,ux,uy,uz,los,axis=0):
    ''' straight line profile through specified axis of ux,uy,uz,los'''
    # Show profile line in separate LOS plot
    plt.figure()
    plt.imshow(los, cmap='bwr')
    if axis == 0:
        plt.axhline(ind,color='k', lw=2)
        # plt.annotate('A',(ind,0))
    else:
        plt.vline(ind,color='k', lw=2)
    plt.title('Profile line')


    # Convert to km and cm for ploting
    # x,y = np.array([x,y]) * 1e-3
    ux,uy,uz,los = np.array([ux,uy,uz,los]) * 1e2

    if axis==1:
        ux,uy,uz,los = [x.T for x in [ux,uy,uz,los]]

    # extract profiles
    ux_p = ux[ind]
    uy_p = uy[ind]
    uz_p = uz[ind]
    los_p = los[ind]

    fig, (ax,ax1,ax2,ax3) = plt.subplots(1,4,sharex=True, sharey=True, figsize=(17,6))
    ax.plot(ux_p,'k.-',lw=2)
    ax.set_title('Ux')
    ax.set_ylabel('Displacement [cm]')
    ax.set_xlabel('Distance [km]')
    ax.text(0,0.9,'A',fontweight='bold',ma='center',transform=ax.transAxes)
    ax.text(0.9,0.9,'B',fontweight='bold',ma='center',transform=ax.transAxes)

    ax1.plot(uy_p,'k.-',lw=2)
    ax1.set_title('Uy')

    ax2.plot(uz_p,'k.-',lw=2)
    ax2.set_title('Uz')

    ax3.plot(los_p,'k.-',lw=2)
    ax3.set_title('LOS')

    for a in (ax,ax1,ax2,ax3):
        a.axhline(0,color='gray',linestyle='--')
        a.axvline(ind,color='gray',linestyle='--')
        a.grid(True)

    # ratio of uz to ur
    # NOTE: just print max ratio
    print('Ur/Uz = ', ux_p.max() / uz_p.max())
    '''
    plt.figure()
    plt.plot(uz_p/ux_p,'k.-',lw=2)
    plt.xlabel('Distance [km]')
    plt.ylabel('Uz/Ur Ratio')
    plt.title('Vertical vs. Horizontal Displacements')
    plt.show()
    '''