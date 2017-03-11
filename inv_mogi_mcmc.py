# -*- coding: utf-8 -*-
"""
Mogi inversion for volcano-geodesy analytic models using Markov Chain Monte-Carlo

Author: Mario Angarita & Scott Henderson
Date: 12/12/2016
"""

import sys
import pyproj
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import util
import mogi
import scipy.optimize
from collections import OrderedDict
import pandas as pd
import pymc

#Reading the UNWRAPPED interferogram in plane coordinates
path = 'ruiz_desc_look.utm.phs.vrt'
data, extent, meta = util.load_rasterio(path)
los = data[1]/100

# Set up grid (pixel center coordinates)
r,c = np.indices(los.shape)
X, Y = util.world2rc(c,r,meta['affine'], inverse=True)

# Load Incidence
data, extent, meta = util.load_rasterio('los_ruiz_desc_look.utm.vrt')
incidence = data[0]

# Load Heading
heading=data[1]

#Run Varres before running this script if the image has a very high resolution
#Result of the quadtree made with Varres
f=open('newtest1.txt')
validos=f.readlines()
valids=np.zeros(los.shape,dtype=bool)
for i in validos:
    valids[int(i.split()[2]),int(i.split()[1])]=True
valid = np.isfinite(los) & np.isfinite(incidence) & np.isfinite(heading)
data = los[valid&valids]

#data_error is the standard deviation
data_error = np.ones_like(data)*0.05*data.max()
xargs = [X[valid&valids],Y[valid&valids],incidence[valid&valids],heading[valid&valids]]

def model(xargs,data): 
#Longitude and latitude guess for the dislocation center
    lat_cen = 0.79
    lon_cen = -77.94
    p = pyproj.Proj(proj='utm',zone=18,south=False,ellps='WGS84')
    xceni,yceni = p(lon_cen,lat_cen)
#Distribution for every parameter
    xcen=pymc.Normal('xcen', mu=xceni, tau=(1e2)**-2)
    ycen=pymc.Normal('ycen', mu=yceni, tau=(1e2)**-2)
    depth=pymc.Normal('depth', mu=10.2e3, tau=(1e3)**-2)
    dV=pymc.Normal('dV', mu=42.5e6, tau=(1e3)**-2)

#Deformation is the deterministic variable
    @pymc.deterministic(plot=False)
    def defo(xargs=xargs,xcen=xcen,ycen=ycen,depth=depth,dV=dV):
    	return mogi.invert(xargs=xargs,xcen=xcen,ycen=ycen,depth=depth,dV=dV)
#Probability distribution
    z = pymc.Normal('z', mu=defo, tau=1.0/data_error**2, value=data, observed=True)
    return locals()

#Making MCMC model
MDL = pymc.MCMC(model(xargs,data))

#Choosing Metropolis-Hastings for step method
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.xcen)
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.ycen)
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.depth)
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.dV)

#Number of runs
MDL.sample(20000, 10000, 10) 

#Writing the output in a csv file
MDL.write_csv("mcmc_mogi.csv", variables=["xcen","ycen","depth", "dV"])

#Listing the parameters output
means,error=np.loadtxt('mcmc_mogi.csv',delimiter=', ',usecols=(1,3),unpack=True,skiprows=1)
xoff1=float(means[0])
yoff1=float(means[1])
depth1=float(means[2])
dV1=float(means[3])

resultDict = OrderedDict([('xcen',xoff1),
            ('ycen',yoff1),
            ('d', depth1), #m
            ('dV', dV1)
            ])

# Run Mogi model with the inversion parameters as input
ux,uy,uz = mogi.forward(X,Y,**resultDict)

dataVec = np.dstack([ux, uy, uz])

# Get LOS transform
cart2los = util.get_cart2los(incidence,heading)
losfinal = np.sum(dataVec * cart2los, axis=2)

#Plot of Mogi model with the inversion parameters
plt.figure(1)
extent = [X.min(), X.max(), Y.min(), Y.max()]
plt.imshow(losfinal, extent=extent, cmap='jet',clim=(-np.nanmax(los),np.nanmax(los)))
plt.xlabel('Eastings(m)')
plt.ylabel('Northings(m)')
plt.title('Inversion')
cb=plt.colorbar()
cb.set_label('LOS Deformation(cm)')
plt.grid(True)
plt.savefig('inversion_ok_mcmc.png')

#Plot of residual Mogi model 
plt.figure(2)
extent = [X.min(), X.max(), Y.min(), Y.max()]
residual=los-losfinal
plt.imshow(residual, extent=extent, cmap='jet',clim=(-np.nanmax(los),np.nanmax(los)))
plt.xlabel('Eastings(m)')
plt.ylabel('Northings(m)')
plt.title('Residual')
cb=plt.colorbar()
cb.set_label('LOS Deformation(cm)')
plt.grid(True)
plt.savefig('residual_ok_mcmc.png')

#Plot of observed data
plt.figure(3)
extent = [X.min(), X.max(), Y.min(), Y.max()]
plt.imshow(los, extent=extent, cmap='jet',clim=(-np.nanmax(los),np.nanmax(los)))
plt.xlabel('Eastings(m)')
plt.ylabel('Northings(m)')
plt.title('Observed')
cb=plt.colorbar()
cb.set_label('LOS Deformation(cm)')
plt.grid(True)
plt.savefig('observed.png')

#Plotting the histograms
pymc.Matplot.plot(MDL)
