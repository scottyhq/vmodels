# -*- coding: utf-8 -*-
"""
Okada inversion for volcano-geodesy analytic models using Markov Chain Monte-Carlo

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
import okada
import scipy.optimize
from collections import OrderedDict
import pandas as pd
import pymc

#Reading the UNWRAPPED interferogram in plane coordinates
path = 'stack_rmg2.unw.utm.vrt'
data, extent, meta = util.load_rasterio(path)
los = data[0]

# Set up grid (pixel center coordinates)
r,c = np.indices(los.shape)
X, Y = util.world2rc(c,r,meta['affine'], inverse=True)

# Load Incidence
data, extent, meta = util.load_rasterio('incidencen.utm.tif')
incidence = data[0]

# Load Heading
data, extent, meta = util.load_rasterio('headingn.utm.tif')
heading = data[0]

#Run Varres before running this script if the image has a very high resolution
#Result of the quadtree made with Varres
f=open('newtest1.txt')
validos=f.readlines()
valids=np.zeros(los.shape,dtype=bool)
for i in validos:
    valids[int(i.split()[2]),int(i.split()[1])]=True
valid = np.isfinite(los) & np.isfinite(incidence) & np.isfinite(heading)
data=los[valids&valid]

#data_error is the standard deviation
data_error = np.ones_like(data)*0.05*data.max()
xargs = [X[valids&valid],Y[valids&valid],incidence[valids&valid],heading[valids&valid]]

def model(xargs,data): 

#Longitude and latitude guess for the dislocation center
    lat_cen = 0.79
    lon_cen = -77.94
    p = pyproj.Proj(proj='utm',zone=18,south=False,ellps='WGS84')
    xceni,yceni = p(lon_cen,lat_cen)

#Distribution for every parameter
    xoff=pymc.Normal('xoff', mu=173e3, tau=(2e3)**-2)
    yoff=pymc.Normal('yoff', mu=87.5e3, tau=(2.5e3)**-2)
    depth=pymc.Normal('depth', mu=1.5e3, tau=(0.2e3)**-2)
    dip=pymc.Normal('dip', mu=50, tau=(1)**-2)
    length=pymc.Normal('length', mu=3e3, tau=(1e3)**-2)
    width=pymc.Normal('width', mu=2e3, tau=(1e3)**-2)
    slip=pymc.Normal('slip', mu=0.6, tau=(0.1)**-2)
    strike=pymc.Normal('strike', mu=189, tau=(5)**-2)
    rake=pymc.Normal('rake', mu=123, tau=(5)**-2)

#Deformation is the deterministic variable
    @pymc.deterministic(plot=False)
    def defo(xargs=xargs,xoff=xoff,yoff=yoff,depth=depth,dip=dip,length=length,width=width,slip=slip,strike=strike,rake=rake):
    	return okada.invert(xargs=xargs,xoff=xoff,yoff=yoff,depth=depth,dip=dip,length=length,width=width,slip=slip,strike=strike,rake=rake)

#Probability distribution   	
    z = pymc.Normal('z', mu=defo, tau=1.0/data_error**2, value=data, observed=True)
    return locals()

#Making MCMC model
MDL = pymc.MCMC(model(xargs,data))

#Choosing Metropolis-Hastings for step method
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.xoff)
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.yoff)
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.depth)
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.dip)
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.length)
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.width)
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.slip)
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.strike)
MDL.use_step_method(pymc.AdaptiveMetropolis,MDL.rake)

#Number of runs
MDL.sample(200000, 100000, 10) 

#Writing the output in a csv file
MDL.write_csv("mcmc_okada.csv", variables=["xoff", "yoff", "depth", "dip", "width", "length", "slip", "strike", "rake"])

#Listing the parameters output
means,error=np.loadtxt('mcmc_okada.csv',delimiter=', ',usecols=(1,3),unpack=True,skiprows=1)
xoff1=float(means[0])
yoff1=float(means[1])
depth1=float(means[2])
dip1=float(means[3])
width1=float(means[4])
length1=float(means[5])
slip1=float(means[6])
strike1=float(means[7])
rake1=float(means[8])

#Listing the parameters errors
dxoff1=float(error[0])
dyoff1=float(error[1])
ddepth1=float(error[2])
ddip1=float(error[3])
dwidth1=float(error[4])
dlength1=float(error[5])
dslip1=float(error[6])
dstrike1=float(error[7])
drake1=float(error[8])


resultDict = OrderedDict([('xoff',xoff1),
            ('yoff',yoff1),
            ('depth', depth1), #m
            ('length', length1), #m^3
            ('width', width1),
            ('opening', 0),
            ('slip', slip1),
            ('dip', dip1),
            ('strike', strike1),
            ('rake', rake1)
            ])

# Run Okada model with the inversion parameters as input
ux,uy,uz = okada.forward(X,Y,**resultDict)
dataVec = np.dstack([ux, uy, uz])

# Get LOS transform
cart2los = util.get_cart2los(incidence,heading)
losfinal = np.sum(dataVec * cart2los, axis=2)

#Plot of Okada model with the inversion parameters
plt.figure(1)
extent = [X.min(), X.max(), Y.min(), Y.max()]
plt.imshow(losfinal*100, extent=extent, cmap='jet',clim=(-np.nanmax(los*100),np.nanmax(los*100)))
plt.xlabel('Eastings(m)')
plt.ylabel('Northings(m)')
plt.title('Inversion')
cb=plt.colorbar()
cb.set_label('LOS Deformation(cm)')
plt.grid(True)
plt.savefig('inversion_ok_mcmc.png')

#Plot of residual Okada model 
plt.figure(2)
extent = [X.min(), X.max(), Y.min(), Y.max()]
residual=los-losfinal
plt.imshow(residual*100, extent=extent, cmap='jet',clim=(-np.nanmax(los*100),np.nanmax(los*100)))
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
