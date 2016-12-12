# -*- coding: utf-8 -*-
"""
Okada inversion for volcano-geodesy analytic models using least squares

Author: Mario Angarita & Scott Henderson
Date: 12/12/2016
"""

import sys
import pyproj
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import util
import okada
from collections import OrderedDict
import pandas as pd

#Reading the UNWRAPPED interferogram in plane coordinates
path = 'stack_rmg2.unw.utm.vrt'
data, extent, meta = util.load_rasterio(path)
los = data[0]
los[los==0] = np.nan

#Longitude and latitude guess for the dislocation center
lat_cen = 0.79
lon_cen = -77.94
p = pyproj.Proj(proj='utm',zone=18,south=False,ellps='WGS84')
xcen,ycen = p(lon_cen,lat_cen)

# Set up grid (pixel center coordinates)
r,c = np.indices(los.shape)
X, Y = util.world2rc(c,r,meta['affine'], inverse=True)

# Load Incidence
data, extent, meta = util.load_rasterio('incidencen.utm.tif')
incidence = data[0]

#Load Heading
data, extent, meta = util.load_rasterio('headingn.utm.tif')
heading = data[0]

#Initial parameters
params = OrderedDict([('xoff',xcen),
            ('yoff',ycen),
            ('depth', 2.4e3), #m
            ('length', 3.4e3), #m^3
            ('width', 5e3),
            ('opening', 0.0),
            ('slip', 0.6),
            ('dip', 50.0),
            ('strike', 213.0),
            ('rake', 151.0)
            ])

# Run Okada model with the initial parameters as input
ux,uy,uz = okada.forward(X,Y,**params)
dataVec = np.dstack([ux, uy, uz])

# Get LOS transform
cart2los = util.get_cart2los2(incidence,heading)
losguess = np.sum(dataVec * cart2los, axis=2)

#Initial guess for the least squares inversion
initial_guess = list(params.values())

#Run Varres before running this script if the image has a very high resolution
#Result of the quadtree made with Varres
f=open('newtest1.txt')
validos=f.readlines()
valids=np.zeros(los.shape,dtype=bool)
for i in validos:
    valids[int(i.split()[2]),int(i.split()[1])]=True
valid = np.isfinite(los) & np.isfinite(incidence) & np.isfinite(heading)
data=los[valids&valid]
xargs = [X[valids&valid],Y[valids&valid],incidence[valids&valid],heading[valids&valid]]
losguess1=losguess[valids&valid]

#curve_fit offers three options for least squares inversion (trf, dogleg and lm) for trf and dogleg is necessary to define some bounds for the parameters
#bounds=[[172000,80000,0.1e3,0,2e3,0.5e3,0,0,0],[174000,90000,5e3,60,10e3,7e3,10,360,360]]

#Bounded:
#popt,pcov = scipy.optimize.curve_fit(vmodels.okada.invert, xargs, data, p0=initial_guess, method='trf',bounds=bounds)

#Unbounded:
popt,pcov,info,mes,flag = scipy.optimize.curve_fit(vmodels.okada.invert, xargs, data, p0=initial_guess,full_output=True)

resultDict = OrderedDict([('xoff',popt[0]),
            ('yoff',popt[1]),
            ('depth', popt[2]), #m
            ('length', popt[4]), #m^3
            ('width', popt[5]),
            ('opening', 0),
            ('slip', popt[6]),
            ('dip', popt[3]),
            ('strike', popt[7]),
            ('rake', popt[8])
            ])

rmse=np.sqrt(np.sum( (data-info['fvec'])**2 ) / data.size)
print('Inversion Info:')
print(meta['crs'])
print('Number of Function Calls= {nfev}'.format(**info))
print('RMSE [m]= {:3e}'.format(rmse))

print('Best Fit Parameters:')
print('Easting= {:.4e} [m]'.format(popt[0]))
print('Northing= {:.4e} [m]'.format(popt[1]))
print('Depth= {:.4e} [m]'.format(popt[2]))
print('Opening= {:.4e} [m]'.format(0))
print('Dip= {:.4e} [degrees]'.format(popt[3]))
print('Length= {:.4e} [m]'.format(popt[4]))
print('Width= {:.4e} [m]'.format(popt[5]))
print('Slip= {:.4e} [m]'.format(popt[6]))
print('Strike= {:.4e} [degrees]'.format(popt[7]))
print('Rake= {:.4e} [degrees]'.format(popt[8]))
print('nu= 0.25')

perr = np.sqrt(np.diag(pcov))
print('Parameter 1 std Error:')
print('Easting= {:.4e} [m]'.format(perr[0]))
print('Northing= {:.4e} [m]'.format(perr[1]))
print('Depth= {:.4e} [m]'.format(perr[2]))
print('Opening= {:.4e} [m]'.format(0))
print('Dip= {:.4e} [degrees]'.format(perr[3]))
print('Length= {:.4e} [m]'.format(perr[4]))
print('Width= {:.4e} [m]'.format(perr[5]))
print('Slip= {:.4e} [m]'.format(perr[6]))
print('Strike= {:.4e} [degrees]'.format(perr[7]))
print('Rake= {:.4e} [degrees]'.format(perr[8]))
print('nu= 0.25')

#Listing the parameters 
resultDict = OrderedDict([('xoff',popt[0]),
            ('yoff',popt[1]),
            ('depth', popt[2]), #m
            ('length', popt[4]), #m^3
            ('width', popt[5]),
            ('opening', 0),
            ('slip', popt[6]),
            ('dip', popt[3]),
            ('strike', popt[7]),
            ('rake', popt[8])
            ])

# Run Okada model with the inversion parameters as input
ux,uy,uz = okada.forward(X,Y,**resultDict)
dataVec = np.dstack([ux, uy, uz]) 

# Get LOS transform
cart2los = util.get_cart2los(incidence,heading)
losfinal = np.sum(dataVec * cart2los, axis=2)

#Writing the results in a csv file
df = pd.DataFrame(resultDict, index=[0])
df.to_csv('mogi_inversion.csv', index=False)

#Plot of Okada model with the initial guess
plt.figure(1)
extent = [X.min(), X.max(), Y.min(), Y.max()]
plt.imshow(losguess*100, extent=extent, cmap='jet',clim=(-np.nanmax(los*100),np.nanmax(los*100)))
plt.xlabel('Eastings(m)')
plt.ylabel('Northings(m)')
plt.title('Model')
cb=plt.colorbar()
cb.set_label('LOS Deformation(cm)')
plt.grid(True)
plt.savefig('directo_ok.png')

#Plot of observed data
plt.figure(2)
extent = [X.min(), X.max(), Y.min(), Y.max()]
plt.imshow(los*100, extent=extent, cmap='jet',clim=(-np.nanmax(los*100),np.nanmax(los*100)))
plt.xlabel('Eastings(m)')
plt.ylabel('Northings(m)')
plt.title('Observed')
cb=plt.colorbar()
cb.set_label('LOS Deformation(cm)')
plt.grid(True)
plt.savefig('observed.png')

#Plot of Okada model with the inversion parameters
plt.figure(3)
extent = [X.min(), X.max(), Y.min(), Y.max()]
plt.imshow(losfinal*100, extent=extent, cmap='jet',clim=(-np.nanmax(los*100),np.nanmax(los*100)))
plt.xlabel('Eastings(m)')
plt.ylabel('Northings(m)')
plt.title('Inversion')
cb=plt.colorbar()
cb.set_label('LOS Deformation(cm)')
plt.grid(True)
plt.savefig('inversion_ok.png')

#Plot of residual Okada model 
plt.figure(4)
extent = [X.min(), X.max(), Y.min(), Y.max()]
residual=los-losfinal
plt.imshow(residual*100, extent=extent, cmap='jet',clim=(-np.nanmax(los*100),np.nanmax(los*100)))
plt.xlabel('Eastings(m)')
plt.ylabel('Northings(m)')
plt.title('Residual')
cb=plt.colorbar()
cb.set_label('LOS Deformation(cm)')
plt.grid(True)
plt.savefig('residual_ok.png')



