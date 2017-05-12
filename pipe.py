"""
Prolate spheroid in limit of Vertical pipe

Bonnacorso & Davis 1999
See also Lisowski Analytical Volcano Deformation Source Models
(Ch8 in Vocano Deformation, Dzurisin 2007)
"""
import numpy as np
from . import util

# =====================
# Forward Models
# =====================

def forward(x,y,xcen=0,ycen=0,d=5e3,l=5e3,r=1e3,dP=1e6,mu=1e9,nu=0.25):
    """
    Calculates surface deformation based on a closed pipe
    Reference: Lisowski & Dzuisin eq 8.21

    Args:
    ------------------
    x: x-coordinate grid (m)
    y: y-coordinate grid (m)

    Kwargs:
    -----------------
    xcen: y-offset of point source epicenter (m)
    ycen: y-offset of point source epicenter (m)
    d: depth to top of pipe (m)
    l: length of pip (m) (major axis diameter)
    r: source radius (m) (semiminor axis radius, a>>b=c)
    dP: change in Pressure (Pa)
    mu: shear modulus (Pa)
    nu: poisson's ratio for medium

    Returns:
    -------
    (ux, uy, uz)

    Examples:
    --------
    """
    # Center coordinate grid on point source
    x = x - xcen
    y = y - ycen
    # Convert to surface cylindrical coordinates
    th, rho = util.cart2pol(x,y)
    c1 = d #match notation of Lisowski
    c2 = d + l

    R1 = np.hypot(r,c1)
    R2 = np.hypot(r,c2)
    A = (dP*r**2)/(4*mu)
    T1 = (c1/R1)**3 + 2*c1*(-3+5*nu)/R1 + (5*c2**3*(1-2*nu) - (2*c2*rho**2*(-3 + 5*nu)))/R2**3
    ux = A*T1*x/rho**2
    uy = A*T1*y/rho**2

    T2 = (c1**2/R1**3) + 2*(-2+5*nu)/R1 +  (c2**2*(3-10*nu) - (2*rho**2*(-2 + 5*nu)))/R2**3
    uz = -A*T2

    return np.array([ux,uy,uz])
