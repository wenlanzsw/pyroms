#!/usr/bin/env python
# encoding: utf-8

"""
grid.py

Copyright (C) 2007, Robert Hetland

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the
   distribution.
3. The name of the author may not be used to endorse or promote
   products derived from this software without specific prior written
   permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.    
"""

from numpy import *
import pylab as pl
import netCDF4_classic as netcdf
import os
import warnings
import pyroms
from datetime import datetime, timedelta
from matplotlib.toolkits.basemap import Basemap
from matplotlib.toolkits.basemap.greatcircle import GreatCircle

def _nn_extrap(a):
    """extrapolate masked points in a 2D array to nearest neighbor"""
    if not isinstance(a, ma.MaskedArray): return a
    aout = a.data
    agood = a[~a.mask]
    ii, jj = indices(a.shape)
    igood = ii[~a.mask]
    jgood = jj[~a.mask]
    for i, j in zip(*where(a.mask)):
        d = (i-igood)**2 + (j-jgood)**2
        aout[i,j] = agood[d.flat==d.min()].mean()
    return aout


class Grid(object):
    
    def __init__(self, x=None, y=None,  \
                       lon=None, lat=None, \
                       pm=None, pn=None, angle=None, \
                       mask=None, proj=None, \
                       f=None, h=None, \
                       theta_s=None, theta_b=None, hc=None, N=None):
        """
        Need to add a bit of documentation here.
        
        Defining parts of the grid -- spatial points are verticies either in
        x/y lat/lon or both.  A projection object may be passed for posterity.
        pm, pn, and angle are all necessary.  For a working roms grid, f and h
        must also be added.
        
        inputs must be arrays or masked arrays.
        """
        self.x_vert = x
        self.y_vert = y
        
        self.lon_vert = lon
        self.lat_vert = lat
        
        self.geographic=False
        if  lon is not None and lat is not None:
            self.geographic=True
            if proj is None:
                self.proj = Basemap(projection='merc', resolution=None, lat_ts=0.0)
            else:
                self.proj = proj
        
        # Define relevant shapes
        if x is not None and y is not None:
            self.Mvert, self.Lvert = self.x_vert.shape
        elif lon is not None and lat is not None:
            self.Mvert, self.Lvert = self.lon_vert.shape
        else:
            raise Error, 'Must define either x/y or lat/lon grid verticies'
        
        self.Mp, self.Lp = self.Mvert-1, self.Lvert-1
        self.M, self.L = self.Mvert-2, self.Lvert-2
        self.Mm, self.Lm = self.Mvert-3, self.Lvert-3
        
        self.pm = pm
        self.pn = pn
        self.angle = angle
        
        if mask is not None:
            self.mask_rho = mask
        else:
            self.mask_rho = ones((self.Mp, self.Lp), dtype='d')
        
        if isinstance(self.x_vert, ma.MaskedArray):
            mask = (self.x_vert.mask[:-1,:-1] | self.x_vert.mask[1:,:-1] | \
                    self.x_vert.mask[:-1,1:] | self.x_vert.mask[1:,1:])
            self.mask_rho = asarray(~(~bool_(self.mask_rho) | mask), dtype='d')
        
        if isinstance(self.y_vert, ma.MaskedArray):
            mask = (self.y_vert.mask[:-1,:-1] | self.y_vert.mask[1:,:-1] | \
                    self.y_vert.mask[:-1,1:] | self.y_vert.mask[1:,1:])
            self.mask_rho = asarray(~(~bool_(self.mask_rho) | mask), dtype='d')
        
        if isinstance(self.lon_vert, ma.MaskedArray):
            mask = (self.lon_vert.mask[:-1,:-1] | self.lon_vert.mask[1:,:-1] | \
                    self.lon_vert.mask[:-1,1:] | self.lon_vert.mask[1:,1:])
            self.mask_rho = asarray(~(~bool_(self.mask_rho) | mask), dtype='d')
        
        if isinstance(self.lat_vert, ma.MaskedArray):
            mask = (self.lat_vert.mask[:-1,:-1] | self.lat_vert.mask[1:,:-1] | \
                    self.lat_vert.mask[:-1,1:] | self.lat_vert.mask[1:,1:])
            self.mask_rho = asarray(~(~bool_(self.mask_rho) | mask), dtype='d')
        
        if (self.x_vert is None or self.y_vert is None) and self.geographic:
            self.calc_projection()
        
        if self.pn is None or self.pm is None or self.angle is None:
            self.calc_metrics()
        
        if self.geographic and f is None:
            self.f = 2 * 7.29e-5 * cos(self.lat_rho * pi / 180.)
        else:
            self.f = asarray(f, dtype='d')
            
        self.h = h
        
        self.theta_s = theta_s
        self.theta_b = theta_b
        self.hc = hc
        self.N = N
    
    def _get_mask_u(self):
        return self.mask_rho[:,1:]*self.mask_rho[:,:-1]
    
    def _get_mask_v(self):
        return self.mask_rho[1:,:]*self.mask_rho[:-1,:]
    
    def _get_mask_psi(self):
        return self.mask_rho[1:,1:]*self.mask_rho[:-1,1:]* \
               self.mask_rho[1:,:-1]*self.mask_rho[:-1,:-1]
    
    def _get_x_rho(self):
        if self.x_vert is None or self.y_vert is None: return
        return 0.25*(self.x_vert[1:,1:]+self.x_vert[1:,:-1]+ \
                     self.x_vert[:-1,1:]+self.x_vert[:-1,:-1])
    
    def _get_y_rho(self):
        if self.x_vert is None or self.y_vert is None: return
        return 0.25*(self.y_vert[1:,1:]+self.y_vert[1:,:-1]+ \
                     self.y_vert[:-1,1:]+self.y_vert[:-1,:-1])
    
    def _get_x_u(self):
        if self.x_vert is None or self.y_vert is None: return
        return 0.5*(self.x_vert[:-1,1:-1] + self.x_vert[1:,1:-1])
    
    def _get_y_u(self):
        if self.x_vert is None or self.y_vert is None: return
        return 0.5*(self.y_vert[:-1,1:-1] + self.y_vert[1:,1:-1])
    
    def _get_x_v(self):
        if self.x_vert is None or self.y_vert is None: return
        return 0.5*(self.x_vert[1:-1,:-1] + self.x_vert[1:-1,1:])
    
    def _get_y_v(self):
        if self.x_vert is None or self.y_vert is None: return
        return 0.5*(self.y_vert[1:-1,:-1] + self.y_vert[1:-1,1:])
    
    def _get_x_psi(self):
        if self.x_vert is None or self.y_vert is None: return
        return self.x_vert[1:-1,1:-1]
    
    def _get_y_psi(self):
        if self.x_vert is None or self.y_vert is None: return
        return self.y_vert[1:-1,1:-1]
    
    def _get_lon_rho(self):
        if self.lon_vert is None or self.lat_vert is None: return
        return 0.25*(self.lon_vert[1:,1:]+self.lon_vert[1:,:-1]+ \
                     self.lon_vert[:-1,1:]+self.lon_vert[:-1,:-1])
    
    def _get_lat_rho(self):
        if self.lon_vert is None or self.lat_vert is None: return
        return 0.25*(self.lat_vert[1:,1:]+self.lat_vert[1:,:-1]+ \
                     self.lat_vert[:-1,1:]+self.lat_vert[:-1,:-1])
    
    def _get_lon_u(self):
        if self.lon_vert is None or self.lat_vert is None: return
        return 0.5*(self.lon_vert[:-1,1:-1] + self.lon_vert[1:,1:-1])
    
    def _get_lat_u(self):
        if self.lon_vert is None or self.lat_vert is None: return
        return 0.5*(self.lat_vert[:-1,1:-1] + self.lat_vert[1:,1:-1])
    
    def _get_lon_v(self):
        if self.lon_vert is None or self.lat_vert is None: return
        return 0.5*(self.lon_vert[1:-1,:-1] + self.lon_vert[1:-1,1:])
    
    def _get_lat_v(self):
        if self.lon_vert is None or self.lat_vert is None: return
        return 0.5*(self.lat_vert[1:-1,:-1] + self.lat_vert[1:-1,1:])
    
    def _get_lon_psi(self):
        if self.lon_vert is None or self.lat_vert is None: return
        return self.lon_vert[1:-1,1:-1]
    
    def _get_lat_psi(self):
        if self.lon_vert is None or self.lat_vert is None: return
        return self.lat_vert[1:-1,1:-1]
    
    
    def _get_sc_w(self):
        if None in (self.theta_s, self.theta_b, self.hc, self.N): return
        return mgrid[-1.0:0.0:1j*(self.N+1)]
    
    def _get_sc_r(self):
        if None in (self.theta_s, self.theta_b, self.hc, self.N): return
        sc_w = mgrid[-1.0:0.0:1j*(self.N+1)]
        return 0.5*(sc_w[1:]+sc_w[:-1])
    
    def _get_Cs_r(self):
        if None in (self.theta_s, self.theta_b, self.hc, self.N): return
        if self.theta_s == 0.0: return self._get_sc_r()
        return (1-self.theta_b)*sinh(self.theta_s*self._get_sc_r())/ \
               sinh(self.theta_s)+0.5*self.theta_b \
               *(tanh(self.theta_s*(self._get_sc_r()+0.5)) \
               - tanh(0.5*self.theta_s))/tanh(0.5*self.theta_s)
    
    def _get_Cs_w(self):
        if None in (self.theta_s, self.theta_b, self.hc, self.N): return
        if self.theta_s == 0.0: return self._get_sc_w()
        return (1-self.theta_b)*sinh(self.theta_s*self._get_sc_w())/ \
               sinh(self.theta_s)+0.5*self.theta_b \
               *(tanh(self.theta_s*(self._get_sc_w()+0.5)) \
                - tanh(0.5*self.theta_s))/tanh(0.5*self.theta_s)
    
    
    def calc_metrics(self):
        'Calculates pm, pn, dndx, dmde, and angle from x_vert and y_vert'
        if self.geographic:
            gc_dist = vectorize(lambda lon1, lat1, lon2, lat2: \
                              GreatCircle(6378137.0, 6356752.3142, \
                                          lon1, lat1, lon2, lat2).distance)
            lon_temp = 0.5*(self.lon_vert[1:,:]+self.lon_vert[:-1,:])
            lat_temp = 0.5*(self.lat_vert[1:,:]+self.lat_vert[:-1,:])
            if isinstance(lat_temp, ma.MaskedArray): lat_temp = lat_temp.filled(0.0)
            if isinstance(lon_temp, ma.MaskedArray): lon_temp = lon_temp.filled(0.0)
            self.pm = 1.0 / gc_dist(lon_temp[:,1:],  lat_temp[:,1:], \
                                    lon_temp[:,:-1], lat_temp[:,:-1])
            lon_temp = 0.5*(self.lon_vert[:,1:]+self.lon_vert[:,:-1])
            lat_temp = 0.5*(self.lat_vert[:,1:]+self.lat_vert[:,:-1])
            if isinstance(lat_temp, ma.MaskedArray): lat_temp = lat_temp.filled(0.0)
            if isinstance(lon_temp, ma.MaskedArray): lon_temp = lon_temp.filled(0.0)
            self.pn = 1.0 / gc_dist(lon_temp[1:,:],  lat_temp[1:,:], \
                                    lon_temp[:-1,:], lat_temp[:-1,:])
        else:
            x_temp = 0.5*(self.x_vert[1:,:]+self.x_vert[:-1,:])
            y_temp = 0.5*(self.y_vert[1:,:]+self.y_vert[:-1,:])
            self.pm = 1.0 / sqrt(diff(x_temp, axis=1)**2 + diff(y_temp, axis=1)**2)
            x_temp = 0.5*(self.x_vert[:,1:]+self.x_vert[:,:-1])
            y_temp = 0.5*(self.y_vert[:,1:]+self.y_vert[:,:-1])
            self.pn = 1.0 / sqrt(diff(x_temp, axis=0)**2 + diff(y_temp, axis=0)**2)
        
        if any(~isfinite(self.pm)) or any(~isfinite(self.pm)):
             self.pm = ma.masked_where(~isfinite(self.pm), self.pm)
             self.pn = ma.masked_where(~isfinite(self.pn), self.pn)
        
        if isinstance(self.pn, ma.MaskedArray):
            self.dndx = ma.zeros((self.Mp, self.Lp), dtype='d')
        else:
            self.dndx = zeros((self.Mp, self.Lp), dtype='d')
        
        if isinstance(self.pm, ma.MaskedArray):
            self.dmde = ma.zeros((self.Mp, self.Lp), dtype='d')
        else:
            self.dmde = zeros((self.Mp, self.Lp), dtype='d')
        
        self.dndx[1:-1,1:-1] = 0.5*(1.0/self.pn[1:-1,2:] - 1.0/self.pn[1:-1,:-2])
        self.dmde[1:-1,1:-1] = 0.5*(1.0/self.pm[2:,1:-1] - 1.0/self.pm[:-2,1:-1])
        
        if self.x_vert is None or self.y_vert is None:
            self.calc_projection()
        
        self.angle = arctan2(diff(0.5*(self.y_vert[1:,:]+self.y_vert[:-1,:])), \
                           diff(0.5*(self.x_vert[1:,:]+self.x_vert[:-1,:])))
    
    def calc_projection(self, proj=None):        
        if isinstance(self.lat_vert, ma.MaskedArray):
            mask_lat = self.lat_vert.mask 
            lat_temp = self.lat_vert.filled(0.0)
        else:
            lat_temp = self.lat_vert
        
        if isinstance(self.lon_vert, ma.MaskedArray): 
            mask_lon = self.lon_vert.mask
            lon_temp = self.lon_vert.filled(0.0)
        else:
            lon_temp = self.lon_vert
        
        self.x_vert, self.y_vert = self.proj(lon_temp, lat_temp)
        
        if isinstance(self.lon_vert, ma.MaskedArray):
            self.x_vert = ma.masked_array(self.x_vert, mask=mask_lon)
        
        if isinstance(self.lat_vert, ma.MaskedArray):
            self.y_vert = ma.masked_array(self.y_vert, mask=mask_lat)
    
    def calc_orthogonality(self):
        '''
        Calculate orthogonality error in radiens
        '''
        z = self.x_vert + 1j*self.y_vert
        du = diff(z, axis=1); du = (du/abs(du))[:-1,:]
        dv = diff(z, axis=0); dv = (dv/abs(dv))[:,:-1]
        ang1 = arccos(du.real*dv.real + du.imag*dv.imag)
        du = diff(z, axis=1); du = (du/abs(du))[1:,:]
        dv = diff(z, axis=0); dv = (dv/abs(dv))[:,:-1]
        ang2 = arccos(du.real*dv.real + du.imag*dv.imag)
        du = diff(z, axis=1); du = (du/abs(du))[:-1,:]
        dv = diff(z, axis=0); dv = (dv/abs(dv))[:,1:]
        ang3 = arccos(du.real*dv.real + du.imag*dv.imag)
        du = diff(z, axis=1); du = (du/abs(du))[1:,:]
        dv = diff(z, axis=0); dv = (dv/abs(dv))[:,1:]
        ang4 = arccos(du.real*dv.real + du.imag*dv.imag)
        ang = mean([abs(ang1), abs(ang2), abs(ang3), abs(ang4)], axis=0)
        ang = (ang-pi/2.0)
        return ang
    
    def maskpoly(self, polyverts, inverse=False, geographic=None):
        """
        Mask Cartesian points contained within the polygons contained in the list 'polygons'.
        
        A cell is masked if the cell center (x_rho, y_rho) is within the polygon.
        Other sub-masks (mask_u, mask_v, and mask_psi) are updated automatically.
        """
        if geographic is None:
            geographic = self.geographic
        
        mask = self.mask_rho
        
        if inverse:
            mask = asarray(~bool_(mask), dtype='d')
        
        iwater = mask == 1.0
        if geographic:
            x_wet = self._get_lon_rho()[iwater]
            y_wet = self._get_lat_rho()[iwater]
        else:
            x_wet = self._get_x_rho()[iwater]
            y_wet = self._get_y_rho()[iwater]
        
        mask_wet = mask[iwater]
        
        inside = pyroms.Polygon(polyverts).inside(zip(x_wet, y_wet))
        
        if any(inside):
            mask_wet[inside] = 0.0
            mask[iwater] = mask_wet
            if inverse:
                mask = asarray(~bool_(a), dtype='d')
            self.mask_rho = mask
    
    def write_roms_grid(self, filename='ocean_grd.nc', full_output=True):
        
        Mp, Lp = self.x_rho.shape
        M, L = self.x_psi.shape
        
        xl = self.x_rho[self.mask_rho==1.0].ptp()
        el = self.y_rho[self.mask_rho==1.0].ptp()
        
        # Write ROMS grid to file
        nc = netcdf.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
        nc.Description = 'ROMS grid'
        nc.Author = 'pyroms.gridgen'
        nc.Created = datetime.now().isoformat()
        nc.type = 'ROMS GRD file'
        
        nc.createDimension('xi_rho', Lp)
        nc.createDimension('xi_u', L)
        nc.createDimension('xi_v', Lp)
        nc.createDimension('xi_psi', L)
        
        nc.createDimension('eta_rho', Mp)
        nc.createDimension('eta_u', Mp)
        nc.createDimension('eta_v', M)
        nc.createDimension('eta_psi', M)
        
        nc.createDimension('xi_vert', Lp+1)
        nc.createDimension('eta_vert', Mp+1)
        
        nc.createVariable('xl', 'f8', ())
        nc.variables['xl'].units = 'meters'
        nc.variables['xl'] = xl
        
        nc.createVariable('el', 'f8', ())
        nc.variables['el'].units = 'meters'
        nc.variables['el'] = el
        
        nc.createVariable('spherical', 'S1', ())
        nc.variables['spherical'] = 'F'
        
        def write_nc_var(var, name, dimensions, units=None):
            nc.createVariable(name, 'f8', dimensions)
            if units is not None:
                nc.variables[name].units = units
            nc.variables[name][:] = _nn_extrap(var)
        
        write_nc_var(self.pm, 'pm', ('eta_rho', 'xi_rho'), 'meters-1')
        write_nc_var(self.pn, 'pn', ('eta_rho', 'xi_rho'), 'meters-1')
        write_nc_var(self.dmde, 'dmde', ('eta_rho', 'xi_rho'))
        write_nc_var(self.dndx, 'dndx', ('eta_rho', 'xi_rho'))
        write_nc_var(self.angle, 'angle', ('eta_rho', 'xi_rho'))
        
        write_nc_var(self.f, 'f', ('eta_rho', 'xi_rho'), 'seconds-1')
        write_nc_var(self.h, 'h', ('eta_rho', 'xi_rho'), 'meters')
        
        write_nc_var(self.mask_rho, 'mask_rho', ('eta_rho', 'xi_rho'))
        write_nc_var(self.mask_u, 'mask_u', ('eta_u', 'xi_u'))
        write_nc_var(self.mask_v, 'mask_v', ('eta_v', 'xi_v'))
        write_nc_var(self.mask_psi, 'mask_psi', ('eta_psi', 'xi_psi'))
        
        if full_output:
            write_nc_var(self.x_vert, 'x_vert', ('eta_vert', 'xi_vert'), 'meters')
            write_nc_var(self.y_vert, 'y_vert', ('eta_vert', 'xi_vert'), 'meters')
            write_nc_var(self.x_rho, 'x_rho', ('eta_rho', 'xi_rho'), 'meters')
            write_nc_var(self.y_rho, 'y_rho', ('eta_rho', 'xi_rho'), 'meters')
            write_nc_var(self.x_u, 'x_u', ('eta_u', 'xi_u'), 'meters')
            write_nc_var(self.y_u, 'y_u', ('eta_u', 'xi_u'), 'meters')
            write_nc_var(self.x_v, 'x_v', ('eta_v', 'xi_v'), 'meters')
            write_nc_var(self.y_v, 'y_v', ('eta_v', 'xi_v'), 'meters')
            write_nc_var(self.x_psi, 'x_psi', ('eta_psi', 'xi_psi'), 'meters')
            write_nc_var(self.y_psi, 'y_psi', ('eta_psi', 'xi_psi'), 'meters')
            if self.geographic:
                write_nc_var(self.lon_vert, 'lon_vert', ('eta_vert', 'xi_vert'), 'meters')
                write_nc_var(self.lat_vert, 'lat_vert', ('eta_vert', 'xi_vert'), 'meters')
                write_nc_var(self.lon_rho, 'lon_rho', ('eta_rho', 'xi_rho'), 'meters')
                write_nc_var(self.lat_rho, 'lat_rho', ('eta_rho', 'xi_rho'), 'meters')
                write_nc_var(self.lon_u, 'lon_u', ('eta_u', 'xi_u'), 'meters')
                write_nc_var(self.lat_u, 'lat_u', ('eta_u', 'xi_u'), 'meters')
                write_nc_var(self.lon_v, 'lon_v', ('eta_v', 'xi_v'), 'meters')
                write_nc_var(self.lat_v, 'lat_v', ('eta_v', 'xi_v'), 'meters')
                write_nc_var(self.lon_psi, 'lon_psi', ('eta_psi', 'xi_psi'), 'meters')
                write_nc_var(self.lat_psi, 'lat_psi', ('eta_psi', 'xi_psi'), 'meters')
        
        nc.close()
    
    
    def calc_z_w(self, zeta=None):
        if None in (self.theta_s, self.theta_b, self.hc, self.N): return
        if self.h is not None and any(self.h<self.hc):
            raise ValueError, 'hc is less than minimum depth %f' % self.h.min
        sc_w = self._get_sc_w()[:, newaxis, newaxis]
        Cs_w = self._get_Cs_w()[:, newaxis, newaxis]
        z_w = (sc_w-Cs_w)*self.hc + Cs_w*self.h
        if zeta is not None:
            z_w = zeta*(1.0 + z_w/self.h)
        return z_w
    
    def calc_z_r(self, zeta=None):
        if None in (self.theta_s, self.theta_b, self.hc, self.N): return
        if self.h is not None and any(self.h<self.hc):
            raise ValueError, 'hc is less than minimum depth %f' % h.min
        sc_r = self._get_sc_r()[:, newaxis, newaxis]
        Cs_r = self._get_Cs_r()[:, newaxis, newaxis]
        z_r = (sc_r-Cs_r)*self.hc + Cs_r*self.h
        if zeta is not None:
            z_r = zeta*(1.0 + z_r/self.h)
        return z_r
    
    
    x = property(lambda self: self.x_vert)
    y = property(lambda self: self.y_vert)
    lon = property(lambda self: self.lon_vert)
    lat = property(lambda self: self.lat_vert)
    mask = property(lambda self: self.mask_rho)
    
    mask_u   = property(_get_mask_u)
    mask_v   = property(_get_mask_v)
    mask_psi = property(_get_mask_psi)
    x_rho = property(_get_x_rho)
    x_u   = property(_get_x_u)
    x_v   = property(_get_x_v)
    x_psi = property(_get_x_psi)
    y_rho = property(_get_y_rho)
    y_u   = property(_get_y_u)
    y_v   = property(_get_y_v)
    y_psi = property(_get_y_psi)
    lon_rho = property(_get_lon_rho)
    lon_u   = property(_get_lon_u)
    lon_v   = property(_get_lon_v)
    lon_psi = property(_get_lon_psi)
    lat_rho = property(_get_lat_rho)
    lat_u   = property(_get_lat_u)
    lat_v   = property(_get_lat_v)
    lat_psi = property(_get_lat_psi)
    
    sc_w = property(_get_sc_w)
    sc_r = property(_get_sc_r)
    Cs_w = property(_get_Cs_w)
    Cs_r = property(_get_Cs_r)
    z_r = property(calc_z_r)
    z_w = property(calc_z_w)


def gridgen(xbry, ybry, beta, shape, focus=None, ul_idx=0, \
            geographic=False, proj=None, \
            nnodes=14, precision=1.0e-12, newton=1, thin=1, \
            checksimplepoly=1, nppe=3, \
            verbose=False, windows=False):
    """
    create grid object from gridgen code by Pavel Sakov.  See:
        http://www.marine.csiro.au/~sakov/
    
    grid = Gridgen(xbry, ybry, beta, shape, focus=None, star=0,
                 nnodes=14, precision=1.0e-12, newton=1, thin=1,
                 checksimplepoly=1, nppe=3, verbose=False)
    
    input:
        xbry, ybry, beta, ul_idx=0:
            the input polygon represents the boundaries of the domain
            to be calculated.  The polygon should be defined counter-
            clockwise (positive).  Beta represent the locations where the
            boundary makes a 90 turn to the left (beta=1) or the right
            (beta=-1).  The index 'ul_idx' is the index of the upper left
            corner of the domain.  There may be more than one posibility,
            but the output grid will not change.  See the gridgen documentation
            for more details (this is the 'starred' entry in the bry file)
        geographic:
            If geographic is True (default is False), the coordinates are
            projected using the Basemap instance proj.  The grid instance
            that is created is a geographic grid, with both lon/lat and
            x/y defined.
        proj:
            Is the Basemap instance used to project the geographic coordinates
            to cartesian coordinates before calling gridgen.  The default is
                    proj = Basemap(projection='merc', lat_ts=0.0)
        shape:
            the number of points in the grid (ny, nx).  When creating a
            ROMS grid, note that Lm = nx-3 and Mm = ny-3.
        focus:
            The focus function returns values on a grid between 0 and 1,
            given inputs as a uniform unit grid.  E.g.,
                    x, y = mgrid[0:1:100j, 0:1:100j]
                    xf, yf = focus(x, y)
            The resulting grid defines where the output will have increased
            resolution.  If focus=None, the output grid will be uniformly
            distributed.
        other keyword arguments:
            follow the parameter definitions in the gridgen prm file.
    
    returns:
        ROMS Grid instance.
    """
    xbry = asarray(xbry)
    ybry = asarray(ybry)
    beta = asarray(beta)
    
    if proj is None:
        proj = Basemap(projection='merc', resolution=None, lat_ts=0.0)
        
    if geographic:
        xbry, ybry = proj(xbry, ybry)
    
    assert beta.sum() == 4.0, 'sum of beta must be 4.0'
    star = ['']*len(xbry)
    star[ul_idx] = '*'
    ny = shape[0]
    nx = shape[1]
    
    # tempnam gives security warning
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    prmname = os.tempnam()    # prm
    bryname = os.tempnam()    # input
    focusname = os.tempnam()  # grid
    rectname = os.tempnam()   # rectangle
    sigmaname = os.tempnam()  # sigmas
    gridname = os.tempnam()   # output
    
    # Write BRY file (input polygon)
    f = open(bryname, 'w')
    f.writelines(["%f %f %d%s\n" % (x, y, b, s) \
                  for x, y, b, s \
                  in zip(xbry, ybry, beta, star)])
    f.close()
    
    # Write FOCUS file
    if focus is not None:
        y, x = mgrid[0:1:ny*1j, 0:1:nx*1j]
        xfocus, yfocus = focus(x, y)
        f = open(focusname, 'w')
        f.writelines(["%f %f\n" % (x, y) \
                      for x, y \
                      in zip(xfocus.flat, yfocus.flat)])
        f.close()
    
    # Write PRM file
    f = open(prmname, 'w')
    f.write('input %s\n' % bryname)
    f.write('output %s\n' % gridname)
    if focus is not None:
        f.write('grid %s\n' % focusname)
    f.write('nx %d\n' % nx)
    f.write('ny %d\n' % ny)
    f.write('precision %g\n' % precision)
    f.write('thin %d\n' % thin)
    f.write('checksimplepoly %d\n' % checksimplepoly)
    f.write('newton %d\n' % newton)
    f.write('sigmas %s\n' % sigmaname)
    f.write('rectangle %s\n' % rectname)
    f.write('nppe %d\n' % nppe)
    f.close()
    
    # Run gridgen
    if verbose:
        verbstr = '-v'
    else:
        verbstr = ''
    os.system('gridgen %s %s' % (verbstr, prmname))
    
    # Read in grid file
    if windows == False:
        xp, yp = pl.load(gridname, comments="#", unpack=True)
    else:
        x=[]; y=[]
        gridfile = open(gridname)
        for line in gridfile:
            if '#' in line:
                x.append(nan)
                y.append(nan)
            else:
                data = line.split()
                x.append(float(data[0]))
                y.append(float(data[1]))
        xp = asarray(x)
        yp = asarray(y)
        gridfile.close()
    xp = xp.reshape(ny, nx)
    yp = yp.reshape(ny, nx)
    # remove temporary files
    
    try:
        [os.remove(file) \
         for file \
         in (prmname, bryname, focusname, rectname, sigmaname, gridname)]
    except:
        pass
    
    if any(isnan(xp)) or any(isnan(yp)):
        xp = ma.masked_where(isnan(xp), xp)
        yp = ma.masked_where(isnan(yp), yp)
    
    if geographic:
        if isinstance(xp, ma.MaskedArray): xp=xp.filled(nan)
        if isinstance(yp, ma.MaskedArray): yp=yp.filled(nan)
        lon, lat = proj(xp, yp, inverse=True)
        lon = ma.masked_where(isnan(lon), lon)
        lat = ma.masked_where(isnan(lat), lat)
        return Grid(lon=lon, lat=lat, proj=proj)
    else:
        return Grid(x=xp, y=yp)
    

def rho_to_vert(xr, yr, pm, pn, ang):
    Mp, Lp = xr.shape
    x = empty((Mp+1, Lp+1), dtype='d')
    y = empty((Mp+1, Lp+1), dtype='d')
    x[1:-1, 1:-1] = 0.25*(xr[1:,1:]+xr[1:,:-1]+xr[:-1,1:]+xr[:-1,:-1])
    y[1:-1, 1:-1] = 0.25*(yr[1:,1:]+yr[1:,:-1]+yr[:-1,1:]+yr[:-1,:-1])
    
    # east side
    theta = 0.5*(ang[:-1,-1]+ang[1:,-1])
    dx = 0.5*(1.0/pm[:-1,-1]+1.0/pm[1:,-1])
    dy = 0.5*(1.0/pn[:-1,-1]+1.0/pn[1:,-1])
    x[1:-1,-1] = x[1:-1,-2] + dx*cos(theta)
    y[1:-1,-1] = y[1:-1,-2] + dx*sin(theta)
    
    # west side
    theta = 0.5*(ang[:-1,0]+ang[1:,0])
    dx = 0.5*(1.0/pm[:-1,0]+1.0/pm[1:,0])
    dy = 0.5*(1.0/pn[:-1,0]+1.0/pn[1:,0])
    x[1:-1,0] = x[1:-1,1] - dx*cos(theta)
    y[1:-1,0] = y[1:-1,1] - dx*sin(theta)
    
    # north side
    theta = 0.5*(ang[-1,:-1]+ang[-1,1:])
    dx = 0.5*(1.0/pm[-1,:-1]+1.0/pm[-1,1:])
    dy = 0.5*(1.0/pn[-1,:-1]+1.0/pn[-1,1:])
    x[-1,1:-1] = x[-2,1:-1] - dy*sin(theta)
    y[-1,1:-1] = y[-2,1:-1] + dy*cos(theta)
    
    # here we are now going to the south side..
    theta = 0.5*(ang[0,:-1]+ang[0,1:])
    dx = 0.5*(1.0/pm[0,:-1]+1.0/pm[0,1:])
    dy = 0.5*(1.0/pn[0,:-1]+1.0/pn[0,1:])
    x[0,1:-1] = x[1,1:-1] + dy*sin(theta)
    y[0,1:-1] = y[1,1:-1] - dy*cos(theta)
    
    #Corners
    x[0,0] = 4.0*xr[0,0]-x[1,0]-x[0,1]-x[1,1]
    x[-1,0] = 4.0*xr[-1,0]-x[-2,0]-x[-1,1]-x[-2,1]
    x[0,-1] = 4.0*xr[0,-1]-x[0,-2]-x[1,-1]-x[1,-2]
    x[-1,-1] = 4.0*xr[-1,-1]-x[-2,-2]-x[-2,-1]-x[-1,-2]
    
    y[0,0] = 4.0*yr[0,0]-y[1,0]-y[0,1]-y[1,1]
    y[-1,0] = 4.0*yr[-1,0]-y[-2,0]-y[-1,1]-y[-2,1]
    y[0,-1] = 4.0*yr[0,-1]-y[0,-2]-y[1,-1]-y[1,-2]
    y[-1,-1] = 4.0*yr[-1,-1]-y[-2,-2]-y[-2,-1]-y[-1,-2]
    
    return x, y

def nc_grid(nc):
    '''
    Return grid class based on ROMS GRD file, or other output file with grid
    information.  The NetCDF file must contain either 'vert' veriables, or the
    verticies will be calculated with 'rho' and angle points.
    '''
    nc = pyroms.Dataset(nc)
    
    varlist = ['h', 'f', 'pm', 'pn', 'angle', 'theta_s', 'theta_b', 'hc']
    variables={}
    for var in varlist:
        try: variables[var] = nc.variables[var][:]
        except: variables[var] = None
    
    try: variables['mask']=nc.variables['mask_rho'][:]
    except: variables['mask']=None
    
    try: variables['N']=len(nc.dimensions['N'])
    except: variables['N']=None
    
    if 'x_vert' in nc.variables.keys() and 'y_vert' in nc.variables.keys():
        x = nc.variables['x_vert'][:]
        y = nc.variables['y_vert'][:]
        
        if any(isnan(x)): x = ma.masked_where(isnan(x), x)
        if any(isnan(y)): y = ma.masked_where(isnan(y), y)
        
        try: lon=nc.variables['lon_vert'][:]
        except: lon=None
        
        try: lat=nc.variables['lat_vert'][:]
        except: lat=None
        
        return Grid(x=x, y=y, lon=lon, lat=lat, **variables)
        
    elif 'lon_vert' in nc.variables.keys() and 'lat_vert' in nc.variables.keys():
        lon=nc.variables['lon_vert'][:]
        lat=nc.variables['lat_vert'][:]
        return Grid(lon=lon, lat=lat, **variables)
    
    else:
        try:
            xr = nc.variables['x_rho'][:]
            yr = nc.variables['y_rho'][:]
            pm = nc.variables['pm'][:]
            pn = nc.variables['pn'][:]
        except:
            raise ValueError, 'NetCDF file must contain x_rho, y_rho, pm, and pn'
        
        ang = variables['angle']
        if ang is None:
            ang = zeros(xr.shape, dtype='d')
        
        x, y = rho_to_vert(xr, yr, pm, pn, ang)
        
        if 'x_psi' in nc.variables.keys() and 'y_psi' in nc.variables.keys():
            xp = nc.variables['x_psi'][:]
            yp = nc.variables['y_psi'][:]
            x[1:-1, 1:-1] = xp
            y[1:-1, 1:-1] = yp
        
        return Grid(x=x, y=y, **variables)


def test_gridgen():
    xbry = array([10.,  5.,  5., 0.5,  0.,  5.,  5., 10.])
    ybry = array([10., 10.,  6., 7.5,  6.,  4.,  0.,  0.])
    beta = array([ 1.,  1., -1., 1. ,  1., -1.,  1.,  1.])
    
    def focus(x, y, xo=0.55, yo=0.45):
        xf = tan((x - xo)*2.0)
        yf = tan((y - yo)*2.0)
        xf -= xf.min()
        xf /= xf.max()
        yf -= yf.min()
        yf /= yf.max()
        return xf, yf
    
    # Run Gridgen
    grid = gridgen(xbry, ybry, beta, (36, 36), focus=focus, ul_idx=1)
    
    # Islands for masking
    xmask = array([8., 7., 6., 7.])
    ymask = array([4., 5., 4., 3.])
    grid.maskpoly(zip(xmask, ymask))
    
    return grid

def test_make_grid():
    yv, xv = mgrid[0:20, 0:20:0.5]
    grd = Grid(x=xv, y=yv)
    
    def rot2d(x, y, ang):
        'rotate vectors by geometric angle'
        xr = x*cos(ang) - y*sin(ang)
        yr = x*sin(ang) + y*cos(ang)
        return xr, yr
    
    xvr, yvr = rot2d(xv, yv, pi/4)
    
    grd_r = Grid(x=xvr, y=yvr)
    
    print 'pn match? ', allclose(grd.pn, grd_r.pn)
    print 'pm match? ', allclose(grd.pn, grd_r.pn)
    
    dx = 1.0/grd_r.pm
    print 'dx min, max, and mean = ', dx.min(), dx.max(), dx.mean()
    dy = 1.0/grd_r.pn
    print 'dy min, max, and mean = ', dy.min(), dy.max(), dy.mean()
    
    latv, lonv = mgrid[15.0:32.0, -100.0:-80.0]
    lonv[0:5, 0:5] = nan
    lonv=ma.masked_where(isnan(lonv), lonv)
    grd_geo = Grid(lon=lonv, lat=latv)
    dx = 1.0/grd_geo.pm
    print 'dx min, max, and mean = ', dx.min(), dx.max(), dx.mean()
    dy = 1.0/grd_geo.pn
    print 'dy min, max, and mean = ', dy.min(), dy.max(), dy.mean()
    print grd_geo.f

def test_masking():
    xi, yi = 20*random.rand(2, 100)
    hi = exp( -(xi-10.)**2/10.0 -(yi-10.)**2/5.0)
    h = grid.extrap_xy_to_grid(xi, yi, hi)
    pl.contour(grid.x_rho, grid.y_rho, h)
    print h
    pl.show()

def test_make_masked_grid():
    yv, xv = mgrid[0:20, 0:20:0.5]
    mask = ones(xv.shape, 'd')
    mask[0:10, 0:20] = 0.0
    xv = ma.masked_where(mask==0, xv)
    yv = ma.masked_where(mask==0, yv)
    return make_cart_grid(xv, yv)

def test_grid_3d():
    yv, xv = mgrid[0:10, 0:10:0.5]
    grd = Grid(x=xv, y=yv)
    grd.h = 11.0 + 19*random.rand(*grd.x_rho.shape)
    grd.theta_s = 5.0
    grd.theta_b = 1.0
    grd.hc = 5.0
    grd.N = 20
    print 'z_w =', grd.z_w[:,5,5]
    print 'h = ', grd.h[5,5]

def test_nc_grid():
    grd = nc_grid('/Users/rob/Projects/Merrimack/Grid/merrimack_large_grd.nc')
    print grd.__dict__.keys()
    
    grd = nc_grid('/Users/rob/Models/roms/roms-3.0/ocean_his.nc')
    print grd.__dict__.keys()

def test_write_roms_grid():
    """test write_roms_grid method of Grid class"""
    y, x = mgrid[0:1:100j, 0:1:100j]
    grd = Grid(x, y)
    grd.f = 1e-4
    grd.h = 10.0
    grd.write_roms_grid('cart_test.nc')
    print ' ### wrote cart_test.nc'
    
    lat, lon = mgrid[43:45:100j, -68:-70:100j]
    grdg = Grid(lon=lon, lat=lat)
    grdg.h = 10.
    grdg.f = 1.0e-4
    grdg.write_roms_grid('geo_test.nc')
    print ' ### wrote geo_test.nc'
    

def test_rho_to_vert():
    yv, xv = mgrid[0:20, 0:20:0.5]
    def rot2d(x, y, ang):
        'rotate vectors by geometric angle'
        xr = x*cos(ang) - y*sin(ang)
        yr = x*sin(ang) + y*cos(ang)
        return xr, yr
    
    print 'Verticies calculated from rho points: '
    for ang in arange(0, 2*pi, pi/8):
        xvr, yvr = rot2d(xv, yv, pi/4)
        grd = Grid(x=xvr, y=yvr)
        x, y = rho_to_vert(grd.x_rho, grd.y_rho, grd.pm, grd.pn, grd.angle)
        print '    All values close for rotation of %f?  %s' % \
              (ang, allclose(x, xvr) and allclose(y, yvr))
    
    # pl.plot(xvr, yvr, '-k')
    # pl.plot(xvr.T, yvr.T, '-k')
    # pl.plot(x, y, '-r')
    # pl.plot(x.T, y.T, '-r')
    # pl.show()


if __name__ == '__main__':
    test_write_roms_grid()
    test_make_grid()
    test_gridgen()
    test_grid_3d()
    test_nc_grid()
    test_rho_to_vert()

