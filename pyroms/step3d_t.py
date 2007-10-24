#!/usr/bin/env python
# encoding: utf-8
"""
step3d_t.py

Created by Rob Hetland on 2007-10-24.
Copyright (c) 2007 Texas A&M Univsersity. All rights reserved.
Release under MIT license.
"""

from numpy import *
import pylab as pl
import pyroms
from scipy.interpolate import interp1d

import _step3d_t

class Step3d_t(object):
    """docstring for Step3d_t"""
    def __init__(self, nc, grd=None):
        self.nc = pyroms.Dataset(nc)
        self.dt = nc.variables['dt'][:]
        if grd is None:
            self.grd = pyroms.nc_grid(nc)
        else:
            self.grd = grd
        self.rmask = self.grd.mask
        self.pm = self.grd.pm
        self.pn = self.grd.pn
        self.t = pyroms.ocean_time(nc)
        self.z_w = pyroms.nc_depths(nc, grid='w')
    
    def dyn_step(self, tidx):
        """docstring for fname"""
        tidx=slice(tidx,tidx+2)
        tstart = self.t[tidx][0]
        tend = self.t[tidx][-1]
        
        z_wi = self.z_w[tidx]
        ui = self.nc.variables['u'][tidx]
        vi = self.nc.variables['v'][tidx]
        AKti = self.nc.variables['AKt'][tidx]
        
        N = ceil((tend-tstart)/self.dt)
        dt = (tend-tstart)/N
        
        for ittr, time in enumerate(linspace(tstart, tend, N)):
            w0 = (tend-time)/(tend-tstart)
            w1 = 1.0 - w0
            self.trc = _step3d_t.step3d_t(dt, 
                              self.rmask.T, self.pm.T, self.pn.T, 
                              (w0*z_wi[0] + w1*z_wi[1]).T,
                              (w0*AKti[0] + w1*AKti[1]).T,
                              (w0*ui[0] + w1*ui[1]).T,
                              (w0*vi[0] + w1*vi[1]).T,
                              self.trc.T).T
            print '[%d/%d] %9.4f' % (ittr, N, time)
        
    def static_step(self, tidx, N):
        """docstring for fname"""
        
        z_w = self.z_w[tidx]
        u = self.nc.variables['u'][tidx]
        v = self.nc.variables['v'][tidx]
        AKt = self.nc.variables['AKt'][tidx]
        
        for ittr, time in enumerate(linspace(self.dt, self.dt*N, N)):
            self.trc = _step3d_t.step3d_t(self.dt,
                              self.rmask.T, self.pm.T, self.pn.T, 
                              z_w.T, AKt.T, u.T, v.T, self.trc.T).T
            print '[%d/%d] %9.4f' % (ittr, N, time)


if __name__ == '__main__':
    nc = pyroms.MFDataset('/Volumes/Hetland_merrimack/field_2007/steady/ocean_1750*.nc')
    step = Step3d_t(nc)
    step.dt=5.0
    y, x = mgrid[0:258, 0:514]
    step.trc = sin(x*8*pi/x.max())[newaxis,:,:] * ones((30,), 'd')[:, newaxis, newaxis]
    step.static_step(17, 500)



