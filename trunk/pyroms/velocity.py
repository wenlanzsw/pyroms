#!/usr/bin/env python
# encoding: utf-8
"""
velocity class.py

Created by Rob Hetland on 2007-01-13.
Copyright (c) 2007 Texas A&M Univsersity. All rights reserved.
"""

import pyroms
from numpy import *

class velocity (object):
    
    def __init__(self, nc, grid='psi'):
        self.nc = pyroms.Dataset(nc)
        self.grid = grid
        self.u = self.nc.variables['u']
        self.v = self.nc.variables['v']
        try:
            self.ang = self.nc.variables['angle'][:,:]
        except:
            self.ang = None
    
    def __getitem__(self, elem):
        u = self.u[elem]
        v = self.v[elem]
        u, v = pyroms.shrink(u, v)
        if self.ang is not None:
            ang = pyroms.shrink(self.ang, u.shape)
            u, v = pyroms.rot2d(u, v, ang)
        
        if self.grid=='rho':
            shp = list(u.shape)
            shp[-1] -= 1
            shp[-2] -= 1
            shpr = list(u.shape)
            shpr[-1] += 1
            shpr[-2] += 1
            ur = zeros(shpr)
            vr = zeros(shpr)
            ur[...,1:-1,1:-1] = pyroms.shrink(u, shp)
            vr[...,1:-1,1:-1] = pyroms.shrink(v, shp)
            return ur, vr
        else:
            return u, v


if __name__ == '__main__':
    nc = pyroms.Dataset('/Users/rob/Models/roms/inertial/ocean_his.nc')
    vel = velocity(nc)
    u, v = vel[-1, -1, :]
    def rot2d(x, y, ang):
        'rotate vectors by geometric angle'
        xr = x*cos(ang) - y*sin(ang)
        yr = x*sin(ang) + y*cos(ang)
        return xr, yr
    
    
