#!/usr/bin/env python
# encoding: utf-8
"""
velocity class.py

Created by Rob Hetland on 2007-01-13.
Copyright (c) 2007 Texas A&M Univsersity. All rights reserved.
"""

import roms
from numpy import *

class velocity (object):
    
    def __init__(self, nc):
        self.nc = roms.Dataset(nc)
        self.u = self.nc.variables['u']
        self.v = self.nc.variables['v']
        try:
            self.ang = self.nc.variables['angle'][:,:]
        except:
            self.ang = None
    
    def __getitem__(self, elem):
        u = self.u[elem]
        v = self.v[elem]
        u, v = roms.shrink(u, v)
        if self.ang is not None:
            ang = roms.shrink(self.ang, u.shape)
            return rot2d(u, v, ang)
        else:
            return u, v



if __name__ == '__main__':
    nc = roms.Dataset('/Users/rob/Models/roms/inertial/ocean_his.nc')
    vel = roms_vel(nc)
    u, v = vel[-1, -1, :]
    def rot2d(x, y, ang):
        'rotate vectors by geometric angle'
        xr = x*cos(ang) - y*sin(ang)
        yr = x*sin(ang) + y*cos(ang)
        return xr, yr
    
    
