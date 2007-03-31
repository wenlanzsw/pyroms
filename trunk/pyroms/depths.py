#!/usr/bin/env python
# encoding: utf-8
"""
depths.py

Created by Rob Hetland on 2007-03-30.
Copyright (c) 2007 Texas A&M Univsersity. All rights reserved.
Release under MIT license.
"""

from numpy import *
from pyroms import Dataset

class Depths(object):
    """docstring for Depths"""
    def __init__(self, hc, h, theta_b, theta_s, N, nc=None):
        self.hc = hc
        self.h = h
        self.theta_b = theta_b
        self.theta_s = theta_s
        self.N = N
        self._nc = nc
    
    def _get_sc_w(self):
        return linspace(-1.0, 0.0, self.N+1)
    
    def _get_sc_r(self):
        sc_w = self._get_sc_w()
        return 0.5*(sc_w[1:]+sc_w[:-1])
    
    def _get_Cs_r(self):
        return (1-self.theta_b)*sinh(self.theta_s*self.sc_r)/sinh(self.theta_s) \
                + 0.5*self.theta_b*(tanh(self.theta_s*(self.sc_r+0.5))-tanh(0.5*self.theta_s))\
                  /tanh(0.5*self.theta_s)
    
    def _get_Cs_w(self):
        return (1-self.theta_b)*sinh(self.theta_s*self.sc_w)/sinh(self.theta_s) \
                + 0.5*self.theta_b*(tanh(self.theta_s*(self.sc_w+0.5))-tanh(0.5*self.theta_s))\
                  /tanh(0.5*self.theta_s)
    
    sc_w = property(_get_sc_w)
    sc_r = property(_get_sc_r)
    Cs_w = property(_get_Cs_w)
    Cs_r = property(_get_Cs_r)
    
    def rho(self, tidx=None, zeta=None):
        """docstring for r"""
        if tidx is None and zeta is None:   # default is zeta=0
            zeta = zeros(self.h.shape, 'd')
        elif tidx is not None:              # read in tidx values of zeta
            zeta = self._nc.variables['zeta'][tidx]
        
        hdim = ndim(self.h)
        assert zeta.shape[-hdim:] == self.h.shape, \
            'zeta and h must have the same spatial dimensions'
        
        if ndim(self.h) == ndim(zeta):       # Assure a time-dimension exists
            zeta = zeta[newaxis, :]
        
        ti = zeta.shape[0]
        z = empty((ti, self.N) + self.h.shape, 'd')
        for n in range(ti):
            for  k in range(self.N):
                z0=(self.sc_r[k]-self.Cs_r[k])*self.hc + self.Cs_r[k]*self.h;
                z[n,k,:] = z0 + zeta[n,:]*(1.0 + z0/self.h);
        
        return(squeeze(z))
    
    def w(self, tidx=None, zeta=None):
        """docstring for r"""
        if tidx is None and zeta is None:   # default is zeta=0
            zeta = zeros(self.h.shape, 'd')
        elif tidx is not None:              # read in tidx values of zeta
            zeta = self._nc.variables['zeta'][tidx]
        
        hdim = ndim(self.h)
        assert zeta.shape[-hdim:] == self.h.shape, \
            'zeta and h must have the same spatial dimensions'
        
        if ndim(self.h) == ndim(zeta):       # Assure a time-dimension exists
            zeta = zeta[newaxis, :]
        
        ti = zeta.shape[0]
        z = empty((ti, self.N+1) + self.h.shape, 'd')
        for n in range(ti):
            for  k in range(self.N+1):
                z0=(self.sc_w[k]-self.Cs_w[k])*self.hc + self.Cs_w[k]*self.h;
                z[n,k,:] = z0 + zeta[n,:]*(1.0 + z0/self.h);
        
        return(squeeze(z))
    
    def dz_r(self, tidx=None, zeta=None):
        """docstring for dz_r"""
        return diff(self.rho(tidx, zeta), axis=0)
    
    def dz_w(self, tidx=None, zeta=None):
        """docstring for dz_w"""
        return diff(self.w(tidx, zeta), axis=0)
    

def nc_depths(nc, h=None):
    nc = Dataset(nc)
    hc = nc.variables['hc'][:]
    if h is None:
        h = nc.variables['h'][:]
    theta_b = nc.variables['theta_b'][:]
    theta_s = nc.variables['theta_s'][:]
    N = len(nc.dimensions['N'])
    return Depths(hc, h, theta_b, theta_s, N, nc)
