# multicdf class
#
# Copyright (c) 2006 Robert Hetland
#
# TODO:
#   Make sure that when there is a single index in a slice out of a file, add 
#   a new dimension to that variable so concatination works
#

from numpy import *
import netCDF4, glob

class multivar (object):
    
    def __init__(self, name, ncobjects):
        self._name = name
        self._nc = ncobjects
        self._master = ncobjects[0]
        self._lengths = []
        for nc in self._nc:
            self._lengths.append(nc.variables[self._name]._shape()[0])
        
    def ncattrs(self):
        return self._master.variables[self._name].ncattrs()
    
    def __getattr__(self, name):
        if name[0] != '_':
            return self._master.variables[self._name].__getattribute__(name)
        else:
            return self.__dict__[name]
    
    def _broadcast_index(self, lengths, elem):
        idxa = array(concatenate(map(arange, lengths))[elem])
        idx = []
        for n in range(len(lengths)):
            idx.append(zeros(lengths[n])+n)
        idxn = concatenate(idx)[elem]
        if isinstance(elem, int):
            elems = len(lengths)*[None]
            elems[idxn] = int(idxa)
        else:
            elems = []
            for n in range(len(lengths)):
                idxtmp = idxa[idxn==n]
                if len(idxtmp)==0:
                    elems.append(None)
                else:
                    elems.append(slice(min(idxtmp),max(idxtmp)+1,elem.step))
        return elems
    
    def __getitem__(self, elem):
        if isinstance(elem,slice):   # Singleton slices
            elem = [elem]
        else:                        # Multi-dimension slice
            elem = list(elem)
        # broadcast first (time) index over files
        idx = self._broadcast_index(self._lengths,elem[0])
        # get the requested time slice from each file
        varout = []; maxdim = 1
        for ncidx, nc in enumerate(self._nc):
            elem[0] = idx[ncidx]
            if elem[0] != None:
                vartmp = asarray(nc.variables[self._name][tuple(elem)])
                maxdim = max(maxdim, len(vartmp.shape))
                varout.append(vartmp)
        # make sure that the dimensions are similar for concatenation
        for ncidx, var in enumerate(varout):
            if (len(var.shape) < maxdim) and var is not None:
                if var.shape == ():
                    varout[ncidx] = varout[ncidx][NewAxis]
                else:
                    varout[ncidx] = varout[ncidx][NewAxis]
        return concatenate(varout)
    

class multicdf (object):
    
    def __init__(self, files):
        self._nc = []
        for file in files:
            try:
                self._nc.append(netCDF4.Dataset(file))
            except:
                raise Exception, 'error opening %s' % file
        self._master = self._nc[0]
        # The dimensions is just a copy of the dimension list in the master file
        self.dimensions = self._master.dimensions
        # Find the first unlimited dimension (there should be only one)
        for key in self.dimensions:
            if self.dimensions[key].isunlimited():
                self.unlimdim = key
                break
        # We need a copy of self._master.variables, since it is modified later
        self.variables = self._master.variables.copy()
        # Change all the variables with unlimited dimension to a multivar class
        for key in self.variables:
            if self.unlimdim in self.variables[key].dimensions:
                self.variables[key] = multivar(key, self._nc)
    
    def ncattrs(self):
        return self._master.ncattrs()
    
    def __getattr__(self, name):
        """Return attributes from the master file.  Does not yet include
           methds like self.ncattr(), but does return file attributes."""
        return self._master.__getattribute__(name)
    

def Dataset(ncfile):
    """Return an appropriate netcdf object given a file string, a list of files
       (returns a multicdf instance), or a netcdf object (returns itself)."""
    if isinstance(ncfile,str):
        nc = netCDF4.Dataset(ncfile, 'r')
    elif isinstance(ncfile, list) or isinstance(ncfile, tuple):
        nc = multicdf(ncfile)
    elif isinstance(ncfile, netCDF4.Dataset) or isinstance(ncfile, multicdf):
        nc = ncfile
    else:
        raise TypeError, 'type %s not supported' % type(ncfile)
    return nc


