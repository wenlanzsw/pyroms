# encoding: utf-8

# Dataset
# Copyright (c) 2006 Robert Hetland

import netCDF4
import netCDF4_classic
import MFDataset
import roms

def Dataset(ncfile):
    """Return an appropriate netcdf object given a file string, a list of files
       (returns a multicdf instance), or a netcdf object (returns itself)."""
    if isinstance(ncfile,str):
        nc = netCDF4.Dataset(ncfile, 'r')
    elif isinstance(ncfile, list) or isinstance(ncfile, tuple):
        nc = MFDataset.Dataset(sorted(ncfile))
    elif isinstance(ncfile, netCDF4.Dataset) or \
         isinstance(ncfile, MFDataset.MFDataset.Dataset) or \
         isinstance(ncfile, netCDF4_classic.Dataset):
        nc = ncfile
    else:
        raise TypeError, 'type %s not supported' % type(ncfile)
    return nc

