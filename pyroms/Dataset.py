# encoding: utf-8

# Dataset
# Copyright (c) 2006 Robert Hetland

try:
    import netCDF4
    import netCDF4_classic
    import MFDataset
    
    def Dataset(ncfile):
        """Return an appropriate netcdf object given a file string, a list of files
           (returns a multicdf instance), or a netcdf object (returns itself)."""
        if isinstance(ncfile, str):
            return netCDF4.Dataset(ncfile, 'r')
        elif isinstance(ncfile, list) or isinstance(ncfile, tuple):
            return MFDataset.Dataset(sorted(ncfile))
        elif isinstance(ncfile, netCDF4.Dataset) or \
             isinstance(ncfile, MFDataset.MFDataset.Dataset) or \
             isinstance(ncfile, netCDF4_classic.Dataset):
            return ncfile
        else:
            raise TypeError, 'type %s not supported' % type(ncfile)
except:
    import pupynere
    import warnings
    
    warnings.warn('netCDF4 not found -- using pupynere.')
    
    def Dataset(ncfile):
        if isinstance(ncfile, str):
            return pupynere.NetCDFFile(ncfile)
        elif isinstance(ncfile, pupynere.NetCDFFile):
            return ncfile
        else:
            raise TypeError, 'type %s not supported' % type(ncfile)


if __name__ == '__main__':
    nc = Dataset('/Users/rob/Models/roms/roms-3.0/ocean_his.nc')
    nc1 = Dataset(nc)
    print nc.variables.keys()
    print nc1.variables.keys()
    
    print nc1 is nc
        