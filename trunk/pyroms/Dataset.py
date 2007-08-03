# encoding: utf-8

# Dataset
# Copyright (c) 2006 Robert Hetland

try:
    import netCDF4
    import netCDF4
    try:
        import MFnetCDF4
    except:
        import MFnetCDF4_classic as MFnetCDF4
    
    def Dataset(ncfile):
        """Return an appropriate netcdf object:
                netCDF4 object given a file string
                MFnetCDF4 object given a list of files
            
            A netCDF4 or MFnetCDF4 object returns itself."""
        if isinstance(ncfile, str):
            return netCDF4.Dataset(ncfile, 'r')
        elif isinstance(ncfile, list) or isinstance(ncfile, tuple):
            return MFnetCDF4.Dataset(sorted(ncfile))
        elif hasattr(ncfile, 'variables'):  # accept any oject with a variables attribute
            assert isinstance(ncfile.variables, dict), \
                   'variables attribute must be a dictionary'
            return ncfile
        else:
            raise TypeError, 'type %s not supported' % type(ncfile)

    def MFDataset(ncfile):
        """Return an MFnetCDF4 object given a string or list.  A string is expanded
           with wildcards using glob.  A netCDF4 or MFnetCDF4 object returns itself."""
        if isinstance(ncfile, str):
            return MFnetCDF4.Dataset(ncfile)
        elif isinstance(ncfile, list) or isinstance(ncfile, tuple):
            return MFnetCDF4.Dataset(sorted(ncfile))
        elif hasattr(ncfile, 'variables'):  # accept any oject with a variables attribute
            assert isinstance(ncfile.variables, dict), \
                   'variables attribute must be a dictionary'
            return ncfile
        else:
            raise TypeError, 'type %s not supported' % type(ncfile)
            return MFnetCDF4.Dataset(files)

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
        