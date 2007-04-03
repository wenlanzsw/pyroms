"""Tools to work with the Regional Ocean Modeling System (ROMS).

Based on:
 + NumPy (http://numpy.scipy.org)
 + matplotlib with the basemap toolkit (http://matplotlib.sourceforge.net)
 + netCDF4 (http://www.cdc.noaa.gov/people/jeffrey.s.whitaker/python/netCDF4.html)

Contains:
 + grid
    + Grid - a class for holding horizontal and vertical ROMS grid information
    + nc_grid - make an instance of Grid based on a netcdf file
    + gridgen - create an instance of Grid using the gridgen program
 + roms_time - retrieve time information from a netcdf file
 + roms_tools - functions for working with roms coordinates
"""

classifiers = """\
Development Status :: beta
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: MIT
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries :: Python Modules
"""


from numpy.distutils.core import Extension

ext1 = Extension(name = '_iso',
                 sources = ['pyroms/iso.f'])

doclines = __doc__.split("\n")

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "pyroms",
          version = '0.7.0',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          author = "Robert Hetland",
          author_email = "hetland@tamu.edu",
          url = "http://code.google.com/p/pyroms/",
          packages = ['pyroms'],
          license = 'MIT',
          platforms = ["any"],
          ext_modules = [ext1,],
          classifiers = filter(None, classifiers.split("\n")),
          )
    