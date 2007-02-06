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

from distutils.core import setup

doclines = __doc__.split("\n")

setup(name = "froms",
      version = '0.1.0',
      description = doclines[0],
      long_description = "\n".join(doclines[2:]),
      author = "Robert Hetland",
      author_email = "hetland@tamu.edu",
      url = "http://pong.tamu.edu/~rob/python",
      download_url = "http://pong.tamu.edu/~rob/python/roms.tar.gz",
      packages = ['pyroms'],
      license = 'MIT',
      platforms = ["any"],
      classifiers = filter(None, classifiers.split("\n")),
     )
