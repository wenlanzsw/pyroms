"""Tools to work with the Regional Ocean Modeling System (ROMS).

Based on:
 + NumPy (http://numpy.scipy.org)
 + netCDF4 (http://www.cdc.noaa.gov/people/jeffrey.s.whitaker/python/netCDF4.html)

Contains:
 + multicdf - a wrapper for netCDF4 that concatenates multiple netCDF files
 + cdftime - retrieve time information from a netcdf file
 + roms_tools - functions for working with roms coordinates
"""

classifiers = """\
Development Status :: alpha
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: BSD
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
      packages = ['roms'],
      license = 'BSD',
      platforms = ["any"],
      classifiers = filter(None, classifiers.split("\n")),
     )
