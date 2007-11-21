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
 + ocean_time - retrieve time information from a netcdf file
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

iso = Extension(name = '_iso',
                sources = ['pyroms/iso.f'])

step3d_t = Extension(name = '_step3d_t',
                     sources = ['pyroms/step3d_t.f'])

delaunay = Extension(name = '_delaunay',
                     sources=["pyroms/delaunay/_delaunay.cpp",
                              "pyroms/delaunay/VoronoiDiagramGenerator.cpp",
                              "pyroms/delaunay/delaunay_utils.cpp",
                              "pyroms/delaunay/natneighbors.cpp"])

doclines = __doc__.split("\n")

gshhs_datafiles = ['gshhs-data/gshhs_c.b', 
                   'gshhs-data/gshhs_l.b',
                   'gshhs-data/gshhs_i.b',
                   'gshhs-data/gshhs_h.b',
                   'gshhs-data/gshhs_f.b']

package_data = {'pyroms': gshhs_datafiles}

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "pyroms",
          version = '0.7.0',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          author = "Robert Hetland",
          author_email = "hetland@tamu.edu",
          url = "http://code.google.com/p/pyroms/",
          packages = ['pyroms', 'pyroms/delaunay'],
          license = 'MIT',
          platforms = ["any"],
          ext_modules = [iso, step3d_t, delaunay],
          classifiers = filter(None, classifiers.split("\n")),
          package_data = package_data,
          )
    