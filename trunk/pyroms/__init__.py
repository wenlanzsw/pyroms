# PYROMS class
# (c) Rob Hetland, 2005
''' 
PyROMS is a suite of tools for working with ROMS (the Regional Ocean Modeling System).

Required modules:
    NumPy, netcdf4-python

Optional modules (mainly for the Grid class):
    SciPy, matplotlib

Copyright (C) 2007, Robert Hetland, released under a BSD-style license.
'''

from Dataset import Dataset, MFDataset
from ocean_time import ocean_time
from velocity import nc_velocity

from roms_tools import zatr, zatw, scoordr, scoordw, isoslice, shrink, gc_dist, \
                       rot2d, zslice, iso_integrate, surface, N2, nc_N2, \
                       nc_gls_dissipation, nc_curl, nc_div, nc_pstrain, transect

from ocean import rho_stp, o2_sat
from polygeom import Polygeom
from polyclick import PolyClick
from grid import Grid, gridgen, nc_grid
from depths import Depths, nc_depths
from greatcircle import GreatCircle
from gshhs import gshhs
from step3d_t import step3d_t

__authors__ = ['Robert Hetland <hetland@tamu.edu>']
__version__ = '0.4.5'
