"""
Tools for working with ROMS model input and output.
"""

from numpy import *
from pyroms import Dataset
from matplotlib.toolkits.basemap.greatcircle import GreatCircle

gc_dist = vectorize(lambda lon1, lat1, lon2, lat2: \
                    GreatCircle(6378137.0, 6356752.3142, \
                                lon1, lat1, lon2, lat2).distance, \
                                doc= \
'''
Usage: dist = roms.gc_dist(lon1, lat1, lon2, lat2)
       All inputs can be arrays, and obey array broadcasting.
''')

def nanmask(a):
    return ma.masked_where(isnan(a), a)

def zatr(ncfile, time=False):
    """
    ZATR finds z at rho points (positive up, zero at rest surface) 
    If zeta = True for all times with the calculated value of zeta, 
    otherwise, for one time calculated with zeta = 0 (default)
    
    ncfile - NetCDF file to use (a ROMS history file or netcdf object).
    """
    nc = Dataset(ncfile)
    h = nc.variables['h'][:]
    h = atleast_2d(h)
    hc = nc.variables['hc'][:]
    Cs_r = nc.variables['Cs_r'][:]
    try:
        sc_r = nc.variables['sc_r'][:]  # roms-2.x format
    except:
        sc_r = nc.variables['s_rho'][:]  # roms-3.x format
    N = len(nc.dimensions['N'])
    if time==False:
        zeta = zeros(h.shape, 'd')
        zeta = zeta[newaxis, :]
    else:
        zeta = nc.variables['zeta'][:]
    ti = zeta.shape[0]
    z = empty((ti, N) + h.shape, 'd')
    for n in arange(ti):
        for  k in arange(N):
            z0=(sc_r[k]-Cs_r[k])*hc + Cs_r[k]*h;
            z[n,k,:] = z0 + zeta[n,:]*(1.0 + z0/h);
    return(squeeze(z))

def zatw(ncfile, time=False):
    """
    ZATW finds z at w-points (positive up, zero at rest surface) 
    If zeta = True for all times with the calculated value of zeta, 
    otherwise, for one time calculated with zeta = 0 (default)
    
    ncfile - NetCDF file to use (a ROMS history file or netcdf object).
    """
    nc = Dataset(ncfile)
    h = nc.variables['h'][:]
    hc = nc.variables['hc'][:]
    Cs_w = nc.variables['Cs_w'][:]
    try:
        sc_w = nc.variables['sc_w'][:]  # roms-2.x format
    except:
        sc_w = nc.variables['s_w'][:]  # roms-3.x format
    if time==False:
        zeta = zeros(h.shape, 'd')
        zeta = zeta[newaxis, :]
    else:
        zeta = nc.variables['zeta'][:]
    N = len(nc.dimensions['N']) + 1
    ti = zeta.shape[0]
    z = empty((ti, N) + h.shape, 'd')
    for n in arange(ti):
        for  k in arange(N):
            z0=(sc_w[k]-Cs_w[k])*hc + Cs_w[k]*h;
            z[n,k,:] = z0 + zeta[n,:]*(1.0 + z0/h);
    return(squeeze(z))

def scoordr(h,hc,theta_b,theta_s,N):
    """
    z=scoordr(h,hc,theta_b,theta_s,N)
    SCOORDR finds z at rho points (positive up, zero at rest surface) 
     h = array of depths (e.g., from grd file)
     hc = critical depth
     theta_b = surface/bottom focusing parameter
     theta_s = strength of focusing parameter
     N number of vertical rho-points
    """
    sc_w = arange(-1, 1./N, 1./N, dtype='d')
    sc_r = 0.5*(sc_w[1:]+sc_w[:-1])
    Cs_r = (1-theta_b)*sinh(theta_s*sc_r)/sinh(theta_s)\
          +0.5*theta_b\
          *(tanh(theta_s*(sc_r+0.5))-tanh(0.5*theta_s))/tanh(0.5*theta_s)
    z_r = empty((N,) + h.shape, dtype='d')
    for  k in arange(N):
        z_r[k,:]=(sc_r[k]-Cs_r[k])*hc + Cs_r[k]*h
    return(squeeze(z_r))

def scoordw(h,hc,theta_b,theta_s,N):
    """
    z=scoordr(h,hc,theta_b,theta_s,N)
    SCOORDW finds z at rho points (positive up, zero at rest surface) 
     h = array of depths (e.g., from grd file)
     hc = critical depth
     theta_b = surface/bottom focusing parameter
     theta_s = strength of focusing parameter
     N number of vertical w-points (one more than rho-points)
    """
    sc_w = arange(-1, 1./N, 1./N, dtype='d')
    Cs_w = (1-theta_b)*sinh(theta_s*sc_w)/sinh(theta_s)\
          +0.5*theta_b\
          *(tanh(theta_s*(sc_w+0.5))-tanh(0.5*theta_s))/tanh(0.5*theta_s)
    z_w = empty((N,) + h.shape, dtype='d')
    for  k in arange(N):
        z_w[k,:]=(sc_w[k]-Cs_w[k])*hc + Cs_w[k]*h
    return(squeeze(z_w))

def isoslice(var,prop,isoval=0,masking=True):
    """
    result = isoslice(variable,property[,isoval=0])
    
    result is a a projection of variable at property == isoval in the first
    nonsingleton dimension.  In the case when there is more than one zero
    crossing, the results are averaged.
    
    EXAMPLE:
    Assume two three dimensional variable, s (salt), z (depth), and
    u (velicity), all on the same 3D grid.  x and y are the horizontal 
    positions, with the same horizontal dimensions as z (the 3D depth 
    field).  Here, assume the vertical dimension, that will be projected,
    is the first.  
    
    s_at_m5  = isoslice(s,z,-5);        # s at z == -5
    h_at_s30 = isoslice(z,s,30);       # z at s == 30
    u_at_s30 = isoslice(u,s,30);       # u at s == 30
    """
    if (len(var.squeeze().shape)<=2):
        raise ValueError, 'variable must have at least two dimensions'
    if not prop.shape == var.shape:
        raise ValueError, 'dimension of var and prop must be identical'
    prop=prop-isoval
    sz = shape(var)
    var = var.reshape(sz[0],-1)
    prop = prop.reshape(sz[0],-1)
    #find zero-crossings (zc == 1)
    zc =  where( (prop[:-1,:]*prop[1:,:])<0 ,1., 0.)
    varl = var[:-1,:]*zc
    varh = var[1:,:]*zc
    propl = prop[:-1,:]*zc
    proph = prop[1:,:]*zc
    result = varl - propl*(varh-varl)/(proph-propl)
    result = where(zc==1., result, 0.)
    szc = zc.sum(axis=0)
    szc = where(szc==0., 1, szc)
    result = result.sum(axis=0)/szc
    if masking:
        result = ma.masked_where(zc.sum(axis=0)==0, result)
        if all(result.mask):
            raise Warning, 'property==%f out of range (%f, %f)' % \
                           (isoval, (prop+isoval).min(), (prop+isoval).max())
    result = result.reshape(sz[1:])
    return(result)

def shrink(a,b):
    """
    Return array(s) shrunk to fit a specified shape (or each other).
    
    a = shrink(array, shape)
        array is an numpy ndarray, and shape is a tuple (e.g., from array.shape).
        a is the input array shrunk such that its maximum dimensions are given
        by shape.  If shape has more dimensions than array, the last dimensions
        of shape are fit.
    
    as, bs = shrink(a, b)
        If the second argument is also an array, both a and b are shrunk to the
        dimensions of each other.  The input arrays must have the same number of
        dimensions, and the resulting arrays will have the same shape.
    
    Example:
        >>> shrink(rand(10, 10), (5, 9, 18)).shape
        (9, 10)
        >>> map(shape, shrink(rand(10, 10, 10), rand(5, 9, 18)))        
        [(5, 9, 10), (5, 9, 10)]      
    """
    if isinstance(b, ndarray):
        if not len(a.shape) == len(b.shape):
            raise Exception, 'input arrays must have the same number of dimensions'
        a = shrink(a,b.shape)
        b = shrink(b,a.shape)
        return (a, b)
    if isinstance(b, int):
        b = (b,)
    if len(a.shape) == 1:                   # consider a 1D array as a special case
        dim = b[-1]
        while a.shape[0] > dim:             # only shrink a
            if (dim - a.shape[0]) >= 2:     # trim off edges evenly
                a = a[1:-1]
            else:                           # or average adjacent cells
                a = 0.5*(a[1:] + a[:-1])     
        return a
    for dim_idx in range(-(len(a.shape)),0):
        dim = b[dim_idx]
        a = a.swapaxes(0,dim_idx)           # put working dim first
        while a.shape[0] > dim:             # only shrink a
            if (a.shape[0] - dim) >= 2:     # trim off edges evenly
                a = a[1:-1,:]
            if (a.shape[0] - dim) == 1:      # or average adjacent cells
                a = 0.5*(a[1:,:] + a[:-1,:])
        a = a.swapaxes(0,dim_idx)           # swap working dim back
    return a

def arg_nearest(x, xo, scale=None):
    """
    idx = arg_nearest(x, xo, scale=None)
    returns a tuple of the indices of the N-dimensional arrays closest to
    the N-dimensional vector, xo.  Scale can be used to normalize arrays
    if that is an issue.  All input arrays must have the same dimension,
    but can use singleton dimensions (e.g., newaxis) as placeholders.
    """
    assert len(x) == len(xo), 'must have %d %d-dimensional arrays' % (len(xo), len(xo))
    for dim, p in enumerate(x):
        assert len(p.shape) == len(xo), 'dimensions of array %d must be %d' % (dim+1, len(xo))
    if scale:
        assert len(x) == len(scale), 'dimensions of array and scale must be identical'
        x = [p*q for (p, q) in zip(x, scale)]
        xo = [p*q for (p, q) in zip(xo, scale)]
    q = reduce(add, [(y-yo)**2 for (y, yo) in zip(x, xo)])
    return where(q == q.min())


if __name__ == '__main__':
    '''Some simple tests'''
    a = rand(5, 5)
    print a
    print '-='*40
    print extrapolate(a)
    print '-='*40
    print extrapolate(a, bc='constant')
    print '-='*40
    print extrapolate(a, axis=1)
    print '-='*40
    print extrapolate(a, axis=1, bc='constant')

if 0:
    '''More tests...'''
    a = arange(40.)
    a = a.reshape(5,-1)
    print ' ### original matrix'
    print a
    print ' ### interior points'
    print shrink(a,(3,6))
    print ' ### avg in first dim'
    print shrink(a,(4, 8))
    print ' ### avg in second dim'
    print shrink(a,(5, 7))
    
    print '+-'*20, ' testing shapes'
    print shrink(rand(10, 10, 10),(10,12,9)).shape == (10,10,9)
    print shrink(rand(100),99).shape == (99,)
    a, b = shrink(rand(10, 12, 15), rand(12, 10, 15))
    print a.shape == (10,10,15), a.shape == b.shape
    print '+-'*20, ' testing arrays'
    a = rand(10,12,14)
    a_test = 0.5*(a[1:,:]+a[:-1,:])
    print all(a_test == shrink(a,a_test.shape))
    a_test = a_test[:,1:-1,:]
    print all(a_test == shrink(a,a_test.shape))
    a_test = a_test[:,:,2:-2]
    print all(a_test == shrink(a,a_test.shape))
    print '+-'*20, ' testing 1D arrays'
    a = rand(100)
    b = rand(99)
    ab = 0.5*(a[1:]+a[:-1])
    aa, bb = shrink(a, b)
    print all(aa == ab), aa.shape == bb.shape
    aa = shrink(a, 99)
    print all(aa == ab), aa.shape == bb.shape

if 0:
    '''Really, you'd think it would work by now...'''
    import pylab as pl
    z,y,x = mgrid[0:1:30j,0:1:100j,0:1:100j]
    p = sqrt((x-0.5)**2 + (y-0.5)**2 + (z-1.0)**2)
    # kind of cool
    pz = isoslice(p,z,0.5)
    # more cool (p does not have a value of 0.5 everywhere)
    pp = isoslice(z,p,0.5)
    fig = pl.figure(1, figsize=(8,4))
    ax = fig.add_subplot(121)
    ax.pcolor(x[1], y[1], pz, shading='flat')
    ax=fig.add_subplot(122)
    ax.pcolor(x[1], y[1], pp, shading='flat')
    pl.show()

    
