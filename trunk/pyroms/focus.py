'A class for creating focus elements for use with gridgen'

# Copyright R. Hetland on 2007-11-26.  
# All rights reserved.

# Version 1.0.0 on 2007-11-26
#   Initial version of Focus and FocusPoint classes.

from numpy import *
from scipy.special import erf


class FocusPoint(object):
    """
    Return a transformed, uniform grid, focused around point xo, yo, with a focusing
    factor of focus, and x and y extent given by Rx and Ry.  The region of focusing
    will be approximately Gausian, and the resolution will be increased by approximately
    the value of factor.  To achive focusing on a line in the x- or y-direction, use a
    large value for R in the desired direction; typically a value of 10.0 is good enough
    to have focusing only in one direction.
    """
    def __init__(self, xo, yo, factor=2.0, Rx=0.1, Ry=None):
        
        if Ry is None: Ry = Rx
        
        self.xo = xo
        self.yo = yo
        self.factor = factor
        self.Rx = Rx
        self.Ry = Ry
    
    def __call__(self, x, y):
        
        x = asarray(x)
        y = asarray(y)
        assert not any(x>1.0) or not any(x<0.0) or not any(y>1.0) or not any(x<0.0), \
                'x and y must both be within the range [0, 1].'
        
        alpha = 1.0 - 1.0/self.factor
        def xf(x, y):
            return x - 0.5*sqrt(pi)*self.Rx*alpha*exp(-(y-self.yo)**2/self.Ry**2)*erf((x-self.xo)/self.Rx)
        
        def yf(x, y):
            return y - 0.5*sqrt(pi)*self.Ry*alpha*exp(-(x-self.xo)**2/self.Rx**2)*erf((y-self.yo)/self.Ry)
        
        xf0 = xf(0.0, y); xf1 = xf(1.0, y)
        yf0 = yf(x, 0.0); yf1 = yf(x, 1.0)
        
        return (xf(x, y)-xf0)/(xf1-xf0), (yf(x, y)-yf0)/(yf1-yf0)


class Focus(object):
    """
    foc = Focus(xo, yo, factor=2.0, Rx=0.1, Ry=Rx)
    
    Return a transformed, uniform grid, focused around point xo, yo, with a focusing
    factor of focus, and x and y extent given by Rx and Ry.  The region of focusing
    will be approximately Gausian, and the resolution will be increased by approximately
    the value of factor.  To achive focusing on a line in the x- or y-direction, use a
    large value for R in the desired direction; typically a value of 10.0 is good enough
    to have focusing only in one direction.
    
    Calls to the object return transformed coordinates:
        xf, yf = foc(x, y)
    where x and y must be within [0, 1], and are typically a uniform, normalized grid.
    
    Additional focus points may be added with the add_focus_point method:
        foc.add_focus_point(xo, yo, factor=2.0, Rx=0.1, Ry=Rx)
    subsequent calls to foc will result in a transformed grid with all the focus points
    applied in sequence.
    
    EXAMPLE:
    
    foc = Focus(0.7, 0.0, factor=2.0, Rx=0.2, Ry=10)
    foc.add_focus_point(0.2, 0.7, factor=3.0, Rx=0.2)
    
    y, x = mgrid[0:1:50j, 0:1:70j]
    xf, yf = foc(x, y)
    """
    def __init__(self, xo, yo, factor=2.0, Rx=0.1, Ry=None):
        self._focuspoints = []
        self.add_focus_point(xo, yo, factor, Rx, Ry)
    
    def add_focus_point(self, xo, yo, factor=2.0, Rx=0.1, Ry=None):
        """docstring for add_point"""
        self._focuspoints.append(FocusPoint(xo, yo, factor, Rx, Ry))
    
    def __call__(self, x, y):
        """docstring for __call__"""
        for focuspoint in self._focuspoints:
            x, y = focuspoint(x, y)
        return x, y


if __name__ == '__main__':
    import pylab as pl
    
    y, x = mgrid[0:1:50j, 0:1:70j]
    
    factor = 2.0
    foc = Focus(0.7, 0.0, factor=factor, Rx=0.2, Ry=10)
    foc.add_focus_point(0.2, 0.7, factor=factor, Rx=0.2)
    # foc.add_focus_point(0.7, 0.7, factor=factor, Rx=0.2)
    foc.add_focus_point(0.2, 0.2, factor=factor, Rx=0.2)
    xf, yf = foc(x, y)
    
    dx = diff(xf, axis=-1)
    print 'Focusing factor in x:', dx.max()/dx.min()
    
    dy = diff(yf, axis=0)
    print 'Focusing factor in y:', dy.max()/dy.min()
    
    print 'should both be approximately: ', factor
    
    pl.plot(xf, yf, '-k')
    pl.plot(xf.T, yf.T, '-k')
    pl.show()