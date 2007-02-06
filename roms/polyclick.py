#!/usr/bin/env python
# encoding: utf-8
"""
polyclick.py

Created by Rob Hetland on 2007-01-21.
Copyright (c) 2007 Texas A&M Univsersity. All rights reserved.
"""

from numpy import *
import pylab as pl
from matplotlib.patches import Polygon, PolygonInteractor

class PolyClick(object):
    """
    p = PolyClick()
    
    Interactive polygon creation, and editing.  Create an instance, and start defining the
    initial polynomial.  At this point, vertecies are blue.  Switch to editing mode by hitting
    return.  This changes the vertecies to red, and you can modify them as per PolygonInteractor.
    (i.e., hit 'i' to insert a point, 'd' to delete a point, and 't' to toggle the vertex visibility)
    
    Data are always available in p.x, p.y, and p.verts attributes.
    """
    
    
    def _on_key(self, event):
        if event.key is None:
            # Get rid of line, and old mpl connections.
            self._line.set_visible(False)
            self._ax.figure.canvas.mpl_disconnect(self._key_id)
            self._ax.figure.canvas.mpl_disconnect(self._click_id)
            
            # start a polygon intteractor
            verts = zip(self._xdata, self._ydata)
            poly = Polygon(verts, alpha=0.2, fc='k')
            self._ax.add_patch(poly)
            self._pi = PolygonInteractor(poly)
            self._ax.add_line(self._pi.line)
            self.get_xdata = self._pi.line.get_xdata
            self.get_ydata = self._pi.line.get_ydata
    
    def _on_click(self, event):
        self._xdata.append(event.xdata)
        self._ydata.append(event.ydata)
        self._line.set_data(self._xdata, self._ydata)
        self._ax.figure.canvas.draw_idle()
    
    def __init__(self, xdata=[], ydata=[], ax=None):
        
        if ax is None: ax = pl.gca()
        self._ax = ax
        
        self._xdata = list(xdata); self._ydata = list(ydata)
        self._pi = None
        self._line = pl.Line2D(self._xdata,self._ydata,marker='o', markerfacecolor='b')
        self._ax.add_line(self._line)
        
        self._key_id = self._ax.figure.canvas.mpl_connect('key_press_event', self._on_key)
        self._click_id = self._ax.figure.canvas.mpl_connect('button_press_event', self._on_click)
    
    def get_xdata(self):
        if self._pi is None:
            return self._xdata
        else:
            return self._pi.line.get_xdata() 
    
    def get_ydata(self):
        if self._pi is None:
            return self._ydata
        else:
            return self._pi.line.get_ydata() 

    def get_verts(self):
        return zip(self.get_xdata(), self.get_ydata())

    x = property(get_xdata)
    y = property(get_ydata)
    verts = property(get_verts)



if __name__ == '__main__':
    pass

