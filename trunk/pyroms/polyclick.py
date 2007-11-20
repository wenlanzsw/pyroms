#!/usr/bin/env python
# encoding: utf-8
"""
polyclick.py

Created by Rob Hetland on 2007-01-21.
Copyright (c) 2007 Texas A&M Univsersity. All rights reserved.
"""

from numpy import *
import pylab as pl
from matplotlib.artist import Artist
from matplotlib.patches import Polygon, CirclePolygon
from matplotlib.lines import Line2D
from matplotlib.numerix import sqrt, nonzero, equal, asarray, dot, Float
from matplotlib.numerix.mlab import amin
from matplotlib.mlab import dist_point_to_segment


class BoundaryInteractor:
    """
    An polygon editor, modified from the matplotlib example to include beta values for use with griggen.

    Key-bindings

      't' toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them

      'd' delete the vertex under point

      'i' insert a vertex at point.  You must be within epsilon of the
          line connecting two existing vertices

      'p' plus -- mark a location as a positive (counterclockwise) corner
      
      'm' minus -- mark a location as a negative (clockwise) corner
      
    """

    showverts = True
    showbetas = True
    epsilon = 5  # max pixel distance to count as a vertex hit

    def __init__(self, ax, poly, beta=None):
        if poly.figure is None:
            raise RuntimeError('You must first add the polygon to a figure or canvas before defining the interactor')
        self.ax = ax
        self.canvas = poly.figure.canvas
        self.poly = poly

        x, y = asarray(zip(*self.poly.xy))
        self.line = Line2D(x,y,marker='o', markerfacecolor='k', animated=True)
        self.ax.add_artist(self.line)

        if beta is None:
            self.beta = zeros(len(x))
        else:
            self.beta = asarray(beta)
        
        self.pline = Line2D([], [], marker='^', markersize=20,
                            markerfacecolor='g', animated=True, lw=0)
        self.ax.add_artist(self.pline)
        if any(self.beta==1.0):
            self.pline.set_data(x[self.beta==1.0], y[self.beta==1.0])
        
        self.mline = Line2D([], [], marker='v', markersize=20,
                            markerfacecolor='r', animated=True, lw=0)
        self.ax.add_artist(self.mline)
        if any(self.beta==-1.0):
            self.pline.set_data(x[self.beta==-1.0], y[self.beta==-1.0])
        
        
        cid = self.poly.add_callback(self.poly_changed)
        self._ind = None # the active vert

        self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        
    
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.poly)
        if any(self.beta==1.0):
            self.ax.draw_artist(self.pline)
        if any(self.beta==-1.0):
            self.ax.draw_artist(self.mline)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)

    def poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, poly)
        self.line.set_visible(vis)  # don't use the poly visibility state


    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'
        x, y = zip(*self.poly.xy)

        # display coords
        xt, yt = self.poly.get_transform().numerix_x_y(x, y)
        d = sqrt((xt-event.x)**2 + (yt-event.y)**2)
        indseq = nonzero(equal(d, amin(d)))
        ind = indseq[0]

        if d[ind]>=self.epsilon:
            ind = None

        return ind

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts: return
        if event.inaxes==None: return
        if event.button != 1: return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts: return
        if event.button != 1: return
        self._ind = None

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes: return
        if event.key=='t':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts: self._ind = None
        if event.key=='T':
            self.showbetas = not self.showbetas
            self.pline.set_visible(self.showbetas)
            self.mline.set_visible(self.showbetas)
        elif event.key=='d':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                self.poly.xy = [tup for i,tup in enumerate(self.poly.xy) if i!=ind]
                self.beta = asarray([b for i, b in enumerate(self.beta) if i!=ind])
                self.line.set_data(zip(*self.poly.xy))
        elif event.key=='p':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                x, y = asarray(zip(*self.poly.xy))
                self.beta[ind] = 1.0
                self.pline.set_data(x[self.beta==1.0], y[self.beta==1.0])
                if any(self.beta==1.0):
                    self.ax.draw_artist(self.pline)
                if any(self.beta==-1.0):
                    self.ax.draw_artist(self.mline)
        elif event.key=='m':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                x, y = asarray(zip(*self.poly.xy))
                self.beta[ind] = -1.0
                self.mline.set_data(x[self.beta==-1.0], y[self.beta==-1.0])
                if any(self.beta==1.0):
                    self.ax.draw_artist(self.pline)
                if any(self.beta==-1.0):
                    self.ax.draw_artist(self.mline)
        elif event.key=='z':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                x, y = asarray(zip(*self.poly.xy))
                self.beta[ind] = 0.0
                if any(self.beta==1.0):
                    self.pline.set_data(x[self.beta==1.0], y[self.beta==1.0])
                else:
                    self.pline.set_data([], [])
                if any(self.beta==-1):
                    self.mline.set_data(x[self.beta==-1.0], y[self.beta==-1.0])
                else:
                    self.pline.set_data([], [])
                if any(self.beta==1.0):
                    self.ax.draw_artist(self.pline)
                if any(self.beta==-1.0):
                    self.ax.draw_artist(self.mline)
        elif event.key=='i':
            xys = self.poly.get_transform().seq_xy_tups(self.poly.xy)
            p = event.x, event.y # display coords
            for i in range(len(xys)-1):
                s0 = xys[i]
                s1 = xys[i+1]
                d = dist_point_to_segment(p, s0, s1)
                if d<=self.epsilon:
                    self.poly.xy.insert(i+1, (event.xdata, event.ydata))
                    self.line.set_data(zip(*self.poly.xy))
                    self.beta = list(self.beta)
                    self.beta.insert(i+1, 0)
                    self.beta = asarray(self.beta)
                    break


        self.canvas.draw()

    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showverts: return
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        x,y = event.xdata, event.ydata
        self.poly.xy[self._ind] = x,y
        
        x, y = asarray(zip(*self.poly.xy))
        self.line.set_data(x, y)

        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.poly)
        if any(self.beta==1.0):
            self.pline.set_data(x[self.beta==1.0], y[self.beta==1.0])
            self.ax.draw_artist(self.pline)
        if any(self.beta==-1.0):
            self.mline.set_data(x[self.beta==-1.0], y[self.beta==-1.0])
            self.ax.draw_artist(self.mline)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)


class BoundaryClick(object):
    """
    p = PolyClick()
    
    Interactive polygon creation, and editing.  Create an instance, and start defining the
    initial polynomial.  At this point, vertecies are blue.  Switch to editing mode by hitting
    return.  This changes the vertecies to red, and you can modify them as per PolygonInteractor.
    (i.e., hit 'i' to insert a point, 'd' to delete a point, and 't' to toggle the vertex visibility)
    
    Data are always available in p.x, p.y, and p.verts attributes.
    """
    
    def _boundary_interactor(self):
        """Send polyclick thread to boundary interactor"""
        # Get rid of line, and old mpl connections.
        self._line.set_visible(False)
        self._ax.figure.canvas.mpl_disconnect(self._key_id)
        self._ax.figure.canvas.mpl_disconnect(self._click_id)
    
        # start a polygon intteractor
        verts = zip(self._xdata, self._ydata)
        poly = Polygon(verts, alpha=0.2, fc='k', animated=True)
        self._ax.add_patch(poly)
        self._pi = BoundaryInteractor(self._ax, poly, beta=self._beta)
        # self._ax.add_line(self._pi.line)
    
    def _on_key(self, event):
        if event.key is 'enter':
            self._boundary_interactor()
    
    def _on_click(self, event):
        self._xdata.append(event.xdata)
        self._ydata.append(event.ydata)
        self._line.set_data(self._xdata, self._ydata)
        self._ax.figure.canvas.draw_idle()
    
    def __init__(self, xdata=[], ydata=[], beta=None, ax=None):
        
        if ax is None: ax = pl.gca()
        self._ax = ax
        
        self._xdata = list(xdata);
        self._ydata = list(ydata)
        self._beta = beta
        self._pi = None
        self._line = pl.Line2D(self._xdata,self._ydata,marker='o', markerfacecolor='b')
        self._ax.add_line(self._line)
        
        self._key_id = self._ax.figure.canvas.mpl_connect('key_press_event', self._on_key)
        self._click_id = self._ax.figure.canvas.mpl_connect('button_press_event', self._on_click)
        if len(xdata) > 0 and len(ydata) > 0:
            self._boundary_interactor()
    
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
    
    def get_beta(self):
        if self._pi is None:
            return self._beta
        else:
            return self._pi.beta 

    def get_verts(self):
        return zip(self.get_xdata(), self.get_ydata())

    x = property(get_xdata)
    y = property(get_ydata)
    beta = property(get_beta)
    verts = property(get_verts)



if __name__ == '__main__':
    p = BoundaryClick()
    pl.show()

