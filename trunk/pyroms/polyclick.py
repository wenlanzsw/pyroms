'Interactive polygon generation.'

# Copyright R. Hetland on 2007-11-25.  
# All rights reserved.

# Version 1.0.0 on 2007-11-25
#    Initial (new) version based on version 1.0.0 of BoundaryClick

from numpy import *
import pylab as pl
import cPickle
from pyroms import Polygeom

from matplotlib.artist import Artist
from matplotlib.patches import Polygon, CirclePolygon
from matplotlib.lines import Line2D
from matplotlib.numerix import sqrt, nonzero, equal, asarray, dot, Float
from matplotlib.numerix.mlab import amin
from matplotlib.mlab import dist_point_to_segment


class PolyClick(object):
    """
    bry = PolyClick(verts, ax=gca())
    
    If verts are not given, an interactive polygon creation session is started.
    Switch to editing mode by hitting return. This changes the vertecies to black. The
    points may be moved with the mouse, or modified with the key commands listed below.
    Current data are always available in bry.x, bry.y and bry.verts.
    
    Key commands:
        
        enter : switch to grid editing mode
        
        t : toggle visibility of verticies
        d : delete a vertex
        i : insert a vertex at a point on the polygon line
        
    Methods:
    
        bry.write_bry(bry_file='bry.pickle)
            Write the current boundary informtion (bry.verts) to a cPickle
            file bry_file.
        
        bry.read_bry(bry_file='bry.pickle)
            Read in boundary informtion (verts) from the cPickle file bry_file.
        
    """
    
    _showverts = True
    _epsilon = 5  # max pixel distance to count as a vertex hit
    
    def _init_boundary_interactor(self):
        """Send polyclick thread to boundary interactor"""
        # Get rid old mpl connections.
        
        self._ax.figure.canvas.mpl_disconnect(self._key_id)
        self._ax.figure.canvas.mpl_disconnect(self._click_id)
        
        # Shade the selected region in a polygon
        self._poly = Polygon(self.verts, alpha=0.2, fc='k', animated=True)
        self._ax.add_patch(self._poly)
        
        # modify the line properties
        self._line.set_markerfacecolor('k')
        
        # get the canvas and connect the callback events
        cid = self._poly.add_callback(self._poly_changed)
        self._ind = None # the active vert
        
        self._canvas.mpl_connect('draw_event', self._draw_callback)
        self._canvas.mpl_connect('button_press_event', self._button_press_callback)
        self._canvas.mpl_connect('key_press_event', self._key_press_callback)
        self._canvas.mpl_connect('button_release_event', self._button_release_callback)
        self._canvas.mpl_connect('motion_notify_event', self._motion_notify_callback)
    
    
    def _draw_callback(self, event):
        self._background = self._canvas.copy_from_bbox(self._ax.bbox)
        self._ax.draw_artist(self._poly)
        self._ax.draw_artist(self._line)
        self._canvas.blit(self._ax.bbox)
    
    def _poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self._line.get_visible()
        Artist.update_from(self._line, poly)
        self._line.set_visible(vis)  # don't use the poly visibility state
    
    def _get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'
        x, y = zip(*self._poly.xy)
        
        # display coords
        xt, yt = self._poly.get_transform().numerix_x_y(x, y)
        d = sqrt((xt-event.x)**2 + (yt-event.y)**2)
        indseq = nonzero(equal(d, amin(d)))
        ind = indseq[0]
        
        if d[ind]>=self._epsilon:
            ind = None
        
        return ind
    
    def _button_press_callback(self, event):
        'whenever a mouse button is pressed'
        # if not self._showverts: return
        if event.inaxes==None: return
        if event.button != 1: return
        self._ind = self._get_ind_under_point(event)
    
    def _button_release_callback(self, event):
        'whenever a mouse button is released'
        # if not self._showverts: return
        if event.button != 1: return
        self._ind = None
    
    def _key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes: return
        if event.key=='t':
            self._showverts = not self._showverts
        elif event.key=='d':
            ind = self._get_ind_under_point(event)
            if ind is not None:
                self._poly.xy = [tup for i,tup in enumerate(self._poly.xy) if i!=ind]
                self._line.set_data(zip(*self._poly.xy))
        elif event.key=='i':
            xys = self._poly.get_transform().seq_xy_tups(self._poly.xy)
            p = event.x, event.y # display coords
            for i in range(len(xys)-1):
                s0 = xys[i]
                s1 = xys[i+1]
                d = dist_point_to_segment(p, s0, s1)
                if d<=self._epsilon:
                    self._poly.xy.insert(i+1, (event.xdata, event.ydata))
                    self._line.set_data(zip(*self._poly.xy))
                    break
            s0 = xys[-1]
            s1 = xys[0]
            d = dist_point_to_segment(p, s0, s1)
            if d<=self._epsilon:
                self._poly.xy.append((event.xdata, event.ydata))
                self._line.set_data(zip(*self._poly.xy))
        
        self._draw_callback(event)
        self._canvas.draw()
    
    def _motion_notify_callback(self, event):
        'on mouse movement'
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        x,y = event.xdata, event.ydata
        self._poly.xy[self._ind] = x,y
        
        x, y = zip(*self._poly.xy)
        self._line.set_data(x, y)
        
        self._canvas.restore_region(self._background)
        self._ax.draw_artist(self._poly)
        self._ax.draw_artist(self._line)
        self._canvas.blit(self._ax.bbox)
    
    def _on_return(self, event):
        if event.key is 'enter':
            self._init_boundary_interactor()
    
    def _on_click(self, event):
        self.x.append(event.xdata)
        self.y.append(event.ydata)
        self._line.set_data(self.x, self.y)
        
        self._background = self._canvas.copy_from_bbox(self._ax.bbox)
        self._ax.draw_artist(self._line)
        self._canvas.blit(self._ax.bbox)
    
    def __init__(self, *args, **kwargs):
        
        ax = kwargs.pop('ax', None)
        if ax is None: ax = pl.gca()
        self._ax = ax
        
        assert len(args) <=2, 'Input is either verticies (1 arg), or x and y (2 args)'
        
        if args is not ():
            if isinstance(args[0], str):
                verts = load(x)
                x, y = zip(*verts)
                x = list(x); y = list(y)
            elif len(args) == 1:  # assume verticies input
                verts = args[0]
                x, y = zip(*verts)
                x = list(x); y = list(y)
            else:               # x and y inputs
                x = list(args[0])
                y = list(args[1])
                assert len(x) == len(y), 'x and y must have the same length.'
        else: # no input
            x = []; y = []
        
        self._line = pl.Line2D(x, y, marker='o', markerfacecolor='orange', animated=True)
        self._ax.add_line(self._line)
        
        self._canvas = self._line.figure.canvas        
        
        self._key_id = self._ax.figure.canvas.mpl_connect('key_press_event', self._on_return)
        self._click_id = self._ax.figure.canvas.mpl_connect('button_press_event', self._on_click)
        if len(x) > 0:
            self._init_boundary_interactor()
    
    def write_poly(self, poly_file):
        f = open(bry_file, 'wb')
        cPickle.dump(self.verts, f)
        f.close()
    
    def read_poly(self, poly_file):
        verts = load(bry_file)
        x, y = zip(*verts)
        self._line.set_data(x, y)
        if hasattr(self, '_poly'):
            self._poly.xy = verts
            self._draw_callback(None)
            self._canvas.draw()
    
    def inside(self, *args):
        """docstring for inside"""
        assert len(args) > 0, 'must provide either verticies or x and y.'
        if len(args) == 1:  # assume verticies input
            verts = args[0]
        else:               # x and y inputs
            x = list(args[0])
            y = list(args[1])
            assert len(x) == len(y), 'x and y must have the same length.'
            verts = zip(x, y)
        return self.geom.inside(verts)
        
    
    def get_verts(self):return zip(self.x, self.y)
    verts = property(get_verts)    
    
    def get_xdata(self): return self._line.get_xdata()
    x = property(get_xdata)
    
    def get_ydata(self): return self._line.get_ydata()
    y = property(get_ydata)
    
    def _get_geom(self): return Polygeom(self.verts)
    geom = property(_get_geom)
    
    def _get_area(self): return self.geom.area
    area = property(_get_area)
    
    def _get_centroid(self): return self.geom.centroid
    centroid = property(_get_centroid)


if __name__ == '__main__':
    x, y = random.rand(2, 100)
    pl.plot(x, y, '.')
    
    p = PolyClick()
    pl.show()
    
    

