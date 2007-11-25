'Interactive boundary generation for use with gridgen.'

# Created by Rob Hetland on 2007-01-21.
# Copyright (c) 2007 Texas A&M Univsersity. All rights reserved.

# Version 1.0.0 on 2007-11-25
#   Removed second BoundaryInteractor class, and added interative gridgen
#   capabilities.

from numpy import *
import pylab as pl
import pyroms
import copy
import cPickle

from matplotlib.artist import Artist
from matplotlib.patches import Polygon, CirclePolygon
from matplotlib.lines import Line2D
from matplotlib.numerix import sqrt, nonzero, equal, asarray, dot, Float
from matplotlib.numerix.mlab import amin
from matplotlib.mlab import dist_point_to_segment


class BoundaryClick(object):
    """
    bry = BoundaryClick(x=[], y=[], beta=None, ax=gca(), **gridgen_options)
    
    If x, y and beta are not given, and interactive polygon creation session is
    started. The initial boundary is sketched out by clicking the points (in a
    counterclockwise manner, usually starting in the upper lefthand corner of the
    boundary). At this point the verticies are marked by orange circles.
    
    Switch to editing mode by hitting return. This changes the vertecies to black. At
    this point, beta values may be created to define the corners of the grid (see
    below). Current data are always available in bry.x, bry.y (combined in p.verts)
    and bry.beta attributes.
    
    Key commands:
        
        enter : switch to grid editing mode
        
        t : toggle visibility of verticies
        d : delete a vertex
        i : insert a vertex at a point on the polygon line
        
        p : define vertex as beta=1 (a Positive turn, marked with green triangle)
        m : define vertex as beta=1 (a Negative turn, marked with red triangle)
        z : define vertex as beta=0 (no corner, marked with green triangle)
        
        G : generate grid from the current boundary using gridgen
        T : toggle visability of the current grid
        R : remove the gridlines from the figure
    
    Methods:
    
        bry.write_bry(bry_file)
            Write the current boundary informtion (bry.x, bry.y, bry.beta) to a cPickle
            file bry_file.
        
        bry.read_bry(bry_file)
            Read in boundary informtion (x, y, beta) from the cPickle file bry_file.
        
        bry.remove_grid()  
            Remove gridlines from axes.
        
    """
    
    _showverts = True
    _showbetas = True
    _showgrid = True
    _epsilon = 5  # max pixel distance to count as a vertex hit
    
    def _update_beta_lines(self):
        """Update the m/pline by finding the (x,y) points where self.beta== -/+ 1"""
        x = self._line.get_xdata()
        y = self._line.get_ydata()
        
        xp = [x[n] for n in range(len(x)) if self.beta[n]==1]
        yp = [y[n] for n in range(len(y)) if self.beta[n]==1]
        self._pline.set_data(xp, yp)
        
        xm = [x[n] for n in range(len(x)) if self.beta[n]==-1]
        ym = [y[n] for n in range(len(y)) if self.beta[n]==-1]
        self._mline.set_data(xm, ym)
        
        xz = [x[n] for n in range(len(x)) if self.beta[n]==0]
        yz = [y[n] for n in range(len(y)) if self.beta[n]==0]
        self._zline.set_data(xz, yz)
        
        if len(x)-1 < self.gridgen_options['ul_idx']:
            self.gridgen_options['ul_idx'] = len(x)-1
        xs = x[self.gridgen_options['ul_idx']]
        ys = y[self.gridgen_options['ul_idx']]
        self._sline.set_data(xs, ys)
    
    def remove_grid(self):
        """Remove a generated grid from the BoundaryClick figure"""
        if hasattr(self, '_gridlines'):
            for line in self._gridlines:
                line.remove()
            delattr(self, '_gridlines')
    
    def _init_boundary_interactor(self):
        """Send polyclick thread to boundary interactor"""
        # Get rid old mpl connections.
        self._ax.figure.canvas.mpl_disconnect(self._key_id)
        self._ax.figure.canvas.mpl_disconnect(self._click_id)
        
        # Shade the selected region in a polygon
        self._poly = Polygon(self.verts, alpha=0.2, fc='k', animated=True)
        self._ax.add_patch(self._poly)
        
        # change line to animated
        # self._line.set_animated(True)
        # self._line.set_markerfacecolor('k')
        self._line.set_marker('None')
        
        # Link in the two lines that will show the beta values
        # pline for positive turns, mline for negative (minus) turns.
        self._pline = Line2D([], [], marker='^', ms=12, mfc='g', animated=True, lw=0)
        self._mline = Line2D([], [], marker='v', ms=12, mfc='r', animated=True, lw=0)
        self._zline = Line2D([], [], marker='o', mfc='k', animated=True, lw=0)
        self._sline = Line2D([], [], marker='s', mfc='k', animated=True, lw=0)
        
        self._update_beta_lines()
        self._ax.add_artist(self._pline)
        self._ax.add_artist(self._mline)
        self._ax.add_artist(self._zline)
        self._ax.add_artist(self._sline)
        
        # get the canvas and connect the callback events
        cid = self._poly.add_callback(self._poly_changed)
        self._ind = None # the active vert
        
        self._canvas.mpl_connect('draw_event', self._draw_callback)
        self._canvas.mpl_connect('button_press_event', self._button_press_callback)
        self._canvas.mpl_connect('key_press_event', self._key_press_callback)
        self._canvas.mpl_connect('button_release_event', self._button_release_callback)
        self._canvas.mpl_connect('motion_notify_event', self._motion_notify_callback)
    
    def _on_return(self, event):
        if event.key is 'enter':
            self._init_boundary_interactor()
    
    def _on_click(self, event):
        self.x.append(event.xdata)
        self.y.append(event.ydata)
        self.beta.append(0)
        self._line.set_data(self.x, self.y)
        
        self._background = self._canvas.copy_from_bbox(self._ax.bbox)
        self._ax.draw_artist(self._line)
        self._canvas.blit(self._ax.bbox)
    
    
    def _draw_callback(self, event):
        self._background = self._canvas.copy_from_bbox(self._ax.bbox)
        self._ax.draw_artist(self._poly)
        self._ax.draw_artist(self._pline)
        self._ax.draw_artist(self._mline)
        self._ax.draw_artist(self._zline)
        self._ax.draw_artist(self._sline)
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
        # if event.key=='t':
        #     self._showverts = not self._showverts
        #     self._line.set_visible(self._showverts)
        #     if not self._showverts: self._ind = None
        if event.key=='t':
            self._showbetas = not self._showbetas
            self._pline.set_visible(self._showbetas)
            self._mline.set_visible(self._showbetas)
            self._zline.set_visible(self._showbetas)
            self._sline.set_visible(self._showbetas)
        elif event.key=='d':
            ind = self._get_ind_under_point(event)
            if ind is not None:
                self._poly.xy = [tup for i,tup in enumerate(self._poly.xy) if i!=ind]
                self._line.set_data(zip(*self._poly.xy))
        elif event.key=='p':
            ind = self._get_ind_under_point(event)
            if ind is not None:
                self.beta[ind] = 1.0
        elif event.key=='m':
            ind = self._get_ind_under_point(event)
            if ind is not None:
                self.beta[ind] = -1.0
        elif event.key=='z':
            ind = self._get_ind_under_point(event)
            if ind is not None:
                self.beta[ind] = 0.0
        elif event.key=='s':
            ind = self._get_ind_under_point(event)
            if ind is not None:
                self.gridgen_options['ul_idx'] = ind
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
                    self.beta.insert(i+1, 0)
                    break
            s0 = xys[-1]
            s1 = xys[0]
            d = dist_point_to_segment(p, s0, s1)
            if d<=self._epsilon:
                self._poly.xy.append((event.xdata, event.ydata))
                self._line.set_data(zip(*self._poly.xy))
                self.beta.append(0)
        elif event.key=='G':
            options = copy.deepcopy(self.gridgen_options)
            shp = options.pop('shp')
            self.grd = pyroms.gridgen(self.x, self.y, self.beta, shp, **options)
            self.remove_grid()
            self._showgrid = True
            # axlim = self._ax.axis()
            gridlineprops = {'linestyle':'-', 'color':'k', 'lw':0.25, 'alpha':0.5}
            self._gridlines = []
            for line in self._ax._get_lines(*(self.grd.x, self.grd.y), **gridlineprops):
                self._ax.add_line(line)
                self._gridlines.append(line)
            for line in self._ax._get_lines(*(self.grd.x.T, self.grd.y.T), **gridlineprops):
                self._ax.add_line(line)
                self._gridlines.append(line)
            # self._ax.axis(axlim)
        elif event.key=='D':
            self.remove_grid()
        elif event.key=='T':
            self._showgrid = not self._showgrid
            if hasattr(self, '_lines'):
                for line in self._lines:
                    line.set_visible(self._showgrid)
            if hasattr(self, '_linesT'):
                for line in self._linesT:
                    line.set_visible(self._showgrid)
        
        self._update_beta_lines()
        self._draw_callback(event)
        self._canvas.draw()
    
    def _motion_notify_callback(self, event):
        'on mouse movement'
        # if not self._showverts: return
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        x,y = event.xdata, event.ydata
        self._poly.xy[self._ind] = x,y
        
        x, y = zip(*self._poly.xy)
        self._line.set_data(x, y)
        self._update_beta_lines()
        
        self._canvas.restore_region(self._background)
        self._ax.draw_artist(self._poly)
        self._ax.draw_artist(self._pline)
        self._ax.draw_artist(self._mline)
        self._ax.draw_artist(self._zline)
        self._ax.draw_artist(self._sline)
        self._ax.draw_artist(self._line)
        self._canvas.blit(self._ax.bbox)
    
    
    def __init__(self, x=[], y=[], beta=None, ax=None, **gridgen_options):
        
        if isinstance(x, str):
            x, y, beta = load(x)
        
        if ax is None: ax = pl.gca()
        self._ax = ax

        # Set default gridgen option, and copy over specified options.
        self.gridgen_options = {'ul_idx': 0, 
                                'shp': (32, 32),
                                'verbose': True}
        for key, value in gridgen_options.iteritems():
            self.gridgen_options[key] = gridgen_options[key]
        
        x = list(x); y = list(y)
        assert len(x)==len(y), 'arrays must be equal length'
        
        if beta is None:
            self.beta = [0 for xi in x]
        else:
            assert len(x)==len(beta), 'beta must have same length as x and y'
            self.beta = list(beta)
        
        self._line = pl.Line2D(x, y, marker='o', markerfacecolor='orange', animated=True)
        self._ax.add_line(self._line)
        
        self._canvas = self._line.figure.canvas        
        
        self._key_id = self._ax.figure.canvas.mpl_connect('key_press_event', self._on_return)
        self._click_id = self._ax.figure.canvas.mpl_connect('button_press_event', self._on_click)
        if len(x) > 0:
            self._init_boundary_interactor()
    
    def write_bry(self, bry_file):
        f = open(bry_file, 'wb')
        cPickle.dump((self.x, self.y, self.beta), f)
        f.close()
    
    def read_bry(self, bry_file):
        x, y, self.beta = load(bry_file)
        self._line.set_data(x, y)
        if hasattr(self, '_poly'):
            self._poly.xy = zip(x, y)
            self._update_beta_lines()
            self._draw_callback(None)
            self._canvas.draw()
    
    def _get_verts(self):return zip(self.x, self.y)
    verts = property(_get_verts)    
    def get_xdata(self): return self._line.get_xdata()
    x = property(get_xdata)
    def get_ydata(self): return self._line.get_ydata()
    y = property(get_ydata)
    


if __name__ == '__main__':
    from matplotlib.toolkits.basemap import Basemap
    from matplotlib.ticker import FuncFormatter
    
    m = Basemap(projection='lcc', 
                resolution='i',
                llcrnrlon=5.0,
                llcrnrlat= 52.,
                urcrnrlon=35.5,
                urcrnrlat=68.0,
                lat_0=15.00, 
                lon_0=0.0,
                suppress_ticks=False)
    
    fig=pl.figure()
    # background color will be used for 'wet' areas.
    fig.add_axes([0.1,0.1,0.8,0.8],axisbg='azure')
    m.drawcoastlines()
    m.fillcontinents(color='tan')
    
    x, y, beta = zip(*[(241782.65384467551, 1019981.3539730886,   1.0),
                       (263546.12512432877,  686274.79435173667,  0),
                       (452162.87621465814,  207478.42619936226,  1.0),
                       (1431519.0837990609,  432367.62942244718,  1.0),
                       (1409755.6125194072,  816855.62202965701, -1.0),
                       (1612881.344462839,   816855.62202965701,  1.0),
                       (1620135.83488939,   1194089.1242103155,  -1.0),
                       (2048150.7700559068, 1302906.4806085825,   1.0),
                       (1982860.3562169466, 1520541.1934051164,   1.0),
                       (1358974.1795335496, 1447996.2891396051,  -1.0),
                       (1344465.1986804469, 1716412.4349219969,   0),
                       (1620135.83488939,   2318535.1403257404,   1.0),
                       (1199375.3901494243, 2383825.5541647007,   1.0),
                       (967231.6964997882,  1825229.7913202639,   0),
                       (916450.26351393061, 1266634.0284758268,   0),
                       (851159.84967497038,  715292.75605794112, -1.0),
                       (691561.06029084534,  693529.28477828787, -1.0),
                       (510198.79962706729, 1128798.7103713555,   1.0)])
    
    p = BoundaryClick(x=x, y=y, beta=beta)
    pl.show()

