
from pylab import figure;
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy
from numpy import ndarray, ones, zeros, fabs, linspace
from math import pow, sqrt, fabs, acos, pi
import weld_stitch.weld._weld_geometry as wg

def make_image_coordinates(eb, fixed_coordinate_value=0.0, plane='x-z', haz=1.0):

    # param fixed_coordinate_value: scalar on axis normal to input plane
    #    if 'x-y'==plane, then fixed_coordinate_value should be a 'z' coordinate
    #    if 'x-z'==plane, then fixed_coordinate_value should be a 'y' coordinate
    #    if 'y-z'==plane, then fixed_coordinate_value should be a 'x' coordinate

    num_cells_short_axis=10
    if num_cells_short_axis <= 0 or int is not type(num_cells_short_axis):
        raise TypeError("'num_cells_short_axis' must be an 'integer' >=0")

    t=eb.plate_thickness
    a=eb.pool_width/2.0
    b=eb.pool_length/2.0
    x_range=(-(a+2*haz),(a+2*haz))
    y_range=(-(b+2*haz),(b+2*haz))
    z_range=(0,-(t))

    extent=[]
    c=fixed_coordinate_value
    if 'x-y'==plane:
        c=fixed_coordinate_value
        extent.append(x_range)
        extent.append(y_range)
        extent.append((c,c))
    elif 'x-z'==plane:
        extent.append(x_range)
        extent.append((c,c))
        extent.append(z_range)
    elif 'y-z'==plane:
        extent.append((c,c))
        extent.append(y_range)
        extent.append(z_range)
    else:
        raise TypeError("'plane' must be a string with one of the values: 'x-y', 'x-z', 'y-z'")

    # Compute non-zero axis sizes
    size=[]
    for ex in extent:
        s=fabs(ex[1]-ex[0])
        if s>0.0: size.append(s)

    # find mininum axis size that is greater than 0.0
    lmin=min(size)

    #: h grid spacing
    h=lmin/num_cells_short_axis

    axes=[]
    n=[]
    for i,ex in enumerate(extent):
        # size of axis
        s=fabs(ex[1]-ex[0])
        num_cells_axis=int(s/h)
        if 0==num_cells_axis: 
            n.append(1)
            ar=numpy.ndarray(shape=(1,),dtype=numpy.float64)
            ar.fill(ex[0])
            axes.append(ar)

        else:
            n.append(num_cells_axis)
            axes.append(linspace(ex[0],ex[1],num_cells_axis))

    nx=n[0];ny=n[1];nz=n[2];
    if num_cells_short_axis==nx*ny*nz:
        raise ValueError("More than one axis has zero extent")

    # Used for '2d' image plot
    # dim_2d
    # extent_2d
    dim_2d=None
    extent_2d=None
    if 1==nx:
        dim_2d=(ny,nz)
        extent_2d=(extent[1][0],extent[1][1],extent[2][0],extent[2][1])
    elif 1==ny:
        dim_2d=(nx,nz)
        extent_2d=(extent[0][0],extent[0][1],extent[2][0],extent[2][1])
    else:
        dim_2d=(nx,ny)
        extent_2d=(extent[0][0],extent[0][1],extent[1][0],extent[1][1])

    num_points=nx*ny*nz
    xyz=numpy.ndarray(shape=(num_points,3),dtype=numpy.float64)
    ax=axes[0];ay=axes[1];az=axes[2];
    for k in range(nz):
        z=az[k]
        for j in range(ny):
            y=ay[j]
            for i in range(nx):
                x=ax[i]
                p=k*nx*ny+j*nx+i
                xyz[p,0]=x;
                xyz[p,1]=y;
                xyz[p,2]=z;
    return dim_2d, extent_2d, xyz

def make_distance_image(eb,haz=50.0,plane='x-y',fixed_coordinate_value=0):
    dim_2d,extent_2d,xyz=make_image_coordinates(eb,fixed_coordinate_value,plane,haz)
    num_points=xyz.shape[0]
    d=ndarray(shape=(num_points,),dtype=numpy.float64)
    for p in range(num_points):
        d[p]=eb.distance(xyz[p,:])

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.imshow(d.reshape(dim_2d[1],dim_2d[0]),aspect='equal',extent=extent_2d)
    # turn axes 'on' or 'off'
    # maybe more for 'white' space around image
    # see http://stackoverflow.com/questions/9295026/matplotlib-plots-removing-axis-legends-and-white-spaces
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.show()
    return fig

def make_contour_image(eb,haz=50.0,plane='x-y',fixed_coordinate_value=0):
    dim_2d,extent_2d,xyz=make_image_coordinates(eb,fixed_coordinate_value,plane,haz)
    num_points=xyz.shape[0]
    d=ndarray(shape=(num_points,),dtype=numpy.float64)
    for p in range(num_points):
        d[p]=eb.distance(xyz[p,:])

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.imshow(d.reshape(dim_2d[1],dim_2d[0]),origin='upper',aspect='equal',extent=extent_2d)
    levels=numpy.arange(0,60,10)
    CS=ax.contour(d.reshape(dim_2d[1],dim_2d[0]),levels,colors='k',origin='upper',linewidths=2,aspect='equal',extent=extent_2d)

    # turn axes 'on' or 'off'
    # maybe more for 'white' space around image
    # see http://stackoverflow.com/questions/9295026/matplotlib-plots-removing-axis-legends-and-white-spaces
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.show()
    return fig

class EllipseBezier(object):

    def __init__(self,width=100.0,length=150.0,thickness=35.0,alpha=.5,beta=0.75):
        """
        Creates 'ellipse' tensor product with Bezier curve
        """
        self._pool_width=width
        self._pool_length=length
        self._plate_thickness=thickness
        self._alpha=alpha
        self._beta=beta
        d={}
        d['pool_width']=width
        d['pool_length']=length
        d['plate_thickness']=thickness
        d['scale_param']=alpha
        d['interpolate_param']=beta
        self._capsule=wg._get_ellipse(**d)

    def distance(self,xyz):
        if not type(xyz) is numpy.ndarray or (3,)!=xyz.shape:
            raise TypeError("Expected ndarray with shape=(3,), dtype=numpy.float64")
        return wg._distance(self._capsule,xyz)

    def get_pool_width(self):
        return self._pool_width

    def get_pool_length(self):
        return self._pool_length

    def get_plate_thickness(self):
        return self._plate_thickness

    pool_width=property(get_pool_width)
    pool_length=property(get_pool_length)
    plate_thickness=property(get_plate_thickness)

class Ellipse(object):

    def __init__(self,a=1.0,b=2.0):
        """
        Creates 'ellipse' with sizes (a,b) along x-axis and y-axis respectively
        """
        self._a=a
        self._b=b

    def closest_point(self,x,y):

        a=self.a
        b=self.b
        a2=pow(a,2)
        b2=pow(b,2)
        ax=a*x;
        by=b*y;

        def f(t):
            c1=ax/(a2+t)
            c2=by/(b2+t)
            return pow(c1,2)+pow(c2,2)-1.0

        def df(t):
            c1=-2.0*pow(ax,2)/pow(a2+t,3)
            c2=-2.0*pow(by,2)/pow(b2+t,3)
            return c1+c2

        def coordinates(t):
            x0=a2*x/(a2+t)
            y0=b2*y/(b2+t)
            return x0,y0

        def distance(x0,y0):
            dx=fabs(x0-x);
            dy=fabs(y0-y);
            if 0<dx and dx>=dy: return dx*sqrt(1.0+pow(dy/dx,2))
            elif 0<dy and dy>=dx: return dy*sqrt(1.0+pow(dx/dy,2))
            return 0.0

        def invert(x0,y0):
            """Given points x0,y0 on surface of ellipse, invert 
            for theta parameter that describes ellipse as a curve in space"""
            if x0/a+1.0<0.0: 
                ac=M_PI
            elif x0/a-1.0>0.0:
                ac=0
            else: 
                ac=acos(x0/a)
            # handles quadrants I and II
            if y0>=0.0: return ac

            # handles quadrants III and IV
            return 2*pi-ac

        def solve():
            max_iter=20
            tolerance=1.0e-15
            iter=0
            print("Iteration #   Residual\n----------   ---------\n")
            s="{0:3d}       {1:12.8e}"
            t=0
            r=f(t)
            print(s.format(iter,r))
            while r>tolerance:
                t=t-r/df(t)
                r=f(t)
                iter+=1
                print(s.format(iter,r))
            return t

        t=solve()
        x0,y0=coordinates(t)
        d=distance(x0,y0)
        theta=invert(x0,y0)

        s="Closest point: {0:12.8e}, {1:12.8e}; distance = {2:12.8e}, angle theta (radians) = {3:12.8e} "
        print(s.format(x0,y0,d,theta))

        return x0,y0,d


    def get_a(self):
        return self._a

    def get_b(self):
        return self._b

    a=property(get_a)
    b=property(get_b)
