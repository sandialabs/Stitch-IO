# Copyright 2019 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525
# with NTESS, the U.S. Government retains certain rights in this software.
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
# 
# For more information, contact Jay Lofstead (gflofst@sandeia.gov) or
# John Mitchell (jamitch@sandia.gov) for more information.

from math import sqrt, fabs
import numpy
from numpy import linspace, ndarray
from matplotlib import pyplot as plt
import stitch.am_weld.potts._teardrop as td

def get_demo_control_points(case="20"):
    n=4
    p=numpy.zeros(shape=(n+1,2),dtype=numpy.float64)
    if "20"==case:
        a=[-6.35,0.00]
        b=[a[0],1.50]
        c=[-3.25,3.00]
        d=[0.0,2.75]
        e=[d[0],0.0]
    elif "25"==case:
        a=[-7.35,0.0]
        b=[a[0],1.0]
        c=[-3.75,2.6]
        d=[0.0,2.8]
        e=[d[0],0.0]
    elif "30"==case:
        a=[-7.35,0.0]
        b=[a[0],0.7]
        c=[-3.75,2.0]
        d=[0.0,2.475]
        e=[d[0],0.0]

    p[0,0]= a[0];  p[0,1]= a[1];
    p[1,0]= b[0];  p[1,1]= b[1];
    p[2,0]= c[0];  p[2,1]= c[1];
    p[3,0]= d[0];  p[3,1]= d[1];
    p[4,0]= e[0];  p[4,1]= e[1];
    return p

def make_image_coordinates(teardrop2d_capsule, plate_thickness=0.0, fixed_coordinate_value=0.0, plane='x-z', haz=1.0):

    # param fixed_coordinate_value: scalar on axis normal to input plane
    #    if 'x-y'==plane, then fixed_coordinate_value should be a 'z' coordinate
    #    if 'x-z'==plane, then fixed_coordinate_value should be a 'y' coordinate
    #    if 'y-z'==plane, then fixed_coordinate_value should be a 'x' coordinate

    num_cells_short_axis=200
    if num_cells_short_axis <= 0 or int is not type(num_cells_short_axis):
        raise TypeError("'num_cells_short_axis' must be an 'integer' >=0")

    t=plate_thickness
    pool_width,pool_length=td._get_teardrop2d_size(teardrop2d_capsule)
    a=pool_width/2.0
    b=pool_length/2.0
    x_range=(-(b+2*haz),(b+2*haz))
    y_range=(-(a+2*haz),(a+2*haz))
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

def make_distance_image(case="30",W=100.0,plane='x-y',with_contours=True,fixed_coordinate_value=0):
    # Retrieve pool and 'Teardrop2D' capsule
    capsule,td=get_unit_pool(case,W=W)
    
    W,L=td._get_teardrop2d_size(capsule)
    haz=W/2
    plate_thickness=0.0
    dim_2d,extent_2d,xyz=make_image_coordinates(capsule,plate_thickness,fixed_coordinate_value,plane,haz)
    num_points=xyz.shape[0]
    d=ndarray(shape=(num_points,),dtype=numpy.float64)
    s="CPP violation at x={0:10.4e}, y={1:10.4e}"
    for p in range(num_points):
        x=xyz[p,0]
        y=xyz[p,1]
        # Teardrop is symmetric about 'x-axis'; use absolute value of 'y'
        u,distance=td._get_teardrop2d_cpp(x,fabs(y),capsule);
        if 0>u or 1<u: print(s.format(x,y))
        d[p]=distance

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.imshow(d.reshape(dim_2d[1],dim_2d[0]),origin='lower',aspect='equal',interpolation='none',extent=extent_2d)
    if with_contours is True:
        ax.imshow(d.reshape(dim_2d[1],dim_2d[0]),origin='lower',aspect='equal',extent=extent_2d)
        levels=numpy.arange(0,60,10)
        levels=numpy.arange(0,2*haz,haz/2)
        CS=ax.contour(d.reshape(dim_2d[1],dim_2d[0]),levels,colors='k',origin='lower',linewidths=2,extent=extent_2d)
    # turn axes 'on' or 'off'
    # maybe more for 'white' space around image
    # see http://stackoverflow.com/questions/9295026/matplotlib-plots-removing-axis-legends-and-white-spaces
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.show()

def get_unit_pool(case="20",W=1.0):
    control_points=get_demo_control_points(case=case)
    capsule=td._get_teardrop2d(control_points,W)
    # Pool size
    W0,L=td._get_teardrop2d_size(capsule)
    s="Pool width W={0:6.4f}, pool length L={1:6.4f}, L/W={2:6.4f}"
    print(s.format(W0,L,L/W0))
    # Pool position
    xp=td._get_teardrop2d_position(capsule)
    s="Pool position xp={0:6.4f}"
    print(s.format(xp))
    return capsule,td

def plot_curve(case="20",pool_width=100.0):
    capsule,td=get_unit_pool(case,W=pool_width)
    label={"20":"I","25":"II","30":"III"}
    thickness={"20":3.5,"25":3,"30":2}
    # Draw pool
    num_points_curve=100
    U=linspace(0.0,1.0,num_points_curve)
    c=td._teardrop2d_compute_curve(capsule,U)
    fontsize=24
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_title("Teardrop case "+label[case],fontsize=fontsize)
    #ax.annotate("Case "+label[case],xycoords='figure fraction', xy=(0.5,.5), \
    #        textcoords='figure fraction', xytext=(0.5275,.65), horizontalalignment='center', fontsize=fontsize)
    ax.tick_params(axis='both', which='major',labelsize=fontsize)
    ax.set_xlabel("x-axis (sites)", fontsize=fontsize)
    ax.set_ylabel("y-axis (sites)", fontsize=fontsize)
    # Plot curve
    t=thickness[case]
    ax.plot(c[0,:],c[1,:],color='b',marker='None',linewidth=t)
    # SYMMETRY: Plot curve for y<0
    ax.plot(c[0,:],-c[1,:],color='b',marker='None',linewidth=t)
    ax.axis('equal')
    fig.subplots_adjust(left=0.2,bottom=.17,top=0.92)
    fig.show()

def demo_predictor(x,y,case="20",with_tangent_arrow=True):
    capsule,td=get_unit_pool(case)
    # Draw pool
    num_points_curve=100
    U=linspace(0.0,1.0,num_points_curve)
    c=td._teardrop2d_compute_curve(capsule,U)
    print("c.shape=",str(c.shape))
    fig=plt.figure()
    ax=fig.add_subplot(111)
    # Plot curve
    ax.plot(c[0,:],c[1,:],color='b',marker='None')
    # SYMMETRY: Plot curve for y<0
    ax.plot(c[0,:],-c[1,:],color='b',marker='None')
    # Plot input point
    ax.plot(x,y,color='r',linestyle='None',marker='x')

    # Define function that bounds cp
    def bound_cp(du,umin,umax):
        un=0.5*(umin+umax)
        while umax-umin > du:
            rho,rho_1=td._teardrop2d_tangent_to_curve(un,capsule)
            dr=numpy.fromiter([x-rho[0],y-rho[1]],dtype=numpy.float64)
            dt=numpy.dot(dr,rho_1)
            if dt>0:
                # Point is north of un
                umin=un
            else:
                # Point is south of un
                umax=un
            un=0.5*(umin+umax)
        return umin,umax

    # Compute y-elevation of maximum pool width
    um=td._teardrop2d_parametric_coordinate_at_maximum_pool_width(capsule)
    U=numpy.fromiter([um],dtype=numpy.float64)
    c=td._teardrop2d_compute_curve(capsule,U)
    xm=c[0,0]
    minlim=0.5
    maxlim=1.0

    # Identify quadrant
    if x>xm:
        umin=um;
        umax=1
    else:
        umin=0
        umax=um

    # Compute bounding points on curve
    du=.025
    umin,umax=bound_cp(du,umin,umax)
    U=numpy.fromiter([umin,umax],dtype=numpy.float64)
    c=td._teardrop2d_compute_curve(capsule,U)

    # Plot bounding points on curve
    ax.plot(c[0,:],c[1,:],color='g',marker='o')

    # Compute cpp onto line spanning two points
    delta=numpy.fromiter([c[0,1]-c[0,0],c[1,1]-c[1,0]],dtype=numpy.float64)
    dr=numpy.fromiter([x-c[0,0],y-c[1,0]],dtype=numpy.float64)
    t=numpy.dot(delta,dr)/numpy.dot(delta,delta)
    s=(1-t)*umin+t*umax
    pstr="Computed t={0:4.2f}, umin={1:4.2f}, umax={2:4.2f}, s={3:4.2f}"
    print(pstr.format(t,umin,umax,s))
    if 0 > t or t > 1: 
        print("OVERRIDE t,umin={0:4.2f}, minlim={1:4.2f}, umax={2:4.2f}, maxlim={3:4.2f}".format(umin,minlim,umax,maxlim))
        if t>1 and umax<maxlim: 
            t=1
        else: 
            if t<0 and umin>minlim: t=0
            else: t=0.5
        s=(1-t)*umin+t*umax
        print(pstr.format(t,umin,umax,s))


    # Plot closet point prediction
    rho,rho_1=td._teardrop2d_tangent_to_curve(s,capsule)
    arrow_props={"fc":"k","ec":"k","head_width":0.005,"head_length":0.005}
    # Estimated normal
    ax.arrow(x,y,-x+rho[0],-y+rho[1],**arrow_props)
    n=numpy.fromiter([rho_1[1],-rho_1[0]],dtype=numpy.float64)
    dr=numpy.fromiter([x-rho[0],y-rho[1]],dtype=numpy.float64)
    dot=numpy.dot(n,dr)
    if dot>=0:
        print("INTERIOR point given.")
    else:
        print("EXTERIOR point given.")

    if True==with_tangent_arrow:
        # Tangent at 'rho'
        m=sqrt(rho_1[0]*rho_1[0]+rho_1[1]*rho_1[1])
        ax.arrow(rho[0],rho[1],rho_1[0]/m,rho_1[1]/m,**arrow_props)

    ax.plot(x,y,color='r',marker='x')

    ax.axis('equal')
    fig.show()
    return s

def teardrop2d_cpp(x,y,u0,case="20"):
    """
    To run this function, first run above function:
    u0=demo_predictor(x,y,case="20",with_tangent_arrow=True)

    This provides a good initial guess value for parametric coordinate 
    of closest point projection; then use the predicted value 'u0'
    to finalize the cpp using Newton iterations.
    teardrop_cpp(x,y,u0)
    """
    capsule,td=get_unit_pool(case)
    u=u0; 
    j11,rho,dr,rho_1,rho_11=td._get_teardrop2d_iterate(x,y,u,capsule);
    r1=rho_1.dot(dr)
    resid=fabs(r1)
    iter_max=10
    tolerance=1.0e-15
    iter=0
    header= "Iter #   U        r1         j11       dr[0]       dr[1]   \n"
    header+="---------------------------------------------------------"
    s="{0:4d}    {1:5.3f},  {2:10.4e}   {3:7.4f}   {4:7.4f}    {5:7.4f} "
    print(header)
    while iter<=iter_max and resid>tolerance:
        print(s.format(iter,u,r1,j11,dr[0],dr[1]))
        u+=-r1/j11
        j11,rho,dr,rho_1,rho_11=td._get_teardrop2d_iterate(x,y,u,capsule);
        r1=rho_1.dot(dr)
        resid=fabs(r1)
        iter+=1

def teardrop2d_parametric_coordinate_at_maximum_pool_width(case="20"):
    """
    """
    capsule,td=get_unit_pool(case)
    u0=0.25
    u=u0; 
    # x,y are arbitrary and do not play a role in this calculation;
    # Only used here because of interface to capsule has them. Only 
    # need rho_1 and rho_11
    x=0.0;y=0.0;
    j11,rho,dr,rho_1,rho_11=td._get_teardrop2d_iterate(x,y,u,capsule);
    f=rho_1[0]
    j=rho_11[0]
    resid=fabs(f)
    iter_max=10
    tolerance=1.0e-15
    iter=0
    header= "Iter #   U        resid  \n"
    header+="---------------------------"
    s="{0:4d}    {1:5.3f},  {2:10.4e}   "
    print(header)
    while iter<=iter_max and resid>tolerance:
        print(s.format(iter,u,resid))
        u+=-f/j
        j11,rho,dr,rho_1,rho_11=td._get_teardrop2d_iterate(x,y,u,capsule);
        f=rho_1[0]
        j=rho_11[0]
        resid=fabs(f)
        iter+=1

    umax=td._teardrop2d_parametric_coordinate_at_maximum_pool_width(capsule)
    s="Parametric coordinate u={0:5.3f}"
    print(s.format(umax))

def plot_cpp_functional(x,y,umin,umax,case="20"):
    capsule,td=get_unit_pool(case)
    # Tweak umin and umax to look at interval of interest for input x and y
    # Nice example:  plot_cpp_functional(0.0,.5,.1,.4,case="20")
    a=umin
    b=umax
    num_points_curve=100
    U=linspace(a,b,num_points_curve)
    fc=numpy.zeros(shape=U.shape,dtype=numpy.float64)
    for i,u in enumerate(U):
        j11,rho,dr,rho_1,rho_11=td._get_teardrop2d_iterate(x,y,u,capsule);
        fc[i]=-rho_1.dot(dr)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(U,fc)
    fig.show()


def make_capsule_demo_plot(case="20",pool_width=1.0):
    capsule,td=get_unit_pool(case,W=pool_width)
    num_points_curve=100
    U=linspace(0.0,1.0,num_points_curve)
    c=td._teardrop2d_compute_curve(capsule,U)
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_aspect('equal')
    ax.plot(c[0,:],c[1,:])
    fig.show()
    return ax,fig


def make_demo_plot(case="20"):
    """
    In constrast with above functionality, this function 
    uses python defined BernsteinPolynomials below to 
    plot pool shape.
    """
    control_points=get_demo_control_points(case=case)
    num_points_curve=100
    polynomial_order=4;
    bp=BernsteinPolynomials(polynomial_order)
    c=bp.compute_curve(num_points_curve,control_points)
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(c[0,:],c[1,:])
    ax.axis('equal')
    fig.show()
    return ax,fig


class BernsteinPolynomials(object):
    def __init__(self,n):
        self._n=n
        self._b=numpy.ndarray(shape=(n+1,),dtype=numpy.float64)

    def eval(self,u):
        n=self._n
        b=self._b
        b[0]=1.0
        u1=1-u;
        for j in range(1,n+1):
            bj=0.0
            for k in range(0,j):
                t=b[k]
                b[k]=bj+u1*t
                bj=u*t;
            b[j]=bj;
        return b

    def compute_curve(self,num_points_curve,control_points):
        
        # Parametric values range[0,1] used to compute curve for 
        #   input control_points
        U=linspace(0.0,1.0,num_points_curve)
        n=self._n
        # Control points must have shape=(n+1,2)

        # Curve
        # 1st row is x-coordinates
        # 2nd row is y-coordinates
        c=numpy.ndarray(shape=(2,num_points_curve),dtype=numpy.float64)

        for i,u in enumerate(U):

            # Evaluate array of Bernstein polynomials for parametric coordinate 'u'
            b=self.eval(u)
            # Initialize curve value at 'u'
            c[0,i]=c[1,i]=0.0
            
            for j in range(n+1):
                x=control_points[j,0]
                y=control_points[j,1]
                c[0,i]+=x*b[j]
                c[1,i]+=y*b[j]

        return c

