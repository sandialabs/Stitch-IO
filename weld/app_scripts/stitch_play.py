import getopt
import os, sys, json
import numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.animation as animation
from weld_stitch.weld import bezier_teardrop as bt
from weld_stitch.app_scripts import haz_box as hb

dtype=hb.dtype

def usage():
    print("""

    Name: 

    Description: Reads json file for stitch weld parameters; Writes two files related to computational 
                 volume box location and size; 1) potts model to initialize microstructures; 2) weld model 
                 for simulation on computational volume.


    Example usage:
    Options:
        -h            Print this message.
        -i            json file containing model parameters for computational volume
    """)

def make_rectangles(haz_box,t0,pool_position,alpha,pool_speed,start_up=False):
    tnp1,yp_np1,potts_init_cv,cv=haz_box.get_cv_boxes(t0,pool_position,alpha,pool_speed,start_up)

    # Create rectangles
    # potts_init
    xlo=potts_init_cv[0,0]; xhi=potts_init_cv[0,1];
    ylo=potts_init_cv[1,0]; yhi=potts_init_cv[1,1];
    h=(xhi-xlo); w=(yhi-ylo);
    pr=Rectangle((ylo,xlo),w,h);
    # weld
    xlo=cv[0,0]; xhi=cv[0,1];
    ylo=cv[1,0]; yhi=cv[1,1];
    h=(xhi-xlo); w=(yhi-ylo);
    wr=Rectangle((ylo,xlo),w,h);
    return tnp1,yp_np1,pr,wr

class StitchWeldOptions(object):
    def __init__(self, opts, args):
		
        # check for options
        jason_params_file=None
        yp=None
        t0=None
        start_up=False
        s=""
        for o, a in opts:
            if o=="-i":
                json_params_file=a
                w="json input file specified = {0:s}\n".format(json_params_file)
                s+=w
            elif o=="--yp":
                yp=float(a)
                w="Input weld pool position = {0:3.1f}\n".format(yp)
                s+=w
            elif o=="--t0":
                t0=float(a)
                w="SPPARKS simulation start time = {0:3.1f}\n".format(t0)
                s+=w
            elif o=="--startup":
                w="Startup flag set\n"
                start_up=True
                s+=w
            else:
                print("ERROR: invalid option(s)")
                usage()
                sys.exit(2)

        self._json_params_file=json_params_file
        self._pool_position=yp
        self._t0=t0
        self._start_up=start_up
        self._options = s
        pass

    def print_options(self):
        print(self._options)

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:",longopts=["yp=","t0=","startup"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    # Get inputs
    swo = StitchWeldOptions(opts,args)
    return swo


def create_patch_list(potts_rectangles,color='white'):
    patches=[]
    for r in potts_rectangles:
       w=r.get_width();h=r.get_height();xy=r.get_xy()
       p=Rectangle(xy,w,h,facecolor=color,alpha=1.0,zorder=0)
       patches.append(p)
    return patches

def print_patch_list(potts_rectangles,weld_rectangles):
    s="POTTS: width=%d, xlo=%4d, xhi=%4d, ylo=%4d, yhi=%d"
    t="WELD : width=%d, xlo=%4d, xhi=%4d, ylo=%4d, yhi=%d"
    for u,v in zip(potts_rectangles,weld_rectangles):
       w=u.get_width();h=u.get_height();xy=u.get_xy()
       print(s%(w,xy[0],xy[0]+w,xy[1],xy[1]+h))
       w=v.get_width();h=v.get_height();xy=v.get_xy()
       print(t%(w,xy[0],xy[0]+w,xy[1],xy[1]+h))
    pass


# python stitch_play.py -istitch_weld.in --yp=0.0 --t0=0.0
if __name__=='__main__':
    # Process command line
    swo=main()

    # Open parameters file; load json
    inpfile=open(swo._json_params_file,'r')
    json_params=json.load(inpfile)
    inpfile.close()

    # Read main parameters section from json
    params=json_params["parameters"]
    case=params["case"]
    pool_width=params["pool_width"]
    plate_thickness=params["plate_thickness"]
    haz=params["haz"]
    haz_cushion=params["haz_cushion"]
    haz_box=hb.HazBox(case,pool_width,plate_thickness,haz,haz_cushion)

    # Initial starting time
    t0=swo._t0
    # Initial pool position
    pool_position=swo._pool_position

    # total weld distance
    weld_distance=params["weld_distance"]

    alpha=params["alpha"]
    weld_speed=params["weld_speed"]

    float_str="variable {0:s} equal {1:14.1f}\n"
    start_up=True
    yn=pool_position
    tn=t0
    distance_traveled=0.0
    potts_rectangles=[]
    weld_rectangles=[]
    stop=0

    # Store sequence of pool positions for latter use
    box=haz_box.get_haz_box(yn); dx=(box[0,1]-box[0,0])/2; 
    yp=[]
    yp.append(yn)

    # Compute pool curve
    capsule,td=bt.get_unit_pool(case,W=pool_width)
    num_points_curve=100
    U=numpy.linspace(0.0,1.0,num_points_curve)
    curve=td._teardrop2d_compute_curve(capsule,U)

    # NOTE
    # yp list has one more entry than potts_rectangles because it is
    #   initialized above with initial pool position; See above yp.append(yn)
    while distance_traveled<weld_distance:
        tnp1,ynp1,pr,wr=make_rectangles(haz_box,tn,yn,alpha,weld_speed,start_up=start_up)
        pr.set_facecolor('black')
        potts_rectangles.append(pr)
        weld_rectangles.append(wr)
        start_up=False
        distance_traveled+=(ynp1-yn)
        tn=tnp1
        yn=ynp1
        yp.append(yn)

    #print_patch_list(potts_rectangles,weld_rectangles)

    # Create sequence of images
    num_rect=len(potts_rectangles)

    for s in range(0,num_rect):

        # Must create new array of patches each time 
        patches=create_patch_list(potts_rectangles,color='white')
        for r in range(0,s):
            patches[r].set_facecolor('black')

        # Set this stage colors
        c=plt.cm.gnuplot2(.8)
        e=patches[s]
        e.set_hatch('/')
        e.set_facecolor(c)

        ## Create a figure and add patches
        fig,ax=plt.subplots(1)
        for r in range(0,s+1):
            ax.add_patch(patches[r])

        # Plot pool curve on image
        # Curve needs to be rotated for consistency with 'stitch_play' 
        #    Stitch play x-axis is 'y-axis' in the teardrop
        # Also need to translate pool along axis according to its position.
        # Also need to translate pool vertically along y-axis to center in 
        #    stitch_play images
        _y=-curve[0,:]+dx
        y=curve[0,:]+dx; x=curve[1,:]+yp[s]
        ax.plot(x,y,color='b')
        ax.plot(x,_y,color='b')
        ax.set_xlim(0,1600)
        ax.set_aspect('equal')
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        # Create file name with leading 0's
        frame_num_str='0'*(6-len(str(2*s)))+str(2*s)
        frame='frame_'+frame_num_str+'.png'
        fig.savefig(frame)
        plt.close(fig)

        # Second image for stage advances pool over hatched
        #   area above.  Otherwise, this image is identical
        #   to above.
        # Must create new array of patches each time 
        patches=create_patch_list(weld_rectangles)
        for r in range(0,s+1):
            patches[r].set_facecolor('black')

        ## Create a figure and add patches
        fig,ax=plt.subplots(1)
        for r in range(0,s+1):
            ax.add_patch(patches[r])


        # Plot pool curve on image
        # Curve needs to be rotated for consistency with 'stitch_play' 
        #    Stitch play x-axis is 'y-axis' in the teardrop
        # Also need to translate pool along axis according to its position.
        # Also need to translate pool vertically along y-axis to center in 
        #    stitch_play images
        _y=-curve[0,:]+dx
        y=curve[0,:]+dx; x=curve[1,:]+yp[s+1]
        ax.plot(x,y,color='b')
        ax.plot(x,_y,color='b')
        ax.set_xlim(0,1600)
        ax.set_aspect('equal')
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        frame_num_str='0'*(6-len(str(2*s+1)))+str(2*s+1)
        frame='frame_'+frame_num_str+'.png'
        fig.savefig(frame)
        plt.close(fig)
