import getopt
import os, sys, json
import numpy
from weld_stitch.weld import bezier_teardrop as bt

# Use 32Bit integers
dtype=numpy.int32

def convert_bb_to_libstitch_block(bb):
    xlo=bb[0,0]; xhi=bb[0,1];
    ylo=bb[1,0]; yhi=bb[1,1];
    zlo=bb[2,0]; zhi=bb[2,1];
    if not numpy.dtype(dtype) is bb.dtype:
        raise TypeError("Input block must be ndarray dtype="+str(dtype))
    box=numpy.fromiter([xlo,xhi,ylo,yhi,zlo,zhi],dtype=bb.dtype).reshape(2,3,order='F')
    return box

def compute_num_sites(bb):
    """
    input bb: type=array(shape=(3,2), dtype=numpy.int32)
              represents bounding box for input to spparks for use 
              in the 'region' command
    output Q: type=integer representing number of spparks sites 
    """
    if not numpy.dtype(dtype) is bb.dtype:
        raise TypeError("Input block must be ndarray dtype="+str(dtype))
    Q=1
    for i in range(3):
        Q*=bb[i,1]-bb[i,0]
    return Q

def get_bb_box_to_spparks_parameters_str(bb):
    """
    Write spparks 'vars' to file; for use in parameterizing spparks runs.

    input fileop: type=open('filename','w'); previously opened file
                  writing.
    input bb: type=array(shape=(3,2), dtype=numpy.int32)
              represents bounding box for input to spparks for use 
              in the 'region' command.
    output: x0,x1,y0,y1,z0,z1 parameter definitions used for 
            var substitution in spparks run.
    output Q: type=integer representing number of spparks sites.
    """
    if not numpy.dtype(dtype) is bb.dtype:
        raise TypeError("Input block must be ndarray dtype="+str(dtype))
    # Number of SPPARKS sites in simulation
    Q=compute_num_sites(bb)
    x_str="variable {0:s} equal {1:14d}\n"
    bb_str=x_str.format('Q',Q)
    bb_str+=x_str.format('X0',bb[0,0])
    bb_str+=x_str.format('X1',bb[0,1])
    bb_str+=x_str.format('Y0',bb[1,0])
    bb_str+=x_str.format('Y1',bb[1,1])
    bb_str+=x_str.format('Z0',bb[2,0])
    bb_str+=x_str.format('Z1',bb[2,1])
    return bb_str

def pretty_string_bb_box(bb):
    """
    Use this to print string representation of bb.

    input bb: type=array(shape=(3,2), dtype=numpy.int32)
              represents bounding box for input to spparks for use 
              in the 'region' command.

    Following strings are concatenated for output.
    output: string='x0,x1,y0,y1,z0,z1 = %d, %d, %d, %d, %d, %d\n'
    output: string = 'Number of sites' %d
    """
    if not numpy.dtype(dtype) is bb.dtype:
        raise TypeError("Input block must be ndarray dtype="+str(dtype))
    Q=compute_num_sites(bb)
    x_str="Number of sites = {0:14d}\n"
    bb_str=x_str.format(Q)
    x_str="x0,x1,y0,y1,z0,z1={0:d},{1:d},{2:d},{3:d},{4:d},{5:d}\n"
    bb_str+=x_str.format(bb[0,0],bb[0,1],bb[1,0],bb[1,1],bb[2,0],bb[2,1])
    return bb_str

class HazBox(object):

    def __init__(self,case="20",pool_width=102.7,plate_thickness=30,haz=10, haz_cushion=2):

        # Work with 32Bit integers
        self._state_dtype=dtype

        # These lengths are in 'lattice' units used to build simulation box in spparks.
        # Pool width is generally a floating point number but using 
        #   length scale units of the lattice.
        # Plate thickness must be an integer and is used to describe 
        #       the exact number of lattice sites required for plate thickness.
        # The lattice length scale parameter 'c' is computed
        #       c=(actual plate thickness measured in microns)/ plate_thickness
        # Then actual/physical pool width can be calculated if multiplied by 'c'.
        self._plate_thickness=plate_thickness
        self._haz=haz
        self._haz_cushion=haz_cushion

        # Computed pool_length is a floating point number
        # pool_width is a floating point number
        self._capsule,td=bt.get_unit_pool(case="20",W=pool_width)
        self._pool_width,self._pool_length=td._get_teardrop2d_size(self._capsule)

        # Pool position is geometric average along axis of weld
        # see 'get_pool_position' in teardrop.h

    def get_haz_box(self,pool_position):
        """
        input pool_position: type=floating point position of pool along y-axis
        output bounding box: type=numpy.array(shape=(3,),dtype=numpy.int32); lower left corner, upper right corner
                         associated with heat affected zone.
        """
        # Leading and trailing edges of pool along y-axis
        # pool moves along y-axis
        ymax=pool_position+self._pool_length/2
        ymin=pool_position-self._pool_length/2
        # Build a bounding box for HAZ around pool
        bb=numpy.zeros(shape=(3,2),dtype=self._state_dtype)
        # Lower left corner of box
        bb[0,0]=0
        bb[1,0]=numpy.floor(ymin)-(self._haz+self._haz_cushion)
        bb[2,0]=0;# bb bottom is always at z=0

        # Upper right corner of box
        bb[0,1]=numpy.ceil(self._pool_width)+2*(self._haz+self._haz_cushion)
        bb[1,1]=numpy.ceil(ymax)+(self._haz+self._haz_cushion)
        bb[2,1]=self._plate_thickness # top of box is plate thickness

        return bb


    def write_weld_parameters(self,fileop,start_time,vp,pool_position,cv_bb,delta_t,T):
        """
        param bb: type=numpy array, shape=3,2, spparks box: lower left, upper right
        param T: type=float, spparks simulation temperature
        """
        # T: SPPARKS simulation temperature
        float_str="variable {0:s} equal {1:4.2f}\n"
        o_str =float_str.format('T',T)
        # T0: start time for this CV
        # DT: simulation duration for this CV
        float_str="variable {0:s} equal {1:14.1f}\n"
        o_str+=float_str.format('T0',start_time)
        o_str+=float_str.format('DT',delta_t)
        o_str+=float_str.format('VEL',vp)
        # YP: starting pool position
        # WP: pool width
        # HAZ: heat affected zone
        o_str+=float_str.format('YP',pool_position)
        o_str+=float_str.format('WP',self._pool_width)
        x_str="variable {0:s} equal {1:14d}\n"
        o_str+=x_str.format('HAZ',self._haz)
        fileop.write(o_str)
        write_bb_box_to_spparks_parameters_file(fileop,cv_bb)
        pass

    def write_potts_parameters(self,fileop,bb,T):
        """
        param bb: type=numpy array, shape=3,2, spparks box: lower left, upper right
        param T: type=float, spparks simulation temperature
        """
        # T: SPPARKS simulation temperature
        float_str="variable {0:s} equal {1:4.2f}\n"
        o_str =float_str.format('T',T)
        # Write box 
        fileop.write(o_str)
        write_bb_box_to_spparks_parameters_file(fileop,bb)
        pass


    def compute_alpha(self,pool_position,alpha_0,pool_speed,dS):
        # SEE Peridynamics AM notebook May 15 2018
        # alpha_0 is an input parameter;
        # It is used to scale up the size of the moving haz.  By scaling the
        # HAZ, the distance traveled by the pool is shorter or longer depending
        # upon the magnitude of alpha.
        # However, alpha is constrained to produce travel distances such
        # that for the input pool_speed the required time of travel dt
        # is a multiple 'n' of the 'spparks_output_frequency dS'
        # dt = n * dS

        # Pool speed alias
        v=pool_speed

        # Step n pool position
        yp_n=pool_position
        haz_n=self.get_haz_box(yp_n)

        # Length of moving haz
        L0=haz_n[1,1]-haz_n[1,0]
        if(L0<=0.0):
            raise RuntimeError("Length of HAZ evaluated to < 0.0 which is a runtime error. Check inputs.")
        n0=int((alpha_0*L0)/(v*dS))
        i=0
        alpha_i=0.0
        while alpha_i<alpha_0:
            i+=1
            alpha_i=i*dS*v/L0
            dt=i*dS
            d=alpha_i*L0
            # DEBUG print
            # print("i=%d,alpha_i=%3.1f,d=alpha_i*L0=%5.1f,d/v=dt=%5.2f;dt=i*ds=%5.2f"%(i,alpha_i,d,d/v,i*dS))
            # END DEBUG print
        alpha=alpha_i
        return alpha

    def get_cv_boxes(self,t0,dS,pool_position,alpha_0,pool_speed,start_up=False):

        # Initial time for this CV
        tn=t0

        # Alias pool speed
        vp=pool_speed

        # Step n pool position
        yp_n=pool_position
        haz_n=self.get_haz_box(yp_n)

        # Length of moving haz
        hl=haz_n[1,1]-haz_n[1,0]


        # first cv is slightly different
        if start_up:
            # Length of CV controlled by user-specified 'alpha' multiplier
            # However, alpha is constrained
            alpha=self.compute_alpha(pool_position,alpha_0,pool_speed,dS)
            cvl=alpha*hl

            # First CV travels a bit further than subsequent CV
            distance=cvl

            # Duration of time spent in one CV
            dt=distance/vp

            # Step 'n+1' pool position
            yp_np1=yp_n+distance
            haz_np1=self.get_haz_box(yp_np1)

            # Initial CV consists of upper right corner of haz_np1 
            #   but with lower left corner of haz_n for x,z components and upper right
            #   of haz_n for y-component
            cv=numpy.copy(haz_n)
            cv[:,0]=haz_n[:,0]; cv[1,0]=haz_n[1,1];
            cv[:,1]=haz_np1[:,1]

            # Bounding box for 'potts_init'
            potts_init_cv=cv

        # all subsequent cv
        else:
            # Length of CV controlled by user-specified 'alpha' multiplier
            # However, alpha is constrained; Also it is shorter for these cases.
            alpha=self.compute_alpha(pool_position,alpha_0-1.0,pool_speed,dS)
            # Distance traveled within this CV
            distance=alpha*hl

            # Duration of time spent in one CV
            dt=distance/vp

            # Next weld pool location
            yp_np1=yp_n+distance

            # Final haz location for this CV
            haz_np1=self.get_haz_box(yp_np1)

            # CV consists of upper right corner of haz_np1 and lower left corner of haz_np
            cv=numpy.copy(haz_n)
            cv[:,1]=haz_np1[:,1]

            # CV for 'potts_init' consists of upper right corner of haz_np1 
            #   but with lower left corner of haz_n for x,z components and upper right
            #   of haz_n for y-component
            potts_init_cv=numpy.copy(haz_np1)
            potts_init_cv[:,0]=haz_n[:,0]; potts_init_cv[1,0]=haz_n[1,1];

        # Update starting time for this CV
        tnp1=tn+dt
        return tnp1,yp_np1,potts_init_cv,cv

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

def write_spparks_stage(json_params_file,t0,pool_position,start_up):
    # Open parameters file; load json
    inpfile=open(json_params_file,'r')
    json_params=json.load(inpfile)
    inpfile.close()

    # Read main output section from json
    output=json_params["output"]
    prefix=output['prefix']
    dS=output['spparks_dump_dt']
    dump_abs_tol=output["spparks_dump_abs_tol"]

    # Open log file
    #logfilename=json_params["log"]
    #sys.stdout=open(logfilename,'a')

    # Read main parameters section from json
    params=json_params["parameters"]
    case=params["case"]
    pool_width=params["pool_width"]
    plate_thickness=params["plate_thickness"]
    haz=params["haz"]
    haz_cushion=params["haz_cushion"]

    # DEBUG print
    #p_str="Case=%s, pool_width=%5.1f, plate_thickness=%3d, HAZ=%3d, haz_cushion=%2d"
    #print(p_str%(case,pool_width,plate_thickness,haz,haz_cushion))
    # END DEBUG print
    haz_box=HazBox(case,pool_width,plate_thickness,haz,haz_cushion)

    alpha=params["alpha"]
    if alpha<=1.0:
        e_str="Invalid scaling parameter alpha; alpha<=1.0 is invalid; input alpha={0:15.5e}".format(alpha)
        raise ValueError(e_str)
    weld_speed=params["weld_speed"]
    tnp1,ynp1,potts_init_cv,cv=haz_box.get_cv_boxes(t0,dS,pool_position,alpha,weld_speed,start_up=start_up)

    def make_weld_parameters_output_str():
        # T: SPPARKS simulation temperature
        T=params["simulation_temperature"]
        float_str="variable {0:s} equal {1:4.2f}\n"
        o_str =float_str.format('T',T)
        # T0: start time for this CV
        # DT: simulation duration for this CV
        float_str="variable {0:s} equal {1:14.1f}\n"
        #o_str+=float_str.format('DT',tnp1-t0)
        o_str+=float_str.format('DS',dS)
        o_str+=float_str.format('T0_plus_EPS',t0+.5)
        o_str+=float_str.format('TNP1',tnp1)
        float_str="variable {0:s} equal {1:5.3f}\n"
        o_str+=float_str.format('ABS_TOL',dump_abs_tol)
        # VEL: simulation weld speed
        o_str+=float_str.format('VEL',weld_speed)
        # YP: starting pool position
        # WP: pool width
        # HAZ: heat affected zone
        o_str+=float_str.format('YP',pool_position)
        o_str+=float_str.format('WP',pool_width)
        x_str="variable {0:s} equal {1:14d}\n"
        o_str+=x_str.format('HAZ',haz)
        o_str+=get_bb_box_to_spparks_parameters_str(cv)
        # Name of parameters file
        params_file=prefix+"_weld.params.txt"
        # Write string to parameters file
        wf=open(params_file,'w')
        wf.write(o_str)
        wf.close()

    def make_potts_parameters_output_str():
        # T: SPPARKS simulation temperature
        T=params["simulation_temperature"]
        float_str="variable {0:s} equal {1:4.2f}\n"
        o_str =float_str.format('T',T)
        o_str+=get_bb_box_to_spparks_parameters_str(potts_init_cv)
        # Name of parameters file
        params_file=prefix+"_potts.params.txt"
        # Write string to parameters file
        wf=open(params_file,'w')
        wf.write(o_str)
        wf.close()

    # Print parameters files
    make_weld_parameters_output_str()
    make_potts_parameters_output_str()

    return tnp1,ynp1


# python haz_box.py -istitch_weld.in --yp=0.0 --t0=<float value>
if __name__=='__main__':
    # Capture print statements to a log file
    logfilename="stitch.log"
    sys.stdout=open(logfilename,'a')

    # Process command line
    swo=main()
    swo.print_options()

    json_params_file=swo._json_params_file
    t0=swo._t0
    pool_position=swo._pool_position
    start_up=swo._start_up

    tnp1,ynp1=write_spparks_stage(json_params_file,t0,pool_position,start_up)

    # Switch print to standard out; this allows bash script to
    #  capture tnp1, ynp1
    sys.stdout=sys.__stdout__
    print(tnp1,ynp1)
