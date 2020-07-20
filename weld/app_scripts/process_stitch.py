import numpy
from mpi4py import MPI
from matplotlib import pyplot as plt
from stitch.libstitch import libstitch
from stitch import time_equivalence

# Stitch filename
stitch_filename="stitch_weld.st"
# open stitch file
rc,fid=libstitch.open(MPI.COMM_WORLD,stitch_filename)
# Get all times stored in file
rc,times=libstitch.get_times(fid)

# Times in stitch file
print("Times read=\n%s"%(times,))

absolute_tolerance=1.0e-9
relative_tolerance=1,0e-15
eq=time_equivalence.Equivalence(absolute_tolerance,relative_tolerance)

# These are parameters
bx=[0,127]
by0=[86,356]
bz=[0,30]
d=0.0

# Frequency of output
dt=2
vel=9
db=vel*dt

# Layer
layer=0

# Create a moving box
start=[86,356]
middle=[814,1076]
end=[1354,1616]
whole=[86,1616]
bb,prefix=bx+whole+bz,"whole_"
bb,prefix=bx+start+bz,"start_"
bb,prefix=bx+middle+bz,"middle_"
bb,prefix=bx+end+bz,"end_"
# start images
# ls start_@([1-9]).png start_@([1][0-9]).png start_@([2][0-5]).png

for n,t in enumerate(times):
    # Debugging hack
    #if t>35.0:
    #    break
    block=numpy.fromiter(bb,dtype=numpy.int32).reshape(2,3,order='F')
    (rc, state, u_state, read_flag) = libstitch.read_block (fid, t, block)
    print("time t=%7.1f, read_flag=%d, distance traveled=%7.1f, \n\tblock=\n%s"%(t,read_flag,d,block,))
    # Increment distance traveled
    #d+=db
    # Increment box
    #byn[0]+=db; byn[1]+=db
    #bb=bx+byn+bz
    # NOTE 'zfill VERY necessary so that files are naturally order
    #    on linux/bash command line
    outfile=prefix+str(n).zfill(2)+".png"
    #print("output image file name = %s"%(outfile,))
    fig=plt.figure()
    ax=fig.subplots(1)
    cmap=plt.get_cmap("gist_rainbow")
    ax.imshow(state[:,:,layer],interpolation=None,cmap=cmap)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.savefig(outfile,bbox_inches='tight')
    plt.close(fig)
