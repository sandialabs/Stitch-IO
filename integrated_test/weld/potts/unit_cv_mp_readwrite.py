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

import unittest
import numpy
import sys
import argparse
from math import fabs,pi,cos,sin,log
from numpy import array, ndarray, dot, transpose
from  mpi4py import MPI
from stitch.libstitch import libstitch

def pretty_string_block(block):
    """
    Use this to print string representation of block.

    input bb: type=array(shape=(2,3), dtype=numpy.int32)
              represents bounding box for input to spparks for use
              in the 'region' command.

    Following strings are concatenated for output.
    output: string = 'Number of sites' %d
    output: string='x0,x1,y0,y1,z0,z1 = %d, %d, %d, %d, %d, %d\n'
    """
    x_str="x0,x1,y0,y1,z0,z1={:d},{:d},{:d},{:d},{:d},{:d}\n"
    t=get_block_tuple(block)
    block_str=x_str.format(*t)
    return block_str

def get_block_tuple(block):
    """
    input block: type numpy.ndarray, shape=(2,3), order='Fortran'
    output x0,x1,y0,y1,z0,z1: type=tuple

    Example
    -------
    input:
        box=numpy.fromiter([0,127, 86,344,0,30],dtype=numpy.int32).reshape(2,3,order='F')
    output:
        (0,127,86,344,0,30)
    """
    shape=block.shape
    if 2!=shape[0] or 3!=shape[1]:
        raise ValueError("Input 'block' must be ndarray shape=2,3")
    if not numpy.dtype(numpy.int32) is block.dtype:
        raise TypeError("Input block must be ndarray dtype=numpy.32")

    t=(block[0,0],block[1,0],block[0,1],block[1,1],block[0,2],block[1,2])
    return t

def get_block_width(block):
    """
    input block: type numpy.ndarray, shape=(2,3), order='Fortran'
    output nx,ny,nz: type tuple, length=3, values=nx,ny,nz

    Example input:
        block=numpy.fromiter([0,127, 86,344,0,30],dtype=numpy.int32).reshape(2,3,order='F')
    """
    shape=block.shape
    if 2!=shape[0] or 3!=shape[1]:
        raise ValueError("Input 'block' must be ndarray shape=2,3")
    if not numpy.dtype(numpy.int32) is block.dtype and not numpy.dtype(numpy.int64) is block.dtype and not numpy.dtype(numpy.float64) is block.dtype:
        raise TypeError("Input block must be ndarray dtype=numpy.32")

    nx=block[1,0]-block[0,0];
    ny=block[1,1]-block[0,1];
    nz=block[1,2]-block[0,2];
    return nx,ny,nz

def get_Q(box):
    nx,ny,nz=get_block_width(box)
    return nx*ny*nz


class unit_fixed_size_cv_write_read(unittest.TestCase):
    NPX = 1
    NPY = 1
    NPZ = 1

    def setUp(self):
        unittest.TestCase.setUp(self)
        name='unit_fixed_size_cv_write_read'
        # 'suffix' 'st' for 'stitch file'
        self._fname=name+'.st'
        self._comm=MPI.COMM_WORLD

        self._state_dtype=numpy.int32
        # x1, x2, y1, y2, z1, z2
        b1=numpy.fromiter([0,127, 86,344,0,30],dtype=self._state_dtype).reshape(2,3,order='F')
        self._b1=b1
        nx,ny,nz=get_block_width(b1)
        self._nx=nx;
        self._ny=ny;
        self._nz=nz;
        self._Q=nx*ny*nz
        (rc,self._file) = libstitch.open (self._comm, self._fname);
	(rc, field_id) = libstitch.create_field (self._file, 'spin', 1, 1, -1)
        pass

    def write(self, timestamp, block, state):
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        libstitch.write_block (self._file, field_id, timestamp, block, state);

    def read(self, timestamp, block, state):
        # Read data
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, state, is_new_time) = libstitch.read_block (self._file, field_id, timestamp, block)
        return state

    def test_write_read(self):
        self.assertTrue(self._Q==982980)
        nx,ny,nz=get_block_width(self._b1)
        state=numpy.ndarray(shape=(nx,ny,nz),dtype=self._state_dtype,order='F')
        t0=0.0
        # Assign state
        for k in range(nz):
            # Assign values according to 'z' elevation
            v=k;
            state[:,:,k]=v

        b1_str=pretty_string_block(self._b1)
        self.write(t0, self._b1, state)

        # read
        trial_state=self.read(t0, self._b1, state)

        # Assert that trial was correctly read
        self.assertTrue((trial_state == state).all())
        pass

    def tearDown(self):
        # close the file
        libstitch.close (self._file);
        pass

class unit_variable_size_cv_write_read(unittest.TestCase):
    NPX = 1
    NPY = 1
    NPZ = 1
    ND = 4  # amount of data per process in each dimension. Increment is by half of this

    def new_timestamp(self, time):
        # make a new timestamp for a given time
        (rc, timestamp) = libstitch.create_timestamp (self._file, time)
        return timestamp;

    def setUp(self):
        unittest.TestCase.setUp(self)
        name='unit_variable_size_cv_write_read'

        # MPI comm, rank, and size
        self._comm=MPI.COMM_SELF
        self._pcomm=MPI.COMM_WORLD
        self._rank = self._pcomm.Get_rank()
        self._size = self._pcomm.Get_size()
        # 'suffix' 'st' for 'stitch file'
        self._fname=name+str(self._rank)+'.st'
        self._pfname='p' + name + '.st'
        # number of processes in each dimension
        self._npx = self.NPX
        self._npy = self.NPY
        self._npz = self.NPZ
        # amount of data per process in each dimension
        self._nd = self.ND
        # what the position for this process is in x,y,z
        self._posx = self._rank % self._npx
        self._posy = (self._rank // self._npx) % self._npy
        self._posz = self._rank // (self._npx * self._npy)
        #print ('rank ', self._rank, ' posx ', self._posx, ' posy ', self._posy, ' posz ', self._posz)
        # the offset from (0,0,0) for this process
        self._offx = self._posx * self._nd
        self._offy = self._posy * self._nd
        self._offz = self._posz * self._nd
        #print ('rank ', self._rank, ' offx ', self._offx, ' offy ', self._offy, ' offz ', self._offz)

        # open a single file for typical tests
        (rc,self._file) = libstitch.open (self._comm, self._fname);
	(rc,self._field_id) = libstitch.create_field (self._pfile, 'spin', 1, 1, -1)

        # open a single file for parallel specific tests
        (rc,self._pfile) = libstitch.open (self._pcomm, self._pfname);
	(rc,self._pfield_id) = libstitch.create_field (self._pfile, 'spin', 1, 1, -1)

        self._state_dtype=numpy.int32
        # Define 3 distinct blocks and then a union block
        # x1, x2, y1, y2, z1, z2
        # Blocks are 'non-overlapping'
        self._b1=numpy.fromiter([0,2,0,2,0,2],dtype=self._state_dtype).reshape(2,3,order='F')
        self._b2=numpy.fromiter([2,4,0,2,0,2],dtype=self._state_dtype).reshape(2,3,order='F')
        self._b3=numpy.fromiter([4,6,0,2,0,2],dtype=self._state_dtype).reshape(2,3,order='F')
        # This is a union of all 3 blocks b1 b2 b3
        self._bu=numpy.fromiter([0,6,0,2,0,2],dtype=self._state_dtype).reshape(2,3,order='F')

        # Define state which will be written to each block
        self._s1=  numpy.ones(shape=get_block_width(self._b1),dtype=self._state_dtype,order='F')
        self._s2=2*numpy.ones(shape=get_block_width(self._b2),dtype=self._state_dtype,order='F')
        self._s3=3*numpy.ones(shape=get_block_width(self._b3),dtype=self._state_dtype,order='F')

        self._state_dtype_i32=numpy.int32
        self._state_dtype_i64=numpy.int64
        self._state_dtype_f64=numpy.float64
        # multi-field test area
        self._b4=numpy.fromiter([6,8,0,2,0,2],dtype=self._state_dtype_i32).reshape(2,3,order='F')
        # multi-field test blocks
        self._s4=4*numpy.ones(shape=get_block_width(self._b4),dtype=self._state_dtype_i32,order='F')
        self._s5=5*numpy.ones(shape=get_block_width(self._b4),dtype=self._state_dtype_i64,order='F')
        self._s6=6.0*numpy.ones(shape=get_block_width(self._b4),dtype=self._state_dtype_f64,order='F')

        # retrieve or define a set of times; 
        #(rc,times) = libstitch.get_times (self._file) #  don't do
        self._t0=0.0
        self._t1=1.0
        self._t2=2.0
        self._t3=3.0

        # setup for the parallel tests
        # for the first iteration
        self._pb1 = numpy.fromiter ([self._offx, self._offx + self._nd, self._offy, self._offy + self._nd, self._offz, self._offz + self._nd],dtype=self._state_dtype).reshape(2,3,order='F')
        #print ('rank: ', self._rank, ' ', self._pb1)

        # for second iteration, shift over in x by per proc size + 50%
        x_increment = self._nd // 2 + (self._nd * (self._npx - 1))
        #print (x_increment)
        self._pb2 = numpy.fromiter ([self._offx + x_increment, self._offx + self._nd + x_increment, self._offy, self._offy + self._nd, self._offz, self._offz + self._nd],dtype=self._state_dtype).reshape(2,3,order='F')
        #print ('rank: ', self._rank, ' ', self._pb2)

        # this is the whole space. This is used only for a single proc for reading and verifying
        self._pbu = numpy.fromiter ([0, (self._npx - 1) * self._nd + x_increment, 0, self._npy * self._nd, 0, self._npz * self._nd],dtype=self._state_dtype).reshape (2, 3, order = 'F')
        #print ('rank: ', self._rank, ' ', self._pbu)
        pass

    def tearDown(self):
        # close the file
        libstitch.close (self._file);
        libstitch.close (self._pfile);
        pass


    # NOTE: for these tests to work correctly, they must be run in order of occurrence 
    #   as they are in this file; otherwise the latter test could fail.
    # One way to apparently achieve this is to add alphabetical orders to 
    #   at beginning of each test name.
    def test_a_write_read_b1_at_t1(self):
        # Write/read b1 and 'assert'
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        libstitch.write_block (self._file, field_id, self._t1, self._b1, self._s1);
        (rc, trial_s1, is_new_time) = libstitch.read_block (self._file, field_id, self._t1, self._b1)
        self.assertTrue((trial_s1 == self._s1).all())
        pass

    def test_b_write_read_b1_at_t1(self):
        # Write/read b1 and 'assert'
        # libstitch.write_block (self._file, self._t0, self._b1, self._s1);
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, trial_s1, is_new_time) = libstitch.read_block (self._file, field_id, self._t1, self._b1)
        self.assertTrue((trial_s1 == self._s1).all())
        pass

    def test_c_read_on_block_b1_at_nonexistent_prior_time(self):
        # This time 't0' is prior to time 't0' previously written;
        #    So read should get all '-1' values
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, trial_m1, is_new_time) = libstitch.read_block (self._file, field_id, self._t0, self._b1)
        m1=-1*numpy.ones(shape=get_block_width(self._b1),dtype=self._state_dtype,order='F')
        self.assertTrue((trial_m1 == m1).all())
        pass

    def test_d_write_read_b2_at_t2(self):
        # Write/read b2 and 'assert'
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        libstitch.write_block (self._file, field_id, self._t2, self._b2, self._s2);
        (rc, trial_s2, is_new_time) = libstitch.read_block (self._file, field_id, self._t2, self._b2)
        self.assertTrue((trial_s2 == self._s2).all())
        pass

    def test_e_write_read_b3_at_t3(self):
        # Write/read b3 and 'assert'
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        libstitch.write_block (self._file, field_id, self._t3, self._b3, self._s3);
        (rc, trial_s3, is_new_time) = libstitch.read_block (self._file, field_id, self._t3, self._b3)
        self.assertTrue((trial_s3 == self._s3).all())
        pass

    def test_f_read_write_on_union(self):
        # Create su on block union 'bu'
        # It is already correct for s1 @ t3
        # Assign for s2 @ t3
        # Assign for s3 @ t3
        su=numpy.ones(shape=get_block_width(self._bu),dtype=self._state_dtype,order='F')
        su[2:4,:,:]=2
        su[4:6,:,:]=3
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, trial_su, is_new_time) = libstitch.read_block (self._file, field_id, self._t3, self._bu)
        # FAILS
        # print ('su')
        # print (su)
        # print ('trial_su')
        # print (trial_su)
        self.assertTrue((trial_su == su).all())
        # Now write union at later time
        t4 = 4.0
        su[:,:,:]=4
        libstitch.write_block (self._file, field_id, t4, self._bu, su);
        #print ('su ', t4)
        ##print (su)
        # Assert 
        (rc, trial_s4, is_new_time) = libstitch.read_block (self._file, field_id, t4, self._bu)
        #print ('trial_s4 ', t4)
        #print (trial_s4)
        self.assertTrue((trial_s4 == su).all())
        pass

    def test_g_read_blocks_at_later_time(self):
        # re-run existing test
        self.test_c_read_on_block_b1_at_nonexistent_prior_time()

        # Re-read original blocks at t3 and assert values
        # b1
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, trial_s1, is_new_time) = libstitch.read_block (self._file, field_id, self._t3, self._b1)
        self.assertTrue((trial_s1 == self._s1).all())
        # b2
        (rc, trial_s2, is_new_time) = libstitch.read_block (self._file, field_id, self._t3, self._b2)
        self.assertTrue((trial_s2 == self._s2).all())
        # b3
        (rc, trial_s3, is_new_time) = libstitch.read_block (self._file, field_id, self._t3, self._b3)
        self.assertTrue((trial_s3 == self._s3).all())
        pass

    def test_h_parallel_write_read (self):
        # write from each process
        su = numpy.ones (shape = get_block_width (self._pb1), dtype = self._state_dtype, order = 'F')
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._pfile, 'spin')
        libstitch.write_block (self._pfile, field_id, self._t0, self._pb1, su)
        # read in from one process to see if it is correct
        x_increment = self._nd // 2 + (self._nd * self._npx)
        whole_su = numpy.empty (shape = get_block_width (self._pbu), dtype = self._state_dtype, order = 'F')
        if (self._rank == 0):
            whole_su.fill (-1)
            #print (self._npx * self._nd)
            #print ((self._npx - 1) * self._nd + x_increment)
            whole_su [0:(self._npx * self._nd),:,:] = 1
            #whole_su [(self._npx * self._nd):((self._npx - 1) * self._nd + x_increment),:,:] = 1
            (rc, trial_su, is_new_time) = libstitch.read_block (self._pfile, field_id, self._t0, self._pbu)
            #print ('whole_su ', whole_su)
            #print ('trial_su ', trial_su)
            self.assertTrue ((trial_su == whole_su).all ())
        else:
            libstitch.read_block (self._pfile, field_id, -1, self._pbu) # read is collective
        self._pcomm.Barrier ()

        # shift and write again
        su1 = numpy.empty (shape = get_block_width (self._pb2), dtype = self._state_dtype, order = 'F')
        su1.fill (2)
        libstitch.write_block (self._pfile, field_id, self._t1, self._pb2, su1)
        # read in from one process to see if it is correct
        if (self._rank == 0):
            x_min = (self._npx -1) * self._nd + (self._nd // 2)
            x_max = x_min + (self._npx * self._nd)
            #print ('x_min ', x_min, ' x_max ', x_max)
            whole_su [x_min:x_max, : , :] = 2
            (rc, trial_su1, is_new_time) = libstitch.read_block (self._pfile, field_id, self._t1, self._pbu)
            #print ('whole_su ')
            #print (whole_su)
            #print ('trial_su1 ')
            #print (trial_su1)
            self.assertTrue ((trial_su1 == whole_su).all ())
        else:
            libstitch.read_block (self._pfile, field_id, -1, self._pbu) # read is collective
        self._pcomm.Barrier ()
        pass

    def test_i_multifield_write_read (self):
        # create three fields (one of each type), write each, read each
        (rc, field_id_i32) = libstitch.create_field (self._file, 'field_int32', 1, 1, -1)
        (rc, field_id_i64) = libstitch.create_field (self._file, 'field_int64', 2, 1, -1)
        (rc, field_id_f64) = libstitch.create_field (self._file, 'field_float64', 3, 1, -1.0)

        # setup for getting the block where data is new?
        libstitch.write_block (self._file, field_id_i32, self._t1, self._b4, self._s4)
        libstitch.write_block (self._file, field_id_i64, self._t1, self._b4, self._s5)
        libstitch.write_block (self._file, field_id_f64, self._t1, self._b4, self._s6)

        (rc, trial_su1_i32, is_new_time) = libstitch.read_block (self._file, field_id_i32, self._t1, self._b4)
        (rc, trial_su1_i64, is_new_time) = libstitch.read_block (self._file, field_id_i64, self._t1, self._b4)
        (rc, trial_su1_f64, is_new_time) = libstitch.read_block (self._file, field_id_f64, self._t1, self._b4)

        self.assertTrue ((trial_su1_i32 == self._s4).all ())
        self.assertTrue ((trial_su1_i64 == self._s5).all ())
        self.assertTrue ((trial_su1_f64 == self._s6).all ())
        pass

if __name__ == "__main__":
        parser = argparse.ArgumentParser (description='npx npy npz for decomposition')
        parser.add_argument ('nps', metavar='N', type=int, nargs=3, help='procs per dimension')
        parser.add_argument ('unittest_args', nargs='*')
        args = parser.parse_args ()
        # print (args)
        sys.argv[1:] = args.unittest_args

        unit_fixed_size_cv_write_read.NPX = args.nps [0]
        unit_fixed_size_cv_write_read.NPY = args.nps [1]
        unit_fixed_size_cv_write_read.NPZ = args.nps [2]

        unit_variable_size_cv_write_read.NPX = args.nps [0]
        unit_variable_size_cv_write_read.NPY = args.nps [1]
        unit_variable_size_cv_write_read.NPZ = args.nps [2]

        unittest.main()
