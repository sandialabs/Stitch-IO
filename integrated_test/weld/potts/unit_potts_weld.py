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
from math import fabs,pi,cos,sin,log
import json
from numpy import array, ndarray, dot, transpose
from mpi4py import MPI
from stitch.libstitch import libstitch
from stitch.weld.app_scripts import haz_box as hb
from stitch import time_equivalence
import unit_cv_readwrite as ucvr

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
    x_str="x0,x1,y0,y1,z0,z1={0:d},{1:d},{2:d},{3:d},{4:d},{5:d}\n"
    block_str=x_str.format(block[0,0],block[1,0],block[0,1],block[1,1],block[0,2],block[1,2])
    return block_str

class unit_potts_weld(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        name='unit_potts_weld'

        # Create 'stitch' file
        # 'suffix' 'st' for 'stitch file'
        self._fname=name+'.st'
        self._state_dtype=numpy.int32
        self._comm = MPI.COMM_WORLD
        (self._rc,self._fileid) = libstitch.open (self._comm, self._fname)

        # Open parameters file; load json
        json_params_file=name+".in"
        inpfile=open(json_params_file,'r')
        # Read json parameter file
        json_params=json.load(inpfile)
        inpfile.close()

        yp=0.0
        t0=0.0
        start_up=True
        self._json_params_file=json_params_file
        self._pool_position=yp
        self._t0=t0
        self._start_up=start_up
        self._absolute_tolerance = 1.0e-9
        self._relative_tolerance = 1.0e-15
        self._no_value_present = -1

        # Read main output section from json
        output=json_params["output"]
        self._prefix=output['prefix']
        self._dS=output['spparks_dump_dt']
        self._dump_abs_tol=output["spparks_dump_abs_tol"]

        # Read main parameters section from json
        params=json_params["parameters"]
        self._alpha=params["alpha"]
        self._weld_speed=params["weld_speed"]
        self._weld_distance=params["weld_distance"]
        case=params["case"]
        pool_width=params["pool_width"]
        plate_thickness=params["plate_thickness"]
        haz=params["haz"]
        haz_cushion=params["haz_cushion"]
        self._haz_box=hb.HazBox(case,pool_width,plate_thickness,haz,haz_cushion)
        rc = libstitch.set_parameters (self._fileid, self._absolute_tolerance, self._relative_tolerance, self._no_value_present)
        # print("Case=%s, pool_width=%5.1f, plate_thickness=%3d, HAZ=%3d, haz_cushion=%2d"%(case,pool_width,plate_thickness,haz,haz_cushion))
        pass

    def tearDown(self):
        # close the file
        libstitch.close (self._fileid);
        pass
    
    def test_run_weld(self):
        # Time eqivalence object
        absolute_tolerance=self._absolute_tolerance
        relative_tolerance=self._relative_tolerance
        eq=time_equivalence.Equivalence(absolute_tolerance,relative_tolerance)

        # Setup for stitch simulation
        t0=self._t0
        pool_position=self._pool_position
        alpha=self._alpha
        weld_speed=self._weld_speed
        start_up=self._start_up
        haz_box=self._haz_box
        weld_distance=self._weld_distance
        dS=self._dS
        distance_traveled=0.0
        yn=pool_position
        tn=t0
        times=[tn]
        # Run stitch simulation
        stage=0
        # 1 = STITCH_INT32. Need to create an Enum and work that over.
        (rc, field_id) = libstitch.query_field (self._fileid, 'spin', 1, 1, -1)
        while distance_traveled<weld_distance:
            tnp1,ynp1,potts_init_cv,cv=haz_box.get_cv_boxes(tn,dS,yn,alpha,weld_speed,start_up=start_up)

            #print("tn=%5.1f, potts_init_bb: %s"%(tn,hb.pretty_string_bb_box(potts_init_cv),))
            #print("tn=%5.1f, weld_cv: %s"%(tn,hb.pretty_string_bb_box(cv),))

            # Write potts init microstructure
            potts_value=66
            #print (potts_init_cv)
            potts_block=hb.convert_bb_to_libstitch_block(potts_init_cv)
            potts_block_width=ucvr.get_block_width(potts_block)
            potts_state=potts_value*numpy.ones(shape=potts_block_width,dtype=self._state_dtype,order='F')
            (rc, new_time) = libstitch.write_block (self._fileid, field_id, t0, potts_block, potts_state);

            # Read potts data and assert correctness
            (rc, trial_state, new_time) = libstitch.read_block (self._fileid, field_id, t0, potts_block)
            #print ('potts_block')
            #print (potts_block)
            #print ('trial_state')
            #print (trial_state)
            self.assertTrue((trial_state == potts_state).all())

            # Initialize 'weld' state on weld block; should reflect initial state as written above 
            #   on potts_block over the intersection of 'potts_block' with 'weld_block'
            weld_block=hb.convert_bb_to_libstitch_block(cv)
            weld_block_width=ucvr.get_block_width(weld_block)
            #print (tn)
            (rc, weld_initial_state, new_time) = libstitch.read_block (self._fileid, field_id, tn, weld_block)
            # Assert that 'start time t0' and 'last time tn' are what is returned by libstitch
            (rc, _abs_tol, _rel_tol, _t0, _tn)=libstitch.get_parameters(self._fileid)
            #print("Read stage=%d, _t0=%6.3f, _tn=%6.3f, expected t0=%6.3f, tn=%6.3f\n"%(stage,_t0,_tn,t0,tn,))
            stage+=1
            self.assertTrue(eq.equivalence(_tn,tn))
            self.assertTrue(eq.equivalence(_t0,t0))

            w1x=weld_block_width[0]-potts_block_width[0]
            w1y=weld_block_width[1]-potts_block_width[1]
            w1z=weld_block_width[2]-potts_block_width[2]
            #print("Size of non-overlap block=%d,%d,%d"%(w1x,w1y,w1z,))

            # Weld is moving along y-axis; boxes overlap fully in x and z; therefore 0==w1x, 0==w1y
            weld_value=99
            self.assertTrue(0==w1x)
            self.assertTrue(0==w1z)
            for k in range(0,weld_block_width[2]):
                for j in range(0,w1y):
                    for i in range(0,weld_block_width[0]):
                        self.assertTrue(weld_value==weld_initial_state[i,j,k])

            # Assert potts init section just written on this iteration
            for k in range(0,weld_block_width[2]):
                for j in range(w1y,weld_block_width[1]):
                    for i in range(0,weld_block_width[0]):
                        self.assertTrue(potts_value==weld_initial_state[i,j,k])

            # Now write weld data
            weld_state=weld_value*numpy.ones(shape=weld_block_width,dtype=self._state_dtype,order='F')
            times.append(tnp1)
            (rc, new_time) = libstitch.write_block (self._fileid, field_id, tnp1, weld_block, weld_state);

            start_up=False
            distance_traveled+=(ynp1-yn)
            tn=tnp1
            yn=ynp1

        # Assert 'timestamps' and 'times'
        rc,stitch_times=libstitch.get_times(self._fileid)
        self.assertTrue(len(times)==len(stitch_times))
        for i,(t,) in enumerate(zip(times)):
            self.assertTrue(fabs(t-stitch_times[i])<absolute_tolerance)
            self.assertTrue(eq.equivalence(stitch_times[i],t))

if __name__ == "__main__":
    unittest.main()
