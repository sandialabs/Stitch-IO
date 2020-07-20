import unittest
import numpy
from math import fabs,pi,cos,sin,log
import json
from numpy import array, ndarray, dot, transpose
from stitch.libstitch import libstitch
from stitch.weld.app_scripts import haz_box as hb
from stitch import time_equivalence

class unit_haz_box(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        name='unit_haz_box'

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

        # Read main output section from json
        output=json_params["output"]
        self._prefix=output['prefix']
        self._dS=output['spparks_dump_dt']
        self._dump_abs_tol=output["spparks_dump_abs_tol"]

        # Read main parameters section from json
        params=json_params["parameters"]
        self._alpha_0=params["alpha"]
        self._weld_speed=params["weld_speed"]
        self._weld_distance=params["weld_distance"]
        case=params["case"]
        pool_width=params["pool_width"]
        plate_thickness=params["plate_thickness"]
        haz=params["haz"]
        haz_cushion=params["haz_cushion"]
        self._haz_box=hb.HazBox(case,pool_width,plate_thickness,haz,haz_cushion)
        pass

    def test_alpha(self):
        haz_box=self._haz_box
        ypn=self._pool_position
        alpha_0=self._alpha_0
        vp=self._weld_speed
        dS=self._dS
        alpha=haz_box.compute_alpha(ypn,alpha_0,vp,dS)
        print("Input alpha_0=%3.1f, alpha=%3.1f"%(alpha_0,alpha,))

    def test_play(self):
        pass


if __name__ == "__main__":
    unittest.main()
