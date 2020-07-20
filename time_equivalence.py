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

from math import fabs

def eqivalence(t,t0,abs_tol=1.0e-9,rel_tol=1.0e-15):
    """
    Compare an input value t0 with an existing value t from database.
    Returns True if t and t0 are numerically equivalent.

    INPUTS
    ------
    t: floating point value already existing in database
    t0: input floating point value
    abs_tol: absolute tolerance
    relative_tol: relative tolerance

    OUTPUT
    ------
    boolean
    True if t and t0 are to be taken as numerically equivalent
    False if t and t0 are NOT numerically equivalent
    """
    if t0>=0:
        if fabs(t-t0)<abs_tol+t0*rel_tol:
            return True
        else:
            return False
    else:
        t_str="{0:20.10e}".format(t0)
        raise TypeError("Input time value must be greater than or equal to 0.0; input value="+t_str)

class Equivalence(object):
    def __init__(self,abs_tol=1.0e-9,rel_tol=1.0e-15):
        self._absolute_tolerance=abs_tol
        self._relative_tolerance=rel_tol

    def get_absolute_tolerance(self):
        return self._absolute_tolerance
    abs_tol=property(get_absolute_tolerance,doc="absolute tolerance used for eqivalence between to time values")

    def get_relative_tolerance(self):
        return self._relative_tolerance
    rel_tol=property(get_relative_tolerance,doc="relative tolerance used for eqivalence between to time values")

    def equivalence(self,t,t0):
        return eqivalence(t,t0,abs_tol=self.abs_tol,rel_tol=self.rel_tol)
