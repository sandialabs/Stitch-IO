#!/usr/bin/env python3
# encoding: utf-8
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

import os
from distutils.core import setup, Extension
import numpy
import mpi4py

mpi_compile_args = os.popen("mpicc --showme:compile").read().strip().split(' ')
mpi_link_args    = os.popen("mpicc --showme:link").read().strip().split(' ')

# the c++ extension module
stitch_mod = Extension ('libstitch'
                       ,['stitchmodule.c', 'stitch.c', 'sqlite3.c']
                       ,include_dirs=[numpy.get_include (), mpi4py.get_include ()]
                       ,extra_compile_args = mpi_compile_args
                       ,extra_link_args    = mpi_link_args
                       )

setup (name = 'libstitch'
      ,version='0.1.0'
      ,description='Python interface (module libstitch) to c-based stitch io system.'
      ,ext_modules=[stitch_mod]
      ,include_dirs=[numpy.get_include (), mpi4py.get_include ()]
      )

