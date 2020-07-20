/* Copyright 2019 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525
 * with NTESS, the U.S. Government retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * For more information, contact Jay Lofstead (gflofst@sandeia.gov) or
 * John Mitchell (jamitch@sandia.gov) for more information.
 */ 
#ifndef TYPEDEFS_H
#define TYPEDEFS_H

/*
 * Following two includes must be at the top 
 */
#define NPY_NO_DEPRECATED_API NPY_1_8_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

/*
 * numpy 'enumerated types'
 * MUST BE EXACT correspondence between the enumerated types
 * and the c-type names
 */

// enumerated types
#define NUMPY_ORDINAL NPY_UINT32
#define NUMPY_REAL NPY_FLOAT64
// c-type names
typedef npy_uint32 ordinal; 
typedef npy_float64 real;

#endif
