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
