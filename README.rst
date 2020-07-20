Additive-IO
===========

by Jay Lofstead, John Mitchell, and Steve Plimpton

This project aims to develop an new, efficient IO library optimized for the
style of operation found in the SPPARKS kinetic monte carlo application
(http://spparks.sandia.gov)

The initial design point is using a key-value store for each active region
written with a metadata engine to track all of the basic information and a
generative naming system for finding data later.

The tests and Python interface all assume a current Python 3.x implementation.

The current build setup is fragile, but works. Here is how to make it work:

1. Go to the libstitch directory
2. 'make stitch_module' will create a serial version of the library for Python use and install it in the specified installation directory.
3. 'make stitch.lib' builds a parallel C interface suitable for linking with anything.
4. 'make stitch_test' builds the parallel C test harness that tests the API, but not the data.

To test the library correctness (serial interface):

1. build the Python version (serial)
2. go into integrated_test/weld/potts
3. make clean ; make

This should run clean.

To test the C API correctness (parallel):

1. build the stitch_test test harness
2. mpiexec -np 2 stitch_test 2 1 1

This should just output a bunch of 3-d box corner pairs without any other output. If this happens, it is working correctly.

It is possible to build a serial C library or a Parallel Python library.

To build a serial C library:

1. edit the Makefile to remove the -DSTITCH_PARALLEL from both stitch_test and stitch.lib.
2. make stitch.lib ; make stitch_test
3. stitch_test

To build a parallel Python library:

1. Figure out your mpi4py import statement.
2. edit setup.py to change the mpi4py import to be what is correct for you.
3. Add -DSTITCH_PARALLEL to the extra_compiler_args list.
4. Make stitch_module
5. go to integrated_tests/weld/potts
6. edit the Makefile
7. comment the serial test case out
8. uncomment the two parallel test cases
9. make
10. This should pass all tests.
