Stitch IO
=========

by Jay Lofstead (gflofst@sandia.gov) and John Mitchell (jamitch@sandia.gov)
with consulting by Steve Plimpton

This project aims to develop a new, efficient IO library optimized for the
style of operation found in the SPPARKS kinetic monte carlo applications
(http://spparks.sandia.gov).  This novel library was motivated by additive 
manufacturing applications where material is incrementally added to a part; 
using the 'Stitch' libary, data/material is appended to a file much in the same 
way material is incrementally added to a part; in this way, simulations of the 
additive manufacturing process are conducted locally, ie only on that portion 
of the build which is active.  This approach to IO is very efficient computationally 
since only the active portion of build is computed on; furthermore, the output 
database is extremely compact since state is only written to the file when 
it changes.


Build 
-----

There are two different entities which may be built.

#. libstitch.a; for use with c/c++ applications; can be built for serial
   or parallel calculations.

#. python module; must have python 3, numpy, scipy.  

Below, STITCH denotes the directory where 'Stitch' is cloned.

Build the 'stitch' python module
--------------------------------

The 'stitch' python module allows for creating, reading and writing 'stitch'
files.  It is particularly useful for reading 'stitch' files written by
'spparks'.

The python module should be built with gnu gcc.  It may be possible to do
otherwise but that has not been tested.  In general, the same compiler used to
build python, numpy, and scipy should be used to build the stitch module.

Python 3, numpy, and scipy are required for building the 'stitch' python
module.  This documentation was developed with the following versions.

* gnu gcc 9.2.0
* Python 3.6.9
* numpy '1.19.2'
* scipy '1.5.2'

Other versions for all of the above will probably work provided 
python 3 is used.

The python module is generally built for serial applications that post process
stitch files; this is the default and is recommended.  gcc or other compatible
compilers are used to build the python module.  The python module may be built
with an mpi compiler but will remain a serial python interface to a stitch
file.  

Change to stitch directory 

::

   % cd $STITCH

In the root project directory, create a setup.cfg file appropriate for your
system. Look at the other examples for guidance. This is important for the
python build setup.  The file 'setup.flamer.cfg' and others are included in
this directory; copy one of these to setup.cfg and edit appropriately; specify
where the python module should be installed -- then add this path to your
PYTHONPATH.

Build and install the python module.

::

    % python setup.py build
    % python setup.py install

Some build 'warnings' may occur which are generally harmless.  On the other
hand 'errors' are not good and must be addressed.  

Verify install and PYTHONPATH are correctly set by launching python and
importing the stitch module. Make sure not to do this in the source directory.

If the following command succeeds, python module build and install is
probably good.  Start an interactive python session by launching python 
and import the stitch module:

::

    % python
    % from stitch.libstitch import libstitch

As final verification of python module correctness, run the unit test; remove
all stitch files prior to running the test -- otherwise the test will fail.

:: 

    % cd $STITCH/integrated_test/weld/potts/
    % rm *.st


Run unit test

::

    % python unit_cv_readwrite.py


Screen output from test should look something like:

::

   .........
   ----------------------------------------------------------------------
   Ran 9 tests in 9.383s
   
   OK


The python module may be built for parallel python applications 
using mpi4py -- this works but is not as well
exercised and tested and is not recommended as its still an area of
development.  To build a parallel Python module:

* Build and/or install mpi4py; make sure import mpi4py works
* Add -DSTITCH_PARALLEL to the extra_compiler_args list in libstitch/setup.py
* build and test python module as described for the serial python module build
* See Makefile on howto run parallel python tests
* Should pass all tests.


Build the 'stitch' library for c/c++ application development
------------------------------------------------------------

::

    % cd $STITCH/libstitch


It is possible to build a serial C library or a Parallel Python library.  

Edit path to c compiler by appropriately defining CC in the 'Makefile'.  Use
the same mpi compiler to build the stitch library as will be used for the
parallel application.  Use -DSTITCH_PARALLEL macro for parallel applications
such as the SPPARKS Kinetic Monte Carlo simulation framework.  To build a
serial C library edit the Makefile and remove -DSTITCH_PARALLEL from both
stitch_test and stitch.lib.

Build libstitch.a suitable for linking with c/c++ applications; 

::

    % make stitch.lib  

Build 'stitch_test' written in C; this test exercises the c/c++ parallel API;

::

    % make stitch_test 

To test the C API correctness (parallel):

::

    % mpiexec -np 4 stitch_test

Depending on the mode selected, it will either test (default), generate output
to the console, or write to files to configure a new testing dataset. For the
console output, it should just output a bunch of 3-d box corner pairs and the
values for proc 0. If this happens and it looks right, it is probably working
correctly. The code is configured to be adaptable, but is currently assuming a
2x2x1 process layout per timestep with each process having a 2x2x3 chunk of the
domain. It increments in x and y by 3 and z by 1. Since there are no
incorporated automated testing, it can be scaled and visually tested by
adapting the internal parameters. There are four sets of relevant parameters:
the global (total) domain size, the process count in each dimension, how big
the decomposition is in each dimension, and how far to move when incrementing
in that dimension.


TROUBLESHOOTING
---------------

* Stitch's file format (SQLite3) relies on the file system to properly handle
  POSIX locks to manage concurrency. GPFS mostly to completely does this
  properly. Lustre does not do this at all. To work around this issue, the
  parallel version uses token passing to serialize access. This is created with
  -DSTITCH_PARALLEL. If there are issues related to this, please contact us for
  an additional workaround.

* The library is written to handle concurrency by retrying. It may fail, but it
  should not. There are a lot of error messages that include line numbers,
  function names and other information. Please include that info with any bug
  report.
