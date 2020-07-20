
Build and Install 'weld_stitch' Package
=======================================

Requires
++++++++

* python 3; has been exercised using 3.4.3
* numpy; has been exercised using 1.9.2
* scipy; has been exercised using 0.14.1
* matplotlib; has been exercised using 1.4.3
* ipython; not absolutely required but...; used with 3.2.1
* c++ compiler that can support c++11; has been exercised with gcc 5.3.1; 

Build
+++++

* cd 'weld_stitch'
* 'python setup.py build'

Install
+++++++

* cd 'weld_stitch'
* edit setup.cfg (don't check in your edits please)
* 'python setup.py install'

Using 'weld_stitch'
++++++++++++++++++

* add 'weld_stitch' to your  PYTHONPATH 

.. code-block:: bash

        # WELD_STITCH
        WELD_STITCH=$HOME/local/weld_stitch
        PYTHONPATH=$WELD_STITCH:$PYTHONPATH


* launch 'ipython' and import python interface to elements of teardrop weld pool shape


::

        [ jamitch @ flamer :]ipython
        Python 3.4.3 (default, Aug  9 2016, 15:36:17) 
        Type "copyright", "credits" or "license" for more information.

        IPython 3.2.1 -- An enhanced Interactive Python.
        ?         -> Introduction and overview of IPython's features.
        %quickref -> Quick reference.
        help      -> Python's own help system.
        object?   -> Details about 'object', use 'object??' for extra details.

        In [1]: import weld_stitch.weld.bezier_teardrop as bt
