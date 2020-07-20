/*
 * On the MAC, gcc49 has a sensitivity to the order in which these 
 * headers are included.  The first 2 headers below must be at 
 * the top otherwise on MAC with gcc49, the build fails.
 */
#include "femanica_typedefs.h"
#include "weld_geometry.h"
#include <string>
#include <cstddef>
#include <math.h>
#include <stdio.h>

extern "C" { 


using std::string;

static ordinal num_instances=0;
typedef weld::pool_shape::EllipticBezier EllipticBezier;

void destructor_callback(PyObject *capsule){
   const char *name=PyCapsule_GetName(capsule);
   EllipticBezier *e = (EllipticBezier*)PyCapsule_GetPointer(capsule, name);
   printf("Capsule delete: %s \n",name);

   // delete
   delete e;
   free((void*)name);
}

void* get_capsule_pointer(PyObject *capsule){
   /*
    * Get capsule name
    */
   void *ptr=nullptr;
   const char *capsule_name=nullptr;
   capsule_name=PyCapsule_GetName(capsule);
   if (nullptr==capsule_name) goto fail;

   /*
    * Get object pointer
    */
   ptr=PyCapsule_GetPointer(capsule,capsule_name);
   if (nullptr==ptr) goto fail;

   return ptr;

   fail:
      printf("Failed to 'get_capsule_pointer'");
      return nullptr;
}


static 
PyObject* _get_ellipse(PyObject *self, PyObject *args, PyObject *keywds){

   real a,b,T,alpha,beta;
   char *capsule_name=nullptr;
   PyObject *capsule=nullptr;
   EllipticBezier *ellipse=nullptr;

   const char *kwlist[] = {"pool_width", "pool_length", 
                           "plate_thickness","scale_param","interpolate_param",
                           nullptr};

   if (!PyArg_ParseTupleAndKeywords(args, keywds, "ddddd", 
            const_cast<char **>(kwlist), &a,&b,&T,&alpha,&beta)) {
      return nullptr;
   }

   {
      // Allocate 'capsule_name'
      const char *name="EllipticBezier_";
      // Increment number of objects allocated;
      // Use num_instances as part of unique 'capsule_name'
      num_instances+=1;
      // Buffer to store string representation of unique identifier 
      char _num_instances[33];
      sprintf(_num_instances,"%d",num_instances);
      // add '1' for terminating null character
      size_t name_len=strlen(name)+strlen(_num_instances)+1;
      capsule_name=(char*)calloc(name_len,sizeof(char));
      if(0==capsule_name) {
         printf("Capsule name allocation FAILED: %s \n",name);
         goto fail;
      }
      strcat(capsule_name,name);
      strcat(capsule_name,_num_instances);
   }

   {
      ellipse=new EllipticBezier(a,b,T,alpha,beta);
      if(nullptr==ellipse){
         printf("\'EllipticBezier\' allocation %s FAILED\n",capsule_name);
         free(capsule_name);
         goto fail;
      }
      // If successful, creates 'capsule' with 'capsule->ob_refcnt=1'
      capsule=(PyObject*)PyCapsule_New((void *)ellipse,capsule_name,destructor_callback);
      if (NULL==capsule) {
         // Delete tree since capsule construction failed
         // free pointer to capsule name
         delete ellipse;
         free(capsule_name);
         goto fail;
      }
   }
   return capsule;

   fail:
      Py_XDECREF(capsule);
      Py_INCREF(Py_None);
      return Py_None;
}

static 
PyObject* _distance(PyObject *self, PyObject *args){

   real d;
   PyArrayObject *xyz_array=nullptr;
   PyObject *xyz_object=nullptr, *capsule=nullptr;


   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"OO",&capsule,&xyz_object))
   { return nullptr; }

   // Creates a new reference which must be decremented at end
   // unless you are returning or keeping
   xyz_array = (PyArrayObject*)PyArray_FROM_OTF(xyz_object,NUMPY_REAL,NPY_ARRAY_CARRAY);
   if (nullptr==xyz_array || 1<PyArray_NDIM(xyz_array) || 3!=PyArray_DIMS(xyz_array)[0]) {
      const char *m="Input 'xyz' array has nullptr or shape invalid; expected shape=(3,)";
      PyErr_SetString(PyExc_ValueError,m);
      PyErr_PrintEx(0); 
      goto fail;
   }

   {
      const real *xyz=(real *)PyArray_DATA(xyz_array);
      const EllipticBezier *eb=(const EllipticBezier*)get_capsule_pointer(capsule);
      d=eb->distance(xyz);
      return Py_BuildValue("d",d);
   }

   fail:
      Py_XDECREF(xyz_array);
      Py_INCREF(Py_None);
      return Py_None;
}

static PyMethodDef _weld_geometry_methods[] = {
   {"_get_ellipse",  
     (PyCFunction)_get_ellipse,
     METH_VARARGS|METH_KEYWORDS,
     "Creates object representing weld pool surface.\n\
      Args \n\
      -----\n\
      real pool_width: dimension of minor elliptical axis (top surface)\n\
      real pool_length: dimension of major elliptical axis (top surface)\n\
      real plate_thickness: dimension of plate thickness\n\
      real scale_param: used to create bottom ellipse by scaling \n\
                        ellipse on top surface given by (pool_width, pool_length)\n\
      real interpolate_param: used to create bottom ellipse by scaling \n\
                            ellipse on top surface given by (pool_width, pool_length)\n\
      returns \n\
      -------\n\
      PyCapsule storing c++ pointer to EllipseBezier object\n"},
   {"_distance",  
     (PyCFunction)_distance,
      METH_VARARGS,
     "Computes distance to 'EllipseBezier surface.\n\
      Args \n\
      -----\n\
      param PyCapsule: storing c++ pointer to EllipseBezier object\n\
      param ndarray: dtype=numpy.float64,shape=(3,); coordinates of point\n\
                     outside of pool for which closest point projection\n\
                     is calculated.\n\
      returns \n\
      -------\n\
      Real scalar distance to surface; If point inside pool return value=-1.0\n"},
   {nullptr, nullptr, 0, nullptr}        /* Sentinel */
};

static struct PyModuleDef _weld_geometrymodule = {
   PyModuleDef_HEAD_INIT,
   "_weld_geometry",
   NULL,
   -1,
   _weld_geometry_methods
};

PyMODINIT_FUNC
PyInit__weld_geometry(void)
{
   PyObject *m=NULL;
   m=PyModule_Create(&_weld_geometrymodule);
	import_array();
   return m;
};

};
