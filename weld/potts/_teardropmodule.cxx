/*
 * On the MAC, gcc49 has a sensitivity to the order in which these 
 * headers are included.  The first 2 headers below must be at 
 * the top otherwise on MAC with gcc49, the build fails.
 */
#include "femanica_typedefs.h"
#include "teardrop.h"
#include <string>
#include <string>
#include <cstddef>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <tuple>

//extern "C" { 

using std::string;
using std::vector;
static ordinal num_instances=0;
typedef weld::pool_shape::BezierCurve BezierCurve;
using weld::pool_shape::Teardrop2D;

real dot(const real* y1, const real *y2) { 
return y1[0]*y2[0]+y1[1]*y2[1];
}

template<class value_type>
void destructor_callback(PyObject *capsule){
   const char *name=PyCapsule_GetName(capsule);
   value_type *b = (value_type*)PyCapsule_GetPointer(capsule, name);
   //printf("Capsule delete: %s \n",name);

   // delete
   delete b;
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
PyObject* _teardrop2d_parametric_coordinate_at_maximum_pool_width(PyObject *self, PyObject *args){
   real umax;
   PyObject *capsule=nullptr;
   PyObject *result=nullptr;

   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"O",&capsule))
   { return nullptr; }

   {
     const Teardrop2D *t=(const Teardrop2D*)get_capsule_pointer(capsule);
     umax=t->get_parametric_coordinate_at_maximum_pool_width();
     result=Py_BuildValue("d",umax);
     if(nullptr==result) goto fail; 
     return result;

   fail:
      Py_XDECREF(result);
      Py_INCREF(Py_None);
      return Py_None;
   }
}

static 
PyObject* _teardrop2d_tangent_to_curve(PyObject *self, PyObject *args){
   real u;
   PyObject *capsule=nullptr;
   PyArrayObject *rho_array=nullptr, *rho_1_array=nullptr;
   PyObject *result=nullptr;

   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"dO",&u,&capsule))
   { return nullptr; }

   {
     const Teardrop2D *t=(const Teardrop2D*)get_capsule_pointer(capsule);
     // Compute point on curve at 'u' as well as tangent
     real rho[]={0,0};
     real rho_1[]={0,0};
     t->compute_curve(u,rho,rho_1);

     const npy_intp dim=2;
     npy_intp num_dim=1;
     npy_intp dims[] = { dim };
     // New reference
     rho_array = (PyArrayObject*)PyArray_SimpleNew(num_dim, dims,NUMPY_REAL);
     rho_1_array = (PyArrayObject*)PyArray_SimpleNew(num_dim, dims,NUMPY_REAL);
     if (nullptr==rho_array || nullptr==rho_1_array) goto fail;

     {
        real *_rho=(real*)PyArray_DATA(rho_array);
        real *_rho_1=(real*)PyArray_DATA(rho_1_array);
        for(npy_intp i=0;i<dim;i++){
           _rho[i]=rho[i];
           _rho_1[i]=rho_1[i];
        }
     }

     result=Py_BuildValue("NN",rho_array,rho_1_array);
     if(nullptr==result) goto fail; 
     return result;

   fail:
      Py_XDECREF(result);
      Py_XDECREF(rho_array);
      Py_XDECREF(rho_1_array);
      Py_INCREF(Py_None);
      return Py_None;
   }
}

static 
PyObject* _teardrop2d_predictor(PyObject *self, PyObject *args){
   real x,y;
   PyObject *capsule=nullptr;
   PyArrayObject *rho_array=nullptr, *rho_1_array=nullptr;
   PyObject *result=nullptr;

   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"ddO",&x,&y,&capsule))
   { return nullptr; }

   {
     const Teardrop2D *t=(const Teardrop2D*)get_capsule_pointer(capsule);
     // Run predictor
     real du=0.025;
     real u=t->predict_cpp_parametric_coordinate(x,y,du);
     // Compute point on curve at 'u'
     real rho[]={0,0};
     real rho_1[]={0,0};
     t->compute_curve(u,rho,rho_1);

     const npy_intp dim=2;
     npy_intp num_dim=1;
     npy_intp dims[] = { dim };
     // New reference
     rho_array = (PyArrayObject*)PyArray_SimpleNew(num_dim, dims,NUMPY_REAL);
     rho_1_array = (PyArrayObject*)PyArray_SimpleNew(num_dim, dims,NUMPY_REAL);
     if (nullptr==rho_array || nullptr==rho_1_array) goto fail;

     {
        real *_rho=(real*)PyArray_DATA(rho_array);
        real *_rho_1=(real*)PyArray_DATA(rho_1_array);
        for(npy_intp i=0;i<dim;i++){
           _rho[i]=rho[i];
           _rho_1[i]=rho_1[i];
        }
     }

     result=Py_BuildValue("NN",rho_array,rho_1_array);
     if(nullptr==result) goto fail; 
     return result;

   fail:
      Py_XDECREF(result);
      Py_XDECREF(rho_array);
      Py_XDECREF(rho_1_array);
      Py_INCREF(Py_None);
      return Py_None;
   }
}

static 
PyObject* _get_teardrop2d_iterate(PyObject *self, PyObject *args){

   real x,y,u;
   PyObject *capsule=nullptr;
   PyArrayObject *rho_array=nullptr, *dr_array=nullptr, *rho_1_array=nullptr;
   PyArrayObject *rho_11_array=nullptr;
   PyObject *result=nullptr;

   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"dddO",&x,&y,&u,&capsule))
   { return nullptr; }

   {
     const Teardrop2D *t=(const Teardrop2D*)get_capsule_pointer(capsule);
     // Compute point on curve at 'u' as well as tangent and curvature
     real rho[]={0,0};
     real rho_1[]={0,0};
     real rho_11[]={0,0};
     t->compute_curve(u,rho,rho_1,rho_11);
     // Compute vector from input point to curve
     real dr[]={x-rho[0],y-rho[1]};
     // Jacobian
     real j11=dot(rho_11,dr)-dot(rho_1,rho_1);

     const npy_intp dim=2;
     npy_intp num_dim=1;
     npy_intp dims[] = { dim };
     // New reference
     rho_array = (PyArrayObject*)PyArray_SimpleNew(num_dim, dims,NUMPY_REAL);
     dr_array = (PyArrayObject*)PyArray_SimpleNew(num_dim, dims,NUMPY_REAL);
     rho_1_array = (PyArrayObject*)PyArray_SimpleNew(num_dim, dims,NUMPY_REAL);
     rho_11_array = (PyArrayObject*)PyArray_SimpleNew(num_dim, dims,NUMPY_REAL);
     if (nullptr==rho_array || nullptr==dr_array 
		            || nullptr==rho_1_array || nullptr==rho_11_array) goto fail;

     {
        real *_rho=(real*)PyArray_DATA(rho_array);
        real *_dr=(real*)PyArray_DATA(dr_array);
        real *_rho_1=(real*)PyArray_DATA(rho_1_array);
        real *_rho_11=(real*)PyArray_DATA(rho_11_array);
        for(npy_intp i=0;i<dim;i++){
           _rho[i]=rho[i];
           _dr[i]=dr[i];
           _rho_1[i]=rho_1[i];
           _rho_11[i]=rho_11[i];
        }
     }

     result=Py_BuildValue("dNNNN",j11,rho_array,dr_array,rho_1_array,rho_11_array);
     if(nullptr==result) goto fail; 
     return result;

  }


   fail:
      Py_XDECREF(result);
      Py_XDECREF(rho_array);
      Py_XDECREF(dr_array);
      Py_XDECREF(rho_1_array);
      Py_XDECREF(rho_11_array);
      Py_INCREF(Py_None);
      return Py_None;
}


static 
PyObject* _get_teardrop2d(PyObject *self, PyObject *args){

   real pool_width;
   PyArrayObject *control_points_array=nullptr;
   PyObject *control_points_object=nullptr;
   char *capsule_name=nullptr;
   PyObject *capsule=nullptr;
   Teardrop2D *teardrop2d=nullptr;
   npy_intp num_control_points=0;
   npy_intp dim=0;

   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"Od",&control_points_object,&pool_width))
   { return nullptr; }

   // Creates a new reference which must be decremented at end
   // unless you are returning or keeping
   control_points_array = (PyArrayObject*)PyArray_FROM_OTF(control_points_object,NUMPY_REAL,NPY_ARRAY_CARRAY);
   if (nullptr==control_points_array || 2!=PyArray_NDIM(control_points_array)) {
      const char *m="Input 'control_points' array has nullptr or shape invalid";
      PyErr_SetString(PyExc_ValueError,m);
      PyErr_PrintEx(0); 
      goto fail;
   }
   // expect shape=(num_control_points,dim=2 or 3)
   num_control_points=PyArray_DIMS(control_points_array)[0];
   dim=PyArray_DIMS(control_points_array)[1];
   if (2!=dim && 3!=dim) {
      const char *m="Input 'control_points' must have shape=(num_control_points,dim) where dim=2 or 3.";
      PyErr_SetString(PyExc_ValueError,m);
      PyErr_PrintEx(0); 
      goto fail;
   }
   {
      // Allocate 'capsule_name'
      const char *name="Teardrop2D_";
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
      // Create std:vector of control points
      vector< vector<double> > cp(num_control_points);
      // Pointer to control points passed in as numpy array
      real *pxyz=(real*)PyArray_DATA(control_points_array);
      for(npy_intp i=0;i<num_control_points;i++,pxyz+=dim){
         vector<double> p(dim);
         for(npy_intp j=0;j<dim;j++){
            p[j]=pxyz[j];
         }
         cp[i]=p;
      }
      bool normalize_pool_position=true;
      teardrop2d=new Teardrop2D(cp,pool_width,normalize_pool_position);
      // Store 'teardrop2d' in python capsule
      if(nullptr==teardrop2d){
         printf("\'Teardrop2D\' allocation %s FAILED\n",capsule_name);
         free(capsule_name);
         goto fail;
      }
      // If successful, creates 'capsule' with 'capsule->ob_refcnt=1'
      capsule=(PyObject*)PyCapsule_New((void *)teardrop2d,capsule_name,destructor_callback<Teardrop2D>);
      if (NULL==capsule) {
         // Delete tree since capsule construction failed
         // free pointer to capsule name
         delete teardrop2d;
         free(capsule_name);
         goto fail;
      }
   }
   return capsule;

   fail:
      Py_XDECREF(control_points_array);
      Py_INCREF(Py_None);
      return Py_None;
}

static 
PyObject* _get_teardrop2d_cpp(PyObject *self, PyObject *args){

   real x,y,u,distance;
   PyObject *capsule=nullptr;
   PyObject *result=nullptr;

   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"ddO",&x,&y,&capsule))
   { return nullptr; }

   {
     const Teardrop2D *t=(const Teardrop2D*)get_capsule_pointer(capsule);
     std::tie(u,distance)=t->closest_point(x,y);
   }

   result=Py_BuildValue("dd",u,distance);
   if(nullptr==result) goto fail; 
   return result;

   fail:
      Py_XDECREF(result);
      Py_INCREF(Py_None);
      return Py_None;
}

static 
PyObject* _get_teardrop2d_size(PyObject *self, PyObject *args){

   real width,length;
   PyObject *capsule=nullptr;
   PyObject *result=nullptr;

   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"O",&capsule))
   { return nullptr; }

   {
     const Teardrop2D *t=(const Teardrop2D*)get_capsule_pointer(capsule);
     width=t->get_pool_width();
     length=t->get_pool_length();
   }

   result=Py_BuildValue("dd",width,length);
   if(nullptr==result) goto fail; 
   return result;

   fail:
      Py_XDECREF(result);
      Py_INCREF(Py_None);
      return Py_None;
}

static 
PyObject* _get_teardrop2d_position(PyObject *self, PyObject *args){

   real position;
   PyObject *capsule=nullptr;
   PyObject *result=nullptr;

   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"O",&capsule))
   { return nullptr; }

   {
     const Teardrop2D *t=(const Teardrop2D*)get_capsule_pointer(capsule);
     position=t->get_pool_position();
   }

   result=Py_BuildValue("d",position);
   if(nullptr==result) goto fail; 
   return result;

   fail:
      Py_XDECREF(result);
      Py_INCREF(Py_None);
      return Py_None;
}

static 
PyObject* _get_beziercurve(PyObject *self, PyObject *args){

   PyArrayObject *control_points_array=nullptr;
   PyObject *control_points_object=nullptr;
   char *capsule_name=nullptr;
   PyObject *capsule=nullptr;
   BezierCurve *bezier_curve=nullptr;
   npy_intp num_control_points=0;
   npy_intp dim=0;

   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"O",&control_points_object))
   { return nullptr; }

   // Creates a new reference which must be decremented at end
   // unless you are returning or keeping
   control_points_array = (PyArrayObject*)PyArray_FROM_OTF(control_points_object,NUMPY_REAL,NPY_ARRAY_CARRAY);
   if (nullptr==control_points_array || 2!=PyArray_NDIM(control_points_array)) {
      const char *m="Input 'control_points' array has nullptr or shape invalid";
      PyErr_SetString(PyExc_ValueError,m);
      PyErr_PrintEx(0); 
      goto fail;
   }
   // expect shape=(num_control_points,dim=2 or 3)
   num_control_points=PyArray_DIMS(control_points_array)[0];
   dim=PyArray_DIMS(control_points_array)[1];
   if (2!=dim && 3!=dim) {
      const char *m="Input 'control_points' must have shape=(num_control_points,dim) where dim=2 or 3.";
      PyErr_SetString(PyExc_ValueError,m);
      PyErr_PrintEx(0); 
      goto fail;
   }
   {
      // Allocate 'capsule_name'
      const char *name="BezierCurve_";
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
      // Create std:vector of control points
      vector< vector<double> > cp(num_control_points);
      // Pointer to control points passed in as numpy array
      real *pxyz=(real*)PyArray_DATA(control_points_array);
      for(npy_intp i=0;i<num_control_points;i++,pxyz+=dim){
         vector<double> p(dim);
         for(npy_intp j=0;j<dim;j++){
            p[j]=pxyz[j];
         }
         cp[i]=p;
      }
      bezier_curve=new BezierCurve(cp);
      // Store 'bezier_curve' in python capsule
      if(nullptr==bezier_curve){
         printf("\'BezierCurve\' allocation %s FAILED\n",capsule_name);
         free(capsule_name);
         goto fail;
      }
      // If successful, creates 'capsule' with 'capsule->ob_refcnt=1'
      capsule=(PyObject*)PyCapsule_New((void *)bezier_curve,capsule_name,destructor_callback<BezierCurve>);
      if (NULL==capsule) {
         // Delete tree since capsule construction failed
         // free pointer to capsule name
         delete bezier_curve;
         free(capsule_name);
         goto fail;
      }
   }
   return capsule;

   fail:
      Py_XDECREF(control_points_array);
      Py_INCREF(Py_None);
      return Py_None;
}

static 
PyObject* _teardrop2d_compute_curve(PyObject *self, PyObject *args){

   PyArrayObject *u_array=nullptr, *curve_array=nullptr;
   PyObject *u_object=nullptr, *capsule=nullptr;

   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"OO",&capsule,&u_object))
   { return nullptr; }

   // Creates a new reference which must be decremented at end
   // unless you are returning or keeping
   u_array = (PyArrayObject*)PyArray_FROM_OTF(u_object,NUMPY_REAL,NPY_ARRAY_IN_ARRAY);
   if (nullptr==u_array || 1!=PyArray_NDIM(u_array) || 0>=PyArray_DIMS(u_array)[0]) {
      const char *m="Input 'u' array has nullptr or shape invalid; expected shape=(num_points_curve,)";
      PyErr_SetString(PyExc_ValueError,m);
      PyErr_PrintEx(0); 
      goto fail;
   }
   {
      const Teardrop2D *o=(const Teardrop2D*)get_capsule_pointer(capsule);
      const npy_intp num_points_curve=PyArray_DIMS(u_array)[0];
      const npy_intp dim=2;
      npy_intp num_dim=2;
      npy_intp dims[] = { dim, num_points_curve };
      
      // New reference
      curve_array = (PyArrayObject*)PyArray_SimpleNew(num_dim, dims,NUMPY_REAL);
      if (nullptr==curve_array) goto fail;

      real *curve=(real*)PyArray_DATA(curve_array);
      real *u=(real*)PyArray_DATA(u_array);
      if(nullptr==curve || nullptr==u) goto fail;

      // Evaluate curve
      real p[dim];
      for(npy_intp i=0;i<num_points_curve;i++,u++){
         o->compute_curve(*u,p);
         for(npy_intp offset=0,j=0;j<dim;j++,offset+=num_points_curve){
            npy_intp ptr=offset+i;
            curve[ptr]=p[j];
         }
      }

      return (PyObject*) curve_array;
   }

   fail:
      Py_XDECREF(u_array);
      Py_XDECREF(curve_array);
      Py_INCREF(Py_None);
      return Py_None;
}

static 
PyObject* _beziercurve_compute_curve(PyObject *self, PyObject *args){

   PyArrayObject *u_array=nullptr, *curve_array=nullptr;
   PyObject *u_object=nullptr, *capsule=nullptr;

   /*
    * Receive c-object pointer that has not had its reference 
    * count incremented; Hence it is borrowed and hence I 
    * generally don't have to do anything with these pointers.
    */
   if (!PyArg_ParseTuple(args,"OO",&capsule,&u_object))
   { return nullptr; }

   // Creates a new reference which must be decremented at end
   // unless you are returning or keeping
   u_array = (PyArrayObject*)PyArray_FROM_OTF(u_object,NUMPY_REAL,NPY_ARRAY_IN_ARRAY);
   if (nullptr==u_array || 1!=PyArray_NDIM(u_array) || 0>=PyArray_DIMS(u_array)[0]) {
      const char *m="Input 'u' array has nullptr or shape invalid; expected shape=(num_points_curve,)";
      PyErr_SetString(PyExc_ValueError,m);
      PyErr_PrintEx(0); 
      goto fail;
   }
   {
      const BezierCurve *o=(const BezierCurve*)get_capsule_pointer(capsule);
      const npy_intp num_points_curve=PyArray_DIMS(u_array)[0];
      const npy_intp dim=o->curve_spatial_dimension();
      npy_intp num_dim=2;
      npy_intp dims[] = { dim, num_points_curve };
      
      // New reference
      curve_array = (PyArrayObject*)PyArray_SimpleNew(num_dim, dims,NUMPY_REAL);
      if (nullptr==curve_array) goto fail;

      real *curve=(real*)PyArray_DATA(curve_array);
      real *u=(real*)PyArray_DATA(u_array);
      if(nullptr==curve || nullptr==u) goto fail;

      // Evaluate curve
      real p[dim];
      for(npy_intp i=0;i<num_points_curve;i++,u++){
         o->compute_curve(*u,p);
         for(npy_intp offset=0,j=0;j<dim;j++,offset+=num_points_curve){
            npy_intp ptr=offset+i;
            curve[ptr]=p[j];
         }
      }

      return (PyObject*) curve_array;
   }

   fail:
      Py_XDECREF(u_array);
      Py_XDECREF(curve_array);
      Py_INCREF(Py_None);
      return Py_None;
}

static PyMethodDef _teardrop_methods[] = {
   {"_teardrop2d_parametric_coordinate_at_maximum_pool_width",  
     (PyCFunction)_teardrop2d_parametric_coordinate_at_maximum_pool_width,
     METH_VARARGS|METH_KEYWORDS,
     "Computes parametric coordinate (Right side of pool) at location of maximum pool width.\n\
      Args \n\
      -----\n\
      capsule: PyCapsule storing c++ pointer to Teardrop2D object\n\
      returns \n\
      -------\n\
      u: parametric coordinate at location of maximum pool width.\n"},
   {"_teardrop2d_tangent_to_curve",  
     (PyCFunction)_teardrop2d_tangent_to_curve,
     METH_VARARGS|METH_KEYWORDS,
     "Computes point on curve and tangent to curve for input parametric coordinate.\n\
      Args \n\
      -----\n\
      real u: Parametric coordinate in domain [0.0,1.0].\n\
      capsule: PyCapsule storing c++ pointer to Teardrop2D object\n\
      returns \n\
      -------\n\
      rho: spatial coordinates of curve at u, ndarrray(shape=(2,), dtype=numpy.float64)\n\
      rho_1: ndarrray(shape=(2,), dtype=numpy.float64) tangent to curve at u.\n"},
   {"_teardrop2d_predictor",  
     (PyCFunction)_teardrop2d_predictor,
     METH_VARARGS|METH_KEYWORDS,
     "Predicts closest point on curve; also computes tangent at estimated closest point..\n\
      Args \n\
      -----\n\
      real x: x coordinate of arbitrary point.\n\
      real x: y coordinate of arbitrary point.\n\
      capsule: PyCapsule storing c++ pointer to Teardrop2D object\n\
      returns \n\
      -------\n\
      rho: spatial coordinates of curve at estimated closed point, ndarrray(shape=(2,), dtype=numpy.float64)\n\
      rho_1: ndarrray(shape=(2,), dtype=numpy.float64) tangent to curve at estimated closest point.\n"},
   {"_get_teardrop2d_iterate",  
     (PyCFunction)_get_teardrop2d_iterate,
     METH_VARARGS|METH_KEYWORDS,
     "Computes point on curve, tangent to curve and vector from (x,y) to point on curve.\n\
      Args \n\
      -----\n\
      real x: x coordinate of arbitrary point.\n\
      real x: y coordinate of arbitrary point.\n\
      real u: [0.0,1.0] parametric point where curve is evaluated.\n\
      capsule: PyCapsule storing c++ pointer to Teardrop2D object\n\
      returns \n\
      -------\n\
      j11: real scale Newton Raphson jacobian used in closest point projection.\n\
      rho: ndarrray(shape=(2,), dtype=numpy.float64)\n\
      dr: (x,y)-rho: ndarrray(shape=(2,), dtype=numpy.float64)\n\
      rho_1: ndarrray(shape=(2,), dtype=numpy.float64) tangent to curve at 'u'\n\
      rho_11: ndarrray(shape=(2,), dtype=numpy.float64) derivative of rho_1 at 'u.'\n"},
   {"_get_beziercurve",  
     (PyCFunction)_get_beziercurve,
     METH_VARARGS|METH_KEYWORDS,
     "Creates object representing Bezier curve.\n\
      Args \n\
      -----\n\
      real control_points: Points used to define curve in 2 or 3 dimensions.\n\
      returns \n\
      -------\n\
      PyCapsule storing c++ pointer to BezierCurve object\n"},
   {"_get_teardrop2d",  
     (PyCFunction)_get_teardrop2d,
     METH_VARARGS|METH_KEYWORDS,
     "Creates object representing teardrop curve in 2D.\n\
      Args \n\
      -----\n\
      real control_points: Points used to define curve in 2 dimensions.\n\
      returns \n\
      -------\n\
      PyCapsule storing c++ pointer to Teardrop2D object\n"},
   {"_get_teardrop2d_cpp",  
     (PyCFunction)_get_teardrop2d_cpp,
     METH_VARARGS|METH_KEYWORDS,
     "Computes closest point on curve for input point (x,y).\n\
      Args \n\
      -----\n\
      real x: x coordinate of arbitrary point.\n\
      real x: y coordinate of arbitrary point.\n\
      returns \n\
      -------\n\
      u: double parametric coordinate for closest point on curve.\n\
      distance: double distance value to curve from input point (x,y).\n"},
   {"_get_teardrop2d_size",  
     (PyCFunction)_get_teardrop2d_size,
     METH_VARARGS|METH_KEYWORDS,
     "Computes pool (width,length).\n\
      Args \n\
      -----\n\
      python capsule object storing pointer to 'Teardrop2D' object.\n\
      returns \n\
      -------\n\
      Tuple storing following two double values:\n\
      width: value for pool width.\n\
      length: value for pool length.\n"},
   {"_get_teardrop2d_position",  
     (PyCFunction)_get_teardrop2d_position,
     METH_VARARGS|METH_KEYWORDS,
     "Computes and returns y-component of pool position.\n\
      Args \n\
      -----\n\
      python capsule object storing pointer to 'Teardrop2D' object.\n\
      returns \n\
      -------\n\
      double value for y-component of pool postion.\n"},
   {"_beziercurve_compute_curve",  
     (PyCFunction)_beziercurve_compute_curve,
     METH_VARARGS|METH_KEYWORDS,
     "Computes spatial coordinates of curve.\n\
      Args \n\
      -----\n\
      param u (type=numpy.nddarray(shape=num_points,): Curve is evaluated at given points in [0.0,1.0].\n\
      param PyCapsule: stores BezierCurve c++ allocated object which does computations.\n\
      returns \n\
      -------\n\
      numpy.ndarray(shape=(num_points_curve,dim)): Curve in 2D or 3D depending \n\
              upon original input construction control points.\n"},
   {"_teardrop2d_compute_curve",  
     (PyCFunction)_teardrop2d_compute_curve,
     METH_VARARGS|METH_KEYWORDS,
     "Computes spatial coordinates of curve.\n\
      Args \n\
      -----\n\
      param u (type=numpy.nddarray(shape=num_points,): Curve is evaluated at given points in [0.0,1.0].\n\
      param PyCapsule: stores Teardrop2D c++ allocated object which does computations.\n\
      returns \n\
      -------\n\
      numpy.ndarray(shape=(num_points_curve,dim)): Curve in 2D or 3D depending \n\
              upon original input Teardrop2D construction control points.\n"},
   {nullptr, nullptr, 0, nullptr}        /* Sentinel */
};

static struct PyModuleDef _teardropmodule = {
   PyModuleDef_HEAD_INIT,
   "_teardrop",
   NULL,
   -1,
   _teardrop_methods
};

PyMODINIT_FUNC
PyInit__teardrop(void)
{
   PyObject *m=NULL;
   m=PyModule_Create(&_teardropmodule);
	import_array();
   return m;
};

