#include "Python.h"
#include "QuadFringe.hh"
#include "G4ThreeVector.hh"
#include <iostream>
using namespace std;

QuadFringe::QuadFringe()
{
  //cout << "Initializing QuadFringe." << endl;
  Py_Initialize();
  pName = PyUnicode_FromString("fringe");
  pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  if (pModule != NULL) {
    pFunc = PyObject_GetAttrString(pModule, "fringe");
    Py_XDECREF(pFunc);
    Py_DECREF(pModule);
  }
  else {
    PyErr_Print();
    fprintf(stderr, "Failed to load \"%s\"\n", "fringe");
  }

}

/////////////////////////////////////////////////////////////////////////

QuadFringe::~QuadFringe()
{

  Py_Finalize();
  //cout << "Deleting QuadFringe!" << endl;

}

////////////////////////////////////////////////////////////////////////
// 
// 

G4ThreeVector QuadFringe::GetFieldValue(double x, double y, double z, double maxfield, double length, double radius){
  //cout << "Getting field value." << endl;
  int   i = 0;
  int   myargc = 5;
  float default_radius = 0.15;
  float radius_scale   = default_radius / radius;
  double myargv[5] = {x * radius_scale, y * radius_scale, z, maxfield, length};
  
  //cout << "Determining if function is callable." << endl;
  //cout << pFunc << endl;
  if (pFunc && PyCallable_Check(pFunc)) {
    //cout << "Making array of arguements." << endl;
    pArgs = PyTuple_New(myargc);
    for (i = 0; i < myargc; ++i) {
      pValue = PyFloat_FromDouble(myargv[i]);
      if (!pValue) {
	Py_DECREF(pArgs);
	Py_DECREF(pModule);
	fprintf(stderr, "Cannot convert argument\n");
      }
      /* pValue reference stolen here: */
      PyTuple_SetItem(pArgs, i, pValue);
    }
    //cout << "Making the call." << endl;
    pValue = PyObject_CallObject(pFunc, pArgs);
    Py_DECREF(pArgs);
    if (pValue != NULL) {
      //cout << PyFloat_AsDouble(pValue) << endl;                                                                                                                                                                
      //cout << PyTuple_GET_SIZE(pValue) << endl;
      PyObject* Bx = PyTuple_GetItem(pValue, 0);
      PyObject* By = PyTuple_GetItem(pValue, 1);
      PyObject* Bz = PyTuple_GetItem(pValue, 2);
      G4ThreeVector bfield(PyFloat_AsDouble(Bx) * tesla , PyFloat_AsDouble(By) * tesla , PyFloat_AsDouble(Bz) * tesla ); //turn on for fringe fields
      //G4ThreeVector bfield( maxfield / radius * y , maxfield / radius * x , 0.0); //turn on for no fringe fields
            
      //cout << PyFloat_AsDouble(Bx) << " ";
      //cout << PyFloat_AsDouble(By) << " ";
      //cout << PyFloat_AsDouble(Bz) << endl;
      Py_DECREF(pValue);
      return bfield;
    }
    else {
      Py_DECREF(pFunc);
      Py_DECREF(pModule);
      PyErr_Print();
      fprintf(stderr,"Call failed\n");
    }
  }
  else {
    if (PyErr_Occurred())
      PyErr_Print();
    fprintf(stderr, "Cannot find function \"%s\"\n", "fringe");
  }
  
}
