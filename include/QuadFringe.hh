#ifndef QUADFRINGE_HH
#define QUADFRINGE_HH

#include "Python.h"
#include "G4ThreeVector.hh"
class QuadFringe{

 public: // with description
  
  QuadFringe();
  ~QuadFringe();
  G4ThreeVector GetFieldValue(double, double, double, double, double, double);
  
 private:
  PyObject *pName, *pModule, *pDict, *pFunc;
  PyObject *pArgs, *pValue;

};
#endif

