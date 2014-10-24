//HRS Exp, v 1.0, 2010/12/26 , By Jixie Zhang
// ********************************************************************
//
//
// $Id: HRSSteppingVerbose.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
//
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose.
//   It shows how to extract informations during the tracking of a particle.
//
//..............................................................................
//..............................................................................

class HRSSteppingVerbose;

#ifndef HRSSteppingVerbose_h
#define HRSSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

//..............................................................................

class HRSSteppingVerbose : public G4SteppingVerbose 
{
 public:
   
  HRSSteppingVerbose();
 ~HRSSteppingVerbose();

  void StepInfo();
  void TrackingStarted();

};

//..............................................................................

#endif
