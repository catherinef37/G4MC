// ********************************************************************
//
// $Id: HRSPhysicsList.hh,v 1.0, 2010/12/26   HRS Exp $
// --------------------------------------------------------------

#ifndef HRSPhysicsList_h
#define HRSPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class HRSPhysicsList: public G4VModularPhysicsList
{
public:
  HRSPhysicsList();
  virtual ~HRSPhysicsList();

public:
  // SetCuts()
  virtual void SetCuts();


};


#endif



