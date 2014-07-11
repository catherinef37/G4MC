//HRS Exp, v 1.0, 2010/12/26 , By Jixie Zhang
// ********************************************************************
//
// $Id: HRSGeneralPhysics.hh,v 1.0, 2010/12/26   HRS Exp $
// --------------------------------------------------------------
//

#ifndef HRSGeneralPhysics_h
#define HRSGeneralPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4Decay.hh"

class HRSGeneralPhysics : public G4VPhysicsConstructor
{
  public:
    HRSGeneralPhysics(const G4String& name = "general");
    virtual ~HRSGeneralPhysics();

  public:
    // This method will be invoked in the Construct() method.
    // each particle type will be instantiated
    virtual void ConstructParticle();

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    virtual void ConstructProcess();

  protected:
    G4Decay fDecayProcess;
};


#endif








