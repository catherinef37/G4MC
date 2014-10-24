// ********************************************************************
//
// $Id: HRSMuonPhysics.hh,v 1.0, 2010/12/26   HRS Exp $
// --------------------------------------------------------------

#ifndef HRSMuonPhysics_h
#define HRSMuonPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuIonisation.hh"
#include "G4hIonisation.hh"

class HRSMuonPhysics : public G4VPhysicsConstructor
{
  public:
    HRSMuonPhysics(const G4String& name="muon");
    virtual ~HRSMuonPhysics();

  public:
    // This method will be invoked in the Construct() method.
    // each particle type will be instantiated
    virtual void ConstructParticle();

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    virtual void ConstructProcess();

  protected:

};


#endif

