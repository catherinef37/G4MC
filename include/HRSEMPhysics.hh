// ********************************************************************
//
// $Id: HRSEMPhysics.hh,v 1.0, 2010/12/26   HRS Exp $
// --------------------------------------------------------------

#ifndef HRSEMPhysics_h
#define HRSEMPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4CoulombScattering.hh"
#include "G4eplusAnnihilation.hh"

class HRSEMPhysics : public G4VPhysicsConstructor
{
  public:
    HRSEMPhysics(const G4String& name ="EM");
    virtual ~HRSEMPhysics();

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
