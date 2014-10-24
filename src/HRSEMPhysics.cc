// ********************************************************************
//
// $Id: HRSEMPhysics.cc,v 1.0, 2010/12/26  HRS Exp $
// --------------------------------------------------------------

#include "HRSEMPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>

//#define EMPHYSICS_NO_BREM true

HRSEMPhysics::HRSEMPhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{
}

HRSEMPhysics::~HRSEMPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"


// e- e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"


void HRSEMPhysics::ConstructParticle()
{
  // gamma
  G4Gamma::GammaDefinition();

  // electron
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}


#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"

void HRSEMPhysics::ConstructProcess()
{
   G4ProcessManager * pManager = 0;
 
   //Gamma
   pManager = G4Gamma::Gamma()->GetProcessManager();
   pManager->AddDiscreteProcess(new G4GammaConversion());
   pManager->AddDiscreteProcess(new G4ComptonScattering());
   pManager->AddDiscreteProcess(new G4PhotoElectricEffect());  


   //Electorn
   pManager = G4Electron::Electron()->GetProcessManager();
   int pUseMultipleScatering=1;
   if(pUseMultipleScatering)
   {
	   G4eMultipleScattering* msc = new G4eMultipleScattering();
	   //msc->AddEmModel(0, new G4UrbanMscModel93());
	   msc->SetStepLimitType(fUseDistanceToBoundary);
	   pManager->AddProcess(msc,                   -1, 1, 1);
   }
   else
   {
	   //test if this single scatering will return none zero NIEL
	   //conclusion: the NIEL always zero!! I looked into the source code,
	   //the value is not set by any physics process 
	   pManager->AddProcess(new G4CoulombScattering());
   }
   G4eIonisation* eIoni = new G4eIonisation();
   eIoni->SetStepFunction(0.2, 100*um);      
   pManager->AddProcess(eIoni,                 -1, 2, 2);

   int pUseBrem=1;
   if(pUseBrem)
   {
	   pManager->AddProcess(new G4eBremsstrahlung, -1,3, 3);
   }
   //add step limiter, by jixie for Geant4.7.0 version and up
   pManager->AddProcess(new G4StepLimiter(),      -1, -1,4);
	

   //Positron
   pManager = G4Positron::Positron()->GetProcessManager();
   G4VProcess* theeplusMultipleScattering = new G4eMultipleScattering();
   G4VProcess* theeplusIonisation         = new G4eIonisation();
   G4VProcess* theeplusBremsstrahlung     = new G4eBremsstrahlung();
   G4VProcess* theeplusAnnihilation       = new G4eplusAnnihilation();
  
   pManager->AddProcess(theeplusMultipleScattering);
   pManager->AddProcess(theeplusIonisation);
   pManager->AddProcess(theeplusBremsstrahlung);
   pManager->AddProcess(theeplusAnnihilation);

   //
   // set ordering for AtRestDoIt
   pManager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
   //
   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(theeplusIonisation,         idxAlongStep,2);
   pManager->SetProcessOrdering(theeplusBremsstrahlung,     idxAlongStep,3);

   //
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(theeplusIonisation,         idxPostStep,2);
   pManager->SetProcessOrdering(theeplusBremsstrahlung,     idxPostStep,3);
   pManager->SetProcessOrdering(theeplusAnnihilation,       idxPostStep,4);

}
