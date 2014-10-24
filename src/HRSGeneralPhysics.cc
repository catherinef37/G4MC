// ********************************************************************
//
// $Id: HRSGeneralPhysics.cc,v 1.0, 2010/12/26   HRS Exp $
// --------------------------------------------------------------

#include "HRSGeneralPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>

HRSGeneralPhysics::HRSGeneralPhysics(const G4String& name)
:  G4VPhysicsConstructor(name)
{
}

HRSGeneralPhysics::~HRSGeneralPhysics()
{
}

#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

void HRSGeneralPhysics::ConstructParticle()
{
	// In Alphabetical Order 

	//  Construct all barions
	G4BaryonConstructor* baryonConstructor = new G4BaryonConstructor();
	baryonConstructor -> ConstructParticle();
	delete baryonConstructor;

	// Construct all bosons (including geantinos)
	G4BosonConstructor* bosonConstructor = new G4BosonConstructor();
	bosonConstructor -> ConstructParticle();
	delete bosonConstructor;

	// Construct all ions 
	G4IonConstructor* ionConstructor = new G4IonConstructor();
	ionConstructor -> ConstructParticle();
	delete ionConstructor;

	// Construct all leptons 
	G4LeptonConstructor* leptonConstructor = new G4LeptonConstructor();
	leptonConstructor -> ConstructParticle();
	delete leptonConstructor;

	// Construct all mesons
	G4MesonConstructor* mesonConstructor = new G4MesonConstructor();
	mesonConstructor -> ConstructParticle();
	delete mesonConstructor;

	//  Construct  resonaces and quarks
	G4ShortLivedConstructor* shortLivedConstructor = new G4ShortLivedConstructor();
	shortLivedConstructor -> ConstructParticle();
	delete shortLivedConstructor;

}

#include "G4Decay.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

void HRSGeneralPhysics::ConstructProcess()
{
	// Add Decay Process
	G4Decay* theDecayProcess = new G4Decay();  
	theParticleIterator->reset();
	while( (*theParticleIterator)() ){
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		if (theDecayProcess->IsApplicable(*particle)) {
			pmanager ->AddProcess(theDecayProcess);
			// set ordering for PostStepDoIt and AtRestDoIt
			pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
			pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
		}
	}
}


