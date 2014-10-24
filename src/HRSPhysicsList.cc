// ********************************************************************
//
// $Id: HRSPhysicsList.cc,v 1.0, 2010/12/26   HRS Exp $
// --------------------------------------------------------------
#include "UsageManager.hh"
#include "HRSPhysicsList.hh"

#include "HRSGeneralPhysics.hh"
#include "HRSEMPhysics.hh"
#include "HRSMuonPhysics.hh"
#include "HRSHadronPhysics.hh"
#include "HRSIonPhysics.hh"
#include "G4EmExtraPhysics.hh"


#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip>
#include "G4VUserPhysicsList.hh"

//http://geant4.cern.ch/support/proc_mod_catalog/physics_lists/referencePL.shtml

#include "G4QStoppingPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "HadronPhysicsQGSP_BERT.hh"
#include "HadronPhysicsQGSP_BERT_HP.hh"
#include "HadronPhysicsQGSP_BIC.hh"
#include "HadronPhysicsQGSP_BIC_HP.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "ElectroNuclearPhysics.hh"

HRSPhysicsList::HRSPhysicsList():  G4VModularPhysicsList()
{

	G4cout << "" << G4endl;
	G4cout << "You are using the HRSPhysicsList." << G4endl;
	G4cout << "Full set of particles (barions bosons and mesons) will be created and" << G4endl;
	G4cout << "Standard EM Physics and Low & High Energy parameterized models will be applied." << G4endl;
	G4cout << "HRSPhysicsList is optimized for robustness and not for any particular usage." << G4endl;
	G4cout << "For the hadronic physics, educated guesses of physics list are prepared for " << G4endl;
	G4cout << "various use cases. When you will start REAL calculations for your own interest," << G4endl;
	G4cout << "please consider the usage of hadronic_lists instead of HRSPhysicsLists." << G4endl;
	G4cout << "More information can also be found from the Geant4 HyperNews." << G4endl;
	G4cout << "http://geant4-hn.slac.stanford.edu:5090/Geant4-HyperNews/index" << G4endl;
	G4cout << "http://geant4.cern.ch/support/proc_mod_catalog/physics_lists/referencePL.shtml" << G4endl;
	G4cout << "" << G4endl;

	defaultCutValue = 0.7*mm; //by jixie
	SetVerboseLevel(1);  //LEVEL 1 WILL PRINT WARNINGS

	// in Geant4.9.x, G4VUserPhysicsList::AddTransportation(); has been called 
	// in G4VModularPhysicsList::ConstructProcess(), no need to invoke it again here
	//G4VUserPhysicsList::AddTransportation();
	
	// General Physics
	RegisterPhysics( new HRSGeneralPhysics("general") );

	// Generate physics, also create all Particles
	//RegisterPhysics(  new G4DecayPhysics("decays"));

	// EM Physics
	RegisterPhysics( new HRSEMPhysics("standard EM"));

	// EM physics
	//RegisterPhysics( new G4EmStandardPhysics());

	// Electro-Nuclear Physics
	RegisterPhysics(  new ElectroNuclearPhysics("Electro-Nuclear"));


	// Muon Physics
	RegisterPhysics(  new HRSMuonPhysics("muon"));

	int verboseLevel=1;
	//By Jixie: G4EmExtraPhysics is not working in this file, it caused a segmentation fault
	//if used with my EM physics and general physics
	//I need to figure it out later
	//Without this one, eletron-nucleon interaction is off, then no neutron will be created
	//RegisterPhysics( new G4EmExtraPhysics(verboseLevel));


	// Ion Physics
	RegisterPhysics( new HRSIonPhysics("ion"));

	// Hadron Physics
	extern UsageManager* gConfig;
	G4String pPhysicsModel=gConfig->GetArgument("PhysicsModel");
	if(pPhysicsModel=="QGSP_BERT")
	{
		RegisterPhysics( new G4QStoppingPhysics(verboseLevel));
		RegisterPhysics( new G4HadronElasticPhysics(verboseLevel) );
		RegisterPhysics( new HadronPhysicsQGSP_BERT("hadron"));
		//This cut will force neutron to die in 10 micro seconds
		RegisterPhysics( new G4NeutronTrackingCut(verboseLevel));
	}
	else if(pPhysicsModel=="QGSP_BERT_HP")
	{
		RegisterPhysics( new G4QStoppingPhysics(verboseLevel));
		RegisterPhysics( new G4HadronElasticPhysicsHP(verboseLevel) );
		RegisterPhysics( new HadronPhysicsQGSP_BERT_HP("hadron"));
		RegisterPhysics( new G4NeutronTrackingCut(verboseLevel));
	}
	else if(pPhysicsModel=="QGSP_BIC")
	{
		RegisterPhysics( new G4QStoppingPhysics(verboseLevel));
		RegisterPhysics( new G4HadronElasticPhysics(verboseLevel) );
		RegisterPhysics( new HadronPhysicsQGSP_BIC("hadron"));
		RegisterPhysics( new G4NeutronTrackingCut(verboseLevel));
	}
	else if(pPhysicsModel=="QGSP_BIC_HP")
	{
		RegisterPhysics( new G4QStoppingPhysics(verboseLevel));
		RegisterPhysics( new G4HadronElasticPhysicsHP(verboseLevel) );
		RegisterPhysics( new HadronPhysicsQGSP_BIC_HP("hadron"));
		RegisterPhysics( new G4NeutronTrackingCut(verboseLevel));
	}
	else
	{
		RegisterPhysics( new HRSHadronPhysics("hadron"));
	}
}

HRSPhysicsList::~HRSPhysicsList()
{
}

void HRSPhysicsList::SetCuts()
{
	//  " G4VUserPhysicsList::SetCutsWithDefault" method sets
	//   the default cut value for all particle types
	SetCutsWithDefault();
}

