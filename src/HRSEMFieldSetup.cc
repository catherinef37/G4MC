// ********************************************************************
//
// $Id: HRSEMField.hh,v 1.0, 2010/12/26 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//
//   User Field class Setup implementation.
//
#include "HRSEMFieldSetup.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4UniformMagField.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ios.hh"

//////////////////////////////////////////////////////////////////////////
HRSEMFieldSetup* HRSEMFieldSetup::fHRSEMFieldSetup=0;
HRSEMFieldSetup* HRSEMFieldSetup::GetHRSEMFieldSetup()
{ 
	if(!fHRSEMFieldSetup)  
	{
		G4cout<<"HRSEMFieldSetup is not initialized yet...exit...\n";
		exit(-99);
	}
	return fHRSEMFieldSetup; 
}

//////////////////////////////////////////////////////////////////////////
//
HRSEMFieldSetup::HRSEMFieldSetup()
: fChordFinder(0), fStepper(0), fIntgrDriver(0)
{
	fHRSEMFieldSetup=this;

	//global EM field
	fEMfield = new HRSEMField();
	messenger = new HRSEMFieldSetupMessenger(this) ;
	fEquation = new G4EqMagElectricField(fEMfield);
	fMinStep  = 0.001*mm ; // minimal step of 1 miron, default is 0.01 mm
	fStepperType = 4 ;     // ClassicalRK4 -- the default stepper

	fFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
	fChordFinder = 0;   //will be set in UpdateField()
	UpdateField();


	//Local field  FZB1
	fMagFieldFZB1 = new G4UniformMagField(G4ThreeVector(0,0,0));
	fEquationFZB1 = new G4Mag_UsualEqRhs(fMagFieldFZB1);	
	fStepperFZB1 = new G4ClassicalRK4(fEquationFZB1);
	fLocalFieldManagerFZB1 = new G4FieldManager();
	fChordFinderFZB1 = 0;
	UpdateFieldFZB1();


	//Local field  FZB2
	fMagFieldFZB2 = new G4UniformMagField(G4ThreeVector(0,0,0));
	fEquationFZB2 = new G4Mag_UsualEqRhs(fMagFieldFZB2);	
	fStepperFZB2 = new G4ClassicalRK4(fEquationFZB2);
	fLocalFieldManagerFZB2 = new G4FieldManager();
	fChordFinderFZB2 = 0;
	UpdateFieldFZB2();

}

/////////////////////////////////////////////////////////////////////////////////
//

HRSEMFieldSetup::~HRSEMFieldSetup()
{
	if(fChordFinder) delete fChordFinder;
	if(fStepper)     delete fStepper;
	if(fEquation)    delete fEquation;
	if(fEMfield)     delete fEMfield;
	if(messenger)    delete messenger;
}

/////////////////////////////////////////////////////////////////////////////
//
// Register this field to 'global' Field Manager and
// Create Stepper and Chord Finder with predefined type, minstep (resp.)
//

void HRSEMFieldSetup::UpdateField()
{
	SetStepper();
	G4cout<<"HRSEMFieldSetup:: The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fFieldManager->SetDetectorField(fEMfield);

	if(fChordFinder) delete fChordFinder;
	fIntgrDriver = new G4MagInt_Driver(fMinStep,fStepper,fStepper->GetNumberOfVariables());
	fChordFinder = new G4ChordFinder(fIntgrDriver);
	fFieldManager->SetChordFinder( fChordFinder );
}


/////////////////////////////////////////////////////////////////////////////
void HRSEMFieldSetup::UpdateFieldFZB1()
{
	G4cout<<"HRSEMFieldSetup:: The minimal step for FZB1 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZB1->SetDetectorField(fMagFieldFZB1);

	if(fChordFinderFZB1) delete fChordFinderFZB1;
	fChordFinderFZB1 = new G4ChordFinder((G4MagneticField*) fMagFieldFZB1, fMinStep, fStepperFZB1);
	fLocalFieldManagerFZB1->SetChordFinder( fChordFinderFZB1 );
}

/////////////////////////////////////////////////////////////////////////////
void HRSEMFieldSetup::UpdateFieldFZB2()
{
	G4cout<<"HRSEMFieldSetup:: The minimal step for FZB2 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZB2->SetDetectorField(fMagFieldFZB2);

	if(fChordFinderFZB2) delete fChordFinderFZB2;
	fChordFinderFZB2 = new G4ChordFinder((G4MagneticField*) fMagFieldFZB2, fMinStep, fStepperFZB2);
	fLocalFieldManagerFZB2->SetChordFinder( fChordFinderFZB2 );
}

/////////////////////////////////////////////////////////////////////////////
void HRSEMFieldSetup::SetBField3VFZB1(G4ThreeVector fieldVector)
{
	if(fMagFieldFZB1) delete fMagFieldFZB1;

	if(fieldVector != G4ThreeVector(0.,0.,0.))
	{ 
		fMagFieldFZB1 = new  G4UniformMagField(fieldVector);
	}
	else 
	{
		// If the new field's value is Zero, then setting the pointer to zero ensures 
		// that it is not used for propagation.
		fMagFieldFZB1 = 0; 
	}
	//call UpdateFieldFZB1() to do the update or the next 2 lines can do the same job 
	fLocalFieldManagerFZB1->SetDetectorField(fMagFieldFZB1);
	fEquationFZB1->SetFieldObj( fMagFieldFZB1 );
	//UpdateFieldFZB1();   //No need this  line if the above 2 lines used
}


/////////////////////////////////////////////////////////////////////////////
void HRSEMFieldSetup::SetBField3VFZB2(G4ThreeVector fieldVector)
{
	if(fMagFieldFZB2) delete fMagFieldFZB2;

	if(fieldVector != G4ThreeVector(0.,0.,0.))
	{ 
		fMagFieldFZB2 = new  G4UniformMagField(fieldVector);
	}
	else 
	{
		// If the new field's value is Zero, then setting the pointer to zero ensures 
		// that it is not used for propagation.
		fMagFieldFZB2 = 0; 
	}
	//call UpdateFieldFZB2() to do the update or the next 2 lines can do the same job 
	fLocalFieldManagerFZB2->SetDetectorField(fMagFieldFZB2);
	fEquationFZB2->SetFieldObj( fMagFieldFZB2 );
	//UpdateFieldFZB2();   //No need this  line if the above 2 lines used
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

void HRSEMFieldSetup::SetStepper()
{
	G4int nvar = 8;

	if(fStepper) delete fStepper;

	switch ( fStepperType )
	{
	case 0:
		fStepper = new G4ExplicitEuler( fEquation, nvar );
		G4cout<<"HRSEMFieldSetup:: G4ExplicitEuler is calledS"<<G4endl;
		break;
	case 1:
		fStepper = new G4ImplicitEuler( fEquation, nvar );
		G4cout<<"HRSEMFieldSetup:: G4ImplicitEuler is called"<<G4endl;
		break;
	case 2:
		fStepper = new G4SimpleRunge( fEquation, nvar );
		G4cout<<"HRSEMFieldSetup:: G4SimpleRunge is called"<<G4endl;
		break;
	case 3:
		fStepper = new G4SimpleHeum( fEquation, nvar );
		G4cout<<"HRSEMFieldSetup:: G4SimpleHeum is called"<<G4endl;
		break;
	case 4:
		fStepper = new G4ClassicalRK4( fEquation, nvar );
		G4cout<<"HRSEMFieldSetup:: G4ClassicalRK4 (default) is called"<<G4endl;
		break;
	case 5:
		fStepper = new G4CashKarpRKF45( fEquation, nvar );
		G4cout<<"HRSEMFieldSetup:: G4CashKarpRKF45 is called"<<G4endl;
		break;

	//The following not working for electric field
	case 6:
		fStepper = 0 ; // new G4RKG3_Stepper( fMagEquation );
		G4cout<<"HRSEMFieldSetup:: G4RKG3_Stepper is not currently working for Electric Field"<<G4endl;
		break;
	case 7:
		fStepper = 0 ; // new G4HelixExplicitEuler( fMagEquation );
		G4cout<<"HRSEMFieldSetup:: G4HelixExplicitEuler is not valid for Electric Field"<<G4endl;
		break;
	case 8:
		fStepper = 0 ; //  new G4HelixImplicitEuler( fMagEquation );
		G4cout<<"HRSEMFieldSetup:: G4HelixImplicitEuler is not valid for Electric Field"<<G4endl;
		break;
	case 9:
		fStepper = 0 ; //  new G4HelixSimpleRunge( fMagEquation );
		G4cout<<"HRSEMFieldSetup:: G4HelixSimpleRunge is not valid for Electric Field"<<G4endl;
		break;
	default: fStepper = 0;
	}
}


///////////////////////////////////////////////////////////////////////////////
