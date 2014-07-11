// ********************************************************************
//
// $Id: HRSEMFieldSetup.hh,v 1.0, 2010/12/26   HRS Exp $
//
//    A class for control of the Electromagnetic Field of the detector.
//  The field for this case is reading from class HRSEMFiled.
//
//  It is simply a 'setup' class that creates the field and other necessary parts
// ********************************************************************

#ifndef HRSEMFieldSetup_H
#define HRSEMFieldSetup_H 1

#include "HRSEMField.hh"
#include "HRSEMFieldSetupMessenger.hh"

class G4FieldManager;
class G4ChordFinder;
class G4EquationOfMotion;
class G4Mag_UsualEqRhs;
class G4Mag_EqRhs;
class G4EqMagElectricField;
class G4MagIntegratorStepper;
class G4MagInt_Driver; 
class G4UniformMagField;

class HRSEMFieldSetup 
{
public:
	//Static method which returns the singleton pointer of this class.
	static HRSEMFieldSetup* GetHRSEMFieldSetup();

private:
	static HRSEMFieldSetup* fHRSEMFieldSetup;

public: 
	HRSEMFieldSetup() ;         
	~HRSEMFieldSetup() ;  

	void SetStepper();
	void UpdateField();

	inline  void SetStepperType( G4int val) { fStepperType = val ; }
	inline  G4int GetStepperType() {return fStepperType; }

	inline void SetMinStep(G4double val) { fMinStep = val ; }
	inline G4double GetMinStep() { return fMinStep ; }
	
	G4FieldManager* GetFieldManager(){return fFieldManager;}

	//Local field  FZB1
	void UpdateFieldFZB1();
	void SetBField3VFZB1(G4ThreeVector fieldVector);
	G4FieldManager* GetFieldManagerFZB1(){return fLocalFieldManagerFZB1;}
	
	//Local field  FZB2
	void UpdateFieldFZB2();
	void SetBField3VFZB2(G4ThreeVector fieldVector);
	G4FieldManager* GetFieldManagerFZB2(){return fLocalFieldManagerFZB2;}

private:
	HRSEMFieldSetupMessenger*   messenger;
	HRSEMField*                 fEMfield; 
	G4FieldManager*             fFieldManager;
	G4ChordFinder*              fChordFinder ;
	G4EqMagElectricField*       fEquation ;
	G4MagIntegratorStepper*     fStepper ;
	G4MagInt_Driver*            fIntgrDriver;

	G4int                       fStepperType ;
	G4double                    fMinStep ;

	//for local field at FZB1 and FZB2
	G4UniformMagField*          fMagFieldFZB1 ;
	G4Mag_UsualEqRhs*           fEquationFZB1 ;
	G4ChordFinder*              fChordFinderFZB1 ;
	G4MagIntegratorStepper*     fStepperFZB1 ;
	G4FieldManager*             fLocalFieldManagerFZB1;

	G4UniformMagField*          fMagFieldFZB2 ;
	G4Mag_UsualEqRhs*           fEquationFZB2 ;
	G4ChordFinder*              fChordFinderFZB2 ;
	G4MagIntegratorStepper*     fStepperFZB2 ;
	G4FieldManager*             fLocalFieldManagerFZB2;
};

#endif
