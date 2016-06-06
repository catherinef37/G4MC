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
#include "BField_Quad.hh"
#include "BField_Dipole_Fringe.hh"
#include "BField_Dipole.hh"

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
	void UpdateFieldFZBL1();
	void SetBField3VFZBL1(G4double fieldGradient);
	G4FieldManager* GetFieldManagerFZBL1(){return fLocalFieldManagerFZBL1;}
	void UpdateFieldFZBR1();
	void SetBField3VFZBR1(G4double fieldGradient);
	G4FieldManager* GetFieldManagerFZBR1(){return fLocalFieldManagerFZBR1;}
	
	//Local field  FZB2
	void UpdateFieldFZBL2();
	void SetBField3VFZBL2(G4double fieldGradient);
	G4FieldManager* GetFieldManagerFZBL2(){return fLocalFieldManagerFZBL2;}
	void UpdateFieldFZBR2();
	void SetBField3VFZBR2(G4double fieldGradient);
	G4FieldManager* GetFieldManagerFZBR2(){return fLocalFieldManagerFZBR2;}

	//Local field  FZB3	
        void UpdateFieldFZBL3();
	void SetBField3VFZBL3(G4double fbend);
	G4FieldManager* GetFieldManagerFZBL3(){return fLocalFieldManagerFZBL3;}
        void UpdateFieldFZBR3();
	void SetBField3VFZBR3(G4double fbend);
	G4FieldManager* GetFieldManagerFZBR3(){return fLocalFieldManagerFZBR3;}

	//Local field  FZB4
	void UpdateFieldFZBL4();
	void SetBField3VFZBL4(G4double fieldGradient);
	G4FieldManager* GetFieldManagerFZBL4(){return fLocalFieldManagerFZBL4;}
	void UpdateFieldFZBR4();
	void SetBField3VFZBR4(G4double fieldGradient);
	G4FieldManager* GetFieldManagerFZBR4(){return fLocalFieldManagerFZBR4;}

private:
  HRSEMFieldSetupMessenger*   messenger;
  HRSEMField*                 fEMfield; 
  G4FieldManager*             fFieldManager;
  G4ChordFinder*              fChordFinder ;
  G4EqMagElectricField*       fEquation ;
  G4MagIntegratorStepper*     fStepper ;
  G4MagInt_Driver*            fIntgrDriver;
  G4MagInt_Driver*            fIntgrDriverFZBL1;
  G4MagInt_Driver*            fIntgrDriverFZBR1;
  G4MagInt_Driver*            fIntgrDriverFZBL2;
  G4MagInt_Driver*            fIntgrDriverFZBR2;
  G4MagInt_Driver*            fIntgrDriverFZBL3;
  G4MagInt_Driver*            fIntgrDriverFZBR3;
  G4MagInt_Driver*            fIntgrDriverFZBL4;
  G4MagInt_Driver*            fIntgrDriverFZBR4;
  
  G4int                       fStepperType ;
  G4double                    fMinStep ;
  
  G4int                       mSnakeModel;
  G4double                    mLHRSMomentum;
  G4double                    mRHRSMomentum;
  G4double                    mLHRSAngle;
  G4double                    mRHRSAngle;
  G4double                    mLSeptumAngle;
  G4double                    mRSeptumAngle;

  //for local field at FZB1 and FZB2
  BField_Quad*                fMagFieldFZBL1 ;
  G4Mag_UsualEqRhs*           fEquationFZBL1 ;
  G4ChordFinder*              fChordFinderFZBL1 ;
  G4MagIntegratorStepper*     fStepperFZBL1 ;
  G4FieldManager*             fLocalFieldManagerFZBL1;
  BField_Quad*                fMagFieldFZBR1 ;
  G4Mag_UsualEqRhs*           fEquationFZBR1 ;
  G4ChordFinder*              fChordFinderFZBR1 ;
  G4MagIntegratorStepper*     fStepperFZBR1 ;
  G4FieldManager*             fLocalFieldManagerFZBR1;
  
  BField_Quad*                fMagFieldFZBL2 ;
  G4Mag_UsualEqRhs*           fEquationFZBL2 ;
  G4ChordFinder*              fChordFinderFZBL2 ;
  G4MagIntegratorStepper*     fStepperFZBL2 ;
  G4FieldManager*             fLocalFieldManagerFZBL2;
  BField_Quad*                fMagFieldFZBR2 ;
  G4Mag_UsualEqRhs*           fEquationFZBR2 ;
  G4ChordFinder*              fChordFinderFZBR2 ;
  G4MagIntegratorStepper*     fStepperFZBR2 ;
  G4FieldManager*             fLocalFieldManagerFZBR2;
    
  BField_Dipole*              fMagFieldFZBL3 ;
  G4Mag_UsualEqRhs*           fEquationFZBL3 ;
  G4ChordFinder*              fChordFinderFZBL3 ;
  G4MagIntegratorStepper*     fStepperFZBL3 ;
  G4FieldManager*             fLocalFieldManagerFZBL3;
  BField_Dipole*              fMagFieldFZBR3 ;
  G4Mag_UsualEqRhs*           fEquationFZBR3 ;
  G4ChordFinder*              fChordFinderFZBR3 ;
  G4MagIntegratorStepper*     fStepperFZBR3 ;
  G4FieldManager*             fLocalFieldManagerFZBR3;
  
  BField_Quad*                fMagFieldFZBL4 ;
  G4Mag_UsualEqRhs*           fEquationFZBL4 ;
  G4ChordFinder*              fChordFinderFZBL4 ;
  G4MagIntegratorStepper*     fStepperFZBL4 ;
  G4FieldManager*             fLocalFieldManagerFZBL4;
  BField_Quad*                fMagFieldFZBR4 ;
  G4Mag_UsualEqRhs*           fEquationFZBR4 ;
  G4ChordFinder*              fChordFinderFZBR4 ;
  G4MagIntegratorStepper*     fStepperFZBR4 ;
  G4FieldManager*             fLocalFieldManagerFZBR4;
  
};

#endif
