//HRS Exp, v 1.0, 2010/12/26 , By Jixie Zhang
// ********************************************************************
//
// $Id: HRSEMFieldSetupMessenger.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//............................................................................

#ifndef HRSEMFieldSetupMessenger_h
#define HRSEMFieldSetupMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3VectorAndUnit;

class HRSEMFieldSetup;

class HRSEMFieldSetupMessenger: public G4UImessenger
{
public:
	HRSEMFieldSetupMessenger(HRSEMFieldSetup* pEM);
	~HRSEMFieldSetupMessenger();

	void SetNewValue(G4UIcommand*, G4String );

	//private:

	HRSEMFieldSetup*           fEMFieldSetup;
	G4UIdirectory*             HRSFieldDir;
	G4UIcmdWithAnInteger*      StepperCmd;
	G4UIcmdWithADoubleAndUnit* MinStepCmd;
	G4UIcmdWithoutParameter*   UpdateCmd;
	G4UIcmdWith3VectorAndUnit *BField3VFZB1Cmd;
	G4UIcmdWith3VectorAndUnit *BField3VFZB2Cmd;
	G4UIcmdWith3VectorAndUnit *BField3VFZB3Cmd;
  G4UIcmdWith3VectorAndUnit *BField3VFZB3enCmd;
	G4UIcmdWith3VectorAndUnit *BField3VFZB4Cmd;

};

#endif
