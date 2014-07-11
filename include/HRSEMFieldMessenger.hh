// ********************************************************************
//
// $Id: HRSEMFieldMessenger.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//............................................................................

#ifndef HRSEMFieldMessenger_h
#define HRSEMFieldMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4String.hh"

#include "HRSEMField.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

class HRSEMFieldMessenger: public G4UImessenger
{
public:
	HRSEMFieldMessenger(HRSEMField* pEMfield);
	~HRSEMFieldMessenger();

	void SetNewValue(G4UIcommand*, G4String);
	G4String GetCurrentValue(G4UIcommand* command);

private:
	HRSEMField*              target;

	G4UIcmdWithADoubleAndUnit* ErDCCmd;
	G4UIcmdWithADoubleAndUnit* ErInnerCmd;
	G4UIcmdWith3VectorAndUnit* BField3VCmd;
	G4UIcmdWith3VectorAndUnit* EField3VCmd;

};

#endif

