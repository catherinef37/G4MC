//HRS Exp, v 1.0, 2010/12/26 , By Jixie Zhang
// ********************************************************************
//
// $Id: HRSEventActionMessenger.hh,v 1.0, 2010/12/26   HRS Exp $
// --------------------------------------------------------------
//
#ifndef HRSEventActionMessenger_h
#define HRSEventActionMessenger_h 1

class HRSEventAction;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"
#include "globals.hh"

class HRSEventActionMessenger: public G4UImessenger
{
public:
	HRSEventActionMessenger(HRSEventAction* mpga);
	~HRSEventActionMessenger();

public:
	void SetNewValue(G4UIcommand * command,G4String newValues);
	G4String GetCurrentValue(G4UIcommand * command);

private:
	HRSEventAction* target;

private: //commands
	G4UIcmdWithAnInteger*  verboseCmd;

};

#endif


