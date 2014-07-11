//HRS Exp, v 1.0, 2010/12/26 , By Jixie Zhang
// ********************************************************************
//
// $Id: HRSSteppingActionMessenger.hh,v 1.0, 2010/12/26   HRS Exp $
// --------------------------------------------------------------
//
#ifndef HRSSteppingActionMessenger_h
#define HRSSteppingActionMessenger_h 1

class G4UIdirectory;
class HRSSteppingAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

#include "G4UImessenger.hh"
#include "globals.hh"

class HRSSteppingActionMessenger: public G4UImessenger
{
  public:
    HRSSteppingActionMessenger(HRSSteppingAction* mpga);
    ~HRSSteppingActionMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    HRSSteppingAction* target;

  private: //commands
	G4UIdirectory*             stepDir;
    G4UIcmdWithAnInteger*      verboseCmd;
	G4UIcmdWithAString*        add2PrintListCmd;
	G4UIcmdWithoutParameter*   emptyPrintListCmd;
};

#endif


