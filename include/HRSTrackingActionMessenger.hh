//HRS Exp, v 1.0, 2010/12/26 , By Jixie Zhang
// ********************************************************************
//
// $Id: HRSTrackingActionMessenger.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
#ifndef HRSTrackingActionMessenger_h
#define HRSTrackingActionMessenger_h 1

class HRSTrackingAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"
#include "globals.hh"

class HRSTrackingActionMessenger: public G4UImessenger
{
  public:
    HRSTrackingActionMessenger(HRSTrackingAction* mtarget);
    ~HRSTrackingActionMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    HRSTrackingAction * target;

  private: //commands
    G4UIcmdWithAnInteger*	trackIdCmd;
    G4UIcmdWithAnInteger*	noSecondaryCmd;
    G4UIcmdWithAString*		particleCmd;
};

#endif


