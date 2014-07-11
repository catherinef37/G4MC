//HRS Exp, v 1.0, 2010/12/26 , By Jixie Zhang
// ********************************************************************
//
// $Id: HRSEventActionMessenger.cc,v 1.0, 2010/12/26   HRS Exp $
// --------------------------------------------------------------
//

#include "HRSEventActionMessenger.hh"
#include "HRSEventAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"

HRSEventActionMessenger::HRSEventActionMessenger(HRSEventAction * mpga)
:target(mpga)
{
  verboseCmd = new G4UIcmdWithAnInteger("/mydet/verbose",this);
  verboseCmd->SetGuidance("Verbose level for event action. The following will be displayed:");
  verboseCmd->SetGuidance("Level 0: None.");
  verboseCmd->SetGuidance("Level 1: Print total number of hits in each fired detector.");
  verboseCmd->SetGuidance("Level 2: all the above and 1) at the beginning of event, local event id, ");
  verboseCmd->SetGuidance("         vertex and momenta of primary particles; 2) print each hits");
  verboseCmd->SetGuidance("Level 3: all the above and at the end of event, print the whole trajectory for each track.");
  verboseCmd->SetGuidance("Level 4: all the above and at the end of event, print PDGcode Momentum TrackID for each primary and its daughters.");
  verboseCmd->SetParameterName("eventlevel",true);
  verboseCmd->SetRange("eventlevel>=0");
  verboseCmd->SetDefaultValue(0);
}

HRSEventActionMessenger::~HRSEventActionMessenger()
{
  delete verboseCmd;
}

void HRSEventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==verboseCmd )
  { target->SetVerbose(verboseCmd->GetNewIntValue(newValue)); }
}

G4String HRSEventActionMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==verboseCmd )
  { cv = verboseCmd->ConvertToString(target->GetVerbose()); }

  return cv;
}
