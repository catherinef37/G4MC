// ********************************************************************
//
// $Id: HRSSteppingActionMessenger.cc,v 1.0, 2010/12/26   HRS Exp $
// --------------------------------------------------------------
//

#include "HRSSteppingActionMessenger.hh"
#include "HRSSteppingAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"
#include "G4ios.hh"

HRSSteppingActionMessenger::HRSSteppingActionMessenger(HRSSteppingAction * mpga)
:target(mpga)
{
  stepDir = new G4UIdirectory("/stepping/");
  stepDir->SetGuidance("Jixie's stepping control commands");

  verboseCmd = new G4UIcmdWithAnInteger("/stepping/verbose",this);
  verboseCmd->SetGuidance("Verbose level for each Stepping.");
  verboseCmd->SetGuidance("Maximum number of step is 1024 and the track will be killed whenever it hits the virtual boundary: r>=116mm || fabs(z)>=151mm");
  verboseCmd->SetGuidance("Level 1: doesn't print VolumnName and ProcessName,");
  verboseCmd->SetGuidance("Level 2: if verboselevel>=2, print all, including VolumnName and ProcessName");
  verboseCmd->SetGuidance("Level 3: print all; track will be killed if r> 5000.0mm");
  verboseCmd->SetGuidance("Level 4: print all; track will be killed if r> 2500.0mm");
  verboseCmd->SetGuidance("Level 5: print all; track will be killed if r> 800.0mm");
  verboseCmd->SetGuidance("Level 6: print all; but only print steps inside those physicalVolumns in the printlist \n");
  verboseCmd->SetParameterName("stepverbose",false);
  verboseCmd->SetRange("stepverbose>=0");
  verboseCmd->SetDefaultValue(2);

  
  //cmd to add phys vol into print list
  emptyPrintListCmd = new G4UIcmdWithoutParameter("/stepping/emptyPrintList",this);
  emptyPrintListCmd->SetGuidance("Empty the print list"); 
  emptyPrintListCmd->SetGuidance("This cmd works only if stepping verbose==6"); 

  add2PrintListCmd = new G4UIcmdWithAString("/stepping/add2PrintList",this);
  add2PrintListCmd->SetGuidance("Add 1 physical volumn to the print list"); 
  add2PrintListCmd->SetGuidance("The cmd works only if stepping verbose==6"); 
  add2PrintListCmd->SetGuidance("To add all physical volumns to the print list, use \"all\" as argument."); 
  add2PrintListCmd->SetParameterName("physvol",false);

}

HRSSteppingActionMessenger::~HRSSteppingActionMessenger()
{
  delete stepDir;
  delete verboseCmd;
  delete emptyPrintListCmd;
  delete add2PrintListCmd;
}

void HRSSteppingActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==verboseCmd )
  { 
	  target->SetVerbose(verboseCmd->GetNewIntValue(newValue)); 
  }

  if( command==add2PrintListCmd )
  { 
	  target->Add2PrintList(newValue); 
  }
  
  if( command==emptyPrintListCmd )
  { 
	  target->EmptyPrintList(); 
  }

}

G4String HRSSteppingActionMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==verboseCmd )
  { cv = command->ConvertToString(target->GetVerbose()); }

  return cv;
}
