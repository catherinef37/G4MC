//HRS Exp, v 1.0, 2010/12/26 , By Jixie Zhang
// ********************************************************************
//
// $Id: HRSEMFieldSetupMessenger.cc,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
// 

#include "HRSEMFieldSetup.hh"
#include "HRSEMFieldSetupMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
//////////////////////////////////////////////////////////////////////////////

HRSEMFieldSetupMessenger::HRSEMFieldSetupMessenger(HRSEMFieldSetup* pmsg)
  :fEMFieldSetup(pmsg)
{ 
  HRSFieldDir = new G4UIdirectory("/field/");
  HRSFieldDir->SetGuidance("HRS field tracking control. You can change StepperType and minimum step size");

  UpdateCmd = new G4UIcmdWithoutParameter("/field/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);

  StepperCmd = new G4UIcmdWithAnInteger("/field/setStepperType",this);
  StepperCmd->SetGuidance("Select stepper type for electric field");
  StepperCmd->SetParameterName("choice",true);
  StepperCmd->SetDefaultValue(4);
  StepperCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MinStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setMinStep",this);  
  MinStepCmd->SetGuidance("Define minimal step");
  MinStepCmd->SetGuidance("Magnetic field will be in Z direction.");
  MinStepCmd->SetParameterName("min step",false,false);
  MinStepCmd->SetDefaultUnit("mm");
  MinStepCmd->AvailableForStates(G4State_Idle);  
  
  BFieldFZBL1Cmd = new G4UIcmdWithADoubleAndUnit("/field/setFZBL1Field",this);  
  BFieldFZBL1Cmd->SetGuidance("Set the quad    magnetic field for FZBL1");
  BFieldFZBL1Cmd->SetParameterName("Field",false);
  BFieldFZBL1Cmd->SetDefaultUnit("tesla");  
  BFieldFZBL1Cmd->SetUnitCandidates("tesla kilogauss gause");
  BFieldFZBL1Cmd->AvailableForStates(G4State_Idle);  

  BFieldFZBR1Cmd = new G4UIcmdWithADoubleAndUnit("/field/setFZBR1Field",this);  
  BFieldFZBR1Cmd->SetGuidance("Set the quad    magnetic field for FZBR1");
  BFieldFZBR1Cmd->SetParameterName("Field",false);
  BFieldFZBR1Cmd->SetDefaultUnit("tesla");  
  BFieldFZBR1Cmd->SetUnitCandidates("tesla kilogauss gause");
  BFieldFZBR1Cmd->AvailableForStates(G4State_Idle);  

  BFieldFZBL2Cmd = new G4UIcmdWithADoubleAndUnit("/field/setFZBL2Field",this);  
  BFieldFZBL2Cmd->SetGuidance("Set the quad    magnetic field for FZBL2");
  BFieldFZBL2Cmd->SetParameterName("Field",false);
  BFieldFZBL2Cmd->SetDefaultUnit("tesla");  
  BFieldFZBL2Cmd->SetUnitCandidates("tesla kilogauss gause");
  BFieldFZBL2Cmd->AvailableForStates(G4State_Idle);  
  
  BFieldFZBR2Cmd = new G4UIcmdWithADoubleAndUnit("/field/setFZBR2Field",this);  
  BFieldFZBR2Cmd->SetGuidance("Set the quad    magnetic field for FZBR2");
  BFieldFZBR2Cmd->SetParameterName("Field",false);
  BFieldFZBR2Cmd->SetDefaultUnit("tesla");  
  BFieldFZBR2Cmd->SetUnitCandidates("tesla kilogauss gause");
  BFieldFZBR2Cmd->AvailableForStates(G4State_Idle);  
/*
 * Dipole field not yet adjustable.
  BFieldFZBL3Cmd = new G4UIcmdWithADoubleAndUnit("/field/setFZBL3Field",this);  
  BFieldFZBL3Cmd->SetGuidance("Set the uniform magnetic field for FZBL3");
  BFieldFZBL3Cmd->SetParameterName("Field",false);
  BFieldFZBL3Cmd->SetDefaultUnit("tesla");  
  BFieldFZBL3Cmd->SetUnitCandidates("tesla kilogauss gause");
  BFieldFZBL3Cmd->AvailableForStates(G4State_Idle);  
  
  BFieldFZBR3Cmd = new G4UIcmdWithADoubleAndUnit("/field/setFZBR3Field",this);  
  BFieldFZBR3Cmd->SetGuidance("Set the uniform magnetic field for FZBR3");
  BFieldFZBR3Cmd->SetParameterName("Field",false);
  BFieldFZBR3Cmd->SetDefaultUnit("tesla");  
  BFieldFZBR3Cmd->SetUnitCandidates("tesla kilogauss gause");
  BFieldFZBR3Cmd->AvailableForStates(G4State_Idle);  
*/
  BFieldFZBL4Cmd = new G4UIcmdWithADoubleAndUnit("/field/setFZBL4Field",this);  
  BFieldFZBL4Cmd->SetGuidance("Set the quad    magnetic field for FZBL4");
  BFieldFZBL4Cmd->SetParameterName("Field",false);
  BFieldFZBL4Cmd->SetDefaultUnit("tesla");  
  BFieldFZBL4Cmd->SetUnitCandidates("tesla kilogauss gause");
  BFieldFZBL4Cmd->AvailableForStates(G4State_Idle);  
  
  BFieldFZBR4Cmd = new G4UIcmdWithADoubleAndUnit("/field/setFZBR4Field",this);  
  BFieldFZBR4Cmd->SetGuidance("Set the quad    magnetic field for FZBR4");
  BFieldFZBR4Cmd->SetParameterName("Field",false);
  BFieldFZBR4Cmd->SetDefaultUnit("tesla");  
  BFieldFZBR4Cmd->SetUnitCandidates("tesla kilogauss gause");
  BFieldFZBR4Cmd->AvailableForStates(G4State_Idle);  
}

///////////////////////////////////////////////////////////////////////////////

HRSEMFieldSetupMessenger::~HRSEMFieldSetupMessenger()
{
  delete StepperCmd;
  delete UpdateCmd;
  delete MinStepCmd;
  delete HRSFieldDir;
  delete BFieldFZBL1Cmd;
  delete BFieldFZBR1Cmd;
  delete BFieldFZBL2Cmd;
  delete BFieldFZBR2Cmd;
  //delete BFieldFZBL3Cmd;
  //delete BFieldFZBR3Cmd;
  delete BFieldFZBL4Cmd;
  delete BFieldFZBR4Cmd;
 
}

////////////////////////////////////////////////////////////////////////////
//
//

void HRSEMFieldSetupMessenger::SetNewValue( G4UIcommand* command, G4String newValue)
{ 
  if( command == StepperCmd )
  { 
    fEMFieldSetup->SetStepperType(StepperCmd->GetNewIntValue(newValue));
  }  
  else if( command == UpdateCmd )
  { 
    fEMFieldSetup->UpdateField(); 
  }
  else if( command == MinStepCmd )
  { 
    fEMFieldSetup->SetMinStep(MinStepCmd->GetNewDoubleValue(newValue));
  }

  else if( command == BFieldFZBL1Cmd )
  { 
    fEMFieldSetup->SetBFieldFZBL1(BFieldFZBL1Cmd->GetNewDoubleValue(newValue));
  }
  else if( command == BFieldFZBR1Cmd )
  { 
    fEMFieldSetup->SetBFieldFZBR1(BFieldFZBR1Cmd->GetNewDoubleValue(newValue));
  }

  else if( command == BFieldFZBL2Cmd )
  { 
    fEMFieldSetup->SetBFieldFZBL2(BFieldFZBL2Cmd->GetNewDoubleValue(newValue));
  }
  else if( command == BFieldFZBR2Cmd )
  { 
    fEMFieldSetup->SetBFieldFZBR2(BFieldFZBR2Cmd->GetNewDoubleValue(newValue));
  }

/*  else if( command == BFieldFZBL3Cmd )
  { 
    fEMFieldSetup->SetBFieldFZBL3(BFieldFZBL3Cmd->GetNewDoubleValue(newValue));
  }
  else if( command == BFieldFZBR3Cmd )
  { 
    fEMFieldSetup->SetBFieldFZBR3(BFieldFZBR3Cmd->GetNewDoubleValue(newValue));
  }
*/
  else if( command == BFieldFZBL4Cmd )
  { 
    fEMFieldSetup->SetBFieldFZBL4(BFieldFZBL4Cmd->GetNewDoubleValue(newValue));
  }
  else if( command == BFieldFZBR4Cmd )
  { 
    fEMFieldSetup->SetBFieldFZBR4(BFieldFZBR4Cmd->GetNewDoubleValue(newValue));
  }
}

//
/////////////////////////////////////////////////////////////////////////
