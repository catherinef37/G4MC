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
  /*
  BField3VFZB1Cmd = new G4UIcmdWith3VectorAndUnit("/field/setFZB1Field3V",this);  
  BField3VFZB1Cmd->SetGuidance("Set the quad    magnetic field for FZB1");
  BField3VFZB1Cmd->SetParameterName("Bx","By","Bz",false);
  BField3VFZB1Cmd->SetDefaultUnit("tesla");  
  BField3VFZB1Cmd->SetUnitCandidates("tesla kilogauss gause");
  BField3VFZB1Cmd->AvailableForStates(G4State_Idle);  

  BField3VFZB2Cmd = new G4UIcmdWith3VectorAndUnit("/field/setFZB2Field3V",this);  
  BField3VFZB2Cmd->SetGuidance("Set the quad    magnetic field for FZB2");
  BField3VFZB2Cmd->SetParameterName("Bx","By","Bz",false);
  BField3VFZB2Cmd->SetDefaultUnit("tesla");  
  BField3VFZB2Cmd->SetUnitCandidates("tesla kilogauss gause");
  BField3VFZB2Cmd->AvailableForStates(G4State_Idle);  
  */
  BField3VFZB3Cmd = new G4UIcmdWith3VectorAndUnit("/field/setFZB3Field3V",this);  
  BField3VFZB3Cmd->SetGuidance("Set the uniform magnetic field for FZB3");
  BField3VFZB3Cmd->SetParameterName("Bx","By","Bz",false);
  BField3VFZB3Cmd->SetDefaultUnit("tesla");  
  BField3VFZB3Cmd->SetUnitCandidates("tesla kilogauss gause");
  BField3VFZB3Cmd->AvailableForStates(G4State_Idle);  
  /*
  BField3VFZB4Cmd = new G4UIcmdWith3VectorAndUnit("/field/setFZB4Field3V",this);  
  BField3VFZB4Cmd->SetGuidance("Set the quad    magnetic field for FZB4");
  BField3VFZB4Cmd->SetParameterName("Bx","By","Bz",false);
  BField3VFZB4Cmd->SetDefaultUnit("tesla");  
  BField3VFZB4Cmd->SetUnitCandidates("tesla kilogauss gause");
  BField3VFZB4Cmd->AvailableForStates(G4State_Idle);  
  */
}

///////////////////////////////////////////////////////////////////////////////

HRSEMFieldSetupMessenger::~HRSEMFieldSetupMessenger()
{
  delete StepperCmd;
  delete UpdateCmd;
  delete MinStepCmd;
  delete HRSFieldDir;
  //delete BField3VFZB1Cmd;
  //delete BField3VFZB2Cmd;
  delete BField3VFZB3Cmd;
  //delete BField3VFZB4Cmd;
 
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
  //else if( command == BField3VFZB1Cmd )
    //{ 
    //fEMFieldSetup->SetBField3VFZB1(BField3VFZB1Cmd->GetNew3VectorValue(newValue));
    //}
  //else if( command == BField3VFZB2Cmd )
  //{ 
  //fEMFieldSetup->SetBField3VFZB2(BField3VFZB2Cmd->GetNew3VectorValue(newValue));
  //}
  //else if( command == BField3VFZB3Cmd )
    //{ 
    //fEMFieldSetup->SetBField3VFZB3(BField3VFZB3Cmd->GetNew3VectorValue(newValue));
    //}
  //else if( command == BField3VFZB4Cmd )
  //{ 
  //fEMFieldSetup->SetBField3VFZB4(BField3VFZB4Cmd->GetNew3VectorValue(newValue));
  //}
}

//
/////////////////////////////////////////////////////////////////////////
