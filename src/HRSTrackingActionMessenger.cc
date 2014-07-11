// ********************************************************************
//Jixie Zhang
// $Id: HRSTrackingActionMessenger.cc,v 1.0, 2010/12/26  HRS Exp $
// --------------------------------------------------------------
//

#include "HRSTrackingActionMessenger.hh"
#include "HRSTrackingAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"

HRSTrackingActionMessenger::HRSTrackingActionMessenger(HRSTrackingAction * mtarget)
:target(mtarget)
{
  trackIdCmd = new G4UIcmdWithAnInteger("/tracking/thisTrackOnly",this);
  trackIdCmd->SetGuidance("Output only this track into the ascii file and screen");
  trackIdCmd->SetGuidance("if outtrackid<0 then output all tracks; If trackid=0 no tracks will be printed");
  trackIdCmd->SetGuidance("You may use together with these cmds /tracking/particleOnly and /tracking/noSecondary to control the ascii output");
  trackIdCmd->SetParameterName("outtrackid",false);
  trackIdCmd->SetRange("outtrackid>=-1");
  trackIdCmd->SetDefaultValue(-1);

  noSecondaryCmd = new G4UIcmdWithAnInteger("/tracking/noSecondary",this);
  noSecondaryCmd->SetGuidance("Output secondary tracks or not. if noSecondary=1 then no secondary tracks will be printed");
  noSecondaryCmd->SetGuidance("You may use together with these cmds /tracking/thisTrackOnly and /tracking/particleOnly to control the ascii output");
  noSecondaryCmd->SetParameterName("noSecondary",true);
  noSecondaryCmd->SetRange("noSecondary>=0");
  noSecondaryCmd->SetDefaultValue(1);
  
  particleCmd = new G4UIcmdWithAString("/tracking/particleOnly",this);
  particleCmd->SetGuidance("Output this type of particle only. if particleName=\"all\" or particleName=\"-1\", all particle will be printed");
  particleCmd->SetGuidance("You may use together with these cmds /tracking/thisTrackOnly and /tracking/noSecondary to control the ascii output");
  particleCmd->SetParameterName("particleName",false);
}

HRSTrackingActionMessenger::~HRSTrackingActionMessenger()
{
  delete trackIdCmd;
  delete noSecondaryCmd;
  delete particleCmd;
}

void HRSTrackingActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==trackIdCmd )
  { target->SetOutTrackId(trackIdCmd->GetNewIntValue(newValue)); }
  else if( command==noSecondaryCmd )
  { target->SetNoSecondary(noSecondaryCmd->GetNewIntValue(newValue)); }
  else if( command==particleCmd )
  {//verified that this cmd is valid
	  if(newValue=="all" || newValue=="-1") 
	  {
		target->SetOutParticleName( "all" );
	  }
	  else
	  {
		  G4ParticleDefinition* pd = 
			  G4ParticleTable::GetParticleTable()->FindParticle(newValue);
		  if(pd != 0) target->SetOutParticleName( newValue );    
	  }
  }
}

G4String HRSTrackingActionMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==trackIdCmd )
  { cv = trackIdCmd->ConvertToString(target->GetOutTrackId()); }
  else if( command==noSecondaryCmd )
  { cv = trackIdCmd->ConvertToString(target->GetNoSecondary()); }
  else if( command==particleCmd )
  { cv = particleCmd->ConvertToString(target->GetOutParticleName()); }
  return cv;
}
