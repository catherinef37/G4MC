//
// $Id: HRSTrackingAction.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
#include "HRSTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "HRSTrackInformation.hh"
#include "HRSTrackingActionMessenger.hh"
#include "HRSTrajectory.hh"
#include "UsageManager.hh"
extern UsageManager* gConfig;

HRSTrackingAction::HRSTrackingAction()
{
	iOutTrackId=-1; //output all tracks by default
	iNoSecondary=1;
	gConfig->GetArgument("NoSecondary",iNoSecondary);

	strOutParticle="all"; //output all particles by default
	messenger=new HRSTrackingActionMessenger(this);

	G4cout<<"HRSTrackingAction() construction done!"<<G4endl;
}

HRSTrackingAction::~HRSTrackingAction()
{
	delete 	messenger;
	G4cout<<"delete HRSTrackingAction ... done!"<<G4endl;
}

void HRSTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
	//For primary particles, even if I don't want the output of the step, I should not kill 
	//these tracks here because I want them to be processed till out of the detector
	//For secondary particles, I will kill them if iNoSecondary is set.
	//control tracking behavious
	//all of these will be done in HRSSteppingAction, What I need to do is just 
	//set iTrackOutputStatus to -1 in HRSTrackingInformation

	//G4cout<<"*******debug******HRSTrackingAction::iOutTrackId="<<iOutTrackId<<"***\n";
	G4int iTrackOutputStatus=1;
	//trackid control
	if(iOutTrackId<0) //output all tracks
	{	
		;	
	}
	else if(iOutTrackId==0) //output no track
	{
		iTrackOutputStatus=0;
	}
	else if(iOutTrackId!=aTrack->GetTrackID())
	{//still process this track but no step output	
		iTrackOutputStatus=0;
	}

	//G4cout<<"*******debug******HRSTrackingAction::strOutParticle="<<strOutParticle<<"***\n";
	//track particle control
	if(strOutParticle=="all" || strOutParticle=="-1") //output all particle tracks
	{	
		;	
	}
	else if(strOutParticle!=aTrack->GetDefinition()->GetParticleName()) 	
	{//output this particle's track only	
		iTrackOutputStatus=0;		
	}

	// Create trajectory only for primaries, these trajectories can be used at EndofEvent
	if(aTrack->GetParentID()==0)
	{ 
		fpTrackingManager->SetStoreTrajectory(1); 
	}
	else 
	{ 	
		if(iNoSecondary)
		{			
			fpTrackingManager->SetStoreTrajectory(0);
			//2 ways to stop 2ndary			
			//(1)stop this track here or
			fpTrackingManager->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
			return;
			//(2)will invoke aTrack->SetTrackStatus(fStopAndKill) again in stepping action
			//iTrackOutputStatus=-1;
		}
	}

	if(aTrack->GetUserInformation()==0)
	{
		HRSTrackInformation* anInfo = new HRSTrackInformation(aTrack);					
		G4Track* theTrack = (G4Track*)aTrack;
		theTrack->SetUserInformation(anInfo);
		/*
		//but we can also do it in this 2 ways
		(1)fpTrackingManager->SetUserTrackInformation(anInfo); 
		(2)fpTrackingManager->GetTrack()->SetUserInformation(anInfo);		
		*/
		anInfo->SetStepOutputStatus(iTrackOutputStatus);
		//do not print this header
		//if(iTrackOutputStatus>0) anInfo->Print(aTrack);  //print track header
	}
	else
	{
		HRSTrackInformation* anInfo = (HRSTrackInformation* )(aTrack->GetUserInformation());
		anInfo->SetStepOutputStatus(iTrackOutputStatus);
	}

	fpTrackingManager->SetTrajectory(new HRSTrajectory(aTrack));
}

void HRSTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
	G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
	if(secondaries)
	{
		HRSTrackInformation* info = (HRSTrackInformation*)(aTrack->GetUserInformation());
		size_t nSeco = secondaries->size();
		if(nSeco>0)
		{
			for(size_t i=0;i<nSeco;i++)
			{
				HRSTrackInformation* infoNew = new HRSTrackInformation(info);
				(*secondaries)[i]->SetUserInformation(infoNew);
			}
		}
	}
}
