// ********************************************************************
//
// $Id: HRSRunAction.cc,v 1.0, 2010/12/26 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//..............................................................................

#include "G4Timer.hh"
#include "HRSRunAction.hh"
#include "G4Run.hh"
#include <stdio.h>
#include "HRSRootTree.hh"
#include "HRSGlobal.hh"

extern  HRSRootTree* gHRSTree;
//..............................................................................

HRSRunAction::HRSRunAction()
{
	timer = new G4Timer;
	bRunActionVerbose=true;
	G4cout<<"HRSRunAction() construction done!"<<G4endl;
}

//..............................................................................

HRSRunAction::~HRSRunAction()
{
	delete timer;
	G4cout<<"delete HRSRunAction ... done!"<<G4endl;
}

//..............................................................................

void HRSRunAction::BeginOfRunAction(const G4Run* aRun)
{
	if(bRunActionVerbose && gHRSTree->GetTotalEvtID()>100*aRun->GetRunID() )
	{
		char strLog[255];
		sprintf(strLog,"### Run %d start...\n",aRun->GetRunID());
		WriteLog(strLog);
		timer->Start();
	}
}

//..............................................................................

void HRSRunAction::EndOfRunAction(const G4Run* aRun)
{
	if(bRunActionVerbose && gHRSTree->GetTotalEvtID()>100*aRun->GetRunID() )
	{
		timer->Stop();
		G4cout << "### Number of event = " << aRun->GetNumberOfEvent()
			<< " " << *timer << G4endl;
		char strLog[255];
		sprintf(strLog,"### Processed %d events in run %d. Total= %d  root= %d \n",
			aRun->GetNumberOfEvent(),aRun->GetRunID(),gHRSTree->GetTotalEvtID(),
			gHRSTree->GetRootEvtID());    
		WriteLog(strLog);
	}
}

//..............................................................................
