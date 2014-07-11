// ********************************************************************
//
// $Id: HRSMC.cc,v 1.05, 2013/10/06  HRS Exp $
// --------------------------------------------------------------

#include <stdlib.h>
#include "HRSRootTree.hh"
#include "G4RunManager.hh"

#include "HRSUIExecutive.hh"
#include "G4UImanager.hh"

#include "HRSPhysicsList.hh"		//Jixie's physics model
#include "PhysicsList.hh"

#include "HRSMaterial.hh"
#include "HRSDetectorConstruction.hh"
#include "HRSPrimaryGeneratorAction.hh"

#include "HRSRunAction.hh"
#include "HRSEventAction.hh"
#include "HRSTrackingAction.hh"
#include "HRSSteppingAction.hh"
#include "HRSSteppingVerbose.hh"

#include "UsageManager.hh"
#include "HRSGlobal.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//global variables
UsageManager* gConfig=0; //will be initialized in main() before all other classes
HRSRootTree* gHRSTree=0; //will be initialized in main() after G4RunManager start up before user-action start

int main(int argc,char** argv)
{
	////////////////////////////////////////////////////////////////////
	//take care of the arguments
	////////////////////////////////////////////////////////////////////
	string LogFileName("G2P_G4Sim_Log.txt");
	gConfig=new UsageManager("HRSUsage.ini",argc,argv,LogFileName.c_str());
	gConfig->ReadFile("Detector.ini");
	//gConfig->PrintOpt();
	//gConfig->PrintParamMap(); 

	////////////////////////////////////////////////////////////////////
	//my Stepping verbose output class
	////////////////////////////////////////////////////////////////////
	G4VSteppingVerbose::SetInstance(new HRSSteppingVerbose);

	////////////////////////////////////////////////////////////////////
	// RunManager construction
	////////////////////////////////////////////////////////////////////
	G4RunManager* runManager = new G4RunManager;


	////////////////////////////////////////////////////////////////////
	// Visualization manager construction	
	////////////////////////////////////////////////////////////////////
#ifdef G4VIS_USE
	G4VisManager* visManager = new G4VisExecutive();
	//default verbose level= warning, 
	visManager->Initialize();
#endif

	////////////////////////////////////////////////////////////////////
	// mandatory user initialization classes
	////////////////////////////////////////////////////////////////////
	//initial materials
	HRSMaterial* pMaterialManager=new HRSMaterial();
	runManager->SetUserInitialization(new HRSDetectorConstruction);

	//inital physics list
	int pUseJixieModel=0;
	gConfig->GetArgument("UseJixieModel",pUseJixieModel);
	if(pUseJixieModel)
	{
		//Use Jixie's model
		runManager->SetUserInitialization(new HRSPhysicsList);
	}
	else
	{
		//or use this one which is from example/extended/hadon/had01
		G4String pPhysicsModel=gConfig->GetArgument("PhysicsModel");
		runManager->SetUserInitialization(new PhysicsList(pPhysicsModel));
	}

	// initialize Geant4 kernel
	runManager->Initialize();
	// mandatory user action class
	runManager->SetUserAction(new HRSPrimaryGeneratorAction);

	// The analysis manager
	int pRunNumber=1;
	gConfig->GetArgument("RunNumber",pRunNumber);
	gHRSTree=new HRSRootTree(pRunNumber); 
	//if(pRunNumber<0) gHRSTree->SetRunNumber(abs(pRunNumber)); 

	////////////////////////////////////////////////////////////////////
	// optional user action classes, 
	////////////////////////////////////////////////////////////////////
	// runManager will delete them during closing
	runManager->SetUserAction(new HRSRunAction);
	runManager->SetUserAction(new HRSEventAction);
	runManager->SetUserAction(new HRSTrackingAction);
	runManager->SetUserAction(new HRSSteppingAction);

	////////////////////////////////////////////////////////////////////
	// User interactions
	//////////////////////////////////////////////////////////////////// 
	// execute macro files if provided through arguments
	// Define (G)UI for interactive mode
	// The UImanager will be handled by Ranmanager,no need to delete it manually 
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	char cmd[255],tmpstr[255];
	//take care of cmd files then the mcin file
	int pNofMac=0;
	gConfig->GetArgument("NofMac",pNofMac); 
	std::string MacFile;
	for(int i=0;i<pNofMac;i++)
	{
		sprintf(tmpstr,"MacFile%d",i+1);
		MacFile=gConfig->GetArgument(tmpstr);
		if(CheckFile(MacFile))
		{
			sprintf(cmd,"/control/execute %s",MacFile.c_str()); 
			G4cout<<"Trying to execute cmd: "<<cmd<<G4endl;
			UImanager->ApplyCommand(cmd);
		}
	}

	//By Jixie: Might change the trigger for this part later
	//gConfig->GetArgument("UseRootNtuple",pUseRootNtuple); 
	//if(UseRootNtuple) //root ntuple mode	
	string pPrimaryEngine1=gConfig->GetArgument("PrimaryEngine1"); 
	if(pPrimaryEngine1=="RootNtuple") //root ntuple mode	
	{		
		int pTrigNum=-1;
		gConfig->GetArgument("TrigNum1",pTrigNum);
		if(pTrigNum<=0) pTrigNum=9999999;
		sprintf(cmd,"/run/beamOn %d",pTrigNum); 
		G4cout<<"Trying to execute cmd: "<<cmd<<G4endl;
		UImanager->ApplyCommand(cmd);
	}

	////////////////////////////////////////////////////////////////////
	// interaction mode
	////////////////////////////////////////////////////////////////////
	int pInteractiveMode=0;
	gConfig->GetArgument("InteractiveMode",pInteractiveMode); 
	if(pInteractiveMode)
	{
		int pUseGui=0;
		gConfig->GetArgument("UseGui",pUseGui);
		HRSUIExecutive *pUI = new HRSUIExecutive(argc, argv, pUseGui);
		if(pUI->IsGUI())
		{
			//put your extra cmd here;
			if(CheckFile("gui.mac"))
			{
				UImanager->ApplyCommand("/control/execute gui.mac");
			}
		}
		pUI->SessionStart();
	}

	////////////////////////////////////////////////////////////////////
	//free the memory
	////////////////////////////////////////////////////////////////////
	if(gHRSTree) delete gHRSTree;

#ifdef G4VIS_USE
	delete visManager;
	G4cout << "Vis manager deleted..." << G4endl;
#endif

	// Free the memory: user actions, physics_list and detector_description are
	//                  owned and deleted by the run manager, so they should not
	//                  be deleted in the main() program !

	delete pMaterialManager;

	delete runManager;
	G4cout << "Run manager deleted..." << G4endl;

	return 0;
}

