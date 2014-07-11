// ********************************************************************
//
// $Id: HRSEventAction.cc,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
#include "math.h"
#include "HRSEventAction.hh"
#include "HRSEventActionMessenger.hh"

#ifdef G4ANALYSIS_USE
#include "HRSAnalysisManager.hh"
#endif // G4ANALYSIS_USE

#include "HRSTrajectory.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4PrimaryParticle.hh"
#include "G4ThreeVector.hh"
#include "HRSStdHit.hh"
#include "HRSDCHit.hh"
#include "HRSCalorimeterHit.hh" 
#include "HRSRootTree.hh"
#include "UsageManager.hh"

//I am not using this enum list any longer
//remember to add your SD name into this emum list
enum SensitiveDetectorList {  
	FZB1VD=1,FZB2VD=2,
	sieveSlit=11,septumWindow=12,
	thirdArmShielding=30, thirdArmSC1=31, thirdArmSC2=32, 
	virtualDetector=40,
	BBSC1=51,BBSC2=52,	
	VETO=60,NDSC1=61,NDSC2=62,NDSC3=63,NDSC4=64,	
	CREXUpBlock=70,CREXTarget=71,CREXDownBlock=72,
	virtualBoundary=99,
	RTPC=200
};

extern  HRSRootTree* gHRSTree;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//verbose level: 
//0:  nothing will be printed
//1:  HRSStdHit::Print() will be called
//2:  Event header will be printed, as well as all primary particle info (print togather). 
//3:  HRSTrajectory::ShowTrajectory() will be called;
//4:  HRSEventAction::PrintPrimary() will be call;

HRSEventAction::HRSEventAction()
{
	verboseLevel = 0;
	messenger = new HRSEventActionMessenger(this);

	G4cout<<"HRSEventAction() construction done!"<<G4endl;
}

HRSEventAction::~HRSEventAction()
{
	delete messenger;
	G4cout<<"delete HRSEventAction ... done!"<<G4endl;
}

void HRSEventAction::BeginOfEventAction(const G4Event* evt)
{
	for(int i=0;i<MaxStoredTrjN;i++) mStoreTrackIdList[i]=false;
	if (verboseLevel<2) return;

	//print the following info
	//>>>Begin of Event       1:  Vertex=(10.000,0.000,15.000)mm
	//*********************************************************************************************
	//>>> Track 1: e- Theta=6.000deg    Phi=-90.000deg   Momentum(-0.000,-99.302,944.796)=950.000MeV
	//>>> Track 2: e- Theta=6.000deg    Phi=-90.000deg   Momentum(-0.000,-99.302,944.796)=950.000MeV
	//*********************************************************************************************
	//mm and MeV are the default units
	G4int prec = G4cout.precision(3);
	G4ThreeVector pPosition = evt->GetPrimaryVertex(0)->GetPosition();
	G4cout<<"\n>>>Begin of Event "<<setw(6)<<evt->GetEventID() 
		<< std::fixed<<":  Vertex=" << pPosition<< "mm"<<G4endl;
	G4cout<< "***********************************************************************************************\n";
	for(int i=0;i<evt->GetNumberOfPrimaryVertex();i++)
	{
		G4PrimaryParticle* primary = evt->GetPrimaryVertex(i)->GetPrimary(0);
		G4ThreeVector pMomentum = primary->GetMomentum(); 
		G4cout<<">>>Track "<< std::setw(3)<< i+1 <<":  "
			<<std::setw(10)<< primary->GetG4code()->GetParticleName()
			<< "   Theta=" << pMomentum.theta()/deg<< "deg "
			<< "   Phi=" << pMomentum.phi()/deg<< "deg "
			<< std::fixed<< "  Momentum=" << pMomentum
			<< "="<< std::fixed<< pMomentum.mag()<< "MeV "<< G4endl;
	}
	G4cout<< "***********************************************************************************************\n";
	G4cout.precision(prec);
}

void HRSEventAction::ProcessDCSDHits(G4HCofThisEvent* HCE)
{
	if(!HCE) return;  //no hits collection found 

	int nHC = HCE->GetNumberOfCollections();

	G4String colName;
	HRSDCHitsCollection* theHC = 0;

	for(int ihc=0;ihc<nHC;ihc++)
	{
		//only one way to get a hit collection
		//G4VHitsCollection* G4HCofThisEvent::GetHC(G4int sdid);
		//But there are two ways to get a hit collection id (sdid)
		//
		//usually the SDID came from the following 2 ways
		//	G4SDManager* SDman = G4SDManager::GetSDMpointer();
		//  G4int G4SDManager::GetCollectionID(G4String colName);
		//	G4int G4SDManager::GetCollectionID(G4VHitsCollection * aHC);
		//  These two methods return the ID number of the sensitive detector.

		//Just for debug
		//it has been checked by the next 3 lines that ihc is equal to SDID
		//in other words, SDID is the array index of each SD, which totally depends
		//on how the code is organized. If one changes the detector construction
		//by adding or deleting a SD somewhere, the SDID will change.
		//
		//G4int theSDID = SDman->GetCollectionID(theHC);
		//G4cout<<"HRSEventAction::ProcessDCHits():  ihc="<<ihc<<"   SDID="<<theSDID<<G4endl;
		//G4cout<<"HRSEventAction::ProcessDCHits():  SDName="<<theHC->GetSDname()
		//	<<" colName="<<theHC->GetName()<<G4endl;

		int theSDID = ihc;
		colName = HCE->GetHC(ihc)->GetName();

		bool bUseUserSDIDTable=false;
		if(bUseUserSDIDTable)
		{
			//G4String SDName = HCE->GetHC(ihc)->GetSDname();
			//need to add int GetUserSDID(const char *SDName);
			//theSDID = GetUserSDID(SDName.c_str());
		}


		if(colName=="DCColl")
		{			
			theHC = (HRSDCHitsCollection*) HCE->GetHC(ihc);

			HRSDCHit *aHit=0; 
			G4ThreeVector tmpV3;

			G4int n_hit = theHC->entries();

			if(verboseLevel>=1) 
			{
				if(n_hit>0) G4cout<< (*theHC)[0]->GetPhysV()->GetName()<< " has " << n_hit << " hits." << G4endl;
			}

			for(int ih=0;ih<n_hit;ih++)
			{
				aHit=(*theHC)[ih];
				if(verboseLevel>=2) aHit->Print();

				if(aHit->GetTrackId()<MaxStoredTrjN)
				{
					mStoreTrackIdList[aHit->GetParentTrackId()]=true;
					mStoreTrackIdList[aHit->GetTrackId()]=true;
				}

				int ii=gHRSTree->SD_N;
				if(ii<MaxSDHit)
				{
					gHRSTree->SD_Id[ii]=theSDID*1.0E06+aHit->GetId();   //one hit one interaction
					gHRSTree->SD_Pid[ii]=aHit->GetPdgid(); 
					gHRSTree->SD_Tid[ii]=aHit->GetTrackId(); 
					gHRSTree->SD_ParentTid[ii]=aHit->GetParentTrackId(); 
					gHRSTree->SD_T[ii]=aHit->GetTime()/ns; 
					tmpV3=aHit->GetInPos();
					gHRSTree->SD_X[ii]=tmpV3.x()/mm;
					gHRSTree->SD_Y[ii]=tmpV3.y()/mm; 
					gHRSTree->SD_Z[ii]=tmpV3.z()/mm; 
					gHRSTree->SD_Edep[ii]=aHit->GetEdep()/MeV;
					gHRSTree->SD_NonIonEdep[ii]=aHit->GetNonIonEdep()/MeV;
					tmpV3=aHit->GetInMom();
					gHRSTree->SD_P[ii]=tmpV3.mag()/GeV;  
					gHRSTree->SD_Theta[ii]=tmpV3.theta(); 
					gHRSTree->SD_Phi[ii]=tmpV3.phi(); 
					gHRSTree->SD_Pout[ii]=aHit->GetOutMom().mag()/GeV; 

					gHRSTree->SD_N++; 
				}
			}

		} //end of for ih 
	} //end of ihc
}


void HRSEventAction::ProcessStdSDHits(G4HCofThisEvent* HCE)
{
	if(!HCE) return;  //no hits collection found 

	int nHC = HCE->GetNumberOfCollections();

	G4String colName, SDName;

	HRSStdHitsCollection* theStdHC = 0;

	for(int ihc=0;ihc<nHC;ihc++)
	{
		//only one way to get a hit collection
		//G4VHitsCollection* G4HCofThisEvent::GetHC(G4int sdid);
		//But there are two ways to get a hit collection id (sdid)
		//
		//usually the SDID came from the following 2 ways
		//  G4SDManager* SDman = G4SDManager::GetSDMpointer();
		//  G4int G4SDManager::GetCollectionID(G4String colName);
		//	G4int G4SDManager::GetCollectionID(G4VHitsCollection * aHC);
		//  These two methods return the ID number of the sensitive detector.

		//Just for debug
		//it has been checked by the next 3 lines that ihc is equal to SDID
		//in other words, SDID is the array index of each SD, which totally depends
		//on how the code is organized. If one changes the detector construction
		//by adding or deleting a SD somewhere, the SDID will change.
		//
		//G4int theSDID = SDman->GetCollectionID(theHC);
		//G4cout<<"HRSEventAction::ProcessDCHits():  ihc="<<ihc<<"   SDID="<<theSDID<<G4endl;
		//G4cout<<"HRSEventAction::ProcessDCHits():  SDName="<<theHC->GetSDname()
		//	<<" colName="<<theHC->GetName()<<G4endl;

		int theSDID = ihc;
		colName = HCE->GetHC(ihc)->GetName();
		SDName = HCE->GetHC(ihc)->GetSDname();

		bool bUseUserSDIDTable=false;
		if(bUseUserSDIDTable)
		{
			//need to add int GetUserSDID(const char *SDName);
			//theSDID = GetUserSDID(SDName.c_str());
		}

		if(colName=="StdColl")
		{
			theStdHC = (HRSStdHitsCollection*) HCE->GetHC(ihc);
			HRSStdHit *aHit=0; 
			G4ThreeVector tmpV3;

			G4int n_hit = theStdHC->entries();

			if(verboseLevel>=1) 
			{
				if(n_hit>0) G4cout<< (*theStdHC)[0]->GetPhysV()->GetName()<< " has " << n_hit << " hits." << G4endl;
			}

			for(int ih=0;ih<n_hit;ih++)
			{
				aHit=(*theStdHC)[ih];
				if(verboseLevel>=2) aHit->Print();

				if(aHit->GetTrackId()<MaxStoredTrjN)
				{
					mStoreTrackIdList[aHit->GetParentTrackId()]=true;
					mStoreTrackIdList[aHit->GetTrackId()]=true;
				}

				if(SDName=="thirdArmSC1")
				{
					int ii=gHRSTree->TA1_N;
					if(ii<MaxSDHit)
					{
						gHRSTree->TA1_Pid[ii]=aHit->GetPdgid(); 
						gHRSTree->TA1_Tid[ii]=aHit->GetTrackId(); 
						gHRSTree->TA1_ParentTid[ii]=aHit->GetParentTrackId(); 
						gHRSTree->TA1_T[ii]=aHit->GetTime()/ns; 
						tmpV3=aHit->GetInPos();
						gHRSTree->TA1_X[ii]=tmpV3.x()/mm;
						gHRSTree->TA1_Y[ii]=tmpV3.y()/mm; 
						gHRSTree->TA1_Z[ii]=tmpV3.z()/mm; 
						gHRSTree->TA1_Edep[ii]=aHit->GetEdep()/MeV;
						gHRSTree->TA1_NonIonEdep[ii]=aHit->GetNonIonEdep()/MeV;
						tmpV3=aHit->GetInMom();
						gHRSTree->TA1_P[ii]=tmpV3.mag()/GeV;  
						gHRSTree->TA1_Theta[ii]=tmpV3.theta(); 
						gHRSTree->TA1_Phi[ii]=tmpV3.phi(); 
						gHRSTree->TA1_Pout[ii]=aHit->GetOutMom().mag()/GeV;  

						gHRSTree->TA1_N++; 
					}
				}
				else if(SDName=="thirdArmSC2")
				{
					int ii=gHRSTree->TA2_N;
					if(ii<MaxSDHit)
					{
						gHRSTree->TA2_Pid[ii]=aHit->GetPdgid(); 
						gHRSTree->TA2_Tid[ii]=aHit->GetTrackId(); 
						gHRSTree->TA2_ParentTid[ii]=aHit->GetParentTrackId(); 
						gHRSTree->TA2_T[ii]=aHit->GetTime()/ns; 
						tmpV3=aHit->GetInPos();
						gHRSTree->TA2_X[ii]=tmpV3.x()/mm;
						gHRSTree->TA2_Y[ii]=tmpV3.y()/mm; 
						gHRSTree->TA2_Z[ii]=tmpV3.z()/mm; 
						gHRSTree->TA2_Edep[ii]=aHit->GetEdep()/MeV;
						gHRSTree->TA2_NonIonEdep[ii]=aHit->GetNonIonEdep()/MeV;
						tmpV3=aHit->GetInMom();
						gHRSTree->TA2_P[ii]=tmpV3.mag()/GeV;  
						gHRSTree->TA2_Theta[ii]=tmpV3.theta(); 
						gHRSTree->TA2_Phi[ii]=tmpV3.phi(); 
						gHRSTree->TA2_Pout[ii]=aHit->GetOutMom().mag()/GeV; 

						gHRSTree->TA2_N++; 
					}
				} //end of TA2
				else
				{
					int ii=gHRSTree->SD_N;
					if(ii<MaxSDHit)
					{
						gHRSTree->SD_Id[ii]=theSDID*1.0E06+aHit->GetId();   //one hit one track
						gHRSTree->SD_Pid[ii]=aHit->GetPdgid(); 
						gHRSTree->SD_Tid[ii]=aHit->GetTrackId(); 
						gHRSTree->SD_ParentTid[ii]=aHit->GetParentTrackId(); 
						gHRSTree->SD_T[ii]=aHit->GetTime()/ns; 
						tmpV3=aHit->GetInPos();
						gHRSTree->SD_X[ii]=tmpV3.x()/mm;
						gHRSTree->SD_Y[ii]=tmpV3.y()/mm; 
						gHRSTree->SD_Z[ii]=tmpV3.z()/mm; 
						gHRSTree->SD_Edep[ii]=aHit->GetEdep()/MeV;
						gHRSTree->SD_NonIonEdep[ii]=aHit->GetNonIonEdep()/MeV;
						tmpV3=aHit->GetInMom();
						gHRSTree->SD_P[ii]=tmpV3.mag()/GeV;  
						gHRSTree->SD_Theta[ii]=tmpV3.theta(); 
						gHRSTree->SD_Phi[ii]=tmpV3.phi(); 
						gHRSTree->SD_Pout[ii]=aHit->GetOutMom().mag()/GeV; 

						gHRSTree->SD_N++; 
					}
				}
			} //end of for ih 

		} //end of  if(colName=="StdColl")

	} //end of ihc
}


void HRSEventAction::ProcessCalorimeterSDHits(G4HCofThisEvent* HCE)
{
	if(!HCE) return;  //no hits collection found 

	int nHC = HCE->GetNumberOfCollections();

	G4String colName;
	HRSCalorimeterHitsCollection* theHC = 0;

	for(int ihc=0;ihc<nHC;ihc++)
	{
		//only one way to get a hit collection
		//G4VHitsCollection* G4HCofThisEvent::GetHC(G4int sdid);
		//But there are two ways to get a hit collection id (sdid)
		//
		//usually the SDID came from the following 2 ways
		//	G4SDManager* SDman = G4SDManager::GetSDMpointer();
		//  G4int G4SDManager::GetCollectionID(G4String colName);
		//	G4int G4SDManager::GetCollectionID(G4VHitsCollection * aHC);
		//  These two methods return the ID number of the sensitive detector.

		//Just for debug
		//it has been checked by the next 3 lines that ihc is equal to SDID
		//in other words, SDID is the array index of each SD, which totally depends
		//on how the code is organized. If one changes the detector construction
		//by adding or deleting a SD somewhere, the SDID will change.
		//
		//G4int theSDID = SDman->GetCollectionID(theHC);
		//G4cout<<"HRSEventAction::ProcessDCHits():  ihc="<<ihc<<"   SDID="<<theSDID<<G4endl;
		//G4cout<<"HRSEventAction::ProcessDCHits():  SDName="<<theHC->GetSDname()
		//	<<" colName="<<theHC->GetName()<<G4endl;

		int theSDID = ihc;
		colName = HCE->GetHC(ihc)->GetName();

		bool bUseUserSDIDTable=false;
		if(bUseUserSDIDTable)
		{
			//G4String SDName = HCE->GetHC(ihc)->GetSDname();
			//need to add int GetUserSDID(const char *SDName);
			//theSDID = GetUserSDID(SDName.c_str());
		}


		if(colName=="CalorimeterColl")
		{			
			theHC = (HRSCalorimeterHitsCollection*) HCE->GetHC(ihc);

			HRSCalorimeterHit *aHit=0; 
			G4ThreeVector tmpV3;

			G4int n_hit = theHC->entries();

			if(verboseLevel>=1) 
			{
				if(n_hit>0) G4cout<< (*theHC)[0]->GetPhysV()->GetName()<< " has " << n_hit << " hits." << G4endl;
			}

			for(int ih=0;ih<n_hit;ih++)
			{
				aHit=(*theHC)[ih];
				if(verboseLevel>=2) aHit->Print();

				if(aHit->GetTrackId()<MaxStoredTrjN)
				{
					mStoreTrackIdList[aHit->GetParentTrackId()]=true;
					mStoreTrackIdList[aHit->GetTrackId()]=true;
				}

				int ii=gHRSTree->SD_N;
				if(ii<MaxSDHit)
				{
					gHRSTree->SD_Id[ii]=theSDID;   //one hit one event, no need to multiple 1.0E6
					gHRSTree->SD_Pid[ii]=aHit->GetPdgid(); 
					gHRSTree->SD_Tid[ii]=aHit->GetTrackId(); 
					gHRSTree->SD_ParentTid[ii]=aHit->GetParentTrackId(); 
					gHRSTree->SD_T[ii]=aHit->GetTime()/ns; 
					tmpV3=aHit->GetInPos();
					gHRSTree->SD_X[ii]=tmpV3.x()/mm;
					gHRSTree->SD_Y[ii]=tmpV3.y()/mm; 
					gHRSTree->SD_Z[ii]=tmpV3.z()/mm; 
					gHRSTree->SD_Edep[ii]=aHit->GetEdep()/MeV;
					gHRSTree->SD_NonIonEdep[ii]=aHit->GetNonIonEdep()/MeV;
					tmpV3=aHit->GetInMom();
					gHRSTree->SD_P[ii]=tmpV3.mag()/GeV;  
					gHRSTree->SD_Theta[ii]=tmpV3.theta(); 
					gHRSTree->SD_Phi[ii]=tmpV3.phi(); 
					gHRSTree->SD_Pout[ii]=aHit->GetOutMom().mag()/GeV; 

					gHRSTree->SD_N++; 
				}
			}

		} //end of for ih 
	} //end of ihc
}


void HRSEventAction::ProcessSDHits(G4HCofThisEvent* HCE)
{
	if(!HCE) return;  //no hits collection found 

	int nHC = HCE->GetNumberOfCollections();

	G4String colName, SDName;

	HRSDCHitsCollection* theDCHC = 0;
	HRSStdHitsCollection* theStdHC = 0;
	HRSCalorimeterHitsCollection* theECHC = 0;

	for(int ihc=0;ihc<nHC;ihc++)
	{
		//only one way to get a hit collection
		//G4VHitsCollection* G4HCofThisEvent::GetHC(G4int sdid);
		//But there are two ways to get a hit collection id (sdid)
		//
		//usually the SDID came from the following 2 ways
		//  G4SDManager* SDman = G4SDManager::GetSDMpointer();
		//  G4int G4SDManager::GetCollectionID(G4String colName);
		//	G4int G4SDManager::GetCollectionID(G4VHitsCollection * aHC);
		//  These two methods return the ID number of the sensitive detector.

		//Just for debug
		//it has been checked by the next 3 lines that ihc is equal to SDID
		//in other words, SDID is the array index of each SD, which totally depends
		//on how the code is organized. If one changes the detector construction
		//by adding or deleting a SD somewhere, the SDID will change.
		//
		//G4int theSDID = SDman->GetCollectionID(theHC);
		//G4cout<<"HRSEventAction::ProcessDCHits():  ihc="<<ihc<<"   SDID="<<theSDID<<G4endl;
		//G4cout<<"HRSEventAction::ProcessDCHits():  SDName="<<theHC->GetSDname()
		//	<<" colName="<<theHC->GetName()<<G4endl;


		int theSDID = ihc;
		colName = HCE->GetHC(ihc)->GetName();
		SDName = HCE->GetHC(ihc)->GetSDname();

		//this block allows you to set SDID to a given number
		bool bUseUserSDIDTable=false;
		if(bUseUserSDIDTable)
		{
			//need to add int GetUserSDID(const char *SDName);
			//theSDID = GetUserSDID(SDName.c_str());
		}


		if(colName=="DCColl")
		{
			theDCHC = (HRSDCHitsCollection*) HCE->GetHC(ihc);
			HRSDCHit *aHit=0; 
			G4ThreeVector tmpV3;

			G4int n_hit = theDCHC->entries();

			if(verboseLevel>=1) 
			{
				if(n_hit>0) G4cout<< (*theDCHC)[0]->GetPhysV()->GetName()<< " has " << n_hit << " hits." << G4endl;
			}

			for(int ih=0;ih<n_hit;ih++)
			{
				aHit=(*theDCHC)[ih];
				if(verboseLevel>=2) aHit->Print();

				if(aHit->GetTrackId()<MaxStoredTrjN)
				{
					mStoreTrackIdList[aHit->GetParentTrackId()]=true;
					mStoreTrackIdList[aHit->GetTrackId()]=true;
				}

				int ii=gHRSTree->SD_N;
				if(ii<MaxSDHit)
				{
					gHRSTree->SD_Id[ii]=theSDID*1.0E06+aHit->GetId();   //one hit one interaction
					gHRSTree->SD_Pid[ii]=aHit->GetPdgid(); 
					gHRSTree->SD_Tid[ii]=aHit->GetTrackId(); 
					gHRSTree->SD_ParentTid[ii]=aHit->GetParentTrackId(); 
					gHRSTree->SD_T[ii]=aHit->GetTime()/ns; 
					tmpV3=aHit->GetInPos();
					gHRSTree->SD_X[ii]=tmpV3.x()/mm;
					gHRSTree->SD_Y[ii]=tmpV3.y()/mm; 
					gHRSTree->SD_Z[ii]=tmpV3.z()/mm; 
					gHRSTree->SD_Edep[ii]=aHit->GetEdep()/MeV;
					gHRSTree->SD_NonIonEdep[ii]=aHit->GetNonIonEdep()/MeV;
					tmpV3=aHit->GetInMom();
					gHRSTree->SD_P[ii]=tmpV3.mag()/GeV;  
					gHRSTree->SD_Theta[ii]=tmpV3.theta(); 
					gHRSTree->SD_Phi[ii]=tmpV3.phi(); 
					gHRSTree->SD_Pout[ii]=aHit->GetOutMom().mag()/GeV; 

					gHRSTree->SD_N++; 
				}
			} //end of for ih 

		} //end of  if(colName=="DCColl")
		else if(colName=="StdColl")
		{
			theStdHC = (HRSStdHitsCollection*) HCE->GetHC(ihc);
			HRSStdHit *aHit=0; 
			G4ThreeVector tmpV3;

			G4int n_hit = theStdHC->entries();

			if(verboseLevel>=1) 
			{
				if(n_hit>0) G4cout<< (*theStdHC)[0]->GetPhysV()->GetName()<< " has " << n_hit << " hits." << G4endl;
			}

			for(int ih=0;ih<n_hit;ih++)
			{
				aHit=(*theStdHC)[ih];
				if(verboseLevel>=2) aHit->Print();

				if(aHit->GetTrackId()<MaxStoredTrjN)
				{
					mStoreTrackIdList[aHit->GetParentTrackId()]=true;
					mStoreTrackIdList[aHit->GetTrackId()]=true;
				}

				if(SDName=="thirdArmSC1")
				{
					int ii=gHRSTree->TA1_N;
					if(ii<MaxSDHit)
					{
						gHRSTree->TA1_Pid[ii]=aHit->GetPdgid(); 
						gHRSTree->TA1_Tid[ii]=aHit->GetTrackId(); 
						gHRSTree->TA1_ParentTid[ii]=aHit->GetParentTrackId(); 
						gHRSTree->TA1_T[ii]=aHit->GetTime()/ns; 
						tmpV3=aHit->GetInPos();
						gHRSTree->TA1_X[ii]=tmpV3.x()/mm;
						gHRSTree->TA1_Y[ii]=tmpV3.y()/mm; 
						gHRSTree->TA1_Z[ii]=tmpV3.z()/mm; 
						gHRSTree->TA1_Edep[ii]=aHit->GetEdep()/MeV;
						gHRSTree->TA1_NonIonEdep[ii]=aHit->GetNonIonEdep()/MeV;
						tmpV3=aHit->GetInMom();
						gHRSTree->TA1_P[ii]=tmpV3.mag()/GeV;  
						gHRSTree->TA1_Theta[ii]=tmpV3.theta(); 
						gHRSTree->TA1_Phi[ii]=tmpV3.phi(); 
						gHRSTree->TA1_Pout[ii]=aHit->GetOutMom().mag()/GeV;  

						gHRSTree->TA1_N++; 
					}
				}
				else if(SDName=="thirdArmSC2")
				{
					int ii=gHRSTree->TA2_N;
					if(ii<MaxSDHit)
					{
						gHRSTree->TA2_Pid[ii]=aHit->GetPdgid(); 
						gHRSTree->TA2_Tid[ii]=aHit->GetTrackId(); 
						gHRSTree->TA2_ParentTid[ii]=aHit->GetParentTrackId(); 
						gHRSTree->TA2_T[ii]=aHit->GetTime()/ns; 
						tmpV3=aHit->GetInPos();
						gHRSTree->TA2_X[ii]=tmpV3.x()/mm;
						gHRSTree->TA2_Y[ii]=tmpV3.y()/mm; 
						gHRSTree->TA2_Z[ii]=tmpV3.z()/mm; 
						gHRSTree->TA2_Edep[ii]=aHit->GetEdep()/MeV;
						gHRSTree->TA2_NonIonEdep[ii]=aHit->GetNonIonEdep()/MeV;
						tmpV3=aHit->GetInMom();
						gHRSTree->TA2_P[ii]=tmpV3.mag()/GeV;  
						gHRSTree->TA2_Theta[ii]=tmpV3.theta(); 
						gHRSTree->TA2_Phi[ii]=tmpV3.phi(); 
						gHRSTree->TA2_Pout[ii]=aHit->GetOutMom().mag()/GeV; 

						gHRSTree->TA2_N++; 
					}
				} //end of TA2
				else
				{
					int ii=gHRSTree->SD_N;
					if(ii<MaxSDHit)
					{
						gHRSTree->SD_Id[ii]=theSDID*1.0E06+aHit->GetId();   //one hit one track
						gHRSTree->SD_Pid[ii]=aHit->GetPdgid(); 
						gHRSTree->SD_Tid[ii]=aHit->GetTrackId(); 
						gHRSTree->SD_ParentTid[ii]=aHit->GetParentTrackId(); 
						gHRSTree->SD_T[ii]=aHit->GetTime()/ns; 
						tmpV3=aHit->GetInPos();
						gHRSTree->SD_X[ii]=tmpV3.x()/mm;
						gHRSTree->SD_Y[ii]=tmpV3.y()/mm; 
						gHRSTree->SD_Z[ii]=tmpV3.z()/mm; 
						gHRSTree->SD_Edep[ii]=aHit->GetEdep()/MeV;
						gHRSTree->SD_NonIonEdep[ii]=aHit->GetNonIonEdep()/MeV;
						tmpV3=aHit->GetInMom();
						gHRSTree->SD_P[ii]=tmpV3.mag()/GeV;  
						gHRSTree->SD_Theta[ii]=tmpV3.theta(); 
						gHRSTree->SD_Phi[ii]=tmpV3.phi(); 
						gHRSTree->SD_Pout[ii]=aHit->GetOutMom().mag()/GeV; 

						gHRSTree->SD_N++; 
					}
				}
			} //end of for ih 

		} //end of  if(colName=="StdColl")
		else if(colName=="CalorimeterColl")
		{
			theECHC = (HRSCalorimeterHitsCollection*) HCE->GetHC(ihc);
			HRSCalorimeterHit *aHit=0; 
			G4ThreeVector tmpV3;

			G4int n_hit = theECHC->entries();

			if(verboseLevel>=1) 
			{
				if(n_hit>0) G4cout<< (*theECHC)[0]->GetPhysV()->GetName()<< " has " << n_hit << " hits." << G4endl;
			}

			for(int ih=0;ih<n_hit;ih++)
			{
				aHit=(*theECHC)[ih];
				if(verboseLevel>=2) aHit->Print();

				if(aHit->GetTrackId()<MaxStoredTrjN)
				{
					mStoreTrackIdList[aHit->GetParentTrackId()]=true;
					mStoreTrackIdList[aHit->GetTrackId()]=true;
				}

				int ii=gHRSTree->SD_N;
				if(ii<MaxSDHit)
				{
					gHRSTree->SD_Id[ii]=theSDID;   //one hit one event, no need to multiple 1.0E6
					gHRSTree->SD_Pid[ii]=aHit->GetPdgid(); 
					gHRSTree->SD_Tid[ii]=aHit->GetTrackId(); 
					gHRSTree->SD_ParentTid[ii]=aHit->GetParentTrackId(); 
					gHRSTree->SD_T[ii]=aHit->GetTime()/ns; 
					tmpV3=aHit->GetInPos();
					gHRSTree->SD_X[ii]=tmpV3.x()/mm;
					gHRSTree->SD_Y[ii]=tmpV3.y()/mm; 
					gHRSTree->SD_Z[ii]=tmpV3.z()/mm; 
					gHRSTree->SD_Edep[ii]=aHit->GetEdep()/MeV;
					gHRSTree->SD_NonIonEdep[ii]=aHit->GetNonIonEdep()/MeV;
					tmpV3=aHit->GetInMom();
					gHRSTree->SD_P[ii]=tmpV3.mag()/GeV;  
					gHRSTree->SD_Theta[ii]=tmpV3.theta(); 
					gHRSTree->SD_Phi[ii]=tmpV3.phi(); 
					gHRSTree->SD_Pout[ii]=aHit->GetOutMom().mag()/GeV; 

					gHRSTree->SD_N++; 
				}
			} //end of for ih 

		} //end of  if(colName=="CalorimeterColl")

	} //end of ihc
}


void HRSEventAction::EndOfEventAction(const G4Event* evt)
{
	if(!gHRSTree->iNoDetectorResponse)
	{
		UsageManager* pConfig=UsageManager::GetUsageManager(); 

		///////////////////////////////////////////////////////////////
		G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
		
		if (HCE) ProcessSDHits(HCE);

		///////////////////////////////////////////////////////////////////////
		// get number of stored trajectories
		G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
		G4int n_trajectories = 0;
		if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
		// extract the trajectories and store them into root ntuple
		G4ThreeVector tmpV3;
		HRSTrajectory* aTrj=0;
		//////////////////////////
		int pStoreSecondary=0;
		pConfig->GetArgument("StoreSecondary",pStoreSecondary);
		//0 mean do not store, 1 means store all , 2 means store only those tracks who fired SD 
		//primary tracks will always be stored
		int pStoreTrajectory=0;
		pConfig->GetArgument("StoreTrajectory",pStoreTrajectory);

		//////////////////////////////////////////////////////////
		//Jixie: The G4 store traj in a tree structure: first come last out
		//It processes one primary after another. In each primary,the trajectories
		//are stored according to their created time.
		//I have to put the order in the D tree right

		for(G4int i=n_trajectories-1; i>=0; i--)
		{
			aTrj = (HRSTrajectory*)((*(evt->GetTrajectoryContainer()))[i]);

			//store only those tracks have hits and their parent tracks, 
			//all of them have been flaged in array mStoreTrackIdList[]
			
			//make sure number of flaged tracks will not exceed the 
			//maximum buffer mStoreTrackIdList 
			//Note that the TrajectoryContainer does not keep traj in the trackid order
			int trackid = aTrj->GetTrackID();
			int parentid = aTrj->GetParentID();
			if(trackid >= MaxStoredTrjN) continue;

			if( parentid==0 || pStoreSecondary==1 ||
				(pStoreSecondary>1 && mStoreTrackIdList[trackid]) )
			{
				//store a traj into root tree
				//make sure number of flaged tracks will not exceed the 
				//maximum buffer of the tree
				int ii=gHRSTree->T_N;
				if(ii>=MaxParticle) break;

				gHRSTree->T_Pid[ii]=aTrj->GetPDGEncoding(); 
				gHRSTree->T_Tid[ii]=trackid; 
				gHRSTree->T_ParentTid[ii]=parentid; 
				gHRSTree->T_T[ii]=aTrj->GetInitialTime()/ns; 
				tmpV3=((G4TrajectoryPoint*)(aTrj->GetPoint(0)))->GetPosition();
				gHRSTree->T_X[ii]=tmpV3.x()/mm;
				gHRSTree->T_Y[ii]=tmpV3.y()/mm; 
				gHRSTree->T_Z[ii]=tmpV3.z()/mm; 
				tmpV3=aTrj->GetInitialMomentum();
				gHRSTree->T_P[ii]=tmpV3.mag()/GeV;  
				gHRSTree->T_Theta[ii]=tmpV3.theta(); 
				gHRSTree->T_Phi[ii]=tmpV3.phi(); 

				if(pStoreTrajectory)
				{
					gHRSTree->T_StepN[ii]=aTrj->GetPointEntries(); 
					for(int ss=0;ss<aTrj->GetPointEntries();ss++)
					{
						if(ss>=MaxTrackHit) break;
						tmpV3=((G4TrajectoryPoint*)(aTrj->GetPoint(ss)))->GetPosition();
						gHRSTree->T_StepX[ii][ss]=tmpV3.x()/mm;
						gHRSTree->T_StepY[ii][ss]=tmpV3.y()/mm; 
						gHRSTree->T_StepZ[ii][ss]=tmpV3.z()/mm; 
					}
				}

				gHRSTree->T_N++; 
			}
			//finish storing a traj into root tree

			//print the trajectories out if verbose>=3
			if (verboseLevel>=3)
			{
				aTrj->ShowTrajectory();
			}
		}

		//print primary particle info if verbose>=4
		if (verboseLevel>=4)
		{
			G4cout << G4endl;
			G4cout << "Primary particles -------------------------------------------------" << G4endl;
			G4int n_vertex = evt->GetNumberOfPrimaryVertex();
			for(G4int iv=0;iv<n_vertex;iv++)
			{
				G4PrimaryVertex* pv = evt->GetPrimaryVertex(iv);
				G4cout << G4endl;
				G4cout << "Primary vertex = ("
					<< pv->GetX0()/mm<<", "<<pv->GetY0()<<", "<<pv->GetZ0()
					<< ") [mm]  at time = " << (pv->GetT0())/ns << " [ns]" << G4endl;
				G4PrimaryParticle* pp = pv->GetPrimary();
				while(pp)
				{
					PrintPrimary(pp,0);
					pp = pp->GetNext();
				}
			}
		}

	}  //end of NoDetectorResponse


	//fill SD table into config tree
	if(!gHRSTree->GetConfigTreeFilledFlag())
	{		
		G4SDManager* SDman = G4SDManager::GetSDMpointer();
		G4HCtable* HCtable = SDman->GetHCtable();
		int nHC = HCtable->entries();
		gHRSTree->SetSDNum(nHC);
		for(int ihc=0;ihc<nHC;ihc++)
		{
			gHRSTree->SetSDName(ihc,(HCtable->GetSDname(ihc)).c_str());
		}
	}

	//Do root tree stuff
	gHRSTree->DoRootTree();
}


void HRSEventAction::PrintPrimary(G4PrimaryParticle* pp,G4int ind)
{
	for(G4int ii=0;ii<=ind;ii++) { G4cout << "  "; }
	G4cout << " PDGcode=" << pp->GetPDGcode() << " ";
	if(pp->GetG4code()!=0) { G4cout << "(" << pp->GetG4code()->GetParticleName() << ")"; }
	else { G4cout << "is not defined in G4"; }
	G4cout << " Momentum=" << pp->GetMomentum()/GeV << " [GeV] ";
	if(pp->GetTrackID()<0) { G4cout << G4endl; }
	else { G4cout << ">>> TrackID=" << pp->GetTrackID() << G4endl; }

	G4PrimaryParticle* daughter = pp->GetDaughter();
	while(daughter)
	{
		PrintPrimary(daughter,ind+1);
		daughter = daughter->GetNext();
	}
}

