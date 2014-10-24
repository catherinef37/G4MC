// ********************************************************************
//
// $Id: HRSSteppingAction.cc,v3.1 2008/3/16 HRS Exp $
//
//..............................................................................
#include <iomanip>
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4String.hh"
#include "HRSSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "HRSSteppingActionMessenger.hh"
#include "HRSTrackInformation.hh"
#include "HRSRootTree.hh"
#include "UsageManager.hh"
//..............................................................................

//#define STEPPING_DEBUG 2
extern UsageManager* gConfig;
extern HRSRootTree* gHRSTree;
extern bool CreateFileName(char *instr,char *outstr,int type);
HRSSteppingAction::HRSSteppingAction()
{
	PrintHeadFlag=true;
	verboseLevel=2;
	vPrintList.push_back("all"); 

	int pBookTrees=0, pBookHistos=0, pBookTxt=0;
	gConfig->GetParameter("BookTrees",pBookTrees);
	gConfig->GetParameter("BookHistos",pBookHistos);
	gConfig->GetParameter("BookTxt",pBookTxt);

	CreateTxt=(pBookTxt>0);
	CreateRootNt=(pBookTrees>0 || pBookHistos>0);

	messenger = new HRSSteppingActionMessenger(this);

	G4cout<<"HRSSteppingAction() construction done!"<<G4endl;
}

//..............................................................................

HRSSteppingAction::~HRSSteppingAction()
{
	if (CreateTxt)
	{
		OutTxt.close();
	}
	vPrintList.clear();
	delete messenger;
	G4cout<<"delete HRSSteppingAction ... done!"<<G4endl;
}
//..............................................................................

void HRSSteppingAction::InitOutTxt()
{
	if (!CreateTxt) return;

	char strFileName[200], strRawName[200];
	
	std::string pOutFileName=gConfig->GetArgument("OutFileName");
	if(gHRSTree)
	{
		pOutFileName=gHRSTree->GetRootFileName();
	}
	if(pOutFileName.length()>3) 
	{
		string tmpStr(pOutFileName);
		size_t pos=tmpStr.rfind(".root");
		if(pos!=string::npos) tmpStr.replace(pos,5,".txt"); 
		strcpy(strFileName,tmpStr.c_str());
	}
	else
	{
		//sprintf(strRawName,"G4Sim_txt_%02d.txt",gHRSTree->iRunNumber);
		//CreateFileName(strRawName,strFileName);
		//in order to make the txt output has the same name with the root output, I check 
		//the existance of the ntuple root file instead of txt file
		sprintf(strRawName,"G4Sim_nt_%02d.root",gHRSTree->iRunNumber);
		CreateFileName(strRawName,strFileName,2);
	}
	OutTxt.open (strFileName, fstream::out | fstream::app);
}

void HRSSteppingAction::UserSteppingAction(const G4Step* theStep)
{
	//control stuff
	G4Track * theTrack = theStep->GetTrack();

	//max step num 1024, if more than this number, do not pass to root tree
	if(theTrack->GetCurrentStepNumber()<MaxStepPerTrack)
	{
		if(CreateRootNt) FillRootArray(theStep);
	}
	else
	{
		theTrack->SetTrackStatus(fStopAndKill);
		PrintHeadFlag=true;
		return;
	}  
	
	///////////////////////////////////////////////////////////////////////
	// check the status from HRSTrackInformation, 
	HRSTrackInformation *anInfo=(HRSTrackInformation *)(theTrack->GetUserInformation());
	G4int status=-5;
	status=anInfo->GetStepOutputStatus();
	//G4cout<<"*******debug******HRSSteppingAction::GetStepOutputStatus() status="<<status<<"***\n";
	if(status==-1)
	{
		//this part take care of the 2ndary particles already
		//When NoSecondary==1, all stutus of 2ndaries will be set to -1 in TrackingAction
		theTrack->SetTrackStatus(fStopAndKill);
		PrintHeadFlag=true;
		return;
	}
	else if(status==0) return;

	///////////////////////////////////////////////////////////////////////
	if(verboseLevel>=6)
	{
		//only print those in the print list
		G4String thePhysName=theTrack->GetVolume()->GetName();
		
		bool bDoPrint=false;
		for(size_t ii=0;ii<vPrintList.size();ii++)
		{
			if( thePhysName == vPrintList[ii] || vPrintList[ii]=="all" 
				|| vPrintList[ii]=="All" || vPrintList[ii]=="ALL")
			{
				bDoPrint=true;
				break;
			}
		}
		if(!PrintHeadFlag && !bDoPrint) return;
	}
	DoPrint(theStep);

	
	///////////////////////////////////////////////////////////////////////
	//check if it is alive, FillRootArray() will set this track stutus as fStopAndKill
	//I have to put this after calling FillRootArray() otherwise the intial variables of some events 
	//will not be set 
	if(theTrack->GetTrackStatus()!=fAlive) 
	{
		PrintHeadFlag=true;
		return;
	}

	///////////////////////////////////////////////////////////////////////
	//determine if this track need to be killed or not
	G4double xx=theTrack->GetPosition().x()/mm;
	G4double yy=theTrack->GetPosition().y()/mm;
	//G4double zz=theTrack->GetPosition().z()/mm;

	/*
	//By Jixie: 
	verboseLevel description: Maximum step=1024
	//kill this track whenever it hits the virtual boundary:  r=116mm fabs(z)<151mm
	Level 1:  print all except VolumnName and ProcessName, 
	Level 2:  if verboseLevel>1, VolumnName and ProcessName will be printed
	Level 3:  print all; track will be killed if r> 5000.0mm
	Level 4:  print all; track will be killed if r> 2500mm
	Level 5:  print all; track will be killed if r> 800mm
	Level 6:  print all; but only print in those physicalVolumn in the printlist
	*/

	if(verboseLevel==3)
	{
		if (sqrt(xx*xx+yy*yy)>=5000.0*mm)
		{
			theTrack->SetTrackStatus(fStopAndKill);
			PrintHeadFlag=true;
			return;
		}
	}
	else if(verboseLevel==4)
	{
		if (sqrt(xx*xx+yy*yy)>=2500.0*mm)
		{
			theTrack->SetTrackStatus(fStopAndKill);
			PrintHeadFlag=true;
			return;
		}
	}
	else if(verboseLevel==5)
	{
		if (sqrt(xx*xx+yy*yy)>=800.0*mm)
		{
			theTrack->SetTrackStatus(fStopAndKill);
			PrintHeadFlag=true;
			return;
		}
	}
	
}

//..............................................................................

void HRSSteppingAction::DoPrint(const G4Step* theStep)
{
	if(CreateTxt)
	{
		//check if this file is open, if not just open it
		if(! OutTxt.is_open()) this->InitOutTxt();
		if(PrintHeadFlag)
		{
			PrintHead(theStep,OutTxt);
			if(verboseLevel>0) PrintHeadFlag=true;
		}
		PrintStep(theStep,OutTxt);
	}

	if(verboseLevel>0)
	{
		if(PrintHeadFlag) PrintHead(theStep,G4cout);
		PrintStep(theStep,G4cout);
	}
}


//..............................................................................

void HRSSteppingAction::PrintHead(const G4Step* theStep,ostream& pOut)
{
	G4Track * theTrack = theStep->GetTrack();
	if(theTrack->GetCurrentStepNumber()>1) return;

	G4int prec = pOut.precision(4);
	pOut.setf(ios::fixed);
	pOut<< G4endl;
	pOut<< "*******************************************************"
		<< "*********************************************************"<< G4endl;
	pOut << "* G4Track Information: "
		<< "  Particle = "  << theTrack->GetDefinition()->GetParticleName()
		<< ",  Track ID = " << theTrack->GetTrackID()
		<< ",  Parent ID = "<< theTrack->GetParentID()
		<< ",  Momentum = " << theStep->GetPreStepPoint()->GetMomentum().mag()/MeV<< "MeV "
		<< G4endl;
	pOut<< "*******************************************************"
		<< "*********************************************************"<< G4endl;
	pOut<< std::setw( 4) << "Step#"      << " "
		<< std::setw( 7) << "X(mm)"      << " "		//shift to left by 1
		<< std::setw( 8) << "Y(mm)"      << " "
		<< std::setw( 8) << "Z(mm)"      << " "
		<< std::setw( 9) << "Ekin(MeV)"  << " "		
		<< std::setw( 8) << "dE(KeV)"    << " "		//shift to left by 1
		<< std::setw( 9) << "Mom(MeV)"   << " " 
		<< std::setw( 7) << "Theta(deg)" << " " 
		<< std::setw( 6) << "Phi(deg)"   << " "		//shift to left by 1
		<< std::setw( 9) << "StepL(mm)"  << " " 
		<< std::setw( 9) << "TrackL(mm)" << " "
		<< std::setw( 9) << "Radlen(m)"  << " ";	//shift to left by 2
	if(verboseLevel>=2)
	{
		pOut<< std::setw(12) << "ProcName" << " "   //shift to left by 2
			<< std::setw(18) << "NextVolume";
	}
	pOut<< G4endl;

	//print the initial step
	G4StepPoint* thePrePoint = theStep->GetPreStepPoint();

	pOut.precision(2);
	pOut<< std::setw(4) << 0 << " "
		<< std::setw(8) << thePrePoint->GetPosition().x()/mm   << " "
		<< std::setw(8) << thePrePoint->GetPosition().y()/mm   << " "
		<< std::setw(8) << thePrePoint->GetPosition().z()/mm   << " ";
	pOut.precision(3);
	pOut<< std::setw(9) << thePrePoint->GetKineticEnergy()/MeV << " "
		<< std::setw(9) << 0 << " "
		<< std::setw(9) << thePrePoint->GetMomentum().mag()/MeV << " ";
	pOut.precision(2);
	pOut<< std::setw(7) << thePrePoint->GetMomentum().theta()/deg << " "
		<< std::setw(7) << thePrePoint->GetMomentum().phi()/deg << " "
		<< std::setw(9) << 0 << " "
		<< std::setw(9) << 0 << " ";

	pOut.precision(4); 
	if(theTrack->GetMaterial()->GetRadlen()/m > 10000.0)
	{
		pOut.unsetf(ios::fixed);
		pOut.setf(ios::scientific); 
		pOut<< std::setw(11) << theTrack->GetMaterial()->GetRadlen()/m<<" ";
		pOut.unsetf(ios::scientific);
	}
	else
	{		
		pOut<< std::setw(11) << theTrack->GetMaterial()->GetRadlen()/m<<" ";
	}

	if(verboseLevel>=2)
	{
		pOut<< std::setw(14)<< "initStep"<<"  "
			<< std::setw(18)<< thePrePoint->GetPhysicalVolume()->GetName();
	}
	pOut << G4endl;
	pOut.precision(prec);

	PrintHeadFlag=false;
}
//..............................................................................

void HRSSteppingAction::PrintStep(const G4Step* theStep,ostream& pOut)
{	
	G4Track * theTrack = theStep->GetTrack();
	pOut.setf( ios::fixed ); 

	G4int prec = pOut.precision(2);
	pOut<< std::setw(4) << theTrack->GetCurrentStepNumber() << " "
		<< std::setw(8) << theTrack->GetPosition().x()/mm<< " "
		<< std::setw(8) << theTrack->GetPosition().y()/mm<< " "
		<< std::setw(8) << theTrack->GetPosition().z()/mm<< " ";

	pOut.precision(3); 
	pOut<< std::setw(9) << theTrack->GetKineticEnergy()/MeV<< " "
		<< std::setw(9) << theStep->GetTotalEnergyDeposit()/keV<< " "
		<< std::setw(9) << theTrack->GetMomentum().mag()/MeV<< " ";

	pOut.precision(2); 
	pOut<< std::setw(7) << theTrack->GetMomentum().theta()/deg<< " "
		<< std::setw(7) << theTrack->GetMomentum().phi()/deg<< " "
		<< std::setw(9) << theStep->GetStepLength()/mm<< " "
		<< std::setw(9) << theTrack->GetTrackLength()/mm<< " ";

	//some material has very large radlen, use scientific notation for it if larger than 10^5
	pOut.precision(4); 
	if(theTrack->GetMaterial()->GetRadlen()/m > 10000.0)
	{
		pOut.unsetf(ios::fixed);
		pOut.setf(ios::scientific); 
		pOut<< std::setw(11) << theTrack->GetMaterial()->GetRadlen()/m<<" ";
		pOut.unsetf(ios::scientific);
	}
	else
	{		
		pOut<< std::setw(11) << theTrack->GetMaterial()->GetRadlen()/m<<" ";
	}

	if(verboseLevel>=2)
	{	
		if(theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL)
		{
			pOut<< std::setw(14)
				<< theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
				<< "  "; 
		}
		else
		{
			pOut<< std::setw(14)<<"  userLimit"<<"  ";
		}

		if( theTrack->GetNextVolume() != 0 )
		{
			pOut << std::setw(18) << theTrack->GetVolume()->GetName();
		}
		else
		{
			pOut << std::setw(18) << "OutOfWorld";
		}

	}
	pOut<< G4endl;
	pOut.precision(prec);

}
//..............................................................................

void HRSSteppingAction::FillRootArray(const G4Step* theStep)
{
	G4Track * theTrack = theStep->GetTrack();
	
	string strPhysVolName = theTrack->GetVolume()->GetName();
	//for sencondary tracks, do not record any steps, kill it
	//if it touch vb or absorber
	if(theTrack->GetParentID()!=0) 
	{		
		if ( strPhysVolName.find("virtualBoundaryPhys") != string::npos || 
			strPhysVolName.find("absorberPhys") != string::npos )
		{
			//hit the absorber, kill this track
			theTrack->SetTrackStatus(fStopAndKill);
			PrintHeadFlag=true;
		}
		return;
	}

	int trackid=theTrack->GetTrackID();
	if(trackid>MaxPrimaryNum) return;

	//FILL THE VERTEX, only at the first step
	G4StepPoint* thePrePoint=theStep->GetPreStepPoint();
	MyTrack* RootTrack = gHRSTree->track[trackid-1];
	if(theTrack->GetCurrentStepNumber()<=1)
	{
		//cout<<"Event "<<gHRSTree->GetTotalEvtID()<<": " 
		//	<<"HRSSteppingAction::FillRootArray() is called to initialize the 1st step."<<endl;

		RootTrack->PdgId=theTrack->GetDefinition()->GetPDGEncoding();

		RootTrack->TrackId=trackid;
		RootTrack->X0=thePrePoint->GetPosition().x()/mm;
		RootTrack->Y0=thePrePoint->GetPosition().y()/mm;
		RootTrack->Z0=thePrePoint->GetPosition().z()/mm;
		RootTrack->P0=thePrePoint->GetMomentum().mag()/GeV;
		RootTrack->Theta0=thePrePoint->GetMomentum().getTheta()/rad;
		RootTrack->Phi0=thePrePoint->GetMomentum().getPhi()/rad;

#ifdef STEPPING_DEBUG
		if(STEPPING_DEBUG>=1)
		{
			G4cout<<" Debug: trackid="<<trackid<<" Initial step: particle="
				<<theTrack->GetDefinition()->GetParticleName()
				<<" PDGEncoding="<<RootTrack->PdgId<<G4endl;
		}
#endif

		gHRSTree->material[trackid-1] = theTrack->GetLogicalVolumeAtVertex()->GetMaterial();		
	}

	G4double xx=theTrack->GetPosition().x()/mm;
	G4double yy=theTrack->GetPosition().y()/mm;
	G4double zz=theTrack->GetPosition().z()/mm;
	G4double ekin=theTrack->GetKineticEnergy()/GeV;
	G4double de=theStep->GetTotalEnergyDeposit()/keV;
	G4double dl=theStep->GetStepLength()/mm;
	G4double tl=theTrack->GetTrackLength()/mm;
	G4ThreeVector v3p=theTrack->GetMomentum();

	//fill virtual boundary information
	std::size_t foundVB = strPhysVolName.find("virtualBoundaryPhys");
	if ( foundVB != string::npos )
	{
		//hit the virtualBoundary, record this step then kill this track
		RootTrack->Xvb=xx;
		RootTrack->Yvb=yy;
		RootTrack->Zvb=zz;
		RootTrack->Pvb=v3p.mag()/GeV;
		RootTrack->Thetavb=v3p.theta()/rad;
		RootTrack->Phivb=v3p.phi()/rad;
		RootTrack->VBName=strPhysVolName.substr(foundVB);
		//kill this track whenever it hit the virtual boundary
		theTrack->SetTrackStatus(fStopAndKill);
		PrintHeadFlag=true;
	}
	else if ( strPhysVolName.find("absorberPhys") != string::npos )
	{
		//hit the absorber, kill this track
		theTrack->SetTrackStatus(fStopAndKill);
		PrintHeadFlag=true;
	}

	int i=RootTrack->StepNum;
	if(i>=MaxStepPerTrack) return;
	RootTrack->StepX[i]=xx;
	RootTrack->StepY[i]=yy;
	RootTrack->StepZ[i]=zz;
	RootTrack->StepEkin[i] = ekin;
	RootTrack->StepdE[i] = de;
	RootTrack->StepL[i] = dl;
	
	double tmpField[6]={0,0,0,0,0,0};
	//By Jixie: In case there is a local field, this is the correct way
	G4FieldManager *theFieldManager = theTrack->GetVolume()->GetLogicalVolume()->GetFieldManager();   
	if(!theFieldManager) theFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
	if(theFieldManager)
	{
		double tmpPoint[4]={xx*mm,yy*mm,zz*mm,0};
		theFieldManager->GetDetectorField()->GetFieldValue(tmpPoint,tmpField);
	}
	RootTrack->StepBx[i] = tmpField[0]/tesla;
	RootTrack->StepBy[i] = tmpField[1]/tesla;
	RootTrack->StepBz[i] = tmpField[2]/tesla;
	
	G4ThreeVector VL = v3p.unit() * theStep->GetStepLength()/m;
	G4ThreeVector VB(tmpField[0]/tesla,tmpField[1]/tesla,tmpField[2]/tesla);
	G4ThreeVector VF = VL.cross(VB);
	RootTrack->TrackBdLx += VF.x(); 
	RootTrack->TrackBdLy += VF.y(); 
	RootTrack->TrackBdLz += VF.z(); 

	RootTrack->StepNum +=1;

	//add radiation length and density of each step into the ntuple 
	// Radiation length:     00191   G4double         GetRadlen()    
	G4double pRadlen= theTrack->GetMaterial()->GetRadlen()/mm;
	G4double pDensity= theTrack->GetMaterial()->GetDensity()/(g/cm3);
	RootTrack->StepTL[i] = tl;
	RootTrack->StepRadlen[i] = pRadlen;
	RootTrack->StepDsty[i] = pDensity;
	//total radiation length for this track
	if(pRadlen!=0)  RootTrack->TrackRadlen += dl/pRadlen;
}
//..............................................................................

