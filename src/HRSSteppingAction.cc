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
#include "HRSTransform_TCSNHCS.hh"
//..............................................................................

//#define STEPPING_DEBUG 2
extern UsageManager* gConfig;
extern HRSRootTree* gHRSTree;
extern bool CreateFileName(char *instr,char *outstr,int type);
HRSSteppingAction::HRSSteppingAction()
{
	PrintHeadFlag=true;
	//verboseLevel=2;
	verboseLevel=0;
	vPrintList.push_back("all"); 

	int pBookTrees=0, pBookHistos=0, pBookTxt=0;
	gConfig->GetParameter("BookTrees",pBookTrees);
	gConfig->GetParameter("BookHistos",pBookHistos);
	gConfig->GetParameter("BookTxt",pBookTxt);

	CreateTxt=(pBookTxt>0);
	CreateRootNt=(pBookTrees>0 || pBookHistos>0);

	messenger = new HRSSteppingActionMessenger(this);

	//G4cout<<"HRSSteppingAction() construction done!"<<G4endl;
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
	//G4cout<<"delete HRSSteppingAction ... done!"<<G4endl;
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
  //cout << "WE ARE IN THE USER STEPPING ACTION DOOOOOOOOOOOOOOOOOOOOP!" << endl;
	//control stuff
	G4Track * theTrack = theStep->GetTrack();

	//max step num 1024, if more than this number, do not pass to root tree
	if(theTrack->GetCurrentStepNumber()<MaxStepPerTrack)
	{
		if(CreateRootNt) FillRootArray(theStep);
	}
	else
	{
	  //cout << "killing place 1" << endl;
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
	  //cout << "killing place 2" << endl;
	  theTrack->SetTrackStatus(fStopAndKill);
	  PrintHeadFlag=true;
	  return;
	}
	else if(status==0) return;
	//cout << "We have passed the status check." << endl;
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

	//cout << "THE ";
	
	///////////////////////////////////////////////////////////////////////
	//check if it is alive, FillRootArray() will set this track stutus as fStopAndKill
	//I have to put this after calling FillRootArray() otherwise the intial variables of some events 
	//will not be set 
	if(theTrack->GetTrackStatus()!=fAlive) 
	{
	  //cout << "The track is not alive." << endl;
		PrintHeadFlag=true;
		return;
	}
	//cout << "QUICK ";
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
	//cout << "Vervose level: " << verboseLevel << endl;
	if(verboseLevel==3)
	{
		if (sqrt(xx*xx+yy*yy)>=5000.0*mm)
		{
		  //cout << "killing place 3" << endl;
			theTrack->SetTrackStatus(fStopAndKill);
			PrintHeadFlag=true;
			return;
		}
	}
	else if(verboseLevel==4)
	{
		if (sqrt(xx*xx+yy*yy)>=2500.0*mm)
		{
		  //cout << "killing place 4" << endl;
			theTrack->SetTrackStatus(fStopAndKill);
			PrintHeadFlag=true;
			return;
		}
	}
	else if(verboseLevel==5)
	{
		if (sqrt(xx*xx+yy*yy)>=800.0*mm)
		{
		  //cout << "killing place 5" << endl;
			theTrack->SetTrackStatus(fStopAndKill);
			PrintHeadFlag=true;
			return;
		}
	}
	//cout << "FOX!" << endl;
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
		<< std::setw( 7) << "B_X(?)"     << " "		//shift to left by 1
		<< std::setw( 8) << "B_Y(?)"     << " "
		<< std::setw( 8) << "B_Z(?)"     << " "
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
	    << std::setw(8) << thePrePoint->GetPosition().z()/mm   << " "
	    << std::setw(8) << 0.0                                 << " "
	    << std::setw(8) << 0.0                                 << " "
	    << std::setw(8) << 0.0                                 << " ";
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

	int trackid=theTrack->GetTrackID();
	MyTrack* RootTrack = gHRSTree->track[trackid-1];
	int i=RootTrack->StepNum;
	G4int prec = pOut.precision(2);
	pOut<< std::setw(4) << theTrack->GetCurrentStepNumber() << " "
	    << std::setw(8) << theTrack->GetPosition().x()/mm<< " "
	    << std::setw(8) << theTrack->GetPosition().y()/mm<< " "
	    << std::setw(8) << theTrack->GetPosition().z()/mm<< " ";
	pOut.precision(3); 
	pOut<< std::setw(8) << RootTrack->StepBx[i-1]/mm<< " "
	    << std::setw(8) << RootTrack->StepBy[i-1]/mm<< " "
	    << std::setw(8) << RootTrack->StepBz[i-1]/mm<< " ";

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
  G4double mLHRSAngle;
  G4double mRHRSAngle;
  gConfig->GetParameter("LHRSAngle",mLHRSAngle);         //deg
  mLHRSAngle*=deg;
  gConfig->GetParameter("RHRSAngle",mRHRSAngle);         //deg
  //  G4cout << mRHRSAngle << G4endl;
  mRHRSAngle*=deg;
  //  G4cout << mRHRSAngle << G4endl;

  G4Track * theTrack = theStep->GetTrack();
  
	string strPhysVolName = theTrack->GetVolume()->GetName();
	//for sencondary tracks, do not record any steps, kill it
	//if it touch vb or absorber
	if(theTrack->GetParentID()!=0) 
	{		
		if ( strPhysVolName.find("virtualBoundaryPhys") != string::npos || 
			strPhysVolName.find("absorberPhys") != string::npos )
		{
		  cout << "killing place 6" << endl;
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
		RootTrack->X0=thePrePoint->GetPosition().x()/mm / 1000.;//force to be in meters
		RootTrack->Y0=thePrePoint->GetPosition().y()/mm / 1000.;
		RootTrack->Z0=thePrePoint->GetPosition().z()/mm / 1000.;
		RootTrack->P0=thePrePoint->GetMomentum().mag()/GeV;
		RootTrack->Theta0=thePrePoint->GetMomentum().getTheta()/rad;
		RootTrack->Phi0=thePrePoint->GetMomentum().getPhi()/rad;

		//transform directly to transport coordinates.

		G4double dir_l[3] = { cos( RootTrack->Phi0 ) * sin( RootTrack->Theta0 ),
				      sin( RootTrack->Phi0 ) * sin( RootTrack->Theta0 ),
				      cos( RootTrack->Theta0 ) };
		G4double sep_ang  = 5. * deg; // have to make +/- for L/R HRS.
		if( thePrePoint->GetMomentum().x() > 0 )
		  sep_ang *= -1.;
		G4double rot_y[3][3] = { { cos( sep_ang ), 0., sin( sep_ang )},
					 { 0.            , 1., 0.            },
					 {-sin( sep_ang ), 0., cos( sep_ang )} } ;
		G4double num_th = 0., num_ph = 0., denom = 0.;//numerators and denominator
		for( G4int i = 0; i < 3; i++ ){
		  num_ph += rot_y[0][i] * dir_l[i];
		  num_th += rot_y[1][i] * dir_l[i];
		  denom  += rot_y[2][i] * dir_l[i];
		}
		RootTrack->Theta0_tr = -num_th / denom;
		RootTrack->Phi0_tr   =  num_ph / denom;

		RootTrack->X0_tr     = -RootTrack->Y0;
		RootTrack->Y0_tr     =  RootTrack->X0 * cos( sep_ang );
		RootTrack->Z0_tr     =  RootTrack->Z0;

		/*
		G4cout << "Hall coordinates at target"                                                 << G4endl;
		G4cout << RootTrack->X0        << " " << RootTrack->Y0      << " " << RootTrack->Z0    << G4endl;
		G4cout << RootTrack->Theta0    << " " << RootTrack->Phi0    <<                            G4endl;
		G4cout << "Trns coordinates at target"                                                 << G4endl;
		G4cout << RootTrack->X0_tr     << " " << RootTrack->Y0_tr   << " " << RootTrack->Z0_tr << G4endl;
		G4cout << RootTrack->Theta0_tr << " " << RootTrack->Phi0_tr <<                            G4endl;
		*/

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
	//cout << strPhysVolName << endl;
	//fill virtual boundary information
	//double rot_theta = 12.5 * deg;
	double lrot_theta = -mLHRSAngle;
	double rrot_theta = -mRHRSAngle;
	
	double rot_phi   = 45.0 * deg;

	double lrot_mat_y[3][3] = { { cos(lrot_theta), 0.0,           sin(lrot_theta)},
				    { 0.0,             1.0,           0.0            },
				    {-sin(lrot_theta), 0.0,           cos(lrot_theta)} };
	double rrot_mat_y[3][3] = { { cos(rrot_theta), 0.0,           sin(rrot_theta)},
				    { 0.0,             1.0,           0.0            },
				    {-sin(rrot_theta), 0.0,           cos(rrot_theta)} };
	double rot_mat_x[3][3]  = { { 1.0,             0.0,           0.0            },
				    { 0.0,             cos(rot_phi), -sin(rot_phi)   },
				    { 0.0,             sin(rot_phi),  cos(rot_phi)   } };

	G4double pQ3eny =  3.5853101 * m;
        G4double pQ3enz = 17.0257042 * m;
	G4double pdexy  =  2.4603032 * m;
        G4double pdexz  = 15.9006973 * m;

	//double subtract1[3]    = {0.0, 0.0, 9954.7 + 8400. * tan( 22.5 * deg ) };//in mm,  SEAMUS
	//double subtract1[3]    = {0.0, 0.0, 9961. + 8400. * tan( 22.5 * deg ) };//in mm,  SNAKE
	//double subtract1[3]    = {pQ3enz * sin( mRHRSAngle ), pQ3eny, pQ3enz * cos( mRHRSAngle )};//in mm,  SNAKE, but based off of Q3 instead of dipole "ELBOW"
	G4double subtract1[3]    = {0., pQ3eny, pQ3enz };//in mm,  SNAKE, but based off of Q3 instead of dipole "ELBOW"
	//double subtract1[3]    = {0.0, 0.0, 9960. + 8400. * tan( 22.5 * deg ) };//in mm,  NIM
	//double subtract2[3]    = {pdexz  * sin( mRHRSAngle ), pdexy , pdexz  * cos( mRHRSAngle )};//in mm
	G4double subtract2[3]    = {0., pdexy , pdexz  };//in mm
	//G4cout << subtract2[0] << " " << subtract2[1] << " " << subtract2[2] << " " << mRHRSAngle << G4endl;
	std::size_t foundVB1 = strPhysVolName.find("virtualBoundaryPhys_vb");
	std::size_t foundVB2 = strPhysVolName.find("virtualBoundaryPhys_LHRSQ1Win");
	std::size_t foundVB3 = strPhysVolName.find("virtualBoundaryPhys_LHRS");
	std::size_t foundVB4 = strPhysVolName.find("virtualBoundaryPhys_RHRSQ1Win");
	std::size_t foundVB5 = strPhysVolName.find("virtualBoundaryPhys_RHRS");
	std::size_t foundVB = strPhysVolName.find("");
	if ( (foundVB = foundVB1)  != string::npos ||
	     (foundVB = foundVB2) != string::npos ||
	     (foundVB = foundVB3) != string::npos ||
	     (foundVB = foundVB4) != string::npos ||
	     (foundVB = foundVB5) != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;

	  RootTrack->X_vb=xx;
	  RootTrack->Y_vb=yy;
	  RootTrack->Z_vb=zz;
	  RootTrack->P_vb=v3p.mag()/GeV;
	  RootTrack->Theta_vb=v3p.theta()/rad;
	  RootTrack->Phi_vb=v3p.phi()/rad;
	  RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	  int mSnakeModel;
	  gConfig->GetArgument("SnakeModel",mSnakeModel);
	  if(mSnakeModel != 49 && mSnakeModel < 51){
	    theTrack->SetTrackStatus(fStopAndKill);
	    PrintHeadFlag=true;
	  }
	}else if ( strPhysVolName.find("virtualBoundaryPhys_sen") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_sen=xx;
	  RootTrack->Y_sen=yy;
	  RootTrack->Z_sen=zz;
	  RootTrack->Theta_sen=v3p.theta()/rad;
	  RootTrack->Phi_sen=v3p.phi()/rad;
	  RootTrack->X_sen_tr=xx / 1000.;
	  RootTrack->Y_sen_tr=yy / 1000.;
	  RootTrack->Z_sen_tr=zz / 1000.;
	  RootTrack->P_sen_tr=v3p.mag()/GeV;
	  //Transform::P_HCS2TCS(RootTrack->Theta_sen, RootTrack->Phi_sen, -( RootTrack->X_sen < 0. ? -5. * deg : 5. * deg ),
	  //RootTrack->Theta_sen_tr , RootTrack->Phi_sen_tr);
	  RootTrack->Theta_sen_tr = -RootTrack->Theta_sen_tr;
	  RootTrack->Phi_sen_tr = -RootTrack->Phi_sen_tr;
	  //RootTrack->Theta_sen_tr=v3p.theta()/rad;
	  //RootTrack->Phi_sen_tr=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	}else if ( strPhysVolName.find("virtualBoundaryPhys_sm") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the sm virtual boundary!" << endl;
	  
	  RootTrack->X_sm=xx;
	  RootTrack->Y_sm=yy;
	  RootTrack->Z_sm=zz;
	  RootTrack->Theta_sm=v3p.theta()/rad;
	  RootTrack->Phi_sm=v3p.phi()/rad;
	  RootTrack->X_sm_tr=xx / 1000.;
	  RootTrack->Y_sm_tr=yy / 1000.;
	  RootTrack->Z_sm_tr=zz / 1000.;
	  RootTrack->P_sm_tr=v3p.mag()/GeV;
	  Transform::P_HCS2TCS(RootTrack->Theta_sm, RootTrack->Phi_sm, -( RootTrack->X_sm < 0. ? -5. * deg : 5. * deg ),
		    RootTrack->Theta_sm_tr , RootTrack->Phi_sm_tr);
	  RootTrack->Theta_sm_tr = -RootTrack->Theta_sm_tr;
	  RootTrack->Phi_sm_tr = -RootTrack->Phi_sm_tr;
	  //RootTrack->Theta_sm_tr=v3p.theta()/rad;
	  //RootTrack->Phi_sm_tr=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	}else if ( strPhysVolName.find("virtualBoundaryPhys_sex") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_sex=xx;
	  RootTrack->Y_sex=yy;
	  RootTrack->Z_sex=zz;
	  RootTrack->Theta_sex=v3p.theta()/rad;
	  RootTrack->Phi_sex=v3p.phi()/rad;
	  RootTrack->X_sex_tr=xx / 1000.;
	  RootTrack->Y_sex_tr=yy / 1000.;
	  RootTrack->Z_sex_tr=zz / 1000.;
	  RootTrack->P_sex_tr=v3p.mag()/GeV;
	  Transform::P_HCS2TCS(RootTrack->Theta_sex, RootTrack->Phi_sex, -( RootTrack->X_sex < 0. ? -5. * deg : 5. * deg ),
			       RootTrack->Theta_sex_tr , RootTrack->Phi_sex_tr);
	  RootTrack->Theta_sex_tr = -RootTrack->Theta_sex_tr;
	  RootTrack->Phi_sex_tr = -RootTrack->Phi_sex_tr;
	  //RootTrack->Theta_sex_tr=v3p.theta()/rad;
	  //RootTrack->Phi_sex_tr=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	}else if ( strPhysVolName.find("virtualBoundaryPhys_coil") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_coil=xx;
	  RootTrack->Y_coil=yy;
	  RootTrack->Z_coil=zz;
	  RootTrack->Theta_coil=v3p.theta()/rad;
	  RootTrack->Phi_coil=v3p.phi()/rad;
	  RootTrack->X_coil_tr=xx / 1000.;
	  RootTrack->Y_coil_tr=yy / 1000.;
	  RootTrack->Z_coil_tr=zz / 1000.;
	  RootTrack->P_coil_tr=v3p.mag()/GeV;
	  Transform::P_HCS2TCS(RootTrack->Theta_coil, RootTrack->Phi_coil, -( RootTrack->X_coil < 0. ? -5. * deg : 5. * deg ),
			       RootTrack->Theta_coil_tr , RootTrack->Phi_coil_tr);
	  RootTrack->Theta_coil_tr = -RootTrack->Theta_coil_tr;
	  RootTrack->Phi_coil_tr = -RootTrack->Phi_coil_tr;
	  //RootTrack->Theta_coil_tr=v3p.theta()/rad;
	  //RootTrack->Phi_coil_tr=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	}else if ( strPhysVolName.find("virtualBoundaryPhys_mid") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_mid=xx;
	  RootTrack->Y_mid=yy;
	  RootTrack->Z_mid=zz;
	  RootTrack->Theta_mid=v3p.theta()/rad;
	  RootTrack->Phi_mid=v3p.phi()/rad;
	  RootTrack->X_mid_tr=xx / 1000.;
	  RootTrack->Y_mid_tr=yy / 1000.;
	  RootTrack->Z_mid_tr=zz / 1000.;
	  RootTrack->P_mid_tr=v3p.mag()/GeV;
	  Transform::P_HCS2TCS(RootTrack->Theta_mid, RootTrack->Phi_mid, -( RootTrack->X_mid < 0. ? -5. * deg : 5. * deg ),
			       RootTrack->Theta_mid_tr , RootTrack->Phi_mid_tr);
	  RootTrack->Theta_mid_tr = -RootTrack->Theta_mid_tr;
	  RootTrack->Phi_mid_tr = -RootTrack->Phi_mid_tr;
	  //RootTrack->Theta_mid_tr=v3p.theta()/rad;
	  //RootTrack->Phi_mid_tr=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	}else if ( strPhysVolName.find("virtualBoundaryPhys_col") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_col=xx;
	  RootTrack->Y_col=yy;
	  RootTrack->Z_col=zz;
	  RootTrack->Theta_col=v3p.theta()/rad;
	  RootTrack->Phi_col=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_col > 0. ){
		final [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		final [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }
	    }
	  }
	  RootTrack->X_col_tr = -final[1] / 1000.;
          RootTrack->Y_col_tr =  final[0] / 1000.;
          RootTrack->Z_col_tr =  final[2] / 1000.;
	  RootTrack->P_col_tr=v3p.mag()/GeV;
	  Transform::single_rot(v3p.x(), v3p.y(), v3p.z(), -( RootTrack->X_col > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_col_tr, RootTrack->Phi_col_tr);
	  //Transform::P_HCS2TCS(RootTrack->Theta_col, RootTrack->Phi_col, ( RootTrack->X_col < 0. ? mLHRSAngle : mRHRSAngle ),
	  //RootTrack->Theta_col_tr, RootTrack->Phi_col_tr);
	  //RootTrack->Theta_col_tr = -tan(RootTrack->Theta_col_tr);
	  //RootTrack->Phi_col_tr   = -tan(RootTrack->Phi_col_tr);
	  
	}else if ( strPhysVolName.find("virtualBoundaryPhys_q1en") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_q1en=xx;
	  RootTrack->Y_q1en=yy;
	  RootTrack->Z_q1en=zz;
	  RootTrack->Theta_q1en=v3p.theta()/rad;
	  RootTrack->Phi_q1en=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_q1en > 0. ){
		final [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		final [iii] += rrot_mat_y[iii][jjj] * initial[jjj];	       
	      }
	      //cout << final[iii] << endl;
	    }
	  }
	  RootTrack->X_q1en_tr = -final[1] / 1000.;
          RootTrack->Y_q1en_tr =  final[0] / 1000.;
          RootTrack->Z_q1en_tr =  final[2] / 1000.;
	  RootTrack->P_q1en_tr=v3p.mag()/GeV;
	  //cout << "raw " << xx << " " << yy << " " << zz << endl;
	  //cout << "fin " << final[0] << " " << final[1] << " " << final[2] << endl;
	  Transform::single_rot(v3p.x(), v3p.y(), v3p.z(), -( RootTrack->X_q1en > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_q1en_tr, RootTrack->Phi_q1en_tr);
	  //Transform::P_HCS2TCS(RootTrack->Theta_q1en, RootTrack->Phi_q1en, ( RootTrack->X_q1en < 0. ? mLHRSAngle : mRHRSAngle ),
	  //RootTrack->Theta_q1en_tr, RootTrack->Phi_q1en_tr);
	  //RootTrack->Theta_q1en_tr = -tan(RootTrack->Theta_q1en_tr);
	  //RootTrack->Phi_q1en_tr   = -tan(RootTrack->Phi_q1en_tr);

	}else if ( strPhysVolName.find("virtualBoundaryPhys_q1ex") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_q1ex=xx;
	  RootTrack->Y_q1ex=yy;
	  RootTrack->Z_q1ex=zz;
	  RootTrack->Theta_q1ex=v3p.theta()/rad;
	  RootTrack->Phi_q1ex=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_q1ex > 0.) {
		final [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		final [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }
	    }
	  }
	  RootTrack->X_q1ex_tr = -final[1] / 1000.;
          RootTrack->Y_q1ex_tr =  final[0] / 1000.;
          RootTrack->Z_q1ex_tr =  final[2] / 1000.;
	  Transform::single_rot(v3p.x(), v3p.y(), v3p.z(), -( RootTrack->X_q1ex > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_q1ex_tr, RootTrack->Phi_q1ex_tr);
	  //Transform::P_HCS2TCS(RootTrack->Theta_q1ex, RootTrack->Phi_q1ex, ( RootTrack->X_q1ex < 0. ? mLHRSAngle : mRHRSAngle ),
	  //RootTrack->Theta_q1ex_tr, RootTrack->Phi_q1ex_tr);
	  //RootTrack->Theta_q1ex_tr = -tan(RootTrack->Theta_q1ex_tr);
	  //RootTrack->Phi_q1ex_tr   = -tan(RootTrack->Phi_q1ex_tr);

	}else if ( strPhysVolName.find("virtualBoundaryPhys_q2en") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_q2en=xx;
	  RootTrack->Y_q2en=yy;
	  RootTrack->Z_q2en=zz;
	  RootTrack->Theta_q2en=v3p.theta()/rad;
	  RootTrack->Phi_q2en=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_q2en > 0. ){
		final [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		final [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }
	    }
	  }
	  RootTrack->X_q2en_tr = -final[1] / 1000.;
          RootTrack->Y_q2en_tr =  final[0] / 1000.;
          RootTrack->Z_q2en_tr =  final[2] / 1000.;
	  RootTrack->P_q2en_tr =v3p.mag()/GeV;
	  Transform::single_rot(v3p.x(), v3p.y(), v3p.z(), -( RootTrack->X_q2en > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_q2en_tr, RootTrack->Phi_q2en_tr);
	  //Transform::P_HCS2TCS(RootTrack->Theta_q2en, RootTrack->Phi_q2en, ( RootTrack->X_q2en < 0. ? mLHRSAngle : mRHRSAngle ),
	  //RootTrack->Theta_q2en_tr, RootTrack->Phi_q2en_tr);
	  //RootTrack->Theta_q2en_tr = -tan(RootTrack->Theta_q2en_tr);
	  //RootTrack->Phi_q2en_tr   = -tan(RootTrack->Phi_q2en_tr);

	}else if ( strPhysVolName.find("virtualBoundaryPhys_q2ex") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_q2ex=xx;
	  RootTrack->Y_q2ex=yy;
	  RootTrack->Z_q2ex=zz;
	  RootTrack->Theta_q2ex=v3p.theta()/rad;
	  RootTrack->Phi_q2ex=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_q2en > 0. ){
		final [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		final [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }
	    }
	  }
	  RootTrack->X_q2ex_tr = -final[1] / 1000.;
          RootTrack->Y_q2ex_tr =  final[0] / 1000.;
          RootTrack->Z_q2ex_tr =  final[2] / 1000.;
	  RootTrack->P_q2ex_tr =v3p.mag()/GeV;
	  
	  Transform::single_rot(v3p.x(), v3p.y(), v3p.z(), -( RootTrack->X_q2ex > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_q2ex_tr, RootTrack->Phi_q2ex_tr);
	  //Transform::P_HCS2TCS(RootTrack->Theta_q2ex, RootTrack->Phi_q2ex, ( RootTrack->X_q2ex < 0. ? mLHRSAngle : mRHRSAngle ),
	  //RootTrack->Theta_q2ex_tr, RootTrack->Phi_q2ex_tr);
	  //RootTrack->Theta_q2ex_tr = -tan(RootTrack->Theta_q2ex_tr);
	  //RootTrack->Phi_q2ex_tr   = -tan(RootTrack->Phi_q2ex_tr);
	  
	}else if ( strPhysVolName.find("virtualBoundaryPhys_den") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_den=xx;
	  RootTrack->Y_den=yy;
	  RootTrack->Z_den=zz;
	  RootTrack->Theta_den=v3p.theta()/rad;
	  RootTrack->Phi_den=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_den > 0. ){
		final [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		final [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }
	    }
	  }
	  RootTrack->X_den_tr = -final[1] / 1000.;
          RootTrack->Y_den_tr =  final[0] / 1000.;
          RootTrack->Z_den_tr =  final[2] / 1000.;
	  RootTrack->P_den_tr=v3p.mag()/GeV;
	  Transform::single_rot(v3p.x(), v3p.y(), v3p.z(), -( RootTrack->X_den > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_den_tr, RootTrack->Phi_den_tr);
	  //Transform::P_HCS2TCS(RootTrack->Theta_den, RootTrack->Phi_den, ( RootTrack->X_den < 0. ? mLHRSAngle : mRHRSAngle ),
	  //RootTrack->Theta_den_tr, RootTrack->Phi_den_tr);
	  //RootTrack->Theta_den_tr = -tan(RootTrack->Theta_den_tr);
	  //RootTrack->Phi_den_tr   = -tan(RootTrack->Phi_den_tr);
	  
	}else if ( strPhysVolName.find("virtualBoundaryPhys_dex") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_dex=xx;
	  RootTrack->Y_dex=yy;
	  RootTrack->Z_dex=zz;
	  RootTrack->Theta_dex=v3p.theta()/rad;
	  RootTrack->Phi_dex=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double interim[3]  = { 0.0, 0.0, 0.0 };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  //G4cout << subtract1[0] << " " << subtract1[1] << " " << subtract1[2] << G4endl;
	  //G4cout << initial[0]   << " " << initial  [1] << " " << initial  [2] << G4endl;
	  
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_dex > 0. ){
		interim [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		interim [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }	      
	    }
	    //interim [iii] -= subtract1[iii];
	  }
	  interim[0] -= subtract2[0];
	  interim[1] -= subtract2[1];
	  interim[2] -= subtract2[2];
	  
	  //G4cout << interim[0]   << " " << interim  [1] << " " << interim  [2] << G4endl;
	  for( int iii = 0; iii < 3; iii++ ){//c
	    for( int jjj = 0; jjj < 3; jjj++ ){//a
	      final [iii] += rot_mat_x[iii][jjj] * interim [jjj];
	    }
	  }
	  //G4cout << final  [0]   << " " << final    [1] << " " << final    [2] << G4endl;
	  RootTrack->X_dex_tr = -final[1] / 1000.;
          RootTrack->Y_dex_tr =  final[0] / 1000.;
          RootTrack->Z_dex_tr =  final[2] / 1000.;
	  RootTrack->P_dex_tr =v3p.mag()/GeV;
	  Transform::double_rot(v3p.x(), v3p.y(), v3p.z(), pi / 4., -( RootTrack->X_dex > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_dex_tr, RootTrack->Phi_dex_tr);
	  //Transform::P_HCS2TCS(RootTrack->Theta_dex, RootTrack->Phi_dex, -( RootTrack->X_dex < 0. ? mLHRSAngle : mRHRSAngle ),
	  //RootTrack->Theta_dex_tr, RootTrack->Phi_dex_tr);
	  //RootTrack->Theta_dex_tr = -tan(RootTrack->Theta_dex_tr - 45.0 * deg );
	  //RootTrack->Phi_dex_tr   = -tan(RootTrack->Phi_dex_tr );
	  //G4cout << RootTrack->X_dex_tr << " " << RootTrack->Y_dex_tr << " " << RootTrack->Z_dex_tr << G4endl;

	}else if ( strPhysVolName.find("virtualBoundaryPhys_q3en") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_q3en=xx;
	  RootTrack->Y_q3en=yy;
	  RootTrack->Z_q3en=zz;
	  RootTrack->Theta_q3en=v3p.theta()/rad;
	  RootTrack->Phi_q3en=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double interim[3]  = { 0.0, 0.0, 0.0 };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_q3en > 0.){
		interim [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		interim [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }
	    }
	    //interim [iii] -= subtract1[iii];
	  }
	  interim[0] -= subtract1[0];
	  interim[1] -= subtract1[1];
	  interim[2] -= subtract1[2];
	  
	  for( int iii = 0; iii < 3; iii++ ){//c
	    for( int jjj = 0; jjj < 3; jjj++ ){//a
	      final [iii] += rot_mat_x[iii][jjj] * interim [jjj];
	    }
	  }
	  RootTrack->X_q3en_tr = -final[1] / 1000.;
          RootTrack->Y_q3en_tr =  final[0] / 1000.;
          RootTrack->Z_q3en_tr =  final[2] / 1000.;
	  RootTrack->P_q3en_tr =v3p.mag()/GeV;

	  Transform::double_rot(v3p.x(), v3p.y(), v3p.z(), pi / 4., -( RootTrack->X_q3en > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_q3en_tr, RootTrack->Phi_q3en_tr);
	  //Transform::P_HCS2TCS(RootTrack->Theta_q3en, RootTrack->Phi_q3en, -( RootTrack->X_q3en < 0. ? mLHRSAngle : mRHRSAngle ),
	  //RootTrack->Theta_q3en_tr, RootTrack->Phi_q3en_tr);
	  //RootTrack->Theta_q3en_tr = -tan(RootTrack->Theta_q3en_tr - 45.0 * deg);
	  //RootTrack->Phi_q3en_tr   = -tan(RootTrack->Phi_q3en_tr );
	  //G4cout << RootTrack->X_q3en_tr << " " << RootTrack->Y_q3en_tr << " " << RootTrack->Z_q3en_tr << G4endl;
	}else if ( strPhysVolName.find("virtualBoundaryPhys_q3ex") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_q3ex=xx;
	  RootTrack->Y_q3ex=yy;
	  RootTrack->Z_q3ex=zz;
	  RootTrack->Theta_q3ex=v3p.theta()/rad;
	  RootTrack->Phi_q3ex=v3p.phi()/rad;
	  //RootTrack->VBName=strPhysVolName.substr(foundVB);
	  //kill this track whenever it hit the virtual boundary
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double interim[3]  = { 0.0, 0.0, 0.0 };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_q3ex > 0.){
		interim [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		interim [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }
	    }
	    //interim [iii] -= subtract1[iii];
	  }
	  interim[0] -= subtract1[0];
	  interim[1] -= subtract1[1];
	  interim[2] -= subtract1[2];
	  
	  for( int iii = 0; iii < 3; iii++ ){//c
	    for( int jjj = 0; jjj < 3; jjj++ ){//a
	      final [iii] += rot_mat_x[iii][jjj] * interim [jjj];
	    }
	  }
	  RootTrack->X_q3ex_tr = -final[1] / 1000.;
          RootTrack->Y_q3ex_tr =  final[0] / 1000.;
          RootTrack->Z_q3ex_tr =  final[2] / 1000.;
	  RootTrack->P_q3ex_tr =v3p.mag()/GeV;
	  Transform::double_rot(v3p.x(), v3p.y(), v3p.z(), pi / 4., -( RootTrack->X_q3ex > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_q3ex_tr, RootTrack->Phi_q3ex_tr);
	  //Transform::P_HCS2TCS(RootTrack->Theta_q3ex, RootTrack->Phi_q3ex, -( RootTrack->X_q3ex < 0. ? mLHRSAngle : mRHRSAngle ),
	  //RootTrack->Theta_q3ex_tr, RootTrack->Phi_q3ex_tr);
	  //RootTrack->Theta_q3ex_tr = -tan(RootTrack->Theta_q3ex_tr - 45.0 * deg);
	  //RootTrack->Phi_q3ex_tr   = -tan(RootTrack->Phi_q3ex_tr );

	}else if ( strPhysVolName.find("virtualBoundaryPhys_vdc") != string::npos ){
	  //hit the virtualBoundary, record this step then kill this track
	  //cout << "Recording at the virtual boundary!" << endl;
	  
	  RootTrack->X_vdc=xx;
	  RootTrack->Y_vdc=yy;
	  RootTrack->Z_vdc=zz;
	  RootTrack->P_vdc=v3p.mag()/GeV;
	  RootTrack->Theta_vdc=v3p.theta()/rad;
	  RootTrack->Phi_vdc=v3p.phi()/rad;
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double interim[3]  = { 0.0, 0.0, 0.0 };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_vdc > 0. ){
		interim [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		interim [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }
	    }
	  }
	  interim[0] -= subtract1[0];
	  interim[1] -= subtract1[1];
	  interim[2] -= subtract1[2];
	  
	  for( int iii = 0; iii < 3; iii++ ){//c
	    for( int jjj = 0; jjj < 3; jjj++ ){//a
	      final [iii] += rot_mat_x[iii][jjj] * interim [jjj];
	    }
	  }
	  RootTrack->X_vdc_tr = -final[1] / 1000.;
          RootTrack->Y_vdc_tr =  final[0] / 1000.;
          RootTrack->Z_vdc_tr =  final[2] / 1000.;
	  RootTrack->P_vdc_tr =  v3p.mag()/GeV;
	  Transform::double_rot(v3p.x(), v3p.y(), v3p.z(), pi / 4., -( RootTrack->X_vdc > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_vdc_tr, RootTrack->Phi_vdc_tr);
	}else if ( strPhysVolName.find("virtualBoundaryPhys_qz1") != string::npos ){
	  
	  RootTrack->X_qz1=xx;
	  RootTrack->Y_qz1=yy;
	  RootTrack->Z_qz1=zz;
	  RootTrack->P_qz1=v3p.mag()/GeV;
	  RootTrack->Theta_qz1=v3p.theta()/rad;
	  RootTrack->Phi_qz1=v3p.phi()/rad;
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double interim[3]  = { 0.0, 0.0, 0.0 };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_qz1 > 0. ){
		interim [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		interim [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }
	    }
	  }
	  interim[0] -= subtract1[0];
	  interim[1] -= subtract1[1];
	  interim[2] -= subtract1[2];
	  
	  for( int iii = 0; iii < 3; iii++ ){//c
	    for( int jjj = 0; jjj < 3; jjj++ ){//a
	      final [iii] += rot_mat_x[iii][jjj] * interim [jjj];
	    }
	  }
	  RootTrack->X_qz1_tr = -final[1] / 1000.;
          RootTrack->Y_qz1_tr =  final[0] / 1000.;
          RootTrack->Z_qz1_tr =  final[2] / 1000.;
	  RootTrack->P_qz1_tr =v3p.mag()/GeV;
	  Transform::double_rot(v3p.x(), v3p.y(), v3p.z(), pi / 4., -( RootTrack->X_qz1 > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_qz1_tr, RootTrack->Phi_qz1_tr);
	
	}else if ( strPhysVolName.find("virtualBoundaryPhys_qz2") != string::npos ){
	  
	  RootTrack->X_qz2=xx;
	  RootTrack->Y_qz2=yy;
	  RootTrack->Z_qz2=zz;
	  RootTrack->P_qz2=v3p.mag()/GeV;
	  RootTrack->Theta_qz2=v3p.theta()/rad;
	  RootTrack->Phi_qz2=v3p.phi()/rad;
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double interim[3]  = { 0.0, 0.0, 0.0 };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_qz2 > 0. ){
		interim [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		interim [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }
	    }
	  }
	  interim[0] -= subtract1[0];
	  interim[1] -= subtract1[1];
	  interim[2] -= subtract1[2];
	  
	  for( int iii = 0; iii < 3; iii++ ){//c
	    for( int jjj = 0; jjj < 3; jjj++ ){//a
	      final [iii] += rot_mat_x[iii][jjj] * interim [jjj];
	    }
	  }
	  RootTrack->X_qz2_tr = -final[1] / 1000.;
          RootTrack->Y_qz2_tr =  final[0] / 1000.;
          RootTrack->Z_qz2_tr =  final[2] / 1000.;
	  RootTrack->P_qz2_tr =v3p.mag()/GeV;
	  Transform::double_rot(v3p.x(), v3p.y(), v3p.z(), pi / 4., -( RootTrack->X_qz2 > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_qz2_tr, RootTrack->Phi_qz2_tr);
	
	}else if ( strPhysVolName.find("virtualBoundaryPhys_fp")  != string::npos ){
	  RootTrack->X_fp=xx;
	  RootTrack->Y_fp=yy;
	  RootTrack->Z_fp=zz;
	  RootTrack->X_vb=xx;
	  RootTrack->Y_vb=yy;
	  RootTrack->Z_vb=zz;
	  RootTrack->P_vb=v3p.mag()/GeV;
	  RootTrack->Theta_fp=v3p.theta()/rad;
	  RootTrack->Phi_fp=v3p.phi()/rad;
	  RootTrack->Theta_vb=v3p.theta()/rad;
	  RootTrack->Phi_vb=v3p.phi()/rad;
	  G4double initial[3]  = { xx , yy , zz  };
	  G4double interim[3]  = { 0.0, 0.0, 0.0 };
	  G4double final  [3]  = { 0.0, 0.0, 0.0 };
	  for( int iii = 0; iii < 3; iii++ ){//a
	    for( int jjj = 0; jjj < 3; jjj++ ){//b
	      if( RootTrack->X_fp > 0. ){
		interim [iii] += lrot_mat_y[iii][jjj] * initial[jjj];
	      }else{
		interim [iii] += rrot_mat_y[iii][jjj] * initial[jjj];
	      }
	    }
	  }
	  interim[0] -= subtract1[0];
	  interim[1] -= subtract1[1];
	  interim[2] -= subtract1[2];
	  
	  for( int iii = 0; iii < 3; iii++ ){//c
	    for( int jjj = 0; jjj < 3; jjj++ ){//a
	      final [iii] += rot_mat_x[iii][jjj] * interim [jjj];
	    }
	  }
	  RootTrack->X_fp_tr = -final[1] / 1000.;
          RootTrack->Y_fp_tr =  final[0] / 1000.;
          RootTrack->Z_fp_tr =  final[2] / 1000.;
	  RootTrack->X_vb_tr = -final[1] / 1000.;
          RootTrack->Y_vb_tr =  final[0] / 1000.;
          RootTrack->Z_vb_tr =  final[2] / 1000.;
	  RootTrack->P_fp_tr =v3p.mag()/GeV;
	  //RootTrack->P_vb_tr =v3p.mag()/GeV;

	  //Transform::P_HCS2TCS(RootTrack->Theta_fp,    RootTrack->Phi_fp, -( RootTrack->X_fp < 0. ? mLHRSAngle : mRHRSAngle ),
	  //RootTrack->Theta_fp_tr, RootTrack->Phi_fp_tr);
	  Transform::double_rot(v3p.x(), v3p.y(), v3p.z(), pi / 4., -( RootTrack->X_fp > 0. ? mLHRSAngle : mRHRSAngle ), RootTrack->Theta_fp_tr, RootTrack->Phi_fp_tr);
	  //RootTrack->Theta_fp_tr = -tan(RootTrack->Theta_fp_tr - 45.0 * deg);
	  //RootTrack->Phi_fp_tr   = -tan(RootTrack->Phi_fp_tr );
	  theTrack->SetTrackStatus(fStopAndKill);
	}else if ( strPhysVolName.find("PaulPhys") != string::npos ){
	  theTrack->SetTrackStatus(fStopAndKill);
	}else if ( strPhysVolName.find("LQ1vacPhys") != string::npos ){
	  theTrack->SetTrackStatus(fStopAndKill);
	}else if ( strPhysVolName.find("LQ2vacPhys") != string::npos ){
	  theTrack->SetTrackStatus(fStopAndKill);
	}else if ( strPhysVolName.find("LQ2vac2Phys") != string::npos ){
	  theTrack->SetTrackStatus(fStopAndKill);
	  //}else if ( strPhysVolName.find("LQ3vacPhys") != string::npos ){
	  //theTrack->SetTrackStatus(fStopAndKill);
	  //}else if ( strPhysVolName.find("LQ3vac2Phys") != string::npos ){
	  //theTrack->SetTrackStatus(fStopAndKill);
	  //}else if ( strPhysVolName.find("LQ3Phys") != string::npos ){
	  //theTrack->SetTrackStatus(fStopAndKill);
	}else if ( strPhysVolName.find("RQ1vacPhys") != string::npos ){
	  theTrack->SetTrackStatus(fStopAndKill);
	}else if ( strPhysVolName.find("RQ2vacPhys") != string::npos ){
	  theTrack->SetTrackStatus(fStopAndKill);
	}else if ( strPhysVolName.find("RQ2vac2Phys") != string::npos ){
	  theTrack->SetTrackStatus(fStopAndKill);
	  //}else if ( strPhysVolName.find("RQ3vacPhys") != string::npos ){
	  //theTrack->SetTrackStatus(fStopAndKill);
	  //}else if ( strPhysVolName.find("RQ3vac2Phys") != string::npos ){
	  //theTrack->SetTrackStatus(fStopAndKill);
	  //}else if ( strPhysVolName.find("RQ3Phys") != string::npos ){
	  //theTrack->SetTrackStatus(fStopAndKill);
	}else if ( strPhysVolName.find("absorberPhys") != string::npos ){
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
	RootTrack->StepBx[i] = tmpField[0] / tesla ;
	RootTrack->StepBy[i] = tmpField[1] / tesla ;
	RootTrack->StepBz[i] = tmpField[2] / tesla ;

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

