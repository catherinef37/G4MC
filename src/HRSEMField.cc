// ********************************************************************
//
// $Id: HRSEMField.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//   User Field class Setup implementation.
//
//  
#include "HRSEMField.hh"
#include "HRSEMFieldMessenger.hh"
#include "UsageManager.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

HRSEMField::HRSEMField():mBField_Helm(0), mBField_Septum(0), mBField_SBS(0)
{
	messenger = new HRSEMFieldMessenger(this); 

	UsageManager* gConfig=UsageManager::GetUsageManager();

	//currently target field is always on
	//determine whether helm field is needed
	string pTargetFieldIni=gConfig->GetArgument("TargetFieldIni");
	string pTargetFieldMap=gConfig->GetArgument("TargetFieldMap");
	mBField_Helm = new BField_Helm(pTargetFieldIni.c_str() ,pTargetFieldMap.c_str());

	//determine whether septum field is needed
	//if HRS is built, septum will be loaded
	//There is a potential problem here: if one set SetupLHRS=SetupRHRS=0 in Detector.ini
	//and later on turn them on in G2P or CREX, this class will only response to Detector.ini
	//but not other config files, therefore the septum field will not setup properly
	int pSetupLHRS=0,pSetupRHRS=0;
	gConfig->GetParameter("SetupLHRS",pSetupLHRS); 
	gConfig->GetParameter("SetupRHRS",pSetupRHRS);
	if(pSetupLHRS || pSetupRHRS)
	{
		double pLHRSMomentum,pRHRSMomentum;
		gConfig->GetArgument("LHRSMomentum",pLHRSMomentum); 
		gConfig->GetArgument("RHRSMomentum",pRHRSMomentum);
		string pSeptumFieldIni=gConfig->GetArgument("SeptumFieldIni");
		string pSeptumFieldMap=gConfig->GetArgument("SeptumFieldMap");
		mBField_Septum = new BField_Septum(pLHRSMomentum,pRHRSMomentum,
			pSeptumFieldIni.c_str(),pSeptumFieldMap.c_str());
	}

	int pSetupSuperBigBite=0;
	gConfig->GetParameter("SetupSuperBigBite",pSetupSuperBigBite); 
	if(pSetupSuperBigBite)
	{
		mBField_SBS = new BField_SBS();
	}

	ErDC = 0 *kilovolt/cm; 
	ErInner = 0 *kilovolt/cm; 

	bUseUniformEField=true;
	bUseUniformBField=false;
	EField3V.set(0,0,0);
	BField3V.set(0,0,0);
}

//////////////////////////////////////////////////////////////////////////
//
//  Deconstructors:
HRSEMField::~HRSEMField()
{
	delete messenger;
	if(mBField_Helm)   delete mBField_Helm;
	if(mBField_Septum) delete mBField_Septum;
	if(mBField_SBS)    delete mBField_SBS;
}


////////////////////////////////////////////////////////////////////////////
//input Point[4] (x,y,z,t) 
//
inline void HRSEMField::GetFieldValue(const G4double Point[4],G4double *Bfield) const
{  
	//////////////////////////////////////////////////////////
	//get BField
	if(this->bUseUniformBField) 
	{
		Bfield[0]=BField3V.x();
		Bfield[1]=BField3V.y();
		Bfield[2]=BField3V.z();
	}
	else
	{
		double pB[3],pPos[3]={Point[0]/cm,Point[1]/cm,Point[2]/cm};  //turn into cm
		for(int i=0;i<3;i++) Bfield[i]=0.0;  //reset

		//target field read from map
		if(mBField_Helm)
		{
			for(int i=0;i<3;i++) pB[i]=0.0;  //reset
			if (! mBField_Helm->IsUniformField() )  mBField_Helm->GetBField(pPos,pB); 
			else  mBField_Helm->GetUniformField(pB); 
			for(int i=0;i<3;i++) Bfield[i]+=pB[i]*tesla;
		}

		//septum field read from map
		if(mBField_Septum)
		{
			for(int i=0;i<3;i++) pB[i]=0.0;  //reset
			if (! mBField_Septum->IsUniformField() )  mBField_Septum->GetBField(pPos,pB); 
			else  mBField_Septum->GetUniformField(pB); 
			for(int i=0;i<3;i++) Bfield[i]+=pB[i]*tesla;
		}

		//SBS field read from map
		if(mBField_SBS)
		{
			for(int i=0;i<3;i++) pB[i]=0.0;  //reset
			if (! mBField_SBS->IsUniformField() )  mBField_SBS->GetBField(pPos,pB); 
			else  mBField_SBS->GetUniformField(pB); 
			for(int i=0;i<3;i++) Bfield[i]+=pB[i]*tesla;
		}
	}

	//////////////////////////////////////////////////////////
	//get EFiled,
	//g2p has no Efiled, 

	double  *Efield=&Bfield[3];
	if(this->bUseUniformEField) 
	{
		Efield[0]=EField3V.x();
		Efield[1]=EField3V.y();
		Efield[2]=EField3V.z();
	}
	else
	{
		for(int i=0;i<3;i++) Efield[i]=0.;
		//the following is an example to get Efield, it rovides 2 different fields
		//in 2 regeon
		//G4double pR= sqrt(sqr(Point[0])+sqr(Point[1]));
		////apply region condition 
		//if(fabs(Point[2])<105.0*mm && pR<60.0*mm && pR>30.00607*mm)
		//{           
		//	Efield[0]=-1.*ErDC* Point[0]/pR; 
		//	Efield[1]=-1.*ErDC* Point[1]/pR; 
		//	Efield[2]=0.;
		//}
		//else if(fabs(Point[2])<105.0*mm && pR<30.0*mm && pR>20.00607*mm)
		//{       
		//	Efield[0]=ErInner* Point[0]/pR; 
		//	Efield[1]=ErInner* Point[1]/pR; 
		//	Efield[2]=0.;    
		//}
		//else
		//{ 
		//	for(int i=0;i<3;i++) Efield[i]=0.; 
		//}
	}
}
