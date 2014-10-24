// ********************************************************************
//
// $Id: HRSFZMagField.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//   User Field class Setup implementation.
//
//  
#include "HRSFZMagField.hh"
//#include "HRSFZMagFieldMessenger.hh"
#include "UsageManager.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

static const double inch = 2.54*cm;
HRSFZMagField::HRSFZMagField()
{
	//messenger = new HRSFZMagFieldMessenger(this); 
	
	UsageManager* gConfig=UsageManager::GetUsageManager();

	bIsFZB1Uniform=true;
	bIsFZB2Uniform=true;

	int pUseDefaultFZB1=1,pUseDefaultFZB2=1;
	gConfig->GetArgument("UseDefaultFZB1",pUseDefaultFZB1);
	gConfig->GetArgument("UseDefaultFZB2",pUseDefaultFZB2); 
	
	gConfig->GetParameter("PivotXOffset",mPivotXOffset); mPivotXOffset*=mm;
	gConfig->GetParameter("PivotYOffset",mPivotYOffset); mPivotYOffset*=mm; 
	gConfig->GetParameter("PivotZOffset",mPivotZOffset); mPivotZOffset*=mm;

	double mFZB1TiltedAngle=0,mFZB2TiltedAngle=0;

	double mFZB1PosX=0.0*cm+mPivotXOffset;
	double mFZB1PosY=0.0*cm+mPivotYOffset;
	double mFZB1PosZ=-622*cm+mPivotZOffset;

	double mFZB2PosX=0.0*cm+mPivotXOffset;
	double mFZB2PosY=0.0*cm+mPivotYOffset;
	double mFZB2PosZ=-266*cm+mPivotZOffset;

	double mFZB1Bx=0,mFZB1By=0,mFZB1Bz=0;
	double mFZB2Bx=0,mFZB2By=0,mFZB2Bz=0;
	///////////////////////////////////////////////////////
	//Here is all the default settings. will be used based on 
	//the beam energy, HRS angle and target field setting

	//judge the whcih setting it is
	
	double pLHRSAngle=5.69*deg,pLSeptumAngle=12.5*deg;
	gConfig->GetParameter("LHRSAngle",pLHRSAngle);  pLHRSAngle*=deg;
	gConfig->GetParameter("LSeptumAngle",pLSeptumAngle);  pLSeptumAngle*=deg;
	bool bIsSeptaIn=(fabs(pLHRSAngle-pLSeptumAngle)>0.1*deg)?true:false;

	double pBeamE_GeV=2.257;
	gConfig->GetArgument("BeamEnergy",pBeamE_GeV);   //in GeV

	double pHelm_CurrentRatio=1.0;
	gConfig->GetParameter("Helm_CurrentRatio",pHelm_CurrentRatio);

	int pIsG2P1OrGEP2=0;
	double pHeam_Rot1Angle=270*deg;
	gConfig->GetParameter("Helm_RotAngle1",pHeam_Rot1Angle);  pHeam_Rot1Angle*=deg;
	if(fabs(pHeam_Rot1Angle-270*deg)<5) pIsG2P1OrGEP2=1;
	else if(fabs(pHeam_Rot1Angle-354*deg)<5) pIsG2P1OrGEP2=2;

	/*B field in unit of KG
	                         B1     B2
1.159 GeV/c:
  2.5T  90deg  HAdump     -3.02  7.02
  2.5T  20deg  Local      -1.02  2.38

 1.706 GeV/c
    20 deg  2.5T  HAdump  -1.02   2.38
    90 deg  2.5T  HAdump  -2.99   6.98

2.257 GeV/c:
  2.5T  20deg  HAdump    -1.02   2.38
  2.5T  90deg  Local     -1.47   3.45
  5.01T 20deg  Local     -2.04   4.77
  5.01T 90deg  Local     -3.44   8.1

3.355 GeV/c

  5.01T 90deg  Local     -2.94   6.92

WITHOUT SEPTA (TARGET AT PIVOT)

1.159 GeV/c
  2.5T  90deg  Local     -1.97   3.95
2.257 GeV/c
  5.01T 90deg  Local     -3.94   7.92
3.355 GeV/c
  5.01T 90deg  Local     -3.94   7.92
  */
	////////////////////////////////////////////////////////
	if(pUseDefaultFZB1)
	{
		//HRS=12.5 degree, TargetFiled=2.5T, Beam=1.159GeV
		if(!bIsSeptaIn && fabs(pHelm_CurrentRatio)<0.7 && fabs(pBeamE_GeV-1.159)<0.3) 
		{
			if(pIsG2P1OrGEP2==1)
			{
				mFZB1TiltedAngle = -3.02*deg;
				mFZB1PosY       = -4*cm+mPivotYOffset;
				mFZB1Bx          = 0.197*tesla;
			}
			else if(pIsG2P1OrGEP2==2)
			{
				mFZB1TiltedAngle = -3.02*deg;
				mFZB1PosY       = -4*cm+mPivotYOffset;
				mFZB1Bx          = 0.197*tesla;
			}

		}
		else if(!bIsSeptaIn && fabs(pHelm_CurrentRatio)>0.9 && fabs(pBeamE_GeV-2.257)<0.3) 
		{	
			if(pIsG2P1OrGEP2==1)
			{
				mFZB1TiltedAngle = -1.55*deg;
				mFZB1PosY       = -4*cm+mPivotYOffset;
				mFZB1Bx          = 0.394*tesla;
			}
			else if(pIsG2P1OrGEP2==2)
			{	
				mFZB1TiltedAngle = -1.55*deg;
				mFZB1PosY       = -4*cm+mPivotYOffset;
				mFZB1Bx          = 0.394*tesla;
			}

		}
		else if(!bIsSeptaIn && fabs(pHelm_CurrentRatio)>0.9 && fabs(pBeamE_GeV-3.359)<0.3)
		{	
			if(pIsG2P1OrGEP2==1)
			{
				mFZB1TiltedAngle = -1.04*deg;
				mFZB1PosY       = -4*cm+mPivotYOffset;
				mFZB1Bx          = 0.394*tesla;
			}
			else if(pIsG2P1OrGEP2==2)
			{
				mFZB1TiltedAngle = -1.04*deg;
				mFZB1PosY       = -4*cm+mPivotYOffset;
				mFZB1Bx          = 0.394*tesla;
			}
		}
		else if(true)
		{
			mFZB1TiltedAngle = -3.02*deg;
			mFZB1PosY       = -4*cm+mPivotYOffset;
			mFZB1Bx          = 0.197*tesla;
		}

		//I have to updated these parameters in the map so they can be written into the config tree
		gConfig->SetArgument("FZB1TiltedAngle",mFZB1TiltedAngle/rad); 
		gConfig->SetArgument("FZB1PosX",mFZB1PosX/mm); 
		gConfig->SetArgument("FZB1PosY",mFZB1PosY/mm); 
		gConfig->SetArgument("FZB1PosZ",mFZB1PosZ/mm); 
		gConfig->SetArgument("FZB1Bx",mFZB1Bx/tesla); 
		gConfig->SetArgument("FZB1By",mFZB1By/tesla); 
		gConfig->SetArgument("FZB1Bz",mFZB1Bz/tesla); 
	}
	else
	{
		//reading the value from input arguments of option -FZB1 or -ChicaneMagnet1
		gConfig->GetArgument("FZB1TiltedAngle",mFZB1TiltedAngle); mFZB1TiltedAngle*=deg;
		gConfig->GetArgument("FZB1PosX",mFZB1PosX); mFZB1PosX*=mm;
		gConfig->GetArgument("FZB1PosY",mFZB1PosY); mFZB1PosY*=mm;
		gConfig->GetArgument("FZB1PosZ",mFZB1PosZ); mFZB1PosZ*=mm;
		gConfig->GetArgument("FZB1Bx",mFZB1Bx); mFZB1Bx*=tesla;
		gConfig->GetArgument("FZB1By",mFZB1By); mFZB1By*=tesla;
		gConfig->GetArgument("FZB1Bz",mFZB1Bz); mFZB1Bz*=tesla;
	}

	////////////////////////////////////////////////////////
	if(pUseDefaultFZB2)
	{
		//HRS=12.5 degree, TargetFiled=2.5T, Beam=1.159GeV
		if(!bIsSeptaIn && fabs(pHelm_CurrentRatio)>0.9 && fabs(pBeamE_GeV-1.159)<0.3)
		{
			if(pIsG2P1OrGEP2==1)
			{
				mFZB2TiltedAngle = 7.06*deg;
				mFZB2PosY       = -24*cm+mPivotYOffset;
				mFZB2Bx          = -0.395*tesla;
			}
			else if(pIsG2P1OrGEP2==2)
			{
				mFZB2TiltedAngle = 7.06*deg;
				mFZB2PosY       = -24*cm+mPivotYOffset;
				mFZB2Bx          = -0.395*tesla;
			}
		}
		else if(!bIsSeptaIn && fabs(pHelm_CurrentRatio)>0.9 && fabs(pBeamE_GeV-2.257)<0.3) 
		{
			if(pIsG2P1OrGEP2==1)
			{
				mFZB2TiltedAngle = 3.62*deg;
				mFZB2PosY       = -24*cm+mPivotYOffset;
				mFZB2Bx          = -0.792*tesla;
			}
			else if(pIsG2P1OrGEP2==2)
			{
				mFZB2TiltedAngle = 3.62*deg;
				mFZB2PosY       = -24*cm+mPivotYOffset;
				mFZB2Bx          = -0.792*tesla;
			}
		}
		else if(!bIsSeptaIn && fabs(pHelm_CurrentRatio)>0.9 && fabs(pBeamE_GeV-3.359)<0.3)
		{
			if(pIsG2P1OrGEP2==1)
			{
				mFZB2TiltedAngle = 2.44*deg;
				mFZB2PosY       = -12*cm+mPivotYOffset;
				mFZB2Bx          = -0.792*tesla;
			}
			else if(pIsG2P1OrGEP2==2)
			{
				mFZB2TiltedAngle = 2.44*deg;
				mFZB2PosY       = -12*cm+mPivotYOffset;
				mFZB2Bx          = -0.792*tesla;
			}
		}
		else if(true)
		{
			mFZB2TiltedAngle = 7.06*deg;
			mFZB2PosY       = -24*cm+mPivotYOffset;
			mFZB2Bx          = -0.395*tesla;
		}

		//I have to updated these parameters in the map so they can be written into the config tree
		gConfig->SetArgument("FZB2TiltedAngle",mFZB2TiltedAngle/rad); 
		gConfig->SetArgument("FZB2PosX",mFZB2PosX/mm); 
		gConfig->SetArgument("FZB2PosY",mFZB2PosY/mm); 
		gConfig->SetArgument("FZB2PosZ",mFZB2PosZ/mm); 
		gConfig->SetArgument("FZB2Bx",mFZB2Bx/tesla); 
		gConfig->SetArgument("FZB2By",mFZB2By/tesla); 
		gConfig->SetArgument("FZB2Bz",mFZB2Bz/tesla); 
	}
	else
	{
		//reading the value from input arguments of option -FZB2 or -ChicaneMagnet2
		gConfig->GetArgument("FZB2TiltedAngle",mFZB2TiltedAngle); mFZB2TiltedAngle*=deg;
		gConfig->GetArgument("FZB2PosX",mFZB2PosX); mFZB2PosX*=mm;
		gConfig->GetArgument("FZB2PosY",mFZB2PosY); mFZB2PosY*=mm;
		gConfig->GetArgument("FZB2PosZ",mFZB2PosZ); mFZB2PosZ*=mm;
		gConfig->GetArgument("FZB2Bx",mFZB2Bx); mFZB2Bx*=tesla;
		gConfig->GetArgument("FZB2By",mFZB2By); mFZB2By*=tesla;
		gConfig->GetArgument("FZB2Bz",mFZB2Bz); mFZB2Bz*=tesla;
	}
	////////////////////////////////////////////////////////

	//now build the vectors
	 mFZB1Field3V.set(mFZB1Bx,mFZB1By,mFZB1Bz);
	 mFZB2Field3V.set(mFZB2Bx,mFZB2By,mFZB2Bz);

	//some const number to descripe the size of field container
	double mFieldHalfX,mFieldHalfY,mFieldHalfZ;
	mFieldHalfX=1.66*inch/2;
	mFieldHalfY=11.5*inch/2;
	mFieldHalfZ=76.34*inch/2;
	
	//the half after rotation
	mFZB1HalfY=fabs(mFieldHalfY/cos(mFZB1TiltedAngle));
	mFZB2HalfY=fabs(mFieldHalfY/cos(mFZB2TiltedAngle));


	//the following is the max after rotation
	double pR=sqrt(mFieldHalfY*mFieldHalfY+mFieldHalfZ*mFieldHalfZ);
	double pAlpha= atan(mFieldHalfY/mFieldHalfZ);
	mFZB1MaxX=mFZB2MaxX=mFieldHalfX;
	mFZB1MaxY=pR*sin(fabs(mFZB1TiltedAngle)+pAlpha);
	mFZB1MaxZ=pR*cos(fabs(mFZB1TiltedAngle)-pAlpha);
	mFZB2MaxY=pR*sin(fabs(mFZB2TiltedAngle)+pAlpha);
	mFZB2MaxZ=pR*cos(fabs(mFZB2TiltedAngle)-pAlpha);

}

//////////////////////////////////////////////////////////////////////////
//
//  Deconstructors:
HRSFZMagField::~HRSFZMagField()
{
	//delete messenger;
}


////////////////////////////////////////////////////////////////////////////
//input pos[4] (x,y,z,t) 
//
inline void HRSFZMagField::GetFieldValue(const G4double pos[4],G4double *Bfield) const
{  
	//reset 
	Bfield[0]=Bfield[1]=Bfield[2]=0;
	
	double x=pos[1]-mPivotXOffset, y=pos[2]-mPivotYOffset, z=pos[3]-mPivotZOffset;

	//judge its location, to save speed, I just 
	bool bIsInFZB1=false;
	bool bIsInFZB2=false;

	//check z to tell which magnet it is 
	if(fabs(z-mFZB1PosZ) < mFZB1MaxZ) bIsInFZB1=true; 
	else if(fabs(z-mFZB2PosZ) < mFZB2MaxZ) bIsInFZB2=true; 
	if(!bIsInFZB1 && !bIsInFZB2)  return; 

	//check x and y
	//after rotation, y=tan(mFZB1TiltedAngle)*(z-mFZB1PosZ)+mFZB1PosY; 
	//x is not changed,  |x|<1.66*inch, the length in z is 76.34 inch
	if (bIsInFZB1)
	{
		if (fabs(x-mFZB1PosX)>mFZB1MaxX)  return; 
		double pMidYAtThisZ=tan(mFZB1TiltedAngle)*(z-mFZB1PosZ)+mFZB1PosY;
		if(fabs(y-pMidYAtThisZ)>mFZB1HalfY) return;
	}
	else if (bIsInFZB2)
	{
		if (fabs(x-mFZB2PosX)>mFZB2MaxX)  return; 
		double pMidYAtThisZ=tan(mFZB2TiltedAngle)*(z-mFZB2PosZ)+mFZB2PosY;
		if(fabs(y-pMidYAtThisZ)>mFZB2HalfY) return;
	}

	//////////////////////////////////////////////////////////
	//now ready to set the field 
	if (bIsInFZB1)
	{
		if(bIsFZB1Uniform) 
		{
			Bfield[0]=mFZB1Field3V.x();
			Bfield[1]=mFZB1Field3V.y();
			Bfield[2]=mFZB1Field3V.z();
		}
	}
	else if (bIsInFZB2)
	{
		if(bIsFZB2Uniform) 
		{
			Bfield[0]=mFZB2Field3V.x();
			Bfield[1]=mFZB2Field3V.y();
			Bfield[2]=mFZB2Field3V.z();
		}
	}
	return;
}
