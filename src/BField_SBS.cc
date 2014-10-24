// ********************************************************************
// $Id: BField_Tosca_SBS.cc,v 3.0, 2011/1/19  G2P Exp $
// Implementation of the BField_Tosca_SBS class.
//
//////////////////////////////////////////////////////////////////////

#include "BField_SBS.hh"
#include "UsageManager.hh"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include  <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"


//#define CREATE_MAP_NTUPLE 1
//#define BFIELD_TOSCA_SBS_DEBUG 0

#ifdef BFIELD_TOSCA_SBS_DEBUG
#include "GlobalDebuger.hh"
#endif


using namespace std;

BField_Tosca_SBS* BField_Tosca_SBS::fInstance=0;
BField_Tosca_SBS* BField_Tosca_SBS::GetInstance()
{ 
	if(!fInstance)  new BField_Tosca_SBS();
	return fInstance; 
}



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BField_Tosca_SBS::BField_Tosca_SBS(double pMomentum,const char *inifile,const char *mapfile)
{
#ifdef BFIELD_TOSCA_SBS_DEBUG
	if(BFIELD_TOSCA_SBS_DEBUG > Global_Debug_Level)
		SetGlobalDebugLevel("BField_Tosca_SBS::BField_Tosca_SBS()", (int)BFIELD_TOSCA_SBS_DEBUG);
#endif

	fInstance=this;

	this->ReadIni(inifile);

	mDoShift=(mOrigin[0]*mOrigin[0]+mOrigin[1]*mOrigin[1]+mOrigin[2]*mOrigin[2]>0)?true:false;
	//mDoRotation=(fabs(mRotAngle[0])+fabs(mRotAngle[1])+fabs(mRotAngle[2])>0)?true:false;

	mDoRotation=false;
	// Default constructor. Gives a unit matrix
	CLHEP::HepRotation pRot[3];

	for(int i=0;i<3;i++)
	{
		if(mRotAxis[i]==0 || fabs(mRotAngle[i])<1.0E-10) continue;
		mDoRotation=true;
		if(mRotAxis[i]==1) pRot[i].rotateX(mRotAngle[i]);
		else if(mRotAxis[i]==2) pRot[i].rotateY(mRotAngle[i]);
		else if(mRotAxis[i]==3) pRot[i].rotateZ(mRotAngle[i]);
	}

	mRotL2F=new CLHEP::HepRotation();
	mRotF2L=new CLHEP::HepRotation();
	//By Jixie; I have to debug this expression, what order shall it be, 210 or 012?
	*mRotL2F=pRot[2]*pRot[1]*pRot[0];
	*mRotF2L=mRotL2F->inverse(); 
#ifdef BFIELD_TOSCA_SBS_DEBUG
	double rad2deg=180./acos(-1.0);
	if(Global_Debug_Level>=1)
	{
		cout<<"\nCLHEP Lab2Field EulerAngles: phi="<<mRotL2F->getPhi()*rad2deg
			<<"  theta="<<mRotL2F->getTheta()*rad2deg
			<<"  psi="<<mRotL2F->getPsi()*rad2deg<<endl;
	}
#endif

	mXNum=int((mXmax-mXmin)/mStepX)+1;
	if (mXNum<2)
	{
		mXNum=3; //set the minimum to 3
		mXmax=mXmin+3*mStepX;
	}
	mYNum=int((mYmax-mYmin)/mStepY)+1;
	if (mYNum<2)
	{
		mYNum=3; //set the minimum to 3
		mYmax=mYmin+3*mStepY;
	}
	mZNum=int((mZmax-mZmin)/mStepZ)+1;
	if (mZNum<2)
	{
		mZNum=3; //set the minimum to 3
		mZmax=mZmin+3*mStepZ;
	}

	///////////allocate start//////////
	int i,j,k,l;
	mBField=new double ***[mXNum];
	for(i=0;i<mXNum;i++)
	{
		mBField[i]=new double **[mYNum];
		for (j=0;j<mYNum;j++)
		{
			mBField[i][j]=new double *[mZNum];
			for	(k=0;k<mZNum;k++)
			{
				mBField[i][j][k]=new double [mNPara];
				//initial the array
				for (l=0;l<mNPara;l++)    mBField[i][j][k][l]=0.0;
			}
		}
	}
	/////////////allocate end////////////

	if(mUseUniformB!=1) this->ReadMap(mapfile);

	//update the Current Ratio if necessary
	if(pMomentum>0.0 && fabs(mCurrentRatio-pMomentum/mDefaultMomentum)>1.0E-04)
	{
		mCurrentRatio=pMomentum/mDefaultMomentum;
		cout<<"\n##BField_Tosca_SBS::BField_Tosca_SBS() Septum field current ratio is set to "<<mCurrentRatio<<endl;
		
		UsageManager *pConfig=UsageManager::GetUsageManager();
		char tmpStr[20];
		sprintf(tmpStr,"%.6f",mCurrentRatio);
		pConfig->SetParameter("Tosca_CurrentRatio",tmpStr);
	}
}

BField_Tosca_SBS::~BField_Tosca_SBS()
{
	int i,j,k;
	for(i=0;i<mXNum;i++)
	{
		for (j=0;j<mYNum;j++)
		{
			for (k=0;k<mZNum;k++)
			{
				delete [] mBField[i][j][k];
			}
			delete [] mBField[i][j];
		}
		delete [] mBField[i];
	}
	delete mRotL2F;
	delete mRotF2L;
}

/////////////////////////////////////////////////////////////////////
bool BField_Tosca_SBS::ReadIni(const char *filename)
{
	double deg2rad=acos(-1.0)/180.;
	//By Jixie: I am not use this routine to read ini file any longer
	//I prefer to use UsageManager::ReadFile
	//
	UsageManager *pConfig=UsageManager::GetUsageManager();
	bool ret=pConfig->ReadFile(filename);

	pConfig->GetParameter("Tosca_SBS_UseUniformB",	mUseUniformB);
	pConfig->GetParameter("Tosca_SBS_UniformB_x",	mUniformB[0]);
	pConfig->GetParameter("Tosca_SBS_UniformB_y",	mUniformB[1]);
	pConfig->GetParameter("Tosca_SBS_UniformB_z",	mUniformB[2]);

	pConfig->GetParameter("Tosca_SBS_FieldUnit",	mFieldUnit);
	pConfig->GetParameter("Tosca_SBS_PosUnit",		mPosUnit);
	pConfig->GetParameter("Tosca_SBS_FirstDataLine",mFirstDataLine);
	pConfig->GetParameter("Tosca_SBS_NPara",		mNPara);
	pConfig->GetParameter("Tosca_SBS_FlipX",		mFlip[0]);
	pConfig->GetParameter("Tosca_SBS_FlipY",		mFlip[1]);
	pConfig->GetParameter("Tosca_SBS_FlipZ",		mFlip[2]);
	pConfig->GetParameter("Tosca_SBS_Direction",	mDirection);

	pConfig->GetParameter("Tosca_SBS_Xmin",			mXmin);
	pConfig->GetParameter("Tosca_SBS_Xmax",			mXmax);
	pConfig->GetParameter("Tosca_SBS_Ymin",			mYmin);
	pConfig->GetParameter("Tosca_SBS_Ymax",			mYmax);
	pConfig->GetParameter("Tosca_SBS_Zmin",			mZmin);
	pConfig->GetParameter("Tosca_SBS_Zmax",			mZmax);
	pConfig->GetParameter("Tosca_SBS_StepX",		mStepX);
	pConfig->GetParameter("Tosca_SBS_StepY",		mStepY);
	pConfig->GetParameter("Tosca_SBS_StepZ",		mStepZ);
	pConfig->GetParameter("Tosca_SBS_InterpolateOutOfRange", mInterpolateOutOfRange);

	pConfig->GetParameter("Tosca_SBS_OriginX",		mOrigin[0]);
	pConfig->GetParameter("Tosca_SBS_OriginY",		mOrigin[1]);
	pConfig->GetParameter("Tosca_SBS_OriginZ",		mOrigin[2]);
	pConfig->GetParameter("Tosca_SBS_RotAxis1",		mRotAxis[0]);
	pConfig->GetParameter("Tosca_SBS_RotAxis2",		mRotAxis[1]);
	pConfig->GetParameter("Tosca_SBS_RotAxis3",		mRotAxis[2]);
	pConfig->GetParameter("Tosca_SBS_RotAngle1",	mRotAngle[0]); mRotAngle[0]*=deg2rad;
	pConfig->GetParameter("Tosca_SBS_RotAngle2",	mRotAngle[1]); mRotAngle[1]*=deg2rad;
	pConfig->GetParameter("Tosca_SBS_RotAngle3",	mRotAngle[2]); mRotAngle[2]*=deg2rad;

	pConfig->GetParameter("Tosca_SBS_DefaultMomentum",mDefaultMomentum);
	pConfig->GetParameter("Tosca_SBS_CurrentRatio",	mCurrentRatio);	

#ifdef BFIELD_TOSCA_SBS_DEBUG
	pConfig->PrintParamMap(); 
#endif
	return ret;

}

/////////////////////////////////////////////////////////////////////
bool BField_Tosca_SBS::ReadMap(const char *filename)
{
	char strLog[1024];
	sprintf(strLog,"BField_Tosca_SBS::ReadMap() is loading field map %s......\n",filename);
	UsageManager::WriteLog(strLog);

	ifstream ins;
	int indexX=0,indexY,indexZ=0,col=0;
	double tempLine[10];
	ins.open(filename);
	if (ins.fail())
	{
		sprintf(strLog,"***ERROR! Can not open field map %s...exit!***\n",filename);
		UsageManager::WriteLog(strLog);		
		exit(-1);
		return false;
	}

	//eat the first mFirstDataLine-1 lines, which is the header block
	int nLine2Eat=mFirstDataLine-1;   //7 lines only for merged.table
	char tempname[256];
	for(int i=0;i<nLine2Eat;i++) ins.getline (tempname,256);

	while (!ins.eof())
	{
		ins.getline(tempname,256);				
		//check if it is an empty line	
		if(strlen(tempname) < size_t(2*mNPara-1)) continue;

		istringstream s1(tempname);
		for(col=0;col<mNPara;col++)    s1>>tempLine[col];

		//check for X, Y and Z
		if (tempLine[0]>=mXmin && tempLine[0]<=mXmax && tempLine[1]>=mYmin && tempLine[1]<=mYmax &&
			tempLine[2]>=mZmin && tempLine[2]<=mZmax)
		{//store the value

			//in case there is an empty line, z=r=Btot=0.0
			if( sqrt(tempLine[0]*tempLine[0]+tempLine[1]*tempLine[1]+tempLine[2]*tempLine[2]) < 1.0E-8 && 
				sqrt(tempLine[3]*tempLine[3]+tempLine[4]*tempLine[4]+tempLine[5]*tempLine[5]) < 1.0E-8 ) 
			{
				cout<<"***Warning: ZERO field at  x="<<tempLine[0]<<"  y="<<tempLine[1]
				<<"  z="<<tempLine[2]<<endl;
				cout<<"There could be an empty line in the map "<<filename<<endl;
			}
			else 
			{
				indexX=int((tempLine[0]-mXmin)/mStepX);
				indexY=int((tempLine[1]-mYmin)/mStepY);
				indexZ=int((tempLine[2]-mZmin)/mStepZ);
				for(col=0;col<mNPara;col++)
				{
					mBField[indexX][indexY][indexZ][col]=tempLine[col]*mPosUnit;		   //change unit to cm
					if(col>=3 && col<=5) mBField[indexX][indexY][indexZ][col]*=mFieldUnit; //change unit to tesla
				}
			}
		}
	}
	ins.close();

#ifdef CREATE_MAP_NTUPLE
	char rootfile[255];
	double x=0.0,y=0.0,z=0.0,Bx=0.0,By=0.0,Bz=0.0,r=0.0,Br=0.0,Btot=0.0;
	sprintf(rootfile,"%s.root",filename);
	TFile *file=new TFile(rootfile,"RECREATE");
	TTree *field=new TTree("field","field map");
	field->Branch("x",&x,"x/D");
	field->Branch("y",&y,"y/D");
	field->Branch("r",&r,"r/D");
	field->Branch("z",&z,"z/D");
	field->Branch("Bx",&Bx,"Bx/D");
	field->Branch("By",&By,"By/D");
	field->Branch("Br",&Br,"Br/D");
	field->Branch("Bz",&Bz,"Bz/D");
	field->Branch("Btot",&Btot,"Btot/D");


#ifdef BFIELD_TOSCA_SBS_DEBUG
	if(Global_Debug_Level>=4)
	{
		printf("The Magnetic field Map is:\n");		
		printf("       x        y        z       Bx       By       Bz        r       Br     Btot\n");
	}
#endif
	for (indexX=0;indexX<mXNum;indexX++)
	{
		for (indexY=0;indexY<mYNum;indexY++)
		{
			for (indexZ=0;indexZ<mZNum;indexZ++)
			{
				x=mBField[indexX][indexY][indexZ][0];
				y=mBField[indexX][indexY][indexZ][1];
				r=sqrt(x*x+y*y);
				z=mBField[indexX][indexY][indexZ][2];
				Bx=mBField[indexX][indexY][indexZ][3];
				By=mBField[indexX][indexY][indexZ][4];  
				Bz=mBField[indexX][indexY][indexZ][5];				
				Br=sqrt(Bx*Bx+By*By);
				Btot=sqrt(Bx*Bx+By*By+Bz*Bz);

				if(fabs(r)<1.0E-8 && fabs(Btot)<1.0E-8) 
				{
					cout<<"***Warning: ZERO field at  x="<<x<<"  y="<<y<<"  z="<<z<<endl;
					cout<<"***Warning: There is an empty line in the map buffer..."<<endl;
				}
				else field->Fill();

#ifdef BFIELD_TOSCA_SBS_DEBUG
				if(Global_Debug_Level>=4) 
				{					
					printf("%8.3f %8.3f %8.3f %8.6f %8.6f %8.6f %8.3f %8.6f %8.6f\n",
						x,y,z,Bx,By,Bz,r,Br,Btot);

				}
#endif
				//reset 
				x=y=z=Bx=By=Bz=r=Br=Btot=0.0;
			}
		}
	}
	file->Write();
	file->Close();
	cout << "Close root file " << rootfile <<endl;  
	file->Delete();
#endif

#ifdef BFIELD_TOSCA_SBS_DEBUG
	if(Global_Debug_Level>=4)
	{
		printf("The Magnetic field Map is:\n");
		printf("      x         y        z        Bx        By        Bz\n");
		for (indexX=0;indexX<mXNum;indexX++)
		{
			for (indexY=0;indexY<mYNum;indexY++)
			{
				for (indexZ=0;indexZ<mZNum;indexZ++)
				{					
					for(col=0;col<3;col++)  printf(" %8.2f ",mBField[indexX][indexY][indexZ][col]);
					for(col=3;col<mNPara;col++)  printf(" %8.6f ",mBField[indexX][indexY][indexZ][col]);
					printf("\n");
				}
			}
		}
	}
#endif

	return true;
}

/////////////////////////////////////////////////////////////////////
bool BField_Tosca_SBS::Interpolation(double Pos[3],double B[3],int n)
{/*////////////////////////////////////////
 //function: calculate the nth order interpolation
 //1)found out B[X0-Xn][Y0-Yn][Z0-Zn], (n+1)x(n+1)x(n+1) matrix, 
 //2)then interpolate by X, found 3x(n+1)x(n+1) matrix
 //3)then interpolate by Y, found 3x(n+1) matrix 
 //4)Finally interpolate by Z to get Bx,By,Bz 
 //Input:
 // Pos[3]: the position in x,y,z coordinate in cm,origin variable
 // n : calculate the Nth order
 //Output:
 // B[3]: the magnetic field in B_x,B_y,B_z in tesla
 *//////////////////////////////////////////

#ifdef BFIELD_TOSCA_SBS_DEBUG
	//I do not want to check these in release version, just make sure you provide the right parameters
	if(n<1) n=1; 
	if(n>4) n=4; 

	if (n>=mXNum || n>=mYNum|| n>=mZNum)
	{
		printf("\n**Error! Too few points to finish %dth order interpolation!**",n);
		printf("**Please change the config file to load more data points or use lower order Number!**\n");
		return false;
	}
#endif

	//if not provide the whole space then need to do flips
	double x=Pos[0],y=Pos[1],z=Pos[2];
	if(mFlip[0]!=0) x=fabs(Pos[0]);
	if(mFlip[1]!=0) y=fabs(Pos[1]);
	if(mFlip[2]!=0) z=fabs(Pos[2]);

	int StartIndexX=0,StartIndexY=0,StartIndexZ=0;  //the first point to do the interpolation
	StartIndexX=int((x-mBField[0][0][0][0])/mStepX);
	StartIndexY=int((y-mBField[0][0][0][1])/mStepY);
	StartIndexZ=int((z-mBField[0][0][0][2])/mStepZ);

#ifdef BFIELD_TOSCA_SBS_DEBUG
	if (StartIndexX<0 || StartIndexX>=mXNum)
	{
		if(BFIELD_TOSCA_SBS_DEBUG>=1)
		printf("\n**Warning!!! X=%.4f is out of range [%.4f,%.4f) !!!**\n",
			x,mBField[0][0][0][0],mBField[mXNum-1][0][0][0]);
		B[0]=B[1]=B[2]=0.0;return false;
	}
	if (StartIndexY<0 || StartIndexY>=mYNum)
	{
		if(BFIELD_TOSCA_SBS_DEBUG>=1)
		printf("\n**Warning!!! Y=%.4f is out of range [%.4f,%.4f) !!!**\n",
			y,mBField[0][0][0][1],mBField[0][mYNum-1][0][1]);
		B[0]=B[1]=B[2]=0.0;return false;
	}
	if (StartIndexZ<0 || StartIndexZ>=mZNum)
	{
		if(BFIELD_TOSCA_SBS_DEBUG>=1)
		printf("\n**Warning!!! Z=%.4f is out of range [%.4f,%.4f) !!!**\n",
			z,mBField[0][0][0][2],mBField[0][0][mZNum-1][2]);
		B[0]=B[1]=B[2]=0.0;return false;
	}
#endif

	if(mInterpolateOutOfRange==0)
	{
		if ((StartIndexX<0 || StartIndexX>=mXNum) || 
			(StartIndexY<0 || StartIndexY>=mYNum) || 
			(StartIndexZ<0 || StartIndexZ>=mZNum) )
		{
			B[0]=B[1]=B[2]=0.0;return false;
		}
	}

	//choose a best position to interpolate, considering the lower and upper limits
	if (StartIndexX>mXNum-1-n) StartIndexX=mXNum-1-n;
	else if (StartIndexX<0) StartIndexX=0;
	if (StartIndexY>mYNum-1-n) StartIndexY=mYNum-1-n;
	else if (StartIndexY<0) StartIndexY=0;
	if (StartIndexZ>mZNum-1-n) StartIndexZ=mZNum-1-n;
	else if (StartIndexZ<0) StartIndexZ=0;

	// interpolate by X first
	//put the n+1=5 to avoid new and delete, this can save a lot of time
	double temp,tempBX[3][5][5],tempBY[3][5];

	int iX,iY,iZ,ii,jj,kk,ll;
	for (ll=0;ll<=n;ll++)
	{
		iZ=StartIndexZ+ll;
		for (ii=0;ii<=n;ii++)
		{
			//initial
			tempBX[0][ii][ll]=0.;
			tempBX[1][ii][ll]=0.;
			tempBX[2][ii][ll]=0.;
			iY=StartIndexY+ii;		
			for (kk=0;kk<=n;kk++)
			{
				temp=1.;
				for (jj=0;jj<=n;jj++)
				{
					iX=StartIndexX+jj;
					if(jj!=kk)
					{
						temp*=(x-mBField[iX][iY][iZ][0])/(mBField[StartIndexX+kk][iY][iZ][0]-mBField[iX][iY][iZ][0]);
					}
				}
				tempBX[0][ii][ll]+=temp*mBField[StartIndexX+kk][iY][iZ][3];
				tempBX[1][ii][ll]+=temp*mBField[StartIndexX+kk][iY][iZ][4];
				tempBX[2][ii][ll]+=temp*mBField[StartIndexX+kk][iY][iZ][5];
			}
		}
	}
	// interpolate by Y 
	iX=StartIndexX;
	for (ii=0;ii<=n;ii++)
	{
		//initial
		for (int i=0;i<3;i++) tempBY[i][ii]=0.;
		iZ=StartIndexZ+ii;
		for (kk=0;kk<=n;kk++)
		{
			temp=1.;
			for (jj=0;jj<=n;jj++)
			{
				iY=StartIndexY+jj;
				if(jj!=kk)
				{
					temp*=(y-mBField[iX][iY][iZ][1])/(mBField[iX][StartIndexY+kk][iZ][1]-mBField[iX][iY][iZ][1]);
				}
			}
			tempBY[0][ii]+=temp*tempBX[0][kk][ii];			
			tempBY[1][ii]+=temp*tempBX[1][kk][ii];
			tempBY[2][ii]+=temp*tempBX[2][kk][ii];
		}
	}
	// interpolate by Z 
	//initial
	iX=StartIndexX;
	iY=StartIndexY;
	for (int i=0;i<3;i++) B[i]=0.;
	for (kk=0;kk<=n;kk++)
	{
		for(int i=0;i<3;i++ )temp=1.;
		for (jj=0;jj<=n;jj++)
		{
			iZ=StartIndexZ+jj;
			if(jj!=kk)
			{
				temp*=(z-mBField[iX][iY][iZ][2])/(mBField[iX][iY][StartIndexZ+kk][2]-mBField[iX][iY][iZ][2]);
			}
		}
		B[0]+=temp*tempBY[0][kk];			
		B[1]+=temp*tempBY[1][kk];
		B[2]+=temp*tempBY[2][kk];
	}
	return  true;
}

void BField_Tosca_SBS::Rotate_Lab2Field(const double LabP[3],double FieldP[3])
{
	Hep3Vector pFieldP(LabP[0],LabP[1],LabP[2]);
	pFieldP.transform(*mRotL2F);
	FieldP[0]=pFieldP.x();
	FieldP[1]=pFieldP.y();
	FieldP[2]=pFieldP.z();
}

void BField_Tosca_SBS::Transform_Lab2Field(const double LabP[3],double FieldP[3])
{
	for(int i=0;i<3;i++) FieldP[i]=LabP[i]-mOrigin[i];
	if(mDoRotation) Rotate_Lab2Field(FieldP,FieldP);
}


void BField_Tosca_SBS::Rotate_Field2Lab(const double FieldP[3],double LabP[3])
{
	//the following 3 lines do the same thing, but the 3rd line will also change pFieldP
	//pLabP=(*mRotF2L)(pFieldP);
	//pLabP=mRotF2L->operator ()(pFieldP);
	//pLabP=pFieldP.transform(*mRotF2L);
	Hep3Vector pLabP(FieldP[0],FieldP[1],FieldP[2]);
	pLabP.transform(*mRotF2L);
	LabP[0]=pLabP.x();
	LabP[1]=pLabP.y();
	LabP[2]=pLabP.z();
}

void BField_Tosca_SBS::Transform_Field2Lab(const double FieldP[3],double LabP[3])
{
	Rotate_Field2Lab(FieldP,LabP);
	LabP[0]+=mOrigin[0];
	LabP[1]+=mOrigin[1];
	LabP[2]+=mOrigin[2];
}

////////////////////////////////////////////////////////////////////////////////////
void BField_Tosca_SBS::DoFlip(double pPos[], double flag[])
{
	/*
	##################################
	#In most cases the map only cover 1 or 2 sectors of the whole space (in order to minimize map size)
	#We need to use the sysmmetry to flip the sign. Here it shows how to setup the field direction sysmetry
	#Tosca_FlipX=-1 means anti-sysmetry in Bx, the field direction in -x is in the opposite to that in +x  
	#Tosca_FlipX= 0 means NO need to flip sign in Bx, the map already covers the whole range along the X axis  
	#Tosca_FlipX= 1 means mirror sysmetry, the field direction in -x is in the same direction as in +x  
	#
	##################################
	#Most magnet is a dipole (or helmholtz coil) or a solinoid, which will be axial sysmmetrical
	#For a RZ map, the field alwaya pointing to z direction: FlipZ=0 or 1, FlipX=FlipY=-1;
	###For a RZ map, assuming only 1st sector is in the map, (x>0,y>0 and z>0) then we should
	# if (z>=0) {if(x<0) Bx*=-1; if(y<0) By*=-1}
	# else if(z<0) {if(x>0) Bx*=-1; if(y>0) By*=-1}
	# or we can write in this way:  if(z*x<0) Bx*=-1; if(z*y<0) By*=-1}
	#Tosca_FlipX=-1;
	#Tosca_FlipY=-1;
	#Tosca_FlipZ=1;
	#Tosca_Direction=3;
	###For a RZ map, assuming only 1st and 5th sector are in the map, (x>0,y>0 and all z) then we should
	# if (z>=0) {if(x<0) Bx*=-1; if(y<0) By*=-1;}
	# else if(z<0) {if(x<0) Bx*=-1; if(y<0) By*=-1;}
	# or we can write in this way:  if(x<0) Bx*=-1; if(y<0) By*=-1;}
	#Tosca_FlipX=-1;
	#Tosca_FlipY=-1;
	#Tosca_FlipZ=0;
	#Tosca_Direction=3;
	###For a RZ map, assuming 4 sectors are in the map, (full x,y but half z) then we should
	# if(z<0) {Bx*=-1; By*=-1;}
	#Tosca_FlipX=0;
	#Tosca_FlipY=0;
	#Tosca_FlipZ=1;
	#Tosca_Direction=3;
	###For a Z axials map, assuming 8 sectors are in the map, (full x,y and z) then we should
	# not do any flipping
	#Tosca_FlipX=0;
	#Tosca_FlipY=0;
	#Tosca_FlipZ=0;
	#Tosca_Direction=3;
	##################################
	*/

	flag[0]=flag[1]=flag[2]=1.0;
	
	if(mFlip[0]==0 && mFlip[1]==0 && mFlip[2]==0) 
	{
		//the map already cover all space 
		;
	}
	else 
	{ 
		//need to flip 1 or 2 directions
		/////////////////////////////////////////////////////////
		//Field always in x direction 
		if(mDirection==1)
		{
			if(mFlip[0]==1)
			{
				if(mFlip[1] && mFlip[2])
				{
					//+/-1 +/-1 1
					if (pPos[0]*pPos[1]<0) {flag[1]*=(mFlip[1]<0)?-1.0:1.0;}
					if (pPos[0]*pPos[2]<0) {flag[2]*=(mFlip[2]<0)?-1.0:1.0;}
				}
				else if(mFlip[1]==0 && mFlip[2]==0)
				{		
					//0 0 1
					if (pPos[1]<0) {flag[1]*=(mFlip[1]<0)?-1.0:1.0;}
					if (pPos[2]<0) {flag[2]*=(mFlip[2]<0)?-1.0:1.0;}
				}
				else
				{
					cout<<"Do not know how to flip this field map! Wrong parameters in field ini file\n";
				}
			}
			else if(mFlip[0]==0)
			{
				if(mFlip[1] && mFlip[2])
				{
					//+/-1 +/-1 1
					if (pPos[1]<0) {flag[1]*=(mFlip[1]<0)?-1.0:1.0;}
					if (pPos[2]<0) {flag[2]*=(mFlip[2]<0)?-1.0:1.0;}
				}
				else if(mFlip[1]==0 && mFlip[2]==0)
				{		
					flag[0]=flag[1]=flag[2]=1.0;
				}
				else
				{
					cout<<"Do not know how to flip this field map! Wrong parameters in field ini file\n";
				}
			} 
		}
		//Field always in y direction  
		else if(mDirection==2)
		{
			if(mFlip[1]==1)
			{
				if(mFlip[0] && mFlip[2])
				{
					//+/-1 +/-1 1
					if (pPos[1]*pPos[0]<0) {flag[0]*=(mFlip[0]<0)?-1.0:1.0;}
					if (pPos[1]*pPos[2]<0) {flag[2]*=(mFlip[2]<0)?-1.0:1.0;}
				}
				else if(mFlip[0]==0 && mFlip[2]==0)
				{		
					//0 0 1
					if (pPos[0]<0) {flag[0]*=(mFlip[0]<0)?-1.0:1.0;}
					if (pPos[2]<0) {flag[2]*=(mFlip[2]<0)?-1.0:1.0;}
				}
				else
				{
					cout<<"Do not know how to flip this field map! Wrong parameters in field ini file\n";
				}
			}
			else if(mFlip[1]==0)
			{
				if(mFlip[0] && mFlip[2])
				{
					//+/-1 +/-1 1
					if (pPos[0]<0) {flag[0]*=(mFlip[0]<0)?-1.0:1.0;}
					if (pPos[2]<0) {flag[2]*=(mFlip[2]<0)?-1.0:1.0;}
				}
				else if(mFlip[0]==0 && mFlip[2]==0)
				{		
					flag[0]=flag[1]=flag[2]=1.0;
				}
				else
				{
					cout<<"Do not know how to flip this field map! Wrong parameters in field ini file\n";
				}
			}
		}
		//Field always in z direction  
		else if(mDirection==3)
		{
			if(mFlip[2]==1)
			{
				if(mFlip[0] && mFlip[1])
				{
				//+/-1 +/-1 1
					if (pPos[2]*pPos[0]<0) {flag[0]*=(mFlip[0]<0)?-1.0:1.0;}
					if (pPos[2]*pPos[1]<0) {flag[1]*=(mFlip[1]<0)?-1.0:1.0;}
				}
				else if(mFlip[0]==0 && mFlip[1]==0)
				{		
				//0 0 1
					if (pPos[0]<0) {flag[0]*=(mFlip[0]<0)?-1.0:1.0;}
					if (pPos[1]<0) {flag[1]*=(mFlip[1]<0)?-1.0:1.0;}
				}
				else
				{
					cout<<"Do not know how to flip this field map! Wrong parameters in field ini file\n";
				}
			}
			else if(mFlip[2]==0)
			{
				if(mFlip[0] && mFlip[1])
				{
					//+/-1 +/-1 1
					if (pPos[0]<0) {flag[0]*=(mFlip[0]<0)?-1.0:1.0;}
					if (pPos[1]<0) {flag[1]*=(mFlip[1]<0)?-1.0:1.0;}
				}
				else if(mFlip[0]==0 && mFlip[1]==0)
				{		
					flag[0]=flag[1]=flag[2]=1.0;
				}
				else
				{
					cout<<"Do not know how to flip this field map! Wrong parameters in field ini file\n";
				}
			} 
		}//end of z map
	}
}

/////////////////////////////////////////////////////////////////////
bool BField_Tosca_SBS::GetBField(double Pos[3],double B[3])
{//input x,y,z in centimeter, return B field in Tesla
	int i;

	if(mUseUniformB==1)
	{
		for (i=0;i<3;i++) B[i]=mUniformB[i];
		return true;
	}

	if(fabs(mCurrentRatio)<1.0E-10)
	{
		for (i=0;i<3;i++) B[i]=0.0;
		return true;
	}

	double pPos[3],pB[3]={0,0,0},flag[3]={1.0,1.0,1.0};
	//shift and rotate the origin to the field coordinate
	if(mDoShift || mDoRotation) Transform_Lab2Field(Pos,pPos);
	else
	{
		for (i=0;i<3;i++) pPos[i]=Pos[i];
	}

	//do flips if necessary
	DoFlip(pPos,flag);
	
	if(!Interpolation(pPos,pB,1)) {B[0]=B[1]=B[2]=0.0; return false;}

#ifdef BFIELD_TOSCA_SBS_DEBUG
	if(Global_Debug_Level>=3)
	{
		printf("Input position(x,y,z) in cm=(%f, %f, %f):==>\n",Pos[0],Pos[1],Pos[2]);
		printf("Map position(x,y,z) in cm=(%f, %f, %f);\n",pPos[0],pPos[1],pPos[2]);
		printf("The raw magnetic field in Tesla without apply Tosca_CurrentRatio:\n");
		printf("(B_x=%7.5f, B_y=%7.5f, B_z=%7.5f);\n",
			pB[0]*=flag[0],pB[1]*=flag[1],pB[2]*=flag[2]);
	}
#endif

	//apply real current ratio and the flag
	for (i=0;i<3;i++) B[i]=pB[i]*mCurrentRatio*flag[i];

#ifdef BFIELD_TOSCA_SBS_DEBUG
	if(Global_Debug_Level>=2)
	{
		printf("The magnetic field in Tesla after apply Tosca_CurrentRatio: \n");
		printf("(B_x=%7.5f, B_y=%7.5f, B_z=%7.5f);\n",B[0],B[1],B[2]);
	}
#endif

	//rotate back to the Hall coordinate
	if(mDoRotation) 
	{
		Rotate_Field2Lab(B,B);

#ifdef BFIELD_TOSCA_SBS_DEBUG
		if(Global_Debug_Level>=2)
		{
			printf("The magnetic field in Tesla after apply Rotation: \n");
			printf("(B_x=%7.5f, B_y=%7.5f, B_z=%7.5f);\n",B[0],B[1],B[2]);
		}
#endif
	}

	return true;
}

/////////////////////////////////////////////////////////////////////
bool BField_Tosca_SBS::GetBField(float fPos[3],float fB[3])
{//input x,y,z in centimeter
	bool status=false;
	double dPos[3],dB[3];
	int i;
	for ( i=0;i<3;i++)
	{
		dPos[i]=(double) fPos[i];
	}
	status=GetBField(dPos,dB);

	for ( i=0;i<3;i++)
	{
		fB[i]=(float) dB[i];
	}
	return status;
}

/////////////////////////////////////////////////////////////////////
//
//void CreateNtuple(const char *rootfile="field.root",int verbose=0,
//				  double xmin=-50.0,double xmax=50.0, double xstep=1.0,
//				  double ymin=-15.0,double ymax=15.0, double ystep=1.0,
//				  double zmin=-100.0,double zmax=100.0, double zstep=1.0);
//
//Create root ntuple: if verbose != 0 will print steps out 
void BField_Tosca_SBS::CreateNtuple(const char *rootfile,int verbose,double xmin,double xmax,double xstep,
								 double ymin,double ymax,double ystep,double zmin,double zmax,double zstep)
{
	double x,y,r,z,Bx,By,Br,Bz,Btot;

	TFile *file=new TFile(rootfile,"RECREATE");
	TTree *field=new TTree("field","field map");
	field->Branch("x",&x,"x/D");
	field->Branch("y",&y,"y/D");
	field->Branch("r",&r,"r/D");
	field->Branch("z",&z,"z/D");
	field->Branch("Bx",&Bx,"Bx/D");
	field->Branch("By",&By,"By/D");
	field->Branch("Br",&Br,"Br/D");
	field->Branch("Bz",&Bz,"Bz/D");
	field->Branch("Btot",&Btot,"Btot/D");

	if(verbose)
	{
		printf("The Magnetic field Map is:\n");
		printf("       x         y        z       Bx       By       Bz        r        Br        Btot\n");
	}

	int pXNum=(int)ceil((xmax-xmin)/xstep)+1;
	int pYNum=(int)ceil((ymax-ymin)/ystep)+1;
	int pZNum=(int)ceil((zmax-zmin)/zstep)+1;

	double pPos[3],pB[3];

	for (int ix=0;ix<pXNum;ix++)
	{
		pPos[0]=x=xmin+xstep*ix;
		for (int iy=0;iy<pYNum;iy++)
		{
			pPos[1]=y=ymin+ystep*iy;
			r=sqrt(x*x+y*y);
			for (int iz=0;iz<pZNum;iz++)
			{
				pPos[2]=z=zmin+zstep*iz;

				//Get The Field value
				GetBField(pPos,pB);
				Bx=pB[0]; By=pB[1]; Bz=pB[2];  //in Tesla

				Br=sqrt(Bx*Bx+By*By);
				Btot=sqrt(Bx*Bx+By*By+Bz*Bz);

				if(fabs(Btot)<1.0E-13) 
				{
					cout<<"***Warning: Zero field at  x="<<x<<"  y="<<y<<"  z="<<z<<endl;
				}

				field->Fill();

				if(verbose)
				{
					printf("%8.3f %8.3f %8.3f %8.6f %8.6f %8.6f %8.3f %8.6f %8.6f\n",
						x,y,z,Bx,By,Bz,r,Br,Btot);
				}
			}
		}
	}

	file->Write();
	file->Close();
	file->Delete();

}
/////////////////////////////////////////////////////////////////////

