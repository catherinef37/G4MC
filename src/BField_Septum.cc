// ********************************************************************
// $Id: BField_Septum.cc,v 3.0, 2011/1/19  G2P Exp $
// Implementation of the BField_Septum class.
//
//////////////////////////////////////////////////////////////////////

#include "BField_Septum.hh"
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
//#define BFIELD_SEPTUM_DEBUG 0

#ifdef BFIELD_SEPTUM_DEBUG
#include "GlobalDebuger.hh"
#endif


using namespace std;

BField_Septum* BField_Septum::fInstance=0;
BField_Septum* BField_Septum::GetInstance()
{ 
	if(!fInstance)  
	{
		//new BField_Septum();
		cout<<"BField_Septum is not nitialized yet...exit...\n";
		exit(-99);
	}
	return fInstance; 
}



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BField_Septum::BField_Septum(double pMomentumL,double pMomentumR,
							 const char *inifile,const char *mapfile)
{
#ifdef BFIELD_SEPTUM_DEBUG
	if(BFIELD_SEPTUM_DEBUG > Global_Debug_Level)
		SetGlobalDebugLevel("BField_Septum::BField_Septum()", (int)BFIELD_SEPTUM_DEBUG);
#endif

	fInstance=this;

	this->ReadIni(inifile);
	if(fabs(mDefaultMomentumL)<1.0E-5 || fabs(mDefaultMomentumR)<1.0E-5) 
	{
		cerr<<"\n##DefaultMomentumL is invaid ("<<mDefaultMomentumL<<", "<<mDefaultMomentumR
			<<") in BField_Septum.ini. I quit... \n";
		exit(-1);
	}

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
	*mRotL2F=pRot[2]*pRot[1]*pRot[0];
	*mRotF2L=mRotL2F->inverse(); 
#ifdef BFIELD_SEPTUM_DEBUG
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
	if(fabs(pMomentumL)>1.0E-5 || fabs(pMomentumR)>1.0E-5) 
		SetMomentum(pMomentumL,pMomentumR);
}

BField_Septum::~BField_Septum()
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
void   BField_Septum::SetMomentum(double pMomentumL,double pMomentumR)
{ 
	SetMomentumL(pMomentumL);
	SetMomentumR(pMomentumR);
}

void   BField_Septum::SetMomentumL(double pMomentumL)
{ 
	//update the Current Ratio if necessary
	double RatioL=pMomentumL/mDefaultMomentumL;
	if(fabs(RatioL)<0.00001) RatioL=0.0;
	SetCurrentRatioL(RatioL);
}

void   BField_Septum::SetMomentumR(double pMomentumR)
{ 
	//update the Current Ratio if necessary
	double RatioR=pMomentumR/mDefaultMomentumR;
	if(fabs(RatioR)<0.00001) RatioR=0.0;
	SetCurrentRatioR(RatioR);
}


/////////////////////////////////////////////////////////////////////
void   BField_Septum::SetCurrentRatio(double valL, double valR)
{ 
	SetCurrentRatioL(valL);
	SetCurrentRatioR(valR);
}	

void   BField_Septum::SetCurrentRatioL(double valL)
{ 
	if(fabs(mCurrentRatioL-valL)>1.0E-05)
	{
		mCurrentRatioL=valL;
		cout<<"\n##BField_Septum::SetCurrentRatioL(rL): Left septum current ratio is set to "
			<<mCurrentRatioL<<endl;
		
		UsageManager *pConfig=UsageManager::GetUsageManager();
		//all these 3 methods work
		//char tmpStr[20];
		//sprintf(tmpStr,"%.6f",mCurrentRatioL);
		//std::string theStr(tmpStr);
		//pConfig->SetParameter("Septum_CurrentRatioL",tmpStr);  
		//pConfig->SetParameter("Septum_CurrentRatioL",theStr);
		pConfig->SetParameter("Septum_CurrentRatioL",mCurrentRatioL);
	}
	
}
void   BField_Septum::SetCurrentRatioR(double valR)
{ 
	if(fabs(mCurrentRatioR-valR)>1.0E-05)
	{
		mCurrentRatioR=valR;
		cout<<"\n##BField_Septum::SetCurrentRatioR(rR): Right septum current ratio is set to "
			<<mCurrentRatioR<<endl;
		
		UsageManager *pConfig=UsageManager::GetUsageManager();
		pConfig->SetParameter("Septum_CurrentRatioR",mCurrentRatioR);
	}
}


/////////////////////////////////////////////////////////////////////
bool BField_Septum::ReadIni(const char *filename)
{
	double deg2rad=acos(-1.0)/180.;
	//By Jixie: I am not use this routine to read ini file any longer
	//I prefer to use UsageManager::ReadFile
	//
	UsageManager *pConfig=UsageManager::GetUsageManager();
	bool ret=pConfig->ReadFile(filename);

	pConfig->GetParameter("Septum_UseUniformB",mUseUniformB);
	pConfig->GetParameter("Septum_UniformB_x", mUniformB[0]);
	pConfig->GetParameter("Septum_UniformB_y", mUniformB[1]);
	pConfig->GetParameter("Septum_UniformB_z", mUniformB[2]);
	
	pConfig->GetParameter("Septum_FieldUnit",	mFieldUnit);
	pConfig->GetParameter("Septum_FirstDataLine",mFirstDataLine);
	pConfig->GetParameter("Septum_NPara",		mNPara);

	pConfig->GetParameter("Septum_Xmin",		mXmin);
	pConfig->GetParameter("Septum_Xmax",		mXmax);
	pConfig->GetParameter("Septum_Ymin",		mYmin);
	pConfig->GetParameter("Septum_Ymax",		mYmax);
	pConfig->GetParameter("Septum_Zmin",		mZmin);
	pConfig->GetParameter("Septum_Zmax",		mZmax);
	pConfig->GetParameter("Septum_StepX",		mStepX);
	pConfig->GetParameter("Septum_StepY",		mStepY);
	pConfig->GetParameter("Septum_StepZ",		mStepZ);
	pConfig->GetParameter("Septum_InterpolateOutOfRange", mInterpolateOutOfRange);

	pConfig->GetParameter("Septum_OriginX",		mOrigin[0]);
	pConfig->GetParameter("Septum_OriginY",		mOrigin[1]);
	pConfig->GetParameter("Septum_OriginZ",		mOrigin[2]);
	pConfig->GetParameter("Septum_RotAxis1",	mRotAxis[0]);
	pConfig->GetParameter("Septum_RotAxis2",	mRotAxis[1]);
	pConfig->GetParameter("Septum_RotAxis3",	mRotAxis[2]);
	pConfig->GetParameter("Septum_RotAngle1",	mRotAngle[0]); mRotAngle[0]*=deg2rad;
	pConfig->GetParameter("Septum_RotAngle2",	mRotAngle[1]); mRotAngle[1]*=deg2rad;
	pConfig->GetParameter("Septum_RotAngle3",	mRotAngle[2]); mRotAngle[2]*=deg2rad;

	pConfig->GetParameter("Septum_DefaultMomentumL",mDefaultMomentumL);
	pConfig->GetParameter("Septum_CurrentRatioL",mCurrentRatioL);	
	pConfig->GetParameter("Septum_DefaultMomentumR",mDefaultMomentumR);
	pConfig->GetParameter("Septum_CurrentRatioR",mCurrentRatioR);	

#ifdef BFIELD_SEPTUM_DEBUG
	pConfig->PrintParamMap(); 
#endif
	return ret;

	//The following is the original code
	//please keep it this way
	FILE *ini;
	char ch[]="=;\n";
	char *name,*value,line[100];
	char *pDest;
	int iPos=-1;

	if((ini=fopen(filename,"r"))==NULL)
	{
		printf("***Error! Can not open configure file \"%s\"!",filename);
		return false;
	}

	printf("\nThe magnetic field configuration is:\n");
	while(!feof(ini))
	{
		fgets(line,100,ini);
		/* Search forward. for the '#' to skip the comment*/
		pDest = strchr( line, '#' );
		iPos = pDest - line + 1;
		//in Linux, if not found '#', iPos==1073747457//
		if(iPos>0 && iPos<100) continue;
		name=strtok(line,ch);
		value=strtok(0,ch);
		//show the confif info

		if(name&&value) printf("%15s = %s\n",name,value);
		else printf("read %s error\n",filename);

		if (strcmp(name,"UseUniformB")==0)				mUseUniformB=atoi(value);
		else if (strcmp(name,"Septum_UniformB_x")==0)	mUniformB[0]=atof(value);
		else if (strcmp(name,"Septum_UniformB_y")==0)	mUniformB[1]=atof(value);
		else if (strcmp(name,"Septum_UniformB_z")==0)	mUniformB[2]=atof(value);
		
		else if (strcmp(name,"Septum_FieldUnit")==0)	mFieldUnit=atof(value);
		else if (strcmp(name,"Septum_FirstDataLine")==0)mFirstDataLine=atoi(value);
		else if (strcmp(name,"Septum_NPara")==0)		mNPara=atoi(value);

		else if (strcmp(name,"Septum_Xmin")==0)			mXmin=atof(value);
		else if (strcmp(name,"Septum_Xmax")==0)			mXmax=atof(value);
		else if (strcmp(name,"Septum_Ymin")==0)			mYmin=atof(value);
		else if (strcmp(name,"Septum_Ymax")==0)			mYmax=atof(value);
		else if (strcmp(name,"Septum_Zmin")==0)			mZmin=atof(value);
		else if (strcmp(name,"Septum_Zmax")==0)			mZmax=atof(value);
		else if (strcmp(name,"Septum_StepX")==0)		mStepX=atof(value);
		else if (strcmp(name,"Septum_StepY")==0)		mStepY=atof(value);
		else if (strcmp(name,"Septum_StepZ")==0)		mStepZ=atof(value);
		else if (strcmp(name,"Septum_InterpolateOutOfRange")==0) mInterpolateOutOfRange=atoi(value);

		else if (strcmp(name,"Septum_OriginX")==0)		mOrigin[0]=atof(value);
		else if (strcmp(name,"Septum_OriginY")==0)		mOrigin[1]=atof(value);
		else if (strcmp(name,"Septum_OriginZ")==0)		mOrigin[2]=atof(value);
		else if (strcmp(name,"Septum_RotAxis1")==0)		mRotAxis[0]=atoi(value);
		else if (strcmp(name,"Septum_RotAxis2")==0)		mRotAxis[1]=atoi(value);
		else if (strcmp(name,"Septum_RotAxis3")==0)		mRotAxis[2]=atoi(value);
		else if (strcmp(name,"Septum_RotAngle1")==0)	mRotAngle[0]=atof(value)*deg2rad;
		else if (strcmp(name,"Septum_RotAngle2")==0)	mRotAngle[1]=atof(value)*deg2rad;
		else if (strcmp(name,"Septum_RotAngle3")==0)	mRotAngle[2]=atof(value)*deg2rad;

		else if (strcmp(name,"Septum_DefaultMomentumL")==0)	mDefaultMomentumL=atof(value);
		else if (strcmp(name,"Septum_CurrentRatioL")==0)	mCurrentRatioL=atof(value);
		else if (strcmp(name,"Septum_DefaultMomentumR")==0)	mDefaultMomentumR=atof(value);
		else if (strcmp(name,"Septum_CurrentRatioR")==0)	mCurrentRatioR=atof(value);
		else continue;
	}
	fclose(ini);

	if(mNPara<6) mNPara=6; //mNPara should not less than 6 columns
	return true;

}

/////////////////////////////////////////////////////////////////////
bool BField_Septum::ReadMap(const char *filename)
{
	char strLog[1024];
	sprintf(strLog,"BField_Septum::ReadMap() is loading field map %s......\n",filename);
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

	//eat the first 15 lines
	int nLine2Eat=mFirstDataLine-1;
	int LineNum=nLine2Eat;
	char tempname[256];
	for(int i=0;i<nLine2Eat;i++) ins.getline (tempname,256);

	while (!ins.eof())
	{
		ins.getline(tempname,256);	
		LineNum++;
		//check if it is an empty line	
		if(strlen(tempname) < size_t(2*mNPara-1)) continue;

		istringstream s1(tempname);
		for(col=0;col<mNPara;col++)    s1>>tempLine[col];

		//check for X, Y and Z
		if (tempLine[0]>=mXmin && tempLine[0]<=mXmax && tempLine[1]>=mYmin && tempLine[1]<=mYmax &&
			tempLine[2]>=mZmin && tempLine[2]<=mZmax)
		{//store the value

			//in case there is an empty line, r=Btot=0.0
			if( sqrt(tempLine[0]*tempLine[0]+tempLine[1]*tempLine[1]+tempLine[2]*tempLine[2]) < 1.0E-8 && 
				sqrt(tempLine[3]*tempLine[3]+tempLine[4]*tempLine[4]+tempLine[5]*tempLine[5]) < 1.0E-8 ) 
			{
				cout<<"***Warning: ZERO field at  x="<<tempLine[0]<<"  y="<<tempLine[1]
				<<"  z="<<tempLine[2]<<endl;
				cout<<"***Line "<<LineNum<<" could be an empty line in the map "<<filename<<endl;
			}
			else 
			{
				indexX=int((tempLine[0]-mXmin)/mStepX);
				indexY=int((tempLine[1]-mYmin)/mStepY);
				indexZ=int((tempLine[2]-mZmin)/mStepZ);
				for(col=0;col<mNPara;col++)
				{
					mBField[indexX][indexY][indexZ][col]=tempLine[col];
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


#ifdef BFIELD_SEPTUM_DEBUG
	if(Global_Debug_Level>=4)
	{
		printf("The Magnetic field Map is:\n");		
		printf("       x         y        z       Bx       By       Bz        r        Br        Btot\n");
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

				if(fabs(r)<1.0E-8 && fabs(z)<1.0E-08 && fabs(Btot)<1.0E-8) 
				{
					cout<<"***Warning: ZERO field at  x="<<x<<"  y="<<y<<"  z="<<z<<endl;
					cout<<"***Warning: There is an empty line in the map buffer..."<<endl;
				}
				else field->Fill();

#ifdef BFIELD_SEPTUM_DEBUG
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

#ifdef BFIELD_SEPTUM_DEBUG
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
bool BField_Septum::Interpolation(double Pos[3],double B[3],int n)
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

#ifdef BFIELD_SEPTUM_DEBU
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

	//just flip by z
	double x=Pos[0],y=Pos[1],z=Pos[2];

	int StartIndexX=0,StartIndexY=0,StartIndexZ=0;  //the first point to do the interpolation
	StartIndexX=int((x-mBField[0][0][0][0])/mStepX);
	StartIndexY=int((y-mBField[0][0][0][1])/mStepY);
	StartIndexZ=int((z-mBField[0][0][0][2])/mStepZ);


#ifdef BFIELD_SEPTUM_DEBUG
	if (StartIndexX<0 || StartIndexX>=mXNum)
	{
		if(BFIELD_SEPTUM_DEBUG>=1)
		printf("\n**Warning!!! X=%.4f is out of range [%.4f,%.4f) !!!**\n",
			x,mBField[0][0][0][0],mBField[mXNum-1][0][0][0]);
		B[0]=B[1]=B[2]=0.0;return false;
	}
	if (StartIndexY<0 || StartIndexY>=mYNum)
	{
		if(BFIELD_SEPTUM_DEBUG>=1)
		printf("\n**Warning!!! Y=%.4f is out of range [%.4f,%.4f) !!!**\n",
			y,mBField[0][0][0][1],mBField[0][mYNum-1][0][1]);
		B[0]=B[1]=B[2]=0.0;return false;
	}
	if (StartIndexZ<0 || StartIndexZ>=mZNum)
	{
		if(BFIELD_SEPTUM_DEBUG>=1)
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

void BField_Septum::Rotate_Lab2Field(const double LabP[3],double FieldP[3])
{
	Hep3Vector pFieldP(LabP[0],LabP[1],LabP[2]);
	pFieldP.transform(*mRotL2F);
	FieldP[0]=pFieldP.x();
	FieldP[1]=pFieldP.y();
	FieldP[2]=pFieldP.z();
}

void BField_Septum::Transform_Lab2Field(const double LabP[3],double FieldP[3])
{
	for(int i=0;i<3;i++) FieldP[i]=LabP[i]-mOrigin[i];
	if(mDoRotation) Rotate_Lab2Field(FieldP,FieldP);
}


void BField_Septum::Rotate_Field2Lab(const double FieldP[3],double LabP[3])
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

void BField_Septum::Transform_Field2Lab(const double FieldP[3],double LabP[3])
{
	Rotate_Field2Lab(FieldP,LabP);
	LabP[0]+=mOrigin[0];
	LabP[1]+=mOrigin[1];
	LabP[2]+=mOrigin[2];
}


/////////////////////////////////////////////////////////////////////
bool BField_Septum::GetBField(double Pos[3],double B[3])
{//input x,y,z in centimeter, return B field in Tesla
	int i;

	if(mUseUniformB==1)
	{
		for (i=0;i<3;i++) B[i]=mUniformB[i];
		return true;
	}

	double pPos[3],pB[3]={0,0,0},flag[3]={1.0,1.0,1.0};
	//shift and rotate the origin to the field coordinate
	if(mDoShift || mDoRotation) Transform_Lab2Field(Pos,pPos);
	else
	{
		for (i=0;i<3;i++) pPos[i]=Pos[i];
	}

	if( (fabs(mCurrentRatioL)<1.0E-8 && pPos[0]>=0) || 
		(fabs(mCurrentRatioR)<1.0E-8 && pPos[0]<0) )
	{
		for (i=0;i<3;i++) B[i]=0.0;
		return true;
	}

	//the map provide fields for the whole range of x and y, but only half of z
	//field is y direction. Need to flip for -z 
	//if z<0, flip Bz. Bx does not need to flip 
	if (pPos[2]<0) {flag[2]*=-1.0;}

	//Note that the map is only covers z>0, I have to take the absolute values
	double pPosAtMap[]={pPos[0],pPos[1],fabs(pPos[2])};
	if(!Interpolation(pPosAtMap,pB,1)) {B[0]=B[1]=B[2]=0.0; return false;}

#ifdef BFIELD_SEPTUM_DEBUG
	if(Global_Debug_Level>=3)
	{
		printf("Input position(x,y,z) in cm=(%f, %f, %f):==>\n",Pos[0],Pos[1],Pos[2]);
		printf("Map position(x,y,z) in cm=(%f, %f, %f);\n",pPos[0],pPos[1],pPos[2]);
		printf("The raw magnetic field in Tesla without apply Septum_CurrentRatio:\n");
		printf("(B_x=%7.5f, B_y=%7.5f, B_z=%7.5f);\n",
			pB[0]*=flag[0],pB[1]*=flag[1],pB[2]*=flag[2]);
	}
#endif

	//apply real current ratio and the flag
	for (i=0;i<3;i++) 
	{
		if(pPos[0]>=0) B[i]=pB[i]*mCurrentRatioL*flag[i];
		else B[i]=pB[i]*mCurrentRatioR*flag[i];
	}
#ifdef BFIELD_SEPTUM_DEBUG
	if(Global_Debug_Level>=2)
	{
		printf("The magnetic field in Tesla after apply Septum_CurrentRatio: \n");
		printf("(B_x=%7.5f, B_y=%7.5f, B_z=%7.5f);\n",B[0],B[1],B[2]);
	}
#endif

	//rotate back to the Hall coordinate
	if(mDoRotation) 
	{
		Rotate_Field2Lab(B,B);

#ifdef BFIELD_SEPTUM_DEBUG
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
bool BField_Septum::GetBField(float fPos[3],float fB[3])
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
//void CreateNtuple(const char *rootfile="septumfield.root",int verbose=0,
//				  double xmin=-43.0,double xmax=43.0, double xstep=1.0,
//				  double ymin=-12.0,double ymax=12.0, double ystep=1.0,
//				  double zmin=-100.0,double zmax=100.0, double zstep=1.0);
//
//Create root ntuple: if verbose != 0 will print steps out 
void BField_Septum::CreateNtuple(const char *rootfile,int verbose,double xmin,double xmax,double xstep,
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

