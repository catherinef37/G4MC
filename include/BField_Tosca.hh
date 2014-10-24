// BField_Tosca.h: interface for the BField_Tosca class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(BFIELD_TOSCA)
#define BFIELD_TOSCA
#include <math.h>
#include <CLHEP/Vector/Rotation.h>
#include <CLHEP/Vector/ThreeVector.h>
using namespace CLHEP;

class BField_Tosca
{
public:
	//Static method which returns the singleton pointer of this class.
	static BField_Tosca* GetInstance();
private:
	static BField_Tosca* fInstance;

public:
	BField_Tosca(double pMomentum=0,const char *inifile="BField_Tosca.ini",
		const char *mapfile="tosca.table");
	virtual ~BField_Tosca();
	bool GetBField(double Pos[3],double B[3]);
	bool GetBField(float fPos[3],float fB[3]);

private:
	bool ReadIni(const char *filename);
	bool ReadMap(const char *filename);
	//"target_5.1T_0deg_1cm.table" for SANE target field
	bool Interpolation(double Pos[3],double B[3],int n=2);

	void DoFlip(double pPos[], double flag[]);

public:

	void   Rotate_Field2Lab(const double FieldP[3],double LabP[3]);
	void   Transform_Field2Lab(const double FieldP[3],double LabP[3]);
	void   Rotate_Lab2Field(const double LabP[3],double FieldP[3]);
	void   Transform_Lab2Field(const double LabP[3],double FieldP[3]);
	double GetCurrentRatio(){ return mCurrentRatio;}
	void   SetCurrentRatio(double val){ mCurrentRatio=val;}
	bool   IsUniformField(){ return (mUseUniformB==0)?false:true;}
	void   GetUniformField(double pB[]){pB[0]=mUniformB[0];pB[1]=mUniformB[1];pB[2]=mUniformB[2];}

	CLHEP::HepRotation *GetRotation_L2F() {return mRotL2F;}
	void GetEulerAngles_L2F(double &phi, double &theta, double &psi) 
	{		
		phi=mRotL2F->getPhi();theta=mRotL2F->getTheta();psi=mRotL2F->getPsi(); 
	}

	//Create root ntuple: if verbose != 0 will print steps out 
	void CreateNtuple(const char *rootfile="septumfield.root",int verbose=0,
		double xmin=-50.0,double xmax=50.0, double xstep=1.0,
		double ymin=-15.0,double ymax=15.0, double ystep=1.0,
		double zmin=-100.0,double zmax=100.0, double zstep=1.0);

private:
	//dynalic allocated the 4 dimensional array: double mBField[indexX][indexY][indexZ][Variables];
	//mBField[indexX][indexY][indexZ][0] x
	//mBField[indexX][indexY][indexZ][1] y
	//mBField[indexX][indexY][indexZ][2] z
	//mBField[indexX][indexY][indexZ][3] Bx
	//mBField[indexX][indexY][indexZ][4] By
	//mBField[indexX][indexY][indexZ][5] Bz
	double ****mBField;

	//parameters from ini file
	int    mUseUniformB;
	double mUniformB[3];
	double mCurrentRatio;
	double mStepX,mStepY,mStepZ;
	double mXmin,mXmax;
	double mYmin,mYmax;
	double mZmin,mZmax;
	int	   mInterpolateOutOfRange;
	double mOrigin[3]; 
	double mFieldUnit,mPosUnit;
	int    mFirstDataLine;
	int    mNPara;
	int    mFlip[3], mDirection;
	int	   mRotAxis[3];
	double mRotAngle[3];
	double mDefaultMomentum;

	bool   mDoShift, mDoRotation;
	int    mZNum,mXNum,mYNum;
	double mStepR;
	CLHEP::HepRotation *mRotL2F;
	CLHEP::HepRotation *mRotF2L;

};

#endif // !defined(BFIELD_TOSCA)

