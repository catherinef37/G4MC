// This file defines a class HRSTransport.
// This class is used in G2PSim class as HRS model.
// The definition of variables and the list of available models can be found in
//+the comments in the body
// The active model is chosen during constructing.
// 
// History:
//   Jan 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Add correction function.
//   Feb 2013, J.X. Zhang, Add VC support
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <map>

#include "TROOT.h"
#include "TObject.h"

#include "HRSTransBase.hh"
//#include "HRSTransSTD.hh"
#include "hamcPREXTrans.hh"
#include "hamcCREXTrans.hh"
//#include "G2PTrans400016.hh"
//#include "G2PTrans484816.hh"
//#include "G2PTrans484816R00.hh"
//#include "GDHTransSTD.hh"
//#include "GDHTransLargeX0.hh"

#include "HRSTransport.hh"
#include "TMath.h"
#include <iostream>
using namespace std;

//#define DEBUG_HRS_FORWARD
//#define DEBUG_HRS_BACKWARD

#ifdef BuildRootClass
ClassImp(HRSTransport);
#endif

const double kDEG = 3.14159265358979323846/180.0;

HRSTransport::HRSTransport()
    :iModelIndex(0), bIsLeftArm(true), fHRSAngle(5.767*kDEG),
     fModelAngle(5.767*kDEG), pModel(NULL)
{
    mModel.clear();
    mModelIndex.clear();

    RegisterModel();
}

HRSTransport::HRSTransport(const char* name)
    :iModelIndex(0), bIsLeftArm(true), fHRSAngle(5.767*kDEG),
     fModelAngle(5.767*kDEG), pModel(NULL)
{
    mModel.clear();
    mModelIndex.clear();

    RegisterModel();
    pModel = mModel[mModelIndex[name]];
    iModelIndex = mModelIndex[name];
    fModelAngle = pModel->GetAngle();
}

HRSTransport::HRSTransport(int setting)
    :iModelIndex(0), bIsLeftArm(true), pModel(NULL)
{
    mModel.clear();
    mModelIndex.clear();

    RegisterModel();
    pModel = mModel[setting];
    iModelIndex = setting;
    fModelAngle = pModel->GetAngle();
}

HRSTransport::~HRSTransport()
{
    map<int, HRSTransBase*>::iterator it = mModel.begin();
    while (it!=mModel.end()) {
        delete it->second;
        it++;
    }
    mModel.clear();
    mModelIndex.clear();
}

///////////////////////////////////////////////////////////////////////////
// Transport particles through HRS using SNAKE model
// Use iModelIndex to identify which SNAKE model to be used
// 1: No Septa, 12.5 deg
// 11: 484816 with shim, 5.65 deg, 3 cm raster, by JJL 
// 12: 403216 with shim, 5.65 deg, SNAKE Model not ready yet 
// 13: 400016 with shim, 5.65 deg, 3 cm raster, by Min
// Index > 10 means test
// 18: 484816 with shim, 5.76 deg, no raster, by Min
// 20: GDH exp with small X0 version
// 21: GDH exp with large X0 version
// May add more HRS packages later
///////////////////////////////////////////////////////////////////////////

void HRSTransport::RegisterModel()
{
    HRSTransBase* temp;
    /*
    temp = new HRSTransSTD();
    mModelIndex["STD"] = 1;
    mModel[1] = temp;

    temp = new G2PTrans484816();
    mModelIndex["484816"] = 11;
    mModel[11] = temp;

    temp = new G2PTrans400016();
    mModelIndex["400016"] = 13;
    mModel[13] = temp;

    temp = new G2PTrans484816R00();
    mModelIndex["484816R00"] = 18;
    mModel[18] = temp;

    temp = new GDHTransSTD();
    mModelIndex["GDHSTD"] = 20;
    mModel[20] = temp;

    temp = new GDHTransLargeX0();
    mModelIndex["GDHLargeX0"] = 21;
    mModel[21] = temp;
    */

    temp = new hamcPREXTrans();
    mModelIndex["PREX"] = 47;
    mModel[47] = temp;

    temp = new hamcCREXTrans();
    mModelIndex["CREX"] = 48;
    mModel[48] = temp;

    temp = new hamcPREXTrans();
    mModelIndex["PREXFIELD"] = 49;
    mModel[49] = temp;

    temp = new hamcPREXTrans();
    mModelIndex["PREXSEPTUM"] = 50;
    mModel[50] = temp;

    temp = new hamcPREXTrans();
    mModelIndex["PREXDATA"] = 51;
    mModel[51] = temp;

    temp = new hamcPREXTrans();
    mModelIndex["gmp"] = 52;
    mModel[52] = temp;

    temp = new hamcPREXTrans();
    mModelIndex["CREXfield"] = 53;
    mModel[53] = temp;

    temp = new hamcPREXTrans();
    mModelIndex["PREXFIELD4DEGREES"] = 54;
    mModel[54] = temp;

    temp = new hamcPREXTrans();
    mModelIndex["supercrex"] = 55;
    mModel[55] = temp;

}

void HRSTransport::ChangeModel(int setting)
{
    bool found=false;
    map<string, int>::iterator it = mModelIndex.begin();
    while (it!=mModelIndex.end()) {
        if(setting == it->second) {
            found=true; 
            break;
        }
        it++;
    }
    
    if (!found) {
        printf("Can not find SNAKE model %d, no change.\n", setting);
        exit(-1);
    }
    else{
        pModel = mModel[setting];
        iModelIndex = setting;
        fModelAngle = pModel->GetAngle();
        printf("Current SNAKE model is \"%s\", setting=%d\n", (it->first).c_str(), it->second);
    }
}

void HRSTransport::Acceptance_Check(double* vec, double* x_check, double* theta_check, double* y_check, double* phi_check, int* acc_bool){
  pModel->Acceptance(vec, x_check, theta_check, y_check, phi_check, acc_bool);
}

void HRSTransport::Acceptance_Check_C(double* vec, double* x_check, double* theta_check, double* y_check, double* phi_check,int* acc_bool){
  pModel->Acceptance_C(vec, x_check, theta_check, y_check, phi_check, acc_bool);
}

bool HRSTransport::Forward(const double* V5_tg, double* V5_fp)
{
  // Definition of variables
  // V5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
  // V5_fp = {x_fg, theta_fp, y_fp, phi_fp, delta@tg};
  // delta does not change
  double V5[5];
  
  V5[0] = V5_tg[0];
  V5[1] = tan(V5_tg[1]);
  V5[2] = V5_tg[2];
  V5[3] = tan(V5_tg[3]);
  V5[4] = V5_tg[4];
  
  
  /*
    V5[0] = 0;
    V5[1] = 0;
    V5[2] = 0;
    V5[3] = 0;//tan( 7.5 * 3.141592654 / 180.);
    V5[4] = 0;
  *///THIS WAS JUST AN IDIOT TEST
#ifdef DEBUG_HRS_FORWARD
  cout << "in:  " << V5_tg[0] << " " << V5_tg[1] << " " << V5_tg[2] << " " << V5_tg[3] << " " << V5_tg[4] << endl;
#endif
  
#ifdef DEBUG_HRS_FORWARD
  printf("HRSTransport: %e\t%e\t%e\t%e\t%e\n", V5[0], V5[1], V5[2], V5[3], V5[4]);
#endif
  
  bool bGoodParticle=false;
  
  if (bIsLeftArm) {
    //cout << "In left arm." << endl;
    if( iModelIndex == 47 || iModelIndex == 48 || iModelIndex == 50 ){
      pModel->CoordsCorrection(0., V5);
      //}if( iModelIndex == 50 ){
      //do nothing!!
    }else{
      pModel->CoordsCorrection(fHRSAngle-fModelAngle, V5);
    }
#ifdef DEBUG_HRS_FORWARD
    printf("HRSTransport: %e\t%e\t%e\t%e\t%e\n", V5[0], V5[1], V5[2], V5[3], V5[4]);
#endif
    //cout << fHRSAngle << " " << fModelAngle << endl;
    bGoodParticle = pModel->TransLeftHRS(V5);
    //pModel->FPCorrLeft(V5_tg, V5);
  }
  else {
    //cout << "In right arm." << endl;
    if( iModelIndex == 47 || iModelIndex == 48 || iModelIndex == 50){
      pModel->CoordsCorrection( 0, V5);
      //}if( iModelIndex == 50 ){
      //do nothing!!
    }else{
      pModel->CoordsCorrection(fHRSAngle+fModelAngle, V5);
    }
    
#ifdef DEBUG_HRS_FORWARD
    printf("HRSTransport: %e\t%e\t%e\t%e\t%e\n", V5[0], V5[1], V5[2], V5[3], V5[4]);
#endif
    //cout << fHRSAngle << " " << fModelAngle << endl;
    if( iModelIndex != 50 ){
      bGoodParticle = pModel->TransRightHRS(V5);
    }else{
      bGoodParticle = pModel->TransRightHRS_C(V5);
    }
    //pModel->FPCorrRight(V5_tg, V5);
  }
  
  V5_fp[0] = V5[0];
  V5_fp[1] = atan(V5[1]);
  V5_fp[2] = V5[2];
  V5_fp[3] = atan(V5[3]);
  V5_fp[4] = V5[4];
#ifdef DEBUG_HRS_FORWARD
  cout << "out: " << V5_fp[0] << " " << V5_fp[1] << " " << V5_fp[2] << " " << V5_fp[3] << " " << V5_fp[4] << endl;
#endif
  //cout << "And bGoodParticle is " << bGoodParticle << endl;
  return bGoodParticle;
}

bool HRSTransport::Backward(const double* V5_fp, double* V5_tg)
{
    // Definition of variables
    // V5_fp = {x_fp, theta_fp, y_fp, phi_fp, x_tg};
    // V5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
    // x_tg does not change
    
    double V5[5];
    
    V5[0] = V5_fp[0];
    V5[1] = tan(V5_fp[1]);
    V5[2] = V5_fp[2];
    V5[3] = tan(V5_fp[3]);
    V5[4] = V5_fp[4];

#ifdef DEBUG_HRS_BACKWARD
    printf("HRSTransport: %e\t%e\t%e\t%e\t%e\n", V5[0], V5[1], V5[2], V5[3], V5[4]);
#endif
    
    if (bIsLeftArm) {
        pModel->ReconLeftHRS(V5);
        pModel->CoordsCorrection(fModelAngle-fHRSAngle, V5);
    }
    else {
        pModel->ReconRightHRS(V5);
        pModel->CoordsCorrection(-fModelAngle-fHRSAngle, V5);
    }

    V5_tg[0] = V5[0];
    V5_tg[1] = atan(V5[1]);
    V5_tg[2] = V5[2];
    V5_tg[3] = atan(V5[3]);
    V5_tg[4] = V5[4];

    bool bGoodParticle = false;

    if (V5_tg[4]<1.0) bGoodParticle = true;
    
    return bGoodParticle;
}
