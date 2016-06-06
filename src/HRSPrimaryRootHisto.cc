// ********************************************************************
//
// $Id: HRSPrimaryRootHisto.cc,v 1.0, 2008/8/29  HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
// --------------------------------------------------------------
//
#include "HRSPrimaryRootHisto.hh"
#include "TRandom.h"

#include "GlobalDebuger.hh"
//#define HRSPRIMARYROOTHISTO_DEBUG 3

HRSPrimaryRootHisto::HRSPrimaryRootHisto()
{
#ifdef HRSPRIMARYROOTHISTO_DEBUG
    if(Global_Debug_Level < HRSPRIMARYROOTHISTO_DEBUG) 
    {
        Global_Debug_Level = HRSPRIMARYROOTHISTO_DEBUG;
        cout<<"HRSPrimaryRootHisto(): Set Global_Debug_Level to "<<Global_Debug_Level<<endl;
    }
#endif
    gRandom->SetSeed(0);
    bHistoOpened=false;
    mBeamE=0.0;
    mTargetMass=0.0;
    //no need to initialize it since we don't know what kind of event the user will use
    //Initialize();

    //cout<<"HRSPrimaryRootHisto() construction done!"<<endl;
}

HRSPrimaryRootHisto::~HRSPrimaryRootHisto()
{
    if(mFile) {mFile->Close();mFile->Delete();} 
    //cout<<"delete HRSPrimaryRootHisto ... done!"<<endl;
}

void HRSPrimaryRootHisto::Initialize()
{ //beam energy and target came from the histo itself

    mFile=TFile::Open("PrimaryEvent.root");
    if (!mFile) 
    {
        cout<<"\nError!!! Cann't load PrimaryEvent.root, File not found!"<<endl;
        return;
    }
    //you need to have 7 histo: 

    //get histogram
    //TH1F *hE,*hT;
    //TH2F *h2ThetaPhi_el,*h2WprimeQ2_Inc,*h2WprimeQ2_Exc,*h2ThetaZ_rtpc;
    //TH3F *h3P_rtpc,*h3P_el,*h3P_pi,*h3P_pr; //Phi(x) Theta(y) P(z) 
    if(mFile) 
    {
        hE = (TH1F*)mFile->Get("hE");
        hT = (TH1F*)mFile->Get("hT");
        h2WprimeQ2_Inc = (TH2F*)mFile->Get("h2WprimeQ2_Inc");
        h2WprimeQ2_Exc = (TH2F*)mFile->Get("h2WprimeQ2_Exc");
        h2ThetaPhi_el = (TH2F*)mFile->Get("h2ThetaPhi_el");
        h2ThetaZ_rtpc = (TH2F*)mFile->Get("h2ThetaZ_rtpc");
        h3P_rtpc = (TH3F*)mFile->Get("h3P_rtpc");
        h3P_el = (TH3F*)mFile->Get("h3P_el");
        h3P_pi = (TH3F*)mFile->Get("h3P_pi");
        h3P_pr = (TH3F*)mFile->Get("h3P_pr");
        bHistoOpened=true;
    }
    if(hE) mBeamE=hE->GetMean();
    if(hT) mTargetMass=hT->GetMean();

    cout<<"HRSPrimaryRootHisto:  Beam Energy="<<mBeamE<<" TargetMass="<<mTargetMass<<endl;

    if (mBeamE<0.9) cout<<"\nError!!! Beam Energy is too low!! mBeamE="<<mBeamE<<endl;
    if (mTargetMass<0.9) cout<<"\nError!!! Target Mass is too small!! mTargetMass="<<mTargetMass<<endl;

}

//e(D,e'p_rtpc)X
void HRSPrimaryRootHisto::GetInclusiveEvent(TVector3 &Vertex, TVector3 &V3P_rtpc, TVector3 &V3P_el)
{  
    double deg2rad = 3.14159/180.0;
    const double E=this->mBeamE;	//GeV
    const double M=this->mTargetMass; //GeV

    TLorentzVector  V4El, V4X, V4Ps, V4T, V4Beam;
    V4Beam.SetXYZM(0,0,E,0.000511);
    V4T.SetXYZM(0,0,0,M+0.9396);

    // zero vector can't be stretched, I need to set TVector3 non-zero before using SetMag()
    Vertex.SetXYZ(0,0,1);
    V3P_rtpc.SetXYZ(0,0,1);
    V3P_el.SetXYZ(0,0,1);
    double x=0.0,y=0.0,z=0.0;

    bool bRTPCCut=false;
    const double P_RTPC_Limit = 0.065;
    double W,Q2,P_el_tot,Theta_deg_el,Phi_deg_el;

    if(!bRTPCCut)
    {
        h3P_rtpc->GetRandom3(x,y,z); //Phi(x) Theta(y) P(z) 
        V3P_rtpc.SetMag(z);
        V3P_rtpc.SetTheta(y*deg2rad);
        V3P_rtpc.SetPhi(x*deg2rad);
    }
    else
    {
        //get P_rtpc ,h3P_rtpc is in GeV unit
        //need to verify that this Ptot will pass the thresould  

        //This cut should not be applied in order to do acceptance study
        while ( z < P_RTPC_Limit )
        {	  
            h3P_rtpc->GetRandom3(x,y,z); //Phi(x) Theta(y) P(z) 		 
            if( z > P_RTPC_Limit && y>15.0 && y<165.0 &&  z > RTPC_P_Threshold(y)) 
            {
                V3P_rtpc.SetMag(z);
                V3P_rtpc.SetTheta(y*deg2rad);
                V3P_rtpc.SetPhi(x*deg2rad);
            }
            else
            {
                z=0.0; //try again
            }
        }   
    }  
#ifdef HRSPRIMARYROOTHISTO_DEBUG
    if(Global_Debug_Level>=4)
    {
        cout<<"GetInclusiveEvent() throw V3P_rtpc: Phi="<<x<<"deg  "
            <<"Theta= "<<y<<"deg  Ptot= "<<z<<"GeV  ==> P_threshold="<<RTPC_P_Threshold(y)<<endl;
    }
#endif
    V4Ps.SetXYZM(V3P_rtpc.x(),V3P_rtpc.y(),V3P_rtpc.z(),0.9383); 

    //get vertex 
    x=y=z=0.0;
    GetXbyY(h2ThetaZ_rtpc,V3P_rtpc.Theta()/deg2rad,z);
    Vertex.SetXYZ(x,y,z);

REPEAT0:
    //Get the P for e' 
    P_el_tot=0.0;
    while (P_el_tot < 0.15 * E || P_el_tot > 0.95 * E)
    {
        h2WprimeQ2_Inc->GetRandom2(Q2,W);
        P_el_tot = E - (W*W-M*M-Q2)/(2.0*M); 
#ifdef HRSPRIMARYROOTHISTO_DEBUG
        if(Global_Debug_Level>=4)
        {
            cout<<"GetInclusiveEvent() throw el: W="<<W<<"GeV  "
                <<"Q2= "<<Q2<<"GeV^2  E'=P_el_tot= "<<P_el_tot<<"GeV  "<<endl;
        }
#endif
    } 
    Theta_deg_el=acos(1.0-Q2/(4.0*E*P_el_tot))/deg2rad;
    GetXbyY(h2ThetaPhi_el,Theta_deg_el,Phi_deg_el);  
    V3P_el.SetMag(P_el_tot);
    V3P_el.SetTheta(Theta_deg_el*deg2rad);
    V3P_el.SetPhi(Phi_deg_el*deg2rad);

    V4El.SetXYZM(V3P_el.x(),V3P_el.y(),V3P_el.z(),0.000511); 

    V4X=V4Beam+V4T-V4El-V4Ps;
    if( V4X.Mag()<0.2 || V4X.Energy() < 0.9383) goto REPEAT0;

#ifdef HRSPRIMARYROOTHISTO_DEBUG
    if(Global_Debug_Level>=2)
    {
        cout<<"HRSPrimaryRootHisto::GetInclusiveEvent()  BeamE="<<E<<"  TargetMass="<<M<<"\n"
            <<"  Vertex="<<Vertex<<"  V3P_rtpc="<<V3P_rtpc<<"  V3P_el="<<V3P_el<<endl;
    }
#endif
    return;
}

//e(D,e'pi-p_rtpc)p_clas
void HRSPrimaryRootHisto::GetExclusiveEvent(TVector3 &Vertex, TVector3 &V3P_rtpc, TVector3 &V3P_el,
                                              TVector3 &V3P_pi, TVector3 &V3P_clas)
{
    //for exclusive events: (1) get rtpc monetum, then using Theta Vs Z of rtpc to detemind vertex
    //(2) get scattered electron from Wprime and Q2
    //(3) get pion minus momentum for pi- histo
    //(4) based on missing mass to get the clas proton momenum, relolution P=min(3%,50MeV) 

    double deg2rad = 3.14159/180.0;
    const double E=this->mBeamE; 
    const double M=this->mTargetMass;

    TLorentzVector  V4El, V4Pi, V4P, V4Ps, V4T, V4Beam;
    V4Beam.SetXYZM(0,0,E,0.000511);
    V4T.SetXYZM(0,0,0,M+0.9396);

    // zero vector can't be stretched, I need to set TVector3 non-zero before using SetMag()
    Vertex.SetXYZ(0,0,1);
    V3P_rtpc.SetXYZ(0,0,1);
    V3P_el.SetXYZ(0,0,1);
    V3P_pi.SetXYZ(0,0,1);
    V3P_clas.SetXYZ(0,0,1);
    double x,y,z;


    bool bRTPCCut=false;
    double P_pi_tot=0.0;
    const double P_RTPC_Limit = 0.065;
    double W,Q2,P_el_tot,Theta_deg_el,Phi_deg_el;
    int REPEAT2_counter=0;

REPEAT1:

    if(!bRTPCCut)
    {
        h3P_rtpc->GetRandom3(x,y,z); //Phi(x) Theta(y) P(z) 
        V3P_rtpc.SetMag(z);
        V3P_rtpc.SetTheta(y*deg2rad);
        V3P_rtpc.SetPhi(x*deg2rad);
    }
    else
    {
        //get P_rtpc ,h3P_rtpc is in GeV unit
        //need to verify that this Ptot will pass the thresould  
        x=y=z=0.0;
        while ( z < P_RTPC_Limit )
        {	  
            h3P_rtpc->GetRandom3(x,y,z); //Phi(x) Theta(y) P(z) 	
            if( z > P_RTPC_Limit && y>15.0 && y<165.0 &&  z > RTPC_P_Threshold(y)) 
            {
                V3P_rtpc.SetMag(z);
                V3P_rtpc.SetTheta(y*deg2rad);
                V3P_rtpc.SetPhi(x*deg2rad);
            }
            else
            {
                z=0.0; //try again
            }
        }   
    }  
#ifdef HRSPRIMARYROOTHISTO_DEBUG
    if(Global_Debug_Level>=4)
    {
        cout<<"GetInclusiveEvent() throw V3P_rtpc: Phi="<<x<<"deg  "
            <<"Theta= "<<y<<"deg  Ptot= "<<z<<"GeV  ==> P_threshold="<<RTPC_P_Threshold(y)<<endl;
    }
#endif
    V4Ps.SetXYZM(V3P_rtpc.x(),V3P_rtpc.y(),V3P_rtpc.z(),0.9383); 

    //get vertex 
    x=y=z=0.0;
    GetXbyY(h2ThetaZ_rtpc,V3P_rtpc.Theta()/deg2rad,z);
    Vertex.SetXYZ(x,y,z);

    //Get the P for e' 
    P_el_tot=0.0;
    while (P_el_tot < 0.15 * E || P_el_tot > 0.95 * E)
    {
        h2WprimeQ2_Exc->GetRandom2(Q2,W);
        P_el_tot = E - (W*W-M*M-Q2)/(2.0*M); 
#ifdef HRSPRIMARYROOTHISTO_DEBUG
        if(Global_Debug_Level>=4)
        {
            cout<<"GetExclusiveEvent() throw el: W="<<W<<"GeV  "
                <<"Q2= "<<Q2<<"GeV^2  E'=P_el_tot= "<<P_el_tot<<"GeV  "<<endl;
        }
#endif
    } 
    Theta_deg_el=acos(1.0-Q2/(4.0*E*P_el_tot))/deg2rad;
    GetXbyY(h2ThetaPhi_el,Theta_deg_el,Phi_deg_el);  
    V3P_el.SetMag(P_el_tot);
    V3P_el.SetTheta(Theta_deg_el*deg2rad);
    V3P_el.SetPhi(Phi_deg_el*deg2rad);

    V4El.SetXYZM(V3P_el.x(),V3P_el.y(),V3P_el.z(),0.000511); 

    V4P=V4Beam+V4T-V4El-V4Ps;
    if( V4P.Mag()<0.2 || V4P.Energy() < 0.1396+0.9383) goto REPEAT1;


    REPEAT2_counter=0;
REPEAT2:
    //P_pi
    P_pi_tot=0;
    while (P_pi_tot<0.2)
    {
        h3P_pi->GetRandom3(x,y,z);
#ifdef HRSPRIMARYROOTHISTO_DEBUG
        if(Global_Debug_Level>=4)
        {
            cout<<"GetExclusiveEvent() throw V3P_pi-: Phi="<<x<<"deg  "
                <<"Theta= "<<y<<"deg  Ptot= "<<z<<"GeV  "<<endl;
        }
#endif
        if(z > 0.2) V3P_pi.SetMag(z);
        V3P_pi.SetTheta(y*deg2rad);
        V3P_pi.SetPhi(x*deg2rad);
        P_pi_tot=z;
    }
    V4Pi.SetXYZM(V3P_pi.x(),V3P_pi.y(),V3P_pi.z(),0.1396); 
    REPEAT2_counter++;

    //using missing mass to determind the P_clas
    V4P=V4Beam+V4T-V4El-V4Pi-V4Ps;
    if( V4P.Mag()<0.2 || fabs(V4P.M()-0.9383)>0.05 ) 
    {
        if(REPEAT2_counter<100) goto REPEAT2;
        else goto REPEAT1;
    }
    V3P_clas=V4P.Vect();
    /*
    //smear the total momentum, but Gail recommend not to use this feature
    double Mean=V4P.P();
    double Sigma=(Mean*0.03 < 0.05) ? Mean*0.03 : 0.05;
    TF1 *f=new TF1("Gaussian","gaus",Mean-2.0*Sigma,Mean+2.0*Sigma);
    f->SetParameters(100.0,Mean,Sigma);
    //throw the total momentum
    double P_pr_tot=f->GetRandom();  
    V3P_clas.SetMag(P_pr_tot);
    V3P_clas.SetTheta(V4P.Theta());
    V3P_clas.SetPhi(V4P.Phi());
    #ifdef HRSPRIMARYROOTHISTO_DEBUG
    if(Global_Debug_Level>=3)
    {
    cout<<"GetExclusiveEvent() Missing Mass techniche ==> Mean="<<Mean/GeV<<"GeV  "
    <<"Sigma= "<<Sigma/GeV<<"GeV,  after smearing P_pr_tot= "<<P_pr_tot/GeV<<"GeV  "<<endl;
    }
    #endif
    */

    //h3P_rtpc->GetRandom3(x,y,z);
    //V3P_rtpc.SetXYZ(x,y,z);

#ifdef HRSPRIMARYROOTHISTO_DEBUG
    if(Global_Debug_Level>=2)
    {
        cout<<"HRSPrimaryRootHisto::GetExclusiveEvent()  BeamE="<<E<<"  TargetMass="<<M<<endl
            <<"  Vertex="<<Vertex<<"  V3P_rtpc="<<V3P_rtpc<<"  V3P_el="<<V3P_el<<endl
            <<"  V3P_pi="<<V3P_pi<<"  V3P_clas="<<V3P_clas<<endl;
    }
#endif
    return;
}

void HRSPrimaryRootHisto::GetYbyX(TH2 *h2, double x, double &y)
{
    int bin;
    //get  X distribution
    TH1D *xproj = (TH1D*) h2->ProjectionX("xproj"); 
    bin = xproj->FindBin(x);
    TH1D *yproj = (TH1D*) h2->ProjectionY("yproj",bin-2,bin+2); 
    y=yproj->GetRandom(); 
    delete xproj;
    delete yproj;

    //cout << "x=" << x << " ==>  y="<< y <<endl;
    return;
}



void HRSPrimaryRootHisto::GetXbyY(TH2 *h2, double y, double &x)
{
    int bin;
    //get Y distribution
    TH1D *yproj = (TH1D*) h2->ProjectionY("yproj"); 
    bin = yproj->FindBin(y);
    TH1D *xproj = (TH1D*) h2->ProjectionX("xproj",bin-2,bin+2); 
    x=xproj->GetRandom(); 
    delete xproj;
    delete yproj;

    //cout << "y=" << y << " ==>  x="<< x <<endl;
    return;
}

//input: theta in deg
//return mementum in GeV unit
double HRSPrimaryRootHisto::RTPC_P_Threshold(double theta_deg)
{//P0_MeV=7.47E-07x^4-2.67E-04x^3+3.74E-02x^2-2.44x+128.0, where x is theta_deg

    int i;
    double p0=0.;
    double A[]={7.47E-07,-2.67E-04,3.74E-02,-2.44,128.0};
    for(i=4;i>=0;i--)
    {
        p0+=A[i]*pow(theta_deg,4-i);
        //printf("A[%d]=%.8f P0=%f\n",i,A[i],p0);
    }
    return p0/1000.0;
}

//	



