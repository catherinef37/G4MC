// ********************************************************************
//
// $Id: HRSPrimaryGeneratorAction.cc,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
// --------------------------------------------------------------
//
#include "HRSPrimaryGeneratorAction.hh"
#include "HRSPrimaryGeneratorMessenger.hh"
#include "HRSPrimaryRootEvent.hh"
#include "HRSPrimaryRootHisto.hh"
#include "UsageManager.hh"
#include "GlobalDebuger.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "TString.h"

extern UsageManager* gConfig;


//////////////////////////////////////////////////////////////////////////////////////////
//PARIS deuteron wave function, fit pol8 (<160 MeV) and pol3 (>150 MeV)
//already normalized, peak value 0.042
const double kDWave_min=0.000; 
const double kDWave_max=0.042; 
inline double PARIS_DWave(double p_mev)
{
	/****************************************
	Minimizer is Linear
	Chi2                      =  3.91081e-07
	NDf                       =           31
	p0                        =   0.00125991   +/-   0.000174876 
	p1                        =  -0.00074123   +/-   5.05681e-05 
	p2                        =  0.000202088   +/-   4.66432e-06 
	p3                        = -7.63826e-06   +/-   1.95654e-07 
	p4                        =   1.3416e-07   +/-   4.36081e-09 
	p5                        = -1.32008e-09   +/-   5.50711e-11 
	p6                        =  7.50679e-12   +/-   3.95191e-13 
	p7                        = -2.31136e-14   +/-   1.50081e-15 
	p8                        =  2.98865e-17   +/-   2.34062e-18 
	Info in <TCanvas::Print>: png file dwave1.png has been created

	****************************************
	Minimizer is Linear
	Chi2                      =  5.41139e-07
	NDf                       =           59
	p0                        =    0.0396013   +/-   0.000741135 
	p1                        =  -0.00034627   +/-   8.73537e-06 
	p2                        =  1.03511e-06   +/-   3.29518e-08 
	p3                        = -1.03947e-09   +/-   3.99671e-11 
	****************************************/
	double par_8[]={ 0.00125991, -0.00074123, 0.000202088, -7.63826e-06,  
		1.3416e-07, -1.32008e-09, 7.50679e-12,-2.31136e-14,  
		2.98865e-17 };
	double par_3[]={ 0.0396013, -0.00034627, 1.03511e-06, -1.03947e-09};

	double prob=0.0;
	if(p_mev<155.0) 
	{
		for(int i=0;i<=8;i++) prob+=par_8[i]*pow(p_mev,double (i));
	}
	else
	{
		for(int i=0;i<=3;i++) prob+=par_3[i]*pow(p_mev,double (i));
	}
	return prob;
}

//This routine will return the probability to find photon with energy 
//fracjtion y=E_g/Beam 
inline double BremDistr(double y, double *par)
{
	//from PDG (2010)  EQ 27.27
	double E = par[0];
	return (4.0/3.0-4.0*y/3.0+y*y)/(y*E);
}

namespace TRANSFORM
{
	void P_TCS2HCS(double Theta_tr, double Phi_tr, double EndPlaneTheta_hall, 
		double &Theta_hall, double &Phi_hall);
}

HRSPrimaryGeneratorAction::HRSPrimaryGeneratorAction()
{
	gConfig->GetParameter("TargetXOffset",kTargetXOffset);
	kTargetXOffset*=mm;
	gConfig->GetParameter("TargetYOffset",kTargetYOffset);
	kTargetYOffset*=mm;
	gConfig->GetParameter("TargetZOffset",kTargetZOffset);
	kTargetZOffset*=mm;
	gConfig->GetParameter("LSeptumAngle",kLSeptumAngle);
	kLSeptumAngle*=deg;
	gConfig->GetParameter("RSeptumAngle",kRSeptumAngle);
	kRSeptumAngle*=deg;

	//the following variables defined in BField_HelmXXX.ini, in unit of cm
	kHelmCurrentRatio=1.0;
	gConfig->GetParameter("Helm_CurrentRatio",kHelmCurrentRatio); 
	gConfig->GetParameter("Helm_OriginX",kHelmXOffset);
	kHelmXOffset*=cm;
	gConfig->GetParameter("Helm_OriginY",kHelmYOffset);
	kHelmYOffset*=cm;
	gConfig->GetParameter("Helm_OriginZ",kHelmZOffset);
	kHelmZOffset*=cm;

	HMSMomentum=0.0;
	int pSetupHMS=0;
	gConfig->GetParameter("SetupHMS",pSetupHMS);
	if(pSetupHMS) gConfig->GetParameter("HMSMomentum",HMSMomentum); //GeV

	gConfig->GetArgument("BeamEnergy",beamEnergy);
	beamEnergy*=GeV;
	gConfig->GetArgument("TargetMass",targetMass);
	targetMass*=GeV;
	gConfig->GetArgument("LHRSMomentum",leftHRSMomentum); //GeV
	leftHRSMomentum*=GeV;
	gConfig->GetArgument("RHRSMomentum",rightHRSMomentum); //GeV
	rightHRSMomentum*=GeV;
	gConfig->GetArgument("BeamTiltedAngle",beamTiltedAngle);
	beamTiltedAngle*=deg;

	gConfig->GetArgument("TargetAtomicNumber",targetAtomicNumber); 

	double tmpPoint[4]={kTargetXOffset,kTargetYOffset,kTargetZOffset,0};
	const G4Field *theField = G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField();
	theField->GetFieldValue(tmpPoint,kEMFieldAtTg);

	rasterMode=1;
	//the following has been defined for G2P
	//In scatter chamber coordinate system
	kRasterR=10.0*mm;		//raster size
	kZLow=-15.0*mm+kTargetZOffset;		//target vertex lower limit
	kZHigh=15.0*mm+kTargetZOffset;		//target vertex higher limit
	kZLowTrig=-100.0*cm+kTargetZOffset;	//the trigger to use random vertex 
	kZHighTrig=400.0*cm+kTargetZOffset;	//the trigger to use random vertex 

	///////////////////////////////////////////////////////////
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String thisParticleName;
	electron = particleTable->FindParticle(thisParticleName="e-");
	positron = particleTable->FindParticle(thisParticleName="e+");
	muonplus = particleTable->FindParticle(thisParticleName="mu+");
	pion0 = particleTable->FindParticle(thisParticleName="pi0");
	pionplus = particleTable->FindParticle(thisParticleName="pi+");
	pionminus = particleTable->FindParticle(thisParticleName="pi-");
	kaon0 = particleTable->FindParticle(thisParticleName="kaon0");
	kaonplus = particleTable->FindParticle(thisParticleName="kaon-");
	kaonminus = particleTable->FindParticle(thisParticleName="kaon+");
	proton = particleTable->FindParticle(thisParticleName="proton");

	//set default values
	position3V.set(kTargetXOffset,kTargetYOffset,kTargetZOffset);

	gunZLow=kZLow;
	gunZHigh=kZHigh;
	gunZ=kZLowTrig-10.0*cm;

	gunRLow=0.0*mm;
	gunRHigh=kRasterR;

	particleNum=1; //number of primary particles in one event

	//load the engine from cmd line
	char tmpStr[100];
	for(int ii=0;ii<MaxPrimaryNum;ii++)
	{
		sprintf(tmpStr,"PrimaryEngine%d",ii+1);
		gConfig->GetArgument(tmpStr,primaryEngine[ii]);

		//set the inital value for these 2 arrays, they will be used to throw theta_tr and phi_tr
		//they will be updated by G4 command during run time
		randomizeInTCS[ii]=0;
		detectorAngle[ii]=(ii%2)?kRSeptumAngle:kLSeptumAngle;
	}
	//primaryEngine[0]="HRSElasEl"; this is default

	//defalut value for primaries, 
	//random energy, random theta, random phi,random vertex, raster on
	for(int i=0;i<MaxPrimaryNum;i++)
	{
		momentumLow[i]=0.001*MeV;
		momentumHigh[i]=beamEnergy*0.99999;

		momentum[i] = -0.2*GeV;
		sigmaMomentum[i] = 0.*MeV;
		sigmaTheta[i] = 0.*deg;
		sigmaPhi[i] = 0.*deg;

		outPlaneAngleLow[i]=-90.0*deg;
		outPlaneAngleHigh[i]=90.0*deg;
		inPlaneAngleLow[i]=-10.0*deg;
		inPlaneAngleHigh[i]=10.0*deg;

		thetaLow[i]=0.5*deg;
		thetaHigh[i]=120.*deg;
		thetaAngle[i]=-10.*deg;

		phiLow[i]=0.0*deg;
		phiHigh[i]=360.*deg;
		phiAngle[i]=-365.*deg;

		particle[i]=electron;
		momentum3V[i].set(0.,0.,0.);
		coupleToPrimary[i]=0;
		polarization[i]=0;

		useMom3V[i]=false;
	}

	fixedPointBL3V.set(0.,0.,kTargetZOffset);
	slopeBL3V.set(0.,0.,kTargetZOffset);
	vertexFlag=1;

	//random event generator based on real data
	bUseRootHisto=false;
	mRootHisto=new HRSPrimaryRootHisto();

	//ntuple event generator
	bUseRootEvent=false;
	iRootEventType=4;
	mRootEvent=new HRSPrimaryRootEvent();

	//by default, randomize is set to false
	randomizePrimary = false;
	bUseFastProtonGenerator=false;

	G4int n_particle = 1;  //this is the number of identical particle you want to shoot 
	//at a time, it may not be the number of primary you want
	particleGun  = new HRSParticleGun(n_particle);

	//create a messenger for this class
	gunMessenger = new HRSPrimaryGeneratorMessenger(this);

	particleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
	particleGun->SetParticleDefinition(electron);

	photonEFracMax=0.999;
	photonEFracMin=0.01;
	incidentEnergy=beamEnergy;
	nThrown=helicity=0;

	G4cout<<"HRSPrimaryGeneratorAction() construction done!"<<G4endl;
}

HRSPrimaryGeneratorAction::~HRSPrimaryGeneratorAction()
{
	delete particleGun;
	delete gunMessenger;
	G4cout<<"delete HRSPrimaryGeneratorAction ... done!"<<G4endl;
}

	
//input: y = E_gamma / E_el
double HRSPrimaryGeneratorAction::GetBremPhotonPol(double y)
{
	double Pol_e = 0.9;
	return Pol_e * (4*y-y*y)/(4-4*y+3*y*y);
}


void HRSPrimaryGeneratorAction::GetPosition()
{
	//calaulate the vertex position
	//vertexFlag;  //3:fixed point 3V; 2: Random z along BeamLine, with raster 1: Ramdom RZ
	//by default vertexFlag=1,RLow=0,RHigh=kRasterR=(10.0mm); beam line is along z axis
	//case 1:by default R is zero, beam line is along z axis ==>cylindrical shape
	//case 2:by default the pivot point=(0.,0.,kTargetZOffset); slopeBL3V=(0.,0.,kTargetZOffset);
	//       a linear line go throught one fixed pivot point: X=x0+a(z-z0) Y=y0+b(z-z0)
	//       You can set the beam line to any direction and go throught any pivot point, 
	//       if z is out of range, a random z will be used. otherwise it is a fixed point
	//case 3:use the fixed vertex point, no raster will be used

	if(vertexFlag==3) return;

	G4double tmprr,tmpphiphi,tmpxx=0.,tmpyy=0.,tmpzz=0.;
	tmprr=tmpphiphi=0.0;

	///////////////////////////////////////////////////////
	//use random x,y position when raster is on
	if(rasterMode==3)  //raster mode: 1 is circle, 2 is rectangle, 3 is elipse
	{
		do
		{
			tmpxx = mRand.fRand();
			tmpyy = mRand.fRand();
		}while(sqrt(tmpxx*tmpxx+tmpyy*tmpyy)>1.0);

		tmpxx*=gunRLow;
		tmpyy*=gunRHigh;
		if(mRand.fRand()>0.5) tmpxx = -1.0*tmpxx;
		if(mRand.fRand()>0.5) tmpyy = -1.0*tmpyy;

	}
	else if(rasterMode==2)  // 2 is rectangle, gunRLow is half X, gunRHigh is half Y
	{
		tmpxx = mRand.fRand()*gunRLow;
		tmpyy = mRand.fRand()*gunRHigh; 

		if(mRand.fRand()>0.5) tmpxx = -1.0*tmpxx;
		if(mRand.fRand()>0.5) tmpyy = -1.0*tmpyy;
	}
	else
	{
		//This model gives uniform distribution
		//it becomes slow if (gunRHigh-gunRLow)/gunRLow is small 
		int tmpmax=0, tmpfound=1;
		do
		{
			tmpxx = mRand.fRand()*gunRHigh;
			tmpyy = mRand.fRand()*gunRHigh;
			tmprr = sqrt(tmpxx*tmpxx+tmpyy*tmpyy);
			tmpmax++;
			if(tmpmax>100) {tmpfound=0;break;}
		}while(tmprr<gunRLow || tmprr>gunRHigh);

		if(tmpfound)
		{
			if(mRand.fRand()>0.5) tmpxx = -1.0*tmpxx;
			if(mRand.fRand()>0.5) tmpyy = -1.0*tmpyy;
		}
		else
		{
			//By Jixie@20110915: This model did not provide a uniform distribution
			//it is only uniform in R, but runs very fast
			tmprr = mRand.fRand(gunRLow,gunRHigh);
			tmpphiphi = 2.0*3.1415926535*mRand.fRand();
			tmpxx = tmprr*cos(tmpphiphi);
			tmpyy = tmprr*sin(tmpphiphi);
		}
	}

	///////////////////////////////////////////////////////

	if(vertexFlag==2)
	{
		tmpzz = slopeBL3V.z()*mm;
		if (tmpzz<kZLowTrig || tmpzz>kZHighTrig)
		{//use random z position
			//tmpzz = G4UniformRand()*(gunZHigh-gunZLow)+gunZLow;
			tmpzz = mRand.fRand(gunZLow,gunZHigh);
		}

		double tmpxx0=fixedPointBL3V.x()+slopeBL3V.x()*(tmpzz-fixedPointBL3V.z());
		double tmpyy0=fixedPointBL3V.y()+slopeBL3V.y()*(tmpzz-fixedPointBL3V.z());

		//apply the raster and offset
		tmpxx+=tmpxx0+kTargetXOffset;
		tmpyy+=tmpyy0+kTargetYOffset;
		//G4cout<<"debug: Use slopeBL3V="<<slopeBL3V<<"mm"
		//	<<" ==> vertex=("<<tmpxx<<", "<<tmpyy<<", "<<tmpyy<<")" <<G4endl;
	}
	else  if(vertexFlag==1)
	{
		tmpzz = gunZ;
		if (tmpzz<kZLowTrig || tmpzz>kZHighTrig)
		{//use random z position
			//tmpzz = G4UniformRand()*(gunZHigh-gunZLow)+gunZLow;
			tmpzz = mRand.fRand(gunZLow,gunZHigh);
		}

		//apply the raster and offset
		tmpxx+=kTargetXOffset;
		tmpyy+=kTargetYOffset;
	}
	position3V.set(tmpxx,tmpyy,tmpzz);

}

void HRSPrimaryGeneratorAction::GetMomentum(int i)
{
	//get the momentum,
	//if useMom3V[i] is true, use 3-vector momentum3V[i] directorly
	//otherwise get the momentum from total momentum,theta and phi commands
	//Be sure the folowing randomize trigger:
	//  in case momentum[i]<0, generate total P randomly
	//  in case thetaAngle[i]<0, generate theta randomly
	//  in case phiAngle[i]<-360, generate phi randomly

	G4double pPtot=momentum3V[i].mag();
	if(useMom3V[i] && pPtot>=keV) return;

	//comparing string is very slow ...... need to improve this
	if(primaryEngine[i]=="Uniform")			 {/*do thing, jump over this if block;*/;}
	else if(primaryEngine[i]=="HRSElasEl")           {HRSElasElectronEngine(i);return;}
	else if(primaryEngine[i]=="HRSElasNucleus") {HRSElasNucleusEngine(i); return;}
	else if(primaryEngine[i]=="HRSQuasiElasEl") {HRSQuasiElasElectronEngine(i); return;}
	else if(primaryEngine[i]=="HRSQuasiElasNucleon") {HRSQuasiElasNucleonEngine(i); return;}
	else if(primaryEngine[i]=="Compton")     {ComptonEngine(i); return;}
	else if(primaryEngine[i]=="TwoBody")     {TwoBodyEngine(i); return;}
	else if(primaryEngine[i]=="BonusProton") {BoNuSProtonEngine(i); return;}
	else if(primaryEngine[i]=="FastProton")  {FastProtonEngine(i); return;}
	else if(primaryEngine[i]=="RootHisto")   {HRSHistoEngine(i); return;}
	else if(primaryEngine[i]=="RootNtuple")  {HRSNtupleEngine(i); return;}
	else if(primaryEngine[i]=="H90UserFit")  {H90UserFitEngine(i); return;}
	else if(primaryEngine[i]=="BdL")         {BdLEngine(i); return;}
	else 
	{
		G4cout<<"Warning: Unknown event generator \""<<primaryEngine[i]
		<<"\". Assuming 'Uniform' !"<<endl; 
	}

	//Uniform Engine in HCS or TCS
	G4double pTheta=0.0,pPhi=0.0;
	if(randomizeInTCS[i]>0) RandTCSThetaPhi(i,pTheta,pPhi);
	else RandHCSThetaPhi(i,pTheta,pPhi);

	if (momentum[i]<=0.*GeV)
	{//use random total momentum; Let low<=total momentum<high
		pPtot = mRand.fRand(momentumLow[i],momentumHigh[i]);
	}
	else
	{//use specify total momentum and smear it with sigma
		pPtot = mRand.fGaus(momentum[i],sigmaMomentum[i]);
	}

	//G4cout<<"Debug GetMomentum(): pTheta="<<pTheta/deg<<"deg  pPhi="<<pPhi/deg<<"deg \n";
	momentum3V[i].setRThetaPhi(pPtot,pTheta,pPhi);
}


//  Note that, even if the particle gun shoots more than one particles at one invokation of
// GeneratePrimaryVertex() method, all particles have the same physical quantities. If the
// user wants to shoot two particles with different momentum, position, etc., invoke
// GeneratePrimaryVertex() method twice and set quantities on demand to the particle gun.
//
void HRSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//I have to reset the incidentEnergy for each event. incidentEnergy is mainly for 
	//comptom engine, or one can also use it to store the radiated incoming energy of an electron
	incidentEnergy=beamEnergy;
	helicity=(nThrown++)%4;

	if(bUseRootEvent)
	{
		this->GenPrimariesFromHisto(anEvent);
		return;
	}

	//get vertex
	this->GetPosition();  //the result will be stored at postion3V

	for(int i=0;i<particleNum;i++)
	{
		effectiveTheta[i]=-10.0;
		if(randomizePrimary)
		{
			G4int ii = mRand.iRand(0, 11);
			switch(ii)
			{
			case 0:
				particle[i] = positron;
				break;
			case 1:
				particle[i] = electron;
				break;
			case 2:
				particle[i] = muonplus;
				break;
			case 3:
				particle[i] = pion0;
				break;
			case 4:
				particle[i] = pionplus;
				break;
			case 5:
				particle[i] = pionminus;
				break;
			case 6:
				particle[i] = kaon0;
				break;
			case 7:
				particle[i] = kaonplus;
				break;
			case 8:
				particle[i] = kaonminus;
				break;
			case 9:
				particle[i] = proton;
				break;
			default:
				particle[i] = proton;
				break;
			}
		}
		else
		{
			//particle[i] have been set by user cmd, if no user cmd is invoked, 
			//use the default paricle electron

			//this will use the last particle or particle specified in gun command 
			//like "/gun/particle alpha"
			if(!particle[i]) particle[i] = particleGun->GetParticleDefinition();
		}

		//use the flat spectator particle generator for the first particle 
		if(i==0 && bUseFastProtonGenerator) FastProtonEngine(i);
		else GetMomentum(i);  //the result will be stored at momentum3V[i]

		//particleGun->SetParticlePolarization(G4ThreeVector)
		particleGun->SetParticlePosition(position3V);
		particleGun->SetParticleDefinition(particle[i]);

		particleGun->SetParticleMomentum(momentum3V[i]);

		/*
		//Or use the following 2 methods together to set the momentum
		//inline void SetParticleEnergy(G4double aKineticEnergy)
		//inline void SetParticleMomentumDirection(G4ParticleMomentum aMomentumDirection)
		G4double mass =  particle[i]->GetPDGMass();
		G4double ekin=std::sqrt(momentum3V[i].mag2()+mass*mass)-mass;
		particleGun->SetParticleEnergy(ekin);
		particleGun->SetParticleMomentumDirection(momentum3V[i].unit());
		*/

		particleGun->GeneratePrimaryVertex(anEvent);
	}
}


//Fast spectator proton generator
//if user has specified the momentum using 3V command, e.g.
//  $>/mydet/particle1/momentum3V 10.7 50.0 20.9 MeV 
//then do nothing. otherwise, generate flat theta, phi acorrding to z, then
//generate flat momenta between momentumLow[0],momentumHigh[0], which are
//specified by user.
void HRSPrimaryGeneratorAction::FastProtonEngine(int index)
{
	G4double dPtotal=momentum3V[index].mag();
	if(useMom3V[index])  return;

	G4double dTheta,dPhi;
	double z0_mm=position3V.z()/mm;
	if (thetaAngle[index]/deg<0.)
	{//use random theta angle//Let low<=theta<high
		//the inner R and outer R of drift region is 30.mm and 60.mm,
		//I take the average 45. to be the exit window low boundary at the end plate
		const double PI=4.0*atan(1.0);
		G4double tmpThetaHigh=PI-atan(45.0/(z0_mm+105.));
		G4double tmpThetaLow=atan(45.0/(-z0_mm+105.0));
		dTheta = mRand.fRand(tmpThetaLow,tmpThetaHigh);
	}
	else
	{//use specify theta angle
		dTheta=mRand.fGaus(thetaAngle[index],sigmaTheta[index]);
	}

	//phi distribution should be flat
	if (phiAngle[index]/deg<-360.)
	{//use random phi angle
		dPhi = mRand.fRand(0.,360.) *deg;
	}
	else
	{//use specify phi angle
		dPhi = mRand.fGaus(phiAngle[index],sigmaPhi[index]);
	}

	//use random total momentum; Let low<=total momentum<high
	//dPtotal = mRand.fRand(momentumLow[index],momentumHigh[index]);

	//throw Ps with dwave distribution between (0,300) MeV
	dPtotal = mRand.VNReject(momentumLow[index]/MeV,momentumHigh[index]/MeV,
		kDWave_min,kDWave_max,PARIS_DWave)*MeV;

	momentum3V[index].setRThetaPhi(dPtotal,dTheta,dPhi);
}

//same as fase proton but require it passes the BoNuS threshold
void HRSPrimaryGeneratorAction::BoNuSProtonEngine(int index)
{
	G4double dPtotal=momentum3V[index].mag();
	if(useMom3V[index] && dPtotal>=65.*MeV)  return;

	G4double dTheta,dPhi;
	double z0_mm=position3V.z()/mm;
	if (thetaAngle[index]/deg<0.)
	{//use random theta angle//Let low<=theta<high
		//the inner R and outer R of drift region is 30.mm and 60.mm,
		//I take the average 45. to be the exit window low boundary at the end plate
		const double PI=3.1415926535*rad;
		G4double tmpThetaHigh=PI-atan(45.0/(z0_mm+105.));
		G4double tmpThetaLow=atan(45.0/(-z0_mm+105.0));
		dTheta = mRand.fRand(tmpThetaLow,tmpThetaHigh);
	}
	else
	{//use specify theta angle
		dTheta=mRand.fGaus(thetaAngle[index],sigmaTheta[index]);
	}

	//phi distribution should be flat
	if (phiAngle[index]/deg<-360.)
	{//use random phi angle
		dPhi = mRand.fRand(0.,360.) *deg;
	}
	else
	{//use specify phi angle
		dPhi = mRand.fGaus(phiAngle[index],sigmaPhi[index]);
	}

	//use random total momentum; Let Low(or Threshold)<=momentum<high
	double pBoNuSThreshold=0.,xx=dTheta/deg;
	double A[]={7.47E-07,-2.67E-04,3.74E-02,-2.44,128.0};
	for(int ii=4;ii>=0;ii--) pBoNuSThreshold+=A[ii]*pow(xx,4-ii); 
	pBoNuSThreshold*=MeV;

	double pMin=(pBoNuSThreshold<momentumLow[index])?pBoNuSThreshold:momentumLow[index];
	dPtotal = mRand.fRand(pMin,momentumHigh[index]);

	momentum3V[index].setRThetaPhi(dPtotal,dTheta,dPhi);
}

void HRSPrimaryGeneratorAction::RandTCSThetaPhi(int i,double &theTheta,double &thePhi)
{
	double pTheta0_tr = mRand.fRand(outPlaneAngleLow[i],outPlaneAngleHigh[i]);
	double pPhi0_tr = mRand.fRand(inPlaneAngleLow[i],inPlaneAngleHigh[i]);

	//do elipse shape if randomizeInTCS==2
	if(randomizeInTCS[i]==2)
	{
		double X0=(outPlaneAngleLow[i]+outPlaneAngleHigh[i])/2.0;
		double Y0=(inPlaneAngleLow[i]+inPlaneAngleHigh[i])/2.0;
		double A=(outPlaneAngleHigh[i]-outPlaneAngleLow[i])/2.0;
		double B=(inPlaneAngleHigh[i]-inPlaneAngleLow[i])/2.0;
		double R2=pow((pTheta0_tr-X0)/A,2.0) + pow((pPhi0_tr-Y0)/B,2.0);
		while(R2>1.0)
		{
			pTheta0_tr = mRand.fRand(outPlaneAngleLow[i],outPlaneAngleHigh[i]);
			pPhi0_tr = mRand.fRand(inPlaneAngleLow[i],inPlaneAngleHigh[i]);

			X0=(outPlaneAngleLow[i]+outPlaneAngleHigh[i])/2.0;
			Y0=(inPlaneAngleLow[i]+inPlaneAngleHigh[i])/2.0;
			A=(outPlaneAngleHigh[i]-outPlaneAngleLow[i])/2.0;
			B=(inPlaneAngleHigh[i]-inPlaneAngleLow[i])/2.0;
			R2=pow((pTheta0_tr-X0)/A,2.0) + pow((pPhi0_tr-Y0)/B,2.0);
		}
	}

	double pTheta=detectorAngle[i],pPhi=0;

	//calculate these angle directly, which will be faster
	pTheta=acos(cos(detectorAngle[i]+pPhi0_tr)*cos(pTheta0_tr));
	pPhi=atan( -tan(pTheta0_tr) / sin(detectorAngle[i]+pPhi0_tr) );
	//only good for left arm, need to +180 for right arm, why?
	if(sin(detectorAngle[i])<0) pPhi=180.*deg+pPhi;
	//or call the routines
	//TRANSFORM::P_TCS2HCS(pTheta0_tr,pPhi0_tr,detectorAngle[i],pTheta,pPhi);

	//now set the angle only if the random trigger is on
	if (thetaAngle[i]/deg<0.) theTheta = pTheta;
	else theTheta= mRand.fGaus(thetaAngle[i],sigmaTheta[i]);

	if (phiAngle[i]/deg<-360.) thePhi = pPhi;
	else thePhi = mRand.fGaus(phiAngle[i],sigmaPhi[i]);

	//G4cout<<"Theta0_tr ="<<pTheta0_tr/mrad<<"mrad  Phi0_tr ="<<pPhi0_tr/mrad
	//	<<"Theta0_lab="<<pTheta/mrad<<"mrad  Phi0_lab="<<pPhi/mrad<<G4endl;
}

void HRSPrimaryGeneratorAction::RandHCSThetaPhi(int i,double &pTheta,double &pPhi)
{
	if (thetaAngle[i]/deg<0.)
	{
		//use random theta angle,Let low<=theta<high
		//By Jixie:  we should randomize in costheta, not theta
		pTheta = acos(mRand.fRand(cos(thetaLow[i]),cos(thetaHigh[i])));
	}
	else
	{
		//use specify theta angle, gaussian distribution
		pTheta = mRand.fGaus(thetaAngle[i],sigmaTheta[i]);
	}

	if (phiAngle[i]/deg<-360.)
	{
		//use random phi angle, have the same unit as phiHigh phiLow
		pPhi = mRand.fRand(phiLow[i],phiHigh[i]);
	}
	else
	{
		//use specify phi angle, gaussian distribution
		pPhi = mRand.fGaus(phiAngle[i],sigmaPhi[i]);
	}

}

void HRSPrimaryGeneratorAction::H90UserFitEngine(int index)
{
	double pPtot=0.0;
	if (momentum[index]<=0.*GeV)
	{
		//use random total momentum; Let low<=total momentum<high
		pPtot = mRand.fRand(momentumLow[index],momentumHigh[index]);
	}
	else
	{
		//use specify total momentum and smear it with sigma
		pPtot = mRand.fGaus(momentum[index],sigmaMomentum[index]);
	}

	double x=pPtot/GeV,p[4];
	double pTheta,pPhi;
	if(randomizeInTCS[index]>0) 
	{
		double pMeanTheta0_tr=0.0,pMeanPhi0_tr=0.0;
		/*  		
		for kTargetYOffset = 0 only
		6 degrees, fit to Theta0_tr vs P0 f(x) = p0/x+p1+p2*x  
		NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
		1  p0         -2.23271e-001  1.46173e-003  2.85477e-006 -1.49338e-006
		2  p1         -2.78692e-003  2.52935e-003  2.38083e-006  2.93668e-005
		3  p2          9.41369e-004  8.30292e-004  1.14022e-006  9.40355e-005

		12.5 degrees, fit to Theta0_tr vs P0 f(x) = p0/x+p1+p2*x  
		NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
		1  p0         -1.88639e-001  1.40148e-003  3.63280e-006 -4.49860e-006
		2  p1         -2.80966e-002  2.43511e-003  3.06206e-006  2.30655e-005
		3  p2          8.16754e-003  8.02532e-004  1.47375e-006  7.81053e-005
		*/
		if(fabs(cos(detectorAngle[index])-cos(6.0*deg))<0.005)
		{
			//6 degree Theta0_tr Average, 
			//parameters fitted in unit of degrees
			if(fabs(kTargetYOffset-0.0*cm)<1.5*cm)
			{
				p[0] = -9.44; p[1] = -3.041; p[2] = 0.7901; 
			}
			if(fabs(kTargetYOffset-3.0*cm)<1.5*cm)
			{
				p[0] = -10.16; p[1] = -0.8209; p[2] = 0.5809; 
			}
			if(fabs(kTargetYOffset-6.0*cm)<1.5*cm)
			{
				p[0] = -10.82; p[1] = 1.331; p[2] = 0.3738; 
			}
			if(fabs(kTargetYOffset-9.0*cm)<1.5*cm)
			{
				p[0] = -9.653; p[1] = -1.536; p[2] = 0.625; 
			}
			else
			{
				p[0] = -9.44; p[1] = -3.041; p[2] = 0.7901; 
			}
			pMeanTheta0_tr=(p[0]/x+p[1]+p[2]*x)*deg;
		}
		else if(fabs(cos(detectorAngle[index])-cos(12.5*deg))<0.005)
		{
			//parameters fitted in unit of rad
			p[0] = -1.88639e-001; p[1] = -2.80966e-002; p[2] =  8.16754e-003; 
			pMeanTheta0_tr=p[0]/x+p[1]+p[2]*x;
		}
		else 
		{
			//parameters fitted in unit of rad
			p[0] = -1.88639e-001; p[1] = -2.80966e-002; p[2] =  8.16754e-003; 
			pMeanTheta0_tr=p[0]/x+p[1]+p[2]*x;
		}
		pMeanTheta0_tr*=kHelmCurrentRatio;

		//the range is +/- 80% of the center value
		double pTheta0_trRange=0.08;
		//if(x<0.6) pTheta0_trRange=0.20;
		//else if(x<1.0) pTheta0_trRange=0.15;
		//else if(x<2.0) pTheta0_trRange=0.12;
		//else if(x<3.5) pTheta0_trRange=0.10;
		//else pTheta0_trRange=0.08;

		outPlaneAngleLow[index]=pMeanTheta0_tr-pTheta0_trRange;
		outPlaneAngleHigh[index]=pMeanTheta0_tr+pTheta0_trRange;

		inPlaneAngleLow[index]=pMeanPhi0_tr-0.05;
		inPlaneAngleHigh[index]=pMeanPhi0_tr+0.05;
		RandTCSThetaPhi(index,pTheta,pPhi);
	}
	else 
	{
		//randomize in HCS
		/*  
		6 degrees, fit to Theta0 vs P0  f(x) = p0/x+p1+p2*x
		NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
		1  p0          2.18866e-001  1.32490e-003  2.64764e-006 -1.55922e-005
		2  p1          1.98425e-002  2.15861e-003  1.91189e-006 -1.26731e-005
		3  p2          1.29083e-002  6.74513e-004  8.44859e-007 -1.18584e-005

		12.5 degrees, fit to Theta0 vs P0  f(x) = p0/x+p1+p2*x
		NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
		1  p0          1.60763e-001  1.08128e-003  2.46757e-006  9.67520e-006
		2  p1          1.35127e-001  1.72363e-003  1.70439e-006  2.16783e-006
		3  p2          1.44681e-002  5.37374e-004  7.45904e-007 -1.00022e-005
		*/

		double pMeanTheta0=0.0,pMeanPhi0=0.0;
		if(fabs(cos(detectorAngle[index])-cos(6.0*deg))<0.005)
		{
			//6 degree Theta0 Average
			//parameters fitted in unit of degrees
			if(fabs(kTargetYOffset-0.0*cm)<1.5*cm)
			{
				p[0] = 8.869; p[1] = 4.381; p[2] = -0.05345; 
			}
			if(fabs(kTargetYOffset-3.0*cm)<1.5*cm)
			{
				p[0] = 9.556; p[1] = 2.074; p[2] = 0.4412; 
			}
			if(fabs(kTargetYOffset-6.0*cm)<1.5*cm)
			{
				p[0] = 10.13; p[1] = -0.0769; p[2] = 0.9697; 
			}
			if(fabs(kTargetYOffset-9.0*cm)<1.5*cm)
			{
				p[0] = 8.665; p[1] = -0.1527; p[2] = 1.11; 
			}
			else
			{
				p[0] = 8.869; p[1] = 4.381; p[2] = -0.05345; 
			}
			pMeanTheta0=(p[0]/x+p[1]+p[2]*x)*deg;
		}
		else if(fabs(cos(detectorAngle[index])-cos(12.5*deg))<0.005)
		{
			p[0] = 1.60763e-001;
			p[1] = 1.35127e-001;
			p[2] = 1.44681e-002; 
			p[3] = 0;
			pMeanTheta0=p[0]/x+p[1]+p[2]*x;
		}
		else 
		{
			p[0] = 1.60763e-001;
			p[1] = 1.35127e-001;
			p[2] = 1.44681e-002; 
			p[3] = 0;
			pMeanTheta0=p[0]/x+p[1]+p[2]*x;
		}

		//the range is +/- 80% of the center value
		double pTheta0Range=0.8*pMeanTheta0;
		//if(x<0.6) pTheta0Range=0.20;
		//else if(x<1.0) pTheta0Range=0.15;
		//else if(x<2.0) pTheta0Range=0.12;
		//else if(x<3.5) pTheta0Range=0.10;
		//else pTheta0Range=0.08;
		if(pTheta0Range<0.08) pTheta0Range=0.08;

		pTheta=mRand.fGaus(pMeanTheta0,pTheta0Range);

		/*
		6 degrees, fit to Phi0 vs P0  f(x)=p0/x+p1+p2*x
		NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
		1  p0         -2.25265e-003  3.10550e-003  2.84987e-006 -3.07966e-004
		2  p1          1.55681e+000  1.16818e-002  5.61828e-006 -3.64233e-006
		3  p2         -4.90198e-001  1.13504e-002  5.38191e-006  2.57751e-004
		4  p3          5.74754e-002  2.78736e-003  2.45831e-006  7.60456e-004		
		12.5 degrees, fit to Phi0 vs P0  f(x)=p0/x+p1+p2*x+p3*x*x  
		NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
		1  p0          8.27403e-002  3.87983e-003  5.70135e-006  2.33003e-004
		2  p1          1.13122e+000  1.25472e-002  8.42188e-006  9.67033e-005
		3  p2         -5.31028e-001  1.05687e-002  5.39756e-006 -1.74557e-005
		4  p3          8.06495e-002  2.36553e-003  2.11588e-006 -2.33861e-004
		*/
		if(fabs(cos(detectorAngle[index])-cos(6.0*deg))<0.005)
		{
			//6 degree Phi0 Average
			//parameters fitted in unit of degrees
			if(fabs(kTargetYOffset-0.0*cm)<1.5*cm)
			{
				p[0] = -4.908; p[1] = 98.74; p[2] = -34.83; p[3] = 4.436; 
			}
			if(fabs(kTargetYOffset-3.0*cm)<1.5*cm)
			{
				p[0] = -5.201; p[1] = 103.0; p[2] = -42.83; p[3] = 5.51; 
			}
			if(fabs(kTargetYOffset-6.0*cm)<1.5*cm)
			{
				p[0] = -7.808; p[1] = 115.2; p[2] = -57.99; p[3] = 8.061; 
			}
			if(fabs(kTargetYOffset-9.0*cm)<1.5*cm)
			{
				p[0] = -4.333; p[1] = 114.1; p[2] = -67.1; p[3] = 9.769;
			}
			else
			{
				p[0] = -4.908; p[1] = 98.74; p[2] = -34.83; p[3] = 4.436; 
			}
			pMeanPhi0 = (p[0]/x+p[1]+p[2]*x+p[3]*x*x)*deg;
		}
		else if(fabs(cos(detectorAngle[index])-cos(12.5*deg))<0.005)
		{
			p[0] = 8.27403e-002;
			p[1] = 1.13122e+000;
			p[2] = -5.31028e-001; 
			p[3] = 8.06495e-002;
			pMeanPhi0 = p[0]/x+p[1]+p[2]*x+p[3]*x*x;	
		}
		else
		{
			//for Yoffset=0 6deg
			p[0] = -2.25265e-003;
			p[1] = 1.55681e+000;
			p[2] = -4.90198e-001; 
			p[3] = 5.74754e-002;
			pMeanPhi0 = p[0]/x+p[1]+p[2]*x+p[3]*x*x;
		}

		//judge if it is the right arm, for right arm, phi=180-phi
		if(sin(detectorAngle[index])<0) pMeanPhi0=180.0*deg-pMeanPhi0;
		pPhi=mRand.fGaus(pMeanPhi0,0.3);
	}

	//check if Theta out of range
	if(pTheta<0.*deg) {pTheta*=-1;}
	if(pTheta>180.*deg) {pTheta=fmod(pTheta,180.0*deg);}
	momentum3V[index].setRThetaPhi(pPtot,pTheta,pPhi);
}


void HRSPrimaryGeneratorAction::ComptonEngine(int index)
{	
	int couple2ThisParticle = (coupleToPrimary[index]>0)?coupleToPrimary[index]-1:0;
	if(index>0 && couple2ThisParticle<index)
	{
		//couple with previous electron
		//Beam + rest_N = e' + N'
		//N' = (Beam + rest_N) - e' = (-P_x_e',-P_y_e',E-P_z_e',E+M_N-E')
		momentum3V[index].setX(-momentum3V[couple2ThisParticle].x());
		momentum3V[index].setY(-momentum3V[couple2ThisParticle].y());
		momentum3V[index].setZ(incidentEnergy-momentum3V[couple2ThisParticle].z());
		
		//The following block of code shows that the above method gives the
		//same result as GetComptonMomentumX
		//debug:
		//G4cout<<"ComptonEngine():: old: momentum3V[index]="<<momentum3V[index]<<G4endl;
		//double pTheta=momentum3V[couple2ThisParticle].theta();
		//double kPi = asin(1.0)*2.0;
		//double pPhi=momentum3V[couple2ThisParticle].phi() + kPi;
		//double pMassP = proton_mass_c2;
		//double pE_gamma, pE_p, pTheta_p, pE0=beamEnergy;
		//	
		//GetComptonMomentumX(pE0,pTheta,pMassP,pE_gamma,pE_p,pTheta_p);
		//momentum3V[index].setRThetaPhi(pE_p,pTheta_p,pPhi);
		//G4cout<<"ComptonEngine():: new: momentum3V[index]="<<momentum3V[index]<<G4endl;
	}
	else
	{
		if(incidentEnergy<=0.0 || incidentEnergy>beamEnergy-0.0001)
		{
			//flat distr
			//incidentEnergy = mRand.fRand(photonEFracMin,photonEFracMax) * beamEnergy;
			//brem XS distr
			double par[]={beamEnergy};
			double kBrem_min = BremDistr(1.0,par);
			double kBrem_max = BremDistr(photonEFracMin,par);
			incidentEnergy = mRand.VNReject(photonEFracMin,photonEFracMax,kBrem_min,kBrem_max,
				par,BremDistr) * beamEnergy;
		}

		//for index=0, assuming this is the gamma
		double pTheta,pPhi;

		if(randomizeInTCS[index]>0) RandTCSThetaPhi(index,pTheta,pPhi);
		else RandHCSThetaPhi(index,pTheta,pPhi);
		//For compton engine, gamma will not be affected by the target field
		//and we assume the beam is not tilted. Therefore no need to calculate 
		//the effective angle
		effectiveTheta[index] = pTheta;

		//RCS energy at theta_gamma_lab angle
		double pM_pr = proton_mass_c2;
		double pPtot = 0;

		//if this is a photon, set its polarization value
		if(particle[index]->GetPDGEncoding()==22) 
		{
			double cosTh_gamma = cos(effectiveTheta[index]);
			pPtot =  incidentEnergy / (1.0 + incidentEnergy/pM_pr * (1.0-cosTh_gamma) );
		}
		else
		{
			//I am assuming this is proton all the time
			pPtot = GetElasMomentumX(incidentEnergy, effectiveTheta[index], pM_pr);
		}
	
		momentum3V[index].setRThetaPhi(pPtot,pTheta,pPhi);
	}

	//if this is a photon, set its polarization value
	if(particle[index]->GetPDGEncoding()==22)
	{
		polarization[index] = GetBremPhotonPol(momentum3V[index].mag()/beamEnergy); 
	}
}


//calculate the BdL (tesla.m) from target to 805mm away along ray with angle of pEndPlaneAngle  
void HRSPrimaryGeneratorAction::GetHelmBdL(G4ThreeVector &V3BdL, double pEndPlaneAngle)
{
	G4FieldManager *theFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
	
	double pStepLength=1.0*mm;
	G4ThreeVector V3UnitCentralRay(sin(pEndPlaneAngle),0,cos(pEndPlaneAngle));
	G4ThreeVector tmpV3B(0,0,0);

	double tmpPos[4]={0,0,0,0};
	double tmpField[6]={0,0,0,0,0,0};

	//acummulate BdL in step length of 1 cm along the central ray up to 805 mm
	double pStep0=(position3V.z()-kHelmZOffset)/cos(pEndPlaneAngle);
	//G4cout<<"Z0="<<position3V.z()/cm<<"cm  HelmZ="<<kHelmZOffset/cm<<"cm"<<G4endl;
	
	V3BdL.set(0,0,0);
	for(double step=pStep0;step<805*mm;step+=pStepLength)
	{
		tmpPos[0]=step*sin(pEndPlaneAngle)+kHelmXOffset;	
		tmpPos[1]=kHelmYOffset;
		tmpPos[2]=step*cos(pEndPlaneAngle)+kHelmZOffset;  
		theFieldManager->GetDetectorField()->GetFieldValue(tmpPos,tmpField);  //in tesla
		tmpV3B.set(tmpField[0]/tesla,tmpField[1]/tesla,tmpField[2]/tesla);
		V3BdL += pStepLength/m*V3UnitCentralRay.cross(tmpV3B);
		//G4cout<<"TrackL="<<step/cm<<"cm  Pos Z="<<tmpPos[2]/cm<<"cm  V3B="<<tmpV3B/tesla<<"tesla  V3BdL="<<V3BdL<<G4endl;
	}
	//G4cout<<"Start from Z0="<<position3V.z()/cm<<"cm  V3BdL="<<V3BdL<<G4endl;
	return;
}


void HRSPrimaryGeneratorAction::BdLEngine(int index)
{
	double pPtot=0.0;
	if (momentum[index]<=0.*GeV)
	{
		//use random total momentum; Let low<=total momentum<high
		pPtot = mRand.fRand(momentumLow[index],momentumHigh[index]);
	}
	else
	{
		//use specify total momentum and smear it with sigma
		pPtot = mRand.fGaus(momentum[index],sigmaMomentum[index]);
	}

	//get the BdL for current detector angle and field ratio
	//can speed up if use fitted function:  pBdL_y = 0.667 - 5.0 * dZ_m * kHelmCurrentRatio; 
	G4ThreeVector pV3BdL(0,0,0);
	GetHelmBdL(pV3BdL,detectorAngle[index]);
	

	double x=pPtot/GeV;
	double pTheta,pPhi;
	//will alway throw in TCS

	double pMeanTheta0_tr=0.0,pMeanPhi0_tr=0.0;
	//dAngle = 0.3 *BdL / P;

	pMeanTheta0_tr = -0.3*pV3BdL.y()/x;
	pMeanPhi0_tr   =  0.3*pV3BdL.x()/x;

	//the range is +/- 80mrad of the center value
	double pTheta0_tr_Range=0.08,pPhi0_tr_Range=0.05;

	outPlaneAngleLow[index]=pMeanTheta0_tr-pTheta0_tr_Range;
	outPlaneAngleHigh[index]=pMeanTheta0_tr+pTheta0_tr_Range;

	inPlaneAngleLow[index]=pMeanPhi0_tr-pPhi0_tr_Range;
	inPlaneAngleHigh[index]=pMeanPhi0_tr+pPhi0_tr_Range;
	RandTCSThetaPhi(index,pTheta,pPhi);

	//check if Theta out of range
	if(pTheta<0.*deg) {pTheta*=-1;}
	if(pTheta>180.*deg) {pTheta=fmod(pTheta,180.0*deg);}
	momentum3V[index].setRThetaPhi(pPtot,pTheta,pPhi);
}


//if the beam is tilted, need to get the the effective scattering angle
//this routine will be used only if the beam is tilted and the particle
//is a charged particle
double HRSPrimaryGeneratorAction::GetEffectiveTheta(double pZ, double pTheta, double pPhi)
{
	if(fabs(kEMFieldAtTg[0])<1.0E-6)  return pTheta; 

	//get the integrated BdL
	//assuming uniform field in X direction, 
	//double pBdL = 5.0*kHelmCurrentRatio*(pZ-kTargetZOffset);
	double pBdL = kEMFieldAtTg[0]*(pZ-kHelmZOffset)/m;  //in tesla.m
	double pDeltaTheta=0.3/beamEnergy*pBdL;  //in unit of rad
	double pThetaBeam=beamTiltedAngle-pDeltaTheta;
	double pPhiBeam=90.0*deg;
	if(pThetaBeam<0.0) {pThetaBeam=fabs(pThetaBeam);pPhiBeam*=-1.0;}

	//declared as module public to avoid contruction and deconstruction,
	//G4ThreeVector mBeamV3,mParticleV3;

	mBeamV3.setRThetaPhi(beamEnergy,pThetaBeam,pPhiBeam);  
	mParticleV3.setRThetaPhi(1.0,pTheta,pPhi); 
	double pThetaEff=mBeamV3.angle(mParticleV3);

	return 	pThetaEff;
}

void HRSPrimaryGeneratorAction::HRSElasElectronEngine(int index)
{
	int couple2ThisParticle = (coupleToPrimary[index]>0)?coupleToPrimary[index]-1:0;
	if(index>0 && couple2ThisParticle<index)
	{
		//couple with previous particle, none-electron
		//Beam + rest_N = e' + N'
		//e' = (Beam + rest_N) - N' = (-P_x_N',-P_y_N',E-P_z_N',E+M_N-E_N')
		momentum3V[index].setX(-momentum3V[couple2ThisParticle].x());
		momentum3V[index].setY(-momentum3V[couple2ThisParticle].y());
		momentum3V[index].setZ(beamEnergy-momentum3V[couple2ThisParticle].z());
	}
	else
	{
		G4double pTheta,pPhi;

		if(randomizeInTCS[index]>0) RandTCSThetaPhi(index,pTheta,pPhi);
		else RandHCSThetaPhi(index,pTheta,pPhi);
		effectiveTheta[index]=GetEffectiveTheta(position3V.z(),pTheta,pPhi);

		//E-N elastic energy at theta angle
		double pEprime=beamEnergy/(1.0+ beamEnergy*(1-cos(effectiveTheta[index]))/targetMass );
		double pMass = electron_mass_c2;
		double pPtot=sqrt(pEprime*pEprime-pMass*pMass);
		momentum3V[index].setRThetaPhi(pPtot,pTheta,pPhi);
	}

}

//This is the routine based on Jixie's calculation
//It gives the elastic momentum for X' in  kX->k'X'elastic scattering
//where k can be either electron or photon
double HRSPrimaryGeneratorAction::GetElasMomentumX(double E0,double Th_x,double M_x)
{
	double A = E0 * cos(Th_x) / (E0 + M_x);
	double P_x = 2 * M_x * A / (1 - A*A);
	return P_x;
}

//It turns out that this routine does not work
//It does gives the same elastic momentum as the above one
double HRSPrimaryGeneratorAction::GetElasMomentumX_old(double pE0,double pTheta,double pMtg)
{
	//pMx is the mass of the outgoing Nucleon
	double pPlab=0.0, pMe=electron_mass_c2,pMx=pMtg;
	double W2 = pMtg*pMtg + 2*pE0*pMtg;
	double W = sqrt(W2);
	double Beta_cm_x = pE0 / (pE0+pMtg);
	double Gamma_cm_x = (pE0+pMtg) / W;
	double E_cm_x = (W2 - pMe*pMe + pMx*pMx) / (2.*W);

	double cosTheta = cos(pTheta);
	double coefA = 1. - pow(Beta_cm_x*cosTheta,2.0);
	double coefB = -2.*E_cm_x*Beta_cm_x*cosTheta/Gamma_cm_x;
	double coefC = pMx*pMx - pow(E_cm_x/Gamma_cm_x,2.0);
	double delta = coefB*coefB - 4.*coefA*coefC;
	if(delta < 0.0)  pPlab = -100.;
	else  pPlab = (sqrt(delta)-coefB) / (2.*coefA);
	return pPlab;
}


void HRSPrimaryGeneratorAction::GetComptonMomentumX(double pE0, double pTheta_gamma, double pMtg,
													double &pE_gamma, double &pP_x, double &pTheta_x)
{
	if(pE0<1.0E-8)  return;
	//E'_gamma = Ei / (1.0 + Ei/Mtg * (1.0-cos(Theta_gamma) )
	//V3P' = ( V3K - V3K')
	//P' = sqrt(K^2 + K'^2 -2*K*K'*cos(Theta_gamma))
	//P' * sin(Theta_p) = K' * sin(Theta_gamma)

	double cosTh_gamma = cos(pTheta_gamma);
	pE_gamma =  pE0 / (1.0 + pE0/pMtg * (1.0-cosTh_gamma) );

	pP_x = sqrt( pE0*pE0 + pE_gamma*pE_gamma - 2*pE0*pE_gamma*cosTh_gamma );

	double sinTheta_x = pE_gamma * sin(pTheta_gamma) / pP_x;
	pTheta_x = asin(sinTheta_x);
}


void   HRSPrimaryGeneratorAction::HRSElasNucleusEngine(int index)
{
	int couple2ThisParticle = (coupleToPrimary[index]>0)?coupleToPrimary[index]-1:0;
	if(index>0 && couple2ThisParticle<index)
	{
		//couple with previous electron
		//Beam + rest_N = e' + N'
		//N' = (Beam + rest_N) - e' = (-P_x_e',-P_y_e',E-P_z_e',E+M_N-E')
		momentum3V[index].setX(-momentum3V[couple2ThisParticle].x());
		momentum3V[index].setY(-momentum3V[couple2ThisParticle].y());
		momentum3V[index].setZ(beamEnergy-momentum3V[couple2ThisParticle].z());
	}
	else
	{
		G4double pTheta,pPhi;

		if(randomizeInTCS[index]>0) RandTCSThetaPhi(index,pTheta,pPhi);
		else RandHCSThetaPhi(index,pTheta,pPhi);

		effectiveTheta[index]=GetEffectiveTheta(position3V.z(),pTheta,pPhi);
		//E-N elastic energy at theta angle
		double pPtot=GetElasMomentumX(beamEnergy,effectiveTheta[index],targetMass);
		momentum3V[index].setRThetaPhi(pPtot,pTheta,pPhi);
	}
}

void   HRSPrimaryGeneratorAction::HRSQuasiElasElectronEngine(int index)
{
	G4double pTheta,pPhi;

	if(randomizeInTCS[index]>0) RandTCSThetaPhi(index,pTheta,pPhi);
	else RandHCSThetaPhi(index,pTheta,pPhi);
	effectiveTheta[index]=GetEffectiveTheta(position3V.z(),pTheta,pPhi);

	//e-p or e-n  elastic energy at theta angle
	//double dMassN = neutron_mass_c2;
	double pMassN = proton_mass_c2;
	double pEprime=beamEnergy/(1.0+ beamEnergy*(1-cos(effectiveTheta[index]))/pMassN );
	double pMass = electron_mass_c2;
	double pPtot=sqrt(pEprime*pEprime-pMass*pMass);

	momentum3V[index].setRThetaPhi(pPtot,pTheta,pPhi);
}

void   HRSPrimaryGeneratorAction::HRSQuasiElasNucleonEngine(int index)
{
	int couple2ThisParticle = (coupleToPrimary[index]>0)?coupleToPrimary[index]-1:0;
	if(index>0 && couple2ThisParticle<index)
	{
		//couple with previous electron
		//Beam + rest_N = e' + N'
		//N' = (Beam + rest_N) - e' = (-P_x_e',-P_y_e',E-P_z_e',E+M_N-E')
		momentum3V[index].setX(-momentum3V[couple2ThisParticle].x());
		momentum3V[index].setY(-momentum3V[couple2ThisParticle].y());
		momentum3V[index].setZ(beamEnergy-momentum3V[couple2ThisParticle].z());
	}
	else
	{
		G4double pTheta,pPhi,pPtot;

		if(randomizeInTCS[index]>0) RandTCSThetaPhi(index,pTheta,pPhi);
		else RandHCSThetaPhi(index,pTheta,pPhi);
		effectiveTheta[index]=GetEffectiveTheta(position3V.z(),pTheta,pPhi);

		//E-N elastic energy at theta angle
		//double pMassN = neutron_mass_c2;
		double pMassN = particle[index]->GetPDGMass();

		//Conclusion: Both rutione return the same Ptot 
		//pPtot=GetElasMomentumX_old(beamEnergy,effectiveTheta[index],pMassN);
		//G4cout<<"Old HRSQuasiElasNucleonEngine return  Theta="<<effectiveTheta[index]/deg<<"  Ptot="<<pPtot<<G4endl;

		pPtot=GetElasMomentumX(beamEnergy,effectiveTheta[index],pMassN);
		//G4cout<<"New HRSQuasiElasNucleonEngine return  Theta="<<effectiveTheta[index]/deg<<"  Ptot="<<pPtot<<G4endl;
		momentum3V[index].setRThetaPhi(pPtot,pTheta,pPhi);
	}
}

void   HRSPrimaryGeneratorAction::TwoBodyEngine(int index)
{
	int couple2ThisParticle = (coupleToPrimary[index]>0)?coupleToPrimary[index]-1:0;
	if(index>0 && couple2ThisParticle<index)
	{
		//couple with previous electron
		//Beam + rest_N = e' + N'
		//N' = (Beam + rest_N) - e' = (-P_x_e',-P_y_e',E-P_z_e',E+M_N-E')
		momentum3V[index].setX(-momentum3V[couple2ThisParticle].x());
		momentum3V[index].setY(-momentum3V[couple2ThisParticle].y());
		momentum3V[index].setZ(beamEnergy-momentum3V[couple2ThisParticle].z());

		//debug:
		//G4cout<<"TwoBodyEngine(): Couple primary "<<index+1<<" to parmary "<<couple2ThisParticle+1<<": \n"
		//	<<"\t momentum3V["<<couple2ThisParticle<<"]="<<momentum3V[couple2ThisParticle]<<G4endl
		//	<<"\t momentum3V["<<index<<"]="<<momentum3V[index]<<G4endl
		//	<<"\t Beam="<<beamEnergy<<G4endl;
	}
	else
	{
		if(index==0)
			cout<<"***Error: There is no exist primary pariticle that particle 1 can couple yet.\n"
			<<"You have to specified another engine for particle 1\n";
		else
			cout<<"***Error: do not know which primary particle to couple. Please given the ordinal \n"
			<<"          number using command  /mydet/particle"<<index+1<<"/coupleToPrimary "<<endl;
	}
}

// using root histo to generate event	
int  HRSPrimaryGeneratorAction::GenPrimariesFromHisto(G4Event* anEvent)
{
	//if file is not open yet or the user has specified a new file,then open it
	if(!mRootHisto->bHistoOpened) 
	{
		mRootHisto->Initialize();
	}	

	TVector3 v3Vertex,v3P_rtpc,v3P_el,v3P_pi,v3P_pr;
	G4ThreeVector G4V3Mom, G4V3Vx;
	//get one exclusive pi- evnt
	if (iRootEventType==3 || iRootEventType==4)
	{
		mRootHisto->GetExclusiveEvent(v3Vertex,v3P_rtpc,v3P_el,v3P_pi,v3P_pr);

		G4V3Vx.set(v3Vertex.x()*cm,v3Vertex.y()*cm,v3Vertex.z()*cm);
		//pr_rtpc
		G4V3Mom.set(v3P_rtpc.x()*GeV,v3P_rtpc.y()*GeV,v3P_rtpc.z()*GeV);
		particleGun->SetParticlePosition(G4V3Vx);
		particleGun->SetParticleDefinition(proton);
		particleGun->SetParticleMomentum(G4V3Mom);	 
		particleGun->GeneratePrimaryVertex(anEvent);

		//el
		G4V3Mom.set(v3P_el.x()*GeV,v3P_el.y()*GeV,v3P_el.z()*GeV);
		particleGun->SetParticlePosition(G4V3Vx);
		particleGun->SetParticleDefinition(electron);
		particleGun->SetParticleMomentum(G4V3Mom);	 
		particleGun->GeneratePrimaryVertex(anEvent);

		//pi-
		G4V3Mom.set(v3P_pi.x()*GeV,v3P_pi.y()*GeV,v3P_pi.z()*GeV);
		particleGun->SetParticlePosition(G4V3Mom);
		particleGun->SetParticleDefinition(pionminus);
		particleGun->SetParticleMomentum(G4V3Mom);	 
		particleGun->GeneratePrimaryVertex(anEvent);

		particleNum=3;

		if(iRootEventType==4)
		{
			//pr_clas
			G4V3Mom.set(v3P_pr.x()*GeV,v3P_pr.y()*GeV,v3P_pr.z()*GeV);
			particleGun->SetParticlePosition(G4V3Vx);
			particleGun->SetParticleDefinition(proton);
			particleGun->SetParticleMomentum(G4V3Mom);	 
			particleGun->GeneratePrimaryVertex(anEvent);

			particleNum=4;
		}
	}


	//get one inclusive event
	else if (iRootEventType==1 || iRootEventType==2)
	{
		mRootHisto->GetInclusiveEvent(v3Vertex,v3P_rtpc,v3P_el);

		G4V3Vx.set(v3Vertex.x()*cm,v3Vertex.y()*cm,v3Vertex.z()*cm);
		//pr_rtpc
		G4V3Mom.set(v3P_rtpc.x()*GeV,v3P_rtpc.y()*GeV,v3P_rtpc.z()*GeV);
		particleGun->SetParticlePosition(G4V3Vx);
		particleGun->SetParticleDefinition(proton);
		particleGun->SetParticleMomentum(G4V3Mom);	 
		particleGun->GeneratePrimaryVertex(anEvent);

		particleNum=1;

		if(iRootEventType==2)
		{
			//el
			G4V3Mom.set(v3P_el.x()*GeV,v3P_el.y()*GeV,v3P_el.z()*GeV);
			particleGun->SetParticlePosition(G4V3Vx);
			particleGun->SetParticleDefinition(electron);
			particleGun->SetParticleMomentum(G4V3Mom);	 
			particleGun->GeneratePrimaryVertex(anEvent);

			particleNum=2;
		}

	}

	return particleNum;
}

// randomize events from root histo	
void  HRSPrimaryGeneratorAction::HRSHistoEngine(int index)
{
	index=0;
	//still need to debug
}

// reading events from root file	
void  HRSPrimaryGeneratorAction::HRSNtupleEngine(int index)
{
	//if file is not open yet or the user has specified a new file,then open it

	if(!mRootEvent->mNtuple[index]) 
	{
		char tmpStr[255];
		sprintf(tmpStr,"InRootFileName%d",index+1);
		string pInRootFileName = gConfig->GetArgument(tmpStr);
		sprintf(tmpStr,"InRootTreeName%d",index+1);
		string pInRootTreeName = gConfig->GetArgument(tmpStr);
		int pSkipNum=0,pTrigNum=-1;
		sprintf(tmpStr,"SkipNum%d",index+1);
		gConfig->GetArgument(tmpStr,pSkipNum);
		sprintf(tmpStr,"TrigNum%d",index+1);
		gConfig->GetArgument(tmpStr,pTrigNum);

		mRootEvent->LoadNtuple(index,pInRootFileName.c_str(),pInRootTreeName.c_str(),pTrigNum,pSkipNum);
	}	


	//get vertex and get momentum
	int pdgid=0;
	TVector3 v3x,v3p;

	int ret=mRootEvent->GetParticle(index,pdgid,v3x,v3p);
	if(ret>=0)
	{
		position3V.set(v3x.x()*mm,v3x.y()*mm,v3x.z()*mm);
		momentum3V[index].set(v3p.x()*GeV,v3p.y()*GeV,v3p.z()*GeV);
		if(particle[index]->GetPDGEncoding() != pdgid)
		{
			particle[index]=G4ParticleTable::GetParticleTable()->FindParticle(pdgid);	
		}
	}
	else if(ret==-1)
	{
		//end of file, terminate this run
		G4RunManager* theRunManager=G4RunManager::GetRunManager(); 
		theRunManager->AbortRun();
	}

}


// reading events from root file	
int  HRSPrimaryGeneratorAction::GenPrimariesFromNtuple(G4Event* anEvent)
{
	//if file is not open yet or the user has specified a new file,then open it
	for(int i=0;i<particleNum;i++)
	{
		if(!mRootEvent->mNtuple[i]) 
		{
			char tmpStr[255];
			sprintf(tmpStr,"InRootFileName%d",i+1);
			string pInRootFileName = gConfig->GetArgument(tmpStr);
			sprintf(tmpStr,"InRootTreeName%d",i+1);
			string pInRootTreeName = gConfig->GetArgument(tmpStr);
			int pSkipNum=0,pTrigNum=-1;
			sprintf(tmpStr,"SkipNum%d",i+1);
			gConfig->GetArgument(tmpStr,pSkipNum);
			sprintf(tmpStr,"TrigNum%d",i+1);
			gConfig->GetArgument(tmpStr,pTrigNum);

			mRootEvent->LoadNtuple(i,pInRootFileName.c_str(),pInRootTreeName.c_str(),pTrigNum,pSkipNum);
		}	
	}

	//get vertex and get momentum
	int pdgid=0;
	TVector3 v3x,v3p;
	G4ThreeVector G4V3Mom, G4V3Vx;

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	int pPartNum=0;
	for(int trackid=0;trackid<particleNum;trackid++)
	{
		int ret=mRootEvent->GetParticle(trackid,pdgid,v3x,v3p);

		if(ret>=0)
		{
			G4V3Vx.set(v3x.x()*mm,v3x.y()*mm,v3x.z()*mm);
			G4V3Mom.set(v3p.x()*GeV,v3p.y()*GeV,v3p.z()*GeV);

			particleGun->SetParticleDefinition(particleTable->FindParticle(pdgid));
			particleGun->SetParticlePosition(G4V3Vx);
			particleGun->SetParticleMomentum(G4V3Mom);	 
			particleGun->GeneratePrimaryVertex(anEvent);
			pPartNum++;
		}
	}

	particleNum=pPartNum;
	return particleNum;
}

