
#include "HRSNtReader.hh"


HRSNtReader::HRSNtReader(const char *filename, const char *treename, int loadconfigtree)
{
	mChain=new TChain(treename);
	char fullpath[255];
	sprintf(fullpath,"%s/%s",filename,treename);
	int n=mChain->Add(fullpath);
	if(n>0) 
	{
		if(loadconfigtree) ReadConfig(filename);
		Init();
		mEntries=mChain->GetEntries(); 
		cout<<"HRSNtReader is loading ntuple "<<filename<<"/"<<treename<<endl;
	}
	else
	{
		cout<<"HRSNtReader fail to load ntuple "<<filename<<".  Exit..."<<endl;
		exit(-1);
	}
}

HRSNtReader::~HRSNtReader()
{
	if (!mChain) return;
	delete mChain;
}

Bool_t HRSNtReader::ReadConfig(const char *filename)
{
	TChain *config = new TChain("config");
	char fullpath[255];
	sprintf(fullpath,"%s/config",filename);
	config->Add(fullpath);

	// Set branch addresses.
	config->SetBranchAddress("Run",&Run);
	config->SetBranchAddress("SkimLevel",&SkimLevel);
	config->SetBranchAddress("BookTrees",&BookTrees);
	config->SetBranchAddress("Beam",&Beam);
	config->SetBranchAddress("TargetM",&TargetM);
	config->SetBranchAddress("TargetAtomicNumber",&TargetAtomicNumber);
	config->SetBranchAddress("LHRSAngle",&LHRSAngle);
	config->SetBranchAddress("RHRSAngle",&RHRSAngle);
	config->SetBranchAddress("TargetXOffset",&TargetXOffset);
	config->SetBranchAddress("TargetYOffset",&TargetYOffset);
	config->SetBranchAddress("TargetZOffset",&TargetZOffset);
	config->SetBranchAddress("PivotXOffset",&PivotXOffset);
	config->SetBranchAddress("PivotYOffset",&PivotYOffset);
	config->SetBranchAddress("PivotZOffset",&PivotZOffset);
	config->SetBranchAddress("LHRSMomentum",&LHRSMomentum);
	config->SetBranchAddress("RHRSMomentum",&RHRSMomentum);
	config->SetBranchAddress("UseHelmField",&UseHelmField);
	config->SetBranchAddress("HelmXOffset",&HelmXOffset);
	config->SetBranchAddress("HelmYOffset",&HelmYOffset);
	config->SetBranchAddress("HelmZOffset",&HelmZOffset);
	config->SetBranchAddress("HelmRotAxis1",&HelmRotAxis1);
	config->SetBranchAddress("HelmRotAxis2",&HelmRotAxis2);
	config->SetBranchAddress("HelmRotAxis3",&HelmRotAxis3);
	config->SetBranchAddress("HelmRotAngle1",&HelmRotAngle1);
	config->SetBranchAddress("HelmRotAngle2",&HelmRotAngle2);
	config->SetBranchAddress("HelmRotAngle3",&HelmRotAngle3);
	config->SetBranchAddress("HelmCurrentRatio",&HelmCurrentRatio);
	config->SetBranchAddress("UseSeptumField",&UseSeptumField);
	config->SetBranchAddress("SeptumXOffset",&SeptumXOffset);
	config->SetBranchAddress("SeptumYOffset",&SeptumYOffset);
	config->SetBranchAddress("SeptumZOffset",&SeptumZOffset);
	config->SetBranchAddress("SeptumRotAxis1",&SeptumRotAxis1);
	config->SetBranchAddress("SeptumRotAxis2",&SeptumRotAxis2);
	config->SetBranchAddress("SeptumRotAxis3",&SeptumRotAxis3);
	config->SetBranchAddress("SeptumRotAngle1",&SeptumRotAngle1);
	config->SetBranchAddress("SeptumRotAngle2",&SeptumRotAngle2);
	config->SetBranchAddress("SeptumRotAngle3",&SeptumRotAngle3);
	config->SetBranchAddress("SeptumCurrentRatioL",&SeptumCurrentRatioL);
	config->SetBranchAddress("SeptumCurrentRatioR",&SeptumCurrentRatioR);
	config->SetBranchAddress("BigBiteAngle",&BigBiteAngle);
	config->SetBranchAddress("BigBiteTiltAngle",&BigBiteTiltAngle);
	config->SetBranchAddress("Pivot2BigBiteFace",&Pivot2BigBiteFace);

	Long64_t nentries = config->GetEntries();

	if(nentries==0) return false;
	mIsCombinedTree=false;
	double tmpAtg,tmpE,tmpZOff;
	for (Long64_t i=0; i<nentries;i++) 
	{
		config->GetEntry(i);
		//check if this is a combined tree with various beam energies
		if(i>0 && (fabs(tmpAtg-TargetAtomicNumber)>0.01 || fabs(tmpE-Beam)>0.01 ||
			fabs(tmpZOff-PivotZOffset)>30.0 ) )  
		{
			mIsCombinedTree=true;
			break;
		}
		else
		{
			tmpAtg=TargetAtomicNumber;
			tmpE=Beam;
			tmpZOff=PivotZOffset;
		}
	}
	delete config;

	return mIsCombinedTree;
}

void HRSNtReader::Init()
{
	// Set branch addresses 
	if (!mChain) return;

	mChain->SetBranchAddress("Index",&Index);
	mChain->SetBranchAddress("PdgId",&PdgId);
	mChain->SetBranchAddress("TrackId",&TrackId);
	mChain->SetBranchAddress("TrackClass",&TrackClass);
	mChain->SetBranchAddress("X0",&X0);
	mChain->SetBranchAddress("Y0",&Y0);
	mChain->SetBranchAddress("Z0",&Z0);
	mChain->SetBranchAddress("P0",&P0);
	mChain->SetBranchAddress("Theta0",&Theta0);
	mChain->SetBranchAddress("Phi0",&Phi0);
	mChain->SetBranchAddress("X0_tr",&X0_tr);
	mChain->SetBranchAddress("Y0_tr",&Y0_tr);
	mChain->SetBranchAddress("Z0_tr",&Z0_tr);
	mChain->SetBranchAddress("Theta0_tr",&Theta0_tr);
	mChain->SetBranchAddress("Phi0_tr",&Phi0_tr);
	mChain->SetBranchAddress("Xvb",&Xvb);
	mChain->SetBranchAddress("Yvb",&Yvb);
	mChain->SetBranchAddress("Zvb",&Zvb);
	mChain->SetBranchAddress("Pvb",&Pvb);
	mChain->SetBranchAddress("Thetavb",&Thetavb);
	mChain->SetBranchAddress("Phivb",&Phivb);
	mChain->SetBranchAddress("Xvb_tr",&Xvb_tr);
	mChain->SetBranchAddress("Yvb_tr",&Yvb_tr);
	mChain->SetBranchAddress("Zvb_tr",&Zvb_tr);
	mChain->SetBranchAddress("Thetavb_tr",&Thetavb_tr);
	mChain->SetBranchAddress("Phivb_tr",&Phivb_tr);
	mChain->SetBranchAddress("Xfp_tr",&Xfp_tr);
	mChain->SetBranchAddress("Yfp_tr",&Yfp_tr);
	mChain->SetBranchAddress("Thetafp_tr",&Thetafp_tr);
	mChain->SetBranchAddress("Phifp_tr",&Phifp_tr);
	mChain->SetBranchAddress("Delta",&Delta);
	mChain->SetBranchAddress("Delta_rec",&Delta_rec);
	mChain->SetBranchAddress("X_rec_tr",&X_rec_tr);
	mChain->SetBranchAddress("Y_rec_tr",&Y_rec_tr);
	mChain->SetBranchAddress("Theta_rec_tr",&Theta_rec_tr);
	mChain->SetBranchAddress("Phi_rec_tr",&Phi_rec_tr);
	mChain->SetBranchAddress("Z_rec",&Z_rec);
	mChain->SetBranchAddress("P_rec",&P_rec);
	mChain->SetBranchAddress("Theta_rec",&Theta_rec);
	mChain->SetBranchAddress("Phi_rec",&Phi_rec);
	mChain->SetBranchAddress("TrackRadlen",&TrackRadlen);
	mChain->SetBranchAddress("Theta0Eff",&Theta0Eff);
	mChain->SetBranchAddress("StepNum",&StepNum);
	mChain->SetBranchAddress("StepX",StepX);
	mChain->SetBranchAddress("StepY",StepY);
	mChain->SetBranchAddress("StepZ",StepZ);
	mChain->SetBranchAddress("StepdE",StepdE);
	mChain->SetBranchAddress("StepL",StepL);
	mChain->SetBranchAddress("StepEkin",StepEkin);
	mChain->SetBranchAddress("StepTL",StepTL);
	mChain->SetBranchAddress("StepRadlen",StepRadlen);
	mChain->SetBranchAddress("StepDsty",StepDsty);
	mChain->SetBranchAddress("StepBx",StepBx);
	mChain->SetBranchAddress("StepBy",StepBy);
	mChain->SetBranchAddress("StepBz",StepBz);
	mChain->SetBranchAddress("TrackBdLx",&TrackBdLx);
	mChain->SetBranchAddress("TrackBdLy",&TrackBdLy);
	mChain->SetBranchAddress("TrackBdLz",&TrackBdLz);
	mChain->SetBranchAddress("R0",&R0);
	mChain->SetBranchAddress("A0",&A0);
	mChain->SetBranchAddress("B0",&B0);
	
	if(mChain->GetBranch("ElasXS"))
	{
		mChain->SetBranchAddress("ElasXS",&ElasXS);
		mChain->SetBranchAddress("XS",&XS);
	}
}

Int_t HRSNtReader::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!mChain) return 0;
	return mChain->GetEntry(entry);
}

Int_t HRSNtReader::Cut()
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	if(P0<0) return 0;
	return 1;
}


void HRSNtReader::Loop()
{
   if (mChain == 0) return;

   Long64_t nentries = mChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      nb = mChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
