
//Declaration of leaves types tree
   Int_t           Run;
   Int_t           SkimLevel;
   Int_t           BookTrees;
   Double_t        Beam;
   Double_t        TargetM;
   Double_t        TargetAtomicNumber;
   Double_t        LHRSAngle;
   Double_t        RHRSAngle;
   Double_t        TargetXOffset;
   Double_t        TargetYOffset;
   Double_t        TargetZOffset;
   Double_t        PivotXOffset;
   Double_t        PivotYOffset;
   Double_t        PivotZOffset;
   Double_t        LHRSMomentum;
   Double_t        RHRSMomentum;
   Int_t           UseHelmField;
   Double_t        HelmXOffset;
   Double_t        HelmYOffset;
   Double_t        HelmZOffset;
   Double_t        HelmRotAxis1;
   Double_t        HelmRotAxis2;
   Double_t        HelmRotAxis3;
   Double_t        HelmRotAngle1;
   Double_t        HelmRotAngle2;
   Double_t        HelmRotAngle3;
   Double_t        HelmCurrentRatio=1.0;
   Int_t           UseSeptumField;
   Double_t        SeptumXOffset;
   Double_t        SeptumYOffset;
   Double_t        SeptumZOffset;
   Double_t        SeptumRotAxis1;
   Double_t        SeptumRotAxis2;
   Double_t        SeptumRotAxis3;
   Double_t        SeptumRotAngle1;
   Double_t        SeptumRotAngle2;
   Double_t        SeptumRotAngle3;
   Double_t        SeptumCurrentRatio;
   Double_t        BigBiteAngle;
   Double_t        BigBiteTiltAngle;
   Double_t        Pivot2BigBiteFace;

	bool			bIsCombinedTree;	//to ideneify if this is a combined ntuple



bool ReadConfig()
{
	TTree *config = (TTree*)gDirectory->Get("config");

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
   if(config->GetBranch("HelmCurrentRatio"))
   {
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
   config->SetBranchAddress("SeptumCurrentRatio",&SeptumCurrentRatio);
   config->SetBranchAddress("BigBiteAngle",&BigBiteAngle);
   config->SetBranchAddress("BigBiteTiltAngle",&BigBiteTiltAngle);
   config->SetBranchAddress("Pivot2BigBiteFace",&Pivot2BigBiteFace);
   }
	//     This is the loop skeleton
	//       To read only selected branches, Insert statements like:
	// config->SetBranchStatus("*",0);  // disable all branches
	// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

	Long64_t nentries = config->GetEntries();

	bIsCombinedTree=false;
	if(nentries==0) return false;
	Long64_t nbytes = 0;
	double tmpAtg,tmpE,tmpZOff;
	for (Long64_t i=0; i<nentries;i++) 
	{
		nbytes += config->GetEntry(i);
		//check if this is a combined tree with various beam energies
		if(i>0 && (fabs(tmpAtg-TargetAtomicNumber)>0.01 || fabs(tmpE-Beam)>0.01 ||
		fabs(tmpZOff-PivotZOffset)>30.0 ) )  
		{
			bIsCombinedTree=true;
			break;
		}
		else
		{
			tmpAtg=TargetAtomicNumber;
			tmpE=Beam;
			tmpZOff=PivotZOffset;
		}
	}
	return true;
}


void BeamShift(char *filename="",double EndPlaneZ=640, double Range=10)
{
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(0);

	if(strlen(filename)>5)
	{
		cout<<"trying to open root file "<<filename<<endl;
		TFile* _file0=TFile::Open(filename);
	}
	ReadConfig();

	TF1 *f1 = new TF1("f1","[0]/x+[1]+[2]*x",0.2,3.3);
	f1->SetFillColor(0);
	f1->SetMarkerStyle(20);
	f1->SetMarkerColor(2);
	f1->SetLineColor(2);
	f1->SetLineWidth(3);
	
	char EndPlaneCut[1024];
	sprintf(EndPlaneCut,"abs(StepZ-(%.1f+%.1f))<%.1f", TargetZOffset, EndPlaneZ,Range);
	TCanvas *pCan=new TCanvas("pCan","Check Beam shift",800,600);
	
	pCan->cd();
	track0->Draw("StepY-Y0:P0>>h2Yshift",EndPlaneCut,"prof");
	h2Yshift->SetMarkerStyle(20);
	h2Yshift->SetMarkerColor(4);
	h2Yshift->SetLineColor(4);
	h2Yshift->SetTitle(Form("Y_{shift} for HRS=%.1f^{o}, Z=%.0f mm, Y_{offset}=%.0f mm, R_{field}=%.1f; P0 (GeV); Y(mm)",
		 LHRSAngle*180./3.14159,EndPlaneZ,TargetYOffset,HelmCurrentRatio));
	h2Yshift->Fit(f1,"R","", 0.55,3.4);
	h2Yshift->Draw("prof");
	f1->Draw("same");
	
	TPaveText *pt2 = new TPaveText(0.65,0.315,0.898,0.59,"brNDC");
	pt2->SetBorderSize(1);
	pt2->SetFillColor(0);
    pt2->SetTextAlign(12);
	pt2->AddText("f(x)=p0/x+p1+p2*x");	
	pt2->AddText(Form("p0 = %f",f1->GetParameter(0)));	
	pt2->AddText(Form("p1 = %f",f1->GetParameter(1)));	
	pt2->AddText(Form("p1 = %f",f1->GetParameter(2)));	
	pt2->Draw();

	pCan->SaveAs(Form("BeamShift_Z%.0f_R%.1f.gif",EndPlaneZ,HelmCurrentRatio));
 }
 

void SeptumShift(char *filename="", double Range=0.5)
{
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(0);

	if(strlen(filename)>5)
	{
		cout<<"trying to open root file "<<filename<<endl;
		TFile* _file0=TFile::Open(filename);
	}
	ReadConfig();

	TF1 *f1fit = new TF1("f1fit","[0]+[1]*x",900.0,2000.0);
	f1fit->SetFillColor(0);
	f1fit->SetMarkerStyle(20);
	f1fit->SetMarkerColor(7);
	f1fit->SetLineColor(7);
	f1fit->SetLineWidth(3);
	
	TF1 *f6 = new TF1("f6",Form("tan(6.0*3.14159/180.)*(x-(%.2f))",TargetZOffset),-900.0,1281.0);
	f6->SetFillColor(0);
	f6->SetMarkerStyle(20);
	f6->SetMarkerColor(6);
	f6->SetLineColor(6);
	f6->SetLineWidth(3);
	
	TF1 *f12 = new TF1("f12","tan(12.5*3.14159/180.)*x",0.0,2000.0);
	f12->SetFillColor(0);
	f12->SetMarkerStyle(20);
	f12->SetMarkerColor(2);
	f12->SetLineColor(2);
	f12->SetLineWidth(3);
	
	TF1 *f12fit = new TF1("f12","tan(12.5*3.14159/180.)*(x+[0])",900.0,2000.0);
	f12fit->SetFillColor(0);
	f12fit->SetMarkerStyle(20);
	f12fit->SetMarkerColor(7);
	f12fit->SetLineColor(7);
	f12fit->SetLineWidth(3);
	
	
	char TheCut[1024];
	sprintf(TheCut,Form("abs(Thetavb*180./3.14159-12.5)<%.1f && abs(P0-Pvb)/P0<0.0005",Range));
	
	TCanvas *pCan=new TCanvas("pCan","Check steptum shift",700,920);
	pCan->Divide(1,2);
	
	//Q1 entrance Z=1281
	pCan->cd(1);
	track0->Draw("abs(StepX):StepZ>>h2XZ",TheCut,"prof");
	h2XZ->SetMarkerStyle(20);
	h2XZ->SetMarkerColor(4);
	h2XZ->SetLineColor(4);
	h2XZ->SetTitle(Form("X Vs Z: R_{field}=%.1f, Z_{septum}=%.2f mm; Z (mm); X(mm)",
		SeptumCurrentRatio,SeptumZOffset*10));
	h2XZ->Fit(f12fit,"R","", 1200,2000);
	h2XZ->Draw("prof");
	f6->Draw("same");
	f12->Draw("same");
	
	TPaveText *pt1 = new TPaveText(0.12,0.61,0.50,0.89,"brNDC");
	pt1->SetBorderSize(1);
	pt1->SetFillColor(0);
    pt1->SetTextAlign(12);
	double z0=f12fit->GetParameter(0);
	pt1->AddText("f(z)=tan(12.5*deg)*(z+Z_{off})");	
	pt1->AddText(Form("Z_{off} = %.2f",z0));
	pt1->AddText(Form("X_{off} = %.2f", z0* tan(12.5*3.14159/180))); 
	pt1->AddText(Form("#color[2]{Shift at Q1 plane = %.2f}", z0* sin(12.5*3.14159/180))); 
	pt1->Draw();

	pCan->cd(2);
	track0->Draw("Thetavb*180/3.14159:Pvb>>h2ThetaPvb",TheCut,"prof");
	h2ThetaPvb->SetMarkerStyle(20);
	h2ThetaPvb->SetMarkerColor(4);
	h2ThetaPvb->SetLineColor(4);
	h2ThetaPvb->SetTitle(Form("#theta_{vb} Vs P_{vb}: R_{field}=%.1f, Z_{septum}=%.2f mm; P_{vb} (GeV); #theta_{vb} (mm)",
		SeptumCurrentRatio,SeptumZOffset*10));
	h2ThetaPvb->Fit(f1fit,"R","", 3.02,3.18);
	h2ThetaPvb->GetYaxis()->SetRangeUser(11.5,13.5);
	h2ThetaPvb->Draw("prof");
	f1fit->Draw("same");
	
	TPaveText *pt2 = new TPaveText(0.15,0.14,0.45,0.40,"brNDC");
	pt2->SetBorderSize(1);
	pt2->SetFillColor(0);
    pt2->SetTextAlign(12);
	pt2->AddText("f(x)= p0 + p1*x");	
	double *pp=f1fit->GetParameters();
	pt2->AddText(Form("p0 = %f",pp[0]));
	pt2->AddText(Form("p1 = %f",pp[1]));
	//y=f(x)=p0+p1*x ==> f^{-1}(y)=(y-p0)/p1
	double xx = (12.5-pp[0])/pp[1];
	pt2->AddText(Form("#color[2]{Central Ray P = %.4f}",xx));	
	pt2->Draw();
	
	pCan->SaveAs(Form("SeptumShift_R%.1fZ%.2f.gif",SeptumCurrentRatio,SeptumZOffset));
	pCan->SaveAs(Form("SeptumShift_R%.1fZ%.2f.C",SeptumCurrentRatio,SeptumZOffset));
 }
 
void CheckBeamShift(char *filename="", double Range=10)
{
	BeamShift(filename,640,Range);
	BeamShift(filename,1200,Range);
	SeptumShift();
}

