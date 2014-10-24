
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
	h2Yshift->SetTitle(Form("Y_{shift} for beam, Z=%.0f mm, Y_{offset}=%.0f mm, R_{field}=%.3f; P0 (GeV); Y(mm)",
		 EndPlaneZ,TargetYOffset,HelmCurrentRatio));
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

	pCan->SaveAs(Form("BeamShift_R%.1f_Helm%.0f_Z%.0f.png",HelmCurrentRatio,HelmRotAngle1*180/3.14159,EndPlaneZ));
 }
 

void YShift(char *filename="",double EndPlaneZ=640, double Range=10)
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
	TCanvas *pCan=new TCanvas("pCan","Check Y shift",800,600);
	
	pCan->cd();
	track0->Draw("StepY-Y0:P0>>h2Yshift",EndPlaneCut,"prof");
	h2Yshift->SetMarkerStyle(20);
	h2Yshift->SetMarkerColor(4);
	h2Yshift->SetLineColor(4);
	h2Yshift->SetTitle(Form("Y_{shift} for HRS=%.2f^{o}, Z=%.0f mm, Y_{offset}=%.0f mm, R_{field}=%.3f; P0 (GeV); Y(mm)",
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

	pCan->SaveAs(Form("YShift_HRS%.2f_R%.1f_Helm%.0f_Z%.0f.png",LHRSAngle*180./3.14159,HelmCurrentRatio,HelmRotAngle1*180/3.14159,EndPlaneZ));
 }
 
void SeptumShift(char *filename="", double Range=0.2)
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
	h2XZ->SetTitle(Form("X Vs Z: R_{field}=%.3f, Z_{septum}=%.2f mm; Z (mm); X(mm)",
		SeptumCurrentRatio,SeptumZOffset));
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
	track0->Draw("Theta0*180/3.14159>>hTheta0",TheCut,"");
	double Theta0Mean=hTheta0->GetMean();
	track0->Draw("Z0>>hZ0",TheCut,"");
	double Z0Mean=hZ0->GetMean();
	
	track0->Draw("Thetavb*180/3.14159:Pvb>>h2ThetaPvb",TheCut,"prof");
	h2ThetaPvb->SetMarkerStyle(20);
	h2ThetaPvb->SetMarkerColor(4);
	h2ThetaPvb->SetLineColor(4);
	h2ThetaPvb->SetTitle(Form("#theta_{vb} Vs P_{vb}: R_{field}=%.3f, Z_{septum}=%.2f mm; P_{vb} (GeV); #theta_{vb} (mm)",
		SeptumCurrentRatio,SeptumZOffset));
	h2ThetaPvb->Fit(f1fit,"R","", 2.70,3.100);
	//h2ThetaPvb->Fit(f1fit);
	h2ThetaPvb->GetYaxis()->SetRangeUser(12.4,12.6);
	h2ThetaPvb->Draw("prof");
	f1fit->Draw("same");
	
	TPaveText *pt2 = new TPaveText(0.12,0.12,0.47,0.42,"brNDC");
	pt2->SetBorderSize(1);
	pt2->SetFillColor(0);
    pt2->SetTextAlign(12);
	pt2->AddText(Form("Z_{0}=%.2fmm,  #theta_{0}=%.2f^{o}",Z0Mean,Theta0Mean));
	pt2->AddText("f(x)= p0 + p1*x");	
	double *pp=f1fit->GetParameters();
	pt2->AddText(Form("p0 = %f",pp[0]));
	pt2->AddText(Form("p1 = %f",pp[1]));
	//y=f(x)=p0+p1*x ==> f^{-1}(y)=(y-p0)/p1
	double xx = (12.5-pp[0])/pp[1];
	pt2->AddText(Form("#color[2]{Central Ray P = %.4f}",xx));	
	pt2->Draw();
	
	pCan->SaveAs(Form("SeptumShift_R%.3fZ%.2f.png",SeptumCurrentRatio,SeptumZOffset));
	pCan->SaveAs(Form("SeptumShift_R%.3fZ%.2f.C",SeptumCurrentRatio,SeptumZOffset));
 }
 

void AngleShift(char *filename="", int TrackClassCut=0)
{
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(10);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetStatH(0.2);
	gStyle->SetStatW(0.17);
	//gStyle->SetStatStyle(4000);
	//gStyle->SetStatX(0.43);
	gStyle->SetStatX(0.95);
	gStyle->SetStatY(0.65);
	

	if(strlen(filename)>5)
	{
		cout<<"trying to open root file "<<filename<<endl;
		TFile* _file0=TFile::Open(filename);
	}
	ReadConfig();

	TF1 *f1 = new TF1("f1","[0]/x+[1]+[2]*x",0.3,3.3);
	f1->SetFillColor(0);
	f1->SetMarkerStyle(20);
	f1->SetMarkerColor(4);
	f1->SetLineColor(4);
	f1->SetLineWidth(3);

	TText *text=0;
	char HRSCut[255];
	sprintf(HRSCut,"abs((P0-Pvb)/P0)<0.01 && TrackClass>=%d",TrackClassCut);

	TCanvas *pCan=new TCanvas("pCan","Theta0_tr in G2P",800,600);
	gPad->SetRightMargin(0.05);
	
	track0->Draw("Theta0*180/3.14159>>hTheta0(100,0,20)",HRSCut,"");
	double pTheta0Mean=hTheta0->GetMean();
	bool pIsBeam=(pTheta0Mean<0.5)?true:false;
	
	if(pIsBeam)
	{
	track0->Draw("Thetavb*180/3.14159:P0>>h2ThetavbP0(175,0.,3.50,80,0,80)",HRSCut,"prof");
	h2ThetavbP0->SetTitle(Form("#delta#theta Vs P_{0} for beam electron, R_{helm}=%.1f; P_{0} (GeV); #delta#theta (deg)",HelmCurrentRatio));
	}
	else
	{
	track0->Draw("Thetavb_tr*180/3.14159:P0>>h2ThetavbP0(175,0.,3.50,80,0,80)",HRSCut,"prof");
	h2ThetavbP0->SetTitle(Form("#delta#theta Vs P_{0} for scattered electron, R_{helm}=%.1f; P_{0} (GeV); #delta#theta (deg)",HelmCurrentRatio));
	}
	
	h2ThetavbP0->SetMarkerStyle(20);
	h2ThetavbP0->SetMarkerColor(1);
	h2ThetavbP0->GetYaxis()->SetRangeUser(0.0,25.0);
	h2ThetavbP0->Draw("prof");
	h2ThetavbP0->Fit(f1,"R","", 0.5,3.3);
	
	TPaveText *pt2 = new TPaveText(0.65,0.65,0.95,0.73,"brNDC");
	pt2->SetBorderSize(0);
	pt2->SetFillColor(0);
 //pt2->SetFillStyle(4000);  ///////////////////////////////////////////
	text=pt2->AddText("f(x)=p0/x+p1+p2*x");
	text->SetTextColor(4);
	pt2->Draw();

	pCan->cd();
	if(!pIsBeam)
		pCan->SaveAs(Form("AngleShift_%.2fdeg_R%.1f_Helm%.0f.png",pTheta0Mean,HelmCurrentRatio,HelmRotAngle1*180/3.14159));
	else
		pCan->SaveAs(Form("BeamAngleShift_R%.1f_Helm%.0f.png",HelmCurrentRatio,HelmRotAngle1*180/3.14159));
 }

 
void BeamVerticalShift(const char *filename="")
{
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(10);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetStatH(0.2);
	gStyle->SetStatW(0.17);
	//gStyle->SetStatStyle(4000);
	//gStyle->SetStatX(0.43);
	gStyle->SetStatX(0.95);
	gStyle->SetStatY(0.65);
	

	if(strlen(filename)>5)
	{
		cout<<"trying to open root file "<<filename<<endl;
		TFile* _file0=TFile::Open(filename);
	}
	ReadConfig();

	int TrackClassCut=0;
	char HRSCut[255];
	sprintf(HRSCut,"abs((P0-Pvb)/P0)<0.05 && TrackClass>=%d",TrackClassCut);

	TCanvas *c1=new TCanvas("c1","Theta0_tr in G2P",1000,700);
	gPad->SetRightMargin(0.05);
	
	track0->Draw("Theta0*180/3.14159>>hTheta0(100,0,20)",HRSCut,"");
	double pTheta0Mean=hTheta0->GetMean();
	c1->Clear();
	c1->Divide(2,2);	
	c1->cd(1);
	track0->Draw("StepY:StepZ>>h2YE1","TrackClass>=0 && abs(P0-1.159)/1.159<0.002","*");
	h2YE1->SetTitle(Form("E=1.159GeV, BeamAngle=%.1f^{o}; StepZ (mm);StepY (mm);",pTheta0Mean));
	h2YE1->Draw();
	
	c1->cd(2);
	track0->Draw("StepY:StepZ>>h2YE2","TrackClass>=0 && abs(P0-1.706)/1.706<0.002","*");
	h2YE2->SetTitle(Form("E=1.706GeV, BeamAngle=%.1f^{o}; StepZ (mm);StepY (mm);",pTheta0Mean));
	h2YE2->Draw();

	c1->cd(3);
	track0->Draw("StepY:StepZ>>h2YE3","TrackClass>=0 && abs(P0-2.257)/2.257<0.002","*");
	h2YE3->SetTitle(Form("E=2.257GeV, BeamAngle=%.1f^{o}; StepZ (mm);StepY (mm);",pTheta0Mean));
	h2YE3->Draw();

	c1->cd(4);
	track0->Draw("StepY:StepZ>>h2YE4","TrackClass>=0 && abs(P0-3.355)/3.355<0.002","*");
	h2YE4->SetTitle(Form("E=3.355GeV, BeamAngle=%.1f^{o}; StepZ (mm);StepY (mm);",pTheta0Mean));
	h2YE4->Draw();

	c1->UseCurrentStyle();
	c1->SaveAs(Form("BeamYShift_tilted%.1f.png",pTheta0Mean));
	c1->SaveAs(Form("BeamYShift_tilted%.1f.C",pTheta0Mean));

}
 

void DeltaP(char *filename="", int TrackClassCut=0)
{
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(10);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetStatH(0.2);
	gStyle->SetStatW(0.17);
	//gStyle->SetStatStyle(4000);
	//gStyle->SetStatX(0.43);
	gStyle->SetStatX(0.95);
	gStyle->SetStatY(0.65);
	

	if(strlen(filename)>5)
	{
		cout<<"trying to open root file "<<filename<<endl;
		TFile* _file0=TFile::Open(filename);
	}
	ReadConfig();

	TText *text=0;
	char HRSCut[255];
	sprintf(HRSCut,"abs((P0-Pvb)/P0)<0.1 && TrackClass>=%d",TrackClassCut);

	TCanvas *pCan=new TCanvas("pCan","Theta0_tr in G2P",800,600);
	
	pCan->cd();
	gPad->SetRightMargin(0.05);
	
	track0->Draw("P0>>hP0",HRSCut,"");
	double pP0Mean=hP0->GetMean();
	track0->Draw("Theta0*180/3.14159>>hTheta0",HRSCut,"");
	double pTheta0Mean=hTheta0->GetMean();
	
	track0->Draw("P0-Pvb>>hdP(100,0,0.01)",HRSCut,"");
	hdP->SetTitle(Form("#deltaP for beam, E=%.3fGeV, Angle=%.1f^{o}, R_{helm}=%.1f; #deltaP (GeV/c)",pP0Mean,pTheta0Mean,HelmCurrentRatio));
	
	pCan->SaveAs(Form("BeamELoss_E%.3f_R%.1f_Helm%.0f.png",pP0Mean,HelmCurrentRatio,HelmRotAngle1*180/3.14159));
 }

void CheckBeamShift(char *filename="", double Range=10)
{
	//BeamShift(filename,640,Range);
	//BeamVerticalShift();
	//YShift(filename,640,Range);
	//YShift(filename,790,Range);
	//AngleShift(filename,0);
	SeptumShift(filename);
	//DeltaP(filename);
}

