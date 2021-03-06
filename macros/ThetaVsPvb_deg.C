

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
   Double_t        HelmCurrentRatio;
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

void Phi0VsPvb(char *filename="")
{
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetStatH(0.2);
	gStyle->SetStatW(0.17);
	gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
	

	if(strlen(filename)>5)
	{
		cout<<"trying to open root file "<<filename<<endl;
		TFile* _file0=TFile::Open(filename);
	}
	ReadConfig();

	TF1 *f1 = new TF1("f1","[0]/x+[1]+[2]*x+[3]*x*x",0.2,3.3);
	f1->SetFillColor(0);
	f1->SetMarkerStyle(20);
	f1->SetMarkerColor(2);
	f1->SetLineColor(2);
	f1->SetLineWidth(3);
	
	TF1 *f2 = new TF1("f2","[0]/x+[1]+[2]*x+[3]*x*x",0.2,3.3);
	f2->SetFillColor(0);
	f2->SetMarkerStyle(20);
	f2->SetMarkerColor(4);
	f2->SetLineColor(4);
	f2->SetLineWidth(3);

	TText *text=0;
	char HRSCut[255];
	sprintf(HRSCut,"abs((P0-Pvb)/P0)<0.05 && TrackClass>3");

	TCanvas *pCan=new TCanvas("pCan","Phi0 in G2P",800,600);
	gPad->SetRightMargin(0.05);
	
	track0->Draw("Phi0*180/3.14159:P0>>h2Phi0P0(175,0.,3.50,100,-10,90)",HRSCut,"prof");
	h2Phi0P0->SetMarkerStyle(20);
	h2Phi0P0->SetMarkerColor(1);
	h2Phi0P0->SetTitle(Form("#phi_{0} Vs P_{0} for HRS=%.1f^{o}; P_{0} (GeV); #phi_{0} (deg)",LHRSAngle*180./3.14159));
	h2Phi0P0->GetYaxis()->SetRangeUser(0,90);
	h2Phi0P0->Draw("prof");
	h2Phi0P0->Fit(f2,"R","", 0.3,3.5);
	
	TPaveText *pt2 = new TPaveText(0.3,22,2.0,35,"br");
	pt2->SetBorderSize(0);
	pt2->SetFillColor(0);
	text=pt2->AddText("f(x)=p0/x+p1+p2*x+p3*x*x");
	text->SetTextColor(4);
	pt2->Draw();

	pCan->cd();
	pCan->SaveAs(Form("PhiVsP_HRS%.1f.gif",LHRSAngle*180./3.14159));	
 }
 
void Theta0VsPvb(char *filename="")
{
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetStatH(0.2);
	gStyle->SetStatW(0.17);
	gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
	

	if(strlen(filename)>5)
	{
		cout<<"trying to open root file "<<filename<<endl;
		TFile* _file0=TFile::Open(filename);
	}
	ReadConfig();
	
	TF1 *f1 = new TF1("f1","exp([0]+[1]*x)+[2]+[3]*x",0.2,3.3);
	f1->SetFillColor(0);
	f1->SetMarkerStyle(20);
	f1->SetMarkerColor(2);
	f1->SetLineColor(2);
	f1->SetLineWidth(3);
	
	TF1 *f2 = new TF1("f2","exp([0]+[1]*x)+[2]+[3]*x",0.2,3.3);
	f2->SetFillColor(0);
	f2->SetMarkerStyle(20);
	f2->SetMarkerColor(4);
	f2->SetLineColor(4);
	f2->SetLineWidth(3);

	TText *text=0;
	char HRSCut[255];
	sprintf(HRSCut,"abs((P0-Pvb)/P0)<0.05 && TrackClass>3");

	TCanvas *pCan=new TCanvas("pCan","Theta0 in G2P",800,600);
	gPad->SetRightMargin(0.05);
	
	track0->Draw("Theta0*180/3.14159:P0>>h2Theta0P0(350,0.,3.50,80,0,80)",HRSCut,"prof");
	h2Theta0P0->SetMarkerStyle(20);
	h2Theta0P0->SetMarkerColor(1);
	h2Theta0P0->SetTitle(Form("#theta_{0} Vs P_{0} for HRS=%.1f^{o}; P_{0} (GeV); #theta_{0} (deg)",LHRSAngle*180./3.14159));
	h2Theta0P0->GetYaxis()->SetRangeUser(0,50);
	h2Theta0P0->Draw("prof");
	h2Theta0P0->Fit(f2,"R","", 0.3,3.3);
	
	// TPaveText *pt2 = new TPaveText(1.5,0.35,3.3,0.45,"br");
	// pt2->SetBorderSize(0);
	// pt2->SetFillColor(0);
	// text=pt2->AddText("f(x)=exp(p0+p1*x)+p2+p3*x");
	// text->SetTextColor(4);
	// pt2->Draw();

	TF1 *fQ2P1 = new TF1("fQ2P1","acos(1.0-0.1/(4*1.1*x))*180/3.14159",0.2,3.3);
	fQ2P1->SetMarkerStyle(22);
	fQ2P1->SetMarkerColor(8);
	fQ2P1->SetLineColor(8);
	fQ2P1->SetLineWidth(3);
	TF1 *fQ2P01 = new TF1("fQ2P01","acos(1.0-0.01/(4*1.1*x))*180/3.14159",0.2,3.3);
	fQ2P01->SetMarkerStyle(24);
	fQ2P01->SetMarkerColor(6);
	fQ2P01->SetLineColor(6);
	fQ2P01->SetLineWidth(3);
	TF1 *fQ2P05 = new TF1("fQ2P05","acos(1.0-0.05/(4*1.1*x))*180/3.14159",0.2,3.3);
	fQ2P05->SetMarkerStyle(28);
	fQ2P05->SetMarkerColor(7);
	fQ2P05->SetLineColor(7);
	fQ2P05->SetLineWidth(3);
	
	fQ2P1->Draw("same");
	fQ2P01->Draw("same");
	fQ2P05->Draw("same");
	
   TLegend *leg = new TLegend(0.64,0.35,0.95,0.58,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("f2","f(x)=e^{p0+p1*x}+p2+p3*x","lpf");
   entry->SetFillColor(0);
   entry->SetLineColor(4);
   entry->SetMarkerColor(4);
   entry->SetMarkerStyle(20);
   entry=leg->AddEntry("fQ2P1","Q^{2}=0.1","lpf");
   entry->SetFillColor(0);
   entry->SetLineColor(8);
   entry->SetMarkerColor(8);
   entry->SetMarkerStyle(22);
   entry=leg->AddEntry("fQ2P01","Q^{2}=0.01","lpf");
   entry->SetFillColor(0);
   entry->SetLineColor(6);
   entry->SetMarkerColor(6);
   entry->SetMarkerStyle(24);
   entry=leg->AddEntry("fQ2P05","Q^{2}=0.05","lpf");
   entry->SetFillColor(0);
   entry->SetLineColor(7);
   entry->SetMarkerColor(7);
   entry->SetMarkerStyle(28);
   leg->Draw();
   
	pCan->cd();
	pCan->SaveAs(Form("ThetaVsP_HRS%.1f.gif",LHRSAngle*180./3.14159));	
 }
 
void Theta0VsPvb_newfcn(char *filename="")
{
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetStatH(0.2);
	gStyle->SetStatW(0.17);
	gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
	

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
	
	TF1 *f2 = new TF1("f2","[0]/x+[1]+[2]*x",0.2,3.3);
	f2->SetFillColor(0);
	f2->SetMarkerStyle(20);
	f2->SetMarkerColor(4);
	f2->SetLineColor(4);
	f2->SetLineWidth(3);

	TText *text=0;
	char HRSCut[255];
	sprintf(HRSCut,"abs((P0-Pvb)/P0)<0.05 && TrackClass>3");

	TCanvas *pCan=new TCanvas("pCan","Theta0 in G2P",800,600);
	
	gPad->SetRightMargin(0.05);
	
	track0->Draw("Theta0*180/3.14159:P0>>h2Theta0P0(175,0.,3.50,80,0,80)",HRSCut,"prof");
	h2Theta0P0->SetMarkerStyle(20);
	h2Theta0P0->SetMarkerColor(1);
	h2Theta0P0->SetTitle(Form("#theta_{0} Vs P_{0} for HRS=%.1f^{o}; P_{0} (GeV); #theta_{0} (deg)",LHRSAngle*180./3.14159));
	h2Theta0P0->GetYaxis()->SetRangeUser(0,50);
	h2Theta0P0->Draw("prof");
	h2Theta0P0->Fit(f2,"R","", 0.3,3.3);
	
	// TPaveText *pt2 = new TPaveText(1.5,0.35,3.3,0.45,"br");
	// pt2->SetBorderSize(0);
	// pt2->SetFillColor(0);
	// text=pt2->AddText("f(x)=p0/x+p1+p2*x");
	// text->SetTextColor(4);
	// pt2->Draw();

	TF1 *fQ2P1 = new TF1("fQ2P1","acos(1.0-0.1/(4*1.1*x))*180/3.14159",0.2,3.3);
	fQ2P1->SetMarkerStyle(22);
	fQ2P1->SetMarkerColor(8);
	fQ2P1->SetLineColor(8);
	fQ2P1->SetLineWidth(3);
	TF1 *fQ2P01 = new TF1("fQ2P01","acos(1.0-0.01/(4*1.1*x))*180/3.14159",0.2,3.3);
	fQ2P01->SetMarkerStyle(24);
	fQ2P01->SetMarkerColor(6);
	fQ2P01->SetLineColor(6);
	fQ2P01->SetLineWidth(3);
	TF1 *fQ2P05 = new TF1("fQ2P05","acos(1.0-0.05/(4*1.1*x))*180/3.14159",0.2,3.3);
	fQ2P05->SetMarkerStyle(28);
	fQ2P05->SetMarkerColor(7);
	fQ2P05->SetLineColor(7);
	fQ2P05->SetLineWidth(3);
	
	fQ2P1->Draw("same");
	fQ2P01->Draw("same");
	fQ2P05->Draw("same");
	
   TLegend *leg = new TLegend(0.64,0.35,0.95,0.58,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("f2","f(x)=p0/x+p1+p2*x","lpf");
   entry->SetFillColor(0);
   entry->SetLineColor(4);
   entry->SetMarkerColor(4);
   entry->SetMarkerStyle(20);
   entry=leg->AddEntry("fQ2P1","Q^{2}=0.1","lpf");
   entry->SetFillColor(0);
   entry->SetLineColor(8);
   entry->SetMarkerColor(8);
   entry->SetMarkerStyle(22);
   entry=leg->AddEntry("fQ2P01","Q^{2}=0.01","lpf");
   entry->SetFillColor(0);
   entry->SetLineColor(6);
   entry->SetMarkerColor(6);
   entry->SetMarkerStyle(24);
   entry=leg->AddEntry("fQ2P05","Q^{2}=0.05","lpf");
   entry->SetFillColor(0);
   entry->SetLineColor(7);
   entry->SetMarkerColor(7);
   entry->SetMarkerStyle(28);
   leg->Draw();
   
	pCan->cd();
	pCan->SaveAs(Form("ThetaVsP_HRS%.1f_newfcn.gif",LHRSAngle*180./3.14159));	
 }
 ////////////////////////////////////////////////////////////////
 
void Theta0_trVsPvb(char *filename="")
{
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(10);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetStatH(0.2);
	gStyle->SetStatW(0.17);
	//gStyle->SetStatStyle(4000);
	gStyle->SetStatX(0.43);
	gStyle->SetStatX(0.95);
	gStyle->SetStatY(0.65);
	

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
	
	TF1 *f2 = new TF1("f2","[0]/x+[1]+[2]*x",0.2,3.3);
	f2->SetFillColor(0);
	f2->SetMarkerStyle(20);
	f2->SetMarkerColor(4);
	f2->SetLineColor(4);
	f2->SetLineWidth(3);

	TText *text=0;
	char HRSCut[255];
	sprintf(HRSCut,"abs((P0-Pvb)/P0)<0.05 && TrackClass>3");

	TCanvas *pCan=new TCanvas("pCan","Theta0_tr in G2P",800,600);
	gPad->SetRightMargin(0.05);
	
	track0->Draw("Theta0_tr*180/3.14159:P0>>h2Theta0_trP0(175,0.,3.50,80,0,80)",HRSCut,"prof");
	h2Theta0_trP0->SetMarkerStyle(20);
	h2Theta0_trP0->SetMarkerColor(1);
	h2Theta0_trP0->SetTitle(Form("#theta_{0}^{tr} Vs P_{0} for HRS=%.1f^{o}; P_{0} (GeV); #theta_{0}^{tr} (deg)",LHRSAngle*180./3.14159));
	h2Theta0_trP0->GetYaxis()->SetRangeUser(-45.0,5.0);
	h2Theta0_trP0->Draw("prof");
	h2Theta0_trP0->Fit(f2,"R","", 0.3,3.3);
	
	TPaveText *pt2 = new TPaveText(1.5,-28,3.3,-22,"br");
	pt2->SetBorderSize(0);
	pt2->SetFillColor(0);
	text=pt2->AddText("f(x)=p0/x+p1+p2*x");
	text->SetTextColor(4);
	pt2->Draw();

	pCan->cd();
	pCan->SaveAs(Form("Theta0_trVsP_HRS%.1f.gif",LHRSAngle*180./3.14159));	
 }

 void ThetaVsPvb_deg(const char *filename="")
 {
	 //plot them all
	 Theta0VsPvb(filename);
	 Theta0VsPvb_newfcn(filename);
	 Phi0VsPvb(filename);
	 Theta0_trVsPvb(filename);
 }
 
