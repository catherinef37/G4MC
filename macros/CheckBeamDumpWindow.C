
//Declaration of leaves types tree
   Int_t           Run=-1;
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


CheckBeamDumpWindow(char *filename="")
{
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);

	if(strlen(filename)>5)
	{
		cout<<"trying to open root file "<<filename<<endl;
		TFile* _file0=TFile::Open(filename);
	}
	ReadConfig();

	char EntrancePlaneCut[1024],ExitPlaneCut[1024];
	sprintf(EntrancePlaneCut,"abs(StepZ-(%.1f+640))<2 && TrackClass>3",TargetZOffset);
	sprintf(ExitPlaneCut,"abs(StepZ-(%.1f+790))<2 && TrackClass>3",TargetZOffset);
	
	TCanvas *pCan=new TCanvas("pCan","Check Beam Dump tutnnel size",600,800);
	pCan->Divide(1,2);
	pCan->cd(1);
	
	track0->Draw("StepY:StepX>>h2Entrance",EntrancePlaneCut,"");
	h2Entrance->SetMarkerStyle(4);
	double EntranceMax=0.5* h2Entrance->GetMaximum();
	h2Entrance->SetMaximum(EntranceMax);
	h2Entrance->SetTitle(Form("HRS=%.2fdeg, Z=640+/-2 mm; X (mm); Y(mm)",LHRSAngle*180./3.14159));
	gPad->Update();
	h2Entrance->Draw("contz");
	
	pCan->cd(2);
	track0->Draw("StepY:StepX>>h2Exit",ExitPlaneCut,"");
	h2Exit->SetMarkerStyle(4);
	double ExitMax=0.5* h2Exit->GetMaximum();
	h2Exit->SetMaximum(ExitMax);
	h2Exit->SetTitle(Form("HRS=%.2fdeg, Z=790+/-1 mm; X (mm); Y(mm)",LHRSAngle*180./3.14159));
	h2Exit->Draw("contz");
	gPad->Update();
	
	pCan->cd();
	if(bIsCombinedTree)
		pCan->SaveAs(Form("BeamDump_HRS%.2f_EAll.gif",LHRSAngle*180./3.14159));	
	else
		pCan->SaveAs(Form("BeamDump_HRS%.2f_E%.1f.gif",LHRSAngle*180./3.14159,Beam));
 }
 
 
//input:x1<X_Entrance<x2   y1<Y_Entrance<y2   x3<X_Exit<x4   y3<Y_Exit<y4
void DrawWindow(double x1,double x2,double y1,double y2,double x3,double x4,double y3,double y4,
	char *filename="")
{
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);

	if(strlen(filename)>5)
	{
		cout<<"trying to open root file "<<filename<<endl;
		TFile* _file0=TFile::Open(filename);
	}
	ReadConfig();

	TCanvas *pCan=new TCanvas("pCan","Check Beam Dump tutnnel size",600,800);
	pCan->Divide(1,2);
	pCan->cd(1);
	char EntrancePlaneCut[1024],ExitPlaneCut[1024];
	char EntranceWindowCut[1024],ExitWindowCut[1024];
//	sprintf(EntrancePlaneCut,"abs(StepZ-(%.1f+640))<2 && TrackClass>3",TargetZOffset);
//	sprintf(ExitPlaneCut,"abs(StepZ-(%.1f+790))<2 && TrackClass>3",TargetZOffset);

	sprintf(EntrancePlaneCut,"abs(StepZ-(%.1f+640))<2 && TrackClass>3 && P0>1.0",TargetZOffset);
	sprintf(ExitPlaneCut,"abs(StepZ-(%.1f+790))<2 && TrackClass>3 && P0>1.0",TargetZOffset);
	sprintf(EntranceWindowCut,"(%s) && StepX>%.1f && StepX<%.1f && StepY>%.1f && StepY<%.1f",
		EntrancePlaneCut,x1,x2,y1,y2);
	sprintf(ExitWindowCut,"(%s) && StepX>%.1f && StepX<%.1f && StepY>%.1f && StepY<%.1f",
		ExitPlaneCut,x3,x4,y3,y4);
		
	TLine* LEntranceV1=new TLine(x1,y1,x1,y2);
	TLine* LEntranceV2=new TLine(x2,y1,x2,y2);
	TLine* LEntranceH1=new TLine(x1,y1,x2,y1);
	TLine* LEntranceH2=new TLine(x1,y2,x2,y2);
	LEntranceV1->SetLineColor(2);LEntranceV1->SetLineWidth(3);
	LEntranceV2->SetLineColor(2);LEntranceV2->SetLineWidth(3);
	LEntranceH1->SetLineColor(2);LEntranceH1->SetLineWidth(3);
	LEntranceH2->SetLineColor(2);LEntranceH2->SetLineWidth(3);
	
	TLine* LExitV1=new TLine(x3,y3,x3,y4);
	TLine* LExitV2=new TLine(x4,y3,x4,y4);
	TLine* LExitH1=new TLine(x3,y3,x4,y3);
	TLine* LExitH2=new TLine(x3,y4,x4,y4);
	LExitV1->SetLineColor(2);LExitV1->SetLineWidth(3);
	LExitV2->SetLineColor(2);LExitV2->SetLineWidth(3);
	LExitH1->SetLineColor(2);LExitH1->SetLineWidth(3);
	LExitH2->SetLineColor(2);LExitH2->SetLineWidth(3);
	
	track0->Draw("StepY:StepX>>h2Entrance",EntrancePlaneCut,"");
	h2Entrance->SetMarkerStyle(4);
	double EntranceMax=0.5* h2Entrance->GetMaximum();
	h2Entrance->SetMaximum(EntranceMax);
	h2Entrance->SetTitle(Form("HRS=%.2fdeg, Entrance Plane; X (mm); Y(mm)",LHRSAngle*180./3.14159));
	track0->Draw("StepY:StepX>>h2EntranceWin",EntranceWindowCut,"same");
	//double NEntrance,NEntranceCut;
	//NEntrance=h2Entrance->GetEntries();
	//NEntranceCut=h2EntranceWin->GetEntries();
	
	//redraw
	h2Entrance->GetListOfFunctions()->Add(LEntranceV1);
	h2Entrance->GetListOfFunctions()->Add(LEntranceV2);
	h2Entrance->GetListOfFunctions()->Add(LEntranceH1);
	h2Entrance->GetListOfFunctions()->Add(LEntranceH2);
	h2Entrance->Draw("contz");
	LEntranceV1->Draw("same");	
	LEntranceV2->Draw("same");	
	LEntranceH1->Draw("same");	
	LEntranceH2->Draw("same");	
	
	//show the stats
	//TPaveStats *pPSEn = new TPaveStats(0.10,0.60,0.43,0.90,"brNDC");
	TPaveStats *pPSEn = new TPaveStats(0.30,0.35,0.60,0.60,"brNDC");
	pPSEn->SetName("mywords");
	pPSEn->SetBorderSize(2);
	pPSEn->SetFillColor(0);
	pPSEn->SetTextAlign(12);
    //pPSEn->SetOptStat(110);
	pPSEn->AddText(Form("Raw %.0f",h2Entrance->GetEntries()));
	pPSEn->AddText(Form("Cut %.0f",h2EntranceWin->GetEntries()));
	pPSEn->AddText(Form("%.1f < X < %.1f",x1,x2));
	pPSEn->AddText(Form("%.1f < Y < %.1f",y1,y2));
	//pPSEn->SetFillStyle(4000);
	pPSEn->Draw();
	//h2Entrance->GetListOfFunctions()->Add(pPSEn);
	//pPSEn->SetParent(h2Entrance->GetListOfFunctions());
	
	
	pCan->cd(2);
	track0->Draw("StepY:StepX>>h2Exit",ExitPlaneCut,"");
	h2Exit->SetMarkerStyle(4);
	double ExitMax=0.5* h2Exit->GetMaximum();
	h2Exit->SetMaximum(ExitMax);
	h2Exit->SetTitle(Form("HRS=%.2fdeg, Exit Plane; X (mm); Y(mm)",LHRSAngle*180./3.14159));
	track0->Draw("StepY:StepX>>h2ExitWin",ExitWindowCut,"same");
	
	//redraw
	h2Exit->GetListOfFunctions()->Add(LExitV1);
	h2Exit->GetListOfFunctions()->Add(LExitV2);
	h2Exit->GetListOfFunctions()->Add(LExitH1);
	h2Exit->GetListOfFunctions()->Add(LExitH2);
	h2Exit->Draw("contz");
	LExitV1->Draw("same");	
	LExitV2->Draw("same");	
	LExitH1->Draw("same");	
	LExitH2->Draw("same");	
	//show the stats
	TPaveStats *pPSEx = new TPaveStats(0.30,0.35,0.60,0.60,"brNDC");
	pPSEx->SetName("mywords");
	pPSEx->SetBorderSize(2);
	pPSEx->SetFillColor(0);
	pPSEx->SetTextAlign(12);
    //pPSEx->SetOptStat(110);
	pPSEx->AddText(Form("Raw %7.0f",h2Exit->GetEntries()));
	pPSEx->AddText(Form("Cut %7.0f",h2ExitWin->GetEntries()));
	pPSEx->AddText(Form("%.1f < X < %.1f",x3,x4));
	pPSEx->AddText(Form("%.1f < Y < %.1f",y3,y4));
	//pPSEx->SetFillStyle(4000);
	pPSEx->Draw();
	//h2Exit->GetListOfFunctions()->Add(pPSEx);
	//pPSEx->SetParent(h2Exit->GetListOfFunctions());
	
	
	pCan->cd();
	if(bIsCombinedTree)
		pCan->SaveAs(Form("BeamDump_HRS%.2f_EAll_R%.1f_result.gif",LHRSAngle*180./3.14159,HelmCurrentRatio));	
	else
		pCan->SaveAs(Form("BeamDump_HRS%.2f_E%.1f_R%.1f_result.gif",LHRSAngle*180./3.14159,Beam,HelmCurrentRatio));
 }
 
<<<<<<< .mine
 
void DoPlot(char* filename="")
{
  if(Run<0) ReadConfig();
  if(LHRSAngle*180/3.14159<12)
  {
    DrawWindow(46,87,-43,50,58,106,-53,58,filename);
  }
  else
  {
    DrawWindow(112,169,-54,60,140,208,-66,70,filename);
  }
}

=======
 
void DoPlot(char* filename="")
{
  if(Run<0) ReadConfig();
  if(LHRSAngle*180/3.14159<12)
  {
    DrawWindow(46,87,-43,50,58,106,-53,58,filename);
  }
  else
  {
    DrawWindow(112,169,-54,60,140,208,-66,70,filename);
  }
}
>>>>>>> .r187
