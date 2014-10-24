void XY_3rdArm()
{

	gStyle->SetStatW(0.28);
	gStyle->SetStatH(0.06);
	gStyle->SetOptStat(110);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	track1->AddFriend("track0");
	TCanvas* pCan=new TCanvas("pCan","",1024,768);
	pCan->Divide(3,3);
	//check bigbite hits

	char EndPlaneCut[1024],DetectorCut[1024];
	sprintf(EndPlaneCut,"abs(StepZ-637.5)<40.0 && abs(StepX+172)<60 && abs(StepY+90)<60");
	sprintf(DetectorCut,"TrackClass>0 && abs(StepZ-637.5)<40.0 && abs(StepX+172)<60 && abs(StepY+90)<60");

	pCan->cd(1);
	//track0->Draw("StepY:StepX","","colz");
	//track0->Draw("StepY:StepX>>h2YX_pr","abs(StepX-3750)<400","colz");
	track0->Draw("Yvb:Xvb>>h2YX_pr","abs(Xvb-3700)<350","colz");
	h2YX_pr->SetTitle("Pr: Yvb:Xvb @ 3rd arm; Xvb;Yvb");
	h2YX_pr->Draw("colz");

	pCan->cd(2);
	track0->Draw("Theta0*57.3:Phi0*57.3>>h2ThetaPhi_pr","","colz");
	h2ThetaPhi_pr->SetTitle("Pr: Theta0:Phi0 ;Phi0 (deg);Theta0 (deg)");
	h2ThetaPhi_pr->Draw("colz");

	pCan->cd(3);
	track0->Draw("P0>>hP0");
	hP0->SetTitle("Pr: P0(black) and Pvb(red); P0 or Pvb (GeV)");
	track0->Draw("Pvb>>hPvb","","same");
	hPvb->SetLineColor(2);  
	hP0->Draw();
	hPvb->Draw("same");

	pCan->cd(4);
	//track1->Draw("637.5*tan(Theta0)-StepY >> hDiff",
	//	Form("%s && abs(640*tan(Theta0)-StepY-280)<120",EndPlaneCut),"");
	//hDiff->SetTitle("e: Displacement due to field; Ydiff (mm)");
	TH2D *h2=0;
	for(int ThetaSlice=0; ThetaSlice<13;ThetaSlice+=3)
	{
		pCan->cd(4+ThetaSlice/3);
		track0->Draw(Form("Xvb_tr:Yvb_tr>>h2XY_T%d",ThetaSlice),
		Form("abs(Theta0*180./3.14159-62.5-%d)<1.5 && Pvb>0",ThetaSlice),"colz");
		h2=(TH2D*)gROOT->FindObject(Form("h2XY_T%d",ThetaSlice));
		h2->SetTitle(Form("Yvb_tr:Xvb_tr  %.0f<Theta0<%.0f; Yvb_tr (mm); Xvb_tr (mm)",
			61.0+ThetaSlice,64.0+ThetaSlice));
		h2->Draw("colz");
	}
	
	// pCan->cd(5);
	// track1->Draw("StepY>>hStepY",EndPlaneCut,"");
	// hStepY->SetTitle("e: StepY: All(black) and Pass HRS(red);StepY (mm)");
	// track1->Draw("StepY>>hStepY_HRS",DetectorCut,"same");
	// hStepY_HRS->SetLineColor(2);  

	// pCan->cd(6);
	// track1->Draw("StepX>>hStepX",EndPlaneCut,"");
	// hStepX->SetTitle("e: StepX: All(black) and Pass HRS(red); StepX (mm)");
	// track1->Draw("StepX>>hStepX_HRS",DetectorCut,"same");
	// hStepX_HRS->SetLineColor(2);  

	// pCan->cd(7);
	// track1->Draw("Theta0*57.3:Phi0*57.3>>h2Theta0Phi0",EndPlaneCut,"");
	// h2Theta0Phi0->SetTitle("e: Theta0:Phi0 ; Phi0 (deg);Theta0 (deg)");
	// h2Theta0Phi0->Draw("colz");

	// pCan->cd(8);
	// track1->Draw("Yfp_tr:Xfp_tr>>h2YvbXvb",DetectorCut,"");
	// h2YvbXvb->SetTitle("e: Yfp:Xfp @ focus plane;Xfp (mm);Yfp (mm)");
	// h2YvbXvb->Draw("colz");
	
	pCan->cd(9);
	track0->Draw("Xvb_tr:Yvb_tr>>h2XYAll","Pvb>0","colz");
	h2XYAll->SetTitle("Xvb_tr:Yvb_tr  All Theta; Yvb_tr (mm):Xvb_tr (mm)");
	
	
	pCan->SaveAs("ThirdArm_XY.png");

}

