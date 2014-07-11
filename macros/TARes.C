{
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  
  TCanvas *c1=new TCanvas("c1","",1200,800);
  c1->Divide(4,2);
 
  D->AddFriend("track0");
  c1->cd(3);
  config->Draw("Beam");
  double Beam=htemp->GetMean();
  
  c1->cd(1);
  track0->Draw("-Xvb_tr:-Yvb_tr >> h2Good","Pvb>0","contz");
  h2Good->SetTitle(" Hits on VD; horizontal (mm) ; vertical (mm) ");
  h2Good->Draw("contz");
  
  c1->cd(2);
  track0->Draw("Theta0*57.3 >> ht(40,60,80)","Pvb>0","");

  c1->cd(3);
  track0->Draw("P0 >> hp(100,0.1,0.6)","Pvb>0","");
  hp->SetTitle("P0 (black) and Pvb(red); P0 or Pvb (GeV)");
  track0->Draw("Pvb >> hpvb","Pvb>0","same");
  hpvb->SetLineColor(2);
  
  double dECut1=8.5,dECut2=20.0;
  if(Beam>1.7)
    {
      dECut1=7.5,dECut2=39.0;
    }
  else if(Beam>2.2)
    {
      dECut1=7.5,dECut2=39.0;
    }

  TCut ADC1Cut = Form("TA1_Edep>%.1f && TA1_N>0",dECut1);
  TCut ADC2Cut = Form("TA2_Edep>%.1f && TA2_N>0",dECut2);
 
  c1->cd(4);
  D->Draw("TA1_Edep:TA1_P>>hdEP1(100,0.1,0.6,80,0,40)","TA1_N>0","colz");
  hdEP1->SetTitle(Form("dE plane ADC Cut=%.1fMeV; TA1_P (GeV) ; TA1_Edep (MeV)",dECut1));
  TF1 *f1=new TF1("f1",Form("0*x+%.1f",dECut1),0,1);
  f1->SetLineWidth(2);
  hdEP1->Draw("colz");
  f1->Draw("same");

  c1->cd(5);
  D->Draw("TA2_Edep:TA1_P>>hdEP2(80,0.2,0.6,100,0,100)","TA2_N>0","colz");
  hdEP2->SetTitle(Form("E plane ADC Cut=%.1fMeV; TA1_P (GeV) ; TA2_Edep (MeV)",dECut2));
  TF1 *f2=new TF1("f2",Form("0*x+%.1f",dECut2),0,1);
  f2->SetLineWidth(2);
  hdEP2->Draw("colz");
  f2->Draw("same");


  c1->cd(6);
  D->Draw("TA1_P>>hpin(70,0.2,0.55)","TA1_N>0");
  hpin->SetTitle("#splitline{TA1_P: All(black)}{Penetrated(red) and ADCCut(blue)}; TA1_P (GeV)");
  D->Draw("TA1_P>>hpin_coin","TA1_Pout>0","same");
  hpin_coin->SetLineColor(2);
  D->Draw("TA1_P>>hpin_cut",ADC1Cut && ADC2Cut,"same");
  hpin_cut->SetLineColor(4);


  c1->cd(7);
  D->Draw("TA1_Edep>>hE1","TA1_N>0");
  double pEtot1=hE1->GetMean()*hE1->GetEntries();
  hE1->SetTitle(Form("Integrated Edep=%.1f MeV; TA1_Edep (MeV)",pEtot1));


  c1->cd(8);
  D->Draw("TA2_Edep>>hE2","TA2_N>0");
  double pEtot2=hE2->GetMean()*hE2->GetEntries();
  hE2->SetTitle(Form("Integrated Edep=%.1f MeV; TA2_Edep (MeV)",pEtot2));
  
  
  c1->SaveAs(Form("TARes_E%.3f.png",Beam)); 
}

