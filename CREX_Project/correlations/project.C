project(){
  TFile* infile1 = new TFile("theta.root", "READ");
  infile1->cd();

  TH2F* vb = new TH2F("vb" , "Virtual Boundary", 200, -75, 75, 200, -75, 75);
  TH2F* fp = new TH2F("fp" , "Focal Plane"     , 200, -75, 75, 200, -75, 75);
  TH2F* pr = new TH2F("pr" , "Projection 1 m beyond Focal Plane", 200, -75, 75, 200, -75, 75);
  TH2F* th = new TH2F("th" , "Theta and Y_{fp} Correlation", 200,   3,  7, 200, -75, 75);
  TH2F* pvb= new TH2F("pvb", "Virtual Boundary Momentum", 200, -2, 2, 200, -2, 2);
  vb->GetXaxis()->SetTitle("X_vb");
  vb->GetYaxis()->SetTitle("Y_vb");
  fp->GetXaxis()->SetTitle("X_fp");
  fp->GetYaxis()->SetTitle("Y_fp");
  pr->GetXaxis()->SetTitle("X_pr");
  pr->GetYaxis()->SetTitle("Y_pr");
  th->GetXaxis()->SetTitle("#theta_0");
  th->GetYaxis()->SetTitle("Y_fp");
  pvb->GetXaxis()->SetTitle("p_X_vb");
  pvb->GetYaxis()->SetTitle("p_Y_vb");
  TCanvas* can = new TCanvas("can", "can", 20, 20, 1100, 700);
  can->Divide(3, 2);
  can->cd(1);
  track0->Draw("Yvb_tr:Xvb_tr>>vb", "rate_208Pb", "COLZ");
  can->cd(2);
  track0->Draw("Yfp_tr:Xfp_tr>>fp", "rate_208Pb", "COLZ");
  can->cd(3);
  track0->Draw("Ypr_tr:Xpr_tr>>pr", "rate_208Pb", "COLZ");
  can->cd(4);
  track0->Draw("Thetavb*Pvb : Phivb*Pvb >> pvb", "rate_208Pb", "COLZ");
  can->cd(6);
  track0->Draw("Yfp_tr:Theta0*TMath::RadToDeg()>>th", "", "COLZ");
  th->ProfileX()->Draw("SAME");
  th->ProfileX()->SetMarkerColor(2);
  th->ProfileX()->SetLineColor(2);
  TF1* f1 = new TF1("f1", "[0] + [1] * x", 3, 7);
  th->ProfileX()->Fit("f1", "", "", 3., 7.);

  TFile* infile2 = new TFile("p.root", "READ");
  infile2->cd();

  TH2F* vb2 = new TH2F("vb2", "Virtual Boundary", 200, -75, 75, 200, -75, 75);
  TH2F* fp2 = new TH2F("fp2", "Focal Plane"     , 200, -75, 75, 200, -75, 75);
  TH2F* pr2 = new TH2F("pr2", "Projection 1 m beyond Focal Plane", 200, -75, 75, 200, -75, 75);
  TH2F* th2 = new TH2F("p2" , "P_{0} and X_{fp} Correlation", 200, 1., 1.2, 200, -75, 75);
  vb2->GetXaxis()->SetTitle("X_vb");
  vb2->GetYaxis()->SetTitle("Y_vb");
  fp2->GetXaxis()->SetTitle("X_fp");
  fp2->GetYaxis()->SetTitle("Y_fp");
  pr2->GetXaxis()->SetTitle("X_pr");
  pr2->GetYaxis()->SetTitle("Y_pr");
  p2->GetXaxis()->SetTitle("P_0");
  p2->GetYaxis()->SetTitle("X_fp");
  TCanvas* can2 = new TCanvas("can2", "can2", 40, 40, 1100, 700);
  can2->Divide(3, 2);
  can2->cd(1);
  track0->Draw("Yvb_tr:Xvb_tr>>vb2", "rate_208Pb", "COLZ");
  can2->cd(2);
  track0->Draw("Yfp_tr:Xfp_tr>>fp2", "rate_208Pb", "COLZ");
  can2->cd(3);
  track0->Draw("Ypr_tr:Xpr_tr>>pr2", "rate_208Pb", "COLZ");
  can2->cd(6);
  track0->Draw("Xfp_tr:P0>>p2", "", "COLZ");
  th2->ProfileX()->Draw("SAME");
  th2->ProfileX()->SetMarkerColor(2);
  th2->ProfileX()->SetLineColor(2);
  TF1* f2 = new TF1("f2", "[0] + [1] * x", 2, 2.4);
  th2->ProfileX()->Fit("f2", "", "", 2., 2.25);

  can->SaveAs("theta.png");
  can2->SaveAs("p.png");

}
