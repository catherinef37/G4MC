//this script is used to plot the acceptance for one perticular 
//sensertive detector, one has to provide the sdid, which can be
//found in the config tree. simply type config->Show(0) and then all 
//sdid will be shown
 
void PlotAcc_PhiTheta(const char* key="", int sdid=9)
{
  gStyle->SetOptStat(0);
  //set the draw-option "text" format
  gStyle->SetPaintTextFormat(".0f");

  TCanvas *c1=new TCanvas("c1","",800,600);
  c1->SetRightMargin(0.15);
  c1->cd(0);
  D->Draw("T_Phi[1]*57.3:T_Theta[1]*57.3>>hD(25,5,30,40,-60,60)","","");
  D->Draw("T_Phi[1]*57.3:T_Theta[1]*57.3>>hN(25,5,30,40,-60,60)",Form("SD_Tid==2 && SD_Id==%d",sdid),"");
  TH2D *hD=(TH2D*) gROOT->FindObject("hD");
  TH2D *hN=(TH2D*) gROOT->FindObject("hN");
  TH2D* h2N = hN->Clone("h2N");
  h2N->Divide(hD);
  h2N->SetTitle("RTPC+ABS Acc in \%; #theta (deg) ;#phi (deg) ");
  h2N->Scale(100.0);
  h2N->Draw("colz text");

  char strName[100];
  sprintf(strName,"SBSAcc_PhiTheta%s.png",key);
  c1->SaveAs(strName);
  sprintf(strName,"SBSAcc_PhiTheta%s.C",key);
  c1->SaveAs(strName);
}

void PlotAcc_PTheta(const char* key="", int sdid=9)
{
  gStyle->SetOptStat(0);
  //set the draw-option "text" format
  gStyle->SetPaintTextFormat(".0f");
  

  TCanvas *c2=new TCanvas("c2","",800,600);
  c2->cd(0);
  c2->SetRightMargin(0.15);
 
  D->Draw("T_P[1]:T_Theta[1]*57.3>>hDpt(25,5,30,48,0,12)","","");
  D->Draw("T_P[1]:T_Theta[1]*57.3>>hNpt(25,5,30,48,0,12)",Form("SD_Tid==2 && SD_Id==%d",sdid),"");
  TH2D *hD=(TH2D*) gROOT->FindObject("hDpt");
  TH2D *hN=(TH2D*) gROOT->FindObject("hNpt");
  TH2D* h2N = hN->Clone("h2Npt");
  h2N->Divide(hD);
  h2N->SetTitle("RTPC+ABS Acc in \%; #theta (deg) ; P (GeV/c) ");
  h2N->Scale(100.0);
  h2N->Draw("colz text");

  char strName[100];
  sprintf(strName,"SBSAcc_PTheta%s.png",key);
  c2->SaveAs(strName);
  sprintf(strName,"SBSAcc_PTheta%s.C",key);
  c2->SaveAs(strName);
}


void PlotAcc(const char* key="",int sdid=9)
{
  PlotAcc_PhiTheta(key,sdid);
  PlotAcc_PTheta(key,sdid);
}
