
void SieveStat(char *key="L", int left=1)
{
  TCanvas *c1=new TCanvas("c1","Sieve Rates",1300,800);
  TCanvas *c2=new TCanvas("c2","Sieve Rates",1100,600);

  //this is the perfect bin size, but not good in binning
  //double xbins[8]={-24.4,-18.3,-12.2,-6.1,0.0,4.8,9.6,14.4};
  //this is a work arround
  double xbinsL[8]={-24.0,-18.0,-12.0,-6.0,0.0,4.5,9.0,13.5};
  double xbinsR[8]={-13.5,-9.0,-4.5,0.0,6.0,12.0,18.0,24.0};
  double *xbins=xbinsL;
  if(left==0) xbins=xbinsR;
  //vertical span=13.3mm x 7=93.1mm, from -43.6 to 49.6
  double ybins[8]={-43.6,-30.3,-17.0,-3.7,9.6,22.9,36.2,49.5};

  TH2F *h2;
  h2=new TH2F("h0",Form("Sieve hole XS weighted yield: %s; -Yvb_tr (mm); -Xvb_tr (mm)",key),
	      7,xbins[0]-1.0,xbins[7]+0.5,7,ybins[0]-2.0,ybins[7]+1.0);
  //		    7,-25.5,13.5,7,-45.9,49.5);
  //		    25,-24.0,13.5,21,-43.6,49.5);
  h2->SetMinimum(5);
  h2->SetMarkerSize(1);
  h2->GetXaxis()->SetTitleOffset(0.6);
  h2->GetYaxis()->SetTitleOffset(0.6);

  TLine *pLH[8];
  TLine *pLV[8];
  bool bDrawGrid=false;
  if (bDrawGrid)
    {
      for(int i=0;i<=7;i++)
	{
	  pLH[i]=new TLine(xbins[0],ybins[i],xbins[7],ybins[i]);
	  pLV[i]=new TLine(xbins[i],ybins[0],xbins[i],ybins[7]);
	  pLH[i]->SetLineStyle(3);
	  pLV[i]->SetLineStyle(3);
	}
    }

  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  
  double topmargin=0.092,bottommargin=0.08,leftmargin=0.08,rightmargin=0.02;
  double padheight=1.0-topmargin-bottommargin-0.03;
  double padwidth=1.0-leftmargin-rightmargin-0.03;
  c1->SetTopMargin(topmargin);
  c1->SetLeftMargin(leftmargin);
  c1->SetRightMargin(rightmargin);
  c1->SetBottomMargin(bottommargin);
  c1->SetFrameBorderMode(0);
  
  c1->cd();  h2->Draw();


  TH1F *h1Thetavb_tr[49],*h1Phivb_tr[49];
  TH1F *h1Theta0_tr[49],*h1Phi0_tr[49];
  TH1F *h1Theta0Eff[49],*h1Theta0[49];
  TH1F *h1P0[49],*h1Pvb[49];

  c2->Divide(4,2,0.001,0.001);
  system("mkdir -p stat");

  TPaveText *pPad[49];
  TText *text;
  double XX=(xbins[7]-xbins[0])*1.0;
  double YY=(ybins[7]-ybins[0])*1.0;
  for(int i=0;i<7;i++)
    {
      for(int j=0;j<7;j++)
	{
          int idx=i*7+j;

	  TCut thisCut=Form("(TrackClass>4 && -Xvb_tr>%.3f && -Xvb_tr<=%.3f && -Yvb_tr>%.3f && -Yvb_tr<=%.3f)*ElasXS/1000",ybins[j],ybins[j+1],xbins[i],xbins[i+1]);

	  h1Theta0[idx]=new TH1F(Form("Theta0_%d%d",i,j),Form("Theta0: hole_%d%d",i,j),200,0.0,0.20);
	  h1Theta0Eff[idx]=new TH1F(Form("Theta0Eff_%d%d",i,j),"Theta0Eff",200,0.0,0.20);
	  h1Theta0_tr[idx]=new TH1F(Form("Theta0_tr_%d%d",i,j),"Theta0_tr",400,-0.30,0.10);
	  h1Thetavb_tr[idx]=new TH1F(Form("Thetavb_tr_%d%d",i,j),"Thetavb_tr",100,-0.05,0.05);

	  h1Phi0_tr[idx]=new TH1F(Form("Phi0_tr_%d%d",i,j),"Phi0_tr",100,-0.05,0.05);
	  h1Phivb_tr[idx]=new TH1F(Form("Phivb_tr_%d%d",i,j),"Phivb_tr",100,-0.05,0.05);

	  h1P0[idx]=new TH1F(Form("P0_%d%d",i,j),"P0",6000,0.4,3.4);
	  h1Pvb[idx]=new TH1F(Form("Pvb_%d%d",i,j),"Pvb",6000,0.4,3.4);

	  c2->cd(1);
	  track0->Draw(Form("Theta0>>Theta0_%d%d",i,j),thisCut,"");
	  double pMeanTheta0=h1Theta0[idx]->GetMean();

	  c2->cd(2);
	  track0->Draw(Form("Theta0Eff>>Theta0Eff_%d%d",i,j),thisCut,"");
	  double pMeanTheta0Eff=h1Theta0Eff[idx]->GetMean();

	  c2->cd(3);
	  track0->Draw(Form("Theta0_tr>>Theta0_tr_%d%d",i,j),thisCut,"");
	  double pMeanTheta0_tr=h1Theta0_tr[idx]->GetMean();

	  c2->cd(4);
	  track0->Draw(Form("Thetavb_tr>>Thetavb_tr_%d%d",i,j),thisCut,"");
	  double pMeanThetavb_tr=h1Thetavb_tr[idx]->GetMean();
	  
	  c2->cd(5);
	  track0->Draw(Form("Phi0_tr>>Phi0_tr_%d%d",i,j),thisCut,"");
	  double pMeanPhi0_tr=h1Phi0_tr[idx]->GetMean();

	  c2->cd(6);
	  track0->Draw(Form("Phivb_tr>>Phivb_tr_%d%d",i,j),thisCut,"");
	  double pMeanPhivb_tr=h1Phivb_tr[idx]->GetMean();

	  c2->cd(7);
	  track0->Draw(Form("P0>>P0_%d%d",i,j),thisCut,"");
	  double pMeanP0=h1P0[idx]->GetMean();

	  c2->cd(8);
	  track0->Draw(Form("Pvb>>Pvb_%d%d",i,j),thisCut,"");
	  double pMeanPvb=h1Pvb[idx]->GetMean();

	  int pN=h1P0[idx]->Integral();

	  c2->cd();
	  c2->SaveAs(Form("stat/SieveStat_Hole%d%d_%s.png",i,j,key));

	  if(pN<1) continue;

	  c1->cd();
	  /*
	  pPad[idx]=new TPaveText( (xbins[i]-xbins[0])/XX*padwidth+leftmargin+0.02,
				   (ybins[j]-ybins[0])/YY*padheight+bottommargin+0.02,
				   (xbins[i+1]-xbins[0])/XX*padwidth+leftmargin+0.02,
				   (ybins[j+1]-ybins[0])/YY*padheight+bottommargin+0.02,
				   "brNDC");
	  */
	  pPad[idx]=new TPaveText( xbins[i],ybins[j],xbins[i+1],ybins[j+1],"br");
	  pPad[idx]->SetFillColor(0);
	  pPad[idx]->SetBorderSize(1);
	  text=pPad[idx]->AddText(Form("Entries=#color[2]{%5d}",pN));
	  text=pPad[idx]->AddText(Form("<#theta_{0}>=%.1f, <#theta_{0}^{eff}>=%.1f",
				       pMeanTheta0*57.3,pMeanTheta0Eff*57.3));
	  text=pPad[idx]->AddText(Form("<#theta_{0}^{tr}>=%.1f, <#theta_{vb}^{tr}>=%.1f",
				       pMeanTheta0_tr*57.3,pMeanThetavb_tr*57.3));
	  text=pPad[idx]->AddText(Form("<#phi_{0}^{tr}>=%.1f, <#phi_{vb}^{tr}>=%.1f",
				       pMeanPhi0_tr*57.3,pMeanPhivb_tr*57.3));
	  //text=pPad[idx]->AddText(Form("<P_{0}>=%.4f, <P_{vb}>=%.4f",pMeanP0,pMeanPvb));
	  text=pPad[idx]->AddText(Form("<P_{0}>=%.4f",pMeanP0));
	  text=pPad[idx]->AddText(Form("<P_{vb}>=%.4f",pMeanPvb));

	  pPad[idx]->Draw("same");
	}
    }

  c1->cd(0);
  c1->SaveAs(Form("SieveStat_%s.png",key));
  return;
}


void SieveRate_Core(char *key="L",int left=1)
{
  TCanvas *c1=new TCanvas("c0","Sieve Rates",800,600);

  TCut cut1 = "TrackClass>4 && abs(-Yvb_tr+5)<20 && abs(Xvb_tr)<50";
  TCut cut2 = "(TrackClass>4 && abs(-Yvb_tr+5)<20 && abs(Xvb_tr)<50)*ElasXS/1000"; //Weighted by cross section

  //this is the perfect bin size, but not good in binning
  //double xbins[8]={-24.4,-18.3,-12.2,-6.1,0.0,4.8,9.6,14.4};
  //this is a work arround
  double xbinsL[8]={-24.0,-18.0,-12.0,-6.0,0.0,4.5,9.0,13.5};
  double xbinsR[8]={-13.5,-9.0,-4.5,0.0,6.0,12.0,18.0,24.0};
  double *xbins=xbinsL;
  if(left==0) xbins=xbinsR;

  //vertical span=13.3mm x 7=93.1mm, from -43.6 to 49.6
  double ybins[8]={-43.6,-30.3,-17.0,-3.7,9.6,22.9,36.2,49.5};


  TH2F *h2=new TH2F("h2",Form("Sieve hole XS weighted yield: %s; -Yvb_tr (mm); -Xvb_tr (mm)",key),
		    7,-25.5,13.5,7,-45.9,49.5);
  //		    25,-24.0,13.5,21,-43.6,49.5);
  h2->SetMinimum(5);
  h2->SetMarkerSize(1);
  h2->GetXaxis()->SetTitleOffset(0.6);
  h2->GetYaxis()->SetTitleOffset(0.6);

  TLine *pLH[8];
  TLine *pLV[8];
  bool bDrawGrid=false;
  if (bDrawGrid)
    {
      for(int i=0;i<=7;i++)
	{
	  pLH[i]=new TLine(xbins[0],ybins[i],xbins[7],ybins[i]);
	  pLV[i]=new TLine(xbins[i],ybins[0],xbins[i],ybins[7]);
	  pLH[i]->SetLineStyle(3);
	  pLV[i]->SetLineStyle(3);
	}
    }

  //TPad TPad(const char* name, const char* title, Double_t xlow, Double_t ylow, Double_t xup, 
  //          Double_t yup, Color_t color = -1, Short_t bordersize = -1, Short_t bordermode = -2)
 
  //c1->Range(-28.6875,-55.2375,18.1875,61.1375);
  //c1->Range(xbins[0],ybins[0],xbins[7],ybins[7]);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  //c1->SetLogz();
  c1->SetTopMargin(0.0995);
  c1->SetLeftMargin(0.08);
  c1->SetRightMargin(0.08);
  c1->SetBottomMargin(0.08);
  c1->SetFrameBorderMode(0);
  h2->Draw();

  TH1F *pH1[49];
  TPad *pPad[49];
  double XX=(xbins[7]-xbins[0])*1.0;
  double YY=(ybins[7]-ybins[0])*1.0;
  for(int i=0;i<7;i++)
    {
      for(int j=0;j<7;j++)
	{
          int idx=i*7+j;
	  pPad[idx]=new TPad(Form("pad_%d%d",i+1,j+1),Form("Hole %d%d",i,j),
			   (xbins[i]-xbins[0])/XX*0.82+0.1,(ybins[j]-ybins[0])/YY*0.8+0.1,
			   (xbins[i+1]-xbins[0])/XX*0.82+0.1,(ybins[j+1]-ybins[0])/YY*0.8+0.1,
			   0,1,2);
	  pPad[idx]->Draw("same");
	  pPad[idx]->cd();
	  pPad[idx]->SetLeftMargin(0);
	  pPad[idx]->SetRightMargin(0);
	  pPad[idx]->SetTopMargin(0);
	  pPad[idx]->SetBottomMargin(0);
	  pPad[idx]->SetFrameFillColor(0);
	  pPad[idx]->SetFrameBorderMode(0);
	  pPad[idx]->SetFrameBorderSize(0);

	  TCut thisCut=Form("(TrackClass>4 && -Xvb_tr>%.3f && -Xvb_tr<=%.3f && -Yvb_tr>%.3f && -Yvb_tr<=%.3f)*ElasXS/1000",ybins[j],ybins[j+1],xbins[i],xbins[i+1]);
	  pH1[idx]=new TH1F(Form("h_%d%d",i,j),"",//Form("Hole %d%d",i,j),
			    1,xbins[i],xbins[i+1]);
	  pH1[idx]->SetMarkerSize(14);
	  track0->Draw(Form("-Yvb_tr>>h_%d%d",i,j),thisCut,"text0");

	  c1->cd();
	}
    }

  if (bDrawGrid)
    {
      for(int i=0;i<=7;i++)
	{
	  pLH[i]->Draw("same");
	  pLV[i]->Draw("same");
	}
    }

  c1->cd(0);
  c1->SaveAs(Form("SieveRates_%s.png",key));
  return;
}


void SieveRate0(char *key="_")
{
  TCanvas *c1=new TCanvas("c1","Sieve Rates",1000,700);

  TCut cut1 = "TrackClass>4 && abs(-Yvb_tr+5)<20 && abs(Xvb_tr)<50";
  TCut cut2 = "(TrackClass>4 && abs(-Yvb_tr+5)<20 && abs(Xvb_tr)<50)*ElasXS/1000"; //Weighted by cross section

  //this is the perfect bin size, but not good in binning
  //double xbins[8]={-24.4,-18.3,-12.2,-6.1,0.0,4.8,9.6,14.4};
  //this is a work arround
  double xbins[8]={-24.0,-18.0,-12.0,-6.0,0.0,4.5,9.0,13.5};
  //vertical span=13.3mm x 7=93.1mm, from -43.6 to 49.6
  double ybins[8]={-43.6,-30.3,-17.0,-3.7,9.6,22.9,36.2,49.5};


  TH2F *h2=new TH2F("h2","XS weighted yield for sieve holes; -Yvb_tr (mm); -Xvb_tr (mm)",
		    25,-24,13.5,21,-43.6,49.5);
  h2->SetMinimum(5);
  h2->SetMarkerSize(1);
 
  TLine *pLH[8];
  TLine *pLV[8];
  bool bDrawGrid=true;
  if (bDrawGrid)
    {
      for(int i=0;i<=7;i++)
	{
	  pLH[i]=new TLine(xbins[0],ybins[i],xbins[7],ybins[i]);
	  pLV[i]=new TLine(xbins[i],ybins[0],xbins[i],ybins[7]);
	  pLH[i]->SetLineStyle(3);
	  pLV[i]->SetLineStyle(3);
	}
    }

  track0->Draw("-Xvb_tr:-Yvb_tr>>h2",cut2,"contztext");  
  gPad->SetLogz(1); 

  if (bDrawGrid)
    {
      for(int i=0;i<=7;i++)
	{
	  pLH[i]->Draw("same");
	  pLV[i]->Draw("same");
	}
    }

  c1->cd(0);
  c1->SaveAs(Form("SieveRates0%s.png",key));
  return;
}

void SieveRate1(char *key="_")
{
  TCanvas *c1=new TCanvas("c1","Sieve Rates",1000,700);

  TCut cut1 = "TrackClass>4 && abs(-Yvb_tr+5)<20 && abs(Xvb_tr)<50";
  TCut cut2 = "(TrackClass>4 && abs(-Yvb_tr+5)<20 && abs(Xvb_tr)<50)*ElasXS/1000"; //Weighted by cross section

  //variable bin size histogram did not work for TH2, only TH1
  //  -case 2  xbins!=0
  //   a new histogram is created (you should specify newname).
  //   The parameter is the number of variable size bins in the created histogram.
  //   The array xbins must contain ngroup+1 elements that represent the low-edge
  //   of the bins.
  //   If the original histogram has errors stored (via Sumw2), the resulting
  //   histograms has new errors correctly calculated.
  //
  //   examples: if h1 is an existing TH1F histogram with 100 bins
  //     Double_t xbins[25] = {...} array of low-edges (xbins[25] is the upper edge of last bin
  //     h1->Rebin(24,"hnew",xbins);  //creates a new variable bin size histogram hnew
  
  c1->Divide(2,2,0.01,0.001);
  //this is the perfect bin size, but not good in binning
  //double xbins[8]={-24.4,-18.3,-12.2,-6.1,0.0,4.8,9.6,14.4};
  //this is a work arround
  double xbins[8]={-24.0,-18.0,-12.0,-6.0,0.0,4.0,9.0,13.5};
  //vertical span=13.3mm x 7=93.1mm, from -43.6 to 49.5


  c1->cd(1);
  //track0->Draw("-Yvb_tr>>h0(288,-24.4,14.4)",cut2,"");    
  //track0->Draw("-Yvb_tr>>h0(75,-24,13.5)",cut2,"");
  TH1F *h0=new TH1F("h0","XS weighted yield for sieve holes; -Yvb_tr (mm)",75,-24,13.5);
  track0->Draw("-Yvb_tr>>h0",cut2,"");
  
  c1->cd(2);
  TH1F *h0s=h0->Rebin(7,"h0s",xbins);  
  h0s->SetMarkerSize(3);
  h0s->Draw("hist text0");

  c1->cd(3);
  TH1F *h1=new TH1F("h1","XS weighted yield for sieve holes; -Xvb_tr (mm)",7,-43.6,49.5);
  h1->SetMarkerSize(3);
  track0->Draw("-Xvb_tr>>h1",cut2,"histtext0");

  
  c1->cd(4);
  TH2F *h2=new TH2F("h2","XS weighted yield for sieve holes; -Yvb_tr (mm); -Xvb_tr (mm)",
		    25,-24,13.5,7,-43.6,49.5);
  h2->SetMinimum(5);
  h2->SetMarkerSize(2);
  track0->Draw("-Xvb_tr:-Yvb_tr>>h2",cut2,"colztext");
  gPad->SetLogz(1);   

  c1->cd(0);
  c1->SaveAs(Form("SieveRates1%s.png",key));
  return;
}


void SieveRate2(char *key="_")
{
  TCanvas *c1=new TCanvas("c1","Sieve Rates",1000,700);

  TCut cut1 = "TrackClass>4";
  TCut cut2 = "(TrackClass>4)*ElasXS/1000"; //Weighted by cross section

  c1->Clear();

  
  c1->Divide(1,2,0.001,0.001);
  TPad *pad1, *pad2;
  pad1=(TPad*)c1->cd(1);
  pad1->Divide(2,1,0.0,0.0);
  pad1->cd(1);
  track0->Draw("-Xvb_tr:-Yvb_tr>>h1(4,-24.4,0,7,-43.6,49.5)",cut1,"text");
  h1->SetTitle("Raw sieve slit yield; Horizontal (mm) ; Vertical (mm)");
  h1->SetMarkerSize(3);
  h1->Draw("text");

  pad1->cd(2);
  track0->Draw("-Xvb_tr:-Yvb_tr>>h2(3,0,14.4, 7,-43.6,49.5)",cut1,"text");
  h2->SetTitle("Raw sieve slit yield; Horizontal (mm) ; Vertical (mm)");
  h2->SetMarkerSize(3);
  h2->Draw("text");

  pad2=(TPad*)c1->cd(2);
  pad2->Divide(2,1,0.0,0.0);

  bool doLogz=false;
  pad2->cd(1);  //6.1 * 4  = 24.4 mm, start from 0; 
  track0->Draw("-Xvb_tr:-Yvb_tr>>h3(4,-24.4,0,7,-43.6,49.5)",cut2,"text");
  h3->SetTitle("Elas XS weighted sieve slit yield; Horizontal (mm) ; Vertical (mm)");
  if(h3->GetMaximum()>1000) doLogz=true;
  if(doLogz) gPad->SetLogz(1);
  h3->SetMarkerSize(3);
  h3->Draw("text");

  pad2->cd(2);  //4.8 * 3 = 14.4 mm, start from 0 
  track0->Draw("-Xvb_tr:-Yvb_tr>>h4(3,0,14.4,7,-43.6,49.5)",cut2,"text");
  h4->SetTitle("Elas XS weighted sieve slit yield; Horizontal (mm) ; Vertical (mm)");
  if(doLogz) gPad->SetLogz(1);
  h4->SetMarkerSize(3);
  h4->Draw("text");

  c1->cd(0);
  c1->SaveAs(Form("SieveRates2%s.png",key));

}



void SieveRate(char *key="L",int left=1)
{
  SieveRate_Core(key,left);
  SieveStat(key,left);
}
