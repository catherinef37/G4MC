#include <TH1.h>
#include <TCanvas.h>
#include <TQObject.h>
#include <TString.h>

void PickSieve(char *key="")
{
   // Example of using signal/slot in TCanvas/TPad to get feedback
   // about processed events. 
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
  
  if(gROOT->GetListOfFiles()->GetEntries()<1) {
    cout<<"no root file opened yet, quit ... \n"; 
    return;
  }

  TString filename=_file0->GetName();
  if(strlen(key)<2)
    {
      int start=0;
      if(filename.Last('/')>=0) filename.Remove(0,filename.Last('/')+1);
      filename.Remove(filename.Length()-5,5);
      key=filename.Data();
    }

  TCanvas *pCan=new TCanvas("pCan","Sieve Slit",800,600);
  pCan->Divide(2,1);
  pCan->cd(1);
  TH2F *h2VH=0;
  h2VH=(TH2F*)gROOT->FindObject("h2VH");
  if (!h2VH)
    {
      h2VH=new TH2F("h2VH",
		    Form("%s; -Yvb_tr (mm); -Xvb_tr (mm)",key),
		    50,-24.0,13.5,42,-43.6,49.5);
    }
  track0->Draw("-Xvb_tr:-Yvb_tr >> h2VH","TrackClass>5","contz");

  pCan->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
               "DoEvent(Int_t,Int_t,Int_t,TObject*)");
}

void DoEvent(Int_t event, Int_t x, Int_t y, TObject *selected)
{
  //do not response to mouse movement
  if(event>20) return;
  if(strcmp("h2VH",selected->GetName())!=0) return;

  TCanvas *c = (TCanvas *) gTQSender;
  double px = gPad->AbsPixeltoX(x);
  double p2x = gPad->PadtoX(px);
  double py = gPad->AbsPixeltoY(y);
  double p2y = gPad->PadtoY(py);
  //printf("Canvas %s: event=%d, x=%d, y=%d, selected_type=%s, selected=%s, p2x=%f, p2y=%f\n", c->GetName(),
  //	 event, x, y, selected->IsA()->GetName(),selected->GetName(),px,py);
  
  if(event==11) DrawElipse(x,y);
  else if(event==12) DrawElipse(x,y,3.0,6.0);
  
}

void DrawElipse(double px,double py,double a=2.0,double b=3.5)
{
  TVirtualPad *padsav = gPad;
  double x = gPad->AbsPixeltoX(px);
  double y = gPad->AbsPixeltoY(py);

  if(gROOT->FindObject("CUTG")) delete gROOT->FindObject("CUTG");
  const int kNpt=9;
  TCutG *cutg = new TCutG("CUTG",kNpt);
  cutg->SetVarX("-Yvb_tr");
  cutg->SetVarY("-Xvb_tr");
  cutg->SetTitle("SieveHoleCut");
  cutg->SetFillColor(0);
  cutg->SetMarkerStyle(20);
  cutg->SetLineWidth(2);

  for(int i=0;i<=kNpt;i++)
    {
      double phi=i*4*acos(0)/kNpt;       
      cutg->SetPoint(i,x+a*cos(phi),y+b*sin(phi));
    }
  
  cutg->Draw("");
  gPad->Update();

  TCanvas *c = (TCanvas *) gTQSender;
  UpdateCanvas(c,x,y);
  padsav->cd();
}

void UpdateCanvas(TCanvas *c, double x, double y)
{
  c->cd(2);
  double x0=x,y0=y;
  int holeindex=GetHoleIndex(x,y,x0,y0);
  TH2F *h2Hole=0;
  h2Hole=(TH2F*)gROOT->FindObject("h2Hole");
  if(h2Hole) {/*cout<<"recreate "<<h2Hole->GetName()<<endl;*/delete h2Hole;}
  h2Hole=new TH2F("h2Hole",Form("Hole %02d; -Yvb_tr(mm);-Xvb_tr(mm)",holeindex),
		  8,x0-3.0,x0+3.0,6,y0-6.65,y0+6.65);		   
  h2Hole->SetMarkerSize(2);
  track0->Draw("-Xvb_tr:-Yvb_tr>>h2Hole","CUTG","contztext0");

}

//get the index and return the center position
int GetHoleIndex(double x, double y, double &x0, double &y0)
{
  double xbinsL[8]={-24.0,-18.0,-12.0,-6.0,0.0,4.5,9.0,13.5};
  double xbinsR[8]={-13.5,-9.0,-4.5,0.0,6.0,12.0,18.0,24.0};
  double *xbins=xbinsL;
  
  int left=1;
  //determine if it is left or right arm  
  track0->Draw("Xvb>>hXvb","TrackClass>5");
  
  double MeanXvb=hXvb->GetMean();  
  if(MeanXvb<-40) left=0;
  if(left==0) xbins=xbinsR;
  
  //vertical span=13.3mm x 7=93.1mm, from -43.6 to 49.6
  double ybins[8]={-43.6,-30.3,-17.0,-3.7,9.6,22.9,36.2,49.5};
  int h=-1,v=-1;
  for(int i=0;i<8;i++)
    {
      if(x>=xbins[i] && x<xbins[i+1]) {h=i;x0=(xbins[h]+xbins[h+1])/2;break;} 
    }
  for(int j=0;j<8;j++)
    {
      if(y>=ybins[j] && y<ybins[j+1]) {v=j;y0=(ybins[v]+ybins[v+1])/2;break;} 
    }

  int idx=(h<0 || v<0)?-1:10*h+v; 
  cout<<"x="<<x<<" y="<<y<<" ==> holeindex="<<idx<<endl;
  return idx;
}
