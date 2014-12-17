crazy(){
  ifstream INFILE;
  INFILE.open("parse.dat");
  Double_t B1, B2, B3;
  Double_t BB[25215];
  Double_t ii[25215];
  Int_t    index = 0;
  for( Int_t i = 0; i < 25215; i++ ) {
    INFILE >> B1 >> B2 >> B3;
    //cout << B1 << " " << B2 << " " << B3 << endl;
    Double_t B = TMath::Sqrt( B1 * B1 + B2 * B2 + B3 * B3 );
    BB[index] = B;
    ii[index] = index;
    index++;
  }

  TGraph* gr = new TGraph(25215, ii, BB);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.5);
  gr->Draw();

}

