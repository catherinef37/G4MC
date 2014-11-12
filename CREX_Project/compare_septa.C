compare_septa(){

  TCanvas* can = new TCanvas("can", "can", 20, 20, 1100, 700);
  can->cd();

  TH1F* dummy = new TH1F("compare septa fields", "", 1, -100., 100.);
  dummy->GetYaxis()->SetRangeUser(0., 10.);
  dummy->Draw();

  ifstream INFILE1;
  INFILE1.open("prex_septumfield.dat");
 
  skip_lines(INFILE1, 8);
  TString input = "";

  Float_t crazy1[25215];
  Float_t cray1 [25215];
  Int_t crayi1 = 0;
  Float_t crazy2[665691];
  Float_t cray2 [665691];
  Int_t crayi2 = 0;

  Float_t x = 30.;
  Float_t y = 0.;
  Float_t z = 2.;

  Float_t x_min1 =  -40.0;
  Float_t x_max1 =   40.0;
  Float_t y_min1 = -100.0;
  Float_t y_max1 =  100.0;
  Float_t z_min1 =  -14.0;
  Float_t z_max1 =   14.0;
  Float_t x_step1=    2.0;
  Float_t y_step1=    5.0;
  Float_t z_step1=    2.0;

  Float_t x_min2 =  -42.0;
  Float_t x_max2 =   42.0;
  Float_t y_min2 =   -9.5;
  Float_t y_max2 =    9.5;
  Float_t x_step2=    0.5;
  Float_t y_step2=    0.5;

  Int_t x_index1 = TMath::FloorNint( ( x - x_min1 ) / x_step1 );
  Int_t y_index1 = TMath::FloorNint( ( y - y_min1 ) / y_step1 );
  Int_t z_index1 = TMath::FloorNint( ( z - z_min1 ) / z_step1 );
  Int_t x_index2 = TMath::FloorNint( ( x - x_min2 ) / x_step2 );
  Int_t y_index2 = TMath::FloorNint( ( y - y_min2 ) / y_step2 );

  cout << x_index1 << " " << y_index1 << " " << x_index2 << " " << y_index2 << endl;

  Float_t X1 [41][15][41];
  Float_t Y1 [41][15][41];
  Float_t Z1 [41][15][41];
  Float_t BX1[41][15][41];
  Float_t BY1[41][15][41];
  Float_t BZ1[41][15][41];
  Float_t B1 [41][15][41];
  for( Int_t k = 0; k < 15; k++ ) {
    for( Int_t j = 0; j < 41; j++ ) {
      for( Int_t i = 0; i < 41; i++ ) {
	INFILE1 >> input; 
	X1 [i][k][j] = input.Atof();
	INFILE1 >> input;
	Y1 [i][k][j] = input.Atof();
	INFILE1 >> input;
	Z1 [i][k][j] = input.Atof();
	INFILE1 >> input;
	BX1[i][k][j] = input.Atof();
	INFILE1 >> input;
	BY1[i][k][j] = input.Atof();
	INFILE1 >> input;
	BZ1[i][k][j] = input.Atof();
	B1 [i][k][j] = TMath::Sqrt( BX1[i][k][j] * BX1[i][k][j] + BY1[i][k][j] * BY1[i][k][j] + BZ1[i][k][j] * BZ1[i][k][j] );
	//if( Z1[i][k][j] == 10. && X1[i][k][j] == 10. )
	  //cout << X1 [i][k][j] << " " << Y1 [i][k][j] << " " << Z1 [i][k][j] << " " << B1 [i][k][j] << endl;
	crazy1[crayi1] = B1[i][k][j];
	cray1 [crayi1] = crayi1;
	crayi1++;
      }
    }
  }

  Float_t the_B1[15];
  Float_t the_Z1[15];
  for(Int_t i = 0; i < 15; i++){
    the_B1[i] = B1 [i][y_index1][x_index1];
    the_Z1[i] = Y1 [i][y_index1][x_index1];
  }

  

  INFILE1.close();
  //TGraph* gr1 = new TGraph(41, Z1[x_index1][y_index1], B1[x_index1][y_index1]);
  TGraph* gr1 = new TGraph(41, Y1[x_index1][z_index1], B1[x_index1][z_index1]);
  cout << X1[x_index1][z_index1][0] << " " << Z1[x_index1][z_index1][0] << endl;
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.5);
  gr1->Draw("p");

  ifstream INFILE2;
  INFILE2.open("g2p_septumfield.dat");
  skip_lines(INFILE2, 8);

  Float_t X2 [169][39][101];
  Float_t Y2 [169][39][101];
  Float_t Z2 [169][39][101];
  Float_t BX2[169][39][101];
  Float_t BY2[169][39][101];
  Float_t BZ2[169][39][101];
  Float_t B2 [169][39][101];
  for( Int_t i = 0; i < 169; i++ ) {
    for( Int_t j = 0; j < 39; j++ ) {
      for( Int_t k = 0; k < 101; k++ ) {
	INFILE2 >> input; 
	X2 [i][j][k] = input.Atof();
	INFILE2 >> input;
	Y2 [i][j][k] = input.Atof();
	INFILE2 >> input;
	Z2 [i][j][k] = input.Atof();
	INFILE2 >> input;
	BX2[i][j][k] = input.Atof();
	INFILE2 >> input;
	BY2[i][j][k] = input.Atof();
	INFILE2 >> input;
	BZ2[i][j][k] = input.Atof();
	B2 [i][j][k] = TMath::Sqrt( BX2[i][j][k] * BX2[i][j][k] + BY2[i][j][k] * BY2[i][j][k] + BZ2[i][j][k] * BZ2[i][j][k] ) / 10000;
	crazy2[crayi2] = B2[i][j][k];
	cray2 [crayi2] = crayi2;
	crayi2++;

	//cout << X2 [i][j][k] << " " << Y2 [i][j][k] << " " << Z2 [i][j][k] << " " << B2 [i][j][k] << endl;
      }
    }
  }

  INFILE2.close();
  TGraph* gr2 = new TGraph(101, Z2[x_index2][y_index2], B2[x_index2][y_index2]);
  //for(Int_t i = 0; i < 39; i++){
  //cout << i << endl;
  cout << X2[x_index2][y_index2][0] << " " << Y2[x_index2][y_index2][0] << endl;
  //cout << X2[x_index2][i][0] << " " << Y2[x_index2][i][0] << endl;
  //}
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.5);
  gr2->SetMarkerColor(2);
  gr2->Draw("p");

  TCanvas* craycan = new TCanvas("craycan", "craycan", 40, 40, 1100, 700);
  TGraph* grgr = new TGraph(665691, cray2, crazy2);
  grgr->SetMarkerStyle(20);
  grgr->SetMarkerSize(0.5);
  //grgr->Draw();
  TGraph* grgrgr = new TGraph(25215, cray1, crazy1);
  grgrgr->SetMarkerStyle(20);
  grgrgr->SetMarkerSize(0.5);
  grgrgr->SetMarkerColor(1);
  grgrgr->Draw();

}

void skip_lines(std::istream& pStream, size_t pLines)
{
  std::string s;
  for (; pLines; --pLines)
    std::getline(pStream, s);
}
