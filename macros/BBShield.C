void BBShield(char *grname="PFraction_1cm_Lead.png")
{
gStyle->SetPadGridX(true);
gStyle->SetPadGridY(true);

TCanvas *c1=new TCanvas("c1","Shielding",1000,700);
c1->Divide(4,2);

c1->cd(1);
track0->Draw("Pvb/P0:P0 >> h1","TrackClass>=0","");
h1->SetTitle("Pvb/P0 Vs P0 for gamma;P_{0} (GeV); P_{vb}/P_{0}");
h1->Draw("colz");

c1->cd(2);
track0->Draw("P0 >> hp0all_0(100,0,1.0)","P0>0","");
hp0all_0->SetTitle("Thrown(Black), Accepted Thrown(Red) and Detected(Blue);P_{0} or P_{vb} (GeV)");
track0->Draw("P0 >> hp0_0(100,0,1.0)","Pvb>0","");
hp0_0->SetLineColor(2);
track0->Draw("Pvb >> hpvb_0(100,0,1.0)","Pvb>0","");
hpvb_0->SetLineColor(4);
hp0all_0->Draw();
hp0_0->Draw("same");
hpvb_0->Draw("same");

c1->cd(3);
track1->Draw("Pvb/P0:P0 >> h2","TrackClass>=0","");
h2->SetTitle("Pvb/P0 Vs P0 for Proton;P_{0} (GeV); P_{vb}/P_{0}");
h2->Draw("colz");

c1->cd(4);
track1->Draw("P0 >> hp0all_1(100,0,1.0)","P0>0","");
hp0all_1->SetTitle("Thrown(Black), Accepted Thrown(Red) and Detected(Blue);P_{0} or P_{vb} (GeV)");
track1->Draw("P0 >> hp0_1(100,0,1.0)","Pvb>0","");
hp0_1->SetLineColor(2);
track1->Draw("Pvb >> hpvb_1(100,0,1.0)","Pvb>0","");
hpvb_1->SetLineColor(4);
hp0all_1->Draw();
hp0_1->Draw("same");
hpvb_1->Draw("same");

c1->cd(5);
track2->Draw("Pvb/P0:P0 >> h3","TrackClass>=0","");
h3->SetTitle("Pvb/P0 Vs P0 for #pi^{-};P_{0} (GeV); P_{vb}/P_{0}");
h3->Draw("colz");

c1->cd(6);
track2->Draw("P0 >> hp0all_2(100,0,1.0)","P0>0","");
hp0all_2->SetTitle("Thrown(Black), Accepted Thrown(Red) and Detected(Blue);P_{0} or P_{vb} (GeV)");
track2->Draw("P0 >> hp0_2(100,0,1.0)","Pvb>0","");
hp0_2->SetLineColor(2);
track2->Draw("Pvb >> hpvb_2(100,0,1.0)","Pvb>0","");
hpvb_2->SetLineColor(4);
hp0all_2->Draw();
hp0_2->Draw("same");
hpvb_2->Draw("same");

c1->cd(7);
track3->Draw("Pvb/P0:P0 >> h4","TrackClass>=0","");
h4->SetTitle("Pvb/P0 Vs P0 for neutron;P_{0} (GeV); P_{vb}/P_{0}");
h4->Draw("colz");

c1->cd(8);
track3->Draw("P0 >> hp0all_3(100,0,1.0)","P0>0","");
hp0all_3->SetTitle("Thrown(Black), Accepted Thrown(Red) and Detected(Blue);P_{0} or P_{vb} (GeV)");
track3->Draw("P0 >> hp0_3(100,0,1.0)","Pvb>0","");
hp0_3->SetLineColor(2);
track3->Draw("Pvb >> hpvb_3(100,0,1.0)","Pvb>0","");
hpvb_3->SetLineColor(4);
hp0all_3->Draw();
hp0_3->Draw("same");
hpvb_3->Draw("same");

c1->SaveAs(grname);
}
