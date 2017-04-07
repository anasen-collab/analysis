{
  //TFile *f1 = new TFile("run235_245out_Slope1.root");
  TFile *f1 = new TFile("../../../../root/run1226-9m.root");

  TCanvas *can = new TCanvas("can","can",900,1100);

  can->Divide(1,2);
  TF1 *fun2 = new TF1("fun2","x",0,10000);
  can->cd(1);
  Q3_back_vs_front->Draw("colz");

  can->cd(2);
  //TFile *f2 = new TFile("run235_245out_Step2_012816.root");
  TFile *f2 = new TFile("../../../../root/run1226-9mQ1.root");
  Q3_back_vs_front->Draw("colz");
  //fun2->Draw("same");
}
