{
  TFile *f = new TFile("run235_245out_Slope1.root");

  TCanvas *can = new TCanvas("can","can",900,1100);

  can->Divide(1,2);
  TF1 *fun2 = new TF1("fun2","x",0,10000);
  can->cd(1);
  QQQback_vs_front->Draw("colz");

  can->cd(2);
  TFile *f2 = new TFile("run235_245out_Step2_012816.root");
  QQQback_vs_front->Draw("colz");
  //fun2->Draw("same");
}
