{
  TFile *f = new TFile("run235_245out_Step1_012816.root");
  TCanvas *can = new TCanvas("can","can",900,1100);

  /*
  TH2F *hist = up_vs_down4_back0;
  TH2F *hist2 = up_vs_down4_back1;
  TH2F *hist3 = up_vs_down4_back2;
  TH2F *hist4 = up_vs_down4_back3;
  */

  Double_t x1[5] = { 1550, 8071, 7065, 1000, 1550 };
  Double_t y1[5] = { 1100, 7150, 8500, 2400, 1100 };
  TCutG *cut = new TCutG("cut",5,x1,y1);

  can->Divide(2,2);
  Int_t counter = 1;
  for (Int_t DetNum=0; DetNum<4; DetNum++){
    TH2F *hist = (TH2F*)f->Get(Form("QQQback_vs_front%i",DetNum));
    TF1 *fun = new TF1("fun","x",0,10000);
    can->cd(counter);
    hist->Draw("colz");
    //cut->Draw("same");
    counter++;
    fun->Draw("same");
  }
  
  /*
  TF1 *fun = new TF1("fun","1 - x",0,1);

  QQQback_vs_front->Draw("colz");
  */
}
