{
  TFile *f = new TFile("run_alpha0_282_284_geometry.root");
  TCanvas *can = new TCanvas("can","can",900,1100);

  can->Divide(3,4);
  Int_t counter = 0;
  for (Int_t DetNum=4; DetNum<16; DetNum++){
    counter++;
    TH1F *hist = NULL;
    TH1F *hist = (TH1F*)f->Get(Form("SX3ZposCal%i",DetNum));
    if (hist==NULL){
      continue;
    }
    can->cd(counter);
    //hist->Rebin(4);
    hist->Draw();
    hist->GetXaxis()->SetRange(350,800);
  }
  
  /*
  can->Divide(4,4);
  Int_t counter = 0;
  for (Int_t DetNum=0; DetNum<4; DetNum++){
    for (Int_t FrontN=0; FrontN<4; FrontN++){
      counter++;
      TH1F *hist = NULL;
      TH1F *hist = (TH1F*)f->Get(Form("ZPos5_%i_%i",DetNum,FrontN));
      if (hist==NULL){
	continue;
      }
      can->cd(counter);
      //hist->Rebin(4);
      hist->Draw();
    }
  }
  */


  /*
  TF1 *fun = new TF1("fun","1 - x",0,1);
  can->Divide(1,2);
  can->cd(1);
  down_vs_up->Draw("colz");
  down_vs_up->GetXaxis()->SetRange(0,120);
  down_vs_up->GetYaxis()->SetRange(0,120);
  fun->Draw("same");

  can->cd(2);
  back_vs_front->Draw("colz");
  */

}
