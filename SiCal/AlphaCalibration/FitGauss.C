{
  TFile *f = new TFile ("run_alpha0_282_284.root");

  Double_t range_L = 6000;
  Double_t range_U = 7000;

  TF1 *combo = new TF1("combo","[0]*exp(-0.5*((x-[1])/[2])^2)",range_L,range_U);
  combo->SetParameter(0,500);
  combo->SetParameter(1,6200);
  combo->SetParameter(2,100);

  combo->SetParName(0,"amplitude");
  combo->SetParName(1,"mean");
  combo->SetParName(2,"sigma");
  TH1F *hist = (TH1F*)f->Get("Energy27");

  hist->Fit("combo","","",range_L,range_U);
  //  background->Draw();
  if (combo->GetParameter(2) < 0)
    combo->SetParameter(2,-combo->GetParameter(2));

  Double_t fwhm, area, mean;
  mean = combo->GetParameter(1);
  fwhm = 2*sqrt(2*log(2))*combo->GetParameter(2);
  area = sqrt(2*3.14159265358979)*combo->GetParameter(0)*combo->GetParameter(2)/1;

  //combo->Draw("same");

  //Double_t resolution = 104.3;
  //Double_t fwhm_ns = fwhm/resolution;
  //Double_t fwhm2_ns = fwhm2/resolution;

  cout << "mean 1: " << mean << endl;
  cout << "fwhm 1: " << fwhm << endl;
  cout << "Area 1: " << area << endl;
  //cout << "fwhm 1 (ns): " << fwhm_ns << "     " << "fwhm 2 (ns): " << fwhm2_ns << endl;
  //cout << "Resolution: " << resolution << endl;
}
