//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Relative calibration of Si gains--Method2
////This method works very similar to method 1, except that instead of a predefined cut, the user must draw their own graphical cut
////To do this, in the canvas: View->Toolbar and click on the scissors on the top right hand corner
////Then, select the region around your data by clicking. Double click to close the cut.
////A best fit line should appear through your data
////
////Things to watch out for:
////In the MyFit() funtion, the binning is assumed. This could cause runtime errors (probably won't crash though) if your binning is different than what the program assumes.
////
//// Edited by : John Parker , 2016Jan22
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <fstream>
#include <exception>
#include <TCutG.h>
#include <TVector.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t MyFit(TH2F* hist, TCanvas *can){
  hist->Draw("colz");

  vector<double> x1;
  vector<double> y1;

  TCutG *cut;
  cut = (TCutG*)can->WaitPrimitive("CUTG");
  x1.resize(cut->GetN());
  y1.resize(cut->GetN());

  for(int n=0;n<cut->GetN();n++){
    cut->GetPoint(n,x1[n],y1[n]);
    cout << x1[n] << "\t" << y1[n] << endl;
  }

  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      if ( !cut->IsInside((Double_t)i*20,j*20) ){
	continue;
      }
      for (int k=0; k<hist->GetBinContent(i,j); k++){
	counter++;
      }
    }
  }

  Double_t *x = new Double_t[counter];
  Double_t *y = new Double_t[counter];

  counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      if ( !cut->IsInside((Double_t)i*20,j*20) ){
	continue;
      }
      for (int k=0; k<hist->GetBinContent(i,j); k++){
	x[counter] = i*20;
	y[counter] = j*20;
	counter++;
      }
    }
  }

  TGraph *graph = new TGraph(counter,x,y);
  graph->Draw("*same");

  TF1 *fun2 = new TF1("fun2","[0]*x",0,10000);
  graph->Fit("fun2");
  can->Update();
  can->WaitPrimitive();

  Double_t gain = fun2->GetParameter(0);
  delete x;
  delete y;
  delete graph;
  delete fun2;

  if (gain<0.1){
    gain = 1;
  }

  return gain;
}

void SiRelativeGains_Step1(void)
{
  using namespace std;

  TFile *f1 = new TFile("run235_245out_Step3_012816.root");//front
  if ( !f1->IsOpen() ){
    cout << "Error: Root File Does Not Exist\n";
    exit(EXIT_FAILURE);
  }
  TCanvas *can = new TCanvas("can","can",800,600);

  ofstream outfile;

  outfile.open("X3RelativeGains012816_Step4.dat");

  ifstream infile;
  infile.open("X3RelativeGains012816_Step3.dat");
  Int_t det=0,ch=0;
  Double_t slope[24][12];
  Double_t dummy;
  if (infile.is_open()){
    while (!infile.eof()){
      infile >> det >> ch >> dummy;
      slope[det-4][ch] = dummy;
    }
  }else{
    cout << "Infile not opened\n";
    exit(EXIT_FAILURE);
  }
  infile.close();

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t bad_back[288];
  Int_t count_bad = 0;

  for (Int_t DetNum=27; DetNum<28; DetNum++){
    for (Int_t i=0; i<4; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
    }
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      Int_t BackChNum = 0;
      if ( DetNum==11 ){
	BackChNum = 1;
      }

      TH2F *hist = NULL;
      hist = (TH2F*)f1->Get(Form("back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChNum));
     if (hist==NULL){
	cout << "Histo does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	bad_back[count_bad] = BackChNum;
	count_bad++;
	continue;
      }


      Double_t gain = MyFit(hist,can);
      slope[DetNum-4][FrontChNum+4] = slope[DetNum-4][FrontChNum+4]*gain;
      slope[DetNum-4][FrontChNum+8] = slope[DetNum-4][FrontChNum+8]*gain;

    }
    for (Int_t i=4; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
    }
  }
  outfile.close();
  cout << "List of bad detectos:\n";
  for (int i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << "  " << bad_back[i] << endl;
  }

  delete can;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		 

	    
	     
    
