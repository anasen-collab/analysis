//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Relative calibration of Si gains
////
//// Usage: root -l SiPulser_All.C++ (from the same directory).
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t MyFit(TH2F* hist, TCanvas *can){
  hist->Draw("colz");

  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
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
      for (int k=0; k<hist->GetBinContent(i,j); k++){
	x[counter] = i*100;
	y[counter] = j*100;
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

void SiRelativeGains_Step2(void)
{
  using namespace std;

  TFile *f1 = new TFile("run236out_Step1.root");//front
  
  ofstream outfile;
  outfile.open("X3RelativeGains012216_NoClickStep2.dat");

  ifstream infile;
  infile.open("X3RelativeGains012216_Step1.dat");
  Int_t det=0,ch=0;
  Double_t slope[24][12];
  if (infile.is_open()){
    while (!infile.eof()){
      infile >> det >> ch >> slope[det-4][ch];
    }
  }else{
    cout << "Infile not opened\n";
    exit(EXIT_FAILURE);
  }
  infile.close();

  TCanvas *can = new TCanvas("can","can",800,1100);

  for (Int_t DetNum=4; DetNum<28; DetNum++){
    for (Int_t BackChNum=0; BackChNum<4; BackChNum++){
      if ( DetNum==11 && (BackChNum==0 || BackChNum==2 || BackChNum==3) ){
	continue;
      }
      if (DetNum==18 && BackChNum==3){
	continue;
      }
      if (DetNum==24 && BackChNum==3){
	continue;
      }
      if (DetNum==26 && (BackChNum==1 || BackChNum==2 || BackChNum==3) ){
	continue;
      }
      if (DetNum==27 && (BackChNum==2 || BackChNum==3) ){
	continue;
      }
      //cout << DetNum << "  " << BackChNum << endl;
      TH2F *hist;
      hist = (TH2F*)f1->Get(Form("back_vs_front%i_0_%i",DetNum,BackChNum));

      Double_t gain = MyFit(hist,can);
      slope[DetNum-4][BackChNum] = slope[DetNum-4][BackChNum]/gain;
    }
    for (int i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
    }
  }
  outfile.close();
  delete can;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		 

	    
	     
    
