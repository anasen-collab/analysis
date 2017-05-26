//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////This program fixes the relative gains for the up and down on the PC
//// Usage: root -l PosCal.C
//// Edited by : John Parker , 2016Jan24
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//C++
#include <fstream>
#include <exception>
//ROOT
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TCutG.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t MyFit(TH2F* hist, TCanvas* can) {
  hist->Draw("colz");

  vector<double> x1;
  vector<double> y1;

  TCutG *cut;
  cut = (TCutG*)can->WaitPrimitive("CUTG");
  x1.resize(cut->GetN());
  y1.resize(cut->GetN());
  
  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++) {
    for (int j=1; j<hist->GetNbinsY(); j++) {
      if ( !cut->IsInside((Double_t)i*0.1667,(Double_t)j*0.1667) ) {
	continue;
      }
      for (int k=0; k<hist->GetBinContent(i,j); k++) {
	counter++;
      }
    }
  }
  
  Double_t *x = new Double_t[counter];
  Double_t *y = new Double_t[counter];

  counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++) {
    for (int j=1; j<hist->GetNbinsY(); j++) {
      if ( !cut->IsInside((Double_t)i*0.1667,(Double_t)j*0.1667) ) {
	continue;
      }
      for (int k=0; k<hist->GetBinContent(i,j); k++) {
	x[counter] = (Double_t)i*0.1667;
	y[counter] = (Double_t)j*0.1667;
	counter++;
      }
    }
  }
  
  TGraph *graph = new TGraph(counter,x,y);
  graph->Draw("*same");

  TF1 *fun2 = new TF1("fun2","[0]+[1]*x",0,6000);
  fun2->SetParameter(0,5000);
  fun2->SetParameter(1,-1);
  graph->Fit("fun2");
  can->Update();

  //cout << fun2->GetChisquare() << "  " << fun2->GetChisquare()/counter;

  can->WaitPrimitive();

  Double_t gain = fun2->GetParameter(1);

  delete x;
  delete y;
  delete graph;
  delete fun2;

  if (gain < -3 || gain > -0.1 ){//basically if the slope is very far from one, I don't trust the fit
    gain = -1;
  }
  return gain;
}

void PosCal(void) {
  using namespace std;
  Double_t x[10];
  Double_t y[10];
  string rootfile = "/home/lighthall/anasen/root/track/spacer0t.root";
  TFile *f1 = new TFile(rootfile.c_str());
  if (!f1->IsOpen()){
    cout << "Root file: " << rootfile << " could not be opened.\n";
    exit(EXIT_FAILURE);
  }

  TCanvas *can = new TCanvas("can","can",800,600);

  Double_t slope[24];
  Double_t offset[24];
  Double_t dummy_slope = 0;
  Double_t dummy_offset = 0;
  Int_t wire = 0;
  ifstream infile;
  string dummy;
  infile.open("saves/PCWireCal_init.dat");
  if (infile.is_open()){
    infile >> dummy >> dummy >> dummy;
    while (!infile.eof()){
      infile >> wire >> dummy_slope >> dummy_offset;
      slope[wire] = dummy_slope;
      offset[wire] = dummy_offset;
    }
  }else{
    cout << "Infile not opened\n";
    exit(EXIT_FAILURE);
  }
  infile.close();

  ofstream outfile;
  outfile.open("saves\PCWireCal_170526.dat");
  outfile << "Wire\tSlope\tShift\n";

  TCutG *cut;
  for (Int_t DetNum=0; DetNum<24; DetNum++){
    TH2F *hist = NULL;
    hist = (TH2F*)f1->Get(Form("PCZ_vs_Z%i",DetNum));
    if (hist==NULL){
      outfile << DetNum << "\t" << 1 << "\t" << 0 << endl;
      cout << "Warning: hist with wire number " << DetNum << " does not exist.\n";
	continue;
    }
    hist->Draw("colz");
    cut = (TCutG*)can->WaitPrimitive("CUTG");
	
    for(int n=0;n<cut->GetN()-2;n++){
      cut->GetPoint(n,x[n],y[n]);
      cout << x[n] << "\t" << y[n] << endl;
    }
    TGraph *graph = new TGraph(cut->GetN()-2,x,y);
    hist->Draw("colz");
    graph->Draw("*same");

    TF1 *fun = new TF1("fun","[0] + [1]*x",0,30);
    graph->Fit("fun");
	
    can->Update();
    can->WaitPrimitive();

    slope[DetNum] = slope[DetNum]/fun->GetParameter(1);
    offset[DetNum] = offset[DetNum] - fun->GetParameter(0);

    outfile << DetNum << "\t" << slope[DetNum]*fun->GetParameter(1) << "\t" << offset[DetNum] << endl;

  }  

}




