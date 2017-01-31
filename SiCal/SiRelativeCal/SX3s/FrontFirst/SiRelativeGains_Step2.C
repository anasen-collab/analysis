//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Relative calibration of Si gains
////This step fixes the relative gains of the 4 front strips with respect to a single back strip
////It works by dumping the x-y coordinates of a 2D histogram into a TGraph.
////---The coordinates are weighted by the bin content of the 2D hist
////---Only the coordinates that fall within a predefined cut are accepted. This gets rid of noise that can affect the fit. However, on first iteration, not all of the data will fit neatly into the cut, so it may need to be adjusted for a few detectors. Can use method 2 or method 3 to do this easily.
////It then fits the TGraph with a function of the form m*x+b.
////The two corresponding front gains (up and down) are then multiplied by the slope to get the new gain.
////-------------------------------------------------------------------------------------------------------
////Usage:
////First, histograms should be created for events in which one front strip (up and down) and one back strip fired. Using more complicated multiplicities here will confuse things (if only the up fired, then front != back, which defeats the initial assumption that front == back).
////Whichever gains file was used in the program that creates the histos should be the input file to this program.
////Input histogram file and relative gains file
////Things to watch out for:
////In the MyFit() funtion, the binning is assumed. This could cause runtime errors (probably won't crash though) if your binning is different than what the program assumes.
////
////root -l SiRelativeGains_Step1.C++


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

  Double_t x1[5] = { 2, 12, 9, 0.4, 2 };
  Double_t y1[5] = { 0.8, 9.8, 11, 1.1, 0.8 };
  TCutG *cut = new TCutG("cut",5,x1,y1);
  //cut = (TCutG*)can->WaitPrimitive("CUTG");
  cut->Draw("same");

  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      if ( !cut->IsInside((Double_t)i*0.1,j*0.1) ){
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
      if ( !cut->IsInside((Double_t)i*0.1,j*0.1) ){
	continue;
      }
      for (int k=0; k<hist->GetBinContent(i,j); k++){
	x[counter] = i*0.1;
	y[counter] = j*0.1;
	counter++;
      }
    }
  }

  TGraph *graph = new TGraph(counter,x,y);
  graph->Draw("*same");

  TF1 *fun2 = new TF1("fun2","[0]*x +[1]",0,10000);
  graph->Fit("fun2");
  can->Update();
  can->WaitPrimitive();

  Double_t gain = fun2->GetParameter(0);
  delete x;
  delete y;
  delete graph;
  delete fun2;

  return gain;
}

void SiRelativeGains_Step2(void)
{
  using namespace std;

  TFile *f1 = new TFile("../../../OrganizeRaw_root/run567_051116.root");//front
  if ( !f1->IsOpen() ){
    cout << "Error: Root File Does Not Exist\n";
    exit(EXIT_FAILURE);
  }
  TCanvas *can = new TCanvas("can","can",800,600);

  ofstream outfile;
  outfile.open("X3RelativeGains051116.dat");

  ifstream infile;
  infile.open("X3RelativeGains022316.dat");
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

  for (Int_t DetNum=8; DetNum<9; DetNum++){
    for (Int_t i=0; i<4; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
    }
    if (DetNum==19 || DetNum==24){
      for (Int_t i=4; i<12; i++){
	outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
      }
      continue;
    }
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      Int_t BackChNum = 1;
      if ( DetNum==11 ){
	BackChNum = 1;
      }
      if (DetNum==24 ){
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
  cout << "List of bad detectors:\n";
  for (int i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << "  " << bad_back[i] << endl;
  }

  delete can;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		 

	    
	     
    
