//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Relative calibration of Si gains
////
////This step fixes the relative gains of the 3 remaining back channels to that of a designated front channel
////It works in the same was as step 2
////The back gain is equal to the inverse of the slope of the best fit line.
////
////First, histograms should be created for events in which one front strip (up and down) and one back strip fired. Using more complicated multiplicities here will confuse things (if only the up fired, then front != back, which defeats the initial assumption that front == back).
////Whichever gains file was used in the program that creates the histos should be the input file to this program.
////Input histogram file and relative gains file
////Things to watch out for:
////In the MyFit() funtion, the binning is assumed. This could cause runtime errors (probably won't crash though) if your binning is different than what the program assumes.
////
////root -l SiRelativeGains_Step3.C++
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

  Double_t x1[5] = { 2, 12, 9, 0.4, 2 };
  Double_t y1[5] = { 0.8, 9.8, 11, 1.1, 0.8 };
  TCutG *cut = new TCutG("cut",5,x1,y1);

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

void SiRelativeGains_Step3(void)
{
  using namespace std;

  TFile *f1 = new TFile("../../../OrganizeRaw_root/run567_051116.root");//front
  if ( !f1->IsOpen() ){
    cout << "Error: Root File Does Not Exist\n";
    exit(EXIT_FAILURE);
  }
  TCanvas *can = new TCanvas("can","can",800,600);
  
  ofstream outfile;
  outfile.open("X3RelativeGains051116_Step3.dat");

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

  for (Int_t DetNum=8; DetNum<9; DetNum++) {
    if (DetNum==7 || DetNum==9 || DetNum==19 || DetNum==24) {
      for (Int_t i=0; i<12; i++){
	outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
      }
      continue;
    }
    for (Int_t BackChNum=0; BackChNum<4; BackChNum++){
      Int_t FrontChNum = 1;
      //if (BackChNum==1){
      //continue;
      //}

      TH2F *hist = NULL;
      hist = (TH2F*)f1->Get(Form("back_vs_front%i_back%i",DetNum,BackChNum));
      if (hist==NULL){
	cout << "Histo does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	bad_back[count_bad] = BackChNum;
	count_bad++;
	continue;
      }

      Double_t gain = MyFit(hist,can);
      slope[DetNum-4][BackChNum] = slope[DetNum-4][BackChNum]/gain;
    }
    for (Int_t i=0; i<12; i++){
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

	    
	     
    
