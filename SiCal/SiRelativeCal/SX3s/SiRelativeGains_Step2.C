//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Relative calibration of Si gains
////This step fixes the relative gains of the 4 front strips with respect to a single back strip
////It works by dumping the x-y coordinates of a 2D histogram into a TGraph.
////---The coordinates are weighted by the bin content of the 2D hist
////---Only the coordinates that fall within a predefined cut are accepted. This gets rid of noise that can affect the fit. 
////However, on first iteration, not all of the data will fit neatly into the cut, so it may need to be adjusted for a few detectors. 
////Can use method 2 or method 3 to do this easily.
////It then fits the TGraph with a function of the form m*x+b.
////The two corresponding front gains (up and down) are then multiplied by the slope to get the new gain.
////-------------------------------------------------------------------------------------------------------
////Usage:
////First, histograms should be created for events in which one front strip (up and down) and one back strip fired. 
////Using more complicated multiplicities here will confuse things (if only the up fired, then front != back, 
////which defeats the initial assumption that front == back).
////
////Whichever gains file was used in the program that creates the histos should be the input file to this program.
////Input histogram file and relative gains file
////Things to watch out for:
////In the MyFit() funtion, the binning is assumed. This could cause runtime errors (probably won't crash though) 
////if your binning is different than what the program assumes.
////The above binning issue is probably OK now as it automatically extracts the correct number of bins 
////using the GetXmax() functions below  

////
////root -l SiRelativeGains_Step2.C++
////
//// Edited by : John Parker , 2016Jan22
//// Edited by : Maria Anastasiou, 2016Sept20
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
Double_t MyFit1(TH2F* hist, TCanvas *can){
  hist->Draw("colz");

  //if it is necessary to limit the area of your data that you want to fit see in our Canvas and 
  //input below on the CUT the coordinates of the points that surround this area. 

  //Double_t x1[5] = { 2, 12, 9, 0.4, 2 };
  //Double_t y1[5] = { 0.8, 9.8, 11, 1.1, 0.8 };
  //TCutG *cut = new TCutG("cut",5,x1,y1);
  //cut = (TCutG*)can->WaitPrimitive("CUTG"); //not necessary to uncomment
  //cut->Draw("same");                        //not necessary to uncomment

  Double_t maxbinNumberX = hist->GetXaxis()->GetXmax();
  Double_t maxbinNumberY = hist->GetYaxis()->GetXmax();
  Double_t maxbinX = (maxbinNumberX/hist->GetNbinsX());
  Double_t maxbinY = (maxbinNumberY/hist->GetNbinsY());

  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      //if ( !cut->IsInside((Double_t)i*maxbinX,j*maxbinY) ){
      //	continue;
      //}
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
      //if ( !cut->IsInside((Double_t)i*maxbinX,j*maxbinY) ){
      //	continue;
      // }
      for (int k=0; k<hist->GetBinContent(i,j); k++){
	x[counter] = i*maxbinX;
	y[counter] = j*maxbinY;
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

Double_t MyFit3(TH2F* hist, TCanvas *can){//cut-as-line fit; copied from Old/Clickable_Step3
  hist->Draw("colz");

  Double_t x[10];
  Double_t y[10];

  TCutG *cut;
  cut = (TCutG*)can->WaitPrimitive("CUTG");
      
  for(int n=0;n<cut->GetN()-2;n++){
    cut->GetPoint(n,x[n],y[n]);
    cout << x[n] << "\t" << y[n] << endl;
  }
	
  TGraph *graph = new TGraph(cut->GetN()-2,x,y);
  hist->Draw("colz");
  graph->Draw("*same");
	
  TF1 *fun = new TF1("fun","[0] + [1]*x",0,1);
  graph->Fit("fun");
	
  can->Update();
  can->WaitPrimitive();

  Double_t gain = fun->GetParameter(1);
      
  delete graph;
  delete fun;
      
  return gain;
}

void SiRelativeGains_Step2(void)
{
  using namespace std;

  //input the root file that was created in the Main(Organize) using the X3RelativeGains_Step1.dat file

  TFile *f1 = new TFile("/data0/manasta/OrganizeRaw_files/run930_931_nospacer_X3step1_divideback.root");
  if ( !f1->IsOpen() ){
    cout << "Error: Root File Does Not Exist\n";
    exit(EXIT_FAILURE);
  }
  
  ofstream outfile;
  outfile.open("X3RelativeGains09202016_Step2_all.dat"); //output file

  ifstream infile;
  infile.open("X3RelativeGains09192016_Step1_det27.dat"); //input the .dat file that was created in Step1 process
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

  TCanvas *can = new TCanvas("can","can",800,600);

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t bad_back[288];
  Int_t count_bad = 0;

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      Int_t BackChNum = 0;   // some if-statements that differ between each data set
      if ( DetNum==8 ){
      	BackChNum = 3;
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

      Double_t gain = MyFit1(hist,can);
      slope[DetNum-4][FrontChNum+4] = slope[DetNum-4][FrontChNum+4]*gain;
      slope[DetNum-4][FrontChNum+8] = slope[DetNum-4][FrontChNum+8]*gain;
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
}
