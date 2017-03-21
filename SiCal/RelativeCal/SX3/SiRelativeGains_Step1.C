//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Relative calibration of Si gains
////
////General Usage of all three steps
////
////If you loop over a subset of detectors, only the paramaters for those detectors will be written into the new .dat file
////First, run Organize.cpp with a given X3cal.dat file
////   Input organize.root and the .dat file into this code
////   This code will output another .dat file that you should define
////   root -l SiRelativeGains_Step1.C++
////Second, run Organize.cpp again with the new .dat file
////   Input the new organize.root and new .dat file into Step2.C
////   root -l SiRelativeGains_Step2.C++
////Finally, run Organize.cpp again with the new .dat file
////   Input the new organize.root and new .dat file into Step3.C
////   root -l SiRelativeGains_Step3.C++
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////This program fixes the relative gains for the up and down on the SX3
////How to use: Create Histograms in the Main.C(OrganizeIC.cpp):
////-----Plot down vs up energy for each front channel of each detector
////-----Code reads in histogram of name down_vs_up%i_front%i--e.g. down_vs_up4_front0 (det4, channel0)
////-----It can loop over any range of detectors you want, but if the histogram does not exist, the code will crash and your work will not be saved
////
////
////The program reads in an X3RelativeGains.dat file and outputs a new file with updated coefficients
////Before running, make sure that the root file you are reading in has the right histograms and has been created using the X3RelativeGains.dat file that you are inputting into this code
////
////This file changes the relative gains on the down relative to the up
////Once you have completed this program, rerun Main.C(OrganizeIC.cpp) with this new X3RelativeGains_Step1.dat
////You will input this new root file with this new relative gains file into step 2.
////
//// Edited by : John Parker , 2016Jan22
//// Edited by : Maria Anastasiou, 2016Sept20
////
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
Double_t MyFit(TH2F* hist, TCanvas* can){
  hist->Draw("colz");
  //hist->GetXaxis()->SetRange(0,180);
  //hist->GetYaxis()->SetRange(0,180);
  
  //if it is necessary to limit the area of your data that you want to fit see in our Canvas and 
  //input below on the CUT the coordinates of the points that surround this area. 
 
  //Double_t x1[8] = {190, 240, 240, 13900, 4230, 1700, 165, 190};
  //Double_t y1[8] = {315, 3215, 13100, 780, 420, 170, 136, 315};
  //TCutG *cut = new TCutG("cut",8,x1,y1);
  
  Double_t maxbinNumberX = hist->GetXaxis()->GetXmax();
  Double_t maxbinNumberY = hist->GetYaxis()->GetXmax();
  Double_t maxbinX = (maxbinNumberX/hist->GetNbinsX());
  Double_t maxbinY = (maxbinNumberY/hist->GetNbinsY());

  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      //if ( !cut->IsInside((Double_t)i*maxbinX,(Double_t)j*maxbinY) ){
      //continue;
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
      //if ( !cut->IsInside((Double_t)i*maxbinX,(Double_t)j*maxbinY) ){
      //continue;
      //}
      for (int k=0; k<hist->GetBinContent(i,j); k++){
	x[counter] = (Double_t)i*maxbinX;
	y[counter] = (Double_t)j*maxbinY;
	counter++;
      }
    }
  }

  TGraph *graph = new TGraph(counter,x,y);
  graph->Draw("*same");

  TF1 *fun2 = new TF1("fun2","[0]+[1]*x",0,16000);
  //fun2->SetParameter(0,10);
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

void SiRelativeGains_Step1(void)
{
  using namespace std;

  //if you are just starting the calibration you can use an .dat input file where ALL SLOPES ARE ONE(1) apart from the MASKED CHANNELS which are ZERO(0)

  //input the root file that was created in the Main(Organize) using the .dat file where ALL SLOPES are ONE.

  //TFile *f1 = new TFile("/data0/manasta/OrganizeRaw_files/run930_931_nospacer_X3slope1_divideback.root"); // root file name created with Main(Organize)
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9m.root");

  TCanvas *can = new TCanvas("can","can",800,600);
  
  ofstream outfile;

  Double_t average_slope = 0;
  
  outfile.open("X3RelativeGains_Step1.dat"); //output file name

  ifstream infile;
  infile.open("X3RelativeGains_Slope1.dat"); //input file name
  Int_t det=0,ch=0;
  Double_t dummy_slope = 0;
  Double_t slope[24][12];
  if (infile.is_open()){
    while (!infile.eof()){
      infile >> det >> ch >> dummy_slope;
      slope[det-4][ch] = dummy_slope;
    }
  }else{
    cout << "Infile not opened\n";
    exit(EXIT_FAILURE);
  }
  infile.close();

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t count_bad = 0;

  //You can do the calibration detector by detector looping one detector at a time or loop from DetNum=4 to 27.
  //if you do so change the loop below.

  //histo "down_vs_up_divideBack%i_front_%i" is a new extra histo that was created in the OrganizeIC.cpp for data that are normalized with the BackEnergy
  //if your data is not normalized use histo "down_vs_up%i_front_%i" for this code. 

  for (Int_t DetNum=4; DetNum<28; DetNum++){
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      TH2F *hist = NULL;
      TString hname=Form("down_vs_up%i_f%i",DetNum,FrontChNum);
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << hname << " histogram does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	count_bad++;
	continue;
      }
      
      average_slope = MyFit(hist,can);
      slope[DetNum-4][FrontChNum+8] = -slope[DetNum-4][FrontChNum+8]/average_slope;
    }
    for (Int_t i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
    }
  }
  outfile.close();
  cout << "List of bad detectors:\n";
  for (int i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << endl;
  }
  delete can;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		 

	    
	     
    
