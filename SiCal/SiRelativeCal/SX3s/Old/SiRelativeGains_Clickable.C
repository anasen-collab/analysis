//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////This program fixes the relative gains for the up and down on the SX3
////How to use: Create Histograms in the Main.C:
////-----Plot down vs up energy for each front channel of each detector
////-----Code reads in histogram of name down_vs_up%i_front%i--e.g. down_vs_up4_front0 (det4, channel0)
////-----It can loop over any range of detectors you want, but if the histogram does not exist, the code will crash and your work will not be saved
////
////-----Craching fixed by putting an if-statement when histo is NULL continue //M.Anastasiou
////
////The program reads in an X3RelativeGains.dat file and outputs a new file with updated coefficients
////Before running, make sure that the root file you are reading in has the right histograms and has been created using the X3RelativeGains.dat file that you are inputting into this code
////
////This file changes the relative gains on the down relative to the up
////
////root -l SiRelativeGains_Clickable.C++
////in canvas: View->Toolbar->GraphicalCut (pair of scissors on right)
////the plot should look like a straight line. click along the straight line. when you are done, double click in canvas and a best fit line will appear. The best fit line should follow the data very well, if not you are doing something wrong.
////Because the back gains are not set properly, the line may appear segmented. If so, choose your favorite segment and get a best fit line for that. Do not click in each segment as that will throw your best fit line off. All of the segments for a detector should have the same slope.
////
////After the best fit line appears, double click to move onto the next channel
////
////
////Once you have completed this program, rerun Main.C with this new X3RelativeGains_Step1.dat
////You will input this new root file with this new relative gains file into step 2.

//// Edited by : John Parker , 2016Jan24
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void SiRelativeGains_Clickable(void){

  using namespace std;
  Double_t x[10];
  Double_t y[10];
  TFile *f1 = new TFile("run236out_nocal.root");

  ofstream outfile;
  outfile.open("X3RelativeGains012216_Step1_dummy.dat");

  ifstream infile;
  infile.open("X3RelativeGains012216_slope1.dat"); //input file name
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

  TCanvas *can = new TCanvas("can","can",800,600);
  TCutG *cut;

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t count_bad = 0;

  for (Int_t DetNum=16; DetNum<17; DetNum++){
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      if (DetNum==22 && (FrontChNum==2 || FrontChNum==3) ){
	continue;
      }
      if (DetNum==25 && FrontChNum==2){
	continue;
      }
      TH2F *hist = NULL;
      hist = (TH2F*)f1->Get(Form("down_vs_up%i_front%i",DetNum,FrontChNum));
      if (hist==NULL){
	cout << "Histo does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	count_bad++;
	continue;
      }
      
      hist->GetXaxis()->SetRange(0,120);
      hist->GetYaxis()->SetRange(0,120);
      hist->Draw("colz");
	
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
	
      slope[DetNum-4][FrontChNum+4] = -1./fun->GetParameter(1)*slope[DetNum-4][FrontChNum+4];
    }
  }
  for (Int_t i=0; i<24; i++){
    for (Int_t j=0; j<12; j++){
      outfile << i+4 << "\t" << j << "\t" << slope[i][j] << endl;
    }
  } 
  outfile.close();
  cout << "List of bad detectors:\n";
  for (int i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << endl;
  }

}




