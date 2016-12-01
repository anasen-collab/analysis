//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////This program fixes the back gains in the SX3 relative to a front channel
////
////The root file should have histograms that are plotting back_vs_front
////-----Code reads in histogram of name back_vs_front%i_0_%i--e.g. down_vs_up4_0_1 (det4, front channel0, back channel 1)
////-----If the front channel 0 does not exist for a given detector, choose a different front channel for everything to be relative too.
////-----It can loop over any range of detectors you want, but if the histogram does not exist, the code will crash and your work will not be saved
////
////The program reads in an X3RelativeGains.dat file and outputs a new file with updated coefficients
////Before running, make sure that the root file you are reading in has the right histograms and has been created using the X3RelativeGains.dat file that you are inputting into this code
////
////root -l SiRelativeGains_Clickable_Step2.C++
////in canvas: View->Toolbar->GraphicalCut (pair of scissors on right)
////the plot should look like a straight line. click along the straight line. when you are done, double click in canvas and a best fit line will appear. The best fit line should follow the data very well, if not you are doing something wrong.
////
////After the best fit line appears, double click to move onto the next channel
////
////
////Once you have completed this program, rerun Main.C with this new X3RelativeGains_Step2.dat
////You will input this new root file with this new relative gains file into step 3.
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


void SiRelativeGains_Clickable_Step2(void){

  using namespace std;
  Double_t x[10];
  Double_t y[10];
  TFile *f1 = new TFile("run236out_Step1.root");

  ofstream outfile;
  outfile.open("X3RelativeGains012216_Step2_dummy.dat");


  ifstream infile;
  infile.open("X3RelativeGains012216_Step1_dummy.dat");
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
  }

  TCanvas *can = new TCanvas("can","can",800,600);
  TCutG *cut;

  Double_t new_slope[24][12];
  for (int i=0; i<24; i++){
    for (int j=0; j<12; j++){
      new_slope[i][j] = 1;
    }
  }

  for (Int_t DetNum=16; DetNum<17; DetNum++){
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

      TH2F *hist = (TH2F*)f1->Get(Form("back_vs_front%i_0_%i",DetNum,BackChNum));
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
	
      slope[DetNum-4][BackChNum] = 1./fun->GetParameter(1)*slope[DetNum-4][BackChNum];
	
    }
  }
  for (Int_t i=0; i<24; i++){
    for (Int_t j=0; j<12; j++){
      outfile << i+4 << "\t" << j << "\t" << slope[i][j] << endl;
    }
  }

}




