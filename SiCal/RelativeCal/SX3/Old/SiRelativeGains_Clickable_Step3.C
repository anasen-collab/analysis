//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////This program fixes the front gains vs the back channels.
////----Code reads in histogram of name back_vs_front%i_%i_0--e.g. down_vs_up4_1_0 (det4, front channel1, back channel 0)
////----If step 2 fixed the gains relative to front channel 0, in this code, you should loop over front channels 1,2,3
////-----If the back channel 0 does not exist for a given detector, choose a different front channel for everything to be relative too (in this code, det 11 uses back channel 1).
////-----It can loop over any range of detectors you want, but if the histogram does not exist, the code will crash and your work will not be saved
////
////The program reads in an X3RelativeGains.dat file and outputs a new file with updated coefficients
////Before running, make sure that the root file you are reading in has the right histograms and has been created using the X3RelativeGains.dat file that you are inputting into this code
////
////root -l SiRelativeGains_Clickable_Step3.C++
////in canvas: View->Toolbar->GraphicalCut (pair of scissors on right)
////the plot should look like a straight line. click along the straight line. when you are done, double click in canvas and a best fit line will appear. The best fit line should follow the data very well, if not you are doing something wrong.
////
////After the best fit line appears, double click to move onto the next channel
////
////
////Once you have completed this program, rerun Main.C with this new X3RelativeGains_Step3.dat
////Everything should now be calibrated properly, but check the histograms to make sure.
////If the histograms are not as good as desired, you can repeat this process for any individual detector (all three programs). Just be sure to input the correct RelativeGains.dat file and root file.
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


void SiRelativeGains_Clickable_Step3(void){

  using namespace std;
  Double_t x[10];
  Double_t y[10];
  TFile *f1 = new TFile("run236out_Step2.root");

  ofstream outfile;
  outfile.open("X3RelativeGains012216_Step3_dummy.dat");


  ifstream infile;
  infile.open("X3RelativeGains012216_Step2_dummy.dat");
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

  for (Int_t DetNum=16; DetNum<17; DetNum++){
    for (Int_t FrontChNum=1; FrontChNum<4; FrontChNum++){
      Int_t BackChannelNumber = 0;
      if ( DetNum==11 ){
	BackChannelNumber = 1;
      }
      if (DetNum==22 && (FrontChNum==2 || FrontChNum==3) ){
	continue;
      }
      if (DetNum==25 && FrontChNum==2){
	continue;
      }
      
      TH2F *hist = (TH2F*)f1->Get(Form("back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChannelNumber));
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
	
      slope[DetNum-4][FrontChNum+4] = slope[DetNum-4][FrontChNum+4]*fun->GetParameter(1);
      slope[DetNum-4][FrontChNum+8] = slope[DetNum-4][FrontChNum+8]*fun->GetParameter(1);
	
    }
  }
  for (Int_t i=0; i<24; i++){
    for (Int_t j=0; j<12; j++){
      outfile << i+4 << "\t" << j << "\t" << slope[i][j] << endl;
    }
  }

}




