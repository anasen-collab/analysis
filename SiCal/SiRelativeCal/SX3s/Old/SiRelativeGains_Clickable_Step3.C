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
Double_t MyFit(TH2F* hist, TCanvas *can){
    
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

void SiRelativeGains_Clickable_Step3(void)
{
  using namespace std;
  
  //TFile *f1 = new TFile("run236out_Step2.root");
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-61m.root");//all proton scattering

  ofstream outfile;
  outfile.open("saves/X3RelativeGains012216_Step3_click.dat");

  ifstream infile;
  infile.open("../saves/X3RelativeGains_09182016_Slope1.dat"); //input file name
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

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t bad_back[288];
  Int_t count_bad = 0;

  for (Int_t DetNum=4; DetNum<28; DetNum++){
    for (Int_t FrontChNum=1; FrontChNum<4; FrontChNum++){
      Int_t BackChanNumber = 0;
      if ( DetNum==11 ){
	BackChanNumber = 1;
      }
      if (DetNum==22 && (FrontChNum==2 || FrontChNum==3) ){
	continue;
      }
      if (DetNum==25 && FrontChNum==2){
	continue;
      }
      
      TH2F *hist = NULL;
      hist = (TH2F*)f1->Get(Form("back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChanNumber));
      if (hist==NULL){
	cout << "Histo does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_back[count_bad] = FrontChNum;
	count_bad++;
	continue;
      }
      
      Double_t average_slope = MyFit(hist,can);	
      slope[DetNum-4][FrontChNum+4] = slope[DetNum-4][FrontChNum+4]*average_slope;
      slope[DetNum-4][FrontChNum+8] = slope[DetNum-4][FrontChNum+8]*average_slope;
	
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




