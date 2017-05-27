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
void PosCal(void) {
  using namespace std;
  Double_t x[99];
  Double_t y[99];
  string rootfile = "/home/lighthall/anasen/root/track/spacer0t.root";
  TFile *f1 = new TFile(rootfile.c_str());
  if (!f1->IsOpen()){
    cout << "Root file: " << rootfile << " could not be opened.\n";
    exit(EXIT_FAILURE);
  }

  TCanvas *can = new TCanvas("can","can",800,600);
  can->ToggleEventStatus(); 
  can->ToggleToolBar();
  
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
  outfile.open("saves/PCWireCal_170527_4.dat");
  outfile << "Wire\tSlope\tShift\n";

  TCutG *cut;
  for (Int_t DetNum=0; DetNum<24; DetNum++) {
    TH2F *hist = NULL;
    if(DetNum!=18) continue;
    hist = (TH2F*)f1->Get(Form("PCZ_vs_Zg%i",DetNum));
    if (hist==NULL){
      outfile << DetNum << "\t" << 1 << "\t" << 0 << endl;
      cout << "Warning: hist with wire number " << DetNum << " does not exist.\n";
	continue;
    }
    can->SetLogz();
    hist->Draw("colz");
    hist->GetXaxis()->SetRangeUser(-1.1,1.1);
    hist->GetYaxis()->SetRangeUser(5,27);
    
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
    //can->WaitPrimitive();

    
    printf(" Previous gain   = %f \t Slope  = %f \t New gain  = %f\n",
	   slope[DetNum],fun->GetParameter(1),slope[DetNum]*fun->GetParameter(1));
    printf(" Previous offset   = %f \t offset  = %f \t new offset  = %f\n",
	   offset[DetNum],fun->GetParameter(0),offset[DetNum]+fun->GetParameter(0));

    //slope[DetNum] *= fun->GetParameter(1);
    //offset[DetNum] += fun->GetParameter(0);

    //slope[DetNum] = fun->GetParameter(1);                      
    //offset[DetNum] = fun->GetParameter(0);
    
    outfile << DetNum << "\t"
	    << slope[DetNum]*fun->GetParameter(1) << "\t"
	    << offset[DetNum]+fun->GetParameter(0) << endl;

  }  
}
