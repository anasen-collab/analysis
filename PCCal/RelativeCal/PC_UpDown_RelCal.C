//////////////////////////////////////////////////////////////////////////
// Author : Nabin Rijal, 20170428
//////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////
void PC_UpDown_RelCal(void) {
  using namespace std;

  Double_t x[10];
  Double_t y[10];  

  string rootfile = "/home/lighthall/anasen/root/main/run1036_rel.root";
  TFile *f1 = new TFile(rootfile.c_str());
  if (!f1->IsOpen()) {
    cout << "Root file: " << rootfile << " could not be opened.\n";
    exit(EXIT_FAILURE);
  } 
  ofstream outfile;  
  outfile.open("saves/PC_UD_RelCal_2018_output.dat");
  outfile << "Wire\tSlope\tShift\n";

  TCanvas *pad = new TCanvas("pad","pad",800,600);
  if(!(pad->GetShowEventStatus()))pad->ToggleEventStatus();
  if(!(pad->GetShowToolBar()))pad->ToggleToolBar();
  TCutG *cut;
  for (Int_t WireNum=0; WireNum<24; WireNum++){
    TH2F *hist = NULL;  
    hist = (TH2F*)f1->Get(Form("PC_Down_vs_Up_BeforeCal_Wire%i",WireNum));
    if (hist==NULL){
      outfile << WireNum << "\t" << 1 << "\t" << 0 << endl;
      cout << "Warning: hist with wire number " << WireNum << " does not exist.\n";
	continue;
    } 
    hist->Draw("colz");	
    /*    cut = (TCutG*)pad->WaitPrimitive("CUTG");
	
    for(int n=0;n<cut->GetN()-2;n++){
      cut->GetPoint(n,x[n],y[n]);
      //cout << x[n] << "\t" << y[n] << endl;
    }	
    TGraph *graph = new TGraph(cut->GetN()-2,x,y);
    hist->Draw("colz");
    graph->Draw("*same");
    */
    
    TF1 *f2 = new TF1("f2","[0] + [1]*x",0,1);
    hist->Fit("f2");

    Float_t Slope = f2->GetParameter(1);
    Float_t Shift = f2->GetParameter(0);

    cout<< "Slope = "<<Slope<<"   Shift = "<<Shift<<endl;

    pad->Update();
    //    pad->WaitPrimitive();          
    outfile << WireNum << "\t" << Slope << "\t" << Shift << endl;    
  } 
}
