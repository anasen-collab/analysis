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
#include <iostream>
//#include <exception>

void print_down_vs_up(void) {

   using namespace std;
   
   //string rootfile = "/data0/manasta/OrganizeRaw_files/run924_testnewercals_pulseroffsets_11112016.root";
   string rootfile = "/data0/manasta/OrganizeRaw_files/run924_step2_fixBack_otherdet_11282016.root";
   //string rootfile = "/data0/manasta/OrganizeRaw_files/24Mg/run1233_1236_allnewercals_11172016.root";
   TFile *f1 = new TFile(rootfile.c_str());
   if (!f1->IsOpen()){
     cout << "Root file: " << rootfile << " could not be opened.\n";
     exit(EXIT_FAILURE);
   }

   TCanvas **can = new TCanvas*[24];
   for(Int_t i=0; i<24;++i)
     {
       can[i]=new TCanvas(Form("can%i",i),Form("det%i",i+4),1200,1200);
       can[i]->Divide(2,2);
     }
  
   Int_t DetNum = 0;
   Int_t BackChNum = 0;
   Int_t PlotCount = 1;
   Int_t i=0;
   //Long64_t entries = -1;
   
   TH2F* hist1 = NULL;
   //TH1F* hist2 = NULL;
   //TH1F* hist3 = NULL;
   
   for(DetNum=4;DetNum<28;DetNum++)
     {
       PlotCount=1;
   for(BackChNum=0;BackChNum<4;BackChNum++){     


   TF1 *fun = new TF1("fun","1-x",0,1);
    
   hist1=(TH2F*)f1->Get(Form("down_vs_up_divideBack%i_front%i",DetNum,BackChNum));
  
   if (hist1!=NULL)
     {
   can[i]->cd(PlotCount); PlotCount++;    
   hist1->Draw("colz");
   fun->Draw("same");
   can[i]->Update();
     }
  
   }
   can[i]->SaveAs(Form("det%i.png",DetNum));
   i++;
   
     }
   
   
}
