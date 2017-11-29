#include <fstream>
#include <iostream>

#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TSpectrum.h>

void print_hist_sx3(void) {

   using namespace std;
   
   string rootfile = "/data0/manasta/OrganizeRaw_files/run930_offsetchecks3_10312016.root";
   TFile *f1 = new TFile(rootfile.c_str());
   if (!f1->IsOpen()){
     cout << "Root file: " << rootfile << " could not be opened.\n";
     exit(EXIT_FAILURE);
   }

   TCanvas **can = new TCanvas*[3];
   for(Int_t i=0; i<3;++i)
     {
       can[i]=new TCanvas(Form("can%i",i),Form("can%i",i),1200,1200);
       can[i]->Divide(2,2);
     
     }

   //TTree *input_tree = (TTree*) f1->Get("MainTree");

   Int_t DetNum = 0;
   //Int_t UpChNum = 0;
   Int_t PlotCount = 1;
   Int_t i=0;
   //Long64_t entries = -1;
   
   TH2F* hist1 = NULL;
   //TH1F* hist2 = NULL;
   //TH1F* hist3 = NULL;
   
   for(DetNum=4;DetNum<8;DetNum++)
     {
       PlotCount=1;
          
       hist1 = (TH2F*)f1->Get(Form("sx3offset_back_vs_front%i",DetNum));    
   
       if (hist1!=NULL)
	 {
	   can[i]->cd(PlotCount); PlotCount++;
	   //can[i]->cd(DetNum+1);
	   hist1->Draw("colz");
	   can[i]->Update();
	 }
     
       hist1 = (TH2F*)f1->Get(Form("sx3offset_back_vs_front%i",DetNum+4));    
   
       if (hist1!=NULL)
       	 {
	   can[i]->cd(PlotCount); PlotCount++;
       	   //can[i]->cd(DetNum+1);
       	   hist1->Draw("colz");
       	   can[i]->Update();
	 }

       hist1 = (TH2F*)f1->Get(Form("sx3offset_back_vs_front%i",DetNum+8));    
   
       if (hist1!=NULL)
       	 {
	   can[i]->cd(PlotCount); PlotCount++;
       	   //can[i]->cd(DetNum+1);
       	   hist1->Draw("colz");
       	   can[i]->Update();
	 }
       
   
       i++;
     
     }

     

}
