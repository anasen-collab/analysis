#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
//#include <TF1.h>
#include <TSpectrum.h>
#include <fstream>
#include <iostream>
//#include <exception>
//#include <TCutG.h>



void print_hist_qqq(void)
{

   using namespace std;
   
   string rootfile = "/data0/manasta/OrganizeRaw_files/run930_931_nospacer_copy2NoPcCal.root";
   TFile *f1 = new TFile(rootfile.c_str());
   if (!f1->IsOpen()){
     cout << "Root file: " << rootfile << " could not be opened.\n";
     exit(EXIT_FAILURE);
   }

   TCanvas **can = new TCanvas*[4];
   for(Int_t i=0; i<4;++i)
     {
       can[i]=new TCanvas(Form("can%i",i),Form("det%i",i),1200,1200);
       can[i]->Divide(2,15);
     }
  

   TTree *input_tree = (TTree*) f1->Get("MainTree");

   Int_t DetNum = 0;
   Int_t UpChNum = 0;
   Int_t PlotCount = 1;
   Int_t i=0;
   Long64_t entries = -1;
   
   TH1F* hist1 = NULL;
   //TH1F* hist2 = NULL;
   TH1F* hist3 = NULL;
   
   for(DetNum=0;DetNum<4;DetNum++)
     {
       PlotCount=1;
   for(UpChNum=0;UpChNum<15;UpChNum++){     

   can[i]->cd(PlotCount); PlotCount++;
   hist1 = new TH1F(Form("det%i_up%i",DetNum,UpChNum),Form("det%i_up%i",DetNum,UpChNum),100,0,30);
   input_tree->Project(Form("det%i_up%i",DetNum,UpChNum),"Si.Detector.EnergyUp_Cal",Form("Si.Detector.DetID==%i && Si.Detector.UpChNum==%i",DetNum,UpChNum));

   cout << hist1->GetEntries() << endl;
   if(hist1->GetEntries()!=0)
     {
   hist1->Draw();
     }
   can[i]->Update();

   /*hist2 = new TH1F(Form("det%i_down%i",DetNum,UpChNum),Form("det%i_down%i",DetNum,UpChNum),100,0,30);
   input_tree->Project(Form("det%i_down%i",DetNum,UpChNum),"Si.Detector.EnergyDown_Cal",Form("Si.Detector.DetID==%i && Si.Detector.DownChNum==%i",DetNum,UpChNum));
  

   cout << hist2->GetEntries() << endl;
   if(hist2->GetEntries()!=0)
     {
       hist2->SetLineColor(2);
       hist2->Draw("same");
     }
   can[i]->Update();
   */

   hist3 = new TH1F(Form("det%i_back%i",DetNum,UpChNum),Form("det%i_back%i",DetNum,UpChNum),100,0,30);
   input_tree->Project(Form("det%i_back%i",DetNum,UpChNum),"Si.Detector.EnergyBack_Cal",Form("Si.Detector.DetID==%i && Si.Detector.BackChNum==%i",DetNum,UpChNum));
  
   can[i]->cd(PlotCount); PlotCount++;
   cout << hist3->GetEntries() << endl;
   if(hist3->GetEntries()!=0)
     {
       hist3->Draw();
     }
   can[i]->Update();
   
   cout << "UpChNum " << UpChNum << " entries " << entries << endl;
   
   }
   can[i]->SaveAs(Form("det%i.png",DetNum));
   i++;
   
     }
   
   
}
