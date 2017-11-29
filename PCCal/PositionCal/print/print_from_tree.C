#include <fstream>
#include <iostream>
//#include <exception>

#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
//#include <TF1.h>
#include <TSpectrum.h>

void print_from_tree(void) {

   using namespace std;
   
   string rootfile = "/data0/manasta/OrganizeRaw_files/run924_testnewercals_pulseroffsets_11112016.root";
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
  

   TTree *input_tree = (TTree*) f1->Get("MainTree");

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

   can[i]->cd(PlotCount); PlotCount++;
   //hist1 = new TH1F(Form("det%i_back%i",DetNum,BackChNum),Form("det%i_back%i",DetNum,BackChNum),100,100,130);
   hist1 = new TH2F(Form("det%i_up%i_back%i",DetNum,BackChNum,1),Form("det%i_up%i_back%i",DetNum,BackChNum,1),512,0,16384,200,-1,1);
   //input_tree->Project(Form("det%i_up%i",DetNum,UpChNum),"Si.Hit.HitType",Form("Si.Detector.DetID==%i && Si.Detector.UpChNum==%i",DetNum,UpChNum));

   input_tree->Project(Form("det%i_up%i_back%i",DetNum,BackChNum,1),"((Si.Detector.EDown_Rel-Si.Detector.EUp_Rel)/Si.Detector.EBack_Rel):Si.Detector.EBack_Rel",
		       Form("Si.Detector.DetID==%i && Si.Detector.HitType==111 && Si.Detector.UpChNum==%i && Si.Detector.BackChNum==%i",DetNum,BackChNum,1));


   //MyFill(Form("offset_check%i",Si.det_obj.DetID),512,0,16384,Si.det_obj.EBack_Rel[0],200,-1,1,((Si.det_obj.EDown_Rel[0]-Si.det_obj.EUp_Rel[0])/Si.det_obj.EBack_Rel[0]));

   cout << hist1->GetEntries() << endl;
   if(hist1->GetEntries()!=0)
     {
   hist1->Draw("colz");
     }
   can[i]->Update();

   /* hist2 = new TH1F(Form("det%i_down%i",DetNum,UpChNum),Form("det%i_down%i",DetNum,UpChNum),100,100,130);
   input_tree->Project(Form("det%i_down%i",DetNum,UpChNum),"Si.Hit.HitType",Form("Si.Detector.DetID==%i && Si.Detector.DownChNum==%i",DetNum,UpChNum));
  

   cout << hist2->GetEntries() << endl;
   if(hist2->GetEntries()!=0)
     {
       hist2->SetLineColor(2);
       hist2->Draw("same");
     }
   can[i]->Update();


   hist3 = new TH1F(Form("det%i_back%i",DetNum,UpChNum),Form("det%i_back%i",DetNum,UpChNum),100,100,130);
   input_tree->Project(Form("det%i_back%i",DetNum,UpChNum),"Si.Hit.HitType",Form("Si.Detector.DetID==%i && Si.Detector.BackChNum==%i",DetNum,UpChNum));
  
   can[i]->cd(PlotCount); PlotCount++;
   cout << hist3->GetEntries() << endl;
   if(hist3->GetEntries()!=0)
     {
       hist3->Draw();
     }
   can[i]->Update();
   
   cout << "UpChNum " << UpChNum << " entries " << entries << endl;
   */

   }
   can[i]->SaveAs(Form("det%i.png",DetNum));
   i++;
   
     }
   
   
}
