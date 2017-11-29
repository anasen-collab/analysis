#include <fstream>
#include <iostream>

#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TSpectrum.h>

void print_histos(void) {

   using namespace std;
   
   string rootfile = "/data0/manasta/OrganizeRaw_files/run930_pcposcal_step1_protoncutwwithSiRonQQQ_01272017.root";
   //string rootfile = "/data0/manasta/OrganizeRaw_files/24Mg/run1233_1236_allnewercals_11172016.root";
   TFile *f1 = new TFile(rootfile.c_str());
   if (!f1->IsOpen()){
     cout << "Root file: " << rootfile << " could not be opened.\n";
     exit(EXIT_FAILURE);
   }

   TCanvas *can1 = new TCanvas("can1","can1",1200,1200);
   TCanvas *can2 = new TCanvas("can2","can2",1200,1200);
   TCanvas *can3 = new TCanvas("can3","can3",1200,1200);
   TCanvas *can4 = new TCanvas("can4","can4",1200,1200);
   // TCanvas *can5 = new TCanvas("can5","can5",1200,1200);
   // TCanvas *can6 = new TCanvas("can6","can6",1200,1200);

   can1->Divide(2,3);
   can2->Divide(2,3);
   can3->Divide(2,3);
   can4->Divide(2,3);
   //can5->Divide(2,2);
   //can6->Divide(2,2);
   

  
   
     for (Int_t DetNum=0; DetNum<6; DetNum++){
      
       TH2F *hist1 = NULL;
       
       //TF1 *fun = new TF1("fun","x",0,16384);
       
       //hist1 = (TH2F*)f1->Get(Form("back_vs_front%i",DetNum+4));
       //hist1 = (TH2F*)f1->Get(Form("sx3offset_back_vs_front%i",DetNum+4));
       //hist1 = (TH2F*)f1->Get(Form("offset_check%i",DetNum+4));
       // hist1 = (TH2F*)f1->Get(Form("Si.Hit.Z_linear_vs_Si.Hit.EnergyBack_det%i",DetNum+4));
       //PCZoffset_vs_PCEnergy_%i
        hist1 = (TH2F*)f1->Get(Form("PCZoffset_vs_PCEnergy_%i",DetNum));
	 if (hist1!=NULL){
	   can1->cd(DetNum+1);
	   //hist1->SetXTitle("EnergyFront_Rel");
	   //hist1->SetYTitle("EBack_Rel-EFront_Rel");
	   //hist1->SetXTitle("EnergyBack_Rel");
	   //hist1->SetYTitle("EDown_Rel-EUp_Rel/EBack_Rel");
	   hist1->Draw("colz");
	   //fun->Draw("same");
	   can1->Update();
	 } 

	 
	 //hist1 = (TH2F*)f1->Get(Form("back_vs_front%i",DetNum+10));
	 //hist1 = (TH2F*)f1->Get(Form("sx3offset_back_vs_front%i",DetNum+10));
	 //hist1 = (TH2F*)f1->Get(Form("offset_check%i",DetNum+8));
	 //hist1 = (TH2F*)f1->Get(Form("Si.Hit.Z_linear_vs_Si.Hit.EnergyBack_det%i",DetNum+10));
	 hist1 = (TH2F*)f1->Get(Form("PCZoffset_vs_PCEnergy_%i",DetNum+6));
	 if (hist1!=NULL){
	   can2->cd(DetNum+1);
	   //hist1->SetXTitle("EnergyFront_Rel");
	   //hist1->SetYTitle("EBack_Rel-EFront_Rel");
	   //hist1->SetXTitle("EnergyBack_Rel");
	   //hist1->SetYTitle("EDown_Rel-EUp_Rel/EBack_Rel");
	   hist1->Draw("colz");
	   //fun->Draw("same");
	   can2->Update();
	 } 
	 
	 //hist1 = (TH2F*)f1->Get(Form("back_vs_front%i",DetNum+16));
	 //hist1 = (TH2F*)f1->Get(Form("sx3offset_back_vs_front%i",DetNum+16));
	 //hist1 = (TH2F*)f1->Get(Form("offset_check%i",DetNum+12));
	 //hist1 = (TH2F*)f1->Get(Form("Si.Hit.Z_linear_vs_Si.Hit.EnergyBack_det%i",DetNum+16));
	 hist1 = (TH2F*)f1->Get(Form("PCZoffset_vs_PCEnergy_%i",DetNum+12));
	 if (hist1!=NULL){
	   can3->cd(DetNum+1);
	   //hist1->SetXTitle("EnergyFront_Rel");
	   //hist1->SetYTitle("EBack_Rel-EFront_Rel");
	   //hist1->SetXTitle("EnergyBack_Rel");
	   //hist1->SetYTitle("EDown_Rel-EUp_Rel/EBack_Rel");
	   hist1->Draw("colz");
	   //fun->Draw("same");
	   can3->Update();
	 } 
	
	 
	 //hist1 = (TH2F*)f1->Get(Form("back_vs_front%i",DetNum+22));
	 //hist1 = (TH2F*)f1->Get(Form("sx3offset_back_vs_front%i",DetNum+22));
	 //hist1 = (TH2F*)f1->Get(Form("offset_check%i",DetNum+16));
	 //hist1 = (TH2F*)f1->Get(Form("Si.Hit.Z_linear_vs_Si.Hit.EnergyBack_det%i",DetNum+22));
	 hist1 = (TH2F*)f1->Get(Form("PCZoffset_vs_PCEnergy_%i",DetNum+18));
	 if (hist1!=NULL){
	   can4->cd(DetNum+1);
	   //hist1->SetXTitle("EnergyFront_Rel");
	   //hist1->SetYTitle("EBack_Rel-EFront_Rel");
	   //hist1->SetXTitle("EnergyBack_Rel");
	   //hist1->SetYTitle("EDown_Rel-EUp_Rel/EBack_Rel");
	   hist1->Draw("colz");
	   //fun->Draw("same");
	   can4->Update();
	 } 

	 /*
	 //hist1 = (TH2F*)f1->Get(Form("back_vs_front%i",DetNum+20));
	 hist1 = (TH2F*)f1->Get(Form("sx3offset_back_vs_front_normback%i",DetNum+20));
	 //hist1 = (TH2F*)f1->Get(Form("offset_check%i",DetNum+20));
	 if (hist1!=NULL){
	   can5->cd(DetNum+1);
	   //hist1->SetXTitle("EnergyFront_Rel");
	   //hist1->SetYTitle("EBack_Rel-EFront_Rel");
	   //hist1->SetXTitle("EnergyBack_Rel");
	   //hist1->SetYTitle("EDown_Rel-EUp_Rel/EBack_Rel");
	   hist1->Draw("colz");
	   //fun->Draw("same");
	   can5->Update();
	 } 
	 

	 //hist1 = (TH2F*)f1->Get(Form("back_vs_front%i",DetNum+24));
	 hist1 = (TH2F*)f1->Get(Form("sx3offset_back_vs_front_normback%i",DetNum+24));
	 //hist1 = (TH2F*)f1->Get(Form("offset_check%i",DetNum+24));
	 if (hist1!=NULL){
	   can6->cd(DetNum+1);
	   //hist1->SetXTitle("EnergyFront_Rel");
	   //hist1->SetYTitle("EBack_Rel-EFront_Rel");
	   //hist1->SetXTitle("EnergyBack_Rel");
	   //hist1->SetYTitle("EDown_Rel-EUp_Rel/EBack_Rel");
	   hist1->Draw("colz");
	   //fun->Draw("same");
	   can6->Update();
	   } 
	 */
	  
	 can1->SaveAs("pcposcal_run930_pcposcal_step1_protoncutwwithSiRonQQQ_01272017_1.png");
	 can2->SaveAs("pcposcal_run930_pcposcal_step1_protoncutwwithSiRonQQQ_01272017_2.png");
	 can3->SaveAs("pcposcal_run930_pcposcal_step1_protoncutwwithSiRonQQQ_01272017_3.png");
	 can4->SaveAs("pcposcal_run930_pcposcal_step1_protoncutwwithSiRonQQQ_01272017_4.png");

	 
       }
     
   
    
   
   
}
