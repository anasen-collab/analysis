//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative calibration of Si gains
//
// Output file (e.g."Sipulser_2015Dec13.dat") has the following columns:
// MBID, CBID, ASICs_Channel, ZeroShift(offset), Voltage_per_Ch(slope)
//
// Usage: root -l AlphaCal.C++ (from the same directory).
//
// Edited by : John Parker , 2016Jan22
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
#include <TVector.h>
#include <TLine.h>
#include <TLegend.h>
//Methods
#include "SiAlphaCal.h"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace std;

void AlphaCal(void) {
  const Int_t npeaks = 3;
   
  //for gold at spacer zero, this is what the alpha energy should be
  //see lab book page 9
  Double_t energy[28];
  for (Int_t i=0; i<28; i++) {
    if (i==0 || i==3) {
      energy[i] = 12.03;
    }else if (i==1 || i==2) {
      energy[i] = 12.07;
    }else if (i>3 && i<16) {
      energy[i] = 12.3;
    }else if (i>15) {
      energy[i] = 12.84;
    }
  }

  //TFile *f1 = new TFile("/home2/parker/ANASEN/LSU/ParkerMain_root/run_alpha0_282_284_cal022716.root");
  //TFile *f1 = new TFile("/home2/parker/ANASEN/LSU/ParkerMain_root/run417_cal.root");
  //TFile *f1 = new TFile("/home/lighthall/anasen/root/main/run1255-7mQ2S3_geo_init.root");
  TFile *f1 = new TFile("/home/lighthall/anasen/root/main/run1255-61mQ2S3.root");
  
  if ( !f1->IsOpen() ) {
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
    
  TTree *MainTree = NULL;
  MainTree = (TTree*)f1->Get("MainTree");
  if (MainTree==NULL) {
    cout << "Tree does not exist\n";
    exit(EXIT_FAILURE);
  }
  Gains gains;
  Bool_t isrecal=kTRUE;
  if(isrecal) {
    gains.Load("saves/AlphaCal_170515.edit.dat");
    range=13;
  }
  else {
    gains.Load("saves/AlphaCal_init.dat");
  }
  gains.Save("saves/AlphaCal");
    
  TCanvas *can = new TCanvas("can","can",800,600);
  can->Divide(1,2);
  TH1F *hist = new TH1F("hist","hist",1024*2,0,range);
  TF1 *fit = new TF1("fit","pol1",0,range);
  TF1 *fit2 = new TF1("fit2","[0]*x",0,range);
  Double_t zeroshift = 0;
  Double_t MeVperCh = 0;
  Double_t gain = 0;
  TSpectrum *s = 0;
  TGraph *FitGraph = 0;
  TLegend *leg;

  Float_t Energies[npeaks] = {4.841,6.852,9.854};
  Float_t Energies1[npeaks-1];
  
  for (Int_t DetNum=0; DetNum<28; DetNum++) {
    if (DetNum==1 || DetNum==2) {
      Energies[0] = 4.841;
      Energies[1] = 6.852;
      Energies[2] = 9.852;
    }
    else if (DetNum>3 && DetNum<16) {
      Energies[0] = 4.837;
      Energies[1] = 6.849;
      Energies[2] = 9.846; 
    }
    else if (DetNum>15) {
      Energies[0] = 4.819;
      Energies[1] = 6.826;
      Energies[2] = 9.817; 
    }
    Energies1[0] = Energies[0];
    Energies1[1] = Energies[2];

    Float_t *average_slope;
    can->cd(1);	 
    MainTree->Draw("EnergyBack>>hist",Form("DetID==%i && (HitType==111 || HitType==11)",DetNum),"");
    hist->SetTitle(Form("Det%i",DetNum));  

    if(hist->GetEntries()==0) {
      printf("DetNum %2d: Histogram %s has zero entries.\n", DetNum,hist->GetTitle());
      outfile << DetNum << "\t"<< 0 << "\t" << 1 <<endl;
      continue;
    }
	      
    if(s!=0) {
      delete s;
    }

    if(isrecal)
      hist->GetXaxis()->SetRangeUser(4,11);
    else
      hist->GetXaxis()->SetRangeUser(1500,6000);
    can->Update();
    TSpectrum *s = new TSpectrum();
    Int_t nfound = s->Search(hist,2," ",0.05);//9 and 0.15
   
    if(nfound <(npeaks-1)) {
      printf("DetNum %2d: peaks = %d, ",DetNum,nfound);
      printf("Less than %d peaks found. Aborting.\n",npeaks);
      //outfile << DetNum << "\t"<< 0 << "\t" << 1 <<endl;
      continue;
    }

    if(nfound >npeaks) {
      printf("DetNum %2d: peaks = %d, ",DetNum,nfound);
      printf("greater than %d peaks found. Aborting.\n",npeaks);
      //outfile << DetNum << "\t"<< 0 << "\t" << 1 <<endl;
      continue;
    }
    
    //sort peaks in order of channle number
    Float_t *xpeaks = s->GetPositionX();
    Float_t Temp=0;    
    for(Int_t i=0;i<nfound;i++) {
      for(Int_t j=i;j<nfound;j++) {
	if (xpeaks[j] < xpeaks[i]) {
	  Temp = xpeaks[i];
	  xpeaks[i] = xpeaks[j];
	  xpeaks[j] = Temp;  
	}
      }
    }
    
    if (npeaks==1) {
      if (nfound != 1) {
	printf("DetNum %2d: peaks = %d",DetNum,nfound);
	cout << "No peaks found\n";
	//Double_t median = hist->GetMaximumBin()*0.00266666 + 0.8;
	Double_t median = hist->GetMaximumBin()*0.0012 + 0.8;
	gains.Add(DetNum,median,energy[DetNum]/median);
	TLine *line = new TLine(median,0,median,1000);
	line->Draw("same");
	can->Update();
	continue;
      }
      average_slope = s->GetPositionX();
      gains.Add(DetNum,average_slope[0],energy[DetNum]/average_slope[0]);
      outfile << DetNum << "\t" << 0 << "\t" << gains.old[DetNum] << endl;
      can->Update();
    }

    if(FitGraph!=0) {
      delete FitGraph;
    }

    if(nfound==npeaks) {
      FitGraph = new TGraph(nfound,xpeaks, &(Energies[0]));
    }
    else if(nfound==(npeaks-1)) {
      FitGraph = new TGraph(nfound,xpeaks, &(Energies1[0]));
    }
    can->cd(2);
    FitGraph->Draw("AP*");
    FitGraph->Fit("fit","q");
    zeroshift = fit->GetParameter(0);
    MeVperCh = fit->GetParameter(1);

    // cout << zeroshift << " " << MeVperCh << endl;;
    FitGraph->Fit("fit2","q");
    gain = fit2->GetParameter(0);
    fit2->SetLineColor(3);
    fit2->SetLineStyle(3);
    
    printf("DetNum %2d: peaks = %d, offset = %f\tslope = %f\tgain = %f (%5.2f%% diff)\n",
	   DetNum,nfound,zeroshift,MeVperCh,gain,((MeVperCh-gain)/gain)*100);
    leg = new TLegend(0.1,0.75,0.2,0.9);
    leg->AddEntry(FitGraph,"peaks","p");
    leg->AddEntry(fit,"linear fit (for reference)","l");
    leg->AddEntry(fit2,"scale fit","l");   
    leg->Draw();
    fit->Draw("same");
    fit2->Draw("same");

    can->Update();
    outfile << DetNum << "\t"<< 0 << "\t" << gain <<endl;
  }
  delete FitGraph;
  delete s;
}
