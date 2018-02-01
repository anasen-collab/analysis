//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Geometry calibration for SX3
// 
// Usage: root -l GeomtryCal.C+
//
// Edited by : John Parker , 2016Jan24
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//C++
#include <fstream>
#include <exception>
#include <iomanip>
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
#include <TLine.h>
//Methods
#include "SiGeometry.h"

#define ratio 3.0
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeometryCal(void) {

  using namespace std;
  //TFile *f1 = new TFile("/home/lighthall/anasen/root/main/run1255-7mQ2S3_geo_init.root");//10MeV only
  TFile *f1 = new TFile("/home/lighthall/anasen/root/main/spacer0.root");//in gas
  if ( !f1->IsOpen() ) {
    cout << "Error: Root File Does Not Exist\n";
    exit(EXIT_FAILURE);
  }

  Double_t x[2]={0};
  Double_t y[2]={0};
  Int_t wide=6;
  Int_t prec=3;

  Gains gains;
  gains.Save("saves/X3geometry");

  TCanvas *can = new TCanvas("can","can",800,600);
  BadDetectors bad;
  bad.count=0;
  cout << "DetNum\tFrontCh\tBackCh\tmin\tmax\n";
  for (Int_t DetNum=4; DetNum<ndets+4; DetNum++) {
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++) {
      for (Int_t BackChNum=0; BackChNum<4; BackChNum++) {
	//if(!(DetNum==8)) continue;
	TH1F *hist = NULL;
	TString hname=Form("SX3ZposCal_%i_%i_%i",DetNum,FrontChNum,BackChNum);
	hist = (TH1F*)f1->Get(hname.Data());
	if (hist==NULL) {
	  cout << hname << " histogram does not exist\n";
	  bad.Add(DetNum,FrontChNum,BackChNum);
	  x[0] = -1+0.5*BackChNum;
	  x[1] = x[0]+0.5;
	  outfile << DetNum << "\t" << FrontChNum << "\t" << BackChNum << "\t"
		  << right << fixed <<  setw(wide) << setprecision(prec) << x[0] << "\t"
		  << right << fixed <<  setw(wide) << setprecision(prec) << x[1] << endl;
	  continue;
	}
	//hist->Rebin(2);
	hist->Draw();
	Int_t Nbins = hist->GetNbinsX();
	Int_t maxbin = hist->GetMaximumBin();
	Double_t half_max = hist->GetBinContent(maxbin)/ratio;//seems to be a good point

	Int_t counter = 0;
	
	for (int i=2; i<Nbins; i++) {
	  if (hist->GetBinContent(i)>half_max && counter==0) {//the first time the bin is greater than the half-max, this is the down content
	    x[0] = hist->GetBinLowEdge(i);
	    y[0] = hist->GetBinContent(i);
	    counter++;
	  }
	}

	counter = 0;
	for (int i=Nbins-1; i>0; i--) {
	  if (hist->GetBinContent(i)>half_max && counter==0) {//the first time the bin is greater than the half-max, this is the down content
	    x[1] = hist->GetBinLowEdge(i)+hist->GetBinWidth(i);
	    y[1] = hist->GetBinContent(i);
	    counter++;
	  }
	}

	Float_t width=x[1]-x[0];
	if(width<0.25||width>0.5)
	  printf("ERROR");

	TLine *line = new TLine(x[0],0,x[0],half_max*ratio);
	line->SetLineColor(2);
	line->Draw("same");

	TLine *line2 = new TLine(x[1],0,x[1],half_max*ratio);
	line2->SetLineColor(2);
	line2->Draw("same");

	TLine *line3 = new TLine(x[0],half_max,x[1],half_max);
	line3->SetLineColor(3);
	line3->Draw("same");

	cout << DetNum << "\t" << FrontChNum << "\t" << BackChNum << "\t"
	     << right << fixed <<  setw(wide) << setprecision(prec) << x[0] << "  "
	     << right << fixed <<  setw(wide) << setprecision(prec) << x[1] << endl;
	
	can->Update();
	//can->WaitPrimitive();

	outfile << DetNum << "\t" << FrontChNum << "\t" << BackChNum << "\t"
		<< right << fixed <<  setw(wide) << setprecision(prec) << x[0] << "\t"
		<< right << fixed <<  setw(wide) << setprecision(prec) << x[1] << endl;
      }
    }
  }
  outfile.close();
  bad.Print();
}
