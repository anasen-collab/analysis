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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeometryCal(void) {

  using namespace std;
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-7mQ2S1.root");//10MeV only
  //TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9mQ2S3.root");//in gas
  if ( !f1->IsOpen() ) {
    cout << "Error: Root File Does Not Exist\n";
    exit(EXIT_FAILURE);
  }
  Gains gains;
  gains.Save("saves/X3geometry");

  TCanvas *can = new TCanvas("can","can",800,600);
  BadDetectors bad;
  bad.count=0;
  cout << "DetNum\tFrontCh\tBackCh\tmin\tmax\n";
  for (Int_t DetNum=4; DetNum<ndets; DetNum++) {
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++) {
      for (Int_t BackChNum=0; BackChNum<4; BackChNum++) { 
	TH1F *hist = NULL;
	TString hname=Form("SX3ZposCal_%i_%i_%i",DetNum,FrontChNum,BackChNum);
	hist = (TH1F*)f1->Get(hname.Data());
	if (hist==NULL) {
	  cout << hname << " histogram does not exist\n";
	  bad.Add(DetNum,FrontChNum,BackChNum);
	  outfile << DetNum << "\t" << FrontChNum << "\t" << BackChNum << "\t" << 0 << "\t" << 0 << endl;
	  continue;
	}
	hist->Rebin(2);
	hist->Draw();
	Int_t Nbins = hist->GetNbinsX();
	Int_t maxbin = hist->GetMaximumBin();
	Double_t half_max = hist->GetBinContent(maxbin)/3.0;//seems to be a good point

	Double_t x[10];
	Double_t y[10];
	for (int i=0; i<10; i++) {
	  x[i] = -10;
	  y[i] = 0;
	}
	Int_t counter = 0;
	Double_t x_pos = -1 + 2./Nbins;

	for (int i=1; i<Nbins; i++) {
	  x_pos += 2./Nbins;
	  if (hist->GetBinContent(i)>half_max && counter==0) {//the first time the bin is greater than the half-max, this is the down content
	    x[0] = x_pos;
	    y[0] = hist->GetBinContent(i);
	    counter++;
	  }
	}

	x_pos = 1-2./Nbins;
	counter = 0;
	for (int i=Nbins-1; i>0; i--) {
	  x_pos -= 2./Nbins;
	  if (hist->GetBinContent(i)>half_max && counter==0) {//the first time the bin is greater than the half-max, this is the down content
	    x[1] = x_pos;
	    y[1] = hist->GetBinContent(i);
	    counter++;
	  }
	}

	//TGraph *graph = new TGraph(2,x,y);
	//graph->Draw("*same");
	TLine *line = new TLine(x[0],0,x[0],half_max*3);
	line->Draw("same");

	TLine *line2 = new TLine(x[1],0,x[1],half_max*3);
	line2->Draw("same");

	Int_t wide=6;
	Int_t prec=3;
	
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
