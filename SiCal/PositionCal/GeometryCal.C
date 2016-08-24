//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//// Edited by : John Parker , 2016Jan24
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
#include <TLine.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void GeometryCal(void){

  using namespace std;

  TFile *f1 = new TFile("../../OrganizeRaw_root/run567_051116.root");
  if ( !f1->IsOpen() ){
    cout << "Error: Root File Does Not Exist\n";
    exit(EXIT_FAILURE);
  }
  ofstream outfile;
  outfile.open("X3geometry_051116.dat");

  TCanvas *can = new TCanvas("can","can",800,600);
  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t bad_back[288];
  Int_t count_bad = 0;

  for (Int_t DetNum=8; DetNum<9; DetNum++){
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      for (Int_t BackChNum=0; BackChNum<4; BackChNum++){

	TH1F *hist = NULL;
	hist = (TH1F*)f1->Get(Form("ZPos%i_%i_%i",DetNum,FrontChNum,BackChNum));
	if (hist==NULL){
	  cout << "Histo does not exist\n";
	  bad_det[count_bad] = DetNum;
	  bad_front[count_bad] = FrontChNum;
	  bad_back[count_bad] = BackChNum;
	  count_bad++;
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
	for (int i=0; i<10; i++){
	  x[i] = -10;
	  y[i] = 0;
	}
	Int_t counter = 0;
	Double_t x_pos = -1 + 2./Nbins;

	for (int i=1; i<Nbins; i++){
	  x_pos += 2./Nbins;
	  if (hist->GetBinContent(i)>half_max && counter==0){//the first time the bin is greater than the half-max, this is the down content
	    x[0] = x_pos;
	    y[0] = hist->GetBinContent(i);
	    counter++;
	  }
	}

	x_pos = 1-2./Nbins;
	counter = 0;
	for (int i=Nbins-1; i>0; i--){
	  x_pos -= 2./Nbins;
	  if (hist->GetBinContent(i)>half_max && counter==0){//the first time the bin is greater than the half-max, this is the down content
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
	
	cout << x[0] << "  " << x[1] << endl;

	can->Update();
	can->WaitPrimitive();

	outfile << DetNum << "\t" << FrontChNum << "\t" << BackChNum << "\t" << x[0] << "\t" << x[1] << endl;
      }
    }
  }
  for (int i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << "  " << bad_back[i] << endl;
  }

}




