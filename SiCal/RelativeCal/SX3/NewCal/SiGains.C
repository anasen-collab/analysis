//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Relative calibration of Si gains
////
////
//// Edited by : John Parker , 2016Jan22
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
#include <TVector.h>
#include <TLine.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace std;

void Step1(TTree *tree, TCanvas *can, TH1F *hist, Double_t array[][12]);
void Step2(TTree *tree, TCanvas *can, TH1F *hist, Double_t array[][12]);
void Step3(TTree *tree, TCanvas *can, TH2F *down_vs_up, Double_t array[][12]);

void SiGains(void){

  TFile *f1 = new TFile("/home2/parker/ANASEN/LSU/ParkerMain_root/run410_419_cal030816b.root");
  //TFile *f1 = new TFile("/home2/parker/ANASEN/LSU/ParkerMain_root/run417_cal.root");


  TCanvas *can = new TCanvas("can","can",800,600);

  ofstream outfile;
  outfile.open("X3RelativeGains030816b.dat");

  ifstream infile;
  infile.open("/home2/parker/ANASEN/LSU/CalParamFiles/X3RelativeGains030816.dat");
  Int_t det=0,ch=0;
  Double_t dummy_slope = 0;
  Double_t slope[24][12];
  if (infile.is_open()){
    while (!infile.eof()){
      infile >> det >> ch >> dummy_slope;
      slope[det-4][ch] = dummy_slope;
    }
  }else{
    cout << "Infile not opened\n";
    exit(EXIT_FAILURE);
  }
  infile.close();

  Double_t array[24][12];
  Double_t array_old[24][12];

  for (Int_t i=0; i<24; i++){
    for (Int_t j=0; j<12; j++){
      array[i][j] = 1;
    }
  }
  TTree *tree = NULL;
  tree = (TTree*)f1->Get("MainTree");
  if (tree==NULL){
    cout << "Tree does not exist\n";
    exit(EXIT_FAILURE);
  }
  TH1F *hist = new TH1F("hist","hist",800,0.5,1.5);
  TH2F *down_vs_up = new TH2F("down_vs_up","down_vs_up",300,0,1,300,0,1);

  Step3(tree, can, down_vs_up, array);
  for ( Int_t iterations=0; iterations<2; iterations++ ){
    for (Int_t i=0; i<24; i++){
      for (Int_t j=0; j<12; j++){
	array_old[i][j] = array[i][j];
      }
    }
    Step1(tree, can, hist, array);
    Step2(tree, can, hist, array);
    cout << "Diff after iteration " << iterations << endl;
    for (Int_t i=0; i<24; i++){
      for ( Int_t j=0; j<12; j++ ){
	cout << i+4 << "  " << j << "  " << array_old[i][j] - array[i][j] << endl;
      }
    }
    cout << endl << endl;

  }
  
  
  for (Int_t DetNum=4; DetNum<28; DetNum++){
    for (Int_t ChNum = 0; ChNum<12; ChNum++){
      slope[DetNum-4][ChNum] = slope[DetNum-4][ChNum]*array[DetNum-4][ChNum];
    }
  }
  for (Int_t DetNum=4; DetNum<28; DetNum++){
    for (Int_t ChNum = 0; ChNum<12; ChNum++){
      outfile << DetNum << "\t" << ChNum << "\t" << slope[DetNum-4][ChNum] << endl;
    }
  }
}

void Step1(TTree *tree, TCanvas *can, TH1F *hist, Double_t array[][12]){
  cout << "Step 1" << endl;
  Int_t BackChNum;
  Float_t *average_slope;

  for (Int_t DetNum=4; DetNum<28; DetNum++){
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      BackChNum = 0;
      if ( DetNum==7 || DetNum==9 ){
	continue;
      }
      if ( DetNum==11 || DetNum==18 || DetNum==24 ){
	//BackChNum = 0;
      }
      if ( DetNum==22 && (FrontChNum==2 || FrontChNum==3) ){
	continue;
      }
      if ( DetNum==25 && FrontChNum==2 ){
	continue;
      }
      tree->Draw(Form("EnergyBack*%f/(EnergyDown*%f+EnergyUp*%f)>>hist",array[DetNum-4][BackChNum],array[DetNum-4][FrontChNum+4],array[DetNum-4][FrontChNum+8]),Form("DetID==%i && UpChNum==%i && BackChNum==%i && HitType==111",DetNum,FrontChNum,BackChNum),"");
      hist->SetTitle(Form("Det%i FrontCh%i BackCh%i",DetNum,FrontChNum,BackChNum));  

      TSpectrum *s = new TSpectrum(2);
      Int_t nfound = s->Search(hist,0.02," ",0.8);//9 and 0.15
      //hist->Fit("fun");
      if (nfound != 1){
	cout << DetNum << "  " << FrontChNum << "  " << BackChNum << "  " << nfound << endl;
	cout << "No peaks found\n";
	//Double_t median = hist->GetMaximumBin()*0.00266666 + 0.8;
	Double_t median = hist->GetMaximumBin()*0.00125 + 0.5;
	cout << median << endl;
	array[DetNum-4][FrontChNum+4] *= median;
	array[DetNum-4][FrontChNum+8] *= median;
	cout << array[DetNum-4][FrontChNum+4] << "  " << array[DetNum-4][FrontChNum+8] << endl;
 
	TLine *line = new TLine(median,0,median,1000);
	line->Draw("same");
	can->Update();
	continue;
      }
      average_slope = s->GetPositionX();
      array[DetNum-4][FrontChNum+4] *= average_slope[0];
      array[DetNum-4][FrontChNum+8] *= average_slope[0];

      //      cout << average_slope[0] << endl;
      cout << array[DetNum-4][FrontChNum+4] << "  " << array[DetNum-4][FrontChNum+8] << endl;

      can->Update();
    }
  }

}

void Step2(TTree *tree, TCanvas *can, TH1F *hist, Double_t array[][12]){
  cout << "Step 2" << endl;
  Float_t *average_slope;

  for (Int_t DetNum=4; DetNum<28; DetNum++){
    for (Int_t BackChNum=1; BackChNum<4; BackChNum++){
      Int_t FrontChNum = 0;
      if ( DetNum==7 || DetNum==9 || DetNum==11 || DetNum==26){
	continue;
      }
      if ( DetNum==18 && BackChNum==2 ){
	continue;
      }
      if ( DetNum==24 && (BackChNum==1 || BackChNum==2) ){
	continue;
      }
      if ((DetNum==19 || DetNum==27) && (BackChNum==2 || BackChNum==3) ){
	continue;
      }
      tree->Draw(Form("(EnergyDown*%f+EnergyUp*%f)/(EnergyBack*%f)>>hist",array[DetNum-4][FrontChNum+4],array[DetNum-4][FrontChNum+8],array[DetNum-4][BackChNum]),Form("DetID==%i && UpChNum==%i && BackChNum==%i && HitType==111",DetNum,FrontChNum,BackChNum),"");
      hist->SetTitle(Form("Det%i FrontCh%i BackCh%i",DetNum,FrontChNum,BackChNum));  

      TSpectrum *s = new TSpectrum(2);
      Int_t nfound = s->Search(hist,0.02," ",0.8);//9 and 0.15
      //hist->Fit("fun");
      if (nfound != 1){
	cout << DetNum << "  " << FrontChNum << "  " << BackChNum << "  " << nfound << endl;
	cout << "No peaks found\n";

	//Double_t median = hist->GetMaximumBin()*0.00266666 + 0.8;
	Double_t median = hist->GetMaximumBin()*0.00125 + 0.5;
	cout << median << endl;
	array[DetNum-4][BackChNum] *= median;
	cout << array[DetNum-4][BackChNum] << endl;
 
	TLine *line = new TLine(median,0,median,1000);
	line->Draw("same");
	can->Update();
	continue;
	can->Update();
	continue;
      }
      average_slope = s->GetPositionX();
      array[DetNum-4][BackChNum] *= average_slope[0];

      //      cout << average_slope[0] << endl;
      cout << array[DetNum-4][BackChNum] << endl;
      can->Update();
    }
  }

}
void Step3(TTree *tree, TCanvas *can, TH2F *down_vs_up, Double_t array[][12]){
  cout << "Step 3" << endl;
  Float_t slope[2];
  slope[0] = 1;
  slope[1] = 1;

  for (Int_t DetNum=4; DetNum<28; DetNum++){
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      slope[0] = 1;
      slope[1] = 1;
      Int_t BackChNum = 1;
      if ( DetNum==7 || DetNum==9){
	continue;
      }
      if (DetNum==22 && (FrontChNum==2 || FrontChNum==3) ){
	continue;
      }
      if (DetNum==25 && FrontChNum==2 ){
	continue;
      }
      if ( DetNum==11 || DetNum==24 || DetNum==26){
	BackChNum = 0;
      }
      TF1 *fun = new TF1("fun","1-x",0,1);
      TF1 *fun1 = new TF1("fun1","1.1-x",0,1);
      TF1 *fun2 = new TF1("fun2","1.2-x",0,1);
      TF1 *fun3 = new TF1("fun3","0.9-x",0,1);
      TF1 *fun4 = new TF1("fun4","0.8-x",0,1);
      TF1 *fun5 = new TF1("fun5","0.7-x",0,1);
      TF1 *fun6 = new TF1("fun6","0.6-x",0,1);
      TF1 *fun7 = new TF1("fun7","1.02-x",0,1);
      TF1 *fun8 = new TF1("fun8","0.98-x",0,1);
      fun->SetLineWidth(1);
      fun1->SetLineWidth(1);
      fun2->SetLineWidth(1);
      fun3->SetLineWidth(1);
      fun4->SetLineWidth(1);
      fun5->SetLineWidth(1);
      fun6->SetLineWidth(1);
      fun7->SetLineWidth(1);
      fun8->SetLineWidth(1);

      cout << "Enter correction factor to make data line up with red line\n";
      cout << "When satisfied, type '900'\n";
      while ( slope[1]!=900 ){
	slope[0]=slope[1];
	tree->Draw(Form("EnergyDown*%f/EnergyBack:EnergyUp/EnergyBack>>down_vs_up",slope[1]),Form("DetID==%i && UpChNum==%i && BackChNum==%i && HitType==111",DetNum,FrontChNum,BackChNum),"colz");
	fun->Draw("same");
	fun1->Draw("same");
	fun2->Draw("same");
	fun3->Draw("same");
	fun4->Draw("same");
	fun5->Draw("same");
	fun6->Draw("same");
	down_vs_up->SetTitle(Form("Det%i FrontCh%i BackChNum%i",DetNum,FrontChNum, BackChNum));
	can->Update();
	cin >> slope[1];
	if (slope[1]==900){
	  array[DetNum-4][FrontChNum+4] *= slope[0];
	}
      }

    }
  }

}
