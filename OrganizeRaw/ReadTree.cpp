//Written by John Parker--28 Apr, 2016
//This program 'reads' the organize.root file event by event
//MainTree->Show(2) doesn't work because the tree is defined from vectors and root needs a dictionary to read it
//This program works by looping over all the entries in your tree and stopping and writing out information if a condition is met. It will continue once you type anything and hit enter (if you type quit, it will exit)
//I use the flags that I set in organize.cpp here...
//Usage:
//  rootcint -f organize_dictionary.cxx -c ../Include/organizetree.h LinkDef.h
//  g++ -o ReadTree organize_dictionary.cxx ReadTree.cpp `root-config --cflags --glibs`
//  ./ReadTree

//C and C++ libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>

//ROOT libraries
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TObjArray.h>

using namespace std;

#define ReadTrack

//Detector Libraries
#include "/home/manasta/Desktop/parker_codes/Include/organizetree.h"

bool Track::Si_sort_method(struct Silicon_Event a,struct Silicon_Event b){
  if(a.SiEnergy > b.SiEnergy)
    return 1;
  return 0;
};

bool Track::PC_sort_method(struct PropCounter_Event a,struct PropCounter_Event b){
  if(a.PCEnergy > b.PCEnergy)
    return 1;
  return 0;
};

int main(int argc, char *argv[]){
  TApplication *myapp=new TApplication("myapp",0,0); //Don't know what this does, but libraries won't load without it
  //Long64_t entry = 2;

  int entry = atoi(argv[1]);

  SiHit Si;
  Track Tr;
  Si.ReadDet = 0;
  Si.ReadHit = 0;
  Tr.ReadTrEvent = 0;
  //Si.detsort = 0;

  //PCHit PC;

  TFile *inputFile = new TFile("/home/manasta/Desktop/parker_codes/OrganizeRaw_files/run410_429_050416.root");
  TTree *tree = (TTree*) inputFile->Get("MainTree");
  tree->SetBranchAddress("Si.NSiHits",&Si.NSiHits);
  tree->SetBranchAddress("Si.Detector",&Si.ReadDet);
  tree->SetBranchAddress("Si.Hit",&Si.ReadHit);
  tree->SetBranchAddress("Tr.TrackEvent",&Tr.ReadTrEvent);

  Long64_t nentries = tree->GetEntries();
  for (Int_t global_evt=0; global_evt<nentries; global_evt++){
    Int_t status = tree->GetEvent(global_evt);

    //cout << Si.ReadHit->at(0).HitType << endl;
    //if ( !(Si.ReadDet->at(0).DetID == 4) )
    //continue;
    //if ( Tr.ReadTrEvent->at(0).HitType != 12 )
    //continue;

    if (global_evt != entry ){
      continue;
    }

    cout << "Event Number: " << global_evt << endl;
    cout << "\nNumber of Detectors Fired: " << Si.NSiHits << endl;
    cout << "------------------------------------------------------------\n";
    cout << "Data Organized By Detector\n";
    for (Int_t i=0; i<Si.ReadDet->size(); i++){
      Si.det_place_holder = Si.ReadDet->at(i);
      cout << "DetID: " << Si.det_place_holder.DetID << endl;
      cout << "UpMult: " << Si.det_place_holder.UpMult << endl;
      cout << "DownMult: " << Si.det_place_holder.DownMult << endl;
      cout << "BackMult: " << Si.det_place_holder.BackMult << endl;
      cout << "HitType: " << Si.det_place_holder.HitType << endl;
      
      cout << "Up:\n";
      for (Int_t j=0; j<Si.det_place_holder.UpMult; j++ ){
	cout << "  ChNum: " << Si.det_place_holder.UpChNum.at(j) << endl;
	cout << "  Energy Raw: " << Si.det_place_holder.EnergyUp_Raw.at(j) << endl;
	cout << "  Energy Pulser: " << Si.det_place_holder.EnergyUp_Pulser.at(j) << endl;
	cout << "  Energy Rel: " << Si.det_place_holder.EnergyUp_Rel.at(j) << endl;
	cout << "  Energy Cal: " << Si.det_place_holder.EnergyUp_Cal.at(j) << endl;
	cout << "  Time: " << Si.det_place_holder.TimeUp.at(j) << endl;
	
      }
      cout << "Down:\n";
      for (Int_t j=0; j<Si.det_place_holder.DownMult; j++ ){
	cout << "  ChNum: " << Si.det_place_holder.DownChNum.at(j) << endl;
	cout << "  Energy Raw: " << Si.det_place_holder.EnergyDown_Raw.at(j) << endl;
	cout << "  Energy Pulser: " << Si.det_place_holder.EnergyDown_Pulser.at(j) << endl;
	cout << "  Energy Rel: " << Si.det_place_holder.EnergyDown_Rel.at(j) << endl;
	cout << "  Energy Cal: " << Si.det_place_holder.EnergyDown_Cal.at(j) << endl;
	cout << "  Time: " << Si.det_place_holder.TimeDown.at(j) << endl;
      }
      cout << "BackChNum:\n";
      for (Int_t j=0; j<Si.det_place_holder.BackMult; j++ ){
	cout << "  ChNum: " << Si.det_place_holder.BackChNum.at(j) << endl;
	cout << "  Energy Raw: " << Si.det_place_holder.EnergyBack_Raw.at(j) << endl;
	cout << "  Energy Pulser: " << Si.det_place_holder.EnergyBack_Pulser.at(j) << endl;
	cout << "  Energy Rel: " << Si.det_place_holder.EnergyBack_Rel.at(j) << endl;
	cout << "  Energy Cal: " << Si.det_place_holder.EnergyBack_Cal.at(j) << endl;
	cout << "  Time: " << Si.det_place_holder.TimeBack.at(j) << endl;
      }
      cout << "---------------------------\n";
    }
    cout << "------------------------------------------------------------\n";
    cout << "Data Processed\n";
    for (Int_t i=0; i<Si.ReadHit->size(); i++){
      Si.hit_place_holder = Si.ReadHit->at(i);
      cout << "NHitsInDet: " << Si.hit_place_holder.NHitsInDet << endl;
      cout << "DetID: " << Si.hit_place_holder.DetID << endl;

      cout << "  HitType: " << Si.hit_place_holder.HitType << endl;
      cout << "  Channel Back: " << Si.hit_place_holder.BackChannel << endl;
      cout << "  Channel Front: " << Si.hit_place_holder.FrontChannel << endl;
      cout << "  Energy Back: " << Si.hit_place_holder.EnergyBack << endl;
      cout << "  Energy Front: " << Si.hit_place_holder.EnergyFront << endl;
      cout << "  Energy: " << Si.hit_place_holder.Energy << endl;
      cout << "  Time: " << Si.hit_place_holder.Time << endl;
      cout << "  X: " << Si.hit_place_holder.X << endl;
      cout << "  Z: " << Si.hit_place_holder.Z << endl;
      cout << "  XW: " << Si.hit_place_holder.XW << endl;
      cout << "  YW: " << Si.hit_place_holder.YW << endl;
      cout << "  ZW: " << Si.hit_place_holder.ZW << endl;
      cout << "  RW: " << Si.hit_place_holder.RW << endl;
      cout << "  PhiW: " << Si.hit_place_holder.PhiW << endl;
      cout << "  RFSubtract: " << Si.hit_place_holder.RFSubtract << endl;

      cout << "------------------------------------------------------------\n";

#ifdef ReadTrack
      cout << "Tracked Data\n";
      Tr.track_place_holder = Tr.ReadTrEvent->at(i);
      cout << "DetID: " << Tr.track_place_holder.DetID << endl;
      cout << "WireID: " << Tr.track_place_holder.WireID << endl;
      cout << "  SiZ: " << Tr.track_place_holder.SiZ << endl;
      cout << "  SiR: " << Tr.track_place_holder.SiR << endl;
      cout << "  SiPhi: " << Tr.track_place_holder.SiPhi << endl;
      cout << "  SiEnergy: " << Tr.track_place_holder.SiEnergy << endl;
      cout << "  PCZ: " << Tr.track_place_holder.PCZ << endl;
      cout << "  PCR: " << Tr.track_place_holder.PCR << endl;
      cout << "  PCPhi: " << Tr.track_place_holder.PCPhi << endl;
      cout << "  PCEnergy: " << Tr.track_place_holder.PCEnergy << endl;
      cout << "  IntPoint: " << Tr.track_place_holder.IntPoint << endl;
      cout << "  Theta: " << Tr.track_place_holder.Theta << endl;
      cout << "------------------------------------------------------------\n";
#endif
      
    }
    string dummy;
    cin >> dummy;
    if (dummy == "quit"){
      exit(EXIT_FAILURE);
    }

  }

  return 0;
}
