//Written by John Parker--28 Apr, 2016
//This program 'reads' the ParkerTrack.root file event by event
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

//#define ReadTrack

//Detector Libraries
#include "../include/organizetree.h"

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

  //int entry = atoi(argv[1]);

  Track Tr;
  Tr.ReadTrEvent = 0;
  Tr.ReadSiEvent = 0;
  Tr.ReadPCEvent = 0;

  TFile *inputFile = new TFile("/home2/parker/ANASEN/LSU/ParkerReconstruct_root/run410_429_050416.root");
  TTree *tree = (TTree*) inputFile->Get("MainTree");
  tree->SetBranchAddress("Tr.NTracks",&Tr.NTracks);
  tree->SetBranchAddress("Tr.NTracks1",&Tr.NTracks1);
  tree->SetBranchAddress("Tr.NTracks2",&Tr.NTracks2);
  tree->SetBranchAddress("Tr.NTracks3",&Tr.NTracks3);
  tree->SetBranchAddress("Tr.TrEvent",&Tr.ReadTrEvent);
  tree->SetBranchAddress("Tr.SiEvent",&Tr.ReadSiEvent);
  tree->SetBranchAddress("Tr.PCEvent",&Tr.ReadPCEvent);


  Long64_t nentries = tree->GetEntries();
  for (Int_t global_evt=0; global_evt<nentries; global_evt++){
    Int_t status = tree->GetEvent(global_evt);

    if ( !(Tr.NTracks1 == 2 && Tr.NTracks2==1) )
      continue;
    //if ( !(Tr.NTracks2 == 1) )
    //continue;

    cout << "Event Number: " << global_evt << endl;
    cout << "-----------------------------------------------------------------------------------\n";

    cout << "Number of Tracks: " << Tr.NTracks1 << endl;
    for (Int_t i=0; i<Tr.ReadTrEvent->size(); i++){
      Tr.track_place_holder = Tr.ReadTrEvent->at(i);
      cout << "Si Detector: " << Tr.track_place_holder.DetID << endl;
      cout << "Si Energy:   " << Tr.track_place_holder.SiEnergy << endl;
      cout << "Si Z:        " << Tr.track_place_holder.SiZ << endl;
      cout << "Si R:        " << Tr.track_place_holder.SiR << endl;
      cout << "Si Phi:      " << Tr.track_place_holder.SiPhi << endl;

      cout << "PC Detector: " << Tr.track_place_holder.WireID << endl;
      cout << "PC Energy:   " << Tr.track_place_holder.PCEnergy << endl;
      cout << "PC Z:        " << Tr.track_place_holder.PCZ << endl;
      cout << "PC R:        " << Tr.track_place_holder.PCR << endl;
      cout << "PC Phi:      " << Tr.track_place_holder.PCPhi << endl;
      cout << "--------------------\n";
    }
    cout << "---------------------------------------------\n";
    cout << "Number of Si Hits: " << Tr.NTracks2 << endl;
    for (Int_t i=0; i<Tr.ReadSiEvent->size(); i++){
      Tr.si_place_holder = Tr.ReadSiEvent->at(i);
      cout << "Si Detector: " << Tr.si_place_holder.DetID << endl;
      cout << "Si Energy:   " << Tr.si_place_holder.SiEnergy << endl;
      cout << "Si Z:        " << Tr.si_place_holder.SiZ << endl;
      cout << "Si R:        " << Tr.si_place_holder.SiR << endl;
      cout << "Si Phi:      " << Tr.si_place_holder.SiPhi << endl;
    cout << "-----------------------\n";
    }
    cout << "---------------------------------------------\n";
    cout << "Number of PC Hits: " << Tr.NTracks3 << endl;
    for (Int_t i=0; i<Tr.ReadPCEvent->size(); i++){
      Tr.pc_place_holder = Tr.ReadPCEvent->at(i);
      cout << "PC Wire:     " << Tr.pc_place_holder.WireID << endl;
      cout << "PC Energy:   " << Tr.pc_place_holder.PCEnergy << endl;
      cout << "PC Z:        " << Tr.pc_place_holder.PCZ << endl;
      cout << "PC R:        " << Tr.pc_place_holder.PCR << endl;
      cout << "PC Phi:      " << Tr.pc_place_holder.PCPhi << endl;
      cout << "----------------------\n";
    }

    string dummy;
    cin >> dummy;
    if (dummy == "quit"){
      exit(EXIT_FAILURE);
    }

  }

  return 0;
}
