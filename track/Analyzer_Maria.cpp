//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Goal: To Analyze the Tracking and Reconstruction of the 7Be+d-> p + alpha + alpha Events  
// & to Analyze the Elastic Scattering of Deuterons to Calibrate ANASEN for 7Be+d Experiments..
// & for the other (d,p),(d,alpha)..etc..ANASEN experiments with Gas volume target
//
// //To create a dictionary:
//  rootcint -f tr_dict.cxx -c tree_structure.h LinkDef.h
//
// Usage: g++ -o Analyzer_Maria tr_dict.cxx LookUp.cpp Analyzer_Maria.cpp `root-config --cflags --glibs` -O3
//
// ./Analyzer_Maria  DataList.txt /data0/manasta/OrganizeRaw_files/run###.root  cuts/alphasRun924newPCThres03212017QQQ.root //
//
// Uses Lookup tables instead of doing integration multiple times for Energyloss, 
// Final Energy, Initial Energy & Distance calculation.
//
// Author: Nabin Rijal, 2016 September.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define FillTree
#define FillEdE_cor
#define CheckBasic
#define IsCal

#define DiffIP 2 //cm
#define ConvAngle 180./TMath::Pi() //when multiplied, Converts to Degree from Radian 

#define EdE

//#define PCWireCal
//#define PCPlots
#define MCP_RF_Cut

#define MaxSiHits   500
#define MaxADCHits  500
#define MaxTDCHits  500
#define MaxTracks   100

#define pcr 3.846284509 //3.75+0.096284509; //correction for the centroid Kx applied
#define La 55.0545   //Length of ANASEN gas volume as measured 2/22/2017 with Lagy

/#define gold_pos 27.7495 //cm based on geometry measurements we did with Lagy at 2/22/2017 run930
//#define gold_pos 22.8981   // spacer 1 = all in - 4.8514 cm run932
//#define gold_pos 16.8783 // spacer 2 = all in - 10.8712 cm run934
//#define gold_pos 15.6083 // spacer 3 = all in - 12.1412 cm run936
//#define gold_pos 12.4587 // spacer 4 = all in - 15.2908 cm
#define gold_pos -2.8505 // spacer 7 = all in - 30.6 cm

///////////////////Nuclear Masses ///////////////////////////////////////////////////
//---MARIA nuclear masses //MeV--------------------------------------
#define M_P 938.27197
#define M_alpha 3727.37892
#define M_16O 14895.079
#define M_17F 17692.29961
#define M_18Ne 16767.09917
#define M_21Na 19553.56884
#define M_24Mg 22335.79043
#define M_27Al 25126.49834
//---maria-----------------------------------

//#define BeamE 54.24 // 16O beam for calibration in 18Ne run
//#define BeamE 75.7 // 24Mg beam run
#define BeamE 71.62 // 18Ne beam !!! CORRECTED 5/23/2017 76.55
#define QValue 2.63819 // MeV 18Ne(a,p)21Na 

/////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <map>

//ROOT libraries
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TObjArray.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TRandom1.h>
#include <TList.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "tree_structure.h"
#include "LookUp.h"
#include "/home/manasta/Desktop/parker_codes/Include/ReconstructMaria.h" // so that the Reconstruction process is in separate script
//#include "/home/maria/rayMountPoint/Desktop/parker_codes/Include/ReconstructMaria.h"
//#include "/home/manasta/Desktop/parker_codes/Include/EnergyLoss.h" // used to be the method to use

using namespace std;
////////////////////////////////////////////////////////////////////////////////////
Int_t FindMaxPC(Double_t phi, PCHit& PC);

void MyFill(string name,int binsX, double lowX, double highX, double valueX);

void MyFill(string name,int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY);

TList* fhlist;
std::map<string,TH1*> fhmap;
////////////////////////////////////////////////////////////////////////////////////
bool Track::Tr_Sisort_method(struct TrackEvent a,struct TrackEvent b){
  if(a.SiEnergy > b.SiEnergy)
    return 1;
  return 0;
};
////////////////////////////////////////////////////////////////////////////////////
bool Track::Tr_PCsort_method(struct TrackEvent c,struct TrackEvent d){
  if(c.PCEnergy > d.PCEnergy)
    return 1;
  return 0;
};
////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) { 

  //Don't know what this does, but libraries won't load without it
  TApplication *myapp=new TApplication("myapp",0,0); 

  Int_t numarg=4;
#ifdef IsCal
  numarg=3;
#endif
  if (argc!=numarg) {
    cout << "Error: Wrong Number of Arguments\n";
    exit(EXIT_FAILURE);
  }

  char* file_raw  = new char [300]; // for input .root file
  char* file_cal = new char [300]; // for output .root file

  strcpy( file_raw, argv[1] );
  strcpy( file_cal, argv[2] );

  cout << argv[0] << endl;
  cout << argv[1] << endl;
  cout << argv[2] << endl;
  
#ifdef IsCal
  //////////////// CUTS ////////////////////

  char* file_cut1 = new char[300]; //for allcut
  strcpy( file_cut1, argv[3] );
 
  //He4 cut
  TFile *cut_file1 = new TFile(file_cut1);
  if (!cut_file1->IsOpen()){
    cout << "Cut file1: " << file_cut1 << " could not be opened.\n";
    exit(EXIT_FAILURE);
  }
  TCutG *cut1 = NULL;
  //cut1 = (TCutG*)cut_file1->Get("alphasRun924newPCThres03212017QQQ"); // cut for 160(a,a) run924
  cut1 = (TCutG*)cut_file1->Get("alphas_run778_911_timecut_protoncut_Qvalue_BeamCorrect_05252017"); // cut for 18Ne(a,a) run778_911_timecut
  
  if (cut1 == NULL){
    cout << "Cut1 does not exist\n";
    exit(EXIT_FAILURE);
  }
  
  ///////----------cut2-------------------///////////////////
   
  char* file_cut2 = new char[300]; //for allcut
  strcpy( file_cut2, argv[3] );
    
  //proton cut
  TFile *cut_file2 = new TFile(file_cut2);
  if (!cut_file2->IsOpen()){
    cout << "Cut file2: " << file_cut2 << " could not be opened.\n";
    exit(EXIT_FAILURE);
  }
  TCutG *cut2 = NULL;
  cut2 = (TCutG*)cut_file2->Get("protonsQ3R1_run778_911_05232017");
  
  if (cut2 == NULL){
    cout << "Cut2 does not exist\n";
    exit(EXIT_FAILURE);
  }
cout << argv[3] << endl;
#endif
 
////////////////////////////////////////////////////////////////////////////////////////////////////
  TObjArray *RootObjects = new TObjArray();
  SiHit Si;
  PCHit PC;
  ///CsIHit CsI;
  Track Tr; 
  Int_t RFTime, MCPTime;
  Int_t Old_RFTime,Old_MCPTime;

  Double_t PCGoodEnergy[24];
  
  Si.ReadDet = 0;
  Si.ReadHit = 0;
  PC.ReadHit = 0;
  //CsI.ReadHit = 0;

  //--------------------------------MARIA Eloss---------------------------------------------------
 
  ///////-----------------E_Loss_16O-------------/////////////////////////////////////////////////////////////////////


  //LookUp *E_Loss_16O = new LookUp("/home/maria/rayMountPoint/Desktop/anasen_analysis_software/srim_files/16O_in_HeCO2_377Torr_18Nerun.eloss",M_16O); // when I work from home
  //LookUp *E_Loss_alpha = new LookUp("/home/maria/rayMountPoint/Desktop/anasen_analysis_software/srim_files/He_in_HeCO2_377Torr_18Nerun.eloss",M_alpha); 


  /*  
  LookUp *E_Loss_16O = new LookUp("/home/manasta/Desktop/anasen_analysis_software/srim_files/16O_in_HeCO2_377Torr_18Nerun.eloss",M_16O);  
  E_Loss_16O->InitializeLookupTables(80.0,1200.0,0.02,0.04);
  //EnergyLoss *E_Loss_16O = new EnergyLoss("/home/manasta/Desktop/anasen_analysis_software/srim_files/16O_in_HeCO2_377Torr_18Nerun.eloss",M_16O);

  
  LookUp *E_Loss_alpha = new LookUp("/home/manasta/Desktop/anasen_analysis_software/srim_files/He_in_HeCO2_377Torr_18Nerun.eloss",M_alpha);
  E_Loss_alpha->InitializeLookupTables(50.0,1800.0,0.02,0.04); 
  //EnergyLoss *E_Loss_alpha = new EnergyLoss("/home/manasta/Desktop/anasen_analysis_software/srim_files/He_in_HeCO2_377Torr_18Nerun.eloss",M_alpha); 
  

  //LookUp *E_Loss_proton = new LookUp("/home/manasta/Desktop/anasen_analysis_software/srim_files/H_in_HeCO2_377Torr_18Nerun.eloss",M_P); 
  //E_Loss_proton->InitializeLookupTables(30.0,9000.0,0.02,0.04); 
  
  */

   ///////-----------------E_Loss_18Ne-------------/////////////////////////////////////////////////////////////////////


 
  //LookUp *E_Loss_18Ne = new LookUp("/home/maria/rayMountPoint/Desktop/anasen_analysis_software/srim_files/18Ne_in_HeCO2_377Torr_18Nerun.eloss",M_18Ne);  
  LookUp *E_Loss_18Ne = new LookUp("/home/manasta/Desktop/anasen_analysis_software/srim_files/18Ne_in_HeCO2_377Torr_18Nerun.eloss",M_18Ne);  
  E_Loss_18Ne->InitializeLookupTables(90.0,850.0,0.02,0.04);
  //E_Loss_18Ne->PrintLookupTables();
  //EnergyLoss *E_Loss_18Ne = new EnergyLoss("/home/manasta/Desktop/anasen_analysis_software/srim_files/16O_in_HeCO2_377Torr_18Nerun.eloss",M_18Ne);

  //LookUp *E_Loss_alpha = new LookUp("/home/maria/rayMountPoint/Desktop/anasen_analysis_software/srim_files/He_in_HeCO2_377Torr_18Nerun.eloss",M_alpha);
  LookUp *E_Loss_alpha = new LookUp("/home/manasta/Desktop/anasen_analysis_software/srim_files/He_in_HeCO2_377Torr_18Nerun.eloss",M_alpha);
  E_Loss_alpha->InitializeLookupTables(50.0,1800.0,0.02,0.04); 
  //EnergyLoss *E_Loss_alpha = new EnergyLoss("/home/manasta/Desktop/anasen_analysis_software/srim_files/He_in_HeCO2_377Torr_18Nerun.eloss",M_alpha); 
  
  //LookUp *E_Loss_proton = new LookUp("/home/maria/rayMountPoint/Desktop/anasen_analysis_software/srim_files/H_in_HeCO2_377Torr_18Nerun.eloss",M_P); 
  LookUp *E_Loss_proton = new LookUp("/home/manasta/Desktop/anasen_analysis_software/srim_files/H_in_HeCO2_377Torr_18Nerun.eloss",M_P); 
  E_Loss_proton->InitializeLookupTables(30.0,9000.0,0.02,0.04); 
  
  
  ///////-----------------E_Loss_24Mg-------------/////////////////////////////////////////////////////////////////////
  
  /*
  LookUp *E_Loss_24Mg = new LookUp("/home/manasta/Desktop/anasen_analysis_software/srim_files/24Mg_in_HeCO2_303Torr_24Mgrun.eloss",M_24Mg);
  E_Loss_24Mg->InitializeLookupTables(80.0,70.0,0.02,0.04);
  //getchar();
  //EnergyLoss *E_Loss_16O = new EnergyLoss("/home/manasta/Desktop/anasen_analysis_software/srim_files/16O_in_HeCO2_377Torr_18Nerun.eloss",M_16O);
  LookUp *E_Loss_alpha = new LookUp("/home/manasta/Desktop/anasen_analysis_software/srim_files/He_in_HeCO2_303Torr_24Mgrun.eloss",M_alpha); 
  E_Loss_alpha->InitializeLookupTables(30.0,800.0,0.02,0.04); 
  */
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  TFile *outputfile = new TFile(file_cal,"RECREATE");
  TTree *MainTree = new TTree("MainTree","MainTree");
  //++
  MainTree->Branch("Tr.NTracks",&Tr.NTracks,"NTracks/I");
  MainTree->Branch("Tr.NTracks1",&Tr.NTracks1,"NTracks1/I");
  MainTree->Branch("Tr.NTracks2",&Tr.NTracks2,"NTracks2/I");
  MainTree->Branch("Tr.NTracks3",&Tr.NTracks3,"NTracks3/I");
  MainTree->Branch("Tr.TrEvent",&Tr.TrEvent); 
  MainTree->Branch("Tr.HeavyEvent",&Tr.HeavyEvent);
  MainTree->Branch("RFTime",&RFTime,"RFTime/I");
  MainTree->Branch("MCPTime",&MCPTime,"MCPTime/I");

  RootObjects->Add(MainTree);

  fhlist = new TList;
  RootObjects->Add(fhlist);  
  ////////////////////////////////////////////////////  
  Int_t count_2A_1P=0, count_2A=0, count_1A_1P=0;
  Int_t count_1track=0, count_2tracks=0, count_3tracks =0,count_morethan3tracks=0;
  Int_t count_3He=0;
  Int_t count_3He_NT2;
  ////////////////////////////////////////////////////
  //============  Read the root files from the Data list ==============  
  ifstream inFileList;
  inFileList.open(file_raw);

  if (!inFileList.is_open()){
    cout << "In File List not open\n";
    exit(EXIT_FAILURE);
  }
  
  string rootfile;
  char rootfile_char[200];

  //------------------counters--------------------------//
  /////////////////////////////////////////////////////////
  Int_t counter0 = 0;
  Int_t counter = 0;
  Int_t counter1 = 0;
  Int_t counter2 = 0;
  Int_t counter3 = 0;
  Int_t counter4 = 0;
  Int_t counter5 = 0;
  Int_t counter6 = 0;
  Int_t counter7 = 0;
  Int_t counter8 = 0;
  Int_t counter9 = 0;
  Int_t counter10 = 0;
  Int_t counter11 = 0;
  Int_t counter12 = 0;
  Int_t counter13 = 0;
  Int_t counter14 = 0;
  Int_t counterNeg = 0;
  Int_t counterPos = 0;
  Int_t counterPCMINUSTEN = 0;
  Int_t counterZeroTracks = 0;
  Int_t counterZeroTracks1 = 0;
  Int_t counterZeroTracks2 = 0;
  Int_t counterZeroTracks3 = 0;
  Int_t counterEn1000 = 0;
  Int_t counterEnLessZero = 0;
  Int_t counterEnPC1000 = 0;
  Int_t counterEnPCLessZero = 0;
  
  Int_t counterBeamE_largeQQQ = 0;
  Int_t counterBeamE_largeR1 = 0;
  Int_t counterBeamE_largeR2 = 0;
  Int_t counterBeamE_largeR2_special = 0;
  Int_t counterBeamE_large_HeR2_special = 0; 

  /////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  /*
  ofstream outfileQQQ;
  ofstream outfileR1;
  ofstream outfileR2;
  ofstream outfileR2special;
  ofstream outfile_HeR2special;
  outfileQQQ.open("LargeEnergiesQQQ_run930_newthres.txt");
  outfileR1.open("LargeEnergiesR1_run930_newthres.txt");
  outfileR2.open("LargeEnergiesR2_run930_newthres.txt");
  outfileR2special.open("specialR2_run930_newthres.txt");
  outfile_HeR2special.open("specialHeR2_run930_newthres.txt");
  */  

  while (!inFileList.eof()){//===================loop over all of the incoming root files================
    getline(inFileList,rootfile);
    if (rootfile.empty()){
      //cout << "Root File does not exist: " << rootfile << endl;
      break;
    }
    strcpy(rootfile_char,rootfile.c_str());
    TFile *inputFile = new TFile(rootfile_char);//open root file and make sure it exists-----
    if (!inputFile->IsOpen()){
      cout << "Root file: " << rootfile << " could not be opened.\n";     
      continue;
    }   
    cout << "Processing File: " << rootfile << endl;   

    TTree *raw_tree = (TTree*) inputFile->Get("MainTree");
    raw_tree->SetBranchAddress("Si.NSiHits",&Si.NSiHits);
    raw_tree->SetBranchAddress("Si.Detector",&Si.ReadDet);
    raw_tree->SetBranchAddress("Si.Hit",&Si.ReadHit);
    raw_tree->SetBranchAddress("PC.NPCHits",&PC.NPCHits);
    raw_tree->SetBranchAddress("PC.Hit",&PC.ReadHit);
    //raw_tree->SetBranchAddress("RFTime",&RFTime);
    //raw_tree->SetBranchAddress("MCPTime",&MCPTime);
    raw_tree->SetBranchAddress("RFTime",&Old_RFTime);
    raw_tree->SetBranchAddress("MCPTime",&Old_MCPTime);
 
    Long64_t nentries = raw_tree->GetEntries();
    cout<<"nentries = "<<nentries<<endl;

    Int_t status;
    for (Long64_t i=0; i<nentries; i++){//====================loop over all events=================
      //cout<<" i =  "<<i<<endl;

      status = raw_tree->GetEvent(i);
      std::cout << "\rDone: " << i*100./nentries << "%          " << std::flush;
    
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      ////------------------------- MCP/ RF Time --------------------///////////////////////////////////

      Double_t correct=1.004009623;
      MCPTime = Old_MCPTime;
      RFTime = Old_RFTime;

      MyFill("Timing",600,1,600,fmod((MCPTime*correct-RFTime),546));
      
      //if((fmod((MCPTime*correct-RFTime),546)>26.54 && fmod((MCPTime*correct-RFTime),546)<53.86) || (fmod((MCPTime*correct-RFTime),546)>302.91 && fmod((MCPTime*correct-RFTime),546)<326.21))
      //	MyFill("Timing_Cut",600,1,600,fmod((MCPTime*correct-RFTime),546));

#ifdef MCP_RF_Cut    
      if(MCPTime > 0 && RFTime>0){
	//cout<<"   RFTime  == "<<RFTime<<"   MCPTime  =="<<MCPTime<<endl;

	if( (fmod((MCPTime*correct - RFTime),546)<26.54) || (fmod((MCPTime*correct - RFTime),546)>53.86  && fmod((MCPTime*correct - RFTime),546)<302.91) || fmod((MCPTime*correct - RFTime),546)>326.21 ){
	  //if( (((MCPTime - RFTime)% 538)<60) || (((MCPTime - RFTime)% 538)>110  && ((MCPTime - RFTime)% 538)<325) || ((MCPTime - RFTime)% 538)>380 ){
	  continue;
	}
      }else{
	continue;
      }
#endif
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      Tr.zeroTrack();
      Int_t GoodPC = -1;  
      for(Int_t k=0; k<24;k++) {
	  PCGoodEnergy[k]=0;
	}
      /////////////////////////////////////////////////////////////////////////////////////////////////////    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
      //cout<<"Si.ReadHit->size() = "<<Si.ReadHit->size()<<endl;  
      MyFill("Si_ReadHit_size",500,0,50,Si.ReadHit->size());  

      for (Int_t j=0; j<Si.ReadHit->size(); j++) {//loop over all silicon
	
	Si.hit_obj = Si.ReadHit->at(j);	//if we have a good hit type set the parameters in your new tree

    
	 if ( Si.hit_obj.Energy <= 0 ) {
	      counterNeg++;
	      continue;
	}
	 //else{
	   
	  counterPos++;

	  GoodPC = FindMaxPC(Si.hit_obj.PhiW, PC);

	  if (GoodPC > -1){//if a PC is found do Tracking

	    PC.pc_obj = PC.ReadHit->at(GoodPC);

	    Tr.ZeroTr_obj();

	    Tr.track_obj.TrackType = 1;
	  
	    Tr.track_obj.SiEnergy = Si.hit_obj.Energy;
	    Tr.track_obj.SiPhi = Si.hit_obj.PhiW;	  
	    Tr.track_obj.SiZ = Si.hit_obj.ZW;
	    Tr.track_obj.SiR = Si.hit_obj.RW;
	    Tr.track_obj.DetID = Si.hit_obj.DetID;
	    Tr.track_obj.SiBCh = Si.hit_obj.BackChannel;
	    Tr.track_obj.HitType = Si.hit_obj.HitType;

	    /*
	    if (  Tr.track_obj.SiEnergy <= 0 ){
	      counter++;
	      continue;
	    }
	    else {counter0++;}
	    */
	  
	    Tr.track_obj.PCEnergy = PC.pc_obj.Energy;

	    // check to see if I really cut the low PCEnergies from the Main 4/11/2017

	    if (Tr.track_obj.PCEnergy < 0.003 && Tr.track_obj.PCEnergy>0 ){
	      cout << "PCEnergy " << Tr.track_obj.PCEnergy << endl;
	      exit(1);
	    }

	    /////////////////////////////////////////////

	    Tr.track_obj.PCPhi = PC.pc_obj.PhiW;	  
	    Tr.track_obj.PCZ = PC.pc_obj.ZW; ////////////sos!!!!!!!!!!!!//////////////////////
	    Tr.track_obj.PCR = PC.pc_obj.RW;
	    Tr.track_obj.WireID = PC.pc_obj.WireID;

	    //had a lot of zeros in my PCZ coming from these wires that don't work so need to exclude them

	    //if (Tr.track_obj.WireID == 0 || Tr.track_obj.WireID == 6 || Tr.track_obj.WireID == 16 || Tr.track_obj.WireID == 17) //24Mg data
	    //  Tr.track_obj.PCZ = -10.0;  

	    if ( Tr.track_obj.WireID == 6 || Tr.track_obj.WireID == 16){  // added in 4/11/2017
	        Tr.track_obj.PCZ = -10.0;
		Tr.track_obj.PCEnergy = -10.0;
		continue;// added in 05/01/2017 cause it was still counting the PCEnergies for these wires showing up in PCPhi plots where it should have been empty
	    }
	     
	    // check if I have actual PCZ=0 not just very small, close to zero  // 4/11/2017
	    if(Tr.track_obj.PCZ==0)
	      cout << " i " << i << " " << j << " " << Tr.track_obj.PCZ << endl;


	    PCGoodEnergy[Tr.track_obj.WireID] = Tr.track_obj.PCEnergy;
	    

	    // eliminate Si and Wire from further tracking
	    
	    Si.ReadHit->at(j).Energy = -1000;  // marks the current tracked hit with Energy = -1000 so that it is not counted again
	    //counterPCMINUSTEN++;
	    PC.ReadHit->at(GoodPC).Energy = -10;

	    Tr.NTracks++; //total no. of tracks...Good or bad
	    Tr.NTracks1++;//good tracks...PC & Si both
	  }else{
	    continue;
	  }
	  Tr.TrEvent.push_back(Tr.track_obj);
	  //}  // closes "else" after SiEnergies<0
      }  
      sort( Tr.TrEvent.begin(), Tr.TrEvent.begin()+Tr.NTracks1,Tr.Tr_Sisort_method );    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      for (Int_t k=0; k<Si.ReadHit->size(); k++){//loop over all silicon
	
	Si.hit_obj = Si.ReadHit->at(k);
	/*
	if ( (Si.hit_obj.Energy == -1000) || (Si.hit_obj.Energy <= 0)){ //make sure that the Silicon energy was filled
	  counter8++;
	  continue;
	}else{
	  counter10++;
	}
	*/
	if (Si.hit_obj.Energy == -1000) {
	    counterEn1000++;
	    continue;
	  }
	else if(Si.hit_obj.Energy <= 0) {
	    counterEnLessZero++;
	    continue;
	  }
	else {
	  counter10++;
	
	  Tr.ZeroTr_obj();
	  Tr.track_obj.TrackType = 2;	
	  Tr.track_obj.SiEnergy = Si.hit_obj.Energy;
	  Tr.track_obj.SiPhi = Si.hit_obj.PhiW;
	  Tr.track_obj.SiZ = Si.hit_obj.ZW;	  
	  Tr.track_obj.SiR = Si.hit_obj.RW;	  
	  Tr.track_obj.DetID = Si.hit_obj.DetID;
	  Tr.track_obj.SiBCh = Si.hit_obj.BackChannel;
	  Tr.track_obj.HitType = Si.hit_obj.HitType;

	  Tr.NTracks2++;//Only Si && no PC
	  Tr.NTracks++;//total no. of tracks...Good or bad
	  Tr.TrEvent.push_back(Tr.track_obj);
	}
      }      
      sort( Tr.TrEvent.begin()+Tr.NTracks1, Tr.TrEvent.begin()+Tr.NTracks1+Tr.NTracks2,Tr.Tr_Sisort_method );  
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //cout<<"PC.ReadHit->size() = "<<PC.ReadHit->size()<<endl;
      MyFill("PC_ReadHit_size",500,0,50,PC.ReadHit->size());  

      for (Int_t l=0; l<PC.ReadHit->size(); l++ ){//loop over pc
	
	PC.pc_obj = PC.ReadHit->at(l);
	/*
	if((PC.ReadHit->at(l).Energy == -10) || (PC.ReadHit->at(l).Energy <= 0)){
	  counter14++;
	  continue;
	}else{
	  counter12++;
	}
	*/

	if (PC.ReadHit->at(l).Energy == -10)
	  {
	    counterEnPC1000++;
	    continue;
	  }
	else if(PC.ReadHit->at(l).Energy <= 0)
	  {
	    counterEnPCLessZero++;
	    continue;
	  }
	else{
	  counter12++;
	
	  Tr.ZeroTr_obj();  
	  Tr.track_obj.TrackType = 3;

	  Tr.track_obj.PCEnergy = PC.pc_obj.Energy;

	   if (Tr.track_obj.PCEnergy < 0.003 && Tr.track_obj.PCEnergy>0 ){
	      cout << "PCEnergy " << Tr.track_obj.PCEnergy << endl;
	      exit(1);}

	  Tr.track_obj.WireID = PC.pc_obj.WireID;	
	  Tr.track_obj.PCZ = PC.pc_obj.ZW;
	  Tr.track_obj.PCR = PC.pc_obj.RW;
	  Tr.track_obj.PCPhi = PC.pc_obj.PhiW;

	   //had a lot of zeros in my PCZ coming from these wires that don't work so need to exclude them
	  
	  // if (Tr.track_obj.WireID == 0 || Tr.track_obj.WireID == 6 || Tr.track_obj.WireID == 16 || Tr.track_obj.WireID == 17)  //24Mg data
	  //    Tr.track_obj.PCZ = -10.0;

	  if ( Tr.track_obj.WireID == 6 || Tr.track_obj.WireID == 16){
	        Tr.track_obj.PCZ = -10.0;
		Tr.track_obj.PCEnergy = -10.0;
		continue;
	    }

	 
	  Tr.NTracks3++;//Only PC && no Si
	  Tr.NTracks++;//total no. of tracks...Good && bad both
	  
	  PCGoodEnergy[Tr.track_obj.WireID] = Tr.track_obj.PCEnergy;

	  Tr.TrEvent.push_back(Tr.track_obj);
	}
      }        
      sort(Tr.TrEvent.begin()+Tr.NTracks1+Tr.NTracks2, Tr.TrEvent.end(),Tr.Tr_PCsort_method );      
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //////////////////////////////////////////////////////////////////////////////////////
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //cout<<"Tr.NTracks = "<<Tr.NTracks<<" Tr.NTracks1 = "<<Tr.NTracks1<<" Tr.NTracks2 = "<<Tr.NTracks2<<" Tr.NTracks3 = "<<Tr.NTracks3<<endl;      
    
      //check for neighboring wires cross-talk

      if(Tr.NTracks1==1){
	if(Tr.TrEvent[0].WireID >-1 && Tr.TrEvent[0].WireID < 23)
	  MyFill(Form("PCEnergyPlus_Wire%i_Wire%i",Tr.TrEvent[0].WireID+1,Tr.TrEvent[0].WireID),300,0,0.25,PCGoodEnergy[Tr.TrEvent[0].WireID],300,0,0.25,PCGoodEnergy[Tr.TrEvent[0].WireID+1]);
	
	if(Tr.TrEvent[0].WireID >0 && Tr.TrEvent[0].WireID < 24)
	  MyFill(Form("PCEnergyMinus_Wire%i_Wire%i",Tr.TrEvent[0].WireID-1,Tr.TrEvent[0].WireID),300,0,0.25,PCGoodEnergy[Tr.TrEvent[0].WireID],300,0,0.25,PCGoodEnergy[Tr.TrEvent[0].WireID-1]);

      }
	 

      //////////////test//////////////////////////////

      if(Tr.NTracks1==1 && Tr.TrEvent[0].DetID > 15 && Tr.TrEvent[0].DetID < 28){
	for(Int_t m=0; m<Tr.NTracks3; m++)
	  {
	    MyFill("PCEnergy_NTracks3_R2",1000,-0.1,1,Tr.TrEvent[m+Tr.NTracks1+Tr.NTracks2].PCEnergy);
	  }
      }else if(Tr.NTracks1==1 && Tr.TrEvent[0].DetID>-1 && Tr.TrEvent[0].DetID < 16){
	for(Int_t m=0; m<Tr.NTracks3; m++)
	  {
	    MyFill("PCEnergy_NTracks3_QQQ_R1",1000,-0.1,1,Tr.TrEvent[m+Tr.NTracks1+Tr.NTracks2].PCEnergy);
	  }
      }


      ////////////// checking for the heavy hit &/or cross talk in the wire ///////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////    
#ifdef PCPlots
      Int_t pct;
      for(Int_t pc=0; pc<Tr.NTracks; pc++) {
	//cout<<"   Check 10 " <<Tr.NTracks<<endl;

	//All PCWire vs Energy
	MyFill("WireID_vs_PCEnegy",25,0,24,Tr.TrEvent[pc].WireID,500,0,0.35,Tr.TrEvent[pc].PCEnergy);

	for(Int_t pca=0; pca<Tr.NTracks1; pca++){

	  //PCWire with track in channel 12 & other tracks & non-tracks in other channels
	  MyFill("WireID_mod1_vs_PCEnegy_12",
		 25,0,24,(((Int_t)Tr.TrEvent[pc].WireID-(Int_t)Tr.TrEvent[pca].WireID +12)%24),
		 500,0,0.35,Tr.TrEvent[pc].PCEnergy);
	  MyFill("WireID_mod1_vs_PCEnegy_0" ,
		 25,0,24,(((Int_t)Tr.TrEvent[pc].WireID-(Int_t)Tr.TrEvent[pca].WireID)%24),
		 500,0,0.35,Tr.TrEvent[pc].PCEnergy);
	}
      }

      for(Int_t m=0; m<Tr.NTracks; m++) {
	MyFill(Form("PCWireID_DetID_%i",Tr.TrEvent[m].DetID),25,0,24,Tr.TrEvent[m].WireID); 
	MyFill(Form("PCWireID_DetID2_%i",Tr.TrEvent[0].DetID),25,0,24,Tr.TrEvent[m].WireID); 
      }
#endif     

      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ///////////////////////////////////For the Tracking///////////////////////////////////    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
      //reconstruction variables
      Double_t m = 0, b = 0; 

      for(Int_t p=0; p<Tr.NTracks1;p++){
	
	//CCCCCCCCCCCCCCCCCCCCCCCCCC///Let's put some checks //2016July28 CCCCCCCCCCC	
#ifdef CheckBasic
	if(Tr.TrEvent[p].SiZ < -4.0 || Tr.TrEvent[p].PCZ < -4.0){
	  continue;
	}
	if(Tr.TrEvent[p].PCR>4.0 || Tr.TrEvent[p].PCR<3.5){
	  continue;
	}
	if(Tr.TrEvent[p].SiR>11.0 || Tr.TrEvent[p].SiR<4.0){
	  continue;
	}
#endif
	//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	  
	//if(Tr.TrEvent[p].SiZ > 0 && Tr.TrEvent[p].PCZ <= 0){
	////cout<<" SiZ ="<<Tr.TrEvent[p].SiZ<<" PCZ ="<<Tr.TrEvent[p].PCZ<<endl;
	//}
	//cout<<"SiR = "<<Tr.TrEvent[p].SiR<<" PCR ="<<Tr.TrEvent[p].PCR<<" SiZ ="<<Tr.TrEvent[p].SiZ<<" PCZ ="<<Tr.TrEvent[p].PCZ<<endl;

	m = (Tr.TrEvent[p].PCR-Tr.TrEvent[p].SiR)/(Tr.TrEvent[p].PCZ-Tr.TrEvent[p].SiZ);
	b = (Tr.TrEvent[p].PCR - m*Tr.TrEvent[p].PCZ);
	Tr.TrEvent[p].IntPoint = -b/m;
	////cout<<" IntPoint = "<<Tr.TrEvent[p].IntPoint<<endl;
	//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#ifdef CheckBasic
	if(Tr.TrEvent[p].IntPoint<-5.0 || Tr.TrEvent[p].IntPoint>55.0){
	  continue;
	}
#endif
	//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	/////////////////////Calculating Theta and PathLength/////////////////////////////////////////////
	// Theta is the angle between particle and beam, but our beam points in the negative z-direction 
	if ((Tr.TrEvent[p].IntPoint - Tr.TrEvent[p].SiZ) > 0){
	  // This is forward-scattering, so we have theta b/w 0 and 90 
	  Tr.TrEvent[p].Theta = atan(Tr.TrEvent[p].SiR/(Tr.TrEvent[p].IntPoint - Tr.TrEvent[p].SiZ));
	  Tr.TrEvent[p].PathLength = Tr.TrEvent[p].SiR/sin(Tr.TrEvent[p].Theta);
	 
	  //cout<<" Tr.TrEvent[p].Theta1 =  "<<Tr.TrEvent[p].Theta*ConvAngle<<" Tr.TrEvent[p].PathLength1 = "<<Tr.TrEvent[p].PathLength<<endl;
	}
	else if ((Tr.TrEvent[p].IntPoint - Tr.TrEvent[p].SiZ) < 0){
	  Tr.TrEvent[p].Theta = TMath::Pi() + atan(Tr.TrEvent[p].SiR/(Tr.TrEvent[p].IntPoint - Tr.TrEvent[p].SiZ));
	  Tr.TrEvent[p].PathLength = Tr.TrEvent[p].SiR/sin(Tr.TrEvent[p].Theta);
	
	  if( Tr.TrEvent[p].Theta <(TMath::Pi()/2)  ||  Tr.TrEvent[p].Theta > TMath::Pi())
	    cout << "Theta: " << Tr.TrEvent[p].Theta*ConvAngle << endl;

	  //if(Tr.TrEvent[p].Theta>0){
	  //cout<<" Tr.TrEvent[p].Theta2 =  "<<Tr.TrEvent[p].Theta*ConvAngle<<" Tr.TrEvent[p].PathLength2 = "<<Tr.TrEvent[p].PathLength<<endl;
	  //}
	}
	else{
	  Tr.TrEvent[p].Theta = TMath::Pi()/2;
	  Tr.TrEvent[p].PathLength = Tr.TrEvent[p].SiR;

	  //cout<<" Tr.TrEvent[p].Theta3 =  "<<Tr.TrEvent[p].Theta*ConvAngle<<" Tr.TrEvent[p].PathLength3 = "<<Tr.TrEvent[p].PathLength<<endl;
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////
	//------------------for the BeamEnergy & EnergyLoss----------/////////
	
	if(Tr.TrEvent[p].IntPoint>0.0 && Tr.TrEvent[p].IntPoint<La){

	  Tr.TrEvent[p].EnergyLoss = E_Loss_18Ne->GetEnergyLoss(BeamE,(La-Tr.TrEvent[p].IntPoint));
	  //Tr.TrEvent[p].BeamEnergy = BeamE - Tr.TrEvent[p].EnergyLoss;
	    
	  if((La-Tr.TrEvent[p].IntPoint)>0.0 && (La-Tr.TrEvent[p].IntPoint)<La){

	    Tr.TrEvent[p].BeamEnergy = E_Loss_18Ne->GetLookupEnergy(BeamE,(La-Tr.TrEvent[p].IntPoint)); // with LookUp
	    //Tr.TrEvent[p].BeamEnergy = E_Loss_16O->GetFinalEnergy(BeamE,La-Tr.TrEvent[p].IntPoint,0.1); // used to do with EnergyLoss.cpp

	    /*
	    Double_t E_4He_rxn = 0;
	    E_4He_rxn = E_Loss_alpha->GetLookupEnergy(Tr.TrEvent[p].SiEnergy,(-Tr.TrEvent[p].PathLength));

	      cout <<  Tr.TrEvent[p].DetID  << " " << Tr.TrEvent[p].Theta*180/3.14 << " " << Tr.TrEvent[p].SiEnergy <<
		" " << E_4He_rxn << " " << Tr.TrEvent[p].PathLength << " " << Tr.TrEvent[p].IntPoint << " " << Tr.TrEvent[p].BeamEnergy <<  
	       " " << Tr.TrEvent[p].PCZ << " " << Tr.TrEvent[p].SiZ << " " << Tr.TrEvent[p].SiR << endl;

	    */
	   	    
	    //cout << " La-IntPoint: " << (La-Tr.TrEvent[p].IntPoint) << " BeamEnergy: " << Tr.TrEvent[p].BeamEnergy <<endl;
	  }
	  // else{
	  //cout << " La-IntPoint: " << (La-Tr.TrEvent[p].IntPoint) << endl;
	  //}
	}

	//else{
	  // cout << "\n La-IntPoint: " << (La-Tr.TrEvent[p].IntPoint) << " BeamEnergy: " << Tr.TrEvent[p].BeamEnergy <<endl;
	  // cout << " IntPoint: " << Tr.TrEvent[p].IntPoint << " slope: " << m << 
	  //" PCZ: " << Tr.TrEvent[p].PCZ << " SiZ: " << Tr.TrEvent[p].SiZ << endl;
	  
	  //continue;
	  //}
	      
	////////////////////////////////////////////////////////////////////////////////////////
      }//end of for loop Tracking.
      ////////////////////////////////////////////////////////////////////////////////////////
      

#ifdef PCWireCal      
      Double_t mpc,bpc;
      Double_t mpc1,bpc1;

      for(Int_t s=0; s<Tr.NTracks1;s++){

	if(Tr.NTracks!=Tr.NTracks1)   // added 05/01/2017 since it eliminates some events and makes it look better!
	 continue;

	//if (cut2->IsInside(Tr.TrEvent[s].SiEnergy,Tr.TrEvent[s].PCEnergy*sin(atan(Tr.TrEvent[s].SiR/(gold_pos-Tr.TrEvent[s].SiZ))))){

	//determine PC position from Silicon position and gold position
	mpc = Tr.TrEvent[s].SiR/(Tr.TrEvent[s].SiZ - gold_pos);
	bpc = Tr.TrEvent[s].SiR - mpc*Tr.TrEvent[s].SiZ;
	Tr.TrEvent[s].pcz_ref = (3.846284509-bpc)/mpc;
	//cout<<"mpc1 = "<<mpc1<<"  bpc1 = "<<bpc1<<"  PCZ_Ref = "<<Tr.TrEvent[s].pcz_ref<<endl;


	//Alternative way aka John Parker's way
	mpc1 = (Tr.TrEvent[s].SiZ - gold_pos)/Tr.TrEvent[s].SiR;
	bpc1 = Tr.TrEvent[s].SiZ - mpc1*Tr.TrEvent[s].SiR;
	Tr.TrEvent[s].PCZ_Ref = mpc1*pcr + bpc1;

	//cout<<"mpc = "<<mpc<<"  bpc = "<<bpc<<"  PCZ_Ref = "<<Tr.TrEvent[s].PCZ_Ref<<endl;

	//cout<<" BeamE = "<<BeamE<<"  pcr = "<<pcr<<endl;
	//cout<<"bpc = "<<bpc<<"  PCZ_Ref = "<<Tr.TrEvent[s].PCZ_Ref<<endl;

	MyFill(Form("PCZ_Ref%i",Tr.TrEvent[s].WireID),300,1,30,Tr.TrEvent[s].pcz_ref);

	MyFill(Form("SiZ_vs_SiEnergy_WireID%i",Tr.TrEvent[s].WireID),600,-1,30,Tr.TrEvent[s].SiEnergy,600,-1,25,Tr.TrEvent[s].SiZ);
	if(Tr.TrEvent[s].DetID<28 && Tr.TrEvent[s].DetID>3)
	  MyFill("SiZ_vs_SiEnergy_R1_R2",600,-1,30,Tr.TrEvent[s].SiEnergy,600,-1,25,Tr.TrEvent[s].SiZ);

	////////////////////////////////////////////////////////////////

	MyFill(Form("PCZ_vs_Z_beforeCal%i",Tr.TrEvent[s].WireID),600,-1.5,1.5,Tr.TrEvent[s].PCZ,600,1.0,30.0,Tr.TrEvent[s].pcz_ref); // before PCWIRECAL applied
	MyFill(Form("PCZ_vs_Z_afterCal%i",Tr.TrEvent[s].WireID),600,-1.0,30.0,Tr.TrEvent[s].PCZ,600,-1.0,30.0,Tr.TrEvent[s].pcz_ref); // after PCWIRECAL applied

	if(Tr.TrEvent[s].DetID<16 && Tr.TrEvent[s].DetID>-1){
	  MyFill(Form("PCZ_vs_Z_afterCal_q3_r1%i",Tr.TrEvent[s].WireID),600,-1.0,30.0,Tr.TrEvent[s].PCZ,600,-1.0,30.0,Tr.TrEvent[s].pcz_ref);
	  MyFill("PCZ_vs_Z_afterCal_q3_r1",600,-1.0,30.0,Tr.TrEvent[s].PCZ,600,-1.0,30.0,Tr.TrEvent[s].pcz_ref);
	}

	if(Tr.TrEvent[s].DetID<28 && Tr.TrEvent[s].DetID>3){
	  MyFill(Form("PCZ_vs_Z_afterCal_r1_r2%i",Tr.TrEvent[s].WireID),600,-1.0,30.0,Tr.TrEvent[s].PCZ,600,-1.0,30.0,Tr.TrEvent[s].pcz_ref);
	}
	
	MyFill("PCZ_vs_Z_afterCal_All",600,-1.0,30.0,Tr.TrEvent[s].PCZ,600,-1.0,30.0,Tr.TrEvent[s].pcz_ref);

	///////////////////////////////////////////////////////////////

	MyFill(Form("PCEnergy_vs_PCZ_beforeCal_WireID%i",Tr.TrEvent[s].WireID),600,-1.5,1.5,Tr.TrEvent[s].PCZ,600,-0.1,0.5,Tr.TrEvent[s].PCEnergy);// before PCWIRECAL applied 
	MyFill(Form("PCEnergy_vs_PCZ_afterCal_WireID%i",Tr.TrEvent[s].WireID),600,-10,60,Tr.TrEvent[s].PCZ,600,-0.1,0.5,Tr.TrEvent[s].PCEnergy);// after PCWIRECAL applied
	
	if(Tr.TrEvent[s].DetID<16 && Tr.TrEvent[s].DetID>-1){
	  MyFill(Form("PCEnergy_vs_PCZ_q3_r1_%i",Tr.TrEvent[s].WireID),600,-10,60,Tr.TrEvent[s].PCZ,600,-0.1,0.5,Tr.TrEvent[s].PCEnergy);
	  MyFill("PCEnergy_vs_PCZ_q3_r1",600,-10,60,Tr.TrEvent[s].PCZ,600,-0.1,0.5,Tr.TrEvent[s].PCEnergy);
	}
	
	MyFill("PCEnergy_vs_PCZ_All",600,-10,60,Tr.TrEvent[s].PCZ,600,-0.1,0.5,Tr.TrEvent[s].PCEnergy);
	

	//--------PCOffset = PCZ-PCZ_ref----------------------//////////////////

	MyFill(Form("PCZoffset_vs_PCEnergy_%i",Tr.TrEvent[s].WireID),600,-0.1,0.5,Tr.TrEvent[s].PCEnergy,600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));

	if(Tr.TrEvent[s].DetID<16 && Tr.TrEvent[s].DetID>-1){
	  MyFill(Form("PCZoffset_vs_PCEnergy_q3_r1_%i",Tr.TrEvent[s].WireID),600,-0.1,0.5,Tr.TrEvent[s].PCEnergy,600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));
	  MyFill("PCZoffset_vs_PCEnergy_q3_r1",600,-0.1,0.5,Tr.TrEvent[s].PCEnergy,600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));
	}

	MyFill("PCZoffset_vs_PCEnergy_All",600,-0.1,0.5,Tr.TrEvent[s].PCEnergy,600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));

	/////////////////////////////////////////////////////////////////////////

	MyFill(Form("PCZoffset_vs_PCZ_beforeCal%i",Tr.TrEvent[s].WireID),600,-1.5,1.5,Tr.TrEvent[s].PCZ,600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref)); // before PCWIRECAL applied
	MyFill(Form("PCZoffset_vs_PCZ_afterCal%i",Tr.TrEvent[s].WireID),600,-1,30,Tr.TrEvent[s].PCZ,600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref)); // after PCWIRECAL applied

	if(Tr.TrEvent[s].DetID<16 && Tr.TrEvent[s].DetID>-1){
	  MyFill(Form("PCZoffset_vs_PCZ_q3_r1_%i",Tr.TrEvent[s].WireID),600,-1,30,Tr.TrEvent[s].PCZ,600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));
	  MyFill("PCZoffset_vs_PCZ_q3_r1",600,-1,30,Tr.TrEvent[s].PCZ,600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));
	}

	MyFill("PCZoffset_vs_PCZ_All",600,-1,30,Tr.TrEvent[s].PCZ,600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));
	
	////////////////////////////////////////////////////////////////////////////

	MyFill(Form("PCZoffset_vs_Theta_%i",Tr.TrEvent[s].WireID),600,0,190,Tr.TrEvent[s].Theta*180/TMath::Pi(),600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));

	if(Tr.TrEvent[s].DetID<16 && Tr.TrEvent[s].DetID>-1){
	  MyFill(Form("PCZoffset_vs_Theta_q3_r1_%i",Tr.TrEvent[s].WireID),600,0,190,Tr.TrEvent[s].Theta*180/TMath::Pi(),600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));
	  MyFill("PCZoffset_vs_Theta_q3_r1",600,0,190,Tr.TrEvent[s].Theta*180/TMath::Pi(),600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));
	}
	MyFill("PCZoffset_vs_Theta_All",600,0,190,Tr.TrEvent[s].Theta*180/TMath::Pi(),600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));

	// for all the wires together the Phi_s make more sense
	MyFill("PCZoffset_vs_PCPhi",600,0,360,Tr.TrEvent[s].PCPhi*180/TMath::Pi(),600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));
	MyFill("PCZoffset_vs_SiPhi",600,0,360,Tr.TrEvent[s].SiPhi*180/TMath::Pi(),600,-20,20,(Tr.TrEvent[s].PCZ-Tr.TrEvent[s].pcz_ref));


	//IntPoint Reconstruction
	if(Tr.TrEvent[s].DetID<4 && Tr.TrEvent[s].DetID>-1)
	  MyFill("IntPoint_cut_Q3",600,-10,60,Tr.TrEvent[s].IntPoint);
	if(Tr.TrEvent[s].DetID<16 && Tr.TrEvent[s].DetID>3)
	  MyFill("IntPoint_cut_SX3_R1",600,-10,60,Tr.TrEvent[s].IntPoint);
	if(Tr.TrEvent[s].DetID<28 && Tr.TrEvent[s].DetID>15)
	  MyFill("IntPoint_cut_SX3_R2",600,-10,60,Tr.TrEvent[s].IntPoint);
	

	  

	if(Tr.NTracks1==1){
	if(Tr.TrEvent[0].WireID >-1 && Tr.TrEvent[0].WireID < 23)
	  MyFill(Form("PCEnergyPlus_cut_Wire%i_Wire%i",Tr.TrEvent[0].WireID+1,Tr.TrEvent[0].WireID),300,0,0.25,PCGoodEnergy[Tr.TrEvent[0].WireID],300,0,0.25,PCGoodEnergy[Tr.TrEvent[0].WireID+1]);
	
	if(Tr.TrEvent[0].WireID >0 && Tr.TrEvent[0].WireID < 24)
	  MyFill(Form("PCEnergyMinus_cut_Wire%i_Wire%i",Tr.TrEvent[0].WireID-1,Tr.TrEvent[0].WireID),300,0,0.25,PCGoodEnergy[Tr.TrEvent[0].WireID],300,0,0.25,PCGoodEnergy[Tr.TrEvent[0].WireID-1]);

      }
	//} //end of cut2 loop

      } // end of for(Int_t s=0; s<Tr.NTracks1;s++){ for PCWireCal

#endif

      /////////////////////////////////////////////////////////////////////////////////////
#ifdef FillEdE_cor
      for(Int_t q=0; q<Tr.NTracks1;q++){
	
	MyFill(Form("PCEnergy_vs_PCZ_WireID%i",Tr.TrEvent[q].WireID),600,-10,60,Tr.TrEvent[q].PCZ,600,-0.1,0.5,Tr.TrEvent[q].PCEnergy);
	MyFill("BeamEnergy_vs_IntPoint",100,-20,80,Tr.TrEvent[q].IntPoint,100,0,90,Tr.TrEvent[q].BeamEnergy);

	MyFill("E_de",600,-1,35,Tr.TrEvent[q].SiEnergy,600,-0.01,0.35,Tr.TrEvent[q].PCEnergy);
	MyFill("E_de_corrected_ALL",600,-1,35,Tr.TrEvent[q].SiEnergy,600,-0.01,0.35,Tr.TrEvent[q].PCEnergy *sin(Tr.TrEvent[q].Theta));
	MyFill("InteractionPoint",300,-10,56,Tr.TrEvent[q].IntPoint);	
	MyFill("E_si_vs_Theta",500,0,200,Tr.TrEvent[q].Theta*ConvAngle,500,0,35,Tr.TrEvent[q].SiEnergy);

	if(Tr.TrEvent[q].DetID<4 && Tr.TrEvent[q].DetID>-1){
	  MyFill("E_de_Q3",600,-1,35,Tr.TrEvent[q].SiEnergy,600,-0.01,0.35,Tr.TrEvent[q].PCEnergy);
	  MyFill("E_de_corrected_Q3",600,-1,35,Tr.TrEvent[q].SiEnergy,600,-0.01,0.35,Tr.TrEvent[q].PCEnergy *sin(Tr.TrEvent[q].Theta));
	  MyFill("InteractionPoint_Q3",300,-10,56,Tr.TrEvent[q].IntPoint);
	  MyFill("E_si_vs_Theta_Q3",500,0,200,Tr.TrEvent[q].Theta*ConvAngle,500,0,35,Tr.TrEvent[q].SiEnergy);
	}
	if(Tr.TrEvent[q].DetID<16 && Tr.TrEvent[q].DetID>3){
	  MyFill("E_de_SX3_1",600,-1,35,Tr.TrEvent[q].SiEnergy,600,-0.01,0.35,Tr.TrEvent[q].PCEnergy);
	  MyFill("E_de_corrected_SX3_1",600,-1,35,Tr.TrEvent[q].SiEnergy,600,-0.01,0.35,Tr.TrEvent[q].PCEnergy *sin(Tr.TrEvent[q].Theta));
	  MyFill("InteractionPoint_SX3_1",300,-10,56,Tr.TrEvent[q].IntPoint);
	  MyFill("E_si_vs_Theta_SX3_1",500,0,200,Tr.TrEvent[q].Theta*ConvAngle,500,0,35,Tr.TrEvent[q].SiEnergy);
	}
	if(Tr.TrEvent[q].DetID<28 && Tr.TrEvent[q].DetID>15){
	  
	  MyFill("E_de_SX3_2",600,-1,35,Tr.TrEvent[q].SiEnergy,600,-0.01,0.35,Tr.TrEvent[q].PCEnergy);
	  MyFill("E_de_corrected_SX3_2",600,-1,35,Tr.TrEvent[q].SiEnergy,600,-0.01,0.35,Tr.TrEvent[q].PCEnergy *sin(Tr.TrEvent[q].Theta));
	  MyFill("InteractionPoint_SX3_2",300,-10,56,Tr.TrEvent[q].IntPoint);
	  MyFill("E_si_vs_Theta_SX3_2",500,0,200,Tr.TrEvent[q].Theta*ConvAngle,500,0,35,Tr.TrEvent[q].SiEnergy);
	}

      }
#endif

      ///////////////////////// DO your analysis here.... /////////////////////////////////////// 
      /*
      Float_t E_4He_rxn1 =0.0;
      Float_t E_4He_rxn2 =0.0;
      Float_t Energy_16O_4He=0.0;
      Float_t Theta_16O_4He=0.0;
      */

      Reconstruct Elastic1(M_18Ne,M_alpha,M_P,M_21Na);
      Elastic1.E_Loss_light = E_Loss_proton;
      // Elastic1.SetELossFile("/home/maria/rayMountPoint/Desktop/anasen_analysis_software/srim_files/H_in_HeCO2_377Torr_18Nerun.eloss");
     for(Int_t c=0; c<Tr.NTracks1; c++){
       //if(Tr.NTracks==Tr.NTracks1){
	if((Tr.TrEvent[c].DetID>-1 && Tr.TrEvent[c].DetID<16) && (Tr.TrEvent[c].BeamEnergy>0 && Tr.TrEvent[c].BeamEnergy<20)){
	  if ( cut2->IsInside(Tr.TrEvent[c].SiEnergy,Tr.TrEvent[c].PCEnergy*sin(Tr.TrEvent[c].Theta)) )
	      {
		MyFill("Timing_Cut",600,1,600,fmod((MCPTime*correct-RFTime),546));
		//Elastic1.ReconstructHeavy(Tr,0);
		Elastic1.ReconstructHeavy_Qvalue(QValue,Tr,c);
		std::cout << Tr.NTracks1 << " " << c  << std::endl;
	      }    
	   }
	//}
     }
      Elastic1.E_Loss_light = NULL;
      

      /*
      Reconstruct Elastic2(M_16O,M_alpha,M_alpha,M_16O);
      Elastic2.SetELossFile("/home/manasta/Desktop/anasen_analysis_software/srim_files/He_in_HeCO2_377Torr_18Nerun.eloss");
      Elastic2.ReconstructElastic(Energy_16O_4He, Theta_16O_4He, Tr, 0);
      */


      //---------------Elastic Scattering-----------------------------/////////////////////////////////////
      /*    
      Float_t M_Beam=M_18Ne;
      Float_t E_4He_rxn=0.0;
      Float_t Energy_4He_elastic=0.0;
      Float_t Energy_Beam_4He=0.0;
      Float_t Theta_Beam_4He=0.0;
      

      for(Int_t c=0; c<Tr.NTracks1;c++){
	 
	if (cut1->IsInside(Tr.TrEvent[c].SiEnergy,Tr.TrEvent[c].PCEnergy*sin(Tr.TrEvent[c].Theta))){

	  MyFill("4He_E_si",500,0,35,Tr.TrEvent[c].SiEnergy);
	  MyFill("4He_E_si_vs_Theta",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,35,Tr.TrEvent[c].SiEnergy);
	  MyFill("4He_E_si_vs_IntPoint",500,0,60,Tr.TrEvent[c].IntPoint,500,0,35,Tr.TrEvent[c].SiEnergy);

	  if(Tr.TrEvent[c].DetID<4 && Tr.TrEvent[c].DetID>-1){
	    MyFill("4He_E_si_Q3",500,0,35,Tr.TrEvent[c].SiEnergy);
	    MyFill("4He_E_si_vs_Theta_Q3",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,35,Tr.TrEvent[c].SiEnergy);
	    MyFill("4He_E_si_vs_IntPoint_Q3",500,0,60,Tr.TrEvent[c].IntPoint,500,0,35,Tr.TrEvent[c].SiEnergy);
	  }else if(Tr.TrEvent[c].DetID<16 && Tr.TrEvent[c].DetID>3){
	    MyFill("4He_E_si_SX3_1",500,0,35,Tr.TrEvent[c].SiEnergy);
	    MyFill("4He_E_si_vs_Theta_SX3_1",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,35,Tr.TrEvent[c].SiEnergy);
	    MyFill("4He_E_si_vs_IntPoint_SX3_1",500,0,60,Tr.TrEvent[c].IntPoint,500,0,35,Tr.TrEvent[c].SiEnergy);
	  }else if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15){
	    MyFill("4He_E_si_SX3_2",500,0,35,Tr.TrEvent[c].SiEnergy);
	    MyFill("4He_E_si_vs_Theta_SX3_2",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,35,Tr.TrEvent[c].SiEnergy);
	    //cout << "IntPoint: " << Tr.TrEvent[c].IntPoint << " La-IntPoint: " << (La-Tr.TrEvent[c].IntPoint) << " BeamEnergy: " << Tr.TrEvent[c].BeamEnergy <<endl;	
	    MyFill("4He_E_si_vs_IntPoint_SX3_2",500,0,60,Tr.TrEvent[c].IntPoint,500,0,35,Tr.TrEvent[c].SiEnergy);
	  }


	  
	  if(Tr.TrEvent[c].SiEnergy>0.0 && Tr.TrEvent[c].SiEnergy<40.0 && Tr.TrEvent[c].PathLength>0.0 && Tr.TrEvent[c].PathLength<100.0 &&
	     Tr.TrEvent[c].IntPoint>0.0 && Tr.TrEvent[c].IntPoint<La && Tr.TrEvent[c].BeamEnergy<BeamE){
	   
	    E_4He_rxn = E_Loss_alpha->GetLookupEnergy(Tr.TrEvent[c].SiEnergy,(-Tr.TrEvent[c].PathLength));
	    Tr.TrEvent[c].LightParEnergy = E_4He_rxn;
	    //E_4He_rxn = E_Loss_alpha->GetInitialEnergy(Tr.TrEvent[c].SiEnergy,Tr.TrEvent[c].PathLength,0.1);
	      MyFill("4He_E_rxn",500,0,40,E_4He_rxn);
	      if(Tr.TrEvent[c].DetID<4 && Tr.TrEvent[c].DetID>-1)
		MyFill("4He_E_rxn_q3",500,0,40,E_4He_rxn);
	      if(Tr.TrEvent[c].DetID<16 && Tr.TrEvent[c].DetID>3)
		MyFill("4He_E_rxn_r1",500,0,40,E_4He_rxn);
	      if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15)
		MyFill("4He_E_rxn_r2",500,0,40,E_4He_rxn);


	      //Beam Energy Calculation from the Elastic scattering: 
	      Energy_Beam_4He = ((M_Beam+M_alpha)*(M_Beam+M_alpha)*E_4He_rxn/(4*M_Beam*M_alpha*cos(Tr.TrEvent[c].Theta)*cos(Tr.TrEvent[c].Theta)));
	      Tr.TrEvent[c].BeamQvalue = Energy_Beam_4He;

	      //Beam Theta Calculation from the Elastic scattering:
	      Theta_Beam_4He = asin((sin(Tr.TrEvent[c].Theta)*sqrt(M_alpha/M_Beam))/sqrt((Energy_Beam_4He/E_4He_rxn)-1))*ConvAngle;
	      Tr.TrEvent[c].ThetaQvalue = Theta_Beam_4He;
	      //MyFill("16O_4He_Energy_VS_ThetaQvalue",300,0,20,Theta_16O_4He,300,0,56,Energy_16O_4He);

	      //4He energy calculated from elastic scattering assuming known SRIM beam energy (not from SiEnergy and Lookup)
	      Energy_4He_elastic = (4*M_alpha*M_Beam*Tr.TrEvent[c].BeamEnergy*cos(Tr.TrEvent[c].Theta)*cos(Tr.TrEvent[c].Theta))/((M_Beam+M_alpha)*(M_Beam+M_alpha));
	      Tr.TrEvent[c].HeEnergyQvalue =  Energy_4He_elastic;

	      ///--------------Alpha particles elastic graphs-----------------------///////////////////////////
	       MyFill("4He_Elastic",500,0,40,Energy_4He_elastic);
	      if(Tr.TrEvent[c].DetID<4 && Tr.TrEvent[c].DetID>-1)
		MyFill("4He_Elastic_q3",500,0,40,Energy_4He_elastic);
	      if(Tr.TrEvent[c].DetID<16 && Tr.TrEvent[c].DetID>3)
		MyFill("4He_Elastic_r1",500,0,40,Energy_4He_elastic);
	      if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15)
		MyFill("4He_Elastic_r2",500,0,40,Energy_4He_elastic);


	      if(Tr.TrEvent[c].DetID<4 && Tr.TrEvent[c].DetID>-1){
		MyFill("LightPar_rxn_vs_Elastic_4He_Energy_Q3",600,0,40,Energy_4He_elastic,600,0,40,E_4He_rxn);
	      }else if(Tr.TrEvent[c].DetID<16 && Tr.TrEvent[c].DetID>3){
		MyFill("LightPar_rxn_vs_Elastic_4He_Energy_SX3_1",600,0,40,Energy_4He_elastic,600,0,40,E_4He_rxn);
	      }else if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15){
		MyFill("LightPar_rxn_vs_Elastic_4He_Energy_SX3_2",600,0,40,Energy_4He_elastic,600,0,40,E_4He_rxn);
	      }
		
	      ///--------------Beam particles elastic graphs-----------------------///////////////////////////
	      
	      MyFill("Beam_4He_Energy",1000,0,100,Energy_Beam_4He);
	      MyFill("BeamEnergy_vs_Beam_4He_Energy",1000,0,100,Energy_Beam_4He,1000,0,100,Tr.TrEvent[c].BeamEnergy);

  
	      if(Tr.TrEvent[c].DetID<4 && Tr.TrEvent[c].DetID>-1){
		MyFill("Beam_4He_Energy_Q3",1000,0,100,Energy_Beam_4He);
		MyFill("BeamEnergy_vs_Beam_4He_Energy_Q3",1000,0,100,Energy_Beam_4He,1000,0,100,Tr.TrEvent[c].BeamEnergy);	     
	      }else if(Tr.TrEvent[c].DetID<16 && Tr.TrEvent[c].DetID>3){
		//cout << "IntPoint: " << Tr.TrEvent[c].IntPoint << " La-IntPoint: " << (La-Tr.TrEvent[c].IntPoint) << " BeamEnergy: " << Tr.TrEvent[c].BeamEnergy <<endl;
		MyFill("Beam_4He_Energy_SX3_1",1000,0,100,Energy_Beam_4He);
		MyFill("BeamEnergy_vs_Beam_4He_Energy_SX3_1",1000,0,100,Energy_Beam_4He,1000,0,100,Tr.TrEvent[c].BeamEnergy);	     
	      }else if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15){
		//cout << "IntPoint: " << Tr.TrEvent[c].IntPoint << " La-IntPoint: " << (La-Tr.TrEvent[c].IntPoint) << " BeamEnergy: " 
		//     << Tr.TrEvent[c].BeamEnergy << " BeamQval: " << Energy_16O_4He << endl;
		MyFill("Beam_4He_Energy_SX3_2",1000,0,100,Energy_Beam_4He);
		MyFill("BeamEnergy_vs_Beam_4He_Energy_SX3_2",1000,0,100,Energy_Beam_4He,1000,0,100,Tr.TrEvent[c].BeamEnergy);	     
	      }

	      
	      /////////// counters and outfiles for check //////////////////////////////////////////
	      
	      if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15 && Energy_4He_elastic<12 && Energy_4He_elastic>5 && E_4He_rxn>5 && E_4He_rxn<12){
		counterBeamE_large_HeR2_special++;
		//outfile_HeR2special << " BeamE: " << Energy_16O_4He << " Tracked_Beam: " << Tr.TrEvent[c].BeamEnergy << " Theta_16O: " << Theta_16O_4He << " Theta_He: " << Tr.TrEvent[c].Theta*ConvAngle << endl
		//<< " He_energy: " << E_4He_rxn  << " HeQvalueEnegy: " << Energy_4He_elastic << endl 
		//<< " SiEnergy: " << Tr.TrEvent[c].SiEnergy  << " IntPoint: " << Tr.TrEvent[c].IntPoint << " SiZ: " << Tr.TrEvent[c].SiZ << " PathLength: " << Tr.TrEvent[c].PathLength << endl
		// << " PCZ: " << Tr.TrEvent[c].PCZ  << " PCEnergy: " << Tr.TrEvent[c].PCEnergy << endl;
				 
	      }
	      
	      

	      if(Tr.TrEvent[c].DetID<4 && Tr.TrEvent[c].DetID>-1 && Energy_16O_4He>BeamE){
		counterBeamE_largeQQQ++;
		//outfileQQQ << " BeamE: " << Energy_16O_4He << " Theta_16O: " << Theta_16O_4He << " Theta_He: " << Tr.TrEvent[c].Theta*ConvAngle << " He_energy: " << E_4He_rxn << " SiEnergy: " << Tr.TrEvent[c].SiEnergy  
		//<< " IntPoint: " << Tr.TrEvent[c].IntPoint << " PCZ: " << Tr.TrEvent[c].PCZ << " PCEnergy: " << Tr.TrEvent[c].PCEnergy  
		//<< " SiZ: " << Tr.TrEvent[c].SiZ   << " SiR: " << Tr.TrEvent[c].SiR << " PathLength: " << Tr.TrEvent[c].PathLength << " Tracked_Beam: " << Tr.TrEvent[c].BeamEnergy << endl;
	      }

	      if(Tr.TrEvent[c].DetID<16 && Tr.TrEvent[c].DetID>3 && Energy_16O_4He>BeamE){
		counterBeamE_largeR1++;
		//outfileR1 << " BeamE: " << Energy_16O_4He << " Theta_16O: " << Theta_16O_4He << " Theta_He: " << Tr.TrEvent[c].Theta*ConvAngle << " He_energy: " << E_4He_rxn << " SiEnergy: " << Tr.TrEvent[c].SiEnergy  
		//<< " IntPoint: " << Tr.TrEvent[c].IntPoint << " PCZ: " << Tr.TrEvent[c].PCZ << " PCEnergy: " << Tr.TrEvent[c].PCEnergy  
		//<< " SiZ: " << Tr.TrEvent[c].SiZ << " PathLength: " << Tr.TrEvent[c].PathLength << " Tracked_Beam: " << Tr.TrEvent[c].BeamEnergy << endl;
	      }


	      if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15 && Energy_16O_4He>BeamE){
		counterBeamE_largeR2++;
		//outfileR2 << " BeamE: " << Energy_16O_4He << " Theta_16O: " << Theta_16O_4He << " Theta_He: " << Tr.TrEvent[c].Theta*ConvAngle << " He_energy: " << E_4He_rxn << " SiEnergy: " << Tr.TrEvent[c].SiEnergy  
		//<< " IntPoint: " << Tr.TrEvent[c].IntPoint << " PCZ: " << Tr.TrEvent[c].PCZ  << " PCEnergy: " << Tr.TrEvent[c].PCEnergy 
		//<< " SiZ: " << Tr.TrEvent[c].SiZ << " PathLength: " << Tr.TrEvent[c].PathLength << " Tracked_Beam: " << Tr.TrEvent[c].BeamEnergy << endl;
	      }


	      if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15 && Tr.TrEvent[c].BeamQvalue<38 && Tr.TrEvent[c].BeamQvalue>18 && Tr.TrEvent[c].BeamEnergy>25 && Tr.TrEvent[c].BeamEnergy<38){
		counterBeamE_largeR2_special++;
		//outfileR2special << " BeamE: " << Energy_16O_4He << " Tracked_Beam: " << Tr.TrEvent[c].BeamEnergy << " Theta_16O: " << Theta_16O_4He << " Theta_He: " << Tr.TrEvent[c].Theta*ConvAngle << endl
		//<< " He_energy: " << E_4He_rxn  << " HeQvalueEnegy: " << Energy_4He_elastic << endl 
		//<< " SiEnergy: " << Tr.TrEvent[c].SiEnergy  << " IntPoint: " << Tr.TrEvent[c].IntPoint << " SiZ: " << Tr.TrEvent[c].SiZ << " PathLength: " << Tr.TrEvent[c].PathLength << endl
		// << " PCZ: " << Tr.TrEvent[c].PCZ  << " PCEnergy: " << Tr.TrEvent[c].PCEnergy << endl;
				 
	      }
	      

	      
	  }	 
	}
      }
      */      

      	
	///------------proton cut-----------------------//////////////////////

      /*   	  
      Float_t E_proton_rxn =0.0;
      Float_t Energy_16O_proton=0.0;
      Float_t Theta_16O_proton=0.0;
      Float_t  Energy_proton_elastic=0.0;
      

      for(Int_t c=0; c<Tr.NTracks1;c++){
	
	if (cut2->IsInside(Tr.TrEvent[c].SiEnergy,Tr.TrEvent[c].PCEnergy*sin(Tr.TrEvent[c].Theta))){

	  MyFill("proton_E_si",500,0,35,Tr.TrEvent[c].SiEnergy);
	  MyFill("proton_E_si_vs_Theta",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,35,Tr.TrEvent[c].SiEnergy);
	  MyFill("proton_E_si_vs_IntPoint",500,0,60,Tr.TrEvent[c].IntPoint,500,0,35,Tr.TrEvent[c].SiEnergy);

	  if(Tr.TrEvent[c].DetID<4 && Tr.TrEvent[c].DetID>-1){
	    MyFill("proton_E_si_Q3",500,0,35,Tr.TrEvent[c].SiEnergy);
	    MyFill("proton_E_si_vs_Theta_Q3",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,35,Tr.TrEvent[c].SiEnergy);
	    MyFill("proton_E_si_vs_IntPoint_Q3",500,0,60,Tr.TrEvent[c].IntPoint,500,0,35,Tr.TrEvent[c].SiEnergy);
	  }else if(Tr.TrEvent[c].DetID<16 && Tr.TrEvent[c].DetID>3){
	    MyFill("proton_E_si_SX3_1",500,0,35,Tr.TrEvent[c].SiEnergy);
	    MyFill("proton_E_si_vs_Theta_SX3_1",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,35,Tr.TrEvent[c].SiEnergy);
	    MyFill("proton_E_si_vs_IntPoint_SX3_1",500,0,60,Tr.TrEvent[c].IntPoint,500,0,35,Tr.TrEvent[c].SiEnergy);
	  }else if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15){
	    MyFill("proton_E_si_SX3_2",500,0,35,Tr.TrEvent[c].SiEnergy);
	    MyFill("proton_E_si_vs_Theta_SX3_2",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,35,Tr.TrEvent[c].SiEnergy);
	    MyFill("proton_E_si_vs_IntPoint_SX3_2",500,0,60,Tr.TrEvent[c].IntPoint,500,0,35,Tr.TrEvent[c].SiEnergy);
	  }

	  if(Tr.TrEvent[c].SiEnergy>0.0 && Tr.TrEvent[c].SiEnergy<40.0 && Tr.TrEvent[c].PathLength>0.0 && Tr.TrEvent[c].PathLength<100.0 &&
	     Tr.TrEvent[c].IntPoint>0.0 && Tr.TrEvent[c].IntPoint<La && Tr.TrEvent[c].BeamEnergy<BeamE){
	   
	    E_proton_rxn = E_Loss_proton->GetLookupEnergy(Tr.TrEvent[c].SiEnergy,(-Tr.TrEvent[c].PathLength));
	    Tr.TrEvent[c].LightParEnergy = E_proton_rxn;
	    
	      MyFill("proton_E_rxn",500,0,40,E_proton_rxn);

	      if(Tr.TrEvent[c].DetID<4 && Tr.TrEvent[c].DetID>-1){
		MyFill("proton_E_q3",500,0,40,E_proton_rxn);
		MyFill("InteractionPoint_protoncut_QQQ",300,-10,56,Tr.TrEvent[c].IntPoint);
	      }
	      if(Tr.TrEvent[c].DetID<16 && Tr.TrEvent[c].DetID>3){
		MyFill("proton_E_r1",500,0,40,E_proton_rxn);
		MyFill("InteractionPoint_protoncut_r1",300,-10,56,Tr.TrEvent[c].IntPoint);
	      }
	      if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15){
		MyFill("proton_E_r2",500,0,40,E_proton_rxn);
		MyFill("InteractionPoint__protoncut_r2",300,-10,56,Tr.TrEvent[c].IntPoint);
	      }

	    //16O beam Energy Calculation from the Elastic scattering: 
	      Energy_16O_proton = ((M_16O+M_P)*(M_16O+M_P)*E_proton_rxn/(4*M_16O*M_P*cos(Tr.TrEvent[c].Theta)*cos(Tr.TrEvent[c].Theta)));
	      Tr.TrEvent[c].BeamQvalue = Energy_16O_proton;

	    //16O beam Theta Calculation from the Elastic scattering:
	      Theta_16O_proton = asin((sin(Tr.TrEvent[c].Theta)*sqrt(M_P/M_16O))/sqrt((Energy_16O_proton/E_proton_rxn)-1))*ConvAngle;
	      Tr.TrEvent[c].ThetaQvalue = Theta_16O_proton;
	      //MyFill("16O_4He_Energy_VS_ThetaQvalue",300,0,20,Theta_16O_4He,300,0,56,Energy_16O_4He);

	      Energy_proton_elastic = (4*M_P*M_16O*Tr.TrEvent[c].BeamEnergy*cos(Tr.TrEvent[c].Theta)*cos(Tr.TrEvent[c].Theta))/((M_16O+M_P)*(M_16O+M_P));
	      Tr.TrEvent[c].HeEnergyQvalue =  Energy_proton_elastic;

	  }
	}

      }
      */
      //////////////////////////////////////////////////////////////////////////////////////////// 
#ifdef FillTree
      MainTree->Fill();
#endif
      //////////////////////////////////////////////////////////////////////////////

   if(Si.NSiHits > 0){
     counter1++;
   }
   else
     {counter4++;}
   
   if(PC.NPCHits > 0){
     counter2++;
   }
   else{counter5++;}
   
   

   if( Tr.NTracks > 0) // total number of events with a Track
     counter7+=Tr.NTracks;
   else{counterZeroTracks++;}
   
   if( Tr.NTracks1 > 0)
     counter9+=Tr.NTracks1;
   else{counterZeroTracks1++;}
   

   if( Tr.NTracks2 > 0)
     counter11+=Tr.NTracks2;
   else{counterZeroTracks2++;}

   if( Tr.NTracks3 > 0)
     counter13+=Tr.NTracks3;
   else{counterZeroTracks3++;}
  
   //cout << "Tr.NTracks= " << Tr.NTracks << " Tr.NTracks1= " << Tr.NTracks1 << " Tr.NTracks2= " << Tr.NTracks2 << " Tr.NTracks3= " << Tr.NTracks3 <<endl;

    }
    ////////////////////////////////////////////////////////////////////////////////   
  }
  //////////////////////////////////////////////////////////////////////////////////
  /*
  outfileQQQ.close(); 
  outfileR1.close(); 
  outfileR2.close(); 
  outfileR2special.close();
  outfile_HeR2special.close();
  */

  cout << endl;
  outputfile->cd();
  RootObjects->Write(); 
  cout << "RootObjects are Written" << endl;
  outputfile->Close();
  cout << "Outputfile Closed\n";
  /*
  cout << "positive SiEnergies OUT GoodPC track1= " << counterPos << endl;
  cout << "negative SiEnergies OUT GoodPC track1= " << counterNeg << endl;
  cout << "positive SiEnergies in GoodPC track1= " << counter0 << endl;
  cout << "negative SiEnergies in GoodPC track1= " << counter << endl;
  
  cout << "positive SiEnergies in GoodPC track2= " << counter10 << endl;
  cout << "negative at 1000 SiEnergies in GoodPC track2= " << counterEn1000 << endl;
  cout << "negative SiEnergies in GoodPC track2= " << counterEnLessZero << endl;

  cout << "positive PCEnergies track3= " << counter12 << endl;
  //cout << "negative PCEnergies track3= " << counter14 << endl;
  cout << "negative at -10 pcEnergies in GoodPC track3= " << counterEnPC1000 << endl;
  cout << "negative pcEnergies in GoodPC track3= " << counterEnPCLessZero << endl;
 
  
  cout << "Si.NSiHits= " << counter1 << endl;
  cout << "SiHits miss= " << counter4 << endl;
  cout << "PC.NPCHits= " << counter2 << endl;
  cout << "PC.NPCHits miss= " << counter5 << endl;
  //cout << "pc_minus_10= " << counterPCMINUSTEN << endl;
  
  cout << "Tr.NTracks= " << counter7 << endl;
  cout << "Tr.NTracks miss= " << counterZeroTracks << endl;
  cout << "Tr.NTracks1= " << counter9 << endl;
  cout << "Tr.NTracks1 miss= " << counterZeroTracks1 << endl;
  cout << "Tr.NTracks2= " << counter11 << endl;
  cout << "Tr.NTracks2 miss= " << counterZeroTracks2 << endl;
  cout << "Tr.NTracks3= " << counter13 << endl;
  cout << "Tr.NTracks3 miss= " << counterZeroTracks3 << endl;
 
  
  cout << "Large BeamEs_QQQ: " << counterBeamE_largeQQQ << endl;
  cout << "Large BeamEs_R1: " << counterBeamE_largeR1 << endl;
  cout << "Large BeamEs_R2: " << counterBeamE_largeR2 << endl;
  cout << "Large BeamEs_HeR2special: " << counterBeamE_large_HeR2_special << endl;
  cout << "Large BeamEs_R2special: " << counterBeamE_largeR2_special << endl;
  */

}//end of Main
////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
Float_t phidiff ( Float_t phi1,Float_t phi2)
{
  return (fmodf((fabs(phi1-phi2) + 2*TMath::Pi()), 2*TMath::Pi()));
}
/////////////////////////////////////////////////////////////////////////////////////
/* 
// Finds a maximum PC within a given phi range
// Nabin Rijal, June 2016

Int_t FindMaxPC(Double_t phi, PCHit& PC){
  Int_t GoodPC = -1;
  Double_t MaxPC = -10;
  //Double_t MinPhi = 0.2619;
  Double_t MinPhi = 0.5238;

  for (int k=0; k<PC.NPCHits; k++){//loop over the pc hits
    //if the PC falls in a range of phi then it is possible correlated
    //we find the maximum energy on the pc
    PC.pc_obj = PC.ReadHit->at(k);

    //if (PC.pc_obj.TrackType == 1){
    //continue;
    //}

    if ( (fabs(PC.pc_obj.PhiW-phi) <= MinPhi) || ((2*TMath::Pi() - fabs(PC.pc_obj.PhiW-phi)) <= MinPhi) ) {
      if ( PC.pc_obj.Energy >= MaxPC ){
	MaxPC = PC.pc_obj.Energy;
	GoodPC = k;
      }
    }
  }
  return GoodPC;
}
*/

///////////////////////////////////////////////////////////////////////////////////
// If there are more than one silicon firing within the range of given phi, 
// picks the one with closer phi and assigns the another pc hit to the next silicon.

// Nabin Rijal, September 18, 2016

Int_t FindMaxPC(Double_t phi, PCHit& PC){
 
  Int_t MaxPCindex = -1,NexttoMaxPCindex =-1;
  Double_t MaxPC = -10, NexttoMaxPC = -10;

  Double_t MinPhi = 0.2619; // 15 degree opening for search  
  //Double_t MinPhi = 0.5238;   // 30 degree opening for search 

  for (int k=0; k<PC.NPCHits; k++){//loop over the pc hits 
    //if the PC falls in a range of phi then it is possible correlated
    //we find the maximum energy on the pc
  
    PC.pc_obj = PC.ReadHit->at(k);   

    if(PC.pc_obj.Energy<0)
      continue;
 
    if ( phidiff(PC.pc_obj.PhiW,phi) < MinPhi) {      

      if ( PC.pc_obj.Energy >= MaxPC ){
	NexttoMaxPC = MaxPC;
	NexttoMaxPCindex = MaxPCindex;
	MaxPC = PC.pc_obj.Energy;
	MaxPCindex = k;

      }
      else if (PC.pc_obj.Energy >= NexttoMaxPC ){
	NexttoMaxPC = PC.pc_obj.Energy;
	NexttoMaxPCindex = k;
      }else{
	continue;
      }
    }    
  }  
  if (NexttoMaxPCindex>0){
    // there is an ambiguity: pick the one with smaller DeltaPhi
    if (phidiff( PC.ReadHit->at(MaxPCindex).PhiW,phi) > phidiff(PC.ReadHit->at(NexttoMaxPCindex).PhiW,phi)){  // it was wrong compared to the Analyzer sent (check email)
      return(NexttoMaxPCindex);
    }
    else
      return(MaxPCindex);      
  }  
  return MaxPCindex; 
}
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX){
  try{
    fhmap.at(name)->Fill(valueX);
  } catch(out_of_range e) {
    TH1F* newHist = new TH1F(name.c_str(),name.c_str(),
			     binsX,lowX,highX);
    newHist->Fill(valueX);
    fhlist->Add(newHist);
    fhmap[name] = newHist;
  }
}
////=================================================================================
void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY){

  try{
    //cout << name << endl;
    fhmap.at(name)->Fill(valueX,valueY);
  } catch(std::out_of_range e) {
    TH2F* newHist = new TH2F(name.c_str(),name.c_str(),
			     binsX,lowX,highX,
			     binsY,lowY,highY);
    newHist->Fill(valueX,valueY);
    fhlist->Add(newHist);
    fhmap[name] = newHist;
  }  
}
/////////////////////////////////////////////////////////////////////////////////////
