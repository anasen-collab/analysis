/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This file creates a tree with vectors, objects and members
// See readme.md for general instructions.
//
// Author: Nabin Rijal, John Parker, Ingo Wiedenhover -- 2016 September.
// Edited by : Jon Lighthall, 2016.12
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define MaxEntries (Long64_t) 1e8
#define MaxADCHits  64
#define MaxTDCHits  500

//Set PC Thresholds here
#define PC_Min_threshold 0  //to forbid Noise
#define PC_Max_threshold 4096 //to forbid Overflow
///////////////////////////////////////////////////////// Switches ////////////////////////////////////////////////////////////

// To select the component of the beam, mostly for radioactive beams
// disable while you work with Calibration data & enable while you do data analysis
//#define Time_hists
//#define MCP_RF_Cut

// beam diagnostic histograms
//#define IC_hists
//#define IC_cut

// Energy for each calibration steps can be switched off after calibration
//#define FillTree_Esteps

// Check pulser calibration?
//#define Pulser_ReRun
#define Re_zeo

// Select the histograms for performing calibration or to check calibration
//#define Hist_for_Si_Cal
//#define Hist_for_PC_Cal
//#define ZPosCal
//#define Hist_after_Cal

///////////////////////////////////////////////////// include Libraries ///////////////////////////////////////////////////////
//C/C++
#include <stdexcept>
#include <map>
#include <fstream>
#include <string>
#include <iomanip>

//ROOT
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TList.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TCutG.h>

//Associated header files/methods
#include "ChannelMap.h"
#include "../include/2016_detclass.h"
#include "Silicon_Cluster.h"

using namespace std;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Methods that fills 1D && 2D Histograms
void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY);
void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX);

//Methods to sort Silicon energy in descending orders
bool Track::Tr_Sisort_method(struct TrackEvent a,struct TrackEvent b){
  if(a.SiEnergy > b.SiEnergy)
    return 1;
  return 0;
};

//Methods to sort PC energy in descending orders
bool Track::Tr_PCsort_method(struct TrackEvent a,struct TrackEvent b){
  if(a.PCEnergy > b.PCEnergy)
    return 1;
  return 0;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TList* fhlist;
std::map<string,TH1*> fhmap;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

  if (argc<3) {
    cout << " Error: Wrong number of arguments\n";
    exit(EXIT_FAILURE);
  }
  
  SiHit Si;
  PCHit PC;  
  
  Silicon_Cluster SiSort;
  SiSort.Initialize();

  //read in command line arguments
  char* filename_histout = new char [200];//output root file
  char* filename_callist = new char [200];//input file list
  
  strcpy( filename_callist, argv[1] );
  strcpy( filename_histout, argv[2] );

  TFile *outputFile = new TFile(filename_histout,"RECREATE");
  TTree *MainTree = new TTree("MainTree","MainTree");
  
  MainTree->Branch("Si.NSiHits",&Si.NSiHits,"NSiHits/I");
  MainTree->Branch("Si.Detector",&Si.Detector);
  MainTree->Branch("Si.Hit",&Si.Hit);
  MainTree->Branch("PC.NPCHits",&PC.NPCHits,"NPCHits/I");
  MainTree->Branch("PC.Hit",&PC.Hit); 

  Int_t RFTime,MCPTime;
  MainTree->Branch("RFTime",&RFTime,"RFTime/I");
  MainTree->Branch("MCPTime",&MCPTime,"MCPTime/I");

#ifdef IC_hists 
  Int_t IC_dE,IC_E;
  MainTree->Branch("IC.dE",&IC_dE,"IC_dE/I");
  MainTree->Branch("IC.E",&IC_E,"IC_E/I");
#endif

  TObjArray *RootObjects = new TObjArray();
  RootObjects->Add(MainTree);

  fhlist = new TList;
  RootObjects->Add(fhlist);

  ChannelMap *CMAP;
  CMAP = new ChannelMap();
  cout<<" ============================================================================================"<<endl; 
  //Initialization of the main channel map  
  ////////////////////////////////////////////////////////////////////////////
  
  //initialize 17F
  if(1) {//load calibration files
    CMAP->Init("Param/24Mg_cals/initialize/ASICS_cmap_022716",
	       "Param/17F_cals/Sipulser_2016.07.20offsets_centroid.dat",
	       "Param/17F_cals/AlphaCal_170515.edit.dat",
	       "Param/17F_cals/X3RelativeGains_Step3_170525.dat",
	       "Param/17F_cals/QQQRelativeGains_Step2_170428.dat");
    CMAP->FinalInit("Param/17F_cals/X3FinalFix_Step3_170525.dat","Param/17F_cals/X3geometry_170502.dat");
    CMAP->LoadQ3FinalFix("Param/17F_cals/QQQFinalFix_Step2_170428.dat");
    CMAP->InitPCCalibration("Param/17F_cals/PCpulserCal_zero_2017-11-06.dat");
    //CMAP->InitPCCalibration("Param/17F_cals/PCpulserCal_zero_offset_2017-11-07.dat");
    //CMAP->InitPCWireCal("Param/17F_cals/PCWireCal_170527_average.dat");
    CMAP->InitPCWireCal("Param/initialize/PCWireCal_init.dat");
    CMAP->Init_PC_UD_RelCal("Param/initialize/PC_UD_RelCal_init.dat");
    CMAP->Init_PCWire_RelGain("Param/initialize/PCWire_RelGain_init.dat");
  }
  else {//load trivial calibration, before any cal where all slopes are one and offsets zero
    CMAP->Init("Param/24Mg_cals/initialize/ASICS_cmap_022716",
	       "Param/initialize/Sipulser_init.dat",
	       "Param/initialize/AlphaCalibration_init.dat",
	       "Param/initialize/X3RelativeGains_Slope1.dat",
	       "Param/initialize/QQQRelativeGains_Slope1.dat");
    CMAP->FinalInit("Param/initialize/X3FinalFix_init.dat","Param/initialize/X3geometry_init.dat");
    CMAP->LoadQ3FinalFix("Param/initialize/QQQFinalFix_init.dat");
    CMAP->InitPCCalibration("Param/initialize/PCpulser_init.dat");
    CMAP->InitPCWireCal("Param/initialize/PCWireCal_init.dat");
    CMAP->Init_PC_UD_RelCal("Param/initialize/PC_UD_RelCal_init.dat");
    CMAP->Init_PCWire_RelGain("Param/initialize/PCWire_RelGain_init.dat");
  }
    
  /*
  //intialize 24Mg
  CMAP->Init("Param/24Mg_cals/initialize/ASICS_cmap_022716",
	     "Param/24Mg_cals/initialize/alignchannels_24Mg_11082016_1262.dat",
	     "Param/initialize/AlphaCalibration_init.dat",
  	     "Param/24Mg_cals/initialize/X3RelativeGains_11022016_Slope1.dat",
	     "Param/24Mg_cals/initialize/QQQRelativeGains11022016_Slope1.dat");

  CMAP->Init("Param/24Mg_cals/initialize/ASICS_cmap_022716",
	     "Param/24Mg_cals/initialize/alignchannels_24Mg_11082016_1262.dat",
	     "Param/initialize/AlphaCalibration_init.dat",
	     "Param/24Mg_cals/SX3Rel/X3RelativeGains_11172016_mix.dat",
	     "Param/24Mg_cals/QQQRel/QQQRelativeGains11092016_Step2.dat");

  //initialize 18Ne
  CMAP->Init("Param/18Ne_cals/ASICS_cmap_06292016",
	     "Param/18Ne_cals/alignchannels_10242016.dat",
	     "Param/initialize/AlphaCalibration_init.dat",
  	     "Param/18Ne_cals/SX3Rel/X3RelativeGains_09182016_Slope1.dat",
	     "Param/18Ne_cals/QQQRel/QQQRelativeGains09122016_Slope1.dat");

  CMAP->Init("Param/18Ne_cals/ASICS_cmap_06292016",
	     "Param/18Ne_cals/alignchannels_10242016.dat",
	     "Param/18Ne_cals/AlphaCal_10312016.dat",
  	     "Param/18Ne_cals/SX3Rel/X3RelativeGains_mix_10262016.dat",
	     "Param/18Ne_cals/QQQRel/QQQRelativeGains10252016_Step2_vol2.dat");
            
  //old cal files for 18Ne
  CMAP->Init("Param/18Ne_cals/ASICS_cmap_06292016",
	     "Param/18Ne_cals/alignchannels_09122016.dat",
	     "Param/18Ne_cals/AlphaCal_09222016.dat",
	     "Param/18Ne_cals/SX3Rel/X3RelativeGains_10052016_step3redo_16_23.dat",
	     "Param/18Ne_cals/QQQRel/QQQRelativeGains09182016_Step2.dat");
	     
  CMAP->Init("Param/18Ne_cals/ASICS_cmap_06292016",
	     "Param/18Ne_cals/alignchannels_09122016.dat",
	     "Param/18Ne_cals/AlphaCal_09222016.dat",
	     "Param/18Ne_cals/SX3Rel/X3RelativeGains09222016_Step3.dat",
	     "Param/18Ne_cals/QQQRel/QQQRelativeGains09182016_Step2.dat");
  
  CMAP->Init("Param/18Ne_cals/ASICS_cmap_06292016",
	     "Param/18Ne_cals/alignchannels_09122016.dat",
	     "Param/initialize/AlphaCalibration_init.dat",
	     "Param/18Ne_cals/SX3Rel/X3RelativeGains_10052016_step3redo_16_23.dat",
	     "Param/18Ne_cals/QQQRel/QQQRelativeGains09122016_Slope1.dat");
  */
  
  CMAP->InitWorldCoordinates("Param/17F_cals/WorldCoord_170223.dat");  
  CMAP->InitPCADC("Param/initialize/NewPCMap");  
  cout<<" ============================================================================================"<<endl;
  //------------------------------------------------------------------------------------------
  TFile *inputFile = new TFile(filename_callist);//open root file and make sure it exists---------------------
  if (!inputFile->IsOpen()) {
    cout << "Root file: " << filename_callist << " could not be opened.\n";
    exit(EXIT_FAILURE);
  }
  
  //Create Objects of the Detector Classes
  ASICHit Si_Old;
  CAENHit ADC;
  CAENHit TDC;
  
  //Set Branch Addresses so that they can be accessed ------------------------
  TTree *input_tree = (TTree*) inputFile->Get("DataTree");
  input_tree->SetBranchAddress("Si.Nhits",&Si_Old.Nhits);
  input_tree->SetBranchAddress("Si.MBID",Si_Old.MBID);
  input_tree->SetBranchAddress("Si.CBID",Si_Old.CBID);
  input_tree->SetBranchAddress("Si.ChNum",Si_Old.ChNum);
  input_tree->SetBranchAddress("Si.Energy",Si_Old.Energy);
  input_tree->SetBranchAddress("Si.Time",Si_Old.Time);
  
  input_tree->SetBranchAddress("ADC.Nhits",&ADC.Nhits);
  input_tree->SetBranchAddress("ADC.ID",ADC.ID);
  input_tree->SetBranchAddress("ADC.ChNum",ADC.ChNum);
  input_tree->SetBranchAddress("ADC.Data",ADC.Data);

  input_tree->SetBranchAddress("TDC.Nhits",&TDC.Nhits);
  input_tree->SetBranchAddress("TDC.ID",TDC.ID);
  input_tree->SetBranchAddress("TDC.ChNum",TDC.ChNum);
  input_tree->SetBranchAddress("TDC.Data",TDC.Data);
  //--------------------------------------------------------------------------

  //Initialize Detector Numbers && Channels
  Int_t DN = -1;
  Int_t DetCh=-1;
  
  Double_t ZeroShift=0,VperCh=0;
  //Double_t p0_coeff = 0,p1_coeff = 0, p2_coeff = 0, x0_offset = 0;
  Double_t Gain_Rel=0,Gain_Alpha=0;
  Double_t FinalShift=0;
  
  Double_t SX3Energy[NumSX3][MaxSX3Ch];
  Double_t SX3Time[NumSX3][MaxSX3Ch];
  Double_t SX3Energy_Pulser[NumSX3][MaxSX3Ch];
  Double_t SX3Energy_Rel[NumSX3][MaxSX3Ch];
  Double_t SX3Energy_Cal[NumSX3][MaxSX3Ch];

  Double_t Q3Energy[NumQ3][MaxQ3Ch];
  Double_t Q3Time[NumQ3][MaxQ3Ch];
  Double_t Q3Energy_Pulser[NumQ3][MaxQ3Ch];
  Double_t Q3Energy_Rel[NumQ3][MaxQ3Ch];
  Double_t Q3Energy_Cal[NumQ3][MaxQ3Ch];
  
  //PC Variables
  Int_t Side,WireID;
  Double_t PCDown[NPCWires];
  Double_t PCUp[NPCWires];
  Double_t PCDownVoltage[NPCWires];
  Double_t PCUpVoltage[NPCWires];
  Bool_t ConvTest=kFALSE;
  Double_t Vcal;

  Int_t status = 0;
  Long64_t nentries = input_tree->GetEntries();
  cout << " nentries = " << nentries<<"  in  "<< filename_callist <<endl;
  Bool_t btrunc=kFALSE;
  Long64_t nstep=nentries/MaxEntries;
  //Long64_t ncount=0;

  if(nentries>MaxEntries) {
    cout << " Max entries exceeded! truncating data set from " << nentries << " to " << MaxEntries << " or " << nentries/(nstep) <<endl;
    cout << " processing 1 out of every " << nstep << " entries" <<endl;
    btrunc=kTRUE;
    nentries=MaxEntries;
  }
  Float_t print_step=0.1;
  if(nentries>5e5)
    print_step/=10;
  
  for (Long64_t global_evt=0; global_evt<nentries; global_evt++) {//loop over all entries in tree------
    if(btrunc && global_evt%nstep>0) continue;
    status = input_tree->GetEvent(global_evt);
    if(global_evt%TMath::Nint(nentries*print_step)==0) { cout << endl << "  Done: "
						      << right << fixed << setw(3)
						      << TMath::Nint(global_evt*100./nentries) << "%" << std::flush;
      //cout << " global_evt = " << global_evt << ", ncount = " <<ncount<<endl;
    }
    if(global_evt%TMath::Nint(nentries*print_step/10)==0 && global_evt>0) cout << "." << std::flush;
    //std::cout << "\rDone: " << global_evt*100./nentries << "%          ";// << std::flush;
   
    /////////////////////////////////  CAEN section (PC, IC, CsI,..etc) ////////////////////////////////

    //======================= PC variables are initialized and Filled here.=============================

    PC.zeroPCHit();  

    Int_t Identifier=-1;  
    //************************************************************************************************
    // Identifier provides information on the detector type recorded by ADC.
    // 0 - no identification was possible, this channel is not described in the ADCChannels.xxxx file
    // 1 - Gas Proportional counter data
    // 2 - CsI(Tl) scintillator   
    //*************************************************************************************************
    // Initialize variables named above.
    for (Int_t i=0; i<NPCWires; i++) {
      PCDown[i]  = sqrt(-1);
      PCUp[i]    = sqrt(-1);
      PCDownVoltage[i] = sqrt(-1);
      PCUpVoltage[i]   = sqrt(-1);
    }
    // Make sure your ADC.Nhits is within bounds.
    if (ADC.Nhits>MaxADCHits) {
      printf("MaxADCHits exceeded! %d > %d\n",ADC.Nhits,MaxADCHits);
      ADC.Nhits=MaxADCHits;
    }
    
    for (Int_t n=0; n<ADC.Nhits; n++) {
      // Identify which detector type we have a hit in.
      CMAP->IdentifyADC(ADC.ID[n],ADC.ChNum[n],Identifier);
      switch (Identifier) {
      case 0:
	break;
      case 1:
	CMAP->IdentifyWire(ADC.ID[n],ADC.ChNum[n],WireID,Side);
	ConvTest =  CMAP->ConvertToVoltage(ADC.ID[n],ADC.ChNum[n],ADC.Data[n],Vcal);
	if (Side==1) {
	  PCDown[WireID] = (Double_t)ADC.Data[n];
	  if (ConvTest) PCDownVoltage[WireID] = Vcal;
	}
	if (Side==2) {
	  PCUp[WireID]   = (Double_t)ADC.Data[n];
	  if (ConvTest) PCUpVoltage[WireID]  = Vcal;
	}
	ConvTest = kFALSE;
	break;
      default:
	break;
      }// End switch
    }// End loop over ADC.Nhits
  
    //=================================
    Double_t XWPC,YWPC,ZWPC,RWPC,PhiWPC,PCRelGain, SlopeUD, OffsetUD;
   
    for ( Int_t i=0; i<NPCWires; i++ ) {

      PC.ZeroPC_obj();

      if ( PCDown[i] > PC_Min_threshold && PCUp[i] > PC_Min_threshold && PCDown[i] < PC_Max_threshold && PCUp[i] < PC_Max_threshold ) { 

	PC.pc_obj.WireID = i;
	PC.pc_obj.Down = PCDown[i];
	PC.pc_obj.Up = PCUp[i];
	PC.pc_obj.Sum=PC.pc_obj.Down+PC.pc_obj.Up;
	PC.pc_obj.DownVoltage = PCDownVoltage[i];
	PC.pc_obj.UpVoltage = PCUpVoltage[i];

#ifdef Re_zero
	Float_t rezero=1.53544582426548004e-02;
	PC.pc_obj.DownVoltage += rezero/2; 
	PC.pc_obj.UpVoltage += rezero/2;
#endif
	
	PC.pc_obj.SumVoltage=PC.pc_obj.DownVoltage+PC.pc_obj.UpVoltage;

	Int_t bins=512;
	Float_t vmin=-0.1;
	Float_t vmax=0.2;
	Float_t orange=0.008;
	
#ifdef Hist_for_PC_Cal	
	MyFill(Form("PC_Down_vs_Up_BeforeCal_Wire%i",i),
	       bins,vmin,vmax,PC.pc_obj.UpVoltage,bins,vmin,vmax,PC.pc_obj.DownVoltage);
	MyFill(Form("PC_Offset_vs_Down_BeforeCal_Wire%i",i),
	       bins,vmin,vmax,PC.pc_obj.DownVoltage,
	       bins,-orange,orange,(PC.pc_obj.DownVoltage-PC.pc_obj.UpVoltage));
	MyFill(Form("PC_sum_vs_diff%i",i),
	       bins,-vmax,vmax,(PC.pc_obj.DownVoltage-PC.pc_obj.UpVoltage),
	       bins,vmin,2*vmax,PC.pc_obj.SumVoltage);
#endif

	CMAP->Get_PC_UD_RelCal(i, SlopeUD, OffsetUD);
	PC.pc_obj.DownRel = (PC.pc_obj.DownVoltage/SlopeUD) - OffsetUD;//PCDownRel[i];
	PC.pc_obj.UpRel = PC.pc_obj.UpVoltage;//PCUpRel[i];
	
#ifdef Hist_after_Cal	
	MyFill(Form("PC_Up_vs_Down_AfterCal_Wire%i",i),
	       bins,vmin,vmax,PC.pc_obj.UpRel,bins,vmin,vmax,PC.pc_obj.DownRel);
	MyFill(Form("PC_Offset_vs_Down_AfterCal_Wire%i",i),
	       bins,vmin,vmax,PC.pc_obj.DownVoltage,
	       bins,-orange,orange,(PC.pc_obj.DownVoltage-PC.pc_obj.UpVoltage));
#endif

	if (PC.pc_obj.DownRel>0 && PC.pc_obj.UpRel>0) {
	  PC.pc_obj.Energy = PC.pc_obj.DownRel + PC.pc_obj.UpRel;
	  PC.pc_obj.Z = (PC.pc_obj.UpRel - PC.pc_obj.DownRel)/PC.pc_obj.Energy;	 
	    
	  CMAP->Get_PCWire_RelGain(PC.pc_obj.WireID, PCRelGain);
	  //cout<<"  PCRelGain == "<<PCRelGain<<endl;
	  PC.pc_obj.Energy *= PCRelGain;

	  Float_t zrange=1.2;

#ifdef Hist_for_PC_Cal		  
	  MyFill(Form("PC_sum_vs_Z%i",PC.pc_obj.WireID),
		 bins,-zrange,zrange,PC.pc_obj.Z,
		 bins,vmin,2*vmax,PC.pc_obj.Energy);
	  MyFill(Form("PCZ%i",PC.pc_obj.WireID),
		 bins,-zrange,zrange,PC.pc_obj.Z);
#endif	

	  // calculate world coordinates
	  CMAP->GetPCWorldCoordinates(PC.pc_obj.WireID,PC.pc_obj.Z,XWPC,YWPC,ZWPC,RWPC,PhiWPC);
	  PC.pc_obj.XW = XWPC;
	  PC.pc_obj.YW = YWPC;
	  PC.pc_obj.ZW = ZWPC;
	  PC.pc_obj.RW = RWPC;
	  PC.pc_obj.PhiW = PhiWPC;
	}
	PC.Hit.push_back(PC.pc_obj);
	PC.NPCHits++;
      }
    }//end NPCWires loops
    ////=============================== MCP && RF =====================================================
    if (TDC.Nhits>MaxTDCHits) {
      printf("MaxTDCHits exceeded! %d > %d\n",TDC.Nhits,MaxTDCHits);
      TDC.Nhits=MaxTDCHits;
    }

    RFTime = 0;
    MCPTime = 0;
 
    Double_t slope=0.9852; //slope of MCP vs RF
    Double_t offset=271.58; //peak-to-peak spacing
    Double_t wrap=offset*2;
    Double_t TOF,TOFc,TOFw;
    Int_t tbins=600;

    for (Int_t n=0; n<TDC.Nhits; n++) {     
      if(TDC.ID[n] == 12 && TDC.ChNum[n]==0) {
	RFTime = (Int_t)TDC.Data[n];
	if(RFTime > 0) {
#ifdef Time_hists	  
	  MyFill("Time_RF",1028,0,4096,RFTime);
#endif
	}
      }
      if(TDC.ID[n] == 12 && TDC.ChNum[n]==7) {
	MCPTime = (Int_t)TDC.Data[n];
	if(MCPTime >0) {
#ifdef Time_hists
	  MyFill("Time_MCP",1028,0,4096,MCPTime);
#endif
	}
      }
    }
    if(RFTime >0 && MCPTime >0) {
      TOF=MCPTime-RFTime;//Time-of-flight
      TOFc=MCPTime-slope*RFTime;//corrected TOF
      TOFw=fmod(TOFc+4*offset,offset);//wrapped TOF
#ifdef Time_hists
      MyFill("Time_MCP_vs_RF",512,0,4096,RFTime,512,0,4096,MCPTime);
      MyFill("TOF_vs_RF",512,0,4096,RFTime,512,-4096,4096,TOF);
      MyFill("TOFc_vs_RF",512,0,4096,RFTime,512,-4096,4096,TOFc);      
      MyFill("TOFc",tbins*2,-4096,4096,TOFc);
      MyFill("TOFw",tbins,0,300,TOFw);
      MyFill("TOFw2",tbins,0,600,fmod(TOFc+4*offset,wrap));//wrapped TOF
#endif
	//=========================== MCP - RF Gate =================================================
#ifdef MCP_RF_Cut    
      if(TOFw>120 && TOFw<178) {//inside gate; keep
#ifdef Time_hists
	MyFill("TOF_wrapped_in",tbins,0,300,TOFw);
#endif
      }
      else {//outside gate; exclude
#ifdef Time_hists
	MyFill("TOF_wrapped_out",tbins,0,300,TOFw);
#endif
	continue;
      }
    }
    else {//bad time; exclude
      continue;
#endif
    }
    
#ifdef IC_hists       
    //------------Ion Chamber----------------------
    IC_dE = 0; IC_E = 0;
    for(Int_t n=0; n<ADC.Nhits; n++) {
      if(ADC.ID[n]==3 && ADC.ChNum[n]==24) {
	IC_dE = ADC.Data[n];
	if(IC_dE >0)
	  MyFill("IC_dE",1028,0,4096,IC_dE);
      }
      if(ADC.ID[n]==3 && ADC.ChNum[n]==28) {
	IC_E = (Int_t)ADC.Data[n];
	if(IC_E >0) {
	  MyFill("IC_E",1028,0,4096,IC_E);
	  if(RFTime >0 && MCPTime >0) {
	    MyFill("IC_TOF_vs_ESi",512,0,4096,IC_E,512,0,4096,TOF);
	    MyFill("IC_TOFc_vs_ESi",512,0,4096,IC_E,512,0,4096,TOFc);
	  }
	}
      }
    }
    if(IC_dE >0 && IC_E >0)
      MyFill("IC_EdE",512,0,4096,IC_E,512,0,4096,IC_dE);
    
#ifdef IC_cut
    char* file_cut1 = "/home/lighthall/anasen/root/main/17F_cut.root";
    TFile *cut_file1 = new TFile(file_cut1);
    if (!cut_file1->IsOpen()){
      cout << "Cut file1: " << file_cut1 << " could not be opened.\n";
      exit(EXIT_FAILURE);
    }
    TCutG *cut1 = NULL;
    cut1 = (TCutG*)cut_file1->Get("CUTG");
    
    if (cut1 == NULL){
      cout << "cut1 does not exist\n";
      exit(EXIT_FAILURE);
    }
#endif
#endif     
    
    /////////////////////////////////////////////  ASICS Section //////////////////////////////////////////////////////
    // ==========================  ASICS  variables are initialized here ====================

    Si.zeroSiHit();  

    for (int i=0; i<NumQ3; i++){
      for (int j=0; j<MaxQ3Ch; j++){
	Q3Energy[i][j] = 0;
	Q3Energy_Pulser[i][j] = 0;
	Q3Energy_Rel[i][j] = 0;
	Q3Energy_Cal[i][j] = 0;
	Q3Time[i][j] = 0;
      }
    }
    for (int i=0; i<NumSX3; i++){
      for (int j=0; j<MaxSX3Ch; j++){
	SX3Energy[i][j] = 0;
	SX3Energy_Pulser[i][j] = 0;
	SX3Energy_Rel[i][j] = 0;
	SX3Energy_Cal[i][j] = 0;
	SX3Time[i][j] = 0;
      }
    }
    //////////////////////////////////////////////  Fill ASICS  ///////////////////////////////////////////////////////

    // Make sure there are not too many hits.
    if (Si_Old.Nhits>MaxSiHits) Si_Old.Nhits=MaxSiHits;

    for (Int_t n=0; n<Si_Old.Nhits; n++){//loop over all Si hits--very similar to Main.C from Jeff
      ZeroShift  = 0;
      VperCh     = 1;
      //p0_coeff = 0;
      //p1_coeff = 0;
      //p2_coeff = 1;
      //x0_offset = 0;

      Gain_Rel  = 1;
      Gain_Alpha = 1;
      FinalShift = 0;    

      CMAP->IdentifyDetChan(Si_Old.MBID[n],Si_Old.CBID[n],Si_Old.ChNum[n], DN, DetCh);
      CMAP->AlignASICsChannels(Si_Old.MBID[n],Si_Old.CBID[n],Si_Old.ChNum[n], ZeroShift, VperCh);
      
      //added by M.Anastasiou
      //CMAP->AlignASICsChannels_Quadratic(Si_Old.MBID[n],Si_Old.CBID[n],Si_Old.ChNum[n], p2_coeff, p1_coeff, p0_coeff, x0_offset);

      //Check the Detector Numbers
      if(DN<0 || DN>27){continue;}

     
      if(DN>3){ //For SX3's
	CMAP->GetSX3MeVPerChannel1(DN,DetCh,Gain_Rel);  
	CMAP->GetSX3MeVPerChannel2(DN,DetCh,Gain_Alpha);  
	CMAP->GetSX3FinalEnergyOffsetInMeV(DN,DetCh,FinalShift);

	SX3Energy[DN-4][DetCh] = (Double_t)Si_Old.Energy[n];
	SX3Energy_Pulser[DN-4][DetCh] = SX3Energy[DN-4][DetCh]+ZeroShift/VperCh;
	//SX3Energy_Pulser[DN-4][DetCh] = SX3Energy[DN-4][DetCh]*VperCh+ZeroShift;   //pulser CHECK

	//in case of quadratic fit
	//SX3Energy_Pulser[DN-4][DetCh] = p2_coeff*(pow(SX3Energy[DN-4][DetCh],2)) + p1_coeff*SX3Energy[DN-4][DetCh] + p0_coeff;
	
	SX3Energy_Rel[DN-4][DetCh] = SX3Energy_Pulser[DN-4][DetCh]*Gain_Rel;
	SX3Energy_Rel[DN-4][DetCh] += FinalShift;
	SX3Energy_Cal[DN-4][DetCh] = SX3Energy_Rel[DN-4][DetCh]*Gain_Alpha;

	SX3Time[DN-4][DetCh]   = (Double_t)Si_Old.Time[n];

      }else{ //For Q3's
	CMAP->GetQ3MeVPerChannel1(DN,DetCh,Gain_Rel);
	CMAP->GetQ3MeVPerChannel2(DN,DetCh,Gain_Alpha);
	CMAP->GetQ3FinalEnergyOffsetInMeV(DN,DetCh,FinalShift);

	Q3Energy[DN][DetCh] = (Double_t)Si_Old.Energy[n];
	Q3Energy_Pulser[DN][DetCh] = Q3Energy[DN][DetCh] + ZeroShift/VperCh;
	//Q3Energy_Pulser[DN][DetCh] = Q3Energy[DN][DetCh]*VperCh + ZeroShift;

	//in case of quadratic fit
	//Q3Energy_Pulser[DN][DetCh] = p2_coeff*(pow(Q3Energy[DN][DetCh],2)) + p1_coeff*Q3Energy[DN][DetCh] + p0_coeff;

	Q3Energy_Rel[DN][DetCh] = Q3Energy_Pulser[DN][DetCh]*Gain_Rel;
	Q3Energy_Rel[DN][DetCh] += FinalShift;
	Q3Energy_Cal[DN][DetCh] = Q3Energy_Rel[DN][DetCh]*Gain_Alpha;
	Q3Time[DN][DetCh] = (Double_t)Si_Old.Time[n];
      }
    }// End loop over Si_Nhits-----------------------------------------------------------------------------------------------
    //////////////////////////////////////////// Push back SX3 Detector-members///////////////////////////////////////////////

    for (int i=0; i<NumSX3; i++){
     
      Si.ZeroSi_obj();

      for (int j=0; j<MaxSX3Ch; j++){

	if ( (SX3Energy_Cal[i][j] > 0.0)){

	  if (j<4){//SX3 Back

	    Si.det_obj.BackChNum.push_back(j);
#if defined(FillTree_Esteps) || defined(Hist_for_Si_Cal)
	    Si.det_obj.EBack_Raw.push_back(SX3Energy[i][j]);
	    Si.det_obj.EBack_Pulser.push_back(SX3Energy_Pulser[i][j]);
	    Si.det_obj.EBack_Rel.push_back(SX3Energy_Rel[i][j]);
#endif	   
	    Si.det_obj.EBack_Cal.push_back(SX3Energy_Cal[i][j]);
	    Si.det_obj.TBack.push_back(SX3Time[i][j]);	


	    //Calculate && push back SX3_ZUp & SX3_ZDown
	    /* if(SX3Energy_Cal[i][j+4]>0){
	      Si.det_obj.SX3_ZUp.push_back(1-(2*SX3Energy_Cal[i][j+4]/SX3Energy_Cal[i][j]));
	    }else if(SX3Energy_Cal[i][j+8]>0){	  
	      Si.det_obj.SX3_ZDown.push_back((2*SX3Energy_Cal[i][j+8]/SX3Energy_Cal[i][j])-1);	         
	    }else{
	      break;
	      }*/

	  }else if(j>3 && j<8){//SX3 Front Up
	    Si.det_obj.UpChNum.push_back(j-4);
#if defined(FillTree_Esteps) || defined(Hist_for_Si_Cal)
	    Si.det_obj.EUp_Raw.push_back(SX3Energy[i][j]);
	    Si.det_obj.EUp_Pulser.push_back(SX3Energy_Pulser[i][j]);
	    Si.det_obj.EUp_Rel.push_back(SX3Energy_Rel[i][j]);
#endif	  	   
	    Si.det_obj.EUp_Cal.push_back(SX3Energy_Cal[i][j]);
	    Si.det_obj.TUp.push_back(SX3Time[i][j]);

	  }else if(j>7 && j<12){//SX3 Front Down
	    Si.det_obj.DownChNum.push_back(j-8);
#if defined(FillTree_Esteps) || defined(Hist_for_Si_Cal)
	    Si.det_obj.EDown_Raw.push_back(SX3Energy[i][j]);
	    Si.det_obj.EDown_Pulser.push_back(SX3Energy_Pulser[i][j]);
	    Si.det_obj.EDown_Rel.push_back(SX3Energy_Rel[i][j]);
#endif
	    Si.det_obj.EDown_Cal.push_back(SX3Energy_Cal[i][j]);
	    Si.det_obj.TDown.push_back(SX3Time[i][j]);
	  }
	}	
	//===========================================
      }//end of for (int j=0; j<MaxSX3Ch; j++){
      //=========================================== 

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if ( (Si.det_obj.EDown_Cal.size()!=0 || Si.det_obj.EUp_Cal.size()!=0) && Si.det_obj.EBack_Cal.size()!=0 ) {

	Si.det_obj.UpMult = Si.det_obj.EUp_Cal.size();
	Si.det_obj.DownMult = Si.det_obj.EDown_Cal.size();
	Si.det_obj.BackMult = Si.det_obj.EBack_Cal.size();
	Si.det_obj.DetID = i+4;

	Si.det_obj.HitType = Si.det_obj.BackMult*100 + Si.det_obj.UpMult*10 + Si.det_obj.DownMult;
		
	SiSort.SortSX3(&Si,CMAP);
	Si.Detector.push_back(Si.det_obj);
	

	//////////////////////////////////////////// Fill SX3 Histograms after Calibration////////////////////////////////////////////
#ifdef Hist_after_Cal
	MyFill(Form("back_vs_front_Cal%i",Si.hit_obj.DetID),500,0,30,Si.hit_obj.EnergyFront,500,0,30,Si.hit_obj.EnergyBack);
	//MyFill(Form("back_vs_front_Cal%i_f%i",Si.hit_obj.DetID,Si.hit_obj.FrontChannel),100,0,30,Si.hit_obj.EnergyFront,100,0,30,Si.hit_obj.EnergyBack);
	//MyFill(Form("back_vs_front_Cal%i_b%i",Si.hit_obj.DetID,Si.hit_obj.BackChannel),100,0,30,Si.hit_obj.EnergyFront,100,0,30,Si.hit_obj.EnergyBack);
	//MyFill(Form("back_vs_front_Cal%i_%i_%i",Si.hit_obj.DetID,Si.hit_obj.FrontChannel,Si.hit_obj.BackChannel),100,0,30,Si.hit_obj.EnergyFront,100,0,30,Si.hit_obj.EnergyBack);
	if(Si.det_obj.HitType ==111){//Requires both Up and Down signal	
	  MyFill(Form("sx3offset_back_vs_front_Cal%i",Si.hit_obj.DetID),960,0,30,Si.hit_obj.EnergyBack,640,-10,10,(Si.hit_obj.EnergyBack-Si.hit_obj.EnergyFront));
	  //MyFill(Form("down_vs_up_Cal%i",Si.det_obj.DetID),100,0,30,Si.det_obj.EUp_Cal[0],100,0,30,Si.det_obj.EDown_Cal[0]);
	  //MyFill(Form("down_vs_up_Cal%i_f%i",Si.det_obj.DetID,Si.det_obj.UpChNum[0]),100,0,30,Si.det_obj.EUp_Cal[0],100,0,30,Si.det_obj.EDown_Cal[0]);
	}
#endif
	//////////////////////////////////////////// Fill SX3 Histograms for Calibration////////////////////////////////////////////
#ifdef Hist_for_Si_Cal
	Int_t udmax=4*4096/3;
	Int_t fbmax=4*4096/3;
	Int_t bins=512;
	if(Si.det_obj.HitType ==111) {//Requires both Up and Down signal	----- Down vs Up histo needs it //Back vs front will be simpler

	  // Step 1 RelCal/U-D, all energies changed to E_Rel
	  MyFill(Form("down_vs_up%i_%i_%i",Si.det_obj.DetID,Si.det_obj.UpChNum[0],Si.det_obj.BackChNum[0]),
		 bins,0,udmax, Si.det_obj.EUp_Rel[0],bins,0,udmax,Si.det_obj.EDown_Rel[0]);
	  MyFill(Form("down_vs_up_divB%i_%i_%i",Si.det_obj.DetID,Si.det_obj.UpChNum[0],Si.det_obj.BackChNum[0]),
		 bins,0,1.3, Si.det_obj.EUp_Rel[0]/Si.det_obj.EBack_Rel[0],
		 bins,0,1.3,Si.det_obj.EDown_Rel[0]/Si.det_obj.EBack_Rel[0]);
	  MyFill(Form("down_vs_up%i_f%i",Si.det_obj.DetID,Si.det_obj.UpChNum[0]),
		 bins,0,udmax, Si.det_obj.EUp_Rel[0],bins,0,udmax,Si.det_obj.EDown_Rel[0]);
	  MyFill(Form("down_vs_up_divB%i_f%i",Si.det_obj.DetID,Si.det_obj.UpChNum[0]),
		 bins,0,1.3,(Si.det_obj.EUp_Rel[0]/Si.det_obj.EBack_Rel[0]),
		 bins,0,1.3,(Si.det_obj.EDown_Rel[0]/Si.det_obj.EBack_Rel[0]));

	  // Step 2 RelCal//F-B //Condition: RelGain Cal from Up-Down is applied
	  MyFill(Form("back_vs_front%i_%i_%i",Si.det_obj.DetID,Si.det_obj.UpChNum[0],Si.det_obj.BackChNum[0]),
		 bins,0,fbmax,Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0],
		 bins,0,fbmax,Si.det_obj.EBack_Rel[0]);

	  // Step 3 RelCal//F-B //RelGain Cal from Step 2 is applied
	  MyFill(Form("back_vs_front%i_b%i",Si.det_obj.DetID,Si.det_obj.BackChNum[0]),bins,0,fbmax,Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0],bins,0,fbmax,Si.det_obj.EBack_Rel[0]);
	  
	  //// just for checking histograms per detector
	  MyFill(Form("back_vs_front%i",Si.det_obj.DetID),bins,0,fbmax,Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0],bins,0,fbmax,Si.det_obj.EBack_Rel[0]);
	  MyFill(Form("down_vs_up%i",Si.det_obj.DetID),600,0,6000,Si.det_obj.EUp_Rel[0],600,0,6000,Si.det_obj.EDown_Rel[0]);
	  Int_t obins=300;
	  Int_t omax=400;
	  
	  // check offset
  	  //MyFill(Form("front_vs_offset%i",Si.det_obj.DetID),
	  // 	 obins,-omax,omax,(Si.det_obj.EBack_Rel[0]-(Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0])),
	  // 	 bins,0,fbmax,Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0]
	  // 	 );
	  
	  MyFill(Form("back_vs_offset_divB%i",Si.det_obj.DetID),
		 obins,-0.2,0.2,((Si.det_obj.EBack_Rel[0]-(Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0]))/Si.det_obj.EBack_Rel[0]),
		 bins,0,fbmax,Si.det_obj.EBack_Rel[0]
		 );

	  MyFill(Form("back_vs_offset%i_%i_%i",Si.det_obj.DetID,Si.det_obj.UpChNum[0],Si.det_obj.BackChNum[0]),
		 obins,-omax,omax,(Si.det_obj.EBack_Rel[0]-(Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0])),
		 bins,0,fbmax,Si.det_obj.EBack_Rel[0]);
	  
	  MyFill(Form("back_vs_offset%i" ,Si.det_obj.DetID),
		 obins,-omax,omax,(Si.det_obj.EBack_Rel[0]-(Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0])),
		 bins,0,fbmax,Si.det_obj.EBack_Rel[0]
		 );

	  // position
	  //MyFill(Form("back_vs_pos2%i",Si.det_obj.DetID),//position over [-1,1]
	  // 	 obins,-1,1,(Si.det_obj.EDown_Rel[0]-Si.det_obj.EUp_Rel[0])/Si.det_obj.EBack_Rel[0],
	  // 	 //obins,-1,1,(Si.det_obj.EDown_Rel[0]-Si.det_obj.EUp_Rel[0])/(Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0])
	  // 	 bins,0,fbmax,Si.det_obj.EBack_Rel[0]
	  // 	 );

	  MyFill(Form("back_vs_pos%i",Si.det_obj.DetID),//position over [0,1]
		 obins,-0.1,1.1,(1./2)*(1+(Si.det_obj.EDown_Rel[0]-Si.det_obj.EUp_Rel[0])/Si.det_obj.EBack_Rel[0]),
		 bins,0,fbmax,Si.det_obj.EBack_Rel[0]
		 );
	}  
#endif	 	  
	//////////////////////////////////////////// Fill SX3 Histograms for Z-Position Calibration//////////////////////////////////
#ifdef ZPosCal
	///////////////////  ZPosCal from the raw data from the Detector  ///////////////////
#ifdef Hist_for_PC_cal
	if(Si.det_obj.HitType ==111) {//Requires both Up and Down signal 
	  if(Si.det_obj.EUp_Cal[0] > 0 && (Si.det_obj.EUp_Cal[0] >= Si.det_obj.EDown_Cal[0])) {
	    //if(Si.det_obj.SX3_ZUp[0] >= Si.det_obj.SX3_ZDown[0]) {	  
	    MyFill(Form("SX3Zpos_%i_%i_%i",Si.det_obj.DetID,Si.det_obj.UpChNum[0],Si.det_obj.BackChNum[0]),100,-1,1,Si.det_obj.SX3_ZUp[0]);
	  }else if(Si.det_obj.EDown_Cal[0] >0 && (Si.det_obj.EUp_Cal[0] < Si.det_obj.EDown_Cal[0])) {
	    //}else if(Si.det_obj.SX3_ZUp[0]< Si.det_obj.SX3_ZDown[0]) {	  
	    MyFill(Form("SX3Zpos_%i_%i_%i",Si.det_obj.DetID,Si.det_obj.DownChNum[0],Si.det_obj.BackChNum[0]),100,-1,1,Si.det_obj.SX3_ZDown[0]);
	  }else{
	    cout<<"1:  it shouldn't happen "<<endl;
	    break;
	  } 
	}	
#endif

	/////////////////// ZPosCal from the Processed data from the Hit  ///////////////////

	//if(Si.det_obj.HitType ==111) {
	if(Si.hit_obj.ZUp_Dummy <= 1.0 && Si.hit_obj.ZUp_Dummy >= -1.0 ) {
	  MyFill(Form("SX3ZposCal_%i_%i_%i",Si.hit_obj.DetID,Si.hit_obj.FrontChannel,Si.hit_obj.BackChannel),600,-1,1,Si.hit_obj.ZUp_Dummy);
	  MyFill(Form("SX3ZposCal_%i_f%i",Si.hit_obj.DetID,Si.hit_obj.FrontChannel),600,-1,1,Si.hit_obj.ZUp_Dummy);
	  //MyFill(Form("SX3ZposCal_%i_b%i",Si.hit_obj.DetID,Si.hit_obj.BackChannel),600,-1,1,Si.hit_obj.ZUp_Dummy);
	  MyFill(Form("SX3ZposCal_%i",Si.hit_obj.DetID),600,-1,1,Si.hit_obj.ZUp_Dummy);
	}
	if(Si.hit_obj.ZDown_Dummy <= 1.0 && Si.hit_obj.ZDown_Dummy >= -1.0 ) {
	  MyFill(Form("SX3ZposCal_%i_%i_%i",Si.hit_obj.DetID,Si.hit_obj.FrontChannel,Si.hit_obj.BackChannel),600,-1,1,Si.hit_obj.ZDown_Dummy);
	  MyFill(Form("SX3ZposCal_%i_f%i",Si.hit_obj.DetID,Si.hit_obj.FrontChannel),600,-1,1,Si.hit_obj.ZDown_Dummy);
	  //MyFill(Form("SX3ZposCal_%i_b%i",Si.hit_obj.DetID,Si.hit_obj.BackChannel),600,-1,1,Si.hit_obj.ZDown_Dummy);//show all front for given back
	  MyFill(Form("SX3ZposCal_%i",Si.hit_obj.DetID),600,-1,1,Si.hit_obj.ZDown_Dummy);
	}
	//}
#endif 		
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef Hist_for_PC_Cal	
	for ( Int_t hits=0; hits<PC.NPCHits; hits++ ) {
	  MyFill("PCPhi_vs_SiPhi_SX3",500,0,8,Si.hit_obj.PhiW,500,0,8,PC.Hit[hits].PhiW);
	}
#endif
	
      }//end ((down||up)&&back)
      else if ( Si.det_obj.EUp_Cal.size()!=0 || Si.det_obj.EDown_Cal.size()!=0 || Si.det_obj.EBack_Cal.size()!=0 ) {
#ifdef Pulser_ReRun
	//if only either of Up, Down or Back is fired in SX3, continue, unless it is a pulser check
	Si.det_obj.UpMult = Si.det_obj.EUp_Cal.size();
	Si.det_obj.DownMult = Si.det_obj.EDown_Cal.size();
	Si.det_obj.BackMult = Si.det_obj.EBack_Cal.size();
	  
	Si.det_obj.DetID = i+4;
	Si.det_obj.HitType = Si.det_obj.BackMult*100 + Si.det_obj.UpMult*10 + Si.det_obj.DownMult;

	Si.Detector.push_back(Si.det_obj);
	Si.NSiHits++;
#else
	continue;
#endif
      }

      //============================================================ 
    }//end of for(int i=0; i<NumSX3; i++){
    ////////////////////////////////////// Push back Q3 Detector-members ////////////////////////////////////////////////////
    for (int i=0; i<NumQ3; i++){

      Si.ZeroSi_obj();

      for (int j=0; j<MaxQ3Ch; j++){

	if ( (Q3Energy_Cal[i][j] > 0)){

	  if (j<16){//Back Channels of Q3

	    Si.det_obj.BackChNum.push_back(j);

#if defined(FillTree_Esteps) || defined(Hist_for_Si_Cal)
	    Si.det_obj.EBack_Raw.push_back(Q3Energy[i][j]);
	    Si.det_obj.EBack_Pulser.push_back(Q3Energy_Pulser[i][j]);
	    Si.det_obj.EBack_Rel.push_back(Q3Energy_Rel[i][j]);
#endif
	    Si.det_obj.EBack_Cal.push_back(Q3Energy_Cal[i][j]);
	    Si.det_obj.TBack.push_back(Q3Time[i][j]);

	  }

	  else if(j>15){//Front Channels of Q3 ,

	    Si.det_obj.FrontChNum.push_back(j-16); 

#if defined(FillTree_Esteps) || defined(Hist_for_Si_Cal)
	    Si.det_obj.EFront_Raw.push_back(Q3Energy[i][j]);
	    Si.det_obj.EFront_Pulser.push_back(Q3Energy_Pulser[i][j]);
	    Si.det_obj.EFront_Rel.push_back(Q3Energy_Rel[i][j]);
#endif
	    Si.det_obj.EFront_Cal.push_back(Q3Energy_Cal[i][j]);
	    Si.det_obj.TFront.push_back(Q3Time[i][j]);
	  }
	}
	//===========================================
      }//end of for(int j=0; j<MaxQ3Ch; j++){
      //===========================================
      Si.det_obj.FrontMult = Si.det_obj.EFront_Cal.size();
      Si.det_obj.BackMult = Si.det_obj.EBack_Cal.size();


      Si.det_obj.DetID = i;
      Si.det_obj.HitType = Si.det_obj.BackMult*10 + Si.det_obj.FrontMult;

      if ( Si.det_obj.EFront_Cal.size()!=0 && Si.det_obj.EBack_Cal.size()!=0 ){	

	SiSort.SortQ3(&Si,CMAP);
	Si.Detector.push_back(Si.det_obj);
	
	/////////////////////////////////////////  Fill Q3 Histograms after Calibration  ///////////////////////////////////////
#ifdef Hist_after_Cal
	if(Si.hit_obj.HitType ==11){
	  MyFill(Form("back_vs_front_Cal%i",Si.hit_obj.DetID),500,0,30,Si.hit_obj.EnergyFront,500,0,30,Si.hit_obj.EnergyBack);
	  MyFill(Form("Q3_offset_back_vs_front_Cal%i",Si.hit_obj.DetID),960,0,30,Si.hit_obj.EnergyBack,340,-10,10,(Si.hit_obj.EnergyBack-Si.hit_obj.EnergyFront));
	  //MyFill(Form("back_vs_front_Cal%i_f%i",Si.hit_obj.DetID,Si.hit_obj.FrontChannel),100,0,30,Si.hit_obj.EnergyFront,100,0,30,Si.hit_obj.EnergyBack);
	  //MyFill(Form("back_vs_front_Cal%i_b%i",Si.hit_obj.DetID,Si.hit_obj.BackChannel),100,0,30,Si.hit_obj.EnergyFront,100,0,30,Si.hit_obj.EnergyBack);
	  //MyFill(Form("back_vs_front_Cal%i_%i_%i",Si.hit_obj.DetID,Si.hit_obj.FrontChannel,Si.hit_obj.BackChannel),100,0,30,Si.hit_obj.EnergyFront,100,0,30,Si.hit_obj.EnergyBack);
	}
#endif
	////////////////////////////////////////  Fill Q3 Histograms for Calibration  ///////////////////////////////////////////
#ifdef Hist_for_Si_Cal
	Int_t fbmax=4*4096/3;
	Int_t bins=512;
	if(Si.det_obj.HitType ==11){//Just to make it simple.
	  //Step 1 RelCal//F-B	 //No need to give it a name Q3_ as we have DetID but since we are Calibration Differently for Q3 && SX3. 
	  MyFill(Form("Q3_back_vs_front%i_%i_%i",Si.det_obj.DetID,Si.det_obj.FrontChNum[0],Si.det_obj.BackChNum[0]),bins,0,fbmax,Si.det_obj.EFront_Rel[0],bins,0,fbmax,Si.det_obj.EBack_Rel[0]);
	  
	  //Step 2 RelCal//F-B //RelGain Cal from Step 1 is applied
	  MyFill(Form("Q3_back_vs_front%i_b%i",Si.det_obj.DetID,Si.det_obj.BackChNum[0]),bins,0,fbmax,Si.det_obj.EFront_Rel[0],bins,0,fbmax,Si.det_obj.EBack_Rel[0]);
	  //Just for check
	  MyFill(Form("Q3_back_vs_front%i",Si.det_obj.DetID),bins,0,fbmax,Si.det_obj.EFront_Rel[0],bins,0,fbmax,Si.det_obj.EBack_Rel[0]);

	  //check offset
  	  MyFill(Form("Q3_back_vs_offset%i",Si.det_obj.DetID),
		 200,-400,400,(Si.det_obj.EBack_Rel[0]-Si.det_obj.EFront_Rel[0]),
		 bins,0,fbmax,Si.det_obj.EBack_Rel[0]);
	}
#endif
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef Hist_for_PC_Cal
	for ( Int_t hits=0; hits<PC.NPCHits; hits++ ) {
	  MyFill("PCPhi_vs_SiPhi_Q3",500,0,8,Si.hit_obj.PhiW,500,0,8,PC.Hit[hits].PhiW);
	}
#endif	
	
      }else if ( Si.det_obj.EFront_Cal.size()!=0 || Si.det_obj.EBack_Cal.size()!=0 ){
	//if only either of Front or Back is fired in Q3, continue, unless it is a pulser check
#ifdef Pulser_ReRun
	Si.Detector.push_back(Si.det_obj);
	Si.NSiHits++;
#else
	continue;
#endif
      }
      //============================================================ 
    }//end of for(int i=0; i<NumQ3; i++){
    //============================================================ 
#ifdef IC_hists
    if(IC_E > 0) {
      MainTree->Fill();
    }
#else
    if (Si.NSiHits > 0) {
    //if (PC.NPCHits > 0) {    
      MainTree->Fill();
    }
#endif
    //============================================================
    //ncount++;
  }//end of for (global_evt=0; global_evt<nentries; global_evt++){
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  cout << endl << " Changing to output file... ";
  outputFile->cd();
  cout << filename_histout  << endl;
  cout << " Writing objects... ";
  RootObjects->Write();
  cout << "RootObjects are Written" << endl;
  outputFile->Close();
  cout << " Outputfile Closed\n";
  for(int i=0;i<1;i++) {//print beeps at end of program
    printf(" beep!\a\n");
    sleep(1);
  }
  return 0;
}
//end of Main()
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MyFill(string name, int binsX, double lowX, double highX, double valueX){
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
//===========================================================================
void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY){
  try{
    fhmap.at(name)->Fill(valueX,valueY);
  } catch(out_of_range e) {
    TH2F* newHist = new TH2F(name.c_str(),name.c_str(),
			     binsX,lowX,highX,
			     binsY,lowY,highY);
    newHist->Fill(valueX,valueY);
    fhlist->Add(newHist);
    fhmap[name] = newHist;  
  }  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
