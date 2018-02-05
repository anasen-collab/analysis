/////////////////////////////////////////////////////////////////////////////////////
// Class ChannelMap describes ANASEN channel map.
// Calculates X,Y,Z,Phi,...etc. for Q3 && SX3 in the forward annular & 2 barrel configurations detectors.
// Loads all the calibration coefficients from files when called in Main.
// Calculation solely depends on the files loaded.
//
// Usage: Include in the Main.cpp
//
// Authors: Nabin Rijal, John Parker, Ingo Wiedenhover, DSG, JAB, KTM, JCB.
//
// Last modified //2016 August
/////////////////////////////////////////////////////////////////////////////////////

#ifndef ChannelMap_h
#define ChannelMap_h

#include <TROOT.h>
#include <TMath.h>
#include <TRandom3.h>

#include <iostream>
#include <fstream>
#include <string>

#define NumDet 28
#define MaxChNum 500
#define NumSX3 24
#define NumQ3 4
#define MaxSX3Ch 12
#define MaxQ3Ch 32

#define WireNum 24
#define MaxADCChs 200
#define MaxADC 5
#define MaxADCCh 32

#define doprint kFALSE

using namespace std;
///////////////////////////////////////////////////////////////////////////////////
class ChannelMap {
  
  Int_t MBID[MaxChNum],CID[MaxChNum],ASICs_Ch[MaxChNum];
  Int_t Detector[MaxChNum],Det_Ch[MaxChNum];
  
  Int_t MBID_Align[MaxChNum],CID_Align[MaxChNum],ASICs_Ch_Align[MaxChNum];
  Int_t TotalNumberOfChannels,TotalAlignedASICsChannels,NumberOfSX3AlphaCalibrated;
  Int_t NumberOfQ3AlphaCalibrated,NumberOfQ3RelativeSlopes,NumberOfAlphaCalibrated;
  Int_t NumberOfSX3RelativeSlopes,NumberOfADCChannels;
  
  Int_t ADC[MaxADCChs];
  Int_t Channel[MaxADCChs];
  Int_t DetTypeID[MaxADCChs];
  Int_t Parameter1[MaxADCChs];
  Int_t Parameter2[MaxADCChs];
  Int_t BackChNum;
  
  Double_t PCSlope[WireNum], PCShift[WireNum];
  Double_t PCWire_RelGain[WireNum];
  
  string Comment[MaxChNum];
  Double_t SiGains[NumDet];
  Double_t SiOffsets[NumDet];
  Double_t SX3RelativeSlope[NumSX3][MaxSX3Ch],SX3FinalFix[NumSX3][MaxSX3Ch];
  Double_t Q3RelativeSlope[NumQ3][MaxQ3Ch];
  Double_t Q3FinalFix[NumQ3][MaxQ3Ch];
  Double_t EdgeU, EdgeD;
  Double_t zerosh[MaxChNum],vperch[MaxChNum];

  //added by M.Anastasiou 10/20/2016
  Double_t a[MaxChNum],b[MaxChNum],c[MaxChNum],q0[MaxChNum];

  Double_t ZOffset[NumSX3],XAt0[NumSX3],XAt4[NumSX3],YAt0[NumSX3],YAt4[NumSX3];

  Double_t PCPulser_YOffset[MaxADC][MaxADCCh];
  Double_t PCPulser_Slope[MaxADC][MaxADCCh];

  //------------------------PC Relative Gains---------------------added 05/05/2017------------------//
  Double_t PC_UD_Slope[WireNum];
  Double_t PC_UD_Offset[WireNum];

  TRandom3 *Randomm;
  
 public:
  
  Double_t EdgeUp[NumSX3][4][MaxSX3Ch], EdgeDown[NumSX3][4][MaxSX3Ch];

  //////////////////// Constructor /////////////////////////////////////
  ChannelMap() {
    TotalNumberOfChannels=0;
    TotalAlignedASICsChannels=0;
    NumberOfSX3AlphaCalibrated=0;
    NumberOfQ3AlphaCalibrated=0;
    NumberOfAlphaCalibrated=0;
    NumberOfSX3RelativeSlopes=0;
    NumberOfQ3RelativeSlopes=0;
    NumberOfADCChannels=0;

    Randomm = new TRandom3();
  };
  //////////////// Destructor //////////////////////////////////////////
  ~ChannelMap() {
    
    delete Randomm;
  };
  ///////////  Load && Initializes Maps and Calibration Files //////////////////////

  int LoadASICsChannelMapFile(const char* ASICsChannelMapFilename);
  int LoadASICsPulserAlignment(const char* ASICsPulserFilename);

  //added by M.Anastasiou 10/20/2016
  //int LoadASICsPulserAlignment_Quadratic (const char* ASICsPulserFilename);

  int LoadSiGains(const char* SiGainsFilename);
  int LoadSX3RelativeSlopes(const char* SX3SlopeFilename);
  int LoadQ3RelativeSlopes(const char* Q3SlopeFilename);  
 
  int Init(const char* ASICsChannelMapFilename, const char* ASICsPulserFilename, const char* SiGainsFilename, 
	   const char* SX3RelativeSlopeFilename, const char* Q3RelativeSlopeFilename);
  int FinalInit(const char* FinalAdjustments, const char* SX3Geo);
  int LoadQ3FinalFix(const char* Q3FinalFixFilename);
  int InitWorldCoordinates(const char* WorldCoordinatesFilename);   

  int InitPCADC(const char* PCMapFilename);
  int InitPCCalibration(const char* PCCalibrationFilename);
  int InitPCWireCal(const char* PCWireCalFilename);
  int Init_PCWire_RelGain(const char* PCWire_RelGain_Filename);

  int LoadFail(const char*);

  /////////////////// Inline Functions //////////////////////////////////////////////

  int GetNumSX3Det() {return NumberOfSX3AlphaCalibrated;};
  int GetNumQ3Det() {return NumberOfQ3AlphaCalibrated;};
  int GetNumDet() {return NumberOfAlphaCalibrated;}
  int GetNumberOfADCChannels() {return NumberOfADCChannels;};
  int GetNumberOfMappedASICsChannels() {return TotalNumberOfChannels;};
  int GetNumberOfASICsAligned() {return TotalAlignedASICsChannels;};  
  int GetNumberOfAlphaCalibrated() {return NumberOfAlphaCalibrated;};
  int GetNumberOfSX3RelativeSlopes() {return NumberOfSX3RelativeSlopes;};
  int GetNumberOfQ3RelativeSlopes() {return NumberOfQ3RelativeSlopes;};

  /////////////////////// Use Coefficients && Calculations //////////////////////////
  
  void IdentifyADC(Int_t ADCid, Int_t CHid, Int_t& DetType);
  void IdentifyWire(Int_t ADCid, Int_t CHid, Int_t& WireID, Int_t& Side);
  Bool_t ConvertToVoltage(Int_t ADCid, Int_t CHid, Int_t PCData, Double_t& Vcal);
  void Get_PCWire_RelGain(Int_t WireID,Double_t& PCRelGain); 
  void GetPCWorldCoordinates(Int_t wireid, Double_t zpos, Double_t& xw, Double_t& yw, Double_t& zw, Double_t& rw, Double_t& phiw);
 
  
  void IdentifyDetChan(Int_t mb_id, Int_t chip_id, Int_t asic_ch, Int_t& det, Int_t& det_ch);  
  void AlignASICsChannels(Int_t mb_id_al, Int_t chip_id_al, Int_t asic_ch_al, Double_t& zero, Double_t& gain);  

  //added by M.Anastasiou 10/20/2016 
  //void AlignASICsChannels_Quadratic(Int_t mb_id_al, Int_t chip_id_al, Int_t asic_ch_al, Double_t& alpha, Double_t& beta, Double_t& gamma, Double_t& Z_shift);  

  
  void GetSX3MeVPerChannel1(Int_t DN, Int_t DetCh, Double_t& slope);
  void GetSX3MeVPerChannel2(Int_t DN, Int_t DetCh, Double_t& slope);
  void GetSX3FinalEnergyOffsetInMeV(Int_t DNum, Int_t ChNum, Double_t& zshift);
  void GetQ3MeVPerChannel1(Int_t DN, Int_t DetCh, Double_t& slope);
  void GetQ3MeVPerChannel2(Int_t DN, Int_t DetCh, Double_t& slope);  
  void GetQ3FinalEnergyOffsetInMeV(Int_t DNum, Int_t ChNum, Double_t& zshift);

  void GetQ3WorldCoordinates(Int_t DID, Double_t SiX, Double_t SiY, Double_t& WSiX, Double_t& WSiY, Double_t& WSiR, Double_t& WSiPhi );
  void GetSX3WorldCoordinates(Int_t DID, Double_t SiX, Double_t SiZ, Double_t& WSiX, Double_t& WSiY, Double_t& WSiZ, Double_t& WSiR, Double_t& WSiPhi);
  void PosCal(Int_t DNum, Int_t StripNum, Int_t ChNum, Double_t FinalZPos, Double_t& FinalZPosCal);
  
  void GetZeroShift(Int_t det, Int_t det_ch, Double_t& zero, Double_t& slope);
  Double_t GetRelLinCoeff(Int_t det, Int_t det_ch) { return SX3RelativeSlope[det-4][det_ch]; };  
  void IdentifyMbChipChan(Int_t det, Int_t det_ch,Int_t &mb_id,Int_t &chip_id,Int_t &asic_ch);

  //------------------------PC Relative Gains---------------------added 05/05/2017------------------//
  int Init_PC_UD_RelCal(const char* PC_UD_RelCal_Filename);
  void Get_PC_UD_RelCal(Int_t WireID, Double_t& SlopeUD, Double_t& OffsetUD);
};

////////////////////////////////////////////////////////////////////////////////////////////////////
int ChannelMap::LoadFail(const char* filename) {
  cout << "Cannot open file " << filename <<endl;
  cout << "Macro terminated abnormally!" << endl;
  exit(EXIT_FAILURE);
  return 0;
}

//------------------------------------------------------------------------------------------------//
int ChannelMap::LoadASICsChannelMapFile (const char* ASICsChannelMapFilename) {
  ifstream channelmapfile;
  string line;
  
  channelmapfile.open(ASICsChannelMapFilename);
  
  cout << ASICsChannelMapFilename << endl;
	
  if (channelmapfile.is_open()) {
    cout << "The channel map file " << ASICsChannelMapFilename << " opened successfully." << endl;
    getline (channelmapfile,line);//Skips the first line in ASICsChannelMapFilename.
    if(doprint) cout<<"line = "<<line<<endl;
    Int_t i=0;

    while (!channelmapfile.eof()) {
      channelmapfile >> MBID[i] >> CID[i] >> ASICs_Ch[i] >> Detector[i] >> Det_Ch[i] >> Comment[i];
      //if(doprint) cout << MBID[i] << "   " << CID[i] << "   " << endl;
      i++;
    }
    TotalNumberOfChannels = i-1;
  }
  else LoadFail(ASICsChannelMapFilename);
  channelmapfile.close();
  return 1;
}

//------------------------------------------------------------------------------------------------//
int ChannelMap::LoadASICsPulserAlignment (const char* ASICsPulserFilename) {
  ifstream alignchan;
  string line;

  alignchan.open(ASICsPulserFilename);
  if (alignchan.is_open()) {
    cout << "The channel allignment file " << ASICsPulserFilename << " opened successfully." << endl;
    getline (alignchan,line);//Skips the first line in ASICsPulserFilename.
    if(doprint) cout<<"line = "<<line<<endl;
    Int_t j=0;
    while (!alignchan.eof()) {
      alignchan >>  MBID_Align[j] >> CID_Align[j] >> ASICs_Ch_Align[j] >> zerosh[j] >> vperch[j];
      j++;
    }
    TotalAlignedASICsChannels = j-1;
  }
  else LoadFail(ASICsPulserFilename);
  return 1;
}

//------------------------------------------------------------------------------------------------//
/*int ChannelMap::LoadASICsPulserAlignment_Quadratic (const char* ASICsPulserFilename) {
  ifstream alignchan;
  string line;

  alignchan.open(ASICsPulserFilename);
  if (alignchan.is_open()) {
    cout << "The channel allignment file " << ASICsPulserFilename << " opened successfully." << endl;
    //getline (alignchan,line);//Skips the first line in ASICsPulserFilename.
    if(doprint) cout<<"line = "<<line<<endl;
    /*    Int_t j=0;
	  while (!alignchan.eof()) {
	  alignchan >>  MBID_Align[j] >> CID_Align[j] >> ASICs_Ch_Align[j] >> a[j] >> b[j] >> c[j] >> q0[j];
	  j++;
	  }
	  TotalAlignedASICsChannels = j-1;
	  }
	  else LoadFail(TotalAlignedASICsChannels);
  	  return 1;
	  }*/

//------------------------------------------------------------------------------------------------//
int ChannelMap::LoadSiGains(const char* SiGainsFilename) {
  ///////////////////////////////////////////////////////////
  // Read the energyslopes file and store data in the array
  // ESlope[DN].
  // NOTE: In the energy slope file there is unnecessary
  //       information in column 2.  Therefore we read it in
  //       as 'Blank' and forget about it.
  ///////////////////////////////////////////////////////////
  
  ifstream SiGainsFile;
  string line;
  Int_t i,j;

  SiGainsFile.open(SiGainsFilename);
  if (SiGainsFile.is_open()) {
    cout << "The Si Gains file " << SiGainsFilename << " opened successfully." << endl;

    getline (SiGainsFile,line);//Skips the first line in SiGainsFilename.
    if(doprint) cout<<"line = "<<line<<endl;
    i=0; 
    j=0;
    Int_t Dnum = 0;
    Double_t ChNum = 0.0;
    Double_t dummy = 0.0;
    
    while (!SiGainsFile.eof()) {
      if(Dnum>=4 && Dnum<=27) {i++;}
      else if(Dnum<=3) {j++;}
      else {
      	cout << Dnum << "not read" <<endl;
	continue;
	}
      SiGainsFile >>  Dnum >> ChNum >> dummy;
      SiGains[Dnum] = dummy;
      SiOffsets[Dnum]= ChNum;
      //printf("  SiGains[%d] = %f\n",Dnum,SiGains[Dnum]);
    }
    
    NumberOfSX3AlphaCalibrated = i-1; // The Minus one accounts for the end of line character
    NumberOfQ3AlphaCalibrated = j;
    NumberOfAlphaCalibrated = NumberOfSX3AlphaCalibrated + NumberOfQ3AlphaCalibrated;
    if(doprint)printf(" %d + %d = %d calibrated\n",NumberOfSX3AlphaCalibrated,NumberOfQ3AlphaCalibrated,NumberOfAlphaCalibrated);
  }
  else LoadFail(SiGainsFilename);
  SiGainsFile.close();
  return 1;
}

//------------------------------------------------------------------------------------------------//
int ChannelMap::LoadSX3RelativeSlopes(const char* SX3RelativeSlopeFilename) {
  ifstream SX3RelativeSlopeFile;
  string line;

  SX3RelativeSlopeFile.open(SX3RelativeSlopeFilename);  
  
  if (SX3RelativeSlopeFile.is_open()) {
    cout << "The SX3 slopes file " << SX3RelativeSlopeFilename << " opened successfully." << endl;

    getline (SX3RelativeSlopeFile,line);//Skips the first line in SX3RelativeSlopeFilename
    if(doprint) cout<<"line = "<<line<<endl;
    Int_t m=0;
    Int_t DnumSX3 = 0;
    Int_t SX3ChNum = 0;
    Double_t dummy = 0.0;
    
    while (!SX3RelativeSlopeFile.eof()) {
      SX3RelativeSlopeFile >>  DnumSX3 >> SX3ChNum >> dummy;
      SX3RelativeSlope[DnumSX3-4][SX3ChNum] = dummy;
      m++;
    }
    NumberOfSX3RelativeSlopes = m-1;
    
  }
  else LoadFail(SX3RelativeSlopeFilename);
  SX3RelativeSlopeFile.close();
  return 1;
}

//------------------------------------------------------------------------------------------------//
int ChannelMap::LoadQ3RelativeSlopes(const char* Q3RelativeSlopeFilename) {
  ifstream Q3RelativeSlopeFile;
  string line;

  Q3RelativeSlopeFile.open(Q3RelativeSlopeFilename);  
  
  if (Q3RelativeSlopeFile.is_open()) {
    cout << "The Q3 slopes file " << Q3RelativeSlopeFilename << " opened successfully." << endl;

    getline(Q3RelativeSlopeFile,line);//Skips the first line in Q3RelativeSlopeFilename
    if(doprint) cout<<"line = "<<line<<endl;
    Int_t m=0;
    Int_t DnumQ3 = 0;
    Int_t Q3ChNum = 0;
    Double_t dummy = 0.0;
    
    while (!Q3RelativeSlopeFile.eof()) {
      Q3RelativeSlopeFile >>  DnumQ3 >> Q3ChNum >> dummy;
      Q3RelativeSlope[DnumQ3][Q3ChNum] = dummy;
      m++;
    }
    NumberOfQ3RelativeSlopes = m-1;
    
  }
  else LoadFail(Q3RelativeSlopeFilename);
  Q3RelativeSlopeFile.close();
  return 1;
}

//------------------------------------------------------------------------------------------------//
int ChannelMap::Init(const char* ASICsChannelMapFilename, 
		     const char* ASICsPulserFilename,
                     const char* SiGainsFilename, 
		     const char* SX3RelativeSlopeFilename, 
		     const char* Q3RelativeSlopeFilename) {
  Int_t status = 0;
  status = LoadASICsChannelMapFile(ASICsChannelMapFilename);
  if(status == 1) {
    status = LoadASICsPulserAlignment(ASICsPulserFilename);
    //status = LoadASICsPulserAlignment_Quadratic(ASICsPulserFilename);
  }
  else LoadFail(ASICsChannelMapFilename);
  
  if(status == 1) {
    status = LoadSiGains(SiGainsFilename);
  }
  else LoadFail(ASICsPulserFilename);
  
  if(status == 1) {
    status = LoadSX3RelativeSlopes(SX3RelativeSlopeFilename);
  }
  else LoadFail(SiGainsFilename);
  
  if(status == 1) {
    status = LoadQ3RelativeSlopes(Q3RelativeSlopeFilename);
  }
  else LoadFail(SX3RelativeSlopeFilename);
  
  if(status == 0) LoadFail(Q3RelativeSlopeFilename);
  
  return status;
}

//------------------------------------------------------------------------------------------------//
// Description: Fine tuning parameters for better energy calibration
// The file contains detector number, channel number and zero shift in MeV
// By default, the zero shift is zero. Only those channels that have non-zero
// shift need to be described in the final fix file.

int ChannelMap::FinalInit(const char* FinalFixFilename, const char* SX3GeoFilename) {
  ifstream finalfix;
  ifstream x3geo;
  string line1, line2;

  finalfix.open(FinalFixFilename);
  x3geo.open(SX3GeoFilename);

  Int_t DNum,ChNum,StripNum;
  Double_t Zero;
  
  for (Int_t i=0; i<NumSX3; i++) {
    for (Int_t c=0; c<MaxSX3Ch; c++) {
      SX3FinalFix[i][c] = 0;
    }
  }
  //---------------------------------------------------------------
  if (finalfix.is_open()) {
    cout << "File with final fix zero shifts " << FinalFixFilename;
    cout << " opened successfully." << endl;

    getline(finalfix,line1);//Skips the first line in FinalFixFilename
    if(doprint) cout<<"line = "<<line1<<endl;

    while (!finalfix.eof()) {
      finalfix >> DNum >> ChNum >> Zero;
      SX3FinalFix[DNum-4][ChNum] = Zero;
    }

  }
  else LoadFail(FinalFixFilename);
  //---------------------------------------------------------------  
  
  if(x3geo.is_open()) {
    cout << "File with SX3 Geometry " << SX3GeoFilename;
    cout << " opened successfully." << endl;

    getline(x3geo,line2);//Skips the first line in SX3GeoFilename
    if(doprint) cout<<"line = "<<line2<<endl;

    while (!x3geo.eof()) {
      x3geo >> DNum >> StripNum >> BackChNum >> EdgeD >> EdgeU;
      EdgeDown[DNum-4][StripNum][BackChNum] = EdgeD;
      EdgeUp[DNum-4][StripNum][BackChNum] = EdgeU;
    }
  }
  else LoadFail(SX3GeoFilename);
  //---------------------------------------------------------------
  finalfix.close();
  x3geo.close();
  return 1;
}

//------------------------------------------------------------------------------------------------//
int ChannelMap::LoadQ3FinalFix(const char* Q3FinalFixFilename) {
  ifstream Q3finalfix;
  string line;
  Q3finalfix.open(Q3FinalFixFilename);
  Double_t Zero = 0;
  Int_t DNum,ChNum;
  
  for (Int_t i=0; i<4; i++) {
    for (Int_t c=0; c<32; c++) {
      Q3FinalFix[i][c] = 0;
    }
  }
  //------------------------------------
  if (Q3finalfix.is_open()) {
    cout << "File with Q3 final fix zero shifts " << Q3FinalFixFilename;
    cout << " opened successfully." << endl;

    getline(Q3finalfix,line);//Skips the first line in Q3FinalFixFilename
    if(doprint) cout<<"line = "<<line<<endl;

    while (!Q3finalfix.eof()) {
      Q3finalfix >> DNum >> ChNum >> Zero;
      Q3FinalFix[DNum][ChNum] = Zero;
    }
  }
  else LoadFail(Q3FinalFixFilename);
  return 1;  
}

//------------------------------------------------------------------------------------------------//
int ChannelMap::InitWorldCoordinates(const char* WorldCoordinatesFilename) {
  ifstream worldCfile;
  worldCfile.open(WorldCoordinatesFilename);
  
  Int_t Dnum;
  Double_t Zoff,XMin,XMax,YMin,YMax;
  string DumComment,line;
  
  // The ZOffset is given with respect to the face of the forward Si detectors.
  // Note that if initial ZOffset is negative then the detector is considered as missing.
  // If YAt0 and YAt4 are both more than +50.0 cm then this is the forward detector.
  
  for (Int_t i = 0; i<NumSX3; i++) {
    ZOffset[i]  = -100.0;
    XAt0[i]     = 0;
    XAt4[i]     = 0;
    YAt0[i]     = 0;
    YAt4[i]     = 0;
  }
  
  if (worldCfile.is_open()) {
    cout << "File with world coordinates " << WorldCoordinatesFilename;
    cout << " opened successfully." << endl;
    getline (worldCfile,line);//Skips the first line in WorldCoordinatesFilename.
    if(doprint) cout<<"line = "<<line<<endl;
   
    while (!worldCfile.eof()) {
   // worldCfile >> Dnum >> Zoff >> XMin >> XMax >> YMin >> YMax >> DumComment;
      worldCfile >> Dnum >> Zoff >> XMin >> XMax >> YMin >> YMax ;

      ZOffset[Dnum-4]  = Zoff;
      XAt0[Dnum-4]     = XMin;
      XAt4[Dnum-4]     = XMax;
      YAt0[Dnum-4]     = YMin;
      YAt4[Dnum-4]     = YMax;
    }

    if(doprint)
      for (Int_t i = 0; i<NumSX3; i++) {
	cout << "\t" << i + 4 << "\t" << ZOffset[i] << "\t " <<  XAt0[i] << "\t " <<  XAt4[i] << "\t " << YAt0[i] << "\t " << YAt4[i] << endl;
      }
    
  }
  else LoadFail(WorldCoordinatesFilename);
  worldCfile.close();
  return 1;
}

//------------------------------------------------------------------------------------------------//
int ChannelMap::InitPCADC(const char* PCMapFilename) {
  ifstream pcmap;
  string line;
  Int_t DumADCID,DumCh,DumWireID,DumSide;
  
  pcmap.open(PCMapFilename);
  
  if (pcmap.is_open()) {
    cout << "Channel map for proportional counter " << PCMapFilename << " opened sucessfully. " << endl;
    getline (pcmap,line);//Skips the first line in PCMapFilename.
    if(doprint) cout<<"line = "<<line<<endl;
    while (!pcmap.eof()) {
      pcmap >> DumADCID >> DumCh >> DumWireID >> DumSide ;

      ADC[NumberOfADCChannels] = DumADCID;
      Channel[NumberOfADCChannels] = DumCh;
      DetTypeID[NumberOfADCChannels] = 1;
 
      Parameter1[NumberOfADCChannels] = DumWireID;
      Parameter2[NumberOfADCChannels] = DumSide;
      NumberOfADCChannels++;
    }
    NumberOfADCChannels--;
  }
  else LoadFail(PCMapFilename);
  pcmap.close();
  return 1;
}

//------------------------------------------------------------------------------------------------//
int ChannelMap::InitPCCalibration(const char* PCCalibrationFilename) {
  ifstream pccal;
  string line;
  Int_t adcid,chnum;
  Double_t dum_Slope,dum_Yoff;
  
  pccal.open(PCCalibrationFilename);
  
  for (Int_t k=0; k<MaxADC; k++) {
    for (Int_t i=0; i<MaxADCCh; i++) {
      PCPulser_YOffset[k][i] = 0;
      PCPulser_Slope[k][i] = 0;
    }
  }
  
  if (pccal.is_open()) {
    cout << "Pulser Calibration for proportional counter " << PCCalibrationFilename << " opened successfully. " << endl;
    getline (pccal,line);//Skips the first line in PCCalibrationFilename
    if(doprint) cout<<"line = "<<line<<endl;

    while (!pccal.eof()) {
      pccal >> adcid >> chnum >> dum_Yoff >> dum_Slope;
     
      PCPulser_YOffset[adcid][chnum] = dum_Yoff;
      PCPulser_Slope[adcid][chnum] = dum_Slope;
    }
  }
  else LoadFail(PCCalibrationFilename);

  /* for (Int_t k=2; k<4; k++) {
    for (Int_t i=0; i<MaxADCCh; i++) {
      if(k==3 && i>15) continue;
      cout << k << "\t "<< i << "\t"<< PCPulser_YOffset[k][i] << "\t"<< PCPulser_Slope[k][i] << endl;
    }
    }*/
  
  return 1;
}

//------------------------------------------------------------------------------------------------//
int ChannelMap::InitPCWireCal(const char* PCWireCalFilename) {
  ifstream pcwirecal;
  string line;
  Double_t pcslopedum, pcshiftdum;
  Int_t WireNumber=0;
  
  pcwirecal.open(PCWireCalFilename);
  
  for (Int_t i=0; i<WireNum; i++) {
    PCSlope[i]=0.0;
    PCShift[i]=0.0;
  }
  
  if (pcwirecal.is_open()) {
    cout << "Alpha Calibration for proportional counter wires ";
    cout << PCWireCalFilename << " opened successfully. " << endl;
    getline (pcwirecal,line);//Skips the first line in PCWireCalFilename
    if(doprint) cout<<"line = "<<line<<endl;

    while (!pcwirecal.eof()) {
      pcwirecal >> WireNumber >> pcslopedum >> pcshiftdum;
      PCSlope[WireNumber] = pcslopedum;
      PCShift[WireNumber] = pcshiftdum;
    }
  }
  else LoadFail(PCWireCalFilename);
  return 1;
}

//------------------------------------------------------------------------------------------------//
int ChannelMap::Init_PCWire_RelGain(const char* PCWire_RelGain_Filename) {
  ifstream pcwire_rel;
  string line1;
  Double_t rel_dummy;
  Int_t WireID=0;

  pcwire_rel.open(PCWire_RelGain_Filename);

  for (Int_t i=0; i<WireNum; i++) {   
    PCWire_RelGain[i]=1;
  }

  if (pcwire_rel.is_open()) {
    cout << "Relative Calibration for PC Wires ";
    cout << PCWire_RelGain_Filename << " opened successfully. " << endl;

    getline (pcwire_rel,line1);//Skips the first line in PCWire_RelGain_Filename.
    if(doprint) cout<<"line = "<<line1<<endl;

    while (!pcwire_rel.eof()) {
      pcwire_rel >> WireID >> rel_dummy;
      PCWire_RelGain[WireID] = rel_dummy;
      if(doprint) printf("Gain for wire %d is %f\n",WireID,PCWire_RelGain[WireID]);
    }
  }
  else LoadFail(PCWire_RelGain_Filename);
  return 1;
}

//------------------------------------------------------------------------------------------------//
void ChannelMap::IdentifyADC(Int_t ADCid, Int_t CHid, Int_t& DetType) {
  DetType = 0;
  for (Int_t i=0; i<NumberOfADCChannels; i++) {
    if (ADCid == ADC[i] && CHid == Channel[i])  DetType=DetTypeID[i];
  }  
}

//------------------------------------------------------------------------------------------------//
void ChannelMap::IdentifyWire(Int_t ADCid, Int_t CHid, Int_t& WireID, Int_t& Side) {
  for (Int_t i=0; i<NumberOfADCChannels; i++) {
    if (ADCid == ADC[i] && CHid == Channel[i]) {
      WireID = Parameter1[i];
      Side   = Parameter2[i];
    }
  }  
}

//------------------------------------------------------------------------------------------------//
Bool_t ChannelMap::ConvertToVoltage(Int_t ADCid, Int_t CHid, Int_t PCData, Double_t& Vcal) {  
  Vcal = sqrt(-1);
  if (PCData > 0) {
    Vcal = PCData*PCPulser_Slope[ADCid][CHid] + PCPulser_YOffset[ADCid][CHid];  
    return kTRUE;
  }
  else
    return kFALSE;
}

//-------------------------------------------------------------------------------------------------//
void ChannelMap::Get_PCWire_RelGain(Int_t WireID,Double_t& PCRelGain) {

  PCRelGain = PCWire_RelGain[WireID];

}

//------------------------------------------------------------------------------------------------//
void ChannelMap::GetPCWorldCoordinates(Int_t wireid, Double_t zpos, Double_t& xw, Double_t& yw, Double_t& zw, Double_t& rw, Double_t& phiw) {
  Randomm->SetSeed();
  Double_t Radius = 3.8463;
  Double_t Angle = (24-(Double_t)wireid+Randomm->Rndm()-0.5)*TMath::TwoPi()/24.0 + TMath::Pi()/2;
  //changed angle on Feb27 to correspond to reality
  xw = Radius*TMath::Cos(Angle);
  yw = Radius*TMath::Sin(Angle);
  zw = zpos*PCSlope[wireid]+PCShift[wireid];
  rw = TMath::Sqrt(xw*xw + yw*yw);
  /*
    if(yw > 0 && xw>0) {
    phiw = TMath::ATan(xw/yw);
    }else if(yw > 0 && xw<0) {
    phiw = TMath::ATan(-xw/yw) + 0.5*TMath::Pi();
    }else if(yw < 0 && xw<0) {
    phiw = TMath::ATan(xw/yw) + TMath::Pi();
    }else if(yw < 0 && xw>0) {
    phiw = TMath::ATan(-xw/yw) + 1.5*TMath::Pi();
    }
  */
  //phiw = Angle-TMath::Pi()/2;
  phiw = Angle;
  if (phiw > 2*TMath::Pi()){
    phiw = phiw - 2*TMath::Pi();
  }
}
//------------------------------------------------------------------------------------------------//
//Description: From the channel map provided in the Init method this function get returns the
//             detector number and detector channel number for a given mother board, chip and
//             chip channel.

void ChannelMap::IdentifyDetChan(Int_t mb_id, Int_t chip_id, Int_t asic_ch, Int_t& det, Int_t& det_ch) {
  for (Int_t i = 0; i<TotalNumberOfChannels; i++) {
    if (MBID[i] == mb_id && CID[i] == chip_id && ASICs_Ch[i] == asic_ch) {
      det = Detector[i];
      det_ch = Det_Ch[i];
      return;
    }
  }
}

//------------------------------------------------------------------------------------------------//
void ChannelMap::AlignASICsChannels(Int_t mb_id, Int_t chip_id, Int_t asic_ch, Double_t& zero, Double_t& gain) {
  for (Int_t i = 0; i<TotalAlignedASICsChannels; i++) {
    if (MBID_Align[i] == mb_id && CID_Align[i] == chip_id && ASICs_Ch_Align[i] == asic_ch) {
      zero = zerosh[i];
      gain = vperch[i];
    }    
  }
}

//------------------------------------------------------------------------------------------------//
/*void ChannelMap::AlignASICsChannels_Quadratic(Int_t mb_id, Int_t chip_id, Int_t asic_ch, Double_t& alpha, Double_t& beta, Double_t& gamma, Double_t& X_shift) {
  for (Int_t i = 0; i<TotalAlignedASICsChannels; i++)
  {
  if (MBID_Align[i] == mb_id && CID_Align[i] == chip_id && ASICs_Ch_Align[i] == asic_ch)
  {
  alpha = a[i];
  beta = b[i];
  gamma = c[i];
  X_shift = q0[i];
  }    
  }
  }*/

//------------------------------------------------------------------------------------------------//
void ChannelMap::GetSX3MeVPerChannel1(Int_t DN, Int_t DetCh, Double_t& slope) {
  slope = SX3RelativeSlope[DN-4][DetCh];
}

//------------------------------------------------------------------------------------------------//
void ChannelMap::GetSX3MeVPerChannel2(Int_t DN, Int_t DetCh, Double_t& slope) {
  slope = SiGains[DN];
}

//------------------------------------------------------------------------------------------------//
void ChannelMap::GetSX3FinalEnergyOffsetInMeV(Int_t DNum, Int_t ChNum, Double_t& zshift) {
  zshift = SX3FinalFix[DNum-4][ChNum];
}

//------------------------------------------------------------------------------------------------//
void ChannelMap::GetQ3MeVPerChannel1(Int_t DN, Int_t DetCh, Double_t& slope) {
  slope = Q3RelativeSlope[DN][DetCh];
}

//------------------------------------------------------------------------------------------------//
void ChannelMap::GetQ3MeVPerChannel2(Int_t DN, Int_t DetCh, Double_t& slope) {
  slope = SiGains[DN];
}

//------------------------------------------------------------------------------------------------//
void ChannelMap::GetQ3FinalEnergyOffsetInMeV(Int_t DNum, Int_t ChNum, Double_t& zshift) {
  zshift = Q3FinalFix[DNum][ChNum];
}

//------------------------------------------------------------------------------------------------//
void ChannelMap::GetQ3WorldCoordinates(Int_t DID, Double_t SiX, Double_t SiY, Double_t& WSiX, Double_t& WSiY, Double_t& WSiR, Double_t& WSiPhi) {
  //Double_t Theta = TMath::Pi()*(DID+1)/2; //Rotate detectors CCW by an angle Theta
 
  Double_t Theta = 0;
  if (DID==0) {
    Theta = 0;
  }else if (DID==1) {
    Theta = 3*TMath::Pi()/2;
  }else if (DID==2) {
    Theta = TMath::Pi();
  }else if (DID==3) {
    Theta = TMath::Pi()/2;
  }
 
  WSiX = SiX*TMath::Cos(Theta)-SiY*TMath::Sin(Theta);
  WSiY = SiX*TMath::Sin(Theta)+SiY*TMath::Cos(Theta);
  //WSiY = SiX*TMath::Cos(Theta)-SiY*TMath::Sin(Theta);
  //WSiX = SiX*TMath::Sin(Theta)+SiY*TMath::Cos(Theta);

  WSiR = TMath::Sqrt(WSiX*WSiX + WSiY*WSiY);

  //Phi calculated to go from 0 to 2*pi--0 is along positive x
  //in lab frame, x is North, z is upstream, and y is up

  if(WSiY >= 0 && WSiX>=0) {
    WSiPhi = TMath::ATan(WSiY/WSiX);
  }else if(WSiY >= 0 && WSiX<0) {
    WSiPhi = TMath::ATan(WSiY/WSiX) + TMath::Pi();
  }else if(WSiY < 0 && WSiX<0) {
    WSiPhi = TMath::ATan(WSiY/WSiX) + TMath::Pi();
  }else if(WSiY < 0 && WSiX>=0) {
    WSiPhi = TMath::ATan(WSiY/WSiX) + 2*TMath::Pi();
  }
  if(doprint) cout << WSiX << "  " << WSiY << "  " << WSiPhi << "  " << endl;
}

//------------------------------------------------------------------------------------------------//
void ChannelMap::GetSX3WorldCoordinates(Int_t DID, Double_t SiX, Double_t SiZ, Double_t& WSiX,
					Double_t& WSiY, Double_t& WSiZ, Double_t& WSiR, Double_t& WSiPhi) {
  if (ZOffset[DID-4]>-0.1 && YAt0[DID-4]<50. && YAt0[DID-4]<50.) {
   
    //WSiZ = ZOffset[DID-4] + 7.5 - SiZ;
    WSiZ = ZOffset[DID-4] + SiZ;
    if(doprint) cout << "Z: " << SiZ << "   Z Offset: " << ZOffset[DID-4] << "   Det ID: " << DID << endl;
    WSiX = (XAt4[DID-4] - XAt0[DID-4])*0.25*SiX + XAt0[DID-4];
    WSiY = (YAt4[DID-4] - YAt0[DID-4])*0.25*SiX + YAt0[DID-4];
    WSiR = TMath::Sqrt(WSiX*WSiX + WSiY*WSiY);
    if(WSiY >= 0 && WSiX>=0) {
      WSiPhi = TMath::ATan(WSiY/WSiX);
    }else if(WSiY >= 0 && WSiX<0) {
      WSiPhi = TMath::ATan(WSiY/WSiX) + TMath::Pi();
    }else if(WSiY < 0 && WSiX<0) {
      WSiPhi = TMath::ATan(WSiY/WSiX) + TMath::Pi();
    }else if(WSiY < 0 && WSiX>=0) {
      WSiPhi = TMath::ATan(WSiY/WSiX) + 2*TMath::Pi();
    }
  }
  else{
    //This detector does not have world coordinates description
    cout << " NO COORDINATES";
    cout << " Det ID: " << DID << "Z: " << SiZ << " Z Offset: " << ZOffset[DID-4] << endl;
    WSiX = 0;
    WSiY = 0;
    WSiZ = 0;
  }
}

//------------------------------------------------------------------------------------------------//
void ChannelMap::PosCal(Int_t DNum, Int_t StripNum, Int_t BChNum, Double_t FinalZPos, Double_t& FinalZPosCal) {
  Double_t EdgeDCal=sqrt(-1);
  Double_t EdgeUCal=sqrt(-1);
  
  switch(BChNum) {
    case 0:
      EdgeDCal = 7.5;
      EdgeUCal = 5.625;
      break;
    case 1:
      EdgeDCal = 5.625;
      EdgeUCal = 3.75;
      break;
    case 2:
      EdgeDCal = 3.75;
      EdgeUCal = 1.875;
      break;
    case 3:
      EdgeDCal = 1.875;
      EdgeUCal = 0.00;
      break;
    default:
      break;
    }
  
  FinalZPosCal = (EdgeDCal-EdgeUCal)/(EdgeDown[DNum-4][StripNum][BChNum]-EdgeUp[DNum-4][StripNum][BChNum])
    *(FinalZPos-EdgeDown[DNum-4][StripNum][BChNum])+EdgeDCal;

  if(doprint)
    if(DNum>15) {
      printf("Det num = %d Z = %f Zcal = %f\n",DNum,FinalZPos,FinalZPosCal);
      printf(" Down = %f Up = %f Diff = %f\n",EdgeDown[DNum-4][StripNum][BChNum],EdgeUp[DNum-4][StripNum][BChNum],(EdgeDown[DNum-4][StripNum][BChNum]-EdgeUp[DNum-4][StripNum][BChNum]));
    }
  
  //Check if Z Position is within the physical limits of the detector 
  //if not return negative value.
  //Caution!!! This option should not be used in the calibration procedures!!!!
  //Comment out the next statement if you are calibrating detectors.

  Double_t zthresh=0.25; //tolerance in cm
  if (FinalZPosCal<(0-zthresh) || FinalZPosCal>(7.5+zthresh)) {
    FinalZPosCal =sqrt(-1);
  }
 
}

//------------------------------------------------------------------------------------------------//
//Description: This function gets the zero-shift from the alignment file by providing the
//             detector and channel number. The coefficients saved in the alignment file
//             are a and b in the equation volts = a*signal + b.  So the get the zero-shift
//             we do z = -b/a (usually b<0 and a>0, thus z>0).  After this, one subtracts 'z'
//             from the measured signal to get the aligned signal.
void ChannelMap::GetZeroShift (Int_t det, Int_t det_ch, Double_t& zero, Double_t& slope) {
  Int_t detector, channel;
  Double_t b,a, ZS=0;
  for (Int_t i = 0; i<TotalAlignedASICsChannels; i++) {
    IdentifyDetChan(MBID_Align[i],CID_Align[i],ASICs_Ch_Align[i],detector, channel);
    if(det == detector && det_ch == channel) {
      zero = zerosh[i];
      slope = vperch[i];
      if(a!=0) ZS = -b/a;
      break;
    }
  }
}

//------------------------------------------------------------------------------------------------//
//Description: From the channel map provided in the Init method this function get returns the
//             motherboard ID, chip number, and chip channel for a detector number and detector
//             channel number.
void ChannelMap::IdentifyMbChipChan(Int_t det, Int_t det_ch,Int_t &mb_id,Int_t &chip_id,Int_t &asic_ch) {
  for (Int_t i = 0; i<TotalNumberOfChannels; i++) {
    if(Detector[i] == det && Det_Ch[i] == det_ch) {
      mb_id = MBID[i];
      chip_id = CID[i];
      asic_ch = ASICs_Ch[i];
      return;
    }
  }
}

//------------------------PC Relative Gains---------------------added 05/05/2017------------------//
int ChannelMap::Init_PC_UD_RelCal(const char* PC_UD_RelCal_Filename) {
  ifstream pc_ud_rel;
  string line11;
  Double_t ud_slope_dummy, ud_offset_dummy;
  Int_t WireID=0;
  
  pc_ud_rel.open(PC_UD_RelCal_Filename);
  
  for (Int_t i=0; i<WireNum; i++) {   
    PC_UD_Slope[i]=0.0;
    PC_UD_Offset[i]=0.0;
  }
  
  if (pc_ud_rel.is_open()) {
    cout << "Relative Calibration for PC Up_Down ";
    cout << PC_UD_RelCal_Filename << " opened successfully. " << endl;

    getline (pc_ud_rel,line11);//Skips the first line in PC_UD_RelCal_Filename.
    if(doprint) cout<<"line = "<<line11<<endl;

    while (!pc_ud_rel.eof()) {
      pc_ud_rel >> WireID >> ud_slope_dummy>> ud_offset_dummy;
      PC_UD_Slope[WireID] = ud_slope_dummy;
      PC_UD_Offset[WireID] = ud_offset_dummy;
    }
  }
  else LoadFail(PC_UD_RelCal_Filename);
  return 1;
}

//------------------------PC Relative Gains---------------------added 05/05/2017------------------//
void ChannelMap::Get_PC_UD_RelCal(Int_t WireID,Double_t& SlopeUD,Double_t& OffsetUD) {
  SlopeUD = PC_UD_Slope[WireID];
  OffsetUD = PC_UD_Offset[WireID];
}
#endif
