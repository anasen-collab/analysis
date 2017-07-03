/////////////////////////////////////////////////////////////////////////////////////
// This class has been generated on
// Feb 16
// This class describes ANASEN channel map
////////////////////////////////////////////////////////////////////////////////////

//  CHANGES ///////////////////////////////////////////////////////////////////////
//
// 2012-03-13 DSG: Changed Init method form void to int and if the files used
//                as arguments do not exist it now sends an error message but
//                does not stop the Selector. It would be nice to make it stop
//                if such error occurs.
//
// 2012-05-03 JAB: Added overloaded Init method that Grisha and I developed to
//                read X3slopes file and energyslopes file.  Added new method
//                called 'GetX3Slope' that gets the relative slope and the
//                energy slope, multiplies the two, and assigns this value to
//                an input variable (Defined as slope in the method).  Also
//                added several variables to accompany these new methods.
//
// 2012-05-03 DSG: Changed all Float_t to Double_t.
//
// 2012-10-03 KTM: Cleaning up and refactoring. Made function and variable
//                 names more descriptive. Added void function
//                 IdentifyMbChipChan(...);
//                 That allows reverse lookup for a ASICs signal
//                 (Detector Number, Channel Number)
//                                ---> (Motherboard ID,Chip ID, Chip Channel)
//
// 2015-12-09 JCB: Change PC from 19->24 wire
/////////////////////////////////////////////////////////////////////////////////////

#ifndef ChannelMap_h
#define ChannelMap_h

#include <TROOT.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <TRandom3.h>

#include <string>

using namespace std;

#define MaxChNum 500
#define NumX3 24
#define NumQQQ3 4
#define NumDet 28
#define MaxX3Ch 12
#define MaxQQQ3Ch 32
#define MaxADCChannels 200
#define MaxADC 5
#define MaxADCCh 32
#define MaxPulserPeaks 17
#define WireNum 24

class ChannelMap {
  
  Int_t MBID[MaxChNum],CID[MaxChNum],ASICs_Ch[MaxChNum];
  Int_t Detector[MaxChNum],Det_Ch[MaxChNum];
  
  Int_t MBID_Align[MaxChNum],CID_Align[MaxChNum],ASICs_Ch_Align[MaxChNum];
  Int_t TotalNumberOfChannels,TotalAlignedASICsChannels,NumberOfX3AlphaCalibrated;
  Int_t NumberOfQQQ3AlphaCalibrated,NumberOfQQQ3RelativeSlopes,NumberOfAlphaCalibrated;
  Int_t NumberOfX3RelativeSlopes,NumberOfADCChannels;
  
  Int_t ADC[MaxADCChannels];
  Int_t Channel[MaxADCChannels];
  Int_t DetTypeID[MaxADCChannels];
  Int_t Parameter1[MaxADCChannels];
  Int_t Parameter2[MaxADCChannels];
  Int_t BackChNum;
  
  Double_t PCSlope[WireNum], PCShift[WireNum];
  
  string Comment[MaxChNum];
  Double_t SiGains[NumDet];
  Double_t X3RelativeSlope[NumX3][MaxX3Ch],X3FinalFix[NumX3][MaxX3Ch];
  Double_t QQQ3RelativeSlope[NumQQQ3][MaxQQQ3Ch];
  Double_t QQQ3FinalFix[NumQQQ3][MaxQQQ3Ch];
  Double_t EdgeU, EdgeD;
  Double_t zerosh[MaxChNum],vperch[MaxChNum];
  Double_t ZOffset[NumX3],XAt0[NumX3],XAt4[NumX3],YAt0[NumX3],YAt4[NumX3];
  //Double_t Voltage[MaxADC][MaxADCCh][MaxPulserPeaks];
  //Double_t PCPulserPeak[MaxADC][MaxADCCh][MaxPulserPeaks];
  Double_t PCPulser_Offset[MaxADC][MaxADCCh];
  Double_t PCPulser_Slope[MaxADC][MaxADCCh];
  TRandom3 *XRandom2;
  TRandom3 *XYRandom;
  
public:
  //Double_t SiGains[NumDet];
  //Double_t X3RelativeSlope[NumX3][MaxX3Ch],X3FinalFix[NumX3][MaxX3Ch];
  //Double_t QQQ3RelativeSlope[NumQQQ3][MaxQQQ3Ch];
  //Double_t QQQ3FinalFix[NumQQQ3][MaxQQQ3Ch];
  Double_t EdgeUp[NumX3][4][MaxX3Ch], EdgeDown[NumX3][4][MaxX3Ch];

  ChannelMap()
  {
    TotalNumberOfChannels=0;
    TotalAlignedASICsChannels=0;
    NumberOfX3AlphaCalibrated=0;
    NumberOfQQQ3AlphaCalibrated=0;
    NumberOfAlphaCalibrated=0;
    NumberOfX3RelativeSlopes=0;
    NumberOfQQQ3RelativeSlopes=0;
    NumberOfADCChannels=0;
    XRandom2 = new TRandom3();
    XYRandom = new TRandom3();
  };
  ~ChannelMap()
  {
    delete XRandom2;
    delete XYRandom;
  };
  int LoadASICsChannelMapFile(const char* ASICsChannelMapFilename);
  int LoadASICsPulserAlignment(const char* ASICsPulserFilename);
  int LoadSiGains(const char* SiGainsFilename);
  int LoadX3RelativeSlopes(const char* X3SlopeFilename);
  int LoadQQQ3RelativeSlopes(const char* QQQ3SlopeFilename);
  int LoadQQQ3FinalFix(const char* QQQ3FinalFixFilename);
  
  int Init(const char* ASICsChannelMapFilename);
  int Init(const char* ASICsChannelMapFilename, const char* ASICsPulserFilename);
  int Init(const char* ASICsChannelMapFilename, const char* ASICsPulserFilename, const char* SiGainsFilename);
  int Init(const char* ASICsChannelMapFilename, const char* ASICsPulserFilename,
           const char* SiGainsFilename, const char* X3RelativeSlopeFilename);
  int Init(const char* ASICsChannelMapFilename, const char* ASICsPulserFilename,
           const char* SiGainsFilename, const char* X3RelativeSlopeFilename, const char* QQQ3RelativeSlopeFilename);
  int FinalInit(const char* FinalAdjustments, const char* X3Geo);
  int InitWorldCoordinates(const char* WorldCoordinatesFilename);
  
  int GetNumX3Det() {return NumberOfX3AlphaCalibrated;};
  int GetNumQQQ3Det() {return NumberOfQQQ3AlphaCalibrated;};
  int GetNumDet() {return NumberOfAlphaCalibrated;}
  int GetNumberOfADCChannels() {return NumberOfADCChannels;};
  int GetNumberOfMappedASICsChannels() {return TotalNumberOfChannels;};
  int GetNumberOfASICsAligned() {return TotalAlignedASICsChannels;};
  int GetNumberOfAlphaCalibrated() {return NumberOfAlphaCalibrated;};
  int GetNumberOfX3RelativeSlopes() {return NumberOfX3RelativeSlopes;};
  int GetNumberOfQQQ3RelativeSlopes() {return NumberOfQQQ3RelativeSlopes;};
  
  int InitPCADC(const char* PCMapFilename);
  int InitPCCalibration(const char* PCCalibrationFilename);
  int InitPCWireCal(const char* PCWireCalFilename);
  
  void IdentifyADC(Int_t ADCid, Int_t CHid, Int_t& DetType);
  void IdentifyWire(Int_t ADCid, Int_t CHid, Int_t& WireID, Int_t& Side);
  int ConvertToVoltage(Int_t ADCid, Int_t CHid, Int_t PCData, Double_t& Vcal);
  
  void IdentifyDetChan(Int_t mb_id, Int_t chip_id, Int_t asic_ch, Int_t& det, Int_t& det_ch);
  void IdentifyMbChipChan(Int_t det, Int_t det_ch,Int_t &mb_id,Int_t &chip_id,Int_t &asic_ch);
  
  void AlignASICsChannels(Int_t mb_id_al, Int_t chip_id_al, Int_t asic_ch_al, Double_t& zero, Double_t& gain);
  
  void GetX3MeVPerChannel1(Int_t DN, Int_t DetCh, Double_t& slope);
  void GetX3MeVPerChannel2(Int_t DN, Int_t DetCh, Double_t& slope);
  void GetQQQ3MeVPerChannel1(Int_t DN, Int_t DetCh, Double_t& slope);
  void GetQQQ3MeVPerChannel2(Int_t DN, Int_t DetCh, Double_t& slope);
  void GetX3FinalEnergyOffsetInMeV(Int_t DNum, Int_t ChNum, Double_t& zshift);
  void GetQQQ3FinalEnergyOffsetInMeV(Int_t DNum, Int_t ChNum, Double_t& zshift);
  void PosCal(Int_t DNum, Int_t StripNum, Int_t ChNum, Double_t FinalZPos, Double_t& FinalZPosCal);
  void GetX3WorldCoordinates(Int_t DID, Double_t SiX, Double_t SiZ, Double_t& WSiX, Double_t& WSiY, Double_t& WSiZ, Double_t& WSiR, Double_t& WSiPhi);
  void GetQQQ3WorldCoordinates(Int_t DID, Double_t SiX, Double_t SiY, Double_t& WSiX, Double_t& WSiY, Double_t& WSiR, Double_t& WSiPhi );
  //void GetPCWorldCoordinates(Int_t wireid, Double_t zpos, Double_t& xw, Double_t& yw, Double_t& zw);
  void GetPCWorldCoordinates(Int_t wireid, Double_t zpos, Double_t& xw, Double_t& yw, Double_t& zw, Double_t& rw, Double_t& phiw);
  void GetZeroShift(Int_t det, Int_t det_ch, Double_t& zero, Double_t& slope);
  Double_t GetRelLinCoeff(Int_t det, Int_t det_ch) 	{ return X3RelativeSlope[det-4][det_ch]; };
  
};

//------------------------------------------------------------------------------------------------//
int ChannelMap::LoadASICsChannelMapFile (const char* ASICsChannelMapFilename)
{
  ifstream channelmapfile;
  string line;
  
  channelmapfile.open(ASICsChannelMapFilename);
  
  cout << ASICsChannelMapFilename << endl;
	
  if (channelmapfile.is_open())
  {
    cout << "The channel map file " << ASICsChannelMapFilename << " opened successfully." << endl;
    getline (channelmapfile,line);
    Int_t i=0;
    string dummy_det;//Added by JJPIV
    while (!channelmapfile.eof())
    {
      channelmapfile >> MBID[i] >> CID[i] >> ASICs_Ch[i] >> Detector[i] >> Det_Ch[i] >> Comment[i];
      //cout << MBID[i] << "   " << CID[i] << "   " << dummy_det << endl;
      //break;
      //if statement added by JJPIV
      //Nabin's map has an extra column labled DetType that Jeff's map/code doesn't have
      //we read in that column into a dummy variable
      //if the det type is an SX3, we increment it by 4 to match with Jeff's map
      //cout << dummy_det << endl;
      /*
      if (dummy_det == "SX3"){
	//cout << "SX3 hit\n";
	Detector[i] += 4;
      }
      */
      i++;
    }
    TotalNumberOfChannels = i-1;
  }
  else
  {
    cout << "*** ERROR: File " << ASICsChannelMapFilename << " could not be opened!" << endl;
    return 0;
  }
  channelmapfile.close();
  return 1;
}

//------------------------------------------------------------------------------------------------//

int ChannelMap::LoadASICsPulserAlignment (const char* ASICsPulserFilename)
{
  ifstream alignchan;
  alignchan.open(ASICsPulserFilename);
  if (alignchan.is_open())
  {
    cout << "The channel allignment file " << ASICsPulserFilename << " opened successfully." << endl;
    Int_t j=0;
    while (!alignchan.eof())
    {
      alignchan >>  MBID_Align[j] >> CID_Align[j] >> ASICs_Ch_Align[j] >> zerosh[j] >> vperch[j];
      j++;
    }
    TotalAlignedASICsChannels = j-1;
  }
  else
  {
    cout << "*** ERROR: File " << TotalAlignedASICsChannels << " could not be opened!" << endl;
    return 0;
  }
  
  return 1;
}

//------------------------------------------------------------------------------------------------//

int ChannelMap::LoadSiGains(const char* SiGainsFilename)
{
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
  if (SiGainsFile.is_open())
  {
    cout << "The Si Gains file " << SiGainsFilename << " opened successfully." << endl;
    i=0;
    j=0;
    Int_t Dnum = 0;
    Int_t ChNum = 0;
    Double_t dummy = 0.0;
    
    while (!SiGainsFile.eof())
    {
      if(Dnum>=4 && Dnum<=27) {i++;}
      else if (Dnum<=3) {j++;}
      else {continue;}
      SiGainsFile >>  Dnum >> ChNum >> dummy;
      SiGains[Dnum] = dummy;
    }
    
    NumberOfX3AlphaCalibrated = i-1; // The Minus one accounts for the end of line character
    NumberOfQQQ3AlphaCalibrated = j;
    NumberOfAlphaCalibrated = NumberOfX3AlphaCalibrated + NumberOfQQQ3AlphaCalibrated;
  }
  else
  {
    cout << "\t*** ERROR: File " << SiGainsFilename << " could not be opened!" << endl;
    return 0;
  }
  
  SiGainsFile.close();
  
  return 1;
}

//------------------------------------------------------------------------------------------------//

int ChannelMap::LoadX3RelativeSlopes(const char* X3RelativeSlopeFilename)
{
  ifstream X3RelativeSlopeFile;
  string line;
  X3RelativeSlopeFile.open(X3RelativeSlopeFilename);
  
  
  if (X3RelativeSlopeFile.is_open())
  {
    cout << "The X3 slopes file " << X3RelativeSlopeFilename << " opened successfully." << endl;
    Int_t m=0;
    Int_t DnumX3 = 0;
    Int_t X3ChNum = 0;
    Double_t dummy = 0.0;
    
    while (!X3RelativeSlopeFile.eof())
    {
      X3RelativeSlopeFile >>  DnumX3 >> X3ChNum >> dummy;
      X3RelativeSlope[DnumX3-4][X3ChNum] = dummy;
      m++;
    }
    NumberOfX3RelativeSlopes = m-1;
    
  }
  else
  {
    cout << "Cannot open X3 slopes file " << X3RelativeSlopeFilename << endl;
    cout << "Macro terminated abnormally !!!" << endl;
    return 0;
  }
  
  X3RelativeSlopeFile.close();
  
  return 1;
}

//------------------------------------------------------------------------------------------------//

int ChannelMap::LoadQQQ3RelativeSlopes(const char* QQQ3RelativeSlopeFilename)
{
  ifstream QQQ3RelativeSlopeFile;
  string line;
  QQQ3RelativeSlopeFile.open(QQQ3RelativeSlopeFilename);
  
  
  if (QQQ3RelativeSlopeFile.is_open())
  {
    cout << "The QQQ3 slopes file " << QQQ3RelativeSlopeFilename << " opened successfully." << endl;
    Int_t m=0;
    Int_t DnumQQQ3 = 0;
    Int_t QQQ3ChNum = 0;
    Double_t dummy = 0.0;
    
    while (!QQQ3RelativeSlopeFile.eof())
    {
      QQQ3RelativeSlopeFile >>  DnumQQQ3 >> QQQ3ChNum >> dummy;
      QQQ3RelativeSlope[DnumQQQ3][QQQ3ChNum] = dummy;
      m++;
    }
    NumberOfQQQ3RelativeSlopes = m-1;
    
  }
  else
  {
    cout << "Cannot open QQQ3 slopes file " << QQQ3RelativeSlopeFilename << endl;
    cout << "Macro terminated abnormally !!!" << endl;
    return 0;
  }
  
  QQQ3RelativeSlopeFile.close();
  
  return 1;
  
}

//------------------------------------------------------------------------------------------------//

int ChannelMap::LoadQQQ3FinalFix(const char* QQQ3FinalFixFilename)
{
  ifstream QQQ3FinalFixFile;
  QQQ3FinalFixFile.open(QQQ3FinalFixFilename);
  Double_t Zero = 0;
  Int_t DNum,ChNum;
  
  for (Int_t i=0; i<4; i++)
  {
    for (Int_t c=0; c<32; c++)
    {
      QQQ3FinalFix[i][c] = 0;
    }
  }
  
  if (QQQ3FinalFixFile.is_open())
  {
    cout << "File with QQQ3 final fix zero shifts " << QQQ3FinalFixFilename;
    cout << " opened successfully." << endl;
  }
  else
  {
    cout << "Cannot open QQQ3 final fix zero shifts file " << QQQ3FinalFixFilename << endl;
    cout << "Macro terminated abnormally!!! " << endl;
    return 0;
  }
  
  while (!QQQ3FinalFixFile.eof())
  {
    QQQ3FinalFixFile >> DNum >> ChNum >> Zero;
    QQQ3FinalFix[DNum][ChNum] = Zero;
  }
  
  return 1;
  
}

//------------------------------------------------------------------------------------------------//

//Description: Initialize only with the channel map.  This is meant to be used in
//             codes where the channel alignment has not yet been made.
int ChannelMap::Init (const char* ASICsChannelMapFilename)
{
  Int_t status = LoadASICsChannelMapFile(ASICsChannelMapFilename);
  return status;
}

//------------------------------------------------------------------------------------------------//

//Description: This is meant to be used when the zero-shifts
//             and the alignment coefficients have already been obtaied.  Most likely
//             the user would want to get the energy calibration (linear) coefficients
//             if he/she is using this function.

int ChannelMap::Init (const char* ASICsChannelMapFilename, const char* ASICsPulserFilename)
{
  Int_t status = 0;
  status = LoadASICsChannelMapFile(ASICsChannelMapFilename);
  
  if(status == 1)
  {
    status = LoadASICsPulserAlignment(ASICsPulserFilename);
  }
  else
  {
    cout << "LoadASICsChannelMapFile Failed" << endl;
    cout << "Macro terminated abnormally !!!" << endl;
    return 0;
  }
  
  return status;
}

//------------------------------------------------------------------------------------------------//

int ChannelMap::Init (const char* ASICsChannelMapFilename, const char* ASICsPulserFilename, const char* SiGainsFilename)
{
  Int_t status = 0;
  status = LoadASICsChannelMapFile(ASICsChannelMapFilename);
  if(status == 1)
  {
    status = LoadASICsPulserAlignment(ASICsPulserFilename);
  }
  else
  {
    cout << "LoadASICsChannelMapFile Failed" << endl;
    cout << "Macro terminated abnormally !!!" << endl;
    return 0;
  }
  
  if(status == 1)
  {
    status = LoadSiGains(SiGainsFilename);
  }
  else
  {
    cout << "LoadASICsPulserAlignment Failed" << endl;
    cout << "Macro terminated abnormally !!!" << endl;
    return 0;
  }
  
  return status;
}

//------------------------------------------------------------------------------------------------//

//Description:
//			In this Init method you need 4 files: the channelmap file, the alignchannels file
//			(produced from pulser alignment), the energyslopes file
//			(currently produced by hand from alpha data using only the reference channel for each detector),
//			and the X3slope file (produced from the X3ChargeAlignment.C code).  Therefore, you will be using this
//          version of Init only if you have all of these things and are ready to look at energy calibrated data
//			in the X3's.

int ChannelMap::Init(const char* ASICsChannelMapFilename, const char* ASICsPulserFilename,
                     const char* SiGainsFilename, const char* X3RelativeSlopeFilename)
{
  Int_t status = 0;
  status = LoadASICsChannelMapFile(ASICsChannelMapFilename);
  if(status == 1)
  {
    status = LoadASICsPulserAlignment(ASICsPulserFilename);
  }
  else
  {
    cout << "LoadASICsChannelMapFile Failed" << endl;
    cout << "Macro terminated abnormally !!!" << endl;
    return 0;
  }
  
  if(status == 1)
  {
    status = LoadSiGains(SiGainsFilename);
  }
  else
  {
    cout << "LoadASICsPulserAlignment Failed" << endl;
    cout << "Macro terminated abnormally !!!" << endl;
    return 0;
  }
  
  if(status == 1)
  {
    status = LoadX3RelativeSlopes(X3RelativeSlopeFilename);
  }
  else
  {
    cout << "LoadSiGains Failed" << endl;
    cout << "Macro terminated abnormally !!!" << endl;
  }
  
  return status;
}

//------------------------------------------------------------------------------------------------//

//Description:


int ChannelMap::Init(const char* ASICsChannelMapFilename, const char* ASICsPulserFilename,
                     const char* SiGainsFilename, const char* X3RelativeSlopeFilename, const char* QQQ3RelativeSlopeFilename)
{
  Int_t status = 0;
  status = LoadASICsChannelMapFile(ASICsChannelMapFilename);
  if(status == 1)
  {
    status = LoadASICsPulserAlignment(ASICsPulserFilename);
  }
  else
  {
    cout << "LoadASICsChannelMapFile Failed" << endl;
    cout << "Macro terminated abnormally !!!" << endl;
    return 0;
  }
  
  if(status == 1)
  {
    status = LoadSiGains(SiGainsFilename);
  }
  else
  {
    cout << "LoadASICsPulserAlignment Failed" << endl;
    cout << "Macro terminated abnormally !!!" << endl;
    return 0;
  }
  
  if(status == 1)
  {
    status = LoadX3RelativeSlopes(X3RelativeSlopeFilename);
  }
  else
  {
    cout << "LoadSiGains Failed" << endl;
    cout << "Macro terminated abnormally !!!" << endl;
  }
  
  if(status == 1)
  {
    status = LoadQQQ3RelativeSlopes(QQQ3RelativeSlopeFilename);
  }
  else
  {
    cout << "LoadX3RelativeSlopes Failed" << endl;
    cout << "Macro terminated abnormally !!!" << endl;
  }
  
  if(status == 0)
  {
    cout << "LoadQQQ3RelativeSlopes Faild" << endl;
    cout << "Macro terminated abnormally !!!" << endl;
  }
  
  return status;
}

//------------------------------------------------------------------------------------------------//

// Description: Fine tuning parameters for better energy calibration
// The file contains detector number, channel number and zero shift in MeV
// By default, the zero shift is zero. Only those channels that have non-zero
// shift need to be described in the final fix file.

int ChannelMap::FinalInit(const char* FinalFixFilename, const char* X3GeoFilename)
{
  ifstream finalfix;
  ifstream x3geo;
  finalfix.open(FinalFixFilename);
  x3geo.open(X3GeoFilename);
  Int_t DNum,ChNum,StripNum;
  Double_t Zero;
  
  for (Int_t i=0; i<NumX3; i++)
  {
    for (Int_t c=0; c<MaxX3Ch; c++)
    {
      X3FinalFix[i][c] = 0;
    }
  }
  
  if (finalfix.is_open())
  {
    cout << "File with final fix zero shifts " << FinalFixFilename;
    cout << " opened successfully." << endl;
  }
  else
  {
    cout << "Cannot open final fix zero shifts file " << FinalFixFilename << endl;
    cout << "Macro terminated abnormally!!! " << endl;
    return 0;
  }
  
  while (!finalfix.eof())
  {
    finalfix >> DNum >> ChNum >> Zero;
    X3FinalFix[DNum-4][ChNum] = Zero;
  }
  
  
  if(x3geo.is_open())
  {
    cout << "File with X3 Geometry " << X3GeoFilename;
    cout << " opened successfully." << endl;
    while (!x3geo.eof())
    {
      x3geo >> DNum >> StripNum >> BackChNum >> EdgeD >> EdgeU;
      EdgeDown[DNum-4][StripNum][BackChNum] = EdgeD;
      EdgeUp[DNum-4][StripNum][BackChNum] = EdgeU;
    }
    //cout << "X3 geometry adjustment using file " << X3GeoFilename << " was performed." << endl;
  }
  else
  {
    cout << "*** ERROR: File " << x3geo << " could not be opened!" << endl;
    cout << "Macro terminated abnormally!!! " << endl;
    return 0;
  }
  finalfix.close();
  x3geo.close();
  return 1;
}

//------------------------------------------------------------------------------------------------//

int ChannelMap::InitWorldCoordinates(const char* WorldCoordinatesFilename)
{
  ifstream worldCfile;
  worldCfile.open(WorldCoordinatesFilename);
  
  Int_t Dnum;
  Double_t Zoff,XMin,XMax,YMin,YMax;
  string DumComment,line;
  
  // The ZOffset is given with respect to the face of the forward Si detectors.
  // Note that if initial ZOffset is negative then the detector is considered as missing.
  // If YAt0 and YAt4 are both more than +50.0 cm then this is the forward detector.
  
  for (Int_t i = 0; i<NumX3; i++)
  {
    ZOffset[i]  = -100.0;
    XAt0[i]     = 0;
    XAt4[i]     = 0;
    YAt0[i]     = 0;
    YAt4[i]     = 0;
  }
  
  
  if (worldCfile.is_open())
  {
    cout << "File with world coordinates " << WorldCoordinatesFilename;
    cout << " opened successfully." << endl;
    getline (worldCfile,line);
    // cout << line << endl;
    while (!worldCfile.eof())
    {
//      worldCfile >> Dnum >> Zoff >> XMin >> XMax >> YMin >> YMax >> DumComment;
//  Removed comment from input file - JCB 12/10/15
        worldCfile >> Dnum >> Zoff >> XMin >> XMax >> YMin >> YMax ;
      ZOffset[Dnum-4]  = Zoff;
      XAt0[Dnum-4]     = XMin;
      XAt4[Dnum-4]     = XMax;
      YAt0[Dnum-4]     = YMin;
      YAt4[Dnum-4]     = YMax;
    }
    
  }
  else
  {
    cout << "Cannot open file with world coordinates " << WorldCoordinatesFilename << endl;
    cout << "Macro terminated abnormally!!! " << endl;
    return 0;
  }
  worldCfile.close();
  return 1;
}

//------------------------------------------------------------------------------------------------//

int ChannelMap::InitPCADC(const char* PCMapFilename)
{
  ifstream pcmap;
  string line;
  Int_t DumADCID,DumCh,DumWireID,DumSide;
  
  pcmap.open(PCMapFilename);
  
  if (pcmap.is_open())
  {
    cout << "Channel map for proportional counter " << PCMapFilename << " opened sucessfully. " << endl;
    getline (pcmap,line);
    while (!pcmap.eof())
    {
      pcmap >> DumADCID >> DumCh >> DumWireID >> DumSide ;
      ADC[NumberOfADCChannels] = DumADCID;
      Channel[NumberOfADCChannels] = DumCh;
      DetTypeID[NumberOfADCChannels] = 1;
      //cout << "Goody" << DumADCID << DumCh << endl;
      Parameter1[NumberOfADCChannels] = DumWireID;
      Parameter2[NumberOfADCChannels] = DumSide;
      NumberOfADCChannels++;
    }
    NumberOfADCChannels--;
  }
  else
  {
    cout << "Cannot open file with Proportional Counter Map " << PCMapFilename << endl;
    cout << "Macro terminated abnormally!!! " << endl;
    return 0;
  }
  pcmap.close();
  return 1;
}

//------------------------------------------------------------------------------------------------//

int ChannelMap::InitPCCalibration(const char* PCCalibrationFilename)
{
  ifstream pccal;
  string line;
  Int_t adcid,chnum;
  Double_t dumvoltage,dumpeak;
  
  pccal.open(PCCalibrationFilename);
  
  for (Int_t k=0; k<MaxADC; k++){
    for (Int_t i=0; i<MaxADCCh; i++){
      for (Int_t j=0; j<MaxPulserPeaks; j++){
	//Voltage[k][i][j] = 0.;
	//PCPulserPeak[k][i][j] = 0;
      }
      PCPulser_Offset[k][i] = 0;
      PCPulser_Slope[k][i] = 0;
    }
  }
  
  if (pccal.is_open())
  {
    
    cout << "Pulser Calibration for proportional counter " << PCCalibrationFilename << " opened successfully. " << endl;
    getline (pccal,line);
    while (!pccal.eof())
    {
      pccal >> adcid >> chnum >> dumvoltage >> dumpeak;
      //Int_t i = 0;
      //while (Voltage[adcid][chnum][i]>0) i++;
      //Voltage[adcid][chnum][i] = dumvoltage;
      //PCPulserPeak[adcid][chnum][i] = dumpeak;
      PCPulser_Offset[adcid][chnum] = dumvoltage;
      PCPulser_Slope[adcid][chnum] = dumpeak;
    }
  }
  else
  {
    cout << "\n**************************************** " << endl;
    cout << "ERROR: Pulser Calibration file: " << PCCalibrationFilename << " failed to open!" << endl;
    cout << "****************************************\n" << endl;
    cout << "Pulser Calibration for proportional counter " << endl;
  }
  
  return 1;
}

//------------------------------------------------------------------------------------------------//

int ChannelMap::InitPCWireCal(const char* PCWireCalFilename)
{
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
    getline (pcwirecal,line);
    while (!pcwirecal.eof()) {
      pcwirecal >> WireNumber >> pcslopedum >> pcshiftdum;
      PCSlope[WireNumber] = pcslopedum;
      PCShift[WireNumber]	= pcshiftdum;
    }
  }
  else cout << "DIDN'T OPEN FILE!" << endl;
  
  return 1;
}

//------------------------------------------------------------------------------------------------//

void ChannelMap::IdentifyADC(Int_t ADCid, Int_t CHid, Int_t& DetType)
{
  DetType = 0;
  for (Int_t i=0; i<NumberOfADCChannels; i++)
  {
    if (ADCid == ADC[i] && CHid == Channel[i])  DetType=DetTypeID[i];
  }
  
}

//------------------------------------------------------------------------------------------------//

void ChannelMap::IdentifyWire(Int_t ADCid, Int_t CHid, Int_t& WireID, Int_t& Side)
{
  for (Int_t i=0; i<NumberOfADCChannels; i++)
  {
    if (ADCid == ADC[i] && CHid == Channel[i])
    {
      WireID = Parameter1[i];
      Side   = Parameter2[i];
    }
  }
  
}

//------------------------------------------------------------------------------------------------//

Int_t ChannelMap::ConvertToVoltage(Int_t ADCid, Int_t CHid, Int_t PCData, Double_t& Vcal)
{
  XRandom2->SetSeed();
  Int_t i=0;
  Vcal = 0;
  if (PCData > 0){

    Vcal = PCData*PCPulser_Slope[ADCid][CHid] + PCPulser_Offset[ADCid][CHid];

    /*
    while ( ((Double_t)PCData>PCPulserPeak[ADCid][CHid][i]) && i<MaxPulserPeaks )
    {
      i++;
    }
    
    if (i>0)
    {
      Vcal = (PCData-PCPulserPeak[ADCid][CHid][i-1]+XRandom2->Rndm()-0.5)/
	    (PCPulserPeak[ADCid][CHid][i]-PCPulserPeak[ADCid][CHid][i-1]+XRandom2->Rndm()-0.5)*
	    (Voltage[ADCid][CHid][i]-Voltage[ADCid][CHid][i-1])+
	    Voltage[ADCid][CHid][i-1];
    }
    else
    {
      Vcal = (Voltage[ADCid][CHid][1]-Voltage[ADCid][CHid][0])/
	    (PCPulserPeak[ADCid][CHid][1]-PCPulserPeak[ADCid][CHid][0]+XRandom2->Rndm()-0.5)*
	    (PCData - PCPulserPeak[ADCid][CHid][1]+XRandom2->Rndm()-0.5) + Voltage[ADCid][CHid][1];
      return 1;
    }
    */
  }
  return 1;
}

//------------------------------------------------------------------------------------------------//

//Description: From the channel map provided in the Init method this function get returns the
//             detector number and detector channel number for a given mother board, chip and
//             chip channel.

void ChannelMap::IdentifyDetChan(Int_t mb_id, Int_t chip_id, Int_t asic_ch, Int_t& det, Int_t& det_ch)
{
  for (Int_t i = 0; i<TotalNumberOfChannels; i++)
  {
    if (MBID[i] == mb_id && CID[i] == chip_id && ASICs_Ch[i] == asic_ch)
    {
      det = Detector[i];
      det_ch = Det_Ch[i];
      return;
    }
  }
}

//------------------------------------------------------------------------------------------------//

//Description: From the channel map provided in the Init method this function get returns the
//             motherboard ID, chip number, and chip channel for a detector number and detector
//             channel number.

void ChannelMap::IdentifyMbChipChan(Int_t det, Int_t det_ch,Int_t &mb_id,Int_t &chip_id,Int_t &asic_ch)
{
  for (Int_t i = 0; i<TotalNumberOfChannels; i++)
  {
    if(Detector[i] == det && Det_Ch[i] == det_ch)
    {
      mb_id = MBID[i];
      chip_id = CID[i];
      asic_ch = ASICs_Ch[i];
      return;
    }
  }
}

//------------------------------------------------------------------------------------------------//

void ChannelMap::AlignASICsChannels(Int_t mb_id, Int_t chip_id, Int_t asic_ch, Double_t& zero, Double_t& gain)
{
  for (Int_t i = 0; i<TotalAlignedASICsChannels; i++)
  {
    if (MBID_Align[i] == mb_id && CID_Align[i] == chip_id && ASICs_Ch_Align[i] == asic_ch)
    {
      zero = zerosh[i];
      gain = vperch[i];
    }    
  }
}

//------------------------------------------------------------------------------------------------//

//Description: This function gets the zero-shift from the alignment file by providing the
//             detector and channel number. The coefficients saved in the alignment file
//             are a and b in the equation volts = a*signal + b.  So the get the zero-shift
//             we do z = -b/a (usually b<0 and a>0, thus z>0).  After this, one subtracts 'z'
//             from the measured signal to get the aligned signal.
void ChannelMap::GetZeroShift (Int_t det, Int_t det_ch, Double_t& zero, Double_t& slope)
{
  Int_t detector, channel;
  Double_t b,a, ZS=0;
  for (Int_t i = 0; i<TotalAlignedASICsChannels; i++)
  {
    IdentifyDetChan(MBID_Align[i],CID_Align[i],ASICs_Ch_Align[i],detector, channel);
    if(det == detector && det_ch == channel)
    {
      zero = zerosh[i];
      slope = vperch[i];
      if(a!=0) ZS = -b/a;
      break;
    }
  }
  
  //return ZS;
}

//------------------------------------------------------------------------------------------------//

void ChannelMap::PosCal(Int_t DNum, Int_t StripNum, Int_t BChNum, Double_t FinalZPos, Double_t& FinalZPosCal)
{
  
  Double_t EdgeDCal=-100.0;
  Double_t EdgeUCal=100.00;
  
  switch(BChNum)
  {
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

  //if (DNum==5){
    //cout << EdgeDCal << "  " << EdgeUCal << "  " << EdgeDown[DNum-4][StripNum][BackChNum] << " " << EdgeUp[DNum-4][StripNum][BackChNum] << " " << FinalZPos << " " << FinalZPosCal << endl;
//cout << FinalZDown << "  " << FinalZDownCal << "  " << Si.DownChNum[Si.NSiHits] << "  " << Si.BackChNum[Si.NSiHits] << endl;
  //}
  
  //Check if Z Position is within the physical limits of the detector (with 1 mm tolerance),
  //if not return negative value.
  //Caution!!! This option should not be used in the calibration procedures!!!!
  //Comment out the next statement if you are calibrating detectors.
  
  if (FinalZPosCal<-0.1 || FinalZPosCal>7.6) {FinalZPosCal =-10;}
  /*
  if (DNum==4){
    //cout << StripNum << " " << BackChNum << " " << EdgeDown[DNum][StripNum][BackChNum] << " " << EdgeUp[DNum][StripNum][BackChNum] << " " << FinalZPos << " " << FinalZPosCal << endl;
  }
  */
}

//------------------------------------------------------------------------------------------------//

void ChannelMap::GetX3WorldCoordinates(Int_t DID, Double_t SiX, Double_t SiZ, Double_t& WSiX,
                                       Double_t& WSiY, Double_t& WSiZ, Double_t& WSiR, Double_t& WSiPhi)
{
  if (ZOffset[DID-4]>-0.1 && YAt0[DID-4]<50. && YAt0[DID-4]<50.){
    //This is barrel detector
    //WSiZ = ZOffset[DID-4] + 7.5 - SiZ;
    WSiZ = ZOffset[DID-4] + SiZ;//JJPIV changed on 2-11-2016.
    //cout << "Z: " << SiZ << "   Z Offset: " << ZOffset[DID-4] << "   Det ID: " << DID << endl;
    WSiX = (XAt4[DID-4] - XAt0[DID-4])*0.25*SiX + XAt0[DID-4];
    WSiY = (YAt4[DID-4] - YAt0[DID-4])*0.25*SiX + YAt0[DID-4];
    WSiR = TMath::Sqrt(WSiX*WSiX + WSiY*WSiY);
    if(WSiY >= 0 && WSiX>=0){
      WSiPhi = TMath::ATan(WSiY/WSiX);
    }else if(WSiY >= 0 && WSiX<0){
      WSiPhi = TMath::ATan(WSiY/WSiX) + TMath::Pi();
    }else if(WSiY < 0 && WSiX<0){
      WSiPhi = TMath::ATan(WSiY/WSiX) + TMath::Pi();
    }else if(WSiY < 0 && WSiX>=0){
      WSiPhi = TMath::ATan(WSiY/WSiX) + 2*TMath::Pi();
    }
    /*
    if (DID==4 && WSiZ !=23.7){
      cout << DID << "   " << WSiX << "   " << WSiY << "   " << WSiZ << endl;
    }
    */
  }else{
    //This detector does not have world coordinates description
    WSiX = 0;
    WSiY = 0;
    WSiZ = 0;
  }
}

//------------------------------------------------------------------------------------------------//

void ChannelMap::GetQQQ3WorldCoordinates(Int_t DID, Double_t SiX, Double_t SiY, Double_t& WSiX, Double_t& WSiY, Double_t& WSiR, Double_t& WSiPhi)
{
  //Double_t Theta = TMath::Pi()*(DID+1)/2; //Rotate detectors CCW by an angle Theta
  //changed on Feb27 to reflect reality--JJPIV
  Double_t Theta = 0;
  if (DID==0){
    Theta = 0;
  }else if (DID==1){
    Theta = 3*TMath::Pi()/2;
  }else if (DID==2){
    Theta = TMath::Pi();
  }else if (DID==3){
    Theta = TMath::Pi()/2;
  }
  //cout << endl << endl;
  //cout << DID << endl;
  //cout << SiX << "  " << SiY << "  " << Theta << endl;
  WSiX = SiX*TMath::Cos(Theta)-SiY*TMath::Sin(Theta);
  WSiY = SiX*TMath::Sin(Theta)+SiY*TMath::Cos(Theta);
  //WSiY = SiX*TMath::Cos(Theta)-SiY*TMath::Sin(Theta);
  //WSiX = SiX*TMath::Sin(Theta)+SiY*TMath::Cos(Theta);

  WSiR = TMath::Sqrt(WSiX*WSiX + WSiY*WSiY);
  //JJPIV--added 6Jan2016
  //Phi calculated to go from 0 to 2*pi--0 is along positive x
  //in lab frame, x is North, z is upstream, and y is up

  if(WSiY >= 0 && WSiX>=0){
    WSiPhi = TMath::ATan(WSiY/WSiX);
  }else if(WSiY >= 0 && WSiX<0){
    WSiPhi = TMath::ATan(WSiY/WSiX) + TMath::Pi();
  }else if(WSiY < 0 && WSiX<0){
    WSiPhi = TMath::ATan(WSiY/WSiX) + TMath::Pi();
  }else if(WSiY < 0 && WSiX>=0){
    WSiPhi = TMath::ATan(WSiY/WSiX) + 2*TMath::Pi();
  }
  //cout << WSiX << "  " << WSiY << "  " << WSiPhi << "  " << endl;

}

//------------------------------------------------------------------------------------------------//

void ChannelMap::GetPCWorldCoordinates(Int_t wireid, Double_t zpos, Double_t& xw, Double_t& yw, Double_t& zw, Double_t& rw, Double_t& phiw)
{
  XYRandom->SetSeed();
  // Changed to new PC geometry  Radius 3.0->3.75
    //  Double_t Radius = 3.0;
  Double_t Radius = 3.75;
  Double_t Angle = (24-(Double_t)wireid+XYRandom->Rndm()-0.5)*TMath::TwoPi()/24.0 + TMath::Pi()/2;
  //changed angle on Feb27 to correspond to reality
  xw = Radius*TMath::Cos(Angle);
  yw = Radius*TMath::Sin(Angle);
  zw = zpos*PCSlope[wireid]+PCShift[wireid];
  rw = TMath::Sqrt(xw*xw + yw*yw);
  /*
  if(yw > 0 && xw>0){
    phiw = TMath::ATan(xw/yw);
  }else if(yw > 0 && xw<0){
    phiw = TMath::ATan(-xw/yw) + 0.5*TMath::Pi();
  }else if(yw < 0 && xw<0){
    phiw = TMath::ATan(xw/yw) + TMath::Pi();
  }else if(yw < 0 && xw>0){
    phiw = TMath::ATan(-xw/yw) + 1.5*TMath::Pi();
  }
  */
  //phiw = Angle-TMath::Pi()/2;
  phiw = Angle;//Change by JJPIV on Feb8
  if (phiw > 2*TMath::Pi()){
    phiw = phiw - 2*TMath::Pi();
  }
}

//------------------------------------------------------------------------------------------------//

void ChannelMap::GetX3MeVPerChannel1(Int_t DN, Int_t DetCh, Double_t& slope)
{
  slope = X3RelativeSlope[DN-4][DetCh];
}
void ChannelMap::GetX3MeVPerChannel2(Int_t DN, Int_t DetCh, Double_t& slope)
{
  slope = SiGains[DN];
}
//------------------------------------------------------------------------------------------------//

void ChannelMap::GetQQQ3MeVPerChannel1(Int_t DN, Int_t DetCh, Double_t& slope)
{
  slope = QQQ3RelativeSlope[DN][DetCh];
}
void ChannelMap::GetQQQ3MeVPerChannel2(Int_t DN, Int_t DetCh, Double_t& slope)
{
  slope = SiGains[DN];
}

//------------------------------------------------------------------------------------------------//

void ChannelMap::GetX3FinalEnergyOffsetInMeV(Int_t DNum, Int_t ChNum, Double_t& zshift)
{
  zshift = X3FinalFix[DNum-4][ChNum];
}

//------------------------------------------------------------------------------------------------//

void ChannelMap::GetQQQ3FinalEnergyOffsetInMeV(Int_t DNum, Int_t ChNum, Double_t& zshift)
{
  zshift = QQQ3FinalFix[DNum][ChNum];
}

//------------------------------------------------------------------------------------------------//
#endif
