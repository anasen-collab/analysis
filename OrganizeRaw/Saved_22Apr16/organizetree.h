//////////////////////////////////////////////////////////////////////////////////////
#include <TROOT.h>
#include <vector>
using namespace std;

// 12/09/15 Modified for 24 wires in PC
//2016/01/04 Modified for different X, Y, Z for Silicon & PC
//////////////////////////////////////////////////////////////////////////////////////
#define MaxPCHits 24
#define NPCWires  24
#define MaxSiHits   500
//////////////////////////////////////////////////////////////////////////////////////
//gROOT->ProcessLine("#include <vector>");

class SiHit {
  // This class is for Si hit
 public:
  SiHit(){};
  /*
  Int_t NSiHits;
  Int_t NHitsInDet[MaxSiHits];
  Int_t DetID[MaxSiHits][2];
  Int_t UpMult[MaxSiHits];
  Int_t DownMult[MaxSiHits];
  Int_t BackMult[MaxSiHits];
  Int_t HitType[MaxSiHits][2];

  Int_t UpChNum[MaxSiHits][16];
  Int_t DownChNum[MaxSiHits][16];
  Int_t BackChNum[MaxSiHits][16];

  Double_t EnergyUp_Raw[MaxSiHits][16];
  Double_t EnergyDown_Raw[MaxSiHits][16];
  Double_t EnergyBack_Raw[MaxSiHits][16];
  Double_t EnergyUp_Pulser[MaxSiHits][16];
  Double_t EnergyDown_Pulser[MaxSiHits][16];
  Double_t EnergyBack_Pulser[MaxSiHits][16];
  Double_t EnergyUp_Rel[MaxSiHits][16];
  Double_t EnergyDown_Rel[MaxSiHits][16];
  Double_t EnergyBack_Rel[MaxSiHits][16];
  Double_t EnergyUp_Cal[MaxSiHits][16];
  Double_t EnergyDown_Cal[MaxSiHits][16];
  Double_t EnergyBack_Cal[MaxSiHits][16];
  */
  Int_t NSiHits;
  vector<Int_t> NHitsInDet;
  vector<vector<Int_t> > DetID;

  vector<Int_t> UpMult;
  vector<Int_t> DownMult;
  vector<Int_t> BackMult;

  vector<vector<Int_t> > HitType;
  vector<vector<Int_t> > UpChNum;
  vector<vector<Int_t> > DownChNum;
  vector<vector<Int_t> > BackChNum;

  vector<vector<Double_t> > EnergyUp_Raw;
  vector<vector<Double_t> > EnergyDown_Raw;
  vector<vector<Double_t> > EnergyBack_Raw;
  vector<vector<Double_t> > EnergyUp_Pulser;
  vector<vector<Double_t> > EnergyDown_Pulser;
  vector<vector<Double_t> > EnergyBack_Pulser;
  vector<vector<Double_t> > EnergyUp_Rel;
  vector<vector<Double_t> > EnergyDown_Rel;
  vector<vector<Double_t> > EnergyBack_Rel;
  vector<vector<Double_t> > EnergyUp_Cal;
  vector<vector<Double_t> > EnergyDown_Cal;
  vector<vector<Double_t> > EnergyBack_Cal;

  vector<vector<Double_t> > TimeUp;
  vector<vector<Double_t> > TimeDown;
  vector<vector<Double_t> > TimeBack;

  vector<vector<Double_t> Energy;
  vector<vector<Double_t> Time;
  vector<vector<Double_t> X, Y, Z;

  vector<vector<Double_t> XW, YW, ZW;
  vector<vector<Double_t> RW, PhiW;
  vector<vector<Double_t> RFSubtract;
  /*
  Double_t TimeUp[MaxSiHits][16];
  Double_t TimeDown[MaxSiHits][16];
  Double_t TimeBack[MaxSiHits][16];

  Double_t Energy[MaxSiHits][2];
  Double_t Time[MaxSiHits][2];
  Double_t X[MaxSiHits][2],Y[MaxSiHits][2],Z[MaxSiHits][2];
  Double_t XW[MaxSiHits][2],YW[MaxSiHits][2],ZW[MaxSiHits][2];
  Double_t RW[MaxSiHits][2], PhiW[MaxSiHits][2];
  Double_t RFSubtract[MaxSiHits][2];
  */
  void zeroSiHit()
  {
    NSiHits = 0;
    EnergyUp_Raw.clear();
    for ( Int_t i=0; i<MaxSiHits; i++ ) {
      NHitsInDet[i] = 0;
      UpMult[i] = 0;
      DownMult[i] = 0;
      BackMult[i] = 0;
      for (Int_t j=0; j<2; j++){
	DetID[i][j] = -1;
	HitType[i][j] = 0;
	X[i][j]  = 0;
	Y[i][j]  = 0;
	Z[i][j]  = 0;
	XW[i][j] = 0;
	YW[i][j] = 0;
	ZW[i][j] = 0;
	RW[i][j] = 0;
	PhiW[i][j] = 0;
	Energy[i][j] = 0;
	Time[i][j] = 0;
	RFSubtract[i][j] = 0;
      }
      for ( Int_t j=0; j<16; j++ ){
	UpChNum[i][j] = -1;
	DownChNum[i][j] = -1;
	BackChNum[i][j] = -1;
	EnergyBack_Raw[i][j] = 0;
	//EnergyUp_Raw[i][j] = 0;
	EnergyDown_Raw[i][j] = 0;
	EnergyBack_Pulser[i][j] = 0;
	EnergyUp_Pulser[i][j] = 0;
	EnergyDown_Pulser[i][j] = 0;
	EnergyBack_Rel[i][j] = 0;
	EnergyUp_Rel[i][j] = 0;
	EnergyDown_Rel[i][j] = 0;
	EnergyBack_Cal[i][j] = 0;
	EnergyUp_Cal[i][j] = 0;
	EnergyDown_Cal[i][j] = 0;
	TimeBack[i][j] = 0;
	TimeUp[i][j] = 0;
	TimeDown[i][j] = 0;
      }
    }
  };
};
//////////////////////////////////////////////////////////////////////////////////////
class PCHit {
  // This class is for Proportional Counter hit
public:
  PCHit(){};
  Int_t NPCHits,WireID[MaxPCHits];
  Double_t Down[MaxPCHits],Up[MaxPCHits];
  Double_t DownVoltage[MaxPCHits], UpVoltage[MaxPCHits];
  Double_t Zp[MaxPCHits];
  Double_t XWp[MaxPCHits],YWp[MaxPCHits],ZWp[MaxPCHits];
  Double_t RWp[MaxPCHits], PhiWp[MaxPCHits];

  void zeroPCHit(){
    NPCHits = 0;
    for (Int_t i=0; i<MaxPCHits; i++) {
      WireID[i] = -1;
      Down[i] = 0;
      Up[i]   = 0;
      DownVoltage[i] = 0;
      UpVoltage[i] = 0;
      Zp[i] = 0;
      XWp[i] = 0;
      YWp[i] = 0;
      ZWp[i] = 0;
      RWp[i] = 0;
      PhiWp[i] = 0;
    }
  };
};

class SETTINGS {
public:
  SETTINGS(){};

  struct set_values {
    Double_t X3Pulser_Offset[24][12];
    Double_t X3Pulser_Slope[24][12];
    Double_t X3RelativeGains[24][12];
    Double_t X3FinalFix[24][12];
    Double_t X3Geometry_Up[24][4][4];
    Double_t X3Geometry_Down[24][4][4];
    Double_t QQQPulser_Offset[4][32];
    Double_t QQQPulser_Slope[4][32];
    Double_t QQQRelativeGains[4][32];
    Double_t QQQFinalFix[4][32];
    Double_t SiAlphaCal[28];

  }set_val;


};
