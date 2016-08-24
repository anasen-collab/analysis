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

  vector<vector<Double_t> > Energy;
  vector<vector<Double_t> > Time;
  vector<vector<Double_t> > X, Y, Z;

  vector<vector<Double_t> > XW, YW, ZW;
  vector<vector<Double_t> > RW, PhiW;
  vector<vector<Double_t> > RFSubtract;

  void zeroSiHit()
  {
    NSiHits = 0;
    NHitsInDet.clear();
    DetID.clear();

    UpMult.clear();
    DownMult.clear();
    BackMult.clear();

    HitType.clear();
    UpChNum.clear();
    DownChNum.clear();
    BackChNum.clear();

    EnergyUp_Raw.clear();
    EnergyDown_Raw.clear();
    EnergyBack_Raw.clear();
    EnergyUp_Pulser.clear();
    EnergyDown_Pulser.clear();
    EnergyBack_Pulser.clear();
    EnergyUp_Rel.clear();
    EnergyDown_Rel.clear();
    EnergyBack_Rel.clear();
    EnergyUp_Cal.clear();
    EnergyDown_Cal.clear();
    EnergyBack_Cal.clear();
    
    TimeUp.clear();
    TimeDown.clear();
    TimeBack.clear();
    
    Energy.clear();
    Time.clear();
    X.clear();
    Y.clear();
    Z.clear();
    
    XW.clear();
    YW.clear();
    ZW.clear();
    RW.clear();
    PhiW.clear();
    RFSubtract.clear();
    
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
