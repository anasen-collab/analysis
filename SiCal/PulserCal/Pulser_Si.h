//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  3 13:51:08 2012 by ROOT version 5.34/00
// from TTree DataTree/DataTree
// found on file: /data0/ANASEN08_Oxygen/ROOT_Files/CalibrationData/Pulser_Si_r1327.root
//////////////////////////////////////////////////////////

#ifndef Pulser_Si_h
#define Pulser_Si_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Pulser_Si : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           Si_Nhits;
   Int_t           Si_MBID[360];   //[SiNhits]
   Int_t           Si_CBID[360];   //[SiNhits]
   Int_t           Si_ChNum[360];   //[SiNhits]
   Int_t           Si_Energy[360];   //[SiNhits]
   Int_t           Si_Time[360];   //[SiNhits]
   Int_t           ADC_Nhits;
   Int_t           ADC_ID[116];   //[ADCNhits]
   Int_t           ADC_ChNum[116];   //[ADCNhits]
   Int_t           ADC_Data[116];   //[ADCNhits]
   Int_t           TDC_Nhits;
   Int_t           TDC_ID[20];   //[TDCNhits]
   Int_t           TDC_ChNum[20];   //[TDCNhits]
   Int_t           TDC_Data[20];   //[TDCNhits]

   // List of branches
   TBranch        *b_SiNhits;   //!
   TBranch        *b_Si_MBID;   //!
   TBranch        *b_Si_CBID;   //!
   TBranch        *b_Si_ChNum;   //!
   TBranch        *b_Si_Energy;   //!
   TBranch        *b_Si_Time;   //!
   TBranch        *b_ADCNhits;   //!
   TBranch        *b_ADC_ID;   //!
   TBranch        *b_ADC_ChNum;   //!
   TBranch        *b_ADC_Data;   //!
   TBranch        *b_TDCNhits;   //!
   TBranch        *b_TDC_ID;   //!
   TBranch        *b_TDC_ChNum;   //!
   TBranch        *b_TDC_Data;   //!

   Pulser_Si(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Pulser_Si() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Pulser_Si,0);
};

#endif

#ifdef Pulser_Si_cxx
void Pulser_Si::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Si.Nhits", &Si_Nhits, &b_SiNhits);
   fChain->SetBranchAddress("Si.MBID", Si_MBID, &b_Si_MBID);
   fChain->SetBranchAddress("Si.CBID", Si_CBID, &b_Si_CBID);
   fChain->SetBranchAddress("Si.ChNum", Si_ChNum, &b_Si_ChNum);
   fChain->SetBranchAddress("Si.Energy", Si_Energy, &b_Si_Energy);
   fChain->SetBranchAddress("Si.Time", Si_Time, &b_Si_Time);
   fChain->SetBranchAddress("ADC.Nhits", &ADC_Nhits, &b_ADCNhits);
   fChain->SetBranchAddress("ADC.ID", ADC_ID, &b_ADC_ID);
   fChain->SetBranchAddress("ADC.ChNum", ADC_ChNum, &b_ADC_ChNum);
   fChain->SetBranchAddress("ADC.Data", ADC_Data, &b_ADC_Data);
   fChain->SetBranchAddress("TDC.Nhits", &TDC_Nhits, &b_TDCNhits);
   fChain->SetBranchAddress("TDC.ID", TDC_ID, &b_TDC_ID);
   fChain->SetBranchAddress("TDC.ChNum", TDC_ChNum, &b_TDC_ChNum);
   fChain->SetBranchAddress("TDC.Data", TDC_Data, &b_TDC_Data);
}

Bool_t Pulser_Si::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Pulser_Si_cxx
