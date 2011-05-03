//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec  7 22:11:40 2010 by ROOT version 5.26/00
// from TTree data/CMSSW Quarkonia J/psi Polarization+Trigger Tree
// found on file: TTree_pol_noTriggerFilter_Run2010A-Nov4ReReco_v1_06Dec2010.root
//////////////////////////////////////////////////////////

#ifndef MCTnPEff_h
#define MCTnPEff_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"

class MCTnPEff {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   TLorentzVector  *onia, *muPos, *muNeg;
   TLorentzVector  *onia_Gen, *muPos_Gen, *muNeg_Gen;
   Int_t           eventNb;
   Int_t           runNb;
   Int_t           lumiBlock;
   Int_t           nPriVtx;
   Int_t           JpsiType;
   Int_t           JpsiCharge;
   Double_t        Jpsict;
   Double_t        JpsictErr;
   Double_t        JpsiVprob;
   Double_t        JpsiDistM1;
   Double_t        JpsiDphiM1;
   Double_t        JpsiDrM1;
   Double_t        JpsiDistM2;
   Double_t        JpsiDphiM2;
   Double_t        JpsiDrM2;
   Int_t           HLT_Mu0_Track0_Jpsi;
   Int_t           HLT_Mu3_Track0_Jpsi;
   Int_t           HLT_Mu5_Track0_Jpsi;
   Int_t           HLT_Mu0_TkMu0_Jpsi;
   Int_t           HLT_Mu3_TkMu0_Jpsi;
   Int_t           HLT_Mu5_TkMu0_Jpsi;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi;
   Int_t           HLT_Mu3_TkMu0_OST_Jpsi;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1;
   Int_t           HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2;
   Int_t           HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3;
   Int_t           HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3;
   Int_t           HLT_Mu3_Track3_Jpsi;
   Int_t           HLT_Mu3_Track3_Jpsi_v2;
   Int_t           HLT_Mu3_Track5_Jpsi_v1;
   Int_t           HLT_Mu3_Track5_Jpsi_v2;
   Int_t           HLT_DoubleMu0;
   Int_t           HLT_DoubleMu0_Quarkonium_v1;
   Int_t           HLT_DoubleMu0_Quarkonium_LS_v1;
   Int_t           HLT_L1DoubleMuOpen;
   Int_t           HLT_L1DoubleMuOpen_Tight;
   Int_t           HLT_DoubleMu3;
   Int_t           HLT_Mu3;
   Int_t           HLT_Mu5;
   Int_t           HLT_Mu7;
   Int_t           HLT_Mu9;
   Int_t           HLT_Mu11;

   // List of branches
   TBranch        *b_eventNb;   //!
   TBranch        *b_runNb;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_nPriVtx;   //!
   TBranch        *b_JpsiType;   //!
   TBranch        *b_JpsiCharge;   //!
   TBranch        *b_Jpsict;   //!
   TBranch        *b_JpsictErr;   //!
   TBranch        *b_JpsiVprob;   //!
   TBranch        *b_JpsiDistM1;
   TBranch        *b_JpsiDphiM1;
   TBranch        *b_JpsiDrM1;
   TBranch        *b_JpsiDistM2;
   TBranch        *b_JpsiDphiM2;
   TBranch        *b_JpsiDrM2;
   TBranch        *b_HLT_Mu0_Track0_Jpsi;   //!
   TBranch        *b_HLT_Mu3_Track0_Jpsi;   //!
   TBranch        *b_HLT_Mu5_Track0_Jpsi;   //!
   TBranch        *b_HLT_Mu0_TkMu0_Jpsi;   //!
   TBranch        *b_HLT_Mu3_TkMu0_Jpsi;   //!
   TBranch        *b_HLT_Mu5_TkMu0_Jpsi;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi;   //!
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1;   //!
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2;   //!
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3;   //!
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi_v2;   //!
   TBranch        *b_HLT_Mu3_Track5_Jpsi_v1;   //!
   TBranch        *b_HLT_Mu3_Track5_Jpsi_v2;   //!
   TBranch        *b_HLT_DoubleMu0;   //!
   TBranch        *b_HLT_DoubleMu0_Quarkonium_v1;   //!
   TBranch        *b_HLT_DoubleMu0_Quarkonium_LS_v1;   //!
   TBranch        *b_HLT_L1DoubleMuOpen;   //!
   TBranch        *b_HLT_L1DoubleMuOpen_Tight;   //!
   TBranch        *b_HLT_DoubleMu3;   //!
   TBranch        *b_HLT_Mu3;   //!
   TBranch        *b_HLT_Mu5;   //!
   TBranch        *b_HLT_Mu7;   //!
   TBranch        *b_HLT_Mu9;   //!
   TBranch        *b_HLT_Mu11;   //!

   MCTnPEff(TTree *tree=0);
   virtual ~MCTnPEff();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t effSample, Char_t *trigLabel);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MCTnPEff_cxx
MCTnPEff::MCTnPEff(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TTree_pol_noTriggerFilter_Run2010A-Nov4ReReco_v1_06Dec2010.root");
      if (!f) {
         f = new TFile("TTree_pol_noTriggerFilter_Run2010A-Nov4ReReco_v1_06Dec2010.root");
      }
      tree = (TTree*)gDirectory->Get("data");

   }
   Init(tree);
}

MCTnPEff::~MCTnPEff()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MCTnPEff::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MCTnPEff::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MCTnPEff::Init(TTree *tree)
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
   fCurrent = -1;
   //fChain->SetMakeClass(1); //H: don't use this with TLorentzVectors!!!

   onia = 0;
   muNeg = 0;
   muPos = 0;
   onia_Gen = 0;
   muNeg_Gen = 0;
   muPos_Gen = 0;

   fChain->SetBranchAddress("JpsiP", &onia);
   fChain->SetBranchAddress("muNegP", &muNeg);
   fChain->SetBranchAddress("muPosP", &muPos);
   fChain->SetBranchAddress("JpsiP_Gen", &onia_Gen);
   fChain->SetBranchAddress("muNegP_Gen", &muNeg_Gen);
   fChain->SetBranchAddress("muPosP_Gen", &muPos_Gen);

   fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("nPriVtx", &nPriVtx, &b_nPriVtx);
   fChain->SetBranchAddress("JpsiType", &JpsiType, &b_JpsiType);
   fChain->SetBranchAddress("JpsiCharge", &JpsiCharge, &b_JpsiCharge);
   fChain->SetBranchAddress("Jpsict", &Jpsict, &b_Jpsict);
   fChain->SetBranchAddress("JpsictErr", &JpsictErr, &b_JpsictErr);
   fChain->SetBranchAddress("JpsiVprob", &JpsiVprob, &b_JpsiVprob);
   fChain->SetBranchAddress("JpsiDistM1", &JpsiDistM1, &b_JpsiDistM1);
   fChain->SetBranchAddress("JpsiDphiM1", &JpsiDphiM1, &b_JpsiDphiM1);
   fChain->SetBranchAddress("JpsiDrM1", &JpsiDrM1, &b_JpsiDrM1);
   fChain->SetBranchAddress("JpsiDistM2", &JpsiDistM2, &b_JpsiDistM2);
   fChain->SetBranchAddress("JpsiDphiM2", &JpsiDphiM2, &b_JpsiDphiM2);
   fChain->SetBranchAddress("JpsiDrM2", &JpsiDrM2, &b_JpsiDrM2);
   fChain->SetBranchAddress("HLT_Mu0_Track0_Jpsi", &HLT_Mu0_Track0_Jpsi, &b_HLT_Mu0_Track0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_Track0_Jpsi", &HLT_Mu3_Track0_Jpsi, &b_HLT_Mu3_Track0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu5_Track0_Jpsi", &HLT_Mu5_Track0_Jpsi, &b_HLT_Mu5_Track0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_Jpsi", &HLT_Mu0_TkMu0_Jpsi, &b_HLT_Mu0_TkMu0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_Jpsi", &HLT_Mu3_TkMu0_Jpsi, &b_HLT_Mu3_TkMu0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_Jpsi", &HLT_Mu5_TkMu0_Jpsi, &b_HLT_Mu5_TkMu0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi", &HLT_Mu0_TkMu0_OST_Jpsi, &b_HLT_Mu0_TkMu0_OST_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi", &HLT_Mu3_TkMu0_OST_Jpsi, &b_HLT_Mu3_TkMu0_OST_Jpsi);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi", &HLT_Mu5_TkMu0_OST_Jpsi, &b_HLT_Mu5_TkMu0_OST_Jpsi);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1", &HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1, &b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1", &HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1, &b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1", &HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1, &b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2", &HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2, &b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2", &HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2, &b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2", &HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2, &b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3", &HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3, &b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3", &HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3, &b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3", &HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3, &b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi", &HLT_Mu3_Track3_Jpsi, &b_HLT_Mu3_Track3_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_v2", &HLT_Mu3_Track3_Jpsi_v2, &b_HLT_Mu3_Track3_Jpsi_v2);
   fChain->SetBranchAddress("HLT_Mu3_Track5_Jpsi_v1", &HLT_Mu3_Track5_Jpsi_v1, &b_HLT_Mu3_Track5_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Mu3_Track5_Jpsi_v2", &HLT_Mu3_Track5_Jpsi_v2, &b_HLT_Mu3_Track5_Jpsi_v2);
   fChain->SetBranchAddress("HLT_DoubleMu0", &HLT_DoubleMu0, &b_HLT_DoubleMu0);
   fChain->SetBranchAddress("HLT_DoubleMu0_Quarkonium_v1", &HLT_DoubleMu0_Quarkonium_v1, &b_HLT_DoubleMu0_Quarkonium_v1);
   fChain->SetBranchAddress("HLT_DoubleMu0_Quarkonium_LS_v1", &HLT_DoubleMu0_Quarkonium_LS_v1, &b_HLT_DoubleMu0_Quarkonium_LS_v1);
   fChain->SetBranchAddress("HLT_L1DoubleMuOpen", &HLT_L1DoubleMuOpen, &b_HLT_L1DoubleMuOpen);
   fChain->SetBranchAddress("HLT_L1DoubleMuOpen_Tight", &HLT_L1DoubleMuOpen_Tight, &b_HLT_L1DoubleMuOpen_Tight);
   fChain->SetBranchAddress("HLT_DoubleMu3", &HLT_DoubleMu3, &b_HLT_DoubleMu3);
   fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
   fChain->SetBranchAddress("HLT_Mu5", &HLT_Mu5, &b_HLT_Mu5);
   fChain->SetBranchAddress("HLT_Mu7", &HLT_Mu7, &b_HLT_Mu7);
   fChain->SetBranchAddress("HLT_Mu9", &HLT_Mu9, &b_HLT_Mu9);
   fChain->SetBranchAddress("HLT_Mu11", &HLT_Mu11, &b_HLT_Mu11);
   Notify();
}

Bool_t MCTnPEff::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MCTnPEff::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MCTnPEff::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MCTnPEff_cxx
