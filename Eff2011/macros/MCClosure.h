//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec  7 22:11:40 2010 by ROOT version 5.26/00
// from TTree data/CMSSW Quarkonia J/psi Polarization+Trigger Tree
// found on file: TTree_pol_noTriggerFilter_Run2010A-Nov4ReReco_v1_06Dec2010.root
//////////////////////////////////////////////////////////

#ifndef MCClosure_h
#define MCClosure_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"

class MCClosure {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   TLorentzVector  *onia, *muPos, *muNeg;
   TLorentzVector  *onia_Gen, *muPos_Gen, *muNeg_Gen;

   // Declaration of leaf types
   Int_t           eventNb;
   Int_t           runNb;
   Int_t           lumiBlock;
   Int_t           nPriVtx;
   Double_t        vertexWeight;
   Double_t        sumPTPV;
   Int_t           countTksOfPV;
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
   Int_t           HLT_Dimuon5_Upsilon_Barrel_v3;
   Int_t           HLT_Dimuon5_Upsilon_Barrel_v3_PreScale;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v3;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v3_PreScale;
   Int_t           HLT_Dimuon7_Jpsi_X_Barrel_v3;
   Int_t           HLT_Dimuon7_Jpsi_X_Barrel_v3_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v1;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v1_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v1;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v1_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_v3;
   Int_t           HLT_Dimuon0_Jpsi_v3_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_v3;
   Int_t           HLT_Dimuon0_Upsilon_v3_PreScale;
   Int_t           HLT_Dimuon5_Upsilon_Barrel_v5;
   Int_t           HLT_Dimuon5_Upsilon_Barrel_v5_PreScale;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v5;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v5_PreScale;
   Int_t           HLT_Dimuon7_Jpsi_X_Barrel_v5;
   Int_t           HLT_Dimuon7_Jpsi_X_Barrel_v5_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v6;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v6_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v6;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v6_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_v5;
   Int_t           HLT_Dimuon0_Jpsi_v5_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_v5;
   Int_t           HLT_Dimuon0_Upsilon_v5_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_NoVertexing_v2;
   Int_t           HLT_Dimuon0_Jpsi_NoVertexing_v2_PreScale;
   Int_t           HLT_Dimuon7_Upsilon_Barrel_v1;
   Int_t           HLT_Dimuon7_Upsilon_Barrel_v1_PreScale;
   Int_t           HLT_Dimuon9_Upsilon_Barrel_v1;
   Int_t           HLT_Dimuon9_Upsilon_Barrel_v1_PreScale;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v6;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v6_PreScale;
   Int_t           HLT_Dimuon13_Jpsi_Barrel_v1;
   Int_t           HLT_Dimuon13_Jpsi_Barrel_v1_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v7;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v7_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v7;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v7_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_v6;
   Int_t           HLT_Dimuon0_Jpsi_v6_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_v6;
   Int_t           HLT_Dimuon0_Upsilon_v6_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_NoVertexing_v3;
   Int_t           HLT_Dimuon0_Jpsi_NoVertexing_v3_PreScale;
   Int_t           HLT_Dimuon7_PsiPrime_v3;
   Int_t           HLT_Dimuon7_PsiPrime_v3_PreScale;
   Int_t           HLT_Dimuon7_PsiPrime_v5;
   Int_t           HLT_Dimuon7_PsiPrime_v5_PreScale;
   Int_t           HLT_Dimuon9_PsiPrime_v1;
   Int_t           HLT_Dimuon9_PsiPrime_v1_PreScale;
   Int_t           HLT_Dimuon11_PsiPrime_v1;
   Int_t           HLT_Dimuon11_PsiPrime_v1_PreScale;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v8;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v8_PreScale;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v9;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v9_PreScale;
   Int_t           HLT_Mu5_Track2_Jpsi_v9;
   Int_t           HLT_Mu5_Track2_Jpsi_v9_PreScale;
   Int_t           HLT_Mu7_Track7_Jpsi_v10;
   Int_t           HLT_Mu7_Track7_Jpsi_v10_PreScale;
   Int_t           MCType;
   Double_t        Jpsict_Gen;

   // List of branches
   TBranch        *b_eventNb;   //!
   TBranch        *b_runNb;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_nPriVtx;   //!
   TBranch        *b_vertexWeight;   //!
   TBranch        *b_sumPTPV;   //!
   TBranch        *b_countTksOfPV;   //!
   TBranch        *b_JpsiType;   //!
   TBranch        *b_JpsiCharge;   //!
   TBranch        *b_Jpsict;   //!
   TBranch        *b_JpsictErr;   //!
   TBranch        *b_JpsiVprob;   //!
   TBranch        *b_JpsiDistM1;   //!
   TBranch        *b_JpsiDphiM1;   //!
   TBranch        *b_JpsiDrM1;   //!
   TBranch        *b_JpsiDistM2;   //!
   TBranch        *b_JpsiDphiM2;   //!
   TBranch        *b_JpsiDrM2;   //!
   TBranch        *b_HLT_Dimuon5_Upsilon_Barrel_v3;   //!
   TBranch        *b_HLT_Dimuon5_Upsilon_Barrel_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v3;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_X_Barrel_v3;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_X_Barrel_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v1;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v1;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v3;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v3;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon5_Upsilon_Barrel_v5;   //!
   TBranch        *b_HLT_Dimuon5_Upsilon_Barrel_v5_PreScale;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v5;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v5_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_X_Barrel_v5;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_X_Barrel_v5_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v6;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v6_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v6;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v6_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v5;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v5_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v5;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v5_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing_v2;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_Upsilon_Barrel_v1;   //!
   TBranch        *b_HLT_Dimuon7_Upsilon_Barrel_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon9_Upsilon_Barrel_v1;   //!
   TBranch        *b_HLT_Dimuon9_Upsilon_Barrel_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v6;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v6_PreScale;   //!
   TBranch        *b_HLT_Dimuon13_Jpsi_Barrel_v1;   //!
   TBranch        *b_HLT_Dimuon13_Jpsi_Barrel_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v7;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v7_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v7;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v7_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v6;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v6_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v6;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v6_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing_v3;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_PsiPrime_v3;   //!
   TBranch        *b_HLT_Dimuon7_PsiPrime_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_PsiPrime_v5;   //!
   TBranch        *b_HLT_Dimuon7_PsiPrime_v5_PreScale;   //!
   TBranch        *b_HLT_Dimuon9_PsiPrime_v1;   //!
   TBranch        *b_HLT_Dimuon9_PsiPrime_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon11_PsiPrime_v1;   //!
   TBranch        *b_HLT_Dimuon11_PsiPrime_v1_PreScale;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v8;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v8_PreScale;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v9;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v9_PreScale;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v9;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v9_PreScale;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v10;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v10_PreScale;   //!
   TBranch        *b_MCType;   //!
   TBranch        *b_Jpsict_Gen;   //!


   MCClosure(TTree *tree=0);
   virtual ~MCClosure();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t effSample, Char_t *trigLabel, Bool_t rejectCowboys, Bool_t use2DGraph);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MCClosure_cxx
MCClosure::MCClosure(TTree *tree)
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

MCClosure::~MCClosure()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MCClosure::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MCClosure::LoadTree(Long64_t entry)
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

void MCClosure::Init(TTree *tree)
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
   fChain->SetBranchAddress("vertexWeight", &vertexWeight, &b_vertexWeight);
   fChain->SetBranchAddress("sumPTPV", &sumPTPV, &b_sumPTPV);
   fChain->SetBranchAddress("countTksOfPV", &countTksOfPV, &b_countTksOfPV);
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
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v3", &HLT_Dimuon5_Upsilon_Barrel_v3, &b_HLT_Dimuon5_Upsilon_Barrel_v3);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v3_PreScale", &HLT_Dimuon5_Upsilon_Barrel_v3_PreScale, &b_HLT_Dimuon5_Upsilon_Barrel_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v3", &HLT_Dimuon10_Jpsi_Barrel_v3, &b_HLT_Dimuon10_Jpsi_Barrel_v3);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v3_PreScale", &HLT_Dimuon10_Jpsi_Barrel_v3_PreScale, &b_HLT_Dimuon10_Jpsi_Barrel_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel_v3", &HLT_Dimuon7_Jpsi_X_Barrel_v3, &b_HLT_Dimuon7_Jpsi_X_Barrel_v3);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel_v3_PreScale", &HLT_Dimuon7_Jpsi_X_Barrel_v3_PreScale, &b_HLT_Dimuon7_Jpsi_X_Barrel_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v1", &HLT_Dimuon0_Jpsi_Muon_v1, &b_HLT_Dimuon0_Jpsi_Muon_v1);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v1_PreScale", &HLT_Dimuon0_Jpsi_Muon_v1_PreScale, &b_HLT_Dimuon0_Jpsi_Muon_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v1", &HLT_Dimuon0_Upsilon_Muon_v1, &b_HLT_Dimuon0_Upsilon_Muon_v1);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v1_PreScale", &HLT_Dimuon0_Upsilon_Muon_v1_PreScale, &b_HLT_Dimuon0_Upsilon_Muon_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v3", &HLT_Dimuon0_Jpsi_v3, &b_HLT_Dimuon0_Jpsi_v3);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v3_PreScale", &HLT_Dimuon0_Jpsi_v3_PreScale, &b_HLT_Dimuon0_Jpsi_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v3", &HLT_Dimuon0_Upsilon_v3, &b_HLT_Dimuon0_Upsilon_v3);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v3_PreScale", &HLT_Dimuon0_Upsilon_v3_PreScale, &b_HLT_Dimuon0_Upsilon_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v5", &HLT_Dimuon5_Upsilon_Barrel_v5, &b_HLT_Dimuon5_Upsilon_Barrel_v5);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v5_PreScale", &HLT_Dimuon5_Upsilon_Barrel_v5_PreScale, &b_HLT_Dimuon5_Upsilon_Barrel_v5_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v5", &HLT_Dimuon10_Jpsi_Barrel_v5, &b_HLT_Dimuon10_Jpsi_Barrel_v5);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v5_PreScale", &HLT_Dimuon10_Jpsi_Barrel_v5_PreScale, &b_HLT_Dimuon10_Jpsi_Barrel_v5_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel_v5", &HLT_Dimuon7_Jpsi_X_Barrel_v5, &b_HLT_Dimuon7_Jpsi_X_Barrel_v5);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel_v5_PreScale", &HLT_Dimuon7_Jpsi_X_Barrel_v5_PreScale, &b_HLT_Dimuon7_Jpsi_X_Barrel_v5_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v6", &HLT_Dimuon0_Jpsi_Muon_v6, &b_HLT_Dimuon0_Jpsi_Muon_v6);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v6_PreScale", &HLT_Dimuon0_Jpsi_Muon_v6_PreScale, &b_HLT_Dimuon0_Jpsi_Muon_v6_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v6", &HLT_Dimuon0_Upsilon_Muon_v6, &b_HLT_Dimuon0_Upsilon_Muon_v6);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v6_PreScale", &HLT_Dimuon0_Upsilon_Muon_v6_PreScale, &b_HLT_Dimuon0_Upsilon_Muon_v6_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v5", &HLT_Dimuon0_Jpsi_v5, &b_HLT_Dimuon0_Jpsi_v5);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v5_PreScale", &HLT_Dimuon0_Jpsi_v5_PreScale, &b_HLT_Dimuon0_Jpsi_v5_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v5", &HLT_Dimuon0_Upsilon_v5, &b_HLT_Dimuon0_Upsilon_v5);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v5_PreScale", &HLT_Dimuon0_Upsilon_v5_PreScale, &b_HLT_Dimuon0_Upsilon_v5_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_v2", &HLT_Dimuon0_Jpsi_NoVertexing_v2, &b_HLT_Dimuon0_Jpsi_NoVertexing_v2);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_v2_PreScale", &HLT_Dimuon0_Jpsi_NoVertexing_v2_PreScale, &b_HLT_Dimuon0_Jpsi_NoVertexing_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_Upsilon_Barrel_v1", &HLT_Dimuon7_Upsilon_Barrel_v1, &b_HLT_Dimuon7_Upsilon_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon7_Upsilon_Barrel_v1_PreScale", &HLT_Dimuon7_Upsilon_Barrel_v1_PreScale, &b_HLT_Dimuon7_Upsilon_Barrel_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon9_Upsilon_Barrel_v1", &HLT_Dimuon9_Upsilon_Barrel_v1, &b_HLT_Dimuon9_Upsilon_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon9_Upsilon_Barrel_v1_PreScale", &HLT_Dimuon9_Upsilon_Barrel_v1_PreScale, &b_HLT_Dimuon9_Upsilon_Barrel_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v6", &HLT_Dimuon10_Jpsi_Barrel_v6, &b_HLT_Dimuon10_Jpsi_Barrel_v6);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v6_PreScale", &HLT_Dimuon10_Jpsi_Barrel_v6_PreScale, &b_HLT_Dimuon10_Jpsi_Barrel_v6_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon13_Jpsi_Barrel_v1", &HLT_Dimuon13_Jpsi_Barrel_v1, &b_HLT_Dimuon13_Jpsi_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon13_Jpsi_Barrel_v1_PreScale", &HLT_Dimuon13_Jpsi_Barrel_v1_PreScale, &b_HLT_Dimuon13_Jpsi_Barrel_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v7", &HLT_Dimuon0_Jpsi_Muon_v7, &b_HLT_Dimuon0_Jpsi_Muon_v7);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v7_PreScale", &HLT_Dimuon0_Jpsi_Muon_v7_PreScale, &b_HLT_Dimuon0_Jpsi_Muon_v7_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v7", &HLT_Dimuon0_Upsilon_Muon_v7, &b_HLT_Dimuon0_Upsilon_Muon_v7);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v7_PreScale", &HLT_Dimuon0_Upsilon_Muon_v7_PreScale, &b_HLT_Dimuon0_Upsilon_Muon_v7_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v6", &HLT_Dimuon0_Jpsi_v6, &b_HLT_Dimuon0_Jpsi_v6);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v6_PreScale", &HLT_Dimuon0_Jpsi_v6_PreScale, &b_HLT_Dimuon0_Jpsi_v6_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v6", &HLT_Dimuon0_Upsilon_v6, &b_HLT_Dimuon0_Upsilon_v6);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v6_PreScale", &HLT_Dimuon0_Upsilon_v6_PreScale, &b_HLT_Dimuon0_Upsilon_v6_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_v3", &HLT_Dimuon0_Jpsi_NoVertexing_v3, &b_HLT_Dimuon0_Jpsi_NoVertexing_v3);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_v3_PreScale", &HLT_Dimuon0_Jpsi_NoVertexing_v3_PreScale, &b_HLT_Dimuon0_Jpsi_NoVertexing_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v3", &HLT_Dimuon7_PsiPrime_v3, &b_HLT_Dimuon7_PsiPrime_v3);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v3_PreScale", &HLT_Dimuon7_PsiPrime_v3_PreScale, &b_HLT_Dimuon7_PsiPrime_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v5", &HLT_Dimuon7_PsiPrime_v5, &b_HLT_Dimuon7_PsiPrime_v5);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v5_PreScale", &HLT_Dimuon7_PsiPrime_v5_PreScale, &b_HLT_Dimuon7_PsiPrime_v5_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon9_PsiPrime_v1", &HLT_Dimuon9_PsiPrime_v1, &b_HLT_Dimuon9_PsiPrime_v1);
   fChain->SetBranchAddress("HLT_Dimuon9_PsiPrime_v1_PreScale", &HLT_Dimuon9_PsiPrime_v1_PreScale, &b_HLT_Dimuon9_PsiPrime_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon11_PsiPrime_v1", &HLT_Dimuon11_PsiPrime_v1, &b_HLT_Dimuon11_PsiPrime_v1);
   fChain->SetBranchAddress("HLT_Dimuon11_PsiPrime_v1_PreScale", &HLT_Dimuon11_PsiPrime_v1_PreScale, &b_HLT_Dimuon11_PsiPrime_v1_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v8", &HLT_Mu5_L2Mu2_Jpsi_v8, &b_HLT_Mu5_L2Mu2_Jpsi_v8);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v8_PreScale", &HLT_Mu5_L2Mu2_Jpsi_v8_PreScale, &b_HLT_Mu5_L2Mu2_Jpsi_v8_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v9", &HLT_Mu5_L2Mu2_Jpsi_v9, &b_HLT_Mu5_L2Mu2_Jpsi_v9);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v9_PreScale", &HLT_Mu5_L2Mu2_Jpsi_v9_PreScale, &b_HLT_Mu5_L2Mu2_Jpsi_v9_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v9", &HLT_Mu5_Track2_Jpsi_v9, &b_HLT_Mu5_Track2_Jpsi_v9);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v9_PreScale", &HLT_Mu5_Track2_Jpsi_v9_PreScale, &b_HLT_Mu5_Track2_Jpsi_v9_PreScale);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v10", &HLT_Mu7_Track7_Jpsi_v10, &b_HLT_Mu7_Track7_Jpsi_v10);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v10_PreScale", &HLT_Mu7_Track7_Jpsi_v10_PreScale, &b_HLT_Mu7_Track7_Jpsi_v10_PreScale);
   fChain->SetBranchAddress("MCType", &MCType, &b_MCType);
   fChain->SetBranchAddress("Jpsict_Gen", &Jpsict_Gen, &b_Jpsict_Gen);

   Notify();
}

Bool_t MCClosure::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MCClosure::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MCClosure::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MCClosure_cxx
