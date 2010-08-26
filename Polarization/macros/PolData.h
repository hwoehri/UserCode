//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 28 14:54:21 2010 by ROOT version 5.26/00
// from TTree data/new data
// found on file: tree_promptJPsi_v4.root
//////////////////////////////////////////////////////////

#ifndef PolData_h
#define PolData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class PolData {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        JpsiMass;
   Double_t        JpsiPt;
   Double_t        JpsiRap;
   Double_t        JpsiPx;
   Double_t        JpsiPy;
   Double_t        JpsiPz;
   Double_t        Jpsict;
   Double_t        JpsictErr;
   Int_t           JpsiType_idx;
   Char_t          JpsiType_lbl[3];
   /* Int_t           JpsiPtType_idx; */
   /* Char_t          JpsiPtType_lbl[4]; */
   /* Int_t           JpsiRapType_idx; */
   /* Char_t          JpsiRapType_lbl[5]; */
   Double_t        muPosPx;
   Double_t        muPosPy;
   Double_t        muPosPz;
   Double_t        muNegPx;
   Double_t        muNegPy;
   Double_t        muNegPz;
   Double_t        JpsiMass_Gen;
   Double_t        JpsiPt_Gen;
   Double_t        JpsiRap_Gen;
   Double_t        JpsiPx_Gen;
   Double_t        JpsiPy_Gen;
   Double_t        JpsiPz_Gen;
   Double_t        muPosPx_Gen;
   Double_t        muPosPy_Gen;
   Double_t        muPosPz_Gen;
   Double_t        muNegPx_Gen;
   Double_t        muNegPy_Gen;
   Double_t        muNegPz_Gen;
   Double_t        eventNb;
   Double_t        runNb;
   Double_t        lumiBlock;
   /* Int_t           type_L1DoubleMuOpen_idx; */
   /* Char_t          type_L1DoubleMuOpen_lbl[9]; */
   /* Int_t           type_L2DoubleMu0_idx; */
   /* Char_t          type_L2DoubleMu0_lbl[9]; */
   /* Int_t           type_Mu0_Track0_Jpsi_idx; */
   /* Char_t          type_Mu0_Track0_Jpsi_lbl[9]; */
   /* Int_t           type_Mu3_Track0_Jpsi_idx; */
   /* Char_t          type_Mu3_Track0_Jpsi_lbl[9]; */
   /* Int_t           type_Mu5_Track0_Jpsi_idx; */
   /* Char_t          type_Mu5_Track0_Jpsi_lbl[9]; */
   /* Int_t           type_DoubleMu0_idx; */
   /* Char_t          type_DoubleMu0_lbl[9]; */
   /* Int_t           type_DoubleMu3_idx; */
   /* Char_t          type_DoubleMu3_lbl[9]; */
   Int_t           MCType_idx;
   Char_t          MCType_lbl[3];

   // List of branches
   TBranch        *b_JpsiMass;   //!
   TBranch        *b_JpsiPt;   //!
   TBranch        *b_JpsiRap;   //!
   TBranch        *b_JpsiPx;   //!
   TBranch        *b_JpsiPy;   //!
   TBranch        *b_JpsiPz;   //!
   TBranch        *b_Jpsict;   //!
   TBranch        *b_JpsictErr;   //!
   TBranch        *b_JpsiType_idx;   //!
   TBranch        *b_JpsiType_lbl;   //!
   /* TBranch        *b_JpsiPtType_idx;   //! */
   /* TBranch        *b_JpsiPtType_lbl;   //! */
   /* TBranch        *b_JpsiRapType_idx;   //! */
   /* TBranch        *b_JpsiRapType_lbl;   //! */
   TBranch        *b_muPosPx;   //!
   TBranch        *b_muPosPy;   //!
   TBranch        *b_muPosPz;   //!
   TBranch        *b_muNegPx;   //!
   TBranch        *b_muNegPy;   //!
   TBranch        *b_muNegPz;   //!
   TBranch        *b_JpsiMass_Gen;   //!
   TBranch        *b_JpsiPt_Gen;   //!
   TBranch        *b_JpsiRap_Gen;   //!
   TBranch        *b_JpsiPx_Gen;   //!
   TBranch        *b_JpsiPy_Gen;   //!
   TBranch        *b_JpsiPz_Gen;   //!
   TBranch        *b_muPosPx_Gen;   //!
   TBranch        *b_muPosPy_Gen;   //!
   TBranch        *b_muPosPz_Gen;   //!
   TBranch        *b_muNegPx_Gen;   //!
   TBranch        *b_muNegPy_Gen;   //!
   TBranch        *b_muNegPz_Gen;   //!
   TBranch        *b_eventNb;   //!
   TBranch        *b_runNb;   //!
   TBranch        *b_lumiBlock;   //!
   /* TBranch        *b_type_L1DoubleMuOpen_idx;   //! */
   /* TBranch        *b_type_L1DoubleMuOpen_lbl;   //! */
   /* TBranch        *b_type_L2DoubleMu0_idx;   //! */
   /* TBranch        *b_type_L2DoubleMu0_lbl;   //! */
   /* TBranch        *b_type_Mu0_Track0_Jpsi_idx;   //! */
   /* TBranch        *b_type_Mu0_Track0_Jpsi_lbl;   //! */
   /* TBranch        *b_type_Mu3_Track0_Jpsi_idx;   //! */
   /* TBranch        *b_type_Mu3_Track0_Jpsi_lbl;   //! */
   /* TBranch        *b_type_Mu5_Track0_Jpsi_idx;   //! */
   /* TBranch        *b_type_Mu5_Track0_Jpsi_lbl;   //! */
   /* TBranch        *b_type_DoubleMu0_idx;   //! */
   /* TBranch        *b_type_DoubleMu0_lbl;   //! */
   /* TBranch        *b_type_DoubleMu3_idx;   //! */
   /* TBranch        *b_type_DoubleMu3_lbl;   //! */
   TBranch        *b_MCType_idx;   //!
   TBranch        *b_MCType_lbl;   //!

   PolData(TTree *tree=0);
   virtual ~PolData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t selDimuType);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PolData_cxx
PolData::PolData(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/hermine/CMS/Work/Polarization/Florian/25Aug2010/Spring10_PromptJPsi_Histos_TEST_0_10.root");
      if (!f) {
         f = new TFile("/home/hermine/CMS/Work/Polarization/Florian/25Aug2010/Spring10_PromptJPsi_Histos_TEST_0_10.root");
      }
      tree = (TTree*)gDirectory->Get("data");

   }
   Init(tree);
}

PolData::~PolData()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PolData::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PolData::LoadTree(Long64_t entry)
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

void PolData::Init(TTree *tree)
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
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("JpsiMass", &JpsiMass, &b_JpsiMass);
   fChain->SetBranchAddress("JpsiPt", &JpsiPt, &b_JpsiPt);
   fChain->SetBranchAddress("JpsiRap", &JpsiRap, &b_JpsiRap);
   fChain->SetBranchAddress("JpsiPx", &JpsiPx, &b_JpsiPx);
   fChain->SetBranchAddress("JpsiPy", &JpsiPy, &b_JpsiPy);
   fChain->SetBranchAddress("JpsiPz", &JpsiPz, &b_JpsiPz);
   fChain->SetBranchAddress("Jpsict", &Jpsict, &b_Jpsict);
   fChain->SetBranchAddress("JpsictErr", &JpsictErr, &b_JpsictErr);
   fChain->SetBranchAddress("JpsiType_idx", &JpsiType_idx, &b_JpsiType_idx);
   fChain->SetBranchAddress("JpsiType_lbl", JpsiType_lbl, &b_JpsiType_lbl);
   /* fChain->SetBranchAddress("JpsiPtType_idx", &JpsiPtType_idx, &b_JpsiPtType_idx); */
   /* fChain->SetBranchAddress("JpsiPtType_lbl", JpsiPtType_lbl, &b_JpsiPtType_lbl); */
   /* fChain->SetBranchAddress("JpsiRapType_idx", &JpsiRapType_idx, &b_JpsiRapType_idx); */
   /* fChain->SetBranchAddress("JpsiRapType_lbl", JpsiRapType_lbl, &b_JpsiRapType_lbl); */
   fChain->SetBranchAddress("muPosPx", &muPosPx, &b_muPosPx);
   fChain->SetBranchAddress("muPosPy", &muPosPy, &b_muPosPy);
   fChain->SetBranchAddress("muPosPz", &muPosPz, &b_muPosPz);
   fChain->SetBranchAddress("muNegPx", &muNegPx, &b_muNegPx);
   fChain->SetBranchAddress("muNegPy", &muNegPy, &b_muNegPy);
   fChain->SetBranchAddress("muNegPz", &muNegPz, &b_muNegPz);
   fChain->SetBranchAddress("JpsiMass_Gen", &JpsiMass_Gen, &b_JpsiMass_Gen);
   fChain->SetBranchAddress("JpsiPt_Gen", &JpsiPt_Gen, &b_JpsiPt_Gen);
   fChain->SetBranchAddress("JpsiRap_Gen", &JpsiRap_Gen, &b_JpsiRap_Gen);
   fChain->SetBranchAddress("JpsiPx_Gen", &JpsiPx_Gen, &b_JpsiPx_Gen);
   fChain->SetBranchAddress("JpsiPy_Gen", &JpsiPy_Gen, &b_JpsiPy_Gen);
   fChain->SetBranchAddress("JpsiPz_Gen", &JpsiPz_Gen, &b_JpsiPz_Gen);
   fChain->SetBranchAddress("muPosPx_Gen", &muPosPx_Gen, &b_muPosPx_Gen);
   fChain->SetBranchAddress("muPosPy_Gen", &muPosPy_Gen, &b_muPosPy_Gen);
   fChain->SetBranchAddress("muPosPz_Gen", &muPosPz_Gen, &b_muPosPz_Gen);
   fChain->SetBranchAddress("muNegPx_Gen", &muNegPx_Gen, &b_muNegPx_Gen);
   fChain->SetBranchAddress("muNegPy_Gen", &muNegPy_Gen, &b_muNegPy_Gen);
   fChain->SetBranchAddress("muNegPz_Gen", &muNegPz_Gen, &b_muNegPz_Gen);
   /* fChain->SetBranchAddress("type_L1DoubleMuOpen_idx", &type_L1DoubleMuOpen_idx, &b_type_L1DoubleMuOpen_idx); */
   /* fChain->SetBranchAddress("type_L1DoubleMuOpen_lbl", type_L1DoubleMuOpen_lbl, &b_type_L1DoubleMuOpen_lbl); */
   /* fChain->SetBranchAddress("type_L2DoubleMu0_idx", &type_L2DoubleMu0_idx, &b_type_L2DoubleMu0_idx); */
   /* fChain->SetBranchAddress("type_L2DoubleMu0_lbl", type_L2DoubleMu0_lbl, &b_type_L2DoubleMu0_lbl); */
   /* fChain->SetBranchAddress("type_Mu0_Track0_Jpsi_idx", &type_Mu0_Track0_Jpsi_idx, &b_type_Mu0_Track0_Jpsi_idx); */
   /* fChain->SetBranchAddress("type_Mu0_Track0_Jpsi_lbl", type_Mu0_Track0_Jpsi_lbl, &b_type_Mu0_Track0_Jpsi_lbl); */
   /* fChain->SetBranchAddress("type_Mu3_Track0_Jpsi_idx", &type_Mu3_Track0_Jpsi_idx, &b_type_Mu3_Track0_Jpsi_idx); */
   /* fChain->SetBranchAddress("type_Mu3_Track0_Jpsi_lbl", type_Mu3_Track0_Jpsi_lbl, &b_type_Mu3_Track0_Jpsi_lbl); */
   /* fChain->SetBranchAddress("type_Mu5_Track0_Jpsi_idx", &type_Mu5_Track0_Jpsi_idx, &b_type_Mu5_Track0_Jpsi_idx); */
   /* fChain->SetBranchAddress("type_Mu5_Track0_Jpsi_lbl", type_Mu5_Track0_Jpsi_lbl, &b_type_Mu5_Track0_Jpsi_lbl); */
   /* fChain->SetBranchAddress("type_DoubleMu0_idx", &type_DoubleMu0_idx, &b_type_DoubleMu0_idx); */
   /* fChain->SetBranchAddress("type_DoubleMu0_lbl", type_DoubleMu0_lbl, &b_type_DoubleMu0_lbl); */
   /* fChain->SetBranchAddress("type_DoubleMu3_idx", &type_DoubleMu3_idx, &b_type_DoubleMu3_idx); */
   /* fChain->SetBranchAddress("type_DoubleMu3_lbl", type_DoubleMu3_lbl, &b_type_DoubleMu3_lbl); */
   fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("MCType_idx", &MCType_idx, &b_MCType_idx);
   fChain->SetBranchAddress("MCType_lbl", MCType_lbl, &b_MCType_lbl);
   Notify();
}

Bool_t PolData::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PolData::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PolData::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PolData_cxx
