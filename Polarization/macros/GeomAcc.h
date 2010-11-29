//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 29 14:53:03 2010 by ROOT version 5.26/00
// from TTree QQbarGenTree/generator level info for acceptance calculation
// found on file: jpsiGun_Tree.root
//////////////////////////////////////////////////////////

#ifndef GeomAcc_h
#define GeomAcc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class GeomAcc {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        onia_Px;
   Double_t        onia_Py;
   Double_t        onia_Pz;
   Double_t        muPos_Px;
   Double_t        muPos_Py;
   Double_t        muPos_Pz;
   Double_t        muNeg_Px;
   Double_t        muNeg_Py;
   Double_t        muNeg_Pz;

   // List of branches
   TBranch        *b_onia_Px;   //!
   TBranch        *b_onia_Py;   //!
   TBranch        *b_onia_Pz;   //!
   TBranch        *b_muPos_Px;   //!
   TBranch        *b_muPos_Py;   //!
   TBranch        *b_muPos_Pz;   //!
   TBranch        *b_muNeg_Px;   //!
   TBranch        *b_muNeg_Py;   //!
   TBranch        *b_muNeg_Pz;   //!

   GeomAcc(TTree *tree=0);
   virtual ~GeomAcc();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GeomAcc_cxx
GeomAcc::GeomAcc(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("jpsiGun_Tree.root");
      if (!f) {
         f = new TFile("jpsiGun_Tree.root");
      }
      tree = (TTree*)gDirectory->Get("QQbarGenTree");

   }
   Init(tree);
}

GeomAcc::~GeomAcc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GeomAcc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GeomAcc::LoadTree(Long64_t entry)
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

void GeomAcc::Init(TTree *tree)
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

   fChain->SetBranchAddress("onia_Px", &onia_Px, &b_onia_Px);
   fChain->SetBranchAddress("onia_Py", &onia_Py, &b_onia_Py);
   fChain->SetBranchAddress("onia_Pz", &onia_Pz, &b_onia_Pz);
   fChain->SetBranchAddress("muPos_Px", &muPos_Px, &b_muPos_Px);
   fChain->SetBranchAddress("muPos_Py", &muPos_Py, &b_muPos_Py);
   fChain->SetBranchAddress("muPos_Pz", &muPos_Pz, &b_muPos_Pz);
   fChain->SetBranchAddress("muNeg_Px", &muNeg_Px, &b_muNeg_Px);
   fChain->SetBranchAddress("muNeg_Py", &muNeg_Py, &b_muNeg_Py);
   fChain->SetBranchAddress("muNeg_Pz", &muNeg_Pz, &b_muNeg_Pz);
   Notify();
}

Bool_t GeomAcc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GeomAcc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GeomAcc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GeomAcc_cxx
