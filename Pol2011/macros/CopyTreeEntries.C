#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"
#include "TLorentzVector.h"

void CopyTreeEntries(Int_t iRapBin = 1, Int_t iPTBin = 1,
		     Char_t *fileNameIn = "RootFiles/selEvents_data_Ups_2Aug2011.root"){

  Char_t name[100];  
  Char_t fileNameOut[100];
  sprintf(fileNameOut, "RootFiles/data_Ups_rap%d_pT%d.root", iPTBin, iRapBin);
  printf("storing events in file %s\n", fileNameOut);

  TFile *fin = new TFile(fileNameIn);
  TLorentzVector *lepP;
  TLorentzVector *lepN;
  TTree *treeIn = (TTree *) gDirectory->Get("selectedData");

  TFile *fOut = new TFile(fileNameOut, "UPDATE");
  TTree *treeOut = treeIn->CloneTree(0);

  lepP = 0; lepN = 0;
  treeIn->SetBranchAddress("lepP", &lepP);
  treeIn->SetBranchAddress("lepN", &lepN);
  TLorentzVector *onia = new TLorentzVector();
  for(int iEn = 0; iEn < treeIn->GetEntries(); iEn++){ 
    Long64_t iEntry = treeIn->LoadTree(iEn);
    treeIn->GetEntry(iEntry);
    if(iEn % 100000 == 0)
      printf("entry %d out of %d\n", iEntry, treeIn->GetEntries());
    *onia = *(lepP) + *(lepN);
    if(onia->Pt() > onia::pTRange[iRapBin][iPTBin-1] && onia->Pt() < onia::pTRange[iRapBin][iPTBin] &&
       TMath::Abs(onia->Rapidity()) > onia::rapForPTRange[iRapBin-1] && TMath::Abs(onia->Rapidity()) < onia::rapForPTRange[iRapBin]){
      //printf("pT %1.3f, |y| %1.3f\n", onia->Pt(), TMath::Abs(onia->Rapidity()));
      treeOut->Fill();
    }
  }
  fOut->cd();
  treeOut->Write();
  fOut->Close();
}
