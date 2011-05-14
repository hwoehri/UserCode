#include "../interface/rootIncludes.inc"
#include "DataTnPTriggerEff.C"
#include "TChain.h"

void LoadEfficiencies(Int_t iEff, Int_t iEffSample);
void BookHistos(Char_t *oniaLabel);
void WriteHistos();
//======================================
void runDataTnPTriggerEff(Char_t *fileNameOut = "DataTnPTriggerEff_HLTDoubleMu0_RunA_14May2011.root",
			  Char_t *trigLabel = "HLT_DoubleMu0", //"HLT_DoubleMu0", "HLT_Mu0_TkMu0_OST_Jpsi" ... universal name for all of its instances
			  Int_t effSample = DATA, //DATA, MC, MCTRUTH
			  Char_t *oniaLabel = "J/#psi" //"J/#psi", "#psi'", "Ups(1S)", "Ups(2S)", "Ups(3S)"
		 ){

  // TFile *fIn = new TFile(fileNameIn);
  // TTree *treeData = (TTree*)fIn->Get("data");
  TChain *treeData = new TChain("data");
  treeData->Add("/Users/hwoehri/CMS/Work/TriggerEfficiency/OrthogonalPDs2010/Ilse/12April2011/Jet_Nov4th_merged_OniaMuMuv7.root");
  treeData->Add("/Users/hwoehri/CMS/Work/TriggerEfficiency/OrthogonalPDs2010/Ilse/12April2011/JetMET_Nov4th_merged_OniaMuMuv7.root");
  treeData->Add("/Users/hwoehri/CMS/Work/TriggerEfficiency/OrthogonalPDs2010/Roberto/10May2011/TTree_MuOnia_MultiJet2010.root");
  treeData->Add("/Users/hwoehri/CMS/Work/TriggerEfficiency/OrthogonalPDs2010/Roberto/10May2011/TTree_MuOnia_BTau2010.root");

  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  printf("initializing tree\n");
  DataTnPTriggerEff tree(treeData);  printf("...done\n");


  for(int iEff = 0; iEff < kNbEff; iEff++)
    LoadEfficiencies(iEff, effSample);
  printf("efficiencies loaded\n");

  BookHistos(oniaLabel);

  tree.Loop(effSample, trigLabel);

  fOut->cd();
  printf("writing out the histos\n");
  WriteHistos();
  fOut->Close();
}

//==============================================================
void LoadEfficiencies(Int_t iEff, Int_t iEffSample){

  TFile *fIn = new TFile(effFileNames[iEff]);
  Char_t name[100];
  //store the central value of the efficiencies
  sprintf(name, "hEff_%s_central", effSampleName[iEffSample]);
  hMuEff[iEff][iEffSample][CENTRAL] = (TH2D *) gDirectory->Get(name);
  sprintf(name, "h%s_%s", effName[iEff], effSampleName[iEffSample]);
  hMuEff[iEff][iEffSample][CENTRAL]->SetName(name);
  printf("%s, histo %p\n", effName[iEff], hMuEff[iEff][iEffSample][CENTRAL]->GetName());

  //load also the pT differential efficiencies from a fit
  if(iEff == L1L2Eff || iEff == L3Eff){
    for(int iEta = 0; iEta < kNbEtaBins; iEta++){
      sprintf(name, "fit%s_%s_pt_etaBin%d", effName[iEff], effSampleName[iEffSample], iEta);
      fMuEff_pT[iEff][iEffSample][iEta] = (TF1 *) gDirectory->Get(name);
    }
  }
}

//==========================================
void BookHistos(Char_t *oniaLabel){

  Char_t name[100], title[100];
  Int_t nMass = 40;
  Double_t massMin = 2.6, massMax = 3.4; //GeV
  //histos taking together +y and -y:
  for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < kNbPTBins[iRapBin]+1; iPTBin++){
      //RECO events
      sprintf(name, "reco_M_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";M [GeV]");
      reco_M[iPTBin][iRapBin] = new TH1D(name, title, nMass, massMin, massMax);
      reco_M[iPTBin][iRapBin]->Sumw2();
      //histo for passing trigger flag
      sprintf(name, "passedTrigger_M_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";M [GeV]");
      passedTrigger_M[iPTBin][iRapBin] = new TH1D(name, title, nMass, massMin, massMax);
      passedTrigger_M[iPTBin][iRapBin]->Sumw2();
      //histo for failing trigger flag
      sprintf(name, "failedTrigger_M_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";M [GeV]");
      failedTrigger_M[iPTBin][iRapBin] = new TH1D(name, title, nMass, massMin, massMax);
      failedTrigger_M[iPTBin][iRapBin]->Sumw2();
    }
  }
}

//==========================================
void WriteHistos(){

  for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < kNbPTBins[iRapBin]+1; iPTBin++){
      reco_M[iPTBin][iRapBin]->Write();
      passedTrigger_M[iPTBin][iRapBin]->Write();
      failedTrigger_M[iPTBin][iRapBin]->Write();
    }
  }
}
