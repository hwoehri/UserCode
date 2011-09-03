#ifndef __CINT__
#endif

#include "TChain.h"
#include "../interface/rootIncludes.inc"
#include "PolData.C"

void BookHistosReco();
void WriteHistosReco(Char_t *fNameOut);
//==========================================
void runData(Char_t *fNameOut = "RootFiles/selEvents_data_Ups_1Aug2011.root",
	     Bool_t rejectCowboys = kTRUE
	     ){
  
  TChain *chain = new TChain("data");
  chain->Add("/Users/hwoehri/CMS/Work/Data2011/May10ReReco/TTree_Onia2MuMu_V8_May10ReReco_v1.root");
  chain->Add("/Users/hwoehri/CMS/Work/Data2011/PromptReco/TTree_Onia2MuMu_V8_PromptReco_v4.root");

  TTree *tree = chain;
  
  TFile *fOut = fOut = new TFile(fNameOut, "RECREATE");

  PolData treeReco(tree);
  BookHistosReco();
  printf("after booking of histo\n");
  Int_t selDimuType = 4; //0...only GG, 1... only GT, 2... only TT, 3...GG+GT, 4...GG+GT+TT
  treeReco.Loop(selDimuType, rejectCowboys);
  printf("writing out the histograms\n");
  WriteHistosReco(fNameOut);

  fOut->Close();
}
//==========================================
void BookHistosReco(){

  //mass
  Int_t nBinsMass = 320;
  Double_t massMin = 8.4, massMax = 11.6;
  //pt
  Int_t nBinsPt = 1000;
  Double_t pTMin = 0., pTMaxOnia = 100.;
  //rap
  Int_t nBinsRap = 100;
  Double_t rapMin = -2.5, rapMax = 2.5;

  Char_t name[100], title[300];
  //statistics
  Reco_StatEv = new TH1F("Reco_StatEv", "", 12, 0., 12.);

  //reconstruction variables for the Onia
  for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){
      //Mass:
      sprintf(name, "Reco_Onia_mass_rap%d_pT%d", iRapBin, iPTBin);
      sprintf(title, ";M [GeV/c^{2}]");
      Reco_Onia_mass[iPTBin][iRapBin] = new TH1F(name, title, nBinsMass, massMin, massMax);
      Reco_Onia_mass[iPTBin][iRapBin]->Sumw2();
    }
  }

  sprintf(name, "Reco_Onia_rap_pt");
  sprintf(title, ";y(#mu#mu);p_{T}^{#mu#mu} [GeV/c]");
  Reco_Onia_rap_pT = new TH2F(name, title, nBinsRap,rapMin,rapMax, nBinsPt,pTMin,pTMaxOnia);
  Reco_Onia_rap_pT->Sumw2();

  for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
    //pT
    sprintf(name, "Reco_Onia_pt_rap%d", iRapBin);
    sprintf(title, ";p_{T}^{#mu#mu} [GeV/c]");
    Reco_Onia_pt[iRapBin]  = new TH1F(name, title, nBinsPt,pTMin,pTMaxOnia);
    Reco_Onia_pt[iRapBin]->Sumw2();
  }
  for(int iPTBin = 0; iPTBin < onia::kNbPTMaxBins+1; iPTBin++){
    //rap
    sprintf(name, "Reco_Onia_rap_pT%d", iPTBin);
    sprintf(title, ";y(#mu#mu)");
    Reco_Onia_rap[iPTBin]  = new TH1F(name, title, nBinsRap, rapMin,rapMax);
    Reco_Onia_rap[iPTBin]->Sumw2();
  }

  //prepare the branches for the output tree
  treeOut = new TTree ("selectedData", "selected events");
  lepP = new TLorentzVector();
  lepN = new TLorentzVector();
  treeOut->Branch("lepP", "TLorentzVector", &lepP);
  treeOut->Branch("lepN", "TLorentzVector", &lepN);

  // //polarization histos:
  // for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

  //   //pT and y double differential pol histos:
  //   for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
  //     for(int iPTBin = 0; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){
  // 	//2D histo:
  // 	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d", onia::frameLabel[iFrame], iPTBin, iRapBin);
  // 	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", onia::frameLabel[iFrame], onia::frameLabel[iFrame]);
  // 	Reco2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2F(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, 
  // 								   onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
  // 	Reco2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
  //     }
  //   }
  //   //pT and y double differential pol histos: FWD and BWD rapidities separately
  //   for(int iRapBin = 0; iRapBin < 2*onia::kNbRapBins; iRapBin++){
  //     Int_t matchRapBin = fabs(onia::kNbRapForPTBins - iRapBin);
  //     if(iRapBin >= onia::kNbRapForPTBins) matchRapBin += 1;
  //     for(int iPTBin = 0; iPTBin < onia::kNbPTBins[matchRapBin]+1; iPTBin++){
  // 	sprintf(name, "Reco_Onia_cosTh_%s_pT%d_rap%d_FWDBWD", onia::frameLabel[iFrame], iPTBin, iRapBin);
  // 	sprintf(title, ";cos#theta_{%s}", onia::frameLabel[iFrame]);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cosThPol] = new TH1F(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cosThPol]->Sumw2();
  // 	//
  // 	sprintf(name, "Reco_Onia_phi_%s_pT%d_rap%d_FWDBWD", onia::frameLabel[iFrame], iPTBin, iRapBin);
  // 	sprintf(title, ";#phi_{%s} [deg]", onia::frameLabel[iFrame]);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::phiPol] = new TH1F(name, title, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::phiPol]->Sumw2();
  // 	//
  // 	sprintf(name, "Reco_Onia_cos2Phi_%s_pT%d_rap%d_FWDBWD", onia::frameLabel[iFrame], iPTBin, iRapBin);
  // 	sprintf(title, ";cos(2#phi_{%s})", onia::frameLabel[iFrame]);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cos2PhiPol] = new TH1F(name, title, onia::kNbBinsCos2Phi, onia::cos2PhiMin, onia::cos2PhiMax);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cos2PhiPol]->Sumw2();
  // 	//2D histo:
  // 	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d_FWDBWD", onia::frameLabel[iFrame], iPTBin, iRapBin);
  // 	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", onia::frameLabel[iFrame], onia::frameLabel[iFrame]);
  // 	Reco2D_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin] = new TH2F(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, 
  // 								   onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
  // 	Reco2D_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin]->Sumw2();
  //     }
  //   }
  // }

  // for(int iRapBin = 1; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 1; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){
  //     //checking the rotation angle between HX and CS:
  //     sprintf(name, "Reco_hDelta_pT%d_rap%d", iPTBin, iRapBin);
  //     hDelta[iPTBin][iRapBin] = new TH1F(name, ";#delta(HX --> CS) [#circ]", 180, 0., 180.);
  //     hDelta[iPTBin][iRapBin]->Sumw2();
  //     sprintf(name, "Reco_hSin2Delta_pT%d_rap%d", iPTBin, iRapBin);
  //     hSin2Delta[iPTBin][iRapBin] = new TH1F(name, ";sin^{2}#delta(HX --> CS)", 100, 0., 1.);
  //     hSin2Delta[iPTBin][iRapBin]->Sumw2();
  //   }
  // }
}

//==========================================
void WriteHistosReco(Char_t *fNameOut){

  treeOut->Write();

  Reco_StatEv->Write();

  for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){  
      Reco_Onia_mass[iPTBin][iRapBin]->Write();
    }
  }

  Reco_Onia_rap_pT->Write();
  for(int iPTBin = 0; iPTBin < onia::kNbPTMaxBins+1; iPTBin++)
    Reco_Onia_rap[iPTBin]->Write();
  for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++)
    Reco_Onia_pt[iRapBin]->Write();


  // //polarization histos: Reco
  // for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
  //   for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
  //     for(int iPTBin = 0; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){
  // 	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][onia::cosThPol]->Write();
  // 	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][onia::phiPol]->Write();
  // 	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][onia::cos2PhiPol]->Write();
  // 	Reco2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
  //     }
  //   }
  // }
  // for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
  //   for(int iRapBin = 0; iRapBin < 2*onia::kNbRapBins; iRapBin++){
  //     Int_t matchRapBin = fabs(onia::kNbRapForPTBins - iRapBin);
  //     if(iRapBin >= onia::kNbRapForPTBins) matchRapBin += 1;
  //     for(int iPTBin = 0; iPTBin < onia::kNbPTBins[matchRapBin]+1; iPTBin++){
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cosThPol]->Write();
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::phiPol]->Write();
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cos2PhiPol]->Write();
  // 	Reco2D_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin]->Write();
  //     }
  //   }
  // }

  // for(int iRapBin = 1; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 1; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){
  //     hDelta[iPTBin][iRapBin]->Write();
  //     hSin2Delta[iPTBin][iRapBin]->Write();
  //   }
  // }

  // for(int iCut = 0; iCut < 4; iCut++){
  //   Reco_muHLT_pT_eta[iCut]->Write();
  //   Reco_muHLT_p_eta[iCut]->Write();
  //   Reco_muTM_pT_eta[iCut]->Write();
  //   Reco_muTM_p_eta[iCut]->Write();
  // }
}
