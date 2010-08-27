#ifndef __CINT__
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#endif

#include "../interface/rootIncludes.inc"
#include "PolMC.C"

void BookHistosMC(Char_t *oniaLabel);
void WriteHistosMC(Char_t *fNameOut);
//==========================================
//usage: root 'runMC.C+' or
//       root 'runMC.C+("pol_MC_HLT_Mu0Track0Jpsi.root", kTRUE)' (e.g.) 
//==========================================
void runMC(Char_t *fNameOut = "pol_MC_HLT_Mu0Track0Jpsi.root",
	   Char_t *fNameIn = "/home/hermine/CMS/Work/Polarization/Florian/25Aug2010/RooDataSet_Spring10_PromptJPsi_TEST_0_10.root",
	   Char_t *nameDataSet = "data", //"data" or "recoData"
	   Int_t selDimuType = 3, //0...only GG, 1... only GT, 2... only TT, 3...GG+GT, 4...GG+GT+TT
	   Char_t *oniaLabel = "J/#psi"){//"Ups(1S)"

  TFile *fIn = new TFile(fNameIn);
  RooDataSet* ds = (RooDataSet*)fIn->Get(nameDataSet);
  TTree *treeData = (TTree*)ds->tree();

  // TTree *treeData = (TTree*)fIn->Get("data");

  TFile *fOut = new TFile(fNameOut, "RECREATE");

  PolMC treeMC(treeData);
  BookHistosMC(oniaLabel);
  treeMC.Loop(selDimuType);
  WriteHistosMC(fNameOut);

  fOut->Close();
}

//==========================================
void BookHistosMC(Char_t *oniaLabel){

  //mass
  Int_t nBinsMass = 80;
  Double_t massMin = 8.0, massMax = 12.0;
  if(strncmp(oniaLabel, "J/#psi", 6) == 0){
    nBinsMass = 80;
    massMin = 2.7;
    massMax = 3.5;
  }
  //pt
  Int_t nBinsPt = 100, nBinsPtGamma = 30;
  Double_t pTMin = 0., pTMaxOnia = 30., pTMaxGamma = 3.0;
  //energy
  Int_t nBinsE = 50, nBinsEOnia = 200;
  Double_t enMin = 0., enMax = 5., enMaxOnia = 50.;
  //phi
  Int_t nBinsPhi = 157; //314
  Double_t phiMin = -3.14, phiMax = 3.14;
  //rap
  Int_t nBinsRap = 100;
  Double_t rapMin = -2.5, rapMax = 2.5;
  //pseudo-rap
  Int_t nBinsEtaGamma = 120;
  Double_t etaMinGamma = -3.0, etaMaxGamma = 3.0;

  Char_t name[100], title[300];
  //statistics
  hGen_StatEv = new TH1F("hGen_StatEv", "", 12, 0., 12.);

  //generated variables for the Onia
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
    for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++){
      //Mass:
      sprintf(name, "hGen_Onia_mass_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";M [GeV/c^{2}]");
      hGen_Onia_mass[iPTBin][iRapBin] = new TH1F(name, title, nBinsMass,massMin,massMax);
      hGen_Onia_mass[iPTBin][iRapBin]->Sumw2();
      //phi-lab:
      sprintf(name, "hGen_Onia_phi_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";%s #phi (rad)", oniaLabel);
      hGen_Onia_phi[iPTBin][iRapBin] = new TH1F(name, title, nBinsPhi,phiMin,phiMax);
      hGen_Onia_phi[iPTBin][iRapBin]->Sumw2();
    }
    //eta:
    sprintf(name, "hGen_Onia_eta_pT%d", iPTBin);
    sprintf(title, ";%s #eta", oniaLabel);
    hGen_Onia_eta[iPTBin] = new TH1F(name, title,nBinsEtaGamma,etaMinGamma,etaMaxGamma);
    hGen_Onia_eta[iPTBin]->Sumw2();
    //rap
    sprintf(name, "hGen_Onia_rap_pT%d", iPTBin);
    sprintf(title, ";y(%s)", oniaLabel);
    hGen_Onia_rap[iPTBin] = new TH1F(name, title, nBinsRap,rapMin,rapMax);
    hGen_Onia_rap[iPTBin]->Sumw2();
  }
  for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++){
    //pT
    sprintf(name, "hGen_Onia_pt_rap%d", iRapBin);
    sprintf(title, ";%s p_{T} [GeV/c]", oniaLabel);
    hGen_Onia_pt[iRapBin]  = new TH1F(name, title, nBinsPt,pTMin,pTMaxOnia);
    hGen_Onia_pt[iRapBin]->Sumw2();
  }

  //2D Onia histos:
  sprintf(name, "hGen_Onia_rap_pt");
  sprintf(title, ";%s y;p_{T} [GeV/c]", oniaLabel);
  hGen_Onia_rap_pT = new TH2F(name, title, nBinsRap,rapMin,rapMax, nBinsPt,pTMin,pTMaxOnia);
  hGen_Onia_rap_pT->Sumw2();


  //debugging histos (single Muons):
  for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
    for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){

      sprintf(name, "hGen_mupl_pt_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;p_{T}(#mu^{+})[GeV/c]",
          rapForPTRange[iRapBin-1], oniaLabel, rapForPTRange[iRapBin],
          pTRange[iPTBin-1], oniaLabel, pTRange[iPTBin]);
      hGen_mupl_pt[iPTBin][iRapBin]  = new TH1F(name, title,nBinsPt,pTMin,pTMaxOnia);
      sprintf(name,"hGen_mupl_eta_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;#eta(#mu^{+})",
          rapForPTRange[iRapBin-1], oniaLabel, rapForPTRange[iRapBin],
          pTRange[iPTBin-1], oniaLabel, pTRange[iPTBin]);
      hGen_mupl_eta[iPTBin][iRapBin] = new TH1F(name,title,nBinsRap,rapMin,rapMax);
      sprintf(name,"hGen_mupl_phi_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;#phi(#mu^{+})",
          rapForPTRange[iRapBin-1], oniaLabel, rapForPTRange[iRapBin],
          pTRange[iPTBin-1], oniaLabel, pTRange[iPTBin]);
      hGen_mupl_phi[iPTBin][iRapBin] = new TH1F(name,title, nBinsPhi,phiMin,phiMax);

      sprintf(name,"hGen_mumi_pt_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;p_{T}(#mu^{-}) [GeV/c]",
          rapForPTRange[iRapBin-1], oniaLabel, rapForPTRange[iRapBin],
          pTRange[iPTBin-1], oniaLabel, pTRange[iPTBin]);
      hGen_mumi_pt[iPTBin][iRapBin]  = new TH1F(name, title,nBinsPt,pTMin,pTMaxOnia);
      sprintf(name,"hGen_mumi_eta_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;#eta(#mu^{-})",
          rapForPTRange[iRapBin-1], oniaLabel, rapForPTRange[iRapBin],
          pTRange[iPTBin-1], oniaLabel, pTRange[iPTBin]);
      hGen_mumi_eta[iPTBin][iRapBin] = new TH1F(name,title,nBinsRap,rapMin,rapMax);
      sprintf(name,"hGen_mumi_phi_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;#phi(#mu^{-})",
          rapForPTRange[iRapBin-1], oniaLabel, rapForPTRange[iRapBin],
          pTRange[iPTBin-1], oniaLabel, pTRange[iPTBin]);
      hGen_mumi_phi[iPTBin][iRapBin] = new TH1F(name,title, nBinsPhi,phiMin,phiMax);

      hGen_mupl_pt[iPTBin][iRapBin] ->Sumw2();
      hGen_mupl_eta[iPTBin][iRapBin]->Sumw2();
      hGen_mupl_phi[iPTBin][iRapBin]->Sumw2();

      hGen_mumi_pt[iPTBin][iRapBin] ->Sumw2();
      hGen_mumi_eta[iPTBin][iRapBin]->Sumw2();
      hGen_mumi_phi[iPTBin][iRapBin]->Sumw2();

      sprintf(name, "hGen_hPhiPos_PhiNeg_pT%d_rap%d" , iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < |p_{T}(%s)| < %1.1f GeV/c;#phi(#mu^{-});#phi(#mu^{+})",
          rapForPTRange[iRapBin-1], oniaLabel, rapForPTRange[iRapBin],
          pTRange[iPTBin-1], oniaLabel, pTRange[iPTBin]);
      hPhiPos_PhiNeg[iPTBin][iRapBin] = new TH2F(name, title, 60,-180.,180., 60,-180.,180.);
      sprintf(name, "hGen_hPtPos_PtNeg_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < |p_{T}(%s)| < %1.1f GeV/c;p_{T}(#mu^{-});p_{T}(#mu^{+})",
          rapForPTRange[iRapBin-1], oniaLabel, rapForPTRange[iRapBin],
          pTRange[iPTBin-1], oniaLabel, pTRange[iPTBin]);
      hPtPos_PtNeg[iPTBin][iRapBin] = new TH2F(name, title, 20, 0., 10., 20, 0., 10.);
      sprintf(name, "hGen_hEtaPos_EtaNeg_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < |p_{T}(%s)| < %1.1f GeV/c;#eta(#mu^{-});#eta(#mu^{+})",
          rapForPTRange[iRapBin-1], oniaLabel, rapForPTRange[iRapBin],
          pTRange[iPTBin-1], oniaLabel, pTRange[iPTBin]);
      hEtaPos_EtaNeg[iPTBin][iRapBin] = new TH2F(name, title, 24, -2.4, 2.4, 24, -2.4, 2.4);
    }
    for(int iRapBin = 1; iRapBin < 2*kNbRapBins+1; iRapBin++){
      sprintf(name, "hGen_hDeltaPhi_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < y(%s) < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;#phi(#mu^{+}) - #phi(#mu^{-})",
          rapRange[iRapBin-1], oniaLabel, rapRange[iRapBin],
          pTRange[iPTBin-1], oniaLabel, pTRange[iPTBin]);
      hDeltaPhi[iPTBin][iRapBin] = new TH1F(name, title, 96, -1.6, 1.6);
      hDeltaPhi[iPTBin][iRapBin]->Sumw2();
    }
  }

  //polarization histos:
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
      sprintf(name, "hGen_Onia_cosTh_%s_pT%d", frameLabel[iFrame], iPTBin);
      sprintf(title, ";cos#theta_{%s}", frameLabel[iFrame]);
      hGen_Onia_pol_pT[iFrame][iPTBin][cosThPol] = new TH1F(name, title, nBinsCosT, cosTMin, cosTMax);
      hGen_Onia_pol_pT[iFrame][iPTBin][cosThPol]->Sumw2();
      //
      sprintf(name, "hGen_Onia_phi_%s_pT%d", frameLabel[iFrame], iPTBin);
      sprintf(title, ";#phi_{%s} [deg]", frameLabel[iFrame]);
      hGen_Onia_pol_pT[iFrame][iPTBin][phiPol] = new TH1F(name, title, nBinsPhiPol, phiPolMin, phiPolMax);
      hGen_Onia_pol_pT[iFrame][iPTBin][phiPol]->Sumw2();
      //
      sprintf(name, "hGen_Onia_cos2Phi_%s_pT%d", frameLabel[iFrame], iPTBin);
      sprintf(title, ";cos(2#phi_{%s})", frameLabel[iFrame]);
      hGen_Onia_pol_pT[iFrame][iPTBin][cos2PhiPol] = new TH1F(name, title, nBinsCos2Phi, cos2PhiMin, cos2PhiMax);
      hGen_Onia_pol_pT[iFrame][iPTBin][cos2PhiPol]->Sumw2();
      //2D histo:
      sprintf(name, "hGen2D_Onia_%s_pT%d", frameLabel[iFrame], iPTBin);
      sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", frameLabel[iFrame], frameLabel[iFrame]);
      hGen2D_Onia_pol_pT[iFrame][iPTBin] = new TH2F(name, title, nBinsCosT, cosTMin, cosTMax,
                            nBinsPhiPol, phiPolMin, phiPolMax);
      hGen2D_Onia_pol_pT[iFrame][iPTBin]->Sumw2();
    }
    for(int iRapBin = 0; iRapBin < 2*kNbRapBins+1; iRapBin++){
      sprintf(name, "hGen_Onia_cosTh_%s_rap%d", frameLabel[iFrame], iRapBin);
      sprintf(title, ";cos#theta_{%s}", frameLabel[iFrame]);
      hGen_Onia_pol_rap[iFrame][iRapBin][cosThPol] = new TH1F(name, title, nBinsCosT, cosTMin, cosTMax);
      hGen_Onia_pol_rap[iFrame][iRapBin][cosThPol]->Sumw2();
      //
      sprintf(name, "hGen_Onia_phi_%s_rap%d", frameLabel[iFrame], iRapBin);
      sprintf(title, ";#phi_{%s} [deg]", frameLabel[iFrame]);
      hGen_Onia_pol_rap[iFrame][iRapBin][phiPol] = new TH1F(name, title, nBinsPhiPol, phiPolMin, phiPolMax);
      hGen_Onia_pol_rap[iFrame][iRapBin][phiPol]->Sumw2();
      //
      sprintf(name, "hGen_Onia_cos2Phi_%s_rap%d", frameLabel[iFrame], iRapBin);
      sprintf(title, ";cos(2#phi_{%s})", frameLabel[iFrame]);
      hGen_Onia_pol_rap[iFrame][iRapBin][cos2PhiPol] = new TH1F(name, title, nBinsCos2Phi, cos2PhiMin, cos2PhiMax);
      hGen_Onia_pol_rap[iFrame][iRapBin][cos2PhiPol]->Sumw2();
      //2D histo:
      sprintf(name, "hGen2D_Onia_%s_rap%d", frameLabel[iFrame], iRapBin);
      sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", frameLabel[iFrame], frameLabel[iFrame]);
      hGen2D_Onia_pol_rap[iFrame][iRapBin] = new TH2F(name, title, nBinsCosT, cosTMin, cosTMax,
                             nBinsPhiPol, phiPolMin, phiPolMax);
      hGen2D_Onia_pol_rap[iFrame][iRapBin]->Sumw2();
    }
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
    sprintf(name, "hGen_Onia_cosTh_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
    sprintf(title, ";cos#theta_{%s}", frameLabel[iFrame]);
    hGen_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol] = new TH1F(name, title, nBinsCosT, cosTMin, cosTMax);
    hGen_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->Sumw2();
    //
    sprintf(name, "hGen_Onia_phi_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
    sprintf(title, ";#phi_{%s} [deg]", frameLabel[iFrame]);
    hGen_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol] = new TH1F(name, title, nBinsPhiPol, phiPolMin, phiPolMax);
    hGen_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->Sumw2();
    //
    sprintf(name, "hGen_Onia_cos2Phi_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
    sprintf(title, ";cos(2#phi_{%s})", frameLabel[iFrame]);
    hGen_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][cos2PhiPol] = new TH1F(name, title, nBinsCos2Phi, cos2PhiMin, cos2PhiMax);
    hGen_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][cos2PhiPol]->Sumw2();
    //2D histo:
    sprintf(name, "hGen2D_Onia_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
    sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", frameLabel[iFrame], frameLabel[iFrame]);
    hGen2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2F(name, title, nBinsCosT, cosTMin, cosTMax,
                                   nBinsPhiPol, phiPolMin, phiPolMax);
    hGen2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();

      }
    }
  }

  for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
    for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
      //checking the rotation angle between HX and CS:
      sprintf(name, "hGen_hDelta_pT%d_rap%d", iPTBin, iRapBin);
//       hDelta[iPTBin][iRapBin] = new TH1F(name, ";#delta(HX --> CS) [rad]", 64, 0., 3.2);
      hDelta[iPTBin][iRapBin] = new TH1F(name, ";#delta(HX --> CS) [#circ]", 180, 0., 180.);
      hDelta[iPTBin][iRapBin]->Sumw2();
      sprintf(name, "hGen_hSin2Delta_pT%d_rap%d", iPTBin, iRapBin);
      hSin2Delta[iPTBin][iRapBin] = new TH1F(name, ";sin^{2}#delta(HX --> CS)", 100, 0., 1.);
      hSin2Delta[iPTBin][iRapBin]->Sumw2();
    }
  }
}

//==========================================
void WriteHistosMC(Char_t *fNameOut){

  hGen_StatEv->Write();

  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++)
    for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++)
      hGen_Onia_mass[iPTBin][iRapBin]->Write();
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++)
    for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++)
      hGen_Onia_phi[iPTBin][iRapBin]->Write();

  for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++)
    hGen_Onia_pt[iRapBin]->Write();
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++)
    hGen_Onia_eta[iPTBin] ->Write();
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++)
    hGen_Onia_rap[iPTBin] ->Write();

  hGen_Onia_rap_pT->Write();

  for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
     for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
       //debugging histos (single muons):
       hGen_mupl_pt[iPTBin][iRapBin] ->Write();
       hGen_mupl_eta[iPTBin][iRapBin]->Write();
       hGen_mupl_phi[iPTBin][iRapBin]->Write();

       hGen_mumi_pt[iPTBin][iRapBin] ->Write();
       hGen_mumi_eta[iPTBin][iRapBin]->Write();
       hGen_mumi_phi[iPTBin][iRapBin]->Write();

       hPhiPos_PhiNeg[iPTBin][iRapBin]->Write();
       hPtPos_PtNeg[iPTBin][iRapBin]->Write();
       hEtaPos_EtaNeg[iPTBin][iRapBin]->Write();
     }
     for(int iRapBin = 1; iRapBin < 2*kNbRapBins+1; iRapBin++)
       hDeltaPhi[iPTBin][iRapBin]->Write();
  }

  //polarization histos: hGen
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
      hGen_Onia_pol_pT[iFrame][iPTBin][cosThPol]->Write();
      hGen_Onia_pol_pT[iFrame][iPTBin][phiPol]->Write();
      hGen_Onia_pol_pT[iFrame][iPTBin][cos2PhiPol]->Write();
      hGen2D_Onia_pol_pT[iFrame][iPTBin]->Write();
    }
    for(int iRapBin = 0; iRapBin < 2*kNbRapBins+1; iRapBin++){
      hGen_Onia_pol_rap[iFrame][iRapBin][cosThPol]->Write();
      hGen_Onia_pol_rap[iFrame][iRapBin][phiPol]->Write();
      hGen_Onia_pol_rap[iFrame][iRapBin][cos2PhiPol]->Write();
      hGen2D_Onia_pol_rap[iFrame][iRapBin]->Write();
    }
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
    hGen_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->Write();
    hGen_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->Write();
    hGen_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][cos2PhiPol]->Write();
    hGen2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
    for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
      hDelta[iPTBin][iRapBin]->Write();
      hSin2Delta[iPTBin][iRapBin]->Write();
    }
  }
}



