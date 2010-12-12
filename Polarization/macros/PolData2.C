#define PolData2_cxx
#include "PolData2.h"
#include "calcPol.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>

// Int_t eventNb, runNb, lumiBlock;
// Int_t nPriVtx;
// Int_t JpsiType, JpsiCharge, Jpsict, JpsictErr, JpsiVprob;
// Int_t HLT_Mu0_Track0_Jpsi, HLT_Mu3_Track0_Jpsi, HLT_Mu5_Track0_Jpsi;
// Int_t HLT_Mu0_TkMu0_Jpsi, HLT_Mu3_TkMu0_Jpsi, HLT_Mu5_TkMu0_Jpsi;
// Int_t HLT_Mu0_TkMu0_OST_Jpsi, HLT_Mu3_TkMu0_OST_Jpsi
//some statistics
TH1F *Reco_StatEv;

//histos at reco level for mu+, mu-, gamma, Onia
TH1F *Reco_mupl_pt[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_mumi_pt[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_mupl_eta[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_mumi_eta[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_mupl_phi[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_mumi_phi[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_Onia_mass[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_Onia_pt[jpsi::kNbRapForPTBins+1];
TH1F *Reco_Onia_phi[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2F *Reco_Onia_rap_pT;
//debugging histos:
TH2F *hPhiPos_PhiNeg[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2F *hPtPos_PtNeg[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2F *hEtaPos_EtaNeg[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *hDeltaPhi[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];

//polarization histos:
TH1F *Reco_Onia_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1][jpsi::kNbPolVar];
TH2F *Reco2D_Onia_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
//same pol histos, but split into FWD and BWD rapidity:
TH1F *Reco_Onia_pol_pT_rap_FWDBWD[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][2*jpsi::kNbRapBins][jpsi::kNbPolVar];
TH2F *Reco2D_Onia_pol_pT_rap_FWDBWD[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][2*jpsi::kNbRapBins];

//
TH1F *hDelta[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *hSin2Delta[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];

Double_t CalcPolWeight(Double_t pf_onia_P, Double_t thisCosTh);
//==============================================
void PolData2::Loop(Int_t selDimuType, Bool_t writeOutEvents)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t cutAtRecEvent = nentries;
  Long64_t countRecEvent = 0;
  Long64_t nb = 0;
  printf("number of entries = %d\n", (Int_t) nentries);
  FILE *fOutputTextFile;
  if(writeOutEvents)
    fOutputTextFile = fopen("jPsiCandidates.txt", "write");

  //loop over the events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if(jentry % 100000 == 0) printf("event %d\n", (Int_t) jentry);

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    Reco_StatEv->Fill(0.5);//count all events

    //check the trigger flag: 0... no trigger, 1 ... triggered+matched, 2 ... triggered
    //for a full list of accessible triggers, check out "PolData2.h"
    //and https://espace.cern.ch/cms-quarkonia/onia-polarization/L1%20%20HLT/unprescaledTriggersVsRun.aspx
    //for the run ranges per HLT path during unprescaled running periods
    if(runNb >= 140116 && runNb <= 144114){
      if(HLT_Mu0_TkMu0_Jpsi != 1) continue;
    }
    else if(runNb >= 146428 && runNb <= 148058){
      if(HLT_Mu0_TkMu0_OST_Jpsi != 1) continue;
    }
    else if(runNb >= 148819 && runNb <= 149182){ 
      if(HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2 != 1) continue;    
    }
    else if(runNb >= 149291 && runNb <= 149442){
      if(HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3 != 1) continue;    
    }
    else{
      printf("rejecting events in run %d\n", runNb);
      continue;
    }
    //for usage of the HLT_DoubleMu0* paths:
//     if(runNb >= 133446 && runNb <= 147116){
//       if(HLT_DoubleMu0 != 1)
// 	continue;
//     }
//     else if(runNb >= 147196 && runNb <= 149442){
//       if(HLT_DoubleMu0_Quarkonium_v1 != 1)
// 	continue;
//     }
//     else{
//       printf("rejecting events in run %d\n", runNb);
//     continue;
//   }

    Reco_StatEv->Fill(1.5);

    //reject processing of events where the dimuon type (GG, GT or TT)
    //does not correspond to the chosen one
    if(selDimuType < 3 && JpsiType != selDimuType)
      continue;
    else if(selDimuType == 3 && JpsiType > 1) //only GG or GT
      continue;

    Reco_StatEv->Fill(2.5);//count all events

    Double_t etaMuPos = muPos->PseudoRapidity();
    Double_t etaMuNeg = muNeg->PseudoRapidity();
    Double_t pTMuPos = muPos->Pt();
    Double_t pTMuNeg = muNeg->Pt();
    Double_t pMuPos = muPos->P();
    Double_t pMuNeg = muNeg->P();

    //take muons only within a certain eta range
    if((fabs(etaMuPos) < jpsi::etaPS[0] && pTMuPos < jpsi::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuPos) > jpsi::etaPS[0] && fabs(etaMuPos) < jpsi::etaPS[1] && pMuPos < jpsi::pMuMin[1]) ||
       (fabs(etaMuPos) > jpsi::etaPS[1] && fabs(etaMuPos) < jpsi::etaPS[2] && pTMuPos < jpsi::pTMuMin[1]))
      continue;
    //(b) on the negative muon
    if((fabs(etaMuNeg) < jpsi::etaPS[0] && pTMuNeg < jpsi::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuNeg) > jpsi::etaPS[0] && fabs(etaMuNeg) < jpsi::etaPS[1] && pMuNeg < jpsi::pMuMin[1]) ||
       (fabs(etaMuNeg) > jpsi::etaPS[1] && fabs(etaMuNeg) < jpsi::etaPS[2] && pTMuNeg < jpsi::pTMuMin[1]))
      continue;

    Reco_StatEv->Fill(4.5);

    Double_t onia_mass = onia->M();
    Double_t onia_pt = onia->Pt();
    Double_t onia_P = onia->P();
    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();
    Double_t onia_phi = onia->Phi();
    Double_t onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);

    Int_t rapIndex = -1;
    for(int iRap = 0; iRap < 2*jpsi::kNbRapBins; iRap++){
      if(onia_rap > jpsi::rapRange[iRap] && onia_rap < jpsi::rapRange[iRap+1]){
    	rapIndex = iRap;
    	break;
      }
    }
    Int_t rapForPTIndex = -1;
    for(int iRap = 0; iRap < jpsi::kNbRapForPTBins; iRap++){
      if(TMath::Abs(onia_rap) > jpsi::rapForPTRange[iRap] && 
	 TMath::Abs(onia_rap) < jpsi::rapForPTRange[iRap+1]){
	 rapForPTIndex = iRap+1;
	break;
      }
    }
    Int_t pTIndex = -1;
    for(int iPT = 0; iPT < jpsi::kNbPTBins[rapForPTIndex]; iPT++){
      if(onia_pt > jpsi::pTRange[rapForPTIndex][iPT] && onia_pt < jpsi::pTRange[rapForPTIndex][iPT+1]){
	pTIndex = iPT+1;
	break;
      }
    }
    Int_t rapIntegratedPTIndex = -1;
    for(int iPT = 0; iPT < jpsi::kNbPTBins[0]; iPT++){
      if(onia_pt > jpsi::pTRange[0][iPT] && onia_pt < jpsi::pTRange[0][iPT+1]){
	rapIntegratedPTIndex = iPT+1;
	break;
      }
    }
    if(rapIndex < 0){
      // printf("rapIndex %d, rap(onia) = %f\n", rapIndex, onia_rap);
      continue;
    }
    if(rapForPTIndex < 1){
      // printf("rapForPTIndex %d, rap(onia) = %f\n", rapForPTIndex, onia_rap);
      continue;
    }
    if(pTIndex < 1){
      // printf("pTIndex %d, pT(onia) = %f\n", pTIndex, onia_pt);
      continue;
    }

    Reco_StatEv->Fill(5.5);

    //select events with a cut on the lifetime to reject NP J/psis:
    if(Jpsict > jpsi::JpsiCtauMax)
      continue;

    Reco_StatEv->Fill(6.5);

    //select events within a narrow mass window around the J/psi
    //(rapidity dependence of the resolution --> different mass windows)
    //values obtained from Gaussian fits in "plotMass.C"
    Double_t jPsiMassMin = jpsi::polMassJpsi[rapForPTIndex] - jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
    Double_t jPsiMassMax = jpsi::polMassJpsi[rapForPTIndex] + jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
    if(onia_mass < jPsiMassMin || onia_mass > jPsiMassMax)
      continue;

    Reco_StatEv->Fill(7.5);

    if(TMath::Abs(onia_rap) > jpsi::rapYPS)
      continue;

    Reco_StatEv->Fill(8.5);

    // //test: reject all Cowboy dimuons
    // if((muPos->Phi() - muNeg->Phi()) > 0)
    //   continue;

    // Reco_StatEv->Fill(9.5);

    //remaining of the events will be used for the analysis
    countRecEvent++;
    if(writeOutEvents)
      fprintf(fOutputTextFile, "%d\t%d\t%d\n", (Int_t) runNb, (Int_t) lumiBlock, (Int_t)  eventNb);

    //fill mass, phi, pt, eta and rap distributions
    //a) all bins
    Reco_Onia_mass[0][0]->Fill(onia_mass);
    Reco_Onia_mass[rapIntegratedPTIndex][0]->Fill(onia_mass);
    Reco_Onia_mass[0][rapForPTIndex]->Fill(onia_mass);
    //
    if(onia_mass > (jpsi::polMassJpsi[rapForPTIndex] - 2.*jpsi::sigmaMassJpsi[rapForPTIndex]) &&
       onia_mass < (jpsi::polMassJpsi[rapForPTIndex] + 2.*jpsi::sigmaMassJpsi[rapForPTIndex])){
	 Reco_Onia_phi[0][0]->Fill(onia_phi);
	 Reco_Onia_phi[rapIntegratedPTIndex][0]->Fill(onia_phi);
	 Reco_Onia_phi[0][rapForPTIndex]->Fill(onia_phi);
	 //      
	 Reco_Onia_pt[0]->Fill(onia_pt);
    }
    //b) individual pT and rap bins:
    Reco_Onia_mass[pTIndex][rapForPTIndex]->Fill(onia_mass);
    if(onia_mass > (jpsi::polMassJpsi[rapForPTIndex] - 2.*jpsi::sigmaMassJpsi[rapForPTIndex]) &&
       onia_mass < (jpsi::polMassJpsi[rapForPTIndex] + 2.*jpsi::sigmaMassJpsi[rapForPTIndex])){

      Reco_Onia_phi[pTIndex][rapForPTIndex]->Fill(onia_phi);
      Reco_Onia_pt[rapForPTIndex]->Fill(onia_pt);

      Reco_Onia_rap_pT->Fill(onia_rap, onia_pt);
    }

    //=====================
    calcPol(*muPos, *muNeg);
    //=====================

    //===================================================
    //calculate delta, the angle between the CS and HX frame
    //Formula from: EJP C69 (2010) 657
    Double_t deltaHXToCS = TMath::ACos(onia_mass * onia->Pz() / (onia_mT * onia_P));
    //     Double_t deltaCSToHX = -deltaHXToCS;
    Double_t sin2Delta = pow((onia_pt * onia->Energy() / (onia_P * onia_mT)),2);
    //sin2Delta does not change sign when going from HX-->CS or vice versa
    hDelta[pTIndex][rapForPTIndex]->Fill(deltaHXToCS * 180./TMath::Pi());
    hSin2Delta[pTIndex][rapForPTIndex]->Fill(sin2Delta);
    //===================================================

    Double_t deltaPhi = muPos->Phi() - muNeg->Phi();
    if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
    else if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();

    //debugging histos
    hPhiPos_PhiNeg[pTIndex][rapForPTIndex]->Fill(180./TMath::Pi() * muNeg->Phi(), 180./TMath::Pi() * muPos->Phi());
    hPtPos_PtNeg[pTIndex][rapForPTIndex]->Fill(muNeg->Pt(), muPos->Pt());
    hEtaPos_EtaNeg[pTIndex][rapForPTIndex]->Fill(muNeg->PseudoRapidity(), muPos->PseudoRapidity());
    hDeltaPhi[pTIndex][rapForPTIndex]->Fill(deltaPhi);

    Reco_mupl_pt[pTIndex][rapForPTIndex]->Fill(muPos->Pt());
    Reco_mupl_eta[pTIndex][rapForPTIndex]->Fill(muPos->PseudoRapidity());
    Reco_mupl_phi[pTIndex][rapForPTIndex]->Fill(muPos->Phi());
      
    Reco_mumi_pt[pTIndex][rapForPTIndex]->Fill(muNeg->Pt());
    Reco_mumi_eta[pTIndex][rapForPTIndex]->Fill(muNeg->PseudoRapidity());
    Reco_mumi_phi[pTIndex][rapForPTIndex]->Fill(muNeg->Phi());

    //fill the histos for all the different frames
    for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){

      thisCosPhi[iFrame] = TMath::Cos(2.*thisPhi_rad[iFrame]);

      Double_t weight = CalcPolWeight(onia_P, thisCosTh[iFrame]);

      //3) polariztion histos - pT and rap Bin
      //all pT and rapidities
      Reco2D_Onia_pol_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      if(rapIntegratedPTIndex > 0)
	Reco2D_Onia_pol_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      if(rapForPTIndex > 0)
	Reco2D_Onia_pol_pT_rap[iFrame][0][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      if(pTIndex > 0 && rapForPTIndex > 0){
	Reco_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex][jpsi::cosThPol]->Fill(thisCosTh[iFrame], weight);
	Reco_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex][jpsi::phiPol]->Fill(thisPhi[iFrame], weight);
	Reco_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex][jpsi::cos2PhiPol]->Fill(thisCosPhi[iFrame], weight);
	Reco2D_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      }
      //4) polarization histos - pT and rap Bin for FWD and BWD rapidities, separately 
      if(pTIndex > 0 && rapIndex >= 0){
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][pTIndex][rapIndex][jpsi::cosThPol]->Fill(thisCosTh[iFrame], weight);
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][pTIndex][rapIndex][jpsi::phiPol]->Fill(thisPhi[iFrame], weight);
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][pTIndex][rapIndex][jpsi::cos2PhiPol]->Fill(thisCosPhi[iFrame], weight);
	Reco2D_Onia_pol_pT_rap_FWDBWD[iFrame][pTIndex][rapIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      }
    }
  }//loop over entries

  if(writeOutEvents)
    fclose(fOutputTextFile);

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}

//==========================================================
Double_t CalcPolWeight(Double_t onia_P, Double_t thisCosTh){

  Double_t lambdaTh = 0.;
  Double_t weight = 1. + lambdaTh * TMath::Power(thisCosTh, 2.);
  return weight;
}
