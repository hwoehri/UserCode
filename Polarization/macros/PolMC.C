#define PolMC_cxx
#include "PolMC.h"
#include "calcPol.C"

#include "TLorentzVector.h"
#include "TClonesArray.h"
#include <TH2.h>
#include <TCanvas.h>

//some statistics
TH1F *hGen_StatEv;

//histos at gen level for mu+, mu-, gamma, Onia
TH1F *hGen_mupl_pt[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH1F *hGen_mumi_pt[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH1F *hGen_mupl_eta[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH1F *hGen_mumi_eta[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH1F *hGen_mupl_phi[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH1F *hGen_mumi_phi[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH1F *hGen_Onia_mass[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH1F *hGen_Onia_pt[jpsi::kNbRapForPTBins+1];
TH1F *hGen_Onia_eta[jpsi::kNbPTBins+1];
TH1F *hGen_Onia_rap[jpsi::kNbPTBins+1];
TH1F *hGen_Onia_phi[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH2F *hGen_Onia_rap_pT;
//debugging histos:
TH2F *hPhiPos_PhiNeg[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH2F *hPtPos_PtNeg[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH2F *hEtaPos_EtaNeg[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH1F *hDeltaPhi[jpsi::kNbPTBins+1][2*jpsi::kNbRapBins+1];

//polarization histos:
TH1F *hGen_Onia_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbPolVar];
TH1F *hGen_Onia_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1][jpsi::kNbPolVar];
TH1F *hGen_Onia_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1][jpsi::kNbPolVar];
TH2F *hGen2D_Onia_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1];
TH2F *hGen2D_Onia_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1];
TH2F *hGen2D_Onia_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];

TH1F *hDelta[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH1F *hSin2Delta[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];

Double_t CalcPolWeight(Double_t pf_onia_P, Double_t thisCosTh);
//==============================================
void PolMC::Loop(Int_t selDimuType)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t countGenEvent = 0;
  Long64_t nb = 0;
  printf("number of entries = %d\n", (Int_t) nentries);

  //loop over the events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if(jentry % 100000 == 0) printf("event %d\n", (Int_t) jentry);

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    //protection against dummy events:
    if(muPosPx_Gen < -9900.)
      continue;

    hGen_StatEv->Fill(0.5);//count all events
    //do NOT select on the dimuon type (only a RECO variable)
    hGen_StatEv->Fill(1.5);//count all events

    Double_t enMuPos = sqrt(muPosPx_Gen*muPosPx_Gen + muPosPy_Gen*muPosPy_Gen + muPosPz_Gen*muPosPz_Gen + jpsi::muMass*jpsi::muMass);
    Double_t enMuNeg = sqrt(muNegPx_Gen*muNegPx_Gen + muNegPy_Gen*muNegPy_Gen + muNegPz_Gen*muNegPz_Gen + jpsi::muMass*jpsi::muMass);
    TLorentzVector *muPos = new TLorentzVector();
    TLorentzVector *muNeg = new TLorentzVector();
    muPos->SetPxPyPzE(muPosPx_Gen, muPosPy_Gen, muPosPz_Gen, enMuPos);
    muNeg->SetPxPyPzE(muNegPx_Gen, muNegPy_Gen, muNegPz_Gen, enMuNeg);
    Double_t etaMuPos = muPos->PseudoRapidity();
    Double_t etaMuNeg = muNeg->PseudoRapidity();
    Double_t pTMuPos = muPos->Pt();
    Double_t pTMuNeg = muNeg->Pt();

    //take muons only within a certain eta range
    if(TMath::Abs(etaMuPos) > jpsi::etaPS || TMath::Abs(etaMuNeg) > jpsi::etaPS){
      // printf("eta(pos. muon) = %f, eta(neg. muon) = %f\n", etaMuPos, etaMuNeg);
      continue;
    }
    hGen_StatEv->Fill(2.5);//count all events

    if(pTMuPos < jpsi::pTMuMin && pTMuNeg < jpsi::pTMuMin){
      // printf("pT(pos. muon) = %f, pT(neg. muon) = %f\n", pTMuPos, pTMuNeg);
      continue;
    }
    hGen_StatEv->Fill(3.5);

    //test according to Gavin's proposal:
    //if any of the two muons is within 1.4 < eta < 1.6 AND
    //the two muons are close in eta (deltaEta < 0.2)
    //reject the dimuon (no matter whether it is "Seagull" or 
    //"Cowboy"):
//     if(TMath::Abs(etaMuPos - etaMuNeg) < 0.2 &&
//        ((etaMuPos > 1.4 && etaMuPos < 1.6) || (etaMuNeg > 1.4 && etaMuNeg < 1.6))){
//       printf("rejecting the event!\n");
//       continue;
//     }

    //build the invariant mass, pt, ... of the two muons
    TLorentzVector *onia = new TLorentzVector();
    *onia = *(muPos) + *(muNeg);

    Double_t onia_mass = onia->M();
    Double_t onia_pt = onia->Pt();
    Double_t onia_P = onia->P();
    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();
    Double_t onia_phi = onia->Phi();
    Double_t onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);

    Int_t pTIndex = -1;
    for(int iPT = 0; iPT < jpsi::kNbPTBins; iPT++){
      if(onia_pt > jpsi::pTRange[iPT] && onia_pt < jpsi::pTRange[iPT+1]){
	pTIndex = iPT+1;
	break;
      }
    }
    Int_t rapIndex = -1;
    for(int iRap = 0; iRap < 2*jpsi::kNbRapBins; iRap++){
      if(onia_rap > jpsi::rapRange[iRap] && onia_rap < jpsi::rapRange[iRap+1]){
	rapIndex = iRap+1;
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

    if(pTIndex < 1){
      printf("pTIndex %d, pT(onia) = %f\n", pTIndex, onia_pt);
      continue;
    }
    if(rapIndex < 1){
      printf("rapIndex %d, rap(onia) = %f\n", rapIndex, onia_rap);
      continue;
    }
    if(rapForPTIndex < 1){
      printf("rapForPTIndex %d, rap(onia) = %f\n", rapForPTIndex, onia_rap);
      continue;
    }

    hGen_StatEv->Fill(4.5);

    if(TMath::Abs(onia_rap) > jpsi::rapYPS)
      continue;

    hGen_StatEv->Fill(7.5);

    //remaining of the events will be used for the analysis
    countGenEvent++;

    //fill mass, phi, pt, eta and rap distributions
    //a) all bins
    hGen_Onia_mass[0][0]->Fill(onia_mass);
    hGen_Onia_phi[0][0]->Fill(onia_phi);
    hGen_Onia_mass[pTIndex][0]->Fill(onia_mass);
    hGen_Onia_phi[pTIndex][0]->Fill(onia_phi);
    hGen_Onia_mass[0][rapForPTIndex]->Fill(onia_mass);
    hGen_Onia_phi[0][rapForPTIndex]->Fill(onia_phi);
      
    hGen_Onia_pt[0]->Fill(onia_pt);
    hGen_Onia_eta[0]->Fill(onia_eta);
    hGen_Onia_rap[0]->Fill(onia_rap);
    //b) individual pT and rap bins:
    hGen_Onia_mass[pTIndex][rapForPTIndex]->Fill(onia_mass);
    hGen_Onia_phi[pTIndex][rapForPTIndex]->Fill(onia_phi);
    hGen_Onia_pt[rapForPTIndex]->Fill(onia_pt);
    hGen_Onia_eta[pTIndex]->Fill(onia_eta);
    hGen_Onia_rap[pTIndex]->Fill(onia_rap);

    hGen_Onia_rap_pT->Fill(onia_rap, onia_pt);

    //=====================
    calcPol(*muPos, *muNeg);
    //=====================

    //test:
//     calcPol(*muNeg, *muPos);
//     //H: test:
//     if(jentry%2 == 0)
//       calcPol(*muPos, *muNeg);
//     else
//       calcPol(*muNeg, *muPos);

    //===================================================
    //calculate delta, the angle between the CS and HX frame
    //Formula from EPJC paper
    Double_t deltaHXToCS = TMath::ACos(onia_mass * onia->Pz() / (onia_mT * onia_P));
    //     Double_t deltaCSToHX = -deltaHXToCS;
    Double_t sin2Delta = pow((onia_pt * onia->Energy() / (onia_P * onia_mT)),2);
    //sin2Delta does not change sign when going from HX-->CS or vice versa
    hDelta[pTIndex][rapForPTIndex]->Fill(deltaHXToCS * 180./TMath::Pi());
    hSin2Delta[pTIndex][rapForPTIndex]->Fill(sin2Delta);
    //===================================================

    Double_t deltaPhi = muPos->Phi() - muNeg->Phi();
    if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
    else if(deltaPhi > TMath::Pi()) deltaPhi = 2.*TMath::Pi() - deltaPhi;

    //debugging histos
    hPhiPos_PhiNeg[pTIndex][rapForPTIndex]->Fill(180./TMath::Pi() * muNeg->Phi(), 180./TMath::Pi() * muPos->Phi());
    hPtPos_PtNeg[pTIndex][rapForPTIndex]->Fill(muNeg->Pt(), muPos->Pt());
    hEtaPos_EtaNeg[pTIndex][rapForPTIndex]->Fill(muNeg->PseudoRapidity(), muPos->PseudoRapidity());
    hDeltaPhi[pTIndex][rapIndex]->Fill(deltaPhi);

    hGen_mupl_pt[pTIndex][rapForPTIndex]->Fill(muPos->Pt());
    hGen_mupl_eta[pTIndex][rapForPTIndex]->Fill(muPos->PseudoRapidity());
    hGen_mupl_phi[pTIndex][rapForPTIndex]->Fill(muPos->Phi());
      
    hGen_mumi_pt[pTIndex][rapForPTIndex]->Fill(muNeg->Pt());
    hGen_mumi_eta[pTIndex][rapForPTIndex]->Fill(muNeg->PseudoRapidity());
    hGen_mumi_phi[pTIndex][rapForPTIndex]->Fill(muNeg->Phi());

    //fill the histos for all the different frames
    for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){

      thisCosPhi[iFrame] = TMath::Cos(2.*thisPhi_rad[iFrame]);

      Double_t weight = CalcPolWeight(onia_P, thisCosTh[iFrame]);

      //1a) polariztion histos - all pT
      hGen_Onia_pol_pT[iFrame][0][jpsi::cosThPol]->Fill(thisCosTh[iFrame], weight);
      hGen_Onia_pol_pT[iFrame][0][jpsi::phiPol]->Fill(thisPhi[iFrame], weight);
      hGen_Onia_pol_pT[iFrame][0][jpsi::cos2PhiPol]->Fill(thisCosPhi[iFrame], weight);
      hGen2D_Onia_pol_pT[iFrame][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);

      //1b) polariztion histos - pT Bin
      if(pTIndex > 0){
	hGen_Onia_pol_pT[iFrame][pTIndex][jpsi::cosThPol]->Fill(thisCosTh[iFrame], weight);
	hGen_Onia_pol_pT[iFrame][pTIndex][jpsi::phiPol]->Fill(thisPhi[iFrame], weight);
	hGen_Onia_pol_pT[iFrame][pTIndex][jpsi::cos2PhiPol]->Fill(thisCosPhi[iFrame], weight);
	hGen2D_Onia_pol_pT[iFrame][pTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      }

      //2a) polariztion histos - all Rap
      hGen_Onia_pol_rap[iFrame][0][jpsi::cosThPol]->Fill(thisCosTh[iFrame], weight);
      hGen_Onia_pol_rap[iFrame][0][jpsi::phiPol]->Fill(thisPhi[iFrame], weight);
      hGen_Onia_pol_rap[iFrame][0][jpsi::cos2PhiPol]->Fill(thisCosPhi[iFrame], weight);
      hGen2D_Onia_pol_rap[iFrame][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);

      //2b) polariztion histos - rap Bin
      if(rapIndex > 0){
	hGen_Onia_pol_rap[iFrame][rapIndex][jpsi::cosThPol]->Fill(thisCosTh[iFrame], weight);
	hGen_Onia_pol_rap[iFrame][rapIndex][jpsi::phiPol]->Fill(thisPhi[iFrame], weight);
	hGen_Onia_pol_rap[iFrame][rapIndex][jpsi::cos2PhiPol]->Fill(thisCosPhi[iFrame], weight);
	hGen2D_Onia_pol_rap[iFrame][rapIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      }

      //3) polariztion histos - pT and rap Bin
      //all pT and rapidities
      hGen2D_Onia_pol_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      if(pTIndex > 0)
	hGen2D_Onia_pol_pT_rap[iFrame][pTIndex][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      if(rapForPTIndex > 0)
	hGen2D_Onia_pol_pT_rap[iFrame][0][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      if(pTIndex > 0 && rapForPTIndex > 0){
	hGen_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex][jpsi::cosThPol]->Fill(thisCosTh[iFrame], weight);
	hGen_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex][jpsi::phiPol]->Fill(thisPhi[iFrame], weight);
	hGen_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex][jpsi::cos2PhiPol]->Fill(thisCosPhi[iFrame], weight);
	hGen2D_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      }
    }

    delete muPos;
    delete muNeg;
    delete onia;
  }//loop over entries

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countGenEvent, (Int_t) nentries);
}

//==========================================================
Double_t CalcPolWeight(Double_t onia_P, Double_t thisCosTh){

//   Double_t const p0 = 5.;
//   Double_t const kappa = 0.6;
//  //pol. acc. to PRL:
//   Double_t lambdaTh = 1. - TMath::Power(2., 1.-TMath::Power(onia_P/p0, kappa));


  Double_t lambdaTh = -1.;
//   Double_t weight = 1. + lambdaTh * TMath::Power(thisCosTh, 2.);
  Double_t weight = 1.;
  return weight;
}
