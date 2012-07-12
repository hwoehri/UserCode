#define MCTruthEff_cxx
#include "MCTruthEff.h"
#include "calcPol.C"
#include "isMuonInAcceptance.C"
#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>
#include "TEfficiency.h"

//=================================================
//1D histos versus J/psi pT, phi and y; cosTheta and phi [kNbFrames]
//=================================================
//histos for neg. and pos. rapidity separately:
TEfficiency *recoEff_pT, *recoEff_y, *recoEff_phi, *recoEff_cosTheta[eff::kNbFrames], *recoEff_phiPol[eff::kNbFrames];
TEfficiency *trigEff_pT, *trigEff_y, *trigEff_phi, *trigEff_cosTheta[eff::kNbFrames], *trigEff_phiPol[eff::kNbFrames];
TEfficiency *totEff_pT, *totEff_y, *totEff_phi, *totEff_cosTheta[eff::kNbFrames], *totEff_phiPol[eff::kNbFrames];
//=================================================
//2D histos versus J/psi pT and y
//=================================================
//histos for neg. and pos. rapidity separately:
TEfficiency *recoEff2D_pT_rapNP, *recoEff2D_pT_rap, *recoEff2D_cosTheta_phiPol[eff::kNbFrames];
TEfficiency *trigEff2D_pT_rapNP, *trigEff2D_pT_rap, *trigEff2D_cosTheta_phiPol[eff::kNbFrames];
TEfficiency *totEff2D_pT_rapNP, *totEff2D_pT_rap, *totEff2D_cosTheta_phiPol[eff::kNbFrames];
//=================================================
//2D polarization histos, for various J/psi pT and y
//=================================================
// //histos for neg. and pos. rapidity separately:
TEfficiency *recoEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
//histos taking together +y and -y:
TEfficiency *recoEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
//histos taking together +y and -y and 4-folding in phi
TEfficiency *recoEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TEfficiency *recoEff_phiPol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff_phiPol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff_phiPol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *recoEff_cosTheta_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff_cosTheta_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff_cosTheta_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

//=================================================
//deltaR, deltaPhi vs deltaEta for various J/psi pT and y
//=================================================
// TEfficiency *recoEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *trigEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *totEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

// TEfficiency *recoEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *trigEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *totEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

// TEfficiency *recoEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *trigEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *totEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

// TEfficiency *recoEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *trigEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *totEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

// TEfficiency *recoEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *trigEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *totEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

// TEfficiency *recoEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *trigEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TEfficiency *totEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

enum {LOOSE,TIGHT};//set of muon fiducial cuts
enum {JPSI, PSIP, UPS1S, UPS2S, UPS3S};
Double_t pTMin = 10.; //5 GeV/c for Upsilon; 10 GeV/c for J/psi

//==============================================
void MCTruthEff::Loop(Int_t resonance, Bool_t rejectCowboys, Int_t useSoftMuons)
{
  if (fChain == 0) return;

  if(resonance == JPSI)
    pTMin = 10.;
  else if(resonance == UPS1S || resonance == UPS2S || resonance == UPS3S)
    pTMin = 5.;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nb = 0;
  Long64_t countRecEvent = 0;
  printf("number of entries = %d\n", (Int_t) nentries);

  //loop over the events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  //for (Long64_t jentry=0; jentry<50;jentry++) {

    //    if(jentry % 100000 == 0) printf("event %d\n", (Int_t) jentry);
    printf("event %d\n", (Int_t) jentry);
 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);


    if(onia_Gen->Pt() > 990. || muPos_Gen->Pt() > 990. || muNeg_Gen->Pt() > 990.){
      printf("found a generated event w/o proper entries: onia=%1.3f, muPos=%1.3f, muNeg=%1.3f\n", 
	     onia_Gen->Pt(), muPos_Gen->Pt(), muNeg_Gen->Pt());
      continue;
    }

    //===================================================================
    //1.) check all generated (and filtered) events and 
    //fill the corresponding histograms
    //===================================================================

    Double_t etaMuPos_Gen = muPos_Gen->PseudoRapidity();
    Double_t etaMuNeg_Gen = muNeg_Gen->PseudoRapidity();
    Double_t pTMuPos_Gen = muPos_Gen->Pt();
    Double_t pTMuNeg_Gen = muNeg_Gen->Pt();
    Double_t pMuPos_Gen = muPos_Gen->P();
    Double_t pMuNeg_Gen = muNeg_Gen->P();
    Double_t phiMuPos_Gen = muPos_Gen->Phi();
    Double_t phiMuNeg_Gen = muNeg_Gen->Phi();

    Double_t deltaPhi = phiMuNeg_Gen - phiMuPos_Gen;
    if(rejectCowboys){
      if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
      else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
      if(deltaPhi < 0.) //reject cowboys
	continue;
    }

    // // if(!(isMuonInAcceptance(LOOSE, pTMuPos_Gen, etaMuPos_Gen) && isMuonInAcceptance(LOOSE, pTMuNeg_Gen, etaMuNeg_Gen)))
    // if(!(isMuonInAcceptance(TIGHT, pTMuPos_Gen, etaMuPos_Gen) && isMuonInAcceptance(TIGHT, pTMuNeg_Gen, etaMuNeg_Gen)))
    //   continue;
    // // if(pTMuPos_Gen < 2.5 || pTMuNeg_Gen < 2.5)
    // //   continue;
    // if(fabs(etaMuPos_Gen) > 1.6 || fabs(etaMuNeg_Gen) > 1.6)
    //   continue;

    Bool_t decisionPos = kFALSE, decisionNeg = kFALSE;
    //positive muon
    if(TMath::Abs(etaMuPos_Gen)<1.2 && pTMuPos_Gen>4.5) decisionPos=kTRUE;
    if(TMath::Abs(etaMuPos_Gen)>1.2 && TMath::Abs(etaMuPos_Gen)<1.4 && pTMuPos_Gen>3.5) decisionPos=kTRUE;
    if(TMath::Abs(etaMuPos_Gen)>1.4 && TMath::Abs(etaMuPos_Gen)<1.6 && pTMuPos_Gen>3.) decisionPos=kTRUE;
    //negative muon
    if(TMath::Abs(etaMuNeg_Gen)<1.2 && pTMuNeg_Gen>4.5) decisionNeg=kTRUE;
    if(TMath::Abs(etaMuNeg_Gen)>1.2 && TMath::Abs(etaMuNeg_Gen)<1.4 && pTMuNeg_Gen>3.5) decisionNeg=kTRUE;
    if(TMath::Abs(etaMuNeg_Gen)>1.4 && TMath::Abs(etaMuNeg_Gen)<1.6 && pTMuNeg_Gen>3.) decisionNeg=kTRUE;

    if(!decisionPos || !decisionNeg)
      continue;

    Double_t onia_Gen_mass = onia_Gen->M();
    Double_t onia_Gen_pt = onia_Gen->Pt();
    Double_t onia_Gen_P = onia_Gen->P();
    Double_t onia_Gen_eta = onia_Gen->PseudoRapidity();
    Double_t onia_Gen_rap = onia_Gen->Rapidity();
    Double_t onia_Gen_phi = onia_Gen->Phi();
    Double_t onia_Gen_mT = sqrt(onia_Gen_mass*onia_Gen_mass + onia_Gen_pt*onia_Gen_pt);

    if(fabs(onia_Gen_rap) > eff::rapMax)
      continue;

    Int_t rapIndex_Gen = -1;
    for(int iRap = 0; iRap < 2*eff::kNbRapBins; iRap++){
      if(onia_Gen_rap > eff::rapRange[iRap] && onia_Gen_rap < eff::rapRange[iRap+1]){
	rapIndex_Gen = iRap;
	break;
      }
    }
    Int_t rapForPTIndex_Gen = -1;
    for(int iRap = 0; iRap < eff::kNbRapForPTBins; iRap++){
      if(TMath::Abs(onia_Gen_rap) > eff::rapForPTRange[iRap] && 
	 TMath::Abs(onia_Gen_rap) < eff::rapForPTRange[iRap+1]){
	rapForPTIndex_Gen = iRap+1;
	break;
      }
    }
    Int_t pTIndex_Gen = -1;
    for(int iPT = 0; iPT < eff::kNbPTBins[rapForPTIndex_Gen]; iPT++){
      if(onia_Gen_pt > eff::pTRange[rapForPTIndex_Gen][iPT] && onia_Gen_pt < eff::pTRange[rapForPTIndex_Gen][iPT+1]){
	pTIndex_Gen = iPT+1;
	break;
      }
    }
    Int_t rapIntegratedPTIndex_Gen = -1;
    for(int iPT = 0; iPT < eff::kNbPTBins[0]; iPT++){
      if(onia_Gen_pt > eff::pTRange[0][iPT] && onia_Gen_pt < eff::pTRange[0][iPT+1]){
	rapIntegratedPTIndex_Gen = iPT+1;
	break;
      }
    }
    if(rapIndex_Gen < 0){
      // printf("rapIndex_Gen %d, rap(onia) = %f\n", rapIndex_Gen, onia_rap);
      continue;
    }
    if(rapForPTIndex_Gen < 1){
      //printf("rapForPTIndex_Gen %d, rap(onia) = %f\n", rapForPTIndex_Gen, onia_rap);
      continue;
    }
    if(pTIndex_Gen < 1){
      // printf("pTIndex_Gen %d, pT(onia) = %f\n", pTIndex_Gen, onia_pt);
      continue;
    }

    // if(rapForPTIndex_Gen == 1 && (fabs(etaMuPos_Gen) > 0.8 || fabs(etaMuNeg_Gen) > 0.8))
    //   continue; //fiducial cut for rap bin 1

    //==============================
    calcPol(*muPos_Gen, *muNeg_Gen);
    //==============================

    //===============================
    //set up the trigger logic:
    Int_t trigValue = -5;
    if(resonance == JPSI)
    //if(strncmp("HLT_Dimuon10_Jpsi_Barrel", trigLabel, 24) == 0)
    //HLT_Dimuon10_Jpsi_Barrel_v3... 1.4E33
    //HLT_Dimuon10_Jpsi_Barrel_v5... no cowboys
    //HLT_Dimuon10_Jpsi_Barrel_v6... L1DoubleMu0_HighQ
    //HLT_Dimuon13_Jpsi_Barrel_v1... L1DoubleMu0_HighQ
      trigValue = HLT_Dimuon10_Jpsi_Barrel_v3; 
    else if(resonance == UPS1S || resonance == UPS2S || resonance == UPS3S)
      //else if(strncmp("HLT_Dimuon5_Upsilon_Barrel", trigLabel, 26) == 0)
      //HLT_Dimuon5_Upsilon_Barrel_v3... 1.4E33
      //HLT_Dimuon5_Upsilon_Barrel_v5... no cowboys
      //HLT_Dimuon7_Upsilon_Barrel_v1... L1DoubleMu0_HighQ
      //HLT_Dimuon9_Upsilon_Barrel_v1... L1DoubleMu0_HighQ
      trigValue = HLT_Dimuon5_Upsilon_Barrel_v3; 
    else{
      printf("chosen trigger path not a valid option!\n");
      exit(0);
    }
    Double_t etaMuPos, etaMuNeg, pTMuPos, pTMuNeg;
    decisionPos = kFALSE;
    decisionNeg = kFALSE;

    if(onia->Pt() < 990.){
      etaMuPos = muPos->PseudoRapidity();
      etaMuNeg = muNeg->PseudoRapidity();
      pTMuPos = muPos->Pt();
      pTMuNeg = muNeg->Pt();
      //positive muon
      if(TMath::Abs(etaMuPos)<1.2 && pTMuPos>4.5) decisionPos=kTRUE;
      if(TMath::Abs(etaMuPos)>1.2 && TMath::Abs(etaMuPos)<1.4 && pTMuPos>3.5) decisionPos=kTRUE;
      if(TMath::Abs(etaMuPos)>1.4 && TMath::Abs(etaMuPos)<1.6 && pTMuPos>3.) decisionPos=kTRUE;
      //negative muon
      if(TMath::Abs(etaMuNeg)<1.2 && pTMuNeg>4.5) decisionNeg=kTRUE;
      if(TMath::Abs(etaMuNeg)>1.2 && TMath::Abs(etaMuNeg)<1.4 && pTMuNeg>3.5) decisionNeg=kTRUE;
      if(TMath::Abs(etaMuNeg)>1.4 && TMath::Abs(etaMuNeg)<1.6 && pTMuNeg>3.) decisionNeg=kTRUE;
    }

    Bool_t recoPassed = kFALSE, totPassed = kFALSE;
    if(onia->Pt() < 990. && 
       fabs(onia->Rapidity()) < eff::rapMax && 
       JpsiVprob > 0.01 &&
       // isMuonInAcceptance(TIGHT, pTMuPos, etaMuPos) &&
       // isMuonInAcceptance(TIGHT, pTMuNeg, etaMuNeg) &&
       // fabs(muPos->Eta()) < 1.6 &&
       // fabs(muNeg->Eta()) < 1.6
       decisionPos && decisionNeg){      
      recoPassed = kTRUE;
      //if(trigValue == 1) //H: disable trigger condition!
	totPassed = kTRUE;
    }

    if(useSoftMuons > 0){//need to apply the TMOST and "sanity checks"
      //Bool_t tightSelection = kFALSE;
      if(ismuPosTMOneStationTight && ismuNegTMOneStationTight){ //not applied anylonger in the TTree
	if(muPosPglobalMuonHits > 0 && muNegPglobalMuonHits > 0 &&
	   muPosPglobalchi2 < 20. && muNegPglobalchi2 < 20.){
	  recoPassed = kTRUE;
	  totPassed = kTRUE;
	}
	// tightSelection = kTRUE;
      }
      // if(!tightSelection){ //correct the decision in case we work with "tight muons"
      // 	recoPassed = kFALSE;
      // 	totPassed = kFALSE;
      // }
    }
    //    if(onia->Pt() < 990. && fabs(onia->Rapidity()) < eff::rapMax && JpsiVprob > 0.01 && trigValue == 1) totPassed = kTRUE;

    // printf("onia pT=%1.3f, y=%1.3f, muPos pT=%1.3f, y=%1.3f, muNeg pT=%1.3f, y=%1.3f; recoPassed=%d, trigValue=%d, totPassed=%d\n",
    // 	   onia_Gen_pt, fabs(onia_Gen_rap), pTMuPos_Gen, etaMuPos_Gen, pTMuNeg_Gen, etaMuNeg_Gen, recoPassed, trigValue, totPassed);

    recoEff_pT->Fill(recoPassed, onia_Gen_pt);
    recoEff2D_pT_rapNP->Fill(recoPassed, onia_Gen_rap, onia_Gen_pt);
    recoEff2D_pT_rap->Fill(recoPassed, fabs(onia_Gen_rap), onia_Gen_pt);

    totEff_pT->Fill(totPassed, onia_Gen_pt);
    totEff2D_pT_rapNP->Fill(totPassed, onia_Gen_rap, onia_Gen_pt);
    totEff2D_pT_rap->Fill(totPassed, fabs(onia_Gen_rap), onia_Gen_pt);

    // if(onia_Gen_pt > 10. && onia_Gen_pt < 25. && onia_Gen_rap > -1.2 && onia_Gen_rap < 1.2){
    if(onia_Gen_pt > pTMin){
      recoEff_y->Fill(recoPassed, onia_Gen_rap);
      recoEff_phi->Fill(recoPassed, onia_Gen_phi);
      totEff_y->Fill(totPassed, onia_Gen_rap);
      totEff_phi->Fill(totPassed, onia_Gen_phi);
    }

    //fill the eff. histos for all the different frames
    for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){

      //fill some of the histos for the rho factor in the
      //phase space of the real data only
      //if(onia_Gen_pt > 10. && onia_Gen_pt < 25. && onia_Gen_rap > -1.2 && onia_Gen_rap < 1.2){
      if(onia_Gen_pt > pTMin){
      	recoEff_cosTheta[iFrame]->Fill(recoPassed, thisCosTh[iFrame]);
      	recoEff_phiPol[iFrame]->Fill(recoPassed, thisPhi[iFrame]);
      	recoEff2D_cosTheta_phiPol[iFrame]->Fill(recoPassed, thisCosTh[iFrame], thisPhi[iFrame]);
	
      	totEff_cosTheta[iFrame]->Fill(totPassed, thisCosTh[iFrame]);
      	totEff_phiPol[iFrame]->Fill(totPassed, thisPhi[iFrame]);
      	totEff2D_cosTheta_phiPol[iFrame]->Fill(totPassed, thisCosTh[iFrame], thisPhi[iFrame]);

      	recoEff_cosTheta_pT_rap[iFrame][0][0]->Fill(recoPassed, thisCosTh[iFrame]);
      	recoEff_phiPol_pT_rap[iFrame][0][0]->Fill(recoPassed, thisPhi[iFrame]);

     	totEff_cosTheta_pT_rap[iFrame][0][0]->Fill(totPassed, thisCosTh[iFrame]);
      	totEff_phiPol_pT_rap[iFrame][0][0]->Fill(totPassed, thisPhi[iFrame]);
      }

      //histos for neg. and pos. rapidity separately:
      if(rapIndex_Gen >= 0){
	  recoEff2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(recoPassed, thisCosTh[iFrame], thisPhi[iFrame]);
	  totEff2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(totPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(pTIndex_Gen > 0 && rapIndex_Gen >= 0){
	recoEff2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(recoPassed, thisCosTh[iFrame], thisPhi[iFrame]);
	totEff2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(totPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }

      //histos taking together +y and -y
      recoEff2D_pol_pT_rap[iFrame][0][0]->Fill(recoPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      totEff2D_pol_pT_rap[iFrame][0][0]->Fill(totPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      if(rapIntegratedPTIndex_Gen > 0){
	recoEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(recoPassed, thisCosTh[iFrame], thisPhi[iFrame]);
	totEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(totPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(rapForPTIndex_Gen > 0){
	recoEff2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(recoPassed, thisCosTh[iFrame], thisPhi[iFrame]);
	totEff2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(totPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0){
	recoEff2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, thisCosTh[iFrame], thisPhi[iFrame]);
	totEff2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, thisCosTh[iFrame], thisPhi[iFrame]);

      	recoEff_cosTheta_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, thisCosTh[iFrame]);
      	recoEff_phiPol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, thisPhi[iFrame]);

     	totEff_cosTheta_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, thisCosTh[iFrame]);
      	totEff_phiPol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, thisPhi[iFrame]);
      }
	  
      //histos taking together +y and -y and phi 4-folding
      Double_t phiFolded = thisPhi[iFrame];
      Double_t thetaAdjusted = thisCosTh[iFrame];
      if(thisPhi[iFrame] > -90. && thisPhi[iFrame] < 0.)
	phiFolded *= -1;
      else if(thisPhi[iFrame] > 90 && thisPhi[iFrame] < 180){
	phiFolded = 180. - thisPhi[iFrame];
	thetaAdjusted *= -1;
      }
      else if(thisPhi[iFrame] > -180. && thisPhi[iFrame] < -90.){
	phiFolded = 180. + thisPhi[iFrame];
	thetaAdjusted *= -1;
      }
      recoEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(recoPassed, thetaAdjusted, phiFolded);
      totEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(totPassed, thetaAdjusted, phiFolded);
      if(rapIntegratedPTIndex_Gen > 0){
	recoEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(recoPassed, thetaAdjusted, phiFolded);
	totEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(totPassed, thetaAdjusted, phiFolded);
      }
      if(rapForPTIndex_Gen > 0){
	recoEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(recoPassed, thetaAdjusted, phiFolded);
	totEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(totPassed, thetaAdjusted, phiFolded);
      }
      if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0){
	recoEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, thetaAdjusted, phiFolded);
	totEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, thetaAdjusted, phiFolded);
      }
    }

    // //fill the series of correlation histos for lab-frame
    // Double_t deltaEta = fabs(etaMuPos_Gen - etaMuNeg_Gen);
    // Double_t deltaPhi = phiMuPos_Gen - phiMuNeg_Gen;
    // if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
    // else if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
    // Double_t deltaPhiDeg = deltaPhi * 180./TMath::Pi();
    // Double_t deltaR = sqrt(pow(deltaEta,2) + pow(deltaPhi,2));

    // //convert the JpsiDphiM2 variable from rad to deg:
    // Double_t JpsiDetaM2 = sqrt(pow(JpsiDrM2,2) - pow(JpsiDphiM2,2));
    // JpsiDphiM2 *= 180./TMath::Pi();

    // //deltaR
    // recoEff_deltaR_pT_rap[0][0]->Fill(recoPassed, deltaR);
    // totEff_deltaR_pT_rap[0][0]->Fill(totPassed, deltaR);
    // //deltaRM2
    // recoEff_deltaRM2_pT_rap[0][0]->Fill(recoPassed, JpsiDrM2);
    // totEff_deltaRM2_pT_rap[0][0]->Fill(totPassed, JpsiDrM2);
    // //deltaPhiM2
    // recoEff_deltaPhiM2_pT_rap[0][0]->Fill(recoPassed, JpsiDphiM2);
    // totEff_deltaPhiM2_pT_rap[0][0]->Fill(totPassed, JpsiDphiM2);
    // //deltaEtaM2
    // recoEff_deltaEtaM2_pT_rap[0][0]->Fill(recoPassed, JpsiDetaM2);
    // totEff_deltaEtaM2_pT_rap[0][0]->Fill(totPassed, JpsiDetaM2);
    // //distM2
    // recoEff_distM2_pT_rap[0][0]->Fill(recoPassed, JpsiDistM2);
    // totEff_distM2_pT_rap[0][0]->Fill(totPassed, JpsiDistM2);
    // //deltaPhi vs deltaEta
    // recoEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(recoPassed, deltaEta, deltaPhiDeg);
    // totEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(totPassed, deltaEta, deltaPhiDeg);
    // if(rapIntegratedPTIndex_Gen > 0){
    //   //deltaR
    //   recoEff_deltaR_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, deltaR);
    //   totEff_deltaR_pT_rap[pTIndex_Gen][0]->Fill(totPassed, deltaR);
    //   //deltaRM2
    //   recoEff_deltaRM2_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, JpsiDrM2);
    //   totEff_deltaRM2_pT_rap[pTIndex_Gen][0]->Fill(totPassed, JpsiDrM2);
    //   //deltaPhiM2
    //   recoEff_deltaPhiM2_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, JpsiDphiM2);
    //   totEff_deltaPhiM2_pT_rap[pTIndex_Gen][0]->Fill(totPassed, JpsiDphiM2);
    //   //deltaEtaM2
    //   recoEff_deltaEtaM2_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, JpsiDetaM2);
    //   totEff_deltaEtaM2_pT_rap[pTIndex_Gen][0]->Fill(totPassed, JpsiDetaM2);
    //   //distM2
    //   recoEff_distM2_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, JpsiDistM2);
    //   totEff_distM2_pT_rap[pTIndex_Gen][0]->Fill(totPassed, JpsiDistM2);
    //   //deltaPhi vs deltaEta
    //   recoEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, deltaEta, deltaPhiDeg);
    //   totEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][0]->Fill(totPassed, deltaEta, deltaPhiDeg);
    // }
    // if(rapForPTIndex_Gen > 0){
    //   //deltaR
    //   recoEff_deltaR_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, deltaR);
    //   totEff_deltaR_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, deltaR);
    //   //deltaRM2
    //   recoEff_deltaRM2_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDrM2);
    //   totEff_deltaRM2_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, JpsiDrM2);
    //   //deltaPhiM2
    //   recoEff_deltaPhiM2_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDphiM2);
    //   totEff_deltaPhiM2_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, JpsiDphiM2);
    //   //deltaEtaM2
    //   recoEff_deltaEtaM2_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDetaM2);
    //   totEff_deltaEtaM2_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, JpsiDetaM2);
    //   //distM2
    //   recoEff_distM2_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDistM2);
    //   totEff_distM2_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, JpsiDistM2);
    //   //deltaPhi vs deltaEta
    //   recoEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, deltaEta, deltaPhiDeg);
    //   totEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, deltaEta, deltaPhiDeg);
    // }
    // if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0){
    //   //deltaR
    //   recoEff_deltaR_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, deltaR);
    //   totEff_deltaR_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, deltaR);
    //   //deltaRM2
    //   recoEff_deltaRM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDrM2);
    //   totEff_deltaRM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, JpsiDrM2);
    //   //deltaPhiM2
    //   recoEff_deltaPhiM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDphiM2);
    //   totEff_deltaPhiM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, JpsiDphiM2);
    //   //deltaEtaM2
    //   recoEff_deltaEtaM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDetaM2);
    //   totEff_deltaEtaM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, JpsiDetaM2);
    //   //distM2
    //   recoEff_distM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDistM2);
    //   totEff_distM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, JpsiDistM2);
    //   //deltaPhi vs deltaEta
    //   recoEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, deltaEta, deltaPhiDeg);
    //   totEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, deltaEta, deltaPhiDeg);
    // }

    //continue only to process RECO info if
    //there is really a RECO dimuon in the event:
    if(onia->Pt() > 990.)
      continue;
    if(JpsiVprob < 0.01) //vertex probability cut not applied anylonger
      continue;
  
    //take muons only within a certain eta range
    // Double_t etaMuPos = muPos->PseudoRapidity();
    // Double_t etaMuNeg = muNeg->PseudoRapidity();
    // Double_t pTMuPos = muPos->Pt();
    // Double_t pTMuNeg = muNeg->Pt();

    // // // if(!(isMuonInAcceptance(LOOSE, pTMuPos, etaMuPos) && isMuonInAcceptance(LOOSE, pTMuNeg, etaMuNeg)))
    // if(!(isMuonInAcceptance(TIGHT, pTMuPos, etaMuPos) && isMuonInAcceptance(TIGHT, pTMuNeg, etaMuNeg)))
    //   continue;
    // // if(pTMuPos < 2.5 || pTMuNeg < 2.5)
    // //   continue;
    // if(fabs(etaMuPos) > 1.6 || fabs(etaMuNeg) > 1.6)
    //   continue;

    if(!decisionPos || !decisionNeg)
      continue;


    if(fabs(onia->Rapidity()) > eff::rapMax) 
      continue;

    Double_t onia_mass = onia->M();
    Double_t onia_pt = onia->Pt();
    Double_t onia_P = onia->P();
    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();
    Double_t onia_phi = onia->Phi();
    Double_t onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);
    
    Bool_t trigPassed = kFALSE;
    if(trigValue == 1) trigPassed = kTRUE;

    // printf("trigPassed=%d\n", trigValue);

    trigEff_pT->Fill(trigPassed, onia_pt);
    trigEff2D_pT_rapNP->Fill(trigPassed, onia_rap, onia_pt);
    trigEff2D_pT_rap->Fill(trigPassed, fabs(onia_rap), onia_pt);

    //if(onia_Gen_pt > 10. && onia_Gen_pt < 25. && onia_Gen_rap > -1.2 && onia_Gen_rap < 1.2){
    if(onia_Gen_pt > pTMin){
      trigEff_y->Fill(trigPassed, onia_rap);
      trigEff_phi->Fill(trigPassed, onia_phi);
    }

    Int_t rapIndex = -1;
    for(int iRap = 0; iRap < 2*eff::kNbRapBins; iRap++){
      if(onia_rap > eff::rapRange[iRap] && onia_rap < eff::rapRange[iRap+1]){
	rapIndex = iRap;
	break;
      }
    }
    Int_t rapForPTIndex = -1;
    for(int iRap = 0; iRap < eff::kNbRapForPTBins; iRap++){
      if(TMath::Abs(onia_rap) > eff::rapForPTRange[iRap] && 
	 TMath::Abs(onia_rap) < eff::rapForPTRange[iRap+1]){
	rapForPTIndex = iRap+1;
	break;
      }
    }
    Int_t pTIndex = -1;
    for(int iPT = 0; iPT < eff::kNbPTBins[rapForPTIndex]; iPT++){
      if(onia_pt > eff::pTRange[rapForPTIndex][iPT] && onia_pt < eff::pTRange[rapForPTIndex][iPT+1]){
	pTIndex = iPT+1;
	break;
      }
    }
    Int_t rapIntegratedPTIndex = -1;
    for(int iPT = 0; iPT < eff::kNbPTBins[0]; iPT++){
      if(onia_pt > eff::pTRange[0][iPT] && onia_pt < eff::pTRange[0][iPT+1]){
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

    // if(rapForPTIndex == 1 && (fabs(etaMuPos) > 0.8 || fabs(etaMuNeg) > 0.8))
    //   continue; //fiducial cut for rap bin 1

    //fill the eff. histos for all the different frames
    for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){

      //fill some of the histos for the rho factor in the
      //phase space of the real data only
      if(onia_Gen_pt > pTMin){
	//if(onia_Gen_pt > 10. && onia_Gen_pt < 25. && onia_Gen_rap > -1.2 && onia_Gen_rap < 1.2){

      	trigEff_cosTheta[iFrame]->Fill(trigPassed, thisCosTh[iFrame]);
      	trigEff_phiPol[iFrame]->Fill(trigPassed, thisPhi[iFrame]);
      	trigEff2D_cosTheta_phiPol[iFrame]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);	

	trigEff_phiPol_pT_rap[iFrame][0][0]->Fill(trigPassed, thisPhi[iFrame]);
	trigEff_cosTheta_pT_rap[iFrame][0][0]->Fill(trigPassed, thisCosTh[iFrame]);
      }

      //histos for neg. and pos. rapidity separately:
      if(rapIndex >= 0){
	  trigEff2D_pol_pT_rapNP[iFrame][0][rapIndex]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(pTIndex > 0 && rapIndex >= 0){
	trigEff2D_pol_pT_rapNP[iFrame][pTIndex][rapIndex]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
	
      //histos taking together +y and -y
      trigEff2D_pol_pT_rap[iFrame][0][0]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      if(rapIntegratedPTIndex > 0){
	trigEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(rapForPTIndex > 0){
	trigEff2D_pol_pT_rap[iFrame][0][rapForPTIndex]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(pTIndex > 0 && rapForPTIndex > 0){
	trigEff2D_pol_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
	trigEff_phiPol_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(trigPassed, thisPhi[iFrame]);
	trigEff_cosTheta_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(trigPassed, thisCosTh[iFrame]);
      }
	  
      //histos taking together +y and -y and phi 4-folding
      Double_t phiFolded = thisPhi[iFrame];
      Double_t thetaAdjusted = thisCosTh[iFrame];
      if(thisPhi[iFrame] > -90. && thisPhi[iFrame] < 0.)
	phiFolded *= -1;
      else if(thisPhi[iFrame] > 90 && thisPhi[iFrame] < 180){
	phiFolded = 180. - thisPhi[iFrame];
	thetaAdjusted *= -1;
      }
      else if(thisPhi[iFrame] > -180. && thisPhi[iFrame] < -90.){
	phiFolded = 180. + thisPhi[iFrame];
	thetaAdjusted *= -1;
      }
      trigEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(trigPassed, thetaAdjusted, phiFolded);
      if(rapIntegratedPTIndex > 0){
	trigEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex][0]->Fill(trigPassed, thetaAdjusted, phiFolded);
      }
      if(rapForPTIndex > 0){
	trigEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex]->Fill(trigPassed, thetaAdjusted, phiFolded);
      }
      if(pTIndex > 0 && rapForPTIndex > 0){
	trigEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex][rapForPTIndex]->Fill(trigPassed, thetaAdjusted, phiFolded);
      }
    }

    // //deltaR
    // trigEff_deltaR_pT_rap[0][0]->Fill(trigPassed, deltaR);
    // //deltaRM2
    // trigEff_deltaRM2_pT_rap[0][0]->Fill(trigPassed, JpsiDrM2);
    // //deltaPhiM2
    // trigEff_deltaPhiM2_pT_rap[0][0]->Fill(trigPassed, JpsiDphiM2);
    // //deltaEtaM2
    // trigEff_deltaEtaM2_pT_rap[0][0]->Fill(trigPassed, JpsiDetaM2);
    // //distM2
    // trigEff_distM2_pT_rap[0][0]->Fill(trigPassed, JpsiDistM2);
    // //deltaPhi vs deltaEta
    // trigEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(trigPassed, deltaEta, deltaPhiDeg);
    // if(rapIntegratedPTIndex > 0){
    //   //deltaR
    //   trigEff_deltaR_pT_rap[pTIndex][0]->Fill(trigPassed, deltaR);
    //   //deltaRM2
    //   trigEff_deltaRM2_pT_rap[pTIndex][0]->Fill(trigPassed, JpsiDrM2);
    //   //deltaPhiM2
    //   trigEff_deltaPhiM2_pT_rap[pTIndex][0]->Fill(trigPassed, JpsiDphiM2);
    //   //deltaEtaM2
    //   trigEff_deltaEtaM2_pT_rap[pTIndex][0]->Fill(trigPassed, JpsiDetaM2);
    //   //distM2
    //   trigEff_distM2_pT_rap[pTIndex][0]->Fill(trigPassed, JpsiDistM2);
    //   //deltaPhi vs deltaEta
    //   trigEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex][0]->Fill(trigPassed, deltaEta, deltaPhiDeg);
    // }
    // if(rapForPTIndex > 0){
    //   //deltaR
    //   trigEff_deltaR_pT_rap[0][rapForPTIndex]->Fill(trigPassed, deltaR);
    //   //deltaRM2
    //   trigEff_deltaRM2_pT_rap[0][rapForPTIndex]->Fill(trigPassed, JpsiDrM2);
    //   //deltaPhiM2
    //   trigEff_deltaPhiM2_pT_rap[0][rapForPTIndex]->Fill(trigPassed, JpsiDphiM2);
    //   //deltaEtaM2
    //   trigEff_deltaEtaM2_pT_rap[0][rapForPTIndex]->Fill(trigPassed, JpsiDetaM2);
    //   //distM2
    //   trigEff_distM2_pT_rap[0][rapForPTIndex]->Fill(trigPassed, JpsiDistM2);
    //   //deltaPhi vs deltaEta
    //   trigEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex]->Fill(trigPassed, deltaEta, deltaPhiDeg);
    // }
    // if(pTIndex > 0 && rapForPTIndex > 0){
    //   //deltaR
    //   trigEff_deltaR_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, deltaR);
    //   //deltaRM2
    //   trigEff_deltaRM2_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, JpsiDrM2);
    //   //deltaPhiM2
    //   trigEff_deltaPhiM2_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, JpsiDphiM2);
    //   //deltaEtaM2
    //   trigEff_deltaEtaM2_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, JpsiDetaM2);
    //   //distM2
    //   trigEff_distM2_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, JpsiDistM2);
    //   //deltaPhi vs deltaEta
    //   trigEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, deltaEta, deltaPhiDeg);
    // }
  }//loop over entries

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}

