#define MCTruthEff_cxx
#include "MCTruthEff.h"
#include "calcPol.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>
#include "TEfficiency.h"

//=================================================
//1D histos versus J/psi pT, phi and y
//=================================================
//histos for neg. and pos. rapidity separately:
TEfficiency *recoEff_pT, *recoEff_y, *recoEff_phi;
TEfficiency *trigEff_pT, *trigEff_y, *trigEff_phi;
TEfficiency *totEff_pT, *totEff_y, *totEff_phi;
//=================================================
//2D histos versus J/psi pT and y
//=================================================
//histos for neg. and pos. rapidity separately:
TEfficiency *recoEff2D_pT_rapNP, *recoEff2D_pT_rap;
TEfficiency *trigEff2D_pT_rapNP, *trigEff2D_pT_rap;
TEfficiency *totEff2D_pT_rapNP, *totEff2D_pT_rap;
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

Double_t massMuOnia;
//==============================================
void MCTruthEff::Loop(Int_t selDimuType, Char_t *trigLabel)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nb = 0;
  Long64_t countRecEvent = 0;
  printf("number of entries = %d\n", (Int_t) nentries);

  //loop over the events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if(jentry % 100000 == 0) printf("event %d\n", (Int_t) jentry);

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

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

    //take muons only within a certain eta range
    if((fabs(etaMuPos_Gen) < eff::etaPS[0] && pTMuPos_Gen < eff::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuPos_Gen) > eff::etaPS[0] && fabs(etaMuPos_Gen) < eff::etaPS[1] && pMuPos_Gen < eff::pMuMin[1]) ||
       (fabs(etaMuPos_Gen) > eff::etaPS[1] && fabs(etaMuPos_Gen) < eff::etaPS[2] && pTMuPos_Gen < eff::pTMuMin[2]))
      continue;
    //(b) on the negative muon
    if((fabs(etaMuNeg_Gen) < eff::etaPS[0] && pTMuNeg_Gen < eff::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuNeg_Gen) > eff::etaPS[0] && fabs(etaMuNeg_Gen) < eff::etaPS[1] && pMuNeg_Gen < eff::pMuMin[1]) ||
       (fabs(etaMuNeg_Gen) > eff::etaPS[1] && fabs(etaMuNeg_Gen) < eff::etaPS[2] && pTMuNeg_Gen < eff::pTMuMin[2]))
      continue;

    Double_t onia_Gen_mass = onia_Gen->M();
    Double_t onia_Gen_pt = onia_Gen->Pt();
    Double_t onia_Gen_P = onia_Gen->P();
    Double_t onia_Gen_eta = onia_Gen->PseudoRapidity();
    Double_t onia_Gen_rap = onia_Gen->Rapidity();
    Double_t onia_Gen_phi = onia_Gen->Phi();
    Double_t onia_Gen_mT = sqrt(onia_Gen_mass*onia_Gen_mass + onia_Gen_pt*onia_Gen_pt);
    
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
      // printf("rapForPTIndex_Gen %d, rap(onia) = %f\n", rapForPTIndex_Gen, onia_rap);
      continue;
    }
    if(pTIndex_Gen < 1){
      // printf("pTIndex_Gen %d, pT(onia) = %f\n", pTIndex_Gen, onia_pt);
      continue;
    }

    //==============================
    calcPol(*muPos_Gen, *muNeg_Gen);
    //==============================

    Int_t trigValue;
    //process events further that pass the trigger under study
    if(strncmp("HLT_Mu0_TkMu0_OST_Jpsi", trigLabel, 22) == 0)
      trigValue = HLT_Mu0_TkMu0_OST_Jpsi;
    else if(strncmp("HLT_DoubleMu0", trigLabel, 13) == 0)
      trigValue = HLT_DoubleMu0;
    else{
      printf("chosen trigger path, %s, not a valid option!\n", trigLabel);
      exit(0);
    }

    Bool_t dimuTypeFlag = kTRUE;
    if(selDimuType < 3 && JpsiType != selDimuType)
      dimuTypeFlag = kFALSE;
    else if(selDimuType == 3 && JpsiType > 1) //only GG or GT
      dimuTypeFlag = kFALSE;

    Bool_t recoPassed = kFALSE, totPassed = kFALSE;
    if(onia->Pt() < 990. && dimuTypeFlag) recoPassed = kTRUE;
    if(onia->Pt() < 990. && dimuTypeFlag && trigValue == 1) totPassed = kTRUE;
    
    recoEff_pT->Fill(recoPassed, onia_Gen_pt);
    recoEff_y->Fill(recoPassed, onia_Gen_rap);
    recoEff_phi->Fill(recoPassed, onia_Gen_phi);
    recoEff2D_pT_rapNP->Fill(recoPassed, onia_Gen_rap, onia_Gen_pt);
    recoEff2D_pT_rap->Fill(recoPassed, fabs(onia_Gen_rap), onia_Gen_pt);

    totEff_pT->Fill(totPassed, onia_Gen_pt);
    totEff_y->Fill(totPassed, onia_Gen_rap);
    totEff_phi->Fill(totPassed, onia_Gen_phi);
    totEff2D_pT_rapNP->Fill(totPassed, onia_Gen_rap, onia_Gen_pt);
    totEff2D_pT_rap->Fill(totPassed, fabs(onia_Gen_rap), onia_Gen_pt);

    //fill the eff. histos for all the different frames
    for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){

	//histos for neg. and pos. rapidity separately:
      if(rapIndex_Gen >= 0){
	  recoEff2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(recoPassed, thisCosTh[iFrame], thisPhi[iFrame]);
	  totEff2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(recoPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(pTIndex_Gen > 0 && rapIndex_Gen >= 0){
	recoEff2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(recoPassed, thisCosTh[iFrame], thisPhi[iFrame]);
	totEff2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(recoPassed, thisCosTh[iFrame], thisPhi[iFrame]);
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

    //continue only to process RECO info if
    //there is really a RECO dimuon in the event:
    if(onia->Pt() > 990.)
      continue;
    if(selDimuType < 3 && JpsiType != selDimuType)
      continue;
    else if(selDimuType == 3 && JpsiType > 1) //only GG or GT
      continue;

    Bool_t trigPassed = kFALSE;
    if(trigValue == 1) trigPassed = kTRUE;

    trigEff_pT->Fill(trigPassed, onia_Gen_pt);
    trigEff_y->Fill(trigPassed, onia_Gen_rap);
    trigEff_phi->Fill(trigPassed, onia_Gen_phi);
    trigEff2D_pT_rapNP->Fill(trigPassed, onia_Gen_rap, onia_Gen_pt);
    trigEff2D_pT_rap->Fill(trigPassed, fabs(onia_Gen_rap), onia_Gen_pt);

    //fill the eff. histos for all the different frames
    for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){

	//histos for neg. and pos. rapidity separately:
      if(rapIndex_Gen >= 0){
	  trigEff2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(pTIndex_Gen > 0 && rapIndex_Gen >= 0){
	trigEff2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
	
      //histos taking together +y and -y
      trigEff2D_pol_pT_rap[iFrame][0][0]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      if(rapIntegratedPTIndex_Gen > 0){
	trigEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(rapForPTIndex_Gen > 0){
	trigEff2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0){
	trigEff2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(trigPassed, thisCosTh[iFrame], thisPhi[iFrame]);
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
      if(rapIntegratedPTIndex_Gen > 0){
	trigEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(trigPassed, thetaAdjusted, phiFolded);
      }
      if(rapForPTIndex_Gen > 0){
	trigEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(trigPassed, thetaAdjusted, phiFolded);
      }
      if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0){
	trigEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(trigPassed, thetaAdjusted, phiFolded);
      }
    }
  }//loop over entries

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}

