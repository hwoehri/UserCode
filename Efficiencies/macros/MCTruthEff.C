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
//=================================================
//deltaR, deltaPhi vs deltaEta for various J/psi pT and y
//=================================================
TEfficiency *recoEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TEfficiency *recoEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TEfficiency *recoEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TEfficiency *recoEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TEfficiency *recoEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TEfficiency *recoEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
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
//   for (Long64_t jentry=0; jentry<10000;jentry++) {

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
    Double_t phiMuPos_Gen = muPos_Gen->Phi();
    Double_t phiMuNeg_Gen = muNeg_Gen->Phi();

    //take muons only within a certain eta range
    //(a) positive muon
    if((fabs(etaMuPos_Gen) < eff::etaPS[0] && pTMuPos_Gen < eff::pTMuMin[0]) || //mid-rapidity cut
       // (fabs(etaMuPos_Gen) > eff::etaPS[0] && fabs(etaMuPos_Gen) < eff::etaPS[1] && pMuPos_Gen < eff::pMuMin[1]) ||
       (fabs(etaMuPos_Gen) > eff::etaPS[0] && fabs(etaMuPos_Gen) < eff::etaPS[1] && pTMuPos_Gen < eff::pTMuMin[1]) ||
       (fabs(etaMuPos_Gen) > eff::etaPS[1] && fabs(etaMuPos_Gen) < eff::etaPS[2] && pTMuPos_Gen < eff::pTMuMin[2]))
      continue;
    //(b) negative muon
    if((fabs(etaMuNeg_Gen) < eff::etaPS[0] && pTMuNeg_Gen < eff::pTMuMin[0]) || //mid-rapidity cut
       // (fabs(etaMuNeg_Gen) > eff::etaPS[0] && fabs(etaMuNeg_Gen) < eff::etaPS[1] && pMuNeg_Gen < eff::pMuMin[1]) ||
       (fabs(etaMuNeg_Gen) > eff::etaPS[0] && fabs(etaMuNeg_Gen) < eff::etaPS[1] && pTMuNeg_Gen < eff::pTMuMin[1]) ||
       (fabs(etaMuNeg_Gen) > eff::etaPS[1] && fabs(etaMuNeg_Gen) < eff::etaPS[2] && pTMuNeg_Gen < eff::pTMuMin[2]))
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

    //===============================
    //set up the trigger logic:
    Int_t trigValue = -5;
    if(strncmp("HLT_Mu0_TkMu0_OST_Jpsi", trigLabel, 22) == 0){
      trigValue = HLT_Mu0_TkMu0_OST_Jpsi; //0... not matched; 3... fired, but at least one not matched; 1, 2, -2 fired+both matched
      //accept the event only if the higher pT muon is the HLT muon
      if(trigValue == 1){
	if(pTMuPos_Gen < pTMuNeg_Gen)
	  trigValue = 0;
	else
	  trigValue = 1;
      }
      else if(trigValue == -1){
	if(pTMuNeg_Gen < pTMuNeg_Gen)
	  trigValue = 0;
	else
	  trigValue = 1;
      }
      else if(trigValue == 2)
	trigValue = 1;
    }
    else if(strncmp("HLT_DoubleMu0", trigLabel, 13) == 0)
      trigValue = HLT_DoubleMu0; //0... not matched, 1... fired+matched, 3... only fired
    else{
      printf("chosen trigger path, %s, not a valid option!\n", trigLabel);
      exit(0);
    }
    // //trigValue must NOT be 0 and not 3 --> set to 1 if this is the case
    // if(trigValue == 1 || trigValue == -1 || trigValue == 2)
    //   trigValue = 1;
    // else
    //   trigValue = 0;

    // ---> For the asymmetric triggers:
    // 0 : event not firing the corresponding trigger
    // 1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, POSITIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
    // -1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, NEGATIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
    // 2 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, both matched to the tighter HLT object (usually a L3 muon)
    // 3 : event firing the corresponding trigger but at least one RECO muon is not matched to the HLT objects 

    Bool_t dimuTypeFlag = kTRUE;
    if(selDimuType < 3 && JpsiType != selDimuType)
      dimuTypeFlag = kFALSE;
    else if(selDimuType == 3 && JpsiType > 1) //only GG or GT
      dimuTypeFlag = kFALSE;

    Bool_t recoPassed = kFALSE, totPassed = kFALSE;
    if(onia->Pt() < 990. && fabs(onia->Rapidity()) < eff::rapMax && dimuTypeFlag && JpsiVprob > 0.01) recoPassed = kTRUE;
    if(onia->Pt() < 990. && fabs(onia->Rapidity()) < eff::rapMax && dimuTypeFlag && JpsiVprob > 0.01 && trigValue == 1) totPassed = kTRUE;
    
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

    //fill the series of correlation histos for lab-frame
    Double_t deltaEta = fabs(etaMuPos_Gen - etaMuNeg_Gen);
    Double_t deltaPhi = phiMuPos_Gen - phiMuNeg_Gen;
    if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
    else if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
    Double_t deltaPhiDeg = deltaPhi * 180./TMath::Pi();
    Double_t deltaR = sqrt(pow(deltaEta,2) + pow(deltaPhi,2));

    //convert the JpsiDphiM2 variable from rad to deg:
    Double_t JpsiDetaM2 = sqrt(pow(JpsiDrM2,2) - pow(JpsiDphiM2,2));
    JpsiDphiM2 *= 180./TMath::Pi();

    //deltaR
    recoEff_deltaR_pT_rap[0][0]->Fill(recoPassed, deltaR);
    totEff_deltaR_pT_rap[0][0]->Fill(totPassed, deltaR);
    //deltaRM2
    recoEff_deltaRM2_pT_rap[0][0]->Fill(recoPassed, JpsiDrM2);
    totEff_deltaRM2_pT_rap[0][0]->Fill(totPassed, JpsiDrM2);
    //deltaPhiM2
    recoEff_deltaPhiM2_pT_rap[0][0]->Fill(recoPassed, JpsiDphiM2);
    totEff_deltaPhiM2_pT_rap[0][0]->Fill(totPassed, JpsiDphiM2);
    //deltaEtaM2
    recoEff_deltaEtaM2_pT_rap[0][0]->Fill(recoPassed, JpsiDetaM2);
    totEff_deltaEtaM2_pT_rap[0][0]->Fill(totPassed, JpsiDetaM2);
    //distM2
    recoEff_distM2_pT_rap[0][0]->Fill(recoPassed, JpsiDistM2);
    totEff_distM2_pT_rap[0][0]->Fill(totPassed, JpsiDistM2);
    //deltaPhi vs deltaEta
    recoEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(recoPassed, deltaEta, deltaPhiDeg);
    totEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(totPassed, deltaEta, deltaPhiDeg);
    if(rapIntegratedPTIndex_Gen > 0){
      //deltaR
      recoEff_deltaR_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, deltaR);
      totEff_deltaR_pT_rap[pTIndex_Gen][0]->Fill(totPassed, deltaR);
      //deltaRM2
      recoEff_deltaRM2_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, JpsiDrM2);
      totEff_deltaRM2_pT_rap[pTIndex_Gen][0]->Fill(totPassed, JpsiDrM2);
      //deltaPhiM2
      recoEff_deltaPhiM2_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, JpsiDphiM2);
      totEff_deltaPhiM2_pT_rap[pTIndex_Gen][0]->Fill(totPassed, JpsiDphiM2);
      //deltaEtaM2
      recoEff_deltaEtaM2_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, JpsiDetaM2);
      totEff_deltaEtaM2_pT_rap[pTIndex_Gen][0]->Fill(totPassed, JpsiDetaM2);
      //distM2
      recoEff_distM2_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, JpsiDistM2);
      totEff_distM2_pT_rap[pTIndex_Gen][0]->Fill(totPassed, JpsiDistM2);
      //deltaPhi vs deltaEta
      recoEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][0]->Fill(recoPassed, deltaEta, deltaPhiDeg);
      totEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][0]->Fill(totPassed, deltaEta, deltaPhiDeg);
    }
    if(rapForPTIndex_Gen > 0){
      //deltaR
      recoEff_deltaR_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, deltaR);
      totEff_deltaR_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, deltaR);
      //deltaRM2
      recoEff_deltaRM2_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDrM2);
      totEff_deltaRM2_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, JpsiDrM2);
      //deltaPhiM2
      recoEff_deltaPhiM2_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDphiM2);
      totEff_deltaPhiM2_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, JpsiDphiM2);
      //deltaEtaM2
      recoEff_deltaEtaM2_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDetaM2);
      totEff_deltaEtaM2_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, JpsiDetaM2);
      //distM2
      recoEff_distM2_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDistM2);
      totEff_distM2_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, JpsiDistM2);
      //deltaPhi vs deltaEta
      recoEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex_Gen]->Fill(recoPassed, deltaEta, deltaPhiDeg);
      totEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex_Gen]->Fill(totPassed, deltaEta, deltaPhiDeg);
    }
    if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0){
      //deltaR
      recoEff_deltaR_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, deltaR);
      totEff_deltaR_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, deltaR);
      //deltaRM2
      recoEff_deltaRM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDrM2);
      totEff_deltaRM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, JpsiDrM2);
      //deltaPhiM2
      recoEff_deltaPhiM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDphiM2);
      totEff_deltaPhiM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, JpsiDphiM2);
      //deltaEtaM2
      recoEff_deltaEtaM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDetaM2);
      totEff_deltaEtaM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, JpsiDetaM2);
      //distM2
      recoEff_distM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, JpsiDistM2);
      totEff_distM2_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, JpsiDistM2);
      //deltaPhi vs deltaEta
      recoEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(recoPassed, deltaEta, deltaPhiDeg);
      totEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(totPassed, deltaEta, deltaPhiDeg);
    }

    //continue only to process RECO info if
    //there is really a RECO dimuon in the event:
    if(onia->Pt() > 990.)
      continue;
    if(selDimuType < 3 && JpsiType != selDimuType)
      continue;
    else if(selDimuType == 3 && JpsiType > 1) //only GG or GT
      continue;
    if(JpsiVprob < 0.01) //vertex probability cut not applied anylonger
      continue;
  
    //take muons only within a certain eta range
    Double_t etaMuPos = muPos->PseudoRapidity();
    Double_t etaMuNeg = muNeg->PseudoRapidity();
    Double_t pTMuPos = muPos->Pt();
    Double_t pTMuNeg = muNeg->Pt();
    Double_t pMuPos = muPos->P();
    Double_t pMuNeg = muNeg->P();
    //(a) positive muon
    if((fabs(etaMuPos) < eff::etaPS[0] && pTMuPos < eff::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuPos) > eff::etaPS[0] && fabs(etaMuPos) < eff::etaPS[1] && pMuPos < eff::pMuMin[1]) ||
       (fabs(etaMuPos) > eff::etaPS[1] && fabs(etaMuPos) < eff::etaPS[2] && pTMuPos < eff::pTMuMin[2]))
      continue;
    //(b) negative muon
    if((fabs(etaMuNeg) < eff::etaPS[0] && pTMuNeg < eff::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuNeg) > eff::etaPS[0] && fabs(etaMuNeg) < eff::etaPS[1] && pMuNeg < eff::pMuMin[1]) ||
       (fabs(etaMuNeg) > eff::etaPS[1] && fabs(etaMuNeg) < eff::etaPS[2] && pTMuNeg < eff::pTMuMin[2]))
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

    trigEff_pT->Fill(trigPassed, onia_pt);
    trigEff_y->Fill(trigPassed, onia_rap);
    trigEff_phi->Fill(trigPassed, onia_phi);
    trigEff2D_pT_rapNP->Fill(trigPassed, onia_rap, onia_pt);
    trigEff2D_pT_rap->Fill(trigPassed, fabs(onia_rap), onia_pt);

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

    //fill the eff. histos for all the different frames
    for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){

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

    //deltaR
    trigEff_deltaR_pT_rap[0][0]->Fill(trigPassed, deltaR);
    //deltaRM2
    trigEff_deltaRM2_pT_rap[0][0]->Fill(trigPassed, JpsiDrM2);
    //deltaPhiM2
    trigEff_deltaPhiM2_pT_rap[0][0]->Fill(trigPassed, JpsiDphiM2);
    //deltaEtaM2
    trigEff_deltaEtaM2_pT_rap[0][0]->Fill(trigPassed, JpsiDetaM2);
    //distM2
    trigEff_distM2_pT_rap[0][0]->Fill(trigPassed, JpsiDistM2);
    //deltaPhi vs deltaEta
    trigEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(trigPassed, deltaEta, deltaPhiDeg);
    if(rapIntegratedPTIndex > 0){
      //deltaR
      trigEff_deltaR_pT_rap[pTIndex][0]->Fill(trigPassed, deltaR);
      //deltaRM2
      trigEff_deltaRM2_pT_rap[pTIndex][0]->Fill(trigPassed, JpsiDrM2);
      //deltaPhiM2
      trigEff_deltaPhiM2_pT_rap[pTIndex][0]->Fill(trigPassed, JpsiDphiM2);
      //deltaEtaM2
      trigEff_deltaEtaM2_pT_rap[pTIndex][0]->Fill(trigPassed, JpsiDetaM2);
      //distM2
      trigEff_distM2_pT_rap[pTIndex][0]->Fill(trigPassed, JpsiDistM2);
      //deltaPhi vs deltaEta
      trigEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex][0]->Fill(trigPassed, deltaEta, deltaPhiDeg);
    }
    if(rapForPTIndex > 0){
      //deltaR
      trigEff_deltaR_pT_rap[0][rapForPTIndex]->Fill(trigPassed, deltaR);
      //deltaRM2
      trigEff_deltaRM2_pT_rap[0][rapForPTIndex]->Fill(trigPassed, JpsiDrM2);
      //deltaPhiM2
      trigEff_deltaPhiM2_pT_rap[0][rapForPTIndex]->Fill(trigPassed, JpsiDphiM2);
      //deltaEtaM2
      trigEff_deltaEtaM2_pT_rap[0][rapForPTIndex]->Fill(trigPassed, JpsiDetaM2);
      //distM2
      trigEff_distM2_pT_rap[0][rapForPTIndex]->Fill(trigPassed, JpsiDistM2);
      //deltaPhi vs deltaEta
      trigEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex]->Fill(trigPassed, deltaEta, deltaPhiDeg);
    }
    if(pTIndex > 0 && rapForPTIndex > 0){
      //deltaR
      trigEff_deltaR_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, deltaR);
      //deltaRM2
      trigEff_deltaRM2_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, JpsiDrM2);
      //deltaPhiM2
      trigEff_deltaPhiM2_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, JpsiDphiM2);
      //deltaEtaM2
      trigEff_deltaEtaM2_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, JpsiDetaM2);
      //distM2
      trigEff_distM2_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, JpsiDistM2);
      //deltaPhi vs deltaEta
      trigEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex][rapForPTIndex]->Fill(trigPassed, deltaEta, deltaPhiDeg);
    }
  }//loop over entries

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}

