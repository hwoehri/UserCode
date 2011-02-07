#define RecoEff_cxx
#include "RecoEff.h"
#include "calcPol.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>

//some statistics
//TH1F *Reco_StatEv;

//histos at reco level for mu+, mu-, gamma, Onia
// TH1D *Reco_mupl_pt[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
// TH1D *Reco_mumi_pt[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
// TH1D *Reco_mupl_eta[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
// TH1D *Reco_mumi_eta[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
// TH1D *Reco_mupl_phi[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
// TH1D *Reco_mumi_phi[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1D *Reco_Onia_mass[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1D *Reco_Onia_pt[jpsi::kNbRapForPTBins+1];
TH1D *Reco_Onia_phi[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2D *Reco_Onia_rap_pT;
//debugging histos:
// TH2D *hPhiPos_PhiNeg[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
// TH2D *hPtPos_PtNeg[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
// TH2D *hEtaPos_EtaNeg[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
// TH1D *hDeltaPhi[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];

//================================================
//2D polarization histos, differential in pT and y
//================================================
//histos for neg. and pos. rapidity separately:
TH2D *hGen2D_pol_pT_rapNP[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2D *Reco2D_pol_pT_rapNP[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
//histos taking together +y and -y:
TH2D *hGen2D_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2D *Reco2D_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
//histos taking together +y and -y and 4-folding in phi
TH2D *hGen2D_pol_pT_rap_phiFolded[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2D *Reco2D_pol_pT_rap_phiFolded[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];


// TH1F *Reco_Onia_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1][jpsi::kNbPolVar];
// //same pol histos, but split into FWD and BWD rapidity:
// TH1F *Reco_Onia_pol_pT_rap_FWDBWD[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][2*jpsi::kNbRapBins][jpsi::kNbPolVar];

//
// TH1F *hDelta[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
// TH1F *hSin2Delta[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];

Double_t massMuOnia;

Double_t CalcPolWeight(Double_t thisCosTh);
//==============================================
void RecoEff::Loop(Int_t selDimuType)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t cutAtRecEvent = nentries;
  Long64_t countRecEvent = 0;
  Long64_t nb = 0;
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

    Bool_t genMuonsInPS = kTRUE;
    Double_t etaMuPos_Gen = muPos_Gen->PseudoRapidity();
    Double_t etaMuNeg_Gen = muNeg_Gen->PseudoRapidity();
    Double_t pTMuPos_Gen = muPos_Gen->Pt();
    Double_t pTMuNeg_Gen = muNeg_Gen->Pt();
    Double_t pMuPos_Gen = muPos_Gen->P();
    Double_t pMuNeg_Gen = muNeg_Gen->P();

    //take muons only within a certain eta range
    if((fabs(etaMuPos_Gen) < jpsi::etaPS[0] && pTMuPos_Gen < jpsi::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuPos_Gen) > jpsi::etaPS[0] && fabs(etaMuPos_Gen) < jpsi::etaPS[1] && pMuPos_Gen < jpsi::pMuMin[1]) ||
       (fabs(etaMuPos_Gen) > jpsi::etaPS[1] && fabs(etaMuPos_Gen) < jpsi::etaPS[2] && pTMuPos_Gen < jpsi::pTMuMin[2]))
      //      genMuonsInPS = kFALSE;
      continue;
    //(b) on the negative muon
    if((fabs(etaMuNeg_Gen) < jpsi::etaPS[0] && pTMuNeg_Gen < jpsi::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuNeg_Gen) > jpsi::etaPS[0] && fabs(etaMuNeg_Gen) < jpsi::etaPS[1] && pMuNeg_Gen < jpsi::pMuMin[1]) ||
       (fabs(etaMuNeg_Gen) > jpsi::etaPS[1] && fabs(etaMuNeg_Gen) < jpsi::etaPS[2] && pTMuNeg_Gen < jpsi::pTMuMin[2]))
//       genMuonsInPS = kFALSE;
      continue;

//     if(genMuonsInPS){ //fill the generator level info only if the muons are emitted in the PS

      Double_t onia_Gen_mass = onia_Gen->M();
      Double_t onia_Gen_pt = onia_Gen->Pt();
      Double_t onia_Gen_P = onia_Gen->P();
      Double_t onia_Gen_eta = onia_Gen->PseudoRapidity();
      Double_t onia_Gen_rap = onia_Gen->Rapidity();
      Double_t onia_Gen_phi = onia_Gen->Phi();
      Double_t onia_Gen_mT = sqrt(onia_Gen_mass*onia_Gen_mass + onia_Gen_pt*onia_Gen_pt);

      Int_t rapIndex_Gen = -1;
      for(int iRap = 0; iRap < 2*jpsi::kNbRapBins; iRap++){
	if(onia_Gen_rap > jpsi::rapRange[iRap] && onia_Gen_rap < jpsi::rapRange[iRap+1]){
	  rapIndex_Gen = iRap;
	  break;
	}
      }
      Int_t rapForPTIndex_Gen = -1;
      for(int iRap = 0; iRap < jpsi::kNbRapForPTBins; iRap++){
	if(TMath::Abs(onia_Gen_rap) > jpsi::rapForPTRange[iRap] && 
	   TMath::Abs(onia_Gen_rap) < jpsi::rapForPTRange[iRap+1]){
	  rapForPTIndex_Gen = iRap+1;
	  break;
	}
      }
      Int_t pTIndex_Gen = -1;
      for(int iPT = 0; iPT < jpsi::kNbPTBins[rapForPTIndex_Gen]; iPT++){
	if(onia_Gen_pt > jpsi::pTRange[rapForPTIndex_Gen][iPT] && onia_Gen_pt < jpsi::pTRange[rapForPTIndex_Gen][iPT+1]){
	  pTIndex_Gen = iPT+1;
	  break;
	}
      }
      Int_t rapIntegratedPTIndex_Gen = -1;
      for(int iPT = 0; iPT < jpsi::kNbPTBins[0]; iPT++){
	if(onia_Gen_pt > jpsi::pTRange[0][iPT] && onia_Gen_pt < jpsi::pTRange[0][iPT+1]){
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
    
      //=====================
      calcPol(*muPos_Gen, *muNeg_Gen);
      //=====================

      //fill the generator level histograms
      for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){

// 	Double_t weight = CalcPolWeight(thisCosTh[iFrame]);
	Double_t weight = 1.;

	//histos for neg. and pos. rapidity separately:
	if(rapIndex_Gen >= 0)
	  hGen2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(pTIndex_Gen > 0 && rapIndex_Gen >= 0)
	  hGen2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	
	//histos taking together +y and -y
	hGen2D_pol_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapIntegratedPTIndex_Gen > 0)
	  hGen2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapForPTIndex_Gen > 0)
	  hGen2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0)
	  hGen2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	  
	//histos taking together +y and -y and phi-4-folding
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
	hGen2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(thetaAdjusted, phiFolded, weight);
	if(rapIntegratedPTIndex_Gen > 0)
	  hGen2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thetaAdjusted, phiFolded, weight);
	if(rapForPTIndex_Gen > 0)
	  hGen2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded, weight);
	if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0)
	  hGen2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded, weight);
      }
//     }
    //===================================================================
    //end of Generator level block
    //===================================================================

    //continue only to process RECO info if
    //there is really a RECO dimuon in the event:
    if(onia->Pt() > 990.) 
      continue;

    //reject processing of events where the dimuon type (GG, GT or TT)
    //does not correspond to the chosen one
    if(selDimuType < 3 && JpsiType != selDimuType)
      continue;
    else if(selDimuType == 3 && JpsiType > 1) //only GG or GT
      continue;

//     Reco_StatEv->Fill(2.5);//count all events

    Double_t etaMuPos = muPos->PseudoRapidity();
    Double_t etaMuNeg = muNeg->PseudoRapidity();
    Double_t pTMuPos = muPos->Pt();
    Double_t pTMuNeg = muNeg->Pt();
    Double_t pMuPos = muPos->P();
    Double_t pMuNeg = muNeg->P();

    //take muons only within a certain eta range
    if((fabs(etaMuPos) < jpsi::etaPS[0] && pTMuPos < jpsi::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuPos) > jpsi::etaPS[0] && fabs(etaMuPos) < jpsi::etaPS[1] && pMuPos < jpsi::pMuMin[1]) ||
       (fabs(etaMuPos) > jpsi::etaPS[1] && fabs(etaMuPos) < jpsi::etaPS[2] && pTMuPos < jpsi::pTMuMin[2]))
      continue;
    //(b) on the negative muon
    if((fabs(etaMuNeg) < jpsi::etaPS[0] && pTMuNeg < jpsi::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuNeg) > jpsi::etaPS[0] && fabs(etaMuNeg) < jpsi::etaPS[1] && pMuNeg < jpsi::pMuMin[1]) ||
       (fabs(etaMuNeg) > jpsi::etaPS[1] && fabs(etaMuNeg) < jpsi::etaPS[2] && pTMuNeg < jpsi::pTMuMin[2]))
      continue;

//     Reco_StatEv->Fill(4.5);

    Double_t onia_mass = onia->M();
    Double_t onia_pt = onia->Pt();
    Double_t onia_P = onia->P();
    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();
    Double_t onia_phi = onia->Phi();
    Double_t onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);

//     Int_t rapIndex = -1;
//     for(int iRap = 0; iRap < 2*jpsi::kNbRapBins; iRap++){
//       if(onia_rap > jpsi::rapRange[iRap] && onia_rap < jpsi::rapRange[iRap+1]){
//     	rapIndex = iRap;
//     	break;
//       }
//     }
//     Int_t rapForPTIndex = -1;
//     for(int iRap = 0; iRap < jpsi::kNbRapForPTBins; iRap++){
//       if(TMath::Abs(onia_rap) > jpsi::rapForPTRange[iRap] && 
// 	 TMath::Abs(onia_rap) < jpsi::rapForPTRange[iRap+1]){
// 	 rapForPTIndex = iRap+1;
// 	break;
//       }
//     }
//     Int_t pTIndex = -1;
//     for(int iPT = 0; iPT < jpsi::kNbPTBins[rapForPTIndex]; iPT++){
//       if(onia_pt > jpsi::pTRange[rapForPTIndex][iPT] && onia_pt < jpsi::pTRange[rapForPTIndex][iPT+1]){
// 	pTIndex = iPT+1;
// 	break;
//       }
//     }
//     Int_t rapIntegratedPTIndex = -1;
//     for(int iPT = 0; iPT < jpsi::kNbPTBins[0]; iPT++){
//       if(onia_pt > jpsi::pTRange[0][iPT] && onia_pt < jpsi::pTRange[0][iPT+1]){
// 	rapIntegratedPTIndex = iPT+1;
// 	break;
//       }
//     }
//     if(rapIndex < 0){
//       // printf("rapIndex %d, rap(onia) = %f\n", rapIndex, onia_rap);
//       continue;
//     }
//     if(rapForPTIndex < 1){
//       // printf("rapForPTIndex %d, rap(onia) = %f\n", rapForPTIndex, onia_rap);
//       continue;
//     }
//     if(pTIndex < 1){
//       // printf("pTIndex %d, pT(onia) = %f\n", pTIndex, onia_pt);
//       continue;
//     }

//     Reco_StatEv->Fill(5.5);

//     //select events with a cut on the lifetime to reject NP J/psis:
//     if(Jpsict > jpsi::JpsiCtauMax)
//       continue;

//     Reco_StatEv->Fill(6.5);

//     //select events within a narrow mass window around the J/psi
//     //(rapidity dependence of the resolution --> different mass windows)
//     //values obtained from Gaussian fits in "plotMass.C"
//     Double_t jPsiMassMin = jpsi::polMassJpsi[rapForPTIndex] - jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
//     Double_t jPsiMassMax = jpsi::polMassJpsi[rapForPTIndex] + jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
//     if(onia_mass < jPsiMassMin || onia_mass > jPsiMassMax)
//       continue;

//     Reco_StatEv->Fill(7.5);

    if(TMath::Abs(onia_rap) > jpsi::rapYPS)
      continue;

//     Reco_StatEv->Fill(8.5);

    // //test: reject all Cowboy dimuons
    // if((muPos->Phi() - muNeg->Phi()) > 0)
    //   continue;

    // Reco_StatEv->Fill(9.5);

    //remaining of the events will be used for the analysis
    countRecEvent++;

    //fill mass, phi, pt, eta and rap distributions
    //a) all bins
    Reco_Onia_mass[0][0]->Fill(onia_mass);
    Reco_Onia_mass[rapIntegratedPTIndex_Gen][0]->Fill(onia_mass);
    Reco_Onia_mass[0][rapForPTIndex_Gen]->Fill(onia_mass);
    //
    if(onia_mass > (jpsi::polMassJpsi[rapForPTIndex_Gen] - 2.*jpsi::sigmaMassJpsi[rapForPTIndex_Gen]) &&
       onia_mass < (jpsi::polMassJpsi[rapForPTIndex_Gen] + 2.*jpsi::sigmaMassJpsi[rapForPTIndex_Gen])){
	 Reco_Onia_phi[0][0]->Fill(onia_phi);
	 Reco_Onia_phi[rapIntegratedPTIndex_Gen][0]->Fill(onia_phi);
	 Reco_Onia_phi[0][rapForPTIndex_Gen]->Fill(onia_phi);
	 //      
	 Reco_Onia_pt[0]->Fill(onia_pt);
    }
    //b) individual pT and rap bins:
    Reco_Onia_mass[pTIndex_Gen][rapForPTIndex_Gen]->Fill(onia_mass);
    if(onia_mass > (jpsi::polMassJpsi[rapForPTIndex_Gen] - 2.*jpsi::sigmaMassJpsi[rapForPTIndex_Gen]) &&
       onia_mass < (jpsi::polMassJpsi[rapForPTIndex_Gen] + 2.*jpsi::sigmaMassJpsi[rapForPTIndex_Gen])){

      Reco_Onia_phi[pTIndex_Gen][rapForPTIndex_Gen]->Fill(onia_phi);
      Reco_Onia_pt[rapForPTIndex_Gen]->Fill(onia_pt);

      Reco_Onia_rap_pT->Fill(onia_rap, onia_pt);
    }

    //=====================
//     calcPol(*muPos, *muNeg);
//      calcPol(*muPos_Gen, *muNeg_Gen);
    //=====================

//     //===================================================
//     //calculate delta, the angle between the CS and HX frame
//     //Formula from: EJP C69 (2010) 657
//     Double_t deltaHXToCS = TMath::ACos(onia_mass * onia->Pz() / (onia_mT * onia_P));
//     //     Double_t deltaCSToHX = -deltaHXToCS;
//     Double_t sin2Delta = pow((onia_pt * onia->Energy() / (onia_P * onia_mT)),2);
//     //sin2Delta does not change sign when going from HX-->CS or vice versa
//     hDelta[pTIndex_Gen][rapForPTIndex_Gen]->Fill(deltaHXToCS * 180./TMath::Pi());
//     hSin2Delta[pTIndex_Gen][rapForPTIndex_Gen]->Fill(sin2Delta);
//     //===================================================

//     Double_t deltaPhi = muPos->Phi() - muNeg->Phi();
//     if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
//     else if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();

//     //debugging histos
//     hPhiPos_PhiNeg[pTIndex_Gen][rapForPTIndex_Gen]->Fill(180./TMath::Pi() * muNeg->Phi(), 180./TMath::Pi() * muPos->Phi());
//     hPtPos_PtNeg[pTIndex_Gen][rapForPTIndex_Gen]->Fill(muNeg->Pt(), muPos->Pt());
//     hEtaPos_EtaNeg[pTIndex_Gen][rapForPTIndex_Gen]->Fill(muNeg->PseudoRapidity(), muPos->PseudoRapidity());
//     hDeltaPhi[pTIndex_Gen][rapForPTIndex_Gen]->Fill(deltaPhi);

//     Reco_mupl_pt[pTIndex_Gen][rapForPTIndex_Gen]->Fill(muPos->Pt());
//     Reco_mupl_eta[pTIndex_Gen][rapForPTIndex_Gen]->Fill(muPos->PseudoRapidity());
//     Reco_mupl_phi[pTIndex_Gen][rapForPTIndex_Gen]->Fill(muPos->Phi());
      
//     Reco_mumi_pt[pTIndex_Gen][rapForPTIndex_Gen]->Fill(muNeg->Pt());
//     Reco_mumi_eta[pTIndex_Gen][rapForPTIndex_Gen]->Fill(muNeg->PseudoRapidity());
//     Reco_mumi_phi[pTIndex_Gen][rapForPTIndex_Gen]->Fill(muNeg->Phi());

    //fill the histos for all the different frames
    for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){

// 	Double_t weight = CalcPolWeight(thisCosTh[iFrame]);
      Double_t weight = 1.;

	//histos for neg. and pos. rapidity separately:
	if(rapIndex_Gen >= 0)
	  Reco2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(pTIndex_Gen > 0 && rapIndex_Gen >= 0)
	  Reco2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	
	//histos taking together +y and -y
	Reco2D_pol_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapIntegratedPTIndex_Gen > 0)
	  Reco2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapForPTIndex_Gen > 0)
	  Reco2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0)
	  Reco2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	  
	//histos taking together +y and -y and phi-4-folding
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
	Reco2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(thetaAdjusted, phiFolded, weight);
	if(rapIntegratedPTIndex_Gen > 0)
	  Reco2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thetaAdjusted, phiFolded, weight);
	if(rapForPTIndex_Gen > 0)
	  Reco2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded, weight);
	if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0)
	  Reco2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded, weight);
      }

  }//loop over entries

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}

//==========================================================
Double_t CalcPolWeight(Double_t thisCosTh){

  Double_t lambdaTh = 0.;
  Double_t weight = 1. + lambdaTh * TMath::Power(thisCosTh, 2.);
  return weight;
}
