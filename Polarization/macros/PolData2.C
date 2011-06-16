#define PolData2_cxx
#include "PolData2.h"
#include "calcPol.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>

//some statistics
TH1F *Reco_StatEv;

//histos at reco level for mu+, mu-, gamma, Onia
TH1F *Reco_mupl_pt[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_mumi_pt[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_mupl_eta[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_mumi_eta[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_mupl_phi[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_mumi_phi[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2F *Reco_mupl_eta_pT, *Reco_mumi_eta_pT;
TH2F *Reco_mupl_eta_p, *Reco_mumi_eta_p;
TH1F *Reco_Onia_mass[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_Onia_pt[jpsi::kNbRapForPTBins+1];
TH1F *Reco_Onia_phi[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1F *Reco_Onia_rap[jpsi::kNbPTMaxBins+1];
TH2F *Reco_Onia_rap_pT;
//check the lower pT and p of the single muons vs eta
TH2F *Reco_muHLT_pT_eta[4], *Reco_muHLT_p_eta[4]; //higher pT muon: all, after acc. cuts, triggered; using the "correct" muon as the HLT muon
TH2F *Reco_muTM_pT_eta[4], *Reco_muTM_p_eta[4];//lower pT muon: all, after acc. cuts, triggered; using the "correct" muon as the HLT muon
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
enum {LOOSE,TIGHT};//set of muon fiducial cuts
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

    //if we process MC, we must ensure that we only consider
    //reconstructed events
    if(onia->Pt() > 990.)
      continue;

    if(JpsiVprob < 0.01)
      continue;

    Reco_StatEv->Fill(0.5);//count all events

    //check the trigger flag: 0... no trigger, 1 ... triggered+matched, 3 ... triggered (HLT_DoubleMu0)
    //for a full list of accessible triggers, check out "PolData2.h"
    //and https://espace.cern.ch/cms-quarkonia/onia-polarization/L1%20%20HLT/unprescaledTriggersVsRun.aspx
    //for the run ranges per HLT path during unprescaled running periods

    // ---> For the asymmetric triggers:
    // 0 : event not firing the corresponding trigger
    // 1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, POSITIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
    // -1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, NEGATIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
    // 2 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, both matched to the tighter HLT object (usually a L3 muon)
    // 3 : event firing the corresponding trigger but at least one RECO muon is not matched to the HLT objects 
    Int_t trigDecision = -99;
    if(runNb >= 140116 && runNb <= 144114)
      trigDecision = HLT_Mu0_TkMu0_Jpsi;
    else if(runNb >= 146428 && runNb <= 148058)
      trigDecision = HLT_Mu0_TkMu0_OST_Jpsi;
    else if(runNb >= 148819 && runNb <= 149182)
      trigDecision = HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2;
    else if(runNb >= 149291 && runNb <= 149442)
      trigDecision = HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3;

    //for usage of the HLT_DoubleMu0* paths:
//     if(runNb >= 133446 && runNb <= 147116)
//       trigDecision = HLT_DoubleMu0;
//     else if(runNb >= 147196 && runNb <= 149442)
//       trigDecision = HLT_DoubleMu0_Quarkonium_v1;
    
    // //alternatively, for MC use:
    // trigDecision = HLT_DoubleMu0;
//     trigDecision = HLT_Mu0_TkMu0_OST_Jpsi;

    //reject processing of events where the dimuon type (GG, GT or TT)
    //does not correspond to the chosen one
    if(selDimuType < 3 && JpsiType != selDimuType)
      continue;
    else if(selDimuType == 3 && JpsiType > 1) //only GG or GT
      continue;

    Reco_StatEv->Fill(1.5);//count all events

    Double_t onia_mass = onia->M();
    Double_t onia_pt = onia->Pt();
    Double_t onia_P = onia->P();
    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();
    Double_t onia_phi = onia->Phi();
    Double_t onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);

    if(TMath::Abs(onia_rap) > jpsi::rapYPS)
      continue;
    Reco_StatEv->Fill(2.5);

    Double_t etaMuPos = muPos->PseudoRapidity();
    Double_t etaMuNeg = muNeg->PseudoRapidity();
    Double_t pTMuPos = muPos->Pt();
    Double_t pTMuNeg = muNeg->Pt();
    Double_t pMuPos = muPos->P();
    Double_t pMuNeg = muNeg->P();

    //select events with a cut on the lifetime to reject NP J/psis:
    if(Jpsict > jpsi::JpsiCtauMax)
      continue;
    Reco_StatEv->Fill(3.5);

    Double_t jPsiMassMin = jpsi::polMassJpsi[0] - jpsi::nSigMass*jpsi::sigmaMassJpsi[0];
    Double_t jPsiMassMax = jpsi::polMassJpsi[0] + jpsi::nSigMass*jpsi::sigmaMassJpsi[0];
    if(onia_mass < jPsiMassMin || onia_mass > jPsiMassMax)
      continue;

    Reco_StatEv->Fill(4.5);

    Reco_mupl_eta_pT->Fill(etaMuPos, pTMuPos);
    Reco_mumi_eta_pT->Fill(etaMuNeg, pTMuNeg);
    Reco_mupl_eta_p->Fill(etaMuPos, pMuPos);
    Reco_mumi_eta_p->Fill(etaMuNeg, pMuNeg);

    if(pTMuPos > pTMuNeg){
      Reco_muHLT_pT_eta[0]->Fill(fabs(etaMuPos), pTMuPos);
      Reco_muHLT_p_eta[0]->Fill(fabs(etaMuPos), pMuPos);
      Reco_muTM_pT_eta[0]->Fill(fabs(etaMuNeg), pTMuNeg);
      Reco_muTM_p_eta[0]->Fill(fabs(etaMuNeg), pMuNeg);
    }
    else{
      Reco_muHLT_pT_eta[0]->Fill(fabs(etaMuNeg), pTMuNeg);
      Reco_muHLT_p_eta[0]->Fill(fabs(etaMuNeg), pMuNeg);
      Reco_muTM_pT_eta[0]->Fill(fabs(etaMuPos), pTMuPos);
      Reco_muTM_p_eta[0]->Fill(fabs(etaMuPos), pMuPos);
    }

    // //take muons only within a certain eta range
    // if((fabs(etaMuPos) < jpsi::etaPS[0] && pTMuPos < jpsi::pTMuMin[0]) || //mid-rapidity cut
    //    (fabs(etaMuPos) > jpsi::etaPS[0] && fabs(etaMuPos) < jpsi::etaPS[1] && pMuPos < jpsi::pMuMin[1]) ||
    //    (fabs(etaMuPos) > jpsi::etaPS[1] && fabs(etaMuPos) < jpsi::etaPS[2] && pTMuPos < jpsi::pTMuMin[2]))
    //   continue;
    // //(b) on the negative muon
    // if((fabs(etaMuNeg) < jpsi::etaPS[0] && pTMuNeg < jpsi::pTMuMin[0]) || //mid-rapidity cut
    //    (fabs(etaMuNeg) > jpsi::etaPS[0] && fabs(etaMuNeg) < jpsi::etaPS[1] && pMuNeg < jpsi::pMuMin[1]) ||
    //    (fabs(etaMuNeg) > jpsi::etaPS[1] && fabs(etaMuNeg) < jpsi::etaPS[2] && pTMuNeg < jpsi::pTMuMin[2]))
    //   continue;

    //apply different fiducial cuts for the muon matched to the HLT leg or the TM leg
    //for simplicity we decided to apply the stronger cuts to the higher pT muon and not
    //strictly to the HLT muon
    Bool_t muonsInAcc = kFALSE;
    if(pTMuPos > pTMuNeg){
      muonsInAcc = areMuonsInAcceptance(LOOSE, pTMuPos, etaMuPos, pTMuNeg, etaMuNeg); //code needs to be adjusted for the DoubleMu0 trigger!!!!
      //muonsInAcc = areMuonsInAcceptance(TIGHT, pTMuPos, etaMuPos, pTMuNeg, etaMuNeg); //code needs to be adjusted for the DoubleMu0 trigger!!!!
      if(muonsInAcc){
	Reco_muHLT_pT_eta[1]->Fill(fabs(etaMuPos), pTMuPos);
	Reco_muHLT_p_eta[1]->Fill(fabs(etaMuPos), pMuPos);
	Reco_muTM_pT_eta[1]->Fill(fabs(etaMuNeg), pTMuNeg);
	Reco_muTM_p_eta[1]->Fill(fabs(etaMuNeg), pMuNeg);
      }
    }
    else{
      muonsInAcc = areMuonsInAcceptance(LOOSE, pTMuNeg, etaMuNeg, pTMuPos, etaMuPos); //code needs to be adjusted for the DoubleMu0 trigger!!!!
      //muonsInAcc = areMuonsInAcceptance(TIGHT, pTMuNeg, etaMuNeg, pTMuPos, etaMuPos); //code needs to be adjusted for the DoubleMu0 trigger!!!!
      if(muonsInAcc){
	Reco_muHLT_pT_eta[1]->Fill(fabs(etaMuNeg), pTMuNeg);
	Reco_muHLT_p_eta[1]->Fill(fabs(etaMuNeg), pMuNeg);
	Reco_muTM_pT_eta[1]->Fill(fabs(etaMuPos), pTMuPos);
	Reco_muTM_p_eta[1]->Fill(fabs(etaMuPos), pMuPos);
      }
    }
    if(!muonsInAcc)
      continue;

   Reco_StatEv->Fill(5.5);

    //check the trigger
    if(trigDecision != 1 &&  trigDecision != -1 && trigDecision != 2){
      //       printf("rejecting events in run %d\n", runNb);
      continue;
    }
    if(pTMuPos > pTMuNeg){
      Reco_muHLT_pT_eta[2]->Fill(fabs(etaMuPos), pTMuPos);
      Reco_muHLT_p_eta[2]->Fill(fabs(etaMuPos), pMuPos);
      Reco_muTM_pT_eta[2]->Fill(fabs(etaMuNeg), pTMuNeg);
      Reco_muTM_p_eta[2]->Fill(fabs(etaMuNeg), pMuNeg);
    }
    else{
      Reco_muHLT_pT_eta[2]->Fill(fabs(etaMuNeg), pTMuNeg);
      Reco_muHLT_p_eta[2]->Fill(fabs(etaMuNeg), pMuNeg);
      Reco_muTM_pT_eta[2]->Fill(fabs(etaMuPos), pTMuPos);
      Reco_muTM_p_eta[2]->Fill(fabs(etaMuPos), pMuPos);
    }

    Reco_StatEv->Fill(6.5);

    //definitions done for the low pT J/psi trigger...
    if(trigDecision == 1){
      Reco_muHLT_pT_eta[3]->Fill(fabs(etaMuPos), pTMuPos);
      Reco_muHLT_p_eta[3]->Fill(fabs(etaMuPos), pMuPos);
      Reco_muTM_pT_eta[3]->Fill(fabs(etaMuNeg), pTMuNeg);
      Reco_muTM_p_eta[3]->Fill(fabs(etaMuNeg), pMuNeg);
    }
    else if(trigDecision == -1){
      Reco_muHLT_pT_eta[3]->Fill(fabs(etaMuNeg), pTMuNeg);
      Reco_muHLT_p_eta[3]->Fill(fabs(etaMuNeg), pMuNeg);
      Reco_muTM_pT_eta[3]->Fill(fabs(etaMuPos), pTMuPos);
      Reco_muTM_p_eta[3]->Fill(fabs(etaMuPos), pMuPos);
    }
    else if(trigDecision == 1){
      Reco_muTM_pT_eta[3]->Fill(fabs(etaMuNeg), pTMuNeg);
      Reco_muTM_p_eta[3]->Fill(fabs(etaMuNeg), pMuNeg);
      Reco_muTM_pT_eta[3]->Fill(fabs(etaMuPos), pTMuPos);
      Reco_muTM_p_eta[3]->Fill(fabs(etaMuPos), pMuPos);

      Reco_muHLT_pT_eta[3]->Fill(fabs(etaMuNeg), pTMuNeg);
      Reco_muHLT_p_eta[3]->Fill(fabs(etaMuNeg), pMuNeg);
      Reco_muHLT_pT_eta[3]->Fill(fabs(etaMuPos), pTMuPos);
      Reco_muHLT_p_eta[3]->Fill(fabs(etaMuPos), pMuPos);
    }

    //set up the pT and y indices
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
    Int_t rapIndex = -1;
    for(int iRap = 0; iRap < 2*jpsi::kNbRapBins; iRap++){
      if(onia_rap > jpsi::rapRange[iRap] && onia_rap < jpsi::rapRange[iRap+1]){
    	rapIndex = iRap;
    	break;
      }
    }

    //select events within a narrow mass window around the J/psi
    //(rapidity dependence of the resolution --> different mass windows)
    //values obtained from Gaussian fits in "plotMass.C"
    jPsiMassMin = jpsi::polMassJpsi[rapForPTIndex] - jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
    jPsiMassMax = jpsi::polMassJpsi[rapForPTIndex] + jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
    if(onia_mass < jPsiMassMin || onia_mass > jPsiMassMax)
      continue;

    Reco_StatEv->Fill(7.5);

    //fill mass, phi, pt, eta and rap distributions
    //for all events, before rejecting events with pT > 30 GeV/c etc.
    //a) all bins
    Reco_Onia_rap_pT->Fill(onia_rap, onia_pt);
    Reco_Onia_mass[0][0]->Fill(onia_mass);
    
    Reco_Onia_rap[0]->Fill(onia_rap);
    Reco_Onia_pt[0]->Fill(onia_pt);
    Reco_Onia_phi[0][0]->Fill(onia_phi);

    Reco_mupl_pt  [0][0]->Fill(muPos->Pt());
    Reco_mupl_eta [0][0]->Fill(muPos->PseudoRapidity());
    Reco_mupl_phi [0][0]->Fill(muPos->Phi());
    Reco_mumi_pt [0][0]->Fill(muNeg->Pt());
    Reco_mumi_eta[0][0]->Fill(muNeg->PseudoRapidity());
    Reco_mumi_phi[0][0]->Fill(muNeg->Phi());

    if(rapIntegratedPTIndex >= 0){
      Reco_Onia_mass[rapIntegratedPTIndex][0]->Fill(onia_mass);
      Reco_Onia_phi[rapIntegratedPTIndex][0]->Fill(onia_phi);
      Reco_Onia_rap[rapIntegratedPTIndex]->Fill(onia_rap);

      Reco_mupl_pt [rapIntegratedPTIndex][0]->Fill(muPos->Pt());
      Reco_mupl_eta[rapIntegratedPTIndex][0]->Fill(muPos->PseudoRapidity());
      Reco_mupl_phi[rapIntegratedPTIndex][0]->Fill(muPos->Phi());
      Reco_mumi_pt [rapIntegratedPTIndex][0]->Fill(muNeg->Pt());
      Reco_mumi_eta[rapIntegratedPTIndex][0]->Fill(muNeg->PseudoRapidity());
      Reco_mumi_phi[rapIntegratedPTIndex][0]->Fill(muNeg->Phi());
    }
    if(rapForPTIndex > 0){

      Reco_Onia_mass[0][rapForPTIndex]->Fill(onia_mass);
      Reco_Onia_pt[rapForPTIndex]->Fill(onia_pt);
      Reco_Onia_phi[0][rapForPTIndex]->Fill(onia_phi);
      //
      Reco_mupl_pt [0][rapForPTIndex]->Fill(muPos->Pt());
      Reco_mupl_eta[0][rapForPTIndex]->Fill(muPos->PseudoRapidity());
      Reco_mupl_phi[0][rapForPTIndex]->Fill(muPos->Phi());
      Reco_mumi_pt [0][rapForPTIndex]->Fill(muNeg->Pt());
      Reco_mumi_eta[0][rapForPTIndex]->Fill(muNeg->PseudoRapidity());
      Reco_mumi_phi[0][rapForPTIndex]->Fill(muNeg->Phi());
    }

    //continue only if we have events within the bins we are interested in
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

    Reco_StatEv->Fill(8.5);

    //b) individual pT and rap bins:
    Reco_Onia_mass[pTIndex][rapForPTIndex]->Fill(onia_mass);
    Reco_Onia_phi[pTIndex][rapForPTIndex]->Fill(onia_phi);

    // //test: reject all Cowboy dimuons
    // if((muPos->Phi() - muNeg->Phi()) > 0)
    //   continue;

    // Reco_StatEv->Fill(8.5);

    //remaining of the events will be used for the analysis
    countRecEvent++;
    if(writeOutEvents)
      fprintf(fOutputTextFile, "%d\t%d\t%d\n", (Int_t) runNb, (Int_t) lumiBlock, (Int_t)  eventNb);

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
      Reco_Onia_pol_pT_rap[iFrame][0][0][jpsi::cosThPol]->Fill(thisCosTh[iFrame], weight);
      Reco_Onia_pol_pT_rap[iFrame][0][0][jpsi::phiPol]->Fill(thisPhi[iFrame], weight);
      Reco_Onia_pol_pT_rap[iFrame][0][0][jpsi::cos2PhiPol]->Fill(thisCosPhi[iFrame], weight);

      if(rapIntegratedPTIndex > 0){
	Reco2D_Onia_pol_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	Reco_Onia_pol_pT_rap[iFrame][rapIntegratedPTIndex][0][jpsi::cosThPol]->Fill(thisCosTh[iFrame], weight);
	Reco_Onia_pol_pT_rap[iFrame][rapIntegratedPTIndex][0][jpsi::phiPol]->Fill(thisPhi[iFrame], weight);
	Reco_Onia_pol_pT_rap[iFrame][rapIntegratedPTIndex][0][jpsi::cos2PhiPol]->Fill(thisCosPhi[iFrame], weight);
      }
      if(rapForPTIndex > 0){
	Reco2D_Onia_pol_pT_rap[iFrame][0][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	Reco_Onia_pol_pT_rap[iFrame][0][rapForPTIndex][jpsi::cosThPol]->Fill(thisCosTh[iFrame], weight);
	Reco_Onia_pol_pT_rap[iFrame][0][rapForPTIndex][jpsi::phiPol]->Fill(thisPhi[iFrame], weight);
	Reco_Onia_pol_pT_rap[iFrame][0][rapForPTIndex][jpsi::cos2PhiPol]->Fill(thisCosPhi[iFrame], weight);
      }
      if(pTIndex > 0 && rapForPTIndex > 0){
	Reco2D_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	Reco_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex][jpsi::cosThPol]->Fill(thisCosTh[iFrame], weight);
	Reco_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex][jpsi::phiPol]->Fill(thisPhi[iFrame], weight);
	Reco_Onia_pol_pT_rap[iFrame][pTIndex][rapForPTIndex][jpsi::cos2PhiPol]->Fill(thisCosPhi[iFrame], weight);
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
Double_t PolData2::CalcPolWeight(Double_t onia_P, Double_t thisCosTh){

  Double_t lambdaTh = 0.;
  Double_t weight = 1. + lambdaTh * TMath::Power(thisCosTh, 2.);
  return weight;
}

//==========================================================
Bool_t PolData2::areMuonsInAcceptance(Int_t iCut, Double_t pTHigh, Double_t etaHigh, Double_t pTLow, Double_t etaLow){

  Double_t etaBorderHLT[2][4] = {{0., 1.1, 1.4, 2.4}, {0., 1.2, 1.3, 2.2}}; //LOOSE, TIGHT cuts
  Double_t pTBorderHLT[2][4] = {{4.6, 4.0, 2.75, 2.0}, {5.2, 4.7, 3.3, 3.0}};
  Double_t etaBorderTM[2][4] = {{0., 1.1, 1.65, 2.4}, {0., 1.2, 1.6, 2.25}};
  Double_t pTBorderTM[2][4] = {{3.5, 3.5, 1.2, 1.2}, {3.8, 3.8, 1.75, 1.75}};
  
  Double_t minPT_HLT, minPT_TM;
  Bool_t decision = kFALSE;

  //loop over higher pT muon
  for(int iEta = 0; iEta < 3; iEta++){
    if(fabs(etaHigh) > etaBorderHLT[iCut][iEta] && fabs(etaHigh) < etaBorderHLT[iCut][iEta+1]){
      minPT_HLT = (pTBorderHLT[iCut][iEta+1]-pTBorderHLT[iCut][iEta]) / (etaBorderHLT[iCut][iEta+1]-etaBorderHLT[iCut][iEta]) * (fabs(etaHigh) - etaBorderHLT[iCut][iEta]) + pTBorderHLT[iCut][iEta];
      break;
    }
    else if(fabs(etaHigh) > etaBorderHLT[iCut][3])
      minPT_HLT = 1000.; //reject all events with |eta| > 2.4 (or 2.2, ...)
  }
  //loop over lower pT muon
  for(int iEta = 0; iEta < 3; iEta++){
    if(fabs(etaLow) > etaBorderTM[iCut][iEta] && fabs(etaLow) < etaBorderTM[iCut][iEta+1]){
      minPT_TM = (pTBorderTM[iCut][iEta+1]-pTBorderTM[iCut][iEta]) / (etaBorderTM[iCut][iEta+1]-etaBorderTM[iCut][iEta]) * (fabs(etaLow) - etaBorderTM[iCut][iEta]) + pTBorderTM[iCut][iEta];
      break;
    }
    else if(fabs(etaLow) > etaBorderTM[iCut][3])
      minPT_TM = 1000.; //reject all events with |eta| > 2.4 (or 2.25, ...)
  }
  if(pTHigh > minPT_HLT && pTLow > minPT_TM)
    decision = kTRUE;

  return decision;
}
