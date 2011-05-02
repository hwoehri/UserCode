#define MCTnPEff_cxx
#include "MCTnPEff.h"
#include "calcPol.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>

Int_t const kNbEff = 7;
Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/LucaPerrozzi/28March2011/MuonTrackingEff_28March2011.root",
				"/Users/hwoehri/CMS/Work/TnP/Zongchang/29March2011/MuonIDEff_29March2011_fitted.root",
				"/Users/hwoehri/CMS/Work/TnP/Xianyou/25March2011/MuonQualEff_25March2011_fitted.root",
				"/Users/hwoehri/CMS/Work/TnP/Francesco/2April2011/L1L2_DMu0_TriggerEfficiencies_RUNA_distM2gt120_7April2011_fitted.root",
				"/Users/hwoehri/CMS/Work/TnP/Luigi/5April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM1gt150_7April2011_fitted.root",
				"/Users/hwoehri/CMS/Work/TnP/Herbert/28March2011/HLT_Track_Mu0TkMu0_Run1.root",
				"/Users/hwoehri/CMS/Work/TnP/Ilse/24March2011/HLTMuonTrack_Mu0_TkMu0_TM.root"};
//old versions:
//"/Users/hwoehri/CMS/Work/TnP/Zongchang/19March2011/MuonIDEff_19March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Francesco/23March2011/L1L2_DMu0_TriggerEfficiencies_23March2011.root", //mix of RunA and B
//"/Users/hwoehri/CMS/Work/TnP/Luigi/25March2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_25March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/LucaPerrozzi/17March2011/MuonTrackingEff_17March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Xianyou/19March2011/MuonQualEff_19March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Francesco/18March2011/L1L2_DoubleMu0_TriggerEfficiencies_18March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Luigi/19March2011/L3_DoubleMu0_TriggerEfficiencies_19March2011.root"};
enum {TrkEff, MuIDEff, MuQualEff, L1L2Eff, L3Eff, Trk_TkMu0, Mu_TkMu0};
Char_t *effName[kNbEff] = {"TrkEff", "MuIDEff", "MuQualEff", "L1L2Eff", "L3Eff", "Trk_TkMu0Eff", "Mu_TkMu0Eff"};
Int_t const kNbEffSample = 3;
enum {DATA, MC, MCTRUTH};
Char_t *effSampleName[kNbEffSample] = {"DATA", "MC", "MCTRUTH"};
//
enum {CENTRAL, UPPER, LOWER};
TH2D *hMuEff[kNbEff][kNbEffSample][3];
Int_t const kNbEtaBins = 8;
Float_t binsEta[kNbEtaBins] = {0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4}; //needs manual adjustment!!!
TF1 *fMuEff_pT[kNbEff][kNbEffSample][kNbEtaBins];
//=================================================
//1D histos versus J/psi pT, phi and y
//=================================================
//histos for neg. and pos. rapidity separately:
TH1D *hGen_pT, *hGen_y, *hGen_phi;
TH1D *recoEff_pT, *recoEff_y, *recoEff_phi;
TH1D *trigEff_pT, *trigEff_y, *trigEff_phi;
TH1D *totEff_pT, *totEff_y, *totEff_phi;
//=================================================
//2D histos versus J/psi pT and y
//=================================================
//histos for neg. and pos. rapidity separately:
TH2D *hGen2D_pT_rapNP, *hGen2D_pT_rap;
TH2D *recoEff2D_pT_rapNP, *recoEff2D_pT_rap;
TH2D *trigEff2D_pT_rapNP, *trigEff2D_pT_rap;
TH2D *totEff2D_pT_rapNP, *totEff2D_pT_rap;
//=================================================
//2D polarization histos, for various J/psi pT and y
//=================================================
// //histos for neg. and pos. rapidity separately:
TH2D *hGen2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH2D *hGen2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH2D *hGen2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TH2D *recoEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH2D *trigEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH2D *totEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
//histos taking together +y and -y:
TH2D *recoEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH2D *trigEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH2D *totEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
//histos taking together +y and -y and 4-folding in phi
TH2D *recoEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH2D *trigEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH2D *totEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
//=================================================
//deltaR, deltaPhi vs deltaEta for various J/psi pT and y
//=================================================
TH1D *hGen_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *recoEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *trigEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *totEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TH2D *hGen2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH2D *recoEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH2D *trigEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH2D *totEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];


Double_t GetEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt);
Double_t GetEfficiency_FromParametrization(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt);
//==============================================
void MCTnPEff::Loop(Int_t effSample, Char_t *trigLabel)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nb = 0;
  Long64_t countRecEvent = 0;
  printf("number of entries = %d\n", (Int_t) nentries);

  //loop over the events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//     for (Long64_t jentry=0; jentry<10000;jentry++) {

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
    if((fabs(etaMuPos_Gen) < eff::etaPS[0] && pTMuPos_Gen < eff::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuPos_Gen) > eff::etaPS[0] && fabs(etaMuPos_Gen) < eff::etaPS[1] && pMuPos_Gen < eff::pMuMin[1]) ||
       (fabs(etaMuPos_Gen) > eff::etaPS[1] && fabs(etaMuPos_Gen) < eff::etaPS[2] && pTMuPos_Gen < eff::pTMuMin[2]))
      continue;
    //(b) on the negative muon
    if((fabs(etaMuNeg_Gen) < eff::etaPS[0] && pTMuNeg_Gen < eff::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuNeg_Gen) > eff::etaPS[0] && fabs(etaMuNeg_Gen) < eff::etaPS[1] && pMuNeg_Gen < eff::pMuMin[1]) ||
       (fabs(etaMuNeg_Gen) > eff::etaPS[1] && fabs(etaMuNeg_Gen) < eff::etaPS[2] && pTMuNeg_Gen < eff::pTMuMin[2]))
      continue;

    if(JpsiVprob < 0.01)
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

    // Double_t epsTrack_Pos  = GetEfficiency(TrkEff, effSample, etaMuPos_Gen, pTMuPos_Gen); 
    // Double_t epsMuonID_Pos = GetEfficiency(MuIDEff, effSample, etaMuPos_Gen, pTMuPos_Gen);
    // Double_t epsQual_Pos   = GetEfficiency(MuQualEff, effSample, etaMuPos_Gen, pTMuPos_Gen);
    // Double_t epsTrack_Neg  = GetEfficiency(TrkEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
    // Double_t epsMuonID_Neg = GetEfficiency(MuIDEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
    // Double_t epsQual_Neg   = GetEfficiency(MuQualEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);

    Double_t epsTrack_Pos  = GetEfficiency(TrkEff, effSample, etaMuPos_Gen, pTMuPos_Gen); 
    Double_t epsTrack_Neg  = GetEfficiency(TrkEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
    //
    Double_t epsMuonID_Pos = GetEfficiency_FromParametrization(MuIDEff, effSample, etaMuPos_Gen, pTMuPos_Gen);
    Double_t epsQual_Pos   = GetEfficiency_FromParametrization(MuQualEff, effSample, etaMuPos_Gen, pTMuPos_Gen);
    Double_t epsMuonID_Neg = GetEfficiency_FromParametrization(MuIDEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
    Double_t epsQual_Neg   = GetEfficiency_FromParametrization(MuQualEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
    //
    Double_t epsL1L2Trig_Pos, epsL3Trig_Pos, epsL1L2Trig_Neg, epsL3Trig_Neg;
    if(strncmp("HLT_DoubleMu0", trigLabel, 13) == 0){
//       epsL1L2Trig_Pos = GetEfficiency(L1L2Eff, effSample, etaMuPos_Gen, pTMuPos_Gen);
//       epsL3Trig_Pos   = GetEfficiency(L3Eff, effSample, etaMuPos_Gen, pTMuPos_Gen);
//       epsL1L2Trig_Neg = GetEfficiency(L1L2Eff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
//       epsL3Trig_Neg   = GetEfficiency(L3Eff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
      epsL1L2Trig_Pos = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuPos_Gen, pTMuPos_Gen);
      epsL3Trig_Pos   = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuPos_Gen, pTMuPos_Gen);
      epsL1L2Trig_Neg = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
      epsL3Trig_Neg   = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
    }
    else if(strncmp("HLT_Mu0_TkMu0_OST_Jpsi", trigLabel, 22) == 0){
      if(pTMuPos_Gen > pTMuNeg_Gen){//approximate requirement; needs fixing
	epsL1L2Trig_Pos  = GetEfficiency(L1L2Eff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	epsL3Trig_Pos    = GetEfficiency(L3Eff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	epsL1L2Trig_Neg  = GetEfficiency(Mu_TkMu0, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	epsL3Trig_Neg    = GetEfficiency(Trk_TkMu0, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
      }
      else{
	epsL1L2Trig_Neg  = GetEfficiency(L1L2Eff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	epsL3Trig_Neg    = GetEfficiency(L3Eff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	epsL1L2Trig_Pos  = GetEfficiency(Mu_TkMu0, effSample, etaMuPos_Gen, pTMuPos_Gen);
	epsL3Trig_Pos    = GetEfficiency(Trk_TkMu0, effSample, etaMuPos_Gen, pTMuPos_Gen);
      }
    }

    Double_t recoEff_Pos = epsTrack_Pos * epsMuonID_Pos * epsQual_Pos;
    Double_t recoEff_Neg = epsTrack_Neg * epsMuonID_Neg * epsQual_Neg;
    Double_t recoEff = recoEff_Pos * recoEff_Neg;
    Double_t trigEff = epsL1L2Trig_Pos * epsL3Trig_Pos * epsL1L2Trig_Neg * epsL3Trig_Neg;
    Double_t totEff = trigEff * recoEff;
    
    hGen_pT->Fill(onia_Gen_pt);
    hGen_y->Fill(onia_Gen_rap);
    hGen_phi->Fill(onia_Gen_phi);
    hGen2D_pT_rapNP->Fill(onia_Gen_rap, onia_Gen_pt);
    hGen2D_pT_rap->Fill(fabs(onia_Gen_rap), onia_Gen_pt);

    recoEff_pT->Fill(onia_Gen_pt, recoEff);
    recoEff_y->Fill(onia_Gen_rap, recoEff);
    recoEff_phi->Fill(onia_Gen_phi, recoEff);
    recoEff2D_pT_rapNP->Fill(onia_Gen_rap, onia_Gen_pt, recoEff);
    recoEff2D_pT_rap->Fill(fabs(onia_Gen_rap), onia_Gen_pt, recoEff);

    trigEff_pT->Fill(onia_Gen_pt, trigEff);
    trigEff_y->Fill(onia_Gen_rap, trigEff);
    trigEff_phi->Fill(onia_Gen_phi, trigEff);
    trigEff2D_pT_rapNP->Fill(onia_Gen_rap, onia_Gen_pt, trigEff);
    trigEff2D_pT_rap->Fill(fabs(onia_Gen_rap), onia_Gen_pt, trigEff);

    totEff_pT->Fill(onia_Gen_pt, totEff);
    totEff_y->Fill(onia_Gen_rap, totEff);
    totEff_phi->Fill(onia_Gen_phi, totEff);
    totEff2D_pT_rapNP->Fill(onia_Gen_rap, onia_Gen_pt, totEff);
    totEff2D_pT_rap->Fill(fabs(onia_Gen_rap), onia_Gen_pt, totEff);

    //fill the eff. histos for all the different frames
    for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){

	//histos for neg. and pos. rapidity separately:
      if(rapIndex_Gen >= 0){
	  hGen2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
	  recoEff2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], recoEff);
	  trigEff2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], trigEff);
	  totEff2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], totEff);
      }
      if(pTIndex_Gen > 0 && rapIndex_Gen >= 0){
	hGen2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
	recoEff2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], recoEff);
	trigEff2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], trigEff);
	totEff2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], totEff);
      }
	
      //histos taking together +y and -y
      hGen2D_pol_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
      recoEff2D_pol_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], recoEff);
      trigEff2D_pol_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], trigEff);
      totEff2D_pol_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], totEff);
      if(rapIntegratedPTIndex_Gen > 0){
	hGen2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
	recoEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], recoEff);
	trigEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], trigEff);
	totEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], totEff);
      }
      if(rapForPTIndex_Gen > 0){
	hGen2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
	recoEff2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], recoEff);
	trigEff2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], trigEff);
	totEff2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], totEff);
      }
      if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0){
	hGen2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
	recoEff2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], recoEff);
	trigEff2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], trigEff);
	totEff2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thisCosTh[iFrame], thisPhi[iFrame], totEff);
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
      hGen2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(thetaAdjusted, phiFolded);
      recoEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(thetaAdjusted, phiFolded, recoEff);
      trigEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(thetaAdjusted, phiFolded, trigEff);
      totEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(thetaAdjusted, phiFolded, totEff);
      if(rapIntegratedPTIndex_Gen > 0){
	hGen2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thetaAdjusted, phiFolded);
	recoEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thetaAdjusted, phiFolded, recoEff);
	trigEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thetaAdjusted, phiFolded, trigEff);
	totEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(thetaAdjusted, phiFolded, totEff);
      }
      if(rapForPTIndex_Gen > 0){
	hGen2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded);
	recoEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded, recoEff);
	trigEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded, trigEff);
	totEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded, totEff);
      }
      if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0){
	hGen2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded);
	recoEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded, recoEff);
	trigEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded, trigEff);
	totEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(thetaAdjusted, phiFolded, totEff);
      }
    }

    //fill the series of correlation histos for lab-frame
    Double_t deltaEta = fabs(etaMuPos_Gen - etaMuNeg_Gen);
    Double_t deltaPhi = phiMuPos_Gen - phiMuNeg_Gen;
    if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
    else if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
    Double_t deltaPhiDeg = deltaPhi * 180./TMath::Pi();
    Double_t deltaR = sqrt(pow(deltaEta,2) + pow(deltaPhi,2));

    hGen_deltaR_pT_rap[0][0]->Fill(deltaR);
    recoEff_deltaR_pT_rap[0][0]->Fill(deltaR, recoEff);
    trigEff_deltaR_pT_rap[0][0]->Fill(deltaR, trigEff);
    totEff_deltaR_pT_rap[0][0]->Fill(deltaR, totEff);

    hGen2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(deltaEta, deltaPhiDeg);
    recoEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(deltaEta, deltaPhiDeg, recoEff);
    trigEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(deltaEta, deltaPhiDeg, trigEff);
    totEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(deltaEta, deltaPhiDeg, totEff);

    if(rapIntegratedPTIndex_Gen > 0){
      hGen_deltaR_pT_rap[rapIntegratedPTIndex_Gen][0]->Fill(deltaR);
      recoEff_deltaR_pT_rap[rapIntegratedPTIndex_Gen][0]->Fill(deltaR, recoEff);
      trigEff_deltaR_pT_rap[rapIntegratedPTIndex_Gen][0]->Fill(deltaR, trigEff);
      totEff_deltaR_pT_rap[rapIntegratedPTIndex_Gen][0]->Fill(deltaR, totEff);

      hGen2D_deltaPhiVsDeltaEta_pT_rap[rapIntegratedPTIndex_Gen][0]->Fill(deltaEta, deltaPhiDeg);
      recoEff2D_deltaPhiVsDeltaEta_pT_rap[rapIntegratedPTIndex_Gen][0]->Fill(deltaEta, deltaPhiDeg, recoEff);
      trigEff2D_deltaPhiVsDeltaEta_pT_rap[rapIntegratedPTIndex_Gen][0]->Fill(deltaEta, deltaPhiDeg, trigEff);
      totEff2D_deltaPhiVsDeltaEta_pT_rap[rapIntegratedPTIndex_Gen][0]->Fill(deltaEta, deltaPhiDeg, totEff);
    }
    if(rapForPTIndex_Gen > 0){
      hGen_deltaR_pT_rap[0][rapForPTIndex_Gen]->Fill(deltaR);
      recoEff_deltaR_pT_rap[0][rapForPTIndex_Gen]->Fill(deltaR, recoEff);
      trigEff_deltaR_pT_rap[0][rapForPTIndex_Gen]->Fill(deltaR, trigEff);
      totEff_deltaR_pT_rap[0][rapForPTIndex_Gen]->Fill(deltaR, totEff);

      hGen2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex_Gen]->Fill(deltaEta, deltaPhiDeg);
      recoEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex_Gen]->Fill(deltaEta, deltaPhiDeg, recoEff);
      trigEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex_Gen]->Fill(deltaEta, deltaPhiDeg, trigEff);
      totEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex_Gen]->Fill(deltaEta, deltaPhiDeg, totEff);
    }
    if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0){
      hGen_deltaR_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(deltaR);
      recoEff_deltaR_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(deltaR, recoEff);
      trigEff_deltaR_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(deltaR, trigEff);
      totEff_deltaR_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(deltaR, totEff);

      hGen2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(deltaEta, deltaPhiDeg);
      recoEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(deltaEta, deltaPhiDeg, recoEff);
      trigEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(deltaEta, deltaPhiDeg, trigEff);
      totEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex_Gen][rapForPTIndex_Gen]->Fill(deltaEta, deltaPhiDeg, totEff);
    }

  }//loop over entries

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}

//==============================================================
Double_t GetEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt){

  if(fabs(eta) > 2.4) return 0.; //no acceptance beyond 2.4
//   if(pt > 20.) return 1.; //efficiencies beyond 20 GeV/c not reliable from T&P, but are always 100%
  if(iEff == TrkEff && fabs(eta) > 2.1 && pt > 10.) return 1.; //Data contains no value

  if(iEff >= kNbEff)
    printf("%d not a valid efficiency!!!\n", iEff);
  if(iEffSample >= kNbEffSample)
    printf("%d not a valide efficiency sample!!!\n", iEffSample);

  Int_t binX = hMuEff[iEff][iEffSample][CENTRAL]->GetXaxis()->FindBin(fabs(eta));
  Int_t binY = hMuEff[iEff][iEffSample][CENTRAL]->GetYaxis()->FindBin(pt);
  Double_t eff = hMuEff[iEff][iEffSample][CENTRAL]->GetBinContent(binX, binY);

//   printf("%s, efficiency for |eta|=%1.3f and pT=%1.2f GeV/c is %1.3f\n",
// 	 effName[iEff], fabs(eta), pt, eff);
  return eff;
}

//==============================================================
Double_t GetEfficiency_FromParametrization(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt){

  if(fabs(eta) > 2.4) return 0.; //no acceptance beyond 2.4

  if(iEff >= kNbEff)
    printf("%d not a valid efficiency!!!\n", iEff);
  if(iEffSample >= kNbEffSample)
    printf("%d not a valide efficiency sample!!!\n", iEffSample);

  Int_t etaBin = -1;
  for(int iEta = 0; iEta < kNbEtaBins; iEta++){
    if(fabs(eta) < binsEta[iEta]){
      etaBin = iEta;
      break;
    }
  }

  Double_t eff = fMuEff_pT[iEff][iEffSample][etaBin]->Eval(pt);
  if(eff > 1.) eff = 1.;
  else if(eff < 0.) eff = 0.;

//   printf("%s: eta %f --> mapped into bin %d... efficiency for pt = %f GeV/c is %f\n", effName[iEff], eta, etaBin, pt, eff);

  return eff;

}
