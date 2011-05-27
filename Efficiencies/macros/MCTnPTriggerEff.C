#define MCTnPTriggerEff_cxx
#include "MCTnPTriggerEff.h"
#include "calcPol.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>
#include "TRandom.h"
#include "TMath.h"
//#include "TF1.h"
#include "TEfficiency.h"

Int_t const kNbEff = 9;
Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/LucaPerrozzi/28March2011/MuonTrackingEff_28March2011.root",
				"/Users/hwoehri/CMS/Work/TnP/Zongchang/29March2011/MuonIDEff_29March2011_fitted.root",
				"/Users/hwoehri/CMS/Work/TnP/Xianyou/25March2011/MuonQualEff_25March2011_fitted.root",
				"/Users/hwoehri/CMS/Work/TnP/Francesco/2April2011/L1L2_DMu0_TriggerEfficiencies_RUNA_distM2gt120_7April2011_fitted.root",
				//"/Users/hwoehri/CMS/Work/TnP/Francesco/2April2011/L1L2_DMu0_TriggerEfficiencies_RUNB_distM2gt120_7April2011_fitted.root",
				"/Users/hwoehri/CMS/Work/TnP/Luigi/17April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM2gt120_17April2011_fitted.root",
				// "/Users/hwoehri/CMS/Work/TnP/Herbert/13May2011/HLT_Track_Mu0TkMu0_Run1_13May2011.root",
				"/Users/hwoehri/CMS/Work/TnP/Herbert/13May2011/HLT_Track_Mu0TkMu0_Run1_13May2011_MCCorrected.root",
				"/Users/hwoehri/CMS/Work/TnP/Ilse/29March2011/HLTMuonTrack_Mu0_TkMu0_TM_runB1_29March2011.root",
				//"/Users/hwoehri/CMS/CMSSW/hWoehri/Efficiencies/macros/TMEffConditionalToHLT_MC_20May2011.root"};
				"/Users/hwoehri/CMS/Work/TnP/Herbert/24May2011/MuX_L2Mu0_L3Mu0.root",
				"/Users/hwoehri/CMS/Work/TnP/Ilse/23May2011/TkMu0L1L2-pt-abseta-runA-distM2gt120_modified.root"};
//old versions:
//"/Users/hwoehri/CMS/Work/TnP/Zongchang/19March2011/MuonIDEff_19March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Francesco/23March2011/L1L2_DMu0_TriggerEfficiencies_23March2011.root", //mix of RunA and B
//"/Users/hwoehri/CMS/Work/TnP/Luigi/25March2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_25March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Luigi/5April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM1gt150_7April2011_fitted.root",
//"/Users/hwoehri/CMS/Work/TnP/LucaPerrozzi/17March2011/MuonTrackingEff_17March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Xianyou/19March2011/MuonQualEff_19March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Francesco/18March2011/L1L2_DoubleMu0_TriggerEfficiencies_18March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Luigi/19March2011/L3_DoubleMu0_TriggerEfficiencies_19March2011.root"};
//"/Users/hwoehri/CMS/Work/TnP/Herbert/28March2011/HLT_Track_Mu0TkMu0_Run1.root",
enum {TrkEff, MuIDEff, MuQualEff, L1L2Eff, L3Eff, Trk_TkMu0, Mu_TkMu0, TMandHLT1Eff, TMandHLT2Eff};
Char_t *effName[kNbEff] = {"TrkEff", "MuIDEff", "MuQualEff", "L1L2Eff", "L3Eff", "Trk_TkMu0Eff", "Mu_TkMu0Eff", "Track_TkMu0andL1L2", "Mu_TkMu0andL1L2"};
Int_t const kNbEffSample = 3;
enum {DATA, MC, MCTRUTH};
Char_t *effSampleName[kNbEffSample] = {"DATA", "MC", "MCTRUTH"};
//
enum {CENTRAL, UPPER, LOWER};
TH2D *hMuEff[kNbEff][kNbEffSample][3];
Int_t const kNbMaxEtaBins = 14;
Int_t const kNbEtaBins[kNbEff] = {0,8,8,8,14,0,0,0,8};//number of TGraphs with parametrized pT diff. efficiencies
Float_t binsEta[kNbEff][kNbMaxEtaBins] = {{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},//needed to fetch the correct eta bin for the fitted pT turn on curves
					  {0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4, 0.,0.,0.,0.,0.,0.},//MuID
					  {0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4, 0.,0.,0.,0.,0.,0.},//MuQual
					  {0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4, 0.,0.,0.,0.,0.,0.},//L1L2
					  {0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4},//L3
					  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},//Trk_TkMu0
					  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},//Mu_TkMu0
					  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},//TMandHLT1Eff
					  {0.2, 0.3, 0.6, 0.8, 1.2, 1.6, 2.1, 2.4, 0.,0.,0.,0.,0.,0.}};//TMandHLT2Eff
TF1 *fMuEff_pT[kNbEff][kNbEffSample][kNbMaxEtaBins];
//=================================================
//1D histos versus J/psi pT, phi and y
//=================================================
//histos for neg. and pos. rapidity separately:
// TH1D *hGen_pT, *hGen_y, *hGen_phi;
// TH1D *recoEff_pT, *recoEff_y, *recoEff_phi;
// TH1D *trigEff_pT, *trigEff_y, *trigEff_phi;
// TH1D *totEff_pT, *totEff_y, *totEff_phi;
TEfficiency *trigEff_pT, *trigEff_y, *trigEff_phi;
//=================================================
//2D histos versus J/psi pT and y
//=================================================
//histos for neg. and pos. rapidity separately:
// TH2D *hGen2D_pT_rapNP, *hGen2D_pT_rap;
// TH2D *recoEff2D_pT_rapNP, *recoEff2D_pT_rap;
// TH2D *trigEff2D_pT_rapNP, *trigEff2D_pT_rap;
// TH2D *totEff2D_pT_rapNP, *totEff2D_pT_rap;
TEfficiency *trigEff2D_pT_rap;
TEfficiency *trigEff2D_pT_rapNP;
//=================================================
//2D polarization histos, for various J/psi pT and y
//=================================================
// //histos for neg. and pos. rapidity separately:
// TH2D *hGen2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH2D *hGen2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH2D *hGen2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

// TH2D *recoEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH2D *trigEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH2D *totEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// //histos taking together +y and -y:
// TH2D *recoEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH2D *trigEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH2D *totEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// //histos taking together +y and -y and 4-folding in phi
// TH2D *recoEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH2D *trigEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH2D *totEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

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

TH1D *hGen_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *recoEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *trigEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *totEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TH1D *hGen_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *recoEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *trigEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *totEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TH1D *hGen_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *recoEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *trigEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *totEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TH1D *hGen_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *recoEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *trigEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TH1D *totEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

Double_t GetEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt);
Double_t GetEfficiency_FromParametrization(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt);
//==============================================
void MCTnPTriggerEff::Loop(Int_t effSample, Char_t *trigLabel)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nb = 0;
  Long64_t countRecEvent = 0;
  printf("number of entries = %d\n", (Int_t) nentries);

  Bool_t incrementTrig;
  Double_t epsL1L2Trig_Pos, epsL3Trig_Pos, epsL1L2Trig_Neg, epsL3Trig_Neg;
  Double_t eps1TMandHLT_Pos, eps1TMandHLT_Neg;
  Double_t eps2TMandHLT_Pos, eps2TMandHLT_Neg;
  Double_t epsTMandHLT_Pos, epsTMandHLT_Neg;
  Double_t trigEff = 0.;

  //loop over the events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  // for (Long64_t jentry=0; jentry<500;jentry++) {

    if(jentry % 100000 == 0) printf("event %d\n", (Int_t) jentry);
    //printf("event %d\n", (Int_t) jentry);

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    //check that this is a RECO event
    if(onia->Pt() > 990.)
      continue;

    Double_t etaMuPos = muPos->PseudoRapidity();
    Double_t etaMuNeg = muNeg->PseudoRapidity();
    Double_t pTMuPos = muPos->Pt();
    Double_t pTMuNeg = muNeg->Pt();
    Double_t pMuPos = muPos->P();
    Double_t pMuNeg = muNeg->P();
    Double_t phiMuPos = muPos->Phi();
    Double_t phiMuNeg = muNeg->Phi();

    //take muons only within a certain eta range
    if((fabs(etaMuPos) < eff::etaPS[0] && pTMuPos < eff::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuPos) > eff::etaPS[0] && fabs(etaMuPos) < eff::etaPS[1] && pMuPos < eff::pMuMin[1]) ||
       (fabs(etaMuPos) > eff::etaPS[1] && fabs(etaMuPos) < eff::etaPS[2] && pTMuPos < eff::pTMuMin[2]))
      continue;
    //(b) on the negative muon
    if((fabs(etaMuNeg) < eff::etaPS[0] && pTMuNeg < eff::pTMuMin[0]) || //mid-rapidity cut
       (fabs(etaMuNeg) > eff::etaPS[0] && fabs(etaMuNeg) < eff::etaPS[1] && pMuNeg < eff::pMuMin[1]) ||
       (fabs(etaMuNeg) > eff::etaPS[1] && fabs(etaMuNeg) < eff::etaPS[2] && pTMuNeg < eff::pTMuMin[2]))
      continue;

    if(fabs(onia->Rapidity()) > eff::rapMax)
      continue;
    if(JpsiVprob > 0.01) 
      continue;

    Double_t onia_mass = onia->M();
    Double_t onia_pt = onia->Pt();
    Double_t onia_P = onia->P();
    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();
    Double_t onia_phi = onia->Phi();
    Double_t onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);
    
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

    //==============================
    calcPol(*muPos, *muNeg);
    //==============================

    Bool_t usePTFit = kTRUE; //alternative to the T&P histograms use fitted pT differential efficiency

    // ---> For the asymmetric triggers:
    // 0 : event not firing the corresponding trigger
    // 1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, POSITIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
    // -1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, NEGATIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
    // 2 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, both matched to the tighter HLT object (usually a L3 muon)
    // 3 : event firing the corresponding trigger but at least one RECO muon is not matched to the HLT objects 

    if(strncmp("HLT_DoubleMu0", trigLabel, 13) == 0){//DoubleMu0 trigger
      if(usePTFit){
	epsL1L2Trig_Pos = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuPos, pTMuPos);
	epsL3Trig_Pos   = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuPos, pTMuPos);
	epsL1L2Trig_Neg = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuNeg, pTMuNeg);
	epsL3Trig_Neg   = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuNeg, pTMuNeg);
      }
      else{
	epsL1L2Trig_Pos = GetEfficiency(L1L2Eff, effSample, etaMuPos, pTMuPos);
	epsL3Trig_Pos   = GetEfficiency(L3Eff, effSample, etaMuPos, pTMuPos);
	epsL1L2Trig_Neg = GetEfficiency(L1L2Eff, effSample, etaMuNeg, pTMuNeg);
	epsL3Trig_Neg   = GetEfficiency(L3Eff, effSample, etaMuNeg, pTMuNeg);
      }
      trigEff = epsL1L2Trig_Pos * epsL3Trig_Pos * epsL1L2Trig_Neg * epsL3Trig_Neg;
    }
    else if(strncmp("HLT_Mu0_TkMu0_OST_Jpsi", trigLabel, 22) == 0){//any of the "low pT J/psi triggers" --> steering via the input filename

      eps1TMandHLT_Pos = GetEfficiency(TMandHLT1Eff, effSample, etaMuPos, pTMuPos);
      eps1TMandHLT_Neg = GetEfficiency(TMandHLT1Eff, effSample, etaMuNeg, pTMuNeg);

      if(usePTFit){
	epsL1L2Trig_Pos  = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuPos, pTMuPos);
	epsL3Trig_Pos    = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuPos, pTMuPos);
	epsL1L2Trig_Neg  = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuNeg, pTMuNeg);
	epsL3Trig_Neg    = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuNeg, pTMuNeg);

	eps2TMandHLT_Pos = GetEfficiency_FromParametrization(TMandHLT2Eff, effSample, etaMuPos, pTMuPos);
	eps2TMandHLT_Neg = GetEfficiency_FromParametrization(TMandHLT2Eff, effSample, etaMuNeg, pTMuNeg);
      }
      else{
	epsL1L2Trig_Pos  = GetEfficiency(L1L2Eff, effSample, etaMuPos, pTMuPos);
	epsL3Trig_Pos    = GetEfficiency(L3Eff, effSample, etaMuPos, pTMuPos);
	epsL1L2Trig_Neg  = GetEfficiency(L1L2Eff, effSample, etaMuNeg, pTMuNeg);
	epsL3Trig_Neg    = GetEfficiency(L3Eff, effSample, etaMuNeg, pTMuNeg);

	eps2TMandHLT_Pos = GetEfficiency(TMandHLT2Eff, effSample, etaMuPos, pTMuPos);
	eps2TMandHLT_Neg = GetEfficiency(TMandHLT2Eff, effSample, etaMuNeg, pTMuNeg);
      }

      epsTMandHLT_Pos = eps1TMandHLT_Pos * eps2TMandHLT_Pos;
      epsTMandHLT_Neg = eps1TMandHLT_Neg * eps2TMandHLT_Neg;


      Double_t epsTM1_Neg  = GetEfficiency(Trk_TkMu0, effSample, etaMuNeg, pTMuNeg);
      Double_t epsTM2_Neg  = GetEfficiency(Mu_TkMu0, effSample, etaMuNeg, pTMuNeg);
      
      Double_t epsHLT_Pos = epsL1L2Trig_Pos * epsL3Trig_Pos;
      Double_t epsHLT_Neg = epsL1L2Trig_Neg * epsL3Trig_Neg;
      Double_t epsTM_Neg = epsTM1_Neg *epsTM2_Neg;
      
      Double_t epsTM1_Pos = GetEfficiency(Trk_TkMu0, effSample, etaMuPos, pTMuPos);
      Double_t epsTM2_Pos = GetEfficiency(Mu_TkMu0, effSample, etaMuPos, pTMuPos);
      Double_t epsTM_Pos = epsTM1_Pos *epsTM2_Pos;

      trigEff = epsHLT_Pos * (epsTM_Neg - epsTMandHLT_Neg) + TMath::Max(0., (epsHLT_Pos - epsTMandHLT_Pos)) * epsTMandHLT_Neg; //last version
      //unconditioned (pos.) HLT muon and exclusive (neg.) TM + exclusive (pos.) HLT and incluseve (neg.) HLT 
      trigEff += epsHLT_Neg * (epsTM_Pos - epsTMandHLT_Pos) + TMath::Max(0., (epsHLT_Neg - epsTMandHLT_Neg)) * epsTMandHLT_Pos; //last version
      //mirror image of the above (exchanging the charge)
      trigEff += epsTMandHLT_Pos * epsTMandHLT_Neg; 
      //combination of 2 inclusve HLT and TM 

      // if((epsTM_Neg - epsTMandHLT_Neg) < 0.)//this component is always okay
      // 	printf("negative value for (epsTM_Neg - epsTMandHLT_Neg) = %1.3f - %1.3f = %f\n", epsTM_Neg, epsTMandHLT_Neg, epsTM_Neg - epsTMandHLT_Neg);
      // if((epsHLT_Pos - epsTMandHLT_Pos) < 0.)
      // 	printf("negative value for (epsHLT_Pos - epsTMandHLT_Pos) = %1.3f - %1.3f = %f\n", epsHLT_Pos, epsTMandHLT_Pos, epsHLT_Pos - epsTMandHLT_Pos);
      // if((epsTM_Pos - epsTMandHLT_Pos) < 0.)//this component is always okay
      // 	printf("negative value for (epsTM_Pos - epsTMandHLT_Pos) = %1.3f - %1.3f = %f\n", epsTM_Pos, epsTMandHLT_Pos, epsTM_Pos - epsTMandHLT_Pos);
      // if((epsHLT_Neg - epsTMandHLT_Neg) < 0.)
      // 	printf("negative value for (epsHLT_Neg - epsTMandHLT_Neg): %1.3f - %1.3f = %f\n", epsHLT_Neg, epsTMandHLT_Neg, epsHLT_Neg - epsTMandHLT_Neg);


    	if(trigEff < 0.)
    	  trigEff = 0.;
    	else if(trigEff > 1.)
    	  trigEff = 1.;
    }//end of low pT J/psi trigger efficiency

    //the indiv. histograms will be filled depending on the
    //assigned probability
    Double_t randNb = gRandom->Uniform();

    if(trigEff > randNb)
      incrementTrig = kTRUE;
    else
      incrementTrig = kFALSE;

    // hGen_pT->Fill(onia_pt);
    // hGen_y->Fill(onia_rap);
    // hGen_phi->Fill(onia_phi);
    // hGen2D_pT_rapNP->Fill(onia_rap, onia_pt);
    // hGen2D_pT_rap->Fill(fabs(onia_rap), onia_pt);

    // if(incrementTrig){
    //   trigEff_pT->Fill(onia_pt);
    //   trigEff_y->Fill(onia_rap);
    //   trigEff_phi->Fill(onia_phi);
    //   trigEff2D_pT_rapNP->Fill(onia_rap, onia_pt);
    //   trigEff2D_pT_rap->Fill(fabs(onia_rap), onia_pt);
    // 	//      }
    //   }
    // }
    // if(incrementTot){
    //   totEff_pT->Fill(onia_pt);
    //   totEff_y->Fill(onia_rap);
    //   totEff_phi->Fill(onia_phi);
    //   totEff2D_pT_rapNP->Fill(onia_rap, onia_pt);
    //   totEff2D_pT_rap->Fill(fabs(onia_rap), onia_pt);
    // }

    trigEff_pT->Fill(incrementTrig, onia_pt);
    trigEff_y->Fill(incrementTrig, onia_rap);
    trigEff_phi->Fill(incrementTrig, onia_phi);
      
    trigEff2D_pT_rap->Fill(incrementTrig, fabs(onia_rap), onia_pt);
    trigEff2D_pT_rapNP->Fill(incrementTrig, onia_rap, onia_pt);

    //fill the eff. histos for all the different frames
    for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){

      //histos for neg. and pos. rapidity separately:
      if(rapIndex >= 0)
	trigEff2D_pol_pT_rapNP[iFrame][0][rapIndex]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
      if(pTIndex > 0 && rapIndex >= 0)
	trigEff2D_pol_pT_rapNP[iFrame][pTIndex][rapIndex]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	
      //histos taking together +y and -y
      trigEff2D_pol_pT_rap[iFrame][0][0]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
      trigEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
      if(rapForPTIndex > 0)
	trigEff2D_pol_pT_rap[iFrame][0][rapForPTIndex]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
      if(pTIndex > 0 && rapForPTIndex > 0)
	trigEff2D_pol_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	  
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

      trigEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(incrementTrig, thetaAdjusted, phiFolded);
      if(rapIntegratedPTIndex > 0)
	trigEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex][0]->Fill(incrementTrig, thetaAdjusted, phiFolded);
      if(rapForPTIndex > 0)
	trigEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex]->Fill(incrementTrig, thetaAdjusted, phiFolded);
      if(pTIndex > 0 && rapForPTIndex > 0)
	trigEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex][rapForPTIndex]->Fill(incrementTrig, thetaAdjusted, phiFolded);
    }

    // //fill the series of correlation histos for lab-frame
    // Double_t deltaEta = fabs(etaMuPos - etaMuNeg);
    // Double_t deltaPhi = phiMuPos - phiMuNeg;
    // if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
    // else if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
    // Double_t deltaPhiDeg = deltaPhi * 180./TMath::Pi();
    // Double_t deltaR = sqrt(pow(deltaEta,2) + pow(deltaPhi,2));

    // //convert the deltaPhi variable from the TTree from rad to deg:
    // Double_t JpsiDetaM2 = sqrt(pow(JpsiDrM2,2) - pow(JpsiDphiM2,2));
    // JpsiDphiM2 *= 180./TMath::Pi();


    // hGen_deltaR_pT_rap[0][0]->Fill(deltaR);
    // hGen_deltaRM2_pT_rap[0][0]->Fill(JpsiDrM2);
    // hGen_deltaPhiM2_pT_rap[0][0]->Fill(JpsiDphiM2);
    // hGen_deltaEtaM2_pT_rap[0][0]->Fill(JpsiDetaM2);
    // hGen_distM2_pT_rap[0][0]->Fill(JpsiDistM2);
    // if(incrementReco){
    //   recoEff_deltaR_pT_rap[0][0]->Fill(deltaR);
    //   recoEff_deltaRM2_pT_rap[0][0]->Fill(JpsiDrM2);
    //   recoEff_deltaPhiM2_pT_rap[0][0]->Fill(JpsiDphiM2);
    //   recoEff_deltaEtaM2_pT_rap[0][0]->Fill(JpsiDetaM2);
    //   recoEff_distM2_pT_rap[0][0]->Fill(JpsiDistM2);
    //   //    }
    //   if(incrementTrig){
    // 	trigEff_deltaR_pT_rap[0][0]->Fill(deltaR);
    // 	trigEff_deltaRM2_pT_rap[0][0]->Fill(JpsiDrM2);
    // 	trigEff_deltaPhiM2_pT_rap[0][0]->Fill(JpsiDphiM2);
    // 	trigEff_deltaEtaM2_pT_rap[0][0]->Fill(JpsiDetaM2);
    // 	trigEff_distM2_pT_rap[0][0]->Fill(JpsiDistM2);
    //   }
    // }
    // if(incrementTot){
    //   totEff_deltaR_pT_rap[0][0]->Fill(deltaR);
    //   totEff_deltaRM2_pT_rap[0][0]->Fill(JpsiDrM2);
    //   totEff_deltaPhiM2_pT_rap[0][0]->Fill(JpsiDphiM2);
    //   totEff_deltaEtaM2_pT_rap[0][0]->Fill(JpsiDetaM2);
    //   totEff_distM2_pT_rap[0][0]->Fill(JpsiDistM2);
    // }
    // hGen2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(deltaEta, deltaPhiDeg);
    // if(incrementReco){
    //   recoEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(deltaEta, deltaPhiDeg);
    //   if(incrementTrig)
    // 	trigEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(deltaEta, deltaPhiDeg);
    // }
    // if(incrementTot)
    //   totEff2D_deltaPhiVsDeltaEta_pT_rap[0][0]->Fill(deltaEta, deltaPhiDeg);

    // if(rapIntegratedPTIndex > 0){
    //   hGen_deltaR_pT_rap[rapIntegratedPTIndex][0]->Fill(deltaR);
    //   hGen_deltaRM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDrM2);
    //   hGen_deltaPhiM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDphiM2);
    //   hGen_deltaEtaM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDetaM2);
    //   hGen_distM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDistM2);
    //   if(incrementReco){
    // 	recoEff_deltaR_pT_rap[rapIntegratedPTIndex][0]->Fill(deltaR);
    // 	recoEff_deltaRM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDrM2);
    // 	recoEff_deltaPhiM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDphiM2);
    // 	recoEff_deltaEtaM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDetaM2);
    // 	recoEff_distM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDistM2);
    //   // }
    // 	if(incrementTrig){
    // 	  trigEff_deltaR_pT_rap[rapIntegratedPTIndex][0]->Fill(deltaR);
    // 	  trigEff_deltaRM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDrM2);
    // 	  trigEff_deltaPhiM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDphiM2);
    // 	  trigEff_deltaEtaM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDetaM2);
    // 	  trigEff_distM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDistM2);
    // 	}
    //   }
    //   if(incrementTot){
    // 	totEff_deltaR_pT_rap[rapIntegratedPTIndex][0]->Fill(deltaR);
    // 	totEff_deltaRM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDrM2);
    // 	totEff_deltaPhiM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDphiM2);
    // 	totEff_deltaEtaM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDetaM2);
    // 	totEff_distM2_pT_rap[rapIntegratedPTIndex][0]->Fill(JpsiDistM2);
    //   }

    //   hGen2D_deltaPhiVsDeltaEta_pT_rap[rapIntegratedPTIndex][0]->Fill(deltaEta, deltaPhiDeg);
    //   if(incrementReco){
    // 	recoEff2D_deltaPhiVsDeltaEta_pT_rap[rapIntegratedPTIndex][0]->Fill(deltaEta, deltaPhiDeg);
    // 	if(incrementTrig)
    // 	  trigEff2D_deltaPhiVsDeltaEta_pT_rap[rapIntegratedPTIndex][0]->Fill(deltaEta, deltaPhiDeg);
    //   }
    //   if(incrementTot)
    // 	totEff2D_deltaPhiVsDeltaEta_pT_rap[rapIntegratedPTIndex][0]->Fill(deltaEta, deltaPhiDeg);
    // }
    // if(rapForPTIndex > 0){
    //   hGen_deltaR_pT_rap[0][rapForPTIndex]->Fill(deltaR);
    //   hGen_deltaRM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDrM2);
    //   hGen_deltaPhiM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDphiM2);
    //   hGen_deltaEtaM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDetaM2);
    //   hGen_distM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDistM2);
    //   if(incrementReco){
    // 	recoEff_deltaR_pT_rap[0][rapForPTIndex]->Fill(deltaR);
    // 	recoEff_deltaRM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDrM2);
    // 	recoEff_deltaPhiM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDphiM2);
    // 	recoEff_deltaEtaM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDetaM2);
    // 	recoEff_distM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDistM2);
    //   // }
    // 	if(incrementTrig){
    // 	  trigEff_deltaR_pT_rap[0][rapForPTIndex]->Fill(deltaR);
    // 	  trigEff_deltaRM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDrM2);
    // 	  trigEff_deltaPhiM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDphiM2);
    // 	  trigEff_deltaEtaM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDetaM2);
    // 	  trigEff_distM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDistM2);
    // 	}
    //   }
    //   if(incrementTot){
    // 	totEff_deltaR_pT_rap[0][rapForPTIndex]->Fill(deltaR);
    // 	totEff_deltaRM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDrM2);
    // 	totEff_deltaPhiM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDphiM2);
    // 	totEff_deltaEtaM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDetaM2);
    // 	totEff_distM2_pT_rap[0][rapForPTIndex]->Fill(JpsiDistM2);
    //   }

    //   hGen2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex]->Fill(deltaEta, deltaPhiDeg);
    //   if(incrementReco){
    // 	recoEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex]->Fill(deltaEta, deltaPhiDeg);
    // 	if(incrementTrig)
    // 	  trigEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex]->Fill(deltaEta, deltaPhiDeg);
    //   }
    //   if(incrementTot)
    // 	totEff2D_deltaPhiVsDeltaEta_pT_rap[0][rapForPTIndex]->Fill(deltaEta, deltaPhiDeg);
    // }
    // if(pTIndex > 0 && rapForPTIndex > 0){
    //   hGen_deltaR_pT_rap[pTIndex][rapForPTIndex]->Fill(deltaR);
    //   hGen_deltaRM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDrM2);
    //   hGen_deltaPhiM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDphiM2);
    //   hGen_deltaEtaM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDetaM2);
    //   hGen_distM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDistM2);
    //   if(incrementReco){
    // 	recoEff_deltaR_pT_rap[pTIndex][rapForPTIndex]->Fill(deltaR);
    // 	recoEff_deltaRM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDrM2);
    // 	recoEff_deltaPhiM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDphiM2);
    // 	recoEff_deltaEtaM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDetaM2);
    // 	recoEff_distM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDistM2);
    // 	//}
    // 	if(incrementTrig){
    // 	  trigEff_deltaR_pT_rap[pTIndex][rapForPTIndex]->Fill(deltaR);
    // 	  trigEff_deltaRM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDrM2);
    // 	  trigEff_deltaPhiM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDphiM2);
    // 	  trigEff_deltaEtaM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDetaM2);
    // 	  trigEff_distM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDistM2);
    // 	}
    //   }
    //   if(incrementTot){
    // 	totEff_deltaR_pT_rap[pTIndex][rapForPTIndex]->Fill(deltaR);
    // 	totEff_deltaRM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDrM2);
    // 	totEff_deltaPhiM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDphiM2);
    // 	totEff_deltaEtaM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDetaM2);
    // 	totEff_distM2_pT_rap[pTIndex][rapForPTIndex]->Fill(JpsiDistM2);
    //   }

    //   hGen2D_deltaPhiVsDeltaEta_pT_rap[pTIndex][rapForPTIndex]->Fill(deltaEta, deltaPhiDeg);
    //   if(incrementReco){
    // 	recoEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex][rapForPTIndex]->Fill(deltaEta, deltaPhiDeg);
    // 	if(incrementTrig)
    // 	  trigEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex][rapForPTIndex]->Fill(deltaEta, deltaPhiDeg);
    //   }
    //   if(incrementTot)
    // 	totEff2D_deltaPhiVsDeltaEta_pT_rap[pTIndex][rapForPTIndex]->Fill(deltaEta, deltaPhiDeg);
    // }

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
    printf("%d not a valid efficiency sample!!!\n", iEffSample);

  Int_t binX = hMuEff[iEff][iEffSample][CENTRAL]->GetXaxis()->FindBin(fabs(eta));
  Int_t binY = hMuEff[iEff][iEffSample][CENTRAL]->GetYaxis()->FindBin(pt);
  Double_t eff = hMuEff[iEff][iEffSample][CENTRAL]->GetBinContent(binX, binY);

  // printf("%s, efficiency for |eta|=%1.3f and pT=%1.2f GeV/c is %1.3f\n",
  // 	 effName[iEff], fabs(eta), pt, eff);

  if(eff > 1.) eff = 1.;
  else if(eff < 0.) eff = 0.;

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
  // for(int iEta = 0; iEta < kNbEtaBins; iEta++){
  for(int iEta = 0; iEta < kNbEtaBins[iEff]; iEta++){
    if(fabs(eta) < binsEta[iEff][iEta]){
      etaBin = iEta;
      break;
    }
  }

  Double_t eff = fMuEff_pT[iEff][iEffSample][etaBin]->Eval(pt);
  if(eff > 1.) eff = 1.;
  else if(eff < 0.) eff = 0.;

  //  printf("%s: eta %f --> mapped into bin %d... efficiency for pt = %f GeV/c is %f\n", effName[iEff], eta, etaBin, pt, eff);

  return eff;

}
