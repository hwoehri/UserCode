#define TnPDimuonEff_flatGen_cxx
#include "TnPDimuonEff_flatGen.h"
#include "calcPol.C"
#include "/Users/hwoehri/CMS/CMSSW/hWoehri/Polarization/macros/areMuonsInAcceptance.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>
#include "TRandom.h"
#include "TMath.h"
//#include "TF1.h"
#include "TEfficiency.h"

Int_t const kNbEff = 8;
Char_t *effFileNames[kNbEff] = {//(1)
                                "/Users/hwoehri/CMS/Work/TnP/LucaPerrozzi/28March2011/MuonTrackingEff_28March2011.root",
				//(2)
				"/Users/hwoehri/CMS/Work/TnP/Zongchang/29March2011/MuonIDEff_29March2011_fitted.root",
				//(3)
				"/Users/hwoehri/CMS/Work/TnP/Xianyou/25March2011/MuonQualEff_25March2011_fitted.root",
				//(4)
				"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L1L2EfficiencyMu0_Mu3_Track3_pt_abseta_seagull_looseCuts_runA.root", 
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L1L2EfficiencyMu0_Mu3_Track3_pt_abseta_seagull_looseCuts_runB.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L1L2EfficiencyMu0_Mu3_Track3_pt_abseta_seagull_tightCuts_runA.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L1L2EfficiencyMu0_Mu3_Track3_pt_abseta_seagull_tightCuts_runB.root",
				//(5)
				"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L3DoubleMu0_seagull_looseCuts_runA.root", 
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L3DoubleMu0_seagull_looseCuts_runB.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L3DoubleMu0_seagull_tightCuts_runA.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L3DoubleMu0_seagull_tightCuts_runB.root",
				//(6)
				"/Users/hwoehri/CMS/Work/TnP/Ilse/23June2011/Track_MuXTrack0_seagull_looseCuts_runAB1andB2.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/23June2011/Track_MuXTrack0_seagull_looseCuts_runB3.root",
				//(7)
				"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_looseCuts_runA.root", 
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_looseCuts_runB1.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_looseCuts_runB2.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_tightCuts_runA.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_tightCuts_runB1.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_tightCuts_runB2.root",
				//(8)
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/29March2011/HLTMuonTrack_Mu0_TkMu0_TM_runB1_29March2011.root",
				//"/Users/hwoehri/CMS/CMSSW/hWoehri/Efficiencies/macros/TMEffConditionalToHLT_MC_20May2011.root"};
				"/Users/hwoehri/CMS/Work/TnP/Ilse/23May2011/TkMu0L1L2-pt-abseta-runA-distM2gt120.root"};
//old versions:
//"/Users/hwoehri/CMS/Work/TnP/Herbert/13May2011/HLT_Track_Mu0TkMu0_Run1_13May2011_MCCorrected.root",
// "/Users/hwoehri/CMS/Work/TnP/Herbert/13May2011/HLT_Track_Mu0TkMu0_Run1_13May2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Luigi/17April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM2gt120_17April2011_fitted.root",
//"/Users/hwoehri/CMS/Work/TnP/Francesco/2April2011/L1L2_DMu0_TriggerEfficiencies_RUNA_distM2gt120_7April2011_fitted.root",
//"/Users/hwoehri/CMS/Work/TnP/Francesco/2April2011/L1L2_DMu0_TriggerEfficiencies_RUNB_distM2gt120_7April2011_fitted.root",
//"/Users/hwoehri/CMS/Work/TnP/Zongchang/19March2011/MuonIDEff_19March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Francesco/23March2011/L1L2_DMu0_TriggerEfficiencies_23March2011.root", //mix of RunA and B
//"/Users/hwoehri/CMS/Work/TnP/Luigi/25March2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_25March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Luigi/5April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM1gt150_7April2011_fitted.root",
//"/Users/hwoehri/CMS/Work/TnP/LucaPerrozzi/17March2011/MuonTrackingEff_17March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Xianyou/19March2011/MuonQualEff_19March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Francesco/18March2011/L1L2_DoubleMu0_TriggerEfficiencies_18March2011.root",
//"/Users/hwoehri/CMS/Work/TnP/Luigi/19March2011/L3_DoubleMu0_TriggerEfficiencies_19March2011.root"};
//"/Users/hwoehri/CMS/Work/TnP/Herbert/28March2011/HLT_Track_Mu0TkMu0_Run1.root",
enum {TrkEff, MuIDEff, MuQualEff, L1L2Eff, L3Eff, Trk_TkMu0, Mu_TkMu0, TMandHLT2};
Char_t *effName[kNbEff] = {"TrkEff", "MuIDEff", "MuQualEff", "L1L2Eff", "L3Eff", "Trk_TkMu0Eff", "Mu_TkMu0Eff", "TMandHLT2"};
Int_t const kNbEffSample = 3;
enum {DATA, MC, MCTRUTH};
Char_t *effSampleName[kNbEffSample] = {"DATA", "MC", "MCTRUTH"};
//
enum {CENTRAL, UPPER, LOWER};
TH2D *hMuEff[kNbEff][kNbEffSample][3];
Int_t const kNbMaxEtaBins = 14;
Int_t const kNbEtaBins[kNbEff] = {0,8,8,8,14,0,0,0};//number of TGraphs with parametrized pT diff. efficiencies
Float_t binsEta[kNbEff][kNbMaxEtaBins] = {{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},//needed to fetch the correct eta bin for the fitted pT turn on curves
					  {0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4,0.,0.,0.,0.,0.,0.},//MuID
					  {0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4,0.,0.,0.,0.,0.,0.},//MuQual
					  {0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4,0.,0.,0.,0.,0.,0.},//L1L2
					  {0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4},//L3
					  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},//Trk_TkMu0
					  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},//Mu_TkMu0
					  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},};//TMgivenHLT
TF1 *fMuEff_pT[kNbEff][kNbEffSample][kNbMaxEtaBins];
enum {LOOSE,TIGHT};//set of muon fiducial cuts
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
// TH1D *hGen_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *recoEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *trigEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *totEff_deltaR_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

// TH2D *hGen2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH2D *recoEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH2D *trigEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH2D *totEff2D_deltaPhiVsDeltaEta_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

// TH1D *hGen_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *recoEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *trigEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *totEff_deltaRM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

// TH1D *hGen_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *recoEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *trigEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *totEff_deltaPhiM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

// TH1D *hGen_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *recoEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *trigEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *totEff_deltaEtaM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

// TH1D *hGen_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *recoEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *trigEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
// TH1D *totEff_distM2_pT_rap[eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

Double_t GetEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt);
Double_t GetEfficiency_FromParametrization(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt);
//==============================================
void TnPDimuonEff_flatGen::Loop(Int_t effSample, Char_t *trigLabel, Bool_t rejectCowboys, Int_t muSelectionCuts)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nb = 0;
  Long64_t countRecEvent = 0;
  printf("number of entries = %d\n", (Int_t) nentries);

   //assign a weight according to a given pT and y distribution (taken from PYTHIA)
   //PYTHIA distributions fit with the macro "projectPTRap_Pythia.C(kTRUE)"
   //pT distribtuions do NOT show any rapidity dependence, while for pT > 10 GeV/c
   //there seems to be some rapidity dependence which is neglected in this parameterization
  //    TF1 *fPT = new TF1("fPT", "[0]*x*pow(1.+(1./([1]-2.))*x*x/[2],-[1])", 0., 50.);
  //    fPT->FixParameter(0, 0.61554); fPT->SetParName(0, "norm");
  //    fPT->FixParameter(1, 4.05847); fPT->SetParName(1, "beta");
  //    fPT->FixParameter(2, 19.93); fPT->SetParName(2, "<pT2> [GeV2]");
  TF1 *fRap = new TF1("fRap", "pol2", -2.5, 2.5);
  fRap->FixParameter(0, 0.8632);
  fRap->FixParameter(1, -7.311e-4);
  fRap->FixParameter(2, -3.041e-2);

  //alternatively, the pT distribution can be taken from the
  //preliminary ATLAS data shown in Moriond, 17th April
  TF1 *fPT = new TF1("fPT", "[0]*x*pow(1.+(1./([1]-2.))*x*x/[2],-[1])", 0., 50.);
  fPT->FixParameter(0, 1.); fPT->SetParName(0, "norm");
  fPT->FixParameter(1, 3.87004); fPT->SetParName(1, "beta");
  fPT->FixParameter(2, 15.5821); fPT->SetParName(2, "<pT2> [GeV2]");
  
  //alternatively, the pT distribution can be taken from CASCADE:
  //    TF1 *fPT = new TF1("fPT", "[0]*x*pow(1.+(1./([1]-2.))*x*x/[2],-[1])", 0., 50.);
  //    fPT->FixParameter(0, 1.18390); fPT->SetParName(0, "norm");
  //    fPT->FixParameter(1, 2.70764); fPT->SetParName(1, "beta");
  //    fPT->FixParameter(2, 12.0053); fPT->SetParName(2, "<pT2> [GeV2]");
  
   //for NP J/psi's use the distribution taken from the 
   //preliminary ATLAS data shown in Moriond, 17th April
  //    TF1 *fPT = new TF1("fPT", "[0]*x*pow(1.+(1./([1]-2.))*x*x/[2],-[1])", 0., 50.);
  //    fPT->FixParameter(0, 0.9); fPT->SetParName(0, "norm");
  //    fPT->FixParameter(1, 2.98180); fPT->SetParName(1, "beta");
  //    fPT->FixParameter(2, 28.8277); fPT->SetParName(2, "<pT2> [GeV2]");

  Bool_t incrementReco, incrementTrig, incrementTot;

  //loop over the events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if(jentry % 100000 == 0) printf("event %d\n", (Int_t) jentry);
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    //=================================================
    //get the individual muons that are affected by FSR
    Double_t enMuPos = sqrt(muPos_Px*muPos_Px + muPos_Py*muPos_Py + muPos_Pz*muPos_Pz + eff::muMass*eff::muMass);
    Double_t enMuNeg = sqrt(muNeg_Px*muNeg_Px + muNeg_Py*muNeg_Py + muNeg_Pz*muNeg_Pz + eff::muMass*eff::muMass);
    TLorentzVector *muPos = new TLorentzVector();
    TLorentzVector *muNeg = new TLorentzVector();
    muPos->SetPxPyPzE(muPos_Px, muPos_Py, muPos_Pz, enMuPos);
    muNeg->SetPxPyPzE(muNeg_Px, muNeg_Py, muNeg_Pz, enMuNeg);
    Double_t etaMuPos = muPos->PseudoRapidity();
    Double_t etaMuNeg = muNeg->PseudoRapidity();
    Double_t pTMuPos = muPos->Pt();
    Double_t pTMuNeg = muNeg->Pt();
    Double_t pMuPos = muPos->P();
    Double_t pMuNeg = muNeg->P();

    //1.) apply the single muon kinematical cuts:
    if(rejectCowboys)
      if((muNeg->Phi() - muPos->Phi()) < 0.)
	continue;
    
    //apply different fiducial cuts for the muon matched to the HLT leg or the TM leg
    //for simplicity we decided to apply the stronger cuts to the higher pT muon and not
    //strictly to the HLT muon
    Bool_t muonsInAcc = kFALSE;
    if(pTMuPos > pTMuNeg)
      muonsInAcc = areMuonsInAcceptance(muSelectionCuts, pTMuPos, etaMuPos, pTMuNeg, etaMuNeg); //code needs to be adjusted for the DoubleMu0 trigger!!!!
    else
      muonsInAcc = areMuonsInAcceptance(muSelectionCuts, pTMuNeg, etaMuNeg, pTMuPos, etaMuPos); //code needs to be adjusted for the DoubleMu0 trigger!!!!
    if(!muonsInAcc)
      continue;

    //reject furthermore all events in which one of the muons has a pT smaller than 3 GeV/c
    if(pTMuPos < 3.0 || pTMuNeg < 3.0)
      continue;

    //build the invariant mass, pt, ... of the two muons
    TLorentzVector *oniaFSR = new TLorentzVector();
    *oniaFSR = *(muPos) + *(muNeg);
    Double_t rapOniaFSR = oniaFSR->Rapidity();
    Double_t pTOniaFSR = oniaFSR->Pt();
    Double_t massOniaFSR = oniaFSR->M();
    Double_t phiOniaFSR = oniaFSR->Phi();

    //2.) build the dimuon fiducial cuts:
    if(fabs(rapOniaFSR) > eff::rapMax)
      continue;

    Int_t rapIndex = -1;
    for(int iRap = 0; iRap < 2*eff::kNbRapBins; iRap++){
      if(rapOniaFSR > eff::rapRange[iRap] && rapOniaFSR < eff::rapRange[iRap+1]){
	rapIndex = iRap;
	break;
      }
    }
    Int_t rapForPTIndex = -1;
    for(int iRap = 0; iRap < eff::kNbRapForPTBins; iRap++){
      if(TMath::Abs(rapOniaFSR) > eff::rapForPTRange[iRap] && 
	 TMath::Abs(rapOniaFSR) < eff::rapForPTRange[iRap+1]){
	rapForPTIndex = iRap+1;
	break;
      }
    }
    Int_t pTIndex = -1;
    for(int iPT = 0; iPT < eff::kNbPTBins[rapForPTIndex]; iPT++){
      if(pTOniaFSR > eff::pTRange[rapForPTIndex][iPT] && pTOniaFSR < eff::pTRange[rapForPTIndex][iPT+1]){
	pTIndex = iPT+1;
	break;
      }
    }
    Int_t rapIntegratedPTIndex = -1;
    for(int iPT = 0; iPT < eff::kNbPTBins[0]; iPT++){
      if(pTOniaFSR > eff::pTRange[0][iPT] && pTOniaFSR < eff::pTRange[0][iPT+1]){
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
    // Bool_t usePTFit = kTRUE; //alternative to the T&P histograms use fitted pT differential efficiency
    Bool_t usePTFit = kFALSE; //alternative to the T&P histograms use fitted pT differential efficiency

    Double_t epsTrack_Pos  = GetEfficiency(TrkEff, effSample, etaMuPos, pTMuPos); 
    Double_t epsTrack_Neg  = GetEfficiency(TrkEff, effSample, etaMuNeg, pTMuNeg); 

    Double_t epsMuonID_Pos, epsQual_Pos, epsMuonID_Neg, epsQual_Neg;
    if(usePTFit){
      epsMuonID_Pos = GetEfficiency_FromParametrization(MuIDEff, effSample, etaMuPos, pTMuPos);
      epsQual_Pos   = GetEfficiency_FromParametrization(MuQualEff, effSample, etaMuPos, pTMuPos);
      epsMuonID_Neg = GetEfficiency_FromParametrization(MuIDEff, effSample, etaMuNeg, pTMuNeg);
      epsQual_Neg   = GetEfficiency_FromParametrization(MuQualEff, effSample, etaMuNeg, pTMuNeg);
    }
    else{
      epsMuonID_Pos = GetEfficiency(MuIDEff, effSample, etaMuPos, pTMuPos);
      epsQual_Pos   = GetEfficiency(MuQualEff, effSample, etaMuPos, pTMuPos);
      epsMuonID_Neg = GetEfficiency(MuIDEff, effSample, etaMuNeg, pTMuNeg);
      epsQual_Neg   = GetEfficiency(MuQualEff, effSample, etaMuNeg, pTMuNeg);
    }
    Double_t recoEff_Pos = epsTrack_Pos * epsMuonID_Pos * epsQual_Pos;
    Double_t recoEff_Neg = epsTrack_Neg * epsMuonID_Neg * epsQual_Neg;
    Double_t recoEff = recoEff_Pos * recoEff_Neg;

    //the indiv. histograms will be filled depending on the
    //assigned probability
    Double_t randNb = gRandom->Uniform();
    if(recoEff > randNb)
      incrementReco = kTRUE;
    else
      incrementReco = kFALSE;

    //printf("recoEff %1.3f --> accepted = %d\n", recoEff, incrementReco);

    Double_t epsL1L2Trig_Pos, epsL3Trig_Pos, epsL1L2Trig_Neg, epsL3Trig_Neg;
    Double_t epsTMandHLT_Pos, epsTMandHLT_Neg;
    Double_t trigEff = 0.;
    if(incrementReco){//check the trigger efficiency only for those events that pass RECO

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

	//the following code is used from 27th of May 2011 onwards:
	if(usePTFit){
	  epsL1L2Trig_Pos  = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuPos, pTMuPos);
	  epsL3Trig_Pos    = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuPos, pTMuPos);
	  epsL1L2Trig_Neg  = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuNeg, pTMuNeg);
	  epsL3Trig_Neg    = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuNeg, pTMuNeg);
	}
	else{
	  epsL1L2Trig_Pos  = GetEfficiency(L1L2Eff, effSample, etaMuPos, pTMuPos);
	  epsL3Trig_Pos    = GetEfficiency(L3Eff, effSample, etaMuPos, pTMuPos);
	  epsL1L2Trig_Neg  = GetEfficiency(L1L2Eff, effSample, etaMuNeg, pTMuNeg);
	  epsL3Trig_Neg    = GetEfficiency(L3Eff, effSample, etaMuNeg, pTMuNeg);
	}
	Double_t epsTM1_Neg  = GetEfficiency(Trk_TkMu0, effSample, etaMuNeg, pTMuNeg);
	Double_t epsTM2_Neg  = GetEfficiency(Mu_TkMu0, effSample, etaMuNeg, pTMuNeg);
	Double_t epsTM1_Pos = GetEfficiency(Trk_TkMu0, effSample, etaMuPos, pTMuPos);
	Double_t epsTM2_Pos = GetEfficiency(Mu_TkMu0, effSample, etaMuPos, pTMuPos);

	if(pTMuPos > pTMuNeg)
	  trigEff = epsL1L2Trig_Pos * epsL3Trig_Pos * epsTM1_Neg * epsTM2_Neg;
	else
	  trigEff = epsL1L2Trig_Neg * epsL3Trig_Neg * epsTM1_Pos * epsTM2_Pos;
      }//end of low pT J/psi trigger efficiency
    }

    if(trigEff > randNb)
      incrementTrig = kTRUE;
    else
      incrementTrig = kFALSE;

    //printf("trigEff %f --> accepted %d\n", trigEff, incrementTrig);

    Double_t totEff = trigEff * recoEff;
    if(totEff > randNb) 
      incrementTot = kTRUE;
    else
      incrementTot = kFALSE;

    //printf("totEff %f --> accepted %d\n", totEff, incrementTot);

    //fill the pT, y, phi 1D and 2D efficiency histos
    recoEff_pT->Fill(incrementReco, pTOniaFSR);
    recoEff_y->Fill(incrementReco, rapOniaFSR);
    recoEff_phi->Fill(incrementReco, phiOniaFSR);
    recoEff2D_pT_rapNP->Fill(incrementReco, rapOniaFSR, pTOniaFSR);
    recoEff2D_pT_rap->Fill(incrementReco, fabs(rapOniaFSR), pTOniaFSR);

    trigEff_pT->Fill(incrementTrig, pTOniaFSR);
    trigEff_y->Fill(incrementTrig, rapOniaFSR);
    trigEff_phi->Fill(incrementTrig, phiOniaFSR);
    trigEff2D_pT_rapNP->Fill(incrementTrig, rapOniaFSR, pTOniaFSR);
    trigEff2D_pT_rap->Fill(incrementTrig, fabs(rapOniaFSR), pTOniaFSR);

    totEff_pT->Fill(incrementTot, pTOniaFSR);
    totEff_y->Fill(incrementTot, rapOniaFSR);
    totEff_phi->Fill(incrementTot, phiOniaFSR);
    totEff2D_pT_rapNP->Fill(incrementTot, rapOniaFSR, pTOniaFSR);
    totEff2D_pT_rap->Fill(incrementTot, fabs(rapOniaFSR), pTOniaFSR);

    //fill the eff. histos for all the different frames
    for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){

      //histos for neg. and pos. rapidity separately:
      if(rapIndex >= 0){
	recoEff2D_pol_pT_rapNP[iFrame][0][rapIndex]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	trigEff2D_pol_pT_rapNP[iFrame][0][rapIndex]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	totEff2D_pol_pT_rapNP[iFrame][0][rapIndex]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(pTIndex > 0 && rapIndex >= 0){
	recoEff2D_pol_pT_rapNP[iFrame][pTIndex][rapIndex]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	trigEff2D_pol_pT_rapNP[iFrame][pTIndex][rapIndex]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	totEff2D_pol_pT_rapNP[iFrame][pTIndex][rapIndex]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
      }
	
      //histos taking together +y and -y
      recoEff2D_pol_pT_rap[iFrame][0][0]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
      trigEff2D_pol_pT_rap[iFrame][0][0]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
      totEff2D_pol_pT_rap[iFrame][0][0]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
      if(rapIntegratedPTIndex > 0){
	recoEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	trigEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	totEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(rapForPTIndex > 0){
	recoEff2D_pol_pT_rap[iFrame][0][rapForPTIndex]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	trigEff2D_pol_pT_rap[iFrame][0][rapForPTIndex]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	totEff2D_pol_pT_rap[iFrame][0][rapForPTIndex]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
      }
      if(pTIndex > 0 && rapForPTIndex > 0){
	recoEff2D_pol_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	trigEff2D_pol_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	totEff2D_pol_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
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

      recoEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(incrementReco, thetaAdjusted, phiFolded);
      trigEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(incrementTrig, thetaAdjusted, phiFolded);
      totEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(incrementTot, thetaAdjusted, phiFolded);

      if(rapIntegratedPTIndex > 0){
	recoEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex][0]->Fill(incrementReco, thetaAdjusted, phiFolded);
	trigEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex][0]->Fill(incrementTrig, thetaAdjusted, phiFolded);
	totEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex][0]->Fill(incrementTot, thetaAdjusted, phiFolded);
      }
      if(rapForPTIndex > 0){
	recoEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex]->Fill(incrementReco, thetaAdjusted, phiFolded);
	trigEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex]->Fill(incrementTrig, thetaAdjusted, phiFolded);
	totEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex]->Fill(incrementTot, thetaAdjusted, phiFolded);
      }
      if(pTIndex > 0 && rapForPTIndex > 0){
	recoEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex][rapForPTIndex]->Fill(incrementReco, thetaAdjusted, phiFolded);
	trigEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex][rapForPTIndex]->Fill(incrementTrig, thetaAdjusted, phiFolded);
	totEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex][rapForPTIndex]->Fill(incrementTot, thetaAdjusted, phiFolded);
      }
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
