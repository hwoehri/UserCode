#define MCTnPEff_cxx
#include "MCTnPEff.h"
#include "calcPol.C"
#include "isMuonInAcceptance.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>
#include "TRandom.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TGraph2D.h"

//switch here between a total single muon efficiency:
Char_t *dimuEffName = "DimuVertexEff";
Char_t *dimuEffFileName = "/Users/hwoehri/CMS/Work/TnP/2011/Linlin/1Dec2011/DimuonVtxEff_Dimuon0Jpsi_cosTheta_Phi_TrkCuts80_CS_01Dec2011.root";
enum {SingleMuEff};
Bool_t useIndivEff = kFALSE;
Int_t const kNbEff = 1;
Char_t *effName[kNbEff] = {"SingleMuEff"};
Char_t *effFileNames[kNbEff] = {"singleMuTruthEff_16July2012_40GeVrap1_2pT100GeV_FineBins200MeV.root"};//singleMu MCTruth efficiency
//Char_t *effFileNames[kNbEff] = {"singleMuTruthEff_18Jan2012_40GeVrap1_2pT100GeV_EtaCut_FineBins200MeV.root"};//singleMu MCTruth efficiency
//Char_t *effFileNames[kNbEff] = {"singleMuTruthEff_17Jan2012_40GeVrap1_2pT100GeV_EtaCut_FineBins100MeVand0p1Eta.root"};//singleMu MCTruth efficiency
//Char_t *effFileNames[kNbEff] = {"singleMuTruthEff_17Jan2012_40GeVrap1_2pT100GeV_EtaCut_FineBins100MeV.root"};//singleMu MCTruth efficiency
//Char_t *effFileNames[kNbEff] = {"singleMuTruthEff_17Jan2012_40GeVrap1_2pT100GeV_EtaCut_FineBins250MeV.root"};//singleMu MCTruth efficiency
//Char_t *effFileNames[kNbEff] = {"singleMuTruthEff_17Jan2012_40GeVrap1_2pT100GeV_EtaCut.root"};//singleMu MCTruth efficiency
//Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/8Dec2011/EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Trk80Cuts_19Nov2011.root"};
//Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/6Feb2012/EfficiencyProduct_combinedMC_Trk80Cuts_6Feb2011.root"};
//Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/6Feb2012/singleMuonEff_combinedMC_Trk80Cuts_6Feb2011.root"};//SingleMuEff a la Matt
//Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Linlin/7Dec2011/singleMuonEfficiency_ProbeTrackMatched_data_mc_pt_abseta_tracker80Cuts_7Dec2011.root"};//SingleMuEff a la Matt
//Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Linlin/2Dec2011/SingleMuEff_Dimuon10Jpsi_ProbeTrackMatched_data_mc_pt_abseta_tracker80Cuts_02Dec2011.root"};//SingleMuEff a la Matt (L1*L2*L3 trigger eff only)
//Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/7Dec2011/EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Trk80Cuts_19Nov2011.root"};
//Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/23Nov2011/EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Trk80Cuts_19Nov2011.root"};

//... or by using the inividual single muon efficiencies
// Bool_t useIndivEff = kTRUE;
// Int_t const kNbEff = 4;
// Char_t *effName[kNbEff] = {"MuIDEff", "MuQualEff", "L1L2Eff", "L3Eff"};
// Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/30Nov2011/MuonID_pt_abseta_runA_inclMay10_Trk80Cuts_19Nov2011_corrected.root",
// 				"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/30Nov2011/MuonQual_pt_abseta_runA_inclMay10_Trk80Cuts_19Nov2011_corrected.root",
// 				"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/30Nov2011/L1L2Dimuon0Jpsi_pt_abseta_seagulls_run1_Trk80Cuts_19Nov2011_corrected.root",
// 				"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/30Nov2011/L3Dimuon0Jpsi_pt_abseta_seagulls_run1_Trk80Cuts_19Nov2011_corrected.root"};

//Tracker50 muon cuts:
// Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/1Nov2011/MuonID_pt_abseta_TrkCuts_20Oct2011_corrected.root",
// 				"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/1Nov2011/MuonQual_pt_abseta_TrkCuts_20Oct2011_corrected.root",
// 				"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/1Nov2011/L1L2Dimuon0Jpsi_pt_abseta_seagulls_run1_TrkCuts_20Oct2011_corrected.root",
// 				"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/1Nov2011/L3Dimuon0Jpsi_pt_abseta_seagulls_TrkCuts_20Oct2011_corrected.root"};
enum {MuIDEff, MuQualEff, L1L2Eff, L3Eff};
Int_t const kNbEffSample = 3;
enum {DATA, MC, MCTRUTH};
Char_t *effSampleName[kNbEffSample] = {"DATA", "MC", "MCTRUTH"};
//
enum {JPSI, PSIP, UPS1S, UPS2S, UPS3S};
//
enum {CENTRAL, UPPER, LOWER};
TH2D *hMuEff[kNbEff][kNbEffSample][3];
TH2D *hMuEff_smoothed[kNbEff][kNbEffSample];
TH2D *hDiMuEff[kNbEffSample][3];
TEfficiency *tMuEff[kNbEff][kNbEffSample];
Int_t const kNbMaxEtaBins = 10;
Int_t const kNbEtaBins[kNbEff] = {10};//number of TGraphs with parametrized pT diff. efficiencies
// Float_t binsEta[kNbEff][kNbMaxEtaBins+1] = {{0, 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.1, 2.4},
// 					    {0, 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.1, 2.4},
// 					    {0, 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.1, 2.4},
// 					    {0, 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.1, 2.4}};//needed to fetch the correct eta bin for the fitted pT 
//turn on curves
Float_t binsEta[kNbEff][kNbMaxEtaBins+1] = {{0, 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.1, 2.4}};//needed to fetch the correct eta bin for the fitted pT turn on curves

//Tracker50 muon cuts:
// Int_t const kNbEtaBins[kNbEff] = {8};//number of TGraphs with parametrized pT diff. efficiencies
// Float_t binsEta[kNbEff][kNbMaxEtaBins+1] = {{0, 0.2, 0.3, 0.6, 0.8, 1.2, 1.6, 2.1, 2.4},
// 					    {0, 0.2, 0.3, 0.6, 0.8, 1.2, 1.6, 2.1, 2.4},
// 					    {0, 0.2, 0.3, 0.6, 0.8, 1.2, 1.6, 2.1, 2.4},
// 					    {0, 0.2, 0.3, 0.6, 0.8, 1.2, 1.6, 2.1, 2.4}};//needed to fetch the correct eta bin for the fitted pT turn on curves

TGraph2D *gEff2D[kNbEff][kNbEffSample];
TF1 *fMuEff_pT[kNbEff][kNbEffSample][kNbMaxEtaBins];
enum {LOOSE,TIGHT};//set of muon fiducial cuts
//=================================================
//1D histos versus J/psi pT, phi and y
//=================================================
//histos for neg. and pos. rapidity separately:
// TH1D *hGen_pT, *hGen_y, *hGen_phi;
// TH1D *recoEff_pT, *recoEff_y, *recoEff_phi;
// TH1D *trigEff_pT, *trigEff_y, *trigEff_phi;
// TH1D *totEff_pT, *totEff_y, *totEff_phi;
TEfficiency *recoEff_pT, *recoEff_y, *recoEff_phi;
TEfficiency *trigEff_pT, *trigEff_y, *trigEff_phi;
TEfficiency *totEff_pT, *totEff_y, *totEff_phi;
//=================================================
//2D histos versus J/psi pT and y
//=================================================
//histos for neg. and pos. rapidity separately:
// TH2D *hGen2D_pT_rapNP, *hGen2D_pT_rap;
// TH2D *recoEff2D_pT_rapNP, *recoEff2D_pT_rap;
// TH2D *trigEff2D_pT_rapNP, *trigEff2D_pT_rap;
// TH2D *totEff2D_pT_rapNP, *totEff2D_pT_rap;
TEfficiency *trigEff2D_pT_rap, *recoEff2D_pT_rap, *totEff2D_pT_rap;
TEfficiency *trigEff2D_pT_rapNP, *recoEff2D_pT_rapNP, *totEff2D_pT_rapNP;
TEfficiency *recoEff_phiPol[eff::kNbFrames], *recoEff_cosTheta[eff::kNbFrames];
TEfficiency *recoEff2D_cosTheta_phiPol[eff::kNbFrames];
TEfficiency *trigEff_phiPol[eff::kNbFrames], *trigEff_cosTheta[eff::kNbFrames];
TEfficiency *trigEff2D_cosTheta_phiPol[eff::kNbFrames];
TEfficiency *totEff_phiPol[eff::kNbFrames], *totEff_cosTheta[eff::kNbFrames];
TEfficiency *totEff2D_cosTheta_phiPol[eff::kNbFrames];
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
TEfficiency *recoEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][2*eff::kNbRapBins+1];
TEfficiency *recoEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *recoEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TEfficiency *recoEff_phiPol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *recoEff_cosTheta_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TEfficiency *trigEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][2*eff::kNbRapBins+1];
TEfficiency *trigEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TEfficiency *trigEff_phiPol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *trigEff_cosTheta_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TEfficiency *totEff2D_pol_pT_rapNP[eff::kNbFrames][eff::kNbPTMaxBins+1][2*eff::kNbRapBins+1];
TEfficiency *totEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff2D_pol_pT_rap_phiFolded[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

TEfficiency *totEff_phiPol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *totEff_cosTheta_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
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

TH2F *hCorrRECO;
Double_t pTMin = 5.;//5... Upsilon; 10... Jpsi
Double_t GetDimuEfficiency(Int_t iEffSample, Double_t cosTheta, Double_t phi);
Double_t GetEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt, Bool_t use2DGraph = kFALSE);
Double_t GetTEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt);
Double_t GetEfficiency_FromParametrization(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt);
//==============================================
void MCTnPEff::Loop(Int_t effSample, Int_t resonance, Bool_t rejectCowboys, Bool_t use2DGraph, Bool_t useTEfficiency, Int_t useEventNTimes)
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

  Bool_t incrementReco, incrementTrig, incrementTot;

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
    Double_t phiMuPos_Gen = muPos_Gen->Phi();
    Double_t phiMuNeg_Gen = muNeg_Gen->Phi();

    Double_t deltaPhi = phiMuNeg_Gen - phiMuPos_Gen;
    if(rejectCowboys){
      if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
      else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
      if(deltaPhi < 0.) //reject cowboys
	continue;
    }

    // //if(!(isMuonInAcceptance(LOOSE, pTMuPos_Gen, etaMuPos_Gen) && isMuonInAcceptance(LOOSE, pTMuNeg_Gen, etaMuNeg_Gen)))
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

    // Bool_t usePTFit = kTRUE; //alternative to the T&P histograms use fitted pT differential efficiency
    Bool_t usePTFit = kFALSE; //alternative to the T&P histograms use fitted pT differential efficiency

    Double_t randNb = 0.;
    Double_t totEff = 0.;

    for(int iN = 0; iN < useEventNTimes; iN++){

      randNb = gRandom->Uniform();
      if(!useIndivEff){
      
	totEff = 0.99*0.99; //tracking efficiency
	//totEff = 1.0 * 1.0;
      
	if(usePTFit){
	  if(!useTEfficiency){
	    totEff *= GetEfficiency_FromParametrization(SingleMuEff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	    totEff *= GetEfficiency_FromParametrization(SingleMuEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	  }
	}
	else{
	  if(!useTEfficiency){
	    totEff *= GetEfficiency(SingleMuEff, effSample, etaMuPos_Gen, pTMuPos_Gen, use2DGraph);
	    totEff *= GetEfficiency(SingleMuEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen, use2DGraph);
	  }
	  else{
	    totEff *= GetTEfficiency(SingleMuEff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	    totEff *= GetTEfficiency(SingleMuEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	  }
	}
	// //dimuon vertexing cut NOT included in the (single) muon quality cuts...
	// if(JpsiVprob < 0.01)
	// 	totEff = 0.;
      }
      else{
	Double_t epsTrack_Pos = 0.99;
	Double_t epsTrack_Neg = 0.99;
	// Double_t epsTrack_Pos  = 1.0;
	// Double_t epsTrack_Neg = 1.0;
	Double_t epsMuonID_Pos, epsQual_Pos, epsMuonID_Neg, epsQual_Neg;
	if(usePTFit){
	  epsMuonID_Pos = GetEfficiency_FromParametrization(MuIDEff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	  epsQual_Pos   = GetEfficiency_FromParametrization(MuQualEff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	  epsMuonID_Neg = GetEfficiency_FromParametrization(MuIDEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	  epsQual_Neg   = GetEfficiency_FromParametrization(MuQualEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	}
	else{
	  epsMuonID_Pos = GetEfficiency(MuIDEff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	  epsQual_Pos   = GetEfficiency(MuQualEff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	  epsMuonID_Neg = GetEfficiency(MuIDEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	  epsQual_Neg   = GetEfficiency(MuQualEff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	}
	Double_t recoEff_Pos = epsTrack_Pos * epsMuonID_Pos * epsQual_Pos;
	Double_t recoEff_Neg = epsTrack_Neg * epsMuonID_Neg * epsQual_Neg;
	Double_t recoEff = recoEff_Pos * recoEff_Neg;
	
	// //dimuon vertexing cut NOT included in the (single) muon quality cuts...
	// if(JpsiVprob < 0.01)
	// 	recoEff = 0.;
	
	//the indiv. histograms will be filled depending on the
	//assigned probability
	if(recoEff > randNb)
	  incrementReco = kTRUE;
	else
	  incrementReco = kFALSE;

	//test: check whether this particular event was really reconstructed...
	Bool_t recoPassed = kFALSE;
	//      if(onia->Pt() < 990. && fabs(onia->Rapidity()) < eff::rapMax && JpsiVprob > 0.01) recoPassed = kTRUE;
	if(onia->Pt() < 990. && fabs(onia->Rapidity()) < eff::rapMax) recoPassed = kTRUE;
	hCorrRECO->Fill(recoPassed, incrementReco);

	Double_t epsL1L2Trig_Pos, epsL3Trig_Pos, epsL1L2Trig_Neg, epsL3Trig_Neg;
	Double_t trigEff = 0.;
	// if(recoPassed){//check the trigger efficiency only for those events that pass RECO
	incrementTrig = kFALSE;
	if(incrementReco){

	  //if(strncmp("HLT_Dimuon10_Jpsi_Barrel", trigLabel, 24) == 0 || strncmp("HLT_Dimuon5_Upsilon_Barrel", trigLabel, 26) == 0){
	  //HLT_Dimuon10_Jpsi_Barrel_v3... 1.4E33
	  //HLT_Dimuon10_Jpsi_Barrel_v5... no cowboys
	  //HLT_Dimuon10_Jpsi_Barrel_v6... L1DoubleMu0_HighQ
	  //HLT_Dimuon13_Jpsi_Barrel_v1... L1DoubleMu0_HighQ
	  //HLT_Dimuon5_Upsilon_Barrel_v3... 1.4E33
	  //HLT_Dimuon5_Upsilon_Barrel_v5... no cowboys
	  //HLT_Dimuon7_Upsilon_Barrel_v1... L1DoubleMu0_HighQ
	  //HLT_Dimuon9_Upsilon_Barrel_v1... L1DoubleMu0_HighQ

	  if(usePTFit){
	    epsL1L2Trig_Pos = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	    epsL3Trig_Pos   = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	    epsL1L2Trig_Neg = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	    epsL3Trig_Neg   = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	  }
	  else{
	    epsL1L2Trig_Pos = GetEfficiency(L1L2Eff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	    epsL3Trig_Pos   = GetEfficiency(L3Eff, effSample, etaMuPos_Gen, pTMuPos_Gen);
	    epsL1L2Trig_Neg = GetEfficiency(L1L2Eff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	    epsL3Trig_Neg   = GetEfficiency(L3Eff, effSample, etaMuNeg_Gen, pTMuNeg_Gen);
	  }
	  trigEff = epsL1L2Trig_Pos * epsL3Trig_Pos * epsL1L2Trig_Neg * epsL3Trig_Neg;
	  //}
	}
	if(trigEff > randNb)
	  incrementTrig = kTRUE;
	else
	  incrementTrig = kFALSE;

	totEff = trigEff * recoEff;
      }//useIndivEff
      // //account for the dimuon vertexing efficiency:
      // totEff *= GetDimuEfficiency(effSample, thisCosTh[eff::CS], thisPhi[eff::CS]);

      if(totEff > randNb) 
	incrementTot = kTRUE;
      else
	incrementTot = kFALSE;

      if(useIndivEff){
	recoEff_pT->Fill(incrementReco, onia_Gen_pt);
	// if(onia_Gen_pt > 10. && onia_Gen_pt < 25. && onia_Gen_rap > -2.0 && onia_Gen_rap < 2.0){
	if(onia_Gen_pt > pTMin){
	  recoEff_y->Fill(incrementReco, onia_Gen_rap);
	  recoEff_phi->Fill(incrementReco, onia_Gen_phi);
	}
	recoEff2D_pT_rapNP->Fill(incrementReco, onia_Gen_rap, onia_Gen_pt);
	recoEff2D_pT_rap->Fill(incrementReco, fabs(onia_Gen_rap), onia_Gen_pt);
    
	//      if(incrementReco){ //calculate the trigger efficiency only for events that pass RECO
	trigEff_pT->Fill(incrementTrig, onia_Gen_pt);
	//if(onia_Gen_pt > 10. && onia_Gen_pt < 25. && onia_Gen_rap > -2.0 && onia_Gen_rap < 2.0){
	if(onia_Gen_pt > pTMin){
	  trigEff_y->Fill(incrementTrig, onia_Gen_rap);
	  trigEff_phi->Fill(incrementTrig, onia_Gen_phi);
	}
	trigEff2D_pT_rapNP->Fill(incrementTrig, onia_Gen_rap, onia_Gen_pt);
	trigEff2D_pT_rap->Fill(incrementTrig, fabs(onia_Gen_rap), onia_Gen_pt);
	//}
      }

      totEff_pT->Fill(incrementTot, onia_Gen_pt);

      //if(onia_Gen_pt > 10. && onia_Gen_pt < 25. && onia_Gen_rap > -2.0 && onia_Gen_rap < 2.0){
      if(onia_Gen_pt > pTMin){
	totEff_y->Fill(incrementTot, onia_Gen_rap);
	totEff_phi->Fill(incrementTot, onia_Gen_phi);
      }
      totEff2D_pT_rapNP->Fill(incrementTot, onia_Gen_rap, onia_Gen_pt);
      totEff2D_pT_rap->Fill(incrementTot, fabs(onia_Gen_rap), onia_Gen_pt);
      
      //fill the eff. histos for all the different frames
      for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){

	//if(onia_Gen_pt > 10. && onia_Gen_pt < 25. && onia_Gen_rap > -2.0 && onia_Gen_rap < 2.0){
	if(onia_Gen_pt > pTMin){
	  if(useIndivEff){
	    recoEff_cosTheta[iFrame]->Fill(incrementReco, thisCosTh[iFrame]);
	    recoEff_phiPol[iFrame]->Fill(incrementReco, thisPhi[iFrame]);
	    recoEff2D_cosTheta_phiPol[iFrame]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	    
	    recoEff_phiPol_pT_rap[iFrame][0][0]->Fill(incrementReco, thisPhi[iFrame]);
	    recoEff_cosTheta_pT_rap[iFrame][0][0]->Fill(incrementReco, thisCosTh[iFrame]);
	  
	    // if(incrementReco){
	    trigEff_cosTheta[iFrame]->Fill(incrementTrig, thisCosTh[iFrame]);
	    trigEff_phiPol[iFrame]->Fill(incrementTrig, thisPhi[iFrame]);
	    trigEff2D_cosTheta_phiPol[iFrame]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	    
	    trigEff_phiPol_pT_rap[iFrame][0][0]->Fill(incrementTrig, thisPhi[iFrame]);
	    trigEff_cosTheta_pT_rap[iFrame][0][0]->Fill(incrementTrig, thisCosTh[iFrame]);
	    // }
	  }

	  totEff_cosTheta[iFrame]->Fill(incrementTot, thisCosTh[iFrame]);
	  totEff_phiPol[iFrame]->Fill(incrementTot, thisPhi[iFrame]);
	  totEff2D_cosTheta_phiPol[iFrame]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
	  
	  totEff_phiPol_pT_rap[iFrame][0][0]->Fill(incrementTot, thisPhi[iFrame]);
	  totEff_cosTheta_pT_rap[iFrame][0][0]->Fill(incrementTot, thisCosTh[iFrame]);
	}
	//histos for neg. and pos. rapidity separately:
	if(rapIndex_Gen >= 0){

	  if(useIndivEff){
	    recoEff2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	    // if(incrementReco)
	    trigEff2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	  }
	  totEff2D_pol_pT_rapNP[iFrame][0][rapIndex_Gen]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
	}
	if(pTIndex_Gen > 0 && rapIndex_Gen >= 0){
	  
	  if(useIndivEff){
	    recoEff2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	    // if(incrementReco)
	    trigEff2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	  }
	  totEff2D_pol_pT_rapNP[iFrame][pTIndex_Gen][rapIndex_Gen]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
	}
	
	//histos taking together +y and -y
	if(useIndivEff){
	  recoEff2D_pol_pT_rap[iFrame][0][0]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	  // if(incrementReco)
	  trigEff2D_pol_pT_rap[iFrame][0][0]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	}
	totEff2D_pol_pT_rap[iFrame][0][0]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
	if(rapIntegratedPTIndex_Gen > 0){
	  if(useIndivEff){
	    recoEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	    // if(incrementReco)
	    trigEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	  }
	  totEff2D_pol_pT_rap[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
	}
	if(rapForPTIndex_Gen > 0){
	  if(useIndivEff){
	    recoEff2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	    // if(incrementReco)
	    trigEff2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	  }
	  totEff2D_pol_pT_rap[iFrame][0][rapForPTIndex_Gen]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
	}
	if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0){
	  if(useIndivEff){
	    recoEff2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementReco, thisCosTh[iFrame], thisPhi[iFrame]);
	    recoEff_phiPol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementReco, thisPhi[iFrame]);
	    recoEff_cosTheta_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementReco, thisCosTh[iFrame]);
	    
	    // if(incrementReco)
	    trigEff2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementTrig, thisCosTh[iFrame], thisPhi[iFrame]);
	    
	    trigEff_phiPol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementTrig, thisPhi[iFrame]);
	    trigEff_cosTheta_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementTrig, thisCosTh[iFrame]);
	  }
	  totEff2D_pol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementTot, thisCosTh[iFrame], thisPhi[iFrame]);
	  
	  totEff_phiPol_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementTot, thisPhi[iFrame]);
	  totEff_cosTheta_pT_rap[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementTot, thisCosTh[iFrame]);
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
	if(useIndivEff){
	  recoEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(incrementReco, thetaAdjusted, phiFolded);
	  // if(incrementReco)
	  trigEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(incrementTrig, thetaAdjusted, phiFolded);
	}
	totEff2D_pol_pT_rap_phiFolded[iFrame][0][0]->Fill(incrementTot, thetaAdjusted, phiFolded);
	if(rapIntegratedPTIndex_Gen > 0){
	  if(useIndivEff){
	    recoEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(incrementReco, thetaAdjusted, phiFolded);
	    // if(incrementReco)
	    trigEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(incrementTrig, thetaAdjusted, phiFolded);
	  }
	  totEff2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex_Gen][0]->Fill(incrementTot, thetaAdjusted, phiFolded);
	}
	if(rapForPTIndex_Gen > 0){
	  if(useIndivEff){
	    recoEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(incrementReco, thetaAdjusted, phiFolded);
	    // if(incrementReco)
	    trigEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(incrementTrig, thetaAdjusted, phiFolded);
	  }
	  totEff2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex_Gen]->Fill(incrementTot, thetaAdjusted, phiFolded);
	}
	if(pTIndex_Gen > 0 && rapForPTIndex_Gen > 0){
	  if(useIndivEff){
	    recoEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementReco, thetaAdjusted, phiFolded);
	    // if(incrementReco)
	    trigEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementTrig, thetaAdjusted, phiFolded);
	  }
	  totEff2D_pol_pT_rap_phiFolded[iFrame][pTIndex_Gen][rapForPTIndex_Gen]->Fill(incrementTot, thetaAdjusted, phiFolded);
	} 
      }//frame
    }//useEventNTimes
  }//loop over entries

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}

//==============================================================
Double_t GetEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt, Bool_t use2DGraph){

  if(fabs(eta) > 2.4) return 0.; //no acceptance beyond 2.4

  if(iEff >= kNbEff)
    printf("<GetEfficiency> %d not a valid efficiency!!!\n", iEff);
  if(iEffSample >= kNbEffSample)
    printf("<GetEfficiency> %d not a valid efficiency sample!!!\n", iEffSample);

  Int_t binX, binY;
  Double_t eff;
  if(use2DGraph){
    binX = hMuEff_smoothed[iEff][iEffSample]->GetXaxis()->FindBin(fabs(eta));
    binY = hMuEff_smoothed[iEff][iEffSample]->GetYaxis()->FindBin(pt);
    eff = hMuEff_smoothed[iEff][iEffSample]->GetBinContent(binX, binY);
  }
  else{
    binX = hMuEff[iEff][iEffSample][CENTRAL]->GetXaxis()->FindBin(fabs(eta));
    binY = hMuEff[iEff][iEffSample][CENTRAL]->GetYaxis()->FindBin(pt);
    eff = hMuEff[iEff][iEffSample][CENTRAL]->GetBinContent(binX, binY);

    // if(iEffSample == MC && binX == 2 && binY == 14)
    //   eff = 0.96;
  }
  // printf("%s, efficiency for |eta|=%1.3f and pT=%1.2f GeV/c is %1.3f\n",
  // 	 effName[iEff], fabs(eta), pt, eff);

  //eff += 0.03; //upscale efficiency artificially
  //eff *= 1.03; //upscale efficiency artificially


  if(eff > 1.) eff = 1.;
  else if(eff < 0.) eff = 0.;

  return eff;
}

//==============================================================
Double_t GetTEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt){

  if(fabs(eta) > 2.4) return 0.; //no acceptance beyond 2.4

  if(iEff >= kNbEff)
    printf("<GetTEfficiency> %d not a valid efficiency!!!\n", iEff);
  if(iEffSample != MCTRUTH)
    printf("<GetTEfficiency> %d not a valid efficiency sample!!!\n", iEffSample);

  Int_t binX, binY;
  Double_t eff;
  Int_t globalBin = tMuEff[iEff][iEffSample]->FindFixBin(fabs(eta), pt);
  eff = tMuEff[iEff][iEffSample]->GetEfficiency(globalBin);

  if(eff > 1.) eff = 1.;
  else if(eff < 0.) eff = 0.;

  return eff;
}

//==============================================================
Double_t GetDimuEfficiency(Int_t iEffSample, Double_t cosTheta, Double_t phi){

  if(fabs(cosTheta) > 1) return 0.;
  if(phi < -180. || phi > 180.) return 0.;

  if(iEffSample >= kNbEffSample)
    printf("<GetDimuEfficiency> %d not a valid efficiency sample!!!\n", iEffSample);

  Int_t binX = hDiMuEff[iEffSample][CENTRAL]->GetXaxis()->FindBin(cosTheta);
  Int_t binY = hDiMuEff[iEffSample][CENTRAL]->GetYaxis()->FindBin(phi);
  Double_t eff = hDiMuEff[iEffSample][CENTRAL]->GetBinContent(binX, binY);
  
  if(eff > 1.) eff = 1.;
  else if(eff < 0.) eff = 0.;

  return eff;
}

//==============================================================
Double_t GetEfficiency_FromParametrization(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt){

  if(fabs(eta) > 2.4) return 0.; //no acceptance beyond 2.4

  if(iEff >= kNbEff)
    printf("<GetEfficiency_FromParametrization> %d not a valid efficiency!!!\n", iEff);
  if(iEffSample >= kNbEffSample)
    printf("<GetEfficiency_FromParametrization> %d not a valide efficiency sample!!!\n", iEffSample);

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
