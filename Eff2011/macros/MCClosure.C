#define MCClosure_cxx
#include "MCClosure.h"
#include "calcPol.C"
#include "isMuonInAcceptance.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>
#include "TRandom.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TGraph2D.h"

//load the rho factor
Char_t *effTypeName[3] = {"reco", "trig", "tot"};
enum {RHO_RECO, RHO_TRIG, RHO_TOT};
//Char_t *rhoFFileName = "rhoFactor_SingleMuEff_noDimuVtxEffCorr_noJpsiVprobCut_15Dec2011.root";
//Char_t *rhoFFileName = "rhoFactor_ProdSingleMuEff_noDimuVtxEffCorr_noJpsiVprobCut_15Dec2011.root";
Char_t *rhoFFileName = "rhoFactor_ProdSingleMuEff_16Dec2011.root";
Int_t const kNbMaxFrame = 3;
TH2D *hRho_pol[3][kNbMaxFrame][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];

//switch here between a total single muon efficiency:
Char_t *dimuEffName = "DimuVertexEff";
Char_t *dimuEffFileName = "/Users/hwoehri/CMS/Work/TnP/2011/Linlin/1Dec2011/DimuonVtxEff_Dimuon0Jpsi_cosTheta_Phi_TrkCuts80_CS_01Dec2011.root";
//
enum {SingleMuEff};
Bool_t useIndivEff = kFALSE;
Int_t const kNbEff = 1;
Char_t *effName[kNbEff] = {"SingleMuEff"};
Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/8Dec2011/EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Trk80Cuts_19Nov2011.root"};
//Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Linlin/7Dec2011/singleMuonEfficiency_ProbeTrackMatched_data_mc_pt_abseta_tracker80Cuts_7Dec2011.root"};//SingleMuEff a la Matt
//Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Linlin/2Dec2011/SingleMuEff_Dimuon10Jpsi_ProbeTrackMatched_data_mc_pt_abseta_tracker80Cuts_02Dec2011.root"};//SingleMuEff a la Matt (L1*L2*L3 trigger eff only)
//Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/7Dec2011/EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Trk80Cuts_19Nov2011.root"};
//Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/23Nov2011/EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Trk80Cuts_19Nov2011.root"};

Char_t *effFileName_Func[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/14Dec2011/fitProd_Trk80Cuts_14Dec2011.root"};
//Char_t *effFileName_Func[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/2011/Ilse/14Dec2011/fitOne_Trk80Cuts_14Dec2011.root"};


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
Int_t const kNbEffSample = 2;
enum {DATA, MC};
Char_t *effSampleName[kNbEffSample] = {"DATA", "MC"};
//
enum {CENTRAL, UPPER, LOWER};
TH2D *hMuEff[kNbEff][kNbEffSample][3];
TH2D *hMuEff_smoothed[kNbEff][kNbEffSample];
TH2D *hDiMuEff[kNbEffSample][3];
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

Int_t const kNbSteps = 6; //0... all generated; 1... generated+Acc; 2... GEN+ACC+RECO, 3... GEN+ACC+RECO+TRIG, 4... eff-corrected[GEN+ACC+RECO+TRIG], 5... eff-corrected[GEN+ACC+RECO+TRIG] using RECO variables
TH2D *hPTRap[kNbSteps];
TH1D *hCosTheta[eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1][kNbSteps];
TH1D *hPhi[eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1][kNbSteps];
Char_t *stepName[kNbSteps] = {"GEN", "GEN_ACC", "GEN_ACC_RECO", "GEN_ACC_RECO_TRIG", "effCorr", "effCorr_RECO"};
enum {GEN, GEN_ACC, GEN_ACC_RECO, GEN_ACC_RECO_TRIG, EFFCORR, EFFCORR_RECO};

Double_t GetDimuEfficiency(Int_t iEffSample, Double_t cosTheta, Double_t phi);
Double_t GetEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt, Bool_t use2DGraph = kFALSE);
Double_t GetEfficiency_FromParametrization(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt);
Double_t GetTotSingleMuEff(Bool_t usePTFit, Bool_t use2DGraph, Int_t effSample, Double_t etaMuPos, Double_t pTMuPos, Double_t etaMuNeg, Double_t pTMuNeg, Char_t *trigLabel);
Double_t GetRhoFactor(Int_t index, Int_t iFrame, Double_t cosTheta, Double_t phi, Int_t iRap, Int_t iPT);
//==============================================
void MCClosure::Loop(Int_t effSample, Char_t *trigLabel, Bool_t rejectCowboys, Bool_t useRhoFactor, Bool_t use2DGraph)
{
  if (fChain == 0) return;

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
    Double_t phiMuPos_Gen = muPos_Gen->Phi();
    Double_t phiMuNeg_Gen = muNeg_Gen->Phi();

    Double_t onia_Gen_mass = onia_Gen->M();
    Double_t onia_Gen_pt = onia_Gen->Pt();
    Double_t onia_Gen_P = onia_Gen->P();
    Double_t onia_Gen_eta = onia_Gen->PseudoRapidity();
    Double_t onia_Gen_rap = onia_Gen->Rapidity();
    Double_t onia_Gen_phi = onia_Gen->Phi();
    
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

    hPTRap[GEN]->Fill(fabs(onia_Gen_rap), onia_Gen_pt);
    calcPol(*muPos_Gen, *muNeg_Gen);
    for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
      hCosTheta[iFrame][rapForPTIndex_Gen][pTIndex_Gen][GEN]->Fill(thisCosTh[iFrame]);
      hPhi[iFrame][rapForPTIndex_Gen][pTIndex_Gen][GEN]->Fill(thisPhi[iFrame]);
    }

    if(rejectCowboys)
      if((phiMuNeg_Gen - phiMuPos_Gen) < 0.)
	continue;
    
    //if(!(isMuonInAcceptance(LOOSE, pTMuPos_Gen, etaMuPos_Gen) && isMuonInAcceptance(LOOSE, pTMuNeg_Gen, etaMuNeg_Gen)))
    if(!(isMuonInAcceptance(TIGHT, pTMuPos_Gen, etaMuPos_Gen) && isMuonInAcceptance(TIGHT, pTMuNeg_Gen, etaMuNeg_Gen)))
      continue;
    if(pTMuPos_Gen < 2.5 || pTMuNeg_Gen < 2.5)
      continue;
    if(fabs(etaMuPos_Gen) > 1.6 || fabs(etaMuNeg_Gen) > 1.6)
      continue;

    hPTRap[GEN_ACC]->Fill(fabs(onia_Gen_rap), onia_Gen_pt);
    for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
      hCosTheta[iFrame][rapForPTIndex_Gen][pTIndex_Gen][GEN_ACC]->Fill(thisCosTh[iFrame]);
      hPhi[iFrame][rapForPTIndex_Gen][pTIndex_Gen][GEN_ACC]->Fill(thisPhi[iFrame]);
    }

    //===============================
    //set up the trigger logic:
    Int_t trigValue = -5;
    if(strncmp("HLT_Dimuon10_Jpsi_Barrel", trigLabel, 24) == 0)
    //HLT_Dimuon10_Jpsi_Barrel_v3... 1.4E33
    //HLT_Dimuon10_Jpsi_Barrel_v5... no cowboys
    //HLT_Dimuon10_Jpsi_Barrel_v6... L1DoubleMu0_HighQ
    //HLT_Dimuon13_Jpsi_Barrel_v1... L1DoubleMu0_HighQ
      trigValue = HLT_Dimuon10_Jpsi_Barrel_v3; 
    else{
      printf("chosen trigger path, %s, not a valid option!\n", trigLabel);
      exit(0);
    }
    
    if(onia->Pt() < 990. && 
       fabs(onia->Rapidity()) < eff::rapMax && 
       JpsiVprob > 0.01 &&
       isMuonInAcceptance(TIGHT, muPos->Pt(), muPos->Eta()) &&
       isMuonInAcceptance(TIGHT, muNeg->Pt(), muNeg->Eta()) &&
       muPos->Pt() > 2.5 &&
       muNeg->Pt() > 2.5 &&
       fabs(muPos->Eta()) < 1.6 &&
       fabs(muNeg->Eta()) < 1.6){
      
      hPTRap[GEN_ACC_RECO]->Fill(fabs(onia_Gen_rap), onia_Gen_pt);
      for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
       	hCosTheta[iFrame][rapForPTIndex_Gen][pTIndex_Gen][GEN_ACC_RECO]->Fill(thisCosTh[iFrame]);
	hPhi[iFrame][rapForPTIndex_Gen][pTIndex_Gen][GEN_ACC_RECO]->Fill(thisPhi[iFrame]);
      }
    }
    else
      continue;
    
    if(trigValue == 1){
      hPTRap[GEN_ACC_RECO_TRIG]->Fill(fabs(onia_Gen_rap), onia_Gen_pt);
      for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
      	hCosTheta[iFrame][rapForPTIndex_Gen][pTIndex_Gen][GEN_ACC_RECO_TRIG]->Fill(thisCosTh[iFrame]);
       	hPhi[iFrame][rapForPTIndex_Gen][pTIndex_Gen][GEN_ACC_RECO_TRIG]->Fill(thisPhi[iFrame]);
      }
    }
    else
      continue;

    //calculate the efficiencies using GEN variables:
    //Bool_t usePTFit = kTRUE; //alternative to the T&P histograms use fitted pT differential efficiency
    Bool_t usePTFit = kFALSE; //alternative to the T&P histograms use fitted pT differential efficiency
    Double_t totEff = GetTotSingleMuEff(usePTFit, use2DGraph, effSample, etaMuPos_Gen, pTMuPos_Gen, etaMuNeg_Gen, pTMuNeg_Gen, trigLabel);
    //totEff *= GetDimuEfficiency(effSample, thisCosTh[eff::CS], thisPhi[eff::CS]);

    if(useRhoFactor && pTIndex_Gen > 5){
      Double_t rho[eff::kNbFrames];
      for(int iFrame = 0; iFrame < kNbMaxFrame; iFrame++){
	rho[iFrame] = GetRhoFactor(RHO_TOT, iFrame, thisCosTh[iFrame], thisPhi[iFrame], rapForPTIndex_Gen, pTIndex_Gen);

	if(rho[iFrame] > 1e-5){
	  hCosTheta[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR]->Fill(thisCosTh[iFrame], 1./(totEff*rho[iFrame]));
	  hPhi[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR]->Fill(thisPhi[iFrame], 1./(totEff*rho[iFrame]));
	}
	else{
	  hCosTheta[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR]->Fill(thisCosTh[iFrame], 1./(totEff*1e-5));
	  hPhi[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR]->Fill(thisPhi[iFrame], 1./(totEff*1e-5));
	}
      }
      if(rho[eff::CS] > 1e-5)
	hPTRap[EFFCORR]->Fill(fabs(onia_Gen_rap), onia_Gen_pt, 1./(totEff*rho[eff::CS]));
      else
	hPTRap[EFFCORR]->Fill(fabs(onia_Gen_rap), onia_Gen_pt, 1./(totEff*1e-5));
    }
    else{
      hPTRap[EFFCORR]->Fill(fabs(onia_Gen_rap), onia_Gen_pt, 1./totEff);
      for(int iFrame = 0; iFrame < kNbMaxFrame; iFrame++){
	hCosTheta[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR]->Fill(thisCosTh[iFrame], 1./totEff);
	hPhi[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR]->Fill(thisPhi[iFrame], 1./totEff);
      }
    }
    //calculate the efficiencies using the RECO variables:
    totEff = GetTotSingleMuEff(usePTFit, use2DGraph, effSample, muPos->Eta(), muPos->Pt(), muNeg->Eta(), muNeg->Pt(), trigLabel);
    calcPol(*muPos, *muNeg);
    //totEff *= GetDimuEfficiency(effSample, thisCosTh[eff::CS], thisPhi[eff::CS]);
    if(useRhoFactor && pTIndex_Gen > 5){
      Double_t rho[eff::kNbFrames];
      for(int iFrame = 0; iFrame < kNbMaxFrame; iFrame++){
	rho[iFrame] = GetRhoFactor(RHO_TOT, iFrame, thisCosTh[iFrame], thisPhi[iFrame], rapForPTIndex_Gen, pTIndex_Gen);

	if(rho[iFrame] > 1e-5){
	  hCosTheta[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR_RECO]->Fill(thisCosTh[iFrame], 1./(totEff*rho[iFrame]));
	  hPhi[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR_RECO]->Fill(thisPhi[iFrame], 1./(totEff*rho[iFrame]));
	}
	else{
	  hCosTheta[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR_RECO]->Fill(thisCosTh[iFrame], 1./(totEff*1e-5));
	  hPhi[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR_RECO]->Fill(thisPhi[iFrame], 1./(totEff*1e-5));
	}
      }
      if(rho[eff::CS] > 1e-5)
	hPTRap[EFFCORR_RECO]->Fill(fabs(onia_Gen_rap), onia_Gen_pt, 1./(totEff*rho[eff::CS]));
      else
	hPTRap[EFFCORR_RECO]->Fill(fabs(onia_Gen_rap), onia_Gen_pt, 1./(totEff*1e-5));
    }
    else{
      hPTRap[EFFCORR_RECO]->Fill(fabs(onia_Gen_rap), onia_Gen_pt, 1./totEff);
      for(int iFrame = 0; iFrame < kNbMaxFrame; iFrame++){
	hCosTheta[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR_RECO]->Fill(thisCosTh[iFrame], 1./totEff);
	hPhi[iFrame][rapForPTIndex_Gen][pTIndex_Gen][EFFCORR_RECO]->Fill(thisPhi[iFrame], 1./totEff);
      }
    }
  }//loop over entries

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}

//==============================================================
Double_t GetEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt, Bool_t use2DGraph){

  if(fabs(eta) > 2.4) return 0.; //no acceptance beyond 2.4

  if(iEff >= kNbEff)
    printf("%d not a valid efficiency!!!\n", iEff);
  if(iEffSample >= kNbEffSample)
    printf("%d not a valid efficiency sample!!!\n", iEffSample);

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
Double_t GetDimuEfficiency(Int_t iEffSample, Double_t cosTheta, Double_t phi){

  if(fabs(cosTheta) > 1) return 0.;
  if(phi < -180. || phi > 180.) return 0.;

  if(iEffSample >= kNbEffSample)
    printf("%d not a valid efficiency sample!!!\n", iEffSample);

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

//============================================
Double_t GetTotSingleMuEff(Bool_t usePTFit, Bool_t use2DGraph,
			   Int_t effSample, 
			   Double_t etaMuPos, Double_t pTMuPos, 
			   Double_t etaMuNeg, Double_t pTMuNeg,
			   Char_t *trigLabel){

  Double_t totEff = 0.;
  if(!useIndivEff){

    totEff = 0.99*0.99; //tracking efficiency
    //      totEff = 1.0 * 1.0;

    if(usePTFit){	
      totEff *= GetEfficiency_FromParametrization(SingleMuEff, effSample, etaMuPos, pTMuPos);
      totEff *= GetEfficiency_FromParametrization(SingleMuEff, effSample, etaMuNeg, pTMuNeg);
    }
    else{
      totEff *= GetEfficiency(SingleMuEff, effSample, etaMuPos, pTMuPos, use2DGraph);
      totEff *= GetEfficiency(SingleMuEff, effSample, etaMuNeg, pTMuNeg, use2DGraph);
    }
    // //dimuon vertexing cut NOT included in the (single) muon quality cuts...
    // if(JpsiVprob < 0.01)
    // 	totEff = 0.;
  }
  else{
    Double_t epsTrack_Pos  = 0.99;
    Double_t epsTrack_Neg = 0.99;
    // Double_t epsTrack_Pos  = 1.0;
    // Double_t epsTrack_Neg = 1.0;
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
    
    //dimuon vertexing cut NOT included in the (single) muon quality cuts...
    // if(JpsiVprob < 0.01)
    // 	recoEff = 0.;
    
    Double_t epsL1L2Trig_Pos, epsL3Trig_Pos, epsL1L2Trig_Neg, epsL3Trig_Neg;
    Double_t trigEff = 0.;
    
    //if(incrementReco){

      if(strncmp("HLT_Dimuon10_Jpsi_Barrel", trigLabel, 24) == 0){
	//HLT_Dimuon10_Jpsi_Barrel_v3... 1.4E33
	//HLT_Dimuon10_Jpsi_Barrel_v5... no cowboys
	//HLT_Dimuon10_Jpsi_Barrel_v6... L1DoubleMu0_HighQ
	//HLT_Dimuon13_Jpsi_Barrel_v1... L1DoubleMu0_HighQ
	
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
      //    }
    totEff = trigEff * recoEff;
  }

  return totEff;
}

//=======================================
Double_t GetRhoFactor(Int_t index, Int_t iFrame, Double_t cosTheta, Double_t phi, Int_t iRap, Int_t iPT){

  // printf("index %d, frame %d, cosTheta %1.3f, phi %1.3f, rap %d, pT %d\n",
  // 	 index, iFrame, cosTheta, phi, iRap, iPT);

  if(fabs(cosTheta) > 1.) return 0.;
  if(phi < -180. || phi > 180.) return 0.;

  Int_t binX = hRho_pol[index][iFrame][iRap][iPT]->GetXaxis()->FindBin(cosTheta);
  Int_t binY = hRho_pol[index][iFrame][iRap][iPT]->GetYaxis()->FindBin(phi);
  Double_t rho = hRho_pol[index][iFrame][iRap][iPT]->GetBinContent(binX, binY);

  return rho;
}
