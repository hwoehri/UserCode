#define DataTnPTriggerEff_cxx
#include "DataTnPTriggerEff.h"
#include "calcPol.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>
#include "TRandom.h"
#include "TF1.h"


Int_t const kNbEff = 4;
Char_t *effFileNames[kNbEff] = {//"/Users/hwoehri/CMS/Work/TnP/Francesco/2April2011/L1L2_DMu0_TriggerEfficiencies_RUNA_distM2gt120_7April2011_fitted.root",//RunA
                                "/Users/hwoehri/CMS/Work/TnP/Francesco/2April2011/L1L2_DMu0_TriggerEfficiencies_RUNB_distM2gt120_7April2011_fitted.root",//RunB
				//"/Users/hwoehri/CMS/Work/TnP/Luigi/5April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010A_distM1gt150_7April2011_fitted.root",//RunA
				"/Users/hwoehri/CMS/Work/TnP/Luigi/5April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM1gt150_7April2011_fitted.root",//RunB
				"/Users/hwoehri/CMS/Work/TnP/Herbert/28March2011/HLT_Track_Mu0TkMu0_Run1.root",
				//change input file for different instances of the low pT J/psi trigger
				"/Users/hwoehri/CMS/Work/TnP/Ilse/24March2011/HLTMuonTrack_Mu0_TkMu0_TM.root"};
enum {L1L2Eff, L3Eff, Trk_TkMu0, Mu_TkMu0};
Char_t *effName[kNbEff] = {"L1L2Eff", "L3Eff", "Trk_TkMu0Eff", "Mu_TkMu0Eff"};
Int_t const kNbEffSample = 3;
enum {DATA, MC, MCTRUTH};
Char_t *effSampleName[kNbEffSample] = {"DATA", "MC", "MCTRUTH"};
//
enum {CENTRAL, UPPER, LOWER};
TH2D *hMuEff[kNbEff][kNbEffSample][3];
Int_t const kNbEtaBins = 8;
Float_t binsEta[kNbEtaBins] = {0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4}; //needs manual adjustment!!!
TF1 *fMuEff_pT[kNbEff][kNbEffSample][kNbEtaBins];

Int_t const kNbPTMaxBins = 4;
Int_t const kNbRapForPTBins = 2;
Int_t const kNbPTBins[kNbRapForPTBins+1] = {4,4,4};//all y, y1, y2
Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {{0., 10., 20., 30., 500.},//all y
						       {0., 10., 20., 30., 500.},//y1
						       {0., 10., 20., 30., 500.}};//y2
Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 1.2, 2.4};

//=================================================
//1D histos versus dimuon mass
//=================================================
TH1D *reco_M[kNbPTMaxBins+1][kNbRapForPTBins+1];
TH1D *passedTrigger_M[kNbPTMaxBins+1][kNbRapForPTBins+1];
TH1D *failedTrigger_M[kNbPTMaxBins+1][kNbRapForPTBins+1];

Double_t GetEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt);
Double_t GetEfficiency_FromParametrization(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt);
//==============================================
void DataTnPTriggerEff::Loop(Int_t effSample, Char_t *trigLabel)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nb = 0;
  Long64_t countRecEvent = 0;
  printf("number of entries = %d\n", (Int_t) nentries);

  Bool_t incrementTrig, incrementTot;

  //loop over the events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    // for (Long64_t jentry=0; jentry<1000000;jentry++) {

    if(jentry % 100000 == 0) printf("event %d\n", (Int_t) jentry);
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    //make sure we have a reconstructed event:
    if(onia->Pt() > 990.) continue;

    //take only events within the following mass window 2.6 < M < 3.4 GeV
    if(onia->M() < 2.6 || onia->M() > 3.4) continue;

    //===================================================================
    //1.) check all reconstructed (and filtered) events and 
    //fill the corresponding histograms
    //===================================================================

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

    if(JpsiVprob < 0.01)
      continue;

    Double_t onia_mass = onia->M();
    Double_t onia_pt = onia->Pt();
    Double_t onia_P = onia->P();
    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();
    Double_t onia_phi = onia->Phi();
    Double_t onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);

    if(fabs(onia_rap) > eff::rapMax)
      continue;
    
    Int_t rapForPTIndex = -1;
    for(int iRap = 0; iRap < kNbRapForPTBins; iRap++){
      if(TMath::Abs(onia_rap) > rapForPTRange[iRap] && 
	 TMath::Abs(onia_rap) < rapForPTRange[iRap+1]){
	rapForPTIndex = iRap+1;
	break;
      }
    }
    Int_t pTIndex = -1;
    for(int iPT = 0; iPT < kNbPTBins[rapForPTIndex]; iPT++){
      if(onia_pt > pTRange[rapForPTIndex][iPT] && onia_pt < pTRange[rapForPTIndex][iPT+1]){
	pTIndex = iPT+1;
	break;
      }
    }
    Int_t rapIntegratedPTIndex = -1;
    for(int iPT = 0; iPT < kNbPTBins[0]; iPT++){
      if(onia_pt > pTRange[0][iPT] && onia_pt < pTRange[0][iPT+1]){
	rapIntegratedPTIndex = iPT+1;
	break;
      }
    }
    if(rapForPTIndex < 1 || pTIndex < 1 || rapIntegratedPTIndex < 1){
      printf("onia y = %1.3f, pT = %1.3f; rapForPTIndex %d, pTIndex %d, rapIntegratedPTIndex %d\n",
	     onia_rap, onia_pt, rapForPTIndex, pTIndex, rapIntegratedPTIndex);
    }
    // //==============================
    // calcPol(*muPos, *muNeg);
    // //==============================
     // ---> For the asymmetric triggers:
      // 0 : event not firing the corresponding trigger
      // 1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, POSITIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
      // -1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, NEGATIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
      // 2 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, both matched to the tighter HLT object (usually a L3 muon)
      // 3 : event firing the corresponding trigger but at least one RECO muon is not matched to the HLT objects 

    Double_t epsL1L2Trig_Pos, epsL3Trig_Pos, epsL1L2Trig_Neg, epsL3Trig_Neg;
    if(strncmp("HLT_DoubleMu0", trigLabel, 13) == 0){
      epsL1L2Trig_Pos = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuPos, pTMuPos);
      epsL3Trig_Pos   = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuPos, pTMuPos);
      epsL1L2Trig_Neg = GetEfficiency_FromParametrization(L1L2Eff, effSample, etaMuNeg, pTMuNeg);
      epsL3Trig_Neg   = GetEfficiency_FromParametrization(L3Eff, effSample, etaMuNeg, pTMuNeg);
    }
    else if(strncmp("HLT_Mu0_TkMu0_OST_Jpsi", trigLabel, 22) == 0){ //
      if(HLT_Mu0_TkMu0_OST_Jpsi == 1){//pos. muon matched to L3 muon, neg. muon to tk-muon
	epsL1L2Trig_Pos  = GetEfficiency(L1L2Eff, effSample, etaMuPos, pTMuPos);
	epsL3Trig_Pos    = GetEfficiency(L3Eff, effSample, etaMuPos, pTMuPos);
	epsL1L2Trig_Neg  = GetEfficiency(Mu_TkMu0, effSample, etaMuNeg, pTMuNeg);
	epsL3Trig_Neg    = GetEfficiency(Trk_TkMu0, effSample, etaMuNeg, pTMuNeg);
      }
      else if(HLT_Mu0_TkMu0_OST_Jpsi == -1){//neg. muon matched to L3 muon, pos. muon to tk-muon
	epsL1L2Trig_Neg  = GetEfficiency(L1L2Eff, effSample, etaMuNeg, pTMuNeg);
	epsL3Trig_Neg    = GetEfficiency(L3Eff, effSample, etaMuNeg, pTMuNeg);
	epsL1L2Trig_Pos  = GetEfficiency(Mu_TkMu0, effSample, etaMuPos, pTMuPos);
	epsL3Trig_Pos    = GetEfficiency(Trk_TkMu0, effSample, etaMuPos, pTMuPos);
      }
      else if(HLT_Mu0_TkMu0_OST_Jpsi == 2){//both matched to L3 muon, as well as tk-muon
	//assign "randomly" the efficiency of the L3 muon 
	//either to the pos or to the neg muon
	if(gRandom->Uniform() > 0.5){
	  epsL1L2Trig_Pos  = GetEfficiency(Mu_TkMu0, effSample, etaMuPos, pTMuPos);
	  epsL3Trig_Pos    = GetEfficiency(Trk_TkMu0, effSample, etaMuPos, pTMuPos);

	  epsL1L2Trig_Neg  = GetEfficiency(L1L2Eff, effSample, etaMuNeg, pTMuNeg);
	  epsL3Trig_Neg    = GetEfficiency(L3Eff, effSample, etaMuNeg, pTMuNeg);
	}
	else{
	  epsL1L2Trig_Pos  = GetEfficiency(L1L2Eff, effSample, etaMuPos, pTMuPos);
	  epsL3Trig_Pos    = GetEfficiency(L3Eff, effSample, etaMuPos, pTMuPos);

	  epsL1L2Trig_Neg  = GetEfficiency(Mu_TkMu0, effSample, etaMuNeg, pTMuNeg);
	  epsL3Trig_Neg    = GetEfficiency(Trk_TkMu0, effSample, etaMuNeg, pTMuNeg);
	}
      }
    }

    Double_t trigEff = epsL1L2Trig_Pos * epsL3Trig_Pos * epsL1L2Trig_Neg * epsL3Trig_Neg;
    // if(trigEff < 0.1){
    //   printf("muPos: pT %1.3f, eta %1.3f, muNeg: pT %1.3f, eta %1.3f; effPos: L1 %1.3f, L3 %1.3f, effNeg: L1 %1.3f, L3 %1.3f --> trigEff %1.3f\n",
    // 	     pTMuPos, etaMuPos, pTMuNeg, etaMuNeg, epsL1L2Trig_Pos, epsL3Trig_Pos, epsL1L2Trig_Neg, epsL3Trig_Neg, trigEff);
    // }

    //the indiv. histograms will be filled depending on the
    //assigned probability (given by the trigger efficiency)
    //--> this way the histograms have Possonian entries
    //and the option "B" when dividing them can be used
    Double_t randNb = gRandom->Uniform();
    if(trigEff > randNb) incrementTrig = kTRUE;
    else incrementTrig = kFALSE;

    reco_M[0][0]->Fill(onia_mass);
    if(incrementTrig)//passed trigger condition
      passedTrigger_M[0][0]->Fill(onia_mass);
    else//failed trigger condition
      failedTrigger_M[0][0]->Fill(onia_mass);


    if(rapForPTIndex < 1){
      // printf("rapForPTIndex %d, rap(onia) = %f\n", rapForPTIndex, onia_rap);
      continue;
    }
    if(pTIndex < 1){
      // printf("pTIndex %d, pT(onia) = %f\n", pTIndex, onia_pt);
      continue;
    }




    if(rapIntegratedPTIndex > 0){
      reco_M[rapIntegratedPTIndex][0]->Fill(onia_mass);
      if(incrementTrig)//passed trigger condition
	passedTrigger_M[rapIntegratedPTIndex][0]->Fill(onia_mass);
      else//failed trigger condition
	failedTrigger_M[rapIntegratedPTIndex][0]->Fill(onia_mass);
    }
    if(rapForPTIndex > 0){
      reco_M[0][rapForPTIndex]->Fill(onia_mass);
      if(incrementTrig)//passed trigger condition
	passedTrigger_M[0][rapForPTIndex]->Fill(onia_mass);
      else//failed trigger condition
	failedTrigger_M[0][rapForPTIndex]->Fill(onia_mass);
    }
    // printf("pTIndex %f, rapForPTIndex %f\n", pTIndex, rapForPTIndex);
    if(pTIndex > 0 && rapForPTIndex > 0){
      reco_M[pTIndex][rapForPTIndex]->Fill(onia_mass);
      if(incrementTrig)//passed trigger condition
	passedTrigger_M[pTIndex][rapForPTIndex]->Fill(onia_mass);
      else//failed trigger condition
	failedTrigger_M[pTIndex][rapForPTIndex]->Fill(onia_mass);
    }
  }//loop over entries

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}

//==============================================================
Double_t GetEfficiency(Int_t iEff, Int_t iEffSample, Double_t eta, Double_t pt){

  if(fabs(eta) > 2.4) return 0.; //no acceptance beyond 2.4
//   if(pt > 20.) return 1.; //efficiencies beyond 20 GeV/c not reliable from T&P, but are always 100%

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
