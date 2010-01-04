#define ProjectQQ_cxx
#include "ProjectQQ.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"

#include "commonVar.inc"

//dimuon related histograms
TH1F *hReco_mass[kNbCharge][kNbCath][kNbCuts];
TH1F *hReco_pT[kNbDimuonSet][kNbCharge][kNbCath][kNbCuts];
TH1F *hReco_rap[kNbDimuonSet][kNbCharge][kNbCath][kNbCuts];
TH1F *hRecoBest_mass[kNbCharge][kNbCath];
TH1F *hRecoBest_pT[kNbDimuonSet][kNbCharge][kNbCath];
TH1F *hRecoBest_rap[kNbDimuonSet][kNbCharge][kNbCath];
TH1F *hMult_QQ[kNbCharge][kNbCath][kNbCuts];
TH2F *hReco_eta1_eta2[kNbDimuonSet][kNbCharge][kNbCath+1][kNbCuts];
TH1F *hRecoNoBest;

//single muon histos for two dimuon mass windows:
TH1F *hMuon_eta[kNbDimuonSet][kNbCharge][kNbCath];
TH1F *hMuon_pT[kNbDimuonSet][kNbCharge][kNbCath];
TH2F *hMuon_pT_eta[kNbDimuonSet][kNbCharge][kNbMuCath];
// quality histogams
TH1F *hD0[kNbDimuonSet][kNbCharge][kNbMuCath];
TH1F *hDz[kNbDimuonSet][kNbCharge][kNbMuCath];
TH1F *hChi2GlobalFit[kNbDimuonSet][kNbCharge];
TH1F *hChi2TrackerFit[kNbDimuonSet][kNbCharge][kNbMuCath];
TH1F *hTM2DCompatibilityTight[kNbDimuonSet][kNbCharge];
TH1F *hTMLastStationOptimizedLowPtLoose[kNbDimuonSet][kNbCharge];
TH1F *hTrackerMuonArbitrated[kNbDimuonSet][kNbCharge];//tracker only
TH1F *hCaloCompatibility[kNbDimuonSet][kNbCharge][kNbMuCath];
TH1F *hNHitsSilicon[kNbDimuonSet][kNbCharge][kNbMuCath];
TH1F *hDimuVtxProb[kNbDimuonSet][kNbCharge][kNbCath];
//=======================================
void ProjectQQ::Loop(Bool_t removeQQ, Bool_t matchMC, Bool_t printGoodEvents)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  //    nentries = 10000;

  Long64_t nb = 0;

  Int_t countAllRecoQQEvents = 0, countNotFoundBestQQEvents = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
 
    if(jentry % 1000 == 0) printf("event %d/%d\n", (Int_t) jentry, (Int_t) nentries);

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    if(removeQQ && Mc_QQ_size > 0) continue; //Mc_QQ_size contains the nb. of MC J/psi's

    Int_t chargeID;
    Int_t nComb[kNbCharge][kNbCath][kNbCuts];
    for(int iCharge = 0; iCharge < kNbCharge; iCharge++)
      for(int iCat = 0; iCat < kNbCath; iCat++)
	for(int iCut = 0; iCut < kNbCuts; iCut++)
	  nComb[iCharge][iCat][iCut] = 0;

    //loop over all quarkonia
    for(int iQQ = 0; iQQ < Reco_QQ_size; iQQ++){

      chargeID = 0;
      if(Reco_QQ_sign[iQQ] < 0) chargeID = 1;
      else if(Reco_QQ_sign[iQQ] > 0) chargeID = 2;

      TLorentzVector *Reco_QQ = (TLorentzVector *) Reco_QQ_4mom->At(iQQ);
      Double_t mass_Reco_QQ = Reco_QQ->M();
      Double_t pT_Reco_QQ = Reco_QQ->Pt();
      Double_t rap_Reco_QQ = Reco_QQ->Rapidity();

      //get the individual muons:
      TLorentzVector *Reco_QQ_posMu, *Reco_QQ_negMu;
      Bool_t etaCorr = kFALSE;
      if(Reco_QQ_type[iQQ] == 0){ //global-global
	if(Reco_QQ_mupl[iQQ] < Reco_mu_glb_size && Reco_QQ_mumi[iQQ] < Reco_mu_glb_size){
	  Reco_QQ_posMu = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_mupl[iQQ]);
	  Reco_QQ_negMu = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_mumi[iQQ]);
	  etaCorr = kTRUE;
	}
      }
      else if(Reco_QQ_type[iQQ] == 2){ //tracker-tracker
	if(Reco_QQ_mupl[iQQ] < Reco_mu_trk_size && Reco_QQ_mumi[iQQ] < Reco_mu_trk_size){
	  Reco_QQ_posMu = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_mupl[iQQ]);
	  Reco_QQ_negMu = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_mumi[iQQ]);
	  etaCorr = kTRUE;
	}
      }

      Bool_t isMatched[2];//only relevant for MC 
      for (int iCut=0; iCut<kNbCuts; iCut++)
	{
	  //start with the MC matching:
// 	  if (iCut == 0 && matchMC){
 	  if (matchMC){

	    if(Mc_QQmupl_indx[0] > 1 || Mc_QQmumi_indx[0] > 1) continue; //sanity check

	    isMatched[0] = kFALSE;
	    isMatched[1] = kFALSE;

	    switch (Reco_QQ_type[iQQ]) {
	    case 0: // global-global
	      //first MC muon
	      if(Reco_QQ_muhpt[iQQ] < Reco_mu_glb_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iQQ] < Reco_mu_glb_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    case 1: // global-tracker
	      //first MC muon
	      if(Reco_QQ_muhpt[iQQ] < Reco_mu_glb_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iQQ] < Reco_mu_trk_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    case 2: // tracker-tracker
	      //first MC muon
	      if(Reco_QQ_muhpt[iQQ] < Reco_mu_trk_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iQQ] < Reco_mu_trk_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    case 3: // global-calo
	      //first MC muon
	      if(Reco_QQ_muhpt[iQQ] < Reco_mu_glb_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iQQ] < Reco_mu_cal_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    case 4: // tracker-calo
	      //first MC muon
	      if(Reco_QQ_muhpt[iQQ] < Reco_mu_trk_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iQQ] < Reco_mu_cal_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    case 5: // calo-calo
	      //first MC muon
	      if(Reco_QQ_muhpt[iQQ] < Reco_mu_cal_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_muhpt[iQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iQQ] < Reco_mu_cal_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    default:
	      break;
	    }	    
	    //check the MC matching before continuing...
	    if(isMatched[0] == kFALSE || isMatched[1] == kFALSE)
	      continue;
	  }

	  if (iCut>0) {
	    // distance cut
	    switch (Reco_QQ_type[iQQ]) {
	    case 0: // global-global
	      if ( !(fabs(Reco_mu_glb_d0[Reco_QQ_muhpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_glb_dz[Reco_QQ_muhpt[iQQ]]) < MAX_dz_trk &&
		     fabs(Reco_mu_glb_d0[Reco_QQ_mulpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_glb_dz[Reco_QQ_mulpt[iQQ]]) < MAX_dz_trk) )
		continue;
	      break;
	    case 1: // global-tracker
	      if ( !(fabs(Reco_mu_glb_d0[Reco_QQ_muhpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_glb_dz[Reco_QQ_muhpt[iQQ]]) < MAX_dz_trk &&
		     fabs(Reco_mu_trk_d0[Reco_QQ_mulpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_trk_dz[Reco_QQ_mulpt[iQQ]]) < MAX_dz_trk) )
		continue;
	      break;
	    case 2: // tracker-tracker
	      if ( !(fabs(Reco_mu_trk_d0[Reco_QQ_muhpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_trk_dz[Reco_QQ_muhpt[iQQ]]) < MAX_dz_trk &&
		     fabs(Reco_mu_trk_d0[Reco_QQ_mulpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_trk_dz[Reco_QQ_mulpt[iQQ]]) < MAX_dz_trk) )
		continue;
	      break;
	    case 3: // global-calo
	      if ( !(fabs(Reco_mu_glb_d0[Reco_QQ_muhpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_glb_dz[Reco_QQ_muhpt[iQQ]]) < MAX_dz_trk &&
		     fabs(Reco_mu_cal_d0[Reco_QQ_mulpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_cal_dz[Reco_QQ_mulpt[iQQ]]) < MAX_dz_trk) )
		continue;
	      break;
	    case 4: // tracker-calo
	      if ( !(fabs(Reco_mu_trk_d0[Reco_QQ_muhpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_trk_dz[Reco_QQ_muhpt[iQQ]]) < MAX_dz_trk &&
		     fabs(Reco_mu_cal_d0[Reco_QQ_mulpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_cal_dz[Reco_QQ_mulpt[iQQ]]) < MAX_dz_trk) )
		continue;
	      break;
	    case 5: // calo-calo
	      if ( !(fabs(Reco_mu_cal_d0[Reco_QQ_muhpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_cal_dz[Reco_QQ_muhpt[iQQ]]) < MAX_dz_trk &&
		     fabs(Reco_mu_cal_d0[Reco_QQ_mulpt[iQQ]]) < MAX_d0_trk &&
		     fabs(Reco_mu_cal_dz[Reco_QQ_mulpt[iQQ]]) < MAX_dz_trk) )
		continue;
	      break;
	    default:
	      break;
	    }
	  }
	    
	  if (iCut>1) {
	    // for tracker muons: TM2DCompatibilityTight || TMLastStationOptimizedLowPtTight 
	    // for calo muons: calo compatibility
	    switch (Reco_QQ_type[iQQ]) {
	    case 0: // global-global
	      break;
	    case 1: // global-tracker
	      if ((Reco_mu_trk_PIDmask[Reco_QQ_mulpt[iQQ]] & (int) pow(2,5))/(int) pow(2,5) == 0 &&
		  (Reco_mu_trk_PIDmask[Reco_QQ_mulpt[iQQ]] & (int) pow(2,8))/(int) pow(2,8) == 0)
		continue;
	      break;
	    case 2: // tracker-tracker
	      if ( ((Reco_mu_trk_PIDmask[Reco_QQ_muhpt[iQQ]] & (int) pow(2,5))/(int) pow(2,5) == 0 &&
		    (Reco_mu_trk_PIDmask[Reco_QQ_muhpt[iQQ]] & (int) pow(2,8))/(int) pow(2,8) == 0) ||
		   ((Reco_mu_trk_PIDmask[Reco_QQ_mulpt[iQQ]] & (int) pow(2,5))/(int) pow(2,5) == 0 &&
		    (Reco_mu_trk_PIDmask[Reco_QQ_mulpt[iQQ]] & (int) pow(2,8))/(int) pow(2,8) == 0) )
		continue;
	      break;
	    case 3: // global-calo
	      if ( !(Reco_mu_cal_caloComp[Reco_QQ_mulpt[iQQ]] > MIN_caloComp) )
		continue;
	      break;
	    case 4: // tracker-calo
	      if ( ((Reco_mu_trk_PIDmask[Reco_QQ_muhpt[iQQ]] & (int) pow(2,5))/(int) pow(2,5) == 0 &&
		    (Reco_mu_trk_PIDmask[Reco_QQ_muhpt[iQQ]] & (int) pow(2,8))/(int) pow(2,8) == 0) ||
		   Reco_mu_cal_caloComp[Reco_QQ_mulpt[iQQ]] < MIN_caloComp)
		continue;
	      break;
	    case 5: // calo-calo
	      if ( !(Reco_mu_cal_caloComp[Reco_QQ_muhpt[iQQ]] > MIN_caloComp &&
		     Reco_mu_cal_caloComp[Reco_QQ_mulpt[iQQ]] > MIN_caloComp) )
		continue;
	      break;
	    default:
	      break;
	    }
	  }

	  if (iCut>2) {
	    // for global: chi2/ndf (global fit)
	    // for tracker and calo: chi2/ndf (track fit) && #hits
	    switch (Reco_QQ_type[iQQ]) {
	    case 0: // global-global
	      if ( !(Reco_mu_glb_normChi2[Reco_QQ_muhpt[iQQ]] < MAX_normchi2_glb &&
		     Reco_mu_glb_normChi2[Reco_QQ_mulpt[iQQ]] < MAX_normchi2_glb) )
		continue;
	      break;
	    case 1: // global-tracker
	      if ( !(Reco_mu_glb_normChi2[Reco_QQ_muhpt[iQQ]] < MAX_normchi2_glb &&
		     Reco_mu_trk_normChi2[Reco_QQ_mulpt[iQQ]] < MAX_normchi2_trk &&
		     Reco_mu_trk_nhitstrack[Reco_QQ_mulpt[iQQ]] > MIN_nhits_trk) )
		continue;
	      break;
	    case 2: // tracker-tracker
	      if ( !(Reco_mu_trk_normChi2[Reco_QQ_muhpt[iQQ]] < MAX_normchi2_trk &&
		     Reco_mu_trk_nhitstrack[Reco_QQ_muhpt[iQQ]] > MIN_nhits_trk &&		       
		     Reco_mu_trk_normChi2[Reco_QQ_mulpt[iQQ]] < MAX_normchi2_trk &&
		     Reco_mu_trk_nhitstrack[Reco_QQ_mulpt[iQQ]] > MIN_nhits_trk) )
		continue;
	      break;
	    case 3: // global-calo
	      if ( !(Reco_mu_glb_normChi2[Reco_QQ_muhpt[iQQ]] < MAX_normchi2_glb &&
		     Reco_mu_cal_normChi2[Reco_QQ_mulpt[iQQ]] < MAX_normchi2_cal &&
		     Reco_mu_cal_nhitstrack[Reco_QQ_mulpt[iQQ]] > MIN_nhits_trk) )
		continue;
	      break;
	    case 4: // tracker-calo
	      if ( !(Reco_mu_trk_normChi2[Reco_QQ_muhpt[iQQ]] < MAX_normchi2_trk &&
		     Reco_mu_trk_nhitstrack[Reco_QQ_muhpt[iQQ]] > MIN_nhits_trk &&
		     Reco_mu_cal_normChi2[Reco_QQ_mulpt[iQQ]] < MAX_normchi2_cal &&
		     Reco_mu_cal_nhitstrack[Reco_QQ_mulpt[iQQ]] > MIN_nhits_trk) )
		continue;
	      break;
	    case 5: // calo-calo
	      if ( !(Reco_mu_cal_normChi2[Reco_QQ_muhpt[iQQ]] < MAX_normchi2_cal &&
		     Reco_mu_cal_nhitstrack[Reco_QQ_muhpt[iQQ]] > MIN_nhits_trk &&
		     Reco_mu_cal_normChi2[Reco_QQ_mulpt[iQQ]] < MAX_normchi2_cal &&
		     Reco_mu_cal_nhitstrack[Reco_QQ_mulpt[iQQ]] > MIN_nhits_trk) )
		continue;
	      break;
	    default:
	      break;
	    }
	  }
	    
	  if (iCut>3) {
	    // dimuon vertex probability
	    if (Reco_QQ_probChi2[iQQ]<= MIN_vtxprob_jpsi)
	      continue;
	    else
	      if(printGoodEvents && Reco_QQ_type[iQQ] < 3 && Reco_QQ_sign[iQQ] == 0){
		printf("run %d, lumiSec %d event %d contains a dimuon of type %d with mass %1.2f GeV, pT %1.2f GeV, prob(vertex-chi2) = %1.3e and a cTau of %1.3e\n", 
		       runNb, lumiBlock, eventNb, Reco_QQ_type[iQQ], mass_Reco_QQ, pT_Reco_QQ, Reco_QQ_probChi2[iQQ], Reco_QQ_ctau[iQQ]);
// 		printf("contains a dimuon of type %d with mass %1.2f, pT %1.2f, prob(vertex-chi2) = %1.3e and a cTau of %1.3e\n", 
// 		       Reco_QQ_type[iQQ], mass_Reco_QQ, pT_Reco_QQ, Reco_QQ_probChi2[iQQ], Reco_QQ_ctau[iQQ]);
	      }
	  }

	  nComb[chargeID][Reco_QQ_type[iQQ]][iCut]++;

	  //mass histograms
	  hReco_mass[chargeID][Reco_QQ_type[iQQ]][iCut]->Fill(mass_Reco_QQ);
	  if(chargeID == 1 || chargeID == 2) //LS histo
	    hReco_mass[3][Reco_QQ_type[iQQ]][iCut]->Fill(mass_Reco_QQ);
	    
	  //fill pT and rap histos only for the signal:
	  for(int iMassW = 0; iMassW < kNbDimuonSet; iMassW++){
	    if(mass_Reco_QQ > massMIN[iMassW] && mass_Reco_QQ < massMAX[iMassW]){
	      //pT histograms
	      hReco_pT[iMassW][chargeID][Reco_QQ_type[iQQ]][iCut]->Fill(pT_Reco_QQ);
	      if(chargeID == 1 || chargeID == 2) //LS histo
		hReco_pT[iMassW][3][Reco_QQ_type[iQQ]][iCut]->Fill(pT_Reco_QQ);
	      //rap histograms
	      hReco_rap[iMassW][chargeID][Reco_QQ_type[iQQ]][iCut]->Fill(rap_Reco_QQ);
	      if(chargeID == 1 || chargeID == 2) //LS histo
		hReco_rap[iMassW][3][Reco_QQ_type[iQQ]][iCut]->Fill(rap_Reco_QQ);
	    
	      //correlating the single muons:
	      if(etaCorr)
		hReco_eta1_eta2[iMassW][chargeID][Reco_QQ_type[iQQ]][iCut]->Fill(Reco_QQ_posMu->Eta(),Reco_QQ_negMu->Eta());
    
	      // fill single muon quantities, after each successive cut
	      // prepare single muon indices to avoid double counting
	      int glMuID[10000], trMuID[10000],caMuID[10000];
	      for (int i=0; i<10000; i++){
		glMuID[i] = -1;
		trMuID[i] = -1;
		caMuID[i] = -1;
	      }

	      //1.) before first cut
	      if (iCut==0) { // comp, before d0 and dz
		switch (Reco_QQ_type[iQQ]) {
		case 0: // global+global
		  if(glMuID[Reco_QQ_muhpt[iQQ]]!=1){
		    hD0[iMassW][chargeID][0]->Fill(Reco_mu_glb_d0[Reco_QQ_muhpt[iQQ]]);
		    hDz[iMassW][chargeID][0]->Fill(Reco_mu_glb_dz[Reco_QQ_muhpt[iQQ]]);
		    glMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (glMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hD0[iMassW][chargeID][0]->Fill(Reco_mu_glb_d0[Reco_QQ_mulpt[iQQ]]);
		    hDz[iMassW][chargeID][0]->Fill(Reco_mu_glb_dz[Reco_QQ_mulpt[iQQ]]);
		    glMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 1: // global+tracker
		  if(glMuID[Reco_QQ_muhpt[iQQ]]!=1){
		    hD0[iMassW][chargeID][0]->Fill(Reco_mu_glb_d0[Reco_QQ_muhpt[iQQ]]);
		    hDz[iMassW][chargeID][0]->Fill(Reco_mu_glb_dz[Reco_QQ_muhpt[iQQ]]);
		    glMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (trMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hD0[iMassW][chargeID][1]->Fill(Reco_mu_trk_d0[Reco_QQ_mulpt[iQQ]]);
		    hDz[iMassW][chargeID][1]->Fill(Reco_mu_trk_dz[Reco_QQ_mulpt[iQQ]]);
		    trMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 2: // tracker+tracker
		  if (trMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    hD0[iMassW][chargeID][1]->Fill(Reco_mu_trk_d0[Reco_QQ_muhpt[iQQ]]);
		    hDz[iMassW][chargeID][1]->Fill(Reco_mu_trk_dz[Reco_QQ_muhpt[iQQ]]);
		    trMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (trMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hD0[iMassW][chargeID][1]->Fill(Reco_mu_trk_d0[Reco_QQ_mulpt[iQQ]]);
		    hDz[iMassW][chargeID][1]->Fill(Reco_mu_trk_dz[Reco_QQ_mulpt[iQQ]]);
		    trMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 3: // global+calo
		  if(glMuID[Reco_QQ_muhpt[iQQ]]!=1){
		    hD0[iMassW][chargeID][0]->Fill(Reco_mu_glb_d0[Reco_QQ_muhpt[iQQ]]);
		    hDz[iMassW][chargeID][0]->Fill(Reco_mu_glb_dz[Reco_QQ_muhpt[iQQ]]);
		    glMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (caMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hD0[iMassW][chargeID][2]->Fill(Reco_mu_cal_d0[Reco_QQ_mulpt[iQQ]]);
		    hDz[iMassW][chargeID][2]->Fill(Reco_mu_cal_dz[Reco_QQ_mulpt[iQQ]]);
		    caMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 4: // tracker+calo
		  if (trMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    hD0[iMassW][chargeID][1]->Fill(Reco_mu_trk_d0[Reco_QQ_muhpt[iQQ]]);
		    hDz[iMassW][chargeID][1]->Fill(Reco_mu_trk_dz[Reco_QQ_muhpt[iQQ]]);
		    trMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (caMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hD0[iMassW][chargeID][2]->Fill(Reco_mu_cal_d0[Reco_QQ_mulpt[iQQ]]);
		    hDz[iMassW][chargeID][2]->Fill(Reco_mu_cal_dz[Reco_QQ_mulpt[iQQ]]);
		    caMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 5: // calo+calo
		  if (caMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    hD0[iMassW][chargeID][2]->Fill(Reco_mu_cal_d0[Reco_QQ_muhpt[iQQ]]);
		    hDz[iMassW][chargeID][2]->Fill(Reco_mu_cal_dz[Reco_QQ_muhpt[iQQ]]);
		    caMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (caMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hD0[iMassW][chargeID][2]->Fill(Reco_mu_cal_d0[Reco_QQ_mulpt[iQQ]]);
		    hDz[iMassW][chargeID][2]->Fill(Reco_mu_cal_dz[Reco_QQ_mulpt[iQQ]]);
		    caMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		default:
		  break;
		}
	      }

	      //2.) after d0 and dz cut
	      if (iCut==1) {
		switch (Reco_QQ_type[iQQ]) {
		case 0: // global+global
		  break;
		case 1: // global+tracker
		  if (trMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hTM2DCompatibilityTight[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_mulpt[iQQ]] & (int) pow(2,5))/(int) pow(2,5));
		    hTMLastStationOptimizedLowPtLoose[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_mulpt[iQQ]] & (int) pow(2,8))/(int) pow(2,8));
		    hTrackerMuonArbitrated[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_mulpt[iQQ]] & (int)pow(2,1))/(int)pow(2,1));
		    hCaloCompatibility[iMassW][chargeID][1]->Fill(Reco_mu_trk_caloComp[Reco_QQ_mulpt[iQQ]]);
		    trMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 2: // tracker+tracker
		  if (trMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    hTM2DCompatibilityTight[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_muhpt[iQQ]] & (int) pow(2,5))/(int) pow(2,5));
		    hTMLastStationOptimizedLowPtLoose[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_muhpt[iQQ]] & (int) pow(2,8))/(int) pow(2,8));
		    hTrackerMuonArbitrated[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_muhpt[iQQ]] & (int)pow(2,1))/(int)pow(2,1));

		    hCaloCompatibility[iMassW][chargeID][1]->Fill(Reco_mu_trk_caloComp[Reco_QQ_muhpt[iQQ]]);
		    trMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (trMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hTM2DCompatibilityTight[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_mulpt[iQQ]] & (int) pow(2,5))/(int) pow(2,5));
		    hTMLastStationOptimizedLowPtLoose[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_mulpt[iQQ]] & (int) pow(2,8))/(int) pow(2,8));
		    hTrackerMuonArbitrated[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_mulpt[iQQ]] & (int)pow(2,1))/(int)pow(2,1));

		    hCaloCompatibility[iMassW][chargeID][1]->Fill(Reco_mu_trk_caloComp[Reco_QQ_mulpt[iQQ]]);
		    trMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 3: // global+calo
		  if (caMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    if (Reco_mu_cal_normChi2[Reco_QQ_mulpt[iQQ]] < MAX_normchi2_cal &&
			Reco_mu_cal_nhitstrack[Reco_QQ_mulpt[iQQ]] > MIN_nhits_trk)
		      hCaloCompatibility[iMassW][chargeID][2]->Fill(Reco_mu_cal_caloComp[Reco_QQ_mulpt[iQQ]]);
		    caMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 4: // tracker+calo
		  if (trMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    hTM2DCompatibilityTight[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_muhpt[iQQ]] & (int) pow(2,5))/(int) pow(2,5));
		    hTMLastStationOptimizedLowPtLoose[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_muhpt[iQQ]] & (int) pow(2,8))/(int) pow(2,8));
		    hTrackerMuonArbitrated[iMassW][chargeID]->Fill((Reco_mu_trk_PIDmask[Reco_QQ_muhpt[iQQ]] & (int)pow(2,1))/(int)pow(2,1));
		    hCaloCompatibility[iMassW][chargeID][1]->Fill(Reco_mu_trk_caloComp[Reco_QQ_muhpt[iQQ]]);
		    trMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (caMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    if (Reco_mu_cal_normChi2[Reco_QQ_mulpt[iQQ]] < MAX_normchi2_cal &&
			Reco_mu_cal_nhitstrack[Reco_QQ_mulpt[iQQ]] > MIN_nhits_trk)
		      hCaloCompatibility[iMassW][chargeID][2]->Fill(Reco_mu_cal_caloComp[Reco_QQ_mulpt[iQQ]]);
		    caMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 5: // calo+calo
		  if (caMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    if (Reco_mu_cal_normChi2[Reco_QQ_muhpt[iQQ]] < MAX_normchi2_cal &&
			Reco_mu_cal_nhitstrack[Reco_QQ_muhpt[iQQ]] > MIN_nhits_trk) {
		      hCaloCompatibility[iMassW][chargeID][2]->Fill(Reco_mu_cal_caloComp[Reco_QQ_muhpt[iQQ]]);
		      caMuID[Reco_QQ_muhpt[iQQ]]=1;
		    }
		  }
		  if (caMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    if (Reco_mu_cal_normChi2[Reco_QQ_mulpt[iQQ]] < MAX_normchi2_cal &&
			Reco_mu_cal_nhitstrack[Reco_QQ_mulpt[iQQ]] > MIN_nhits_trk) {
		      hCaloCompatibility[iMassW][chargeID][2]->Fill(Reco_mu_cal_caloComp[Reco_QQ_mulpt[iQQ]]);
		      caMuID[Reco_QQ_mulpt[iQQ]]=1;
		    }
		  }
		  break;
		default:
		  break;
		}
	      }
	      
	      //3.) after hCaloComp, TM2DCompatibilityTight, TMLastStationOptimizedLowPtTight
	      if (iCut==2) {
		switch (Reco_QQ_type[iQQ]) {
		case 0: // global+global
		  if (glMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    hChi2GlobalFit[iMassW][chargeID]->Fill(Reco_mu_glb_normChi2[Reco_QQ_muhpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][0]->Fill(Reco_mu_glb_nhitstrack[Reco_QQ_muhpt[iQQ]]);
		    glMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (glMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hChi2GlobalFit[iMassW][chargeID]->Fill(Reco_mu_glb_normChi2[Reco_QQ_mulpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][0]->Fill(Reco_mu_glb_nhitstrack[Reco_QQ_mulpt[iQQ]]);
		    glMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 1: // global+tracker
		  if (glMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    hChi2GlobalFit[iMassW][chargeID]->Fill(Reco_mu_glb_normChi2[Reco_QQ_muhpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][0]->Fill(Reco_mu_glb_nhitstrack[Reco_QQ_muhpt[iQQ]]);
		    glMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  if (trMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hChi2TrackerFit[iMassW][chargeID][1]->Fill(Reco_mu_trk_normChi2[Reco_QQ_mulpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][1]->Fill(Reco_mu_trk_nhitstrack[Reco_QQ_mulpt[iQQ]]);
		    trMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 2: // tracker+tracker
		  if (trMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    hChi2TrackerFit[iMassW][chargeID][1]->Fill(Reco_mu_trk_normChi2[Reco_QQ_muhpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][1]->Fill(Reco_mu_trk_nhitstrack[Reco_QQ_muhpt[iQQ]]);
		    trMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (trMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hChi2TrackerFit[iMassW][chargeID][1]->Fill(Reco_mu_trk_normChi2[Reco_QQ_mulpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][1]->Fill(Reco_mu_trk_nhitstrack[Reco_QQ_mulpt[iQQ]]);
		    trMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 3: // global+calo
		  if (glMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    hChi2GlobalFit[iMassW][chargeID]->Fill(Reco_mu_glb_normChi2[Reco_QQ_muhpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][0]->Fill(Reco_mu_glb_nhitstrack[Reco_QQ_muhpt[iQQ]]);
		    glMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (caMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hChi2TrackerFit[iMassW][chargeID][2]->Fill(Reco_mu_cal_normChi2[Reco_QQ_mulpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][2]->Fill(Reco_mu_cal_nhitstrack[Reco_QQ_mulpt[iQQ]]);
		    caMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 4: // tracker+calo
		  if (trMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    hChi2TrackerFit[iMassW][chargeID][1]->Fill(Reco_mu_trk_normChi2[Reco_QQ_muhpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][1]->Fill(Reco_mu_trk_nhitstrack[Reco_QQ_muhpt[iQQ]]);
		    trMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (caMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hChi2TrackerFit[iMassW][chargeID][2]->Fill(Reco_mu_cal_normChi2[Reco_QQ_mulpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][2]->Fill(Reco_mu_cal_nhitstrack[Reco_QQ_mulpt[iQQ]]);
		    caMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		case 5: // calo+calo
		  if (caMuID[Reco_QQ_muhpt[iQQ]]!=1) {
		    hChi2TrackerFit[iMassW][chargeID][2]->Fill(Reco_mu_cal_normChi2[Reco_QQ_muhpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][2]->Fill(Reco_mu_cal_nhitstrack[Reco_QQ_muhpt[iQQ]]);
		    caMuID[Reco_QQ_muhpt[iQQ]]=1;
		  }
		  if (caMuID[Reco_QQ_mulpt[iQQ]]!=1) {
		    hChi2TrackerFit[iMassW][chargeID][2]->Fill(Reco_mu_cal_normChi2[Reco_QQ_mulpt[iQQ]]);
		    hNHitsSilicon[iMassW][chargeID][2]->Fill(Reco_mu_cal_nhitstrack[Reco_QQ_mulpt[iQQ]]);
		    caMuID[Reco_QQ_mulpt[iQQ]]=1;
		  }
		  break;
		default:
		  break;
		}
	      }
	      
	      //4.) after chi2 and #hits
	      if (iCut==3) { // dimuon vertex prob, after comp
		hDimuVtxProb[iMassW][chargeID][Reco_QQ_type[iQQ]]->Fill(Reco_QQ_probChi2[iQQ]);
		if(chargeID == 1 || chargeID == 2) //LS histo
		  hDimuVtxProb[iMassW][3][Reco_QQ_type[iQQ]]->Fill(Reco_QQ_probChi2[iQQ]);
	      }
	    } // fill for signal
	  }//different mass windows
	}//cuts   
    }//#Reco_QQ combinations

    for(int iCharge = 0; iCharge < kNbCharge; iCharge++)
      for(int iCat = 0; iCat < kNbCath; iCat++)
	for(int iCut = 0; iCut < kNbCuts; iCut++)
	  if(nComb[iCharge][iCat][iCut] > 0)
	    hMult_QQ[chargeID][iCat][iCut]->Fill(nComb[iCharge][iCat][iCut]);
    for(int iCat = 0; iCat < kNbCath; iCat++)
      for(int iCut = 0; iCut < kNbCuts; iCut++)
	hMult_QQ[3][iCat][iCut]->Fill(nComb[1][iCat][iCut]+nComb[2][iCat][iCut]);

    //now fill the histos for the "best" QQ pair only
    Int_t iDBestQQ = -1;
    if(Reco_QQ_size > 0){
      countAllRecoQQEvents++;
      iDBestQQ = theBestQQ();
      // 	if(iDBestQQ >= 0){
      // 	  printf("best QQ is %d (out of %d) and is of type %s\n",
      // 		 iDBestQQ, Reco_QQ_size, oniaCatName[Reco_QQ_type[iDBestQQ]]);
      // 	}
      // 	else{
      // 	  printf("best QQ not identified... %d total QQbar\n", Reco_QQ_size);
      // 	  countNotFoundBestQQEvents++;
      // 	}
      if(iDBestQQ < 0){
	countNotFoundBestQQEvents++;
	for(int iRec = 0; iRec < Reco_QQ_size; iRec++){
	    
	  if(Reco_QQ_sign[iRec] != 0) continue;
	  //we count ALL reco QQpairs in the event:
	  hRecoNoBest->Fill(Reco_QQ_type[iRec]);
	}
      }
      else{

	TLorentzVector *Reco_QQ = (TLorentzVector *) Reco_QQ_4mom->At(iDBestQQ);
	Double_t mass_Reco_QQ = Reco_QQ->M();
	Double_t pT_Reco_QQ = Reco_QQ->Pt();
	Double_t rap_Reco_QQ = Reco_QQ->Rapidity();

	chargeID = 0;
	if(Reco_QQ_sign[iDBestQQ] < 0) chargeID = 1;
	else if(Reco_QQ_sign[iDBestQQ] > 0) chargeID = 2;

	Bool_t isMatched[2] = {kTRUE, kTRUE};//only relevant for MC; for data must be set to kTRUE
	if (matchMC){

	  isMatched[0] = kFALSE;
	  isMatched[1] = kFALSE;

	  if(Mc_QQmupl_indx[0] < 2 && Mc_QQmumi_indx[0] < 2){ //sanity check

	    switch (Reco_QQ_type[iDBestQQ]) {
	    case 0: // global-global
	      //first MC muon
	      if(Reco_QQ_muhpt[iDBestQQ] < Reco_mu_glb_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iDBestQQ] < Reco_mu_glb_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    case 1: // global-tracker
	      //first MC muon
	      if(Reco_QQ_muhpt[iDBestQQ] < Reco_mu_glb_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iDBestQQ] < Reco_mu_trk_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    case 2: // tracker-tracker
	      //first MC muon
	      if(Reco_QQ_muhpt[iDBestQQ] < Reco_mu_trk_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iDBestQQ] < Reco_mu_trk_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    case 3: // global-calo
	      //first MC muon
	      if(Reco_QQ_muhpt[iDBestQQ] < Reco_mu_glb_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iDBestQQ] < Reco_mu_cal_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    case 4: // tracker-calo
	      //first MC muon
	      if(Reco_QQ_muhpt[iDBestQQ] < Reco_mu_trk_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iDBestQQ] < Reco_mu_cal_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    case 5: // calo-calo
	      //first MC muon
	      if(Reco_QQ_muhpt[iDBestQQ] < Reco_mu_cal_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iDBestQQ] < Reco_mu_cal_size){//sanity check
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
	      }
	      break;
	    default:
	      break;
	    }	    
	  }
	}
	//check the MC matching before continuing...
	if(isMatched[0] == kTRUE && isMatched[1] == kTRUE){

	  //mass histograms
	  hRecoBest_mass[chargeID][Reco_QQ_type[iDBestQQ]]->Fill(mass_Reco_QQ);
	  if(chargeID == 1 || chargeID == 2) //LS histo
	    hRecoBest_mass[3][Reco_QQ_type[iDBestQQ]]->Fill(mass_Reco_QQ);

	  //fill pT and rap histos only for the signal:
	  for(int iMassW = 0; iMassW < kNbDimuonSet; iMassW++){
	    if(mass_Reco_QQ > massMIN[iMassW] && mass_Reco_QQ < massMAX[iMassW]){
	      //pT histograms
	      hRecoBest_pT[iMassW][chargeID][Reco_QQ_type[iDBestQQ]]->Fill(pT_Reco_QQ);
	      if(chargeID == 1 || chargeID == 2) //LS histo
		hRecoBest_pT[iMassW][3][Reco_QQ_type[iDBestQQ]]->Fill(pT_Reco_QQ);
	      //rap histograms
	      hRecoBest_rap[iMassW][chargeID][Reco_QQ_type[iDBestQQ]]->Fill(rap_Reco_QQ);
	      if(chargeID == 1 || chargeID == 2) //LS histo
		hRecoBest_rap[iMassW][3][Reco_QQ_type[iDBestQQ]]->Fill(rap_Reco_QQ);
	    }
	    //===============================================================
	    //single muon histograms for best dimuon within 3.0 < M < 3.2 GeV
	    //and for dimuons with M > 2 GeV
	    //===============================================================
	    TLorentzVector *Reco_mu_1, *Reco_mu_2;
	    Int_t iMuCat[2] = {-1,-1};
	    if(Reco_QQ_type[iDBestQQ] == 0){//gl-gl
	      Reco_mu_1 = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
	      Reco_mu_2 = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
	      iMuCat[0] = 0; iMuCat[1] = 0;
	    }
	    else if(Reco_QQ_type[iDBestQQ] == 1){//gl-tr
	      Reco_mu_1 = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
	      Reco_mu_2 = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
	      iMuCat[0] = 0; iMuCat[1] = 1;
	    }
	    else if(Reco_QQ_type[iDBestQQ] == 2){//tr-tr
	      Reco_mu_1 = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
	      Reco_mu_2 = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
	      iMuCat[0] = 1; iMuCat[1] = 1;
	    }
	    else if(Reco_QQ_type[iDBestQQ] == 3){//gl-calo
	      Reco_mu_1 = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
	      Reco_mu_2 = (TLorentzVector *) Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
	      iMuCat[0] = 0; iMuCat[1] = 2;
	    }
	    else if(Reco_QQ_type[iDBestQQ] == 4){//tr-calo
	      Reco_mu_1 = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
	      Reco_mu_2 = (TLorentzVector *) Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
	      iMuCat[0] = 1; iMuCat[1] = 2;
	    }
	    else if(Reco_QQ_type[iDBestQQ] == 5){//calo-calo
	      Reco_mu_1 = (TLorentzVector *) Reco_mu_cal_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
	      Reco_mu_2 = (TLorentzVector *) Reco_mu_cal_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
	      iMuCat[0] = 2; iMuCat[1] = 2;
	    }
	    
	    Double_t etaMu[2] = {Reco_mu_1->Eta(), Reco_mu_2->Eta()};
	    Double_t pTMu[2]  = {Reco_mu_1->Pt(),	 Reco_mu_2->Pt()};
	    for(int iMu = 0; iMu < 2; iMu++){
	      hMuon_eta[iMassW][chargeID][iMuCat[iMu]]->Fill(etaMu[iMu]);
	      hMuon_pT[iMassW][chargeID][iMuCat[iMu]]->Fill(pTMu[iMu]);
	      hMuon_pT_eta[iMassW][chargeID][iMuCat[iMu]]->Fill(etaMu[iMu], pTMu[iMu]);
	    }
	  }//mass windows
	}//matched: relevant for MC
      }//best QQ
    }//Reco_QQ_size > 0
  }//events

  if(countAllRecoQQEvents > 0){
    hRecoNoBest->Scale(1./countAllRecoQQEvents);
    printf("\n\n The number of events w/o identification of the best QQbar is %d out of %d, i.e. %1.1f %%\n", countNotFoundBestQQEvents, countAllRecoQQEvents, 100.*countNotFoundBestQQEvents/countAllRecoQQEvents);
  }
  }

//=============================
int ProjectQQ::theBestQQ() {
    
  int theBest = -1;
  float thehighestPt = -1.;
 
  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 0 ) {

      int thehptMu = Reco_QQ_muhpt[iqq];   if (thehptMu >= Reco_mu_glb_size) continue;
      int thelptMu = Reco_QQ_mulpt[iqq];   if (thelptMu >= Reco_mu_glb_size) continue;
      if (Reco_QQ_probChi2[iqq] > MIN_vtxprob_jpsi && 
	  //Reco_mu_glb_nhitstrack[thehptMu] > MIN_nhits_trk && 
	  Reco_mu_glb_normChi2[thehptMu] < MAX_normchi2_glb && 
	  //(((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel) || ((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thehptMu] == 1)) && 
	  fabs(Reco_mu_glb_d0[thehptMu]) < MAX_d0_trk && 
	  fabs(Reco_mu_glb_dz[thehptMu]) < MAX_dz_trk && 
	  //Reco_mu_glb_nhitstrack[thelptMu] > MIN_nhits_trk && 
	  Reco_mu_glb_normChi2[thelptMu] < MAX_normchi2_glb &&
	  //(((Reco_mu_glb_nhitsPixB[thelptMu] + Reco_mu_glb_nhitsPixE[thelptMu]) > MIN_nhits_pixel) || ((Reco_mu_glb_nhitsPixB[thelptMu] + Reco_mu_glb_nhitsPixE[thelptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thelptMu] == 1)) &&
	  fabs(Reco_mu_glb_d0[thelptMu]) < MAX_d0_trk && 
	  fabs(Reco_mu_glb_dz[thelptMu]) < MAX_dz_trk
	  ) {
	return iqq;
      }
    }
  }

  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 1 ) {
      
      int thehptMu = Reco_QQ_muhpt[iqq];  if (thehptMu >= Reco_mu_glb_size) continue;
      int thelptMu = Reco_QQ_mulpt[iqq];  if (thelptMu >= Reco_mu_trk_size) continue;

      if ( Reco_QQ_probChi2[iqq] > MIN_vtxprob_jpsi &&
	   //Reco_mu_glb_nhitstrack[thehptMu] > MIN_nhits_trk && 
	   Reco_mu_glb_normChi2[thehptMu] < MAX_normchi2_glb &&
	   //(((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel) || ((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thehptMu] == 1)) &&
	   fabs(Reco_mu_glb_d0[thehptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_glb_dz[thehptMu]) < MAX_dz_trk && 
	   Reco_mu_trk_nhitstrack[thelptMu] > MIN_nhits_trk && 
	   ((Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,5))/(int)pow(2,5) > 0 || (Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,8))/(int)pow(2,8) > 0) &&
	   //(Reco_mu_trk_nhitsPixB[thelptMu] + Reco_mu_trk_nhitsPixE[thelptMu]) > MIN_nhits_pixel &&
	   Reco_mu_trk_normChi2[thelptMu] < MAX_normchi2_trk &&
	   fabs(Reco_mu_trk_d0[thelptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_trk_dz[thelptMu]) < MAX_dz_trk) {
	
        TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(thelptMu);
        if (theTrMumom->Perp() > thehighestPt) {
	  thehighestPt = theTrMumom->Perp();
          theBest = iqq;
	}
      }
    }    
  }
  
  if (theBest >= 0) return theBest;

  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 2 ) {
      
      int thehptMu = Reco_QQ_muhpt[iqq];  if (thehptMu >= Reco_mu_trk_size) continue;
      int thelptMu = Reco_QQ_mulpt[iqq];  if (thelptMu >= Reco_mu_trk_size) continue;

      if ( Reco_QQ_probChi2[iqq] > MIN_vtxprob_jpsi &&
	   Reco_mu_trk_nhitstrack[thehptMu] > MIN_nhits_trk && 
	   ((Reco_mu_trk_PIDmask[thehptMu] & (int)pow(2,5))/(int)pow(2,5) > 0 || (Reco_mu_trk_PIDmask[thehptMu] & (int)pow(2,8))/(int)pow(2,8) > 0) &&
	   //(Reco_mu_trk_nhitsPixB[thehptMu] + Reco_mu_trk_nhitsPixE[thehptMu]) > MIN_nhits_pixel &&
	   Reco_mu_trk_normChi2[thehptMu] < MAX_normchi2_trk &&
	   fabs(Reco_mu_trk_d0[thehptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_trk_dz[thehptMu]) < MAX_dz_trk &&
	   Reco_mu_trk_nhitstrack[thelptMu] > MIN_nhits_trk && 
	   ((Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,5))/(int)pow(2,5) > 0 || (Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,8))/(int)pow(2,8) > 0) &&
	   //(Reco_mu_trk_nhitsPixB[thelptMu] + Reco_mu_trk_nhitsPixE[thelptMu]) > MIN_nhits_pixel &&
	   Reco_mu_trk_normChi2[thelptMu] < MAX_normchi2_trk &&
	   fabs(Reco_mu_trk_d0[thelptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_trk_dz[thelptMu]) < MAX_dz_trk) {
	
        TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(thehptMu);
        if (theTrMumom->Perp() > thehighestPt) {
	  thehighestPt = theTrMumom->Perp();
          theBest = iqq;
	}
      }
    }    
  }
  
  if (theBest >= 0) return theBest;

  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 3 ) {

      int thehptMu = Reco_QQ_muhpt[iqq];
      int thelptMu = Reco_QQ_mulpt[iqq];
      if (thelptMu >= Reco_mu_cal_size) {
	// cout << "Non deve succedere! cmIndex = " << thelptMu+1 << " cmSize = " << Reco_mu_cal_size << endl;
	continue;
      }

      if ( //Reco_mu_glb_nhitstrack[thehptMu] > MIN_nhits_trk && 
	   Reco_mu_glb_normChi2[thehptMu] < MAX_normchi2_glb &&
	   Reco_mu_cal_nhitstrack[thelptMu] > MIN_nhits_trk && 
	   Reco_mu_cal_normChi2[thelptMu] < MAX_normchi2_cal && 
	   Reco_mu_cal_caloComp[thelptMu] > MIN_caloComp) {
	
        TLorentzVector *theCaMumom = (TLorentzVector*)Reco_mu_cal_4mom->At(thelptMu);
        if (theCaMumom->Perp() > thehighestPt) {
	  thehighestPt = theCaMumom->Perp();
          theBest = iqq;
	}
      }
    }    
  }
  
  return theBest;
}
//=========================================
Double_t ProjectQQ::deltaR(TLorentzVector* t, TLorentzVector* u){

  return sqrt(pow(t->Eta()-u->Eta(),2) +pow(PhiInRange(t->Phi()-u->Phi()),2));
}
//=========================================
double ProjectQQ::PhiInRange(double phi){

      double phiout = phi;

      if( phiout > 2*M_PI || phiout < -2*M_PI) {
            phiout = fmod( phiout, 2*M_PI);
      }
      if (phiout <= -M_PI) phiout += 2*M_PI;
      else if (phiout >  M_PI) phiout -= 2*M_PI;

      return phiout;
}
