#define GeomAcc_cxx
#include "GeomAcc.h"
#include "calcPol.C"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TRandom.h"

using namespace std;

Double_t massMuOnia;

TH1D *hMass[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH1D *hMass_Smeared[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];

TH2D *hGenPtRap, *hGenCutPtRap;

//================================================
//2D polarization histos, differential in pT and y
//================================================
//histos for neg. and pos. rapidity separately:
TH2D *hGen_pT_rapNP[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2D *hGenCut_pT_rapNP[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
//histos taking together +y and -y:
TH2D *hGen_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2D *hGenCut_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
//histos taking together +y and -y and 4-folding in phi
TH2D *hGen_pT_rap_phiFolded[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2D *hGenCut_pT_rap_phiFolded[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];

enum {LOOSE,TIGHT};//set of muon fiducial cuts
//================================================
Double_t smearValues[27] = {1.27851,
			    0.00679612,
			    0.,
			    0.00630424,
			    0.000557345,
			    0.,
			    0.0145714,
			    0.860521,
			    0.000990502,
			    -0.00023,
			    -0.0005,
			    0.00107,
			    5.80784e-05,
			    0.000101772,
			    1.96263,
			    2.18953,
			    0.0257171,
			    0.,
			    0.160361,
			    2.26713,
			    0.,0.,0.,0.,0.,0.,0.};

//================================
void GeomAcc::Loop(Bool_t smearing, Bool_t rejectCowboys){
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

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

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     if(jentry%10000 == 0) printf("event %d / %d\n", (Int_t) jentry, (Int_t) nentries);
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //====================================
      //get the *generated* onia (w/o FSR!)
      Double_t enOniaGen = sqrt(onia_Px*onia_Px + onia_Py*onia_Py + onia_Pz*onia_Pz + massMuOnia*massMuOnia);
      TLorentzVector *oniaGen = new TLorentzVector();
      oniaGen->SetPxPyPzE(onia_Px, onia_Py, onia_Pz, enOniaGen);
      Double_t rapOniaGen = oniaGen->Rapidity();
      Double_t pTOniaGen = oniaGen->Pt();
      Double_t massOniaGen = oniaGen->M();
      //====================================

      //=================================================
      //get the individual muons that are affected by FSR
      Double_t enMuPos = sqrt(muPos_Px*muPos_Px + muPos_Py*muPos_Py + muPos_Pz*muPos_Pz + jpsi::muMass*jpsi::muMass);
      Double_t enMuNeg = sqrt(muNeg_Px*muNeg_Px + muNeg_Py*muNeg_Py + muNeg_Pz*muNeg_Pz + jpsi::muMass*jpsi::muMass);
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

      //build the invariant mass, pt, ... of the two muons
      TLorentzVector *oniaFSR = new TLorentzVector();
      *oniaFSR = *(muPos) + *(muNeg);
      Double_t rapOniaFSR = oniaFSR->Rapidity();
      Double_t pTOniaFSR = oniaFSR->Pt();
      Double_t massOniaFSR = oniaFSR->M();

      Int_t rapIndex = -1;
      for(int iRap = 0; iRap < 2*jpsi::kNbRapBins; iRap++){
	if(rapOniaFSR > jpsi::rapRange[iRap] && rapOniaFSR < jpsi::rapRange[iRap+1]){
	  rapIndex = iRap;
	  break;
	}
      }
      Int_t rapForPTIndex = -1;
      for(int iRap = 0; iRap < jpsi::kNbRapForPTBins; iRap++){
	if(TMath::Abs(rapOniaFSR) > jpsi::rapForPTRange[iRap] && 
	   TMath::Abs(rapOniaFSR) < jpsi::rapForPTRange[iRap+1]){
	  rapForPTIndex = iRap+1;
	  break;
	}
      }
      Int_t pTIndex = -1;
      for(int iPT = 0; iPT < jpsi::kNbPTBins[rapForPTIndex]; iPT++){
	if(pTOniaFSR > jpsi::pTRange[rapForPTIndex][iPT] && pTOniaFSR < jpsi::pTRange[rapForPTIndex][iPT+1]){
	  pTIndex = iPT+1;
	  break;
	}
      }
      Int_t rapIntegratedPTIndex = -1;
      for(int iPT = 0; iPT < jpsi::kNbPTBins[0]; iPT++){
	if(pTOniaFSR > jpsi::pTRange[0][iPT] && pTOniaFSR < jpsi::pTRange[0][iPT+1]){
	  rapIntegratedPTIndex = iPT+1;
	  break;
	}
      }

      //select events within a narrow mass window around the J/psi
      //(rapidity dependence of the resolution --> different mass windows)
      //values obtained from Gaussian fits in "plotMass.C"
      Double_t jPsiMassMin = jpsi::polMassJpsi[rapForPTIndex] - jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
      Double_t jPsiMassMax = jpsi::polMassJpsi[rapForPTIndex] + jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
      if(massOniaFSR < jPsiMassMin || massOniaFSR > jPsiMassMax)
	continue;

      if(rejectCowboys)
	if((muNeg->Phi() - muPos->Phi()) < 0.)
	  continue;

      calcPol(*muPos, *muNeg);  //get the polarization variables at the generation level (with FSR muons)

      //assign a weight according to a given pT and y distribution (taken from PYTHIA)
      //PYTHIA distributions fit with the macro "projectPTRap_Pythia.C(kTRUE)"
      //pT distribtuions do NOT show any rapidity dependence, while for pT > 10 GeV/c
      //there seems to be some rapidity dependence which is neglected in this parameterization
      Double_t pTweight = fPT->Eval(pTOniaFSR);
      Double_t rapweight = fRap->Eval(rapOniaFSR);
      Double_t kinWeight = pTweight*rapweight;
//       printf("pT weight %1.3f, y weight %1.3f\n", pTweight, rapweight);
      
      //fill the generator histos only if both muons are within |eta| < 2.4
//       if(fabs(etaMuPos) < jpsi::etaPS[2] && fabs(etaMuNeg) < jpsi::etaPS[2]){ //H: 5March2011 we should only limit the J/psi's kinematics!
      if(fabs(rapOniaFSR) < jpsi::rapYPS){

	hGenPtRap->Fill(rapOniaFSR, pTOniaFSR, kinWeight);

	if(pTIndex >= 0 && rapForPTIndex >= 0)
	  hMass[pTIndex][rapForPTIndex]->Fill(massOniaFSR, kinWeight);
	if(rapForPTIndex > 0)
	  hMass[0][rapForPTIndex]->Fill(massOniaFSR, kinWeight);
	if(pTIndex >= 0)
	  hMass[pTIndex][0]->Fill(massOniaFSR, kinWeight);
	hMass[0][0]->Fill(massOniaFSR, kinWeight);

	for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){

	  Double_t weight = CalcPolWeight(thisCosTh[iFrame]);
	  weight *= kinWeight;

	  //histos for neg. and pos. rapidity separately:
	  if(rapIndex >= 0)
	    hGen_pT_rapNP[iFrame][0][rapIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	  if(pTIndex > 0 && rapIndex >= 0)
	    hGen_pT_rapNP[iFrame][pTIndex][rapIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	
	  //histos taking together +y and -y
	  hGen_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	  if(rapIntegratedPTIndex > 0)
	    hGen_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	  if(rapForPTIndex > 0)
	    hGen_pT_rap[iFrame][0][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	  if(pTIndex > 0 && rapForPTIndex > 0)
	    hGen_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	  
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
	  hGen_pT_rap_phiFolded[iFrame][0][0]->Fill(thetaAdjusted, phiFolded, weight);
	  if(rapIntegratedPTIndex > 0)
	    hGen_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex][0]->Fill(thetaAdjusted, phiFolded, weight);
	  if(rapForPTIndex > 0)
	    hGen_pT_rap_phiFolded[iFrame][0][rapForPTIndex]->Fill(thetaAdjusted, phiFolded, weight);
	  if(pTIndex > 0 && rapForPTIndex > 0)
	    hGen_pT_rap_phiFolded[iFrame][pTIndex][rapForPTIndex]->Fill(thetaAdjusted, phiFolded, weight);
	}
      }

      //now, smear the two muons with the detector resolution
      //==================================================
      TLorentzVector *muPosSmeared, *muNegSmeared;
      if(smearing){
	muPosSmeared = ApplySmearing(muPos, smearValues);
	muNegSmeared  = ApplySmearing(muNeg, smearValues);
      }
      else{
	muPosSmeared = new TLorentzVector(*muPos);
	muNegSmeared = new TLorentzVector(*muNeg);
      }
      
      Double_t etaMuPosSmeared = muPosSmeared->Eta();
      Double_t etaMuNegSmeared = muNegSmeared->Eta();
      Double_t pTMuPosSmeared = muPosSmeared->Pt();
      Double_t pTMuNegSmeared = muNegSmeared->Pt();
      Double_t pMuPosSmeared = muPosSmeared->P();
      Double_t pMuNegSmeared = muNegSmeared->P();

      // //continue only if both muons are within |eta| < 2.4:
      // if(fabs(etaMuPosSmeared) > jpsi::etaPS[2] || fabs(etaMuNegSmeared) > jpsi::etaPS[2])
      // 	continue;

      //old definition of acceptance
      // //apply the kinematical phase space cut:
      // //(a) on the positive muon
      // if((fabs(etaMuPosSmeared) < jpsi::etaPS[0] && pTMuPosSmeared < jpsi::pTMuMin[0]) || //mid-rapidity cut
      // 	 (fabs(etaMuPosSmeared) > jpsi::etaPS[0] && fabs(etaMuPosSmeared) < jpsi::etaPS[1] && pMuPosSmeared < jpsi::pMuMin[1]) ||
      // 	 (fabs(etaMuPosSmeared) > jpsi::etaPS[1] && fabs(etaMuPosSmeared) < jpsi::etaPS[2] && pTMuPosSmeared < jpsi::pTMuMin[2]))
      // 	continue;
      // //(b) on the negative muon
      // if((fabs(etaMuNegSmeared) < jpsi::etaPS[0] && pTMuNegSmeared < jpsi::pTMuMin[0]) || //mid-rapidity cut
      // 	 (fabs(etaMuNegSmeared) > jpsi::etaPS[0] && fabs(etaMuNegSmeared) < jpsi::etaPS[1] && pMuNegSmeared < jpsi::pMuMin[1]) ||
      // 	 (fabs(etaMuNegSmeared) > jpsi::etaPS[1] && fabs(etaMuNegSmeared) < jpsi::etaPS[2] && pTMuNegSmeared < jpsi::pTMuMin[2]))
      // 	continue;
      
      Bool_t muonsInAcc = kFALSE;
      if(pTMuPosSmeared > pTMuNegSmeared){
	muonsInAcc = areMuonsInAcceptance(LOOSE, pTMuPosSmeared, etaMuPosSmeared, pTMuNegSmeared, etaMuNegSmeared); //code needs to be adjusted for the DoubleMu0 trigger!!!!
	//muonsInAcc = areMuonsInAcceptance(TIGHT, pTMuPosSmeared, etaMuPosSmeared, pTMuNegSmeared, etaMuNegSmeared); //code needs to be adjusted for the DoubleMu0 trigger!!!!
      }
      else{
	muonsInAcc = areMuonsInAcceptance(LOOSE, pTMuNegSmeared, etaMuNegSmeared, pTMuPosSmeared, etaMuPosSmeared); //code needs to be adjusted for the DoubleMu0 trigger!!!!
	//muonsInAcc = areMuonsInAcceptance(TIGHT, pTMuNegSmeared, etaMuNegSmeared, pTMuPosSmeared, etaMuPosSmeared); //code needs to be adjusted for the DoubleMu0 trigger!!!!
      }
      if(!muonsInAcc)
	continue;

      //reject furthermore all events in which one of the muons has a pT smaller than 3 GeV/c
      if(pTMuPosSmeared < 3.0 || pTMuNegSmeared < 3.0)
	continue;

      //build the invariant mass, pt, ... of the two muons
      //includes FSR and smearing (if switched on)
      TLorentzVector *oniaSim = new TLorentzVector();
      *oniaSim = *(muPosSmeared) + *(muNegSmeared);

      Double_t oniaSim_mass = oniaSim->M();
      Double_t oniaSim_pt = oniaSim->Pt();
      Double_t oniaSim_P = oniaSim->P();
      Double_t oniaSim_eta = oniaSim->PseudoRapidity();
      Double_t oniaSim_rap = oniaSim->Rapidity();
      Double_t oniaSim_phi = oniaSim->Phi();
      Double_t oniaSim_mT = sqrt(oniaSim_mass*oniaSim_mass + oniaSim_pt*oniaSim_pt);

      hMass_Smeared[0][0]->Fill(oniaSim_mass, kinWeight); //others will be filled after all indices are set

      rapIndex = -1;
      for(int iRap = 0; iRap < 2*jpsi::kNbRapBins; iRap++){
	if(oniaSim_rap > jpsi::rapRange[iRap] && oniaSim_rap < jpsi::rapRange[iRap+1]){
	  rapIndex = iRap;
	  break;
	}
      }
      rapForPTIndex = -1;
      for(int iRap = 0; iRap < jpsi::kNbRapForPTBins; iRap++){
	if(TMath::Abs(oniaSim_rap) > jpsi::rapForPTRange[iRap] && 
	   TMath::Abs(oniaSim_rap) < jpsi::rapForPTRange[iRap+1]){
	  rapForPTIndex = iRap+1;
	  break;
	}
      }
      pTIndex = -1;
      for(int iPT = 0; iPT < jpsi::kNbPTBins[rapForPTIndex]; iPT++){
	if(oniaSim_pt > jpsi::pTRange[rapForPTIndex][iPT] && oniaSim_pt < jpsi::pTRange[rapForPTIndex][iPT+1]){
	  pTIndex = iPT+1;
	  break;
	}
      }
      rapIntegratedPTIndex = -1;
      for(int iPT = 0; iPT < jpsi::kNbPTBins[0]; iPT++){
	if(oniaSim_pt > jpsi::pTRange[0][iPT] && oniaSim_pt < jpsi::pTRange[0][iPT+1]){
	  rapIntegratedPTIndex = iPT+1;
	  break;
	}
      }
      if(rapIndex < 0){
	// printf("rapIndex %d, rap(oniaSim) = %f\n", rapIndex, oniaSim_rap);
	continue;
      }
      if(rapForPTIndex < 1){
	// printf("rapForPTIndex %d, rap(oniaSim) = %f\n", rapForPTIndex, oniaSim_rap);
	continue;
      }
      if(pTIndex < 1){
	// printf("pTIndex %d, pT(oniaSim) = %f\n", pTIndex, oniaSim_pt);
	continue;
      }

      hGenCutPtRap->Fill(oniaSim_rap, oniaSim_pt, kinWeight);

      hMass_Smeared[pTIndex][rapForPTIndex]->Fill(oniaSim_mass, kinWeight);
      hMass_Smeared[0][rapForPTIndex]->Fill(oniaSim_mass, kinWeight);
      hMass_Smeared[pTIndex][0]->Fill(oniaSim_mass, kinWeight);

      //select events within a narrow mass window around the J/psi
      //(rapidity dependence of the resolution --> different mass windows)
      //values obtained from Gaussian fits in "plotMass.C"
      jPsiMassMin = jpsi::polMassJpsi[rapForPTIndex] - jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
      jPsiMassMax = jpsi::polMassJpsi[rapForPTIndex] + jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
      if(oniaSim_mass < jPsiMassMin || oniaSim_mass > jPsiMassMax)
	continue;

      //=====================================
      calcPol(*muPosSmeared, *muNegSmeared);
      //=====================================

      for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){

	Double_t weight = CalcPolWeight(thisCosTh[iFrame]);
	weight *= kinWeight;

	//histos for neg. and pos. rapidity separately:
	if(rapIndex >= 0)
	  hGenCut_pT_rapNP[iFrame][0][rapIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(pTIndex > 0 && rapIndex >= 0)
	  hGenCut_pT_rapNP[iFrame][pTIndex][rapIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);

	//histos taking together +y and -y
	hGenCut_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapIntegratedPTIndex > 0)
	  hGenCut_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapForPTIndex > 0)
	  hGenCut_pT_rap[iFrame][0][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(pTIndex > 0 && rapForPTIndex > 0)
	  hGenCut_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);

	//histos taking together +y and -y and phi-4-folding
	Double_t phiFolded = thisPhi[iFrame];
	Double_t thetaAdjusted = thisCosTh[iFrame];
	if(thisPhi[iFrame] > -90. && thisPhi[iFrame] < 0.)
	  phiFolded *= -1;
	else if(thisPhi[iFrame] > 90. && thisPhi[iFrame] < 180.){
	  phiFolded = 180. - thisPhi[iFrame];
	  thetaAdjusted *= -1;
	}
	else if(thisPhi[iFrame] > -180. && thisPhi[iFrame] < -90.){
	  phiFolded = 180. + thisPhi[iFrame];
	  thetaAdjusted *= -1;
	}
	hGenCut_pT_rap_phiFolded[iFrame][0][0]->Fill(thetaAdjusted, phiFolded, weight);
	if(rapIntegratedPTIndex > 0)
	  hGenCut_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex][0]->Fill(thetaAdjusted, phiFolded, weight);
	if(rapForPTIndex > 0)
	  hGenCut_pT_rap_phiFolded[iFrame][0][rapForPTIndex]->Fill(thetaAdjusted, phiFolded, weight);
	if(pTIndex > 0 && rapForPTIndex > 0)
	  hGenCut_pT_rap_phiFolded[iFrame][pTIndex][rapForPTIndex]->Fill(thetaAdjusted, phiFolded, weight);

// 	//debugging:
// 	if(rapForPTIndex == 1 && pTIndex == 4 && iFrame == jpsi::CS){
// 	  if(thetaAdjusted > 0.6 && thetaAdjusted < 0.9 && phiFolded < 50.){
// 	    printf("pT %1.2f, y %1.2f, M %1.3f, cosThetaCS %1.3f, phiCS %1.3f --> posMu: pT %1.3f, eta %1.3f, phi %1.3f, negMu: pT %1.3f, eta %1.3f, phi %1.3f\n",
// 		   oniaSim_pt, oniaSim_rap, oniaSim_mass, thetaAdjusted, phiFolded,
// 		   muPosSmeared->Pt(), muPosSmeared->Eta(), muPosSmeared->Phi(),
// 		   muNegSmeared->Pt(), muNegSmeared->Eta(), muNegSmeared->Phi());
// 	  }
// 	}
      }

      //delete the created LorentzVectors:
      delete oniaGen;
      delete muPos;
      delete muNeg; 
      delete muPosSmeared;
      delete muNegSmeared; 
      delete oniaFSR;
      delete oniaSim;
   }
   delete fRap;
   delete fPT;
}

//==========================================================
Double_t GeomAcc::CalcPolWeight(Double_t cosTh){

//   Double_t lambdaTh = -1.;
  Double_t lambdaTh = 0.;
  Double_t weight = 1. + lambdaTh * TMath::Power(cosTh, 2.);
  return weight;
}

//==========================================================
TLorentzVector* GeomAcc::ApplySmearing(TLorentzVector *vec, const double *parval){

  TLorentzVector *vecSmeared;
  Double_t pT = vec->Pt();
  Double_t eta = vec->Eta();
  Double_t phi = vec->Phi();
  Double_t pTot = vec->P();
  Double_t theta = vec->Theta();
  Double_t pTRes = sigmaPt(pT,  eta, parval);
  Double_t pTSmeared = gRandom->Gaus(pT, pTRes);
//         printf("pT %1.3f, eta %1.3f --> pT resolution: %1.3f --> smeared pT: %1.3f\n", 
//   	     pT, eta, pTRes, pTSmeared);

  Double_t resCotgTh = sigmaCotgTh(pT,  eta, parval);
  Double_t CotgThSmeared = gRandom->Gaus(1./TMath::Tan(theta), resCotgTh);
  Double_t ThetaSmeared = TMath::ATan(1./CotgThSmeared);
  Double_t etaSmeared = -TMath::Log(TMath::Tan(fabs(ThetaSmeared)/2.));
  if(eta < 0. && etaSmeared > 0.) etaSmeared *= -1.;
//   printf("pT %1.3f, eta %1.3f --> res in CotgTh %1.3f --> cotgThSmeared: %1.3f --> theta %1.3f --> thetaSmeared %1.3f --> Eta smeared of %1.3f\n",
// 	 pT, eta, resCotgTh, CotgThSmeared, theta, ThetaSmeared, etaSmeared);


//   printf("pT %1.3f GeV/c, eta = %1.3f: pT changed by %1.3f %%, eta by %1.3f%%\n", 
// 	 pT, eta, 100.* pTSmeared / pT, 100.*etaSmeared / eta);

	 
//   Double_t pLSmeared = sqrt(pTot*pTot - pow(pTSmeared,2));
//   Double_t thetaSmeared = TMath::ATan(pTSmeared / pLSmeared);
//   Double_t etaSmeared = -TMath::Log(TMath::Tan(thetaSmeared / 2.));
//   if(eta < 0. && etaSmeared > 0.)
//     etaSmeared *= -1.;
// //   printf("eta %1.3f --> smeared eta %1.3f\n", eta, etaSmeared);

  vecSmeared = new TLorentzVector();
  vecSmeared->SetPtEtaPhiM(pTSmeared, etaSmeared, phi, jpsi::muMass);
  return vecSmeared;
}

//==========================================================
//taken from 
//http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/MuonAnalysis/MomentumScaleCalibration/interface/Functions.h?revision=1.52&view=markup
//according to the mail of Marco De Mattia on 02/12/2010
double GeomAcc::sigmaPt(const double &pt, const double &eta, const double *parval){

  double fabsEta = fabs(eta);
  
  double ptPart = parval[13]*pt;
  if( fabsEta > 2.0 ) {
    ptPart = parval[22]*pt + parval[23]*pt*pt;
  }
  else if( fabsEta > 1.4 ) {
    ptPart = parval[20]*pt + parval[21]*pt*pt;
  }
  if(fabsEta<parval[0]) {
    return( ptPart + parval[1] + parval[2]*fabsEta + parval[3]*eta*eta );
  }
  // Return a line connecting the two parabolas
  else if( fabsEta < parval[14] ) {
    double x_1 = parval[0];
    double y_1 = parval[1] + parval[2]*parval[0] + 
      parval[3]*parval[0]*parval[0];
    double x_2 = parval[14];
    double y_2 = parval[4] + parval[5]*fabs((parval[14]-parval[7])) + 
      parval[6]*(parval[14]-parval[7])*(parval[14]-parval[7]);
    return( (fabsEta - x_1)*(y_2 - y_1)/(x_2 - x_1) + y_1 );
  }
  else if( fabsEta < parval[15] ) {
    return( ptPart + parval[4] + parval[5]*fabs(fabsEta-parval[7]) + 
	    parval[6]*(fabsEta-parval[7])*(fabsEta-parval[7]) );
  }
  else {
    return( ptPart + parval[16] + parval[17]*fabs(fabsEta-parval[19]) + 
	    parval[18]*(fabsEta-parval[19])*(fabsEta-parval[19]) );
  }
}

//===================================================
double GeomAcc::sigmaCotgTh(const double & pt, const double & eta, const double *parval) {

  double fabsEta = fabs(eta);
  double value = parval[8] + parval[9]*fabsEta + parval[10]*eta*eta + parval[11]*fabsEta*fabsEta*fabsEta;
  if( value > 0 ) {
    return( value );
  }
  return 0;
}

//==========================================================
Bool_t GeomAcc::areMuonsInAcceptance(Int_t iCut, Double_t pTHigh, Double_t etaHigh, Double_t pTLow, Double_t etaLow){

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
