#define RecoEff_cxx
#include "RecoEff.h"
#include "calcPol.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>
#include "TRandom.h"

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
Double_t massMu = 0.105658367;
Double_t CalcPolWeight(Double_t thisCosTh);
//==============================================
void RecoEff::Loop(Int_t selDimuType, Bool_t useMCWeight)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t cutAtRecEvent = nentries;
  Long64_t countRecEvent = 0, countTrigEvent = 0;
  Long64_t nb = 0;
  printf("<TrigEff> number of entries = %d\n", (Int_t) nentries);

  //apply a weight as the ratio of pT distributions:
  //1.) default Pythia (parametrization)
  TF1 *fPT = new TF1("fPT", "[0]*x*pow(1.+(1./([1]-2.))*x*x/[2],-[1])", 0., 50.);
  fPT->FixParameter(0, 0.61554); fPT->SetParName(0, "norm");
  fPT->FixParameter(1, 4.05847); fPT->SetParName(1, "beta");
  fPT->FixParameter(2, 19.93); fPT->SetParName(2, "<pT2> [GeV2]");

  //alternatively, the pT distribution can be taken from the
  //a CASCADE calculation:
  TF1 *fPTCASCADE = new TF1("fPTCASCADE", "[0]*x*pow(1.+(1./([1]-2.))*x*x/[2],-[1])", 0., 50.);
  fPTCASCADE->FixParameter(0, 1.); fPTCASCADE->SetParName(0, "norm");
  fPTCASCADE->FixParameter(1, 3.87004); fPTCASCADE->SetParName(1, "beta");
  fPTCASCADE->FixParameter(2, 15.5821); fPTCASCADE->SetParName(2, "<pT2> [GeV2]");
    
  //alternatively, the pT distribution can be taken from the
  //preliminary ATLAS data shown in Moriond, 17th April
  TF1 *fPTATLAS = new TF1("fPTATLAS", "[0]*x*pow(1.+(1./([1]-2.))*x*x/[2],-[1])", 0., 50.);
  fPTATLAS->FixParameter(0, 1.); fPTATLAS->SetParName(0, "norm");
  fPTATLAS->FixParameter(1, 3.87004); fPTATLAS->SetParName(1, "beta");
  fPTATLAS->FixParameter(2, 15.5821); fPTATLAS->SetParName(2, "<pT2> [GeV2]");


  //loop over the events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if(jentry % 100000 == 0) printf("<TrigEff> event %d\n", (Int_t) jentry);

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    //continue only if there is really a RECO dimuon in the event:
    if(onia->Pt() > 990.) 
      continue;

    //===================================================================
    //1.) fill the reconstructed variables in our phase space
    //===================================================================

    Double_t etaMuPos = muPos->PseudoRapidity();
    Double_t etaMuNeg = muNeg->PseudoRapidity();
    Double_t pTMuPos = muPos->Pt();
    Double_t pTMuNeg = muNeg->Pt();
    Double_t pMuPos = muPos->P();
    Double_t pMuNeg = muNeg->P();
    Double_t phiMuPos = muPos->Phi();
    Double_t phiMuNeg = muNeg->Phi();
    Double_t deltaPhi = phiMuPos - phiMuNeg;

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

    Double_t onia_mass = onia->M();
    Double_t onia_pt = onia->Pt();
    Double_t onia_P = onia->P();
    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();
    Double_t onia_phi = onia->Phi();
    Double_t onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);

    Int_t rapIndex = -1;
    for(int iRap = 0; iRap < 2*jpsi::kNbRapBins; iRap++){
      if(onia_rap > jpsi::rapRange[iRap] && onia_rap < jpsi::rapRange[iRap+1]){
    	rapIndex = iRap;
    	break;
      }
    }
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

    //reject processing of events where the dimuon type (GG, GT or TT)
    //does not correspond to the chosen one
    if(selDimuType < 3 && JpsiType != selDimuType)
      continue;
    else if(selDimuType == 3 && JpsiType > 1) //only GG or GT
      continue;

    if(TMath::Abs(onia_rap) > jpsi::rapYPS)
      continue;

    //select events within a narrow mass window around the J/psi
    //(rapidity dependence of the resolution --> different mass windows)
    //values obtained from Gaussian fits in "plotMass.C"
    Double_t jPsiMassMin = jpsi::polMassJpsi[rapForPTIndex] - jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
    Double_t jPsiMassMax = jpsi::polMassJpsi[rapForPTIndex] + jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
    if(onia_mass < jPsiMassMin || onia_mass > jPsiMassMax)
      continue;

    countRecEvent++;

    //=====================
    calcPol(*muPos, *muNeg);
    //=====================

    Double_t weight = 1.;
    if(useMCWeight)
      weight *= MCweight;

//     weight *= fPTCASCADE->Eval(onia_pt) / fPT->Eval(onia_pt);
    weight *= fPTATLAS->Eval(onia_pt) / fPT->Eval(onia_pt);

    //fill the  histograms after reconstruction
    for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){

      //weight = CalcPolWeight(thisCosTh[iFrame]);

      //histos for neg. and pos. rapidity separately:
      if(rapIndex >= 0)
	hGen2D_pol_pT_rapNP[iFrame][0][rapIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      if(pTIndex > 0 && rapIndex >= 0)
	hGen2D_pol_pT_rapNP[iFrame][pTIndex][rapIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	
      //histos taking together +y and -y
      hGen2D_pol_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      if(rapIntegratedPTIndex > 0)
	hGen2D_pol_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      if(rapForPTIndex > 0)
	hGen2D_pol_pT_rap[iFrame][0][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      if(pTIndex > 0 && rapForPTIndex > 0)
	hGen2D_pol_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	  
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
      if(rapIntegratedPTIndex > 0)
	hGen2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex][0]->Fill(thetaAdjusted, phiFolded, weight);
      if(rapForPTIndex > 0)
	hGen2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex]->Fill(thetaAdjusted, phiFolded, weight);
      if(pTIndex > 0 && rapForPTIndex > 0)
	hGen2D_pol_pT_rap_phiFolded[iFrame][pTIndex][rapForPTIndex]->Fill(thetaAdjusted, phiFolded, weight);
    }


    //apply the trigger flag (Fall10 samples have only the HLT_Mu0_TkMu0_OST_Jpsi implemented):
//     if(HLT_Mu0_TkMu0_OST_Jpsi != 1) continue;
    if(HLT_DoubleMu0 != 1) continue;

    //=====================
    calcPol(*muPos, *muNeg);
    //=====================

    //remaining of the events will be used for the analysis
    countTrigEvent++;

    //fill mass, phi, pt, eta and rap distributions
    //a) all bins
    Reco_Onia_mass[0][0]->Fill(onia_mass);
    Reco_Onia_mass[rapIntegratedPTIndex][0]->Fill(onia_mass);
    Reco_Onia_mass[0][rapForPTIndex]->Fill(onia_mass);
    //
    if(onia_mass > (jpsi::polMassJpsi[rapForPTIndex] - 2.*jpsi::sigmaMassJpsi[rapForPTIndex]) &&
       onia_mass < (jpsi::polMassJpsi[rapForPTIndex] + 2.*jpsi::sigmaMassJpsi[rapForPTIndex])){
	 Reco_Onia_phi[0][0]->Fill(onia_phi);
	 Reco_Onia_phi[rapIntegratedPTIndex][0]->Fill(onia_phi);
	 Reco_Onia_phi[0][rapForPTIndex]->Fill(onia_phi);
	 //      
	 Reco_Onia_pt[0]->Fill(onia_pt);
    }
    //b) individual pT and rap bins:
    Reco_Onia_mass[pTIndex][rapForPTIndex]->Fill(onia_mass);
    if(onia_mass > (jpsi::polMassJpsi[rapForPTIndex] - 2.*jpsi::sigmaMassJpsi[rapForPTIndex]) &&
       onia_mass < (jpsi::polMassJpsi[rapForPTIndex] + 2.*jpsi::sigmaMassJpsi[rapForPTIndex])){

      Reco_Onia_phi[pTIndex][rapForPTIndex]->Fill(onia_phi);
      Reco_Onia_pt[rapForPTIndex]->Fill(onia_pt);

      Reco_Onia_rap_pT->Fill(onia_rap, onia_pt);
    }

//     //===================================================
//     //calculate delta, the angle between the CS and HX frame
//     //Formula from: EJP C69 (2010) 657
//     Double_t deltaHXToCS = TMath::ACos(onia_mass * onia->Pz() / (onia_mT * onia_P));
//     //     Double_t deltaCSToHX = -deltaHXToCS;
//     Double_t sin2Delta = pow((onia_pt * onia->Energy() / (onia_P * onia_mT)),2);
//     //sin2Delta does not change sign when going from HX-->CS or vice versa
//     hDelta[pTIndex][rapForPTIndex]->Fill(deltaHXToCS * 180./TMath::Pi());
//     hSin2Delta[pTIndex][rapForPTIndex]->Fill(sin2Delta);
//     //===================================================

//     Double_t deltaPhi = muPos->Phi() - muNeg->Phi();
//     if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
//     else if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();

//     //debugging histos
//     hPhiPos_PhiNeg[pTIndex][rapForPTIndex]->Fill(180./TMath::Pi() * muNeg->Phi(), 180./TMath::Pi() * muPos->Phi());
//     hPtPos_PtNeg[pTIndex][rapForPTIndex]->Fill(muNeg->Pt(), muPos->Pt());
//     hEtaPos_EtaNeg[pTIndex][rapForPTIndex]->Fill(muNeg->PseudoRapidity(), muPos->PseudoRapidity());
//     hDeltaPhi[pTIndex][rapForPTIndex]->Fill(deltaPhi);

//     Reco_mupl_pt[pTIndex][rapForPTIndex]->Fill(muPos->Pt());
//     Reco_mupl_eta[pTIndex][rapForPTIndex]->Fill(muPos->PseudoRapidity());
//     Reco_mupl_phi[pTIndex][rapForPTIndex]->Fill(muPos->Phi());
      
//     Reco_mumi_pt[pTIndex][rapForPTIndex]->Fill(muNeg->Pt());
//     Reco_mumi_eta[pTIndex][rapForPTIndex]->Fill(muNeg->PseudoRapidity());
//     Reco_mumi_phi[pTIndex][rapForPTIndex]->Fill(muNeg->Phi());



    //fill the histos for all the different frames
    for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){

      // 	weight = CalcPolWeight(thisCosTh[iFrame]);

	//histos for neg. and pos. rapidity separately:
	if(rapIndex >= 0)
	  Reco2D_pol_pT_rapNP[iFrame][0][rapIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(pTIndex > 0 && rapIndex >= 0)
	  Reco2D_pol_pT_rapNP[iFrame][pTIndex][rapIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	
	//histos taking together +y and -y
	Reco2D_pol_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapIntegratedPTIndex > 0)
	  Reco2D_pol_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapForPTIndex > 0)
	  Reco2D_pol_pT_rap[iFrame][0][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(pTIndex > 0 && rapForPTIndex > 0)
	  Reco2D_pol_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	  
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
	if(rapIntegratedPTIndex > 0)
	  Reco2D_pol_pT_rap_phiFolded[iFrame][rapIntegratedPTIndex][0]->Fill(thetaAdjusted, phiFolded, weight);
	if(rapForPTIndex > 0)
	  Reco2D_pol_pT_rap_phiFolded[iFrame][0][rapForPTIndex]->Fill(thetaAdjusted, phiFolded, weight);
	if(pTIndex > 0 && rapForPTIndex > 0)
	  Reco2D_pol_pT_rap_phiFolded[iFrame][pTIndex][rapForPTIndex]->Fill(thetaAdjusted, phiFolded, weight);
      }

  }//loop over entries

  printf("nb. of rec. and triggered events are %d and %d, of a total of %d events\n", 
	 (Int_t) countRecEvent, (Int_t) countTrigEvent, (Int_t) nentries);
}

//==========================================================
Double_t CalcPolWeight(Double_t cosTh){

  Double_t lambdaTh = 0.;
  Double_t weight = 1. + lambdaTh * TMath::Power(cosTh, 2.);
  return weight;
}
