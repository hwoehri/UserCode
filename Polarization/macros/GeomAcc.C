#define GeomAcc_cxx
#include "GeomAcc.h"
#include "calcPol.C"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"

using namespace std;

Double_t massMuOnia;

TH1D *hMass_Smeared[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2D *hGen_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TH2D *hGenCut_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];

Double_t CalcPolWeight(Double_t thisCosTh);
//================================
void GeomAcc::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     if(jentry%10000 == 0) printf("event %d\n", (Int_t) jentry);
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //get the *generated* onia
      Double_t enOnia = sqrt(onia_Px*onia_Px + onia_Py*onia_Py + onia_Pz*onia_Pz + massMuOnia*massMuOnia);
      TLorentzVector *onia4mom = new TLorentzVector();
      onia4mom->SetPxPyPzE(onia_Px, onia_Py, onia_Pz, enOnia);
      Double_t rapOnia = onia4mom->Rapidity();
      Double_t pTOnia = onia4mom->Pt();

      //get the individual muons that can be affected by FSR
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

      //=====================
      calcPol(*muPos, *muNeg);
      //=====================

      Int_t rapIndex = -1;
      for(int iRap = 0; iRap < 2*jpsi::kNbRapBins; iRap++){
	if(rapOnia > jpsi::rapRange[iRap] && rapOnia < jpsi::rapRange[iRap+1]){
	  rapIndex = iRap;
	  break;
	}
      }
      Int_t rapForPTIndex = -1;
      for(int iRap = 0; iRap < jpsi::kNbRapForPTBins; iRap++){
	if(TMath::Abs(rapOnia) > jpsi::rapForPTRange[iRap] && 
	   TMath::Abs(rapOnia) < jpsi::rapForPTRange[iRap+1]){
	  rapForPTIndex = iRap+1;
	  break;
	}
      }
      Int_t pTIndex = -1;
      for(int iPT = 0; iPT < jpsi::kNbPTBins[rapForPTIndex]; iPT++){
	if(pTOnia > jpsi::pTRange[rapForPTIndex][iPT] && pTOnia < jpsi::pTRange[rapForPTIndex][iPT+1]){
	  pTIndex = iPT+1;
	  break;
	}
      }
      Int_t rapIntegratedPTIndex = -1;
      for(int iPT = 0; iPT < jpsi::kNbPTBins[0]; iPT++){
	if(pTOnia > jpsi::pTRange[0][iPT] && pTOnia < jpsi::pTRange[0][iPT+1]){
	  rapIntegratedPTIndex = iPT+1;
	  break;
	}
      }

      for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){

	Double_t weight = CalcPolWeight(thisCosTh[iFrame]);

	hGen_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapIntegratedPTIndex > 0)
	  hGen_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapForPTIndex > 0)
	  hGen_pT_rap[iFrame][0][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(pTIndex > 0 && rapForPTIndex > 0)
	  hGen_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      }

      //now, smear the two muons with the detector resolution
      //H: --> to be implemented
      //==================================================

      //apply the kinematical phase space cut:
      //(a) on the positive muon
      if((fabs(etaMuPos) < jpsi::etaPS[0] && pTMuPos < jpsi::pTMuMin[0]) || //mid-rapidity cut
	 (fabs(etaMuPos) > jpsi::etaPS[0] && fabs(etaMuPos) < jpsi::etaPS[1] && pMuPos < jpsi::pMuMin[1]) ||
	 (fabs(etaMuPos) > jpsi::etaPS[1] && fabs(etaMuPos) < jpsi::etaPS[2] && pTMuPos < jpsi::pTMuMin[1]))
	continue;
      //(b) on the negative muon
      if((fabs(etaMuNeg) < jpsi::etaPS[0] && pTMuNeg < jpsi::pTMuMin[0]) || //mid-rapidity cut
	 (fabs(etaMuNeg) > jpsi::etaPS[0] && fabs(etaMuNeg) < jpsi::etaPS[1] && pMuNeg < jpsi::pMuMin[1]) ||
	 (fabs(etaMuNeg) > jpsi::etaPS[1] && fabs(etaMuNeg) < jpsi::etaPS[2] && pTMuNeg < jpsi::pTMuMin[1]))
	continue;

      //build the invariant mass, pt, ... of the two muons
      TLorentzVector *oniaSim = new TLorentzVector();
      *oniaSim = *(muPos) + *(muNeg);

      Double_t oniaSim_mass = oniaSim->M();
      Double_t oniaSim_pt = oniaSim->Pt();
      Double_t oniaSim_P = oniaSim->P();
      Double_t oniaSim_eta = oniaSim->PseudoRapidity();
      Double_t oniaSim_rap = oniaSim->Rapidity();
      Double_t oniaSim_phi = oniaSim->Phi();
      Double_t oniaSim_mT = sqrt(oniaSim_mass*oniaSim_mass + oniaSim_pt*oniaSim_pt);

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

      hMass_Smeared[pTIndex][rapForPTIndex]->Fill(oniaSim_mass);
      hMass_Smeared[0][rapForPTIndex]->Fill(oniaSim_mass);
      hMass_Smeared[pTIndex][0]->Fill(oniaSim_mass);
      hMass_Smeared[0][0]->Fill(oniaSim_mass);

      //select events within a narrow mass window around the J/psi
      //(rapidity dependence of the resolution --> different mass windows)
      //values obtained from Gaussian fits in "plotMass.C"
      Double_t jPsiMassMin = jpsi::polMassJpsi[rapForPTIndex] - jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
      Double_t jPsiMassMax = jpsi::polMassJpsi[rapForPTIndex] + jpsi::nSigMass*jpsi::sigmaMassJpsi[rapForPTIndex];
      if(oniaSim_mass < jPsiMassMin || oniaSim_mass > jPsiMassMax)
	continue;

      //=====================
      calcPol(*muPos, *muNeg);
      //=====================
      for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){

	Double_t weight = CalcPolWeight(thisCosTh[iFrame]);

	hGenCut_pT_rap[iFrame][0][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapIntegratedPTIndex > 0)
	  hGenCut_pT_rap[iFrame][rapIntegratedPTIndex][0]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(rapForPTIndex > 0)
	  hGenCut_pT_rap[iFrame][0][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
	if(pTIndex > 0 && rapForPTIndex > 0)
	  hGenCut_pT_rap[iFrame][pTIndex][rapForPTIndex]->Fill(thisCosTh[iFrame], thisPhi[iFrame], weight);
      }
   }

}

//==========================================================
Double_t CalcPolWeight(Double_t thisCosTh){

//   Double_t lambdaTh = -1.;
  Double_t lambdaTh = 0.;
  Double_t weight = 1. + lambdaTh * TMath::Power(thisCosTh, 2.);
  return weight;
}
