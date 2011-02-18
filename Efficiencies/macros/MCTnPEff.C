#define MCTnPEff_cxx
#include "MCTnPEff.h"
#include "calcPol.C"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>

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

Double_t massMuOnia;

Double_t GetTrackingEffPos(Double_t etaMuPos, Double_t pTMuPos);
Double_t GetMuonIDEffPos(Double_t etaMuPos, Double_t pTMuPos);
Double_t GetMuQualEffPos(Double_t etaMuPos, Double_t pTMuPos);
Double_t GetTrackingEffNeg(Double_t etaMuNeg, Double_t pTMuNeg);
Double_t GetMuonIDEffNeg(Double_t etaMuNeg, Double_t pTMuNeg);
Double_t GetMuQualEffNeg(Double_t etaMuNeg, Double_t pTMuNeg);
Double_t GetL1L2TrigEffNeg(Double_t etaMuNeg, Double_t pTMuNeg);
Double_t GetL3TrigEffNeg(Double_t etaMuNeg, Double_t pTMuNeg);
Double_t GetTrackTrigEffNeg(Double_t etaMuNeg, Double_t pTMuNeg);
Double_t GetTkMuTrigEffNeg(Double_t etaMuNeg, Double_t pTMuNeg);
Double_t GetL1L2TrigEffPos(Double_t etaMuPos, Double_t pTMuPos);
Double_t GetL3TrigEffPos(Double_t etaMuPos, Double_t pTMuPos);
Double_t GetTrackTrigEffPos(Double_t etaMuPos, Double_t pTMuPos);
Double_t GetTkMuTrigEffPos(Double_t etaMuPos, Double_t pTMuPos);
//==============================================
void MCTnPEff::Loop(Char_t *trigLabel)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nb = 0;
  Long64_t countRecEvent = 0;
  printf("number of entries = %d\n", (Int_t) nentries);

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

    Double_t epsTrack_Pos  = GetTrackingEffPos(etaMuPos_Gen, pTMuPos_Gen);
    Double_t epsMuonID_Pos = GetMuonIDEffPos(etaMuPos_Gen, pTMuPos_Gen);
    Double_t epsQual_Pos   = GetMuQualEffPos(etaMuPos_Gen, pTMuPos_Gen);
    Double_t epsTrack_Neg  = GetTrackingEffNeg(etaMuNeg_Gen, pTMuNeg_Gen);
    Double_t epsMuonID_Neg = GetMuonIDEffNeg(etaMuNeg_Gen, pTMuNeg_Gen);
    Double_t epsQual_Neg   = GetMuQualEffNeg(etaMuNeg_Gen, pTMuNeg_Gen);
    //check which muon is matched to the HLT-muon:
    Double_t epsL1L2Trig, epsL3Trig, epsTrackTrig, epsTkMuTrig;
    if(1){//dummy requirement; needs fixing
      epsL1L2Trig  = GetL1L2TrigEffPos(etaMuPos_Gen, pTMuPos_Gen);
      epsL3Trig    = GetL3TrigEffPos(etaMuPos_Gen, pTMuPos_Gen);
      epsTrackTrig = GetTrackTrigEffNeg(etaMuNeg_Gen, pTMuNeg_Gen);
      epsTkMuTrig  = GetTkMuTrigEffNeg(etaMuNeg_Gen, pTMuNeg_Gen);
    }
    else{
      epsL1L2Trig = GetL1L2TrigEffNeg(etaMuNeg_Gen, pTMuNeg_Gen);
      epsL3Trig   = GetL3TrigEffNeg(etaMuNeg_Gen, pTMuNeg_Gen);
      epsTrackTrig = GetTrackTrigEffPos(etaMuPos_Gen, pTMuPos_Gen);
      epsTkMuTrig  = GetTkMuTrigEffPos(etaMuPos_Gen, pTMuPos_Gen);
    }


    Double_t recoEff_Pos = epsTrack_Pos * epsMuonID_Pos * epsQual_Pos;
    Double_t recoEff_Neg = epsTrack_Neg * epsMuonID_Neg * epsQual_Neg;
    Double_t recoEff = recoEff_Pos * recoEff_Neg;
    Double_t trigEff = epsL1L2Trig * epsL3Trig * epsTrackTrig * epsTkMuTrig;
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

  }//loop over entries

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}

//==============================================================
Double_t GetTrackingEffPos(Double_t etaMuPos, Double_t pTMuPos){
  return 1.;
}
//==============================================================
Double_t GetMuonIDEffPos(Double_t etaMuPos, Double_t pTMuPos){
  return 1.;
}
//==============================================================
Double_t GetMuQualEffPos(Double_t etaMuPos, Double_t pTMuPos){
  return 1.;
}
//==============================================================
Double_t GetTrackingEffNeg(Double_t etaMuNeg, Double_t pTMuNeg){
  return 1.;
}
//==============================================================
Double_t GetMuonIDEffNeg(Double_t etaMuNeg, Double_t pTMuNeg){
  return 1.;
}
//==============================================================
Double_t GetMuQualEffNeg(Double_t etaMuNeg, Double_t pTMuNeg){
  return 1.;
}
//==============================================================
Double_t GetL1L2TrigEffNeg(Double_t etaMuNeg, Double_t pTMuNeg){
  return 1.;
}
//==============================================================
Double_t GetL3TrigEffNeg(Double_t etaMuNeg, Double_t pTMuNeg){
  return 1.;
}
//==============================================================
Double_t GetTrackTrigEffNeg(Double_t etaMuNeg, Double_t pTMuNeg){
  return 1.;
}
//==============================================================
Double_t GetTkMuTrigEffNeg(Double_t etaMuNeg, Double_t pTMuNeg){
  return 1.;
}
//==============================================================
Double_t GetL1L2TrigEffPos(Double_t etaMuPos, Double_t pTMuPos){
  return 1.;
}
//==============================================================
Double_t GetL3TrigEffPos(Double_t etaMuPos, Double_t pTMuPos){
  return 1.;
}
//==============================================================
Double_t GetTrackTrigEffPos(Double_t etaMuPos, Double_t pTMuPos){
  return 1.;
}
//==============================================================
Double_t GetTkMuTrigEffPos(Double_t etaMuPos, Double_t pTMuPos){
  return 1.;
}



