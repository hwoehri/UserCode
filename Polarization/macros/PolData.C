#define PolData_cxx
#include "PolData.h"
#include "../interface/commonVar.h"

#include "TLorentzVector.h"
#include "TClonesArray.h"
#include <TH2.h>
#include <TCanvas.h>
//
Double_t thisCosTh_CS, thisPhi_CS, thisPhi_CS_rad, thisCosPhi_CS;
Double_t thisCosTh_HX, thisPhi_HX, thisPhi_HX_rad, thisCosPhi_HX;

//some statistics
TH1F *Reco_StatEv;

//histos at reco level for mu+, mu-, gamma, Onia
TH1F *Reco_mupl_pt[kNbPTBins+1][kNbRapForPTBins+1];
TH1F *Reco_mumi_pt[kNbPTBins+1][kNbRapForPTBins+1];
TH1F *Reco_mupl_eta[kNbPTBins+1][kNbRapForPTBins+1];
TH1F *Reco_mumi_eta[kNbPTBins+1][kNbRapForPTBins+1];
TH1F *Reco_mupl_phi[kNbPTBins+1][kNbRapForPTBins+1];
TH1F *Reco_mumi_phi[kNbPTBins+1][kNbRapForPTBins+1];
TH1F *Reco_Onia_mass[kNbPTBins+1][kNbRapForPTBins+1];
TH1F *Reco_Onia_pt[kNbRapForPTBins+1];
TH1F *Reco_Onia_eta[kNbPTBins+1];
TH1F *Reco_Onia_rap[kNbPTBins+1];
TH1F *Reco_Onia_phi[kNbPTBins+1][kNbRapForPTBins+1];
TH2F *Reco_Onia_rap_pT;
//debugging histos:
TH2F *hPhiPos_PhiNeg[kNbPTBins+1][kNbRapForPTBins+1];
TH2F *hPtPos_PtNeg[kNbPTBins+1][kNbRapForPTBins+1];
TH2F *hEtaPos_EtaNeg[kNbPTBins+1][kNbRapForPTBins+1];
TH1F *hDeltaPhi[kNbPTBins+1][2*kNbRapBins+1];

//polarization histos:
TH1F *Reco_Onia_pol_pT[kNbFrames][kNbPTBins+1][kNbPolVar];
TH1F *Reco_Onia_pol_rap[kNbFrames][2*kNbRapBins+1][kNbPolVar];
TH1F *Reco_Onia_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1][kNbPolVar];
TH2F *Reco2D_Onia_pol_pT[kNbFrames][kNbPTBins+1];
TH2F *Reco2D_Onia_pol_rap[kNbFrames][2*kNbRapBins+1];
TH2F *Reco2D_Onia_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1];

TH1F *hDelta[kNbPTBins+1][kNbRapForPTBins+1];
TH1F *hSin2Delta[kNbPTBins+1][kNbRapForPTBins+1];


void calcPol(TLorentzVector muplus_LAB, TLorentzVector muminus_LAB);
Double_t CalcPolWeight(Double_t pf_onia_P, Double_t thisCosTh_CS);

//==============================================
void PolData::Loop(Int_t selDimuType, Bool_t writeOutEvents)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t cutAtRecEvent = nentries;
  Long64_t countRecEvent = 0;
  Long64_t nb = 0;
  printf("number of entries = %d\n", (Int_t) nentries);
  FILE *fOutputTextFile;
  if(writeOutEvents)
    fOutputTextFile = fopen("jPsiCandidates.txt", "write");

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if(jentry % 100000 == 0) printf("event %d\n", (Int_t) jentry);

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    //fill the RECO variables only
    //--> reject all events that contain only GEN information:
    if(JpsiPx < -9900.)
      continue;

    Reco_StatEv->Fill(0.5);//count all events

    //reject processing of events where the dimuon type (GG, GT or TT)
    //does not correspond to the chosen one
    if(selDimuType < 3 && JpsiType_idx != selDimuType)
      continue;
    else if(selDimuType == 3 && JpsiType_idx > 1) //only GG or GT
      continue;

    Reco_StatEv->Fill(1.5);//count all events

    Double_t muMass = 0.105658;
    Double_t enMuPos = sqrt(muPosPx*muPosPx + muPosPy*muPosPy + muPosPz*muPosPz + muMass*muMass);
    Double_t enMuNeg = sqrt(muNegPx*muNegPx + muNegPy*muNegPy + muNegPz*muNegPz + muMass*muMass);
    TLorentzVector *muPos = new TLorentzVector();
    TLorentzVector *muNeg = new TLorentzVector();
    muPos->SetPxPyPzE(muPosPx, muPosPy, muPosPz, enMuPos);
    muNeg->SetPxPyPzE(muNegPx, muNegPy, muNegPz, enMuNeg);
    Double_t etaMuPos = muPos->PseudoRapidity();
    Double_t etaMuNeg = muNeg->PseudoRapidity();
    Double_t pTMuPos = muPos->Pt();
    Double_t pTMuNeg = muNeg->Pt();

    //test of fiducial area:
    if(TMath::Abs(etaMuPos) > etaPS || TMath::Abs(etaMuNeg) > etaPS){
      // printf("eta(pos. muon) = %f, eta(neg. muon) = %f\n", etaMuPos, etaMuNeg);
      continue;
    }

    Reco_StatEv->Fill(2.5);//count all events

    if(pTMuPos < pTMuMin && pTMuNeg < pTMuMin){
      // printf("pT(pos. muon) = %f, pT(neg. muon) = %f\n", pTMuPos, pTMuNeg);
      continue;
    }

    //test according to Gavin's proposal:
    //if any of the two muons is within 1.4 < eta < 1.6 AND
    //the two muons are close in eta (deltaEta < 0.2)
    //reject the dimuon (no matter whether it is "Seagull" or 
    //"Cowboy"):
//     if(TMath::Abs(etaMuPos - etaMuNeg) < 0.2 &&
//        ((etaMuPos > 1.4 && etaMuPos < 1.6) || (etaMuNeg > 1.4 && etaMuNeg < 1.6))){
//       printf("rejecting the event!\n");
//       continue;
//     }

    Reco_StatEv->Fill(3.5);

    calcPol(*muPos, *muNeg);
    //test:
//     calcPol(*muNeg, *muPos);
//     //H: test:
//     if(jentry%2 == 0)
//       calcPol(*muPos, *muNeg);
//     else
//       calcPol(*muNeg, *muPos);

    thisCosPhi_CS = TMath::Cos(2.*thisPhi_CS_rad);
    thisCosPhi_HX = TMath::Cos(2.*thisPhi_HX_rad);
    
    //build the invariant mass, pt, ... of the two muons
    TLorentzVector *onia = new TLorentzVector();
    *onia = *(muPos) + *(muNeg);

    Double_t onia_mass = onia->M();
    Double_t onia_pt = onia->Pt();
    Double_t onia_P = onia->P();
    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();
    Double_t onia_phi = onia->Phi();
    Double_t onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);

//     Double_t weight = CalcPolWeight(onia_P, thisCosTh_CS);
    Double_t weight = CalcPolWeight(onia_P, thisCosTh_HX);

    Int_t pTIndex = -1;
    for(int iPT = 0; iPT < kNbPTBins; iPT++){
      if(onia_pt > pTRange[iPT] && onia_pt < pTRange[iPT+1]){
	pTIndex = iPT+1;
	break;
      }
    }
    Int_t rapIndex = -1;
    for(int iRap = 0; iRap < 2*kNbRapBins; iRap++){
      if(onia_rap > rapRange[iRap] && onia_rap < rapRange[iRap+1]){
	rapIndex = iRap+1;
	break;
      }
    }
    Int_t rapForPTIndex = -1;
    for(int iRap = 0; iRap < kNbRapForPTBins; iRap++){
      if(TMath::Abs(onia_rap) > rapForPTRange[iRap] && 
	 TMath::Abs(onia_rap) < rapForPTRange[iRap+1]){
	 rapForPTIndex = iRap+1;
	break;
      }
    }

    if(pTIndex < 1){
      printf("pTIndex %d, pT(onia) = %f\n", pTIndex, onia_pt);
      continue;
    }
    if(rapIndex < 1){
      printf("rapIndex %d, rap(onia) = %f\n", rapIndex, onia_rap);
      continue;
    }
    if(rapForPTIndex < 1){
      printf("rapForPTIndex %d, rap(onia) = %f\n", rapForPTIndex, onia_rap);
      continue;
    }

    Reco_StatEv->Fill(4.5);

    //select events with a cut on the lifetime to reject NP J/psis:
    if(Jpsict > JpsiCtauMax)
      continue;

    Reco_StatEv->Fill(5.5);

    //select events within a narrow mass window around the J/psi
    //(rapidity dependence of the resolution --> different mass windows)
    //values obtained from Gaussian fits in "plotMass.C"
    if(JpsiMass < JpsiMassMin[rapForPTIndex-1] || 
       JpsiMass > JpsiMassMax[rapForPTIndex-1])
      continue;

    Reco_StatEv->Fill(6.5);

    countRecEvent++;

    if(writeOutEvents)
      fprintf(fOutputTextFile, "%d\t%d\t%d\n", (Int_t) runNb, (Int_t) lumiBlock, (Int_t)  eventNb);

    //===================================================
    //calculate delta, the angle between the CS and HX frame
    //Formula from EPJC paper
    Double_t deltaHXToCS = TMath::ACos(onia_mass * onia->Pz() / (onia_mT * onia_P));
//     Double_t deltaCSToHX = -deltaHXToCS;
    Double_t sin2Delta = pow((onia_pt * onia->Energy() / (onia_P * onia_mT)),2);
    //sin2Delta does not change sign when going from HX-->CS or vice versa
    hDelta[pTIndex][rapForPTIndex]->Fill(deltaHXToCS * 180./TMath::Pi());
    hSin2Delta[pTIndex][rapForPTIndex]->Fill(sin2Delta);
    //===================================================


    Double_t deltaPhi = muPos->Phi() - muNeg->Phi();
    if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
    else if(deltaPhi > TMath::Pi()) deltaPhi = 2.*TMath::Pi() - deltaPhi;

    //debugging histos
    hPhiPos_PhiNeg[pTIndex][rapForPTIndex]->Fill(180./TMath::Pi() * muNeg->Phi(), 180./TMath::Pi() * muPos->Phi());
    hPtPos_PtNeg[pTIndex][rapForPTIndex]->Fill(muNeg->Pt(), muPos->Pt());
    hEtaPos_EtaNeg[pTIndex][rapForPTIndex]->Fill(muNeg->PseudoRapidity(), muPos->PseudoRapidity());
    hDeltaPhi[pTIndex][rapIndex]->Fill(deltaPhi);

    Reco_mupl_pt[pTIndex][rapForPTIndex]->Fill(muPos->Pt());
    Reco_mupl_eta[pTIndex][rapForPTIndex]->Fill(muPos->PseudoRapidity());
    Reco_mupl_phi[pTIndex][rapForPTIndex]->Fill(muPos->Phi());
      
    Reco_mumi_pt[pTIndex][rapForPTIndex]->Fill(muNeg->Pt());
    Reco_mumi_eta[pTIndex][rapForPTIndex]->Fill(muNeg->PseudoRapidity());
    Reco_mumi_phi[pTIndex][rapForPTIndex]->Fill(muNeg->Phi());
	
    if(TMath::Abs(onia_rap) < rapYPS){

      //all bins
      Reco_Onia_mass[0][0]->Fill(onia_mass, weight);
      Reco_Onia_phi[0][0]->Fill(onia_phi, weight);
      Reco_Onia_mass[pTIndex][0]->Fill(onia_mass, weight);
      Reco_Onia_phi[pTIndex][0]->Fill(onia_phi, weight);
      Reco_Onia_mass[0][rapForPTIndex]->Fill(onia_mass, weight);
      Reco_Onia_phi[0][rapForPTIndex]->Fill(onia_phi, weight);
      
      Reco_Onia_pt[0]->Fill(onia_pt, weight);
      Reco_Onia_eta[0]->Fill(onia_eta, weight);
      Reco_Onia_rap[0]->Fill(onia_rap, weight);
      //individual pT and rap bins:
      Reco_Onia_mass[pTIndex][rapForPTIndex]->Fill(onia_mass, weight);
      Reco_Onia_phi[pTIndex][rapForPTIndex]->Fill(onia_phi, weight);
      Reco_Onia_pt[rapForPTIndex]->Fill(onia_pt, weight);
      Reco_Onia_eta[pTIndex]->Fill(onia_eta, weight);
      Reco_Onia_rap[pTIndex]->Fill(onia_rap, weight);

      Reco_Onia_rap_pT->Fill(onia_rap, onia_pt, weight);

      //1a) polariztion histos - all pT
      //CS frame
      Reco_Onia_pol_pT[CS][0][cosThPol]->Fill(thisCosTh_CS, weight);
      Reco_Onia_pol_pT[CS][0][phiPol]->Fill(thisPhi_CS, weight);
      Reco_Onia_pol_pT[CS][0][cos2PhiPol]->Fill(thisCosPhi_CS, weight);
      Reco2D_Onia_pol_pT[CS][0]->Fill(thisCosTh_CS, thisPhi_CS, weight);
      //HX frame
      Reco_Onia_pol_pT[HX][0][cosThPol]->Fill(thisCosTh_HX, weight);
      Reco_Onia_pol_pT[HX][0][phiPol]->Fill(thisPhi_HX, weight);
      Reco_Onia_pol_pT[HX][0][cos2PhiPol]->Fill(thisCosPhi_HX, weight);
      Reco2D_Onia_pol_pT[HX][0]->Fill(thisCosTh_HX, thisPhi_HX, weight);

      //1b) polariztion histos - pT Bin
      if(pTIndex > 0){
	//CS frame
	Reco_Onia_pol_pT[CS][pTIndex][cosThPol]->Fill(thisCosTh_CS, weight);
	Reco_Onia_pol_pT[CS][pTIndex][phiPol]->Fill(thisPhi_CS, weight);
	Reco_Onia_pol_pT[CS][pTIndex][cos2PhiPol]->Fill(thisCosPhi_CS, weight);
	Reco2D_Onia_pol_pT[CS][pTIndex]->Fill(thisCosTh_CS, thisPhi_CS, weight);
	//HX frame
	Reco_Onia_pol_pT[HX][pTIndex][cosThPol]->Fill(thisCosTh_HX, weight);
	Reco_Onia_pol_pT[HX][pTIndex][phiPol]->Fill(thisPhi_HX, weight);
	Reco_Onia_pol_pT[HX][pTIndex][cos2PhiPol]->Fill(thisCosPhi_HX, weight);
	Reco2D_Onia_pol_pT[HX][pTIndex]->Fill(thisCosTh_HX, thisPhi_HX, weight);
      }

      //2a) polariztion histos - all Rap
      //CS frame
      Reco_Onia_pol_rap[CS][0][cosThPol]->Fill(thisCosTh_CS, weight);
      Reco_Onia_pol_rap[CS][0][phiPol]->Fill(thisPhi_CS, weight);
      Reco_Onia_pol_rap[CS][0][cos2PhiPol]->Fill(thisCosPhi_CS, weight);
      Reco2D_Onia_pol_rap[CS][0]->Fill(thisCosTh_CS, thisPhi_CS, weight);
      //HX frame
      Reco_Onia_pol_rap[HX][0][cosThPol]->Fill(thisCosTh_HX, weight);
      Reco_Onia_pol_rap[HX][0][phiPol]->Fill(thisPhi_HX, weight);
      Reco_Onia_pol_rap[HX][0][cos2PhiPol]->Fill(thisCosPhi_HX, weight);
      Reco2D_Onia_pol_rap[HX][0]->Fill(thisCosTh_HX, thisPhi_HX, weight);

      //2b) polariztion histos - rap Bin
      if(rapIndex > 0){
	//CS frame
	Reco_Onia_pol_rap[CS][rapIndex][cosThPol]->Fill(thisCosTh_CS, weight);
	Reco_Onia_pol_rap[CS][rapIndex][phiPol]->Fill(thisPhi_CS, weight);
	Reco_Onia_pol_rap[CS][rapIndex][cos2PhiPol]->Fill(thisCosPhi_CS, weight);
	Reco2D_Onia_pol_rap[CS][rapIndex]->Fill(thisCosTh_CS, thisPhi_CS, weight);
	//HX frame
	Reco_Onia_pol_rap[HX][rapIndex][cosThPol]->Fill(thisCosTh_HX, weight);
	Reco_Onia_pol_rap[HX][rapIndex][phiPol]->Fill(thisPhi_HX, weight);
	Reco_Onia_pol_rap[HX][rapIndex][cos2PhiPol]->Fill(thisCosPhi_HX, weight);
	Reco2D_Onia_pol_rap[HX][rapIndex]->Fill(thisCosTh_HX, thisPhi_HX, weight);
      }

      //3) polariztion histos - pT and rap Bin
      //all pT and rapidities
      Reco2D_Onia_pol_pT_rap[CS][0][0]->Fill(thisCosTh_CS, thisPhi_CS, weight);
      Reco2D_Onia_pol_pT_rap[HX][0][0]->Fill(thisCosTh_HX, thisPhi_HX, weight);
      if(pTIndex > 0){
	Reco2D_Onia_pol_pT_rap[CS][pTIndex][0]->Fill(thisCosTh_CS, thisPhi_CS, weight);
	Reco2D_Onia_pol_pT_rap[HX][pTIndex][0]->Fill(thisCosTh_HX, thisPhi_HX, weight);
      }
      if(rapForPTIndex > 0){
	Reco2D_Onia_pol_pT_rap[CS][0][rapForPTIndex]->Fill(thisCosTh_CS, thisPhi_CS, weight);
	Reco2D_Onia_pol_pT_rap[HX][0][rapForPTIndex]->Fill(thisCosTh_HX, thisPhi_HX, weight);
      }
      if(pTIndex > 0 && rapForPTIndex > 0){
	//CS frame
	Reco_Onia_pol_pT_rap[CS][pTIndex][rapForPTIndex][cosThPol]->Fill(thisCosTh_CS, weight);
	Reco_Onia_pol_pT_rap[CS][pTIndex][rapForPTIndex][phiPol]->Fill(thisPhi_CS, weight);
	Reco_Onia_pol_pT_rap[CS][pTIndex][rapForPTIndex][cos2PhiPol]->Fill(thisCosPhi_CS, weight);
	Reco2D_Onia_pol_pT_rap[CS][pTIndex][rapForPTIndex]->Fill(thisCosTh_CS, thisPhi_CS, weight);
	//HX frame
	Reco_Onia_pol_pT_rap[HX][pTIndex][rapForPTIndex][cosThPol]->Fill(thisCosTh_HX, weight);
	Reco_Onia_pol_pT_rap[HX][pTIndex][rapForPTIndex][phiPol]->Fill(thisPhi_HX, weight);
	Reco_Onia_pol_pT_rap[HX][pTIndex][rapForPTIndex][cos2PhiPol]->Fill(thisCosPhi_HX, weight);
	Reco2D_Onia_pol_pT_rap[HX][pTIndex][rapForPTIndex]->Fill(thisCosTh_HX, thisPhi_HX, weight);
      }
    }

    delete muPos;
    delete muNeg;
    delete onia;
  }//loop over entries

  if(writeOutEvents)
    fclose(fOutputTextFile);

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}

//=========================================
void calcPol(TLorentzVector muplus_LAB, 
	     TLorentzVector muminus_LAB){
  
  TLorentzVector qqbar_LAB = muplus_LAB + muminus_LAB;
  // calculation of decay angular parameters

  // boost beams and positive muon into the q-qbar rest frame:
  TVector3 LAB_to_QQBAR = -qqbar_LAB.BoostVector();

  TLorentzVector beam1_QQBAR = beam1_LAB;
  beam1_QQBAR.Boost( LAB_to_QQBAR );

  TLorentzVector beam2_QQBAR = beam2_LAB;
  beam2_QQBAR.Boost( LAB_to_QQBAR );

  TLorentzVector muplus_QQBAR = muplus_LAB;
  muplus_QQBAR.Boost( LAB_to_QQBAR );

  // reference directions in the Jpsi rest frame:

  TVector3 beam1_direction     = beam1_QQBAR.Vect().Unit();
  TVector3 beam2_direction     = beam2_QQBAR.Vect().Unit();
  TVector3 qqbar_direction     = qqbar_LAB.Vect().Unit();
  TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();

  // all polarization frames have the same Y axis = the normal to the plane formed by
  // the directions of the colliding hadrons
  TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

  /////////////////////////////////////////////////////////////////////
  // CS frame

  TVector3 newZaxis = beam1_beam2_bisect;
  TVector3 newYaxis = Yaxis;
  TVector3 newXaxis = newYaxis.Cross( newZaxis );

  TRotation rotation;
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();   // transforms coordinates from the "xyz" system
  // to the "new" (rotated) system having the polarization axis
  // as z axis

  TVector3 muplus_QQBAR_rotated(muplus_QQBAR.Vect());
  
  muplus_QQBAR_rotated.Transform( rotation );
      
  thisCosTh_CS = muplus_QQBAR_rotated.CosTheta();

  thisPhi_CS_rad = muplus_QQBAR_rotated.Phi();
  thisPhi_CS = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
//   if ( thisPhi_CS < 0. ) thisPhi_CS = 360. + thisPhi_CS;      // phi defined in degrees from 0 to 360
  thisPhi_CS += 180.;

  /////////////////////////////////////////////////////////////////////
  // HELICITY frame

  newZaxis = qqbar_direction;
  newYaxis = Yaxis;
  newXaxis = newYaxis.Cross( newZaxis );

  rotation.SetToIdentity();
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();

  muplus_QQBAR_rotated = muplus_QQBAR.Vect();

  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh_HX = muplus_QQBAR_rotated.CosTheta();

  thisPhi_HX_rad = muplus_QQBAR_rotated.Phi();
  thisPhi_HX = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
//   if ( thisPhi_HX < 0. ) thisPhi_HX = 360. + thisPhi_HX;
  thisPhi_HX += 180.;
}

//==========================================================
Double_t CalcPolWeight(Double_t onia_P, Double_t thisCosTh_CS){

//   Double_t const p0 = 5.;
//   Double_t const kappa = 0.6;
//  //pol. acc. to PRL:
//   Double_t lambdaTh = 1. - TMath::Power(2., 1.-TMath::Power(onia_P/p0, kappa));


  Double_t lambdaTh = -1.;
//   Double_t weight = 1. + lambdaTh * TMath::Power(thisCosTh_CS, 2.);
  Double_t weight = 1.;
  return weight;
}
