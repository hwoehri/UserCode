#include "../interface/rootIncludes.inc"
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "isMuonInAcceptance.C"
#include "calcPol.C"

TChain *data;
TLorentzVector *pMu = 0, *pMu_Fixed, *pMu_Gen = 0, *pMu_Gen_Fixed;
TLorentzVector *pMu_Tk = 0, *pMu_Tk_Fixed = 0;
Int_t HLT_Dimuon10_Jpsi_Barrel_v3, HLT_Dimuon10_Jpsi_Barrel_v6;
Int_t HLT_Dimuon0_Jpsi_v3, HLT_Dimuon0_Jpsi_NoVertexing_v3;

Int_t const kNbpT = 17;
Double_t pTBins[kNbpT+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 10., 15., 20., 30., 50.};
// Int_t const kNbpT = 47;
// Double_t pTBins[kNbpT+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 46., 48., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100};
// Int_t const kNbpT = 65;
// Double_t pTBins[kNbpT+1] = {2.0, 2.2, 2.4, 2.6, 2.8,   
// 			    3.0, 3.2, 3.4, 3.6, 3.8,   
// 			    4.0, 4.2, 4.4, 4.6, 4.8,   
// 			    5.0, 5.2, 5.4, 5.6, 5.8,   
// 			    6.0, 6.2, 6.4, 6.6, 6.8,   
// 			    7.0, 7.25, 7.5, 8.0, 9.0, 10., 11., 12., 13., 14., 
// 			    15., 16., 17., 18., 19., 20., 22., 24., 26., 28., 
// 			    30., 32., 34., 36., 38., 40., 42., 44., 46., 48., 
// 			    50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.};
Int_t const kNbEta = 10;
Double_t etaBins[kNbEta+1] = {0, 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.1, 2.4};
// Int_t const kNbEta = 18;
// Double_t etaBins[kNbEta+1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2.1, 2.4};
TEfficiency *recoEff_pT[kNbEta], *recoEff_eta[kNbpT], *recoEff_phi;
TEfficiency *trigEff_pT[kNbEta], *trigEff_eta[kNbpT], *trigEff_phi;
TEfficiency *totEff_pT[kNbEta], *totEff_eta[kNbpT], *totEff_phi;
TEfficiency *recoEff_pT_eta, *trigEff_pT_eta, *totEff_pT_eta;

TGraphAsymmErrors *gtotEff_pT[kNbEta], *gtotEff_eta[kNbpT];

TH1D *yTimesX_pT[kNbEta], *y_pT[kNbEta];
TH1D *yTimesX_eta[kNbpT], *y_eta[kNbpT];

enum {LOOSE,TIGHT};//set of muon fiducial cuts

void SetBranches(Bool_t startFromTrack);
void FillHistos(Bool_t startFromTrack);
void BookHistos();
void WriteHistos();
void CalcWeightedAverage();
//==============================
void singleMuTruthEff(Char_t *fileNameOut = "singleMuTruthEff_25June2012.root",
		      Bool_t startFromTrack = kTRUE){

  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/onia2MuMu_tree.root";
  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed20GeVrap1MuFree_11Jan2012.root";
  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1MuFree_11Jan2012.root";
  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1_MuFree2PT100GeV_12Jan2012.root";
  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1_MuFree2PT100GeV_New_12Jan2012.root";
  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1_MuFree2PT100GeV_New_L1MatchingOnly_12Jan2012.root";
  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1_MuFree2PT100GeV_New_L2Matching_12Jan2012.root";

  data = new TChain("data");
  if(startFromTrack){
    data->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/March2012/TTree_Onia2MuMu_muFixed40GeVrap1_MuFree2PT20GeV_RecoTracks_24March2012.root");
    data->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/March2012/TTree_Onia2MuMu_muFixed40GeVrap1_MuFree2PT100GeV_RecoTracks_24March2012.root");
  }
  else{
    data->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1_MuFree2PT100GeV_New_12Jan2012.root");
    data->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1_MuFree2PT20GeV_New_17Jan2012.root");
  }
  cout << data->GetEntries() << endl;

  SetBranches(startFromTrack);

  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  BookHistos();
  FillHistos(startFromTrack);
  CalcWeightedAverage();
  WriteHistos();
  fOut->Close();
}
//===============================
void FillHistos(Bool_t startFromTrack){

  Long64_t nbEv = data->GetEntries();
  for(Long64_t iEv = 0; iEv < nbEv; iEv++){
  //for(Long64_t iEv = 0; iEv < 600000; iEv++){

    data->GetEntry(iEv);
    if(iEv % 100000 == 0) 
      printf("event %d out of %d\n", (Int_t) iEv, (Int_t) nbEv);

    Double_t ptGen = pMu_Gen->Pt();
    Double_t etaGen = pMu_Gen->Eta();
    Double_t phiGen = pMu_Gen->Phi();

    // if(!(isMuonInAcceptance(TIGHT, ptGen, etaGen)) || !(isMuonInAcceptance(TIGHT, pMu_Gen_Fixed->Pt(), pMu_Gen_Fixed->Eta()))) //take only muons in the fiducial area
    //   continue;
    // if(etaGen > 1.6 || pMu_Gen_Fixed->Eta() > 1.6)
    //   continue;

    if(pMu_Fixed->Pt() > 990.) //do not introduce inefficiencies, because of the fixed muon
      continue; 

    if(startFromTrack){//reject all events w/o a reco track
      if(pMu_Tk->Pt() > 990. || pMu_Tk_Fixed->Pt() > 990.){
	//printf("rejecting event:  %1.3f (fixed), %1.3f (free)\n", pMu->Pt(), pMu_Fixed->Pt());
	continue;
      }
    }
    //printf("pT of recoTrack: %1.3f (fixed), %1.3f (free)\n", pMu->Pt(), pMu_Fixed->Pt());

    Double_t pt = pMu->Pt();
    Bool_t isRECO = kFALSE;
    if(pt < 990.)
      isRECO = kTRUE;

    Bool_t isTRIG = kFALSE;
    //trigger coding: 0... not fired; -2... only neg. mu fired, 2... only pos. mu fired, 1... both muons fired; 
    //3... trigger fired in event but none of the muons could be trigger matched (maybe do not pass the selection criteria)
    TLorentzVector *dimu_Gen = 0;
    dimu_Gen = &(*pMu_Gen + *pMu_Fixed);
    Double_t dimuRap = dimu_Gen->Rapidity();
    Double_t pTGen = dimu_Gen->Pt();
    // if(pTGen < 10. || fabs(dimuRap) > 1.25) //needed for the Dimuon10 trigger flag
    //   continue;

    Double_t phiMuNeg_Gen = pMu_Gen->Phi();//negative muon
    Double_t phiMuPos_Gen = pMu_Fixed->Phi();//positive muon

    Double_t deltaPhi = phiMuNeg_Gen - phiMuPos_Gen;
    if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
    else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
    if(deltaPhi < 0.) //reject cowboys
      continue;

    // if(HLT_Dimuon10_Jpsi_Barrel_v3 == 1 || HLT_Dimuon10_Jpsi_Barrel_v3 == -2) //1.4E33
    //if(HLT_Dimuon10_Jpsi_Barrel_v6 == 1 || HLT_Dimuon10_Jpsi_Barrel_v6 == -2) //L1DoubleMu0_HighQ, no cowboys
    //if(HLT_Dimuon0_Jpsi_v3 == 1 || HLT_Dimuon0_Jpsi_v3 == -2) 
    //if(HLT_Dimuon0_Jpsi_NoVertexing_v3 == 1 || HLT_Dimuon0_Jpsi_NoVertexing_v3 == -2) 
    if(HLT_Dimuon0_Jpsi_NoVertexing_v3 == 1) 
      isTRIG = kTRUE;

    Bool_t isUseful = kFALSE;
    if(isRECO && isTRIG)
      isUseful = kTRUE;

    Int_t thisPTBin = -1, thisEtaBin = -1;
    for(int iPT = 0; iPT < kNbpT; iPT++){
      if(ptGen > pTBins[iPT] && ptGen < pTBins[iPT+1]){
	thisPTBin = iPT;
	break;
      }
    }
    for(int iEta = 0; iEta < kNbEta; iEta++){
      if(fabs(etaGen) > etaBins[iEta] && fabs(etaGen) <= etaBins[iEta+1]){
	thisEtaBin = iEta;
	break;
      }
    }
    if(thisEtaBin >= 0){
      recoEff_pT[thisEtaBin]->Fill(isRECO, ptGen);
      trigEff_pT[thisEtaBin]->Fill(isTRIG, ptGen);
      totEff_pT[thisEtaBin]->Fill(isUseful, ptGen);
      if(isUseful){
	yTimesX_pT[thisEtaBin]->Fill(ptGen, ptGen);
	y_pT[thisEtaBin]->Fill(ptGen);
      }
    }
    if(thisPTBin >= 0){
      recoEff_eta[thisPTBin]->Fill(isRECO, fabs(etaGen));
      trigEff_eta[thisPTBin]->Fill(isTRIG, fabs(etaGen));
      totEff_eta[thisPTBin]->Fill(isUseful, fabs(etaGen));
      if(isUseful){
	yTimesX_eta[thisPTBin]->Fill(fabs(etaGen), fabs(etaGen));
	y_eta[thisPTBin]->Fill(fabs(etaGen));
      }
    }

    recoEff_phi->Fill(isRECO, phiGen);
    trigEff_phi->Fill(isTRIG, phiGen);
    totEff_phi->Fill(isUseful, phiGen);

    recoEff_pT_eta->Fill(isRECO, etaGen, ptGen);
    trigEff_pT_eta->Fill(isTRIG, etaGen, ptGen);
    totEff_pT_eta->Fill(isUseful, etaGen, ptGen);
  }
}
//===============================
void SetBranches(Bool_t startFromTrack){

  data->SetBranchAddress("muNegP_Gen", &pMu_Gen);//negative muon
  data->SetBranchAddress("muPosP_Gen", &pMu_Gen_Fixed); //positive muon

  if(startFromTrack){//take the reconstructed track
    data->SetBranchAddress("muNegP_tk", &pMu_Tk);
    data->SetBranchAddress("muPosP_tk", &pMu_Tk_Fixed);
  }
  // else{
  data->SetBranchAddress("muNegP", &pMu);
  data->SetBranchAddress("muPosP", &pMu_Fixed);
  // }
  data->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v3", &HLT_Dimuon10_Jpsi_Barrel_v3);
  data->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v6", &HLT_Dimuon10_Jpsi_Barrel_v6);
  data->SetBranchAddress("HLT_Dimuon0_Jpsi_v3", &HLT_Dimuon0_Jpsi_v3);
  data->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_v3", &HLT_Dimuon0_Jpsi_NoVertexing_v3);
}

//===============================
void BookHistos(){

  Char_t name[100];
  Int_t const kNbPhi = 10;
  Double_t phi[kNbPhi+1];
  for(int iphi = 0; iphi <= kNbPhi; iphi++)
    phi[iphi] = -TMath::Pi() + (2.*TMath::Pi() / kNbPhi * iphi);

  recoEff_phi = new TEfficiency("recoEff_phi", ";#phi", kNbPhi, phi); 
  trigEff_phi = new TEfficiency("trigEff_phi", ";#phi", kNbPhi, phi); 
  totEff_phi = new TEfficiency("totEff_phi", ";#phi", kNbPhi, phi); 

  for(int iEta = 0; iEta < kNbEta; iEta++){
    sprintf(name, "recoEff_MCTRUTH_PT_AETA%d", iEta);
    recoEff_pT[iEta] = new TEfficiency(name, ";p_{T} [GeV/c]", kNbpT, pTBins); 
    sprintf(name, "trigEff_MCTRUTH_PT_AETA%d", iEta);
    trigEff_pT[iEta] = new TEfficiency(name, ";p_{T} [GeV/c]", kNbpT, pTBins); 
    sprintf(name, "totEff_MCTRUTH_PT_AETA%d", iEta);
    totEff_pT[iEta] = new TEfficiency(name, ";p_{T} [GeV/c]", kNbpT, pTBins); 
  }
  for(int iPT = 0; iPT < kNbpT; iPT++){
    sprintf(name, "recoEff_MCTRUTH_AETA_PT%d", iPT);
    recoEff_eta[iPT] = new TEfficiency(name, ";#eta", kNbEta, etaBins); 
    sprintf(name, "trigEff_MCTRUTH_AETA_PT%d", iPT);
    trigEff_eta[iPT] = new TEfficiency(name, ";#eta", kNbEta, etaBins); 
    sprintf(name, "totEff_MCTRUTH_AETA_PT%d", iPT);
    totEff_eta[iPT] = new TEfficiency(name, ";#eta", kNbEta, etaBins); 
  }
  
  recoEff_pT_eta = new TEfficiency("recoEff_MCTRUTH_pT_eta", ";#eta; p_{T} [GeV/c]", kNbEta, etaBins, kNbpT, pTBins);
  trigEff_pT_eta = new TEfficiency("trigEff_MCTRUTH_pT_eta", ";#eta; p_{T} [GeV/c]", kNbEta, etaBins, kNbpT, pTBins);
  totEff_pT_eta = new TEfficiency("totEff_MCTRUTH_pT_eta", ";#eta; p_{T} [GeV/c]", kNbEta, etaBins, kNbpT, pTBins);

  //book some histos to keep track of the "centre-of-gravity"
  for(int iEta = 0; iEta < kNbEta; iEta++){
    sprintf(name, "yTimesX_PT_AETA%d", iEta);
    yTimesX_pT[iEta] = new TH1D(name, ";p_{T} [GeV/c]", kNbpT, pTBins);
    sprintf(name, "y_PT_AETA%d", iEta);
    y_pT[iEta] = new TH1D(name, ";p_{T} [GeV/c]", kNbpT, pTBins);
  }
  for(int iPT = 0; iPT < kNbpT; iPT++){
    sprintf(name, "yTimesX_AETA_PT%d", iPT);
    yTimesX_eta[iPT] = new TH1D(name, ";#eta", kNbEta, etaBins); 
    sprintf(name, "y_AETA_PT%d", iPT);
    y_eta[iPT] = new TH1D(name, ";#eta", kNbEta, etaBins); 
  }
}

//===============================
void CalcWeightedAverage(){

  printf("calculating the <pT> from all events\n");

  Char_t name[100];
  for(int iEta = 0; iEta < kNbEta; iEta++){
    Double_t avPT[kNbpT], errL_avPT[kNbpT], errR_avPT[kNbpT];
    Double_t eff[kNbpT], errEff_pos[kNbpT], errEff_neg[kNbpT];
    printf("etaBin %d\n", iEta);
    for(int iBin = 1; iBin < yTimesX_pT[iEta]->GetNbinsX(); iBin++){
      avPT[iBin-1] = yTimesX_pT[iEta]->GetBinContent(iBin) / y_pT[iEta]->GetBinContent(iBin);
      errL_avPT[iBin-1] = avPT[iBin-1] - yTimesX_pT[iEta]->GetBinLowEdge(iBin);
      errR_avPT[iBin-1] =  yTimesX_pT[iEta]->GetBinLowEdge(iBin + 1) - avPT[iBin-1];

      eff[iBin-1] = totEff_pT[iEta]->GetEfficiency(iBin);
      errEff_neg[iBin-1] = totEff_pT[iEta]->GetEfficiencyErrorLow(iBin);
      errEff_pos[iBin-1] = totEff_pT[iEta]->GetEfficiencyErrorUp(iBin);

      printf("%1.3f < %1.3f < %1.3f --> eff = %1.3f - %1.3f + %1.3f\n", 
	     yTimesX_pT[iEta]->GetBinLowEdge(iBin), avPT[iBin-1], 
	     yTimesX_pT[iEta]->GetBinLowEdge(iBin + 1),
	     eff[iBin-1], errEff_neg[iBin-1], errEff_pos[iBin-1]);
    }
    gtotEff_pT[iEta] = new TGraphAsymmErrors(yTimesX_pT[iEta]->GetNbinsX(), avPT, eff,
					     errL_avPT, errR_avPT,
					     errEff_neg, errEff_pos);
    sprintf(name, "gtotEff_MCTRUTH_PT_AETA%d", iEta);
    gtotEff_pT[iEta]->SetName(name);
  }
}

//===============================
void WriteHistos(){

  for(int iEta = 0; iEta < kNbEta; iEta++){
    recoEff_pT[iEta]->Write();
    trigEff_pT[iEta]->Write();
    totEff_pT[iEta]->Write();
    gtotEff_pT[iEta]->Write();
  }
  for(int iPT = 0; iPT < kNbpT; iPT++){
    recoEff_eta[iPT]->Write();
    trigEff_eta[iPT]->Write();
    totEff_eta[iPT]->Write();
  }
  recoEff_phi->Write();
  trigEff_phi->Write();
  totEff_phi->Write();

  recoEff_pT_eta->Write();
  trigEff_pT_eta->Write();
  totEff_pT_eta->Write();
  
}
