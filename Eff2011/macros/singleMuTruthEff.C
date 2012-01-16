#include "../interface/rootIncludes.inc"
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "isMuonInAcceptance.C"
#include "calcPol.C"

TTree *data;
TLorentzVector *pMu = 0, *pMu_Fixed, *pMu_Gen = 0, *pMu_Gen_Fixed;
Int_t HLT_Dimuon10_Jpsi_Barrel_v3, HLT_Dimuon0_Jpsi_v3;

//Int_t const kNbpT = 14;
//Double_t pTBins[kNbpT+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10., 15., 20., 50.};
Int_t const kNbpT = 47;
Double_t pTBins[kNbpT+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 46., 48., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100};
Int_t const kNbEta = 10;
Double_t etaBins[kNbEta+1] = {0, 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.1, 2.4};
TEfficiency *recoEff_pT[kNbEta], *recoEff_eta[kNbpT], *recoEff_phi;
TEfficiency *trigEff_pT[kNbEta], *trigEff_eta[kNbpT], *trigEff_phi;
TEfficiency *totEff_pT[kNbEta], *totEff_eta[kNbpT], *totEff_phi;
TEfficiency *recoEff_pT_eta, *trigEff_pT_eta, *totEff_pT_eta;;

enum {LOOSE,TIGHT};//set of muon fiducial cuts

void SetBranches(Char_t *inputFile);
void FillHistos();
void BookHistos();
void WriteHistos();
//==============================
void singleMuTruthEff(Char_t *fileNameOut = "singleMuTruthEff_9Jan2012.root"){

  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/onia2MuMu_tree.root";
  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed20GeVrap1MuFree_11Jan2012.root";
  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1MuFree_11Jan2012.root";
  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1_MuFree2PT100GeV_12Jan2012.root";
  Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1_MuFree2PT100GeV_New_12Jan2012.root";
  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1_MuFree2PT100GeV_New_L1MatchingOnly_12Jan2012.root";
  //Char_t *inputFile = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/January2012/TTree_Onia2MuMu_v10_muFixed40GeVrap1_MuFree2PT100GeV_New_L2Matching_12Jan2012.root";

  SetBranches(inputFile);

  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  BookHistos();
  FillHistos();
  WriteHistos();
  fOut->Close();
}
//===============================
void FillHistos(){

  Long64_t nbEv = data->GetEntries();
  for(Long64_t iEv = 0; iEv < nbEv; iEv++){

    data->GetEntry(iEv);
    if(iEv % 100000 == 0) 
      printf("event %d out of %d\n", (Int_t) iEv, (Int_t) nbEv);

    Double_t ptGen = pMu_Gen->Pt();
    Double_t etaGen = pMu_Gen->Eta();
    Double_t phiGen = pMu_Gen->Phi();

    if(!(isMuonInAcceptance(TIGHT, ptGen, etaGen)) && !(isMuonInAcceptance(TIGHT, pMu_Gen_Fixed->Pt(), pMu_Gen_Fixed->Eta()))) //take only muons in the fiducial area
      continue;
    if(etaGen > 1.6 || pMu_Gen_Fixed->Eta() > 1.6)
      continue;

    if(pMu_Fixed->Pt() > 990.) //do not introduce inefficiencies, because of the fixed muon
      continue; 

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
    if(pTGen < 10. || fabs(dimuRap) > 1.3)
      continue;
    if(HLT_Dimuon10_Jpsi_Barrel_v3 == 1 || HLT_Dimuon10_Jpsi_Barrel_v3 == -2) 
      //if(HLT_Dimuon0_Jpsi_v3 == 1 || HLT_Dimuon0_Jpsi_v3 == -2) 
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
    }
    if(thisPTBin >= 0){
      recoEff_eta[thisPTBin]->Fill(isRECO, fabs(etaGen));
      trigEff_eta[thisPTBin]->Fill(isTRIG, fabs(etaGen));
      totEff_eta[thisPTBin]->Fill(isUseful, fabs(etaGen));
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
void SetBranches(Char_t *inputFile){

  TFile *f = new TFile(inputFile);
  data = (TTree *) gDirectory->Get("data");

  data->SetBranchAddress("muNegP_Gen", &pMu_Gen);
  data->SetBranchAddress("muPosP_Gen", &pMu_Gen_Fixed);
  data->SetBranchAddress("muNegP", &pMu);
  data->SetBranchAddress("muPosP", &pMu_Fixed);
  data->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v3", &HLT_Dimuon10_Jpsi_Barrel_v3);
  data->SetBranchAddress("HLT_Dimuon0_Jpsi_v3", &HLT_Dimuon0_Jpsi_v3);
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
  for(int iPT = 0; iPT <= kNbpT; iPT++){
    sprintf(name, "recoEff_MCTRUTH_AETA_PT%d", iPT);
    recoEff_eta[iPT] = new TEfficiency(name, ";#eta", kNbEta, etaBins); 
    sprintf(name, "trigEff_MCTRUTH_AETA_PT%d", iPT);
    trigEff_eta[iPT] = new TEfficiency(name, ";#eta", kNbEta, etaBins); 
    sprintf(name, "totEff_MCTRUTH_AETA_PT%d", iPT);
    totEff_eta[iPT] = new TEfficiency(name, ";#eta", kNbEta, etaBins); 
  }

  recoEff_pT_eta = new TEfficiency(name, ";#eta; p_{T} [GeV/c]", kNbEta, etaBins, kNbpT, pTBins);
  trigEff_pT_eta = new TEfficiency(name, ";#eta; p_{T} [GeV/c]", kNbEta, etaBins, kNbpT, pTBins);
  totEff_pT_eta = new TEfficiency(name, ";#eta; p_{T} [GeV/c]", kNbEta, etaBins, kNbpT, pTBins);
}
//===============================
void WriteHistos(){

  for(int iEta = 0; iEta < kNbEta; iEta++){
    recoEff_pT[iEta]->Write();
    trigEff_pT[iEta]->Write();
    totEff_pT[iEta]->Write();
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
