#include "../interface/rootIncludes.inc"
#include "TEfficiency.h"

//tracker80 bins:
Int_t const kNbpT = 14;
Double_t pTBins[kNbpT+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10., 15., 20., 50.};
Int_t const kNbEta = 10;
Double_t etaBins[kNbEta+1] = {0, 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.1, 2.4};

Int_t const kNbType = 3;
enum{TNP, TNPALLINONE, MCTRUTH};
TGraphAsymmErrors *gEff[2][kNbEta];
TEfficiency *tEff[kNbEta];
TH1F *hRatio[2][kNbEta];

void LoadEfficiencies();
void DrawEff(Int_t iEta, Bool_t drawTNP, Bool_t drawAllInOne);
void MakeRatio(Int_t iEta, Int_t index);
void DrawRatio(Int_t iEta);
//==============================
void compareEff(){

  LoadEfficiencies();
  // for(int iEta = 0; iEta < kNbEta-1; iEta++)
  //   DrawEff(iEta, kTRUE, kTRUE);

  // for(int iEta = 0; iEta < kNbEta-1; iEta++)
  //   DrawEff(iEta, kTRUE, kFALSE);

  // for(int iEta = 0; iEta < kNbEta-1; iEta++)
  //   DrawEff(iEta, kFALSE, kTRUE);

  for(int iEta = 0; iEta < kNbEta-1; iEta++){
    MakeRatio(iEta, TNP);
    MakeRatio(iEta, TNPALLINONE);
    DrawRatio(iEta);
  }
}
//==============================
void DrawRatio(Int_t iEta){

  Char_t name[100];
  sprintf(name, "c2_%d", iEta);
  TCanvas *c2 = new TCanvas(name, name);
  TH1F *hFrame2 = gPad->DrawFrame(0., 0.5, 20.1, 1.5);
  hFrame2->SetXTitle("p_{T} [GeV/c]");
  hFrame2->SetYTitle("T&P / Truth");

  hRatio[TNP][iEta]->SetMarkerStyle(20);
  hRatio[TNP][iEta]->SetMarkerColor(1);
  hRatio[TNP][iEta]->Draw("p same");

  hRatio[TNPALLINONE][iEta]->SetMarkerStyle(24);
  hRatio[TNPALLINONE][iEta]->SetMarkerColor(4);
  hRatio[TNPALLINONE][iEta]->Draw("p same");

  TLine *line2 = new TLine(0., 1., 20.1, 1.);
  line2->SetLineStyle(3); line2->Draw();

  sprintf(name, "%1.1f < |#eta| < %1.1f", etaBins[iEta], etaBins[iEta+1]);
  TLegend *leg2 = new TLegend(0.5991379,0.7436441,0.9022989,0.8813559, name);
  leg2->AddEntry(hRatio[TNPALLINONE][iEta], "all in one MC T&P", "pl");
  leg2->AddEntry(hRatio[TNP][iEta], "product MC T&P", "pl");
  leg2->SetTextSize(0.04); leg2->SetBorderSize(0);
  leg2->SetFillColor(0); leg2->Draw();

  sprintf(name, "Figures/ratioEff_Eta%d.pdf", iEta);
  c2->Print(name);
}

//==============================
void MakeRatio(Int_t iEta, Int_t index){

  Char_t name[100];
  sprintf(name, "hRatio_%d_eta%d", index, iEta);
  hRatio[index][iEta] = new TH1F(name, ";p_{T} [GeV/c]", kNbpT, pTBins);
  Double_t pT, eff, effTruth;
  Int_t thisBin = -1;
  for(int iPT = 0; iPT < kNbpT; iPT++){
    gEff[index][iEta]->GetPoint(iPT, pT, eff);
    // for(int iPTBin = 0; iPTBin < kNbpT; iPTBin++){
    //   if(pT > pTBins[iPTBin] && pT < pTBins[iPTBin+1]){
    // 	thisBin = iPTBin;
    // 	break;
    //   }
    // }

    thisBin = tEff[iEta]->FindFixBin(pT);
    effTruth = tEff[iEta]->GetEfficiency(thisBin);
    printf("pTBin %d, eff %1.3f, effTruth %1.3f\n", thisBin, eff, effTruth);
    hRatio[index][iEta]->Fill(pT, eff / effTruth);
  }
}

//==============================
void DrawEff(Int_t iEta, Bool_t drawTNP, Bool_t drawAllInOne){

  Char_t name[100];
  sprintf(name, "c_iEta%d", iEta); 
  TCanvas *c1 = new TCanvas(name);
  TH1F *hFrame1 = gPad->DrawFrame(0., 0., 20.1, 1.1);
  hFrame1->SetXTitle("p_{T} [GeV/c]");
  hFrame1->SetYTitle("single muon efficiency");

  if(drawTNP)
    gEff[TNP][iEta]->Draw("p same");
  if(drawAllInOne)
    gEff[TNPALLINONE][iEta]->Draw("p same");

  tEff[iEta]->Draw("p same");

  TLine *line = new TLine(0., 1., 20.1, 1.);
  line->SetLineStyle(3);
  line->Draw();

  sprintf(name, "%1.1f < |#eta| < %1.1f", etaBins[iEta], etaBins[iEta+1]);
  TLegend *tex1 = new TLegend(0.5847701,0.1758475,0.9353448,0.3834746, name);
  if(drawTNP)
    tex1->AddEntry(gEff[TNP][iEta], "product MC T&P", "pl");
  if(drawAllInOne)
    tex1->AddEntry(gEff[TNPALLINONE][iEta], "all in one MC T&P", "pl");
  tex1->AddEntry(tEff[iEta], "MC truth", "pl");
  tex1->SetFillColor(0); tex1->SetTextSize(0.04);
  tex1->SetBorderSize(0.);
  tex1->Draw(); 

  if(drawTNP && drawAllInOne)
    sprintf(name, "Figures/tnpEff_allInOneEff_truthEff_Eta%d.pdf", iEta);
  else if(drawTNP)
    sprintf(name, "Figures/tnpEff_truthEff_Eta%d.pdf", iEta);
  else if(drawAllInOne)
    sprintf(name, "Figures/allInOneEff_truthEff_Eta%d.pdf", iEta);
  c1->Print(name);
}

//==============================
void LoadEfficiencies(){

  TFile *fTnP = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Ilse/8Dec2011/EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Trk80Cuts_19Nov2011.root");
  Char_t name[100];
  for(int iEta = 0; iEta < kNbEta; iEta++){
    sprintf(name, "gEff_MC_PT_AETA%d", iEta);
    gEff[TNP][iEta] = (TGraphAsymmErrors *) fTnP->Get(name);
    gEff[TNP][iEta]->SetMarkerStyle(20);
    gEff[TNP][iEta]->SetMarkerColor(1);
    gEff[TNP][iEta]->SetLineColor(1);
  }

  TFile *fTnPAllInOne = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Linlin/7Dec2011/singleMuonEfficiency_ProbeTrackMatched_data_mc_pt_abseta_tracker80Cuts_7Dec2011.root");
  for(int iEta = 0; iEta < kNbEta-1; iEta++){
    sprintf(name, "gEff_MC_PT_ABSETA%d", iEta);
    gEff[TNPALLINONE][iEta] = (TGraphAsymmErrors *) fTnPAllInOne->Get(name);
    gEff[TNPALLINONE][iEta]->SetMarkerStyle(24);
    gEff[TNPALLINONE][iEta]->SetMarkerColor(4);
    gEff[TNPALLINONE][iEta]->SetLineColor(4);
  }

  TFile *fTruth = new TFile("singleMuTruthEff_12Jan2012_40GeVrap1_2pT100GeV_New_EtaCut.root");
  for(int iEta = 0; iEta < kNbEta; iEta++){
    sprintf(name, "totEff_MCTRUTH_PT_AETA%d", iEta);
    tEff[iEta] = (TEfficiency *) fTruth->Get(name);
    tEff[iEta]->SetMarkerStyle(21);
    tEff[iEta]->SetMarkerSize(0.4);
    tEff[iEta]->SetMarkerColor(2);
    tEff[iEta]->SetLineColor(2);
  }
}
