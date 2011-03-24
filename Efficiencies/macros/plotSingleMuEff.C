#include "../interface/rootIncludes.inc"
#include "TGraphAsymmErrors.h"


Int_t const kNbEff = 6;
Char_t *effFileNames[kNbEff] = {"/home/hermine/CMS/Work/TnP/LucaPerrozzi/17March2011/MuonTrackingEff_17March2011.root",
				"/home/hermine/CMS/Work/TnP/Zongchang/19March2011/MuonIDEff_19March2011.root",
				"/home/hermine/CMS/Work/TnP/Xianyou/19March2011/MuonQualEff_19March2011.root",
 				"/home/hermine/CMS/Work/TnP/Francesco/23March2011/L1L2_DMu0_TriggerEfficiencies_23March2011.root",
				"/home/hermine/CMS/Work/TnP/Luigi/19March2011/L3_DoubleMu0_TriggerEfficiencies_19March2011.root",
				"/home/hermine/CMS/Work/TnP/Ilse/24March2011/HLTMuonTrack_Mu0_TkMu0_TM.root"};
//old version of scripts:
// 				"/home/hermine/CMS/Work/TnP/Francesco/18March2011/L1L2_DoubleMu0_TriggerEfficiencies_18March2011.root",
enum {TrkEff, MuIDEff, MuQualEff, L1L2Eff, L3Eff};
Char_t *effName[kNbEff] = {"TrkEff", "MuIDEff", "MuQualEff", "L1L2Eff", "L3Eff", "Mu-TkMu0"};
Char_t *effNameLong[kNbEff] = {"offline tracking efficiency", "muon identification efficiency", 
				"efficiency of muon quality cuts", "L1-L2 trigger efficiency (Mu0)", 
				"L3 trigger efficiency (Mu0)", "muon trigger efficency of TkMu0"};
Int_t const kNbEffSample = 3;
enum {DATA, MC, MCTRUTH};
Char_t *effSampleName[kNbEffSample] = {"DATA", "MC", "MCTRUTH"};
//
enum {CENTRAL, UPPER, LOWER};
Char_t *valName[3] = {"central", "lower", "upper"};
TH2D *hMuEff[kNbEff][kNbEffSample][3];
Int_t const kNbEta = 8;
Int_t const kNbPT = 12;
TGraphAsymmErrors *gEff_pT[kNbEff][kNbEffSample][kNbEta];
TGraphAsymmErrors *gEff_eta[kNbEff][kNbEffSample][kNbPT];
//
Int_t const marker[kNbEffSample] = {20, 25, 21};
Int_t const colourEta[kNbPT] = {1,2,3,4,kYellow-6,6,7,8,9,10,11};
Int_t const colourPT[kNbEta] = {1,2,3,4,kYellow-6,6,7,8};
Float_t binsEta[kNbEta] = {0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4};
Float_t binsPT[kNbPT+1] = {0.7, 1., 1.25, 1.5, 1.75, 2., 2.5, 3.5, 4., 6., 10., 20., 30.};
void LoadEfficiencies(Int_t iEff, Int_t iEffSample);
void PlotEff_DataVsMC_PT(Int_t iEff, Int_t iEta);
void PlotEff_DataVsMC_Eta(Int_t iEff, Int_t iPT);
//==============================================================
void plotSingleMuEff(){

//   for(int iEff = 0; iEff < kNbEff; iEff++){
//   for(int iEff = 5; iEff < kNbEff; iEff++){//Muon trigger eff. of TkMu0 (Ilse)
  for(int iEff = 3; iEff < 4; iEff++){//L1/L2 trigger efficiency
    for(int iEffSample = 0; iEffSample < kNbEffSample; iEffSample++) //DATA, MC, MCTRUTH
      LoadEfficiencies(iEff, iEffSample);

    for(int iEta = 0; iEta < kNbEta; iEta++){
      PlotEff_DataVsMC_PT(iEff, iEta);
    }
    for(int iPT = 0; iPT < kNbPT; iPT++){
      PlotEff_DataVsMC_Eta(iEff, iPT);
    }
  }
}

//==============================================================
void PlotEff_DataVsMC_PT(Int_t iEff, Int_t iEta){

  Char_t name[100];
  sprintf(name, "c1_%s_pT_etaBin%d", effName[iEff], iEta);
  TCanvas *c1 = new TCanvas(name, name);
  Double_t minY = 0, maxY = 1.2;
  if(iEff == 5){minY = 0.6, maxY = 1.1;}
  TH1F *hFrame1 = gPad->DrawFrame(0., minY, 31., maxY);
  hFrame1->SetXTitle("p_{T}(#mu) [GeV]");

  gEff_pT[iEff][MC][iEta]->Draw("p same");   gEff_pT[iEff][MC][iEta]->SetLineColor(2); gEff_pT[iEff][MC][iEta]->SetMarkerColor(2);
  gEff_pT[iEff][MCTRUTH][iEta]->Draw("p same"); gEff_pT[iEff][MCTRUTH][iEta]->SetLineColor(4); gEff_pT[iEff][MCTRUTH][iEta]->SetMarkerColor(4);
  gEff_pT[iEff][DATA][iEta]->Draw("p same"); gEff_pT[iEff][DATA][iEta]->SetLineColor(1); gEff_pT[iEff][DATA][iEta]->SetMarkerColor(1);

  TLatex *tex1;
  if(iEff == 5)
    tex1 = new TLatex(2., 0.95*maxY, effNameLong[iEff]);
  else
    tex1 = new TLatex(2., 0.9*maxY, effNameLong[iEff]);
  tex1->SetTextSize(0.045); tex1->Draw();
  
  TLine *line1 = new TLine(0., 1., 31., 1.);
  line1->SetLineStyle(3); line1->Draw();

  if(iEta == 0)
    sprintf(name, "|#eta(#mu)| < %1.1f", binsEta[iEta]);
  else
    sprintf(name, "%1.1f < |#eta(#mu)| < %1.1f", binsEta[iEta-1], binsEta[iEta]);
  TLegend *leg = new TLegend(0.7456897,0.1610169,0.9454023,0.3622881, name);
  leg->AddEntry(gEff_pT[iEff][DATA][iEta], "Data T&P", "pl");
  leg->AddEntry(gEff_pT[iEff][MC][iEta], "MC T&P", "pl");
  leg->AddEntry(gEff_pT[iEff][MCTRUTH][iEta], "MC Truth", "pl");
  leg->SetFillColor(0); leg->SetTextSize(0.04); 
  leg->SetBorderSize(0); leg->Draw();

  sprintf(name, "Figures/%s_pT_DataVsMC_etaBin%d.pdf", effName[iEff], iEta);  c1->Print(name);

}

//==============================================================
void PlotEff_DataVsMC_Eta(Int_t iEff, Int_t iPT){

  Char_t name[100];
  sprintf(name, "c1_%s_eta_pTBin%d", effName[iEff], iPT);
  TCanvas *c1 = new TCanvas(name, name);
  Double_t minY = 0, maxY = 1.2;
  if(iEff == 5){minY = 0.6, maxY = 1.1;}
  TH1F *hFrame1 = gPad->DrawFrame(0., minY, 2.5, maxY);
  hFrame1->SetXTitle("|#eta(#mu)|");

  gEff_eta[iEff][MC][iPT]->Draw("p same");   gEff_eta[iEff][MC][iPT]->SetLineColor(2); gEff_eta[iEff][MC][iPT]->SetMarkerColor(2);
  gEff_eta[iEff][MCTRUTH][iPT]->Draw("p same"); gEff_eta[iEff][MCTRUTH][iPT]->SetLineColor(4); gEff_eta[iEff][MCTRUTH][iPT]->SetMarkerColor(4);
  gEff_eta[iEff][DATA][iPT]->Draw("p same"); gEff_eta[iEff][DATA][iPT]->SetLineColor(1); gEff_eta[iEff][DATA][iPT]->SetMarkerColor(1);

  TLatex *tex1;
  if(iEff == 5) 
    tex1 = new TLatex(0.2, 0.95*maxY, effNameLong[iEff]);
  else
    tex1 = new TLatex(0.2, 0.9*maxY, effNameLong[iEff]);
  tex1->SetTextSize(0.045); tex1->Draw();
  
  TLine *line1 = new TLine(0., 1., 2.5, 1.);
  line1->SetLineStyle(3); line1->Draw();

  sprintf(name, "%1.1f < p_{T}(#mu) < %1.1f [GeV]", binsPT[iPT], binsPT[iPT+1]);
  TLegend *leg = new TLegend(0.30,0.1652542,0.50,0.3665254, name);
  leg->AddEntry(gEff_eta[iEff][DATA][iPT], "Data T&P", "pl");
  leg->AddEntry(gEff_eta[iEff][MC][iPT], "MC T&P", "pl");
  leg->AddEntry(gEff_eta[iEff][MCTRUTH][iPT], "MC Truth", "pl");
  leg->SetFillColor(0); leg->SetTextSize(0.04); 
  leg->SetBorderSize(0); leg->Draw();

  sprintf(name, "Figures/%s_eta_DataVsMC_pTBin%d.pdf", effName[iEff], iPT);  c1->Print(name);

}

//==============================================================
void LoadEfficiencies(Int_t iEff, Int_t iEffSample){

  TFile *fIn = new TFile(effFileNames[iEff]);
  Char_t name[100];
  for(int iPT = 0; iPT < kNbPT; iPT++){
    sprintf(name, "gEff_%s_AETA_PT%d", effSampleName[iEffSample], iPT);
    gEff_eta[iEff][iEffSample][iPT] = (TGraphAsymmErrors *) gDirectory->Get(name);
    printf("%s, %s eta differential for pT bin %d is: %p\n", effName[iEff], effSampleName[iEffSample], iPT, gEff_eta[iEff][iEffSample][iPT]);
    gEff_eta[iEff][iEffSample][iPT]->SetMarkerStyle(marker[iEffSample]);
    gEff_eta[iEff][iEffSample][iPT]->SetMarkerColor(colourEta[iPT]);
    gEff_eta[iEff][iEffSample][iPT]->SetLineColor(colourEta[iPT]);
  }
  for(int iEta = 0; iEta < kNbEta; iEta++){
    sprintf(name, "gEff_%s_PT_AETA%d", effSampleName[iEffSample], iEta);
    gEff_pT[iEff][iEffSample][iEta] = (TGraphAsymmErrors *) gDirectory->Get(name);
    printf("%s, %s pT differential for eta bin %d is: %p\n", effName[iEff], effSampleName[iEffSample], iEta, gEff_pT[iEff][iEffSample][iEta]);
    gEff_pT[iEff][iEffSample][iEta]->SetMarkerStyle(marker[iEffSample]);
    gEff_pT[iEff][iEffSample][iEta]->SetMarkerColor(colourPT[iEta]);
    gEff_pT[iEff][iEffSample][iEta]->SetLineColor(colourPT[iEta]);
  }
}
