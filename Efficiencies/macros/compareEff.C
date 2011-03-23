#include "/home/hermine/Work/rootIncludes.inc"
#include "TEfficiency.h"

Char_t const *trigName = "HLT_DoubleMu0";

Int_t const kNbEff = 3;
enum {RECO, TRIG, TOT};
Char_t *effName[kNbEff] = {"reco", "trig", "tot"};

Int_t const kNbVar = 3;
Char_t *varName[kNbVar] = {"pT", "y", "phi"};
enum {PT,Y,PHI};
TEfficiency *gEffMCTruth1D[kNbEff][kNbVar];
TH1D *hEffMCTruth1D[kNbEff][kNbVar];
enum {ABSETA, ETA};
TEfficiency *gEffMCTruth2D[kNbEff][2];
TH2D *hEffMCTruth2D[kNbEff][2];

TH1D *hEffTP[kNbEff][kNbVar];
TH2D *hEffTP2D[kNbEff][2];

void LoadMCTruthEff(Char_t *fileNameIn, Int_t iEff);
void LoadTPEfficiencies(Char_t *fileNameIn, Int_t iEff);
void Plot1DEff(Int_t iEff, Int_t iVar);
void PlotRatio1D(Int_t iEff, Int_t iVar);
void Plot2DEff(Int_t iEff, Int_t iRap);
void PlotRatio2D(Int_t iEff, Int_t iRap);
//=======================
void compareEff(Char_t *fileNameMCTruth = "MCTruthEff_HLTDoubleMu0_22March2011.root",
// 		Char_t *fileNameTPEff = "MCTnPEff_SingleMuMCTruth_HLTDoubleMu0_21March2011.root"){
		Char_t *fileNameTPEff = "MCTnPEff_HLTDoubleMu0_22March2011.root"){

  for(int iEff = 0; iEff < kNbEff; iEff++)
    LoadMCTruthEff(fileNameMCTruth, iEff);
    
  for(int iEff = 0; iEff < kNbEff; iEff++)
    LoadTPEfficiencies(fileNameTPEff, iEff);

  for(int iEff = 0; iEff < kNbEff; iEff++)
    for(int iVar = 0; iVar < kNbVar; iVar++)
      Plot1DEff(iEff, iVar);

  for(int iEff = 0; iEff < kNbEff; iEff++)
    for(int iVar = 0; iVar < kNbVar; iVar++)
      PlotRatio1D(iEff, iVar);

  for(int iEff = 0; iEff < kNbEff; iEff++)
    for(int iRap = 0; iRap < 2; iRap++)
      Plot2DEff(iEff, iRap);

  for(int iEff = 0; iEff < kNbEff; iEff++)
    for(int iRap = 0; iRap < 2; iRap++)
      PlotRatio2D(iEff, iRap);

}

//========================
void PlotRatio1D(Int_t iEff, Int_t iVar){

  Char_t name[100];
  sprintf(name, "c1_%sEff_%s", effName[iEff], varName[iVar]);
  TCanvas *c1 = new TCanvas(name, name);

  sprintf(name, "hRatio_%sEff_%s", effName[iEff], varName[iVar]);
  TH1D *hRatio = (TH1D *) hEffTP[iEff][iVar]->Clone(name);
  hRatio->Divide(hEffMCTruth1D[iEff][iVar]);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMinimum(0.);
  hRatio->SetMaximum(1.5);
  hRatio->Draw("p");

  sprintf(name, "Figures/ratio_%sEff_%s.pdf", effName[iEff], varName[iVar]);
  c1->Print(name);
}

//=====================================
void Plot1DEff(Int_t iEff, Int_t iVar){

  gStyle->SetOptStat(0);

  Char_t name[100];
  sprintf(name, "c2_%sEff_%s", effName[iEff], varName[iVar]);
  TCanvas *c1 = new TCanvas(name, name);

  hEffTP[iEff][iVar]->SetMinimum(0.);
  hEffTP[iEff][iVar]->SetMaximum(1.);
  hEffTP[iEff][iVar]->Draw("p");
  hEffTP[iEff][iVar]->SetMarkerStyle(20);
  gEffMCTruth1D[iEff][iVar]->Draw("p same");
  gEffMCTruth1D[iEff][iVar]->SetMarkerStyle(24);

  if(iVar == 0){
    if(iEff == 1)
      sprintf(name, "eff. for %s", trigName);
    else
      sprintf(name, "%s efficiency", effName[iEff]);
    TLegend *leg1 = new TLegend(0.65,0.1610169,0.8,0.3601695, name);
    leg1->AddEntry(gEffMCTruth1D[iEff][iVar], "MC truth", "p");
    leg1->AddEntry(hEffTP[iEff][iVar], "MC T&P", "p");
    leg1->SetFillColor(0);  leg1->SetBorderSize(0);
    leg1->SetTextSize(0.04);
    leg1->Draw("same");
  }

  sprintf(name, "Figures/eff1D_%sEff_%s.pdf", effName[iEff], varName[iVar]);
  c1->Print(name);
}

//=====================================
void Plot2DEff(Int_t iEff, Int_t iRap){

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("5.2f");
  gStyle->SetPadRightMargin(0.1);

  Char_t name[100];
  //plot the TnP MC histogram
  sprintf(name, "MCTnP_%sEff_%d", effName[iEff], iRap);
  TCanvas *c1 = new TCanvas(name, name);

  hEffTP2D[iEff][iRap]->SetMinimum(0.);
  hEffTP2D[iEff][iRap]->SetMaximum(1.);
  hEffTP2D[iEff][iRap]->Draw("colz text");

  sprintf(name, "Figures/eff2D_TnPMC_%sEff_rap%d.pdf", effName[iEff], iRap);
  c1->Print(name);


  //plot the MC truth histogram
  sprintf(name, "MCTruth_%sEff_%d", effName[iEff], iRap);
  TCanvas *c2 = new TCanvas(name, name);

  gEffMCTruth2D[iEff][iRap]->Draw("colz text");

  sprintf(name, "Figures/eff2D_MCTruth_%sEff_rap%d.pdf", effName[iEff], iRap);
  c2->Print(name);
}

//========================
void PlotRatio2D(Int_t iEff, Int_t iRap){

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("5.2f");

  Char_t name[100];
  sprintf(name, "ratio2D_%sEff_%d", effName[iEff], iRap);
  TCanvas *c1 = new TCanvas(name, name);

  sprintf(name, "hRatio2D_%sEff_rap%d", effName[iEff], iRap);
  TH2D *hRatio = (TH2D *) hEffTP2D[iEff][iRap]->Clone(name);
  hRatio->Divide(hEffMCTruth2D[iEff][iRap]);
  hRatio->SetMinimum(0.);
  hRatio->SetMaximum(1.5);
  hRatio->Draw("colz text");

  sprintf(name, "Figures/ratio2D_%sEff_rap%d.pdf", effName[iEff], iRap);
  c1->Print(name);
}

//========================
void LoadTPEfficiencies(Char_t *fileNameIn, Int_t iEff){

  TFile *fIn = new TFile(fileNameIn);
  Char_t name[100];
  for(int iVar = 0; iVar < kNbVar; iVar++){
    sprintf(name, "%sEff_%s", effName[iEff], varName[iVar]);
    hEffTP[iEff][iVar] = (TH1D *) gDirectory->Get(name);
    sprintf(name, "%sEffMCTruth_%s", effName[iEff], varName[iVar]);
    hEffTP[iEff][iVar]->SetName(name);
  }

  sprintf(name, "%sEff2D_pT_rap", effName[iEff]);
  hEffTP2D[iEff][ABSETA] = (TH2D *) gDirectory->Get(name);
  sprintf(name, "%sEffMCTruth2D_pT_rap", effName[iEff]);
  hEffTP2D[iEff][ABSETA]->SetName(name);

  sprintf(name, "%sEff2D_pT_rapNP", effName[iEff]);
  hEffTP2D[iEff][ETA] = (TH2D *) gDirectory->Get(name);
  sprintf(name, "%sEffMCTruth2D_pT_rapNP", effName[iEff]);
  hEffTP2D[iEff][ETA]->SetName(name);

}

//========================
void LoadMCTruthEff(Char_t *fileNameIn, Int_t iEff){

  TFile *fIn = new TFile(fileNameIn);
  Char_t name[100];
  TH1D *hPassed, *hTot;
  for(int iVar = 0; iVar < kNbVar; iVar++){
    sprintf(name, "%sEff_%s", effName[iEff], varName[iVar]);
    gEffMCTruth1D[iEff][iVar] = (TEfficiency *) gDirectory->Get(name);
    sprintf(name, "%sEffMCTruth_%s", effName[iEff], varName[iVar]);
    gEffMCTruth1D[iEff][iVar]->SetName(name);

    //copy the values into a histogram:
    hPassed = (TH1D *) gEffMCTruth1D[iEff][iVar]->GetPassedHistogram();
    hTot = (TH1D *) gEffMCTruth1D[iEff][iVar]->GetTotalHistogram();
    sprintf(name, "h%sEffMCTruth_%s", effName[iEff], varName[iVar]);
    hEffMCTruth1D[iEff][iVar] = (TH1D *) hPassed->Clone(name);
    hEffMCTruth1D[iEff][iVar]->Divide(hTot);
  }

  sprintf(name, "%sEff2D_pT_rap", effName[iEff]);
  gEffMCTruth2D[iEff][ABSETA] = (TEfficiency *) gDirectory->Get(name);
  sprintf(name, "%sEffMCTruth2D_pT_rap", effName[iEff]);
  gEffMCTruth2D[iEff][ABSETA]->SetName(name);

  sprintf(name, "%sEff2D_pT_rapNP", effName[iEff]);
  gEffMCTruth2D[iEff][ETA] = (TEfficiency *) gDirectory->Get(name);
  sprintf(name, "%sEffMCTruth2D_pT_rapNP", effName[iEff]);
  gEffMCTruth2D[iEff][ETA]->SetName(name);
  
  //copy the values into a histogram:
  TH2D *hPassed2D, *hTot2D;
  hPassed2D = (TH2D *) gEffMCTruth2D[iEff][ABSETA]->GetPassedHistogram();
  hTot2D = (TH2D *) gEffMCTruth2D[iEff][ABSETA]->GetTotalHistogram();
  sprintf(name, "h%sEffMCTruth2D_pT_rap", effName[iEff]);
  hEffMCTruth2D[iEff][ABSETA] = (TH2D *) hPassed2D->Clone(name);
  hEffMCTruth2D[iEff][ABSETA]->Divide(hTot2D);

  hPassed2D = (TH2D *) gEffMCTruth2D[iEff][ETA]->GetPassedHistogram();
  hTot2D = (TH2D *) gEffMCTruth2D[iEff][ETA]->GetTotalHistogram();
  sprintf(name, "h%sEffMCTruth2D_pT_rapNP", effName[iEff]);
  hEffMCTruth2D[iEff][ETA] = (TH2D *) hPassed2D->Clone(name);
  hEffMCTruth2D[iEff][ETA]->Divide(hTot2D);

}
