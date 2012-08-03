#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"
#include "TEfficiency.h"

TEfficiency *totEff2D_pol_pT_rap[eff::kNbFrames][eff::kNbPTMaxBins][eff::kNbRapForPTBins+1];
TH2D *hRelError[eff::kNbFrames][eff::kNbPTMaxBins][eff::kNbRapForPTBins+1];
void LoadEff(Char_t *fName);
void PlotEff(Int_t iFrame, Int_t iRapBin, Int_t iPTBin);
void PlotError(Int_t iFrame, Int_t iRapBin, Int_t iPTBin);
//============================
void plotDimuAccAndEff(Char_t *fileNameIn = "MCTruthEffAndAcc_UPS1S_HLTDimuon7UpsilonBarrelV1_1sigma_2Aug2012.root"){
  //Char_t *fileNameIn = "MCTruthEffAndAcc_HLTDimuon10JpsiBarrel_2Aug2012.root"){

  LoadEff(fileNameIn);

  gStyle->SetPaintTextFormat("5.2g");
  //for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
  for(int iFrame = 0; iFrame < 3; iFrame++)
    for(int iRapBin = 1; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
      for(int iPTBin = 1; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
	PlotEff(iFrame, iRapBin, iPTBin);	

  for(int iFrame = 0; iFrame < 3; iFrame++)
    for(int iRapBin = 1; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
      for(int iPTBin = 1; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
	PlotError(iFrame, iRapBin, iPTBin);	

}

//============================
void PlotError(Int_t iFrame, Int_t iRapBin, Int_t iPTBin){

  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.14);
  Char_t name[100];
  sprintf(name, "c2_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
  TCanvas *c2 = new TCanvas(name, name, 700, 700);
  hRelError[iFrame][iPTBin][iRapBin]->Draw("colz text");

  sprintf(name, "Figures/relUncertainty_dimuTruthAccAndEff_%s_rap%d_pT%d.pdf", eff::frameLabel[iFrame], iRapBin, iPTBin);
  c2->Print(name);

}

//============================
void PlotEff(Int_t iFrame, Int_t iRapBin, Int_t iPTBin){

  gStyle->SetPadRightMargin(0.14);
  Char_t name[100];
  sprintf(name, "c1_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
  TCanvas *c1 = new TCanvas(name, name, 700, 700);
  //TH1F *hFrame1 = gPad->DrawFrame(0., 0., 50., 1.1);
  totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Draw("colz text");

  sprintf(name, "Figures/dimuTruthAccAndEff_%s_rap%d_pT%d.pdf", eff::frameLabel[iFrame], iRapBin, iPTBin);
  c1->Print(name);

}

//============================
void LoadEff(Char_t *fName){

  Char_t name[100], title[100];
  TFile *fIn = new TFile(fName);
  Double_t nEntries;
  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    //for(int iFrame = 0; iFrame < 3; iFrame++){
    for(int iRapBin = 1; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 1; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

	sprintf(name, "totEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TEfficiency *) gDirectory->Get(name);
	printf("frame %d, rap %d, pT %d, %p\n", iFrame, iRapBin, iPTBin, totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]);
	// totEff2D_pol_pT_rap[iFrame][iRapBin][iPTBin]->SetMarkerColor(index+1);
	// totEff2D_pol_pT_rap[iFrame][iRapBin][iPTBin]->SetLineColor(index+1);
	sprintf(name, "hRelErr_%s", totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetName());
	sprintf(title, ";%s;%s", totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetXaxis()->GetTitle(), totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetYaxis()->GetTitle());
	hRelError[iFrame][iPTBin][iRapBin] = new TH2D(name, title, totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetNbinsX(), -1., 1., totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetNbinsY(), -180., 180);
	for(int iBinX = 1; iBinX <= hRelError[iFrame][iPTBin][iRapBin]->GetNbinsX(); iBinX++){
	  for(int iBinY = 1; iBinY <= hRelError[iFrame][iPTBin][iRapBin]->GetNbinsY(); iBinY++){
	    nEntries = totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetBinContent(iBinX, iBinY);
	    if(nEntries > 0)
	      hRelError[iFrame][iPTBin][iRapBin]->SetBinContent(iBinX, iBinY, 100./sqrt(nEntries));
	  }
	}
      }
    }
  }
  printf("after loading eff\n");
}
