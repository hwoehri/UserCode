#include "rootIncludes.inc"
#include "commonVar.h"
#include "RooPlot.h"
#include "TROOT.h"


//mass and lifetime
TGraphAsymmErrors *gMass[kNbFrames][kNbRap][kNbPTMaxBins];
TGraphAsymmErrors *gMass_Fit[kNbComp][kNbFrames][kNbRap][kNbPTMaxBins];
Double_t maxMass[kNbFrames][kNbRap][kNbPTMaxBins];
Double_t chi2Mass[kNbFrames][kNbRap][kNbPTMaxBins];

TGraphAsymmErrors *gLifetime[kNbSigBG][kNbFrames][kNbRap][kNbPTMaxBins];
TGraphAsymmErrors *gLifetime_Fit[kNbSigBG][kNbComp][kNbFrames][kNbRap][kNbPTMaxBins];
Double_t maxLifetime[kNbSigBG][kNbFrames][kNbRap][kNbPTMaxBins];
Double_t chi2Lifetime[kNbSigBG][kNbFrames][kNbRap][kNbPTMaxBins];

//polarization variables
TGraphAsymmErrors *gPhi[kNbSpecies][kNbFrames][kNbRap][kNbPTMaxBins];
TGraphAsymmErrors *gCosTh[kNbSpecies][kNbFrames][kNbRap][kNbPTMaxBins];
TGraphAsymmErrors *gPhi_Fit[kNbSpecies][kNbComp][kNbFrames][kNbRap][kNbPTMaxBins];
TGraphAsymmErrors *gCosTh_Fit[kNbSpecies][kNbComp][kNbFrames][kNbRap][kNbPTMaxBins];

//polarization variables
TGraphAsymmErrors *gPhiBG[kNbBG][kNbFrames][kNbRap][kNbPTMaxBins];
TGraphAsymmErrors *gCosThBG[kNbBG][kNbFrames][kNbRap][kNbPTMaxBins];
TGraphAsymmErrors *gPhi_FitBG[kNbBG][kNbComp-1][kNbFrames][kNbRap][kNbPTMaxBins];
TGraphAsymmErrors *gCosTh_FitBG[kNbBG][kNbComp-1][kNbFrames][kNbRap][kNbPTMaxBins];
TH1F *hPhiBG[kNbBG][kNbFrames][kNbRap][kNbPTMaxBins];
TH1F *hCosThBG[kNbBG][kNbFrames][kNbRap][kNbPTMaxBins];

//contour plots
TGraphAsymmErrors *gContour[kNbVarComb][kNbSpecies][kNbFrames][kNbRap][kNbPTMaxBins];
TLine *gContour_line1[kNbVarComb][kNbSpecies][kNbFrames][kNbRap][kNbPTMaxBins];
TLine *gContour_line2[kNbVarComb][kNbSpecies][kNbFrames][kNbRap][kNbPTMaxBins];
TLine *gContour_marker[kNbVarComb][kNbSpecies][kNbFrames][kNbRap][kNbPTMaxBins];

Int_t const kNbVar = 2;
Double_t maximumPol[kNbVar][kNbSpecies][kNbFrames][kNbRap][kNbPTMaxBins];
Double_t chi2[kNbVar][kNbSpecies][kNbFrames][kNbRap][kNbPTMaxBins];
Double_t maximumPolBG[kNbVar][kNbBG][kNbFrames][kNbRap][kNbPTMaxBins];

void GetPolGraphs(Int_t iFrame, Int_t iRap, Int_t iPT);
void PlotPolGraphs(Int_t iFrame, Int_t iRap, Int_t iPT, Int_t iSpecies);
void GetMassLifetime(Int_t iFrame, Int_t iRap, Int_t iPT);
void PlotMassLifetime(Int_t iFrame, Int_t iRap, Int_t iPT);
void GetContour(Int_t iFrame, Int_t iRap, Int_t iPT);
void PlotContour(Int_t iFrame, Int_t iRap, Int_t iPT, Int_t iSpecies);
void GetPolGraphsLRSidebands(Int_t iFrame, Int_t iRap, Int_t iPT);
void PlotPolGraphsLRSidebands(Int_t iFrame, Int_t iRap, Int_t iPT, Int_t iSideband);
//=======================
void plotPolFits(){

  //produces the full suite of plots in the signal region
  for(int iRap = 1; iRap <= 2; iRap++){
    for(int iFrame = CS; iFrame <= HX; iFrame++){
      for(int iPT = 5; iPT <= 8; iPT++){

	if(iRap == 1 && (iPT == 5 || iPT == 8)) continue;
	if(iRap == 2 && iPT == 8 && iFrame == HX) continue;

	//a) mass and lifetime distributions
	GetMassLifetime(iFrame, iRap, iPT);
	PlotMassLifetime(iFrame, iRap, iPT);
	//b) cosTheta and phi distributions
	GetPolGraphs(iFrame, iRap, iPT);
	PlotPolGraphs(iFrame, iRap, iPT, P);  PlotPolGraphs(iFrame, iRap, iPT, NP);
	//c) 68% C.L. contour plots for phi vs. cosTheta
	GetContour(iFrame, iRap, iPT);
	PlotContour(iFrame, iRap, iPT, P);  PlotContour(iFrame, iRap, iPT, NP);
      }
    }
  }

  //produces the "pedagogical" cosTheta and phi figures of the side band regions
  for(int iRap = 1; iRap <= 2; iRap++){
    for(int iFrame = CS; iFrame <= HX; iFrame++){
      for(int iPT = 5; iPT <= 8; iPT++){

	if(iRap == 1 && (iPT == 5 || iPT == 8)) continue;
	if(iRap == 2 && iPT == 8 && iFrame == HX) continue;

	GetPolGraphsLRSidebands(iFrame, iRap, iPT);
	PlotPolGraphsLRSidebands(iFrame, iRap, iPT, L); PlotPolGraphsLRSidebands(iFrame, iRap, iPT, R); 
      }
    }
  }
}

//=======================
void PlotContour(Int_t iFrame, Int_t iRap, Int_t iPT, Int_t iSpecies){

  Char_t name[100];
  gStyle->SetPadLeftMargin(0.16);

  //===================================
  //Lambda_theta vs Lambda_phi
  //===================================
  Double_t minXLine = gContour_line1[TH_PHI][iSpecies][iFrame][iRap][iPT]->GetX1();
  Double_t maxXLine = gContour_line1[TH_PHI][iSpecies][iFrame][iRap][iPT]->GetX2();
  Double_t minYLine = gContour_line2[TH_PHI][iSpecies][iFrame][iRap][iPT]->GetY1();
  Double_t maxYLine = gContour_line2[TH_PHI][iSpecies][iFrame][iRap][iPT]->GetY2();

  Double_t minX = -0.3, maxX = 0.3;
  Double_t minY = -0.3, maxY = 0.3;
  if(minXLine < minX || maxXLine > maxX || minYLine < minY || maxYLine > maxY){ 
    minX = minXLine - 0.1; maxX = minXLine + 0.5; 
    minY = minYLine - 0.1; maxY = minYLine + 0.5;
  }

  sprintf(name, "LthVsLPhi_%s_%s_rap%d_pT%d", speciesLabel[iSpecies], frameLabel[iFrame], iRap, iPT);
  TCanvas *c1 = new TCanvas(name, name, 500, 500);
  TH1F *hFrame1 = gPad->DrawFrame(minX, minY, maxX, maxY);
  sprintf(name, "#lambda_{#phi}^{%s} (%s)", speciesLabel[iSpecies], frameLabel[iFrame]);   
  hFrame1->SetYTitle(name);
  sprintf(name, "#lambda_{#theta}^{%s} (%s)", speciesLabel[iSpecies], frameLabel[iFrame]);   
  hFrame1->SetXTitle(name);
  hFrame1->SetTitleOffset(1.9, "y");

  gContour[TH_PHI][iSpecies][iFrame][iRap][iPT]->Draw("l same");
  gContour_line1[TH_PHI][iSpecies][iFrame][iRap][iPT]->Draw("l same");
  gContour_line2[TH_PHI][iSpecies][iFrame][iRap][iPT]->Draw("l same");
  gContour_marker[TH_PHI][iSpecies][iFrame][iRap][iPT]->Draw("p same");

  if(iRap == 1)
    sprintf(name, "|y| < %1.1f", rapForPTRange[iRap]);
  else 
    sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[iRap-1], rapForPTRange[iRap]);
  TLatex *tex1 = new TLatex(minX+0.02, minY+0.05, name);
  tex1->SetTextSize(0.04); tex1->Draw();
  sprintf(name, "%1.1f < p_{T} < %1.1f GeV", pTRange[iRap][iPT-1], pTRange[iRap][iPT]);
  tex1->DrawLatex(minX+0.02, minY+0.02, name);

  TLine *line1 = new TLine(minX, 0., maxX, 0.);
  line1->Draw(); line1->SetLineStyle(3);
  line1->DrawLine(0., minY, 0., maxY);

  TF1 *func1 = new TF1("lambdaTilde0", "-x[0]/3.", minX, maxX);
  func1->SetLineColor(4); func1->SetLineWidth(1);
  func1->Draw("same");
//   if(maxX < 0.25){
//     TLatex *tex1a = new TLatex(0.13, -0.1, "#tilde{#lambda} = 0");
//     tex1a->SetTextColor(4); tex1a->Draw();
//   }

  sprintf(name, "Figures/contour_LPhiVsLTheta_%s_%s_rap%d_pT%d.pdf", speciesLabel[iSpecies], frameLabel[iFrame], iRap, iPT);
  c1->Print(name);
}

//=======================
void PlotMassLifetime(Int_t iFrame, Int_t iRap, Int_t iPT){

  Char_t name[100];

  //===================================
  //Mass
  //===================================
  sprintf(name, "mass_%s_rap%d_pT%d", frameLabel[iFrame], iRap, iPT);
  TCanvas *c1 = new TCanvas(name, name);
  TH1F *hFrame1 = gPad->DrawFrame(2.7, 0., 3.5, maxMass[iFrame][iRap][iPT]);
  sprintf(name, "Mass [GeV]");
  hFrame1->SetXTitle(name);

  gMass[iFrame][iRap][iPT]->Draw("p same");
  gMass_Fit[P][iFrame][iRap][iPT]->Draw("l same"); //P
  gMass_Fit[NP][iFrame][iRap][iPT]->Draw("l same"); //NP
  gMass_Fit[BKG][iFrame][iRap][iPT]->Draw("l same"); //BG
  gMass_Fit[TOT][iFrame][iRap][iPT]->Draw("l same"); //Tot

  if(iRap == 1)
    sprintf(name, "|y| < %1.1f", rapForPTRange[iRap]);
  else 
    sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[iRap-1], rapForPTRange[iRap]);
  TLatex *tex1 = new TLatex(3.32, 0.6*maxMass[iFrame][iRap][iPT], name);
  tex1->SetTextSize(0.04); tex1->Draw();
  sprintf(name, "%1.1f < p_{T} < %1.1f GeV", pTRange[iRap][iPT-1], pTRange[iRap][iPT]);
  tex1->DrawLatex(3.24, 0.53*maxMass[iFrame][iRap][iPT], name);

//   sprintf(name, "#chi^{2}/ndf = %1.2f", chi2Mass[iFrame][iRap][iPT]);
//   tex1->DrawLatex(2.73, 0.9*maxMass[iFrame][iRap][iPT], name);

  TLegend *leg1 = new TLegend(0.6752874,0.7394068,0.96,0.9385593);
  leg1->AddEntry(gMass_Fit[TOT][iFrame][iRap][iPT], "sum", "l");
  leg1->AddEntry(gMass_Fit[P][iFrame][iRap][iPT], "prompt J/#psi", "l");
  leg1->AddEntry(gMass_Fit[NP][iFrame][iRap][iPT], "non-prompt J/#psi", "l");
  leg1->AddEntry(gMass_Fit[BKG][iFrame][iRap][iPT], "background", "l");
  leg1->SetTextSize(0.04); leg1->SetFillColor(0);
  leg1->SetBorderSize(0);  leg1->Draw();

  sprintf(name, "Figures/mass_%s_rap%d_pT%d.pdf", frameLabel[iFrame], iRap, iPT);
  c1->Print(name);

  //===================================
  //Lifetime in the BG region
  //===================================
  sprintf(name, "lifetime_BG_%s_rap%d_pT%d", frameLabel[iFrame], iRap, iPT);
  TCanvas *c2 = new TCanvas(name, name);
  gPad->SetLogy();
  TH1F *hFrame2 = gPad->DrawFrame(-1.0, 0.5, 2.5, 2.5*maxLifetime[BG][iFrame][iRap][iPT]);
  sprintf(name, "l_{J/#psi} [mm]");
  hFrame2->SetXTitle(name);

  gLifetime[BG][iFrame][iRap][iPT]->Draw("p same");
//   gLifetime_Fit[BG][P][iFrame][iRap][iPT]->Draw("l same"); //P
//   gLifetime_Fit[BG][NP][iFrame][iRap][iPT]->Draw("l same"); //NP
//   gLifetime_Fit[BG][BKG][iFrame][iRap][iPT]->Draw("l same"); //BG
  gLifetime_Fit[BG][TOT][iFrame][iRap][iPT]->Draw("l same"); //Tot

  if(iRap == 1)
    sprintf(name, "|y| < %1.1f", rapForPTRange[iRap]);
  else 
    sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[iRap-1], rapForPTRange[iRap]);
  TLatex *tex2 = new TLatex(1.75, 0.12*maxLifetime[BG][iFrame][iRap][iPT], name);
  tex2->SetTextSize(0.04); tex2->Draw();
  sprintf(name, "%1.1f < p_{T} < %1.1f GeV", pTRange[iRap][iPT-1], pTRange[iRap][iPT]);
  tex2->DrawLatex(1.4, 0.07*maxLifetime[BG][iFrame][iRap][iPT], name);
  tex2->DrawLatex(1.54, 0.2*maxLifetime[BG][iFrame][iRap][iPT], "mass side bands");

//   sprintf(name, "#chi^{2}/ndf = %1.2f", chi2Lifetime[BG][iFrame][iRap][iPT]);
//   tex2->DrawLatex(-0.9, 1.2*maxLifetime[BG][iFrame][iRap][iPT], name);

//   TLegend *leg2 = new TLegend(0.6752874,0.7394068,0.9798851,0.9385593);
//   leg2->AddEntry(gLifetime_Fit[S][TOT][iFrame][iRap][iPT], "sum", "l");
//   leg2->AddEntry(gLifetime_Fit[S][P][iFrame][iRap][iPT], "prompt J/#psi", "l");
//   leg2->AddEntry(gLifetime_Fit[S][NP][iFrame][iRap][iPT], "non-prompt J/#psi", "l");
//   leg2->AddEntry(gLifetime_Fit[S][BKG][iFrame][iRap][iPT], "background", "l");
//   leg2->SetTextSize(0.04); leg2->SetFillColor(0);
//   leg2->SetBorderSize(0);  leg2->Draw();

  sprintf(name, "Figures/lifetime_BGRegion_%s_rap%d_pT%d.pdf", frameLabel[iFrame], iRap, iPT);
  c2->Print(name);

  //===================================
  //Lifetime in the signal region
  //===================================
  sprintf(name, "lifetime_signal_%s_rap%d_pT%d", frameLabel[iFrame], iRap, iPT);
  TCanvas *c3 = new TCanvas(name, name);
  gPad->SetLogy();
  TH1F *hFrame3 = gPad->DrawFrame(-1.0, 0.5, 2.5, 2.5*maxLifetime[S][iFrame][iRap][iPT]);
  sprintf(name, "l_{J/#psi} [mm]");
  hFrame3->SetXTitle(name);

  gLifetime[S][iFrame][iRap][iPT]->Draw("p same");
  gLifetime_Fit[S][P][iFrame][iRap][iPT]->Draw("l same"); //P
  gLifetime_Fit[S][NP][iFrame][iRap][iPT]->Draw("l same"); //NP
  gLifetime_Fit[S][BKG][iFrame][iRap][iPT]->Draw("l same"); //BG
  gLifetime_Fit[S][TOT][iFrame][iRap][iPT]->Draw("l same"); //Tot

  if(iRap == 1)
    sprintf(name, "|y| < %1.1f", rapForPTRange[iRap]);
  else 
    sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[iRap-1], rapForPTRange[iRap]);
  TLatex *tex3 = new TLatex(1.75, 0.06*maxLifetime[S][iFrame][iRap][iPT], name);
  tex3->SetTextSize(0.04); tex3->Draw();
  sprintf(name, "%1.1f < p_{T} < %1.1f GeV", pTRange[iRap][iPT-1], pTRange[iRap][iPT]);
  tex3->DrawLatex(1.4, 0.03*maxLifetime[S][iFrame][iRap][iPT], name);
  tex3->DrawLatex(1.57, 0.01*maxLifetime[S][iFrame][iRap][iPT], "J/#psi mass region");

//   sprintf(name, "#chi^{2}/ndf = %1.2f", chi2Lifetime[S][iFrame][iRap][iPT]);
//   tex3->DrawLatex(-0.9, 0.7*maxLifetime[S][iFrame][iRap][iPT], name);

  TLegend *leg3 = new TLegend(0.6752874,0.7394068,0.9798851,0.9385593);
  leg3->AddEntry(gLifetime_Fit[S][TOT][iFrame][iRap][iPT], "sum", "l");
  leg3->AddEntry(gLifetime_Fit[S][P][iFrame][iRap][iPT], "prompt J/#psi", "l");
  leg3->AddEntry(gLifetime_Fit[S][NP][iFrame][iRap][iPT], "non-prompt J/#psi", "l");
  leg3->AddEntry(gLifetime_Fit[S][BKG][iFrame][iRap][iPT], "background", "l");
  leg3->SetTextSize(0.04); leg3->SetFillColor(0);
  leg3->SetBorderSize(0);  leg3->Draw();

  sprintf(name, "Figures/lifetime_SignalRegion_%s_rap%d_pT%d.pdf", frameLabel[iFrame], iRap, iPT);
  c3->Print(name);

}

//=======================
void PlotPolGraphs(Int_t iFrame, Int_t iRap, Int_t iPT, Int_t iSpecies){

  Char_t name[100];
  //CosTheta
  sprintf(name, "c1_cosTheta_%s_%s_rap%d_pT%d", speciesLabel[iSpecies], frameLabel[iFrame], iRap, iPT);
  TCanvas *c1 = new TCanvas(name, name);
  TH1F *hFrame1 = gPad->DrawFrame(-1., 0., 1., maximumPol[0][iSpecies][iFrame][iRap][iPT]);
  sprintf(name, "cos#theta_{%s}", frameLabel[iFrame]);
  hFrame1->SetXTitle(name);

  gCosTh[iSpecies][iFrame][iRap][iPT]->Draw("p same");
  gCosTh_Fit[iSpecies][FIT][iFrame][iRap][iPT]->Draw("l same");
  gCosTh_Fit[iSpecies][PONE][iFrame][iRap][iPT]->Draw("l same");
  gCosTh_Fit[iSpecies][ZERO][iFrame][iRap][iPT]->Draw("l same");
  gCosTh_Fit[iSpecies][MONE][iFrame][iRap][iPT]->Draw("l same");

  if(iRap == 1)
    sprintf(name, "|y| < %1.1f", rapForPTRange[iRap]);
  else 
    sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[iRap-1], rapForPTRange[iRap]);
  TLatex *tex1 = new TLatex(-0.93, 0.9*maximumPol[0][iSpecies][iFrame][iRap][iPT], name);
  tex1->SetTextSize(0.04); tex1->Draw();
  sprintf(name, "%1.1f < p_{T} < %1.1f GeV", pTRange[iRap][iPT-1], pTRange[iRap][iPT]);
  tex1->DrawLatex(-0.93, 0.83*maximumPol[0][iSpecies][iFrame][iRap][iPT], name);
  if(iSpecies == P)
    sprintf(name, "l_{J/#psi} < 0.1 mm");
  else if(iSpecies == NP)
    sprintf(name, "l_{J/#psi} > 0.1 mm");
  tex1->DrawLatex(-0.93, 0.76*maximumPol[0][iSpecies][iFrame][iRap][iPT], name);

//   sprintf(name, "#chi^{2}/ndf = %1.2f", chi2[0][iSpecies][iFrame][iRap][iPT]);
//   tex1->DrawLatex(0.6, 0.9*maximumPol[0][iSpecies][iFrame][iRap][iPT], name);

  sprintf(name, "Figures/cosTheta_%s_%s_rap%d_pT%d.pdf", speciesLabel[iSpecies], frameLabel[iFrame], iRap, iPT);
  c1->Print(name);

  //==========================================
  //Phi
  //==========================================
  sprintf(name, "c2_phi_%s_%s_rap%d_pT%d", speciesLabel[iSpecies], frameLabel[iFrame], iRap, iPT);
  TCanvas *c2 = new TCanvas(name, name);
  TH1F *hFrame2 = gPad->DrawFrame(0., 0., 90., maximumPol[1][iSpecies][iFrame][iRap][iPT]);
  sprintf(name, "#phi_{%s}", frameLabel[iFrame]);
  hFrame2->SetXTitle(name);

  gPhi[iSpecies][iFrame][iRap][iPT]->Draw("p same");
  gPhi_Fit[iSpecies][FIT][iFrame][iRap][iPT]->Draw("l same");
  gPhi_Fit[iSpecies][PONE][iFrame][iRap][iPT]->Draw("l same");
  gPhi_Fit[iSpecies][ZERO][iFrame][iRap][iPT]->Draw("l same");
  gPhi_Fit[iSpecies][MONE][iFrame][iRap][iPT]->Draw("l same");

  if(iRap == 1)
    sprintf(name, "|y| < %1.1f", rapForPTRange[iRap]);
  else 
    sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[iRap-1], rapForPTRange[iRap]);
  TLatex *tex2;
  if(iFrame == HX){
    tex2 = new TLatex(3., 0.21*maximumPol[1][iSpecies][iFrame][iRap][iPT], name);
    tex2->SetTextSize(0.04); tex2->Draw();
    sprintf(name, "%1.1f < p_{T} < %1.1f GeV", pTRange[iRap][iPT-1], pTRange[iRap][iPT]);
    tex2->DrawLatex(3., 0.14*maximumPol[1][iSpecies][iFrame][iRap][iPT], name);
    if(iSpecies == P)
      sprintf(name, "l_{J/#psi} < 0.1 mm");
    else if(iSpecies == NP)
      sprintf(name, "l_{J/#psi} > 0.1 mm");
    tex2->DrawLatex(3., 0.07*maximumPol[1][iSpecies][iFrame][iRap][iPT], name);
  }
  else{ //CS
    tex2 = new TLatex(3., 0.9*maximumPol[1][iSpecies][iFrame][iRap][iPT], name);
    tex2->SetTextSize(0.04); tex2->Draw();
    sprintf(name, "%1.1f < p_{T} < %1.1f GeV", pTRange[iRap][iPT-1], pTRange[iRap][iPT]);
    tex2->DrawLatex(3., 0.83*maximumPol[1][iSpecies][iFrame][iRap][iPT], name);
    if(iSpecies == P)
      sprintf(name, "l_{J/#psi} < 0.1 mm");
    else if(iSpecies == NP)
      sprintf(name, "l_{J/#psi} > 0.1 mm");
    tex2->DrawLatex(3., 0.76*maximumPol[1][iSpecies][iFrame][iRap][iPT], name);
  }
//   sprintf(name, "#chi^{2}/ndf = %1.2f", chi2[1][iSpecies][iFrame][iRap][iPT]);
//   tex2->DrawLatex(70., 0.9*maximumPol[1][iSpecies][iFrame][iRap][iPT], name);

  sprintf(name, "Figures/phi_%s_%s_rap%d_pT%d.pdf", speciesLabel[iSpecies], frameLabel[iFrame], iRap, iPT);
  c2->Print(name);
}

//===================================================================
void PlotPolGraphsLRSidebands(Int_t iFrame, Int_t iRap, Int_t iPT, Int_t iSideBand){

  Char_t name[100];
  //CosTheta
  sprintf(name, "c1_cosTheta_%s_%s_rap%d_pT%d", BGName[iSideBand], frameLabel[iFrame], iRap, iPT);
  TCanvas *c1 = new TCanvas(name, name);
  TH1F *hFrame1 = gPad->DrawFrame(-1., 0., 1., maximumPolBG[0][iSideBand][iFrame][iRap][iPT]);
  sprintf(name, "cos#theta_{%s}", frameLabel[iFrame]);
  hFrame1->SetXTitle(name);

  gCosThBG[iSideBand][iFrame][iRap][iPT]->Draw("p same");
//   hCosThBG[iSideBand][iFrame][iRap][iPT]->Draw("p same");
  gCosTh_FitBG[iSideBand][PONE-1][iFrame][iRap][iPT]->Draw("l same");
  gCosTh_FitBG[iSideBand][ZERO-1][iFrame][iRap][iPT]->Draw("l same");
  gCosTh_FitBG[iSideBand][MONE-1][iFrame][iRap][iPT]->Draw("l same");

  if(iRap == 1)
    sprintf(name, "|y| < %1.1f", rapForPTRange[iRap]);
  else 
    sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[iRap-1], rapForPTRange[iRap]);
  TLatex *tex1 = new TLatex(-0.93, 0.9*maximumPolBG[0][iSideBand][iFrame][iRap][iPT], name);
  tex1->SetTextSize(0.04); tex1->Draw();
  sprintf(name, "%1.1f < p_{T} < %1.1f GeV", pTRange[iRap][iPT-1], pTRange[iRap][iPT]);
  tex1->DrawLatex(-0.93, 0.83*maximumPolBG[0][iSideBand][iFrame][iRap][iPT], name);
  if(iSideBand == L)
    sprintf(name, "L mass sideband");
  else if(iSideBand == R)
    sprintf(name, "R mass sideband");
  tex1->DrawLatex(-0.93, 0.76*maximumPolBG[0][iSideBand][iFrame][iRap][iPT], name);

//   sprintf(name, "#chi^{2}/ndf = %1.2f", chi2[0][iSideBand][iFrame][iRap][iPT]);
//   tex1->DrawLatex(0.6, 0.9*maximumPolBG[0][iSideBand][iFrame][iRap][iPT], name);

  sprintf(name, "Figures/cosTheta_%s_%s_rap%d_pT%d.pdf", BGName[iSideBand], frameLabel[iFrame], iRap, iPT);
  c1->Print(name);

  //==========================================
  //Phi
  //==========================================
  sprintf(name, "c2_phi_%s_%s_rap%d_pT%d", BGName[iSideBand], frameLabel[iFrame], iRap, iPT);
  TCanvas *c2 = new TCanvas(name, name);
  TH1F *hFrame2 = gPad->DrawFrame(0., 0., 90., maximumPolBG[1][iSideBand][iFrame][iRap][iPT]);
  sprintf(name, "#phi_{%s}", frameLabel[iFrame]);
  hFrame2->SetXTitle(name);

  gPhiBG[iSideBand][iFrame][iRap][iPT]->Draw("p same");
  //hPhiBG[iSideBand][iFrame][iRap][iPT]->Draw("p same");
  gPhi_FitBG[iSideBand][PONE-1][iFrame][iRap][iPT]->Draw("l same");
  gPhi_FitBG[iSideBand][ZERO-1][iFrame][iRap][iPT]->Draw("l same");
  gPhi_FitBG[iSideBand][MONE-1][iFrame][iRap][iPT]->Draw("l same");

  if(iRap == 1)
    sprintf(name, "|y| < %1.1f", rapForPTRange[iRap]);
  else 
    sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[iRap-1], rapForPTRange[iRap]);
  TLatex *tex2;
//   if(iFrame == HX){
//     tex2 = new TLatex(3., 0.21*maximumPolBG[1][iSideBand][iFrame][iRap][iPT], name);
//     tex2->SetTextSize(0.04); tex2->Draw();
//     sprintf(name, "%1.1f < p_{T} < %1.1f GeV", pTRange[iRap][iPT-1], pTRange[iRap][iPT]);
//     tex2->DrawLatex(3., 0.14*maximumPolBG[1][iSideBand][iFrame][iRap][iPT], name);
//     if(iSideBand == P)
//       sprintf(name, "l_{J/#psi} < 0.1 mm");
//     else if(iSideBand == NP)
//       sprintf(name, "l_{J/#psi} > 0.1 mm");
//     tex2->DrawLatex(3., 0.07*maximumPolBG[1][iSideBand][iFrame][iRap][iPT], name);
//   }
//   else{ //CS
    tex2 = new TLatex(3., 0.9*maximumPolBG[1][iSideBand][iFrame][iRap][iPT], name);
    tex2->SetTextSize(0.04); tex2->Draw();
    sprintf(name, "%1.1f < p_{T} < %1.1f GeV", pTRange[iRap][iPT-1], pTRange[iRap][iPT]);
    tex2->DrawLatex(3., 0.83*maximumPolBG[1][iSideBand][iFrame][iRap][iPT], name);
    if(iSideBand == L)
      sprintf(name, "L mass sideband");
    else if(iSideBand == R)
      sprintf(name, "R mass sideband");
    tex2->DrawLatex(3., 0.76*maximumPolBG[1][iSideBand][iFrame][iRap][iPT], name);
//   }

  sprintf(name, "Figures/phi_%s_%s_rap%d_pT%d.pdf", BGName[iSideBand], frameLabel[iFrame], iRap, iPT);
  c2->Print(name);
}

//=======================
void GetPolGraphs(Int_t iFrame, Int_t iRap, Int_t iPT){

  Char_t name[100];
  sprintf(name, "/home/hermine/CMS/Work/Polarization/PlotsForPaper/RootFiles/pedagogical/fitFrame%s_%d_%d-smoothed.root", frameLabel[iFrame], iRap, iPT, frameLabel[iFrame]);
  TFile *f = new TFile(name);

  RooPlot *tempCosTh[2], *tempPhi[2];

  //theta for prompt
//   sprintf(name, "polCosTh%sPrompt_%d_%d", frameLabel[iFrame], iRap, iPT);
  sprintf(name, "polCosTh%sPrompt_peda_%d_%d", frameLabel[iFrame], iRap, iPT);
  tempCosTh[P] = (RooPlot *) gDirectory->Get(name);
  gCosTh[P][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[P]->getObject(0); //Data
  gCosTh_Fit[P][FIT][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[P]->getObject(1); //Fit
  gCosTh_Fit[P][PONE][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[P]->getObject(2); //lambda_theta = +1
  gCosTh_Fit[P][ZERO][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[P]->getObject(3); //lambda_theta = 0
  gCosTh_Fit[P][MONE][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[P]->getObject(4); //lambda_theta = -1
  gCosTh_Fit[P][PONE][iFrame][iRap][iPT]->SetLineWidth(2);
  gCosTh_Fit[P][ZERO][iFrame][iRap][iPT]->SetLineWidth(2);
  gCosTh_Fit[P][MONE][iFrame][iRap][iPT]->SetLineWidth(2);

  //theta for non-prompt
  sprintf(name, "polCosTh%sNonPrompt_peda_%d_%d", frameLabel[iFrame], iRap, iPT);
  tempCosTh[NP] = (RooPlot *) gDirectory->Get(name);
  gCosTh[NP][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[NP]->getObject(0); //Data
  gCosTh_Fit[NP][FIT][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[NP]->getObject(1); //Fit
  gCosTh_Fit[NP][PONE][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[NP]->getObject(2); //lambda_theta = +1
  gCosTh_Fit[NP][ZERO][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[NP]->getObject(3); //lambda_theta = 0
  gCosTh_Fit[NP][MONE][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[NP]->getObject(4); //lambda_theta = -1
  gCosTh_Fit[NP][PONE][iFrame][iRap][iPT]->SetLineWidth(2);
  gCosTh_Fit[NP][ZERO][iFrame][iRap][iPT]->SetLineWidth(2);
  gCosTh_Fit[NP][MONE][iFrame][iRap][iPT]->SetLineWidth(2);

  //
  //phi for prompt
  sprintf(name, "polPhi%sPrompt_peda_%d_%d", frameLabel[iFrame], iRap, iPT);
  tempPhi[P] = (RooPlot *) gDirectory->Get(name);
  gPhi[P][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[P]->getObject(0); //Data
  gPhi_Fit[P][FIT][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[P]->getObject(1); //Fit
  gPhi_Fit[P][PONE][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[P]->getObject(2); //lambda_theta = +1
  gPhi_Fit[P][ZERO][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[P]->getObject(3); //lambda_theta = 0
  gPhi_Fit[P][MONE][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[P]->getObject(4); //lambda_theta = -1
  gPhi_Fit[P][PONE][iFrame][iRap][iPT]->SetLineWidth(2);
  gPhi_Fit[P][ZERO][iFrame][iRap][iPT]->SetLineWidth(2);
  gPhi_Fit[P][MONE][iFrame][iRap][iPT]->SetLineWidth(2);

  //phi for non-prompt
  sprintf(name, "polPhi%sNonPrompt_peda_%d_%d", frameLabel[iFrame], iRap, iPT);
  tempPhi[NP] = (RooPlot *) gDirectory->Get(name);
  gPhi[NP][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[NP]->getObject(0); //Data
  gPhi_Fit[NP][FIT][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[NP]->getObject(1); //Fit
  gPhi_Fit[NP][PONE][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[NP]->getObject(2); //lambda_theta = +1
  gPhi_Fit[NP][ZERO][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[NP]->getObject(3); //lambda_theta = 0
  gPhi_Fit[NP][MONE][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[NP]->getObject(4); //lambda_theta = -1
  gPhi_Fit[NP][PONE][iFrame][iRap][iPT]->SetLineWidth(2);
  gPhi_Fit[NP][ZERO][iFrame][iRap][iPT]->SetLineWidth(2);
  gPhi_Fit[NP][MONE][iFrame][iRap][iPT]->SetLineWidth(2);

  //
  maximumPol[0][P][iFrame][iRap][iPT] = tempCosTh[P]->GetMaximum();
  chi2[0][P][iFrame][iRap][iPT] = tempCosTh[P]->chiSquare(tempCosTh[P]->nameOf(1), tempCosTh[P]->nameOf(0));
  printf("%s, rap %d, pT %d --> chi2/ndf (cosTheta fit P) = %1.3f\n", frameLabel[iFrame], iRap, iPT, chi2[0][P][iFrame][iRap][iPT]);

  maximumPol[0][NP][iFrame][iRap][iPT] = tempCosTh[NP]->GetMaximum();
  chi2[0][NP][iFrame][iRap][iPT] = tempCosTh[NP]->chiSquare(tempCosTh[NP]->nameOf(1), tempCosTh[NP]->nameOf(0));
  printf("%s, rap %d, pT %d --> chi2/ndf (cosTheta fit NP) = %1.3f\n", frameLabel[iFrame], iRap, iPT, chi2[0][NP][iFrame][iRap][iPT]);
  //
  maximumPol[1][P][iFrame][iRap][iPT] = tempPhi[P]->GetMaximum();
  chi2[1][P][iFrame][iRap][iPT] = tempPhi[P]->chiSquare(tempPhi[P]->nameOf(1), tempPhi[P]->nameOf(0));
  printf("%s, rap %d, pT %d --> chi2/ndf (Phi fit P) = %1.3f\n", frameLabel[iFrame], iRap, iPT, chi2[1][P][iFrame][iRap][iPT]);

  maximumPol[1][NP][iFrame][iRap][iPT] = tempPhi[NP]->GetMaximum();
  chi2[1][NP][iFrame][iRap][iPT] = tempPhi[NP]->chiSquare(tempPhi[NP]->nameOf(1), tempPhi[NP]->nameOf(0));
  printf("%s, rap %d, pT %d --> chi2/ndf (Phi fit NP) = %1.3f\n", frameLabel[iFrame], iRap, iPT, chi2[1][NP][iFrame][iRap][iPT]);
}

//================================================
void GetPolGraphsLRSidebands(Int_t iFrame, Int_t iRap, Int_t iPT){

  Char_t name[100];
  sprintf(name, "/home/hermine/CMS/Work/Polarization/PlotsForPaper/RootFiles/pedagogical/fitFrame%s_%d_%d-smoothed.root", frameLabel[iFrame], iRap, iPT);
  TFile *f = new TFile(name);

  RooPlot *tempCosTh[2], *tempPhi[2];

  //theta for L
  sprintf(name, "polCosTh%s%s_peda_%d_%d", frameLabel[iFrame], BGName[L], iRap, iPT);
  tempCosTh[L] = (RooPlot *) gDirectory->Get(name);
  gCosThBG[L][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[L]->getObject(0); //Data
  hCosThBG[L][iFrame][iRap][iPT] = (TH1F *) tempCosTh[L]->getHist(tempCosTh[L]->nameOf(0));
  hCosThBG[L][iFrame][iRap][iPT]->SetMarkerStyle(20);
  hCosThBG[L][iFrame][iRap][iPT]->Rebin(6);
  gCosTh_FitBG[L][PONE-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[L]->getObject(1); //lambda_theta = +1
  gCosTh_FitBG[L][ZERO-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[L]->getObject(2); //lambda_theta = 0
  gCosTh_FitBG[L][MONE-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[L]->getObject(3); //lambda_theta = -1
  gCosTh_FitBG[L][PONE-1][iFrame][iRap][iPT]->SetLineWidth(2);
  gCosTh_FitBG[L][ZERO-1][iFrame][iRap][iPT]->SetLineWidth(2);
  gCosTh_FitBG[L][MONE-1][iFrame][iRap][iPT]->SetLineWidth(2);

  //theta for R
  sprintf(name, "polCosTh%s%s_peda_%d_%d", frameLabel[iFrame], BGName[R], iRap, iPT);
  tempCosTh[R] = (RooPlot *) gDirectory->Get(name);
  gCosThBG[R][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[R]->getObject(0); //Data
  hCosThBG[R][iFrame][iRap][iPT] = (TH1F *) tempCosTh[R]->getHist(tempCosTh[R]->nameOf(0));
  hCosThBG[R][iFrame][iRap][iPT]->Rebin(4);
  gCosTh_FitBG[R][PONE-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[R]->getObject(1); //lambda_theta = +1
  gCosTh_FitBG[R][ZERO-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[R]->getObject(2); //lambda_theta = 0
  gCosTh_FitBG[R][MONE-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempCosTh[R]->getObject(3); //lambda_theta = -1
  gCosTh_FitBG[R][PONE-1][iFrame][iRap][iPT]->SetLineWidth(2);
  gCosTh_FitBG[R][ZERO-1][iFrame][iRap][iPT]->SetLineWidth(2);
  gCosTh_FitBG[R][MONE-1][iFrame][iRap][iPT]->SetLineWidth(2);
  //
  //phi for L
  sprintf(name, "polPhi%s%s_peda_%d_%d", frameLabel[iFrame], BGName[L], iRap, iPT);
  tempPhi[L] = (RooPlot *) gDirectory->Get(name);
  gPhiBG[L][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[L]->getObject(0); //Data
  hPhiBG[L][iFrame][iRap][iPT] = (TH1F *) tempPhi[L]->getHist(tempPhi[L]->nameOf(0));
  hPhiBG[L][iFrame][iRap][iPT]->Rebin(4);
  gPhi_FitBG[L][PONE-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[L]->getObject(1); //lambda_theta = +1
  gPhi_FitBG[L][ZERO-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[L]->getObject(2); //lambda_theta = 0
  gPhi_FitBG[L][MONE-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[L]->getObject(3); //lambda_theta = -1
  gPhi_FitBG[L][PONE-1][iFrame][iRap][iPT]->SetLineWidth(2);
  gPhi_FitBG[L][ZERO-1][iFrame][iRap][iPT]->SetLineWidth(2);
  gPhi_FitBG[L][MONE-1][iFrame][iRap][iPT]->SetLineWidth(2);

  //phi for R
  sprintf(name, "polPhi%s%s_peda_%d_%d", frameLabel[iFrame], BGName[R], iRap, iPT);
  tempPhi[R] = (RooPlot *) gDirectory->Get(name);
  gPhiBG[R][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[R]->getObject(0); //Data
  hPhiBG[R][iFrame][iRap][iPT] = (TH1F *) tempPhi[R]->getHist(tempPhi[R]->nameOf(0));
  hPhiBG[R][iFrame][iRap][iPT]->Rebin(2);
  gPhi_FitBG[R][PONE-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[R]->getObject(1); //lambda_theta = +1
  gPhi_FitBG[R][ZERO-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[R]->getObject(2); //lambda_theta = 0
  gPhi_FitBG[R][MONE-1][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempPhi[R]->getObject(3); //lambda_theta = -1
  gPhi_FitBG[R][PONE-1][iFrame][iRap][iPT]->SetLineWidth(2);
  gPhi_FitBG[R][ZERO-1][iFrame][iRap][iPT]->SetLineWidth(2);
  gPhi_FitBG[R][MONE-1][iFrame][iRap][iPT]->SetLineWidth(2);

  //
  maximumPolBG[0][L][iFrame][iRap][iPT] = tempCosTh[L]->GetMaximum();
  maximumPolBG[0][R][iFrame][iRap][iPT] = tempCosTh[R]->GetMaximum();
  maximumPolBG[1][L][iFrame][iRap][iPT] = tempPhi[L]->GetMaximum();
  maximumPolBG[1][R][iFrame][iRap][iPT] = tempPhi[R]->GetMaximum();
}

//=======================
void GetMassLifetime(Int_t iFrame, Int_t iRap, Int_t iPT){

  Char_t name[100];
  sprintf(name, "/home/hermine/CMS/Work/Polarization/PlotsForPaper/RootFiles/species/fitFrame%s_%d_%d-%s.root", frameLabel[iFrame], iRap, iPT, frameLabel[iFrame]);
  TFile *f = new TFile(name);

  RooPlot *tempMass, *tempLifetime;

  //mass
  sprintf(name, "mass_plot_rap%d_pt%d", iRap, iPT);
  tempMass = (RooPlot *) gDirectory->Get(name);
  gMass[iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempMass->getObject(0); //Data
  gMass_Fit[TOT][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempMass->getObject(1); //Fit-Total
  gMass_Fit[P][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempMass->getObject(4); //Fit-P
  gMass_Fit[NP][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempMass->getObject(3); //Fit-NP
  gMass_Fit[BKG][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempMass->getObject(2); //Fit-BG
  maxMass[iFrame][iRap][iPT] = tempMass->GetMaximum();
  chi2Mass[iFrame][iRap][iPT] = tempMass->chiSquare(tempMass->nameOf(1), tempMass->nameOf(0));
  printf("%s, rap %d, pT %d --> chi2/ndf (mass fit) = %1.3f\n", frameLabel[iFrame], iRap, iPT, chi2Mass[iFrame][iRap][iPT]);

  //lifetime - signal region
  sprintf(name, "ctausig_plot_rap%d_pt%d", iRap, iPT);
  tempLifetime = (RooPlot *) gDirectory->Get(name);
  gLifetime[S][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempLifetime->getObject(0); //Data
  gLifetime_Fit[S][TOT][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempLifetime->getObject(1); //Fit-Total
  gLifetime_Fit[S][P][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempLifetime->getObject(2); //Fit-P
  gLifetime_Fit[S][NP][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempLifetime->getObject(3); //Fit-NP
  gLifetime_Fit[S][BKG][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempLifetime->getObject(4); //Fit-BG
  maxLifetime[S][iFrame][iRap][iPT] = tempLifetime->GetMaximum();
  chi2Lifetime[S][iFrame][iRap][iPT] = tempLifetime->chiSquare(tempLifetime->nameOf(1), tempLifetime->nameOf(0));
  printf("%s, rap %d, pT %d --> chi2/ndf (lifetime-signal fit) = %1.3f\n", frameLabel[iFrame], iRap, iPT, chi2Lifetime[S][iFrame][iRap][iPT]);

  //lifetime - BG region
  sprintf(name, "ctaubkg_plot_rap%d_pt%d", iRap, iPT);
  tempLifetime = (RooPlot *) gDirectory->Get(name);
  gLifetime[BG][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempLifetime->getObject(0); //Data
  gLifetime_Fit[BG][TOT][iFrame][iRap][iPT] = (TGraphAsymmErrors *) tempLifetime->getObject(1); //Fit-Total
  maxLifetime[BG][iFrame][iRap][iPT] = tempLifetime->GetMaximum();
  chi2Lifetime[BG][iFrame][iRap][iPT] = tempLifetime->chiSquare(tempLifetime->nameOf(1), tempLifetime->nameOf(0));
  printf("%s, rap %d, pT %d --> chi2/ndf (lifetime-BG fit) = %1.3f\n", frameLabel[iFrame], iRap, iPT, chi2Lifetime[BG][iFrame][iRap][iPT]);
  //
}

//=======================
void GetContour(Int_t iFrame, Int_t iRap, Int_t iPT){

  Char_t name[100];
  sprintf(name, "/home/hermine/CMS/Work/Polarization/PlotsForPaper/RootFiles/species/fitFrame%s_%d_%d-%s.root", frameLabel[iFrame], iRap, iPT, frameLabel[iFrame]);
  TFile *f = new TFile(name);

  RooPlot *contour[kNbVarComb][kNbSpecies];
  //1.) PROMPT J/PSI'S
  //theta vs phi
  sprintf(name, "th_ph_p_rap%d_pt%d", iRap, iPT);
  contour[TH_PHI][P] = (RooPlot *) gDirectory->Get(name);
  gContour[TH_PHI][P][iFrame][iRap][iPT] = (TGraphAsymmErrors *) contour[TH_PHI][P]->getObject(0); //Ellipse
  gContour_line1[TH_PHI][P][iFrame][iRap][iPT] = (TLine *) contour[TH_PHI][P]->getObject(1); //line
  gContour_line2[TH_PHI][P][iFrame][iRap][iPT] = (TLine *) contour[TH_PHI][P]->getObject(2); //line
  gContour_marker[TH_PHI][P][iFrame][iRap][iPT] = (TLine *) contour[TH_PHI][P]->getObject(3); //marker
  //theta vs thetaPhi
  sprintf(name, "th_thph_p_rap%d_pt%d", iRap, iPT);
  contour[TH_THETAPHI][P] = (RooPlot *) gDirectory->Get(name);
  gContour[TH_THETAPHI][P][iFrame][iRap][iPT] = (TGraphAsymmErrors *) contour[TH_THETAPHI][P]->getObject(0); //Ellipse
  gContour_line1[TH_THETAPHI][P][iFrame][iRap][iPT] = (TLine *) contour[TH_THETAPHI][P]->getObject(1); //line
  gContour_line2[TH_THETAPHI][P][iFrame][iRap][iPT] = (TLine *) contour[TH_THETAPHI][P]->getObject(2); //line
  gContour_marker[TH_THETAPHI][P][iFrame][iRap][iPT] = (TLine *) contour[TH_THETAPHI][P]->getObject(3); //marker
  //thetaPhi vs phi
  sprintf(name, "thph_ph_p_rap%d_pt%d", iRap, iPT);
  contour[THETAPHI_PHI][P] = (RooPlot *) gDirectory->Get(name);
  gContour[THETAPHI_PHI][P][iFrame][iRap][iPT] = (TGraphAsymmErrors *) contour[THETAPHI_PHI][P]->getObject(0); //Ellipse
  gContour_line1[THETAPHI_PHI][P][iFrame][iRap][iPT] = (TLine *) contour[THETAPHI_PHI][P]->getObject(1); //line
  gContour_line2[THETAPHI_PHI][P][iFrame][iRap][iPT] = (TLine *) contour[THETAPHI_PHI][P]->getObject(2); //line
  gContour_marker[THETAPHI_PHI][P][iFrame][iRap][iPT] = (TLine *) contour[THETAPHI_PHI][P]->getObject(3); //marker


  //2.) NON-PROMPT J/PSI'S
  //theta vs phi
  sprintf(name, "th_ph_np_rap%d_pt%d", iRap, iPT);
  contour[TH_PHI][NP] = (RooPlot *) gDirectory->Get(name);
  gContour[TH_PHI][NP][iFrame][iRap][iPT] = (TGraphAsymmErrors *) contour[TH_PHI][NP]->getObject(0); //Ellipse
  gContour_line1[TH_PHI][NP][iFrame][iRap][iPT] = (TLine *) contour[TH_PHI][NP]->getObject(1); //line
  gContour_line2[TH_PHI][NP][iFrame][iRap][iPT] = (TLine *) contour[TH_PHI][NP]->getObject(2); //line
  gContour_marker[TH_PHI][NP][iFrame][iRap][iPT] = (TLine *) contour[TH_PHI][NP]->getObject(3); //marker
  //theta vs thetaPhi
  sprintf(name, "th_thph_np_rap%d_pt%d", iRap, iPT);
  contour[TH_THETAPHI][NP] = (RooPlot *) gDirectory->Get(name);
  gContour[TH_THETAPHI][NP][iFrame][iRap][iPT] = (TGraphAsymmErrors *) contour[TH_THETAPHI][NP]->getObject(0); //Ellipse
  gContour_line1[TH_THETAPHI][NP][iFrame][iRap][iPT] = (TLine *) contour[TH_THETAPHI][NP]->getObject(1); //line
  gContour_line2[TH_THETAPHI][NP][iFrame][iRap][iPT] = (TLine *) contour[TH_THETAPHI][NP]->getObject(2); //line
  gContour_marker[TH_THETAPHI][NP][iFrame][iRap][iPT] = (TLine *) contour[TH_THETAPHI][NP]->getObject(3); //marker
  //thetaPhi vs phi
  sprintf(name, "th_thph_np_rap%d_pt%d", iRap, iPT);
  contour[THETAPHI_PHI][NP] = (RooPlot *) gDirectory->Get(name);
  gContour[THETAPHI_PHI][NP][iFrame][iRap][iPT] = (TGraphAsymmErrors *) contour[THETAPHI_PHI][NP]->getObject(0); //Ellipse
  gContour_line1[THETAPHI_PHI][NP][iFrame][iRap][iPT] = (TLine *) contour[THETAPHI_PHI][NP]->getObject(1); //line
  gContour_line2[THETAPHI_PHI][NP][iFrame][iRap][iPT] = (TLine *) contour[THETAPHI_PHI][NP]->getObject(2); //line
  gContour_marker[THETAPHI_PHI][NP][iFrame][iRap][iPT] = (TLine *) contour[THETAPHI_PHI][NP]->getObject(3); //marker
}
