#include "../interface/rootIncludes.inc"
#include "TGraphAsymmErrors.h"
#include "TMath.h"

Int_t const kNbEff = 9;
Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/LucaPerrozzi/17March2011/MuonTrackingEff_17March2011.root",

				"/Users/hwoehri/CMS/Work/TnP/Zongchang/29March2011/MuonIDEff_29March2011_fitted.root",

				"/Users/hwoehri/CMS/Work/TnP/Xianyou/25March2011/MuonQualEff_25March2011_fitted.root",

				"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L1L2EfficiencyMu0_Mu3_Track3_pt_abseta_seagull_looseCuts_runA.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L1L2EfficiencyMu0_Mu3_Track3_pt_abseta_seagull_looseCuts_runB.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L1L2EfficiencyMu0_Mu3_Track3_pt_abseta_seagull_tightCuts_runA.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L1L2EfficiencyMu0_Mu3_Track3_pt_abseta_seagull_tightCuts_runB.root",

				"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L3DoubleMu0_seagull_looseCuts_runA.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L3DoubleMu0_seagull_looseCuts_runB.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L3DoubleMu0_seagull_tightCuts_runA.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/22June2011/L3DoubleMu0_seagull_tightCuts_runB.root",

				"/Users/hwoehri/CMS/Work/TnP/Herbert/22June2011/MuX_L2Mu0_3etabinstight.root", //

				//"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_looseCuts_runA.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_looseCuts_runB1.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_looseCuts_runB2.root",
				"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_tightCuts_runA.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_tightCuts_runB1.root",
				//"/Users/hwoehri/CMS/Work/TnP/Ilse/20June2011/TkMu_Mu3_Track3_pt_abseta_seagull_tightCuts_runB2.root",

				"/Users/hwoehri/CMS/Work/TnP/Herbert/24May2011/MuX_L2Mu0_L3Mu0.root",

				"/Users/hwoehri/CMS/Work/TnP/Ilse/23May2011/TkMu0L1L2-pt-abseta-runA-distM2gt120_modified.root"};
//old version of scripts:
//  				"/Users/hwoehri/CMS/Work/TnP/Luigi/25March2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_25March2011_fitted.root",
//				"/Users/hwoehri/CMS/Work/TnP/Luigi/5April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010A_distM1gt150_7April2011_fitted.root",
//				"/Users/hwoehri/CMS/Work/TnP/Luigi/5April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM1gt150_7April2011_fitted.root",
				//"/Users/hwoehri/CMS/Work/TnP/Francesco/2April2011/L1L2_DMu0_TriggerEfficiencies_RUNB_distM2gt120_7April2011_fitted.root",
  				//"/Users/hwoehri/CMS/Work/TnP/Francesco/2April2011/L1L2_DMu0_TriggerEfficiencies_RUNA_distM2gt120_7April2011_fitted.root",
//				"/Users/hwoehri/CMS/Work/TnP/Luigi/17April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010A_distM2gt120_17April2011.root",
//				"/Users/hwoehri/CMS/Work/TnP/Luigi/17April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM2gt120_17April2011.root",
//				"/Users/hwoehri/CMS/Work/TnP/Luigi/17April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010A_distM2gt120_17April2011_fitted.root",
//				"/Users/hwoehri/CMS/Work/TnP/Luigi/17April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM2gt120_17April2011_fitted.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Herbert/28March2011/HLT_Track_Mu0TkMu0_Run1.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Herbert/29March2011/HLT_Track_Mu0TkMu0_Run_Tightv3.root",
				//"/Users/hwoehri/CMS/Work/TnP/Herbert/13May2011/HLT_Track_Mu0TkMu0_Run1_13May2011.root", //
//  				"/Users/hwoehri/CMS/Work/TnP/Francesco/23March2011/L1L2_DMu0_TriggerEfficiencies_23March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Zongchang/19March2011/MuonIDEff_19March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Francesco/18March2011/L1L2_DoubleMu0_TriggerEfficiencies_18March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Luigi/25March2011/L3_DoubleMu0_TriggerEfficiencies_RunABmixed_25March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Luigi/19March2011/L3_DoubleMu0_TriggerEfficiencies_19March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Xianyou/19March2011/MuonQualEff_19March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Ilse/24March2011/HLTMuonTrack_Mu0_TkMu0_TM.root"};
// 				"/Users/hwoehri/CMS/Work/TnP/Ilse/29March2011/HLTMuonTrack_Mu0_TkMu0_TM_runA_29March2011.root"};
//				"/Users/hwoehri/CMS/Work/TnP/Ilse/29March2011/HLTMuonTrack_Mu0_TkMu0_TM_runB1_29March2011.root",

enum {TrkEff, MuIDEff, MuQualEff, L1L2Eff, L3Eff, TkMu_Track, TkMu_Mu, TkMuL1L2, TkMuL3};
Char_t *effName[kNbEff] = {"TrkEff", "MuIDEff", "MuQualEff", "L1L2Eff", "L3Eff", "Trk_TkMu0Eff", "Mu_TkMu0Eff", "Track_TkMu0AndL3Mu0", "Mu_TkMu0andL1L2"};
Char_t *effNameLong[kNbEff] = {"offline tracking efficiency", "muon identification efficiency", 
				"efficiency of muon quality cuts", "L1-L2 trigger efficiency (Mu0)", 
			       "L3 trigger efficiency (Mu0)", "tracking efficiency of TkMu0",
			       "muon trigger efficency of TkMu0",
			       "tracking efficiency of TkMu0 and L3 of Mu0",
			       "muon trigger efficiency of TkMu0 and L1-L2 combined"};
Int_t const kNbEffSample = 3;
enum {DATA, MC, MCTRUTH};
Char_t *effSampleName[kNbEffSample] = {"DATA", "MC", "MCTRUTH"};
//
enum {CENTRAL, UPPER, LOWER};
Char_t *valName[3] = {"central", "lower", "upper"};
TH2D *hMuEff[kNbEff][kNbEffSample][3];
Int_t const kNbEta = 8;
//Int_t const kNbEta = 3; //tracking eff of TkMu
//Int_t const kNbEta = 14; //new L1/L2 and L3 efficiencies
//Int_t const kNbPT = 12; //"old" PT binning
//Int_t const kNbPT = 11; //"new" PT binning
//Int_t const kNbPT = 15; //"new" L1/L2 and L3 PT binning
//Int_t const kNbPT = 11; //TkMu_Mu, 20June 2011: loose cuts
//Int_t const kNbPT = 9; //TkMu_Mu, 20June 2011: tight cuts
Int_t const kNbPT = 8; //L1*L2, 20June 2011: loose/tight cuts
TGraphAsymmErrors *gEff_pT[kNbEff][kNbEffSample][kNbEta];
TGraphAsymmErrors *gEff_eta[kNbEff][kNbEffSample][kNbPT];
//
TF1 *fitPT, *fitPT_PlateauCorrected;

Int_t const marker[kNbEffSample] = {20, 25, 21};
// Float_t binsEta[kNbEta] = {1.2, 1.6,  2.4}; //Track-TkMu0
// Float_t binsEta[kNbEta] = {0.2, 0.3, 0.6, 0.8, 1.2, 1.6, 2.1, 2.4};
//Float_t binsEta[kNbEta] = {0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4};//new L1/L2 and L3 binning
//Float_t binsPT[kNbPT+1] = {0.7, 1., 1.25, 1.5, 1.75, 2., 2.5, 3.5, 4., 6., 10., 20., 30.}; //old binning
// Float_t binsPT[kNbPT+1] = {0.8, 1.25, 1.5, 1.75, 2., 2.5, 3.3, 4., 6., 10., 20., 30.};//standard binning
//Float_t binsPT[kNbPT+1] = {0.8, 1.25, 1.5, 1.75, 2., 2.25, 2.5, 3, 3.3,  4.,  5., 6., 7., 10., 20., 30.};//new L1/L2 and L3 binning
//June 2011:
//Float_t binsEta[kNbEta] = {0.2, 0.3, 0.6, 0.8, 1.1, 1.65, 2.1, 2.4};//TkMu_Mu, 20 June 2011: loose cuts
//Float_t binsPT[kNbPT+1] = {1.2, 1.5, 1.75, 2., 2.5, 3., 3.5, 4., 6., 10., 20., 30.};//TkMu_Mu, 20 June 2011: loose cuts
//Float_t binsEta[kNbEta] = {0.2, 0.3, 0.6, 0.8, 1.2, 1.6, 2.1, 2.25};//TkMu_Mu, 20 June 2011: tight cuts
//Float_t binsPT[kNbPT+1] = {1.75, 2., 2.5, 3., 3.8, 4.5, 6., 10., 20., 30.};//TkMu, 20 June 2011: tight cuts

Float_t binsEta[kNbEta] = {0.2, 0.3, 0.6, 0.8, 1.1, 1.4, 2.1, 2.4};//L1*L2, L3, 22 June: loose cuts
Float_t binsPT[kNbPT+1] = {2., 2.75, 3., 4.0, 4.6, 6., 10., 20., 30.};//L1*L2, L3, 20 June 2011: loose cuts
// Float_t binsEta[kNbEta] = {0.2, 0.3, 0.6, 0.8, 1.2, 1.3, 1.6, 2.2};//L1*L2, L3, 22 June: tight cuts
// Float_t binsPT[kNbPT+1] = {3., 3.3, 4.0, 4.7, 5.2, 6., 10., 20., 30.};//L1*L2, L3, 20 June 2011: tight cuts

Double_t kNbPTMax = 30.; //upper edge of calculations (needed to constructe the large fitting histogram)
TGraphAsymmErrors *gEffPT_allEta;
Double_t pTMin = 3.; //fit efficiency curves only from "pTMin" onwards
// Int_t etaIndicesToMerge[] = {0, 2, 3, 4};//MC
//Int_t etaIndicesToMerge[] = {1};//MC
Int_t etaIndicesToMerge[] = {5, 6, 7};//MC
void LoadEfficiencies(Int_t iEff, Int_t iEffSample);
void PlotEff_DataVsMC_PT(Int_t iEff, Int_t iEta);
void PlotEff_DataVsMC_Eta(Int_t iEff, Int_t iPT);
Double_t erfunction(Double_t *x, Double_t *par);
Double_t erfunction2(Double_t *x, Double_t *par);
Double_t erfunction_corr(Double_t *x, Double_t *par);
void ParametrizePT(Int_t iEff, Int_t iEffSample, Int_t iEta, Bool_t saveFit, Bool_t plateauCorrection);
void ConcatenateGraphs(Int_t iEff, Int_t iEffSample, Int_t nBinsToMerge, Double_t pTMin);
Double_t fitConcatenatedGraph(Double_t *x, Double_t *par);
void FitConcatenatedGraph(Int_t nBins, Int_t iEff, Int_t iEffSample);
Double_t calcChi2(TGraphAsymmErrors *graph, TF1 *func, Int_t fitPar);
//==============================================================
void plotSingleMuEff(Bool_t saveFit = kFALSE,
		     Bool_t plateauCorrection = kFALSE){

  //for(int iEff = 0; iEff < kNbEff; iEff++){
  //for(int iEff = 1; iEff < 2; iEff++){//MuonID efficiency (Zongchang)
  //for(int iEff = 2; iEff < 3; iEff++){//Muon Quality efficiency (Xianyou)
  //for(int iEff = 3; iEff < 4; iEff++){//L1/L2 trigger efficiency (Francesco)
   for(int iEff = 4; iEff < 5; iEff++){//L3 trigger efficiency (Luigi)
  // for(int iEff = 5; iEff < 6; iEff++){//Tracking trigger eff. of TkMu0 (Herbert)
  // for(int iEff = 6; iEff < 7; iEff++){//Muon trigger eff. of TkMu0 (Ilse)
  //for(int iEff = 7; iEff < 8; iEff++){//Tracking trigger eff. of TkMu0 + L3 (Mu0) (Herbert)
  // for(int iEff = 8; iEff < 9; iEff++){//Muon trigger eff. of TkMu0 + L1/L2 (Ilse)

    for(int iEffSample = 0; iEffSample < kNbEffSample; iEffSample++) //DATA, MC, MCTRUTH
      LoadEfficiencies(iEff, iEffSample);
    
    // Int_t maxEta = kNbEta;
    // for(int iEta = 0; iEta < maxEta; iEta++)
    //   PlotEff_DataVsMC_PT(iEff, iEta);

    // for(int iPT = 0; iPT < kNbPT; iPT++)
    //   PlotEff_DataVsMC_Eta(iEff, iPT);
 
    Int_t theEffSample = MCTRUTH;

    // printf("nBinsToMerge = %d\n", sizeof(etaIndicesToMerge) / sizeof(Int_t));
    ConcatenateGraphs(iEff, theEffSample, sizeof(etaIndicesToMerge) / sizeof(Int_t), pTMin);
    FitConcatenatedGraph(sizeof(etaIndicesToMerge) / sizeof(Int_t), iEff, theEffSample);

    // // fitPT = new TF1("fitPT", erfunction, 0., 100., 3);
    // // fitPT->SetParName(0, "xPos-Offset");
    // // fitPT->SetParName(1, "xPos-Plateau");
    // // fitPT->SetParName(2, "yPos-Plateau");
    // fitPT = new TF1("fitPT", erfunction2, 0., 100., 3);
    // fitPT->SetParName(0, "norm");
    // fitPT->SetParName(1, "xPos-Offset");
    // fitPT->SetParName(2, "rising-slope");
    // // fitPT_PlateauCorrected = new TF1("fitPT_PlateauCorrected", erfunction_corr, 0., 100., 3);
    // // fitPT_PlateauCorrected->SetParName(0, "xPos-Offset");
    // // fitPT_PlateauCorrected->SetParName(1, "xPos-Plateau");
    // // fitPT_PlateauCorrected->SetParName(2, "yPos-Plateau");

    // // for(int iEffSample = 1; iEffSample < 2; iEffSample++){
    // //   for(int iEta = 0; iEta < maxEta; iEta++){
    // // 	ParametrizePT(iEff, iEffSample, iEta, saveFit, plateauCorrection);
    // //   }
    // // }
  }
}

//==============================================================
void PlotEff_DataVsMC_PT(Int_t iEff, Int_t iEta){

  printf("will be plotting %s efficiency versus pT for eta slice %d\n", effName[iEff], iEta);
  Char_t name[100];
  sprintf(name, "c1_%s_pT_etaBin%d", effName[iEff], iEta);
  TCanvas *c1 = new TCanvas(name, name);
  Double_t minY = 0, maxY = 1.2;
  if(iEff == 6){minY = 0.6, maxY = 1.1;}
  TH1F *hFrame1 = gPad->DrawFrame(0., minY, 31., maxY);
  hFrame1->SetXTitle("p_{T}(#mu) GeV");
  
  gEff_pT[iEff][MC][iEta]->Draw("p same");   gEff_pT[iEff][MC][iEta]->SetLineColor(2); gEff_pT[iEff][MC][iEta]->SetMarkerColor(2);
  if(iEff != TkMu_Track)
    gEff_pT[iEff][MCTRUTH][iEta]->Draw("p same"); gEff_pT[iEff][MCTRUTH][iEta]->SetLineColor(4); gEff_pT[iEff][MCTRUTH][iEta]->SetMarkerColor(4);
  gEff_pT[iEff][DATA][iEta]->Draw("p same"); gEff_pT[iEff][DATA][iEta]->SetLineColor(1); gEff_pT[iEff][DATA][iEta]->SetMarkerColor(1);

  TLatex *tex1;
  if(iEff == 6)
    tex1 = new TLatex(2., 0.95*maxY, effNameLong[iEff]);
  else
    tex1 = new TLatex(2., 0.9*maxY, effNameLong[iEff]);
  tex1->SetTextSize(0.045); tex1->Draw();
  
  TLine *line1 = new TLine(0., 1., 31., 1.);
  line1->SetLineStyle(3); line1->Draw();

  if(iEff != 5){
    if(iEta == 0)
      sprintf(name, "|#eta(#mu)| < %1.1f", binsEta[iEta]);
    else
      sprintf(name, "%1.2f < |#eta(#mu)| < %1.2f", binsEta[iEta-1], binsEta[iEta]);
  }
  else{
    if(iEta == 0)
      sprintf(name, "|#eta(#mu)| < 1.2");
    else if(iEta == 1)
      sprintf(name, "1.2 < |#eta(#mu)| < 1.6");
    else if(iEta == 2)
      sprintf(name, "1.6 < |#eta(#mu)| < 2.4");
  }
  TLegend *leg = new TLegend(0.7456897,0.1610169,0.9454023,0.3622881, name);
  leg->AddEntry(gEff_pT[iEff][DATA][iEta], "Data T&P", "pl");
  leg->AddEntry(gEff_pT[iEff][MC][iEta], "MC T&P", "pl");
  if(iEff != TkMu_Track)
    leg->AddEntry(gEff_pT[iEff][MCTRUTH][iEta], "MC Truth", "pl");
  leg->SetFillColor(0); leg->SetTextSize(0.04); 
  leg->SetBorderSize(0); leg->Draw();

  sprintf(name, "Figures/%s_pT_DataVsMC_etaBin%d.pdf", effName[iEff], iEta);  c1->Print(name);
}

//==============================================================
void PlotEff_DataVsMC_Eta(Int_t iEff, Int_t iPT){

  Char_t name[100];
  sprintf(name, "c2_%s_eta_pTBin%d", effName[iEff], iPT);
  TCanvas *c1 = new TCanvas(name, name);
  Double_t minY = 0, maxY = 1.2;
  if(iEff == 6){minY = 0.6, maxY = 1.1;}
  TH1F *hFrame1 = gPad->DrawFrame(0., minY, 2.5, maxY);
  hFrame1->SetXTitle("|#eta(#mu)|");

  gEff_eta[iEff][MC][iPT]->Draw("p same");   gEff_eta[iEff][MC][iPT]->SetLineColor(2); gEff_eta[iEff][MC][iPT]->SetMarkerColor(2);
  if(iEff != TkMu_Track)
    gEff_eta[iEff][MCTRUTH][iPT]->Draw("p same"); gEff_eta[iEff][MCTRUTH][iPT]->SetLineColor(4); gEff_eta[iEff][MCTRUTH][iPT]->SetMarkerColor(4);
  gEff_eta[iEff][DATA][iPT]->Draw("p same"); gEff_eta[iEff][DATA][iPT]->SetLineColor(1); gEff_eta[iEff][DATA][iPT]->SetMarkerColor(1);

  TLatex *tex1;
  if(iEff == 6) 
    tex1 = new TLatex(0.2, 0.95*maxY, effNameLong[iEff]);
  else
    tex1 = new TLatex(0.2, 0.9*maxY, effNameLong[iEff]);
  tex1->SetTextSize(0.045); tex1->Draw();
  
  TLine *line1 = new TLine(0., 1., 2.5, 1.);
  line1->SetLineStyle(3); line1->Draw();

  sprintf(name, "%1.1f < p_{T}(#mu) < %1.1f GeV", binsPT[iPT], binsPT[iPT+1]);
  TLegend *leg = new TLegend(0.30,0.1652542,0.50,0.3665254, name);
  leg->AddEntry(gEff_eta[iEff][DATA][iPT], "Data T&P", "pl");
  leg->AddEntry(gEff_eta[iEff][MC][iPT], "MC T&P", "pl");
  if(iEff != TkMu_Track)
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
    // gEff_eta[iEff][iEffSample][iPT]->SetMarkerColor(colourEta[iPT]);
    // gEff_eta[iEff][iEffSample][iPT]->SetLineColor(colourEta[iPT]);
  }
  Int_t maxEta = kNbEta;
  for(int iEta = 0; iEta < maxEta; iEta++){
    sprintf(name, "gEff_%s_PT_AETA%d", effSampleName[iEffSample], iEta);
    gEff_pT[iEff][iEffSample][iEta] = (TGraphAsymmErrors *) gDirectory->Get(name);
    printf("%s, %s pT differential for eta bin %d is: %p\n", effName[iEff], effSampleName[iEffSample], iEta, gEff_pT[iEff][iEffSample][iEta]);
    gEff_pT[iEff][iEffSample][iEta]->SetMarkerStyle(marker[iEffSample]);
    // gEff_pT[iEff][iEffSample][iEta]->SetMarkerColor(colourPT[iEta]);
    // gEff_pT[iEff][iEffSample][iEta]->SetLineColor(colourPT[iEta]);
  }
}


//==========================================
void ParametrizePT(Int_t iEff, Int_t iEffSample, Int_t iEta, Bool_t saveFit, Bool_t plateauCorrection){

  Char_t name[100];
  sprintf(name, "c3_%s_%s_pT_etaBin%d", effName[iEff], effSampleName[iEffSample], iEta);
  TCanvas *c1 = new TCanvas(name, name);
  Double_t minY = 0, maxY = 1.2;
  if(iEff == 6){minY = 0.6, maxY = 1.1;}
  TH1F *hFrame1 = gPad->DrawFrame(0., minY, 30., maxY);
  hFrame1->SetXTitle("p_{T} [GeV]");

  gEff_pT[iEff][iEffSample][iEta]->Draw("p same"); gEff_pT[iEff][iEffSample][iEta]->SetLineColor(1); gEff_pT[iEff][iEffSample][iEta]->SetMarkerColor(1);
  fitPT->SetParameters(1., 3., 0.8);
  fitPT->SetParLimits(0., 0.9, 1.1);
  fitPT->SetParLimits(1, 0., 5.);
  fitPT->SetParLimits(2, 0., 2.);
  //following is useful for the 14 bins in Eta
  // if(iEta <= 5)
  //   gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT", "", "", 3., 30.);
  // else if(iEta <= 8)
  //   gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT", "", "", 2., 30.);
  // else if(iEta <= 10)
  //   gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT", "", "", 1.5, 30.);
  // else 
  //   gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT", "", "", 0.5, 30.);

  //following is useful for the 8 bins in eta:
  if(iEta <= 4)
    gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT", "", "", 3., 30.);
  else if(iEta <= 5)
    gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT", "", "", 1.8, 30.);
  else 
    gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT", "", "", 0.5, 30.);

  TF1 *fitFunc;
  
  sprintf(name, "fit%s_%s_pt_etaBin%d", effName[iEff], effSampleName[iEffSample], iEta);
  fitFunc = new TF1(name, erfunction2, 0., 100., 3);
  fitPT = gEff_pT[iEff][iEffSample][iEta]->GetFunction("fitPT");
  for(int iPar = 0; iPar < 3; iPar++)
    fitFunc->FixParameter(iPar, fitPT->GetParameter(iPar));

  fitFunc->SetName(name);
  fitFunc->SetLineWidth(1);
  fitFunc->SetLineColor(2);
  fitFunc->Draw("same");
  
  if(saveFit){
    TFile *fOut = new TFile(effFileNames[iEff], "UPDATE");
    fitFunc->Write();
    fOut->Close();
  }

  TLatex *tex1;
  if(iEff == 6) 
    tex1 = new TLatex(2., 0.95*maxY, effNameLong[iEff]);
  else
    tex1 = new TLatex(2., 0.9*maxY, effNameLong[iEff]);
  tex1->SetTextSize(0.045); tex1->Draw();

  if(iEff != 5){
    if(iEta == 0)
      sprintf(name, "|#eta(#mu)| < %1.1f", binsEta[iEta]);
    else
      sprintf(name, "%1.2f < |#eta(#mu)| < %1.2f", binsEta[iEta-1], binsEta[iEta]);
  }
  else{
    if(iEta == 0)
      sprintf(name, "|#eta(#mu)| < 1.2");
    else if(iEta == 1)
      sprintf(name, "1.2 < |#eta(#mu)| < 1.6");
    else if(iEta == 2)
      sprintf(name, "1.6 < |#eta(#mu)| < 2.4");
  }

  TLine *line1 = new TLine(0., 1., 30., 1.);
  line1->SetLineStyle(3); line1->Draw();

  TLegend *leg = new TLegend(0.7456897,0.2610169,0.9454023,0.3622881, name);
  leg->AddEntry(gEff_pT[iEff][iEffSample][iEta], effSampleName[iEffSample], "pl");
  leg->SetFillColor(0); leg->SetTextSize(0.04); 
  leg->SetBorderSize(0); leg->Draw();

  sprintf(name, "Figures/%s_pT_fitted_etaBin%d.pdf", effName[iEff], iEta);
  c1->Print(name);
}


//==========================================
void ConcatenateGraphs(Int_t iEff, Int_t iEffSample, Int_t nBinsToMerge, Double_t pTMin){

  Double_t pT[300], eff[300], errPTLow[300], errPTHigh[300], errEffLow[300], errEffHigh[300];
  Double_t pTOfRunningGraph, effOfRunningGraph;
  Int_t index = 0;
  Int_t maxPT;
  for(int iBin = 0; iBin < nBinsToMerge; iBin++){
    printf("taking pT values of eta bin %d\n", etaIndicesToMerge[iBin]);
    maxPT = TMath::Min(kNbPT, gEff_pT[iEff][iEffSample][etaIndicesToMerge[iBin]]->GetN());
    for(int iPT = 0; iPT < maxPT; iPT++){

      gEff_pT[iEff][iEffSample][etaIndicesToMerge[iBin]]->GetPoint(iPT, pTOfRunningGraph, effOfRunningGraph);      
      if(pTOfRunningGraph < pTMin) //discard points below a value of "pTMin"
	continue;

      printf("index has value %d (pT %f)\n", index, pTOfRunningGraph);
      pT[index] = pTOfRunningGraph + kNbPTMax*iBin;
      eff[index] = effOfRunningGraph;
      errPTLow[index] = gEff_pT[iEff][iEffSample][etaIndicesToMerge[iBin]]->GetErrorXlow(iPT);
      errPTHigh[index] = gEff_pT[iEff][iEffSample][etaIndicesToMerge[iBin]]->GetErrorXhigh(iPT);
      errEffLow[index] = gEff_pT[iEff][iEffSample][etaIndicesToMerge[iBin]]->GetErrorYlow(iPT);
      errEffHigh[index] = gEff_pT[iEff][iEffSample][etaIndicesToMerge[iBin]]->GetErrorYhigh(iPT);
      index++;
    }
  }

  //make a large concatenated TGraph:
  gEffPT_allEta = new TGraphAsymmErrors(index, pT, eff,  errPTLow, errPTHigh, errEffLow, errEffHigh);
  gEffPT_allEta->SetMarkerStyle(marker[iEffSample]);
  gEffPT_allEta->SetMarkerSize(0.7);
}

//==========================================
Double_t erfunction2(Double_t *x, Double_t *par){

  Double_t value;
  Double_t norm = par[0];
  value = norm/2. *(1. + TMath::Erf((x[0] - par[1])/(sqrt(2.)*par[2])));
  // if(value > 1.) 
  //   return 1.;
  // if(value < 0.)
  //   return 0.;

  return value;
}

//==========================================
Double_t erfunction(Double_t *x, Double_t *par){

  Double_t value;
  value = par[2] + TMath::Erf(par[0] + x[0]/ par[1]);
  if(value > 1.) 
    return 1.;
  // if(value < 0.)
  //   return 0.;

  return value;
}

//==========================================
Double_t erfunction_corr(Double_t *x, Double_t *par){

  Double_t value;
  value = par[2] + TMath::Erf(par[0] + x[0]/par[1]);

  if(value < 1.0 && x[0] > (par[1]+par[0]) && x[0] < 20.){
    Double_t k = (value-1.) / ((par[1]+par[0]) - 20.);
    Double_t d = 1.-20.*k;
    value = k*x[0] + d;
  }
  else if(x[0] > 20.)
    value = 1.;

  if(value < 0.)
    value = 0.;
  // if(value > 1.)
  //   value = 1.;

  return value;
}

//==========================================
Double_t fitConcatenatedGraph(Double_t *x, Double_t *par){

  Int_t iBin = ((Int_t) x[0]) / kNbPTMax;
  Double_t pT = x[0] - iBin * kNbPTMax;
  // printf("original pT %1.3f, calculated pT %1.3f\n", x[0], pT);

  Double_t norm = par[0];
  Double_t value;
  if(pT < pTMin)
    value = 0.;
  else
    value= norm/2. *(1. + TMath::Erf((pT - par[1])/(sqrt(2.)*par[2])));

  return value;
}

//==========================================
void FitConcatenatedGraph(Int_t nBins, Int_t iEff, Int_t iEffSample){

  Char_t name[100];
  sprintf(name, "cFit_%s_%s", effName[iEff], effNameLong[iEffSample]);
  TCanvas *c = new TCanvas(name, "fit of concatenated graph", 1200, 500);
  Double_t minX = 0.7;
  TH1F *hFrame = gPad->DrawFrame(0., minX, kNbPTMax*nBins+1., 1.1);
  hFrame->SetYTitle("Efficiency");
  hFrame->SetXTitle("iEta * pT [GeV/c]");

  gEffPT_allEta->Draw("psame");

  fitPT = new TF1("fitPT", fitConcatenatedGraph, 0., kNbPTMax*nBins, 3);
  fitPT->SetParName(0, "norm");
  fitPT->SetParName(1, "xPos-Offset");
  fitPT->SetParName(2, "rising-slope");
  //MIND: if a parameter is fixed, NDF of "calcChi2" needs to be adjusted!!!
  fitPT->SetParameters(1., 3., 0.8);
  // fitPT->SetParLimits(0., 0.9, 1.1);
  // fitPT->SetParLimits(1, 0., 5.);
  // fitPT->SetParLimits(2, 0., 2.);
  fitPT->FixParameter(0, 1.);
  gEffPT_allEta->Fit("fitPT", "0", "", pTMin, kNbPTMax*nBins);

  TF1 *fitFunc;  
  sprintf(name, "fit%s_%s_pT", effName[iEff], effSampleName[iEffSample]);
  fitFunc = new TF1(name, fitConcatenatedGraph, 0., kNbPTMax*nBins, 3);
  fitPT = gEffPT_allEta->GetFunction("fitPT");
  Double_t globalChi2 = fitPT->GetChisquare();
  Int_t ndf = fitPT->GetNDF();
  for(int iPar = 0; iPar < 3; iPar++)
    fitFunc->FixParameter(iPar, fitPT->GetParameter(iPar));

  fitFunc->SetName(name);
  fitFunc->SetLineWidth(1);
  fitFunc->SetLineColor(2);
  fitFunc->Draw("same");
 
  TLatex *tex1 = new TLatex(0.05*kNbPTMax*nBins, 1.05, effNameLong[iEff]);
  tex1->SetTextSize(0.045); tex1->Draw();
  sprintf(name, "#chi^{2} / ndf = %1.2f / %d = %1.2f\n", globalChi2, ndf, globalChi2/ndf);
  tex1->DrawLatex(0.7*kNbPTMax*nBins, 1.02, name);

  for(int iBin = 0; iBin < nBins; iBin++){
    if(etaIndicesToMerge[iBin] == 0)
      sprintf(name, "|eta| < %1.2f", binsEta[etaIndicesToMerge[iBin]]);
    else
      sprintf(name, "%1.2f < |eta| < %1.2f", binsEta[etaIndicesToMerge[iBin]-1], binsEta[etaIndicesToMerge[iBin]]);
    tex1->DrawLatex(iBin*kNbPTMax+5., minX*1.05, name);

    sprintf(name, "#chi^{2}/ndf = %1.3f", calcChi2(gEff_pT[iEff][iEffSample][etaIndicesToMerge[iBin]], fitFunc, 2));
    tex1->DrawLatex(iBin*kNbPTMax+5., minX*1.01, name);
  }

  TLine *line1 = new TLine(0., 1., kNbPTMax*nBins+1, 1.);
  line1->SetLineStyle(3); line1->Draw();
  for(int iBin = 0; iBin < nBins; iBin++)
    line1->DrawLine(0., iBin*kNbPTMax, 1.1,  iBin*kNbPTMax);  

  // sprintf(name, "Figures/%s_
  // c->Print(name);
}

//==========================================
Double_t calcChi2(TGraphAsymmErrors *graph, TF1 *func, Int_t fitPar){

  Double_t pTOfRunningGraph, effOfRunningGraph;
  Double_t errEffLow, errEffHigh;
  Double_t chi2 = 0.;
  Int_t nBins = 0;
  for(int iPT = 0; iPT < graph->GetN(); iPT++){
    graph->GetPoint(iPT, pTOfRunningGraph, effOfRunningGraph);      

    if(pTOfRunningGraph < pTMin) //discard points below a value of "pTMin"
      continue;
    errEffLow = graph->GetErrorYlow(iPT);
    errEffHigh = graph->GetErrorYhigh(iPT);

    nBins++;
    printf("measured eff %1.3f, fitted eff %1.3f, errEffLow %1.3f, errEffHigh %1.3f\n",
	   effOfRunningGraph, func->Eval(pTOfRunningGraph), errEffLow, errEffHigh);
    if(errEffHigh < 1e-6)
      chi2 += pow(effOfRunningGraph - func->Eval(pTOfRunningGraph), 2)/(errEffLow*errEffLow);
    else if(errEffLow < 1e-6)
      chi2 += pow(effOfRunningGraph - func->Eval(pTOfRunningGraph), 2)/(errEffHigh*errEffHigh);
    else
      chi2 += pow(effOfRunningGraph - func->Eval(pTOfRunningGraph), 2)/(errEffLow*errEffHigh);
  }
  chi2 /= (Double_t) (nBins - fitPar);
  return chi2;
}
