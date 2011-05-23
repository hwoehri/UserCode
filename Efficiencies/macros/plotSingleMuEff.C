#include "../interface/rootIncludes.inc"
#include "TGraphAsymmErrors.h"
#include "TMath.h"

Int_t const kNbEff = 7;
Char_t *effFileNames[kNbEff] = {"/Users/hwoehri/CMS/Work/TnP/LucaPerrozzi/17March2011/MuonTrackingEff_17March2011.root",
				"/Users/hwoehri/CMS/Work/TnP/Zongchang/29March2011/MuonIDEff_29March2011_fitted.root",
				"/Users/hwoehri/CMS/Work/TnP/Xianyou/25March2011/MuonQualEff_25March2011_fitted.root",
				//"/Users/hwoehri/CMS/Work/TnP/Francesco/2April2011/L1L2_DMu0_TriggerEfficiencies_RUNB_distM2gt120_7April2011_fitted.root",
  				"/Users/hwoehri/CMS/Work/TnP/Francesco/2April2011/L1L2_DMu0_TriggerEfficiencies_RUNA_distM2gt120_7April2011_fitted.root",
//				"/Users/hwoehri/CMS/Work/TnP/Luigi/17April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010A_distM2gt120_17April2011.root",
//				"/Users/hwoehri/CMS/Work/TnP/Luigi/17April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM2gt120_17April2011.root",
				"/Users/hwoehri/CMS/Work/TnP/Luigi/17April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010A_distM2gt120_17April2011_fitted.root",
//				"/Users/hwoehri/CMS/Work/TnP/Luigi/17April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM2gt120_17April2011_fitted.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Herbert/28March2011/HLT_Track_Mu0TkMu0_Run1.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Herbert/29March2011/HLT_Track_Mu0TkMu0_Run_Tightv3.root",
				"/Users/hwoehri/CMS/Work/TnP/Herbert/13May2011/HLT_Track_Mu0TkMu0_Run1_13May2011.root", //
// 				"/Users/hwoehri/CMS/Work/TnP/Ilse/29March2011/HLTMuonTrack_Mu0_TkMu0_TM_runA_29March2011.root"};
				"/Users/hwoehri/CMS/Work/TnP/Ilse/29March2011/HLTMuonTrack_Mu0_TkMu0_TM_runB1_29March2011.root"};
//old version of scripts:
//  				"/Users/hwoehri/CMS/Work/TnP/Luigi/25March2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_25March2011_fitted.root",
//				"/Users/hwoehri/CMS/Work/TnP/Luigi/5April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010A_distM1gt150_7April2011_fitted.root",
//				"/Users/hwoehri/CMS/Work/TnP/Luigi/5April2011/L3_DoubleMu0_TriggerEfficiencies_Run2010B_distM1gt150_7April2011_fitted.root",
//  				"/Users/hwoehri/CMS/Work/TnP/Francesco/23March2011/L1L2_DMu0_TriggerEfficiencies_23March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Zongchang/19March2011/MuonIDEff_19March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Francesco/18March2011/L1L2_DoubleMu0_TriggerEfficiencies_18March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Luigi/25March2011/L3_DoubleMu0_TriggerEfficiencies_RunABmixed_25March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Luigi/19March2011/L3_DoubleMu0_TriggerEfficiencies_19March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Xianyou/19March2011/MuonQualEff_19March2011.root",
// 				"/Users/hwoehri/CMS/Work/TnP/Ilse/24March2011/HLTMuonTrack_Mu0_TkMu0_TM.root"};
enum {TrkEff, MuIDEff, MuQualEff, L1L2Eff, L3Eff};
Char_t *effName[kNbEff] = {"TrkEff", "MuIDEff", "MuQualEff", "L1L2Eff", "L3Eff", "Track-TkMu0", "Mu-TkMu0"};
Char_t *effNameLong[kNbEff] = {"offline tracking efficiency", "muon identification efficiency", 
				"efficiency of muon quality cuts", "L1-L2 trigger efficiency (Mu0)", 
			       "L3 trigger efficiency (Mu0)", "tracking efficiency of TkMu0",
			       "muon trigger efficency of TkMu0"};
Int_t const kNbEffSample = 3;
enum {DATA, MC, MCTRUTH};
Char_t *effSampleName[kNbEffSample] = {"DATA", "MC", "MCTRUTH"};
//
enum {CENTRAL, UPPER, LOWER};
Char_t *valName[3] = {"central", "lower", "upper"};
TH2D *hMuEff[kNbEff][kNbEffSample][3];
// Int_t const kNbEta = 3;//Track-TkMu0
//Int_t const kNbEta = 8;
Int_t const kNbEta = 14; //new L1/L2 and L3 efficiencies
// Int_t const kNbPT = 12; //"old" PT binning
// Int_t const kNbPT = 11; //"new" PT binning
Int_t const kNbPT = 15; //"new" L1/L2 and L3 PT binning
TGraphAsymmErrors *gEff_pT[kNbEff][kNbEffSample][kNbEta];
TGraphAsymmErrors *gEff_eta[kNbEff][kNbEffSample][kNbPT];
//
TF1 *fitPT, *fitPT_PlateauCorrected;

Int_t const marker[kNbEffSample] = {20, 25, 21};
Int_t const colourEta[kNbPT] = {1,2,3,4,kYellow-6,6,7,8,9,10,11};
//Int_t const colourPT[kNbEta] = {1,2,3,4,kYellow-6,6,7,8};
Int_t const colourPT[kNbEta] = {1,2,4};
//Float_t binsEta[kNbEta] = {1.2, 1.6,  2.4}; //Track-TkMu0
//Float_t binsEta[kNbEta] = {0.2, 0.3, 0.6, 0.8, 1.2, 1.6, 2.1, 2.4};
Float_t binsEta[kNbEta] = {0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4};//new L1/L2 and L3 binning
//Float_t binsPT[kNbPT+1] = {0.7, 1., 1.25, 1.5, 1.75, 2., 2.5, 3.5, 4., 6., 10., 20., 30.}; //old binning
//Float_t binsPT[kNbPT+1] = {0.8, 1.25, 1.5, 1.75, 2., 2.5, 3.3, 4., 6., 10., 20., 30.};//standard binning
Float_t binsPT[kNbPT+1] = {0.8, 1.25, 1.5, 1.75, 2., 2.25, 2.5, 3, 3.3,  4.,  5., 6., 7., 10., 20., 30.};//new L1/L2 and L3 binning

void LoadEfficiencies(Int_t iEff, Int_t iEffSample);
void PlotEff_DataVsMC_PT(Int_t iEff, Int_t iEta);
void PlotEff_DataVsMC_Eta(Int_t iEff, Int_t iPT);
Double_t erfunction(Double_t *x, Double_t *par);
Double_t erfunction2(Double_t *x, Double_t *par);
Double_t erfunction_corr(Double_t *x, Double_t *par);
void ParametrizePT(Int_t iEff, Int_t iEffSample, Int_t iEta, Bool_t saveFit, Bool_t plateauCorrection);
//==============================================================
void plotSingleMuEff(Bool_t saveFit = kFALSE,
		     Bool_t plateauCorrection = kFALSE){

  //for(int iEff = 0; iEff < kNbEff; iEff++){
  //for(int iEff = 1; iEff < 2; iEff++){//MuonID efficiency (Zongchang)
  //for(int iEff = 2; iEff < 3; iEff++){//Muon Quality efficiency (Xianyou)
  // for(int iEff = 3; iEff < 4; iEff++){//L1/L2 trigger efficiency (Francesco)
  for(int iEff = 4; iEff < 5; iEff++){//L3 trigger efficiency (Luigi)
  //for(int iEff = 5; iEff < 6; iEff++){//Tracking trigger eff. of TkMu0 (Herbert)
    //for(int iEff = 6; iEff < 7; iEff++){//Muon trigger eff. of TkMu0 (Ilse)
    for(int iEffSample = 0; iEffSample < kNbEffSample; iEffSample++) //DATA, MC, MCTRUTH

      LoadEfficiencies(iEff, iEffSample);
    
    Int_t maxEta = kNbEta;
    for(int iEta = 0; iEta < maxEta; iEta++)
      PlotEff_DataVsMC_PT(iEff, iEta);

    for(int iPT = 0; iPT < kNbPT; iPT++)
      PlotEff_DataVsMC_Eta(iEff, iPT);
 

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

    // for(int iEffSample = 0; iEffSample < 1; iEffSample++){
    //   for(int iEta = 0; iEta < maxEta; iEta++){
    // 	ParametrizePT(iEff, iEffSample, iEta, saveFit, plateauCorrection);
    //   }
    // }
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
  
  // gEff_pT[iEff][MC][iEta]->Draw("p same");   gEff_pT[iEff][MC][iEta]->SetLineColor(2); gEff_pT[iEff][MC][iEta]->SetMarkerColor(2);
  // gEff_pT[iEff][MCTRUTH][iEta]->Draw("p same"); gEff_pT[iEff][MCTRUTH][iEta]->SetLineColor(4); gEff_pT[iEff][MCTRUTH][iEta]->SetMarkerColor(4);
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
      sprintf(name, "%1.1f < |#eta(#mu)| < %1.1f", binsEta[iEta-1], binsEta[iEta]);
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
  // leg->AddEntry(gEff_pT[iEff][MC][iEta], "MC T&P", "pl");
  // leg->AddEntry(gEff_pT[iEff][MCTRUTH][iEta], "MC Truth", "pl");
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

  // gEff_eta[iEff][MC][iPT]->Draw("p same");   gEff_eta[iEff][MC][iPT]->SetLineColor(2); gEff_eta[iEff][MC][iPT]->SetMarkerColor(2);
  // gEff_eta[iEff][MCTRUTH][iPT]->Draw("p same"); gEff_eta[iEff][MCTRUTH][iPT]->SetLineColor(4); gEff_eta[iEff][MCTRUTH][iPT]->SetMarkerColor(4);
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
  // leg->AddEntry(gEff_eta[iEff][MC][iPT], "MC T&P", "pl");
  // leg->AddEntry(gEff_eta[iEff][MCTRUTH][iPT], "MC Truth", "pl");
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
  Int_t maxEta = kNbEta;
  for(int iEta = 0; iEta < maxEta; iEta++){
    sprintf(name, "gEff_%s_PT_AETA%d", effSampleName[iEffSample], iEta);
    gEff_pT[iEff][iEffSample][iEta] = (TGraphAsymmErrors *) gDirectory->Get(name);
    printf("%s, %s pT differential for eta bin %d is: %p\n", effName[iEff], effSampleName[iEffSample], iEta, gEff_pT[iEff][iEffSample][iEta]);
    gEff_pT[iEff][iEffSample][iEta]->SetMarkerStyle(marker[iEffSample]);
    gEff_pT[iEff][iEffSample][iEta]->SetMarkerColor(colourPT[iEta]);
    gEff_pT[iEff][iEffSample][iEta]->SetLineColor(colourPT[iEta]);
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
  // if((iEta == 1 || iEta == 3) || !plateauCorrection){
  // if(!plateauCorrection){
    fitPT->SetParameters(1., 3., 0.8);
    fitPT->SetParLimits(0., 0.9, 1.1);
    fitPT->SetParLimits(1, 0., 5.);
    fitPT->SetParLimits(2, 0., 2.);
    if(iEta <= 5)
      gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT", "", "", 3., 30.);
    else if(iEta <= 8)
      gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT", "", "", 2., 30.);
    else if(iEta <= 10)
      gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT", "", "", 1.5, 30.);
    else 
      gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT", "", "", 0.5, 30.);
  // }
  // else{
  //   fitPT_PlateauCorrected->SetParameter(0, -0.6);
  //   fitPT_PlateauCorrected->SetParameter(1, 7.);
  //   fitPT_PlateauCorrected->SetParameter(2, -0.01);
  //   if(iEta > 4)
  //     gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT_PlateauCorrected", "", "", 0., 20.);
  //   else if(iEta == 3) //L3 trigger efficiency
  //     gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT_PlateauCorrected", "", "", 2., 10.);
  //   else
  //   gEff_pT[iEff][iEffSample][iEta]->Fit("fitPT_PlateauCorrected", "", "", 3., 20.);
  // }

  TF1 *fitFunc;
  
  // if(!plateauCorrection){
    sprintf(name, "fit%s_%s_pt_etaBin%d", effName[iEff], effSampleName[iEffSample], iEta);
    fitFunc = new TF1(name, erfunction2, 0., 100., 3);
    fitPT = gEff_pT[iEff][iEffSample][iEta]->GetFunction("fitPT");
    for(int iPar = 0; iPar < 3; iPar++)
      fitFunc->FixParameter(iPar, fitPT->GetParameter(iPar));
  // }
  // else{
  //   sprintf(name, "fit%s_%s_pt_etaBin%d", effName[iEff], effSampleName[iEffSample], iEta);
  //   fitFunc = new TF1(name, erfunction_corr, 0., 100., 3);
  //   fitPT_PlateauCorrected = gEff_pT[iEff][iEffSample][iEta]->GetFunction("fitPT_PlateauCorrected");
  //   for(int iPar = 0; iPar < 3; iPar++){
  //     fitFunc->FixParameter(iPar, fitPT_PlateauCorrected->GetParameter(iPar));
  //     printf("setting parameter %d to %1.3f\n", iPar, fitPT_PlateauCorrected->GetParameter(iPar));
  //   }
  // }
  fitFunc->SetName(name);
  fitFunc->SetLineWidth(1);
  fitFunc->SetLineColor(2);
  fitFunc->Draw("same");
  //fitFunc->SetRange(0., 100.);
  
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
      sprintf(name, "%1.1f < |#eta(#mu)| < %1.1f", binsEta[iEta-1], binsEta[iEta]);
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
