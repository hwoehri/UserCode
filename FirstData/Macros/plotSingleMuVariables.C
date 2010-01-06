#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include "commonVar.inc"

const Int_t kNbSandBG = 3;
const Char_t *sBGName[kNbSandBG] = {"signal", "BG", "data"};
const Char_t *sMuonType[kNbMuCath] = {"global", "tracker", "calo"};
Int_t marker[kNbSandBG] = {20, 24};
Int_t colour[kNbMuCath] = {1, 2, 4};
Int_t colourDimu[kNbCath] = {1, 2, 4, 3, 6, 7};

//single muon histos for two dimuon mass windows:
//const Int_t kNbDimuonSet = 2; //[0]..3.0 < M < 3.2; [1]... M > 2 GeV
TH1F *hMuon_eta[kNbSandBG][kNbDimuonSet][kNbCharge][kNbCath];
TH1F *hMuon_pT[kNbSandBG][kNbDimuonSet][kNbCharge][kNbCath];
TH2F *hMuon_pT_eta[kNbSandBG][kNbDimuonSet][kNbCharge][kNbMuCath];
//quality histograms
TH1F *hD0[kNbSandBG][kNbDimuonSet][kNbCharge][kNbMuCath];
TH1F *hDz[kNbSandBG][kNbDimuonSet][kNbCharge][kNbMuCath];
TH1F *hChi2GlobalFit[kNbSandBG][kNbDimuonSet][kNbCharge];//global only
TH1F *hChi2TrackerFit[kNbSandBG][kNbDimuonSet][kNbCharge][kNbMuCath];
TH1F *hTrackerMuonArbitrated[kNbSandBG][kNbDimuonSet][kNbCharge];//tracker only
TH1F *hTM2DCompatibilityTight[kNbSandBG][kNbDimuonSet][kNbCharge];//tracker only
TH1F *hTMLastStationOptimizedLowPtLoose[kNbSandBG][kNbDimuonSet][kNbCharge];//tracker only
TH1F *hCaloCompatibility[kNbSandBG][kNbDimuonSet][kNbCharge][kNbMuCath];
TH1F *hDimuVtxProb[kNbSandBG][kNbDimuonSet][kNbCharge][kNbCath];
TH1F *hNHitsSilicon[kNbSandBG][kNbDimuonSet][kNbCharge][kNbMuCath];

void PlotHistos(Int_t iSet, Int_t iCharge);
void LoadHistos(Char_t *fileNameIn, Int_t iS);
//================================
void plotSingleMuVariables(Char_t *fileNameIn="histos_jPsi_2360GeV_STARTUP.root",//MC signal
			   //			   Char_t *fileNameIn2="histos_MB_2360_STARTUP.root", //MC BG
			   Char_t *fileNameIn2 = "histos_ppMuXLoose_2360GeV_STARTUP.root",
                           Char_t *fileNameIn3="histos_data09_2360GeV.root",//data
			   Int_t iSet = 1){ //[0] 3.0<M<3.2, [1] M > 2 GeV, [2]... M < 4 GeV

  LoadHistos(fileNameIn, 0);
  printf("first set of histos loaded\n");
  LoadHistos(fileNameIn2, 1);
  printf("second set of histos loaded\n");
  LoadHistos(fileNameIn3, 2);
  printf("third set of histos loaded\n");

  Int_t iCharge = 0;//OS only
  PlotHistos(iSet, iCharge);
}
//=================================
void PlotHistos(Int_t iSet, Int_t iCharge){

  gStyle->SetOptStat(0);
  Char_t name[100], title[100];

  //=========================================
  //global muons:
  //=========================================
  TCanvas *c1 = new TCanvas("c1", "global chi2: signal and BG");
  gPad->SetLogy();
  TH1F *hFrame1 = gPad->DrawFrame(0., 1e-5, 10., 1.);
  hFrame1->SetXTitle(hChi2GlobalFit[0][iSet][iCharge]->GetXaxis()->GetTitle());
  hChi2GlobalFit[0][iSet][iCharge]->SetLineColor(4);
  hChi2GlobalFit[0][iSet][iCharge]->DrawNormalized("same");
  hChi2GlobalFit[1][iSet][iCharge]->SetLineColor(2);
  hChi2GlobalFit[1][iSet][iCharge]->DrawNormalized("same");
  hChi2GlobalFit[2][iSet][iCharge]->SetMarkerStyle(20);
  hChi2GlobalFit[2][iSet][iCharge]->DrawNormalized("histsame");
  TLatex *tex1 = new TLatex(1.0229, 0.434, "MC global muons from background muon pairs of M > 2 GeV");
  tex1->SetTextSize(0.04); tex1->Draw(); tex1->SetTextColor(2);
  TLatex *tex1a = new TLatex(4.809, 0.202,"MC global muons from J/#psi decays");
  tex1a->SetTextSize(0.04); tex1a->Draw(); tex1a->SetTextColor(4);
  TLatex *tex1b = new TLatex(9.0709, 0.0954,"Data");
  tex1b->SetTextSize(0.04); tex1b->Draw(); tex1b->SetTextColor(1);

  c1->Print("Figures/chi2GlobalFit_signalBGAndData.gif");

  //d0 (of all muons)
  TCanvas *c1a[kNbMuCath];
  TH1F *hFrame1a[kNbMuCath];
  TLatex *tex1aa[kNbMuCath], *tex1aaa[kNbMuCath], *tex1ba[kNbMuCath];
  for(int iCat = 0; iCat < kNbMuCath; iCat++){
    sprintf(name, "c1a_%d", iCat);
    c1a[iCat] = new TCanvas(name, "d0: signal and BG");
    gPad->SetLogy();
    hFrame1a[iCat] = gPad->DrawFrame(0., 1e-5, 5., 1.);
    hFrame1a[iCat]->SetXTitle(hD0[0][iSet][iCharge][iCat]->GetXaxis()->GetTitle());
    hD0[0][iSet][iCharge][iCat]->SetLineColor(4);
    hD0[0][iSet][iCharge][iCat]->DrawNormalized("same");
    hD0[1][iSet][iCharge][iCat]->SetLineColor(2);
    hD0[1][iSet][iCharge][iCat]->DrawNormalized("same");
    hD0[2][iSet][iCharge][iCat]->SetMarkerStyle(20);
    hD0[2][iSet][iCharge][iCat]->DrawNormalized("histsame");
    
    sprintf(name, "MC %s muons from background muon pairs of M > 2 GeV", sMuonType[iCat]);
    tex1aa[iCat] = new TLatex(0.46, 0.434, name);
    tex1aa[iCat]->SetTextSize(0.04); tex1aa[iCat]->Draw(); tex1aa[iCat]->SetTextColor(2);
    sprintf(name, "MC %s muons from J/#psi decays", sMuonType[iCat]);
    tex1aaa[iCat] = new TLatex(2.33, 0.202, name);
    tex1aaa[iCat]->SetTextSize(0.04); tex1aaa[iCat]->Draw(); tex1aaa[iCat]->SetTextColor(4);
    tex1ba[iCat] = new TLatex(4.52, 0.0954,"Data");
    tex1ba[iCat]->SetTextSize(0.04); tex1ba[iCat]->Draw(); tex1ba[iCat]->SetTextColor(1);
    
    sprintf(name, "d0_%s_signalBGAndData.gif", sMuonType[iCat]);
    c1a[iCat]->Print(name);
  }
 
  //dz (of all muons)
  TCanvas *c1b[kNbMuCath];
  TH1F *hFrame1b[kNbMuCath];
  TLatex *tex100[kNbMuCath], *tex100a[kNbMuCath], *tex100b[kNbMuCath];
  for(int iCat = 0; iCat < kNbMuCath; iCat++){
    sprintf(name, "c1b_%d", iCat);
   c1b[iCat] = new TCanvas(name, "dz: signal and BG");
   gPad->SetLogy();
   hFrame1b[iCat] = gPad->DrawFrame(0., 1e-4, 20., 1.);
   hFrame1b[iCat]->SetXTitle(hDz[0][iSet][iCharge][iCat]->GetXaxis()->GetTitle());
   hDz[0][iSet][iCharge][iCat]->SetLineColor(4);
   hDz[0][iSet][iCharge][iCat]->DrawNormalized("same");
   hDz[1][iSet][iCharge][iCat]->SetLineColor(2);
   hDz[1][iSet][iCharge][iCat]->DrawNormalized("same");
   hDz[2][iSet][iCharge][iCat]->SetMarkerStyle(20);
   hDz[2][iSet][iCharge][iCat]->DrawNormalized("histsame");

   sprintf(name, "MC %s muons from background muon pairs of M > 2 GeV", sMuonType[iCat]);
   tex100[iCat] = new TLatex(1.7, 0.434, name);
   tex100[iCat]->SetTextSize(0.04); tex100[iCat]->Draw(); tex100[iCat]->SetTextColor(2);
   sprintf(name, "MC %s muons from J/#psi decays", sMuonType[iCat]);
   tex100a[iCat] = new TLatex(9.3, 0.202, name);
   tex100a[iCat]->SetTextSize(0.04); tex100a[iCat]->Draw(); tex100a[iCat]->SetTextColor(4);
   tex100b[iCat] = new TLatex(18., 0.0954,"Data");
   tex100b[iCat]->SetTextSize(0.04); tex100b[iCat]->Draw(); tex100b[iCat]->SetTextColor(1);
   
   sprintf(name, "dz_%s_signalBGAndData.gif", sMuonType[iCat]);
   c1b[iCat]->Print(name);
  }

  //=========================================
  //tracker muons
  //=========================================
  TCanvas *c2[kNbMuCath];
  TH1F *hFrame2[kNbMuCath];
  TLatex *tex2[kNbMuCath], *tex2a[kNbMuCath], *tex2b[kNbMuCath];
  for(int iCat = 1; iCat < kNbMuCath; iCat++){
    sprintf(name, "c2_%d", iCat);
    sprintf(title, "tracker chi2 for mu %d", iCat);
    c2[iCat] = new TCanvas(name, title);
    hFrame2[iCat] = gPad->DrawFrame(0., 1e-5, 15., 1.);
    hFrame2[iCat]->SetXTitle(hChi2TrackerFit[0][iSet][iCharge][iCat]->GetXaxis()->GetTitle());
    gPad->SetLogy();
    hChi2TrackerFit[0][iSet][iCharge][iCat]->SetLineColor(4);//signal
    hChi2TrackerFit[0][iSet][iCharge][iCat]->DrawNormalized("same");
    hChi2TrackerFit[1][iSet][iCharge][iCat]->SetLineColor(2);//BG
    hChi2TrackerFit[1][iSet][iCharge][iCat]->DrawNormalized("same");
    hChi2TrackerFit[2][iSet][iCharge][iCat]->SetMarkerStyle(20);
    hChi2TrackerFit[2][iSet][iCharge][iCat]->DrawNormalized("histsame");

    sprintf(name, "MC %s muons from background muon pairs of M > 2 GeV", sMuonType[iCat]);
    tex2[iCat] = new TLatex(1.331643,0.4740476, name);
    tex2[iCat]->SetTextSize(0.04); tex2[iCat]->Draw(); tex2[iCat]->SetTextColor(2);
    sprintf(name, "MC %s muons from J/#psi decays", sMuonType[iCat]);
    tex2a[iCat] = new TLatex(6.935091,0.1963064, name);
    tex2a[iCat]->SetTextSize(0.04); tex2a[iCat]->Draw(); tex2a[iCat]->SetTextColor(4);
    tex2b[iCat] = new TLatex(13.3499, 0.0837,"Data");
    tex2b[iCat]->SetTextSize(0.04); tex2b[iCat]->Draw(); tex2b[iCat]->SetTextColor(1);

    sprintf(name, "Figures/hChi2TrackerFit_%s_signalBGAndData.gif", sMuonType[iCat]);
    c2[iCat]->Print(name);
  }

  //=========================================
  TCanvas *c14 = new TCanvas("c14", "hNHitsSilicon: signal and BG (tracker)");
  TH1F *hFrame14 = gPad->DrawFrame(0., 0., 32, 0.12);
  hFrame14->SetXTitle("#hits in silicon tracker");
  hNHitsSilicon[0][iSet][iCharge][1]->SetLineColor(4);
  hNHitsSilicon[0][iSet][iCharge][1]->DrawNormalized("same");//signal
  hNHitsSilicon[1][iSet][iCharge][1]->SetLineColor(2);//BG
  hNHitsSilicon[1][iSet][iCharge][1]->DrawNormalized("same");//BG
  hNHitsSilicon[2][iSet][iCharge][1]->SetMarkerStyle(20);//data
  hNHitsSilicon[2][iSet][iCharge][1]->SetLineColor(1);//data
  hNHitsSilicon[2][iSet][iCharge][1]->DrawNormalized("same");//data
  TLatex *tex14 = new TLatex(1., 0.11, "MC tracker muons from background muon pairs of M > 2 GeV");
  tex14->SetTextSize(0.04); tex14->Draw(); tex14->SetTextColor(2);
  TLatex *tex14a = new TLatex(1., 0.10,"MC tracker muons from J/#psi decays");
  tex14a->SetTextSize(0.04); tex14a->Draw(); tex14a->SetTextColor(4);
  TLatex *tex14b = new TLatex(1., 0.09,"Data");
  tex14b->SetTextSize(0.04); tex14b->Draw(); tex14b->SetTextColor(1);

  c14->Print("Figures/hNHitsSilicon_tracker_signalBGAndData.pdf");
  c14->Print("Figures/hNHitsSilicon_tracker_signalBGAndData.gif");

  //=========================================
  TCanvas *c15 = new TCanvas("c15", "hNHitsSilicon: signal and BG (calo)");
  TH1F *hFrame15 = gPad->DrawFrame(0., 0., 32, 0.2);
  hFrame15->SetXTitle("#hits in silicon tracker");
  hNHitsSilicon[0][iSet][iCharge][2]->SetLineColor(4);
  hNHitsSilicon[0][iSet][iCharge][2]->DrawNormalized("same");//signal
  hNHitsSilicon[1][iSet][iCharge][2]->SetLineColor(2);
  hNHitsSilicon[1][iSet][iCharge][2]->DrawNormalized("same");//BG
  hNHitsSilicon[2][iSet][iCharge][2]->SetMarkerStyle(20);
  hNHitsSilicon[2][iSet][iCharge][2]->DrawNormalized("same");//data

  TLatex *tex15 = new TLatex(1., 0.18, "MC calo muons from background muon pairs of M > 2 GeV");
  tex15->SetTextSize(0.04); tex15->Draw(); tex15->SetTextColor(2);
  TLatex *tex15a = new TLatex(1., 0.165,"MC calo muons from J/#psi decays");
  tex15a->SetTextSize(0.04); tex15a->Draw(); tex15a->SetTextColor(4);
  TLatex *tex15b = new TLatex(1., 0.15,"Data");
  tex15b->SetTextSize(0.04); tex15b->Draw(); tex15b->SetTextColor(1);

  c15->Print("Figures/hNHitsSilicon_calo_signalBGAndData.pdf");
  c15->Print("Figures/hNHitsSilicon_calo_signalBGAndData.gif");

  //=========================================
  //calo muons
  //=========================================
  TCanvas *c8 = new TCanvas("c8", "hCaloCompatibility: signal and BG");
  TH1F *hFrame8 = gPad->DrawFrame(0.4, 0., 1., 0.12);
  hFrame8->SetXTitle("CaloCompatibility");
  hCaloCompatibility[0][iSet][iCharge][2]->SetLineColor(4);
  hCaloCompatibility[0][iSet][iCharge][2]->DrawNormalized("same");//calo-signal
  hCaloCompatibility[1][iSet][iCharge][2]->SetLineColor(2);
  hCaloCompatibility[1][iSet][iCharge][2]->DrawNormalized("same");//calo-BG
  hCaloCompatibility[2][iSet][iCharge][2]->SetLineColor(1);
  hCaloCompatibility[2][iSet][iCharge][2]->SetMarkerStyle(20);
  hCaloCompatibility[2][iSet][iCharge][2]->DrawNormalized("psame");//calo-data

  TLatex *tex8 = new TLatex(0.42, 0.11, "MC calo muons from background muon pairs of M > 2 GeV");
  tex8->SetTextSize(0.04); tex8->Draw(); tex8->SetTextColor(2);
  TLatex *tex8a = new TLatex(0.42, 0.10,"MC calo muons from J/#psi decays");
  tex8a->SetTextSize(0.04); tex8a->Draw(); tex8a->SetTextColor(4);
  TLatex *tex8b = new TLatex(0.42, 0.09,"Data");
  tex8b->SetTextSize(0.04); tex8b->Draw(); tex8b->SetTextColor(1);

  c8->Print("Figures/hCaloCompatibility_calo_signalBGAndData.pdf");
  c8->Print("Figures/hCaloCompatibility_calo_signalBGAndData.gif");

  //=========================================
  //signal, background and data: tracker-tracker
  TCanvas *c17a = new TCanvas("c17a", "hDimuVtxProb: signal and BG (tracker muons)");
  TH1F *hFrame17a = gPad->DrawFrame(5e-5, 5e-5, 1., 1.3);
  hFrame17a->SetXTitle(hDimuVtxProb[1][iSet][iCharge][2]->GetXaxis()->GetTitle());
  gPad->SetLogx();  gPad->SetLogy();
  hDimuVtxProb[1][iSet][iCharge][2]->SetLineColor(2);
  hDimuVtxProb[1][iSet][iCharge][2]->DrawNormalized("same");//BG
  hDimuVtxProb[0][iSet][iCharge][2]->SetLineColor(4);//signal
  hDimuVtxProb[0][iSet][iCharge][2]->DrawNormalized("same");//signal
  hDimuVtxProb[2][iSet][iCharge][2]->SetLineColor(1);//data
  hDimuVtxProb[2][iSet][iCharge][2]->SetMarkerStyle(20);//
  hDimuVtxProb[2][iSet][iCharge][2]->DrawNormalized("same");//

  sprintf(name, "Figures/dimuVtxProb_%s_signalBGAndData.gif", oniaCatName[2]);
  c17a->Print(name);

}

//=================================
void LoadHistos(Char_t *fileNameIn, Int_t iS){

  Char_t name[100], title[300];
  TFile *fIn = new TFile(fileNameIn);

  for(int iSet = 0; iSet < kNbDimuonSet; iSet++){
//     for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
    for(int iCharge = 0; iCharge < 1; iCharge++){
      for(int iCat = 0; iCat < kNbMuCath; iCat++){

	//chi2 of global fit
	if(iCat == 0){//global muons
	  sprintf(name, "hChi2GlobalFit_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	  hChi2GlobalFit[iS][iSet][iCharge] = (TH1F *) gDirectory->Get(name);
	  sprintf(name, "%s_%s", hChi2GlobalFit[iS][iSet][iCharge]->GetName(), sBGName[iS]);
	  hChi2GlobalFit[iS][iSet][iCharge]->SetName(name);
	}
	//d0
	sprintf(name, "hD0_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hD0[iS][iSet][iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	sprintf(name, "%s_%s", hD0[iS][iSet][iCharge][iCat]->GetName(), sBGName[iS]);
	hD0[iS][iSet][iCharge][iCat]->SetName(name);
	//dz
	sprintf(name, "hDz_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hDz[iS][iSet][iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	sprintf(name, "%s_%s", hDz[iS][iSet][iCharge][iCat]->GetName(), sBGName[iS]);
	hDz[iS][iSet][iCharge][iCat]->SetName(name);

	if(iCat > 0){
	  //chi2 of tracker fit
	  sprintf(name, "hChi2TrackerFit_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	  hChi2TrackerFit[iS][iSet][iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	  sprintf(name, "%s_%s", hChi2TrackerFit[iS][iSet][iCharge][iCat]->GetName(), sBGName[iS]);
	  hChi2TrackerFit[iS][iSet][iCharge][iCat]->SetName(name);
	}
	//#hits in silicon tracker
	sprintf(name, "hNHitsSilicon_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hNHitsSilicon[iS][iSet][iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	sprintf(name, "%s_%s", hNHitsSilicon[iS][iSet][iCharge][iCat]->GetName(), sBGName[iS]);
	hNHitsSilicon[iS][iSet][iCharge][iCat]->SetName(name);
	hNHitsSilicon[iS][iSet][iCharge][iCat]->SetAxisRange(0., 34.);
	if(iS == 2) 
	  hNHitsSilicon[iS][iSet][iCharge][iCat]->Sumw2();

	if(iCat == 1){//tracker muons
	  
	  //hTM2DCompatibilityTight
	  sprintf(name, "hTM2DCompatibilityTight_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	  hTM2DCompatibilityTight[iS][iSet][iCharge] = (TH1F *) gDirectory->Get(name);
	  sprintf(name, "%s_%s", hTM2DCompatibilityTight[iS][iSet][iCharge]->GetName(), sBGName[iS]);
	  hTM2DCompatibilityTight[iS][iSet][iCharge]->SetName(name);
	  hTM2DCompatibilityTight[iS][iSet][iCharge]->SetMarkerStyle(marker[iS]);
	  //hTMLastStationOptimizedLowPtLoose
	  sprintf(name, "hTMLastStationOptimizedLowPtLoose_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	  hTMLastStationOptimizedLowPtLoose[iS][iSet][iCharge] = (TH1F *) gDirectory->Get(name);
	  sprintf(name, "%s_%s", hTMLastStationOptimizedLowPtLoose[iS][iSet][iCharge]->GetName(), sBGName[iS]);
	  hTMLastStationOptimizedLowPtLoose[iS][iSet][iCharge]->SetName(name);
	  hTMLastStationOptimizedLowPtLoose[iS][iSet][iCharge]->SetMarkerStyle(marker[iS]);
	}
	//calo compatibility
	sprintf(name, "hCaloCompatibility_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hCaloCompatibility[iS][iSet][iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	sprintf(name, "%s_%s", hCaloCompatibility[iS][iSet][iCharge][iCat]->GetName(), sBGName[iS]);
	hCaloCompatibility[iS][iSet][iCharge][iCat]->SetName(name);
	if(iS == 2) 
	  hCaloCompatibility[iS][iSet][iCharge][iCat]->Sumw2();

	//eta
	sprintf(name, "hMuon_eta_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hMuon_eta[iS][iSet][iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	sprintf(name, "%s_%s", hMuon_eta[iS][iSet][iCharge][iCat]->GetName(), sBGName[iS]);
	hMuon_eta[iS][iSet][iCharge][iCat]->SetName(name);
	//pT
	sprintf(name, "hMuon_pT_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hMuon_pT[iS][iSet][iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	sprintf(name, "%s_%s", hMuon_pT[iS][iSet][iCharge][iCat]->GetName(), sBGName[iS]);
	hMuon_pT[iS][iSet][iCharge][iCat]->SetName(name);
	//pT vs eta
	sprintf(name, "hMuon_pTvsEta_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hMuon_pT_eta[iS][iSet][iCharge][iCat] = (TH2F *) gDirectory->Get(name);
	sprintf(name, "%s_%s", hMuon_pT_eta[iS][iSet][iCharge][iCat]->GetName(), sBGName[iS]);
	hMuon_pT_eta[iS][iSet][iCharge][iCat]->SetName(name);
      }//muon category
    }
  }
  
  //   for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
  for(int iCharge = 0; iCharge < 1; iCharge++){
    for(int iCat = 0; iCat < kNbCath; iCat++){
      for(int iSet = 0; iSet < kNbDimuonSet; iSet++){
	//dimuon vertexing probability
	sprintf(name, "hDimuVtxProb_%d_%s_%s", iSet, chargeName[iCharge], oniaCatName[iCat]);
	hDimuVtxProb[iS][iSet][iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	sprintf(name, "%s_%s", hDimuVtxProb[iS][iSet][iCharge][iCat]->GetName(), sBGName[iS]);
	hDimuVtxProb[iS][iSet][iCharge][iCat]->SetName(name);
      }
    }
  }
}
