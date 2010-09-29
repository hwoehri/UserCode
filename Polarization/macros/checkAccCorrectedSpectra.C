#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

#include "TGraphAsymmErrors.h"
#include "TH1D.h"

Double_t minAcc = 0.001;//default: 0.1
Double_t fracMaxAcc = 0.1;
Double_t fracMax = 0.01;
Double_t minBinContent = 5.;

Int_t nbBinsCosTheta, nbBinsPhi, nBins1D; //will be filled in Get1DHistoFrom2D

//polarization histos from Data:
TH1D *Reco_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbPolVar];
TH1D *Reco_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1][jpsi::kNbPolVar];
TH1D *Reco_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1][jpsi::kNbPolVar];
TH2D *Reco2D_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1];
TH2D *Reco2D_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1];
TH2D *Reco2D_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];

//acceptance histos
TGraphAsymmErrors *gAcc_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbPolVar];
TGraphAsymmErrors *gAcc_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1][jpsi::kNbPolVar];
TGraphAsymmErrors *gAcc_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1][jpsi::kNbPolVar];
TH2D *hAcc2D_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1];
TH2D *hAcc2D_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1];
TH2D *hAcc2D_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TGraphAsymmErrors *gAcc2D_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];

//data corrected for acc:
TH1D *hData_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbPolVar];
TH1D *hData_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1][jpsi::kNbPolVar];
TH1D *hData_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1][jpsi::kNbPolVar];
TH2D *hData2D_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1];
TH2D *hData2D_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1];
TH2D *hData2D_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];

//mega histo containing the N phi bins appended:
TH1D *hData1D_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];

TH2D *thisHist;
TH2D *thisReco, *thisAcc;
TH1D *thisHist1D;

void ReadData(Char_t *fileNameData, Int_t rebinCosTh, Int_t rebinPhi);
void PlotUncorrData(Int_t iFrame, Char_t *polTag);
void PlotUncorrData2D(Int_t iFrame, Char_t *polTag, Char_t *oniaLabel);
void PlotCorrData2D(Int_t iFrame, Char_t *polTag, Char_t *oniaLabel);
void PlotUnCorr2D_OneByOne(Int_t iFrame, Char_t *label, Char_t *oniaLabel, Int_t rapBin, Int_t pTBin);
void PlotCorr2D_OneByOne(Int_t iFrame, Char_t *label, Char_t *oniaLabel, Int_t rapBin, Int_t pTBin);
void ReadAccHistos(Char_t *fileNameMC, Int_t rebinCosTh, Int_t rebinPhi, Bool_t normalising);
void CorrectForAcc(Char_t *polTag);
void Get1DHistoFrom2D();
void PlotHistos(Int_t iFrame, Char_t *polTag);
void PlotAll(Char_t *polTag);
void PlotRapIntegrated(Char_t *label, Int_t pTMin, Int_t pTMax);
void SaveCorrData(Char_t *polTag);
void AccCorrect(TH1D *hist, TGraphAsymmErrors *gAcc);
//======================================
//usage: root checkAccCorrectedSpectra.C+ or
//       root 'checkAccCorrectedSpectra.C+("HLT_Mu0Track0Jpsi", 1, 1)' (e.g.)
//======================================
void checkAccCorrectedSpectra(Char_t *hltTag = "HLT_Mu0TkMu0Jpsi",
			      Char_t *tag = "27Sep2010",
			      Int_t rebinCosTh = 2, //histos will be rebinned by "rebinCosTh"
			      Int_t rebinPhi = 2, //histos will be rebinned by "rebinPhi"
			      Bool_t normalising = kTRUE, //normalising Acc maps to 1
			      Char_t *oniaLabel = "J/#psi"){

  Char_t fileNameMC[200];
  //sprintf(fileNameMC, "/home/hermine/CMS/Work/Polarization/Florian/23Aug2010/files/accHistos_%s.root", hltTag);
  //sprintf(fileNameMC, "accHistos_HLT_Mu0Track0Jpsi_cut_30Aug2010.root");
  sprintf(fileNameMC, "accHistos_P_all_HLT_Mu0Track0Jpsi_26Sep2010.root");

  Char_t label[100];
  sprintf(label, "%s_%s", hltTag, tag);

  Char_t fileNameData[200];
  sprintf(fileNameData, "pol_data_all_%s.root", label);

  ReadData(fileNameData, rebinCosTh, rebinPhi);
  rebinCosTh = 1; rebinPhi = 1;
  ReadAccHistos(fileNameMC, rebinCosTh, rebinPhi, normalising);
  CorrectForAcc(hltTag); //calls internally "Get1DHistoFrom2D" to fill
                         //the histo hData1D_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];

  // //1D histos (cosTheta and phi):
  // // PlotUncorrData(jpsi::CS, hltTag); //integrated in rapidity
  // // PlotUncorrData(jpsi::HX, hltTag); //integrated in rapidity
  // // PlotHistos(jpsi::CS, hltTag); //1D plots are not acceptance corrected!!!
  // // PlotHistos(jpsi::HX, hltTag); //1D plots are not acceptance corrected!!!

  PlotAll(hltTag);
  // Int_t pTMin = 3, pTMax = 5;
  // PlotRapIntegrated(hltTag, pTMin, pTMax);

  // //2D histos (cosTheta and phi):
  // PlotUncorrData2D(jpsi::CS, hltTag, oniaLabel);
  // PlotUncorrData2D(jpsi::HX, hltTag, oniaLabel);
  // PlotCorrData2D(jpsi::CS, hltTag, oniaLabel);
  // PlotCorrData2D(jpsi::HX, hltTag, oniaLabel);
  // Int_t pTBin = 4;
  // for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
  //   PlotUnCorr2D_OneByOne(jpsi::CS, hltTag, oniaLabel, iRapBin, pTBin);
  // PlotUnCorr2D_OneByOne(jpsi::HX, hltTag, oniaLabel, iRapBin, pTBin);
  //   PlotCorr2D_OneByOne(jpsi::CS, hltTag, oniaLabel, iRapBin, pTBin);
  //   PlotCorr2D_OneByOne(jpsi::HX, hltTag, oniaLabel, iRapBin, pTBin);
  // }
}

//======================================
void Get1DHistoFrom2D(){

  Char_t name[100];
  Double_t content, error;
  Double_t contentReco, acc, maxAcc;
  Int_t binID = 0;
  nbBinsCosTheta = hData2D_pol_pT_rap[0][1][1]->GetNbinsX();
  nbBinsPhi = hData2D_pol_pT_rap[0][1][1]->GetNbinsY();
  nBins1D = nbBinsCosTheta * nbBinsPhi;

  printf("cosTheta consists of %d bins and phi of %d\n", nbBinsCosTheta, nbBinsPhi);
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
	//
	sprintf(name, "hData1D_Onia_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH1D(name, "", nBins1D, 0., (Double_t) nBins1D);
	binID = 0;
// 	printf("\n\n\n%s, pT %d, rap %d\n", jpsi::frameLabel[iFrame], iPTBin, iRapBin);

	maxAcc = hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetMaximum();
	for(int iBinY = 1; iBinY <= hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY(); iBinY++){
	  for(int iBinX = 1; iBinX <= hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(); iBinX++){
	    binID++;
// 	    printf("mapping binCosTh = %d, binPhi = %d into bin %d\n", iBinX, iBinY, binID);
	    content = hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iBinX,iBinY);
	    error = hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinError(iBinX,iBinY);

	    //flag the cells that have little acceptance or few binContents:
	    contentReco = Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iBinX,iBinY);
	    acc = hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iBinX,iBinY);

// 	    if(contentReco < minBinContent || acc < fracMaxAcc*maxAcc){
 	    // if(contentReco < minBinContent || acc < minAcc){
	    //   content = -1.;
	    //   error = 1.;
	    //   hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(binID, content);
	    //   hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinError(binID, error);
	    // }
	    // else{
	      hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(binID, content);
	      hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinError(binID, error);
	    // }
	  }
	}
      }
    }
  }
}

//======================================
void PlotUncorrData(Int_t iFrame, Char_t *label){

  gStyle->SetOptStat(10);
  Char_t name[100], title[100];
  //===========================================
  //cosTheta for different pT bins
  //===========================================
  sprintf(name, "c20CosTh_%s", jpsi::frameLabel[iFrame]);
  TCanvas *c20CosTh = new TCanvas(name, "cosTheta for pT bins", 900, 700);
  TH1F *hFrame20 = gPad->DrawFrame(-1., 0., 1., 1.3*Reco_pol_pT[iFrame][0][jpsi::cosThPol]->GetMaximum());
  sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);
  hFrame20->SetXTitle(name);
  sprintf(name, "Acc #times dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);
  hFrame20->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//   for(int iPTBin = 0; iPTBin < 3; iPTBin++){
    Reco_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->Draw("psame");
  }
  TLegend *leg20a = new TLegend(0.65,0.735119,0.9866071,0.985119);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    if(iPTBin == 0) sprintf(name, "all p_{T}");
    else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
    else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
    leg20a->AddEntry(Reco_pol_pT[iFrame][iPTBin][jpsi::cosThPol], name, "p");
  }
  leg20a->SetTextSize(0.035); leg20a->SetFillColor(0);
  leg20a->Draw();
 
  sprintf(name, "Figures/reco_cosTheta_%s_%s_pTBins.eps", jpsi::frameLabel[iFrame], label);  c20CosTh->Print(name);
  sprintf(name, "Figures/reco_cosTheta_%s_%s_pTBins.pdf", jpsi::frameLabel[iFrame], label);  c20CosTh->Print(name);
  sprintf(name, "Figures/reco_cosTheta_%s_%s_pTBins.gif", jpsi::frameLabel[iFrame], label);  c20CosTh->Print(name);

  //===========================================
  //phi for different pT bins
  //===========================================
  sprintf(name, "c20Phi_%s", jpsi::frameLabel[iFrame]);
  TCanvas *c20Phi = new TCanvas(name, "phi for pT bins", 900, 700);
  TH1F *hFrame20b = gPad->DrawFrame(0., 0., 360., 1.3*Reco_pol_pT[iFrame][0][jpsi::phiPol]->GetMaximum());
  sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);
  hFrame20b->SetXTitle(name);
  sprintf(name, "Acc #times dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);
  hFrame20b->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//   for(int iPTBin = 0; iPTBin < 3; iPTBin++){
    Reco_pol_pT[iFrame][iPTBin][jpsi::phiPol]->Draw("psame");
  }
//   TLegend *leg20b = new TLegend(0.1417411,0.1339286,0.5178571,0.3839286);
//   for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//     if(iPTBin == 0) sprintf(name, "all p_{T}");
//     else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
//     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
//     leg20b->AddEntry(Reco_pol_pT[iFrame][iPTBin][jpsi::phiPol], name, "p");
//   }
//   leg20b->SetTextSize(0.035); leg20b->SetFillColor(0);
//   leg20b->Draw();

  sprintf(name, "Figures/reco_phi_%s_%s_pTBins.eps", jpsi::frameLabel[iFrame], label);  c20Phi->Print(name);
  sprintf(name, "Figures/reco_phi_%s_%s_pTBins.pdf", jpsi::frameLabel[iFrame], label);  c20Phi->Print(name);
  sprintf(name, "Figures/reco_phi_%s_%s_pTBins.gif", jpsi::frameLabel[iFrame], label);  c20Phi->Print(name);

}

//=================================
void PlotUncorrData2D(Int_t iFrame, Char_t *label, Char_t *oniaLabel){

  Char_t name[100], title[100];
  TCanvas *c2D[jpsi::kNbRapForPTBins+1];
  for(int iRap = 1; iRap <= jpsi::kNbRapForPTBins; iRap++){
    sprintf(name, "c2D_%s_rap%d", jpsi::frameLabel[iFrame], iRap);
    sprintf(title, "phi vs cosTheta for pT bins, rap = %d (%s)", iRap, jpsi::frameLabel[iFrame]);
    c2D[iRap] = new TCanvas(name, title, 1000, 700);
    c2D[iRap]->Divide(3,2);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      c2D[iRap]->cd(iPTBin);
      if(iPTBin == jpsi::kNbPTBins) 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", oniaLabel, jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPTBin-1]);
      else 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", oniaLabel, jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
      Reco2D_pol_pT_rap[iFrame][iPTBin][iRap]->SetTitle(name);
//       Reco2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("colz");
      Reco2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("lego2");
    }
    sprintf(name, "Figures/reco2D_%s_%s_rap%d_pTBins.eps", jpsi::frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
    sprintf(name, "Figures/reco2D_%s_%s_rap%d_pTBins.pdf", jpsi::frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
    sprintf(name, "Figures/reco2D_%s_%s_rap%d_pTBins.png", jpsi::frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
  }
}



// //=================================
void PlotUnCorr2D_OneByOne(Int_t iFrame, Char_t *label, Char_t *oniaLabel, Int_t rapBin, Int_t pTBin){

  //  gStyle->SetOptStat(0);

  Char_t name[100], title[100];
  TCanvas *c2D[jpsi::kNbRapForPTBins+1];
  for(int iRap = rapBin; iRap < rapBin+1; iRap++){
    sprintf(name, "c2D_2_%s_rap%d", jpsi::frameLabel[iFrame], iRap);
    sprintf(title, "phi vs cosTheta for pT bins, rap = %d (%s)", iRap, jpsi::frameLabel[iFrame]);
    c2D[iRap] = new TCanvas(name, title, 1000, 700);
    for(int iPTBin = pTBin; iPTBin < pTBin+1; iPTBin++){
      if(iPTBin == jpsi::kNbPTBins) 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", oniaLabel, jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPTBin-1]);
      else 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", oniaLabel, jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
      Reco2D_pol_pT_rap[iFrame][iPTBin][iRap]->SetTitle(name);
      //       Reco2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("colz");
      Reco2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("lego2");
      sprintf(name, "Figures/reco2D_%s_%s_rap%d_pT%d.eps", jpsi::frameLabel[iFrame], label, iRap, iPTBin);  c2D[iRap]->Print(name);
      sprintf(name, "Figures/reco2D_%s_%s_rap%d_pT%d.pdf", jpsi::frameLabel[iFrame], label, iRap, iPTBin);  c2D[iRap]->Print(name);
      sprintf(name, "Figures/reco2D_%s_%s_rap%d_pT%d.png", jpsi::frameLabel[iFrame], label, iRap, iPTBin);  c2D[iRap]->Print(name);
    }
  }
}

//======================================
void PlotHistos(Int_t iFrame, Char_t *label){

  Char_t name[100];
  //===========================================
  //cosTheta for different pT bins
  //===========================================
  sprintf(name, "c10CosTh_%s", jpsi::frameLabel[iFrame]);
  TCanvas *c10CosTh = new TCanvas(name, "cosTheta for pT bins", 900, 700);
  TH1F *hFrame10 = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT[iFrame][0][jpsi::cosThPol]->GetMaximum());
  sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);
  hFrame10->SetXTitle(name);
  sprintf(name, "dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);
  hFrame10->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//   for(int iPTBin = 0; iPTBin < 3; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->Draw("psame");
  }
  TLegend *leg10a = new TLegend(0.6104911,0.735119,0.9866071,0.985119);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    if(iPTBin == 0) sprintf(name, "all p_{T}");
    else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
    else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
    leg10a->AddEntry(hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol], name, "p");
  }
  leg10a->SetTextSize(0.035); leg10a->SetFillColor(0);
  leg10a->Draw();

 
  sprintf(name, "Figures/dataCorr_cosTheta_%s_%s_pTBins.eps", jpsi::frameLabel[iFrame], label);  c10CosTh->Print(name);
  sprintf(name, "Figures/dataCorr_cosTheta_%s_%s_pTBins.pdf", jpsi::frameLabel[iFrame], label);  c10CosTh->Print(name);
  sprintf(name, "Figures/dataCorr_cosTheta_%s_%s_pTBins.gif", jpsi::frameLabel[iFrame], label);  c10CosTh->Print(name);

  //===========================================
  //phi for different pT bins
  //===========================================
  sprintf(name, "c10Phi_%s", jpsi::frameLabel[iFrame]);
  TCanvas *c10Phi = new TCanvas(name, "phi for pT bins", 900, 700);
  TH1F *hFrame10b = gPad->DrawFrame(0., 0., 360., 1.3*hData_pol_pT[iFrame][0][jpsi::phiPol]->GetMaximum());
  sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);
  hFrame10b->SetXTitle(name);
  sprintf(name, "dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);
  hFrame10b->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//   for(int iPTBin = 0; iPTBin < 3; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][jpsi::phiPol]->Draw("psame");
  }
//   TLegend *leg10b = new TLegend(0.1417411,0.1339286,0.5178571,0.3839286);
//   for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//     if(iPTBin == 0) sprintf(name, "all p_{T}");
//     else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
//     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
//     leg10b->AddEntry(hData_pol_pT[iFrame][iPTBin][jpsi::phiPol], name, "p");
//   }
//   leg10b->SetTextSize(0.035); leg10b->SetFillColor(0);
//   leg10b->Draw();

  sprintf(name, "Figures/dataCorr_phi_%s_%s_pTBins.eps", jpsi::frameLabel[iFrame], label);  c10Phi->Print(name);
  sprintf(name, "Figures/dataCorr_phi_%s_%s_pTBins.pdf", jpsi::frameLabel[iFrame], label);  c10Phi->Print(name);
  sprintf(name, "Figures/dataCorr_phi_%s_%s_pTBins.gif", jpsi::frameLabel[iFrame], label);  c10Phi->Print(name);
}

//======================================
void PlotAll(Char_t *label){

  Char_t name[100];
  //=====================================
  //rapidity integrated spectra:
  //raw data
  //=====================================
  sprintf(name, "c9All_%s", label);
  TCanvas *c9All = new TCanvas(name, "raw spectra for pT bins (all y)", 1000, 700);
  c9All->Divide(2,2);
  c9All->cd(1);
  Int_t iFrame = 0;
  TH1F *hFrame9 = gPad->DrawFrame(-1., 0., 1., 1.3*Reco_pol_pT[iFrame][0][jpsi::cosThPol]->GetMaximum());
  sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);  hFrame9->SetXTitle(name);
  sprintf(name, "Acc #times dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);  hFrame9->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    Reco_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->Draw("psame");
  }
  sprintf(name, "|y| < %1.1f", jpsi::rapForPTRange[jpsi::kNbRapForPTBins]);
  TLatex *tex9 = new TLatex(-0.9, 1.1*Reco_pol_pT[iFrame][0][jpsi::cosThPol]->GetMaximum(), name);
  tex9->SetTextSize(0.06); tex9->Draw();

  c9All->cd(2);
  TH1F *hFrame9b = gPad->DrawFrame(0., 0., 360., 1.3*Reco_pol_pT[iFrame][0][jpsi::phiPol]->GetMaximum());
  sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);  hFrame9b->SetXTitle(name);
  sprintf(name, "Acc #times dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);  hFrame9b->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    Reco_pol_pT[iFrame][iPTBin][jpsi::phiPol]->Draw("psame");
  }
  TLegend *leg9a = new TLegend(0.6171352,0.6922123,0.9936412,0.9929315);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    if(iPTBin == 0) sprintf(name, "all p_{T}");
    else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
    else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
    leg9a->AddEntry(Reco_pol_pT[iFrame][iPTBin][jpsi::cosThPol], name, "p");
  }
  leg9a->SetTextSize(0.035); leg9a->SetFillColor(0);
  leg9a->Draw();

  c9All->cd(3);
  iFrame = 1;
  TH1F *hFrame9c = gPad->DrawFrame(-1., 0., 1., 1.3*Reco_pol_pT[iFrame][0][jpsi::cosThPol]->GetMaximum());
  sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);  hFrame9c->SetXTitle(name);
  sprintf(name, "Acc #times dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);  hFrame9c->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    Reco_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->Draw("psame");
  }
  c9All->cd(4);
  TH1F *hFrame9d = gPad->DrawFrame(0., 0., 360., 1.3*Reco_pol_pT[iFrame][0][jpsi::phiPol]->GetMaximum());
  sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);  hFrame9d->SetXTitle(name);
  sprintf(name, "Acc #times dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);  hFrame9d->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    Reco_pol_pT[iFrame][iPTBin][jpsi::phiPol]->Draw("psame");
  }

  sprintf(name, "Figures/dataUncorr_%s_pTBins.eps", label);  c9All->Print(name);
  sprintf(name, "Figures/dataUncorr_%s_pTBins.pdf", label);  c9All->Print(name);
  sprintf(name, "Figures/dataUncorr_%s_pTBins.gif", label);  c9All->Print(name);

  //=====================================
  //rapidity integrated spectra:
  //data corrected for acceptance
  //=====================================
  sprintf(name, "c10All_%s", label);
  TCanvas *c10All = new TCanvas(name, "corrected spectra for pT bins (all y)", 1000, 700);
  c10All->Divide(2,2);
  c10All->cd(1);
  TH1F *hFrame10 = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT[iFrame][0][jpsi::cosThPol]->GetMaximum());
  sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);  hFrame10->SetXTitle(name);
  sprintf(name, "dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);  hFrame10->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->Draw("psame");
  }
  sprintf(name, "|y| < %1.1f", jpsi::rapForPTRange[jpsi::kNbRapForPTBins]);
  TLatex *tex10 = new TLatex(-0.9, 1.1*hData_pol_pT[iFrame][0][jpsi::cosThPol]->GetMaximum(), name);
  tex10->SetTextSize(0.06); tex10->Draw();
  c10All->cd(2);
  TH1F *hFrame10b = gPad->DrawFrame(0., 0., 360., 1.3*hData_pol_pT[iFrame][0][jpsi::phiPol]->GetMaximum());
  sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);  hFrame10b->SetXTitle(name);
  sprintf(name, "dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);  hFrame10b->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][jpsi::phiPol]->Draw("psame");
  }
  TLegend *leg10a = new TLegend(0.6171352,0.6922123,0.9936412,0.9929315);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    if(iPTBin == 0) sprintf(name, "all p_{T}");
    else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
    else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
    leg10a->AddEntry(hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol], name, "p");
  }
  leg10a->SetTextSize(0.035); leg10a->SetFillColor(0);
  leg10a->Draw();

  c10All->cd(3);
  iFrame = 1;
  TH1F *hFrame10c = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT[iFrame][0][jpsi::cosThPol]->GetMaximum());
  sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);  hFrame10c->SetXTitle(name);
  sprintf(name, "dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);  hFrame10c->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->Draw("psame");
  }
  c10All->cd(4);
  TH1F *hFrame10d = gPad->DrawFrame(0., 0., 360., 1.3*hData_pol_pT[iFrame][0][jpsi::phiPol]->GetMaximum());
  sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);  hFrame10d->SetXTitle(name);
  sprintf(name, "dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);  hFrame10d->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][jpsi::phiPol]->Draw("psame");
  }

  sprintf(name, "Figures/dataCorr_%s_pTBins.eps", label);  c10All->Print(name);
  sprintf(name, "Figures/dataCorr_%s_pTBins.pdf", label);  c10All->Print(name);
  sprintf(name, "Figures/dataCorr_%s_pTBins.gif", label);  c10All->Print(name);

  //=====================================
  //rapidity differential spectra
  //before acceptance correction
  //=====================================
  Char_t title[100];
  TCanvas *c11All[jpsi::kNbRapForPTBins+1];
  Double_t maxY = 0;

  for(int iRap = 1; iRap < jpsi::kNbRapForPTBins+1; iRap++){
    sprintf(name, "c11All_%s_rap%d", label, iRap);
    sprintf(title, "raw spectra for pT bins (rap %d)", iRap);
    c11All[iRap] = new TCanvas(name, title, 1000, 700);
    c11All[iRap]->Divide(2,2);
    c11All[iRap]->cd(1);

    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->GetMaximum();

    TH1F *hFrame10 = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
    sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::CS]);  hFrame10->SetXTitle(name);
    sprintf(name, "Acc #times dN/d(cos#theta_{%s})", jpsi::frameLabel[jpsi::CS]);  hFrame10->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

    TLegend *leg10a = new TLegend(0.7,0.6922123,0.9936412,0.9929315);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      if(iPTBin == 0) sprintf(name, "all p_{T}");
      else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
      else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
      leg10a->AddEntry(Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol], name, "p");
    }
    leg10a->SetTextSize(0.035); leg10a->SetFillColor(0);
    leg10a->Draw();

    if(iRap == 1)
      sprintf(name, "|y| < %1.1f", jpsi::rapForPTRange[iRap]);
    else
      sprintf(name, "%1.1f < |y| < %1.1f", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);

    TLatex *tex11a = new TLatex(-0.9, 1.1*maxY, name);
    tex11a->SetTextSize(0.06); tex11a->Draw();
    
    c11All[iRap]->cd(2);
    // gPad->SetLogy(); 
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->GetMaximum();

    TH1F *hFrame11b = gPad->DrawFrame(0., 0., 360., 1.3*maxY);
    sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::CS]);  hFrame11b->SetXTitle(name);
    sprintf(name, "Acc #times dN/d#phi_{%s}", jpsi::frameLabel[jpsi::CS]);  hFrame11b->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->Draw("psame");

    c11All[iRap]->cd(3);
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->GetMaximum();
    TH1F *hFrame11c = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
    sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::HX]);  hFrame11c->SetXTitle(name);
    sprintf(name, "Acc #times dN/d(cos#theta_{%s})", jpsi::frameLabel[jpsi::HX]);  hFrame11c->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

    c11All[iRap]->cd(4);
    // gPad->SetLogy(); 
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->GetMaximum();
    TH1F *hFrame11d = gPad->DrawFrame(0., 0., 360., 1.3*maxY);
    sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::HX]);  hFrame11d->SetXTitle(name);
    sprintf(name, "Acc #times dN/d#phi_{%s}", jpsi::frameLabel[jpsi::HX]);  hFrame11d->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->Draw("psame");
    }

    sprintf(name, "Figures/dataUncorr_%s_rap%d_pTBins.eps", label, iRap);  c11All[iRap]->Print(name);
    sprintf(name, "Figures/dataUncorr_%s_rap%d_pTBins.pdf", label, iRap);  c11All[iRap]->Print(name);
    sprintf(name, "Figures/dataUncorr_%s_rap%d_pTBins.gif", label, iRap);  c11All[iRap]->Print(name);
  }

  //=====================================
  //rapidity differential spectra
  //before acceptance correction:
  //each distribution on a different canvas
  //=====================================
  Char_t rapLabel[100];
  TCanvas *c11Alla[jpsi::kNbRapForPTBins+1][4];
  for(int iRap = 1; iRap < jpsi::kNbRapForPTBins+1; iRap++){

    TLegend *leg10a = new TLegend(0.7,0.6922123,0.9936412,0.9929315);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      if(iPTBin == 0) sprintf(name, "all p_{T}");
      else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
      else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
      leg10a->AddEntry(Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol], name, "p");
    }
    leg10a->SetTextSize(0.035); leg10a->SetFillColor(0);

    if(iRap == 1)
      sprintf(rapLabel, "|y| < %1.1f", jpsi::rapForPTRange[iRap]);
    else
      sprintf(rapLabel, "%1.1f < |y| < %1.1f", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);

    TLatex *tex11a = new TLatex(-0.9, 1.1*maxY, rapLabel);
    tex11a->SetTextSize(0.06);
    tex11a->Draw();

    //1.) cosTheta in the CS frame
    sprintf(name, "c11Alla_%s_rap%d", label, iRap);
    sprintf(title, "raw spectra for pT bins (rap %d)", iRap);
    c11Alla[iRap][0] = new TCanvas(name, title);

    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->GetMaximum();

    TH1F *hFrame10 = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
    sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::CS]);  hFrame10->SetXTitle(name);
    sprintf(name, "Acc #times dN/d(cos#theta_{%s})", jpsi::frameLabel[jpsi::CS]);  hFrame10->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

    leg10a->DrawClone();
    tex11a->DrawLatex(-0.9, 1.1*maxY, rapLabel);

    sprintf(name, "Figures/dataUncorr_%s_%s_cosTh_rap%d_pTBins.eps", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c11Alla[iRap][0]->Print(name);
    sprintf(name, "Figures/dataUncorr_%s_%s_cosTh_rap%d_pTBins.pdf", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c11Alla[iRap][0]->Print(name);
    sprintf(name, "Figures/dataUncorr_%s_%s_cosTh_rap%d_pTBins.gif", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c11Alla[iRap][0]->Print(name);

    //
    //2.) phi in the CS frame
    c11Alla[iRap][1] = new TCanvas(name, title);    
    // gPad->SetLogy(); 
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->GetMaximum();

    TH1F *hFrame11b = gPad->DrawFrame(0., 0., 360., 1.3*maxY);
    sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::CS]);  hFrame11b->SetXTitle(name);
    sprintf(name, "Acc #times dN/d#phi_{%s}", jpsi::frameLabel[jpsi::CS]);  hFrame11b->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->Draw("psame");

    leg10a->DrawClone();
    tex11a->DrawLatex(15., 1.1*maxY, rapLabel);

    sprintf(name, "Figures/dataUncorr_%s_%s_phi_rap%d_pTBins.eps", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c11Alla[iRap][1]->Print(name);	    
    sprintf(name, "Figures/dataUncorr_%s_%s_phi_rap%d_pTBins.pdf", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c11Alla[iRap][1]->Print(name);	    
    sprintf(name, "Figures/dataUncorr_%s_%s_phi_rap%d_pTBins.gif", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c11Alla[iRap][1]->Print(name);

    //3.) cosTheta in the HX frame
    c11Alla[iRap][2] = new TCanvas(name, title);    
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->GetMaximum();
    TH1F *hFrame11c = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
    sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::HX]);  hFrame11c->SetXTitle(name);
    sprintf(name, "Acc #times dN/d(cos#theta_{%s})", jpsi::frameLabel[jpsi::HX]);  hFrame11c->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

    leg10a->DrawClone();
    tex11a->DrawLatex(-0.9, 1.1*maxY, rapLabel);

    sprintf(name, "Figures/dataUncorr_%s_%s_cosTh_rap%d_pTBins.eps", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c11Alla[iRap][2]->Print(name);	    
    sprintf(name, "Figures/dataUncorr_%s_%s_cosTh_rap%d_pTBins.pdf", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c11Alla[iRap][2]->Print(name);	    
    sprintf(name, "Figures/dataUncorr_%s_%s_cosTh_rap%d_pTBins.gif", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c11Alla[iRap][2]->Print(name);

    //4.) phi in the HX frame
    c11Alla[iRap][3] = new TCanvas(name, title);    
    // gPad->SetLogy(); 
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->GetMaximum();
    TH1F *hFrame11d = gPad->DrawFrame(0., 0., 360., 1.3*maxY);
    sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::HX]);  hFrame11d->SetXTitle(name);
    sprintf(name, "Acc #times dN/d#phi_{%s}", jpsi::frameLabel[jpsi::HX]);  hFrame11d->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->Draw("psame");
    }

    leg10a->DrawClone();
    tex11a->DrawLatex(15., 1.1*maxY, rapLabel);

    sprintf(name, "Figures/dataUncorr_%s_%s_phi_rap%d_pTBins.eps", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c11Alla[iRap][3]->Print(name);
    sprintf(name, "Figures/dataUncorr_%s_%s_phi_rap%d_pTBins.pdf", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c11Alla[iRap][3]->Print(name);
    sprintf(name, "Figures/dataUncorr_%s_%s_phi_rap%d_pTBins.gif", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c11Alla[iRap][3]->Print(name);
  }

  //=============================================
  //rapidity differential cosTheta and phispectra
  //after acceptance correction
  //=============================================
  TCanvas *c12All[jpsi::kNbRapForPTBins+1];
  for(int iRap = 1; iRap < jpsi::kNbRapForPTBins+1; iRap++){
    sprintf(name, "c12All_%s_rap%d", label, iRap);
    sprintf(title, "corrected spectra for pT bins (rap %d)", iRap);
    c12All[iRap] = new TCanvas(name, title, 1000, 700);
    c12All[iRap]->Divide(2,2);
    c12All[iRap]->cd(1);
    Int_t iFrame = 0;
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::cosThPol]->GetMaximum() > maxY)
	maxY = hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::cosThPol]->GetMaximum();
    TH1F *hFrame12 = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
    sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);  hFrame12->SetXTitle(name);
    sprintf(name, "dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);  hFrame12->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

    TLegend *leg12a = new TLegend(0.7,0.6922123,0.9936412,0.9929315);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      if(iPTBin == 0) sprintf(name, "all p_{T}");
      else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
      else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
      leg12a->AddEntry(hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::cosThPol], name, "p");
    }
    leg12a->SetTextSize(0.035); leg12a->SetFillColor(0);
    leg12a->Draw();

    if(iRap == 1)
      sprintf(name, "|y| < %1.1f", jpsi::rapForPTRange[iRap]);
    else
      sprintf(name, "%1.1f < |y| < %1.1f", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);

    TLatex *tex12a = new TLatex(-0.9, 1.1*maxY, name);
    tex12a->SetTextSize(0.06); tex12a->Draw();
    
    c12All[iRap]->cd(2);
    // gPad->SetLogy();
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::phiPol]->GetMaximum() > maxY)
	maxY = hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::phiPol]->GetMaximum();

    TH1F *hFrame12b = gPad->DrawFrame(0., 0., 360., 1.3*maxY);
    sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);  hFrame12b->SetXTitle(name);
    sprintf(name, "dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);  hFrame12b->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::phiPol]->Draw("psame");

    c12All[iRap]->cd(3);
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::cosThPol]->GetMaximum() > maxY)
	maxY = hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::cosThPol]->GetMaximum();
    TH1F *hFrame12c = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
    sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);  hFrame12c->SetXTitle(name);
    sprintf(name, "dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);  hFrame12c->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

    c12All[iRap]->cd(4);
    // gPad->SetLogy(); 
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::phiPol]->GetMaximum() > maxY)
	maxY = hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::phiPol]->GetMaximum();
    TH1F *hFrame12d = gPad->DrawFrame(0., 0., 360., 1.3*maxY);
    sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);  hFrame12d->SetXTitle(name);
    sprintf(name, "dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);  hFrame12d->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::phiPol]->Draw("psame");
    }

    sprintf(name, "Figures/dataCorr_%s_rap%d_pTBins.eps", label, iRap);  c12All[iRap]->Print(name);
    sprintf(name, "Figures/dataCorr_%s_rap%d_pTBins.pdf", label, iRap);  c12All[iRap]->Print(name);
    sprintf(name, "Figures/dataCorr_%s_rap%d_pTBins.gif", label, iRap);  c12All[iRap]->Print(name);
  }

  //=====================================
  //rapidity differential spectra
  //after acceptance correction:
  //each distribution on a different canvas
  //=====================================
  TCanvas *c12Alla[jpsi::kNbRapForPTBins+1][4];
  for(int iRap = 1; iRap < jpsi::kNbRapForPTBins+1; iRap++){

    TLegend *leg10a = new TLegend(0.7,0.6922123,0.9936412,0.9929315);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      if(iPTBin == 0) sprintf(name, "all p_{T}");
      else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
      else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
      leg10a->AddEntry(hData_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol], name, "p");
    }
    leg10a->SetTextSize(0.035); leg10a->SetFillColor(0);

    if(iRap == 1)
      sprintf(rapLabel, "|y| < %1.1f", jpsi::rapForPTRange[iRap]);
    else
      sprintf(rapLabel, "%1.1f < |y| < %1.1f", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);
    TLatex *tex11a = new TLatex(-0.9, 1.1*maxY, rapLabel);
    tex11a->SetTextSize(0.06); 

    //1.) cosTheta in the CS frame
    sprintf(name, "c12Alla_%s_rap%d", label, iRap);
    sprintf(title, "corrected spectra for pT bins (rap %d)", iRap);
    c12Alla[iRap][0] = new TCanvas(name, title);

    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(hData_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->GetMaximum() > maxY)
	maxY = hData_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->GetMaximum();

    TH1F *hFrame10 = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
    sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::CS]);  hFrame10->SetXTitle(name);
    sprintf(name, "dN/d(cos#theta_{%s})", jpsi::frameLabel[jpsi::CS]);  hFrame10->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      hData_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

    leg10a->DrawClone();
    tex11a->DrawLatex(-0.9, 1.1*maxY, rapLabel);

    sprintf(name, "Figures/dataCorr_%s_%s_cosTh_rap%d_pTBins.eps", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c12Alla[iRap][0]->Print(name);
    sprintf(name, "Figures/dataCorr_%s_%s_cosTh_rap%d_pTBins.pdf", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c12Alla[iRap][0]->Print(name);
    sprintf(name, "Figures/dataCorr_%s_%s_cosTh_rap%d_pTBins.gif", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c12Alla[iRap][0]->Print(name);

    //
    //2.) phi in the CS frame
    c12Alla[iRap][1] = new TCanvas(name, title);    
    // gPad->SetLogy(); 
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(hData_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->GetMaximum() > maxY)
	maxY = hData_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->GetMaximum();

    TH1F *hFrame11b = gPad->DrawFrame(0., 0., 360., 1.3*maxY);
    sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::CS]);  hFrame11b->SetXTitle(name);
    sprintf(name, "dN/d#phi_{%s}", jpsi::frameLabel[jpsi::CS]);  hFrame11b->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      hData_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->Draw("psame");

    leg10a->DrawClone();
    tex11a->DrawLatex(15., 1.1*maxY, rapLabel);

    sprintf(name, "Figures/dataCorr_%s_%s_phi_rap%d_pTBins.eps", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c12Alla[iRap][1]->Print(name);	    
    sprintf(name, "Figures/dataCorr_%s_%s_phi_rap%d_pTBins.pdf", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c12Alla[iRap][1]->Print(name);	    
    sprintf(name, "Figures/dataCorr_%s_%s_phi_rap%d_pTBins.gif", label, jpsi::frameLabel[jpsi::CS], iRap);  
    c12Alla[iRap][1]->Print(name);

    //3.) cosTheta in the HX frame
    c12Alla[iRap][2] = new TCanvas(name, title);    
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(hData_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->GetMaximum() > maxY)
	maxY = hData_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->GetMaximum();
    TH1F *hFrame11c = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
    sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::HX]);  hFrame11c->SetXTitle(name);
    sprintf(name, "dN/d(cos#theta_{%s})", jpsi::frameLabel[jpsi::HX]);  hFrame11c->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      hData_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

    leg10a->DrawClone();
    tex11a->DrawLatex(-0.9, 1.1*maxY, rapLabel);

    sprintf(name, "Figures/dataCorr_%s_%s_cosTh_rap%d_pTBins.eps", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c12Alla[iRap][2]->Print(name);	    
    sprintf(name, "Figures/dataCorr_%s_%s_cosTh_rap%d_pTBins.pdf", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c12Alla[iRap][2]->Print(name);	    
    sprintf(name, "Figures/dataCorr_%s_%s_cosTh_rap%d_pTBins.gif", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c12Alla[iRap][2]->Print(name);

    //4.) phi in the HX frame
    c12Alla[iRap][3] = new TCanvas(name, title);    
    // gPad->SetLogy(); 
    maxY = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
      if(hData_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->GetMaximum() > maxY)
	maxY = hData_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->GetMaximum();
    TH1F *hFrame11d = gPad->DrawFrame(0., 0., 360., 1.3*maxY);
    sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::HX]);  hFrame11d->SetXTitle(name);
    sprintf(name, "dN/d#phi_{%s}", jpsi::frameLabel[jpsi::HX]);  hFrame11d->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      hData_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->Draw("psame");
    }

    leg10a->DrawClone();
    tex11a->DrawLatex(15., 1.1*maxY, rapLabel);

    sprintf(name, "Figures/dataCorr_%s_%s_phi_rap%d_pTBins.eps", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c12Alla[iRap][3]->Print(name);
    sprintf(name, "Figures/dataCorr_%s_%s_phi_rap%d_pTBins.pdf", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c12Alla[iRap][3]->Print(name);
    sprintf(name, "Figures/dataCorr_%s_%s_phi_rap%d_pTBins.gif", label, jpsi::frameLabel[jpsi::HX], iRap);  
    c12Alla[iRap][3]->Print(name);
  }
}

//=====================================
void PlotRapIntegrated(Char_t *label, Int_t pTMin, Int_t pTMax){

  Char_t name[100];
  //=====================================
  //rapidity integrated spectra:
  //raw data
  //=====================================
  Double_t maxY = 0;

  sprintf(name, "c9All_%s", label);
  TCanvas *c9All = new TCanvas(name, "raw spectra for some pT bins (all y)", 1000, 700);
  c9All->Divide(2,2);
  c9All->cd(1);
  maxY = Reco_pol_pT[jpsi::CS][1][jpsi::cosThPol]->GetMaximum();
  TH1F *hFrame9 = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
  sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::CS]);  hFrame9->SetXTitle(name);
  sprintf(name, "Acc #times dN/d(cos#theta_{%s})", jpsi::frameLabel[jpsi::CS]);  hFrame9->SetYTitle(name);
  for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++){
    Reco_pol_pT[jpsi::CS][iPTBin][jpsi::cosThPol]->Draw("psame");
  }
  sprintf(name, "|y| < %1.1f", jpsi::rapForPTRange[jpsi::kNbRapForPTBins]);
  TLatex *tex9 = new TLatex(-0.9, 1.1*maxY, name);
  tex9->SetTextSize(0.06); tex9->Draw();

  c9All->cd(2);
  maxY = Reco_pol_pT[jpsi::CS][2][jpsi::phiPol]->GetMaximum();
  TH1F *hFrame9b = gPad->DrawFrame(0., 0., 360., 1.3*maxY);
  sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::CS]);  hFrame9b->SetXTitle(name);
  sprintf(name, "Acc #times dN/d#phi_{%s}", jpsi::frameLabel[jpsi::CS]);  hFrame9b->SetYTitle(name);
  for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++){
    Reco_pol_pT[jpsi::CS][iPTBin][jpsi::phiPol]->Draw("psame");
  }
  TLegend *leg9a = new TLegend(0.6171352,0.75,0.9936412,0.9929315);
  for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++){
    if(iPTBin == 0) sprintf(name, "all p_{T}");
    else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
    else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
    leg9a->AddEntry(Reco_pol_pT[jpsi::CS][iPTBin][jpsi::cosThPol], name, "p");
  }
  leg9a->SetTextSize(0.035); leg9a->SetFillColor(0);
  leg9a->Draw();

  c9All->cd(3);
  maxY = Reco_pol_pT[jpsi::HX][2][jpsi::cosThPol]->GetMaximum();
  TH1F *hFrame9c = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
  sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::HX]);  hFrame9c->SetXTitle(name);
  sprintf(name, "Acc #times dN/d(cos#theta_{%s})", jpsi::frameLabel[jpsi::HX]);  hFrame9c->SetYTitle(name);
  for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++){
    Reco_pol_pT[jpsi::HX][iPTBin][jpsi::cosThPol]->Draw("psame");
  }
  c9All->cd(4);
  maxY = Reco_pol_pT[jpsi::HX][2][jpsi::phiPol]->GetMaximum();
  TH1F *hFrame9d = gPad->DrawFrame(0., 0., 360., 1.3*maxY);
  sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::HX]);  hFrame9d->SetXTitle(name);
  sprintf(name, "Acc #times dN/d#phi_{%s}", jpsi::frameLabel[jpsi::HX]);  hFrame9d->SetYTitle(name);
  for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++){
    Reco_pol_pT[jpsi::HX][iPTBin][jpsi::phiPol]->Draw("psame");
  }

  sprintf(name, "Figures/dataUncorr_%s_pT%dTo%d.eps",     label, pTMin, pTMax);  c9All->Print(name);
  sprintf(name, "Figures/dataUncorr_%s_pT%dTo%d.pdf", label, pTMin, pTMax);  c9All->Print(name);
  sprintf(name, "Figures/dataUncorr_%s_pT%dTo%d.gif", label, pTMin, pTMax);  c9All->Print(name);

  // //=====================================
  // //rapidity integrated spectra:
  // //data corrected for acceptance
  // //=====================================
  // sprintf(name, "c10All_%s", label);
  // TCanvas *c10All = new TCanvas(name, "corrected spectra for pT bins (all y)", 1000, 700);
  // c10All->Divide(2,2);
  // c10All->cd(1);
  // TH1F *hFrame10 = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT[iFrame][0][jpsi::cosThPol]->GetMaximum());
  // sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);  hFrame10->SetXTitle(name);
  // sprintf(name, "dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);  hFrame10->SetYTitle(name);
  // for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  //   hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->Draw("psame");
  // }
  // sprintf(name, "|y| < %1.1f", jpsi::rapForPTRange[jpsi::kNbRapForPTBins]);
  // TLatex *tex10 = new TLatex(-0.9, 1.1*hData_pol_pT[iFrame][0][jpsi::cosThPol]->GetMaximum(), name);
  // tex10->SetTextSize(0.06); tex10->Draw();
  // c10All->cd(2);
  // TH1F *hFrame10b = gPad->DrawFrame(0., 0., 360., 1.3*hData_pol_pT[iFrame][0][jpsi::phiPol]->GetMaximum());
  // sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);  hFrame10b->SetXTitle(name);
  // sprintf(name, "dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);  hFrame10b->SetYTitle(name);
  // for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  //   hData_pol_pT[iFrame][iPTBin][jpsi::phiPol]->Draw("psame");
  // }
  // TLegend *leg10a = new TLegend(0.6171352,0.6922123,0.9936412,0.9929315);
  // for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  //   if(iPTBin == 0) sprintf(name, "all p_{T}");
  //   else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
  //   else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
  //   leg10a->AddEntry(hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol], name, "p");
  // }
  // leg10a->SetTextSize(0.035); leg10a->SetFillColor(0);
  // leg10a->Draw();

  // c10All->cd(3);
  // iFrame = 1;
  // TH1F *hFrame10c = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT[iFrame][0][jpsi::cosThPol]->GetMaximum());
  // sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);  hFrame10c->SetXTitle(name);
  // sprintf(name, "dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);  hFrame10c->SetYTitle(name);
  // for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  //   hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->Draw("psame");
  // }
  // c10All->cd(4);
  // TH1F *hFrame10d = gPad->DrawFrame(0., 0., 360., 1.3*hData_pol_pT[iFrame][0][jpsi::phiPol]->GetMaximum());
  // sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);  hFrame10d->SetXTitle(name);
  // sprintf(name, "dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);  hFrame10d->SetYTitle(name);
  // for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  //   hData_pol_pT[iFrame][iPTBin][jpsi::phiPol]->Draw("psame");
  // }

  // sprintf(name, "Figures/dataCorr_%s_pTBins.eps", label);  c10All->Print(name);
  // sprintf(name, "Figures/dataCorr_%s_pTBins.pdf", label);  c10All->Print(name);
  // sprintf(name, "Figures/dataCorr_%s_pTBins.gif", label);  c10All->Print(name);

  //=====================================
  //rapidity differential spectra
  //before acceptance correction
  //=====================================
  Char_t title[100];
  TCanvas *c11All[jpsi::kNbRapForPTBins+1];
  Int_t maxIndex[jpsi::kNbRapForPTBins+1] = {1, 4, 3, 3};
  for(int iRap = 1; iRap < jpsi::kNbRapForPTBins+1; iRap++){

    maxY = 0;
    for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->GetMaximum();

    sprintf(name, "c11All_%s_rap%d", label, iRap);
    sprintf(title, "raw spectra for pT bins (rap %d)", iRap);
    c11All[iRap] = new TCanvas(name, title, 1000, 700);
    c11All[iRap]->Divide(2,2);
    c11All[iRap]->cd(1);
    TH1F *hFrame10 = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
    sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::CS]);  hFrame10->SetXTitle(name);
    sprintf(name, "Acc #times dN/d(cos#theta_{%s})", jpsi::frameLabel[jpsi::CS]);  hFrame10->SetYTitle(name);
    for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++)
      Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

    TLegend *leg10a = new TLegend(0.7,0.75,0.9936412,0.9929315);
    for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++){
      if(iPTBin == 0) sprintf(name, "all p_{T}");
      else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
      else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
      leg10a->AddEntry(Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::cosThPol], name, "p");
    }
    leg10a->SetTextSize(0.035); leg10a->SetFillColor(0);
    leg10a->Draw();

    if(iRap == 1)
      sprintf(name, "|y| < %1.1f", jpsi::rapForPTRange[iRap]);
    else
      sprintf(name, "%1.1f < |y| < %1.1f", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);

    TLatex *tex11a = new TLatex(-0.9, 1.1*maxY, name);
    tex11a->SetTextSize(0.06); tex11a->Draw();
    
    c11All[iRap]->cd(2);
    // gPad->SetLogy(); 
    //    gPad->SetGridy();
    maxY = 0;
    for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->GetMaximum();

    TH1F *hFrame11b = gPad->DrawFrame(0., 1., 360., 1.3*maxY);
    sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::CS]);  hFrame11b->SetXTitle(name);
    sprintf(name, "Acc #times dN/d#phi_{%s}", jpsi::frameLabel[jpsi::CS]);  hFrame11b->SetYTitle(name);
    for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++)
      Reco_pol_pT_rap[jpsi::CS][iPTBin][iRap][jpsi::phiPol]->Draw("psame");

    c11All[iRap]->cd(3);
    maxY = 0;
    for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->GetMaximum();
    TH1F *hFrame11c = gPad->DrawFrame(-1., 0., 1., 1.3*maxY);
    sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::HX]);  hFrame11c->SetXTitle(name);
    sprintf(name, "Acc #times dN/d(cos#theta_{%s})", jpsi::frameLabel[jpsi::HX]);  hFrame11c->SetYTitle(name);
    for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++)
      Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

    c11All[iRap]->cd(4);
    // gPad->SetLogy(); 
    // gPad->SetGridy();
    maxY = 0;
    for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++)
      if(Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->GetMaximum() > maxY)
	maxY = Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->GetMaximum();
    TH1F *hFrame11d = gPad->DrawFrame(0., 1., 360., 1.3*maxY);
    sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::HX]);  hFrame11d->SetXTitle(name);
    sprintf(name, "Acc #times dN/d#phi_{%s}", jpsi::frameLabel[jpsi::HX]);  hFrame11d->SetYTitle(name);
    for(int iPTBin = pTMin; iPTBin < pTMax+1; iPTBin++){
      Reco_pol_pT_rap[jpsi::HX][iPTBin][iRap][jpsi::phiPol]->Draw("psame");
    }

    sprintf(name, "Figures/dataUncorr_%s_rap%d_pT%dTo%d.eps", label, iRap, pTMin, pTMax);  c11All[iRap]->Print(name);
    sprintf(name, "Figures/dataUncorr_%s_rap%d_pT%dTo%d.pdf", label, iRap, pTMin, pTMax);  c11All[iRap]->Print(name);
    sprintf(name, "Figures/dataUncorr_%s_rap%d_pT%dTo%d.gif", label, iRap, pTMin, pTMax);  c11All[iRap]->Print(name);
  }

  // //=====================================
  // //rapidity differential spectra:
  // //corrected for acceptance
  // //=====================================
  // TCanvas *c12All[jpsi::kNbRapForPTBins+1];
  // for(int iRap = 1; iRap < jpsi::kNbRapForPTBins+1; iRap++){
  //   sprintf(name, "c12All_%s_rap%d", label, iRap);
  //   sprintf(title, "corrected spectra for pT bins (rap %d)", iRap);
  //   c12All[iRap] = new TCanvas(name, title, 1000, 700);
  //   c12All[iRap]->Divide(2,2);
  //   c12All[iRap]->cd(1);
  //   Int_t iFrame = 0;
  //   TH1F *hFrame12 = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT_rap[iFrame][maxIndex[iRap]][iRap][jpsi::cosThPol]->GetMaximum());
  //   sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);  hFrame12->SetXTitle(name);
  //   sprintf(name, "dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);  hFrame12->SetYTitle(name);
  //   for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
  //     hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

  //   TLegend *leg12a = new TLegend(0.7,0.6922123,0.9936412,0.9929315);
  //   for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  //     if(iPTBin == 0) sprintf(name, "all p_{T}");
  //     else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
  //     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
  //     leg12a->AddEntry(hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::cosThPol], name, "p");
  //   }
  //   leg12a->SetTextSize(0.035); leg12a->SetFillColor(0);
  //   leg12a->Draw();

  //   if(iRap == 1)
  //     sprintf(name, "|y| < %1.1f", jpsi::rapForPTRange[iRap]);
  //   else
  //     sprintf(name, "%1.1f < |y| < %1.1f", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);

  //   TLatex *tex12a = new TLatex(-0.9, 1.1*hData_pol_pT_rap[iFrame][maxIndex[iRap]][iRap][jpsi::cosThPol]->GetMaximum(), name);
  //   tex12a->SetTextSize(0.06); tex12a->Draw();
    
  //   c12All[iRap]->cd(2);
  //   // gPad->SetLogy(); 
  //   gPad->SetGridy();
  //   TH1F *hFrame12b = gPad->DrawFrame(0., 1., 360., 1.3*hData_pol_pT_rap[iFrame][maxIndex[iRap]][iRap][jpsi::phiPol]->GetMaximum());
  //   sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);  hFrame12b->SetXTitle(name);
  //   sprintf(name, "dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);  hFrame12b->SetYTitle(name);
  //   for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
  //     hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::phiPol]->Draw("psame");

  //   c12All[iRap]->cd(3);
  //   iFrame = 1;
  //   TH1F *hFrame12c = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT_rap[iFrame][maxIndex[iRap]][iRap][jpsi::cosThPol]->GetMaximum());
  //   sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);  hFrame12c->SetXTitle(name);
  //   sprintf(name, "dN/d(cos#theta_{%s})", jpsi::frameLabel[iFrame]);  hFrame12c->SetYTitle(name);
  //   for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)
  //     hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::cosThPol]->Draw("psame");

  //   c12All[iRap]->cd(4);
  //   // gPad->SetLogy(); 
  //   gPad->SetGridy();
  //   TH1F *hFrame12d = gPad->DrawFrame(0., 1., 360., 1.3*hData_pol_pT_rap[iFrame][maxIndex[iRap]][iRap][jpsi::phiPol]->GetMaximum());
  //   sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);  hFrame12d->SetXTitle(name);
  //   sprintf(name, "dN/d#phi_{%s}", jpsi::frameLabel[iFrame]);  hFrame12d->SetYTitle(name);
  //   for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  //     hData_pol_pT_rap[iFrame][iPTBin][iRap][jpsi::phiPol]->Draw("psame");
  //   }

  //   sprintf(name, "Figures/dataCorr_%s_rap%d_pTBins.eps", label, iRap);  c12All[iRap]->Print(name);
  //   sprintf(name, "Figures/dataCorr_%s_rap%d_pTBins.pdf", label, iRap);  c12All[iRap]->Print(name);
  //   sprintf(name, "Figures/dataCorr_%s_rap%d_pTBins.gif", label, iRap);  c12All[iRap]->Print(name);
  // }
}

//======================================
void PlotCorrData2D(Int_t iFrame, Char_t *label, Char_t *oniaLabel){

  gStyle->SetOptStat(0);
  Char_t name[100], title[100];
  TCanvas *c2D[jpsi::kNbRapForPTBins+1];
  for(int iRap = 1; iRap <= jpsi::kNbRapForPTBins; iRap++){
    sprintf(name, "c2DAfter_%s_rap%d", jpsi::frameLabel[iFrame], iRap);
    sprintf(title, "Acc corrected phi vs cosTheta for pT bins, rap = %d (%s)", iRap, jpsi::frameLabel[iFrame]);
    c2D[iRap] = new TCanvas(name, title, 1000, 700);
    c2D[iRap]->Divide(3,2);
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      c2D[iRap]->cd(iPTBin);
      if(iPTBin == jpsi::kNbPTBins) 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", oniaLabel, jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPTBin-1]);
      else 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", oniaLabel, jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
      hData2D_pol_pT_rap[iFrame][iPTBin][iRap]->SetTitle(name);
//       hData2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("colz");
      hData2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("lego2");
    }
    sprintf(name, "Figures/dataCorr2D_%s_%s_rap%d_pTBins.eps", jpsi::frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
    sprintf(name, "Figures/dataCorr2D_%s_%s_rap%d_pTBins.pdf", jpsi::frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
    sprintf(name, "Figures/dataCorr2D_%s_%s_rap%d_pTBins.png", jpsi::frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
  }
}

//======================================
void PlotCorr2D_OneByOne(Int_t iFrame, Char_t *label, Char_t *oniaLabel, Int_t rapBin, Int_t pTBin){

  gStyle->SetOptStat(0);
  Char_t name[100], title[100];
  TCanvas *c2D[jpsi::kNbRapForPTBins+1];
  for(int iRap = rapBin; iRap < rapBin+1; iRap++){
    sprintf(name, "c2DAfter2_%s_rap%d", jpsi::frameLabel[iFrame], iRap);
    sprintf(title, "Acc corrected phi vs cosTheta for pT bin %d, rap = %d (%s)", pTBin, iRap, jpsi::frameLabel[iFrame]);
    c2D[iRap] = new TCanvas(name, title, 1000, 700);
    for(int iPTBin = pTBin; iPTBin < pTBin+1; iPTBin++){
      c2D[iRap]->cd(iPTBin);
      if(iPTBin == jpsi::kNbPTBins) 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", oniaLabel, jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPTBin-1]);
      else 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", oniaLabel, jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
      hData2D_pol_pT_rap[iFrame][iPTBin][iRap]->SetTitle(name);
//       hData2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("colz");
      hData2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("lego2");
    }
    sprintf(name, "Figures/dataCorr2D_%s_%s_rap%d_pT%d.eps", jpsi::frameLabel[iFrame], label, iRap, pTBin);  c2D[iRap]->Print(name);
    sprintf(name, "Figures/dataCorr2D_%s_%s_rap%d_pT%d.pdf", jpsi::frameLabel[iFrame], label, iRap, pTBin);  c2D[iRap]->Print(name);
    sprintf(name, "Figures/dataCorr2D_%s_%s_rap%d_pT%d.gif", jpsi::frameLabel[iFrame], label, iRap, pTBin);  c2D[iRap]->Print(name);
  }
}

//======================================
void CorrectForAcc(Char_t *label){

  Char_t name[100];
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){

      sprintf(name, "hData_Onia_cosTh_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
      hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol] = (TH1D *) Reco_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->Clone(name);
      // hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->Divide(hAcc_pol_pT[iFrame][iPTBin][jpsi::cosThPol]);
      printf("frame %d, pT %d, acceptance correcting 1D now\n", iFrame, iPTBin);
      AccCorrect(hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol], gAcc_pol_pT[iFrame][iPTBin][jpsi::cosThPol]);
      hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
      hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
      hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);

      sprintf(name, "hData_Onia_phi_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
      hData_pol_pT[iFrame][iPTBin][jpsi::phiPol] = (TH1D *) Reco_pol_pT[iFrame][iPTBin][jpsi::phiPol]->Clone(name);
      // hData_pol_pT[iFrame][iPTBin][jpsi::phiPol]->Divide(hAcc_pol_pT[iFrame][iPTBin][jpsi::phiPol]);
      printf("frame %d, pT %d, acceptance correcting 1D now\n", iFrame, iPTBin);
      AccCorrect(hData_pol_pT[iFrame][iPTBin][jpsi::phiPol], gAcc_pol_pT[iFrame][iPTBin][jpsi::phiPol]);
      hData_pol_pT[iFrame][iPTBin][jpsi::phiPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
      hData_pol_pT[iFrame][iPTBin][jpsi::phiPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
      hData_pol_pT[iFrame][iPTBin][jpsi::phiPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
      
      sprintf(name, "hData2D_Onia_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
      hData2D_pol_pT[iFrame][iPTBin] = (TH2D *) Reco2D_pol_pT[iFrame][iPTBin]->Clone(name);
      hData2D_pol_pT[iFrame][iPTBin]->Divide(hAcc2D_pol_pT[iFrame][iPTBin]);
      hData2D_pol_pT[iFrame][iPTBin]->SetLineColor(jpsi::colour_pT[iPTBin]);
      hData2D_pol_pT[iFrame][iPTBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
      hData2D_pol_pT[iFrame][iPTBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    }
    for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
      sprintf(name, "hData_Onia_cosTh_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
      hData_pol_rap[iFrame][iRapBin][jpsi::cosThPol] = (TH1D *) Reco_pol_rap[iFrame][iRapBin][jpsi::cosThPol]->Clone(name);
      // hData_pol_rap[iFrame][iRapBin][jpsi::cosThPol]->Divide(hAcc_pol_rap[iFrame][iRapBin][jpsi::cosThPol]);
      printf("frame %d, rap %d, acceptance correcting 1D now\n", iFrame, iRapBin);
      AccCorrect(hData_pol_rap[iFrame][iRapBin][jpsi::cosThPol], gAcc_pol_pT[iFrame][iRapBin][jpsi::cosThPol]);
      hData_pol_rap[iFrame][iRapBin][jpsi::cosThPol]->SetLineColor(jpsi::colour_rap[iRapBin]);
      hData_pol_rap[iFrame][iRapBin][jpsi::cosThPol]->SetMarkerColor(jpsi::colour_rap[iRapBin]);
      hData_pol_rap[iFrame][iRapBin][jpsi::cosThPol]->SetMarkerStyle(jpsi::marker_rap[iRapBin]);
      
      sprintf(name, "hData_Onia_phi_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
      hData_pol_rap[iFrame][iRapBin][jpsi::phiPol] = (TH1D *) Reco_pol_rap[iFrame][iRapBin][jpsi::phiPol]->Clone(name);
      // hData_pol_rap[iFrame][iRapBin][jpsi::phiPol]->Divide(hAcc_pol_rap[iFrame][iRapBin][jpsi::phiPol]);
      printf("frame %d, rap %d, acceptance correcting 1D now\n", iFrame, iRapBin);
      AccCorrect(hData_pol_rap[iFrame][iRapBin][jpsi::phiPol], gAcc_pol_pT[iFrame][iRapBin][jpsi::phiPol]);
      hData_pol_rap[iFrame][iRapBin][jpsi::phiPol]->SetLineColor(jpsi::colour_rap[iRapBin]);
      hData_pol_rap[iFrame][iRapBin][jpsi::phiPol]->SetMarkerColor(jpsi::colour_rap[iRapBin]);
      hData_pol_rap[iFrame][iRapBin][jpsi::phiPol]->SetMarkerStyle(jpsi::marker_rap[iRapBin]);
      
      sprintf(name, "hData2D_Onia_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
      hData2D_pol_rap[iFrame][iRapBin] = (TH2D *) Reco2D_pol_rap[iFrame][iRapBin]->Clone(name);
      hData2D_pol_rap[iFrame][iRapBin]->Divide(hAcc2D_pol_rap[iFrame][iRapBin]);
      hData2D_pol_rap[iFrame][iRapBin]->SetLineColor(jpsi::colour_rap[iRapBin]);
      hData2D_pol_rap[iFrame][iRapBin]->SetMarkerColor(jpsi::colour_rap[iRapBin]);
      hData2D_pol_rap[iFrame][iRapBin]->SetMarkerStyle(jpsi::marker_rap[iRapBin]);
    }   
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
	sprintf(name, "hData_Onia_cosTh_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol] = (TH1D *) Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->Clone(name);
	//hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->Divide(hAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]);
	printf("frame %d, pT %d, rap %d acceptance correcting 1D now\n", iFrame, iPTBin, iRapBin);
	AccCorrect(hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol], gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
	
	sprintf(name, "hData_Onia_phi_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol] = (TH1D *) Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->Clone(name);
	//hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->Divide(hAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]);
	printf("frame %d, pT %d, rap %d acceptance correcting 1D now\n", iFrame, iPTBin, iRapBin);
	AccCorrect(hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol], gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);

	sprintf(name, "hData2D_Onia_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Clone(name);
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Divide(hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]);//H: can be done in the way below (commented)
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetLineColor(jpsi::colour_pT[iPTBin]);
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
      }
    }
  }
  
  // Double_t binContentReco, binErrorReco;
  // Double_t acc, binContentData, binErrorData;
  // Double_t content, error, maxAcc;

  // //acceptance correct the 2D histogram
  // for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
  //   for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  //     for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
  // 	for(int iBinX = 1; iBinX <= hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(); iBinX++){
  // 	  for(int iBinY = 1; iBinY <= hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY(); iBinY++){
	    
  // 	    binContentReco = Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iBinX, iBinY);
  // 	    binErrorReco = Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinError(iBinX, iBinY);
  // 	    acc = hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iBinX, iBinY);

  // 	    if(acc > 0.){
  // 	      binErrorData = binErrorReco / acc;
  // 	      binContentData = binContentReco / acc;
  // 	    }
  // 	    else{
  // 	      binErrorData = 1.;
  // 	      binContentData = 0.;
  // 	    }
  // 	    hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(iBinX, iBinY, binContentData);
  // 	    hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinError(iBinX, iBinY, binErrorData);
  // 	  }
  // 	}
  //     }
  //   }
  // }
  Get1DHistoFrom2D();
  SaveCorrData(label); //save the histos here BEFORE we set bin contents to zero for further processing!

}

//======================================
void ReadData(Char_t *fileNameData, Int_t rebinCosTh, Int_t rebinPhi){

  TFile *fIn = new TFile(fileNameData);
  
  Char_t name[100];
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      sprintf(name, "Reco_Onia_cosTh_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
      Reco_pol_pT[iFrame][iPTBin][jpsi::cosThPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
      sprintf(name, "Reco_Onia_phi_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
      Reco_pol_pT[iFrame][iPTBin][jpsi::phiPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_pT[iFrame][iPTBin][jpsi::phiPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][jpsi::phiPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][jpsi::phiPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
      sprintf(name, "Reco2D_Onia_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
      Reco2D_pol_pT[iFrame][iPTBin] = (TH2D *) gDirectory->Get(name);
      Reco2D_pol_pT[iFrame][iPTBin]->SetLineColor(jpsi::colour_pT[iPTBin]);
      Reco2D_pol_pT[iFrame][iPTBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
      Reco2D_pol_pT[iFrame][iPTBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
      if(rebinCosTh > 1 || rebinPhi > 1)
	Reco2D_pol_pT[iFrame][iPTBin]->Rebin2D(rebinCosTh, rebinPhi);
    }
    for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
      sprintf(name, "Reco_Onia_cosTh_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
      Reco_pol_rap[iFrame][iRapBin][jpsi::cosThPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_rap[iFrame][iRapBin][jpsi::cosThPol]->SetLineColor(jpsi::colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][jpsi::cosThPol]->SetMarkerColor(jpsi::colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][jpsi::cosThPol]->SetMarkerStyle(jpsi::marker_rap[iRapBin]);
      sprintf(name, "Reco_Onia_phi_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
      Reco_pol_rap[iFrame][iRapBin][jpsi::phiPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_rap[iFrame][iRapBin][jpsi::phiPol]->SetLineColor(jpsi::colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][jpsi::phiPol]->SetMarkerColor(jpsi::colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][jpsi::phiPol]->SetMarkerStyle(jpsi::marker_rap[iRapBin]);
      sprintf(name, "Reco2D_Onia_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
      Reco2D_pol_rap[iFrame][iRapBin] = (TH2D *) gDirectory->Get(name);
      Reco2D_pol_rap[iFrame][iRapBin]->SetLineColor(jpsi::colour_rap[iRapBin]);
      Reco2D_pol_rap[iFrame][iRapBin]->SetMarkerColor(jpsi::colour_rap[iRapBin]);
      Reco2D_pol_rap[iFrame][iRapBin]->SetMarkerStyle(jpsi::marker_rap[iRapBin]);
      if(rebinCosTh > 1 || rebinPhi > 1)
	Reco2D_pol_rap[iFrame][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
    }
    Int_t totEv = 0;
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
	sprintf(name, "Reco_Onia_cosTh_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol] = (TH1D *) gDirectory->Get(name);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
	sprintf(name, "Reco_Onia_phi_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol] = (TH1D *) gDirectory->Get(name);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) gDirectory->Get(name);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetLineColor(jpsi::colour_pT[iPTBin]);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
	printf("%s: pTbin %d, rapBin %d --> statistics in data: %1.1f\n",
	       jpsi::frameLabel[iFrame], iPTBin, iRapBin, Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Integral());

	totEv += Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Integral();
	if(rebinCosTh > 1 || rebinPhi > 1)
	  Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
      }
    }
    printf("%s: total nb. of events: %d\n\n\n", jpsi::frameLabel[iFrame], totEv);
    
  }
}

//===================================
void ReadAccHistos(Char_t *fileNameMC, Int_t rebinCosTh, Int_t rebinPhi, Bool_t normalising){

  printf("reading in the acceptance histograms\n");
  TFile *fInMC = new TFile(fileNameMC);

  Char_t name[100];
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      sprintf(name, "gAcc_Onia_cosTh_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
      gAcc_pol_pT[iFrame][iPTBin][jpsi::cosThPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
      // printf("frame %d, pT %d, histo %p\n", iFrame, iPTBin, gAcc_pol_pT[iFrame][iPTBin][jpsi::cosThPol]);
      //
      sprintf(name, "gAcc_Onia_phi_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
      gAcc_pol_pT[iFrame][iPTBin][jpsi::phiPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
      // printf("frame %d, pT %d, histo %p\n", iFrame, iPTBin, gAcc_pol_pT[iFrame][iPTBin][jpsi::phiPol]);
      //
      sprintf(name, "hAcc2D_Onia_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
      hAcc2D_pol_pT[iFrame][iPTBin] = (TH2D *) gDirectory->Get(name);
      // printf("frame %d, pT %d, histo %p\n", iFrame, iPTBin, hAcc2D_pol_pT[iFrame][iPTBin]);
      if(rebinCosTh > 1 || rebinPhi > 1)
	hAcc2D_pol_pT[iFrame][iPTBin]->Rebin2D(rebinCosTh, rebinPhi);
    }
    for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
      sprintf(name, "gAcc_Onia_cosTh_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
      gAcc_pol_rap[iFrame][iRapBin][jpsi::cosThPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
      // printf("frame %d, rap %d, histo %p\n", iFrame, iRapBin, gAcc_pol_rap[iFrame][iRapBin][jpsi::cosThPol]);
      //
      sprintf(name, "gAcc_Onia_phi_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
      gAcc_pol_rap[iFrame][iRapBin][jpsi::phiPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
      // printf("frame %d, rap %d, histo %p\n", iFrame, iRapBin, gAcc_pol_rap[iFrame][iRapBin][jpsi::phiPol]);
      //
      sprintf(name, "hAcc2D_Onia_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
      hAcc2D_pol_rap[iFrame][iRapBin] = (TH2D *) gDirectory->Get(name);
      // printf("frame %d, rap %d, histo %p\n", iFrame, iRapBin, hAcc2D_pol_rap[iFrame][iRapBin]);
      if(rebinCosTh > 1 || rebinPhi > 1)
	hAcc2D_pol_rap[iFrame][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
    }
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
	sprintf(name, "gAcc_Onia_cosTh_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
	// printf("frame %d, pT %d, rap %d, histo %p\n", iFrame, iPTBin, iRapBin, gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]);
	//
	sprintf(name, "gAcc_Onia_phi_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
	// printf("frame %d, pT %d, rap %d, histo %p\n", iFrame, iPTBin, iRapBin, gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]);
	//
	sprintf(name, "hAcc2D_Onia_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) gDirectory->Get(name);
	// printf("frame %d, pT %d, rap %d, histo %p\n", iFrame, iPTBin, iRapBin, hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]);
	if(rebinCosTh > 1 || rebinPhi > 1)
	  hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
	if(normalising && hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Integral() > 0)
	  hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Scale(1./hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Integral());
      }
    }
  }
}

//============================================
void SaveCorrData(Char_t *label){

  Char_t name[100];
  sprintf(name, "polarization_corrData_%s.root", label);
  TFile *fOut = new TFile(name, "RECREATE");
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      hData_pol_pT[iFrame][iPTBin][jpsi::cosThPol]->Write();
      hData_pol_pT[iFrame][iPTBin][jpsi::phiPol]->Write();
      hData2D_pol_pT[iFrame][iPTBin]->Write();
    }
    for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
      hData_pol_rap[iFrame][iRapBin][jpsi::cosThPol]->Write();
      hData_pol_rap[iFrame][iRapBin][jpsi::phiPol]->Write();
      hData2D_pol_rap[iFrame][iRapBin]->Write();
    }
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->Write();
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->Write();
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }
}

//============================================
void AccCorrect(TH1D *hist, TGraphAsymmErrors *gAcc){

  // Int_t nBinsHist = hist->GetNbinsX();
  // Int_t nBinsGraph = gAcc->GetN();
  // Double_t x, acc;
  // printf("#bins(hist) = %d, #bins(graph) = %d\n", nBinsHist, nBinsGraph);
  // if(nBinsHist != nBinsGraph){
  //   printf("cannot acceptance correct because the TGraph and the histo have differnt number of entries: %d vs %d\n", 
  //   	   nBinsGraph, nBinsHist);
  //   // for(int iP = 0; iP < nBinsGraph; iP++){
  //   //   gAcc->GetPoint(iP, x, acc);
  //   //   printf("bin %d sits at x = %f (TGraph) and %f (histo) \n", iP, x, hist->GetBinCenter(iP+1));
  //   // }
  //   // printf("last bin of histo at x = %f\n", hist->GetBinCenter( hist->GetNbinsX()));
  //   return;
  // }

  Double_t x, acc;
  Double_t content, binErr;
  for(int iBin = 0; iBin < gAcc->GetN(); iBin++){

    gAcc->GetPoint(iBin, x, acc);
    Int_t binHist = hist->GetXaxis()->FindBin(x);
    content = hist->GetBinContent(binHist);
    binErr = hist->GetBinError(binHist);

    if(acc > 0.){
      hist->SetBinContent(binHist, content / acc);
      hist->SetBinError(binHist, binErr / acc);
    }
  }
}
