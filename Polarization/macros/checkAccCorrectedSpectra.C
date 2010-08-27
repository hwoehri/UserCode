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
TH1D *Reco_pol_pT[kNbFrames][kNbPTBins+1][kNbPolVar];
TH1D *Reco_pol_rap[kNbFrames][2*kNbRapBins+1][kNbPolVar];
TH1D *Reco_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1][kNbPolVar];
TH2D *Reco2D_pol_pT[kNbFrames][kNbPTBins+1];
TH2D *Reco2D_pol_rap[kNbFrames][2*kNbRapBins+1];
TH2D *Reco2D_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1];

//acceptance histos
TGraphAsymmErrors *gAcc_pol_pT[kNbFrames][kNbPTBins+1][kNbPolVar];
TGraphAsymmErrors *gAcc_pol_rap[kNbFrames][2*kNbRapBins+1][kNbPolVar];
TGraphAsymmErrors *gAcc_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1][kNbPolVar];
TH2D *hAcc2D_pol_pT[kNbFrames][kNbPTBins+1];
TH2D *hAcc2D_pol_rap[kNbFrames][2*kNbRapBins+1];
TH2D *hAcc2D_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1];
TGraphAsymmErrors *gAcc2D_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1];

//data corrected for acc:
TH1D *hData_pol_pT[kNbFrames][kNbPTBins+1][kNbPolVar];
TH1D *hData_pol_rap[kNbFrames][2*kNbRapBins+1][kNbPolVar];
TH1D *hData_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1][kNbPolVar];
TH2D *hData2D_pol_pT[kNbFrames][kNbPTBins+1];
TH2D *hData2D_pol_rap[kNbFrames][2*kNbRapBins+1];
TH2D *hData2D_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1];

//mega histo containing the N phi bins appended:
TH1D *hData1D_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1];

TH2D *thisHist;
TH2D *thisReco, *thisAcc;
TH1D *thisHist1D;

void ReadData(Char_t *fileNameData, Int_t rebinCosTh, Int_t rebinPhi);
void PlotUncorrData(Int_t iFrame, Char_t *polTag);
void PlotUncorrData2D(Int_t iFrame, Char_t *polTag, Char_t *oniaLabel);
void PlotCorrData2D(Int_t iFrame, Char_t *polTag, Char_t *oniaLabel);
void ReadAccHistos(Char_t *fileNameMC);
void CorrectForAcc(Char_t *polTag);
void Get1DHistoFrom2D();
void PlotHistos(Int_t iFrame, Char_t *polTag);
void PlotAll(Char_t *polTag);
void SaveCorrData(Char_t *polTag);
void AccCorrect(TH1D *hist, TGraphAsymmErrors *gAcc);
//======================================
//usage: root checkAccCorrectedSpectra.C+ or
//       root 'checkAccCorrectedSpectra.C+("HLT_Mu0Track0Jpsi", 1, 1)' (e.g.)
//======================================
void checkAccCorrectedSpectra(Char_t *hltTag = "HLT_Mu0Track0Jpsi",
			      Int_t rebinCosTh = 2, //histos will be rebinned by "rebinCosTh"
			      Int_t rebinPhi = 2, //histos will be rebinned by "rebinPhi"
			      Char_t *oniaLabel = "J/#psi"){

  Char_t fileNameMC[200];
  sprintf(fileNameMC, "/home/hermine/CMS/Work/Polarization/Florian/23Aug2010/files/accHistos_%s.root", hltTag);

  Char_t label[100];
  sprintf(label, "%s", hltTag);

  Char_t fileNameData[200];
  sprintf(fileNameData, "pol_data_HLT_Mu0Track0Jpsi.root", label);

  ReadData(fileNameData, rebinCosTh, rebinPhi);
  ReadAccHistos(fileNameMC);
  CorrectForAcc(label); //calls internally "Get1DHistoFrom2D" to fill
                         //the histo hData1D_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1];

  //1D histos (cosTheta and phi):
  PlotUncorrData(CS, label);
  PlotUncorrData(HX, label);
  PlotHistos(CS, label);
  PlotHistos(HX, label);
  PlotAll(label);

  //2D histos (cosTheta and phi):
  PlotUncorrData2D(CS, label, oniaLabel);
  PlotUncorrData2D(HX, label, oniaLabel);
  PlotCorrData2D(CS, label, oniaLabel);
  PlotCorrData2D(HX, label, oniaLabel);

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
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
	//
	sprintf(name, "hData1D_Onia_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH1D(name, "", nBins1D, 0., (Double_t) nBins1D);
	binID = 0;
// 	printf("\n\n\n%s, pT %d, rap %d\n", frameLabel[iFrame], iPTBin, iRapBin);

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
  sprintf(name, "c20CosTh_%s", frameLabel[iFrame]);
  TCanvas *c20CosTh = new TCanvas(name, "cosTheta for pT bins", 900, 700);
  TH1F *hFrame20 = gPad->DrawFrame(-1., 0., 1., 1.3*Reco_pol_pT[iFrame][0][cosThPol]->GetMaximum());
  sprintf(name, "cos#theta_{%s}", frameLabel[iFrame]);
  hFrame20->SetXTitle(name);
  sprintf(name, "Acc #times dN/d(cos#theta_{%s})", frameLabel[iFrame]);
  hFrame20->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
//   for(int iPTBin = 0; iPTBin < 3; iPTBin++){
    Reco_pol_pT[iFrame][iPTBin][cosThPol]->Draw("psame");
  }
  TLegend *leg20a = new TLegend(0.65,0.735119,0.9866071,0.985119);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
    if(iPTBin == 0) sprintf(name, "all p_{T}");
    else if(iPTBin == kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", pTRange[iPTBin-1]);
    else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", pTRange[iPTBin-1], pTRange[iPTBin]);
    leg20a->AddEntry(Reco_pol_pT[iFrame][iPTBin][cosThPol], name, "p");
  }
  leg20a->SetTextSize(0.035); leg20a->SetFillColor(0);
  leg20a->Draw();
 
  sprintf(name, "Figures/reco_cosTheta_%s_%s_pTBins.eps", frameLabel[iFrame], label);  c20CosTh->Print(name);
  sprintf(name, "Figures/reco_cosTheta_%s_%s_pTBins.pdf", frameLabel[iFrame], label);  c20CosTh->Print(name);
  sprintf(name, "Figures/reco_cosTheta_%s_%s_pTBins.gif", frameLabel[iFrame], label);  c20CosTh->Print(name);

  //===========================================
  //phi for different pT bins
  //===========================================
  sprintf(name, "c20Phi_%s", frameLabel[iFrame]);
  TCanvas *c20Phi = new TCanvas(name, "phi for pT bins", 900, 700);
  TH1F *hFrame20b = gPad->DrawFrame(0., 0., 360., 1.3*Reco_pol_pT[iFrame][0][phiPol]->GetMaximum());
  sprintf(name, "#phi_{%s} [deg]", frameLabel[iFrame]);
  hFrame20b->SetXTitle(name);
  sprintf(name, "Acc #times dN/d#phi_{%s}", frameLabel[iFrame]);
  hFrame20b->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
//   for(int iPTBin = 0; iPTBin < 3; iPTBin++){
    Reco_pol_pT[iFrame][iPTBin][phiPol]->Draw("psame");
  }
//   TLegend *leg20b = new TLegend(0.1417411,0.1339286,0.5178571,0.3839286);
//   for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
//     if(iPTBin == 0) sprintf(name, "all p_{T}");
//     else if(iPTBin == kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", pTRange[iPTBin-1]);
//     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", pTRange[iPTBin-1], pTRange[iPTBin]);
//     leg20b->AddEntry(Reco_pol_pT[iFrame][iPTBin][phiPol], name, "p");
//   }
//   leg20b->SetTextSize(0.035); leg20b->SetFillColor(0);
//   leg20b->Draw();

  sprintf(name, "Figures/reco_phi_%s_%s_pTBins.eps", frameLabel[iFrame], label);  c20Phi->Print(name);
  sprintf(name, "Figures/reco_phi_%s_%s_pTBins.pdf", frameLabel[iFrame], label);  c20Phi->Print(name);
  sprintf(name, "Figures/reco_phi_%s_%s_pTBins.gif", frameLabel[iFrame], label);  c20Phi->Print(name);

}

//=================================
void PlotUncorrData2D(Int_t iFrame, Char_t *label, Char_t *oniaLabel){

  Char_t name[100], title[100];
  TCanvas *c2D[kNbRapForPTBins+1];
  for(int iRap = 1; iRap <= kNbRapForPTBins; iRap++){
    sprintf(name, "c2D_%s_rap%d", frameLabel[iFrame], iRap);
    sprintf(title, "phi vs cosTheta for pT bins, rap = %d (%s)", iRap, frameLabel[iFrame]);
    c2D[iRap] = new TCanvas(name, title, 1000, 700);
    c2D[iRap]->Divide(3,2);
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      c2D[iRap]->cd(iPTBin);
      if(iPTBin == kNbPTBins) 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", oniaLabel, rapForPTRange[iRap-1], rapForPTRange[iRap], pTRange[iPTBin-1]);
      else 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", oniaLabel, rapForPTRange[iRap-1], rapForPTRange[iRap], pTRange[iPTBin-1], pTRange[iPTBin]);
      Reco2D_pol_pT_rap[iFrame][iPTBin][iRap]->SetTitle(name);
//       Reco2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("colz");
      Reco2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("lego2");
    }
    sprintf(name, "Figures/reco2D_%s_%s_rap%d_pTBins.eps", frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
    sprintf(name, "Figures/reco2D_%s_%s_rap%d_pTBins.pdf", frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
    sprintf(name, "Figures/reco2D_%s_%s_rap%d_pTBins.png", frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
  }
}

//======================================
void PlotHistos(Int_t iFrame, Char_t *label){

  Char_t name[100];
  //===========================================
  //cosTheta for different pT bins
  //===========================================
  sprintf(name, "c10CosTh_%s", frameLabel[iFrame]);
  TCanvas *c10CosTh = new TCanvas(name, "cosTheta for pT bins", 900, 700);
  TH1F *hFrame10 = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT[iFrame][0][cosThPol]->GetMaximum());
  sprintf(name, "cos#theta_{%s}", frameLabel[iFrame]);
  hFrame10->SetXTitle(name);
  sprintf(name, "dN/d(cos#theta_{%s})", frameLabel[iFrame]);
  hFrame10->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
//   for(int iPTBin = 0; iPTBin < 3; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][cosThPol]->Draw("psame");
  }
  TLegend *leg10a = new TLegend(0.6104911,0.735119,0.9866071,0.985119);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
    if(iPTBin == 0) sprintf(name, "all p_{T}");
    else if(iPTBin == kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", pTRange[iPTBin-1]);
    else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", pTRange[iPTBin-1], pTRange[iPTBin]);
    leg10a->AddEntry(hData_pol_pT[iFrame][iPTBin][cosThPol], name, "p");
  }
  leg10a->SetTextSize(0.035); leg10a->SetFillColor(0);
  leg10a->Draw();

 
  sprintf(name, "Figures/dataCorr_cosTheta_%s_%s_pTBins.eps", frameLabel[iFrame], label);  c10CosTh->Print(name);
  sprintf(name, "Figures/dataCorr_cosTheta_%s_%s_pTBins.pdf", frameLabel[iFrame], label);  c10CosTh->Print(name);
  sprintf(name, "Figures/dataCorr_cosTheta_%s_%s_pTBins.gif", frameLabel[iFrame], label);  c10CosTh->Print(name);

  //===========================================
  //phi for different pT bins
  //===========================================
  sprintf(name, "c10Phi_%s", frameLabel[iFrame]);
  TCanvas *c10Phi = new TCanvas(name, "phi for pT bins", 900, 700);
  TH1F *hFrame10b = gPad->DrawFrame(0., 0., 360., 1.3*hData_pol_pT[iFrame][0][phiPol]->GetMaximum());
  sprintf(name, "#phi_{%s} [deg]", frameLabel[iFrame]);
  hFrame10b->SetXTitle(name);
  sprintf(name, "dN/d#phi_{%s}", frameLabel[iFrame]);
  hFrame10b->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
//   for(int iPTBin = 0; iPTBin < 3; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][phiPol]->Draw("psame");
  }
//   TLegend *leg10b = new TLegend(0.1417411,0.1339286,0.5178571,0.3839286);
//   for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
//     if(iPTBin == 0) sprintf(name, "all p_{T}");
//     else if(iPTBin == kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", pTRange[iPTBin-1]);
//     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", pTRange[iPTBin-1], pTRange[iPTBin]);
//     leg10b->AddEntry(hData_pol_pT[iFrame][iPTBin][phiPol], name, "p");
//   }
//   leg10b->SetTextSize(0.035); leg10b->SetFillColor(0);
//   leg10b->Draw();

  sprintf(name, "Figures/dataCorr_phi_%s_%s_pTBins.eps", frameLabel[iFrame], label);  c10Phi->Print(name);
  sprintf(name, "Figures/dataCorr_phi_%s_%s_pTBins.pdf", frameLabel[iFrame], label);  c10Phi->Print(name);
  sprintf(name, "Figures/dataCorr_phi_%s_%s_pTBins.gif", frameLabel[iFrame], label);  c10Phi->Print(name);
}

//======================================
void PlotAll(Char_t *label){

  Char_t name[100];
  //=====================================
  //rapidity integrated spectra
  //=====================================
  sprintf(name, "c10All_%s", label);
  TCanvas *c10All = new TCanvas(name, "corrected spectra for pT bins", 1000, 700);
  c10All->Divide(2,2);
  c10All->cd(1);
  Int_t iFrame = 0;
  TH1F *hFrame10 = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT[iFrame][0][cosThPol]->GetMaximum());
  sprintf(name, "cos#theta_{%s}", frameLabel[iFrame]);  hFrame10->SetXTitle(name);
  sprintf(name, "dN/d(cos#theta_{%s})", frameLabel[iFrame]);  hFrame10->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][cosThPol]->Draw("psame");
  }
  c10All->cd(2);
  TH1F *hFrame10b = gPad->DrawFrame(0., 0., 360., 1.3*hData_pol_pT[iFrame][0][phiPol]->GetMaximum());
  sprintf(name, "#phi_{%s} [deg]", frameLabel[iFrame]);  hFrame10b->SetXTitle(name);
  sprintf(name, "dN/d#phi_{%s}", frameLabel[iFrame]);  hFrame10b->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][phiPol]->Draw("psame");
  }
  TLegend *leg10a = new TLegend(0.6171352,0.6922123,0.9936412,0.9929315);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
    if(iPTBin == 0) sprintf(name, "all p_{T}");
    else if(iPTBin == kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", pTRange[iPTBin-1]);
    else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", pTRange[iPTBin-1], pTRange[iPTBin]);
    leg10a->AddEntry(hData_pol_pT[iFrame][iPTBin][cosThPol], name, "p");
  }
  leg10a->SetTextSize(0.035); leg10a->SetFillColor(0);
  leg10a->Draw();

  c10All->cd(3);
  iFrame = 1;
  TH1F *hFrame10c = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT[iFrame][0][cosThPol]->GetMaximum());
  sprintf(name, "cos#theta_{%s}", frameLabel[iFrame]);  hFrame10c->SetXTitle(name);
  sprintf(name, "dN/d(cos#theta_{%s})", frameLabel[iFrame]);  hFrame10c->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][cosThPol]->Draw("psame");
  }
  c10All->cd(4);
  TH1F *hFrame10d = gPad->DrawFrame(0., 0., 360., 1.3*hData_pol_pT[iFrame][0][phiPol]->GetMaximum());
  sprintf(name, "#phi_{%s} [deg]", frameLabel[iFrame]);  hFrame10d->SetXTitle(name);
  sprintf(name, "dN/d#phi_{%s}", frameLabel[iFrame]);  hFrame10d->SetYTitle(name);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
    hData_pol_pT[iFrame][iPTBin][phiPol]->Draw("psame");
  }

  sprintf(name, "Figures/dataCorr_%s_pTBins.eps", label);  c10All->Print(name);
  sprintf(name, "Figures/dataCorr_%s_pTBins.pdf", label);  c10All->Print(name);
  sprintf(name, "Figures/dataCorr_%s_pTBins.gif", label);  c10All->Print(name);

  //=====================================
  //rapidity differential spectra
  //=====================================
  Char_t title[100];
  TCanvas *c11All[kNbRapForPTBins+1];
  Int_t maxIndex[kNbRapForPTBins+1] = {1, 4, 3, 2};
  for(int iRap = 1; iRap < kNbRapForPTBins+1; iRap++){
    sprintf(name, "c11All_%s_rap%d", label, iRap);
    sprintf(title, "corrected spectra for pT bins (rap %d)", iRap);
    c11All[iRap] = new TCanvas(name, title, 1000, 700);
    c11All[iRap]->Divide(2,2);
    c11All[iRap]->cd(1);
    Int_t iFrame = 0;
    TH1F *hFrame10 = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT_rap[iFrame][maxIndex[iRap]][iRap][cosThPol]->GetMaximum());
    sprintf(name, "cos#theta_{%s}", frameLabel[iFrame]);  hFrame10->SetXTitle(name);
    sprintf(name, "dN/d(cos#theta_{%s})", frameLabel[iFrame]);  hFrame10->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++)
      hData_pol_pT_rap[iFrame][iPTBin][iRap][cosThPol]->Draw("psame");

    TLegend *leg10a = new TLegend(0.7,0.6922123,0.9936412,0.9929315);
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      if(iPTBin == 0) sprintf(name, "all p_{T}");
      else if(iPTBin == kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", pTRange[iPTBin-1]);
      else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", pTRange[iPTBin-1], pTRange[iPTBin]);
      leg10a->AddEntry(hData_pol_pT_rap[iFrame][iPTBin][iRap][cosThPol], name, "p");
    }
    leg10a->SetTextSize(0.035); leg10a->SetFillColor(0);
    leg10a->Draw();

    if(iRap == 1)
      sprintf(name, "|y| < %1.1f", rapForPTRange[iRap]);
    else
      sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[iRap-1], rapForPTRange[iRap]);

    TLatex *tex10a = new TLatex(-0.9, 1.1*hData_pol_pT_rap[iFrame][maxIndex[iRap]][iRap][cosThPol]->GetMaximum(), name);
    tex10a->SetTextSize(0.06); tex10a->Draw();
    
    c11All[iRap]->cd(2);
    // gPad->SetLogy(); 
    gPad->SetGridy();
    TH1F *hFrame10b = gPad->DrawFrame(0., 1., 360., 1.3*hData_pol_pT_rap[iFrame][maxIndex[iRap]][iRap][phiPol]->GetMaximum());
    sprintf(name, "#phi_{%s} [deg]", frameLabel[iFrame]);  hFrame10b->SetXTitle(name);
    sprintf(name, "dN/d#phi_{%s}", frameLabel[iFrame]);  hFrame10b->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++)
      hData_pol_pT_rap[iFrame][iPTBin][iRap][phiPol]->Draw("psame");

    c11All[iRap]->cd(3);
    iFrame = 1;
    TH1F *hFrame10c = gPad->DrawFrame(-1., 0., 1., 1.3*hData_pol_pT_rap[iFrame][maxIndex[iRap]][iRap][cosThPol]->GetMaximum());
    sprintf(name, "cos#theta_{%s}", frameLabel[iFrame]);  hFrame10c->SetXTitle(name);
    sprintf(name, "dN/d(cos#theta_{%s})", frameLabel[iFrame]);  hFrame10c->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++)
      hData_pol_pT_rap[iFrame][iPTBin][iRap][cosThPol]->Draw("psame");

    c11All[iRap]->cd(4);
    // gPad->SetLogy(); 
    gPad->SetGridy();
    TH1F *hFrame10d = gPad->DrawFrame(0., 1., 360., 1.3*hData_pol_pT_rap[iFrame][maxIndex[iRap]][iRap][phiPol]->GetMaximum());
    sprintf(name, "#phi_{%s} [deg]", frameLabel[iFrame]);  hFrame10d->SetXTitle(name);
    sprintf(name, "dN/d#phi_{%s}", frameLabel[iFrame]);  hFrame10d->SetYTitle(name);
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      hData_pol_pT_rap[iFrame][iPTBin][iRap][phiPol]->Draw("psame");
    }

    sprintf(name, "Figures/dataCorr_%s_rap%d_pTBins.eps", label, iRap);  c11All[iRap]->Print(name);
    sprintf(name, "Figures/dataCorr_%s_rap%d_pTBins.pdf", label, iRap);  c11All[iRap]->Print(name);
    sprintf(name, "Figures/dataCorr_%s_rap%d_pTBins.gif", label, iRap);  c11All[iRap]->Print(name);
  }
}

//======================================
void PlotCorrData2D(Int_t iFrame, Char_t *label, Char_t *oniaLabel){

  gStyle->SetOptStat(0);
  Char_t name[100], title[100];
  TCanvas *c2D[kNbRapForPTBins+1];
  for(int iRap = 1; iRap <= kNbRapForPTBins; iRap++){
    sprintf(name, "c2DAfter_%s_rap%d", frameLabel[iFrame], iRap);
    sprintf(title, "Acc corrected phi vs cosTheta for pT bins, rap = %d (%s)", iRap, frameLabel[iFrame]);
    c2D[iRap] = new TCanvas(name, title, 1000, 700);
    c2D[iRap]->Divide(3,2);
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      c2D[iRap]->cd(iPTBin);
      if(iPTBin == kNbPTBins) 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", oniaLabel, rapForPTRange[iRap-1], rapForPTRange[iRap], pTRange[iPTBin-1]);
      else 
	sprintf(name, "%s: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", oniaLabel, rapForPTRange[iRap-1], rapForPTRange[iRap], pTRange[iPTBin-1], pTRange[iPTBin]);
      hData2D_pol_pT_rap[iFrame][iPTBin][iRap]->SetTitle(name);
//       hData2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("colz");
      hData2D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("lego2");
    }
    sprintf(name, "Figures/dataCorr2D_%s_%s_rap%d_pTBins.eps", frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
    sprintf(name, "Figures/dataCorr2D_%s_%s_rap%d_pTBins.pdf", frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
    sprintf(name, "Figures/dataCorr2D_%s_%s_rap%d_pTBins.gif", frameLabel[iFrame], label, iRap);  c2D[iRap]->Print(name);
  }
}

//======================================
void CorrectForAcc(Char_t *label){

  Char_t name[100];
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){

      sprintf(name, "hData_Onia_cosTh_%s_pT%d", frameLabel[iFrame], iPTBin);
      hData_pol_pT[iFrame][iPTBin][cosThPol] = (TH1D *) Reco_pol_pT[iFrame][iPTBin][cosThPol]->Clone(name);
      // hData_pol_pT[iFrame][iPTBin][cosThPol]->Divide(hAcc_pol_pT[iFrame][iPTBin][cosThPol]);
      printf("frame %d, pT %d, acceptance correcting 1D now\n", iFrame, iPTBin);
      AccCorrect(hData_pol_pT[iFrame][iPTBin][cosThPol], gAcc_pol_pT[iFrame][iPTBin][cosThPol]);
      hData_pol_pT[iFrame][iPTBin][cosThPol]->SetLineColor(colour_pT[iPTBin]);
      hData_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerColor(colour_pT[iPTBin]);
      hData_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerStyle(marker_pT[iPTBin]);

      sprintf(name, "hData_Onia_phi_%s_pT%d", frameLabel[iFrame], iPTBin);
      hData_pol_pT[iFrame][iPTBin][phiPol] = (TH1D *) Reco_pol_pT[iFrame][iPTBin][phiPol]->Clone(name);
      // hData_pol_pT[iFrame][iPTBin][phiPol]->Divide(hAcc_pol_pT[iFrame][iPTBin][phiPol]);
      printf("frame %d, pT %d, acceptance correcting 1D now\n", iFrame, iPTBin);
      AccCorrect(hData_pol_pT[iFrame][iPTBin][phiPol], gAcc_pol_pT[iFrame][iPTBin][phiPol]);
      hData_pol_pT[iFrame][iPTBin][phiPol]->SetLineColor(colour_pT[iPTBin]);
      hData_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerColor(colour_pT[iPTBin]);
      hData_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerStyle(marker_pT[iPTBin]);
      
      sprintf(name, "hData2D_Onia_%s_pT%d", frameLabel[iFrame], iPTBin);
      hData2D_pol_pT[iFrame][iPTBin] = (TH2D *) Reco2D_pol_pT[iFrame][iPTBin]->Clone(name);
      hData2D_pol_pT[iFrame][iPTBin]->Divide(hAcc2D_pol_pT[iFrame][iPTBin]);
      hData2D_pol_pT[iFrame][iPTBin]->SetLineColor(colour_pT[iPTBin]);
      hData2D_pol_pT[iFrame][iPTBin]->SetMarkerColor(colour_pT[iPTBin]);
      hData2D_pol_pT[iFrame][iPTBin]->SetMarkerStyle(marker_pT[iPTBin]);
    }
    for(int iRapBin = 0; iRapBin < 2*kNbRapBins+1; iRapBin++){
      sprintf(name, "hData_Onia_cosTh_%s_rap%d", frameLabel[iFrame], iRapBin);
      hData_pol_rap[iFrame][iRapBin][cosThPol] = (TH1D *) Reco_pol_rap[iFrame][iRapBin][cosThPol]->Clone(name);
      // hData_pol_rap[iFrame][iRapBin][cosThPol]->Divide(hAcc_pol_rap[iFrame][iRapBin][cosThPol]);
      printf("frame %d, rap %d, acceptance correcting 1D now\n", iFrame, iRapBin);
      AccCorrect(hData_pol_rap[iFrame][iRapBin][cosThPol], gAcc_pol_pT[iFrame][iRapBin][cosThPol]);
      hData_pol_rap[iFrame][iRapBin][cosThPol]->SetLineColor(colour_rap[iRapBin]);
      hData_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerColor(colour_rap[iRapBin]);
      hData_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerStyle(marker_rap[iRapBin]);
      
      sprintf(name, "hData_Onia_phi_%s_rap%d", frameLabel[iFrame], iRapBin);
      hData_pol_rap[iFrame][iRapBin][phiPol] = (TH1D *) Reco_pol_rap[iFrame][iRapBin][phiPol]->Clone(name);
      // hData_pol_rap[iFrame][iRapBin][phiPol]->Divide(hAcc_pol_rap[iFrame][iRapBin][phiPol]);
      printf("frame %d, rap %d, acceptance correcting 1D now\n", iFrame, iRapBin);
      AccCorrect(hData_pol_rap[iFrame][iRapBin][phiPol], gAcc_pol_pT[iFrame][iRapBin][phiPol]);
      hData_pol_rap[iFrame][iRapBin][phiPol]->SetLineColor(colour_rap[iRapBin]);
      hData_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerColor(colour_rap[iRapBin]);
      hData_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerStyle(marker_rap[iRapBin]);
      
      sprintf(name, "hData2D_Onia_%s_rap%d", frameLabel[iFrame], iRapBin);
      hData2D_pol_rap[iFrame][iRapBin] = (TH2D *) Reco2D_pol_rap[iFrame][iRapBin]->Clone(name);
      hData2D_pol_rap[iFrame][iRapBin]->Divide(hAcc2D_pol_rap[iFrame][iRapBin]);
      hData2D_pol_rap[iFrame][iRapBin]->SetLineColor(colour_rap[iRapBin]);
      hData2D_pol_rap[iFrame][iRapBin]->SetMarkerColor(colour_rap[iRapBin]);
      hData2D_pol_rap[iFrame][iRapBin]->SetMarkerStyle(marker_rap[iRapBin]);
    }   
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
	sprintf(name, "hData_Onia_cosTh_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol] = (TH1D *) Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->Clone(name);
	//hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->Divide(hAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]);
	printf("frame %d, pT %d, rap %d acceptance correcting 1D now\n", iFrame, iPTBin, iRapBin);
	AccCorrect(hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol], gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetLineColor(colour_pT[iPTBin]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerColor(colour_pT[iPTBin]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerStyle(marker_pT[iPTBin]);
	
	sprintf(name, "hData_Onia_phi_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol] = (TH1D *) Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->Clone(name);
	//hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->Divide(hAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]);
	printf("frame %d, pT %d, rap %d acceptance correcting 1D now\n", iFrame, iPTBin, iRapBin);
	AccCorrect(hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol], gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetLineColor(colour_pT[iPTBin]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerColor(colour_pT[iPTBin]);
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerStyle(marker_pT[iPTBin]);

	sprintf(name, "hData2D_Onia_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Clone(name);
	// 	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Divide(hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]);//H: will be done below
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetLineColor(colour_pT[iPTBin]);
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerColor(colour_pT[iPTBin]);
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerStyle(marker_pT[iPTBin]);
      }
    }
  }
  
  Double_t binContentReco, binErrorReco;
  Double_t acc, binContentData, binErrorData;
  Double_t content, error, maxAcc;

  //acceptance correct the 2D histogram
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
	for(int iBinX = 1; iBinX <= hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(); iBinX++){
	  for(int iBinY = 1; iBinY <= hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY(); iBinY++){
	    
	    binContentReco = Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iBinX, iBinY);
	    binErrorReco = Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinError(iBinX, iBinY);
	    acc = hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iBinX, iBinY);

	    if(acc > 0.){
	      binErrorData = binErrorReco / acc;
	      binContentData = binContentReco / acc;
	    }
	    else{
	      binErrorData = 1.;
	      binContentData = 0.;
	    }
	    hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(iBinX, iBinY, binContentData);
	    hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinError(iBinX, iBinY, binErrorData);
	  }
	}
      }
    }
  }
  Get1DHistoFrom2D();
  SaveCorrData(label); //save the histos here BEFORE we set bin contents to zero for further processing!

}

//======================================
void ReadData(Char_t *fileNameData, Int_t rebinCosTh, Int_t rebinPhi){

  TFile *fIn = new TFile(fileNameData);
  
  Char_t name[100];
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
      sprintf(name, "Reco_Onia_cosTh_%s_pT%d", frameLabel[iFrame], iPTBin);
      Reco_pol_pT[iFrame][iPTBin][cosThPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_pT[iFrame][iPTBin][cosThPol]->SetLineColor(colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerColor(colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerStyle(marker_pT[iPTBin]);
      sprintf(name, "Reco_Onia_phi_%s_pT%d", frameLabel[iFrame], iPTBin);
      Reco_pol_pT[iFrame][iPTBin][phiPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_pT[iFrame][iPTBin][phiPol]->SetLineColor(colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerColor(colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerStyle(marker_pT[iPTBin]);
      sprintf(name, "Reco2D_Onia_%s_pT%d", frameLabel[iFrame], iPTBin);
      Reco2D_pol_pT[iFrame][iPTBin] = (TH2D *) gDirectory->Get(name);
      Reco2D_pol_pT[iFrame][iPTBin]->SetLineColor(colour_pT[iPTBin]);
      Reco2D_pol_pT[iFrame][iPTBin]->SetMarkerColor(colour_pT[iPTBin]);
      Reco2D_pol_pT[iFrame][iPTBin]->SetMarkerStyle(marker_pT[iPTBin]);
      if(rebinCosTh > 1 || rebinPhi > 1)
	Reco2D_pol_pT[iFrame][iPTBin]->Rebin2D(rebinCosTh, rebinPhi);
    }
    for(int iRapBin = 0; iRapBin < 2*kNbRapBins+1; iRapBin++){
      sprintf(name, "Reco_Onia_cosTh_%s_rap%d", frameLabel[iFrame], iRapBin);
      Reco_pol_rap[iFrame][iRapBin][cosThPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_rap[iFrame][iRapBin][cosThPol]->SetLineColor(colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerColor(colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerStyle(marker_rap[iRapBin]);
      sprintf(name, "Reco_Onia_phi_%s_rap%d", frameLabel[iFrame], iRapBin);
      Reco_pol_rap[iFrame][iRapBin][phiPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_rap[iFrame][iRapBin][phiPol]->SetLineColor(colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerColor(colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerStyle(marker_rap[iRapBin]);
      sprintf(name, "Reco2D_Onia_%s_rap%d", frameLabel[iFrame], iRapBin);
      Reco2D_pol_rap[iFrame][iRapBin] = (TH2D *) gDirectory->Get(name);
      Reco2D_pol_rap[iFrame][iRapBin]->SetLineColor(colour_rap[iRapBin]);
      Reco2D_pol_rap[iFrame][iRapBin]->SetMarkerColor(colour_rap[iRapBin]);
      Reco2D_pol_rap[iFrame][iRapBin]->SetMarkerStyle(marker_rap[iRapBin]);
      if(rebinCosTh > 1 || rebinPhi > 1)
	Reco2D_pol_rap[iFrame][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
    }
    Int_t totEv = 0;
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
	sprintf(name, "Reco_Onia_cosTh_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol] = (TH1D *) gDirectory->Get(name);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetLineColor(colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerColor(colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerStyle(marker_pT[iPTBin]);
	sprintf(name, "Reco_Onia_phi_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol] = (TH1D *) gDirectory->Get(name);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetLineColor(colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerColor(colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerStyle(marker_pT[iPTBin]);
	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) gDirectory->Get(name);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetLineColor(colour_pT[iPTBin]);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerColor(colour_pT[iPTBin]);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerStyle(marker_pT[iPTBin]);
	printf("%s: pTbin %d, rapBin %d --> statistics in data: %1.1f\n",
	       frameLabel[iFrame], iPTBin, iRapBin, Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Integral());

	totEv += Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Integral();
	if(rebinCosTh > 1 || rebinPhi > 1)
	  Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
      }
    }
    printf("%s: total nb. of events: %d\n\n\n", frameLabel[iFrame], totEv);
    
  }
}

//===================================
void ReadAccHistos(Char_t *fileNameMC){

  printf("reading in the acceptance histograms\n");
  TFile *fInMC = new TFile(fileNameMC);

  Char_t name[100];
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
      sprintf(name, "gAcc_Onia_cosTh_%s_pT%d", frameLabel[iFrame], iPTBin);
      gAcc_pol_pT[iFrame][iPTBin][cosThPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
      // printf("frame %d, pT %d, histo %p\n", iFrame, iPTBin, gAcc_pol_pT[iFrame][iPTBin][cosThPol]);
      //
      sprintf(name, "gAcc_Onia_phi_%s_pT%d", frameLabel[iFrame], iPTBin);
      gAcc_pol_pT[iFrame][iPTBin][phiPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
      // printf("frame %d, pT %d, histo %p\n", iFrame, iPTBin, gAcc_pol_pT[iFrame][iPTBin][phiPol]);
      //
      sprintf(name, "hAcc2D_Onia_%s_pT%d", frameLabel[iFrame], iPTBin);
      hAcc2D_pol_pT[iFrame][iPTBin] = (TH2D *) gDirectory->Get(name);
      // printf("frame %d, pT %d, histo %p\n", iFrame, iPTBin, hAcc2D_pol_pT[iFrame][iPTBin]);
    }
    for(int iRapBin = 0; iRapBin < 2*kNbRapBins+1; iRapBin++){
      sprintf(name, "gAcc_Onia_cosTh_%s_rap%d", frameLabel[iFrame], iRapBin);
      gAcc_pol_rap[iFrame][iRapBin][cosThPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
      // printf("frame %d, rap %d, histo %p\n", iFrame, iRapBin, gAcc_pol_rap[iFrame][iRapBin][cosThPol]);
      //
      sprintf(name, "gAcc_Onia_phi_%s_rap%d", frameLabel[iFrame], iRapBin);
      gAcc_pol_rap[iFrame][iRapBin][phiPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
      // printf("frame %d, rap %d, histo %p\n", iFrame, iRapBin, gAcc_pol_rap[iFrame][iRapBin][phiPol]);
      //
      sprintf(name, "hAcc2D_Onia_%s_rap%d", frameLabel[iFrame], iRapBin);
      hAcc2D_pol_rap[iFrame][iRapBin] = (TH2D *) gDirectory->Get(name);
      // printf("frame %d, rap %d, histo %p\n", iFrame, iRapBin, hAcc2D_pol_rap[iFrame][iRapBin]);
      //
    }
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
	sprintf(name, "gAcc_Onia_cosTh_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
	// printf("frame %d, pT %d, rap %d, histo %p\n", iFrame, iPTBin, iRapBin, gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]);
	//
	sprintf(name, "gAcc_Onia_phi_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol] = (TGraphAsymmErrors *) gDirectory->Get(name);
	// printf("frame %d, pT %d, rap %d, histo %p\n", iFrame, iPTBin, iRapBin, gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]);
	//
	sprintf(name, "hAcc2D_Onia_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) gDirectory->Get(name);
	// printf("frame %d, pT %d, rap %d, histo %p\n", iFrame, iPTBin, iRapBin, hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]);
      }
    }
  }
}

//============================================
void SaveCorrData(Char_t *label){

  Char_t name[100];
  sprintf(name, "polarization_corrData_%s.root", label);
  TFile *fOut = new TFile(name, "RECREATE");
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
      hData_pol_pT[iFrame][iPTBin][cosThPol]->Write();
      hData_pol_pT[iFrame][iPTBin][phiPol]->Write();
      hData2D_pol_pT[iFrame][iPTBin]->Write();
    }
    for(int iRapBin = 0; iRapBin < 2*kNbRapBins+1; iRapBin++){
      hData_pol_rap[iFrame][iRapBin][cosThPol]->Write();
      hData_pol_rap[iFrame][iRapBin][phiPol]->Write();
      hData2D_pol_rap[iFrame][iRapBin]->Write();
    }
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->Write();
	hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->Write();
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
