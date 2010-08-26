#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

TH1F *hMass[kNbPTBins+1][kNbRapForPTBins+1];
Double_t nSigma = 2.5;
void ReadInHistos(Char_t *fileNameIn);
void FitJPsi(Int_t iPTBin, Int_t iRapBin);
Double_t CalcSigOvBG(TF1 *func);
//====================================
//usage: root plotMass.C+ or
//       root 'plotMass.C+("myHistoInputFile.root")'
//====================================
void plotMass(Char_t *fileNameIn = "pol_data_HLT_Mu0Track0Jpsi.root"){

  ReadInHistos(fileNameIn);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++)
    for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++)
      FitJPsi(iPTBin, iRapBin);

}

//==========================================
void FitJPsi(Int_t thePT, Int_t theRap){

  printf("will fit pT %d and rap %d\n", thePT, theRap);
  gStyle->SetOptFit(kTRUE);
  gStyle->SetOptStat(kFALSE);
  Char_t name[100];

  Bool_t performFit = kTRUE;
  if(hMass[thePT][theRap]->Integral() < 30)
    performFit = kFALSE;

  TLatex *tex1;
  Double_t sigOvBG;
  sprintf(name, "c1_pT%d_rap%d", thePT, theRap);
  TCanvas *c1 = new TCanvas(name, "fits");
  //H: MIND that if the bin width is not 10 MeV anymore, the
  //scaling factor of 0.1 must be adjusted!!!
  hMass[thePT][theRap]->Rebin(2);
  TF1 *fRECO;
  if(performFit){

    sprintf(name, "fit_pT%d_rap%d", thePT, theRap);
    fRECO = new TF1(name, "0.2*pol1(0) + 0.2*gaus(2)", 2.8, 3.4);
    fRECO->SetParameters(25.,-0.01,100,3.097,0.04);
    //     fRECO->FixParameter(4, 0.035);
    fRECO->SetParLimits(4, 0., 1.);
  
    fRECO->SetParNames("d", "k", "N_{sig}", "#mu", "#sigma_{M}");
    // fRECO->SetLineColor(colour);
    //     hMass[thePT][theRap]->SetAxisRange(2.6, 3.6);
    hMass[thePT][theRap]->Fit(fRECO, "", "", 2.9, 3.3);
    fRECO = hMass[thePT][theRap]->GetFunction(name);
  }
  else
    hMass[thePT][theRap]->Draw("");

  //add some text
  if(theRap == 0 && thePT == 0) sprintf(name, "|y| < 2.3, all p_{T}");
  else if(theRap == 0) sprintf(name, "|y| < 2.3, %1.1f < p_{T} < %1.1f", pTRange[thePT-1], pTRange[thePT]);
  else if(theRap > 0 && thePT == 0) sprintf(name, "%1.1f < |y| < %1.1f, all p_{T}", rapForPTRange[theRap-1], rapForPTRange[theRap]);
  else if(theRap > 0) sprintf(name, "%1.1f < |y| < %1.1f, %1.1f < p_{T} < %1.1f", 
			      rapForPTRange[theRap-1], rapForPTRange[theRap], pTRange[thePT-1], pTRange[thePT]);

  tex1 = new TLatex(2.73, 0.9*hMass[thePT][theRap]->GetMaximum(), name);
  tex1->Draw();
  tex1->DrawLatex(2.73, 0.8*hMass[thePT][theRap]->GetMaximum(), "HLT_Mu0_TkMu0_Jpsi");

  if(performFit){
    sigOvBG = CalcSigOvBG(fRECO);
    sprintf(name, "S / B (#pm %1.1f #sigma_{M}) = %1.1f", nSigma, sigOvBG);
    tex1->DrawLatex(2.73, 0.7*hMass[thePT][theRap]->GetMaximum(), name);
  }
  //
  sprintf(name, "Figures/fitData_Jpsi_pT%d_rap%d.eps", thePT, theRap); c1->Print(name);
  sprintf(name, "Figures/fitData_Jpsi_pT%d_rap%d.gif", thePT, theRap); c1->Print(name);
  sprintf(name, "Figures/fitData_Jpsi_pT%d_rap%d.pdf", thePT, theRap); c1->Print(name);
}

//====================================
void ReadInHistos(Char_t *fileNameIn){

  TFile *fOut = new TFile(fileNameIn);
  Char_t name[100];
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
    for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++){
      //Mass:
      sprintf(name, "Reco_Onia_mass_pT%d_rap%d", iPTBin, iRapBin);
      hMass[iPTBin][iRapBin] = (TH1F *) gDirectory->Get(name);
      printf("pT %d, rap %d, histo %p\n", iPTBin, iRapBin, hMass[iPTBin][iRapBin]);
    }
  }
}

//=============================
Double_t CalcSigOvBG(TF1 *func){

  TF1 *gaus = new TF1("gaus", "gaus(0)", 2., 4.);
  Double_t mean = func->GetParameter(3);
  Double_t sigma = func->GetParameter(4);
  gaus->FixParameter(0, func->GetParameter(2));
  gaus->FixParameter(1, mean);
  gaus->FixParameter(2, sigma);

  Double_t signal = gaus->Integral(mean - nSigma*sigma, mean + nSigma*sigma);
  printf("signal = %f (within %f and %f GeV)\n", signal, mean - nSigma*sigma, mean + nSigma*sigma);

  TF1 *pol1 = new TF1("bg", "pol1", 2., 4.);
  pol1->FixParameter(0, func->GetParameter(0));
  pol1->FixParameter(1, func->GetParameter(1));

  Double_t bgContr = pol1->Integral(mean - nSigma*sigma, mean + nSigma*sigma);
  printf("bg = %f\n", bgContr);

  delete gaus;
  delete pol1;

  return signal / bgContr;
}
