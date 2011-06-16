#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"
#include "TPaveStats.h"

// Double_t fitRangeMin[jpsi::kNbRapForPTBins+1] = {2.8, 2.9, 2.85, 2.75};
// Double_t fitRangeMax[jpsi::kNbRapForPTBins+1] = {3.4, 3.3, 3.35, 3.45};
Double_t fitRangeMin[jpsi::kNbRapForPTBins+1] = {2.6, 2.6, 2.6, 2.6, 2.6};
Double_t fitRangeMax[jpsi::kNbRapForPTBins+1] = {3.5, 3.5, 3.5, 3.5, 3.5};
Double_t theChi2[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
Int_t theNDF[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
Double_t theMean[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
Double_t theMeanErr[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
Double_t theSigma[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
Double_t theSigmaErr[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
Double_t theNSig[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
Double_t theNSigErr[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
Double_t theFracBG[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
Double_t theFracBGErr[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];

TH1F *hMass[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TF1 *fRECO[jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
Double_t binWidth = 0.01; //10 MeV
Double_t nSigma = 2.0; //2.5;
void ReadInHistos(Char_t *fileNameIn);
void FitJPsi(Int_t iPTBin, Int_t iRapBin, Int_t fitFunc, Char_t *hltTag);
Double_t fitGaussLin(Double_t *x, Double_t *par);
Double_t fitGaussExp(Double_t *x, Double_t *par);
Double_t fitGaussPol2(Double_t *x, Double_t *par);
Double_t fitCBLin(Double_t *x, Double_t *par);
Double_t CalcSigOvBG(TF1 *func, Int_t fitFunc);
void PrintFitPar(Char_t *fileNameOut);
//====================================
//usage: root plotMass.C+ or
//       root 'plotMass.C+("myHistoInputFile.root")'
//====================================
void plotMass(Char_t *fileNameIn = "pol_data_HLT_Mu0TkMu0Jpsi.root",
	      Char_t *fileFitResultsOut = "Results/dimuMassFitPar.txt",
	      Char_t *hltTag = "HLT_Mu0TkMu0Jpsi",
	      Int_t fitFunc = 0 //[0]...Gauss+Lin, [1]... Gauss+Exp, [2]...Gauss+Pol2, [2]...CB+Lin
	      ){

  ReadInHistos(fileNameIn);
  for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++)
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++)

  // for(int iPTBin = 3; iPTBin < 4; iPTBin++)
  //   for(int iRapBin = 1; iRapBin < 2; iRapBin++)
      FitJPsi(iPTBin, iRapBin, fitFunc, hltTag);
  PrintFitPar(fileFitResultsOut);
}

//==========================================
void FitJPsi(Int_t thePT, Int_t theRap, Int_t fitFunc, Char_t *hltTag){

  printf("will fit pT %d and rap %d\n", thePT, theRap);
  gStyle->SetOptFit(kTRUE);
  gStyle->SetOptStat("ne");
  // gStyle->SetOptStat(0);
  TPaveStats *st;
  Char_t name[100];

  Bool_t performFit = kTRUE;
  Int_t binMin = hMass[thePT][theRap]->GetXaxis()->FindBin(3.0);
  Int_t binMax = hMass[thePT][theRap]->GetXaxis()->FindBin(3.2);
  if(hMass[thePT][theRap]->Integral(binMin, binMax) < 30)
    performFit = kFALSE;

  TLatex *tex1;
  Double_t sigOvBG;
  sprintf(name, "c1_pT%d_rap%d", thePT, theRap);
  TCanvas *c1 = new TCanvas(name, "fits");
  TH1F *hFrame;
  if(performFit){

    sprintf(name, "fit_pT%d_rap%d", thePT, theRap);
    if(fitFunc == 0){
      printf("fitting with a Gaussian + a linear\n");
      fRECO[thePT][theRap] = new TF1(name, fitGaussLin, fitRangeMin[theRap], fitRangeMax[theRap], 5);
      fRECO[thePT][theRap]->SetParameters(10000, 3.097, 0.04, 250.,-0.1);
      fRECO[thePT][theRap]->SetParLimits(2, 0., 1.);  
      fRECO[thePT][theRap]->SetParNames("N_{sig}", "#mu", "#sigma_{M}", "d", "k");
      // // fRECO[thePT][theRap]->SetLineColor(colour);
      // //     hMass[thePT][theRap]->SetAxisRange(2.6, 3.6);
    }
    if(fitFunc == 1){
      printf("fitting with a Gaussian + an exp\n");
      fRECO[thePT][theRap] = new TF1(name, fitGaussExp, fitRangeMin[theRap], fitRangeMax[theRap], 5);
      fRECO[thePT][theRap]->SetParameters(10000, 3.097, 0.04, 250.,-0.1);
      fRECO[thePT][theRap]->SetParLimits(2, 0., 1.);  
      fRECO[thePT][theRap]->SetParNames("N_{sig}", "#mu", "#sigma_{M}", "d", "k");
      // // fRECO[thePT][theRap]->SetLineColor(colour);
      // //     hMass[thePT][theRap]->SetAxisRange(2.6, 3.6);
    }
    else if(fitFunc == 2){
      printf("fitting with a Gaussian + a Pol2\n");
      fRECO[thePT][theRap] = new TF1(name, fitGaussPol2, fitRangeMin[theRap], fitRangeMax[theRap], 6);
      fRECO[thePT][theRap]->SetParameters(10000, 3.097, 0.04, 250.,-0.1, 1.);
      fRECO[thePT][theRap]->SetParLimits(2, 0., 1.);  
      fRECO[thePT][theRap]->SetParNames("N_{sig}", "#mu", "#sigma_{M}", "p0", "p1", "p2");
      // // fRECO[thePT][theRap]->SetLineColor(colour);
      // //     hMass[thePT][theRap]->SetAxisRange(2.6, 3.6);
    }
    else if(fitFunc == 3){
      fRECO[thePT][theRap] = new TF1(name, fitCBLin, fitRangeMin[theRap], fitRangeMax[theRap], 7);
      fRECO[thePT][theRap]->SetParameters(250.,10., 10000, 3.097, 0.04, 1.6, 5.7);
      fRECO[thePT][theRap]->SetParLimits(4, 0., 1.);
      // fRECO[thePT][theRap]->FixParameter(5, 1.6);
      // fRECO[thePT][theRap]->FixParameter(6, 5.7);
      fRECO[thePT][theRap]->SetParNames("d", "k", "N", "#mu", "#sigma_{M}", "#alpha", "n");
    }

    binWidth = hMass[thePT][theRap]->GetBinWidth(1);
    printf("binWidth is %f\n", binWidth);

    hMass[thePT][theRap]->Fit(fRECO[thePT][theRap], "+0", "", 
			      fitRangeMin[theRap], fitRangeMax[theRap]);    

    fRECO[thePT][theRap] = hMass[thePT][theRap]->GetFunction(name);
    theChi2[thePT][theRap] = fRECO[thePT][theRap]->GetChisquare();
    theNDF[thePT][theRap] = fRECO[thePT][theRap]->GetNDF();

    theMean[thePT][theRap] = fRECO[thePT][theRap]->GetParameter(1);
    theMeanErr[thePT][theRap] = fRECO[thePT][theRap]->GetParError(1);
    theSigma[thePT][theRap] = fRECO[thePT][theRap]->GetParameter(2);
    theSigmaErr[thePT][theRap] = fRECO[thePT][theRap]->GetParError(2);
    theNSig[thePT][theRap] = fRECO[thePT][theRap]->GetParameter(0);
    theNSigErr[thePT][theRap] = fRECO[thePT][theRap]->GetParError(0);

    hMass[thePT][theRap]->Draw();
    gPad->Update();

    hFrame = gPad->DrawFrame(2.5, 0.8*hMass[thePT][theRap]->GetMinimum()+0.5, 4.1, 5.*hMass[thePT][theRap]->GetMaximum());
    gPad->SetLogy();
    hFrame->SetXTitle(hMass[thePT][theRap]->GetXaxis()->GetTitle());
    hFrame->Draw();

    hMass[thePT][theRap]->Draw("same");
    fRECO[thePT][theRap]->Draw("same");
    st = (TPaveStats *) hMass[thePT][theRap]->FindObject("stats");
    printf("statistics box: %p\n", st);
    st->SetOptStat(11);
    st->SetOptFit(1);
    st->Paint();
  }
  else
    hMass[thePT][theRap]->Draw("");

  //add some text
  if(theRap == 0 && thePT == 0) sprintf(name, "|y| < 2.4, all p_{T}");
  else if(theRap == 0) sprintf(name, "|y| < 2.4, %1.1f < p_{T} < %1.1f", jpsi::pTRange[theRap][thePT-1], jpsi::pTRange[theRap][thePT]);
  else if(theRap > 1 && thePT == 0) sprintf(name, "%1.1f < |y| < %1.1f, all p_{T}", jpsi::rapForPTRange[theRap-1], jpsi::rapForPTRange[theRap]);
  else if(theRap == 1 && thePT == 0) sprintf(name, "|y| < %1.1f, all p_{T}", jpsi::rapForPTRange[theRap]);
  else if(theRap == 1) sprintf(name, "|y| < %1.1f, %1.1f < p_{T} < %1.1f", jpsi::rapForPTRange[theRap], jpsi::pTRange[theRap][thePT-1], jpsi::pTRange[theRap][thePT]);
  else if(theRap > 1)  sprintf(name, "%1.1f < |y| < %1.1f, %1.1f < p_{T} < %1.1f", 
	    jpsi::rapForPTRange[theRap-1], jpsi::rapForPTRange[theRap], jpsi::pTRange[theRap][thePT-1], jpsi::pTRange[theRap][thePT]);

  //w/o Psi':
  // tex1 = new TLatex(2.73, 0.9*hMass[thePT][theRap]->GetMaximum(), name);
  // tex1->Draw();
  // tex1->DrawLatex(2.73, 0.8*hMass[thePT][theRap]->GetMaximum(), hltTag);

  //with Psi' and log scale
  tex1 = new TLatex(2.52, 3.*hMass[thePT][theRap]->GetMaximum(), name);
  tex1->Draw();
  tex1->DrawLatex(2.52, 1.9*hMass[thePT][theRap]->GetMaximum(), hltTag);

  if(performFit){
    sigOvBG = CalcSigOvBG(fRECO[thePT][theRap], fitFunc);
    sprintf(name, "S / B (#pm %1.1f #sigma_{M}) = %1.1f", nSigma, sigOvBG);
    // tex1->DrawLatex(2.73, 0.7*hMass[thePT][theRap]->GetMaximum(), name);
    tex1->DrawLatex(2.52, 1.3*hMass[thePT][theRap]->GetMaximum(), name);
    
    theFracBG[thePT][theRap] = 1./(sigOvBG + 1.);
  }
  
  //add the lines where we will base the mass cuts:
  Double_t jPsiMassMin = jpsi::polMassJpsi[theRap] - jpsi::nSigMass*jpsi::sigmaMassJpsi[theRap];
  Double_t jPsiMassMax = jpsi::polMassJpsi[theRap] + jpsi::nSigMass*jpsi::sigmaMassJpsi[theRap];
  TLine *line = new TLine(jPsiMassMin, 0., jPsiMassMin, 0.6*hMass[thePT][theRap]->GetMaximum());
  //line->SetLineStyle(3);
  line->Draw();
  line->DrawLine(jPsiMassMax, 0., jPsiMassMax, 0.6*hMass[thePT][theRap]->GetMaximum());

  sprintf(name, "Figures/fitData_Jpsi_pT%d_rap%d_%s.eps", thePT, theRap, hltTag); c1->Print(name);
  sprintf(name, "Figures/fitData_Jpsi_pT%d_rap%d_%s.gif", thePT, theRap, hltTag); c1->Print(name);
  sprintf(name, "Figures/fitData_Jpsi_pT%d_rap%d_%s.pdf", thePT, theRap, hltTag); c1->Print(name);
}

//====================================
void ReadInHistos(Char_t *fileNameIn){

  TFile *fIn = new TFile(fileNameIn);
  Char_t name[100];
  printf("reading file %s\n", fileNameIn);
  for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
      //Mass:
      sprintf(name, "Reco_Onia_mass_pT%d_rap%d", iPTBin, iRapBin);
      hMass[iPTBin][iRapBin] = (TH1F *) gDirectory->Get(name);
      printf("pT %d, rap %d, histo %p has %f entries and integral %f\n", 
	     iPTBin, iRapBin, hMass[iPTBin][iRapBin], hMass[iPTBin][iRapBin]->GetEntries(), hMass[iPTBin][iRapBin]->Integral());
    }
  }
}

//=============================
Double_t CalcSigOvBG(TF1 *func, Int_t fitFunc){

  TF1 *gaus = new TF1("gaus", "gaus(0)", 2., 4.);
  Double_t mean = func->GetParameter(1);
  Double_t sigma = func->GetParameter(2);

  gaus->FixParameter(0, func->GetParameter(0) / (TMath::Sqrt(2.*TMath::Pi()) * sigma));
  gaus->FixParameter(1, mean);
  gaus->FixParameter(2, sigma);

  Double_t signal = gaus->Integral(mean - nSigma*sigma, mean + nSigma*sigma);
  printf("signal = %f (within %f and %f GeV)\n", signal, mean - nSigma*sigma, mean + nSigma*sigma);

  TF1 *bgFunc;
  if(fitFunc == 0){
    bgFunc = new TF1("bg", "pol1", 2., 4.);
    bgFunc->FixParameter(0, func->GetParameter(3));
    bgFunc->FixParameter(1, func->GetParameter(4));
  }
  if(fitFunc == 1){
    bgFunc = new TF1("bg", "expo", 2., 4.);
    bgFunc->FixParameter(0, func->GetParameter(3));
    bgFunc->FixParameter(1, func->GetParameter(4));
  }
  else if(fitFunc == 2){
    bgFunc = new TF1("bg", "pol2", 2., 4.);
    bgFunc->FixParameter(0, func->GetParameter(3));
    bgFunc->FixParameter(1, func->GetParameter(4));
    bgFunc->FixParameter(2, func->GetParameter(4));
  }

  Double_t bgContr = bgFunc->Integral(mean - nSigma*sigma, mean + nSigma*sigma);
  printf("bg = %f\n", bgContr);

  delete gaus;
  delete bgFunc;

  return signal / bgContr;
}

//=========================
Double_t fitGaussLin(Double_t *x, Double_t *par){

  Double_t pol1 = par[3] + x[0]*par[4];
  Double_t mean = par[1];
  Double_t sigma = par[2];

  Double_t gauss = par[0] / (TMath::Sqrt(2.*TMath::Pi()) * sigma) * TMath::Exp(-(pow(x[0]-mean,2)/(2.*sigma*sigma)));
  Double_t result = pol1 + gauss;
  result *= binWidth; //correct for the bin width
  return result;
}

//=========================
Double_t fitGaussExp(Double_t *x, Double_t *par){

  Double_t exp = par[3]*TMath::Exp(-par[4]);
  Double_t mean = par[1];
  Double_t sigma = par[2];

  Double_t gauss = par[0] / (TMath::Sqrt(2.*TMath::Pi()) * sigma) * TMath::Exp(-(pow(x[0]-mean,2)/(2.*sigma*sigma)));
  Double_t result = exp + gauss;
  result *= binWidth; //correct for the bin width
  return result;
}

//=========================
Double_t fitGaussPol2(Double_t *x, Double_t *par){

  Double_t pol2 = par[3] + x[0]*par[4] + x[0]*x[0]*par[5];
  Double_t mean = par[1];
  Double_t sigma = par[2];

  Double_t gauss = par[0] / (TMath::Sqrt(2.*TMath::Pi()) * sigma) * TMath::Exp(-(pow(x[0]-mean,2)/(2.*sigma*sigma)));
  Double_t result = pol2 + gauss;
  result *= binWidth; //correct for the bin width
  return result;
}

//=========================
Double_t fitCBLin(Double_t *x, Double_t *par){

  Double_t pol1 = par[0] + x[0]*par[1];
  Double_t norm = par[2];
  Double_t mean = par[3];
  Double_t sigma = par[4];
  Double_t alpha = par[5];
  Double_t n = par[6];

  Double_t CB;
  if(((x[0] - mean)/sigma) > -alpha)
    CB = TMath::Exp(-(pow(x[0]-mean,2)/(2.*sigma*sigma)));
  else{
    Double_t A = pow(n / fabs(alpha),n) * TMath::Exp(-alpha*alpha/2.);
    Double_t B = n / fabs(alpha) - fabs(alpha);
    CB = A*pow(B - (x[0]-mean)/sigma, -n);
  }
  Double_t result = pol1 + norm*CB;
  result *= binWidth; //correct for the bin width
  return result;
}

//=========================
void PrintFitPar(Char_t *fileNameOut){

  FILE *fOut = fopen(fileNameOut, "write");

  fprintf(fOut, "mass resolution [MeV]:\n");
  fprintf(fOut, "======================\n\n");
  for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++){

    if(iRap == 0)
      fprintf(fOut, "pT [GeV/c]\t   all y\n");
    else if(iRap == 1)
      fprintf(fOut, "pT [GeV/c]       |y| < %1.1f\n", jpsi::rapForPTRange[1]);
    else
      fprintf(fOut, "pT [GeV/c]      %1.1f<|y|<%1.1f\n", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);
    fprintf(fOut, "---------------------------\n");

    fprintf(fOut, "all pT   \t%5.1f+-%3.1f\n", 1000.*theSigma[0][iRap], 1000.*theSigmaErr[0][iRap]);
    fprintf(fOut, "---------------------------\n");
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[iRap]+1; iPTBin++){
      fprintf(fOut, "%4.1f-%4.1f\t%5.1f+-%3.1f\n",
	      jpsi::pTRange[iRap][iPTBin-1], jpsi::pTRange[iRap][iPTBin],
	      1000.*theSigma[iPTBin][iRap], 1000.*theSigmaErr[iPTBin][iRap]);
    }
    fprintf(fOut, "---------------------------\n\n");
  }

  fprintf(fOut, "\n\nJ/psi mass [MeV]:\n");
  fprintf(fOut, "==================\n\n");
  for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++){

    if(iRap == 0)
      fprintf(fOut, "pT [GeV/c]\t   all y\n");
    else if(iRap == 1)
      fprintf(fOut, "pT [GeV/c]       |y| < %1.1f\n", jpsi::rapForPTRange[1]);
    else
      fprintf(fOut, "pT [GeV/c]      %1.1f<|y|<%1.1f\n", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);
    fprintf(fOut, "---------------------------\n");

    fprintf(fOut, "all pT   \t%5.1f+-%3.1f\n", 1000.*theMean[0][iRap], 1000.*theMeanErr[0][iRap]);
    fprintf(fOut, "---------------------------\n");
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[iRap]+1; iPTBin++){
      fprintf(fOut, "%4.1f-%4.1f\t%5.1f+-%3.1f\n",
	      jpsi::pTRange[iRap][iPTBin-1], jpsi::pTRange[iRap][iPTBin],
	      1000.*theMean[iPTBin][iRap], 1000.*theMeanErr[iPTBin][iRap]);
    }
    fprintf(fOut, "---------------------------\n\n");
  }

  fprintf(fOut, "\n\nnumber of fitted J/psi's and fraction of BG:\n");
  fprintf(fOut, "================================================\n\n");
  for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++){

    if(iRap == 0)
      fprintf(fOut, "pT [GeV/c]\t   all y\t BG/tot\n");
    else if(iRap == 1)
      fprintf(fOut, "pT [GeV/c]       |y| < %1.1f\t BG/tot\n", jpsi::rapForPTRange[1]);
    else
      fprintf(fOut, "pT [GeV/c]      %1.1f<|y|<%1.1f\t BG/tot\n", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);
    fprintf(fOut, "----------------------------------------\n");

    fprintf(fOut, "all pT   \t%5.1f+-%3.1f\n", theNSig[0][iRap], theNSigErr[0][iRap]);
    fprintf(fOut, "----------------------------------------\n");
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[iRap]+1; iPTBin++){
      fprintf(fOut, "%4.1f-%4.1f\t%6.0f+-%4.0f\t%1.3f\n",
	      jpsi::pTRange[iRap][iPTBin-1], jpsi::pTRange[iRap][iPTBin],	      
	      theNSig[iPTBin][iRap], theNSigErr[iPTBin][iRap], theFracBG[iPTBin][iRap]);
    }
    fprintf(fOut, "-----------------------------------------\n\n");
  }

  // fprintf(fOut, "\n\nJ/psi mass [MeV]:\n");
  // fprintf(fOut, "==================\n\n");
  // fprintf(fOut, "pT [GeV/c]\t   all y\t  |y|<%1.1f       %1.1f<|y|<%1.1f     %1.1f<|y|<%1.1f     %1.1f<|y|<%1.1f\n",
  // 	  jpsi::rapForPTRange[1], 
  // 	  jpsi::rapForPTRange[1], jpsi::rapForPTRange[2],
  // 	  jpsi::rapForPTRange[2], jpsi::rapForPTRange[3], 
  // 	  jpsi::rapForPTRange[3], jpsi::rapForPTRange[4]);
  // fprintf(fOut, "-------------------------------------------------------------------------------------------\n");
  // fprintf(fOut, "all pT   \t%5.1f+-%3.1f\t%5.1f+-%3.1f\t%5.1f+-%3.1f\t%5.1f+-%3.1f\t%5.1f+-%3.1f\n",
  // 	  1000.*theMean[0][0], 1000.*theMeanErr[0][0], 
  // 	  1000.*theMean[0][1], 1000.*theMeanErr[0][1], 
  // 	  1000.*theMean[0][2], 1000.*theMeanErr[0][2], 
  // 	  1000.*theMean[0][3], 1000.*theMeanErr[0][3],
  // 	  1000.*theMean[0][4], 1000.*theMeanErr[0][4]);
  // fprintf(fOut, "-------------------------------------------------------------------------------------------\n");
  // for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[1]+1; iPTBin++){
  //   fprintf(fOut, "%4.1f-%4.1f\t%5.1f+-%3.1f\t%5.1f+-%3.1f\t%5.1f+-%3.1f\t%5.1f+-%3.1f\t%5.1f+-%3.1f\n",
  // 	    jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin],
  // 	    1000.*theMean[iPTBin][0], 1000.*theMeanErr[iPTBin][0], 
  // 	    1000.*theMean[iPTBin][1], 1000.*theMeanErr[iPTBin][1], 
  // 	    1000.*theMean[iPTBin][2], 1000.*theMeanErr[iPTBin][2], 
  // 	    1000.*theMean[iPTBin][3], 1000.*theMeanErr[iPTBin][3],
  // 	    1000.*theMean[iPTBin][4], 1000.*theMeanErr[iPTBin][4]);
  // }

  // fprintf(fOut, "\n\nnumber of fitted J/psi's:\n");
  // fprintf(fOut, "=============================\n\n");
  // fprintf(fOut, "pT [GeV/c]\t      all y\t   |y|<%1.1f       %1.1f<|y|<%1.1f     %1.1f<|y|<%1.1f     %1.1f<|y|<%1.1f\n",
  // 	  jpsi::rapForPTRange[1], 
  // 	  jpsi::rapForPTRange[1], jpsi::rapForPTRange[2],
  // 	  jpsi::rapForPTRange[2], jpsi::rapForPTRange[3], 
  // 	  jpsi::rapForPTRange[3], jpsi::rapForPTRange[4]);
  // fprintf(fOut, "-------------------------------------------------------------------------------------------\n");
  // fprintf(fOut, "all pT   \t%6.0f+-%4.0f\t%6.0f+-%4.0f\t%6.0f+-%4.0f\t%6.0f+-%4.0f\t%6.0f+-%4.0f\n",
  // 	  theNSig[0][0], theNSigErr[0][0], 
  // 	  theNSig[0][1], theNSigErr[0][1], 
  // 	  theNSig[0][2], theNSigErr[0][2], 
  // 	  theNSig[0][3], theNSigErr[0][3],
  // 	  theNSig[0][4], theNSigErr[0][4]);
  // fprintf(fOut, "-------------------------------------------------------------------------------------------\n");
  // for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[1]+1; iPTBin++){
  //   fprintf(fOut, "%4.1f-%4.1f\t%6.0f+-%4.0f\t%6.0f+-%4.0f\t%6.0f+-%4.0f\t%6.0f+-%4.0f\t%6.0f+-%4.0f\n",
  // 	    jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin],
  // 	    theNSig[iPTBin][0], theNSigErr[iPTBin][0], 
  // 	    theNSig[iPTBin][1], theNSigErr[iPTBin][1], 
  // 	    theNSig[iPTBin][2], theNSigErr[iPTBin][2], 
  // 	    theNSig[iPTBin][3], theNSigErr[iPTBin][3],
  // 	    theNSig[iPTBin][4], theNSigErr[iPTBin][4]);
  // }

  // fprintf(fOut, "\n\nquality of the fit (chi2/ndf):\n");
  // fprintf(fOut, "==================================\n\n");
  // fprintf(fOut, "pT [GeV/c]   \tall y\t\t|y|<%1.1f \t %1.1f<|y|<%1.1f\t %1.1f<|y|<%1.1f\t %1.1f<|y|<%1.1f\n",
  // 	  jpsi::rapForPTRange[1], 
  // 	  jpsi::rapForPTRange[1], jpsi::rapForPTRange[2],
  // 	  jpsi::rapForPTRange[2], jpsi::rapForPTRange[3], 
  // 	  jpsi::rapForPTRange[3], jpsi::rapForPTRange[4]);
  // fprintf(fOut, "-------------------------------------------------------------------------------------------\n");
  // fprintf(fOut, "all pT   \t%5.1f /%3d\t%5.1f /%3d\t%5.1f /%3d\t%5.1f /%3d\t%5.1f /%3d\n",
  // 	  theChi2[0][0], theNDF[0][0], theChi2[0][1], theNDF[0][1], 
  // 	  theChi2[0][2], theNDF[0][2], theChi2[0][3], theNDF[0][3], 
  // 	  theChi2[0][4], theNDF[0][4]);
  // fprintf(fOut, "-------------------------------------------------------------------------------------------\n");
  // for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[1]+1; iPTBin++){
  //   fprintf(fOut, "%4.1f-%4.1f\t%5.1f /%3d\t%5.1f /%3d\t%5.1f /%3d\t%5.1f /%3d\t%5.1f /%3d\n",
  // 	    jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin],
  // 	    theChi2[iPTBin][0], theNDF[iPTBin][0], theChi2[iPTBin][1], theNDF[iPTBin][1], 
  // 	    theChi2[iPTBin][2], theNDF[iPTBin][2], theChi2[iPTBin][3], theNDF[iPTBin][3],
  // 	    theChi2[iPTBin][4], theNDF[iPTBin][4]);

  // }

  fclose(fOut);
}
