#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

Double_t fitRangeMin[kNbRapForPTBins+1] = {2.8, 2.9, 2.85, 2.75};
Double_t fitRangeMax[kNbRapForPTBins+1] = {3.4, 3.3, 3.35, 3.45};
Double_t theChi2[kNbPTBins+1][kNbRapForPTBins+1];
Int_t theNDF[kNbPTBins+1][kNbRapForPTBins+1];
Double_t theMean[kNbPTBins+1][kNbRapForPTBins+1];
Double_t theSigma[kNbPTBins+1][kNbRapForPTBins+1];
Double_t theNSig[kNbPTBins+1][kNbRapForPTBins+1];

TH1F *hMass[kNbPTBins+1][kNbRapForPTBins+1];
TF1 *fRECO[kNbPTBins+1][kNbRapForPTBins+1];
Double_t binWidth = 0.01; //10 MeV
Double_t nSigma = 2.5;
void ReadInHistos(Char_t *fileNameIn);
void FitJPsi(Int_t iPTBin, Int_t iRapBin, Int_t fitFunc);
Double_t fitGaussLin(Double_t *x, Double_t *par);
Double_t fitCBLin(Double_t *x, Double_t *par);
Double_t CalcSigOvBG(TF1 *func);
void PrintFitPar(Char_t *fileNameOut);
//====================================
//usage: root plotMass.C+ or
//       root 'plotMass.C+("myHistoInputFile.root")'
//====================================
void plotMass(Char_t *fileNameIn = "pol_data_HLT_Mu0TkMu0Jpsi.root",
	      Char_t *fileFitResultsOut = "Results/dimuMassFitPar.txt",
	      Int_t fitFunc = 0 //[0]...Gauss+Lin, [1]...CB+Lin
	      ){

  ReadInHistos(fileNameIn);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++)
    for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++)
      FitJPsi(iPTBin, iRapBin, fitFunc);
  PrintFitPar(fileFitResultsOut);
}

//==========================================
void FitJPsi(Int_t thePT, Int_t theRap, Int_t fitFunc){

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
  if(performFit){

    sprintf(name, "fit_pT%d_rap%d", thePT, theRap);
    if(fitFunc == 0){
      fRECO[thePT][theRap] = new TF1(name, fitGaussLin, 2.8, 3.4, 5);
      fRECO[thePT][theRap]->SetParameters(250.,-0.1, 10000, 3.097, 0.04);
      fRECO[thePT][theRap]->SetParLimits(4, 0., 1.);  
      fRECO[thePT][theRap]->SetParNames("d", "k", "N_{sig}", "#mu", "#sigma_{M}");
      // // fRECO[thePT][theRap]->SetLineColor(colour);
      // //     hMass[thePT][theRap]->SetAxisRange(2.6, 3.6);
    }
    else if(fitFunc == 1){
      fRECO[thePT][theRap] = new TF1(name, fitCBLin, 2.8, 3.4, 7);
      fRECO[thePT][theRap]->SetParameters(250.,10., 10000, 3.097, 0.04, 1.6, 5.7);
      fRECO[thePT][theRap]->SetParLimits(4, 0., 1.);
      // fRECO[thePT][theRap]->FixParameter(5, 1.6);
      // fRECO[thePT][theRap]->FixParameter(6, 5.7);
      fRECO[thePT][theRap]->SetParNames("d", "k", "N", "#mu", "#sigma_{M}", "#alpha", "n");
    }

    binWidth = hMass[thePT][theRap]->GetBinWidth(1);
    printf("binWidth is %f\n", binWidth);

    hMass[thePT][theRap]->Fit(fRECO[thePT][theRap], "+", "", fitRangeMin[theRap], fitRangeMax[theRap]);
    printf("name of function is %s\n", name);
    fRECO[thePT][theRap] = hMass[thePT][theRap]->GetFunction(name);
    theChi2[thePT][theRap] = fRECO[thePT][theRap]->GetChisquare();
    theNDF[thePT][theRap] = fRECO[thePT][theRap]->GetNDF();
    theMean[thePT][theRap] = fRECO[thePT][theRap]->GetParameter(3);
    theSigma[thePT][theRap] = fRECO[thePT][theRap]->GetParameter(4);
    printf("sigma = %e\n", theSigma[thePT][theRap]);
    theNSig[thePT][theRap] = fRECO[thePT][theRap]->GetParameter(2);
  }
  else
    hMass[thePT][theRap]->Draw("");

  //add some text
  if(theRap == 0 && thePT == 0) sprintf(name, "|y| < 2.3, all p_{T}");
  else if(theRap == 0) sprintf(name, "|y| < 2.3, %1.1f < p_{T} < %1.1f", pTRange[thePT-1], pTRange[thePT]);
  else if(theRap > 1 && thePT == 0) sprintf(name, "%1.1f < |y| < %1.1f, all p_{T}", rapForPTRange[theRap-1], rapForPTRange[theRap]);
  else if(theRap == 1 && thePT == 0) sprintf(name, "|y| < %1.1f, all p_{T}", rapForPTRange[theRap]);
  else if(theRap == 1) sprintf(name, "|y| < %1.1f, %1.1f < p_{T} < %1.1f", rapForPTRange[theRap], pTRange[thePT-1], pTRange[thePT]);
  else if(theRap > 1)  sprintf(name, "%1.1f < |y| < %1.1f, %1.1f < p_{T} < %1.1f", 
	    rapForPTRange[theRap-1], rapForPTRange[theRap], pTRange[thePT-1], pTRange[thePT]);

  tex1 = new TLatex(2.73, 0.9*hMass[thePT][theRap]->GetMaximum(), name);
  tex1->Draw();
  tex1->DrawLatex(2.73, 0.8*hMass[thePT][theRap]->GetMaximum(), "HLT_Mu0_TkMu0_Jpsi");

  if(performFit){
    sigOvBG = CalcSigOvBG(fRECO[thePT][theRap]);
    sprintf(name, "S / B (#pm %1.1f #sigma_{M}) = %1.1f", nSigma, sigOvBG);
    tex1->DrawLatex(2.73, 0.7*hMass[thePT][theRap]->GetMaximum(), name);
  }
  
  sprintf(name, "Figures/fitData_Jpsi_pT%d_rap%d.eps", thePT, theRap); c1->Print(name);
  sprintf(name, "Figures/fitData_Jpsi_pT%d_rap%d.gif", thePT, theRap); c1->Print(name);
  sprintf(name, "Figures/fitData_Jpsi_pT%d_rap%d.pdf", thePT, theRap); c1->Print(name);
}

//====================================
void ReadInHistos(Char_t *fileNameIn){

  TFile *fIn = new TFile(fileNameIn);
  Char_t name[100];
  printf("reading file %s\n", fileNameIn);
  for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
    for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++){
      //Mass:
      sprintf(name, "Reco_Onia_mass_pT%d_rap%d", iPTBin, iRapBin);
      hMass[iPTBin][iRapBin] = (TH1F *) gDirectory->Get(name);
      printf("pT %d, rap %d, histo %p has %f entries and integral %f\n", 
	     iPTBin, iRapBin, hMass[iPTBin][iRapBin], hMass[iPTBin][iRapBin]->GetEntries(), hMass[iPTBin][iRapBin]->Integral());
    }
  }
}

//=============================
Double_t CalcSigOvBG(TF1 *func){

  TF1 *gaus = new TF1("gaus", "gaus(0)", 2., 4.);
  Double_t mean = func->GetParameter(3);
  Double_t sigma = func->GetParameter(4);

  gaus->FixParameter(0, func->GetParameter(2) / (TMath::Sqrt(2.*TMath::Pi()) * sigma));
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

//=========================
Double_t fitGaussLin(Double_t *x, Double_t *par){

  Double_t pol1 = par[0] + x[0]*par[1];
  Double_t mean = par[3];
  Double_t sigma = par[4];

  Double_t gauss = par[2] / (TMath::Sqrt(2.*TMath::Pi()) * sigma) * TMath::Exp(-(pow(x[0]-mean,2)/(2.*sigma*sigma)));
  Double_t result = pol1 + gauss;
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
  fprintf(fOut, "pT [GeV/c]\tall y\t|y|<%1.1f  %1.1f<|y|<%1.1f %1.1f<|y|<%1.1f\n",
	  rapForPTRange[1], rapForPTRange[1], rapForPTRange[2],
	  rapForPTRange[2], rapForPTRange[3]);
  fprintf(fOut, "------------------------------------------------------\n");
  fprintf(fOut, "all pT   \t%5.1f\t%5.1f\t%5.1f\t\t%5.1f\n",
	  1000.*theSigma[0][0], 1000.*theSigma[0][1], 1000.*theSigma[0][2], 1000.*theSigma[0][3]);
  fprintf(fOut, "------------------------------------------------------\n");
  for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
    fprintf(fOut, "%4.1f-%4.1f\t%5.1f\t%5.1f\t%5.1f\t\t%5.1f\n",
	    pTRange[iPTBin-1], pTRange[iPTBin],
	    1000.*theSigma[iPTBin][0], 1000.*theSigma[iPTBin][1], 1000.*theSigma[iPTBin][2], 1000.*theSigma[iPTBin][3]);
  }


  fprintf(fOut, "\n\nJ/psi mass [MeV]:\n");
  fprintf(fOut, "==================\n\n");
  fprintf(fOut, "pT [GeV/c]\tall y\t|y|<%1.1f  %1.1f<|y|<%1.1f %1.1f<|y|<%1.1f\n",
	  rapForPTRange[1], rapForPTRange[1], rapForPTRange[2],
	  rapForPTRange[2], rapForPTRange[3]);
  fprintf(fOut, "------------------------------------------------------\n");
  fprintf(fOut, "all pT   \t%5.1f\t%5.1f\t%5.1f\t\t%5.1f\n",
	  1000.*theMean[0][0], 1000.*theMean[0][1], 1000.*theMean[0][2], 1000.*theMean[0][3]);
  fprintf(fOut, "------------------------------------------------------\n");
  for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
    fprintf(fOut, "%4.1f-%4.1f\t%5.1f\t%5.1f\t%5.1f\t\t%5.1f\n",
	    pTRange[iPTBin-1], pTRange[iPTBin],
	    1000.*theMean[iPTBin][0], 1000.*theMean[iPTBin][1], 
	    1000.*theMean[iPTBin][2], 1000.*theMean[iPTBin][3]);
  }

  fprintf(fOut, "\n\nnumber of fitted J/psi's:\n");
  fprintf(fOut, "=============================\n\n");
  fprintf(fOut, "pT [GeV/c]\tall y\t|y|<%1.1f  %1.1f<|y|<%1.1f %1.1f<|y|<%1.1f\n",
	  rapForPTRange[1], rapForPTRange[1], rapForPTRange[2],
	  rapForPTRange[2], rapForPTRange[3]);
  fprintf(fOut, "-----------------------------------------------------------\n");
  fprintf(fOut, "all pT   \t%6.0f\t%6.0f\t%6.0f\t\t%6.0f\n",
	  theNSig[0][0], theNSig[0][1], theNSig[0][2], theNSig[0][3]);
  fprintf(fOut, "-----------------------------------------------------------\n");
  for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
    fprintf(fOut, "%4.1f-%4.1f\t%6.0f\t%6.0f\t%6.0f\t\t%6.0f\n",
	    pTRange[iPTBin-1], pTRange[iPTBin],
	    theNSig[iPTBin][0], theNSig[iPTBin][1], theNSig[iPTBin][2], theNSig[iPTBin][3]);
  }

  fprintf(fOut, "\n\nquality of the fit (chi2/ndf):\n");
  fprintf(fOut, "==================================\n\n");
  fprintf(fOut, "pT [GeV/c]   \tall y\t\t|y|<%1.1f \t %1.1f<|y|<%1.1f\t %1.1f<|y|<%1.1f\n",
	  rapForPTRange[1], rapForPTRange[1], rapForPTRange[2],
	  rapForPTRange[2], rapForPTRange[3]);
  fprintf(fOut, "----------------------------------------------------------------------------\n");
  fprintf(fOut, "all pT   \t%5.1f /%3d\t%5.1f /%3d\t%5.1f /%3d\t%5.1f /%3d\n",
	  theChi2[0][0], theNDF[0][0], theChi2[0][1], theNDF[0][1], 
	  theChi2[0][2], theNDF[0][2], theChi2[0][3], theNDF[0][3]);
  fprintf(fOut, "----------------------------------------------------------------------------\n");
  for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
    fprintf(fOut, "%4.1f-%4.1f\t%5.1f /%3d\t%5.1f /%3d\t%5.1f /%3d\t%5.1f /%3d\n",
	    pTRange[iPTBin-1], pTRange[iPTBin],
	    theChi2[iPTBin][0], theNDF[iPTBin][0], theChi2[iPTBin][1], theNDF[iPTBin][1], 
	    theChi2[iPTBin][2], theNDF[iPTBin][2], theChi2[iPTBin][3], theNDF[iPTBin][3]);

  }

  fclose(fOut);
}
