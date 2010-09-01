#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"
#include "TGraphAsymmErrors.h"

Double_t maxCosTh = 1.0;
Double_t minPhi = 0., maxPhi = 360.;
Double_t minAcc = 0.05;//default: 0.1
Double_t fracMaxAcc = 0.1;
Double_t minBinContent = 5.;

Double_t fracMax = 0.01;

Int_t nbBinsCosTheta, nbBinsPhi, nBins1D; //will be filled in Get1DHistoFrom2D

//polarization histos from Data:
TH1D *Reco_pol_pT[kNbFrames][kNbPTBins+1][kNbPolVar];
TH1D *Reco_pol_rap[kNbFrames][2*kNbRapBins+1][kNbPolVar];
TH1D *Reco_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1][kNbPolVar];
TH2D *Reco2D_pol_pT[kNbFrames][kNbPTBins+1];
TH2D *Reco2D_pol_rap[kNbFrames][2*kNbRapBins+1];
TH2D *Reco2D_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1];

//acceptance histos
TH1D *hAcc_pol_pT[kNbFrames][kNbPTBins+1][kNbPolVar];
TH1D *hAcc_pol_rap[kNbFrames][2*kNbRapBins+1][kNbPolVar];
TH1D *hAcc_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1][kNbPolVar];
TH2D *hAcc2D_pol_pT[kNbFrames][kNbPTBins+1];
TH2D *hAcc2D_pol_rap[kNbFrames][2*kNbRapBins+1];
TH2D *hAcc2D_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1];

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

//parameters of 1D fits from the 2D projection histo:
Double_t norm1D_rap_pT[kNbRapForPTBins+1][kNbPTBins+1];
Double_t lambda1D_rap_pT[kNbRapForPTBins+1][kNbPTBins+1][3];
Double_t errLambda1D_rap_pT[kNbRapForPTBins+1][kNbPTBins+1][3];
Double_t frameInv[kNbRapForPTBins+1][kNbPTBins+1];
Double_t chi2_1D_rap_pT[kNbRapForPTBins+1][kNbPTBins+1];
Double_t chi2_1D_rap_global[kNbRapForPTBins+1], probChi2_1D_rap_global;
Int_t NDF_1D_rap_pT[kNbRapForPTBins+1][kNbPTBins+1], NDF_1D_rap_global[kNbRapForPTBins+1];
Bool_t fit1DConverged[kNbRapForPTBins+1][kNbPTBins+1];
TF1 *f1D_rap_pT[kNbFrames][kNbRapForPTBins+1][kNbPTBins+1];
TF1 *plotF1D_rap_pT[kNbFrames][kNbRapForPTBins+1][kNbPTBins+1];

TGraphAsymmErrors *gLambdaTh1D_rap[kNbFrames][kNbRapForPTBins+1];
TGraphAsymmErrors *gLambdaPhi1D_rap[kNbFrames][kNbRapForPTBins+1];
TGraphAsymmErrors *gLambdaThPhi1D_rap[kNbFrames][kNbRapForPTBins+1];
TGraph *gFrameInv1D_rap[kNbFrames][kNbRapForPTBins+1];

void ReadData(Char_t *fileNameData, Int_t rebinCosTh, Int_t rebinPhi);
void ReadAccHistos(Char_t *fileNameMC, Bool_t normalise);
void CorrectForAcc(Char_t *polTag);
void PlotUncorrData(Int_t iFrame, Char_t *polTag);
void PlotUncorrData2D(Int_t iFrame, Char_t *polTag);
void PlotCorrData2D(Int_t iFrame, Char_t *polTag);
void Get1DHistoFrom2D();
void Fit1DHisto(Int_t iFrame, Char_t *polTag, Int_t rap=0, Int_t pT=0, Bool_t displayOn=kFALSE);
void PlotFitResults1D(TH2D *hHist, TF1 *fFit, Int_t iFrame, Int_t iRap, Int_t iPT, Int_t binPhi, Char_t *label);
void SaveCorrData(Char_t *polTag);
void Write1DFitResults(Char_t *fOutName);
Double_t fit1D(Double_t *x, Double_t *par);
//======================================
void fitPolPar(Char_t *label = "data_HLT_Mu0TkMu0Jpsi_cut120_1Sep2010",
	       Int_t pT = 4,
	       Int_t rap = 3,
	       Bool_t normalise = kTRUE, //normalises the acceptance histograms to 1
	       Bool_t displayOn = kFALSE,
	       Int_t rebinCosTh = 2,
	       Int_t rebinPhi = 2){

  Char_t fileNameMC[200];
  sprintf(fileNameMC, "accHistos_HLT_Mu0Track0Jpsi_cut120_1Sep2010.root", label);

  Char_t fileNameData[200], fOutName1D[200];
  sprintf(fileNameData, "pol_%s.root", label);
  sprintf(fOutName1D, "Results/fitPar1D_lambda_%s.root", label);

  ReadData(fileNameData, rebinCosTh, rebinPhi);
//   PlotUncorrData(CS, label);
//   PlotUncorrData(HX, label);
  ReadAccHistos(fileNameMC, normalise);
  CorrectForAcc(label); //calls internally "Get1DHistoFrom2D" to fill
                         //the histo hData1D_pol_pT_rap[kNbFrames][kNbPTBins+1][kNbRapForPTBins+1];

  //1D analysis projected from 2D:
  Fit1DHisto(CS, label, rap, pT, displayOn);
  Fit1DHisto(HX, label, rap, pT, displayOn);

  if(rap == 0 && pT == 0)
    Write1DFitResults(fOutName1D);  //shall be run only when looping over all rap and pT bins

  //2D analysis: (should not be used)
//   PlotUncorrData2D(CS, label);
//   PlotUncorrData2D(HX, label);
//   PlotCorrData2D(CS, label);
//   PlotCorrData2D(HX, label);
}

//======================================
Double_t fit1D(Double_t *x, Double_t *par){

  if(x[0] < 0. || x[0] > (Double_t) nBins1D) return 0.;

  // Double_t binsCosTheta = 40.; //H: change the "20" manually if a different binsize is taken
  Int_t index = (Int_t) x[0];
  Int_t iPhi = 1 + Int_t(x[0] / nbBinsCosTheta); 
  Double_t phi = nbBinsCosTheta * (Float_t)iPhi - (nbBinsCosTheta/2.);
  Int_t iCosTh = 1+(Int_t)(x[0]- nbBinsCosTheta *(iPhi-1));
  Double_t cosTh = -1.05 + (Float_t)iCosTh / (nbBinsCosTheta/2.);

  // printf("index %d, phi %f (bin %d), cosTh %f (bin %d)\n", 
  // 	 index, phi, iPhi, cosTh, iCosTh);

  Double_t cosTh2 = cosTh*cosTh;

  Double_t norm = par[0];
  Double_t lambdaTh = par[1];
  Double_t lambdaPhi = par[2];
  Double_t lambdaThPhi = par[3];  
 
//   Int_t iBin = thisHist1D->GetXaxis()->FindBin(x[0]);
  Int_t iBin = (Int_t)x[0] + 1;
  Double_t content = thisHist1D->GetBinContent(iBin);
//   if(content < 1.0){
  if(content < 0.0){//reject all bins that had no events and those that were flagged as -1.

    // printf("rejecting point because content is -1\n");
    TF1::RejectPoint();
    return 0.;
  }

 Double_t conv = TMath::Pi() / 180.;
 Double_t termTheta = 1.+lambdaTh*cosTh2;
 Double_t termPhi = lambdaPhi*(1. - cosTh2)*TMath::Cos(2.*phi*conv);
 Double_t termThetaPhi = lambdaThPhi*2.*cosTh*sqrt(1.-cosTh2)*cos(phi*conv);
 Double_t funcValue = termTheta + termPhi + termThetaPhi;
 return norm*funcValue;
}


//======================================
void Fit1DHisto(Int_t iFrame, Char_t *label, Int_t rap, Int_t pT, Bool_t displayOn){

  gStyle->SetOptFit(kTRUE);
  gStyle->SetOptStat(kFALSE);

  Char_t name[200], title[200];
  Char_t myLabel[200];
  //==============================================
  //fitting phi and cosTheta
  //==============================================
  TCanvas *c1DFits[kNbRapForPTBins][kNbPTBins+1];
  Int_t fitResult;

  Int_t rapMin = rap, rapMax = rap;
  if(rap == 0){ rapMin = 1; rapMax = kNbRapForPTBins;}

  Int_t pTMin = pT, pTMax = pT+1;
  if(pT == 0){ pTMin = 1; pTMax = kNbPTBins;}

  //preparing an output file:
  if(rap == 0 && pT == 0)
    sprintf(name, "Results/fitResults_%s_%s.txt", frameLabel[iFrame], label);
  else
    sprintf(name, "Results/fitResults_%s_rap%d_pT%d_%s.txt", frameLabel[iFrame], rap, pT, label);
  FILE *fResults = fopen(name, "write");
  printf("opened text file %s to write out the results\n", name);

  TLine *lineP = new TLine(0., 0., 0., hData1D_pol_pT_rap[iFrame][1][1]->GetMaximum());
  lineP->SetLineStyle(3);
  lineP->SetLineColor(4);
  Int_t nBCosTh = hData2D_pol_pT_rap[iFrame][1][1]->GetNbinsX();
  Int_t nBPhi = hData2D_pol_pT_rap[iFrame][1][1]->GetNbinsY();
  Int_t nB = nBCosTh * nBPhi;

  //for(int iRap = 1; iRap <= kNbRapForPTBins; iRap++){
  for(int iRap = rapMin; iRap <= rapMax; iRap++){
    fprintf(fResults, "\n");
    //for(int iPTBin = 1; iPTBin < kNbPTBins; iPTBin++){//exclude last pTBin
    for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
      sprintf(name, "c1DFits_rap%d_pT%d_%s", iRap, iPTBin, frameLabel[iFrame]);    
      c1DFits[iRap][iPTBin] = new TCanvas(name, name, 1200, 700);

      hData1D_pol_pT_rap[iFrame][iPTBin][iRap]->SetMinimum(0.);
      hData1D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("colz");
//        hData1D_pol_pT_rap[iFrame][iPTBin][iRap]->Draw("lego2");

      //fit function:
      sprintf(name, "fit1D_rap%d_pT%d_%s", iRap, iPTBin, frameLabel[iFrame]);
      f1D_rap_pT[iFrame][iRap][iPTBin] = new TF1(name, fit1D, 0., 360., 4);
      f1D_rap_pT[iFrame][iRap][iPTBin]->SetParNames("norm", "lambda_theta", "lambda_phi", "lambda_thetaPhi");
      f1D_rap_pT[iFrame][iRap][iPTBin]->SetParameters(0.2, 0.0, 0.0, 0.0); 
      // f1D_rap_pT[iFrame][iRap][iPTBin]->SetParameters(200., 0.0, 0.0, 0.0); 
//       f1D_rap_pT[iFrame][iRap][iPTBin]->SetParLimits(1, -1.0, 1.0);
//       f1D_rap_pT[iFrame][iRap][iPTBin]->SetParLimits(2, -1.0, 1.0); 
      f1D_rap_pT[iFrame][iRap][iPTBin]->FixParameter(3, 0.); 

      printf("\n\n\nfitting histogram for %s, rap %d, pT %d\n", frameLabel[iFrame], iRap, iPTBin);
      thisHist1D = hData1D_pol_pT_rap[iFrame][iPTBin][iRap];
      thisHist1D->Fit(name, "+");
//       fitResult = thisHist1D->Fit(name, "EM", "", -maxCosTh, maxCosTh);
      fitResult = thisHist1D->Fit(name, "E");

      if(iRap == 1){
      	sprintf(myLabel, "%s: |y| < %1.1f,   %1.1f < p_{T} < %1.1f GeV/c", 
		frameLabel[iFrame], rapForPTRange[iRap], pTRange[iPTBin-1], pTRange[iPTBin]);
      }
      else{
      	sprintf(myLabel, "%s: %1.1f < |y| < %1.1f,   %1.1f < p_{T} < %1.1f GeV/c", 
      		frameLabel[iFrame], rapForPTRange[iRap-1], rapForPTRange[iRap], 
		pTRange[iPTBin-1], pTRange[iPTBin]);
      }
      TLatex *tex1 = new TLatex(10., 0.95*hData1D_pol_pT_rap[iFrame][iPTBin][iRap]->GetMaximum(), myLabel);
      //TLatex *tex1 = new TLatex(10., 0.9*thisHist1D->GetMaximum(), myLabel);
      tex1->SetTextSize(0.05); tex1->Draw("same");


      if(fitResult != 0){
	printf("fit did not converge: %d\n", fitResult);
	fit1DConverged[iRap][iPTBin] = kFALSE;
      }
      else
	fit1DConverged[iRap][iPTBin] = kTRUE;

      if(fit1DConverged[iRap][iPTBin]){
 	f1D_rap_pT[iFrame][iRap][iPTBin] = (TF1 *) thisHist1D->GetFunction(name);
	chi2_1D_rap_pT[iRap][iPTBin] = f1D_rap_pT[iFrame][iRap][iPTBin]->GetChisquare();
	NDF_1D_rap_pT[iRap][iPTBin] = f1D_rap_pT[iFrame][iRap][iPTBin]->GetNDF();
	if(chi2_1D_rap_pT[iRap][iPTBin] > 1000000.)
	  continue;

	norm1D_rap_pT[iRap][iPTBin] = f1D_rap_pT[iFrame][iRap][iPTBin]->GetParameter(0);
	lambda1D_rap_pT[iRap][iPTBin][0] = f1D_rap_pT[iFrame][iRap][iPTBin]->GetParameter(1);
	lambda1D_rap_pT[iRap][iPTBin][1] = f1D_rap_pT[iFrame][iRap][iPTBin]->GetParameter(2);
	lambda1D_rap_pT[iRap][iPTBin][2] = f1D_rap_pT[iFrame][iRap][iPTBin]->GetParameter(3);
	errLambda1D_rap_pT[iRap][iPTBin][0] = f1D_rap_pT[iFrame][iRap][iPTBin]->GetParError(1);
	errLambda1D_rap_pT[iRap][iPTBin][1] = f1D_rap_pT[iFrame][iRap][iPTBin]->GetParError(2);
	errLambda1D_rap_pT[iRap][iPTBin][2] = f1D_rap_pT[iFrame][iRap][iPTBin]->GetParError(3);

	frameInv[iRap][iPTBin] = (1. + lambda1D_rap_pT[iRap][iPTBin][0] + 2.*lambda1D_rap_pT[iRap][iPTBin][1]) / (3.+lambda1D_rap_pT[iRap][iPTBin][0]);
	printf("frameInvariant is %f; lambda_theta = %1.3f, lambda_phi = %f\n", 
	       frameInv[iRap][iPTBin], lambda1D_rap_pT[iRap][iPTBin][0], lambda1D_rap_pT[iRap][iPTBin][1]);

	//       f1D_rap_pT[iFrame][iRap][iPTBin]->Draw("same");

	sprintf(name, "plotFit1D_rap%d_pT%d_%s", iRap, iPTBin, frameLabel[iFrame]);
	plotF1D_rap_pT[iFrame][iRap][iPTBin] = new TF1(name, fit1D, 0., 360., 4);
	plotF1D_rap_pT[iFrame][iRap][iPTBin]->FixParameter(0, norm1D_rap_pT[iRap][iPTBin]);
	plotF1D_rap_pT[iFrame][iRap][iPTBin]->FixParameter(1, lambda1D_rap_pT[iRap][iPTBin][0]);
	plotF1D_rap_pT[iFrame][iRap][iPTBin]->FixParameter(2, lambda1D_rap_pT[iRap][iPTBin][1]);
	plotF1D_rap_pT[iFrame][iRap][iPTBin]->FixParameter(3, lambda1D_rap_pT[iRap][iPTBin][2]);
  	// plotF1D_rap_pT[iFrame][iRap][iPTBin]->Draw("csame");

	for(int iL = 0; iL < nB; iL++)
	  lineP->DrawLine(nBCosTh*iL, 0., nBCosTh*iL, hData1D_pol_pT_rap[iFrame][iPTBin][iRap]->GetMaximum());

	if(displayOn){
 	  for(int iPhi = 0; iPhi < nBPhi; iPhi++){
	    printf("preparing 1D plot for phi bin %d\n", iPhi);
	    PlotFitResults1D(hData2D_pol_pT_rap[iFrame][iPTBin][iRap], plotF1D_rap_pT[iFrame][iRap][iPTBin], iFrame, iRap, iPTBin, iPhi, label);
	    printf("done\n");
	  }
// 	  PlotFitResults1D(hData2D_pol_pT_rap[iFrame][iPTBin][iRap], plotF1D_rap_pT[iFrame][iRap][iPTBin], iPTBin, iRap, 3);
	}

	//store the results in an output file:	
	if(chi2_1D_rap_pT[iRap][iPTBin] < 100000.){
	  printf("storing results in outputfile\n");
	  fprintf(fResults, "%s, rap %d, pT %d: theta=%1.3f +- %1.3f, phi=%1.3f +- %1.3f, thetaPhi=%1.3f +- %1.3f, F=%1.3f, chi2/ndf = %1.3f / %d = %1.3f\n", frameLabel[iFrame], iRap, iPTBin, lambda1D_rap_pT[iRap][iPTBin][0], errLambda1D_rap_pT[iRap][iPTBin][0], lambda1D_rap_pT[iRap][iPTBin][1], errLambda1D_rap_pT[iRap][iPTBin][1], lambda1D_rap_pT[iRap][iPTBin][2], errLambda1D_rap_pT[iRap][iPTBin][2], frameInv[iRap][iPTBin], chi2_1D_rap_pT[iRap][iPTBin], NDF_1D_rap_pT[iRap][iPTBin], chi2_1D_rap_pT[iRap][iPTBin]/NDF_1D_rap_pT[iRap][iPTBin]);
	  printf("done\n");
	}
	c1DFits[iRap][iPTBin]->Update();
	sprintf(name, "Figures/FitResults1D_%s/fits1D_%s_rap%d_pT%d.gif", label, frameLabel[iFrame], iRap, iPTBin);
	c1DFits[iRap][iPTBin]->Print(name);
	sprintf(name, "Figures/FitResults1D_%s/fits1D_%s_rap%d_pT%d.pdf", label, frameLabel[iFrame], iRap, iPTBin);
	c1DFits[iRap][iPTBin]->Print(name);
	sprintf(name, "Figures/FitResults1D_%s/fits1D_%s_rap%d_pT%d.eps", label, frameLabel[iFrame], iRap, iPTBin);
	c1DFits[iRap][iPTBin]->Print(name);
      }
    }

    //calculate the global chi2:
    chi2_1D_rap_global[iRap] = 0.; NDF_1D_rap_global[iRap] = 0;
    for(int iPTBin = 1; iPTBin < kNbPTBins; iPTBin++){//exclude last PTBin
      if(fit1DConverged[iRap][iPTBin]){
	chi2_1D_rap_global[iRap] += chi2_1D_rap_pT[iRap][iPTBin];
	NDF_1D_rap_global[iRap] += NDF_1D_rap_pT[iRap][iPTBin];
      }
    }
    probChi2_1D_rap_global = TMath::Prob(chi2_1D_rap_global[iRap], NDF_1D_rap_global[iRap]);
    chi2_1D_rap_global[iRap] /= NDF_1D_rap_global[iRap];
    printf("global (norm.) chi2 = %1.3f and probability %1.3f\n", chi2_1D_rap_global[iRap], probChi2_1D_rap_global);
    //
    Double_t lambda_pT_forPlot[3][kNbPTBins], errLambda_pT_forPlot[3][kNbPTBins];
    Double_t errL[kNbPTBins], errR[kNbPTBins];
    Double_t frameInv_pT_forPlot[kNbPTBins];

    for(int iPTBin = 1; iPTBin < kNbPTBins; iPTBin++){//exclude last PTBin
      for(int iPar = 0; iPar < 3; iPar++){
	if(fit1DConverged[iRap][iPTBin] && chi2_1D_rap_pT[iRap][iPTBin] < 100000.){
	  lambda_pT_forPlot[iPar][iPTBin-1] = lambda1D_rap_pT[iRap][iPTBin][iPar];
	  errLambda_pT_forPlot[iPar][iPTBin-1] = errLambda1D_rap_pT[iRap][iPTBin][iPar];
	  errL[iPTBin-1] = pTWCentre_rap[iRap-1][iPTBin-1] - pTRange[iPTBin-1];
	  errR[iPTBin-1] = pTRange[iPTBin] - pTWCentre_rap[iRap-1][iPTBin-1];
	  frameInv_pT_forPlot[iPTBin-1] = frameInv[iRap][iPTBin];
	}
	else{
	  lambda_pT_forPlot[iPar][iPTBin-1] = 1000.;
	  errLambda_pT_forPlot[iPar][iPTBin-1] = 0.;
	  errL[iPTBin-1] = pTWCentre_rap[iRap-1][iPTBin-1] - pTRange[iPTBin-1];
	  errR[iPTBin-1] = pTRange[iPTBin] - pTWCentre_rap[iRap-1][iPTBin-1];
	  frameInv_pT_forPlot[iPTBin-1] = 1000.;
	}
      }
    }
    //TGraph for lambdaTheta
    gLambdaTh1D_rap[iFrame][iRap] = new TGraphAsymmErrors(kNbPTBins-1, pTWCentre_rap[iRap-1], lambda_pT_forPlot[0], errL, errR, errLambda_pT_forPlot[0], errLambda_pT_forPlot[0]);
    sprintf(name, "gLambdaTh1D_%s_rap%d", frameLabel[iFrame], iRap);
    gLambdaTh1D_rap[iFrame][iRap]->SetName(name);
    gLambdaTh1D_rap[iFrame][iRap]->SetMarkerStyle(marker_rapForPTBins[iRap]);
    gLambdaTh1D_rap[iFrame][iRap]->SetMarkerColor(colour_rapForPTBins[iRap]);
    gLambdaTh1D_rap[iFrame][iRap]->SetLineColor(colour_rapForPTBins[iRap]);

    //TGraph for lambdaPhi
    gLambdaPhi1D_rap[iFrame][iRap] = new TGraphAsymmErrors(kNbPTBins-1, pTWCentre_rap[iRap-1], lambda_pT_forPlot[1], errL, errR, errLambda_pT_forPlot[1], errLambda_pT_forPlot[1]);
    sprintf(name, "gLambdaPhi1D_%s_rap%d", frameLabel[iFrame], iRap);
    gLambdaPhi1D_rap[iFrame][iRap]->SetName(name);
    gLambdaPhi1D_rap[iFrame][iRap]->SetMarkerStyle(marker_rapForPTBins[iRap]);
    gLambdaPhi1D_rap[iFrame][iRap]->SetMarkerColor(colour_rapForPTBins[iRap]);
    gLambdaPhi1D_rap[iFrame][iRap]->SetLineColor(colour_rapForPTBins[iRap]);

    //TGraph for lambdaThetaPhi
    gLambdaThPhi1D_rap[iFrame][iRap] = new TGraphAsymmErrors(kNbPTBins-1, pTWCentre_rap[iRap-1], lambda_pT_forPlot[2], errL, errR, errLambda_pT_forPlot[2], errLambda_pT_forPlot[2]);
    sprintf(name, "gLambdaThPhi1D_%s_rap%d", frameLabel[iFrame], iRap);
    gLambdaThPhi1D_rap[iFrame][iRap]->SetName(name);
    gLambdaThPhi1D_rap[iFrame][iRap]->SetMarkerStyle(marker_rapForPTBins[iRap]);
    gLambdaThPhi1D_rap[iFrame][iRap]->SetMarkerColor(colour_rapForPTBins[iRap]);
    gLambdaThPhi1D_rap[iFrame][iRap]->SetLineColor(colour_rapForPTBins[iRap]);

    //TGraph for frameInvariant F
    gFrameInv1D_rap[iFrame][iRap] = new TGraph(kNbPTBins-1, pTWCentre_rap[iRap-1], frameInv_pT_forPlot);
    sprintf(name, "gFrameInv1D_%s_rap%d", frameLabel[iFrame], iRap);
    gFrameInv1D_rap[iFrame][iRap]->SetName(name);
    gFrameInv1D_rap[iFrame][iRap]->SetMarkerStyle(marker_rapForPTBins[iRap]);
    gFrameInv1D_rap[iFrame][iRap]->SetMarkerColor(colour_rapForPTBins[iRap]);
    gFrameInv1D_rap[iFrame][iRap]->SetLineColor(colour_rapForPTBins[iRap]);
  }

  fclose(fResults);
  printf("after closing output file\n");

  if(rap != 0 || pT != 0) return;

  for(int iRap = 1; iRap <= kNbRapForPTBins; iRap++){
//    for(int iRap = 1; iRap <= 2; iRap++){
    for(int iPT = 1; iPT < kNbPTBins; iPT++){//exclude last pT bin

      if(fit1DConverged[iRap][iPT]){
	printf("%10f %10f %10f, F = %10f, chi2 %1.3f, ndf = %1.3f, chi2/ndf = %1.3f\n", 
	       lambda1D_rap_pT[iRap][iPT][0],
	       lambda1D_rap_pT[iRap][iPT][1],
	       lambda1D_rap_pT[iRap][iPT][2],
	       frameInv[iRap][iPT],
	       chi2_1D_rap_pT[iRap][iPT],
	       NDF_1D_rap_pT[iRap][iPT],
	       chi2_1D_rap_pT[iRap][iPT] / NDF_1D_rap_pT[iRap][iPT]);
      }
    }
  }

  //============================================
  //lambda_theta vs pT  (3 rap bins)
  //============================================
  sprintf(name, "cLambdaTh1D_%s", frameLabel[iFrame]);
  TCanvas *cLambdaTh1D = new TCanvas(name, name);
  TH1F *hFrameLTh1D = gPad->DrawFrame(0., -1.2, 30., 1.2);
  hFrameLTh1D->SetXTitle("p_{T} [GeV/c]");
  sprintf(name, "#lambda_{#theta}^{%s}", frameLabel[iFrame]);
  hFrameLTh1D->SetYTitle(name);

  gLambdaTh1D_rap[iFrame][1]->Draw("psame");
  gLambdaTh1D_rap[iFrame][2]->Draw("psame");
  gLambdaTh1D_rap[iFrame][3]->Draw("psame");

  TLine *line = new TLine(0., 0., 30., 0.);
  line->SetLineStyle(3);  line->Draw();
  line->DrawLine(0., 1., 30., 1.);
  line->DrawLine(0., -1., 30., -1.);

  sprintf(name, "Figures/lambdaTh1D_%s_%s.eps", frameLabel[iFrame], label);  cLambdaTh1D->Print(name);
  sprintf(name, "Figures/lambdaTh1D_%s_%s.pdf", frameLabel[iFrame], label);  cLambdaTh1D->Print(name);
  sprintf(name, "Figures/lambdaTh1D_%s_%s.gif", frameLabel[iFrame], label);  cLambdaTh1D->Print(name);

  //============================================
  //lambda_phi vs pT (3 rap bins)
  //============================================
  sprintf(name, "cLambdaPhi1D_%s", frameLabel[iFrame]);
  TCanvas *cLambdaPhi1D_pT = new TCanvas(name, name);
  TH1F *hFrameLPhi1D_pT = gPad->DrawFrame(0., -1.2, 30., 1.4);
  hFrameLPhi1D_pT->SetXTitle("p_{T} [GeV/c]");
  sprintf(name, "#lambda_{#phi}^{%s}", frameLabel[iFrame]);
  hFrameLPhi1D_pT->SetYTitle(name);

  gLambdaPhi1D_rap[iFrame][1]->Draw("psame");
  gLambdaPhi1D_rap[iFrame][2]->Draw("psame");
  gLambdaPhi1D_rap[iFrame][3]->Draw("psame");

  TLine *line2 = new TLine(0., 0., 30., 0.);
  line2->SetLineStyle(3);  line2->Draw();
  line2->DrawLine(0., 0.33333, 30., 0.333333);
  line2->DrawLine(0., -1., 30., -1.);

  sprintf(name, "Figures/lambdaPhi1D_%s_%s.eps", frameLabel[iFrame], label);  cLambdaPhi1D_pT->Print(name);
  sprintf(name, "Figures/lambdaPhi1D_%s_%s.pdf", frameLabel[iFrame], label);  cLambdaPhi1D_pT->Print(name);
  sprintf(name, "Figures/lambdaPhi1D_%s_%s.gif", frameLabel[iFrame], label);  cLambdaPhi1D_pT->Print(name);

  //============================================
  //lambda_ThetaPhi vs pT (3 rap bins)
  //============================================
  sprintf(name, "cLambdaThPhi1D_%s", frameLabel[iFrame]);
  TCanvas *cLambdaThPhi1D_pT = new TCanvas(name, name);
  TH1F *hFrameLThPhi1D_pT = gPad->DrawFrame(0., -1.4, 30., 1.4);
  hFrameLThPhi1D_pT->SetXTitle("p_{T} [GeV/c]");
  sprintf(name, "#lambda_{#theta#phi}^{%s}", frameLabel[iFrame]);
  hFrameLThPhi1D_pT->SetYTitle(name);

  gLambdaThPhi1D_rap[iFrame][1]->Draw("psame");
  gLambdaThPhi1D_rap[iFrame][2]->Draw("psame");
  gLambdaThPhi1D_rap[iFrame][3]->Draw("psame");

  TLine *line3 = new TLine(0., 0., 30., 0.);
  line3->SetLineStyle(3);  line3->Draw();
  line3->DrawLine(0., 0.33333, 30., 0.333333);
  line3->DrawLine(0., -1., 30., -1.);

  sprintf(name, "Figures/lambdaThPhi1D_%s_%s.eps", frameLabel[iFrame], label);  cLambdaThPhi1D_pT->Print(name);
  sprintf(name, "Figures/lambdaThPhi1D_%s_%s.pdf", frameLabel[iFrame], label);  cLambdaThPhi1D_pT->Print(name);
  sprintf(name, "Figures/lambdaThPhi1D_%s_%s.gif", frameLabel[iFrame], label);  cLambdaThPhi1D_pT->Print(name);

  //============================================
  //frame invariant F
  //============================================
  sprintf(name, "cFrameInvF1D_%s", frameLabel[iFrame]);
  TCanvas *cFrameInv1D_pT = new TCanvas(name, name);
  TH1F *hFrameInvF1D_pT = gPad->DrawFrame(0., -1.2, 30., 1.2);
  hFrameInvF1D_pT->SetXTitle("p_{T} [GeV/c]");
  sprintf(name, "frame Invariant F for %s", frameLabel[iFrame]);
  hFrameInvF1D_pT->SetYTitle(name);

  gFrameInv1D_rap[iFrame][1]->Draw("psame");
  gFrameInv1D_rap[iFrame][2]->Draw("psame");
  gFrameInv1D_rap[iFrame][3]->Draw("psame");

  TLine *line4 = new TLine(0., 0., 30., 0.);
  line4->SetLineStyle(3);  line4->Draw();
  line4->DrawLine(0., 1., 30., 1.);
  line4->DrawLine(0., -1., 30., -1.);

  sprintf(name, "Figures/frameInvF1D_%s_%s.eps", frameLabel[iFrame], label);  cFrameInv1D_pT->Print(name);
  sprintf(name, "Figures/frameInvF1D_%s_%s.pdf", frameLabel[iFrame], label);  cFrameInv1D_pT->Print(name);
  sprintf(name, "Figures/frameInvF1D_%s_%s.gif", frameLabel[iFrame], label);  cFrameInv1D_pT->Print(name);
}


//======================================
void Get1DHistoFrom2D(){

  Char_t name[200];
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
	    acc = hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iBinX,iBinY) / hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetMaximum();

// 	    if(contentReco < minBinContent || acc < fracMaxAcc*maxAcc){
 	    if(contentReco < minBinContent || acc < minAcc){
	      // if(iPTBin == 4 && iRapBin == 1){
	      // 	printf("bin %d and %d: content in RECO histo %f, acc %f --> setting bin content of %d to -1\n", 
	      // 	       iBinY, iBinX, contentReco, acc, binID);
	      // }
	      content = -1.;
	      error = 10000.;
	      hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(binID, content);
	      hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinError(binID, error);
	    }
	    else{
	      // if(iPTBin == 4 && iRapBin == 1){
	      // 	printf("bin %d and %d: content in RECO histo %f, acc %f --> setting bin content of %d to %f\n", 
	      // 	       iBinY, iBinX, contentReco, acc, binID, content);
	      // }

	      hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(binID, content);
	      hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinError(binID, error);
	    }
	  }
	}
      }
    }
  }
}

//======================================
void PlotFitResults1D(TH2D *hHist, TF1 *fFit, Int_t iFrame, Int_t thisRap, Int_t thisPT, Int_t binPhi, Char_t *label){

  printf("<PlotFitResults1D>: begining\n");
  if(binPhi > hHist->GetNbinsY() || binPhi < 0){
    printf("<PlotFitResults1D> rejecting phi bin %d\n", binPhi);
    return;
  }
  Char_t name[300], title[300];
  //project the histos into 1D histos:
  TH1D *hCosTh[nBinsCosT], *hPhi[nBinsPhiPol]; //maximum bin size if not rebinned
  for(int iBinX = 0; iBinX < hHist->GetNbinsX(); iBinX++){
    sprintf(name, "%s_phi%d", hHist->GetName(), iBinX);
    hPhi[iBinX] = (TH1D *) hHist->ProjectionY(name, iBinX+1, iBinX+1);
    // hPhi[iBinX]->SetLineColor(iBinX);
  }
  for(int iBinY = 0; iBinY < hHist->GetNbinsY(); iBinY++){
    sprintf(name, "%s_cosTh%d", hHist->GetName(), iBinY);
    hCosTh[iBinY] = (TH1D *) hHist->ProjectionX(name, iBinY+1, iBinY+1);
    // hCosTh[iBinY]->SetLineColor(iBinY);
  }

  Double_t phiMin = hHist->GetYaxis()->GetBinLowEdge(binPhi+1);
  Double_t phiMax = hHist->GetYaxis()->GetBinLowEdge(binPhi+2);
  Double_t deltaPhi = (phiMax - phiMin);
  // printf("<PlotFitResults1D> phi slice %d: want to loop from %1.2f to %1.2f in phi\n", binPhi, phiMin, phiMax);

  //store the function values for this phi slice in a histogram:
  sprintf(name, "hPlot_%s_rap%d_pT%d_phi%d", frameLabel[iFrame], thisRap, thisPT, binPhi);
  TH1F *hPlot = new TH1F(name, "", nBinsCosT, -1., 1.);
  hPlot->SetLineWidth(2);
  Double_t content;
  Double_t binCentre, myX;
  for(int iBin = 1; iBin <= nBinsCosT; iBin++){
    binCentre = hPlot->GetBinCenter(iBin);
    myX = (binCentre + 1)*10. + binPhi*deltaPhi;
    content = fFit->Eval(myX);
    // printf("<PlotFitResults1D> bin %d corresponds to phi = %1.2f and function has value %f\n", iBin, myX, content);
    hPlot->SetBinContent(iBin, content);
  }
  //plot only were the data sit:
  Double_t cosThMinPlot = -1., cosThMaxPlot = 1.;
  Bool_t minFound = kFALSE;
  for(int iBin = 1; iBin <= hCosTh[binPhi]->GetNbinsX(); iBin++){
    content = hCosTh[binPhi]->GetBinContent(iBin);
    if(content > 0){
      if(!minFound){
	cosThMinPlot = hCosTh[binPhi]->GetBinLowEdge(iBin);
	minFound = kTRUE;
      }
      cosThMaxPlot = hCosTh[binPhi]->GetBinLowEdge(iBin+1);
    }
  }
  printf("phi %d: function will be drawn from cosTheta %1.3f to %1.3f\n", 
	 binPhi, cosThMinPlot, cosThMaxPlot);
  Double_t eps = 0.005;
  hPlot->SetAxisRange(cosThMinPlot+eps, cosThMaxPlot-eps);

  // if(hCosTh[binPhi]->Integral() < 1){
  //   printf("<PlotFitResults1D> integral is %f\n", hCosTh[binPhi]->Integral());
  //   return;
  // }

  sprintf(name, "canv_%s_%d", hHist->GetName(), binPhi);
  TCanvas *cCosTh = new TCanvas(name, "", 700, 500);
  TH1F *hFrameCosTh = gPad->DrawFrame(-1., 0., 1., 1.5*hCosTh[binPhi]->GetMaximum());
  // TH1F *hFrameCosTh = gPad->DrawFrame(-1., 0., 1., 10000.);
  hFrameCosTh->SetXTitle(hHist->GetXaxis()->GetTitle());
  Double_t xmin, xmax, ymin, ymax;

  hCosTh[binPhi]->SetLineColor(1);
  hCosTh[binPhi]->SetMarkerColor(1);
  if(binPhi < hHist->GetNbinsY() / 2){
    hCosTh[binPhi]->SetMarkerStyle(20);
    hCosTh[binPhi]->SetMarkerColor(2);
  }
  else{
    hCosTh[binPhi]->SetMarkerStyle(25);
    hCosTh[binPhi]->SetMarkerColor(4);
  }
  hCosTh[binPhi]->Draw("psame");

  hPlot->Draw("csame");

  printf("after drawing\n");

  if(thisRap == 1){
    sprintf(name, "|y| < %1.1f,   %1.1f < p_{T} < %1.1f GeV/c,   %1.0f < #phi < %1.0f", 
	    rapForPTRange[thisRap], pTRange[thisPT-1], pTRange[thisPT], phiMin, phiMax);
  }
  else{
    sprintf(name, "%1.1f < |y| < %1.1f,   %1.1f < p_{T} < %1.1f GeV/c,   %1.0f < #phi < %1.0f", 
	    rapForPTRange[thisRap-1], rapForPTRange[thisRap], pTRange[thisPT-1], pTRange[thisPT],
	    phiMin, phiMax);
  }
  TLatex *phiL = new TLatex(-0.93, 1.4*hCosTh[binPhi]->GetMaximum(), name);
  phiL->SetTextSize(0.05); phiL->Draw();

  sprintf(name, "Figures/FitResults1D_%s/fits1D_1DRep_%s_rapBin%d_pTBin%d_phiBin%d.eps", label, frameLabel[iFrame], thisRap, thisPT, binPhi);  cCosTh->Print(name);
  sprintf(name, "Figures/FitResults1D_%s/fits1D_1DRep_%s_rapBin%d_pTBin%d_phiBin%d.pdf", label, frameLabel[iFrame], thisRap, thisPT, binPhi);  cCosTh->Print(name);
  sprintf(name, "Figures/FitResults1D_%s/fits1D_1DRep_%s_rapBin%d_pTBin%d_phiBin%d.gif", label, frameLabel[iFrame], thisRap, thisPT, binPhi);  cCosTh->Print(name);
  printf("<PlotFitResults1D>: end\n");
}

//======================================
void PlotUncorrData(Int_t iFrame, Char_t *label){

  gStyle->SetOptStat(10);
  Char_t name[200], title[200];
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
void PlotUncorrData2D(Int_t iFrame, Char_t *label){

  Char_t name[200], title[200];
  TCanvas *c2D[kNbRapForPTBins+1];
  for(int iRap = 1; iRap <= kNbRapForPTBins; iRap++){
    sprintf(name, "c2D_%s_rap%d", frameLabel[iFrame], iRap);
    sprintf(title, "phi vs cosTheta for pT bins, rap = %d (%s)", iRap, frameLabel[iFrame]);
    c2D[iRap] = new TCanvas(name, title, 1000, 700);
    c2D[iRap]->Divide(3,2);
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      c2D[iRap]->cd(iPTBin);
      if(iPTBin == kNbPTBins) 
	sprintf(name, "#Upsilon(1S): %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", rapForPTRange[iRap-1], rapForPTRange[iRap], pTRange[iPTBin-1]);
      else 
	sprintf(name, "#Upsilon(1S): %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", rapForPTRange[iRap-1], rapForPTRange[iRap], pTRange[iPTBin-1], pTRange[iPTBin]);
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
void PlotCorrData2D(Int_t iFrame, Char_t *label){

  gStyle->SetOptStat(0);
  Char_t name[200], title[200];
  TCanvas *c2D[kNbRapForPTBins+1];
  for(int iRap = 1; iRap <= kNbRapForPTBins; iRap++){
    sprintf(name, "c2DAfter_%s_rap%d", frameLabel[iFrame], iRap);
    sprintf(title, "Acc corrected phi vs cosTheta for pT bins, rap = %d (%s)", iRap, frameLabel[iFrame]);
    c2D[iRap] = new TCanvas(name, title, 1000, 700);
    c2D[iRap]->Divide(3,2);
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      c2D[iRap]->cd(iPTBin);
      if(iPTBin == kNbPTBins) 
	sprintf(name, "#Upsilon(1S): %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", rapForPTRange[iRap-1], rapForPTRange[iRap], pTRange[iPTBin-1]);
      else 
	sprintf(name, "#Upsilon(1S): %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", rapForPTRange[iRap-1], rapForPTRange[iRap], pTRange[iPTBin-1], pTRange[iPTBin]);
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

  printf("\n\ncorrecting for acceptance\n\n\n");

  Char_t name[200];
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    // for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){

    //   sprintf(name, "hData_Onia_cosTh_%s_pT%d", frameLabel[iFrame], iPTBin);
    //   hData_pol_pT[iFrame][iPTBin][cosThPol] = (TH1D *) Reco_pol_pT[iFrame][iPTBin][cosThPol]->Clone(name);
    //   hData_pol_pT[iFrame][iPTBin][cosThPol]->Divide(hAcc_pol_pT[iFrame][iPTBin][cosThPol]);
    //   hData_pol_pT[iFrame][iPTBin][cosThPol]->SetLineColor(colour_pT[iPTBin]);
    //   hData_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerColor(colour_pT[iPTBin]);
    //   hData_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerStyle(marker_pT[iPTBin]);

    //   sprintf(name, "hData_Onia_phi_%s_pT%d", frameLabel[iFrame], iPTBin);
    //   hData_pol_pT[iFrame][iPTBin][phiPol] = (TH1D *) Reco_pol_pT[iFrame][iPTBin][phiPol]->Clone(name);
    //   hData_pol_pT[iFrame][iPTBin][phiPol]->Divide(hAcc_pol_pT[iFrame][iPTBin][phiPol]);
    //   hData_pol_pT[iFrame][iPTBin][phiPol]->SetLineColor(colour_pT[iPTBin]);
    //   hData_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerColor(colour_pT[iPTBin]);
    //   hData_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerStyle(marker_pT[iPTBin]);

    //   sprintf(name, "hData2D_Onia_%s_pT%d", frameLabel[iFrame], iPTBin);
    //   hData2D_pol_pT[iFrame][iPTBin] = (TH2D *) Reco2D_pol_pT[iFrame][iPTBin]->Clone(name);
    //   hData2D_pol_pT[iFrame][iPTBin]->Divide(hAcc2D_pol_pT[iFrame][iPTBin]);
    //   hData2D_pol_pT[iFrame][iPTBin]->SetLineColor(colour_pT[iPTBin]);
    //   hData2D_pol_pT[iFrame][iPTBin]->SetMarkerColor(colour_pT[iPTBin]);
    //   hData2D_pol_pT[iFrame][iPTBin]->SetMarkerStyle(marker_pT[iPTBin]);
    // }
    // for(int iRapBin = 0; iRapBin < 2*kNbRapBins+1; iRapBin++) {
    //   sprintf(name, "hData_Onia_cosTh_%s_rap%d", frameLabel[iFrame], iRapBin);
    //   hData_pol_rap[iFrame][iRapBin][cosThPol] = (TH1D *) Reco_pol_rap[iFrame][iRapBin][cosThPol]->Clone(name);
    //   hData_pol_rap[iFrame][iRapBin][cosThPol]->Divide(hAcc_pol_rap[iFrame][iRapBin][cosThPol]);
    //   hData_pol_rap[iFrame][iRapBin][cosThPol]->SetLineColor(colour_rap[iRapBin]);
    //   hData_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerColor(colour_rap[iRapBin]);
    //   hData_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerStyle(marker_rap[iRapBin]);

    //   sprintf(name, "hData_Onia_phi_%s_rap%d", frameLabel[iFrame], iRapBin);
    //   hData_pol_rap[iFrame][iRapBin][phiPol] = (TH1D *) Reco_pol_rap[iFrame][iRapBin][phiPol]->Clone(name);
    //   hData_pol_rap[iFrame][iRapBin][phiPol]->Divide(hAcc_pol_rap[iFrame][iRapBin][phiPol]);
    //   hData_pol_rap[iFrame][iRapBin][phiPol]->SetLineColor(colour_rap[iRapBin]);
    //   hData_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerColor(colour_rap[iRapBin]);
    //   hData_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerStyle(marker_rap[iRapBin]);

    //   sprintf(name, "hData2D_Onia_%s_rap%d", frameLabel[iFrame], iRapBin);
    //   hData2D_pol_rap[iFrame][iRapBin] = (TH2D *) Reco2D_pol_rap[iFrame][iRapBin]->Clone(name);
    //   hData2D_pol_rap[iFrame][iRapBin]->Divide(hAcc2D_pol_rap[iFrame][iRapBin]);
    //   hData2D_pol_rap[iFrame][iRapBin]->SetLineColor(colour_rap[iRapBin]);
    //   hData2D_pol_rap[iFrame][iRapBin]->SetMarkerColor(colour_rap[iRapBin]);
    //   hData2D_pol_rap[iFrame][iRapBin]->SetMarkerStyle(marker_rap[iRapBin]);
    // }
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
	// sprintf(name, "hData_Onia_cosTh_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol] = (TH1D *) Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->Clone(name);
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->Divide(hAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]);
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetLineColor(colour_pT[iPTBin]);
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerColor(colour_pT[iPTBin]);
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerStyle(marker_pT[iPTBin]);

	// sprintf(name, "hData_Onia_phi_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol] = (TH1D *) Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->Clone(name);
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->Divide(hAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]);
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetLineColor(colour_pT[iPTBin]);
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerColor(colour_pT[iPTBin]);
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerStyle(marker_pT[iPTBin]);

	sprintf(name, "hData2D_Onia_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Clone(name);
	// printf("data have %d cosTheta and %d phi bins\n", 
	//        hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(), 
	//        hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY());
	// printf("while acceptance histo has %d and %d\n", 
	//        hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(), 
	//        hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY());
	if(iPTBin == 4 && iRapBin == 1){
	  printf("bin 6 and 7 of RECO is %f\n", Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(6,7));
	}
 	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Divide(hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]);//H: will be done below
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetLineColor(colour_pT[iPTBin]);
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerColor(colour_pT[iPTBin]);
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerStyle(marker_pT[iPTBin]);
      }
    }
  }

  Double_t binContentReco, binErrorReco;
  Double_t acc, binContentData, binErrorData;
  Double_t content, error, maxAcc;

  // //do NOT introduce any error on the "data" due to the acceptance correction:
  // for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
  //   for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
  //     for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
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

  printf("before Get1DHistoFrom2D\n");
  Get1DHistoFrom2D();
  printf("after Get1DHistoFrom2D\n");
  SaveCorrData(label); //save the histos here BEFORE we set bin contents to zero for further processing!
  printf("after SaveCorrData\n");
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){

	maxAcc = hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetMaximum();
	// if(iPTBin == 4 && iRapBin == 1){
	//   printf("pT %d, rap %d, data have %d and %d bins\n", iPTBin, iRapBin, 
	// 	 hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(), 
	// 	 hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY());
	// }
	for(int iBinX = 1; iBinX <= hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(); iBinX++){
	  for(int iBinY = 1; iBinY <= hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY(); iBinY++){

	    //in case the bin content of the reco histo was below "minBinContent"
	    //set the value in the data histogram to -1
	    content = Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iBinX, iBinY);
	    if(content < minBinContent){
	      // if(iPTBin == 4 && iRapBin == 1){
	      // 	printf("reco2D has %f entries in bin %d and %d --> setting content to -1\n", content, iBinX, iBinY);
	      // }
	      content = -1.;
	      error = 10000.;
	      hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(iBinX, iBinY, content);
	      hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinError(iBinX, iBinY, error);
	    }

	    //in case the acceptance  was less than "fracMaxAcc"s lower than the maximum
	    //set the value in the data histogram to -1
	    acc = hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iBinX, iBinY) / 
	      hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetMaximum();
// 	    if(acc < fracMaxAcc*maxAcc){
 	    if(acc < minAcc){
	      content = -1.;
	      error = 10000.;
	      hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(iBinX, iBinY, content);
	      hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinError(iBinX, iBinY, error);
	    }
	  }
	}
      }
    }
  }
}

//======================================
void ReadData(Char_t *fileNameData, Int_t rebinCosTh, Int_t rebinPhi){

  TFile *fIn = new TFile(fileNameData);
  
  Char_t name[200];
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
      sprintf(name, "Reco_Onia_cosTh_%s_pT%d", frameLabel[iFrame], iPTBin);
      Reco_pol_pT[iFrame][iPTBin][cosThPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_pT[iFrame][iPTBin][cosThPol]->SetLineColor(colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerColor(colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerStyle(marker_pT[iPTBin]);
      // printf("frame %s, pT bin %d, histo %p\n", frameLabel[iFrame], iPTBin, Reco_pol_pT[iFrame][iPTBin][cosThPol]);
      //
      sprintf(name, "Reco_Onia_phi_%s_pT%d", frameLabel[iFrame], iPTBin);
      Reco_pol_pT[iFrame][iPTBin][phiPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_pT[iFrame][iPTBin][phiPol]->SetLineColor(colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerColor(colour_pT[iPTBin]);
      Reco_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerStyle(marker_pT[iPTBin]);
      // printf("frame %s, pT bin %d, histo %p\n", frameLabel[iFrame], iPTBin, Reco_pol_pT[iFrame][iPTBin][phiPol]);
      //
      sprintf(name, "Reco2D_Onia_%s_pT%d", frameLabel[iFrame], iPTBin);
      Reco2D_pol_pT[iFrame][iPTBin] = (TH2D *) gDirectory->Get(name);
      Reco2D_pol_pT[iFrame][iPTBin]->SetLineColor(colour_pT[iPTBin]);
      Reco2D_pol_pT[iFrame][iPTBin]->SetMarkerColor(colour_pT[iPTBin]);
      Reco2D_pol_pT[iFrame][iPTBin]->SetMarkerStyle(marker_pT[iPTBin]);
      // printf("frame %s, pT bin %d, histo %p\n", frameLabel[iFrame], iPTBin, Reco2D_pol_pT[iFrame][iPTBin]);
      if(rebinCosTh > 1 || rebinPhi > 1)
	Reco2D_pol_pT[iFrame][iPTBin]->Rebin2D(rebinCosTh, rebinPhi);
    }
    for(int iRapBin = 0; iRapBin < 2*kNbRapBins+1; iRapBin++){
      sprintf(name, "Reco_Onia_cosTh_%s_rap%d", frameLabel[iFrame], iRapBin);
      Reco_pol_rap[iFrame][iRapBin][cosThPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_rap[iFrame][iRapBin][cosThPol]->SetLineColor(colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerColor(colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerStyle(marker_rap[iRapBin]);
      // printf("frame %s, rap bin %d, histo %p\n", frameLabel[iFrame], iRapBin, Reco_pol_rap[iFrame][iRapBin][cosThPol]);
      //
      sprintf(name, "Reco_Onia_phi_%s_rap%d", frameLabel[iFrame], iRapBin);
      Reco_pol_rap[iFrame][iRapBin][phiPol] = (TH1D *) gDirectory->Get(name);
      Reco_pol_rap[iFrame][iRapBin][phiPol]->SetLineColor(colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerColor(colour_rap[iRapBin]);
      Reco_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerStyle(marker_rap[iRapBin]);
      // printf("frame %s, rap bin %d, histo %p\n", frameLabel[iFrame], iRapBin, Reco_pol_rap[iFrame][iRapBin][phiPol]);
      //
      sprintf(name, "Reco2D_Onia_%s_rap%d", frameLabel[iFrame], iRapBin);
      Reco2D_pol_rap[iFrame][iRapBin] = (TH2D *) gDirectory->Get(name);
      Reco2D_pol_rap[iFrame][iRapBin]->SetLineColor(colour_rap[iRapBin]);
      Reco2D_pol_rap[iFrame][iRapBin]->SetMarkerColor(colour_rap[iRapBin]);
      Reco2D_pol_rap[iFrame][iRapBin]->SetMarkerStyle(marker_rap[iRapBin]);
      // printf("frame %s, rap bin %d, histo %p\n", frameLabel[iFrame], iRapBin, Reco2D_pol_rap[iFrame][iRapBin]);
      if(rebinCosTh > 1 || rebinPhi > 1)
	Reco2D_pol_rap[iFrame][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
    }
    Int_t totEv = 0;
    printf("here\n");
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
	sprintf(name, "Reco_Onia_cosTh_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol] = (TH1D *) gDirectory->Get(name);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetLineColor(colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerColor(colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerStyle(marker_pT[iPTBin]);
	// printf("frame %s, pT %d, rap %d, histo %p\n", frameLabel[iFrame], iPTBin, iRapBin, Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]);
	//
	sprintf(name, "Reco_Onia_phi_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol] = (TH1D *) gDirectory->Get(name);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetLineColor(colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerColor(colour_pT[iPTBin]);
	Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerStyle(marker_pT[iPTBin]);
	// printf("frame %s, pT %d, rap %d, histo %p\n", frameLabel[iFrame], iPTBin, iRapBin, Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]);
	//
	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) gDirectory->Get(name);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetLineColor(colour_pT[iPTBin]);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerColor(colour_pT[iPTBin]);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerStyle(marker_pT[iPTBin]);
	// printf("frame %s, pT %d, rap %d, histo %p\n", frameLabel[iFrame], iPTBin, iRapBin, Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]);
	printf("%s: pTbin %d, rapBin %d --> statistics in data: %1.1f\n",
	       frameLabel[iFrame], iPTBin, iRapBin, Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Integral());

	totEv += Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Integral();
	if(rebinCosTh > 1 || rebinPhi > 1){
	  printf("rebinning 2D histo by %d and %d\n", rebinCosTh, rebinPhi);
	  Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
	  printf("#bins are now %d and %d\n",
		 Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(),
		 Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY());

	  if(iPTBin == 4 && iRapBin == 1){
	    for(int iBinX = 1; iBinX < Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(); iBinX++){
	      for(int iBinY = 1; iBinY < Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY(); iBinY++){
	  	printf("bin %d and %d --> binContent %f\n", iBinX, iBinY, 
	  	       Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iBinX, iBinY));
	      }
	    }
	  }
	}
      }
    }
    printf("%s: total nb. of events: %d\n\n\n", frameLabel[iFrame], totEv);
  }
}

//===================================
void ReadAccHistos(Char_t *fileNameMC, Bool_t normalise){

  printf("\n\n loading acceptance histos\n\n\n");
  TFile *fInMC = new TFile(fileNameMC);

  Char_t name[200];
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
      sprintf(name, "hAcc_Onia_cosTh_%s_pT%d", frameLabel[iFrame], iPTBin);
      hAcc_pol_pT[iFrame][iPTBin][cosThPol] = (TH1D *) gDirectory->Get(name);
      // printf("%p\n", hAcc_pol_pT[iFrame][iPTBin][cosThPol]);
      sprintf(name, "hAcc_Onia_phi_%s_pT%d", frameLabel[iFrame], iPTBin);      
      hAcc_pol_pT[iFrame][iPTBin][phiPol] = (TH1D *) gDirectory->Get(name);
      // printf("%p\n", hAcc_pol_pT[iFrame][iPTBin][phiPol]);
      sprintf(name, "hAcc2D_Onia_%s_pT%d", frameLabel[iFrame], iPTBin);
      hAcc2D_pol_pT[iFrame][iPTBin] = (TH2D *) gDirectory->Get(name);
      // printf("%p\n", hAcc2D_pol_pT[iFrame][iPTBin]);
    }
    for(int iRapBin = 0; iRapBin < 2*kNbRapBins+1; iRapBin++){
      sprintf(name, "hAcc_Onia_cosTh_%s_rap%d", frameLabel[iFrame], iRapBin);
      hAcc_pol_rap[iFrame][iRapBin][cosThPol] = (TH1D *) gDirectory->Get(name);
      // printf("%p\n", hAcc_pol_rap[iFrame][iRapBin][cosThPol]);
      sprintf(name, "hAcc_Onia_phi_%s_rap%d", frameLabel[iFrame], iRapBin);
      hAcc_pol_rap[iFrame][iRapBin][phiPol] = (TH1D *) gDirectory->Get(name);
      // printf("%p\n", hAcc_pol_rap[iFrame][iRapBin][phiPol]);
      sprintf(name, "hAcc2D_Onia_%s_rap%d", frameLabel[iFrame], iRapBin);
      hAcc2D_pol_rap[iFrame][iRapBin] = (TH2D *) gDirectory->Get(name);
      // printf("%p\n", hAcc2D_pol_rap[iFrame][iRapBin]);
    }
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
	sprintf(name, "hAcc_Onia_cosTh_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	hAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol] = (TH1D *) gDirectory->Get(name);
	// printf("%p\n", hAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]);
	sprintf(name, "hAcc_Onia_phi_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	hAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol] = (TH1D *) gDirectory->Get(name);
	// printf("%p\n", hAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]);
	sprintf(name, "hAcc2D_Onia_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) gDirectory->Get(name);
	// printf("%p\n", hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]);
	if(normalise)
	  hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Scale(1./hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Integral());
      }
    }
  }
}

//============================================
void Write1DFitResults(Char_t *fOutName){

  TFile *fOut = new TFile(fOutName, "RECREATE");
  //graphs from 1D fits:
  for(int iRap = 1; iRap <= kNbRapForPTBins; iRap++){
    gLambdaTh1D_rap[CS][iRap]->Write();
    gLambdaTh1D_rap[HX][iRap]->Write();
    gLambdaPhi1D_rap[CS][iRap]->Write();
    gLambdaPhi1D_rap[HX][iRap]->Write();
    gLambdaThPhi1D_rap[CS][iRap]->Write();
    gLambdaThPhi1D_rap[HX][iRap]->Write();
    gFrameInv1D_rap[CS][iRap]->Write();
    gFrameInv1D_rap[HX][iRap]->Write();
  }
  fOut->Close();
}

//============================================
void SaveCorrData(Char_t *label){

  Char_t name[200];
  sprintf(name, "polarization_corrData_%s.root", label);
  TFile *fOut = new TFile(name, "RECREATE");
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    // for(int iPTBin = 0; iPTBin < kNbPTBins+1; iPTBin++){
    //   hData_pol_pT[iFrame][iPTBin][cosThPol]->Write();
    //   hData_pol_pT[iFrame][iPTBin][phiPol]->Write();
    //   hData2D_pol_pT[iFrame][iPTBin]->Write();
    // }
    // for(int iRapBin = 0; iRapBin < 2*kNbRapBins+1; iRapBin++){
    //   hData_pol_rap[iFrame][iRapBin][cosThPol]->Write();
    //   hData_pol_rap[iFrame][iRapBin][phiPol]->Write();
    //   hData2D_pol_rap[iFrame][iRapBin]->Write();
    // }
    for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
      for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->Write();
	// hData_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->Write();
	hData2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	hData1D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }
}
