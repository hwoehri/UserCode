#include "/Users/hwoehri/Work/rootIncludes.inc"
#include "DataTnPTriggerEff.C"

Double_t binWidth = 0.01; //10 MeV
TF1 *fitFunc[kNbRapForPTBins+1][kNbPTMaxBins+1];
TF1 *fitFuncPassed[kNbRapForPTBins+1][kNbPTMaxBins+1];
//Int_t colour[kNbRapForPTBins+1] = {1, 1, 2};
Int_t colour[kNbRapForPTBins+1] = {1, 1, 2, 3, 4};
Double_t nTot[kNbRapForPTBins+1][kNbPTMaxBins+1];
Double_t errNTot[kNbRapForPTBins+1][kNbPTMaxBins+1];
Double_t nPassed[kNbRapForPTBins+1][kNbPTMaxBins+1];
Double_t errNPassed[kNbRapForPTBins+1][kNbPTMaxBins+1];
Double_t eps[kNbRapForPTBins+1][kNbPTMaxBins+1];
Double_t errEps[kNbRapForPTBins+1][kNbPTMaxBins+1];
TGraphAsymmErrors *gEff[kNbRapForPTBins+1];

void LoadHistos(Char_t *fileNameIn);
void FitHistos(Int_t iRapBin, Int_t iPTBin);
Double_t fitGaussExp(Double_t *x, Double_t *par);
Double_t fitCBExp(Double_t *x, Double_t *par);
void MakeEffGraph(Int_t iRapBin);
void PlotEffGraph(Char_t *hltName);
//===================================================
void calcEfficiency(Char_t *fileNameIn = "DataTnPTriggerEff_HLTDoubleMu0_RunB_15May2011.root",
		    //Char_t *fileNameIn = "DataTnPTriggerEff_HLTMu0TkMu0OSTJpsi_14May2011.root",){
		    Char_t *hltName = "HLT_DoubleMu0_RunB"){

  LoadHistos(fileNameIn);

  for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < kNbPTBins[iRapBin]+1; iPTBin++){
      binWidth = reco_M[iPTBin][iRapBin]->GetBinWidth(1);
      printf("binWidth = %f\n", binWidth);
      FitHistos(iRapBin, iPTBin);
    }
    MakeEffGraph(iRapBin);
  }
  PlotEffGraph(hltName);
}
//===================================================
void PlotEffGraph(Char_t *hltName){

  Char_t name[100];
  TCanvas *c1 = new TCanvas("cEff", "pT differential efficiency");
  TH1F *hFrame1 = gPad->DrawFrame(0., 0., 60., 1.2);
  hFrame1->SetXTitle("p_{T} [GeV/c]");
  hFrame1->SetYTitle("#varepsilon^{trigger}(#mu_{pos}) #upoint #varepsilon^{trigger}(#mu_{neg})");

  for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
    gEff[iRapBin]->Draw("p same");
  }
  TLegend *leg = new TLegend(0.6695402,0.1377119,0.9698276,0.2923729);
  for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
    if(iRapBin == 1)
      sprintf(name, "|y| < %1.2f", rapForPTRange[iRapBin]);
    else
      sprintf(name, "%1.2f < |y| < %1.2f", rapForPTRange[iRapBin-1], rapForPTRange[iRapBin]);
    leg->AddEntry(gEff[iRapBin], name, "lp");
  }
  leg->SetFillColor(0); 
  leg->SetBorderSize(0); leg->Draw();

  TLine *line = new TLine(0., 1., 60., 1.);
  line->SetLineStyle(3); line->Draw();

  TLatex *tex = new TLatex(2., 1.05, hltName);
  tex->SetTextSize(0.04); tex->Draw();

  sprintf(name, "Figures/triggerEff_TnP_orthogonalPDs_%s.pdf", hltName);
  c1->Print(name);
}

//===================================================
void MakeEffGraph(Int_t iRapBin){

  Double_t pT[kNbPTMaxBins], pTR[kNbPTMaxBins], pTL[kNbPTMaxBins];
  Double_t errEpsNeg[kNbPTMaxBins], errEpsPos[kNbPTMaxBins];
  for(int iPT = 0; iPT < kNbPTBins[iRapBin]; iPT++){
    pT[iPT] = pTRange[iRapBin][iPT] + 0.5*(pTRange[iRapBin][iPT+1] - pTRange[iRapBin][iPT]);
    pTL[iPT] = pT[iPT] - pTRange[iRapBin][iPT];
    pTR[iPT] = pTRange[iRapBin][iPT+1] - pT[iPT];
    errEpsNeg[iPT] = errEps[iRapBin][iPT];
    errEpsPos[iPT] = errEps[iRapBin][iPT];
    //correct the error...
    if((eps[iRapBin][iPT] + errEps[iRapBin][iPT]) > 1.)
      errEpsPos[iPT] = 1. - eps[iRapBin][iPT];
    if((eps[iRapBin][iPT] - errEps[iRapBin][iPT]) < 0.)
      errEpsNeg[iPT] = eps[iRapBin][iPT];

    printf("pT %1.2f, eff %1.3f +- %1.3f\n", pT[iPT], eps[iRapBin][iPT],  errEps[iRapBin][iPT]);
  }
  pT[kNbPTBins[iRapBin]-1] = 45.; pTL[kNbPTBins[iRapBin]-1] = 15., pTR[kNbPTBins[iRapBin]-1] = 15.;

  gEff[iRapBin] = new TGraphAsymmErrors(kNbPTBins[iRapBin], pT, eps[iRapBin], pTL, pTR, errEpsNeg, errEpsPos);
  gEff[iRapBin]->SetMarkerStyle(20);
  gEff[iRapBin]->SetMarkerColor(colour[iRapBin]);
  gEff[iRapBin]->SetLineColor(colour[iRapBin]);

}
//===================================================
void FitHistos(Int_t iRapBin, Int_t iPTBin){

  if(reco_M[iPTBin][iRapBin]->GetEntries() < 200.){//too little statistics
    nPassed[iRapBin][iPTBin] = 0;
    nTot[iRapBin][iPTBin] = 0;
    eps[iRapBin][iPTBin-1] = 0.;
    errEps[iRapBin][iPTBin-1] = 0.;
    return;
  }

  Char_t name[100];
  // gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(kTRUE);
  gStyle->SetStatX(0.5); gStyle->SetStatY(0.94);
  sprintf(name, "cTot_rap%d_pT%d\n", iRapBin, iPTBin);
  TCanvas *cTot = new TCanvas(name, name);

  sprintf(name, "fitFunc_rap%d_pt%d", iRapBin, iPTBin);
  fitFunc[iPTBin][iRapBin] = new TF1(name, fitCBExp, 2.6, 3.4, 7);
  fitFunc[iPTBin][iRapBin]->SetParameters(5000., 3.097, 0.03, 2., 2., 1e4, 1.);
  fitFunc[iPTBin][iRapBin]->SetParNames("N_{sig}", "#mu", "#sigma_{M}", "#alpha", "n", "N_{exp}", "tau");
  fitFunc[iPTBin][iRapBin]->SetParLimits(1, 3.0, 3.1); //mean
  fitFunc[iPTBin][iRapBin]->SetParLimits(2, 0.02, 0.1); //sigma
  fitFunc[iPTBin][iRapBin]->SetParLimits(3, 1., 10.); //alpha (determines how long the tail continues)
  fitFunc[iPTBin][iRapBin]->SetParLimits(4, 1., 2.); //n (determines how fast the hight of the tail goes to zero)

  reco_M[iPTBin][iRapBin]->SetMinimum(0.);
  reco_M[iPTBin][iRapBin]->Fit(name, "Q");
  fitFunc[iPTBin][iRapBin] = reco_M[iPTBin][iRapBin]->GetFunction(name);
  nTot[iRapBin][iPTBin] = fitFunc[iPTBin][iRapBin]->GetParameter(0);
  errNTot[iRapBin][iPTBin] = fitFunc[iPTBin][iRapBin]->GetParError(0);
  Double_t sigma = fitFunc[iPTBin][iRapBin]->GetParameter(2);
  Double_t alpha= fitFunc[iPTBin][iRapBin]->GetParameter(3);
  Double_t n = fitFunc[iPTBin][iRapBin]->GetParameter(4);

  sprintf(name, "Figures/nTot_rap%d_pT%d.pdf", iRapBin, iPTBin);
  cTot->Print(name);
  //=====================================================
  sprintf(name, "cPassed_rap%d_pT%d\n", iRapBin, iPTBin);
  TCanvas *cPassed = new TCanvas(name, name);

  sprintf(name, "fitFuncPassed_rap%d_pt%d", iRapBin, iPTBin);
  fitFuncPassed[iPTBin][iRapBin] = new TF1(name, fitCBExp, 2.6, 3.4, 7);
  fitFuncPassed[iPTBin][iRapBin]->SetParNames("N_{sig}", "#mu", "#sigma_{M}", "#alpha", "n", "N_{exp}", "tau");

  fitFuncPassed[iPTBin][iRapBin]->SetParameter(0, 0.5*fitFunc[iPTBin][iRapBin]->GetParameter(0));
  //impose the fit parameters
  fitFuncPassed[iPTBin][iRapBin]->FixParameter(1, fitFunc[iPTBin][iRapBin]->GetParameter(1));
  fitFuncPassed[iPTBin][iRapBin]->FixParameter(2, sigma);
  fitFuncPassed[iPTBin][iRapBin]->FixParameter(3, alpha);
  fitFuncPassed[iPTBin][iRapBin]->FixParameter(4, n);
  fitFuncPassed[iPTBin][iRapBin]->SetParameter(5, 0.5*fitFunc[iPTBin][iRapBin]->GetParameter(5));
  fitFuncPassed[iPTBin][iRapBin]->SetParameter(6, fitFunc[iPTBin][iRapBin]->GetParameter(6));

  passedTrigger_M[iPTBin][iRapBin]->SetMinimum(0.);
  passedTrigger_M[iPTBin][iRapBin]->Fit(name, "Q");
  fitFuncPassed[iPTBin][iRapBin] = passedTrigger_M[iPTBin][iRapBin]->GetFunction(name);
  nPassed[iRapBin][iPTBin] = fitFuncPassed[iPTBin][iRapBin]->GetParameter(0);
  errNPassed[iRapBin][iPTBin] = fitFuncPassed[iPTBin][iRapBin]->GetParError(0);

  sprintf(name, "Figures/nPassed_rap%d_pT%d.pdf", iRapBin, iPTBin);
  cPassed->Print(name);

  if(nTot[iRapBin][iPTBin] > 0.){
    if(iPTBin != 0)//save in the first indices the pT differential values
      eps[iRapBin][iPTBin-1] = nPassed[iRapBin][iPTBin] / nTot[iRapBin][iPTBin];
    else// and in the last bin the pT integrated value
      eps[iRapBin][kNbPTBins[iRapBin]-1] = nPassed[iRapBin][iPTBin] / nTot[iRapBin][iPTBin];

    if(nPassed[iRapBin][iPTBin] > 0.){
      if(iPTBin != 0){
	errEps[iRapBin][iPTBin-1] = eps[iRapBin][iPTBin-1]*(1.-eps[iRapBin][iPTBin-1]) / nPassed[iRapBin][iPTBin];
	errEps[iRapBin][iPTBin-1] = sqrt(errEps[iRapBin][iPTBin-1]);
      }
      else{
	errEps[iRapBin][kNbPTBins[iRapBin]-1] = eps[iRapBin][kNbPTBins[iRapBin]-1]*(1.-eps[iRapBin][kNbPTBins[iRapBin]-1]) / nPassed[iRapBin][kNbPTBins[iRapBin]-1];
	errEps[iRapBin][kNbPTBins[iRapBin]-1] = sqrt(errEps[iRapBin][kNbPTBins[iRapBin]-1]);
      }
    }
  }
  else{
    if(iPTBin != 0){
      eps[iRapBin][iPTBin-1] = 0.;
      errEps[iRapBin][iPTBin-1] = 0.;
    }
    else{
      eps[iRapBin][kNbPTBins[iRapBin]-1] = 0.;
      errEps[iRapBin][kNbPTBins[iRapBin]-1] = 0.;
    }
  }
  if(iPTBin != 0){
    printf("Npassed = %1.3f, Ntot = %1.3f --> efficiency = %1.3f +- %1.3f\n", 
	   nPassed[iRapBin][iPTBin], nTot[iRapBin][iPTBin],
	   eps[iRapBin][iPTBin-1], errEps[iRapBin][iPTBin-1]);
  }
}
//===================================================
void LoadHistos(Char_t *fileNameIn){

  TFile *fIn = new TFile(fileNameIn);

  Char_t name[100];
  for(int iRapBin = 0; iRapBin < kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < kNbPTBins[iRapBin]+1; iPTBin++){
      sprintf(name, "reco_M_pT%d_rap%d", iPTBin, iRapBin);
      reco_M[iPTBin][iRapBin] = (TH1D *) gDirectory->Get(name);

      sprintf(name, "passedTrigger_M_pT%d_rap%d", iPTBin, iRapBin);
      passedTrigger_M[iPTBin][iRapBin] = (TH1D *) gDirectory->Get(name);
      //      printf("%p\n", passedTrigger_M[iPTBin][iRapBin]);
    }
  }
}
//=========================
Double_t fitGaussExp(Double_t *x, Double_t *par){

  Double_t myExp = par[3]*TMath::Exp(-par[4]*x[0]);
  Double_t mean = par[1];
  Double_t sigma = par[2];

  Double_t gauss = par[0] / (TMath::Sqrt(2.*TMath::Pi()) * sigma) * TMath::Exp(-(pow(x[0]-mean,2)/(2.*sigma*sigma)));
  Double_t result = myExp + gauss;
  result *= binWidth; //correct for the bin width
  return result;
}
//=========================
Double_t fitCBExp(Double_t *x, Double_t *par){

  Double_t norm = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];
  Double_t alpha = par[3];
  Double_t n = par[4];
  Double_t myExp = par[5]*TMath::Exp(-par[6]*x[0]); 

  // Double_t CB;
  // if(((x[0] - mean)/sigma) > -alpha)
  //   CB = TMath::Exp(-(pow(x[0]-mean,2)/(2.*sigma*sigma)));
  // else{
  //   Double_t A = pow(n / fabs(alpha),n) * TMath::Exp(-alpha*alpha/2.);
  //   Double_t B = n / fabs(alpha) - fabs(alpha);
  //   CB = A*pow(B - (x[0]-mean)/sigma, -n);
  // }

  Double_t CB;
  Double_t t = (x[0] - mean)/sigma;
  if(alpha < 0) t = -t;
  Double_t absAlpha = fabs(alpha);
  if(t >= -absAlpha)
    CB = TMath::Exp(-0.5*t*t);
  else{
    Double_t a = pow(n/absAlpha,n)*TMath::Exp(-0.5*absAlpha*absAlpha);
    Double_t b = n/absAlpha - absAlpha;

    CB = a/pow(b - t, n);
  }
  Double_t result = myExp + norm / (TMath::Sqrt(2.*TMath::Pi()) * sigma) * CB;
  result *= binWidth; //correct for the bin width
  return result;
}
