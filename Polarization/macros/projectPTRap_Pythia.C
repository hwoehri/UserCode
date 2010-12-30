#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

Int_t colour[jpsi::kNbRapForPTBins+1] = {1,2,3,4,1,6};
//====================================================
void projectPTRap_Pythia(Bool_t plotNorm = kTRUE,
			 Char_t *fNameIn = "pythiaDistributions_JpsiWithFSR_7TeV.root"){


  TFile *fIn = new TFile(fNameIn);
  TH2D *hPtRap = (TH2D *) gDirectory->Get("hPtRap");
  
  TH1D *hPt[jpsi::kNbRapBins+1];
  TH1D *hRap[jpsi::kNbPTMaxBins+1];
  TH1D *hPtNorm[jpsi::kNbRapBins+1];
  TH1D *hRapNorm[jpsi::kNbPTMaxBins+1];


  Char_t name[100];
  for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++){
    sprintf(name, "hPt_rap%d", iRap);
    hPt[iRap] = (TH1D *) gDirectory->Get(name);
    hPt[iRap]->SetLineColor(colour[iRap]);
    sprintf(name, "%s_Clone", hPt[iRap]->GetName());
    hPtNorm[iRap] = (TH1D *) hPt[iRap]->Clone(name);
    hPtNorm[iRap]->Scale(4./(hPtNorm[iRap]->Integral() * hPtNorm[iRap]->GetBinWidth(1)));
  }
  for(int iPt = 0; iPt < jpsi::kNbPTMaxBins+1; iPt++){
    sprintf(name, "hRap_pT%d", iPt);
    hRap[iPt] = (TH1D *) gDirectory->Get(name);
    hRap[iPt]->SetLineColor(1+iPt);
    sprintf(name, "%s_Clone", hRap[iPt]->GetName());
    hRapNorm[iPt] = (TH1D *) hRap[iPt]->Clone(name);
    hRapNorm[iPt]->Scale(4./(hRapNorm[iPt]->Integral() * hRapNorm[iPt]->GetBinWidth(1)));
  }

  //================================================
  //pT distributions
  //================================================
  TCanvas *cPt = new TCanvas("cPt", "pT distributions");
  TH1F *hFramePt;
  if(plotNorm)
    hFramePt = gPad->DrawFrame(0., 1e-4, 30., 1.);
  else
    hFramePt = gPad->DrawFrame(0., 1., 30., 2.*hPt[0]->GetMaximum());
  hFramePt->SetXTitle("p_{T} [GeV/c]");
  gPad->SetLogy();
  //==============================================================
  //fitting the curves with: pT*[1 + (1/beta-2)*pT2/<pT2>]^(-beta)
  //fit parameters: beta and the av. pT^2, <pT2> 
  //==============================================================
  TF1  *fitPT[jpsi::kNbRapForPTBins+1];
  Double_t beta[jpsi::kNbRapForPTBins+1], pt2[jpsi::kNbRapForPTBins+1];
  for(int iRap  = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++){
    if(plotNorm)
      hPtNorm[iRap]->Draw("same");
    else
      hPt[iRap]->Draw("same");
    sprintf(name, "fitPT_rap%d", iRap);
    fitPT[iRap] = new TF1(name, "[0]*x*pow(1.+(1./([1]-2.))*x*x/[2],-[1])", 5., 30.);
    fitPT[iRap]->SetParName(0, "norm");
    fitPT[iRap]->SetParName(1, "beta");
    fitPT[iRap]->SetParName(2, "<pT2> [GeV2]");
    fitPT[iRap]->SetParameter(0, 1000.);
    fitPT[iRap]->SetParameter(1, 4.);
    fitPT[iRap]->SetParameter(2, 20.);

    if(plotNorm){
      hPtNorm[iRap]->Fit(name, "0+");
      fitPT[iRap] = hPtNorm[iRap]->GetFunction(name);
    }
    else{
      hPt[iRap]->Fit(name, "0+");
      fitPT[iRap] = hPt[iRap]->GetFunction(name);
    }
    fitPT[iRap]->SetLineWidth(1);
    fitPT[iRap]->SetLineColor(colour[iRap]);
    fitPT[iRap]->Draw("same");
    beta[iRap] = fitPT[iRap]->GetParameter(1);
    pt2[iRap] = fitPT[iRap]->GetParameter(2);
  }

  TLegend *legPT = new TLegend(0.5043103,0.6461864,0.9784483,0.9427966);
  for(int iRap  = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++){
    if(iRap == 0) sprintf(name, "all y, #beta = %1.2f, <p_{T}^{2}> = %1.1f GeV^{2}", beta[iRap], pt2[iRap]);
    else if(iRap == 1) sprintf(name, "|y| < %1.1f, #beta = %1.2f, <p_{T}^{2}> = %1.1f GeV^{2}", jpsi::rapForPTRange[iRap], beta[iRap], pt2[iRap]);
    else sprintf(name, "%1.1f < |y| < %1.1f, #beta = %1.2f, <p_{T}^{2}> = %1.1f GeV^{2}", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], beta[iRap], pt2[iRap]);
    legPT->AddEntry(hPt[iRap], name, "l");
  }
  legPT->SetTextSize(0.03); legPT->SetBorderSize(0);
  legPT->SetFillColor(0); legPT->Draw();

  TLatex *tex = new TLatex(2., 5., "PYTHIA6 with FSR");
  tex->SetTextSize(0.04); tex->Draw();
  tex->DrawLatex(2., 2., "prompt J/#psi @ 7 TeV");

  cPt->Print("Figures/pythia6_JpsiWithFSR_pT.gif");

  //=================================================
  //rapidity distributions
  //=================================================
  TCanvas *cRap = new TCanvas("cRap", "rapidity distributions");
  TH1F *hFrameRap;
  if(plotNorm)
    hFrameRap = gPad->DrawFrame(-2.5, 1e-1, 2.5, 2.);
  else{
    hFrameRap = gPad->DrawFrame(-2.5, 1., 2.5, 5.*hRap[0]->GetMaximum());
    gPad->SetLogy();
  }
  hFrameRap->SetXTitle("y");


  Double_t par1[jpsi::kNbPTMaxBins+1], par2[jpsi::kNbPTMaxBins+1];
  TF1 *fRap[jpsi::kNbPTMaxBins+1];
  for(int iPT = 0; iPT < jpsi::kNbPTMaxBins+1; iPT++){
    if(plotNorm)
      hRapNorm[iPT]->Draw("same");
    else
      hRap[iPT]->Draw("same");

    sprintf(name, "fRap_pT%d", iPT);
    fRap[iPT] = new TF1(name, "pol2", -2.5, 2.5);
    if(plotNorm){
      hRapNorm[iPT]->Fit(name, "0+");
      fRap[iPT] = hRapNorm[iPT]->GetFunction(name);
    }
    else{
      hRap[iPT]->Fit(name, "0+");
      fRap[iPT] = hRap[iPT]->GetFunction(name);
    }
    fRap[iPT]->SetLineColor(iPT+1);
    if(iPT > 0)
      fRap[iPT]->SetLineWidth(1);
    else
      fRap[iPT]->SetLineWidth(3);
    fRap[iPT]->Draw("same");
    par1[iPT] = fRap[iPT]->GetParameter(1);
    par2[iPT] = fRap[iPT]->GetParameter(2);
    printf("pT %d, par0 %1.3e, par1 %1.3e, par2 %1.3e\n", iPT, 
	   fRap[iPT]->GetParameter(0), par1[iPT], par2[iPT]);
  }

//   TLegend *legRap = new TLegend(0.5043103,0.6461864,0.9784483,0.9427966);
//   for(int iPt  = 0; iPt < jpsi::kNbPTMaxBins+1; iPt++){
//     if(iPt == 0) sprintf(name, "all p_{T}, p_{1} = %1.2f, p_{2} = %1.1f", par1[iPt], par2[iPt]);
//     else if(iPt == 1) sprintf(name, "p_{T} < %1.1f GeV/c, p_{1} = %1.2f, p_{2} = %1.1f", jpsi::pTRange[0][iPt], par1[iPt], par2[iPt]);
//     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c, p_{1} = %1.2f, p_{2} = %1.1f", jpsi::pTRange[0][iPt-1], jpsi::pTRange[0][iPt], par1[iPt], par2[iPt]);
//     legRap->AddEntry(hRap[iPt], name, "l");
//   }
//   legRap->SetTextSize(0.03); legRap->SetBorderSize(0);
//   legRap->SetFillColor(0);   legRap->Draw();

  TLatex *texRap = new TLatex(-2., 2.5, "PYTHIA6 with FSR");
  texRap->SetTextSize(0.04); texRap->Draw();
  texRap->DrawLatex(-2., 1.2, "prompt J/#psi @ 7 TeV");

  cRap->Print("Figures/pythia6_JpsiWithFSR_rap.gif");
}
