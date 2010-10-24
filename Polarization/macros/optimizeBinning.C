#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"
#include "TMath.h"

Double_t binRange[100];
Double_t weightedCentre[100];
void Draw2DHist(TH2F *h2D_Pt_Rap, Char_t *hltTag);
void calcBinWidthPT(Int_t minNb, TH1F *hRecoPT);
//==========================================
void optimizeBinning(Int_t minNb = 10000, //min. nb. entries per pT bin
		     Char_t *fileNameIn = "pol_data_all_HLT_Mu0TkMu0Jpsi_2sigma_21Oct2010.root",
		     Char_t *hltTag = "HLT_Mu0TkMu0_Jpsi"){

  TFile *fIn = new TFile(fileNameIn);
  TH2F *h2D_Pt_Rap = (TH2F *) gDirectory->Get("Reco_Onia_rap_pt");
  Draw2DHist(h2D_Pt_Rap, hltTag);

  TH1F *hReco[jpsi::kNbRapForPTBins+1];
  Char_t name[100];
  for(int iRap = 0; iRap <= jpsi::kNbRapForPTBins; iRap++){
    sprintf(name, "Reco_Onia_pt_rap%d", iRap);
    hReco[iRap] = (TH1F*) gDirectory->Get(name);
    if(iRap == 0) printf("rapidity: |y| < %1.1f\n", jpsi::rapForPTRange[jpsi::kNbRapForPTBins]);
    else if(iRap == 1) printf("rapidity: |y| < %1.1f\n", jpsi::rapForPTRange[iRap]);
    else printf("rapidity: %1.1f < |y| < %1.1f\n", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);
    calcBinWidthPT(minNb, hReco[iRap]);
  }

  // hReco[1] = (TH1F*) gDirectory->Get("Reco_Onia_rap_pT0");
  TCanvas *c1 = new TCanvas("c1", "pT distributions");
  gPad->SetLogy();
  TH1F *hFrame = gPad->DrawFrame(0., 0.5, 30., 2.*hReco[0]->GetMaximum());
  hFrame->SetXTitle("p_{T} [GeV/c]");
  for(int iRap = 0; iRap <= jpsi::kNbRapForPTBins; iRap++){
    hReco[iRap]->SetLineColor(jpsi::colour_rapForPTBins[iRap]);
    hReco[iRap]->Draw("same");
  }
  TLegend *leg = new TLegend(0.7,0.66,0.9755747,0.9449153, "M(J/#psi) #pm 2 #upoint #sigma_{M}");
  leg->AddEntry(hReco[0], "all rapidities", "l");
  for(int iRap = jpsi::kNbRapForPTBins; iRap > 1; iRap--){
    sprintf(name, "%1.1f < |y| < %1.1f\n", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);
    leg->AddEntry(hReco[iRap], name, "l");
  }
  sprintf(name, "|y| < %1.1f\n", jpsi::rapForPTRange[1]);
  leg->AddEntry(hReco[1], name, "l");
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  TLatex *tex = new TLatex(7., 1., hltTag);
  tex->SetTextSize(0.04); tex->Draw();

  sprintf(name, "Figures/pTdist_%s.gif", hltTag);
  c1->Print(name);
}

//========================================
void Draw2DHist(TH2F *h2D_Pt_Rap, Char_t *hltTag){

  Char_t name[100];
  TCanvas *c2 = new TCanvas("c2", "2D distribution");
  //  TH1F *hFrame = gPad->DrawFrame(
  h2D_Pt_Rap->Draw("box");
  TLine *lineV = new TLine();
  lineV->SetLineColor(4); lineV->SetLineWidth(2);
  for(int iRap = 0; iRap < 2*jpsi::kNbRapBins+1; iRap++)
    lineV->DrawLine(jpsi::rapRange[iRap], 0., jpsi::rapRange[iRap], 30.);

  TLine *lineH = new TLine();
  lineH->SetLineColor(4); lineH->SetLineWidth(2);
  for(int iRap = 1; iRap < jpsi::kNbRapForPTBins+1; iRap++){
    for(int iPT = 1; iPT < jpsi::kNbPTBins[iRap]+1; iPT++){
      lineH->DrawLine(jpsi::rapForPTRange[iRap-1], jpsi::pTRange[iRap][iPT], jpsi::rapForPTRange[iRap], jpsi::pTRange[iRap][iPT]);
      lineH->DrawLine(-jpsi::rapForPTRange[iRap-1], jpsi::pTRange[iRap][iPT], -jpsi::rapForPTRange[iRap], jpsi::pTRange[iRap][iPT]);
    }
  }

  TLatex *tex = new TLatex(-2.3, 27., hltTag);
  tex->SetTextSize(0.03);
  tex->Draw();
  sprintf(name, "Figures/ptVsRap_%s.gif", hltTag);  c2->Print(name);
  sprintf(name, "Figures/ptVsRap_%s.pdf", hltTag);  c2->Print(name);
}

//========================================
void calcBinWidthPT(Int_t minNb, TH1F *hRecoPT){

  Double_t tot, integral;
  Int_t nBins;
  Double_t nbPerBin;

  tot = hRecoPT->GetEntries();
  integral = hRecoPT->Integral();

  nBins = ((Int_t) tot) / minNb;

  nbPerBin = integral / (Double_t) nBins;

  printf("nb. entries %1.1f --> we can make %d bins with %1.1f entries\n", tot, nBins, nbPerBin);

  Double_t runSum = 0.;
  Int_t count = 1;
  //1.) define the bin ranges:
  binRange[0] = hRecoPT->GetBinLowEdge(1);
  for(int iBin = 1; iBin <= hRecoPT->GetNbinsX(); iBin++){

    runSum += hRecoPT->GetBinContent(iBin);
    if(runSum > count*nbPerBin){
      binRange[count] = hRecoPT->GetBinLowEdge(iBin);
      count++;
    }
  }
  binRange[count] = hRecoPT->GetBinLowEdge(hRecoPT->GetNbinsX()+1);
  printf("for pT the optimal bins are\n"); 
  for(int iBin = 0; iBin <= count; iBin++){
    printf("%f, ", binRange[iBin]);
  }
  printf("\n\n");
}
