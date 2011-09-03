#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"
#include "TMath.h"

Double_t binRange[100];
Double_t weightedCentre[100];
Int_t colour[6] = {1,1,2,4,3,6};
void Draw2DHist(TH2F *h2D_Pt_Rap, Char_t *hltTag);
void calcBinWidthPT(Int_t minNb, TH1F *hRecoPT);

//==========================================
void optimizeBinning(Int_t minNb = 6000, //min. nb. entries per pT bin
  		     Char_t *fileNameIn = "RootFiles/selEvents_data_Ups1Sonly_noRapCut_noCowboys_2Sep2011.root",
		     Char_t *hltTag = "8.4 < M < 11.6 GeV"){

  TFile *fIn = new TFile(fileNameIn);
  TH2F *h2D_Pt_Rap = (TH2F *) gDirectory->Get("Reco_Onia_rap_pt");
  Draw2DHist(h2D_Pt_Rap, hltTag);

  gStyle->SetPadRightMargin(0.02);

  TH1F *hReco[onia::kNbRapForPTBins+1];
  Char_t name[100];
  // for(int iRap = 0; iRap <= onia::kNbRapForPTBins; iRap++){
  for(int iRap = 0; iRap <= 2; iRap++){
    sprintf(name, "Reco_Onia_pt_rap%d", iRap);
    hReco[iRap] = (TH1F*) gDirectory->Get(name);
    hReco[iRap]->Rebin(5); //[100 MeV bins]
    if(iRap == 0) printf("rapidity: |y| < %1.1f\n", onia::rapForPTRange[onia::kNbRapForPTBins]);
    else if(iRap == 1) printf("rapidity: |y| < %1.1f\n", onia::rapForPTRange[iRap]);
    else printf("rapidity: %1.1f < |y| < %1.1f\n", onia::rapForPTRange[iRap-1], onia::rapForPTRange[iRap]);
    calcBinWidthPT(minNb, hReco[iRap]);
  }

  //==================================
  //all pT curves
  //==================================
  TCanvas *c1 = new TCanvas("c1", "pT distributions");
  gPad->SetLogy();
  TH1F *hFrame = gPad->DrawFrame(0., 0.5, 60., 2.*hReco[0]->GetMaximum());
  hFrame->SetXTitle("p_{T} [GeV/c]");
  // for(int iRap = 0; iRap <= onia::kNbRapForPTBins; iRap++){
  for(int iRap = 0; iRap <= 2; iRap++){
    hReco[iRap]->SetLineColor(colour[iRap]);
    hReco[iRap]->SetMarkerColor(colour[iRap]);
    hReco[iRap]->SetLineWidth(2);
    hReco[iRap]->Draw("same");
  }
  TLegend *leg = new TLegend(0.7,0.66,0.9755747,0.9449153, "M(#Upsilon(1S)) #pm 2 #upoint #sigma_{M}");
  leg->AddEntry(hReco[0], "all rapidities", "l");
  // for(int iRap = onia::kNbRapForPTBins; iRap > 1; iRap--){
  for(int iRap = 2; iRap > 1; iRap--){
    sprintf(name, "%1.1f < |y| < %1.1f\n", onia::rapForPTRange[iRap-1], onia::rapForPTRange[iRap]);
    leg->AddEntry(hReco[iRap], name, "l");
  }
  sprintf(name, "|y| < %1.1f\n", onia::rapForPTRange[1]);
  leg->AddEntry(hReco[1], name, "l");
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  TLatex *tex = new TLatex(7., 1., hltTag);
  tex->SetTextSize(0.04); tex->Draw();

  //  sprintf(name, "Figures/pTdist_%s.gif", hltTag);  c1->Print(name);
//   sprintf(name, "Figures/pTdist_%s.pdf", hltTag);  c1->Print(name);
  sprintf(name, "Figures/pTdist_Ups1S.gif");  c1->Print(name);
  sprintf(name, "Figures/pTdist_Ups1S.pdf");  c1->Print(name);

  //==================================
  //selected curves
  //==================================
  TCanvas *c2 = new TCanvas("c2", "pT distributions", 600, 600);
  gPad->SetLogy();
//   TH1F *hFrame2 = gPad->DrawFrame(0., 2., 30., 2.*hReco[0]->GetMaximum());
  TH1F *hFrame2 = gPad->DrawFrame(0., 2., 60., 2.*hReco[0]->GetMaximum());
  hFrame2->SetXTitle("p_{T} [GeV]");
  hFrame2->SetYTitle("events / 500 MeV");
  hFrame2->GetYaxis()->SetTitleOffset(1.6);

  // for(int iRap = 1; iRap < onia::kNbRapForPTBins; iRap++){
  for(int iRap = 1; iRap <= 2; iRap++){
    hReco[iRap]->Draw("same");
  }
  TLegend *leg2 = new TLegend(0.7,0.6,0.97,0.8, "M(#Upsilon(1S)) #pm 2 #upoint #sigma_{M}");
  // for(int iRap = onia::kNbRapForPTBins-1; iRap > 1; iRap--){
  for(int iRap = 2; iRap > 1; iRap--){
    sprintf(name, "%1.1f < |y| < %1.1f\n", onia::rapForPTRange[iRap-1], onia::rapForPTRange[iRap]);
    leg2->AddEntry(hReco[iRap], name, "l");
  }
  sprintf(name, "|y| < %1.1f\n", onia::rapForPTRange[1]);
  leg2->AddEntry(hReco[1], name, "l");
  leg2->SetTextSize(0.035);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->Draw();

  TLatex *tex2 = new TLatex(7., 3., hltTag);
  tex2->SetTextSize(0.04); tex2->Draw();
  // TLatex *tex2 = new TLatex(7.5, 1.2*hReco[0]->GetMaximum(), "CMS Preliminary,  #sqrt{s} = 7 TeV,  L = 39.6 pb^{-1}");
  // tex2->SetTextSize(0.04); tex2->Draw();

  // sprintf(name, "Figures/pTdist_%s_ANnote.gif", hltTag);  c2->Print(name);
  // sprintf(name, "Figures/pTdist_%s_ANnote.pdf", hltTag);  c2->Print(name);
  sprintf(name, "Figures/pTdist2_Ups1S.gif");  c2->Print(name);
  sprintf(name, "Figures/pTdist2_Ups1S.pdf");  c2->Print(name);
}

//========================================
void Draw2DHist(TH2F *h2D_Pt_Rap, Char_t *hltTag){

  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.12);
  Char_t name[100];
  TCanvas *c2 = new TCanvas("c2a", "2D distribution", 600, 600);
  gPad->SetLogz();
  //  TH1F *hFrame = gPad->DrawFrame(
//   h2D_Pt_Rap->Draw("box");
  h2D_Pt_Rap->GetYaxis()->SetTitleOffset(1.6);
  h2D_Pt_Rap->SetAxisRange(0,30.,"y");
  h2D_Pt_Rap->SetMaximum(4000.);
  h2D_Pt_Rap->Draw("colz");
  h2D_Pt_Rap->SetXTitle("y");
  h2D_Pt_Rap->SetYTitle("p_{T} [GeV]");

  TLine *lineV = new TLine();
  lineV->SetLineColor(0); lineV->SetLineWidth(2);
  for(int iRap = 0; iRap < 2*onia::kNbRapBins+1; iRap++)
    lineV->DrawLine(onia::rapRange[iRap], 0., onia::rapRange[iRap], 30.);

  TLine *lineH = new TLine();
  lineH->SetLineColor(0); lineH->SetLineWidth(2);
  for(int iRap = 1; iRap < onia::kNbRapForPTBins+1; iRap++){
    for(int iPT = 1; iPT < onia::kNbPTBins[iRap]+1; iPT++){
      lineH->DrawLine(onia::rapForPTRange[iRap-1], onia::pTRange[iRap][iPT], onia::rapForPTRange[iRap], onia::pTRange[iRap][iPT]);
      lineH->DrawLine(-onia::rapForPTRange[iRap-1], onia::pTRange[iRap][iPT], -onia::rapForPTRange[iRap], onia::pTRange[iRap][iPT]);
    }
  }

//   TLatex *tex = new TLatex(-2.3, 27., hltTag);
//   tex->SetTextSize(0.045);
//   tex->Draw();

  // sprintf(name, "Figures/ptVsRap_%s.gif", hltTag);  c2->Print(name);
  // sprintf(name, "Figures/ptVsRap_%s.pdf", hltTag);  c2->Print(name);
  sprintf(name, "Figures/ptVsRap_Upsilons.gif");  c2->Print(name);
  sprintf(name, "Figures/ptVsRap_Upsilons.pdf");  c2->Print(name);
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

  Double_t runSum[100];
  for(int iBin = 0; iBin < 100; iBin++)
    runSum[iBin] = 0.;
  
  Int_t count = 1;
  //1.) define the bin ranges:
  //  binRange[0] = hRecoPT->GetBinLowEdge(1);
  Int_t binMin = hRecoPT->GetXaxis()->FindBin(5.);
  binRange[0] = hRecoPT->GetBinLowEdge(binMin);
  for(int iBin = binMin; iBin <= hRecoPT->GetNbinsX(); iBin++){

    runSum[count-1] += hRecoPT->GetBinContent(iBin);
    // printf("histoBin %d, pTBin %d, runSum %1.3f\n", iBin, count-1, runSum[count-1]);
    if(runSum[count-1] > nbPerBin){

      binRange[count] = hRecoPT->GetBinLowEdge(iBin);
      count++;
    }
  }
  binRange[count] = hRecoPT->GetBinLowEdge(hRecoPT->GetNbinsX()+1);
  printf("for pT the optimal bins are\n"); 
  for(int iBin = 0; iBin <= count; iBin++){
    printf("%1.1f, ", binRange[iBin]);
  }
  printf("\n\n");
}
