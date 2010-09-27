#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

TH1F *hDeltaPhi[jpsi::kNbPTBins+1][2*jpsi::kNbRapBins+1];

void PlotHistos(Int_t iPTBin, Int_t iRapBin, Int_t oniaID);
void ReadHistos(Char_t *fileNameIn);
Char_t *oniaLabel[100] = {"J/#psi", "Ups(1S)"};
//==============================
void plotDeltaPhi(Char_t *fileNameIn = "pol_data_all_HLT_Mu0TkMu0Jpsi_27Sep2010.root", 
		  Int_t oniaID = 0){ //0... J/psi, 1...Ups(1S)

  ReadHistos(fileNameIn);

  for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++)
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++)  
      PlotHistos(iPTBin, iRapBin, oniaID);

}
//==============================
void PlotHistos(Int_t iPTBin, Int_t iRapBin, Int_t oniaID){

  Char_t name[100];
  sprintf(name, "c1_pT%d_y%d", iPTBin, iRapBin);
  TCanvas *c1 = new TCanvas(name, name);
  TH1F *hFrame = gPad->DrawFrame(-TMath::Pi()/2., 0., TMath::Pi()/2., 1.4*hDeltaPhi[iPTBin][iRapBin]->GetMaximum());
  hFrame->SetXTitle(hDeltaPhi[iPTBin][iRapBin]->GetXaxis()->GetTitle());
  hDeltaPhi[iPTBin][iRapBin]->Draw("same");

  Int_t nBins = hDeltaPhi[iPTBin][iRapBin]->GetNbinsX();
  Double_t all = hDeltaPhi[iPTBin][iRapBin]->Integral();
  // printf("total % bins = %d, seagulls integrated from bin 1 to %d; coyboys from %d to %d\n",
  // 	 nBins, nBins/2, nBins/2+1, nBins);
  Double_t seagulls = hDeltaPhi[iPTBin][iRapBin]->Integral(1, nBins/2);
  Double_t cowboys = hDeltaPhi[iPTBin][iRapBin]->Integral(nBins/2+1, nBins);

  if(all > 0){
    cowboys /= all;
    seagulls /= all;
  }

  if(all > 10){
    if(iRapBin == 1){
      sprintf(name, "|y(%s)| < %1.1f, %1.1f < p_{T} < %1.1f GeV/c", 
	      oniaLabel[oniaID], jpsi::rapForPTRange[iRapBin], 
	      jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
    }
    else{
      sprintf(name, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T} < %1.1f GeV/c", 
	      jpsi::rapForPTRange[iRapBin-1], 
	      oniaLabel[oniaID], jpsi::rapForPTRange[iRapBin],
	      jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
    }
    TLatex *tex = new TLatex(-1.4, 1.3*hDeltaPhi[iPTBin][iRapBin]->GetMaximum(), name);
    tex->SetTextSize(0.04);  tex->Draw();
    sprintf(name, "fraction of cowboys, #phi(#mu^{+})-#phi(#mu^{-}) > 0: %1.1f %%\n", 100.*cowboys);
    tex->DrawLatex(-1.4, 1.2*hDeltaPhi[iPTBin][iRapBin]->GetMaximum(), name);
  }

  sprintf(name, "Figures/deltaPhiLab-pT%d-rap%d.gif", iPTBin, iRapBin);
  c1->Print(name);
  sprintf(name, "Figures/deltaPhiLab-pT%d-rap%d.pdf", iPTBin, iRapBin);
  c1->Print(name);
}

//==============================
void ReadHistos(Char_t *fileNameIn){

  TFile *fIn = new TFile(fileNameIn);
  Char_t name[100];

  for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    for(int iRapBin = 1; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
      sprintf(name, "Reco_hDeltaPhi_pT%d_rap%d", iPTBin, iRapBin);
      hDeltaPhi[iPTBin][iRapBin] = (TH1F *) gDirectory->Get(name);
    }
  }
}
