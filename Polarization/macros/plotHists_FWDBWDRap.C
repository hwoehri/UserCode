#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

Int_t const kNbTrigger = 3;
//data:
//Char_t *trigLabel[kNbTrigger] = {"HLT_Mu0TkMu0Jpsi", "HLT_DoubleMu0"};
//MC:
Char_t *trigLabel[kNbTrigger] = {"RECO", "HLT_Mu0TkMu0_OST_Jpsi", "HLT_DoubleMu0"};
Double_t markerSize = 0.7;
Int_t markerStyle[kNbTrigger] = {24, 20, 21};

//reconstructed histos
TH1D *Reco_pol_pT_rap[kNbTrigger][jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][2*jpsi::kNbRapBins][jpsi::kNbPolVar];
TH2D *Reco2D_pol_pT_rap[kNbTrigger][jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][2*jpsi::kNbRapBins];

void ReadInRecHistos(Char_t *fileNameIn, Int_t iTrigger, Bool_t normalize, Int_t rebinCosTh, Int_t rebinPhi);
void PlotRecHistos(Int_t pT, Int_t rap, Int_t iFrame);
//===============================
void plotHists_FWDBWDRap(Bool_t normalize = kFALSE,
			 Int_t rebinCosTh = 1, 
			 Int_t rebinPhi = 1){

  Char_t name[100];
  Char_t fileNameIn[100];
  for(int iTr = 0; iTr < kNbTrigger; iTr++){
    sprintf(fileNameIn, "pol_MC_P_all_%s_rapDiff_2sigma_5Nov2010.root", trigLabel[iTr]);

    printf("reading in histos from file %s\n", fileNameIn);
    ReadInRecHistos(fileNameIn, iTr, normalize, rebinCosTh, rebinPhi);
  }
  //
  for(int iFrame = 0; iFrame < 2; iFrame++){
    for(int iRap = 0; iRap < 2*jpsi::kNbRapBins; iRap++){
      Int_t matchRapBin = fabs(jpsi::kNbRapForPTBins - iRap);
      if(iRap >= jpsi::kNbRapForPTBins) matchRapBin += 1;
      for(int iPT = 1; iPT < jpsi::kNbPTBins[matchRapBin]+1; iPT++){
	PlotRecHistos(iPT, iRap, iFrame);
	printf("will be plotting rap %d and pT bin %d\n", iRap, iPT);
      }
    }
  }
}

//===============================
void PlotRecHistos(Int_t iPTBin, Int_t iRapBin, Int_t iFrame){

  gStyle->SetOptStat(0);

  Int_t matchRapBin = fabs(jpsi::kNbRapForPTBins - iRapBin);
  if(iRapBin >= jpsi::kNbRapForPTBins) matchRapBin += 1;

  Char_t label[100], pTLabel[100], rapLabel[100];
  if(iPTBin == 0)
    sprintf(pTLabel, "all p_{T}");
  else if(iPTBin == 1)
    sprintf(pTLabel, "p_{T} < %1.1f", jpsi::pTRange[matchRapBin][iPTBin]);
  else
    sprintf(pTLabel, "%1.1f < p_{T} < %1.1f", jpsi::pTRange[matchRapBin][iPTBin-1], jpsi::pTRange[matchRapBin][iPTBin]);

  sprintf(rapLabel, "%1.1f < y < %1.1f", jpsi::rapRange[iRapBin], jpsi::rapRange[iRapBin+1]);

  sprintf(label, "%s, %s", rapLabel, pTLabel);

  Double_t maximum = Reco2D_pol_pT_rap[0][iFrame][iPTBin][iRapBin]->GetMaximum();

  Char_t name[100];
  sprintf(name, "c1_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
  TCanvas *c1 = new TCanvas(name, name, 1000, 700);
  c1->Divide(2,2);
  TH1F *hFrame1[4];

  c1->cd(1);
  gPad->SetRightMargin(0.1);
  Reco2D_pol_pT_rap[0][iFrame][iPTBin][iRapBin]->SetMaximum(maximum);
  Reco2D_pol_pT_rap[0][iFrame][iPTBin][iRapBin]->Draw("colz");

  TLatex *tex11 = new TLatex(-0.95, 320., trigLabel[0]);
  tex11->SetTextSize(0.05); tex11->Draw();

  //reconstructed cosT
  c1->cd(2);
//   gPad->SetLogy();
//   hFrame1[1] = gPad->DrawFrame(-1., 0., 1., 1.4*Reco_pol_pT_rap[0][iFrame][iPTBin][iRapBin][jpsi::cosThPol]->GetMaximum());
  hFrame1[1] = gPad->DrawFrame(-1., 1., 1., 1.4*Reco_pol_pT_rap[0][iFrame][iPTBin][iRapBin][jpsi::cosThPol]->GetMaximum());
  sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);
  hFrame1[1]->SetXTitle(name);

  for(int iTr = 0; iTr < kNbTrigger; iTr++){
    Reco_pol_pT_rap[iTr][iFrame][iPTBin][iRapBin][jpsi::cosThPol]->SetLineColor(iTr+1);
    Reco_pol_pT_rap[iTr][iFrame][iPTBin][iRapBin][jpsi::cosThPol]->SetMarkerColor(iTr+1);
    Reco_pol_pT_rap[iTr][iFrame][iPTBin][iRapBin][jpsi::cosThPol]->Draw("psame");
  }

  TLatex *tex1 = new TLatex(-0.9, 1.25*Reco_pol_pT_rap[0][iFrame][iPTBin][iRapBin][jpsi::cosThPol]->GetMaximum(), label);
  tex1->SetTextSize(0.05); tex1->Draw();

  TLegend *leg1 = new TLegend(0.55,0.7573165,0.9957329,0.9805308);
  for(int iTr = 0; iTr < kNbTrigger; iTr++){
sprintf(name, "%s: %1.0f", trigLabel[iTr], Reco2D_pol_pT_rap[iTr][iFrame][iPTBin][iRapBin]->Integral());
    leg1->AddEntry(Reco_pol_pT_rap[iTr][iFrame][iPTBin][iRapBin][jpsi::cosThPol], name, "p");
  }
  leg1->SetTextSize(0.05); leg1->SetFillColor(0); leg1->Draw();

  c1->cd(3);
  gPad->SetRightMargin(0.1);
  Reco2D_pol_pT_rap[1][iFrame][iPTBin][iRapBin]->Draw("colz");
  Reco2D_pol_pT_rap[1][iFrame][iPTBin][iRapBin]->SetMaximum(maximum);
  TLatex *tex13 = new TLatex(-0.95, 320., trigLabel[1]);
  tex13->SetTextSize(0.05); tex13->Draw();

  //reconstructed cosT
  c1->cd(4);
  hFrame1[3] = gPad->DrawFrame(0, 0., 360, 1.2*Reco_pol_pT_rap[0][iFrame][iPTBin][iRapBin][jpsi::phiPol]->GetMaximum());
  sprintf(name, "#phi_{%s}", jpsi::frameLabel[iFrame]);
  hFrame1[3]->SetXTitle(name);
  for(int iTr = 0; iTr < kNbTrigger; iTr++){
    Reco_pol_pT_rap[iTr][iFrame][iPTBin][iRapBin][jpsi::phiPol]->SetLineColor(iTr+1);
    Reco_pol_pT_rap[iTr][iFrame][iPTBin][iRapBin][jpsi::phiPol]->SetMarkerColor(iTr+1);
    Reco_pol_pT_rap[iTr][iFrame][iPTBin][iRapBin][jpsi::phiPol]->Draw("psame");
  }
  sprintf(name, "Figures/recoTrig_FWDBWD_pT%d_rap%d_%s.eps", iPTBin, iRapBin, jpsi::frameLabel[iFrame]);  c1->Print(name);
  sprintf(name, "Figures/recoTrig_FWDBWD_pT%d_rap%d_%s.gif", iPTBin, iRapBin, jpsi::frameLabel[iFrame]);  c1->Print(name);
  sprintf(name, "Figures/recoTrig_FWDBWD_pT%d_rap%d_%s.pdf", iPTBin, iRapBin, jpsi::frameLabel[iFrame]);  c1->Print(name);
}

//===============================
void ReadInRecHistos(Char_t *fileNameIn, Int_t iTrigger, Bool_t normalize, Int_t rebinCosTh, Int_t rebinPhi){

  TFile *fIn = new TFile(fileNameIn);
  Char_t name[100];

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins; iRapBin++){
      Int_t matchRapBin = fabs(jpsi::kNbRapForPTBins - iRapBin);
      if(iRapBin >= jpsi::kNbRapForPTBins) matchRapBin += 1;
      for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[matchRapBin]+1; iPTBin++){

	//generated
	sprintf(name, "Reco_Onia_cosTh_%s_pT%d_rap%d_FWDBWD", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	Reco_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin][jpsi::cosThPol] = (TH1D *) gDirectory->Get(name);
	sprintf(name, "%s_Trigger%d", Reco_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin][jpsi::cosThPol]->GetName(), iTrigger);
  	Reco_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin][jpsi::cosThPol]->SetMarkerStyle(markerStyle[iTrigger]);  
	Reco_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin][jpsi::cosThPol]->SetMarkerSize(markerSize);          
	//
	sprintf(name, "Reco_Onia_phi_%s_pT%d_rap%d_FWDBWD", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	Reco_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin][jpsi::phiPol] = (TH1D *) gDirectory->Get(name);
	sprintf(name, "%s_Trigger%d", Reco_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin][jpsi::phiPol]->GetName(), iTrigger);
	Reco_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin][jpsi::phiPol]->SetName(name);
 	Reco_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin][jpsi::phiPol]->SetMarkerStyle(markerStyle[iTrigger]); 
	Reco_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin][jpsi::phiPol]->SetMarkerSize(markerSize);         
	//
	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d_FWDBWD", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	Reco2D_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin] = (TH2D *) gDirectory->Get(name);
	sprintf(name, "%s_Trigger%d", Reco2D_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin]->GetName(), iTrigger);
	Reco2D_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin]->SetName(name);
	Reco2D_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin]->SetMarkerSize(markerSize);         
	if(rebinCosTh > 1 || rebinPhi > 1)
	  Reco2D_pol_pT_rap[iTrigger][iFrame][iPTBin][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
      }
    }
  }
}
