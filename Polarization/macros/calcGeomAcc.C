#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

#include "TStyle.h"

Double_t markerSize = 0.7;

//generated histos
TH2D *hGen_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
//generated histos with kinematical cuts (and smeared)
TH2D *hGenCut_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];

//acceptance histos
TH2D *hAcc2D_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];
TGraphAsymmErrors *gAcc2D_pT_rap[jpsi::kNbFrames][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1];

void ReadInHistos(Char_t *fileNameIn, Char_t *histoName, Int_t rebinCosTh, Int_t rebinPhi);
void CalcAcceptance();
void Plot2DAcceptance(Int_t iFrame);
void Plot2DAccOneByOne(Int_t iFrame, Int_t iRap, Int_t iPT);
void WriteAccHistos(Char_t *fileNameOut);
void CopyHistGraph(TH1D *hist, TGraphAsymmErrors*);
void CalcAcc(TH2D *hhGenCut, TH2D *hGen, TGraphAsymmErrors *graph);
//================================================================
void calcGeomAcc(Char_t *tag = "WithFSR_6Dec2010",
		 Char_t *histoName = "phiFolded_",
		 //Char_t *tag = "29Nov2010",
		 Int_t rebinCosTh = 1,
		 Int_t rebinPhi = 1){

  Char_t name[200];
  Char_t fileNameIn[200], fileNameOut[200];
  Char_t *processTag = "P_all";
  sprintf(fileNameIn,  "geomAcc_%s.root", tag);
  sprintf(fileNameOut, "geomAccHistos_%s.root", tag);

  ReadInHistos(fileNameIn, histoName, rebinCosTh, rebinPhi);

  CalcAcceptance();

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++)
    Plot2DAcceptance(iFrame);

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++)
    for(int iRap = 1; iRap <= jpsi::kNbRapForPTBins; iRap++)
      for(int iPT = 1; iPT <= jpsi::kNbPTBins[iRap]; iPT++)
	Plot2DAccOneByOne(iFrame, iRap, iPT);

  WriteAccHistos(fileNameOut);
}


//===============================
void Plot2DAcceptance(Int_t iFrame){

  gStyle->SetTitleW(1);

  Char_t name[100];
  Char_t title[100];
//   //=============================================
//   //1.) 2D acceptance for different pT bins
//   //=============================================
//   sprintf(name, "c1_2D_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c1_2D = new TCanvas(name, "acc vs cosTheta and phi for PT Bins", 1000, 700);
//   c1_2D->Divide(3,2);
//   TH2F *hFrame1a[jpsi::kNbPTBins+1];
//   TLatex *tex1[jpsi::kNbPTBins+1];
//   for(int iPT = 1; iPT <= jpsi::kNbPTBins; iPT++){

//     c1_2D->cd(iPT);

//     if(iPT == 0) sprintf(name, "J/#psi: all p_{T}");
//     else if(iPT == jpsi::kNbPTBins) sprintf(name, "J/#psi: p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPT-1]);
//     else sprintf(name, "J/#psi: %1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPT-1], jpsi::pTRange[iPT]);
//     hAcc2D_pT[iFrame][iPT]->SetTitle(name);
//     hAcc2D_pT[iFrame][iPT]->Draw("colz");//"colz" or "cont"
// //     hAcc2D_pT[iFrame][iPT]->Draw("lego2");//"colz" or "cont"
//   }
  
//   sprintf(name, "Figures/acceptance2D_%s_%s_pTBins.eps", jpsi::frameLabel[iFrame], hltTag);  c1_2D->Print(name);
//   sprintf(name, "Figures/acceptance2D_%s_%s_pTBins.pdf", jpsi::frameLabel[iFrame], hltTag);  c1_2D->Print(name);
//   sprintf(name, "Figures/acceptance2D_%s_%s_pTBins.gif", jpsi::frameLabel[iFrame], hltTag);  c1_2D->Print(name);

//   //========================================
//   //2.) 2D acceptance for different rap bins
//   //========================================
//   sprintf(name, "c2_2D_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c2_2D = new TCanvas(name, "acc vs cosTheta and phi for rap Bins", 1000, 700);
//   c2_2D->Divide(3,2);
//   TH2F *hFrame2a[jpsi::kNbRapBins+1]; 
//   for(int iRapBin = 1; iRapBin <= jpsi::kNbRapBins; iRapBin++){

//     c2_2D->cd(iRapBin);
//     if(iRapBin == 0) sprintf(name, "all y");
//     else sprintf(name, "J/#psi: %1.1f < y < %1.1f GeV/c", jpsi::rapRange[iRapBin-1], jpsi::rapRange[iRapBin]);
//     hAcc2D_rap[iFrame][iRapBin]->SetTitle(name);
//     hAcc2D_rap[iFrame][iRapBin]->Draw("colz"); //"colz" or "cont"
// //     hAcc2D_rap[iFrame][iRapBin]->Draw("lego2"); //"colz" or "cont"
//   }

//   sprintf(name, "Figures/acceptance2D_%s_%s_rapBins.eps", jpsi::frameLabel[iFrame], hltTag);  c2_2D->Print(name);
//   sprintf(name, "Figures/acceptance2D_%s_%s_rapBins.pdf", jpsi::frameLabel[iFrame], hltTag);  c2_2D->Print(name);
//   sprintf(name, "Figures/acceptance2D_%s_%s_rapBins.gif", jpsi::frameLabel[iFrame], hltTag);  c2_2D->Print(name);

  //=====================================================
  //3.) cosTheta acceptance for different rap and pT bins
  //=====================================================
  TCanvas *c3_2D[jpsi::kNbRapForPTBins];
  for(int iRap = 1; iRap <= jpsi::kNbRapForPTBins; iRap++){

    sprintf(name, "c3_2D_%s_rap%d", jpsi::frameLabel[iFrame], iRap);
    sprintf(title, "acc vs cosTheta and phi for pT; rap Bin %d", iRap);
    c3_2D[iRap] = new TCanvas(name, title, 1000, 700);
    c3_2D[iRap]->Divide(3,2);

    for(int iPT = 1; iPT <= jpsi::kNbPTBins[iRap]; iPT++){

      c3_2D[iRap]->cd(iPT);
      if(iPT == 0) 
	sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);
      else if(iPT == jpsi::kNbPTBins[iRap]) 
	sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iRap][iPT-1]);
      else 
	sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iRap][iPT-1], jpsi::pTRange[iRap][iPT]);
      hAcc2D_pT_rap[iFrame][iPT][iRap]->SetTitle(name);
      hAcc2D_pT_rap[iFrame][iPT][iRap]->Draw("colz");//"colz" or "cont"
//       hAcc2D_pT_rap[iFrame][iPT][iRap]->Draw("lego2");//"colz" or "cont"
    }

    sprintf(name, "Figures/geomAcc2D_%s_TBins_rap%d.eps", jpsi::frameLabel[iFrame], iRap); c3_2D[iRap]->Print(name);
    sprintf(name, "Figures/geomAcc2D_%s_TBins_rap%d.pdf", jpsi::frameLabel[iFrame], iRap); c3_2D[iRap]->Print(name);
    sprintf(name, "Figures/geomAcc2D_%s_TBins_rap%d.gif", jpsi::frameLabel[iFrame], iRap); c3_2D[iRap]->Print(name);
  }
}

//===============================
void Plot2DAccOneByOne(Int_t iFrame, Int_t iRap, Int_t iPT){

  gStyle->SetTitleW(1);

  // hAcc2D_pT_rap[iFrame][iPT][iRap]->SetMaximum(1.); //set the scale to 100 %

  Char_t name[100];
  Char_t title[100];

  //=========================================================
  //phi vs. cosTheta acceptance for different rap and pT bins
  //=========================================================
  TCanvas *c30_2D;
  sprintf(name, "c30_2D_%s_rap%d_pT%d", jpsi::frameLabel[iFrame], iRap, iPT);
  sprintf(title, "acc vs cosTheta and phi for pT %d rap Bin %d", iPT, iRap);
  c30_2D = new TCanvas(name, title, 500, 500);
  // TH1F *hFrame30 = gPad->DrawFrame(-1., 0., 1., 360.);
  // hFrame30->SetXTitle("hAcc2D_pT_rap[iFrame][iPT][iRap]->GetXTitle());
  // hFrame30->SetYTitle("hAcc2D_pT_rap[iFrame][iPT][iRap]->GetYTitle());

  if(iPT == 0) 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);
  else if(iPT == jpsi::kNbPTBins[iRap]) 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iRap][iPT-1]);
  else 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iRap][iPT-1], jpsi::pTRange[iRap][iPT]);
  hAcc2D_pT_rap[iFrame][iPT][iRap]->SetTitle(name);
  hAcc2D_pT_rap[iFrame][iPT][iRap]->Draw("colz");//"colz" or "cont"
  //       hAcc2D_pT_rap[iFrame][iPT][iRap]->Draw("lego2");//"colz" or "cont"

  sprintf(name, "Figures/geomAcc2D_%s_rap%d_pT%d.eps", jpsi::frameLabel[iFrame], iRap, iPT); 
  c30_2D->Print(name);
  sprintf(name, "Figures/geomAcc2D_%s_rap%d_pT%d.pdf", jpsi::frameLabel[iFrame], iRap, iPT); 
  c30_2D->Print(name);
  sprintf(name, "Figures/geomAcc2D_%s_rap%d_pT%d.gif", jpsi::frameLabel[iFrame], iRap, iPT); 
  c30_2D->Print(name);
}

//===============================
void CalcAcceptance(){

  Char_t name[100];

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){

// 	printf("preparing 2D differential acceptance (pT bin %d, rap bin %d)\n", iPTBin, iRapBin);
// 	gAcc2D_pT_rap[iFrame][iPTBin][iRapBin] = new TGraphAsymmErrors();
// 	sprintf(name, "gAcc2D_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
// 	gAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->SetName(name);
// 	Calc2DAcc(Reco2D_pT_rap[iFrame][iPTBin][iRapBin], 
// 		  hGen2D_pT_rap[iFrame][iPTBin][iRapBin],
// 		  gAcc2D_pT_rap[iFrame][iPTBin][iRapBin]);
      }
    }
    for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){

	sprintf(name, "hAcc2D_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	hAcc2D_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) hGenCut_pT_rap[iFrame][iPTBin][iRapBin]->Clone(name);
	Double_t nEntries;
// 	for(int iX = 1; iX <= hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(); iX++){
// 	  for(int iY = 1; iY <= hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY(); iY++){
// 	    nEntries = hGenCut_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iX, iY);
// 	    in case a bin has less than N events, set the acceptance to 0
// 	    and increase the bin error to something very large
// 	    if(nEntries < minEntriesPerBin){
// 	      hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(iX, iY, 0.);
// 	      hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->SetBinError(iX, iY, 100.);
// 	    }
// 	  }
// 	}

	hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->Divide(hGen_pT_rap[iFrame][iPTBin][iRapBin]);
	// hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->SetLineColor(jpsi::colour_pT[iPTBin]);
	// hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
	// hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
	hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->SetLineColor(jpsi::colour_rapForPTBins[iPTBin]);
	hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerColor(jpsi::colour_rapForPTBins[iPTBin]);
	hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerStyle(jpsi::marker_rapForPTBins[iPTBin]);
	hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerSize(markerSize);
	printf("%s: pT %d, rap %d: max. acceptance: %1.3e\n", jpsi::frameLabel[iFrame], iPTBin, iRapBin,
	       hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->GetMaximum());

	// //in case a bin content is higher than 1: set it back to 1
	// Double_t content;
	// for(int iX = 1; iX <= hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(); iX++){
	//   for(int iY = 1; iY <= hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY(); iY++){
	//     content = hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iX, iY);
	//     if(content > 1.)
	//       hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(iX, iY, 1.);
	//   }
	// }
      }
    }
  }
}

//===============================
void ReadInHistos(Char_t *fileNameIn, Char_t *histoName, Int_t rebinCosTh, Int_t rebinPhi){

  TFile *fIn = new TFile(fileNameIn);
  Char_t name[100];

  Int_t totEntries[jpsi::kNbFrames];
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    totEntries[iFrame] = 0;
    for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	sprintf(name, "hGen_%s%s_pT%d_rap%d", histoName, jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	hGen_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) gDirectory->Get(name);
	// hGen_pT_rap[iFrame][iPTBin][iRapBin]->SetLineColor(jpsi::colour_pT[iPTBin]);   
	// hGen_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]); 
	// hGen_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]); 
	hGen_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerSize(markerSize);         
	if(rebinCosTh > 1 || rebinPhi > 1)
	  hGen_pT_rap[iFrame][iPTBin][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
	printf("%s: pTbin %d, rapBin %d --> statistics in  GEN: %1.1f\n",
	       jpsi::frameLabel[iFrame], iPTBin, iRapBin, 
	       hGen_pT_rap[iFrame][iPTBin][iRapBin]->GetEntries());

	//reconstructed
	sprintf(name, "hGenCut_%s%s_pT%d_rap%d", histoName, jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	hGenCut_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) gDirectory->Get(name);
	// hGenCut_pT_rap[iFrame][iPTBin][iRapBin]->SetLineColor(jpsi::colour_pT[iPTBin]);   
	// hGenCut_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]); 
	// hGenCut_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]); 
	hGenCut_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerSize(markerSize);         
	printf("%s: pTbin %d, rapBin %d --> statistics in RECO: %1.1f\n",
	       jpsi::frameLabel[iFrame], iPTBin, iRapBin, 
	       hGenCut_pT_rap[iFrame][iPTBin][iRapBin]->GetEntries());

	totEntries[iFrame] += hGenCut_pT_rap[iFrame][iPTBin][iRapBin]->GetEntries();
	if(rebinCosTh > 1 || rebinPhi > 1)
	  hGenCut_pT_rap[iFrame][iPTBin][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
      }
    }
    printf("%s: integrated statistics in RECO is %d\n", 
	   jpsi::frameLabel[iFrame], totEntries[iFrame]);
  }
}

//===============================
void WriteAccHistos(Char_t *fileNameOut){

  printf("<WriteAccHistos> writing out the histograms\n");

  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
//     for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
//       for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){

// 	gAcc_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->Write();
// 	gAcc_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->Write();
// 	gAcc_pT_rap[iFrame][iPTBin][iRapBin]->Write();
//       }
//     }
    for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	hAcc2D_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }
  fOut->Close();
}

//=======================================================
void CalcAcc(TH2D *hhGenCut, TH2D *hGen, TGraphAsymmErrors *graph){

  Int_t nBinsX1 = hhGenCut->GetNbinsX();
  Int_t nBinsX2 = hGen->GetNbinsX();
  Int_t nBinsY1 = hhGenCut->GetNbinsY();
  Int_t nBinsY2 = hGen->GetNbinsY();
  if(nBinsX1 != nBinsX2 && nBinsY1 != nBinsY2) return;

  //build 1D histos containing the slices of the  histo appended
  Char_t name[100];
  sprintf(name, "%s_1D", hGen->GetName());
  TH1D *hGen1D = new TH1D(name, "", nBinsX1*nBinsY1, 0., nBinsX1*nBinsY1);
  sprintf(name, "%s_1D", hhGenCut->GetName());
  TH1D *hhGenCut1D = new TH1D(name, "", nBinsX1*nBinsY1, 0., nBinsX1*nBinsY1);
  Double_t gen, rec;
  for(int iY = 1; iY <= nBinsY1; iY++){
    for(int iX = 1; iX <= nBinsX1; iX++){
      gen = hGen->GetBinContent(iX, iY);
      rec = hhGenCut->GetBinContent(iX, iY);

      hGen1D->SetBinContent((iY-1)*nBinsX1 + iX, gen);
      hhGenCut1D->SetBinContent((iY-1)*nBinsX1 + iX, rec);
    }
  }

  //use Bayes divide to get the errors
  graph->BayesDivide(hhGenCut1D, hGen1D);
  
  // //***use FeldmanCousinsBinomialInterval/ClopperPearsonBinomialInterval to get the errors
  // //BinomialInterval::FeldmanCousinsBinomialInterval fc;
  // FeldmanCousinsBinomialInterval fc;
  // //ClopperPearsonBinomialInterval cp;

  // //alpha = 1 - CL
  // const double alpha = (1-0.682);
  // fc.init(alpha);
  // //cp.init(alpha);

  // Int_t binX = hGen1D->GetNbinsX();

  // for(int iX = 1; iX <= binX; iX++){
  //     gen = hGen1D->GetBinContent(iX);
  //     rec = hhGenCut1D->GetBinContent(iX);

  //     double acc = rec / gen;
  //     fc.calculate(rec, gen);
  //     //cp.calculate(rec, gen);
  //     double errorLow = acc - fc.lower();
  //     double errorHigh = fc.upper() - acc;

  //     Double_t xVal =  hGen1D->GetBinCenter(iX);

  //     graph->SetPoint(iX,xVal,acc);
  //     graph->SetPointEYhigh(iX, errorHigh);
  //     graph->SetPointEYlow(iX, errorLow);
  // }

}

//=========================================
void CopyHistGraph(TH1D *hist, TGraphAsymmErrors* graph){

  graph->Set(hist->GetNbinsX());

  Double_t centre, content, error;
  Double_t integral = hist->GetEntries();
  for(int iBin = 1; iBin <= hist->GetNbinsX(); iBin++){
    centre = hist->GetBinCenter(iBin);
    content = hist->GetBinContent(iBin);
    error = hist->GetBinError(iBin);
    graph->SetPoint(iBin-1, centre, content/integral);
    graph->SetPointError(iBin-1, 
		  hist->GetBinCenter(iBin)-hist->GetBinLowEdge(iBin),
		  hist->GetBinLowEdge(iBin+1)-hist->GetBinCenter(iBin),
		  error/integral, error/integral);
  }
}
