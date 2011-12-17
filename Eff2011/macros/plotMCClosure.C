#include "../interface/commonVar.h"
#include "../interface/rootIncludes.inc"

Int_t const kNbSteps = 6; //0... all generated; 1... generated+Acc; 2... GEN+ACC+RECO, 3... GEN+ACC+RECO+TRIG, 4... eff-corrected[GEN+ACC+RECO+TRIG]
TH2D *hPTRap[kNbSteps];
TH2D *hPTRap_Ratio[kNbSteps];
Char_t *stepName[kNbSteps] = {"GEN", "GEN_ACC", "GEN_ACC_RECO", "GEN_ACC_RECO_TRIG", "effCorr", "effCorr_RECO"};
enum {GEN, GEN_ACC, GEN_ACC_RECO, GEN_ACC_RECO_TRIG, EFFCORR, EFFCORR_RECO};
TH1D *hCosTheta[eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1][kNbSteps];
TH1D *hPhi[eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1][kNbSteps];
TH1D *hCosTheta_Ratio[eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1][kNbSteps];
TH1D *hPhi_Ratio[eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1][kNbSteps];

void LoadHistos(Char_t *fileNameIn);
void BuildRatio(Int_t iNumerator, Int_t iDenom);
void PlotHistos();
void PlotRatioHistos(Int_t iFrame, Int_t iRapBin, Int_t iPTBin, Int_t iNumerator, Int_t iDenom);
//===============================================
void plotMCClosure(Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_16Dec2011_ProdSingleMuEff_woRhoFactor.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_16Dec2011_ProdSingleMuEff_withRhoFactor.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_15Dec2011_ProdSingleMuEff_withRhoFactor.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_14Dec2011_ProdSingleMuEff_noDimuVtxEffCorrection.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_14Dec2011_SingleMuEff_noDimuVtxEffCorrection.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_14Dec2011_ProdSingleMuEff_fitted.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_14Dec2011_SingleMuEff.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_14Dec2011_ProdSingleMuEff.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_14Dec2011_SingleMuEff_pTMin4c5GeV.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_14Dec2011_ProdSingleMuEff_pTMin4c5GeV.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_14Dec2011_ProdSingleMuEff.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_14Dec2011_ProdSingleMuEff.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_13Dec2011_ProdSingleMuEff_pTMin2c5GeV.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_13Dec2011_ProdSingleMuEff_pTMin4c5GeV.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_13Dec2011_SingleMuEff.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_13Dec2011_ProdSingleMuEff.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_13Dec2011_withTrackingEff.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_13Dec2011_noDimuVertexEff.root"){
  //Char_t *fileNameIn = "MCClosure_HLTDimuon10JpsiBarrel_13Dec2011.root"){

  LoadHistos(fileNameIn);
  PlotHistos();

  Int_t iFrame = 1;
  for(int iRapBin = 1; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
    for(int iPTBin = 6; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
      PlotRatioHistos(iFrame, iRapBin, iPTBin, EFFCORR, GEN_ACC);

}

//===============================================
void PlotRatioHistos(Int_t iFrame, Int_t iRapBin, Int_t iPTBin, Int_t iNumerator, Int_t iDenom){

  //CosTheta distribution
  Char_t name[100];
  sprintf(name, "hCosTheta_Ratio_%d_rap%d_pT%d_%s_%s", iFrame, iRapBin, iPTBin, stepName[iNumerator], stepName[iDenom]);
  hCosTheta_Ratio[iFrame][iRapBin][iPTBin][iNumerator] = (TH1D *) hCosTheta[iFrame][iRapBin][iPTBin][iNumerator]->Clone(name);
  hCosTheta_Ratio[iFrame][iRapBin][iPTBin][iNumerator]->Sumw2();
  hCosTheta_Ratio[iFrame][iRapBin][iPTBin][iNumerator]->Divide(hCosTheta[iFrame][iRapBin][iPTBin][iDenom]);

  sprintf(name, "c_CosTheta_Ratio_%d_rap%d_pT%d_%s_%s", iFrame, iRapBin, iPTBin, stepName[iNumerator], stepName[iDenom]);
  TCanvas *c1 = new TCanvas(name);
  TH1F *hFrame1 = gPad->DrawFrame(-1., 0.8, 1., 1.2);
  hFrame1->SetXTitle(hCosTheta[iFrame][iRapBin][iPTBin][iNumerator]->GetXaxis()->GetTitle());
  hFrame1->SetYTitle("1/#varepsilon #upoint (RECO + TRIG) / GEN");
  hCosTheta_Ratio[iFrame][iRapBin][iPTBin][iNumerator]->SetMarkerStyle(20);
  hCosTheta_Ratio[iFrame][iRapBin][iPTBin][iNumerator]->Draw("p same");

  sprintf(name, "%1.1f < |y| < %1.2f, %1.1f < p_{T} < %1.1f GeV/c", 
	  eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin],
	  eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);
  TLatex *tex1 = new TLatex(-0.9, 1.15, name);
  tex1->SetTextSize(0.04);
  tex1->Draw();

  TLine *line = new TLine(-1., 1., 1., 1.);
  line->SetLineStyle(3); line->Draw();

  if(iNumerator == EFFCORR)
    sprintf(name, "Figures/MCClosure_cosTheta_%s_rap%d_pT%d.pdf", eff::frameLabel[iFrame], iRapBin, iPTBin);
  else if(iNumerator == EFFCORR_RECO)
    sprintf(name, "Figures/MCClosure_smearing_cosTheta_%s_rap%d_pT%d.pdf", eff::frameLabel[iFrame], iRapBin, iPTBin);
  c1->Print(name);



  //phi distribution
  sprintf(name, "hPhi_Ratio_%d_rap%d_pT%d_%s_%s", iFrame, iRapBin, iPTBin, stepName[iNumerator], stepName[iDenom]);
  hPhi_Ratio[iFrame][iRapBin][iPTBin][iNumerator] = (TH1D *) hPhi[iFrame][iRapBin][iPTBin][iNumerator]->Clone(name);
  hPhi_Ratio[iFrame][iRapBin][iPTBin][iNumerator]->Sumw2();
  hPhi_Ratio[iFrame][iRapBin][iPTBin][iNumerator]->Divide(hPhi[iFrame][iRapBin][iPTBin][iDenom]);

  sprintf(name, "c_Phi_Ratio_%d_rap%d_pT%d_%s_%s", iFrame, iRapBin, iPTBin, stepName[iNumerator], stepName[iDenom]);
  TCanvas *c2 = new TCanvas(name);
  TH1F *hFrame2 = gPad->DrawFrame(-180., 0.8, 180., 1.2);
  hFrame2->SetXTitle(hPhi[iFrame][iRapBin][iPTBin][iNumerator]->GetXaxis()->GetTitle());
  hFrame2->SetYTitle("1/#varepsilon #upoint (RECO + TRIG) / GEN");
  hPhi_Ratio[iFrame][iRapBin][iPTBin][iNumerator]->SetMarkerStyle(20);
  hPhi_Ratio[iFrame][iRapBin][iPTBin][iNumerator]->Draw("p same");

  sprintf(name, "%1.1f < |y| < %1.2f, %1.1f < p_{T} < %1.1f GeV/c", 
	  eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin],
	  eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);
  TLatex *tex2 = new TLatex(-170., 1.15, name);
  tex2->SetTextSize(0.04);
  tex2->Draw();

  TLine *line2 = new TLine(-180., 1., 180., 1.);
  line2->SetLineStyle(3); line2->Draw();

  if(iNumerator == EFFCORR)
    sprintf(name, "Figures/MCClosure_phi_%s_rap%d_pT%d.pdf", eff::frameLabel[iFrame], iRapBin, iPTBin);
  else if(iNumerator == EFFCORR_RECO)
    sprintf(name, "Figures/MCClosure_smearing_phi_%s_rap%d_pT%d.pdf", eff::frameLabel[iFrame], iRapBin, iPTBin);
  c2->Print(name);
}

//===============================================
void PlotHistos(){

  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPaintTextFormat("5.2f");

  TCanvas *c1 = new TCanvas("c1", "", 700, 700);
  c1->Divide(2,3);
  for(int iStep = 0; iStep < kNbSteps; iStep++){
    c1->cd(iStep+1);
    hPTRap[iStep]->Draw("colz text");
  }

  BuildRatio(EFFCORR, GEN_ACC);
  TCanvas *c2 = new TCanvas("c2", "");
  gPad->SetLogy();  hPTRap_Ratio[EFFCORR]->GetYaxis()->SetMoreLogLabels();
  hPTRap_Ratio[EFFCORR]->GetYaxis()->SetRangeUser(9.999, 50);
  hPTRap_Ratio[EFFCORR]->SetMinimum(0.8);
  hPTRap_Ratio[EFFCORR]->SetMaximum(1.2);
  hPTRap_Ratio[EFFCORR]->SetMarkerSize(2.);
  hPTRap_Ratio[EFFCORR]->Draw("colz text");

  c2->Print("Figures/MCClosure.pdf");

  BuildRatio(EFFCORR_RECO, GEN_ACC);
  TCanvas *c3 = new TCanvas("c3", "");
  gPad->SetLogy();  hPTRap_Ratio[EFFCORR_RECO]->GetYaxis()->SetMoreLogLabels();
  hPTRap_Ratio[EFFCORR_RECO]->GetYaxis()->SetRangeUser(9.999, 50);
  hPTRap_Ratio[EFFCORR_RECO]->SetMinimum(0.8);
  hPTRap_Ratio[EFFCORR_RECO]->SetMaximum(1.2);
  hPTRap_Ratio[EFFCORR_RECO]->SetMarkerSize(2.);
  hPTRap_Ratio[EFFCORR_RECO]->Draw("colz text");

  c3->Print("Figures/MCClosure_smearing.pdf");

}

//===============================================
void BuildRatio(Int_t iNumerator, Int_t iDenom){

  Char_t name[100];
  sprintf(name, "hRatio_%s_%s", stepName[iNumerator], stepName[iDenom]);
  hPTRap_Ratio[iNumerator] = (TH2D *) hPTRap[iNumerator]->Clone(name);
  hPTRap_Ratio[iNumerator]->Sumw2();
  printf("will be dividing %d by %d\n", iNumerator, iDenom);
  hPTRap_Ratio[iNumerator]->Divide(hPTRap[iDenom]);
  
}

//===============================================
void LoadHistos(Char_t *fileNameIn){

  Char_t name[100];
  TFile *fIn = new TFile(fileNameIn);
  for(int iStep = 0; iStep < kNbSteps; iStep++){
    sprintf(name, "hPTRap_%s", stepName[iStep]);
    hPTRap[iStep] = (TH2D *) gDirectory->Get(name);
    printf("%s --> %p\n", name, hPTRap[iStep]);
  }

  for(int iStep = 0; iStep < kNbSteps; iStep++){
    for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
      for(int iRapBin = 0; iRapBin < eff::kNbRapBins+1; iRapBin++){
    	for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
   	  sprintf(name, "hCosTheta_%s_rap%d_pT%d_step%d", eff::frameLabel[iFrame], iRapBin, iPTBin, iStep);
	  hCosTheta[iFrame][iRapBin][iPTBin][iStep] = (TH1D *) gDirectory->Get(name);

	  sprintf(name, "hPhi_%s_rap%d_pT%d_step%d", eff::frameLabel[iFrame], iRapBin, iPTBin, iStep);
	  hPhi[iFrame][iRapBin][iPTBin][iStep] = (TH1D *) gDirectory->Get(name);
	}
      }
    }
  }
}


