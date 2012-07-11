#include "../interface/rootIncludes.inc"
#include "TEfficiency.h"

//tracker80 bins:
// Int_t const kNbpTTNP = 52;
// Double_t pTBinsTNP[kNbpTTNP+1] = {2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.2, 4.4, 4.6, 4.8, 5., 5.2, 5.4, 5.6, 5.8, 6., 6.2, 6.4, 6.6, 6.8, 7., 7.2, 7.4, 7.6, 7.8, 8., 8.2, 8.4, 8.6, 8.8, 9, 9.2, 9.4, 9.6, 9.8, 10., 10.2, 10.4, 10.6, 10.8, 11, 11.2, 11.4, 11.6, 11.8, 12, 13., 14, 15};
Int_t const kNbpTTNP = 17;
Double_t pTBinsTNP[kNbpTTNP+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 10., 15., 20., 30., 50.};
Int_t const kNbpTTNPALLINONE = 14;
Double_t pTBinsTNPALLINONE[kNbpTTNPALLINONE+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10., 15., 20., 50.};
Int_t const kNbEta = 10;
Double_t etaBins[kNbEta+1] = {0, 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.1, 2.4};
// Int_t const kNbEta = 18;
// Double_t etaBins[kNbEta+1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2.1, 2.4};

Double_t minDeltaPT = 0.5; //set the minDeltaPT in which the T&P eff. were studied; is crucial for the extrapolation

Int_t const kNbType = 3;
enum{TNP, TNPALLINONE, MCTRUTH};
TGraphAsymmErrors *gEff[3][kNbEta];
TEfficiency *tEff[kNbEta]; //MC truth stored as TEfficiency

TH1F *hRatio[2][kNbEta];
enum{DATA,MC};
Char_t *effSampleName[2] = {"DATA", "MC"};

void LoadEfficiencies(Int_t effSample);
void DrawEff(Int_t iEta, Bool_t drawTNP, Bool_t drawAllInOne, Int_t effSample);
void MakeRatio(Int_t iEta, Int_t index, Bool_t useGraph_MCTRUTH, Bool_t useExtrapolation);
void DrawRatio(Int_t iEta, Int_t effSample);
//==============================
void compareEff(Int_t effSample = 1, //0...DATA, 1... MC
		Bool_t useExtrapolation = kTRUE
		){

  LoadEfficiencies(effSample);

  for(int iEta = 0; iEta < kNbEta-1; iEta++)
    DrawEff(iEta, kFALSE, kFALSE, effSample); //MCtruth alone

  for(int iEta = 0; iEta < kNbEta-1; iEta++)
    DrawEff(iEta, kTRUE, kTRUE, effSample);

  for(int iEta = 0; iEta < kNbEta-1; iEta++)
    DrawEff(iEta, kTRUE, kFALSE, effSample);

  for(int iEta = 0; iEta < kNbEta-1; iEta++)
    DrawEff(iEta, kFALSE, kTRUE, effSample);

  //Bool_t useGraph_MCTRUTH = kFALSE;
  Bool_t useGraph_MCTRUTH = kTRUE;
  for(int iEta = 0; iEta < kNbEta-1; iEta++){
  //for(int iEta = 0; iEta < 1; iEta++){
    MakeRatio(iEta, TNP, useGraph_MCTRUTH, useExtrapolation);
    MakeRatio(iEta, TNPALLINONE, useGraph_MCTRUTH, useExtrapolation);
    DrawRatio(iEta, effSample);
  }
}
//==============================
void DrawRatio(Int_t iEta, Int_t effSample){

  Char_t name[100];
  sprintf(name, "c2_%d", iEta);
  TCanvas *c2 = new TCanvas(name, name);
  TH1F *hFrame2 = gPad->DrawFrame(0., 0.5, 20.1, 1.5);
  hFrame2->SetXTitle("p_{T} [GeV/c]");
  hFrame2->SetYTitle("T&P / Truth");

  hRatio[TNP][iEta]->SetMarkerStyle(20);
  hRatio[TNP][iEta]->SetMarkerColor(1);
  hRatio[TNP][iEta]->Draw("p same");

  // hRatio[TNPALLINONE][iEta]->SetMarkerStyle(24);
  // hRatio[TNPALLINONE][iEta]->SetMarkerColor(4);
  // hRatio[TNPALLINONE][iEta]->Draw("p same");

  TLine *line2 = new TLine(0., 1., 20.1, 1.);
  line2->SetLineStyle(3); line2->Draw();

  sprintf(name, "%1.1f < |#eta| < %1.1f", etaBins[iEta], etaBins[iEta+1]);
  TLegend *leg2 = new TLegend(0.5991379,0.7436441,0.9022989,0.8813559, name);
  if(effSample == MC){
    // leg2->AddEntry(hRatio[TNPALLINONE][iEta], "all in one MC T&P", "pl");
    // leg2->AddEntry(hRatio[TNP][iEta], "product MC T&P", "pl");
    leg2->AddEntry(hRatio[TNPALLINONE][iEta], "inclusive MC T&P", "pl");
    leg2->AddEntry(hRatio[TNP][iEta], "factorized MC T&P", "pl");
  }
  else if(effSample == DATA){
    leg2->AddEntry(hRatio[TNPALLINONE][iEta], "all in one Data T&P", "pl");
    leg2->AddEntry(hRatio[TNP][iEta], "product Data T&P", "pl");
  }
  leg2->SetTextSize(0.04); leg2->SetBorderSize(0);
  leg2->SetFillColor(0); leg2->Draw();

  sprintf(name, "Figures/ratioEff_Eta%d.pdf", iEta);
  c2->Print(name);
}

//==============================
void MakeRatio(Int_t iEta, Int_t index, Bool_t useGraph_MCTRUTH, Bool_t useExtrapolation){

  printf("index %d\n", index);

  Char_t name[100];
  sprintf(name, "hRatio_%d_eta%d", index, iEta);
  // Int_t maxPT = kNbpTTNP;
  Int_t maxPT = gEff[index][iEta]->GetN();
  if(index == 0 || index == 1){
    hRatio[index][iEta] = new TH1F(name, ";p_{T} [GeV/c]", kNbpTTNP, pTBinsTNP);
     // maxPT = kNbpTTNP;
  }
  // else if(index == 1){
  //   hRatio[index][iEta] = new TH1F(name, ";p_{T} [GeV/c]", kNbpTTNPALLINONE, pTBinsTNPALLINONE);
  //   maxPT = kNbpTTNPALLINONE;
  // }
  Double_t pT=0, eff=0, effTruth=0;
  Double_t pT1=0, eff1=0, pTTruthCentre=0;
  Double_t pT2=0, eff2=0;
  Double_t effTruthMinusOne=0, effTruthPlusOne=0, pTMinusOne=0, pTPlusOne=0;
  Int_t thisBin = -1;
  Double_t k=0, d=0, effTruthExtrap=0;
  Int_t binMin = 0;
  Bool_t binMinAssigned = kFALSE;
  for(int iPT = 0; iPT < maxPT; iPT++){

    gEff[index][iEta]->GetPoint(iPT, pT, eff); //pT ... centre-of-gravity
    binMinAssigned = kFALSE;

    if(eff < 0.05 )
       continue;

    if(pT > 19.9) //H: temporary protection for missing data
      continue;

    // for(int iPTBin = 0; iPTBin < kNbpT; iPTBin++){
    //   if(pT > pTBins[iPTBin] && pT < pTBins[iPTBin+1]){
    // 	thisBin = iPTBin;
    // 	break;
    //   }
    // }

    //uncomment if no extrapolation is wanted...
    if(!useExtrapolation){
      thisBin = tEff[iEta]->FindFixBin(pT);
      effTruth = tEff[iEta]->GetEfficiency(thisBin);
      printf("pTBin %d, eff %1.3f, effTruth %1.3f\n", thisBin, eff, effTruth);
      hRatio[index][iEta]->Fill(pT, eff / effTruth);
    }
    else{
      if(!useGraph_MCTRUTH){

	thisBin = tEff[iEta]->FindFixBin(pT);
	effTruth = tEff[iEta]->GetEfficiency(thisBin);
	pTTruthCentre = tEff[iEta]->GetPassedHistogram()->GetBinCenter(thisBin);
	if(pT < pTTruthCentre && thisBin > 0){
	  effTruthMinusOne = tEff[iEta]->GetEfficiency(thisBin-1);
	  pTMinusOne = tEff[iEta]->GetPassedHistogram()->GetBinCenter(thisBin-1);
	}
	else{
	  effTruthPlusOne = tEff[iEta]->GetEfficiency(thisBin+1);
	  pTPlusOne = tEff[iEta]->GetPassedHistogram()->GetBinCenter(thisBin+1);
	}
      }
      else{
	
	Double_t pTOld = 1000.;
	for(int iPoint = 0; iPoint < gEff[MCTRUTH][iEta]->GetN(); iPoint++){
	  
	  gEff[MCTRUTH][iEta]->GetPoint(iPoint, pTTruthCentre, effTruth);
	  if(effTruth > 0.001 && binMinAssigned == kFALSE){
	    binMin = iPoint;
	    binMinAssigned = kTRUE;
	  }

	  if(fabs(pT - pTTruthCentre) < fabs(pTOld - pT)){
	    thisBin = iPoint;
	    printf("pTBin %d, pT %1.3f, pTOld %1.3f, pTTruthCentre %1.3f (point %d), thisBin %d\n",
		   iPT, pT, pTOld, pTTruthCentre, iPoint, thisBin);

	    // if(fabs(pT - pTTruthCentre) < minDeltaPT)
	    //   break;
	    pTOld = pTTruthCentre;
	  }
	  if(effTruth < 0.05)
	     continue;

	}
	//printf("thisBin %d\n", thisBin);
	// thisBin = tEff[iEta]->FindFixBin(pT);
	gEff[MCTRUTH][iEta]->GetPoint(thisBin, pTTruthCentre, effTruth);
	if(pT < pTTruthCentre && thisBin > binMin)
	  gEff[MCTRUTH][iEta]->GetPoint(thisBin-1, pTMinusOne, effTruthMinusOne);
	else
	  gEff[MCTRUTH][iEta]->GetPoint(thisBin+1, pTPlusOne, effTruthPlusOne);
      }
      //      if(pT < pTTruthCentre && thisBin > binMin){
	 
	 if(pT < pTTruthCentre){
	   eff1 = effTruthMinusOne;
	   pT1 = pTMinusOne;
	   pT2 = pTTruthCentre;
	   eff2 = effTruth;
	 }
	 else{
	   eff1 = effTruth;
	   pT1 = pTTruthCentre;
	   pT2 = pTPlusOne;
	   eff2 = effTruthPlusOne;
	 }


      k = (eff2 - eff1) / (pT2 - pT1);
      d = eff1 - k*pT1;

      //      printf("eff2 %1.3f, eff1 %1.3f, pT2 %1.3f, pT1 %1.3f\n", eff2, eff1, pT2, pT1);
      effTruthExtrap = k * pT + d;
      hRatio[index][iEta]->Fill(pT, eff / effTruthExtrap);
      printf("pT(T&P)=%1.3f corresponds to pTBin %d, eff %1.3f, effTruthExtrap %1.3f\n", 
	     pT, thisBin, eff, effTruthExtrap);
    }
  }
}

//==============================
void DrawEff(Int_t iEta, Bool_t drawTNP, Bool_t drawAllInOne, Int_t effSample){

  Char_t name[100];
  sprintf(name, "c_iEta%d_%d_%d_%d", iEta, drawTNP, drawAllInOne, effSample); 
  TCanvas *c1 = new TCanvas(name);
  TH1F *hFrame1 = gPad->DrawFrame(0., 0., 20.1, 1.1);
  hFrame1->SetXTitle("p_{T} [GeV/c]");
  hFrame1->SetYTitle("single muon efficiency");

  if(drawTNP)
    gEff[TNP][iEta]->Draw("p same");
  if(drawAllInOne)
    gEff[TNPALLINONE][iEta]->Draw("p same");

  tEff[iEta]->Draw("p same");

  TLine *line = new TLine(0., 1., 20.1, 1.);
  line->SetLineStyle(3);
  line->Draw();

  sprintf(name, "%1.1f < |#eta| < %1.1f", etaBins[iEta], etaBins[iEta+1]);
  TLegend *tex1 = new TLegend(0.5847701,0.1758475,0.9353448,0.3834746, name);
  if(effSample == MC){
    if(drawTNP)
      tex1->AddEntry(gEff[TNP][iEta], "factorized MC T&P", "pl");
      // tex1->AddEntry(gEff[TNP][iEta], "product MC T&P", "pl");
    if(drawAllInOne)
      tex1->AddEntry(gEff[TNPALLINONE][iEta], "inclusive MC T&P", "pl");
      // tex1->AddEntry(gEff[TNPALLINONE][iEta], "all in one MC T&P", "pl");
  }
  else if(effSample == DATA){
    if(drawTNP)
      tex1->AddEntry(gEff[TNP][iEta], "factorized Data T&P", "pl");
      //      tex1->AddEntry(gEff[TNP][iEta], "product Data T&P", "pl");
    if(drawAllInOne)
      tex1->AddEntry(gEff[TNPALLINONE][iEta], "inclusive Data T&P", "pl");
      // tex1->AddEntry(gEff[TNPALLINONE][iEta], "all in one Data T&P", "pl");
  }
  tex1->AddEntry(tEff[iEta], "MC truth", "pl");
  tex1->SetFillColor(0); tex1->SetTextSize(0.04);
  tex1->SetBorderSize(0.);
  tex1->Draw(); 

  if(drawTNP && drawAllInOne)
    sprintf(name, "Figures/tnpEff_allInOneEff_truthEff_Eta%d.pdf", iEta);
  else if(drawTNP)
    sprintf(name, "Figures/tnpEff_truthEff_Eta%d.pdf", iEta);
  else if(drawAllInOne)
    sprintf(name, "Figures/allInOneEff_truthEff_Eta%d.pdf", iEta);
  else if(!drawTNP && !drawAllInOne)
    sprintf(name, "Figures/truthEff_Eta%d.pdf", iEta);
  c1->Print(name);
}

//==============================
void LoadEfficiencies(Int_t effSample){

  Char_t name[100];

  //TFile *fTnP = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Ilse/8Dec2011/EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Trk80Cuts_19Nov2011.root");
  //TFile *fTnP = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Ilse/2Feb2012/EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Trk80Cuts_29Jan2012.root");
  //TFile *fTnP = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Ilse/6Feb2012/EfficiencyProduct_combinedMC_Trk80Cuts_6Feb2011.root");
  //TFile *fTnP = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Ilse/10Feb2012/EfficiencyProductDimuon0Jpsi_fineMC_Trk80Cuts_3Feb2012.root");
  //TFile *fTnP = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Linlin/6March2012/EfficiencyProductDimuon0Jpsi_combinedMC_fineBins_Trk80Cuts_6Mar2011.root");//TMOneStationTight moved from MuonQual to MuonID
  //TFile *fTnP = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Linlin/13March2012/EfficiencyProductDimuon0Jpsi_combined_DATA_MC_fineBins_Trk80Cuts_13Mar2012.root");
  TFile *fTnP = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Ilse/5June2012-TnP-WOCuts/FactorizedEfficiency_5June2012.root");
  for(int iEta = 0; iEta < kNbEta; iEta++){
    sprintf(name, "gEff_%s_PT_AETA%d", effSampleName[effSample], iEta);
    gEff[TNP][iEta] = (TGraphAsymmErrors *) fTnP->Get(name);
    gEff[TNP][iEta]->SetMarkerStyle(20);
    gEff[TNP][iEta]->SetMarkerColor(1);
    gEff[TNP][iEta]->SetLineColor(1);
  }
  printf("product of efficiencies read in\n");

  //TFile *fTnPAllInOne = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Linlin/7Dec2011/singleMuonEfficiency_ProbeTrackMatched_data_mc_pt_abseta_tracker80Cuts_7Dec2011.root");
  //TFile *fTnPAllInOne = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Ilse/11Feb2012/singleMuonEff_pt_abseta_seagulls_run1_fineMC_Trk80Cuts_3Feb2012_corrected.root");//very fine pT bins
  //TFile *fTnPAllInOne = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Ilse/6Feb2012/singleMuonEff_combinedMC_Trk80Cuts_6Feb2011.root");
  //TFile *fTnPAllInOne = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Ilse/6Feb2012/singleMuonEff_pt_abseta_seagulls_run1_Trk80Cuts_3Feb2012_corrected.root");
  //TFile *fTnPAllInOne = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Linlin/12March2012/Inclusive_SingleMu_Eff_data_oldMC_Trk80Cuts_12Mar2012.root");
  TFile *fTnPAllInOne = new TFile("/Users/hwoehri/CMS/Work/TnP/2011/Linlin/14March2012/EfficiencyInclusive_SingleMu_combined_DATA_MC_Trk80Cuts_14Mar2012.root");

  for(int iEta = 0; iEta < kNbEta-1; iEta++){
    //sprintf(name, "gEff_%s_PT_ABSETA%d", effSampleName[effSample], iEta);
    sprintf(name, "gEff_%s_PT_AETA%d", effSampleName[effSample], iEta);
    gEff[TNPALLINONE][iEta] = (TGraphAsymmErrors *) fTnPAllInOne->Get(name);
    gEff[TNPALLINONE][iEta]->SetMarkerStyle(24);
    gEff[TNPALLINONE][iEta]->SetMarkerColor(4);
    gEff[TNPALLINONE][iEta]->SetLineColor(4);
  }

  printf("in-one-go efficiencies read in\n");

  //TFile *fTruth = new TFile("singleMuTruthEff_12Jan2012_40GeVrap1_2pT100GeV_New_EtaCut.root");
  //TFile *fTruth = new TFile("singleMuTruthEff_17Jan2012_40GeVrap1_2pT100GeV_EtaCut.root");
  //TFile *fTruth = new TFile("singleMuTruthEff_17Jan2012_40GeVrap1_2pT100GeV_EtaCut_FineBins100MeV.root");
  //TFile *fTruth = new TFile("singleMuTruthEff_17Jan2012_40GeVrap1_2pT100GeV_EtaCut_FineBins100MeVand0p1Eta.root");
  //TFile *fTruth = new TFile("singleMuTruthEff_18Jan2012_40GeVrap1_2pT100GeV_EtaCut_FineBins200MeV.root");
  //TFile *fTruth = new TFile("singleMuTruthEff_7March2012_40GeVrap1_TnPBins.root");
  //TFile *fTruth = new TFile("singleMuTruthEff_22March2012_40GeVrap1_FineBins200MeV_binWeighting.root");
  //TFile *fTruth = new TFile("singleMuTruthEff_24March2012_40GeVrap1_FineBins200MeV_RecoTracks.root");
  //TFile *fTruth = new TFile("singleMuTruthEff_5June2012_40GeVrap1_14PTBins_RecoTracks_noKinemCuts.root");
  //TFile *fTruth = new TFile("singleMuTruthEff_7June2012_40GeVrap1_17PTBins_noKinemCuts_Dimuon0JpsiNoVertexing.root");
  //TFile *fTruth = new TFile("singleMuTruthEff_7June2012_40GeVrap1_17PTBins_RecoTracks_noKinemCuts_Dimuon0JpsiNoVertexing.root");
  TFile *fTruth = new TFile("singleMuTruthEff_7June2012_40GeVrap1_17PTBins_RecoTracks_noKinemCuts_Dimuon0JpsiNoVertexing_2.root");

  for(int iEta = 0; iEta < kNbEta; iEta++){
    sprintf(name, "totEff_MCTRUTH_PT_AETA%d", iEta);
    tEff[iEta] = (TEfficiency *) fTruth->Get(name);
    tEff[iEta]->SetMarkerStyle(21);
    tEff[iEta]->SetMarkerSize(0.4);
    tEff[iEta]->SetMarkerColor(2);
    tEff[iEta]->SetLineColor(2);

    //works only for files from 22/03/2012 onwards...
    sprintf(name, "gtotEff_MCTRUTH_PT_AETA%d", iEta);
    gEff[MCTRUTH][iEta] = (TGraphAsymmErrors *) fTruth->Get(name);
    gEff[MCTRUTH][iEta]->SetMarkerStyle(21);
    gEff[MCTRUTH][iEta]->SetMarkerSize(0.4);
    gEff[MCTRUTH][iEta]->SetMarkerColor(2);
    gEff[MCTRUTH][iEta]->SetLineColor(2);
  }
}
