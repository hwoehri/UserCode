#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"
#include "TEfficiency.h"

//Char_t const *trigName = "HLT_DoubleMu10_Jpsi_Barrel";
Char_t const *trigName = "HLT_DoubleMu5_Upsilon_Barrel";
Int_t const kNbMaxFrame = 3;
// Int_t const kNbEff = 3;
// enum {RECO,TRIG,TOT};
// Char_t *effName[kNbEff] = {"reco", "trig", "tot"};
// Int_t const kNbEff = 1;
// enum {TRIG};
// Char_t *effName[kNbEff] = {"trig"};
Int_t const kNbEff = 1;
enum {TOT};
Char_t *effName[kNbEff] = {"tot"};

Int_t const kNbVar = 3;
Char_t *varName[kNbVar] = {"pT", "y", "phi"};
enum {PT,Y,PHI};
Char_t *effTypeLabel[2] = {"TnP", "Truth"};
enum {TnP, TRUTH};
//MC Truth efficiencies:
TEfficiency *gEffMCTruth1D[kNbEff][kNbVar];
TH1D *hEffMCTruth1D[kNbEff][kNbVar];
enum {ABSETA, ETA};
TEfficiency *gEffMCTruth2D[kNbEff][2];
TH2D *hEffMCTruth2D[kNbEff][2];

//MC T&P efficiencies
TEfficiency *gEffMCTP1D[kNbEff][kNbVar];
TH1D *hEffMCTP1D[kNbEff][kNbVar];
TEfficiency *gEffMCTP2D[kNbEff][2];
TH2D *hEffMCTP2D[kNbEff][2];

//cosTheta and phi maps
TEfficiency *gEff_pol_MCTruth[kNbEff][eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH2D *hEff_pol_MCTruth[kNbEff][eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TEfficiency *gEff_pol_MCTP[kNbEff][eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH2D *hEff_pol_MCTP[kNbEff][eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH2D *hRho_pol[kNbEff][eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];

//1D cosTheta and phi
TEfficiency *gEffMCTP_phiPol[2][kNbEff][eff::kNbFrames];
TEfficiency *gEffMCTP_cosTheta[2][kNbEff][eff::kNbFrames];
TEfficiency *gEffMCTP2D_cosTheta_phiPol[2][kNbEff][eff::kNbFrames];
TEfficiency *gEffMCTP_phiPol_pT_rap[2][kNbEff][eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
TEfficiency *gEffMCTP_cosTheta_pT_rap[2][kNbEff][eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];

//deltaR and deltaPhi vs deltaEta:
TEfficiency *gEff_deltaR_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH1D *hEff_deltaR_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH1D *hEff_deltaR_MCTP[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TEfficiency *gEff_deltaPhiVsDeltaEta_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH2D *hEff_deltaPhiVsDeltaEta_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH2D *hEff_deltaPhiVsDeltaEta_MCTP[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TEfficiency *gEff_deltaRM2_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH1D *hEff_deltaRM2_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH1D *hEff_deltaRM2_MCTP[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TEfficiency *gEff_deltaPhiM2_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH1D *hEff_deltaPhiM2_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH1D *hEff_deltaPhiM2_MCTP[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TEfficiency *gEff_deltaEtaM2_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH1D *hEff_deltaEtaM2_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH1D *hEff_deltaEtaM2_MCTP[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TEfficiency *gEff_distM2_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH1D *hEff_distM2_MCTruth[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH1D *hEff_distM2_MCTP[kNbEff][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];

//helper histograms:
TH1D *hPassed[kNbVar], *hTot[kNbVar];
TH2D *hPassed2D, *hTot2D;
TH1D *hPassed1D, *hTot1D;
TH2D *hPassed2D_2, *hTot2D_2;
TH2D *hPassed2D_pol[kNbEff][eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];
TH2D *hTot2D_pol[kNbEff][eff::kNbFrames][eff::kNbRapForPTBins+1][eff::kNbPTMaxBins+1];

void LoadMCTruthEff(Char_t *fileNameIn, Int_t iEff);
void LoadTPEfficiencies(Char_t *fileNameIn, Int_t iEff);
void LoadPolEfficiencies(Char_t *fileName, Int_t iEff, Int_t iEffType);
void Plot1DEff(Int_t iEff, Int_t iVar);
void PlotRatio1D(Int_t iEff, Int_t iVar);
void Plot2DEff(Int_t iEff, Int_t iRap);
void PlotRatio2D(Int_t iEff, Int_t iRap);
void PlotRhoPol(Int_t iEff, Int_t iFrame);
void PlotEffPol(Int_t iEff, Int_t iFrame);
void PlotEffPol_pT_rap(Int_t iEff, Int_t iFrame, Int_t iRapBin, Int_t iPTBin);
void PlotRhoPol_pT_rap(Int_t iEff, Int_t iFrame, Int_t iRapBin, Int_t iPTBin);
void Plot2DMaps_CosThetaPhi(Int_t iEff, Int_t iFrame, Int_t iRapBin, Int_t iPTBin);
void PlotRatio2D_CosThetaPhi(Int_t iEff, Int_t iFrame, Int_t iRapBin, Int_t iPTBin);
void Plot2DMaps_DeltaPhiDeltaEta(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
void PlotRatio2D_DeltaPhiDeltaEta(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
void PlotMaps_DeltaR(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
void PlotMaps_DeltaRM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
void PlotMaps_DeltaPhiM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
void PlotMaps_DeltaEtaM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
void PlotMaps_DistM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
void PlotRatio_DeltaR(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
void PlotRatio_DeltaRM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
void PlotRatio_DeltaPhiM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
void PlotRatio_DeltaEtaM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
void PlotRatio_DistM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin);
//=======================
//WARNING: when changing from J/psi to Upsilon, adjust pTBinMin and pTBinMax!
void compareEff_CosThetaPhi(Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon5UpsilonBarrel_13March2012_newPTBins.root",
			    Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon5UpsilonBarrel_13March2012_MCTruthEff_FineBins200MeV_newPTBins.root"){
			    // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon5UpsilonBarrel_23Jan2012.root", 
			    // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon5UpsilonBarrel_7Feb2012_ProdSingleMuEff.root"){
			    // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon5UpsilonBarrel_23Jan2012.root", 
			    // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon5UpsilonBarrel_6Feb2012_SingleMuEff.root"){
			    // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_23Jan2012.root", 
			    // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_23Jan2012_MCTruthEff_FineBins200MeV.root"){
			    // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon5UpsilonBarrel_23Jan2012.root", 
			    // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon5UpsilonBarrel_23Jan2012_MCTruthEff_FineBins200MeV.root"){
			    // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon5UpsilonBarrel_22Jan2012.root", 
			    // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon5UpsilonBarrel_22Jan2012_MCTruthEff_FineBins200MeV.root"){
			    // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon5UpsilonBarrel_21Jan2012.root", 
			    // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon5UpsilonBarrel_21Jan2012_MCTruthEff_FineBins200MeV.root"){
			    // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon5UpsilonBarrel_20Jan2012.root", 
			    // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon5UpsilonBarrel_20Jan2012_MCTruthEff_FineBins200MeV.root"){
			    // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon5UpsilonBarrel_18Jan2011.root", 
			    // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon5UpsilonBarrel_18Jan2011_MCTruthEff_FineBins100MeV.root"){
			    // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_18Jan2011.root", 
			    // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_18Jan2011_MCTruthEff_FineBins100MeV.root"){
			    // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_17Jan2011.root", 
			    // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_17Jan2011_MCTruthEff_FineBins100MeVand0p1Eta.root"){
                            //Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_16Dec2011.root", 
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_17Jan2011_MCTruthEff_FineBins100MeV.root"){
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_17Jan2011_MCTruthEff.root"){
                            //Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon5UpsilonBarrel_19Dec2011.root",
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon5UpsilonBarrel_19Dec2011_ProdSingleMuEff.root"){
			    //Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_17Dec2011_VeryTightPSCuts.root")}
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_18Dec2011_ProdSingleMuEff_DataEff.root"){
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_17Dec2011_ProdSingleMuEff_VeryTightPSCuts.root"){
                            //Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_16Dec2011.root",
                            //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_17Dec2011_ProdSingleMuEff_noTrackingIneff.root"){
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_17Dec2011_SingleMuEff.root"){
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_16Dec2011_ProdSingleMuEff.root"){
			    //Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_12Dec2011.root",
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_SingleMuEff_noDimuVtxEffCorr_noJpsiVprobCut_15Dec2011.root"){	
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_ProdSingleMuEff_noDimuVtxEffCorr_noJpsiVprobCut_15Dec2011.root"){	
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_SingleMuEff_DimuVertexEff_13Dec2011_noTrackingInEff.root"){	
                            //Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_12Dec2011.root",
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_4SingleMuEff_DimuVertexEff_12Dec2011_recoEffonly.root"){	
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_4SingleMuEff_DimuVertexEff_12Dec2011_trigEffonly.root"){	
  //  Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_12Dec2011_singleMuEta1p4or0p8.root",
                            //Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_12Dec2011_singleMuEta1p4.root",
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_ProdSingleMu_DimuVertexEff_12Dec2011_noTrackingIneff_singleMuEta1p4.root"){
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_ProdSingleMu_DimuVertexEff_12Dec2011_noTrackingIneff_singleMuEta1p4.root"){
  //  Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_7Dec2011_2.root",
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_SingleMuEff_DimuVertexEff_12Dec2011_noTrackingIneff_centralTimes1c03.root"){//SingleMuEff a la Matt R.
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_ProdSingleMu_DimuVertexEff_12Dec2011_noTrackingIneff_centralTimes1c03.root"){
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_ProdSingleMu_DimuVertexEff_12Dec2011_noTrackingIneff_centralPlus3Perc.root"){
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_ProdSingleMu_DimuVertexEff_12Dec2011_noTrackingIneff.root"){
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_SingleMuEff_DimuVertexEff_12Dec2011_noTrackingIneff.root"){//SingleMuEff a la Matt R.
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_SingleMuEff_DimuVertexEff_12Dec2011.root"){//SingleMuEff a la Matt R.
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_ProdSingleMu_DimuVertexEff_12Dec2011.root"){
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_ProdSingleMu_DimuVertexEff_12Dec2011_centralPlus5Perc.root"){
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_ProdSingleMu_DimuVertexEff_8Dec2011_interpolated2D.root"){
  //			    Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_ProdSingleMu_DimuVertexEff_8Dec2011_interpolated.root"){
			    //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_SingleMuEff_DimuVertexEff_7Dec2011_2.root"){//SingleMuEff a la Matt R.
  //Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_30Nov2011.root",
  //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_ProdSingleMu_6Dec2011.root"){//SingleMuEff a la Matt R.
  // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_30Nov2011.root",
  // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_L1L2L3MuEff_DimuVertexEff_6Dec2011.root"){//SingleMuEff a la Matt R., adding the online dimuon vertex ineff.
  // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_30Nov2011.root",
  // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_L1L2L3MuEff_5Dec2011.root"){//SingleMuEff a la Matt R. (trigger only)
  // Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_30Nov2011.root",
  // Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_30Nov2011.root"){
  //Char_t *fileNameMCTruth = "MCTruthEff_HLTDimuon10JpsiBarrel_12Nov2011.root",
  //Char_t *fileNameTPEff = "MCTnPEff_HLTDimuon10JpsiBarrel_12Nov2011.root"){

  //(1) loading the efficiencies
  printf("loading MCTruth efficiencies\n");
  for(int iEff = 0; iEff < kNbEff; iEff++)
    LoadMCTruthEff(fileNameMCTruth, iEff);

  printf("loading MCTnP efficiencies\n");

  for(int iEff = 0; iEff < kNbEff; iEff++)
    LoadTPEfficiencies(fileNameTPEff, iEff);

  printf("loading polarization related efficiencies\n");

  // Int_t pTBinMin[4] = {6, 6, 6, 10};//Jpsi
  // Int_t pTBinMax[4] = {13, 13, 13, 14};//Jpsi
  // Int_t pTBinMin[4] = {1, 1, 1, 10}; //Ups, until 13 March
  // Int_t pTBinMax[4] = {13, 13, 13, 14}; //Ups, until 13 March
  Int_t pTBinMin[4] = {1, 1, 1, 10}; //Ups
  Int_t pTBinMax[4] = {12, 12, 12, 12}; //Ups
  
  for(int iEff = 0; iEff < kNbEff; iEff++){
  //for(int iEff = 2; iEff < kNbEff; iEff++){
    LoadPolEfficiencies(fileNameMCTruth, iEff, TRUTH);
    LoadPolEfficiencies(fileNameTPEff,   iEff, TnP);
    printf("first pass done\n");

    for(int iFrame = 0; iFrame < kNbMaxFrame; iFrame++){
      PlotEffPol(iEff, iFrame);
      PlotRhoPol(iEff, iFrame);

      for(int iRapBin = 1; iRapBin <= 2; iRapBin++){
      	for(int iPTBin = pTBinMin[iRapBin-1]; iPTBin <= pTBinMax[iRapBin-1]; iPTBin++){
      	  PlotEffPol_pT_rap(iEff, iFrame, iRapBin, iPTBin);
      	  PlotRhoPol_pT_rap(iEff, iFrame, iRapBin, iPTBin);
      	}
      }
    }
  }

  // //(2) plotting the 1D efficiencies: pT, eta, phi
  // for(int iEff = 0; iEff < kNbEff; iEff++)
  //   for(int iVar = 0; iVar < kNbVar; iVar++)
  //     Plot1DEff(iEff, iVar);

  // for(int iEff = 0; iEff < kNbEff; iEff++)
  //   for(int iVar = 0; iVar < kNbVar; iVar++)
  //     PlotRatio1D(iEff, iVar);

  // //(3) plotting the pT vs eta differential efficiencies
  // for(int iEff = 0; iEff < kNbEff; iEff++)
  //   for(int iRap = 0; iRap < 2; iRap++)
  //     Plot2DEff(iEff, iRap);

  // for(int iEff = 0; iEff < kNbEff; iEff++)
  //   for(int iRap = 0; iRap < 2; iRap++)
  //     PlotRatio2D(iEff, iRap);

  // //(4) plot efficiency and ratio of CosTheta vs phi
  // for(int iEff = 0; iEff < kNbEff; iEff++)
  //   for(int iFrame = 0; iFrame < kNbMaxFrame; iFrame++)
  //     for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  //   	//for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // 	for(int iPTBin = 6; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // 	  Plot2DMaps_CosThetaPhi(iEff, iFrame, iRapBin, iPTBin);

  TFile *fOut = new TFile("rhoFactor.root", "RECREATE");
  for(int iEff = 0; iEff < kNbEff; iEff++){
    for(int iFrame = 0; iFrame < kNbMaxFrame; iFrame++){
      for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
	for(int iPTBin = pTBinMin[iRapBin-1]; iPTBin <= pTBinMax[iRapBin-1]; iPTBin++){
  	  PlotRatio2D_CosThetaPhi(iEff, iFrame, iRapBin, iPTBin);
	  if(hRho_pol[iEff][iFrame][iRapBin][iPTBin])
  	    hRho_pol[iEff][iFrame][iRapBin][iPTBin]->Write();
  	}
      }
    }
  }
  fOut->Close();

  // // //(5) plot efficiency and ratio vs deltaEta / deltaPhi
  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	Plot2DMaps_DeltaPhiDeltaEta(iEff, iRapBin, iPTBin);

  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	PlotRatio2D_DeltaPhiDeltaEta(iEff, iRapBin, iPTBin);

  // // //(6) plot efficiency and ratio vs deltaR
  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	PlotMaps_DeltaR(iEff, iRapBin, iPTBin);

  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	PlotRatio_DeltaR(iEff, iRapBin, iPTBin);

  // // //(7) plot efficiency and ratio vs deltaRM2
  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	PlotMaps_DeltaRM2(iEff, iRapBin, iPTBin);

  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	PlotRatio_DeltaRM2(iEff, iRapBin, iPTBin);

  // // //(8) plot efficiency and ratio vs deltaPhiM2
  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	PlotMaps_DeltaPhiM2(iEff, iRapBin, iPTBin);

  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	PlotRatio_DeltaPhiM2(iEff, iRapBin, iPTBin);

  // // //(9) plot efficiency and ratio vs deltaEtaM2
  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	PlotMaps_DeltaEtaM2(iEff, iRapBin, iPTBin);

  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	PlotRatio_DeltaEtaM2(iEff, iRapBin, iPTBin);

  // // //(10) plot efficiency and ratio vs distM2
  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	PlotMaps_DistM2(iEff, iRapBin, iPTBin);

  // // for(int iEff = 0; iEff < kNbEff; iEff++)
  // //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)
  // //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++)
  // // 	PlotRatio_DistM2(iEff, iRapBin, iPTBin);
}


//========================
void PlotRatio1D(Int_t iEff, Int_t iVar){

  Char_t name[100];
  sprintf(name, "c1_%sEff_%s", effName[iEff], varName[iVar]);
  TCanvas *c1 = new TCanvas(name, name);

  Double_t minX[kNbVar] = {10., -1.2, -3.145};
  Double_t maxX[kNbVar] = {50., 1.2, 3.145};

  sprintf(name, "hRatio_%sEff_%s", effName[iEff], varName[iVar]);
  TH1D *hRatio = (TH1D *) hEffMCTruth1D[iEff][iVar]->Clone(name);
  hRatio->Divide(hEffMCTP1D[iEff][iVar]);
  
  hRatio->SetTitle("");
  hRatio->SetMarkerStyle(20);
  hRatio->SetMinimum(0.7);
  hRatio->SetMaximum(1.3);
  hRatio->SetAxisRange(minX[iVar], maxX[iVar]);
  hRatio->Draw("p");

  TLine *line = new TLine(minX[iVar], 1., maxX[iVar], 1.);
  line->SetLineStyle(3); line->Draw();

  sprintf(name, "Figures/ratio_%sEff_%s.pdf", effName[iEff], varName[iVar]);
  c1->Print(name);
}

//=====================================
void Plot1DEff(Int_t iEff, Int_t iVar){

  gStyle->SetOptStat(0);

  Double_t minX[kNbVar] = {10., -1.2, -3.145};
  Double_t maxX[kNbVar] = {30., 1.2, 3.145};

  Char_t name[100];
  sprintf(name, "c2_%sEff_%s", effName[iEff], varName[iVar]);
  TCanvas *c1 = new TCanvas(name, name);

  hEffMCTP1D[iEff][iVar]->SetMinimum(0.);
  hEffMCTP1D[iEff][iVar]->SetMaximum(1.);
  hEffMCTP1D[iEff][iVar]->SetAxisRange(minX[iVar], maxX[iVar]);
  hEffMCTP1D[iEff][iVar]->Draw("p");  
  hEffMCTP1D[iEff][iVar]->SetMarkerStyle(20);
  gEffMCTruth1D[iEff][iVar]->Draw("p same");
  gEffMCTruth1D[iEff][iVar]->SetMarkerStyle(24);

  if(iVar == 0){
    if(iEff == 1)
      sprintf(name, "%s", trigName);
    else
      sprintf(name, "%s", trigName);
    TLegend *leg1 = new TLegend(0.55,0.1610169,0.8,0.3601695, name);
    leg1->AddEntry(gEffMCTruth1D[iEff][iVar], "MC truth", "p");
    leg1->AddEntry(hEffMCTP1D[iEff][iVar], "MC T&P", "p");
    leg1->SetFillColor(0);  leg1->SetBorderSize(0);
    leg1->SetTextSize(0.04);
    leg1->Draw("same");
  }

  sprintf(name, "Figures/eff1D_%sEff_%s.pdf", effName[iEff], varName[iVar]);
  c1->Print(name);
}

//=====================================
void Plot2DEff(Int_t iEff, Int_t iRap){

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("5.2f");
  gStyle->SetPadRightMargin(0.1);

  Char_t name[100];
  //plot the TnP MC histogram
  sprintf(name, "MCTnP_%sEff_%d", effName[iEff], iRap);
  TCanvas *c1 = new TCanvas(name, name);

  hEffMCTP2D[iEff][iRap]->SetMinimum(0.);
  hEffMCTP2D[iEff][iRap]->SetMaximum(1.);
  hEffMCTP2D[iEff][iRap]->SetMarkerSize(1.3);
  hEffMCTP2D[iEff][iRap]->GetYaxis()->SetRangeUser(10., 50.);
  hEffMCTP2D[iEff][iRap]->Draw("colz text e");

  sprintf(name, "Figures/eff2D_TnPMC_%sEff_rap%d.pdf", effName[iEff], iRap);
  c1->Print(name);


  //plot the MC truth histogram
  sprintf(name, "MCTruth_%sEff_%d", effName[iEff], iRap);
  TCanvas *c2 = new TCanvas(name, name);

  //gEffMCTruth2D[iEff][iRap]->Draw("colz text e");
  hEffMCTruth2D[iEff][iRap]->GetYaxis()->SetRangeUser(10., 50.);
  hEffMCTruth2D[iEff][iRap]->Draw("colz text e");

  sprintf(name, "Figures/eff2D_MCTruth_%sEff_rap%d.pdf", effName[iEff], iRap);
  c2->Print(name);
}

//========================
void Plot2DMaps_CosThetaPhi(Int_t iEff, Int_t iFrame, Int_t iRapBin, Int_t iPTBin) {

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("5.2f");
  gStyle->SetPadRightMargin(0.12);

  Char_t name[100];
  //plot MC truth
  sprintf(name, "pol_MCTruth_%sEff_%s_rap%d_pT%d", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  TCanvas *c1 = new TCanvas(name, name, 600, 600);

  if(iPTBin == 0) 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
  else if(iPTBin == eff::kNbPTBins[iRapBin]+1)
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
  else 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);


  gEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->SetTitle(name);
  // gEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->SetMinimum(0.);
  // gEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->SetMaximum(1.);
  gEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->Draw("colz text e");

  sprintf(name, "Figures/%sEff_MCTruth_%s_rap%d_pT%d.pdf", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  c1->Print(name);

  //plot MC T&P
  sprintf(name, "pol_MCTP_%sEff_%s_rap%d_pT%d", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  TCanvas *c2 = new TCanvas(name, name, 600, 600);

  if(iPTBin == 0) 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
  else if(iPTBin == eff::kNbPTBins[iRapBin]+1)//H: 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
  else 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);

  gEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->SetTitle(name);
//   gEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->SetMinimum(0.);
//   gEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->SetMaximum(1.);
  gEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->Draw("colz text e");

  sprintf(name, "Figures/%sEff_MCTP_%s_rap%d_pT%d.pdf", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  c2->Print(name);

}

//========================
void PlotRatio2D_CosThetaPhi(Int_t iEff, Int_t iFrame, Int_t iRapBin, Int_t iPTBin) {

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("5.2f");
  gStyle->SetPadRightMargin(0.12);

  Char_t name[100];
  sprintf(name, "ratio_pol_%sEff_%s_rap%d_pT%d", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  TCanvas *c1 = new TCanvas(name, name, 600, 600);

//   printf("numerator %p\n", hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]);
  sprintf(name, "hRho_pol_%sEff_%s_rap%d_pT%d", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  hRho_pol[iEff][iFrame][iRapBin][iPTBin] = (TH2D *) hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->Clone(name);
//   printf("numerator %p\n", hEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]);
  hRho_pol[iEff][iFrame][iRapBin][iPTBin]->Divide(hEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]);
  hRho_pol[iEff][iFrame][iRapBin][iPTBin]->SetTitle(name);
  hRho_pol[iEff][iFrame][iRapBin][iPTBin]->SetMinimum(0.75);
  hRho_pol[iEff][iFrame][iRapBin][iPTBin]->SetMaximum(1.25);

  if(iPTBin == 0) 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
  else if(iPTBin == eff::kNbPTBins[iRapBin]+1)//H: 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
  else 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);

  hRho_pol[iEff][iFrame][iRapBin][iPTBin]->SetTitle(name);
  hRho_pol[iEff][iFrame][iRapBin][iPTBin]->SetMarkerSize(1.);
  hRho_pol[iEff][iFrame][iRapBin][iPTBin]->Draw("colz text e");

  sprintf(name, "Figures/ratio_pol_%sEff_%s_rap%d_pT%d.pdf", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  c1->Print(name);

}

//========================
void PlotRatio2D(Int_t iEff, Int_t iRap){

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("5.2f");

  Char_t name[100];
  sprintf(name, "ratio2D_%sEff_%d", effName[iEff], iRap);
  TCanvas *c1 = new TCanvas(name, name);

  sprintf(name, "hRatio2D_%sEff_rap%d", effName[iEff], iRap);
  TH2D *hRatio = (TH2D *) hEffMCTruth2D[iEff][iRap]->Clone(name);
  hRatio->Divide(hEffMCTP2D[iEff][iRap]);
  hRatio->SetTitle("");
  hRatio->SetMinimum(0.7);
  hRatio->SetMaximum(1.3);
  hRatio->SetMarkerSize(2.0);

  for(int iBinX = 1; iBinX <= hRatio->GetNbinsX(); iBinX++){
    for(int iBinY = 1; iBinY <= hRatio->GetNbinsY(); iBinY++){

      if(hEffMCTruth2D[iEff][iRap]->GetBinContent(iBinX, iBinY) < 0.1)
	hRatio->SetBinContent(iBinX, iBinY, 0.);
	hRatio->SetBinError(iBinX, iBinY, 0.);
    }
  }

  hRatio->GetYaxis()->SetRangeUser(10., 50.);
  hRatio->Draw("colz text");

  sprintf(name, "Figures/ratio2D_%sEff_rap%d.pdf", effName[iEff], iRap);
  c1->Print(name);
}


//========================
void Plot2DMaps_DeltaPhiDeltaEta(Int_t iEff, Int_t iRapBin, Int_t iPTBin) {

//   gStyle->SetOptStat(0);
//   gStyle->SetPaintTextFormat("5.2f");
//   gStyle->SetPadRightMargin(0.12);

//   Char_t name[100];
//   //plot MC truth
//   sprintf(name, "deltaPhiDeltaEta_MCTruth_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c1 = new TCanvas(name, name, 600, 600);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);


//   hEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin]->SetTitle(name);
//   hEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin]->SetMinimum(0.);
//   hEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin]->SetMaximum(1.);
//   hEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin]->Draw("colz text e");

//   sprintf(name, "Figures/%sEff_deltaPhiVsDeltaEta_MCTruth_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c1->Print(name);

//   //plot MC T&P
//   sprintf(name, "deltaPhiDeltaEta_MCTP_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c2 = new TCanvas(name, name, 600, 600);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)//H: 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);

//   hEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]->SetTitle(name);
// //   hEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]->SetMinimum(0.);
// //   hEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]->SetMaximum(1.);
//   hEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]->Draw("colz text e");

//   sprintf(name, "Figures/%sEff_deltaPhiVsDeltaEta_MCTP_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c2->Print(name);

  }

//========================
void PlotMaps_DeltaR(Int_t iEff, Int_t iRapBin, Int_t iPTBin){

  gStyle->SetOptStat(0);
  // gStyle->SetPaintTextFormat("5.2f");
  // gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadTopMargin(0.1);

  Char_t name[100];
  //plot MC truth
  sprintf(name, "deltaR_MCTruth_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
  TCanvas *c1 = new TCanvas(name, name);

  if(iPTBin == 0) 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
  else if(iPTBin == eff::kNbPTBins[iRapBin]+1)
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
  else 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);
  hEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->SetTitle(name);
  sprintf(name, "#varepsilon(%s)", effName[iEff]);
  hEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->SetYTitle(name);
  hEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->SetMinimum(0.);
  hEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->SetMaximum(1.);
  hEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->SetMarkerStyle(24);
  hEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->Draw("p");

  hEff_deltaR_MCTP[iEff][iRapBin][iPTBin]->SetMarkerStyle(20);
  hEff_deltaR_MCTP[iEff][iRapBin][iPTBin]->Draw("p same");

  TLegend *leg = new TLegend(0.8, 0.8, 0.96, 0.94);
  leg->AddEntry(hEff_deltaR_MCTruth[iEff][iRapBin][iPTBin], "MC Truth", "p");
  leg->AddEntry(hEff_deltaR_MCTP[iEff][iRapBin][iPTBin], "MC T&P", "p");
  leg->SetFillColor(0); leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->Draw();

  sprintf(name, "Figures/%sEff_deltaR_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
  c1->Print(name);

}

//========================
void PlotMaps_DeltaRM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin) {

//   gStyle->SetOptStat(0);
//   // gStyle->SetPaintTextFormat("5.2f");
//   // gStyle->SetPadRightMargin(0.12);
//   gStyle->SetPadTopMargin(0.1);

//   Char_t name[100];
//   //plot MC truth
//   sprintf(name, "deltaRM2_MCTruth_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c1 = new TCanvas(name, name);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);
//   hEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->SetTitle(name);
//   sprintf(name, "#varepsilon(%s)", effName[iEff]);
//   hEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->SetYTitle(name);
//   hEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->SetMinimum(0.);
//   hEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->SetMaximum(1.);
//   hEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->SetMarkerStyle(24);
//   hEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->Draw("p");

//   hEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin]->SetMarkerStyle(20);
//   hEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin]->Draw("p same");

//   TLegend *leg = new TLegend(0.8, 0.8, 0.96, 0.94);
//   leg->AddEntry(hEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin], "MC Truth", "p");
//   leg->AddEntry(hEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin], "MC T&P", "p");
//   leg->SetFillColor(0); leg->SetTextSize(0.04);
//   leg->SetBorderSize(0);
//   leg->Draw();

//   sprintf(name, "Figures/%sEff_deltaRM2_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c1->Print(name);

}

//========================
void PlotMaps_DeltaPhiM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin) {

//   gStyle->SetOptStat(0);
//   // gStyle->SetPaintTextFormat("5.2f");
//   // gStyle->SetPadRightMargin(0.12);
//   gStyle->SetPadTopMargin(0.1);

//   Char_t name[100];
//   //plot MC truth
//   sprintf(name, "deltaPhiM2_MCTruth_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c1 = new TCanvas(name, name);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);
//   hEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->SetTitle(name);
//   sprintf(name, "#varepsilon(%s)", effName[iEff]);
//   hEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->SetYTitle(name);
//   hEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->SetMinimum(0.);
//   hEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->SetMaximum(1.);
//   hEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->SetMarkerStyle(24);
//   hEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->Draw("p");

//   hEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin]->SetMarkerStyle(20);
//   hEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin]->Draw("p same");

//   TLegend *leg = new TLegend(0.8, 0.8, 0.96, 0.94);
//   leg->AddEntry(hEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin], "MC Truth", "p");
//   leg->AddEntry(hEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin], "MC T&P", "p");
//   leg->SetFillColor(0); leg->SetTextSize(0.04);
//   leg->SetBorderSize(0);
//   leg->Draw();

//   sprintf(name, "Figures/%sEff_deltaPhiM2_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c1->Print(name);

}

//========================
void PlotMaps_DeltaEtaM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin) {

//   gStyle->SetOptStat(0);
//   // gStyle->SetPaintTextFormat("5.2f");
//   // gStyle->SetPadRightMargin(0.12);
//   gStyle->SetPadTopMargin(0.1);

//   Char_t name[100];
//   //plot MC truth
//   sprintf(name, "deltaEtaM2_MCTruth_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c1 = new TCanvas(name, name);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);
//   hEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->SetTitle(name);
//   sprintf(name, "#varepsilon(%s)", effName[iEff]);
//   hEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->SetYTitle(name);
//   hEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->SetMinimum(0.);
//   hEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->SetMaximum(1.);
//   hEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->SetMarkerStyle(24);
//   hEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->Draw("p");

//   hEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin]->SetMarkerStyle(20);
//   hEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin]->Draw("p same");

//   TLegend *leg = new TLegend(0.8, 0.8, 0.96, 0.94);
//   leg->AddEntry(hEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin], "MC Truth", "p");
//   leg->AddEntry(hEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin], "MC T&P", "p");
//   leg->SetFillColor(0); leg->SetTextSize(0.04);
//   leg->SetBorderSize(0);
//   leg->Draw();

//   sprintf(name, "Figures/%sEff_deltaEtaM2_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c1->Print(name);

}

//========================
void PlotMaps_DistM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin) {

//   gStyle->SetOptStat(0);
//   // gStyle->SetPaintTextFormat("5.2f");
//   // gStyle->SetPadRightMargin(0.12);
//   gStyle->SetPadTopMargin(0.1);

//   Char_t name[100];
//   //plot MC truth
//   sprintf(name, "distM2_MCTruth_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c1 = new TCanvas(name, name);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);
//   hEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->SetTitle(name);
//   sprintf(name, "#varepsilon(%s)", effName[iEff]);
//   hEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->SetYTitle(name);
//   hEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->SetMinimum(0.);
//   hEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->SetMaximum(1.);
//   hEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->SetMarkerStyle(24);
//   hEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->Draw("p");

//   hEff_distM2_MCTP[iEff][iRapBin][iPTBin]->SetMarkerStyle(20);
//   hEff_distM2_MCTP[iEff][iRapBin][iPTBin]->Draw("p same");

//   TLegend *leg = new TLegend(0.8, 0.8, 0.96, 0.94);
//   leg->AddEntry(hEff_distM2_MCTruth[iEff][iRapBin][iPTBin], "MC Truth", "p");
//   leg->AddEntry(hEff_distM2_MCTP[iEff][iRapBin][iPTBin], "MC T&P", "p");
//   leg->SetFillColor(0); leg->SetTextSize(0.04);
//   leg->SetBorderSize(0);
//   leg->Draw();

//   sprintf(name, "Figures/%sEff_distM2_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c1->Print(name);

}

//========================
void PlotRatio2D_DeltaPhiDeltaEta(Int_t iEff, Int_t iRapBin, Int_t iPTBin) {

//   gStyle->SetOptStat(0);
//   gStyle->SetPaintTextFormat("5.2f");
//   gStyle->SetPadRightMargin(0.12);

//   Char_t name[100];
//   sprintf(name, "ratio_DeltaPhiDeltaEta_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c1 = new TCanvas(name, name, 600, 600);

// //   printf("numerator %p\n", hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]);
//   sprintf(name, "hRatio_deltaPhiVsDeltaEta_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TH2D *hRatio = (TH2D *) hEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin]->Clone(name);

//   hRatio->Divide(hEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]);
//   hRatio->SetTitle(name);
//   hRatio->SetMinimum(0.75);
//   hRatio->SetMaximum(1.25);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)//H: 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);

//   hRatio->SetTitle(name);
//   hRatio->SetMarkerSize(1.);
//   hRatio->Draw("colz text");

// //   for(int iBinX = 1; iBinX <= hRatio->GetNbinsX(); iBinX++){
// //     for(int iBinY = 1; iBinY <= hRatio->GetNbinsY(); iBinY++){

// //       if(hEffMCTruth2D[iEff][iRapBin]->GetBinContent(iBinX, iBinY) < 0.05)
// // 	hRatio->SetBinContent(iBinX, iBinY, 0.);
// // 	hRatio->SetBinError(iBinX, iBinY, 0.);
// //     }
// //   }

//   sprintf(name, "Figures/ratio_deltaPhiVsDeltaEta_%sEff_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c1->Print(name);

}

//========================
void PlotRatio_DeltaR(Int_t iEff, Int_t iRapBin, Int_t iPTBin) {

//   gStyle->SetOptStat(0);
//   // gStyle->SetPaintTextFormat("5.2f");
//   // gStyle->SetPadRightMargin(0.12);

//   Char_t name[100];
//   sprintf(name, "ratio_DeltaR_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c1 = new TCanvas(name, name);

// //   printf("numerator %p\n", hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]);
//   sprintf(name, "hRatio_deltaR_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TH2D *hRatio = (TH2D *) hEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->Clone(name);
//   hRatio->Divide(hEff_deltaR_MCTP[iEff][iRapBin][iPTBin]);
//   hRatio->SetTitle(name);
//   hRatio->SetMinimum(0.5);
//   hRatio->SetMaximum(1.5);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)//H: 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);

//   hRatio->SetTitle(name);
//   hRatio->SetYTitle("#rho = MCTruth / MC T&P");
//   hRatio->SetMarkerSize(1.);
//   hRatio->Draw("p");

//   TLine *line = new TLine(0., 1., 3.14, 1.);
//   line->SetLineStyle(3); line->Draw();

//   sprintf(name, "Figures/ratio_deltaR_%sEff_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c1->Print(name);

}

//========================
void PlotRatio_DeltaRM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin){

//   gStyle->SetOptStat(0);
//   // gStyle->SetPaintTextFormat("5.2f");
//   // gStyle->SetPadRightMargin(0.12);

//   Char_t name[100];
//   sprintf(name, "ratio_DeltaRM2_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c1 = new TCanvas(name, name);

// //   printf("numerator %p\n", hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]);
//   sprintf(name, "hRatio_deltaRM2_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TH2D *hRatio = (TH2D *) hEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->Clone(name);
//   hRatio->Divide(hEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin]);
//   hRatio->SetTitle(name);
//   hRatio->SetMinimum(0.5);
//   hRatio->SetMaximum(1.5);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)//H: 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);

//   hRatio->SetTitle(name);
//   hRatio->SetYTitle("#rho = MCTruth / MC T&P");
//   hRatio->SetMarkerSize(1.);
//   hRatio->Draw("p");

//   TLine *line = new TLine(0., 1., 3.14, 1.);
//   line->SetLineStyle(3); line->Draw();

//   sprintf(name, "Figures/ratio_deltaRM2_%sEff_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c1->Print(name);

}

//========================
void PlotRatio_DeltaPhiM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin) {

//   gStyle->SetOptStat(0);
//   // gStyle->SetPaintTextFormat("5.2f");
//   // gStyle->SetPadRightMargin(0.12);

//   Char_t name[100];
//   sprintf(name, "ratio_DeltaPhiM2_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c1 = new TCanvas(name, name);

// //   printf("numerator %p\n", hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]);
//   sprintf(name, "hRatio_deltaPhiM2_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TH2D *hRatio = (TH2D *) hEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->Clone(name);
//   hRatio->Divide(hEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin]);
//   hRatio->SetTitle(name);
//   hRatio->SetMinimum(0.5);
//   hRatio->SetMaximum(1.5);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)//H: 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);

//   hRatio->SetTitle(name);
//   hRatio->SetYTitle("#rho = MCTruth / MC T&P");
//   hRatio->SetMarkerSize(1.);
//   hRatio->Draw("p");

//   TLine *line = new TLine(0., 1., 180., 1.);
//   line->SetLineStyle(3); line->Draw();

//   sprintf(name, "Figures/ratio_deltaPhiM2_%sEff_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c1->Print(name);

}

//========================
void PlotRatio_DeltaEtaM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin) {

//   gStyle->SetOptStat(0);
//   // gStyle->SetPaintTextFormat("5.2f");
//   // gStyle->SetPadRightMargin(0.12);

//   Char_t name[100];
//   sprintf(name, "ratio_DeltaEtaM2_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c1 = new TCanvas(name, name);

// //   printf("numerator %p\n", hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]);
//   sprintf(name, "hRatio_deltaEtaM2_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TH2D *hRatio = (TH2D *) hEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->Clone(name);
//   hRatio->Divide(hEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin]);
//   hRatio->SetTitle(name);
//   hRatio->SetMinimum(0.5);
//   hRatio->SetMaximum(1.5);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)//H: 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);

//   hRatio->SetTitle(name);
//   hRatio->SetYTitle("#rho = MCTruth / MC T&P");
//   hRatio->SetMarkerSize(1.);
//   hRatio->Draw("p");

//   TLine *line = new TLine(0., 1., 2.5, 1.);
//   line->SetLineStyle(3); line->Draw();

//   sprintf(name, "Figures/ratio_deltaEtaM2_%sEff_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c1->Print(name);

}

//========================
void PlotRatio_DistM2(Int_t iEff, Int_t iRapBin, Int_t iPTBin) {

//   gStyle->SetOptStat(0);
//   // gStyle->SetPaintTextFormat("5.2f");
//   // gStyle->SetPadRightMargin(0.12);

//   Char_t name[100];
//   sprintf(name, "ratio_DistM2_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TCanvas *c1 = new TCanvas(name, name);

// //   printf("numerator %p\n", hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]);
//   sprintf(name, "hRatio_distM2_%sEff_rap%d_pT%d", effName[iEff], iRapBin, iPTBin);
//   TH2D *hRatio = (TH2D *) hEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->Clone(name);
//   hRatio->Divide(hEff_distM2_MCTP[iEff][iRapBin][iPTBin]);
//   hRatio->SetTitle(name);
//   hRatio->SetMinimum(0.5);
//   hRatio->SetMaximum(1.5);

//   if(iPTBin == 0) 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T}", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin]);
//   else if(iPTBin == eff::kNbPTBins[iRapBin]+1)//H: 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c\n", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1]);
//   else 
//     sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c", eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);

//   hRatio->SetTitle(name);
//   hRatio->SetYTitle("#rho = MCTruth / MC T&P");
//   hRatio->SetMarkerSize(1.);
//   hRatio->Draw("p");

//   TLine *line = new TLine(0., 1., 200., 1.);
//   line->SetLineStyle(3); line->Draw();

//   sprintf(name, "Figures/ratio_distM2_%sEff_rap%d_pT%d.pdf", effName[iEff], iRapBin, iPTBin);
//   c1->Print(name);

}

//========================
void LoadTPEfficiencies(Char_t *fileNameIn, Int_t iEff){

  TFile *fIn = new TFile(fileNameIn);
  printf("reading file %s for %s eff\n", fileNameIn, effName[iEff]);

  Char_t name[100];
  for(int iVar = 0; iVar < kNbVar; iVar++){
    sprintf(name, "%sEff_%s", effName[iEff], varName[iVar]);
    printf("taking %s ...\t", name);
    gEffMCTP1D[iEff][iVar] = (TEfficiency *) gDirectory->Get(name);
    printf("%p\n", gEffMCTP1D[iEff][iVar]);
    sprintf(name, "%sEffMCTP_%s", effName[iEff], varName[iVar]);
    gEffMCTP1D[iEff][iVar]->SetName(name);

    //copy the values into a histogram:
    hPassed[iVar] = (TH1D *) gEffMCTP1D[iEff][iVar]->GetPassedHistogram();
    hPassed[iVar]->Sumw2(); 
    printf("%s has %d passed entries\n", name, hPassed[iVar]->GetEntries());
    hTot[iVar] = (TH1D *) gEffMCTP1D[iEff][iVar]->GetTotalHistogram();
    hTot[iVar]->Sumw2();
    sprintf(name, "h%sEffMCTP_%s", effName[iEff], varName[iVar]);
    hEffMCTP1D[iEff][iVar] = (TH1D *) hPassed[iVar]->Clone(name);
    hEffMCTP1D[iEff][iVar]->Divide(hPassed[iVar], hTot[iVar], 1., 1., "B");
    hEffMCTP1D[iEff][iVar]->SetTitle("");
  }

  sprintf(name, "%sEff2D_pT_rap", effName[iEff]);
  printf("taking %s ...\t", name);
  gEffMCTP2D[iEff][ABSETA] = (TEfficiency *) gDirectory->Get(name);
  printf("%p\n", gEffMCTP2D[iEff][ABSETA]);
  sprintf(name, "%sEffMCTP2D_pT_rap", effName[iEff]);
  gEffMCTP2D[iEff][ABSETA]->SetName(name);
  //copy the values into a histogram:
  // TH2D *hPassed2D, *hTot2D;
  hPassed2D = (TH2D *) gEffMCTP2D[iEff][ABSETA]->GetPassedHistogram();
  hPassed2D->Sumw2();
  printf("%s has %d passed entries\n", name, hPassed2D->GetEntries());
  hTot2D = (TH2D *) gEffMCTP2D[iEff][ABSETA]->GetTotalHistogram();
  hTot2D->Sumw2();
  sprintf(name, "h%sEffMCTP2D_pT_rap", effName[iEff]);
  hEffMCTP2D[iEff][ABSETA] = (TH2D *) hPassed2D->Clone(name);
  hEffMCTP2D[iEff][ABSETA]->Divide(hPassed2D, hTot2D, 1., 1., "B");
  hEffMCTP2D[iEff][ABSETA]->SetTitle("");
  hEffMCTP2D[iEff][ABSETA]->SetMaximum(1.);

  sprintf(name, "%sEff2D_pT_rapNP", effName[iEff]);
  gEffMCTP2D[iEff][ETA] = (TEfficiency *) gDirectory->Get(name);
  sprintf(name, "%sEffMCTP2D_pT_rapNP", effName[iEff]);
  gEffMCTP2D[iEff][ETA]->SetName(name);

  // TH2D *hPassed2D_2, *hTot2D_2;
  hPassed2D_2 = (TH2D *) gEffMCTP2D[iEff][ETA]->GetPassedHistogram();
  hPassed2D_2->Sumw2();
  printf("%s has %d passed entries\n", name, hPassed2D_2->GetEntries());
  hTot2D_2 = (TH2D *) gEffMCTP2D[iEff][ETA]->GetTotalHistogram();
  hTot2D_2->Sumw2();
  sprintf(name, "h%sEffMCTP2D_pT_rapNP", effName[iEff]);
  hEffMCTP2D[iEff][ETA] = (TH2D *) hPassed2D_2->Clone(name);
  hEffMCTP2D[iEff][ETA]->Divide(hPassed2D_2, hTot2D_2, 1., 1., "B");
  hEffMCTP2D[iEff][ETA]->SetTitle("");
  hEffMCTP2D[iEff][ETA]->SetMaximum(1.);

  //load the cosTheta and phi 2D maps:
  for(int iFrame = 0; iFrame < kNbMaxFrame; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

  	// //4-folding in phi
        // sprintf(name, "%sEff2D_Onia_phiFolded_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
  	// gEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  	// sprintf(name, "%sEffMCTP_cosThetaPhi_phiFolded_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
  	// gEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->SetName(name);

  	//no folding in phi and cosTheta
        sprintf(name, "%sEff2D_Onia_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
  	gEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  	sprintf(name, "%sEffMCTP_cosThetaPhi_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
  	gEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->SetName(name);

  	hPassed2D_pol[iEff][iFrame][iRapBin][iPTBin] = (TH2D *) gEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->GetPassedHistogram();
  	hPassed2D_pol[iEff][iFrame][iRapBin][iPTBin]->Sumw2();
	printf("%s has %d passed entries\n", name, hPassed2D_pol[iEff][iFrame][iRapBin][iPTBin]->GetEntries());
  	hTot2D_pol[iEff][iFrame][iRapBin][iPTBin] = (TH2D *) gEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->GetTotalHistogram();
  	hTot2D_pol[iEff][iFrame][iRapBin][iPTBin]->Sumw2();
  	sprintf(name, "h%sEffMCTP_cosThetaPhi_phiFolded_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
  	hEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin] = (TH2D *) hPassed2D_pol[iEff][iFrame][iRapBin][iPTBin]->Clone(name);
  	hEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->Divide(hPassed2D_pol[iEff][iFrame][iRapBin][iPTBin], hTot2D_pol[iEff][iFrame][iRapBin][iPTBin], 1., 1., "B");
  	hEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->SetTitle("");
  	hEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->SetMaximum(1.);

  	// printf("eff %s, frame %s, rap %d, pT %d has %1.0f and %1.0f entries\n", 
  	//        effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin, 
  	//        hPassed2D_pol[iEff][iFrame][iRapBin][iPTBin]->GetEntries(),
  	//        hTot2D_pol[iEff][iFrame][iRapBin][iPTBin]->GetEntries());
      }
    }
  }

  // for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
      
  //     //deltaPhiVsDeltaEta cut:
  //     sprintf(name, "%sEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_deltaPhiVsDeltaEta_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed2D = (TH2D *) gEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed2D->Sumw2();
  //     hTot2D = (TH2D *) gEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot2D->Sumw2();
  //     sprintf(name, "h%sEffMCTP_deltaPhiVsDeltaEta_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin] = (TH2D *) hPassed2D->Clone(name);
  //     hEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]->Divide(hPassed2D, hTot2D, 1., 1., "B");
  //     hEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]->SetTitle("");
  //     hEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]->SetMaximum(1.);

  //     //deltaR cut:
  //     sprintf(name, "%sEff_deltaR_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaR_MCTP[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_deltaR_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaR_MCTP[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed1D = (TH1D *) gEff_deltaR_MCTP[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed1D->Sumw2();
  //     hTot1D = (TH1D *) gEff_deltaR_MCTP[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot1D->Sumw2();
  //     sprintf(name, "h%sEffMCTP_deltaR_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaR_MCTP[iEff][iRapBin][iPTBin] = (TH1D *) hPassed1D->Clone(name);
  //     hEff_deltaR_MCTP[iEff][iRapBin][iPTBin]->Divide(hPassed1D, hTot1D, 1., 1., "B");
  //     hEff_deltaR_MCTP[iEff][iRapBin][iPTBin]->SetTitle("");

  //     //deltaRM2 cut:
  //     sprintf(name, "%sEff_deltaRM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_deltaRM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed1D = (TH1D *) gEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed1D->Sumw2();
  //     hTot1D = (TH1D *) gEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot1D->Sumw2();
  //     sprintf(name, "h%sEffMCTP_deltaRM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin] = (TH1D *) hPassed1D->Clone(name);
  //     hEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin]->Divide(hPassed1D, hTot1D, 1., 1., "B");
  //     hEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin]->SetTitle("");

  //     //deltaPhiM2 cut:
  //     sprintf(name, "%sEff_deltaPhiM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_deltaPhiM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed1D = (TH1D *) gEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed1D->Sumw2();
  //     hTot1D = (TH1D *) gEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot1D->Sumw2();
  //     sprintf(name, "h%sEffMCTP_deltaPhiM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin] = (TH1D *) hPassed1D->Clone(name);
  //     hEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin]->Divide(hPassed1D, hTot1D, 1., 1., "B");
  //     hEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin]->SetTitle("");

  //     //deltaEtaM2 cut:
  //     sprintf(name, "%sEff_deltaEtaM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_deltaEtaM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed1D = (TH1D *) gEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed1D->Sumw2();
  //     hTot1D = (TH1D *) gEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot1D->Sumw2();
  //     sprintf(name, "h%sEffMCTP_deltaEtaM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin] = (TH1D *) hPassed1D->Clone(name);
  //     hEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin]->Divide(hPassed1D, hTot1D, 1., 1., "B");
  //     hEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin]->SetTitle("");

  //     //distM2 cut:
  //     sprintf(name, "%sEff_distM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_distM2_MCTP[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_distM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_distM2_MCTP[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed1D = (TH1D *) gEff_distM2_MCTP[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed1D->Sumw2();
  //     hTot1D = (TH1D *) gEff_distM2_MCTP[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot1D->Sumw2();
  //     sprintf(name, "h%sEffMCTP_distM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_distM2_MCTP[iEff][iRapBin][iPTBin] = (TH1D *) hPassed1D->Clone(name);
  //     hEff_distM2_MCTP[iEff][iRapBin][iPTBin]->Divide(hPassed1D, hTot1D, 1., 1., "B");
  //     hEff_distM2_MCTP[iEff][iRapBin][iPTBin]->SetTitle("");
  //   }
  // }







  //old version using histograms

  // TFile *fIn = new TFile(fileNameIn);
  // Char_t name[100];
  // for(int iVar = 0; iVar < kNbVar; iVar++){
  //   sprintf(name, "%sEff_%s", effName[iEff], varName[iVar]);
  //   hEffTP[iEff][iVar] = (TH1D *) gDirectory->Get(name);
  //   sprintf(name, "%sEffMCTP_%s", effName[iEff], varName[iVar]);
  //   hEffTP[iEff][iVar]->SetName(name);
  // }

  // sprintf(name, "%sEff2D_pT_rap", effName[iEff]);
  // hEffTP2D[iEff][ABSETA] = (TH2D *) gDirectory->Get(name);
  // sprintf(name, "%sEffMCTP2D_pT_rap", effName[iEff]);
  // hEffTP2D[iEff][ABSETA]->SetName(name);

  // sprintf(name, "%sEff2D_pT_rapNP", effName[iEff]);
  // hEffTP2D[iEff][ETA] = (TH2D *) gDirectory->Get(name);
  // sprintf(name, "%sEffMCTP2D_pT_rapNP", effName[iEff]);
  // hEffTP2D[iEff][ETA]->SetName(name);

  // //load the cosTheta and phi 2D maps:
  // for(int iFrame = 0; iFrame < kNbMaxFrame; iFrame++){
  //   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
  //     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

  // 	//4-folding in phi                                                                                                        
  //       sprintf(name, "%sEff2D_Onia_phiFolded_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
  // 	hEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin] = (TH2D *) gDirectory->Get(name);
  // 	sprintf(name, "%sEffMCTP_cosThetaPhi_phiFolded_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
  // 	hEff_pol_MCTP[iEff][iFrame][iRapBin][iPTBin]->SetName(name);
  //     }
  //   }
  // }

  // for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
  //     //deltaPhiVsDeltaEta
  //     sprintf(name, "%sEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin] = (TH2D *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_deltaPhiVsDeltaEta_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaPhiVsDeltaEta_MCTP[iEff][iRapBin][iPTBin]->SetName(name);
  //     //deltaR
  //     sprintf(name, "%sEff_deltaR_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaR_MCTP[iEff][iRapBin][iPTBin] = (TH1D *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_deltaR_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaR_MCTP[iEff][iRapBin][iPTBin]->SetName(name);
  //     //deltaRM2
  //     sprintf(name, "%sEff_deltaRM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin] = (TH1D *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_deltaRM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaRM2_MCTP[iEff][iRapBin][iPTBin]->SetName(name);
  //     //deltaPhiM2
  //     sprintf(name, "%sEff_deltaPhiM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin] = (TH1D *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_deltaPhiM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaPhiM2_MCTP[iEff][iRapBin][iPTBin]->SetName(name);
  //     //deltaEtaM2
  //     sprintf(name, "%sEff_deltaEtaM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin] = (TH1D *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_deltaEtaM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaEtaM2_MCTP[iEff][iRapBin][iPTBin]->SetName(name);
  //     //distM2
  //     sprintf(name, "%sEff_distM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_distM2_MCTP[iEff][iRapBin][iPTBin] = (TH1D *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTP_distM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_distM2_MCTP[iEff][iRapBin][iPTBin]->SetName(name);
  //   }
  // }
}

//========================
void LoadMCTruthEff(Char_t *fileNameIn, Int_t iEff){

  TFile *fIn = new TFile(fileNameIn);
  printf("reading file %s for %s eff\n", fileNameIn, effName[iEff]);

  Char_t name[100];
  // TH1D *hPassed[kNbVar], *hTot[kNbVar];
  for(int iVar = 0; iVar < kNbVar; iVar++){
    sprintf(name, "%sEff_%s", effName[iEff], varName[iVar]);
    gEffMCTruth1D[iEff][iVar] = (TEfficiency *) gDirectory->Get(name);
    sprintf(name, "%sEffMCTruth_%s", effName[iEff], varName[iVar]);
    gEffMCTruth1D[iEff][iVar]->SetName(name);

    //copy the values into a histogram:
    hPassed[iVar] = (TH1D *) gEffMCTruth1D[iEff][iVar]->GetPassedHistogram();
    hPassed[iVar]->Sumw2(); 
    hTot[iVar] = (TH1D *) gEffMCTruth1D[iEff][iVar]->GetTotalHistogram();
    hTot[iVar]->Sumw2();
    sprintf(name, "h%sEffMCTruth_%s", effName[iEff], varName[iVar]);
    hEffMCTruth1D[iEff][iVar] = (TH1D *) hPassed[iVar]->Clone(name);
    hEffMCTruth1D[iEff][iVar]->Divide(hPassed[iVar], hTot[iVar], 1., 1., "B");
    hEffMCTruth1D[iEff][iVar]->SetTitle("");
  }
  printf("1D histos done...\n");

  sprintf(name, "%sEff2D_pT_rap", effName[iEff]);
  gEffMCTruth2D[iEff][ABSETA] = (TEfficiency *) gDirectory->Get(name);
  sprintf(name, "%sEffMCTruth2D_pT_rap", effName[iEff]);
  gEffMCTruth2D[iEff][ABSETA]->SetName(name);

  sprintf(name, "%sEff2D_pT_rapNP", effName[iEff]);
  gEffMCTruth2D[iEff][ETA] = (TEfficiency *) gDirectory->Get(name);
  sprintf(name, "%sEffMCTruth2D_pT_rapNP", effName[iEff]);
  gEffMCTruth2D[iEff][ETA]->SetName(name);

  //copy the values into a histogram:
  // TH2D *hPassed2D, *hTot2D;
  hPassed2D = (TH2D *) gEffMCTruth2D[iEff][ABSETA]->GetPassedHistogram();
  hPassed2D->Sumw2();
  hTot2D = (TH2D *) gEffMCTruth2D[iEff][ABSETA]->GetTotalHistogram();
  hTot2D->Sumw2();
  sprintf(name, "h%sEffMCTruth2D_pT_rap", effName[iEff]);
  hEffMCTruth2D[iEff][ABSETA] = (TH2D *) hPassed2D->Clone(name);
  hEffMCTruth2D[iEff][ABSETA]->Divide(hPassed2D, hTot2D, 1., 1., "B");
  hEffMCTruth2D[iEff][ABSETA]->SetTitle("");
  hEffMCTruth2D[iEff][ABSETA]->SetMaximum(1.);

  // TH2D *hPassed2D_2, *hTot2D_2;
  hPassed2D_2 = (TH2D *) gEffMCTruth2D[iEff][ETA]->GetPassedHistogram();
  hPassed2D_2->Sumw2();
  hTot2D_2 = (TH2D *) gEffMCTruth2D[iEff][ETA]->GetTotalHistogram();
  hTot2D_2->Sumw2();
  sprintf(name, "h%sEffMCTruth2D_pT_rapNP", effName[iEff]);
  hEffMCTruth2D[iEff][ETA] = (TH2D *) hPassed2D_2->Clone(name);
  hEffMCTruth2D[iEff][ETA]->Divide(hPassed2D_2, hTot2D_2, 1., 1., "B");
  hEffMCTruth2D[iEff][ETA]->SetTitle("");
  hEffMCTruth2D[iEff][ETA]->SetMaximum(1.);

  printf("2D histos done... \n");

  // TH2D *hPassed2D_pol[kNbEff][eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
  // TH2D *hTot2D_pol[kNbEff][eff::kNbFrames][eff::kNbPTMaxBins+1][eff::kNbRapForPTBins+1];
  //load the cosTheta and phi 2D maps:
  for(int iFrame = 0; iFrame < kNbMaxFrame; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

	// //4-folding in phi
        // sprintf(name, "%sEff2D_Onia_phiFolded_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
	// gEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
	// sprintf(name, "%sEffMCTruth_cosThetaPhi_phiFolded_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
	// gEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->SetName(name);

 	//no folding in phi and cosTheta
        sprintf(name, "%sEff2D_Onia_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
	gEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
	printf("got %s\n", name);
	sprintf(name, "%sEffMCTruth_cosThetaPhi_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
	gEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->SetName(name);

	hPassed2D_pol[iEff][iFrame][iRapBin][iPTBin] = (TH2D *) gEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->GetPassedHistogram();
	hPassed2D_pol[iEff][iFrame][iRapBin][iPTBin]->Sumw2();
	hTot2D_pol[iEff][iFrame][iRapBin][iPTBin] = (TH2D *) gEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->GetTotalHistogram();
	hTot2D_pol[iEff][iFrame][iRapBin][iPTBin]->Sumw2();
	sprintf(name, "h%sEffMCTruth_cosThetaPhi_phiFolded_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
	hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin] = (TH2D *) hPassed2D_pol[iEff][iFrame][iRapBin][iPTBin]->Clone(name);
	hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->Divide(hPassed2D_pol[iEff][iFrame][iRapBin][iPTBin], hTot2D_pol[iEff][iFrame][iRapBin][iPTBin], 1., 1., "B");
	hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->SetTitle("");
	hEff_pol_MCTruth[iEff][iFrame][iRapBin][iPTBin]->SetMaximum(1.);

	// printf("eff %s, frame %s, rap %d, pT %d has %1.0f and %1.0f entries\n", 
	//        effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin, 
	//        hPassed2D_pol[iEff][iFrame][iRapBin][iPTBin]->GetEntries(),
	//        hTot2D_pol[iEff][iFrame][iRapBin][iPTBin]->GetEntries());
      }
    }
  }

  printf("2D histos for different frames, rap and pT bins done\n");

  //for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++)// {
  //   for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
      
  //     //deltaPhiVsDeltaEta cut:
  //     sprintf(name, "%sEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTruth_deltaPhiVsDeltaEta_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed2D = (TH2D *) gEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed2D->Sumw2();
  //     hTot2D = (TH2D *) gEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot2D->Sumw2();
  //     sprintf(name, "h%sEffMCTruth_deltaPhiVsDeltaEta_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin] = (TH2D *) hPassed2D->Clone(name);
  //     hEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin]->Divide(hPassed2D, hTot2D, 1., 1., "B");
  //     hEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin]->SetTitle("");
  //     hEff_deltaPhiVsDeltaEta_MCTruth[iEff][iRapBin][iPTBin]->SetMaximum(1.);

  //     //deltaR cut:
  //     sprintf(name, "%sEff_deltaR_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaR_MCTruth[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTruth_deltaR_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed1D = (TH1D *) gEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed1D->Sumw2();
  //     hTot1D = (TH1D *) gEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot1D->Sumw2();
  //     sprintf(name, "h%sEffMCTruth_deltaR_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaR_MCTruth[iEff][iRapBin][iPTBin] = (TH1D *) hPassed1D->Clone(name);
  //     hEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->Divide(hPassed1D, hTot1D, 1., 1., "B");
  //     hEff_deltaR_MCTruth[iEff][iRapBin][iPTBin]->SetTitle("");

  //     //deltaRM2 cut:
  //     sprintf(name, "%sEff_deltaRM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTruth_deltaRM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed1D = (TH1D *) gEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed1D->Sumw2();
  //     hTot1D = (TH1D *) gEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot1D->Sumw2();
  //     sprintf(name, "h%sEffMCTruth_deltaRM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin] = (TH1D *) hPassed1D->Clone(name);
  //     hEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->Divide(hPassed1D, hTot1D, 1., 1., "B");
  //     hEff_deltaRM2_MCTruth[iEff][iRapBin][iPTBin]->SetTitle("");

  //     //deltaPhiM2 cut:
  //     sprintf(name, "%sEff_deltaPhiM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTruth_deltaPhiM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed1D = (TH1D *) gEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed1D->Sumw2();
  //     hTot1D = (TH1D *) gEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot1D->Sumw2();
  //     sprintf(name, "h%sEffMCTruth_deltaPhiM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin] = (TH1D *) hPassed1D->Clone(name);
  //     hEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->Divide(hPassed1D, hTot1D, 1., 1., "B");
  //     hEff_deltaPhiM2_MCTruth[iEff][iRapBin][iPTBin]->SetTitle("");

  //     //deltaEtaM2 cut:
  //     sprintf(name, "%sEff_deltaEtaM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTruth_deltaEtaM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed1D = (TH1D *) gEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed1D->Sumw2();
  //     hTot1D = (TH1D *) gEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot1D->Sumw2();
  //     sprintf(name, "h%sEffMCTruth_deltaEtaM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin] = (TH1D *) hPassed1D->Clone(name);
  //     hEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->Divide(hPassed1D, hTot1D, 1., 1., "B");
  //     hEff_deltaEtaM2_MCTruth[iEff][iRapBin][iPTBin]->SetTitle("");

  //     //distM2 cut:
  //     sprintf(name, "%sEff_distM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_distM2_MCTruth[iEff][iRapBin][iPTBin] = (TEfficiency *) gDirectory->Get(name);
  //     sprintf(name, "%sEffMCTruth_distM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     gEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->SetName(name);

  //     hPassed1D = (TH1D *) gEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->GetPassedHistogram();
  //     hPassed1D->Sumw2();
  //     hTot1D = (TH1D *) gEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->GetTotalHistogram();
  //     hTot1D->Sumw2();
  //     sprintf(name, "h%sEffMCTruth_distM2_pT%d_rap%d", effName[iEff], iPTBin, iRapBin);
  //     hEff_distM2_MCTruth[iEff][iRapBin][iPTBin] = (TH1D *) hPassed1D->Clone(name);
  //     hEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->Divide(hPassed1D, hTot1D, 1., 1., "B");
  //     hEff_distM2_MCTruth[iEff][iRapBin][iPTBin]->SetTitle("");
  //   }
  // }
}

//===========================================
void LoadPolEfficiencies(Char_t *fileName, Int_t iEff, Int_t iEffType){

  TFile *f = new TFile(fileName);

  Char_t name[100];
  for(int iFrame = 0; iFrame < kNbMaxFrame; iFrame++){
    //1D versus phi
    sprintf(name, "%sEff_phiPol_%s", effName[iEff], eff::frameLabel[iFrame]);
    gEffMCTP_phiPol[iEffType][iEff][iFrame] = (TEfficiency *) gDirectory->Get(name);
    sprintf(name, "%sEff_phiPol_%s_%s", effName[iEff], eff::frameLabel[iFrame], effTypeLabel[iEffType]);
    gEffMCTP_phiPol[iEffType][iEff][iFrame]->SetName(name);
    //1D versus cosTheta
    sprintf(name, "%sEff_cosTheta_%s", effName[iEff], eff::frameLabel[iFrame]);
    gEffMCTP_cosTheta[iEffType][iEff][iFrame] = (TEfficiency *) gDirectory->Get(name);
    sprintf(name, "%sEff_cosTheta_%s_%s", effName[iEff], eff::frameLabel[iFrame], effTypeLabel[iEffType]);
    gEffMCTP_cosTheta[iEffType][iEff][iFrame]->SetName(name);
    printf("%s has %1.0f / %1.0f entries\n", gEffMCTP_cosTheta[iEffType][iEff][iFrame]->GetName(),
	   gEffMCTP_cosTheta[iEffType][iEff][iFrame]->GetPassedHistogram()->GetEntries(),
	   gEffMCTP_cosTheta[iEffType][iEff][iFrame]->GetTotalHistogram()->GetEntries());

    //2D versus cosTheta
    sprintf(name, "%sEff_phiVsCosTheta_%s", effName[iEff], eff::frameLabel[iFrame]);
    gEffMCTP2D_cosTheta_phiPol[iEffType][iEff][iFrame] = (TEfficiency *) gDirectory->Get(name);
    sprintf(name, "%sEff_phiVsCosTheta_%s_%s", effName[iEff], eff::frameLabel[iFrame], effTypeLabel[iEffType]);
    gEffMCTP2D_cosTheta_phiPol[iEffType][iEff][iFrame]->SetName(name);

    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

    	sprintf(name, "%sEff_phiPol_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
    	gEffMCTP_phiPol_pT_rap[iEffType][iEff][iFrame][iPTBin][iRapBin] = (TEfficiency *) gDirectory->Get(name);
    	sprintf(name, "%sEff_phiPol_%s_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], effTypeLabel[iEffType], iPTBin, iRapBin);
    	gEffMCTP_phiPol_pT_rap[iEffType][iEff][iFrame][iPTBin][iRapBin]->SetName(name);

    	sprintf(name, "%sEff_cosTheta_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
    	gEffMCTP_cosTheta_pT_rap[iEffType][iEff][iFrame][iPTBin][iRapBin] = (TEfficiency *) gDirectory->Get(name);
    	sprintf(name, "%sEff_cosTheta_%s_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], effTypeLabel[iEffType], iPTBin, iRapBin);
    	gEffMCTP_cosTheta_pT_rap[iEffType][iEff][iFrame][iPTBin][iRapBin]->SetName(name);
      }
    }
  }
}

//============================================
void PlotEffPol(Int_t iEff, Int_t iFrame){

  Char_t name[100];

  sprintf(name, "cPol_%s_phi_%s", effName[iEff], eff::frameLabel[iFrame]);
  TCanvas *c1 = new TCanvas(name, name);
  TH1F *hFrame1 = gPad->DrawFrame(-180, 0., 180., 1.1);
  hFrame1->SetXTitle(gEffMCTP_phiPol[TRUTH][iEff][iFrame]->GetPassedHistogram()->GetXaxis()->GetTitle());
  sprintf(name, "#varepsilon_{%s}", effName[iEff]);
  hFrame1->SetYTitle(name);
  gEffMCTP_phiPol[TRUTH][iEff][iFrame]->SetMarkerStyle(24);
  gEffMCTP_phiPol[TRUTH][iEff][iFrame]->Draw("p same");
  gEffMCTP_phiPol[TnP][iEff][iFrame]->SetMarkerStyle(20);
  gEffMCTP_phiPol[TnP][iEff][iFrame]->Draw("p same");

  TLegend *leg1 = new TLegend(0.7212644,0.1716102,0.9252874,0.3368644);
  leg1->AddEntry(gEffMCTP_phiPol[TRUTH][iEff][iFrame], "MC Truth", "p");
  leg1->AddEntry(gEffMCTP_phiPol[TnP][iEff][iFrame], "MC T&P", "p");
  leg1->SetFillColor(0); leg1->SetBorderSize(0);
  leg1->SetTextSize(0.04); leg1->Draw();
  
  TLine *line1 = new TLine(-180., 1., 180., 1.);
  line1->SetLineStyle(3); line1->Draw();

  sprintf(name, "Figures/%sEff_phi_%s.pdf", effName[iEff], eff::frameLabel[iFrame]);
  c1->Print(name);

  //================================
  sprintf(name, "cPol_%s_cosTheta_%s", effName[iEff], eff::frameLabel[iFrame]);
  TCanvas *c2 = new TCanvas(name, name);
  TH1F *hFrame2 = gPad->DrawFrame(-1, 0., 1., 1.1);
  hFrame2->SetXTitle(gEffMCTP_cosTheta[TRUTH][iEff][iFrame]->GetPassedHistogram()->GetXaxis()->GetTitle());
  sprintf(name, "#varepsilon_{%s}", effName[iEff]);
  hFrame2->SetYTitle(name);
  gEffMCTP_cosTheta[TRUTH][iEff][iFrame]->SetMarkerStyle(24);
  gEffMCTP_cosTheta[TRUTH][iEff][iFrame]->Draw("p same");
  gEffMCTP_cosTheta[TnP][iEff][iFrame]->SetMarkerStyle(20);
  gEffMCTP_cosTheta[TnP][iEff][iFrame]->Draw("p same");

  TLegend *leg2 = new TLegend(0.7212644,0.1716102,0.9252874,0.3368644);
  leg2->AddEntry(gEffMCTP_cosTheta[TRUTH][iEff][iFrame], "MC Truth", "p");
  leg2->AddEntry(gEffMCTP_cosTheta[TnP][iEff][iFrame], "MC T&P", "p");
  leg2->SetFillColor(0); leg2->SetBorderSize(0);
  leg2->SetTextSize(0.04); leg2->Draw();
  
  TLine *line2 = new TLine(-1., 1., 1., 1.);
  line2->SetLineStyle(3); line2->Draw();

  sprintf(name, "Figures/%sEff_cosTheta_%s.pdf", effName[iEff], eff::frameLabel[iFrame]);
  c2->Print(name);

  //=========================================
  sprintf(name, "c3TRUTH_cosTheta_phi_%sEff_%s", effName[iEff], eff::frameLabel[iFrame]);
  TCanvas *c3 = new TCanvas(name, name, 700, 700);
  gStyle->SetPaintTextFormat("5.2f");
  gEffMCTP2D_cosTheta_phiPol[TRUTH][iEff][iFrame]->Draw("text colz e");
  sprintf(name, "Figures/%sEff_MCTRUTH_cosThetaVsPhi_%s.pdf", effName[iEff], eff::frameLabel[iFrame]);
  c3->Print(name);

  //=========================================
  sprintf(name, "c4TnP_cosTheta_phi_%sEff_%s", effName[iEff], eff::frameLabel[iFrame]);
  TCanvas *c4 = new TCanvas(name, name, 700, 700);
  gStyle->SetPaintTextFormat("5.2f");
  gEffMCTP2D_cosTheta_phiPol[TnP][iEff][iFrame]->Draw("text colz e");
  sprintf(name, "Figures/%sEff_TnP_cosThetaVsPhi_%s.pdf", effName[iEff], eff::frameLabel[iFrame]);
  c4->Print(name);
  
}

//============================================
void PlotEffPol_pT_rap(Int_t iEff, Int_t iFrame, Int_t iRapBin, Int_t iPTBin){

  Char_t name[100];

  sprintf(name, "cPol_%s_phi_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
  TCanvas *c1 = new TCanvas(name, name);
  TH1F *hFrame1 = gPad->DrawFrame(-180, 0., 180., 1.1);
  hFrame1->SetXTitle(gEffMCTP_phiPol_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetXaxis()->GetTitle());
  sprintf(name, "#varepsilon_{%s}", effName[iEff]);
  hFrame1->SetYTitle(name);
  gEffMCTP_phiPol_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->SetMarkerStyle(24);
  gEffMCTP_phiPol_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->Draw("p same");
  gEffMCTP_phiPol_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin]->SetMarkerStyle(20);
  gEffMCTP_phiPol_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin]->Draw("p same");

  TLegend *leg1 = new TLegend(0.7212644,0.1716102,0.9252874,0.3368644);
  leg1->AddEntry(gEffMCTP_phiPol_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin], "MC Truth", "p");
  leg1->AddEntry(gEffMCTP_phiPol_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin], "MC T&P", "p");
  leg1->SetFillColor(0); leg1->SetBorderSize(0);
  leg1->SetTextSize(0.04); leg1->Draw();
  
  TLine *line1 = new TLine(-180., 1., 180., 1.);
  line1->SetLineStyle(3); line1->Draw();

  sprintf(name, "%1.1f < |y| < %1.1f, %1.0f < p_{T} < %1.0f GeV/c", 
	  eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], 
	  eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);
  TLatex *tex1 = new TLatex(-150., 0.1, name);
  tex1->SetTextSize(0.04); tex1->Draw();

  sprintf(name, "Figures/%sEff_phi_%s_rap%d_pT%d.pdf", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  c1->Print(name);

  //================================
  sprintf(name, "cPol_%s_cosTheta_%s_pT%d_rap%d", effName[iEff], eff::frameLabel[iFrame], iPTBin, iRapBin);
  TCanvas *c2 = new TCanvas(name, name);
  TH1F *hFrame2 = gPad->DrawFrame(-1, 0., 1., 1.1);
  hFrame2->SetXTitle(gEffMCTP_cosTheta_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetXaxis()->GetTitle());
  sprintf(name, "#varepsilon_{%s}", effName[iEff]);
  hFrame2->SetYTitle(name);
  gEffMCTP_cosTheta_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->SetMarkerStyle(24);
  gEffMCTP_cosTheta_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->Draw("p same");
  gEffMCTP_cosTheta_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin]->SetMarkerStyle(20);
  gEffMCTP_cosTheta_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin]->Draw("p same");

  TLegend *leg2 = new TLegend(0.7212644,0.1716102,0.9252874,0.3368644);
  leg2->AddEntry(gEffMCTP_cosTheta_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin], "MC Truth", "p");
  leg2->AddEntry(gEffMCTP_cosTheta_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin], "MC T&P", "p");
  leg2->SetFillColor(0); leg2->SetBorderSize(0);
  leg2->SetTextSize(0.04); leg2->Draw();
  
  TLine *line2 = new TLine(-1., 1., 1., 1.);
  line2->SetLineStyle(3); line2->Draw();

  sprintf(name, "%1.1f < |y| < %1.1f, %1.0f < p_{T} < %1.0f GeV/c", 
	  eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], 
	  eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);
  TLatex *tex2 = new TLatex(-0.9, 0.1, name);
  tex2->SetTextSize(0.04); tex2->Draw();

  sprintf(name, "Figures/%sEff_cosTheta_%s_rap%d_pT%d.pdf", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  c2->Print(name);
}


//============================================
void PlotRhoPol(Int_t iEff, Int_t iFrame){

  Char_t name[100];

  sprintf(name, "cRho_%s_phi_%s", effName[iEff], eff::frameLabel[iFrame]);
  TCanvas *c1 = new TCanvas(name, name);
  TH1F *hFrame1 = gPad->DrawFrame(-180, 0.8, 180., 1.2);
  hFrame1->SetXTitle(gEffMCTP_phiPol[TRUTH][iEff][iFrame]->GetPassedHistogram()->GetXaxis()->GetTitle());
  sprintf(name, "#rho (%s efficiency)", effName[iEff]);
  hFrame1->SetYTitle(name);
  // TH1D *hPassed1 = (TH1D*) gEffMCTP_phiPol[TRUTH][iEff][iFrame]->GetPassedHistogram();
  // TH1D *hTotal1 = (TH1D*) gEffMCTP_phiPol[TRUTH][iEff][iFrame]->GetTotalHistogram();
  // TH1D *hPassed2 =  (TH1D*) gEffMCTP_phiPol[TnP][iEff][iFrame]->GetPassedHistogram();
  // TH1D *hTotal2 =  (TH1D*) gEffMCTP_phiPol[TnP][iEff][iFrame]->GetTotalHistogram();
  // hPassed2->Divide(hPassed2,hTotal2,1.,1.,"b");
  // hPassed1->Divide(hPassed2);

  Double_t Np1, Np2, Nf1, Nf2;
  Double_t errNp1, errNp2, errNf1, errNf2;
  Double_t eff1, eff2, errEff1, errEff2, rho, errRho;
  sprintf(name, "rho_%s_phi_%s", effName[iEff], eff::frameLabel[iFrame]);
  TH1D *hRho = (TH1D *) ((TH1D *) gEffMCTP_phiPol[TRUTH][iEff][iFrame]->GetPassedHistogram())->Clone(name);
  hRho->Reset();
  for(int iBinsX = 1; iBinsX <= hRho->GetNbinsX(); iBinsX++){

    Np1 = gEffMCTP_phiPol[TRUTH][iEff][iFrame]->GetPassedHistogram()->GetBinContent(iBinsX);
    Np2 = gEffMCTP_phiPol[TnP][iEff][iFrame]->GetPassedHistogram()->GetBinContent(iBinsX);
    Nf1 = gEffMCTP_phiPol[TRUTH][iEff][iFrame]->GetTotalHistogram()->GetBinContent(iBinsX) - Np1;
    Nf2 = gEffMCTP_phiPol[TnP][iEff][iFrame]->GetTotalHistogram()->GetBinContent(iBinsX) - Np2;
    errNp1 = sqrt(Np1);    errNp2 = sqrt(Np2);
    errNf1 = sqrt(Nf1);    errNf2 = sqrt(Nf2);

    // printf("Np1 = %1.3f +- %1.3f, Np2 = %1.3f +- %1.3f, Nf1 = %1.3f +- %1.3f, Nf2 = %1.3f +- %1.3f\n",
    // 	   Np1, errNp1, Np2, errNp2, Nf1, errNf1, Nf2, errNf2);
    if((Np1+Nf1) > 0 && (Np2+Nf2) > 0 && Np1*errNp1 > 0 && Np2*errNp2 > 0){
      errEff1 = Np1/pow(Np1+Nf1,2) * sqrt(pow(Nf1/Np1*errNp1,2) + pow(errNf1,2));
      errEff2 = Np2/pow(Np2+Nf2,2) * sqrt(pow(Nf2/Np2*errNp2,2) + pow(errNf2,2));

      eff1 = Np1/(Np1+Nf1);
      eff2 = Np2/(Np2+Nf2);
      rho =  eff1 / eff2;

      errRho = rho * sqrt(pow(errEff1/eff1,2) + pow(errEff2/eff2,2));
      printf("eff1 = %1.3e +- %1.3e; eff2 = %1.3e +- %1.3e; rho = %1.3f +- %1.3f\n",
	     eff1, errEff1, eff2, errEff2, rho, errRho);
    }
    else{
      rho = 0.;
      errRho = 0.;
    }
    // hPassed1->SetBinError(iBinsX, errRho);
    hRho->SetBinContent(iBinsX, rho);
    hRho->SetBinError(iBinsX, errRho);
   }
  

  hRho->SetMarkerStyle(20);
  hRho->Draw("p same");

  TLine *line1 = new TLine(-180., 1., 180., 1.);
  line1->SetLineStyle(3); line1->Draw();

  sprintf(name, "Figures/rho_%sEff_phi_%s.pdf", effName[iEff], eff::frameLabel[iFrame]);
  c1->Print(name);

  //=========================================
  sprintf(name, "cRho_%s_cosTheta_%s", effName[iEff], eff::frameLabel[iFrame]);
  TCanvas *c2 = new TCanvas(name, name);
  TH1F *hFrame2 = gPad->DrawFrame(-1., 0.8, 1., 1.2);
  hFrame2->SetXTitle(gEffMCTP_cosTheta[TRUTH][iEff][iFrame]->GetPassedHistogram()->GetXaxis()->GetTitle());
  sprintf(name, "#rho (%s efficiency)", effName[iEff]);
  hFrame2->SetYTitle(name);
  // TH1D *hPassed1 = (TH1D*) gEffMCTP_cosTheta[TRUTH][iEff][iFrame]->GetPassedHistogram();
  // TH1D *hTotal1 = (TH1D*) gEffMCTP_cosTheta[TRUTH][iEff][iFrame]->GetTotalHistogram();
  // hPassed1->Divide(hPassed1,hTotal1,1.,1.,"b");
  // TH1D *hPassed2 =  (TH1D*) gEffMCTP_cosTheta[TnP][iEff][iFrame]->GetPassedHistogram();
  // TH1D *hTotal2 =  (TH1D*) gEffMCTP_cosTheta[TnP][iEff][iFrame]->GetTotalHistogram();
  // hPassed2->Divide(hPassed2,hTotal2,1.,1.,"b");
  // hPassed1->Divide(hPassed2);

  // hPassed1->SetMarkerStyle(20);
  // hPassed1->Draw("p same");

  TH1D *hRho2 = (TH1D *) ((TH1D *) gEffMCTP_cosTheta[TRUTH][iEff][iFrame]->GetPassedHistogram())->Clone(name);
  hRho2->Reset();
  for(int iBinsX = 1; iBinsX <= hRho2->GetNbinsX(); iBinsX++){

    Np1 = gEffMCTP_cosTheta[TRUTH][iEff][iFrame]->GetPassedHistogram()->GetBinContent(iBinsX);
    Np2 = gEffMCTP_cosTheta[TnP][iEff][iFrame]->GetPassedHistogram()->GetBinContent(iBinsX);
    Nf1 = gEffMCTP_cosTheta[TRUTH][iEff][iFrame]->GetTotalHistogram()->GetBinContent(iBinsX) - Np1;
    Nf2 = gEffMCTP_cosTheta[TnP][iEff][iFrame]->GetTotalHistogram()->GetBinContent(iBinsX) - Np2;

    errNp1 = sqrt(Np1);    errNp2 = sqrt(Np2);
    errNf1 = sqrt(Nf1);    errNf2 = sqrt(Nf2);

    // printf("Np1 = %1.3f +- %1.3f, Np2 = %1.3f +- %1.3f, Nf1 = %1.3f +- %1.3f, Nf2 = %1.3f +- %1.3f\n",
    // 	   Np1, errNp1, Np2, errNp2, Nf1, errNf1, Nf2, errNf2);
    if((Np1+Nf1) > 0 && (Np2+Nf2) > 0 && Np1*errNp1 > 0 && Np2*errNp2 > 0){
      errEff1 = Np1/pow(Np1+Nf1,2) * sqrt(pow(Nf1/Np1*errNp1,2) + pow(errNf1,2));
      errEff2 = Np2/pow(Np2+Nf2,2) * sqrt(pow(Nf2/Np2*errNp2,2) + pow(errNf2,2));

      eff1 = Np1/(Np1+Nf1);
      eff2 = Np2/(Np2+Nf2);
      rho =  eff1 / eff2;

      errRho = rho * sqrt(pow(errEff1/eff1,2) + pow(errEff2/eff2,2));
      printf("eff1 = %1.3e +- %1.3e; eff2 = %1.3e +- %1.3e; rho = %1.3f +- %1.3f\n",
	     eff1, errEff1, eff2, errEff2, rho, errRho);
    }
    else{
      rho = 0.;
      errRho = 0.;
    }

    // hPassed1->SetBinError(iBinsX, errRho);
    hRho2->SetBinContent(iBinsX, rho);
    hRho2->SetBinError(iBinsX, errRho);
   }

  hRho2->SetMarkerStyle(20);
  hRho2->Draw("p same");

  TLine *line2 = new TLine(-1., 1., 1., 1.);
  line2->SetLineStyle(3); line2->Draw();

  sprintf(name, "Figures/rho_%sEff_cosTheta_%s.pdf", effName[iEff], eff::frameLabel[iFrame]);
  c2->Print(name);

}


//============================================
void PlotRhoPol_pT_rap(Int_t iEff, Int_t iFrame, Int_t iRapBin, Int_t iPTBin){

  Char_t name[100];

  //=========================================
  sprintf(name, "cRho_%s_%s_rap%d_pt%d_cosTheta", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  TCanvas *c2 = new TCanvas(name, name);
  TH1F *hFrame2 = gPad->DrawFrame(-1., 0.8, 1., 1.2);
  hFrame2->SetXTitle(gEffMCTP_cosTheta_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetXaxis()->GetTitle());
  sprintf(name, "#rho (%s efficiency)", effName[iEff]);
  hFrame2->SetYTitle(name);
  // TH1D *hPassed1 = (TH1D*) gEffMCTP_cosTheta_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram();
  // TH1D *hTotal1 = (TH1D*) gEffMCTP_cosTheta_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetTotalHistogram();
  // hPassed1->Divide(hPassed1,hTotal1,1.,1.,"b");
  // TH1D *hPassed2 =  (TH1D*) gEffMCTP_cosTheta_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram();
  // TH1D *hTotal2 =  (TH1D*) gEffMCTP_cosTheta_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin]->GetTotalHistogram();
  // hPassed2->Divide(hPassed2,hTotal2,1.,1.,"b");
  // hPassed1->Divide(hPassed2);

  // hPassed1->SetMarkerStyle(20);
  // hPassed1->Draw("p same");

  Double_t Np1, Np2, Nf1, Nf2;
  Double_t errNp1, errNp2, errNf1, errNf2;
  Double_t eff1, eff2, errEff1, errEff2, rho, errRho;
  sprintf(name, "rho_%s_%s_rap%d_pt%d_cosTheta", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  TH1D *hRho2 = (TH1D *) ((TH1D *) gEffMCTP_cosTheta_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram())->Clone(name);
  hRho2->Reset();
  for(int iBinsX = 1; iBinsX <= hRho2->GetNbinsX(); iBinsX++){

    Np1 = gEffMCTP_cosTheta_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetBinContent(iBinsX);
    Np2 = gEffMCTP_cosTheta_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetBinContent(iBinsX);
    Nf1 = gEffMCTP_cosTheta_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetTotalHistogram()->GetBinContent(iBinsX) - Np1;
    Nf2 = gEffMCTP_cosTheta_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin]->GetTotalHistogram()->GetBinContent(iBinsX) - Np2;
    errNp1 = sqrt(Np1);    errNp2 = sqrt(Np2);
    errNf1 = sqrt(Nf1);    errNf2 = sqrt(Nf2);

    // printf("Np1 = %1.3f +- %1.3f, Np2 = %1.3f +- %1.3f, Nf1 = %1.3f +- %1.3f, Nf2 = %1.3f +- %1.3f\n",
    // 	   Np1, errNp1, Np2, errNp2, Nf1, errNf1, Nf2, errNf2);
    if((Np1+Nf1) > 0 && (Np2+Nf2) > 0){
      errEff1 = Np1/pow(Np1+Nf1,2) * sqrt(pow(Nf1/Np1*errNp1,2) + pow(errNf1,2));
      errEff2 = Np2/pow(Np2+Nf2,2) * sqrt(pow(Nf2/Np2*errNp2,2) + pow(errNf2,2));

      eff1 = Np1/(Np1+Nf1);
      eff2 = Np2/(Np2+Nf2);
      if(eff1 > 0 && eff2 > 0){
	rho =  eff1 / eff2;

	errRho = rho * sqrt(pow(errEff1/eff1,2) + pow(errEff2/eff2,2));
	printf("eff1 = %1.3e +- %1.3e; eff2 = %1.3e +- %1.3e; rho = %1.3f +- %1.3f\n",
	       eff1, errEff1, eff2, errEff2, rho, errRho);
      }
      else{
	rho = 0.;
	errRho = 0.;
      }
    }
    else{
      rho = 0.;
      errRho = 0.;
    }

    // hPassed1->SetBinError(iBinsX, errRho);
    hRho2->SetBinContent(iBinsX, rho);
    hRho2->SetBinError(iBinsX, errRho);
   }
  
  hRho2->SetMarkerStyle(20);
  hRho2->Draw("p same");

   TF1 *fit = new TF1("fit", "[0]  + (1+[1]*x*x)", -1, 1);
  //TF1 *fit = new TF1("fit", "[0] + (1./(3.+[1]) * (1+[1]*x*x))", -1, 1);
  fit->SetParameter(0,1.);
  fit->SetParameter(1,0.3);
  fit->SetLineWidth(1);
  hRho2->Fit("fit", "0+");
  fit = hRho2->GetFunction("fit");
  fit->Draw("same");
  Double_t lambdaTheta = fit->GetParameter(1);

  sprintf(name, "#lambda_{#theta}^{eff} = %1.2f", lambdaTheta);
  TLatex *tex2 = new TLatex(-0.9, 0.82, name);
  tex2->SetTextSize(0.05); tex2->Draw();
  sprintf(name, "%1.1f < |y| < %1.1f, %1.0f < p_{T} < %1.0f GeV/c", 
	  eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], 
	  eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);
  tex2->DrawLatex(-0.9, 1.17, name);

  TLine *line2 = new TLine(-1., 1., 1., 1.);
  line2->SetLineStyle(3); line2->Draw();

  sprintf(name, "Figures/rho_%sEff_cosTheta_%s_rap%d_pt%d.pdf", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  c2->Print(name);

  //======================================================

  sprintf(name, "cRho_%s_%s_rap%d_pt%d_phi", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  TCanvas *c1 = new TCanvas(name, name);
  TH1F *hFrame1 = gPad->DrawFrame(-180, 0.8, 180., 1.2);
  hFrame1->SetXTitle(gEffMCTP_phiPol_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetXaxis()->GetTitle());
  sprintf(name, "#rho (%s efficiency)", effName[iEff]);
  hFrame1->SetYTitle(name);
  // TH1D *hPassed1 = (TH1D*) gEffMCTP_phiPol[TRUTH][iEff][iFrame]->GetPassedHistogram();
  // TH1D *hTotal1 = (TH1D*) gEffMCTP_phiPol[TRUTH][iEff][iFrame]->GetTotalHistogram();
  // TH1D *hPassed2 =  (TH1D*) gEffMCTP_phiPol[TnP][iEff][iFrame]->GetPassedHistogram();
  // TH1D *hTotal2 =  (TH1D*) gEffMCTP_phiPol[TnP][iEff][iFrame]->GetTotalHistogram();
  // hPassed2->Divide(hPassed2,hTotal2,1.,1.,"b");
  // hPassed1->Divide(hPassed2);

  sprintf(name, "rho_%s_%s_rap%d_pt%d_phi", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  TH1D *hRho = (TH1D *) ((TH1D *) gEffMCTP_phiPol_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram())->Clone(name);
  hRho->Reset();
  for(int iBinsX = 1; iBinsX <= hRho->GetNbinsX(); iBinsX++){

    Np1 = gEffMCTP_phiPol_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetBinContent(iBinsX);
    Np2 = gEffMCTP_phiPol_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin]->GetPassedHistogram()->GetBinContent(iBinsX);
    Nf1 = gEffMCTP_phiPol_pT_rap[TRUTH][iEff][iFrame][iPTBin][iRapBin]->GetTotalHistogram()->GetBinContent(iBinsX) - Np1;
    Nf2 = gEffMCTP_phiPol_pT_rap[TnP][iEff][iFrame][iPTBin][iRapBin]->GetTotalHistogram()->GetBinContent(iBinsX) - Np2;
    errNp1 = sqrt(Np1);    errNp2 = sqrt(Np2);
    errNf1 = sqrt(Nf1);    errNf2 = sqrt(Nf2);

    if((Np1+Nf1) > 0 && (Np2+Nf2) > 0){
      // printf("Np1 = %1.3f +- %1.3f, Np2 = %1.3f +- %1.3f, Nf1 = %1.3f +- %1.3f, Nf2 = %1.3f +- %1.3f\n",
      // 	   Np1, errNp1, Np2, errNp2, Nf1, errNf1, Nf2, errNf2);
      errEff1 = Np1/pow(Np1+Nf1,2) * sqrt(pow(Nf1/Np1*errNp1,2) + pow(errNf1,2));
      errEff2 = Np2/pow(Np2+Nf2,2) * sqrt(pow(Nf2/Np2*errNp2,2) + pow(errNf2,2));

      eff1 = Np1/(Np1+Nf1);
      eff2 = Np2/(Np2+Nf2);
      if(eff1 > 0 && eff2 > 0){
	rho =  eff1 / eff2;

	errRho = rho * sqrt(pow(errEff1/eff1,2) + pow(errEff2/eff2,2));
	printf("eff1 = %1.3e +- %1.3e; eff2 = %1.3e +- %1.3e; rho = %1.3f +- %1.3f\n",
	       eff1, errEff1, eff2, errEff2, rho, errRho);
      }
      else{
	rho = 0.;
	errRho = 0.;
      }
    }
    else{
      rho = 0.;
      errRho = 0.;
    }
    // hPassed1->SetBinError(iBinsX, errRho);
    hRho->SetBinContent(iBinsX, rho);
    hRho->SetBinError(iBinsX, errRho);
   }
  

  hRho->SetMarkerStyle(20);
  hRho->Draw("p same");


  TF1 *fit2 = new TF1("fit2", "[0]  + (1.+((2.*[1])/(3.+[2])*cos(2.*x/180.*3.1415926535897931)))", -180., 180.);
  fit2->SetParameter(0,-1.);
  fit2->SetParameter(1,0.3); //lambda_phi
  fit2->FixParameter(2,lambdaTheta);
  fit2->SetLineWidth(1);
  hRho->Fit("fit2", "0+");
  fit2 = hRho->GetFunction("fit2");
  fit2->Draw("same");
  Double_t lambdaPhi = fit2->GetParameter(1);

  sprintf(name, "#lambda_{#phi}^{eff} = %1.2f", lambdaPhi);
  TLatex *tex1 = new TLatex(-170., 0.82, name);
  tex1->SetTextSize(0.05); tex1->Draw();
  sprintf(name, "%1.1f < |y| < %1.1f, %1.0f < p_{T} < %1.0f GeV/c", 
	  eff::rapForPTRange[iRapBin-1], eff::rapForPTRange[iRapBin], 
	  eff::pTRange[iRapBin][iPTBin-1], eff::pTRange[iRapBin][iPTBin]);
  tex1->DrawLatex(-170., 1.15, name);

  TLine *line1 = new TLine(-180., 1., 180., 1.);
  line1->SetLineStyle(3); line1->Draw();

  sprintf(name, "Figures/rho_%sEff_phi_%s_rap%d_pt%d.pdf", effName[iEff], eff::frameLabel[iFrame], iRapBin, iPTBin);
  c1->Print(name);

}
