#include "../interface/rootIncludes.inc"
#include "MCTruthEffAndAcc.C"

void BookHistos();
void WriteHistos();
//======================================
void runMCTruthEffAndAcc(Char_t *fileNameOut = "MCTruthEffAndAcc_HLTDimuon10JpsiBarrel_2Aug2012.root",
			 Int_t nSigma = 1,
			 Int_t resonance = UPS1S, //JPSI=0, PSIP=1, UPS1S=2, UPS2=3, UPS3S=4
			 Int_t trigTime = 1,//0...HLT_Dimuon5_Upsilon_Barrel_v3, 1...HLT_Dimuon5_Upsilon_Barrel_v5, 2...HLT_Dimuon7_Upsilon_Barrel_v1, 3...HLT_Dimuon9_Upsilon_Barrel_v1
			 Int_t useSoftMuons = 1, //1 = "soft"
			 Bool_t rejectCowboys = kTRUE
			 ){
// void runMCTruthEffAndAcc(Char_t *fileNameOut = "MCTruthEffAndAcc_HLTDimuon5UpsilonBarrel_19Dec2011.root",
// 		   Char_t *trigLabel = "HLT_Dimuon5_Upsilon_Barrel_v3",
// 		   Bool_t rejectCowboys = kTRUE
// 		   ){

  TChain *treeData = new TChain("data");
  if(resonance == JPSI){
    printf("\n\npreparing TTrees for J/psi processing\n");
    treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_JpsiLowPt_PGun_HLT1E33_3E33.root");
    treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_Jpsi_Pt9_21_PGun_HLT_1E33_3E33.root");
    treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_JpsiLowPt_PGun_HLT_1E33_3E33_2.root");
  }
  else if(resonance == UPS1S){
    printf("\n\npreparing TTrees for Ups(1S) processing...");
    if(useSoftMuons == 0){
      printf("original trees\n");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_1.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_2.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_3.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_4.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_5.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_6.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_7.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_8.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_9.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_10.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_11.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_12.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_13.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_14.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/December2011/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33_15.root");
    }
    else{
      printf("trees where TMOST flag needs to be applied by hand\n");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt30_rap1p3_HLT_1E33_3E33_20June2012_0.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt30_rap1p3_HLT_1E33_3E33_20June2012_1.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt30_rap1p3_HLT_1E33_3E33_20June2012_2.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt30_rap1p3_HLT_1E33_3E33_20June2012_3.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt30_rap1p3_HLT_1E33_3E33_20June2012_4.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt30_rap1p3_HLT_1E33_3E33_20June2012_5.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt30_rap1p3_HLT_1E33_3E33_20June2012_6.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt30_rap1p3_HLT_1E33_3E33_20June2012_7.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt30_rap1p3_HLT_1E33_3E33_20June2012_8.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt30_rap1p3_HLT_1E33_3E33_20June2012_9.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt15_rap1p3_HLT_1E33_3E33_20June2012_10.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_0pt15_rap1p3_HLT_1E33_3E33_20June2012_11.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_4pt15_rap1p3_HLT_1E33_3E33_20June2012_12.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_4pt15_rap1p3_HLT_1E33_3E33_20June2012_13.root");
      treeData->Add("/Users/hwoehri/CMS/Work/Data2011/FlatGen/June2012/TTree_Ups_PGun_4pt15_rap1p3_HLT_1E33_3E33_20June2012_14.root");
    }
  }
  else if(resonance == UPS2S){
    printf("preparing TTrees for Ups(2S) processing\n");
    //treeData->Add("/Users/hwoehri/CMS/Work/Data2011/Gen_DataPTDistribution/TTree_Onia2MuMu_Upsi2s_measuredPt.root"); //do NOT use (does not contain all GEN events)
    treeData->Add("/Users/hwoehri/CMS/Work/Data2011/Gen_DataPTDistribution/TTree_Onia2MuMu_Upsi2s_measuredPt_2April2012.root");
  }
  else if(resonance == UPS3S){
    printf("preparing TTrees for Ups(3S) processing\n");
    //treeData->Add("/Users/hwoehri/CMS/Work/Data2011/Gen_DataPTDistribution/TTree_Onia2MuMu_Upsi3s_measuredPt.root"); //do NOT use (does not contain all GEN events)
    treeData->Add("/Users/hwoehri/CMS/Work/Data2011/Gen_DataPTDistribution/TTree_Onia2MuMu_Upsi3s_measuredPt_2April2012.root");
  }
  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  printf("initializing tree\n");
  MCTruthEffAndAcc tree(treeData);
  printf("...done\n");  printf("booking the histos now\n");
  BookHistos();
  printf("...done\n");

  tree.Loop(resonance, trigTime, rejectCowboys, nSigma, useSoftMuons);
  WriteHistos();
  fOut->Close();
}

//==========================================
void BookHistos(){

  Char_t name[200], title[200];

  //======================================================
  //1D distributions: 
  //======================================================
  recoEff_pT = new TEfficiency("recoEff_pT", ";p_{T} [GeV/c]", eff::kNbPTBins[1], eff::pTRange[1]); 
  recoEff_y = new TEfficiency("recoEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); 
  recoEff_phi = new TEfficiency("recoEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); 
  //
  trigEff_pT = new TEfficiency("trigEff_pT", ";p_{T} [GeV/c]", eff::kNbPTBins[1], eff::pTRange[1]); 
  trigEff_y = new TEfficiency("trigEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); 
  trigEff_phi = new TEfficiency("trigEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); 
  //
  totEff_pT = new TEfficiency("totEff_pT", ";p_{T} [GeV/c]", eff::kNbPTBins[1], eff::pTRange[1]); 
  totEff_y = new TEfficiency("totEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); 
  totEff_phi = new TEfficiency("totEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); 

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){

    //1D distributions for phi in the polarization frames
    sprintf(name, "recoEff_phiPol_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";#phi_{%s}", eff::frameLabel[iFrame]);
    recoEff_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 
    sprintf(name, "trigEff_phiPol_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";#phi_{%s}", eff::frameLabel[iFrame]);
    trigEff_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 
    sprintf(name, "totEff_phiPol_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";#phi_{%s}", eff::frameLabel[iFrame]);
    totEff_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 

    //1D distributions for cosTheta in the polarization frames
    sprintf(name, "recoEff_cosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s}", eff::frameLabel[iFrame]);
    recoEff_cosTheta[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax); 
    sprintf(name, "trigEff_cosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s}", eff::frameLabel[iFrame]);
    trigEff_cosTheta[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax); 
    sprintf(name, "totEff_cosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s}", eff::frameLabel[iFrame]);
    totEff_cosTheta[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax); 

    // printf("cosTheta done, %p %p %p\n", recoEff_cosTheta[iFrame], trigEff_cosTheta[iFrame], totEff_cosTheta[iFrame]);

    //2D distributions for cosTheta vs phi in the polarization frames
    sprintf(name, "recoEff_phiVsCosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s};#phi_{%s}", eff::frameLabel[iFrame],eff::frameLabel[iFrame]);
    recoEff2D_cosTheta_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
    							eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 

    sprintf(name, "trigEff_phiVsCosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s};#phi_{%s}", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
    trigEff2D_cosTheta_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
    							eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 
    sprintf(name, "totEff_phiVsCosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s};#phi_{%s}", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
    totEff2D_cosTheta_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
    						       eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 
  }

  //======================================================
  //2D distributions: 
  //======================================================
  recoEff2D_pT_rapNP = new TEfficiency("recoEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTBins[1], eff::pTRange[1]);
  recoEff2D_pT_rap = new TEfficiency("recoEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTBins[1], eff::pTRange[1]);
  trigEff2D_pT_rapNP = new TEfficiency("trigEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTBins[1], eff::pTRange[1]);
  trigEff2D_pT_rap = new TEfficiency("trigEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTBins[1], eff::pTRange[1]);
  totEff2D_pT_rapNP = new TEfficiency("totEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTBins[1], eff::pTRange[1]);
  totEff2D_pT_rap = new TEfficiency("totEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTBins[1], eff::pTRange[1]);

  //======================================================
  //polarization histograms:
  //======================================================

  //histos for neg. and pos. rapidity separately:
  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*eff::kNbRapBins; iRapBin++){
      Int_t rapIndex = eff::kNbRapBins - iRapBin;
      if(iRapBin >= eff::kNbRapBins) rapIndex = iRapBin - eff::kNbRapBins + 1;
//       printf("iRapBin %d uses definitions of rapID %d\n", iRapBin, rapIndex);

      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[rapIndex]+1; iPTBin++){

	//all events
	sprintf(name, "totEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT2D, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol2D, eff::phiPolMin, eff::phiPolMax);
	//RECO events
	sprintf(name, "recoEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT2D, eff::cosTMin, eff::cosTMax,
								eff::kNbBinsPhiPol2D, eff::phiPolMin, eff::phiPolMax);

	//RECO + triggered events
	sprintf(name, "trigEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT2D, eff::cosTMin, eff::cosTMax,
								eff::kNbBinsPhiPol2D, eff::phiPolMin, eff::phiPolMax);
      }
    }
  }
  //histos taking together +y and -y:
  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	//(1) no folding in phi and cosTheta
	//all events
	sprintf(name, "totEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT2D, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol2D, eff::phiPolMin, eff::phiPolMax);
	//RECO events
	sprintf(name, "recoEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT2D, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol2D, eff::phiPolMin, eff::phiPolMax);
	//RECO + trigEff events
	sprintf(name, "trigEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT2D, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol2D, eff::phiPolMin, eff::phiPolMax);

	//(2) 4-folding in phi
	//all events
	sprintf(name, "totEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT/2, eff::cosTMin, eff::cosTMax,
								  eff::kNbBinsPhiPol, 0., 90.);
	//RECO events
	sprintf(name, "recoEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT/2, eff::cosTMin, eff::cosTMax,
								     eff::kNbBinsPhiPol, 0., 90.);

	//RECO + trigEff events
	sprintf(name, "trigEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT/2, eff::cosTMin, eff::cosTMax,
								     eff::kNbBinsPhiPol, 0., 90.);

	//add the 1D histos
	sprintf(name, "recoEff_phiPol_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";#phi_{%s}", eff::frameLabel[iFrame]);
	recoEff_phiPol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);

	sprintf(name, "recoEff_cosTheta_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s}", eff::frameLabel[iFrame]);
	recoEff_cosTheta_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax);

	sprintf(name, "trigEff_phiPol_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";#phi_{%s}", eff::frameLabel[iFrame]);
	trigEff_phiPol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);

	sprintf(name, "trigEff_cosTheta_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s}", eff::frameLabel[iFrame]);
	trigEff_cosTheta_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax);

	sprintf(name, "totEff_phiPol_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";#phi_{%s}", eff::frameLabel[iFrame]);
	totEff_phiPol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);

	sprintf(name, "totEff_cosTheta_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s}", eff::frameLabel[iFrame]);
	totEff_cosTheta_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax);

      }
    }
  }
}

//==========================================
void WriteHistos(){

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*eff::kNbRapBins; iRapBin++){
      Int_t rapIndex = eff::kNbRapBins - iRapBin;
      if(iRapBin >= eff::kNbRapBins) rapIndex = iRapBin - eff::kNbRapBins + 1;
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[rapIndex]+1; iPTBin++){
	totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Write();
	recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Write();
	trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
  //for(int iFrame = 0; iFrame < 2; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    //for(int iFrame = 0; iFrame < 2; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
	recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  //efficiencies
  recoEff_pT->Write();   recoEff_y->Write();   recoEff_phi->Write();
  trigEff_pT->Write();   trigEff_y->Write();   trigEff_phi->Write();
  totEff_pT->Write();   totEff_y->Write();   totEff_phi->Write();
  recoEff2D_pT_rapNP->Write();  recoEff2D_pT_rap->Write();
  trigEff2D_pT_rapNP->Write();  trigEff2D_pT_rap->Write();
  totEff2D_pT_rapNP->Write();  totEff2D_pT_rap->Write();

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    recoEff_phiPol[iFrame]->Write();
    trigEff_phiPol[iFrame]->Write();
    totEff_phiPol[iFrame]->Write();

    recoEff_cosTheta[iFrame]->Write();
    trigEff_cosTheta[iFrame]->Write();
    totEff_cosTheta[iFrame]->Write();
    
    recoEff2D_cosTheta_phiPol[iFrame]->Write();
    trigEff2D_cosTheta_phiPol[iFrame]->Write();
    totEff2D_cosTheta_phiPol[iFrame]->Write();
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	recoEff_phiPol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	trigEff_phiPol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	totEff_phiPol_pT_rap[iFrame][iPTBin][iRapBin]->Write();

	recoEff_cosTheta_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	trigEff_cosTheta_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	totEff_cosTheta_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }
}
