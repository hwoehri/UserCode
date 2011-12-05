#include "../interface/rootIncludes.inc"
#include "MCTnPEff.C"
#include "TEfficiency.h"

void LoadEfficiencies(Int_t iEff, Int_t iEffSample);
void BookHistos();
void WriteHistos();
//======================================
void runMCTnPEff(Char_t *fileNameOut = "MCTnPEff_HLTDimuon10JpsiBarrel_12Nov2011.root",
		 Int_t effSample = MC, //DATA, MC, MCTRUTH
		 Bool_t rejectCowboys = kTRUE,
		 Char_t *trigLabel = "HLT_Dimuon10_Jpsi_Barrel_v3", 
		 Char_t *fileNameIn = "/Users/hwoehri/CMS/Work/Data2011/FlatGen/PGun_HLT1E33_3E33_TTree.root"
		 ){

  TFile *fIn = new TFile(fileNameIn);
  TTree *treeData = (TTree*)fIn->Get("data");

  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  printf("initializing tree\n");
  MCTnPEff tree(treeData);  printf("...done\n");

  for(int iEff = 0; iEff < kNbEff; iEff++)
    LoadEfficiencies(iEff, effSample);
  printf("efficiencies loaded\n");

  BookHistos();

  tree.Loop(effSample, trigLabel, rejectCowboys);

  fOut->cd();
  WriteHistos();
  fOut->Close();
}

//==============================================================
void LoadEfficiencies(Int_t iEff, Int_t iEffSample){

  TFile *fIn = new TFile(effFileNames[iEff]);
  Char_t name[100];
  // if(iEff < TMgivenHLT){
    //store the central value of the efficiencies
    sprintf(name, "hEff_%s_central", effSampleName[iEffSample]);
    hMuEff[iEff][iEffSample][CENTRAL] = (TH2D *) gDirectory->Get(name);
    sprintf(name, "h%s_%s", effName[iEff], effSampleName[iEffSample]);
    hMuEff[iEff][iEffSample][CENTRAL]->SetName(name);
    printf("%s, histo %p\n", effName[iEff], hMuEff[iEff][iEffSample][CENTRAL]->GetName());
  // }
  // else{//eff for a TM given an HLT muon
  //   sprintf(name, "trigEff2D_pT_eta_BOTH");
  //   TEfficiency *eff = (TEfficiency *) gDirectory->Get(name);
  //   TH2D *hPassed = (TH2D *) eff->GetPassedHistogram();
  //   TH2D *hTot = (TH2D *) eff->GetTotalHistogram();
  //   sprintf(name, "h%s_%s", effName[iEff], effSampleName[iEffSample]);
  //   hMuEff[iEff][iEffSample][CENTRAL] = (TH2D *) hPassed->Clone(name);
  //   hMuEff[iEff][iEffSample][CENTRAL]->Divide(hTot);
  //   printf("%s, histo %p\n", effName[iEff], hMuEff[iEff][iEffSample][CENTRAL]->GetName());
  // }

  //load also the pT differential efficiencies from a fit
  // if(iEff == MuIDEff || iEff == MuQualEff || iEff == L1L2Eff || iEff == L3Eff){
  //   for(int iEta = 0; iEta < kNbEtaBins[iEff]; iEta++){
  //     sprintf(name, "fit%s_%s_pt_etaBin%d", effName[iEff], effSampleName[iEffSample], iEta);
  //     fMuEff_pT[iEff][iEffSample][iEta] = (TF1 *) gDirectory->Get(name);
  //   }
  // }
}

//==========================================
void BookHistos(){

  Char_t name[100], title[100];

  //======================================================
  //1D distributions: 
  //======================================================
  // hGen_pT = new TH1D("hGen_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); hGen_pT->Sumw2();
  // hGen_y = new TH1D("hGen_y", ";y", 2*eff::kNbRapBins, eff::rapRange); hGen_y->Sumw2();
  // hGen_phi = new TH1D("hGen_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); hGen_phi->Sumw2();
  // //
  // recoEff_pT = new TH1D("recoEff_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); recoEff_pT->Sumw2();
  // recoEff_y = new TH1D("recoEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); recoEff_y->Sumw2();
  // recoEff_phi = new TH1D("recoEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); recoEff_phi->Sumw2();
  // //
  // trigEff_pT = new TH1D("trigEff_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); trigEff_pT->Sumw2();
  // trigEff_y = new TH1D("trigEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); trigEff_y->Sumw2();
  // trigEff_phi = new TH1D("trigEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); trigEff_phi->Sumw2();
  // //
  // totEff_pT = new TH1D("totEff_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); totEff_pT->Sumw2();
  // totEff_y = new TH1D("totEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); totEff_y->Sumw2();
  // totEff_phi = new TH1D("totEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); totEff_phi->Sumw2();

  recoEff_pT  = new TEfficiency("recoEff_pT", ";p_{T} [GeV/c]", eff::kNbPTBins[1], eff::pTRange[1]); 
  recoEff_y   = new TEfficiency("recoEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); 
  recoEff_phi = new TEfficiency("recoEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); 

  trigEff_pT  = new TEfficiency("trigEff_pT", ";p_{T} [GeV/c]", eff::kNbPTBins[1], eff::pTRange[1]); 
  trigEff_y   = new TEfficiency("trigEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); 
  trigEff_phi = new TEfficiency("trigEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); 

  totEff_pT  = new TEfficiency("totEff_pT", ";p_{T} [GeV/c]", eff::kNbPTBins[1], eff::pTRange[1]); 
  totEff_y   = new TEfficiency("totEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); 
  totEff_phi = new TEfficiency("totEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); 

  //======================================================
  //2D distributions: 
  //======================================================
  // hGen2D_pT_rapNP = new TH2D("hGen2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  // hGen2D_pT_rap = new TH2D("hGen2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  // hGen2D_pT_rapNP->Sumw2(); hGen2D_pT_rap->Sumw2();
  // //
  // recoEff2D_pT_rapNP = new TH2D("recoEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  // recoEff2D_pT_rap = new TH2D("recoEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  // trigEff2D_pT_rapNP = new TH2D("trigEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  // trigEff2D_pT_rap = new TH2D("trigEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  // totEff2D_pT_rapNP = new TH2D("totEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  // totEff2D_pT_rap = new TH2D("totEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  // recoEff2D_pT_rapNP->Sumw2();  recoEff2D_pT_rap->Sumw2();  
  // trigEff2D_pT_rapNP->Sumw2();  trigEff2D_pT_rap->Sumw2();  
  // totEff2D_pT_rapNP->Sumw2(); totEff2D_pT_rap->Sumw2();

  recoEff2D_pT_rap = new TEfficiency("recoEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, 
  				     eff::kNbPTBins[1], eff::pTRange[1]);
  recoEff2D_pT_rapNP = new TEfficiency("recoEff2D_pT_rapNP", ";|y|;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, 
  				       eff::kNbPTBins[1], eff::pTRange[1]);

  trigEff2D_pT_rap = new TEfficiency("trigEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, 
  				     eff::kNbPTBins[1], eff::pTRange[1]);
  trigEff2D_pT_rapNP = new TEfficiency("trigEff2D_pT_rapNP", ";|y|;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, 
  				       eff::kNbPTBins[1], eff::pTRange[1]);

  totEff2D_pT_rap = new TEfficiency("totEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, 
  				     eff::kNbPTBins[1], eff::pTRange[1]);
  totEff2D_pT_rapNP = new TEfficiency("totEff2D_pT_rapNP", ";|y|;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, 
  				       eff::kNbPTBins[1], eff::pTRange[1]);
  //======================================================
  //polarization histograms:
  //======================================================

  //histos for neg. and pos. rapidity separately:
  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){

    sprintf(name, "recoEff_phiPol_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";#phi_{%s}", eff::frameLabel[iFrame]);
    recoEff_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 

    sprintf(name, "recoEff_cosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s}", eff::frameLabel[iFrame]);
    recoEff_cosTheta[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax); 

    sprintf(name, "recoEff_phiVsCosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s};#phi_{%s}", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
    recoEff2D_cosTheta_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
    							eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 

    sprintf(name, "trigEff_phiPol_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";#phi_{%s}", eff::frameLabel[iFrame]);
    trigEff_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 

    sprintf(name, "trigEff_cosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s}", eff::frameLabel[iFrame]);
    trigEff_cosTheta[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax); 

    sprintf(name, "trigEff_phiVsCosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s};#phi_{%s}", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
    trigEff2D_cosTheta_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
    							eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 

    sprintf(name, "totEff_phiPol_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";#phi_{%s}", eff::frameLabel[iFrame]);
    totEff_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 

    sprintf(name, "totEff_cosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s}", eff::frameLabel[iFrame]);
    totEff_cosTheta[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax); 

    sprintf(name, "totEff_phiVsCosTheta_%s", eff::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s};#phi_{%s}", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
    totEff2D_cosTheta_phiPol[iFrame] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
    							eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax); 
  }

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*eff::kNbRapBins; iRapBin++){
      Int_t rapIndex = eff::kNbRapBins - iRapBin;
      if(iRapBin >= eff::kNbRapBins) rapIndex = iRapBin - eff::kNbRapBins + 1;
//       printf("iRapBin %d uses definitions of rapID %d\n", iRapBin, rapIndex);

      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[rapIndex]+1; iPTBin++){

	sprintf(name, "recoEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);

	sprintf(name, "trigEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);

	sprintf(name, "totEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);

	// //all events
	// sprintf(name, "hGen2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax, 
	// 							 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	// hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();

	// //total efficiency
	// sprintf(name, "totEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
	// 							 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	// totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();

	// //RECO events
	// sprintf(name, "recoEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
	// 							eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	// recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();

	// //RECO + triggered events
	// sprintf(name, "trigEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
	// 							eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	// trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();
      }
    }
  }
  //histos taking together +y and -y:
  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

	//(1) no folding in phi and cosTheta
	sprintf(name, "recoEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);

	sprintf(name, "trigEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);

	sprintf(name, "totEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	//(2) 4-folding in phi
	sprintf(name, "recoEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT/2, eff::cosTMin, eff::cosTMax,
								     eff::kNbBinsPhiPol, 0., 90.);

	sprintf(name, "trigEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT/2, eff::cosTMin, eff::cosTMax,
								     eff::kNbBinsPhiPol, 0., 90.);

	sprintf(name, "totEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT/2, eff::cosTMin, eff::cosTMax,
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

	// //(1) no folding in phi and cosTheta
	// //all events
	// sprintf(name, "hGen2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
	// 							 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	// hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();

	// //total efficiency
	// sprintf(name, "totEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
	// 							 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	// totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
	// //RECO events
	// sprintf(name, "recoEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
	// 							 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	// recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
	// //RECO + trigEff events
	// sprintf(name, "trigEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
	// 							 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	// trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
	// //(2) 4-folding in phi
	// //all events
	// sprintf(name, "hGen2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
	// 							  eff::kNbBinsPhiPol, 0., 90.);
	// hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();
	// //total efficiencies
	// sprintf(name, "totEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
	// 							  eff::kNbBinsPhiPol, 0., 90.);
	// totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();
	// //RECO events
	// sprintf(name, "recoEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
	// 							     eff::kNbBinsPhiPol, 0., 90.);
	// recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();
	// //RECO + trigEff events
	// sprintf(name, "trigEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	// sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	// trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
	// 							     eff::kNbBinsPhiPol, 0., 90.);
 	// trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();
      }
    }
  }


  // //efficiencies vs. lab coordinates
  // Int_t nBinsDeltaEta = 10, nBinsDeltaPhi = 9;
  // Double_t deltaEtaMin = 0., deltaEtaMax = 2.5;
  // Double_t deltaPhiMin = 0., deltaPhiMax = 180.;
  // Int_t nBinsDeltaR = 32;
  // Double_t deltaRMin = 0., deltaRMax = 3.2;
  // Int_t nBinsDistM2 = 20;
  // Double_t distMin = 0., distMax = 200.;
  // for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

  //     //1.) as a function of deltaPhi vs deltaEta
  //     //all events
  //     sprintf(name, "hGen2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
  //     hGen2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
  //     hGen2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //total efficiencies
  //     sprintf(name, "totEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
  //     totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
  //     totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO events
  //     sprintf(name, "recoEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
  //     recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
  //     recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO + TRIG events
  //     sprintf(name, "trigEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
  //     trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
  //     trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //2.) as a function of deltaR
  //     //all events
  //     sprintf(name, "hGen_deltaR_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#DeltaR");
  //     hGen_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
  //     hGen_deltaR_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //total efficiencies
  //     sprintf(name, "totEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#DeltaR");
  //     totEff_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
  //     totEff_deltaR_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO events
  //     sprintf(name, "recoEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#DeltaR");
  //     recoEff_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
  //     recoEff_deltaR_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO + TRIG events
  //     sprintf(name, "trigEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#DeltaR");
  //     trigEff_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
  //     trigEff_deltaR_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //3.) as a function of deltaRM2
  //     //all events
  //     sprintf(name, "hGen_deltaRM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#DeltaR(M2)");
  //     hGen_deltaRM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
  //     hGen_deltaRM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //total efficiencies
  //     sprintf(name, "totEff_deltaRM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#DeltaR(M2)");
  //     totEff_deltaRM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
  //     totEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO events
  //     sprintf(name, "recoEff_deltaRM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#DeltaR(M2)");
  //     recoEff_deltaRM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
  //     recoEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO + TRIG events
  //     sprintf(name, "trigEff_deltaRM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#DeltaR(M2)");
  //     trigEff_deltaRM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
  //     trigEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //4.) as a function of deltaPhiM2
  //     //all events
  //     sprintf(name, "hGen_deltaPhiM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#phi(M2)");
  //     hGen_deltaPhiM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
  //     hGen_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //total efficiencies
  //     sprintf(name, "totEff_deltaPhiM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#phi(M2)");
  //     totEff_deltaPhiM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
  //     totEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO events
  //     sprintf(name, "recoEff_deltaPhiM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#phi(M2)");
  //     recoEff_deltaPhiM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
  //     recoEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO + TRIG events
  //     sprintf(name, "trigEff_deltaPhiM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#phi(M2)");
  //     trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
  //     trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //5.) as a function of deltaEtaM2
  //     //all events
  //     sprintf(name, "hGen_deltaEtaM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#eta(M2)");
  //     hGen_deltaEtaM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaEta, deltaEtaMin, deltaEtaMax);
  //     hGen_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //total efficiencies
  //     sprintf(name, "totEff_deltaEtaM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#eta(M2)");
  //     totEff_deltaEtaM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaEta, deltaEtaMin, deltaEtaMax);
  //     totEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO events
  //     sprintf(name, "recoEff_deltaEtaM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#eta(M2)");
  //     recoEff_deltaEtaM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaEta, deltaEtaMin, deltaEtaMax);
  //     recoEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO + TRIG events
  //     sprintf(name, "trigEff_deltaEtaM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";#Delta#phi(M2)");
  //     trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaEta, deltaEtaMin, deltaEtaMax);
  //     trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //6.) as a function of distM2
  //     //all events
  //     sprintf(name, "hGen_distM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";dist(M2) [cm]");
  //     hGen_distM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDistM2, distMin, distMax);
  //     hGen_distM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //total efficiencies
  //     sprintf(name, "totEff_distM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";dist(M2) [cm]");
  //     totEff_distM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDistM2, distMin, distMax);
  //     totEff_distM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO events
  //     sprintf(name, "recoEff_distM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";dist(M2) [cm]");
  //     recoEff_distM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDistM2, distMin, distMax);
  //     recoEff_distM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //     //RECO + TRIG events
  //     sprintf(name, "trigEff_distM2_pT%d_rap%d", iPTBin, iRapBin);
  //     sprintf(title, ";dist(M2) [cm]");
  //     trigEff_distM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDistM2, distMin, distMax);
  //     trigEff_distM2_pT_rap[iPTBin][iRapBin]->Sumw2();

  //   }
  // }

  hCorrRECO = new TH2F("hCorrRECO", "", 2, 0., 2., 2, 0., 2.);
}

//==========================================
void WriteHistos(){

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
  //for(int iFrame = 0; iFrame < 3; iFrame++){
    recoEff_phiPol[iFrame]->Write();
    recoEff_cosTheta[iFrame]->Write();
    recoEff2D_cosTheta_phiPol[iFrame]->Write();
    trigEff_phiPol[iFrame]->Write();
    trigEff_cosTheta[iFrame]->Write();
    trigEff2D_cosTheta_phiPol[iFrame]->Write();
    totEff_phiPol[iFrame]->Write();
    totEff_cosTheta[iFrame]->Write();
    totEff2D_cosTheta_phiPol[iFrame]->Write();
  }

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    //for(int iFrame = 0; iFrame < 3; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*eff::kNbRapBins; iRapBin++){
      Int_t rapIndex = eff::kNbRapBins - iRapBin;
      if(iRapBin >= eff::kNbRapBins) rapIndex = iRapBin - eff::kNbRapBins + 1;
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[rapIndex]+1; iPTBin++){
	recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Write();
	trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Write();
	totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    //for(int iFrame = 0; iFrame < 3; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    //for(int iFrame = 0; iFrame < 3; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
	totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    //for(int iFrame = 0; iFrame < 3; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	recoEff_phiPol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	recoEff_cosTheta_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	trigEff_phiPol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	trigEff_cosTheta_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	totEff_phiPol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	totEff_cosTheta_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    //for(int iFrame = 0; iFrame < 3; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
	totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  // for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
  //     recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Write();
  //     trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Write();
  //     totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Write();
  //   }
  // }
  // for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
  //     totEff_deltaR_pT_rap[iPTBin][iRapBin]->Write();
  //     recoEff_deltaR_pT_rap[iPTBin][iRapBin]->Write();
  //     trigEff_deltaR_pT_rap[iPTBin][iRapBin]->Write();
  //   }
  // }
  // for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
  //     totEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Write();
  //     recoEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Write();
  //     trigEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Write();
  //   }
  // }
  // for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
  //     totEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Write();
  //     recoEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Write();
  //     trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Write();
  //   }
  // }

  // for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
  //     totEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Write();
  //     recoEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Write();
  //     trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Write();
  //   }
  // }

  // for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
  //     totEff_distM2_pT_rap[iPTBin][iRapBin]->Write();
  //     recoEff_distM2_pT_rap[iPTBin][iRapBin]->Write();
  //     trigEff_distM2_pT_rap[iPTBin][iRapBin]->Write();
  //   }
  // }

  //efficiencies
  recoEff_pT->Write();   recoEff_y->Write();   recoEff_phi->Write();
  trigEff_pT->Write();   trigEff_y->Write();   trigEff_phi->Write();
  totEff_pT->Write();   totEff_y->Write();   totEff_phi->Write();
  recoEff2D_pT_rapNP->Write();  recoEff2D_pT_rap->Write();
  trigEff2D_pT_rapNP->Write();  trigEff2D_pT_rap->Write();
  totEff2D_pT_rapNP->Write();  totEff2D_pT_rap->Write();
  
  hCorrRECO->Write();
}
