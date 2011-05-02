#include "../interface/rootIncludes.inc"
#include "MCTnPEff.C"

void LoadEfficiencies(Int_t iEff, Int_t iEffSample);
void BookHistos(Char_t *oniaLabel);
void WriteHistos();
void DivideHistos();
//======================================
void runMCTnPEff(Char_t *fileNameOut = "MCTnPEff_20March2011.root",
		 Int_t effSample = MC, //DATA, MC, MCTRUTH
		 Char_t *trigLabel = "HLT_DoubleMu0", //"HLT_DoubleMu0", "HLT_Mu0_TkMu0_OST_Jpsi",
		 Char_t *fileNameIn = "/Users/hwoehri/CMS/CMSSW/hWoehri/Polarization/macros/JPsiToMuMu_Fall10-START38_V12-HLTrereco-WithAllMCEvents.root",
		 //Char_t *fileNameIn = "/Users/hwoehri/CMS/Work/Polarization/Ilse/3April2011/JPsiToMuMu_pol_Fall10_noDimuVtxCut_1April2011.root",
		 // 		 Char_t *fileNameIn = "/home/hermine/CMS/CMSSW/hWoehri/Polarization/macros/JPsiToMuMu_Fall10-START38_V12-v1-Onia2MuMu-v6-WithAllMCEvents_merged.root",
		 Char_t *oniaLabel = "J/#psi" //"J/#psi", "#psi'", "Ups(1S)", "Ups(2S)", "Ups(3S)"
		 ){

  TFile *fIn = new TFile(fileNameIn);
  TTree *treeData = (TTree*)fIn->Get("data");

  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  printf("initializing tree\n");
  MCTnPEff tree(treeData);  printf("...done\n");

  for(int iEff = 0; iEff < kNbEff; iEff++)
    LoadEfficiencies(iEff, effSample);
  printf("efficiencies loaded\n");

  BookHistos(oniaLabel);

  tree.Loop(effSample, trigLabel);

  DivideHistos();
  fOut->cd();
  WriteHistos();
  fOut->Close();
}

//==============================================================
void LoadEfficiencies(Int_t iEff, Int_t iEffSample){

  TFile *fIn = new TFile(effFileNames[iEff]);
  Char_t name[100];
  //store the central value of the efficiencies
  sprintf(name, "hEff_%s_central", effSampleName[iEffSample]);
  hMuEff[iEff][iEffSample][CENTRAL] = (TH2D *) gDirectory->Get(name);
  sprintf(name, "h%s_%s", effName[iEff], effSampleName[iEffSample]);
  hMuEff[iEff][iEffSample][CENTRAL]->SetName(name);
  printf("%s, histo %p\n", effName[iEff], hMuEff[iEff][iEffSample][CENTRAL]->GetName());

  //load also the pT differential efficiencies from a fit
  if(iEff == MuIDEff || iEff == MuQualEff || iEff == L1L2Eff || iEff == L3Eff){
    for(int iEta = 0; iEta < kNbEtaBins; iEta++){
      sprintf(name, "fit%s_%s_pt_etaBin%d", effName[iEff], effSampleName[iEffSample], iEta);
      fMuEff_pT[iEff][iEffSample][iEta] = (TF1 *) gDirectory->Get(name);
    }
  }
}

//==========================================
void BookHistos(Char_t *oniaLabel){

//   //initialise the "massMuOnia" variable
//   std::map<std::string, double> fOniaToMass;
//   fOniaToMass["J/#psi"] = 3.0969;
//   fOniaToMass["#psi'"] = 3.686;
//   fOniaToMass["Ups(1S)"] = 9.460;
//   fOniaToMass["Ups(2S)"] = 10.023;
//   fOniaToMass["Ups(3S)"] = 10.355;

//   std::map<std::string, double>::iterator iter = fOniaToMass.find(oniaLabel);
//   if(iter != fOniaToMass.end())
//     massMuOnia = iter->second;
//   printf("will use a mass of %1.3f GeV for the %s\n", massMuOnia, oniaLabel);
//   //===================================

  Char_t name[100], title[100];

  //======================================================
  //1D distributions: 
  //======================================================
  hGen_pT = new TH1D("hGen_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); hGen_pT->Sumw2();
  hGen_y = new TH1D("hGen_y", ";y", 2*eff::kNbRapBins, eff::rapRange); hGen_y->Sumw2();
  hGen_phi = new TH1D("hGen_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); hGen_phi->Sumw2();
  //
  recoEff_pT = new TH1D("recoEff_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); recoEff_pT->Sumw2();
  recoEff_y = new TH1D("recoEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); recoEff_y->Sumw2();
  recoEff_phi = new TH1D("recoEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); recoEff_phi->Sumw2();
  //
  trigEff_pT = new TH1D("trigEff_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); trigEff_pT->Sumw2();
  trigEff_y = new TH1D("trigEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); trigEff_y->Sumw2();
  trigEff_phi = new TH1D("trigEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); trigEff_phi->Sumw2();
  //
  totEff_pT = new TH1D("totEff_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); totEff_pT->Sumw2();
  totEff_y = new TH1D("totEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); totEff_y->Sumw2();
  totEff_phi = new TH1D("totEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); totEff_phi->Sumw2();

  //======================================================
  //2D distributions: 
  //======================================================
  hGen2D_pT_rapNP = new TH2D("hGen2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  hGen2D_pT_rap = new TH2D("hGen2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  hGen2D_pT_rapNP->Sumw2(); hGen2D_pT_rap->Sumw2();
  //
  recoEff2D_pT_rapNP = new TH2D("recoEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  recoEff2D_pT_rap = new TH2D("recoEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  trigEff2D_pT_rapNP = new TH2D("trigEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  trigEff2D_pT_rap = new TH2D("trigEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  totEff2D_pT_rapNP = new TH2D("totEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  totEff2D_pT_rap = new TH2D("totEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  recoEff2D_pT_rapNP->Sumw2();  recoEff2D_pT_rap->Sumw2();  trigEff2D_pT_rapNP->Sumw2();
  trigEff2D_pT_rap->Sumw2();  totEff2D_pT_rap->Sumw2();

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
	sprintf(name, "hGen2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax, 
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();

	//total efficiency
	sprintf(name, "totEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();

	//RECO events
	sprintf(name, "recoEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();

	//RECO + triggered events
	sprintf(name, "trigEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();
      }
    }
  }
  //histos taking together +y and -y:
  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	//(1) no folding in phi and cosTheta
	//all events
	sprintf(name, "hGen2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();

	//total efficiency
	sprintf(name, "totEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
	//RECO events
	sprintf(name, "recoEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
	//RECO + trigEff events
	sprintf(name, "trigEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
	//(2) 4-folding in phi
	//all events
	sprintf(name, "hGen2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								  eff::kNbBinsPhiPol, 0., 90.);
	hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();
	//total efficiencies
	sprintf(name, "totEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								  eff::kNbBinsPhiPol, 0., 90.);
	totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();
	//RECO events
	sprintf(name, "recoEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								     eff::kNbBinsPhiPol, 0., 90.);
	recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();
	//RECO + trigEff events
	sprintf(name, "trigEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								     eff::kNbBinsPhiPol, 0., 90.);
 	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();
     }
    }
  }


  //efficiencies vs. lab coordinates
  Int_t nBinsDeltaEta = 10, nBinsDeltaPhi = 9;
  Double_t deltaEtaMin = 0., deltaEtaMax = 2.5;
  Double_t deltaPhiMin = 0., deltaPhiMax = 180.;
  Int_t nBinsDeltaR = 32;
  Double_t deltaRMin = 0., deltaRMax = 3.2;
  for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

      //1.) as a function of deltaPhi vs deltaEta
      //all events
      sprintf(name, "hGen2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
      hGen2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);

      //total efficiencies
      sprintf(name, "totEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
      totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);

      //RECO events
      sprintf(name, "recoEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
      recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);

      //RECO + TRIG events
      sprintf(name, "trigEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
      trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);

      //2.) as a function of deltaR
      //all events
      sprintf(name, "hGen_deltaR_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#DeltaR");
      hGen_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);

      //total efficiencies
      sprintf(name, "totEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#DeltaR");
      totEff_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);

      //RECO events
      sprintf(name, "recoEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#DeltaR");
      recoEff_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);

      //RECO + TRIG events
      sprintf(name, "trigEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#DeltaR");
      trigEff_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
    }
  }
}

//==========================================
void DivideHistos(){

  recoEff_pT->Divide(recoEff_pT, hGen_pT, 1., 1., "B"); //binomial errors
  recoEff_y->Divide(recoEff_y, hGen_y, 1., 1., "B");
  recoEff_phi->Divide(recoEff_phi, hGen_phi, 1., 1., "B");

  trigEff_pT->Divide(trigEff_pT, hGen_pT, 1., 1., "B");
  trigEff_y->Divide(trigEff_y, hGen_y, 1., 1., "B");
  trigEff_phi->Divide(trigEff_phi, hGen_phi, 1., 1., "B");

  totEff_pT->Divide(totEff_pT, hGen_pT, 1., 1., "B");
  totEff_y->Divide(totEff_y, hGen_y, 1., 1., "B");
  totEff_phi->Divide(totEff_phi, hGen_phi, 1., 1., "B");

  recoEff2D_pT_rapNP->Divide(recoEff2D_pT_rapNP, hGen2D_pT_rapNP, 1., 1., "B");
  recoEff2D_pT_rap->Divide(recoEff2D_pT_rap, hGen2D_pT_rap);
  trigEff2D_pT_rapNP->Divide(trigEff2D_pT_rapNP, hGen2D_pT_rapNP, 1., 1., "B");
  trigEff2D_pT_rap->Divide(trigEff2D_pT_rap, hGen2D_pT_rap);
  totEff2D_pT_rapNP->Divide(totEff2D_pT_rapNP, hGen2D_pT_rapNP, 1., 1., "B");
  totEff2D_pT_rap->Divide(totEff2D_pT_rap, hGen2D_pT_rap);

  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*eff::kNbRapBins; iRapBin++){
      Int_t rapIndex = eff::kNbRapBins - iRapBin;
      if(iRapBin >= eff::kNbRapBins) rapIndex = iRapBin - eff::kNbRapBins + 1;
//       printf("iRapBin %d uses definitions of rapID %d\n", iRapBin, rapIndex);

      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[rapIndex]+1; iPTBin++){
	recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Divide(recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], 1., 1., "B");
	trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Divide(trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], 1., 1., "B");
	totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Divide(totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], 1., 1., "B");
      }
    }

    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Divide(recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin], 1., 1., "B");
	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Divide(trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin], 1., 1., "B");
	totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Divide(totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin], 1., 1., "B");

	recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Divide(recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], 1., 1., "B");
	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Divide(trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], 1., 1., "B");
	totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Divide(totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], 1., 1., "B");
      }
    }
  }

  for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

      recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Divide(recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], hGen2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], 1., 1., "B");
      trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Divide(trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], hGen2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], 1., 1., "B");
      totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Divide(totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], hGen2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], 1., 1., "B");

      recoEff_deltaR_pT_rap[iPTBin][iRapBin]->Divide(recoEff_deltaR_pT_rap[iPTBin][iRapBin], hGen_deltaR_pT_rap[iPTBin][iRapBin], 1., 1., "B");
      trigEff_deltaR_pT_rap[iPTBin][iRapBin]->Divide(trigEff_deltaR_pT_rap[iPTBin][iRapBin], hGen_deltaR_pT_rap[iPTBin][iRapBin], 1., 1., "B");
      totEff_deltaR_pT_rap[iPTBin][iRapBin]->Divide(totEff_deltaR_pT_rap[iPTBin][iRapBin], hGen_deltaR_pT_rap[iPTBin][iRapBin], 1., 1., "B");
    }
  }
}

//==========================================
void WriteHistos(){

//   for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
  for(int iFrame = 0; iFrame < 2; iFrame++){
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

//   for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
  for(int iFrame = 0; iFrame < 2; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

//   for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
  for(int iFrame = 0; iFrame < 2; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
	recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }
  for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
      totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Write();
      recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Write();
      trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Write();
    }
  }
  for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
      totEff_deltaR_pT_rap[iPTBin][iRapBin]->Write();
      recoEff_deltaR_pT_rap[iPTBin][iRapBin]->Write();
      trigEff_deltaR_pT_rap[iPTBin][iRapBin]->Write();
    }
  }

  //efficiencies (in histo form)
  recoEff_pT->Write();   recoEff_y->Write();   recoEff_phi->Write();
  trigEff_pT->Write();   trigEff_y->Write();   trigEff_phi->Write();
  totEff_pT->Write();   totEff_y->Write();   totEff_phi->Write();
  recoEff2D_pT_rapNP->Write();  recoEff2D_pT_rap->Write();
  trigEff2D_pT_rapNP->Write();  trigEff2D_pT_rap->Write();
  totEff2D_pT_rapNP->Write();  totEff2D_pT_rap->Write();

}
