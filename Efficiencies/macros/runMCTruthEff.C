#include "../interface/rootIncludes.inc"
#include "MCTruthEff.C"

void BookHistos(Char_t *oniaLabel);
void WriteHistos();
//======================================
void runMCTruthEff(Char_t *fileNameOut = "MCTruthEff_HLTMu0TkMu0OSTJpsi_18May2011.root",
		   Char_t *trigLabel = "HLT_Mu0_TkMu0_OST_Jpsi", //"HLT_DoubleMu0", "HLT_Mu0_TkMu0_OST_Jpsi",
		   Bool_t rejectCowboys = kTRUE,
		   Char_t *fileNameIn = "/Users/hwoehri/CMS/CMSSW/hWoehri/Polarization/macros/JPsiToMuMu_Fall10-START38_V12-HLTrereco-WithAllMCEvents.root", 
		   //Char_t *fileNameIn = "/Users/hwoehri/CMS/Work/Polarization/Ilse/3April2011/JPsiToMuMu_pol_Fall10_noDimuVtxCut_1April2011.root",
		   //Char_t *fileNameIn = "/Users/hwoehri/CMS/CMSSW/hWoehri/Polarization/macros/JPsiToMuMu_Fall10-START38_V12-v1-Onia2MuMu-v6-WithAllMCEvents_merged.root",
		   Int_t selDimuType = 4, //0...only GG, 1... only GT, 2... only TT, 3...GG+GT, 4...GG+GT+TT		
		   Char_t *oniaLabel = "J/#psi" //"J/#psi", "#psi'", "Ups(1S)", "Ups(2S)", "Ups(3S)"
		   ){

  Char_t *nameDataSet = "data";
  TFile *fIn = new TFile(fileNameIn);
  TTree *treeData = (TTree*)fIn->Get(nameDataSet);

  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  printf("initializing tree\n");
  MCTruthEff tree(treeData);
  printf("...done\n");
  BookHistos(oniaLabel);

  tree.Loop(selDimuType, trigLabel, rejectCowboys);
  WriteHistos();
  fOut->Close();
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
  recoEff_pT = new TEfficiency("recoEff_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); 
  recoEff_y = new TEfficiency("recoEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); 
  recoEff_phi = new TEfficiency("recoEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); 
  //
  trigEff_pT = new TEfficiency("trigEff_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); 
  trigEff_y = new TEfficiency("trigEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); 
  trigEff_phi = new TEfficiency("trigEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); 
  //
  totEff_pT = new TEfficiency("totEff_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); 
  totEff_y = new TEfficiency("totEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); 
  totEff_phi = new TEfficiency("totEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); 

  //======================================================
  //2D distributions: 
  //======================================================
  recoEff2D_pT_rapNP = new TEfficiency("recoEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  recoEff2D_pT_rap = new TEfficiency("recoEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  trigEff2D_pT_rapNP = new TEfficiency("trigEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  trigEff2D_pT_rap = new TEfficiency("trigEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  totEff2D_pT_rapNP = new TEfficiency("totEff2D_pT_rapNP", ";y;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  totEff2D_pT_rap = new TEfficiency("totEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);

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
	totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	//RECO events
	sprintf(name, "recoEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);

	//RECO + triggered events
	sprintf(name, "trigEff2D_Onia_%s_pT%d_NPrap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
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
	totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	//RECO events
	sprintf(name, "recoEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	//RECO + trigEff events
	sprintf(name, "trigEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
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
      }
    }
  }

  //efficiencies vs. lab coordinates
  Int_t nBinsDeltaEta = 10, nBinsDeltaPhi = 9;
  Double_t deltaEtaMin = 0., deltaEtaMax = 2.5;
  Double_t deltaPhiMin = 0., deltaPhiMax = 180.;
  Int_t nBinsDeltaR = 32;
  Double_t deltaRMin = 0., deltaRMax = 3.2;
  Int_t nBinsDistM2 = 20;
  Double_t distMin = 0., distMax = 200.;
  for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

      //1.) as a function of deltaPhi vs deltaEta
      //all events
      sprintf(name, "totEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
      totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);

      //RECO events
      sprintf(name, "recoEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
      recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);

      //RECO + TRIG events
      sprintf(name, "trigEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
      trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);

      //2.) as a function of deltaR
      //all events
      sprintf(name, "totEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#DeltaR");
      totEff_deltaR_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDeltaR, deltaRMin, deltaRMax);

      //RECO events
      sprintf(name, "recoEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#DeltaR");
      recoEff_deltaR_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDeltaR, deltaRMin, deltaRMax);

      //RECO + TRIG events
      sprintf(name, "trigEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#DeltaR");
      trigEff_deltaR_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDeltaR, deltaRMin, deltaRMax);

      //3.) as a function of deltaRM2
      //all events
      sprintf(name, "totEff_deltaRM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#DeltaR(M2)");
      totEff_deltaRM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDeltaR, deltaRMin, deltaRMax);

      //RECO events
      sprintf(name, "recoEff_deltaRM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#DeltaR(M2)");
      recoEff_deltaRM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDeltaR, deltaRMin, deltaRMax);

      //RECO + TRIG events
      sprintf(name, "trigEff_deltaRM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#DeltaR(M2)");
      trigEff_deltaRM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDeltaR, deltaRMin, deltaRMax);

      //4.) as a function of deltaPhiM2
      //all events
      sprintf(name, "totEff_deltaPhiM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#phi(M2)");
      totEff_deltaPhiM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, 2*nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);

      //RECO events
      sprintf(name, "recoEff_deltaPhiM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#phi(M2)");
      recoEff_deltaPhiM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, 2*nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);

      //RECO + TRIG events
      sprintf(name, "trigEff_deltaPhiM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#phi(M2)");
      trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, 2*nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);

      //5.) as a function of deltaEtaM2
      //all events
      sprintf(name, "totEff_deltaEtaM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#eta(M2)");
      totEff_deltaEtaM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, 2*nBinsDeltaEta, deltaEtaMin, deltaEtaMax);

      //RECO events
      sprintf(name, "recoEff_deltaEtaM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#eta(M2)");
      recoEff_deltaEtaM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, 2*nBinsDeltaEta, deltaEtaMin, deltaEtaMax);

      //RECO + TRIG events
      sprintf(name, "trigEff_deltaEtaM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";#Delta#eta(M2)");
      trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, 2*nBinsDeltaEta, deltaEtaMin, deltaEtaMax);

      //6.) as a function of distM2
      //all events
      sprintf(name, "totEff_distM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";dist(M2) [cm]");
      totEff_distM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDistM2, distMin, distMax);

      //RECO events
      sprintf(name, "recoEff_distM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";dist(M2) [cm]");
      recoEff_distM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDistM2, distMin, distMax);

      //RECO + TRIG events
      sprintf(name, "trigEff_distM2_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";dist(M2) [cm]");
      trigEff_distM2_pT_rap[iPTBin][iRapBin] = new TEfficiency(name, title, nBinsDistM2, distMin, distMax);
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
  for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
      totEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Write();
      recoEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Write();
      trigEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Write();
    }
  }
  for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
      totEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Write();
      recoEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Write();
      trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Write();
    }
  }
  for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
      totEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Write();
      recoEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Write();
      trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Write();
    }
  }
  for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
      totEff_distM2_pT_rap[iPTBin][iRapBin]->Write();
      recoEff_distM2_pT_rap[iPTBin][iRapBin]->Write();
      trigEff_distM2_pT_rap[iPTBin][iRapBin]->Write();
    }
  }


  //efficiencies
  recoEff_pT->Write();   recoEff_y->Write();   recoEff_phi->Write();
  trigEff_pT->Write();   trigEff_y->Write();   trigEff_phi->Write();
  totEff_pT->Write();   totEff_y->Write();   totEff_phi->Write();
  recoEff2D_pT_rapNP->Write();  recoEff2D_pT_rap->Write();
  trigEff2D_pT_rapNP->Write();  trigEff2D_pT_rap->Write();
  totEff2D_pT_rapNP->Write();  totEff2D_pT_rap->Write();
}
