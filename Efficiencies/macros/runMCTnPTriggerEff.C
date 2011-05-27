#include "../interface/rootIncludes.inc"
#include "MCTnPTriggerEff.C"

void LoadEfficiencies(Int_t iEff, Int_t iEffSample);
void BookHistos(Char_t *oniaLabel);
void WriteHistos();
void DivideHistos();
//======================================
void runMCTnPTriggerEff(Char_t *fileNameOut = "MCTnPTriggerEff_HLTMu0TkMu0OSTJpsi_23May2011.root",
			Int_t effSample = MC, //DATA, MC, MCTRUTH
			Char_t *trigLabel = "HLT_Mu0_TkMu0_OST_Jpsi", //"HLT_DoubleMu0", "HLT_Mu0_TkMu0_OST_Jpsi",
			Char_t *fileNameIn = "/Users/hwoehri/CMS/CMSSW/hWoehri/Polarization/macros/JPsiToMuMu_Fall10-START38_V12-HLTrereco-WithAllMCEvents.root",
			//Char_t *fileNameIn = "/Users/hwoehri/CMS/Work/Polarization/Ilse/3April2011/JPsiToMuMu_pol_Fall10_noDimuVtxCut_1April2011.root",
			// 		 Char_t *fileNameIn = "/home/hermine/CMS/CMSSW/hWoehri/Polarization/macros/JPsiToMuMu_Fall10-START38_V12-v1-Onia2MuMu-v6-WithAllMCEvents_merged.root",
			Char_t *oniaLabel = "J/#psi" //"J/#psi", "#psi'", "Ups(1S)", "Ups(2S)", "Ups(3S)"
		 ){

  TFile *fIn = new TFile(fileNameIn);
  TTree *treeData = (TTree*)fIn->Get("data");

  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  printf("initializing tree\n");
  MCTnPTriggerEff tree(treeData);  printf("...done\n");

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
  if(iEff == MuIDEff || iEff == MuQualEff || iEff == L1L2Eff || iEff == L3Eff || iEff == TMandHLT2Eff){
    for(int iEta = 0; iEta < kNbEtaBins[iEff]; iEta++){
      sprintf(name, "fit%s_%s_pt_etaBin%d", effName[iEff], effSampleName[iEffSample], iEta);
      fMuEff_pT[iEff][iEffSample][iEta] = (TF1 *) gDirectory->Get(name);
      printf("fitfunction for %s, etaBin %d: %p\n", effName[iEff], iEta, fMuEff_pT[iEff][iEffSample][iEta]);
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

  // //======================================================
  // //1D distributions: 
  // //======================================================
  trigEff_pT = new TEfficiency("trigEff_pT", ";p_{T} [GeV/c]", eff::kNbPTMaxBins, eff::pTRange[0]); 
  trigEff_y = new TEfficiency("trigEff_y", ";y", 2*eff::kNbRapForPTBins, eff::rapRange); 
  trigEff_phi = new TEfficiency("trigEff_phi", ";#phi", eff::nBinsPhi1D, eff::phiMin, eff::phiMax); 

  // //======================================================
  // //2D distributions: 
  // //======================================================
  trigEff2D_pT_rap = new TEfficiency("trigEff2D_pT_rap", ";|y|;p_{T} [GeV/c]", eff::kNbRapForPTBins, eff::rapForPTRange, eff::kNbPTMaxBins, eff::pTRange[0]);
  trigEff2D_pT_rapNP = new TEfficiency("trigEff2D_pT_rapNP", ";|y|;p_{T} [GeV/c]", 2*eff::kNbRapForPTBins, eff::rapRange, eff::kNbPTMaxBins, eff::pTRange[0]);
//   //======================================================
//   //polarization histograms:
//   //======================================================

  //histos for neg. and pos. rapidity separately:
  for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*eff::kNbRapBins; iRapBin++){
      Int_t rapIndex = eff::kNbRapBins - iRapBin;
      if(iRapBin >= eff::kNbRapBins) rapIndex = iRapBin - eff::kNbRapBins + 1;
//       printf("iRapBin %d uses definitions of rapID %d\n", iRapBin, rapIndex);

      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[rapIndex]+1; iPTBin++){

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
	sprintf(name, "trigEff2D_Onia_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT, eff::cosTMin, eff::cosTMax,
								 eff::kNbBinsPhiPol, eff::phiPolMin, eff::phiPolMax);
	//(2) 4-folding in phi
	sprintf(name, "trigEff2D_Onia_phiFolded_%s_pT%d_rap%d", eff::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", eff::frameLabel[iFrame], eff::frameLabel[iFrame]);
	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TEfficiency(name, title, eff::kNbBinsCosT/2, eff::cosTMin, eff::cosTMax,
								     eff::kNbBinsPhiPol, 0., 90.);
     }
    }
  }


//   //efficiencies vs. lab coordinates
//   Int_t nBinsDeltaEta = 10, nBinsDeltaPhi = 9;
//   Double_t deltaEtaMin = 0., deltaEtaMax = 2.5;
//   Double_t deltaPhiMin = 0., deltaPhiMax = 180.;
//   Int_t nBinsDeltaR = 32;
//   Double_t deltaRMin = 0., deltaRMax = 3.2;
//   Int_t nBinsDistM2 = 20;
//   Double_t distMin = 0., distMax = 200.;
//   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
//     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

//       //1.) as a function of deltaPhi vs deltaEta
//       //all events
//       sprintf(name, "hGen2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
//       hGen2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
//       hGen2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //total efficiencies
//       sprintf(name, "totEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
//       totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
//       totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO events
//       sprintf(name, "recoEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
//       recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
//       recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO + TRIG events
//       sprintf(name, "trigEff2D_deltaPhiVsDeltaEta_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#eta;#Delta#phi [deg]");
//       trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin] = new TH2D(name, title, nBinsDeltaEta, deltaEtaMin, deltaEtaMax, nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
//       trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //2.) as a function of deltaR
//       //all events
//       sprintf(name, "hGen_deltaR_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#DeltaR");
//       hGen_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
//       hGen_deltaR_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //total efficiencies
//       sprintf(name, "totEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#DeltaR");
//       totEff_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
//       totEff_deltaR_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO events
//       sprintf(name, "recoEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#DeltaR");
//       recoEff_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
//       recoEff_deltaR_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO + TRIG events
//       sprintf(name, "trigEff_deltaR_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#DeltaR");
//       trigEff_deltaR_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
//       trigEff_deltaR_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //3.) as a function of deltaRM2
//       //all events
//       sprintf(name, "hGen_deltaRM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#DeltaR(M2)");
//       hGen_deltaRM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
//       hGen_deltaRM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //total efficiencies
//       sprintf(name, "totEff_deltaRM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#DeltaR(M2)");
//       totEff_deltaRM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
//       totEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO events
//       sprintf(name, "recoEff_deltaRM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#DeltaR(M2)");
//       recoEff_deltaRM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
//       recoEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO + TRIG events
//       sprintf(name, "trigEff_deltaRM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#DeltaR(M2)");
//       trigEff_deltaRM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDeltaR, deltaRMin, deltaRMax);
//       trigEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //4.) as a function of deltaPhiM2
//       //all events
//       sprintf(name, "hGen_deltaPhiM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#phi(M2)");
//       hGen_deltaPhiM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
//       hGen_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //total efficiencies
//       sprintf(name, "totEff_deltaPhiM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#phi(M2)");
//       totEff_deltaPhiM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
//       totEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO events
//       sprintf(name, "recoEff_deltaPhiM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#phi(M2)");
//       recoEff_deltaPhiM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
//       recoEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO + TRIG events
//       sprintf(name, "trigEff_deltaPhiM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#phi(M2)");
//       trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaPhi, deltaPhiMin, deltaPhiMax);
//       trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //5.) as a function of deltaEtaM2
//       //all events
//       sprintf(name, "hGen_deltaEtaM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#eta(M2)");
//       hGen_deltaEtaM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaEta, deltaEtaMin, deltaEtaMax);
//       hGen_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //total efficiencies
//       sprintf(name, "totEff_deltaEtaM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#eta(M2)");
//       totEff_deltaEtaM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaEta, deltaEtaMin, deltaEtaMax);
//       totEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO events
//       sprintf(name, "recoEff_deltaEtaM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#eta(M2)");
//       recoEff_deltaEtaM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaEta, deltaEtaMin, deltaEtaMax);
//       recoEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO + TRIG events
//       sprintf(name, "trigEff_deltaEtaM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";#Delta#phi(M2)");
//       trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, 2*nBinsDeltaEta, deltaEtaMin, deltaEtaMax);
//       trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //6.) as a function of distM2
//       //all events
//       sprintf(name, "hGen_distM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";dist(M2) [cm]");
//       hGen_distM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDistM2, distMin, distMax);
//       hGen_distM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //total efficiencies
//       sprintf(name, "totEff_distM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";dist(M2) [cm]");
//       totEff_distM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDistM2, distMin, distMax);
//       totEff_distM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO events
//       sprintf(name, "recoEff_distM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";dist(M2) [cm]");
//       recoEff_distM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDistM2, distMin, distMax);
//       recoEff_distM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//       //RECO + TRIG events
//       sprintf(name, "trigEff_distM2_pT%d_rap%d", iPTBin, iRapBin);
//       sprintf(title, ";dist(M2) [cm]");
//       trigEff_distM2_pT_rap[iPTBin][iRapBin] = new TH1D(name, title, nBinsDistM2, distMin, distMax);
//       trigEff_distM2_pT_rap[iPTBin][iRapBin]->Sumw2();

//     }
//   }

}

//==========================================
void DivideHistos(){

//   //given that the "reco" histograms were filled with Possonian errors,
//   //we can calculate the binomial errors
//   //1.) ID histograms
//   trigEff_pT->Divide(trigEff_pT, recoEff_pT, 1., 1., "B");
//   trigEff_y->Divide(trigEff_y, recoEff_y, 1., 1., "B");
//   trigEff_phi->Divide(trigEff_phi, recoEff_phi, 1., 1., "B");

//   recoEff_pT->Divide(recoEff_pT, hGen_pT, 1., 1., "B"); 
//   recoEff_y->Divide(recoEff_y, hGen_y, 1., 1., "B");
//   recoEff_phi->Divide(recoEff_phi, hGen_phi, 1., 1., "B");

//   // trigEff_pT->Divide(trigEff_pT, hGen_pT, 1., 1., "B");
//   // trigEff_y->Divide(trigEff_y, hGen_y, 1., 1., "B");
//   // trigEff_phi->Divide(trigEff_phi, hGen_phi, 1., 1., "B");

//   totEff_pT->Divide(totEff_pT, hGen_pT, 1., 1., "B");
//   totEff_y->Divide(totEff_y, hGen_y, 1., 1., "B");
//   totEff_phi->Divide(totEff_phi, hGen_phi, 1., 1., "B");

//   //2.) 2D histograms
//   trigEff2D_pT_rapNP->Divide(trigEff2D_pT_rapNP, recoEff2D_pT_rapNP, 1., 1., "B");
//   trigEff2D_pT_rap->Divide(trigEff2D_pT_rap, recoEff2D_pT_rap);

//   recoEff2D_pT_rapNP->Divide(recoEff2D_pT_rapNP, hGen2D_pT_rapNP, 1., 1., "B");
//   recoEff2D_pT_rap->Divide(recoEff2D_pT_rap, hGen2D_pT_rap);
//   // trigEff2D_pT_rapNP->Divide(trigEff2D_pT_rapNP, hGen2D_pT_rapNP, 1., 1., "B");
//   // trigEff2D_pT_rap->Divide(trigEff2D_pT_rap, hGen2D_pT_rap);
//   totEff2D_pT_rapNP->Divide(totEff2D_pT_rapNP, hGen2D_pT_rapNP, 1., 1., "B");
//   totEff2D_pT_rap->Divide(totEff2D_pT_rap, hGen2D_pT_rap);

//   //3.) polarization histograms:
//   for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
//     for(int iRapBin = 0; iRapBin < 2*eff::kNbRapBins; iRapBin++){
//       Int_t rapIndex = eff::kNbRapBins - iRapBin;
//       if(iRapBin >= eff::kNbRapBins) rapIndex = iRapBin - eff::kNbRapBins + 1;
// //       printf("iRapBin %d uses definitions of rapID %d\n", iRapBin, rapIndex);

//       for(int iPTBin = 0; iPTBin < eff::kNbPTBins[rapIndex]+1; iPTBin++){
// 	trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Divide(trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], 1., 1., "B");
// 	recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Divide(recoEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], 1., 1., "B");
// 	//	trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Divide(trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], 1., 1., "B");
// 	totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Divide(totEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin], 1., 1., "B");
//       }
//     }
    
//     for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
//       for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
// 	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Divide(trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin], recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin], 1., 1., "B");
// 	recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Divide(recoEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin], 1., 1., "B");
// 	// trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Divide(trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin], 1., 1., "B");
// 	totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Divide(totEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin], 1., 1., "B");
      
// 	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Divide(trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], 1., 1., "B");
// 	recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Divide(recoEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], 1., 1., "B");
// 	// trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Divide(trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], 1., 1., "B");
// 	totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Divide(totEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin], 1., 1., "B");
//       }
//     }
//   }

//   //lab-frame variables
//   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
//     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){

//       trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Divide(trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Divide(recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], hGen2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       // trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Divide(trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], hGen2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Divide(totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], hGen2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin], 1., 1., "B");

//       trigEff_deltaR_pT_rap[iPTBin][iRapBin]->Divide(trigEff_deltaR_pT_rap[iPTBin][iRapBin], recoEff_deltaR_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       recoEff_deltaR_pT_rap[iPTBin][iRapBin]->Divide(recoEff_deltaR_pT_rap[iPTBin][iRapBin], hGen_deltaR_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       // trigEff_deltaR_pT_rap[iPTBin][iRapBin]->Divide(trigEff_deltaR_pT_rap[iPTBin][iRapBin], hGen_deltaR_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       totEff_deltaR_pT_rap[iPTBin][iRapBin]->Divide(totEff_deltaR_pT_rap[iPTBin][iRapBin], hGen_deltaR_pT_rap[iPTBin][iRapBin], 1., 1., "B");

//       trigEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Divide(trigEff_deltaRM2_pT_rap[iPTBin][iRapBin], recoEff_deltaRM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       recoEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Divide(recoEff_deltaRM2_pT_rap[iPTBin][iRapBin], hGen_deltaRM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       //trigEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Divide(trigEff_deltaRM2_pT_rap[iPTBin][iRapBin], hGen_deltaRM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       totEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Divide(totEff_deltaRM2_pT_rap[iPTBin][iRapBin], hGen_deltaRM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");

//       trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Divide(trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin], recoEff_deltaPhiM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       recoEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Divide(recoEff_deltaPhiM2_pT_rap[iPTBin][iRapBin], hGen_deltaPhiM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       //      trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Divide(trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin], hGen_deltaPhiM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       totEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Divide(totEff_deltaPhiM2_pT_rap[iPTBin][iRapBin], hGen_deltaPhiM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");

//       trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Divide(trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin], recoEff_deltaEtaM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       recoEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Divide(recoEff_deltaEtaM2_pT_rap[iPTBin][iRapBin], hGen_deltaEtaM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       // trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Divide(trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin], hGen_deltaEtaM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       totEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Divide(totEff_deltaEtaM2_pT_rap[iPTBin][iRapBin], hGen_deltaEtaM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");

//       trigEff_distM2_pT_rap[iPTBin][iRapBin]->Divide(trigEff_distM2_pT_rap[iPTBin][iRapBin], recoEff_distM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       recoEff_distM2_pT_rap[iPTBin][iRapBin]->Divide(recoEff_distM2_pT_rap[iPTBin][iRapBin], hGen_distM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       // trigEff_distM2_pT_rap[iPTBin][iRapBin]->Divide(trigEff_distM2_pT_rap[iPTBin][iRapBin], hGen_distM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//       totEff_distM2_pT_rap[iPTBin][iRapBin]->Divide(totEff_distM2_pT_rap[iPTBin][iRapBin], hGen_distM2_pT_rap[iPTBin][iRapBin], 1., 1., "B");
//     }
//   }
}

//==========================================
void WriteHistos(){

//   for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
  for(int iFrame = 0; iFrame < 2; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*eff::kNbRapBins; iRapBin++){
      Int_t rapIndex = eff::kNbRapBins - iRapBin;
      if(iRapBin >= eff::kNbRapBins) rapIndex = iRapBin - eff::kNbRapBins + 1;
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[rapIndex]+1; iPTBin++){
	trigEff2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

//   for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
  for(int iFrame = 0; iFrame < 2; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	trigEff2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

//   for(int iFrame = 0; iFrame < eff::kNbFrames; iFrame++){
  for(int iFrame = 0; iFrame < 2; iFrame++){
    for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
	trigEff2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

//   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
//     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
//       totEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Write();
//       recoEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Write();
//       trigEff2D_deltaPhiVsDeltaEta_pT_rap[iPTBin][iRapBin]->Write();
//     }
//   }
//   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
//     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
//       totEff_deltaR_pT_rap[iPTBin][iRapBin]->Write();
//       recoEff_deltaR_pT_rap[iPTBin][iRapBin]->Write();
//       trigEff_deltaR_pT_rap[iPTBin][iRapBin]->Write();
//     }
//   }
//   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
//     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
//       totEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Write();
//       recoEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Write();
//       trigEff_deltaRM2_pT_rap[iPTBin][iRapBin]->Write();
//     }
//   }
//   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
//     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
//       totEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Write();
//       recoEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Write();
//       trigEff_deltaPhiM2_pT_rap[iPTBin][iRapBin]->Write();
//     }
//   }

//   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
//     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
//       totEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Write();
//       recoEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Write();
//       trigEff_deltaEtaM2_pT_rap[iPTBin][iRapBin]->Write();
//     }
//   }

//   for(int iRapBin = 0; iRapBin < eff::kNbRapForPTBins+1; iRapBin++){
//     for(int iPTBin = 0; iPTBin < eff::kNbPTBins[iRapBin]+1; iPTBin++){
//       totEff_distM2_pT_rap[iPTBin][iRapBin]->Write();
//       recoEff_distM2_pT_rap[iPTBin][iRapBin]->Write();
//       trigEff_distM2_pT_rap[iPTBin][iRapBin]->Write();
//     }
//   }

//   //efficiencies (in histo form)
//   recoEff_pT->Write();   recoEff_y->Write();   recoEff_phi->Write();
//   trigEff_pT->Write();   trigEff_y->Write();   trigEff_phi->Write();
//   totEff_pT->Write();   totEff_y->Write();   totEff_phi->Write();
//   recoEff2D_pT_rapNP->Write();  recoEff2D_pT_rap->Write();
//   trigEff2D_pT_rapNP->Write();  trigEff2D_pT_rap->Write();
//   totEff2D_pT_rapNP->Write();  totEff2D_pT_rap->Write();
  
  trigEff_pT->Write();   trigEff_y->Write();   trigEff_phi->Write();
  trigEff2D_pT_rap->Write();
  trigEff2D_pT_rapNP->Write();

}
