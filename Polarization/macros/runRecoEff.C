#include "../interface/rootIncludes.inc"
//#include "RecoEff.C"
#include "TrigEff.C"

void BookHistos(Char_t *oniaLabel);
void WriteHistos();
//======================================
void runRecoEff(Char_t *fileNameOut = "recoEff_7Feb2011.root",
		Bool_t useMCWeight = kFALSE, //used for NP sample which is composed of 4 simulations
		Int_t selDimuType = 4, //0...only GG, 1... only GT, 2... only TT, 3...GG+GT, 4...GG+GT+TT		
		Char_t *nameDataSet = "data", //"data" or "recoData"
 		Char_t *fileNameIn = "JPsiToMuMu_Fall10-START38_V12-v1-Onia2MuMu-v6-WithAllMCEvents_merged.root",
//  		Char_t *fileNameIn = "TTree_NPmerged_Bx_toJpsi_toMuMu.root",
		Char_t *oniaLabel = "J/#psi" //"J/#psi", "#psi'", "Ups(1S)", "Ups(2S)", "Ups(3S)"
		){

  TFile *fIn = new TFile(fileNameIn);
  TTree *treeData = (TTree*)fIn->Get(nameDataSet);

  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  printf("initializing tree\n");
  RecoEff tree(treeData);
  printf("...done\n");
  BookHistos(oniaLabel);

  tree.Loop(selDimuType, useMCWeight);
  WriteHistos();
  fOut->Close();
}

//==========================================
void BookHistos(Char_t *oniaLabel){

  //initialise the "massMuOnia" variable
  std::map<std::string, double> fOniaToMass;
  fOniaToMass["J/#psi"] = 3.0969;
  fOniaToMass["#psi'"] = 3.686;
  fOniaToMass["Ups(1S)"] = 9.460;
  fOniaToMass["Ups(2S)"] = 10.023;
  fOniaToMass["Ups(3S)"] = 10.355;

  std::map<std::string, double>::iterator iter = fOniaToMass.find(oniaLabel);
  if(iter != fOniaToMass.end())
    massMuOnia = iter->second;
  printf("will use a mass of %1.3f GeV for the %s\n", massMuOnia, oniaLabel);
  //===================================

  //mass
  Int_t nBinsMass = 80;
  Double_t massMin = 8.0, massMax = 12.0;
  if(strncmp(oniaLabel, "J/#psi", 6) == 0){
    nBinsMass = 160;
    massMin = 2.5;
    massMax = 4.1;
  }
  //pt
  Int_t nBinsPt = 500;
  Double_t pTMin = 0., pTMaxOnia = 50.;
  //rap
  Int_t nBinsRap = 100;
  Double_t rapMin = -2.5, rapMax = 2.5;

//   hGenPtRap = new TH2D("hGenPtRap", "", 50, -2.5, 2.5, 500, 0., 50.);
//   hGenCutPtRap = new TH2D("hGenCutPtRap", "", 50, -2.5, 2.5, 500, 0., 50.);

  Char_t name[100], title[100];
  for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	//
      //Mass:
      sprintf(name, "Reco_Onia_mass_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";M [GeV/c^{2}]");
      Reco_Onia_mass[iPTBin][iRapBin] = new TH1D(name, title, nBinsMass, massMin, massMax);
      Reco_Onia_mass[iPTBin][iRapBin]->Sumw2();
      //phi-lab:
      sprintf(name, "Reco_Onia_phi_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";%s #phi (rad)", oniaLabel);
      Reco_Onia_phi[iPTBin][iRapBin] = new TH1D(name, title, jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
      Reco_Onia_phi[iPTBin][iRapBin]->Sumw2();
    }
  }
  for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    //pT
    sprintf(name, "Reco_Onia_pt_rap%d", iRapBin);
    sprintf(title, ";%s p_{T} [GeV/c]", oniaLabel);
    Reco_Onia_pt[iRapBin]  = new TH1D(name, title, nBinsPt,pTMin,pTMaxOnia);
    Reco_Onia_pt[iRapBin]->Sumw2();
  }
  //2D Onia histos:
  sprintf(name, "Reco_Onia_rap_pt");
  sprintf(title, ";%s y;p_{T} [GeV/c]", oniaLabel);
  Reco_Onia_rap_pT = new TH2D(name, title, nBinsRap,rapMin,rapMax, nBinsPt,pTMin,pTMaxOnia);
  Reco_Onia_rap_pT->Sumw2();


  //histos for neg. and pos. rapidity separately:
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins; iRapBin++){
      Int_t rapIndex = jpsi::kNbRapBins - iRapBin;
      if(iRapBin >= jpsi::kNbRapBins) rapIndex = iRapBin - jpsi::kNbRapBins + 1;
//       printf("iRapBin %d uses definitions of rapID %d\n", iRapBin, rapIndex);

      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[rapIndex]+1; iPTBin++){

	//all events
	sprintf(name, "hGen2D_Onia_%s_pT%d_NPrap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								 jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();
	//accepted events
	sprintf(name, "Reco2D_Onia_%s_pT%d_NPrap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	Reco2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	Reco2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();
      }
    }
  }
  //histos taking together +y and -y:
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	//(1) no folding in phi and cosTheta
	//all events
	sprintf(name, "hGen2D_Onia_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								 jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
	//accepted events
	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								 jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
	//(2) 4-folding in phi
	//all events
	sprintf(name, "hGen2D_Onia_phiFolded_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								  jpsi::kNbBinsPhiPol/4, 0., 90.);
	hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();
	//accepted events
	sprintf(name, "Reco2D_Onia_phiFolded_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	Reco2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								     jpsi::kNbBinsPhiPol/4, 0., 90.);
	Reco2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();
      }
    }
  }
}

//==========================================
void WriteHistos(){

//   hGenPtRap->Write();
//   hGenCutPtRap->Write();

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins; iRapBin++){
      Int_t rapIndex = jpsi::kNbRapBins - iRapBin;
      if(iRapBin >= jpsi::kNbRapBins) rapIndex = iRapBin - jpsi::kNbRapBins + 1;
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[rapIndex]+1; iPTBin++){
	hGen2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Write();
	Reco2D_pol_pT_rapNP[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	hGen2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
	Reco2D_pol_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }
}
