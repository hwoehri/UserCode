#include "../interface/rootIncludes.inc"
#include "GeomAcc.C"

void BookHistos(Char_t *oniaLabel);
void WriteHistos();
//======================================
void runGeomAcc(Char_t *fileNameOut = "geomAcc_WithFSR_kinWeighted.root",
		Char_t *fileNameIn = "jpsiGun_WithFSR_Tree.root",
		Bool_t applySmearing = kTRUE,
		Bool_t rejectCowboys = kTRUE,
		//Char_t *fileNameOut = "geomAcc.root",
		//		Char_t *fileNameIn = "jpsiGun_Tree.root",
		Char_t *oniaLabel = "J/#psi" //"J/#psi", "#psi'", "Ups(1S)", "Ups(2S)", "Ups(3S)"
		){

  TFile *fIn = new TFile(fileNameIn);
  TTree *treeData = (TTree*)fIn->Get("QQbarGenTree");

  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  printf("initializing tree\n");
  GeomAcc tree(treeData);
  printf("...done\n");
  BookHistos(oniaLabel);

  tree.Loop(applySmearing, rejectCowboys);
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

  hGenPtRap = new TH2D("hGenPtRap", "", 50, -2.5, 2.5, 500, 0., 50.);
  hGenCutPtRap = new TH2D("hGenCutPtRap", "", 50, -2.5, 2.5, 500, 0., 50.);

  Char_t name[100], title[100];
  //histos for neg. and pos. rapidity separately:
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins; iRapBin++){
      Int_t rapIndex = jpsi::kNbRapBins - iRapBin;
      if(iRapBin >= jpsi::kNbRapBins) rapIndex = iRapBin - jpsi::kNbRapBins + 1;
//       printf("iRapBin %d uses definitions of rapID %d\n", iRapBin, rapIndex);

      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[rapIndex]+1; iPTBin++){

	//all events
	sprintf(name, "hGen_%s_pT%d_NPrap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	hGen_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								 jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	hGen_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();
	//accepted events
	sprintf(name, "hGenCut_%s_pT%d_NPrap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	hGenCut_pT_rapNP[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								 jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	hGenCut_pT_rapNP[iFrame][iPTBin][iRapBin]->Sumw2();
      }
    }
  }
  //histos taking together +y and -y:
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	//(1) no folding in phi and cosTheta
	//all events
	sprintf(name, "hGen_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	hGen_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								 jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	hGen_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
	//accepted events
	sprintf(name, "hGenCut_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	hGenCut_pT_rap[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								 jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	hGenCut_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
	//(2) 4-folding in phi
	//all events
	sprintf(name, "hGen_phiFolded_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	hGen_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								  jpsi::kNbBinsPhiPol/4, 0., 90.);
	hGen_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();
	//accepted events
	sprintf(name, "hGenCut_phiFolded_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	hGenCut_pT_rap_phiFolded[iFrame][iPTBin][iRapBin] = new TH2D(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								     jpsi::kNbBinsPhiPol/4, 0., 90.);
	hGenCut_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Sumw2();

      }
    }
  }
  for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	//
	sprintf(name, "hMass_pT%d_rap%d", iPTBin, iRapBin);
	sprintf(title, ";Mass [GeV]");
	hMass[iPTBin][iRapBin] = new TH1D(name, title, nBinsMass, massMin, massMax);
	hMass[iPTBin][iRapBin]->Sumw2();
	//
	sprintf(name, "hMass_Smeared_pT%d_rap%d", iPTBin, iRapBin);
	sprintf(title, ";Mass [GeV]");
	hMass_Smeared[iPTBin][iRapBin] = new TH1D(name, title, nBinsMass, massMin, massMax);
	hMass_Smeared[iPTBin][iRapBin]->Sumw2();
    }
  }
}

//==========================================
void WriteHistos(){

  hGenPtRap->Write();
  hGenCutPtRap->Write();

  for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
      //
      hMass[iPTBin][iRapBin]->Write();
      hMass_Smeared[iPTBin][iRapBin]->Write();
    }
  }

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins; iRapBin++){
      Int_t rapIndex = jpsi::kNbRapBins - iRapBin;
      if(iRapBin >= jpsi::kNbRapBins) rapIndex = iRapBin - jpsi::kNbRapBins + 1;
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[rapIndex]+1; iPTBin++){
	hGen_pT_rapNP[iFrame][iPTBin][iRapBin]->Write();
	hGenCut_pT_rapNP[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	hGen_pT_rap[iFrame][iPTBin][iRapBin]->Write();
	hGenCut_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	hGen_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
	hGenCut_pT_rap_phiFolded[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }
}
