#include "ProjectQQ.C"

void BookHistos();
void WriteHistos(Bool_t normalize, Double_t intLumi, Char_t *fNameOut);
//==========================================
void runProjectQQ(Char_t *fNameOut = "histos_jpsi_MC_900GeV_STARTUP.root", //fileName where the histos will be stored
		  Bool_t matchMC = kFALSE,//set to kTRUE for J/psi MC sample and to kFALSE for BG MC and real data
		  Bool_t removeQQ = kFALSE, //for MB/ppMuXLoose samples set to kTRUE; for J/psi MC and real data: set to kFALSE
		  Bool_t normalize = kFALSE, //set to true *only* if you want to normalise the histos; in standard mode: should be kFALSE
		  Double_t intLumi = 16.675//only relevant if "normalise"=kTRUE; 2.36 TeV: 13557/813=16.675; 900 GeV: 6724./404=16.644
		  ){

  ProjectQQ tree;

  BookHistos();

  tree.Loop(removeQQ, matchMC);
  WriteHistos(normalize, intLumi, fNameOut);

}
//==========================================
void BookHistos(){

  //mass
//   Int_t nBinsMass = 20;
  Int_t nBinsMass = 1200;
  Double_t massMin = 0.0, massMax = 12.0;
  //pt
  Int_t nBinsPT = 20;
  Double_t pTMin = 0., pTMax = 10.;
  //rap
  Int_t nBinsRap = 60;
  Double_t rapMin = -3.0, rapMax = 3.0;

  double ptbins[] = {0.0, 0.5, 1.0, 1.5,
                     2.0, 2.5, 3.0, 3.5,
                     4.0, 4.5, 5.0, 5.5,
                     6.0, 6.5, 7.0, 7.5,
                     8.0, 8.5, 9.0, 9.5,
                     10.0, 10.5, 11.0, 11.5,
                     12.0, 12.5, 13.0, 13.5,
                     14.0, 14.5, 15.0, 15.5,
                     16.0, 17.0, 18.0, 19.0,
                     20.0, 22.0, 24.0, 26.0,
                     28.0, 30.0, 35.0, 40.0,
                     50.0, 70.0, 100.};
  int nbinspt = sizeof(ptbins)/sizeof(double)-1;


  Char_t name[100], title[300];
  for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
    for(int iCat = 0; iCat < kNbCath; iCat++){
      for (int iCut = 0; iCut < kNbCuts; iCut++) {

	//mass histos
	sprintf(name, "hReco_mass_%s_%s_%s", chargeName[iCharge], oniaCatName[iCat], muCutName[iCut]);
	sprintf(title, ";M [GeV/c^{2}];dN/dM [per %1.0f MeV/c^{2}]", 
		1000.* (massMax-massMin)/nBinsMass);
	hReco_mass[iCharge][iCat][iCut] = new TH1F(name, title, nBinsMass, massMin, massMax);
	hReco_mass[iCharge][iCat][iCut]->Sumw2();
	
	for(int iSet = 0; iSet < kNbDimuonSet; iSet++){
	  //pT histos
	  sprintf(name, "hReco_pT_%d_%s_%s_%s", iSet, chargeName[iCharge], oniaCatName[iCat], muCutName[iCut]);
	  if(iSet == 0)
	    sprintf(title, ";p_{T}(J/#psi) [GeV/c];dN/dp_{T} [per %1.0f MeV/c]", 1000.* (pTMax-pTMin)/nBinsPT);
	  else
	    sprintf(title, ";p_{T} [GeV/c];dN/dp_{T} [per %1.0f MeV/c]", 1000.* (pTMax-pTMin)/nBinsPT);
	  hReco_pT[iSet][iCharge][iCat][iCut] = new TH1F(name, title, nBinsPT, pTMin, pTMax);
	  hReco_pT[iSet][iCharge][iCat][iCut]->Sumw2();
	  //rap histos
	  sprintf(name, "hReco_rap_%d_%s_%s_%s", iSet, chargeName[iCharge], oniaCatName[iCat], muCutName[iCut]);
	  if(iSet == 0)
	    sprintf(title, ";y(J/#psi);dN/dy");
	  else
	    sprintf(title, ";y;dN/dy");
	  hReco_rap[iSet][iCharge][iCat][iCut] = new TH1F(name, title, nBinsRap, rapMin, rapMax);
	  hReco_rap[iSet][iCharge][iCat][iCut]->Sumw2();

	  //eta(mu1) vs eta(mu2)
	  sprintf(name, "hReco_eta1_eta2_%d_%s_%s_%s", iSet, chargeName[iCharge], oniaCatName[iCat], muCutName[iCut]);
	  hReco_eta1_eta2[iSet][iCharge][iCat][iCut] = new TH2F(name, ";#eta(#mu^{+});#eta(#mu^{-})", nBinsRap, rapMin, rapMax, nBinsRap, rapMin, rapMax);

	  if(iSet == 0){
	    //multiplicity
	    sprintf(name, "hMult_QQ_%s_%s_%s", chargeName[iCharge], oniaCatName[iCat], muCutName[iCut]);
	    hMult_QQ[iCharge][iCat][iCut] = new TH1F(name, ";multiplicity", 10, 0., 10.);
	  }
	  if(iCut == 0){

	    //dimuon vertexing probability
	    sprintf(name, "hDimuVtxProb_%d_%s_%s", iSet, chargeName[iCharge], oniaCatName[iCat]);
	    hDimuVtxProb[iSet][iCharge][iCat] = new TH1F(name, ";#chi^{2} prob. of dimuon vertex", 2000, 0., 1.);
	  }
	} //iSet
      } //iCut

	//mass histos: best QQbar
      sprintf(name, "hRecoBest_mass_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
      sprintf(title, ";M [GeV/c^{2}];dN/dM [per %1.0f MeV/c^{2}]", 
	      1000.* (massMax-massMin)/nBinsMass);
      hRecoBest_mass[iCharge][iCat] = new TH1F(name, title, nBinsMass, massMin, massMax);
      hRecoBest_mass[iCharge][iCat]->Sumw2();

      for(int iSet = 0; iSet < kNbDimuonSet; iSet++){

	//pT histos: best QQbar
	sprintf(name, "hRecoBest_pT_%d_%s_%s", iSet, chargeName[iCharge], oniaCatName[iCat]);
	sprintf(title, ";p_{T}(J/#psi) [GeV/c];dN/dp_{T} [per %1.0f MeV/c]", 
		1000.* (pTMax-pTMin)/nBinsPT);
	hRecoBest_pT[iSet][iCharge][iCat] = new TH1F(name, title, nBinsPT, pTMin, pTMax);
	hRecoBest_pT[iSet][iCharge][iCat]->Sumw2();
	//rap histos: best QQbar
	sprintf(name, "hRecoBest_rap_%d_%s_%s", iSet, chargeName[iCharge], oniaCatName[iCat]);
	sprintf(title, ";y(J/#psi);dN/dy");
	hRecoBest_rap[iSet][iCharge][iCat] = new TH1F(name, title, nBinsRap, rapMin, rapMax);
	hRecoBest_rap[iSet][iCharge][iCat]->Sumw2();
      }
    } // iCat
  } // iCharge

  //counter for events where no "bestQQ" could be found:
  sprintf(name, "hRecoNoBest");
  hRecoNoBest = new TH1F(name, ";J/#psi category", kNbCath, 0., kNbCath);

  //single muon histograms:
  for(int iSet = 0; iSet < kNbDimuonSet; iSet++){
    for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
      for(int iCat = 0; iCat < kNbMuCath; iCat++){

	//eta
	sprintf(name, "hMuon_eta_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hMuon_eta[iSet][iCharge][iCat] = new TH1F(name, ";#eta(#mu)", nBinsRap, rapMin, rapMax);
	hMuon_eta[iSet][iCharge][iCat]->Sumw2();

	//pT
	sprintf(name, "hMuon_pT_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hMuon_pT[iSet][iCharge][iCat] = new TH1F(name, ";p_{T}(#mu) [GeV/c]", nBinsPT, pTMin, pTMax);
	hMuon_pT[iSet][iCharge][iCat]->Sumw2();

	//pT vs eta
	sprintf(name, "hMuon_pTvsEta_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hMuon_pT_eta[iSet][iCharge][iCat] = new TH2F(name, ";#eta(#mu);p_{T}(#mu) [GeV/c]", nBinsRap, rapMin, rapMax, nBinsPT, pTMin, pTMax);
	hMuon_pT_eta[iSet][iCharge][iCat]->Sumw2();

	//chi2 of global fit
	if(iCat == 0){//global muons
	  sprintf(name, "hChi2GlobalFit_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	  hChi2GlobalFit[iSet][iCharge] = new TH1F(name, ";#chi^{2}/ndf (global fit)", 200, 0, 20.);
	}
	//chi2 of tracker fit
	sprintf(name, "hChi2TrackerFit_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hChi2TrackerFit[iSet][iCharge][iCat] = new TH1F(name, ";#chi^{2}/ndf (tracker fit)", 200, 0, 20.);

	//#hits in silicon tracker
	sprintf(name, "hNHitsSilicon_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hNHitsSilicon[iSet][iCharge][iCat] = new TH1F(name, ";#hits in silicon tracker", 40, 0, 40.);

	if(iCat == 1){//tracker muons
	  //hTM2DCompatibilityTight
	  sprintf(name, "hTM2DCompatibilityTight_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	  hTM2DCompatibilityTight[iSet][iCharge] = new TH1F(name, ";TM2DCompatibilityTight", 2, 0, 2);
	  //hTMLastStationOptimizedLowPtLoose
	  sprintf(name, "hTMLastStationOptimizedLowPtLoose_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	  hTMLastStationOptimizedLowPtLoose[iSet][iCharge] = new TH1F(name, ";TMLastStationOptimizedLowPtLoose", 2, 0, 2);
	  //hTrackerMuonArbitrated
	  sprintf(name, "hTrackerMuonArbitrated_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	  hTrackerMuonArbitrated[iSet][iCharge] = new TH1F(name, ";TrackerMuonArbitrated", 2, 0, 2);
	}
	//calo compatibility
	sprintf(name, "hCaloCompatibility_%d_%s_%s", iSet, chargeName[iCharge], muCatName[iCat]);
	hCaloCompatibility[iSet][iCharge][iCat] = new TH1F(name, ";CaloCompatibility", 100, 0., 1.);
      }//cat
    }//charge
  }//set
}

//==========================================
void WriteHistos(Bool_t normalize, Double_t intLumi, Char_t *fNameOut){

  TFile *fOut = new TFile(fNameOut, "RECREATE");
  for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
    for(int iCat = 0; iCat < kNbCath; iCat++){
      for(int iCut = 0; iCut < kNbCuts; iCut++){
	//normalise the histo to 1 nb-1:
	if(normalize){
	  hReco_mass[iCharge][iCat][iCut]->Scale(1./intLumi);
	  for(int iSet = 0; iSet < kNbDimuonSet; iSet++){
	    hReco_pT[iSet][iCharge][iCat][iCut]->Scale(1./intLumi);
	    hReco_rap[iSet][iCharge][iCat][iCut]->Scale(1./intLumi);
	  }
	}
	hReco_mass[iCharge][iCat][iCut]->Write();
	for(int iSet = 0; iSet < kNbDimuonSet; iSet++){
	  hReco_pT[iSet][iCharge][iCat][iCut]->Write();
	  hReco_rap[iSet][iCharge][iCat][iCut]->Write();
	  if(iCat == 0 || iCat == 2)
	    hReco_eta1_eta2[iSet][iCharge][iCat][iCut]->Write();
	}
	hMult_QQ[iCharge][iCat][iCut]->Write();
      }//cut
      
      //normalise the histo to 1 nb-1:
      if(normalize){
	hRecoBest_mass[iCharge][iCat]->Scale(1./intLumi);
	for(int iSet = 0; iSet < kNbDimuonSet; iSet++){
	  hRecoBest_pT[iSet][iCharge][iCat]->Scale(1./intLumi);
	  hRecoBest_rap[iSet][iCharge][iCat]->Scale(1./intLumi);
	}
      }
      //mass
      hRecoBest_mass[iCharge][iCat]->Write();
      for(int iSet = 0; iSet < kNbDimuonSet; iSet++){
	hRecoBest_pT[iSet][iCharge][iCat]->Write();
	hRecoBest_rap[iSet][iCharge][iCat]->Write();
      }

      for(int iSet = 0; iSet < kNbDimuonSet; iSet++)
	if(iCharge == 0)
	  hDimuVtxProb[iSet][iCharge][iCat]->Write();
    }
  }

  hRecoNoBest->Write();


  for(int iSet = 0; iSet < kNbDimuonSet; iSet++){
//     for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
    for(int iCharge = 0; iCharge < 1; iCharge++){//only muons from OS pairs are filled 
      for(int iCat = 0; iCat < kNbMuCath; iCat++){

	hMuon_eta[iSet][iCharge][iCat]->Write();
	hMuon_pT[iSet][iCharge][iCat]->Write();
	hMuon_pT_eta[iSet][iCharge][iCat]->Write();

	if(iCat == 0)
	  hChi2GlobalFit[iSet][iCharge]->Write();
	if(iCat == 1){
	  hTM2DCompatibilityTight[iSet][iCharge]->Write();
	  hTMLastStationOptimizedLowPtLoose[iSet][iCharge]->Write();
	  hTrackerMuonArbitrated[iSet][iCharge]->Write();
	}
	if(iCat > 0)
	  hChi2TrackerFit[iSet][iCharge][iCat]->Write();

	hNHitsSilicon[iSet][iCharge][iCat]->Write();
	hCaloCompatibility[iSet][iCharge][iCat]->Write();
      }
    }
  }

  fOut->Close();
}
