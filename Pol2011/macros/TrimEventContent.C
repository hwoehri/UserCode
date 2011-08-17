#include "../interface/rootIncludes.inc"
#include "calcPol.C"

enum{L,R};
const Int_t kNbSpecies = 3;
enum{UPS1S, UPS2S, UPS3S};

void TrimEventContent(Int_t iRapBin = 1, 
		      Int_t iPTBin = 1,
		      Double_t fracL = 0.5, Double_t nSigma = 2., 
		      Int_t nUpsState=0//[0]... 1S, [1]... 2S, [2]... 3S
		      ){


  printf("\n\n\nfracL = %1.1f, nSigma = %1.1f, iState = %d, rap %d, pT %d\n", fracL, nSigma, nUpsState, iRapBin, iPTBin);
  Char_t name[100], title[100];  
  Char_t fileNameIn[100];
  sprintf(fileNameIn, "RootFiles/data_Ups_rap%d_pT%d.root", iRapBin, iPTBin);
  //==============================
  //read inputs from input file:
  TFile *fIn = new TFile(fileNameIn);
  TLorentzVector *lepP;
  TLorentzVector *lepN;
  TTree *treeIn = (TTree *) gDirectory->Get("selectedData");
  TH2D *hCosThetaPhiLR[onia::kNbFrames][2];
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    sprintf(name, "hCosThetaPhi_%s_L", onia::frameLabel[iFrame]);
    hCosThetaPhiLR[iFrame][L] = (TH2D *) gDirectory->Get(name);
    sprintf(name, "hCosThetaPhi_%s_R", onia::frameLabel[iFrame]);
    hCosThetaPhiLR[iFrame][R] = (TH2D *) gDirectory->Get(name);
  }
  //==============================

  //definition of output variables 
  Char_t fileNameOut[100];
  sprintf(fileNameOut, "RootFiles/data_%dSUps_rap%d_pT%d_%1.0fSigma.root", nUpsState+1, iRapBin, iPTBin, nSigma);
  TFile *fOut = new TFile(fileNameOut, "RECREATE");
  gStyle->SetPadRightMargin(0.2);
  TTree *treeOut = treeIn->CloneTree(0);
  // treeOut->SetName("data");
  TH2D *hCosThetaPhiBG[onia::kNbFrames];
  TH2D *hCosThetaPhiSignal[onia::kNbFrames];
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    //book the histo for the signal
    sprintf(name, "total_%s", onia::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", onia::frameLabel[iFrame], onia::frameLabel[iFrame]);
    hCosThetaPhiSignal[iFrame] = new TH2D(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, 
					  onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hCosThetaPhiSignal[iFrame]->Sumw2();
    //copy the L and R sideband histos into one output BG histogram
    hCosThetaPhiLR[iFrame][L]->Scale(fracL/hCosThetaPhiLR[iFrame][L]->Integral());
    hCosThetaPhiLR[iFrame][R]->Scale((1.-fracL)/hCosThetaPhiLR[iFrame][R]->Integral());
    sprintf(name, "background_%s", onia::frameLabel[iFrame]);
    hCosThetaPhiBG[iFrame] = (TH2D *) hCosThetaPhiLR[iFrame][L]->Clone(name);
    hCosThetaPhiBG[iFrame]->Add(hCosThetaPhiLR[iFrame][R]);
  }

  //======================================================
  //reading fit parameters to establish signal mass window
  //======================================================
  fIn->cd();
  TTree *treeFitPar = (TTree *) gDirectory->Get("massFitParameters");
  TF1 *fUps[kNbSpecies], *fBG = 0;
  fUps[0] = 0, fUps[1] = 0, fUps[2] = 0;
  treeFitPar->SetBranchAddress("fUps1S", &fUps[0]);
  treeFitPar->SetBranchAddress("fUps2S", &fUps[1]);
  treeFitPar->SetBranchAddress("fUps3S", &fUps[2]);
  treeFitPar->SetBranchAddress("fBG", &fBG);
  treeFitPar->LoadTree(0);
  treeFitPar->GetEntry(0);

  Double_t mass[kNbSpecies], sigma[kNbSpecies];
  for(int iState = 0; iState < kNbSpecies; iState++){
    mass[iState] = fUps[iState]->GetParameter(1);
    sigma[iState] = fUps[iState]->GetParameter(2);
  }
  printf("1S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS1S], sigma[UPS1S]);
  printf("2S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS2S], sigma[UPS2S]);
  printf("3S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS3S], sigma[UPS3S]);
  Double_t poleMass = mass[nUpsState], massMin, massMax;
  massMin = poleMass - nSigma*sigma[nUpsState];
  massMax = poleMass + nSigma*sigma[nUpsState];

  printf("--> signal mass window: %1.3f < M < %1.3f GeV\n", massMin, massMax);

  //calculate the fraction of BG under the signal mass window
  Double_t nBG = fBG->Integral(massMin, massMax);
  Double_t nSignal = fUps[nUpsState]->Integral(massMin, massMax);
  Double_t fracBG = nBG / (nBG + nSignal);
  sprintf(name, ";;fraction of BG in %1.1f sigma window", nSigma);
  TH1D *hFracBG = new TH1D("backgroundFraction", name, 1, 0., 1.);
  hFracBG->SetBinContent(1, fracBG);

  lepP = 0; lepN = 0;
  treeIn->SetBranchAddress("lepP", &lepP);
  treeIn->SetBranchAddress("lepN", &lepN);
  TLorentzVector *onia = new TLorentzVector();

  Double_t onia_mass;
  Int_t index;
  for(int iEn = 0; iEn < treeIn->GetEntries(); iEn++){ 
    Long64_t iEntry = treeIn->LoadTree(iEn);
    treeIn->GetEntry(iEntry);
    if(iEn % 10000 == 0)
      cout << "entry " << iEntry << " out of " << treeIn->GetEntries() << endl;

    *onia = *(lepP) + *(lepN);
    onia_mass = onia->M();
    if(onia_mass > massMin && onia_mass < massMax){
      treeOut->Fill(); //stores TLorenzVectors of the two muons

      //store now the cosTheta and phi distributions of the signal window:
      calcPol(*lepP, *lepN);

      for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
	hCosThetaPhiSignal[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
    }
  }



  //write the output
  fOut->cd();
  treeOut->Write();

  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    hCosThetaPhiBG[iFrame]->Write();
    hCosThetaPhiSignal[iFrame]->Write();
  }
  hFracBG->Write();

  fOut->Close();
  fIn->Close();
}
