#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

#include "TStyle.h"

// #include "PhysicsTools/RooStatsCms/interface/ClopperPearsonBinomialInterval.h"
// #include "PhysicsTools/RooStatsCms/interface/FeldmanCousinsBinomialInterval.h"

Double_t markerSize = 0.7;
Int_t const kNbSpecies = 4; //B0, Bp, Bs, LambdaB
Char_t *speciesName[kNbSpecies+1] = {"B0", "Bp", "Bs", "LambdaB", "NP"};
//the individual generations need to added with correct weights
//the downscaling factors (0.894 etc) refer to the fact that only 90% of
//the full MC sample was used to produce the acceptance maps
//to leave 10% for a MC "pseudoData" sample...
Double_t weights[kNbSpecies] = {10./(2.*44.46*0.894), 10./(2.*36.7*0.898), 10./(2.*33.5*0.9), 10./(2.*156.63*0.899)};
//generated histos
// TH1D *hGen_pt[jpsi::kNbRapForPTBins+1];
// TH1D *hGen_rap[jpsi::kNbPTBins+1];
// TH1D *hGen_phi[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
// TH1D *hGen_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbPolVar];//[all,pT1,pT2,pT...]
// TH1D *hGen_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1][jpsi::kNbPolVar];
// TH1D *hGen_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1][jpsi::kNbPolVar];
// TH2D *hGen2D_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1];
// TH2D *hGen2D_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1];
TH2D *hGen2D_pol_pT_rap[kNbSpecies+1][jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];

//reconstructed histos
// TH1D *Reco_pt[jpsi::kNbRapForPTBins+1];
// TH1D *Reco_rap[jpsi::kNbPTBins+1];
// TH1D *Reco_phi[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
// TH1D *Reco_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbPolVar];
// TH1D *Reco_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1][jpsi::kNbPolVar];
// TH1D *Reco_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1][jpsi::kNbPolVar];
// TH2D *Reco2D_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1];
// TH2D *Reco2D_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1];
TH2D *Reco2D_pol_pT_rap[kNbSpecies+1][jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];

//acceptance histos
TGraphAsymmErrors *gAcc_pt[jpsi::kNbRapForPTBins+1];
TGraphAsymmErrors *gAcc_rap[jpsi::kNbPTBins+1];
TGraphAsymmErrors *gAcc_phi[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
// TGraphAsymmErrors *gAcc_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbPolVar];
// TGraphAsymmErrors *gAcc_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1][jpsi::kNbPolVar];
// TGraphAsymmErrors *gAcc_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1][jpsi::kNbPolVar];
// TH2D *hAcc2D_pol_pT[jpsi::kNbFrames][jpsi::kNbPTBins+1];
// TH2D *hAcc2D_pol_rap[jpsi::kNbFrames][2*jpsi::kNbRapBins+1];
TH2D *hAcc2D_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
TH2D *hAccErr2D_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
// TGraphAsymmErrors *gAcc2D_pol_pT_rap[jpsi::kNbFrames][jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];

void ReadInHistos(Char_t *fileNameIn, Int_t iSpecies, Int_t rebinCosTh, Int_t rebinPhi);
void AddSpecies();
void PlotGenRec(Int_t iFrame, Bool_t skipAll);
void PlotAllGenRec(Int_t iFrame);
void PlotAllRec(Bool_t skipAll);
void CalcAcceptance();
void PlotAcceptance(Int_t iFrame);
void Plot2DAcceptance(Int_t iFrame, Char_t *hltTag);
void Plot2DAccOneByOne(Int_t iFrame, Char_t *hltTag, Int_t iRap, Int_t iPT);
void Plot2DAccErrOneByOne(Int_t iFrame, Char_t *hltTag, Int_t iRap, Int_t iPT);
void PlotKinVarAcceptance();
void WriteAccHistos(Char_t *fileNameOut);
void PlotAllAcceptances();
void Calc2DAcc(TH2D *hReco, TH2D *hGen, TGraphAsymmErrors *graph);

//================================================================
//usage: root getAcceptance.C+ or
//       'root getAcceptance.C+("HLT_Mu0Track0Jpsi", 2, 2)' (e.g.)
//================================================================
void getNPAcceptance(Char_t *hltTag = "HLT_Mu0Track0Jpsi", 
		     Char_t *lastLabel = "", //optional: if filename contains another string at the end
		     Int_t rebinCosTh = 2,
		     Int_t rebinPhi = 2){

  Char_t name[200];
  Char_t fileNameIn[4][200], fileNameOut[200];

  sprintf(fileNameIn[0],  "pol_MC_B0_all_%s_%s.root", hltTag, lastLabel);
  sprintf(fileNameIn[1],  "pol_MC_Bp_all_%s_%s.root", hltTag, lastLabel);
  sprintf(fileNameIn[2],  "pol_MC_Bs_all_%s_%s.root", hltTag, lastLabel);
  sprintf(fileNameIn[3],  "pol_MC_LambdaB_all_%s_%s.root", hltTag, lastLabel);
  sprintf(fileNameOut, "accHistos_NP_all_%s_%s.root", hltTag, lastLabel);

  for(int iSpecies = 0; iSpecies < kNbSpecies; iSpecies++)
    ReadInHistos(fileNameIn[iSpecies], iSpecies, rebinCosTh, rebinPhi);
  AddSpecies();

  // Bool_t skipAll = kTRUE; //if kTRUE --> do not show the "all pT" curve
//   PlotGenRec(CS, skipAll);
//   PlotGenRec(HX, skipAll);
//   PlotAllGenRec(CS);
//   PlotAllGenRec(HX);
//   PlotAllRec(skipAll);
  CalcAcceptance();

  // PlotKinVarAcceptance();
//   PlotAcceptance(CS);
// //   PlotAcceptance(HX);
// //   PlotAllAcceptances();
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++)
    Plot2DAcceptance(iFrame, hltTag);

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++)
    for(int iRap = 1; iRap <= jpsi::kNbRapForPTBins; iRap++)
      for(int iPT = 1; iPT <= jpsi::kNbPTBins; iPT++)
	Plot2DAccOneByOne(iFrame, hltTag, iRap, iPT);

  // for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++)
  //   for(int iRap = 1; iRap <= jpsi::kNbRapForPTBins; iRap++)
  //     for(int iPT = 1; iPT <= jpsi::kNbPTBins; iPT++){
  // 	Plot2DAccErrOneByOne(iFrame, hltTag, iRap, iPT);

  WriteAccHistos(fileNameOut);
}

//===============================
void PlotGenRec(Int_t iFrame, Bool_t skipAll){

//   Int_t pTMin = 0, pTMax = jpsi::kNbPTBins; //jpsi::kNbPTBins: last one will be suppressed
//   Int_t indexMax = 0;
//   if(skipAll){
//     pTMin = 1;
//     indexMax = 2;
//   }
  
//   Char_t name[100];
//   //=======================================
//   //1.) generated cosTheta for different pT bins
//   //=======================================
//   sprintf(name, "c1CosThGen_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c1CosThGen = new TCanvas(name, "generated cosTheta for PT Bins", 700, 500);
//   TH1F *hFrame1a = gPad->DrawFrame(-1, 0., 1., 1.7*hGen_pol_pT[iFrame][indexMax][cosThPol]->GetMaximum());
//   sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);
//   hFrame1a->SetXTitle(name);
//   hFrame1a->SetYTitle("dN/d(cos#theta)_{generated}");
//   for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
//     hGen_pol_pT[iFrame][iPTBin][cosThPol]->Draw("histsame");    
//   }
//   TLegend *leg1a = new TLegend(0.65,0.735119,0.9866071,0.985119);
//   for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
//     if(iPTBin == 0) sprintf(name, "all p_{T}");
//     else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
//     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
//     leg1a->AddEntry(hGen_pol_pT[iFrame][iPTBin][cosThPol], name, "pl");
//   }
//   leg1a->SetTextSize(0.035); leg1a->SetFillColor(0);
//   leg1a->Draw();

//   sprintf(name, "Figures/generated_cosTheta_%s_pTBins.eps", jpsi::frameLabel[iFrame]);  c1CosThGen->Print(name);
//   sprintf(name, "Figures/generated_cosTheta_%s_pTBins.pdf", jpsi::frameLabel[iFrame]);  c1CosThGen->Print(name);
//   sprintf(name, "Figures/generated_cosTheta_%s_pTBins.gif", jpsi::frameLabel[iFrame]);  c1CosThGen->Print(name);

//   //=======================================
//   //2.) generated phi for different pT bins
//   //=======================================
//   sprintf(name, "c1PhiGen_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c1PhiGen = new TCanvas(name, "generated phi for PT Bins", 700, 500);
//   TH1F *hFrame1b = gPad->DrawFrame(0., 0., 360., 1.2*hGen_pol_pT[iFrame][indexMax][phiPol]->GetMaximum());
//   sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);
//   hFrame1b->SetXTitle(name);
//   hFrame1b->SetYTitle("dN/d#phi_{generated}");
//   for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
//     hGen_pol_pT[iFrame][iPTBin][phiPol]->Draw("histsame");    
//   }
// //   TLegend *leg1b = new TLegend(0.65,0.735119,0.9866071,0.985119);
// //   for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
// //     if(iPTBin == 0) sprintf(name, "all p_{T}");
// //     else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
// //     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
// //     leg1b->AddEntry(hGen_pol_pT[iFrame][iPTBin][phiPol], name, "pl");
// //   }
// //   leg1a->SetTextSize(0.035); leg1a->SetFillColor(0);
// //   leg1a->Draw();

//   sprintf(name, "Figures/generated_phi_%s_pTBins.eps", jpsi::frameLabel[iFrame]);  c1PhiGen->Print(name);
//   sprintf(name, "Figures/generated_phi_%s_pTBins.pdf", jpsi::frameLabel[iFrame]);  c1PhiGen->Print(name);
//   sprintf(name, "Figures/generated_phi_%s_pTBins.gif", jpsi::frameLabel[iFrame]);  c1PhiGen->Print(name);

//   //=======================================
//   //3.) reconstructed cosTheta for different pT bins
//   //=======================================
//   sprintf(name, "c1CosThRec_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c1CosThRec = new TCanvas(name, "reconstructed cosTheta for PT Bins", 900, 700);
//   TH1F *hFrame3a = gPad->DrawFrame(-1, 0., 1, 1.2*Reco_pol_pT[iFrame][indexMax][cosThPol]->GetMaximum());
//   sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);
//   hFrame3a->SetXTitle(name);
//   hFrame3a->SetYTitle("dN/d(cos#theta)_{reconstructed}");
//   for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
//     Reco_pol_pT[iFrame][iPTBin][cosThPol]->Draw("psame");    
//   }
//   TLegend *leg3a = new TLegend(0.65,0.735119,0.9866071,0.985119);
//   for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//     if(iPTBin == 0) sprintf(name, "all p_{T}");
//     else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
//     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
//     leg3a->AddEntry(Reco_pol_pT[iFrame][iPTBin][cosThPol], name, "pl");
//   }
//   leg3a->SetTextSize(0.035); leg3a->SetFillColor(0);
//   leg3a->Draw();

//   sprintf(name, "Figures/reconstructed_cosTheta_%s_pTBins.eps", jpsi::frameLabel[iFrame]);  c1CosThRec->Print(name);
//   sprintf(name, "Figures/reconstructed_cosTheta_%s_pTBins.pdf", jpsi::frameLabel[iFrame]);  c1CosThRec->Print(name);
//   sprintf(name, "Figures/reconstructed_cosTheta_%s_pTBins.gif", jpsi::frameLabel[iFrame]);  c1CosThRec->Print(name);

//   //=======================================
//   //4.) reconstructed phi for different pT bins
//   //=======================================
//   sprintf(name, "c1PhiRec_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c1PhiRec = new TCanvas(name, "reconstructed phi for PT Bins", 900, 700);
//   TH1F *hFrame3b = gPad->DrawFrame(0., 0., 360., 1.2*Reco_pol_pT[iFrame][indexMax][phiPol]->GetMaximum());
//   sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);
//   hFrame3b->SetXTitle(name);
//   hFrame3b->SetYTitle("dN/d#phi_{reconstructed}");
//   for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
//     Reco_pol_pT[iFrame][iPTBin][phiPol]->Draw("psame");    
//   }
// //   TLegend *leg3b = new TLegend(0.65,0.735119,0.9866071,0.985119);
// //   for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
// //     if(iPTBin == 0) sprintf(name, "all p_{T}");
// //     else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
// //     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
// //     leg3b->AddEntry(Reco_pol_pT[iFrame][iPTBin][phiPol], name, "pl");
// //   }
// //   leg1a->SetTextSize(0.035); leg1a->SetFillColor(0);
// //   leg1a->Draw();

//   sprintf(name, "Figures/reconstructed_phi_%s_pTBins.eps", jpsi::frameLabel[iFrame]);  c1PhiRec->Print(name);
//   sprintf(name, "Figures/reconstructed_phi_%s_pTBins.pdf", jpsi::frameLabel[iFrame]);  c1PhiRec->Print(name);
//   sprintf(name, "Figures/reconstructed_phi_%s_pTBins.gif", jpsi::frameLabel[iFrame]);  c1PhiRec->Print(name);

}

//===============================
void PlotAllRec(Bool_t skipAll){

  // Int_t pTMin = 0, pTMax = jpsi::kNbPTBins; //jpsi::kNbPTBins: last one will be suppressed
  // Int_t indexMax = 0;
  // if(skipAll){
  //   pTMin = 1;
  //   indexMax = 2;
  // }

  // Char_t name[100];
  // sprintf(name, "c1Recs");
  // TCanvas *c1CosThRec = new TCanvas(name, "reconstructed CS and HX for PT Bins", 1000, 700);
  // c1CosThRec->Divide(2,2);
  // TH1F *hFrame3a[4];
  // //HX: cosTheta
  // c1CosThRec->cd(1);
  // hFrame3a[0] = gPad->DrawFrame(-1, 0., 1, 1.2*Reco_pol_pT[jpsi::HX][indexMax][cosThPol]->GetMaximum());
  // sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::HX]);
  // hFrame3a[0]->SetXTitle(name);
  // hFrame3a[0]->SetYTitle("dN/d(cos#theta)_{reconstructed}");
  // for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
  //   Reco_pol_pT[jpsi::HX][iPTBin][cosThPol]->Draw("psame");    
  // }
  // TLegend *leg3a = new TLegend(0.65,0.735119,0.9866071,0.985119);
  // for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
  //   if(iPTBin == 0) sprintf(name, "all p_{T}");
  //   else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
  //   else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
  //   leg3a->AddEntry(Reco_pol_pT[jpsi::HX][iPTBin][cosThPol], name, "pl");
  // }
  // leg3a->SetTextSize(0.035); leg3a->SetFillColor(0);
  // leg3a->Draw();

  // //HX: phi
  // c1CosThRec->cd(2);
  // hFrame3a[1] = gPad->DrawFrame(0., 0., 360., 1.2*Reco_pol_pT[jpsi::HX][indexMax][phiPol]->GetMaximum());
  // sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::HX]);
  // hFrame3a[1]->SetXTitle(name);
  // hFrame3a[1]->SetYTitle("dN/d#phi_{reconstructed}");
  // for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
  //   Reco_pol_pT[jpsi::HX][iPTBin][phiPol]->Draw("psame");    
  // }
  // //CS: cosTheta
  // c1CosThRec->cd(3);
  // hFrame3a[2] = gPad->DrawFrame(-1, 0., 1, 1.2*Reco_pol_pT[jpsi::CS][indexMax][cosThPol]->GetMaximum());
  // sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::CS]);
  // hFrame3a[2]->SetXTitle(name);
  // hFrame3a[2]->SetYTitle("dN/d(cos#theta)_{reconstructed}");
  // for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
  //   Reco_pol_pT[jpsi::CS][iPTBin][cosThPol]->Draw("psame");    
  // }
  // //CS: phi
  // c1CosThRec->cd(4);
  // hFrame3a[3] = gPad->DrawFrame(0., 0., 360., 1.2*Reco_pol_pT[jpsi::CS][indexMax][phiPol]->GetMaximum());
  // sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[jpsi::CS]);
  // hFrame3a[3]->SetXTitle(name);
  // hFrame3a[3]->SetYTitle("dN/d#phi_{reconstructed}");
  // for(int iPTBin = pTMin; iPTBin < pTMax; iPTBin++){
  //   Reco_pol_pT[jpsi::CS][iPTBin][phiPol]->Draw("psame");    
  // }

  // sprintf(name, "Figures/rec_pTBins.gif"); c1CosThRec->Print(name);
  // sprintf(name, "Figures/rec_pTBins.eps"); c1CosThRec->Print(name);
  // sprintf(name, "Figures/rec_pTBins.pdf"); c1CosThRec->Print(name);

}


//===============================
void PlotAcceptance(Int_t iFrame){

//   Char_t name[100];
//   //=============================================
//   //1.) cosTheta acceptance for different pT bins
//   //=============================================
//   sprintf(name, "c1CosTh_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c1CosTh = new TCanvas(name, "acc vs cosTheta for PT Bins", 900, 700);
// //   TH1F *hFrame1a = gPad->DrawFrame(-1, 0., 1, 1.2*gAcc_pol_pT[iFrame][jpsi::kNbPTBins][cosThPol]->GetMaximum());
//   TH1F *hFrame1a = gPad->DrawFrame(-1, 0., 1, 1.0);
//   sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);
//   hFrame1a->SetXTitle(name);
//   hFrame1a->SetYTitle("Acceptance");
//   for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//      gAcc_pol_pT[iFrame][iPTBin][cosThPol]->Draw("psame");    
//   }
//   TLegend *leg1a = new TLegend(0.65,0.735119,0.9866071,0.985119);
//   for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//     if(iPTBin == 0) sprintf(name, "all p_{T}");
//     else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
//     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
//     leg1a->AddEntry(gAcc_pol_pT[iFrame][iPTBin][cosThPol], name, "pl");
//   }
//   leg1a->SetTextSize(0.035); leg1a->SetFillColor(0);
//   leg1a->Draw();

//   sprintf(name, "Figures/acceptance_cosTheta_%s_pTBins.eps", jpsi::frameLabel[iFrame]);  c1CosTh->Print(name);
//   sprintf(name, "Figures/acceptance_cosTheta_%s_pTBins.pdf", jpsi::frameLabel[iFrame]);  c1CosTh->Print(name);
//   sprintf(name, "Figures/acceptance_cosTheta_%s_pTBins.gif", jpsi::frameLabel[iFrame]);  c1CosTh->Print(name);

//   //=======================================
//   //2.) phi acceptance for different pT bins
//   //=======================================
//   sprintf(name, "c1Phi_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c1Phi = new TCanvas(name, "acc vs phi for PT Bins", 900, 700);
// //   TH1F *hFrame1b = gPad->DrawFrame(0., 0., 360., 1.2*gAcc_pol_pT[iFrame][0][phiPol]->GetMaximum());
//   TH1F *hFrame1b = gPad->DrawFrame(0., 0., 360., 1.0);
//   sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);
//   hFrame1b->SetXTitle(name);
//   hFrame1b->SetYTitle("Acceptance");
//   for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//     gAcc_pol_pT[iFrame][iPTBin][phiPol]->Draw("psame");    
//   }
//   TLegend *leg1b = new TLegend(0.65,0.735119,0.9866071,0.985119);
//   for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//     if(iPTBin == 0) sprintf(name, "all p_{T}");
//     else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
//     else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
//     leg1b->AddEntry(gAcc_pol_pT[iFrame][iPTBin][phiPol], name, "pl");
//   }
//   leg1a->SetTextSize(0.035); leg1a->SetFillColor(0);
//   leg1a->Draw();

//   sprintf(name, "Figures/acceptance_phi_%s_pTBins.eps", jpsi::frameLabel[iFrame]);  c1Phi->Print(name);
//   sprintf(name, "Figures/acceptance_phi_%s_pTBins.pdf", jpsi::frameLabel[iFrame]);  c1Phi->Print(name);
//   sprintf(name, "Figures/acceptance_phi_%s_pTBins.gif", jpsi::frameLabel[iFrame]);  c1Phi->Print(name);

//   //=======================================
//   //3.) cosTheta acceptance for different rap bins
//   //=======================================
//   sprintf(name, "c2CosTh_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c2CosTh = new TCanvas(name, "acc vs cosTheta for rap Bins", 900, 700);
// //   TH1F *hFrame2a = gPad->DrawFrame(-1, 0., 1, 1.2*gAcc_pol_rap[iFrame][0][cosThPol]->GetMaximum());
//   TH1F *hFrame2a = gPad->DrawFrame(-1, 0., 1, 1.0);
//   sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);
//   hFrame2a->SetXTitle(name);
//   hFrame2a->SetYTitle("Acceptance");
//   for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
//     gAcc_pol_rap[iFrame][iRapBin][cosThPol]->Draw("psame");    
//   }
//   TLegend *leg2a = new TLegend(0.75,0.65,0.9866071,0.985119);
//   for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
//     if(iRapBin == 0) sprintf(name, "all y");
//     else sprintf(name, "%1.1f < y < %1.1f", rapRange[iRapBin-1], rapRange[iRapBin]);
//     leg2a->AddEntry(gAcc_pol_rap[iFrame][iRapBin][cosThPol], name, "pl");
//   }
//   leg2a->SetTextSize(0.035); leg2a->SetFillColor(0);
//   leg2a->Draw();

//   sprintf(name, "Figures/acceptance_cosTheta_%s_rapBins.eps", jpsi::frameLabel[iFrame]);  c2CosTh->Print(name);
//   sprintf(name, "Figures/acceptance_cosTheta_%s_rapBins.pdf", jpsi::frameLabel[iFrame]);  c2CosTh->Print(name);
//   sprintf(name, "Figures/acceptance_cosTheta_%s_rapBins.gif", jpsi::frameLabel[iFrame]);  c2CosTh->Print(name);

//   //=======================================
//   //4.) phi acceptance for different rap bins
//   //=======================================
//   sprintf(name, "c2Phi_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c2Phi = new TCanvas(name, "acc vs phi for rap Bins", 900, 700);
// //   TH1F *hFrame2b = gPad->DrawFrame(0., 0., 360., 1.2*gAcc_pol_rap[iFrame][0][phiPol]->GetMaximum());
//   TH1F *hFrame2b = gPad->DrawFrame(0., 0., 360., 1.0);
//   sprintf(name, "#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);
//   hFrame2b->SetXTitle(name);
//   hFrame2b->SetYTitle("Acceptance");
//   for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
//     gAcc_pol_rap[iFrame][iRapBin][phiPol]->Draw("psame");    
//   }
//   TLegend *leg2b = new TLegend(0.75,0.735119,0.9866071,0.985119);
//   for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
//     if(iRapBin == 0) sprintf(name, "all y");
//     else sprintf(name, "%1.1f < y < %1.1f", rapRange[iRapBin-1], rapRange[iRapBin]);
//     leg2b->AddEntry(gAcc_pol_rap[iFrame][iRapBin][phiPol], name, "pl");
//   }
//   leg2a->SetTextSize(0.035); leg2a->SetFillColor(0);
//   leg2a->Draw();

//   sprintf(name, "Figures/acceptance_phi_%s_rapBins.eps", jpsi::frameLabel[iFrame]);  c2Phi->Print(name);
//   sprintf(name, "Figures/acceptance_phi_%s_rapBins.pdf", jpsi::frameLabel[iFrame]);  c2Phi->Print(name);
//   sprintf(name, "Figures/acceptance_phi_%s_rapBins.gif", jpsi::frameLabel[iFrame]);  c2Phi->Print(name);

//   //=====================================================
//   //5.) cosTheta acceptance for different rap and pT bins
//   //=====================================================
//   sprintf(name, "c3CosTh_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c3CosTh = new TCanvas(name, "acc vs cosTheta for pT and rap Bins", 900, 700);
// //   TH1F *hFrame3a = gPad->DrawFrame(-1, 0., 1, 1.2*gAcc_pol_rap[iFrame][0][cosThPol]->GetMaximum());
//   c3CosTh->Divide(1,jpsi::kNbRapForPTBins);
//   TH1F *hFrame3a[jpsi::kNbRapForPTBins];
//   for(int iRap = 1; iRap < jpsi::kNbRapForPTBins+1; iRap++){

//     c3CosTh->cd(iRap);
//     hFrame3a[iRap] = gPad->DrawFrame(-1, 0., 1, 1.0);
// //     hFrame3a[iRap] = gPad->DrawFrame(-1, 0., 1, 1.2*gAcc_pol_pT_rap[iFrame][1][iRap][cosThPol]->GetMaximum());
//     sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);
//     hFrame3a[iRap]->SetXTitle(name);
//     hFrame3a[iRap]->SetYTitle("Acceptance");
//     for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//        gAcc_pol_pT_rap[iFrame][iPTBin][iRap][cosThPol]->Draw("psame");
//     }
//     sprintf(name, "%1.2f < y < %1.2f", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);
//     TLegend *leg3a = new TLegend(0.88,0.55,0.9866071,0.985119, name);
//     for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//       if(iPTBin == 0) sprintf(name, "all p_{T}");
//       else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
//       else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
//       leg3a->AddEntry(gAcc_pol_pT_rap[iFrame][iPTBin][iRap][phiPol], name, "pl");
//     }
//     leg3a->SetTextSize(0.035); leg3a->SetFillColor(0);  leg3a->Draw();
//   }
//   sprintf(name, "Figures/acceptance_cosTheta_%s_pTBins_rap.eps", jpsi::frameLabel[iFrame]);  c3CosTh->Print(name);
//   sprintf(name, "Figures/acceptance_cosTheta_%s_pTBins_rap.pdf", jpsi::frameLabel[iFrame]);  c3CosTh->Print(name);
//   sprintf(name, "Figures/acceptance_cosTheta_%s_pTBins_rap.gif", jpsi::frameLabel[iFrame]);  c3CosTh->Print(name);

//   //=====================================================
//   //6.) phi acceptance for different rap and pT bins
//   //=====================================================
//   sprintf(name, "c3Phi_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c3Phi = new TCanvas(name, "acc vs phi for pT and rap Bins", 900, 700);
//   c3Phi->Divide(1,jpsi::kNbRapForPTBins);
//   TH1F *hFrame3b[jpsi::kNbRapForPTBins];
//   for(int iRap = 1; iRap < jpsi::kNbRapForPTBins+1; iRap++){

//     c3Phi->cd(iRap);
//     hFrame3b[iRap] = gPad->DrawFrame(0., 0., 360., 1.0);
// //      hFrame3b[iRap] = gPad->DrawFrame(0., 0., 360., 1.2*gAcc_pol_pT_rap[iFrame][1][iRap][phiPol]->GetMaximum() );
//     sprintf(name, "#phi_{%s}", jpsi::frameLabel[iFrame]);
//     hFrame3b[iRap]->SetXTitle(name);
//     hFrame3b[iRap]->SetYTitle("Acceptance");
//     for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//       gAcc_pol_pT_rap[iFrame][iPTBin][iRap][phiPol]->Draw("psame");    
//     }
//     sprintf(name, "%1.2f < y < %1.2f", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap]);
//     TLegend *leg3b = new TLegend(0.88,0.65,0.9866071,0.985119, name);
//     for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
//       if(iPTBin == 0) sprintf(name, "all p_{T}");
//       else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
//       else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
//       leg3b->AddEntry(gAcc_pol_pT_rap[iFrame][iPTBin][iRap][phiPol], name, "pl");
//     }
//     leg3b->SetTextSize(0.035); leg3b->SetFillColor(0);  leg3b->Draw();
//   }
//   sprintf(name, "Figures/acceptance_phi_%s_pTBins_rap.eps", jpsi::frameLabel[iFrame]);  c3Phi->Print(name);
//   sprintf(name, "Figures/acceptance_phi_%s_pTBins_rap.pdf", jpsi::frameLabel[iFrame]);  c3Phi->Print(name);
//   sprintf(name, "Figures/acceptance_phi_%s_pTBins_rap.gif", jpsi::frameLabel[iFrame]);  c3Phi->Print(name);

}

//===============================
void PlotAllAcceptances(){

  // Char_t name[100];
  // //===========================================================
  // //cosTheta and phi acceptance for CS and HX for pT dependence
  // //===========================================================
  // sprintf(name, "c1All_PT");
  // TCanvas *c1CosTh = new TCanvas(name, "acc vs cosTheta for PT Bins", 1000, 700);
  // c1CosTh->Divide(2,2);
  // TH1F *hFrame1a[4];
  // for(int i = 0; i < 4; i++){
  //   c1CosTh->cd(i+1);
  //   if(i == 0 || i == 2) hFrame1a[i] = gPad->DrawFrame(-1, 0., 1, 1.0);
  //   else if(i == 1 || i == 3) hFrame1a[i] = gPad->DrawFrame(0., 0., 360., 1.0);
  //   if(i == 0) sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::CS]);
  //   else if(i == 1) sprintf(name, "#phi_{%s}", jpsi::frameLabel[jpsi::CS]);
  //   else if(i == 2) sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::HX]);
  //   else if(i == 3) sprintf(name, "#phi_{%s}", jpsi::frameLabel[jpsi::HX]);
  //   hFrame1a[i]->SetXTitle(name);
  //   hFrame1a[i]->SetYTitle("Acceptance");
  //   for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  //     if(i == 0) gAcc_pol_pT[jpsi::CS][iPTBin][cosThPol]->Draw("psame");    
  //     else if(i == 1) gAcc_pol_pT[jpsi::CS][iPTBin][phiPol]->Draw("psame");    
  //     else if(i == 2) gAcc_pol_pT[jpsi::HX][iPTBin][cosThPol]->Draw("psame");    
  //     else if(i == 3) gAcc_pol_pT[jpsi::HX][iPTBin][phiPol]->Draw("psame");    
  //   }
  //   if(i == 3){
  //     TLegend *leg1a = new TLegend(0.65,0.735119,0.9866071,0.985119);
  //     for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  // 	if(iPTBin == 0) sprintf(name, "all p_{T}");
  // 	else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
  // 	else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
  // 	leg1a->AddEntry(gAcc_pol_pT[jpsi::HX][iPTBin][phiPol], name, "pl");
  //     }
  //     leg1a->SetTextSize(0.035); leg1a->SetFillColor(0);
  //     leg1a->Draw();
  //   }
  // }
  // sprintf(name, "Figures/acceptances_pTBins.eps");  c1CosTh->Print(name);
  // sprintf(name, "Figures/acceptances_pTBins.pdf");  c1CosTh->Print(name);
  // sprintf(name, "Figures/acceptances_pTBins.gif");  c1CosTh->Print(name);

  // //===========================================================
  // //cosTheta and phi acceptance for CS and HX for rap dependence
  // //===========================================================
  // sprintf(name, "c2All_Rap");
  // TCanvas *c2CosTh = new TCanvas(name, "acc vs cosTheta for Rap Bins", 1000, 700);
  // c2CosTh->Divide(2,2);
  // TH1F *hFrame2a[4];
  // for(int i = 0; i < 4; i++){
  //   c2CosTh->cd(i+1);
  //   if(i == 0 || i == 2) hFrame2a[i] = gPad->DrawFrame(-1, 0., 1, 1.0);
  //   else if(i == 1 || i == 3) hFrame2a[i] = gPad->DrawFrame(0., 0., 360., 1.0);
  //   if(i == 0) sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::CS]);
  //   else if(i == 1) sprintf(name, "#phi_{%s}", jpsi::frameLabel[jpsi::CS]);
  //   else if(i == 2) sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[jpsi::HX]);
  //   else if(i == 3) sprintf(name, "#phi_{%s}", jpsi::frameLabel[jpsi::HX]);
  //   hFrame2a[i]->SetXTitle(name);
  //   hFrame2a[i]->SetYTitle("Acceptance");
  //   for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
  //     if(i == 0) gAcc_pol_rap[jpsi::CS][iRapBin][cosThPol]->Draw("psame");    
  //     else if(i == 1) gAcc_pol_rap[jpsi::CS][iRapBin][phiPol]->Draw("psame");    
  //     else if(i == 2) gAcc_pol_rap[jpsi::HX][iRapBin][cosThPol]->Draw("psame");    
  //     else if(i == 3) gAcc_pol_rap[jpsi::HX][iRapBin][phiPol]->Draw("psame");    
  //   }
  //   if(i == 3){
  //     TLegend *leg2a = new TLegend(0.75,0.65,0.9866071,0.985119);
  //     for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
  // 	if(iRapBin == 0) sprintf(name, "all y");
  // 	else sprintf(name, "%1.1f < y < %1.1f", rapRange[iRapBin-1], rapRange[iRapBin]);
  // 	leg2a->AddEntry(gAcc_pol_rap[jpsi::HX][iRapBin][phiPol], name, "pl");
  //     }
  //     leg2a->SetTextSize(0.035); leg2a->SetFillColor(0);
  //     leg2a->Draw();
  //   }
  // }
  // sprintf(name, "Figures/acceptances_rapBins.eps");  c2CosTh->Print(name);
  // sprintf(name, "Figures/acceptances_rapBins.pdf");  c2CosTh->Print(name);
  // sprintf(name, "Figures/acceptances_rapBins.gif");  c2CosTh->Print(name);
}

//===============================
void Plot2DAcceptance(Int_t iFrame, Char_t *hltTag){

  gStyle->SetOptStat(0);
  gStyle->SetTitleW(1);

  Char_t name[100];
  Char_t title[100];
//   //=============================================
//   //1.) 2D acceptance for different pT bins
//   //=============================================
//   sprintf(name, "c1_2D_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c1_2D = new TCanvas(name, "acc vs cosTheta and phi for PT Bins", 1000, 700);
//   c1_2D->Divide(3,2);
//   TH2F *hFrame1a[jpsi::kNbPTBins+1];
//   TLatex *tex1[jpsi::kNbPTBins+1];
//   for(int iPT = 1; iPT <= jpsi::kNbPTBins; iPT++){

//     c1_2D->cd(iPT);

//     if(iPT == 0) sprintf(name, "J/#psi: all p_{T}");
//     else if(iPT == jpsi::kNbPTBins) sprintf(name, "J/#psi: p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPT-1]);
//     else sprintf(name, "J/#psi: %1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPT-1], jpsi::pTRange[iPT]);
//     hAcc2D_pol_pT[iFrame][iPT]->SetTitle(name);
//     hAcc2D_pol_pT[iFrame][iPT]->Draw("colz");//"colz" or "cont"
// //     hAcc2D_pol_pT[iFrame][iPT]->Draw("lego2");//"colz" or "cont"
//   }
  
//   sprintf(name, "Figures/acceptance2D_%s_%s_pTBins.eps", jpsi::frameLabel[iFrame], hltTag);  c1_2D->Print(name);
//   sprintf(name, "Figures/acceptance2D_%s_%s_pTBins.pdf", jpsi::frameLabel[iFrame], hltTag);  c1_2D->Print(name);
//   sprintf(name, "Figures/acceptance2D_%s_%s_pTBins.gif", jpsi::frameLabel[iFrame], hltTag);  c1_2D->Print(name);

//   //========================================
//   //2.) 2D acceptance for different rap bins
//   //========================================
//   sprintf(name, "c2_2D_%s", jpsi::frameLabel[iFrame]);
//   TCanvas *c2_2D = new TCanvas(name, "acc vs cosTheta and phi for rap Bins", 1000, 700);
//   c2_2D->Divide(3,2);
//   TH2F *hFrame2a[jpsi::kNbRapBins+1]; 
//   for(int iRapBin = 1; iRapBin <= jpsi::kNbRapBins; iRapBin++){

//     c2_2D->cd(iRapBin);
//     if(iRapBin == 0) sprintf(name, "all y");
//     else sprintf(name, "J/#psi: %1.1f < y < %1.1f GeV/c", rapRange[iRapBin-1], rapRange[iRapBin]);
//     hAcc2D_pol_rap[iFrame][iRapBin]->SetTitle(name);
//     hAcc2D_pol_rap[iFrame][iRapBin]->Draw("colz"); //"colz" or "cont"
// //     hAcc2D_pol_rap[iFrame][iRapBin]->Draw("lego2"); //"colz" or "cont"
//   }

//   sprintf(name, "Figures/acceptance2D_%s_%s_rapBins.eps", jpsi::frameLabel[iFrame], hltTag);  c2_2D->Print(name);
//   sprintf(name, "Figures/acceptance2D_%s_%s_rapBins.pdf", jpsi::frameLabel[iFrame], hltTag);  c2_2D->Print(name);
//   sprintf(name, "Figures/acceptance2D_%s_%s_rapBins.gif", jpsi::frameLabel[iFrame], hltTag);  c2_2D->Print(name);

  //=====================================================
  //3.) cosTheta acceptance for different rap and pT bins
  //=====================================================
  TCanvas *c3_2D[jpsi::kNbRapForPTBins];
  for(int iRap = 1; iRap <= jpsi::kNbRapForPTBins; iRap++){

    sprintf(name, "c3_2D_%s_rap%d", jpsi::frameLabel[iFrame], iRap);
    sprintf(title, "acc vs cosTheta and phi for pT; rap Bin %d", iRap);
    c3_2D[iRap] = new TCanvas(name, title, 1000, 700);
    c3_2D[iRap]->Divide(3,2);

    for(int iPT = 1; iPT <= jpsi::kNbPTBins; iPT++){

      c3_2D[iRap]->cd(iPT);
      if(iPT == 0) 
	sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T} (%s)", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], hltTag);
      else if(iPT == jpsi::kNbPTBins) 
	sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c (%s)\n", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPT-1], hltTag);
      else 
	sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c (%s)", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPT-1], jpsi::pTRange[iPT], hltTag);
      hAcc2D_pol_pT_rap[iFrame][iPT][iRap]->SetTitle(name);
      hAcc2D_pol_pT_rap[iFrame][iPT][iRap]->Draw("colz");//"colz" or "cont"
//       hAcc2D_pol_pT_rap[iFrame][iPT][iRap]->Draw("lego2");//"colz" or "cont"
    }

    sprintf(name, "Figures/acceptance2D_%s_%s_%s_pTBins_rap%d.eps", speciesName[kNbSpecies], jpsi::frameLabel[iFrame], hltTag, iRap); c3_2D[iRap]->Print(name);
    sprintf(name, "Figures/acceptance2D_%s_%s_%s_pTBins_rap%d.pdf", speciesName[kNbSpecies], jpsi::frameLabel[iFrame], hltTag, iRap); c3_2D[iRap]->Print(name);
    sprintf(name, "Figures/acceptance2D_%s_%s_%s_pTBins_rap%d.gif", speciesName[kNbSpecies], jpsi::frameLabel[iFrame], hltTag, iRap); c3_2D[iRap]->Print(name);
  }
}

//===============================
void Plot2DAccOneByOne(Int_t iFrame, Char_t *hltTag, Int_t iRap, Int_t iPT){

  gStyle->SetOptStat(0);
  gStyle->SetTitleW(1);

  // hAcc2D_pol_pT_rap[iFrame][iPT][iRap]->SetMaximum(1.); //set the scale to 100 %

  Char_t name[100];
  Char_t title[100];

  //=========================================================
  //phi vs. cosTheta acceptance for different rap and pT bins
  //=========================================================
  TCanvas *c30_2D;
  sprintf(name, "c30_2D_%s_rap%d_pT%d", jpsi::frameLabel[iFrame], iRap, iPT);
  sprintf(title, "acc vs cosTheta and phi for pT %d rap Bin %d", iPT, iRap);
  c30_2D = new TCanvas(name, title, 500, 500);
  // TH1F *hFrame30 = gPad->DrawFrame(-1., 0., 1., 360.);
  // hFrame30->SetXTitle("hAcc2D_pol_pT_rap[iFrame][iPT][iRap]->GetXTitle());
  // hFrame30->SetYTitle("hAcc2D_pol_pT_rap[iFrame][iPT][iRap]->GetYTitle());

  if(iPT == 0) 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T} (%s)", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], hltTag);
  else if(iPT == jpsi::kNbPTBins) 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c (%s)\n", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPT-1], hltTag);
  else 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c (%s)", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPT-1], jpsi::pTRange[iPT], hltTag);
  hAcc2D_pol_pT_rap[iFrame][iPT][iRap]->SetTitle(name);
  hAcc2D_pol_pT_rap[iFrame][iPT][iRap]->Draw("colz");//"colz" or "cont"
  //       hAcc2D_pol_pT_rap[iFrame][iPT][iRap]->Draw("lego2");//"colz" or "cont"

  sprintf(name, "Figures/acceptance2D_%s_%s_%s_rap%d_pT%d.eps", speciesName[kNbSpecies], jpsi::frameLabel[iFrame], hltTag, iRap, iPT); 
  c30_2D->Print(name);
  sprintf(name, "Figures/acceptance2D_%s_%s_%s_rap%d_pT%d.pdf", speciesName[kNbSpecies], jpsi::frameLabel[iFrame], hltTag, iRap, iPT); 
  c30_2D->Print(name);
  sprintf(name, "Figures/acceptance2D_%s_%s_%s_rap%d_pT%d.gif", speciesName[kNbSpecies], jpsi::frameLabel[iFrame], hltTag, iRap, iPT); 
  c30_2D->Print(name);
}

//===============================
void Plot2DAccErrOneByOne(Int_t iFrame, Char_t *hltTag, Int_t iRap, Int_t iPT){

  gStyle->SetTitleW(1);
  gStyle->SetOptStat(0);

  // hAcc2D_pol_pT_rap[iFrame][iPT][iRap]->SetMaximum(1.); //set the scale to 100 %

  Char_t name[100];
  Char_t title[100];

  //==================================================================
  //phi vs. cosTheta error on acceptance for different rap and pT bins
  //==================================================================
  TCanvas *c31_2D;
  sprintf(name, "c31_2D_%s_rap%d_pT%d", jpsi::frameLabel[iFrame], iRap, iPT);
  sprintf(title, "acc error vs cosTheta and phi for pT %d rap Bin %d", iPT, iRap);
  c31_2D = new TCanvas(name, title, 500, 500);
  // TH1F *hFrame31 = gPad->DrawFrame(-1., 0., 1., 360.);
  // hFrame31->SetXTitle("hAccErr2D_pol_pT_rap[iFrame][iPT][iRap]->GetXTitle());
  // hFrame31->SetYTitle("hAccErr2D_pol_pT_rap[iFrame][iPT][iRap]->GetYTitle());

  if(iPT == 0) 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, all p_{T} (%s)", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], hltTag);
  else if(iPT == jpsi::kNbPTBins) 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, p_{T} > %1.1f GeV/c (%s)\n", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPT-1], hltTag);
  else 
    sprintf(name, "J/#psi: %1.2f <|y|< %1.2f, %1.1f < p_{T} < %1.1f GeV/c (%s)", jpsi::rapForPTRange[iRap-1], jpsi::rapForPTRange[iRap], jpsi::pTRange[iPT-1], jpsi::pTRange[iPT], hltTag);
  hAccErr2D_pol_pT_rap[iFrame][iPT][iRap]->SetTitle(name);
  hAccErr2D_pol_pT_rap[iFrame][iPT][iRap]->SetZTitle("rel. stat. err. on Acc*Eff");
  hAccErr2D_pol_pT_rap[iFrame][iPT][iRap]->SetMaximum(0.25);
  hAccErr2D_pol_pT_rap[iFrame][iPT][iRap]->SetTitleOffset(1.5, "z");
  //hAccErr2D_pol_pT_rap[iFrame][iPT][iRap]->Draw("colz");//"colz" or "cont"
  hAccErr2D_pol_pT_rap[iFrame][iPT][iRap]->Draw("lego2");//"colz" or "cont"

  sprintf(name, "Figures/accErr2D_%s_%s_%s_rap%d_pT%d.eps", speciesName[kNbSpecies], jpsi::frameLabel[iFrame], hltTag, iRap, iPT); 
  c31_2D->Print(name);		 
  sprintf(name, "Figures/accErr2D_%s_%s_%s_rap%d_pT%d.pdf", speciesName[kNbSpecies], jpsi::frameLabel[iFrame], hltTag, iRap, iPT); 
  c31_2D->Print(name);		 
  sprintf(name, "Figures/accErr2D_%s_%s_%s_rap%d_pT%d.gif", speciesName[kNbSpecies], jpsi::frameLabel[iFrame], hltTag, iRap, iPT); 
  c31_2D->Print(name);
}

//===============================
void PlotAllGenRec(Int_t iFrame){

  gStyle->SetOptStat(0);

  Char_t name[100];
  // //===============================================================
  // //cosTheta and phi, generated and reconstructed for pT dependence
  // //===============================================================
  // sprintf(name, "c1AllGenRec_PT_%s", jpsi::frameLabel[iFrame]);
  // TCanvas *c1CosTh = new TCanvas(name, "generated and reco histos for PT Bins", 1000, 700);
  // c1CosTh->Divide(2,2);
  // TH1F *hFrame1a[4];
  // for(int i = 0; i < 4; i++){
  //   c1CosTh->cd(i+1);
  //   if(i == 0){
  //     sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);
  //     hFrame1a[i] = gPad->DrawFrame(-1, 0., 1, 1.2*hGen_pol_pT[iFrame][0][cosThPol]->GetMaximum());
  //     hFrame1a[i]->SetYTitle("Generated");
  //   }
  //   else if(i == 1){
  //     sprintf(name, "#phi_{%s}", jpsi::frameLabel[iFrame]);
  //     hFrame1a[i] = gPad->DrawFrame(0., 0., 360., 1.2*hGen_pol_pT[iFrame][0][phiPol]->GetMaximum());
  //     hFrame1a[i]->SetYTitle("Generated");
  //   }
  //   else if(i == 2){
  //     gPad->SetLogy();
  //     sprintf(name, "cos#theta_{%s}", jpsi::frameLabel[iFrame]);
  //     hFrame1a[i] = gPad->DrawFrame(-1, 100., 1, 3.*Reco_pol_pT[iFrame][0][cosThPol]->GetMaximum());
  //     hFrame1a[i]->SetYTitle("Reconstructed");
  //   }
  //   else if(i == 3){
  //     gPad->SetLogy();
  //     sprintf(name, "#phi_{%s}", jpsi::frameLabel[iFrame]);
  //     hFrame1a[i] = gPad->DrawFrame(0., 100., 360., 3.*Reco_pol_pT[iFrame][0][phiPol]->GetMaximum());
  //     hFrame1a[i]->SetYTitle("Reconstructed");
  //   }
  //   hFrame1a[i]->SetXTitle(name);
  //   for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  //     if(i == 0) hGen_pol_pT[iFrame][iPTBin][cosThPol]->Draw("psame");    
  //     else if(i == 1) hGen_pol_pT[iFrame][iPTBin][phiPol]->Draw("psame");    
  //     else if(i == 2) Reco_pol_pT[iFrame][iPTBin][cosThPol]->Draw("psame");    
  //     else if(i == 3) Reco_pol_pT[iFrame][iPTBin][phiPol]->Draw("psame");    
  //   }
  //   if(i == 1){
  //     TLegend *leg1a = new TLegend(0.6129518,0.4317956,0.9497155,0.7449157);
  //     for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
  // 	if(iPTBin == 0) sprintf(name, "all p_{T}");
  // 	else if(iPTBin == jpsi::kNbPTBins) sprintf(name, "p_{T} > %1.1f GeV/c\n", jpsi::pTRange[iPTBin-1]);
  // 	else sprintf(name, "%1.1f < p_{T} < %1.1f GeV/c", jpsi::pTRange[iPTBin-1], jpsi::pTRange[iPTBin]);
  // 	leg1a->AddEntry(Reco_pol_pT[iFrame][iPTBin][phiPol], name, "pl");
  //     }
  //     leg1a->SetTextSize(0.035); leg1a->SetFillColor(0);
  //     leg1a->Draw();
  //   }
  // }
  // sprintf(name, "Figures/genReco_pTBins_%s.eps", jpsi::frameLabel[iFrame]); c1CosTh->Print(name);
  // sprintf(name, "Figures/genReco_pTBins_%s.pdf", jpsi::frameLabel[iFrame]); c1CosTh->Print(name);
  // sprintf(name, "Figures/genReco_pTBins_%s.gif", jpsi::frameLabel[iFrame]); c1CosTh->Print(name);

}
//===============================
void PlotKinVarAcceptance(){

  // Char_t name[100];
  // TCanvas *cpT = new TCanvas("cpT", "pT acceptance");
  // TH1F *hFramePT = gPad->DrawFrame(0., 0., 30., 1.);
  // hFramePT->SetXTitle("p_{T}(J/#psi) [GeV/c]");
  // hFramePT->SetYTitle("Acceptance");
  // for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++)
  //   gAcc_pt[iRap]->Draw("psame");
  // cpT->Print("Figures/acceptance_pT.eps");
  // cpT->Print("Figures/acceptance_pT.pdf");
  // cpT->Print("Figures/acceptance_pT.gif");

  // TCanvas *cRap = new TCanvas("cRap", "Rap acceptance");
  // TH1F *hFrameRap = gPad->DrawFrame(-2.5, 0., 2.5, 1.);
  // hFrameRap->SetXTitle("y(J/#psi)");
  // hFrameRap->SetYTitle("Acceptance");
  // for(int iPT = 0; iPT < jpsi::kNbPTBins+1; iPT++)
  //   gAcc_rap[iPT]->Draw("psame");
  // cRap->Print("Figures/acceptance_y.eps");
  // cRap->Print("Figures/acceptance_y.pdf");
  // cRap->Print("Figures/acceptance_y.gif");

  // TCanvas *cPhi[jpsi::kNbRapForPTBins+1];
  // for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++){
  //   sprintf(name, "cPhi_rap%d", iRap);
  //   cPhi[iRap] = new TCanvas(name, "Phi acceptance");
  //   TH1F *hFramePhi = gPad->DrawFrame(-3.14, 0., 3.14, 1.);
  //   hFramePhi->SetXTitle("#phi(J/#psi)");
  //   hFramePhi->SetYTitle("Acceptance");
  //   for(int iPT = 0; iPT < jpsi::kNbPTBins+1; iPT++)
  //     gAcc_phi[iPT][iRap]->Draw("psame");

  //   sprintf(name, "Figures/acceptance_phi_rap%d.eps", iRap); cPhi[iRap]->Print(name);
  //   sprintf(name, "Figures/acceptance_phi_rap%d.pdf", iRap); cPhi[iRap]->Print(name);
  //   sprintf(name, "Figures/acceptance_phi_rap%d.gif", iRap); cPhi[iRap]->Print(name);
  // }
}

//===============================
void CalcAcceptance(){

  // Bool_t dividing = kTRUE;
  Bool_t dividing = kFALSE;
  Char_t name[100];

  // printf("preparing pT differential acceptance\n");
  // for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++){
  //   gAcc_pt[iRap] = new TGraphAsymmErrors();
  //   gAcc_pt[iRap]->BayesDivide(Reco_pt[iRap], hGen_pt[iRap]);  
  //   sprintf(name, "gAcc_pt_rap%d", iRap);
  //   gAcc_pt[iRap]->SetName(name);
  // }
  // printf("preparing rap differential acceptance\n");
  // for(int iPT = 0; iPT < jpsi::kNbPTBins+1; iPT++){
  //   gAcc_rap[iPT] = new TGraphAsymmErrors();
  //   gAcc_rap[iPT]->BayesDivide(Reco_rap[iPT], hGen_rap[iPT]); 
  //   sprintf(name, "gAcc_rap_pT%d", iPT);
  //   gAcc_rap[iPT]->SetName(name);
  // }
  // printf("preparing phi differential acceptance\n");
  // for(int iPT = 0; iPT < jpsi::kNbPTBins+1; iPT++){
  //   for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++){

  //     gAcc_phi[iPT][iRap] = new TGraphAsymmErrors();
  //     Reco_phi[iPT][iRap]->Rebin(4); 
  //     hGen_phi[iPT][iRap]->Rebin(4);
  //     gAcc_phi[iPT][iRap]->BayesDivide(Reco_phi[iPT][iRap], hGen_phi[iPT][iRap]); 
  //     sprintf(name, "gAcc_phi_pT%d_rap%d", iPT, iRap);
  //     gAcc_phi[iPT][iRap]->SetName(name);
  //   }
  // }

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    // for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){

    //   printf("preparing cosTheta 1D differential acceptance (pT bin %d)\n", iPTBin);
    //   gAcc_pol_pT[iFrame][iPTBin][cosThPol] = new TGraphAsymmErrors();
    //   gAcc_pol_pT[iFrame][iPTBin][cosThPol]->BayesDivide(Reco_pol_pT[iFrame][iPTBin][cosThPol],
    //   							 hGen_pol_pT[iFrame][iPTBin][cosThPol]);
    //   sprintf(name, "gAcc_Onia_cosTh_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
    //   gAcc_pol_pT[iFrame][iPTBin][cosThPol]->SetName(name);
    //   gAcc_pol_pT[iFrame][iPTBin][cosThPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
    //   gAcc_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
    //   gAcc_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    //   gAcc_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerSize(markerSize);
    //   //

    //   printf("preparing phi 1D differential acceptance (pT bin %d)\n", iPTBin);
    //   gAcc_pol_pT[iFrame][iPTBin][phiPol] = new TGraphAsymmErrors();
    //   gAcc_pol_pT[iFrame][iPTBin][phiPol]->BayesDivide(Reco_pol_pT[iFrame][iPTBin][phiPol],
    //   						       hGen_pol_pT[iFrame][iPTBin][phiPol]);
    //   sprintf(name, "gAcc_Onia_phi_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
    //   gAcc_pol_pT[iFrame][iPTBin][phiPol]->SetName(name);
    //   gAcc_pol_pT[iFrame][iPTBin][phiPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
    //   gAcc_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
    //   gAcc_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    //   gAcc_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerSize(markerSize);
    //   //
    //   sprintf(name, "hAcc2D_Onia_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
    //   hAcc2D_pol_pT[iFrame][iPTBin] = (TH2D *) Reco2D_pol_pT[iFrame][iPTBin]->Clone(name);
    //   if(dividing)
    //   	hAcc2D_pol_pT[iFrame][iPTBin]->Divide(hGen2D_pol_pT[iFrame][iPTBin]);
    //   hAcc2D_pol_pT[iFrame][iPTBin]->SetLineColor(jpsi::colour_pT[iPTBin]);
    //   hAcc2D_pol_pT[iFrame][iPTBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
    //   hAcc2D_pol_pT[iFrame][iPTBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    //   hAcc2D_pol_pT[iFrame][iPTBin]->SetMarkerSize(markerSize);
    // }
    // for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){

    //   printf("preparing cosTheta 1D differential acceptance (rap bin %d)\n", iRapBin);
    //   gAcc_pol_rap[iFrame][iRapBin][cosThPol] = new TGraphAsymmErrors();
    //   gAcc_pol_rap[iFrame][iRapBin][cosThPol]->BayesDivide(Reco_pol_rap[iFrame][iRapBin][cosThPol],
    //   							   hGen_pol_rap[iFrame][iRapBin][cosThPol]);
    //   sprintf(name, "gAcc_Onia_cosTh_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
    //   gAcc_pol_rap[iFrame][iRapBin][cosThPol]->SetName(name);
    //   gAcc_pol_rap[iFrame][iRapBin][cosThPol]->SetLineColor(jpsi::colour_rap[iRapBin]);
    //   gAcc_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerColor(jpsi::colour_rap[iRapBin]);
    //   gAcc_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerStyle(jpsi::marker_rap[iRapBin]);
    //   gAcc_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerSize(markerSize);           
    //   //
    //   printf("preparing phi 1D differential acceptance (rap bin %d)\n", iRapBin);
    //   gAcc_pol_rap[iFrame][iRapBin][phiPol] = new TGraphAsymmErrors();
    //   gAcc_pol_rap[iFrame][iRapBin][phiPol]->BayesDivide(Reco_pol_rap[iFrame][iRapBin][phiPol],
    //   							 hGen_pol_rap[iFrame][iRapBin][phiPol]);
    //   sprintf(name, "gAcc_Onia_phi_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
    //   gAcc_pol_rap[iFrame][iRapBin][phiPol]->SetName(name);
    //   gAcc_pol_rap[iFrame][iRapBin][phiPol]->SetLineColor(jpsi::colour_rap[iRapBin]);
    //   gAcc_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerColor(jpsi::colour_rap[iRapBin]);
    //   gAcc_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerStyle(jpsi::marker_rap[iRapBin]);
    //   gAcc_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerSize(markerSize);
    //   //H: implement by hand!!!
    //   sprintf(name, "hAcc2D_Onia_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
    //   hAcc2D_pol_rap[iFrame][iRapBin] = (TH2D *) Reco2D_pol_rap[iFrame][iRapBin]->Clone(name);
    //   if(dividing)
    //   	hAcc2D_pol_rap[iFrame][iRapBin]->Divide(hGen2D_pol_rap[iFrame][iRapBin]);
    //   hAcc2D_pol_rap[iFrame][iRapBin]->SetLineColor(jpsi::colour_rap[iRapBin]);
    //   hAcc2D_pol_rap[iFrame][iRapBin]->SetMarkerColor(jpsi::colour_rap[iRapBin]);
    //   hAcc2D_pol_rap[iFrame][iRapBin]->SetMarkerStyle(jpsi::marker_rap[iRapBin]);
    //   hAcc2D_pol_rap[iFrame][iRapBin]->SetMarkerSize(markerSize);
    // }
    // for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    //   for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){

    // 	printf("preparing cosTheta 1D differential acceptance (pT bin %d, rap bin %d)\n", iPTBin, iRapBin);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol] = new TGraphAsymmErrors();
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->BayesDivide(Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol],
    // 									hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]);
    // 	sprintf(name, "gAcc_Onia_cosTh_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetName(name);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerSize(markerSize);         
    // 	//

    // 	printf("preparing phi 1D differential acceptance (pT bin %d, rap bin %d)\n", iPTBin, iRapBin);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol] = new TGraphAsymmErrors();
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->BayesDivide(Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol],
    // 								      hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]);
    // 	sprintf(name, "gAcc_Onia_phi_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetName(name);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerSize(markerSize);

    // 	//H: implement by hand!!!
    // 	printf("preparing 2D differential acceptance (pT bin %d, rap bin %d)\n", iPTBin, iRapBin);
    // 	gAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TGraphAsymmErrors();
    // 	sprintf(name, "gAcc2D_Onia_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
    // 	gAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetName(name);
    // 	Calc2DAcc(Reco2D_pol_pT_rap[iFrame][iPTBin][iRapBin], 
    // 		  hGen2D_pol_pT_rap[iFrame][iPTBin][iRapBin],
    // 		  gAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]);
    //   }
    // }
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){

	sprintf(name, "hAcc2D_Onia_%s_pT%d_rap%d_%s", jpsi::frameLabel[iFrame], iPTBin, iRapBin, speciesName[kNbSpecies]);
	hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) Reco2D_pol_pT_rap[kNbSpecies][iFrame][iPTBin][iRapBin]->Clone(name);
	//prepare a histogram containing the errors of the acceptance:
	sprintf(name, "hAccErr2D_Onia_%s_pT%d_rap%d_%s", jpsi::frameLabel[iFrame], iPTBin, iRapBin, speciesName[kNbSpecies]);
	hAccErr2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2D *) Reco2D_pol_pT_rap[0][iFrame][iPTBin][iRapBin]->Clone(name);
	hAccErr2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Reset();

	Double_t nEntries;
	for(int iX = 1; iX <= hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(); iX++){
	  for(int iY = 1; iY <= hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY(); iY++){
	    nEntries = 0.;
	    for(int iSpecies = 0; iSpecies < kNbSpecies; iSpecies++)
	      //make the check on the *unscaled* histograms
	      nEntries += Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->GetBinContent(iX, iY);
	    //in case a bin has less than N events, set the acceptance to 0
	    //and increase the bin error to something very large
	    // if(nEntries < minEntriesPerBin){
	    //   hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(iX, iY, 0.);
	    //   hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinError(iX, iY, 100.);
	    // }
	    if(!dividing){//pure Possonian error:
	      if(nEntries > 0.)
		hAccErr2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(iX,iY,1./sqrt(nEntries));
	    }
	  }
	}
	if(dividing)
	  hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Divide(hGen2D_pol_pT_rap[kNbSpecies][iFrame][iPTBin][iRapBin]);
	hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetLineColor(jpsi::colour_pT[iPTBin]);
	hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
	hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
	hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetMarkerSize(markerSize);
	printf("%s: pT %d, rap %d: max. acceptance: %1.3e\n", jpsi::frameLabel[iFrame], iPTBin, iRapBin,
	       hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetMaximum());

	// //in case a bin content is higher than 1: set it back to 1
	// Double_t content;
	// for(int iX = 1; iX <= hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsX(); iX++){
	//   for(int iY = 1; iY <= hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetNbinsY(); iY++){
	//     content = hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetBinContent(iX, iY);
	//     if(content > 1.)
	//       hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->SetBinContent(iX, iY, 1.);
	//   }
	// }
      }
    }
  }
}

//===============================
void ReadInHistos(Char_t *fileNameIn, Int_t iSpecies, Int_t rebinCosTh, Int_t rebinPhi){

  TFile *fIn = new TFile(fileNameIn);
  Char_t name[100];

  // for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++){
  //   sprintf(name, "hGen_Onia_pt_rap%d", iRap);
  //   hGen_pt[iRap] = (TH1D *)  gDirectory->Get(name);
  //   sprintf(name, "Reco_Onia_pt_rap%d", iRap);
  //   Reco_pt[iRap] = (TH1D *)  gDirectory->Get(name);
  // }
  // for(int iPT = 0; iPT < jpsi::kNbPTBins+1; iPT++){
  //   sprintf(name, "hGen_Onia_rap_pT%d", iPT);
  //   hGen_rap[iPT] = (TH1D *) gDirectory->Get(name);
  //   sprintf(name, "Reco_Onia_rap_pT%d", iPT);
  //   Reco_rap[iPT] = (TH1D *) gDirectory->Get(name);
  // }
  // for(int iPT = 0; iPT < jpsi::kNbPTBins+1; iPT++){
  //   for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++){
  //     sprintf(name, "hGen_Onia_phi_pT%d_rap%d", iPT, iRap);
  //     hGen_phi[iPT][iRap] = (TH1D *) gDirectory->Get(name);
  //     sprintf(name, "Reco_Onia_phi_pT%d_rap%d", iPT, iRap);
  //     Reco_phi[iPT][iRap] = (TH1D *) gDirectory->Get(name);
  //   }
  // }

  Int_t totEntries[jpsi::kNbFrames];
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    totEntries[iFrame] = 0;
    // for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    //   //generated
    //   sprintf(name, "hGen_Onia_cosTh_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
    //   hGen_pol_pT[iFrame][iPTBin][cosThPol] = (TH1D *) gDirectory->Get(name);
    //   hGen_pol_pT[iFrame][iPTBin][cosThPol]->SetLineColor(jpsi::colour_pT[iPTBin]);
    //   hGen_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
    //   hGen_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    //   hGen_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerSize(markerSize);        
    //   sprintf(name, "hGen_Onia_phi_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
    //   hGen_pol_pT[iFrame][iPTBin][phiPol] = (TH1D *) gDirectory->Get(name);
    //   hGen_pol_pT[iFrame][iPTBin][phiPol]->SetLineColor(jpsi::colour_pT[iPTBin]);  
    //   hGen_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
    //   hGen_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    //   hGen_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerSize(markerSize);        
    //   sprintf(name, "hGen2D_Onia_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
    //   hGen2D_pol_pT[iFrame][iPTBin] = (TH2D *) gDirectory->Get(name);
    //   hGen2D_pol_pT[iFrame][iPTBin]->SetLineColor(jpsi::colour_pT[iPTBin]);  
    //   hGen2D_pol_pT[iFrame][iPTBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
    //   hGen2D_pol_pT[iFrame][iPTBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    //   hGen2D_pol_pT[iFrame][iPTBin]->SetMarkerSize(markerSize);        
    //   if(rebinCosTh > 1 || rebinPhi > 1)
    // 	hGen2D_pol_pT[iFrame][iPTBin]->Rebin2D(rebinCosTh, rebinPhi);

    //   //reco
    //   sprintf(name, "Reco_Onia_cosTh_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
    //   Reco_pol_pT[iFrame][iPTBin][cosThPol] = (TH1D *) gDirectory->Get(name);
    //   Reco_pol_pT[iFrame][iPTBin][cosThPol]->SetLineColor(jpsi::colour_pT[iPTBin]);  
    //   Reco_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
    //   Reco_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    //   Reco_pol_pT[iFrame][iPTBin][cosThPol]->SetMarkerSize(markerSize);        
    //   sprintf(name, "Reco_Onia_phi_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
    //   Reco_pol_pT[iFrame][iPTBin][phiPol] = (TH1D *) gDirectory->Get(name);
    //   Reco_pol_pT[iFrame][iPTBin][phiPol]->SetLineColor(jpsi::colour_pT[iPTBin]);  
    //   Reco_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
    //   Reco_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    //   Reco_pol_pT[iFrame][iPTBin][phiPol]->SetMarkerSize(markerSize);
    //   sprintf(name, "Reco2D_Onia_%s_pT%d", jpsi::frameLabel[iFrame], iPTBin);
    //   Reco2D_pol_pT[iFrame][iPTBin] = (TH2D *) gDirectory->Get(name);
    //   Reco2D_pol_pT[iFrame][iPTBin]->SetLineColor(jpsi::colour_pT[iPTBin]);  
    //   Reco2D_pol_pT[iFrame][iPTBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]);
    //   Reco2D_pol_pT[iFrame][iPTBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);
    //   Reco2D_pol_pT[iFrame][iPTBin]->SetMarkerSize(markerSize);
    //   if(rebinCosTh > 1 || rebinPhi > 1)
    // 	Reco2D_pol_pT[iFrame][iPTBin]->Rebin2D(rebinCosTh, rebinPhi);
    // }
    // for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
    //   //generated
    //   sprintf(name, "hGen_Onia_cosTh_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
    //   hGen_pol_rap[iFrame][iRapBin][cosThPol] = (TH1D *) gDirectory->Get(name);
    //   hGen_pol_rap[iFrame][iRapBin][cosThPol]->SetLineColor(jpsi::colour_rap[iRapBin]);   
    //   hGen_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerColor(jpsi::colour_rap[iRapBin]); 
    //   hGen_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerStyle(jpsi::marker_rap[iRapBin]); 
    //   hGen_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerSize(markerSize);           
    //   sprintf(name, "hGen_Onia_phi_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
    //   hGen_pol_rap[iFrame][iRapBin][phiPol] = (TH1D *) gDirectory->Get(name);
    //   hGen_pol_rap[iFrame][iRapBin][phiPol]->SetLineColor(jpsi::colour_rap[iRapBin]);   
    //   hGen_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerColor(jpsi::colour_rap[iRapBin]); 
    //   hGen_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerStyle(jpsi::marker_rap[iRapBin]); 
    //   hGen_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerSize(markerSize);           
    //   sprintf(name, "hGen2D_Onia_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
    //   hGen2D_pol_rap[iFrame][iRapBin] = (TH2D *) gDirectory->Get(name);
    //   hGen2D_pol_rap[iFrame][iRapBin]->SetLineColor(jpsi::colour_rap[iRapBin]);   
    //   hGen2D_pol_rap[iFrame][iRapBin]->SetMarkerColor(jpsi::colour_rap[iRapBin]); 
    //   hGen2D_pol_rap[iFrame][iRapBin]->SetMarkerStyle(jpsi::marker_rap[iRapBin]); 
    //   hGen2D_pol_rap[iFrame][iRapBin]->SetMarkerSize(markerSize);           
    //   if(rebinCosTh > 1 || rebinPhi > 1)
    // 	hGen2D_pol_rap[iFrame][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
    //   //reconstructed
    //   sprintf(name, "Reco_Onia_cosTh_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
    //   Reco_pol_rap[iFrame][iRapBin][cosThPol] = (TH1D *) gDirectory->Get(name);
    //   Reco_pol_rap[iFrame][iRapBin][cosThPol]->SetLineColor(jpsi::colour_rap[iRapBin]);   
    //   Reco_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerColor(jpsi::colour_rap[iRapBin]); 
    //   Reco_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerStyle(jpsi::marker_rap[iRapBin]); 
    //   Reco_pol_rap[iFrame][iRapBin][cosThPol]->SetMarkerSize(markerSize);           
    //   sprintf(name, "Reco_Onia_phi_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
    //   Reco_pol_rap[iFrame][iRapBin][phiPol] = (TH1D *) gDirectory->Get(name);
    //   Reco_pol_rap[iFrame][iRapBin][phiPol]->SetLineColor(jpsi::colour_rap[iRapBin]);   
    //   Reco_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerColor(jpsi::colour_rap[iRapBin]); 
    //   Reco_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerStyle(jpsi::marker_rap[iRapBin]); 
    //   Reco_pol_rap[iFrame][iRapBin][phiPol]->SetMarkerSize(markerSize);           
    //   sprintf(name, "Reco2D_Onia_%s_rap%d", jpsi::frameLabel[iFrame], iRapBin);
    //   Reco2D_pol_rap[iFrame][iRapBin] = (TH2D *) gDirectory->Get(name);
    //   Reco2D_pol_rap[iFrame][iRapBin]->SetLineColor(jpsi::colour_rap[iRapBin]);   
    //   Reco2D_pol_rap[iFrame][iRapBin]->SetMarkerColor(jpsi::colour_rap[iRapBin]); 
    //   Reco2D_pol_rap[iFrame][iRapBin]->SetMarkerStyle(jpsi::marker_rap[iRapBin]); 
    //   Reco2D_pol_rap[iFrame][iRapBin]->SetMarkerSize(markerSize);           
    //   if(rebinCosTh > 1 || rebinPhi > 1)
    // 	Reco2D_pol_rap[iFrame][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
    // }
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
	//generated
	// sprintf(name, "hGen_Onia_cosTh_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	// hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol] = (TH1D *) gDirectory->Get(name);
	// hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetLineColor(jpsi::colour_pT[iPTBin]);   
	// hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]);  
	// hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]);  
	// hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerSize(markerSize);          
	// sprintf(name, "hGen_Onia_phi_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	// hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol] = (TH1D *) gDirectory->Get(name);
	// hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetLineColor(jpsi::colour_pT[iPTBin]);   
	// hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]); 
	// hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]); 
	// hGen_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerSize(markerSize);         
	sprintf(name, "hGen2D_Onia_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	hGen2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin] = (TH2D *) gDirectory->Get(name);
	sprintf(name, "%s_%s", hGen2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->GetName(), speciesName[iSpecies]);
	hGen2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->SetName(name);
	hGen2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->SetLineColor(jpsi::colour_pT[iPTBin]);   
	hGen2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]); 
	hGen2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]); 
	hGen2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->SetMarkerSize(markerSize);         
	if(rebinCosTh > 1 || rebinPhi > 1)
	  hGen2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
	printf("%s: pTbin %d, rapBin %d --> statistics in  GEN: %1.1f\n",
	       jpsi::frameLabel[iFrame], iPTBin, iRapBin, 
	       hGen2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->GetEntries());

	//reconstructed
	// sprintf(name, "Reco_Onia_cosTh_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	// Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol] = (TH1D *) gDirectory->Get(name);
	// Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetLineColor(jpsi::colour_pT[iPTBin]);   
	// Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]); 
	// Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]); 
	// Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->SetMarkerSize(markerSize);         
	// sprintf(name, "Reco_Onia_phi_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	// Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol] = (TH1D *) gDirectory->Get(name);
	// Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetLineColor(jpsi::colour_pT[iPTBin]);   
	// Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerColor(jpsi::colour_pT[iPTBin]); 
	// Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerStyle(jpsi::marker_pT[iPTBin]); 
	// Reco_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->SetMarkerSize(markerSize);         
	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin] = (TH2D *) gDirectory->Get(name);
	sprintf(name, "%s_%s", Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->GetName(), speciesName[iSpecies]);
	Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->SetName(name);
	Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->SetLineColor(jpsi::colour_pT[iPTBin]);   
	Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->SetMarkerColor(jpsi::colour_pT[iPTBin]); 
	Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->SetMarkerStyle(jpsi::marker_pT[iPTBin]); 
	Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->SetMarkerSize(markerSize);         
	printf("%s: pTbin %d, rapBin %d --> statistics in RECO: %1.1f\n",
	       jpsi::frameLabel[iFrame], iPTBin, iRapBin, 
	       Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->GetEntries());

	totEntries[iFrame] += Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->GetEntries();
	if(rebinCosTh > 1 || rebinPhi > 1)
	  Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
      }
    }
    printf("%s: integrated statistics in RECO is %d\n", 
	   jpsi::frameLabel[iFrame], totEntries[iFrame]);
  }
}

//===============================
void WriteAccHistos(Char_t *fileNameOut){

  printf("<WriteAccHistos> writing out the histograms\n");

  TFile *fOut = new TFile(fileNameOut, "RECREATE");

  // for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++)
  //   gAcc_pt[iRap]->Write();
  // for(int iPT = 0; iPT < jpsi::kNbPTBins+1; iPT++)
  //   gAcc_rap[iPT]->Write();
  // for(int iPT = 0; iPT < jpsi::kNbPTBins+1; iPT++)
  //   for(int iRap = 0; iRap < jpsi::kNbRapForPTBins+1; iRap++)
  //     gAcc_phi[iPT][iRap]->Write();

  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    // for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    //   gAcc_pol_pT[iFrame][iPTBin][cosThPol]->Write();
    //   gAcc_pol_pT[iFrame][iPTBin][phiPol]->Write();
    //   hAcc2D_pol_pT[iFrame][iPTBin]->Write();
    // }
    // for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++){
    //   gAcc_pol_rap[iFrame][iRapBin][cosThPol]->Write();
    //   gAcc_pol_rap[iFrame][iRapBin][phiPol]->Write();
    //   hAcc2D_pol_rap[iFrame][iRapBin]->Write();
    // }
    // for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
    //   for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][cosThPol]->Write();
    // 	gAcc_pol_pT_rap[iFrame][iPTBin][iRapBin][phiPol]->Write();
    // 	gAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
    //   }
    // }
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
	hAccErr2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
	hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }
  fOut->Close();
}

//=======================================================
void Calc2DAcc(TH2D *hReco, TH2D *hGen, TGraphAsymmErrors *graph){

  Int_t nBinsX1 = hReco->GetNbinsX();
  Int_t nBinsX2 = hGen->GetNbinsX();
  Int_t nBinsY1 = hReco->GetNbinsY();
  Int_t nBinsY2 = hGen->GetNbinsY();
  if(nBinsX1 != nBinsX2 && nBinsY1 != nBinsY2) return;

  //build 1D histos containing the slices of the 2D histo appended
  Char_t name[100];
  sprintf(name, "%s_1D", hGen->GetName());
  TH1D *hGen1D = new TH1D(name, "", nBinsX1*nBinsY1, 0., nBinsX1*nBinsY1);
  sprintf(name, "%s_1D", hReco->GetName());
  TH1D *hReco1D = new TH1D(name, "", nBinsX1*nBinsY1, 0., nBinsX1*nBinsY1);
  Double_t gen, rec;
  for(int iY = 1; iY <= nBinsY1; iY++){
    for(int iX = 1; iX <= nBinsX1; iX++){
      gen = hGen->GetBinContent(iX, iY);
      rec = hReco->GetBinContent(iX, iY);

      hGen1D->SetBinContent((iY-1)*nBinsX1 + iX, gen);
      hReco1D->SetBinContent((iY-1)*nBinsX1 + iX, rec);
    }
  }

  //use Bayes divide to get the errors
  graph->BayesDivide(hReco1D, hGen1D);
  
  // //***use FeldmanCousinsBinomialInterval/ClopperPearsonBinomialInterval to get the errors
  // //BinomialInterval::FeldmanCousinsBinomialInterval fc;
  // FeldmanCousinsBinomialInterval fc;
  // //ClopperPearsonBinomialInterval cp;

  // //alpha = 1 - CL
  // const double alpha = (1-0.682);
  // fc.init(alpha);
  // //cp.init(alpha);

  // Int_t binX = hGen1D->GetNbinsX();

  // for(int iX = 1; iX <= binX; iX++){
  //     gen = hGen1D->GetBinContent(iX);
  //     rec = hReco1D->GetBinContent(iX);

  //     double acc = rec / gen;
  //     fc.calculate(rec, gen);
  //     //cp.calculate(rec, gen);
  //     double errorLow = acc - fc.lower();
  //     double errorHigh = fc.upper() - acc;

  //     Double_t xVal =  hGen1D->GetBinCenter(iX);

  //     graph->SetPoint(iX,xVal,acc);
  //     graph->SetPointEYhigh(iX, errorHigh);
  //     graph->SetPointEYlow(iX, errorLow);
  // }

}

//======================================
void AddSpecies(){

  Char_t name[100];
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins+1; iPTBin++){
      for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
	//generated histograms:
	sprintf(name, "hGen2D_Onia_%s_pT%d_rap%d_%s", jpsi::frameLabel[iFrame], iPTBin, iRapBin, speciesName[kNbSpecies]);
	hGen2D_pol_pT_rap[kNbSpecies][iFrame][iPTBin][iRapBin] = (TH2D *) hGen2D_pol_pT_rap[0][iFrame][iPTBin][iRapBin]->Clone(name);
	hGen2D_pol_pT_rap[kNbSpecies][iFrame][iPTBin][iRapBin]->Scale(weights[0]);

	for(int iSpecies = 1; iSpecies < kNbSpecies; iSpecies++)
	  hGen2D_pol_pT_rap[kNbSpecies][iFrame][iPTBin][iRapBin]->Add(hGen2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin], weights[iSpecies]);

	//reconstructed histograms:
	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d_%s", jpsi::frameLabel[iFrame], iPTBin, iRapBin, speciesName[kNbSpecies]);
	Reco2D_pol_pT_rap[kNbSpecies][iFrame][iPTBin][iRapBin] = (TH2D *) Reco2D_pol_pT_rap[0][iFrame][iPTBin][iRapBin]->Clone(name);
	Reco2D_pol_pT_rap[kNbSpecies][iFrame][iPTBin][iRapBin]->Scale(weights[0]);

	for(int iSpecies = 1; iSpecies < kNbSpecies; iSpecies++)
	  Reco2D_pol_pT_rap[kNbSpecies][iFrame][iPTBin][iRapBin]->Add(Reco2D_pol_pT_rap[iSpecies][iFrame][iPTBin][iRapBin], weights[iSpecies]);
      }
    }
  }
}
