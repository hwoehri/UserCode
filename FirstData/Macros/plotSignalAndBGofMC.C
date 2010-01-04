#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "THStack.h"
#include "TLegend.h"

#include "commonVar.inc"

TH1F *hReco_mass[kNbCharge][kNbCath+2];
TH1F *hReco_pT[kNbCharge][kNbCath];
TH1F *hReco_rap[kNbCharge][kNbCath];
TH1F *hRecoBG_mass[kNbCharge][kNbCath+2];
TH1F *hRecoSum_mass[kNbCharge][kNbCath+2];//[kNbCath] contains sum of gl-gl and gl-tr
TH1F *hReco_eta1_eta2[kNbCharge][kNbCath];

Double_t massMinJPsi = 3.0, massMaxJPsi = 3.2;
// Double_t massMinJPsi = 2.9, massMaxJPsi = 3.3;

Char_t *catNames[kNbCath+2] = {"global-global", "global-tracker", "tracker-tracker",
			       "global-calo", "tracker-calo", "calo-calo", 
			       "gl-gl and gl-tr", "gl-gl, gl-tr and tr-tr"};

Double_t integral[kNbCath+2], errPossonian[kNbCath+2]; //signal counts
Double_t integralBG[kNbCath+2], errPossonianBG[kNbCath+2]; //signal counts
Double_t sigOvBG[kNbCath+2], errSigOvBG[kNbCath+2];

Int_t colour[kNbCath] = {1,2,4,3,6,7};

void LoadHistos(Char_t *fileNameIn, Bool_t takeBest);
void LoadBGHistos(Char_t *fileNameIn, Bool_t takeBest, Double_t lumi);
void Normalize(Bool_t normalize, Double_t intLumi);
void AddSignalBG();
void PlotHistos(Char_t *enLabel, Char_t *fileNameIn);
void PlotHistosWithBG(Char_t *enLabel, Char_t *fileNameIn);
Double_t gausexp(Double_t *x, Double_t *par);
void PrintStat();
//=======================
void plotSignalAndBGofMC(Char_t *fileNameIn = "histos_jPsi_2360GeV_STARTUP", //don't put ".root" at the end!
			 Char_t *enLabel = "2360 GeV",
			 Bool_t takeBest = kTRUE, //uses the histograms after all cuts
			 Bool_t normalize = kTRUE, //normalizes the histos to 1 nb-1 
			 Double_t intLumi = 261.//jPsi: DESIGN:  2.36 TeV: 266.; 900 GeV: 263.
			 //jPsi: STARTUP: 2.36 TeV: 261.; 900 GeV: 384.
			 ){

  LoadHistos(fileNameIn, takeBest);
  Normalize(normalize, intLumi);
  LoadBGHistos("histos_ppMuXLoose_2360GeV_STARTUP.root", takeBest, 5.5);

  PlotHistos(enLabel, fileNameIn);
  PrintStat();
}

//=======================
void PlotHistosWithBG(Char_t *enLabel, Char_t *fileNameIn){

  TCanvas *c1MassBG[kNbCath];
  TLatex *texMass[kNbCath];
  Char_t name[100];
  Double_t xPrint = 2.1, scale;
  Double_t normBG[kNbCath+1], slopeBG[kNbCath+1];
  TF1* gauss2[kNbCath+1];
  TF1* expBG[kNbCath+1];
  Int_t iCatToPlot[4] = {0,1,2,kNbCath};
  for(int iCat = 0; iCat < 4; iCat++){

    sprintf(name, "c1MassBG_%d", iCat);
    c1MassBG[iCat] = new TCanvas(name, "mass with BG");
    hRecoSum_mass[0][iCatToPlot[iCat]]->SetYTitle("dN/dM [per 100 MeV]");
    if(iCat == 3) 
      hRecoSum_mass[0][iCatToPlot[iCat]]->SetMaximum(1.7*hRecoSum_mass[0][iCatToPlot[iCat]]->GetMaximum());
    hRecoSum_mass[0][iCatToPlot[iCat]]->Draw();
    
    if(iCat > 1){
      //fit exponential and a Gauss:
      sprintf(name, "gauss2_cat%d", iCat);
      gauss2[iCatToPlot[iCat]] = new TF1(name,gausexp,2.0,4.0,5);
      gauss2[iCatToPlot[iCat]]->SetParNames("Signal","Mean","Sigma","Background","ExpSlope");
      gauss2[iCatToPlot[iCat]]->SetParameters(hRecoSum_mass[0][iCatToPlot[iCat]]->GetBinContent(11),3.1,0.04,
					      hRecoSum_mass[0][iCatToPlot[iCat]]->GetBinContent(1),-0.1);
      gauss2[iCatToPlot[iCat]]->SetParLimits(0,0.,3.*hRecoSum_mass[0][iCatToPlot[iCat]]->GetBinContent(11));
      gauss2[iCatToPlot[iCat]]->SetParLimits(1,2.9,3.2); //Mean
      gauss2[iCatToPlot[iCat]]->SetParLimits(2,0.,0.1); //Sigma
      //     gauss2[iCatToPlot[iCat]]->SetParLimits(3,0.,1000.); //BG
      gauss2[iCatToPlot[iCat]]->SetParLimits(4,-5.0,0.);//slope
      
      gauss2[iCatToPlot[iCat]]->FixParameter(1,3.096);
      gauss2[iCatToPlot[iCat]]->FixParameter(2,0.03);
      gauss2[iCatToPlot[iCat]]->SetLineWidth(1);    
      hRecoSum_mass[0][iCatToPlot[iCat]]->Fit(name,"","",2.0,3.9);
      //end fit
      gauss2[iCatToPlot[iCat]] = hRecoSum_mass[0][iCatToPlot[iCat]]->GetFunction(name);
      normBG[iCatToPlot[iCat]] =  gauss2[iCatToPlot[iCat]]->GetParameter(3);
      slopeBG[iCatToPlot[iCat]] =  gauss2[iCatToPlot[iCat]]->GetParameter(4);
      sprintf(name, "BG_%d", iCatToPlot[iCat]);
      expBG[iCatToPlot[iCat]] = new TF1(name, "[0]*exp([1]*x)", 2.0,4.0);
      expBG[iCatToPlot[iCat]]->FixParameter(0, normBG[iCatToPlot[iCat]]);
      expBG[iCatToPlot[iCat]]->FixParameter(1, slopeBG[iCatToPlot[iCat]]);
      expBG[iCatToPlot[iCat]]->SetLineColor(4);
      expBG[iCatToPlot[iCat]]->SetLineWidth(1);
      expBG[iCatToPlot[iCat]]->Draw("same");
    }

    if(iCat == 3) scale = 0.8;
    else scale = 1.0;

    texMass[iCat] = new TLatex(3.7, scale*1.*hRecoSum_mass[0][iCatToPlot[iCat]]->GetMaximum(), enLabel);
    texMass[iCat]->SetTextSize(0.04); texMass[iCat]->Draw();

    if(iCat == 2)
      xPrint = 2.4;
    else 
      xPrint = 2.1;
    texMass[iCat]->DrawLatex(xPrint, scale*1.0*hRecoSum_mass[0][iCatToPlot[iCat]]->GetMaximum(), catNames[iCatToPlot[iCat]]);
    sprintf(name, "#J/#psi per nb^{-1} :  %1.0f #pm %1.0f", 
	    integral[iCatToPlot[iCat]], errPossonian[iCatToPlot[iCat]]);
    texMass[iCat]->DrawLatex(xPrint, scale*0.9*hRecoSum_mass[0][iCatToPlot[iCat]]->GetMaximum(), name);

    if(iCatToPlot[iCat] == kNbCath)
      sprintf(name, "Figures/%s_withBG_mass_glglAndgltr.gif", fileNameIn);
    else
      sprintf(name, "Figures/%s_withBG_mass_%s.gif", fileNameIn, oniaCatName[iCatToPlot[iCat]]);
    c1MassBG[iCat]->Print(name);
    //pdf
    if(iCatToPlot[iCat] == kNbCath)
      sprintf(name, "Figures/%s_withBG_mass_glglAndgltr.pdf", fileNameIn);
    else
      sprintf(name, "Figures/%s_withBG_mass_%s.pdf", fileNameIn, oniaCatName[iCatToPlot[iCat]]);
    c1MassBG[iCat]->Print(name);
  }
}
//=======================
void Normalize(Bool_t normalize, Double_t intLumi){

  for(int iCat = 0; iCat < kNbCath+2; iCat++){
    Int_t binMin = hReco_mass[0][iCat]->GetXaxis()->FindBin(massMinJPsi);
    Int_t binMax = hReco_mass[0][iCat]->GetXaxis()->FindBin(massMaxJPsi) - 1;
    printf("taking integral for bins %d to %d\n", binMin, binMax);
    integral[iCat] = hReco_mass[0][iCat]->Integral(binMin, binMax);
    errPossonian[iCat] = sqrt(integral[iCat]);
  }
//   //add statistics from gl-gl and gl-tr
//   integral[kNbCath] = integral[0] + integral[1];
//   errPossonian[kNbCath] = sqrt(integral[kNbCath]);
//   //add statistics from gl-gl, gl-tr and tr-tr
//   integral[kNbCath+1] = integral[0] + integral[1] + integral[2];
//   errPossonian[kNbCath+1] = sqrt(integral[kNbCath]);

  for(int iCat = 0; iCat < kNbCath+2; iCat++){
    if(normalize){
      integral[iCat] /= intLumi;
      errPossonian[iCat]/= intLumi;
      //normalize now to the statistics available in the BG:
//       integral[iCat] *= 0.21;
//       errPossonian[iCat] *= 0.21;
    }
  }
  if(normalize){
//     for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
     for(int iCharge = 0; iCharge < 1; iCharge++){
      for(int iCat = 0; iCat < kNbCath+2; iCat++){
	hReco_mass[iCharge][iCat]->Scale(1./intLumi);
	if(iCat >= kNbCath) continue;
	hReco_pT[iCharge][iCat]->Scale(1./intLumi);
	hReco_rap[iCharge][iCat]->Scale(1./intLumi);
// 	hReco_mass[iCharge][iCat]->Scale(0.21/intLumi);
// 	hReco_pT[iCharge][iCat]->Scale(0.21/intLumi);
// 	hReco_rap[iCharge][iCat]->Scale(0.21/intLumi);
      }
    }
  }
}


//=======================
void PlotHistos(Char_t *enLabel, Char_t *fileNameIn){

  Bool_t plotBG = kTRUE;

  gStyle->SetOptStat(kFALSE);
//   gStyle->SetOptFit(kTRUE);

  Char_t name[100];

  //plot the mass histograms:
  TCanvas *c1Mass[kNbCath+2];
  TLatex *texMass[kNbCath+2];
  for(int iCat = 0; iCat < kNbCath+2; iCat++){

    if(iCat > 2 && iCat < kNbCath) continue;
    sprintf(name, "c1Mass_%d", iCat);
    c1Mass[iCat] = new TCanvas(name, "mass");

//     if(iCat == 2 && plotBG)
//       hReco_mass[0][iCat]->SetMaximum(5.);
    if(iCat == kNbCath+1 && plotBG)
      hReco_mass[0][iCat]->SetMaximum(12.);
    hReco_mass[0][iCat]->Draw();

    hReco_mass[0][iCat]->Fit("gaus", "", "", 3.05, 3.15);
//     (hReco_mass[0][iCat]->GetFunction("gaus"))->SetLineWidth(1);
    if(iCat < kNbCath)
      (hReco_mass[0][iCat]->GetFunction("gaus"))->SetLineColor(colour[iCat]);

    texMass[iCat] = new TLatex(3.7, 0.65*hReco_mass[0][iCat]->GetMaximum(), enLabel);
    texMass[iCat]->SetTextSize(0.04); texMass[iCat]->Draw();

    texMass[iCat]->DrawLatex(2.1, 0.6*hReco_mass[0][iCat]->GetMaximum(), catNames[iCat]);
//     sprintf(name, "#J/#psi per L_{int}=1 nb^{-1}:  %1.0f #pm %1.0f", 
// 	    integral[iCat], errPossonian[iCat]);
//     sprintf(name, "#J/#psi per L_{int}=1 nb^{-1}:  %1.0f", integral[iCat]);
     sprintf(name, "#J/#psi per nb^{-1}:  %1.0f", integral[iCat]);
    texMass[iCat]->DrawLatex(2.1, 0.5*hReco_mass[0][iCat]->GetMaximum(), name);

    if(plotBG)
      hRecoBG_mass[0][iCat]->Draw("same");


    if(iCat < kNbCath){
      sprintf(name, "Figures/%s_mass_%s.gif", fileNameIn, oniaCatName[iCat]);
      c1Mass[iCat]->Print(name);
      sprintf(name, "Figures/%s_mass_%s.pdf", fileNameIn, oniaCatName[iCat]);
      c1Mass[iCat]->Print(name);
    }
    else if(iCat == kNbCath){
      sprintf(name, "Figures/%s_mass_GGandGT.gif", fileNameIn);
      c1Mass[iCat]->Print(name);
      sprintf(name, "Figures/%s_mass_GGandGT.pdf", fileNameIn);
      c1Mass[iCat]->Print(name);
    }
    else if(iCat == kNbCath+1){
      sprintf(name, "Figures/%s_mass_all.gif", fileNameIn);
      c1Mass[iCat]->Print(name);
      sprintf(name, "Figures/%s_mass_all.pdf", fileNameIn);
      c1Mass[iCat]->Print(name);
    }
  }

  //plot the GG and GT combination with its BG AND
  //the GG with its BG overlayed:
  //plot the mass histograms:
  TCanvas *c1Massa[kNbCath+2];
  TLatex *texMassa[kNbCath+2];
  for(int iCat = kNbCath+1; iCat < kNbCath+2; iCat++){

    sprintf(name, "c1Massa_%d", iCat);
    c1Massa[iCat] = new TCanvas(name, "mass special");

    //gl-gl and gl-tr and tr-tr
    hReco_mass[0][iCat]->SetLineColor(4);
    (hReco_mass[0][iCat]->GetFunction("gaus"))->SetLineColor(4);
    hReco_mass[0][iCat]->Draw();
    //gl-gl and gl-tr
    hReco_mass[0][kNbCath]->SetLineColor(2);
    (hReco_mass[0][kNbCath]->GetFunction("gaus"))->SetLineColor(2);
    hReco_mass[0][kNbCath]->Draw("same");
    //gl-gl
//     hReco_mass[0][0]->SetLineColor(4);
//     (hReco_mass[0][0]->GetFunction("gaus"))->SetLineColor(4);
    hReco_mass[0][0]->Draw("same");

//     hReco_mass[0][iCat]->Fit("gaus", "", "", 3.05, 3.15);
//     hReco_mass[0][iCat]->SetLineWidth(1);

    texMassa[iCat] = new TLatex(3.6, 0.65*hReco_mass[0][iCat]->GetMaximum(), enLabel);
    texMassa[iCat]->SetTextSize(0.04); texMassa[iCat]->Draw();

    texMassa[iCat]->DrawLatex(2.45, 11., catNames[iCat]);
//     sprintf(name, "#J/#psi per L_{int}=1 nb^{-1}:  %1.0f #pm %1.0f", 
// 	    integral[iCat], errPossonian[iCat]);
//     sprintf(name, "#J/#psi per L_{int}=1 nb^{-1}:  %1.0f", integral[iCat]);
     sprintf(name, "#J/#psi per nb^{-1}:  %1.0f", integral[iCat]);
    texMassa[iCat]->DrawLatex(2.45, 9.8, name);

    if(plotBG){
      //gl-gl and gl-tr and tr-tr
      hRecoBG_mass[0][iCat]->SetLineColor(4);
      hRecoBG_mass[0][iCat]->Draw("same");
      //gl-gl and gl-tr
      hRecoBG_mass[0][kNbCath]->SetLineColor(2);
      hRecoBG_mass[0][kNbCath]->Draw("same");
      //gl-gl
      hRecoBG_mass[0][0]->SetLineColor(1);
      hRecoBG_mass[0][0]->Draw("same");
    }
    sprintf(name, "Figures/%s_mass_GGandGT_withGG.gif", fileNameIn);
    c1Massa[iCat]->Print(name);
    sprintf(name, "Figures/%s_mass_GGandGT_withGG.pdf", fileNameIn);
    c1Massa[iCat]->Print(name);
    TLatex *tex; 

    if(iCat == kNbCath+1){
      tex = new TLatex(3.231686,11.05442,"CMS simulation 2009");
      tex->Draw();
      tex = new TLatex(3.625803,6.712477,"per nb^{-1}");
   tex->SetLineWidth(2);
   tex->Draw();
//       tex = new TLatex(2.331331,10.70125,"global-global");
//    tex->SetTextColor(4);
//    tex->SetTextSize(0.035);
//    tex->SetLineWidth(2);
//    tex->Draw();
//       tex = new TLatex(2.368704,9.966102,"plus global-tracker");
//    tex->SetTextColor(4);
//    tex->SetTextSize(0.035);
//    tex->SetLineWidth(2);
//    tex->Draw();
      tex = new TLatex(2.399282,8.320604,"plus tracker-tracker");
   tex->SetTextSize(0.04);
   tex->SetTextColor(4);
   tex->SetLineWidth(2);
   tex->Draw();
//       tex = new TLatex(2.05273,4.08495,"global-global");
//    tex->SetTextSize(0.035);
//    tex->SetTextColor(2);
//    tex->SetLineWidth(2);
//    tex->Draw();
      tex = new TLatex(2.02555,3.946498,"plus global-tracker");
   tex->SetTextSize(0.04);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(2.035742,0.5930161,"global-global");
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();

    }
  }


  //plot the pT histograms:
  TCanvas *c1Pt;
  TLatex *texPt;
  TH1F *hFramePt;
  sprintf(name, "c1Pt");
  c1Pt = new TCanvas(name, "pt");
  hFramePt = gPad->DrawFrame(0., 0., 8., 1.2*hReco_pT[0][1]->GetMaximum());
  hFramePt->SetXTitle(hReco_pT[0][1]->GetXaxis()->GetTitle());
  hFramePt->SetYTitle(hReco_pT[0][1]->GetYaxis()->GetTitle());
//   hFramePt->SetYTitle("dN/dp_{T} [per 400 MeV/c]");
  TLegend *legPt = new TLegend(0.72,0.684322,0.9899425,0.9618644, enLabel);

  for(int iCat = 0; iCat < 3; iCat++){
    hReco_pT[0][iCat]->SetLineWidth(2);
    hReco_pT[0][iCat]->SetLineColor(colour[iCat]);
    hReco_pT[0][iCat]->Draw("same");
    legPt->AddEntry(hReco_pT[0][iCat], catNames[iCat], "l");
  }
  legPt->SetTextSize(0.04); legPt->SetFillColor(0);
  legPt->Draw();
  sprintf(name, "Figures/%s_pT.gif", fileNameIn);
  c1Pt->Print(name);
  sprintf(name, "Figures/%s_pT.pdf", fileNameIn);
  c1Pt->Print(name);

  //get one stacked histogram:
  THStack *hStackPT = new THStack("hStackPT", "");
  hReco_pT[0][0]->SetFillColor(1);
  hStackPT->Add(hReco_pT[0][0]);
  hReco_pT[0][1]->SetFillColor(2);
  hStackPT->Add(hReco_pT[0][1]);
  hReco_pT[0][2]->SetFillColor(4);
  hStackPT->Add(hReco_pT[0][2]);

  TCanvas *c1PTStack = new TCanvas("c1PTStack", "pT stacked");
  TH1F *hFramePT = gPad->DrawFrame(0., 0., 8., 1.1*hStackPT->GetMaximum());
  hFramePT->SetXTitle("p_{T} (J/#psi) [GeV/c]");
//   hFramePT->SetYTitle("dN/dp_{T}");
  hFramePT->SetYTitle(hReco_pT[0][1]->GetYaxis()->GetTitle());
  hStackPT->Draw("histsame");

  TLegend *legPTSt = new TLegend(0.72,0.684322,0.9899425,0.9618644, enLabel);
  for(int iCat = 2; iCat >= 0; iCat--)
    legPTSt->AddEntry(hReco_pT[0][iCat], catNames[iCat], "f");
  legPTSt->SetTextSize(0.04); legPTSt->SetFillColor(0);
  legPTSt->Draw();
  sprintf(name, "Figures/%s_pTStacked.gif", fileNameIn);
  c1PTStack->Print(name);
  sprintf(name, "Figures/%s_pTStacked.pdf", fileNameIn);
  c1PTStack->Print(name);


  //plot the rap histograms:
  TCanvas *c1Rap;
  TLatex *texRap;
  TH1F *hFrameRap;
  sprintf(name, "c1Rap");
  c1Rap = new TCanvas(name, "rap");
  hFrameRap = gPad->DrawFrame(-2.5, 0., 2.5, 1.2*hReco_rap[0][1]->GetMaximum());
  hFrameRap->SetXTitle(hReco_rap[0][1]->GetXaxis()->GetTitle());
  hFrameRap->SetYTitle(hReco_rap[0][1]->GetYaxis()->GetTitle());
//    hFrameRap->SetYTitle("dN/dy");
  TLegend *legRap = new TLegend(0.5474138,0.6567797,0.8175287,0.93432, enLabel);

  for(int iCat = 0; iCat < 3; iCat++){
    hReco_rap[0][iCat]->SetLineWidth(2);
    hReco_rap[0][iCat]->SetLineColor(colour[iCat]);
    hReco_rap[0][iCat]->Draw("same");
    legRap->AddEntry(hReco_rap[0][iCat], catNames[iCat], "l");
  }
  legRap->SetTextSize(0.04); legRap->SetFillColor(0);
  legRap->Draw();
  sprintf(name, "Figures/%s_rap.gif", fileNameIn);
  c1Rap->Print(name);
  sprintf(name, "Figures/%s_rap.pdf", fileNameIn);
  c1Rap->Print(name);

  //get one stacked histogram:
  THStack *hStack = new THStack("hStack", "");
  hReco_rap[0][0]->SetFillColor(1);
  hStack->Add(hReco_rap[0][0]);
  hReco_rap[0][1]->SetFillColor(2);
  hStack->Add(hReco_rap[0][1]);
  hReco_rap[0][2]->SetFillColor(4);
  hStack->Add(hReco_rap[0][2]);

  TCanvas *c1RapStack = new TCanvas("c1RapStack", "rap stacked");
  TH1F *hFrame = gPad->DrawFrame(-2.5, 0., 2.5, 1.1*hStack->GetMaximum());
  hFrame->SetXTitle("y(J/#psi)");
  hFrame->SetYTitle("dN/dy");
  hStack->Draw("histsame");

  TLegend *legRapSt = new TLegend(0.5474138,0.6567797,0.8175287,0.93432, enLabel);
  for(int iCat = 2; iCat >= 0; iCat--)
    legRapSt->AddEntry(hReco_rap[0][iCat], catNames[iCat], "f");
  legRapSt->SetTextSize(0.04); legRapSt->SetFillColor(0);
  legRapSt->Draw();
  sprintf(name, "Figures/%s_rapStacked.gif", fileNameIn);
  c1RapStack->Print(name);
  sprintf(name, "Figures/%s_rapStacked.pdf", fileNameIn);
  c1RapStack->Print(name);

}

//=======================
void LoadHistos(Char_t *fileNameIn, Bool_t takeBest){

  Char_t name[100], title[300];
  sprintf(name, "%s.root", fileNameIn);
  TFile *fIn = new TFile(name);
//   for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
  for(int iCharge = 0; iCharge < 1; iCharge++){
    for(int iCat = 0; iCat < kNbCath; iCat++){
      if(takeBest){
	//mass histos
	sprintf(name, "hRecoBest_mass_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
	hReco_mass[iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	hReco_mass[iCharge][iCat]->SetLineColor(colour[iCat]);
// 	hReco_mass[iCharge][iCat]->SetLineWidth(2);
	hReco_mass[iCharge][iCat]->SetAxisRange(2.,4.);
	//pT histos
	sprintf(name, "hRecoBest_pT_0_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
	hReco_pT[iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	//rap histos
	sprintf(name, "hRecoBest_rap_0_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
	hReco_rap[iCharge][iCat] = (TH1F *) gDirectory->Get(name);
      }
      else{
	//mass histos
	sprintf(name, "hReco_mass_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
	hReco_mass[iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	hReco_mass[iCharge][iCat]->SetLineColor(colour[iCat]);
	hReco_mass[iCharge][iCat]->SetAxisRange(2.,4.);
// 	hReco_mass[iCharge][iCat]->SetLineWidth(2);
	//pT histos
	sprintf(name, "hReco_pT_0_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
	hReco_pT[iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	//rap histos
	sprintf(name, "hReco_rap_0_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
	hReco_rap[iCharge][iCat] = (TH1F *) gDirectory->Get(name);
      }
    }
    //add the mass histos:
    sprintf(name, "hReco_mass_%s_GGandGT", chargeName[iCharge]);
    hReco_mass[iCharge][kNbCath] = (TH1F *) hReco_mass[iCharge][0]->Clone(name);
    hReco_mass[iCharge][kNbCath]->Add(hReco_mass[iCharge][1]);

    //add the mass histos:
    sprintf(name, "hReco_mass_%s_all", chargeName[iCharge]);
    hReco_mass[iCharge][kNbCath+1] = (TH1F *) hReco_mass[iCharge][0]->Clone(name);
    hReco_mass[iCharge][kNbCath+1]->Add(hReco_mass[iCharge][1]);
    hReco_mass[iCharge][kNbCath+1]->Add(hReco_mass[iCharge][2]);
  }
}

//=======================
void LoadBGHistos(Char_t *fileNameIn, Bool_t takeBest, Double_t lumi){

  Char_t name[100], title[300];
  TFile *fIn = new TFile(fileNameIn);
  Int_t rebinBy = 5;
//   for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
  for(int iCharge = 0; iCharge < 1; iCharge++){
    for(int iCat = 0; iCat < kNbCath; iCat++){

      if(takeBest){
	sprintf(name, "hRecoBest_mass_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
	hRecoBG_mass[iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	sprintf(name, "hRecoBestBG_mass_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
	hRecoBG_mass[iCharge][iCat]->SetName(name);
      }
      else{
	sprintf(name, "hReco_mass_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
	hRecoBG_mass[iCharge][iCat] = (TH1F *) gDirectory->Get(name);
	sprintf(name, "hRecoBG_mass_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
	hRecoBG_mass[iCharge][iCat]->SetName(name);
      }
      hRecoBG_mass[iCharge][iCat]->SetLineColor(colour[iCat]);
      hRecoBG_mass[iCharge][iCat]->SetLineWidth(2);
    }
    //add the mass histos:
    sprintf(name, "hRecoBG_mass_%s_GGandGT", chargeName[iCharge]);
    hRecoBG_mass[iCharge][kNbCath] = (TH1F *) hRecoBG_mass[iCharge][0]->Clone(name);
    hRecoBG_mass[iCharge][kNbCath]->Add(hRecoBG_mass[iCharge][1]);
    //add the mass histos:
    sprintf(name, "hRecoBG_mass_%s_all", chargeName[iCharge]);
    hRecoBG_mass[iCharge][kNbCath+1] = (TH1F *) hRecoBG_mass[iCharge][0]->Clone(name);
    hRecoBG_mass[iCharge][kNbCath+1]->Add(hRecoBG_mass[iCharge][1]);
    hRecoBG_mass[iCharge][kNbCath+1]->Add(hRecoBG_mass[iCharge][2]);
  }
  //save the yields in the J/psi mass window:
  for(int iCat = 0; iCat < kNbCath+2; iCat++){
    Int_t binMin = hRecoBG_mass[0][iCat]->GetXaxis()->FindBin(massMinJPsi);
    Int_t binMax = hRecoBG_mass[0][iCat]->GetXaxis()->FindBin(massMaxJPsi) - 1;
    printf("BG: taking integral for bins %d to %d\n", binMin, binMax);
    integralBG[iCat] = hRecoBG_mass[0][iCat]->Integral(binMin, binMax);
    errPossonianBG[iCat] = sqrt(integralBG[iCat]);
  }
  //normalise to 1 nb-1 and 10 MeV bin size
//   for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
   for(int iCharge = 0; iCharge < 1; iCharge++){
    for(int iCat = 0; iCat < kNbCath+2; iCat++){
      hRecoBG_mass[iCharge][iCat]->Rebin(rebinBy);
      //histos correspond to Lint = 5.5 nb-1
      //--> normalise to Lint = 1 nb-1
      hRecoBG_mass[iCharge][iCat]->Scale(1./(lumi*rebinBy)); 
      if(iCharge == 0){
	integralBG[iCat] /= lumi;
	errPossonianBG[iCat] /= lumi;
      }
    }
  }
}

//=======================
void AddSignalBG(){

  Char_t name[100], title[300];
//   for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
  for(int iCharge = 0; iCharge < 1; iCharge++){
    for(int iCat = 0; iCat < kNbCath; iCat++){

      sprintf(name, "hRecoSum_mass_%s_%s", chargeName[iCharge], oniaCatName[iCat]);
      hRecoSum_mass[iCharge][iCat] = (TH1F *) hReco_mass[iCharge][iCat]->Clone(name);
      hRecoSum_mass[iCharge][iCat]->Add(hRecoBG_mass[iCharge][iCat]);
    }
    //add now the statistics from the gl-gl and gl-tr
    hRecoSum_mass[iCharge][kNbCath] = (TH1F *) hRecoSum_mass[iCharge][0]->Clone("hRecoSum_mass_OS_glglAndgltr");
    hRecoSum_mass[iCharge][kNbCath]->Add(hRecoSum_mass[iCharge][1]);
  }
}
//========================================
void PrintStat(){

  printf("integral within %1.1f < M < %1.1f GeV\n", massMinJPsi, massMaxJPsi);
  for(int iCat = 0; iCat < kNbCath+2; iCat++){
    if(integralBG[iCat] > 0){
      sigOvBG[iCat] = integral[iCat] / integralBG[iCat];
      errSigOvBG[iCat] = sqrt(pow(errPossonian[iCat] / integral[iCat], 2) + 
			      pow(errPossonianBG[iCat] / integralBG[iCat], 2));
      errSigOvBG[iCat] *= sigOvBG[iCat];

      printf("%29s: signal = %1.0f +- %1.0f; BG = %1.0f +- %1.0f; --> S / BG = %1.0f +- %1.0f\n",
	     catNames[iCat], integral[iCat], errPossonian[iCat], 
	     integralBG[iCat], errPossonianBG[iCat],
	     sigOvBG[iCat], errSigOvBG[iCat]);
    }
  }
}

//========================================
Double_t gausexp(Double_t *x, Double_t *par) {
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2)) +
    par[3]*exp(par[4]*x[0]);
 
}

