#include "rootIncludes.inc"
#include "commonVar.h"

Int_t myColour[kNbRap] = {2,4};
Int_t myMarkerStyle[kNbRap] = {20, 25};
Int_t myStars[kNbFrames] = {29, 30};

TGraphAsymmErrors *gLTh[kNbSpecies][kNbFrames][kNbRap];
TGraphAsymmErrors *gLPhi[kNbSpecies][kNbFrames][kNbRap];
TGraphAsymmErrors *gLThPhi[kNbSpecies][kNbFrames][kNbRap];
TGraphAsymmErrors *gLTilde[kNbSpecies][kNbFrames][kNbRap];

TGraphAsymmErrors *gLTh_CDF_P;

void ReadGraphs(Int_t iFrame, Int_t iRap);
void PlotGraphsPandNP(Int_t iFrame, Int_t iRap);
void PlotGraphsRap(Int_t iFrame, Int_t iSpecies);
void PlotLambdaTilde(Int_t iRap, Int_t iSpecies);
void LoadCDF();
void PlotTogetherWithCDF();
void PrintPars();
//====================================
void plotPolPar(){

  for(int iRap = 0; iRap < kNbRap; iRap++){
    ReadGraphs(HX, iRap);
    ReadGraphs(CS, iRap);
  }
  //remove a faulty bin (might not be needed in the future...)
  //|y| < 0.9, 10-15 GeV/c
  for(int iSpecies = 0; iSpecies < kNbSpecies; iSpecies++){
    for(int iFrame = CS; iFrame <= CS; iFrame++){
      gLTh[iSpecies][iFrame][0]->RemovePoint(0);
      gLPhi[iSpecies][iFrame][0]->RemovePoint(0);
      gLThPhi[iSpecies][iFrame][0]->RemovePoint(0);
      gLTilde[iSpecies][iFrame][0]->RemovePoint(0);
    }
  }
  //remove many more points:
  for(int iSpecies = 0; iSpecies < kNbSpecies; iSpecies++){
    for(int iFrame = CS; iFrame < kNbFrames; iFrame++){
      gLTh[iSpecies][iFrame][1]->RemovePoint(0);
      gLPhi[iSpecies][iFrame][1]->RemovePoint(0);
      gLThPhi[iSpecies][iFrame][1]->RemovePoint(0);
      gLTilde[iSpecies][iFrame][1]->RemovePoint(0);

//       gLTh[iSpecies][iFrame][1]->RemovePoint(0);
//       gLPhi[iSpecies][iFrame][1]->RemovePoint(0);
//       gLThPhi[iSpecies][iFrame][1]->RemovePoint(0);
//       gLTilde[iSpecies][iFrame][1]->RemovePoint(0);

      gLTh[iSpecies][iFrame][1]->RemovePoint(2);
      gLPhi[iSpecies][iFrame][1]->RemovePoint(2);
      gLThPhi[iSpecies][iFrame][1]->RemovePoint(2);
      gLTilde[iSpecies][iFrame][1]->RemovePoint(2);
    }
  }


  PrintPars();

  LoadCDF();

//   for(int iRap = 1; iRap <= kNbRap; iRap++){
//     PlotGraphsPandNP(CS, iRap);
//     PlotGraphsPandNP(HX, iRap);
//   }

  PlotGraphsRap(CS, P);
  PlotGraphsRap(HX, P);
  PlotGraphsRap(CS, NP);
  PlotGraphsRap(HX, NP);


  for(int iRap = 0; iRap < kNbRap; iRap++){
    PlotLambdaTilde(iRap, P);
    PlotLambdaTilde(iRap, NP);
  }
  
  PlotTogetherWithCDF();
}

//====================================
void PlotTogetherWithCDF(){

  Char_t name[100];
  gStyle->SetTitleOffset(0.9, "y");

  TCanvas *cCDF = new TCanvas("cCDF");
  TH1F *hFrame4 = gPad->DrawFrame(0., -1., 30.5, 1.);
  hFrame4->SetXTitle("p_{T} [GeV]");
  hFrame4->GetYaxis()->SetTitleSize(0.06);
  sprintf(name, "#lambda_{#theta}^{HX}");
  hFrame4->SetYTitle(name);

  gLTh[P][HX][0]->Draw("p same");
  gLTh_CDF_P->Draw("p same");

  TLine *line4 = new TLine(0., 0., 30.5, 0.);
  line4->SetLineStyle(3);  line4->Draw();

  TLegend *leg5 = new TLegend(0.1408046,0.1334746,0.33,0.28);
  leg5->AddEntry(gLTh[P][HX][0], "CMS Preliminary,  #sqrt{s} = 7 TeV,  L = 15.8 pb^{-1},  |y| < 0.9", "p");
  leg5->AddEntry(gLTh_CDF_P, "CDF,  #sqrt{s} = 1.96 TeV,  L = 800 pb^{-1},  |y| < 0.6", "p");
  leg5->SetFillColor(0);  leg5->SetTextSize(0.04);
  leg5->SetBorderSize(0);
  leg5->Draw();

//   TLatex *tex5 = new TLatex(21.6, -0.76, "CMS preliminary");
//   tex5->SetTextSize(0.045);  tex5->Draw();
//   tex5->DrawLatex(18.6, -0.9, "#sqrt{s} = 7 TeV, L = 15.8 pb^{-1}");

  sprintf(name, "Figures/lambdaTheta_CDF.gif");
  cCDF->Print(name);
  sprintf(name, "Figures/lambdaTheta_CDF.pdf");
  cCDF->Print(name);

}

//====================================
void PlotLambdaTilde(Int_t iRap, Int_t iSpecies){

  Char_t name[100];
  //===============================================
  //lambda_thTilde
  //===============================================
  sprintf(name, "LambdaThTilde_%s_%d", speciesLabel[iSpecies], iRap);
  TCanvas *c4 = new TCanvas(name, name);
  TH1F *hFrame4 = gPad->DrawFrame(0., -1., 30.5, 1.);
  hFrame4->SetXTitle("p_{T} [GeV]");
  hFrame4->GetYaxis()->SetTitleSize(0.06);
  sprintf(name, "#tilde{#lambda}");
  hFrame4->SetYTitle(name);

  gLTilde[iSpecies][CS][iRap]->Draw("p same");
  gLTilde[iSpecies][HX][iRap]->Draw("p same");

  TLine *line4 = new TLine(0., 0., 30.5, 0.);
  line4->SetLineStyle(3);  line4->Draw();

  if(iSpecies == 0)
    sprintf(name, "prompt J/#psi, %s", rapLabel[iRap]);
  else if(iSpecies == 1)
    sprintf(name, "non-prompt J/#psi, %s", rapLabel[iRap]);
  TLegend *leg4 = new TLegend(0.1408046,0.1334746,0.375,0.3241525, name);
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++)
    leg4->AddEntry(gLTilde[iSpecies][iFrame][iRap], frameLabel[iFrame], "p");
  leg4->SetFillColor(0);  leg4->SetTextSize(0.04);
  leg4->SetBorderSize(0);
  leg4->Draw();

  TLatex *tex4 = new TLatex(21.6, -0.76, "CMS preliminary");
  tex4->SetTextSize(0.045);  tex4->Draw();
  tex4->DrawLatex(18.6, -0.9, "#sqrt{s} = 7 TeV, L = 15.8 pb^{-1}");

  sprintf(name, "Figures/lambdaTilde_%s_rap%d.gif", speciesLabel[iSpecies], iRap+1);
  c4->Print(name);
  sprintf(name, "Figures/lambdaTilde_%s_rap%d.pdf", speciesLabel[iSpecies], iRap+1);
  c4->Print(name);

}

//====================================
void PlotGraphsRap(Int_t iFrame, Int_t iSpecies){

  gStyle->SetTitleOffset(0.9, "y");
  Char_t name[100];

  //===============================================
  //lambda_theta
  //===============================================
  sprintf(name, "LambdaTheta_%s_%s", frameLabel[iFrame], speciesLabel[iSpecies]);
  TCanvas *c1 = new TCanvas(name, name);
  TH1F *hFrame1 = gPad->DrawFrame(0., -1., 30.5, 1.);
  hFrame1->SetXTitle("p_{T} [GeV]");
  hFrame1->GetYaxis()->SetTitleSize(0.06);
  sprintf(name, "#lambda_{#theta}^{%s}", frameLabel[iFrame]);
  hFrame1->SetYTitle(name);


  for(int iRap = 0; iRap < kNbRap; iRap++)
    gLTh[iSpecies][iFrame][iRap]->Draw("p same");

  TLine *line1 = new TLine(0., 0., 30.5, 0.);
  line1->SetLineStyle(3);
  line1->Draw();

  if(iSpecies == 0)
    sprintf(name, "prompt J/#psi");
  else if(iSpecies == 1)
    sprintf(name, "non-prompt J/#psi");
  TLegend *leg1 = new TLegend(0.1408046,0.1334746,0.375,0.3241525, name);
  for(int iRap = 0; iRap < kNbRap; iRap++)
    leg1->AddEntry(gLTh[iSpecies][iFrame][iRap], rapLabel[iRap], "p");
  leg1->SetFillColor(0);  leg1->SetTextSize(0.04);
  leg1->SetBorderSize(0);
  leg1->Draw();

  TLatex *tex1 = new TLatex(21.6, 0.8, "CMS preliminary");
  tex1->SetTextSize(0.045);  tex1->Draw();
  tex1->DrawLatex(18.6, 0.62, "#sqrt{s} = 7 TeV, L = 15.8 pb^{-1}");

  sprintf(name, "Figures/lambdaTh_%s_%s.gif", speciesLabel[iSpecies], frameLabel[iFrame]);
  c1->Print(name);
  sprintf(name, "Figures/lambdaTh_%s_%s.pdf", speciesLabel[iSpecies], frameLabel[iFrame]);
  c1->Print(name);

  //===============================================
  //lambda_phi
  //===============================================
  sprintf(name, "LambdaPhi_%s_%s", frameLabel[iFrame], speciesLabel[iSpecies]);
  TCanvas *c2 = new TCanvas(name, name);
  TH1F *hFrame2 = gPad->DrawFrame(0., -1., 30.5, 1.);
  hFrame2->SetXTitle("p_{T} [GeV]");
  hFrame2->GetYaxis()->SetTitleSize(0.06);
  sprintf(name, "#lambda_{#phi}^{%s}", frameLabel[iFrame]);
  hFrame2->SetYTitle(name);

  for(int iRap = 0; iRap < kNbRap; iRap++)
    gLPhi[iSpecies][iFrame][iRap]->Draw("p same");

  TLine *line2 = new TLine(0., 0., 30.5, 0.);
  line2->SetLineStyle(3);
  line2->Draw();

  if(iSpecies == 0)
    sprintf(name, "prompt J/#psi");
  else if(iSpecies == 1)
    sprintf(name, "non-prompt J/#psi");
  TLegend *leg2 = new TLegend(0.1408046,0.1334746,0.375,0.3241525, name);
  for(int iRap = 0; iRap < kNbRap; iRap++)
    leg2->AddEntry(gLPhi[iSpecies][iFrame][iRap], rapLabel[iRap], "p");
  leg2->SetFillColor(0);  leg2->SetTextSize(0.04);
  leg2->SetBorderSize(0);
  leg2->Draw();

  TLatex *tex2 = new TLatex(21.6, 0.8, "CMS preliminary");
  tex2->SetTextSize(0.045);  tex2->Draw();
  tex2->DrawLatex(18.6, 0.62, "#sqrt{s} = 7 TeV, L = 15.8 pb^{-1}");

  sprintf(name, "Figures/lambdaPhi_%s_%s.gif", speciesLabel[iSpecies], frameLabel[iFrame]);
  c2->Print(name);
  sprintf(name, "Figures/lambdaPhi_%s_%s.pdf", speciesLabel[iSpecies], frameLabel[iFrame]);
  c2->Print(name);


  //===============================================
  //lambda_thetaPhi
  //===============================================
  sprintf(name, "LambdaThetaPhi_%s_%s", frameLabel[iFrame], speciesLabel[iSpecies]);
  TCanvas *c3 = new TCanvas(name, name);
  TH1F *hFrame3 = gPad->DrawFrame(0., -1./sqrt(2.), 30.5, 1./sqrt(2.));
  hFrame3->SetXTitle("p_{T} [GeV]");
  hFrame3->GetYaxis()->SetTitleSize(0.06);
  sprintf(name, "#lambda_{#theta#phi}^{%s}", frameLabel[iFrame]);
  hFrame3->SetYTitle(name);

  for(int iRap = 0; iRap < kNbRap; iRap++)
    gLThPhi[iSpecies][iFrame][iRap]->Draw("p same");

  TLine *line3 = new TLine(0., 0., 30.5, 0.);
  line3->SetLineStyle(3);
  line3->Draw();

  if(iSpecies == 0)
    sprintf(name, "prompt J/#psi");
  else if(iSpecies == 1)
    sprintf(name, "non-prompt J/#psi");
  TLegend *leg3 = new TLegend(0.1408046,0.1334746,0.375,0.3241525, name);
  for(int iRap = 0; iRap < kNbRap; iRap++)
    leg3->AddEntry(gLThPhi[iSpecies][iFrame][iRap], rapLabel[iRap], "p");
  leg3->SetFillColor(0);  leg3->SetTextSize(0.04);
  leg3->SetBorderSize(0);
  leg3->Draw();

  TLatex *tex3 = new TLatex(21.6, 0.55, "CMS preliminary");
  tex3->SetTextSize(0.045);  tex3->Draw();
  tex3->DrawLatex(18.6, 0.42, "#sqrt{s} = 7 TeV, L = 15.8 pb^{-1}");

  sprintf(name, "Figures/lambdaThPhi_%s_%s.gif", speciesLabel[iSpecies], frameLabel[iFrame]);
  c3->Print(name);
  sprintf(name, "Figures/lambdaThPhi_%s_%s.pdf", speciesLabel[iSpecies], frameLabel[iFrame]);
  c3->Print(name);

}


//====================================
void PlotGraphsPandNP(Int_t iFrame, Int_t iRap){

  gStyle->SetTitleOffset(0.9, "y");
  Char_t name[100];

  //===============================================
  //lambda_theta
  //===============================================
  sprintf(name, "c1_%s_rap%d", frameLabel[iFrame], iRap);
  TCanvas *c1 = new TCanvas(name, name);
  TH1F *hFrame1 = gPad->DrawFrame(0., -1., 30.5, 1.);
  hFrame1->SetXTitle("p_{T} [GeV]");
  hFrame1->GetYaxis()->SetTitleSize(0.06);
  sprintf(name, "#lambda_{#theta}^{%s}", frameLabel[iFrame]);
  hFrame1->SetYTitle(name);


  gLTh[P][iFrame][iRap]->Draw("p same");
  gLTh[NP][iFrame][iRap]->Draw("p same");

  TLine *line1 = new TLine(0., 0., 30.5, 0.);
  line1->SetLineStyle(3);
  line1->Draw();

  if(iRap == 1)
    sprintf(name, "  |y| < 0.9");
  else if(iRap == 2)
    sprintf(name, "  0.9 < |y| < 1.2");
  TLegend *leg1 = new TLegend(0.1408046,0.1334746,0.375,0.3241525, name);
  leg1->AddEntry(gLTh[P][iFrame][iRap], "prompt J/#psi", "p");
  leg1->AddEntry(gLTh[NP][iFrame][iRap], "non-prompt J/#psi", "p");
  leg1->SetFillColor(0);  leg1->SetTextSize(0.04);
  leg1->SetBorderSize(0);
  leg1->Draw();

  TLatex *tex1 = new TLatex(21.6, 0.8, "CMS preliminary");
  tex1->SetTextSize(0.045);  tex1->Draw();
  tex1->DrawLatex(18.6, 0.62, "#sqrt{s} = 7 TeV, L = 15.8 pb^{-1}");

  sprintf(name, "Figures/lambdaTh_PandNP_%s_rap%d.gif", frameLabel[iFrame], iRap);
  c1->Print(name);
  sprintf(name, "Figures/lambdaTh_PandNP_%s_rap%d.pdf", frameLabel[iFrame], iRap);
  c1->Print(name);

  //===============================================
  //lambda_theta_tilde
  //===============================================
  sprintf(name, "c2_%s_rap%d", frameLabel[iFrame], iRap);
  TCanvas *c2 = new TCanvas(name, name);
  TH1F *hFrame2 = gPad->DrawFrame(0., -1., 30.5, 1.);
  hFrame2->SetXTitle("p_{T} [GeV]");
  hFrame2->GetYaxis()->SetTitleSize(0.06);
  sprintf(name, "#tilde{#lambda}^{ %s}", frameLabel[iFrame]);
  hFrame2->SetYTitle(name);


  gLTilde[P][iFrame][iRap]->Draw("p same");
  gLTilde[NP][iFrame][iRap]->Draw("p same");

  TLine *line2 = new TLine(0., 0., 30.5, 0.);
  line2->SetLineStyle(3);
  line2->Draw();

  if(iRap == 1)
    sprintf(name, "  |y| < 0.9");
  else if(iRap == 2)
    sprintf(name, "  0.9 < |y| < 1.2");
  TLegend *leg2 = new TLegend(0.1408046,0.1334746,0.375,0.3241525, name);
  leg2->AddEntry(gLTilde[P][iFrame][iRap], "prompt J/#psi", "p");
  leg2->AddEntry(gLTilde[NP][iFrame][iRap], "non-prompt J/#psi", "p");
  leg2->SetFillColor(0);  leg2->SetTextSize(0.04);
  leg2->SetBorderSize(0);
  leg2->Draw();

  TLatex *tex2 = new TLatex(21.6, -0.6, "CMS preliminary");
  tex2->SetTextSize(0.045);  tex2->Draw();
  tex2->DrawLatex(18.6, -0.78, "#sqrt{s} = 7 TeV, L = 15.8 pb^{-1}");

  sprintf(name, "Figures/lambdaThTilde_PandNP_rap%d.gif", iRap);
  c2->Print(name);
  sprintf(name, "Figures/lambdaThTilde_PandNP_rap%d.pdf", iRap);
  c2->Print(name);
}

//====================================
void ReadGraphs(Int_t iFrame, Int_t iRap){

  Char_t name[100];
  sprintf(name, "/home/hermine/CMS/Work/Polarization/PlotsForPaper/PolVsPt%s_rap%d.root", frameLabel[iFrame], iRap+1);
  TFile *fIn = new TFile(name);
  //1.) get the parameters for P J/psi:
  //a.) lambda_theta
  gLTh[P][iFrame][iRap] = (TGraphAsymmErrors *) gDirectory->Get("lth_pt_p");
  sprintf(name, "lth_P_%s_rap%d", frameLabel[iFrame], iRap);
  gLTh[P][iFrame][iRap]->SetName(name);
  gLTh[P][iFrame][iRap]->SetMarkerStyle(myMarkerStyle[iRap]);
  gLTh[P][iFrame][iRap]->SetMarkerColor(myColour[iRap]);
  gLTh[P][iFrame][iRap]->SetMarkerSize(1.2);
  //b.) lambda_phi
  gLPhi[P][iFrame][iRap] = (TGraphAsymmErrors *) gDirectory->Get("lphi_pt_p");
  sprintf(name, "lphi_P_%s_rap%d", frameLabel[iFrame], iRap);
  gLPhi[P][iFrame][iRap]->SetName(name);
  gLPhi[P][iFrame][iRap]->SetMarkerStyle(myMarkerStyle[iRap]);
  gLPhi[P][iFrame][iRap]->SetMarkerColor(myColour[iRap]);
  gLTh[P][iFrame][iRap]->SetMarkerSize(1.2);
  //c.) lambda_thetaPhi
  gLThPhi[P][iFrame][iRap] = (TGraphAsymmErrors *) gDirectory->Get("lthphi_pt_p");
  sprintf(name, "lthphi_P_%s_rap%d", frameLabel[iFrame], iRap);
  gLThPhi[P][iFrame][iRap]->SetName(name);
  gLThPhi[P][iFrame][iRap]->SetMarkerStyle(myMarkerStyle[iRap]);
  gLThPhi[P][iFrame][iRap]->SetMarkerColor(myColour[iRap]);
  //d.) lambda_tilde
  gLTilde[P][iFrame][iRap] = (TGraphAsymmErrors *) gDirectory->Get("lthtilde_pt_p");
  sprintf(name, "lthtilde_P_%s_rap%d", frameLabel[iFrame], iRap);
  gLTilde[P][iFrame][iRap]->SetName(name);
  gLTilde[P][iFrame][iRap]->SetMarkerStyle(myStars[iFrame]);
  gLTilde[P][iFrame][iRap]->SetMarkerSize(1.4);

  //2.) get the parameters for NP J/psi:
  //a.) lambda_theta
  gLTh[NP][iFrame][iRap] = (TGraphAsymmErrors *) gDirectory->Get("lth_pt_np");
  sprintf(name, "lth_NP_%s_rap%d", frameLabel[iFrame], iRap);
  gLTh[NP][iFrame][iRap]->SetName(name);
  gLTh[NP][iFrame][iRap]->SetMarkerStyle(myMarkerStyle[iRap]);
  gLTh[NP][iFrame][iRap]->SetMarkerColor(myColour[iRap]);
  //b.) lambda_phi
  gLPhi[NP][iFrame][iRap] = (TGraphAsymmErrors *) gDirectory->Get("lphi_pt_np");
  sprintf(name, "lphi_NP_%s_rap%d", frameLabel[iFrame], iRap);
  gLPhi[NP][iFrame][iRap]->SetName(name);
  gLPhi[NP][iFrame][iRap]->SetMarkerStyle(myMarkerStyle[iRap]);
  gLPhi[NP][iFrame][iRap]->SetMarkerColor(myColour[iRap]);
  //c.) lambda_thetaPhi
  gLThPhi[NP][iFrame][iRap] = (TGraphAsymmErrors *) gDirectory->Get("lthphi_pt_np");
  sprintf(name, "lthphi_NP_%s_rap%d", frameLabel[iFrame], iRap);
  gLThPhi[NP][iFrame][iRap]->SetName(name);
  gLThPhi[NP][iFrame][iRap]->SetMarkerStyle(myMarkerStyle[iRap]);
  gLThPhi[NP][iFrame][iRap]->SetMarkerColor(myColour[iRap]);
  //d.) lambda_tilde
  gLTilde[NP][iFrame][iRap] = (TGraphAsymmErrors *) gDirectory->Get("lthtilde_pt_np");
  sprintf(name, "lthtilde_NP_%s_rap%d", frameLabel[iFrame], iRap);
  gLTilde[NP][iFrame][iRap]->SetName(name);
  gLTilde[NP][iFrame][iRap]->SetMarkerStyle(myStars[iFrame]);
  gLTilde[NP][iFrame][iRap]->SetMarkerSize(1.4);

  fIn->Close();
}

//==============================
void LoadCDF(){

 double  pT_CDF[6]  = {  5.5,  6.5,  7.8,  10.1,  13.7,  20.0  };
 double D1pT_CDF[6] = {  0.5,  0.5,  0.8,   1.1,   1.7,   3.0  };
 double D2pT_CDF[6] = {  0.5,  0.5,  1.2,   1.9,   3.3,  10.0  };

 double lth_CDF[6] = {-0.004, -0.015, -0.077, -0.094, -0.140, -0.187};
 double Dlth_CDF[6] = {0.0304, 0.0297, 0.0264, 0.0289, 0.0436, 0.0907 };

 gLTh_CDF_P = new TGraphAsymmErrors(6, pT_CDF, lth_CDF, D1pT_CDF, D2pT_CDF, Dlth_CDF, Dlth_CDF );
 gLTh_CDF_P->SetMarkerStyle(24);

}

//=============================
void PrintPars(){
  
  for(int iSpecies = 0; iSpecies < kNbSpecies; iSpecies++){
    printf("%s\n", speciesLabel[iSpecies]);
    for(int iFrame = kNbFrames-1; iFrame >= 0; iFrame--){
      printf("%s\n", frameLabel[iFrame]);
      for(int iRap = 0; iRap < kNbRap; iRap++){

	Int_t nPoints = gLTh[iSpecies][iFrame][iRap]->GetN();
	Double_t x, y[4];
	Double_t errXlow, errXhigh;
	Double_t errYlow[4], errYhigh[4];
	for(int iP = 0; iP < nPoints; iP++){

	  //lambda_theta
	  gLTh[iSpecies][iFrame][iRap]->GetPoint(iP,x,y[0]);
	  errXlow = gLTh[iSpecies][iFrame][iRap]->GetErrorXlow(iP);
	  errXhigh = gLTh[iSpecies][iFrame][iRap]->GetErrorXhigh(iP);
	  errYlow[0] = gLTh[iSpecies][iFrame][iRap]->GetErrorYlow(iP);
	  errYhigh[0] = gLTh[iSpecies][iFrame][iRap]->GetErrorYhigh(iP);
	  //lambda_phi
	  gLPhi[iSpecies][iFrame][iRap]->GetPoint(iP,x,y[1]);
	  errYlow[1] = gLPhi[iSpecies][iFrame][iRap]->GetErrorYlow(iP);
	  errYhigh[1] = gLPhi[iSpecies][iFrame][iRap]->GetErrorYhigh(iP);
	  //lambda_thetaPhi
	  gLThPhi[iSpecies][iFrame][iRap]->GetPoint(iP,x,y[2]);
	  errYlow[2] = gLThPhi[iSpecies][iFrame][iRap]->GetErrorYlow(iP);
	  errYhigh[2] = gLThPhi[iSpecies][iFrame][iRap]->GetErrorYhigh(iP);

	  printf("%1.1f--%1.1f & %1.0f--%1.0f & $%1.1f^{+%1.1f}_{-%1.1f}$ & $%1.2f^{+%1.2f}_{-%1.2f}$ & $%1.2f^{+%1.2f}_{-%1.2f}$ & $%1.2f^{+%1.2f}_{-%1.2f}$\n",
		 rapForPTRange[iRap], rapForPTRange[iRap+1], 
		 x-errXlow, x+errXhigh, //pTmin-pTmax
		 x, errXhigh, errXlow,  //pTmean and +-
		 y[0], errYhigh[0], errYlow[0], 
		 y[1], errYhigh[1], errYlow[1], 
		 y[2], errYhigh[2], errYlow[2]);
	}
      }
    }
  }
}
