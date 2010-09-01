#include "/home/hermine/CMS/CMSSW/hWoehri/Polarization/interface/commonVar.h"
#include "/home/hermine/Work/rootIncludes.inc"

TGraphAsymmErrors *gLambdaTh1D_rap[kNbFrames][kNbRapForPTBins+1];
TGraphAsymmErrors *gLambdaPhi1D_rap[kNbFrames][kNbRapForPTBins+1];
TGraphAsymmErrors *gLambdaThPhi1D_rap[kNbFrames][kNbRapForPTBins+1];
TGraph *gFrameInv1D_rap[kNbFrames][kNbRapForPTBins+1];


void LoadGraphs(Char_t *label);
void FinalizeGraphsForPlotting();
void PlotGraphs(Int_t iFrame, Char_t *label);
void PlotAll(Char_t *label, Char_t *dataTypeLabel);
//=================================
void plotLambdas(Char_t *label = "pseudoData_HLT_Mu0Track0Jpsi_cut120_1Sep2010",
		 Char_t *dataTypeLabel = "MC"){

  LoadGraphs(label);
  FinalizeGraphsForPlotting(); //removes non-trustworthy bins (needs costumizing)

  PlotAll(label, dataTypeLabel);
  PlotGraphs(CS, label);
  PlotGraphs(HX, label);
}

//=================================
void PlotAll(Char_t *label, Char_t *dataTypeLabel){

  Char_t name[200];
  Int_t iFrame = 0; //CS
  gStyle->SetTitleSize(0.08, "y");
  gStyle->SetTitleOffset(0.6, "y");

  sprintf(name, "cLambdas");
  TCanvas *cLambdas = new TCanvas(name, name, 1000, 700);
  cLambdas->Divide(2,2);
  cLambdas->cd(1);
  //============================================
  //lambda_theta vs pT  (3 rap bins): CS
  //============================================
  TH1F *hFrameLTh1D = gPad->DrawFrame(0., -1.2, 30., 1.2);
  hFrameLTh1D->SetXTitle("p_{T} [GeV/c]");
  sprintf(name, "#lambda_{#theta}");
  hFrameLTh1D->SetYTitle(name);

  // gLambdaTh1D_rap[iFrame][1]->Draw("psame");
  gLambdaTh1D_rap[iFrame][2]->Draw("psame");
  gLambdaTh1D_rap[iFrame][3]->Draw("psame");
  
  TLine *line = new TLine(0., 0., 30., 0.);
  line->SetLineStyle(3);  line->Draw();
  line->DrawLine(0., 1., 30., 1.);
  line->DrawLine(0., -1., 30., -1.);

  TLatex *tex = new TLatex(2., 0.7, frameLabel[iFrame]);
  tex->SetTextSize(0.06); tex->Draw();
  tex->DrawLatex(24., -0.8, dataTypeLabel);

  cLambdas->cd(2);
  //============================================
  //lambda_phi vs pT (3 rap bins)
  //============================================
  TH1F *hFrameLPhi1D_pT = gPad->DrawFrame(0., -1.2, 30., 1.2);
  hFrameLPhi1D_pT->SetXTitle("p_{T} [GeV/c]");
  sprintf(name, "#lambda_{#phi}");
  hFrameLPhi1D_pT->SetYTitle(name);

  // gLambdaPhi1D_rap[iFrame][1]->Draw("psame");
  gLambdaPhi1D_rap[iFrame][2]->Draw("psame");
  gLambdaPhi1D_rap[iFrame][3]->Draw("psame");

  TLine *line2 = new TLine(0., 0., 30., 0.);
  line2->SetLineStyle(3);  line2->Draw();
  line2->DrawLine(0., 1., 30., 1.);
  line2->DrawLine(0., -1., 30., -1.);

  TLatex *tex2 = new TLatex(2., 0.7, frameLabel[iFrame]);
  tex2->SetTextSize(0.06); tex2->Draw();

  TLegend *leg = new TLegend(0.65,0.75,0.97,0.9);
  // sprintf(name, "|y| < %1.1f", rapForPTRange[1]);
  // leg->AddEntry(gLambdaTh1D_rap[iFrame][1], name, "p");
  sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[1], rapForPTRange[2]);
  leg->AddEntry(gLambdaTh1D_rap[iFrame][2], name, "p");
  sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[2], rapForPTRange[3]);
  leg->AddEntry(gLambdaTh1D_rap[iFrame][3], name, "p");
  // leg->SetBorderSize(0);
  leg->SetTextSize(0.05); leg->SetFillColor(0); leg->Draw();


  cLambdas->cd(3);
  iFrame = 1; //HX
  //============================================
  //lambda_theta vs pT  (3 rap bins): CS
  //============================================
  TH1F *hFrameLTh1D_2 = gPad->DrawFrame(0., -1.2, 30., 1.2);
  hFrameLTh1D_2->SetXTitle("p_{T} [GeV/c]");
  sprintf(name, "#lambda_{#theta}");
  hFrameLTh1D_2->SetYTitle(name);
  
  // gLambdaTh1D_rap[iFrame][1]->Draw("psame");
  gLambdaTh1D_rap[iFrame][2]->Draw("psame");
  gLambdaTh1D_rap[iFrame][3]->Draw("psame");
  
  TLine *line_2 = new TLine(0., 0., 30., 0.);
  line_2->SetLineStyle(3);  line_2->Draw();
  line_2->DrawLine(0., 1., 30., 1.);
  line_2->DrawLine(0., -1., 30., -1.);

  TLatex *tex3 = new TLatex(2., 0.7, frameLabel[iFrame]);
  tex3->SetTextSize(0.06); tex3->Draw();

  // TLegend *leg_2 = new TLegend(0.6609195,0.8,0.97,0.94);
  // // sprintf(name, "|y| < %1.1f", rapForPTRange[1]);
  // // leg_2->AddEntry(gLambdaTh1D_rap[iFrame][1], name, "p");
  // sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[1], rapForPTRange[2]);
  // leg_2->AddEntry(gLambdaTh1D_rap[iFrame][2], name, "p");
  // sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[2], rapForPTRange[3]);
  // leg_2->AddEntry(gLambdaTh1D_rap[iFrame][3], name, "p");
  // leg_2->SetBorderSize(0);
  // leg_2->SetTextSize(0.05); leg_2->SetFillColor(0); leg_2->Draw();

  cLambdas->cd(4);
  //============================================
  //lambda_phi vs pT (3 rap bins)
  //============================================
  sprintf(name, "cLambdaPhi1D_2_%s", frameLabel[iFrame]);
  TH1F *hFrameLPhi1D_pT_2 = gPad->DrawFrame(0., -1.2, 30., 1.2);
  hFrameLPhi1D_pT_2->SetXTitle("p_{T} [GeV/c]");
  sprintf(name, "#lambda_{#phi}");
  hFrameLPhi1D_pT_2->SetYTitle(name);

  // gLambdaPhi1D_rap[iFrame][1]->Draw("psame");
  gLambdaPhi1D_rap[iFrame][2]->Draw("psame");
  gLambdaPhi1D_rap[iFrame][3]->Draw("psame");

  TLine *line2_2 = new TLine(0., 0., 30., 0.);
  line2_2->SetLineStyle(3);  line2->Draw();
  line2_2->DrawLine(0., 1., 30., 1.);
  line2_2->DrawLine(0., -1., 30., -1.);

  TLatex *tex4 = new TLatex(2., 0.7, frameLabel[iFrame]);
  tex4->SetTextSize(0.06); tex4->Draw();

  sprintf(name, "Figures/lambdas_%s_%s.eps", frameLabel[iFrame], label); cLambdas->Print(name);
  sprintf(name, "Figures/lambdas_%s_%s.gif", frameLabel[iFrame], label); cLambdas->Print(name);
  sprintf(name, "Figures/lambdas_%s_%s.pdf", frameLabel[iFrame], label); cLambdas->Print(name);
}

//=================================
void PlotGraphs(Int_t iFrame, Char_t *label){

  Char_t name[200];
  //============================================
  //lambda_theta vs pT  (3 rap bins)
  //============================================
  sprintf(name, "cLambdaTh1D_%s", frameLabel[iFrame]);
  TCanvas *cLambdaTh1D = new TCanvas(name, name);
  TH1F *hFrameLTh1D = gPad->DrawFrame(0., -1.2, 30., 1.2);
  hFrameLTh1D->SetXTitle("p_{T} [GeV/c]");
  sprintf(name, "#lambda_{#theta}^{%s}", frameLabel[iFrame]);
  hFrameLTh1D->SetYTitle(name);
  
  // gLambdaTh1D_rap[iFrame][1]->Draw("psame");
  gLambdaTh1D_rap[iFrame][2]->Draw("psame");
  gLambdaTh1D_rap[iFrame][3]->Draw("psame");
  
  TLine *line = new TLine(0., 0., 30., 0.);
  line->SetLineStyle(3);  line->Draw();
  line->DrawLine(0., 1., 30., 1.);
  line->DrawLine(0., -1., 30., -1.);

  // TLegend *leg = new TLegend(0.6609195,0.8,0.97,0.94);
  // // sprintf(name, "|y| < %1.1f", rapForPTRange[1]);
  // // leg->AddEntry(gLambdaTh1D_rap[iFrame][1], name, "p");
  // sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[1], rapForPTRange[2]);
  // leg->AddEntry(gLambdaTh1D_rap[iFrame][2], name, "p");
  // sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[2], rapForPTRange[3]);
  // leg->AddEntry(gLambdaTh1D_rap[iFrame][3], name, "p");
  // leg->SetBorderSize(0);
  // leg->SetTextSize(0.05); leg->SetFillColor(0); leg->Draw();

  sprintf(name, "Figures/lambdaTh1D_%s_%s.eps", frameLabel[iFrame], label);  cLambdaTh1D->Print(name);
  sprintf(name, "Figures/lambdaTh1D_%s_%s.pdf", frameLabel[iFrame], label);  cLambdaTh1D->Print(name);
  sprintf(name, "Figures/lambdaTh1D_%s_%s.gif", frameLabel[iFrame], label);  cLambdaTh1D->Print(name);
  
  //============================================
  //lambda_phi vs pT (3 rap bins)
  //============================================
  sprintf(name, "cLambdaPhi1D_%s", frameLabel[iFrame]);
  TCanvas *cLambdaPhi1D_pT = new TCanvas(name, name);
  TH1F *hFrameLPhi1D_pT = gPad->DrawFrame(0., -1.2, 30., 1.2);
  hFrameLPhi1D_pT->SetXTitle("p_{T} [GeV/c]");
  sprintf(name, "#lambda_{#phi}^{%s}", frameLabel[iFrame]);
  hFrameLPhi1D_pT->SetYTitle(name);

  // gLambdaPhi1D_rap[iFrame][1]->Draw("psame");
  gLambdaPhi1D_rap[iFrame][2]->Draw("psame");
  gLambdaPhi1D_rap[iFrame][3]->Draw("psame");

  TLine *line2 = new TLine(0., 0., 30., 0.);
  line2->SetLineStyle(3);  line2->Draw();
  line2->DrawLine(0., 1., 30., 1.);
  line2->DrawLine(0., -1., 30., -1.);

  sprintf(name, "Figures/lambdaPhi1D_%s_%s.eps", frameLabel[iFrame], label);  cLambdaPhi1D_pT->Print(name);
  sprintf(name, "Figures/lambdaPhi1D_%s_%s.pdf", frameLabel[iFrame], label);  cLambdaPhi1D_pT->Print(name);
  sprintf(name, "Figures/lambdaPhi1D_%s_%s.gif", frameLabel[iFrame], label);  cLambdaPhi1D_pT->Print(name);

  // //============================================
  // //lambda_ThetaPhi vs pT (3 rap bins)
  // //============================================
  // sprintf(name, "cLambdaThPhi1D_%s", frameLabel[iFrame]);
  // TCanvas *cLambdaThPhi1D_pT = new TCanvas(name, name);
  // TH1F *hFrameLThPhi1D_pT = gPad->DrawFrame(0., -1.4, 30., 1.2);
  // hFrameLThPhi1D_pT->SetXTitle("p_{T} [GeV/c]");
  // sprintf(name, "#lambda_{#theta#phi}^{%s}", frameLabel[iFrame]);
  // hFrameLThPhi1D_pT->SetYTitle(name);

  // gLambdaThPhi1D_rap[iFrame][1]->Draw("psame");
  // gLambdaThPhi1D_rap[iFrame][2]->Draw("psame");
  // gLambdaThPhi1D_rap[iFrame][3]->Draw("psame");

  // TLine *line3 = new TLine(0., 0., 30., 0.);
  // line3->SetLineStyle(3);  line3->Draw();
  // line3->DrawLine(0., 1., 30., 1.);
  // line3->DrawLine(0., -1., 30., -1.);

  // sprintf(name, "Figures/lambdaThPhi1D_%s_%s.eps", frameLabel[iFrame], label);  cLambdaThPhi1D_pT->Print(name);
  // sprintf(name, "Figures/lambdaThPhi1D_%s_%s.pdf", frameLabel[iFrame], label);  cLambdaThPhi1D_pT->Print(name);
  // sprintf(name, "Figures/lambdaThPhi1D_%s_%s.gif", frameLabel[iFrame], label);  cLambdaThPhi1D_pT->Print(name);

  //============================================
  //frame invariant F
  //============================================
  sprintf(name, "cFrameInvF1D_%s", frameLabel[iFrame]);
  TCanvas *cFrameInv1D_pT = new TCanvas(name, name);
  TH1F *hFrameInvF1D_pT = gPad->DrawFrame(0., -0.2, 30., 1.2);
  hFrameInvF1D_pT->SetXTitle("p_{T} [GeV/c]");
  sprintf(name, "frame Invariant F for %s", frameLabel[iFrame]);
  hFrameInvF1D_pT->SetYTitle(name);

  // gFrameInv1D_rap[iFrame][1]->Draw("psame");
  gFrameInv1D_rap[iFrame][2]->Draw("psame");
  gFrameInv1D_rap[iFrame][3]->Draw("psame");

  TLine *line4 = new TLine(0., 0., 30., 0.);
  line4->SetLineStyle(3);  line4->Draw();
  line4->DrawLine(0., 1., 30., 1.);

  sprintf(name, "Figures/frameInvF1D_%s_%s.eps", frameLabel[iFrame], label);  cFrameInv1D_pT->Print(name);
  sprintf(name, "Figures/frameInvF1D_%s_%s.pdf", frameLabel[iFrame], label);  cFrameInv1D_pT->Print(name);
  sprintf(name, "Figures/frameInvF1D_%s_%s.gif", frameLabel[iFrame], label);  cFrameInv1D_pT->Print(name);
}



//=================================
void LoadGraphs(Char_t *label){

  TFile *fIn;
  Char_t name[200], fileName[200];
  sprintf(fileName, "Results/fitPar1D_lambda_%s.root", label);
  printf("reading TGraphs from file %s\n", fileName);
  fIn = new TFile(fileName);
  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iRap = 1; iRap < kNbRapForPTBins+1; iRap++){
      sprintf(name, "gLambdaTh1D_%s_rap%d", frameLabel[iFrame], iRap);
      gLambdaTh1D_rap[iFrame][iRap] = (TGraphAsymmErrors *) gDirectory->Get(name);
      sprintf(name, "gLambdaPhi1D_%s_rap%d", frameLabel[iFrame], iRap);
      gLambdaPhi1D_rap[iFrame][iRap] = (TGraphAsymmErrors *) gDirectory->Get(name);
      sprintf(name, "gLambdaThPhi1D_%s_rap%d", frameLabel[iFrame], iRap);
      gLambdaThPhi1D_rap[iFrame][iRap] = (TGraphAsymmErrors *) gDirectory->Get(name);
      sprintf(name, "gFrameInv1D_%s_rap%d", frameLabel[iFrame], iRap);
      gFrameInv1D_rap[iFrame][iRap] = (TGraph *) gDirectory->Get(name);
    }
  }
}

//=================================
void FinalizeGraphsForPlotting(){

  for(int iFrame = 0; iFrame < kNbFrames; iFrame++){
    for(int iRap = 1; iRap < kNbRapForPTBins+1; iRap++){
      gLambdaTh1D_rap[iFrame][iRap]->RemovePoint(0);
      gLambdaPhi1D_rap[iFrame][iRap]->RemovePoint(0);
      gLambdaThPhi1D_rap[iFrame][iRap]->RemovePoint(0);
      gFrameInv1D_rap[iFrame][iRap]->RemovePoint(0);
    }
  }
}
