#include "rootIncludes.inc"
#include "commonVar.h"

Int_t const nBinsPT[2] = {2, 4}; //needs manual adjustment!
TGraphErrors *gFrac[kNbSpecies][kNbComp-1][kNbRap];
Float_t NEv[kNbSpecies][kNbRap][kNbPTMaxBins];

void FillGraphs();
void PlotGraphs(Int_t iRegion, Int_t iRap);
//=============================
void plotFractions(){

  FillGraphs();
  for(int iRap = 0; iRap < kNbRap; iRap++){
    PlotGraphs(NP, iRap);  PlotGraphs(P, iRap);
  }
}

//=============================
void PlotGraphs(Int_t iRegion, Int_t iRap){

  Char_t name[100];
  sprintf(name, "fractions_%s_rap%d", compName[iRegion], iRap);
  TCanvas *c1 = new TCanvas(name, name);
  TH1F *hFrame1 = gPad->DrawFrame(0., -0.1, 30.5, 1.1);
  hFrame1->SetXTitle("p_{T} [GeV]");
  hFrame1->SetYTitle("fraction");
  gFrac[iRegion][P][iRap]->Draw("p same");
  gFrac[iRegion][NP][iRap]->Draw("p same");
  gFrac[iRegion][BKG][iRap]->Draw("p same");

  //print the total statistics in [k]
  TLatex *tex1a = new TLatex(0., 0., "");
  tex1a->SetTextSize(0.04); tex1a->Draw();
  for(int iPT = 0; iPT < nBinsPT[iRap]; iPT++){
    sprintf(name, "%1.1f k", NEv[iRegion][iRap][iPT] / 1000.);
    Double_t pT, frac;
    gFrac[iRegion][P][iRap]->GetPoint(iPT, pT, frac);
    tex1a->DrawLatex(0.9*pT, 1.03, name);
  }

  if(iRap == 0)
    sprintf(name, "|y| < %1.1f", rapForPTRange[iRap+1]);
  else if(iRap == 1)
    sprintf(name, "%1.1f < |y| < %1.1f", rapForPTRange[iRap], rapForPTRange[iRap+1]);
  TLatex *tex1 = new TLatex(2., 0.7, name);
  tex1->SetTextSize(0.045), tex1->Draw();
  if(iRegion == P)
    sprintf(name, "l_{J/#psi} < 0.1 mm");
  else
    sprintf(name, "l_{J/#psi} > 0.1 mm");
  tex1->DrawLatex(2., 0.62, name);

  TLegend *leg1 = new TLegend(0.1551724,0.42,0.4554598,0.6);
  if(iRegion == P){
    leg1->AddEntry(gFrac[iRegion][P][iRap], "prompt J/#psi", "p");
    leg1->AddEntry(gFrac[iRegion][NP][iRap], "non-prompt J/#psi", "p");
    leg1->AddEntry(gFrac[iRegion][BKG][iRap], "background", "p");
  }
  else{
    leg1->AddEntry(gFrac[iRegion][NP][iRap], "non-prompt J/#psi", "p");
    leg1->AddEntry(gFrac[iRegion][BKG][iRap], "background", "p");
    leg1->AddEntry(gFrac[iRegion][P][iRap], "prompt J/#psi", "p");
  }
  leg1->SetTextSize(0.04); leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->Draw();

  TLine *line1 = new TLine(0., 0., 30.5, 0.);
  line1->SetLineStyle(3); line1->Draw();
  line1->DrawLine(0.,1.,30.5,1.);

  sprintf(name, "Figures/fractions_%sRegion_rap%d.pdf", compName[iRegion], iRap+1);
  c1->Print(name);
}

//=============================
void FillGraphs(){

  FILE *fIn = fopen("fractions.txt", "read");

  Char_t line[1000];
  Float_t fraction[kNbSpecies][kNbComp-1][kNbRap][kNbPTMaxBins];
  Float_t eFraction[kNbSpecies][kNbComp-1][kNbRap][kNbPTMaxBins];

  Float_t pT[kNbPTMaxBins], errPT[kNbPTMaxBins];
  for(int iRap = 0; iRap < kNbRap; iRap++){
    for(int iPT = 0; iPT < nBinsPT[iRap]; iPT++){
      fgets(line, sizeof(line), fIn); //comment

      //get the total number of events:
      fgets(line, sizeof(line), fIn);
      sscanf(line, "%*s %*s %*s %*f %*s %f", &NEv[0][iRap][iPT]);
      fgets(line, sizeof(line), fIn);
      sscanf(line, "%*s %*s %*s %*f %*s %f", &NEv[1][iRap][iPT]);
      
      //first the components in the P dominated region
      fgets(line, sizeof(line), fIn);
      sscanf(line, "%*s %*s %f %*s %f", &fraction[0][P][iRap][iPT], &eFraction[0][P][iRap][iPT]);
      fgets(line, sizeof(line), fIn);
      sscanf(line, "%*s %*s %f %*s %f", &fraction[0][NP][iRap][iPT], &eFraction[0][NP][iRap][iPT]);
      fgets(line, sizeof(line), fIn);
      sscanf(line, "%*s %*s %f %*s %f", &fraction[0][BKG][iRap][iPT], &eFraction[0][BKG][iRap][iPT]);
      //then in the NP dominated region
      fgets(line, sizeof(line), fIn);
      sscanf(line, "%*s %*s %f %*s %f", &fraction[1][P][iRap][iPT], &eFraction[1][P][iRap][iPT]);
      fgets(line, sizeof(line), fIn);
      sscanf(line, "%*s %*s %f %*s %f", &fraction[1][NP][iRap][iPT], &eFraction[1][NP][iRap][iPT]);
      fgets(line, sizeof(line), fIn);
      sscanf(line, "%*s %*s %f %*s %f", &fraction[1][BKG][iRap][iPT], &eFraction[1][BKG][iRap][iPT]);

      fgets(line, sizeof(line), fIn); //empty line

    }
    //set up the pT value and the error:
    if(iRap == 0){
      pT[0] = 17.5; pT[1] = 25;
      errPT[0] = 2.5; errPT[1] = 5.;
    }
    else if(iRap == 1){
      pT[0] = 9.; pT[1] = 12.5; pT[2] = 17.5; pT[3] = 25;
      errPT[0] = 1.; errPT[1] = 2.5; errPT[2] = 2.5; errPT[3] = 5.;
    }
    for(int iRegions = 0; iRegions < kNbSpecies; iRegions++){
      gFrac[iRegions][P][iRap] = new TGraphErrors(nBinsPT[iRap], pT, fraction[iRegions][P][iRap], errPT, eFraction[iRegions][P][iRap]);
      gFrac[iRegions][NP][iRap] = new TGraphErrors(nBinsPT[iRap], pT, fraction[iRegions][NP][iRap], errPT, eFraction[iRegions][NP][iRap]);
      gFrac[iRegions][BKG][iRap] = new TGraphErrors(nBinsPT[iRap], pT, fraction[iRegions][BKG][iRap], errPT, eFraction[iRegions][BKG][iRap]);
      gFrac[iRegions][P][iRap]->SetMarkerStyle(20); gFrac[iRegions][P][iRap]->SetMarkerColor(4); gFrac[iRegions][P][iRap]->SetLineColor(4); 
      gFrac[iRegions][NP][iRap]->SetMarkerStyle(25); gFrac[iRegions][NP][iRap]->SetMarkerColor(2); gFrac[iRegions][NP][iRap]->SetLineColor(2);
      gFrac[iRegions][BKG][iRap]->SetMarkerStyle(30); gFrac[iRegions][BKG][iRap]->SetMarkerColor(kMagenta+2); gFrac[iRegions][BKG][iRap]->SetLineColor(kMagenta+2);
    }
  }
  fclose(fIn);
}
