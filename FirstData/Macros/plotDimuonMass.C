#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TStyle.h"
#include "TLatex.h"

#include "commonVar.inc"
Int_t colour[kNbCharge] = {1,2,3,2};
TH1F *hReco_mass[kNbCharge][kNbCath][kNbCuts];
Char_t *oniaLabel[3] = {"global-global", "global-tracker", "tracker-tracker"};

void LoadData(Char_t *fileNameIn);
void PlotData(Char_t *enLabel);
//=======================================
void plotDimuonMass(Char_t *fileNameIn = "histos_data09_2360GeV.root", 
		    Char_t *enLabel = "2360"){

  LoadData(fileNameIn);
  PlotData(enLabel);

}
//=======================================
void PlotData(Char_t *enLabel){

  gStyle->SetOptStat(0);

  TCanvas *c1[kNbCuts];
  TLatex *tex1[kNbCath][kNbCuts];
  Char_t name[100], title[100];
  for(int iCut = 0; iCut < kNbCuts; iCut++){
    sprintf(name, "c1_cut%d", iCut);
    c1[iCut] = new TCanvas(name, "", 1100, 400);
    c1[iCut]->Divide(3,1);
    for(int iCat = 0; iCat < 3; iCat++){
      c1[iCut]->cd(iCat+1);
      hReco_mass[0][iCat][iCut]->SetMaximum(2.*hReco_mass[0][iCat][iCut]->GetMaximum());
      hReco_mass[0][iCat][iCut]->Draw("hist");
      hReco_mass[3][iCat][iCut]->Draw("histsame");
      tex1[iCat][iCut] = new TLatex(3.3, 0.8*hReco_mass[0][iCat][iCut]->GetMaximum(), oniaLabel[iCat]);
      tex1[iCat][iCut]->SetTextSize(0.05);
      tex1[iCat][iCut]->Draw();
      sprintf(name, "%s GeV", enLabel);
      tex1[iCat][iCut]->DrawLatex(3.3, 0.87*hReco_mass[0][iCat][iCut]->GetMaximum(), name);
      tex1[iCat][iCut]->DrawLatex(3.3, 0.7*hReco_mass[0][iCat][iCut]->GetMaximum(), muCutName[iCut]);
    }
    sprintf(name, "Figures/dimuMass_%sGeV_%s.gif", muCutName[iCut], enLabel);
    c1[iCut]->Print(name);
  }
}
//=======================================
void LoadData(Char_t *fileNameIn){

  TFile *fIn = new TFile(fileNameIn, "r");
  Char_t name[100];
  for(int iCharge = 0; iCharge < kNbCharge; iCharge++){
    for(int iCat = 0; iCat < 3; iCat++){
      for(int iCut = 0; iCut < kNbCuts; iCut++){
	sprintf(name, "hReco_mass_%s_%s_%s", chargeName[iCharge], oniaCatName[iCat], muCutName[iCut]);
	hReco_mass[iCharge][iCat][iCut] = (TH1F *) gDirectory->Get(name);
	hReco_mass[iCharge][iCat][iCut]->Rebin(5);
	hReco_mass[iCharge][iCat][iCut]->SetYTitle("dN/dM [per 50 MeV/c^{2}]");
	hReco_mass[iCharge][iCat][iCut]->SetAxisRange(1., 5.);
	hReco_mass[iCharge][iCat][iCut]->SetLineColor(colour[iCharge]);
      }
    }
  }
}
