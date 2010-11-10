#include <string.h>
#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"
#include "TPaveStats.h"
#include "TFile.h"
#include "TKey.h"


// Declaring functions
void ReadInHistos(Char_t*);
void PlotHistos(Int_t iRapBin, Int_t iPTBin);
void ReadInHistos(Char_t*);

// Global variables
Int_t const kNbMassBins = 2;              // Number of mass regions we consider (2 sidebands for now)
char massLabel[kNbMassBins] = {'L','R'};  // Name of the sidebands
Int_t const hLifetimeRebin = -1;          // -1 if we don't want to rebin (use positive value if bins are too little populated)
Double_t kMassMin = 2.5;
Double_t kMassMax = 4.1;


TH1F *hLifetime[kNbMassBins][jpsi::kNbPTMaxBins+1][jpsi::kNbRapForPTBins+1]; // Array of pointers to lifetime histograms



//====================================
void ReadInHistos(Char_t *fileNameIn){

   cout << "### EXECUTING: ReadInHistos" << endl;
   TFile *fIn = new TFile(fileNameIn);
   Char_t name[100];
   printf("reading file %s\n", fileNameIn);
   for(int iMass = 0; iMass < kNbMassBins; iMass++){
      for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
         for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
            // lifetime histograms
            sprintf(name, "Sideband_Lifetime_%c_pT%d_rap%d", massLabel[iMass], iPTBin, iRapBin);
            hLifetime[iMass][iPTBin][iRapBin] = (TH1F *) gDirectory->Get(name);
            printf("-> %c sideband, pT %d, rap %d has %f entries and integral %f\n", massLabel[iMass], iPTBin, iRapBin, hLifetime[iMass][iPTBin][iRapBin]->GetEntries(), hLifetime[iMass][iPTBin][iRapBin]->Integral() );
	    //hLifetime[iMass][iPTBin][iRapBin]->Rebin(2);
            
            // costheta and phi histograms (later)
         }
      }
   }
   cout << "### EXITING: ReadInHistos" << endl;
}



//====================================
void PlotHistos(Int_t iPTBin, Int_t iRapBin){
   Char_t name[256];
   
   Float_t max_y = hLifetime[0][iPTBin][iRapBin]->GetMaximum(); // gets maximum height in both histograms
   if (max_y < hLifetime[1][iPTBin][iRapBin]->GetMaximum()) 
      max_y = hLifetime[1][iPTBin][iRapBin]->GetMaximum();
   
   if(max_y < 1) return; // we don't want empty histograms

   sprintf(name, "c1_rap%d_pT%d", iRapBin, iPTBin);
   TCanvas *c1 = new TCanvas(name, name);
   c1->Draw();
   gPad->SetLogy();
   TH1F *hFrame = gPad->DrawFrame(kMassMin, 1, kMassMax, 0.9*max_y);
   hFrame->SetXTitle(hLifetime[0][iPTBin][iRapBin]->GetXaxis()->GetTitle());
   
   hLifetime[0][iPTBin][iRapBin]->SetMarkerColor(2);
   hLifetime[1][iPTBin][iRapBin]->SetMarkerColor(4);
   hLifetime[0][iPTBin][iRapBin]->SetMarkerStyle(20);
   hLifetime[1][iPTBin][iRapBin]->SetMarkerStyle(20);
   hLifetime[0][iPTBin][iRapBin]->SetMarkerSize(1.0);
   hLifetime[1][iPTBin][iRapBin]->SetMarkerSize(1.0);
   
   // if (hLifetimeRebin > 0) {
   //    hLifetime[0][iPTBin][iRapBin]->Rebin(hLifetimeRebin);
   //    hLifetime[1][iPTBin][iRapBin]->Rebin(hLifetimeRebin);
   // }
   
   hLifetime[0][iPTBin][iRapBin]->Draw("SAME"); //L
   hLifetime[1][iPTBin][iRapBin]->Draw("SAME"); //R

   //add some text
   if(iRapBin == 0 && iPTBin == 0) sprintf(name, "|y| < %1.1f, all p_{T}", jpsi::rapForPTRange[jpsi::kNbRapForPTBins]);
   else if(iRapBin == 0 && iPTBin > 0) sprintf(name, "|y| < %1.1f, %1.1f < p_{T} < %1.1f", jpsi::rapForPTRange[jpsi::kNbRapForPTBins], jpsi::pTRange[iRapBin][iPTBin-1], jpsi::pTRange[iRapBin][iPTBin]);
   else if(iRapBin > 0 && iPTBin == 0) sprintf(name, "%1.1f < |y| < %1.1f, all p_{T}", jpsi::rapForPTRange[iRapBin-1], jpsi::rapForPTRange[iRapBin]);
   
//    else if(iRapBin == 1 && iPTBin == 0) sprintf(name, "|y| < %1.1f, all p_{T}", jpsi::rapForPTRange[iRapBin]);
//    else if(iRapBin == 1) sprintf(name, "|y| < %1.1f, %1.1f < p_{T} < %1.1f", jpsi::rapForPTRange[iRapBin], jpsi::pTRange[iRapBin][iPTBin-1], jpsi::pTRange[iRapBin][iPTBin]);
//    else if(iRapBin > 1)  sprintf(name, "%1.1f < |y| < %1.1f, %1.1f < p_{T} < %1.1f", 
// 				jpsi::rapForPTRange[iRapBin-1], jpsi::rapForPTRange[iRapBin], jpsi::pTRange[iRapBin][iPTBin-1], jpsi::pTRange[iRapBin][iPTBin]);
   
   TLatex *tex = new TLatex(-0.8, 3*max_y, name);
   tex->SetTextSize(0.04); tex->Draw();

   TLegend *leg = new TLegend(0.6853448,0.7775424,0.9741379,0.940678);
   leg->AddEntry(hLifetime[0][iPTBin][iRapBin], "Left mass window", "p");
   leg->AddEntry(hLifetime[1][iPTBin][iRapBin], "Right mass window", "p");
   leg->SetTextSize(0.035);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->Draw();

   sprintf(name, "Figures/lifetimeBG_rap%d_pT%d.gif", iRapBin, iPTBin);
   c1->Print(name);
   sprintf(name, "Figures/lifetimeBG_rap%d_pT%d.pdf", iRapBin, iPTBin);
   c1->Print(name);
   // sprintf(name, "Figures/lifetimeBG_rap%d_pT%d.root", iRapBin, iPTBin);
   // c1->Print(name);
}



//====================================
void plotBGDist(Char_t *fileNameIn = "/afs/cern.ch/user/c/calligar/scratch0/data/pol_data_testoutput.root"){
   cout << "EXECUTING: plotBGDist" << endl;
   ReadInHistos(fileNameIn);
   for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
         PlotHistos(iRapBin, iPTBin);
      }
   }
}
