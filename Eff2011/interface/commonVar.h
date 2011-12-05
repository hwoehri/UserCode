#include "TLorentzVector.h"
#include "TMath.h"

namespace eff{

  // beam energy in GeV
  const double pbeam = 3500.;
  // masses
  const double Mprot = 0.9382720;
  const double Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot );
  const TLorentzVector beam1_LAB( 0., 0., pbeam, Ebeam );
  const TLorentzVector beam2_LAB( 0., 0., -pbeam, Ebeam );
  const double muMass = 0.105658;
  //rap bins
  Int_t const kNbRapForPTBins = 2;
  Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.6, 1.2};
  Int_t const kNbPTMaxBins = 13;
  Int_t const kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, 12, 12};//all y, y1, y2
  Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
    {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 15., 17., 20., 30.0, 50.},//all rapidities
    {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 13.0, 15., 17., 20., 30.0, 50.},//0-0.6
    {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 13.0, 15., 17., 20., 30.0, 50.}};//0.6-1.2
  //the following values are dummy and need to be filled out at end of analysis
  Double_t pTWCentre_rap[kNbRapForPTBins+1][kNbPTMaxBins+1] = 
    {{5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 14., 18.5, 25., 40., 75.},
     {5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.0, 14., 18.5, 25., 60.},
     {5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.0, 14., 18.5, 25., 60.}};

  //number of reference frames
  enum {CS, HX, PHX, sGJ, GJ1, GJ2};
  Int_t const kNbFrames = 6;
  Char_t *frameLabel[kNbFrames] = {"CS", "HX", "PHX", "sGJ", "GJ1", "GJ2"};
/*   Int_t const kNbFrames = 2; */
/*   Char_t *frameLabel[kNbFrames] = {"CS", "HX"}; */

  //cosTheta 
  Int_t const kNbBinsCosT = 10;
  Double_t cosTMin = -1., cosTMax = 1.;
  //phi for pol. 
   Int_t const kNbBinsPhiPol = 9;
  Double_t phiPolMin = -180., phiPolMax = 180.;

  //pt
  const Int_t nBinsPt1D = 19;
  Double_t pT1D[nBinsPt1D+1] = {0., 0.5, 0.7, 1., 1.25, 1.5, 1.75, 2., 2.5, 3.5, 4., 5., 6., 8., 10., 12., 15., 20., 25., 30.};
  const Int_t nBinsPt2D = 16;
  Double_t pT2D[nBinsPt2D+1] = {0., 0.5, 0.7, 1., 1.25, 1.5, 1.75, 2., 2.5, 3.5, 4., 5., 6., 10., 15., 20., 30.};
  //rap
  const Int_t nBinsRap1D_NP = 16;
  Double_t rap1D_NP[nBinsRap1D_NP+1] = {-2.4, -2.1, -1.5, -1.2, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4};
  const Int_t nBinsRap1D = 8;
  Double_t rap1D[nBinsRap1D+1] = {0.0, 0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4};
  const Int_t nBinsRap2D_NP = 15;
  Double_t rap2D_NP[nBinsRap2D_NP+1] = {-2.4, -2.1, -1.5, -1.2, -0.8, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4};
  const Int_t nBinsRap2D = 8;
  Double_t rap2D[nBinsRap2D+1] = {0.0, 0.2, 0.3, 0.6, 0.8, 1.2, 1.5, 2.1, 2.4};
  //phi-lab
  const Int_t nBinsPhi1D = 18;
  const Int_t nBinsPhi2D = 9;
  Double_t phiMin = -TMath::Pi(), phiMax = TMath::Pi();

  //study the negative and positive rapidity sides separately
  Int_t const kNbRapBins = kNbRapForPTBins;
  /* Double_t rapRange[2*kNbRapBins+1] = {-2.3, -1.9, -1.5, -0.9, 0., 0.9, 1.5, 1.9, 2.3}; */
  /* Double_t rapRange[2*kNbRapBins+1] = {-2.4, -2.1, -1.6, -1.2, -0.9, 0., 0.9, 1.2, 1.6, 2.1, 2.4}; */
  Double_t rapRange[2*kNbRapBins+1] = {-1.2, -0.6, 0., 0.6, 1.2};
  
  //phase space limiting cuts:
  Double_t rapMax = 1.2;
}
