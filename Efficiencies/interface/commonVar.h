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
  Int_t const kNbRapForPTBins = 5;
  /* Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.9, 1.5, 1.9, 2.3}; */
  Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.9, 1.2, 1.6, 2.1, 2.4};
  //pT bins (optimized to have at least 10.000 entries per bin)
  Int_t const kNbPTMaxBins = 12;
  Int_t const kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, 7,8,9,12,12};//all y, y1, y2, y3, y4, y5
  Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
    {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.},//all rapidities
    {0., 6., 7., 8., 10., 15., 20., 30.},//mid-rap
    {0., 4., 6., 7., 8., 10., 15., 20., 30.},
    {0., 4., 5., 6., 7., 8., 10., 15., 20., 30.},
    {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.},
    {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.}};//most forward
  //the following values are dummy and need to be filled out at end of analysis
  Double_t pTWCentre_rap[kNbRapForPTBins+1][23] = 
    {{0.5, 1.25, 1.75, 1.95, 2.5, 2.8, 3.1, 3.4, 3.9, 4.3, 4.6, 5.0, 5.5, 6.0, 6.5, 7.2, 7.8, 8.5, 9.6, 11.0, 14., 27.},
     {7., 12.},
     {4.5, 6.0, 7.3, 15.},
     {1.0, 2.5, 3.5, 4.5, 5.5, 6.7, 10.},
     {0.5, 1.5, 2.5, 3.5, 5.0, 13.},
     {0.5, 1.5, 2.5, 3.5, 5.0, 13.}};
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
  Double_t rapRange[2*kNbRapBins+1] = {-2.4, -2.1, -1.6, -1.2, -0.9, 0., 0.9, 1.2, 1.6, 2.1, 2.4};
  
  //phase space limiting cuts:
  Int_t const kNbEtaRegions = 3;
  Double_t etaPS[kNbEtaRegions] = {1.3, 2.2, 2.4}; //pseudo-rap cuts for muons
  Double_t pTMuMin[kNbEtaRegions] = {3.3, 0., 0.8};
  Double_t pMuMin[kNbEtaRegions] = {0., 2.9, 0.};

}
