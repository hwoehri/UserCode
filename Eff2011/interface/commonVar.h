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
  Int_t const kNbPTMaxBins = 10;
  Int_t const kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, 10, 10};//all y, y1, y2
  Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
    {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 16., 20., 30.0, 50.},//all rapidities
    {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 16., 20., 30.0, 50.},//0-0.6
    {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 16., 20., 30.0, 50.}};//0.6-1.2
  //the following values are dummy and need to be filled out at end of analysis
  Double_t pTWCentre_rap[kNbRapForPTBins+1][kNbPTMaxBins+1] = 
    {{5.5, 6.5, 7.5, 8.5, 9.5, 11.5, 14., 18.5, 25., 40.},
     {5.5, 6.5, 7.5, 8.5, 9.5, 12.0, 14., 18.5, 25., 40.},
     {5.5, 6.5, 7.5, 8.5, 9.5, 12.0, 14., 18.5, 25., 40.}};

  //number of reference frames
  enum {CS, HX, PHX, sGJ, GJ1, GJ2};
  Int_t const kNbFrames = 6;
  Char_t *frameLabel[kNbFrames] = {"CS", "HX", "PHX", "sGJ", "GJ1", "GJ2"};
/*   Int_t const kNbFrames = 2; */
/*   Char_t *frameLabel[kNbFrames] = {"CS", "HX"}; */

  //cosTheta 
  Int_t const kNbBinsCosT = 20;
  Int_t const kNbBinsCosT2D = 10;
  Double_t cosTMin = -1., cosTMax = 1.;
  //phi for pol.
  //Int_t const kNbBinsPhiPol2D = 9;
  Int_t const kNbBinsPhiPol2D = 12;
  Int_t const kNbBinsPhiPol = 18;
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
  Double_t rapRange[2*kNbRapBins+1] = {-1.2, -0.6, 0., 0.6, 1.2};
  
  //phase space limiting cuts:
  Double_t rapMax = 1.2;
  //fitted polemass and widths, for the UPS(1S):
  Double_t polMassOnia[kNbRapForPTBins+1] = {9.45, 9.455, 9.454};//[all rap, rap bin 1-n]
  Double_t sigmaMassOnia[kNbRapForPTBins+1] = {0.080, 0.070, 0.095};//[all rap, rap bin 1-n]

  Double_t massMinUps1S_1sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 9.38372, 9.38435, 9.38425, 9.38487, 9.38326},
    {0,0,0,0,0, 9.35644, 9.3574, 9.35555, 9.35522, 9.34896}};
  Double_t massMaxUps1S_1sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 9.518, 9.51711, 9.51826, 9.51701, 9.51941},
    {0,0,0,0,0, 9.53728, 9.53853, 9.5394, 9.54169, 9.54538}};

  Double_t massMinUps2S_1sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 9.94218, 9.94285, 9.94275, 9.9434, 9.94169},
    {0,0,0,0,0, 9.91328, 9.91429, 9.91234, 9.91198, 9.90535}};
  Double_t massMaxUps2S_1sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 10.0844, 10.0835, 10.0847, 10.0834, 10.086},
    {0,0,0,0,0, 10.1049, 10.1062, 10.1071, 10.1096, 10.1135}};

  Double_t massMinUps3S_1sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 10.2715, 10.2722, 10.2721, 10.2728, 10.271},
    {0,0,0,0,0, 10.2416, 10.2427, 10.2407, 10.2403, 10.2335}};
  Double_t massMaxUps3S_1sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 10.4185, 10.4175, 10.4188, 10.4174, 10.42},
    {0,0,0,0,0, 10.4396, 10.441, 10.4419, 10.4444, 10.4485}};

  //=============================

  Double_t massMinUps1S_3sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 9.24945, 9.2516, 9.25025, 9.25273, 9.24711},
    {0,0,0,0,0, 9.1756, 9.17626, 9.17171, 9.16874, 9.15253}};
  Double_t massMaxUps1S_3sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 9.65227, 9.64987, 9.65227, 9.64915, 9.65557},
    {0,0,0,0,0, 9.71811, 9.71966, 9.72048, 9.72149, 9.72017}};

  Double_t massMinUps2S_3sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 9.79992, 9.80219, 9.80076, 9.80339, 9.79744},
    {0,0,0,0,0, 9.72168, 9.72238, 9.72048, 9.72149, 9.72017}};
  Double_t massMaxUps2S_3sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 10.1765, 10.1763, 10.1769, 10.1765, 10.177},
    {0,0,0,0,0, 10.1721, 10.1733, 10.1728, 10.1739, 10.1725}};

  Double_t massMinUps3S_3sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 10.1765, 10.1763, 10.1769, 10.1765, 10.177},
    {0,0,0,0,0, 10.1721, 10.1733, 10.1728, 10.1739, 10.1725}};
  Double_t massMaxUps3S_3sigma[kNbRapForPTBins+1][kNbPTMaxBins] = {
    {0,0,0,0,0, 0,0,0,0,0},
    {0,0,0,0,0, 10.5655, 10.5628, 10.5655, 10.562, 10.5691},
    {0,0,0,0,0, 10.6375, 10.6392, 10.6432, 10.6485, 10.6635}};

}
