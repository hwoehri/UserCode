#include "TLorentzVector.h"
#include "TMath.h"

namespace onia{

  // beam energy in GeV
  const double pbeam = 3500.;
  // masses
  const double Mprot = 0.9382720;
  const double Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot );
  const TLorentzVector beam1_LAB( 0., 0., pbeam, Ebeam );
  const TLorentzVector beam2_LAB( 0., 0., -pbeam, Ebeam );
  const double muMass = 0.105658;
  Int_t const kNbRapForPTBins = 4;
  Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.6, 1.2, 1.8, 2.4}; //relevant for psi' study
  //study the negative and positive rapidity sides separately
  Int_t const kNbRapBins = kNbRapForPTBins;
  Double_t rapRange[2*kNbRapBins+1] = {-2.4, -1.8, -1.2, -0.6, -0., 0.6, 1.2, 1.8, 2.4};

  Int_t const kNbPTMaxBins = 14;
  Int_t const kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, 12, 12, 12, 12};//all y, yi
  Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
    {5., 6., 7., 8., 9., 10., 11., 13., 15., 17., 20., 25., 30., 50., 100.},//all rapidities
    {5., 6., 7., 8., 9., 10., 11., 13., 15., 17., 20., 30., 100.},//0.0-0.6
    {5., 6., 7., 8., 9., 10., 11., 13., 15., 17., 20., 30., 100.},//0.6-1.2
    {5., 6., 7., 8., 9., 10., 11., 13., 15., 17., 20., 30., 100.},//1.2-1.8
    {5., 6., 7., 8., 9., 10., 11., 13., 15., 17., 20., 30., 100.}};//1.8-2.4
  //number of reference frames
  Int_t const kNbFrames = 6;
  Char_t *frameLabel[kNbFrames] = {"CS", "HX", "PHX", "sGJ", "GJ1", "GJ2"};
  enum {CS, HX, PHX, sGJ, GJ1, GJ2};
  //polarization variables
  Int_t const kNbPolVar = 2; //cosTheta, phi
  enum {cosThPol,phiPol};
  //cosTheta
  Int_t const kNbBinsCosT = 10;
  Double_t cosTMin = -1., cosTMax = 1.;
  //phi for pol. 
   Int_t const kNbBinsPhiPol = 10;
  Double_t phiPolMin = -180., phiPolMax = 180.;
  Double_t massMinL = 8.6; //left mass edge of the BG mass window
  Double_t massMaxR = 11.4; //right mass edge of the BG mass window
  Double_t nSigmaL = 4.; //nSigma to define the right edge of the L BG mass window
  Double_t nSigmaR = 3.; //nSigma to define the left edge of the R BG mass window

  //phase space limiting cuts:
  Double_t rapYPS = 2.4;
  Double_t oniaCtauMax = 1000.; //effectively no cut on lifetime
  Double_t nSigMass = 1000.;
  /* Double_t nSigMass = 2.; */
  Double_t polMassOnia[kNbRapForPTBins+1] = {3.6832, 3.6837, 3.6837, 3.6825, 3.6805};//[all rap, rap bin 1-n]
  Double_t sigmaMassOnia[kNbRapForPTBins+1] = {0.036, 0.025, 0.036, 0.047, 0.053};//[all rap, rap bin 1-n]

  //some make up to use the same colour and marker for each pT and rapidity bin
  //in every plotting macro:
  Int_t colour_pT[kNbPTMaxBins+1] = {1, 2, 3, 4, 6, 7, 8, 49, 38, 46, 12, 40};
  Int_t marker_pT[kNbPTMaxBins+1] = {20, 21, 25, 22, 23, 26, 27, 28, 29, 30, 20, 20};
  
  Int_t colour_rapForPTBins[kNbRapForPTBins+1] = {1, 30, 4, 2, 3};
  Int_t marker_rapForPTBins[kNbRapForPTBins+1] = {20, 21, 25, 20, 22};

}
