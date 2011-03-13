Int_t const kNbFrames = 2;
enum {CS, HX};
Char_t *frameLabel[kNbFrames] = {"CS", "HX"};
Int_t const kNbRap = 2;
Char_t *rapLabel[kNbRap] = {"|y| < 0.9", "0.9 < |y| < 1.2"};
//
Int_t const kNbSpecies = 2;
Char_t *speciesLabel[kNbFrames] = {"P", "NP"};
//
Int_t const kNbComp = 4;
Char_t *compName[kNbComp] = {"P", "NP", "BG", "Tot"};
enum {P, NP, BKG, TOT};
enum {FIT,PONE,ZERO,MONE};
//
Int_t const kNbSigBG = 2;
Char_t *sigBGName[kNbSigBG] = {"S", "BG"};
enum {S, BG};
Int_t const kNbBG = 2;
Char_t *BGName[kNbBG] = {"L", "R"};
enum {L, R};
//
Int_t const kNbVars = 4;
enum {LTH, LPHI, LTHPHI, LTILDE};
Char_t *varName[kNbVars] = {"lambdaTheta", "lambdaPhi", "lambdaThetaPhi", "lambdaTilde"};
//
Int_t const kNbVarComb = 3;
enum {TH_PHI, TH_THETAPHI, THETAPHI_PHI};


Int_t const kNbRapForPTBins = 5;
Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.9, 1.2, 1.6, 2.1, 2.4};
Int_t const kNbPTMaxBins = 12;
Int_t const kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, 7,8,9,12,12};//all y, y1, y2, y3, y4, y5
Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
 {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.},//all rapidities
 {0., 6., 7., 8., 10., 15., 20., 30.},//mid-rap
 {0., 4., 6., 7., 8., 10., 15., 20., 30.},
 {0., 4., 5., 6., 7., 8., 10., 15., 20., 30.},
 {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.},
 {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.}};//most forward

