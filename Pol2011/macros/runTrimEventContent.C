#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

//====================================
void runTrimEventContent(Double_t fracL = 0.5, 
			 Double_t nSigma = 2.
			 ){

  gROOT->ProcessLine(".L TrimEventContent.C+");

  for(int iState = 0; iState < 3; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 1; iPT < max; iPT++){
	TrimEventContent(iRap, iPT, fracL, nSigma, iState);
      }
    }
  }
}
