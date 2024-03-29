#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

//====================================
void runMassFit(Double_t nSigma = 2., //sigma_mass cut for preparation of figures
		Char_t *fileNameIn = "RootFiles/3Sep2011_noCowboys/selEvents_data_Ups_noCowboys_3Sep2011.root"){

  gROOT->ProcessLine(".L upsilon_2StepFit.C+");

  for(int iRap = 0; iRap <= onia::kNbRapForPTBins; iRap++){
    Int_t max = onia::kNbPTBins[iRap]+1;
    for(int iPT = 0; iPT < max; iPT++){
      upsilon_2StepFit(iRap, iPT, nSigma, fileNameIn);
    }
  }
}
