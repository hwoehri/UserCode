#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

//====================================
void runCopyTreeEntries(Char_t *fileNameIn = "RootFiles/3Sep2011_noCowboys/selEvents_data_Ups_noCowboys_3Sep2011.root"){

  gROOT->ProcessLine(".L CopyTreeEntries.C+");

  for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
    Int_t max = onia::kNbPTBins[iRap]+1;
    for(int iPT = 1; iPT < max; iPT++){
      CopyTreeEntries(iRap, iPT, fileNameIn);
    }
  }
}
