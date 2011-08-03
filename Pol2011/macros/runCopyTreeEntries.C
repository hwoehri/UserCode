#include "../interface/rootIncludes.inc"
#include "../interface/commonVar.h"

//====================================
void runCopyTreeEntries(Char_t *fileNameIn = "RootFiles/selEvents_data_Ups_2Aug2011.root"){

  gROOT->ProcessLine(".L CopyTreeEntries.C+");

  for(int iRap = 1; iRap <= 2; iRap++){
    Int_t max = onia::kNbPTBins[iRap]+1;
    for(int iPT = 0; iPT < max; iPT++){
      CopyTreeEntries(iRap, iPT, fileNameIn);
    }
  }
}
