//==========================================================
Bool_t areMuonsInAcceptance(Int_t iCut, Double_t pTHigh, Double_t etaHigh, Double_t pTLow, Double_t etaLow){

  Double_t etaBorderHLT[2][4] = {{0., 1.1, 1.4, 2.4}, {0., 1.2, 1.3, 2.2}}; //LOOSE, TIGHT cuts
  Double_t pTBorderHLT[2][4] = {{4.6, 4.0, 2.75, 2.0}, {5.2, 4.7, 3.3, 3.0}};
  Double_t etaBorderTM[2][4] = {{0., 1.1, 1.65, 2.4}, {0., 1.2, 1.6, 2.25}};
  Double_t pTBorderTM[2][4] = {{3.5, 3.5, 1.2, 1.2}, {3.8, 3.8, 1.75, 1.75}};
  
  Double_t minPT_HLT, minPT_TM;
  Bool_t decision = kFALSE;

  //loop over higher pT muon
  for(int iEta = 0; iEta < 3; iEta++){
    if(fabs(etaHigh) > etaBorderHLT[iCut][iEta] && fabs(etaHigh) < etaBorderHLT[iCut][iEta+1]){
      minPT_HLT = (pTBorderHLT[iCut][iEta+1]-pTBorderHLT[iCut][iEta]) / (etaBorderHLT[iCut][iEta+1]-etaBorderHLT[iCut][iEta]) * (fabs(etaHigh) - etaBorderHLT[iCut][iEta]) + pTBorderHLT[iCut][iEta];
      break;
    }
    else if(fabs(etaHigh) > etaBorderHLT[iCut][3])
      minPT_HLT = 1000.; //reject all events with |eta| > 2.4 (or 2.2, ...)
  }
  //loop over lower pT muon
  for(int iEta = 0; iEta < 3; iEta++){
    if(fabs(etaLow) > etaBorderTM[iCut][iEta] && fabs(etaLow) < etaBorderTM[iCut][iEta+1]){
      minPT_TM = (pTBorderTM[iCut][iEta+1]-pTBorderTM[iCut][iEta]) / (etaBorderTM[iCut][iEta+1]-etaBorderTM[iCut][iEta]) * (fabs(etaLow) - etaBorderTM[iCut][iEta]) + pTBorderTM[iCut][iEta];
      break;
    }
    else if(fabs(etaLow) > etaBorderTM[iCut][3])
      minPT_TM = 1000.; //reject all events with |eta| > 2.4 (or 2.25, ...)
  }
  if(pTHigh > minPT_HLT && pTLow > minPT_TM)
    decision = kTRUE;

  return decision;
}
