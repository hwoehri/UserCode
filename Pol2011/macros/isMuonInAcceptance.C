//==========================================================
Bool_t isMuonInAcceptance(Int_t iCut, Double_t pT, Double_t eta){

  Double_t etaBorderHLT[2][4] = {{0., 1.1, 1.4, 2.4}, {0., 1.2, 1.3, 2.2}}; //LOOSE, TIGHT cuts
  Double_t pTBorderHLT[2][4] = {{4.6, 4.0, 2.75, 2.0}, {5.2, 4.7, 3.3, 3.0}};
  
  Double_t minPT_HLT;
  Bool_t decision = kFALSE;

  //loop over higher pT muon
  for(int iEta = 0; iEta < 3; iEta++){
    if(fabs(eta) > etaBorderHLT[iCut][iEta] && fabs(eta) < etaBorderHLT[iCut][iEta+1]){
      minPT_HLT = (pTBorderHLT[iCut][iEta+1]-pTBorderHLT[iCut][iEta]) / (etaBorderHLT[iCut][iEta+1]-etaBorderHLT[iCut][iEta]) * (fabs(eta) - etaBorderHLT[iCut][iEta]) + pTBorderHLT[iCut][iEta];
      break;
    }
    else if(fabs(eta) > etaBorderHLT[iCut][3])
      minPT_HLT = 1000.; //reject all events with |eta| > 2.4 (or 2.2, ...)
  }


  if(pT > minPT_HLT)
    decision = kTRUE;
  
  return decision;
}
