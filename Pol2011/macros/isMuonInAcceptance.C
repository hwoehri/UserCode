//==========================================================
Bool_t isMuonInAcceptance(Int_t iCut, Double_t pT, Double_t eta){

  //inclusive tracker muons + efficiency triggers only go down to 2 GeV/c:
  Double_t etaBorder[2][4] = {{0., 1.2, 1.6, 2.1}, {0., 1.2, 1.6, 2.1}}; //tracker50, tracker80
  Double_t pTBorder[2][4] = {{3.5, 3.5, 2.5, 2.5}, {3.8, 3.8, 2.5, 2.5}};
  // Double_t etaBorder[2][4] = {{0., 1.2, 1.6, 2.4}, {0., 1.2, 1.6, 2.1}}; //tracker50, tracker80
  // Double_t pTBorder[2][4] = {{3.5, 3.5, 2.0, 2.0}, {3.8, 3.8, 2.0, 2.0}};

  // //inclusive tracker muons:
  // Double_t etaBorder[2][4] = {{0., 1.1, 1.65, 2.4}, {0., 1.2, 1.6, 2.25}}; //LOOSE, TIGHT cuts
  // Double_t pTBorder[2][4] = {{3.5, 3.5, 1.2, 1.2}, {3.8, 3.8, 1.75, 1.75}};

  //exclusive global muons
  // Double_t etaBorder[2][4] = {{0., 1.1, 1.4, 2.4}, {0., 1.2, 1.3, 2.2}}; //LOOSE, TIGHT cuts
  // Double_t pTBorder[2][4] = {{4.6, 4.0, 2.75, 2.0}, {5.2, 4.7, 3.3, 3.0}};
  
  Double_t minPT;
  Bool_t decision = kFALSE;

  //loop over higher pT muon
  for(int iEta = 0; iEta < 3; iEta++){
    if(fabs(eta) > etaBorder[iCut][iEta] && fabs(eta) < etaBorder[iCut][iEta+1]){
      minPT = (pTBorder[iCut][iEta+1]-pTBorder[iCut][iEta]) / (etaBorder[iCut][iEta+1]-etaBorder[iCut][iEta]) * (fabs(eta) - etaBorder[iCut][iEta]) + pTBorder[iCut][iEta];
      break;
    }
    else if(fabs(eta) > etaBorder[iCut][3])
      minPT = 1000.; //reject all events with |eta| > 2.4 (or 2.2, ...)
  }


  if(pT > minPT)
    decision = kTRUE;
  
  return decision;
}
