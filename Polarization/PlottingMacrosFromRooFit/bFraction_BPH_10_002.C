// #include "TGraphAsymmErrors.h"
// #include "TGraphErrors.h"
//BPH-10-002
TGraphAsymmErrors *gB_mid, *gB_inter, *gB_FWD;
TGraphErrors *gB_mid_syst, *gB_inter_syst, *gB_FWD_syst;
TGraphErrors *gB_mid_tot, *gB_inter_tot, *gB_FWD_tot;

//==============================
void bFraction_BPH_10_002(){

//   Int_t const kNbRapB = 3;
//   Double_t rapRange[kNbRapB+1] = {0., 1.2, 1.6, 2.4};
  Int_t const kNbpTMid = 3;
  Int_t const kNbpTInter = 4;
  Int_t const kNbpTFWD = 8;
  
  Int_t colourCMS[] = {2, 4, 1};
  Int_t markerCMS[] = {24, 25, 27};
  Double_t markerSizeCMS[] = {0.7, 0.7, 1.2};

  //================================
  //BPH-10-002, 314 nb-1
  //================================
  //mid-rapidity, |y| < 1.2
  Double_t bFraction_mid[kNbpTMid] = {0.179, 0.257, 0.395};
  Double_t bFrac_Stat_mid[kNbpTMid] = {0.024, 0.015, 0.018};
  Double_t bFrac_Syst_mid[kNbpTMid] = {0.015, 0.014, 0.005};
  Double_t bFrac_tot_mid[kNbpTMid];
  Double_t pTMid[kNbpTMid] = {5.46, 8.14, 13.5};
  Double_t pTMid_Err[kNbpTMid] = {0.56, 0.97, 3.53};
  Double_t pTRange_mid[kNbpTMid+1] = {4.5,6.5,10.,30.};
  Double_t pTL_mid[kNbpTMid], pTR_mid[kNbpTMid];
  Double_t errPTSyst_mid[kNbpTMid];
  for(int ipT = 0; ipT < kNbpTMid; ipT++){
    // pTL_mid[ipT] = 0.;
    // pTR_mid[ipT] = 0.;
    pTMid_Err[ipT] = 0.; //don't show the RMS of the data
    pTL_mid[ipT] = pTMid_Err[ipT];
    pTR_mid[ipT] = pTMid_Err[ipT];
    errPTSyst_mid[ipT] = 0.1;
    bFrac_tot_mid[ipT] = sqrt(pow(bFrac_Stat_mid[ipT],2.) + pow(bFrac_Syst_mid[ipT],2.));
  }
  gB_mid = new TGraphAsymmErrors(kNbpTMid, pTMid, bFraction_mid, pTL_mid, pTR_mid, bFrac_Stat_mid, bFrac_Stat_mid);
  gB_mid->SetMarkerStyle(markerCMS[0]);
  gB_mid->SetMarkerColor(colourCMS[0]);
  gB_mid->SetLineColor(colourCMS[0]);
  gB_mid->SetLineWidth(2);
  gB_mid->SetMarkerSize(markerSizeCMS[0]);
  gB_mid_syst = new TGraphErrors(kNbpTMid, pTMid, bFraction_mid, errPTSyst_mid, bFrac_Syst_mid);
  // gB_mid_syst->SetLineColor(colourCMS[0]);
  gB_mid_syst->SetLineColor(colourCMS[0]);
  gB_mid_syst->SetFillStyle(0);

  gB_mid_tot = new TGraphErrors(kNbpTMid, pTMid, bFraction_mid, pTMid_Err, bFrac_tot_mid);
  gB_mid_tot->SetLineColor(colourCMS[0]);
  gB_mid_tot->SetLineWidth(2);
  gB_mid_tot->SetMarkerColor(colourCMS[0]);
  gB_mid_tot->SetMarkerStyle(markerCMS[0]);
  gB_mid_tot->SetMarkerSize(markerSizeCMS[0]);

  //remove the first data point:
  gB_mid->RemovePoint(0);
  gB_mid_syst->RemovePoint(0);
  gB_mid_tot->RemovePoint(0);

  //intermediate-rapidity, 1.2 < |y| < 1.6
  Double_t bFraction_inter[kNbpTInter] = {0.146, 0.180, 0.203, 0.360};
  Double_t bFrac_Stat_inter[kNbpTInter] = {0.021, 0.017, 0.017, 0.031};
  Double_t bFrac_Syst_inter[kNbpTInter] = {0.028, 0.019, 0.014, 0.016};
  Double_t bFrac_tot_inter[kNbpTInter];
  Double_t pTInter[kNbpTInter] = {3.27, 5.48, 7.89, 12.96};
  Double_t pTInter_Err[kNbpTInter] = {0.75, 0.55, 0.93, 3.06};
  Double_t pTRange_inter[kNbpTInter+1] = {2.0,4.5,6.5,10.,30.};
  Double_t pTL_inter[kNbpTInter], pTR_inter[kNbpTInter];
  Double_t errPTSyst_inter[kNbpTInter];
  for(int ipT = 0; ipT < kNbpTInter; ipT++){
    // pTL_inter[ipT] = 0.;
    // pTR_inter[ipT] = 0.;
    pTInter_Err[ipT] = 0.; //don't show the RMS of the data
    pTL_inter[ipT] = pTInter_Err[ipT];
    pTR_inter[ipT] = pTInter_Err[ipT];
    errPTSyst_inter[ipT] = 0.1;
    bFrac_tot_inter[ipT] = sqrt(pow(bFrac_Stat_inter[ipT],2.) + pow(bFrac_Syst_inter[ipT],2.));
  }
  gB_inter = new TGraphAsymmErrors(kNbpTInter, pTInter, bFraction_inter, pTL_inter, pTR_inter, bFrac_Stat_inter, bFrac_Stat_inter);
  gB_inter->SetMarkerStyle(markerCMS[1]);
  gB_inter->SetMarkerColor(colourCMS[1]);
  gB_inter->SetLineColor(colourCMS[1]);
  gB_inter->SetMarkerSize(markerSizeCMS[1]);
  gB_inter->SetLineWidth(2);
  gB_inter_syst = new TGraphErrors(kNbpTInter, pTInter, bFraction_inter, errPTSyst_inter, bFrac_Syst_inter);
  gB_inter_syst->SetLineColor(colourCMS[1]);
  gB_inter_syst->SetFillStyle(0);

  gB_inter_tot = new TGraphErrors(kNbpTInter, pTInter, bFraction_inter, pTInter_Err, bFrac_tot_inter);
  gB_inter_tot->SetLineColor(colourCMS[1]);
  gB_inter_tot->SetLineWidth(2);
  gB_inter_tot->SetMarkerColor(colourCMS[1]);
  gB_inter_tot->SetMarkerStyle(markerCMS[1]);
  gB_inter_tot->SetMarkerSize(markerSizeCMS[1]);

  //forward-rapidity, 1.6 < |y| < 2.4
  Double_t bFraction_FWD[kNbpTFWD] = {0.057, 0.087, 0.113, 0.139, 0.160, 0.177, 0.235, 0.374};
  Double_t bFrac_Stat_FWD[kNbpTFWD] = {0.021, 0.014, 0.013, 0.014, 0.014, 0.012, 0.016, 0.031};
  Double_t bFrac_Syst_FWD[kNbpTFWD] = {0.042, 0.022, 0.020, 0.010, 0.013, 0.012, 0.012, 0.008};
  Double_t bFrac_tot_FWD[kNbpTFWD];
  Double_t pT_FWD[kNbpTFWD] = {0.79, 1.60, 2.35, 3.10, 3.96, 5.35, 7.86, 13.11};
  Double_t pT_FWD_Err[kNbpTFWD] = {0.29, 0.21, 0.22, 0.21, 0.29, 0.57, 0.97, 3.23};
  Double_t pTRange_FWD[kNbpTFWD+1] = {0., 1.25, 2., 2.75, 3.5, 4.5, 6.5, 10., 30.};
  Double_t pTL_FWD[kNbpTFWD], pTR_FWD[kNbpTFWD];
  Double_t errPTSyst_FWD[kNbpTFWD];
  for(int ipT = 0; ipT < kNbpTFWD; ipT++){
    // pTL_FWD[ipT] = 0.;
    // pTR_FWD[ipT] = 0.;
    pT_FWD_Err[ipT] = 0.; //don't show the RMS of the data
    pTL_FWD[ipT] = pT_FWD_Err[ipT];
    pTR_FWD[ipT] = pT_FWD_Err[ipT];
    errPTSyst_FWD[ipT] = 0.1;
    bFrac_tot_FWD[ipT] = sqrt(pow(bFrac_Stat_FWD[ipT],2.) + pow(bFrac_Syst_FWD[ipT],2.));
  }
  gB_FWD = new TGraphAsymmErrors(kNbpTFWD, pT_FWD, bFraction_FWD, pTL_FWD, pTR_FWD, bFrac_Stat_FWD, bFrac_Stat_FWD);
  gB_FWD->SetMarkerStyle(markerCMS[2]);
  gB_FWD->SetMarkerColor(colourCMS[2]);
  gB_FWD->SetLineColor(colourCMS[2]);
  gB_FWD->SetMarkerSize(markerSizeCMS[2]);
  gB_FWD->SetLineWidth(2);
  gB_FWD_syst = new TGraphErrors(kNbpTFWD, pT_FWD, bFraction_FWD, errPTSyst_FWD, bFrac_Syst_FWD);
  gB_FWD_syst->SetLineColor(colourCMS[2]);
  gB_FWD_syst->SetFillStyle(0);
  gB_FWD_tot = new TGraphErrors(kNbpTFWD, pT_FWD, bFraction_FWD, pT_FWD_Err, bFrac_tot_FWD);
  gB_FWD_tot->SetLineColor(colourCMS[2]);
  gB_FWD_tot->SetLineWidth(2);
  gB_FWD_tot->SetMarkerColor(colourCMS[2]);
  gB_FWD_tot->SetMarkerStyle(markerCMS[2]);
  gB_FWD_tot->SetMarkerSize(markerSizeCMS[2]);
}
