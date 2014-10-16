#include "TFile.h"
#include "TH2F.h"


TFile * f;

//   0 --> ele + muo
//   1 --> ele
//   2 --> muo

TH2F * h_den_b_0[3];
TH2F * h_den_c_0[3];
TH2F * h_den_l_0[3];
TH2F * h_numL_b_0[3];
TH2F * h_numL_c_0[3];
TH2F * h_numL_l_0[3];
TH2F * h_numT_b_0[3];
TH2F * h_numT_c_0[3];
TH2F * h_numT_l_0[3];

TH2F * h_effL_b_0[3];
TH2F * h_effL_c_0[3];
TH2F * h_effL_l_0[3];
TH2F * h_effT_b_0[3];
TH2F * h_effT_c_0[3];
TH2F * h_effT_l_0[3];

TH2F * h_den_b[3];
TH2F * h_den_c[3];
TH2F * h_den_l[3];
TH2F * h_numL_b[3];
TH2F * h_numL_c[3];
TH2F * h_numL_l[3];
TH2F * h_numT_b[3];
TH2F * h_numT_c[3];
TH2F * h_numT_l[3];

TH2F * h_effL_b[3];
TH2F * h_effL_c[3];
TH2F * h_effL_l[3];
TH2F * h_effT_b[3];
TH2F * h_effT_c[3];
TH2F * h_effT_l[3];


TH2F * histo_rebin(TH2F *h, const int nx, const Double_t *xbins, const int ny, const Double_t *ybins){

  // NB: the errors of the rebinned histogram are wrong!

  TH2F * h_new = new TH2F("h_new",h->GetTitle(),nx,xbins,ny,ybins);
  h_new->Sumw2();

  TAxis *xaxis = h->GetXaxis();
  TAxis *yaxis = h->GetYaxis();
  
  for (int jbin=1; jbin<=yaxis->GetNbins(); ++jbin){
    for (int ibin=1; ibin<=xaxis->GetNbins(); ++ibin){
      h_new->Fill(xaxis->GetBinCenter(ibin),yaxis->GetBinCenter(jbin),h->GetBinContent(ibin,jbin));
    }
  }
  
  return (TH2F*) h_new->Clone();
  
}


void btag_eff_calc(const TString flavour="b", const TString tagger="CSVL",const bool rebin=false){

  const TString path = "/gpfs/cms/users/casarsa/analysis/Zc/work/output/v10/";
  const TString filename = "DYJetsToLL.root";


  f = new TFile(path+filename);


  // --- Read ee and mm histigrams:
  h_numL_b_0[1] = (TH2F*) f->Get("anaEle/hc_CSVL_eff_b");
  h_numL_c_0[1] = (TH2F*) f->Get("anaEle/hc_CSVL_eff_c");
  h_numL_l_0[1] = (TH2F*) f->Get("anaEle/hc_CSVL_eff_l");
  h_numL_b_0[2] = (TH2F*) f->Get("anaMuo/hc_CSVL_eff_b");
  h_numL_c_0[2] = (TH2F*) f->Get("anaMuo/hc_CSVL_eff_c");
  h_numL_l_0[2] = (TH2F*) f->Get("anaMuo/hc_CSVL_eff_l");

  h_numT_b_0[1] = (TH2F*) f->Get("anaEle/hc_CSVT_eff_b");
  h_numT_c_0[1] = (TH2F*) f->Get("anaEle/hc_CSVT_eff_c");
  h_numT_l_0[1] = (TH2F*) f->Get("anaEle/hc_CSVT_eff_l");
  h_numT_b_0[2] = (TH2F*) f->Get("anaMuo/hc_CSVT_eff_b");
  h_numT_c_0[2] = (TH2F*) f->Get("anaMuo/hc_CSVT_eff_c");
  h_numT_l_0[2] = (TH2F*) f->Get("anaMuo/hc_CSVT_eff_l");

  h_den_b_0[1] = (TH2F*) f->Get("anaEle/h_eff_b");
  h_den_c_0[1] = (TH2F*) f->Get("anaEle/h_eff_c");
  h_den_l_0[1] = (TH2F*) f->Get("anaEle/h_eff_l");
  h_den_b_0[2] = (TH2F*) f->Get("anaMuo/h_eff_b");
  h_den_c_0[2] = (TH2F*) f->Get("anaMuo/h_eff_c");
  h_den_l_0[2] = (TH2F*) f->Get("anaMuo/h_eff_l");


  // --- Rebin the histograms
  
  if (rebin) {

    // original binning:
    // 20., 30., 40., 50., 60., 70., 80., 100., 120., 160., 210., 260., 320., 400., 500., 600., 800.

    // for b jets
    const Double_t xbins[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 800};
    const Double_t ybins[] = {-2.4, -1.5, -1.0, -0.5, 0., 0.5, 1.0, 1.5, 2.4};
    const int nx = 9;
    const int ny = 8;

    // for c jets
    //const Double_t xbins[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 800};
    //const Double_t ybins[] = {-2.4, -1., -0., 1., 2.4};
    //const int nx = 9;
    //const int ny = 4;

    
    // for light jets CSVL
    //const Double_t xbins[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 800};
    //const Double_t ybins[] = {-2.4, -1.5, -1.0, -0.5, 0., 0.5, 1.0, 1.5, 2.4};
    //const int nx = 11;
    //const int ny = 8;

    
    // for light jets CSVT
    //const Double_t xbins[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 800};
    //const Double_t ybins[] = {-2.4, 2.4};
    //const int nx = 10;
    //const int ny = 1;


    for (int i=1;i<3; ++i) {
      h_numL_b[i] = histo_rebin( h_numL_b_0[i], nx, xbins, ny, ybins);
      h_numL_c[i] = histo_rebin( h_numL_c_0[i], nx, xbins, ny, ybins);
      h_numL_l[i] = histo_rebin( h_numL_l_0[i], nx, xbins, ny, ybins);
      h_numT_b[i] = histo_rebin( h_numT_b_0[i], nx, xbins, ny, ybins);
      h_numT_c[i] = histo_rebin( h_numT_c_0[i], nx, xbins, ny, ybins);
      h_numT_l[i] = histo_rebin( h_numT_l_0[i], nx, xbins, ny, ybins);
      h_den_b[i]  = histo_rebin( h_den_b_0[i],  nx, xbins, ny, ybins);
      h_den_c[i]  = histo_rebin( h_den_c_0[i],  nx, xbins, ny, ybins); 
      h_den_l[i]  = histo_rebin( h_den_l_0[i],  nx, xbins, ny, ybins); 
    }

  }
  else {

    for (int i=1;i<3; ++i) {
      h_numL_b[i] = h_numL_b_0[i];
      h_numL_c[i] = h_numL_c_0[i];
      h_numL_l[i] = h_numL_l_0[i];
      h_numT_b[i] = h_numT_b_0[i];
      h_numT_c[i] = h_numT_c_0[i];
      h_numT_l[i] = h_numT_l_0[i];
      h_den_b[i]  = h_den_b_0[i];
      h_den_c[i]  = h_den_c_0[i]; 
      h_den_l[i]  = h_den_l_0[i]; 
    }

  }

  
  // --- Sum ee and mm histograms:
  h_numL_b[0] = (TH2F*) h_numL_b[1]->Clone("h_numL_b");
  h_numL_b[0]->Add(h_numL_b[2]);
  h_numL_c[0] = (TH2F*) h_numL_c[1]->Clone("h_numL_c");
  h_numL_c[0]->Add(h_numL_c[2]);
  h_numL_l[0] = (TH2F*) h_numL_l[1]->Clone("h_numL_l");
  h_numL_l[0]->Add(h_numL_l[2]);

  h_numT_b[0] = (TH2F*) h_numT_b[1]->Clone("h_numT_b");
  h_numT_b[0]->Add(h_numT_b[2]);
  h_numT_c[0] = (TH2F*) h_numT_c[1]->Clone("h_numT_c");
  h_numT_c[0]->Add(h_numT_c[2]);
  h_numT_l[0] = (TH2F*) h_numT_l[1]->Clone("h_numT_l");
  h_numT_l[0]->Add(h_numT_l[2]);

  h_den_b[0] = (TH2F*) h_den_b[1]->Clone("h_den_b");
  h_den_b[0]->Add(h_den_b[2]);
  h_den_c[0] = (TH2F*) h_den_c[1]->Clone("h_den_c");
  h_den_c[0]->Add(h_den_c[2]);
  h_den_l[0] = (TH2F*) h_den_l[1]->Clone("h_den_l");
  h_den_l[0]->Add(h_den_l[2]);



  // --- Calculate ratios:
  for (int i=0; i<3; ++i) {

    h_effL_b[i] = (TH2F*) h_numL_b[i]->Clone("h_effL_b");
    h_effL_c[i] = (TH2F*) h_numL_c[i]->Clone("h_effL_c");
    h_effL_l[i] = (TH2F*) h_numL_l[i]->Clone("h_effL_l");
    h_effT_b[i] = (TH2F*) h_numT_b[i]->Clone("h_effT_b");
    h_effT_c[i] = (TH2F*) h_numT_c[i]->Clone("h_effT_c");
    h_effT_l[i] = (TH2F*) h_numT_l[i]->Clone("h_effT_l");

    h_effL_b[i]->Reset();
    h_effL_c[i]->Reset();
    h_effL_l[i]->Reset();
    h_effT_b[i]->Reset();
    h_effT_c[i]->Reset();
    h_effT_l[i]->Reset();

    h_effL_b[i]->Divide(h_numL_b[i],h_den_b[i],1.,1.,"B");
    h_effL_c[i]->Divide(h_numL_c[i],h_den_c[i],1.,1.,"B");
    h_effL_l[i]->Divide(h_numL_l[i],h_den_l[i],1.,1.,"B");
    h_effT_b[i]->Divide(h_numT_b[i],h_den_b[i],1.,1.,"B");
    h_effT_c[i]->Divide(h_numT_c[i],h_den_c[i],1.,1.,"B");
    h_effT_l[i]->Divide(h_numT_l[i],h_den_l[i],1.,1.,"B");

  }


  // --- Printout the CSV efficiencies: 

  TH2F * h;

  if ( flavour=="b" && tagger=="CSVL" )
    h = h_effL_b[0];
  else if ( flavour=="c" && tagger=="CSVL" )
    h = h_effL_c[0];
  else if ( flavour=="l" && tagger=="CSVL" )
    h = h_effL_l[0];
  else if ( flavour=="b" && tagger=="CSVT" )
    h = h_effT_b[0];
  else if ( flavour=="c" && tagger=="CSVT" )
    h = h_effT_c[0];
  else if ( flavour=="l" && tagger=="CSVT" )
    h = h_effT_l[0];
  else {
    cout << "*** Wrong arguments! ***" << endl;
    return;
  }
//  if ( flavour=="b" && tagger=="CSVL" )
//    h = h_numL_b[0];
//  else if ( flavour=="c" && tagger=="CSVL" )
//    h = h_numL_c[0];
//  else if ( flavour=="l" && tagger=="CSVL" )
//    h = h_numL_l[0];
//  else if ( flavour=="b" && tagger=="CSVT" )
//    h = h_numT_b[0];
//  else if ( flavour=="c" && tagger=="CSVT" )
//    h = h_numT_c[0];
//  else if ( flavour=="l" && tagger=="CSVT" )
//    h = h_numT_l[0];
//  else {
//    cout << "*** Wrong arguments! ***" << endl;
//    return;
//  }

  for (unsigned int ipt=1; ipt<=h_effL_b[0]->GetNbinsX(); ++ipt){
    for (unsigned int ieta=1; ieta<=h_effL_b[0]->GetNbinsY(); ++ieta){

      cout << h->GetXaxis()->GetBinLowEdge(ipt) << "\t"
	   << h->GetXaxis()->GetBinLowEdge(ipt) + h->GetXaxis()->GetBinWidth(ipt) << "\t"
	   << h->GetYaxis()->GetBinLowEdge(ieta) << "\t"
	   << h->GetYaxis()->GetBinLowEdge(ieta) + h->GetYaxis()->GetBinWidth(ieta) << "\t"
	   << h->GetBinContent(ipt,ieta) << "\t"
	   << h->GetBinError(ipt,ieta) << "\t"
	   << h->GetBinError(ipt,ieta) 
	   << endl;
    }
  }

  
}
