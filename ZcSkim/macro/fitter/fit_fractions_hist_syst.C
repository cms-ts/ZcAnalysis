///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Notes: 
//       - Prior to compiling this macro, the RooFit include path has to be added:
//         gSystem->AddIncludePath("-I/cvmfs/cms.cern.ch/slc5_amd64_gcc462/lcg/roofit/5.34.02-cms/include")
//         .L fit_fractions.C+
//
//       - In order to use RooMyChebychev, RooMyChebychev.cxx has to be compiled and loaded.
//
// Usage:
//     flag = 0 --> fit the data
//     flag = 1 --> fit the DY MC 
//     flag = 2 --> fit a toy MC
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include <fstream>

#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h" 
#include "RooNLLVar.h"

#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1.h"
#include "TLatex.h"

#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TRandom3.h"


// ========================================================================================

void MultiGaus(const TVectorD& parMeans, const TMatrixDSym& covMatrix, TVectorD& genPars)
{

  TRandom3 rnd(0);

  int nPars = parMeans.GetNrows();
  if(nPars <= 0) {
    Error("MultiGaus", "Must have >0 pars");
    return;
  }
  if(covMatrix.GetNrows() != nPars) {
    Error("MultiGaus", "parMeans.GetNrows() != covMatrix.GetNrows()");
    return;
  }
 
  // Check that covMatrix is symmetric
  for(int iRow = 0; iRow < nPars; iRow++) {
    for(int iCol = iRow; iCol < nPars; iCol++) {
      if(covMatrix(iRow, iCol) != covMatrix(iCol, iRow)) {
        Error("MultiGaus", "malformed cov matrix at row %d, col %d", iRow, iCol);
        return;
      }
    }
  }

  genPars.ResizeTo(nPars);

  TMatrixDSymEigen eigenvariances(covMatrix);
  
  TMatrixD V = eigenvariances.GetEigenVectors();

  TVectorD rotParMeans = V * parMeans;

  for(int iPar = 0; iPar < nPars; iPar++) {
    double variance = eigenvariances.GetEigenValues()[iPar];
    // check for positive-definiteness of covMatrix
    if(variance < 0) {
      Error("MultiGaus", "Got a negative eigenvariance (%f) on iPar = %d", variance, iPar);
    }
    genPars[iPar] = rnd.Gaus(rotParMeans[iPar], sqrt(variance));
  }

  V.Invert();
  
  genPars = V * genPars;

}

// ========================================================================================


//  HISTPDF = 1  --> Use the histograms as pdf's.
#define HISTPDF 0

using namespace RooFit;


TCanvas* canvas   = NULL;
RooPlot* frame    = NULL;
RooPlot* b_frame  = NULL;
RooPlot* c_frame  = NULL;
RooPlot* l_frame  = NULL;

TLatex * text[4];

TFile * f     = NULL;
TFile * f_pdf = NULL;
TH1F  * h     = NULL;

TCanvas * c_bpdf = NULL;
TCanvas * c_cpdf = NULL;
TCanvas * c_lpdf = NULL;

TCanvas* c_fc       = NULL;
TCanvas* c_chi2     = NULL;
TCanvas* c_chi2_pdf = NULL;
TCanvas* c_prob     = NULL;
TCanvas* c_bpars    = NULL;
TCanvas* c_cpars    = NULL;
TCanvas* c_lpars    = NULL;

TH1F h_fc("h_fc","c fraction",100.,0.,0.5);
TH1F h_chi2("h_chi2","chi2",100.,0.,10.);
TH1F h_prob("h_prob","prob",100.,0.,1.);

TH1F h_chi2_b("h_chi2_b","b pdf chi2",100.,0.,100.);
TH1F h_chi2_c("h_chi2_c","c pdf chi2",100.,0.,100.);
TH1F h_chi2_l("h_chi2_l","dusg pdf chi2",100.,0.,100.);

TH1F h_ks_b("h_ks_b","b pdf KS prob",100.,0.,1.);
TH1F h_ks_c("h_ks_c","c pdf KS prob",100.,0.,1.);
TH1F h_ks_l("h_ks_l","dusg pdf KS prob",100.,0.,1.);

TH1F h_b1("h_b1","b1",100.,2.,7.);
TH1F h_b2("h_b2","b2",100.,-7.,-4.);
TH1F h_b3("h_b3","b3",100.,-7.,-4.);
TH1F h_b4("h_b4","b4",100.,1.,2.5);
TH1F h_b5("h_b5","b5",100.,0.3,0.7);

TH1F h_c1("h_c1","c1",100.,1.,4.);
TH1F h_c2("h_c2","c2",100.,-6.,-2.);
TH1F h_c3("h_c3","c3",100.,-4.,-1.);
TH1F h_c4("h_c4","c4",100.,0.5,2.);
TH1F h_c5("h_c5","c5",100.,0.3,0.6);

TH1F h_l1("h_l1","l1",100.,0.45,0.8);
TH1F h_l2("h_l2","l2",100.,0.1,0.4);
TH1F h_l3("h_l3","l3",100.,0.5,2.);
TH1F h_l4("h_l4","l4",100.,0.,0.4);
TH1F h_l5("h_l5","l5",100.,1.5,3.);
TH1F h_l6("h_l6","l6",100.,-6.,-3.);


const int n_toys   = 10000.;
const int n_events = 32000.;


const TString path  = "/gpfs/cms/users/casarsa/analysis/Zc/work/output/v09/";
const TString hname = "hc_svxM_corr";


void fit_fractions_hist_syst(const Int_t n_loop = 1, const Int_t flag=2, const TString channel="ee", 
			     const Bool_t debug=false, const int rebin=4)
{

  int generated_events = 23513 - 1266;
  if ( channel == "mm" ) 
    generated_events = 32719 - 1715;

  // ========================================================================================
  //  Define the fit observable
  // ========================================================================================
  
  RooRealVar x("x","SVX mass",0.,10.);



  // ========================================================================================
  //  Define the fit parameters
  // ========================================================================================

  double f_c_0 = 0.242451;
  double f_l_0 = 0.494342;
  if ( channel=="mm" ){
    f_c_0 = 0.272683;
    f_l_0 = 0.493308;
  }

  RooRealVar f_c("f_c","c-jet fraction",    f_c_0,0.,1.);
  RooRealVar f_l("f_l","dusg-jet fraction", f_l_0,0.,1.);
 


  // ========================================================================================
  //  Build the Likelihood
  // ========================================================================================

#if HISTPDF>0
  
  TString dir = ( channel=="ee" ? "anaEle/" : ("anaMuo/") );


  TFile * f_pdf = new TFile(path + "DYJetsToLL.root");
  TH1F * h_b = (TH1F*) f_pdf->Get(dir + hname + "_b")->Clone("h_b");
  TH1F * h_c = (TH1F*) f_pdf->Get(dir + hname + "_c")->Clone("h_c");
  TH1F * h_l = (TH1F*) f_pdf->Get(dir + hname + "_l")->Clone("h_l");

  if ( rebin>0 ) {
    h_b->Rebin(rebin);
    h_c->Rebin(rebin);
    h_l->Rebin(rebin);
  }

  const int nSmooth = 0;

  // --- b-jet pdf:
  RooHistPdf * Pdf_b = new RooHistPdf("Pdf_b","b pdf", x, 
				      RooDataHist("h_b","",RooArgSet(x),h_b),
				      nSmooth);

  // --- c-jet pdf:
  RooHistPdf * Pdf_c = new RooHistPdf("Pdf_c","c pdf",x,
				      RooDataHist("h_c","",RooArgSet(x),h_c),
				      nSmooth);
  // --- dusg-jet pdf:
  RooHistPdf * Pdf_l = new RooHistPdf("Pdf_l","l pdf",x,
				      RooDataHist("h_l","",RooArgSet(x),h_l),
				      nSmooth);

#else  

  // --- b-jet pdf:

  RooRealVar b0("b0","b0",  3.68511e+02);
  RooRealVar b1("b1","b1",  4.18671e+00);
  RooRealVar b2("b2","b2", -5.43157e+00);
  RooRealVar b3("b3","b3", -5.37153e+00);
  RooRealVar b4("b4","b4",  1.68118e+00);
  RooRealVar b5("b5","b5",  4.74296e-01);
  if (channel=="mm"){
    b0.setVal( 2.15313e+03);
    b1.setVal( 4.95404e+00);
    b2.setVal(-6.02701e+00);
    b3.setVal(-6.03382e+00);
    b4.setVal( 1.85650e+00);
    b5.setVal( 4.98217e-01);
  }

  RooArgList b_pdf_par(b0,b1,b2,b3,b4,b5,x);
  RooGenericPdf * Pdf_b = new RooGenericPdf("Pdf_b","b pdf",
					    "@0*(TMath::Power(@6+@1,@2)+TMath::Exp(@3*@6))*(1.+TMath::Erf((@6-@4)/@5))",
					    b_pdf_par);


  RooRealVar b0_0("b0_0","b0_0",  3.68511e+02);
  RooRealVar b1_0("b1_0","b1_0",  4.18671e+00);
  RooRealVar b2_0("b2_0","b2_0", -5.43157e+00);
  RooRealVar b3_0("b3_0","b3_0", -5.37153e+00);
  RooRealVar b4_0("b4_0","b4_0",  1.68118e+00);
  RooRealVar b5_0("b5_0","b5_0",  4.74296e-01);
  if (channel=="mm"){
    b0_0.setVal( 2.15313e+03);
    b1_0.setVal( 4.95404e+00);
    b2_0.setVal(-6.02701e+00);
    b3_0.setVal(-6.03382e+00);
    b4_0.setVal( 1.85650e+00);
    b5_0.setVal( 4.98217e-01);
  }

  RooArgList b_pdf_par_0(b0_0,b1_0,b2_0,b3_0,b4_0,b5_0,x);
  RooGenericPdf * Pdf_b_0 = new RooGenericPdf("Pdf_b","b pdf",
					      "@0*(TMath::Power(@6+@1,@2)+TMath::Exp(@3*@6))*(1.+TMath::Erf((@6-@4)/@5))",
					      b_pdf_par_0);

  // --- c-jet  pdf:


  RooRealVar c0("c0","c0",  2.82240e+00);
  RooRealVar c1("c1","c1",  2.55373e+00);
  RooRealVar c2("c2","c2", -4.04638e+00);
  RooRealVar c3("c3","c3", -2.60516e+00);
  RooRealVar c4("c4","c4",  1.23763e+00);
  RooRealVar c5("c5","c5",  4.27637e-01);
  if (channel=="mm"){
    c0.setVal( 2.40900e+01);
    c1.setVal( 3.45622e+00);
    c2.setVal(-4.75556e+00);
    c3.setVal(-3.84626e+00);
    c4.setVal( 1.40510e+00);
    c5.setVal( 4.41045e-01);
  }

  RooArgList c_pdf_par(c0,c1,c2,c3,c4,c5,x);
  RooGenericPdf * Pdf_c = new RooGenericPdf("Pdf_c","c pdf",
					    "@0*(TMath::Power(@6+@1,@2)+TMath::Exp(@3*@6))*(1.+TMath::Erf((@6-@4)/@5))",
					    c_pdf_par);

  RooRealVar c0_0("c0_0","c0_0",  2.82240e+00);
  RooRealVar c1_0("c1_0","c1_0",  2.55373e+00);
  RooRealVar c2_0("c2_0","c2_0", -4.04638e+00);
  RooRealVar c3_0("c3_0","c3_0", -2.60516e+00);
  RooRealVar c4_0("c4_0","c4_0",  1.23763e+00);
  RooRealVar c5_0("c5_0","c5_0",  4.27637e-01);
  if (channel=="mm"){
    c0_0.setVal( 2.40900e+01);
    c1_0.setVal( 3.45622e+00);
    c2_0.setVal(-4.75556e+00);
    c3_0.setVal(-3.84626e+00);
    c4_0.setVal( 1.40510e+00);
    c5_0.setVal( 4.41045e-01);
  }

  RooArgList c_pdf_par_0(c0_0,c1_0,c2_0,c3_0,c4_0,c5_0,x);
  RooGenericPdf * Pdf_c_0 = new RooGenericPdf("Pdf_c_0","c pdf",
					      "@0*(TMath::Power(@6+@1,@2)+TMath::Exp(@3*@6))*(1.+TMath::Erf((@6-@4)/@5))",
					      c_pdf_par_0);


  // --- dusg-jet pdf:

  RooRealVar l0("l0","l0",  6.31062e+00);
  RooRealVar l1("l1","l1",  6.12407e-01);
  RooRealVar l2("l2","l2",  2.67040e-01);
  RooRealVar l3("l3","l3",  1.13059e+00);
  RooRealVar l4("l4","l4",  1.80036e-01);
  RooRealVar l5("l5","l5",  2.22999e+00);
  RooRealVar l6("l6","l6", -4.48947e+00);
   if (channel=="mm"){
   l0.setVal( 3.65986e+01);
   l1.setVal( 5.83322e-01);
   l2.setVal( 2.51782e-01);
   l3.setVal( 1.10027e+00);
   l4.setVal( 2.16042e-01);
   l5.setVal( 3.00909e+00);
   l6.setVal(-5.10485e+00);
  }
  
  RooArgList l_pdf_par(l0,l1,l2,l3,l4,l5,l6,x);
  RooGenericPdf * Pdf_l = new RooGenericPdf("Pdf_l","l pdf",
  					    "@0*(2.+TMath::Erf((@7-@1)/@2)+TMath::Erf((@7-@3)/@4))*TMath::Power(@7+@5,@6)",
  					    l_pdf_par);

  RooRealVar l0_0("l0_0","l0_0",  6.31062e+00);
  RooRealVar l1_0("l1_0","l1_0",  6.12407e-01);
  RooRealVar l2_0("l2_0","l2_0",  2.67040e-01);
  RooRealVar l3_0("l3_0","l3_0",  1.13059e+00);
  RooRealVar l4_0("l4_0","l4_0",  1.80036e-01);
  RooRealVar l5_0("l5_0","l5_0",  2.22999e+00);
  RooRealVar l6_0("l6_0","l6_0", -4.48947e+00);
  if (channel=="mm"){
   l0_0.setVal( 3.65986e+01);
   l1_0.setVal( 5.83322e-01);
   l2_0.setVal( 2.51782e-01);
   l3_0.setVal( 1.10027e+00);
   l4_0.setVal( 2.16042e-01);
   l5_0.setVal( 3.00909e+00);
   l6_0.setVal(-5.10485e+00);
  }
   
  RooArgList l_pdf_par_0(l0_0,l1_0,l2_0,l3_0,l4_0,l5_0,l6_0,x);
  RooGenericPdf * Pdf_l_0 = new RooGenericPdf("Pdf_l_0","l pdf",
					      "@0*(2.+TMath::Erf((@7-@1)/@2)+TMath::Erf((@7-@3)/@4))*TMath::Power(@7+@5,@6)",
					      l_pdf_par_0);

  if ( debug ){
    frame = x.frame();
    Pdf_b->plotOn(frame); 
    frame->Draw();
    //return;
  } 


#endif


  // --- build the likelihood:
  
  RooAddPdf model("model","Likelihood",RooArgList(*Pdf_c,*Pdf_l,*Pdf_b),RooArgList(f_c,f_l));
  RooAddPdf model_0("model_0","Likelihood",RooArgList(*Pdf_c_0,*Pdf_l_0,*Pdf_b_0),RooArgList(f_c,f_l));


  // ========================================================================================
  //  Read the data
  // ========================================================================================


  RooDataHist *dataHist = NULL;

  if ( flag==0 || flag==1 || flag==2 ){
    
    // --- Fit the data
    if ( flag==0 ){

      if ( channel == "ee" ){
	//f = new TFile(path + "DoubleElectron_2012_merge.root");
	//h = (TH1F*) f->Get("anaEle/"+hname)->Clone("h");
	f = new TFile("data_noBkg.root");
	h = (TH1F*) f->Get("hc_ee")->Clone("h");
      }
      else if ( channel == "mm" ){
	//f = new TFile(path + "DoubleMu_2012_merge.root");
	//h = (TH1F*) f->Get("anaMuo/"+hname)->Clone("h");
	f = new TFile("data_noBkg.root");
	h = (TH1F*) f->Get("hc_mm")->Clone("h");
      }
      else {
	std::cout << "\n*** ERROR: wrong channel name\n" << std::endl;
	return;
      }

      if ( rebin>0 )
	h->Rebin(rebin);

      dataHist = new RooDataHist("dataHist","x",RooArgSet(x),h);
      
    }

    // --- Fit the DY sample
    else if ( flag==1 ){
      
      f = new TFile(path + "DYJetsToLL.root");
      
      if ( channel == "ee" ){
	h = (TH1F*) f->Get("anaEle/"+hname)->Clone("h");
      }
      else if ( channel == "mm" ){
	h = (TH1F*) f->Get("anaMuo/"+hname)->Clone("h");
      }
      else {
	std::cout << "\n*** ERROR: wrong channel name\n" << std::endl;
	return;
      }
      
      if ( rebin>0 )
	h->Rebin(rebin);

      dataHist = new RooDataHist("dataHist","x",RooArgSet(x),h);
      
    }
    
    // --- Generate a toy MC
    else if ( flag==2 ){
      
      x.setBins(50);
      dataHist = model_0.generateBinned(RooArgSet(x), generated_events);
   
    }

    if ( debug ){
      RooFitResult * fitRes = model_0.fitTo(*dataHist,Save(),SumW2Error(kFALSE),Minimizer("Minuit2")); 
      frame = x.frame();
      dataHist->plotOn(frame);
      model_0.plotOn(frame);
      model_0.plotOn(frame, Components(*Pdf_b_0),LineStyle(kDashed),LineColor(kRed)) ;
      model_0.plotOn(frame, Components(*Pdf_c_0),LineStyle(kDashed),LineColor(kBlue)) ;
      model_0.plotOn(frame, Components(*Pdf_l_0),LineStyle(kDashed),LineColor(kBlack)) ;
      canvas = new TCanvas("canvas","canvas");
      frame->Draw();
      //return;
    }

    // --- Set the parameters of MultiGaus

    // b flavour

    int nDim_h = 5;
    int nDim_l = 6;

    TVectorD parMeans_b(nDim_h);
    parMeans_b(0) = b1.getVal();
    parMeans_b(1) = b2.getVal();
    parMeans_b(2) = b3.getVal();
    parMeans_b(3) = b4.getVal();
    parMeans_b(4) = b5.getVal();
    TMatrixDSym covMatrix_b(nDim_h);
    covMatrix_b(0,0) = 1.671e+00;        covMatrix_b(0,1) = -1.254e+00;       covMatrix_b(0,2) = -1.296e+00;       covMatrix_b(0,3) =  4.291e-01;       covMatrix_b(0,4) =  6.217e-02;
    covMatrix_b(1,0) = covMatrix_b(0,1); covMatrix_b(1,1) =  1.011e+00;       covMatrix_b(1,2) =  1.045e+00;       covMatrix_b(1,3) = -3.317e-01;       covMatrix_b(1,4) = -4.667e-02;
    covMatrix_b(2,0) = covMatrix_b(0,2); covMatrix_b(2,1) = covMatrix_b(1,2); covMatrix_b(2,2) =  1.104e+00;       covMatrix_b(2,3) = -3.351e-01;       covMatrix_b(2,4) = -4.523e-02;
    covMatrix_b(3,0) = covMatrix_b(0,3); covMatrix_b(3,1) = covMatrix_b(1,3); covMatrix_b(3,2) = covMatrix_b(2,3); covMatrix_b(3,3) =  1.143e-01;       covMatrix_b(3,4) =  1.704e-02;
    covMatrix_b(4,0) = covMatrix_b(0,4); covMatrix_b(4,1) = covMatrix_b(1,4); covMatrix_b(4,2) = covMatrix_b(2,4); covMatrix_b(4,3) = covMatrix_b(3,4); covMatrix_b(4,4) =  2.746e-03;
    covMatrix_b.Print();
    TVectorD genPars_b(nDim_l);
    

    // c flavour

    TVectorD parMeans_c(nDim_h);
    parMeans_c(0) = c1.getVal();
    parMeans_c(1) = c2.getVal();
    parMeans_c(2) = c3.getVal();
    parMeans_c(3) = c4.getVal();
    parMeans_c(4) = c5.getVal();
    TMatrixDSym covMatrix_c(nDim_h);
    covMatrix_c(0,0) = 2.401e-01;        covMatrix_c(0,1) = 1.875e-01;        covMatrix_c(0,2) = 1.673e-01;        covMatrix_c(0,3) = -1.856e-02;       covMatrix_c(0,4) = -1.120e-03;
    covMatrix_c(1,0) = covMatrix_c(0,1); covMatrix_c(1,1) = 3.539e-01;        covMatrix_c(1,2) = 3.859e-01;        covMatrix_c(1,3) = -5.883e-02;       covMatrix_c(1,4) = -8.237e-03;
    covMatrix_c(2,0) = covMatrix_c(0,2); covMatrix_c(2,1) = covMatrix_c(1,2); covMatrix_c(2,2) = 4.353e-01;        covMatrix_c(2,3) = -6.739e-02;       covMatrix_c(2,4) = -9.588e-03;
    covMatrix_c(3,0) = covMatrix_c(0,3); covMatrix_c(3,1) = covMatrix_c(1,3); covMatrix_c(3,2) = covMatrix_c(2,3); covMatrix_c(3,3) =  1.112e-02;       covMatrix_c(3,4) =  1.785e-03;
    covMatrix_c(4,0) = covMatrix_c(0,4); covMatrix_c(4,1) = covMatrix_c(1,4); covMatrix_c(4,2) = covMatrix_c(2,4); covMatrix_c(4,3) = covMatrix_c(3,4); covMatrix_c(4,4) =  3.735e-04;
    covMatrix_c.Print();
    TVectorD genPars_c(nDim_l);


    // light flavour

    TVectorD parMeans_l(nDim_l);
    parMeans_l(0) = l1.getVal();
    parMeans_l(1) = l2.getVal();
    parMeans_l(2) = l3.getVal();
    parMeans_l(3) = l4.getVal();
    parMeans_l(4) = l5.getVal();
    parMeans_l(5) = l6.getVal();
    TMatrixDSym covMatrix_l(nDim_l);
    covMatrix_l(0,0) = 7.090e-04;        covMatrix_l(0,1) = 3.629e-04;        covMatrix_l(0,2) = 6.320e-04;        covMatrix_l(0,3) = -7.600e-05;       covMatrix_l(0,4) = -1.336e-02;       covMatrix_l(0,5) =  1.077e-02;
    covMatrix_l(1,0) = covMatrix_l(0,1); covMatrix_l(1,1) = 2.661e-04;        covMatrix_l(1,2) = 2.863e-04;        covMatrix_l(1,3) =  1.026e-05;       covMatrix_l(1,4) = -5.617e-03;       covMatrix_l(1,5) =  4.414e-03;
    covMatrix_l(2,0) = covMatrix_l(0,2); covMatrix_l(2,1) = covMatrix_l(1,2); covMatrix_l(2,2) = 1.198e-03;        covMatrix_l(2,3) = -2.218e-05;       covMatrix_l(2,4) = -1.200e-02;       covMatrix_l(2,5) =  9.095e-03;
    covMatrix_l(3,0) = covMatrix_l(0,3); covMatrix_l(3,1) = covMatrix_l(1,3); covMatrix_l(3,2) = covMatrix_l(2,3); covMatrix_l(3,3) =  1.214e-03;       covMatrix_l(3,4) =  3.792e-03;       covMatrix_l(3,5) = -3.238e-03;
    covMatrix_l(4,0) = covMatrix_l(0,4); covMatrix_l(4,1) = covMatrix_l(1,4); covMatrix_l(4,2) = covMatrix_l(2,4); covMatrix_l(4,3) = covMatrix_l(3,4); covMatrix_l(4,4) =  3.927e-01;       covMatrix_l(4,5) = -3.457e-01;
    covMatrix_l(5,0) = covMatrix_l(0,5); covMatrix_l(5,1) = covMatrix_l(1,5); covMatrix_l(5,2) = covMatrix_l(2,5); covMatrix_l(5,3) = covMatrix_l(3,5); covMatrix_l(5,4) = covMatrix_l(4,5); covMatrix_l(5,5) =  3.127e-01;
    covMatrix_l.Print();
    TVectorD genPars_l(nDim_l);
 
    b_frame = x.frame();
    c_frame = x.frame();
    l_frame = x.frame();

    TH1* h0_b = Pdf_b_0->createHistogram("h0_b",x); 
    h0_b->SetEntries((int) (1.-f_c_0-f_l_0)*generated_events);
    h0_b->Scale((1.-f_c_0-f_l_0)*generated_events/h0_b->Integral());
    TH1* h0_c = Pdf_c_0->createHistogram("h0_c",x); 
    h0_c->SetEntries((int) f_c_0*generated_events);
    h0_c->Scale(f_c_0*generated_events/h0_c->Integral());
    TH1* h0_l = Pdf_l_0->createHistogram("h0_l",x); 
    h0_l->SetEntries((int) f_l_0*generated_events);
    h0_l->Scale(f_l_0*generated_events/h0_l->Integral());

    for (int iloop=0; iloop<n_loop; ++iloop ){

      if ( iloop % 100 == 0 ) 
	cout << " >>>>>>>>>>>> Iteration # " << iloop << endl;


      f_l.setVal(f_l_0);
      f_c.setVal(f_c_0);

      double chi2 = 99999.;
      double prob = -999.;

      double chi2_b = 9999.;
      double chi2_c = 9999.;
      double chi2_l = 9999.;

      double ks_b = -999.;
      double ks_c = -999.;
      double ks_l = -999.;

      //while ( chi2_b>0.5 || std::isnan(chi2_b) ){
      while ( ks_b<0.1 ){

	MultiGaus(parMeans_b, covMatrix_b, genPars_b);
	//genPars_b.Print();

	b1.setVal(genPars_b[0]);
	b2.setVal(genPars_b[1]);
	b3.setVal(genPars_b[2]);
	b4.setVal(genPars_b[3]);
	b5.setVal(genPars_b[4]);
      

	TH1* h1_b = Pdf_b->createHistogram("h1_b",x); 
	h1_b->SetEntries((int) (1.-f_c_0-f_l_0)*generated_events);
	h1_b->Scale((1.-f_c_0-f_l_0)*generated_events/h1_b->Integral());
	//h1_b->Draw("SAME");
	chi2_b = h1_b->Chi2Test(h0_b,"WW CHI2/NDF");
	h_chi2_b.Fill(chi2_b);
	ks_b = h1_b->KolmogorovTest(h0_b);
	h_ks_b.Fill(ks_b);

	delete h1_b;
      }

      //while ( chi2_c>0.5 || std::isnan(chi2_c) ){
      while ( ks_c<0.1 ){

	MultiGaus(parMeans_c, covMatrix_c, genPars_c);
	//genPars_c.Print();

	c1.setVal(genPars_c[0]);
	c2.setVal(genPars_c[1]);
	c3.setVal(genPars_c[2]);
	c4.setVal(genPars_c[3]);
	c5.setVal(genPars_c[4]);

	TH1* h1_c = Pdf_c->createHistogram("h1_c",x); 
	h1_c->SetEntries((int) f_c_0*generated_events);
	h1_c->Scale(f_c_0*generated_events/h1_c->Integral());
	//h1_c->Draw("SAME");
	chi2_c = h1_c->Chi2Test(h0_c,"WW CHI2/NDF");
	h_chi2_c.Fill(chi2_c);
	ks_c = h1_c->KolmogorovTest(h0_c);
	h_ks_c.Fill(ks_c);

	delete h1_c;
      }


      //while ( chi2_l>0.5 || std::isnan(chi2_l) ){
      while ( ks_l<0.1 ){

	MultiGaus(parMeans_l, covMatrix_l, genPars_l);
	//genPars_l.Print();

	l1.setVal(genPars_l[0]);
	l2.setVal(genPars_l[1]);
	l3.setVal(genPars_l[2]);
	l4.setVal(genPars_l[3]);
	l5.setVal(genPars_l[4]);
	l6.setVal(genPars_l[5]);

	TH1* h1_l = Pdf_l->createHistogram("h1_l",x); 
	h1_l->SetEntries((int) f_l_0*generated_events);
	h1_l->Scale(f_l_0*generated_events/h1_l->Integral());
	//h1_l->Draw("SAME");
	chi2_l = h1_l->Chi2Test(h0_l,"WW CHI2/NDF");
	h_chi2_l.Fill(chi2_l);
	ks_l = h1_l->KolmogorovTest(h0_l);
	h_ks_l.Fill(ks_l);

	delete h1_l;
      }

      RooFitResult * fitRes = model.fitTo(*dataHist,Save(),SumW2Error(kFALSE),Minimizer("Minuit2"),PrintLevel(1),Verbose(0));

      if (fitRes->status() != 0 ) continue;
      
      frame = x.frame();
      dataHist->plotOn(frame);
      model.plotOn(frame);
       
      //chi2 = frame->chiSquare("model_Norm[x]","h_dataHist",2);
      chi2 = frame->chiSquare(2);
      prob = TMath::Prob(chi2,2);
      
      h_chi2.Fill(chi2);
      h_prob.Fill(prob);
      

      //if  ( prob < 0.6 ) continue;
      if ( chi2 > 5. ) continue;
      if ( f_c.getVal()<0.01 ||
      	   f_l.getVal()<0.01 ||
	   (1.-f_c.getVal()-f_l.getVal())<0.01 ) continue;
      
      if ( f_c.getVal()>0.99 ||
      	   f_l.getVal()>0.99 ||
      	   (1.-f_c.getVal()-f_l.getVal())>0.99 ) continue;

      
      Pdf_b->plotOn(b_frame);
      Pdf_c->plotOn(c_frame);
      Pdf_l->plotOn(l_frame);
      
      
      h_b1.Fill(genPars_b[0]);
      h_b2.Fill(genPars_b[1]);
      h_b3.Fill(genPars_b[2]);
      h_b4.Fill(genPars_b[3]);
      h_b5.Fill(genPars_b[4]);
      
      h_c1.Fill(genPars_c[0]);
      h_c2.Fill(genPars_c[1]);
      h_c3.Fill(genPars_c[2]);
      h_c4.Fill(genPars_c[3]);
      h_c5.Fill(genPars_c[4]);
      
      h_l1.Fill(genPars_l[0]);
      h_l2.Fill(genPars_l[1]);
      h_l3.Fill(genPars_l[2]);
      h_l4.Fill(genPars_l[3]);
      h_l5.Fill(genPars_l[4]);
      h_l6.Fill(genPars_l[5]);

      if ( debug ) {
	model.plotOn(frame, Components(*Pdf_b),LineStyle(kDashed),LineColor(kRed)) ;
	model.plotOn(frame, Components(*Pdf_c),LineStyle(kDashed),LineColor(kBlue)) ;
	model.plotOn(frame, Components(*Pdf_l),LineStyle(kDashed),LineColor(kBlack)) ;
	canvas = new TCanvas("canvas","canvas");
	frame->Draw();
      }

      h_fc.Fill(f_c.getVal());

    }

    Pdf_b_0->plotOn(b_frame,LineColor(kRed),LineWidth(1));
    Pdf_c_0->plotOn(c_frame,LineColor(kRed),LineWidth(1));
    Pdf_l_0->plotOn(l_frame,LineColor(kRed),LineWidth(1));

    c_chi2 = new TCanvas();
    h_chi2.Draw();
    c_prob = new TCanvas();
    h_prob.Draw();

    c_chi2_pdf = new TCanvas("c_chi2_pdf", "Pdf's chi2", 1024, 768);
    c_chi2_pdf->Divide(3,2);
    c_chi2_pdf->cd(1);
    h_chi2_b.Draw();
    c_chi2_pdf->cd(2);
    h_chi2_c.Draw();
    c_chi2_pdf->cd(3);
    h_chi2_l.Draw();
    c_chi2_pdf->cd(4);
    h_ks_b.Draw();
    c_chi2_pdf->cd(5);
    h_ks_c.Draw();
    c_chi2_pdf->cd(6);
    h_ks_l.Draw();

    c_bpars = new TCanvas("c_bpars", "generated b parameters: distributions", 1024, 768);
    c_bpars->Divide(3,2);
    c_bpars->cd(1);
    h_b1.Draw();
    c_bpars->cd(2);
    h_b2.Draw();
    c_bpars->cd(3);
    h_b3.Draw();
    c_bpars->cd(4);
    h_b4.Draw();
    c_bpars->cd(5);
    h_b5.Draw();

    c_cpars = new TCanvas("c_cpars", "generated c parameters: distributions", 1024, 768);
    c_cpars->Divide(3,2);
    c_cpars->cd(1);
    h_c1.Draw();
    c_cpars->cd(2);
    h_c2.Draw();
    c_cpars->cd(3);
    h_c3.Draw();
    c_cpars->cd(4);
    h_c4.Draw();
    c_cpars->cd(5);
    h_c5.Draw();

    c_lpars = new TCanvas("l_lpars", "generated dusg parameters: distributions", 1024, 768);
    c_lpars->Divide(3,2);
    c_lpars->cd(1);
    h_l1.Draw();
    c_lpars->cd(2);
    h_l2.Draw();
    c_lpars->cd(3);
    h_l3.Draw();
    c_lpars->cd(4);
    h_l4.Draw();
    c_lpars->cd(5);
    h_l5.Draw();
    c_lpars->cd(6);
    h_l6.Draw();

    c_bpdf = new TCanvas();
    b_frame->Draw();
    c_cpdf = new TCanvas();
    c_frame->Draw();
    c_lpdf = new TCanvas();
    l_frame->Draw();

    c_fc = new TCanvas();
    h_fc.Draw();
   
    c_bpdf ->SaveAs("figs/c_bpdf.png");
    c_cpdf ->SaveAs("figs/c_cpdf.png");
    c_lpdf ->SaveAs("figs/c_lpdf.png");
    c_fc   ->SaveAs("figs/c_fc.png");
    c_chi2 ->SaveAs("figs/c_chi2.png");
    c_prob ->SaveAs("figs/c_prob.png");
    c_bpars->SaveAs("figs/c_bpars.png");
    c_cpars->SaveAs("figs/c_cpars.png");
    c_lpars->SaveAs("figs/c_lpars.png");
    c_chi2_pdf->SaveAs("figs/c_chi2_pdf.png");
    
    c_bpdf ->SaveAs("figs/c_bpdf.root");
    c_cpdf ->SaveAs("figs/c_cpdf.root");
    c_lpdf ->SaveAs("figs/c_lpdf.root");
    c_fc   ->SaveAs("figs/c_fc.root");
    c_chi2 ->SaveAs("figs/c_chi2.root");
    c_prob ->SaveAs("figs/c_prob.root");
    c_bpars->SaveAs("figs/c_bpars.root");
    c_cpars->SaveAs("figs/c_cpars.root");
    c_lpars->SaveAs("figs/c_lpars.root");
    c_chi2_pdf->SaveAs("figs/c_chi2_pdf.root");

  }

}
