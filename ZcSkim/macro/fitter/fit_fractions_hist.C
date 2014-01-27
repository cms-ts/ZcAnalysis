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
//     flag = 3 --> toy MC study
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

//#ifndef __CINT__
//#include "RooGlobalFunc.h"
//#endif
//#include "RooRealVar.h"
//#include "RooArgList.h"
//#include "RooArgSet.h"
//#include "RooChebychev.h"
//#include "RooAddPdf.h"
//#include "RooDataSet.h"
//#include "RooDataHist.h"
//#include "RooMCStudy.h"
//#include "RooPlot.h"
//#include "RooFitResult.h"
//#include "RooGenericPdf.h" 
//#include "RooNLLVar.h"
//
//#include "TCanvas.h"
//#include "TFile.h"
//#include "TH1F.h"
//
//#include <fstream>


#define HISTPDF 0

using namespace RooFit;


TCanvas* canvas = NULL;
RooPlot* frame  = NULL;

TLatex * text[4];

TFile * f     = NULL;
TFile * f_pdf = NULL;
TH1F  * h     = NULL;

const int n_toys   = 10000.;
const int n_events = 185000.;


const TString path  = "/gpfs/cms/users/casarsa/analysis/Zc/work/output/v01/";
const TString hname = "hc_BJP";


void fit_fractions_hist(const Int_t flag=0, const TString channel="ee", 
			const Bool_t verbose=false, const Bool_t rebin=0)
{


  // ========================================================================================
  //  Define the fit observable
  // ========================================================================================
  
  RooRealVar x("x","BJP",0.,10.);



  // ========================================================================================
  //  Define the fit parameters
  // ========================================================================================

  RooRealVar f_b("f_b","b-jet fraction",0.07,0.,1.);
  RooRealVar f_l("f_l","dusg-jet fraction",0.82,0.,1.);
 
  //f_b.setConstant(kTRUE);



  // ========================================================================================
  //  Build the Likelihood
  // ========================================================================================

#if HISTPDF>0
  
  TFile * f_pdf = new TFile(path + "DYJetsToLL.root");
  TH1F * h_b = (TH1F*) f_pdf->Get("anaEle/"+ hname + "_ee_b")->Clone("h_b");
  TH1F * h_c = (TH1F*) f_pdf->Get("anaEle/"+ hname + "_ee_c")->Clone("h_c");
  TH1F * h_l = (TH1F*) f_pdf->Get("anaEle/"+ hname + "_ee_l")->Clone("h_l");

  if ( rebin ) {
    h_b->Rebin();
    h_c->Rebin();
    h_l->Rebin();
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

  //RooRealVar b0("b0","b0", 1.49696e+03);
  //RooRealVar b1("b1","b1", 8.97614e-01);
  //RooRealVar b2("b2","b2", 1.59427e+00);
  //RooRealVar b3("b3","b3", 1.44409e+00);
  //
  //RooArgList b_pdf_par(b0,b1,b2,b3,x);
  //RooGenericPdf * Pdf_b = new RooGenericPdf("Pdf_b","b pdf","@0*TMath::Exp(-@1*@4)*(1.+TMath::Erf((@4-@2)/@3))",b_pdf_par);


  RooRealVar b0 ("b0", "b0",   4.47870e+03);
  RooRealVar b1 ("b1", "b1",   4.47874e+03);
  RooRealVar b2 ("b2", "b2",  -8.78794e+02);
  RooRealVar b3 ("b3", "b3",  -5.54399e+03);
  RooRealVar b4 ("b4", "b4",  -4.37121e+03);
  RooRealVar b5 ("b5", "b5",   1.08433e+03);
  RooRealVar b6 ("b6", "b6",   5.55864e+03);
  RooRealVar b7 ("b7", "b7",   7.13993e+03);
  RooRealVar b8 ("b8", "b8",   5.49716e+03);
  RooRealVar b9 ("b9", "b9",   2.61302e+03);
  RooRealVar b10("b10","b10",  9.24420e+02);

  RooArgList b_pdf_par(b0,b1,b2,b3,b4,b5,b6,b7,b8);
  b_pdf_par.add(b9);  
  b_pdf_par.add(b10);  

  RooMyChebychev *  Pdf_b = new RooMyChebychev("Pdf_b", "b pdf", x, b_pdf_par);
  


  // --- c-jet  pdf:

  //RooRealVar c0("c0","c0", 1.23326e+04);
  //RooRealVar c1("c1","c1", 1.55897e+00);
  //RooRealVar c2("c2","c2", 1.55481e+00);
  //RooRealVar c3("c3","c3", 1.21464e+00);
  //
  //RooArgList c_pdf_par(c0,c1,c2,c3,x);
  //RooGenericPdf * Pdf_c = new RooGenericPdf("Pdf_c","c pdf","@0*TMath::Exp(-@1*@4)*(1.+TMath::Erf((@4-@2)/@3))",c_pdf_par);


  RooRealVar c0 ("c0", "c0", 8.47143e-01);
  RooRealVar c1 ("c1", "c1", 1.49919e+00);
  RooRealVar c2 ("c2", "c2", 1.15742e+00);
  RooRealVar c3 ("c3", "c3", 8.51517e-01);
  RooRealVar c4 ("c4", "c4", 8.52124e-01);
  RooRealVar c5 ("c5", "c5", 1.06918e+00);
  RooRealVar c6 ("c6", "c6", 1.19842e+00);
  RooRealVar c7 ("c7", "c7", 1.04439e+00);
  RooRealVar c8 ("c8", "c8", 6.48985e-01);
  RooRealVar c9 ("c9", "c9", 2.70436e-01);
  RooRealVar c10("c10","c10",6.29495e-02);
  //RooRealVar c11("c11","c11",);
  //RooRealVar c12("c12","c12",);


  RooArgList c_pdf_par(c0,c1,c2,c3,c4,c5,c6,c7,c8);
  c_pdf_par.add(c9);  
  c_pdf_par.add(c10);  
  //c_pdf_par.add(c11);  
  //c_pdf_par.add(c12);  

  RooMyChebychev *  Pdf_c = new RooMyChebychev("Pdf_c", "c pdf", x, c_pdf_par);


  // --- dusg-jet pdf:

  RooRealVar l0("l0","l0", 1.21171e+05);
  RooRealVar l1("l1","l1", 1.53030e+00);
  RooRealVar l2("l2","l2", 7.25955e-01);
  RooRealVar l3("l3","l3", 8.11109e-01);
  
  RooArgList l_pdf_par(l0,l1,l2,l3,x);
  RooGenericPdf * Pdf_l = new RooGenericPdf("Pdf_l","dusg pdf","@0*TMath::Exp(-@1*@4)*(1.+TMath::Erf((@4-@2)/@3))",l_pdf_par);

  //RooRealVar l0 ("l0", "l0",  -7.37000e+03);
  //RooRealVar l1 ("l1", "l1",  -2.64181e+04);
  //RooRealVar l2 ("l2", "l2",   9.76844e+03);
  //RooRealVar l3 ("l3", "l3",   1.01456e+03);
  //RooRealVar l4 ("l4", "l4",  -1.11664e+03);
  //RooRealVar l5 ("l5", "l5",  -9.17895e+03);
  //RooRealVar l6 ("l6", "l6",  -4.81038e+03);
  //RooRealVar l7 ("l7", "l7",   9.99850e+03);
  //RooRealVar l8 ("l8", "l8",   2.99891e+03);
  //RooRealVar l9 ("l9", "l9",  -4.28644e+03);
  //RooRealVar l10("l10","l10", -1.38560e+04);
  //RooRealVar l11("l11","l11", -4.44506e+03);
  //RooRealVar l12("l12","l12",  1.03354e+04);
  //RooRealVar l13("l13","l13",  1.85403e+04);
  //RooRealVar l14("l14","l14",  1.72235e+04);
  //RooRealVar l15("l15","l15",  8.43691e+03);
  //RooRealVar l16("l16","l16",  2.96895e+03);
  ////RooRealVar l17("l17","l17",  5.54591e+03);
  ////RooRealVar l18("l18","l18",  6.54835e+03);
  ////RooRealVar l19("l19","l19",  3.21628e+03);
  ////RooRealVar l20("l20","l20",  1.28814e+03);
  //
  //RooArgList l_pdf_par(l0,l1,l2,l3,l4,l5,l6,l7,l8);
  //l_pdf_par.add(l9);  
  //l_pdf_par.add(l10);  
  //l_pdf_par.add(l11);  
  //l_pdf_par.add(l12);  
  //l_pdf_par.add(l13);  
  //l_pdf_par.add(l14);  
  //l_pdf_par.add(l15);  
  //l_pdf_par.add(l16);  
  ////l_pdf_par.add(l17);  
  ////l_pdf_par.add(l18);  
  ////l_pdf_par.add(l19);  
  ////l_pdf_par.add(l20);  
  //
  //RooMyChebychev *  Pdf_l = new RooMyChebychev("Pdf_l", "l pdf", x, l_pdf_par);


  //fr = x.frame();
  //Pdf_l->plotOn(fr);
  //fr->Draw();
  //return;



#endif


  // --- build the likelihood:
  
  RooAddPdf model("model","Likelihood",RooArgList(*Pdf_b,*Pdf_l,*Pdf_c),RooArgList(f_b,f_l));



  // ========================================================================================
  //  Read the data
  // ========================================================================================


  RooDataHist *dataHist = NULL;

  if ( flag==0 || flag==1 || flag==2 ){
    
    // --- Fit the data
    if ( flag==0 ){

      if ( channel == "ee" ){
	f = new TFile(path + "DoubleElectron_2012_merge.root");
	h = (TH1F*) f->Get("anaEle/"+hname+"_ee")->Clone("h");
      }
      else if ( channel == "mm" ){
	f = new TFile(path + "DoubleMu_2012_merge.root");
	h = (TH1F*) f->Get("anaMuo/"+hname+"_mm")->Clone("h");
      }
      else {
	std::cout << "\n*** ERROR: wrong channel name\n" << std::endl;
	return;
      }

      if ( rebin )
	h->Rebin();

      dataHist = new RooDataHist("dataHist","x",RooArgSet(x),h);
      
    }

    // --- Fit the DY sample
    else if ( flag==1 ){
      
      f = new TFile(path + "DYJetsToLL.root");
      
      if ( channel == "ee" ){
	h = (TH1F*) f->Get("anaEle/"+hname+"_ee")->Clone("h");
      }
      else if ( channel == "mm" ){
	h = (TH1F*) f->Get("anaMuo/"+hname+"_mm")->Clone("h");
      }
      else {
	std::cout << "\n*** ERROR: wrong channel name\n" << std::endl;
	return;
      }
      
      if ( rebin )
	h->Rebin();

      dataHist = new RooDataHist("dataHist","x",RooArgSet(x),h);
      
    }
    
    // --- Generate a toy MC
    else if ( flag==2 ){
      
      RooDataSet * data = model.generate(x,10000);
      
      x.setBins(100) ;
      dataHist = new RooDataHist("dataHist","x",RooArgSet(x),*data) ;
      
    }


    // ========================================================================================
    //  Perform the fit
    // ========================================================================================

    RooFitResult * fitRes = model.fitTo(*dataHist,Save()); 


    // ========================================================================================
    //  Plot the data and the composite PDF overlaid
    // ========================================================================================

    frame = x.frame();
    dataHist->plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(*Pdf_b),LineStyle(kDashed),LineColor(kRed)) ;
    model.plotOn(frame, Components(*Pdf_c),LineStyle(kDashed),LineColor(kBlue)) ;
    model.plotOn(frame, Components(*Pdf_l),LineStyle(kDashed),LineColor(kBlack)) ;

    canvas = new TCanvas("canvas","canvas");
    frame->Draw();

    TLegend *leg = new TLegend(0.55, 0.55, 0.875, 0.875);
    //leg->AddEntry(frame->findObject("h_dataHist"),"Data","LP");
    leg->AddEntry(frame->findObject("h_dataHist"),"DY MC","LP");
    leg->AddEntry(frame->findObject("model_Norm[x]"),"Fit projection","L");
    leg->AddEntry(frame->findObject("model_Norm[x]_Comp[Pdf_b]"),"b component","L");
    leg->AddEntry(frame->findObject("model_Norm[x]_Comp[Pdf_c]"),"c component","L");
    leg->AddEntry(frame->findObject("model_Norm[x]_Comp[Pdf_l]"),"dusg component","L");
    leg->Draw();

    double f_c =  1.-f_b.getVal()-f_l.getVal();
    double f_c_err =  sqrt(f_b.getError()*f_b.getError()+f_l.getError()*f_l.getError()+
			   2.*fitRes->correlation(f_b,f_l)*f_b.getError()*f_l.getError());

    Double_t chi2 = frame.chiSquare("model_Norm[x]","h_dataHist",2);


    text[0] = new TLatex(4., 1500., Form("f_{b}    = %4.3f #pm %4.3f", f_b.getVal(),  f_b.getError()));
    text[0]->Draw();
    text[1] = new TLatex(4., 1000., Form("f_{dusg} = %4.3f #pm %4.3f", f_l.getVal(),  f_l.getError()));
    text[1]->Draw();
    text[2] = new TLatex(4., 500., Form("f_{c}    = %4.3f #pm %4.3f", f_c,  f_c_err));
    text[2]->SetTextColor(4);
    text[2]->Draw();
    text[3] = new TLatex(2.5, 3500., Form("#chi^{2}/d.o.f.  = %4.2f", chi2));
    text[3]->Draw();

    if ( verbose ) {

      // Plot the likelihood scans
      RooNLLVar nll("nll","nll",model,*dataHist) ;

      RooPlot* fbscan_frame = f_b.frame(Range(0.,0.3),Title("-log(L) scan vs f_{b}"));
      nll.plotOn(fbscan_frame,PrintEvalErrors(0),ShiftToZero(),LineColor(kRed),Precision(1e-4));
      //fbscan_frame->SetMaximum(15);
      //fbscan_frame->SetMinimum(0);

      TCanvas * fbscan_canvas = new TCanvas("fbscan_canvas","-log(L) scan vs fb");
      fbscan_frame->Draw();

      RooPlot* flscan_frame = f_l.frame(Range(0.4,1.),Title("-log(L) scan vs f_{l}"));
      nll.plotOn(flscan_frame,PrintEvalErrors(0),ShiftToZero(),LineColor(kRed),Precision(1e-4));
      //fbscan_frame->SetMaximum(15);
      //fbscan_frame->SetMinimum(0);

      TCanvas * flscan_canvas = new TCanvas("flscan_canvas","-log(L) scan vs fl");
      flscan_frame->Draw();

    }


    // ========================================================================================
    //  Print the fit result
    // ========================================================================================

    cout << endl;
    cout << "Channel = " << channel << endl;
    cout << endl;
    dataHist->Print();
    cout << endl;

    //fitRes->Print();

    cout << "f_b = " << f_b.getVal() << " +/- " << f_b.getError() << endl;
    cout << "f_l = " << f_l.getVal() << " +/- " << f_l.getError() << endl;
    cout << "f_c = " 
	 << 1.-f_b.getVal()-f_l.getVal() << " +/- " 
	 << sqrt(f_b.getError()*f_b.getError()+f_l.getError()*f_l.getError()+
		 2.*fitRes->correlation(f_b,f_l)*f_b.getError()*f_l.getError()) 
	 << endl;
    cout << endl;

  }

  // ========================================================================================
  //  Toy MC study
  // ========================================================================================
  else if ( flag==3 ) {

    // Create manager
    RooMCStudy* mcstudy = new RooMCStudy(model,x,Binned(kTRUE),Silence(),Extended(),
					 FitOptions(Save(kTRUE),PrintEvalErrors(0)));

    // Generate and fit  events
    mcstudy->generateAndFit(n_toys,n_events);

    // Plot the distrution of the pull
    RooPlot * fbpull_frame = mcstudy->plotPull(f_b,-3.,3.,40,kTRUE);
    TCanvas * fbpull_canvas = new TCanvas("fbpull_canvas","fb pull");
    fbpull_frame->Draw();
    

    // Plot the distrution of the pull
    RooPlot * flpull_frame = mcstudy->plotPull(f_l,-3.,3.,40,kTRUE);
    TCanvas * flpull_canvas = new TCanvas("flpull_canvas","fl pull");
    flpull_frame->Draw();

    if ( verbose ) {

      // Plot the distrution of the value 
      RooPlot * fb_frame = f_b.frame(0.,0.5); 
      mcstudy->plotParamOn(fb_frame); 
      TCanvas * fb_canvas = new TCanvas("fb_canvas","fb");
      fb_frame->Draw();

      // Plot the distrution of the error 
      RooPlot * fberr_frame = mcstudy->plotError(f_b,0.,0.1);
      TCanvas * fberr_canvas = new TCanvas("fberr_canvas","fb error");
      fberr_frame->Draw();

      // Plot the distribution of the NLL
      RooPlot * nll_frame = mcstudy->plotNLL();
      TCanvas * nll_canvas = new TCanvas("nll_canvas","NLL");
      nll_frame->Draw();

    }

  }
  else {
    cout << "\n*** ERROR: wrong flag value\n" << endl;
    return;
  }

}
