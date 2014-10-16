#include <iostream>
#include <iomanip>
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLine.h"

using namespace std;


TFile * f[9];
TCanvas * c[3];
TH1F * h[3][11];
THStack * hs[3];
TH1F * h_sum[3];
TH1F * h_sig[3];
TH1F * h_bkg[3];
TH1F * h_ratio[3];



// =============================================================================================================
//  Constants for MC normalization:

const double Lumi2012[] = { 19789.,  // dielectron path    
			    19751.,  // di muon path
			    19780.}; // electron-muon path

const double mc_lumi[] = { 31200./57709905.,    // Wj
			   54.838/10000431.,    // WW
			   33.21/10000283.,     // WZ
			   8.059/9799908.,      // ZZ
			   225.197/6923750.,    // ttbar
			   3503.71/30459503. }; // DY

// =============================================================================================================


void make_plots(const TString h_name="h_M", const bool log_scale=true, const bool bkg_subtraction=false, 
		const int rebin=0, const bool save_plots=false, const double x_low=-999., const double x_up=-999.){


  // ----------------------------------------------------------------------------------
  //  root files

  const int channels_to_plot = 3;

  const TString version = "v12";

  const TString path    = "/gpfs/cms/users/casarsa/analysis/Zc/work/output/" + version;

  const TString dir[] = { "anaEle/",
			  "anaMuo/",
			  "anaEleMuo/" };

  const TString channel[] = { "ee", "mm", "em"};

  vector <TString> file_names;
  // --- data
  file_names.push_back("DoubleElectron_2012_merge.root");
  file_names.push_back("DoubleMu_2012_merge.root");
  file_names.push_back("MuEG_2012_merge.root");
  // --- MC
  file_names.push_back("Wj.root");
  file_names.push_back("WW.root");
  file_names.push_back("WZ.root");
  file_names.push_back("ZZ.root");
  file_names.push_back("TTbar.root");
  file_names.push_back("DYJetsToLL.root");
  //file_names.push_back("DYJetsToLL_PtW_hc.root");



  // ----------------------------------------------------------------------------------
  //  get the data histos from root files

  for (unsigned int ifile=0; ifile<3; ++ifile){

    f[ifile] = new TFile(path + "/" + file_names[ifile]);
    h[ifile][0] =  (TH1F*) f[ifile]->Get(dir[ifile] + h_name)->Clone(); 

    if ( rebin>0 )
      h[ifile][0]->Rebin(rebin);

    // --- histo for all MC channels
    h_sum[ifile] = (TH1F*)  h[ifile][0]->Clone(Form("h_sum_%d",ifile));
    h_sum[ifile]->Reset();

    // --- histo for the MC signal channels
    h_sig[ifile] = (TH1F*)  h[ifile][0]->Clone(Form("h_sig_%d",ifile));
    h_sig[ifile]->Reset();

    // --- histo for the MC bkg channels
    h_bkg[ifile] = (TH1F*)  h[ifile][0]->Clone(Form("h_bkg_%d",ifile));
    h_bkg[ifile]->Reset();
    
  }



  // ----------------------------------------------------------------------------------
  //  get the MC histos from root files

  for (unsigned int ifile=3; ifile<file_names.size(); ++ifile){

    f[ifile] = new TFile(path + "/" + file_names[ifile]);
    
    for (int icha=0;icha<3; ++icha){

      h[icha][ifile-2] = (TH1F*) f[ifile]->Get(dir[icha] + h_name)->Clone(); 

      if (rebin>0)
	h[icha][ifile-2]->Rebin(rebin);

      // --- kludge for the DY file      
      if ( file_names[ifile].Contains("DYJetsToLL") ) {

	h[icha][ifile-1] = (TH1F*) f[ifile]->Get(dir[icha] + h_name + "_tt")->Clone(); 
	h[icha][ifile  ] = (TH1F*) f[ifile]->Get(dir[icha] + h_name + "_b")->Clone(); 
	h[icha][ifile+1] = (TH1F*) f[ifile]->Get(dir[icha] + h_name + "_c")->Clone(); 
	h[icha][ifile+2] = (TH1F*) f[ifile]->Get(dir[icha] + h_name + "_l")->Clone(); 

	if (rebin>0){
	  h[icha][ifile-1]->Rebin(rebin);
	  h[icha][ifile  ]->Rebin(rebin);
	  h[icha][ifile+1]->Rebin(rebin);
	  h[icha][ifile+2]->Rebin(rebin);

	}
      }

      // --- scale MC histograms to data luminosity
      h[icha][ifile-2]->Scale(Lumi2012[icha]*mc_lumi[ifile-3]);

      if ( file_names[ifile].Contains("DYJetsToLL") ) {
	h[icha][ifile-1]->Scale(Lumi2012[icha]*mc_lumi[ifile-3]);
 	h[icha][ifile  ]->Scale(Lumi2012[icha]*mc_lumi[ifile-3]);
 	h[icha][ifile+1]->Scale(Lumi2012[icha]*mc_lumi[ifile-3]);
 	h[icha][ifile+2]->Scale(Lumi2012[icha]*mc_lumi[ifile-3]);
      }
      

      // --- sum histograms
      h_sum[icha]->Add(h[icha][ifile-2]);

      if ( file_names[ifile].Contains("DYJetsToLL") ) {

	h_sum[icha]->Add(h[icha][ifile-1]);
	h_bkg[icha]->Add(h[icha][ifile-1]);

	h_sig[icha]->Add(h[icha][ifile]);
	h_sig[icha]->Add(h[icha][ifile+1]);
	h_sig[icha]->Add(h[icha][ifile+2]);

      }
      else {
	h_bkg[icha]->Add(h[icha][ifile-2]);
      }

    }

  }



  // ----------------------------------------------------------------------------------
  //  printout the yields
  
  cout << "                   ee          mm           em" << endl;
  cout << "---------------------------------------------------" << endl;
  cout <<  " Z+jets    " << setw(12) << h[0][6]->Integral()  << " " 
       << setw(12) << h[1][6]->Integral() << " " 
       << setw(12) << h[2][6]->Integral() << endl;
  cout <<  "      b    " << setw(12) << h[0][8]->Integral()  << " " 
       << setw(12) << h[1][8]->Integral() << " " 
       << setw(12) << h[2][8]->Integral() << endl;
  cout <<  "      c    " << setw(12) << h[0][9]->Integral()  << " " 
       << setw(12) << h[1][9]->Integral() << " " 
       << setw(12) << h[2][9]->Integral() << endl;
  cout <<  "      dusg " << setw(12) << h[0][10]->Integral() << " " 
       << setw(12) << h[1][10]->Integral()<< " " 
       << setw(12) << h[2][10]->Integral()<< endl;
  cout <<  " Z->tautau " << setw(12) << h[0][7]->Integral()  << " " 
       << setw(12) << h[1][7]->Integral() << " " 
       << setw(12) << h[2][7]->Integral() << endl;
  cout <<  " ttbar     " << setw(12) << h[0][5]->Integral()  << " " 
       << setw(12) << h[1][5]->Integral() << " "
       << setw(12) << h[2][5]->Integral() << endl;
  cout <<  " ZZ        " << setw(12) << h[0][4]->Integral()  << " " 
       << setw(12) << h[1][4]->Integral() << " " 
       << setw(12) << h[2][4]->Integral() << endl;
  cout <<  " WZ        " << setw(12) << h[0][3]->Integral()  << " " 
       << setw(12) << h[1][3]->Integral() << " " 
       << setw(12) << h[2][3]->Integral() << endl;
  cout <<  " WW        " << setw(12) << h[0][2]->Integral()  << " " 
       << setw(12) << h[1][2]->Integral() << " " 
       << setw(12) << h[2][2]->Integral() << endl;
  cout <<  " W+jets    " << setw(12) << h[0][1]->Integral()  << " " 
       << setw(12) << h[1][1]->Integral() << " " 
       << setw(12) << h[2][1]->Integral() << endl;
  cout << "---------------------------------------------------" << endl;
  cout <<  " Total bkg " << setw(12) << h_bkg[0]->Integral() << " " 
       << setw(12) << h_bkg[1]->Integral() << " " 
       << setw(12) << h_bkg[2]->Integral() << endl;
  cout <<  " Total MC  " << setw(12) << h_sum[0]->Integral() << " " 
       << setw(12) << h_sum[1]->Integral() << " " << setw(12) 
       << h_sum[2]->Integral() << endl;
  cout << "---------------------------------------------------" << endl;
  cout <<  " Data      " << setw(12) << h[0][0]->Integral()  << " " 
       << setw(12) << h[1][0]->Integral()  << " " 
       << setw(12) << h[2][0]->Integral()  << endl;



  // ----------------------------------------------------------------------------------
  //  make the stacks and plot

  for ( int icha=0; icha<channels_to_plot; ++icha ){

    if ( c[icha] ) 
      delete  c[icha];

    // -- set colors
    h[icha][1]->SetFillColor(kYellow-7);
    h[icha][2]->SetFillColor(kOrange-3);
    h[icha][3]->SetFillColor(kPink+8);
    h[icha][4]->SetFillColor(kAzure+2);
    h[icha][5]->SetFillColor(kRed+1);
    h[icha][6]->SetFillColor(kGreen+2);
    h[icha][7]->SetFillColor(kCyan+1);
    h[icha][8]->SetFillColor(kBlue-8);
    h[icha][9]->SetFillColor(kViolet+2);
    h[icha][10]->SetFillColor(kMagenta+1);


    hs[icha] = new THStack("hs_0",h[icha][0]->GetTitle());
    if ( bkg_subtraction ){
      hs[icha]->Add(h[icha][8]);
      hs[icha]->Add(h[icha][9]);
      hs[icha]->Add(h[icha][10]);
    }
    else {
      hs[icha]->Add(h[icha][2]);
      hs[icha]->Add(h[icha][3]);
      hs[icha]->Add(h[icha][4]);
      hs[icha]->Add(h[icha][5]);
      hs[icha]->Add(h[icha][7]);
      hs[icha]->Add(h[icha][1]);
      hs[icha]->Add(h[icha][6]);
    }


    c[icha] = new TCanvas(Form("c_%d",icha),Form("c_%d",icha), 800, 600);
    c[icha]->cd();

    TPad *pad_1 = new TPad("pad_1","pad_1",0.0,0.3,1.0,1.0);
    pad_1->SetBottomMargin(0.001);
    pad_1->Draw();
    pad_1->cd();
    if (log_scale)
      pad_1->SetLogy(kTRUE);

    hs[icha]->SetMinimum(0.1);
    hs[icha]->Draw("HIST");
    hs[icha]->GetXaxis()->SetTitle(h_name.Data());


    if ( x_low != x_up ){
      int ilow = (int) (x_low-hs[icha]->GetXaxis()->GetBinLowEdge(1))/hs[icha]->GetXaxis()->GetBinWidth(1) +1 ;
      int iup  = (int) (x_up -hs[icha]->GetXaxis()->GetBinLowEdge(1))/hs[icha]->GetXaxis()->GetBinWidth(1);

      hs[icha]->GetXaxis()->SetRange(ilow,iup);
      h[icha][0]->GetXaxis()->SetRange(ilow,iup);
    }

    h[icha][0]->SetLineColor(kBlack);
    h[icha][0]->SetMarkerStyle(20);
    h[icha][0]->SetMarkerColor(kBlack);
    if ( bkg_subtraction ){
      h[icha][0]->Add(h_bkg[icha],-1.);
      h[icha][0]->Draw("EP,SAME");
    }
    else {
      h[icha][0]->Draw("EP,SAME");
    }


    TLegend *leg = new TLegend(0.69, 0.52, 0.96, 0.88);
    leg->SetBorderSize(0);
    leg->SetEntrySeparation(0.01);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);

    if (bkg_subtraction) {
      leg->AddEntry(h[icha][0],"bkg-subtracted data","p");
      leg->AddEntry(h[icha][8],"Z+b","f");
      leg->AddEntry(h[icha][9],"Z+c","f");
      leg->AddEntry(h[icha][10],"Z+dusg","f");
    } else {
      leg->AddEntry(h[icha][0],"Z+jets","p");
      leg->AddEntry(h[icha][5],"t#bar{t}","f");
      leg->AddEntry(h[icha][2],"WW","f");
      leg->AddEntry(h[icha][3],"WZ","f");
      leg->AddEntry(h[icha][4],"ZZ","f");
      leg->AddEntry(h[icha][7],"Z#rightarrow #tau#tau","f");
      leg->AddEntry(h[icha][1],"W+jets","f");
      leg->AddEntry(h[icha][6],"DY","f");
    }
    leg->Draw();


    pad_1->Update();
    c[icha]->Update();
  
    c[icha]->cd();
  
    TPad *pad_2 = new TPad("pad_2","pad_2",0,0,1,0.3);
    pad_2->SetTopMargin(0.);
    pad_2->SetBottomMargin(0.3);
    pad_2->Draw();
    pad_2->cd();


    h_ratio[icha] = (TH1F*) h[icha][0]->Clone("h_ratio");
    if ( bkg_subtraction )
      h_ratio[icha]->Divide(h_sig[icha]);
    else
      h_ratio[icha]->Divide(h_sum[icha]);

    h_ratio[icha]->SetTitle("");
    h_ratio[icha]->SetStats(0);

    h_ratio[icha]->GetXaxis()->SetTitleOffset(0.9);
    h_ratio[icha]->GetXaxis()->SetTitleSize(0.1);
    h_ratio[icha]->GetXaxis()->SetLabelFont(42);
    h_ratio[icha]->GetXaxis()->SetLabelSize(0.08);
    h_ratio[icha]->GetXaxis()->SetTitleFont(42);
    h_ratio[icha]->GetYaxis()->SetTitle("Data/MC ratio");
    h_ratio[icha]->GetYaxis()->SetNdivisions(505);
    h_ratio[icha]->GetYaxis()->SetTitleSize(0.09);
    h_ratio[icha]->GetYaxis()->SetLabelSize(0.08);
    h_ratio[icha]->GetYaxis()->SetRangeUser(0.5, 1.5);
    h_ratio[icha]->GetYaxis()->SetTitleOffset(0.4);

    h_ratio[icha]->SetMarkerStyle(20);
    h_ratio[icha]->SetMarkerColor(kBlack);

    TLine *OLine = new TLine(h_ratio[icha]->GetXaxis()->GetXmin(),1.,h_ratio[icha]->GetXaxis()->GetXmax(),1.);
    OLine->SetLineColor(kRed);
    OLine->SetLineWidth(1);


    if ( x_low != x_up ){
      int ilow = (int) (x_low-hs[icha]->GetXaxis()->GetBinLowEdge(1))/hs[icha]->GetXaxis()->GetBinWidth(1) +1 ;
      int iup  = (int) (x_up -hs[icha]->GetXaxis()->GetBinLowEdge(1))/hs[icha]->GetXaxis()->GetBinWidth(1);
      
      h_ratio[icha]->GetXaxis()->SetRange(ilow,iup);
      OLine->SetX1(x_low);
      OLine->SetX2(x_up);
    }

    h_ratio[icha]->Draw("P");

    OLine->Draw();

    if (save_plots)
      c[icha]->SaveAs( h_name + "_" + channel[icha] + ".pdf");


  }

}

