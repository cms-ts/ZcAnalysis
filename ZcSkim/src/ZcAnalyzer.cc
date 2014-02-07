// -*- C++ -*-
//
// Package: ZcAnalyzer
// Class: ZcAnalyzer
//
/**\class ZcAnalyzer ZcAnalyzer.cc ZcAnalysis/ZcAnalyzer/src/ZcAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]

*/
//
// Original Authors: Alessandra Trevisan and Massimo Casarsa,
//                   inspired by ZbAnalysis/ZbAnalyzer/src/ZbAnalyzer.cc
// Created: Thu Sep 12 11:03:00 CET 2013
// $Id: ZcAnalyzer.cc,v 1.111 2013/07/20 07:23:57 dellaric Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <string>


// CMSSW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"


// root include files
#include <Math/VectorUtil.h>
#include "TH1F.h"
#include "TH2F.h"


// user include files
#include "table.h"


//
// class declaration
//

class ZcAnalyzer : public edm::EDProducer {

public:

  explicit ZcAnalyzer (const edm::ParameterSet &);
  ~ZcAnalyzer ();

private:

  virtual void beginJob ();
  virtual void produce (edm::Event &, const edm::EventSetup &);
  virtual void endJob ();

  virtual void beginRun (edm::Run const &, edm::EventSetup const &);
  virtual void endRun (edm::Run const &, edm::EventSetup const &);
  virtual void beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);
  virtual void endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);


  double ctagWeight(bool isMC, std::vector< pat::Jet >& jets) {
    
    if (isMC == false) return 1.;

    double weight = 1.;

    for (unsigned int i=0; i<jets.size(); i++) {
      
      double discrCSV = jets[i].bDiscriminator("combinedSecondaryVertexBJetTags");

      if ( fabs(jets[i].partonFlavour()) == 5 ) {
	
	double SFL  = BtSFL_->Val(jets[i].pt(), jets[i].eta());
	double SFT  = BtSFT_->Val(jets[i].pt(), jets[i].eta());
	double effL = BtEffL_->Val(jets[i].pt(), jets[i].eta());
	double effT = BtEffT_->Val(jets[i].pt(), jets[i].eta());
      
	if ( discrCSV>0.244 && discrCSV<0.898 )
	  weight *= (SFL*effL-SFT*effT)/(effL-effT);
	else if ( discrCSV<=0.244 )
	  weight *= (1.-SFL*effL)/(1.-effL);
	else if ( discrCSV>=0.898 )
	  weight *= SFT;
	
      }
      else if ( fabs(jets[i].partonFlavour()) == 4 ) {
	
	double SFL  = BtSFL_->Val(jets[i].pt(), jets[i].eta());
	double SFT  = BtSFT_->Val(jets[i].pt(), jets[i].eta());
	double effL = CtEffL_->Val(jets[i].pt(), jets[i].eta());
	double effT = CtEffT_->Val(jets[i].pt(), jets[i].eta());
      
	if ( discrCSV>0.244 && discrCSV<0.898 )
	  weight *= (SFL*effL-SFT*effT)/(effL-effT);
	else if ( discrCSV<=0.244 )
	  weight *= (1.-SFL*effL)/(1.-effL);
	else if ( discrCSV>=0.898 )
	  weight *= SFT;
	  
      }
      else {
	
	double SFL  = LtSFL_->Val(jets[i].pt(), fabs(jets[i].eta()));
	double SFT  = LtSFT_->Val(jets[i].pt(), jets[i].eta());
	double effL = LtEffL_->Val(jets[i].pt(), jets[i].eta());
	double effT = LtEffT_->Val(jets[i].pt(), jets[i].eta());
      
	if ( discrCSV>0.244 && discrCSV<0.898 )
	  weight *= (SFL*effL-SFT*effT)/(effL-effT);
	else if ( discrCSV<=0.244 )
	  weight *= (1.-SFL*effL)/(1.-effL);
	else if ( discrCSV>=0.898 )
	  weight *= SFT;
	  
      }


    }

    return weight;

  };


  double jetResolutionCorrection(double jetEta, double jetPt, double jetPtGen, int syst) {

    if (jetPt <= 0 || jetPtGen <= 0) return jetPt;

    double correctionFactor[5]     = {1.052, 1.057, 1.096, 1.134, 1.288};
    double correctionFactorUp[5]   = {1.115, 1.114, 1.161, 1.228, 1.488};
    double correctionFactorDown[5] = {0.990, 1.001, 1.032, 1.042, 1.089};

    int index = 0;

    if (fabs(jetEta) <= 0.5) index = 0;
    if (fabs(jetEta) > 0.5 && fabs(jetEta) <= 1.1) index = 1;
    if (fabs(jetEta) > 1.1 && fabs(jetEta) <= 1.7) index = 2;
    if (fabs(jetEta) > 1.7 && fabs(jetEta) <= 2.3) index = 3;
    if (fabs(jetEta) > 2.3 && fabs(jetEta) <= 3.0) index = 4;

    double jetPtNew = jetPt;

    if (syst ==  0) jetPtNew = jetPtGen + correctionFactor[index]     * (jetPt-jetPtGen);
    if (syst == +1) jetPtNew = jetPtGen + correctionFactorUp[index]   * (jetPt-jetPtGen);
    if (syst == -1) jetPtNew = jetPtGen + correctionFactorDown[index] * (jetPt-jetPtGen);

    return fmax(0.0, jetPtNew/jetPt);

  };



  // ----------member data ---------------------------

  std::string pileup_;
  std::string pileupMC_;
  std::string pileupDT_;
  std::string lepton_;
  double par_;
  double par2_;
  bool usePartonFlavour_;
  bool pcut_;
  bool useDeltaR_;
  std::string path_;
  unsigned int icut_;

  JetCorrectionUncertainty *jetCorrectionUncertainty_;
  edm::LumiReWeighting LumiWeights_;

  table* ElSF_;
  table* ElSF2_;

  table* MuSF_;
  table* MuSF2_;
  table* BtSFT_;
  table* BtSFL_;
  table* LtSFT_;
  table* LtSFL_;
  table* BtEffT_;
  table* BtEffL_;
  table* CtEffT_;
  table* CtEffL_;
  table* LtEffT_;
  table* LtEffL_;

  int    n_ztautau[2];
  double w_ztautau[2];
  double w2_ztautau[2];

  int    n_events[2];
  double w_events[2];
  double w2_events[2];


  // ================================================================================================
  //  Histogram declaration


  // --- Z+jets

  TH1F* h_M;
  TH1F* h_M_b;
  TH1F* h_M_c;
  TH1F* h_M_l;
  TH1F* h_M_tt;

  TH1F* h_Pt;
  TH1F* h_Pt_b;
  TH1F* h_Pt_c;
  TH1F* h_Pt_l;
  TH1F* h_Pt_tt;

  TH1F* h_njet;
  TH1F* h_njet_b;
  TH1F* h_njet_c;
  TH1F* h_njet_l;
  TH1F* h_njet_tt;

  TH1F* h_MET;
  TH1F* h_MET_b;
  TH1F* h_MET_c;
  TH1F* h_MET_l;
  TH1F* h_MET_tt;

  TH1F* h_sMET;
  TH1F* h_sMET_b;
  TH1F* h_sMET_c;
  TH1F* h_sMET_l;
  TH1F* h_sMET_tt;

  TH1F* h_CSV;
  TH1F* h_CSV_b;
  TH1F* h_CSV_c;
  TH1F* h_CSV_l;
  TH1F* h_CSV_tt;

  TH1F* h_CSV_all;
  TH1F* h_CSV_all_b;
  TH1F* h_CSV_all_c;
  TH1F* h_CSV_all_l;
  TH1F* h_CSV_all_tt;

  TH1F* h_BJP;
  TH1F* h_BJP_b;
  TH1F* h_BJP_c;
  TH1F* h_BJP_l;
  TH1F* h_BJP_tt;

  TH2F* h_eff_b;
  TH2F* h_eff_c;
  TH2F* h_eff_l;


  // --- Z+c

  TH1F* hc_M;
  TH1F* hc_M_b;
  TH1F* hc_M_c;
  TH1F* hc_M_l;
  TH1F* hc_M_tt;

  TH1F* hc_Pt;
  TH1F* hc_Pt_b;
  TH1F* hc_Pt_c;
  TH1F* hc_Pt_l;
  TH1F* hc_Pt_tt;

  TH1F* hc_njet;
  TH1F* hc_njet_b;
  TH1F* hc_njet_c;
  TH1F* hc_njet_l;
  TH1F* hc_njet_tt;

  TH1F* hc_ncjet;
  TH1F* hc_ncjet_b;
  TH1F* hc_ncjet_c;
  TH1F* hc_ncjet_l;
  TH1F* hc_ncjet_tt;

  TH1F* hc_Pt_cjet1;
  TH1F* hc_Pt_cjet1_b;
  TH1F* hc_Pt_cjet1_c;
  TH1F* hc_Pt_cjet1_l;
  TH1F* hc_Pt_cjet1_tt;

  TH1F* hc_MET;
  TH1F* hc_MET_b;
  TH1F* hc_MET_c;
  TH1F* hc_MET_l;
  TH1F* hc_MET_tt;

  TH1F* hc_sMET;
  TH1F* hc_sMET_b;
  TH1F* hc_sMET_c;
  TH1F* hc_sMET_l;
  TH1F* hc_sMET_tt;

  TH1F* hc_CSV;
  TH1F* hc_CSV_b;
  TH1F* hc_CSV_c;
  TH1F* hc_CSV_l;
  TH1F* hc_CSV_tt;

  TH1F* hc_CSV_ctag;
  TH1F* hc_CSV_ctag_b;
  TH1F* hc_CSV_ctag_c;
  TH1F* hc_CSV_ctag_l;
  TH1F* hc_CSV_ctag_tt;

  TH1F* hc_CSV_all;
  TH1F* hc_CSV_all_b;
  TH1F* hc_CSV_all_c;
  TH1F* hc_CSV_all_l;
  TH1F* hc_CSV_all_tt;

  TH1F* hc_BJP;
  TH1F* hc_BJP_b;
  TH1F* hc_BJP_c;
  TH1F* hc_BJP_l;
  TH1F* hc_BJP_tt;

  TH1F* hc_BJP_ctag;
  TH1F* hc_BJP_ctag_b;
  TH1F* hc_BJP_ctag_c;
  TH1F* hc_BJP_ctag_l;
  TH1F* hc_BJP_ctag_tt;

  TH1F* hc_JPB;
  TH1F* hc_JPB_b;
  TH1F* hc_JPB_c;
  TH1F* hc_JPB_l;
  TH1F* hc_JPB_tt;

  TH1F* hc_CHP;
  TH1F* hc_CHP_b;
  TH1F* hc_CHP_c;
  TH1F* hc_CHP_l;
  TH1F* hc_CHP_tt;

  TH1F* hc_CHE;
  TH1F* hc_CHE_b;
  TH1F* hc_CHE_c;
  TH1F* hc_CHE_l;
  TH1F* hc_CHE_tt;

  TH1F* hc_svxM;
  TH1F* hc_svxM_b;
  TH1F* hc_svxM_c;
  TH1F* hc_svxM_l;
  TH1F* hc_svxM_tt;

  TH1F* hc_svxM_corr;
  TH1F* hc_svxM_corr_b;
  TH1F* hc_svxM_corr_c;
  TH1F* hc_svxM_corr_l;
  TH1F* hc_svxM_corr_tt;

  TH1F* hc_svxEfr;
  TH1F* hc_svxEfr_b;
  TH1F* hc_svxEfr_c;
  TH1F* hc_svxEfr_l;
  TH1F* hc_svxEfr_tt;

  TH2F* hc_CSVL_eff_b;
  TH2F* hc_CSVL_eff_c;
  TH2F* hc_CSVL_eff_l;
  TH2F* hc_CSVT_eff_b;
  TH2F* hc_CSVT_eff_c;
  TH2F* hc_CSVT_eff_l;

};

using namespace  pat;

//
// constants, enums and typedefs
//


//
// static data member definitions
//


//
// constructors and destructor
//
ZcAnalyzer::ZcAnalyzer (const edm::ParameterSet & iConfig) {

  pileupMC_ = iConfig.getUntrackedParameter < std::string > ("pileupMC", "S10");
  pileupDT_ = iConfig.getUntrackedParameter < std::string > ("pileupDT", "");
  lepton_ = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");
  par_ = iConfig.getUntrackedParameter <double> ("JEC", 0);
  par2_ = iConfig.getUntrackedParameter <double> ("JER", 0);
  path_ = iConfig.getUntrackedParameter < std::string > ("path", "/gpfs/cms/users/casarsa/analysis/Zc/work/SF");
  icut_ = iConfig.getUntrackedParameter <unsigned int> ("icut", 0);
  usePartonFlavour_ = iConfig.getUntrackedParameter <bool> ("usePartonFlavour", false);
  pcut_ = iConfig.getUntrackedParameter <bool> ("pcut", false);
  useDeltaR_ = iConfig.getUntrackedParameter <bool> ("useDeltaR", false);

  // now do what ever initialization is needed
  edm::Service < TFileService > fs;


  std::string channel = pileupDT_;
  if (  pileupDT_ == "mm" )
    channel = "#mu#mu"; 
  else if (  pileupDT_ == "em" )
    channel = "e#mu"; 


  // ------------------------------------------------------------------------------------------------
  //  Z+jets histos


  // dilepton mass
  std::string htitle = Form("%s - mass;M_{%s} [GeV/c^{2}]", channel.data(), channel.data());
  h_M    = fs->make < TH1F > ("h_M", htitle.data(), 50, 71., 111.);
  h_M    ->Sumw2();
  htitle = Form("%s - mass (b);M_{%s} [GeV/c^{2}]", channel.data(), channel.data());
  h_M_b  = fs->make < TH1F > ("h_M_b", htitle.data(), 50, 71., 111.);
  h_M_b  ->Sumw2();
  htitle = Form("%s - mass (c);M_{%s} [GeV/c^{2}]", channel.data(), channel.data());
  h_M_c  = fs->make < TH1F > ("h_M_c", htitle.data(), 50, 71., 111.);
  h_M_c  ->Sumw2();
  htitle = Form("%s - mass (dusg);M_{%s} [GeV/c^{2}]", channel.data(), channel.data());
  h_M_l  = fs->make < TH1F > ("h_M_l", htitle.data(), 50, 71., 111.);
  h_M_l  ->Sumw2();
  htitle = Form("%s - mass (Z#rightarrow#tau#tau);M_{%s} [GeV/c^{2}]", channel.data(), channel.data());
  h_M_tt = fs->make < TH1F > ("h_M_tt", htitle.data(), 50, 71., 111.);
  h_M_tt ->Sumw2();

  // dilepton PT
  htitle  = Form("%s - P_{T};P_{T}^{%s} [GeV/c]", channel.data(), channel.data());
  h_Pt    = fs->make < TH1F > ("h_Pt",    htitle.data(), 150, 0., 150.);
  h_Pt    ->Sumw2();
  htitle  = Form("%s - P_{T} (b);P_{T}^{%s} [GeV/c]", channel.data(), channel.data());
  h_Pt_b  = fs->make < TH1F > ("h_Pt_b",  htitle.data(), 150, 0., 150.);
  h_Pt_b  ->Sumw2();
  htitle  = Form("%s - P_{T} (c);P_{T}^{%s} [GeV/c]", channel.data(), channel.data());
  h_Pt_c  = fs->make < TH1F > ("h_Pt_c",  htitle.data(), 150, 0., 150.);
  h_Pt_c  ->Sumw2();
  htitle  = Form("%s - P_{T} (dusg);P_{T}^{%s} [GeV/c]", channel.data(), channel.data());
  h_Pt_l  = fs->make < TH1F > ("h_Pt_l",  htitle.data(), 150, 0., 150.);
  h_Pt_l  ->Sumw2();
  htitle  = Form("%s - P_{T} (Z#rightarrow#tau#tau);P_{T}^{%s} [GeV/c]", channel.data(), channel.data());
  h_Pt_tt = fs->make < TH1F > ("h_Pt_tt", htitle.data(), 150, 0., 150.);
  h_Pt_tt ->Sumw2();

  // jet multiplicity
  htitle    = Form("%s - jet multiplicity;N_{jet}", channel.data());
  h_njet    = fs->make < TH1F > ("h_njet",    htitle.data(), 10, 0.5, 10.5);
  h_njet    ->Sumw2();
  htitle    = Form("%s - jet multiplicity (b);N_{jet}", channel.data());
  h_njet_b  = fs->make < TH1F > ("h_njet_b",  htitle.data(), 10, 0.5, 10.5);
  h_njet_b  ->Sumw2();
  htitle    = Form("%s - jet multiplicity (c);N_{jet}", channel.data());
  h_njet_c  = fs->make < TH1F > ("h_njet_c",  htitle.data(), 10, 0.5, 10.5);
  h_njet_c  ->Sumw2();
  htitle    = Form("%s - jet multiplicity (dusg);N_{jet}", channel.data());
  h_njet_l  = fs->make < TH1F > ("h_njet_l",  htitle.data(), 10, 0.5, 10.5);
  h_njet_l  ->Sumw2();
  htitle    = Form("%s - jet multiplicity (Z#rightarrow#tau#tau);N_{jet}", channel.data());
  h_njet_tt = fs->make < TH1F > ("h_njet_tt", htitle.data(), 10, 0.5, 10.5);
  h_njet_tt ->Sumw2();

  // MET
  htitle   = Form("%s - MET;MET [GeV]", channel.data());
  h_MET    = fs->make < TH1F > ("h_MET",    htitle.data(), 100, 0., 250.);
  h_MET    ->Sumw2();
  htitle   = Form("%s - MET (b);MET [GeV]", channel.data());
  h_MET_b  = fs->make < TH1F > ("h_MET_b ", htitle.data(), 100, 0., 250.);
  h_MET_b  ->Sumw2();
  htitle   = Form("%s - MET (c);MET [GeV]", channel.data());
  h_MET_c  = fs->make < TH1F > ("h_MET_c ", htitle.data(), 100, 0., 250.);
  h_MET_c  ->Sumw2();
  htitle   = Form("%s - MET (dusg);MET [GeV]", channel.data());
  h_MET_l  = fs->make < TH1F > ("h_MET_l ", htitle.data(), 100, 0., 250.);
  h_MET_l  ->Sumw2();
  htitle   = Form("%s - MET (Z#rightarrow#tau#tau);MET [GeV]", channel.data());
  h_MET_tt = fs->make < TH1F > ("h_MET_tt", htitle.data(), 100, 0., 250.);
  h_MET_tt ->Sumw2();

  // MET significance
  htitle    = Form("%s - MET significance;MET significance", channel.data());
  h_sMET    = fs->make < TH1F > ("h_sMET",    htitle.data(), 100, 0., 250.);
  h_sMET    ->Sumw2();
  htitle    = Form("%s - MET significance (b);MET significance", channel.data());
  h_sMET_b  = fs->make < TH1F > ("h_sMET_b ", htitle.data(), 100, 0., 250.);
  h_sMET_b  ->Sumw2();
  htitle    = Form("%s - MET significance (c);MET significance", channel.data());
  h_sMET_c  = fs->make < TH1F > ("h_sMET_c ", htitle.data(), 100, 0., 250.);
  h_sMET_c  ->Sumw2();
  htitle    = Form("%s - MET significance (dusg);MET significance", channel.data());
  h_sMET_l  = fs->make < TH1F > ("h_sMET_l ", htitle.data(), 100, 0., 250.);
  h_sMET_l  ->Sumw2();
  htitle    = Form("%s - MET significance (Z#rightarrow#tau#tau);MET significance", channel.data());
  h_sMET_tt = fs->make < TH1F > ("h_sMET_tt", htitle.data(), 100, 0., 250.);
  h_sMET_tt ->Sumw2();

  // CSV discriminant
  htitle   = Form("%s - CSV discriminant;CSV", channel.data());
  h_CSV    = fs->make < TH1F > ("h_CSV", htitle.data(), 100, 0., 1.);
  h_CSV    ->Sumw2();
  htitle   = Form("%s - CSV discriminant (b);CSV", channel.data());
  h_CSV_b  = fs->make < TH1F > ("h_CSV_b", htitle.data(), 100, 0., 1.);
  h_CSV_b  ->Sumw2();
  htitle   = Form("%s - CSV discriminant (c);CSV", channel.data());
  h_CSV_c  = fs->make < TH1F > ("h_CSV_c", htitle.data(), 100, 0., 1.);
  h_CSV_c  ->Sumw2();
  htitle   = Form("%s - CSV discriminant (dusg);CSV", channel.data());
  h_CSV_l  = fs->make < TH1F > ("h_CSV_l", htitle.data(), 100, 0., 1.);
  h_CSV_l  ->Sumw2();
  htitle   = Form("%s - CSV discriminant (Z#rightarrow#tau#tau);CSV", channel.data());
  h_CSV_tt = fs->make < TH1F > ("h_CSV_tt", htitle.data(), 100, 0., 1.);
  h_CSV_tt ->Sumw2();

  htitle       = Form("%s (all jets) - CSV discriminant;CSV", channel.data());
  h_CSV_all    = fs->make < TH1F > ("h_CSV_all", htitle.data(), 100, 0., 1.);
  h_CSV_all    ->Sumw2();
  htitle       = Form("%s (all jets) - CSV discriminant (b);CSV", channel.data());
  h_CSV_all_b  = fs->make < TH1F > ("h_CSV_all_b", htitle.data(), 100, 0., 1.);
  h_CSV_all_b  ->Sumw2();
  htitle       = Form("%s (all jets) - CSV discriminant (c);CSV", channel.data());
  h_CSV_all_c  = fs->make < TH1F > ("h_CSV_all_c", htitle.data(), 100, 0., 1.);
  h_CSV_all_c  ->Sumw2();
  htitle       = Form("%s (all jets) - CSV discriminant (dusg);CSV", channel.data());
  h_CSV_all_l  = fs->make < TH1F > ("h_CSV_all_l", htitle.data(), 100, 0., 1.);
  h_CSV_all_l  ->Sumw2();
  htitle       = Form("%s (all jets) - CSV discriminant (Z#rightarrow#tau#tau);CSV", channel.data());
  h_CSV_all_tt = fs->make < TH1F > ("h_CSV_all_tt", htitle.data(), 100, 0., 1.);
  h_CSV_all_tt ->Sumw2();

  // BJP discriminant
  htitle   = Form("%s - BJP discriminant;BJP", channel.data());
  h_BJP    = fs->make < TH1F > ("h_BJP", htitle.data(), 100, 0., 10.);
  h_BJP    ->Sumw2();
  htitle   = Form("%s - BJP discriminant (b);BJP", channel.data());
  h_BJP_b  = fs->make < TH1F > ("h_BJP_b", htitle.data(), 100, 0., 10.);
  h_BJP_b  ->Sumw2();
  htitle   = Form("%s - BJP discriminant (c);BJP", channel.data());
  h_BJP_c  = fs->make < TH1F > ("h_BJP_c", htitle.data(), 100, 0., 10.);
  h_BJP_c  ->Sumw2();
  htitle   = Form("%s - BJP discriminant (dusg);BJP", channel.data());
  h_BJP_l  = fs->make < TH1F > ("h_BJP_l", htitle.data(), 100, 0., 10.);
  h_BJP_l  ->Sumw2();
  htitle   = Form("%s - BJP discriminant (Z#rightarrow#tau#tau);BJP", channel.data());
  h_BJP_tt = fs->make < TH1F > ("h_BJP_tt", htitle.data(), 100, 0., 10.);
  h_BJP_tt ->Sumw2();

  // counters for b-tagging efficiency calculation
  const float xbins[17] = {20., 30., 40., 50., 60., 70., 80., 100., 120., 160., 210., 260., 320., 400., 500., 600., 800.};
  const float ybins[9]  = {-2.4, -1.5, -1.0, -0.5, 0., 0.5, 1.0, 1.5, 2.4};
  htitle   = Form("%s - before tagging (b)", channel.data());
  h_eff_b = fs->make < TH2F > ("h_eff_b", htitle.data(), 16, xbins, 8, ybins);
  h_eff_b->Sumw2();
  htitle   = Form("%s - before tagging (l)", channel.data());
  h_eff_c = fs->make < TH2F > ("h_eff_c", htitle.data(), 16, xbins, 8, ybins);
  h_eff_c->Sumw2();
  htitle   = Form("%s - before tagging (dusg)", channel.data());
  h_eff_l = fs->make < TH2F > ("h_eff_l", htitle.data(), 16, xbins, 8, ybins);
  h_eff_l->Sumw2();



  // ------------------------------------------------------------------------------------------------
  //  Z+c histos


  // dilepton mass
  htitle = Form("%s - mass;M_{%s} [GeV/c^{2}]", channel.data(), channel.data());
  hc_M    = fs->make < TH1F > ("hc_M", htitle.data(), 50, 71., 111.);
  hc_M    ->Sumw2();
  htitle  = Form("%s - mass (b);M_{%s} [GeV/c^{2}]", channel.data(), channel.data());
  hc_M_b  = fs->make < TH1F > ("hc_M_b", htitle.data(), 50, 71., 111.);
  hc_M_b  ->Sumw2();
  htitle  = Form("%s - mass (c);M_{%s} [GeV/c^{2}]", channel.data(), channel.data());
  hc_M_c  = fs->make < TH1F > ("hc_M_c", htitle.data(), 50, 71., 111.);
  hc_M_c  ->Sumw2();
  htitle  = Form("%s - mass (dusg);M_{%s} [GeV/c^{2}]", channel.data(), channel.data());
  hc_M_l  = fs->make < TH1F > ("hc_M_l", htitle.data(), 50, 71., 111.);
  hc_M_l  ->Sumw2();
  htitle  = Form("%s - mass (Z#rightarrow#tau#tau);M_{%s} [GeV/c^{2}]", channel.data(), channel.data());
  hc_M_tt = fs->make < TH1F > ("hc_M_tt", htitle.data(), 50, 71., 111.);
  hc_M_tt ->Sumw2();

  // dilepton PT
  htitle   = Form("%s - P_{T};P_{T}^{%s} [GeV/c]", channel.data(), channel.data());
  hc_Pt    = fs->make < TH1F > ("hc_Pt",    htitle.data(), 150, 0., 150.);
  hc_Pt    ->Sumw2();
  htitle   = Form("%s - P_{T} (b);P_{T}^{%s} [GeV/c]", channel.data(), channel.data());
  hc_Pt_b  = fs->make < TH1F > ("hc_Pt_b",  htitle.data(), 150, 0., 150.);
  hc_Pt_b  ->Sumw2();
  htitle   = Form("%s - P_{T} (c);P_{T}^{%s} [GeV/c]", channel.data(), channel.data());
  hc_Pt_c  = fs->make < TH1F > ("hc_Pt_c",  htitle.data(), 150, 0., 150.);
  hc_Pt_c  ->Sumw2();
  htitle   = Form("%s - P_{T} (dusg);P_{T}^{%s} [GeV/c]", channel.data(), channel.data());
  hc_Pt_l  = fs->make < TH1F > ("hc_Pt_l",  htitle.data(), 150, 0., 150.);
  hc_Pt_l  ->Sumw2();
  htitle   = Form("%s - P_{T} (Z#rightarrow#tau#tau);P_{T}^{%s} [GeV/c]", channel.data(), channel.data());
  hc_Pt_tt = fs->make < TH1F > ("hc_Pt_tt", htitle.data(), 150, 0., 150.);
  hc_Pt_tt ->Sumw2();

  // jet multiplicity
  htitle     = Form("%s - jet multiplicity;N_{jet}", channel.data());
  hc_njet    = fs->make < TH1F > ("hc_njet",    htitle.data(), 10, 0.5, 10.5);
  hc_njet    ->Sumw2();
  htitle     = Form("%s - jet multiplicity (b);N_{jet}", channel.data());
  hc_njet_b  = fs->make < TH1F > ("hc_njet_b",  htitle.data(), 10, 0.5, 10.5);
  hc_njet_b  ->Sumw2();
  htitle     = Form("%s - jet multiplicity (c);N_{jet}", channel.data());
  hc_njet_c  = fs->make < TH1F > ("hc_njet_c",  htitle.data(), 10, 0.5, 10.5);
  hc_njet_c  ->Sumw2();
  htitle     = Form("%s - jet multiplicity (dusg);N_{jet}", channel.data());
  hc_njet_l  = fs->make < TH1F > ("hc_njet_l",  htitle.data(), 10, 0.5, 10.5);
  hc_njet_l  ->Sumw2();
  htitle     = Form("%s - jet multiplicity (Z#rightarrow#tau#tau);N_{jet}", channel.data());
  hc_njet_tt = fs->make < TH1F > ("hc_njet_tt", htitle.data(), 10, 0.5, 10.5);
  hc_njet_tt ->Sumw2();

  // c-tagged jet multiplicity
  htitle      = Form("%s - c-jet multiplicity;N_{c-jet}", channel.data());
  hc_ncjet    = fs->make < TH1F > ("hc_ncjet",    htitle.data(), 10, 0.5, 10.5);
  hc_ncjet    ->Sumw2();
  htitle      = Form("%s - c-jet multiplicity (b);N_{c-jet}", channel.data());
  hc_ncjet_b  = fs->make < TH1F > ("hc_ncjet_b",  htitle.data(), 10, 0.5, 10.5);
  hc_ncjet_b  ->Sumw2();
  htitle      = Form("%s - c-jet multiplicity (c);N_{c-jet}", channel.data());
  hc_ncjet_c  = fs->make < TH1F > ("hc_ncjet_c",  htitle.data(), 10, 0.5, 10.5);
  hc_ncjet_c  ->Sumw2();
  htitle      = Form("%s - c-jet multiplicity (dusg);N_{c-jet}", channel.data());
  hc_ncjet_l  = fs->make < TH1F > ("hc_ncjet_l",  htitle.data(), 10, 0.5, 10.5);
  hc_ncjet_l  ->Sumw2();
  htitle      = Form("%s - c-jet multiplicity (Z#rightarrow#tau#tau);N_{c-jet}", channel.data());
  hc_ncjet_tt = fs->make < TH1F > ("hc_ncjet_tt", htitle.data(), 10, 0.5, 10.5);
  hc_ncjet_tt ->Sumw2();

  // PT of the leading c-tagged jet
  htitle         = Form("%s - P_{T};P_{T}^{jet} [GeV/c]", channel.data());
  hc_Pt_cjet1    = fs->make < TH1F > ("hc_Pt_cjet1", htitle.data(), 500, 0., 2000.);
  hc_Pt_cjet1    ->Sumw2();
  htitle         = Form("%s - P_{T} (b);P_{T}^{jet} [GeV/c]", channel.data());
  hc_Pt_cjet1_b  = fs->make < TH1F > ("hc_Pt_cjet1_b", htitle.data(), 500, 0., 2000.);
  hc_Pt_cjet1_b  ->Sumw2();
  htitle         = Form("%s - P_{T} (c);P_{T}^{jet} [GeV/c]", channel.data());
  hc_Pt_cjet1_c  = fs->make < TH1F > ("hc_Pt_cjet1_c", htitle.data(), 500, 0., 2000.);
  hc_Pt_cjet1_c  ->Sumw2();
  htitle         = Form("%s - P_{T} (dusg);P_{T}^{jet} [GeV/c]", channel.data());
  hc_Pt_cjet1_l  = fs->make < TH1F > ("hc_Pt_cjet1_l", htitle.data(), 500, 0., 2000.);
  hc_Pt_cjet1_l  ->Sumw2();
  htitle         = Form("%s - P_{T} (Z#rightarrow#tau#tau);P_{T}^{jet} [GeV/c]", channel.data());
  hc_Pt_cjet1_tt = fs->make < TH1F > ("hc_Pt_cjet1_tt", htitle.data(), 500, 0., 2000.);
  hc_Pt_cjet1_tt ->Sumw2();

  // MET
  htitle    = Form("%s - MET;MET [GeV]", channel.data());
  hc_MET    = fs->make < TH1F > ("hc_MET",    htitle.data(), 100, 0., 250.);
  hc_MET    ->Sumw2();
  htitle    = Form("%s - MET (b);MET [GeV]", channel.data());
  hc_MET_b  = fs->make < TH1F > ("hc_MET_b ", htitle.data(), 100, 0., 250.);
  hc_MET_b  ->Sumw2();
  htitle    = Form("%s - MET (c);MET [GeV]", channel.data());
  hc_MET_c  = fs->make < TH1F > ("hc_MET_c ", htitle.data(), 100, 0., 250.);
  hc_MET_c  ->Sumw2();
  htitle    = Form("%s - MET (dusg);MET [GeV]", channel.data());
  hc_MET_l  = fs->make < TH1F > ("hc_MET_l ", htitle.data(), 100, 0., 250.);
  hc_MET_l  ->Sumw2();
  htitle    = Form("%s - MET (Z#rightarrow#tau#tau);MET [GeV]", channel.data());
  hc_MET_tt = fs->make < TH1F > ("hc_MET_tt", htitle.data(), 100, 0., 250.);
  hc_MET_tt ->Sumw2();

  // MET significance
  htitle     = Form("%s - MET significance;MET significance", channel.data());
  hc_sMET    = fs->make < TH1F > ("hc_sMET",    htitle.data(), 100, 0., 250.);
  hc_sMET    ->Sumw2();
  htitle     = Form("%s - MET significance (b);MET significance", channel.data());
  hc_sMET_b  = fs->make < TH1F > ("hc_sMET_b ", htitle.data(), 100, 0., 250.);
  hc_sMET_b  ->Sumw2();
  htitle     = Form("%s - MET significance (c);MET significance", channel.data());
  hc_sMET_c  = fs->make < TH1F > ("hc_sMET_c ", htitle.data(), 100, 0., 250.);
  hc_sMET_c  ->Sumw2();
  htitle     = Form("%s - MET significance (dusg);MET significance", channel.data());
  hc_sMET_l  = fs->make < TH1F > ("hc_sMET_l ", htitle.data(), 100, 0., 250.);
  hc_sMET_l  ->Sumw2();
  htitle     = Form("%s - MET significance (Z#rightarrow#tau#tau);MET significance", channel.data());
  hc_sMET_tt = fs->make < TH1F > ("hc_sMET_tt", htitle.data(), 100, 0., 250.);
  hc_sMET_tt ->Sumw2();

  // CSV discriminant
  htitle    = Form("%s - CSV discriminant;CSV", channel.data());
  hc_CSV    = fs->make < TH1F > ("hc_CSV", htitle.data(), 100, 0., 1.);
  hc_CSV    ->Sumw2();
  htitle    = Form("%s - CSV discriminant (b);CSV", channel.data());
  hc_CSV_b  = fs->make < TH1F > ("hc_CSV_b", htitle.data(), 100, 0., 1.);
  hc_CSV_b  ->Sumw2();
  htitle    = Form("%s - CSV discriminant (c);CSV", channel.data());
  hc_CSV_c  = fs->make < TH1F > ("hc_CSV_c", htitle.data(), 100, 0., 1.);
  hc_CSV_c  ->Sumw2();
  htitle    = Form("%s - CSV discriminant (dusg);CSV", channel.data());
  hc_CSV_l  = fs->make < TH1F > ("hc_CSV_l", htitle.data(), 100, 0., 1.);
  hc_CSV_l  ->Sumw2();
  htitle    = Form("%s - CSV discriminant (Z#rightarrow#tau#tau);CSV", channel.data());
  hc_CSV_tt = fs->make < TH1F > ("hc_CSV_tt", htitle.data(), 100, 0., 1.);
  hc_CSV_tt ->Sumw2();

  htitle        = Form("%s (all jets) - CSV discriminant;CSV", channel.data());
  hc_CSV_all    = fs->make < TH1F > ("hc_CSV_all", htitle.data(), 100, 0., 1.);
  hc_CSV_all    ->Sumw2();
  htitle        = Form("%s (all jets) - CSV discriminant (b);CSV", channel.data());
  hc_CSV_all_b  = fs->make < TH1F > ("hc_CSV_all_b", htitle.data(), 100, 0., 1.);
  hc_CSV_all_b  ->Sumw2();
  htitle        = Form("%s (all jets) - CSV discriminant (c);CSV", channel.data());
  hc_CSV_all_c  = fs->make < TH1F > ("hc_CSV_all_c", htitle.data(), 100, 0., 1.);
  hc_CSV_all_c  ->Sumw2();
  htitle        = Form("%s (all jets) - CSV discriminant (dusg);CSV", channel.data());
  hc_CSV_all_l  = fs->make < TH1F > ("hc_CSV_all_l", htitle.data(), 100, 0., 1.);
  hc_CSV_all_l  ->Sumw2();
  htitle        = Form("%s (all jets) - CSV discriminant (Z#rightarrow#tau#tau);CSV", channel.data());
  hc_CSV_all_tt = fs->make < TH1F > ("hc_CSV_all_tt", htitle.data(), 100, 0., 1.);
  hc_CSV_all_tt ->Sumw2();

  htitle         = Form("%s (c-tagget jets) - CSV discriminant;CSV", channel.data());
  hc_CSV_ctag    = fs->make < TH1F > ("hc_CSV_ctag", htitle.data(), 100, 0., 1.);
  hc_CSV_ctag    ->Sumw2();
  htitle         = Form("%s (c-tagged jets) - CSV discriminant (b);CSV", channel.data());
  hc_CSV_ctag_b  = fs->make < TH1F > ("hc_CSV_ctag_b", htitle.data(), 100, 0., 1.);
  hc_CSV_ctag_b  ->Sumw2();
  htitle         = Form("%s (c-tagged jets) - CSV discriminant (c);CSV", channel.data());
  hc_CSV_ctag_c  = fs->make < TH1F > ("hc_CSV_ctag_c", htitle.data(), 100, 0., 1.);
  hc_CSV_ctag_c  ->Sumw2();
  htitle         = Form("%s (c-tagged jets) - CSV discriminant (dusg);CSV", channel.data());
  hc_CSV_ctag_l  = fs->make < TH1F > ("hc_CSV_ctag_l", htitle.data(), 100, 0., 1.);
  hc_CSV_ctag_l  ->Sumw2();
  htitle         = Form("%s (c-tagged jets) - CSV discriminant (Z#rightarrow#tau#tau);CSV", channel.data());
  hc_CSV_ctag_tt = fs->make < TH1F > ("hc_CSV_ctag_tt", htitle.data(), 100, 0., 1.);
  hc_CSV_ctag_tt ->Sumw2();

  // BJP discriminant
  htitle    = Form("%s - BJP discriminant;BJP", channel.data());
  hc_BJP    = fs->make < TH1F > ("hc_BJP", htitle.data(), 100, 0., 10.);
  hc_BJP    ->Sumw2();
  htitle    = Form("%s - BJP discriminant (b);BJP", channel.data());
  hc_BJP_b  = fs->make < TH1F > ("hc_BJP_b", htitle.data(), 100, 0., 10.);
  hc_BJP_b  ->Sumw2();
  htitle    = Form("%s - BJP discriminant (c);BJP", channel.data());
  hc_BJP_c  = fs->make < TH1F > ("hc_BJP_c", htitle.data(), 100, 0., 10.);
  hc_BJP_c  ->Sumw2();
  htitle    = Form("%s - BJP discriminant (dusg);BJP", channel.data());
  hc_BJP_l  = fs->make < TH1F > ("hc_BJP_l", htitle.data(), 100, 0., 10.);
  hc_BJP_l  ->Sumw2();
  htitle    = Form("%s - BJP discriminant (Z#rightarrow#tau#tau);BJP", channel.data());
  hc_BJP_tt = fs->make < TH1F > ("hc_BJP_tt", htitle.data(), 100, 0., 10.);
  hc_BJP_tt ->Sumw2();

  htitle         = Form("%s (c-tagget jets) - BJP discriminant;BJP", channel.data());
  hc_BJP_ctag    = fs->make < TH1F > ("hc_BJP_ctag", htitle.data(), 100, 0., 10.);
  hc_BJP_ctag    ->Sumw2();
  htitle         = Form("%s (c-tagged jets) - BJP discriminant (b);BJP", channel.data());
  hc_BJP_ctag_b  = fs->make < TH1F > ("hc_BJP_ctag_b", htitle.data(), 100, 0., 10.);
  hc_BJP_ctag_b  ->Sumw2();
  htitle         = Form("%s (c-tagged jets) - BJP discriminant (c);BJP", channel.data());
  hc_BJP_ctag_c  = fs->make < TH1F > ("hc_BJP_ctag_c", htitle.data(), 100, 0., 10.);
  hc_BJP_ctag_c  ->Sumw2();
  htitle         = Form("%s (c-tagged jets) - BJP discriminant (dusg);BJP", channel.data());
  hc_BJP_ctag_l  = fs->make < TH1F > ("hc_BJP_ctag_l", htitle.data(), 100, 0., 10.);
  hc_BJP_ctag_l  ->Sumw2();
  htitle         = Form("%s (c-tagged jets) - BJP discriminant (Z#rightarrow#tau#tau);BJP", channel.data());
  hc_BJP_ctag_tt = fs->make < TH1F > ("hc_BJP_ctag_tt", htitle.data(), 100, 0., 10.);
  hc_BJP_ctag_tt ->Sumw2();

  // JBP discriminant
  htitle    = Form("%s - JPB discriminant;JPB", channel.data());
  hc_JPB    = fs->make < TH1F > ("hc_JPB", htitle.data(), 100, 0., 2.);
  hc_JPB    ->Sumw2();
  htitle    = Form("%s - JPB discriminant (b);JPB", channel.data());
  hc_JPB_b  = fs->make < TH1F > ("hc_JPB_b", htitle.data(), 100, 0., 2.);
  hc_JPB_b  ->Sumw2();
  htitle    = Form("%s - JPB discriminant (c);JPB", channel.data());
  hc_JPB_c  = fs->make < TH1F > ("hc_JPB_c", htitle.data(), 100, 0., 2.);
  hc_JPB_c  ->Sumw2();
  htitle    = Form("%s - JPB discriminant (dusg);JPB", channel.data());
  hc_JPB_l  = fs->make < TH1F > ("hc_JPB_l", htitle.data(), 100, 0., 2.);
  hc_JPB_l  ->Sumw2();
  htitle    = Form("%s - JPB discriminant (Z#rightarrow#tau#tau);JPB", channel.data());
  hc_JPB_tt = fs->make < TH1F > ("hc_JPB_tt", htitle.data(), 100, 0., 2.);
  hc_JPB_tt ->Sumw2();

  // CHP discriminant
  htitle    = Form("%s - CHP discriminant;CHP", channel.data());
  hc_CHP    = fs->make < TH1F > ("hc_CHP", htitle.data(), 100, 0., 20.);
  hc_CHP    ->Sumw2();
  htitle    = Form("%s - CHP discriminant (b);CHP", channel.data());
  hc_CHP_b  = fs->make < TH1F > ("hc_CHP_b", htitle.data(), 100, 0., 20.);
  hc_CHP_b  ->Sumw2();
  htitle    = Form("%s - CHP discriminant (c);CHP", channel.data());
  hc_CHP_c  = fs->make < TH1F > ("hc_CHP_c", htitle.data(), 100, 0., 20.);
  hc_CHP_c  ->Sumw2();
  htitle    = Form("%s - CHP discriminant (dusg);CHP", channel.data());
  hc_CHP_l  = fs->make < TH1F > ("hc_CHP_l", htitle.data(), 100, 0., 20.);
  hc_CHP_l  ->Sumw2();
  htitle    = Form("%s - CHP discriminant (Z#rightarrow#tau#tau);CHP", channel.data());
  hc_CHP_tt = fs->make < TH1F > ("hc_CHP_tt", htitle.data(), 100, 0., 20.);
  hc_CHP_tt ->Sumw2();

  // CHE discriminant
  htitle    = Form("%s - CHE discriminant;CHE", channel.data());
  hc_CHE    = fs->make < TH1F > ("hc_CHE", htitle.data(), 100, 0., 20.);
  hc_CHE    ->Sumw2();
  htitle    = Form("%s - CHE discriminant (b);CHE", channel.data());
  hc_CHE_b  = fs->make < TH1F > ("hc_CHE_b", htitle.data(), 100, 0., 20.);
  hc_CHE_b  ->Sumw2();
  htitle    = Form("%s - CHE discriminant (c);CHE", channel.data());
  hc_CHE_c  = fs->make < TH1F > ("hc_CHE_c", htitle.data(), 100, 0., 20.);
  hc_CHE_c  ->Sumw2();
  htitle    = Form("%s - CHE discriminant (dusg);CHE", channel.data());
  hc_CHE_l  = fs->make < TH1F > ("hc_CHE_l", htitle.data(), 100, 0., 20.);
  hc_CHE_l  ->Sumw2();
  htitle    = Form("%s - CHE discriminant (Z#rightarrow#tau#tau);CHE", channel.data());
  hc_CHE_tt = fs->make < TH1F > ("hc_CHE_tt", htitle.data(), 100, 0., 20.);
  hc_CHE_tt ->Sumw2();

  // Secondary VTX mass
  htitle     = Form("%s - Secondary VTX mass;M_{SVX} [GeV/c^{2}]", channel.data());
  hc_svxM    = fs->make < TH1F > ("hc_svxM",    htitle.data(), 200, 0., 10.);
  hc_svxM    ->Sumw2();
  htitle     = Form("%s - Secondary VTX mass (b);M_{SVX} [GeV/c^{2}]", channel.data());
  hc_svxM_b  = fs->make < TH1F > ("hc_svxM_b",  htitle.data(), 200, 0., 10.);
  hc_svxM_b  ->Sumw2();
  htitle     = Form("%s - Secondary VTX mass (c);M_{SVX} [GeV/c^{2}]", channel.data());
  hc_svxM_c  = fs->make < TH1F > ("hc_svxM_c",  htitle.data(), 200, 0., 10.);
  hc_svxM_c  ->Sumw2();
  htitle     = Form("%s - Secondary VTX mass (dusg);M_{SVX} [GeV/c^{2}]", channel.data());
  hc_svxM_l  = fs->make < TH1F > ("hc_svxM_l",  htitle.data(), 200, 0., 10.);
  hc_svxM_l  ->Sumw2();
  htitle     = Form("%s - Secondary VTX mass (Z#rightarrow#tau#tau);M_{SVX} [GeV/c^{2}]", channel.data());
  hc_svxM_tt = fs->make < TH1F > ("hc_svxM_tt", htitle.data(), 200, 0., 10.);
  hc_svxM_tt ->Sumw2();
  
  htitle          = Form("%s - Secondary VTX mass;M_{SVX} [GeV/c^{2}]", channel.data());
  hc_svxM_corr    = fs->make < TH1F > ("hc_svxM_corr",    htitle.data(), 200, 0., 10.);
  hc_svxM_corr    ->Sumw2();
  htitle          = Form("%s - Secondary VTX mass (b);M_{SVX} [GeV/c^{2}]", channel.data());
  hc_svxM_corr_b  = fs->make < TH1F > ("hc_svxM_corr_b",  htitle.data(), 200, 0., 10.);
  hc_svxM_corr_b  ->Sumw2();
  htitle          = Form("%s - Secondary VTX mass (c);M_{SVX} [GeV/c^{2}]", channel.data());
  hc_svxM_corr_c  = fs->make < TH1F > ("hc_svxM_corr_c",  htitle.data(), 200, 0., 10.);
  hc_svxM_corr_c  ->Sumw2();
  htitle          = Form("%s - Secondary VTX mass (dusg);M_{SVX} [GeV/c^{2}]", channel.data());
  hc_svxM_corr_l  = fs->make < TH1F > ("hc_svxM_corr_l",  htitle.data(), 200, 0., 10.);
  hc_svxM_corr_l  ->Sumw2();
  htitle          = Form("%s - Secondary VTX mass (Z#rightarrow#tau#tau);M_{SVX} [GeV/c^{2}]", channel.data());
  hc_svxM_corr_tt = fs->make < TH1F > ("hc_svxM_corr_tt", htitle.data(), 200, 0., 10.);
  hc_svxM_corr_tt ->Sumw2();
  
  // Secondary VTX energy fraction
  htitle       = Form("%s - Sec. VTX energy fraction;SVX energy fraction", channel.data());
  hc_svxEfr    = fs->make < TH1F > ("hc_svxEfr",    htitle.data(), 100, 0., 1.);
  hc_svxEfr    ->Sumw2();
  htitle       = Form("%s - Sec. VTX energy fraction (b);SVX energy fraction", channel.data());
  hc_svxEfr_b  = fs->make < TH1F > ("hc_svxEfr_b",  htitle.data(), 100, 0., 1.);
  hc_svxEfr_b  ->Sumw2();
  htitle       = Form("%s - Sec. VTX energy fraction (c);SVX energy fraction", channel.data());
  hc_svxEfr_c  = fs->make < TH1F > ("hc_svxEfr_c",  htitle.data(), 100, 0., 1.);
  hc_svxEfr_c  ->Sumw2();
  htitle       = Form("%s - Sec. VTX energy fraction (dusg);SVX energy fraction", channel.data());
  hc_svxEfr_l  = fs->make < TH1F > ("hc_svxEfr_l",  htitle.data(), 100, 0., 1.);
  hc_svxEfr_l  ->Sumw2();
  htitle       = Form("%s - Sec. VTX energy fraction (Z#rightarrow#tau#tau);SVX energy fraction", channel.data());
  hc_svxEfr_tt = fs->make < TH1F > ("hc_svxEfr_tt", htitle.data(), 100, 0., 1.);
  hc_svxEfr_tt ->Sumw2();
  
  // counters for b-tagging efficiency calculation
  htitle        = Form("%s - after CSVL tagging (b)", channel.data());
  hc_CSVL_eff_b = fs->make < TH2F > ("hc_CSVL_eff_b", htitle.data(), 16, xbins, 8, ybins);
  hc_CSVL_eff_b ->Sumw2();
  htitle        = Form("%s - after CSVL tagging (c)", channel.data());
  hc_CSVL_eff_c = fs->make < TH2F > ("hc_CSVL_eff_c", htitle.data(), 16, xbins, 8, ybins);
  hc_CSVL_eff_c ->Sumw2();
  htitle        = Form("%s - after CSVL tagging (dusg)", channel.data());
  hc_CSVL_eff_l = fs->make < TH2F > ("hc_CSVL_eff_l", htitle.data(), 16, xbins, 8, ybins);
  hc_CSVL_eff_l ->Sumw2();

  htitle        = Form("%s - after CSVT tagging (b)", channel.data());
  hc_CSVT_eff_b = fs->make < TH2F > ("hc_CSVT_eff_b", htitle.data(), 16, xbins, 8, ybins);
  hc_CSVT_eff_b ->Sumw2();
  htitle        = Form("%s - after CSVT tagging (c)", channel.data());
  hc_CSVT_eff_c = fs->make < TH2F > ("hc_CSVT_eff_c", htitle.data(), 16, xbins, 8, ybins);
  hc_CSVT_eff_c ->Sumw2();
  htitle        = Form("%s - after CSVT tagging (dusg)", channel.data());
  hc_CSVT_eff_l = fs->make < TH2F > ("hc_CSVT_eff_l", htitle.data(), 16, xbins, 8, ybins);
  hc_CSVT_eff_l ->Sumw2();


}

ZcAnalyzer::~ZcAnalyzer () {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event ------------
void ZcAnalyzer::produce (edm::Event & iEvent, const edm::EventSetup & iSetup) {


  using namespace edm;
  using namespace std;

  bool isMC = false;

  bool ee_event = false;
  bool mm_event = false;
  bool em_event = false;

  bool is_tautau = false;
  unsigned int n_b = 0;
  unsigned int n_c = 0;

  double MyWeight = 1.;



  // ------------------------------------------------------------------------------------------------
  //  Access event information

  // Get electron collection
  edm::Handle < pat::ElectronCollection > electrons;
  iEvent.getByLabel ("matchedElectrons", electrons);
  if (lepton_=="electron+muon") iEvent.getByLabel ("matchedElectronsEM", electrons);


  // Get muon collection
  edm::Handle < pat::MuonCollection > muons;
  iEvent.getByLabel ("matchedMuons", muons);
  if (lepton_=="electron+muon") iEvent.getByLabel ("matchedMuonsEM", muons);


  // Get jet collection
  edm::Handle < vector < pat::Jet > > jets;
  iEvent.getByLabel ("goodJets", jets);


  // Get MET
  edm::Handle < vector < pat::MET > > met;
  iEvent.getByLabel (edm::InputTag ("patMETsPFlow"), met);


  // Get GEN particle collection
  edm::Handle<vector<reco::GenParticle> > genPart;
  iEvent.getByLabel ("genParticles", genPart);


  if (lepton_=="electron" && !electrons.isValid()) return;
  if (lepton_=="muon" && !muons.isValid()) return;
  if (lepton_=="electron+muon" && !electrons.isValid()) return;
  if (lepton_=="electron+muon" && !muons.isValid()) return;

  

  // ------------------------------------------------------------------------------------------------
  //  Pile-up reweighting


  edm::Handle < vector < PileupSummaryInfo > > PupInfo;

  if (iEvent.getByLabel (edm::InputTag ("addPileupInfo"), PupInfo)) {

    isMC = true;

    float Tnpv = -1.;

    for (vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI!=PupInfo->end(); ++PVI) {
      int BX = PVI->getBunchCrossing ();
      if (BX == 0) {
        Tnpv = PVI->getTrueNumInteractions ();
        continue;
      }
    }

    MyWeight = LumiWeights_.weight(Tnpv);

  }
 

  edm::Handle<GenEventInfoProduct> genEventInfoHandle;

  if (iEvent.getByLabel ("generator", genEventInfoHandle)) {

    double mcWeight = genEventInfoHandle->weight();

    MyWeight *= mcWeight;

  }



  // ------------------------------------------------------------------------------------------------
  //  Z reconstruction


  double dilep_mass = -999.;
  double dilep_pt   = -999.;
  math::XYZTLorentzVector dilep_mom;

  // --- electron selection
  vector < pat::Electron > vect_elem;
  vector < pat::Electron > vect_elep;
  for (pat::ElectronCollection::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele) {

    if ( ele->pt()<20. || fabs(ele->eta())>2.4 ) continue; 

    ele->charge()<0. ? vect_elem.push_back(*ele) : vect_elep.push_back(*ele);
    
  }

  // --- muon selection
  vector < pat::Muon > vect_muom;
  vector < pat::Muon > vect_muop;
  for (pat::MuonCollection::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon) {

    if ( muon->pt()<20. || fabs(muon->eta())>2.4 ) continue;

    muon->charge()<0. ? vect_muom.push_back(*muon) : vect_muop.push_back(*muon);

  }


  // --- dielectron channel

  if ( pileupDT_ == "ee" ){

    if ( vect_elem.size()>0. && vect_elep.size()>0. ){

      dilep_mom  = vect_elem[0].p4() + vect_elep[0].p4();
      dilep_mass = dilep_mom.mass();
      dilep_pt   = dilep_mom.pt();
      if ( dilep_mass>71. && dilep_mass<111. ) ee_event = true;

      if (isMC && ee_event) {
	double scalFac_elem = ElSF_ ->Val(vect_elem[0].pt(), vect_elem[0].eta()) * 
	                      ElSF2_->Val(vect_elem[0].pt(), vect_elem[0].eta());
	double scalFac_elep = ElSF_ ->Val(vect_elep[0].pt(), vect_elep[0].eta()) * 
                              ElSF2_->Val(vect_elep[0].pt(), vect_elep[0].eta());
	MyWeight *= scalFac_elem * scalFac_elep;
      }
      
    }
    
  }
  

  // --- dimuon channel

  else if ( pileupDT_ == "mm" ){

    if ( vect_muom.size()>0. && vect_muop.size()>0. ){

      dilep_mom  = vect_muom[0].p4() + vect_muop[0].p4();
      dilep_mass = dilep_mom.mass();
      dilep_pt   = dilep_mom.pt();
      if ( dilep_mass>71. && dilep_mass<111. ) mm_event = true;

      if (isMC && mm_event) {
	double scalFac_muom = MuSF_->Val(vect_muom[0].pt(), vect_muom[0].eta()) * 
	                      sqrt(MuSF2_->Val(fabs(vect_muom[0].eta()), fabs(vect_muop[0].eta())));
	double scalFac_muop = MuSF_->Val(vect_muop[0].pt(), vect_muop[0].eta()) * 
	                      sqrt(MuSF2_->Val(fabs(vect_muom[0].eta()), fabs(vect_muop[0].eta())));
	MyWeight *= scalFac_muom * scalFac_muop;
      }
    }

  }


  // --- electron-muon channel

  else if ( pileupDT_ == "em" ){

    if ( vect_muom.size()>0. && vect_elep.size()>0. ){

      dilep_mom  = vect_muom[0].p4() + vect_elep[0].p4();
      dilep_mass = dilep_mom.mass();
      dilep_pt   = dilep_mom.pt();
      if ( dilep_mass>71. && dilep_mass<111. ) em_event = true;

      if (isMC && em_event) {
	double scalFac_muom = MuSF_->Val(vect_muom[0].pt(), vect_muom[0].eta()) * 
	                      sqrt(MuSF2_->Val(fabs(vect_muom[0].eta()), fabs(vect_muom[0].eta())));
	double scalFac_elep = ElSF_ ->Val(vect_elep[0].pt(), vect_elep[0].eta()) * 
                              ElSF2_->Val(vect_elep[0].pt(), vect_elep[0].eta());
	MyWeight *= scalFac_muom * scalFac_elep;
      }

    }
    else if ( vect_muop.size()>0. && vect_elem.size()>0. ){

      dilep_mom  = vect_muop[0].p4() + vect_elem[0].p4();
      dilep_mass = dilep_mom.mass();
      dilep_pt   = dilep_mom.pt();
      if ( dilep_mass>71. && dilep_mass<111. ) em_event = true;

      if (isMC) {
	double scalFac_muop = MuSF_->Val(vect_muop[0].pt(), vect_muop[0].eta()) * 
	                      sqrt(MuSF2_->Val(fabs(vect_muop[0].eta()), fabs(vect_muop[0].eta())));
	double scalFac_elem = ElSF_ ->Val(vect_elem[0].pt(), vect_elem[0].eta()) * 
	                      ElSF2_->Val(vect_elem[0].pt(), vect_elem[0].eta());
	MyWeight *= scalFac_muop * scalFac_elem;
      }

    }

  }
  else{
    throw cms::Exception("UnknownChannel") << "Unknown analysis channel: " << pileupDT_ << endl;
  }
    


  // ================================================================================================
  //   Keep only Z-->ee, Z-->mm and em events:
  // ================================================================================================
  
  ee_event = ee_event && (lepton_ == "electron");
  mm_event = mm_event && (lepton_ == "muon");
  em_event = em_event && (lepton_ == "electron+muon");
  if ( !ee_event && !mm_event && !em_event ) return;



  // --- Identify Z->tautau and Z+b/c events: 

  if (isMC) {
    for (std::vector <reco::GenParticle>::const_iterator thepart  = genPart->begin(); 
	                                                 thepart != genPart->end(); 
	                                                 ++thepart) {
      if ((int) abs(thepart->pdgId()) == 23) {
        for (UInt_t i=0; i<thepart->numberOfDaughters(); i++){
	  if (abs(thepart->daughter(i)->pdgId()) == 15 && thepart->daughter(i)->status()==3){
	    is_tautau = true;
	  }
        }
      }
      
      if ( thepart->pt() > 25. &&  fabs(thepart->eta())<2.4 ) {
	
 	if ( fabs(thepart->pdgId()) == 5 ) n_b++; 
 	if ( fabs(thepart->pdgId()) == 4 ) n_c++; 

      }

    }
  }



  // ------------------------------------------------------------------------------------------------
  //  Primary vertex selection

  bool vtx_cut = false;

  edm::Handle < vector < reco::Vertex > > vertices;
  iEvent.getByLabel (edm::InputTag ("goodOfflinePrimaryVertices"), vertices);

  if (vertices->size() > 0) {
    const reco::Vertex* theVertex = &(vertices->front());
    if ( theVertex->ndof()                 >  5  &&
         fabs(theVertex->z())              < 24. &&
         fabs(theVertex->position().rho()) <  2. ) vtx_cut = true;
  } 



  // ================================================================================================
  //  Keep only events with a good primary vertex:
  // ================================================================================================

  if ( !vtx_cut ) return;



  // ------------------------------------------------------------------------------------------------
  //  Jet selection

  vector < pat::Jet > vect_jets;
  vector < pat::Jet > vect_cjets;

  vector < bool > has_b;
  vector < bool > has_c;


  for (vector<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet) {

    
    // --- JEC uncertainty
    double jecUnc = 0.;
    if ( !isMC ) {
      jetCorrectionUncertainty_->setJetPt(jet->pt());
      jetCorrectionUncertainty_->setJetEta(jet->eta());
      jecUnc = jetCorrectionUncertainty_->getUncertainty(true);
    }


    // --- JER corrections
    double jerCorr = 1.;
    if (isMC && jet->genJet()) jerCorr = jetResolutionCorrection( jet->eta(), 
								  jet->pt(), 
								  jet->genJet()->pt(), 
								  par2_ );
    
    pat::Jet jetNew = (*jet);
    math::XYZTLorentzVector jetNew_p4 = jetNew.p4();

    jetNew_p4 = jetNew_p4 * (1.0 + jecUnc * par_) * jerCorr;

    jetNew.setP4(jetNew_p4);


    // --- jet selection:
    if ( jetNew.pt()<30. || fabs(jetNew.eta())>2.4  ) continue;
    vect_jets.push_back(jetNew);


    // --- c tagging:
    double discrCSV = jetNew.bDiscriminator("combinedSecondaryVertexBJetTags");


    if ( ( discrCSV>0.244 && discrCSV<0.898 && fabs(jetNew.eta())<2.4 ) ){

      reco::SecondaryVertexTagInfo const * sv = jetNew.tagInfoSecondaryVertex("secondaryVertex");
      if ( sv && sv->nVertices()>0 )
	vect_cjets.push_back(jetNew);

    }


    // --- Fill histograms to calculate the CSV tagging efficiencies
    if ( isMC ) {

      if ( fabs(jetNew.partonFlavour())==5 )
	h_eff_b->Fill(jetNew.pt(),jetNew.eta());
      else if ( fabs(jetNew.partonFlavour())==4 )
	h_eff_c->Fill(jetNew.pt(),jetNew.eta());
      else
	h_eff_l->Fill(jetNew.pt(),jetNew.eta());
      
      if ( discrCSV>0.244 ){
	if ( fabs(jetNew.partonFlavour())==5 )
	  hc_CSVL_eff_b->Fill(jetNew.pt(),jetNew.eta());
	else if ( fabs(jetNew.partonFlavour())==4 )
	  hc_CSVL_eff_c->Fill(jetNew.pt(),jetNew.eta());
	else
	  hc_CSVL_eff_l->Fill(jetNew.pt(),jetNew.eta());
      }
      
      if ( discrCSV>0.898 ){
	if ( fabs(jetNew.partonFlavour())==5 )
	  hc_CSVT_eff_b->Fill(jetNew.pt(),jetNew.eta());
	else if ( fabs(jetNew.partonFlavour())==4 )
	  hc_CSVT_eff_c->Fill(jetNew.pt(),jetNew.eta());
	else
	  hc_CSVT_eff_l->Fill(jetNew.pt(),jetNew.eta());
      }

    }


  } // end loop over jets

  

  // ================================================================================================
  //  Keep only events with at least one reconstructed jet:
  // ================================================================================================

  if ( vect_jets.size()==0 ) return;


  // --- Z Pt reweighting:
  //if ( ee_event ){
  //  
  //  double w_Zpt = 0.940801 +
  //    0.0125043   * dilep_pt -
  //    0.000900183 * pow(dilep_pt,2) +
  //    2.57315e-05 * pow(dilep_pt,3) -
  //    3.25401e-07 * pow(dilep_pt,4) +
  //    1.77911e-09 * pow(dilep_pt,5) -
  //    3.15907e-12 * pow(dilep_pt,6);
  // 
  //  if ( dilep_pt > 80. ) w_Zpt = 1.;
  //
  //  MyWeight *= w_Zpt;
  //
  //  //cout <<  dilep_pt << " " << w_Zpt << endl;
  //
  //
  //}
  //
  //if ( mm_event ){
  //  
  //  double w_Zpt = 0.963693 + 
  //    0.0115836   * dimuo_pt -
  //    0.000900126 * pow(dimuo_pt,2) +
  //    2.8187e-05  * pow(dimuo_pt,3) -
  //    4.10665e-07 * pow(dimuo_pt,4) +
  //    2.74203e-09 * pow(dimuo_pt,5) -
  //    6.6021e-12  * pow(dimuo_pt,6);
  //  
  //  if ( dimuo_pt > 80. ) w_Zpt = 1.;
  //
  //  MyWeight *= w_Zpt;
  //
  //  //cout << dimuo_pt << " " <<  w_Zpt << endl;
  //
  //}




  double scalFac_b = ctagWeight(isMC, vect_jets);


  h_MET->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
  h_sMET->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
    
  if (isMC) {
     
    if (is_tautau) {
      h_MET_tt->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
      h_sMET_tt->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
    }
    else {
      if ( n_b>0 ){
	h_MET_b->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
	h_sMET_b->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
      }
      else if ( n_b==0 && n_c>0 ){
	h_MET_c->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
	h_sMET_c->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
      }
      else if ( n_b==0 && n_c==0 ){
	h_MET_l->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
	h_sMET_l->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
      }
      
    } 


    if ( vect_cjets.size()>0 ) {

      hc_MET->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
      hc_sMET->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);

      if (isMC) {
	
	if (is_tautau) {
	  hc_MET_tt->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	  hc_sMET_tt->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	}
	else {
	  if ( n_b>0 ){
	    hc_MET_b->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	    hc_sMET_b->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	  }
	  else if ( n_b==0 && n_c>0 ){
	    hc_MET_c->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	    hc_sMET_c->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	  }
	  else if ( n_b==0 && n_c==0 ){
	    hc_MET_l->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	    hc_sMET_l->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	  }

	} 
      }

    }
    
  }

  // ================================================================================================
  //  Keep only events with low MET significance
  // ================================================================================================

  if ( !met->empty() && (*met)[0].significance() > 30. ) return; 


  // ------------------------------------------------------------------------------------------------
  //  Z+jet histograms


  if (!is_tautau){

    n_events[0]++;
    w_events[0] += MyWeight;
    w2_events[0] += MyWeight*MyWeight;

    h_M->Fill(dilep_mass, MyWeight);
    h_Pt->Fill(dilep_pt, MyWeight);
    h_njet->Fill(vect_jets.size(), MyWeight);

    h_CSV->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
    h_BJP->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
    for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
      h_CSV_all->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);

  }

  if (isMC) {
     
    if (is_tautau) {

      n_ztautau[0]++;
      w_ztautau[0] += MyWeight;
      w2_ztautau[0] += MyWeight*MyWeight;

      h_M_tt->Fill(dilep_mass, MyWeight);
      h_Pt_tt->Fill(dilep_pt, MyWeight);
      h_njet_tt->Fill(vect_jets.size(), MyWeight);

      h_CSV_tt->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      h_BJP_tt->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
      for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	h_CSV_all_tt->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	
    }
    else {

      if ( n_b>0 ){
      
	h_M_b->Fill(dilep_mass, MyWeight);
	h_Pt_b->Fill(dilep_pt, MyWeight);
	h_njet_b->Fill(vect_jets.size(), MyWeight);

	h_CSV_b->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	h_BJP_b->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
	for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	  h_CSV_all_b->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      
      }
      else if ( n_b==0 && n_c>0 ){
      
	h_M_c->Fill(dilep_mass, MyWeight);
	h_Pt_c->Fill(dilep_pt, MyWeight);
	h_njet_c->Fill(vect_jets.size(), MyWeight);
	
	h_CSV_c->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	h_BJP_c->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
	for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	  h_CSV_all_c->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      
      }
      else if ( n_b==0 && n_c==0 ){
	  
	h_M_l->Fill(dilep_mass, MyWeight);
	h_Pt_l->Fill(dilep_pt, MyWeight);
	h_njet_l->Fill(vect_jets.size(), MyWeight);

	h_CSV_l->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	h_BJP_l->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
	for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	  h_CSV_all_l->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      
      }

    }
    
  }





  // ================================================================================================
  //  Keep only events with at least one reconstructed c-tagged jet:
  // ================================================================================================

  if ( vect_cjets.size()==0 ) return;


  // Get the b-tagging discriminants for the leading c-tagged jet
  double CSV = vect_cjets[0].bDiscriminator("combinedSecondaryVertexBJetTags");
  double BJP = vect_cjets[0].bDiscriminator("jetBProbabilityBJetTags");
  double JPB = vect_cjets[0].bDiscriminator("jetProbabilityBJetTags");
  double CHP = vect_cjets[0].bDiscriminator("trackCountingHighPurBJetTags");
  double CHE = vect_cjets[0].bDiscriminator("trackCountingHighEffBJetTags");


  double SecVtx_mass      = 0.;
  double SecVtx_mass_corr = 0.;
  double SecVtx_efrac     = 0.;

  reco::SecondaryVertexTagInfo const * svTagInfos = vect_cjets[0].tagInfoSecondaryVertex("secondaryVertex");
  reco::TrackIPTagInfo const * ipTagInfos = vect_cjets[0].tagInfoTrackIP("impactParameter");


  if (svTagInfos && svTagInfos->nVertices() > 0) {

    reco::TrackKinematics vertexKinematics;
    reco::TrackKinematics allKinematics;
 
  
  // --- Total momentum of the secondary vertex tracks:

    for (unsigned int isvx=0; isvx<svTagInfos->nVertices(); ++isvx) {

      const reco::Vertex &vertex = svTagInfos->secondaryVertex(isvx);
      bool hasRefittedTracks = vertex.hasRefittedTracks();

      for (reco::Vertex::trackRef_iterator track  = vertex.tracks_begin(); 
	                                   track != vertex.tracks_end(); 
	                                   ++track) {

	double w = vertex.trackWeight(*track);

	if (w < 0.5)
	  continue;
	if (hasRefittedTracks) {
	  reco::Track actualTrack = vertex.refittedTrack(*track);
	  vertexKinematics.add(actualTrack, w);
	} else {
	  vertexKinematics.add(**track, w);
	}
      }
      
    }
    math::XYZTLorentzVector vtxTrkMom = vertexKinematics.weightedVectorSum();
    SecVtx_mass = vtxTrkMom.M();

    // Apply the SVX mass correction:
    GlobalVector dir = svTagInfos->flightDirection(0);
    double vertexPt2 = math::XYZVector(dir.x(), dir.y(), dir.z()).Cross(vtxTrkMom).Mag2() / dir.mag2();
    SecVtx_mass_corr = std::sqrt(SecVtx_mass*SecVtx_mass + vertexPt2) + std::sqrt(vertexPt2);


    // --- Total momentum of all jet tracks:

    const edm::RefVector<reco::TrackCollection> &tracks = ipTagInfos->selectedTracks();
     
    for(edm::RefVector<reco::TrackCollection>::const_iterator track  = tracks.begin(); 
	                                                      track != tracks.end(); 
	                                                      ++track) {
      allKinematics.add(**track);
    
    }
    math::XYZTLorentzVector jetTrkMom = allKinematics.weightedVectorSum();

    // Calculate the SVX/jet energy ratio:
    if (allKinematics.numberOfTracks())
      SecVtx_efrac = vtxTrkMom.E() / jetTrkMom.E();
    else
      SecVtx_efrac = 1.;
  
  }



  // ------------------------------------------------------------------------------------------------
  //  Z+c histograms


  if (!is_tautau){

    n_events[1]++;
    w_events[1]  += MyWeight*scalFac_b;
    w2_events[1] += MyWeight*MyWeight*scalFac_b*scalFac_b;

    hc_M->Fill(dilep_mass, MyWeight*scalFac_b);
    hc_Pt->Fill(dilep_pt, MyWeight*scalFac_b);
    hc_njet->Fill(vect_jets.size(), MyWeight*scalFac_b);
    hc_ncjet->Fill(vect_cjets.size(), MyWeight*scalFac_b);

    hc_Pt_cjet1->Fill(vect_cjets[0].pt(), MyWeight*scalFac_b);
    
    hc_CSV->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
    hc_BJP->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
    hc_CSV_ctag->Fill(CSV, MyWeight*scalFac_b);
    hc_BJP_ctag->Fill(BJP, MyWeight*scalFac_b);
    for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
      hc_CSV_all->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);

    hc_JPB->Fill(JPB, MyWeight*scalFac_b);
    hc_CHP->Fill(CHP, MyWeight*scalFac_b);
    hc_CHE->Fill(CHE, MyWeight*scalFac_b);
    hc_svxM->Fill(SecVtx_mass,  MyWeight*scalFac_b);
    hc_svxM_corr->Fill(SecVtx_mass_corr,  MyWeight*scalFac_b);
    hc_svxEfr->Fill(SecVtx_efrac,  MyWeight*scalFac_b);

  }
      

  if (isMC) {
     
    if (is_tautau) {

      n_ztautau[1]++;
      w_ztautau[1]  += MyWeight*scalFac_b;
      w2_ztautau[1] += MyWeight*MyWeight*scalFac_b;

      hc_M_tt->Fill(dilep_mass, MyWeight*scalFac_b);
      hc_Pt_tt->Fill(dilep_pt, MyWeight*scalFac_b);
      hc_njet_tt->Fill(vect_jets.size(), MyWeight*scalFac_b);
      hc_ncjet_tt->Fill(vect_cjets.size(), MyWeight*scalFac_b);

      hc_Pt_cjet1_tt->Fill(vect_cjets[0].pt(), MyWeight*scalFac_b);

      hc_CSV_tt->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      hc_BJP_tt->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
      hc_CSV_ctag_tt->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      hc_BJP_ctag_tt->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
      for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	hc_CSV_all_tt->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	
      hc_JPB_tt->Fill(JPB, MyWeight*scalFac_b);
      hc_CHP_tt->Fill(CHP, MyWeight*scalFac_b);
      hc_CHE_tt->Fill(CHE, MyWeight*scalFac_b);
      hc_svxM_tt->Fill(SecVtx_mass,  MyWeight*scalFac_b);
      hc_svxM_corr_tt->Fill(SecVtx_mass_corr,  MyWeight*scalFac_b);
      hc_svxEfr_tt->Fill(SecVtx_efrac,  MyWeight*scalFac_b);

    }
    else {

      //if ( n_b>0 && n_c==0 ){
      if ( n_b>0 ){
      
	hc_M_b->Fill(dilep_mass, MyWeight*scalFac_b);
	hc_Pt_b->Fill(dilep_pt, MyWeight*scalFac_b);
	hc_njet_b->Fill(vect_jets.size(), MyWeight*scalFac_b);
	hc_ncjet_b->Fill(vect_cjets.size(), MyWeight*scalFac_b);

	hc_Pt_cjet1_b->Fill(vect_cjets[0].pt(), MyWeight*scalFac_b);

	hc_CSV_b->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	hc_BJP_b->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
	hc_CSV_ctag_b->Fill(CSV, MyWeight*scalFac_b);
	hc_BJP_ctag_b->Fill(BJP, MyWeight*scalFac_b);
	for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	  hc_CSV_all_b->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);

	hc_JPB_b->Fill(JPB, MyWeight*scalFac_b);
	hc_CHP_b->Fill(CHP, MyWeight*scalFac_b);
	hc_CHE_b->Fill(CHE, MyWeight*scalFac_b);
	hc_svxM_b->Fill(SecVtx_mass,  MyWeight*scalFac_b); 
	hc_svxM_corr_b->Fill(SecVtx_mass_corr,  MyWeight*scalFac_b);
	hc_svxEfr_b->Fill(SecVtx_efrac,  MyWeight*scalFac_b);
     
      }
      else if ( n_b==0 && n_c>0 ){
      
	hc_M_c->Fill(dilep_mass, MyWeight*scalFac_b);
	hc_Pt_c->Fill(dilep_pt, MyWeight*scalFac_b);
	hc_njet_c->Fill(vect_jets.size(), MyWeight*scalFac_b);
	hc_ncjet_c->Fill(vect_cjets.size(), MyWeight*scalFac_b);

	hc_Pt_cjet1_c->Fill(vect_cjets[0].pt(), MyWeight*scalFac_b);

	hc_CSV_c->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	hc_BJP_c->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
	hc_CSV_ctag_c->Fill(CSV, MyWeight*scalFac_b);
	hc_BJP_ctag_c->Fill(BJP, MyWeight*scalFac_b);
	for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	  hc_CSV_all_c->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      
	hc_JPB_c->Fill(JPB, MyWeight*scalFac_b);
	hc_CHP_c->Fill(CHP, MyWeight*scalFac_b);
	hc_CHE_c->Fill(CHE, MyWeight*scalFac_b);
	hc_svxM_c->Fill(SecVtx_mass,  MyWeight*scalFac_b);
	hc_svxM_corr_c->Fill(SecVtx_mass_corr,  MyWeight*scalFac_b);
	hc_svxEfr_c->Fill(SecVtx_efrac,  MyWeight*scalFac_b);

      }
      else if ( n_b==0 && n_c==0 ){
      
	hc_M_l->Fill(dilep_mass, MyWeight*scalFac_b);
	hc_Pt_l->Fill(dilep_pt, MyWeight*scalFac_b);
	hc_njet_l->Fill(vect_jets.size(), MyWeight*scalFac_b);
	hc_ncjet_l->Fill(vect_cjets.size(), MyWeight*scalFac_b);

	hc_Pt_cjet1_l->Fill(vect_cjets[0].pt(), MyWeight*scalFac_b);

	hc_CSV_l->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	hc_BJP_l->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
	hc_CSV_ctag_l->Fill(CSV, MyWeight*scalFac_b);
	hc_BJP_ctag_l->Fill(BJP, MyWeight*scalFac_b);
	for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	  hc_CSV_all_l->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      
	hc_JPB_l->Fill(JPB, MyWeight*scalFac_b);
	hc_CHP_l->Fill(CHP, MyWeight*scalFac_b);
	hc_CHE_l->Fill(CHE, MyWeight*scalFac_b);
	hc_svxM_l->Fill(SecVtx_mass,  MyWeight*scalFac_b);
	hc_svxM_corr_l->Fill(SecVtx_mass_corr,  MyWeight*scalFac_b);
	hc_svxEfr_l->Fill(SecVtx_efrac,  MyWeight*scalFac_b);

      }

    }
      
  }



}



// ------------ method called once each job just before starting event loop ------------
void ZcAnalyzer::beginJob () {

  jetCorrectionUncertainty_ = new JetCorrectionUncertainty(path_ + "/" + "Summer13_V4_DATA_Uncertainty_AK5PFchs.txt");
  LumiWeights_ = edm::LumiReWeighting(path_ + "/" + "pileup_" + pileupMC_ + ".root", path_ + "/" + "pileup_2012_" + pileupDT_ + ".root", "pileup", "pileup");

  ElSF_  = new table(path_ + "/" + "ele_eff.txt"  );
  ElSF2_ = new table(path_ + "/" + "ele_eff2.txt" );
  MuSF_  = new table(path_ + "/" + "muon_eff.txt" );
  MuSF2_ = new table(path_ + "/" + "muon_eff2.txt");

  BtSFT_  = new table(path_ + "/" + "btag_CSVT_b_SF.txt" ); // btagging scale factors for CSVT (SFb = SFc) 
  BtSFL_  = new table(path_ + "/" + "btag_CSVL_b_SF.txt" ); // btagging scale factors for CSVL (SFb = SFc) 
  LtSFT_  = new table(path_ + "/" + "btag_CSVT_dusg_SF.txt"); // light flavour scale factors for CSVT
  LtSFL_  = new table(path_ + "/" + "btag_CSVL_dusg_SF.txt"); // light flavour scale factors for CSVL
  BtEffT_ = new table(path_ + "/" + "btag_CSVT_b_eff.txt" ); // b jets efficiencies for CSVT
  BtEffL_ = new table(path_ + "/" + "btag_CSVL_b_eff.txt" ); // b jets efficiencies for CSVL
  CtEffT_ = new table(path_ + "/" + "btag_CSVT_c_eff.txt" ); // c jets efficiencies for CSVT
  CtEffL_ = new table(path_ + "/" + "btag_CSVL_c_eff.txt" ); // c jets efficiencies for CSVL
  LtEffT_ = new table(path_ + "/" + "btag_CSVT_dusg_eff.txt"); // light flavour efficiencies for CSVT
  LtEffL_ = new table(path_ + "/" + "btag_CSVL_dusg_eff.txt"); // light flavour efficiencies for CSVL


  // --- Initialize the counters: 
  for (int i=0; i<2; ++i){
    n_ztautau[i]  = 0 ;
    w_ztautau[i]  = 0.;
    w2_ztautau[i] = 0.;
    n_events[i]   = 0 ;
    w_events[i]   = 0.;
    w2_events[i]  = 0.;
  }

}

// ------------ method called once each job just after ending the event loop ------------
void ZcAnalyzer::endJob () {
  delete jetCorrectionUncertainty_;

  delete ElSF_;
  delete ElSF2_;
  delete MuSF_;
  delete MuSF2_;
  delete BtSFT_;
  delete BtSFL_;
  delete LtSFT_;
  delete LtSFL_;
  delete BtEffT_;
  delete BtEffL_;
  delete CtEffT_;
  delete CtEffL_;
  delete LtEffT_;
  delete LtEffL_;


  // --- Printout the counters:
  std::cout << "======================================================================" << std::endl;
  std::cout << " " << pileupDT_ << " channel:" << std::endl;
  std::cout <<  std::endl;
  std::cout << "Z-->tautau yield (inclusive) = " << n_ztautau[0] << std::endl;
  std::cout << "Z-->tautau yield (tagged)    = " << n_ztautau[1] << std::endl; 
  std::cout << "Z-->tautau weighted yield (inclusive) = "
	    << w_ztautau[0] << " " << w2_ztautau[0] << std::endl;
  std::cout << "Z-->tautau weighted yield (tagged)    = " 
	    << w_ztautau[1] << " " << w2_ztautau[1] << std::endl;

  std::cout <<  std::endl;

  std::cout << "Yield (inclusive) = " << n_events[0] << std::endl;
  std::cout << "Yield (tagged)    = " << n_events[1] << std::endl;
  std::cout << "Weighted yield (inclusive) = " 
	    << w_events[0] << " " << w2_events[0] << std::endl;
  std::cout << "Weighted yield (tagged)    = "  
	    << w_events[1] << " " << w2_events[1] << std::endl;

}

// ------------ method called when starting to processes a run ------------
void ZcAnalyzer::beginRun (edm::Run const &, edm::EventSetup const &) {
}

// ------------ method called when ending the processing of a run ------------
void ZcAnalyzer::endRun (edm::Run const &, edm::EventSetup const &) {
}

// ------------ method called when starting to processes a luminosity block ------------
void ZcAnalyzer::beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &) {
}

// ------------ method called when ending the processing of a luminosity block ------------
void ZcAnalyzer::endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &) {
}

// define this as a plug-in
DEFINE_FWK_MODULE (ZcAnalyzer);
