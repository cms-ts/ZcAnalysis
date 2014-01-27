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
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stddef.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/View.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "PhysicsTools/CandUtils/interface/CandCombiner.h"
#include "CommonTools/Utils/interface/MassRangeSelector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"


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

    if (syst == 0) jetPtNew = jetPtGen + correctionFactor[index] * (jetPt-jetPtGen);
    if (syst == +1) jetPtNew = jetPtGen + correctionFactorUp[index] * (jetPt-jetPtGen);
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

  int    n_ztautau[2][2];
  double w_ztautau[2][2];
  double w2_ztautau[2][2];

  int    n_events[2][2];
  double w_events[2][2];
  double w2_events[2][2];


  // ================================================================================================
  //  Histogram declaration


  // --- Z+jets

  TH1F*     h_M_ee;
  TH1F*     h_M_ee_b;
  TH1F*     h_M_ee_c;
  TH1F*     h_M_ee_l;
  TH1F*     h_M_ee_tt;
  TH1F*     h_M_mm;
  TH1F*     h_M_mm_b;
  TH1F*     h_M_mm_c;
  TH1F*     h_M_mm_l;
  TH1F*     h_M_mm_tt;

  TH1F*     h_njet_ee;
  TH1F*     h_njet_ee_b;
  TH1F*     h_njet_ee_c;
  TH1F*     h_njet_ee_l;
  TH1F*     h_njet_ee_tt;
  TH1F*     h_njet_mm;
  TH1F*     h_njet_mm_b;
  TH1F*     h_njet_mm_c;
  TH1F*     h_njet_mm_l;
  TH1F*     h_njet_mm_tt;

  TH1F*     h_MET_ee;
  TH1F*     h_MET_ee_b;
  TH1F*     h_MET_ee_c;
  TH1F*     h_MET_ee_l;
  TH1F*     h_MET_ee_tt;
  TH1F*     h_MET_mm;
  TH1F*     h_MET_mm_b;
  TH1F*     h_MET_mm_c;
  TH1F*     h_MET_mm_l;
  TH1F*     h_MET_mm_tt;

  TH1F*     h_sMET_ee;
  TH1F*     h_sMET_ee_b;
  TH1F*     h_sMET_ee_c;
  TH1F*     h_sMET_ee_l;
  TH1F*     h_sMET_ee_tt;
  TH1F*     h_sMET_mm;
  TH1F*     h_sMET_mm_b;
  TH1F*     h_sMET_mm_c;
  TH1F*     h_sMET_mm_l;
  TH1F*     h_sMET_mm_tt;

  TH1F*     h_CSV_ee;
  TH1F*     h_CSV_ee_b;
  TH1F*     h_CSV_ee_c;
  TH1F*     h_CSV_ee_b_split;
  TH1F*     h_CSV_ee_c_split;
  TH1F*     h_CSV_ee_l;
  TH1F*     h_CSV_ee_tt;
  TH1F*     h_CSV_mm;
  TH1F*     h_CSV_mm_b;
  TH1F*     h_CSV_mm_c;
  TH1F*     h_CSV_mm_b_split;
  TH1F*     h_CSV_mm_c_split;
  TH1F*     h_CSV_mm_l;
  TH1F*     h_CSV_mm_tt;

  TH1F*     h_BJP_ee;
  TH1F*     h_BJP_ee_b;
  TH1F*     h_BJP_ee_c;
  TH1F*     h_BJP_ee_l;
  TH1F*     h_BJP_ee_tt;
  TH1F*     h_BJP_mm;
  TH1F*     h_BJP_mm_b;
  TH1F*     h_BJP_mm_c;
  TH1F*     h_BJP_mm_l;
  TH1F*     h_BJP_mm_tt;

  TH1F*     h_CSV_all_ee;
  TH1F*     h_CSV_all_ee_b;
  TH1F*     h_CSV_all_ee_c;
  TH1F*     h_CSV_all_ee_l;
  TH1F*     h_CSV_all_ee_tt;
  TH1F*     h_CSV_all_mm;
  TH1F*     h_CSV_all_mm_b;
  TH1F*     h_CSV_all_mm_c;
  TH1F*     h_CSV_all_mm_l;
  TH1F*     h_CSV_all_mm_tt;

  TH2F*     h_eff_b;
  TH2F*     h_eff_c;
  TH2F*     h_eff_l;


  // --- Z+c

  TH1F*     hc_M_ee;
  TH1F*     hc_M_ee_b;
  TH1F*     hc_M_ee_c;
  TH1F*     hc_M_ee_l;
  TH1F*     hc_M_ee_tt;
  TH1F*     hc_M_mm;
  TH1F*     hc_M_mm_b;
  TH1F*     hc_M_mm_c;
  TH1F*     hc_M_mm_l;
  TH1F*     hc_M_mm_tt;

  TH1F*     hc_njet_ee;
  TH1F*     hc_njet_ee_b;
  TH1F*     hc_njet_ee_c;
  TH1F*     hc_njet_ee_l;
  TH1F*     hc_njet_ee_tt;
  TH1F*     hc_njet_mm;
  TH1F*     hc_njet_mm_b;
  TH1F*     hc_njet_mm_c;
  TH1F*     hc_njet_mm_l;
  TH1F*     hc_njet_mm_tt;

  TH1F*     hc_MET_ee;
  TH1F*     hc_MET_ee_b;
  TH1F*     hc_MET_ee_c;
  TH1F*     hc_MET_ee_l;
  TH1F*     hc_MET_ee_tt;
  TH1F*     hc_MET_mm;
  TH1F*     hc_MET_mm_b;
  TH1F*     hc_MET_mm_c;
  TH1F*     hc_MET_mm_l;
  TH1F*     hc_MET_mm_tt;

  TH1F*     hc_sMET_ee;
  TH1F*     hc_sMET_ee_b;
  TH1F*     hc_sMET_ee_c;
  TH1F*     hc_sMET_ee_l;
  TH1F*     hc_sMET_ee_tt;
  TH1F*     hc_sMET_mm;
  TH1F*     hc_sMET_mm_b;
  TH1F*     hc_sMET_mm_c;
  TH1F*     hc_sMET_mm_l;
  TH1F*     hc_sMET_mm_tt;

  TH1F*     hc_CSV_ee;
  TH1F*     hc_CSV_ee_b;
  TH1F*     hc_CSV_ee_c;
  TH1F*     hc_CSV_ee_l;
  TH1F*     hc_CSV_ee_b_split;
  TH1F*     hc_CSV_ee_c_split;
  TH1F*     hc_CSV_ee_tt;
  TH1F*     hc_CSV_mm;
  TH1F*     hc_CSV_mm_b;
  TH1F*     hc_CSV_mm_c;
  TH1F*     hc_CSV_mm_b_split;
  TH1F*     hc_CSV_mm_c_split;
  TH1F*     hc_CSV_mm_l;
  TH1F*     hc_CSV_mm_tt;

  TH1F*     hc_BJP_ee;
  TH1F*     hc_BJP_ee_b;
  TH1F*     hc_BJP_ee_c;
  TH1F*     hc_BJP_ee_l;
  TH1F*     hc_BJP_ee_tt;
  TH1F*     hc_BJP_mm;
  TH1F*     hc_BJP_mm_b;
  TH1F*     hc_BJP_mm_c;
  TH1F*     hc_BJP_mm_l;
  TH1F*     hc_BJP_mm_tt;

  TH1F*     hc_CSV_all_ee;
  TH1F*     hc_CSV_all_ee_b;
  TH1F*     hc_CSV_all_ee_c;
  TH1F*     hc_CSV_all_ee_l;
  TH1F*     hc_CSV_all_ee_tt;
  TH1F*     hc_CSV_all_mm;
  TH1F*     hc_CSV_all_mm_b;
  TH1F*     hc_CSV_all_mm_c;
  TH1F*     hc_CSV_all_mm_l;
  TH1F*     hc_CSV_all_mm_tt;

  TH2F*     hc_CSVL_eff_b;
  TH2F*     hc_CSVL_eff_c;
  TH2F*     hc_CSVL_eff_l;
  TH2F*     hc_CSVT_eff_b;
  TH2F*     hc_CSVT_eff_c;
  TH2F*     hc_CSVT_eff_l;

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
  path_ = iConfig.getUntrackedParameter < std::string > ("path", "/gpfs/cms/users/casarsa/analysis/CMSSW_5_3_11_patch6/src/ZcAnalysis/ZcSkim/SF");
  icut_ = iConfig.getUntrackedParameter <unsigned int> ("icut", 0);
  usePartonFlavour_ = iConfig.getUntrackedParameter <bool> ("usePartonFlavour", false);
  pcut_ = iConfig.getUntrackedParameter <bool> ("pcut", false);
  useDeltaR_ = iConfig.getUntrackedParameter <bool> ("useDeltaR", false);

  // now do what ever initialization is needed
  edm::Service < TFileService > fs;

  h_M_ee    = fs->make < TH1F > ("h_M_ee",    "ee mass;M_{ee} [GeV/c^{2}]", 50, 71., 111.);
  h_M_ee->Sumw2();
  h_M_ee_b  = fs->make < TH1F > ("h_M_ee_b",  "ee mass (b);M_{ee} [GeV/c^{2}]", 50, 71., 111.);
  h_M_ee_b->Sumw2();
  h_M_ee_c  = fs->make < TH1F > ("h_M_ee_c",  "ee mass (c);M_{ee} [GeV/c^{2}]", 50, 71., 111.);
  h_M_ee_c->Sumw2();
  h_M_ee_l  = fs->make < TH1F > ("h_M_ee_l",  "ee mass (dusg);M_{ee} [GeV/c^{2}]", 50, 71., 111.);
  h_M_ee_l->Sumw2();
  h_M_ee_tt = fs->make < TH1F > ("h_M_ee_tt", "ee mass (Z#rightarrow#tau#tau);M_{ee} [GeV/c^{2}]", 50, 71., 111.);
  h_M_ee_tt->Sumw2();

  h_M_mm    = fs->make < TH1F > ("h_M_mm",    "#mu#mu mass;M_{#mu#mu} [GeV/c^{2}]", 50, 71., 111.);
  h_M_mm->Sumw2();
  h_M_mm_b  = fs->make < TH1F > ("h_M_mm_b",  "#mu#mu mass (b);M_{#mu#mu} [GeV/c^{2}]", 50, 71., 111.);
  h_M_mm_b->Sumw2();
  h_M_mm_c  = fs->make < TH1F > ("h_M_mm_c",  "#mu#mu mass (c);M_{#mu#mu} [GeV/c^{2}]", 50, 71., 111.);
  h_M_mm_c->Sumw2();
  h_M_mm_l  = fs->make < TH1F > ("h_M_mm_l",  "#mu#mu mass (dusg);M_{#mu#mu} [GeV/c^{2}]", 50, 71., 111.);
  h_M_mm_l->Sumw2();
  h_M_mm_tt = fs->make < TH1F > ("h_M_mm_tt", "#mu#mu mass (Z#rightarrow#tau#tau);M_{#mu#mu} [GeV/c^{2}]", 50, 71., 111.);
  h_M_mm_tt->Sumw2();

  h_njet_ee    = fs->make < TH1F > ("h_njet_ee",    "ee jet multiplicity;N_{jet}", 10, 0., 10.);
  h_njet_ee_b  = fs->make < TH1F > ("h_njet_ee_b",  "ee jet multiplicity (b);N_{jet}", 10, 0., 10.);
  h_njet_ee_c  = fs->make < TH1F > ("h_njet_ee_c",  "ee jet multiplicity (c);N_{jet}", 10, 0., 10.);
  h_njet_ee_l  = fs->make < TH1F > ("h_njet_ee_l",  "ee jet multiplicity (dusg);N_{jet}", 10, 0., 10.);
  h_njet_ee_tt = fs->make < TH1F > ("h_njet_ee_tt", "ee jet multiplicity (Z#rightarrow#tau#tau);N_{jet}", 10, 0., 10.);

  h_njet_mm    = fs->make < TH1F > ("h_njet_mm",    "#mu#mu jet multiplicity;N_{jet}", 10, 0., 10.);
  h_njet_mm_b  = fs->make < TH1F > ("h_njet_mm_b",  "#mu#mu jet multiplicity (b);N_{jet}", 10, 0., 10.);
  h_njet_mm_c  = fs->make < TH1F > ("h_njet_mm_c",  "#mu#mu jet multiplicity (c);N_{jet}", 10, 0., 10.);
  h_njet_mm_l  = fs->make < TH1F > ("h_njet_mm_l",  "#mu#mu jet multiplicity (dusg);N_{jet}", 10, 0., 10.);
  h_njet_mm_tt = fs->make < TH1F > ("h_njet_mm_tt", "#mu#mu jet multiplicity (Z#rightarrow#tau#tau);N_{jet}", 10, 0., 10.);

  h_MET_ee    = fs->make < TH1F > ("h_MET_ee",    "ee MET;MET [GeV]", 100, 0., 250.);
  h_MET_ee_b  = fs->make < TH1F > ("h_MET_ee_b ", "ee MET (b);MET [GeV]", 100, 0., 250.);
  h_MET_ee_c  = fs->make < TH1F > ("h_MET_ee_c ", "ee MET (c);MET [GeV]", 100, 0., 250.);
  h_MET_ee_l  = fs->make < TH1F > ("h_MET_ee_l ", "ee MET (dusg);MET [GeV]", 100, 0., 250.);
  h_MET_ee_tt = fs->make < TH1F > ("h_MET_ee_tt", "ee MET (Z#rightarrow#tau#tau);MET [GeV]", 100, 0., 250.);

  h_MET_mm    = fs->make < TH1F > ("h_MET_mm",    "#mu#mu MET;MET [GeV]", 100, 0., 250.);
  h_MET_mm_b  = fs->make < TH1F > ("h_MET_mm_b ", "#mu#mu MET (b);MET [GeV]", 100, 0., 250.);
  h_MET_mm_c  = fs->make < TH1F > ("h_MET_mm_c ", "#mu#mu MET (c);MET [GeV]", 100, 0., 250.);
  h_MET_mm_l  = fs->make < TH1F > ("h_MET_mm_l ", "#mu#mu MET (dusg);MET [GeV]", 100, 0., 250.);
  h_MET_mm_tt = fs->make < TH1F > ("h_MET_mm_tt", "#mu#mu MET (Z#rightarrow#tau#tau);MET [GeV]", 100, 0., 250.);

  h_sMET_ee    = fs->make < TH1F > ("h_sMET_ee",    "ee MET significance;MET [GeV]", 100, 0., 250.);
  h_sMET_ee_b  = fs->make < TH1F > ("h_sMET_ee_b ", "ee MET significance (b);MET [GeV]", 100, 0., 250.);
  h_sMET_ee_c  = fs->make < TH1F > ("h_sMET_ee_c ", "ee MET significance (c);MET [GeV]", 100, 0., 250.);
  h_sMET_ee_l  = fs->make < TH1F > ("h_sMET_ee_l ", "ee MET significance (dusg);MET [GeV]", 100, 0., 250.);
  h_sMET_ee_tt = fs->make < TH1F > ("h_sMET_ee_tt", "ee MET significance (Z#rightarrow#tau#tau);MET [GeV]", 100, 0., 250.);

  h_sMET_mm    = fs->make < TH1F > ("h_sMET_mm",    "#mu#mu MET significance;MET [GeV]", 100, 0., 250.);
  h_sMET_mm_b  = fs->make < TH1F > ("h_sMET_mm_b ", "#mu#mu MET significance (b);MET [GeV]", 100, 0., 250.);
  h_sMET_mm_c  = fs->make < TH1F > ("h_sMET_mm_c ", "#mu#mu MET significance (c);MET [GeV]", 100, 0., 250.);
  h_sMET_mm_l  = fs->make < TH1F > ("h_sMET_mm_l ", "#mu#mu MET significance (dusg);MET [GeV]", 100, 0., 250.);
  h_sMET_mm_tt = fs->make < TH1F > ("h_sMET_mm_tt", "#mu#mu MET significance (Z#rightarrow#tau#tau);MET [GeV]", 100, 0., 250.);

  h_CSV_ee    = fs->make < TH1F > ("h_CSV_ee",    "ee CSV discriminant;CSV", 100, 0., 1.);
  h_CSV_ee_b  = fs->make < TH1F > ("h_CSV_ee_b",  "ee CSV discriminant (b);CSV", 100, 0., 1.);
  h_CSV_ee_c  = fs->make < TH1F > ("h_CSV_ee_c",  "ee CSV discriminant (c);CSV", 100, 0., 1.);
  h_CSV_ee_b_split  = fs->make < TH1F > ("h_CSV_ee_b_split",  "ee CSV discriminant (g->bb);CSV", 100, 0., 1.);
  h_CSV_ee_c_split  = fs->make < TH1F > ("h_CSV_ee_c_split",  "ee CSV discriminant (g->cc);CSV", 100, 0., 1.);
  h_CSV_ee_l  = fs->make < TH1F > ("h_CSV_ee_l",  "ee CSV discriminant (dusg);CSV", 100, 0., 1.);
  h_CSV_ee_tt = fs->make < TH1F > ("h_CSV_ee_tt", "ee CSV discriminant (Z#rightarrow#tau#tau);CSV", 100, 0., 1.);

  h_CSV_mm    = fs->make < TH1F > ("h_CSV_mm",    "#mu#mu  CSV discriminant;CSV", 100, 0., 1.);
  h_CSV_mm_b  = fs->make < TH1F > ("h_CSV_mm_b",  "#mu#mu  CSV discriminant (b);CSV", 100, 0., 1.);
  h_CSV_mm_c  = fs->make < TH1F > ("h_CSV_mm_c",  "#mu#mu  CSV discriminant (c);CSV", 100, 0., 1.);
  h_CSV_mm_b_split  = fs->make < TH1F > ("h_CSV_mm_b_split",  "#mu#mu CSV discriminant (g->bb);CSV", 100, 0., 1.);
  h_CSV_mm_c_split  = fs->make < TH1F > ("h_CSV_mm_c_split",  "#mu#mu CSV discriminant (g->cc);CSV", 100, 0., 1.);
  h_CSV_mm_l  = fs->make < TH1F > ("h_CSV_mm_l",  "#mu#mu  CSV discriminant (dusg);CSV", 100, 0., 1.);
  h_CSV_mm_tt = fs->make < TH1F > ("h_CSV_mm_tt", "#mu#mu  CSV discriminant (Z#rightarrow#tau#tau);CSV", 100, 0., 1.);

  h_CSV_all_ee    = fs->make < TH1F > ("h_CSV_all_ee",    "ee CSV discriminant;CSV", 100, 0., 1.);
  h_CSV_all_ee_b  = fs->make < TH1F > ("h_CSV_all_ee_b",  "ee CSV discriminant (b);CSV", 100, 0., 1.);
  h_CSV_all_ee_c  = fs->make < TH1F > ("h_CSV_all_ee_c",  "ee CSV discriminant (c);CSV", 100, 0., 1.);
  h_CSV_all_ee_l  = fs->make < TH1F > ("h_CSV_all_ee_l",  "ee CSV discriminant (dusg);CSV", 100, 0., 1.);
  h_CSV_all_ee_tt = fs->make < TH1F > ("h_CSV_all_ee_tt", "ee CSV discriminant (Z#rightarrow#tau#tau);CSV", 100, 0., 1.);

  h_CSV_all_mm    = fs->make < TH1F > ("h_CSV_all_mm",    "#mu#mu  CSV discriminant;CSV", 100, 0., 1.);
  h_CSV_all_mm_b  = fs->make < TH1F > ("h_CSV_all_mm_b",  "#mu#mu  CSV discriminant (b);CSV", 100, 0., 1.);
  h_CSV_all_mm_c  = fs->make < TH1F > ("h_CSV_all_mm_c",  "#mu#mu  CSV discriminant (c);CSV", 100, 0., 1.);
  h_CSV_all_mm_l  = fs->make < TH1F > ("h_CSV_all_mm_l",  "#mu#mu  CSV discriminant (dusg);CSV", 100, 0., 1.);
  h_CSV_all_mm_tt = fs->make < TH1F > ("h_CSV_all_mm_tt", "#mu#mu  CSV discriminant (Z#rightarrow#tau#tau);CSV", 100, 0., 1.);

  h_BJP_ee    = fs->make < TH1F > ("h_BJP_ee",    "ee BJP discriminant;BJP", 100, 0., 10.);
  h_BJP_ee_b  = fs->make < TH1F > ("h_BJP_ee_b",  "ee BJP discriminant (b);BJP", 100, 0., 10.);
  h_BJP_ee_c  = fs->make < TH1F > ("h_BJP_ee_c",  "ee BJP discriminant (c);BJP", 100, 0., 10.);
  h_BJP_ee_l  = fs->make < TH1F > ("h_BJP_ee_l",  "ee BJP discriminant (dusg);BJP", 100, 0., 10.);
  h_BJP_ee_tt = fs->make < TH1F > ("h_BJP_ee_tt", "ee BJP discriminant (Z#rightarrow#tau#tau);BJP", 100, 0., 10.);

  h_BJP_mm    = fs->make < TH1F > ("h_BJP_mm",    "#mu#mu  BJP discriminant;BJP", 100, 0., 10.);
  h_BJP_mm_b  = fs->make < TH1F > ("h_BJP_mm_b",  "#mu#mu  BJP discriminant (b);BJP", 100, 0., 10.);
  h_BJP_mm_c  = fs->make < TH1F > ("h_BJP_mm_c",  "#mu#mu  BJP discriminant (c);BJP", 100, 0., 10.);
  h_BJP_mm_l  = fs->make < TH1F > ("h_BJP_mm_l",  "#mu#mu  BJP discriminant (dusg);BJP", 100, 0., 10.);
  h_BJP_mm_tt = fs->make < TH1F > ("h_BJP_mm_tt", "#mu#mu  BJP discriminant (Z#rightarrow#tau#tau);BJP", 100, 0., 10.);


  const float xbins[17] = {20., 30., 40., 50., 60., 70., 80., 100., 120., 160., 210., 260., 320., 400., 500., 600., 800.};
  const float ybins[9]  = {-2.4, -1.5, -1.0, -0.5, 0., 0.5, 1.0, 1.5, 2.4};

  h_eff_b = fs->make < TH2F > ("h_eff_b","before tagging (b)", 16, xbins, 8, ybins);
  h_eff_b->Sumw2();
  h_eff_c = fs->make < TH2F > ("h_eff_c","before tagging (c)", 16, xbins, 8, ybins);
  h_eff_c->Sumw2();
  h_eff_l = fs->make < TH2F > ("h_eff_l","before tagging (dusg)", 16, xbins, 8, ybins);
  h_eff_l->Sumw2();

  hc_CSVL_eff_b = fs->make < TH2F > ("hc_CSVL_eff_b","after CSVL tagging (b)", 16, xbins, 8, ybins);
  hc_CSVL_eff_b->Sumw2();
  hc_CSVL_eff_c = fs->make < TH2F > ("hc_CSVL_eff_c","after CSVL tagging (c)", 16, xbins, 8, ybins);
  hc_CSVL_eff_c->Sumw2();
  hc_CSVL_eff_l = fs->make < TH2F > ("hc_CSVL_eff_l","after CSVL tagging (dusg)", 16, xbins, 8, ybins);
  hc_CSVL_eff_l->Sumw2();
  hc_CSVT_eff_b = fs->make < TH2F > ("hc_CSVT_eff_b","after CSVT tagging (b)", 16, xbins, 8, ybins);
  hc_CSVT_eff_b->Sumw2();
  hc_CSVT_eff_c = fs->make < TH2F > ("hc_CSVT_eff_c","after CSVT tagging (c)", 16, xbins, 8, ybins);
  hc_CSVT_eff_c->Sumw2();
  hc_CSVT_eff_l = fs->make < TH2F > ("hc_CSVT_eff_l","after CSVT tagging (dusg)", 16, xbins, 8, ybins);
  hc_CSVT_eff_l->Sumw2();


  hc_M_ee    = fs->make < TH1F > ("hc_M_ee",    "ee mass;M_{ee} [GeV/c^{2}]", 50, 71., 111.);
  hc_M_ee_b  = fs->make < TH1F > ("hc_M_ee_b",  "ee mass (b);M_{ee} [GeV/c^{2}]", 50, 71., 111.);
  hc_M_ee_c  = fs->make < TH1F > ("hc_M_ee_c",  "ee mass (c);M_{ee} [GeV/c^{2}]", 50, 71., 111.);
  hc_M_ee_l  = fs->make < TH1F > ("hc_M_ee_l",  "ee mass (dusg);M_{ee} [GeV/c^{2}]", 50, 71., 111.);
  hc_M_ee_tt = fs->make < TH1F > ("hc_M_ee_tt", "ee mass (Z#rightarrow#tau#tau);M_{ee} [GeV/c^{2}]", 50, 71., 111.);

  hc_M_mm    = fs->make < TH1F > ("hc_M_mm",    "#mu#mu mass;M_{#mu#mu} [GeV/c^{2}]", 50, 71., 111.);
  hc_M_mm_b  = fs->make < TH1F > ("hc_M_mm_b",  "#mu#mu mass (b);M_{#mu#mu} [GeV/c^{2}]", 50, 71., 111.);
  hc_M_mm_c  = fs->make < TH1F > ("hc_M_mm_c",  "#mu#mu mass (c);M_{#mu#mu} [GeV/c^{2}]", 50, 71., 111.);
  hc_M_mm_l  = fs->make < TH1F > ("hc_M_mm_l",  "#mu#mu mass (dusg);M_{#mu#mu} [GeV/c^{2}]", 50, 71., 111.);
  hc_M_mm_tt = fs->make < TH1F > ("hc_M_mm_tt", "#mu#mu mass (Z#rightarrow#tau#tau);M_{#mu#mu} [GeV/c^{2}]", 50, 71., 111.);

  hc_njet_ee    = fs->make < TH1F > ("hc_njet_ee",    "ee jet multiplicity;N_{jet}", 10, 0., 10.);
  hc_njet_ee_b  = fs->make < TH1F > ("hc_njet_ee_b",  "ee jet multiplicity (b);N_{jet}", 10, 0., 10.);
  hc_njet_ee_c  = fs->make < TH1F > ("hc_njet_ee_c",  "ee jet multiplicity (c);N_{jet}", 10, 0., 10.);
  hc_njet_ee_l  = fs->make < TH1F > ("hc_njet_ee_l",  "ee jet multiplicity (dusg);N_{jet}", 10, 0., 10.);
  hc_njet_ee_tt = fs->make < TH1F > ("hc_njet_ee_tt", "ee jet multiplicity (Z#rightarrow#tau#tau);N_{jet}", 10, 0., 10.);

  hc_njet_mm    = fs->make < TH1F > ("hc_njet_mm",    "#mu#mu jet multiplicity;N_{jet}", 10, 0., 10.);
  hc_njet_mm_b  = fs->make < TH1F > ("hc_njet_mm_b",  "#mu#mu jet multiplicity (b);N_{jet}", 10, 0., 10.);
  hc_njet_mm_c  = fs->make < TH1F > ("hc_njet_mm_c",  "#mu#mu jet multiplicity (c);N_{jet}", 10, 0., 10.);
  hc_njet_mm_l  = fs->make < TH1F > ("hc_njet_mm_l",  "#mu#mu jet multiplicity (dusg);N_{jet}", 10, 0., 10.);
  hc_njet_mm_tt = fs->make < TH1F > ("hc_njet_mm_tt", "#mu#mu jet multiplicity (Z#rightarrow#tau#tau);N_{jet}", 10, 0., 10.);

  hc_MET_ee    = fs->make < TH1F > ("hc_MET_ee",    "ee MET;MET [GeV]", 100, 0., 250.);
  hc_MET_ee_b  = fs->make < TH1F > ("hc_MET_ee_b ", "ee MET (b);MET [GeV]", 100, 0., 250.);
  hc_MET_ee_c  = fs->make < TH1F > ("hc_MET_ee_c ", "ee MET (c);MET [GeV]", 100, 0., 250.);
  hc_MET_ee_l  = fs->make < TH1F > ("hc_MET_ee_l ", "ee MET (dusg);MET [GeV]", 100, 0., 250.);
  hc_MET_ee_tt = fs->make < TH1F > ("hc_MET_ee_tt", "ee MET (Z#rightarrow#tau#tau);MET [GeV]", 100, 0., 250.);

  hc_MET_mm    = fs->make < TH1F > ("hc_MET_mm",    "#mu#mu MET;MET [GeV]", 100, 0., 250.);
  hc_MET_mm_b  = fs->make < TH1F > ("hc_MET_mm_b ", "#mu#mu MET (b);MET [GeV]", 100, 0., 250.);
  hc_MET_mm_c  = fs->make < TH1F > ("hc_MET_mm_c ", "#mu#mu MET (c);MET [GeV]", 100, 0., 250.);
  hc_MET_mm_l  = fs->make < TH1F > ("hc_MET_mm_l ", "#mu#mu MET (dusg);MET [GeV]", 100, 0., 250.);
  hc_MET_mm_tt = fs->make < TH1F > ("hc_MET_mm_tt", "#mu#mu MET (Z#rightarrow#tau#tau);MET [GeV]", 100, 0., 250.);

  hc_sMET_ee    = fs->make < TH1F > ("hc_sMET_ee",    "ee MET significance;MET [GeV]", 100, 0., 250.);
  hc_sMET_ee_b  = fs->make < TH1F > ("hc_sMET_ee_b ", "ee MET significance (b);MET [GeV]", 100, 0., 250.);
  hc_sMET_ee_c  = fs->make < TH1F > ("hc_sMET_ee_c ", "ee MET significance (c);MET [GeV]", 100, 0., 250.);
  hc_sMET_ee_l  = fs->make < TH1F > ("hc_sMET_ee_l ", "ee MET significance (dusg);MET [GeV]", 100, 0., 250.);
  hc_sMET_ee_tt = fs->make < TH1F > ("hc_sMET_ee_tt", "ee MET significance (Z#rightarrow#tau#tau);MET [GeV]", 100, 0., 250.);

  hc_sMET_mm    = fs->make < TH1F > ("hc_sMET_mm",    "#mu#mu MET significance;MET [GeV]", 100, 0., 250.);
  hc_sMET_mm_b  = fs->make < TH1F > ("hc_sMET_mm_b ", "#mu#mu MET significance (b);MET [GeV]", 100, 0., 250.);
  hc_sMET_mm_c  = fs->make < TH1F > ("hc_sMET_mm_c ", "#mu#mu MET significance (c);MET [GeV]", 100, 0., 250.);
  hc_sMET_mm_l  = fs->make < TH1F > ("hc_sMET_mm_l ", "#mu#mu MET significance (dusg);MET [GeV]", 100, 0., 250.);
  hc_sMET_mm_tt = fs->make < TH1F > ("hc_sMET_mm_tt", "#mu#mu MET significance (Z#rightarrow#tau#tau);MET [GeV]", 100, 0., 250.);

  hc_CSV_ee    = fs->make < TH1F > ("hc_CSV_ee",    "ee CSV discriminant;CSV", 100, 0., 1.);
  hc_CSV_ee_b  = fs->make < TH1F > ("hc_CSV_ee_b",  "ee CSV discriminant (b);CSV", 100, 0., 1.);
  hc_CSV_ee_c  = fs->make < TH1F > ("hc_CSV_ee_c",  "ee CSV discriminant (c);CSV", 100, 0., 1.);
  hc_CSV_ee_l  = fs->make < TH1F > ("hc_CSV_ee_l",  "ee CSV discriminant (dusg);CSV", 100, 0., 1.);
  hc_CSV_ee_tt = fs->make < TH1F > ("hc_CSV_ee_tt", "ee CSV discriminant (Z#rightarrow#tau#tau);CSV", 100, 0., 1.);

  hc_CSV_mm    = fs->make < TH1F > ("hc_CSV_mm",    "#mu#mu  CSV discriminant;CSV", 100, 0., 1.);
  hc_CSV_mm_b  = fs->make < TH1F > ("hc_CSV_mm_b",  "#mu#mu  CSV discriminant (b);CSV", 100, 0., 1.);
  hc_CSV_mm_c  = fs->make < TH1F > ("hc_CSV_mm_c",  "#mu#mu  CSV discriminant (c);CSV", 100, 0., 1.);
  hc_CSV_mm_l  = fs->make < TH1F > ("hc_CSV_mm_l",  "#mu#mu  CSV discriminant (dusg);CSV", 100, 0., 1.);
  hc_CSV_mm_tt = fs->make < TH1F > ("hc_CSV_mm_tt", "#mu#mu  CSV discriminant (Z#rightarrow#tau#tau);CSV", 100, 0., 1.);

  hc_CSV_all_ee    = fs->make < TH1F > ("hc_CSV_all_ee",    "ee CSV discriminant;CSV", 100, 0., 1.);
  hc_CSV_all_ee_b  = fs->make < TH1F > ("hc_CSV_all_ee_b",  "ee CSV discriminant (b);CSV", 100, 0., 1.);
  hc_CSV_all_ee_c  = fs->make < TH1F > ("hc_CSV_all_ee_c",  "ee CSV discriminant (c);CSV", 100, 0., 1.);
  hc_CSV_ee_b_split  = fs->make < TH1F > ("hc_CSV_ee_b_split",  "ee CSV discriminant (g->bb);CSV", 100, 0., 1.);
  hc_CSV_ee_c_split  = fs->make < TH1F > ("hc_CSV_ee_c_split",  "ee CSV discriminant (g->cc);CSV", 100, 0., 1.);
  hc_CSV_all_ee_l  = fs->make < TH1F > ("hc_CSV_all_ee_l",  "ee CSV discriminant (dusg);CSV", 100, 0., 1.);
  hc_CSV_all_ee_tt = fs->make < TH1F > ("hc_CSV_all_ee_tt", "ee CSV discriminant (Z#rightarrow#tau#tau);CSV", 100, 0., 1.);

  hc_CSV_all_mm    = fs->make < TH1F > ("hc_CSV_all_mm",    "#mu#mu  CSV discriminant;CSV", 100, 0., 1.);
  hc_CSV_all_mm_b  = fs->make < TH1F > ("hc_CSV_all_mm_b",  "#mu#mu  CSV discriminant (b);CSV", 100, 0., 1.);
  hc_CSV_all_mm_c  = fs->make < TH1F > ("hc_CSV_all_mm_c",  "#mu#mu  CSV discriminant (c);CSV", 100, 0., 1.);
  hc_CSV_mm_b_split  = fs->make < TH1F > ("hc_CSV_mm_b_split",  "#mu#mu CSV discriminant (g->bb);CSV", 100, 0., 1.);
  hc_CSV_mm_c_split  = fs->make < TH1F > ("hc_CSV_mm_c_split",  "#mu#mu CSV discriminant (g->cc);CSV", 100, 0., 1.);
  hc_CSV_all_mm_l  = fs->make < TH1F > ("hc_CSV_all_mm_l",  "#mu#mu  CSV discriminant (dusg);CSV", 100, 0., 1.);
  hc_CSV_all_mm_tt = fs->make < TH1F > ("hc_CSV_all_mm_tt", "#mu#mu  CSV discriminant (Z#rightarrow#tau#tau);CSV", 100, 0., 1.);

  hc_BJP_ee    = fs->make < TH1F > ("hc_BJP_ee",    "ee BJP discriminant;BJP", 100, 0., 10.);
  hc_BJP_ee->Sumw2();
  hc_BJP_ee_b  = fs->make < TH1F > ("hc_BJP_ee_b",  "ee BJP discriminant (b);BJP", 100, 0., 10.);
  hc_BJP_ee_b->Sumw2();
  hc_BJP_ee_c  = fs->make < TH1F > ("hc_BJP_ee_c",  "ee BJP discriminant (c);BJP", 100, 0., 10.);
  hc_BJP_ee_c->Sumw2();
  hc_BJP_ee_l  = fs->make < TH1F > ("hc_BJP_ee_l",  "ee BJP discriminant (dusg);BJP", 100, 0., 10.);
  hc_BJP_ee_l->Sumw2();
  hc_BJP_ee_tt = fs->make < TH1F > ("hc_BJP_ee_tt", "ee BJP discriminant (Z#rightarrow#tau#tau);BJP", 100, 0., 10.);
  hc_BJP_ee_tt->Sumw2();

  hc_BJP_mm    = fs->make < TH1F > ("hc_BJP_mm",    "#mu#mu  BJP discriminant;BJP", 100, 0., 10.);
  hc_BJP_mm->Sumw2();
  hc_BJP_mm_b  = fs->make < TH1F > ("hc_BJP_mm_b",  "#mu#mu  BJP discriminant (b);BJP", 100, 0., 10.);
  hc_BJP_mm_b->Sumw2();
  hc_BJP_mm_c  = fs->make < TH1F > ("hc_BJP_mm_c",  "#mu#mu  BJP discriminant (c);BJP", 100, 0., 10.);
  hc_BJP_mm_c->Sumw2();
  hc_BJP_mm_l  = fs->make < TH1F > ("hc_BJP_mm_l",  "#mu#mu  BJP discriminant (dusg);BJP", 100, 0., 10.);
  hc_BJP_mm_l->Sumw2();
  hc_BJP_mm_tt = fs->make < TH1F > ("hc_BJP_mm_tt", "#mu#mu  BJP discriminant (Z#rightarrow#tau#tau);BJP", 100, 0., 10.);
  hc_BJP_mm_tt->Sumw2();

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

  bool is_tautau = false;
  unsigned int n_b = 0;
  unsigned int n_c = 0;
  unsigned int n_b_split = 0;
  unsigned int n_c_split = 0;

  double MyWeight = 1.;



  // ================================================================================================
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

  
  // ================================================================================================
  //  Event weights


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



  // ================================================================================================
  //  Electrons

  vector < pat::Electron > vect_elem;
  vector < pat::Electron > vect_elep;

  for (pat::ElectronCollection::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele) {

    if ( ele->pt()<20. || fabs(ele->eta())>2.4 ) continue; 

    ele->charge()<0. ? vect_elem.push_back(*ele) : vect_elep.push_back(*ele);
    
  } // ele loop


  double diele_mass = 0.;
  math::XYZTLorentzVector z_ee;

  if ( vect_elem.size()>0. && vect_elep.size()>0. ){

    z_ee = vect_elem[0].p4() + vect_elep[0].p4();
    diele_mass = z_ee.mass();
    if ( diele_mass>71. && diele_mass<111. ) ee_event = true;

  }




  // ================================================================================================
  //   Muons
  
  vector < pat::Muon > vect_muom;
  vector < pat::Muon > vect_muop;

  for (pat::MuonCollection::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon) {

    if ( muon->pt()<20. || fabs(muon->eta())>2.4 ) continue;

    muon->charge()<0. ? vect_muom.push_back(*muon) : vect_muop.push_back(*muon);

  } // muon loop


  double dimuo_mass = 0.;
  math::XYZTLorentzVector z_mm;

  if ( vect_muom.size()>0. && vect_muop.size()>0. ){

    z_mm = vect_muom[0].p4() + vect_muop[0].p4();
    dimuo_mass = z_mm.mass();
    if ( dimuo_mass>71. && dimuo_mass<111. ) mm_event = true;

  }

  ee_event = ee_event && (lepton_ == "electron");
  mm_event = mm_event && (lepton_ == "muon");
  


  // Keep only Z-->ee and Z-->mm events:
  if ( !ee_event && !mm_event ) return;



  // --- Identify Z->tautau events: 

  if (isMC) {
    for (std::vector <reco::GenParticle>::const_iterator thepart = genPart->begin(); thepart != genPart->end(); thepart++) {
      if ((int) abs(thepart->pdgId()) == 23) {
        for (UInt_t i=0; i<thepart->numberOfDaughters(); i++){
	  if (abs(thepart->daughter(i)->pdgId()) == 15 && thepart->daughter(i)->status()==3){
	    is_tautau = true;
	  }
        }
      }
    }
  }



  // ================================================================================================
  //  Lepton scale factors

  if (isMC) {
    if (ee_event) {
      double scalFac_elem = ElSF_ ->Val(vect_elem[0].pt(), vect_elem[0].eta()) * 
	                    ElSF2_->Val(vect_elem[0].pt(), vect_elem[0].eta());
      double scalFac_elep = ElSF_ ->Val(vect_elep[0].pt(), vect_elep[0].eta()) * 
                            ElSF2_->Val(vect_elep[0].pt(), vect_elep[0].eta());
      MyWeight *= scalFac_elem * scalFac_elep;
    }
    if (mm_event) {
      double scalFac_muom = MuSF_->Val(vect_muom[0].pt(), vect_muom[0].eta()) * 
	                    sqrt(MuSF2_->Val(fabs(vect_muom[0].eta()), fabs(vect_muop[0].eta())));
      double scalFac_muop = MuSF_->Val(vect_muop[0].pt(), vect_muop[0].eta()) * 
	                    sqrt(MuSF2_->Val(fabs(vect_muom[0].eta()), fabs(vect_muop[0].eta())));
      MyWeight *= scalFac_muom * scalFac_muop;
    }
  }



  // ================================================================================================
  //  Primary vertices

  bool vtx_cut = false;

  edm::Handle < vector < reco::Vertex > > vertices;
  iEvent.getByLabel (edm::InputTag ("goodOfflinePrimaryVertices"), vertices);

  if (vertices->size() > 0) {
    const reco::Vertex* theVertex = &(vertices->front());
    if ( theVertex->ndof()                 >  5 &&
         fabs(theVertex->z())              < 24.&&
         fabs(theVertex->position().rho()) <  2.) vtx_cut = true;
  } 


  // Keep only events with a good primary vertex:
  if ( !vtx_cut ) return;



  // ================================================================================================
  //  Jets

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
    if (isMC && jet->genJet()) jerCorr = jetResolutionCorrection(jet->eta(), jet->pt(), jet->genJet()->pt(), par2_);
    
    pat::Jet jetNew = (*jet);
    math::XYZTLorentzVector jetNew_p4 = jetNew.p4();

    jetNew_p4 = jetNew_p4 * (1.0 + jecUnc * par_) * jerCorr;

    jetNew.setP4(jetNew_p4);


    // --- jet selection:
    if ( jetNew.pt()<30. || fabs(jetNew.eta())>2.4  ) continue;
    vect_jets.push_back(jetNew);


    // --- c tagging:
    double discrCSV = jetNew.bDiscriminator("combinedSecondaryVertexBJetTags");


    if ( ( discrCSV>0.244 && discrCSV<0.898 && fabs(jetNew.eta())<2.4 ) )
      vect_cjets.push_back(jetNew);


    // --- check the MC truth:
    bool hasb = false;
    bool hasc = false;
    if ( isMC && jet->genJet() ) {

      vector <const reco::GenParticle*> genp =  jet->genJet()->getGenConstituents();
      for (unsigned int ipart=0; ipart<genp.size(); ++ipart){

	const reco::Candidate* mot = genp[ipart]->mother(0);

	while (mot){
	  if ( abs(mot->pdgId()) == 5 ) hasb = true;
	  if ( abs(mot->pdgId()) == 4 ) hasc = true;
	  mot = mot->mother(0);
	}

      }


      // --- Fill histograms to calculate the CSV tagging efficiencies
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


    if ( hasb ) {
      if ( fabs(jetNew.partonFlavour())==5 ){
	has_b.push_back(1);
	n_b++;
      }
      else {
	has_b.push_back(2);
	n_b_split++;
      }
    }
    else  has_b.push_back(0);

    if ( hasc ) {
      if ( fabs(jetNew.partonFlavour())==4 ){
	has_c.push_back(1);
	n_c++;
      }
      else {
	has_c.push_back(2);
	n_c_split++;
      }
    }
    else  has_c.push_back(0);

  } // for jet

 

  // Keep only events with at least one reconstructed jet:
  if ( vect_jets.size()==0 ) return;

  double scalFac_b = ctagWeight(isMC, vect_jets);
  


  if ( ee_event ){
    h_MET_ee->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
    h_sMET_ee->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
    
    if (isMC) {
     
      if (is_tautau) {
	h_MET_ee_tt->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
	h_sMET_ee_tt->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
      }
      else {
	if ( n_b>0 ){
	  h_MET_ee_b->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
	  h_sMET_ee_b->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
 	}
	else if ( n_b==0 && n_c>0 ){
	  h_MET_ee_c->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
	  h_sMET_ee_c->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
	}
	else if ( n_b==0 && n_c==0 ){
	  h_MET_ee_l->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
	  h_sMET_ee_l->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
	}

      } 
    }


    if ( vect_cjets.size()>0 ) {

      hc_MET_ee->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
      hc_sMET_ee->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);

      if (isMC) {
	
	if (is_tautau) {
	  hc_MET_ee_tt->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	  hc_sMET_ee_tt->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	}
	else {
	  if ( n_b>0 ){
	    hc_MET_ee_b->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	    hc_sMET_ee_b->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	  }
	  else if ( n_b==0 && n_c>0 ){
	    hc_MET_ee_c->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	    hc_sMET_ee_c->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	  }
	  else if ( n_b==0 && n_c==0 ){
	    hc_MET_ee_l->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	    hc_sMET_ee_l->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	  }

	} 
      }


    }

  }
  else if ( mm_event ){
    h_MET_mm->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
    h_sMET_mm->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);

    if (isMC) {
     
      if (is_tautau) {
	h_MET_mm_tt->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
	h_sMET_mm_tt->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
      }
      else {
	if ( n_b>0 ){
	  h_MET_mm_b->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
	  h_sMET_mm_b->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
 	}
	else if ( n_b==0 && n_c>0 ){
	  h_MET_mm_c->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
	  h_sMET_mm_c->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
	}
	else if ( n_b==0 && n_c==0 ){
	  h_MET_mm_l->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight);
	  h_sMET_mm_l->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight);
	}

      } 
    }


    if ( vect_cjets.size()>0 ) {

      hc_MET_mm->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
      hc_sMET_mm->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);

      if (isMC) {
	
	if (is_tautau) {
	  hc_MET_mm_tt->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	  hc_sMET_mm_tt->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	}
	else {
	  if ( n_b>0 ){
	    hc_MET_mm_b->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	    hc_sMET_mm_b->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	  }
	  else if ( n_b==0 && n_c>0 ){
	    hc_MET_mm_c->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	    hc_sMET_mm_c->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	  }
	  else if ( n_b==0 && n_c==0 ){
	    hc_MET_mm_l->Fill(met->empty() ? 0. : (*met)[0].et(), MyWeight*scalFac_b);
	    hc_sMET_mm_l->Fill(met->empty() ? 0. : (*met)[0].significance(), MyWeight*scalFac_b);
	  }

	} 
      }


    }

  }

  // Keep only events with low MET significance
  if ( !met->empty() && (*met)[0].significance() > 30. ) return; 



  // ================================================================================================
  //  Z+jet histograms


  if ( ee_event ){

    if (!is_tautau){

      n_events[0][0]++;
      w_events[0][0] += MyWeight;
      w2_events[0][0] += MyWeight*MyWeight;

      h_M_ee->Fill(diele_mass, MyWeight);
      h_njet_ee->Fill(vect_jets.size(), MyWeight);

      h_CSV_ee->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      h_BJP_ee->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
      for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	h_CSV_all_ee->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
    }

    if (isMC) {
     
      if (is_tautau) {

	n_ztautau[0][0]++;
	w_ztautau[0][0] += MyWeight;
	w2_ztautau[0][0] += MyWeight*MyWeight;

	h_M_ee_tt->Fill(diele_mass, MyWeight);
	h_njet_ee_tt->Fill(vect_jets.size(), MyWeight);

	h_CSV_ee_tt->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	h_BJP_ee_tt->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
	for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	  h_CSV_all_ee_tt->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	
      }
      else {

	//if ( n_b>0 && n_c==0 ){
	if ( n_b>0 ){
      
	  h_M_ee_b->Fill(diele_mass, MyWeight);
	  h_njet_ee_b->Fill(vect_jets.size(), MyWeight);

	  h_CSV_ee_b->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	  h_BJP_ee_b->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    h_CSV_all_ee_b->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      
	}
	else if ( n_b==0 && n_c>0 ){
      
	  h_M_ee_c->Fill(diele_mass, MyWeight);
	  h_njet_ee_c->Fill(vect_jets.size(), MyWeight);

	  h_CSV_ee_c->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	  h_BJP_ee_c->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    h_CSV_all_ee_c->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      
	}
	else if ( n_b==0 && n_c==0 ){
      
	  h_M_ee_l->Fill(diele_mass, MyWeight);
	  h_njet_ee_l->Fill(vect_jets.size(), MyWeight);

	  h_CSV_ee_l->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	  h_BJP_ee_l->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    h_CSV_all_ee_l->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      
	}

	if ( n_b_split>0 ){
	  h_CSV_ee_b_split->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	}
	if ( n_c_split>0 ){
	  h_CSV_ee_c_split->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	}

      }
      
    }

  }


  if ( mm_event ){
    
    if (!is_tautau){

      n_events[0][1]++;
      w_events[0][1] += MyWeight;
      w2_events[0][1] += MyWeight*MyWeight;

      h_M_mm->Fill(dimuo_mass, MyWeight);
      h_njet_mm->Fill(vect_jets.size(), MyWeight);

      h_CSV_mm->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      h_BJP_mm->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
      for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	h_CSV_all_mm->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);

    }

    if (isMC) {
     
      if (is_tautau) {

	n_ztautau[0][1]++;
	w_ztautau[0][1] += MyWeight;
	w2_ztautau[0][1] += MyWeight*MyWeight;

	h_M_mm_tt->Fill(dimuo_mass, MyWeight);
	h_njet_mm_tt->Fill(vect_jets.size(), MyWeight);

	h_CSV_mm_tt->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	h_BJP_mm_tt->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
	for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	  h_CSV_all_mm_tt->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	
      }
      else {

	//if ( n_b>0 && n_c==0 ){
	if ( n_b>0 ){
      
	  h_M_mm_b->Fill(dimuo_mass, MyWeight);
	  h_njet_mm_b->Fill(vect_jets.size(), MyWeight);

	  h_CSV_mm_b->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	  h_BJP_mm_b->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    h_CSV_all_mm_b->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      
	}
	else if ( n_b==0 && n_c>0 ){
      
	  h_M_mm_c->Fill(dimuo_mass, MyWeight);
	  h_njet_mm_c->Fill(vect_jets.size(), MyWeight);

	  h_CSV_mm_c->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	  h_BJP_mm_c->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    h_CSV_all_mm_c->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      
	}
	else if ( n_b==0 && n_c==0 ){
      
	  h_M_mm_l->Fill(dimuo_mass, MyWeight);
	  h_njet_mm_l->Fill(vect_jets.size(), MyWeight);

	  h_CSV_mm_l->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	  h_BJP_mm_l->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    h_CSV_all_mm_l->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
      
	}

	if ( n_b_split>0 ){
	  h_CSV_mm_b_split->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	}
	if ( n_c_split>0 ){
	  h_CSV_mm_c_split->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight);
	}


      }
      
    }

  }


  // ================================================================================================
  //  Z+c histograms


  // Keep only events with at least one reconstructed c-tagged jet:
  if ( vect_cjets.size()==0 ) return;


  if ( ee_event ){

    if (!is_tautau){

      n_events[1][0]++;
      w_events[1][0] += MyWeight*scalFac_b;
      w2_events[1][0] += MyWeight*MyWeight*scalFac_b*scalFac_b;

      hc_M_ee->Fill(diele_mass, MyWeight*scalFac_b);
      hc_njet_ee->Fill(vect_jets.size(), MyWeight*scalFac_b);

      hc_CSV_ee->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      hc_BJP_ee->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
      for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	hc_CSV_all_ee->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);

    }

    if (isMC) {
     
      if (is_tautau) {

	n_ztautau[1][0]++;
	w_ztautau[1][0] += MyWeight*scalFac_b;
	w2_ztautau[1][0] += MyWeight*MyWeight*scalFac_b;

	hc_M_ee_tt->Fill(diele_mass, MyWeight*scalFac_b);
	hc_njet_ee_tt->Fill(vect_jets.size(), MyWeight*scalFac_b);

	hc_CSV_ee_tt->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	hc_BJP_ee_tt->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
	for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	  hc_CSV_all_ee_tt->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	
      }
      else {

	//if ( n_b>0 && n_c==0 ){
	if ( n_b>0 ){
      
	  hc_M_ee_b->Fill(diele_mass, MyWeight*scalFac_b);
	  hc_njet_ee_b->Fill(vect_jets.size(), MyWeight*scalFac_b);

	  hc_CSV_ee_b->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	  hc_BJP_ee_b->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    hc_CSV_all_ee_b->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      
	}
	else if ( n_b==0 && n_c>0 ){
      
	  hc_M_ee_c->Fill(diele_mass, MyWeight*scalFac_b);
	  hc_njet_ee_c->Fill(vect_jets.size(), MyWeight*scalFac_b);

	  hc_CSV_ee_c->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	  hc_BJP_ee_c->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    hc_CSV_all_ee_c->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      
	}
	else if ( n_b==0 && n_c==0 ){
      
	  hc_M_ee_l->Fill(diele_mass, MyWeight*scalFac_b);
	  hc_njet_ee_l->Fill(vect_jets.size(), MyWeight*scalFac_b);

	  hc_CSV_ee_l->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	  hc_BJP_ee_l->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    hc_CSV_all_ee_l->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      
	}

	if ( n_b_split>0 ){
	  hc_CSV_ee_b_split->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	}
	if ( n_c_split>0 ){
	  hc_CSV_ee_c_split->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	}


      }
      
    }

  }


  if ( mm_event ){

    if (!is_tautau){

      n_events[1][1]++;
      w_events[1][1] += MyWeight*scalFac_b;
      w2_events[1][1] += MyWeight*MyWeight*scalFac_b*scalFac_b;

      hc_M_mm->Fill(dimuo_mass, MyWeight*scalFac_b);
      hc_njet_mm->Fill(vect_jets.size(), MyWeight*scalFac_b);

      hc_CSV_mm->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      hc_BJP_mm->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
      for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	hc_CSV_all_mm->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);

    }

    if (isMC) {
     
      if (is_tautau) {

	n_ztautau[1][1]++;
	w_ztautau[1][1] += MyWeight*scalFac_b;
	w2_ztautau[1][1] += MyWeight*MyWeight*scalFac_b;

	hc_M_mm_tt->Fill(dimuo_mass, MyWeight*scalFac_b);
	hc_njet_mm_tt->Fill(vect_jets.size(), MyWeight*scalFac_b);

	hc_CSV_mm_tt->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	hc_BJP_mm_tt->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
	for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	  hc_CSV_all_mm_tt->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	
      }
      else {

	//if ( n_b>0 && n_c==0 ){
	if ( n_b>0 ){
      
	  hc_M_mm_b->Fill(dimuo_mass, MyWeight*scalFac_b);
	  hc_njet_mm_b->Fill(vect_jets.size(), MyWeight*scalFac_b);

	  hc_CSV_mm_b->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	  hc_BJP_mm_b->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    hc_CSV_all_mm_b->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      
	}
	else if ( n_b==0 && n_c>0 ){
      
	  hc_M_mm_c->Fill(dimuo_mass, MyWeight*scalFac_b);
	  hc_njet_mm_c->Fill(vect_jets.size(), MyWeight*scalFac_b);

	  hc_CSV_mm_c->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	  hc_BJP_mm_c->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    hc_CSV_all_mm_c->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      
	}
	else if ( n_b==0 && n_c==0 ){
      
	  hc_M_mm_l->Fill(dimuo_mass, MyWeight*scalFac_b);
	  hc_njet_mm_l->Fill(vect_jets.size(), MyWeight*scalFac_b);

	  hc_CSV_mm_l->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	  hc_BJP_mm_l->Fill(vect_jets[0].bDiscriminator("jetBProbabilityBJetTags"), MyWeight*scalFac_b);
	  for (unsigned int ijet=0; ijet<vect_jets.size();++ijet)
	    hc_CSV_all_mm_l->Fill(vect_jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
      
	}
	
	if ( n_b_split>0 ){
	  hc_CSV_mm_b_split->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	}
	if ( n_c_split>0 ){
	  hc_CSV_mm_c_split->Fill(vect_jets[0].bDiscriminator("combinedSecondaryVertexBJetTags"), MyWeight*scalFac_b);
	}

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
    for (int j=0; j<2; ++j){
      n_ztautau[i][j]  = 0 ;
      w_ztautau[i][j]  = 0.;
      w2_ztautau[i][j] = 0.;
      n_events[i][j]   = 0 ;
      w_events[i][j]   = 0.;
      w2_events[i][j]  = 0.;
    }
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
  std::cout << "Z-->tautau yield (inclusive) = " 
	    << n_ztautau[0][0] << " " 
	    << n_ztautau[0][1] << std::endl;
  std::cout << "Z-->tautau yield (tagged)    = " 
	    << n_ztautau[1][0] << " " 
	    << n_ztautau[1][1] << std::endl;
  std::cout << "Z-->tautau weighted yield (inclusive) = "
	    << w_ztautau[0][0] << " " << w2_ztautau[0][0] << " " 
	    << w_ztautau[0][1] << " " << w2_ztautau[0][1] << std::endl;
  std::cout << "Z-->tautau weighted yield (tagged)    = " 
	    << w_ztautau[1][0] << " " << w2_ztautau[1][0] << " " 
	    << w_ztautau[1][1] << " " << w2_ztautau[1][1] << std::endl;

  std::cout <<  std::endl;

  std::cout << "Yield (inclusive) = " 
	    << n_events[0][0] << " " 
	    << n_events[0][1] << std::endl;
  std::cout << "Yield (tagged)    = " 
	    << n_events[1][0] << " " 
	    << n_events[1][1] << std::endl;
  std::cout << "Weighted yield (inclusive) = " 
	    << w_events[0][0] << " " << w2_events[0][0] << " " 
	    << w_events[0][1] << " " << w2_events[0][1] << std::endl;
  std::cout << "Weighted yield (tagged)    = "  
	    << w_events[1][0] << " " << w2_events[1][0] << " " 
	    << w_events[1][1] << " " << w2_events[1][1] << std::endl;
  


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
