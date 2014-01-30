import FWCore.ParameterSet.Config as cms

process = cms.Process("ZcAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = cms.untracked.string("WARNING")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
	# replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring(
                #'file:patTuple_1_1_IPe.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1000_1_HGo.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1001_2_cC8.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1002_2_rof.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1003_1_32u.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1004_2_JKa.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1005_2_a1a.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1006_2_smy.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1007_2_Hvj.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1008_1_bEp.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1009_2_PMJ.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1010_2_EMz.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1011_1_dAb.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1012_2_cJT.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1013_2_2yI.root',
                'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DYJetsToLL/patTuple_1014_2_nRc.root',
                #
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/TTbar/patTuple_81_1_C9w.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/TTbar/patTuple_82_2_NMq.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/TTbar/patTuple_83_2_cyK.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/TTbar/patTuple_84_2_FH1.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/TTbar/patTuple_85_1_eTB.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/TTbar/patTuple_86_1_8kl.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/TTbar/patTuple_87_1_VEi.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/TTbar/patTuple_88_2_KvV.root',
                #
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_941_2_xru.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_942_1_D4N.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_943_1_UsT.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_944_2_oTk.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_945_1_gG9.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_946_2_N0Y.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_947_2_CDE.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_948_2_19n.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_949_2_nBu.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_94_1_tD0.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_950_2_Ton.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_951_1_d2t.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_952_2_Ib4.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_953_1_dO2.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_954_1_2Wt.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_955_1_ong.root',
                #'file:/gpfs/grid/srm/cms/store/user/vieri/grid/v11/DoubleElectron_2012D_22Jan13/patTuple_956_1_K5J.root',
        )
)


process.anaEle = cms.EDProducer('ZcAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton = cms.untracked.string("electron")
)

process.anaMuo = cms.EDProducer('ZcAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton = cms.untracked.string("muon")
)

process.anaEleMuo = cms.EDProducer('ZcAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em"),
        lepton = cms.untracked.string("electron+muon")
)


process.TFileService = cms.Service("TFileService",
        fileName = cms.string('ZcAnalyzer.root')
)

process.p = cms.Path(process.anaEle*process.anaMuo*process.anaEleMuo)
