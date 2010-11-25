import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/tmp/hwoehri/jpsiGun_py_GEN.root'
    )
)

process.demo = cms.EDAnalyzer('MakeGenTree'
                              genParticlesLabel = cms.InputTag("generator"),
                              QQbarPDG = cms.int32(443),
                              outputFileName= cms.string("/tmp/hwoehri/jpsiGun_Tree.root")
))


process.p = cms.Path(process.demo)
