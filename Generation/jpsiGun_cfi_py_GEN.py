# Auto generated configuration file
# using: 
# Revision: 1.222.2.6 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/jpsiGun_cfi.py -s GEN --conditions DESIGN_38_V13::All --datatier GEN-SIM-RAW --eventcontent RAWSIM -n 100 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.VtxSmearedRealistic7TeVCollision_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.222.2.6 $'),
    annotation = cms.untracked.string('Configuration/GenProduction/python/jpsiGun_cfi.py nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("EmptySource")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.RandomNumberGeneratorService.generator.initialSeed = 12345

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('/tmp/hwoehri/jpsiGun_cfi_py_GEN.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring('generation_step')
            SelectEvents = cms.vstring('genOniaFilter')
    )
)

# Additional output definition
process.makeGenTree = cms.EDAnalyzer('MakeGenTree',
                                  genParticlesLabel = cms.InputTag("generator"),
                                  QQbarPDG = cms.int32(443),
                                  outputFileName= cms.string("/tmp/hwoehri/jpsiGun_Tree.root")
)

# Other statements
process.GlobalTag.globaltag = 'DESIGN_38_V13::All'
process.generator = cms.EDProducer("Pythia6PtGun",
    PGunParameters = cms.PSet(
        MinPhi = cms.double(-3.14159265359),
        MinPt = cms.double(0.0),
        ParticleID = cms.vint32(443),
        MaxEta = cms.double(6.),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-6.),
        AddAntiParticle = cms.bool(False),
        MaxPt = cms.double(35.0)
    ),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    PythiaParameters = cms.PSet(
        pythiaJpsiDecays = cms.vstring('MDME(858,1)= 0    ! e- e+', 
            'MDME(859,1)= 1    ! mu- mu+', 
            'MDME(860,1)= 0    ! rndmflav        rndmflavbar'),
        parameterSets = cms.vstring('pythiaJpsiDecays')
    )
)

#onia filter
process.oniafilter = cms.EDFilter("PythiaFilter",
    MaxEta = cms.untracked.double(10.),
    MaxRapidity = cms.untracked.double(2.5),
    Status = cms.untracked.int32(2),
    MinRapidity = cms.untracked.double(-2.5),
    MinEta = cms.untracked.double(-10.),
    MinPt = cms.untracked.double(0.0),
    ParticleID = cms.untracked.int32(443)
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genOniaFilter = cms.Path(process.generator*process.oniafilter)
process.genAndWrite = cms.Path(process.generator*process.oniafilter*process.makeGenTree)
process.endjob_step = cms.Path(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.generation_step,process.endjob_step,process.RAWSIMoutput_step)
#process.schedule = cms.Schedule(process.genOniaFilter,process.endjob_step,process.RAWSIMoutput_step)
process.schedule = cms.Schedule(process.genAndWrite,process.endjob_step)

# special treatment in case of production filter sequence
for path in process.paths: 
    getattr(process,path)._seq = process.generator*getattr(process,path)._seq
