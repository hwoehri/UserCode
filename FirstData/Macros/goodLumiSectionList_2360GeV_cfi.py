import FWCore.ParameterSet.Config as cms

# filter on lumisections
goodLumisToProcess = cms.untracked.VLuminosityBlockRange('124120:2-124120:58',
                                                         '124275:2-124275:31')
