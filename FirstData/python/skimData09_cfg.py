import FWCore.ParameterSet.Config as cms
process = cms.Process("Skim")
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/5C95707E-15ED-DE11-A7D6-0015178C646C.root",#run 123970 lumi 13
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/A4FDFAFF-26ED-DE11-90B7-0015178C6740.root", #run 124020 lumi 18
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/52165DFD-19ED-DE11-9679-00151796D6F0.root", #run 124020 lumi 30
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/E8D6B236-25ED-DE11-9826-00151796D688.root",#run124022 lumi 87
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/F65182E5-13ED-DE11-99CC-0024E876841F.root",#run124022 129
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/B872787E-15ED-DE11-A77A-00151796D86C.root",#run124022 163
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/0026565C-2CED-DE11-A3FC-001D0967DD2D.root",#run124023 lumi 40
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/E89FEAFF-23ED-DE11-9F78-0024E8768265.root",#run124030 lumi 7
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/3AB79D38-1FED-DE11-9673-0015178C64BC.root",#run124230 lumi 34
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/223763C6-18ED-DE11-89A4-001D0967D49F.root",#lumi 14&19
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/94AB71E7-1BED-DE11-AE0C-001D0967D314.root" #lumi 52 
    ),
                            eventsToProcess = cms.untracked.VEventRange(
    cms.EventRange(123970,3196130,123970,3196130),#lumi 13
    cms.EventRange(124020,6180262,124020,6180262),#lumi 18
    cms.EventRange(124020,10528700,124020,10528700),#lumi30
    cms.EventRange(124022,13702620,124022,13702620),#lumi 87
    cms.EventRange(124022,29116624,124022,29116624),#lumi 129
    cms.EventRange(124022,41867601,124022,41867601),#lumi 163
    cms.EventRange(124023,3216398,124023,3216398),#lumi 40
    cms.EventRange(124030,2342390,124030,2342390),#lumi 7
    cms.EventRange(124230,11845706,124230,11845706),#lumi 34
    cms.EventRange(124120,4207354,124120,4207354),#lumi14
    cms.EventRange(124120,5686693,124120,5686693),#lumi19
    cms.EventRange(124120,15957076,124120,15957076)#lumi52
    ) )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("dimuon_skim.root"),
                               outputCommands = cms.untracked.vstring("keep *",
                               "drop *_MEtoEDMConverter_*_*"))
process.o = cms.EndPath(process.out)
