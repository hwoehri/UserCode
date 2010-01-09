import FWCore.ParameterSet.Config as cms
process = cms.Process("Skim")
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/5C95707E-15ED-DE11-A7D6-0015178C646C.root",#run 123970 lumi 13
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/12A1BCFD-19ED-DE11-B316-001D0967DF62.root",#run 123970 lumi 18
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/A4FDFAFF-26ED-DE11-90B7-0015178C6740.root",#run 124020 lumi 18
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/52165DFD-19ED-DE11-9679-00151796D6F0.root",#run 124020 lumi 30
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/045412FE-19ED-DE11-987D-001D0967DEEF.root",#run 124020 lumi 47
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/76A9B1E5-1BED-DE11-BABE-0024E86E8D66.root",#run 124020 lumi 70
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/AE9A7219-26ED-DE11-8C62-0015178C0170.root",#run 124020 lumi 72
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/E8D6B236-25ED-DE11-9826-00151796D688.root",#run124022 lumi 87
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/F65182E5-13ED-DE11-99CC-0024E876841F.root",#run124022 lumi 129
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/6EF93111-23ED-DE11-8EEB-0015178C49F8.root",#run124022 lumi 130
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/B872787E-15ED-DE11-A77A-00151796D86C.root",#run124022 lumi 163
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/CA8CDBE7-12ED-DE11-963C-00151796D894.root",#run124022 lumi 182
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/0026565C-2CED-DE11-A3FC-001D0967DD2D.root",#run124023 lumi 40
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/A462B136-25ED-DE11-AED2-0015178C1574.root",#run124023 lumi 47
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/1E3C226E-1AED-DE11-9987-001D0967D6E8.root",#run124023 lumi 56
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/826BF118-26ED-DE11-88A6-001D0967DA49.root",#run124023 lumi 91
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/5A4F5BC4-18ED-DE11-9312-001D0967DAE4.root",#run124024 lumi 18
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/0A0EFE88-10ED-DE11-9349-00151796D674.root",#run124027 lumi 37
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/E89FEAFF-23ED-DE11-9F78-0024E8768265.root",#run124030 lumi 6
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/E89FEAFF-23ED-DE11-9F78-0024E8768265.root",#run124030 lumi 7
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/D6FC56AE-ADEE-DE11-B2F7-0024E87680E7.root",#run124030 lumi 26
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/3AB79D38-1FED-DE11-9673-0015178C64BC.root",#run124230 lumi 34
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/62000402-24ED-DE11-AE55-001D0967D49F.root",#run124230 lumi 40
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/9494BB88-0BED-DE11-837E-001D0967D6C0.root",#run124230 lumi 68
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/A674A9FE-19ED-DE11-8C16-00151785FF78.root",#run124120 lumi 5 (2.36 TeV)
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/223763C6-18ED-DE11-89A4-001D0967D49F.root",#run124120 lumi 14&19 (2.36 TeV)
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/223763C6-18ED-DE11-89A4-001D0967D49F.root",#run124120 lumi 20 (2.36 TeV)
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/B4237151-29ED-DE11-81ED-0015178C1804.root",#run124120 lumi 51 (2.36 TeV)
    "/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/94AB71E7-1BED-DE11-AE0C-001D0967D314.root" #run124120 lumi 52 (2.36 TeV)
    ),
                            eventsToProcess = cms.untracked.VEventRange(
    cms.EventRange(123970,3196130,123970,3196130),#lumi 13
    cms.EventRange(124020,4922604,124020,4922604),#lumi 18
    cms.EventRange(124020,6180262,124020,6180262),#lumi 18
    cms.EventRange(124020,10528700,124020,10528700),#lumi30
    cms.EventRange(124020,17047394,124020,17047394),#lumi47
    cms.EventRange(124022,8126826,124022,8126826),#lumi 70
    cms.EventRange(124022,8834000,124022,8834000),#lumi 72
    cms.EventRange(124022,13702620,124022,13702620),#lumi 87
    cms.EventRange(124022,29116624,124022,29116624),#lumi 129
    cms.EventRange(124022,29692749,124022,29692749),#lumi 130
    cms.EventRange(124022,41867601,124022,41867601),#lumi 163
    cms.EventRange(124022,49182951,124022,49182951),#lumi 182
    cms.EventRange(124023,3216398,124023,3216398),#lumi 40
    cms.EventRange(124023,5717745,124023,5717745),#lumi 47
    cms.EventRange(124023,9160822,124023,9160822),#lumi 56
    cms.EventRange(124023,22468285,124023,22468285),#lumi 91
    cms.EventRange(124024,6436943,124024,6436943),#lumi 18
    cms.EventRange(124027,12941332,124027,12941332),#lumi 37
    cms.EventRange(124030,1985373,124030,1985373),#lumi 6
    cms.EventRange(124030,2342390,124030,2342390),#lumi 7
    cms.EventRange(124230,8638137,124230,8638137),#lumi 26
    cms.EventRange(124230,11845706,124230,11845706),#lumi 34
    cms.EventRange(124230,14098140,124230,14098140),#lumi 40
    cms.EventRange(124230,24268447,124230,24268447),#lumi 68

    cms.EventRange(124120,1320824,124120,1320824),#lumi 5
    cms.EventRange(124120,4207354,124120,4207354),#lumi 14
    cms.EventRange(124120,5686693,124120,5686693),#lumi 19
    cms.EventRange(124120,6094468,124120,6094468),#lumi 20
    cms.EventRange(124120,15823906,124120,15823906),#lumi 51
    cms.EventRange(124120,15957076,124120,15957076)#lumi 52
    ) )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("dimuon_skim_allSignComb.root"),
                               outputCommands = cms.untracked.vstring("keep *",
                               "drop *_MEtoEDMConverter_*_*"))
process.o = cms.EndPath(process.out)
