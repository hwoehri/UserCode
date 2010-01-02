import FWCore.ParameterSet.Config as cms

process = cms.Process('ANA')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("PhysicsTools.HepMCCandAlgos.allMuonsGenParticlesMatch_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

# Input source
process.source = cms.Source("PoolSource",
                            
    fileNames = cms.untracked.vstring(
    # first data
#    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/F682F559-B9ED-DE11-8057-002618943900.root'
#    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/28BD9294-E0F0-DE11-8961-0024E86E8D18.root'
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/FCFEFB8D-E0F0-DE11-A01D-0015178C6B88.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/F8269247-E0F0-DE11-9431-00151796D9C4.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/F4E10D24-E0F0-DE11-870D-0024E8767DAD.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/EAFCDCA2-E0F0-DE11-AF3D-0015178C6B8C.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/E2F75B8B-E0F0-DE11-8435-0024E8768D34.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/DE4AFC37-E0F0-DE11-A291-001D0967DC3D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/DE306D68-E0F0-DE11-989B-00151796C0A0.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/DCDDF98A-E0F0-DE11-BBEB-0015178C4C44.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/DC90039C-E0F0-DE11-930D-0024E87681F0.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/D88AC68D-E0F0-DE11-BF91-0024E87683EB.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/D043639A-E0F0-DE11-9B13-0024E87699C0.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/C8E94A79-E0F0-DE11-B075-0024E8768BFC.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/C0516464-E0F0-DE11-B85F-0024E8768C8B.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/BEE8747D-E0F0-DE11-A3F6-001D0967DC3D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/B891226E-E0F0-DE11-BE6E-0024E87699C0.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/B040DA9B-E0F0-DE11-B40B-0024E8767DAD.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/A2D7E02F-E0F0-DE11-9948-001D0967DA76.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/9C753B29-E0F0-DE11-BDB0-0024E87680A6.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/8E7E2337-E0F0-DE11-8371-001D096B0BDE.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/8E79E01F-E0F0-DE11-8B13-0024E87680DA.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/88B5A966-E0F0-DE11-A1EF-001D0967D8EB.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/88B07018-E0F0-DE11-9EF0-0024E8769999.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/86029E8B-E0F0-DE11-9C7D-00151796D814.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/7E372AAF-E0F0-DE11-9570-0024E8769999.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/70994874-E0F0-DE11-8177-0024E8768826.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/707292A1-E0F0-DE11-9121-0024E876A83B.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/68BC1717-E0F0-DE11-B533-0024E8767D79.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/649DF515-E0F0-DE11-B26F-0024E876A814.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/60DE5988-E0F0-DE11-86F7-0024E8769BA1.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/602C301E-E0F0-DE11-B8FF-0015178C674C.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/5EB80B2C-E0F0-DE11-9DC8-0024E8768412.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/5470C523-E0F0-DE11-B2B4-0015178C64BC.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/52CB5B2F-E0F0-DE11-9564-0024E87680F4.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/403D8262-E0F0-DE11-9DDC-00151796C1E8.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/368D1336-E0F0-DE11-B0AE-001D0967D5A8.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/3668B44B-E0F0-DE11-AB95-00151796D814.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/3286EB68-E0F0-DE11-ADE6-0015178C6B8C.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/28BD9294-E0F0-DE11-8961-0024E86E8D18.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/26F58471-E0F0-DE11-A64A-001D0967DFFD.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/2465C860-E0F0-DE11-97D2-001D0967C649.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/244A0A53-E0F0-DE11-B74C-0024E8767D79.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/22FFC084-E0F0-DE11-8358-0024E8767D79.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/22D2771E-E0F0-DE11-BBD1-0024E87699C0.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/1E832668-E0F0-DE11-8635-0024E8769BA1.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/16FE8F36-E0F0-DE11-AAB2-001D0967DFE4.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/1630A134-E0F0-DE11-8A34-001D0967DC0B.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/0E028D65-E0F0-DE11-9F8A-0024E876A814.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0008/0076986C-E0F0-DE11-B280-0024E86E8D59.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0007/F81487D0-EEEE-DE11-8E88-0024E8768299.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0007/C822227F-F9EE-DE11-94C3-0024E8767D52.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0007/90E667DD-ECEE-DE11-B460-0024E8768272.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0007/8E70550C-F3EE-DE11-B124-0024E8767CEA.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0007/82500F77-F5EE-DE11-B770-0024E8768C64.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0007/78D08CDC-ECEE-DE11-8182-0024E86E8D9A.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0007/1080EF5E-F7EE-DE11-81BE-0024E876A86F.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0007/00EBD30A-F3EE-DE11-BB54-001D0967D517.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/FC659F33-E2EE-DE11-85B9-0024E86E8D8D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/FC5E6C73-E6EE-DE11-9863-001D0967D4A4.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/FC4E836F-E8EE-DE11-A698-0024E8768390.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/FC1656C7-D7EE-DE11-8E19-0024E876A814.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/F8CFE0B7-D9EE-DE11-893F-0024E86E8DB4.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/F8433D57-D3EE-DE11-981B-0024E8769ADE.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/F68D49C0-D9EE-DE11-8DE6-0024E8768C71.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/EAC5876C-E6EE-DE11-9FA5-0024E87687CB.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/E8183AC9-DBEE-DE11-B5B2-001D0967512F.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/E0473C6A-E4EE-DE11-AE29-0024E86E8D8D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/DE022750-D3EE-DE11-93EA-0024E876A807.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/D4485CA6-D9EE-DE11-99F6-0024E87681F0.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/D259446D-E4EE-DE11-A486-0024E876A86F.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/CEDDB885-D7EE-DE11-BBEA-001D0967D580.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/CE8CADE5-DDEE-DE11-A57A-0024E8768D5B.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/CA5CE3BA-D9EE-DE11-A842-001D0968F337.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/C004D0D3-DBEE-DE11-99C2-001D0967C649.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/B8F2E915-E0EE-DE11-B5CD-001D096761A1.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/B2D9EC90-EAEE-DE11-B760-0024E876636C.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/B0FAF1C6-DBEE-DE11-B9A4-001D0967DAC1.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/B0B950B0-D9EE-DE11-B4CB-0024E876A83B.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/ACD451EE-DDEE-DE11-B732-0024E8768072.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/ACBB9A50-D3EE-DE11-98AC-0024E8768BE2.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/AA36E654-D3EE-DE11-9C2F-0024E8768099.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/A68684C1-DBEE-DE11-B5AF-0024E86E8CFE.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/A61A43D3-DBEE-DE11-8098-001D0967DEEF.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/A2B218B7-D9EE-DE11-B6D3-0024E8768265.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/A0C194AB-E8EE-DE11-B0F7-001D0967D03A.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/9C98FF07-E0EE-DE11-8DA7-001D0967D2DD.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/9C7DAB33-E2EE-DE11-AC04-001D0967DEEA.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/964F76D0-DBEE-DE11-8210-001D0967DA99.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/8CA470B5-D9EE-DE11-9E91-001D0967DAC1.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/8CA03223-D1EE-DE11-8230-001D0967D643.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/8C3748E2-DDEE-DE11-A92B-001D0967DAC1.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/887BA3E9-CCEE-DE11-8F0D-0024E87699E7.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/886BDAE3-DDEE-DE11-B9F9-001D0967E04D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/86B4989E-D7EE-DE11-8595-001D0967C1E9.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/864B0381-D5EE-DE11-AA54-0024E87683EB.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/80FA526F-E4EE-DE11-93F9-0024E8768072.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/7CA100C2-D7EE-DE11-BCC3-001D0967D643.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/7AE9F535-E2EE-DE11-9729-0024E876993E.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/7A47E48F-D7EE-DE11-B280-0024E8767D45.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/7862D177-E4EE-DE11-AD8F-001D0967CEBE.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/7052F9EB-DDEE-DE11-B511-0024E86E8D9A.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/6AE7D5C0-D9EE-DE11-902B-0024E8768C8B.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/6ADCB14B-D3EE-DE11-ABBF-0024E8767D38.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/6A1E3B30-E2EE-DE11-8012-001D0967DC3D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/660A1166-D5EE-DE11-8BDF-0024E876A86F.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/647463D6-D7EE-DE11-A117-0024E87683DE.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/62F237E9-DDEE-DE11-9603-001D0967D0A3.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/5C07E14A-D3EE-DE11-B854-0024E8767CEA.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/5A7710BB-D9EE-DE11-B8C8-0024E8768C3D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/54E83F5C-E4EE-DE11-9975-0024E87680E7.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/5001BEDB-DBEE-DE11-B91A-0024E8768258.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/4E7C52BC-D9EE-DE11-BF6D-0024E876A807.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/48D5FE31-E2EE-DE11-A35E-0024E87680B3.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/485A0F5B-D3EE-DE11-9553-0024E876803E.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/480B60BF-DBEE-DE11-8A4F-0024E8768446.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/44454A4D-E2EE-DE11-B920-0024E8768101.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/3CA6BE87-D5EE-DE11-B392-0024E87681F0.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/36DC526E-D5EE-DE11-B0AD-0024E86E8D59.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/321D3314-E0EE-DE11-A1CE-0024E8768867.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/30A0BF94-E8EE-DE11-9AD3-0024E87699A6.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/282ABC75-E8EE-DE11-86F5-0024E87680DA.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/26E6CC21-D1EE-DE11-B601-0024E86E8D66.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/2466002A-D1EE-DE11-AE33-0024E8768C64.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/24404114-E0EE-DE11-865C-0024E876A7D3.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/200D1815-E0EE-DE11-96D0-001D0967DC3D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/1A5972D3-DBEE-DE11-BBD3-001D0967DEEA.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/1A2A5C83-D5EE-DE11-8944-0024E8767D45.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/18C9B3D3-DBEE-DE11-A44D-0024E876A7E0.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/18257394-D7EE-DE11-80EE-0024E8769B05.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/14C0B4B4-D9EE-DE11-8F40-001D0967DEEF.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/1483D8D0-DBEE-DE11-9A7E-001D0967DC92.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/10F36BBF-D9EE-DE11-90B1-0024E87683B7.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/0ECD7750-D3EE-DE11-BD04-0024E876A848.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/0CECD1DF-EAEE-DE11-AD3A-001D0967DBE8.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_336p3_v1/0006/0C56077E-D5EE-DE11-88F4-0024E8768405.root'
#    '/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_341_v1/0004/001D84C8-90ED-DE11-B78E-002618943879.root'
#    '/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/E0D9D3FF-C6ED-DE11-910A-0015178C6A24.root'
#     '/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/C2AB8F41-E3ED-DE11-8F63-00151796D678.root'
#     '/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0102/7405A6D7-17ED-DE11-8C16-001D0967D4A4.root'
#    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/9E264E41-72EA-DE11-ACC1-0024E87680CD.root'
#     '/store/data/BeamCommissioning09/MinimumBias/RECO/Dec14thReReco_v1/0100/46E539F7-C4E9-DE11-921C-001D0967D03A.root'
#    '/store/data/BeamCommissioning09/MinimumBias/RECO/Dec14thReReco_v1/0100/EE7868BA-9CE9-DE11-B52D-0015178C6BF4.root'
    )
)

#process.GlobalTag.globaltag = 'GR09_P_V7::All' #rereco of 14th of Dec
#process.GlobalTag.globaltag = 'GR09_R_34X_V2::All' #rereco of 19th of Dec with CMSSW341
process.GlobalTag.globaltag = 'GR09_R_V5::All' #rereco pf 19th Dec with 336p3
#process.GlobalTag.globaltag = 'GR09_R_V4::All' #rereco pf 14th Dec with 335p4


process.MuonAnalysis = cms.EDAnalyzer("Onia2MuMu",
    OutputFileName       = cms.string('jpsi_data09_900GeV.root'),
    OniaType             = cms.int32(443),
    OniaMaxCat           = cms.int32(6),
    skimOnOniaMaxCat     = cms.bool(False),
    DebugLevel           = cms.int32(0),
    genParticlesLabel    = cms.InputTag("genParticles"),
    # StandAloneMuonsLabel = cms.InputTag("standAloneMuons"),
    GlobalMuonsLabel     = cms.InputTag("globalMuons"),
    MuonsLabel           = cms.InputTag("muons"),
    CaloMuonsLabel       = cms.InputTag("calomuons"),
    BeamSpotLabel        = cms.InputTag("offlineBeamSpot"),
    PrimaryVertexLabel   = cms.InputTag("offlinePrimaryVerticesWithBS"),
    UsePrimaryNoMuons    = cms.bool(True),
    TrackLabel           = cms.InputTag("generalTracks"),
    PhotonLabel          = cms.InputTag("particleFlow"),
    PhotonMinEnergy      = cms.double(2.0),
    triggerEventLabel    = cms.string("hltTriggerSummaryAOD"),
    triggerResultsLabel  = cms.string("TriggerResults"),
    HLTprocessName8e29   = cms.string("HLT"),
    HLTprocessName1e31   = cms.string("HLT"),  
    L1GTReadoutRec       = cms.InputTag("gtDigis"),
    L1MuonLabel          = cms.InputTag("l1extraParticles"),
    PATMuonsLabel        = cms.InputTag("selectedLayer1Muons"),
    # PATPhotonsLabel      = cms.InputTag("selectedLayer1Photons"),
    useKinFit            = cms.bool(False),
    StoreGenFlag         = cms.bool(False),
    StoreHLTFlag         = cms.bool(True),
    StoreL1Flag          = cms.bool(True),
    StoreTrkFlag         = cms.bool(True),
    StorePhotonFlag      = cms.bool(False),
    # StorePFMuonFlag      = cms.bool(True),
    StoreTRKMuonFlag     = cms.bool(True),
    StoreGLBMuonFlag     = cms.bool(True),
    StoreCALMuonFlag     = cms.bool(True),
    StoreBeamSpotFlag    = cms.bool(True),
    StorePriVtxFlag      = cms.bool(True),
    StoreOniaFlag        = cms.bool(True),
    StoreChicFlag        = cms.bool(False),
    StoreBpFlag          = cms.bool(False),                     
    StoreWSOnia          = cms.bool(True),                                  
    StoreOniaRadiation   = cms.bool(False),
    UsingBeamSpot        = cms.bool(False),
    minimumFlag          = cms.bool(False),
    UsingAOD             = cms.bool(False),
    StorePATFlag         = cms.bool(False)
)

# this is for filtering on L1 technical trigger bit
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
# different choices of trigger conditions:
# bsc minbias
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41)')
# bsc minbias and veto on beam halo
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
# bsc minbias in coincidence with bptx
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41)')
# bsc minbias in coinidence with bptx and veto on beam halo
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')


# this is for filtering on HLT physics bit
process.hltPhysicsDeclared = cms.EDFilter("HLTHighLevel",
                                 TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                 HLTPaths = cms.vstring("HLT_PhysicsDeclared"
                                                        ),
                                 eventSetupPathsKey = cms.string(''),
                                 andOr = cms.bool(True),
                                 throw = cms.bool(True)
                                 )


# Beamspot temporary fix
#from CondCore.DBCommon.CondDBSetup_cfi import *
#process.firstCollBeamspot = cms.ESSource("PoolDBESSource",
#                                         CondDBSetup,
#                                         connect = cms.string("frontier://PromptProd/CMS_COND_31X_BEAMSPOT"),
#                                         # connect = cms.string("sqlite_file:/afs/cern.ch/user/y/yumiceva/public/BeamSpotObjects_2009_v1_express@a35f2218-e25f-11de-9d9b-0018f34695d6.db"),
#                                         toGet = cms.VPSet(cms.PSet(record = cms.string("BeamSpotObjectsRcd"),
#                                                                    tag =cms.string("BeamSpotObjects_2009_v3_offline"))
#                                                           )
#                                         )
#process.es_prefer_firstCollBeamspot = cms.ESPrefer("PoolDBESSource","firstCollBeamspot")
#process.load('RecoVertex.BeamSpotProducer.BeamSpot_cfi')

# filter on lumisections
from HeavyFlavorAnalysis.Onia2MuMu.goodLumiSectionList_900GeV_cfi import *
process.source.lumisToProcess = goodLumisToProcess

#process.ana_step = cms.Path(process.hltLevel1GTSeed + process.hltPhysicsDeclared + process.offlineBeamSpot + process.MuonAnalysis)
#process.ana_step = cms.Path(process.hltLevel1GTSeed + process.hltPhysicsDeclared + process.MuonAnalysis)
process.ana_step = cms.Path(process.hltPhysicsDeclared + process.MuonAnalysis)
#process.ana_step = cms.EndPath(process.MuonAnalysis)

# Schedule definition
process.schedule = cms.Schedule(process.ana_step)

