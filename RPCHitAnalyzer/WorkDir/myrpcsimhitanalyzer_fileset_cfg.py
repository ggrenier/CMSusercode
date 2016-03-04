import FWCore.ParameterSet.Config as cms

#inputFile should be provided elsewhere
#define inputFile in extraInfo.py 
from extraInfo import *

hereInputFile='../Data/'+inputFile
hereOutputFile='../Data/histos/histo_simhit_'+inputFile


process = cms.Process("testSimHitAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

#You need these 3 files to get the geometry
#process.load("Geometry.MuonCommonData.muonIdealGeometryXML_upscope_cfi")
#process.load('Configuration.Geometry.GeometryExtended2023Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')

process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'DES23_62_V1::All', '')



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:'+hereInputFile
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(hereOutputFile)
)


process.testee = cms.EDAnalyzer('MyRPCSimHitAnalyzer',
                                RPCSimHitCollection=cms.untracked.InputTag("g4SimHits","MuonRPCHits","SIM")
)


process.p = cms.Path(process.testee)
