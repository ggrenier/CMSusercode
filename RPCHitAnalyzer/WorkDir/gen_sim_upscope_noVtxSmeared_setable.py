# Auto generated configuration file
# using: 
# Revision: 1.13 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleMuPt100_cfi --conditions POSTLS161_V12::All -s GEN,SIM --eventcontent FEVTDEBUG --datatier GEN-SIM -n 25 --geometry ExtendedPostLS1 --fileout SingleMuPt40_RE4_cfi_GEN_SIM.root --no_exec
import FWCore.ParameterSet.Config as cms
process = cms.Process('SIM')

#parameters
gg_a=dict()
gg_a['nevents']=10000
gg_a['Pt']=100
gg_a['Zvtx']=0 #mm
#gg_a['ZvtxSigma']=0.000001
gg_a['ZvtxSigma']=60 #mm
gg_a['etamin']=1.3
gg_a['etamax']=2.6

import getFileName
outputFileName=getFileName.generateFileName(gg_a)


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedPostLS1Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedPostLS1_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_upscope_cfi")

process.VtxSmeared.Phi        = cms.double(0.0)
process.VtxSmeared.BetaStar   = cms.double(70.0)
process.VtxSmeared.Emittance  = cms.double(0.586e-07)
process.VtxSmeared.Alpha      = cms.double(0.0)
process.VtxSmeared.SigmaZ     = cms.double(gg_a['ZvtxSigma']*0.1)
process.VtxSmeared.TimeOffset = cms.double(0.0)
process.VtxSmeared.X0         = cms.double(0.0)
process.VtxSmeared.Y0         = cms.double(0.0)
process.VtxSmeared.Z0         = cms.double(gg_a['Zvtx']*0.1)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(gg_a['nevents'])
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.13 $'),
    annotation = cms.untracked.string('SingleMuPt100_cfi nevts:10000'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string(outputFileName),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PRE_PO61_V1::All', '')

process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(gg_a['Pt']+0.01),
        MinPt = cms.double(gg_a['Pt']-0.01),
        PartID = cms.vint32(-13),
        MaxEta = cms.double(gg_a['etamax']),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(gg_a['etamin']),
        MinPhi = cms.double(-3.14159265359)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single mu pt 100'),
    AddAntiParticle = cms.bool(True),
    firstRun = cms.untracked.uint32(1)
)

# select generated muons and antimuons
process.genMuons = cms.EDFilter("PdgIdCandViewSelector",
    src = cms.InputTag("genParticles"),
    pdgId = cms.vint32( 13, -13 )
)

# Path and EndPath definitions
process.gen_mu_select = cms.Sequence(process.genMuons)
process.generation_step = cms.Path(process.pgen * process.gen_mu_select)
# process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

