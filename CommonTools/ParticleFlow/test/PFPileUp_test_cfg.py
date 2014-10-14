import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python');
options.register ('isMC',False,VarParsing.multiplicity.singleton, VarParsing.varType.int,"option in order to run the analyzer on the MC")
options.register ('isMiniAOD',False,VarParsing.multiplicity.singleton, VarParsing.varType.int,"option in order to run on MiniAOD instead of AOD sim")
options.parseArguments();
print options;

process = cms.Process( "PileUP" )

process.load("FWCore.MessageService.MessageLogger_cfi");
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

### load input file, options, global tag
process.source = cms.Source ("PoolSource");

if options.isMiniAOD :
 process.source.fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU40bx25_POSTLS170_V7-v2/00000/00800BE3-E826-E411-AD01-20CF3019DEE9.root') ;
else:
 process.source.fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU40bx25_POSTLS170_V5-v1/00000/000A6D7D-EB11-E411-9EEC-002590DB9152.root') ;


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents));

if options.isMC == 0:
 process.GlobalTag.globaltag = 'GR_R_71_V4::All';
else :
 process.GlobalTag.globaltag = 'DESIGN71_V5::All';

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.output = cms.OutputModule("PoolOutputModule",
   fileName = cms.untracked.string("output.root"),
   outputCommands = cms.untracked.vstring('drop *',
     'keep *_*_*_PileUP',
     'keep *_*particleFlowPtrs*_*_*',                                       
    ),
    dropMetaData = cms.untracked.string('ALL')
)

process.load('CommonTools.ParticleFlow.pfPileUp_cfi')
if options.isMiniAOD :
 process.path = cms.Path(process.PFPileUpPackedCandidate)
else:
 process.path = cms.Path(process.PFPileUpPFCandidate)

process.EndPath = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.path,process.EndPath)

############################
## Dump the output Python ##
############################
processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
