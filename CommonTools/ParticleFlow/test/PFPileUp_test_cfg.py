import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python');
options.register ('isMC',False,VarParsing.multiplicity.singleton, VarParsing.varType.int,"option in order to run the analyzer on the MC")
options.parseArguments();
print options;

process = cms.Process( "PileUP" )

process.load("FWCore.MessageService.MessageLogger_cfi");
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

### load input file, options, global tag
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring (options.inputFiles));

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
process.path = cms.Path(process.pfPileUp)
process.EndPath = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.path,process.EndPath)

############################
## Dump the output Python ##
############################
processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
