import FWCore.ParameterSet.Config as cms


pfPileUp = cms.EDProducer("PFPileUp",
    PFCandidates = cms.InputTag("particleFlowPtrs"), ## input collection
    Vertices     = cms.InputTag("offlinePrimaryVertices"),  ## vertexes 
    Enable  = cms.bool(True), ## enable the module to run
    ##CHS section
    produceCHS = cms.bool(True), ## produce CHS collections
    checkClosestZVertex = cms.bool(True), ## select vertex z close
    ##soft killer section
    produceSoftKiller = cms.bool(True),
    softKillerParam = cms.PSet( ymax      = cms.double(2.5), 
                                cell_size = cms.double(0.4),
                                applyCHS  = cms.bool(False),
                                clusterInJets = cms.bool(False),
                                jetAlgo   = cms.string("antikt_algorithm"),
                                jetR      = cms.double(0.8), 
                                GhostEtaMax = cms.double(5.0),
                                jetPtCut = cms.double(10.), 
    ),
    ##puppi section
    producePuppi      = cms.bool(False),
    puppiParam = cms.PSet(),
    ##jet cleansing section
    produceJetCleansing = cms.bool(True),
    jetCleansingParam = cms.PSet(   jetAlgo   = cms.string("antikt_algorithm"),
                                    jetR      = cms.double(0.8), 
                                    subjetAlgo = cms.string("kt_algorithm"),
                                    subjetR   = cms.double(0.3), 
                                    GhostEtaMax = cms.double(5.0),
                                    cleansingMethod = cms.string("Gaussian"),
                                    cleansingParam = cms.vdouble(0.617,0.62,0.15,0.22),
                                    applyCHS  = cms.bool(False),
                                    jetPtCut = cms.double(10.), 
    ),
    ##constituent subtraction section
    produceConstituentSubtraction = cms.bool(True),
    constituentSubtractionParam  = cms.PSet(
                                    jetAlgoRho   = cms.string("kt_algorithm"),
                                    jetRRho      = cms.double(0.4),                                    
                                    applyCHS    = cms.bool(False),
                                    clusterInJets = cms.bool(True),
                                    jetAlgo   = cms.string("antikt_algorithm"),
                                    jetR      = cms.double(0.8),                                    
                                    jetPtCut    = cms.double(0.), 
                                    GhostEtaMax = cms.double(5.0),
   ),

   verbose = cms.untracked.bool(False), ## verbosity
)
