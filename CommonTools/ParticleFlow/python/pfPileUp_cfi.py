import FWCore.ParameterSet.Config as cms


PFPileUpPFCandidate = cms.EDProducer("PFPileUpPFCandidate",
    PFCandidates = cms.InputTag("particleFlowPtrs"),       ## input particle flow collection collection
    Vertices     = cms.InputTag("offlinePrimaryVertices"), ## vertexes 
    Enable  = cms.bool(True),                              ## enable the module to run

    ##CHS section
    produceCHS = cms.bool(True),          ## produce CHS collections
    checkClosestZVertex = cms.bool(True), ## select vertex z close

    ##soft killer section
    produceSoftKiller = cms.bool(True),
    softKillerParam = cms.PSet( ymax      = cms.double(5.0), ## maximum extension
                                cell_size = cms.double(0.4), ## patch dimension
                                applyCHS  = cms.bool(True), ## apply CHS --> softkiller only on CHS particles 
                                clusterInJets = cms.bool(False),  ## cluster soft killer products
                                jetAlgo   = cms.string("antikt_algorithm"), ## jet algo for clustering
                                jetR      = cms.double(0.8),  ## jetR dimension
                                GhostEtaMax = cms.double(5.0), ## Eta max for ghosts
                                jetPtCut = cms.double(10.),    ## pt cut on clustered jets
    ),
    ##puppi section
    producePuppi      = cms.bool(True),
    puppiParam = cms.PSet(  applyCHS = cms.untracked.bool (True),  ## give weight zero to PU charged hadrons and 1 to leading vertex
                            useExp   = cms.untracked.bool (False), ## multiply the puppi weights for a probablity of belonging to the vertex proportional to the dZ
                            MinPuppiWeight = cms.untracked.double(0.01), ## cut away particles with weight less than 1%
                            UseDeltaZCut = cms.untracked.bool (True),    ## pfCandidates with dZ > DeltaZCut from the leading vertex are put as pile-up
                            DeltaZCut = cms.untracked.double(0.2),
                            algos = cms.VPSet(),
    ), 

    ##jet cleansing section
    produceJetCleansing = cms.bool(True),
    jetCleansingParam = cms.PSet(   jetAlgo   = cms.string("antikt_algorithm"), ## clustering for cleansing
                                    jetR      = cms.double(0.8),                ## R dimension
                                    subjetAlgo = cms.string("kt_algorithm"),    ## subjet algo
                                    subjetR   = cms.double(0.3),                ## sub jet dimension
                                    GhostEtaMax = cms.double(5.0),
                                    cleansingMethod = cms.string("Gaussian"),   ## type of cleansing: JVF, Linear, Gaussian
                                    cleansingParam = cms.vdouble(0.617,0.62,0.15,0.22), ## parameters
                                    applyCHS  = cms.bool(False), ## apply CHS
                                    jetPtCut = cms.double(10.),  ## jet pt cut
    ),
    ##constituent subtraction section
    produceConstituentSubtraction = cms.bool(True),
    constituentSubtractionParam  = cms.PSet(
                                    jetAlgoRho   = cms.string("kt_algorithm"), ## jet clustering for constituent rho
                                    jetRRho      = cms.double(0.4),            ## dimension                         
                                    applyCHS    = cms.bool(True),             ## apply CHS
                                    clusterInJets = cms.bool(True),            ## cluster in hets
                                    jetAlgo   = cms.string("antikt_algorithm"),
                                    jetR      = cms.double(0.8),                                    
                                    jetPtCut    = cms.double(10.), 
                                    GhostEtaMax = cms.double(5.0),
   ),

   verbose = cms.untracked.bool(False), ## verbosity
)


## puppi algo setup
puppiCentral =  cms.PSet( 
            algoId = cms.untracked.int32(5), #0 is default Puppi
            useCharged = cms.untracked.bool(True),
            applyLowPUCorr = cms.untracked.bool(True),
            combOpt = cms.untracked.int32(0),
            cone = cms.untracked.double(0.3),
            rmsPtMin = cms.untracked.double(0.1),
            rmsScaleFactor = cms.untracked.double(1.0)
)


puppiForward = cms.PSet(
           algoId = cms.untracked.int32(5), #0 is default Puppi
           useCharged = cms.untracked.bool(False),
           applyLowPUCorr = cms.untracked.bool(True),
           combOpt = cms.untracked.int32(0),
           cone = cms.untracked.double(0.3),
           rmsPtMin = cms.untracked.double(0.5),
           rmsScaleFactor = cms.untracked.double(1.0)
)

PFPileUpPFCandidate.puppiParam.algos.insert(False,
 cms.PSet( etaMin = cms.untracked.double(3.0),
           etaMax = cms.untracked.double(10.0),
           ptMin = cms.untracked.double(0.0),
           MinNeutralPt = cms.untracked.double(1.5),
           MinNeutralPtSlope = cms.untracked.double(0.005),
           puppiAlgo = puppiForward)
)

PFPileUpPFCandidate.puppiParam.algos.insert(False,
       cms.PSet(etaMin = cms.untracked.double(2.5),
          etaMax = cms.untracked.double(3.0),
          ptMin = cms.untracked.double(0.0),
          MinNeutralPt = cms.untracked.double(1.0),
          MinNeutralPtSlope = cms.untracked.double(0.005),
          puppiAlgo = puppiForward),
)

PFPileUpPFCandidate.puppiParam.algos.insert(False,
 ## central region inside tracker acceptance
 cms.PSet( etaMin = cms.untracked.double(0.),
           etaMax = cms.untracked.double(2.5),
           ptMin = cms.untracked.double(0.), ## all pfCandidates are taken into account
           MinNeutralPt = cms.untracked.double(0.2), ## only neutrals with pT > 0.2 GeV
           MinNeutralPtSlope = cms.untracked.double(0.02),
           puppiAlgo = puppiCentral),
)

## on MiniAOD
PFPileUpPackedCandidate = cms.EDProducer("PFPileUpPackedCandidate",
    PFCandidates = cms.InputTag("packedPFCandidates"),       ## input particle flow collection collection
    Vertices     = cms.InputTag("offlineSlimmedPrimaryVertices"), ## vertexes 
    Enable  = cms.bool(True),                              ## enable the module to run

    ##CHS section
    produceCHS = cms.bool(True),          ## produce CHS collections
    checkClosestZVertex = cms.bool(True), ## select vertex z close

    ##soft killer section
    produceSoftKiller = cms.bool(True),
    softKillerParam = cms.PSet( ymax      = cms.double(5.0), ## maximum extension
                                cell_size = cms.double(0.4), ## patch dimension
                                applyCHS  = cms.bool(True), ## apply CHS --> softkiller only on CHS particles 
                                clusterInJets = cms.bool(False),  ## cluster soft killer products
                                jetAlgo   = cms.string("antikt_algorithm"), ## jet algo for clustering
                                jetR      = cms.double(0.8),  ## jetR dimension
                                GhostEtaMax = cms.double(5.0), ## Eta max for ghosts
                                jetPtCut = cms.double(10.),    ## pt cut on clustered jets
    ),
    ##puppi section
    producePuppi      = cms.bool(True),
    puppiParam = cms.PSet(  applyCHS = cms.untracked.bool (True),  ## give weight zero to PU charged hadrons and 1 to leading vertex
                            useExp   = cms.untracked.bool (False), ## multiply the puppi weights for a probablity of belonging to the vertex proportional to the dZ
                            MinPuppiWeight = cms.untracked.double(0.01), ## cut away particles with weight less than 1%
                            UseDeltaZCut = cms.untracked.bool (True),    ## pfCandidates with dZ > DeltaZCut from the leading vertex are put as pile-up
                            DeltaZCut = cms.untracked.double(0.2),
                            algos = cms.VPSet(),
    ), 

    ##jet cleansing section
    produceJetCleansing = cms.bool(True),
    jetCleansingParam = cms.PSet(   jetAlgo   = cms.string("antikt_algorithm"), ## clustering for cleansing
                                    jetR      = cms.double(0.8),                ## R dimension
                                    subjetAlgo = cms.string("kt_algorithm"),    ## subjet algo
                                    subjetR   = cms.double(0.3),                ## sub jet dimension
                                    GhostEtaMax = cms.double(5.0),
                                    cleansingMethod = cms.string("Gaussian"),   ## type of cleansing: JVF, Linear, Gaussian
                                    cleansingParam = cms.vdouble(0.617,0.62,0.15,0.22), ## parameters
                                    applyCHS  = cms.bool(False), ## apply CHS
                                    jetPtCut = cms.double(10.),  ## jet pt cut
    ),
    ##constituent subtraction section
    produceConstituentSubtraction = cms.bool(True),
    constituentSubtractionParam  = cms.PSet(
                                    jetAlgoRho   = cms.string("kt_algorithm"), ## jet clustering for constituent rho
                                    jetRRho      = cms.double(0.4),            ## dimension                         
                                    applyCHS    = cms.bool(True),             ## apply CHS
                                    clusterInJets = cms.bool(True),            ## cluster in hets
                                    jetAlgo   = cms.string("antikt_algorithm"),
                                    jetR      = cms.double(0.8),                                    
                                    jetPtCut    = cms.double(10.), 
                                    GhostEtaMax = cms.double(5.0),
   ),

   verbose = cms.untracked.bool(False), ## verbosity
)
 
