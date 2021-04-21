import FWCore.ParameterSet.Config as cms

pfDeepBoostedJetTagInfos = cms.EDProducer('DeepBoostedJetTagInfoProducer',
  ## jet properties
  jet_radius = cms.double(0.8),
  min_jet_pt = cms.double(150),
  max_jet_eta = cms.double(99),
  ## threshold used for the training 
  include_neutrals = cms.bool(True),
  min_pt_for_track_properties = cms.double(-1),
  min_pt_for_pfcandidates = cms.double(-1),
  ## puppi 
  use_puppiP4 = cms.bool(True),
  min_puppi_wgt = cms.double(0.01),
  ## sorting 
  sort_by_sip2dsig = cms.bool(False),
  flip_ip_sign = cms.bool(False),
  sip3dSigMax = cms.double(-1),
  ## HLT vs reco
  use_hlt_features = cms.bool(False),
  ## input collections                                          
  vertices           = cms.InputTag('offlinePrimaryVertices'),
  secondary_vertices = cms.InputTag('inclusiveCandidateSecondaryVertices'),
  pf_candidates      = cms.InputTag('particleFlow'),
  jets               = cms.InputTag('ak8PFJetsPuppi'),
  puppi_value_map    = cms.InputTag('puppi'),
  vertex_associator  = cms.InputTag('primaryVertexAssociation', 'original'),
  mightGet = cms.optional.untracked.vstring
)

pfDeepBoostedJetTagInfosHLT = cms.EDProducer('DeepBoostedJetTagInfoProducer',
  ## jet properties
  jet_radius = cms.double(0.4),
  min_jet_pt = cms.double(25),
  max_jet_eta = cms.double(2.5),
  ## threshold used for the training 
  include_neutrals = cms.bool(True),
  min_pt_for_track_properties = cms.double(0.95),
  min_pt_for_pfcandidates = cms.double(0.1),
  ## puppi 
  use_puppiP4 = cms.bool(False),
  min_puppi_wgt = cms.double(0.01),
  ## sorting 
  sort_by_sip2dsig = cms.bool(False),
  flip_ip_sign = cms.bool(False),
  sip3dSigMax = cms.double(-1),
  ## HLT vs reco
  use_hlt_features = cms.bool(True),
  ## input collections                                          
  vertices           = cms.InputTag('hltVerticesPFFilter'),
  secondary_vertices = cms.InputTag('hltDeepInclusiveMergedVerticesPF'),
  pf_candidates      = cms.InputTag('hltParticleFlow'),
  jets               = cms.InputTag('hltPFJetForBtag'),
  puppi_value_map    = cms.InputTag(''),
  vertex_associator  = cms.InputTag('hltPrimaryVertexAssociation', 'original'),
  mightGet = cms.optional.untracked.vstring
)
