#include "CommonTools/ParticleFlow/plugins/PFMET.h"
#include "CommonTools/ParticleFlow/plugins/Type1PFMET.h"
#include "CommonTools/ParticleFlow/plugins/PFPileUp.h"
#include "CommonTools/ParticleFlow/plugins/PFPileUp.cc"
#include "CommonTools/ParticleFlow/plugins/PFCandidateFwdPtrCollectionFilter.h"
#include "CommonTools/ParticleFlow/plugins/PFJetFwdPtrProducer.h"
#include "CommonTools/ParticleFlow/plugins/PFTauFwdPtrProducer.h"
#include "CommonTools/ParticleFlow/plugins/PFCandidateFromFwdPtrProducer.h"
#include "CommonTools/ParticleFlow/plugins/DeltaBetaWeights.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

DEFINE_FWK_MODULE(PFMET);
DEFINE_FWK_MODULE(Type1PFMET);
typedef PFPileUp<reco::PFCandidate> PFPileUpPFCandidate;
DEFINE_FWK_MODULE(PFPileUpPFCandidate);
typedef PFPileUp<pat::PackedCandidate> PFPileUpPackedCandidate;
DEFINE_FWK_MODULE(PFPileUpPackedCandidate);


DEFINE_FWK_MODULE(PFCandidateFwdPtrCollectionStringFilter);
DEFINE_FWK_MODULE(PFCandidateFwdPtrCollectionPdgIdFilter);
DEFINE_FWK_MODULE(PFJetFwdPtrProducer);
DEFINE_FWK_MODULE(PFTauFwdPtrProducer);
DEFINE_FWK_MODULE(PFCandidateFromFwdPtrProducer);

typedef edm::ProductFromFwdPtrProducer< reco::PFJet >  PFJetFromFwdPtrProducer;
DEFINE_FWK_MODULE(PFJetFromFwdPtrProducer);

DEFINE_FWK_MODULE(DeltaBetaWeights);
