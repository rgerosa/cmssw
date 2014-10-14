#ifndef PhysicsTools_PFCandProducer_PFPileUp_
#define PhysicsTools_PFCandProducer_PFPileUp_

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h"

#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"
#include "CommonTools/ParticleFlow/src/PFPileUpAlgo.cc"

template <class T>
class PFPileUp : public edm::EDProducer {

 public:

  typedef std::vector<edm::FwdPtr<T> >  PFCollection;  // pfCandidates vector inputs --> PFCandidate or Packed ones from SealModule
  typedef edm::View<T>    PFView;        // View from event

  explicit PFPileUp(const edm::ParameterSet&); // default constructor

  ~PFPileUp();

  virtual void produce(edm::Event&, const edm::EventSetup&) override; // produce methdo

 private:

  PFPileUpAlgo<T>    pileUpAlgo_; // object to be called to apply the algorithm

  /// PFCandidates to be analyzed
  edm::EDGetTokenT<PFCollection>           tokenPFCandidates_;
  edm::EDGetTokenT<PFView>                 tokenPFCandidatesView_;
  edm::EDGetTokenT<reco::VertexCollection> tokenVertices_;

  /// enable PFPileUp selection
  bool   enable_;

  /// verbose ?
  bool   verbose_;

  /// use the closest z vertex if a track is not in a vertex
  bool   checkClosestZVertex_;

  edm::ParameterSet  softKillerParam_ ;
  edm::ParameterSet  puppiParam_ ;
  edm::ParameterSet  jetCleansingParam_ ;
  edm::ParameterSet  constituentSubtractionParam_ ;

  bool produceCHS_ ;
  bool produceSoftKiller_ ;
  bool producePuppi_ ;
  bool produceJetCleansing_ ;
  bool produceConstituentSubtraction_ ;}
;

#endif
