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
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"


class PFPileUp : public edm::stream::EDProducer<> {

 public:

  typedef std::vector<edm::FwdPtr<reco::PFCandidate> >  PFCollection; // pfCandidates vector inputs
  typedef edm::View<reco::PFCandidate>                  PFView; // view from event
  typedef std::vector<reco::PFCandidate>                PFCollectionByValue;

  explicit PFPileUp(const edm::ParameterSet&); // default constructor

  ~PFPileUp();

  virtual void produce(edm::Event&, const edm::EventSetup&) override; // produce methdo

 private:

  PFPileUpAlgo    pileUpAlgo_; // object to be called to apply the algorithm

  /// PFCandidates to be analyzed
  edm::EDGetTokenT<PFCollection>   tokenPFCandidates_;
  edm::EDGetTokenT<PFView>         tokenPFCandidatesView_;
  /// vertices
  edm::EDGetTokenT<reco::VertexCollection>   tokenVertices_;

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
  bool produceConstituentSubtraction_ ;


};

#endif
