#ifndef CommonTools_PFCandProducer_PFPileUpAlgo_
#define CommonTools_PFCandProducer_PFPileUpAlgo_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/SoftKiller.hh"
#include "fastjet/contrib/JetCleanser.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/Selector.hh"



typedef std::vector< edm::FwdPtr<reco::PFCandidate> >  PFCollection; // particle flow candidate collection
  
class PFPileUpAlgo { //PF Algorithm for pile-up mitigation 
 
 public:

  PFPileUpAlgo(): checkClosestZVertex_(true), verbose_(false)  {;} // basic constructor, no input pfCandidates and vertex given
    
  PFPileUpAlgo( bool checkClosestZVertex, bool verbose=false): // parsing of boolean informations
   checkClosestZVertex_(checkClosestZVertex), verbose_(verbose) {;}

  PFPileUpAlgo(const PFCollection & pfCandidates, const reco::VertexCollection & vertices, bool checkClosestZVertex, bool verbose=false); // give more input in the constructor

  ~PFPileUpAlgo();

  // set functions 
  inline void setPFCandidatesCollection(const PFCollection & pfCandidates) { pfCandidates_ = pfCandidates; }
  inline void setPFCandidatesFromPU(const PFCollection & pfCandidatesFromPU)  { pfCandidatesFromPU_  = pfCandidatesFromPU; }
  inline void setPFCandidatesFormLV(const PFCollection & pfCandidatesFromVtx) { pfCandidatesFromVtx_ = pfCandidatesFromVtx; }
  inline void setVerticesCollection(const reco::VertexCollection & vertices) { vertices_ = &vertices; }
  inline void setVerbose(bool verbose) { verbose_ = verbose; }
  inline void setCheckClosestZVertex(bool val) { checkClosestZVertex_ = val;}

  // get functions 
  const PFCollection & getPFCandidates()             const {return pfCandidates_;}  
  const PFCollection & getPFCandidatesFromPU()       const {return pfCandidatesFromPU_;}  
  const PFCollection & getPFCandidatesFromVtx()      const {return pfCandidatesFromVtx_;}
  const PFCollection & getSoftKillerPFCandidates()   const {return pfSoftKiller_;}
  const PFCollection & getPuppiPFCandidates()        const {return pfPuppi_;}
  const PFCollection & getJetCleansingPFCandidates() const {return pfJetCleansing_;}
  const PFCollection & getConstituentSubtractionPFCandidates()  const {return pfConstituentSubtraction_;}


  // implement charged hadron subtraction
  void processChargedHadronSubtraction(const PFCollection & pfCandidates, 
                                       const reco::VertexCollection & vertices); // charged hadron subtraction implementation 

  void processChargedHadronSubtraction(); // charged hadron subtraction implementation 


  // implement softKiller pfcandidates
  void processSoftKiller(const PFCollection & pfCandidates, 
                         const reco::VertexCollection & vertices,
                         const edm::ParameterSet& param);

  void processSoftKiller(const edm::ParameterSet& param);

  // implement puppi pfcandidates
  void processPuppi(const PFCollection & pfCandidates, 
                    const reco::VertexCollection & vertices,
                    const edm::ParameterSet& param)  ;

  void processPuppi(const edm::ParameterSet& param)  ;

  // implement pfcandidates after jet cleansing
  void processJetCleansing(const PFCollection & pfCandidates, 
                           const reco::VertexCollection & vertices,
                           const edm::ParameterSet& param);

  void processJetCleansing(const edm::ParameterSet& param);

  // implement pfcandidates after constituent subtraction
  void processConstituentSubtraction(const PFCollection & pfCandidates,
                                     const reco::VertexCollection & vertices,
                                     const edm::ParameterSet& param);

  void processConstituentSubtraction(const edm::ParameterSet& param);


  // identify charged hadron particels from PU vertex
  int chargedHadronVertex(const reco::VertexCollection& vertices, 
 			  const reco::PFCandidate& pfcand ) const;
 

  // convert string to fastjet::JetAlgorithm
  fastjet::JetAlgorithm get_algo(const std::string & algo);

  // convert a PFCollection into a vector of PseudoJet useful for fastjet  
  std::vector<fastjet::PseudoJet> ConvertToPseudoJet(const PFCollection & pfCandidates);

  // convert a a vector of PseudoJet where each pseudo jet is a particle into a PFCollection given the original set of  pfCandidates given before the clustering  
  PFCollection ConvertToPFCollection(const std::vector<fastjet::PseudoJet> & inputJets, const PFCollection & pfParticles);
  // convert a a vector of PseudoJet where each pseudo jet is a jet
  PFCollection ConvertJetToPFCollection(const std::vector<fastjet::PseudoJet> & inputJets, const PFCollection & pfParticles);
 
  // get constituents for cleaning --> take all the particles after clustering, divide them in neutrals, chargedLV and chargedPU using the PFCollections
  void getConstitsForCleansing(const std::vector<fastjet::PseudoJet> & jetParticles, 
                               std::vector<fastjet::PseudoJet> & neutrals,
			       std::vector<fastjet::PseudoJet> & chargedLV, 
                               std::vector<fastjet::PseudoJet> & chargedPU,
			       const PFCollection & eventParticles, const PFCollection & vertexParticles, const PFCollection & pileupParticles);
    
 private  :

  static constexpr double ymax     = 2.5;
  static constexpr double cell_size = 0.4;
  static constexpr double jetR     = 0.4;
  static constexpr double subjetR  = 0.3;
  static constexpr double jetRRho  = 0.4;
  static constexpr double GhostEtaMax  = 5.0;


  /// use the closest z vertex if a track is not in a vertex
  bool   checkClosestZVertex_;    
  /// verbose ?
  bool   verbose_;

  // private collections
  PFCollection pfCandidates_; // original set of pfCandidates
  const reco::VertexCollection* vertices_; // vertex collection

  PFCollection pfCandidatesFromVtx_; // pfCandidates from LV
  PFCollection pfCandidatesFromPU_;  // pfCandidates from PU
  PFCollection pfSoftKiller_ ;       // soft killer candidates
  PFCollection pfPuppi_ ;            // puppi candidates
  PFCollection pfJetCleansing_ ;     // candidates after jet cleansing
  PFCollection pfConstituentSubtraction_ ; // candidates after constituent subtraction

};

#endif
