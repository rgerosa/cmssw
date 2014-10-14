#ifndef CommonTools_PFCandProducer_PFPileUpAlgo_
#define CommonTools_PFCandProducer_PFPileUpAlgo_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "CommonTools/ParticleFlow/interface/PuppiAlgo.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/SoftKiller.hh"
#include "fastjet/contrib/JetCleanser.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/Selector.hh"


#include "TMath.h"
#include "Math/ProbFunc.h"


template <class T>  
class PFPileUpAlgo { //PF Algorithm for pile-up mitigation 
 
 public:

  typedef std::vector<edm::FwdPtr<T> >  CandCollection; // particle flow candidate collection
  typedef std::vector<edm::FwdPtr<reco::PFCandidate> > PFCollection; // particle flow candidate collection
  typedef std::vector<edm::FwdPtr<pat::PackedCandidate> >  PackedCollection; // particle flow candidate collection

  PFPileUpAlgo(): checkClosestZVertex_(true), verbose_(false)  {;} // basic constructor, no input pfCandidates and vertex given
    
  PFPileUpAlgo( bool checkClosestZVertex, bool verbose=false): // parsing of boolean informations
   checkClosestZVertex_(checkClosestZVertex), verbose_(verbose) {;}

  PFPileUpAlgo(const CandCollection & pfCandidates, const reco::VertexCollection & vertices, bool checkClosestZVertex, bool verbose = false); // give more input in the constructor

  ~PFPileUpAlgo();

  // set functions 
  inline void setPFCandidatesCollection(const CandCollection & pfCandidates)    { pfCandidates_ = pfCandidates; }
  inline void setPFCandidatesFromPU(const CandCollection & pfCandidatesFromPU)  { pfCandidatesFromPU_  = pfCandidatesFromPU; }
  inline void setPFCandidatesFormPV(const CandCollection & pfCandidatesFromVtx) { pfCandidatesFromVtx_ = pfCandidatesFromVtx; }
  inline void setPFCandidatesFormPVCharged(const CandCollection & pfCandidatesFromVtxCharged) { pfCandidatesFromVtxCharged_ = pfCandidatesFromVtxCharged; }
  inline void setVerticesCollection(const reco::VertexCollection & vertices)  { vertices_ = &vertices; }
  inline void setVerbose(const bool & verbose)         { verbose_ = verbose; }
  inline void setCheckClosestZVertex(const bool & val) { checkClosestZVertex_ = val;}

  // get functions 
  const CandCollection & getPFCandidates()             const {return pfCandidates_;}  
  const CandCollection & getPFCandidatesFromPU()       const {return pfCandidatesFromPU_;}  
  const CandCollection & getPFCandidatesFromVtx()      const {return pfCandidatesFromVtx_;}
  const CandCollection & getPFCandidatesFromVtxCharged()  const {return pfCandidatesFromVtxCharged_;}
  const CandCollection & getSoftKillerPFCandidates()   const {return pfSoftKiller_;}
  const CandCollection & getPuppiPFCandidates()        const {return pfPuppi_;}
  const CandCollection & getJetCleansingPFCandidates() const {return pfJetCleansing_;}
  const CandCollection & getConstituentSubtractionPFCandidates() const {return pfConstituentSubtraction_;}


  // implement charged hadron subtraction
  void processChargedHadronSubtraction(const CandCollection & pfCandidates, 
                                       const reco::VertexCollection & vertices,
                                       const double & DeltaZCut = 1000); // charged hadron subtraction implementation 

  void processChargedHadronSubtraction(); // charged hadron subtraction implementation 


  // implement softKiller pfcandidates
  void processSoftKiller(const CandCollection & pfCandidates, 
                         const reco::VertexCollection & vertices,
                         const edm::ParameterSet& param);

  void processSoftKiller(const edm::ParameterSet& param);

  // implement puppi pfcandidates
  void processPuppi(const CandCollection & pfCandidates, 
                    const reco::VertexCollection & vertices,
                    const edm::ParameterSet& param)  ;

  void processPuppi(const edm::ParameterSet& param)  ;

  // implement pfcandidates after jet cleansing
  void processJetCleansing(const CandCollection & pfCandidates, 
                           const reco::VertexCollection & vertices,
                           const edm::ParameterSet& param);

  void processJetCleansing(const edm::ParameterSet& param);

  // implement pfcandidates after constituent subtraction
  void processConstituentSubtraction(const CandCollection & pfCandidates,
                                     const reco::VertexCollection & vertices,
                                     const edm::ParameterSet& param);

  void processConstituentSubtraction(const edm::ParameterSet& param);


  // identify charged hadron particels from PU vertex
  int chargedHadronVertex(const reco::VertexCollection& vertices, 
 			  const reco::PFCandidate & pfcand ) const;

  // identify charged hadron particels from PU vertex
  int chargedHadronVertex(const reco::VertexCollection& vertices, 
 			  const pat::PackedCandidate & packedCand ) const;
 

  // convert string to fastjet::JetAlgorithm
  fastjet::JetAlgorithm get_algo(const std::string & algo);

  // convert a CandCollection into a vector of PseudoJet useful for fastjet  
  std::vector<fastjet::PseudoJet> ConvertToPseudoJet(const CandCollection & pfCandidates);

  // convert a a vector of PseudoJet where each pseudo jet is a particle into a CandCollection given the original set of  pfCandidates given before the clustering  
  CandCollection ConvertToCandCollection(const std::vector<fastjet::PseudoJet> & inputJets, const CandCollection & pfParticles);

  // convert a a vector of PseudoJet where each pseudo jet is a jet
  CandCollection ConvertJetToCandCollection(const std::vector<fastjet::PseudoJet> & inputJets, const CandCollection & pfParticles);
 
  // get constituents for cleaning --> take all the particles after clustering, divide them in neutrals, chargedLV and chargedPU using the CandCollections
  void getConstitsForCleansing(const std::vector<fastjet::PseudoJet> & jetParticles, 
                               std::vector<fastjet::PseudoJet> & neutrals,
			       std::vector<fastjet::PseudoJet> & chargedLV, 
                               std::vector<fastjet::PseudoJet> & chargedPU,
			       const CandCollection & eventParticles, const CandCollection & vertexParticles, const CandCollection & pileupParticles);

  // Puppi metric values and RMS and mean value for PU particles
  void getPuppiRMSAvg(const int & iOpt, std::vector<fastjet::PseudoJet> & eventParticles, std::vector<fastjet::PseudoJet> & chargedParticlesPV);
  // get puppi algorithm id as a function of eta and pt
  int getPuppiId(const float & pt,const float & eta);
  // get the chi2 value for the dZ distance of the particles track wrt to the leading vertex
  double getChi2FromdZ(const double & dZ);
  // get the metric value
  double goodVarPuppi(const fastjet::PseudoJet & iPart, std::vector<fastjet::PseudoJet> & iParts, int iOpt, double iRCone);   

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
  CandCollection pfCandidates_;               // original set of pfCandidates
  CandCollection pfCandidatesFromVtx_;        // pfCandidates from LV
  CandCollection pfCandidatesFromVtxCharged_; // pfCandidates from LV charged
  CandCollection pfCandidatesFromPU_;         // pfCandidates from PU
  const reco::VertexCollection* vertices_;  // vertex collection

  CandCollection pfSoftKiller_ ;              // soft killer candidates
  CandCollection pfJetCleansing_ ;            // candidates after jet cleansing
  CandCollection pfConstituentSubtraction_ ;  // candidates after constituent subtraction

  CandCollection pfPuppi_ ;            // puppi candidates
  std::vector<int> particleIndex_;   // index for particles from PV,PU .. etc
  std::vector<PuppiAlgo> fPuppiAlgo_;
  std::vector<double> puppiWeights_ ;
  std::vector<edm::ParameterSet> lAlgos_;  
  std::vector<double> puppiMetricValues_;
};


#endif
