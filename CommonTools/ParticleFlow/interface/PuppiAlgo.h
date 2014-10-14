#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "fastjet/PseudoJet.hh"
#include <vector>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//// Puppi Alogorithm class
class PuppiAlgo { 
 public:
  PuppiAlgo(edm::ParameterSet &iConfig); // constructor from a config parameter set
  ~PuppiAlgo();

  
  void   reset(); // collection reset  
  void   add(const fastjet::PseudoJet & iParticle,  const edm::FwdPtr<reco::PFCandidate> & pfCandidate, const bool & isFromPV, const double & iVal); // add a particle (PseudoJet), an algoID and a Value
  void   add(const fastjet::PseudoJet & iParticle,  const edm::FwdPtr<pat::PackedCandidate> & pfCandidate, const bool & isFromPV, const double & iVal); // add a particle (PseudoJet), an algoID and a Value
  void   computeMedRMS(const double & iPVFrac); // compute RMS given the algo id and  a PVfraction
  double compute(double & iVals, double iChi2); // compute the final weights given a dZ chi2

  double getPtMin(); // get Min pt set in the procedure
  double getEtaMin(); // get min eta
  double getEtaMax(); // get max eta
  int    getAlgoId    (); // get the algo Id
  bool   getIsCharged (); // get if is charged or not
  double getConeSize  (); // get the cone size
  double getNeutralPt (int iNPV); // get neutral pt cut

private:  

  float  fEtaMax_; // eta maximum
  float  fEtaMin_; // eta minimum
  float  fPtMin_ ;  // pt minimum
  double fNeutralPtMin_; // minimum neutral pt
  double fNeutralPtSlope_; // minimum neutral pt slope vs Nvtx
  std::vector<float>  fPupsPU_; 
  std::vector<float>  fPupsPV_;
  int fAlgoId_;  // algo id
  bool fCharged_; // if charge hadron subtraction is applied
  bool fAdjust_;  
  int  fCombId_;
  double fConeSize_;
  double fRMSPtMin_;
  double fRMSScaleFactor_;
  double fRMS_; 
  double fMedian_;
  double fMean_;
  int    fNCount_;
};
