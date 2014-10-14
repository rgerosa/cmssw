#include "../interface/PuppiAlgo.h"
#include "fastjet/internal/base.hh"
#include "Math/QuantFuncMathCore.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "TMath.h"

PuppiAlgo::PuppiAlgo(edm::ParameterSet &iConfig) { 

  fEtaMin_             = iConfig.getUntrackedParameter<double>("etaMin");
  fEtaMax_             = iConfig.getUntrackedParameter<double>("etaMax");
  fPtMin_              = iConfig.getUntrackedParameter<double>("ptMin");
  fNeutralPtMin_       = iConfig.getUntrackedParameter<double>("MinNeutralPt");      // Weighted Neutral Pt Cut
  fNeutralPtSlope_     = iConfig.getUntrackedParameter<double>("MinNeutralPtSlope"); // Slope vs #pv


  edm::ParameterSet lAlgos = iConfig.getParameter<edm::ParameterSet >("puppiAlgo");  // set of parameter set for the different algos to be run
  //User Configurable Puppi  --> import information for each algo
  fAlgoId_        = lAlgos.getUntrackedParameter<int >  ("algoId"); //algoId
  fCharged_       = lAlgos.getUntrackedParameter<bool>  ("useCharged"); //useCharge hadron subtraction
  fAdjust_        = lAlgos.getUntrackedParameter<bool>  ("applyLowPUCorr"); // apply Low PU correction
  fCombId_        = lAlgos.getUntrackedParameter<int>   ("combOpt"); // 0=> add in chi2/1=>Multiply p-values 
  fConeSize_      = lAlgos.getUntrackedParameter<double>("cone"); // Min Pt when computing pt and rms
  fRMSPtMin_      = lAlgos.getUntrackedParameter<double>("rmsPtMin"); // Min Pt when computing pt and rm
  fRMSScaleFactor_ = lAlgos.getUntrackedParameter<double>("rmsScaleFactor"); // Additional Tuning parameter for Jokers
  fRMS_    = 0.;
  fMedian_ = 0.;
  fMean_   = 0.;
  fNCount_ = 0;
  
}

PuppiAlgo::~PuppiAlgo() { 
  fPupsPU_.clear();
  fPupsPV_.clear();
  fAlgoId_  = 0;
  fCharged_ = 0;
  fAdjust_  = 0;
  fCombId_  = 0;
  fConeSize_ = 0;
  fRMSPtMin_ = 0;
  fRMSScaleFactor_= 0;
  fRMS_   = 0;
  fMedian_= 0;
  fMean_  = 0;
  fNCount_= 0;
  
}

// reset not inner parameters but inner values
void PuppiAlgo::reset() { 
  fPupsPU_.clear();
  fPupsPV_.clear();
  fMedian_ =  0; 
  fRMS_    =  0;
  fMean_   =  0;
  fNCount_ =  0;
}

// add a particle
void PuppiAlgo::add(const fastjet::PseudoJet & iParticle, const edm::FwdPtr<reco::PFCandidate> & pfCandidate, const bool & isFromPV, const double & iVal) { 
  if(iParticle.pt() < fRMSPtMin_) return; // if is under the pt cut --> skip it
  if(fCharged_ && pfCandidate->charge() == 0) return; // if you want only charge particles 
  if(fCharged_ && isFromPV) fPupsPV_.push_back(iVal);
  if(fCharged_ && not isFromPV && pfCandidate->charge() <= -2) return;
  fPupsPU_.push_back(iVal); // put the value in the vector .. one vector for each algo, don't want to do a matrix
  fNCount_++; // one more particles in the algo
}

// add a particle
void PuppiAlgo::add(const fastjet::PseudoJet & iParticle, const edm::FwdPtr<pat::PackedCandidate> & packedCandidate, const bool & isFromPV, const double & iVal) { 
  if(iParticle.pt() < fRMSPtMin_) return; // if is under the pt cut --> skip it
  if(fCharged_ && packedCandidate->charge() == 0) return; // if you want only charge particles 
  if(fCharged_ && isFromPV) fPupsPV_.push_back(iVal);
  if(fCharged_ && not isFromPV && packedCandidate->charge() <= -2) return;
  fPupsPU_.push_back(iVal); // put the value in the vector .. one vector for each algo, don't want to do a matrix
  fNCount_++; // one more particles in the algo
}

// compute the median RMS
void PuppiAlgo::computeMedRMS(const double &iPVFrac) { 

  if(fNCount_ == 0) return; // if no particles return

  std::sort(fPupsPU_.begin(),fPupsPU_.begin()+fNCount_); //sort in pt only the part of the vector related to that algo for PU particles

  double lCorr = 1.;
  if(fAdjust_) lCorr *= 1. - iPVFrac; // to adjust algo take the factor as 1-iPVfrac

  int lNum0 = 0;
  for(int i0 = 0; i0 < fNCount_; i0++) { 
    if(fPupsPU_[i0] == 0) lNum0 = i0; 
  }  // drop zero value particles and after sorting we can take the median

  int lNHalfway = lNum0 + int( double( fNCount_-lNum0 )*0.50*lCorr); // adjust the median position by charged fraction 
  fMedian_    = fPupsPU_[lNHalfway]; // take the media for that algo
  double lMed = fMedian_;  //Just to make the readability easier

  int lNRMS = 0;  // compute RMS 
  for(int i0 = 0; i0 < fNCount_; i0++) {
    fMean_ += fPupsPU_[i0]; // sum all the weights
    if(fPupsPU_[i0] == 0) continue;
    if(!fCharged_ && fAdjust_ && fPupsPU_[i0] > lMed) continue; // not charged, and fAdjust and we are greater than median ski
    lNRMS++; // denominatore
    fRMS_ += (fPupsPU_[i0]-lMed)*(fPupsPU_[i0]-lMed); // RMS with respect to the median value
  }

  fMean_/=fNCount_; // mean value
  if(lNRMS > 0) fRMS_/=lNRMS; // division
  if(fRMS_ == 0) fRMS_ = 1e-5; // set default value

  fRMS_ = sqrt(fRMS_); // take the sqrt
  fRMS_ *= fRMSScaleFactor_; // multiply by a scale factor

  if(!fAdjust_) return; // return otherwise adjust the weights

  //Adjust the p-value to correspond to the median
  std::sort(fPupsPV_.begin(),fPupsPV_.end()); // sort PV partucles
  int lNPV = 0;  // particles from PV
  for(unsigned int i0 = 0; i0 < fPupsPV_.size(); i0++){
    if(fPupsPV_[i0] <= lMed ) lNPV++;  // how many less than median from PV
  }
  double lAdjust = 1.5*double(lNPV)/double(fPupsPV_.size()+fNCount_);
  if(lAdjust > 0) fMedian_ -= sqrt(ROOT::Math::chisquared_quantile(lAdjust,1.)*fRMS_);
}

//This code is probably a bit confusing
double PuppiAlgo::compute(double & iVals, double iChi2) {  // compute the chi2 --> used as weight

  if(fAlgoId_ == -1) return 1;

  double lVal  = 0.;
  double lPVal = 1.;
  int    lNDOF = 0; 

  if(fNCount_ == 0) return 1.;   //in the NoPU case return 1.
  if(fCombId_ == 1) {  //Compute the previous p-value so that p-values can be multiplieed
      double pPVal = ROOT::Math::chisquared_cdf(lVal,lNDOF); // compute the chi2 of the weight as it was centered in zero--> probablity
      lPVal *= pPVal; // multiply for each 
      lNDOF = 0; 
      lVal  = 0; 
  }

  double pVal = iVals; // take the value of the single weight

  //Special Check for any algo with log(0) 
  if(fAlgoId_ == 0 && iVals == 0) pVal = fMedian_; // median
  if(fAlgoId_ == 3 && iVals == 0) pVal = fMedian_; // median
  if(fAlgoId_ == 5 && iVals == 0) pVal = fMedian_; // median
  lVal += (pVal-fMedian_)*(fabs(pVal-fMedian_))/fRMS_/fRMS_;  // chi2
  lNDOF++; // increase degree of fredom
  if( iChi2 != 0) lNDOF++;      //Add external Chi2 to first element
  if( iChi2 != 0) lVal+=iChi2;  //Add external Chi2 to first element 

  //Top it off with the last calc
  lPVal *= ROOT::Math::chisquared_cdf(lVal,lNDOF);
  return lPVal;
}

// get methods
double PuppiAlgo::getNeutralPt(int iNPV) { 
  return fNeutralPtMin_ + iNPV * fNeutralPtSlope_;
}

double PuppiAlgo::getPtMin() { 
  return fPtMin_;
}

double PuppiAlgo::getEtaMin() { 
  return fEtaMin_;
}

double PuppiAlgo::getEtaMax() { 
  return fEtaMax_;
}

int PuppiAlgo::getAlgoId() { 
  return fAlgoId_;
}

bool PuppiAlgo::getIsCharged() { 
  return fCharged_;
}

double PuppiAlgo::getConeSize() { 
  return fConeSize_;
}
