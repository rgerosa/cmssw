#ifndef CondFormats_SiStripObjects_PedsFullNoiseAnalysis_H
#define CondFormats_SiStripObjects_PedsFullNoiseAnalysis_H

#include "CondFormats/SiStripObjects/interface/CommissioningAnalysis.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include <boost/cstdint.hpp>
#include <sstream>
#include <vector>

/** 
    @class PedsFullNoiseAnalysis
    @author M. Wingham, R.Bainbridge
    @brief Histogram-based analysis for pedestal run.
*/

class PedsFullNoiseAnalysis : public CommissioningAnalysis {
  
 public:

	// ---------- con(de)structors ----------

    PedsFullNoiseAnalysis( const uint32_t& key );

    PedsFullNoiseAnalysis();

    virtual ~PedsFullNoiseAnalysis() {;}

    friend class PedestalsAlgorithm;
    friend class PedsFullNoiseAlgorithm;

    // ---------- public interface ----------
    
    /** Identifies if analysis is valid or not. */
    bool isValid() const;
    
    // Pedestal, noise and raw noise (128-strip vector per APV)
    inline const VVFloat& peds() const;
    inline const VVFloat& noise() const;
    inline const VVFloat& raw() const;

    // test statistics for each APV (128-strip vector per APV)
    inline const VVFloat& adProbab() const;
    inline const VVFloat& ksProbab() const;
    inline const VVFloat& jbProbab() const;
    inline const VVFloat& chi2Probab() const;

    // Per strip values
    inline const VVFloat& noiseRMS() const; // RMS
    inline const VVFloat& noiseSigmaGaus() const; // from gaus fit
    inline const VVFloat& noiseSignificance() const; // noise significance
    inline const VVFloat& noiseBin84() const;
    inline const VVFloat& noiseSkewness() const;
    inline const VVFloat& noiseKurtosis() const;
    inline const VVFloat& noiseIntegralNsigma() const;
    inline const VVFloat& noiseIntegral() const;
        
    // status for different class of bad or problematic strips
    inline const VVInt& deadStrip() const; 
    inline const VVInt& badStrip() const;
    inline const VVInt& shiftedStrip() const;
    inline const VVInt& lowNoiseStrip() const;
    inline const VVInt& largeNoiseStrip() const;
    inline const VVInt& largeNoiseSignificance() const;
    inline const VVInt& badFitStatus() const;
    inline const VVInt& badADProbab() const;
    inline const VVInt& badKSProbab() const;
    inline const VVInt& badJBProbab() const;
    inline const VVInt& badChi2Probab() const;
    inline const VVInt& badTailStrip() const;
    inline const VVInt& badDoublePeakStrip() const;

    // Mean and rms spread (value per APV)
    inline const VFloat& pedsMean() const;
    inline const VFloat& pedsSpread() const;
    inline const VFloat& noiseMean() const;
    inline const VFloat& noiseSpread() const;
    inline const VFloat& rawMean() const;
    inline const VFloat& rawSpread() const;

    // Max and min values (value per APV)
    inline const VFloat& pedsMax() const;
    inline const VFloat& pedsMin() const; 
    inline const VFloat& noiseMax() const;
    inline const VFloat& noiseMin() const;
    inline const VFloat& rawMax() const;
    inline const VFloat& rawMin() const;

    // ---------- misc ----------

    /** Prints analysis results. */
    void print( std::stringstream&, uint32_t apv_number = 0 );

    /** Overrides base method. */
    void summary( std::stringstream& ) const;

    /** Resets analysis member data. */
    void reset();

    // ---------- private member data ----------

  private:
	
    // VVFloats means: 1 vector per APV, 1 value per strip.

    VVFloat peds_;
    VVFloat noise_;
    VVFloat raw_;
    
    VVFloat adProbab_;
    VVFloat ksProbab_;
    VVFloat jbProbab_;
    VVFloat chi2Probab_;
    VVFloat noiseRMS_;
    VVFloat noiseSigmaGaus_;
    VVFloat noiseSignificance_;
    VVFloat noiseBin84_;
    VVFloat noiseSkewness_;
    VVFloat noiseKurtosis_;
    VVFloat noiseIntegralNsigma_;
    VVFloat noiseIntegral_;

    VVInt deadStrip_;
    VVInt badStrip_;
    VVInt shiftedStrip_;
    VVInt lowNoiseStrip_;
    VVInt largeNoiseStrip_;
    VVInt largeNoiseSignificance_;
    VVInt badFitStatus_;
    VVInt badADProbab_;
    VVInt badKSProbab_;
    VVInt badJBProbab_;
    VVInt badChi2Probab_;
    VVInt badTailStrip_;
    VVInt badDoublePeakStrip_;

    // VFloat: 1 value per APV
    VFloat pedsMean_;
    VFloat pedsSpread_;
    VFloat noiseMean_;
    VFloat noiseSpread_;
    VFloat rawMean_;
    VFloat rawSpread_;
    VFloat pedsMax_;
    VFloat pedsMin_; 
    VFloat noiseMax_;
    VFloat noiseMin_;
    VFloat rawMax_;
    VFloat rawMin_;
    
    // true if legacy histogram naming is used
    bool legacy_;
};

// ---------- Inline methods ----------

const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::peds() const { return peds_; }
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::noise() const { return noise_; }
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::raw() const { return raw_; }

const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::adProbab() const { return adProbab_;}
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::ksProbab() const { return ksProbab_;}
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::jbProbab() const { return jbProbab_;}
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::chi2Probab() const { return chi2Probab_;}

const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::noiseRMS() const { return noiseRMS_;}
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::noiseSigmaGaus() const { return noiseSigmaGaus_;}
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::noiseSignificance() const { return noiseSignificance_;}
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::noiseBin84() const { return noiseBin84_;}
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::noiseSkewness() const { return noiseSkewness_;}
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::noiseKurtosis() const { return noiseKurtosis_;}
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::noiseIntegralNsigma() const { return noiseIntegralNsigma_;}
const PedsFullNoiseAnalysis::VVFloat& PedsFullNoiseAnalysis::noiseIntegral() const { return noiseIntegral_;}

const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::deadStrip() const { return deadStrip_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::badStrip() const { return badStrip_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::shiftedStrip() const { return shiftedStrip_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::lowNoiseStrip() const { return lowNoiseStrip_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::largeNoiseStrip() const { return largeNoiseStrip_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::largeNoiseSignificance() const { return largeNoiseSignificance_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::badFitStatus() const { return badFitStatus_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::badADProbab() const { return badADProbab_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::badKSProbab() const { return badKSProbab_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::badJBProbab() const { return badJBProbab_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::badChi2Probab() const { return badChi2Probab_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::badTailStrip() const { return badTailStrip_;}
const PedsFullNoiseAnalysis::VVInt& PedsFullNoiseAnalysis::badDoublePeakStrip() const { return badDoublePeakStrip_;}
        
const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::pedsMean() const { return pedsMean_; }
const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::pedsSpread() const { return pedsSpread_; }
const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::noiseMean() const { return noiseMean_; }
const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::noiseSpread() const { return noiseSpread_; }
const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::rawMean() const { return rawMean_; }
const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::rawSpread() const { return rawSpread_; }

const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::pedsMax() const { return pedsMax_; }
const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::pedsMin() const { return pedsMin_; } 
const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::noiseMax() const { return noiseMax_; }
const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::noiseMin() const { return noiseMin_; }
const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::rawMax() const { return rawMax_; }
const PedsFullNoiseAnalysis::VFloat& PedsFullNoiseAnalysis::rawMin() const { return rawMin_; }

#endif // CondFormats_SiStripObjects_PedsFullNoiseAnalysis_H


