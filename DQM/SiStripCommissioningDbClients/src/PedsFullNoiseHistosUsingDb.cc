
#include "DQM/SiStripCommissioningDbClients/interface/PedsFullNoiseHistosUsingDb.h"
#include "CondFormats/SiStripObjects/interface/PedsFullNoiseAnalysis.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "DataFormats/SiStripCommon/interface/SiStripFecKey.h"
#include "DataFormats/SiStripCommon/interface/SiStripFedKey.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>

using namespace sistrip;

// -----------------------------------------------------------------------------
/** */
PedsFullNoiseHistosUsingDb::PedsFullNoiseHistosUsingDb( const edm::ParameterSet & pset,
                                                        DQMStore* bei,
                                                        SiStripConfigDb* const db ) 
  : CommissioningHistograms(pset.getParameter<edm::ParameterSet>("PedsFullNoiseParameters"),bei,sistrip::PEDS_FULL_NOISE ),
    CommissioningHistosUsingDb(db,sistrip::PEDS_FULL_NOISE ),
    PedsFullNoiseHistograms( pset.getParameter<edm::ParameterSet>("PedsFullNoiseParameters"),bei )
{
  LogTrace(mlDqmClient_) 
    << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";

  highThreshold_ = this->pset().getParameter<double>("HighThreshold");
  lowThreshold_  = this->pset().getParameter<double>("LowThreshold");

  LogTrace(mlDqmClient_)
    << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]"
    << " Set FED zero suppression high/low threshold to "
    << highThreshold_ << "/" << lowThreshold_;

  disableBadStrips_   = this->pset().getParameter<bool>("DisableBadStrips");
  keepStripsDisabled_ = this->pset().getParameter<bool>("KeepStripsDisabled");
  skipEmptyStrips_    = this->pset().getParameter<bool>("skipEmptyStrips");
  uploadOnlyStripMaskingBit_ = this->pset().getParameter<bool>("UploadOnlyStripMaskingBit");
  
  LogTrace(mlDqmClient_)
    << "[PedestalsHistosUsingDb::" << __func__ << "]"
    << " Disabling strips: " << disableBadStrips_
    << " ; keeping previously disabled strips: " << keepStripsDisabled_
    << " ; skip strips with no data: " << skipEmptyStrips_
    << " ; upload only masking bit: " << uploadOnlyStripMaskingBit_;
    
}

// -----------------------------------------------------------------------------
/** */
PedsFullNoiseHistosUsingDb::~PedsFullNoiseHistosUsingDb() {
  LogTrace(mlDqmClient_) 
    << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]"
    << " Destructing object...";
}

// -----------------------------------------------------------------------------
/** */
void PedsFullNoiseHistosUsingDb::uploadConfigurations() {

  LogTrace(mlDqmClient_) 
    << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]";

  if ( !db() ) {
    edm::LogError(mlDqmClient_) 
      << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]"
      << " NULL pointer to SiStripConfigDb interface!"
      << " Aborting upload...";
    return;
  }
  
  // Update FED descriptions with new peds/noise values
  SiStripConfigDb::FedDescriptionsRange feds = db()->getFedDescriptions(); 
  update( feds );
  if (doUploadConf()) {  // check whether the upload HD config is set to true

    edm::LogVerbatim(mlDqmClient_) 
      << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]"
      << " Uploading pedestals/noise to DB...";

    db()->uploadFedDescriptions(); // change the FED version

    edm::LogVerbatim(mlDqmClient_) 
      << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]"
      << " Completed database upload of " << feds.size() 
      << " FED descriptions!";
  } 
  else {
    edm::LogWarning(mlDqmClient_) 
      << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]"
      << " TEST! No pedestals/noise values will be uploaded to DB...";
  }  
}

// -----------------------------------------------------------------------------
/** */
void PedsFullNoiseHistosUsingDb::update( SiStripConfigDb::FedDescriptionsRange feds ) {
 
  // Iterate through feds and update fed descriptions
  uint16_t updated = 0;
  SiStripConfigDb::FedDescriptionsV::const_iterator ifed;

  for ( ifed = feds.begin(); ifed != feds.end(); ifed++ ) {    

    for ( uint16_t ichan = 0; ichan < sistrip::FEDCH_PER_FED; ichan++ ) {
      // Build FED and FEC keys from the cabling object i.e. checking if there is a connection
      const FedChannelConnection& conn = cabling()->fedConnection( (*ifed)->getFedId(), ichan );
      if ( conn.fecCrate() == sistrip::invalid_ ||
           conn.fecSlot()  == sistrip::invalid_ ||
           conn.fecRing()  == sistrip::invalid_ ||
           conn.ccuAddr()  == sistrip::invalid_ ||
           conn.ccuChan()  == sistrip::invalid_ ||
           conn.lldChannel() == sistrip::invalid_ ) 
	continue; 

      // build the FED and FEC key from the connection object
      SiStripFedKey fed_key( conn.fedId(), 
                             SiStripFedKey::feUnit( conn.fedCh() ),
                             SiStripFedKey::feChan( conn.fedCh() ) );

      SiStripFecKey fec_key( conn.fecCrate(),
                             conn.fecSlot(),
                             conn.fecRing(),
                             conn.ccuAddr(),
                             conn.ccuChan(),
                             conn.lldChannel() );

      // Locate appropriate analysis object --> based on FEC keys cause they are per lldChannel
      Analyses::const_iterator iter = data().find( fec_key.key() );

      // if data are found from the analysis
      if ( iter != data().end() ) {      

        PedsFullNoiseAnalysis* anal = dynamic_cast<PedsFullNoiseAnalysis*>( iter->second );

        if ( !anal ) { 
          edm::LogError(mlDqmClient_)
            << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]"
            << " NULL pointer to analysis object!";
          continue; 
        }

        // Determine the pedestal shift to apply --> this is standard in the pedestal paylaod to avoid loss of signal from common-mode subtraction
        uint32_t pedshift = 127;
        for ( uint16_t iapv = 0; iapv < sistrip::APVS_PER_FEDCH; iapv++ ) {
          uint32_t pedmin = (uint32_t) anal->pedsMin()[iapv];
          pedshift = pedmin < pedshift ? pedmin : pedshift;
        }

        // Iterate through APVs and strips
        for ( uint16_t iapv = 0; iapv < sistrip::APVS_PER_FEDCH; iapv++ ) {
          for ( uint16_t istr = 0; istr < anal->peds()[iapv].size(); istr++ ) { // Loop on the pedestal for each APV

            // get the information on the strip as it was on the db
            Fed9U::Fed9UAddress addr( ichan, iapv, istr );
            Fed9U::Fed9UStripDescription temp = (*ifed)->getFedStrips().getStrip( addr );

            // determine whether we need to disable the strip
            bool disableStrip = false;
	    std::stringstream ss_disable;

	    if(temp.getDisable()) { // strip already disabled in the database
	      ss_disable<<"Already Disabled: "<<conn.fecCrate()
			<<" "<<conn.fecSlot()
			<<" "<<conn.fecRing()
			<<" "<<conn.ccuAddr()
			<<" "<<conn.ccuChan()
			<<" "<<conn.lldChannel()
			<<" "<<iapv*128+istr<<std::endl;

	      if(keepStripsDisabled_) disableStrip = true; // in case one wants to keep them disabled
            }
	    else{
	      
	      // to disable new strips
	      if(disableBadStrips_){
		
		SiStripFedKey fed_key(anal->fedKey());              
              	PedsFullNoiseAnalysis::VInt dead = anal->deadStrip()[iapv];
              	if (not skipEmptyStrips_ and  // if one don't want to skip dead strips
		    find( dead.begin(), dead.end(), istr ) != dead.end() ) {
		  disableStrip = true;
                  ss_disable<<"Disabling Dead Strip: "<<conn.fecCrate()
			    <<" "<<conn.fecSlot()
			    <<" "<<conn.fecRing()
			    <<" "<<conn.ccuAddr()
			    <<" "<<conn.ccuChan()
			    <<" "<<conn.lldChannel()
			    <<" "<<iapv*128+istr<<std::endl;
                }
		
              	PedsFullNoiseAnalysis::VInt badcChan = anal->badStrip()[iapv]; // new feature --> this is the sample of the whole bad strips from the analysis
              	if ( find( badcChan.begin(), badcChan.end(), istr ) != badcChan.end() ) {
		  disableStrip = true;
                  ss_disable<<"Disabling Bad strip: "<<conn.fecCrate()
			    <<" "<<conn.fecSlot()
			    <<" "<<conn.fecRing()
			    <<" "<<conn.ccuAddr()
			    <<" "<<conn.ccuChan()
			    <<" "<<conn.lldChannel()
			    <<" "<<iapv*128+istr<<std::endl;
                }
              }
            }

	    if(edm::isDebugEnabled())
	      LogTrace(mlDqmClient_) << ss_disable.str();
	    
	    uint32_t pedestalVal = 0;
	    uint32_t noiseVal = 0;
	    uint32_t lowThr   = 0;
	    uint32_t highThr  = 0;
	    
	    // download the previous pedestal/noise payload from the DB
	    if(uploadOnlyStripMaskingBit_){ 
	      // understand and implement how to download the payload	      
	      pedestalVal = temp.getPedestal();
	      noiseVal = temp.getNoise();
	      lowThr   = temp.getLowThreshold();
	      highThr  = temp.getHighThreshold();
	    }
	    else{	      
	      pedestalVal = static_cast<uint32_t>( anal->peds()[iapv][istr]-pedshift );
	      noiseVal = static_cast<uint32_t>(anal->noise()[iapv][istr]);
	      lowThr   = lowThreshold_;
	      highThr  = highThreshold_;
	    }
	    
	    
            Fed9U::Fed9UStripDescription data(pedestalVal,highThr,lowThr,noiseVal,disableStrip);

            std::stringstream ss;
            if ( data.getDisable() && edm::isDebugEnabled() ) {
              ss << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]"
                 << " Disabling strip in Fed9UStripDescription object..." << std::endl
                 << " for FED id/channel and APV/strip : "
                 << fed_key.fedId() << "/"
                 << fed_key.fedChannel() << " "
                 << iapv << "/"
                 << istr << std::endl 
                 << " and crate/FEC/ring/CCU/module    : "
                 << fec_key.fecCrate() << "/"
                 << fec_key.fecSlot() << "/"
                 << fec_key.fecRing() << "/"
                 << fec_key.ccuAddr() << "/"
                 << fec_key.ccuChan() << std::endl 
                 << " from ped/noise/high/low/disable  : "
                 << static_cast<uint16_t>( temp.getPedestal() ) << "/" 
                 << static_cast<uint16_t>( temp.getHighThreshold() ) << "/" 
                 << static_cast<uint16_t>( temp.getLowThreshold() ) << "/" 
                 << static_cast<uint16_t>( temp.getNoise() ) << "/" 
                 << static_cast<uint16_t>( temp.getDisable() ) << std::endl;
            }
            (*ifed)->getFedStrips().setStrip( addr, data );
            if ( data.getDisable() && edm::isDebugEnabled() ) {
              ss << " to ped/noise/high/low/disable    : "
                 << static_cast<uint16_t>( data.getPedestal() ) << "/" 
                 << static_cast<uint16_t>( data.getHighThreshold() ) << "/" 
                 << static_cast<uint16_t>( data.getLowThreshold() ) << "/" 
                 << static_cast<uint16_t>( data.getNoise() ) << "/" 
                 << static_cast<uint16_t>( data.getDisable() ) << std::endl;
              LogTrace(mlDqmClient_) << ss.str();
            }	    
          } // end loop on strips
        } // end loop on apvs
        updated++;
      }

      else { // device not found in the analysis	
        if ( deviceIsPresent(fec_key) ) {
          edm::LogWarning(mlDqmClient_) 
            << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]"
            << " Unable to find pedestals/noise for FedKey/Id/Ch: " 
            << hex << setw(8) << setfill('0') << fed_key.key() << dec << "/"
            << (*ifed)->getFedId() << "/"
            << ichan
            << " and device with FEC/slot/ring/CCU/LLD " 
            << fec_key.fecCrate() << "/"
            << fec_key.fecSlot() << "/"
            << fec_key.fecRing() << "/"
            << fec_key.ccuAddr() << "/"
            << fec_key.ccuChan() << "/"
            << fec_key.channel();
        }
      }
    }
  }
  
  edm::LogVerbatim(mlDqmClient_) 
    << "[PedsFullNoiseHistosUsingDb::" << __func__ << "]"
    << " Updated FED pedestals/noise for " 
    << updated << " channels";
}

// -----------------------------------------------------------------------------
/** */
void PedsFullNoiseHistosUsingDb::create( SiStripConfigDb::AnalysisDescriptionsV& desc,
                                         Analysis analysis ) {

  PedsFullNoiseAnalysis* anal = dynamic_cast<PedsFullNoiseAnalysis*>( analysis->second );
  if ( !anal ) { return; }
  
  SiStripFecKey fec_key( anal->fecKey() );
  SiStripFedKey fed_key( anal->fedKey() );
  
  for ( uint16_t iapv = 0; iapv < 2; ++iapv ) {
    // Create description
    PedestalsAnalysisDescription* tmp = NULL;
    /*
    tmp = new PedestalsAnalysisDescription(
					   //// Bad flags for the analysis summary
					   anal->deadStrip()[iapv],
					   anal->badStrip()[iapv],
					   anal->shiftedStrip()[iapv],
					   anal->lowNoiseStrip()[iapv],
					   anal->largeNoiseStrip()[iapv],
					   anal->largeNoiseSignificance()[iapv],
					   anal->badFitStatus()[iapv],
					   anal->badADProbab()[iapv],
					   anal->badKSProbab()[iapv],
					   anal->badJBProbab()[iapv],
					   anal->badChi2Probab()[iapv],
					   anal->badTailStrip()[iapv],
					   anal->badDoublePeakStrip()[iapv],					   
					   ///// Per APV quantities
					   anal->pedsMean()[iapv],
					   anal->pedsSpread()[iapv],
					   anal->noiseMean()[iapv],
					   anal->noiseSpread()[iapv],
					   anal->rawMean()[iapv],
					   anal->rawSpread()[iapv],
					   anal->pedsMax()[iapv], 
					   anal->pedsMin()[iapv], 
					   anal->noiseMax()[iapv],
					   anal->noiseMin()[iapv],
					   anal->rawMax()[iapv],
					   anal->rawMin()[iapv],
					   ///// --> test statistic values
					   anal->adProbab()[iapv],
					   anal->ksProbab()[iapv],
					   anal->jbProbab()[iapv],
					   anal->chi2Probab()[iapv],
					   //// --> Per strip quantities
					   anal->noiseRMS()[iapv],
					   anal->noiseSigmaGaus()[iapv],
					   anal->noiseSignificance()[iapv],
					   anal->noiseBin84()[iapv],
					   anal->residualSkewness()[iapv],
					   anal->residualKurtosis()[iapv],
					   anal->residualIntegralNsigma()[iapv],
					   anal->residualIntegral()[iapv],
					   ///// coordinates
					   fec_key.fecCrate(),
					   fec_key.fecSlot(),
					   fec_key.fecRing(),
					   fec_key.ccuAddr(),
					   fec_key.ccuChan(),
					   SiStripFecKey::i2cAddr( fec_key.lldChan(), !iapv ), 
					   db()->dbParams().partitions().begin()->second.partitionName(),
					   db()->dbParams().partitions().begin()->second.runNumber(),
					   anal->isValid(),
					   "",
					   fed_key.fedId(),
					   fed_key.feUnit(),
					   fed_key.feChan(),
					   fed_key.fedApv()
					   );
    
    */
    // Add comments
    typedef std::vector<std::string> Strings;
    Strings errors = anal->getErrorCodes();
    Strings::const_iterator istr = errors.begin();
    Strings::const_iterator jstr = errors.end();
    for ( ; istr != jstr; ++istr ) { tmp->addComments( *istr ); }
    // Store description
    desc.push_back( tmp );      
  }

}

