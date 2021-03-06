#ifndef FWCore_Framework_stream_EDProducerAdaptorBase_h
#define FWCore_Framework_stream_EDProducerAdaptorBase_h
// -*- C++ -*-
//
// Package:     FWCore/Framework
// Class  :     EDProducerAdaptorBase
// 
/**\class edm::stream::EDProducerAdaptorBase EDProducerAdaptorBase.h "FWCore/Framework/interface/stream/EDProducerAdaptorBase.h"

 Description: [one line class summary]

 Usage:
    <usage>

*/
//
// Original Author:  Chris Jones
//         Created:  Fri, 02 Aug 2013 18:09:15 GMT
//

// system include files

// user include files
#include "FWCore/Framework/interface/stream/ProducingModuleAdaptorBase.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Provenance/interface/ModuleDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/RunIndex.h"
#include "FWCore/Utilities/interface/LuminosityBlockIndex.h"


// forward declarations

namespace edm {
  namespace stream {
    class EDProducerBase;
    class EDProducerAdaptorBase : public ProducingModuleAdaptorBase<EDProducerBase>
    {
      
    public:
      template <typename T> friend class edm::WorkerT;

      EDProducerAdaptorBase();
      
      // ---------- const member functions ---------------------
      
      // ---------- static member functions --------------------
      
      // ---------- member functions ---------------------------
      
      std::string workerType() const { return "WorkerT<EDProducerAdaptorBase>";}
    protected:
      using ProducingModuleAdaptorBase<EDProducerBase>::commit;

    private:
      EDProducerAdaptorBase(const EDProducerAdaptorBase&) =delete; // stop default
      
      const EDProducerAdaptorBase& operator=(const EDProducerAdaptorBase&) =delete; // stop default
      
      bool doEvent(EventPrincipal& ep, EventSetup const& c,
                           CurrentProcessingContext const* cpcp) ;
    };
  }
}

#endif
