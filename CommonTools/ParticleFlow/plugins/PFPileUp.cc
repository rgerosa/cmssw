#include "CommonTools/ParticleFlow/plugins/PFPileUp.h"

using namespace std;
using namespace edm;
using namespace reco;

// constructor from event parameter set
template< class T >
PFPileUp<T>::PFPileUp(const edm::ParameterSet& iConfig) {

  tokenPFCandidates_      = consumes<PFCollection> (iConfig.getParameter<edm::InputTag>("PFCandidates")); // consumes of pfCandidates
  tokenPFCandidatesView_  = mayConsume<PFView>     (iConfig.getParameter<edm::InputTag>("PFCandidates"));

  tokenVertices_ = consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("Vertices"));

  enable_  = iConfig.exists("Enable") ? iConfig.getParameter<bool>("Enable") : true; 
  verbose_ = iConfig.exists("verbose") ? iConfig.getUntrackedParameter<bool>("verbose") : true;
  checkClosestZVertex_ = iConfig.exists("checkClosestZVertex") ? iConfig.getParameter<bool>("checkClosestZVertex") : true;

  // Configure the algo  
  pileUpAlgo_.setVerbose(verbose_);
  pileUpAlgo_.setCheckClosestZVertex(checkClosestZVertex_);

  softKillerParam_   = iConfig.getParameter<edm::ParameterSet>("softKillerParam");
  puppiParam_        = iConfig.getParameter<edm::ParameterSet>("puppiParam");
  jetCleansingParam_ = iConfig.getParameter<edm::ParameterSet>("jetCleansingParam");
  constituentSubtractionParam_ = iConfig.getParameter<edm::ParameterSet>("constituentSubtractionParam");

  produceCHS_                    = iConfig.exists("produceCHS")                    ? iConfig.getParameter<bool>("produceCHS") : true;
  produceSoftKiller_             = iConfig.exists("produceSoftKiller")             ? iConfig.getParameter<bool>("produceSoftKiller") : true;
  producePuppi_                  = iConfig.exists("producePuppi")                  ? iConfig.getParameter<bool>("producePuppi") : true;
  produceJetCleansing_           = iConfig.exists("produceJetCleansing")           ? iConfig.getParameter<bool>("produceJetCleansing") : true;
  produceConstituentSubtraction_ = iConfig.exists("produceConstituentSubtraction") ? iConfig.getParameter<bool>("produceConstituentSubtraction") : true;

  if(produceCHS_)          produces< PFCollection > ("pfCandidatePU"); // produce output CHS
  if(produceCHS_)          produces< PFCollection > ("pfCandidatePV"); // produce output CHS
  if(produceSoftKiller_)   produces< PFCollection > ("pfCandidateSoftKiller"); // produce output soft killer
  if(producePuppi_)        produces< PFCollection > ("pfCandidatePuppi"); // produce output puppi
  if(produceJetCleansing_) produces< PFCollection > ("pfCandidateJetCleansing"); // produce output jet cleansing
  if(produceConstituentSubtraction_) produces< PFCollection > ("pfCandidateConstSubtraction"); // produce output constituent subtraction

}


template< class T >
PFPileUp<T>::~PFPileUp(){} // deconstructor

// basic producer from event and setup
template< class T >
void PFPileUp<T>::produce(edm::Event& iEvent, const edm::EventSetup & iSetup) {

  std::auto_ptr< PFCollection > outputCHSCollectionPV ( new PFCollection );
  std::auto_ptr< PFCollection > outputCHSCollectionPU ( new PFCollection );
  std::auto_ptr< PFCollection > outputSoftKillerCollection ( new PFCollection );
  std::auto_ptr< PFCollection > outputPuppiCollection ( new PFCollection );
  std::auto_ptr< PFCollection > outputJetCleansingCollection ( new PFCollection );
  std::auto_ptr< PFCollection > outputConstituentSubtractionCollection ( new PFCollection );

  if(enable_) {

    // get vertices
    edm::Handle<VertexCollection> vertices;
    iEvent.getByToken(tokenVertices_,vertices);
    // try as they are miniAOD
    edm::Handle<PFCollection> pfCandidates;
    const PFCollection* pfCandidatesRef = 0;
    PFCollection useNoFwdPtrs;

    bool getFromPackedFwdPtr = iEvent.getByToken(tokenPFCandidates_, pfCandidates);
    if ( getFromPackedFwdPtr ) { pfCandidatesRef = pfCandidates.product();
    }
    else {
       edm::Handle<PFView> pfView;
       bool getFromPFView = iEvent.getByToken(tokenPFCandidatesView_, pfView);       
       if(getFromPFView){ 
         for (typename edm::View<T>::const_iterator viewBegin = pfView->begin(), viewEnd = pfView->end(), iview = viewBegin; iview != viewEnd; ++iview ) {            
	   useNoFwdPtrs.push_back(edm::FwdPtr<T>(pfView->ptrAt(iview-viewBegin), pfView->ptrAt(iview-viewBegin)));
	 }
       pfCandidatesRef = &useNoFwdPtrs;
       }
    }
    
    if(pfCandidatesRef == 0) {
      throw cms::Exception("Something went dreadfully wrong with PFPileUp. pfCandidatesRef should never be zero, so this is a logic error.");
    }

    if(produceCHS_){ // if you want to run CHS
      pileUpAlgo_.processChargedHadronSubtraction(*pfCandidatesRef,*vertices.product()); // run CHS alone
      outputCHSCollectionPV->insert(outputCHSCollectionPV->end(),pileUpAlgo_.getPFCandidatesFromVtx().begin(),pileUpAlgo_.getPFCandidatesFromVtx().end()); 
      outputCHSCollectionPU->insert(outputCHSCollectionPU->end(),pileUpAlgo_.getPFCandidatesFromPU().begin(),pileUpAlgo_.getPFCandidatesFromPU().end()); 
    }

    if(produceSoftKiller_){ // if you want to run softkiller
      pileUpAlgo_.processSoftKiller(*pfCandidatesRef,*vertices,softKillerParam_); // run CHS alone
      outputSoftKillerCollection->insert(outputSoftKillerCollection->end(),pileUpAlgo_.getSoftKillerPFCandidates().begin(),pileUpAlgo_.getSoftKillerPFCandidates().end()); 
    }

    if(produceJetCleansing_){ // if you want to run jet cleansing
      pileUpAlgo_.processJetCleansing(*pfCandidatesRef,*vertices,jetCleansingParam_); // run CHS alone
      outputJetCleansingCollection->insert(outputJetCleansingCollection->end(),pileUpAlgo_.getJetCleansingPFCandidates().begin(),pileUpAlgo_.getJetCleansingPFCandidates().end()); 
    }

    if(produceConstituentSubtraction_){ // if you want to run constituent subtraction
      pileUpAlgo_.processConstituentSubtraction(*pfCandidatesRef,*vertices,constituentSubtractionParam_); // run CHS alone
      outputConstituentSubtractionCollection->insert(outputConstituentSubtractionCollection->end(),pileUpAlgo_.getConstituentSubtractionPFCandidates().begin(),pileUpAlgo_.getConstituentSubtractionPFCandidates().end()); 
    }

    if(producePuppi_){ // if you want to run puppi
      pileUpAlgo_.processPuppi(*pfCandidatesRef,*vertices,puppiParam_); // run CHS alone
      outputPuppiCollection->insert(outputPuppiCollection->end(),pileUpAlgo_.getPuppiPFCandidates().begin(),pileUpAlgo_.getPuppiPFCandidates().end()); 
    }


  } // end if enabled

  // outsize of the loop to fill the collection anyway even when disabled
  if(produceCHS_){
    iEvent.put(outputCHSCollectionPV,"pfCandidatePV");
    iEvent.put(outputCHSCollectionPU,"pfCandidatePU");
  }
  
  if(produceSoftKiller_) iEvent.put(outputSoftKillerCollection,"pfCandidateSoftKiller");

  if(producePuppi_) iEvent.put(outputPuppiCollection,"pfCandidatePuppi");

  if(produceJetCleansing_) iEvent.put(outputJetCleansingCollection,"pfCandidateJetCleansing");

  if(produceConstituentSubtraction_) iEvent.put(outputConstituentSubtractionCollection,"pfCandidateConstSubtraction");

}

