#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"

// alternative constructor
template< typename T>
PFPileUpAlgo<T>::PFPileUpAlgo(const CandCollection & pfCandidates, const reco::VertexCollection & vertices, bool checkClosestZVertex, bool verbose){

  checkClosestZVertex_ = checkClosestZVertex;
  verbose_             = verbose;
  vertices_            = &vertices ;
  pfCandidates_        = pfCandidates;
  processChargedHadronSubtraction(pfCandidates_,*vertices_); // process charged hadron subtraction by default

}

template< typename T>
PFPileUpAlgo<T>::~PFPileUpAlgo(){

  pfCandidates_.clear();               // all pf candidates
  pfCandidatesFromVtx_.clear();        // all pf candidates from PV (charged + neutral)
  pfCandidatesFromPU_.clear();         // charge candidates from PU 
  pfCandidatesFromVtxCharged_.clear(); // only charged particles associated to PV

  pfSoftKiller_.clear() ;      // soft killer particles                                                                                                                       

  pfJetCleansing_.clear() ;    // jet cleansing particles                                                                                                                  

  pfConstituentSubtraction_.clear() ; // constituent subtraction particles

  pfPuppi_.clear();            // PUPPI particles                                                                                                                                      
  fPuppiAlgo_.clear() ;        // vector of puppi algorithm objects                                                                                                                
  puppiMetricValues_.clear() ; // metric values
  puppiWeights_.clear() ;      // weights
  lAlgos_.clear();             // algorithms
  particleIndex_.clear();

  if(!vertices_) delete vertices_;
}

template< typename T>
void PFPileUpAlgo<T>::processChargedHadronSubtraction(const CandCollection & pfCandidates, 
			                              const reco::VertexCollection & vertices,
                                                      const double & DeltaZCut){
  
  pfCandidatesFromVtx_.clear(); // clear existing collections
  pfCandidatesFromPU_.clear();  
  pfCandidatesFromVtxCharged_.clear();

  for( unsigned i = 0; i < pfCandidates.size(); i++ ) { // loop on PFcandidates    
     int ivertex = chargedHadronVertex(vertices, *pfCandidates.at(i)); // identify the vertex          
     // no associated vertex, or primary vertex  --> put in the fromVtx list     
     if( ivertex == -1  || ivertex == 0 ) {
      if(verbose_) std::cout<<"VTX "<<i<<" "<<*pfCandidates[i]<<std::endl;
      bool isGood = true ;     
      if(ivertex == 0 and fabs(pfCandidates[i]->vertex().z()-vertices[0].z()) > DeltaZCut) isGood = false;
      if(isGood){
         pfCandidatesFromVtx_.push_back(pfCandidates[i]);
         particleIndex_.push_back(0);
         if(pfCandidates[i]->charge()!=0){ pfCandidatesFromVtxCharged_.push_back(pfCandidates[i]); particleIndex_.back() = 1 ;}      
      }     
      else{     
        pfCandidatesFromPU_.push_back(pfCandidates[i]);
	particleIndex_.push_back(2);
      }
     }
     else{
        pfCandidatesFromPU_.push_back(pfCandidates[i]);
	particleIndex_.push_back(2);
     }
  }
}


// method in prder to run soft killer
template <typename T>
void PFPileUpAlgo<T>::processSoftKiller(const edm::ParameterSet& param){

  processSoftKiller(pfCandidates_,*vertices_,param);

}

template <typename T>
void PFPileUpAlgo<T>::processSoftKiller(const CandCollection & pfCandidates,
				        const reco::VertexCollection & vertices,
				        const edm::ParameterSet& param){

  pfSoftKiller_.clear(); // clean soft killer particle list

  // soft killer definition                
  fastjet::contrib::SoftKiller softKiller ( param.exists("ymax") ? param.getParameter<double>("ymax") : ymax, 
                                            param.exists("cell_size") ? param.getParameter<double>("cell_size") : cell_size); 

  // select the event particle to give as input to softkiller  
  std::vector<fastjet::PseudoJet> eventParticles;
  bool applyCHS = false ;
  if(param.exists("applyCHS")){
    if(param.getParameter<bool>("applyCHS")){ // if CHS has to be applied, cross check that the collection are not empty, in this case pfCandidate vertex association is not re-run
     applyCHS = true ;
     if(pfCandidatesFromVtx_.empty() or pfCandidatesFromPU_.empty()) 
       processChargedHadronSubtraction(pfCandidates,vertices);
     eventParticles = ConvertToPseudoJet(pfCandidatesFromVtx_); // cluster jet on top of CHS candidates
    }
    else eventParticles = ConvertToPseudoJet(pfCandidates); // use full pfCandidates
  }
  else eventParticles = ConvertToPseudoJet(pfCandidates); // use full pfCandidates

  // apply soft killer to particles in the event
  std::vector<fastjet::PseudoJet> softParticles  = softKiller(eventParticles); // take softkiller particles

  if(param.exists("clusterInJets")){ 
    if(param.getParameter<bool>("clusterInJets")){ // if yes, cluster soft particles in jets and take the constituent after a possible pt cut applied on jets

     // define a jet type for softkiller given a algo an dR as softKiller parameter set --> by default used some static info in private part of this class 
     fastjet::JetDefinition jet_def( param.exists("jetAlgo") ? get_algo(param.getParameter<std::string>("jetAlgo")) : fastjet::antikt_algorithm , 
                                   param.exists("jetR")? param.getParameter<double>("jetR") : jetR);

     // area definition with explicit ghosts 
     fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,fastjet::GhostedAreaSpec(fastjet::SelectorAbsRapMax(param.exists("GhostEtaMax") ? param.getParameter<double>("GhostEtaMax") : GhostEtaMax))); 

     fastjet::ClusterSequenceArea pSoftJets(softParticles,jet_def,area_def); // cluster them 
     std::vector<fastjet::PseudoJet> softJets = sorted_by_pt(pSoftJets.inclusive_jets(param.exists("jetPtCut") ? param.getParameter<double>("jetPtCut") : 0)); // take inclusvie jets
     
     if(applyCHS) pfSoftKiller_ = ConvertJetToCandCollection(softJets,pfCandidatesFromVtx_); // give jets as input of converting function and look at Vtx particles
     else         pfSoftKiller_ = ConvertJetToCandCollection(softJets,pfCandidates); // give jets as input of converting function

    }
    else{

      if(applyCHS) pfSoftKiller_ = ConvertToCandCollection(softParticles,pfCandidatesFromVtx_);
      else         pfSoftKiller_ = ConvertToCandCollection(softParticles,pfCandidates); // give particles as input to the converting function
    }    
  }
  else{ if(applyCHS) pfSoftKiller_ = ConvertToCandCollection(softParticles,pfCandidatesFromVtx_);
        else         pfSoftKiller_ = ConvertToCandCollection(softParticles,pfCandidates); // give particles as input to the converting function
  }
}
  
// implement jet cleaning
template <typename T>
void PFPileUpAlgo<T>::processJetCleansing(const edm::ParameterSet& param){

  processJetCleansing(pfCandidates_,*vertices_,param);

}

template <typename T>
void PFPileUpAlgo<T>::processJetCleansing(const CandCollection & pfCandidates,
		 		          const reco::VertexCollection & vertices,
				          const edm::ParameterSet& param){

  pfJetCleansing_.clear();  // clean the PF collection

  // define a benchmarck jet definition for cleansing  
  fastjet::JetDefinition jet_def(param.exists("jetAlgo")? get_algo(param.getParameter<std::string>("jetAlgo")) : fastjet::antikt_algorithm, 
                                 param.exists("jetR") ? param.getParameter<double>("jetR") : jetR);

  // define a benchmarck sub-jet definition for cleansing  
  fastjet::JetDefinition subjet_def(param.exists("subjetAlgo") ? get_algo(param.getParameter<std::string>("subjetAlgo")) : fastjet::kt_algorithm, 
                                    param.exists("subjetR") ? param.getParameter<double>("subjetR") : subjetR);
  
  // throw ghosts
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,fastjet::GhostedAreaSpec(fastjet::SelectorAbsRapMax(param.exists("GhostEtaMax") ? param.getParameter<double>("GhostEtaMax") : GhostEtaMax)));   

  // cleanser definition
  fastjet::contrib::JetCleanser::input_mode projectmode = fastjet::contrib::JetCleanser::input_nc_separate;
  fastjet::contrib::JetCleanser* cleanser = NULL ;

  if(param.exists("cleansingMethod")){
     if(param.getParameter<std::string>("cleansingMethod") == "JVF"){ // JVF method
     cleanser = new fastjet::contrib::JetCleanser(subjet_def,fastjet::contrib::JetCleanser::jvf_cleansing,projectmode);
     if(param.exists("trimCut") and param.getParameter<double>("trimCut") > 0 and  param.exists("nsubjet") and param.getParameter<int>("nsub") > 0) // additional for cuts 
       cleanser->SetTrimming(param.getParameter<double>("trimCut"));
     if(param.exists("trimCut") and param.getParameter<double>("trimCut") < 0 and  param.exists("nsubjet") and param.getParameter<int>("nsub") > 0)
       cleanser->SetFiltering(param.getParameter<double>("nsubjet"));
     if(param.exists("trimCut") and param.getParameter<double>("trimCut") > 0 and  param.exists("nsubjet") and param.getParameter<int>("nsub") > 0)
       cleanser->SetGroomingParameters(param.getParameter<double>("trimCut"),param.getParameter<double>("nsubjet"));  
    }
    else if(param.getParameter<std::string>("cleansingMethod") == "Linear"){ // Linear method
     cleanser = new fastjet::contrib::JetCleanser(subjet_def,fastjet::contrib::JetCleanser::linear_cleansing,projectmode);
     if(param.exists("linearPar")) cleanser->SetLinearParameters(param.getParameter<std::vector<double> >("cleansingParam").at(0)); // if not exsist use the fastjet default
     else cleanser->SetLinearParameters();
     if(param.exists("trimCut") and param.getParameter<double>("trimCut") > 0 and  param.exists("nsubjet") and param.getParameter<int>("nsub") > 0) // additional for cuts
       cleanser->SetTrimming(param.getParameter<double>("trimCut"));
     if(param.exists("trimCut") and param.getParameter<double>("trimCut") < 0 and  param.exists("nsubjet") and param.getParameter<int>("nsub") > 0)
       cleanser->SetFiltering(param.getParameter<double>("nsubjet"));
     if(param.exists("trimCut") and param.getParameter<double>("trimCut") > 0 and  param.exists("nsubjet") and param.getParameter<int>("nsub") > 0)
       cleanser->SetGroomingParameters(param.getParameter<double>("trimCut"),param.getParameter<double>("nsubjet"));
    }
    else if(param.getParameter<std::string>("cleansingMethod") == "Gaussian"){
     cleanser = new fastjet::contrib::JetCleanser(subjet_def,fastjet::contrib::JetCleanser::gaussian_cleansing,projectmode);
     if(param.exists("GaussianParam") and param.getParameter<std::vector<double> >("cleansingParam").size() >=4 ) // if a vdouble don't exist, used fastjet default
       cleanser->SetGaussianParameters(param.getParameter<std::vector<double> >("cleansingParam").at(0),param.getParameter<std::vector<double> >("cleansingParam").at(1),param.getParameter<std::vector<double> >("cleansingParam").at(2),param.getParameter<std::vector<double> >("cleansingParam").at(3));
     else cleanser->SetGaussianParameters();
     if(param.exists("trimCut") and param.getParameter<double>("trimCut") > 0 and  param.exists("nsubjet") and param.getParameter<int>("nsub") > 0) // additional for cuts
      cleanser->SetTrimming(param.getParameter<double>("trimCut"));
     if(param.exists("trimCut") and param.getParameter<double>("trimCut") < 0 and  param.exists("nsubjet") and param.getParameter<int>("nsub") > 0)
      cleanser->SetFiltering(param.getParameter<double>("nsubjet"));
     if(param.exists("trimCut") and param.getParameter<double>("trimCut") > 0 and  param.exists("nsubjet") and param.getParameter<int>("nsub") > 0)
      cleanser->SetGroomingParameters(param.getParameter<double>("trimCut"),param.getParameter<double>("nsubjet"));       
    }
  }
  else{ throw cms::Exception("PFPileUpAlgo --> called cleansing withou parsing parameters"); }
  

  std::vector<fastjet::PseudoJet> cleansedJets;    // jet after cleansing
  std::vector<fastjet::PseudoJet> particleToCluster;    

  if(pfCandidatesFromVtx_.empty() or pfCandidatesFromPU_.empty()) processChargedHadronSubtraction(pfCandidates,vertices) ; // if not available, redo candidate vertex association

  bool applyCHS = false;
  if(param.exists("applyCHS")){
    if(param.getParameter<bool>("applyCHS")){ // apply CHS before cleansing
      applyCHS = true ;
      particleToCluster = ConvertToPseudoJet(pfCandidatesFromVtx_);
    }    
    else particleToCluster = ConvertToPseudoJet(pfCandidates);
  }
  else particleToCluster = ConvertToPseudoJet(pfCandidates);

  fastjet::ClusterSequenceArea pCluster (particleToCluster,jet_def,area_def); // clustering on full PF collection
  std::vector<fastjet::PseudoJet> Jets = fastjet::sorted_by_pt(pCluster.inclusive_jets(param.exists("jetPtCut") ? param.getParameter<double>("jetPtCut") : 0));  

  for(std::vector<fastjet::PseudoJet>::const_iterator itJets = Jets.begin(); itJets != Jets.end(); itJets++){ // loop on jets
    std::vector<fastjet::PseudoJet> neutrals,chargedLV,jetParticles,ghosts,chargedPU;
    fastjet::SelectorIsPureGhost().sift((*itJets).constituents(), ghosts, jetParticles); // filter ghosts
    if(applyCHS) getConstitsForCleansing(jetParticles,neutrals,chargedLV,chargedPU,pfCandidatesFromVtx_,pfCandidatesFromVtx_,pfCandidatesFromPU_); // get constituents
    else getConstitsForCleansing(jetParticles,neutrals,chargedLV,chargedPU,pfCandidates,pfCandidatesFromVtx_,pfCandidatesFromPU_);
    cleansedJets.push_back((*cleanser)(neutrals,chargedLV,chargedPU)); // use cleansing                                                                                                
  }
  
  pfJetCleansing_ = ConvertJetToCandCollection(cleansedJets,pfCandidates);
  delete cleanser ;
  
} 

// method in order to run constituent subtraction
template <typename T>
void PFPileUpAlgo<T>::processConstituentSubtraction(const edm::ParameterSet& param){
  processConstituentSubtraction(pfCandidates_,*vertices_,param);

}

template <typename T>
void PFPileUpAlgo<T>::processConstituentSubtraction(const CandCollection & pfCandidates,
		  		                    const reco::VertexCollection & vertices,
                				    const edm::ParameterSet& param){


  pfConstituentSubtraction_.clear(); // clear the collection

  std::vector<fastjet::PseudoJet> eventParticles;
  
  bool applyCHS = false;
  if(param.exists("applyCHS")){
    if(param.getParameter<bool>("applyCHS")){ // if CHS has to be applied, cross check that the collection are not empty, in this case pfCandidate vertex association is not re-run
     applyCHS = true ;
     if(pfCandidatesFromVtx_.empty() or pfCandidatesFromPU_.empty()) processChargedHadronSubtraction(pfCandidates,vertices) ;
     eventParticles = ConvertToPseudoJet(pfCandidatesFromVtx_); // cluster jet on top of CHS candidates
    }
    else eventParticles = ConvertToPseudoJet(pfCandidates);
  }
  else eventParticles = ConvertToPseudoJet(pfCandidates); // use full pfCandidates

  // jet definition to evaluate rho
  fastjet::JetDefinition jet_def_for_rho (param.exists("jetAlgoRho") ? get_algo(param.getParameter<std::string>("jetAlgoRho")) : fastjet::kt_algorithm,
                                          param.exists("jetRRho") ? param.getParameter<double>("jetRRho") : jetRRho); // jet definition for bkg estimation

  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,fastjet::GhostedAreaSpec(fastjet::SelectorAbsRapMax(param.exists("GhostEtaMax") ? param.getParameter<double>("GhostEtaMax") : GhostEtaMax))); // real ghosts in the PseudoJet list     

  fastjet::Selector rho_range ; // selector for rho range as a function that you want to apply or not CHS
  if(applyCHS) rho_range = fastjet::SelectorAbsRapMax(2.5);  
  else rho_range = fastjet::SelectorAbsRapMax(param.exists("GhostEtaMax") ? param.getParameter<double>("GhostEtaMax") : GhostEtaMax);

  fastjet::JetMedianBackgroundEstimator bge_rhoC(rho_range,jet_def_for_rho, area_def); // define media background estimation
  fastjet::BackgroundJetScalarPtDensity scalarPtDensity; 
  bge_rhoC.set_jet_density_class(&scalarPtDensity); // set particles for evaluate rho
  bge_rhoC.set_particles(eventParticles);
      
  fastjet::contrib::ConstituentSubtractor *const_subtractor = 0; // constituent subtraction objec
  const_subtractor = new fastjet::contrib::ConstituentSubtractor(&bge_rhoC); // use the background
  (*const_subtractor).use_common_bge_for_rho_and_rhom(true);

  std::vector<fastjet::PseudoJet> constituentSubtractionResult ; 
  if(param.exists("clusterInJets")){
    if(param.getParameter<bool>("clusterInJets")){                  
     // define ajet algo, area and cluster sequence
     fastjet::JetDefinition jet_def(param.exists("jetAlgo") ? get_algo(param.getParameter<std::string>("jetAlgo")) : fastjet::antikt_algorithm, 
                                    param.exists("jetR") ? param.getParameter<double>("jetR") : jetR);

     fastjet::ClusterSequenceArea clusetJets (eventParticles,jet_def,area_def);
     std::vector<fastjet::PseudoJet> Jets = fastjet::sorted_by_pt(clusetJets.inclusive_jets(param.exists("jetPtCut") ? param.getParameter<double>("jetPtCut") : 0));  
     for(std::vector<fastjet::PseudoJet>::const_iterator itJet = Jets.begin(); itJet != Jets.end(); itJet++){
      constituentSubtractionResult.push_back((*const_subtractor)((*itJet)));     // apply constituent subtraction
     }
     if(applyCHS) pfConstituentSubtraction_ =  ConvertJetToCandCollection(constituentSubtractionResult,pfCandidatesFromVtx_);  
     else         pfConstituentSubtraction_ =  ConvertJetToCandCollection(constituentSubtractionResult,pfCandidates);
    }
  }
  else{
    for(std::vector<fastjet::PseudoJet>::const_iterator itParticle = eventParticles.begin(); itParticle != eventParticles.end(); itParticle++){
       constituentSubtractionResult.push_back((*const_subtractor)(*itParticle));
    }
    if(applyCHS) pfConstituentSubtraction_ =  ConvertToCandCollection(constituentSubtractionResult,pfCandidatesFromVtx_);  
    else         pfConstituentSubtraction_ =  ConvertToCandCollection(constituentSubtractionResult,pfCandidates);
  }

  delete const_subtractor;

}

// puppi method implementation
template <typename T>
void PFPileUpAlgo<T>::processPuppi(const edm::ParameterSet& param){
  processPuppi(pfCandidates_,*vertices_,param); 
}

template <typename T>
void PFPileUpAlgo<T>::processPuppi(const CandCollection & pfCandidates, // take as input pfCandidates, vertexes and parameter set
		                   const reco::VertexCollection & vertices,
		                   const edm::ParameterSet& param){

  
  pfPuppi_.clear();        // particles flow puppi particles
  fPuppiAlgo_.clear() ;    // vector of puppi algorithm objects                                                                                                         
  puppiWeights_.clear() ;  // puppi weights
  fPuppiAlgo_.clear() ;    // vector of puppi algorithm objects                                                                                                                
  puppiMetricValues_.clear() ;   // values
  lAlgos_.clear();         // algorithms
  
  // build particle collection for puppi event interpretation application
  std::vector<fastjet::PseudoJet> eventParticles,pileUpParticles,primaryVertexChargedParticles;

  if(param.exists("UseDeltaZCut")){
    if(param.getUntrackedParameter<bool>("UseDeltaZCut")){
      pfCandidates_ = pfCandidates;
      processChargedHadronSubtraction(pfCandidates,vertices,param.getUntrackedParameter<double>("DeltaZCut"));
    }
    else{
      pfCandidates_ = pfCandidates;
      if(pfCandidatesFromVtx_.empty() or pfCandidatesFromPU_.empty()) processChargedHadronSubtraction(pfCandidates,vertices) ;
    }
  }
  else{
      pfCandidates_ = pfCandidates;
      if(pfCandidatesFromVtx_.empty() or pfCandidatesFromPU_.empty()) processChargedHadronSubtraction(pfCandidates,vertices) ;
  }  

  eventParticles                = ConvertToPseudoJet(pfCandidates);  // all pfCandidates
  primaryVertexChargedParticles = ConvertToPseudoJet(pfCandidatesFromVtxCharged_); // charged from LV
  pileUpParticles               = ConvertToPseudoJet(pfCandidatesFromPU_); // charged not from LV

  bool applyCHS = false;
  if(param.exists("applyCHS")){
      if(param.getUntrackedParameter<bool>("applyCHS"))
       applyCHS = true;
  }

  double fPuppiWeightCut = 0. ;
  if(param.exists("MinPuppiWeight")){
       fPuppiWeightCut = param.getUntrackedParameter<double>("MinPuppiWeight");
  }

  // start to setup puppi algorithms objects    
  lAlgos_ = param.getParameter<std::vector<edm::ParameterSet> >("algos"); 
  for(size_t iAlgo = 0; iAlgo < lAlgos_.size(); iAlgo++) { // loop on the algorithms : central, forward1 and forward2 usually
    PuppiAlgo pPuppiConfig(lAlgos_[iAlgo]); // create the container
    fPuppiAlgo_.push_back(pPuppiConfig); //create algorithm container as a function of the algo list given as input
    fPuppiAlgo_.back().reset(); // reset
  }

  // run on the algorithms
  for(size_t iAlgo = 0; iAlgo < lAlgos_.size(); iAlgo++) 
    getPuppiRMSAvg(iAlgo,eventParticles,primaryVertexChargedParticles); // for each algo give all the particles

  // run on all the particles
  std::vector<fastjet::PseudoJet> puppiOutput ; // output particles produced by puppi
  int iParticle = 0;
  for(std::vector<fastjet::PseudoJet>::const_iterator itPart = eventParticles.begin(); itPart != eventParticles.end(); itPart++, iParticle++){  // loop on particles
   double pWeight = 1; // default weight
   int pPupId = getPuppiId((*itPart).pt(),(*itPart).eta()); // as a function of eta and pt take the id for the region looking at the algorithm --> only the first one
                                                            // if there is an overlap on the same eta and pt region
   if(pPupId == -1 or pPupId >= fPuppiAlgo_.size()) {
     puppiWeights_.push_back(1); // default value
     continue;
   }

   // fill the p-values                                                                                                                                                                 
   double pChi2   = 0;
   if(param.exists("useExp")){
     if( param.getUntrackedParameter<bool>("useExp")){
       //Compute an Experimental Puppi Weight with delta Z info (very simple example)                                                                                                     
       double dZ = pfCandidates_.at((*itPart).user_index())->vertex().z()-vertices.at(0).z();
       pChi2 = getChi2FromdZ(dZ);					   
       //Now make sure Neutrals are not set                                                                                                                                               
       if(pfCandidates_.at((*itPart).user_index())->charge() == 0) pChi2 = 0;
     }
   }

   pWeight = fPuppiAlgo_[pPupId].compute(puppiMetricValues_[iParticle],pChi2);   // compute the weight for the given puppiId 

   //Apply the CHS weights --> find the particles in the CHS collection  --> this selection can be done at priori
   if(particleIndex_.at(iParticle) == 1 && applyCHS ) pWeight = 1; // charged from LV
   if(particleIndex_.at(iParticle) == 2 && applyCHS ) pWeight = 0; // charged from PU   
   if(std::isnan(pWeight)) std::cerr << "====> Weight is nan : pt " << (*itPart).pt() << " -- eta : " << (*itPart).eta() << " -- Value" << puppiMetricValues_[iParticle] << std::endl;

   //Basic Cuts                                                                                                                                                                           
   if(pWeight                         < fPuppiWeightCut) pWeight = 0;  //==> Elminate the low Weight stuff                                                                                
   if(pWeight*(*itPart).pt() < fPuppiAlgo_[pPupId].getNeutralPt(int(vertices.size())) && particleIndex_.at(iParticle) == 0 ) pWeight = 0;  //threshold cut on the neutral Pt          

   puppiWeights_.push_back(pWeight); // fill the weight

   //Now get rid of the thrown out weights for the particle collection                                                                                                                    
   if(pWeight == 0) continue;
   //Produce                                                                                                                                                                            
   fastjet::PseudoJet curjet( pWeight*(*itPart).px(), pWeight*(*itPart).py(), pWeight*(*itPart).pz(), pWeight*(*itPart).e());
   curjet.set_user_index(iParticle);                                                                                                                              
   puppiOutput.push_back(curjet);   
  }

  pfPuppi_ =  ConvertToCandCollection(puppiOutput,pfCandidates_);  

}

// Method to identify the vertex of a single pfCand
template <typename T>
int PFPileUpAlgo<T>::chargedHadronVertex( const reco::VertexCollection& vertices, const reco::PFCandidate & pfcand) const {
  
  
  if(pfcand.particleId() == reco::PFCandidate::h) {
     auto const & track = pfcand.trackRef(); ;
     size_t  iVertex = 0;
     unsigned int index = 0;
     unsigned int nFoundVertex = 0;
     float bestweight = 0;
     for( auto const & vtx : vertices) { // loop on vertex collection
      float w = vtx.trackWeight(track); //select the vertex for which the track has the highest weight
 	if (w > bestweight){
	  bestweight = w;
	  iVertex = index;
	  nFoundVertex++;
	}
      ++index;
     }

     if (nFoundVertex>0){ // if a vertex is found
      if (nFoundVertex!=1) // mote than one vertex
      edm::LogWarning("TrackOnTwoVertex")<<"a track is shared by at least two verteces. Used to be an assert";
      return iVertex; // return vertex
     }

     // no vertex found with this track. 
     if ( checkClosestZVertex_ ) {
      double dzmin = 10000;
      double ztrack = pfcand.vertex().z();
      bool foundVertex = false;
      index = 0;
      for(auto iv=vertices.begin(); iv!=vertices.end(); ++iv, ++index) {
       double dz = fabs(ztrack - iv->z());
       if(dz < dzmin) { // check closer in z and assigne the particle to it
	dzmin = dz; 
	iVertex = index;
	foundVertex = true;
       }
      }
      if( foundVertex  ) return iVertex;  
     }
     return -1 ;
  }  
  else return -1 ;   

}

// Method to identify the vertex of a single pfCand
template <typename T>
int PFPileUpAlgo<T>::chargedHadronVertex( const reco::VertexCollection& vertices, const pat::PackedCandidate & pfcand) const {

  if(&(*(pfcand.vertexRef())) == 0) return -1 ;  
   
  int pVtxId = pfcand.fromPV();
  if((pVtxId == 3 or pVtxId == 2) and pfcand.charge() == 0) return -1;
  if((pVtxId == 3 or pVtxId == 2) and pfcand.charge() != 0) return 0;
  if((pVtxId == 1 or pVtxId == 0) and pfcand.charge() != 0) return 1;
  if((pVtxId == 1 or pVtxId == 0) and pfcand.charge() == 0) return -1;

  // no vertex found with this track. 
  unsigned int index = 0;
  size_t  iVertex = 0;
  if ( checkClosestZVertex_ ) {
    double dzmin = 10000;
    double ztrack = pfcand.vertex().z();
    bool foundVertex = false;
    index = 0;
    for(auto iv=vertices.begin(); iv!=vertices.end(); ++iv, ++index) {
      double dz = fabs(ztrack - iv->z());
      if(dz < dzmin) { // check closer in z and assigne the particle to it
	dzmin = dz; 
	iVertex = index;
	foundVertex = true;
      }
    }
    if( foundVertex  ) return iVertex;  
  }
  return -1 ;
    
}


template <typename T>
fastjet::JetAlgorithm PFPileUpAlgo<T>::get_algo(const std::string & algo){ // method  to define a jet algorithm in fast jet from a generic string
  fastjet::JetAlgorithm jetalgo;
  if (algo == "kt")      jetalgo = fastjet::kt_algorithm ;
  else if (algo == "ca") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo == "ak") jetalgo = fastjet::antikt_algorithm ;
  else if (algo == "KT") jetalgo = fastjet::kt_algorithm ;
  else if (algo == "CA") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo == "AK") jetalgo = fastjet::antikt_algorithm ;
  else if (algo == "kt_algorithm") jetalgo = fastjet::kt_algorithm ;
  else if (algo == "cambridge_algorithm") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo == "antikt_algorithm") jetalgo = fastjet::antikt_algorithm ;
  else if (algo == "0") jetalgo = fastjet::kt_algorithm ;
  else if (algo == "1") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo == "2") jetalgo = fastjet::antikt_algorithm ;
  else jetalgo = fastjet::antikt_algorithm ;
  return jetalgo;
}

// convert pf candidates to vector of  pseudojets
template <typename T>
std::vector<fastjet::PseudoJet> PFPileUpAlgo<T>::ConvertToPseudoJet(const CandCollection & pfCandidates){ // take a CandCollection and convert it in a list of pseudojets
  std::vector<fastjet::PseudoJet> jetParticles ;
  int ipfCandidate = 0;
  for(typename CandCollection::const_iterator  itPF = pfCandidates.begin(); itPF != pfCandidates.end(); itPF++){
    fastjet::PseudoJet jetParticle ((*itPF)->px(),(*itPF)->py(),(*itPF)->pz(),(*itPF)->energy());
    jetParticle.set_user_index(ipfCandidate);
    jetParticles.push_back(jetParticle);
    ipfCandidate++;
  }
  return jetParticles;
}

template <typename T>
std::vector<edm::FwdPtr<T> > PFPileUpAlgo<T>::ConvertToCandCollection(const std::vector<fastjet::PseudoJet> & inputParticles, const CandCollection & pfParticles){
  CandCollection pfAssociatedParticles ;
  for(std::vector<fastjet::PseudoJet>::const_iterator itPart = inputParticles.begin(); itPart != inputParticles.end(); itPart++){
      if((*itPart).user_index() < 0) continue;
      pfAssociatedParticles.push_back(pfParticles.at((*itPart).user_index()));
  }
  return pfAssociatedParticles;
}

template <typename T>
std::vector<edm::FwdPtr<T> > PFPileUpAlgo<T>::ConvertJetToCandCollection(const std::vector<fastjet::PseudoJet> & inputJets, const CandCollection & pfParticles){
  CandCollection pfAssociatedParticles ;
  for(std::vector<fastjet::PseudoJet>::const_iterator itJets = inputJets.begin(); itJets != inputJets.end(); itJets++){
    std::vector<fastjet::PseudoJet> jetParticles;
    std::vector<fastjet::PseudoJet> ghosts;
    fastjet::SelectorIsPureGhost().sift((*itJets).constituents(), ghosts, jetParticles); // filter ghosts                
    for(std::vector<fastjet::PseudoJet>::const_iterator  itPart = jetParticles.begin(); itPart != jetParticles.end(); itPart++){
      if((*itPart).user_index() < 0) continue;
      pfAssociatedParticles.push_back(pfParticles.at((*itPart).user_index()));
    }
  }
  return pfAssociatedParticles;
}

template <typename T>
void PFPileUpAlgo<T>::getConstitsForCleansing(const std::vector<fastjet::PseudoJet> & jetParticles, std::vector<fastjet::PseudoJet> & neutrals,
                                              std::vector<fastjet::PseudoJet> & chargedLV, std::vector<fastjet::PseudoJet> & chargedPU,
                                              const CandCollection & eventParticles, const CandCollection & vertexParticles, const CandCollection & pileupParticles){

  for(std::vector<fastjet::PseudoJet>::const_iterator itJetConstituent = jetParticles.begin(); itJetConstituent!=jetParticles.end(); ++itJetConstituent){
     int user_index = (*itJetConstituent).user_index(); // particle  user index found in the full PFList
     if(user_index < 0) continue;
     typename CandCollection::const_iterator itCand = pileupParticles.begin();             
     for(; itCand!=pileupParticles.end(); itCand++){
       if( (*itCand)  == eventParticles.at(user_index)) break ;
     }

     if(itCand != pileupParticles.end()){	 
	if((*itCand)->charge() == 0) throw cms::Exception("PFPileUpAlgo --> called cleansing constituents problem with CHS"); 
        else chargedPU.push_back((*itJetConstituent));    
     }
     itCand = vertexParticles.begin();	      
     for(; itCand!=vertexParticles.end(); itCand++){
       if((*itCand) == eventParticles.at(user_index) ) break ;
     }    
     if(itCand != vertexParticles.end()){
      if((*itCand)->charge() == 0) neutrals.push_back((*itJetConstituent));
      else chargedLV.push_back((*itJetConstituent));    
    }
  }
     
  return ;
    
}

template <typename T>
void PFPileUpAlgo<T>::getPuppiRMSAvg(const int & iOpt, std::vector<fastjet::PseudoJet> & eventParticles, std::vector<fastjet::PseudoJet> & chargedParticlesPV) {
  for(size_t iConst = 0; iConst < eventParticles.size(); iConst++ ) { // Loop on the particles
    //double pVal = -1;
    //Calculate the Puppi Algo to use                                                                                                                                             
    int  pPupId = getPuppiId(eventParticles[iConst].pt(),eventParticles[iConst].eta());
    if(pPupId > fPuppiAlgo_.size()) pPupId = -1; // take the algorithm that can be different asaf of eta and pt
    if(pPupId != iOpt) continue;
    if(pPupId == -1) {
        puppiMetricValues_.push_back(-1); 
        continue;
    }
    
    
    //Get the Puppi Sub Algo (given iteration)                                                                                                                                            
    int  pAlgo    = fPuppiAlgo_[pPupId].getAlgoId   ();
    bool pCharged = fPuppiAlgo_[pPupId].getIsCharged();
    double pCone  = fPuppiAlgo_[pPupId].getConeSize ();
    
    //Compute the Puppi Metric
    double pVal = -1 ;                                                                                                                                                            
    if(!pCharged) pVal = goodVarPuppi(eventParticles[iConst],eventParticles,pAlgo,pCone);
    if( pCharged) pVal = goodVarPuppi(eventParticles[iConst],chargedParticlesPV,pAlgo,pCone);
    puppiMetricValues_.push_back(pVal); // fill puppi values

    if(std::isnan(pVal) || std::isinf(pVal)) std::cerr << "====> Value is Nan " << pVal << " == " << eventParticles[iConst].pt() << " -- " << eventParticles[iConst].eta() << std::endl;
    if(std::isnan(pVal) || std::isinf(pVal)) continue;
    bool isFromPV = true ;
    typename CandCollection::const_iterator itCand = pfCandidatesFromPU_.begin();             
    for(; itCand!=pfCandidatesFromPU_.end(); itCand++){
      if( (*itCand)  == pfCandidates_[eventParticles[iConst].user_index()]) break ;
    }
    if(itCand != pfCandidatesFromPU_.end()){	 
	if((*itCand)->charge() == 0) throw cms::Exception("PFPileUpAlgo --> called puppi constituents problem with CHS"); 
        isFromPV = false ;
    }    
    fPuppiAlgo_[pPupId].add(eventParticles[iConst],pfCandidates_[eventParticles[iConst].user_index()],isFromPV,pVal);
  }

  for(size_t iAlgo = 0; iAlgo < fPuppiAlgo_.size(); iAlgo++) 
    fPuppiAlgo_[iAlgo].computeMedRMS(pfCandidatesFromVtxCharged_.size()/(pfCandidatesFromVtxCharged_.size()+pfCandidatesFromPU_.size()));
}

template <typename T>
int PFPileUpAlgo<T>::getPuppiId(const float & pt,const float & eta) {
  int lId = -1;
  for(size_t iAlgo = 0; iAlgo < lAlgos_.size(); iAlgo++) {
    if(fabs(eta) < fPuppiAlgo_[iAlgo].getEtaMin()) continue;
    if(fabs(eta) > fPuppiAlgo_[iAlgo].getEtaMax()) continue;
    if(pt < fPuppiAlgo_[iAlgo].getPtMin()) continue;
    lId = iAlgo;
    break;
  }
  return lId;
}

template <typename T>
double PFPileUpAlgo<T>::goodVarPuppi(const fastjet::PseudoJet & iPart, std::vector<fastjet::PseudoJet> & iParts, int iOpt, double iRCone) {
  double lPup = 0;
  if(iOpt == -1) return 1;
  fastjet::Selector coneSelector = fastjet::SelectorCircle(iRCone); // define a cone region around each particle
  coneSelector.set_reference(iPart); // give the center of  the code
  std::vector<fastjet::PseudoJet> near_particles = coneSelector(iParts); // select the event particles
  for(size_t iParticles = 0; iParticles < near_particles.size(); iParticles++){ // loop on the selected particles
    double pDEta = near_particles[iParticles].eta()-iPart.eta(); // dEta
    double pDPhi = fabs(near_particles[iParticles].phi()-iPart.phi()); // dPhi
    if(pDPhi > 2.*3.14159265-pDPhi) pDPhi = 2.*3.14159265-pDPhi;
    double pDR = sqrt(pDEta*pDEta+pDPhi*pDPhi);
    if(pDR < 0.001) continue; // discard if to close with center
    if(pDR < 0.01)  continue;//pDR = 0.01;
    if(pDR == 0)    continue;
    if(iOpt == 0) lPup += (near_particles[iParticles].pt()/pDR/pDR); // metric = pt/DR^2
    if(iOpt == 1) lPup += near_particles[iParticles].pt(); // metric = pt
    if(iOpt == 2) lPup += (1./pDR)*(1./pDR);  // metric = 1/DR^2
    if(iOpt == 3) lPup += (1./pDR)*(1./pDR); // metric = 1/DR^2
    if(iOpt == 4) lPup += near_particles[iParticles].pt(); // metric = pt
    if(iOpt == 5) lPup += (near_particles[iParticles].pt()/pDR)*(near_particles[iParticles].pt()/pDR); // metric = pt^2/DR^2
  }
  if(iOpt == 0 && lPup != 0) lPup = log(lPup); // take the log
  if(iOpt == 3 && lPup != 0) lPup = log(lPup);
  if(iOpt == 5 && lPup != 0) lPup = log(lPup);
  return lPup;
}

template <typename T>
double PFPileUpAlgo<T>::getChi2FromdZ(const double & dZ) {
  //We need to obtain prob of PU + (1-Prob of LV)
  // Prob(LV) = Gaus(dZ,sigma) where sigma = 1.5mm (its really more like 1mm)
  //double lProbLV = ROOT::Math::normal_cdf_c(fabs(iDZ),0.2)*2.; //*2 is to do it double sided
  //Take iDZ to be corrected by sigma already
  double lProbLV = ROOT::Math::normal_cdf_c(fabs(dZ),1.)*2.; //*2 is to do it double sided
  double lProbPU = 1-lProbLV;
  if(lProbPU <= 0) lProbPU = 1e-16; //Quick Trick to through out infs
  if(lProbPU >= 0) lProbPU = 1-1e-16; //Ditto
  double lChi2PU = TMath::ChisquareQuantile(lProbPU,1);
  lChi2PU*=lChi2PU;
  return lChi2PU;
}
