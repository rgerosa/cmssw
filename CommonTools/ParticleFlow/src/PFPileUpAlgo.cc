#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// alternative constructor
PFPileUpAlgo::PFPileUpAlgo(const PFCollection & pfCandidates, const reco::VertexCollection & vertices, bool checkClosestZVertex, bool verbose){

  checkClosestZVertex_ = checkClosestZVertex;
  verbose_             = verbose;
  vertices_            = &vertices ;
  pfCandidates_        = pfCandidates;
  processChargedHadronSubtraction(pfCandidates_,*vertices_); // process charged hadron subtraction by default

}

PFPileUpAlgo::~PFPileUpAlgo(){

  pfCandidatesFromVtx_.clear();
  pfCandidatesFromPU_.clear();
  pfCandidates_.clear();
  pfSoftKiller_.clear() ;                                                                                                                            
  pfPuppi_.clear();                                                                                                                                           
  pfJetCleansing_.clear() ;                                                                                                                      
  pfConstituentSubtraction_.clear() ;
  delete vertices_ ;

}

void PFPileUpAlgo::processChargedHadronSubtraction(const PFCollection & pfCandidates, 
			                           const reco::VertexCollection & vertices){

  pfCandidatesFromVtx_.clear(); // clear existing collections
  pfCandidatesFromPU_.clear();
  for( unsigned i = 0; i < pfCandidates.size(); i++ ) { // loop on PFcandidates    
    int ivertex;
    switch( pfCandidates.at(i)->particleId() ) {
    case reco::PFCandidate::h:
      ivertex = chargedHadronVertex(vertices, *(pfCandidates.at(i))); // identify the vertex
      break;
    default:
      pfCandidatesFromVtx_.push_back( pfCandidates[i] ); // also leptons, neutrals are taken into account as Vtx particles
      continue;
    }     
    // no associated vertex, or primary vertex  --> put in the fromVtx list
    if( ivertex == -1  || ivertex == 0 ) {
      if(verbose_) std::cout<<"VTX "<<i<<" "<< *(pfCandidates[i])<<std::endl;
      pfCandidatesFromVtx_.push_back( pfCandidates[i] );
    } 
    else { // associated to another vertex
      if(verbose_) std::cout<<"PU  "<<i<<" "<< *(pfCandidates[i])<<std::endl;
      pfCandidatesFromPU_.push_back( pfCandidates[i] );
    }        
  }
}

// method in prder to run soft killer
void PFPileUpAlgo::processSoftKiller(const edm::ParameterSet& param){

  processSoftKiller(pfCandidates_,*vertices_,param);

}

void PFPileUpAlgo::processSoftKiller(const PFCollection & pfCandidates,
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
     
     if(applyCHS) pfSoftKiller_ = ConvertJetToPFCollection(softJets,pfCandidatesFromVtx_); // give jets as input of converting function and look at Vtx particles
     else         pfSoftKiller_ = ConvertJetToPFCollection(softJets,pfCandidates); // give jets as input of converting function

    }
    else{

      if(applyCHS) pfSoftKiller_ = ConvertToPFCollection(softParticles,pfCandidatesFromVtx_);
      else         pfSoftKiller_ = ConvertToPFCollection(softParticles,pfCandidates); // give particles as input to the converting function
    }    
  }
  else{ if(applyCHS) pfSoftKiller_ = ConvertToPFCollection(softParticles,pfCandidatesFromVtx_);
        else         pfSoftKiller_ = ConvertToPFCollection(softParticles,pfCandidates); // give particles as input to the converting function
  }
}
  

// implement jet cleaning
void PFPileUpAlgo::processJetCleansing(const edm::ParameterSet& param){

  processJetCleansing(pfCandidates_,*vertices_,param);

}

void PFPileUpAlgo::processJetCleansing(const PFCollection & pfCandidates,
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
  
  pfJetCleansing_ = ConvertJetToPFCollection(cleansedJets,pfCandidates);
  delete cleanser ;
  
} 

// method in order to run constituent subtraction

void PFPileUpAlgo::processConstituentSubtraction(const edm::ParameterSet& param){
  processConstituentSubtraction(pfCandidates_,*vertices_,param);

}

void PFPileUpAlgo::processConstituentSubtraction(const PFCollection & pfCandidates,
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
     if(applyCHS) pfConstituentSubtraction_ =  ConvertJetToPFCollection(constituentSubtractionResult,pfCandidatesFromVtx_);  
     else         pfConstituentSubtraction_ =  ConvertJetToPFCollection(constituentSubtractionResult,pfCandidates);
    }
  }
  else{
    for(std::vector<fastjet::PseudoJet>::const_iterator itParticle = eventParticles.begin(); itParticle != eventParticles.end(); itParticle++){
       constituentSubtractionResult.push_back((*const_subtractor)(*itParticle));
       if(applyCHS) pfConstituentSubtraction_ =  ConvertToPFCollection(constituentSubtractionResult,pfCandidatesFromVtx_);  
       else         pfConstituentSubtraction_ =  ConvertToPFCollection(constituentSubtractionResult,pfCandidates);
    }
  }

  delete const_subtractor;

}

// puppi method implementation
void PFPileUpAlgo::processPuppi(const edm::ParameterSet& param){
  processPuppi(pfCandidates_,*vertices_,param); 
}

void PFPileUpAlgo::processPuppi(const PFCollection & pfCandidates,
		  const reco::VertexCollection & vertices,
		  const edm::ParameterSet& param){

  std::cout<<" Dummy code "<<std::endl;
  
}

// method to identify the vertex of a single pfCand
int PFPileUpAlgo::chargedHadronVertex( const reco::VertexCollection& vertices, const reco::PFCandidate& pfcand ) const {

  auto const & track = pfcand.trackRef();  // tacke the track    
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
    if( foundVertex ) return iVertex;  
  }

  return -1 ;
}


fastjet::JetAlgorithm PFPileUpAlgo::get_algo(const std::string & algo){ // method  to define a jet algorithm in fast jet from a generic string
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
std::vector<fastjet::PseudoJet> PFPileUpAlgo::ConvertToPseudoJet(const PFCollection & pfCandidates){ // take a PFCollection and convert it in a list of pseudojets
  std::vector<fastjet::PseudoJet> jetParticles ;
  int ipfCandidate = 0;
  for(PFCollection::const_iterator  itPF = pfCandidates.begin(); itPF != pfCandidates.end(); itPF++){
    fastjet::PseudoJet jetParticle ((*itPF)->px(),(*itPF)->py(),(*itPF)->pz(),(*itPF)->energy());
    jetParticle.set_user_index(ipfCandidate);
    jetParticles.push_back(jetParticle);
    ipfCandidate++;
  }
  return jetParticles;
}

PFCollection PFPileUpAlgo::ConvertToPFCollection(const std::vector<fastjet::PseudoJet> & inputParticles, const PFCollection & pfParticles){
  PFCollection pfAssociatedParticles ;
  for(std::vector<fastjet::PseudoJet>::const_iterator itPart = inputParticles.begin(); itPart != inputParticles.end(); itPart++){
      if((*itPart).user_index() < 0) continue;
      pfAssociatedParticles.push_back(pfParticles.at((*itPart).user_index()));
  }
  return pfAssociatedParticles;
}

PFCollection PFPileUpAlgo::ConvertJetToPFCollection(const std::vector<fastjet::PseudoJet> & inputJets, const PFCollection & pfParticles){
  PFCollection pfAssociatedParticles ;
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


void PFPileUpAlgo::getConstitsForCleansing(const std::vector<fastjet::PseudoJet> & jetParticles, std::vector<fastjet::PseudoJet> & neutrals,
                                           std::vector<fastjet::PseudoJet> & chargedLV, std::vector<fastjet::PseudoJet> & chargedPU,
                                           const PFCollection & eventParticles, const PFCollection & vertexParticles, const PFCollection & pileupParticles){

  for(std::vector<fastjet::PseudoJet>::const_iterator itJetConstituent = jetParticles.begin(); itJetConstituent!=jetParticles.end(); ++itJetConstituent){
     int user_index = (*itJetConstituent).user_index(); // particle  user index found in the full PFList
     if(user_index < 0) continue;
     PFCollection::const_iterator itCand = vertexParticles.begin();             
     for(; itCand!=vertexParticles.end(); itCand++){
       if((*itCand) == eventParticles.at(user_index) ) break ;
     }    
     if(itCand != vertexParticles.end()){
      if((*itCand)->charge() == 0) neutrals.push_back((*itJetConstituent));
      else chargedLV.push_back((*itJetConstituent));    
    }
    else{
       itCand = pileupParticles.begin();
       for(; itCand!=pileupParticles.end(); itCand++){
	 if( (*itCand)  == eventParticles.at(user_index)) break ;
       }       
       if(itCand != pileupParticles.end()){	 
	if((*itCand)->charge() == 0) throw cms::Exception("PFPileUpAlgo --> called cleansing constituents problem with CHS"); 
        else chargedPU.push_back((*itJetConstituent));    
      }
    }	      
  }
  return ;
    
}
