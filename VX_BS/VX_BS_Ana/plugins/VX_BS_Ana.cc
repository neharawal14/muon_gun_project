// -*- C++ -*-
//
// Package:    VX_BS/VX_BS_Ana
// Class:      VX_BS_Ana
//
/**\class VX_BS_Ana VX_BS_Ana.cc VX_BS/VX_BS_Ana/plugins/VX_BS_Ana.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Filippo Errico
//         Created:  Fri.clear(); 31 Jul 2020 10:45:38 GMT
//
//

#include <fstream>
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include <vector>

#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "RecoVertex/KalmanVertexFit/interface/SingleTrackVertexConstraint.h"

#include "TMath.h"
#include "TLorentzVector.h"

#include "TTree.h"

#include <tuple>
//
// class declaration
//

// If the analyzer does not use TFileService.clear(); please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


// using reco::Muon;

     

class VX_BS_Ana : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit VX_BS_Ana(const edm::ParameterSet&);
      ~VX_BS_Ana();
      int number_events;
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      void bookPassedEventTree(TString treeName, TTree *tree);      
      TTree *passedEventsTree;
      
      float BS_x, BS_y, BS_z; 
      float BS_xErr, BS_yErr, BS_zErr; 
      float BeamWidth_x, BeamWidth_y;
      float BeamWidth_xErr, BeamWidth_yErr;

      float PV_x, PV_y, PV_z; 
      float PV_xErr, PV_yErr, PV_zErr; 
      float PV_Width_x, PV_Width_y;
      float PV_Width_xErr, PV_Width_yErr;
      int nVtx; 

      float xerr_bs, yerr_bs, xerr_bs_set, yerr_bs_set;
      int number_flo_muons =0 , number_tune_muons =0 , number_reco_muons =0, glb_count =0 ,  flo_count =0, standAlone_count=0, gen_count =0;
      

      std::vector<double> pt_GEN; std::vector<double> eta_GEN; std::vector<double> phi_GEN; std::vector<double> mass_GEN; std::vector<double> charge_GEN; 
      std::vector<double> pt_trk; std::vector<double> eta_trk; std::vector<double> phi_trk; std::vector<double> mass_trk; 
      std::vector<double> dxy_BS_trk; std::vector<double> charge_trk; std::vector<double> d0_trk;
      
 // for innertrack with BS    
      std::vector<double> pt_trk_BS; std::vector<double> eta_trk_BS; std::vector<double> phi_trk_BS; std::vector<double> mass_trk_BS; 
      std::vector<double> dxy_trk_BS_const; 
      std::vector<double> dxy_trk_PV_const;
      std::vector<double> pt_trk_BS_float; std::vector<double> eta_trk_BS_float; std::vector<double> phi_trk_BS_float; std::vector<double> mass_trk_BS_float; 
      std::vector<double> dxy_trk_BS_const_float; 
      std::vector<double> dxy_trk_PV_const_float;

           std::vector<double> pt_trk_BS_1; std::vector<double> eta_trk_BS_1; std::vector<double> phi_trk_BS_1; std::vector<double> mass_trk_BS_1; 
      std::vector<double> dxy_trk_BS_const_1; 
      std::vector<double> dxy_trk_PV_const_1;
      std::vector<double> pt_trk_BS_1_float; std::vector<double> eta_trk_BS_1_float; std::vector<double> phi_trk_BS_1_float; std::vector<double> mass_trk_BS_1_float; 
      std::vector<double> dxy_trk_BS_const_1_float; 
      std::vector<double> dxy_trk_PV_const_1_float;

     
     
      std::vector<double> d0_PV_trk; std::vector<double> d0_BS_trk; std::vector<double> dxy_PV_trk;
      
      std::vector<float> pt_trk_float; std::vector<float> eta_trk_float; std::vector<float> phi_trk_float; std::vector<float> mass_trk_float;
      std::vector<float> dxy_BS_trk_float; std::vector<double> charge_trk_float; std::vector<double> d0_trk_float;
      std::vector<double> d0_PV_trk_float; std::vector<double> d0_BS_trk_float; std::vector<double> dxy_PV_trk_float;
 

      std::vector<double> pt_FLO; std::vector<double> eta_FLO; std::vector<double> phi_FLO; std::vector<double> mass_FLO; std::vector<double> dxy_BS_FLO; std::vector<double> dxy_PV_FLO; 
      std::vector<double> pt_FLO_vtx; std::vector<double> eta_FLO_vtx; std::vector<double> phi_FLO_vtx; std::vector<double> mass_FLO_vtx; std::vector<double> dxy_BS_FLO_vtx;
      std::vector<double> dxy_PV_FLO_vtx;  std::vector<double> dxy_PV_FLO_vtx_float; 
      std::vector<double> charge_FLO; std::vector<double> charge_FLO_float;
      std::vector<double> pt_FLO_single; std::vector<double> eta_FLO_single; std::vector<double> phi_FLO_single; std::vector<double> dxy_BS_FLO_single; 
      std::vector<double> pt_glb; std::vector<double> eta_glb; std::vector<double> phi_glb; std::vector<double> mass_glb; std::vector<double> dxy_BS_glb; std::vector<double> dxy_PV_glb; std::vector<double> charge_glb;
      std::vector<double> pt_glb_vtx; std::vector<double> eta_glb_vtx; std::vector<double> phi_glb_vtx; std::vector<double> mass_glb_vtx; std::vector<double> dxy_BS_glb_vtx; 
      std::vector<double> pt_glb_single; std::vector<double> eta_glb_single; std::vector<double> phi_glb_single; std::vector<double> dxy_BS_glb_single; 
      std::vector<double> pt_stdAlone; std::vector<double> eta_stdAlone; std::vector<double> phi_stdAlone; std::vector<double> dxy_BS_stdAlone; std::vector<double> dxy_PV_stdAlone; std::vector<double> charge_stdAlone; std::vector<double> mass_stdAlone;
      std::vector<double> pt_stdAlone_vtx; std::vector<double> eta_stdAlone_vtx; std::vector<double> phi_stdAlone_vtx; std::vector<double> dxy_BS_stdAlone_vtx; std::vector<double> dxy_PV_stdAlone_vtx; 

      std::vector<float> pt_GEN_float; std::vector<float> eta_GEN_float; std::vector<float> phi_GEN_float; std::vector<float> mass_GEN_float; std::vector<double> charge_GEN_float;
      std::vector<float> pt_FLO_float; std::vector<float> eta_FLO_float; std::vector<float> phi_FLO_float; std::vector<float> mass_FLO_float; std::vector<float> dxy_BS_FLO_float; std::vector<float> dxy_PV_FLO_float; 
      std::vector<float> pt_FLO_vtx_float; std::vector<float> eta_FLO_vtx_float; std::vector<float> phi_FLO_vtx_float; std::vector<float> mass_FLO_vtx_float; std::vector<float> dxy_BS_FLO_vtx_float; 
      std::vector<float> pt_FLO_single_float; std::vector<float> eta_FLO_single_float; std::vector<float> phi_FLO_single_float; std::vector<float> dxy_BS_FLO_single_float; 
      std::vector<float> pt_glb_float; std::vector<float> eta_glb_float; std::vector<float> phi_glb_float; std::vector<float> mass_glb_float; std::vector<float> dxy_BS_glb_float; std::vector<float> dxy_PV_glb_float; std::vector<double> charge_glb_float; 
      std::vector<float> pt_glb_vtx_float; std::vector<float> eta_glb_vtx_float; std::vector<float> phi_glb_vtx_float; std::vector<float> mass_glb_vtx_float; std::vector<float> dxy_BS_glb_vtx_float; 
      std::vector<float> pt_glb_single_float; std::vector<float> eta_glb_single_float; std::vector<float> phi_glb_single_float; std::vector<float> dxy_BS_glb_single_float; 
      std::vector<float> pt_stdAlone_float; std::vector<float> eta_stdAlone_float; std::vector<float> phi_stdAlone_float; std::vector<float> dxy_BS_stdAlone_float; std::vector<float> dxy_PV_stdAlone_float; std::vector<double> charge_stdAlone_float; std::vector<double> mass_stdAlone_float;
      std::vector<float> pt_stdAlone_vtx_float; std::vector<float> eta_stdAlone_vtx_float; std::vector<float> phi_stdAlone_vtx_float; std::vector<float> dxy_BS_stdAlone_vtx_float; std::vector<float> dxy_PV_stdAlone_vtx_float; 

      std::vector<float> dxy_BS_FLO_vtx_second_test_float; std::vector<float> dxy_BS_FLO_vtx_second_test;
      
// for standAlone muons in muon container
       std::vector<double> pt_stdAlone_trk; std::vector<double> eta_stdAlone_trk; std::vector<double> phi_stdAlone_trk; std::vector<double> dxy_BS_stdAlone_trk; std::vector<double> dxy_PV_stdAlone_trk; std::vector<double> charge_stdAlone_trk; std::vector<double> mass_stdAlone_trk;
	     std::vector<float> pt_stdAlone_trk_float; std::vector<float> eta_stdAlone_trk_float; std::vector<float> phi_stdAlone_trk_float; std::vector<float> dxy_BS_stdAlone_trk_float; std::vector<float> dxy_PV_stdAlone_trk_float; std::vector<double> charge_stdAlone_trk_float; std::vector<double> mass_stdAlone_trk_float;

// not tracker muon but global and standAlone muon
       std::vector<double> pt_glb_not_trk; std::vector<double> pt_stdAlone_not_trk;
       std::vector<double> pt_glb_not_trk_float; std::vector<double> pt_stdAlone_not_trk_float;

// for global muons in muon container
        std::vector<double> pt_glb_trk; std::vector<double> eta_glb_trk; std::vector<double> phi_glb_trk; std::vector<double> dxy_BS_glb_trk; std::vector<double> dxy_PV_glb_trk; std::vector<double> charge_glb_trk; std::vector<double> mass_glb_trk;
        std::vector<float> pt_glb_trk_float; std::vector<float> eta_glb_trk_float; std::vector<float> phi_glb_trk_float; std::vector<float> dxy_BS_glb_trk_float; std::vector<float> dxy_PV_glb_trk_float; std::vector<double> charge_glb_trk_float; std::vector<double> mass_glb_trk_float;

        std::vector<double> pt_tuned_float; std::vector<double> eta_tuned_float; std::vector<double> phi_tuned_float; std::vector<double> mass_tuned_float;
        std::vector<double> dxy_BS_tuned_float; std::vector<double> dxy_PV_tuned_float; std::vector<double> charge_tuned_float; std::vector<double> d0_tuned_float;
    
        std::vector<double> pt_tuned; std::vector<double> eta_tuned; std::vector<double> phi_tuned; std::vector<double> mass_tuned;
        std::vector<double> dxy_BS_tuned; std::vector<double> dxy_PV_tuned; std::vector<double> charge_tuned; std::vector<double> d0_tuned;

//variables for beamspot 
       std::vector<float> beam_x; std::vector<float> beam_y; std::vector<float> beam_z; std::vector<float> beam_width_x; std::vector<float> beam_width_y;
       std::vector<float> beam_x_err; std::vector<float> beam_y_err; std::vector<float> beam_z_err; std::vector<float> beam_width_x_err; std::vector<float> beam_width_y_err;

       std::vector<float> beam_x_float; std::vector<float> beam_y_float; std::vector<float> beam_z_float; std::vector<float> beam_width_x_float; std::vector<float> beam_width_y_float;
       std::vector<float> beam_x_err_float; std::vector<float> beam_y_err_float; std::vector<float> beam_z_err_float; std::vector<float> beam_width_x_err_float; std::vector<float> beam_width_y_err_float;

// primary vertex variables
       std::vector<float> pv_x; std::vector<float> pv_y; std::vector<float> pv_z; std::vector<float> pv_width_x; std::vector<float> pv_width_y;
       std::vector<float> pv_x_err; std::vector<float> pv_y_err; std::vector<float> pv_z_err; std::vector<float> pv_width_x_err; std::vector<float> pv_width_y_err;

       std::vector<float> pv_x_float; std::vector<float> pv_y_float; std::vector<float> pv_z_float; std::vector<float> pv_width_x_float; std::vector<float> pv_width_y_float;
       std::vector<float> pv_x_err_float; std::vector<float> pv_y_err_float; std::vector<float> pv_z_err_float; std::vector<float> pv_width_x_err_float; std::vector<float> pv_width_y_err_float;

 // parameters for vx,vy,vz and 1 for dxy also

       std::vector<double> v_x, v_y, v_z; 
       std::vector<double> dxy_thebeamspot_one, dxy_primaryvertex_one; std::vector<double> dxy_thebeamspot_second;
       std::vector<double> dxy_without_argument;

       std::vector<double> v_x_float, v_y_float, v_z_float; 
       std::vector<double> dxy_thebeamspot_one_float, dxy_primaryvertex_one_float; std::vector<double> dxy_thebeamspot_second_float;
       std::vector<double> dxy_without_argument_float;

       // to check the BS error parameter
       std::vector<float> bs_x_error_base, bs_x_error_set , bs_y_error_base, bs_y_error_set;
       std::vector<float> bs_x_error_base_float, bs_x_error_set_float , bs_y_error_base_float, bs_y_error_set_float;


     // ----------member data ---------------------------
      edm::EDGetTokenT<reco::GenParticleCollection> prunedgenParticlesToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
      edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;
      edm::EDGetTokenT<edm::View<reco::Muon> > muonsToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<edm::View<reco::Track> > FLORIDAToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<edm::View<reco::Track> > globalToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<edm::View<reco::Track> > stdAloneToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<edm::View<reco::Track> > stdAloneVtxToken_;  
      //myfile1.close();//used to select what tracks to read from configuration file

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
VX_BS_Ana::VX_BS_Ana(const edm::ParameterSet& iConfig):
  prunedgenParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("prunedgenParticlesSrc"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"))),
  vertexToken_(consumes<edm::View<reco::Vertex> >(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc"))),
  muonsToken_(consumes<edm::View<reco::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc"))),
  FLORIDAToken_(consumes<edm::View<reco::Track> >(iConfig.getUntrackedParameter<edm::InputTag>("FLORIDASrc"))),
  globalToken_(consumes<edm::View<reco::Track> >(iConfig.getUntrackedParameter<edm::InputTag>("globalSrc"))),
  stdAloneToken_(consumes<edm::View<reco::Track> >(iConfig.getUntrackedParameter<edm::InputTag>("stdAloneSrc"))),
  stdAloneVtxToken_(consumes<edm::View<reco::Track> >(iConfig.getUntrackedParameter<edm::InputTag>("stdAloneVtxSrc")))
{
   
   //now do what ever initialization is needed
  passedEventsTree = new TTree("passedEvents","passedEvents");

}


VX_BS_Ana::~VX_BS_Ana()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


float KalmanEnergy(float px, float py, float pz, float mass){

    double E=TMath::Sqrt(px*px + py*py + pz*pz + mass*mass);
    
    return E;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
VX_BS_Ana::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   number_events = number_events+1;
   if(number_events == 126){
    std::cout<<"got till here"<<std::endl; 
   }
   std::cout<<"entering into anlayzer :  event number  :  "<<number_events<<std::endl;   
   edm::ESHandle<TransientTrackBuilder> ttkb; 
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);

   edm::Handle<reco::BeamSpot> beamSpot;
   iEvent.getByToken(beamSpotToken_, beamSpot);
   const reco::BeamSpot BS = *beamSpot;
   
   BS_x =  BS.position().X();
   BS_y =  BS.position().Y();
   BS_z =  BS.position().Z();
   BS_xErr =  BS.x0Error();
   BS_yErr =  BS.y0Error();
   BS_zErr =  BS.z0Error();
   

   BeamWidth_x = BS.BeamWidthX();
   BeamWidth_y = BS.BeamWidthY();
   BeamWidth_xErr = BS.BeamWidthXError();
   BeamWidth_yErr = BS.BeamWidthYError();
   
   xerr_bs = beamSpot->x0Error();
   yerr_bs = beamSpot->y0Error();


  // beamSpot->x0Error() = 0.001;
  // beamSpot->y0Error() = 0.001;
  // beamSpot->z0Error() = 1;

  xerr_bs_set = BS.x0Error();
  yerr_bs_set = BS.y0Error();

   // print out BS parameters for every event
 /* std::cout<<"//////////////////// BS and beamSpot ////////////////////"<<std::endl;
  std::cout<<beamSpot->position().Z()<<" :  "<<BS.position().Z()<<std::endl;
  std::cout<<" x error and y error on BS     "<<std::endl;
  std::cout<<beamSpot->x0Error()<<" :  "<<BS.x0Error()<<std::endl;
  std::cout<<beamSpot->y0Error()<<" :  "<<BS.y0Error()<<std::endl;

  std::cout<<"BS values for event : "<<std::endl;
  std::cout<<" parameters : x y z  : widthx widthy"<<std::endl;
  std::cout<<BS_x<<"  "<<BS_y<<"  "<<BS_z<<"  :  "<<BeamWidth_x<<"  "<<BeamWidth_y<<std::endl;
  std::cout<<"Beam spot x errors  :  "<<xerr_bs<<"  "<<xerr_bs_set<<"   "<<std::endl;
  std::cout<<"Beam spot y errors  :  "<<yerr_bs<<"  "<<yerr_bs_set<<"   "<<std::endl;*/

   Handle<edm::View<reco::Vertex> > vertex;
   iEvent.getByToken(vertexToken_,vertex);
   const reco::Vertex *PV = 0;
//   int theVertex = -1;
 //  std::cout<<"test vertex"<<vertex->size()<<std::endl;

    for (unsigned int i = 0; i < vertex->size(); i++) {
	   PV = &(vertex->at(i));        
     PV_x =  PV->position().X();
     PV_y =  PV->position().Y();
     PV_z =  PV->position().Z();
     PV_xErr =  PV->xError();
     PV_yErr =  PV->yError();
     PV_zErr =  PV->zError();
   
     //PV_Width_x = PV->BeamWidthX();
     //PV_Width_y = PV->BeamWidthY();
     //PV_Width_xErr = PV->BeamWidthXError();
     //PV_Width_yErr = PV->BeamWidthYError();
     //std::cout<<"test vertex position"<<PV->position()<<std::endl;
     //std::cout<<"test vertex position"<<vertex->position()<<std::endl;



	   
	   if (PV->isFake()) continue;
	   if (PV->ndof()<=4 || PV->position().Rho()>2.0 || fabs(PV->position().Z())>24.0) continue;
	   //theVertex=(int)i; 
 //          std::cout<<"test second time vertex position"<<std::endl; //PV->position()<<"   "<<vertex->at(i).Charge()<<std::endl;

//     break;
     } 
   
    //std::cout<<theVertex<<std::endl;
  // if(theVertex == -1) return;
   
  // std::cout<<"vertex thing"<<std::endl;
  // nVtx = vertex->size();
   
   edm::Handle<reco::GenParticleCollection> prunedgenParticles;
   iEvent.getByToken(prunedgenParticlesToken_,  prunedgenParticles);
   
   Handle<edm::View<reco::Muon> > muons;
   iEvent.getByToken(muonsToken_,  muons);
   
   Handle<edm::View<reco::Track> > glbMuon;
   iEvent.getByToken(globalToken_,  glbMuon);
   
   Handle<edm::View<reco::Track> > FLORIDA;
   iEvent.getByToken(FLORIDAToken_,  FLORIDA);

   Handle<edm::View<reco::Track> > stdAlone;
   iEvent.getByToken(stdAloneToken_,  stdAlone);

   Handle<edm::View<reco::Track> > stdAloneVtx;
   iEvent.getByToken(stdAloneVtxToken_,  stdAloneVtx);
  
   pt_GEN.clear(); eta_GEN.clear(); phi_GEN.clear(); mass_GEN.clear(); charge_GEN.clear();
   pt_trk.clear(); eta_trk.clear(); phi_trk.clear(); mass_trk.clear();
   dxy_BS_trk.clear(); charge_trk.clear(); d0_trk.clear();
   dxy_PV_trk.clear(); d0_BS_trk.clear(); d0_PV_trk.clear();

   pt_FLO.clear(); eta_FLO.clear(); phi_FLO.clear(); mass_FLO.clear(); dxy_BS_FLO.clear(); dxy_PV_FLO.clear();
   dxy_BS_FLO_vtx_second_test_float.clear(); dxy_BS_FLO_vtx_second_test.clear();
   charge_FLO.clear(); charge_FLO_float.clear();
   pt_FLO_vtx.clear(); eta_FLO_vtx.clear(); phi_FLO_vtx.clear(); mass_FLO_vtx.clear(); dxy_BS_FLO_vtx.clear();
   //pt_FLO_single.clear(); eta_FLO_single.clear(); phi_FLO_single.clear(); dxy_BS_FLO_single.clear();
   dxy_PV_FLO.clear(); dxy_PV_FLO_vtx.clear();  dxy_PV_FLO_float.clear(); dxy_PV_FLO_vtx_float.clear(); 


   pt_glb_trk.clear(); eta_glb_trk.clear(); phi_glb_trk.clear(); mass_glb_trk.clear(); dxy_BS_glb_trk.clear(); dxy_PV_glb_trk.clear(); charge_glb_trk.clear(); mass_glb_trk.clear();
  /* pt_glb_vtx.clear(); eta_glb_vtx.clear(); phi_glb_vtx.clear(); mass_glb_vtx.clear(); dxy_BS_glb_vtx.clear();
 /  pt_glb_single.clear(); eta_glb_single.clear(); phi_glb_single.clear(); dxy_BS_glb_single.clear();
  */
   pt_stdAlone.clear(); eta_stdAlone.clear(); phi_stdAlone.clear(); dxy_BS_stdAlone.clear(); dxy_PV_stdAlone.clear(); charge_stdAlone.clear(); mass_stdAlone.clear();
  /* pt_stdAlone_vtx.clear(); eta_stdAlone_vtx.clear(); phi_stdAlone_vtx.clear(); dxy_BS_stdAlone_vtx.clear(); dxy_PV_stdAlone_vtx.clear();
   */
   pt_GEN_float.clear(); eta_GEN_float.clear(); phi_GEN_float.clear(); mass_GEN_float.clear(); charge_GEN_float.clear(); 
   pt_trk_float.clear(); eta_trk_float.clear(); phi_trk_float.clear(); mass_trk_float.clear(); 
   dxy_BS_trk_float.clear();  charge_trk_float.clear(); d0_trk_float.clear();
   dxy_PV_trk_float.clear(); d0_BS_trk_float.clear(); d0_PV_trk_float.clear();
  
    pt_trk_BS.clear(); eta_trk_BS.clear(); phi_trk_BS.clear(); mass_trk_BS.clear(); 
    dxy_trk_BS_const.clear(); dxy_trk_PV_const.clear();
    pt_trk_BS_float.clear(); eta_trk_BS_float.clear(); phi_trk_BS_float.clear(); mass_trk_BS_float.clear(); 
    dxy_trk_BS_const_float.clear(); dxy_trk_PV_const_float.clear();
   
    pt_trk_BS_1.clear(); eta_trk_BS_1.clear(); phi_trk_BS_1.clear(); mass_trk_BS_1.clear(); 
    dxy_trk_BS_const_1.clear(); dxy_trk_PV_const_1.clear();
    pt_trk_BS_1_float.clear(); eta_trk_BS_1_float.clear(); phi_trk_BS_1_float.clear(); mass_trk_BS_1_float.clear(); 
    dxy_trk_BS_const_1_float.clear(); dxy_trk_PV_const_1_float.clear();
   

   
    pt_FLO_float.clear(); eta_FLO_float.clear(); phi_FLO_float.clear(); mass_FLO_float.clear(); dxy_BS_FLO_float.clear(); dxy_PV_FLO_float.clear();
   pt_FLO_vtx_float.clear(); eta_FLO_vtx_float.clear(); phi_FLO_vtx_float.clear(); mass_FLO_vtx_float.clear(); dxy_BS_FLO_vtx_float.clear();
   //pt_FLO_single_float.clear(); eta_FLO_single_float.clear(); phi_FLO_single_float.clear(); dxy_BS_FLO_single_float.clear();
  
   pt_glb_trk_float.clear(); eta_glb_trk_float.clear(); phi_glb_trk_float.clear(); mass_glb_trk_float.clear(); dxy_BS_glb_trk_float.clear(); dxy_PV_glb_trk_float.clear(); charge_glb_trk_float.clear();

/*   pt_glb_vtx_float.clear(); eta_glb_vtx_float.clear(); phi_glb_vtx_float.clear(); mass_glb_vtx_float.clear(); dxy_BS_glb_vtx_float.clear();
   pt_glb_single_float.clear(); eta_glb_single_float.clear(); phi_glb_single_float.clear(); dxy_BS_glb_single_float.clear(); 
*/
   pt_stdAlone_float.clear(); eta_stdAlone_float.clear(); phi_stdAlone_float.clear(); dxy_BS_stdAlone_float.clear(); dxy_PV_stdAlone_float.clear(); charge_stdAlone_float.clear(); mass_stdAlone_float.clear();
  /* pt_stdAlone_vtx_float.clear(); eta_stdAlone_vtx_float.clear(); phi_stdAlone_vtx_float.clear(); dxy_BS_stdAlone_vtx_float.clear(); dxy_PV_stdAlone_vtx_float.clear();
*/

// for standAlone tracks from muon container

   pt_stdAlone_trk_float.clear(); eta_stdAlone_trk_float.clear(); phi_stdAlone_trk_float.clear(); dxy_BS_stdAlone_trk_float.clear();  charge_stdAlone_trk_float.clear();
// dxy_PV_stdAlone_trk_float.clear();
   mass_stdAlone_trk_float.clear();
   pt_stdAlone_trk.clear(); eta_stdAlone_trk.clear(); phi_stdAlone_trk.clear(); dxy_BS_stdAlone_trk.clear(); charge_stdAlone_trk.clear(); mass_stdAlone_trk.clear();
// dxy_PV_stdAlone_trk.clear();

// for Beamspot plotting
   pt_glb_not_trk.clear(); pt_stdAlone_not_trk.clear();
   pt_glb_not_trk_float.clear(); pt_stdAlone_not_trk_float.clear();


   pt_tuned_float.clear(); eta_tuned_float.clear(); phi_tuned_float.clear(); mass_tuned_float.clear(); 
   dxy_BS_tuned_float.clear(); dxy_PV_tuned_float.clear();  charge_tuned_float.clear(); d0_tuned_float.clear();
   pt_tuned.clear(); eta_tuned.clear(); phi_tuned.clear(); mass_tuned.clear();
   dxy_BS_tuned.clear(); dxy_PV_tuned.clear();  charge_tuned.clear(); d0_tuned.clear();

   beam_x.clear(); beam_y.clear(); beam_z.clear(); beam_width_x.clear(); beam_width_y.clear(); 
   beam_x_err.clear(); beam_y_err.clear(); beam_z_err.clear(); beam_width_x_err.clear(); beam_width_y_err.clear();

   beam_x_float.clear(); beam_y_float.clear(); beam_z_float.clear(); beam_width_x_float.clear(); beam_width_y_float.clear();
   beam_x_err_float.clear(); beam_y_err_float.clear(); beam_z_err_float.clear(); beam_width_x_err_float.clear(); beam_width_y_err_float.clear();


   pv_x.clear(); pv_y.clear(); pv_z.clear(); pv_width_x.clear(); pv_width_y.clear(); 
   pv_x_err.clear(); pv_y_err.clear(); pv_z_err.clear(); pv_width_x_err.clear(); pv_width_y_err.clear();

   pv_x_float.clear(); pv_y_float.clear(); pv_z_float.clear(); pv_width_x_float.clear(); pv_width_y_float.clear();
   pv_x_err_float.clear(); pv_y_err_float.clear(); pv_z_err_float.clear(); pv_width_x_err_float.clear(); pv_width_y_err_float.clear();

   v_x.clear(); v_y.clear(); v_z.clear(); 
   dxy_thebeamspot_one.clear(); dxy_primaryvertex_one.clear(); dxy_thebeamspot_second.clear();
   v_x_float.clear(); v_y_float.clear(); v_z_float.clear(); 
   dxy_thebeamspot_one_float.clear(); dxy_primaryvertex_one_float.clear(); dxy_thebeamspot_second_float.clear();
 
   dxy_without_argument_float.clear(); dxy_without_argument.clear();
  

   bs_x_error_base.clear(); bs_y_error_base.clear(); bs_x_error_base_float.clear(); bs_y_error_base_float.clear();
   bs_x_error_set.clear(); bs_y_error_set.clear(); bs_x_error_set_float.clear(); bs_y_error_set_float.clear();


   // std::cout<<"before kv fitter"<<std::endl;
    SingleTrackVertexConstraint stvc;
    KalmanVertexFitter KVfitter(true);
    TransientVertex KVertex_BS;	
    //std::cout<<"after kv fitter"<<std::endl;

// beamspot plotting
    beam_x.push_back(BS_x);
    beam_y.push_back(BS_y);
    beam_z.push_back(BS_z);
    beam_width_x.push_back(BeamWidth_x);
    beam_width_y.push_back(BeamWidth_y);
    beam_x_err.push_back(BS_xErr);
    beam_y_err.push_back(BS_yErr);
    beam_z_err.push_back(BS_zErr);
    beam_width_x_err.push_back(BeamWidth_xErr);
    beam_width_y_err.push_back(BeamWidth_yErr);

    bs_x_error_base.push_back(xerr_bs);
    bs_x_error_set.push_back(xerr_bs_set);
    bs_y_error_base.push_back(yerr_bs);
    bs_y_error_set.push_back(yerr_bs_set);

    pv_x.push_back(PV_x);
    pv_y.push_back(PV_y);
    pv_z.push_back(PV_z);
  //  pv_width_x.push_back(PV_Width_x);
  //  pv_width_y.push_back(PV_Width_y);
    pv_x_err.push_back(PV_xErr);
    pv_y_err.push_back(PV_yErr);
    pv_z_err.push_back(PV_zErr);


     std::cout<<" entering muon container : muon size :  "<<muons->size()<<std::endl;
  // Fitting innerTrack   
     std::vector<reco::TransientTrack> ttv_trk;   
  // To count for number of inner Track
     
  //  int count = 0;
  //     double s; 
  
     int tracker_count =0;
     for(unsigned int i = 0; i < muons->size(); i++){
  //  test
  //       std::cout<<"here"<<std::endl;
         /*s = muons->refAt(i)->innerTrack()->pt();     
         if(s){
         count++;
         std::cout<<"inner Track  :  "<<count<<std::endl;
	      }*/
      
      std::cout<<"is muon"<<muons->size()<<std::endl;
 	    if(muons->refAt(i)->isTrackerMuon()){
        number_reco_muons = number_reco_muons +1;
   
        tracker_count++;

  // Fitting inner track to Beamspot      
      	ttv_trk.push_back(ttkb->build(*muons->refAt(i)->innerTrack()));
		    SingleTrackVertexConstraint::BTFtuple a = stvc.constrain(ttkb->build(*muons->refAt(i)->innerTrack()), *beamSpot);
       SingleTrackVertexConstraint::BTFtuple a1 = stvc.constrain(ttkb->build(*muons->refAt(i)->innerTrack()), BS)     ;
		    reco::Track innertrack_vtx_BS = a.get<1>().track();

		    reco::Track innertrack_vtx_BS_1 = a1.get<1>().track();
	  	  TLorentzVector tmp;
        TLorentzVector tmp1;

		    tmp.SetPxPyPzE(innertrack_vtx_BS.px(), innertrack_vtx_BS.py(), innertrack_vtx_BS.pz(), KalmanEnergy(innertrack_vtx_BS.px(), innertrack_vtx_BS.py(), innertrack_vtx_BS.pz(), 0.10565837));
		    pt_trk_BS.push_back(tmp.Pt());
		    eta_trk_BS.push_back(tmp.Eta());
		    phi_trk_BS.push_back(tmp.Phi());
		    mass_trk_BS.push_back(0.10565837);
        dxy_trk_BS_const.push_back(innertrack_vtx_BS.dxy(beamSpot->position()));
        dxy_trk_PV_const.push_back(innertrack_vtx_BS.dxy(PV->position()));
 
        std::cout<<" from 1st fitting beamspot with pointer"<<tmp.Pt()<<"  "<<innertrack_vtx_BS.dxy(beamSpot->position())<<std::endl;

        
    // To check difference btw beamSpot pointer and BS value

        

		    tmp1.SetPxPyPzE(innertrack_vtx_BS_1.px(), innertrack_vtx_BS_1.py(), innertrack_vtx_BS_1.pz(), KalmanEnergy(innertrack_vtx_BS_1.px(), innertrack_vtx_BS_1.py(), innertrack_vtx_BS_1.pz(), 0.10565837));


        std::cout<<" from 2nd fitting beamspot with pointer"<<tmp1.Pt()<<"  "<<innertrack_vtx_BS_1.dxy(BS.position())<<std::endl;
		    pt_trk_BS_1.push_back(tmp1.Pt());
		    eta_trk_BS_1.push_back(tmp1.Eta());
  	    phi_trk_BS_1.push_back(tmp1.Phi());
		    mass_trk_BS_1.push_back(0.10565837);
        dxy_trk_BS_const_1.push_back(innertrack_vtx_BS_1.dxy(BS.position()));
        dxy_trk_PV_const_1.push_back(innertrack_vtx_BS_1.dxy(PV->position()));

     	std::cout<<"a tracker muon  : muon track number  :  tracker count "<<tracker_count<<"pt value   :  "<<muons->refAt(i)->innerTrack()->pt()<<"  number "<<number_reco_muons<<std::endl;

		  	pt_trk.push_back(muons->refAt(i)->innerTrack()->pt());
 			  eta_trk.push_back(muons->refAt(i)->innerTrack()->eta());
 			  phi_trk.push_back(muons->refAt(i)->innerTrack()->phi());
        
        v_x.push_back(muons->refAt(i)->innerTrack()->vx());
        v_y.push_back(muons->refAt(i)->innerTrack()->vy());
        v_z.push_back(muons->refAt(i)->innerTrack()->vz());


 			  dxy_BS_trk.push_back(muons->refAt(i)->innerTrack()->dxy(beamSpot->position()));
        dxy_PV_trk.push_back(muons->refAt(i)->innerTrack()->dxy(PV->position()));
       //         dxy_vertex_test.push_back(muons->refAt(i)->innerTrack()->dxy(vertex->position()));
       //
        dxy_thebeamspot_one.push_back(muons->refAt(i)->innerTrack()->dxy(BS));
        dxy_without_argument.push_back(muons->refAt(i)->innerTrack()->dxy());
        dxy_primaryvertex_one.push_back(muons->refAt(i)->innerTrack()->dxy(PV->position()));
        dxy_thebeamspot_second.push_back(muons->refAt(i)->innerTrack()->dxy(beamSpot->position()));

        
 	 		  mass_trk.push_back(0.10565837);
        charge_trk.push_back(muons->refAt(i)->innerTrack()->charge());
			  d0_trk.push_back(muons->refAt(i)->innerTrack()->d0());

       ttv_trk.push_back(ttkb->build(*muons->refAt(i)->innerTrack()));
       // dxy_thebeamspot_second.push_back(muons->refAt(i)->innerTrack()->dxy(beamSpot->position(vz())));
        
 		
      }

      if(tracker_count ==0){
        std::cout<<"not a tracker muon here you see"<<std::endl;
      }



     	if(muons->refAt(i)->isGlobalMuon()){
        std::cout<<"global pt "<<muons->refAt(i)->globalTrack()->pt()<<std::endl;
			  glb_count++;

      //  std::cout<<"GLOBAL tracker muon  : count   :   "<<glb_count<<std::endl;
      	pt_glb_trk.push_back(muons->refAt(i)->globalTrack()->pt());
 			  eta_glb_trk.push_back(muons->refAt(i)->globalTrack()->eta());
 			  phi_glb_trk.push_back(muons->refAt(i)->globalTrack()->phi());
 			  mass_glb_trk.push_back(0.10565837);
 			  dxy_BS_glb_trk.push_back(muons->refAt(i)->globalTrack()->dxy(beamSpot->position()));
 	      charge_glb_trk.push_back(muons->refAt(i)->globalTrack()->charge());  
 			  dxy_PV_glb_trk.push_back(muons->refAt(i)->globalTrack()->dxy(PV->position()));
      }


      if(muons->refAt(i)->isStandAloneMuon()){
        std::cout<<"standalone pt  "<<muons->refAt(i)->standAloneMuon()->pt()<<std::endl;
      	standAlone_count++;
      //  std::cout<<"STANDALONE  muon  : count  :   "<<standAlone_count<<std::endl;

        pt_stdAlone_trk.push_back(muons->refAt(i)->standAloneMuon()->pt());
        eta_stdAlone_trk.push_back(muons->refAt(i)->standAloneMuon()->eta());
        phi_stdAlone_trk.push_back(muons->refAt(i)->standAloneMuon()->phi());
        mass_stdAlone_trk.push_back(0.10565837);
        dxy_BS_stdAlone_trk.push_back(muons->refAt(i)->standAloneMuon()->dxy(beamSpot->position()));
        charge_stdAlone_trk.push_back(muons->refAt(i)->standAloneMuon()->charge());   
      }

if(!muons->refAt(i)->isTrackerMuon() && muons->refAt(i)->isGlobalMuon())
{
std::cout<<"not a tracker muon but a global muon "<<std::endl;
pt_glb_not_trk.push_back(muons->refAt(i)->globalTrack()->pt());
pt_stdAlone_not_trk.push_back(muons->refAt(i)->standAloneMuon()->pt());

//std::cout<<"pT"<<muons->refAt(i)->globalTrack()->pt()<<"   eta"<<muons->refAt(i)->globalTrack()->eta()<<"  phi"<<muons->refAt(i)->globalTrack()->phi()<<std::endl;

}

if(!muons->refAt(i)->isTrackerMuon()){
std::cout<<"not a tracker muon  : count "<<std::endl; 
 }

//std::cout<<" moving to other muon in the event"<<std::endl;

//std::cout<<" number of tracks :  tracker  "<<tracker_count<<std::endl;
//std::cout<<" number of tracks :  global  "<<glb_count<<std::endl;
//std::cout<<" number of tracks :  standAlone  "<<standAlone_count<<std::endl;



//     }
// 	if(ttv_trk.size() == 2){
// 		KVertex_BS = KVfitter.vertex(ttv_trk, BS);
// 		if(KVertex_BS.hasRefittedTracks()){
// 			for(int i = 0; i < 2; i++){
// 				std::vector <reco::TransientTrack> ttrks_BS = KVertex_BS.refittedTracks();                  			
// 				TLorentzVector tmp;
// 				reco::Track track_vtx_BS = ttrks_BS.at(i).track();
// 				tmp.SetPxPyPzE(track_vtx_BS.px(), track_vtx_BS.py(), track_vtx_BS.pz(), KalmanEnergy(track_vtx_BS.px(), track_vtx_BS.py(), track_vtx_BS.pz(), 0.10565837));
// 				pt_trk_vtx.push_back(tmp.Pt());
// 				eta_trk_vtx.push_back(tmp.Eta());
// 				phi_trk_vtx.push_back(tmp.Phi());
// 				dxy_BS_trk_vtx.push_back(track_vtx_BS.dxy(BS.position()));
// 				dxy_PV_trk_vtx.push_back(track_vtx_PV.dxy(PV.position()));
// 			}
// 		}
// 	}
// 	else{
// 		pt_trk_vtx.push_back(-999);
// 		eta_trk_vtx.push_back(-999);
// 		phi_trk_vtx.push_back(-999);
// 		dxy_BS_trk_vtx.push_back(-999);
// 	}
     }    
/*
std::cout<<" number of tracks :  tracker  "<<tracker_count<<std::endl;
std::cout<<" number of tracks :  global  "<<glb_count<<std::endl;
std::cout<<" number of tracks :  standAlone  "<<standAlone_count<<std::endl;

 std::cout<<"ending muon track"<<std::endl;*/
    
     for(unsigned int i = 0; i < muons->size(); i++){
      std::cout<<"is muon   entering to tuneP  "<<muons->size()<<std::endl;
//        if(muon::isHighPtMuon(reco::Muon & muons, reco::Vertex & vertex)){
      	if(muons->refAt(i)->isGlobalMuon()){
     
        number_tune_muons = number_tune_muons +1;
   
       	std::cout<<"a global muon  : tunning track  "<<"pt value   :  "<<muons->refAt(i)->tunePMuonBestTrack()->pt()<<"  number "<<number_tune_muons<<std::endl;

		  	pt_tuned.push_back(muons->refAt(i)->tunePMuonBestTrack()->pt());
 			  eta_tuned.push_back(muons->refAt(i)->tunePMuonBestTrack()->eta());
 			  phi_tuned.push_back(muons->refAt(i)->tunePMuonBestTrack()->phi());
        
 			  dxy_BS_tuned.push_back(muons->refAt(i)->tunePMuonBestTrack()->dxy(beamSpot->position()));
        dxy_PV_tuned.push_back(muons->refAt(i)->tunePMuonBestTrack()->dxy(PV->position()));
 	 		  mass_tuned.push_back(0.10565837);
        charge_tuned.push_back(muons->refAt(i)->tunePMuonBestTrack()->charge());
			  d0_tuned.push_back(muons->refAt(i)->tunePMuonBestTrack()->d0());
      }
    }

  
     //std::cout<<"\tStart on GEN PARTICLE"<<std::endl;
    reco::GenParticleCollection::const_iterator genPart;
     float pt_value_gen = 5, pt_flo_value ;
    for(genPart = prunedgenParticles->begin(); genPart != prunedgenParticles->end(); ++genPart) {
      gen_count++;
      std::cout<<"gen part here  : pt  "<<genPart->pt()<<std::endl;
	       // if (fabs(genPart->pdgId())!=13){

    //         std::cout<<" this event is not muon  "<<std::endl; 
          //     continue;
          //    }

        //  std::cout<<" entering into genPart loop   "<<k<<"   "<<genPart->pdgId()<<std::endl;
        //  k=k+1;

//        if (fabs(genPart->pdgId())!=23) continue;
//        if (!(genPart->status()==1 || abs(genPart->pdgId())==15)) continue;
//         std::cout<<genPart->pdgId()<<"\t"<<genPart->pt()<<"\t"<<genPart->mass()<<std::endl;
    pt_value_gen = genPart->pt();
		pt_GEN.push_back(genPart->pt());
		eta_GEN.push_back(genPart->eta());
		phi_GEN.push_back(genPart->phi());
		mass_GEN.push_back(genPart->mass());
  	charge_GEN.push_back(genPart->charge());	
// 		std::cout<<genPart->dxyError()<<std::endl;

    }
  
 
    
 // this was commented before : it is for FLO muons   
		//std::cout<<"\tStart FLORIDA "<<std::endl;
 //   std::vector<reco::TransientTrack> ttv_FLO;
    number_flo_muons = number_flo_muons+ FLORIDA->size();
    std::cout<<" florida size  "<<FLORIDA->size()<<std::endl;
    if(FLORIDA->size() ==0){
      std::cout<<" not any event for flo"<<std::endl;
    }
    if(FLORIDA->size() ==1){
      std::cout<<"Flo single size "<<std::endl;
    for(unsigned int i = 0; i < FLORIDA->size(); i++){
      flo_count++;
	//  std::cout<<"entered florida loop "<<FLORIDA->size()<<std::endl;   
//     	if(!FLORIDA->refAt(i)->isGlobalMuon()) continue;
    	


	 /* if(FLORIDA->refAt(i)->pt() < 5) continue;
    if(fabs(FLORIDA->refAt(i)->eta())>2.4) continue;
    if(FLORIDA->refAt(i)->dxy(PV->position()) > 0.5) continue;
    if(FLORIDA->refAt(i)->dz(PV->position()) > 1) continue;    



}*/

//     	if(FLORIDA->refAt(i)->numberOfMatches() > 0) continue;    
    	
//     	double ip = fabs(FLORIDA->refAt(i)->dB(PV));
//     	double ipError = FLORIDA->refAt(i)->edB(PV);
//     	double sip3D = ip/ipError;
//     	if(sip3D > 4) continue;
 //   std::cout<<"   FLO HERE"<<". "<<std::endl;
    pt_flo_value = FLORIDA->refAt(i)->pt();

		pt_FLO.push_back(FLORIDA->refAt(i)->pt());
		eta_FLO.push_back(FLORIDA->refAt(i)->eta());
		phi_FLO.push_back(FLORIDA->refAt(i)->phi());
		mass_FLO.push_back(0.10565837);
		dxy_BS_FLO.push_back(FLORIDA->refAt(i)->dxy(beamSpot->position()));
		dxy_PV_FLO.push_back(FLORIDA->refAt(i)->dxy(PV->position()));
		charge_FLO.push_back(FLORIDA->refAt(i)->charge());
   }
  }
  int i;
   if(FLORIDA->size() == 2){
     if(abs(FLORIDA->refAt(0)->pt() - pt_value_gen) < abs(FLORIDA->refAt(1)->pt() - pt_value_gen)) 
      i = 0;
     else 
      i=1;

    std::cout<<"flo pt with 2 particles  here"<<FLORIDA->refAt(i)->pt()<<std::endl;
    std::cout<<"just printing"<<FLORIDA->refAt(0)->pt()<<"   "<<FLORIDA->refAt(1)->pt()<<std::endl;
    flo_count++;

//      if(FLORIDA->refAt(i)->nberOfMatches() > 0) continue;    
      
//      double ip = fabs(FLORIDA->refAt(i)->dB(PV));
//      double ipError = FLORIDA->refAt(i)->edB(PV);
//      double sip3D = ip/ipError;
//      if(sip3D > 4) continue;
 //   std::cout<<"   FLO HERE"<<". "<<std::endl;
    pt_flo_value = FLORIDA->refAt(i)->pt();

    pt_FLO.push_back(FLORIDA->refAt(i)->pt());
    eta_FLO.push_back(FLORIDA->refAt(i)->eta());
    phi_FLO.push_back(FLORIDA->refAt(i)->phi());
    mass_FLO.push_back(0.10565837);
    dxy_BS_FLO.push_back(FLORIDA->refAt(i)->dxy(beamSpot->position()));
    dxy_PV_FLO.push_back(FLORIDA->refAt(i)->dxy(PV->position()));
    charge_FLO.push_back(FLORIDA->refAt(i)->charge());
	  }

     if(FLORIDA->size() == 3){
       std::cout<<"error here"<<FLORIDA->refAt(0)->pt()<<"   "<<FLORIDA->refAt(1)->pt()<<"  "<<FLORIDA->refAt(2)->pt()<<pt_value_gen<<"  "<<std::endl;
       if(abs(FLORIDA->refAt(0)->pt() - pt_value_gen) < 0.2) 
      i = 0;
     else if(abs(FLORIDA->refAt(1)->pt() - pt_value_gen)<0.2)
      i=1;
    else if(abs(FLORIDA->refAt(2)->pt() - pt_value_gen)<0.2)
      i=2;

  //  std::cout<<"flo pt with 3 particles  here"<<FLORIDA->refAt(i)->pt()<<std::endl;
    std::cout<<"just printing"<<FLORIDA->refAt(0)->pt()<<"   "<<FLORIDA->refAt(1)->pt()<<std::endl;
    flo_count++;

//      if(FLORIDA->refAt(i)->nberOfMatches() > 0) continue;    
      
//      double ip = fabs(FLORIDA->refAt(i)->dB(PV));
//      double ipError = FLORIDA->refAt(i)->edB(PV);
//      double sip3D = ip/ipError;
//      if(sip3D > 4) continue;
 //   std::cout<<"   FLO HERE"<<". "<<std::endl;
    pt_flo_value = FLORIDA->refAt(i)->pt();

    pt_FLO.push_back(FLORIDA->refAt(i)->pt());
    eta_FLO.push_back(FLORIDA->refAt(i)->eta());
    phi_FLO.push_back(FLORIDA->refAt(i)->phi());
    mass_FLO.push_back(0.10565837);
    dxy_BS_FLO.push_back(FLORIDA->refAt(i)->dxy(beamSpot->position()));
    dxy_PV_FLO.push_back(FLORIDA->refAt(i)->dxy(PV->position()));
    charge_FLO.push_back(FLORIDA->refAt(i)->charge());
    }
    if(gen_count != flo_count){
      std::cout<<" gen but not flo  : event number "<<number_events<<"  pt gen : "<<pt_value_gen<<" pt flo : "<<pt_flo_value<<std::endl;
    }
  //	ttv_FLO.push_back(ttkb->build(*FLORIDA->refAt(i)));
			
		//SingleTrackVertexConstraint::BTFtuple a = stvc.constrain(ttkb->build(*FLORIDA->refAt(i)), *beamSpot);
		//reco::Track track_vtx_BS = a.get<1>().track();

	/*	TLorentzVector tmp;
		tmp.SetPxPyPzE(track_vtx_BS.px(), track_vtx_BS.py(), track_vtx_BS.pz(), KalmanEnergy(track_vtx_BS.px(), track_vtx_BS.py(), track_vtx_BS.pz(), 0.10565837));
		pt_FLO_single.push_back(tmp.Pt());
		eta_FLO_single.push_back(tmp.Eta());
		phi_FLO_single.push_back(tmp.Phi());
		dxy_BS_FLO_single.push_back(track_vtx_BS.dxy(BS.position()));*/
	
	

	//std::cout<<"FLO ends"<<std::endl;
// 	KalmanVertexFitter KVfitter(true);
// 	TransientVertex KVertex_BS = KVfitter.vertex(ttv_FLO, BS);
	/*if(ttv_FLO.size() == 2){
		KVertex_BS = KVfitter.vertex(ttv_FLO, BS);
		if(KVertex_BS.hasRefittedTracks()){
			for(int i = 0; i < 2; i++){
				std::vector <reco::TransientTrack> ttrks_BS = KVertex_BS.refittedTracks();                  			
				TLorentzVector tmp;
				reco::Track track_vtx_BS = ttrks_BS.at(i).track();
				tmp.SetPxPyPzE(track_vtx_BS.px(), track_vtx_BS.py(), track_vtx_BS.pz(), KalmanEnergy(track_vtx_BS.px(), track_vtx_BS.py(), track_vtx_BS.pz(), 0.10565837));
				pt_FLO_vtx.push_back(tmp.Pt());
				eta_FLO_vtx.push_back(tmp.Eta());
				phi_FLO_vtx.push_back(tmp.Phi());
				dxy_BS_FLO_vtx.push_back(track_vtx_BS.dxy(BS.position()));
        dxy_BS_FLO_vtx_second_test.push_back(track_vtx_BS.dxy(beamSpot->position()));
        dxy_PV_FLO_vtx.push_back(track_vtx_BS.dxy(PV->position()));

			}
		}
	}
	else{
		pt_FLO_vtx.push_back(-999);
		eta_FLO_vtx.push_back(-999);
		phi_FLO_vtx.push_back(-999);
		dxy_BS_FLO_vtx.push_back(-999);
    dxy_PV_FLO_vtx.push_back(-999);
    dxy_BS_FLO_vtx_second_test.push_back(-999);


	}

 // std::cout<<" FLO completely ends Global track starts"<<std::endl;

    std::vector<reco::TransientTrack> ttv_glb;
	std::cout<<"\tStart Global track"<<std::endl;
    for(unsigned int i = 0; i < glbMuon->size(); i++){
 

    	if(glbMuon->refAt(i)->pt() < 5) continue;
    	if(fabs(glbMuon->refAt(i)->eta())>2.4) continue;
    	if(glbMuon->refAt(i)->dxy(PV->position()) > 0.5) continue;
    	if(glbMuon->refAt(i)->dz(PV->position()) > 1) continue;    

   		pt_glb.push_back(glbMuon->refAt(i)->pt());
		eta_glb.push_back(glbMuon->refAt(i)->eta());
		phi_glb.push_back(glbMuon->refAt(i)->phi());
		mass_glb.push_back(0.10565837);
		dxy_BS_glb.push_back(glbMuon->refAt(i)->dxy(beamSpot->position()));
		dxy_PV_glb.push_back(glbMuon->refAt(i)->dxy(PV->position()));

		ttv_glb.push_back(ttkb->build(*glbMuon->refAt(i)));

		SingleTrackVertexConstraint::BTFtuple a = stvc.constrain(ttkb->build(*glbMuon->refAt(i)), *beamSpot);
		reco::Track track_vtx_BS = a.get<1>().track();

		TLorentzVector tmp;
		tmp.SetPxPyPzE(track_vtx_BS.px(), track_vtx_BS.py(), track_vtx_BS.pz(), KalmanEnergy(track_vtx_BS.px(), track_vtx_BS.py(), track_vtx_BS.pz(), 0.10565837));
		pt_glb_single.push_back(tmp.Pt());
		eta_glb_single.push_back(tmp.Eta());
		phi_glb_single.push_back(tmp.Phi());
		dxy_BS_glb_single.push_back(track_vtx_BS.dxy(BS.position()));

    }
    if(ttv_glb.size() == 2){
		KVertex_BS = KVfitter.vertex(ttv_glb, BS);
		if(KVertex_BS.hasRefittedTracks()){
			for(int i = 0; i < 2; i++){
				std::vector <reco::TransientTrack> ttrks_BS = KVertex_BS.refittedTracks();                  			
				TLorentzVector tmp;
				reco::Track track_vtx_BS = ttrks_BS.at(i).track();
				tmp.SetPxPyPzE(track_vtx_BS.px(), track_vtx_BS.py(), track_vtx_BS.pz(), KalmanEnergy(track_vtx_BS.px(), track_vtx_BS.py(), track_vtx_BS.pz(), 0.10565837));
				pt_glb_vtx.push_back(tmp.Pt());
				eta_glb_vtx.push_back(tmp.Eta());
				phi_glb_vtx.push_back(tmp.Phi());
				dxy_BS_glb_vtx.push_back(track_vtx_BS.dxy(BS.position()));
			}
		}
	}
	else{
		pt_glb_vtx.push_back(-999);
		eta_glb_vtx.push_back(-999);
		phi_glb_vtx.push_back(-999);
		dxy_BS_glb_vtx.push_back(-999);
	}
	
*/
  //  std::cout<<"entering into stdAlone track"<<std::endl;	
	
    for(unsigned int i = 0; i < stdAlone->size(); i++){

		pt_stdAlone.push_back(stdAlone->refAt(i)->pt());
		eta_stdAlone.push_back(stdAlone->refAt(i)->eta());
		phi_stdAlone.push_back(stdAlone->refAt(i)->phi());
		dxy_BS_stdAlone.push_back(stdAlone->refAt(i)->dxy(beamSpot->position()));
    charge_stdAlone.push_back(stdAlone->refAt(i)->charge());
//		dxy_PV_stdAlone.push_back(stdAlone->refAt(i)->dxy(PV->position()));
	
	}
/*
    for(unsigned int i = 0; i < stdAloneVtx->size(); i++){

		pt_stdAlone_vtx.push_back(stdAloneVtx->refAt(i)->pt());
		eta_stdAlone_vtx.push_back(stdAloneVtx->refAt(i)->eta());
		phi_stdAlone_vtx.push_back(stdAloneVtx->refAt(i)->phi());
		dxy_BS_stdAlone_vtx.push_back(stdAloneVtx->refAt(i)->dxy(beamSpot->position()));
		dxy_PV_stdAlone_vtx.push_back(stdAloneVtx->refAt(i)->dxy(PV->position()));
	
	}

*/
	pt_GEN_float.assign(pt_GEN.begin(),pt_GEN.end());
	eta_GEN_float.assign(eta_GEN.begin(),eta_GEN.end());
	phi_GEN_float.assign(phi_GEN.begin(),phi_GEN.end());
	mass_GEN_float.assign(mass_GEN.begin(),mass_GEN.end());
        charge_GEN_float.assign(charge_GEN.begin(),charge_GEN.end());

        pt_trk_float.assign(pt_trk.begin(),pt_trk.end());
        eta_trk_float.assign(eta_trk.begin(),eta_trk.end());
        phi_trk_float.assign(phi_trk.begin(),phi_trk.end());
        mass_trk_float.assign(mass_trk.begin(),mass_trk.end());
        dxy_BS_trk_float.assign(dxy_BS_trk.begin(),dxy_BS_trk.end());
        charge_trk_float.assign(charge_trk.begin(),charge_trk.end());
        d0_trk_float.assign(d0_trk.begin(),d0_trk.end());

        pt_trk_BS_float.assign(pt_trk_BS.begin(),pt_trk_BS.end());
        eta_trk_BS_float.assign(eta_trk_BS.begin(),eta_trk_BS.end());
        phi_trk_BS_float.assign(phi_trk_BS.begin(),phi_trk_BS.end());
        mass_trk_BS_float.assign(mass_trk_BS.begin(),mass_trk_BS.end());
        dxy_trk_BS_const_float.assign(dxy_trk_BS_const.begin(),dxy_trk_BS_const.end());
        dxy_trk_PV_const_float.assign(dxy_trk_PV_const.begin(),dxy_trk_PV_const.end());

        pt_trk_BS_1_float.assign(pt_trk_BS_1.begin(),pt_trk_BS_1.end());
        eta_trk_BS_1_float.assign(eta_trk_BS_1.begin(),eta_trk_BS_1.end());
        phi_trk_BS_1_float.assign(phi_trk_BS_1.begin(),phi_trk_BS_1.end());
        mass_trk_BS_1_float.assign(mass_trk_BS_1.begin(),mass_trk_BS_1.end());
        dxy_trk_BS_const_1_float.assign(dxy_trk_BS_const_1.begin(),dxy_trk_BS_const_1.end());
        dxy_trk_PV_const_1_float.assign(dxy_trk_PV_const_1.begin(),dxy_trk_PV_const_1.end());



         pt_tuned_float.assign(pt_tuned.begin(),pt_tuned.end());
         eta_tuned_float.assign(eta_tuned.begin(),eta_tuned.end());
         phi_tuned_float.assign(phi_tuned.begin(),phi_tuned.end());
         mass_tuned_float.assign(mass_tuned.begin(),mass_tuned.end());
         dxy_BS_tuned_float.assign(dxy_BS_tuned.begin(),dxy_BS_tuned.end());
          dxy_PV_tuned_float.assign(dxy_PV_tuned.begin(),dxy_PV_tuned.end());
         charge_tuned_float.assign(charge_tuned.begin(),charge_tuned.end());
         d0_tuned_float.assign(d0_tuned.begin(),d0_tuned.end());
       
         pt_glb_trk_float.assign(pt_glb_trk.begin(),pt_glb_trk.end());
        eta_glb_trk_float.assign(eta_glb_trk.begin(),eta_glb_trk.end());
        phi_glb_trk_float.assign(phi_glb_trk.begin(),phi_glb_trk.end());
        mass_glb_trk_float.assign(mass_glb_trk.begin(),mass_glb_trk.end());
        dxy_BS_glb_trk_float.assign(dxy_BS_glb_trk.begin(),dxy_BS_glb_trk.end());
        dxy_PV_glb_trk_float.assign(dxy_PV_glb_trk.begin(),dxy_PV_glb_trk.end());
        charge_glb_trk_float.assign(charge_glb_trk.begin(),charge_glb_trk.end());

       pt_glb_not_trk_float.assign(pt_glb_not_trk.begin(), pt_glb_not_trk.end());
       pt_stdAlone_not_trk_float.assign(pt_stdAlone_not_trk.begin(), pt_stdAlone_not_trk.end());

        dxy_PV_trk_float.assign(dxy_PV_trk.begin(),dxy_PV_trk.end());
        d0_BS_trk_float.assign(d0_BS_trk.begin(),d0_BS_trk.end());
   //       d0_PV_trk_float.assign(d0_BS_trk.begin(),d0_PV_trk.end());



	pt_FLO_float.assign(pt_FLO.begin(),pt_FLO.end());
	eta_FLO_float.assign(eta_FLO.begin(),eta_FLO.end());
	phi_FLO_float.assign(phi_FLO.begin(),phi_FLO.end());
	mass_FLO_float.assign(mass_FLO.begin(),mass_FLO.end());
	dxy_BS_FLO_float.assign(dxy_BS_FLO.begin(),dxy_BS_FLO.end());
	dxy_PV_FLO_float.assign(dxy_PV_FLO.begin(),dxy_PV_FLO.end());
  charge_FLO_float.assign(charge_FLO.begin(),charge_FLO.end());



	pt_FLO_vtx_float.assign(pt_FLO_vtx.begin(),pt_FLO_vtx.end());
	eta_FLO_vtx_float.assign(eta_FLO_vtx.begin(),eta_FLO_vtx.end());
	phi_FLO_vtx_float.assign(phi_FLO_vtx.begin(),phi_FLO_vtx.end());
	mass_FLO_vtx_float.assign(mass_FLO_vtx.begin(),mass_FLO_vtx.end());
	dxy_BS_FLO_vtx_float.assign(dxy_BS_FLO_vtx.begin(),dxy_BS_FLO_vtx.end());
  dxy_PV_FLO_vtx_float.assign(dxy_PV_FLO_vtx.begin(),dxy_PV_FLO_vtx.end());
  dxy_BS_FLO_vtx_second_test_float.assign(dxy_BS_FLO_vtx_second_test.begin(),dxy_BS_FLO_vtx_second_test.end());

	/*
  pt_FLO_single_float.assign(pt_FLO_single.begin(),pt_FLO_single.end());
	eta_FLO_single_float.assign(eta_FLO_single.begin(),eta_FLO_single.end());
	phi_FLO_single_float.assign(phi_FLO_single.begin(),phi_FLO_single.end());
	dxy_BS_FLO_single_float.assign(dxy_BS_FLO_single.begin(),dxy_BS_FLO_single.end());

	pt_glb_float.assign(pt_glb.begin(),pt_glb.end());
	eta_glb_float.assign(eta_glb.begin(),eta_glb.end());
	phi_glb_float.assign(phi_glb.begin(),phi_glb.end());
	mass_glb_float.assign(mass_glb.begin(),mass_glb.end());
	dxy_BS_glb_float.assign(dxy_BS_glb.begin(),dxy_BS_glb.end());
	dxy_PV_glb_float.assign(dxy_PV_glb.begin(),dxy_PV_glb.end());

	pt_glb_vtx_float.assign(pt_glb_vtx.begin(),pt_glb_vtx.end());
	eta_glb_vtx_float.assign(eta_glb_vtx.begin(),eta_glb_vtx.end());
	phi_glb_vtx_float.assign(phi_glb_vtx.begin(),phi_glb_vtx.end());
	mass_glb_vtx_float.assign(mass_glb_vtx.begin(),mass_glb_vtx.end());
	dxy_BS_glb_vtx_float.assign(dxy_BS_glb_vtx.begin(),dxy_BS_glb_vtx.end());

	pt_glb_single_float.assign(pt_glb_single.begin(),pt_glb_single.end());
	eta_glb_single_float.assign(eta_glb_single.begin(),eta_glb_single.end());
	phi_glb_single_float.assign(phi_glb_single.begin(),phi_glb_single.end());
	dxy_BS_glb_single_float.assign(dxy_BS_glb_single.begin(),dxy_BS_glb_single.end());
*/
	pt_stdAlone_float.assign(pt_stdAlone.begin(),pt_stdAlone.end());
	eta_stdAlone_float.assign(eta_stdAlone.begin(),eta_stdAlone.end());
	phi_stdAlone_float.assign(phi_stdAlone.begin(),phi_stdAlone.end());
	dxy_BS_stdAlone_float.assign(dxy_BS_stdAlone.begin(),dxy_BS_stdAlone.end());
  charge_stdAlone_float.assign(charge_stdAlone.begin(),charge_stdAlone.end());
//	dxy_PV_stdAlone_float.assign(dxy_PV_stdAlone.begin(),dxy_PV_stdAlone.end());

        pt_stdAlone_trk_float.assign(pt_stdAlone_trk.begin(),pt_stdAlone_trk.end());
        eta_stdAlone_trk_float.assign(eta_stdAlone_trk.begin(),eta_stdAlone_trk.end());
        phi_stdAlone_trk_float.assign(phi_stdAlone_trk.begin(),phi_stdAlone_trk.end());
        dxy_BS_stdAlone_trk_float.assign(dxy_BS_stdAlone_trk.begin(),dxy_BS_stdAlone_trk.end());
        charge_stdAlone_trk_float.assign(charge_stdAlone_trk.begin(),charge_stdAlone_trk.end());

/*	pt_stdAlone_vtx_float.assign(pt_stdAlone_vtx.begin(),pt_stdAlone_vtx.end());
	eta_stdAlone_vtx_float.assign(eta_stdAlone_vtx.begin(),eta_stdAlone_vtx.end());
	phi_stdAlone_vtx_float.assign(phi_stdAlone_vtx.begin(),phi_stdAlone_vtx.end());
	dxy_BS_stdAlone_vtx_float.assign(dxy_BS_stdAlone_vtx.begin(),dxy_BS_stdAlone_vtx.end());
	dxy_PV_stdAlone_vtx_float.assign(dxy_PV_stdAlone_vtx.begin(),dxy_PV_stdAlone_vtx.end());
*/

// beamspot plotting
  
        beam_x_float.assign(beam_x.begin(),beam_x.end());
        beam_y_float.assign(beam_y.begin(),beam_y.end());
        beam_z_float.assign(beam_z.begin(),beam_z.end());
        beam_width_x_float.assign(beam_width_x.begin(),beam_width_x.end());
        beam_width_y_float.assign(beam_width_y.begin(),beam_width_y.end());
        beam_x_err_float.assign(beam_x_err.begin(),beam_x_err.end());
        beam_y_err_float.assign(beam_y_err.begin(),beam_y_err.end());
        beam_z_err_float.assign(beam_z_err.begin(),beam_z_err.end());
        beam_width_x_err_float.assign(beam_width_x_err.begin(),beam_width_x_err.end());
        beam_width_y_err_float.assign(beam_width_y_err.begin(),beam_width_y_err.end());	


        pv_x_float.assign(pv_x.begin(),pv_x.end());
        pv_y_float.assign(pv_y.begin(),pv_y.end());
        pv_z_float.assign(pv_z.begin(),pv_z.end());
     //   pv_width_x_float.assign(pv_width_x.begin(),pv_width_x.end());
     //   pv_width_y_float.assign(pv_width_y.begin(),pv_width_y.end());
        pv_x_err_float.assign(pv_x_err.begin(),pv_x_err.end());
        pv_y_err_float.assign(pv_y_err.begin(),pv_y_err.end());
        pv_z_err_float.assign(pv_z_err.begin(),pv_z_err.end());
     //   pv_width_x_err_float.assign(pv_width_x_err.begin(),pv_width_x_err.end());
     //   pv_width_y_err_float.assign(pv_width_y_err.begin(),pv_width_y_err.end()); 
        v_x_float.assign(v_x.begin(),v_x.end());
        v_y_float.assign(v_y.begin(),v_y.end());
        v_z_float.assign(v_z.begin(),v_z.end());
        dxy_thebeamspot_one_float.assign(dxy_thebeamspot_one.begin(),dxy_thebeamspot_one.end());
        dxy_primaryvertex_one_float.assign(dxy_primaryvertex_one.begin(),dxy_primaryvertex_one.end());
        dxy_thebeamspot_second_float.assign(dxy_thebeamspot_second.begin(),dxy_thebeamspot_second.end());

        dxy_without_argument_float.assign(dxy_without_argument.begin(), dxy_without_argument.end());



        bs_x_error_base_float.assign(bs_x_error_base.begin(), bs_x_error_base.end());
        bs_y_error_base_float.assign(bs_y_error_base.begin(), bs_y_error_base.end());
        bs_x_error_set_float.assign(bs_x_error_set.begin(), bs_x_error_set.end());  
        bs_y_error_set_float.assign(bs_y_error_set.begin(), bs_y_error_set.end());
	

  std::cout<<std::endl;
  std::cout<<std::endl;

  std::cout<<number_flo_muons<<std::endl;
 

  std::cout<<std::endl;
  std::cout<<std::endl;

  passedEventsTree->Fill();     

  if(number_events == 126){
   std::cout<<"here"<<std::endl;
   }


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example", pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
VX_BS_Ana::beginJob()
{
  number_events =0;
	bookPassedEventTree("passedEvents", passedEventsTree);
 // int glb_count =0, gen_count =0, flo_count =0;
}

// ------------ method called once each job just after ending the event loop  ------------
void
VX_BS_Ana::endJob()
{
  std::cout<<"total number of events"<<number_events<<std::endl; 
  std::cout<<" number reco muons "<<number_reco_muons<<std::endl;
  std::cout<<" number reco muons "<<number_tune_muons<<std::endl;

  std::cout<<" global count"<<glb_count<<std::endl;
  std::cout<<" flo count"<<flo_count<<std::endl;
  std::cout<<" gen count"<<gen_count<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VX_BS_Ana::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use,  even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use,  remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", "ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}





void VX_BS_Ana::bookPassedEventTree(TString treeName, TTree *tree)
{     
    
    
    std::cout<<"Iiiiiiiiiii "<<number_flo_muons<<std::endl;
    using namespace edm;
    using namespace std;

//     tree->Branch("Run",&Run,"Run/l");
//     tree->Branch("Event",&Event,"Event/l");
//     tree->Branch("LumiSect",&LumiSect,"LumiSect/l");
   /* tree->Branch("BS_x", &BS_x, "BS_x/F");
    tree->Branch("BS_y", &BS_y, "BS_y/F");
    tree->Branch("BS_z", &BS_z, "BS_z/F");
    tree->Branch("BS_xErr", &BS_xErr, "BS_xErr/F");
    tree->Branch("BS_yErr", &BS_yErr, "BS_yErr/F");
    tree->Branch("BS_zErr", &BS_zErr, "BS_zErr/F");
    tree->Branch("BeamWidth_x", &BeamWidth_x, "BeamWidth_x/F");
    tree->Branch("BeamWidth_y", &BeamWidth_y, "BeamWidth_y/F")f;
    tree->Branch("BeamWidth_xErr", &BeamWidth_xErr, "BeamWidth_xErr/F");
    tree->Branch("BeamWidth_yErr", &BeamWidth_yErr, "BeamWidth_yErr/F");*/
//     tree->Branch("nVtx", &nVtx, "nVtx/I");
    tree->Branch("pt_GEN", &pt_GEN_float);
    tree->Branch("eta_GEN", &eta_GEN_float);
    tree->Branch("phi_GEN", &phi_GEN_float);
    tree->Branch("mass_GEN", &mass_GEN_float);
    tree->Branch("charge_GEN",&charge_GEN_float);
    /*tree->Branch("pt_FLO", &pt_FLO_float);
    tree->Branch("eta_FLO", &eta_FLO_float);
    tree->Branch("phi_FLO", &phi_FLO_float);
    tree->Branch("mass_FLO", &mass_FLO_float);
    tree->Branch("dxy_BS_FLO", &dxy_BS_FLO_float);
    tree->Branch("dxy_PV_FLO", &dxy_PV_FLO_float);
    tree->Branch("pt_FLO_vtx", &pt_FLO_vtx_float);
    tree->Branch("eta_FLO_vtx", &eta_FLO_vtx_float);
    tree->Branch("phi_FLO_vtx", &phi_FLO_vtx_float);
    tree->Branch("mass_FLO_vtx", &mass_FLO_vtx_float);
    tree->Branch("dxy_BS_FLO_vtx", &dxy_BS_FLO_vtx_float);*/
   // tree->Branch("pt_FLO_single", &pt_FLO_single_float);
   //  tree->Branch("eta_FLO_single", &eta_FLO_single_float);
   //tree->Branch("phi_FLO_single", &phi_FLO_single_float);
   // tree->Branch("dxy_BS_FLO_single", &dxy_BS_FLO_single_float);
    tree->Branch("pt_glb", &pt_glb_float);
    tree->Branch("eta_glb", &eta_glb_float);
    tree->Branch("phi_glb", &phi_glb_float);
    tree->Branch("mass_glb", &mass_glb_float);
    tree->Branch("dxy_BS_glb", &dxy_BS_glb_float);
    tree->Branch("dxy_PV_glb", &dxy_PV_glb_float);
    tree->Branch("pt_glb_vtx", &pt_glb_vtx_float);
    tree->Branch("eta_glb_vtx", &eta_glb_vtx_float);
    tree->Branch("phi_glb_vtx", &phi_glb_vtx_float);
    tree->Branch("mass_glb_vtx", &mass_glb_vtx_float);
    tree->Branch("dxy_BS_glb_vtx", &dxy_BS_glb_vtx_float);
    tree->Branch("pt_glb_single", &pt_glb_single_float);
    tree->Branch("eta_glb_single", &eta_glb_single_float);
    tree->Branch("phi_glb_single", &phi_glb_single_float);
    tree->Branch("dxy_BS_glb_single", &dxy_BS_glb_single_float);
    tree->Branch("pt_stdAlone", &pt_stdAlone_float);
    tree->Branch("eta_stdAlone", &eta_stdAlone_float);
    tree->Branch("phi_stdAlone", &phi_stdAlone_float);
    tree->Branch("charge_stdAlone",&charge_stdAlone_float);


    tree->Branch("dxy_BS_stdAlone", &dxy_BS_stdAlone_float);
    tree->Branch("pt_stdAlone_vtx", &pt_stdAlone_vtx_float);
    tree->Branch("eta_stdAlone_vtx", &eta_stdAlone_vtx_float);
    tree->Branch("phi_stdAlone_vtx", &phi_stdAlone_vtx_float);
    tree->Branch("dxy_BS_stdAlone_vtx", &dxy_BS_stdAlone_vtx_float);


    tree->Branch("pt_trk", &pt_trk_float);
    tree->Branch("eta_trk", &eta_trk_float);
    tree->Branch("phi_trk", &phi_trk_float);
    tree->Branch("mass_trk", &mass_trk_float);
    tree->Branch("dxy_BS_trk", &dxy_BS_trk_float);
    tree->Branch("charge_trk",&charge_trk_float);
    tree->Branch("d0_trk", &d0_trk_float);

// Inner track with BS constraint
    tree->Branch("pt_trk_BS", &pt_trk_BS_float);
    tree->Branch("eta_trk_BS", &eta_trk_BS_float);
    tree->Branch("phi_trk_BS", &phi_trk_BS_float);
    tree->Branch("mass_trk_BS", &mass_trk_BS_float);
    tree->Branch("dxy_trk_BS_const", &dxy_trk_BS_const_float);
    tree->Branch("dxy_trk_PV_const", &dxy_trk_PV_const_float);

    tree->Branch("pt_trk_BS_1", &pt_trk_BS_1_float);
    tree->Branch("eta_trk_BS_1", &eta_trk_BS_1_float);
    tree->Branch("phi_trk_BS_1", &phi_trk_BS_1_float);
    tree->Branch("mass_trk_BS_1", &mass_trk_BS_1_float);
    tree->Branch("dxy_trk_BS_const_1", &dxy_trk_BS_const_1_float);
    tree->Branch("dxy_trk_PV_const_1", &dxy_trk_PV_const_1_float);



    tree->Branch("pt_tuned", &pt_tuned_float);
    tree->Branch("eta_tuned", &eta_tuned_float);
    tree->Branch("phi_tuned", &phi_tuned_float);
    tree->Branch("mass_tuned", &mass_tuned_float);
    tree->Branch("dxy_BS_tuned", &dxy_BS_tuned_float);
    tree->Branch("dxy_PV_tuned", &dxy_PV_tuned_float);
    tree->Branch("charge_tuned",&charge_tuned_float);
    tree->Branch("d0_tuned", &d0_tuned_float);


    tree->Branch("pt_glb_trk", &pt_glb_trk_float);
    tree->Branch("eta_glb_trk", &eta_glb_trk_float);
    tree->Branch("phi_glb_trk", &phi_glb_trk_float);
    tree->Branch("mass_glb_trk", &mass_glb_trk_float);
    tree->Branch("dxy_BS_glb_trk", &dxy_BS_glb_trk_float);
    tree->Branch("dxy_PV_glb_trk", &dxy_PV_glb_trk_float);
    tree->Branch("charge_glb_trk",&charge_glb_trk_float);

    tree->Branch("pt_stdAlone_trk", &pt_stdAlone_trk_float);
    tree->Branch("eta_stdAlone_trk", &eta_stdAlone_trk_float);
    tree->Branch("phi_stdAlone_trk", &phi_stdAlone_trk_float);
    tree->Branch("mass_stdAlone_trk", &mass_stdAlone_trk_float);
    tree->Branch("dxy_BS_stdAlone_trk", &dxy_BS_stdAlone_trk_float);
    tree->Branch("charge_stdAlone_trk",&charge_stdAlone_trk_float);

    tree->Branch("pt_global_not_tracker",&pt_glb_not_trk_float);
    tree->Branch("pt_standAlone_not_tracker",&pt_stdAlone_not_trk_float);

    tree->Branch("dxy_PV_trk", &dxy_PV_trk_float);
    tree->Branch("d0_BS_trk", &d0_BS_trk_float);
 //   tree->Branch("d0_PV_trk", &d0_PV_trk_float);

// plotting BS


    tree->Branch("beam_x", &beam_x_float);
    tree->Branch("beam_y", &beam_y_float);
    tree->Branch("beam_z", &beam_z_float);
    tree->Branch("beam_x_err", &beam_x_err_float);
    tree->Branch("beam_y_err", &beam_y_err_float);
    tree->Branch("beam_z_err", &beam_z_err_float);

    tree->Branch("beam_width_x", &beam_width_x_float);
    tree->Branch("beam_width_y", &beam_width_y_float);
    tree->Branch("beam_width_x_err", &beam_width_x_err_float);
    tree->Branch("beam_width_y_err", &beam_width_y_err_float);

    tree->Branch("pv_x", &pv_x_float);
    tree->Branch("pv_y", &pv_y_float);
    tree->Branch("pv_z", &pv_z_float);
    tree->Branch("pv_x_err", &pv_x_err_float);
    tree->Branch("pv_y_err", &pv_y_err_float);
    tree->Branch("pv_z_err", &pv_z_err_float);
 //   tree->Branch("PV_width_x", &pv_width_x_float);
 //   tree->Branch("PV_width_y", &pv_width_y_float);

    tree->Branch("v_x", &v_x_float);
    tree->Branch("v_y", &v_y_float);
    tree->Branch("v_z", &v_z_float);

    tree->Branch("dxy_thebeamspot_one", &dxy_thebeamspot_one_float);
    tree->Branch("dxy_thebeamspot_second", &dxy_thebeamspot_second_float);
    tree->Branch("dxy_primaryvertex_one", &dxy_primaryvertex_one_float);
    tree->Branch("dxy_without_argument", &dxy_without_argument_float);

    //FLO muons

    tree->Branch("pt_FLO_vtx", &pt_FLO_vtx_float);
    tree->Branch("eta_FLO_vtx", &eta_FLO_vtx_float);
    tree->Branch("phi_FLO_vtx", &phi_FLO_vtx_float);
    tree->Branch("dxy_BS_FLO_vtx", &dxy_BS_FLO_vtx_float);
    tree->Branch("dxy_PV_FLO_vtx", &dxy_PV_FLO_vtx_float);
    tree->Branch("dxy_BS_FLO_vtx_second_test", &dxy_BS_FLO_vtx_second_test_float);


    tree->Branch("pt_FLO", &pt_FLO_float);
    tree->Branch("eta_FLO", &eta_FLO_float);
    tree->Branch("phi_FLO", &phi_FLO_float);
    tree->Branch("dxy_BS_FLO", &dxy_BS_FLO_float);
    tree->Branch("dxy_PV_FLO", &dxy_PV_FLO_float);
    tree->Branch("charge_FLO", &charge_FLO_float);


    tree->Branch("bs_x_error_base", &bs_x_error_base_float);
    tree->Branch("bs_y_error_base", &bs_y_error_base_float);
    tree->Branch("bs_x_error_set", &bs_x_error_set_float);
    tree->Branch("bs_y_error_set", &bs_y_error_set_float);


}




//define this as a plug-in
DEFINE_FWK_MODULE(VX_BS_Ana);
