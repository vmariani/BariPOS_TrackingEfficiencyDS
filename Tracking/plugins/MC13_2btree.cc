
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h" 
#include "CLHEP/Matrix/Vector.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

#include <map>
#include <string>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Math/VectorUtil.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include <DataFormats/Math/interface/deltaPhi.h>
#include <iostream>

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include <DataFormats/BeamSpot/interface/BeamSpot.h>
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <algorithm>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Version/interface/GetReleaseVersion.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"

#include "TTree.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <DataFormats/Common/interface/Ref.h>
#include <DataFormats/Common/interface/Ptr.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h>
#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "HepMC/GenVertex.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include <TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h>
#include <TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h>
#include <TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/getRef.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "RecoTracker/TransientTrackingRecHit/interface/ProjectedRecHit2D.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TRecHit1DMomConstraint.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TRecHit2DPosConstraint.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiPixelRecHit.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripMatchedRecHit.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripRecHit1D.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripRecHit2DLocalPos.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "RecoVertex/TrimmedKalmanVertexFinder/interface/ConfigurableTrimmedVertexFinder.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"


using namespace std;
using namespace edm;
using namespace reco;


class MC13_2btree : public edm::EDAnalyzer {
 public:
 explicit MC13_2btree(const edm::ParameterSet&);
 ~MC13_2btree();

 static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
 virtual void beginJob() ;
 virtual void analyze(const edm::Event&, const edm::EventSetup&);
 virtual void endJob() ;

 edm::InputTag assocTags_;

 std::map<std::string,TH1F*> histContainer_;

 std::map<std::string,TH2F*> histContainer2_;

 virtual void beginRun(edm::Run const&, edm::EventSetup const&);
 virtual void endRun(edm::Run const&, edm::EventSetup const&);
 virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
 virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
 int i;
  
 edm::EDGetToken tracks_;
 edm::EDGetToken hVtx_;
 edm::EDGetToken beamSpotHandle_;
 edm::EDGetToken genParticles_;
 edm::EDGetToken hPU_;
 
 TTree *Tree;
 TTree *genTree;
 double Pis_pt_gen, Pis_eta_gen, Pis_phi_gen;
 double Pi_pt_gen, Pi_eta_gen, Pi_phi_gen;
 double K_pt_gen, K_eta_gen, K_phi_gen;
 double DS_pt_gen, DS_eta_gen, DS_phi_gen;
 double D0_pt_gen, D0_eta_gen, D0_phi_gen;
 double Pis_pt, Pis_eta, Pis_phi;
 double Pi_pt, Pi_eta, Pi_phi;
 double K_pt, K_eta, K_phi;
 double DS_pt, DS_eta, DS_phi;
 double D0_pt, D0_eta, D0_phi;
 double CL_vertex;
 double x_p, y_p, z_p;
 double x_s, y_s, z_s;
 double L_abs, L_sigma;
 double cos_phi;
 int vertex_size, num_tracks, truePU, pt_tracks;

 // ----------member data ---------------------------
};

MC13_2btree::MC13_2btree(const edm::ParameterSet& iConfig)
{
 edm::Service<TFileService> fs;

 tracks_= consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
 hVtx_= consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"));
 beamSpotHandle_= consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
 genParticles_= consumes<edm::View<reco::GenParticle>>(edm::InputTag("genParticles"));
 hPU_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"));

 genTree = fs->make<TTree>("Gen", "Gen");

 genTree->Branch("Pis_pt_gen",&Pis_pt_gen,"Pis_pt_gen/D");
 genTree->Branch("Pis_eta_gen",&Pis_eta_gen,"Pis_eta_gen/D");
 genTree->Branch("Pis_phi_gen",&Pis_phi_gen,"Pis_phi_gen/D");

 genTree->Branch("Pi_pt_gen",&Pi_pt_gen,"Pi_pt_gen/D");
 genTree->Branch("Pi_eta_gen",&Pi_eta_gen,"Pi_eta_gen/D");
 genTree->Branch("Pi_phi_gen",&Pi_phi_gen,"Pi_phi_gen/D");

 genTree->Branch("K_pt_gen",&K_pt_gen,"K_pt_gen/D");
 genTree->Branch("K_eta_gen",&K_eta_gen,"K_eta_gen/D");
 genTree->Branch("K_phi_gen",&K_phi_gen,"K_phi_gen/D");
 
 genTree->Branch("D0_pt_gen",&D0_pt_gen,"D0_pt_gen/D");
 genTree->Branch("D0_eta_gen",&D0_eta_gen,"D0_eta_gen/D");
 genTree->Branch("D0_phi_gen",&D0_phi_gen,"D0_phi_gen/D");

 genTree->Branch("DS_pt_gen",&DS_pt_gen,"DS_pt_gen/D");
 genTree->Branch("DS_eta_gen",&DS_eta_gen,"DS_eta_gen/D");
 genTree->Branch("DS_phi_gen",&DS_phi_gen,"DS_phi_gen/D");

 Tree= fs->make<TTree>("Analysis","Analysis");

 Tree->Branch("Pis_pt",&Pis_pt,"Pis_pt/D");
 Tree->Branch("Pis_eta",&Pis_eta,"Pis_eta/D");
 Tree->Branch("Pis_phi",&Pis_phi,"Pis_phi/D");
 
 Tree->Branch("Pi_pt",&Pi_pt,"Pi_pt/D");
 Tree->Branch("Pi_eta",&Pi_eta,"Pi_eta/D");
 Tree->Branch("Pi_phi",&Pi_phi,"Pi_phi/D");
 
 Tree->Branch("K_pt",&K_pt,"K_pt/D");
 Tree->Branch("K_eta",&K_eta,"K_eta/D");
 Tree->Branch("K_phi",&K_phi,"K_phi/D");
 
 Tree->Branch("DS_pt",&DS_pt,"DS_pt/D");
 Tree->Branch("DS_eta",&DS_eta,"DS_eta/D");
 Tree->Branch("DS_phi",&DS_phi,"DS_phi/D");
 
 Tree->Branch("D0_pt",&D0_pt,"D0_pt/D");
 Tree->Branch("D0_eta",&D0_eta,"D0_eta/D");
 Tree->Branch("D0_phi",&D0_phi,"D0_phi/D");
 
 Tree->Branch("CL_vertex",&CL_vertex,"CL_vertex/D");
 
 Tree->Branch("x_p",&x_p,"x_p/D");
 Tree->Branch("y_p",&y_p,"y_p/D");
 Tree->Branch("z_p",&z_p,"z_p/D");
 
 Tree->Branch("x_s",&x_s,"x_s/D");
 Tree->Branch("y_s",&y_s,"y_s/D");
 Tree->Branch("z_s",&z_s,"z_s/D");
 
 Tree->Branch("L_abs",&L_abs,"L_abs/D");
 Tree->Branch("L_sigma",&L_sigma,"L_sigma/D");
 Tree->Branch("cos_phi",&cos_phi,"cos_phi/D");
 Tree->Branch("vertex_size",&vertex_size,"vertex_size/I");
 Tree->Branch("pt_tracks",&pt_tracks,"pt_tracks/D");
 Tree->Branch("truePU",&truePU,"truePU/I");
 Tree->Branch("num_tracks",&num_tracks,"num_tracks/I");
}

MC13_2btree::~MC13_2btree()
{
 

}


void
MC13_2btree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
 using namespace edm;

 edm::Handle<edm::View<reco::GenParticle> > genParticles;
 iEvent.getByToken(genParticles_, genParticles);
   
 edm::Handle<reco::TrackCollection> tracks;
 iEvent.getByToken(tracks_ , tracks );

 reco::BeamSpot beamSpot;
 edm::Handle<reco::BeamSpot> beamSpotHandle;
 iEvent.getByToken(beamSpotHandle_ , beamSpotHandle);

 edm::ESHandle<TransientTrackBuilder> theB;
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  
 edm::ESHandle<TransientTrackBuilder> theB1;
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB1);

 int Dstar = 0, D0 = 0, Pislow = 0, Piplus = 0, Kminus = 0;
 double pi_mass = 0.13957018;
 double K_mass = 0.493677;
 double D0_mass = 1.864841;    
 int id = 0, n = 0, m= 0;
  
 TLorentzVector KvectorGen;
 TLorentzVector PivectorGen;
 TLorentzVector D0vectorGen;
 TLorentzVector PisvectorGen;
 TLorentzVector DSvectorGen;
 
		/*########### reconstruction of the chain at gen level ##########*/

 for(size_t i = 0; i < genParticles->size(); ++ i) {
  const GenParticle & p= (*genParticles)[i];
  id = p.pdgId();
  if (id == 413){Dstar ++;
   m = p.numberOfDaughters();
   if(m == 2){
    Pislow = 0;
    D0 = 0;   
    for(int i = 0; i < m; ++ i) {
     const Candidate * d = p.daughter( i );
     id = d->pdgId();
     if(id == 211){
      Pislow++;
      Pis_pt_gen = d->pt();
      Pis_eta_gen = d->eta();
      Pis_phi_gen = d->phi();
      PisvectorGen.SetPtEtaPhiM(Pis_pt_gen, Pis_eta_gen, Pis_phi_gen, pi_mass);
     }
     else if( id == 421 ){
      D0++;
      n = d->numberOfDaughters();
      if(n == 2){
       Piplus=0; Kminus=0; 
       for(int j = 0; j < n; ++ j) {
        const Candidate * e = d->daughter( j );
        int id2 = e->pdgId();
        if(id2 == 211) {
         Piplus++;
         Pi_pt_gen = e->pt();
         Pi_eta_gen = e->eta();
         Pi_phi_gen = e->phi();
         PivectorGen.SetPtEtaPhiM(Pi_pt_gen, Pi_eta_gen, Pi_phi_gen, pi_mass); 
        }
        if(id2 == -321) {
         Kminus++;
         K_pt_gen = e->pt();
         K_eta_gen = e->eta();
         K_phi_gen = e->phi();
         KvectorGen.SetPtEtaPhiM(K_pt_gen, K_eta_gen, K_phi_gen, K_mass);
        } 
       }//for su due corpi
      }//D0 in due corpi           
     }//id D0
    }//for decadim
   }//prodotti decad D*m
  }//D* meno 
 } //for genparticle

 if(Pislow == 1 && D0 == 1){
  if(Kminus == 1 && Piplus == 1){
   D0vectorGen = PivectorGen + KvectorGen;
   D0_pt_gen = D0vectorGen.Pt();
   DSvectorGen = D0vectorGen + PisvectorGen;
   DS_pt_gen = DSvectorGen.Pt();
   genTree->Fill();

		/*##### info vertices ######*/

   edm::Handle<std::vector<reco::Vertex> > hVtx;
   iEvent.getByToken(hVtx_ , hVtx);
   reco::Vertex primVertex;
   reco::Vertex primVertex_tmp;

   if(hVtx->size() > 0){ 

    /*##### beam spot definition ######*/

    double bsx = 0, bsy = 0, bsz = 0;
    if ( beamSpotHandle.isValid() ){ beamSpot = *beamSpotHandle;
     bsx = beamSpot.x0();
     bsy = beamSpot.y0();
     bsz = beamSpot.z0();
    }
    GlobalPoint BeamSpotGP(bsx, bsy, bsz);

    int sameVertex = 0, ntracks = 0, pttracks = 0;
    double inv_mass1 = 0, sigma_L = 0, CL = 0;
    bool goodtracks=false;
    float v_chi2 = 0;
    double deriv[3];
    double cov_sv[3][3];
    double cov_pv[3][3];
    double xprim = 0, yprim = 0, zprim = 0, xsec = 0, ysec = 0, zsec = 0, dx = 0, dy = 0, dz = 0, px = 0, py = 0, pz = 0, pi = 0, cosphi = 0, L = 0;
    bool goodvertex = false;

    for (unsigned int t = 0; t < hVtx->size(); t++){
     primVertex_tmp = hVtx->at(t);
     if(!primVertex_tmp.isFake()  && primVertex_tmp.isValid() && primVertex_tmp.ndof() > 4 && fabs(primVertex_tmp.z()-bsz) < 10 ){
      primVertex = primVertex_tmp;
      goodvertex = true;
      break;
      }
     }

    if(goodvertex){

     reco::TransientTrack  Ktrans;
     reco::TransientTrack  pitrans;
     std::vector<reco::TransientTrack> tks;
     TLorentzVector D0candidates;
     vector<reco::Track> myvector;
     reco::Track tr1;
     reco::Track tr2; 
     reco::Track goodTrack1;
     reco::Track goodTrack2;
     TLorentzVector Kvector;
     TLorentzVector Pivector;
     TLorentzVector D0vector;
     TLorentzVector Pisvector;
     TLorentzVector DSvector;
     TransientVertex v;

				  //######### tracks preselection ##########

     for (reco::TrackCollection::const_iterator track = tracks->begin();  track != tracks->end();  ++track){
      ntracks ++;
      pttracks = pttracks + track->pt();
      if((track->pt() > 0.5) && ( (track->chi2())/(track->ndof()) <= 2.5) && (track->numberOfValidHits() >= 5) && (track->hitPattern().numberOfValidPixelHits() >= 2) && fabs(track->dxy(primVertex.position())) < 0.1 && fabs(track->dz(primVertex.position())) < 1 && (track->quality(Track::highPurity))) { 
       myvector.push_back(*track);	
      }
     }
      			//############ combinatorial #############
     if (myvector.size() > 1){ 
      for (reco::TrackCollection::const_iterator track = myvector.begin();  track != myvector.end();  ++track) {
       for (reco::TrackCollection::const_iterator track1 = track + 1; track1 != myvector.end(); ++track1){
        if (track->charge() + track1->charge() == 0){	
         goodtracks=false;
         inv_mass1 = 0;
         tr1 = *track;
         tr2 = *track1;
         if ((track->charge()== 1) && (track1->charge()==-1)){
          Kvector.SetPtEtaPhiM(track1->pt(), track1->eta(), track1->phi(), K_mass);
          Pivector.SetPtEtaPhiM(track->pt(), track->eta(), track->phi(), pi_mass);
          D0vector = Kvector + Pivector;	
          inv_mass1 = D0vector.M();
          goodtracks = true;
          goodTrack1 = *track;
          goodTrack2 = *track1;
          Ktrans = (*theB).build(tr1);
          pitrans = (*theB).build(tr2);
          }
         else if ((track->charge()==-1) && (track1->charge()== 1) ){
          Kvector.SetPtEtaPhiM(track->pt(), track->eta(), track->phi(), K_mass);
          Pivector.SetPtEtaPhiM(track1->pt(), track1->eta(), track1->phi(), pi_mass);
          D0vector = Kvector + Pivector;
          inv_mass1 = D0vector.M();
          goodtracks = true;
          goodTrack1 = *track;
          goodTrack2 = *track1;
          Ktrans  = (*theB).build(tr1);
          pitrans = (*theB).build(tr2);
         }
         if (fabs(inv_mass1 - D0_mass) < 0.1){
          v_chi2 = 0;
          CL = 0;
          sameVertex = 0;
          tks.clear();
          tks.push_back(Ktrans); 
          tks.push_back(pitrans);
          if (tks.size() > 1){
           KalmanVertexFitter kalman(true);
           if(goodtracks){ v = kalman.vertex(tks); 
            if(v.isValid()){sameVertex++;
             v_chi2 = v.normalisedChiSquared();
             CL = TMath::Prob(v.totalChiSquared(),(int)v.degreesOfFreedom());
            } 
	       }
          }
       
          if(sameVertex > 0 && CL > 0.01 && goodtracks == true){
           xprim = primVertex.position().x();
           yprim = primVertex.position().y();
           zprim = primVertex.position().z();
           xsec = v.position().x();
           ysec = v.position().y();
           zsec = v.position().z();
           dx = xsec - xprim;
           dy = ysec - yprim;
           dz = zsec - zprim;
           L_abs = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
           px = D0vector.Px();
           py = D0vector.Py();
           pz = D0vector.Pz();
           pi = D0vector.P();
           cosphi = (px*dx + py*dy + pz*dz)/(L_abs*pi);
           L = (px*dx + py*dy + pz*dz)/pi;

           if (cosphi > 0.9){
           //covariance matrix SV
            cov_sv[0][0] = v.positionError().cxx();
            cov_sv[1][0] = v.positionError().cyx();
            cov_sv[2][0] = v.positionError().czx();
            cov_sv[0][1] = cov_sv[1][0];
            cov_sv[1][1] = v.positionError().cyy();
            cov_sv[2][1] = v.positionError().czy();
            cov_sv[0][2] = cov_sv[2][0];
            cov_sv[1][2] = cov_sv[2][1];
            cov_sv[2][2] = v.positionError().czz();
        
           //covariance matrix PV
            cov_pv[0][0] = primVertex.covariance(0,0);
            cov_pv[1][0] = primVertex.covariance(1,0);
            cov_pv[2][0] = primVertex.covariance(2,0);
            cov_pv[0][1] = cov_pv[1][0];
            cov_pv[1][1] = primVertex.covariance(1,1);
            cov_pv[2][1] = primVertex.covariance(2,1);
            cov_pv[0][2] = cov_pv[2][0];
            cov_pv[1][2] = cov_pv[2][1];
            cov_pv[2][2] = primVertex.covariance(2,2);    
            deriv[0] = dx/L_abs;
            deriv[1] = dy/L_abs;
            deriv[2] = dz/L_abs;

            if(D0vector.Pt() > 3){

              // third track selection 

             for(reco::TrackCollection::const_iterator track3 = tracks->begin();  track3 != tracks->end();  ++track3){ 
              if(track3->charge() == 1 && (track3->pt() > 0.3) && (track3->numberOfValidHits() >= 2) && fabs(track3->dxy(primVertex.position())/track3->dxyError()) < 3 && fabs(track3->dz(primVertex.position())/track3->dzError()) < 3 && ((track3->chi2())/(track3->ndof()) <= 3) ) {   
               if(track3->pt() != Pivector.Pt() && track3->eta() != Pivector.Eta()){          
	        Pisvector.SetPtEtaPhiM(track3->pt(), track3->eta(), track3->phi(), pi_mass);
                Pis_pt = Pisvector.Pt();
                DSvector = Pisvector + Kvector + Pivector;

                if (fabs(DSvector.M() - D0vector.M()) < 0.160){
                 if(DSvector.Pt() > 4.){
               
                  Pis_pt=Pisvector.Pt();
                  Pis_eta=Pisvector.Eta();
                  Pis_phi=Pisvector.Phi();
                  Pi_pt=Pivector.Pt();
                  Pi_eta=Pivector.Eta();
                  Pi_phi=Pivector.Phi();

                  K_pt=Kvector.Pt();
                  K_eta=Kvector.Eta();
                  K_phi=Kvector.Phi();
                  DS_pt=DSvector.Pt();
                  DS_eta=DSvector.Eta();
                  DS_phi=DSvector.Phi();
                  D0_pt=D0vector.Pt();
                  D0_eta=D0vector.Eta();
                  D0_phi=D0vector.Phi();

                  CL_vertex = CL;
                  vertex_size = hVtx->size();
                  cos_phi = (px*dx + py*dy + pz*dz)/(L_abs*pi);
                  x_p = xprim;
                  y_p = yprim;
                  z_p = zprim;

                  x_s = xsec;
                  y_s = ysec;
                  z_s = zsec;
                  L = (px*dx + py*dy + pz*dz)/pi;
                  sigma_L = 0;
                  for (int m = 0; m < 3; ++m){
                   for (int n = 0; n < 3; ++n ){
                    sigma_L += deriv[m]*deriv[n]*(cov_pv[m][n] + cov_sv[m][n]);
                   }
                  }
                  sigma_L = sqrt(sigma_L);
                  L_sigma = L/sigma_L;
                  num_tracks = ntracks;
                  pt_tracks = pttracks;

                  //filling tree

                  Tree->Fill();                 
                 
                  } //pt dstar
                 }//diff < 160			 
                }//pislow != pi
               } // third track
              }//for 
             }//D0 pt
            }//coseno
           }// same vertex 
          }//D0 mass
         }//if total charge 0
        }//for track 1 
       }// for track
      }//myvector >1
     }//good vertex
   }//D0
  } //D* 
 }//#vertices > 0
}

// ------------ method called once each job just before starting event loop  ------------
void 
MC13_2btree::beginJob()
{

edm::Service<TFileService> fs;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MC13_2btree::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MC13_2btree::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MC13_2btree::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MC13_2btree::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MC13_2btree::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MC13_2btree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MC13_2btree);
