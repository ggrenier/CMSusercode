#ifndef MyRPCSimHitAnalyzer_HH
#define MyRPCSimHitAnalyzer_HH


// -*- C++ -*-
//
// Package:    MyRPCSimHitAnalyzer
// Class:      MyRPCSimHitAnalyzer
// 
/**\class MyRPCSimHitAnalyzer MyRPCSimHitAnalyzer.cc testGG/MyRPCSimHitAnalyzer/src/MyRPCSimHitAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gerald Grenier,40 4-B16,+41227671558,
//         Created:  Tue Dec 10 12:06:11 CET 2013
// $Id$
//
//




// system include files
#include <memory>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

//Data format
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

//Geometry stuff
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"

#include <Geometry/RPCGeometry/interface/RPCRoll.h>

#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/RPC_MCNavigator.h"

//For plots and histos
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TH1F.h>
#include <TRandom3.h>
#include <TTree.h>

//local for plots and histos
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/RPC_MultiHisto.h"
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/HistoBooker_TFileService.h"

//Local classes
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/ZvertexCalculator.h"
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/ZvertexCalculatorMulti.h"
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/hitDataPair.h"
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/RPC_HitPair_MultiHist.h"




//
// class declaration
//

class MyRPCSimHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyRPCSimHitAnalyzer(const edm::ParameterSet&);
      ~MyRPCSimHitAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      //virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  const BoundPlane & SimHitSurface(const PSimHit& hit);
  GlobalPoint SimHitGlobalPosition(const PSimHit& hit) { return SimHitSurface(hit).toGlobal(hit.entryPoint());}
  double HitRefTime(const PSimHit& hit) {return SimHitGlobalPosition(hit).mag()/RPC::cspeed();}

      // ----------member data ---------------------------

  //histograms
  TH1F *h_nhits;
  TH1F *h_detId, *h_subdetId; //checking histograms
  TH1F *h_detRegion,*h_detStation,*h_detRing;
  TH1F *h_detSector,*h_detSubsector,*h_detRoll;
  TH1F *h_nonMuonHitEnergy,*h_muonHitEnergy;
  TH1F *h_NPVVertex;
  TH1F *h_PVtime;
  TH2F *h_PVTimeVsZpos;
  TH2F *h_GenParticleVertex,*h_GenParticleVertexZoom;
  TH2F *h_GenPartVertexVsHitTimeVertex,*h_GenPartVertexVsHitTimeTriangleVertex;

  RPC_MultiHisto mh_hitPdgId;
  RPC_MultiHisto mh_tofMuon;
  RPC_MultiHisto mh_tofOther;
  enum {NtimeAxes=3};
  RPC_MultiHisto mh_tofRefMuon[NtimeAxes];
  RPC_MultiHisto mh_tofRefOther[NtimeAxes];
  RPC_MultiHisto mh_ZvtxMuon;
  RPC_MultiHisto mh_ZvtxOther;
  RPC_MultiHisto mh_ZvtxDiffMuon;
  RPC_MultiHisto mh_ZvtxDiffOther;

  RPC_HitPair_MultiHist mh_hp_simlevel;
  RPC_HitPair_MultiHist mh_hp_smearT10ps;
  RPC_HitPair_MultiHist mh_hp_smearT50ps;
  RPC_HitPair_MultiHist mh_hp_smearT100ps;
  RPC_HitPair_MultiHist mh_hp_smearT1ns;
  RPC_HitPair_MultiHist mh_hp_smearXY2mm;

  //tree
  TTree *t_francois;
  //variable for the tree branches
  std::vector<int> pidvec,regionvec,stationvec,ringvec;
  std::vector<float> xvec,yvec,zvec,tvec;
  std::vector<float> localdirEtavec,localdirPhivec,localdirThetavec;
  double treePVZ,treePVCT;
  std::vector<float> muonPxvec,muonPyvec,muonPzvec;
  void addHitData(int pidval,hitData &d);
  void clearVec();

  //configurable parameters
  edm::InputTag cfg_RPCSimHitCollection;

  //geometry accessor to be loaded at beginRun
  edm::ESHandle <RPCGeometry> rpcGeo;

  //Access to data
  edm::EDGetTokenT<edm::PSimHitContainer> _theSimHitsToken;
  edm::EDGetTokenT<std::vector<SimTrack> > _simTracksToken;
  edm::EDGetTokenT<std::vector<SimVertex> > _simVerticesToken;



};

#endif
