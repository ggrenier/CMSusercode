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









class hitDataPair
{
public:
  hitDataPair(hitData &hit1, hitData &hit2)
    : d1(hit1), d2(hit2)
  {
    ComputeDistance();
    ComputeTime();
    std::cout << "time in ns (before,after,expect)=("
	      << timeBefore << "," << timeAfter << "," << timeExpect << ") ";
  }

  bool diffEncap() { return d1.region != d2.region;}
  bool lightCone() { return timeAfter==timeExpect;}
  bool lightCone(double timeError) { return fabs(timeAfter-timeExpect)<timeError;}

  ZvertexCalculatorMulti* getVertex()
  {
    return FindZFromPair(d1,d2);
  }

  int station1() {return d1.station;}
  int station2() {return d2.station;}
  
private:
  hitData d1;
  hitData d2;
  double distance;
  double timeBefore;
  double timeAfter;
  double timeExpect;

  void ComputeDistance()
  {
    distance=0;
    distance+=(d1.x-d2.x)*(d1.x-d2.x);
    distance+=(d1.y-d2.y)*(d1.y-d2.y);
    distance+=(d1.z-d2.z)*(d1.z-d2.z);
    distance=sqrt(distance);
  }
  void ComputeTime()
  {
    timeBefore=0;
    timeAfter=0;
    if (d1.t > d2.t)
      {
	timeBefore=d2.t;
	timeAfter=d1.t;
      }
    else
      {
	timeBefore=d1.t;
	timeAfter=d2.t;
      }
    timeExpect=timeBefore+distance/RPC::cspeed();
  }

};





//
// Helper class for plotting histos
//
class RPC_HitPair_MultiHist
{
public:
  enum {DiffEC,LightCone,LightCone10ps,LightCone100ps,LightCone200ps,LightCone1ns,LightCone1BC,Ncases} CASES;
  RPC_HitPair_MultiHist(std::string extensionString="");
  void fill(hitDataPair &hp,float& PVZ);
  void fill(std::vector<hitData>& vechit,float PVZ);
private:
  std::string extension;
  bool pass(unsigned int i,hitDataPair &paire);
  TH1F* statHist;
  TH1F* hist[Ncases];
  TH1F* histSameSideStationPair[(Ncases-1)][6]; //station pair = (1,2) (1,3) (1,4) (2,3) (2,4) (3,4)
  const char* extend(std::string val) {return (val+"_"+extension).c_str();}
  std::string histoTitle(std::string precision) {std::string s("Estimated PV Z - true PV Z (in cm) "); return s+precision;}
  int indice(int station1, int station2);
};

bool RPC_HitPair_MultiHist::pass(unsigned int i,hitDataPair &paire)
{
  switch (i)
    {
    case 0:
      return paire.diffEncap();
    case 1:
      return paire.lightCone();
    case 2:
      return paire.lightCone(0.01);
    case 3:
      return paire.lightCone(0.1);
    case 4:
      return paire.lightCone(0.2);
    case 5:
      return paire.lightCone(1.0);
    case 6:
      return paire.lightCone(25.0);
    }
  return false;
}

RPC_HitPair_MultiHist::RPC_HitPair_MultiHist(std::string extensionString) : extension(extensionString)
{
  std::string names[Ncases]={"DiffEC","LightCone","LightCone10ps","LightCone100ps","LightCone200ps","LightCone1ns","LightCone1BC"};
  edm::Service<TFileService> fs;
  statHist = fs->make<TH1F>(extend("hitPairStatus"),extend("hit pairs quality"),Ncases+1,0,Ncases+2);
  statHist->GetXaxis()->SetBinLabel(1, "All");
  int nbins=1000;
  float xmin=-100;
  float xmax=100;
  for (unsigned int i=0; i<Ncases; i++) 
    {
      statHist->GetXaxis()->SetBinLabel(i+2,names[i].c_str());
      hist[i]=fs->make<TH1F>(extend(names[i]),extend(histoTitle(names[i])),nbins,xmin,xmax);
      if (i==0) continue;
      histSameSideStationPair[i-1][0]=fs->make<TH1F>(extend(names[i]+"_12"),extend(histoTitle(names[i]+" sections 1 and 2")),nbins,xmin,xmax);
      histSameSideStationPair[i-1][1]=fs->make<TH1F>(extend(names[i]+"_13"),extend(histoTitle(names[i]+" sections 1 and 3")),nbins,xmin,xmax);
      histSameSideStationPair[i-1][2]=fs->make<TH1F>(extend(names[i]+"_14"),extend(histoTitle(names[i]+" sections 1 and 4")),nbins,xmin,xmax);
      histSameSideStationPair[i-1][3]=fs->make<TH1F>(extend(names[i]+"_23"),extend(histoTitle(names[i]+" sections 2 and 3")),nbins,xmin,xmax);
      histSameSideStationPair[i-1][4]=fs->make<TH1F>(extend(names[i]+"_24"),extend(histoTitle(names[i]+" sections 2 and 4")),nbins,xmin,xmax);
      histSameSideStationPair[i-1][5]=fs->make<TH1F>(extend(names[i]+"_34"),extend(histoTitle(names[i]+" sections 3 and 4")),nbins,xmin,xmax);
    }
}

int RPC_HitPair_MultiHist::indice(int station1, int station2)
{
  if (station1==1) return station2-2;
  if (station2==1) return station1-2;
  //should be left with (2,3) (2,4) and (3,4)
  if (station1==2) return station2;
  if (station2==2) return station1;
  //should be left with (3,4)
  return 5;
}

void RPC_HitPair_MultiHist::fill(hitDataPair &hp,float& PVZ)
{
  statHist->Fill(1);
  ZvertexCalculatorMulti* Z=hp.getVertex();
  for (unsigned int i=0; i<Ncases; i++)
    {
      if (pass(i,hp)) 
	{
	  statHist->Fill(i+2);
	  if (NULL != Z)
	    {
	      float diffZvtx=Z->meanSolZ()-PVZ;
	      hist[i]->Fill(diffZvtx);
	      if (i!=0) histSameSideStationPair[i-1][indice(hp.station1(),hp.station2())]->Fill(diffZvtx);
	    }
	}
    }
  if (NULL != Z) delete Z;
}

void  RPC_HitPair_MultiHist::fill(std::vector<hitData>& vechit,float PVZ)
{
  for (std::vector<hitData>::iterator it1=vechit.begin(); it1 != vechit.end(); it1++)
    {
      std::vector<hitData>::iterator it2=it1;
      it2++;
      for (; it2 != vechit.end(); it2++)
	{
	  hitDataPair HP(*it1,*it2);
	  fill(HP,PVZ);
	}
    }
}

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
  void addHitData(int pidval,hitData &d)
  {
    pidvec.push_back(pidval);
    regionvec.push_back(d.region);
    stationvec.push_back(d.station);
    ringvec.push_back(d.ring);
    xvec.push_back(d.x);
    yvec.push_back(d.y);
    zvec.push_back(d.z);
    tvec.push_back(d.t);
    localdirEtavec.push_back(d.localdirEta);
    localdirPhivec.push_back(d.localdirPhi);
    localdirThetavec.push_back(d.localdirTheta);
  }
  void clearVec()
  {
    pidvec.clear();
    regionvec.clear();
    stationvec.clear();
    ringvec.clear();
    xvec.clear();
    yvec.clear();
    zvec.clear();
    tvec.clear();
    localdirEtavec.clear();
    localdirPhivec.clear();
    localdirThetavec.clear();
    muonPxvec.clear();
    muonPyvec.clear();
    muonPzvec.clear();
  }

  //configurable parameters
  edm::InputTag cfg_RPCSimHitCollection;

  //geometry accessor to be loaded at beginRun
  edm::ESHandle <RPCGeometry> rpcGeo;

  //Access to data
  edm::EDGetTokenT<edm::PSimHitContainer> _theSimHitsToken;
  edm::EDGetTokenT<std::vector<SimTrack> > _simTracksToken;
  edm::EDGetTokenT<std::vector<SimVertex> > _simVerticesToken;



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
MyRPCSimHitAnalyzer::MyRPCSimHitAnalyzer(const edm::ParameterSet& iConfig) :
  mh_hp_simlevel(),
  mh_hp_smearT10ps("smear10ps"),
  mh_hp_smearT50ps("smear50ps"),
  mh_hp_smearT100ps("smear100ps"),
  mh_hp_smearT1ns("smear1ns"),
  mh_hp_smearXY2mm("smearXY2mm"),
  cfg_RPCSimHitCollection(iConfig.getUntrackedParameter<edm::InputTag>("RPCSimHitCollection",edm::InputTag("g4SimHits","MuonRPCHits","SIM"))),
  //_theSimHitsToken(consumes<edm::PSimHitContainer >( edm::InputTag("g4SimHits","MuonRPCHits","SIM") )),
  _theSimHitsToken(consumes<edm::PSimHitContainer >( cfg_RPCSimHitCollection )),
  _simTracksToken(consumes<std::vector<SimTrack> >(edm::InputTag(cfg_RPCSimHitCollection.label(),"",cfg_RPCSimHitCollection.process()))),
  _simVerticesToken(consumes<std::vector<SimVertex> >(edm::InputTag(cfg_RPCSimHitCollection.label(),"",cfg_RPCSimHitCollection.process())))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  h_nhits = fs->make<TH1F>("Nmuonhits","Number of muon hits",100,0,100);
  h_detId = fs->make<TH1F>("muonHitDetId","Detector Id { Tracker=1,Muon=2,Ecal=3,Hcal=4,Calo=5 }",5,1,6);
  h_subdetId = fs->make<TH1F>("muonHitSubDetId","Muon Detector Sub Id  { DT=1,CSC=2,RPC=3,GEM=4 }",4,1,5);
  h_detRegion = fs->make<TH1F>("muonHitRPCRegion","RPC hit region { 0 for Barrel, +/-1 For +/- Endcap }",3,-1,2);
  h_detStation = fs->make<TH1F>("muonHitECStation","RPC hit EC station",4,1,5);
  h_detRing = fs->make<TH1F>("muonHitECRing","RPC hit EC ring",3,1,4);
  h_detSector = fs->make<TH1F>("muonHitECSector","RPC hit EC sector",6,1,7);
  h_detSubsector = fs->make<TH1F>("muonHitECSubsector","RPC hit EC sub sector",6,1,7);
  h_detRoll = fs->make<TH1F>("muonHitECRoll","RPC hit EC roll (0=wild card)",5,0,5);

  h_nonMuonHitEnergy=fs->make<TH1F>("nonMuonHitEnergy","RPC hit MC part log energy in GeV for part other than muon",100,-10,4);
  h_muonHitEnergy=fs->make<TH1F>("muonHitEnergy","RPC hit MC part log energy in GeV for true muon",100,-10,4);
  h_NPVVertex=fs->make<TH1F>("NPV","Number of found Primary Vertices",100,0,100);
  h_PVtime=fs->make<TH1F>("PVTime","Time of Primary Vertices in ns",200,-1,1);
  h_PVTimeVsZpos=fs->make<TH2F>("PVTimeVsZ","Primary Vertices : time in ns vs Z in cm",100,-30,30,200,-1,1);
  h_GenParticleVertex=fs->make<TH2F>("genVertex","Vertex r vs z for generated particles in cm",100,-30,30,50,0,15);
  h_GenParticleVertexZoom=fs->make<TH2F>("genVertexZoom","Vertex r vs z for generated particles in cm",100,-1,1,50,0,1);
  h_GenPartVertexVsHitTimeVertex=fs->make<TH2F>("genVertexVsHitVertex","Z vertices reconstructed with hit time vs Z vertex in generator (in cm)",100,-30,30,100,-30,30);
  h_GenPartVertexVsHitTimeTriangleVertex=fs->make<TH2F>("genVertexVsHitVertexTriangle","Z vertices reconstructed with hit time using triangle relation vs Z vertex in generator (in cm)",100,-30,30,100,-30,30);

  HistoBooker_TFileService hbooker;
  //Don't call RPC_MultiHisto book function outside this constructor (hbooker pointer will be destroyed)
  mh_hitPdgId.setHistoBooker(&hbooker);
  mh_hitPdgId.book("muonHitPDGid","RPC hit EC PDG ID",20,0,20);
  mh_tofMuon.setHistoBooker(&hbooker);
  mh_tofMuon.book("muonHitTOFtrueMuon","TOF (ns) of EC RPC hits produced by muons",100,0,100);
  mh_tofOther.setHistoBooker(&hbooker);
  mh_tofOther.book("muonHitTOFtrueNotMuon","TOF (ns) of EC RPC hits not produced by muons",100,0,100);


  HistoBooker::Axis axes[NtimeAxes]={HistoBooker::Axis(200,-0.1,0.1),
				     HistoBooker::Axis(200,-1,1),
				     HistoBooker::Axis(200,-100,100)};
  std::string namemuon="RPChitTOFDeltaTrueMuon_";
  std::string namenotmuon="RPChitTOFDeltaTrueNotMuon_";
  for (unsigned int i=0; i<NtimeAxes; i++)
    {
      namemuon+='Z';
      mh_tofRefMuon[i].setHistoBooker(&hbooker);
      mh_tofRefMuon[i].RPC_MultiHistoBase::book(namemuon,"TOF - Straight r/c (ns) of EC RPC hits produced by muons",axes[i]);
      namenotmuon+='Z';
      mh_tofRefOther[i].setHistoBooker(&hbooker);
      mh_tofRefOther[i].RPC_MultiHistoBase::book(namenotmuon,"TOF - Straight r/c (ns) of EC RPC hits not produced by muons",axes[i]);
    }

  mh_ZvtxMuon.setHistoBooker(&hbooker);
  mh_ZvtxMuon.book("EstZVtxPos_muonHitTrueMuon","Estimated position of Z vertex in cm for EC RPC hits produced by muons",300,-30,30);
  mh_ZvtxOther.setHistoBooker(&hbooker);
  mh_ZvtxOther.book("EstZVtxPos_muonHitTrueNotMuon","Estimated position of Z vertex in cm for EC RPC hits not produced by muons",300,-30,30);
  mh_ZvtxDiffMuon.setHistoBooker(&hbooker);
  mh_ZvtxDiffMuon.book("EstZVtxDiffPos_muonHitTrueMuon","Estimated position of Z vertex minus real PV Z in cm for EC RPC hits produced by muons",200,-10,10);
  mh_ZvtxDiffOther.setHistoBooker(&hbooker);
  mh_ZvtxDiffOther.book("EstZVtxDiffPos_muonHitTrueNotMuon","Estimated position of Z vertex minus real PV Z in cm for EC RPC hits not produced by muons",200,-10,10);


  t_francois=fs->make<TTree>("HitTree","CMS RPC hits");
  t_francois->Branch("PVZ",&treePVZ);
  t_francois->Branch("PVCT",&treePVCT);
  t_francois->Branch("pid",&pidvec);
  t_francois->Branch("region",&regionvec);
  t_francois->Branch("station",&stationvec);
  t_francois->Branch("ring",&ringvec);
  t_francois->Branch("x",&xvec);
  t_francois->Branch("y",&yvec);
  t_francois->Branch("z",&zvec);
  t_francois->Branch("t",&tvec);
  t_francois->Branch("localdirEta",&localdirEtavec);
  t_francois->Branch("localdirPhi",&localdirPhivec);
  t_francois->Branch("localdirTheta",&localdirThetavec);
  t_francois->Branch("genMuonPx",&muonPxvec);
  t_francois->Branch("genMuonPy",&muonPyvec);
  t_francois->Branch("genMuonPz",&muonPzvec);


}


MyRPCSimHitAnalyzer::~MyRPCSimHitAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

const BoundPlane & MyRPCSimHitAnalyzer::SimHitSurface(const PSimHit& hit)
{
  RPCDetId rollId(hit.detUnitId());  
  RPCGeomServ rpcsrv(rollId);
  const RPCRoll* rollasociated = rpcGeo->roll(rollId);
  if (0==rollasociated)
    {
      std::cout << "Sorry, I have a problem with the detector geometry....ABORTING" << std::endl;
      abort();
    }
  return rollasociated->surface();
}


// ------------ method called for each event  ------------
void
MyRPCSimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //InputTag tag("g4SimHits","MuonRPCHits","SIM");
    edm::Handle<edm::PSimHitContainer> theSimHitsHandle;
    //iEvent.getByLabel(cfg_RPCSimHitCollection,theSimHitsHandle);
    iEvent.getByToken(_theSimHitsToken, theSimHitsHandle);
    InputTag NoInstanceTag(cfg_RPCSimHitCollection.label(),"",cfg_RPCSimHitCollection.process());

    //access to MC information
    edm::Handle<std::vector<SimTrack> > simTracksHandle;
    edm::Handle<std::vector<SimVertex> > simVerticesHandle;
    //iEvent.getByLabel(NoInstanceTag,simTracksHandle);
    //iEvent.getByLabel(NoInstanceTag,simVerticesHandle);
    iEvent.getByToken(_simTracksToken,simTracksHandle);
    iEvent.getByToken(_simVerticesToken,simVerticesHandle);
    std::set<unsigned int> PV_indices;
    std::set<const SimTrack*> generatedMuons;
    if (simTracksHandle.isValid())
      {
	const std::vector<SimTrack> &simTracks=*simTracksHandle;
	for (std::vector<SimTrack>::const_iterator itt=simTracks.begin(); itt != simTracks.end(); itt++)
	  {
	    if (! itt->noGenpart() && ! itt->noVertex()) PV_indices.insert(itt->vertIndex());
	    if (! itt->noGenpart()) generatedMuons.insert(&(*itt));
	  }
      }
    std::set<const SimVertex*> PrimaryVertices;
    if (simVerticesHandle.isValid())
      {
	const std::vector<SimVertex> &simVertices=*simVerticesHandle;
	for (std::set<unsigned int>::iterator iti=PV_indices.begin(); iti != PV_indices.end(); iti++)
	  {
	    for (std::vector<SimVertex>::const_iterator itv=simVertices.begin(); itv!=simVertices.end(); itv++)
	      {
		if (itv->vertexId()==*iti) {PrimaryVertices.insert(&(*itv)); break;}
	      }
	  }
      }
    for (std::set<const SimVertex*>::iterator itpv=PrimaryVertices.begin(); itpv != PrimaryVertices.end(); itpv++)
      {
	const SimVertex &vertex=*(*itpv);
	h_GenParticleVertex->Fill(vertex.position().Z(),vertex.position().Pt());
	h_GenParticleVertexZoom->Fill(vertex.position().Z(),vertex.position().Pt());
	h_PVtime->Fill(vertex.position().T()*1e9); // time is stored in second
	h_PVTimeVsZpos->Fill(vertex.position().Z(),vertex.position().T()*1e9); 
      }
    double PVZ=3000; //in cm
    double PVCT=3e9; //in cm
    if (! PrimaryVertices.empty()) 
      {
	PVZ=(*PrimaryVertices.begin())->position().Z();
	PVCT=(*PrimaryVertices.begin())->position().T()*1e9*RPC::cspeed();
      }
    h_NPVVertex->Fill(PrimaryVertices.size());
    std::vector<hitData> allMuonHits;
    if (theSimHitsHandle.isValid())
      {
	const edm::PSimHitContainer &theSimHits=*theSimHitsHandle; //this is a std::vector<PSimHit> type
	std::cout << "Got " << theSimHits.size() << " hits" << std::endl;
	h_nhits->Fill(theSimHits.size());

	//loop on PSimHit
	for ( edm::PSimHitContainer::const_iterator iHit=theSimHits.begin(); iHit != theSimHits.end(); iHit++)
	  {
	    //Detector Id
	    DetId simdetid= DetId((*iHit).detUnitId());
	    h_detId->Fill(int(simdetid.det()));
	    if (simdetid.det()==DetId::Muon)  //Muon detector
	      {
		h_subdetId->Fill(simdetid.subdetId());
		if (simdetid.subdetId()==MuonSubdetId::RPC)
		  {
		    RPCDetId rollId(simdetid);
		    int region=rollId.region();
		    h_detRegion->Fill(region);
		    if (region != 0) //RPC endcap
		      {
			int station=rollId.station();
			int ring=rollId.ring();
			h_detStation->Fill(station);
			h_detRing->Fill(ring);
			h_detSector->Fill(rollId.sector());
			h_detSubsector->Fill(rollId.subsector());
			h_detRoll->Fill(rollId.roll());
			mh_hitPdgId.Fill(region,station,ring,abs((*iHit).particleType()));
			if (abs((*iHit).particleType()) == 13 ) mh_tofMuon.Fill(region,station,ring,(*iHit).timeOfFlight());
			else mh_tofOther.Fill(region,station,ring,(*iHit).timeOfFlight());
			double diffreftime=(*iHit).timeOfFlight()-HitRefTime(*iHit);
			for (unsigned int i=0; i<NtimeAxes; i++)
			  {
			    if (abs((*iHit).particleType()) == 13 ) mh_tofRefMuon[i].Fill(region,station,ring,diffreftime);
			    else mh_tofRefOther[i].Fill(region,station,ring,diffreftime);
			  }
			double EstZVertex= (SimHitGlobalPosition(*iHit).mag()-RPC::cspeed()*(*iHit).timeOfFlight())*region;
			hitData d; 
			d.x=SimHitGlobalPosition(*iHit).x();
			d.y=SimHitGlobalPosition(*iHit).y();
			d.r=SimHitGlobalPosition(*iHit).mag();
			d.z=SimHitGlobalPosition(*iHit).z();
			d.t=(*iHit).timeOfFlight();
			d.region=region;
			d.station=station;
			d.ring=ring;
			LocalVector LV=(*iHit).momentumAtEntry();
			d.localdirEta=LV.eta();
			d.localdirPhi=LV.phi();
			d.localdirTheta=LV.theta();
			ZvertexCalculator zvc(d);
			addHitData((*iHit).particleType(),d);
			if (abs((*iHit).particleType()) == 13 ) 
			  {
			    allMuonHits.push_back(d);
			    mh_ZvtxMuon.Fill(region,station,ring,zvc.aMinSolution());
			    mh_ZvtxDiffMuon.Fill(region,station,ring,zvc.aMinSolution()-PVZ);
			    h_GenPartVertexVsHitTimeVertex->Fill(EstZVertex,PVZ);
			    h_GenPartVertexVsHitTimeTriangleVertex->Fill(zvc.aMinSolution(),PVZ);
			    std::cout << "("<<SimHitGlobalPosition(*iHit).x()<<","<<SimHitGlobalPosition(*iHit).y()<<","<< SimHitGlobalPosition(*iHit).z()<<"cm,"<<(*iHit).timeOfFlight()<<"ns) ";
			    std::cout << "Vertex (z,ct)=("<<PVZ<<","<<PVCT<<") ";
			    std::cout << "ZVC : " << zvc.hasSolution() << " " << zvc.aMinSolution() << " " << EstZVertex << std::endl;
			  }
			else 
			  {
			    mh_ZvtxOther.Fill(region,station,ring,zvc.aMinSolution());
			    mh_ZvtxDiffOther.Fill(region,station,ring,zvc.aMinSolution()-PVZ);
			  }
			//access to MC tree
			if (simTracksHandle.isValid() && simVerticesHandle.isValid())
			  {
			    const std::vector<SimTrack> &simTracks=*simTracksHandle;
			    const std::vector<SimVertex> &simVertices=*simVerticesHandle;
			    
			    std::vector<SimTrack>::const_iterator theSimTrack=RPC::getTrack(simTracks,*iHit);
			    if (theSimTrack != simTracks.end())
			      {
				std::cout << "SimTrack : " << *theSimTrack << " parent "; 			    
				const SimTrack& parentTrack=RPC::getParentTrack(simTracks,simVertices,*theSimTrack);
				if (&(*theSimTrack)==&parentTrack) std::cout << "none";
				else  std::cout << parentTrack;
				std::cout << std::endl;
				if (abs((*iHit).particleType()) == 13 ) 
				  h_muonHitEnergy->Fill(log10(theSimTrack->momentum().energy()));
				else
				  h_nonMuonHitEnergy->Fill(log10(theSimTrack->momentum().energy()));
			      }
			    else //(theSimTrack != simTracks.end())
			      std::cout << "PSimHit without SimTrack found" << std::endl;

			  } // simTracksHandle.isValid() && simVerticesHandle.isValid()
				 //std::cout << " Reading the Roll"<<std::endl;

			//const RPCRoll* rollasociated = rpcGeo->roll(rollId);
			//const BoundPlane & RPCSurface = rollasociated->surface(); 
			//GlobalPoint SimHitInGlobal = RPCSurface.toGlobal((*iHit).entryPoint());


		      } //rollId.region() != 0) //RPC endcap
		  } //simdetid.subdetId()==MuonSubdetId::RPC
	      } //simdetid.det()==DetId::Muon
	  } //end loop on PSimHit
	treePVZ=PVZ;
	treePVCT=PVCT;
	for (std::set<const SimTrack*>::const_iterator itgenmu=generatedMuons.begin(); itgenmu!=generatedMuons.end(); ++itgenmu)
	  {
	    const SimTrack &track=**itgenmu;
	    if (abs(track.type())==13)
	      {
		muonPxvec.push_back(track.momentum().Px());
		muonPyvec.push_back(track.momentum().Py());
		muonPzvec.push_back(track.momentum().Pz());
	      }
	  }
	t_francois->Fill();
	clearVec();
      } //theSimHitsHandle.isValid()
    else 
      std::cerr << "Handle is not valid" << std::endl;

    //look at muon hits:
    int count1=0;
    for (std::vector<hitData>::iterator it1=allMuonHits.begin(); it1 != allMuonHits.end(); it1++)
      {
	count1++; int count2=count1;
	std::vector<hitData>::iterator it2=it1;
	it2++;
	for (; it2 != allMuonHits.end(); it2++)
	  {
            count2++;
            std::cout << "Treating "<<count1<<":"<<count2<<" over "<<allMuonHits.size() << std::endl;
	    
	    hitDataPair HP(*it1,*it2);
	    if (HP.diffEncap()) std::cout << "diffEC ";
	    if (HP.lightCone()) std::cout << "lightCone ";
	    if (HP.lightCone(0.01)) std::cout << "lightCone10ps ";
	    ZvertexCalculatorMulti* Z=HP.getVertex();
	    if (NULL==Z)
	      std::cout << "No solution found for Z"<< std::endl;
	    else
	      {
		std::cout << " Zvtx ( T=" << Z->getTimeShift() << "," << " Z=" << Z->meanSolZ() << ") " << std::endl; 
		delete Z;
	      }
	  }//iterator it2
      } //iterator it1
    if (PVZ != 3000)
      mh_hp_simlevel.fill(allMuonHits,PVZ);

    TRandom3 rnd;
    std::vector<hitData> smeared=allMuonHits;
    for (std::vector<hitData>::iterator it=smeared.begin(); it != smeared.end(); it++)
      {
	it->t+=rnd.Gaus(0,0.01);
      }
    mh_hp_smearT10ps.fill(smeared,PVZ);
    
    smeared=allMuonHits;
    for (std::vector<hitData>::iterator it=smeared.begin(); it != smeared.end(); it++)
      {
	it->t+=rnd.Gaus(0,0.05);
      }
    mh_hp_smearT50ps.fill(smeared,PVZ);

    smeared=allMuonHits;
    for (std::vector<hitData>::iterator it=smeared.begin(); it != smeared.end(); it++)
      {
	it->t+=rnd.Gaus(0,0.1);
      }
    mh_hp_smearT100ps.fill(smeared,PVZ);

    smeared=allMuonHits;
    for (std::vector<hitData>::iterator it=smeared.begin(); it != smeared.end(); it++)
      {
	it->t+=rnd.Gaus(0,1.);
      }
    mh_hp_smearT1ns.fill(smeared,PVZ);

    smeared=allMuonHits;
    for (std::vector<hitData>::iterator it=smeared.begin(); it != smeared.end(); it++)
      {
	it->x+=rnd.Gaus(0,0.2);
	it->y+=rnd.Gaus(0,0.2);
      }
    mh_hp_smearXY2mm.fill(smeared,PVZ);

}


// ------------ method called once each job just before starting event loop  ------------
void 
MyRPCSimHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyRPCSimHitAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void 
MyRPCSimHitAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    iSetup.get<MuonGeometryRecord>().get(rpcGeo);
}


// ------------ method called when ending the processing of a run  ------------
/*
void 
MyRPCSimHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyRPCSimHitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyRPCSimHitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyRPCSimHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyRPCSimHitAnalyzer);
