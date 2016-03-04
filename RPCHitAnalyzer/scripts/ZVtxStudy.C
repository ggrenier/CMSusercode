//Pixel geometry
// jusqu'a eta=1.5, on touche les 3 barrels : r=4.3, 7.2 et 11.0 cm
// pour eta entre 1.5 et 1.8, on touche les 2 premiers barrels + un disque a z=34.5 cm
// pour eta entre 1.8 et 2.1; on touche les 2 premiers barells + deux disques : z=34.5 cm et z=46.5 cm
// pour eta entre 2.1 et 2.5, on touche le premier barrel et le deuxieme disque (z=46.5cm)

#define ZVtxStudy_cxx
#include "ZVtxStudy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include <string>
#include <map>
#include <iostream>

#include "computeZvtx.C"
#include "Math/Vector3Dfwd.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

class MyStudy
{
public:
  virtual ~MyStudy() {}
  virtual void StartLoop()=0;
  virtual void ProcessEvent(ZVtxStudy* tree)=0;
  virtual MyStudy* recreate()=0; //send a new MyStudy of the correct type 
  static float cspeed;
};


float MyStudy::cspeed=29.9792458; /* cm/ns */


class myAnalysisSet
{
public:
  typedef float PT;
  void add(float thePt,const char* filename) 
  { 
    if (sets.count(thePt)==1) delete sets[thePt];
    sets[thePt]=new Set(filename);
  }
  void setStudy(MyStudy &s)
  {
    for (std::map<PT,Set*>::iterator it=sets.begin(); it != sets.end(); ++it)
      {
	Set& aSet=*(it->second);
	if (aSet.s != NULL) delete aSet.s;
	aSet.s=s.recreate();
	aSet.zv->Loop(aSet.s);
      }
  }
  void printPt() 
  {
    for(std::map<PT,Set*>::iterator it=sets.begin(); it != sets.end();++it)
      std::cout << it->first << " ";
    std::cout << std::endl;
  }
  ZVtxStudy* getTree(float thePt) { if (sets.count(thePt)==1) return sets[thePt]->zv; else return NULL;}
  MyStudy* getStudy(float thePt) { if (sets.count(thePt)==1) return sets[thePt]->s; else return NULL;}
  std::map<float,MyStudy*> getStudies() 
  {
    std::map<float,MyStudy*> c;
    for(std::map<PT,Set*>::iterator it=sets.begin(); it != sets.end();++it)
      c[it->first]=it->second->s;
    return c;
  }
  unsigned int nStudies() {return sets.size();}
  float * thePts() 
  {
    float * pt=new float[sets.size()];
    int i=0;
    for (std::map<PT,Set*>::iterator it=sets.begin(); it != sets.end(); ++it) pt[i++]=it->first;
    return pt;
  }

  MyStudy ** theStudies() 
  {
    MyStudy ** st=new MyStudy*[sets.size()];
    int i=0;
    for (std::map<PT,Set*>::iterator it=sets.begin(); it != sets.end(); ++it) st[i++]=it->second->s;
    return st;
  }

private:
  struct Set
  {
    ZVtxStudy *zv;
    MyStudy *s;
    ~Set() { delete zv; delete s;}
    Set() : zv(0), s(0) {}
    Set(const char* filename) : zv(new ZVtxStudy(filename)), s(0) {}
  };
  std::map<PT,Set*> sets;
};


class MultiMyStudy : public MyStudy
{
public:
  virtual ~MultiMyStudy() { for (std::vector<MyStudy*>::iterator it=studies.begin(); it!= studies.end(); it++) delete *it;}
  void addStudy(MyStudy *s /*ownership transferred*/,std::string name="NOT SET") {studies.push_back(s); studyNames.push_back(name);}
  void StartLoop() {  for(std::vector<MyStudy*>::iterator it=studies.begin(); it!= studies.end(); it++) (*it)->StartLoop();}
  void ProcessEvent(ZVtxStudy* tree) { for(std::vector<MyStudy*>::iterator it=studies.begin(); it!= studies.end(); it++) (*it)->ProcessEvent(tree);}
  MyStudy* recreate() 
  { 
    MultiMyStudy* m=new MultiMyStudy;
    std::vector<std::string>::iterator itname=studyNames.begin();
    for(std::vector<MyStudy*>::iterator it=studies.begin();  it!= studies.end(); ++it, ++itname) m->addStudy( (*it)->recreate(), *itname);
    return m;
  }
  unsigned int nStudies() {return studies.size();}
  MyStudy* getStudy(unsigned int i) {return studies[i];}
  void printStudyNames()
  {
    for(std::vector<std::string>::iterator itname=studyNames.begin(); itname!=studyNames.end(); ++itname)
      std::cout << *itname << std::endl;
  }
  std::string studyName(unsigned int i) {return studyNames[i];}
private:
  std::vector<MyStudy*> studies;
  std::vector<std::string> studyNames;
};



class HitPairStudyEventManager
{
public:
  HitPairStudyEventManager() {}
  virtual void Reset() {}
  virtual ~HitPairStudyEventManager() {}
  virtual HitPairStudyEventManager* create() {return new HitPairStudyEventManager;}
  virtual void startEChitFilling() {;}
  virtual void endEChitFilling() {;}
  virtual bool fillECplus() {return true;}
  virtual bool fillECminus() {return true;}
  virtual bool clearECplus() {return true;}
  virtual bool clearECminus() {return true;}
  virtual bool noProcess() {return false;}
  virtual bool keepPVinfo() {return false;} //to be called after endEChitFilling
};

class HitPairStudyEventManagerMixECplusWithNextECminus : public HitPairStudyEventManager
{
public:
  HitPairStudyEventManagerMixECplusWithNextECminus() : HitPairStudyEventManager(),
						       ECplusHasBeenFilledInThisEvent(false),
						       ECminusHasBeenFilledInThisEvent(false),
						       ECplusHasBeenFilledInPreviousEvent(false),
						       waitForProcessEvent(false),
						       keepPV(false)
  {}
  void Reset() {ECplusHasBeenFilledInThisEvent=ECminusHasBeenFilledInThisEvent=ECplusHasBeenFilledInPreviousEvent=waitForProcessEvent=false;keepPV=false;}
  HitPairStudyEventManager* create() {return new HitPairStudyEventManagerMixECplusWithNextECminus;}

  void endEChitFilling() 
  {
    if (ECplusHasBeenFilledInPreviousEvent)
      {
	if (ECminusHasBeenFilledInThisEvent)
	  {
	    waitForProcessEvent=false; //event completed
	    ECplusHasBeenFilledInThisEvent=ECminusHasBeenFilledInThisEvent=ECplusHasBeenFilledInPreviousEvent=false; //clear state
	    keepPV=true;
	  }
	else //EC minus missing
	  {
	    waitForProcessEvent=true;
	    ECplusHasBeenFilledInThisEvent=ECminusHasBeenFilledInThisEvent=false;
	    keepPV=false;
	  }
      }
    else
      {
	waitForProcessEvent=true;
	if (ECplusHasBeenFilledInThisEvent)
	  {
	    ECplusHasBeenFilledInPreviousEvent=true;
	    ECplusHasBeenFilledInThisEvent=ECminusHasBeenFilledInThisEvent=false;
	    keepPV=true;
	  }
	else 
	  {
	    ECplusHasBeenFilledInThisEvent=ECminusHasBeenFilledInThisEvent=false;
	    keepPV=false;
	  }
      }
  }
  bool fillECplus() 
  {
    ECplusHasBeenFilledInThisEvent=true;
    return !ECplusHasBeenFilledInPreviousEvent;
  }
  bool fillECminus() 
  {
    if (!ECplusHasBeenFilledInPreviousEvent) return false;
    ECminusHasBeenFilledInThisEvent=true;
    return true;
  }

  bool clearECplus() {return !ECplusHasBeenFilledInPreviousEvent;}
  bool noProcess() {return waitForProcessEvent;}
  bool keepPVinfo() {return keepPV;}

private:
  bool ECplusHasBeenFilledInThisEvent;
  bool ECminusHasBeenFilledInThisEvent;
  bool ECplusHasBeenFilledInPreviousEvent;
  bool waitForProcessEvent;
  bool keepPV;
};


class HitPairStudyEventManagerMixECplusNevents : public HitPairStudyEventManager
{
public:
  HitPairStudyEventManagerMixECplusNevents() : HitPairStudyEventManager(),
					       NFill(0),
					       ECplusHasBeenFilledInThisEvent(false),
					       keepPV(false)
  {}
  void Reset() {NFill=0;ECplusHasBeenFilledInThisEvent=false;keepPV=false;}
  HitPairStudyEventManager* create() {return new HitPairStudyEventManagerMixECplusNevents;}

  virtual void startEChitFilling() {if (NFill==2) NFill=0;}
  virtual void endEChitFilling() 
  {
    if (ECplusHasBeenFilledInThisEvent) {NFill++;keepPV=true;} else keepPV=false;
    ECplusHasBeenFilledInThisEvent=false;
  }

  virtual bool fillECplus() {ECplusHasBeenFilledInThisEvent=true; return true;}
  virtual bool fillECminus() {return false;}

  virtual bool clearECplus() { return (NFill==0);}
  virtual bool noProcess() {return (NFill!=2);}
  virtual bool keepPVinfo() {return keepPV;}

private:
  int NFill;
  bool ECplusHasBeenFilledInThisEvent;
  bool keepPV;
};


class HitSelector
{
public:
  HitSelector() {;}
  virtual ~HitSelector() {;}
  virtual bool keepHit(float /*x*/, float /*y*/, float /*z*/, float /*ct*/, unsigned int /*station*/, unsigned int /*ring*/) {return true;}
};

class HitSelectorInverse : public HitSelector
{
public:
  HitSelectorInverse(HitSelector& h) : HitSelector(), sel(&h) {}
  bool keepHit(float x, float y, float z, float ct, unsigned int station, unsigned int ring) 
  {return ! sel->keepHit(x,y,z,ct,station,ring);}
private:
  HitSelector* sel;
};

class HitSelectorRing : public HitSelector
{
public:
  HitSelectorRing(unsigned int ring=1) : HitSelector(), ringnumber(ring) {}
  bool keepHit(float /*x*/, float /*y*/, float /*z*/, float /*ct*/, unsigned int /*station*/, unsigned int ring) {return ring==ringnumber;}
private:
  unsigned int ringnumber;
};


class HitPairStudy : public MyStudy
{
public:
  HitPairStudy(int nRPChitsfordiffmean=0,float TOFdiffmaxbin=1.0,HitPairStudyEventManager* m=NULL,HitSelector* h=NULL);
  void StartLoop() {eventManager->Reset();}
  void ProcessEvent(ZVtxStudy* tree);
  MyStudy* recreate() {return new HitPairStudy(nRPChitsForMean,TOFdiffmaxbinvalue,eventManager->create(),hitSelect);}
  
  class stationXYZTVector : public ROOT::Math::XYZTVector
  {
  public:
    stationXYZTVector() : ROOT::Math::XYZTVector(), station(0) {}
    stationXYZTVector(double vx, double vy, double vz, double vt, unsigned int st=0,unsigned int ri=0) :   ROOT::Math::XYZTVector(vx,vy,vz,vt), station(st), ring(ri) {}
    unsigned int station;
    unsigned int ring;
  };

  class lesHistos
  {
  public:
    TH1F *nsol;
    TH1F *recoPVZ[3];
    TH1F *diffPVZ[3];
    TH1F *recoPVCT[3];
    TH1F *diffPVCT[3];
    lesHistos(std::string prefix);
    void fill(ZTcalc &a,double PVZ, double PVCT);
    void fill(int i, double Z, double T,double PVZ, double PVCT);
  };
  
  lesHistos sameEC;
  lesHistos diffEC;
  TH1F *nhitPlus;
  TH1F *nhitMinus;
  TH1F *samerecoPVZ4hits;
  TH1F *samediffPVZ4hits;
  TH1F *samerecoPVCT4hits;
  TH1F *samediffPVCT4hits;
  TH1F *samerecoPVZlinearfit;
  TH1F *samediffPVZlinearfit;
  TH1F *diffrecoPVZ8hits;
  TH1F *diffdiffPVZ8hits;
  TH1F *diffrecoPVCT8hits;
  TH1F *diffdiffPVCT8hits;
  TH1F *singlehitdiffPVZwithPVCTknown;
  TH1F *singlehitdiffPVCTwithPVZknownat300microns; //300 microns = bigger than current PVZ resolution 
  TH1F *singlehitdiffPVCTwithPVZknownat10microns;  //10 microns = current silicon vertex pixel resolution

  TH1F *singlehitusingFirstPixelBarrelLayerNsol;
  TH1F *singlehitrecoPVZusingFirstPixelBarrelLayer;
  TH1F *singlehitdiffPVZusingFirstPixelBarrelLayer;

  TH1F *straightLineMuRPCHitRdphi;
  TH1F *straightLineMuRPCHitdR;
  TH1F *straightLineMuRPCHitdT;

  TH1F *correctedTOFdiffsameEC(unsigned int station1,unsigned int station2); 
  TH1F *correctedTOFdiffsameEC_hist[6]; 
  TH1F *correctedTOFdiffsameStation[4];
  TH2F *correctedTOFdiffsameStation3vs4;
  TH1F *correctedTOFdiffsameStationMax3and4;

  TH1F *debugNstoredPVZ;
  TH1F *trueDiffPVZ;
  TH2F *trueDiffPVZvstrueDiffPVCT;
  TH1F *correctedTOFdiffsameStation4_TrueDiffPVZless1mm;
  TH1F *correctedTOFdiffsameStation4_TrueDiffPVZless2mm;
  TH1F *correctedTOFdiffsameStation4_TrueDiffPVZless3mm;
  TH1F *correctedTOFdiffsameStation4_TrueDiffPVZless4mm;
  TH1F *correctedTOFdiffsameStation4_TrueDiffPVZabove4mm;
  TH2F *trueDiffPVZvscorrectedTOFdiffsameStation[4];

protected:
  virtual void fill(std::vector<stationXYZTVector>& EC,float x, float y, float z, float ct,unsigned int station, unsigned int ring);
private:
  std::vector<stationXYZTVector> ECplus;
  std::vector<stationXYZTVector> ECminus;  
  ROOT::Math::XYZTVector trueGenMuplus;
  ROOT::Math::XYZTVector trueGenMuminus;

  //equation of a line : x=x0+s*ax, y=y0+s*ay, z=z0+s*az
  //This is to compute the strigth line position
  //of the muon (4-vector given by dir)
  //which comes form the interaction point
  //that is supposed to be
  // x=0, y=0, z=PVZ
  ROOT::Math::XYZVector CoordinateAtS(double s,  ROOT::Math::XYZTVector& dir) { ROOT::Math::XYZVector origin(0,0,PVZ); return origin+(dir.Vect().Unit()*s);}
  ROOT::Math::XYZVector CoordinateAtR(double r,  ROOT::Math::XYZTVector& dir) { return CoordinateAtS( sqrt((r*r)/dir.Vect().Unit().Perp2()), dir);}
  ROOT::Math::XYZVector CoordinateAtZ(double z,  ROOT::Math::XYZTVector& dir) { return CoordinateAtS( (z-PVZ)/dir.Vect().Unit().Z(), dir);}

  ROOT::Math::XYZTVector extrapolate4D(ROOT::Math::XYZVector& position,ROOT::Math::XYZTVector& RPChit){return ROOT::Math::XYZTVector(position.X(),position.Y(),position.Z(),RPChit.T()-(RPChit.Vect()-position).R());}

  void computeFirstPixelLayer(std::vector<stationXYZTVector>& EC,ROOT::Math::XYZTVector& trueDir);
  void computeFirstPixelLayerPlus() { computeFirstPixelLayer(ECplus,trueGenMuplus);}
  void computeFirstPixelLayerMinus() { computeFirstPixelLayer(ECminus,trueGenMuminus);}

  void computeStraightLineVsRealRPCHit(std::vector<stationXYZTVector>& EC,ROOT::Math::XYZTVector& trueDir);


  void computeSame(std::vector<stationXYZTVector>& EC);
  void computeDiff();
  void addResults(ZTcalc &a, int &ncalc, double &sumZ, double &sumT);
  double computeT0(ROOT::Math::XYZTVector& a, double z0);

protected:

  HitPairStudyEventManager *eventManager;
  HitSelector* hitSelect;
 
  float TOFdiffmaxbinvalue;
  int nRPChitsForMean;

  double PVZ,PVCT;
  std::vector<double> usedPVZ,usedPVCT;
  TRandom3 rnd;
  static HitPairStudyEventManager defaultEventManager;
  static HitSelector defaultHitSelector;
};

HitPairStudyEventManager HitPairStudy::defaultEventManager=HitPairStudyEventManager();
HitSelector HitPairStudy::defaultHitSelector=HitSelector();


class HitPairStudySmearTime : public HitPairStudy
{
public:
  HitPairStudySmearTime(float ts,int nRPChitsfordiffmean=0,float TOFdiffmaxbin=1.0,HitPairStudyEventManager* m=NULL,HitSelector* h=NULL) : HitPairStudy(nRPChitsfordiffmean,TOFdiffmaxbin, m, h), timeSmear(ts) {}
  MyStudy* recreate() {return new HitPairStudySmearTime(timeSmear,nRPChitsForMean,TOFdiffmaxbinvalue,eventManager->create(),hitSelect);}
private:
  void fill(std::vector<stationXYZTVector>& EC,float x, float y, float z, float ct, unsigned int station, unsigned int ring)
  {
    ct+=MyStudy::cspeed*rnd.Gaus(0,timeSmear);
    HitPairStudy::fill(EC,x,y,z,ct,station,ring);
  }
  float timeSmear; //ns
};


class HitPairStudySmearTimeAndStrip : public HitPairStudy
{
public:
  HitPairStudySmearTimeAndStrip(float ts,int nRPChitsfordiffmean=0,float TOFdiffmaxbin=1.0,HitPairStudyEventManager* m=NULL,HitSelector* h=NULL) : HitPairStudy(nRPChitsfordiffmean,TOFdiffmaxbin,m,h), timeSmear(ts) {}
  MyStudy* recreate() {return new HitPairStudySmearTimeAndStrip(timeSmear,nRPChitsForMean,TOFdiffmaxbinvalue,eventManager->create(),hitSelect);}
protected:
  void fill(std::vector<stationXYZTVector>& EC,float x, float y, float z, float ct, unsigned int station, unsigned int ring)
  {
    ct+=MyStudy::cspeed*rnd.Gaus(0,timeSmear);
    float r=sqrt(x*x+y*y);
    float rnew=r+2*MyStudy::cspeed*rnd.Gaus(0,timeSmear);
    HitPairStudy::fill(EC,x*rnew/r,y*rnew/r,z,ct,station,ring);
  }
protected:
  float timeSmear; //ns
};


class HitPairStudyWithMorePairs : public HitPairStudySmearTimeAndStrip
{
public:
  HitPairStudyWithMorePairs(int nextraLayers, float ts, float extralayerZdelta=0.1 /*in cm so this is 1 mm*/, int nRPChitsfordiffmean=0,float TOFdiffmaxbin=1.0,HitPairStudyEventManager* m=NULL,HitSelector* h=NULL) 
    : HitPairStudySmearTimeAndStrip(ts,nRPChitsfordiffmean,TOFdiffmaxbin,m,h), nextra(nextraLayers), ZDelta(extralayerZdelta) {}
  MyStudy* recreate() {return new HitPairStudyWithMorePairs(nextra,timeSmear,ZDelta,nRPChitsForMean,TOFdiffmaxbinvalue,eventManager->create(),hitSelect);}
private:
  void fill(std::vector<stationXYZTVector>& EC,float x, float y, float z, float ct, unsigned int station,unsigned int ring)
  {
    ROOT::Math::XYZTVector PV(0,0,PVZ,PVCT);
    ROOT::Math::XYZTVector hit(x,y,z,ct);
    ROOT::Math::XYZTVector centered=hit-PV;
    for (int i=0; i<=nextra; i++)
      {
	ROOT::Math::XYZTVector newhit=centered*((centered.Z()+i*ZDelta)/centered.Z())+PV;
	HitPairStudySmearTimeAndStrip::fill(EC,newhit.X(),newhit.Y(),newhit.Z(),newhit.T(),station,ring);
      }
  }
  int nextra;
  float ZDelta;
};



TH1F* HitPairStudy::correctedTOFdiffsameEC(unsigned int station1,unsigned int station2)
{
  if (station1==station2) return NULL;
  if (station1>4 || station1==0 || station2>4 || station2==0) return NULL;
  if (station1 > station2)
    {
      int tmp=station1;
      station1=station2;
      station2=tmp;
    }
  if (station1==1) return correctedTOFdiffsameEC_hist[station2-2]; //[0]=1-2 [1]=1-3 [2]=1-4
  if (station1==2) return correctedTOFdiffsameEC_hist[station2]; //[3]=2-3 [4]=2-4
  return correctedTOFdiffsameEC_hist[5]; // [5]=3-4
}

HitPairStudy::HitPairStudy(int nRPChitsfordiffmean,float TOFdiffmaxbin,HitPairStudyEventManager* m,HitSelector* h): sameEC("same"), diffEC("diff"), eventManager(m), hitSelect(h), TOFdiffmaxbinvalue(TOFdiffmaxbin), nRPChitsForMean(nRPChitsfordiffmean)/*, rnd(rand())*/
{
  if (NULL==eventManager) eventManager=&defaultEventManager;
  if (NULL==hitSelect) hitSelect=&defaultHitSelector;
  Int_t nbin=200;
  nhitPlus=new TH1F("nhitECplus","Number of muon hits in EC plus",25,0,25);
  nhitMinus=new TH1F("nhitECminus","Number of muon hits in EC minus",25,0,25);
  samerecoPVZ4hits=new TH1F("samerecoPVZ4hits","sameEC : mean reconstructed PVZ with 4 hits (cm)",100,-100,100);
  samediffPVZ4hits=new TH1F("samediffPVZ4hits","sameEC : mean reconstructed PVZ with 4 hits - true PVZ (cm)",nbin,-20,20);
  samerecoPVCT4hits=new TH1F("samerecoPVCT4hits","sameEC : mean reconstructed PV CT with 4 hits (cm)",100,-100,100);
  samediffPVCT4hits=new TH1F("samediffPVCT4hits","sameEC : mean reconstructed PV CT with 4 hits - true PV CT (cm)",nbin,-20,20);
  diffrecoPVZ8hits=new TH1F("diffrecoPVZ8hits","diffEC : mean reconstructed PVZ with 8 hits (cm)",100,-100,100);
  diffdiffPVZ8hits=new TH1F("diffdiffPVZ8hits","diffEC : mean reconstructed PVZ with 8 hits - true PVZ (cm)",nbin,-20,20);
  diffrecoPVCT8hits=new TH1F("diffrecoPVCT8hits","diffEC : mean reconstructed PV CT with 8 hits (cm)",100,-100,100);
  diffdiffPVCT8hits=new TH1F("diffdiffPVCT8hits","diffEC : mean reconstructed PV CT with 8 hits - true PV CT (cm)",nbin,-20,20);
  singlehitdiffPVZwithPVCTknown=new TH1F("singlediffPVZ","single RPC hits : reco PVZ - true PVZ when using true PV CT (cm)",nbin,-20,20);
  singlehitdiffPVCTwithPVZknownat300microns=new TH1F("singlediffPVCT300mPVZ","single RPC hits : reco PV T - true P CT with PVZ smear 300 microns (ns)",nbin,-1,1);
  singlehitdiffPVCTwithPVZknownat10microns=new TH1F("singlediffPVCT10mPVZ","single RPC hits : reco PV T - true P CT with PVZ smear 10 microns (ns)",nbin,-1,1);
  samerecoPVZlinearfit=new TH1F("samerecoPVZlinear","sameEC : reconstructed PVZ from linear r=az+b (cm)",100,-100,100);
  samediffPVZlinearfit=new TH1F("samediffPVZlinear","sameEC : reconstructed PVZ from linear r=az+b - true PVZ (cm)",nbin,-20,20);
  singlehitusingFirstPixelBarrelLayerNsol=new TH1F("singleHitFirstPixelBarrelNsol","single hit : Number of solution when using first pixel barrel",3,0,3);
  singlehitrecoPVZusingFirstPixelBarrelLayer=new TH1F("singleHitrecoFirstPixelBarrel","single hit : reconstructed PVZ using first pixel barrel (cm)",100,-100,100);
  singlehitdiffPVZusingFirstPixelBarrelLayer=new TH1F("singleHitdiffFirstPixelBarrel","single hit : reconstructed PVZ using first pixel barrel - true PVZ (cm)",nbin,-20,20);
  straightLineMuRPCHitRdphi =  new TH1F("straightLineMuRPCHitRdphi","straigth line true muon extrapolated position - RPC hit position : r*Delta phi (cm)",1000,-500,500);
  straightLineMuRPCHitdR =  new TH1F("straightLineMuRPCHitdR","straigth line true muon extrapolated position - RPC hit position : Delta r (cm)",1000,-50,50);
  straightLineMuRPCHitdT =  new TH1F("straightLineMuRPCHitdT","straigth line true muon extrapolated position - RPC hit position : Delta t (ns)",1000,-1,1);

  std::string stationPair[6]={"1_2","1_3","1_4","2_3","2_4","3_4"};
  for (int i=0; i<6; ++i) correctedTOFdiffsameEC_hist[i]=new TH1F((std::string("corrTOFdiffStations_")+stationPair[i]).c_str(),
								  (std::string("Diff of (TOF-r/C) in ns for hits in station pair : ")+stationPair[i]).c_str(),
								  100,0,TOFdiffmaxbinvalue);
  char chiffre='1';
  for (int i=0; i<4; ++i,++chiffre) 
    {
      correctedTOFdiffsameStation[i]=new TH1F((std::string("corrTOFdiffSameStation_")+chiffre).c_str(),
					      (std::string("Diff of (TOF-r/C) in ns for muon hits in same station : ")+chiffre).c_str(),
					      100,0,TOFdiffmaxbinvalue);
      trueDiffPVZvscorrectedTOFdiffsameStation[i]=new TH2F((std::string("trueDiffPVZvscorrectedTOFdiffsameStation_")+chiffre).c_str(),
							   (std::string("true diff PVZ vs Diff of (TOF-r/C) in ns for muon hits in same station : ")+chiffre).c_str(),
							   100,0,TOFdiffmaxbinvalue,nbin,-20,20);
    }

  correctedTOFdiffsameStation3vs4 = new TH2F("correctedTOFdiffsameStation3vs4","Diff of (TOF-r/C) in ns for muon hits in same station : 3 vs 4",
					     100,0,TOFdiffmaxbinvalue,100,0,TOFdiffmaxbinvalue);
  correctedTOFdiffsameStationMax3and4 = new TH1F("correctedTOFdiffsameStationMax3and4",
						 "Diff of (TOF-r/C) in ns for muon hits in same station : max value for 3 and 4",
						 100,0,TOFdiffmaxbinvalue);

  debugNstoredPVZ = new TH1F("debugNstoredPVZ","DEBUG histo : N stored PVZ values",10,0,10);
  trueDiffPVZ = new TH1F("trueDiffPVZ","true PVZ diff in cm ",400,-20,20);
  trueDiffPVZvstrueDiffPVCT = new TH2F("trueDiffPVZvstruediffPVCT","true PVZ diff in cm vs true PV CT diff in cm",400,-20,20,400,-20,20);
  correctedTOFdiffsameStation4_TrueDiffPVZless1mm = new TH1F("correctedTOFdiffsameStation4_TrueDiffPVZless1mm",
							     "Diff of (TOF-r/C) in ns for muon hits in same station 4 and true PVZ diff less than 1 mm",
							     100,0,TOFdiffmaxbinvalue);
  correctedTOFdiffsameStation4_TrueDiffPVZless2mm = new TH1F("correctedTOFdiffsameStation4_TrueDiffPVZless2mm",
							     "Diff of (TOF-r/C) in ns for muon hits in same station 4 and true PVZ diff less than 2 mm",
							     100,0,TOFdiffmaxbinvalue);
  correctedTOFdiffsameStation4_TrueDiffPVZless3mm = new TH1F("correctedTOFdiffsameStation4_TrueDiffPVZless3mm",
							     "Diff of (TOF-r/C) in ns for muon hits in same station 4 and true PVZ diff less than 3 mm",
							     100,0,TOFdiffmaxbinvalue);
  correctedTOFdiffsameStation4_TrueDiffPVZless4mm = new TH1F("correctedTOFdiffsameStation4_TrueDiffPVZless4mm",
							     "Diff of (TOF-r/C) in ns for muon hits in same station 4 and true PVZ diff less than 4 mm",
							     100,0,TOFdiffmaxbinvalue);
  correctedTOFdiffsameStation4_TrueDiffPVZabove4mm = new TH1F("correctedTOFdiffsameStation4_TrueDiffPVZabove4mm",
							     "Diff of (TOF-r/C) in ns for muon hits in same station 4 and true PVZ diff greater than 4 mm",
							     100,0,TOFdiffmaxbinvalue);
  
}

HitPairStudy::lesHistos::lesHistos(std::string prefix)
{
  Int_t nbin=200;
  nsol=new TH1F((prefix+"Nsol").c_str(),(prefix+"EC : Number of solutions").c_str(),3,0,3);
  recoPVZ[0]=new TH1F((prefix+"recoPVZall").c_str(),(prefix+"EC : reconstructed PVZ (cm) All").c_str(),100,-100,100);
  recoPVZ[1]=new TH1F((prefix+"recoPVZcloseZ").c_str(),(prefix+"EC : reconstructed PVZ (cm) Z closest to 0").c_str(),100,-100,100);
  recoPVZ[2]=new TH1F((prefix+"recoPVZcloseT").c_str(),(prefix+"EC : reconstructed PVZ (cm) CT closest to 0").c_str(),100,-100,100);
  recoPVCT[0]=new TH1F((prefix+"recoPVCTall").c_str(),(prefix+"EC : reconstructed PV CT (cm) All").c_str(),100,-100,100);
  recoPVCT[1]=new TH1F((prefix+"recoPVCTcloseZ").c_str(),(prefix+"EC : reconstructed PV CT (cm) Z closest to 0").c_str(),100,-100,100);  recoPVCT[2]=new TH1F((prefix+"recoPVCTcloseT").c_str(),(prefix+"EC : reconstructed PV CT (cm) CT closest to 0").c_str(),100,-100,100);
  diffPVZ[0]=new TH1F((prefix+"diffPVZall").c_str(),(prefix+"EC : reconstructed PVZ - true PVZ (cm) All").c_str(),nbin,-20,20);
  diffPVZ[1]=new TH1F((prefix+"diffPVZcloseZ").c_str(),(prefix+"EC : reconstructed PVZ - true PVZ (cm) Z closest to 0").c_str(),nbin,-20,20);
  diffPVZ[2]=new TH1F((prefix+"diffPVZcloseT").c_str(),(prefix+"EC : reconstructed PVZ - true PVZ (cm) CT closest to 0").c_str(),nbin,-20,20);
  diffPVCT[0]=new TH1F((prefix+"diffPVCTall").c_str(),(prefix+"EC : reconstructed PVCT - true PVCT (cm) All").c_str(),nbin,-20,20);
  diffPVCT[1]=new TH1F((prefix+"diffPVCTcloseZ").c_str(),(prefix+"EC : reconstructed PVCT - true PVCT (cm) Z closest to 0").c_str(),nbin,-20,20);
  diffPVCT[2]=new TH1F((prefix+"diffPVCTcloseT").c_str(),(prefix+"EC : reconstructed PVCT - true PVCT (cm) T closest to 0").c_str(),nbin,-20,20);
}

void HitPairStudy::lesHistos::fill(ZTcalc &a,double PVZ, double PVCT)
{
  int nsolutions=a.Nsolution();
  nsol->Fill(nsolutions);
  if (nsolutions>=1) fill(0,a.Zsol(0),a.Tsol(0),PVZ,PVCT);
  if (nsolutions==2) 
    {
      fill(0,a.Zsol(1),a.Tsol(1),PVZ,PVCT);
      if (fabs(a.Zsol(0))<fabs(a.Zsol(1)))
	fill(1,a.Zsol(0),a.Tsol(0),PVZ,PVCT);
      else 
	fill(1,a.Zsol(1),a.Tsol(1),PVZ,PVCT);
      if (fabs(a.Tsol(0))<fabs(a.Tsol(1)))
	fill(2,a.Zsol(0),a.Tsol(0),PVZ,PVCT);
      else 
	fill(2,a.Zsol(1),a.Tsol(1),PVZ,PVCT);
    }
  if (nsolutions==1)
    {
      fill(1,a.Zsol(0),a.Tsol(0),PVZ,PVCT);
      fill(2,a.Zsol(0),a.Tsol(0),PVZ,PVCT);
    }
}

void HitPairStudy::lesHistos::fill(int i, double Z, double T,double PVZ, double PVCT)
{
  recoPVZ[i]->Fill(Z);
  diffPVZ[i]->Fill(Z-PVZ);
  recoPVCT[i]->Fill(T);
  diffPVCT[i]->Fill(T-PVCT);
}



void HitPairStudy::addResults(ZTcalc &a, int &ncalc, double &sumZ, double & sumT)
{
  int nsolutions=a.Nsolution();
  if (nsolutions==0) return;
  if (nsolutions==1) 
    {
      ncalc++;
      sumZ+=a.Zsol(0);
      sumT+=a.Tsol(0);
    }
  if (nsolutions==2)
    {
      int i=0;
      if (fabs(a.Tsol(1))<fabs(a.Tsol(0))) i=1;
      ncalc++;
      sumZ+=a.Zsol(i);
      sumT+=a.Tsol(i);
    }
} 



void HitPairStudy::ProcessEvent(ZVtxStudy* tree)
{
  PVZ=tree->PVZ;
  PVCT=tree->PVCT;
  std::vector<float> &genPx=*(tree->genMuonPx);
  std::vector<float> &genPy=*(tree->genMuonPy);
  std::vector<float> &genPz=*(tree->genMuonPz);
  for (unsigned i=0; i<genPz.size();++i)
    {
      ROOT::Math::XYZTVector* toset=&trueGenMuplus;
      if (genPz[i]<0) toset=&trueGenMuminus;
      double mumass=0.10566; //GeV
      toset->SetE(sqrt(genPx[i]*genPx[i]+genPy[i]*genPy[i]+genPz[i]*genPz[i]+mumass*mumass));
      toset->SetPx(genPx[i]);
      toset->SetPy(genPy[i]);
      toset->SetPz(genPz[i]);
    }     

  eventManager->startEChitFilling();
  if (eventManager->clearECplus())  ECplus.clear(); 
  if (eventManager->clearECminus()) ECminus.clear();

  std::vector<int> &pid=*(tree->pid);
  std::vector<int> &region=*(tree->region);
  std::vector<int> &station=*(tree->station);
  std::vector<int> &ring=*(tree->ring);
  for (unsigned int i=0; i<pid.size();i++)
    {
      if (std::abs(pid[i])==13)
	{
	  std::vector<float> &x=*(tree->x);
	  std::vector<float> &y=*(tree->y);
	  std::vector<float> &z=*(tree->z);
	  std::vector<float> &t=*(tree->t);
	  if (hitSelect->keepHit(x[i],y[i],z[i],MyStudy::cspeed*t[i],station[i],ring[i]))
	    {
	      if (region[i]==1 && eventManager->fillECplus()) fill(ECplus,x[i],y[i],z[i],MyStudy::cspeed*t[i],station[i],ring[i]);
	      if (region[i]==-1 && eventManager->fillECminus()) fill(ECminus,x[i],y[i],z[i],MyStudy::cspeed*t[i],station[i],ring[i]);
	    }
	}
    }
  eventManager->endEChitFilling();
  if (eventManager->keepPVinfo())
    {
      usedPVZ.push_back(PVZ);
      usedPVCT.push_back(PVCT);
    }
  if (eventManager->noProcess()) return;
  nhitPlus->Fill(ECplus.size());
  nhitMinus->Fill(ECminus.size());
  computeSame(ECplus);
  computeSame(ECminus);
  computeDiff();
  computeFirstPixelLayerPlus();
  computeFirstPixelLayerMinus();
  computeStraightLineVsRealRPCHit(ECplus,trueGenMuplus);
  computeStraightLineVsRealRPCHit(ECminus,trueGenMuminus);
  debugNstoredPVZ->Fill(usedPVZ.size());
  usedPVZ.clear();
  usedPVCT.clear();
}

void HitPairStudy::fill(std::vector<stationXYZTVector>& EC,float x, float y, float z, float ct,unsigned int station, unsigned int ring)
{
  EC.push_back(stationXYZTVector(x,y,z,ct,station,ring));
}


void HitPairStudy::computeFirstPixelLayer(std::vector<stationXYZTVector>& EC,ROOT::Math::XYZTVector& trueDir)
{
  ROOT::Math::XYZVector pos=CoordinateAtR( 4.3,trueDir  );
  //std::cout << "DEBUG " << pos.Rho() << " ##### " << trueDir.Px() << " " << trueDir.Py() << " " << trueDir.Pz() << std::endl;
  
  for (std::vector<stationXYZTVector>::iterator it=EC.begin(); it != EC.end(); ++it)
    {
      ROOT::Math::XYZTVector pos4D=extrapolate4D( pos , *it );
      //std::cout << "DEBUG 4D " <<  pos4D.Rho() << " " << pos4D.Z() << " " << pos4D.T() << " {{{}}} " << it->Rho() << " " << it->Z() << " " << it->T() << std::endl;
      ZTcalc a( pos4D , *it);
      a.compute();
      singlehitusingFirstPixelBarrelLayerNsol->Fill(a.Nsolution());
      int index=-1;
      if (a.Nsolution()==1) index=0;
      if (a.Nsolution()==2)
	{
	  index=0;
	  if (fabs(a.Tsol(1))<fabs(a.Tsol(0))) index=1;
	}
      if (index != -1) 
	{
	  //std::cout << "DEBUG " << a.Zsol(index) << std::endl;
	  singlehitrecoPVZusingFirstPixelBarrelLayer->Fill(a.Zsol(index));
	  singlehitdiffPVZusingFirstPixelBarrelLayer->Fill(a.Zsol(index)-PVZ);
	}
    }
}


void HitPairStudy::computeStraightLineVsRealRPCHit(std::vector<stationXYZTVector>& EC,ROOT::Math::XYZTVector& trueDir)
{
  for (std::vector<stationXYZTVector>::iterator it=EC.begin(); it != EC.end(); ++it)
    {
      ROOT::Math::XYZVector pos=CoordinateAtZ(it->Z(),trueDir);
      double dPhi = pos.Phi() - it->Phi();
      if (dPhi>ROOT::Math::Pi()) dPhi-=2*ROOT::Math::Pi();
      if (dPhi<-ROOT::Math::Pi()) dPhi+=2*ROOT::Math::Pi();
      double dRho = pos.Rho() - it->Rho();
      double dT = (PVCT + sqrt(pos.Perp2() + (pos.Z()-PVZ)*(pos.Z()-PVZ))) - it->T();
      straightLineMuRPCHitRdphi->Fill(it->Rho()*dPhi);
      straightLineMuRPCHitdR->Fill(dRho);
      straightLineMuRPCHitdT->Fill(dT);
    }
}


double HitPairStudy::computeT0(ROOT::Math::XYZTVector& a, double z0)
{
  double extra=sqrt(a.Perp2()+(a.Z()-z0)*(a.Z()-z0));
  double tplus=a.T()+extra;
  double tminus=a.T()-extra;
  double t=tplus;
  if (fabs(tminus)<fabs(tplus)) t=tminus;
  return t/MyStudy::cspeed;
}


void HitPairStudy::computeSame(std::vector<stationXYZTVector>& EC)
{
  int ncalcul=0;
  double calcZ=0;
  double calcT=0;
  float trueDPVZ=0;
  float trueDPVCT=0;
  if (usedPVZ.size() >= 2) {trueDPVZ=usedPVZ[0]-usedPVZ[1];trueDPVCT=usedPVCT[0]-usedPVCT[1];}
  trueDiffPVZ->Fill(trueDPVZ);
  trueDiffPVZvstrueDiffPVCT->Fill(trueDPVCT,trueDPVZ);
  std::vector<double> memorizeDiffEC3,memorizeDiffEC4;
  for (unsigned int i=0; i<EC.size(); i++)
    {
      double itime=(EC[i].T()-EC[i].R())/MyStudy::cspeed;
      for (unsigned int j=i+1; j<EC.size(); j++)
	{
	  double jtime=(EC[j].T()-EC[j].R())/MyStudy::cspeed;
	  ZTcalc a(EC[i],EC[j]);
	  a.compute();
	  sameEC.fill(a,PVZ,PVCT);
	  addResults(a, ncalcul, calcZ, calcT);
	  double z0=EC[i].Z()-(EC[i].Z()-EC[j].Z())/(1-sqrt(EC[j].Perp2()/EC[i].Perp2()));
	  samerecoPVZlinearfit->Fill(z0);
	  samediffPVZlinearfit->Fill(z0-PVZ);

	  TH1F* histdiffTOF=correctedTOFdiffsameEC(EC[i].station,EC[j].station);
	  if (NULL != histdiffTOF) histdiffTOF->Fill(fabs(itime-jtime));
	  if (EC[i].station==EC[j].station) 
	    {
	      correctedTOFdiffsameStation[EC[i].station-1]->Fill(fabs(itime-jtime));
	      trueDiffPVZvscorrectedTOFdiffsameStation[EC[i].station-1]->Fill(fabs(itime-jtime),trueDPVZ);
	      if (EC[i].station==3) memorizeDiffEC3.push_back(fabs(itime-jtime));
	      if (EC[i].station==4)
		{
		  memorizeDiffEC4.push_back(fabs(itime-jtime));
		  if (fabs(trueDPVZ)<0.1) correctedTOFdiffsameStation4_TrueDiffPVZless1mm->Fill(fabs(itime-jtime));
		  if (fabs(trueDPVZ)<0.2) correctedTOFdiffsameStation4_TrueDiffPVZless2mm->Fill(fabs(itime-jtime));
		  if (fabs(trueDPVZ)<0.3) correctedTOFdiffsameStation4_TrueDiffPVZless3mm->Fill(fabs(itime-jtime));
		  if (fabs(trueDPVZ)<0.4) correctedTOFdiffsameStation4_TrueDiffPVZless4mm->Fill(fabs(itime-jtime));
		  else correctedTOFdiffsameStation4_TrueDiffPVZabove4mm->Fill(fabs(itime-jtime));
		}
	    }
	}
      //compute single hit resolution
      double extra=sqrt((EC[i].T()-PVCT)*(EC[i].T()-PVCT)-EC[i].Perp2());
      double zplus=EC[i].Z()+extra;
      double zminus=EC[i].Z()-extra;
      double z=zplus;
      if (fabs(zminus)<fabs(zplus)) z=zminus;
      singlehitdiffPVZwithPVCTknown->Fill(z-PVZ);

      singlehitdiffPVCTwithPVZknownat300microns->Fill(computeT0(EC[i],PVZ+rnd.Gaus(0,0.03))-PVCT/MyStudy::cspeed);
      singlehitdiffPVCTwithPVZknownat10microns->Fill(computeT0(EC[i],PVZ+rnd.Gaus(0,0.001))-PVCT/MyStudy::cspeed);
    }
  if (memorizeDiffEC3.size()==1 &&  memorizeDiffEC4.size()==1) 
    {
      correctedTOFdiffsameStation3vs4->Fill(memorizeDiffEC4[0],memorizeDiffEC3[0]);
      if (memorizeDiffEC4[0]>memorizeDiffEC3[0])
	correctedTOFdiffsameStationMax3and4->Fill(memorizeDiffEC4[0]);
      else
	correctedTOFdiffsameStationMax3and4->Fill(memorizeDiffEC3[0]);
    }
  if (ncalcul==0) return;
  calcZ=calcZ/ncalcul;
  calcT=calcT/ncalcul;
  samerecoPVZ4hits->Fill(calcZ);
  samediffPVZ4hits->Fill(calcZ-PVZ);
  samerecoPVCT4hits->Fill(calcT);
  samediffPVCT4hits->Fill(calcT-PVCT);
  
}

void HitPairStudy::computeDiff()
{
  int ncalcul=0;
  double calcZ=0;
  double calcT=0;
  for (std::vector<stationXYZTVector>::iterator itplus=ECplus.begin();
       itplus != ECplus.end(); ++itplus)
    for (std::vector<stationXYZTVector>::iterator itminus=ECminus.begin();
       itminus != ECminus.end(); ++itminus)
      {
	ZTcalc a(*itplus,*itminus);
	a.compute();
	diffEC.fill(a,PVZ,PVCT);
	addResults(a, ncalcul, calcZ, calcT);
      }
  if (nRPChitsForMean>=2 && ECplus.size()+ECminus.size()!=std::abs(nRPChitsForMean)) return;
  if (ncalcul==0) return;
  calcZ=calcZ/ncalcul;
  calcT=calcT/ncalcul;
  diffrecoPVZ8hits->Fill(calcZ);
  diffdiffPVZ8hits->Fill(calcZ-PVZ);
  diffrecoPVCT8hits->Fill(calcT);
  diffdiffPVCT8hits->Fill(calcT-PVCT);
}




void ZVtxStudy::Loop(MyStudy* study)
{
//   In a ROOT session, you can do:
//      Root > .L ZVtxStudy.C
//      Root > ZVtxStudy t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   study->StartLoop();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      study->ProcessEvent(this);
   }
}
