#ifndef RPC_MultiHisto_HH
#define RPC_MultiHisto_HH


//NB : this file is both compiled by scram and interpreted by ROOT cint

#include <string>
#include <iostream>

#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "THStack.h"
#include "TLegend.h"

//helper class to retrieve/or book histos
class HistoBooker
{
 public:
  //small helper class : axis
  class Axis
  {
  public:
    Axis(int,float,float);
    int Nbin;
    float minValue,maxValue;
  };

  virtual ~HistoBooker() {}
  TH1* book(std::string name, std::string title, Axis X, Axis *Y=NULL);
  virtual bool YouOwn() {return true;}
 private:
  virtual TH1F* book1D(std::string name, std::string title, Axis X);
  virtual TH2F* book2D(std::string name, std::string title, Axis X,Axis Y);
};
 
class HistoRetriever : public HistoBooker
{
 public:
  HistoRetriever(TDirectory &dir);
  virtual bool YouOwn() {return false;}
 private:
  TDirectory *_dir; //Interactive ROOT doesn't know how to handle References
  virtual TH1F* book1D(std::string name, std::string title, Axis X);
  virtual TH2F* book2D(std::string name, std::string title, Axis X,Axis Y);

  void printWarning(std::string name,TH1* h);
};

class RPC_MultiHistoBase
{
 public:
  enum {NECs=2,Nstations=4,Nrings=3};
  
  RPC_MultiHistoBase(HistoBooker *hb=NULL);
  ~RPC_MultiHistoBase();
  void Write();  

  void setHistoBooker(HistoBooker *booker);
  void book(std::string name, std::string title, HistoBooker::Axis X, HistoBooker::Axis *Y=NULL);

  void setUse(int ecNumber,unsigned stationNumber,unsigned int ringNumber,bool use=true);
  void setUseAll(bool use=true);
  void setUseRing(unsigned int ringNumber,bool use=true);

  THStack* getStack();//ownership transfered to caller
  TH1* getHistoSum(); //ownership transfered to caller 
  TLegend* getLegend(double x1, double y1, double x2, double y2, const char* header="", const char* option="brNDC");

  void setLineAttributes();
  void setFillAttributes();
  void setMarkerAttributes();
  void setAllAttributes();
  void resetLineAttributes();
  void resetFillAttributes();
  void resetMarkerAttributes();

 protected:
  TH1* getRingHistoBase(int ecNumber, unsigned int stationNumber, unsigned int ringNumber);
 
  std::string ringName(unsigned int);
  std::string ringNameShort(unsigned int);
  std::string stationName(unsigned int);
  std::string ECName(int);
  std::string histoNameAppend(int ecNumber,unsigned int stationNumber, unsigned int ringNumber);
  std::string histoTitleAppend(int ecNumber,unsigned int stationNumber, unsigned int ringNumber);
  std::string histoTitleWithoutAppend(std::string appendedTitle,int ecNumber,unsigned int stationNumber, unsigned int ringNumber);
#ifndef __CINT__
  std::string usehitToShortString(); //for stack creation and histo sum
#endif

  unsigned int ecNumberToIndex(int ecNumber);
  int indexToEcNumber(unsigned int index);

  TH1* _hist[NECs][Nstations][Nrings];
  bool _useHist[NECs][Nstations][Nrings];

  std::string _histName;

  HistoBooker *_currentBooker;
  HistoBooker _defaultBooker;

  bool _IOwnHistos;
};


class RPC_MultiHisto : public RPC_MultiHistoBase
{
 public:
  RPC_MultiHisto(HistoBooker *hb=NULL);
  void book(std::string name, std::string title, int nbin, float xmin, float xmax);
  void Fill(int ecNumber,unsigned stationNumber,unsigned int ringNumber,double t);
  TH1F& getRingHisto(int ecNumber, unsigned int stationNumber, unsigned int ringNumber);
};

class RPC_MultiHisto2D : public RPC_MultiHistoBase
{
 public:
  RPC_MultiHisto2D(HistoBooker *hb=NULL);
  void book(std::string name, std::string title, int nbinx, float xmin, float xmax,int nbiny, float ymin, float ymax);
  void Fill(int ecNumber,unsigned stationNumber,unsigned int ringNumber,double x,double y);
  TH2F& getRingHisto(int ecNumber, unsigned int stationNumber, unsigned int ringNumber);
};



inline HistoBooker::Axis::Axis(int n,float min, float max) : Nbin(n),minValue(min),maxValue(max) {}

inline TH1* HistoBooker::book(std::string name, std::string title, Axis X, Axis *Y)
{
  if (Y==NULL) return book1D(name,title,X);
  else         return book2D(name,title,X,*Y);
}

inline TH1F* HistoBooker::book1D(std::string name, std::string title, Axis X)
{
  return new TH1F(name.c_str(), title.c_str(), X.Nbin, X.minValue, X.maxValue);
}

inline TH2F* HistoBooker::book2D(std::string name, std::string title, Axis X,Axis Y)
{
  return new TH2F(name.c_str(), title.c_str(),
                  X.Nbin, X.minValue, X.maxValue ,
                  Y.Nbin, Y.minValue, Y.maxValue);
}


inline HistoRetriever::HistoRetriever(TDirectory &dir) : HistoBooker(), _dir(&dir) {}

inline TH1F* HistoRetriever::book1D(std::string name, std::string title, Axis X)
{
  TH1F *h=NULL; _dir->GetObject(name.c_str(),h);printWarning(name,h);return h;
}

inline TH2F* HistoRetriever::book2D(std::string name, std::string title, Axis X,Axis Y)
{
  TH2F *h=NULL;_dir->GetObject(name.c_str(),h);printWarning(name,h);return h;
}

inline void HistoRetriever::printWarning(std::string name,TH1* h)
{
  if (NULL==h) std::cout << "WARNING : did not retrieved " << name << std::endl;
}


inline void RPC_MultiHistoBase::setHistoBooker(HistoBooker *booker) {_currentBooker=booker; _IOwnHistos=booker->YouOwn();}


inline std::string RPC_MultiHistoBase::ringName(unsigned int i) { std::string A="Ring0"; A[4]+=i; return A;}
inline std::string RPC_MultiHistoBase::ringNameShort(unsigned int i) {std::string A="/0"; A[1]+=i; return A;}
inline std::string RPC_MultiHistoBase::stationName(unsigned int i) {std::string A="RE0"; A[2]+=i; return A;}
inline std::string RPC_MultiHistoBase::ECName(int i) {if (i==-1) return "ECminus"; else return "ECplus";}
inline std::string RPC_MultiHistoBase::histoNameAppend(int ecNumber,unsigned int stationNumber, unsigned int ringNumber)
{ return std::string("_")+stationName(stationNumber)+std::string("_")+ringName(ringNumber)+std::string("_")+ECName(ecNumber); }
inline std::string RPC_MultiHistoBase::histoTitleAppend(int ecNumber,unsigned int stationNumber, unsigned int ringNumber)
{ return std::string(" ")+stationName(stationNumber)+ringNameShort(ringNumber)+std::string(" ")+ECName(ecNumber); }
inline  std::string RPC_MultiHistoBase::histoTitleWithoutAppend(std::string appendedTitle,int ecNumber,unsigned int stationNumber, unsigned int ringNumber)
{ return appendedTitle.substr(0,appendedTitle.rfind(histoTitleAppend(ecNumber,stationNumber,ringNumber)));}

inline unsigned int RPC_MultiHistoBase::ecNumberToIndex(int ecNumber) { return (ecNumber+1)/2;}
inline int RPC_MultiHistoBase::indexToEcNumber(unsigned int index) {return 2*int(index)-1;}


inline TH1* RPC_MultiHistoBase::getRingHistoBase(int ecNumber, unsigned int stationNumber, unsigned int ringNumber)
{
  return _hist[ecNumberToIndex(ecNumber)][stationNumber-1][ringNumber-1];
}

inline void RPC_MultiHistoBase::setUse(int ecNumber,unsigned stationNumber,unsigned int ringNumber,bool use)
{
  _useHist[ecNumberToIndex(ecNumber)][stationNumber-1][ringNumber-1]=use;
}

inline RPC_MultiHisto::RPC_MultiHisto(HistoBooker *hb) : RPC_MultiHistoBase(hb) {}

inline TH1F& RPC_MultiHisto::getRingHisto(int ecNumber, unsigned int stationNumber, unsigned int ringNumber)
{ 
  return *( (TH1F*) getRingHistoBase(ecNumber, stationNumber, ringNumber));
}

inline void RPC_MultiHisto::book(std::string name, std::string title, int nbin, float xmin, float xmax)
{
  HistoBooker::Axis X(nbin,xmin,xmax);
  RPC_MultiHistoBase::book(name,title,X);
}

inline void RPC_MultiHisto::Fill(int ecNumber,unsigned stationNumber,unsigned int ringNumber,double t)
{
  getRingHisto(ecNumber, stationNumber, ringNumber).Fill(t);
}

inline RPC_MultiHisto2D::RPC_MultiHisto2D(HistoBooker *hb) : RPC_MultiHistoBase(hb) {}

inline TH2F& RPC_MultiHisto2D::getRingHisto(int ecNumber, unsigned int stationNumber, unsigned int ringNumber)
{ 
  return *( (TH2F*) getRingHistoBase(ecNumber, stationNumber, ringNumber));
}

inline void RPC_MultiHisto2D::book(std::string name, std::string title, int nbinx, float xmin, float xmax, int nbiny, float ymin, float ymax)
{
  HistoBooker::Axis X(nbinx,xmin,xmax);
  HistoBooker::Axis Y(nbiny,ymin,ymax);
  RPC_MultiHistoBase::book(name,title,X,&Y);
}

inline void RPC_MultiHisto2D::Fill(int ecNumber,unsigned stationNumber,unsigned int ringNumber,double x,double y)
{
  getRingHisto(ecNumber, stationNumber, ringNumber).Fill(x,y);
}


#endif
