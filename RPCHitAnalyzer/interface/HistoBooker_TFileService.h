#ifndef HistoBooker_TFileService_HH
#define HistoBooker_TFileService_HH

#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/RPC_MultiHisto.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class HistoBooker_TFileService : public HistoBooker
{
 public:
  HistoBooker_TFileService() : HistoBooker() {}
  virtual bool YouOwn() {return false;}
 private:
  virtual TH1F* book1D(std::string name, std::string title, Axis X);
  virtual TH2F* book2D(std::string name, std::string title, Axis X,Axis Y);
};

#include <iostream>

inline TH1F* HistoBooker_TFileService::book1D(std::string name, std::string title, Axis X)
{
  edm::Service<TFileService> fs;
  return fs->make<TH1F>(name.c_str(), title.c_str(), X.Nbin, X.minValue, X.maxValue);
}

inline TH2F* HistoBooker_TFileService::book2D(std::string name, std::string title, Axis X,Axis Y)
{
  edm::Service<TFileService> fs;
  return fs->make<TH2F>(name.c_str(), title.c_str(), X.Nbin, X.minValue, X.maxValue, Y.Nbin, Y.minValue, Y.maxValue);
} 

#endif
