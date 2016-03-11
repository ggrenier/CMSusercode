#ifndef RPC_HitPair_MultiHist_HH
#define RPC_HitPair_MultiHist_HH


#include <string>
#include <vector>

//For plots and histos
#include <TH1F.h>

//local for plots and histos
//#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/RPC_MultiHisto.h"
//#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/HistoBooker_TFileService.h"

//Local classes
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/hitDataPair.h"


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
  inline const char* extend(std::string val) {return (val+"_"+extension).c_str();}
  inline std::string histoTitle(std::string precision) {std::string s("Estimated PV Z - true PV Z (in cm) "); return s+precision;}
  int indice(int station1, int station2);
};

#endif
