#ifndef hitDataPair_HH
#define hitDataPair_HH


//Local classes
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/ZvertexCalculatorMulti.h"



class hitDataPair
{
public:
  hitDataPair(hitData &hit1, hitData &hit2);

  inline bool diffEncap() { return d1.region != d2.region;}
  inline bool lightCone() { return timeAfter==timeExpect;}
  inline bool lightCone(double timeError) { return fabs(timeAfter-timeExpect)<timeError;}

  inline ZvertexCalculatorMulti* getVertex() { return FindZFromPair(d1,d2); }

  inline int station1() {return d1.station;}
  inline int station2() {return d2.station;}
  
private:
  hitData d1;
  hitData d2;
  double distance;
  double timeBefore;
  double timeAfter;
  double timeExpect;

  void ComputeDistance();
  void ComputeTime();

};


#endif
