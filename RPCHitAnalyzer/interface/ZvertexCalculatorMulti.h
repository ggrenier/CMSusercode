#ifndef ZvertexCalculatorMulti_HH
#define ZvertexCalculatorMulti_HH


//Local classes
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/ZvertexCalculator.h"

#include <vector>


class ZvertexCalculatorMulti
{
public:
  ZvertexCalculatorMulti() : timeshift(0),hitCalcVec() {}
  void addHit(hitData &d);
  double minZ();
  double meanSolZ();
  double pullSolZ();
  double pullSolZBis();
  void setTimeShift(double t);
  inline double operator()(double t)
  {
    setTimeShift(t);
    double p=pullSolZBis();
    //std::cout << "function call  t= " << t << " answer= " << p << std::endl;
    return p;
  }
  inline double getTimeShift() {return timeshift;}
private:
  double timeshift;
  std::vector<ZvertexCalculator> hitCalcVec;
};

ZvertexCalculatorMulti* FindZFromPair(hitData d1, hitData d2);


#endif


