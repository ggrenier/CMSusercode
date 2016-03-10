
#include <cmath>

#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/ZvertexCalculator.h" 
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/RPC_MCNavigator.h" //For RPC::cspeed


void ZvertexCalculator::setHitData(hitData &hit)
{
  region=hit.region;
  b=hit.t*RPC::cspeed();
  c=hit.r;
  hit_time=hit.t;
  double costheta=hit.z/hit.r;
  cosbeta=computeSignOfCos()*costheta;
  sinbeta=sqrt(1-cosbeta*cosbeta);
  //std::cout << "Zvtx init " << b << " " << c << " " << region << " " << std::endl;
}

void ZvertexCalculator::setTimeShift(double timeshift)
{
  int signBefore=computeSignOfCos();
  b=(hit_time-timeshift)*RPC::cspeed();
  int signAfter=computeSignOfCos();
  int signChange=signBefore*signAfter; 
  cosbeta=signChange*cosbeta; //no need to change sinbeta
}

double ZvertexCalculator::aMinSolution()
{
  if (! hasSolution()) return 3000000;
  double sol=aOneSolution();
  if (hasTwoSolution()) 
    {
      double sol2=aSecondSolution();
      if (fabs(sol2)<fabs(sol)) sol=sol2;
    }
  if ((region==1 && b<c)||(region==-1 && b>c)) return sol;
  return -sol;
}

