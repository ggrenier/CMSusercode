

#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/hitDataPair.h"
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/RPC_MCNavigator.h" //for cspeed

#include <iostream>


hitDataPair::hitDataPair(hitData &hit1, hitData &hit2)
  : d1(hit1), d2(hit2)
{
  ComputeDistance();
  ComputeTime();
  std::cout << "time in ns (before,after,expect)=("
	    << timeBefore << "," << timeAfter << "," << timeExpect << ") ";
}


void hitDataPair::ComputeDistance()
{
  distance=0;
  distance+=(d1.x-d2.x)*(d1.x-d2.x);
  distance+=(d1.y-d2.y)*(d1.y-d2.y);
  distance+=(d1.z-d2.z)*(d1.z-d2.z);
  distance=sqrt(distance);
}

void hitDataPair::ComputeTime()
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



