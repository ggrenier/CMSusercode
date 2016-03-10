



//Local classes
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/ZvertexCalculatorMulti.h"


#include <iostream>

//ROOT minimisation
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"





void ZvertexCalculatorMulti::addHit(hitData &d)
{
  ZvertexCalculator p;
  p.setHitData(d);
  p.setTimeShift(timeshift);
  hitCalcVec.push_back(p);
}

double ZvertexCalculatorMulti::minZ()
{
  double minSol=3e9;
  for (std::vector<ZvertexCalculator>::iterator it=hitCalcVec.begin(); it != hitCalcVec.end(); it++)
    {
      double sol=it->aMinSolution();
	if (minSol>sol) minSol=sol;
    }
  return minSol;
}


double ZvertexCalculatorMulti::meanSolZ()
{
  double mean=0;
  for (std::vector<ZvertexCalculator>::iterator it=hitCalcVec.begin(); it != hitCalcVec.end(); it++)
    mean+=it->aMinSolution();
  return mean/hitCalcVec.size();
}

double ZvertexCalculatorMulti::pullSolZ()
{
  double mean=meanSolZ();
  double shiftTot=0;
  for (std::vector<ZvertexCalculator>::iterator it=hitCalcVec.begin(); it != hitCalcVec.end(); it++)
    {
      double shift=it->aMinSolution()-mean;
      shiftTot+=shift*shift;
    }
  return sqrt(shiftTot);
}

double ZvertexCalculatorMulti::pullSolZBis()
{
  double shiftTot=0;
  for (std::vector<ZvertexCalculator>::iterator it=hitCalcVec.begin(); it != hitCalcVec.end(); it++)
    {
      double thisPos=it->aMinSolution();
      std::vector<ZvertexCalculator>::iterator itb=it; itb++;
      for (;itb != hitCalcVec.end(); itb++)
	{
	  double shift=thisPos-itb->aMinSolution();
	  shiftTot+=shift*shift;
	}
    }
  return sqrt(shiftTot);
}

void ZvertexCalculatorMulti::setTimeShift(double t)
{
  timeshift=t;
  for (std::vector<ZvertexCalculator>::iterator it=hitCalcVec.begin(); it != hitCalcVec.end(); it++)
    it->setTimeShift(t);
}


ZvertexCalculatorMulti* FindZFromPair(hitData d1, hitData d2)
{
  double loopStart=30; //in ns
  
  ZvertexCalculator v1(d1),v2(d2);
  double ts=-loopStart;
  while (ts<=loopStart)
    {
      v1.setTimeShift(ts);
      v2.setTimeShift(ts);
      if (v1.hasSolution() && v2.hasSolution()) break;
      ts+=1;
    }
  double intervalStart=ts;
  ts=loopStart;
  while (ts>=-loopStart)
    {
      v1.setTimeShift(ts);
      v2.setTimeShift(ts);
      if (v1.hasSolution() && v2.hasSolution()) break;
      ts-=1;
    }
  double intervalEnd=ts;
  if (intervalEnd<intervalStart) 
    {
      std::cout << "Can't find non-zero length interval" << std::endl;
      return NULL;
    }

  ZvertexCalculatorMulti *Z=new ZvertexCalculatorMulti() ;
  Z->addHit(d1);
  Z->addHit(d2);

  double ZeroVal=(*Z)(0.0);
  while ((*Z)(intervalStart) < ZeroVal)
    {
      intervalStart-=1.0;
      v1.setTimeShift(intervalStart);
      v2.setTimeShift(intervalStart);
      if ( (!v1.hasSolution()) || (!v2.hasSolution()) || intervalStart<-1000)
	{
	  delete Z;
	  std::cout << "Can't set interval start" << std::endl;
	  return NULL;
	}
    }

  while ((*Z)(intervalEnd) < ZeroVal)
    {
      intervalEnd+=1.0;
      v1.setTimeShift(intervalEnd);
      v2.setTimeShift(intervalEnd);
      if ( (!v1.hasSolution()) || (!v2.hasSolution()) || intervalEnd>1000)
	{
	  delete Z;
	  std::cout << "Can't set interval end" << std::endl;
	  return NULL;
	}
    }

  
  std::cout << " start = " << intervalStart;
  std::cout << " stop = " <<  intervalEnd << std::endl;

  //std::cout<< "Precheck "<<std::endl;
  //for (double v=intervalStart; v<intervalEnd ; v+=1.0) (*Z)(v);

  ROOT::Math::Functor1D func(*Z);

  ROOT::Math::GSLMinimizer1D minBrent;
  minBrent.SetFunction(func,0,intervalStart,intervalEnd); //start minimizing at 0, initial search interval
  minBrent.Minimize(100,0.001,0.01); // max 100 tries or stop if error=0.001 (1 ps) or relative error=0.01 (1%)

  double timeShift=minBrent.XMinimum();
  Z->setTimeShift(timeShift);

  return Z;
}





