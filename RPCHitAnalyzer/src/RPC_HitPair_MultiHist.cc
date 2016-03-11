


//For plots and histos
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TH1F.h>


//Local classes
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/RPC_HitPair_MultiHist.h"






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

