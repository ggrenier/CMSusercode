#ifndef __CINT__
#include "CMSRPCDPGUserCode/RPCHitAnalyzer/interface/RPC_MultiHisto.h"
#else
#include "interface/RPC_MultiHisto.h"
#endif




RPC_MultiHistoBase::RPC_MultiHistoBase(HistoBooker *hb) : _defaultBooker() 
{
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	{
	  _hist[iec][istation][iring]=NULL;
	  _useHist[iec][istation][iring]=true;
	}
  if (hb != NULL) setHistoBooker(hb);
  else setHistoBooker(&_defaultBooker);
}

RPC_MultiHistoBase::~RPC_MultiHistoBase()
{
  if (_IOwnHistos)
    for (unsigned int iec=0; iec<NECs; iec++)
      for (unsigned int istation=0; istation<Nstations; istation++)
	for (unsigned int iring=0; iring<Nrings; iring++)
	  if (_hist[iec][istation][iring]!=NULL)
	    delete _hist[iec][istation][iring];
}


void RPC_MultiHistoBase::Write()
{
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
        _hist[iec][istation][iring]->Write();
}


void RPC_MultiHistoBase::book(std::string name, std::string title, HistoBooker::Axis X, HistoBooker::Axis *Y)
{
  _histName=name;
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	_hist[iec][istation][iring]=_currentBooker->book(name+histoNameAppend(indexToEcNumber(iec),istation+1,iring+1),
							 title+histoTitleAppend(indexToEcNumber(iec),istation+1,iring+1),
							 X,Y);
}



void RPC_MultiHistoBase::setUseAll(bool use)
{
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	_useHist[iec][istation][iring]=use;
}

void RPC_MultiHistoBase::setUseRing(unsigned int ringNumber,bool use)
{
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	if (iring==ringNumber-1)
	  _useHist[iec][istation][iring]=use;      
}

#ifndef __CINT__
std::string RPC_MultiHistoBase::usehitToShortString()
{
  std::string s("_");
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	s+=(_useHist[iec][istation][iring] ? '1' : '0'); 
  return s;
}
#endif

THStack* RPC_MultiHistoBase::getStack()
{
#ifdef __CINT__
  std::string slocal("_");
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	slocal+=(_useHist[iec][istation][iring] ? '1' : '0'); 
#endif
  THStack *stack=NULL;
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	if (_useHist[iec][istation][iring])
	  {
	    if (NULL==stack) 
	      {
		std::string stackName=_histName;
#ifndef __CINT__
		stackName+=usehitToShortString();
#else
		stackName+=slocal;
#endif
		std::string stackTitle=histoTitleWithoutAppend(_hist[iec][istation][iring]->GetTitle(),indexToEcNumber(iec),istation+1,iring+1);
		stackTitle+=" Stack";
		stack=new THStack(stackName.c_str(),stackTitle.c_str());
	      }
	    stack->Add(_hist[iec][istation][iring]);
	  }
  return stack;
}


TH1* RPC_MultiHistoBase::getHistoSum()
{
//NB : I haven't done _hist[iec][istation][iring]->Sumw2(). That may be added if needed
#ifdef __CINT__
  std::string slocal2("_");
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	slocal2+=(_useHist[iec][istation][iring] ? '1' : '0'); 
#endif
  TH1* r=NULL;
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	  if (_useHist[iec][istation][iring])
	    {
	      if (NULL==r)
		{
		  std::string sumName=_histName+"_sum";
#ifndef __CINT__
		  sumName+=usehitToShortString();
#else
		  sumName+=slocal2;
#endif
		  r=(TH1*) _hist[iec][istation][iring]->Clone( sumName.c_str() );
		  r->SetTitle((histoTitleWithoutAppend(r->GetTitle(),indexToEcNumber(iec),istation+1,iring+1)+" Summed").c_str());
		}
	      else r->Add( _hist[iec][istation][iring]);
	    }
  return r; 
}


TLegend* RPC_MultiHistoBase::getLegend(double x1, double y1, double x2, double y2, const char * header, const char *option)
{
  TLegend *leg=new TLegend(x1,y1,x2,y2,header,option);
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
        if (_useHist[iec][istation][iring])
          leg->AddEntry(_hist[iec][istation][iring],histoTitleAppend(indexToEcNumber(iec),istation+1,iring+1).c_str());
  return leg;
}

void RPC_MultiHistoBase::setLineAttributes()
{
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	{
	  _hist[iec][istation][iring]->SetLineColor(istation+1);
	  _hist[iec][istation][iring]->SetLineStyle(iring+1+7*iec);
	}
}

void RPC_MultiHistoBase::setFillAttributes()
{
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	{
	  _hist[iec][istation][iring]->SetFillColor(istation+1);
	  _hist[iec][istation][iring]->SetFillStyle(3004+iring+3*iec);
	}
}

void RPC_MultiHistoBase::setMarkerAttributes()
{
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	{
	  _hist[iec][istation][iring]->SetMarkerColor(istation+1);
	  _hist[iec][istation][iring]->SetMarkerSize(0.5);
	  _hist[iec][istation][iring]->SetMarkerStyle(20+iring+4*iec);
	}
}


void RPC_MultiHistoBase::setAllAttributes()
{
  setLineAttributes();
  setFillAttributes();
  setMarkerAttributes();
}

void RPC_MultiHistoBase::resetLineAttributes()
{
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	_hist[iec][istation][iring]->ResetAttLine();
}

void RPC_MultiHistoBase::resetFillAttributes()
{
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	_hist[iec][istation][iring]->ResetAttFill();
}

void RPC_MultiHistoBase::resetMarkerAttributes()
{
  for (unsigned int iec=0; iec<NECs; iec++)
    for (unsigned int istation=0; istation<Nstations; istation++)
      for (unsigned int iring=0; iring<Nrings; iring++)
	_hist[iec][istation][iring]->ResetAttMarker();
}


