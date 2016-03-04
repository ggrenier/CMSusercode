//It is assumed you have done .L src/RPC_MultiHisto.cc in the parent directory 
// and that you're going to do .L scripts/Utilities.C
//
// this is interpreted code by cint
//

struct dataAccess
{
public:
  TFile *F;
  TDirectory *d;
  HistoRetriever *hr;
  dataAccess(const char* filename, const char *directory="testee");
  dataAccess(TFile* afile, const char *directory="testee");
  void load(const char *directory); 
  bool checkOK();
  void print();
  RPC_MultiHisto* getMultiHist(const char *name);
  TH1F* getTH1F(const char *name);
  TH2F* getTH2F(const char *name);
};


TLegend *DrawTLegendCanvas(TCanvas *c,RPC_MultiHisto *rp)
{
  c->cd();
  TLegend *l=rp->getLegend(0, 0, 1, 1);
  l->Draw();
  return l;
}

dataAccess::dataAccess(const char* filename, const char *directory)
{
  F=TFile::Open(filename);
  load(directory);
}

dataAccess::dataAccess(TFile* afile, const char *directory)
{
  F=afile;
  load(directory);
}

void dataAccess::load(const char *directory)
{
  d=NULL;
  F->GetObject(directory,d);
  if (d) hr=new HistoRetriever(*d);
  else hr=NULL;
  if (checkOK()) print();
  else cout << "problem with loading" << endl;
}

bool dataAccess::checkOK() 
{
  return (NULL != d && NULL != hr);
}

void dataAccess::print()
{
  TList* keyList=d->GetListOfKeys();
  cout << keyList->GetSize() << endl;
  TIter next(keyList);
  TKey* named=NULL;
  TRegexp match("_RE[1-4]_Ring[1-3]_EC");
  while ((named=(TKey*)next()))
    {
      TString sname=named->GetName();
      TString stitle=named->GetTitle();
      if (sname.Contains(match))
	{
	  if (sname.EndsWith("_RE1_Ring1_ECminus"))
	    {
	      sname.Remove(sname.Index("_RE1_Ring1_ECminus"));
              stitle.Remove(stitle.Index(" RE1/1 ECminus"));
	      cout << "RPC_MultiHisto : " << sname <<  " (" << named->GetClassName() <<"*) " << stitle << endl;
	    }
	}
      else
	cout << sname << " (" << named->GetClassName() <<"*) " << named->GetTitle() << endl;
    }
}


RPC_MultiHisto* dataAccess::getMultiHist(const char *name)
{
  RPC_MultiHisto* rhm=new RPC_MultiHisto(hr);
  rhm->book(name,"dummy",0,0,0);
  return rhm;
}

TH1F* dataAccess::getTH1F(const char *name)
{
  TH1F *h1=NULL;
  d->GetObject(name,h1);
  return h1;
}

TH2F* dataAccess::getTH2F(const char *name)
{
  TH2F *h2=NULL;
  d->GetObject(name,h2);
  return h2;
}
