//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 10 10:29:32 2014 by ROOT version 5.34/03
// from TTree HitTree/CMS RPC hits
// found on file: ../../Data/histos/histo_simhit_SingleMu_upscope_Pt100_Zvtx0_ZvtxSigma60_etamax2.6_etamin1.3_10000.root
//////////////////////////////////////////////////////////

#ifndef ZVtxStudy_h
#define ZVtxStudy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
//#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.


class MyStudy;

class ZVtxStudy {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        PVZ;
   Double_t        PVCT;
   std::vector<int>     *pid;
   std::vector<int>     *region;
   std::vector<int>     *station;
   std::vector<int>     *ring;
   std::vector<float>   *x;
   std::vector<float>   *y;
   std::vector<float>   *z;
   std::vector<float>   *t;
   std::vector<float>   *localdirEta;
   std::vector<float>   *localdirPhi;
   std::vector<float>   *localdirTheta;
   std::vector<float>   *genMuonPx;
   std::vector<float>   *genMuonPy;
   std::vector<float>   *genMuonPz;

   // List of branches
   TBranch        *b_PVZ;   //!
   TBranch        *b_PVCT;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_region;   //!
   TBranch        *b_station;   //!
   TBranch        *b_ring;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_t;   //!
   TBranch        *b_localdirEta;   //!
   TBranch        *b_localdirPhi;   //!
   TBranch        *b_localdirTheta;   //!
   TBranch        *b_genMuonPx;   //!
   TBranch        *b_genMuonPy;   //!
   TBranch        *b_genMuonPz;   //!

   ZVtxStudy(TTree *tree=0);
   ZVtxStudy(const char* filename);
   virtual ~ZVtxStudy();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(MyStudy* study);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ZVtxStudy_cxx
ZVtxStudy::ZVtxStudy(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../Data/histos/histo_simhit_SingleMu_upscope_Pt100_Zvtx0_ZvtxSigma60_etamax2.6_etamin1.3_10000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../Data/histos/histo_simhit_SingleMu_upscope_Pt100_Zvtx0_ZvtxSigma60_etamax2.6_etamin1.3_10000.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../../Data/histos/histo_simhit_SingleMu_upscope_Pt100_Zvtx0_ZvtxSigma60_etamax2.6_etamin1.3_10000.root:/testee");
      dir->GetObject("HitTree",tree);

   }
   Init(tree);
}

ZVtxStudy::ZVtxStudy(const char* filename) : fChain(0) 
{
  TTree *tree=0;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f || !f->IsOpen()) {
         f = new TFile(filename);
      }
      TDirectory * dir = (TDirectory*)f->GetDirectory("testee");
      dir->GetObject("HitTree",tree);

   }
   Init(tree);
}


ZVtxStudy::~ZVtxStudy()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ZVtxStudy::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ZVtxStudy::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ZVtxStudy::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pid = 0;
   region = 0;
   station = 0;
   ring = 0;
   x = 0;
   y = 0;
   z = 0;
   t = 0;
   localdirEta = 0;
   localdirPhi = 0;
   localdirTheta = 0;
   genMuonPx = 0;
   genMuonPy = 0;
   genMuonPz = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PVZ", &PVZ, &b_PVZ);
   fChain->SetBranchAddress("PVCT", &PVCT, &b_PVCT);
   fChain->SetBranchAddress("pid", &pid, &b_pid);
   fChain->SetBranchAddress("region", &region, &b_region);
   fChain->SetBranchAddress("station", &station, &b_station);
   fChain->SetBranchAddress("ring", &ring, &b_ring);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("localdirEta", &localdirEta, &b_localdirEta);
   fChain->SetBranchAddress("localdirPhi", &localdirPhi, &b_localdirPhi);
   fChain->SetBranchAddress("localdirTheta", &localdirTheta, &b_localdirTheta);
   fChain->SetBranchAddress("genMuonPx", &genMuonPx, &b_genMuonPx);
   fChain->SetBranchAddress("genMuonPy", &genMuonPy, &b_genMuonPy);
   fChain->SetBranchAddress("genMuonPz", &genMuonPz, &b_genMuonPz);
   Notify();
}

Bool_t ZVtxStudy::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZVtxStudy::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ZVtxStudy::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  if (entry<0) return -1;
   return 1;
}
#endif // #ifdef ZVtxStudy_cxx
