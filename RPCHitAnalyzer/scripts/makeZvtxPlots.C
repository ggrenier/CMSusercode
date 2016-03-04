#include "ZVtxStudy.C"
#include "TGraph.h"
#include "TLegend.h"
#include "TFitResult.h"

void makeZvtxPlots()
{
  //gSystem->Load("ZVtxStudy_C.so");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);


  myAnalysisSet AS;
  
  AS.add(3.0,  "../../Data/histos/histo_simhit_SingleMu_upscope_Pt3_Zvtx0_ZvtxSigma60_etamax2.6_etamin1.3_10000.root");
  AS.add(5.0,  "../../Data/histos/histo_simhit_SingleMu_upscope_Pt5_Zvtx0_ZvtxSigma60_etamax2.6_etamin1.3_10000.root");
  AS.add(10.0, "../../Data/histos/histo_simhit_SingleMu_upscope_Pt10_Zvtx0_ZvtxSigma60_etamax2.6_etamin1.3_10000.root");
  AS.add(20.0, "../../Data/histos/histo_simhit_SingleMu_upscope_Pt20_Zvtx0_ZvtxSigma60_etamax2.6_etamin1.3_10000.root");
  AS.add(50.0, "../../Data/histos/histo_simhit_SingleMu_upscope_Pt50_Zvtx0_ZvtxSigma60_etamax2.6_etamin1.3_10000.root");
  AS.add(100.0,"../../Data/histos/histo_simhit_SingleMu_upscope_Pt100_Zvtx0_ZvtxSigma60_etamax2.6_etamin1.3_10000.root");

  MultiMyStudy stud;
  stud.addStudy(new HitPairStudySmearTimeAndStrip(0.1),"100psTandRsmearing"); //100 ps  
  stud.addStudy(new HitPairStudySmearTimeAndStrip(0.05),"50psTandRsmearing"); //50 ps
  stud.addStudy(new HitPairStudySmearTimeAndStrip(0.02),"20psTandRsmearing"); //20 ps
  stud.addStudy(new HitPairStudySmearTimeAndStrip(0.01),"10psTandRsmearing"); //10 ps
  stud.addStudy(new HitPairStudySmearTimeAndStrip(1),"1nsTandRsmearing"); //1 ns
  stud.addStudy(new HitPairStudyWithMorePairs(1,0.1),"100psTandRsmearing_2HitsPerLayer");


  AS.setStudy(stud);

  unsigned Nfiles=AS.nStudies();
  MultiMyStudy **studies=(MultiMyStudy**)AS.theStudies();
  float* Pts=AS.thePts();

  studies[0]->printStudyNames();


  unsigned int nStudies=studies[0]->nStudies();

  TCanvas c;
  TGraph *PairHitsResolVsPt[nStudies]; float mingrpair=1e8; float maxgrpair=0;
  TGraph *AllHitsResolVsPt[nStudies];
  

  for (unsigned int istud=0; istud<nStudies; istud++)
    {
      c.Clear();    
      c.Divide((Nfiles+1)/2,2);
      PairHitsResolVsPt[istud]=new TGraph(Nfiles);
      PairHitsResolVsPt[istud]->SetLineColor(istud+1);
      PairHitsResolVsPt[istud]->SetMarkerColor(istud+1);
      PairHitsResolVsPt[istud]->SetMarkerStyle(istud+20);
      AllHitsResolVsPt[istud]=new TGraph(Nfiles);
      AllHitsResolVsPt[istud]->SetLineColor(istud+1);
      AllHitsResolVsPt[istud]->SetMarkerColor(istud+1);
      AllHitsResolVsPt[istud]->SetMarkerStyle(istud+20);
      
      for (unsigned int i=0; i<Nfiles; ++i)
	{
	  c.cd(i+1);
	  HitPairStudy* hp=(HitPairStudy*) studies[i]->getStudy(istud);
	  TFitResult& fitp=*(hp->diffEC.diffPVZ[2]->Fit("gaus","S"));
	  hp->diffEC.diffPVZ[2]->Draw();
	  const double *fitval=fitp.GetParams();
	  PairHitsResolVsPt[istud]->SetPoint(i,Pts[i],fitval[2]);
	  if (mingrpair>fitval[2]) mingrpair=fitval[2];
	  if (maxgrpair<fitval[2]) maxgrpair=fitval[2];
	  if (i==0) PairHitsResolVsPt[istud]->SetTitle(hp->diffEC.diffPVZ[2]->GetTitle());
	}
      c.SaveAs((studies[0]->studyName(istud)+"_diffEC2hits.pdf").c_str());

      
      for (unsigned int i=0; i<Nfiles; ++i)
	{
	  c.cd(i+1);
	  HitPairStudy* hp=(HitPairStudy*) studies[i]->getStudy(istud);
	  TFitResult& fitp=*(hp->diffdiffPVZ8hits->Fit("gaus","S"));
	  hp->diffdiffPVZ8hits->Draw();
	  const double *fitval=fitp.GetParams();
	  AllHitsResolVsPt[istud]->SetPoint(i,Pts[i],fitval[2]);	  
	  if (i==0) AllHitsResolVsPt[istud]->SetTitle(hp->diffdiffPVZ8hits->GetTitle());

	  //cout << " Nparameters = " << fitp.NPar();
	  //const double *fitval=fitp.GetParams();
	  //for (int ip=0; ip<fitp.NPar(); ip++) cout << "  fit par "<<ip<<" "<<fitp.GetParameterName(ip)
	  //					<< " value " << fitval[ip];
	  //cout << endl;  
	}

      c.SaveAs((studies[0]->studyName(istud)+"_diffECAllHits.pdf").c_str());

    }
  
  c.Clear();
  c.Divide(2,1);
  c.cd(1);
  TLegend *tl=new TLegend(0.2,0.44,0.99,0.73,"using pair of hit");
  PairHitsResolVsPt[0]->SetMaximum(3.5);
  PairHitsResolVsPt[0]->SetMinimum(0);
  PairHitsResolVsPt[0]->Draw("ALP");
  PairHitsResolVsPt[0]->GetXaxis()->SetTitle("muon Pt (GeV)");
  PairHitsResolVsPt[0]->GetYaxis()->SetTitle("reco PVZ error (cm)");
  tl->AddEntry(PairHitsResolVsPt[0],studies[0]->studyName(0).c_str(),"lp");
  for (unsigned int istud=1; istud<nStudies; istud++)
    {
      if (istud==4) continue;
      PairHitsResolVsPt[istud]->Draw("LP");
      tl->AddEntry(PairHitsResolVsPt[istud],studies[0]->studyName(istud).c_str(),"lp");
    }
  tl->Draw();
  cout << mingrpair << " ---- " << maxgrpair << endl;
  //c.SaveAs("PairHitsResolVsPt.pdf");

  c.cd(2);
  //c.Clear();
  TLegend *tl2=new TLegend(0.2,0.5,0.99,0.8,"using all pairs of hit");
  AllHitsResolVsPt[0]->SetMaximum(3.5);
  AllHitsResolVsPt[0]->SetMinimum(0);
  AllHitsResolVsPt[0]->Draw("ALP");
  AllHitsResolVsPt[0]->GetXaxis()->SetTitle("muon Pt (GeV)");
  AllHitsResolVsPt[0]->GetYaxis()->SetTitle("reco PVZ error (cm)");
  tl2->AddEntry(AllHitsResolVsPt[0],studies[0]->studyName(0).c_str(),"lp");
  for (unsigned int istud=1; istud<nStudies; istud++)
    {
      if (istud==4) continue;
      AllHitsResolVsPt[istud]->Draw("LP");
      tl2->AddEntry(AllHitsResolVsPt[istud],studies[0]->studyName(istud).c_str(),"lp");
    }
  tl2->Draw();
  c.SaveAs("HitsResolVsPt.pdf");

  TGraph *PairHitsResolVsTS[Nfiles];
  float TS[5]={100,50,20,10,1000};
  for (unsigned int ipt=0; ipt<Nfiles; ipt++) 
    {
      PairHitsResolVsTS[ipt]=new TGraph(nStudies-1); //nstudies-1==5
      PairHitsResolVsTS[ipt]->SetLineColor(ipt+1);
      PairHitsResolVsTS[ipt]->SetMarkerColor(ipt+1);
      PairHitsResolVsTS[ipt]->SetMarkerStyle(ipt+20);
      PairHitsResolVsTS[ipt]->SetTitle("using pair of hits");
      
      for (unsigned int its=0; its<5; its++)
	{
	  double* lesResol=PairHitsResolVsPt[its]->GetY();
	  //cout << "pointeur : " << lesResol << endl; 
	  PairHitsResolVsTS[ipt]->SetPoint(its,TS[its],lesResol[ipt]);
	}
    }
  c.Clear();
  PairHitsResolVsTS[0]->Sort();
  PairHitsResolVsTS[0]->Draw("ALP");
  PairHitsResolVsTS[0]->GetXaxis()->SetTitle("RPC hit time smearing (ps)");
  PairHitsResolVsTS[0]->GetYaxis()->SetTitle("reco PVZ error (cm)");
  
  TLegend *tl3=new TLegend(0.7,0.15,0.9,0.6,"muon Pt in GeV");
  TString a; a+=Pts[0];
  tl3->AddEntry(PairHitsResolVsTS[0],a,"lp");
  for (unsigned int ipt=1; ipt<Nfiles; ipt++)
    {
       PairHitsResolVsTS[ipt]->Sort();
       PairHitsResolVsTS[ipt]->Draw("LP");
       a="";
       a+=Pts[ipt];
       tl3->AddEntry(PairHitsResolVsTS[ipt],a,"lp");
    }
  tl3->Draw();
  c.SaveAs("PairHitsResolVsTS.pdf");

  //for (unsigned int ipt=0; ipt<Nfiles; ipt++)
  //PairHitsResolVsTS[ipt].RemovePoint(4);
  
  //c.Update();

  myAnalysisSet ASnoSmear;
  
  ASnoSmear.add(3.0,  "../../Data/histos/histo_simhit_SingleMu_upscope_Pt3_Zvtx0_etamax2.6_etamin1.3_10000.root");
  ASnoSmear.add(5.0,  "../../Data/histos/histo_simhit_SingleMu_upscope_Pt5_Zvtx0_etamax2.6_etamin1.3_10000.root");
  ASnoSmear.add(10.0, "../../Data/histos/histo_simhit_SingleMu_upscope_Pt10_Zvtx0_etamax2.6_etamin1.3_10000.root");
  ASnoSmear.add(20.0, "../../Data/histos/histo_simhit_SingleMu_upscope_Pt20_Zvtx0_etamax2.6_etamin1.3_10000.root");
  ASnoSmear.add(50.0, "../../Data/histos/histo_simhit_SingleMu_upscope_Pt50_Zvtx0_etamax2.6_etamin1.3_10000.root");
  ASnoSmear.add(100.0,"../../Data/histos/histo_simhit_SingleMu_upscope_Pt100_Zvtx0_etamax2.6_etamin1.3_10000.root");

  
  MultiMyStudy studMix;
  studMix.addStudy(new HitPairStudySmearTimeAndStrip(0.01,0,1.0,new HitPairStudyEventManagerMixECplusNevents),"10psTandRsmearing"); //10 ps
  studMix.addStudy(new HitPairStudySmearTimeAndStrip(0.02,0,1.0,new HitPairStudyEventManagerMixECplusNevents),"20psTandRsmearing"); //20 ps
  studMix.addStudy(new HitPairStudySmearTimeAndStrip(0.05,0,1.0,new HitPairStudyEventManagerMixECplusNevents),"50psTandRsmearing"); //50 ps
  studMix.addStudy(new HitPairStudySmearTimeAndStrip(0.1 ,0,1.0,new HitPairStudyEventManagerMixECplusNevents),"100psTandRsmearing"); //100 ps  
  studMix.addStudy(new HitPairStudySmearTime(1,0,1.0,new HitPairStudyEventManagerMixECplusNevents),"1nsTsmearing"); //1 ns


  AS.setStudy(studMix);
  ASnoSmear.setStudy(studMix);
  
  unsigned int nMixStudies=studMix.nStudies();

  studies=(MultiMyStudy**)AS.theStudies(); //retrieve the new study lists.

  MultiMyStudy **studnoSmear=(MultiMyStudy**)ASnoSmear.theStudies();
  unsigned NfilesnoSmear=ASnoSmear.nStudies();
  Pts=ASnoSmear.thePts(); 
  
  gStyle->SetOptStat(0);

  c.Clear();
  c.Divide((Nfiles+1)/2,2);
  for (unsigned int i=0; i<Nfiles; ++i)
    {
      c.cd(i+1);
      HitPairStudy* hp=(HitPairStudy*) studies[i]->getStudy(0);
      hp->correctedTOFdiffsameStation[3]->Draw();
      for (unsigned int istud=1; istud<nMixStudies; istud++)
	{
	  hp=(HitPairStudy*) studies[i]->getStudy(istud);
	  hp->correctedTOFdiffsameStation[3]->SetLineColor(istud+istud/4);
	  hp->correctedTOFdiffsameStation[3]->Draw("SAME");
	}
    }
  c.SaveAs("DiffTOFmuonFromDiffVertex.pdf");

  c.Clear();
  c.Divide((NfilesnoSmear+1)/2,2);
  for (unsigned int i=0; i<NfilesnoSmear; ++i)
    {
      c.cd(i+1);
      HitPairStudy* hp=(HitPairStudy*) studnoSmear[i]->getStudy(0);
      hp->correctedTOFdiffsameStation[3]->Draw();
      for (unsigned int istud=1; istud<nMixStudies; istud++)
	{
	  hp=(HitPairStudy*) studnoSmear[i]->getStudy(istud);
	  hp->correctedTOFdiffsameStation[3]->SetLineColor(istud+istud/4);
	  hp->correctedTOFdiffsameStation[3]->Draw("SAME");
	}
    }
  c.SaveAs("DiffTOFmuonFromSameVertex.pdf");

  c.Clear();
  TLegend *tl4=new TLegend(0.6,0.5,0.9,0.9,"RPC time resolution");

  HitPairStudy* hp=(HitPairStudy*) studies[Nfiles-1]->getStudy(0);
  hp->correctedTOFdiffsameStation[3]->Draw();
  tl4->AddEntry(hp->correctedTOFdiffsameStation[3],studies[Nfiles-1]->studyName(0).c_str(),"lp");

  for (unsigned int istud=1; istud<nMixStudies; istud++)
    {
      hp=(HitPairStudy*) studies[Nfiles-1]->getStudy(istud);
      hp->correctedTOFdiffsameStation[3]->SetLineColor(istud+istud/4);
      hp->correctedTOFdiffsameStation[3]->Draw("SAME");
      tl4->AddEntry(hp->correctedTOFdiffsameStation[3],studies[Nfiles-1]->studyName(istud).c_str(),"lp");
    }
  tl4->Draw();
  c.SaveAs("DiffTOFmuonFromDiffVertexPt100.pdf");

  hp=(HitPairStudy*) studnoSmear[NfilesnoSmear-1]->getStudy(0);
  hp->correctedTOFdiffsameStation[3]->Draw();
  for (unsigned int istud=1; istud<nMixStudies; istud++)
    {
      hp=(HitPairStudy*) studnoSmear[NfilesnoSmear-1]->getStudy(istud);
      hp->correctedTOFdiffsameStation[3]->SetLineColor(istud+istud/4);
      hp->correctedTOFdiffsameStation[3]->Draw("SAME");
    }
  tl4->Draw();
  c.SaveAs("DiffTOFmuonFromSameVertexPt100.pdf");

  //now I assume you have the same Pt for files in AS and ASnoSmear
  TGraph* ROC[Nfiles][nMixStudies];
  std::cout << "Pt \t study             \t RMS_sameVTX  \t  RMS_diffVTX" << std::endl;
  if (Nfiles != NfilesnoSmear)
    cout << "FILE LIST NOT COHERENT" << endl;
  else
    {
      for (unsigned int i=0; i<Nfiles; ++i)
	{
	  c.Clear();
	  c.Divide((nMixStudies+1)/2,2);
	  for (unsigned int istud=0; istud<nMixStudies; istud++)
	    {
	      std::cout << Pts[i] << "\t" << studnoSmear[i]->studyName(istud) << "\t";
	      c.cd(istud+1);
	      TLegend *legloc=new TLegend(0.5,0.7,0.9,0.9,studnoSmear[i]->studyName(istud).c_str());
	      hp=(HitPairStudy*) studnoSmear[i]->getStudy(istud);
	      hp->correctedTOFdiffsameStation[3]->SetLineStyle(2);
	      hp->correctedTOFdiffsameStation[3]->Draw();
	      legloc->AddEntry(hp->correctedTOFdiffsameStation[3],"2 muons from same vertex","lp");
	      std::cout << hp->correctedTOFdiffsameStation[3]->GetRMS() << "\t";
	      double *y_noSmear=hp->correctedTOFdiffsameStation[3]->GetIntegral();
	      hp=(HitPairStudy*) studies[i]->getStudy(istud);
	      hp->correctedTOFdiffsameStation[3]->Draw("SAME");
	      legloc->AddEntry(hp->correctedTOFdiffsameStation[3],"2 muons from diff vertex","lp");
	      std::cout << hp->correctedTOFdiffsameStation[3]->GetRMS() << std::endl;
	      double *x_smear=hp->correctedTOFdiffsameStation[3]->GetIntegral();
	      legloc->Draw();
	      //Double_t *f=hp->correctedTOFdiffsameStation[3]->GetIntegral();
	      //for (int k=0; k<hp->correctedTOFdiffsameStation[3]->GetNbinsX()+2;++k)
	      //std::cout << f[k] << ":";
	      //std::cout<<std::endl;
	      ROC[i][istud]=new TGraph(hp->correctedTOFdiffsameStation[3]->GetNbinsX()+1,x_smear,y_noSmear);
	      ROC[i][istud]->GetXaxis()->SetTitle("Fraction kept hit pairs from muons from different vertices");
	      ROC[i][istud]->GetYaxis()->SetTitle("Fraction kept hit pairs from muons originating from the same vertex");
	      TString graphTitle=studnoSmear[i]->studyName(istud).c_str();
	      graphTitle+="_Pt="; graphTitle+=Pts[i]; graphTitle+="GeV";
	      ROC[i][istud]->SetTitle(graphTitle);
	      ROC[i][istud]->SetLineColor(1+istud+istud/4);
	      ROC[i][istud]->SetLineStyle(i+1);
	    }
	  TString s("DiffTOFmuonPt"); s+=Pts[i]; s+=".pdf";
	  c.SaveAs(s);

	  
	  c.Clear();
	  c.Divide((nMixStudies+1)/2,2);
	  //gStyle->SetOptStat(111111);
	  for (unsigned int istud=0; istud<nMixStudies; istud++)
	    {
	      std::cout << Pts[i] << "\t" << studnoSmear[i]->studyName(istud) << "\t";
	      c.cd(istud+1);
	      TLegend *legloc=new TLegend(0.5,0.7,0.9,0.9,studnoSmear[i]->studyName(istud).c_str());
	      hp=(HitPairStudy*) studies[i]->getStudy(istud);
	      if (istud != 0) hp->correctedTOFdiffsameStation4_TrueDiffPVZless2mm->SetLineColor(istud+istud/4);
	      hp->correctedTOFdiffsameStation4_TrueDiffPVZless2mm->Draw();
	      legloc->AddEntry(hp->correctedTOFdiffsameStation4_TrueDiffPVZless2mm,"2 muons from diff vertex Delta PVZ < 2mm","lp");
	      float nentries=hp->correctedTOFdiffsameStation4_TrueDiffPVZless2mm->GetEntries();
	      hp=(HitPairStudy*) studnoSmear[i]->getStudy(istud);	      
	      /////////////////// WARNING RESCALING HISTO
	      /////////////////// WARNING RESCALING HISTO
	      /////////////////// WARNING RESCALING HISTO
	      /////////////////// WARNING RESCALING HISTO
	      /////////////////// WARNING RESCALING HISTO
	      hp->correctedTOFdiffsameStation[3]->Scale(nentries/hp->correctedTOFdiffsameStation[3]->GetEntries());
	      /////////////////// WARNING RESCALING HISTO
	      /////////////////// WARNING RESCALING HISTO
	      /////////////////// WARNING RESCALING HISTO
	      /////////////////// WARNING RESCALING HISTO
	      hp->correctedTOFdiffsameStation[3]->Draw("SAME");
	      legloc->AddEntry(hp->correctedTOFdiffsameStation[3],"2 muons from same vertex","lp");
	      legloc->Draw();
	    }
	  TString s2("DiffTOF_L1TT_muonPt"); s2+=Pts[i]; s2+=".pdf";
	  c.SaveAs(s2);
	  //gStyle->SetOptStat(0);
	}
    }

  c.Clear();
  ROC[0][0]->Draw("AL");
  for (unsigned int i=0; i<Nfiles; ++i)
    for (unsigned int istud=0; istud<nMixStudies; istud++)
      ROC[i][istud]->Draw();
  
  c.Clear();
  c.Divide((Nfiles+1)/2,2);
  for (unsigned int i=0; i<Nfiles; ++i)
    {
      c.cd(i+1);
      TLegend *tl5=new TLegend(0.4,0.1,0.9,0.4,"RPC time resolution");
      ROC[i][0]->Draw("AL");
      tl5->AddEntry(ROC[i][0],ROC[i][0]->GetTitle(),"lp");
      for (unsigned int istud=1; istud<nMixStudies-1; istud++) //remove 1ns case
	{
	  ROC[i][istud]->Draw("L");
	  tl5->AddEntry(ROC[i][istud],ROC[i][istud]->GetTitle(),"lp");
	}
      tl5->Draw();
    }
  c.SaveAs("ROC_overlayTS.pdf");
   
 
  c.Clear();
  c.Divide((nMixStudies)/2,2); //remove last = 1ns case
  for (unsigned int istud=0; istud<nMixStudies-1; istud++) //remove 1ns case
    {
      c.cd(istud+1);
      TLegend *tl5=new TLegend(0.4,0.1,0.9,0.4,"Muons Pt");
      ROC[0][istud]->Draw("AL");
      tl5->AddEntry(ROC[0][istud],ROC[0][istud]->GetTitle(),"lp");
      for (unsigned int i=1; i<Nfiles; ++i)
	{
	  ROC[i][istud]->Draw("L");
	  tl5->AddEntry(ROC[i][istud],ROC[i][istud]->GetTitle(),"lp");
	}
      tl5->Draw();
    }
  c.SaveAs("ROC_overlayPt.pdf");
}



/*
{
  ZVtxStudy zv;
  HitPairStudy hp;
  zv.Loop(&hp);
  HitPairStudySmearTimeAndStrip hps(0.1);
  zv.Loop(&hp);
  hp.nhitPlus->Draw();
  hp.nhitMinus->SetLineColor(kRed);
  hp.nhitMinus->Draw("SAME");
}
*/

