//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 11 15:37:22 2022 by ROOT version 5.34/36
// from TTree analysis/selected ntcltestle
// found on file: ../../muon_output/analysis_preliminary_484_428757.root
//////////////////////////////////////////////////////////

#ifndef asymmetry_analyzer_h
#define asymmetry_analyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class asymmetry_analyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nmuons;
   Float_t         smddg0;
   Float_t         smpx;
   Float_t         smpy;
   Float_t         smpz;
   Float_t         smrapidity;
   Float_t         smdg0;
   Float_t         smds3;
   Float_t         smtrchi2;
   Float_t         smidchi2;
   Float_t         smxst1;
   Float_t         smxst2;
   Float_t         smxst3;
   Float_t         smyst1;
   Float_t         smyst2;
   Float_t         smyst3;
   Float_t         smidx;
   Float_t         smidy;
   Float_t         smst1px;
   Float_t         smst1py;
   Float_t         smst1pz;
   Float_t         smdcar;
   Float_t         smdcaz;
   Float_t         smx0;
   Float_t         smx0;
   Float_t         smy0;
   Float_t         smz0;
   Int_t           smtrhits;
   Int_t           smidhits;
   Int_t           smntrhits;
   Int_t           smnidhits;
   Bool_t          smmuid1d;
   Bool_t          smmuid1s;
   Bool_t          smcharge;
   Int_t           bbcn;
   Int_t           bbcs;
   Float_t         bbcqn;
   Float_t         bbcqs;
   Float_t         bbcz;
   Float_t         bbczerr;
   Float_t         bbct0;
   Float_t         bbcts;
   Float_t         bbctn;
   Float_t         evtbbcz;
   Float_t         evtbbczerr;
   Float_t         evtvtxx;
   Float_t         evtvtxxerr;
   Float_t         evtvtxy;
   Float_t         evtvtxyerr;
   Float_t         evtvtxz;
   Float_t         evtvtxzerr;
   Float_t         smmaxres_sigma;
   Int_t           smtrackid;
   Int_t           smlastgap;
   Int_t           smhitpattern;
   Int_t           lvl1_trigraw;
   Int_t           lvl1_triglive;
   Int_t           lvl1_trigscaled;
   Int_t           lvl1_clock_cross;
   Int_t           lvl1_rbits;
   Int_t           beamclk;
   Int_t           evtsequence;
   Int_t           evtType;
   Int_t           timestamp;
   Int_t           runnumber;
   Int_t           seg_number;
   Int_t           singleM;
   Int_t           eventnumber;

   // List of branches
   TBranch        *b_nmuons;   //!
   TBranch        *b_smdd0;   //!
   TBranch        *b_smpx;   //!
   TBranch        *b_smpy;   //!
   TBranch        *b_smpz;   //!
   TBranch        *b_smrapidity;   //!
   TBranch        *b_smdg0;   //!
   TBranch        *b_smds3;   //!
   TBranch        *b_smtrchi2;   //!
   TBranch        *b_smidchi2;   //!
   TBranch        *b_smxst1;   //!
   TBranch        *b_smxst2;   //!
   TBranch        *b_smxst3;   //!
   TBranch        *b_smyst1;   //!
   TBranch        *b_smyst2;   //!
   TBranch        *b_smyst3;   //!
   TBranch        *b_smidx;   //!
   TBranch        *b_smidy;   //!
   TBranch        *b_smst1px;   //!
   TBranch        *b_smst1py;   //!
   TBranch        *b_smst1pz;   //!
   TBranch        *b_smdcar;   //!
   TBranch        *b_smdcaz;   //!
   TBranch        *b_smx0;   //!
   TBranch        *b_smx0;   //!
   TBranch        *b_smy0;   //!
   TBranch        *b_smz0;   //!
   TBranch        *b_smtrhits;   //!
   TBranch        *b_smidhits;   //!
   TBranch        *b_smntrhits;   //!
   TBranch        *b_smnidhits;   //!
   TBranch        *b_smmuid1d;   //!
   TBranch        *b_smmuid1s;   //!
   TBranch        *b_smcharge;   //!
   TBranch        *b_bbcn;   //!
   TBranch        *b_bbcs;   //!
   TBranch        *b_bbcqn;   //!
   TBranch        *b_bbcqs;   //!
   TBranch        *b_bbcz;   //!
   TBranch        *b_bbczerr;   //!
   TBranch        *b_bbct0;   //!
   TBranch        *b_bbcts;   //!
   TBranch        *b_bbctn;   //!
   TBranch        *b_evtbbcz;   //!
   TBranch        *b_evtbbczerr;   //!
   TBranch        *b_evtvtxx;   //!
   TBranch        *b_evtvtxxerr;   //!
   TBranch        *b_evtvtxy;   //!
   TBranch        *b_evtvtxyerr;   //!
   TBranch        *b_evtvtxz;   //!
   TBranch        *b_evtvtxzerr;   //!
   TBranch        *b_smmaxres_sigma;   //!
   TBranch        *b_smtrackid;   //!
   TBranch        *b_smlastgap;   //!
   TBranch        *b_smhitpattern;   //!
   TBranch        *b_lvl1_trigraw;   //!
   TBranch        *b_lvl1_triglive;   //!
   TBranch        *b_lvl1_trigscaled;   //!
   TBranch        *b_lvl1_clock_cross;   //!
   TBranch        *b_lvl1_rbits;   //!
   TBranch        *b_beamclk;   //!
   TBranch        *b_evtsequence;   //!
   TBranch        *b_evtType;   //!
   TBranch        *b_timestamp;   //!
   TBranch        *b_runnumber;   //!
   TBranch        *b_seg_number;   //!
   TBranch        *b_singleM;   //!
   TBranch        *b_eventnumber;   //!

   asymmetry_analyzer(TTree *tree=0);
   virtual ~asymmetry_analyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef asymmetry_analyzer_cxx
asymmetry_analyzer::asymmetry_analyzer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../muon_output/analysis_preliminary_484_428757.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../muon_output/analysis_preliminary_484_428757.root");
      }
      f->GetObject("analysis",tree);

   }
   Init(tree);
}

asymmetry_analyzer::~asymmetry_analyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t asymmetry_analyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t asymmetry_analyzer::LoadTree(Long64_t entry)
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

void asymmetry_analyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nmuons", &nmuons, &b_nmuons);
   fChain->SetBranchAddress("smddg0", &smddg0, &b_smdd0);
   fChain->SetBranchAddress("smpx", &smpx, &b_smpx);
   fChain->SetBranchAddress("smpy", &smpy, &b_smpy);
   fChain->SetBranchAddress("smpz", &smpz, &b_smpz);
   fChain->SetBranchAddress("smrapidity", &smrapidity, &b_smrapidity);
   fChain->SetBranchAddress("smdg0", &smdg0, &b_smdg0);
   fChain->SetBranchAddress("smds3", &smds3, &b_smds3);
   fChain->SetBranchAddress("smtrchi2", &smtrchi2, &b_smtrchi2);
   fChain->SetBranchAddress("smidchi2", &smidchi2, &b_smidchi2);
   fChain->SetBranchAddress("smxst1", &smxst1, &b_smxst1);
   fChain->SetBranchAddress("smxst2", &smxst2, &b_smxst2);
   fChain->SetBranchAddress("smxst3", &smxst3, &b_smxst3);
   fChain->SetBranchAddress("smyst1", &smyst1, &b_smyst1);
   fChain->SetBranchAddress("smyst2", &smyst2, &b_smyst2);
   fChain->SetBranchAddress("smyst3", &smyst3, &b_smyst3);
   fChain->SetBranchAddress("smidx", &smidx, &b_smidx);
   fChain->SetBranchAddress("smidy", &smidy, &b_smidy);
   fChain->SetBranchAddress("smst1px", &smst1px, &b_smst1px);
   fChain->SetBranchAddress("smst1py", &smst1py, &b_smst1py);
   fChain->SetBranchAddress("smst1pz", &smst1pz, &b_smst1pz);
   fChain->SetBranchAddress("smdcar", &smdcar, &b_smdcar);
   fChain->SetBranchAddress("smdcaz", &smdcaz, &b_smdcaz);
   fChain->SetBranchAddress("smx0", &smx0, &b_smx0);
//    fChain->SetBranchAddress("smx0", &smx0, &b_smx0);
   fChain->SetBranchAddress("smy0", &smy0, &b_smy0);
   fChain->SetBranchAddress("smz0", &smz0, &b_smz0);
   fChain->SetBranchAddress("smtrhits", &smtrhits, &b_smtrhits);
   fChain->SetBranchAddress("smidhits", &smidhits, &b_smidhits);
   fChain->SetBranchAddress("smntrhits", &smntrhits, &b_smntrhits);
   fChain->SetBranchAddress("smnidhits", &smnidhits, &b_smnidhits);
   fChain->SetBranchAddress("smmuid1d", &smmuid1d, &b_smmuid1d);
   fChain->SetBranchAddress("smmuid1s", &smmuid1s, &b_smmuid1s);
   fChain->SetBranchAddress("smcharge", &smcharge, &b_smcharge);
   fChain->SetBranchAddress("bbcn", &bbcn, &b_bbcn);
   fChain->SetBranchAddress("bbcs", &bbcs, &b_bbcs);
   fChain->SetBranchAddress("bbcqn", &bbcqn, &b_bbcqn);
   fChain->SetBranchAddress("bbcqs", &bbcqs, &b_bbcqs);
   fChain->SetBranchAddress("bbcz", &bbcz, &b_bbcz);
   fChain->SetBranchAddress("bbczerr", &bbczerr, &b_bbczerr);
   fChain->SetBranchAddress("bbct0", &bbct0, &b_bbct0);
   fChain->SetBranchAddress("bbcts", &bbcts, &b_bbcts);
   fChain->SetBranchAddress("bbctn", &bbctn, &b_bbctn);
   fChain->SetBranchAddress("evtbbcz", &evtbbcz, &b_evtbbcz);
   fChain->SetBranchAddress("evtbbczerr", &evtbbczerr, &b_evtbbczerr);
   fChain->SetBranchAddress("evtvtxx", &evtvtxx, &b_evtvtxx);
   fChain->SetBranchAddress("evtvtxxerr", &evtvtxxerr, &b_evtvtxxerr);
   fChain->SetBranchAddress("evtvtxy", &evtvtxy, &b_evtvtxy);
   fChain->SetBranchAddress("evtvtxyerr", &evtvtxyerr, &b_evtvtxyerr);
   fChain->SetBranchAddress("evtvtxz", &evtvtxz, &b_evtvtxz);
   fChain->SetBranchAddress("evtvtxzerr", &evtvtxzerr, &b_evtvtxzerr);
   fChain->SetBranchAddress("smmaxres_sigma", &smmaxres_sigma, &b_smmaxres_sigma);
   fChain->SetBranchAddress("smtrackid", &smtrackid, &b_smtrackid);
   fChain->SetBranchAddress("smlastgap", &smlastgap, &b_smlastgap);
   fChain->SetBranchAddress("smhitpattern", &smhitpattern, &b_smhitpattern);
   fChain->SetBranchAddress("lvl1_trigraw", &lvl1_trigraw, &b_lvl1_trigraw);
   fChain->SetBranchAddress("lvl1_triglive", &lvl1_triglive, &b_lvl1_triglive);
   fChain->SetBranchAddress("lvl1_trigscaled", &lvl1_trigscaled, &b_lvl1_trigscaled);
   fChain->SetBranchAddress("lvl1_clock_cross", &lvl1_clock_cross, &b_lvl1_clock_cross);
   fChain->SetBranchAddress("lvl1_rbits", &lvl1_rbits, &b_lvl1_rbits);
   fChain->SetBranchAddress("beamclk", &beamclk, &b_beamclk);
   fChain->SetBranchAddress("evtsequence", &evtsequence, &b_evtsequence);
   fChain->SetBranchAddress("evtType", &evtType, &b_evtType);
   fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
   fChain->SetBranchAddress("seg_number", &seg_number, &b_seg_number);
   fChain->SetBranchAddress("singleM", &singleM, &b_singleM);
   fChain->SetBranchAddress("eventnumber", &eventnumber, &b_eventnumber);
   Notify();
}

Bool_t asymmetry_analyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void asymmetry_analyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t asymmetry_analyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef asymmetry_analyzer_cxx
