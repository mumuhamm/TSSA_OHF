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
#include <vector>
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TH1.h>
#include <string>
#include <TMath.h>
#include <vector>
#include <TStyle.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TBuffer.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>


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

Float_t             _phist1;
Float_t             _phist2;
Float_t             _phist3;
Float_t             _radst1;
Float_t             _radst2;
Float_t             _radst3;
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
   virtual void     Trig_Efficiency();
   virtual void     Cuts();
   virtual void FillDiagnosticHistos(float weight);
   virtual void SetPzCut();
   virtual void BookHistos();
   virtual void SetpTArrays();
private:
   Float_t  pT =0, pdtheta =0, ref_rad =0, trk_phi =0, xF=0, chi2_trk_vtx=0;
   TVector3 vec_trk;
   TVector3 vec_vtx;
   TVector3 vec_vtx_z;
   std::string dataset;
   std::string filetype;
   std::string cuts;
   std::string triggertype;
   std::string chargetype;
   std::string armtype;
   float weight_by_thrown_pt = 0.0;
   int pTweighting;
   float KaonWeight;
   
   bool MANUAL_PHI_CUTS;
   bool MANUAL_RAD_CUTS;
   
   enum { narm = 2};
   enum { ncharge = 2};
   enum { ngap = 5};
      //enum { nptbin = 17};
      //enum { nmcptbin = 25};
   enum { nptbin = 14};
   enum { nmcptbin = 23};
   enum { nsta = 3};
   enum { npanel = 6};
   enum { norient = 2};
   enum { nhvgroup = 30};
   Int_t _arm ;
   float eff_SG3 = 1.0, eff_MUID1D = 1.0;
   float cut_eta[2];
   float cut_z[2];
   float cut_Pz[narm][2];
   float cut_dPz[narm][2];
   
   
   double varbin_pTarray[nptbin+1];
   double varbin_p_array[nptbin+1];
   double varbin_mc_pTarray[nmcptbin+1];
   
   double varbin_AN_pTarray[10];
   double varbin_AN_pzarray[5];
   //cut variable arrays, assigned by data set
   float cut_mutr_chi2[narm][ngap];
   float cut_pdtheta[narm][ngap];
   
   float cut_pdtheta__;
   float cut_mutr_chi2__;
   float cut_muid_chi2__;
   float cut_road_slope__;
   float cut_dg0__;
   float cut_ddg0__;
   float cut_vertex_chi2__;
   float cut_vertex_rad__;
   
   TFile *cut_file;
   TF1 *fcut_vertex_chi2[narm][ngap];
   TF1 *fcut_dg0[narm][ngap];
   TF1 *fcut_ddg0[narm][ngap];
   TF1 *fcut_vertex_rad[narm][ngap];
   
   TFile *fSD;
   TH1F *hSD_MUID1D[narm];
   TH1F *hSD_SG3_MUID1DH[narm];
   
   TH1F *hcut_vertex_chi2[narm][ngap];
   TH1F *hcut_dg0[narm][ngap];
   TH1F *hcut_ddg0[narm][ngap];
   TH1F *hcut_vertex_rad[narm][ngap];
   
   TFile *ftrig;
   TF1 *ftrig_SG3_pT[narm][ngap];
   TF1 *ftrig_MUID1D_pT[narm];
   
      //! Cut variable histograms
   TH2F *DG0[narm][ngap]; //DG0
   TH2F *DDG0[narm][ngap]; //DDG0
   TH2F *REFRAD[narm][ngap]; //Refrad
   TH2F *MUTR_CHI2[narm][ngap]; //MuTR chi2
   TH2F *VTX_CHI2[narm][ngap]; //vtx chi2
   TH2F *PDTHETA[narm][ngap]; //pdtheta
   
   
      //pT distribution
   TH1F *N_varbin[narm][ngap][ncharge];
   TH1F *N_varbin_AN_pT[narm][ngap][ncharge];
   TH1F *N_varbin_AN_pz[narm][ngap][ncharge];
   TH1F *N_varbin_AN_pz_MUID1D[narm][ngap][ncharge];
   TH1F *N_varbin_AN_pz_SG3MUID1DH[narm][ngap][ncharge];
   TH1F *N_varbin_fake[narm][ngap][ncharge];
   TProfile *avg_pT_varbin_gap[narm][ngap][ncharge];
   TProfile *avg_pT_varbin[narm][ncharge];
   TProfile *avg_p[narm][ngap][ncharge];
   TProfile *avg_pT_AN_pz[narm][ncharge];
   
   TH2F *N_pT_pz[narm][ngap][ngap];
   
   TH1F *N_varbin_pion[narm][ngap][ncharge];
   TH1F *N_varbin_kaon[narm][ngap][ncharge];
   TH1F *N_varbin_kaon0[narm][ngap][ncharge];
   TH1F *N_varbin_proton[narm][ngap][ncharge];
   
   TH2F *YvX_mutr_cut[narm][nsta]; //MuTr station
   TH2F *Rphi_mutr_cut[narm][nsta]; //MuTr station
  
};

#endif

#ifdef asymmetry_analyzer_cxx
asymmetry_analyzer::asymmetry_analyzer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Users/md/Documents/Phenix_HF_Analysis/analysis_preliminary_484_428757.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/Users/md/Documents/Phenix_HF_Analysis/analysis_preliminary_484_428757.root");
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
void asymmetry_analyzer::SetPzCut(){
      
   
   cut_Pz[0][0] = 3.40; //south gap2
   cut_Pz[0][1] = 3.60; //south gap3
   cut_Pz[1][0] = 3.80; //north gap2
   cut_Pz[1][1] = 4.00; //north gap3
   
   cut_dPz[0][0] = 0.20;
   cut_dPz[0][1] = 0.20;
   cut_dPz[1][0] = 0.20;
   cut_dPz[1][1] = 0.20;
}
void asymmetry_analyzer::FillDiagnosticHistos(float weight){
   
   if ( ((smlastgap==3 || smlastgap==2) && fabs(smpz)>(cut_Pz[_arm][smlastgap-2]+1.0*(pT-1.0))) || (smlastgap==4) ) {
        //DG0-DDG0
      if (
          //_rf_slope>cut_road_slope__ //GAP0SLOPE
           ref_rad<cut_vertex_rad__ //REFRAD
          && pdtheta<cut_pdtheta__ //PDTHETA
          && smtrchi2<cut_mutr_chi2__ //MUTRCHI2
          //&& chi2_trk_vtx<cut_vertex_chi2__ //VTXCHI2
          ){
         DG0[_arm][smlastgap]->Fill(pT, smdg0, weight);
         DDG0[_arm][smlastgap]->Fill(pT, smddg0, weight);
      }
      
         //REFRAD
     /* if (
          //_rf_slope>cut_road_slope__ //GAP0SLOPE
          smddg0<cut_ddg0__ //DDG0
          && smdg0<cut_dg0__ //DG0
          && pdtheta<cut_pdtheta__ //PDTHETA
          && smtrchi2<cut_mutr_chi2__ //MUTRCHI2
          && chi2_trk_vtx<cut_vertex_chi2__ //VTXCHI2
          ){
         REFRAD[_arm][smlastgap]->Fill(pT, ref_rad, weight);
      }
      
         //MUTR_CHI2 && MUID_CHI2
      if (
          //_rf_slope>cut_road_slope__ //GAP0SLOPE
           smddg0<cut_ddg0__ //DDG0
          && smdg0<cut_dg0__ //DG0
          && pdtheta<cut_pdtheta__ //PDTHETA
          && ref_rad<cut_vertex_rad__ //REFRAD
          && chi2_trk_vtx<cut_vertex_chi2__ //VTXCHI2
          ){
         MUTR_CHI2[_arm][smlastgap]->Fill(pT, smtrchi2, weight);
      }
      
         //VTXCHI2
     if (
          //_rf_slope>cut_road_slope__ //GAP0SLOPE
           smddg0<cut_ddg0__ //DDG0
          && smdg0<cut_dg0__ //DG0
          && pdtheta<cut_pdtheta__ //PDTHETA
          && smtrchi2<cut_mutr_chi2__ //MUTRCHI2
          && ref_rad<cut_vertex_rad__ //REFRAD
          ){
         VTX_CHI2[_arm][smlastgap]->Fill(pT, chi2_trk_vtx, weight);
         }
      
         //PDTHETA
      if (
          //_rf_slope>cut_road_slope__ //GAP0SLOPE
           smddg0<cut_ddg0__ //DDG0
          && smdg0<cut_dg0__ //DG0
          && smtrchi2<cut_mutr_chi2__ //MUTRCHI2
          && ref_rad<cut_vertex_rad__  //REFRAD
          && chi2_trk_vtx<2.0*cut_vertex_chi2__ //VTXCHI2
          ){
         if ( chi2_trk_vtx<cut_vertex_chi2__ ){
            PDTHETA[_arm][smlastgap]->Fill(pT, pdtheta, weight);
           
         }
      }*/
   }
}
void run_cuts::SetpTArrays(){
   
   
   
   float tmp_pTarray[nptbin+1] = {1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00};
   float tmp_p_array[nptbin+1] = {3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 7.00, 8.00, 9.00, 10.0, 12.0, 14.0, 17.0, 20.0};
   float tmp_mc_pTarray[nmcptbin+1] = {0.80, 1.00, 1.20, 1.40, 1.60,
      1.80, 2.00, 2.20, 2.40, 2.60,
      2.80, 3.00, 3.25, 3.50, 3.75,
      4.00, 4.50, 5.00, 5.50, 6.00,
      7.00, 8.00, 9.00, 10.0
   };
   
   
   for ( int ipt=0; ipt<nptbin+1; ipt++){
      varbin_pTarray[ipt] = tmp_pTarray[ipt];
      varbin_p_array[ipt] = tmp_p_array[ipt];
   }
   
   for ( int ipt=0; ipt<nmcptbin+1; ipt++){
      varbin_mc_pTarray[ipt] = tmp_mc_pTarray[ipt];
   }
   
   float tmp_AN_pTarray[10] = {1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0};
   float tmp_AN_pzarray[5] = {3.5, 5.0, 7.0, 10.0, 20.0};
   
   for (int ipt=0; ipt<10; ipt++) varbin_AN_pTarray[ipt] = tmp_AN_pTarray[ipt];
   for (int ipz=0; ipz<5; ipz++) varbin_AN_pzarray[ipz] = tmp_AN_pzarray[ipz];
   
}

#endif // #ifdef asymmetry_analyzer_cxx
