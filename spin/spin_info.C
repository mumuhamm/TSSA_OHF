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


int muonmultiplicity ;
int clock_candidate;
Long64_t scaler[120];
int bunch;
int run_spin[1000], run_muon[1000];

void spin_info(string output){
   
   
   
  
   TH1F * clock_h = new TH1F("clock_h", "clock_h", 100, 0, 120);
   TH1F * scalar_h = new TH1F("scalar_h", "scalar_h", 100, 0, 120);
   
   std::cout<<"======================================spin tree output======================"<<"\n";
   
   TFile *mufile_spin = new TFile("../spinDB_test.root");
   TTree *mutree_spin = (TTree*)mufile_spin->Get("T");
   int n_entries_spin = mutree_spin->GetEntries();
   std::cout<<" Number of entries for the moment : \t"<< n_entries <<"\n";
   
   int run, fill, bunch;
   float polblue, polyellow;
   int trhits, idhits, eventyield, clock;
   int y_pattern[120], b_pattern[120], cand_run, xshift;
   
   mutree_spin->SetBranchAddress("runnumber",&run);
   mutree_spin->SetBranchAddress("fillnumber",&fill);
   mutree_spin->SetBranchAddress("polblue",&polblue);
   mutree_spin->SetBranchAddress("polyellow",&polyellow);
   mutree_spin->SetBranchAddress("scalerA",&scaler[0]);
   mutree_spin->SetBranchAddress("patternblue",&b_pattern[0]);
   mutree_spin->SetBranchAddress("patternyellow",&y_pattern[0]);
   mutree_spin->SetBranchAddress("xingshift",&xshift);
   
   
   
   
   std::cout<<"======================================candidate tree output======================"<<"\n";
   
   TFile *mufile_data = new TFile("../analysis_preliminary.root");
   TTree *mutree_data = (TTree*)mufile_data->Get("analysis");
   int n_entries_data = mutree_data->GetEntries();
   std::cout<<" Number of entries for the moment : \t"<< n_entries_data <<"\n";
   
   mutree_data->SetBranchAddress("eventnumber",&eventyield);
   mutree_data->SetBranchAddress("lvl1_clock_cross",&clock_candidate);
   mutree_data->SetBranchAddress("runnumber",&cand_run);
   
   for (int ientry = 0; ientry< 2000; ++ientry){//n_entries_data
      
      mutree_data->GetEntry(ientry);
      //if (ientry%100==0) cout << "processing event " << ientry << "/" << n_entries_data <<"\n";
      //std::cout<< " run numer from the candidate tree  :"<<cand_run<< "\n";
     
      for (int jentry = 0; jentry<=n_entries; ++jentry){
         
         mutree_spin->GetEntry(jentry);
         //if (ientry%100==0) cout << "processing event " << ientry << "/" << n_entries <<"\n";
         if(cand_run != run)continue;
         int shifted_clock = (clock_candidate + xshift)%120 ;
         clock_h->Fill(shifted_clock);
         int bluepat = b_pattern[shifted_clock];
         int yellowpat = y_pattern[shifted_clock];
         
         std::cout<< " |  Run number | "<< run<< "\t" << " | Event number | " << ientry <<"\n";
         std::cout<< "| -1 : spin down | +1 spin up |"<<"\n";
         std::cout<< "------------------------------------------------------------"<<"\n";
         std::cout<<" spin pattern in blue beam | " <<bluepat <<"\n"
                  << "spin pattern in yellow beam  |   " <<yellowpat<<"\n";
         std::cout<< "------------------------------------------------------------"<<"\n";
         for(unsigned k =0 ; k<120; ++k){
            Long64_t scaler_bbcvtxcut = scaler[k];
            scalar_h->Fill(scaler_bbcvtxcut);
            
           
           
            
         }
         
         
         
      }
      
      
      
   }
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   
   TCanvas *c2 = new TCanvas();
   clock_h->SetLineColor(kRed);
   scalar_h->SetLineColor(kBlue);
   clock_h->DrawNormalized();
   scalar_h->DrawNormalized("SAME");
   TLegend* legendb = new TLegend(0.15,0.7,0.35,0.9);
   legendb->AddEntry(clock_h,"# of #mu","l");
   legendb->AddEntry(scalar_h,"scaler","l");
   legendb->Draw();
   c2->SaveAs(("/Users/md/Documents/Phenix_HF_Analysis/plot/"+output+".png").c_str());
   
}



