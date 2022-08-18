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

using namespace std;

void cosinemod(){
   
   TH1F *bin1_south_spinup_xfGT0 = new TH1F("bin1_south_spinup_xfGT0", "bin1_south_spinup_xfGT0", 100, 0,7);
   TH1F *bin1_south_spindown_xfGT0 = new TH1F("bin1_south_spindown_xfGT0", "bin1_south_spindown_xfGT0", 100, 0,7);
   TH1F *bin1_north_spinup_xfGT0 = new TH1F("bin1_north_spinup_xfGT0", "bin1_north_spinup_xfGT0", 100, 0,7);
   TH1F *bin1_north_spindown_xfGT0 = new TH1F("TH1F *bin1_north_spindown_xfGT0", "TH1F *bin1_north_spindown_xfGT0", 100, 0,7);
   
   TFile *fIn1 = new TFile("../fit.root");
      //  TFile *fIn1 = new TFile("B0sample.root");
   if (!fIn1){return;}
   TTree* smu = (TTree*)fIn1->Get("fit");
   Int_t n_entries = smu->GetEntries();
   
   std::cout<<n_entries<<"\n";
   float xF, pz, pt;
   bool pt_bin1;
   bool pt_bin2;
   bool pt_bin3;
   bool pt_bin4;
   bool pt_bin5;
   bool pt_bin6;
   
   bool south_cut;
   bool north_cut;
   bool positive_mu;
   bool negative_mu;
   bool Spin_Up;
   bool Spin_Down;
   
   bool xF_positive;
   bool xF_negative;
   
   smu->SetBranchAddress("pt_var",&pt);
   smu->SetBranchAddress("x_F_var",&xF);
   smu->SetBranchAddress("pz_var",&pz);
   smu->SetBranchAddress("pt_bin1",&pt_bin1);
   smu->SetBranchAddress("pt_bin2",&pt_bin2);
   smu->SetBranchAddress("pt_bin3",&pt_bin3);
   smu->SetBranchAddress("pt_bin4",&pt_bin4);
   smu->SetBranchAddress("pt_bin5",&pt_bin5);
   smu->SetBranchAddress("pt_bin6",&pt_bin6);
   smu->SetBranchAddress("south_cut",&south_cut);
   smu->SetBranchAddress("north_cut",&north_cut);
   smu->SetBranchAddress("positive_mu",&positive_mu);
   smu->SetBranchAddress("negative_mu",&negative_mu);
   smu->SetBranchAddress("Spin_Up",&Spin_Up);
   smu->SetBranchAddress("Spin_Down",&Spin_Down);
   smu->SetBranchAddress("xF_positive",&xF_positive);
   smu->SetBranchAddress("xF_negative",&xF_negative);
   
   for (Int_t i=0;i<n_entries;i++){
      smu->GetEntry(i);
      if(pt_bin2 && south_cut && pz < 0 && Spin_Up && xF_positive){bin1_south_spinup_xfGT0->Fill(pt);}
      if(pt_bin2 && south_cut && pz < 0 && Spin_Down && xF_positive){bin1_south_spindown_xfGT0->Fill(pt);}
      if(pt_bin2 && north_cut && pz > 0 && Spin_Up && xF_positive){bin1_north_spinup_xfGT0->Fill(pt);}
      if(pt_bin2 && north_cut && pz > 0 && Spin_Down && xF_positive){bin1_north_spindown_xfGT0->Fill(pt);}
      
      
   }
      
   float asym ;
   asym = (sqrt(bin1_south_spinup_xfGT0->GetEntries()*bin1_north_spindown_xfGT0->GetEntries() ) -sqrt(bin1_south_spindown_xfGT0->GetEntries()*bin1_north_spinup_xfGT0->GetEntries()))/(sqrt(bin1_south_spinup_xfGT0->GetEntries()*bin1_north_spindown_xfGT0->GetEntries() ) +sqrt(bin1_south_spindown_xfGT0->GetEntries()*bin1_north_spinup_xfGT0->GetEntries()));
   std::cout<<asym<<"\n";
   
      TCanvas* c= new TCanvas();
   c->Divide(2,2);
   c->cd(1);bin1_south_spinup_xfGT0->Draw();
   c->cd(2);bin1_south_spindown_xfGT0->Draw();
   c->cd(3);bin1_north_spinup_xfGT0->Draw();
   c->cd(4);bin1_north_spindown_xfGT0->Draw();
   
        
      
      
      
   
}
