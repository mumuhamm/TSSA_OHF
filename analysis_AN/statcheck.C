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

void plot_north(TH1F * hist1 , TH1F * hist2){
   TCanvas *c = new TCanvas();
   gStyle->SetOptTitle(0);
   c->SetFillColor(0);
   c->SetBorderSize(2);
   c->SetLeftMargin(0.1422222);
   c->SetRightMargin(0.04444445);
   c->SetBottomMargin(0.00001);
   c->SetBorderMode(0);
   c->SetBottomMargin(0.15);
   c->SetBorderMode(0);
   hist1->SetLineColor(kRed);
   hist1->SetMarkerStyle(8);
   hist1->Draw();
   hist2->SetLineColor(kBlack);
   hist2->SetMarkerStyle(8);
   hist2->Draw("SAME");
   TLegend* legend = new TLegend(0.15,0.6,0.35,0.8);
   legend->SetHeader("North","L");
   legend->AddEntry(hist1,"#mu^{+}","l");
   legend->AddEntry(hist2,"#mu^{-}","l");
   legend->Draw();
  
}
void plot_south(TH1F * hist1 , TH1F * hist2){
   TCanvas *c = new TCanvas();
   c->SetFillColor(0);
   c->SetBorderSize(2);
   c->SetLeftMargin(0.1422222);
   c->SetRightMargin(0.04444445);
   c->SetBottomMargin(0.00001);
   c->SetBorderMode(0);
   c->SetBottomMargin(0.15);
   c->SetBorderMode(0);
   hist1->SetLineColor(kRed);
   hist1->SetMarkerStyle(8);
   hist1->Draw();
   hist2->SetLineColor(kBlack);
   hist2->SetMarkerStyle(8);
   hist2->Draw("SAME");
   TLegend* legend = new TLegend(0.15,0.6,0.35,0.8);
   legend->SetHeader("South","L");
   legend->AddEntry(hist1,"#mu^{+}","l");
   legend->AddEntry(hist2,"#mu^{-}","l");
   legend->Draw();
   
   
}



void statcheck(){
   
   TH1F *xF_h_n_mp = new TH1F("xF_h_n_mp", "xF_h; x_{F}(=p_{z}/p_z^{max}) (a.u.); Entries", 1000, -0.1, 0.1);
   TH1F *pdtheta_h_n_mp = new TH1F("pdtheta_h_n_mp", "pdtheta_h; p#delta#theta (rad GeV/c); Entries", 100, 0.0, 0.4);
   TH1F *ddg0_h_n_mp = new TH1F("ddg0_h_n_mp", "ddg0_h; DDG0 (cm); Entries", 100, 0, 20);
   TH1F *dg0_h_n_mp = new TH1F("dg0_h_n_mp", "dg0_h; DG0 (Degree); Entries", 100, 0, 50);
   TH1F *phi_h_n_mp = new TH1F("phi_h_n_mp", "phi_h; #phi (rad); Entries", 100, -TMath::Pi()-0.1, TMath::Pi()+0.1);
   TH1F *pt_h_n_mp = new TH1F("pt_h_n_mp", "pt_h; p_{T} (GeV/c); Entries", 100, 0, 7);
   
   TH1F *xF_h_n_mn = new TH1F("xF_h_n_mn", "xF_h; x_{F}(=p_{z}/p_z^{max}) (a.u.); Entries", 1000, -0.1, 0.1);
   TH1F *pdtheta_h_n_mn = new TH1F("pdtheta_h_n_mn", "pdtheta_h; p#delta#theta (rad GeV/c); Entries", 100, 0.0, 0.4);
   TH1F *ddg0_h_n_mn = new TH1F("ddg0_h_n_mn", "ddg0_h; DDG0 (cm); Entries", 100, 0, 20);
   TH1F *dg0_h_n_mn = new TH1F("dg0_h_n_mn", "dg0_h; DG0 (Degree); Entries", 100, 0, 50);
   TH1F *phi_h_n_mn = new TH1F("phi_h_n_mn", "phi_h; #phi (rad); Entries", 100, -TMath::Pi()-0.1, TMath::Pi()+0.1);
   TH1F *pt_h_n_mn = new TH1F("pt_h_n_mn", "pt_h; p_{T} (GeV/c); Entries", 100, 0, 7);
   
   TH1F *xF_h_s_mp = new TH1F("xF_h_s_mp", "xF_h; x_{F}(=p_{z}/p_z^{max}) (a.u.); Entries", 1000, -0.1, 0.1);
   TH1F *pdtheta_h_s_mp = new TH1F("pdtheta_h_s_mp", "pdtheta_h; p#delta#theta (rad GeV/c); Entries", 100, 0.0, 0.4);
   TH1F *ddg0_h_s_mp = new TH1F("ddg0_h_s_mp", "ddg0_h; DDG0 (cm); Entries", 100, 0, 20);
   TH1F *dg0_h_s_mp = new TH1F("dg0_h_s_mp", "dg0_h; DG0 (Degree); Entries", 100, 0, 50);
   TH1F *phi_h_s_mp = new TH1F("phi_h_s_mp", "phi_h; #phi (rad); Entries", 100, -TMath::Pi()-0.1, TMath::Pi()+0.1);
   TH1F *pt_h_s_mp = new TH1F("pt_h_s_mp", "pt_h; p_{T} (GeV/c); Entries", 100, 0, 7);
   
   TH1F *xF_h_s_mn = new TH1F("xF_h_s_mn", "xF_h; x_{F}(=p_{z}/p_z^{max}) (a.u.); Entries", 1000, -0.1, 0.1);
   TH1F *pdtheta_h_s_mn = new TH1F("pdtheta_h_s_mn", "pdtheta_h; p#delta#theta (rad GeV/c); Entries", 100, 0.0, 0.4);
   TH1F *ddg0_h_s_mn = new TH1F("ddg0_h_s_mn", "ddg0_h; DDG0 (cm); Entries", 100, 0, 20);
   TH1F *dg0_h_s_mn = new TH1F("dg0_h_s_mn", "dg0_h; DG0 (Degree); Entries", 100, 0, 50);
   TH1F *phi_h_s_mn = new TH1F("phi_h_s_mn", "phi_h; #phi (rad); Entries", 100, -TMath::Pi()-0.1, TMath::Pi()+0.1);
   TH1F *pt_h_s_mn = new TH1F("pt_h_s_mn", "pt_h; p_{T} (GeV/c); Entries", 100, 0, 7);
   
   TFile *fIn1 = new TFile("/Users/md/Documents/Phenix_HF_Analysis/fitFour.root");
      //  TFile *fIn1 = new TFile("B0sample.root");
   if (!fIn1){return;}
   TTree* smu = (TTree*)fIn1->Get("fit");
   Int_t n_entries = smu->GetEntries();
   
   std::cout<<n_entries<<"\n";
   float x_F_var, pz_var, pt_var, bluebeam_pol_var, yellowbeam_pol_var, pdtheta_var, dg0_var, ddg0_var, phi_var, rapidity_var, idchi2_var, trchi2_var ;
   int bluebeam_spin_pattern, yellowbeam_spin_pattern, run_candidate_var, run_spin_var,trhits_var, idhits_var, lastgap_var, muoncharge_var;
   bool pt_bin1;
   bool pt_bin2;
   bool pt_bin3;
   bool pt_bin4;
   bool pt_bin5;
   bool pt_bin6;
   
   bool south_cut;
   bool north_cut;
   
   
   smu->SetBranchAddress("pt_var",&pt_var);
   smu->SetBranchAddress("x_F_var",&x_F_var);
   smu->SetBranchAddress("pz_var",&pz_var);
   smu->SetBranchAddress("pt_bin1",&pt_bin1);
   smu->SetBranchAddress("pt_bin2",&pt_bin2);
   smu->SetBranchAddress("pt_bin3",&pt_bin3);
   smu->SetBranchAddress("pt_bin4",&pt_bin4);
   smu->SetBranchAddress("pt_bin5",&pt_bin5);
   smu->SetBranchAddress("pt_bin6",&pt_bin6);
   smu->SetBranchAddress("south_cut",&south_cut);
   smu->SetBranchAddress("north_cut",&north_cut);
   smu->SetBranchAddress("bluebeam_pol_var",&bluebeam_pol_var);
   smu->SetBranchAddress("yellowbeam_pol_var",&yellowbeam_pol_var);
   smu->SetBranchAddress("pdtheta_var",&pdtheta_var);
   smu->SetBranchAddress("dg0_var",&dg0_var);
   smu->SetBranchAddress("ddg0_var",&ddg0_var);
   smu->SetBranchAddress("phi_var",&phi_var);
   smu->SetBranchAddress("rapidity_var",&rapidity_var);
   smu->SetBranchAddress("idchi2_var",&idchi2_var);
   smu->SetBranchAddress("trchi2_var",&trchi2_var );
   smu->SetBranchAddress("bluebeam_spin_pattern",&bluebeam_spin_pattern);
   smu->SetBranchAddress("yellowbeam_spin_pattern",&yellowbeam_spin_pattern);
   smu->SetBranchAddress("run_candidate_var",&run_candidate_var);
   smu->SetBranchAddress("run_spin_var",&run_spin_var);
   smu->SetBranchAddress("phi_var",&phi_var);
   smu->SetBranchAddress("trhits_var",&trhits_var);
   smu->SetBranchAddress("idhits_var",&idhits_var);
   smu->SetBranchAddress("lastgap_var",&lastgap_var);
   smu->SetBranchAddress("muoncharge_var",&muoncharge_var);
   
   
   for (Int_t i=0;i<n_entries;i++){
      smu->GetEntry(i);
      
//--------------------------------------------------------------------------Fill histograms - statistics  check
      if(pz_var > 0 && north_cut && lastgap_var < 4 && (pt_var > 0 &&  pt_var < 7)){
         if(muoncharge_var ==1){
            xF_h_n_mp->Fill(x_F_var);
            pdtheta_h_n_mp->Fill(pdtheta_var);
            ddg0_h_n_mp->Fill(ddg0_var);
            dg0_h_n_mp->Fill(dg0_var);
            phi_h_n_mp->Fill(phi_var);
            pt_h_n_mp->Fill(pt_var);
         }
         if(muoncharge_var ==0){
            xF_h_n_mn->Fill(x_F_var);
            pdtheta_h_n_mn->Fill(pdtheta_var);
            ddg0_h_n_mn->Fill(ddg0_var);
            dg0_h_n_mn->Fill(dg0_var);
            phi_h_n_mn->Fill(phi_var);
            pt_h_n_mn->Fill(pt_var);
            
         }
      }
      if (pz_var < 0 && south_cut && lastgap_var < 4 && (pt_var > 0 &&  pt_var < 7)){
         if(muoncharge_var ==1){
            xF_h_s_mp->Fill(x_F_var);
            pdtheta_h_s_mp->Fill(pdtheta_var);
            ddg0_h_s_mp->Fill(ddg0_var);
            dg0_h_s_mp->Fill(dg0_var);
            phi_h_s_mp->Fill(phi_var);
            pt_h_s_mp->Fill(pt_var);
         }
         if(muoncharge_var ==0){
            xF_h_s_mn->Fill(x_F_var);
            pdtheta_h_s_mn->Fill(pdtheta_var);
            ddg0_h_s_mn->Fill(ddg0_var);
            dg0_h_s_mn->Fill(dg0_var);
            phi_h_s_mn->Fill(phi_var);
            pt_h_s_mn->Fill(pt_var);
         }
      }
      
  
      
      
      
      
      
   }

   
   std::cout<< " #========================================================Stat Table================================================================="<<"\n";
   std::cout<< " #========================================================| cuts | ================================================================="<<"\n";
   std::cout<< "These are the cuts \n"
   <<"south_cut = ( trhits > 12 && trchi2 < 10 && idhits > 6 && idchi2 < 5 && ddg0 < 8 && dg0 < 20 && (fabs(rapidity)>1.2 || fabs(rapidity)< 2.0) && pdtheta < 0.2)\n"
   <<"north_cut = ( trhits > 12 && trchi2 < 10 && idhits > 6 && idchi2 < 5 && ddg0 < 8 && dg0 < 10 && (fabs(rapidity)>1.2 || fabs(rapidity)< 2.0) && pdtheta < 0.2)\n"
   <<"For both arm : \t 0< pt< 7, r_ref < 125, lastgap < 4"<<"\n";
   std::cout<< " #========================================================Stats ================================================================="<<"\n";
   std::cout<<"|MuonPlus : xF-N:(S)|\t\t|MuonMinus : xF-N:(S)|"<<"\n";
   std::cout<<"----------------------------------------------"<<"\n";
   std::cout<<xF_h_n_mp->GetEntries()<<":"<<(xF_h_s_mp->GetEntries())<<"\t\t"<<xF_h_n_mn->GetEntries()<<":"<<(xF_h_s_mn->GetEntries())<<"\n";
   std::cout<<"|MuonPlus : pdtheta-N:(S)|\t\t|MuonMinus : pdtheta-N:(S)|"<<"\n";
   std::cout<<"----------------------------------------------"<<"\n";
   std::cout<<pdtheta_h_n_mp->GetEntries()<<":"<<(pdtheta_h_s_mp->GetEntries())<<"\t\t"<<pdtheta_h_n_mn->GetEntries()<<":"<<(pdtheta_h_s_mn->GetEntries())<<"\n";
   std::cout<<"|MuonPlus : DDG0-N:(S)|\t\t|MuonMinus : DDG0-N:(S)|"<<"\n";
   std::cout<<"----------------------------------------------"<<"\n";
   std::cout<<ddg0_h_n_mp->GetEntries()<<":"<<(ddg0_h_s_mp->GetEntries())<<"\t\t"<<ddg0_h_n_mn->GetEntries()<<":"<<(ddg0_h_s_mn->GetEntries())<<"\n";
   std::cout<<"|MuonPlus : DG0-N:(S)|\t\t|MuonMinus : DG0-N:(S)|"<<"\n";
   std::cout<<"----------------------------------------------"<<"\n";
   std::cout<<dg0_h_n_mp->GetEntries()<<":"<<(dg0_h_s_mp->GetEntries())<<"\t\t"<<dg0_h_n_mn->GetEntries()<<":"<<(dg0_h_s_mn->GetEntries())<<"\n";
   std::cout<<"|MuonPlus : PHI-N:(S)|\t\t|MuonMinus : PHI-N:(S)|"<<"\n";
   std::cout<<"----------------------------------------------"<<"\n";
   std::cout<<phi_h_n_mp->GetEntries()<<":"<<(phi_h_s_mp->GetEntries())<<"\t\t"<<phi_h_n_mn->GetEntries()<<":"<<(phi_h_s_mn->GetEntries())<<"\n";
   std::cout<<"|MuonPlus : PT-N:(S)|\t\t|MuonMinus : PT-N:(S)|"<<"\n";
   std::cout<<"----------------------------------------------"<<"\n";
   std::cout<<pt_h_n_mp->GetEntries()<<":"<<(pt_h_s_mp->GetEntries())<<"\t\t"<<pt_h_n_mn->GetEntries()<<":"<<(pt_h_s_mn->GetEntries())<<"\n";
   
   
   plot_north(xF_h_n_mp , xF_h_n_mn);
   plot_north(pdtheta_h_n_mp , pdtheta_h_n_mn);
   plot_north(ddg0_h_n_mp , ddg0_h_n_mn);
   plot_north(dg0_h_n_mp , dg0_h_n_mn);
   plot_north(phi_h_n_mp , phi_h_n_mn);
   plot_north(pt_h_n_mp , pt_h_n_mn);
   
   plot_south(xF_h_s_mp, xF_h_s_mn);
   plot_south(pdtheta_h_s_mp, pdtheta_h_s_mn);
   plot_south(ddg0_h_s_mp, ddg0_h_s_mn);
   plot_south(dg0_h_s_mp, dg0_h_s_mn);
   plot_south(phi_h_s_mp, phi_h_s_mn);
   plot_south(pt_h_s_mp, pt_h_s_mn);
  
   
}
