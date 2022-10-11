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
#include <iomanip>
#include <TStyle.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TBuffer.h>

using namespace std;

Double_t cfunc(Double_t *x, Double_t *par){
  Double_t func;
  Double_t xx = x[0];
  func = par[0]*cos(xx);
  return func;
}


void phi_test(){
   Double_t PHI_SA_YB_PM_Event_R_UP[6][12] = {0};
   Double_t PHI_SA_YB_PM_Event_L_UP[6][12]= {0};
   Double_t PHI_SA_YB_PM_Event_R_DOWN[6][12]= {0};
   Double_t PHI_SA_YB_PM_Event_L_DOWN[6][12]= {0};
   
   
   Double_t epsilonPHI_SA_YB_PM[12] = {0};
   Double_t epsilonPHI_SA_YB_PMError[12] = {0};
   
   Double_t phi_diff[12] ={-(TMath::Pi()-TMath::Pi()/12),-(TMath::Pi()*5/6-TMath::Pi()/12), -(TMath::Pi()*2/3-TMath::Pi()/12),-(TMath::Pi()/2-TMath::Pi()/12), -(TMath::Pi()/3-TMath::Pi()/12), -(TMath::Pi()/6-TMath::Pi()/12),  (0.0+TMath::Pi()/12), (TMath::Pi()/6+TMath::Pi()/12), (TMath::Pi()/3+TMath::Pi()/12), (TMath::Pi()/2+TMath::Pi()/12), (TMath::Pi()*2/3+TMath::Pi()/12), (TMath::Pi()*5/6+TMath::Pi()/12)} ;
   Double_t phi_diffErr[12] ={0};
   
   Double_t PHIBIN[13] = {-TMath::Pi(), -TMath::Pi()*5/6, -TMath::Pi()*2/3,-TMath::Pi()/2, -TMath::Pi()/3, -TMath::Pi()/6, 0, TMath::Pi()/6,TMath::Pi()/3, TMath::Pi()/2, TMath::Pi()*2/3,TMath::Pi()*5/6, TMath::Pi()};
   Int_t  binnumphi = sizeof(PHIBIN)/sizeof(Double_t) - 1;
   std::cout <<" bin number in phi distribution  : "<<binnumphi <<"\n";
   Double_t y=7, hj = 5;
   std::cout<<" sqrt test : squre root of 30 : "<< sqrt(y*hj)<<"\n";
   TH1D * PHI_H = new TH1D("PHI_H", "PHI_H", binnumphi, PHIBIN);
   

TFile *fIn1 = new TFile("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/alibordi/HF_Analysis/samples/fit_tree_runXV.root");

if (!fIn1){return;}
   TTree* smu = (TTree*)fIn1->Get("fit");
   Int_t n_entries = smu->GetEntries();
std::cout<<n_entries<<"\n";
float x_F_var, pz_var, pt_var, bluebeam_pol_var, yellowbeam_pol_var, pdtheta_var, dg0_var, ddg0_var, phi_var, rapidity_var, idchi2_var, trchi2_var ;
   int bluebeam_spin_pattern, yellowbeam_spin_pattern, run_candidate_var, run_spin_var,trhits_var, idhits_var, lastgap_var, muoncharge_var;
   float px_var, r_ref_var;
   Float_t polarization_yellow;
   bool pt_bin1;
   bool pt_bin2;
   bool pt_bin3;
   bool pt_bin4;
   bool pt_bin5;
   bool pt_bin6;
   smu->SetBranchAddress("phi_var",&phi_var);
   smu->SetBranchAddress("yellowbeam_spin_pattern",&yellowbeam_spin_pattern);
   smu->SetBranchAddress("px_var",&px_var);
   smu->SetBranchAddress("pt_bin1",&pt_bin1);
   smu->SetBranchAddress("pt_bin2",&pt_bin2);
   smu->SetBranchAddress("pt_bin3",&pt_bin3);
   smu->SetBranchAddress("pt_bin4",&pt_bin4);
   smu->SetBranchAddress("pt_bin5",&pt_bin5);
   smu->SetBranchAddress("pt_bin6",&pt_bin6);
   smu->SetBranchAddress("muoncharge_var",&muoncharge_var);
   smu->SetBranchAddress("yellowbeam_pol_var",&yellowbeam_pol_var);

for (Int_t i=0;i<n_entries;i++){
      smu->GetEntry(i);
      PHI_H->Fill(phi_var);
      polarization_yellow = yellowbeam_pol_var;
if (muoncharge_var == 1 ){
         for(Int_t j=1;j<=PHI_H->GetNbinsX(); ++j){
            if( yellowbeam_spin_pattern == 1){
               if (px_var < 0){
                  if(pt_bin1)PHI_SA_YB_PM_Event_L_UP[0][j] = PHI_H->GetBinContent(j);
                  if(pt_bin2)PHI_SA_YB_PM_Event_L_UP[1][j] = PHI_H->GetBinContent(j);
                  if(pt_bin3)PHI_SA_YB_PM_Event_L_UP[2][j] = PHI_H->GetBinContent(j);
                  if(pt_bin4)PHI_SA_YB_PM_Event_L_UP[3][j] = PHI_H->GetBinContent(j);
                  if(pt_bin5)PHI_SA_YB_PM_Event_L_UP[4][j] = PHI_H->GetBinContent(j);
                  if(pt_bin6)PHI_SA_YB_PM_Event_L_UP[5][j] = PHI_H->GetBinContent(j);
               }
               if (px_var > 0){
                  if(pt_bin1)PHI_SA_YB_PM_Event_R_UP[0][j] = PHI_H->GetBinContent(j);
                  if(pt_bin2)PHI_SA_YB_PM_Event_R_UP[1][j] = PHI_H->GetBinContent(j);
                  if(pt_bin3)PHI_SA_YB_PM_Event_R_UP[2][j] = PHI_H->GetBinContent(j);
                  if(pt_bin4)PHI_SA_YB_PM_Event_R_UP[3][j] = PHI_H->GetBinContent(j);
                  if(pt_bin5)PHI_SA_YB_PM_Event_R_UP[4][j] = PHI_H->GetBinContent(j);
                  if(pt_bin6)PHI_SA_YB_PM_Event_R_UP[5][j] = PHI_H->GetBinContent(j);
                  
               }
               
            }
            if( yellowbeam_spin_pattern == -1){
               if(px_var < 0){
                  if(pt_bin1)PHI_SA_YB_PM_Event_L_DOWN[0][j] = PHI_H->GetBinContent(j);
                  if(pt_bin2)PHI_SA_YB_PM_Event_L_DOWN[1][j] = PHI_H->GetBinContent(j);
                  if(pt_bin3)PHI_SA_YB_PM_Event_L_DOWN[2][j] = PHI_H->GetBinContent(j);
                  if(pt_bin4)PHI_SA_YB_PM_Event_L_DOWN[3][j] = PHI_H->GetBinContent(j);
                  if(pt_bin5)PHI_SA_YB_PM_Event_L_DOWN[4][j] = PHI_H->GetBinContent(j);
                  if(pt_bin6)PHI_SA_YB_PM_Event_L_DOWN[5][j] = PHI_H->GetBinContent(j);
               }
               if (px_var > 0){
                  if(pt_bin1)PHI_SA_YB_PM_Event_R_DOWN[0][j] = PHI_H->GetBinContent(j);
                  if(pt_bin2)PHI_SA_YB_PM_Event_R_DOWN[1][j] = PHI_H->GetBinContent(j);
                  if(pt_bin3)PHI_SA_YB_PM_Event_R_DOWN[2][j] = PHI_H->GetBinContent(j);
                  if(pt_bin4)PHI_SA_YB_PM_Event_R_DOWN[3][j] = PHI_H->GetBinContent(j);
                  if(pt_bin5)PHI_SA_YB_PM_Event_R_DOWN[4][j] = PHI_H->GetBinContent(j);
                  if(pt_bin6)PHI_SA_YB_PM_Event_R_DOWN[5][j] = PHI_H->GetBinContent(j);
               }
            }


 }
         
      }//muon charge
      
}// event loop
  

TGraphErrors *gr_phi_sa_yb_pm[6];
   auto cm = new TCanvas("cm","Yellow : AN, X :phi",200,10,700,500);
   
   for(Int_t i=0; i<6; ++i){
      for(Int_t j=1; j<=12; ++j){
         
          std::cout<<"pt bin : "<< i<< "\t"<< " phi bin : "<< j<<"\t"<<" L Up : "<<PHI_SA_YB_PM_Event_L_UP[i][j] <<"\n";
          std::cout<<"pt bin : "<< i<< "\t"<< " phi bin : "<< j<<"\t"<<" R Down : "<<PHI_SA_YB_PM_Event_R_DOWN[i][j] <<"\n";
          std::cout<<"pt bin : "<< i<< "\t"<< " phi bin : "<< j<<"\t"<<" L Down : "<<PHI_SA_YB_PM_Event_L_DOWN[i][j] <<"\n";
          std::cout<<"pt bin : "<< i<< "\t"<< " phi bin : "<< j<<"\t"<<" R Up : "<<PHI_SA_YB_PM_Event_R_UP[i][j] <<"\n";
         Float_t var_1 = sqrt(PHI_SA_YB_PM_Event_L_UP[i][j] * PHI_SA_YB_PM_Event_R_DOWN[i][j]);
         std::cout<< " var1 : "<<var_1<<"\n";
         Float_t var_2 = sqrt(PHI_SA_YB_PM_Event_L_DOWN[i][j] * PHI_SA_YB_PM_Event_R_UP[i][j]);
         std::cout<< " var2 : "<<var_2<<"\n";
         epsilonPHI_SA_YB_PM[j] = ((var_1 - var_2)/((var_1 + var_2)*polarization_yellow));
         epsilonPHI_SA_YB_PMError[j] = (sqrt(PHI_SA_YB_PM_Event_L_UP[i][j]*PHI_SA_YB_PM_Event_R_DOWN[i][j]*PHI_SA_YB_PM_Event_L_DOWN[i][j]* PHI_SA_YB_PM_Event_R_UP[i][j])/((var_1+var_2)*(var_1+var_2)))*sqrt(1/PHI_SA_YB_PM_Event_L_UP[i][j] + 1/PHI_SA_YB_PM_Event_R_DOWN[i][j] + 1/PHI_SA_YB_PM_Event_L_DOWN[i][j] + 1/ PHI_SA_YB_PM_Event_R_UP[i][j]);
         
         std::cout<<" raw asymmetry in phi case when jump to next pT bin "<<i<<": \t"<< epsilonPHI_SA_YB_PM[j]<<"\t"<<epsilonPHI_SA_YB_PMError[j]<<"\n";
         
      }
auto fitform  = new TF1("fitform",cfunc,-3.14,3.14,1);
      fitform->SetParameter(0,0.05);
      fitform->SetParLimits(0,-1.0,1.0);
      fitform->SetParName(0,"p0");
      
      gr_phi_sa_yb_pm[i] = new TGraphErrors(binnumphi,phi_diff,epsilonPHI_SA_YB_PM,phi_diffErr,epsilonPHI_SA_YB_PMError);
      gr_phi_sa_yb_pm[i]->SetTitle("x_F < 0 :: Y : AN, X :PHI");
      gr_phi_sa_yb_pm[i]->SetMarkerColor(1);
      gr_phi_sa_yb_pm[i]->SetLineColor(1);
      gr_phi_sa_yb_pm[i]->SetMarkerStyle(8);
      gr_phi_sa_yb_pm[i]->GetXaxis()->SetTitle("#phi");
      gr_phi_sa_yb_pm[i]->GetYaxis()->SetTitle("A_{N}");
      gr_phi_sa_yb_pm[i]->Fit(fitform,"M");
       TFitResultPtr cosinere =    gr_phi_sa_yb_pm[i]->Fit(fitform, "S");
       TMatrixDSym cov_cosine = cosinere->GetCovarianceMatrix();
       Double_t chi2_cosine   = cosinere->Chi2();
       Double_t par0_cosine   = cosinere->Parameter(0);
       Double_t err0_cosine   = cosinere->ParError(0);
       cosinere->Print("V");
      gr_phi_sa_yb_pm[i]->Draw("a p s ; ; 5 s=0.5");
      cm->SaveAs(Form("../../plots/cosine_modulous_%d.pdf",i));
   }
   
   


}
