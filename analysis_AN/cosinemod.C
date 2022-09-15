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
   
   
   
  
   Float_t ntest;
   
   Float_t pt_diff[8]={1.375,1.75,2.25,2.75,3.25,4,6, 8.5};
   Float_t pt_diffErr[8] ={0};
   Float_t N;
   Float_t sinPhi;
   Float_t cosPhi;
   Float_t epsilon[8] = {0};
   Float_t epsilon2[8] = {0};
   Float_t epsilon_PM_NO[8] = {0};
   Float_t epsilon_PM_NO_Error[8] = {0};
   Float_t epsilon_NM_NO[8] = {0};
   Float_t epsilon_NM_NO_Error[8] = {0};
   Float_t epsilon_Error[8] = {0};
   Float_t epsilon2_Error[8] = {0};
   Float_t epsilon_PHI[12] = {0};
   Float_t epsilon_PHIError[12] = {0};
   TH1F * cosmod_h[6];
   for(Int_t i=0; i<6; ++i){ cosmod_h[i]= new TH1F(Form("myhist%d",i), Form("myhist%d",i), 100, -0.2, 0.2);}
   
   Float_t PHI_SA_YB_PM_Event_R_UP[8][12] = {0};
   Float_t PHI_SA_YB_PM_Event_L_UP[8][12]= {0};
   Float_t PHI_SA_YB_PM_Event_R_DOWN[8][12]= {0};
   Float_t PHI_SA_YB_PM_Event_L_DOWN[8][12]= {0};
   
   Float_t PHI_SA_YB_NM_Event_R_UP[8][12] = {0};
   Float_t PHI_SA_YB_NM_Event_L_UP[8][12] = {0};
   Float_t PHI_SA_YB_NM_Event_R_DOWN[8][12] = {0};
   Float_t PHI_SA_YB_NM_Event_L_DOWN[8][12] = {0};
   
   
   
   
   Float_t SA_YB_PM_Event_R_UP[8] = {0};
   Float_t SA_YB_PM_Event_L_UP[8]= {0};
   Float_t SA_YB_PM_Event_R_DOWN[8]= {0};
   Float_t SA_YB_PM_Event_L_DOWN[8]= {0};
   
   Float_t SA_YB_NM_Event_R_UP[8] = {0};
   Float_t SA_YB_NM_Event_L_UP[8] = {0};
   Float_t SA_YB_NM_Event_R_DOWN[8]= {0};
   Float_t SA_YB_NM_Event_L_DOWN[8]= {0};
   
   Float_t NO_BB_PM_Event_R_UP[8] = {0};
   Float_t NO_BB_PM_Event_L_UP[8]= {0};
   Float_t NO_BB_PM_Event_R_DOWN[8]= {0};
   Float_t NO_BB_PM_Event_L_DOWN[8]= {0};
   
   Float_t NO_BB_NM_Event_R_UP[8] = {0};
   Float_t NO_BB_NM_Event_L_UP[8] = {0};
   Float_t NO_BB_NM_Event_R_DOWN[8]= {0};
   Float_t NO_BB_NM_Event_L_DOWN[8]= {0};
   
   Float_t COSPHI_PM[8]={0};
   Float_t COSPHI_NM[8]={0};
   Float_t COSPHI_PHI[12]={0};
   Float_t PHIBIN[13] = {-TMath::Pi(), -TMath::Pi()*5/6, -TMath::Pi()*2/3,-TMath::Pi()/2, -TMath::Pi()/3, -TMath::Pi()/6, 0, TMath::Pi()/6,TMath::Pi()/3, TMath::Pi()/2, TMath::Pi()*2/3,TMath::Pi()*5/6, TMath::Pi()};
   Float_t PTBIN[9] = {1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 5, 7, 10};
   Int_t  binnumpt = sizeof(PTBIN)/sizeof(Float_t) - 1;
   TH1F *PT_H = new TH1F("PT_H", "PT_H", 8, PTBIN);
   TH1F *AN_PT_H = new TH1F("AN_PT_H", "AN_PT_H", 8, PTBIN);
   TH1F * PHI_H = new TH1F("PHI_H", "PHI_H", 12, PHIBIN);
   TFile *fIn1 = new TFile("/Users/md/Documents/Phenix_HF_Analysis/fitFour.root");
   if (!fIn1){return;}
   TTree* smu = (TTree*)fIn1->Get("fit");
   Int_t n_entries = smu->GetEntries();
   
   std::cout<<n_entries<<"\n";
   float x_F_var, pz_var, pt_var, bluebeam_pol_var, yellowbeam_pol_var, pdtheta_var, dg0_var, ddg0_var, phi_var, rapidity_var, idchi2_var, trchi2_var ;
   int bluebeam_spin_pattern, yellowbeam_spin_pattern, run_candidate_var, run_spin_var,trhits_var, idhits_var, lastgap_var, muoncharge_var;
   float px_var, r_ref_var;
   bool pt_bin1;
   bool pt_bin2;
   bool pt_bin3;
   bool pt_bin4;
   bool pt_bin5;
   bool pt_bin6;
   
   bool south_cut;
   bool north_cut;
   Float_t polarization;
   
   smu->SetBranchAddress("pt_var",&pt_var);
   smu->SetBranchAddress("x_F_var",&x_F_var);
   smu->SetBranchAddress("pz_var",&pz_var);
   smu->SetBranchAddress("px_var",&px_var);
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
   smu->SetBranchAddress("r_ref_var",&r_ref_var);
   
   
   for (Int_t i=0;i<n_entries;i++){
      smu->GetEntry(i);
   if(run_candidate_var != run_spin_var)continue;
      polarization = yellowbeam_pol_var;
      sinPhi = sin(phi_var);
      cosPhi = cos(phi_var);
      PHI_H->Fill(phi_var);
      PT_H->Fill(pt_var);
      //std::cout<<PHI_H->GetBinContent(1)<<"\n";
      
  //-------------------------------------------- Studies in x_F < 0 :South : Yellow , North : Blue -------------Run by run study & integration : Calculation on Phi
     // r_ref_var < 125
      if (x_F_var < 0){
         
         //South & Yellow beam
         
         if(south_cut &&  r_ref_var < 180  && px_var > 0 && muoncharge_var ==1){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if (yellowbeam_spin_pattern == 1)SA_YB_PM_Event_R_UP[k] = PT_H->GetBinContent(k);
               if (yellowbeam_spin_pattern == -1)SA_YB_PM_Event_R_DOWN[k] = PT_H->GetBinContent(k);
               for(Int_t j=0;j<PHI_H->GetNbinsX(); ++j){
                  if (yellowbeam_spin_pattern == 1)PHI_SA_YB_PM_Event_R_UP[k][j] = PHI_H->GetBinContent(j);
                  if (yellowbeam_spin_pattern == -1)PHI_SA_YB_PM_Event_R_DOWN[k][j] = PHI_H->GetBinContent(j);
                  
                  std::cout<<" lets get printout the bin wise events in case phi binning : k & j "<<k<<"\t:&:\t"<<j<<"\t"<<PHI_SA_YB_PM_Event_R_UP[k][j]<<"\n";
               }
               
            }
         }
         
         
         if(south_cut && r_ref_var < 180  && px_var < 0 && muoncharge_var ==1){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if( yellowbeam_spin_pattern == 1)SA_YB_PM_Event_L_UP[k] = PT_H->GetBinContent(k);
               if( yellowbeam_spin_pattern == -1)SA_YB_PM_Event_L_DOWN[k] = PT_H->GetBinContent(k);
               for(Int_t j=0;j<PHI_H->GetNbinsX(); ++j){
                  if( yellowbeam_spin_pattern == 1)PHI_SA_YB_PM_Event_L_UP[k][j] = PHI_H->GetBinContent(j);
                  if( yellowbeam_spin_pattern == -1)PHI_SA_YB_PM_Event_L_DOWN[k][j] = PHI_H->GetBinContent(j);
                 
               }
            }
         }
         
         
         
         
         
         if(south_cut &&  r_ref_var < 180  && px_var > 0 && muoncharge_var ==0){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if(yellowbeam_spin_pattern == 1)SA_YB_NM_Event_R_UP[k] = PT_H->GetBinContent(k);
               if(yellowbeam_spin_pattern == -1)SA_YB_NM_Event_R_DOWN[k] = PT_H->GetBinContent(k);
            }
         }
         
        
         if(south_cut && r_ref_var < 180  && px_var < 0 && muoncharge_var ==0){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if(yellowbeam_spin_pattern == 1)SA_YB_NM_Event_L_UP[k] = PT_H->GetBinContent(k);
               if(yellowbeam_spin_pattern == -1)SA_YB_NM_Event_L_DOWN[k] = PT_H->GetBinContent(k);
            }
         }
        
         
        
         //North Blue beam
         
         if(north_cut && r_ref_var < 180  && px_var > 0 && muoncharge_var ==1 ){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if(bluebeam_spin_pattern == 1)NO_BB_PM_Event_R_UP[k] = PT_H->GetBinContent(k);
               if(bluebeam_spin_pattern == -1)NO_BB_PM_Event_R_DOWN[k] = PT_H->GetBinContent(k);
            }
         }
       
        
         if(north_cut && r_ref_var < 180  && px_var < 0 && muoncharge_var ==1 ){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if(bluebeam_spin_pattern == 1)NO_BB_PM_Event_L_UP[k] = PT_H->GetBinContent(k);
               if(bluebeam_spin_pattern == -1)NO_BB_PM_Event_L_DOWN[k] = PT_H->GetBinContent(k);
            }
         }
         
         
         if(north_cut && r_ref_var < 180 && px_var > 0 && muoncharge_var ==0 ){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if( bluebeam_spin_pattern == 1)NO_BB_NM_Event_R_UP[k] = PT_H->GetBinContent(k);
               if( bluebeam_spin_pattern == -1)NO_BB_NM_Event_R_DOWN[k] = PT_H->GetBinContent(k);
            }
         }
         
         
         if(north_cut && r_ref_var < 180  && px_var < 0 && muoncharge_var ==0 ){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if( bluebeam_spin_pattern == 1)NO_BB_NM_Event_L_UP[k] = PT_H->GetBinContent(k);
               if( bluebeam_spin_pattern == -1)NO_BB_NM_Event_L_DOWN[k] = PT_H->GetBinContent(k);
            }
         }
        
         
      }
      
      
   
      
      
      
      
   }//event loop closed
   
   
   
    std::cout<<"-----------------------------------------------------------------------------------------------------------------------------------------------"<<"\n";
   
   for(Int_t j=0; j<8; ++j){
      
      Float_t var_1 = sqrt(SA_YB_PM_Event_L_UP[j] * SA_YB_PM_Event_R_DOWN[j]);
      Float_t var_2 = sqrt(SA_YB_PM_Event_L_DOWN[j] * SA_YB_PM_Event_R_UP[j]);
      epsilon[j] = ((var_1 - var_2)/(var_1+var_2))/(polarization*cosPhi);
      epsilon_Error[j] = (sqrt(SA_YB_PM_Event_L_UP[j]*SA_YB_PM_Event_R_DOWN[j]*SA_YB_PM_Event_L_DOWN[j] * SA_YB_PM_Event_R_UP[j])/((var_1+var_2)*(var_1+var_2)))*sqrt(1/SA_YB_PM_Event_L_UP[j] + 1/SA_YB_PM_Event_R_DOWN[j] + 1/SA_YB_PM_Event_L_DOWN[j] + 1/SA_YB_PM_Event_R_UP[j] );
      
      
      
      Float_t var_3 = sqrt(SA_YB_NM_Event_L_UP[j] * SA_YB_NM_Event_R_DOWN[j]);
      Float_t var_4 = sqrt(SA_YB_NM_Event_L_DOWN[j] * SA_YB_NM_Event_R_UP[j]);
      epsilon2[j] = ((var_3 - var_4)/(var_3+var_4))/(polarization*cosPhi);
      epsilon2_Error[j] = (sqrt(SA_YB_NM_Event_L_UP[j]*SA_YB_NM_Event_R_DOWN[j]*SA_YB_NM_Event_L_DOWN[j] * SA_YB_NM_Event_R_UP[j])/((var_1+var_2)*(var_1+var_2)))*sqrt(1/SA_YB_NM_Event_L_UP[j] + 1/SA_YB_NM_Event_R_DOWN[j] + 1/SA_YB_NM_Event_L_DOWN[j] + 1/SA_YB_NM_Event_R_UP[j] );
      
      
      Float_t var_5 = sqrt(NO_BB_PM_Event_L_UP[j]*NO_BB_PM_Event_R_DOWN[j]);
      Float_t var_6 = sqrt(NO_BB_PM_Event_L_DOWN[j]*NO_BB_PM_Event_R_UP[j]);
      epsilon_PM_NO = ((var_5 - var_6)/(var_5+var_6))/(polarization*cosPhi);
      epsilon_PM_NO_Error = (sqrt(NO_BB_PM_Event_L_UP[j]*NO_BB_PM_Event_R_UP[j]*NO_BB_PM_Event_R_UP[j]*NO_BB_PM_Event_R_DOWN[j])/((var_5+var_6)*(var_5+var_6)))*sqrt(1/NO_BB_PM_Event_L_UP[j] + 1/NO_BB_PM_Event_R_DOWN[j] + 1/NO_BB_PM_Event_L_DOWN[j] + 1/ NO_BB_PM_Event_R_DOWN[j]);
      
      
      Float_t var_7 = sqrt(NO_BB_NM_Event_L_UP[j]*NO_BB_NM_Event_R_DOWN[j]);
      Float_t var_8 = sqrt(NO_BB_NM_Event_L_DOWN[j]*NO_BB_NM_Event_R_UP[j]);
      epsilon_NM_NO = ((var_7 - var_8)/(var_7+var_8))/(polarization*cosPhi);
      epsilon_NM_NO_Error = (sqrt(NO_BB_PM_Event_L_UP[j]*NO_BB_PM_Event_R_UP[j]*NO_BB_PM_Event_R_UP[j]*NO_BB_PM_Event_R_DOWN[j])/((var_5+var_6)*(var_5+var_6)))*sqrt(1/NO_BB_PM_Event_L_UP[j] + 1/NO_BB_PM_Event_R_DOWN[j] + 1/NO_BB_PM_Event_L_DOWN[j] + 1/ NO_BB_PM_Event_R_DOWN[j]);
   }
   
   for(Int_t i=0; i<8; ++i){
      for(Int_t j=0; j<12; ++j){
         
         Float_t var_1 = sqrt(PHI_SA_YB_PM_Event_R_UP[i][j] * PHI_SA_YB_PM_Event_R_DOWN[i][j]);
         Float_t var_2 = sqrt(PHI_SA_YB_PM_Event_L_DOWN[i][j] * PHI_SA_YB_PM_Event_R_UP[i][j]);
         epsilon_PHI[j] = ((var_1 - var_2)/(var_1+var_2))/(polarization*cosPhi);
         std::cout<<"lets take the print out of this fellow raw asymmetry in phi case when jump to next pT bin "<<i<<": \t"<< epsilon_PHI[j] <<"\n";
         
      }
   }
   
   
  
   
  
    
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   TF1 *f1 = new TF1("f1","[0]*cos(x)",1.25,10);
   
   auto gr = new TGraphErrors(binnumpt,pt_diff,epsilon,pt_diffErr,epsilon_Error);
   gr->SetTitle("x_F < 0 :: Y : AN, X :PT");
   gr->SetMarkerColor(4);
   gr->SetLineColor(4);
   gr->SetMarkerStyle(8);
   gr->GetXaxis()->SetTitle("p_{T}");
   gr->GetYaxis()->SetTitle("A_{N}");
   //gr->Fit("f1","q");
   //gr->Draw("AP");
   //gr->Draw("a p s ; ; 5 s=0.5");
   
   
   auto gr1 = new TGraphErrors(binnumpt, pt_diff, epsilon2 , pt_diffErr,epsilon2_Error);
   gr1->SetTitle("x_F < 0 :: Y : AN, X :PT");
   gr1->SetMarkerColor(1);
   gr1->SetLineColor(1);
   gr1->SetMarkerStyle(8);
   gr1->GetXaxis()->SetTitle("p_{T}");
   gr1->GetYaxis()->SetTitle("A_{N}");
      //gr1->Fit("f1","q");
      //gr1->Draw("a p s ; ; 5 s=0.5");
   
   auto c1 = new TCanvas("c1","Y : AN, X :PT",200,10,700,500);
   c1->SetFillColor(0);
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gr);
   mg->Add(gr1);
      //mg->GetXaxis()->SetTitle("p_{T}");
      //mg->GetYaxis()->SetTitle("A_{N}");
   mg->Draw("a p s ; ; 5 s=0.5");
   TLatex T0;
   T0.SetTextColor(2);
   T0.DrawLatexNDC(.4,.30, "x_{F}<0, Yellow, A_{N}=#frac{#epsilon}{P#bulletcos#phi}");
   TLatex T1;
   T1.SetTextColor(4);
   T1.DrawLatexNDC(.6,.20, "#mu^{+}");
   TLatex T2;
   T2.SetTextColor(1);
   T2.DrawLatexNDC(.7,.20, "#mu^{-}");

   
   /*float asym ;
   asym = (sqrt(bin1_south_spinup_xfGT0->GetEntries()*bin1_north_spindown_xfGT0->GetEntries() ) -sqrt(bin1_south_spindown_xfGT0->GetEntries()*bin1_north_spinup_xfGT0->GetEntries()))/(sqrt(bin1_south_spinup_xfGT0->GetEntries()*bin1_north_spindown_xfGT0->GetEntries() ) +sqrt(bin1_south_spindown_xfGT0->GetEntries()*bin1_north_spinup_xfGT0->GetEntries()));
   std::cout<<asym<<"\n";
    
    auto c2 = new TCanvas("c2","Y : AN, X :PT",200,10,700,500);
    c2->SetFillColor(0);
    c2->GetFrame()->SetFillColor(21);
    c2->GetFrame()->SetBorderSize(12);
  */
   
   
   
   
  
   
   
  
}
