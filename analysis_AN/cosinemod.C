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
   Float_t epsilon_SA_BB_PM[8] = {0};
   Float_t epsilon_SA_BB_PMError[8] = {0};
   Float_t epsilon_SA_BB_NM[8] = {0};
   Float_t epsilon_SA_BB_NMError[8] = {0};
   Float_t epsilon_NO_YB_PM[8] = {0};
   Float_t epsilon_NO_YB_PMError[8] = {0};
   Float_t epsilon_NO_YB_NM[8] = {0};
   Float_t epsilon_NO_YB_NMError[8] = {0};
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
   
   Float_t SA_BB_PM_Event_R_UP[8] = {0};
   Float_t SA_BB_PM_Event_L_UP[8]= {0};
   Float_t SA_BB_PM_Event_R_DOWN[8]= {0};
   Float_t SA_BB_PM_Event_L_DOWN[8]= {0};
   
   Float_t SA_BB_NM_Event_R_UP[8] = {0};
   Float_t SA_BB_NM_Event_L_UP[8] = {0};
   Float_t SA_BB_NM_Event_R_DOWN[8]= {0};
   Float_t SA_BB_NM_Event_L_DOWN[8]= {0};
   
   Float_t NO_BB_PM_Event_R_UP[8] = {0};
   Float_t NO_BB_PM_Event_L_UP[8]= {0};
   Float_t NO_BB_PM_Event_R_DOWN[8]= {0};
   Float_t NO_BB_PM_Event_L_DOWN[8]= {0};
   
   Float_t NO_BB_NM_Event_R_UP[8] = {0};
   Float_t NO_BB_NM_Event_L_UP[8] = {0};
   Float_t NO_BB_NM_Event_R_DOWN[8]= {0};
   Float_t NO_BB_NM_Event_L_DOWN[8]= {0};
   
   Float_t NO_YB_PM_Event_R_UP[8] = {0};
   Float_t NO_YB_PM_Event_L_UP[8]= {0};
   Float_t NO_YB_PM_Event_R_DOWN[8]= {0};
   Float_t NO_YB_PM_Event_L_DOWN[8]= {0};
   
   Float_t NO_YB_NM_Event_R_UP[8] = {0};
   Float_t NO_YB_NM_Event_L_UP[8] = {0};
   Float_t NO_YB_NM_Event_R_DOWN[8]= {0};
   Float_t NO_YB_NM_Event_L_DOWN[8]= {0};
   
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
   Float_t polarization_yellow;
   Float_t polarization_blue;
   
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
      polarization_yellow = yellowbeam_pol_var;
      polarization_blue = bluebeam_pol_var;
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
                  
                 // std::cout<<" lets get printout the bin wise events in case phi binning : k & j "<<k<<"\t:&:\t"<<j<<"\t"<<PHI_SA_YB_PM_Event_R_UP[k][j]<<"\n";
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
        
         
      }//x_F < 0 closes
      
      //-------------------------------------------- Studies in x_F > 0 :South : Blue , North : Yellow -------------Run by run study & integration : Calculation on Phi
      
      if(x_F_var > 0 ){
         // South and blue beam
         if(south_cut && r_ref_var < 180 && muoncharge_var ==1 ){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if(bluebeam_spin_pattern == 1 ){
                  
                  if (px_var < 0) SA_BB_PM_Event_L_UP[k]=PT_H->GetBinContent(k);
                  if (px_var > 0) SA_BB_PM_Event_R_UP[k]=PT_H->GetBinContent(k);
               }
               if(bluebeam_spin_pattern == -1 ){

                  if (px_var < 0) SA_BB_PM_Event_L_DOWN[k] = PT_H->GetBinContent(k);
                  if (px_var > 0) SA_BB_PM_Event_R_DOWN[k] = PT_H->GetBinContent(k);
                  
               }
            }
         }
         if(south_cut && r_ref_var < 180 && muoncharge_var == 0 ){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if(bluebeam_spin_pattern == 1 ){
                  
                  if (px_var < 0) SA_BB_NM_Event_L_UP[k]=PT_H->GetBinContent(k);
                  if (px_var > 0) SA_BB_NM_Event_R_UP[k]=PT_H->GetBinContent(k);
               }
               if(bluebeam_spin_pattern == -1 ){
                  
                  if (px_var < 0) SA_BB_NM_Event_L_DOWN[k] = PT_H->GetBinContent(k);
                  if (px_var > 0) SA_BB_NM_Event_R_DOWN[k] = PT_H->GetBinContent(k);
                  
               }
            }
         }
         
         // north and yellow beam
         if(north_cut && r_ref_var < 180 && muoncharge_var == 1 ){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if(yellowbeam_spin_pattern == 1 ){
                  
                  if (px_var < 0) NO_YB_PM_Event_L_UP[k]=PT_H->GetBinContent(k);//std::cout<<"bulldog_y"<<NO_YB_PM_Event_L_UP[k]<<"\n";
                  if (px_var > 0) NO_YB_PM_Event_R_UP[k]=PT_H->GetBinContent(k);//std::cout<<"bulldog_x"<<NO_YB_PM_Event_R_UP[k]<<"\n";
               }
               if(yellowbeam_spin_pattern == -1 ){
                  
                  if (px_var < 0) NO_YB_PM_Event_L_DOWN[k] = PT_H->GetBinContent(k);
                  if (px_var > 0) NO_YB_PM_Event_R_DOWN[k] = PT_H->GetBinContent(k);
                  
               }
            }
         }
         
         if(north_cut && r_ref_var < 180 && muoncharge_var == 0){
            for(Int_t k =0; k< PT_H->GetNbinsX(); ++k){
               if(yellowbeam_spin_pattern == 1 ){
                  
                  if (px_var < 0) NO_YB_NM_Event_L_UP[k]=PT_H->GetBinContent(k);
                  if (px_var > 0) NO_YB_NM_Event_R_UP[k]=PT_H->GetBinContent(k);
               }
               if(yellowbeam_spin_pattern == -1 ){
                  
                  if (px_var < 0) NO_YB_NM_Event_L_DOWN[k] = PT_H->GetBinContent(k);
                  if (px_var > 0) NO_YB_NM_Event_R_DOWN[k] = PT_H->GetBinContent(k);
                  
               }
            }
         }
         
         
         
      }// x_f > 0 closses
   
      
      
      
      
   }//event loop closed
   
   
   
    std::cout<<"-----------------------------------------------------------------------------------------------------------------------------------------------"<<"\n";
   
   for(Int_t j=0; j<8; ++j){
      
      //south & yellow
      Float_t var_1 = sqrt(SA_YB_PM_Event_L_UP[j] * SA_YB_PM_Event_R_DOWN[j]);
      Float_t var_2 = sqrt(SA_YB_PM_Event_L_DOWN[j] * SA_YB_PM_Event_R_UP[j]);
      epsilon[j] = ((var_1 - var_2)/(var_1+var_2))/(polarization_yellow*cosPhi);
      epsilon_Error[j] = (sqrt(SA_YB_PM_Event_L_UP[j]*SA_YB_PM_Event_R_DOWN[j]*SA_YB_PM_Event_L_DOWN[j] * SA_YB_PM_Event_R_UP[j])/((var_1+var_2)*(var_1+var_2)))*sqrt(1/SA_YB_PM_Event_L_UP[j] + 1/SA_YB_PM_Event_R_DOWN[j] + 1/SA_YB_PM_Event_L_DOWN[j] + 1/SA_YB_PM_Event_R_UP[j] );
      
      
      
      Float_t var_3 = sqrt(SA_YB_NM_Event_L_UP[j] * SA_YB_NM_Event_R_DOWN[j]);
      Float_t var_4 = sqrt(SA_YB_NM_Event_L_DOWN[j] * SA_YB_NM_Event_R_UP[j]);
      epsilon2[j] = ((var_3 - var_4)/(var_3+var_4))/(polarization_yellow*cosPhi);
      epsilon2_Error[j] = (sqrt(SA_YB_NM_Event_L_UP[j]*SA_YB_NM_Event_R_DOWN[j]*SA_YB_NM_Event_L_DOWN[j] * SA_YB_NM_Event_R_UP[j])/((var_1+var_2)*(var_1+var_2)))*sqrt(1/SA_YB_NM_Event_L_UP[j] + 1/SA_YB_NM_Event_R_DOWN[j] + 1/SA_YB_NM_Event_L_DOWN[j] + 1/SA_YB_NM_Event_R_UP[j] );
      
      //North & blue
      Float_t var_5 = sqrt(NO_BB_PM_Event_L_UP[j]*NO_BB_PM_Event_R_DOWN[j]);
      Float_t var_6 = sqrt(NO_BB_PM_Event_L_DOWN[j]*NO_BB_PM_Event_R_UP[j]);
      epsilon_PM_NO[j] = ((var_5 - var_6)/(var_5+var_6))/(polarization_blue*cosPhi);
      epsilon_PM_NO_Error[j] = (sqrt(NO_BB_PM_Event_L_UP[j]*NO_BB_PM_Event_R_UP[j]*NO_BB_PM_Event_R_UP[j]*NO_BB_PM_Event_R_DOWN[j])/((var_5+var_6)*(var_5+var_6)))*sqrt(1/NO_BB_PM_Event_L_UP[j] + 1/NO_BB_PM_Event_R_DOWN[j] + 1/NO_BB_PM_Event_L_DOWN[j] + 1/ NO_BB_PM_Event_R_DOWN[j]);
      
      
      Float_t var_7 = sqrt(NO_BB_NM_Event_L_UP[j]*NO_BB_NM_Event_R_DOWN[j]);
      Float_t var_8 = sqrt(NO_BB_NM_Event_L_DOWN[j]*NO_BB_NM_Event_R_UP[j]);
      epsilon_NM_NO[j] = ((var_7 - var_8)/(var_7+var_8))/(polarization_blue*cosPhi);
      epsilon_NM_NO_Error[j] = (sqrt(NO_BB_NM_Event_L_UP[j]*NO_BB_NM_Event_R_UP[j]*NO_BB_NM_Event_R_UP[j]*NO_BB_NM_Event_R_DOWN[j])/((var_7+var_8)*(var_7+var_8)))*sqrt(1/NO_BB_NM_Event_L_UP[j] + 1/NO_BB_NM_Event_R_DOWN[j] + 1/NO_BB_NM_Event_L_DOWN[j] + 1/ NO_BB_NM_Event_R_DOWN[j]);
      
      // South and blue
      
      Float_t var_9 = sqrt(SA_BB_PM_Event_L_UP[j] * SA_BB_PM_Event_R_DOWN[j]);
      Float_t var_10 = sqrt(SA_BB_PM_Event_L_DOWN[j] * SA_BB_PM_Event_R_UP[j]);
      epsilon_SA_BB_PM[j] = ((var_9 - var_10)/(var_9+var_10))/(polarization_blue*cosPhi);
      epsilon_SA_BB_PMError[j] = (sqrt(SA_BB_PM_Event_L_UP[j]*SA_BB_PM_Event_R_DOWN[j]*SA_BB_PM_Event_L_DOWN[j] * SA_BB_PM_Event_R_UP[j])/((var_9+var_10)*(var_9+var_10)))*sqrt(1/SA_BB_PM_Event_L_UP[j] + 1/SA_BB_PM_Event_R_DOWN[j] + 1/SA_BB_PM_Event_L_DOWN[j] + 1/SA_BB_PM_Event_R_UP[j] );
      
      
      
      Float_t var_11 = sqrt(SA_BB_NM_Event_L_UP[j] * SA_BB_NM_Event_R_DOWN[j]);
      Float_t var_12 = sqrt(SA_BB_NM_Event_L_DOWN[j] * SA_BB_NM_Event_R_UP[j]);
      epsilon_SA_BB_NM[j] = ((var_11 - var_12)/(var_11+var_12))/(polarization_blue*cosPhi);
      epsilon_SA_BB_NMError[j] = (sqrt(SA_BB_NM_Event_L_UP[j]*SA_BB_NM_Event_R_DOWN[j]*SA_BB_NM_Event_L_DOWN[j] * SA_BB_NM_Event_R_UP[j])/((var_11+var_12)*(var_11+var_12)))*sqrt(1/SA_BB_NM_Event_L_UP[j] + 1/SA_BB_NM_Event_R_DOWN[j] + 1/SA_BB_NM_Event_L_DOWN[j] + 1/SA_BB_NM_Event_R_UP[j] );
      
      
      
      //North & Yellow
      
      Float_t var_13 = sqrt(NO_YB_PM_Event_L_UP[j] * NO_YB_PM_Event_R_DOWN[j]);
      cout<<var_13<<"\n";
      
      Float_t var_14 = sqrt(NO_YB_PM_Event_L_DOWN[j] * NO_YB_PM_Event_R_UP[j]);
      epsilon_NO_YB_PM[j] = ((var_13 - var_14)/(var_13+var_14))/(polarization_yellow*cosPhi);
      epsilon_NO_YB_PMError[j] = (sqrt(NO_YB_PM_Event_L_UP[j]*NO_YB_PM_Event_R_UP[j]*NO_YB_PM_Event_R_UP[j]*NO_YB_PM_Event_R_DOWN[j])/((var_13+var_14)*(var_13+var_14)))*sqrt(1/NO_YB_PM_Event_L_UP[j] + 1/NO_YB_PM_Event_R_DOWN[j] + 1/NO_YB_PM_Event_L_DOWN[j] + 1/ NO_YB_PM_Event_R_DOWN[j]);
      
     
      Float_t var_15 = sqrt(NO_YB_NM_Event_L_UP[j]*NO_YB_NM_Event_R_DOWN[j]);
      Float_t var_16 = sqrt(NO_YB_NM_Event_L_DOWN[j]*NO_YB_NM_Event_R_UP[j]);
      epsilon_NO_YB_NM[j] = ((var_15 - var_16)/(var_15 + var_16))/(polarization_yellow*cosPhi);
      epsilon_NO_YB_NMError[j] = (sqrt(NO_YB_NM_Event_L_UP[j]*NO_YB_NM_Event_R_UP[j]*NO_YB_NM_Event_R_UP[j]*NO_YB_NM_Event_R_DOWN[j])/((var_15+var_16)*(var_15+var_16)))*sqrt(1/NO_YB_NM_Event_L_UP[j] + 1/NO_YB_NM_Event_R_DOWN[j] + 1/NO_YB_NM_Event_L_DOWN[j] + 1/ NO_YB_NM_Event_R_DOWN[j]);
      
      //std::cout<<"lets take the print out of this fellow raw asymmetry in phi case when jump to next pT bin "<<j<<": \t"<< epsilon_NO_YB_NM[j] <<"\n";
   }
   
   
   
   
   for(Int_t i=0; i<8; ++i){
      for(Int_t j=0; j<12; ++j){
         
         Float_t var_1 = sqrt(PHI_SA_YB_PM_Event_R_UP[i][j] * PHI_SA_YB_PM_Event_R_DOWN[i][j]);
         Float_t var_2 = sqrt(PHI_SA_YB_PM_Event_L_DOWN[i][j] * PHI_SA_YB_PM_Event_R_UP[i][j]);
         epsilon_PHI[j] = ((var_1 - var_2)/(var_1+var_2))/(polarization_yellow*cosPhi);
        // std::cout<<"lets take the print out of this fellow raw asymmetry in phi case when jump to next pT bin "<<i<<": \t"<< epsilon_PHI[j] <<"\n";
         
      }
   }
   
   
 //-----------------------------------------------------------------------------------------------------------------------------------------------"<<"\n";
   
 
   auto gr_no_pm = new TGraphErrors(binnumpt,pt_diff,epsilon_PM_NO,pt_diffErr,epsilon_PM_NO_Error);
   gr_no_pm->SetTitle("x_F < 0 :: Y : AN, X :PT");
   gr_no_pm->SetMarkerColor(4);
   gr_no_pm->SetLineColor(4);
   gr_no_pm->SetMarkerStyle(8);
   gr_no_pm->GetXaxis()->SetTitle("p_{T}");
   gr_no_pm->GetYaxis()->SetTitle("A_{N}");
   
   
   auto gr_no_nm = new TGraphErrors(binnumpt,pt_diff,epsilon_NM_NO,pt_diffErr,epsilon_NM_NO_Error);
   gr_no_nm->SetTitle("x_F < 0 :: Y : AN, X :PT");
   gr_no_nm->SetMarkerColor(1);
   gr_no_nm->SetLineColor(1);
   gr_no_nm->SetMarkerStyle(8);
   gr_no_nm->GetXaxis()->SetTitle("p_{T}");
   gr_no_nm->GetYaxis()->SetTitle("A_{N}");
   
   auto c2 = new TCanvas("c2","Blue : AN, X :PT",200,10,700,500);
   c2->SetFillColor(0);
   c2->GetFrame()->SetFillColor(21);
   c2->GetFrame()->SetBorderSize(12);
   TMultiGraph *mg_no = new TMultiGraph();
   mg_no->Add(gr_no_pm);
   mg_no->Add(gr_no_nm);
   mg_no->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   mg_no->GetYaxis()->SetTitle("A_{N}");
      // Change the axis limits
   gPad->Modified();
   mg_no->GetXaxis()->SetLimits(0.,10.);
   mg_no->Draw("a p s ; ; 5 s=0.5");
   TLatex T0_NO;
   T0_NO.SetTextColor(2);
   T0_NO.DrawLatexNDC(.2,.30, "x_{F}<0, North Arm, Blue, A_{N}=#frac{#epsilon}{P#bulletcos#phi}");
   TLatex T1_NO;
   T1_NO.SetTextColor(4);
   T1_NO.DrawLatexNDC(.6,.20, "#mu^{+}");
   TLatex T2_NO;
   T2_NO.SetTextColor(1);
   T2_NO.DrawLatexNDC(.7,.20, "#mu^{-}");
   
//-----------------------------------------------------------------------------------------------------------------------------------------------"<<"\n";
   auto gr_no_yb_pm = new TGraphErrors(binnumpt,pt_diff,epsilon_NO_YB_PM,pt_diffErr,epsilon_NO_YB_PMError);
   gr_no_yb_pm->SetTitle("x_F < 0 :: Y : AN, X :PT");
   gr_no_yb_pm->SetMarkerColor(4);
   gr_no_yb_pm->SetLineColor(4);
   gr_no_yb_pm->SetMarkerStyle(8);
   gr_no_yb_pm->GetXaxis()->SetTitle("p_{T}");
   gr_no_yb_pm->GetYaxis()->SetTitle("A_{N}");
   
   
   auto gr_no_yb_nm = new TGraphErrors(binnumpt,pt_diff,epsilon_NO_YB_NM,pt_diffErr,epsilon_NO_YB_NMError);
   gr_no_yb_nm->SetTitle("x_F < 0 :: Y : AN, X :PT");
   gr_no_yb_nm->SetMarkerColor(1);
   gr_no_yb_nm->SetLineColor(1);
   gr_no_yb_nm->SetMarkerStyle(8);
   gr_no_yb_nm->GetXaxis()->SetTitle("p_{T}");
   gr_no_yb_nm->GetYaxis()->SetTitle("A_{N}");
   
   auto c4 = new TCanvas("c4","Y : AN, X :PT",200,10,700,500);
   c4->SetFillColor(0);
   c4->GetFrame()->SetFillColor(21);
   c4->GetFrame()->SetBorderSize(12);
   TMultiGraph *mg_no_yb_pm = new TMultiGraph();
   mg_no_yb_pm->Add(gr_no_yb_pm);
   mg_no_yb_pm->Add(gr_no_yb_nm);
   mg_no_yb_pm->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   mg_no_yb_pm->GetYaxis()->SetTitle("A_{N}");
   gPad->Modified();
   mg_no_yb_pm->GetXaxis()->SetLimits(0.,10.);
   mg_no_yb_pm->Draw("a p s ; ; 5 s=0.5");
   TLatex T0_NO_YB;
   T0_NO_YB.SetTextColor(2);
   T0_NO_YB.DrawLatexNDC(.2,.30, "x_{F}> 0, North Arm, Yellow Beam, A_{N}=#frac{#epsilon}{P#bulletcos#phi}");
   TLatex T1_NO_YB;
   T1_NO_YB.SetTextColor(4);
   T1_NO_YB.DrawLatexNDC(.6,.20, "#mu^{+}");
   TLatex T2_NO_YB;
   T2_NO_YB.SetTextColor(1);
   T2_NO_YB.DrawLatexNDC(.7,.20, "#mu^{-}");
   
   
   
   
   
   
   
   
   
   
   
   
   
   TF1 *f1 = new TF1("f1","[0]*cos(x)",1.25,10);
   
   
   auto gr_sa_bb_pm = new TGraphErrors(binnumpt,pt_diff,epsilon_SA_BB_PM,pt_diffErr,epsilon_SA_BB_PMError);
   gr_sa_bb_pm->SetTitle("x_F > 0 :: Y : AN, X :PT");
   gr_sa_bb_pm->SetMarkerColor(4);
   gr_sa_bb_pm->SetLineColor(4);
   gr_sa_bb_pm->SetMarkerStyle(8);
   gr_sa_bb_pm->GetXaxis()->SetTitle("p_{T}");
   gr_sa_bb_pm->GetYaxis()->SetTitle("A_{N}");
   
   
   auto gr_sa_bb_nm = new TGraphErrors(binnumpt,pt_diff,epsilon_SA_BB_NM,pt_diffErr,epsilon_SA_BB_NMError);
   gr_sa_bb_nm->SetTitle("x_F > 0 :: Y : AN, X :PT");
   gr_sa_bb_nm->SetMarkerColor(1);
   gr_sa_bb_nm->SetLineColor(1);
   gr_sa_bb_nm->SetMarkerStyle(8);
   gr_sa_bb_nm->GetXaxis()->SetTitle("p_{T}");
   gr_sa_bb_nm->GetYaxis()->SetTitle("A_{N}");
   
   auto c3 = new TCanvas("c3","Blue : AN, X :PT",200,10,700,500);
   c3->SetFillColor(0);
   c3->GetFrame()->SetFillColor(21);
   c3->GetFrame()->SetBorderSize(12);
   TMultiGraph *mg_sa_bb = new TMultiGraph();
   mg_sa_bb->Add(gr_sa_bb_pm);
   mg_sa_bb->Add(gr_sa_bb_nm);
   mg_sa_bb->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   mg_sa_bb->GetYaxis()->SetTitle("A_{N}");
      // Change the axis limits
   gPad->Modified();
   mg_sa_bb->GetXaxis()->SetLimits(0.,10.);
   mg_sa_bb->Draw("a p s ; ; 5 s=0.5");
   TLatex T0_SA_BB;
   T0_SA_BB.SetTextColor(2);
   T0_SA_BB.DrawLatexNDC(.2,.30, "x_{F}> 0, South Arm, Blue, A_{N}=#frac{#epsilon}{P#bulletcos#phi}");
   TLatex T1_SA_BB;
   T1_SA_BB.SetTextColor(4);
   T1_SA_BB.DrawLatexNDC(.6,.20, "#mu^{+}");
   TLatex T2_SA_BB;
   T2_SA_BB.SetTextColor(1);
   T2_SA_BB.DrawLatexNDC(.7,.20, "#mu^{-}");
   
   //-----------------------------------------------------------------------------------------------------------------------------------------------"<<"\n";
   
   
   
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
   mg->GetXaxis()->SetTitle("p_{T}");
   mg->GetYaxis()->SetTitle("A_{N}");
   gPad->Modified();
   mg->GetXaxis()->SetLimits(0.,10.);
   mg->Draw("a p s ; ; 5 s=0.5");
   TLatex T0;
   T0.SetTextColor(2);
   T0.DrawLatexNDC(.2,.30, "x_{F}<0, South Arm, Yellow, A_{N}=#frac{#epsilon}{P#bulletcos#phi}");
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
