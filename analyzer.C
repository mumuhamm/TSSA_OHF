/* Author : Muhammad Alibordi

 */

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
#include "analyzer.h"

using  namespace std;

void analyzer(string filename){
   
   definition call;
   
   float pt_bins[11] = {0,1,2,3,4,5,6,7,8,9,10};
   
      //---------------------------------------------------------------------------------------------------------------------------South Variable
   
   TH1F *mupt_h_s_mp = new TH1F("mupt_h_s_mp", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 7);//pt_bins
   TH1F *costheta_vtx_h_s_mp = new TH1F("costheta_vtx_h_s_mp","cos#theta_{vtx}; cos(#theta)_{vtx}; Number of events", 100, -1,1);
   TH1F *theta_vtx_h_s_mp = new TH1F("theta_vtx_h_s_mp","#theta_{vtx}; #theta_{vtx} (rad); Number of events", 100, 0, TMath::Pi());
   TH1F *costheta_mutr_h_s_mp = new TH1F("costheta_mutr_h_s_mp","cos#theta_{MuTr}; cos(#theta)_{MuTr}; Number of events", 100, -1,1);
   TH1F *theta_mutr_h_s_mp = new TH1F("theta_mutr_h_s_mp","#theta_{MuTr}; #theta_{MuTr} (rad); Number of events", 100, -TMath::Pi(), TMath::Pi());
   TH1F *delta_theta_h_s_mp = new TH1F("delta_theta_h_s_mp", "#delta#theta = #theta_{MuTr} - #theta_{vtx};#delta#theta (rad);Number of events", 50, -2, 2 );
   TH1F *scaled_dtheta_h_s_mp = new TH1F("scaled_dtheta_h_s_mp", "p#bullet(#theta_{MuTr} - #theta_{vtx});p#bullet(#theta_{MuTr} - #theta_{vtx}) (rad.GeV);Number of events", 100, -25, 25 );
   TH1F *r_ref_h_s_mp = new TH1F("r_ref_h_s_mp", "r_{ref} ; r_{ref} (cm); Number of events", 100, 0, 700);
   TH1F *chi2_trkzvtx_h_s_mp = new TH1F("chi2_trkzvtx_h_s_mp", "#chi^{2} (r_{trk}#to z_{vtx}) ; #chi^{2}; Number of events", 100, 0, 200);
   TH1F *phi_trk_h_s_mp = new TH1F("phi_trk_h_s_mp","#phi_{trk}; #phi_{trk} (rad); Number of events", 100,-TMath::Pi(), TMath::Pi());
   
   
   
   TH1F *mupt_h_s_mn = new TH1F("mupt_h_s_mn", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 7);//pt_bins
   TH1F *costheta_vtx_h_s_mn = new TH1F("costheta_vtx_h_s_mn","cos#theta_{vtx}; cos(#theta)_{vtx}; Number of events", 100, -1,1);
   TH1F *theta_vtx_h_s_mn = new TH1F("theta_vtx_h_s_mn","#theta_{vtx}; #theta_{vtx} (rad); Number of events", 100, 0, TMath::Pi());
   TH1F *costheta_mutr_h_s_mn = new TH1F("costheta_mutr_h_s_mn","cos#theta_{MuTr}; cos(#theta)_{MuTr}; Number of events", 100, -1,1);
   TH1F *theta_mutr_h_s_mn = new TH1F("theta_mutr_h_s_mn","#theta_{MuTr}; #theta_{MuTr} (rad); Number of events", 100, -TMath::Pi(), TMath::Pi());
   TH1F *delta_theta_h_s_mn = new TH1F("delta_theta_h_s_mn", "#delta#theta = #theta_{MuTr} - #theta_{vtx};#delta#theta (rad);Number of events", 50, -2, 2 );
   TH1F *scaled_dtheta_h_s_mn = new TH1F("scaled_dtheta_h_s_mn", "p#bullet(#theta_{MuTr} - #theta_{vtx});p#bullet(#theta_{MuTr} - #theta_{vtx}) (rad.GeV);Number of events", 100, -25, 25 );
   TH1F *r_ref_h_s_mn = new TH1F("r_ref_h_s_mn", "r_{ref} ; r_{ref} (cm); Number of events", 100, 0, 700);
   TH1F *chi2_trkzvtx_h_s_mn = new TH1F("chi2_trkzvtx_h_s_mn", "#chi^{2} (r_{trk}#to z_{vtx}) ; #chi^{2}; Number of events", 100, 0, 200);
   TH1F *phi_trk_h_s_mn = new TH1F("phi_trk_h_s_mn","#phi_{trk}; #phi_{trk} (rad); Number of events", 100,-TMath::Pi(), TMath::Pi());
   
      //---------------------------------------------------------------------------------------------------------------------------North Variable
   
   TH1F *mupt_h_n_mp = new TH1F("mupt_h_n_mp", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 7);//pt_bins
   TH1F *costheta_vtx_h_n_mp = new TH1F("costheta_vtx_h_n_mp","cos#theta_{vtx}; cos(#theta)_{vtx}; Number of events", 100, -1,1);
   TH1F *theta_vtx_h_n_mp = new TH1F("theta_vtx_h_n_mp","#theta_{vtx}; #theta_{vtx} (rad); Number of events", 100, 0, TMath::Pi());
   TH1F *costheta_mutr_h_n_mp = new TH1F("costheta_mutr_h_n_mp","cos#theta_{MuTr}; cos(#theta)_{MuTr}; Number of events", 100, -1,1);
   TH1F *theta_mutr_h_n_mp = new TH1F("theta_mutr_h_n_mp","#theta_{MuTr}; #theta_{MuTr} (rad); Number of events", 100, -TMath::Pi(), TMath::Pi());
   TH1F *delta_theta_h_n_mp = new TH1F("delta_theta_h_n_mp", "#delta#theta = #theta_{MuTr} - #theta_{vtx};#delta#theta (rad);Number of events", 50, -2, 2);
   TH1F *scaled_dtheta_h_n_mp = new TH1F("scaled_dtheta_h_n_mp", "p#bullet(#theta_{MuTr} - #theta_{vtx});p#bullet(#theta_{MuTr} - #theta_{vtx}) (rad.GeV);Number of events", 100, -25, 25 );
   TH1F *r_ref_h_n_mp = new TH1F("r_ref_h_n_mp", "r_{ref} ; r_{ref} (cm); Number of events", 100, 0, 700);
   TH1F *chi2_trkzvtx_h_n_mp = new TH1F("chi2_trkzvtx_h_n_mp", "#chi^{2} (r_{trk}#to z_{vtx}) ; #chi^{2}; Number of events", 100, 0, 200);
   TH1F *phi_trk_h_n_mp = new TH1F("phi_trk_h_n_mp","#phi_{trk}; #phi_{trk} (rad); Number of events", 100,-TMath::Pi(), TMath::Pi());
   
   
   
   TH1F *mupt_h_n_mn = new TH1F("mupt_h_n_mn", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 7);//pt_bins
   TH1F *costheta_vtx_h_n_mn = new TH1F("costheta_vtx_h_n_mn","cos#theta_{vtx}; cos(#theta)_{vtx}; Number of events", 100, -1,1);
   TH1F *theta_vtx_h_n_mn = new TH1F("theta_vtx_h_n_mn","#theta_{vtx}; #theta_{vtx} (rad); Number of events", 100, 0, TMath::Pi());
   TH1F *costheta_mutr_h_n_mn = new TH1F("costheta_mutr_h_n_mn","cos#theta_{MuTr}; cos(#theta)_{MuTr}; Number of events", 100, -1,1);
   TH1F *theta_mutr_h_n_mn = new TH1F("theta_mutr_h_n_mn","#theta_{MuTr}; #theta_{MuTr} (rad); Number of events", 100, -TMath::Pi(), TMath::Pi());
   TH1F *delta_theta_h_n_mn = new TH1F("delta_theta_h_n_mn", "#delta#theta = #theta_{MuTr} - #theta_{vtx};#delta#theta (rad);Number of events", 50, -2, 2);
   TH1F *scaled_dtheta_h_n_mn = new TH1F("scaled_dtheta_h_n_mn", "p#bullet(#theta_{MuTr} - #theta_{vtx});p#bullet(#theta_{MuTr} - #theta_{vtx}) (rad.GeV);Number of events", 100, -25, 25 );
   TH1F *r_ref_h_n_mn = new TH1F("r_ref_h_n_mn", "r_{ref} ; r_{ref} (cm); Number of events", 100, 0, 700);
   TH1F *chi2_trkzvtx_h_n_mn = new TH1F("chi2_trkzvtx_h_n_mn", "#chi^{2} (r_{trk}#to z_{vtx}) ; #chi^{2}; Number of events", 100, 0, 200);
   TH1F *phi_trk_h_n_mn = new TH1F("phi_trk_h_n_mn","#phi_{trk}; #phi_{trk} (rad); Number of events", 100,-TMath::Pi(), TMath::Pi());
   
   
  
   
   TFile *mufile = new TFile((filename).c_str());
   TTree *mutree = (TTree*)mufile->Get("analysis");
   int n_entries = mutree->GetEntries();
   std::cout<<" Number of entries for the moment : \t"<< n_entries <<"\n";
   TLorentzVector *mu_4vec = new TLorentzVector();
   float px, py, pz, pt, pT_cut_val, rapidity, energy, mass, phi, x0, y0, z0, vtx_x, vtx_y, vtx_z, x_st1, y_st1, idx, idy;
   float trchi2, idchi2, dg0, ddg0, phi_trk;
   float  costheta_vtx, theta_vtx, costheta_mutr, theta_mutr, delta_theta, scaled_dtheta;
   float r_ref, chi2_trk_vtx,  sq_norm_z;
   int trhits, idhits;
   bool muoncharge;
   float ipx =0.0, ipy=0.0; //coordinate of IP cosidered as 0,0 , centre of mass frame
   TVector3 vec_trk;
   TVector3 vec_vtx;
   TVector3 vec_mutr;
   TVector3 three_mom;
   TVector3 vec_proj;
   TVector3 vec_vtx_z;
   
   
   mutree->SetBranchAddress("smpx",&px);
   mutree->SetBranchAddress("smpy",&py);
   mutree->SetBranchAddress("smpz",&pz);
   mutree->SetBranchAddress("smrapidity",&rapidity);
   mutree->SetBranchAddress("smtrchi2",&trchi2);
   mutree->SetBranchAddress("smidchi2",&idchi2);
   mutree->SetBranchAddress("smtrhits",&trhits);
   mutree->SetBranchAddress("smidhits",&idhits);
   mutree->SetBranchAddress("smddg0",&ddg0);
   mutree->SetBranchAddress("smdg0",&dg0);
   mutree->SetBranchAddress("smx0",&x0);
   mutree->SetBranchAddress("smy0",&y0);
   mutree->SetBranchAddress("smz0",&z0);
   mutree->SetBranchAddress("evtvtxx",&vtx_x);
   mutree->SetBranchAddress("evtvtxy",&vtx_y);
   mutree->SetBranchAddress("evtvtxz",&vtx_z);
   mutree->SetBranchAddress("smxst1",&x_st1);
   mutree->SetBranchAddress("smyst1",&y_st1);
   mutree->SetBranchAddress("smidx",&idx);
   mutree->SetBranchAddress("smidy",&idy);
   mutree->SetBranchAddress("smcharge",&muoncharge);
   
   
   for (int ientry = 0; ientry<=n_entries; ++ientry){
      
      mutree->GetEntry(ientry);
      if (ientry%10000==0) cout << "processing event " << ientry << "/" << n_entries <<"\n";
         //std::cout<<"muon charge :::     "<<muoncharge<<"\n";
      
      energy = pz*TMath::TanH(rapidity);
      mu_4vec->SetPxPyPzE(px, py, pz, energy);
      three_mom.SetXYZ(px, py, pz);
      vec_trk.SetXYZ(x0,y0,z0);
      vec_vtx.SetXYZ(vtx_x, vtx_y, vtx_z);
      vec_vtx_z.SetXYZ(0, 0, vtx_z);
      vec_mutr.SetXYZ(x_st1, y_st1, 0);
      pt = call.pT(px, py);
      costheta_vtx = call.costheta(vec_trk, vec_vtx);
      costheta_mutr = call.costheta(vec_trk, vec_mutr);
      theta_vtx = call.cosinverse(vec_trk, vec_vtx);
      theta_mutr = call.cosinverse(vec_trk, vec_mutr);
      float delta_theta_inter = theta_mutr - theta_vtx ;
      if(!(delta_theta_inter != delta_theta_inter))
         {
         delta_theta = delta_theta_inter;
         scaled_dtheta = (three_mom.Mag())*delta_theta;
         }//theta cross check
      vec_proj = call.vectorProjection(vec_trk,vec_vtx_z );
      sq_norm_z = vec_proj.Mag2();
      if(!(sq_norm_z !=sq_norm_z)){ chi2_trk_vtx = sq_norm_z;}
      if (idx != -8888 && idy != -8888    ){
         r_ref = sqrt ((ipx-idx)*(ipx-idx) + (ipy-idy)*(ipy-idy));
      }
      
      
      if( trhits > 12 && trchi2 < 10 && idhits > 6 && idchi2 < 5
         && ddg0 < 8 && dg0 < 20 && (fabs(rapidity)>1.2 || fabs(rapidity)< 2.0)
         && (pt > 1.25 || pt < 5.0)){
         
         pT_cut_val = pt;
         phi_trk = vec_trk.Phi();
         //std::cout << " phi value   ::     "<<phi_trk<<"\n";
      }
      
      
         //==================================================================================================================South
      if (pz<0){
         if(muoncharge ==0){
            
            mupt_h_s_mn->Fill(pT_cut_val);
            chi2_trkzvtx_h_s_mn->Fill(chi2_trk_vtx);
            theta_vtx_h_s_mn->Fill(theta_vtx);
            theta_mutr_h_s_mn->Fill(theta_mutr);
            costheta_vtx_h_s_mn->Fill(costheta_vtx);
            costheta_mutr_h_s_mn->Fill(costheta_mutr);
            delta_theta_h_s_mn->Fill(delta_theta);
            scaled_dtheta_h_s_mn->Fill(scaled_dtheta);
            r_ref_h_s_mn->Fill(r_ref);
            phi_trk_h_s_mn->Fill(phi_trk);
         }//negative muo charge
         
         if(muoncharge==1){
            mupt_h_s_mp->Fill(pT_cut_val);
            chi2_trkzvtx_h_s_mp->Fill(chi2_trk_vtx);
            theta_vtx_h_s_mp->Fill(theta_vtx);
            theta_mutr_h_s_mp->Fill(theta_mutr);
            costheta_vtx_h_s_mp->Fill(costheta_vtx);
            costheta_mutr_h_s_mp->Fill(costheta_mutr);
            delta_theta_h_s_mp->Fill(delta_theta);
            scaled_dtheta_h_s_mp->Fill(scaled_dtheta);
            r_ref_h_s_mp->Fill(r_ref);
            phi_trk_h_s_mp->Fill(phi_trk);
         }
         
      }//pz - south
      
         //=================================================================================================================North
      if (pz > 0 ){
         if(muoncharge==0){
            mupt_h_n_mn->Fill(pT_cut_val);
            chi2_trkzvtx_h_n_mn->Fill(chi2_trk_vtx);
            theta_vtx_h_n_mn->Fill(theta_vtx);
            theta_mutr_h_n_mn->Fill(theta_mutr);
            costheta_vtx_h_n_mn->Fill(costheta_vtx);
            costheta_mutr_h_n_mn->Fill(costheta_mutr);
            delta_theta_h_n_mn->Fill(delta_theta);
            scaled_dtheta_h_n_mn->Fill(scaled_dtheta);
            r_ref_h_n_mn->Fill(r_ref);
            phi_trk_h_n_mn->Fill(phi_trk);
         }
         if(muoncharge==1){
            mupt_h_n_mp->Fill(pT_cut_val);
            chi2_trkzvtx_h_n_mp->Fill(chi2_trk_vtx);
            theta_vtx_h_n_mp->Fill(theta_vtx);
            theta_mutr_h_n_mp->Fill(theta_mutr);
            costheta_vtx_h_n_mp->Fill(costheta_vtx);
            costheta_mutr_h_n_mp->Fill(costheta_mutr);
            delta_theta_h_n_mp->Fill(delta_theta);
            scaled_dtheta_h_n_mp->Fill(scaled_dtheta);
            r_ref_h_n_mp->Fill(r_ref);
            phi_trk_h_n_mp->Fill(phi_trk);
         }
         
         
         
      }//pz
      
      
   }//Event loop
   gStyle->SetOptStat(0);
   
   
   call.plot_PT(mupt_h_s_mp, mupt_h_s_mn, mupt_h_n_mp, mupt_h_n_mn, "pTDistribution");
   
   call.plot_south(delta_theta_h_s_mp, delta_theta_h_s_mn, "delta_theta_south");
   call.plot_south(costheta_vtx_h_s_mp, costheta_vtx_h_s_mn, "costheta_vtx_south");
   call.plot_south(theta_vtx_h_s_mp, theta_vtx_h_s_mn, "theta_vtx_south");
   call.plot_south(costheta_mutr_h_s_mp, costheta_mutr_h_s_mn, "costheta_mutr_south");
   call.plot_south(theta_mutr_h_s_mp, theta_mutr_h_s_mn, "theta_mutr_south");
   call.plot_south(scaled_dtheta_h_s_mp, scaled_dtheta_h_s_mn, "scaled_dtheta_south");
   call.plot_south(r_ref_h_s_mp, r_ref_h_s_mn, "r_ref_south");
   call.plot_south(chi2_trkzvtx_h_s_mp, chi2_trkzvtx_h_s_mn, "chi2_trkzvtx_south");
   call.plot_south(phi_trk_h_s_mp, phi_trk_h_s_mn, "phi_trk_south");
   call.plot_trmom_south(mupt_h_s_mp, mupt_h_s_mn, "mupt_south");
   
   call.plot_north(delta_theta_h_n_mp, delta_theta_h_n_mn, "delta_theta_north");
   call.plot_north(costheta_vtx_h_n_mp, costheta_vtx_h_n_mn, "costheta_vtx_north");
   call.plot_north(theta_vtx_h_n_mp, theta_vtx_h_n_mn, "theta_vtx_north");
   call.plot_north(costheta_mutr_h_n_mp, costheta_mutr_h_n_mn, "costheta_mutr_north");
   call.plot_north(theta_mutr_h_n_mp, theta_mutr_h_n_mn, "theta_mutr_north");
   call.plot_north(scaled_dtheta_h_n_mp, scaled_dtheta_h_n_mn, "scaled_dtheta_north");
   call.plot_north(r_ref_h_n_mp, r_ref_h_n_mn, "r_ref_north");
   call.plot_north(chi2_trkzvtx_h_n_mp, chi2_trkzvtx_h_n_mn, "chi2_trkzvtx_north");
   call.plot_north(phi_trk_h_n_mp, phi_trk_h_n_mn, "phi_trk_north");
   call.plot_trmom_north(mupt_h_n_mp, mupt_h_n_mn, "mupt_north");
   
   
   
 
}
