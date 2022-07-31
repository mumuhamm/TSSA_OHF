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
   
   TH1F *mupt_h_s = new TH1F("mupt_h_s", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 7);//pt_bins
   TH1F *costheta_vtx_h_s = new TH1F("costheta_vtx_h_s","cos#theta_{vtx}; cos(#theta)_{vtx}; Number of events", 100, -1,1);
   TH1F *theta_vtx_h_s = new TH1F("theta_vtx_h_s","#theta_{vtx}; #theta_{vtx} (rad); Number of events", 100, 0, TMath::Pi());
   TH1F *costheta_mutr_h_s = new TH1F("costheta_mutr_h_s","cos#theta_{MuTr}; cos(#theta)_{MuTr}; Number of events", 100, -1,1);
   TH1F *theta_mutr_h_s = new TH1F("theta_mutr_h_s","#theta_{MuTr}; #theta_{MuTr} (rad); Number of events", 100, -TMath::Pi(), TMath::Pi());
   TH1F *delta_theta_h_s = new TH1F("delta_theta_h_s", "#delta#theta = #theta_{MuTr} - #theta_{vtx};#delta#theta (rad);Number of events", 50, -2, 2 );
   TH1F *scaled_dtheta_h_s = new TH1F("scaled_dtheta_h_s", "p#bullet(#theta_{MuTr} - #theta_{vtx});p#bullet(#theta_{MuTr} - #theta_{vtx}) (rad.GeV);Number of events", 100, -25, 25 );
   TH1F *r_ref_h_s = new TH1F("r_ref_h_s", "r_{ref} ; r_{ref} (cm); Number of events", 100, 0, 700);
   TH1F *chi2_trkzvtx_h_s = new TH1F("chi2_trkzvtx_h_s", "#chi^{2} (r_{trk}#to z_{vtx}) ; #chi^{2}; Number of events", 100, 0, 200);
   
   
   
   
   //---------------------------------------------------------------------------------------------------------------------------North Variable
   
   TH1F *mupt_h_n = new TH1F("mupt_h_n", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 7);//pt_bins
   TH1F *costheta_vtx_h_n = new TH1F("costheta_vtx_h_n","cos#theta_{vtx}; cos(#theta)_{vtx}; Number of events", 100, -1,1);
   TH1F *theta_vtx_h_n = new TH1F("theta_vtx_h_n","#theta_{vtx}; #theta_{vtx} (rad); Number of events", 100, 0, TMath::Pi());
   TH1F *costheta_mutr_h_n = new TH1F("costheta_mutr_h_n","cos#theta_{MuTr}; cos(#theta)_{MuTr}; Number of events", 100, -1,1);
   TH1F *theta_mutr_h_n = new TH1F("theta_mutr_h_n","#theta_{MuTr}; #theta_{MuTr} (rad); Number of events", 100, -TMath::Pi(), TMath::Pi());
   TH1F *delta_theta_h_n = new TH1F("delta_theta_h_n", "#delta#theta = #theta_{MuTr} - #theta_{vtx};#delta#theta (rad);Number of events", 50, -2, 2);
   TH1F *scaled_dtheta_h_n = new TH1F("scaled_dtheta_h_n", "p#bullet(#theta_{MuTr} - #theta_{vtx});p#bullet(#theta_{MuTr} - #theta_{vtx}) (rad.GeV);Number of events", 100, -25, 25 );
   TH1F *r_ref_h_n = new TH1F("r_ref_h_n", "r_{ref} ; r_{ref} (cm); Number of events", 100, 0, 700);
   TH1F *chi2_trkzvtx_h_n = new TH1F("chi2_trkzvtx_h_n", "#chi^{2} (r_{trk}#to z_{vtx}) ; #chi^{2}; Number of events", 100, 0, 200);
   
   
   
   
   TFile *mufile = new TFile((filename).c_str());
   TTree *mutree = (TTree*)mufile->Get("analysis");
   int n_entries = mutree->GetEntries();
   std::cout<<" Number of entries for the moment : \t"<< n_entries <<"\n";
   TLorentzVector *mu_4vec = new TLorentzVector();
   float px, py, pz, pt, rapidity, energy, mass, phi, x0, y0, z0, vtx_x, vtx_y, vtx_z, x_st1, y_st1, idx, idy;
   float trchi2, idchi2, dg0, ddg0;
   float pt_s, costheta_vtx_s, theta_vtx_s, costheta_mutr_s, theta_mutr_s, delta_theta_s, scaled_dtheta_s;
   float pt_n, costheta_vtx_n, theta_vtx_n, costheta_mutr_n, theta_mutr_n, delta_theta_n, scaled_dtheta_n;
   float r_ref_s, chi2_trk_vtx_s,  sq_norm_z_s;
   float r_ref_n, chi2_trk_vtx_n,  sq_norm_z_n;
   int trhits, idhits;
   float ipx =0.0, ipy=0.0; //coordinate of IP cosidered as 0,0 , centre of mass frame
   TVector3 vec_trk_s;    TVector3 vec_trk_n;
   TVector3 vec_vtx_s;    TVector3 vec_vtx_n;
   TVector3 vec_mutr_s;   TVector3 vec_mutr_n;
   TVector3 three_mom_s;  TVector3 three_mom_n;
   TVector3 vec_proj_s;   TVector3 vec_proj_n;
   TVector3 vec_vtx_z_s;  TVector3 vec_vtx_z_n;
   
   
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
   
   
   for (int i = 0; i<=n_entries; ++i){
      
      mutree->GetEntry(i);
//==================================================================================================================South
      if (pz<0){
         energy = pz*TMath::TanH(rapidity);
         mu_4vec->SetPxPyPzE(px, py, pz, energy);
         three_mom_s.SetXYZ(px, py, pz);
         vec_trk_s.SetXYZ(x0,y0,z0);
         vec_vtx_s.SetXYZ(vtx_x, vtx_y, vtx_z);
         vec_vtx_z_s.SetXYZ(0, 0, vtx_z);
         vec_mutr_s.SetXYZ(x_st1, y_st1, 0);
         pt_s = call.pT(px, py);
         costheta_vtx_s = call.costheta(vec_trk_s, vec_vtx_s);
         costheta_mutr_s = call.costheta(vec_trk_s, vec_mutr_s);
         theta_vtx_s = call.cosinverse(vec_trk_s, vec_vtx_s);
         theta_mutr_s = call.cosinverse(vec_trk_s, vec_mutr_s);
         float delta_theta_inter_s = theta_mutr_s - theta_vtx_s ;
         if(!(delta_theta_inter_s != delta_theta_inter_s))
            {
            delta_theta_s = delta_theta_inter_s;
            scaled_dtheta_s = (three_mom_s.Mag())*delta_theta_s;
            }//theta cross check
         vec_proj_s = call.vectorProjection(vec_trk_s,vec_vtx_z_s );//(vec_trk_s.Dot(vec_vtx_z_s)/vec_vtx_z_s.Mag2())*vec_vtx_z_s;
         sq_norm_z_s = vec_proj_s.Mag2();
         if(!(sq_norm_z_s !=sq_norm_z_s)){ chi2_trk_vtx_s = sq_norm_z_s;}
         if (idx != -8888 && idy != -8888    ){
            r_ref_s = sqrt ((ipx-idx)*(ipx-idx) + (ipy-idy)*(ipy-idy));
            //std::cout<<r_ref_s<<" r reference value "<<"\n";
         }
         if( trhits > 12 && trchi2 < 10 && idhits > 6 && idchi2 < 5
         && ddg0 < 8 && dg0 < 20 && (fabs(rapidity)>1.2 || fabs(rapidity)< 2.0)
         && (pt_s > 1.25 || pt_s < 10.0)){
           
            mupt_h_s->Fill(pt_s);
         }
         //mupt_h->Fill(pt);
         chi2_trkzvtx_h_s->Fill(chi2_trk_vtx_s);
         theta_vtx_h_s->Fill(theta_vtx_s);
         theta_mutr_h_s->Fill(theta_mutr_s);
         costheta_vtx_h_s->Fill(costheta_vtx_s);
         costheta_mutr_h_s->Fill(costheta_mutr_s);
         delta_theta_h_s->Fill(delta_theta_s);
         scaled_dtheta_h_s->Fill(scaled_dtheta_s);
         r_ref_h_s->Fill(r_ref_s);
         
         
      }//pz - south
      
//=================================================================================================================North
      if (pz > 0 ){
         energy = pz*TMath::TanH(rapidity);
         mu_4vec->SetPxPyPzE(px, py, pz, energy);
         three_mom_n.SetXYZ(px, py, pz);
         vec_trk_n.SetXYZ(x0,y0,z0);
         vec_vtx_n.SetXYZ(vtx_x, vtx_y, vtx_z);
         vec_vtx_z_n.SetXYZ(0, 0, vtx_z);
         vec_mutr_n.SetXYZ(x_st1, y_st1, 0);
         pt_n = call.pT(px, py);
         costheta_vtx_n = call.costheta(vec_trk_n, vec_vtx_n);
         costheta_mutr_n = call.costheta(vec_trk_n, vec_mutr_n);
         theta_vtx_n = call.cosinverse(vec_trk_n, vec_vtx_n);
         theta_mutr_n = call.cosinverse(vec_trk_n, vec_mutr_n);
         float delta_theta_inter_n = theta_mutr_n - theta_vtx_n ;
         if(!(delta_theta_inter_n != delta_theta_inter_n))
            {
            delta_theta_n = delta_theta_inter_n;
            scaled_dtheta_n = (three_mom_n.Mag())*delta_theta_n;
            }//theta cross check
         vec_proj_n = call.vectorProjection(vec_trk_n,vec_vtx_z_n );//(vec_trk_n.Dot(vec_vtx_z_n)/vec_vtx_z_n.Mag2())*vec_vtx_z_n;
         sq_norm_z_n = vec_proj_n.Mag2();
         if(!(sq_norm_z_n !=sq_norm_z_n)){ chi2_trk_vtx_n = sq_norm_z_n;}
         if (idx != -8888 && idy != -8888    ){ r_ref_n = sqrt ((ipx-idx)*(ipx-idx) + (ipy-idy)*(ipy-idy));}
         if( trhits > 12 && trchi2 < 10 && idhits > 6 && idchi2 < 5
            && ddg0 < 8 && dg0 < 10 && (fabs(rapidity)>1.2 || fabs(rapidity)< 2.0)
            && (pt_n > 1.25 || pt_n < 10.0)){
            
            mupt_h_n->Fill(pt_n);
         }
            //mupt_h->Fill(pt);
         chi2_trkzvtx_h_n->Fill(chi2_trk_vtx_n);
         theta_vtx_h_n->Fill(theta_vtx_n);
         theta_mutr_h_n->Fill(theta_mutr_n);
         costheta_vtx_h_n->Fill(costheta_vtx_n);
         costheta_mutr_h_n->Fill(costheta_mutr_n);
         delta_theta_h_n->Fill(delta_theta_n);
         scaled_dtheta_h_n->Fill(scaled_dtheta_n);
         r_ref_h_n->Fill(r_ref_n);
         
         
      }//pz
      
      
   }//Event loop
   gStyle->SetOptStat(0);
   
   call.plot(delta_theta_h_s, delta_theta_h_n, "delta_theta");
   call.plot(costheta_vtx_h_s, costheta_vtx_h_n, "costheta_vtx");
   call.plot(theta_vtx_h_s, theta_vtx_h_n, "theta_vtx");
   call.plot(costheta_mutr_h_s, costheta_mutr_h_n, "costheta_mutr");
   call.plot(theta_mutr_h_s, theta_mutr_h_n, "theta_mutr");
   call.plot(scaled_dtheta_h_s, scaled_dtheta_h_n, "scaled_dtheta");
   call.plot(r_ref_h_s, r_ref_h_n, "r_ref");
   call.plot(chi2_trkzvtx_h_s, chi2_trkzvtx_h_n, "chi2_trkzvtx");
   call.plot(mupt_h_s, mupt_h_n, "mupt");
   
   

}
