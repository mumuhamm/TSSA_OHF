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
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

using  namespace std;


#define IS_TYPE_SIGNED(a) ((a-1) < 0)
#define MAX_VALUE_UNSIGNED(a) (((unsigned long long)1 << (sizeof(a) * CHAR_BIT)) - 1)
#define MAX_VALUE_SIGNED(a) (MAX_VALUE_UNSIGNED(a) >> 1)
#define MAX_VALUE(a) (IS_TYPE_SIGNED(a) ? MAX_VALUE_SIGNED(a) : MAX_VALUE_UNSIGNED(a))

void analyzer(string cadidatefile, string spinfile){
   
   definition call;
   
   
   std::cout<<"======================================spin tree output======================"<<"\n";
   
   TFile *mufile_spin = new TFile((spinfile).c_str());
   TTree *mutree_spin = (TTree*)mufile_spin->Get("T");
   int n_entries_spin = mutree_spin->GetEntries();
   std::cout<<" Number of entries for the moment : \t"<< n_entries_spin <<"\n";
   
   int run_spin, bunch_spin;
   float polblue, polyellow;
   int y_pattern[120], b_pattern[120], cand_run, xshift;
   Long64_t scaler[120];
   
   
   mutree_spin->SetBranchAddress("runnumber",&run_spin);
   mutree_spin->SetBranchAddress("polblue",&polblue);
   mutree_spin->SetBranchAddress("polyellow",&polyellow);
   mutree_spin->SetBranchAddress("scalerA",&scaler[0]);
   mutree_spin->SetBranchAddress("patternblue",&b_pattern[0]);
   mutree_spin->SetBranchAddress("patternyellow",&y_pattern[0]);
   mutree_spin->SetBranchAddress("xingshift",&xshift);
   
   
   
   
   std::cout<<"======================================candidate tree output======================"<<"\n";
   
   
   
   
   
   
   
   TH2F *_h2d = new TH2F("_h2d", "_h2d; x( cm); y (cm)", 100, -450, 450, 100, -450, 450);
   TH1F *mupt_h = new TH1F("mupt_h", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 10);//pt_bins
   TH1F *phi_h = new TH1F("phi_h", "phi_h; #phi (rad); Number of events", 100, -TMath::Pi(), TMath::Pi());//pt_bins
   
   TFile *mufile_cadidate = new TFile((cadidatefile).c_str());
   TTree *mutree_cadidate = (TTree*)mufile_cadidate->Get("analysis");
   int n_entries_cadidate = mutree_cadidate->GetEntries();
   std::cout<<" Number of entries for the moment : \t"<< n_entries_cadidate <<"\n";
   TLorentzVector *mu_4vec = new TLorentzVector();
   float px, py, pz, pt, pT_cut_val, rapidity, energy, mass, phi, x0, y0, z0, vtx_x, vtx_y, vtx_z, x_sta1, y_sta1, idx, idy;
   float trchi2, idchi2, dg0, ddg0, phi_trk;
   float  costheta_vtx, theta_vtx, costheta_mutr, theta_mutr, delta_theta, scaled_dtheta;
   float r_ref, chi2_trk_vtx,  sq_norm_z;
   float px_sta1, py_sta1, pz_sta1, x_sta3, y_sta3, pz_max;
   int trhits, idhits, eventYield, clock, run_candidate, lastgap;
   float ctheta, pdtheta, dangle;
   bool muoncharge;
   int clock_candidate;
   float pz_array[n_entries_cadidate];
   pz_max = pz_array[n_entries_cadidate];
   float ipx =0.0, ipy=0.0; //coordinate of IP cosidered as 0,0 , centre of mass frame
   TVector3 vec_trk;
   TVector3 vec_vtx;
   TVector3 vec_mutr;
   TVector3 three_mom;
   TVector3 vec_proj;
   TVector3 vec_vtx_z;
   
   
   mutree_cadidate->SetBranchAddress("smpx",&px);
   mutree_cadidate->SetBranchAddress("smpy",&py);
   mutree_cadidate->SetBranchAddress("smpz",&pz);
   mutree_cadidate->SetBranchAddress("smrapidity",&rapidity);
   mutree_cadidate->SetBranchAddress("smtrchi2",&trchi2);
   mutree_cadidate->SetBranchAddress("smidchi2",&idchi2);
   mutree_cadidate->SetBranchAddress("smtrhits",&trhits);
   mutree_cadidate->SetBranchAddress("smidhits",&idhits);
   mutree_cadidate->SetBranchAddress("smddg0",&ddg0);
   mutree_cadidate->SetBranchAddress("smdg0",&dg0);
   mutree_cadidate->SetBranchAddress("smx0",&x0);
   mutree_cadidate->SetBranchAddress("smy0",&y0);
   mutree_cadidate->SetBranchAddress("smz0",&z0);
   mutree_cadidate->SetBranchAddress("evtvtxx",&vtx_x);
   mutree_cadidate->SetBranchAddress("evtvtxy",&vtx_y);
   mutree_cadidate->SetBranchAddress("evtvtxz",&vtx_z);
   mutree_cadidate->SetBranchAddress("smxst1",&x_sta1);
   mutree_cadidate->SetBranchAddress("smyst1",&y_sta1);
   mutree_cadidate->SetBranchAddress("smidx",&idx);
   mutree_cadidate->SetBranchAddress("smidy",&idy);
   mutree_cadidate->SetBranchAddress("smcharge",&muoncharge);
   mutree_cadidate->SetBranchAddress("eventnumber",&eventYield);
   mutree_cadidate->SetBranchAddress("lvl1_clock_cross",&clock);
   mutree_cadidate->SetBranchAddress("smst1px",&px_sta1);
   mutree_cadidate->SetBranchAddress("smst1py",&py_sta1);
   mutree_cadidate->SetBranchAddress("smst1pz",&pz_sta1);
   mutree_cadidate->SetBranchAddress("smxst3",&x_sta3);
   mutree_cadidate->SetBranchAddress("runnumber",&run_candidate);
   mutree_cadidate->SetBranchAddress("smlastgap",&lastgap);
   mutree_cadidate->SetBranchAddress("lvl1_clock_cross",&clock_candidate);
   
   //================================================================== output tree
   float pz_var =0, rapidity_var=0, pt_var=0, phi_var=0, ddg0_var=0, dg0_var=0, pdtheta_var=0, idchi2_var=0, trchi2_var=0, x_F_var=0;
   int muoncharge_var, trhits_var, idhits_var, lastgap_var, run_candidate_var, run_spin_var;
   int bluebeam_spin_pattern, yellowbeam_spin_pattern;
   float bluebeam_pol_var, yellowbeam_pol_var;
   
   bool pt_bin1 = false;
   bool pt_bin2 = false;
   bool pt_bin3 = false;
   bool pt_bin4 = false;
   bool pt_bin5 = false;
   bool pt_bin6 = false;
   
   bool south_cut = false;
   bool north_cut = false;
   
  
   float ptSpinUp_SouthArm_bin1, ptSpinUp_SouthArm_bin2, ptSpinUp_SouthArm_bin3, ptSpinUp_SouthArm_bin4, ptSpinUp_SouthArm_bin5, ptSpinUp_SouthArm_bin6;
   float ptSpinUp_NorthArm_bin1, ptSpinUp_NorthArm_bin2, ptSpinUp_NorthArm_bin3, ptSpinUp_NorthArm_bin4, ptSpinUp_NorthArm_bin5, ptSpinUp_NorthArm_bin6;
   float ptSpinDown_SouthArm_bin1, ptSpinDown_SouthArm_bin2, ptSpinDown_SouthArm_bin3, ptSpinDown_SouthArm_bin4, ptSpinDown_SouthArm_bin5, ptSpinDown_SouthArm_bin6;
   float ptSpinDown_NorthArm_bin1, ptSpinDown_NorthArm_bin2, ptSpinDown_NorthArm_bin3, ptSpinDown_NorthArm_bin4, ptSpinDown_NorthArm_bin5, ptSpinDown_NorthArm_bin6;
   
   
   TFile *f = new TFile(("../fit_"+output+".root").c_str(),"recreate");
   TTree *fit = new TTree("fit","selected ntcltestle");
   fit->Branch("pz_var",&pz_var,"pz_var/F");
   fit->Branch("rapidity_var",&rapidity_var,"rapidity_var/F");
   fit->Branch("pt_var",&pt_var,"pt_var/F");
   fit->Branch("phi_var",&phi_var,"phi_var/F");
   fit->Branch("ddg0_var",&ddg0_var,"ddg0_var/F");
   fit->Branch("dg0_var",&dg0_var,"dg0_var/F");
   fit->Branch("pdtheta_var",&pdtheta_var,"pdtheta_var/F");
   fit->Branch("idchi2_var",&idchi2_var,"idchi2_var/F");
   fit->Branch("trchi2_var",&trchi2_var,"trchi2_var/F");
   fit->Branch("muoncharge_var",&muoncharge_var,"muoncharge_var/I");
   fit->Branch("trhits_var",&trhits_var,"trhits_var/I");
   fit->Branch("idhits_var",&idhits_var,"idhits_var/I");
   fit->Branch("x_F_var",&x_F_var,"x_F_var/F");
   fit->Branch("lastgap_var",&lastgap_var,"lastgap_var/I");
   fit->Branch("pt_bin1",&pt_bin1,"pt_bin1/O");
   fit->Branch("pt_bin2",&pt_bin2,"pt_bin2/O");
   fit->Branch("pt_bin3",&pt_bin3,"pt_bin3/O");
   fit->Branch("pt_bin4",&pt_bin4,"pt_bin4/O");
   fit->Branch("pt_bin5",&pt_bin5,"pt_bin5/O");
   fit->Branch("pt_bin6",&pt_bin6,"pt_bin6/O");
   fit->Branch("south_cut",&south_cut,"south_cut/O");
   fit->Branch("north_cut",&north_cut,"north_cut/O");
   fit->Branch("bluebeam_spin_pattern",&bluebeam_spin_pattern,"bluebeam_spin_pattern/I");
   fit->Branch("yellowbeam_spin_pattern",&yellowbeam_spin_pattern,"yellowbeam_spin_pattern/I");
   fit->Branch("bluebeam_pol_var",&bluebeam_pol_var,"bluebeam_pol_var/F");
   fit->Branch("yellowbeam_pol_var",&yellowbeam_pol_var,"yellowbeam_pol_var/F");
   fit->Branch("run_candidate_var",&run_candidate_var,"run_candidate_var/I");
   fit->Branch("run_spin_var",&run_spin_var,"run_spin_var/I");
   
  
   
   for (int ientry = 0; ientry<n_entries_cadidate ; ++ientry){//n_entries_cadidate
      
      mutree_cadidate->GetEntry(ientry);
      if (ientry%10000==0) cout << "processing event " << ientry << "/" << n_entries_cadidate<<"\n";
      run_candidate_var = run_candidate;
      
#if 1
      pz_array[ientry] = pz;
      if (pz_array[ientry] > pz_max){pz_max = pz_array[ientry];}
      x_F_var = pz/pz_max;
#endif
      
      _h2d->Fill(x_sta3, y_sta3);
      pt = call.pT(px, py);
      
      
      energy = pz*TMath::TanH(rapidity);
      mu_4vec->SetPxPyPzE(px, py, pz, energy);
      three_mom.SetXYZ(px, py, pz);
      vec_trk.SetXYZ(x0,y0,z0);
      vec_vtx.SetXYZ(vtx_x, vtx_y, vtx_z);
      vec_vtx_z.SetXYZ(0, 0, vtx_z);
      vec_mutr.SetXYZ(x_sta1, y_sta1, 0);
      
      float p = sqrt(px*px+py*py+pz*pz);
      float pSTI = sqrt(px_sta1*px_sta1 + py_sta1*py_sta1 + pz_sta1*pz_sta1);
      ctheta = (px*px_sta1 + py*py_sta1 + pz*pz_sta1) / ( p*pSTI ) ;
      if ( abs(ctheta) < 1.0 ) {
         dangle = acos (ctheta);
      }
      else{
         dangle = 0.0;
      }
      pdtheta = dangle*0.5*(p+pSTI);
      
      vec_proj = call.vectorProjection(vec_trk,vec_vtx_z );
      sq_norm_z = vec_proj.Mag2();
      if(!(sq_norm_z !=sq_norm_z)){ chi2_trk_vtx = sq_norm_z;}
      if (idx != -8888 && idy != -8888    ){
         r_ref = sqrt ((ipx-idx)*(ipx-idx) + (ipy-idy)*(ipy-idy));
      }
      
      
      pt_bin1 = (pt > 1.25 && pt < 1.50);
      pt_bin2 = (pt > 1.50 && pt < 2.00);
      pt_bin3 = (pt > 2.00 && pt < 2.50);
      pt_bin4 = (pt > 2.50 && pt < 3.00);
      pt_bin5 = (pt > 3.00 && pt < 3.50);
      pt_bin6 = (pt > 3.50 && pt < 5.00);
      south_cut = ( trhits > 12 && trchi2 < 10 && idhits > 6 && idchi2 < 5 && ddg0 < 8 && dg0 < 20 && (fabs(rapidity)>1.2 || fabs(rapidity)< 2.0) && pdtheta < 0.2);
      north_cut = ( trhits > 12 && trchi2 < 10 && idhits > 6 && idchi2 < 5 && ddg0 < 8 && dg0 < 10 && (fabs(rapidity)>1.2 || fabs(rapidity)< 2.0) && pdtheta < 0.2);
      
      pz_var = pz;
      pt_var = pt;
      ddg0_var = ddg0;
      dg0_var = dg0;
      rapidity_var = rapidity;
      phi_var = vec_trk.Phi();
      pdtheta_var = pdtheta;
      idchi2_var = idchi2;
      trchi2_var = trchi2;
      idhits_var = idhits;
      trhits_var = trhits;
      muoncharge_var = muoncharge;
      lastgap_var = lastgap;
      
      
      
      
      for (int jentry = 0; jentry<=n_entries_spin; ++jentry){
         
         mutree_spin->GetEntry(jentry);
         run_spin_var = run_spin;
         if(run_candidate != run_spin)continue;
         int shifted_clock = (clock_candidate + xshift)%120 ;
         bluebeam_spin_pattern = b_pattern[shifted_clock] *(-1);
         yellowbeam_spin_pattern = y_pattern[shifted_clock]*(-1);
         bluebeam_pol_var = polblue;
         yellowbeam_pol_var = polyellow;
            
       
      }
      
      
      fit->Fill();
   }//Event loop
   fit->Print();
   _h2d->Write();
   f->Write();
   
   
 
   
    
}
