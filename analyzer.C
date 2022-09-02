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

void analyzer(string filename){
   
   definition call;
   
   
   float pt_bins[11] = {0,1,2,3,4,5,6,7,8,9,10};
   
      //---------------------------------------------------------------------------------------------------------------------------South Variable
   
   TH1F *mupt_h_s_mp = new TH1F("mupt_h_s_mp", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 7);//pt_bins
   TH1F *mudg0_h_s_mp = new TH1F("mudg0_h_s_mp", "mudg0_h_s_mp; DG0_{#mu}; Number of events", 100, 0, 80);
   TH1F *muddg0_h_s_mp = new TH1F("muddg0_h_s_mp", "muddg0_h_s_mp;DDG0_{#mu}; Number of events (a.u.)", 100, 0, 25.0);
   TH1F *costheta_vtx_h_s_mp = new TH1F("costheta_vtx_h_s_mp","cos#theta_{vtx}; cos(#theta)_{vtx}; Number of events", 100, -1,1);
   TH1F *theta_vtx_h_s_mp = new TH1F("theta_vtx_h_s_mp","#theta_{vtx}; #theta_{vtx} (rad); Number of events", 100, 0, TMath::Pi());
   TH1F *costheta_mutr_h_s_mp = new TH1F("costheta_mutr_h_s_mp","cos#theta_{MuTr}; cos(#theta)_{MuTr}; Number of events", 100, -1,1);
   TH1F *theta_mutr_h_s_mp = new TH1F("theta_mutr_h_s_mp","#theta_{MuTr}; #theta_{MuTr} (rad); Number of events", 100, -TMath::Pi(), TMath::Pi());
   TH1F *delta_theta_h_s_mp = new TH1F("delta_theta_h_s_mp", "#delta#theta = #theta_{MuTr} - #theta_{vtx};#delta#theta (rad);Number of events", 50, 0.0, 0.3  );
   TH1F *scaled_dtheta_h_s_mp = new TH1F("scaled_dtheta_h_s_mp", "p#bullet#delta#theta;p#bullet#delta#theta (rad.GeV);Number of events", 100, 0.0, 0.4 );
   TH1F *r_ref_h_s_mp = new TH1F("r_ref_h_s_mp", "r_{ref} ; r_{ref} (cm); Number of events", 100, 0, 700);
   TH1F *chi2_trkzvtx_h_s_mp = new TH1F("chi2_trkzvtx_h_s_mp", "#chi^{2} (r_{trk}#to z_{vtx}) ; #chi^{2}; Number of events", 100, 0, 200);
   TH1F *phi_trk_h_s_mp = new TH1F("phi_trk_h_s_mp","#phi_{trk}; #phi_{trk} (rad); Number of events", 100,-TMath::Pi(), TMath::Pi());
   
   
   
   TH1F *mupt_h_s_mn = new TH1F("mupt_h_s_mn", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 7);//pt_bins
   TH1F *mudg0_h_s_mn = new TH1F("mudg0_h_s_mn", "mudg0_h_s_mn; DG0_{#mu}; Number of events", 100, 0, 80);
   TH1F *muddg0_h_s_mn = new TH1F("muddg0_h_s_mn", "muddg0_h_s_mn;DDG0_{#mu}; Number of events (a.u.)", 100, 0, 25.0);
   TH1F *costheta_vtx_h_s_mn = new TH1F("costheta_vtx_h_s_mn","cos#theta_{vtx}; cos(#theta)_{vtx}; Number of events", 100, -1,1);
   TH1F *theta_vtx_h_s_mn = new TH1F("theta_vtx_h_s_mn","#theta_{vtx}; #theta_{vtx} (rad); Number of events", 100, 0, TMath::Pi());
   TH1F *costheta_mutr_h_s_mn = new TH1F("costheta_mutr_h_s_mn","cos#theta_{MuTr}; cos(#theta)_{MuTr}; Number of events", 100, -1,1);
   TH1F *theta_mutr_h_s_mn = new TH1F("theta_mutr_h_s_mn","#theta_{MuTr}; #theta_{MuTr} (rad); Number of events", 100, -TMath::Pi(), TMath::Pi());
   TH1F *delta_theta_h_s_mn = new TH1F("delta_theta_h_s_mn", "#delta#theta = #theta_{MuTr} - #theta_{vtx};#delta#theta (rad);Number of events", 50, 0.0, 0.3  );
   TH1F *scaled_dtheta_h_s_mn = new TH1F("scaled_dtheta_h_s_mn", "p#bullet#delta#theta;p#bullet#delta#theta (rad.GeV);Number of events", 100, 0.0, 0.4 );
   TH1F *r_ref_h_s_mn = new TH1F("r_ref_h_s_mn", "r_{ref} ; r_{ref} (cm); Number of events", 100, 0, 700);
   TH1F *chi2_trkzvtx_h_s_mn = new TH1F("chi2_trkzvtx_h_s_mn", "#chi^{2} (r_{trk}#to z_{vtx}) ; #chi^{2}; Number of events", 100, 0, 200);
   TH1F *phi_trk_h_s_mn = new TH1F("phi_trk_h_s_mn","#phi_{trk}; #phi_{trk} (rad); Number of events", 100,-TMath::Pi(), TMath::Pi());
   
      //---------------------------------------------------------------------------------------------------------------------------North Variable
   
   TH1F *mupt_h_n_mp = new TH1F("mupt_h_n_mp", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 7);//pt_bins
   TH1F *mudg0_h_n_mp = new TH1F("mudg0_h_n_mp", "mudg0_h_n_mp; DG0_{#mu}; Number of events", 100, 0, 80);
   TH1F *muddg0_h_n_mp = new TH1F("muddg0_h_n_mp", "muddg0_h_n_mp;DDG0_{#mu}; Number of events (a.u.)", 100, 0, 25.0);
   TH1F *costheta_vtx_h_n_mp = new TH1F("costheta_vtx_h_n_mp","cos#theta_{vtx}; cos(#theta)_{vtx}; Number of events", 100, -1,1);
   TH1F *theta_vtx_h_n_mp = new TH1F("theta_vtx_h_n_mp","#theta_{vtx}; #theta_{vtx} (rad); Number of events", 100, 0, TMath::Pi());
   TH1F *costheta_mutr_h_n_mp = new TH1F("costheta_mutr_h_n_mp","cos#theta_{MuTr}; cos(#theta)_{MuTr}; Number of events", 100, -1,1);
   TH1F *theta_mutr_h_n_mp = new TH1F("theta_mutr_h_n_mp","#theta_{MuTr}; #theta_{MuTr} (rad); Number of events", 100, -TMath::Pi(), TMath::Pi());
   TH1F *delta_theta_h_n_mp = new TH1F("delta_theta_h_n_mp", "#delta#theta = #theta_{MuTr} - #theta_{vtx};#delta#theta (rad);Number of events", 50,0.0, 0.3  );
   TH1F *scaled_dtheta_h_n_mp = new TH1F("scaled_dtheta_h_n_mp", "p#bullet#delta#theta;p#bullet#delta#theta (rad.GeV);Number of events", 100, 0.0, 0.4);
   TH1F *r_ref_h_n_mp = new TH1F("r_ref_h_n_mp", "r_{ref} ; r_{ref} (cm); Number of events", 100, 0, 700);
   TH1F *chi2_trkzvtx_h_n_mp = new TH1F("chi2_trkzvtx_h_n_mp", "#chi^{2} (r_{trk}#to z_{vtx}) ; #chi^{2}; Number of events", 100, 0, 200);
   TH1F *phi_trk_h_n_mp = new TH1F("phi_trk_h_n_mp","#phi_{trk}; #phi_{trk} (rad); Number of events", 100,-TMath::Pi(), TMath::Pi());
   
   
   
   TH1F *mupt_h_n_mn = new TH1F("mupt_h_n_mn", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 7);//pt_bins
   TH1F *mudg0_h_n_mn = new TH1F("mudg0_h_n_mn", "mudg0_h_n_mn; DG0_{#mu}; Number of events", 100, 0, 80);
   TH1F *muddg0_h_n_mn = new TH1F("muddg0_h_n_mn", "muddg0_h_n_mn;DDG0_{#mu}; Number of events (a.u.)", 100, 0, 25.0);
   TH1F *costheta_vtx_h_n_mn = new TH1F("costheta_vtx_h_n_mn","cos#theta_{vtx}; cos(#theta)_{vtx}; Number of events", 100, -1,1);
   TH1F *theta_vtx_h_n_mn = new TH1F("theta_vtx_h_n_mn","#theta_{vtx}; #theta_{vtx} (rad); Number of events", 100, 0, TMath::Pi());
   TH1F *costheta_mutr_h_n_mn = new TH1F("costheta_mutr_h_n_mn","cos#theta_{MuTr}; cos(#theta)_{MuTr}; Number of events", 100, -1,1);
   TH1F *theta_mutr_h_n_mn = new TH1F("theta_mutr_h_n_mn","#theta_{MuTr}; #theta_{MuTr} (rad); Number of events", 100, -TMath::Pi(), TMath::Pi());
   TH1F *delta_theta_h_n_mn = new TH1F("delta_theta_h_n_mn", "#delta#theta = #theta_{MuTr} - #theta_{vtx};#delta#theta (rad);Number of events", 50, 0.0, 0.3  );
   TH1F *scaled_dtheta_h_n_mn = new TH1F("scaled_dtheta_h_n_mn", "p#bullet#delta#theta;p#bullet#delta#theta (rad.GeV);Number of events", 100, 0.0, 0.4 );
   TH1F *r_ref_h_n_mn = new TH1F("r_ref_h_n_mn", "r_{ref} ; r_{ref} (cm); Number of events", 100, 0, 700);
   TH1F *chi2_trkzvtx_h_n_mn = new TH1F("chi2_trkzvtx_h_n_mn", "#chi^{2} (r_{trk}#to z_{vtx}) ; #chi^{2}; Number of events", 100, 0, 200);
   TH1F *phi_trk_h_n_mn = new TH1F("phi_trk_h_n_mn","#phi_{trk}; #phi_{trk} (rad); Number of events", 100,-TMath::Pi(), TMath::Pi());
   
   TH1F *xF_s_h = new TH1F("xF_s_h","xF; South x_{F},  p_{z} < 0 ; Number of events", 100,-0.2, 0.0);
   TH1F *xF_n_h = new TH1F("xF_n_h","xF; North x_{F},  p_{z} > 0 ; Number of events", 100,0.0, 0.2);
   TH2F *_h2d = new TH2F("_h2d", "_h2d; x( cm); y (cm)", 100, -450, 450, 100, -450, 450);
   TH1F *_h = new TH1F("_h", "_h", 50, 0, 7);
   TH1F *mupt_h = new TH1F("mupt_h", "mupt_h; p_{T} (GeV); Number of events", 50, 0, 10);//pt_bins
   
   TFile *mufile = new TFile((filename).c_str());
   TTree *mutree = (TTree*)mufile->Get("analysis");
   int n_entries = mutree->GetEntries();
   std::cout<<" Number of entries for the moment : \t"<< n_entries <<"\n";
   TLorentzVector *mu_4vec = new TLorentzVector();
   float px, py, pz, pt, pT_cut_val, rapidity, energy, mass, phi, x0, y0, z0, vtx_x, vtx_y, vtx_z, x_sta1, y_sta1, idx, idy;
   float trchi2, idchi2, dg0, ddg0, phi_trk;
   float  costheta_vtx, theta_vtx, costheta_mutr, theta_mutr, delta_theta, scaled_dtheta;
   float r_ref, chi2_trk_vtx,  sq_norm_z;
   float px_sta1, py_sta1, pz_sta1, x_sta3, y_sta3, pz_max;
   int trhits, idhits, eventYield, clock, runnumber;
   float ctheta, pdtheta, dangle;
   bool muoncharge;
   float pz_array[n_entries];
   pz_max = pz_array[0];
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
   mutree->SetBranchAddress("smxst1",&x_sta1);
   mutree->SetBranchAddress("smyst1",&y_sta1);
   mutree->SetBranchAddress("smidx",&idx);
   mutree->SetBranchAddress("smidy",&idy);
   mutree->SetBranchAddress("smcharge",&muoncharge);
   mutree->SetBranchAddress("eventnumber",&eventYield);
   mutree->SetBranchAddress("lvl1_clock_cross",&clock);
   mutree->SetBranchAddress("smst1px",&px_sta1);
   mutree->SetBranchAddress("smst1py",&py_sta1);
   mutree->SetBranchAddress("smst1pz",&pz_sta1);
   mutree->SetBranchAddress("smxst3",&x_sta3);
   mutree->SetBranchAddress("runnumber",&runnumber);
   
      //================================================================== output tree
   float pz_var =0, rapidity_var=0, pt_var=0, phi_var=0, ddg0_var=0, dg0_var=0, pdtheta_var=0, idchi2_var=0, trchi2_var=0, x_F_var=0;
   int muoncharge_var, trhits_var, idhits_var;
   
   bool pt_bin1 = false;
   bool pt_bin2 = false;
   bool pt_bin3 = false;
   bool pt_bin4 = false;
   bool pt_bin5 = false;
   bool pt_bin6 = false;
   
   bool south_cut = false;
   bool north_cut = false;
   bool positive_mu = false;
   bool negative_mu = false;
   bool Spin_Up = false;
   bool Spin_Down = false;
   
   bool xF_positive = false;
   bool xF_negative = false;
   
   
   float ptSpinUp_SouthArm_bin1, ptSpinUp_SouthArm_bin2, ptSpinUp_SouthArm_bin3, ptSpinUp_SouthArm_bin4, ptSpinUp_SouthArm_bin5, ptSpinUp_SouthArm_bin6;
   float ptSpinUp_NorthArm_bin1, ptSpinUp_NorthArm_bin2, ptSpinUp_NorthArm_bin3, ptSpinUp_NorthArm_bin4, ptSpinUp_NorthArm_bin5, ptSpinUp_NorthArm_bin6;
   float ptSpinDown_SouthArm_bin1, ptSpinDown_SouthArm_bin2, ptSpinDown_SouthArm_bin3, ptSpinDown_SouthArm_bin4, ptSpinDown_SouthArm_bin5, ptSpinDown_SouthArm_bin6;
   float ptSpinDown_NorthArm_bin1, ptSpinDown_NorthArm_bin2, ptSpinDown_NorthArm_bin3, ptSpinDown_NorthArm_bin4, ptSpinDown_NorthArm_bin5, ptSpinDown_NorthArm_bin6;
   
   
   TFile *f = new TFile("../fit.root","recreate");
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
   fit->Branch("pt_bin1",&pt_bin1,"pt_bin1/O");
   fit->Branch("pt_bin2",&pt_bin2,"pt_bin2/O");
   fit->Branch("pt_bin3",&pt_bin3,"pt_bin3/O");
   fit->Branch("pt_bin4",&pt_bin4,"pt_bin4/O");
   fit->Branch("pt_bin5",&pt_bin5,"pt_bin5/O");
   fit->Branch("pt_bin6",&pt_bin6,"pt_bin6/O");
   fit->Branch("south_cut",&south_cut,"south_cut/O");
   fit->Branch("north_cut",&north_cut,"north_cut/O");
   fit->Branch("positive_mu",&positive_mu,"positive_mu/O");
   fit->Branch("negative_mu",&negative_mu,"negative_mu/O");
   fit->Branch("Spin_Up",&Spin_Up,"Spin_Up/O");
   fit->Branch("Spin_Down",&Spin_Down,"Spin_Down/O");
   fit->Branch("xF_positive",&xF_positive,"xF_positive/O");
   fit->Branch("xF_negative",&xF_negative,"xF_negative/O");
   
   
   
   
   
   
   
   
   
   for (int ientry = 0; ientry< n_entries; ++ientry){//
      
      mutree->GetEntry(ientry);
      if (ientry%10000==0) cout << "processing event " << ientry << "/" << n_entries <<"\n";
      
#if 1
      //std::cout<<" run number " << runnumber<<"\n";
      pz_array[ientry] = pz;
      if (pz_array[ientry] > pz_max){pz_max = pz_array[ientry];}
     // std::cout<<"pz value : "<< pz <<"\t"<<"pz max  :"<<pz_max<< "\t the value of xF  : "<< pz/pz_max<<" \n";
      x_F_var = pz/pz_max;
      if (x_F_var > 0){xF_positive=1;}
      if (x_F_var < 0){xF_negative=1;}
#endif
      
      
      _h2d->Fill(x_sta3, y_sta3);
      
      
      pt = sqrt(px*px+ py*py);//call.pT(px, py);
     
         
         
      
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
      
      
      
       pt_bin1 = (pt > 1.25 || pt < 1.50);
       pt_bin2 = (pt > 1.50 || pt < 2.00);
       pt_bin3 = (pt > 2.00 || pt < 2.50);
       pt_bin4 = (pt > 2.50 || pt < 3.00);
       pt_bin5 = (pt > 3.00 || pt < 3.50);
       pt_bin6 = (pt > 3.50 || pt < 5.00);
       south_cut = ( trhits > 12 && trchi2 < 10 && idhits > 6 && idchi2 < 5 && ddg0 < 8 && dg0 < 20 && (fabs(rapidity)>1.2 || fabs(rapidity)< 2.0) && pdtheta < 0.2);
       north_cut = ( trhits > 12 && trchi2 < 10 && idhits > 6 && idchi2 < 5 && ddg0 < 8 && dg0 < 10 && (fabs(rapidity)>1.2 || fabs(rapidity)< 2.0) && pdtheta < 0.2);
       negative_mu = (muoncharge ==0);
       positive_mu = (muoncharge ==1);
       if(py > 0){Spin_Up =1;}
       if(py < 0){Spin_Down =1;}
      
      
      
      
      
      
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
      
      
      
      
      
      
      
      
      if( trhits > 12 && trchi2 < 10 && idhits > 6 && idchi2 < 5
         && ddg0 < 8 && dg0 < 20 && (fabs(rapidity)>1.2 || fabs(rapidity)< 2.0)
         && (pt > 1.25 || pt < 5.0) && pdtheta < 0.2){
         
         pT_cut_val = pt;
         mupt_h->Fill(pT_cut_val);
         phi_trk = vec_trk.Phi();
         
      }
      
      
         //==================================================================================================================South
      if (pz<0){
         xF_s_h->Fill(pz/pz_max);
         if(muoncharge ==0){
            
            mupt_h_s_mn->Fill(pT_cut_val);
            muddg0_h_s_mn->Fill(ddg0);
            mudg0_h_s_mn->Fill(dg0);
            chi2_trkzvtx_h_s_mn->Fill(chi2_trk_vtx);
            delta_theta_h_s_mn->Fill(dangle);
            scaled_dtheta_h_s_mn->Fill(pdtheta);
            r_ref_h_s_mn->Fill(r_ref);
            phi_trk_h_s_mn->Fill(phi_trk);
         }//negative muo charge
         
         if(muoncharge==1){
            mupt_h_s_mp->Fill(pT_cut_val);
            muddg0_h_s_mp->Fill(ddg0);
            mudg0_h_s_mp->Fill(dg0);
            chi2_trkzvtx_h_s_mp->Fill(chi2_trk_vtx);
            delta_theta_h_s_mp->Fill(dangle);
            scaled_dtheta_h_s_mp->Fill(pdtheta);
            r_ref_h_s_mp->Fill(r_ref);
            phi_trk_h_s_mp->Fill(phi_trk);
         }
         
      }//pz - south
      
         //=================================================================================================================North
      if (pz > 0 ){
         xF_n_h->Fill(pz/pz_max);
         if(muoncharge==0){
            mupt_h_n_mn->Fill(pT_cut_val);
            muddg0_h_n_mn->Fill(ddg0);
            mudg0_h_n_mn->Fill(dg0);
            chi2_trkzvtx_h_n_mn->Fill(chi2_trk_vtx);
            delta_theta_h_n_mn->Fill(dangle);
            scaled_dtheta_h_n_mn->Fill(pdtheta);
            r_ref_h_n_mn->Fill(r_ref);
            phi_trk_h_n_mn->Fill(phi_trk);
         }
         if(muoncharge==1){
            mupt_h_n_mp->Fill(pT_cut_val);
            muddg0_h_n_mp->Fill(ddg0);
            mudg0_h_n_mp->Fill(dg0);
            chi2_trkzvtx_h_n_mp->Fill(chi2_trk_vtx);
            delta_theta_h_n_mp->Fill(dangle);
            scaled_dtheta_h_n_mp->Fill(pdtheta);
            r_ref_h_n_mp->Fill(r_ref);
            phi_trk_h_n_mp->Fill(phi_trk);
         }
         
         
         
      }//pz
      
      fit->Fill();
   }//Event loop
   fit->Print();
   f->Write();
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
  
   /*
   call.plot_PT(mupt_h_s_mp, mupt_h_s_mn, mupt_h_n_mp, mupt_h_n_mn, "pTDistribution");
   
   call.plot_south_right(delta_theta_h_s_mp, delta_theta_h_s_mn, "delta_theta_south");
   call.plot_south_right(muddg0_h_s_mp, muddg0_h_s_mn, "ddg0_south");
   call.plot_south_right(mudg0_h_s_mp, mudg0_h_s_mn, "dg0_south");
   call.plot_south_right(scaled_dtheta_h_s_mp, scaled_dtheta_h_s_mn, "scaled_dtheta_south");
   call.plot_south_right(r_ref_h_s_mp, r_ref_h_s_mn, "r_ref_south");
   call.plot_south(chi2_trkzvtx_h_s_mp, chi2_trkzvtx_h_s_mn, "chi2_trkzvtx_south");
   call.plot_south(phi_trk_h_s_mp, phi_trk_h_s_mn, "phi_trk_south");
   call.plot_trmom_south(mupt_h_s_mp, mupt_h_s_mn, "mupt_south");
   
   
   
   
   call.plot_north_right(delta_theta_h_n_mp, delta_theta_h_n_mn, "delta_theta_north");
   call.plot_north_right(muddg0_h_n_mp, muddg0_h_n_mn, "ddg0_north");
   call.plot_north_right(mudg0_h_n_mp, mudg0_h_n_mn, "dg0_north");
   call.plot_north_right(scaled_dtheta_h_n_mp, scaled_dtheta_h_n_mn, "scaled_dtheta_north");
   call.plot_north_right(r_ref_h_n_mp, r_ref_h_n_mn, "r_ref_north");
   call.plot_north(chi2_trkzvtx_h_n_mp, chi2_trkzvtx_h_n_mn, "chi2_trkzvtx_north");
   call.plot_north(phi_trk_h_n_mp, phi_trk_h_n_mn, "phi_trk_north");
   call.plot_trmom_north(mupt_h_n_mp, mupt_h_n_mn, "mupt_north");
   
   
  
   */
   TCanvas* c= new TCanvas();
   _h2d->Draw("colz");
   TCanvas* c1= new TCanvas();
   xF_s_h->Draw();
   TCanvas* c2= new TCanvas();
   xF_n_h->Draw();
   
   TCanvas* c3= new TCanvas();
   std::cout<< " number of the entries of the muon histogram"<< mupt_h->GetEntries()<<"\n";
   mupt_h->Draw();
   std::cout<<mupt_h->GetEntries()<<"\n";
}
