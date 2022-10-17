#define asymmetry_analyzer_cxx
#include "asymmetry_analyzer.h"
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
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>


/*
 Author of this code : Muhammad Alibordi
 In a ROOT session, you can do:
 Root > .L asymmetry_analyzer.C
 Root > asymmetry_analyzer t
 Root > t.GetEntry(12); // Fill t data members with entry number 12
 Root > t.Show();       // Show values of entry 12
 Root > t.Show(16);     // Read and show values of entry 16
 Root > t.Loop();       // Loop on all entries
 */


void asymmetry_analyzer::Loop()
{
   
    
   Trig_Efficiency();
   Cuts();
   
	TH1F* XGEN_h = new TH1F("XGEN_h", "XGEN_h", 100, 0.0, 200);

 
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   Float_t s = 200.0 ; //centre of mass energy






   for(unsigned ibin=1; ibin<=17; ibin++){
      for(unsigned iarm =0;iarm <2; iarm++){
         Float_t cut_dg0 = hcut_dg0[iarm][4]->GetBinContent(ibin);
         Float_t  cut_ddg0 = hcut_ddg0[iarm][4]->GetBinContent(ibin);
         Float_t cut_vertex_rad = hcut_vertex_rad[iarm][4]->GetBinContent(ibin);
         Float_t cut_vertex_chi2 = hcut_vertex_chi2[iarm][4]->GetBinContent(ibin);
         
         std::cout<<"ibin: "<< ibin<<"\t"<<"iarm : "<<iarm<<"\t"<<cut_dg0<<"\t"<<cut_ddg0<<"\t"<<cut_vertex_rad<<"\t"<<cut_vertex_chi2<<"\n";
      }
     
   }

     
   
   
   

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(smpz > 0){_arm = 1;}
      if(smpz < 0){_arm = 0;}
      ref_rad = sqrt(smx0*smx0 + smy0*smy0);
      pT = sqrt(smpx*smpx + smpy*smpy);
      //pT = sqrt(smst1px*smst1px + smst1py*smst1py);
      xF = (2*smpz)/(sqrt(s));
      Float_t track_mom_atvtx = sqrt(smpx*smpx + smpy*smpy + smpz*smpz);
      Float_t track_mom_atsta1 = sqrt (smst1px*smst1px + smst1py*smst1py + smst1pz*smst1pz);
      Float_t ctheta=0, dangle=0;
      ctheta = (smpx*smst1px + smpy*smst1py + smpz*smst1pz) / ( track_mom_atvtx*track_mom_atsta1 ) ;
      if ( abs(ctheta) < 1.0 ) {
         dangle = acos (ctheta);
      }
      else{
         dangle = 0.0;
      }
      pdtheta = dangle*0.5*(track_mom_atvtx+track_mom_atsta1 );
      vec_trk.SetXYZ(smx0,smy0,smz0);
      vec_vtx.SetXYZ(evtvtxx, evtvtxy, evtvtxz);
      vec_vtx_z.SetXYZ(0, 0, evtvtxz);
      TVector3  vec_proj =  (vec_trk.Dot(vec_vtx_z)/vec_vtx_z.Mag2())*vec_vtx_z;//  call.vectorProjection(vec_trk,vec_vtx_z );
      Float_t sq_norm_z = vec_proj.Mag2();
      
      if(!(sq_norm_z !=sq_norm_z)){
         chi2_trk_vtx = sq_norm_z;
         
      }
      
      _phist1 = atan2(smyst1,smxst1);
      _phist2 = atan2(smyst2,smxst2);
      _phist3 = atan2(smyst3,smxst3);
      _radst1 = sqrt(smxst1*smxst1 + smyst1*smyst1);
      _radst2 = sqrt(smxst2*smxst2 + smyst2*smyst2);
      _radst3 = sqrt(smxst3*smxst3 + smyst3*smyst3);

      
      float eff_SG3 = 1.0, eff_MUID1D = 1.0;
      float SD_SG3 = hSD_SG3_MUID1DH[_arm]->GetBinContent(hSD_SG3_MUID1DH[_arm]->FindBin(runnumber));
      float SD_MUID1D = hSD_MUID1D[_arm]->GetBinContent(hSD_MUID1D[_arm]->FindBin(runnumber));
       
      if ( smlastgap==4 ){
         eff_SG3 = ftrig_SG3_pT[_arm][smlastgap]->Eval(pT);
         eff_MUID1D = ftrig_MUID1D_pT[_arm]->Eval(pT);
      }
      else{
         eff_SG3 = ftrig_SG3_pT[_arm][smlastgap]->Eval(pT);
      }
      
      //std::cout<<" trigger : smmuid1d == : "<<smmuid1d<< "\t trigger : smmuid1s == :"<< smmuid1s<< "\t trigger : smcharge == :"<< smcharge <<"\n";
      //cout << "Trig eff, arm : " << _arm << ", lastgap : " << smlastgap << ", eff(SG3_MUID1DH) : " << eff_SG3 << ", eff(MUID1D) : " << eff_MUID1D <<"\n";
      
      if(smmuid1d ==1 ){
         weight_by_thrown_pt = 1./eff_MUID1D*(SD_MUID1D+1);
      }
      if ( smlastgap==2 || smlastgap==3 ){
         weight_by_thrown_pt = 1./eff_SG3*(SD_SG3+1);
      }
      cut_pdtheta__ = cut_pdtheta[_arm][smlastgap];
      cut_mutr_chi2__ = cut_mutr_chi2[_arm][smlastgap];
      cut_muid_chi2__ = 6.0;
      
      cut_dg0__ = hcut_dg0[_arm][smlastgap]->GetBinContent(hcut_dg0[_arm][smlastgap]->FindBin(pT));
      cut_ddg0__ = hcut_ddg0[_arm][smlastgap]->GetBinContent(hcut_ddg0[_arm][smlastgap]->FindBin(pT));
      cut_vertex_rad__ = hcut_vertex_rad[_arm][smlastgap]->GetBinContent(hcut_vertex_rad[_arm][smlastgap]->FindBin(pT));
      cut_vertex_chi2__ = hcut_vertex_chi2[_arm][smlastgap]->GetBinContent(hcut_vertex_chi2[_arm][smlastgap]->FindBin(pT));
     
      
      if ( smlastgap==4 ){
         if (
            // _rf_slope>cut_road_slope__ //1
             ref_rad < cut_vertex_rad__ //2
             && smdg0<cut_dg0__ //3
             && smddg0<cut_ddg0__ //4
             && smtrchi2<cut_mutr_chi2__ //5
             && smidchi2<cut_muid_chi2__ //6
             && pdtheta<cut_pdtheta__ //7
             && chi2_trk_vtx<cut_vertex_chi2__ //8
             ){
            avg_pT_varbin[_arm][smcharge]->Fill(pT, pT, weight_by_thrown_pt);
            avg_pT_varbin_gap[_arm][smlastgap][smcharge]->Fill(pT, pT, weight_by_thrown_pt);
            avg_p[_arm][smlastgap][smcharge]->Fill(pT, sqrt(pT*pT+smpz*smpz), weight_by_thrown_pt);
            avg_pT_AN_pz[_arm][smcharge]->Fill(fabs(smpz), pT, weight_by_thrown_pt);
         }
      }
      
     // std::cout<< " arm :\t"<<_arm<<"\t gap \t "<< smlastgap<< "\t vt dg0 \t"<<cut_dg0__<<  "\t"<< smdg0<<"\n";
      //Fill cut variable distributions
      //std::cout<<pT << "\t"<< cut_dg0__ << "\t"<<  cut_ddg0__<< "\t"<< ref_rad<< "\t"<<cut_vertex_rad__<< "\t"<< cut_vertex_chi2__<<"\n";
     // FillDiagnosticHistos(weight_by_thrown_pt);
      
      XGEN_h->Fill(chi2_trk_vtx);
   }//event loop

   
   
   TCanvas * c = new TCanvas();
   XGEN_h->Draw();
}

void asymmetry_analyzer::Trig_Efficiency()
{
   ftrig = new TFile("/Users/md/Documents/Phenix_HF_Analysis/Run15pp200_trig_eff_func_20170301.root","READ");
   cout << "OPEN trigger efficiency file: " << ftrig->GetName() << endl;
   for (int iarm=0; iarm<narm; iarm++){
      for (int igap=2; igap<ngap; igap++){
         
         ftrig_SG3_pT[iarm][igap] = (TF1*)ftrig->Get(Form("ftrig_eff_pT_SG3_MUID1DH_arm%d_gap%d",iarm,igap));
         
         if ( igap==4 ){
            ftrig_MUID1D_pT[iarm] = (TF1*)ftrig->Get(Form("ftrig_eff_pT_MUID1D_arm%d_gap4",iarm));
         }
         
      }//igap
   }//iarm

   fSD = new TFile("/Users/md/Documents/Phenix_HF_Analysis/Run15pp200_prescale.root","read");
   hSD_MUID1D[0] = (TH1F*)fSD->Get("hSD_MUID1D_S");
   hSD_MUID1D[1] = (TH1F*)fSD->Get("hSD_MUID1D_N");
   hSD_SG3_MUID1DH[0] = (TH1F*)fSD->Get("hSD_SG3_MUID1DH_S");
   hSD_SG3_MUID1DH[1] = (TH1F*)fSD->Get("hSD_SG3_MUID1DH_N");
   
   
}

void asymmetry_analyzer::Cuts()
{
   cut_pdtheta[0][2] = 0.40;
   cut_pdtheta[0][3] = 0.40;
   cut_pdtheta[0][4] = 0.25;
   cut_pdtheta[1][2] = 0.40;
   cut_pdtheta[1][3] = 0.40;
   cut_pdtheta[1][4] = 0.25;
   
   cut_mutr_chi2[0][2] = 15;
   cut_mutr_chi2[0][3] = 15;
   cut_mutr_chi2[0][4] = 15;
   cut_mutr_chi2[1][2] = 20;
   cut_mutr_chi2[1][3] = 20;
   cut_mutr_chi2[1][4] = 20;
   
   
   cut_file = new TFile("/Users/md/Documents/Phenix_HF_Analysis/Run15pp200_HFmu_tight_cut_20170404.root","READ");
   for (int iarm=0; iarm<narm; iarm++){
      for (int igap=2; igap<ngap; igap++){
         hcut_dg0[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_dg0_arm%d_gap%d",iarm,igap));
         hcut_ddg0[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_ddg0_arm%d_gap%d",iarm,igap));
         hcut_vertex_rad[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_refrad_arm%d_gap%d",iarm,igap));
         hcut_vertex_chi2[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_vtx_chi2_arm%d_gap%d",iarm,igap));
      }
   }
   
}

void asymmetry_analyzer::BookHistos(){
   char hname[200];
   for (int iarm=0; iarm<narm; iarm++){
      for (int ich=0; ich<ncharge; ich++){
            // variable bin pT profile
         sprintf(hname, "avg_pT_varbin_arm%d_chg%d", iarm, ich);
         avg_pT_varbin[iarm][ich]= new TProfile(hname, "", nptbin, varbin_pTarray);
         
         sprintf(hname, "avg_pT_AN_pz_arm%d_chg%d", iarm, ich);
         avg_pT_AN_pz[iarm][ich]= new TProfile(hname, "", 4, varbin_AN_pzarray);
         
         for (int igap=2; igap<ngap; igap++){
               // variable bin pT profile
            sprintf(hname, "avg_pT_varbin_arm%d_gap%d_chg%d", iarm, igap, ich);
            avg_pT_varbin_gap[iarm][igap][ich]= new TProfile(hname, "", nptbin, varbin_pTarray);
            
            sprintf(hname, "avg_p_varbin_arm%d_gap%d_chg%d", iarm, igap, ich);
            avg_p[iarm][igap][ich]= new TProfile(hname, "", nptbin, varbin_pTarray);
            
               //gap4 variable bin spectra
            sprintf(hname, "n_varbin_arm%d_gap%d_chg%d", iarm, igap, ich);
            N_varbin[iarm][igap][ich] = new TH1F(hname, "", nptbin, varbin_pTarray);
            N_varbin[iarm][igap][ich]->Sumw2();
            
            sprintf(hname, "n_varbin_AN_pT_arm%d_gap%d_chg%d", iarm, igap, ich);
            N_varbin_AN_pT[iarm][igap][ich] = new TH1F(hname, "", 9, varbin_AN_pTarray);
            N_varbin_AN_pT[iarm][igap][ich]->Sumw2();
            
            sprintf(hname, "n_varbin_AN_pz_arm%d_gap%d_chg%d", iarm, igap, ich);
            N_varbin_AN_pz[iarm][igap][ich] = new TH1F(hname, "", 4, varbin_AN_pzarray);
            N_varbin_AN_pz[iarm][igap][ich]->Sumw2();
            
            sprintf(hname, "n_varbin_AN_pz_MUID1D_arm%d_gap%d_chg%d", iarm, igap, ich);
            N_varbin_AN_pz_MUID1D[iarm][igap][ich] = new TH1F(hname, "", 4, varbin_AN_pzarray);
            N_varbin_AN_pz_MUID1D[iarm][igap][ich]->Sumw2();
            
            sprintf(hname, "n_varbin_AN_pz_SG3MUID1DH_arm%d_gap%d_chg%d", iarm, igap, ich);
            N_varbin_AN_pz_SG3MUID1DH[iarm][igap][ich] = new TH1F(hname, "", 4, varbin_AN_pzarray);
            N_varbin_AN_pz_SG3MUID1DH[iarm][igap][ich]->Sumw2();
            
            sprintf(hname, "n_varbin_fake_arm%d_gap%d_chg%d", iarm, igap, ich);
            N_varbin_fake[iarm][igap][ich] = new TH1F(hname, "", nptbin, varbin_pTarray);
            N_varbin_fake[iarm][igap][ich]->Sumw2();
            
            sprintf(hname, "n_pT_pz_arm%d_gap%d_chg%d", iarm, igap, ich);
            N_pT_pz[iarm][igap][ich] = new TH2F(hname, "", 30, 0, 15, 30, 0, 45);
            N_pT_pz[iarm][igap][ich]->Sumw2();
            }//GAP
         }//CHARGE
      for (int jgap=2; jgap<ngap; jgap++){
            //DG0
         sprintf(hname, "dg0_arm%d_gap%d", iarm, jgap);
         DG0[iarm][jgap] = new TH2F(hname,"",nptbin,varbin_pTarray,80,0,40);
         DG0[iarm][jgap]->Sumw2();
            //DDG0
         sprintf(hname, "ddg0_arm%d_gap%d", iarm, jgap);
         DDG0[iarm][jgap] = new TH2F(hname,"",nptbin,varbin_pTarray,60,0,30);
         DDG0[iarm][jgap]->Sumw2();
            //REFRAD
         sprintf(hname, "refrad_arm%d_gap%d", iarm, jgap);
         REFRAD[iarm][jgap] = new TH2F(hname,"",nptbin,varbin_pTarray,250,0.,250.0);
         REFRAD[iarm][jgap]->Sumw2();
            //MUTR CHI2
         sprintf(hname, "mutr_chi2_arm%d_gap%d", iarm, jgap);
         MUTR_CHI2[iarm][jgap] = new TH2F(hname,"",nptbin,varbin_pTarray,100,0,25);
         MUTR_CHI2[iarm][jgap]->Sumw2();
            //VTX CHI2
         sprintf(hname, "vtx_chi2_arm%d_gap%d",iarm, jgap);
         VTX_CHI2[iarm][jgap] = new TH2F(hname,"",nptbin,varbin_pTarray,400,0,20);
         VTX_CHI2[iarm][jgap]->Sumw2();
         
            //PDTHETA
         sprintf(hname, "pdtheta_arm%d_gap%d", iarm, jgap);
         PDTHETA[iarm][jgap] = new TH2F(hname,hname,nptbin,varbin_pTarray,50,0.0,1.0);
         PDTHETA[iarm][jgap]->Sumw2();
            
        
      }//GAP-2
      }//ARM
}
