#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include "TGraphAsymmErrors.h"
#include "TVirtualFFT.h"
#include "TBinomialEfficiencyFitter.h"
#include "TVectorF.h"
#include "TPaveText.h"
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TEventList.h"
#include "Riostream.h"
#include "string.h"
#include "TList.h"
#include "TDirectory.h"
#include "TCut.h"
#include "TChain.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "RooDataHist.h"
#include <fstream>
using namespace RooFit ;
using namespace std;

void Time_Resolution(){
   ofstream ptr;
   ptr.open("table.dat",ios::out);
   RooDataHist *bindata[16];
   RooGaussian *gauss[16];
   RooFitResult* fitRes[16]; TH2* hcorr[16];
   RooPlot* DelT[16];
   TFile * lfile = new TFile("../timediff_tree_00023342.root");//ltest-00023342.root");
   TTree * ltree = (TTree*)lfile->Get("t");
   Int_t nentries = ltree->GetEntries();
   std::cout<<" Number of events in the tree "<< nentries<<"\n";
   ltree->Print();
   Double_t amplitude_L03C_allpixel_for2D[24][24];
   ltree->SetBranchAddress("amplitude_L03C_allpixel_for2D",&amplitude_L03C_allpixel_for2D[0][0]);
   TH1D * time_r[16];
   Double_t time_diff, amp, x , y, cluster_x, cluster_y ;
   Double_t max_coor_cut,min_coor_cut;
   Double_t amplitude;
   TH2F *hxvsy = new TH2F("hxvsy", "Cluster position x vs cluster position;Cluster position x ;Cluster position y", 24, -50.0,50,24, -50.0,50  );
    auto h2 = new TH2D("h2", "", 24, 0.0, 24, 24, 0.0, 24);
   for (unsigned i =0; i<16; ++i){
      time_r[i]= new TH1D(Form("mhistogram_%d", i), Form("diffrence_time_[%d-(%d-1)-60]", i,i), 100, -500, 500  );
   }
   
   TLeaf* l_tt = (TLeaf*) ltree->GetLeaf("tt");
   TLeaf* l_aa = (TLeaf*) ltree->GetLeaf("aa");
   TLeaf* l_tx = (TLeaf*) ltree->GetLeaf("trk_tx");
   TLeaf* l_ty = (TLeaf*) ltree->GetLeaf("trk_ty");
   TLeaf* l_cluster_L03C_X = (TLeaf*) ltree->GetLeaf("cluster_L03C_X");
   TLeaf* l_cluster_L03C_Y = (TLeaf*) ltree->GetLeaf("cluster_L03C_Y");
   TLeaf* l_amplitude_L03C_allpixel_for2D = (TLeaf*) ltree->GetLeaf("amplitude_L03C_allpixel_for2D");
   
   amplitude_L03C_allpixel_for2D
  
   for(int i = 0; i < ltree->GetEntries(); i++){
      ltree->GetEntry(i);
      x = l_tx->GetValue();
      y = l_ty->GetValue();
      cluster_x = l_cluster_L03C_X->GetValue();
      cluster_y = l_cluster_L03C_Y->GetValue();
      hxvsy->Fill(cluster_x, cluster_y);
      
      max_coor_cut = x + y + 12;
      min_coor_cut = x - y + 12;
      if(abs(max_coor_cut) > 5.0) continue;
      if(abs(min_coor_cut) > 10.0) continue;
      std::cout <<" the value of x and y "<< x << "\t" << y<<"\n";
      for(unsigned j =0; j< 16; ++j){
        // std::cout<<l_tt->GetValue(j)<<"\n";
         if(abs(l_aa->GetValue(j)/l_aa->GetValue(j-1)-1.0) > 0.5)continue;
         if(l_aa->GetValue(j) < 50 && l_aa->GetValue(j) > 250)continue;
         if(l_aa->GetValue(j-1) < 50 && l_aa->GetValue(j-1) > 250)continue;
         time_diff = l_tt->GetValue(j)-l_tt->GetValue(j-1)-60;
         time_r[j]->Fill(time_diff);
         
         
      }
      
      for(unsigned k =0; k<24; ++k){
         for(unsigned l =0 ; l<24; ++l){
            Double_t amp = amplitude_L03C_allpixel_for2D[k][l];
            h2->SetBinContent(k, l, amp);
         }
      }
      
     }//event loop
      
      
   
   RooRealVar *timediff = new RooRealVar("timediff", "timediff", -500, 500);
   RooRealVar *mean= new RooRealVar("mean", "mean", -100, 100);
   RooRealVar *sigma = new RooRealVar("sigma", "sigma",0.,100);
   TCanvas *c = new TCanvas("c", "c",0,0,600,400);
   for(unsigned j =0; j< 16; ++j){
      
      bindata[j] = new RooDataHist("bindata", "bindata", *timediff, Import(*time_r[j]));
      gauss [j] = new RooGaussian("gauss","gauss",*timediff,*mean,*sigma) ;
      fitRes[j] = gauss [j]->fitTo(*bindata[j],Save(),NumCPU(8));
      fitRes[j]->Print("v");
      DelT[j] = timediff->frame(Title("#DeltaT (= t_{j}-t_{j-1}) (ps) "),Bins(40));
      bindata[j]->plotOn(DelT[j],DataError(RooAbsData::SumW2));
      gauss [j]->plotOn(DelT[j]) ;
      gauss [j]->paramOn(DelT[j]);
      DelT[j]->Draw();
       c->Update();
       c->SaveAs(Form("/Users/md/Documents/Phenix_HF_Analysis/lappd_event/time/fittedT_%d.pdf",j));
      std::cout<< time_r[j]->GetRMS() << "\t +-" <<time_r[j]->GetRMSError()<<"\n";
      ptr << "\u0394" <<"T"<<":["<<j<<"-("<<j<<"-1)]"<<"\t"<< time_r[j]->GetRMS()<<"\t"<<"\t"<<"\t" << "\u0394" <<"TErr:     " << time_r[j]->GetRMSError() <<"\n";
      
      
   }
   TCanvas *c1 = new TCanvas();
   for(unsigned j =0; j< 16; ++j){
   
      time_r[j]->Draw();
         //time_r[j]->Fit("gaus", "L");
      c1->SaveAs(Form("/Users/md/Documents/Phenix_HF_Analysis/lappd_event/time/time-diff_%d.png",j));
   }
  
   TCanvas *c8 = new TCanvas();h2->Draw("COLZ"); c8->SaveAs("/Users/md/Documents/Phenix_HF_Analysis/lappd_event/time/xvsy.pdf", "pdf");
   
   


   
}
