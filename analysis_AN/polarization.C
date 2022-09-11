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

using  namespace std;
using namespace RooFit;

void polarization(){
   
   TProfile *fillvspol_b = new TProfile("fillvspol_b","fillvspol;fill;polarization",100,18650,19000,-1, 1 );
   TProfile *fillvspol_y = new TProfile("fillvspol_y","fillvspol;fill;polarization",100,18650,19000,-1, 1 );
   TProfile *runvspol_b = new TProfile("runvspol_b","runvspol;run;polarization",100,4.21E5,4.34E5,-1, 1 );
   TProfile *runvspol_y = new TProfile("runvspol_y","runvspol;run;polarization",100,4.21E5,4.34E5,-1, 1 );
   TProfile *bunchvsscaler = new TProfile("bunchvsscaler","bunchvsscaler;bunch;scaler",120,0,120, 0, 4E9 );
   
   
   int run, fill, bunch;
   Long64_t scalerA[120];
   float polblue, polyellow;
   
   TFile *mufile = new TFile("../spinDB_test.root");
   TTree *mutree = (TTree*)mufile->Get("T");
   int n_entries = mutree->GetEntries();
   std::cout<<" Number of entries for the moment : \t"<< n_entries <<"\n";
   
   mutree->SetBranchAddress("runnumber",&run);
   mutree->SetBranchAddress("fillnumber",&fill);
   mutree->SetBranchAddress("polblue",&polblue);
   mutree->SetBranchAddress("polyellow",&polyellow);
   mutree->SetBranchAddress("scalerA",&scalerA[0]);
   for (int ientry = 0; ientry<=n_entries; ++ientry){
      
      mutree->GetEntry(ientry);
      if (ientry%100==0) cout << "processing event " << ientry << "/" << n_entries <<"\n";
      fillvspol_b->Fill(fill, polblue,1);
      fillvspol_y->Fill(fill, polyellow,1);
      runvspol_y->Fill(run, polyellow,1);
      runvspol_b->Fill(run, polblue,1);
      for(unsigned k =0 ; k<=120; ++k){
         bunch = k;
         Long64_t scaler_bbcvtxcut = scalerA[k];
         std::cout<<" scaler value   :"<<k<<"that was the value of bunch    : "<<scaler_bbcvtxcut<<"\n";
         bunchvsscaler->Fill(k, scaler_bbcvtxcut,1);
      }
      
      
      
   }
   gStyle->SetOptStat(0);
   
   TCanvas *c = new TCanvas();
   fillvspol_b->SetLineColor(kBlue);
   fillvspol_b->SetMarkerColor(kBlue);
   fillvspol_b->SetMarkerStyle(8);
   fillvspol_b->Draw();
   fillvspol_y->SetLineColor(kYellow+1);
   fillvspol_y->SetMarkerColor(kYellow+1);
   fillvspol_y->SetMarkerStyle(8);
   fillvspol_y->Draw("SAME");
   auto legend = new TLegend(0.6,0.6,0.8,0.8);
   legend->SetHeader("Fill","L");
   legend->AddEntry(fillvspol_b,"Blue","l");
   legend->AddEntry(fillvspol_y,"Yellow","l");
   legend->Draw();
   
   
   TCanvas *c1 = new TCanvas();
   runvspol_b->SetLineColor(kBlue);
   runvspol_b->SetMarkerColor(kBlue);
   runvspol_b->SetMarkerStyle(8);
   runvspol_b->Draw();
   runvspol_y->SetLineColor(kYellow+1);
   runvspol_y->SetMarkerColor(kYellow+1);
   runvspol_y->SetMarkerStyle(8);
   runvspol_y->Draw("SAME");
   auto legenda = new TLegend(0.6,0.6,0.8,0.8);
   legenda->SetHeader("Run","L");
   legenda->AddEntry(runvspol_b,"Blue","l");
   legenda->AddEntry(runvspol_y,"Yellow","l");
   legenda->Draw();
   
   TCanvas *c2 = new TCanvas("c2","c2", 1400, 400);
   bunchvsscaler->SetLineColor(kBlue);
   bunchvsscaler->SetFillColor(kOrange+2);
   bunchvsscaler->SetBarWidth(0.4);
   bunchvsscaler->SetBarOffset(0.1);
   bunchvsscaler->SetFillStyle(3005);
   bunchvsscaler->DrawNormalized();
   auto legendb = new TLegend(0.7,0.75,0.87,0.95);
   legendb->SetHeader("scaler=f(bunch)","L");
   legendb->AddEntry(bunchvsscaler,"scaler","l");
   legendb->Draw();
   
}




