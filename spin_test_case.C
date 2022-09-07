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

int muonmultiplicity ;
int clock;
Long64_t scaler[120];
int bunch;

void spin_test_case(string data, string spindb, string output){
   TProfile *bunchvsscaler = new TProfile("bunchvsscaler","bunchvsscaler;bunch;scaler",100,0,120, 0, 4E9 );
   TProfile *mmulvsclock = new TProfile("mmulvsclock","mmulvsclock;bunch;muon multiplicity",120,0,120, 0, 40E6 );
   
   std::cout<<" do ls /direct/phenix+u/alibordi/copy_datacode/ spin files, just the file name--> spindB \n"
   <<" spinDB_test_421815.root, spinDB_test_421815.root  spinDB_test_424227.root  spinDB_test_424896.root  \n"
   << " spinDB_test_421989.root  spinDB_test_424885.root  spinDB_test_431963.root   \n";
   
   std::cout<<" do ls /direct/phenix+u/alibordi/hf_outputs/ spin files, just the file name--> data \n"
            <<" analysis_preliminary_421815.root  \n"
            <<" analysis_preliminary_421989.root   \n"
            <<"  analysis_preliminary_424227.root \n"
            <<"  analysis_preliminary_424885.root  \n"
            <<"  analysis_preliminary_429896.root  \n"
            <<"  analysis_preliminary_431963.root  \n"<<"\n";
            
   
   TFile *mufile = new TFile(("/direct/phenix+u/alibordi/copy_datacode/"+spindb+".root").c_str());
   TTree *mutree = (TTree*)mufile->Get("T");
   int n_entries = mutree->GetEntries();
   std::cout<<" Number of entries for the moment : \t"<< n_entries <<"\n";
   
   int run, fill, bunch;
   float polblue, polyellow;
  int trhits, idhits, eventyield, clock; 
   
   mutree->SetBranchAddress("runnumber",&run);
   mutree->SetBranchAddress("fillnumber",&fill);
   mutree->SetBranchAddress("polblue",&polblue);
   mutree->SetBranchAddress("polyellow",&polyellow);
   mutree->SetBranchAddress("scalerA",&scaler[0]);
   for (int ientry = 0; ientry<=n_entries; ++ientry){
      
      mutree->GetEntry(ientry);
      if (ientry%100==0) cout << "processing event " << ientry << "/" << n_entries <<"\n";
     
      for(unsigned k =0 ; k<120; ++k){
         bunch = k;
         Long64_t scaler_bbcvtxcut = scaler[k];
         //std::cout<<" scaler value   :" <<scaler_bbcvtxcut<<"\n";
         bunchvsscaler->Fill(k, scaler_bbcvtxcut,1);
      }
      
      
      
   }
   
   
   std::cout<<"======================================data tree output======================"<<"\n";
   
   TFile *mufile_data = new TFile(("/direct/phenix+u/alibordi/hf_outputs/"+data+".root").c_str());
   TTree *mutree_data = (TTree*)mufile_data->Get("analysis");
   int n_entries_data = mutree_data->GetEntries();
   std::cout<<" Number of entries for the moment : \t"<< n_entries_data <<"\n";
   
   mutree_data->SetBranchAddress("eventnumber",&eventyield);
   mutree_data->SetBranchAddress("lvl1_clock_cross",&clock);
   
   for (int ientry = 0; ientry< n_entries_data; ++ientry){//
      
      mutree_data->GetEntry(ientry);
      if (ientry%100==0) cout << "processing event " << ientry << "/" << n_entries_data <<"\n";
      //std::cout<<muonmultiplicity <<"\t lolololol "<< (clock +5)%120<<"\n";
      mmulvsclock->Fill((clock +5)%120, eventyield,1);

/*for (unsigned i =0; i<=clock_h->GetNbinsX(); ++i ){
 *          int clock_value = clock_h->GetBinContent(i);
 *                  // mmulvsclock->Fill(i+5, muonmultiplicity,1);
 *                        }*/
   

}
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0); 
   
   TCanvas *c2 = new TCanvas();
   bunchvsscaler->SetLineColor(kRed);
   mmulvsclock->SetLineColor(kBlue);
   bunchvsscaler->DrawNormalized();
   mmulvsclock->DrawNormalized("SAME");
   TLegend* legendb = new TLegend(0.15,0.7,0.35,0.9);
   //legendb->SetHeader(("run#-"+output).c_str(),"C");
   legendb->AddEntry(mmulvsclock,"# of #mu","l");
   legendb->AddEntry(bunchvsscaler,"scaler","l");
   legendb->Draw();
   c2->SaveAs(("/direct/phenix+u/alibordi/copy_datacode/"+output+".png").c_str());
   
   
   
   
   
   
}

