#ifndef ANALYZER_H
#define ANALYZER_H


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
#include "TStyle.h"
using  namespace std;
//using namespace RooFit;

class definition
{
private:
   
public:
   
   float energy(float pz, float rapidity){
      float energy ;
      energy = pz*TMath::TanH(rapidity);
      return energy;
   }
   
   float pT( float px, float py){
      return sqrt(px*px+ py*py);
   }
   float costheta(TVector3 vec_A, TVector3 vec_B){
      float cosineTheta;
      cosineTheta =  (vec_A.Dot( vec_B))/((vec_A.Mag())*(vec_B.Mag()));
      return cosineTheta;
   }
   float cosinverse(TVector3 vec_A, TVector3 vec_B){
      float cosineTheta, Theta;
      cosineTheta =  (vec_A.Dot( vec_B))/((vec_A.Mag())*(vec_B.Mag()));
      if(abs(cosineTheta)<1.0){
      Theta = TMath::ACos(cosineTheta);
      }
      else{
         Theta =0.0;
      }
      return Theta;
   }

   TVector3 vectorProjection(TVector3 vec_A, TVector3 vec_B){	
      TVector3 projection;
      projection = (vec_A.Dot(vec_B)/vec_B.Mag2())*vec_B;
      return projection;
   }


   void plot(TH1F * hist1 , TH1F * hist2, string histname){
      TCanvas *c = new TCanvas();
      hist1->SetLineColor(kRed);
      hist1->SetMarkerStyle(8);
      hist1->Draw();
      hist2->SetLineColor(kBlack);
      hist2->SetMarkerStyle(8);
      hist2->Draw("SAME");
      TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
      legend->AddEntry(hist1,"South","l");
      legend->AddEntry(hist2,"North","l");
      legend->Draw();
      c->SaveAs(("plot/"+histname+".pdf").c_str());
   }


   
   void plot_north(TH1F * hist1 , TH1F * hist2, string histname){
      TCanvas *c = new TCanvas();
      hist1->SetLineColor(kRed);
      hist1->SetMarkerStyle(8);
      hist1->Draw();
      hist2->SetLineColor(kBlack);
      hist2->SetMarkerStyle(8);
      hist2->Draw("SAME");
      TLegend* legend = new TLegend(0.15,0.6,0.35,0.8);
      legend->SetHeader("North","L");
      legend->AddEntry(hist1,"#mu^{+}","l");
      legend->AddEntry(hist2,"#mu^{-}","l");
      legend->Draw();
      c->SaveAs(("../plot/"+histname+".pdf").c_str());
   }
   void plot_north_right(TH1F * hist1 , TH1F * hist2, string histname){
      TCanvas *c = new TCanvas();
      hist1->SetLineColor(kRed);
      hist1->SetMarkerStyle(8);
      hist1->Draw();
      hist2->SetLineColor(kBlack);
      hist2->SetMarkerStyle(8);
      hist2->Draw("SAME");
      TLegend* legend = new TLegend(0.6,0.55,0.8,0.85);
      legend->SetHeader("North","L");
      legend->AddEntry(hist1,"#mu^{+}","l");
      legend->AddEntry(hist2,"#mu^{-}","l");
      legend->Draw();
      c->SaveAs(("../plot/"+histname+".pdf").c_str());
   }
   void plot_south(TH1F * hist1 , TH1F * hist2, string histname){
      TCanvas *c = new TCanvas();
      hist1->SetLineColor(kRed);
      hist1->SetMarkerStyle(8);
      hist1->Draw();
      hist2->SetLineColor(kBlack);
      hist2->SetMarkerStyle(8);
      hist2->Draw("SAME");
      TLegend* legend = new TLegend(0.15,0.6,0.35,0.8);
      legend->SetHeader("South","L");
      legend->AddEntry(hist1,"#mu^{+}","l");
      legend->AddEntry(hist2,"#mu^{-}","l");
      legend->Draw();
      c->SaveAs(("../plot/"+histname+".pdf").c_str());
   }
   void plot_south_right(TH1F * hist1 , TH1F * hist2, string histname){
      TCanvas *c = new TCanvas();
      hist1->SetLineColor(kRed);
      hist1->SetMarkerStyle(8);
      hist1->Draw();
      hist2->SetLineColor(kBlack);
      hist2->SetMarkerStyle(8);
      hist2->Draw("SAME");
      TLegend* legend = new TLegend(0.6,0.55,0.8,0.85);
      legend->SetHeader("South","L");
      legend->AddEntry(hist1,"#mu^{+}","l");
      legend->AddEntry(hist2,"#mu^{-}","l");
      legend->Draw();
      c->SaveAs(("../plot/"+histname+".pdf").c_str());
   }
   
   void plot_PT(TH1F * hist1 , TH1F * hist2, TH1F * hist3, TH1F * hist4, string histname){
      TCanvas *c = new TCanvas();
      hist1->SetLineColor(kRed);
      hist1->SetMarkerStyle(8);
      hist1->Draw();
      hist2->SetLineColor(kBlack);
      hist2->SetMarkerStyle(8);
      hist2->Draw("SAME");
      hist3->SetLineColor(kBlue);
      hist3->SetMarkerStyle(8);
      hist3->Draw("SAME");
      hist4->SetLineColor(kGreen+2);
      hist4->SetMarkerStyle(8);
      hist4->Draw("SAME");
      TLegend* legend = new TLegend(0.6,0.55,0.8,0.85);
      legend->SetHeader("All","C");
      legend->AddEntry(hist1,"South #mu^{+}","l");
      legend->AddEntry(hist2,"South #mu^{-}","l");
      legend->AddEntry(hist3,"North #mu^{+}","l");
      legend->AddEntry(hist4,"North #mu^{-}","l");
      legend->Draw();
      c->SaveAs(("../plot/"+histname+".pdf").c_str());
   }
   
   void plot_trmom_south(TH1F * hist1 , TH1F * hist2, string histname){
      TCanvas *c = new TCanvas();
      hist1->SetLineColor(kRed);
      hist1->SetMarkerStyle(8);
      hist1->Draw();
      hist2->SetLineColor(kBlack);
      hist2->SetMarkerStyle(8);
      hist2->Draw("SAME");
      TLegend* legend = new TLegend(0.6,0.6,0.8,0.8);
      legend->SetHeader("South","L");
      legend->AddEntry(hist1,"#mu^{+}","l");
      legend->AddEntry(hist2,"#mu^{-}","l");
      legend->Draw();
      c->SaveAs(("../plot/"+histname+".pdf").c_str());
   }
   
   void plot_trmom_north(TH1F * hist1 , TH1F * hist2, string histname){
      TCanvas *c = new TCanvas();
      hist1->SetLineColor(kRed);
      hist1->SetMarkerStyle(8);
      hist1->Draw();
      hist2->SetLineColor(kBlack);
      hist2->SetMarkerStyle(8);
      hist2->Draw("SAME");
      TLegend* legend = new TLegend(0.6,0.6,0.8,0.8);
      legend->SetHeader("North","L");
      legend->AddEntry(hist1,"#mu^{+}","l");
      legend->AddEntry(hist2,"#mu^{-}","l");
      legend->Draw();
      c->SaveAs(("../plot/"+histname+".pdf").c_str());
   }
   
   

  
   
   int findMax(int arr[], int low, int high)
   {
   
   
   if (high == low)
      return arr[low];
   int mid = low + (high - low) / 2;
   if(mid==0 && arr[mid]>arr[mid+1])
      {
      return arr[mid];
      }
   
   if (mid < high && arr[mid + 1] < arr[mid] && mid>0 && arr[mid]>arr[mid-1]) {
      return arr[mid];
   }
   if (arr[low] > arr[mid]) {
      return findMax(arr, low, mid - 1);
   }
   else {
      return findMax(arr, mid + 1, high);
   }
   }
   
   /*
    theta_vtx_h_s_mn->Fill(theta_vtx);
    theta_mutr_h_s_mn->Fill(theta_mutr);
    //costheta_vtx_h_s_mn->Fill(costheta_vtx);
    //costheta_mutr_h_s_mn->Fill(costheta_mutr);
    theta_vtx_h_s_mp->Fill(theta_vtx);
    theta_mutr_h_s_mp->Fill(theta_mutr);
    //costheta_vtx_h_s_mp->Fill(costheta_vtx);
    //costheta_mutr_h_s_mp->Fill(costheta_mutr);
    theta_vtx_h_n_mn->Fill(theta_vtx);
    theta_mutr_h_n_mn->Fill(theta_mutr);
    //costheta_vtx_h_n_mn->Fill(costheta_vtx);
    //costheta_mutr_h_n_mn->Fill(costheta_mutr);
    theta_vtx_h_n_mp->Fill(theta_vtx);
    theta_mutr_h_n_mp->Fill(theta_mutr);
    //costheta_vtx_h_n_mp->Fill(costheta_vtx);
    //costheta_mutr_h_n_mp->Fill(costheta_mutr);
    call.plot_north(theta_vtx_h_n_mp, theta_vtx_h_n_mn, "theta_vtx_north");
    call.plot_north(theta_mutr_h_n_mp, theta_mutr_h_n_mn, "theta_mutr_north");
    call.plot_south(theta_mutr_h_s_mp, theta_mutr_h_s_mn, "theta_mutr_south");
    call.plot_south(theta_vtx_h_s_mp, theta_vtx_h_s_mn, "theta_vtx_south");
    call.plot_south(costheta_vtx_h_s_mp, costheta_vtx_h_s_mn, "costheta_vtx_south");
    call.plot_south(costheta_mutr_h_s_mp, costheta_mutr_h_s_mn, "costheta_mutr_south");
    call.plot_north(costheta_vtx_h_n_mp, costheta_vtx_h_n_mn, "costheta_vtx_north");
    call.plot_north(costheta_mutr_h_n_mp, costheta_mutr_h_n_mn, "costheta_mutr_north");
    */
   
};
#endif // PHYSICS_OBSERVALES_H
