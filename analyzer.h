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

using  namespace std;
using namespace RooFit;

class definition
{
private:
   
public:
   void plot(TH1F * hist1 , TH1F * hist2, string histname){
      TCanvas *c = new TCanvas();
      hist1->SetLineColor(kRed);
      hist1->SetMarkerStyle(8);
      hist1->Draw();
      hist2->SetLineColor(kBlack);
      hist2->SetMarkerStyle(8);
      hist2->Draw("SAME");
      auto legend = new TLegend(0.1,0.7,0.48,0.9);
      legend->AddEntry(hist1,"South","l");
      legend->AddEntry(hist2,"North","l");
      legend->Draw();
      c->SaveAs(("plot/"+histname+".pdf").c_str());
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
      Theta = TMath::ACos(cosineTheta);
      return Theta;
   }
   TVector3 vectorProjection(TVector3 vec_A, TVector3 vec_B){
      TVector3 projection;
      projection = (vec_A.Dot(vec_B)/vec_B.Mag2())*vec_B;
      return projection;
   }
   
   
   
};
#endif // PHYSICS_OBSERVALES_H
