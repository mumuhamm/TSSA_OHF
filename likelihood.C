/// \author  - Muhammad Alibordi
// Test of RooKeyPDF has ability to discriminate the background and

#include <TLatex.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include <TString.h>
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooGaussModel.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TCut.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooHist.h"
#include "RooGenericPdf.h"
#include "RooTruthModel.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "RooEffProd.h"
using namespace std;
using namespace RooFit;


void likelihood(){
   
   Int_t nbins = 100;
   TChain* chain_data = new TChain("fit");
   chain_data->Add("../fit.root");
   Int_t nevt = (int)chain_data->GetEntries();
   std::cout<<"Number of total events  : "<<nevt<<"\n";
   
   TChain* chain_data_spin = new TChain("T");
   chain_data_spin->Add("../spinDB_test.root");
   Int_t nevt_spin = (int)chain_data_spin->GetEntries();
   std::cout<<"Number of total events in the spin tree   : "<<nevt_spin<<"\n";
   //Creating a data set which we are going to fit with the variables defined in the PHYSICAL REVIEW D 95, 112001 (2017)
   
   RooRealVar pt_var("pt_var", "p_{T} GeV/c",0.0,7.0);
   RooRealVar phi_var("phi_var", "#phi (rad)",-TMath::Pi(), TMath::Pi());
   RooRealVar pdtheta_var("pdtheta_var", "p#bullet#delta#theta (rad-GeV/c)",0.0, 0.4);
   RooRealVar rapidity_var("rapidity_var", "y (rad)",-TMath::Pi(), TMath::Pi());
   RooRealVar ddg0_var("ddg0_var", "DDG0",0,20.0);
   RooRealVar dg0_var("dg0_var", "DG0",0,50.0);
   RooRealVar pz_var("pz_var", "p_{z} (GeV/c)",-40, 40);
   RooRealVar idchi2_var("idchi2_var", "#chi^{2}_{ID}",0,20.0);
   RooRealVar trchi2_var("trchi2_var", "#chi^{2}_{ID}",0,30.0);
   RooRealVar muoncharge_var("muoncharge_var", "Q_{#mu}",0,2);
   RooRealVar trhits_var("trhits_var", "hits_{TR}",0,30000);
   RooRealVar idhits_var("idhits_var", "hits_{ID}",0,1000);
   
   RooRealVar polyellow("polyellow", "P",0.4,0.8);
   /*RooCategory phi_pol("phi_pol","phi_pol");
   phi_pol.defineType("B",-TMath::Pi()*0.5);
   phi_pol.defineType("Y",TMath::Pi()*0.5);
   phi_pol.defineType("NONE",0);
   */
   RooArgSet s( pz_var, rapidity_var, pt_var, phi_var, ddg0_var, dg0_var, pdtheta_var);
   s.add(RooArgSet(muoncharge_var, trhits_var, idhits_var, idchi2_var, trchi2_var,polyellow/*,phi_pol*/));

   RooDataSet Spindata("Spindata", "Spindata", s, Import(*chain_data_spin));
   
   RooDataSet data("data", "data", s, Import(*chain_data));
   data.append(Spindata);
   data.Print("v");
   
   TCut c1 = "TMath::Abs(rapidity_var)>1.2 || TMath::Abs(rapidity_var)< 2.0";
   TCut c2 = "ddg0_var < 8 && dg0_var < 20";
   TCut c3 = "trhits_var > 12 && trchi2_var < 10";
   TCut c4 = "idhits_var > 6 && idchi2_var < 5";
   TCut c5 = "pdtheta_var < 0.2 && pz_var < 0";
   TCut c6 = "muoncharge_var == 0";
   TCut ptbin_1 = "pt_var > 1.25 || pt_var < 1.50";
   TCut ptbin_2 = "pt_var > 1.50 || pt_var < 2.00";
   TCut ptbin_3 = "pt_var > 2.00 || pt_var < 2.50";
   TCut ptbin_4 = "pt_var > 2.50 || pt_var < 3.00";
   TCut ptbin_5 = "pt_var > 3.00 || pt_var < 3.50";
   TCut ptbin_6 = "pt_var > 3.50 || pt_var < 5.00";
   
   
   TCut signalRegion = c1 && c2 && c3 && c4 && c5 && c6 && ptbin_1 ;
   
   
   RooDataSet *data_SigReg = (RooDataSet*)data.reduce(signalRegion);
   std::cout<< " number of entries in the roodata set   :"<<    data_SigReg->sumEntries()<<"\n";
   
   
   RooRealVar ns("ns","number signal events",0.0,data_SigReg->sumEntries());
   
   RooRealVar AN("AN", "AN",-1,1);
   RooRealVar phi_pol("phi_pol", "phi_pol",-TMath::Pi()*0.5, TMath::Pi()*0.5);
   RooFormulaVar sinpart("sinpart", "sinpart", "sin(TMath::Pi()*0.5-@0)",RooArgSet(phi_var));
   //RooFormulaVar sinpart("sinpart", "sinpart", "sin(@0-@1)",RooArgSet(phi_pol,phi_var));
   //RooFormulaVar sinpart("sinpart", "sinpart", "sin(@0-@1)*sin(@0-@1)",RooArgSet(phi_pol,phi_var));
   RooFormulaVar polpart("polpart", "polpart", "@0*@1",RooArgSet(polyellow,AN));
   RooGenericPdf formula("formula", "1+@0*@1", RooArgSet(sinpart,polpart));
   RooExtendPdf model("model","extended signal p.d.f",formula, ns,"signalRange") ;
   RooFitResult* AsymmResults = model.fitTo(*data_SigReg, Extended(), Save(kTRUE));
   AsymmResults->Print("v");
   gStyle->SetOptStat(0) ;
   gStyle->SetPalette(1) ;
   TH2* hcorr = AsymmResults->correlationHist() ;
   TCanvas* c = new TCanvas("Correlation Matrix","Correlation Matrix",600,600) ;
   hcorr->GetYaxis()->SetTitleOffset(1.4) ; hcorr->Draw("colz");
   
   RooPlot* phi_fit = phi_var.frame(Title("#phi (rad)"),Bins(50));
   data_SigReg->plotOn(phi_fit,DataError(RooAbsData::SumW2));
   model.plotOn(phi_fit) ;
   model.paramOn(phi_fit);
   model.plotOn(phi_fit, LineColor(kBlack), LineWidth(1));
   Double_t chisquare_costheta = phi_fit->chiSquare();
   std::cout<< " chisquare"<< chisquare_costheta << "\n";
   
   TCanvas *cc = new TCanvas("cc", "cc",0,0,800,600);
   phi_fit->Draw();
   
   
   
   
   
   
   
   
   
   /*RooPlot* pt = pt_var.frame(Title("p_{T} GeV/c"),Bins(nbins));
   data_SigReg->plotOn(pt,DataError(RooAbsData::SumW2));
   RooPlot* phi = phi_var.frame(Title("#phi (rad)"),Bins(nbins));
   data_SigReg->plotOn(phi);
   RooPlot* pdtheta = pdtheta_var.frame(Title("p#bullet#delta#theta (rad#bulletGeV/c)"),Bins(nbins));
   data_SigReg->plotOn(pdtheta);
   TCanvas *lz = new TCanvas("c","c", 1200, 400);
   lz->Divide(3,1);
   lz->cd(1);pt->Draw();lz->cd(2);phi->Draw();lz->cd(3);pdtheta->Draw();
  */
}

