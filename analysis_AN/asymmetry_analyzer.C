#define asymmetry_analyzer_cxx
#include "asymmetry_analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void asymmetry_analyzer::Loop()
{
     /*
      In a ROOT session, you can do:
      Root > .L asymmetry_analyzer.C
      Root > asymmetry_analyzer t
      Root > t.GetEntry(12); // Fill t data members with entry number 12
      Root > t.Show();       // Show values of entry 12
      Root > t.Show(16);     // Read and show values of entry 16
      Root > t.Loop();       // Loop on all entries
*/

	TH1F* xf_h = new TH1F("xf_h", "xf_h", 100, -0.1. 0.1); 

 
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   Float_t s = 200 ; //centre of mass energy   








	Float_t  pT =0, pdtheta =0, ref_rad =0, trk_phi =0, xF=0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
	pT = sqrt(smpx*smpx + smpy*smpy);
	xF = (2*smpz)/(sqrt(s));
        xf_h->Fill(xF);     


   }
}
