#include <TTree.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TMath.h>

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

void test_check(TString filename ){

        TFile* file;
	TTree* muon;
	ifstream inrootfiles;
	inrootfiles.open(filename.Data());
	while(!inrootfiles.eof())
		{
			TString tssafiles;
			inrootfiles>>tssafiles;
			if(inrootfiles.eof())break;
			//std::cout<<" my files , technically the data files :"<< tssafiles.Data()<<"\n";
			file = new TFile(tssafiles.Data());
			muon = (TTree*)file->Get("T");	
			Int_t nentries = muon->GetEntries();
         		std::cout<<"No of entries"<<nentries<<"\n";
         		for (Int_t ientry=0; ientry<nentries; ientry++) {
	         		muon->GetEntry(ientry);
        	 		if (ientry%1000==0) cout << "processing event " << ientry << "/" << nentries <<"\n";
 	}

}

}
