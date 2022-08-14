/*

Author : Muhammad Alibordi 
A simple macro , runs in a root session on 
the pico DST files primarily produced by Sanghwa Park
The Kinematc Varibales are concerned only to the Open Muon Heavy Flavour 
decays process for the moment. Gathering the kinematics of TSSA analysis

*/
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


const int NB = 120;

int default_qa;
int runnumber;
int fillnumber;
int xingshift;
int badrun_flag;
int patternblue[NB];
int patternyell[NB];
int badbunch_qa[NB];

long long scaler_bbcvtxcut[NB];
long long scaler_bbcnovtx[NB];
long long scaler_zdcwide[NB];
long long scaler_zdcnarrow[NB];

float b_pol, y_pol;
float b_stat, y_stat;
float b_syst, y_syst;

int pattern_number;
int group;

std::vector<int> base1;
std::vector<int> base2;
std::vector<int> base3;
std::vector<int> base4;
std::vector<int> base1a;
std::vector<int> base2a;
std::vector<int> base3a;
std::vector<int> base4a;
std::vector<int> v_ptb;
std::vector<int> v_pty;

int FindPattern(vector<int> &ptb, vector<int> &pty);
int FindGroup(int pattern);
void DefineBasePatterns();

int verb = 0;










void Variable_histograms(string procnum, string runlist){

	DefineBasePatterns();
	gSystem->Load("libuspin.so");
	SpinDBContent spin_cont;
	SpinDBOutput spin_out("phnxrc");

    TH1F *smddg0_h = new TH1F("smddg0_h", "smddg0_h;DDG0_{#mu}; Number of events (a.u.)", 100, 0, 25.0);
    TH1F *smpx_h = new TH1F("smpx_h", " smpx_h; px_{#mu} (GeV); Number of events", 100, 0, 10);
    TH1F *smpy_h = new TH1F("smpy_h", "smpy_h; py_{#mu} (GeV); Number of events ", 100, 0, 10);
    TH1F *smpz_h = new TH1F("smpz_h", "smpz_h; pz_{#mu} (GeV); Number of events", 100, 0, 30);
    TH1F *smrapidity_h = new TH1F("smrapidity_h", "smrapidity_h; y_{#mu}; Number of events", 100, -5.0, 5.0);
    TH1F *smdg0_h = new TH1F("smdg0_h", "smdg0_h; DG0_{#mu}; Number of events", 100, 0, 80);
    TH1F *smds3_h = new TH1F("smds3_h", "smds3_h; DS3_{#mu}; Number of events ", 100, 0, 80);
    TH1F *smtrchi2_h = new TH1F("smtrchi2_h", "smtrchi2_h; #chi^{2}_{Tr}; Number of events", 100, 0, 30);
    TH1F *smidchi2_h = new TH1F("smidchi2_h", "smidchi2_h; #chi^{2}_{ID}; Number of events", 100, 0,30);
    TH1F *smxst1_h = new TH1F("smxst1_h", "smxst1_h; xst_{1}; Number of events", 100, -200, 200);
    TH1F *smxst2_h = new TH1F("smxst2_h", "smxst2_h; xst_{2}; Number of events", 100, -200, 200);
    TH1F *smxst3_h = new TH1F("smxst3_h", "smxst3_h; xst_{3}; Number of events", 100, -200, 200);
    TH1F *smyst1_h = new TH1F("smyst1_h", "smyst1_h; yst_{1}; Number of events", 100, -200, 200);
    TH1F *smyst2_h = new TH1F("smyst2_h", "smyst2_h; yst_{2}; Number of events", 100, -200, 200);
    TH1F *smyst3_h = new TH1F("smyst3_h", "smyst3_h; yst_{3}; Number of events", 100, -200, 200); 
    TH1F *smidx_h = new TH1F("smidx_h", "smidx_h; idx; Number of events", 100, -300, 300);
    TH1F *smidy_h = new TH1F("smidy_h", "smidy_h; idy; Number of events", 100, -300, 300);
    TH1F *smst1px_h = new TH1F("smsta1px_h", "smsta1px_h; st_{1}^{px}; Number of events", 100, 0, 10);
    TH1F *smst1py_h = new TH1F("smsta1py_h", "smsta1py_h; st_{1}^{py}; Number of events", 100, 0, 10);
    TH1F *smst1pz_h = new TH1F("smsta1pz_h", "smsta1pz_h; st_{1}^{pz}; Number of events", 100, 0, 30);
    TH1F *smdcar_h = new TH1F("smdcar_h", "smdcar_h; DCA_{r} (length); Number of events", 100, 0, 100);
    TH1F *smdcaz_h = new TH1F("smdcaz_h", "smdcaz_h; DCA_{z} (length); Number of events", 100, 0, 100);
    TH1F *smtrhits_h = new TH1F("smtrhits_h", "smtrhits_h; Hit-Multiplicity_{Tr}; Number of events", 100, 10E6, 10E6);
    TH1F *smidhits_h = new TH1F("smidhits_h", "smidhits_h; Hit-Multiplicity_{ID}; Number of events", 100, 0, 10E4);
    TH1F *smntrhits_h = new TH1F("smntrhits_h", "smntrhits_h; Hit-Multiplicity_{nTr}; Number of events", 100, 0, 40);
    TH1F *smnidhits_h = new TH1F("smnidhits_h", "smnidhits_h; Hit-Multiplicity_{nID}; Number of events", 100, 0, 40);
    TH1F *smmuid1d_h = new TH1F("smmuid1d_h", "smmuid1d_h; MUID_{1D}; Number of events", 100, -2,2);
    TH1F *smmuid1s_h = new TH1F("smmuid1s_h", "smmuid1s_h; MUID_{1S}; Number of events", 100, -2,2); 
    TH1F *smcharge_h = new TH1F("smcharge_h", "smcharge_h; Charge_{#mu}; Number of events", 100, -2,2);
    TH1F *smx0_h = new TH1F("smx0_h", "smx0_h; x0_{#mu}; Number of events", 100, 0, 100);
    TH1F *smy0_h = new TH1F("smy0_h", "smy0_h; y0_{#mu}; Number of events", 100, 0, 100);
    TH1F *smz0_h = new TH1F("smz0_h", "smz0_h; z0_{#mu}; Number of events", 100, 0, 100);
    TH1F *smcov_h = new TH1F("smcov_h", "smcov_h; cov_{#mu}; Number of events" , 100, -10, 100);
    TH1F *smx0fvtxmutr_h = new TH1F("smx0fvtxmutr_h", "smx0fvtxmutr_h; x0_{fvtx}^{mutr} (length); Number of events", 100, 0, 10);
    TH1F *smy0fvtxmutr_h = new TH1F("smy0fvtxmutr_h", "smy0fvtxmutr_h; y0_{fvtx}^{mutr} (length); Number of events", 100, 0, 10);
    TH1F *smz0fvtxmutr_h = new TH1F("smz0fvtxmutr_h", "smz0fvtxmutr_h; z0_{fvtx}^{mutr} (length); Number of events", 100, 0, 20);
    TH1F *smpxfvtxmutr_h = new TH1F("smpxfvtxmutr_h", "smpxfvtxmutr_h; px_{fvtx}^{mutr} (GeV); Number of events", 100, 0, 20);
    TH1F *smpyfvtxmutr_h = new TH1F("smpyfvtxmutr_h", "smpyfvtxmutr_h; py_{fvtx}^{mutr} (GeV); Number of events", 100, 0, 20);
    TH1F *smpzfvtxmutr_h = new TH1F("smpzfvtxmutr_h", "smpzfvtxmutr_h; pz_{fvtx}^{mutr} (GeV); Number of events", 100, 0, 30);
    TH1F *smdphifvtx_h = new TH1F("smdphifvtx_h", "smdphifvtx_h; d#phi_{fvtx} ; Number of events", 100, -TMath::Pi(), TMath::Pi() );
    TH1F *smdrfvtx_h = new TH1F("smdrfvtx_h", "smdrfvtx_h; dr_{fvtx}; Number of events", 100, -20, 20);
    TH1F *smdthetafvtx_h = new TH1F("smdthetafvtx_h", "smdthetafvtx_h; d#theta_{fvtx} ; Number of events", 100, -0.2, 0.5*TMath::Pi()); 
    TH1F *smchi2fvtx_h = new TH1F("smchi2fvtx_h", "smchi2fvtx_h; #chi^{2}_{fvtx} ; Number of events", 100, 0, 20);
    TH1F *smclusterssize_h = new TH1F("smclusterssize_h", "smclusterssize_h; clusterssize ; Number of events", 100, 0, 40);
    TH1F *smdcaphi_h = new TH1F("smdcaphi_h", "smdcaphi_h; DCA#phi ; Number of events", 100, -TMath::Pi(), TMath::Pi());
    TH1F *smfvtxtrackid_h = new TH1F("smfvtxtrackid_h", "smfvtxtrackid_h; trackID_{fvtx}; Number of events", 100, -10, 10);
    TH1F *smvtxindex_h = new TH1F("smvtxindex_h", "smvtxindex_h; Vetex-Index; Number of events", 100, 0, 200);
    TH1F *smhitpattern_h = new TH1F("smhitpattern_h", "smhitpattern_h; Hit-Pattern; Number of events", 100, 0, 400);
   

    /*TFile * fn1;
    TTree * singMu; 
    std::ifstream datafile ("input_files.txt");
    if ( datafile.is_open() ) {
		string inputchar;
		while (getline(datafile,inputchar)) {
		std::cout << " data file in each line   :: "<<inputchar <<"\n" ;
    		fn1 = new TFile(("../../data/"+inputchar).c_str());
    		singMu = (TTree*)fn1->Get("T");
}}*/
    //TFile * fn1 = new TFile(("../../data/"+procnum+".root").c_str());
    for (unsigned i = 421815 ; i <= 432008; ++i){
    TFile *fn1 = new TFile(Form("../../../data/%d.root", i));
    if ((!fn1) || (fn1->IsZombie())) continue;
    TTree * singMu = (TTree*)fn1->Get("T");
    TLeaf* l_nsm; //index parameter
    TLeaf* l_smddg0; TLeaf* l_smpx ; TLeaf* l_smpy; TLeaf* l_smpz; TLeaf* l_smrapidity ;
    TLeaf* l_smtrhits; TLeaf* l_smidhits; TLeaf* l_smdg0; TLeaf* l_smds3;
    TLeaf* l_smtrchi2; TLeaf* l_smidchi2 ; TLeaf* l_smntrhits; TLeaf* l_smnidhits; TLeaf* l_smxst1 ;
    TLeaf* l_smxst2; TLeaf* l_smxst3; TLeaf* l_smyst1; TLeaf* l_smyst2; TLeaf* l_smyst3;
    TLeaf* l_smidx; TLeaf* l_smidy ; TLeaf* l_smst1px; TLeaf* l_smst1py; TLeaf* l_smst1pz ;
    TLeaf* l_smmuid1s; TLeaf* l_smmuid1d; TLeaf* l_smdcaz; TLeaf* l_smdcar; TLeaf* l_smcharge;
    TLeaf* l_smx0; TLeaf* l_smy0 ; TLeaf* l_smz0; TLeaf* l_smcov; TLeaf* l_smx0fvtx;
    TLeaf* l_smy0fvtx; TLeaf* l_smz0fvtx; TLeaf* l_smx0fvtxmutr; TLeaf* l_smy0fvtxmutr; 
    TLeaf* l_smz0fvtxmutr; TLeaf* l_smpxfvtxmutr; TLeaf* l_smpyfvtxmutr ; TLeaf* l_smpzfvtxmutr;
    TLeaf* l_smpxfvtx; TLeaf* l_smpyfvtx; TLeaf* l_smpzfvtx; TLeaf * l_smhitpattern;
    TLeaf* l_smnhitsfvtx; TLeaf* l_smdphifvtx; TLeaf* l_smdthetafvtx; TLeaf* l_smdrfvtx; 
    TLeaf* l_smchi2fvtx; TLeaf* l_smfvtxstrip; TLeaf* l_smcovfvtx; TLeaf* l_smcovfvtxmutr; 
    TLeaf* l_smfvtxcharge ; TLeaf* l_smnfvtxcharge; TLeaf* l_smnfvtxtrks; TLeaf* l_smmtootrkindex;
    TLeaf* l_smchi2fvtxmutr; TLeaf* l_smdcaphi; TLeaf* l_smclusterssize; TLeaf* l_smnfvtxtrackletscone; 
    TLeaf* l_smnfvtxtracklets; TLeaf* l_smnfvtxclusterscone; TLeaf* l_smfvtxtrackid ; TLeaf* l_smvtxindex; 
    TLeaf* l_smnmatching; TLeaf* l_hitpattern ; TLeaf* l_smcovmutstal; TLeaf* l_smxfvtxproj; TLeaf* l_smyfvtxproj;
    TLeaf* l_smxmutproj; TLeaf* l_smymutproj; TLeaf* l_smpxfvtxproj; TLeaf* l_smpyfvtxproj; TLeaf* l_drfvtxmutrsigma;
    TLeaf* l_smdphifvtxmutrsigma;
    TLeaf* l_smlastgap; TLeaf* l_smmaxres_sigma; TLeaf* l_smtrackid; 
    TLeaf* l_bbcn; TLeaf* l_bbcs; TLeaf* l_bbcqn; TLeaf*  l_bbcqs; TLeaf* l_bbcz; TLeaf* l_bbczerr;
    TLeaf* l_bbct0; TLeaf* l_bbcts; TLeaf* l_bbctn;
    TLeaf* l_evtbbcz; TLeaf* l_evtbbczerr; TLeaf* l_evtvtxx; TLeaf* l_evtvtxxerr; TLeaf* l_evtvtxy; TLeaf* l_evtvtxyerr;
    TLeaf* l_evtvtxz; TLeaf* l_evtvtxzerr;
    TLeaf* l_lvl1_trigraw; TLeaf* l_lvl1_triglive; TLeaf* l_lvl1_trigscaled; 
    TLeaf* l_lvl1_clock_cross; TLeaf* l_lvl1_rbits; TLeaf* l_beamclk;     





    l_nsm                    = (TLeaf*)singMu->GetLeaf("nSingleMuons");
    l_smddg0                 = (TLeaf*)singMu->GetLeaf("SingleMuons.DDG0");
    l_smpx                   = (TLeaf*)singMu->GetLeaf("SingleMuons.px");
    l_smpy                   = (TLeaf*)singMu->GetLeaf("SingleMuons.py");
    l_smpz                   = (TLeaf*)singMu->GetLeaf("SingleMuons.pz");
    l_smrapidity             = (TLeaf*)singMu->GetLeaf("SingleMuons.rapidity");
    l_smtrhits               = (TLeaf*)singMu->GetLeaf("SingleMuons.trhits");
    l_smidhits               = (TLeaf*)singMu->GetLeaf("SingleMuons.idhits");
    l_smdg0                  = (TLeaf*)singMu->GetLeaf("SingleMuons.DG0");
    l_smds3                  = (TLeaf*)singMu->GetLeaf("SingleMuons.DS3");
    l_smtrchi2               = (TLeaf*)singMu->GetLeaf("SingleMuons.trchi2");
    l_smidchi2               = (TLeaf*)singMu->GetLeaf("SingleMuons.idchi2");
    l_smntrhits              = (TLeaf*)singMu->GetLeaf("SingleMuons.ntrhits");
    l_smnidhits              = (TLeaf*)singMu->GetLeaf("SingleMuons.nidhits");
    l_smxst1                 = (TLeaf*)singMu->GetLeaf("SingleMuons.xst1");
    l_smxst2                 = (TLeaf*)singMu->GetLeaf("SingleMuons.xst2");
    l_smxst3                 = (TLeaf*)singMu->GetLeaf("SingleMuons.xst3");
    l_smyst1                 = (TLeaf*)singMu->GetLeaf("SingleMuons.yst1");
    l_smyst2                 = (TLeaf*)singMu->GetLeaf("SingleMuons.yst2");
    l_smyst3                 = (TLeaf*)singMu->GetLeaf("SingleMuons.yst3");
    l_smidx                  = (TLeaf*)singMu->GetLeaf("SingleMuons.idx");
    l_smidy                  = (TLeaf*)singMu->GetLeaf("SingleMuons.idy");
    l_smst1px                = (TLeaf*)singMu->GetLeaf("SingleMuons.st1px");
    l_smst1py                = (TLeaf*)singMu->GetLeaf("SingleMuons.st1py");
    l_smst1pz                = (TLeaf*)singMu->GetLeaf("SingleMuons.st1pz");
    l_smmuid1s               = (TLeaf*)singMu->GetLeaf("SingleMuons.MUID1S");
    l_smmuid1d               = (TLeaf*)singMu->GetLeaf("SingleMuons.MUID1D");
    l_smdcar                 = (TLeaf*)singMu->GetLeaf("SingleMuons.dca_r");
    l_smdcaz                 = (TLeaf*)singMu->GetLeaf("SingleMuons.dca_z");
    l_smcharge               = (TLeaf*)singMu->GetLeaf("SingleMuons.charge");
    l_smx0                   = (TLeaf*)singMu->GetLeaf("SingleMuons.x0");
    l_smy0                   = (TLeaf*)singMu->GetLeaf("SingleMuons.y0");
    l_smz0                   = (TLeaf*)singMu->GetLeaf("SingleMuons.z0");
    l_smcov                  = (TLeaf*)singMu->GetLeaf("SingleMuons.cov");
//forward vertex information 
    l_smx0fvtx               = (TLeaf*)singMu->GetLeaf("SingleMuons.x0_fvtx");
    l_smy0fvtx               = (TLeaf*)singMu->GetLeaf("SingleMuons.y0_fvtx");
    l_smz0fvtx               = (TLeaf*)singMu->GetLeaf("SingleMuons.z0_fvtx");
    l_smpxfvtx               = (TLeaf*)singMu->GetLeaf("SingleMuons.px_fvtx");
    l_smpyfvtx               = (TLeaf*)singMu->GetLeaf("SingleMuons.py_fvtx");
    l_smpzfvtx               = (TLeaf*)singMu->GetLeaf("SingleMuons.pz_fvtx");
    l_smx0fvtxmutr           = (TLeaf*)singMu->GetLeaf("SingleMuons.x0_fvtxmutr");
    l_smy0fvtxmutr           = (TLeaf*)singMu->GetLeaf("SingleMuons.y0_fvtxmutr");
    l_smz0fvtxmutr           = (TLeaf*)singMu->GetLeaf("SingleMuons.z0_fvtxmutr");
    l_smpxfvtxmutr           = (TLeaf*)singMu->GetLeaf("SingleMuons.px_fvtxmutr");
    l_smpyfvtxmutr           = (TLeaf*)singMu->GetLeaf("SingleMuons.py_fvtxmutr");
    l_smpzfvtxmutr           = (TLeaf*)singMu->GetLeaf("SingleMuons.pz_fvtxmutr");
    l_smnhitsfvtx            = (TLeaf*)singMu->GetLeaf("SingleMuons.nhits_fvtx");
    l_smdphifvtx             = (TLeaf*)singMu->GetLeaf("SingleMuons.dphi_fvtx");
    l_smdthetafvtx           = (TLeaf*)singMu->GetLeaf("SingleMuons.dtheta_fvtx");
    l_smdrfvtx               = (TLeaf*)singMu->GetLeaf("SingleMuons.dr_fvtx");
    l_smchi2fvtx             = (TLeaf*)singMu->GetLeaf("SingleMuons.chi2_fvtx");
    l_smfvtxstrip            = (TLeaf*)singMu->GetLeaf("SingleMuons.fvtx_strip");
    l_smcovfvtx              = (TLeaf*)singMu->GetLeaf("SingleMuons.cov_fvtx");
    l_smcovfvtxmutr          = (TLeaf*)singMu->GetLeaf("SingleMuons.cov_fvtxmutr");
    l_smfvtxcharge           = (TLeaf*)singMu->GetLeaf("SingleMuons.fvtx_charge");
    l_smnfvtxtrks            = (TLeaf*)singMu->GetLeaf("SingleMuons.nfvtx_trks");
    l_smmtootrkindex         = (TLeaf*)singMu->GetLeaf("SingleMuons.mutoo_trk_index");
    l_smchi2fvtxmutr         = (TLeaf*)singMu->GetLeaf("SingleMuons.chi2_fvtxmutr");
    l_smdcaphi               = (TLeaf*)singMu->GetLeaf("SingleMuons.dca_phi");
    l_smclusterssize         = (TLeaf*)singMu->GetLeaf("SingleMuons.clusters_size1");
    l_smnfvtxtrackletscone   = (TLeaf*)singMu->GetLeaf("SingleMuons.nfvtx_tracklets_cone");
    l_smnfvtxtracklets       = (TLeaf*)singMu->GetLeaf("SingleMuons.nfvtx_tracklets");
    l_smnfvtxclusterscone    = (TLeaf*)singMu->GetLeaf("SingleMuons.nfvtx_clusters_cone");
    l_smfvtxtrackid          = (TLeaf*)singMu->GetLeaf("SingleMuons.fvtxtrack_id");
    l_smvtxindex             = (TLeaf*)singMu->GetLeaf("SingleMuons.vtx_index");
    l_smhitpattern           = (TLeaf*)singMu->GetLeaf("SingleMuons.hit_pattern");
    l_smcovmutstal           = (TLeaf*)singMu->GetLeaf("SingleMuons.cov_mutstal");
    l_smxfvtxproj            = (TLeaf*)singMu->GetLeaf("SingleMuons.x_fvtxproj");
    l_smyfvtxproj            = (TLeaf*)singMu->GetLeaf("SingleMuons.y_fvtxproj");
    l_smxmutproj             = (TLeaf*)singMu->GetLeaf("SingleMuons.x_mutproj");
    l_smymutproj             = (TLeaf*)singMu->GetLeaf("SingleMuons.y_mutproj");
    l_smpxfvtxproj           = (TLeaf*)singMu->GetLeaf("SingleMuons.px_fvtxproj");
    l_smpyfvtxproj           = (TLeaf*)singMu->GetLeaf("SingleMuons.py_fvtxproj");
//bbc vertex     
    l_bbcn                   = (TLeaf*)singMu->GetLeaf("bbcn");
    l_bbcs                   = (TLeaf*)singMu->GetLeaf("bbcs");
    l_bbcqn                  = (TLeaf*)singMu->GetLeaf("bbcqn");
    l_bbcqs                  = (TLeaf*)singMu->GetLeaf("bbcqs");
    l_bbcz                   = (TLeaf*)singMu->GetLeaf("bbcz");
    l_bbczerr                = (TLeaf*)singMu->GetLeaf("bbczerr");
    l_bbct0                  = (TLeaf*)singMu->GetLeaf("bbct0");
    l_bbcts                  = (TLeaf*)singMu->GetLeaf("bbcts");
    l_bbctn                  = (TLeaf*)singMu->GetLeaf("bbctn");
    l_evtbbcz		     = (TLeaf*)singMu->GetLeaf("Evt_bbcZ");
    l_evtbbczerr             = (TLeaf*)singMu->GetLeaf("Evt_bbcZ_Err");
    l_evtvtxx                = (TLeaf*)singMu->GetLeaf("Evt_vtxX");
    l_evtvtxxerr             = (TLeaf*)singMu->GetLeaf("Evt_vtxX_Err"); 
    l_evtvtxy                = (TLeaf*)singMu->GetLeaf("Evt_vtxY");
    l_evtvtxyerr             = (TLeaf*)singMu->GetLeaf("Evt_vtxY_Err"); 
    l_evtvtxz                = (TLeaf*)singMu->GetLeaf("Evt_vtxZ");
    l_evtvtxzerr             = (TLeaf*)singMu->GetLeaf("Evt_vtxZ_Err"); 
    l_smmaxres_sigma         = (TLeaf*)singMu->GetLeaf("SingleMuons.maxres_sigma");
    l_smtrackid              = (TLeaf*)singMu->GetLeaf("SingleMuons.track_id");
    l_smlastgap              = (TLeaf*)singMu->GetLeaf("SingleMuons.lastgap");
//trigger and clock 
    l_lvl1_trigraw           = (TLeaf*)singMu->GetLeaf("lvl1_trigraw");
    l_lvl1_triglive          = (TLeaf*)singMu->GetLeaf("lvl1_triglive");
    l_lvl1_trigscaled        = (TLeaf*)singMu->GetLeaf("lvl1_trigscaled");
    l_lvl1_clock_cross       = (TLeaf*)singMu->GetLeaf("lvl1_clock_cross");
    l_lvl1_rbits             = (TLeaf*)singMu->GetLeaf("lvl1_rbits");
    l_beamclk                = (TLeaf*)singMu->GetLeaf("beamclk");







    //------Variable type initialization 
    int nmuons ;
    float smddg0, smpx, smpy, smpz, smrapidity, smdg0, smds3, smtrchi2, smidchi2, smxst1, smxst2, smxst3 , smyst1, smyst2, smyst3, smidx, smidy, smst1px, smst1py, smst1pz;
    float smdcar, smdcaz, smx0, smy0, smz0, smcov;
    float smx0fvtxmutr, smy0fvtxmutr, smz0fvtxmutr, smpxfvtxmutr, smpyfvtxmutr, smpzfvtxmutr, smdphifvtx, smdrfvtx, smdthetafvtx, smchi2fvtx,  smdcaphi;
    float smx0fvtx, smy0fvtx, smz0fvtx, smpxfvtx, smpyfvtx, smpzfvtx, smxfvtxproj , smyfvtxproj, smxmutproj,smymutproj, smpxfvtxproj, smpyfvtxproj;
    float bbcqn, bbcqs, bbcz, bbczerr, bbct0, bbcts, bbctn;
    int smfvtxtrackid, smvtxindex, smhitpattern;   
    int smtrhits , smidhits, smntrhits, smnidhits, smclusterssize, bbcn, bbcs; 
    bool smmuid1d; bool smmuid1s; bool smcharge;
    float smmaxres_sigma, evtbbcz, evtbbczerr, evtvtxx, evtvtxxerr, evtvtxy, evtvtxyerr, evtvtxz, evtvtxzerr; 
    int smlastgap, smtrackid;
    int lvl1_trigraw, lvl1_triglive, lvl1_trigscaled,lvl1_clock_cross, lvl1_rbits, beamclk;


    TFile *f = new TFile(Form("/direct/phenix+u/alibordi/hf_outputs/analysis_preliminary_%d.root",i),"recreate");
    TTree *analysis = new TTree("analysis","selected ntcltestle");
    analysis->Branch("nmuons",&nmuons,"nmuons/I");
    analysis->Branch("smddg0",&smddg0,"smdd0/F");
    analysis->Branch("smpx",&smpx,"smpx/F");
    analysis->Branch("smpy",&smpy,"smpy/F");
    analysis->Branch("smpz",&smpz,"smpz/F");
    analysis->Branch("smrapidity",&smrapidity,"smrapidity/F");
    analysis->Branch("smdg0",&smdg0,"smdg0/F");
    analysis->Branch("smds3",&smds3,"smds3/F");
    analysis->Branch("smtrchi2",&smtrchi2,"smtrchi2/F");
    analysis->Branch("smidchi2",&smidchi2,"smidchi2/F");
    analysis->Branch("smxst1",&smxst1,"smxst1/F");
    analysis->Branch("smxst2",&smxst2,"smxst2/F");
    analysis->Branch("smxst3",&smxst3,"smxst3/F");
    analysis->Branch("smyst1",&smyst1,"smyst1/F");
    analysis->Branch("smyst2",&smyst2,"smyst2/F");
    analysis->Branch("smyst3",&smyst3,"smyst3/F");
    analysis->Branch("smidx",&smidx,"smidx/F");
    analysis->Branch("smidy",&smidy,"smidy/F");
    analysis->Branch("smst1px",&smst1px,"smst1px/F");
    analysis->Branch("smst1py",&smst1py,"smst1py/F");
    analysis->Branch("smst1pz",&smst1pz,"smst1pz/F");
    analysis->Branch("smdcar",&smdcar,"smdcar/F");
    analysis->Branch("smdcaz",&smdcaz,"smdcaz/F");
    analysis->Branch("smx0",&smx0,"smx0/F");
    analysis->Branch("smx0",&smx0,"smx0/F");
    analysis->Branch("smy0",&smy0,"smy0/F");
    analysis->Branch("smz0",&smz0,"smz0/F");
    analysis->Branch("smtrhits",&smtrhits,"smtrhits/I");
    analysis->Branch("smidhits",&smidhits,"smidhits/I");
    analysis->Branch("smntrhits",&smntrhits,"smntrhits/I");
    analysis->Branch("smnidhits",&smnidhits,"smnidhits/I");
    analysis->Branch("smmuid1d",&smmuid1d,"smmuid1d/O");
    analysis->Branch("smmuid1s",&smmuid1s,"smmuid1s/O");
    analysis->Branch("smcharge",&smcharge,"smcharge/O"); 
    analysis->Branch("bbcn",&bbcn,"bbcn/I");
    analysis->Branch("bbcs",&bbcs,"bbcs/I");
    analysis->Branch("bbcqn",&bbcqn,"bbcqn/F");
    analysis->Branch("bbcqs",&bbcqs,"bbcqs/F");
    analysis->Branch("bbcz",&bbcz,"bbcz/F");
    analysis->Branch("bbczerr",&bbczerr,"bbczerr/F");
    analysis->Branch("bbct0",&bbct0,"bbct0/F");
    analysis->Branch("bbcts",&bbcts,"bbcts/F");
    analysis->Branch("bbctn",&bbctn,"bbctn/F");
    analysis->Branch("evtbbcz", &evtbbcz, "evtbbcz/F");
    analysis->Branch("evtbbczerr", &evtbbczerr, "evtbbczerr/F");
    analysis->Branch("evtvtxx", &evtvtxx, "evtvtxx/F");
    analysis->Branch("evtvtxxerr", &evtvtxxerr, "evtvtxxerr/F");
    analysis->Branch("evtvtxy", &evtvtxy, "evtvtxy/F");
    analysis->Branch("evtvtxyerr", &evtvtxyerr, "evtvtxyerr/F");
    analysis->Branch("evtvtxz", &evtvtxz, "evtvtxz/F");
    analysis->Branch("evtvtxzerr", &evtvtxzerr, "evtvtxzerr/F");
    analysis->Branch("smmaxres_sigma", &smmaxres_sigma, "smmaxres_sigma/F");
    analysis->Branch("smtrackid", &smtrackid, "smtrackid/I");
    analysis->Branch("smlastgap", &smlastgap, "smlastgap/I");
    analysis->Branch("smhitpattern", &smhitpattern, "smhitpattern/I");
    analysis->Branch("lvl1_trigraw", &lvl1_trigraw, "lvl1_trigraw/I");
    analysis->Branch("lvl1_triglive", &lvl1_triglive, "lvl1_triglive/I");
    analysis->Branch("lvl1_trigscaled", &lvl1_trigscaled, "lvl1_trigscaled/I");
    analysis->Branch("lvl1_clock_cross", &lvl1_clock_cross, "lvl1_clock_cross/I");
    analysis->Branch("lvl1_rbits", &lvl1_rbits, "lvl1_rbits/I");
    analysis->Branch("beamclk", &beamclk, "beamclk/I");
   
analysis->Branch("runnumber",      &runnumber,         "runnumber/I");
analysis->Branch("fillnumber",     &fillnumber,        "fillnumber/I");
analysis->Branch("qa_level",       &default_qa,        "qa_level/I");
analysis->Branch("xingshift",      &xingshift,         "xingshift/I");
analysis->Branch("badrun_flag",    &badrun_flag,       "badrun_flag/I");
analysis->Branch("polblue",        &b_pol,             "polblue/F");
analysis->Branch("polblue_stat",   &b_stat,            "polblue_stat/F");
analysis->Branch("polblue_sys",    &b_syst,            "polblue_sys/F");
analysis->Branch("polyellow",      &y_pol,             "polyellow/F");
analysis->Branch("polyellow_stat", &y_stat,            "polyellow_stat/F");
analysis->Branch("polyellow_sys",  &y_syst,            "polyellow_sys/F");
analysis->Branch("patternblue",    patternblue,        "patternblue[120]/I");
analysis->Branch("patternyellow",  patternyell,        "patternyellow[120]/I");
analysis->Branch("badbunch_qa",    badbunch_qa,        "badbunch_qa[120]/I");
analysis->Branch("scalerA",        scaler_bbcvtxcut,   "scalerA[120]/L");
analysis->Branch("scalerB",        scaler_bbcnovtx,    "scalerB[120]/L");
analysis->Branch("scalerC",        scaler_zdcwide,     "scalerC[120]/L");
analysis->Branch("scalerD",        scaler_zdcnarrow,   "scalerD[120]/L");
analysis->Branch("pattern",        &pattern_number,    "pattern/I");
analysis->Branch("group",          &group,             "group/I");


//forward vertex information 
/*analysis->Branch("smx0fvtxmutr",&smx0fvtxmutr,"smx0fvtxmutr/F");
analysis->Branch("smy0fvtxmutr",&smy0fvtxmutr,"smy0fvtxmutr/F");
analysis->Branch("smz0fvtxmutr",&smz0fvtxmutr,"smz0fvtxmutr/F");
analysis->Branch("smpxfvtxmutr",&smpxfvtxmutr,"smpxfvtxmutr/F");
analysis->Branch("smpyfvtxmutr",&smpyfvtxmutr,"smpyfvtxmutr/F");
analysis->Branch("smpzfvtxmutr",&smpzfvtxmutr,"smpzfvtxmutr/F");
analysis->Branch("smdcaphi",&smdcaphi,"smdcaphi/F");
analysis->Branch("smdrfvtx",&smdrfvtx,"smdrfvtx/F");
analysis->Branch("smchi2fvtx",&smchi2fvtx,"smchi2fvtx/F");
analysis->Branch("smdthetafvtx",&smdthetafvtx,"smdthetafvtx/F");
analysis->Branch("smdphifvtx",&smdphifvtx,"smdphifvtx/F");
analysis->Branch("smclusterssize",&smclusterssize,"smclusterssize/I");
analysis->Branch("smfvtxtrackid",&smclusterssize,"smfvtxtrackid/I");
analysis->Branch("smx0fvtx",&smx0fvtx,"smx0fvtx/F");
analysis->Branch("smy0fvtx",&smy0fvtx,"smy0fvtx/F");
analysis->Branch("smz0fvtx",&smz0fvtx,"smz0fvtx/F");
analysis->Branch("smpxfvtx",&smpxfvtx,"smpxfvtx/F");
analysis->Branch("smpyfvtx",&smpyfvtx,"smpyfvtx/F");
analysis->Branch("smpzfvtx",&smpzfvtx,"smpzfvtx/F");
*/
//analysis->Branch("smxfvtxproj",&smxfvtxproj,"smxfvtxproj/F");
//analysis->Branch("smyfvtxproj",&smyfvtxproj,"smyfvtxproj/F");
//analysis->Branch("smxmutproj",&smxmutproj,"smxmutproj/F");
//analysis->Branch("smymutproj",&smymutproj,"smymutproj/F");
//analysis->Branch("smpxfvtxproj",&smpxfvtxproj,"smpxfvtxproj/F");
//analysis->Branch("smpyfvtxproj",&smpyfvtxproj,"smpyfvtxproj/F");



  Int_t nentries = singMu->GetEntries(); 
    std::cout<<"No of entries"<<nentries<<"\n";
    for (Int_t ientry=0; ientry<nentries; ientry++) {               
      singMu->GetEntry(ientry);
      if (ientry%100==0) cout << "processing event " << ientry << "/" << nentries <<"\n";
      int nmuons = int(l_nsm->GetValue());
      //cout << nmuons << endl;
      
      // Looop over all single muons for each event
      for(int imuon=0; imuon<nmuons; imuon++) {
	smddg0 = float(l_smddg0->GetValue(imuon));
	smpx = float(l_smpx->GetValue(imuon));
        smpy = float(l_smpy->GetValue(imuon));
        smpz = float(l_smpz->GetValue(imuon));
        smrapidity = float(l_smrapidity->GetValue(imuon));
        smdg0 = float(l_smdg0->GetValue(imuon));
        smds3 = float(l_smds3->GetValue(imuon));
        smtrchi2 = float(l_smtrchi2->GetValue(imuon));
        smidchi2 = float(l_smidchi2->GetValue(imuon));
        smxst1 = float(l_smxst1->GetValue(imuon));
        smxst2 = float(l_smxst2->GetValue(imuon));
        smxst3 = float(l_smxst3->GetValue(imuon));
        smyst1 = float(l_smyst1->GetValue(imuon));
        smyst2 = float(l_smyst2->GetValue(imuon));
        smyst3 = float(l_smyst3->GetValue(imuon));
        smidx = float(l_smidx->GetValue(imuon));
        smidy = float(l_smidy->GetValue(imuon));
        smst1px = float(l_smst1px->GetValue(imuon));
        smst1py = float(l_smst1py->GetValue(imuon));
        smst1pz = float(l_smst1pz->GetValue(imuon));
        smdcar = float(l_smdcar->GetValue(imuon));
        smdcaz = float(l_smdcaz->GetValue(imuon));       
        smtrhits = int(l_smtrhits->GetValue(imuon));
        smidhits = int(l_smidhits->GetValue(imuon));
        smntrhits = int(l_smntrhits->GetValue(imuon));
        smnidhits = int(l_smnidhits->GetValue(imuon));
        smmuid1d = bool(l_smmuid1d->GetValue(imuon));
        smmuid1s = bool(l_smmuid1s->GetValue(imuon));
        smcharge = bool(l_smcharge->GetValue(imuon));
        smx0 = float(l_smx0->GetValue(imuon));
        smy0 = float(l_smy0->GetValue(imuon));
        smz0 = float(l_smz0->GetValue(imuon));
        smcov = float(l_smcov->GetValue(imuon));
	smvtxindex = int(l_smvtxindex->GetValue(imuon));
        smhitpattern = int(l_smhitpattern->GetValue(imuon));
        smmaxres_sigma = float(l_smmaxres_sigma->GetValue(imuon));
        smtrackid = int(l_smtrackid->GetValue(imuon));
        smlastgap = int(l_smlastgap->GetValue(imuon));


        smx0fvtxmutr = float(l_smx0fvtxmutr->GetValue(imuon));
	smy0fvtxmutr = float(l_smy0fvtxmutr->GetValue(imuon));
	smz0fvtxmutr = float(l_smz0fvtxmutr->GetValue(imuon));
	smpxfvtxmutr = float(l_smpxfvtxmutr->GetValue(imuon));
	smpyfvtxmutr = float(l_smpyfvtxmutr->GetValue(imuon));
	smpzfvtxmutr = float(l_smpzfvtxmutr->GetValue(imuon));
	smdcaphi = float(l_smdcaphi->GetValue(imuon));
	smdrfvtx = float(l_smdrfvtx->GetValue(imuon));
	smchi2fvtx = float(l_smchi2fvtx->GetValue(imuon));
	smdthetafvtx = float(l_smdthetafvtx->GetValue(imuon));
	smdphifvtx = float(l_smdphifvtx->GetValue(imuon));
	smclusterssize = float(l_smclusterssize->GetValue(imuon));
	smfvtxtrackid = int(l_smfvtxtrackid->GetValue(imuon));
	smx0fvtx = float(l_smx0fvtx->GetValue(imuon));
	smy0fvtx = float(l_smy0fvtx->GetValue(imuon));
	smz0fvtx = float(l_smz0fvtx->GetValue(imuon));
	smpxfvtx = float(l_smpxfvtx->GetValue(imuon));
	smpyfvtx = float(l_smpyfvtx->GetValue(imuon));
	smpzfvtx = float(l_smpzfvtx->GetValue(imuon));
/*	smxfvtxproj = float(l_smxfvtxproj->GetValue(imuon));
	//std::cout<<smxfvtxproj<<"\n";
	smyfvtxproj = float(l_smyfvtxproj->GetValue(imuon));
	smxmutproj = float(l_smxmutproj->GetValue(imuon));
	smymutproj = float(l_smymutproj->GetValue(imuon));
	smpxfvtxproj = float(l_smpxfvtxproj->GetValue(imuon));
	smpyfvtxproj = float(l_smpyfvtxproj->GetValue(imuon));       
  */      
        smddg0_h->Fill(smddg0);
        smpx_h->Fill(smpx);
        smpy_h->Fill(smpy);
        smpz_h->Fill(smpz);
        smrapidity_h->Fill(smrapidity);
        smdg0_h->Fill(smdg0);
        smds3_h->Fill(smds3);
        smtrchi2_h->Fill(smtrchi2);
        smidchi2_h->Fill(smidchi2);
        smxst1_h->Fill(smxst1);
        smxst2_h->Fill(smxst2);
        smxst3_h->Fill(smxst3);
        smidx_h->Fill(smidx);
        smidy_h->Fill(smidy);
        smst1px_h->Fill(smst1px);
        smst1py_h->Fill(smst1py);
        smst1pz_h->Fill(smst1pz);
        smdcar_h->Fill(smdcar);
        smdcaz_h->Fill(smdcaz);
        smtrhits_h->Fill(smtrhits);
        smidhits_h->Fill(smidhits);
        smntrhits_h->Fill(smntrhits);
        smnidhits_h->Fill(smnidhits);
        smmuid1d_h->Fill(smmuid1d);
        smmuid1s_h->Fill(smmuid1s);
        smcharge_h->Fill(smcharge);
        smx0_h->Fill(smx0);
        smy0_h->Fill(smy0);
        smz0_h->Fill(smz0);
        smcov_h->Fill(smcov);
        smx0fvtxmutr_h->Fill(smx0fvtxmutr);
        smy0fvtxmutr_h->Fill(smy0fvtxmutr);
        smz0fvtxmutr_h->Fill(smz0fvtxmutr);
        smpxfvtxmutr_h->Fill(smpxfvtxmutr);
        smpyfvtxmutr_h->Fill(smpyfvtxmutr);
        smpzfvtxmutr_h->Fill(smpzfvtxmutr);
        smdphifvtx_h->Fill(smdphifvtx);
        smdrfvtx_h->Fill(smdrfvtx);
        smdthetafvtx_h->Fill(smdthetafvtx);
        smchi2fvtx_h->Fill(smchi2fvtx);
        smclusterssize_h->Fill(smclusterssize);
        smdcaphi_h->Fill(smdcaphi);
        smfvtxtrackid_h->Fill(smfvtxtrackid);
        smvtxindex_h->Fill(smvtxindex);
        smhitpattern_h->Fill(smhitpattern);
      } 
      

    bbcn  = int (l_bbcn->GetValue()); 
    bbcs  = int (l_bbcs->GetValue());
    bbcqn  = int (l_bbcqn->GetValue());
    bbcqs  = float (l_bbcqs->GetValue());
    bbcz   = float(l_bbcz->GetValue());
    bbczerr = float(l_bbczerr->GetValue());
    bbct0   = float(l_bbct0->GetValue());
    bbcts   = float(l_bbcts->GetValue());
    bbctn   = float(l_bbctn->GetValue());
    evtbbcz = float(l_evtbbcz->GetValue());
    evtbbczerr = float(l_evtbbczerr->GetValue());
    evtvtxx = float(l_evtvtxx->GetValue());
    evtvtxxerr = float(l_evtvtxxerr->GetValue());
    evtvtxy = float(l_evtvtxy->GetValue());
    evtvtxyerr = float(l_evtvtxyerr->GetValue());
    evtvtxz = float(l_evtvtxz->GetValue());
    evtvtxzerr = float(l_evtvtxzerr->GetValue());
    lvl1_trigraw = int(l_lvl1_trigraw->GetValue());
    lvl1_triglive = int(l_lvl1_triglive->GetValue());
    lvl1_trigscaled = int(l_lvl1_trigscaled->GetValue());
    lvl1_clock_cross = int(l_lvl1_clock_cross->GetValue());
    lvl1_rbits = int(l_lvl1_rbits->GetValue());
    beamclk = int(l_beamclk->GetValue());



     analysis->Fill();


    }// Event loop
    


/*bbcvertexcut->ScalerA, bbcwithoutcut->ScalerB, zdcnarrow->ScalerC, zdcwide-> ScalerD
https://www.phenix.bnl.gov/WWW/offline/wikioff/index.php/GL1P_Information
*/



ifstream ifstr(runlist.c_str());
int run;
while( ifstr >> run )
{
   InitTreeVars();
   
   v_ptb.clear();
   v_pty.clear();
   
   int qa_level  = spin_out.GetDefaultQA(run);
   spin_out.StoreDBContent(run, run);
   int dbout = spin_out.GetDBContentStore(spin_cont, run);
   if(dbout != 1) continue;
   
   default_qa = qa_level;
   runnumber  = spin_cont.GetRunNumber();
   fillnumber = spin_cont.GetFillNumber();
   xingshift  = spin_cont.GetCrossingShift();
   badrun_flag = spin_cont.GetBadRunFlag();
   
   spin_cont.GetPolarizationBlue(1, b_pol, b_stat, b_syst);
   spin_cont.GetPolarizationYellow(1, y_pol, y_stat, y_syst);
   

   for(int ibunch=0; ibunch<NB; ibunch++)
      {
      patternblue[ibunch] = spin_cont.GetSpinPatternBlue(ibunch);
      patternyell[ibunch] = spin_cont.GetSpinPatternYellow(ibunch);
      badbunch_qa[ibunch] = spin_cont.GetBadBunchFlag(ibunch);
      
      scaler_bbcvtxcut[ibunch] = spin_cont.GetScalerBbcVertexCut(ibunch);
      scaler_bbcnovtx[ibunch]  = spin_cont.GetScalerBbcNoCut(ibunch);
      scaler_zdcwide[ibunch]   = spin_cont.GetScalerZdcWide(ibunch);
      scaler_zdcnarrow[ibunch] = spin_cont.GetScalerZdcNarrow(ibunch);
      
         
      if(ibunch < 16)
         {
         v_ptb.push_back(patternblue[ibunch]);
         v_pty.push_back(patternyell[ibunch]);
         }
      }
   
   pattern_number = FindPattern(v_ptb, v_pty);
   group = FindGroup(pattern_number);
   
   analysis->Fill();
   
   if(verb)
      {
      cout << endl;
      cout << "Run: " << runnumber << "\n"
      << "QA level: " << qa_level << "\n"
      << "Fillnumber: " << fillnumber << "\n"
      << "Polarization: " << b_pol << " (blue) " << y_pol << " (yellow) \n"
      << "crossing shift: " << xingshift << endl;
      
      cout << "Spin pattern blue: ";
      for(int i=0; i<NB; i++) cout << patternblue[i] << " ";
      cout << endl;
      cout << "Spin pattern yellow: ";
      for(int i=0; i<NB; i++) cout << patternblue[i] << " ";
      cout << endl;
      }
}



     analysis->Print();
     f->Write();
}
/*
     TCanvas *pl1 = new TCanvas("pl1", "pl1", 1000, 800);
     pl1->Divide(5,3);
     pl1->cd(1);smddg0_h->SetFillColor(kGreen+2);smddg0_h->SetFillStyle(3644);smddg0_h->SetLineColor(kGray+2);smddg0_h->Draw();
     pl1->cd(2);smpx_h->SetFillColor(kGreen+2);smpx_h->SetFillStyle(3644);smpx_h->SetLineColor(kGray+2);smpx_h->Draw();
     pl1->cd(3);smpy_h->SetFillColor(kGreen+2);smpy_h->SetFillStyle(3644);smpy_h->SetLineColor(kGray+2);smpy_h->Draw();
     pl1->cd(4);smpz_h->SetFillColor(kGreen+2);smpz_h->SetFillStyle(3644);smpz_h->SetLineColor(kGray+2);smpz_h->Draw();
     pl1->cd(5);smrapidity_h->SetFillColor(kGreen+2);smrapidity_h->SetFillStyle(3644);smrapidity_h->SetLineColor(kGray+2);smrapidity_h->Draw();
     pl1->cd(6);smdg0_h->SetFillColor(kGreen+2);smdg0_h->SetFillStyle(3644);smdg0_h->SetLineColor(kGray+2);smdg0_h->Draw();
     pl1->cd(7);smds3_h->SetFillColor(kGreen+2);smds3_h->SetFillStyle(3644);smds3_h->SetLineColor(kGray+2);smds3_h->Draw();
     pl1->cd(8);smtrchi2_h->SetFillColor(kGreen+2);smtrchi2_h->SetFillStyle(3644);smtrchi2_h->SetLineColor(kGray+2);smtrchi2_h->Draw();
     pl1->cd(9);smidchi2_h->SetFillColor(kGreen+2);smidchi2_h->SetFillStyle(3644);smidchi2_h->SetLineColor(kGray+2);smidchi2_h->Draw();
     pl1->cd(10);smxst1_h->SetFillColor(kGreen+2);smxst1_h->SetFillStyle(3644);smxst1_h->SetLineColor(kGray+2);smxst1_h->Draw();
     pl1->cd(11);smxst2_h->SetFillColor(kGreen+2);smxst2_h->SetFillStyle(3644);smxst2_h->SetLineColor(kGray+2);smxst2_h->Draw();
     pl1->cd(12);smxst3_h->SetFillColor(kGreen+2);smxst3_h->SetFillStyle(3644);smxst3_h->SetLineColor(kGray+2);smxst3_h->Draw();
     pl1->cd(13);smidx_h->SetFillColor(kGreen+2);smidx_h->SetFillStyle(3644);smidx_h->SetLineColor(kGray+2);smidx_h->Draw();
     pl1->cd(14);smidy_h->SetFillColor(kGreen+2);smidy_h->SetFillStyle(3644);smidy_h->SetLineColor(kGray+2);smidy_h->Draw();
     pl1->cd(15);smst1px_h->SetFillColor(kGreen+2);smst1px_h->SetFillStyle(3644);smst1px_h->SetLineColor(kGray+2);smst1px_h->Draw();
     pl1->SaveAs("/direct/phenix+u/alibordi/hf_outputs/Variable_WOSelection_set1.pdf","pdf");
     TCanvas *pl2 = new TCanvas("pl2", "pl2", 1000, 800);
     pl2->Divide(5,3);
     pl2->cd(1);smst1py_h->SetFillColor(kGreen+2);smst1py_h->SetFillStyle(3644);smst1py_h->SetLineColor(kGray+2);smst1py_h->Draw();
     pl2->cd(2);smst1pz_h->SetFillColor(kGreen+2);smst1pz_h->SetFillStyle(3644);smst1pz_h->SetLineColor(kGray+2);smst1pz_h->Draw();
     pl2->cd(3);smdcar_h->SetFillColor(kGreen+2);smdcar_h->SetFillStyle(3644);smdcar_h->SetLineColor(kGray+2);smdcar_h->Draw();
     pl2->cd(4);smdcaz_h->SetFillColor(kGreen+2);smdcaz_h->SetFillStyle(3644);smdcaz_h->SetLineColor(kGray+2);smdcaz_h->Draw();
     pl2->cd(5);smtrhits_h->SetFillColor(kGreen+2);smtrhits_h->SetFillStyle(3644);smtrhits_h->SetLineColor(kGray+2);smtrhits_h->Draw();
     pl2->cd(6);smidhits_h->SetFillColor(kGreen+2);smidhits_h->SetFillStyle(3644);smidhits_h->SetLineColor(kGray+2);smidhits_h->Draw();
     pl2->cd(7);smntrhits_h->SetFillColor(kGreen+2);smntrhits_h->SetFillStyle(3644);smntrhits_h->SetLineColor(kGray+2);smntrhits_h->Draw();
     pl2->cd(8);smnidhits_h->SetFillColor(kGreen+2);smnidhits_h->SetFillStyle(3644);smnidhits_h->SetLineColor(kGray+2);smnidhits_h->Draw();
     pl2->cd(9);smmuid1s_h->SetFillColor(kGreen+2);smmuid1s_h->SetFillStyle(3644);smmuid1s_h->SetLineColor(kGray+2);smmuid1s_h->Draw();
     pl2->cd(10);smmuid1d_h->SetFillColor(kGreen+2);smmuid1d_h->SetFillStyle(3644);smmuid1d_h->SetLineColor(kGray+2);smmuid1d_h->Draw();
     pl2->cd(11);smcharge_h->SetFillColor(kGreen+2);smcharge_h->SetFillStyle(3644);smcharge_h->SetLineColor(kGray+2);smcharge_h->Draw();
     pl2->cd(12);smx0_h->SetFillColor(kGreen+2);smx0_h->SetFillStyle(3644);smx0_h->SetLineColor(kGray+2);smx0_h->Draw();
     pl2->cd(13);smy0_h->SetFillColor(kGreen+2);smy0_h->SetFillStyle(3644);smy0_h->SetLineColor(kGray+2);smy0_h->Draw();
     pl2->cd(14);smz0_h->SetFillColor(kGreen+2);smz0_h->SetFillStyle(3644);smz0_h->SetLineColor(kGray+2);smz0_h->Draw();
     pl2->cd(15);smcov_h->SetFillColor(kGreen+2);smcov_h->SetFillStyle(3644);smcov_h->SetLineColor(kGray+2);smcov_h->Draw();
     pl2->SaveAs("/direct/phenix+u/alibordi/hf_outputs/Variable_WOSelection_set2.pdf","pdf");


     TCanvas *pl3 = new TCanvas("pl3", "pl3", 1000, 800);
     pl3->Divide(5,3);
     pl3->cd(1);smx0fvtxmutr_h->SetFillColor(kGreen+2);smx0fvtxmutr_h->SetFillStyle(3644);smx0fvtxmutr_h->SetLineColor(kGray+2);smx0fvtxmutr_h->Draw();
     pl3->cd(2);smy0fvtxmutr_h->SetFillColor(kGreen+2);smy0fvtxmutr_h->SetFillStyle(3644);smy0fvtxmutr_h->SetLineColor(kGray+2);smy0fvtxmutr_h->Draw();
     pl3->cd(3);smz0fvtxmutr_h->SetFillColor(kGreen+2);smz0fvtxmutr_h->SetFillStyle(3644);smz0fvtxmutr_h->SetLineColor(kGray+2);smz0fvtxmutr_h->Draw();
     pl3->cd(4);smpxfvtxmutr_h->SetFillColor(kGreen+2);smpxfvtxmutr_h->SetFillStyle(3644);smpxfvtxmutr_h->SetLineColor(kGray+2);smpxfvtxmutr_h->Draw();
     pl3->cd(5);smpyfvtxmutr_h->SetFillColor(kGreen+2);smpyfvtxmutr_h->SetFillStyle(3644);smpyfvtxmutr_h->SetLineColor(kGray+2);smpyfvtxmutr_h->Draw();
     pl3->cd(6);smpzfvtxmutr_h->SetFillColor(kGreen+2);smpzfvtxmutr_h->SetFillStyle(3644);smpzfvtxmutr_h->SetLineColor(kGray+2);smpxfvtxmutr_h->Draw();
     pl3->cd(7);smdphifvtx_h->SetFillColor(kGreen+2);smdphifvtx_h->SetFillStyle(3644);smdphifvtx_h->SetLineColor(kGray+2);smdphifvtx_h->Draw();
     pl3->cd(8);smdrfvtx_h->SetFillColor(kGreen+2);smdrfvtx_h->SetFillStyle(3644);smdrfvtx_h->SetLineColor(kGray+2);smdrfvtx_h->Draw();
     pl3->cd(9);smdthetafvtx_h->SetFillColor(kGreen+2);smdthetafvtx_h->SetFillStyle(3644);smdthetafvtx_h->SetLineColor(kGray+2);smdthetafvtx_h->Draw();
     pl3->cd(10);smchi2fvtx_h->SetFillColor(kGreen+2);smchi2fvtx_h->SetFillStyle(3644);smchi2fvtx_h->SetLineColor(kGray+2);smchi2fvtx_h->Draw();
     pl3->cd(11);smclusterssize_h->SetFillColor(kGreen+2);smclusterssize_h->SetFillStyle(3644);smclusterssize_h->SetLineColor(kGray+2);smclusterssize_h->Draw();
     pl3->cd(12);smdcaphi_h->SetFillColor(kGreen+2);smdcaphi_h->SetFillStyle(3644);smdcaphi_h->SetLineColor(kGray+2);smdcaphi_h->Draw();
     pl3->cd(13);smvtxindex_h->SetFillColor(kGreen+2);smvtxindex_h->SetFillStyle(3644);smvtxindex_h->SetLineColor(kGray+2);smvtxindex_h->Draw();
     pl3->cd(14);smfvtxtrackid_h->SetFillColor(kGreen+2);smfvtxtrackid_h->SetFillStyle(3644);smfvtxtrackid_h->SetLineColor(kGray+2);smfvtxtrackid_h->Draw();
     pl3->cd(15);smhitpattern_h->SetFillColor(kGreen+2);smhitpattern_h->SetFillStyle(3644);smhitpattern_h->SetLineColor(kGray+2);smhitpattern_h->Draw();
     pl3->SaveAs("/direct/phenix+u/alibordi/hf_outputs/Variable_WOSelection_set3.pdf","pdf");

*/


}

int FindGroup(int pattern)
{
   
   
   if( pattern == 1 || pattern == 4 || pattern == 5 || pattern == 8 )
      return 0;
   else if( pattern == 2 || pattern == 3 || pattern == 6 || pattern == 7 )
      return 1;
   else if( pattern == 21 || pattern == 24 || pattern == 25 || pattern == 28 )
      return 2;
   else if( pattern == 22 || pattern == 23 || pattern == 26 || pattern == 27 )
      return 3;
   else
      return -1;
}

void DefineBasePatterns()
{
   
   base1.clear();
   base2.clear();
   base3.clear();
   base4.clear();
   base1a.clear();
   base2a.clear();
   base3a.clear();
   base4a.clear();
   
   int arr_base1[16] = {1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1};
   int arr_base2[16] = {-1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1};
   int arr_base3[16] = {1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1};
   int arr_base4[16] = {-1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1};
   int arr_base1a[16] = {1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1};
   int arr_base2a[16] = {-1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1};
   int arr_base3a[16] = {-1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1};
   int arr_base4a[16] = {1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1};
   
   for(int i=0; i<16; i++)
      {
      base1.push_back(arr_base1[i]);
      base2.push_back(arr_base2[i]);
      base3.push_back(arr_base3[i]);
      base4.push_back(arr_base4[i]);
      base1a.push_back(arr_base1a[i]);
      base2a.push_back(arr_base2a[i]);
      base3a.push_back(arr_base3a[i]);
      base4a.push_back(arr_base4a[i]);
      }
   
}

int FindPattern(vector<int> &ptb, vector<int> &pty)
{
   
   
   if( ptb == base1 && pty == base3 )
      return 1;
   else if( ptb == base2 && pty == base3 )
      return 2;
   else if( ptb == base1 && pty == base4 )
      return 3;
   else if( ptb == base2 && pty == base4 )
      return 4;
   else if( ptb == base3 && pty == base1 )
      return 5;
   else if( ptb == base3 && pty == base2 )
      return 6;
   else if( ptb == base4 && pty == base1 )
      return 7;
   else if( ptb == base4 && pty == base2 )
      return 8;
   else if( ptb == base1a && pty == base3a )
      return 21;
   else if( ptb == base2a && pty == base3a )
      return 22;
   else if( ptb == base1a && pty == base4a )
      return 23;
   else if( ptb == base2a && pty == base4a )
      return 24;
   else if( ptb == base3a && pty == base1a )
      return 25;
   else if( ptb == base3a && pty == base2a )
      return 26;
   else if( ptb == base4a && pty == base1a )
      return 27;
   else if( ptb == base4a && pty == base2a )
      return 28;
   else
      return -1; 
   
}



void InitTreeVars()
{

  default_qa = -999;
  runnumber = -999;
  fillnumber = -999;
  xingshift = -999;
  badrun_flag = -999;
  b_pol = -999;
  y_pol = -999;
  b_stat = -999;
  y_stat= -999;
  b_syst = -999;
  y_syst = -999;

  for(int i=0; i<NB; i++)
    {
      patternblue[i] = -999;
      patternyell[i] = -999;
      badbunch_qa[i] = -999;
      scaler_bbcvtxcut[i] = -999;
      scaler_bbcnovtx[i] = -999;
      scaler_zdcwide[i] = -999;
      scaler_zdcnarrow[i] = -999;
    }
}
