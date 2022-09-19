#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TSpline.h>
#include <TH2.h>
#include <TF1.h>
#include <TH3.h>
#include <TCutG.h>
#include <iostream>
#include <fstream>
#include <TRandom.h>
#include <TRandom3.h>
#include <iomanip>
#include <TMath.h>
#include <TSystem.h>
#include <TVector3.h>
//#include "MuidEffGeom.h"

using namespace std;

class run_cuts{

public:

	TChain		*fChain;

	//Variables
	int           _arm;
	float         _pT;
	float					_pz;
	float         _eta;
	int           _trhits;

	int           _idhits;
	float         _DG0;
	float         _DDG0;
	float         _trchi2;
	float         _idchi2;

	int           _ntrhits;
	int           _nidhits;
	int           _lastgap;
	float         _xst1;
	float         _xst2;

	float         _xst3;
	float         _yst1;
	float         _yst2;
	float         _yst3;
	float					_st1px;

	float					_st1py;
	float					_st1pz;
	float					_phist1;
	float					_phist2;
	float					_phist3;

	float 				_radst1;
	float 				_radst2;
	float 				_radst3;
	int						_halfoctst1;
	int						_halfoctst2;

	int						_halfoctst3;
	int           _charge;
	float         _pdtheta;
	float         _bbcZ;
	int           _runnumber;

	float					_rf_vtxchi2pdf;
	float         _rf_refrad;
	float         _rf_slope;
	float         _rf_gap0x;
	float         _rf_gap0y;

	float         _rf_gap1x;
	float         _rf_gap1y;
	float         _rf_gap2x;
	float         _rf_gap2y;
	float         _rf_gap3x;

	float         _rf_gap3y;
	float         _rf_gap4x;
	float         _rf_gap4y;

	float         _mc_g_pT;
	float         _mc_g_eta;
	float         _mc_g_z;
	int           _mc_g_pid;

	float         _mc_g_px;
	float         _mc_g_py;
	float         _mc_g_pz;

	float         _mc_z;
	int						_mc_pid;
	int						_mc_hits_mutr_true;
	int						_mc_hits_muid_true;

	int						_muid_panel;
	int 					_scaled_trigbit;

	int						_trigcount;

	bool 					_scaled_MUID1DN;
	bool 					_scaled_MUID1DS;
	bool					_scaled_SG3_MUID1DHN;
	bool 					_scaled_SG3_MUID1DHS;

	run_cuts();
	virtual ~run_cuts();
	virtual void InitTree();
	virtual void SetDataset(char *newdataset);
	virtual void SetTriggerType(char *newtriggertype);
	virtual void SetFileType(char *infiletype);
	virtual void SetKinematicCuts(char *Cuts);
	virtual void SetDefaultParameters();
	virtual void SetEtaCut(float eta_low, float eta_high);
	virtual void SetZCut(float z_near, float z_far);
	virtual void SetPzCut();
	virtual void SetpTArrays();
	virtual void SetpTCorrection(int powercorrection);
	virtual void SetKaonCorrection(float kaoncorrection);
	virtual void UsePhiCuts(bool phicuts);
	virtual void UseMutrRadCuts(bool radcuts);
	virtual void PrintTuneWeights(bool printweights);
	virtual void InitializeCuts();
	//virtual void Loop(MuidEffGeom *muigeom);
	virtual void Loop();
	virtual void BookHistos();
	virtual void Load_tuning_weights();
	virtual void FillGap23(float weight, bool cocktail); 
	virtual void FillDiagnosticHistos(float weight);
	//virtual void FillHVGroupEff(float weight);
	virtual bool CheckPhiCut();
	virtual bool CheckRadCut();

private:

	std::string dataset;
	std::string filetype;
	std::string cuts;
	std::string triggertype;
	std::string chargetype;
	std::string armtype;

	int pTweighting;
	float KaonWeight;

	bool MANUAL_PHI_CUTS;
	bool MANUAL_RAD_CUTS;

	enum { narm = 2};
	enum { ncharge = 2};
	enum { ngap = 5};
	//enum { nptbin = 17};
	//enum { nmcptbin = 25};
	enum { nptbin = 14};
	enum { nmcptbin = 23};
	enum { nsta = 3};
	enum { npanel = 6};
	enum { norient = 2};
	enum { nhvgroup = 30};

	float cut_eta[2];
	float cut_z[2];
	float cut_Pz[narm][2];
	float cut_dPz[narm][2];

	double varbin_pTarray[nptbin+1];
	double varbin_p_array[nptbin+1];
	double varbin_mc_pTarray[nmcptbin+1];

	double varbin_AN_pTarray[10];
	double varbin_AN_pzarray[5];

	//cut variable arrays, assigned by data set
	float cut_mutr_chi2[narm][ngap];
	float cut_pdtheta[narm][ngap];

	float cut_pdtheta__;
	float cut_mutr_chi2__;
	float cut_muid_chi2__;
	float cut_road_slope__;
	float cut_dg0__;
	float cut_ddg0__;
	float cut_vertex_chi2__;
	float cut_vertex_rad__;

	TFile *cut_file;
	TF1 *fcut_vertex_chi2[narm][ngap];
	TF1 *fcut_dg0[narm][ngap];
	TF1 *fcut_ddg0[narm][ngap];
	TF1 *fcut_vertex_rad[narm][ngap];

	TFile *fSD;
	TH1F *hSD_MUID1D[narm];
	TH1F *hSD_SG3_MUID1DH[narm];

	TH1F *hcut_vertex_chi2[narm][ngap];
	TH1F *hcut_dg0[narm][ngap];
	TH1F *hcut_ddg0[narm][ngap];
	TH1F *hcut_vertex_rad[narm][ngap];

	float cut_dg0[narm][ngap][nptbin];
	float cut_ddg0[narm][ngap][nptbin];
	float cut_vertex_rad[narm][ngap][nptbin];
	float cut_vertex_chi2[narm][ngap][nptbin];
	
	//! arrays for reweighting
	float tune_weight[narm][ncharge][nmcptbin];

	//! Cut variable histograms
	TH2F *DG0[narm][ngap]; //DG0
	TH2F *DDG0[narm][ngap]; //DDG0
	TH2F *REFRAD[narm][ngap]; //Refrad
	TH2F *MUTR_CHI2[narm][ngap]; //MuTR chi2

	TH2F *VTX_CHI2[narm][ngap]; //vtx chi2
	TH2F *VTX_CHI2_OUT[narm][ngap]; //vtx chi2

	TH2F *PDTHETA[narm][ngap]; //pdtheta
	TH2F *PDTHETA_FAKE[narm][ngap]; //pdtheta
	TH2F *PDTHETA_OUT[narm][ngap]; //pdtheta

	// all cuts
	TH2F *tuning_matrix_allcuts[narm][ngap][ncharge];

	//for cocktail
	TH1F *Ndecay[narm][ngap][ncharge];
	TH1F *Ndecayabs[narm][ngap][ncharge];
	TH1F *Npunchthrough[narm][ngap][ncharge];

	TH1F *Ndecay_AN_pT[narm][ngap][ncharge];
	TH1F *Ndecayabs_AN_pT[narm][ngap][ncharge];
	TH1F *Npunchthrough_AN_pT[narm][ngap][ncharge];

	TH1F *Ndecay_AN_pz[narm][ngap][ncharge];
	TH1F *Ndecayabs_AN_pz[narm][ngap][ncharge];
	TH1F *Npunchthrough_AN_pz[narm][ngap][ncharge];

	//gap4
	//TH1F *inclZ[narm][ncharge][nptbin];
	//TH1F *inclZ_decay[narm][ncharge][nptbin];
	//TH1F *inclZ_punchthrough[narm][ncharge][nptbin];

	TH2F *inclZ[narm][ncharge];
	TH2F *inclZ_decay[narm][ncharge];
	TH2F *inclZ_punchthrough[narm][ncharge];

	//pT distribution
	TH1F *N_varbin[narm][ngap][ncharge];
	TH1F *N_varbin_AN_pT[narm][ngap][ncharge];
	TH1F *N_varbin_AN_pz[narm][ngap][ncharge];
	TH1F *N_varbin_AN_pz_MUID1D[narm][ngap][ncharge];
	TH1F *N_varbin_AN_pz_SG3MUID1DH[narm][ngap][ncharge];
	TH1F *N_varbin_fake[narm][ngap][ncharge];
	TProfile *avg_pT_varbin_gap[narm][ngap][ncharge];
	TProfile *avg_pT_varbin[narm][ncharge];
	TProfile *avg_p[narm][ngap][ncharge];
	TProfile *avg_pT_AN_pz[narm][ncharge];

	TH2F *N_pT_pz[narm][ngap][ngap];

	TH1F *N_varbin_pion[narm][ngap][ncharge];
	TH1F *N_varbin_kaon[narm][ngap][ncharge];
	TH1F *N_varbin_kaon0[narm][ngap][ncharge];
	TH1F *N_varbin_proton[narm][ngap][ncharge];

	TH2F *YvX_mutr_cut[narm][nsta]; //MuTr station
	TH2F *Rphi_mutr_cut[narm][nsta]; //MuTr station

	//MuidEffGeom *muigeom;
	

	//trigger eff
	TFile *ftrig;
	TF1 *ftrig_SG3_pT[narm][ngap];
	TF1 *ftrig_MUID1D_pT[narm];
	
	//charge weight
	TFile *fchg_wt;
	TH1F *hW_KAON[ncharge][10];

};

run_cuts::run_cuts(){

}

run_cuts::~run_cuts(){

	if ( !fChain ) return;
	delete fChain->GetCurrentFile();
}

void run_cuts::SetDataset(char *newdataset){
	cout << "SetDataset" << endl;
	dataset = newdataset;
}

void run_cuts::SetTriggerType(char *newtriggertype){
	cout << "SetTriggerType" << endl;
	triggertype = newtriggertype;
}

void run_cuts::SetFileType(char *infiletype){
	cout << "SetFileType" << endl;
	filetype = infiletype;
	//InitTree();

	SetDefaultParameters();
}

void run_cuts::SetKinematicCuts(char *Cuts){
	cuts = Cuts;
}

void run_cuts::SetDefaultParameters(){
	cout << "SetDefaultParameter" << endl;

	SetEtaCut(1.2, 2.2);
	SetZCut(-30, 30);
	SetPzCut();
	SetpTArrays();
	SetpTCorrection(4);
	SetKaonCorrection(1);
}

void run_cuts::SetEtaCut(float eta_low, float eta_high){
	cut_eta[0] = fabs(eta_low);
	cut_eta[1] = fabs(eta_high);
}

void run_cuts::SetZCut(float z_near, float z_far){
	cut_z[0] = z_near;
	cut_z[1] = z_far;
}

void run_cuts::SetPzCut(){
	//cut_Pz[0][0] = 2.4; //south gap2
	//cut_Pz[0][1] = 2.6; //south gap3
	//cut_Pz[1][0] = 2.6; //north gap2
	//cut_Pz[1][1] = 2.8; //north gap3

	cut_Pz[0][0] = 3.40; //south gap2
	cut_Pz[0][1] = 3.60; //south gap3
	cut_Pz[1][0] = 3.80; //north gap2
	cut_Pz[1][1] = 4.00; //north gap3

	cut_dPz[0][0] = 0.20;
	cut_dPz[0][1] = 0.20;
	cut_dPz[1][0] = 0.20;
	cut_dPz[1][1] = 0.20;
}

void run_cuts::SetpTArrays(){

	/*
	float tmp_pTarray[nptbin+1] = {1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 10.0, 12.0};
	float tmp_p_array[nptbin+1] = {3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 7.00, 8.00, 9.00, 10.0, 12.0, 14.0, 16.0, 19.0, 22.0, 25.0, 30.0};
	float tmp_mc_pTarray[nmcptbin+1] = {0.80, 1.00, 1.20, 1.40, 1.60,
																			1.80, 2.00, 2.20, 2.40, 2.60,
																			2.80, 3.00, 3.25, 3.50, 3.75,
																			4.00, 4.50, 5.00, 5.50, 6.00,
																			7.00, 8.00, 9.00, 10.0, 12.0,
																			15.0};
																			*/

	float tmp_pTarray[nptbin+1] = {1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00};
	float tmp_p_array[nptbin+1] = {3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 7.00, 8.00, 9.00, 10.0, 12.0, 14.0, 17.0, 20.0};
	float tmp_mc_pTarray[nmcptbin+1] = {0.80, 1.00, 1.20, 1.40, 1.60,
																			1.80, 2.00, 2.20, 2.40, 2.60,
																			2.80, 3.00, 3.25, 3.50, 3.75,
																			4.00, 4.50, 5.00, 5.50, 6.00,
																			7.00, 8.00, 9.00, 10.0
																			};


	for ( int ipt=0; ipt<nptbin+1; ipt++){
		varbin_pTarray[ipt] = tmp_pTarray[ipt];
		varbin_p_array[ipt] = tmp_p_array[ipt];
	}

	for ( int ipt=0; ipt<nmcptbin+1; ipt++){
		varbin_mc_pTarray[ipt] = tmp_mc_pTarray[ipt];
	}

	float tmp_AN_pTarray[10] = {1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0};
	float tmp_AN_pzarray[5] = {3.5, 5.0, 7.0, 10.0, 20.0};

	for (int ipt=0; ipt<10; ipt++) varbin_AN_pTarray[ipt] = tmp_AN_pTarray[ipt];
	for (int ipz=0; ipz<5; ipz++) varbin_AN_pzarray[ipz] = tmp_AN_pzarray[ipz];

}

void run_cuts::SetpTCorrection(int powercorrection){
	if ( powercorrection==1 ) pTweighting = 1;
	else if ( powercorrection==3 ) pTweighting = 3;
	else if ( powercorrection==4 ) pTweighting = 4;
	else if ( powercorrection==5 ) pTweighting = 5;
	else{
		cout << "In appropriate power correciton of : pT^ " << powercorrection << endl;
		exit(1);
	}
}

void run_cuts::SetKaonCorrection(float kaoncorrection){
	KaonWeight = kaoncorrection;
	if ( KaonWeight!=1 ){
		cout << "Set Kaon Correction: " << KaonWeight << endl;
	}
}

void run_cuts::UsePhiCuts(bool phicuts){
	if ( phicuts ){
		MANUAL_PHI_CUTS = true;
		cout << "Employing use defined phi cuts to improve MC/DATA matching" << endl;
	}else{
		MANUAL_PHI_CUTS = false;
		cout << "Manual Phi cuts disabled" << endl;
	}
}

void run_cuts::UseMutrRadCuts(bool radcuts){
	if ( radcuts ){
		MANUAL_RAD_CUTS = true;
		cout << "Employing use defined Mutr Rad cuts to improve MC/DATA matching" << endl;
	}else{
		MANUAL_RAD_CUTS = false;
		cout << "Manual Mutr Rad cuts disabled" << endl;
	}
}

void run_cuts::PrintTuneWeights(bool printweights){
	if ( printweights ){
		for (int jpt=0; jpt<nmcptbin; jpt++){
			cout << "MCpT bin : " << jpt << " weights : ";
			cout << " S.neg : "<< setw(6) << setprecision(4) << tune_weight[0][0][jpt];
			cout << " N.neg : "<< setw(6) << setprecision(4) << tune_weight[1][0][jpt];
			cout << " S.pos : "<< setw(6) << setprecision(4) << tune_weight[0][1][jpt];
			cout << " N.pos : "<< setw(6) << setprecision(4) << tune_weight[1][1][jpt];
			cout << endl;
		}
	}
}

void run_cuts::InitTree(){

	string fname;
	ifstream flist;

	fChain = new TChain("ana_tree");

	if ( filetype=="QGSP_BERT" || filetype=="FTFP_BERT" || filetype=="QGSP_BIC" || filetype=="COCKTAIL" ){
		//flist.open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/15.inclusiveHF/list_Run15pp200_hadron_cocktail_QGSP_BERT.lst");
		sprintf((char*)fname.c_str(),"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/15.inclusiveHF/list_Run15pp200_hadron_cocktail_%s.lst",filetype.data());
		cout << "OPEN femtotDST list : " << fname.data() << endl;
		flist.open(fname.data());
		while ( flist >> fname ){
			cout << "OPEN femtoDST named : " << fname.data() << endl;
			fChain->AddFile(fname.data());
		}
	}else if ( filetype=="DATA" ){

		int runnum;
		//flist.open("/direct/phenix+hhj/shlim/work/15.run15/10.runQA/04.goodrun/good_SG3_MUID1DH_BOTH_v0.lst");
		flist.open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/10.runQA/03.yield_run/Run15pp200_goodrun_inclHF_20170301.lst");
		//flist.open("test.lst");

		int _count = 0;
		while ( flist >> runnum ){
			if ( triggertype=="MB" ){
				sprintf((char*)fname.c_str(),"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/06.femtoDST/Run15pp200_femtoDST_MB/Run15pp200MB_femtoDST_%d.root",runnum);
			}else{
				sprintf((char*)fname.c_str(),"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/06.femtoDST/jobdir_Run15pp200MU_pro108_sngmu_00/Run15pp200_pro108_femtoDST_sngmu_%d.root",runnum);
			}

			if ( _count<10 ){
				cout << "#" << _count << ", OPEN femtoDST named : " << fname.data() << endl;
			}

			fChain->AddFile(fname.data());
			_count++;
		}//while

		cout << "Number of good runs: " << _count << endl;

	}else if ( filetype=="MUON" ){
		//fname = "/direct/phenix+hhj/shlim/work/09.run12_pdst/10.pp200_check/Run12pp200_muon_femtoDST.root";
		//fname = "femtoDST.root";
		flist.open("femtodst_muon.lst");
		while ( flist >> fname ){
			cout << "OPEN femtoDST named : " << fname.data() << endl;
			fChain->AddFile(fname.data());
		}
	}//MUON

	//TFile *fin = new TFile(tmp_fname, "read");
	//fChain = (TTree*)fin->Get("ana_tree");

	fChain->SetBranchAddress("_arm",&_arm);
	fChain->SetBranchAddress("_pT",&_pT);
	fChain->SetBranchAddress("_pz",&_pz);
	fChain->SetBranchAddress("_eta",&_eta);
	fChain->SetBranchAddress("_trhits",&_trhits);

	fChain->SetBranchAddress("_idhits",&_idhits);
	fChain->SetBranchAddress("_DG0",&_DG0);
	fChain->SetBranchAddress("_DDG0",&_DDG0);
	fChain->SetBranchAddress("_trchi2",&_trchi2);
	fChain->SetBranchAddress("_idchi2",&_idchi2);

	fChain->SetBranchAddress("_ntrhits",&_ntrhits);
	fChain->SetBranchAddress("_nidhits",&_nidhits);
	fChain->SetBranchAddress("_lastgap",&_lastgap);
	fChain->SetBranchAddress("_xst1",&_xst1);
	fChain->SetBranchAddress("_xst2",&_xst2);

	fChain->SetBranchAddress("_xst3",&_xst3);
	fChain->SetBranchAddress("_yst1",&_yst1);
	fChain->SetBranchAddress("_yst2",&_yst2);
	fChain->SetBranchAddress("_yst3",&_yst3);
	fChain->SetBranchAddress("_st1px",&_st1px);

	fChain->SetBranchAddress("_st1py",&_st1py);
	fChain->SetBranchAddress("_st1pz",&_st1pz);
	//fChain->SetBranchAddress("_phist1",&_phist1);
	//fChain->SetBranchAddress("_phist2",&_phist2);
	//fChain->SetBranchAddress("_phist3",&_phist3);

	//fChain->SetBranchAddress("_radst1",&_radst1);
	//fChain->SetBranchAddress("_radst2",&_radst2);
	//fChain->SetBranchAddress("_radst3",&_radst3);
	//fChain->SetBranchAddress("_halfoctst1",&_halfoctst1);
	//fChain->SetBranchAddress("_halfoctst2",&_halfoctst2);

	//fChain->SetBranchAddress("_halfoctst3",&_halfoctst3);
	fChain->SetBranchAddress("_charge",&_charge);
	fChain->SetBranchAddress("_pdtheta",&_pdtheta);
	fChain->SetBranchAddress("_bbcZ",&_bbcZ);
	fChain->SetBranchAddress("_runnumber",&_runnumber);
	fChain->SetBranchAddress("_scaled_trigbit",&_scaled_trigbit);

	fChain->SetBranchAddress("_rf_vtxchi2pdf",&_rf_vtxchi2pdf);
	fChain->SetBranchAddress("_rf_refrad",&_rf_refrad);
	fChain->SetBranchAddress("_rf_slope",&_rf_slope);
	//fChain->SetBranchAddress("_rf_gap0x",&_rf_gap0x);
	//fChain->SetBranchAddress("_rf_gap0y",&_rf_gap0y);

	//fChain->SetBranchAddress("_rf_gap1x",&_rf_gap1x);
	//fChain->SetBranchAddress("_rf_gap1y",&_rf_gap1y);
	//fChain->SetBranchAddress("_rf_gap2x",&_rf_gap2x);
	//fChain->SetBranchAddress("_rf_gap2y",&_rf_gap2y);
	//fChain->SetBranchAddress("_rf_gap3x",&_rf_gap3x);

	//fChain->SetBranchAddress("_rf_gap3y",&_rf_gap3y);
	//fChain->SetBranchAddress("_rf_gap4x",&_rf_gap4x);
	//fChain->SetBranchAddress("_rf_gap4y",&_rf_gap4y);

	if ( filetype!="DATA" ){
		fChain->SetBranchAddress("_mc_g_pT",&_mc_g_pT);
		fChain->SetBranchAddress("_mc_g_z",&_mc_g_z);
		fChain->SetBranchAddress("_mc_g_pid",&_mc_g_pid);

		fChain->SetBranchAddress("_mc_g_px",&_mc_g_px);
		fChain->SetBranchAddress("_mc_g_py",&_mc_g_py);
		fChain->SetBranchAddress("_mc_g_pz",&_mc_g_pz);

		fChain->SetBranchAddress("_mc_z",&_mc_z);
		fChain->SetBranchAddress("_mc_pid",&_mc_pid);

		fChain->SetBranchAddress("_mc_hits_mutr_true",&_mc_hits_mutr_true);
		fChain->SetBranchAddress("_mc_hits_muid_true",&_mc_hits_muid_true);
	}

	cout << "DONE InitTree()" << endl;
}

void run_cuts::FillDiagnosticHistos(float weight){

	if ( ((_lastgap==3 || _lastgap==2) && fabs(_pz)>(cut_Pz[_arm][_lastgap-2]+1.0*(_pT-1.0))) || (_lastgap==4) ) {

		//float ptot = sqrt(_pT*_pT + _pz*_pz);

		//DG0-DDG0
		if ( 
				_rf_slope>cut_road_slope__ //GAP0SLOPE
				&& _rf_refrad<cut_vertex_rad__ //REFRAD
				&& _pdtheta<cut_pdtheta__ //PDTHETA
				&& _trchi2<cut_mutr_chi2__ //MUTRCHI2
				&& _rf_vtxchi2pdf<cut_vertex_chi2__ //VTXCHI2
				){
			DG0[_arm][_lastgap]->Fill(_pT, _DG0, weight);
			DDG0[_arm][_lastgap]->Fill(_pT, _DDG0, weight);
		}

		//REFRAD
		if ( 
				_rf_slope>cut_road_slope__ //GAP0SLOPE
				&& _DDG0<cut_ddg0__ //DDG0
				&& _DG0<cut_dg0__ //DG0
				&& _pdtheta<cut_pdtheta__ //PDTHETA
				&& _trchi2<cut_mutr_chi2__ //MUTRCHI2
				&& _rf_vtxchi2pdf<cut_vertex_chi2__ //VTXCHI2
				){
			REFRAD[_arm][_lastgap]->Fill(_pT, _rf_refrad, weight);
		}

		//MUTR_CHI2 && MUID_CHI2
		if ( 
				_rf_slope>cut_road_slope__ //GAP0SLOPE
				&& _DDG0<cut_ddg0__ //DDG0
				&& _DG0<cut_dg0__ //DG0
				&& _pdtheta<cut_pdtheta__ //PDTHETA
				&& _rf_refrad<cut_vertex_rad__ //REFRAD
				&& _rf_vtxchi2pdf<cut_vertex_chi2__ //VTXCHI2
				){
			MUTR_CHI2[_arm][_lastgap]->Fill(_pT, _trchi2, weight);
		}

		//VTXCHI2
		if ( 
				_rf_slope>cut_road_slope__ //GAP0SLOPE
				&& _DDG0<cut_ddg0__ //DDG0
				&& _DG0<cut_dg0__ //DG0
				&& _pdtheta<cut_pdtheta__ //PDTHETA
				&& _trchi2<cut_mutr_chi2__ //MUTRCHI2
				&& _rf_refrad<cut_vertex_rad__ //REFRAD
				){
				VTX_CHI2[_arm][_lastgap]->Fill(_pT, _rf_vtxchi2pdf, weight);

				if ( filetype!="DATA" && _mc_g_pT<0.80*_pT ){
					VTX_CHI2_OUT[_arm][_lastgap]->Fill(_pT, _rf_vtxchi2pdf, weight);
				}
		}

		//PDTHETA
		if ( 
				_rf_slope>cut_road_slope__ //GAP0SLOPE
				&& _DDG0<cut_ddg0__ //DDG0
				&& _DG0<cut_dg0__ //DG0
				&& _trchi2<cut_mutr_chi2__ //MUTRCHI2
				&& _rf_refrad<cut_vertex_rad__  //REFRAD
				&& _rf_vtxchi2pdf<2.0*cut_vertex_chi2__ //VTXCHI2
			 ){
			if ( _rf_vtxchi2pdf<cut_vertex_chi2__ ){
				PDTHETA[_arm][_lastgap]->Fill(_pT, _pdtheta, weight);
				if ( filetype!="DATA" && _mc_g_pT<0.80*_pT ){
					PDTHETA_FAKE[_arm][_lastgap]->Fill(_pT, _pdtheta, weight);
				}
			}else{
				PDTHETA_OUT[_arm][_lastgap]->Fill(_pT, _pdtheta, weight);
			}
		}
	}
}

//void run_cuts::FillGap23(int ptbin,float weight,bool cocktail){
void run_cuts::FillGap23(float weight,bool cocktail){

	if ( 
			fabs(_pz)>(cut_Pz[_arm][_lastgap-2]+1.0*(_pT-1.0))
			&& _rf_slope>cut_road_slope__ //1
			&& _rf_refrad<cut_vertex_rad__ //2
			&& _DG0<cut_dg0__ //3
			&& _DDG0<cut_ddg0__ //4
			&& _trchi2<cut_mutr_chi2__ //5
			&& _pdtheta<cut_pdtheta__ //6
			&& _rf_vtxchi2pdf<cut_vertex_chi2__ //7
			){

		avg_pT_varbin[_arm][_charge]->Fill(_pT,_pT,weight);
		avg_pT_varbin_gap[_arm][_lastgap][_charge]->Fill(_pT,_pT,weight);
		avg_p[_arm][_lastgap][_charge]->Fill(_pT,sqrt(_pT*_pT+_pz*_pz),weight);

		N_varbin[_arm][_lastgap][_charge]->Fill(_pT,weight);

		if ( filetype!="DATA" && _mc_g_pT<0.8*_pT ){
			N_varbin_fake[_arm][_lastgap][_charge]->Fill(_pT, weight);
		}

		if ( cocktail ){
			if ( abs(_mc_pid)==13 ){
				if ( fabs(_mc_z)<80.0 ){
					Ndecay[_arm][_lastgap][_charge]->Fill(_pT,weight);
					if ( _pT>1.25 && fabs(_pz)>3.50 ){
						Ndecay_AN_pT[_arm][_lastgap][_charge]->Fill(_pT,weight);
						Ndecay_AN_pz[_arm][_lastgap][_charge]->Fill(fabs(_pz),weight);
					}
				}else{
					Ndecayabs[_arm][_lastgap][_charge]->Fill(_pT,weight);
					if ( _pT>1.25 && fabs(_pz)>3.50 ){
						Ndecayabs_AN_pT[_arm][_lastgap][_charge]->Fill(_pT,weight);
						Ndecayabs_AN_pz[_arm][_lastgap][_charge]->Fill(fabs(_pz),weight);
					}
				}
			}else{
				Npunchthrough[_arm][_lastgap][_charge]->Fill(_pT,weight);
				if ( _pT>1.25 && fabs(_pz)>3.50 ){
					Npunchthrough_AN_pT[_arm][_lastgap][_charge]->Fill(_pT,weight);
					Npunchthrough_AN_pz[_arm][_lastgap][_charge]->Fill(fabs(_pz),weight);
				}
			}

			if ( (_pT-_mc_g_pT)<(0.5*_pT) ){
				tuning_matrix_allcuts[_arm][_lastgap][_charge]->Fill(_mc_g_pT,_pT,weight);
			}

			if ( abs(_mc_g_pid)==211 ){
				N_varbin_pion[_arm][_lastgap][_charge]->Fill(_pT,weight);
			}else if ( abs(_mc_g_pid)==321 ){
				N_varbin_kaon[_arm][_lastgap][_charge]->Fill(_pT,weight);
			}else if ( abs(_mc_g_pid)==2212 ){
				N_varbin_proton[_arm][_lastgap][_charge]->Fill(_pT,weight);
			}

		}//cocktail
	}//quality cuts
}

