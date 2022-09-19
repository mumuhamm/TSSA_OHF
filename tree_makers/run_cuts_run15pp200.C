#include "run_cuts_run15pp200.h"
#include <TStyle.h>
#include <TProfile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <iostream>

#include "/direct/phenix+u/shlim/RunRange.h"

//void run_cuts::Loop(MuidEffGeom *muigeom){
void run_cuts::Loop(){

	InitializeCuts();
	InitTree();

	if ( pTweighting!=1 && pTweighting!=3 && pTweighting!=4 && pTweighting!=5 ){
		cout << "pT weighting scheme not specified, exiting." << endl;
		exit(1);
	}

	cout << "Power correction of pT^ " << pTweighting << " specified." << endl;

	bool COCKTAIL = kFALSE;
	bool DATA = kFALSE;
	bool MUON = kFALSE;
	bool PTCORRECTION = kFALSE;

	cout << "Starting Loop()" << endl;
	cout << "Eta bin from " << cut_eta[0] << " to " << cut_eta[1] << endl;

	cout << " Z cuts set to near-side: " << cut_z[0] 
		 << " and to far-side: " << cut_z[1]
		 << " North arm convention." << endl;



	if ( filetype=="QGSP_BERT" || filetype=="FTFP_BERT" || filetype=="QGSP_BIC" || filetype=="COCKTAIL" ) COCKTAIL = kTRUE;
	if ( filetype=="MUON" ) MUON = kTRUE;
	if ( filetype=="DATA" ) DATA = kTRUE;

	if ( (DATA && COCKTAIL) || (DATA && MUON) ){
		cout << "Data and COCKTAIL are kTRUE simultaneously, exiting ... " << endl;
		exit(1);
	}

	if ( (COCKTAIL) && ( pTweighting==3 || pTweighting==4 || pTweighting==5 ) ){
		cout << "pTweighting = " << pTweighting << "\tDo pT correction!!" << endl;
		PTCORRECTION=kTRUE;
	}

	if ( triggertype=="MUID1D" || triggertype=="SG3_MUID1DH" || triggertype=="MB" || triggertype=="COMBINED" || triggertype=="OTHERS" ){
		cout << "Choose Trigger Type : " << triggertype << endl;
	}else{
		cout << "Choose trigger type, BBCLL1_NOVTX and MUIDLL1, now : " << triggertype << endl;
		exit(1);
	}

	cout << "Kinematic cuts are set to - " << cuts << endl;

	char fname[300];
	sprintf(fname, "runcuts_%s_%s_%s_%s.root",dataset.c_str(),triggertype.c_str(),filetype.c_str(),cuts.c_str());

	TFile *outfile = new TFile(fname, "recreate");
	cout << "OPEN OUTPUT FILE : " << fname << endl;

	BookHistos();
	cout << "DONE BOOK HISTOS!!" << endl;

	Load_tuning_weights();
	if ( COCKTAIL ){
		PrintTuneWeights(1);
	}

	int kpt = -999, mpt = -999;
	float weight_by_thrown_pt = 0.0;

	_trigcount = 0;

	if ( fChain==0 ) return;
	int nentries = fChain->GetEntries();

	cout << nentries << " tracks to be analyzed.." << endl;

	bool _goodflag = kFALSE;

	for (int jentry=0; jentry<nentries; jentry++){

		fChain->GetEntry(jentry);

		if ( (jentry%(nentries/10))==0 ){
			float percent = jentry * 100.0 / nentries;
			//cout << jentry << "\t" << nentries/20 << endl;
			cout << setw(12) << jentry << "\tevnets,\t" << setw(12) << percent << "\tpercent completed..." << endl;
		}

		kpt = -999;
		mpt = -999;

		//////////////////////
		//PT BIN AND ETA BIN//
		//////////////////////
		for (int jj=0; jj<nptbin; jj++){
			if ( _pT>varbin_pTarray[jj] && _pT<=varbin_pTarray[jj+1] ){
				kpt = jj;
				break;
			}
		}

		for (int jj=0; jj<nmcptbin; jj++){
			if ( _mc_g_pT>varbin_mc_pTarray[jj] && _mc_g_pT<=varbin_mc_pTarray[jj+1] ){
				mpt = jj;
				break;
			}
		}

		if ( kpt<0 || kpt>nptbin-1 || isnan(_pT) ){
			continue;
		}
		if ( COCKTAIL && (mpt<0 || mpt>nmcptbin-1) ){
			continue;
		}
		if ( (COCKTAIL||MUON) && _mc_g_pid==-999 ){
			continue;
		}

		if ( COCKTAIL || MUON ){
			if ( (_mc_hits_mutr_true*1.0/_ntrhits)<0.6 ) continue;
			//if ( (float(_mc_hits_muid_true)/_nidhits)<0.6 ) continue;
		}

		/////////////
		//BASIC CUT//
		/////////////
		//if ( COCKTAIL || MUON ) _bbcZ = _mc_g_z; 

		if ( _bbcZ>cut_z[1] || _bbcZ<cut_z[0] ) continue;
		if ( fabs(_eta)<cut_eta[0] || fabs(_eta)>cut_eta[1] ) continue;
		if ( _lastgap<1.5 ) continue;
		if ( _ntrhits<10.5 ) continue;
		if ( _lastgap<4 && fabs(_eta)<1.4 ) continue;


		_phist1 = atan2(_yst1,_xst1);
		_phist2 = atan2(_yst2,_xst2);
		_phist3 = atan2(_yst3,_xst3);
		_radst1 = sqrt(_xst1*_xst1 + _yst1*_yst1);
		_radst2 = sqrt(_xst2*_xst2 + _yst2*_yst2);
		_radst3 = sqrt(_xst3*_xst3 + _yst3*_yst3);

		/*
		if ( _phist1<0 ) _phist1 += 2*TMath::Pi();
		if ( _phist2<0 ) _phist2 += 2*TMath::Pi();
		if ( _phist3<0 ) _phist3 += 2*TMath::Pi();

		_halfoctst1 = int(_phist1/(TMath::Pi()/8.));
		_halfoctst2 = int(_phist2/(TMath::Pi()/8.));
		_halfoctst3 = int(_phist3/(TMath::Pi()/8.));
		*/

		_scaled_MUID1DN = _scaled_MUID1DS = _scaled_SG3_MUID1DHN = _scaled_SG3_MUID1DHS = false;
		if ( DATA ){
			if ( (_scaled_trigbit&0x00400000)>0 ) _scaled_MUID1DN = true;
			if ( (_scaled_trigbit&0x00800000)>0 ) _scaled_MUID1DS = true;
			if ( (_scaled_trigbit&0x01000000)>0 ) _scaled_SG3_MUID1DHN = true;
			if ( (_scaled_trigbit&0x02000000)>0 ) _scaled_SG3_MUID1DHS = true;
		}

		///////////////
		//TRIGGER CUT//
		///////////////
		if ( DATA && triggertype=="SG3_MUID1DH" ){
			if ( _arm==0 && !_scaled_SG3_MUID1DHS ) continue;
			if ( _arm==1 && !_scaled_SG3_MUID1DHN ) continue;
		}

		if ( DATA && triggertype=="MUID1D" ){
			if ( _arm==0 && !_scaled_MUID1DS ) continue;
			if ( _arm==1 && !_scaled_MUID1DN ) continue;
		}

		//Run15pAu200
		if ( DATA && triggertype=="COMBINED" ){
			if ( _arm==0 ){
				if ( _lastgap==2 || _lastgap==3 || (_lastgap==4 && _pT>2.50) ){
					if ( !_scaled_SG3_MUID1DHS ) continue;
				}else{
					if ( !_scaled_MUID1DS ) continue;
				}
			}else{
				if ( _lastgap==2 || _lastgap==3 || (_lastgap==4 && _pT>2.50) ){
					if ( !_scaled_SG3_MUID1DHN ) continue;
				}else{
					if ( !_scaled_MUID1DN ) continue;
				}
			}
		}


		////////////
		//PHI CUTS//
		////////////
		if ( MANUAL_PHI_CUTS ) _goodflag = CheckPhiCut();
		if ( MANUAL_PHI_CUTS && !_goodflag ) continue;
		/////////////////
		//MUTR RAD CUTS//
		/////////////////
		if ( MANUAL_RAD_CUTS ) _goodflag = CheckRadCut();
		if ( MANUAL_RAD_CUTS && !_goodflag ) continue;

		//////////////
		//LOOSE CUTS//
		//////////////
		//if ( COCKTAIL && cuts=="LOOSE" && (_lastgap==2 || _lastgap==3) && (fabs(_mc_z)>40.0 || abs(_mc_pid)==13) ) continue;

		/////////////////
		//PT CORRECTION//
		/////////////////
		weight_by_thrown_pt = 1.0;
		if ( COCKTAIL ){

			TVector3 vec(_mc_g_px,_mc_g_py,_mc_g_pz);
			int etabin = int((fabs(vec.Eta())-1.2)/0.2); 

			float kaon_wt = KaonWeight;
			if ( (etabin>=0 && etabin<8) && abs(_mc_g_pid)==321 ){
				kaon_wt *= hW_KAON[_charge][etabin]->GetBinContent(hW_KAON[_charge][etabin]->FindBin(_mc_g_pT));
			}

			//KAON VARIATION//
			if ( COCKTAIL && (abs(_mc_g_pid)==321 || _mc_g_pid==130 || _mc_g_pid==310) ) weight_by_thrown_pt *= kaon_wt; 
			weight_by_thrown_pt = tune_weight[_arm][_charge][mpt];

		}//COCKTAIL

		if ( (COCKTAIL) && PTCORRECTION ){
			if ( _mc_g_pT<1.0 ){
				weight_by_thrown_pt = weight_by_thrown_pt * 1.0;
			}else{
				if ( pTweighting==3 ) weight_by_thrown_pt = weight_by_thrown_pt / (_mc_g_pT*_mc_g_pT);
				if ( pTweighting==4 ) weight_by_thrown_pt = weight_by_thrown_pt / (_mc_g_pT*_mc_g_pT*_mc_g_pT);
				if ( pTweighting==5 ) weight_by_thrown_pt = weight_by_thrown_pt / (_mc_g_pT*_mc_g_pT*_mc_g_pT*_mc_g_pT);
			}
		} //COCKTAIL && PTCORRECTION

		/////////////////////////////////
		//TRIGGER EFFICIENCY CORRECTION//
		/////////////////////////////////
		if ( DATA ){
			float eff_SG3 = 1.0, eff_MUID1D = 1.0;
			float SD_SG3 = hSD_SG3_MUID1DH[_arm]->GetBinContent(hSD_SG3_MUID1DH[_arm]->FindBin(_runnumber));
			float SD_MUID1D = hSD_MUID1D[_arm]->GetBinContent(hSD_MUID1D[_arm]->FindBin(_runnumber));

			if ( _lastgap==4 ){
				eff_SG3 = ftrig_SG3_pT[_arm][_lastgap]->Eval(_pT); 
				eff_MUID1D = ftrig_MUID1D_pT[_arm]->Eval(_pT); 
			}else{
				eff_SG3 = ftrig_SG3_pT[_arm][_lastgap]->Eval(_pT); 
			}

			if ( _trigcount<20 ){
				cout << "Trig eff, arm : " << _arm << ", lastgap : " << _lastgap << ", eff(SG3_MUID1DH) : " << eff_SG3 << ", eff(MUID1D) : " << eff_MUID1D << endl;
				_trigcount++;
			}

			if ( triggertype=="SG3_MUID1DH" ){
				weight_by_thrown_pt = 1./eff_SG3*(SD_SG3+1);
			}else if ( triggertype=="MUID1D" ){
				weight_by_thrown_pt = 1./eff_MUID1D*(SD_MUID1D+1);
			}else if ( triggertype=="COMBINED" ){
				if ( _lastgap==2 || _lastgap==3 ){
					weight_by_thrown_pt = 1./eff_SG3*(SD_SG3+1);
				}else{
					if ( _pT<2.50 )
						weight_by_thrown_pt = 1./eff_MUID1D*(SD_MUID1D+1);
					else
						weight_by_thrown_pt = 1./eff_SG3*(SD_SG3+1);
				}
			}

			if ( isnan(weight_by_thrown_pt) || isinf(weight_by_thrown_pt) || weight_by_thrown_pt<0 || weight_by_thrown_pt>1000. ){
			//if ( !(weight_by_thrown_pt>0 && weight_by_thrown_pt<1000) ){
				cout << "WRONG WEIGHT: " << weight_by_thrown_pt << ", pT: " << _pT << ", eff(SG3_MUID1DH) : " << eff_SG3 << ", eff(MUID1D) : " << eff_MUID1D << endl;
				continue;
			}
		}//DATA



		//////////////
		//DEFINE CUT//
		//////////////
		//float ptot = sqrt(_pT*_pT + _pz*_pz);
		//cut_pdtheta__ = cut_pdtheta[_arm][_lastgap] + 0.01*_pT;
		cut_pdtheta__ = cut_pdtheta[_arm][_lastgap];
		cut_mutr_chi2__ = cut_mutr_chi2[_arm][_lastgap];
		cut_muid_chi2__ = 6.0;

		cut_dg0__ = hcut_dg0[_arm][_lastgap]->GetBinContent(hcut_dg0[_arm][_lastgap]->FindBin(_pT));
		cut_ddg0__ = hcut_ddg0[_arm][_lastgap]->GetBinContent(hcut_ddg0[_arm][_lastgap]->FindBin(_pT));
		cut_vertex_rad__ = hcut_vertex_rad[_arm][_lastgap]->GetBinContent(hcut_vertex_rad[_arm][_lastgap]->FindBin(_pT));
		cut_vertex_chi2__ = hcut_vertex_chi2[_arm][_lastgap]->GetBinContent(hcut_vertex_chi2[_arm][_lastgap]->FindBin(_pT));
		/*
		if ( _pT<2.0 ){
			cut_vertex_chi2__ *= 1.5;
		}
		*/

		if ( cuts=="LOOSE" ){
			cut_dg0__ *= 1.25;
			cut_ddg0__ *= 1.25;
			cut_vertex_rad__ *= 1.25;
			cut_vertex_chi2__ *= 1.25;
		}

		//if ( (cut_pdtheta__-cut_pdtheta[_arm][_lastgap])>0.2 ) cut_pdtheta__ = cut_pdtheta[_arm][_lastgap] + 0.2;

		//Fill cut variable distributions
		FillDiagnosticHistos(weight_by_thrown_pt);

		////////////////////////////////////////////////////////////
		//Perform analysis cuts and fill final analysis histograms//
		////////////////////////////////////////////////////////////
		if ( _lastgap==2 || _lastgap==3 ){
			FillGap23(weight_by_thrown_pt, COCKTAIL);
		}//lastgap23

		////////////////
		//FILL LASTGAP//
		////////////////
		if ( _lastgap==4 ){
			if ( 
					_rf_slope>cut_road_slope__ //1 
					&& _rf_refrad<cut_vertex_rad__ //2
					&& _DG0<cut_dg0__ //3
					&& _DDG0<cut_ddg0__ //4
					&& _trchi2<cut_mutr_chi2__ //5
					&& _idchi2<cut_muid_chi2__ //6
					&& _pdtheta<cut_pdtheta__ //7
					&& _rf_vtxchi2pdf<cut_vertex_chi2__ //8
				 ){ 

				//both pdtheta and deltaz cut for inclusive yield
				avg_pT_varbin[_arm][_charge]->Fill(_pT, _pT, weight_by_thrown_pt);
				avg_pT_varbin_gap[_arm][_lastgap][_charge]->Fill(_pT, _pT, weight_by_thrown_pt);
				avg_p[_arm][_lastgap][_charge]->Fill(_pT, sqrt(_pT*_pT+_pz*_pz), weight_by_thrown_pt);
				avg_pT_AN_pz[_arm][_charge]->Fill(fabs(_pz), _pT, weight_by_thrown_pt);

				inclZ[_arm][_charge]->Fill(_pT, _bbcZ, weight_by_thrown_pt);
				N_varbin[_arm][_lastgap][_charge]->Fill(_pT, weight_by_thrown_pt);

				N_pT_pz[_arm][_lastgap][_charge]->Fill(_pT, fabs(_pz), weight_by_thrown_pt);

				if ( _pT>1.25 && fabs(_pz)>3.50 ){
					N_varbin_AN_pT[_arm][_lastgap][_charge]->Fill(_pT, weight_by_thrown_pt);
					N_varbin_AN_pz[_arm][_lastgap][_charge]->Fill(fabs(_pz), weight_by_thrown_pt);

					if ( _pT<2.50 ){
						N_varbin_AN_pz_MUID1D[_arm][_lastgap][_charge]->Fill(fabs(_pz), weight_by_thrown_pt);
					}else{
						N_varbin_AN_pz_SG3MUID1DH[_arm][_lastgap][_charge]->Fill(fabs(_pz), weight_by_thrown_pt);
					}
				}

				if ( filetype!="DATA" && _mc_g_pT<0.8*_pT ){
					N_varbin_fake[_arm][_lastgap][_charge]->Fill(_pT, weight_by_thrown_pt);
				}

				if ( COCKTAIL ){
					//if ( (_mc_pid==5 || _mc_pid==6) && (_pT-_mc_g_pT)<(0.5*_pT) ){
					//if ( 1 ){
					if ( abs(_mc_pid)==13 && fabs(_mc_z)<80.0 && (_pT-_mc_g_pT)<(0.5*_pT) ){
						tuning_matrix_allcuts[_arm][_lastgap][_charge]->Fill(_mc_g_pT, _pT, weight_by_thrown_pt);
					}

					if ( abs(_mc_pid)==13 ){
						if ( fabs(_mc_z)<80.0 ){
							Ndecay[_arm][_lastgap][_charge]->Fill(_pT, weight_by_thrown_pt);
							inclZ_decay[_arm][_charge]->Fill(_pT, _bbcZ, weight_by_thrown_pt);
							if ( _pT>1.25 && fabs(_pz)>3.50 ){
								Ndecay_AN_pT[_arm][_lastgap][_charge]->Fill(_pT, weight_by_thrown_pt);
								Ndecay_AN_pz[_arm][_lastgap][_charge]->Fill(fabs(_pz), weight_by_thrown_pt);
							}
						}else{
							Ndecayabs[_arm][_lastgap][_charge]->Fill(_pT, weight_by_thrown_pt);
							inclZ_decay[_arm][_charge]->Fill(_pT, _bbcZ, weight_by_thrown_pt);
							if ( _pT>1.25 && fabs(_pz)>3.50 ){
								Ndecayabs_AN_pT[_arm][_lastgap][_charge]->Fill(_pT, weight_by_thrown_pt);
								Ndecayabs_AN_pz[_arm][_lastgap][_charge]->Fill(fabs(_pz), weight_by_thrown_pt);
							}
						}
					}else{
						Npunchthrough[_arm][_lastgap][_charge]->Fill(_pT, weight_by_thrown_pt);
						inclZ_punchthrough[_arm][_charge]->Fill(_pT, _bbcZ, weight_by_thrown_pt);
						if ( _pT>1.25 && fabs(_pz)>3.50 ){
							Npunchthrough_AN_pT[_arm][_lastgap][_charge]->Fill(_pT, weight_by_thrown_pt);
							Npunchthrough_AN_pz[_arm][_lastgap][_charge]->Fill(fabs(_pz), weight_by_thrown_pt);
						}
					}

				}//COCKTAIL

				if ( 1 ){
					YvX_mutr_cut[_arm][0]->Fill(_xst1, _yst1, weight_by_thrown_pt);
					Rphi_mutr_cut[_arm][0]->Fill(_phist1, _radst1, weight_by_thrown_pt);
					YvX_mutr_cut[_arm][1]->Fill(_xst2, _yst2, weight_by_thrown_pt);
					Rphi_mutr_cut[_arm][1]->Fill(_phist2, _radst2, weight_by_thrown_pt);
					YvX_mutr_cut[_arm][2]->Fill(_xst3, _yst3, weight_by_thrown_pt);
					Rphi_mutr_cut[_arm][2]->Fill(_phist3, _radst3, weight_by_thrown_pt);
				}

			}//quality cut
		}//lastgap4
	}//jentry

	cout << "Writing output file : " << fname << endl;
	outfile->Write();
	outfile->Close();
}

void run_cuts::BookHistos(){

	bool COCKTAIL = kFALSE;
	bool MUON = kFALSE;

	if ( filetype=="QGSP_BERT" || filetype=="FTFP_BERT" || filetype=="QGSP_BIC" || filetype=="COCKTAIL" ) COCKTAIL = kTRUE;
	if ( filetype=="MUON" ) MUON = kTRUE;

	char hname[200];

	//int nbin250 = 28;
	if ( COCKTAIL || MUON ){
		for (int iarm=0; iarm<narm; iarm++){
			for (int igap=2; igap<ngap; igap++){
				for (int ich=0; ich<ncharge; ich++){
					sprintf(hname,"ndecay_arm%d_gap%d_chg%d",iarm,igap,ich);
					Ndecay[iarm][igap][ich] = new TH1F(hname,"",nptbin,varbin_pTarray);
					Ndecay[iarm][igap][ich]->Sumw2();

					sprintf(hname,"ndecayabs_arm%d_gap%d_chg%d",iarm,igap,ich);
					Ndecayabs[iarm][igap][ich] = new TH1F(hname,"",nptbin,varbin_pTarray);
					Ndecayabs[iarm][igap][ich]->Sumw2();

					sprintf(hname,"npunchthrough_arm%d_gap%d_chg%d",iarm,igap,ich);
					Npunchthrough[iarm][igap][ich] = new TH1F(hname,"",nptbin,varbin_pTarray);
					Npunchthrough[iarm][igap][ich]->Sumw2();

					sprintf(hname,"ndecay_AN_pT_arm%d_gap%d_chg%d",iarm,igap,ich);
					Ndecay_AN_pT[iarm][igap][ich] = new TH1F(hname,"",9,varbin_AN_pTarray);
					Ndecay_AN_pT[iarm][igap][ich]->Sumw2();

					sprintf(hname,"ndecayabs_AN_pT_arm%d_gap%d_chg%d",iarm,igap,ich);
					Ndecayabs_AN_pT[iarm][igap][ich] = new TH1F(hname,"",9,varbin_AN_pTarray);
					Ndecayabs_AN_pT[iarm][igap][ich]->Sumw2();

					sprintf(hname,"npunchthrough_AN_pT_arm%d_gap%d_chg%d",iarm,igap,ich);
					Npunchthrough_AN_pT[iarm][igap][ich] = new TH1F(hname,"",9,varbin_AN_pTarray);
					Npunchthrough_AN_pT[iarm][igap][ich]->Sumw2();

					sprintf(hname,"ndecay_AN_pz_arm%d_gap%d_chg%d",iarm,igap,ich);
					Ndecay_AN_pz[iarm][igap][ich] = new TH1F(hname,"",4,varbin_AN_pzarray);
					Ndecay_AN_pz[iarm][igap][ich]->Sumw2();

					sprintf(hname,"ndecayabs_AN_pz_arm%d_gap%d_chg%d",iarm,igap,ich);
					Ndecayabs_AN_pz[iarm][igap][ich] = new TH1F(hname,"",4,varbin_AN_pzarray);
					Ndecayabs_AN_pz[iarm][igap][ich]->Sumw2();

					sprintf(hname,"npunchthrough_AN_pz_arm%d_gap%d_chg%d",iarm,igap,ich);
					Npunchthrough_AN_pz[iarm][igap][ich] = new TH1F(hname,"",4,varbin_AN_pzarray);
					Npunchthrough_AN_pz[iarm][igap][ich]->Sumw2();

					//TUNNING MATRIX
					sprintf(hname, "tuning_matrix_allcuts_arm%d_gap%d_chg%d",iarm,igap,ich);
					tuning_matrix_allcuts[iarm][igap][ich] = new TH2F(hname,"",nmcptbin,varbin_mc_pTarray,nptbin,varbin_pTarray);
				}//CHARGE
			}//GAP
		}//ARM
	}//COCKTAIL||MUON

	float zrange = cut_z[1];
	int nzbin = int(2*cut_z[1] + 0.1);

	for (int iarm=0; iarm<narm; iarm++){
		for (int ich=0; ich<ncharge; ich++){
			// variable bin pT profile
			sprintf(hname, "avg_pT_varbin_arm%d_chg%d", iarm, ich);
			avg_pT_varbin[iarm][ich]= new TProfile(hname, "", nptbin, varbin_pTarray);

			sprintf(hname, "avg_pT_AN_pz_arm%d_chg%d", iarm, ich);
			avg_pT_AN_pz[iarm][ich]= new TProfile(hname, "", 4, varbin_AN_pzarray);

			for (int igap=2; igap<ngap; igap++){
				// variable bin pT profile
				sprintf(hname, "avg_pT_varbin_arm%d_gap%d_chg%d", iarm, igap, ich);
				avg_pT_varbin_gap[iarm][igap][ich]= new TProfile(hname, "", nptbin, varbin_pTarray);

				sprintf(hname, "avg_p_varbin_arm%d_gap%d_chg%d", iarm, igap, ich);
				avg_p[iarm][igap][ich]= new TProfile(hname, "", nptbin, varbin_pTarray);

				//gap4 variable bin spectra
				sprintf(hname, "n_varbin_arm%d_gap%d_chg%d", iarm, igap, ich);
				N_varbin[iarm][igap][ich] = new TH1F(hname, "", nptbin, varbin_pTarray);
				N_varbin[iarm][igap][ich]->Sumw2();

				sprintf(hname, "n_varbin_AN_pT_arm%d_gap%d_chg%d", iarm, igap, ich);
				N_varbin_AN_pT[iarm][igap][ich] = new TH1F(hname, "", 9, varbin_AN_pTarray);
				N_varbin_AN_pT[iarm][igap][ich]->Sumw2();

				sprintf(hname, "n_varbin_AN_pz_arm%d_gap%d_chg%d", iarm, igap, ich);
				N_varbin_AN_pz[iarm][igap][ich] = new TH1F(hname, "", 4, varbin_AN_pzarray);
				N_varbin_AN_pz[iarm][igap][ich]->Sumw2();

				sprintf(hname, "n_varbin_AN_pz_MUID1D_arm%d_gap%d_chg%d", iarm, igap, ich);
				N_varbin_AN_pz_MUID1D[iarm][igap][ich] = new TH1F(hname, "", 4, varbin_AN_pzarray);
				N_varbin_AN_pz_MUID1D[iarm][igap][ich]->Sumw2();

				sprintf(hname, "n_varbin_AN_pz_SG3MUID1DH_arm%d_gap%d_chg%d", iarm, igap, ich);
				N_varbin_AN_pz_SG3MUID1DH[iarm][igap][ich] = new TH1F(hname, "", 4, varbin_AN_pzarray);
				N_varbin_AN_pz_SG3MUID1DH[iarm][igap][ich]->Sumw2();

				sprintf(hname, "n_varbin_fake_arm%d_gap%d_chg%d", iarm, igap, ich);
				N_varbin_fake[iarm][igap][ich] = new TH1F(hname, "", nptbin, varbin_pTarray);
				N_varbin_fake[iarm][igap][ich]->Sumw2();

				sprintf(hname, "n_pT_pz_arm%d_gap%d_chg%d", iarm, igap, ich);
				N_pT_pz[iarm][igap][ich] = new TH2F(hname, "", 30, 0, 15, 30, 0, 45);
				N_pT_pz[iarm][igap][ich]->Sumw2();

				if ( COCKTAIL ){
					sprintf(hname, "n_varbin_pion_arm%d_gap%d_chg%d", iarm, igap, ich);
					N_varbin_pion[iarm][igap][ich] = new TH1F(hname, "", nptbin, varbin_pTarray);
					N_varbin_pion[iarm][igap][ich]->Sumw2();

					sprintf(hname, "n_varbin_kaon_arm%d_gap%d_chg%d", iarm, igap, ich);
					N_varbin_kaon[iarm][igap][ich] = new TH1F(hname, "", nptbin, varbin_pTarray);
					N_varbin_kaon[iarm][igap][ich]->Sumw2();

					sprintf(hname, "n_varbin_kaon0_arm%d_gap%d_chg%d", iarm, igap, ich);
					N_varbin_kaon0[iarm][igap][ich] = new TH1F(hname, "", nptbin, varbin_pTarray);
					N_varbin_kaon0[iarm][igap][ich]->Sumw2();

					sprintf(hname, "n_varbin_proton_arm%d_gap%d_chg%d", iarm, igap, ich);
					N_varbin_proton[iarm][igap][ich] = new TH1F(hname, "", nptbin, varbin_pTarray);
					N_varbin_proton[iarm][igap][ich]->Sumw2();
				}

			}//GAP

			//! inclZ 
			sprintf(hname, "inclZ_arm%d_chg%d", iarm, ich);
			inclZ[iarm][ich] = new TH2F(hname, "", nptbin, varbin_pTarray, nzbin, -zrange, zrange);
			inclZ[iarm][ich]->Sumw2();


			if ( COCKTAIL ){
				sprintf(hname, "inclZ_decay_arm%d_chg%d", iarm, ich);
				inclZ_decay[iarm][ich] = new TH2F(hname, "", nptbin, varbin_pTarray, nzbin, -zrange, zrange);
				inclZ_decay[iarm][ich]->Sumw2();

				sprintf(hname, "inclZ_punchthrough_arm%d_chg%d", iarm, ich);
				inclZ_punchthrough[iarm][ich] = new TH2F(hname, "", nptbin, varbin_pTarray, nzbin, -zrange, zrange);
				inclZ_punchthrough[iarm][ich]->Sumw2();
			}
		}//CHARGE

		for (int ista=0; ista<nsta; ista++){
			int bound = (ista + 1)*160;

			if ( ista==0 ) bound = 140;
			else if ( ista==1 ) bound = 240;
			else if ( ista==2 ) bound = 440;

			sprintf(hname, "YvX_arm%d_sta%d_cut", iarm, ista+1);
			YvX_mutr_cut[iarm][ista] = new TH2F(hname, "", bound, -bound, bound, bound, -bound, bound);

			sprintf(hname, "Rphi_arm%d_sta%d_cut", iarm, ista+1);
			//Rphi_mutr_cut[iarm][ista] = new TH2F(hname, "", 64, -TMath::Pi(), TMath::Pi(), bound/2, 0, bound);
			Rphi_mutr_cut[iarm][ista] = new TH2F(hname, "", 64, 0, 2*TMath::Pi(), bound/2, 0, bound);
		}

		//KINEMATIC CUTS FOR GAP 2, 3, 4
		for (int jgap=2; jgap<ngap; jgap++){
			//DG0
			sprintf(hname, "dg0_arm%d_gap%d", iarm, jgap);
			DG0[iarm][jgap] = new TH2F(hname,"",nptbin,varbin_pTarray,80,0,40);
			DG0[iarm][jgap]->Sumw2();
			//DDG0
			sprintf(hname, "ddg0_arm%d_gap%d", iarm, jgap);
			DDG0[iarm][jgap] = new TH2F(hname,"",nptbin,varbin_pTarray,60,0,30);
			DDG0[iarm][jgap]->Sumw2();
			//REFRAD
			sprintf(hname, "refrad_arm%d_gap%d", iarm, jgap);
			REFRAD[iarm][jgap] = new TH2F(hname,"",nptbin,varbin_pTarray,250,0.,250.0);
			REFRAD[iarm][jgap]->Sumw2();
			//MUTR CHI2
			sprintf(hname, "mutr_chi2_arm%d_gap%d", iarm, jgap);
			MUTR_CHI2[iarm][jgap] = new TH2F(hname,"",nptbin,varbin_pTarray,100,0,25);
			MUTR_CHI2[iarm][jgap]->Sumw2();
			//VTX CHI2
			sprintf(hname, "vtx_chi2_arm%d_gap%d",iarm, jgap);
			VTX_CHI2[iarm][jgap] = new TH2F(hname,"",nptbin,varbin_pTarray,400,0,20);
			VTX_CHI2[iarm][jgap]->Sumw2();
			if ( filetype!="DATA" ){
				sprintf(hname, "vtx_chi2_out_arm%d_gap%d",iarm, jgap);
				VTX_CHI2_OUT[iarm][jgap] = new TH2F(hname,"",nptbin,varbin_pTarray,400,0,20);
				VTX_CHI2_OUT[iarm][jgap]->Sumw2();
			}
			//PDTHETA
			sprintf(hname, "pdtheta_arm%d_gap%d", iarm, jgap);
			PDTHETA[iarm][jgap] = new TH2F(hname,hname,nptbin,varbin_pTarray,50,0.0,1.0);
			PDTHETA[iarm][jgap]->Sumw2();
			//PDTHETA_FAKE
			sprintf(hname, "pdtheta_fake_arm%d_gap%d", iarm, jgap);
			PDTHETA_FAKE[iarm][jgap] = new TH2F(hname,hname,nptbin,varbin_pTarray,50,0.0,1.0);
			PDTHETA_FAKE[iarm][jgap]->Sumw2();
			//PDTHETA_OUT
			sprintf(hname, "pdtheta_out_arm%d_gap%d", iarm, jgap);
			PDTHETA_OUT[iarm][jgap] = new TH2F(hname,hname,nptbin,varbin_pTarray,50,0.0,1.0);
			PDTHETA_OUT[iarm][jgap]->Sumw2();
		}//GAP-2
	}//ARM


}

void run_cuts::InitializeCuts(){

	if ( (cuts=="TIGHT") || (cuts=="LOOSE") || (cuts=="NOCUT") ){ 
		cout << "***inside InitializeCuts*** cuts variable is set to - " << cuts << endl;
	}else{
		cout << "Need to add t->SetKinematicsCuts method" << endl;
		exit(1);
	}

	//double par_dg0[narm][ngap][2];
	//double par_ddg0[narm][ngap][2];
	//double par_vertex_rad[narm][ngap][2];
	//double par_vertex_chi2[narm][ngap][2];

	cut_road_slope__ = 0.15;

	if ( cuts=="LOOSE" ){
		cut_pdtheta[0][2] = 0.35;
		cut_pdtheta[0][3] = 0.35;
		cut_pdtheta[0][4] = 0.25;
		cut_pdtheta[1][2] = 0.35;
		cut_pdtheta[1][3] = 0.35;
		cut_pdtheta[1][4] = 0.25;

		cut_mutr_chi2[0][2] = 20;
		cut_mutr_chi2[0][3] = 20;
		cut_mutr_chi2[0][4] = 20;
		cut_mutr_chi2[1][2] = 20;
		cut_mutr_chi2[1][2] = 20;
		cut_mutr_chi2[1][3] = 20;
		cut_mutr_chi2[1][4] = 20;

		/*
		cut_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/15.inclusiveHF/run_cuts/cut_files/Run15pp200_HFmu_tight_cut.root","READ");
		for (int iarm=0; iarm<narm; iarm++){
			for (int igap=2; igap<ngap; igap++){
				fcut_dg0[iarm][igap] = (TF1*)cut_file->Get(Form("fcut_dg0_arm%d_gap%d",iarm,igap));
				fcut_ddg0[iarm][igap] = (TF1*)cut_file->Get(Form("fcut_ddg0_arm%d_gap%d",iarm,igap));
				fcut_vertex_rad[iarm][igap] = (TF1*)cut_file->Get(Form("fcut_refrad_arm%d_gap%d",iarm,igap));
				fcut_vertex_chi2[iarm][igap] = (TF1*)cut_file->Get(Form("fcut_vtx_chi2_arm%d_gap%d",iarm,igap));
			}
		}
		*/
		cut_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/15.inclusiveHF/run_cuts/cut_files/Run15pp200_HFmu_tight_cut_20170404.root","READ");
		for (int iarm=0; iarm<narm; iarm++){
			for (int igap=2; igap<ngap; igap++){
				hcut_dg0[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_dg0_arm%d_gap%d",iarm,igap));
				hcut_ddg0[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_ddg0_arm%d_gap%d",iarm,igap));
				hcut_vertex_rad[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_refrad_arm%d_gap%d",iarm,igap));
				hcut_vertex_chi2[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_vtx_chi2_arm%d_gap%d",iarm,igap));
			}
		}

	}else if ( cuts=="TIGHT" ){
		cut_pdtheta[0][2] = 0.40;
		cut_pdtheta[0][3] = 0.40;
		cut_pdtheta[0][4] = 0.25;
		cut_pdtheta[1][2] = 0.40;
		cut_pdtheta[1][3] = 0.40;
		cut_pdtheta[1][4] = 0.25;

		cut_mutr_chi2[0][2] = 15;
		cut_mutr_chi2[0][3] = 15;
		cut_mutr_chi2[0][4] = 15;
		cut_mutr_chi2[1][2] = 20;
		cut_mutr_chi2[1][3] = 20;
		cut_mutr_chi2[1][4] = 20;

		cut_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/15.inclusiveHF/run_cuts/cut_files/Run15pp200_HFmu_tight_cut_20170404.root","READ");
		cout << "OPEN CUT FILE: " << cut_file->GetName() << endl;
		for (int iarm=0; iarm<narm; iarm++){
			for (int igap=2; igap<ngap; igap++){
				hcut_dg0[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_dg0_arm%d_gap%d",iarm,igap));
				hcut_ddg0[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_ddg0_arm%d_gap%d",iarm,igap));
				hcut_vertex_rad[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_refrad_arm%d_gap%d",iarm,igap));
				hcut_vertex_chi2[iarm][igap] = (TH1F*)cut_file->Get(Form("hcut_vtx_chi2_arm%d_gap%d",iarm,igap));
			}
		}
	}//TIGHT|LOOSE


	cout << "Initialized cuts for dataset : " << dataset << "  Cut type : " << cuts << endl;

	ftrig = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/10.runQA/01.trig_eff/Run15pp200_trig_eff_func_20170301.root","READ");
	cout << "OPEN trigger efficiency file: " << ftrig->GetName() << endl;
	for (int iarm=0; iarm<narm; iarm++){
		for (int igap=2; igap<ngap; igap++){

			ftrig_SG3_pT[iarm][igap] = (TF1*)ftrig->Get(Form("ftrig_eff_pT_SG3_MUID1DH_arm%d_gap%d",iarm,igap));

			if ( igap==4 ){
				ftrig_MUID1D_pT[iarm] = (TF1*)ftrig->Get(Form("ftrig_eff_pT_MUID1D_arm%d_gap4",iarm));
			}

		}//igap
	}//iarm

	if ( !ftrig_SG3_pT[0][4] || !ftrig_SG3_pT[1][4] || !ftrig_MUID1D_pT[0] || !ftrig_MUID1D_pT[1] ){
		cout << "CAN NOT!! Initialized triggr efficiency file!!" << endl;
		exit(1);
	}else{
		cout << "Initialized triggr efficiency file : " << dataset << "  Trigger type : " << triggertype << endl;
	}

	fSD = new TFile(Form("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/07.genericQA/03.prescale/%s/%s_prescale.root",dataset.c_str(),dataset.c_str()),"read");
	hSD_MUID1D[0] = (TH1F*)fSD->Get("hSD_MUID1D_S");
	hSD_MUID1D[1] = (TH1F*)fSD->Get("hSD_MUID1D_N");
	hSD_SG3_MUID1DH[0] = (TH1F*)fSD->Get("hSD_SG3_MUID1DH_S");
	hSD_SG3_MUID1DH[1] = (TH1F*)fSD->Get("hSD_SG3_MUID1DH_N");

	fchg_wt = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/04.run9_pdst/08.input_study/pp200_kpi_ratio_chg_wt.root","read");
	for (int ichg=0; ichg<ncharge; ichg++){
		for (int ieta=0; ieta<8; ieta++){
			hW_KAON[ichg][ieta] = (TH1F*)fchg_wt->Get(Form("hW_KAON_CHG%d_ETA%d",ichg,11+2*ieta));
			hW_KAON[ichg][ieta]->Print();
		}
	}
}

void run_cuts::Load_tuning_weights(){

	float temp_weight;
	int temp_arm, temp_pt, temp_chg;

	ifstream weightfile;

	cout << "FILETYPE : " << filetype << endl;
	if ( (filetype=="DATA") || (filetype=="MUON") ){
		cout << "Untuned specified, all tuning weights set to 1." << endl;
		for (int iarm=0; iarm<narm; iarm++){
			for (int ich=0; ich<ncharge; ich++){
				for (int ipt=0; ipt<nmcptbin; ipt++){
					tune_weight[iarm][ich][ipt] = 1.0;
				}
			}
		}
	}else {

		char weight_file[200];
		sprintf(weight_file, "weight/Weight_%s.txt",filetype.c_str());

		cout << weight_file << endl;

		weightfile.open(weight_file);

		if ( !weightfile.good() ){
			cout << "Invalid weight file..." << endl;
			exit(1);
		}else{
			cout << "Using the following weight files : " << weight_file << endl;
		}

		while ( weightfile >> temp_arm >> temp_chg >> temp_pt >> temp_weight ){
			tune_weight[temp_arm][temp_chg][temp_pt] = temp_weight;
		}

		weightfile.close();
	}//filetype

	cout << "Applying pT weight for case : " << dataset << "\t" << filetype << endl;

}

bool run_cuts::CheckRadCut(){
	/////////////////
	//MUTR RAD CUTS//
	/////////////////
	bool radgood = kTRUE;

	/*
	if ( _arm==1 && _halfoctst2==9 && (_radst2>140 && _radst2<160) ) radgood = false;
	if ( _arm==1 && _halfoctst2==10 && (_radst2>140 && _radst2<160) ) radgood = false;
	*/

	return radgood;
} 

bool run_cuts::CheckPhiCut(){
	/////////////////
	//MUTR PHI CUTS//
	/////////////////

	//bool phigood = true;

	if ( _arm==1 ){
		if ( _phist1>0.40 && _phist1<1.20 ) return false;
		if ( _phist2>0.40 && _phist2<1.20 ) return false;
		if ( _phist3>0.40 && _phist3<1.20 ) return false;
		if ( _phist1>-1.60 && _phist1<-0.88 ) return false;
		if ( _phist2>-1.60 && _phist2<-0.88 ) return false;
		if ( _phist3>-1.60 && _phist3<-0.88 ) return false;
		if ( _phist1>3.02 ) return false;
		if ( _phist2>3.02 ) return false;
		if ( _phist3>3.02 ) return false;
		if ( _phist1<-2.75 ) return false;
		if ( _phist2<-2.75 ) return false;
		if ( _phist3<-2.75 ) return false;
	}else{
		if ( _phist1>0 && _phist1<0.4 ) return false;
		if ( _phist2>0 && _phist2<0.4 ) return false;
		if ( _phist3>0 && _phist3<0.4 ) return false;
		if ( _phist1>1.55 && _phist1<1.9 ) return false;
		if ( _phist2>1.55 && _phist2<1.9 ) return false;
		if ( _phist3>1.55 && _phist3<1.9 ) return false;
		if ( _phist1>-0.60 && _phist1<-0.40 ) return false;
		if ( _phist2>-0.60 && _phist2<-0.40 ) return false;
		if ( _phist3>-0.60 && _phist3<-0.40 ) return false;
	}

	return true;
}

