/*
Author : Sanghwa Park 
User : Muhammad Alibordi  

NOTE:
  This simple example macro shows how to read contents 
  from spinDB and store them into a rootfile
  
  Takes a run list as input

  It uses uspin library 
  For details check source codes from:
  offline/packages/uspin

  To quickly check list of functions to use
  one can also see from:
  $OFFLINE_MAIN/include/SpinDBContent.hh

 */

const int NB = 120; 

// Output tree variables
int default_qa;
int runnumber;
int fillnumber;
int xingshift;
int badrun_flag;

int patternblue[NB];
int patternyell[NB];
int badbunch_qa[NB];

// GL1p scalers
long long scaler_bbcvtxcut[NB];
long long scaler_bbcnovtx[NB];
long long scaler_zdcwide[NB];
long long scaler_zdcnarrow[NB];

// polarization blue, yellow, stat and sys error
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

void SetOutTree(TTree* tree);
void InitTreeVars();
int FindPattern(vector<int> &ptb, vector<int> &pty);
int FindGroup(int pattern);
void DefineBasePatterns();

int verb = 0;

void SpinExample2(string runlist)
{

  DefineBasePatterns();

  gSystem->Load("libuspin.so");

  TFile* fout = new TFile("spinDB_test.root", "RECREATE");
  TTree* T = new TTree("T", "T");
  SetOutTree(T);

  vector<int> v_ptb;
  vector<int> v_pty;

  SpinDBContent spin_cont;
  SpinDBOutput spin_out("phnxrc");

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

      /*
	Get spin pattern for yellow and blue beams

	GL1P scaler names in the spin DB
	:bbcvertexcut  -> ScalerA
	:bbcwithoutcut -> ScalerB
	:zdcnarrow     -> ScalerC
	:zdcwide       -> ScalerD
	Note that they are **NOT** necessarily what the name say
	To find inputs for each year, visit the wiki page:
	https://www.phenix.bnl.gov/WWW/offline/wikioff/index.php/GL1P_Information
      */    
      
      for(int ibunch=0; ibunch<NB; ibunch++)
	{
	  patternblue[ibunch] = spin_cont.GetSpinPatternBlue(ibunch);
	  patternyell[ibunch] = spin_cont.GetSpinPatternYellow(ibunch);
	  badbunch_qa[ibunch] = spin_cont.GetBadBunchFlag(ibunch);
	  
	  scaler_bbcvtxcut[ibunch] = spin_cont.GetScalerBbcVertexCut(ibunch);
	  scaler_bbcnovtx[ibunch]  = spin_cont.GetScalerBbcNoCut(ibunch);
	  scaler_zdcwide[ibunch]   = spin_cont.GetScalerZdcWide(ibunch);
	  scaler_zdcnarrow[ibunch] = spin_cont.GetScalerZdcNarrow(ibunch);

	  //fill first 16 bunches
	  if(ibunch < 16)
	    {
	      v_ptb.push_back(patternblue[ibunch]);
	      v_pty.push_back(patternyell[ibunch]);
	    }
	}

      pattern_number = FindPattern(v_ptb, v_pty);
      group = FindGroup(pattern_number);

      T->Fill();

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

  fout->cd();
  T->Write();
  fout->Close();

  return;
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

void SetOutTree(TTree* tree)
{
  tree->Branch("runnumber",      &runnumber,         "runnumber/I");
  tree->Branch("fillnumber",     &fillnumber,        "fillnumber/I");
  tree->Branch("qa_level",       &default_qa,        "qa_level/I");
  tree->Branch("xingshift",      &xingshift,         "xingshift/I");
  tree->Branch("badrun_flag",    &badrun_flag,       "badrun_flag/I");
  tree->Branch("polblue",        &b_pol,             "polblue/F");
  tree->Branch("polblue_stat",   &b_stat,            "polblue_stat/F");
  tree->Branch("polblue_sys",    &b_syst,            "polblue_sys/F");
  tree->Branch("polyellow",      &y_pol,             "polyellow/F");
  tree->Branch("polyellow_stat", &y_stat,            "polyellow_stat/F");
  tree->Branch("polyellow_sys",  &y_syst,            "polyellow_sys/F");
  tree->Branch("patternblue",    patternblue,        "patternblue[120]/I");
  tree->Branch("patternyellow",  patternyell,        "patternyellow[120]/I");
  tree->Branch("badbunch_qa",    badbunch_qa,        "badbunch_qa[120]/I");
  tree->Branch("scalerA",        scaler_bbcvtxcut,   "scalerA[120]/L");
  tree->Branch("scalerB",        scaler_bbcnovtx,    "scalerB[120]/L");
  tree->Branch("scalerC",        scaler_zdcwide,     "scalerC[120]/L");
  tree->Branch("scalerD",        scaler_zdcnarrow,   "scalerD[120]/L");
  tree->Branch("pattern",        &pattern_number,    "pattern/I");
  tree->Branch("group",          &group,             "group/I");

  return;

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

  // intended base patterns (gives 8 patterns from 1/2 + 3a/4a combinations)
  /*
    1:  - - + + - - + + 
    2:  + + - - + + - - 
    3a: - - + + + + - - 
    4a: + + - - - - + +

    yuk... root5 doesn't support c++11
  // Define base vectors for first 16 crossings
  // + + - - + + - - + + - - 
  std::vector<int> base1 = {1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1};
  // - - + + - - + + - - + + 
  std::vector<int> base2 = {-1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1};
  // + + + + - - - - + + + + - - 
  std::vector<int> base3 = {1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1};
  // - - - - + + + + - - - - + + 
  std::vector<int> base4 = {-1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1};

  //P1 (B:base1 Y:base3)
  //P2 (B:base2 Y:base3)
  //P3 (B:base1 Y:base4)
  //P4 (B:base2 Y:base4)
  //P5 (B:base3 Y:base1)
  //P6 (B:base3 Y:base2)
  //P7 (B:base4 Y:base1)
  //P8 (B:base4 Y:base2)
  
  // + + - - + + - - 
  std::vector<int> base1a = {1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1};
  // - - + + - - + + 
  std::vector<int> base2a = {-1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1};
  // + + + + - - - -
  std::vector<int> base3a = {-1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1};
  // - - - - + + + +
  std::vector<int> base4a = {1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1};

  //Note: patterns defined based on the table in AN1125, CAD defines P22<->P23
  //P21 (B:base1a Y:base3a)
  //P22 (B:base2a Y:base3a)
  //P23 (B:base1a Y:base4a)
  //P24 (B:base2a Y:base4a)
  //P25 (B:base3a Y:base1a)
  //P26 (B:base3a Y:base2a)
  //P27 (B:base4a Y:base1a)
  //P28 (B:base4a Y:base2a)
  */

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
    return -1; // can't find a matching pattern
  
}

int FindGroup(int pattern)
{

  /* SOOSSOO P1 P4 P5 P8
     OSSOOSS P2 P3 P6 P7
     SSOO    P21 P24 P25 P28
     OOSS    P22 P23 P26 P27
  */
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
