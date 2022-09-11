//We have to include the headers , but any for the moment this is working


void timediff_treemaker(unsigned fnum ){

	auto t1450 = new LAPPD4EIC();
	t1450->ImportDrs4Calibrations("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/alibordi/Test_Beam/mpgd4eic/v1742.db", 0, 5);
	t1450->AddV1742module(2071, 12063);
	t1450->AddV1742module(2000, 12064);
  	t1450->AddV1742module(2010,   106);
  	t1450->AddV1742module(2020, 10906);
  	t1450->AddV1742module(2030, 12067);
  	t1450->AddV1742module(2040,    81);
  	t1450->AddV1742module(2050,    97);
  	t1450->AddV1742module(2060,   120);
 	t1450->AddV1742module(2070,    87);
  	t1450->AddV1742module(2071, 12063);
  	t1450->AddV1742module(2072, 12065);
	t1450->SetStat(10000);
	t1450->SetRcdaqFileNameMask("/gpfs/mnt/gpfs02/eic/TEST.RUNS/2022-FNAL/calibration/calibration_lappd-%08d-0000.evt");
	t1450->SetOutputFileNameMask("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/alibordi/Test_Beam/output_lappd/timediff_tree_%08d.root");
	auto v1742 = t1450->m_V1742s[0];
	t1450->SetPlanaconChannels(v1742->m_Channels[2], v1742->m_Channels[5]);
	t1450->UseCaenAndNaiveOffsetChannelCalibration(2100, 3.72);
        t1450->UseCaenAndNaiveOffsetTriggerCalibration(2200, 1.77);	
	t1450->SetNIMPulseFitParameters(-125.0, -450.0, -300.0);
	t1450->SetNIMPulseBaselineWindow(-70, 30);
	t1450->SetPlanaconPulseSearchThreshold(-25.0);
	t1450->SetPlanaconPulseBaselineWindow(-25, 10);
	t1450->SetPlanaconPulsePeakWindow(20);
	t1450->SetPlanaconPulseFitParameters(0.25, 0.70, 0.50);
	t1450->CorrectBaselineInRawWaveForms(1, 20);
	t1450->SetRegularChannelPeakWindow(-115, 25);
	t1450->SetRegularChannelBaselineWindow(-50, 20);
	t1450->ImportMappingTable("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/alibordi/Test_Beam/mpgd4eic/gerber/output/L03c/L03c.map.v01.root", "TR");
	auto mtx = t1450->m_PixelMatrices[0];
        mtx->SetPitch(4.0);
	t1450->SetPixelChargeThreshold(2.0);
	assert(t1450->m_RcdaqFileNameMask);
	char fname[1024];
	snprintf(fname, 1024-1, t1450->m_RcdaqFileNameMask->Data(), fnum);
	auto it = new DREAMfileEventiterator(fname);
	if (!it) {
     			 printf("Failed to open input file '%s'\n", fname);
      			 exit(1);
                }
	it->DeclareWrappedPacket(3333, -1);
	it->InitializeMappingTables("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/alibordi/Test_Beam/mpgd4eic/mapping/database");
	auto uc = it->GetChamber("D1"); assert(uc);
        auto dc = it->GetChamber("D4"); assert(dc);
	snprintf(fname, 1024-1, t1450->m_OutputFileNameMask->c_str(), fnum);
	auto fout = new TFile(fname, "RECREATE");
	auto tout = new TTree("t", "LAPPD Tree");
	unsigned ic; 
	double tt[16], aa[16], trk_tx, trk_ty, baseline_val[16];//dt, t0, t1,  *taddr[2] = {&t0, &t1};//
	double baseline0, baseline1, *baddr[2] = {&baseline0, &baseline1};
        //double a0, a1, *aaddr[2] = {&a0, &a1}, tx, ty;
	bool ok[16];
	int  bestpeak_position[16], trk_cluster_size[4], trk_hits_first[2][4], trk_hits_second[2][4];
	double trig_t0,trig_t1,trig_t2,trig_t3, *trig_taddr[4] = {&trig_t0, &trig_t1,&trig_t2,&trig_t3};
        double trig_a0, trig_a1, trig_a2, trig_a3, *trig_aaddr[4] = {&trig_a0, &trig_a1,&trig_a2,&trig_a3};






        //=====Trigger====//

        tout->Branch("trig_t0", &trig_t0,"trig_t0/D");
        tout->Branch("trig_t1", &trig_t1,"trig_t1/D");
   	tout->Branch("trig_t2", &trig_t2,"trig_t2/D");
   	tout->Branch("trig_t3", &trig_t3,"trig_t3/D");
  	tout->Branch("trig_a0", &trig_a0,"trig_a0/D");
   	tout->Branch("trig_a1", &trig_a1,"trig_a1/D");
   	tout->Branch("trig_a2", &trig_a2,"trig_a2/D");
   	tout->Branch("trig_a3", &trig_a3,"trig_a3/D");

       //=======Planacon======//

	tout->Branch("tt", tt,"tt[16]/D");
        tout->Branch("ic", &ic,"ic/I");
        tout->Branch("aa",  aa,"aa[16]/D");
        tout->Branch("ok",  ok,"ok[16]/D");
	tout->Branch("baseline_val", &baseline_val, "baseline_val[16]/D");
        tout->Branch("bestpeak_position", &bestpeak_position, "bestpeak_position[16]/I");
        
        //===============DRS4 Cherenkov============//
        
        double charge0_chrnkv, charge1_chrnkv, charge2_chrnkv, charge3_chrnkv;
	double baseline_chrnkv;
	double *chrnkv_chargeaddr[4] = {&charge0_chrnkv, &charge1_chrnkv, &charge2_chrnkv, &charge3_chrnkv};

	tout->Branch("charge0_chrnkv", &charge0_chrnkv,"charge0_chrnkv/D");
	tout->Branch("charge1_chrnkv", &charge1_chrnkv,"charge1_chrnkv/D");
	tout->Branch("charge2_chrnkv", &charge2_chrnkv,"charge2_chrnkv/D");
	tout->Branch("charge3_chrnkv", &charge3_chrnkv,"charge3_chrnkv/D");
	tout->Branch("baseline_chrnkv", &baseline_chrnkv,"baseline_chrnkv/D");
	

        //============LAPPD: L03C board============//

	double planacon_timing_L03C_allpixel[24][24], tdiff_L03C_allpixel[24][24], baseline_L03C_allpixel[24][24], amplitude_L03C_allpixel[24][24];
	double amplitude_L03C_allpixel_for2D[24][24], tracelength_L03C_allpixel[24][24];
	int cluster_L03C_pixelcount;
	double cluster_L03C_amplitude, cluster_L03C_X, cluster_L03C_Y, cluster_L03C_R, cluster_L03C_Time;


	tout->Branch("planacon_timing_L03C_allpixel", &planacon_timing_L03C_allpixel,"planacon_timing_L03C_allpixel[24][24]/D");
	tout->Branch("tdiff_L03C_allpixel", &tdiff_L03C_allpixel,"tdiff_L03C_allpixel[24][24]/D");
	tout->Branch("baseline_L03C_allpixel", &baseline_L03C_allpixel,"baseline_L03C_allpixel[24][24]/D");
	tout->Branch("amplitude_L03C_allpixel", &amplitude_L03C_allpixel,"amplitude_L03C_allpixel[24][24]/D");
	tout->Branch("tracelength_L03C_allpixel", &tracelength_L03C_allpixel,"tracelength_L03C_allpixel[24][24]/D");
	tout->Branch("amplitude_L03C_allpixel_for2D", &amplitude_L03C_allpixel_for2D,"amplitude_L03C_allpixel_for2D[24][24]/D");
	tout->Branch("cluster_L03C_pixelcount", &cluster_L03C_pixelcount,"cluster_L03C_pixelcount/I");
	tout->Branch("cluster_L03C_amplitude", &cluster_L03C_amplitude,"cluster_L03C_amplitude/D");
	tout->Branch("cluster_L03C_X", &cluster_L03C_X,"cluster_L03C_X/D");
	tout->Branch("cluster_L03C_Y", &cluster_L03C_Y,"cluster_L03C_Y/D");
	tout->Branch("cluster_L03C_R", &cluster_L03C_R,"cluster_L03C_R/D");
	tout->Branch("cluster_L03C_Time", &cluster_L03C_Time,"cluster_L03C_Time/D");





        //========Tracker=====================//
	


	tout->Branch("trk_tx", &trk_tx,"trk_tx/D");
        tout->Branch("trk_ty", &trk_ty,"trk_ty/D");
	tout->Branch("trk_cluster_size", &trk_cluster_size, "trk_cluster_size[4]/I");
	tout->Branch("trk_hits_first",&trk_hits_first, "trk_hits_first[2][4]/I");
	tout->Branch("trk_hits_second",&trk_hits_second, "trk_hits_second[2][4]/I");







	std::cout<<"===================== size of the digitizers | " <<t1450->m_V1742s.size()<< " | =========================== "<<"\n";	



	assert(t1450->m_PixelMatrices.size() == 1);
        auto gmtx = t1450->m_PixelMatrices[0];
        unsigned gdimX = gmtx->GetDimX(), gdimY = gmtx->GetDimY();
        std::cout<<" The pixel matrix size : X-coor |"<< gdimX << "| \t and : Y-coor |"<< gdimY<<"\n";


	for(unsigned iv=0; iv<t1450->m_V1742s.size(); iv++) {
					auto v1742 = t1450->m_V1742s[iv];
					for(unsigned ich=0; ich<32; ich++) {
							auto channel = v1742->m_Channels[ich];
							TString hname; hname.Form("wf-v%d-ch%02d", iv, ich);
							channel->m_WFH = new TH1D(hname.Data(), hname.Data(), _1024_, 0, _1024_);
							//tout->Branch( hname.Data(), &channel->m_WFH);
		


                     }//for ich
					for(unsigned tr=0; tr<4; tr++) {
							 auto trigger = v1742->FastTrigger(tr);
							 TString hname; hname.Form("wf-v%d-tr%d", iv, tr);
							 trigger->m_WFH = new TH1D(hname.Data(), hname.Data(), _1024_, 0, _1024_);
							 //tout->Branch( hname.Data(), &trigger->m_WFH);


                     }//for tr
                }//iv

	
	for(unsigned ix=0; ix<gdimX; ix++){
		for(unsigned iy=0; iy<gdimY; iy++) {
   		TString hname; hname.Form("wf-%s-x%02d-y%02d", gmtx->GetName(), ix, iy);
   
  		 auto channel_L03C = gmtx->GetChannel(ix, iy);
   		if (!channel_L03C) continue;
   
   		channel_L03C->m_WFH = new TH1D(hname.Data(), hname.Data(), _1024_, 0, _1024_);
   
   		//tout->Branch( hname.Data(), &channelchannel_L03C->m_WFH);
   
   		channel_L03C->m_WF = new WaveForm((hname + '@').Data());
   		channel_L03C->m_WF->AddSourceHistogram(channel_L03C->m_WFH);
}
}


//===========================================================================================================================Wave form gathering stops here 
//Start run over events and acquire the data 
//1. Trigger 
//2. Planacon
//3. Cherenkov beam line 
//4. LAPPD L03C board 
//5. Tracker : GEM 
//===========================================================================================================================While loop starts 






	Event *evt;
	unsigned evCounter = 0;
	while ( (evt = /*t1450->m_RcdaqIterator->*/it->getNextEvent()) ) {
		if ( evt->getEvtType() != 1) continue;
			printf("\n@@@ event# %6d\n", evt->getEvtSequence());
			t1450->FillV1742s(evt);

//============================================================================================================================Trigger time 


	for(unsigned iv=0; iv<t1450->m_V1742s.size(); iv++) {
          auto v1742 = t1450->m_V1742s[iv];

          for(unsigned tr=0; tr<4; tr++) {
            auto ftrigger = v1742->FastTrigger(tr);
            int i0 = ftrigger->ThresholdCrossing(t1450->m_NIMPulseFitSetPoint);
            //printf("@@@ V1742 %2d, trigger WF %d --> threshold crossing @ sample %4d\n",       iv, tr, i0);
            assert(i0 != -1);

          
              int ifrom = i0 + t1450->m_NIMPulseBaselineWindowOffset;
              int ito = ifrom + int(t1450->m_NIMPulseBaselineWindowWidth) - 1;
              //std::cout<< " What is the value of the ifrom and ito   :"<<ifrom<<"\t"<<ito<<"\n";
	      assert(ifrom >=0 && ito > ifrom);
              double baseline = ftrigger->CalculateBaseline(ifrom, ito);
              //printf("@@@  NIM pulse %d baseline estimate: %7.2f [mV]\n", tr, baseline);
	      int ipeak = ftrigger->SimplePeakSearch(ifrom, ito);
	      double amplitude = ftrigger->m_Data[ipeak] - baseline;
	     // printf("@P@ NIM pulse [%d] peak amplitude: %7.2f [mV]\n", tr, amplitude);
	      
	      *trig_taddr[tr] = ftrigger->m_TimingReference =
                ftrigger->LeadingEdgeRangeFit(t1450->m_NIMPulseFitRange[0] + baseline,
                                              t1450->m_NIMPulseFitRange[1] + baseline,
                                              1, t1450->m_NIMPulseFitSetPoint + baseline);

      	     //printf("@@@   NIM pulse %d precise timing reference: %9.2f [ps]\n", tr, ftrigger->m_TimingReference);           
	     *trig_aaddr[tr] = -amplitude;


         } //for tr
        } //for iv










//======================================================================================================================================For planacon

 
#if 1
	for(unsigned pl=0; pl<16; pl++) {
                ok[pl] = false;
		auto v1742 = t1450->m_V1742s[0];
                auto planacon = v1742->GetChannel(pl);
		double baseline = planacon->CalculateBaseline(0, 10);
               // printf("@P@  Planacon %d baseline assumed : %7.2f [mV]\n", pl, baseline);
                int i0 = planacon->ThresholdCrossing(t1450->m_PlanaconPulseSearchThreshold + baseline);
                //printf("@P@  Planacon %d threshold crossing @ sample %4d\n", pl, i0);
		if (i0 == -1) continue;
		int ifrom = i0 + t1450->m_PlanaconPulseBaselineWindowOffset;
                if (ifrom < 0) continue;//goto _next_event;
                int ito = ifrom + int(t1450->m_PlanaconPulseBaselineWindowWidth) - 1;
		double baseline_pl = planacon->CalculateBaseline(ifrom, ito);
                baseline_val[pl] = baseline_pl;
                //printf("@p@  Planacon %d baseline estimate: %7.2f [mV]\n", pl, baseline_pl);
		int jfrom = i0;
                int jto = jfrom + t1450->m_PlanaconPulsePeakWindowWidth - 1;
                int jpeak = planacon->SimplePeakSearch(jfrom, jto);
		bestpeak_position[pl] = jpeak;
                //printf("@P@  Planacon %d peak position: %4d [ismp]\n", pl, jpeak);

                double amplitude = planacon->m_Data[jpeak] - baseline_pl;
                //printf("@P@  Planacon %d peak amplitude: %7.2f [mV]\n", pl, amplitude);
		aa[pl] = -amplitude;
#if 1
                tt[pl] = planacon->m_TimingReference =
                planacon->LeadingEdgeRangeFit(   t1450->m_PlanaconPulseFitRange[0]*amplitude    + baseline_pl,
                                               t1450->m_PlanaconPulseFitRange[1]*amplitude    + baseline_pl,
                                            1, t1450->m_PlanaconPulseFitSetPoint*amplitude    + baseline_pl);
#else

		*taddr[pl] = planacon->m_TimingReference = 
              planacon->LeadingEdgeSinglePointFit(t1450->m_PlanaconPulseFitSetPoint*amplitude + baseline);
#endif
		//printf("@P@  Planacon %d timing: %7.2f [ps]\n", pl, planacon->m_TimingReference);


#if _NIM_
   {
   int ifrom = i0 + t1450->m_NIMPulseBaselineWindowOffset;
   int ito = ifrom + int(t1450->m_NIMPulseBaselineWindowWidth) - 1;
   assert(ifrom >=0 && ito > ifrom);
   double baseline = planacon->CalculateBaseline(ifrom, ito);
   //printf("@N@  NIM pulse %d baseline estimate: %7.2f [mV]\n", pl, baseline);
   *taddr[pl] = planacon->m_TimingReference =
   planacon->LeadingEdgeRangeFit(t1450->m_NIMPulseFitRange[0]+baseline,
                                 t1450->m_NIMPulseFitRange[1]+baseline,
                                 1, t1450->m_NIMPulseFitSetPoint+baseline);
   //printf("@N@  NIM pulse %d precise timing reference: %9.2f [ps]\n", pl, planacon->m_TimingReference);
   }
#endif

   
#if _OLD_
      
   int ifrom = iNIMref + t1450->m_PlanaconPeakWindowOffset;
   int ito = ifrom + t1450->m_PlanaconPeakWindowWidth - 1;
   int ipeak = planacon->SimplePeakSearch(ifrom, ito);
   //printf("@@@ pl#%d peak sample: %d\n", pl, ipeak);
   if (ipeak < 50) continue;
   
     
   int jfrom = ipeak + t1450->m_PlanaconBaselineWindowOffset;
   int jto = jfrom + t1450->m_PlanaconBaselineWindowWidth - 1;
   double baseline = *baddr[pl] = planacon->CalculateBaseline(jfrom, jto);
   //printf("@@@  Planacon %d baseline estimate: %7.2f [mV]\n", pl, baseline);
   
      
   double amplitude = planacon->m_Data[ipeak] - baseline;
   //printf("@@@  Planacon %d amplitude        : %7.2f [mV]\n", pl, amplitude);
   *aaddr[pl] = -amplitude;
   
   planacon->m_TimingReference = *taddr[pl] =
   planacon->LeadingEdgeRangeFit(t1450->m_PlanaconPulseFitRange[0]*amplitude+baseline,
                                 t1450->m_PlanaconPulseFitRange[1]*amplitude+baseline,
                                 1, t1450->m_PlanaconPulseFitSetPoint*amplitude+baseline);
      
   //printf("@@@  Planacon %d timing reference : %7.2f [ps]\n", pl, planacon->m_TimingReference);
#endif
   
      
   ok[pl] = true;
} //for planacon
#endif




//==========================================================================================================================Cherenkov Beam line 

// Beam line Cherenkov DRS4; recycle part of the V1742 codes; do not mind to
// hardcode the parameters right here;
// Martin's equipment in the DAQ;

const unsigned cdim = 4;
std::vector<V1742channel*> Cherenkov;
for(unsigned ich=0; ich<cdim; ich++)Cherenkov.push_back(new V1742channel(77, ich));
   
{
      
   unsigned eq_chrnkv = 1020;
      
   unsigned blfrom = 30, blto = 50, pfrom = 70, pto = 90;// Baseline and peak evaluation ranges; assume the same for all channels;
   
   Packet *p_chrnkv = evt->getPacket(eq_chrnkv);
   
   if (p_chrnkv) {
     
      
      for(unsigned ich=0; ich<cdim; ich++) {
         auto channel_chrnkv = Cherenkov[ich];
         
         for(unsigned ismp=0; ismp<1024; ismp++) {
            channel_chrnkv->m_Data[ismp] = p_chrnkv->rValue(ismp, ich);
            if (!ich && !ismp)
               printf("@F@ %d %4d -> %7.2f\n", ich, ismp, channel_chrnkv->m_Data[ismp]);
            
            if (channel_chrnkv->m_WFH) channel_chrnkv->m_WFH->SetBinContent(ismp+1, channel_chrnkv->m_Data[ismp]);
         } //for ismp
         
         
          baseline_chrnkv = channel_chrnkv->CalculateBaseline(blfrom, blto);
   double charge_chrnkv   = channel_chrnkv->CalculateCharge(pfrom, pto, baseline_chrnkv);
         
         *chrnkv_chargeaddr[ich] = charge_chrnkv;
      } //for ich
      
      delete p_chrnkv;
   } else
      printf("@@@ No DRS4 eval board data packet found\n");
      
      
}

//==================================================================================================================L03C board




{
   for(unsigned ix=0; ix<gmtx->GetDimX(); ix++) {
      for(unsigned iy=0; iy<gmtx->GetDimY(); iy++) {
         auto channel_L03C = gmtx->GetChannel(ix, iy);
         if (!channel_L03C) continue;
         //amplitude_L03C_allpixel_for2D[ix+1][iy+1] = channel_L03C->m_DataSum;
         auto v1742_L03C = t1450->m_V1742s[channel_L03C->m_V1742]; assert(v1742_L03C);
         if (channel_L03C == t1450->m_Planacon[0] || channel_L03C == t1450->m_Planacon[1]) continue;
         auto t0_L03C = t1450->m_Planacon[1];
         double planacon_timing_L03C = t0_L03C->m_TimingReference - t0_L03C->m_FastTrigger->m_TimingReference;
         planacon_timing_L03C_allpixel[ix][iy] = planacon_timing_L03C;
         double tdiff_L03C = planacon_timing_L03C + channel_L03C->m_FastTrigger->m_TimingReference + 1.0*channel_L03C->GetExternalDelay();
	 tdiff_L03C_allpixel[ix][iy] = tdiff_L03C;
         int ifrom = (int)rint(tdiff_L03C / v1742_L03C->m_NominalSampleWidth) + t1450->m_RegularChannelPeakWindowOffset;
         int ito = ifrom +  t1450->m_RegularChannelPeakWindowWidth - 1;
         printf("@C@ ==========================================Pixel [%2d,%2d] signal window  : %4d .. %4d\n", ix, iy, ifrom, ito);
            
         if (ifrom < 0 || ito > 1023) goto _next_event;
         
         double attenuation_L03C = 1.0;//exp(-channel->GetTraceLength()/1500.); //printf("%f\n", exp(-380./1000.)); exit(0);
        
                                 
         {
         int kpeak = channel_L03C->SimplePeakSearch(ifrom, ito);
         printf("@C@ Pixel [%2d,%2d] peak WF channel   :  %4d (dist %3d), %7.2f [mV]\n", ix, iy, kpeak, kpeak - ifrom, channel_L03C->m_Data[kpeak]);
         
         int kfrom = kpeak + t1450->m_RegularChannelBaselineWindowOffset;
         if (kfrom < 0 ) kfrom = 0;
         int kto = kfrom + t1450->m_RegularChannelBaselineWindowWidth;
         double baseline_L03C = channel_L03C->CalculateBaseline(kfrom, kto);
         baseline_L03C_allpixel[ix][iy] = baseline_L03C;
         printf("@C@ Pixel [%2d,%2d] baseline : %7.2f [mV]\n", ix, iy, baseline_L03C);
         
#if 0
         double amplitude_L03C = channel_L03C->m_DataSum = channel_L03C->CalculateCharge(kpeak-1, kpeak+1, baseline_L03C) /(attenuation_L03C * (3));
         
#else
           
         double amplitude_L03C = channel_L03C->m_DataSum = -(channel_L03C->m_Data[kpeak] - baseline_L03C) / attenuation_L03C;
#endif
         printf("@C@ Pixel [%2d,%2d] amplitude estimate:  %7.2f [mV]\n@Q@\n", ix, iy, amplitude_L03C);
         amplitude_L03C_allpixel[ix][iy] = amplitude_L03C;
         
         
         }
      } //for ix
   } //for iy
   
  
   
    
   auto clusters_L03C = gmtx->GetPixelClusters(t1450->m_PixelChargeThreshold);
   
   unsigned order = 0;
   for(auto cptr: clusters_L03C) {
      auto &cluster_L03C = cptr.second;
      
        
         auto top = cluster_L03C.GetCentralPixel();
         auto channel_L03C = top->GetDigitizerChannel();
         
      cluster_L03C_Time = channel_L03C->m_TimingReference;
      
         
      cluster_L03C.Calculate(gdimX, gdimY, mtx->GetPitch());
      
      printf("@X@ LAPPD cluster: %3d pixel(s), ampl: %7.2f; X = %7.1f [mm], Y = %7.1f [mm] -> R = %7.1f [mm]\n",
             cluster_L03C.GetPixelCount(), cluster_L03C.GetAmplitude(), cluster_L03C.GetX(), cluster_L03C.GetY(),
             sqrt(pow(cluster_L03C.GetX(), 2) + pow(cluster_L03C.GetY(), 2)));
      
      cluster_L03C_pixelcount = cluster_L03C.GetPixelCount();
      cluster_L03C_amplitude = cluster_L03C.GetAmplitude();
      cluster_L03C_X = cluster_L03C.GetX();
      cluster_L03C_Y = cluster_L03C.GetY();
      cluster_L03C_R = sqrt(pow(cluster_L03C.GetX(), 2) + pow(cluster_L03C.GetY(), 2));
      
      
      
      
         
   } //for ptr
   
}





//==================================================================================================================for tracker 

ic = t1450->m_V1742s[0]->m_IndexCell[0];
   
#if 1
    
{
    
   auto p = it->getPacket(evt, 3333);
   if (p) {
      DreamEvent event(it, p);
      
      {
      std::map<unsigned, int> _qhits[2][4];
      
      for(auto entry: event.GetEntries()) {
         unsigned feu = entry.feu(), sample = entry.tsample();
         unsigned channel = entry.channel();
         int ampl = entry.amplitude() - 256;
         switch (feu) {
            case 100:
            {
            auto geo = uc->GetPlaneAndStrip(channel);
            int plane = geo.first, strip = geo.second;
            _qhits[0][plane][channel] += ampl;
            }
               break;
            case 11:
                  
            {
            auto geo = dc->GetPlaneAndStrip(channel);
            int plane = geo.first, strip = geo.second;
            _qhits[1][plane][channel] += ampl;
            }
               break;
            default:
               ;
         } //switch
      } //for entry
        
      for(unsigned ud=0; ud<2; ud++)
         for(unsigned iq=0; iq<4; iq++) {
            unsigned ixy = iq % 2, ieo = (iq/2)%2;
            auto plane = (ud ? dc : uc)->GetPlane(iq);
            std::vector<std::pair<unsigned, int> > hits;
            for(auto qhit: _qhits[ud][iq])
               hits.push_back(std::make_pair(qhit.first, qhit.second));
               for(int i = 0; i < hits.size(); i++)
                     {
			trk_hits_first[ud][iq] = hits[i].first; trk_hits_second[ud][iq] = hits[i].second;
                     }
            auto clusters = plane->GetClusters(hits);
            trk_cluster_size[iq] = (int)clusters.size();
            printf("@D@ %d the cannel %2d cluster(s)\n", iq, (int)clusters.size());
            std::map<double, Cluster*> reordered;
            for(auto &cluster: clusters)
               reordered[cluster.mAmplitude] = &cluster;
            if (!ud && !ieo && reordered.size())
               (ixy ? trk_ty : trk_tx) = reordered.rbegin()->second->mCentroid;
         } //for ud..iq
      }
      it->deletePacket(evt, 3333, p);
   }
   else
      printf("@@@ No DREAM data packet found\n");
      }
#endif







//======================================================================================================================Data aqcuired 



tout->Fill();
_next_event:
evCounter++; if (t1450->m_Stat && evCounter == t1450->m_Stat) break;


}//while 

tout->Write();
fout->Close();
delete fout;
exit(0);


}
