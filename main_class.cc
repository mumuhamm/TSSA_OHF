#include <PHGlobal.h>
#include <ReactionPlaneObject.h>
#include <PHMuoTracksOut.h>
#include <RunHeader.h>
#include <MUTOO.h>
#include <Fun4AllReturnCodes.h>
#include <Fun4AllServer.h>
#include <RunNumberRanges.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <MWG.h>
#include <MWGConsts.h>
#include <MWGVersion.h>
#include <TMuiPseudoLL1Map.h>
#include <mMuiFastRoadFinderPar.h>
#include <TMuiGeometry.hh>
#include <PHTrackIntegratorKF.h>
#include <Tools.h>
#include <TMutTrkMap.h>
#include <TMutTrkPar.hh>
#include <mMutFitVtx.h>
#include <MuonUtil.h>
#include <mMfmMT.h>
#include <MWGVertex.h>
#include <RunHeader.h>
#include <SyncObject.h>
#include <SingleMuon.h>
#include <SingleMuon_v1.h>
#include <SingleMuon_v2.h>
#include <SingleMuon_v3.h>
#include <SingleMuon_v4.h>
#include <SingleMuon_v5.h>
#include <SingleMuon_v6.h>
#include <SingleMuon_v7.h>
#include <SingleMuon_v8.h>
#include <SingleMuon_v9.h>
#include <SingleMuon_v10.h>
#include <SingleMuon_v11.h>
#include <SingleMuon_v12.h>
#include <SingleMuon_v15.h>
#include <SingleMuon_v16.h>
#include <SingleMuon_v17.h>
#include <SingleMuonContainer.h>
#include <SingleMuonContainer_v1.h>
#include <SingleMuonContainer_v2.h>
#include <SingleMuonContainer_v3.h>
#include <SingleMuonContainer_v4.h>
#include <SingleMuonContainer_v5.h>
#include <SingleMuonContainer_v6.h>
#include <SingleMuonContainer_v7.h>
#include <SingleMuonContainer_v8.h>
#include <SingleMuonContainer_v9.h>
#include <SingleMuonContainer_v10.h>
#include <SingleMuonContainer_v11.h>
#include <SingleMuonContainer_v12.h>
#include <SingleMuonContainer_v13.h>
#include <SingleMuonContainer_v14.h>
#include <SingleMuonContainer_v15.h>
#include <SingleMuonContainer_v16.h>
#include <SingleMuonContainer_v17.h>
#include <VtxOut.h>
#include <TRpcTrkMap.h>
#include <TRpcMuoTrkMap.h>
#include <BbcRaw.h>
#include <FvtxConeTracklets.h>
   //#include <TFvtxCompactCoordMap.h>
#include <TFvtxCompactTrkMap.h>
#include <SvxSegmentList.h>
#include <SvxSegment.h>

#include <TH1F.h>
#include "TSpectrum.h"
#include <TF1.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>
#include <PHGslMatrix.h>
#include <TMath.h>

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include "mFillSingleMuonContainer.h"

using namespace std;
typedef PHIODataNode<PHObject> PHObjectNode_t;

   //___________________________________________________________________
int mFillSingleMuonContainer::Init(PHCompositeNode *top_node)
{
      //create output node in the node tree
   PHNodeIterator iter(top_node);
   PHCompositeNode *dstNode
   = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode","DST"));
   
   if(!dstNode) {
      cout << "dstNode not found"<<endl;
   } else {
      cout << "dstNode is found"<<endl;
   }
   
   SingleMuonContainer *muons;
   switch (sm_version)
      {
         case 1: muons = new SingleMuonContainer_v1(); break;
         case 2: muons = new SingleMuonContainer_v2(); break;
         case 3: muons = new SingleMuonContainer_v3(); break;
         case 4: muons = new SingleMuonContainer_v4(); break;
         case 5: muons = new SingleMuonContainer_v5(); break;
         case 6: muons = new SingleMuonContainer_v6(); break;
         case 7: muons = new SingleMuonContainer_v7(); break;
         case 8: muons = new SingleMuonContainer_v8(); break;
         case 9: muons = new SingleMuonContainer_v9(); break;
         case 10: muons = new SingleMuonContainer_v10(); break;
         case 11: muons = new SingleMuonContainer_v11(); break;
         case 12: muons = new SingleMuonContainer_v12(); break;
         case 13: muons = new SingleMuonContainer_v13(); break;
         case 14: muons = new SingleMuonContainer_v14(); break;
         case 15: muons = new SingleMuonContainer_v15(); break;
         case 16: muons = new SingleMuonContainer_v16(); break;
         case 17: muons = new SingleMuonContainer_v17(); break;
         default: muons = new SingleMuonContainer_v17(); break;
      }
   
   if(muons && dstNode) {
      PHObjectNode_t *muonNode =
      new PHIODataNode<PHObject>(muons,"SingleMuonContainer","PHObject");
      dstNode->addNode(muonNode);
      cout << "SingleMuonContainer is added" <<endl;
   } else {
      cout << ThisName << " Init() failed to create output object"<<endl;
      return ABORTRUN;
      
   }
   
   nevents = 0;
   naccepted_N = 0;
   naccepted_S = 0;
   
   if(is_sim)
      {
      set_vtx_vertex_name("SIM");
      cout << "mFillSingleMuonContainer::Init() - Setting vtx name to SIM!" << endl;
      }
   if(!_vtx_top_node.empty())
      {
      Fun4AllServer *se = Fun4AllServer::instance();
      vtxTopNode = se->topNode(_vtx_top_node.c_str());
      cout << "mFillSingleMuonContainer::Init() - Setting vtx top_node to " << _vtx_top_node << endl;
      if(!vtxTopNode) return ABORTRUN;
      }
   
      // used to find if track fired the trigger
   try
   {
   TMutNode<TMuiHitMapO>::find_node(top_node,"TMuiHitMapO_track");
   }
   catch(std::exception& e)
   {
   TMutNode<TMuiHitMapO>::new_node(top_node,"TMuiHitMapO_track");
   }
   
   simMuIDLl1->set_tmuihitmap_name("TMuiHitMapO_track");
   
   return 0;
} //  End of Init()

   //______________________________________________________
int mFillSingleMuonContainer::InitRun(PHCompositeNode *top_node)
{
      // load field for track projections
   static bool initialized = false;
   if(!initialized)
      {
         //      mMfmMT::setMapFileScale(0.98);
      MuonUtil::initialize_database(top_node, true, true);
      initialized = true;
      }
   
      // Make sure global db parameters are set:
   TFvtxGlobalParCntrl::init_run();
   
   MUTOO::PRINT( cout, "mFillSingleMuonContainer::InitRun" );
   TFvtxGlobalParCntrl::print();
   
   
   return 0;
}

   //______________________________________________________
int mFillSingleMuonContainer::process_event(PHCompositeNode *top_node)
{
   
   if(nevents%100000==0){
      cout << "Event: " << nevents << endl;
   }
   nevents++;
   
   int Run_Number = 0;
   
   RunHeader* runh = findNode::getClass<RunHeader>(top_node,"RunHeader");
   if (runh)
      Run_Number = runh->get_RunNumber();
   
      // PHGlobal
   PHGlobal* evt = findNode::getClass<PHGlobal>(top_node, "PHGlobal" );
   if (!evt)
      {
         //      cout << PHWHERE << "mFillSingleMuonContainer:: PHGlobal not in Node Tree" << endl;
      return ABORTRUN;
      }
   
   BbcRaw * bbcraw = findNode::getClass<BbcRaw>(top_node, "BbcRaw" );
   if (!bbcraw)
      {
      cout << PHWHERE << "mFillSingleMuonContainer:: BbcRaw not in Node Tree" << endl;
      return ABORTRUN;
      }
   
      // VTX
   if(_vtx_top_node.empty()) vtx = findNode::getClass<VtxOut>(top_node, "VtxOut" );
   else vtx = findNode::getClass<VtxOut>(vtxTopNode, "VtxOut" );
   if(!vtx){ cout << "Could not find VtxOut on " << vtxTopNode->getName() << endl; return ABORTRUN;}
   
      // retrieve tracks, check if empty
   TMutTrkMap *trk_map = findNode::getClass<TMutTrkMap>(top_node,"TMutTrkMap");
   
      // new framework MWG tracks
   PHMuoTracksOut* muo = findNode::getClass<PHMuoTracksOut>(top_node,"PHMuoTracksOO");
   if (!muo)
      {
         //      cout << "mFillSingleMuonContainer:: PHMuoTracksOO (new framework) not in Node Tree" << endl;
      return DISCARDEVENT;
      }
   if (nevents<2)
      muo->ShutUp();
   
   SingleMuonContainer *muons;
   switch (sm_version)
      {
         case 1: muons = findNode::getClass<SingleMuonContainer_v1>(top_node, "SingleMuonContainer" ); break;
         case 2: muons = findNode::getClass<SingleMuonContainer_v2>(top_node, "SingleMuonContainer" ); break;
         case 3: muons = findNode::getClass<SingleMuonContainer_v3>(top_node, "SingleMuonContainer" ); break;
         case 4: muons = findNode::getClass<SingleMuonContainer_v4>(top_node, "SingleMuonContainer" ); break;
         case 5: muons = findNode::getClass<SingleMuonContainer_v5>(top_node, "SingleMuonContainer" ); break;
         case 6: muons = findNode::getClass<SingleMuonContainer_v6>(top_node, "SingleMuonContainer" ); break;
         case 7: muons = findNode::getClass<SingleMuonContainer_v7>(top_node, "SingleMuonContainer" ); break;
         case 8: muons = findNode::getClass<SingleMuonContainer_v8>(top_node, "SingleMuonContainer" ); break;
         case 9: muons = findNode::getClass<SingleMuonContainer_v9>(top_node, "SingleMuonContainer" ); break;
         case 10: muons = findNode::getClass<SingleMuonContainer_v10>(top_node, "SingleMuonContainer" ); break;
         case 11: muons = findNode::getClass<SingleMuonContainer_v11>(top_node, "SingleMuonContainer" ); break;
         case 12: muons = findNode::getClass<SingleMuonContainer_v12>(top_node, "SingleMuonContainer" ); break;
         case 13: muons = findNode::getClass<SingleMuonContainer_v13>(top_node, "SingleMuonContainer" ); break;
         case 14: muons = findNode::getClass<SingleMuonContainer_v14>(top_node, "SingleMuonContainer" ); break;
         case 15: muons = findNode::getClass<SingleMuonContainer_v15>(top_node, "SingleMuonContainer" ); break;
         case 16: muons = findNode::getClass<SingleMuonContainer_v16>(top_node, "SingleMuonContainer" ); break;
         case 17: muons = findNode::getClass<SingleMuonContainer_v17>(top_node, "SingleMuonContainer" ); break;
         default: muons = findNode::getClass<SingleMuonContainer_v17>(top_node, "SingleMuonContainer" ); break;
      }
   
   if (!muons)
      {
      cout << "mFillSingleMuonContainer:: SingleMuonContainer not in Node Tree" << endl;
      return ABORTRUN;
      }
   muons->Reset();
   
      // Z Vertex and Centrality info
   if (evt){
      muons->set_Evt_bbcZ(evt->getBbcZVertex());
      muons->set_Evt_bbcZ_Err(evt->getBbcZVertexError());
      
      muons->set_Evt_Cent(evt->getCentrality());
   }
      //-- debug
      //  if ( muons->get_Evt_bbcZ() < -200.0 )  cout << " ====>>> muons->get_Evt_bbcZ(" <<  nevents << ") = " <<  muons->get_Evt_bbcZ() << endl;
   
   
   if (bbcraw)
      {
      TH1F *htdcS = new TH1F("htdcS", "htdcS", 50, 0, 3600);
      TH1F *htdcN = new TH1F("htdcN", "htdcN", 50, 0, 3600);
      TSpectrum *specS = new TSpectrum();
      TSpectrum *specN = new TSpectrum();
      
      for (int ipmt=0; ipmt<bbcraw->get_npmt(); ipmt++)
         {
         if (bbcraw->get_Pmt(ipmt)<64)
            htdcS->Fill(bbcraw->get_Tdc0(ipmt));
         else
            htdcN->Fill(bbcraw->get_Tdc0(ipmt));
         }
      
      specS->Search(htdcS, 1, "", 0.1);
      specN->Search(htdcN, 1, "", 0.1);
      
      muons->set_Evt_BbcTdcMean_S(htdcS->GetMean());
      muons->set_Evt_BbcTdcMean_N(htdcN->GetMean());
      muons->set_Evt_BbcTdcRMS_S(htdcS->GetRMS());
      muons->set_Evt_BbcTdcRMS_N(htdcN->GetRMS());
      muons->set_Evt_BbcTdcMode_S(htdcS->GetBinCenter(htdcS->GetMaximumBin()));
      muons->set_Evt_BbcTdcMode_N(htdcN->GetBinCenter(htdcN->GetMaximumBin()));
      muons->set_Evt_BbcTdcnPeaks_S(specS->GetNPeaks());
      muons->set_Evt_BbcTdcnPeaks_N(specN->GetNPeaks());
      
      delete htdcS;
      delete htdcN;
      delete specS;
      delete specN;
      
      }
   
      //Multiplicity   -- identical selection of mFillDiMuonContainer   03/09/2021  MXL
   TFvtxCompactTrkMap* fvtx_trk_map = findNode::getClass<TFvtxCompactTrkMap>(top_node,"TFvtxCompactTrkMap");
   if(sm_version >= 13 && fvtx_trk_map)
      {
      int fvtxMult_N = 0;
      int fvtxMult_S = 0;
      
      float chi2_fvtx_cut  =20;
      float  chi2_fvtx_prob_cut = 0.05;
      int  min_nhits_fvtx =  2;
      
      TFvtxCompactTrkMap::iterator iter( fvtx_trk_map->range() );
      while( TFvtxCompactTrkMap::const_pointer fvtx_ptr = iter.next() )
         {
         if((*fvtx_ptr)->get_chi2_ndf() > 0 && (*fvtx_ptr)->get_chi2_ndf() < chi2_fvtx_cut)  // very loose cuts
            {
            int ndf = 2*(*fvtx_ptr)->get_nhits()-5;
            
            float chi2_prob = TMath::Prob( (*fvtx_ptr)->get_chi2_ndf()*ndf, ndf);
            if (chi2_prob >= chi2_fvtx_prob_cut && (*fvtx_ptr)->get_nhits() >=  min_nhits_fvtx)
               {
               if ((*fvtx_ptr)->get_fvtx_eta() > 0) fvtxMult_N++;
               if ((*fvtx_ptr)->get_fvtx_eta() < 0) fvtxMult_S++;
               }
            } //
         } // while loop
      
      /*
       while( TFvtxCompactTrkMap::const_pointer fvtx_ptr = iter.next() )
       {
       if((*fvtx_ptr)->get_chi2_ndf() > 0 && (*fvtx_ptr)->get_chi2_ndf() < 20)
       {
       if ((*fvtx_ptr)->get_fvtx_eta() > 0) fvtxMult_N++;
       if ((*fvtx_ptr)->get_fvtx_eta() < 0) fvtxMult_S++;
       }
       }
       */
      muons->set_Evt_Mult_FVTXN(fvtxMult_N);
      muons->set_Evt_Mult_FVTXS(fvtxMult_S);
      
         //-- debug
         //      cout << " ====>>>  muons->set_Evt_Mult_FVTXN(fvtxMult_N, S) (Evt=" <<  nevents << ") = " << fvtxMult_N  << "\t" <<fvtxMult_S <<  endl;
      }
      //  cout << "  ====>>>  muons->get_Evt_bbcZ() (Evt=" <<  nevents << ") = " << muons->get_Evt_bbcZ() <<  endl;
   
   
   SvxSegmentList* segmentList = findNode::getClass<SvxSegmentList>(top_node,"SvxSegmentList");
   if(sm_version >= 13 && segmentList)
      {
         //Apply cuts used in SvxPrimVertexFinder
      int svxMult = 0;
      for(int isvx=0; isvx<segmentList->get_nSegments(); ++isvx) {
         SvxSegment *segment = segmentList->get_segment(isvx);
         
         if (isnan(segment->get3Momentum(0))) continue;
         if (isnan(segment->get3Momentum(1))) continue;
         if (isnan(segment->get3Momentum(2))) continue;
         if (isnan(segment->getMomentum())) continue;
         
         if (isnan(segment->getInnerMostProjectedPosition(0))) continue;
         if (isnan(segment->getInnerMostProjectedPosition(1))) continue;
         if (isnan(segment->getInnerMostProjectedPosition(2))) continue;
         
         if ((segment->getChiSq()/segment->getNDF()) <= 0.0) continue;
         if ((segment->getChiSq()/segment->getNDF()) > 6.0) continue;
         
         if ( fabs(segment->getDCA2D())>2.0 ) continue;
         
         svxMult++;
      }
      muons->set_Evt_Mult_SVX(svxMult);
      }
   
   
   if (vtx)
      {
      muons->set_Evt_vtxX(vtx->get_Vertex(_vtx_vertex_name.c_str()).getX());
      muons->set_Evt_vtxY(vtx->get_Vertex(_vtx_vertex_name.c_str()).getY());
      muons->set_Evt_vtxZ(vtx->get_Vertex(_vtx_vertex_name.c_str()).getZ());
         //set vtx vertex error.
         //float vtx1 = vtx->get_VertexError(_vtx_vertex_name.c_str()).getX();
         //float vtx2 = vtx->get_VertexError(_vtx_vertex_name.c_str()).getX();
         //std::cout<<"******************  vtx1= "<<vtx1<<"  vtx2= "<<vtx2<<"  *************"<<std::endl;
      muons->set_Evt_vtxX_Err(vtx->get_VertexError(_vtx_vertex_name.c_str()).getX());
      muons->set_Evt_vtxY_Err(vtx->get_VertexError(_vtx_vertex_name.c_str()).getY());
      muons->set_Evt_vtxZ_Err(vtx->get_VertexError(_vtx_vertex_name.c_str()).getZ());
      
      if (sm_version > 2 )
         {
         muons->set_Evt_fvtxX(vtx->get_Vertex(_fvtx_vertex_names[0].data()).getX());
         muons->set_Evt_fvtxY(vtx->get_Vertex(_fvtx_vertex_names[0].data()).getY());
         muons->set_Evt_fvtxZ(vtx->get_Vertex(_fvtx_vertex_names[0].data()).getZ());
            //set fvtx vertex error.
            //float vtx1 = vtx->get_VertexError(_fvtx_vertex_names[0].data()).getX();
            //float vtx2 = vtx->get_VertexError(_fvtx_vertex_names[0].data()).getY();
            //std::cout<<"******************  vtx1= "<<vtx1<<"  vtx2= "<<vtx2<<"  *************"<<std::endl;
         muons->set_Evt_fvtxX_Err(vtx->get_VertexError(_fvtx_vertex_names[0].data()).getX());
         muons->set_Evt_fvtxY_Err(vtx->get_VertexError(_fvtx_vertex_names[0].data()).getY());
         muons->set_Evt_fvtxZ_Err(vtx->get_VertexError(_fvtx_vertex_names[0].data()).getZ());
         
         muons->set_Evt_fvtxX2(vtx->get_Vertex(_fvtx_vertex_names[1].data()).getX());
         muons->set_Evt_fvtxY2(vtx->get_Vertex(_fvtx_vertex_names[1].data()).getY());
         muons->set_Evt_fvtxZ2(vtx->get_Vertex(_fvtx_vertex_names[1].data()).getZ());
         }
      }
   
   
      //if (fabs(muons->get_Evt_bbcZ())>bbcz_cut && fabs(muons->get_Evt_vtxZ())>bbcz_cut) return DISCARDEVENT; //can't filter out events of nan vertex(?)
   if ( fabs(muons->get_Evt_bbcZ()) > bbcz_cut ) return DISCARDEVENT;
   
      //  if (muons->get_Evt_fvtxZ()<min_fvtxz_cut || muons->get_Evt_fvtxZ()>=max_fvtxz_cut) return ABORTEVENT;
   if (muons->get_Evt_fvtxZ()<min_fvtxz_cut || muons->get_Evt_fvtxZ()>=max_fvtxz_cut) return DISCARDEVENT;  // 03/06/2021 MXL
   
   
   
      //
   
   ReactionPlaneObject* rp = findNode::getClass<ReactionPlaneObject>(top_node, "ReactionPlaneObject");
   if (rp)
      {
         // Check whether rp is implemented to avoid invalid error messages
      const TString v = rp->ClassName();
      if ((v == TString("ReactionPlaneObject")) || (v == TString("ReactionPlaneObjectv4")))
         {
         static bool once = true;
         
         if (once)
            {
            
            cout << "mFillSingleMuonContainer::process_event - WARNING - Incompatible version of ReactionPlaneObject ("
            << v << ") received from the TOP node. Ignore it."
            << endl;
            
            once = false;
            }
         
         }
      else
         {
         muons->set_rx_Ninner(rp->getRXNrp13());
         muons->set_rx_Sinner(rp->getRXNrp10());
         muons->set_rx_Nouter(rp->getRXNrp14());
         muons->set_rx_Souter(rp->getRXNrp11());
         muons->set_rx_N(rp->getRXNrp15());
         muons->set_rx_S(rp->getRXNrp12());
         muons->set_rx_NS(rp->getRXNrp18());
         }
      }
   
      // RpcTrkMap
   rpctrk_map = NULL;
   try
   {
   rpctrk_map = TMutNode<TRpcTrkMap>::find_node(top_node, "TRpcTrkMap");
   }
   catch (exception &e)
   { // MUTOO::TRACE( e.what() );
      static bool once = true;
      
      if (once)
         {
         cout
         << "mFillSingleMuonContainer::process_event - RpcTrkMap not in Node Tree"
         << endl;
         }
      once = false;
      
      rpctrk_map = NULL;
   }
   
      //RPCMuoTrkMap
   rpc_muotrk_map = NULL;
   try
   {
   rpc_muotrk_map = TMutNode<TRpcMuoTrkMap>::find_node(top_node,
                                                       "TRpcMuoTrkMap");
   }
   catch (exception &e)
   {
   static bool once = true;
   
   if (once)
      {
      cout << PHWHERE << "mFillSingleMuonContainer - TRpcMuoTrkMap not in Node Tree"
      << endl;
         // MUTOO::TRACE( e.what() );
      }
   once = false;
   rpc_muotrk_map = NULL;
   }
   
   int npart = muo->get_npart();
   
   SingleMuon* muon;
   
   switch (sm_version)
      {
         case 1: muon = new SingleMuon_v1(); break;
         case 2: muon = new SingleMuon_v2(); break;
         case 3: muon = new SingleMuon_v3(); break;
         case 4: muon = new SingleMuon_v4(); break;
         case 5: muon = new SingleMuon_v5(); break;
         case 6: muon = new SingleMuon_v6(); break;
         case 7: muon = new SingleMuon_v7(); break;
         case 8: muon = new SingleMuon_v8(); break;
         case 9: muon = new SingleMuon_v9(); break;
         case 10: muon = new SingleMuon_v10(); break;
         case 11: muon = new SingleMuon_v11(); break;
         case 12: muon = new SingleMuon_v12(); break;
         case 13: muon = new SingleMuon_v12(); break; // yes that's right v12 here
         case 15: muon = new SingleMuon_v15(); break;
         case 16: muon = new SingleMuon_v16(); break;
         case 17: muon = new SingleMuon_v17(); break;
         default: muon = new SingleMuon_v17(); break;
      }
   
      // Loop over single muons to fill container
   for (int imu=0; imu<npart; imu++)
      {
      muon->Reset();
      
      float px = muo->get_px(0,imu);
      float py = muo->get_py(0,imu);
      float pz = muo->get_pz(0,imu);
      
      muon->set_px(px);
      muon->set_py(py);
      muon->set_pz(pz);
      
      if (fabs(muon->get_p()) < p_cut) continue;
      if (Tools::pT(muon->get_px(), muon->get_py()) < pt_cut ) continue;
      
         // get best road
      int iroad = get_best_road_oo( imu, muo );
      if ( iroad >= 0 )
         muon->set_idhits(muo->get_muIDOOhits(iroad, imu));
      
      for(int igap=4; igap>0; igap--)
         if (muo->is_muIDOOhit( iroad, imu, igap, 0) || muo->is_muIDOOhit( iroad, imu, igap, 1 ))
            {
            muon->set_lastgap(igap);
            break;
            }
      if (muon->get_lastgap() < lastgap_cut) continue;
      if (hadron_cut && muon->get_lastgap()>3) continue;
      
      muon->set_rapidity(Tools::rapidity(MU_MASS,
                                         muon->get_px(),
                                         muon->get_py(),
                                         muon->get_pz()));
      
      muon->set_trhits(muo->get_muTRhits(imu));
      muon->set_DG0(Tools::DG0( muo, imu, iroad ));
      if (muon->get_DG0() > DG0_cut && lastgap_cut>0) continue;
      if (muon->get_p()*muon->get_DG0() > pDG0_cut && lastgap_cut>0) continue;
      muon->set_DDG0(Tools::DDG0( muo, imu, iroad ));
      if (muon->get_DDG0() > DDG0_cut && lastgap_cut>0) continue;
      if (muon->get_p()*muon->get_DDG0() > pDDG0_cut && lastgap_cut>0) continue;
      muon->set_DS3(Tools::DS3( muo, imu, iroad ));
      muon->set_trchi2(muo->get_chisquare(imu));
      if (muon->get_trchi2() > chi2_cut) continue;
      muon->set_idchi2(muo->get_muIDOOchi(iroad,imu));
      if (muon->get_idchi2() > idchi2_cut  && lastgap_cut>0) continue;
      muon->set_ntrhits(Tools::sumbit( muon->get_trhits() ));
      if (muon->get_ntrhits() < ntrhits_cut) continue;
      muon->set_nidhits(Tools::sumbit( muon->get_idhits() ));
      if (muon->get_nidhits() < nidhits_cut  && lastgap_cut>0) continue;
      muon->set_charge(false);
      if (muo->get_charge(imu)>0) muon->set_charge(true);
      muon->set_x0(muo->get_xpos(0,imu));
      muon->set_xst1(muo->get_xpos(1,imu));
      muon->set_xst2(muo->get_xpos(2,imu));
      muon->set_xst3(muo->get_xpos(3,imu));
      muon->set_y0(muo->get_ypos(0,imu));
      muon->set_yst1(muo->get_ypos(1,imu));
      muon->set_yst2(muo->get_ypos(2,imu));
      muon->set_yst3(muo->get_ypos(3,imu));
      muon->set_z0(muo->get_zpos(0,imu));
      if (muon->get_lastgap()>0)
         {
         muon->set_idx(muo->get_muid_hit_x(muon->get_lastgap(),imu));
         muon->set_idy(muo->get_muid_hit_y(muon->get_lastgap(),imu));
         }
      
      TMuiGeometry* muigeom = TMuiGeometry::Geom();
      int arm = (muo->get_pz(0,imu)>0);
         // find direction for the road
      float dx = muo->get_muID_gap0(3,imu);
      float dy = muo->get_muID_gap0(4,imu);
      float dz = 1.0;
      if (arm==0)
         {
         dx *= -1.0;
         dy *= -1.0;
         dz *= -1.0;
         }
      PHVector vroad(dx, dy, dz);
      PHVector idpnt(muo->get_muid_hit_x(0, imu),
                     muo->get_muid_hit_y(0, imu),
                     muigeom->GapZPosition( arm, 0));
      vector<TMuiChannelId> muich_list = muigeom->findTwoPacks( arm, 0, idpnt, vroad);
      if (muich_list.size()>0)
         {
         TMuiChannelId muich = muich_list[0];
         muon->set_idpanel(muich.Panel());
         }
      
      float px_mut = muo->get_px(1,imu);
      float py_mut = muo->get_py(1,imu);
      float pz_mut = muo->get_pz(1,imu);
      muon->set_st1px(px_mut);
      muon->set_st1py(py_mut);
      muon->set_st1pz(pz_mut);
      
      if ( _use_MuIDLl1 )
         muon->set_MUID1D(trigger_LL1( imu, top_node ));
      
         // makes the association of each muon track with one of the vertices
      PHPoint vertex = PHPoint(0,0,evt->getBbcZVertex());
      if (vtx){
         vertex = associate_mut_vertex(muon);
         
            // Use run-averaged x,y position of vertex with slope correction if GlobalPar flag is set to do this:
         if (TFvtxGlobalParCntrl::get_bool_par("beam_use_average_xy")){
            vertex.setX(TFvtxGlobalParCntrl::get_float_par("beam_x_seed") + vertex.getZ() * TFvtxGlobalParCntrl::get_float_par("beam_dxdz"));
            vertex.setY(TFvtxGlobalParCntrl::get_float_par("beam_y_seed") + vertex.getZ() * TFvtxGlobalParCntrl::get_float_par("beam_dydz"));
         }
      }
      
      if (sm_version>4)
         muon->set_clusters_size1(muo->get_clusters_size1(imu));
      
      for (unsigned int i=0; i<5; i++)
         for (unsigned int j=0; j<5; j++)
            muon->set_cov(i, j, muo->get_cov(i,j,imu));
      
      if ( sm_version>11 )
         muon->set_track_id( muo->get_uid(imu) );
      
      if ( trk_map && sm_version>10 )
         {
         TMutTrkMap::const_iterator trk_iter = trk_map->range();
         unsigned long int uid = muo->get_uid(imu);
         while(TMutTrkMap::const_pointer trk_ptr = trk_iter.next())
            {
            unsigned long int icnt = trk_ptr->get()->get_key().get_obj_key();
            if (icnt == uid)
               {
               const TMutTrkPar* trk_par = trk_ptr->get()->get_trk_par_station(MUTOO::Station1);
               for (unsigned int i=0; i<5; i++)
                  for (unsigned int j=0; j<5; j++)
                     muon->set_cov_mutsta1(i,j, trk_par->get_covar(i,j));
               break;
               }
            }
         }
      
      int nhits_fvtx = Tools::sumbit(muo->get_fvtx_hits(imu));
         //      if (nhits_fvtx>=min_nhits_fvtx && sm_version>2 )// nhits not available in run13 production
      if (muo->get_fvtx_dr(imu)>-100 && sm_version>2)
         {
         muon->set_x0_fvtx(muo->get_fvtx_vtx(imu,0));
         muon->set_y0_fvtx(muo->get_fvtx_vtx(imu,1));
         muon->set_z0_fvtx(muo->get_fvtx_vtx(imu,2));
         
         muon->set_px_fvtx(muo->get_fvtx_p(imu,0));
         muon->set_py_fvtx(muo->get_fvtx_p(imu,1));
         muon->set_pz_fvtx(muo->get_fvtx_p(imu,2));
         
         muon->set_x0_fvtxmutr(muo->get_fvtxmutr_vtx(imu,0));
         muon->set_y0_fvtxmutr(muo->get_fvtxmutr_vtx(imu,1));
         muon->set_z0_fvtxmutr(muo->get_fvtxmutr_vtx(imu,2));
         
         muon->set_px_fvtxmutr(muo->get_fvtxmutr_p(imu,0));
         muon->set_py_fvtxmutr(muo->get_fvtxmutr_p(imu,1));
         muon->set_pz_fvtxmutr(muo->get_fvtxmutr_p(imu,2));
         
         muon->set_dphi_fvtx(muo->get_fvtx_dphi(imu));
         muon->set_dtheta_fvtx(muo->get_fvtx_dtheta(imu));
         muon->set_dr_fvtx(muo->get_fvtx_dr(imu));
         muon->set_chi2_fvtx(muo->get_fvtx_chi2(imu));
         muon->set_chi2_fvtxmutr(muo->get_fvtxmutr_chi2(imu));
         
         muon->set_hit_pattern(muo->get_fvtx_hits(imu));
         muon->set_nhits_fvtx(nhits_fvtx);
         
         if (Run_Number > BEGIN_OF_RUN12)
            for (unsigned int i=0; i<5; i++)
               for (unsigned int j=0; j<5; j++)
                  {
                  muon->set_cov_fvtx(i, j, muo->get_fvtx_cov(imu,i,j));
                  muon->set_cov_fvtxmutr(i, j, muo->get_fvtxmutr_cov(imu,i,j));
                  }
         
            // by 06-11-2013 the cluster size variable was buggy and the
            // number of hits from it was always return one.
            // For now the number of hits are taken from number of
            // stations with at least one FVTX coordinates.
            // From production made after 06-11-2013 it will take the
            // total number of FVTX coordinates.
         short nfvtxhits = 0;
         for (int i=0; i<4; i++)
            {
            muon->set_fvtx_strip(i, muo->get_fvtx_global_strip(imu,i));
            muon->set_fvtx_charge(i, muo->get_fvtx_cluster_charge(imu,i));
            if (muo->get_fvtx_cluster_charge(imu,i)>0) nfvtxhits++;
            }
         if (muon->get_nhits_fvtx()==0)
            {
            nfvtxhits = 0;
            short hit_pattern = 0;
            for (int i=0; i<8; i++)
               if (muo->get_fvtx_cluster_size(imu,i)>0)
                  {
                  nfvtxhits++;
                  hit_pattern |= (1 << i);
                  }
            if (nfvtxhits > (short)muon->get_nhits_fvtx())
               muon->set_nhits_fvtx(nfvtxhits);
            
            muon->set_hit_pattern( hit_pattern );
            }
         
         get_dca(vertex, muon);
         }
      if (muon->get_dr_fvtx()<0) // no FVTX matching
         {
         float dca_r, dca_x, dca_y, dca_z;
         Tools::DCA(muo, imu, dca_r, dca_x, dca_y, dca_z);
         muon->set_dca_z( dca_z );
         float pt = sqrt(MUTOO::SQUARE(muon->get_px()) + MUTOO::SQUARE(muon->get_py()));
         muon->set_dca_r( dca_r );
         float dca_phi = (dca_x*muon->get_px() - dca_y*muon->get_py()) / pt;
         muon->set_dca_phi( dca_phi );
         }
      
      if (match_fvtx && muon->get_dr_fvtx()<0) continue;
      
      if (fabs(muon->get_dca_z()) > dca_z_cut) continue;
      
      if (muon->get_pz() > 0) naccepted_N++;
      if (muon->get_pz() < 0) naccepted_S++;
      
      if (sm_version>4)
         muon->set_mutoo_trk_index(imu);
      
         // cone observables
         // require v15 or newer nDST, otherwise directly load from FvtxConeTracklets
      const unsigned int muo_version = MWGVersion::get(muo->ClassName());
      if (muo_version >= 15 && sm_version>2 )
         {
         
         static bool once = true;
         if (once)
            {
            once = false;
            cout <<"mFillSingleMuonContainer::process_event - INFO - "
            <<"Version "<<muo_version<<" of nDST is used in the TOP node. Fill cone observables as normal"
            <<endl;
            }
         
            // v15 or newer nDST
         muon->set_nfvtx_tracklets_cone(muo->get_fvtx_tracklets_cone(imu));
         muon->set_nfvtx_tracklets(muo->get_nfvtx_tracklets(imu));
         muon->set_nfvtx_clusters_cone(muo->get_fvtx_clusters_cone(imu));
         }
      else
         {
            // directly load from FvtxConeTracklets
            // this is intended as a fix for run12 510pp 1st production nDSTs
            // this should not be used in other productions/analysis
         
         static bool search_once = true;
         if (search_once)
            {
            search_once = false;
            
            Fun4AllServer *se = Fun4AllServer::instance();
            
            fvtx_cone = dynamic_cast<FvtxConeTracklets *>(se->getSubsysReco(
                                                                            "FVTXCONETRACKLETS"));
            
            if (fvtx_cone)
               {
               cout <<"mFillSingleMuonContainer::process_event - WARNING - "
               <<"Complete old version of nDST by directly importing cone observables from FvtxConeTracklets"<<endl;
               }
            else
               {
               cout <<"mFillSingleMuonContainer::process_event - WARNING - "
               <<"Cone observables is not available. "
               <<"Please use either v15 or newer nDST or load FvtxConeTracklets in your script"<<endl;
               }
            }
         if (fvtx_cone)
            {
            const PHMuoTracksOut * muo_cone = fvtx_cone->get_muo_local();
            
            if (!muo_cone)
               {
               cout <<"mFillSingleMuonContainer::process_event - Error - "
               <<"cannot read local PHMuoTracksOut from FvtxConeTracklets"<<endl;
               }
            else
               {
               muon->set_nfvtx_tracklets_cone(muo_cone->get_fvtx_tracklets_cone(imu));
               muon->set_nfvtx_tracklets(muo_cone->get_nfvtx_tracklets(imu));
               muon->set_nfvtx_clusters_cone(muo_cone->get_fvtx_clusters_cone(imu));
               }
            }
         }
      
         // Fill Run11 Rpc info
      if (sm_version==5 || sm_version==6)
         {
         Float_t fTracks_Rpcpx  = -9999.;
         Float_t fTracks_Rpcpy  = -9999.;
         Float_t fTracks_Rpcpz  = -9999.;
         Float_t fTracks_RpcDca = -9999. ;
         Float_t fTracks_RpcPt   = -9999.  ;
         
         Float_t fTracks_RpcDcaVt_1  = -9999.;
         Float_t fTracks_Rpctime_1   = -9999.;
         
         if (rpctrk_map)
            {
               //   cout << "after hitmap2 " << endl;
            
            TRpcTrkMap::iterator trk_iter = rpctrk_map->range();
               //cout << "rpctrack hitmap range " << rpctrk_map->size()<< endl;
            while(TRpcTrkMap::pointer trk_ptr = trk_iter.next())
               {
               TMutTrkPar *fTrkPar = (TMutTrkPar *) trk_ptr->get()->get_trk_par_vtx();
                  //TMutTrkPar *fTrkPar = (TMutTrkPar *) trk_ptr->get()->get_trk_par();
                  //int charge_r = fTrkPar->get_charge();
               
               Float_t p2 = fTrkPar->get_px()*fTrkPar->get_px();
               p2 += fTrkPar->get_py()*fTrkPar->get_py();
               p2 += fTrkPar->get_pz()*fTrkPar->get_pz();
               
               fTracks_Rpcpx = fTrkPar->get_px();
               fTracks_Rpcpy = fTrkPar->get_py();
               fTracks_Rpcpz = fTrkPar->get_pz();
               
               
               if ( fabs( muo->get_px(0,imu) -  fTracks_Rpcpx ) < 0.00001 &&
                   fabs( muo->get_py(0,imu) -  fTracks_Rpcpy ) < 0.00001 &&
                   fabs( muo->get_pz(0,imu) -  fTracks_Rpcpz ) < 0.00001)
                  {
                  fTracks_RpcDca   = trk_ptr->get()->get_dca_trk();
                  fTracks_RpcPt    = sqrt((pow(fTracks_Rpcpx,2)+pow(fTracks_Rpcpy,2))/p2); //trk_ptr->get()->get_corr_pT();
                  
                  fTracks_RpcDcaVt_1 = trk_ptr->get()->get_dca_trk_vtx1();
                  fTracks_Rpctime_1  = trk_ptr->get()->get_rpcclus1time_vtx();
                  
                  muon->set_RpcDCA(fTracks_RpcDca);
                  muon->set_RpcpT(fTracks_RpcPt);
                  muon->set_Rpctime((float)trk_ptr->get()->get_rpcclus3time());//fTracks_Rpctime;
                  muon->set_Rpc1DCA(fTracks_RpcDcaVt_1);
                  muon->set_Rpc1time(fTracks_Rpctime_1);
                  }
               }
            } // rpctrk_map
         
         
            // Fill Run12 Rpc info
         if(rpc_muotrk_map)
            {
            TRpcMuoTrkMap::iterator trk_iter = rpc_muotrk_map->range();
               //            cout << "rpcmuotrack map range " << rpc_muotrk_map->size()<< endl;
            while(TRpcMuoTrkMap::pointer trk_ptr = trk_iter.next())
               {
               if ( fabs(muo->get_px(0,imu) - trk_ptr->get()->get_muo_trk_momentum(0) ) < 0.00001 &&
                   fabs(muo->get_py(0,imu) - trk_ptr->get()->get_muo_trk_momentum(1) ) < 0.00001 &&
                   fabs(muo->get_pz(0,imu) - trk_ptr->get()->get_muo_trk_momentum(2) ) < 0.00001)
                  {
                        //               cout << "found matched track " << trk_ptr->get()->get_muo_trk_number() << " " << ipart << endl;
                     muon->set_Rpc1VtxDCA(trk_ptr->get()->get_dca_vtx(0));
                     muon->set_Rpc1VtxTime(trk_ptr->get()->get_rpcclus_time_vtx(0));
                     muon->set_Rpc3VtxDCA (trk_ptr->get()->get_dca_vtx(1));
                     muon->set_Rpc3VtxTime(trk_ptr->get()->get_rpcclus_time_vtx(1));
                     
                     muon->set_Rpc1St3DCA (trk_ptr->get()->get_dca_trk(0,1));
                     muon->set_Rpc1St3Time(trk_ptr->get()->get_rpcclus_time(0,1));
                     muon->set_Rpc3St3DCA (trk_ptr->get()->get_dca_trk(1,1));
                     muon->set_Rpc3St3Time(trk_ptr->get()->get_rpcclus_time(1,1));
                     
                     muon->set_Rpc1St1DCA (trk_ptr->get()->get_dca_trk(0,0));
                     muon->set_Rpc1St1Time(trk_ptr->get()->get_rpcclus_time(0,0));
                     muon->set_Rpc3St1DCA (trk_ptr->get()->get_dca_trk(1,0));
                     muon->set_Rpc3St1Time(trk_ptr->get()->get_rpcclus_time(1,0));
                     
                     muon->set_Rpc1MuIDDCA (trk_ptr->get()->get_dca_muid(0));
                     muon->set_Rpc1MuIDTime(trk_ptr->get()->get_rpcclus_time_muid(0));
                     muon->set_Rpc3MuIDDCA (trk_ptr->get()->get_dca_muid(1));
                     muon->set_Rpc3MuIDTime(trk_ptr->get()->get_rpcclus_time_muid(1));
                     
                  }//matched
               }//trk_ptr
            }//rpc_trk_map
         }
      
      muons->AddSingleMuon(*muon);
      }
   delete muon;
   
   if (muons->get_nSingleMuons() > 0) return 0;
   return DISCARDEVENT;
   
}  // -- end of event processing

   //______________________________________________________
int mFillSingleMuonContainer::End(PHCompositeNode *top_node)
{
   
   MUTOO::PRINT( cout, "mFillSingleMuonContainer::End" );
   
   cout << "events: " << nevents << endl;
   cout << "North Arm accepted muons:  " << naccepted_N << endl;
   cout << "South Arm accepted muons:  " << naccepted_S << endl;
   
   return 0 ;
}

   // --- define user functions for single muons
   //__________________________________________________________
int mFillSingleMuonContainer::get_best_road_oo( int imu, PHMuoTracksOut* muo)
{
   
      // keep track of the "best" accepted road
   int best_accepted_road( -1 );
   double best_accepted_dg0( -1 );
   
   for( int i_road=0; i_road<3; i_road++ )
      {
      
         // check if road is present
      if( !muo->get_muIDOOhits( i_road, imu ) ) continue;
      
         // retrieve track DG0
      double dg0( Tools::DG0( muo, imu, i_road ) );
      
      /*
       check if road DG0 is smaller than others
       it is done for debugging only. The 'true' check
       must be done only for accepted roads
       */
      if( best_accepted_dg0<0 || dg0 < best_accepted_dg0 )
         {
         best_accepted_road = i_road;
         best_accepted_dg0 = dg0;
         }
      }
   
      // return the accepted road which have the smaller DG0
   return best_accepted_road;
   
}

void
mFillSingleMuonContainer::get_dca(PHPoint vertex, TMutTrkPar* fvtx_par,
                                  float& dca_r, float& dca_phi, float& dca_z)
{
   PHTrackIntegratorKF integrator;
   integrator.clear();
   TMutTrkPar extrap_trk_par;
   integrator.initialize(*fvtx_par);
   
   float zref = (fvtx_par->get_pz() > 0) ? 25.0 : -25.0;
   
   integrator.extrapolate ( zref );
   if (integrator.get_error()){
      extrap_trk_par.set_x(fvtx_par->get_x());
      extrap_trk_par.set_y(fvtx_par->get_y());
      extrap_trk_par.set_z(fvtx_par->get_z());
   }
   else{
      integrator.finish( extrap_trk_par );
   }
   
   PHVector vtrack(extrap_trk_par.get_x() - vertex.getX(),
                   extrap_trk_par.get_y() - vertex.getY(),
                   extrap_trk_par.get_z() - vertex.getZ());
   PHLine track(fvtx_par->get_point(), vtrack);
   
   zref = vertex.getZ();
   integrator.initialize(*fvtx_par);
   integrator.extrapolate ( zref );
   if (integrator.get_error()){
      extrap_trk_par.set_x(fvtx_par->get_x());
      extrap_trk_par.set_y(fvtx_par->get_y());
      extrap_trk_par.set_z(fvtx_par->get_z());
   }
   else{
      integrator.finish( extrap_trk_par );
   }
   float dx = extrap_trk_par.get_x() - vertex.getX();
   float dy = extrap_trk_par.get_y() - vertex.getY();
   float R = sqrt(vtrack.getX()*vtrack.getX()+vtrack.getY()*vtrack.getY());
   dca_r = (dx*vtrack.getX() + dy*vtrack.getY())/R;
   dca_phi = (dx*vtrack.getY() - dy*vtrack.getX())/R;
   
   PHPoint CA = PHGeometry::closestApproachLinePoint(track, vertex);
   vtrack.setZ(0);
   dx = CA.getX()-vertex.getX();
   dy = CA.getY()-vertex.getY();
   dca_z = (dx*vtrack.getX() + dy*vtrack.getY())/R;
}

void
mFillSingleMuonContainer::get_dca(PHPoint vertex, SingleMuon* muon, bool swapped)
{
   if (swapped)
      {
      cout << PHWHERE << "this function cannot be called for swapped fvtx tracks" << endl;
      return;
      }
   
      // Extrapolate FVTX track to ~first FVTX station and use this location
      // to decompose into dca_r, dca_phi components rather than the momentum
      // vector at the vertex:
   TMutTrkPar fvtx_par(muon->get_x0_fvtxmutr(),
                       muon->get_y0_fvtxmutr(),
                       muon->get_z0_fvtxmutr(),
                       muon->get_px_fvtxmutr(),
                       muon->get_py_fvtxmutr(),
                       muon->get_pz_fvtxmutr(),
                       static_cast<int>((muon->get_charge()>0)? 1 : -1),
                       muon->get_chi2_fvtxmutr());
   
   for (int i=0; i<TMutTrkPar::COVAR_ROW; i++)
      for (int j=0; j<TMutTrkPar::COVAR_ROW; j++)
         fvtx_par.set_covar(i,j, muon->get_cov_fvtxmutr(i,j));
   
   float dca_r=-999., dca_phi=-999., dca_z=-999.;
   get_dca(vertex, &fvtx_par, dca_r, dca_phi, dca_z);
   
   muon->set_dca_z(dca_z);
   muon->set_dca_r(dca_r);
   muon->set_dca_phi(dca_phi);
}

void
mFillSingleMuonContainer::get_dca(PHPoint vertex, TFvtxCompactTrk* fvtx_track,
                                  float& dca_r, float& dca_phi, float& dca_z)
{
   PHPoint p0 = fvtx_track->get_fvtx_vtx();
   
   float phi = fvtx_track->get_fvtx_phi();
   float theta = fvtx_track->get_fvtx_theta();
   
   float px = cos(phi)*sin(theta);
   float py = sin(phi)*sin(theta);
   float pz = cos(theta);
   
   PHVector vtrack = PHVector(px, py, pz);
   vtrack.normalize();
   
   PHLine track(p0, vtrack);
   
   get_dca(vertex, track, dca_r, dca_phi, dca_z);
}

void
mFillSingleMuonContainer::get_dca(PHPoint vertex, PHLine track, float& dca_r, float& dca_phi, float& dca_z)
{
   float zdir = (track.getDirection().getZ()>0) ? 1.0 : -1.0;
   
   PHVector vtrack = track.getDirection();
   vtrack.setX(vtrack.getX());
   vtrack.setY(vtrack.getY());
   
   PHVector vzref = PHVector(0, 0, zdir);
   vzref.normalize();
   PHLine plzref_norm(vertex, vzref);
   PHPlane plzref = PHPlane(plzref_norm);
   
   PHPoint pdca = PHPoint();
   PHGeometry::intersectionLinePlane(track,plzref,pdca);
   
   float dx0 = pdca.getX() - vertex.getX();
   float dy0 = pdca.getY() - vertex.getY();
   
   float px = vtrack.getX();
   float py = vtrack.getY();
   float pt = sqrt(px*px + py*py);
   
      // dca on the pT direction
   dca_r = (dx0*px + dy0*py)/pt;
      // dca perpendicular to pT direction
   dca_phi = (dx0*py - dy0*px)/pt;
   
   PHLine beam_line = PHLine(vertex, vzref);
   PHPoint CA = PHGeometry::closestApproachLineLine(track, beam_line);
   dca_z = CA.getZ() - vertex.getZ();
}

PHPoint
mFillSingleMuonContainer::associate_mut_vertex(SingleMuon* muon)
{
      // associate track with a vertex
      // First vertex choice is to use the FVTX (which includes VTX). If not available,
      // try SVX_PRECISE, else use BBC. If more than one FVTX vertex available, choose
      // the closest to this track projection:
   
   PHPoint vertex;
   int best_vertex = -1;
   if (is_sim){
      vertex = vtx->get_Vertex("SIM");
      if (sm_version >= 8) muon->set_vtx_index(mFillSingleMuonContainer::SIM);
      return vertex;
   }
      // from now on only real data
   
   string vertex_name[7] = {_fvtx_vertex_names[0].data(),
      _fvtx_vertex_names[1].data(),
      _fvtx_vertex_names[2].data(),
      _fvtx_vertex_names[3].data(),
      "SVX_PRECISE", "BBC"};
   float best_dca = track_fvtxvertex_cut;
      // check if the best vertex match
   for (int ivertex=0; ivertex<4; ivertex++)
      {
      vertex = vtx->get_Vertex(vertex_name[ivertex].c_str());
      if ( !(fabs(vertex.getZ()) < 200) ) continue;
      PHPoint p0(muon->get_x0(), muon->get_y0(), muon->get_z0());
      PHVector vtrack(muon->get_px(), muon->get_py(), muon->get_pz());
      PHLine track(p0, vtrack);
      PHPoint CA = PHGeometry::closestApproachLinePoint(track, vertex);
      float dca_vertex = PHGeometry::distancePointToPoint(CA, vertex);
      if (dca_vertex < best_dca)
         {
         best_dca = dca_vertex;
         best_vertex = ivertex;
         }
      }
   if (best_vertex != -1 && sm_version >= 8)
      {
      muon->set_vtx_index(best_vertex);
      return vtx->get_Vertex(vertex_name[best_vertex].c_str());
      }
   
      // if there is no good FVTX vertex, try SVX vertex
   vertex = vtx->get_Vertex("SVX_PRECISE");
   if ( fabs(vertex.getZ()) < 200 )
      {
      PHPoint p0(muon->get_x0(), muon->get_y0(), muon->get_z0());
      PHVector vtrack(muon->get_px(), muon->get_py(), muon->get_pz());
      PHLine track(p0, vtrack);
      PHPoint CA = PHGeometry::closestApproachLinePoint(track, vertex);
      float dca_vertex = PHGeometry::distancePointToPoint(CA, vertex);
      if ( dca_vertex < best_dca && sm_version>=8 )
         {
         muon->set_vtx_index(mFillSingleMuonContainer::SVX_PRECISE);
         return vertex;
         }
      }
   
      // if vertex detectors failed to provide a vertex, use BBC
   vertex = vtx->get_Vertex("BBC");
   if ( fabs(vertex.getZ()) < 200 )
      {
      PHPoint p0(muon->get_x0(), muon->get_y0(), muon->get_z0());
      PHVector vtrack(muon->get_px(), muon->get_py(), muon->get_pz());
      PHLine track(p0, vtrack);
      PHPoint CA = PHGeometry::closestApproachLinePoint(track, vertex);
      float dca_vertex = PHGeometry::distancePointToPoint(CA, vertex);
      if ( dca_vertex < best_dca && sm_version>=8 )
         {
         muon->set_vtx_index(mFillSingleMuonContainer::BBC);
         return vertex;
         }
      }
      // at last, set a NO vertex association and return whatever vertex available
   if ( sm_version>=8 )
      muon->set_vtx_index(-1);
   return vertex;
}

bool
mFillSingleMuonContainer::trigger_LL1(int imu, PHCompositeNode *top_node)
{
      // this function verify if this track alone can fire the LL1 trigger
   PHMuoTracksOut* muo = findNode::getClass<PHMuoTracksOut>(top_node,"PHMuoTracksOO");
   TMuiHitMapO* muihit_map = TMutNode<TMuiHitMapO>::find_node( top_node, "TMuiHitMapO_track" );
   if (!muihit_map)
      {
      cout << "mSingleMuonContainer::trigger_LL1:: TMuiHitMapO_track not found" << endl;
      }
   TMuiGeometry* muigeom = TMuiGeometry::Geom();
   
   int arm = (muo->get_pz(0,imu)>0);
      // find direction for the road
   float dx = muo->get_muID_gap0(3,imu);
   float dy = muo->get_muID_gap0(4,imu);
   float dz = 1.0;
   if (arm==0)
      {
      dx *= -1.0;
      dy *= -1.0;
      dz *= -1.0;
      }
   PHVector vroad(dx, dy, dz);
   
      // collect road points in each gap and fill TMuiHitMap
   muihit_map->clear();
   for (int igap=0; igap<5; igap++)
      {
      
      PHVector pnt(muo->get_muid_hit_x(igap, imu),
                   muo->get_muid_hit_y(igap, imu),
                   muigeom->GapZPosition( arm, igap));
      vector<TMuiChannelId> muich_list = muigeom->findTwoPacks( arm, igap, pnt, vroad);
      for ( size_t ich=0; ich<muich_list.size(); ich++ )
         {
         TMuiChannelId muich = muich_list[ich];
         muihit_map->insert_new(arm, igap, muich.Panel(), muich.Orient(), muich.TwoPack());
         }
      }
      // run trigger emulator with the hits in the container
   simMuIDLl1->getDataFromMutoo( top_node );
   if (arm==0) return simMuIDLl1->GL1_1Deep_S();
   else return simMuIDLl1->GL1_1Deep_N();
}

/*
 int
 mFillSingleMuonContainer::trigger_2D_LL1(int imu, PHCompositeNode *top_node)
 {
 // this function verify if this track alone can fire the LL1 trigger
 PHMuoTracksOut* muo = findNode::getClass<PHMuoTracksOut>(top_node,"PHMuoTracksOO");
 TMuiHitMapO* muihit_map = TMutNode<TMuiHitMapO>::find_node( top_node, "TMuiHitMapO_track" );
 if (!muihit_map)
 {
 cout << "mSingleMuonContainer::trigger_LL1:: TMuiHitMapO_track not found" << endl;
 }
 TMuiGeometry* muigeom = TMuiGeometry::Geom();
 
 int npart = muo->get_npart();
 for ( int imu2; imu2<npart; imu2++)
 {
 int arm = (muo->get_pz(0,imu)>0);
 // find direction for the road
 float dx = muo->get_muID_gap0(3,imu);
 float dy = muo->get_muID_gap0(4,imu);
 float dz = 1.0;
 if (arm==0)
 {
 dx *= -1.0;
 dy *= -1.0;
 dz *= -1.0;
 }
 PHVector vroad(dx, dy, dz);
 
 // collect road points in each gap and fill TMuiHitMap
 //      muihit_map->clear();
 for (int igap=0; igap<5; igap++)
 {
 
 PHVector pnt(muo->get_muid_hit_x(igap, imu),
 muo->get_muid_hit_y(igap, imu),
 muigeom->GapZPosition( arm, igap));
 vector<TMuiChannelId> muich_list = muigeom->findTwoPacks( arm, igap, pnt, vroad);
 for ( size_t ich=0; ich<muich_list.size(); ich++ )
 {
 TMuiChannelId muich = muich_list[ich];
 muihit_map->insert_new(arm, igap, muich.Panel(), muich.Orient(), muich.TwoPack());
 }
 }
 
 // run trigger emulator with the hits in the container
 simMuIDLl1->getDataFromMutoo( top_node );
 bool fired = false;
 if (arm==0) fired = simMuIDLl1->GL1_2Deep_S();
 else fired = simMuIDLl1->GL1_2Deep_N();
 if ( fired ) return imu2;
 }
 }
 */
