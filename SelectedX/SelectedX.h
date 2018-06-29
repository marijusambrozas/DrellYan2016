#pragma once

#define MaxN 50000
#include <TTree.h>
#include <TTreeCache.h>
#include <TChain.h>
#include <vector>

class SelectedX_t
{
public:
	TChain *chain;

    //Event Informations
    Int_t nVertices;
    Int_t runNum;
    Int_t lumiBlock;
    Int_t evtNum;
    Int_t nPileUp;

    //Trigger variables
    Int_t HLT_ntrig;
    Int_t HLT_trigFired[MaxN];
    vector<string> *HLT_trigName;
    Double_t HLT_trigPt[MaxN];
    Double_t HLT_trigEta[MaxN];
    Double_t HLT_trigPhi[MaxN];

    Double_t GENEvt_weight;

    Double_t Inv_Mass[MaxN];

    //////////////////////////////
    // -- Electron Variables -- //
    //////////////////////////////
    Int_t Nelectrons;
    Double_t Electron_Energy[MaxN];
    Double_t Electron_pT[MaxN];
    Double_t Electron_eta[MaxN];
    Double_t Electron_phi[MaxN];
    Int_t Electron_charge[MaxN];
//    Double_t Electron_gsfpT[MaxN];
//    Double_t Electron_gsfPx[MaxN];
//    Double_t Electron_gsfPy[MaxN];
//    Double_t Electron_gsfPz[MaxN];
//    Double_t Electron_gsfEta[MaxN];
//    Double_t Electron_gsfPhi[MaxN];
//    Double_t Electron_gsfCharge[MaxN];
    Double_t Electron_etaSC[MaxN];
    Double_t Electron_phiSC[MaxN];
//    Double_t Electron_etaWidth[MaxN];
//    Double_t Electron_phiWidth[MaxN];
//    Double_t Electron_dEtaIn[MaxN];
//    Double_t Electron_dPhiIn[MaxN];
//    Double_t Electron_sigmaIEtaIEta[MaxN];
//    Double_t Electron_HoverE[MaxN];
//    Double_t Electron_fbrem[MaxN];
//    Double_t Electron_eOverP[MaxN];
//    Double_t Electron_InvEminusInvP[MaxN];
//    Double_t Electron_dxyVTX[MaxN];
//    Double_t Electron_dzVTX[MaxN];
//    Double_t Electron_dxy[MaxN];
//    Double_t Electron_dz[MaxN];
//    Double_t Electron_dxyBS[MaxN];
//    Double_t Electron_dzBS[MaxN];
//    Double_t Electron_chIso03[MaxN];
//    Double_t Electron_nhIso03[MaxN];
//    Double_t Electron_phIso03[MaxN];
//    Double_t Electron_ChIso03FromPU[MaxN];
//    Int_t Electron_mHits[MaxN];
//    Double_t Electron_EnergySC[MaxN];
//    Double_t Electron_preEnergySC[MaxN];
//    Double_t Electron_rawEnergySC[MaxN];
//    Double_t Electron_etSC[MaxN];
//    Double_t Electron_E15[MaxN];
//    Double_t Electron_E25[MaxN];
//    Double_t Electron_E55[MaxN];
//    Double_t Electron_RelPFIso_dBeta[MaxN];
//    Double_t Electron_RelPFIso_Rho[MaxN];
//    Double_t Electron_r9[MaxN];
//    Int_t Electron_ecalDriven[MaxN];
//    Int_t Electron_passConvVeto[MaxN];

//    Bool_t Electron_passLooseID[MaxN];
//    Bool_t Electron_passMediumID[MaxN];
//    Bool_t Electron_passTightID[MaxN];
//    Bool_t Electron_passMVAID_WP80[MaxN];
//    Bool_t Electron_passMVAID_WP90[MaxN];
//    Bool_t Electron_passHEEPID[MaxN];

    //////////////////////////
    // -- Muon Variables -- //
    //////////////////////////
    // -- Physical Variables -- //
    Double_t Muon_pT[MaxN];
    Double_t Muon_eta[MaxN];
    Double_t Muon_phi[MaxN];
    
    // -- Cut variables -- //
//    Int_t Muon_muonType[MaxN];
//    Double_t Muon_chi2dof[MaxN];
//    Int_t Muon_muonHits[MaxN];
//    Int_t Muon_nSegments[MaxN];
//    Int_t Muon_nMatches[MaxN];
//    Int_t Muon_trackerLayers[MaxN];

//     Int_t Muon_trackerHitsGLB[MaxN];
//    Int_t Muon_pixelHitsGLB[MaxN];
//    Int_t Muon_trackerLayersGLB[MaxN];
    
//    Int_t Muon_pixelHits[MaxN];
//    Double_t Muon_dxyVTX[MaxN];
//    Double_t Muon_dzVTX[MaxN];
//    Double_t Muon_trkiso[MaxN];
//    Int_t isGLBmuon[MaxN];
//    Int_t isPFmuon[MaxN];
//    Int_t isTRKmuon[MaxN];
    Int_t nMuon;
    
    // -- for invariant mass calculation -- //
    Double_t Muon_Px[MaxN];
    Double_t Muon_Py[MaxN];
    Double_t Muon_Pz[MaxN];
    
//    Double_t Muon_dB[MaxN];
    // -- for the muon momentum corrections -- //
    Int_t Muon_charge[MaxN];
    
    //PF information
//    Double_t Muon_PfChargedHadronIsoR04[MaxN];
//    Double_t Muon_PfNeutralHadronIsoR04[MaxN];
//    Double_t Muon_PfGammaIsoR04[MaxN];
//    Double_t Muon_PFSumPUIsoR04[MaxN];

    //Dimuon variables
//    std::vector<double> *CosAngle;
//    std::vector<double> *vtxTrkChi2;
//    std::vector<double> *vtxTrkProb;
//    std::vector<double> *vtxTrkNdof;
//    std::vector<double> *vtxTrkCkt1Pt;
//    std::vector<double> *vtxTrkCkt2Pt;
//    std::vector<double> *vtxTrkDiEChi2;
//    std::vector<double> *vtxTrkDiEProb;
//    std::vector<double> *vtxTrkDiENdof;
//    std::vector<double> *vtxTrkDiE1Pt;
//    std::vector<double> *vtxTrkDiE2Pt;
//    std::vector<double> *vtxTrkEMuChi2;
//    std::vector<double> *vtxTrkEMuProb;
//    std::vector<double> *vtxTrkEMuNdof;
//    std::vector<double> *vtxTrkEMu1Pt;
//    std::vector<double> *vtxTrkEMu2Pt;

//    std::vector<double> *CosAngle_TuneP;
//    std::vector<double> *vtxTrk1Pt_TuneP;
//    std::vector<double> *vtxTrk2Pt_TuneP;
//    std::vector<double> *vtxTrkChi2_TuneP;
//    std::vector<double> *vtxTrkNdof_TuneP;
//    std::vector<double> *vtxTrkProb_TuneP;

    // -- Various Track Information -- //
//    Double_t Muon_Best_pT[MaxN];
//    Double_t Muon_Best_pTError[MaxN];
//    Double_t Muon_Best_Px[MaxN];
//    Double_t Muon_Best_Py[MaxN];
//    Double_t Muon_Best_Pz[MaxN];
//    Double_t Muon_Best_eta[MaxN];
//    Double_t Muon_Best_phi[MaxN];

//    Double_t Muon_Inner_pT[MaxN];
//    Double_t Muon_Inner_pTError[MaxN];
//    Double_t Muon_Inner_Px[MaxN];
//    Double_t Muon_Inner_Py[MaxN];
//    Double_t Muon_Inner_Pz[MaxN];
//    Double_t Muon_Inner_eta[MaxN];
//    Double_t Muon_Inner_phi[MaxN];

//    Double_t Muon_Outer_pT[MaxN];
//    Double_t Muon_Outer_pTError[MaxN];
//    Double_t Muon_Outer_Px[MaxN];
//    Double_t Muon_Outer_Py[MaxN];
//    Double_t Muon_Outer_Pz[MaxN];
//    Double_t Muon_Outer_eta[MaxN];
//    Double_t Muon_Outer_phi[MaxN];

//    Double_t Muon_GLB_pT[MaxN];
//    Double_t Muon_GLB_pTError[MaxN];
//    Double_t Muon_GLB_Px[MaxN];
//    Double_t Muon_GLB_Py[MaxN];
//    Double_t Muon_GLB_Pz[MaxN];
//    Double_t Muon_GLB_eta[MaxN];
//    Double_t Muon_GLB_phi[MaxN];

//    Double_t Muon_TuneP_pT[MaxN];
//    Double_t Muon_TuneP_pTError[MaxN];
//    Double_t Muon_TuneP_Px[MaxN];
//    Double_t Muon_TuneP_Py[MaxN];
//    Double_t Muon_TuneP_Pz[MaxN];
//    Double_t Muon_TuneP_eta[MaxN];
//    Double_t Muon_TuneP_phi[MaxN];

//    Int_t Muon_stationMask[MaxN];
//    Int_t Muon_nMatchesRPCLayers[MaxN];


    // -- Constructor -- //
    NtupleHandle(TChain *chainptr)
    {
    	chain = chainptr;
        Int_t cachesize = 500000000; // -- 500MB -- //
        chain->SetCacheSize(cachesize);
        chain->SetBranchStatus("*", 0);

    	// -- Event Information -- //
    	chain->SetBranchStatus("nVertices", 1);
    	chain->SetBranchStatus("runNum", 1);
        chain->SetBranchStatus("lumiBlock", 1);
        chain->SetBranchStatus("evtNum", 1);
    	chain->SetBranchStatus("nPileUp", 1);

    	chain->SetBranchAddress("nVertices", &nVertices);
    	chain->SetBranchAddress("runNum", &runNum);
        chain->SetBranchAddress("lumiBlock", &lumiBlock);
        chain->SetBranchAddress("evtNum", &evtNum);
    	chain->SetBranchAddress("nPileUp", &nPileUp);

        chain->AddBranchToCache("nVertices", 1);
        chain->AddBranchToCache("runNum", 1);
        chain->AddBranchToCache("lumiBlock", 1);
        chain->AddBranchToCache("evtNum", 1);
        chain->AddBranchToCache("nPileUp", 1);

    	// -- Trigger Information -- //
    	chain->SetBranchStatus("HLT_trigName", 1);
    	chain->SetBranchStatus("HLT_ntrig", 1);
    	chain->SetBranchStatus("HLT_trigFired", 1);
        chain->SetBranchStatus("HLT_trigPt", 1);
    	chain->SetBranchStatus("HLT_trigEta", 1);
    	chain->SetBranchStatus("HLT_trigPhi", 1);

    	chain->SetBranchAddress("HLT_trigName", &HLT_trigName);
    	chain->SetBranchAddress("HLT_trigFired", HLT_trigFired);
    	chain->SetBranchAddress("HLT_ntrig", &HLT_ntrig);
        chain->SetBranchAddress("HLT_trigPt", &HLT_trigPt);
    	chain->SetBranchAddress("HLT_trigEta", &HLT_trigEta);
    	chain->SetBranchAddress("HLT_trigPhi", &HLT_trigPhi);

        chain->AddBranchToCache("HLT_trigName", 1);
        chain->AddBranchToCache("HLT_ntrig", 1);
        chain->AddBranchToCache("HLT_trigFired", 1);
        chain->AddBranchToCache("HLT_trigPt", 1);
        chain->AddBranchToCache("HLT_trigEta", 1);
        chain->AddBranchToCache("HLT_trigPhi", 1);


    	// this->TurnOnBranches_GenLepton();
    	// this->TurnOnBranches_Muon();

        chain->SetBranchStatus("GENEvt_weight", 1);
        chain->SetBranchStatus("Inv_Mass", 1);

        chain->SetBranchAddress("GENEvt_weight", &GENEvt_weight);
        chain->SetBranchAddress("Inv_Mass", &Inv_Mass);

        chain->AddBranchToCache("GENEvt_weight", 1);
        chain->AddBranchToCache("Inv_Mass", 1);

    }

    void Ready()
    {
        chain->StopCacheLearningPhase();
    }

    void TurnOnBranches_Muon()
    {
//    	chain->SetBranchStatus("isPFmuon", 1);
//    	chain->SetBranchStatus("isGLBmuon", 1);
//        chain->SetBranchStatus("isTRKmuon", 1);
//    	chain->SetBranchStatus("CosAngle", 1);
//    	chain->SetBranchStatus("vtxTrkChi2", 1);
//    	chain->SetBranchStatus("vtxTrkProb", 1);
//    	chain->SetBranchStatus("vtxTrkNdof", 1);
//    	chain->SetBranchStatus("vtxTrkCkt1Pt", 1);
//    	chain->SetBranchStatus("vtxTrkCkt2Pt", 1);
//    	chain->SetBranchStatus("vtxTrkDiEChi2", 1);
//    	chain->SetBranchStatus("vtxTrkDiEProb", 1);
//    	chain->SetBranchStatus("vtxTrkDiENdof", 1);
//    	chain->SetBranchStatus("vtxTrkDiE1Pt", 1);
//    	chain->SetBranchStatus("vtxTrkDiE2Pt", 1);
//    	chain->SetBranchStatus("vtxTrkEMuChi2", 1);
//    	chain->SetBranchStatus("vtxTrkEMuProb", 1);
//    	chain->SetBranchStatus("vtxTrkEMuNdof", 1);
//    	chain->SetBranchStatus("vtxTrkEMu1Pt", 1);
//    	chain->SetBranchStatus("vtxTrkEMu2Pt", 1);

//        chain->SetBranchStatus("CosAngle_TuneP", 1);
//        chain->SetBranchStatus("vtxTrk1Pt_TuneP", 1);
//        chain->SetBranchStatus("vtxTrk2Pt_TuneP", 1);
//        chain->SetBranchStatus("vtxTrkChi2_TuneP", 1);
//        chain->SetBranchStatus("vtxTrkNdof_TuneP", 1);
//        chain->SetBranchStatus("vtxTrkProb_TuneP", 1);

        chain->SetBranchStatus("nMuon", 1);
        chain->SetBranchStatus("Muon_pT", 1);
        chain->SetBranchStatus("Muon_eta", 1);
        chain->SetBranchStatus("Muon_phi", 1);
//        chain->SetBranchStatus("Muon_muonType", 1);
//        chain->SetBranchStatus("Muon_chi2dof", 1);
//        chain->SetBranchStatus("Muon_muonHits", 1);
//        chain->SetBranchStatus("Muon_nSegments", 1);
//        chain->SetBranchStatus("Muon_nMatches", 1);
//        chain->SetBranchStatus("Muon_trackerLayers", 1);

        // chain->SetBranchStatus("Muon_trackerHitsGLB", 1);
//        chain->SetBranchStatus("Muon_pixelHitsGLB", 1);
//        chain->SetBranchStatus("Muon_trackerLayersGLB", 1);

//        chain->SetBranchStatus("Muon_pixelHits", 1);
//        chain->SetBranchStatus("Muon_dxyVTX", 1);
//        chain->SetBranchStatus("Muon_dzVTX", 1);
//        chain->SetBranchStatus("Muon_trkiso", 1);
        
        chain->SetBranchStatus("Muon_Px", 1);
        chain->SetBranchStatus("Muon_Py", 1);
        chain->SetBranchStatus("Muon_Pz", 1);

//        chain->SetBranchStatus("Muon_dB", 1);
        
        chain->SetBranchStatus("Muon_charge", 1);
        
//        chain->SetBranchStatus("Muon_PfChargedHadronIsoR04", 1);
//        chain->SetBranchStatus("Muon_PfNeutralHadronIsoR04" ,1);
//        chain->SetBranchStatus("Muon_PfGammaIsoR04", 1);
//        chain->SetBranchStatus("Muon_PFSumPUIsoR04", 1);

//        chain->SetBranchStatus("Muon_Best_pT", 1);
//        chain->SetBranchStatus("Muon_Best_pTError", 1);
//        chain->SetBranchStatus("Muon_Best_Px", 1);
//        chain->SetBranchStatus("Muon_Best_Py", 1);
//        chain->SetBranchStatus("Muon_Best_Pz", 1);
//        chain->SetBranchStatus("Muon_Best_eta", 1);
//        chain->SetBranchStatus("Muon_Best_phi", 1);

//        chain->SetBranchStatus("Muon_Inner_pT", 1);
//        chain->SetBranchStatus("Muon_Inner_pTError", 1);
//        chain->SetBranchStatus("Muon_Inner_eta", 1);
//        chain->SetBranchStatus("Muon_Inner_phi", 1);
//        chain->SetBranchStatus("Muon_Inner_Px", 1);
//        chain->SetBranchStatus("Muon_Inner_Py", 1);
//        chain->SetBranchStatus("Muon_Inner_Pz", 1);

//        chain->SetBranchStatus("Muon_Outer_pT", 1);
//        chain->SetBranchStatus("Muon_Outer_pTError", 1);
//        chain->SetBranchStatus("Muon_Outer_Px", 1);
//        chain->SetBranchStatus("Muon_Outer_Py", 1);
//        chain->SetBranchStatus("Muon_Outer_Pz", 1);
//        chain->SetBranchStatus("Muon_Outer_eta", 1);
//        chain->SetBranchStatus("Muon_Outer_phi", 1);

//        chain->SetBranchStatus("Muon_GLB_pT", 1);
//        chain->SetBranchStatus("Muon_GLB_pTError", 1);
//        chain->SetBranchStatus("Muon_GLB_Px", 1);
//        chain->SetBranchStatus("Muon_GLB_Py", 1);
//        chain->SetBranchStatus("Muon_GLB_Pz", 1);
//        chain->SetBranchStatus("Muon_GLB_eta", 1);
//        chain->SetBranchStatus("Muon_GLB_phi", 1);

//        chain->SetBranchStatus("Muon_TuneP_pT", 1);
//        chain->SetBranchStatus("Muon_TuneP_pTError", 1);
//        chain->SetBranchStatus("Muon_TuneP_eta", 1);
//        chain->SetBranchStatus("Muon_TuneP_phi", 1);
//        chain->SetBranchStatus("Muon_TuneP_Px", 1);
//        chain->SetBranchStatus("Muon_TuneP_Py", 1);
//        chain->SetBranchStatus("Muon_TuneP_Pz", 1);

//        chain->SetBranchStatus("Muon_stationMask", 1);
//        chain->SetBranchStatus("Muon_nMatchesRPCLayers", 1);


    	chain->SetBranchAddress("nMuon", &nMuon);
    	chain->SetBranchAddress("Muon_pT", Muon_pT);
    	chain->SetBranchAddress("Muon_eta", Muon_eta);
    	chain->SetBranchAddress("Muon_phi", Muon_phi);
//    	chain->SetBranchAddress("Muon_muonType", Muon_muonType);
//    	chain->SetBranchAddress("Muon_chi2dof", Muon_chi2dof);
//    	chain->SetBranchAddress("Muon_muonHits", Muon_muonHits);
//    	chain->SetBranchAddress("Muon_nSegments", Muon_nSegments);
//    	chain->SetBranchAddress("Muon_nMatches", Muon_nMatches);
//    	chain->SetBranchAddress("Muon_trackerLayers", Muon_trackerLayers);
        // chain->SetBranchAddress("Muon_trackerHitsGLB", Muon_trackerHitsGLB);
//        chain->SetBranchAddress("Muon_pixelHitsGLB", Muon_pixelHitsGLB);
//        chain->SetBranchAddress("Muon_trackerLayersGLB", Muon_trackerLayersGLB);
//    	chain->SetBranchAddress("Muon_pixelHits", Muon_pixelHits);
//    	chain->SetBranchAddress("Muon_dxyVTX", Muon_dxyVTX);
//    	chain->SetBranchAddress("Muon_dzVTX", Muon_dzVTX);
//    	chain->SetBranchAddress("Muon_trkiso", Muon_trkiso);
//    	chain->SetBranchAddress("isPFmuon", isPFmuon);
//    	chain->SetBranchAddress("isGLBmuon", isGLBmuon);
//        chain->SetBranchAddress("isTRKmuon", isTRKmuon);
    	chain->SetBranchAddress("Muon_Px", Muon_Px );
    	chain->SetBranchAddress("Muon_Py", Muon_Py );
    	chain->SetBranchAddress("Muon_Pz", Muon_Pz );

//        chain->SetBranchAddress("Muon_dB", Muon_dB );
    	
    	chain->SetBranchAddress("Muon_charge", Muon_charge);
    	
//    	chain->SetBranchAddress("Muon_PfChargedHadronIsoR04", Muon_PfChargedHadronIsoR04);
//    	chain->SetBranchAddress("Muon_PfNeutralHadronIsoR04", Muon_PfNeutralHadronIsoR04);
//    	chain->SetBranchAddress("Muon_PfGammaIsoR04", Muon_PfGammaIsoR04);
//        chain->SetBranchAddress("Muon_PFSumPUIsoR04", Muon_PFSumPUIsoR04);

//    	chain->SetBranchAddress("CosAngle", &CosAngle);
//    	chain->SetBranchAddress("vtxTrkChi2", &vtxTrkChi2);
//    	chain->SetBranchAddress("vtxTrkProb", &vtxTrkProb);
//    	chain->SetBranchAddress("vtxTrkNdof", &vtxTrkNdof);
//    	chain->SetBranchAddress("vtxTrkCkt1Pt", &vtxTrkCkt1Pt);
//    	chain->SetBranchAddress("vtxTrkCkt2Pt", &vtxTrkCkt2Pt);
//    	chain->SetBranchAddress("vtxTrkDiEChi2", &vtxTrkDiEChi2);
//    	chain->SetBranchAddress("vtxTrkDiEProb", &vtxTrkDiEProb);
//    	chain->SetBranchAddress("vtxTrkDiENdof", &vtxTrkDiENdof);
//    	chain->SetBranchAddress("vtxTrkDiE1Pt", &vtxTrkDiE1Pt);
//    	chain->SetBranchAddress("vtxTrkDiE2Pt", &vtxTrkDiE2Pt);
//    	chain->SetBranchAddress("vtxTrkEMuChi2", &vtxTrkEMuChi2);
//    	chain->SetBranchAddress("vtxTrkEMuProb", &vtxTrkEMuProb);
//    	chain->SetBranchAddress("vtxTrkEMuNdof", &vtxTrkEMuNdof);
//    	chain->SetBranchAddress("vtxTrkEMu1Pt", &vtxTrkEMu1Pt);
//    	chain->SetBranchAddress("vtxTrkEMu2Pt", &vtxTrkEMu2Pt);

//        chain->SetBranchAddress("CosAngle_TuneP", &CosAngle_TuneP);
//        chain->SetBranchAddress("vtxTrk1Pt_TuneP", &vtxTrk1Pt_TuneP);
//        chain->SetBranchAddress("vtxTrk2Pt_TuneP", &vtxTrk2Pt_TuneP);
//        chain->SetBranchAddress("vtxTrkChi2_TuneP", &vtxTrkChi2_TuneP);
//        chain->SetBranchAddress("vtxTrkNdof_TuneP", &vtxTrkNdof_TuneP);
//        chain->SetBranchAddress("vtxTrkProb_TuneP", &vtxTrkProb_TuneP);
    
//    	chain->SetBranchAddress("Muon_Best_pT", &Muon_Best_pT);
//    	chain->SetBranchAddress("Muon_Best_pTError", &Muon_Best_pTError);
//    	chain->SetBranchAddress("Muon_Best_Px", &Muon_Best_Px);
//    	chain->SetBranchAddress("Muon_Best_Py", &Muon_Best_Py);
//    	chain->SetBranchAddress("Muon_Best_Pz", &Muon_Best_Pz);
//    	chain->SetBranchAddress("Muon_Best_eta", &Muon_Best_eta);
//    	chain->SetBranchAddress("Muon_Best_phi", &Muon_Best_phi);

//    	chain->SetBranchAddress("Muon_Inner_pT", &Muon_Inner_pT);
//    	chain->SetBranchAddress("Muon_Inner_pTError", &Muon_Inner_pTError);
//    	chain->SetBranchAddress("Muon_Inner_Px", &Muon_Inner_Px);
//    	chain->SetBranchAddress("Muon_Inner_Py", &Muon_Inner_Py);
//    	chain->SetBranchAddress("Muon_Inner_Pz", &Muon_Inner_Pz);
//    	chain->SetBranchAddress("Muon_Inner_eta", &Muon_Inner_eta);
//    	chain->SetBranchAddress("Muon_Inner_phi", &Muon_Inner_phi);

//    	chain->SetBranchAddress("Muon_Outer_pT", &Muon_Outer_pT);
//    	chain->SetBranchAddress("Muon_Outer_pTError", &Muon_Outer_pTError);
//    	chain->SetBranchAddress("Muon_Outer_Px", &Muon_Outer_Px);
//    	chain->SetBranchAddress("Muon_Outer_Py", &Muon_Outer_Py);
//    	chain->SetBranchAddress("Muon_Outer_Pz", &Muon_Outer_Pz);
//    	chain->SetBranchAddress("Muon_Outer_eta", &Muon_Outer_eta);
//    	chain->SetBranchAddress("Muon_Outer_phi", &Muon_Outer_phi);

//    	chain->SetBranchAddress("Muon_GLB_pT", &Muon_GLB_pT);
//    	chain->SetBranchAddress("Muon_GLB_pTError", &Muon_GLB_pTError);
//    	chain->SetBranchAddress("Muon_GLB_Px", &Muon_GLB_Px);
//    	chain->SetBranchAddress("Muon_GLB_Py", &Muon_GLB_Py);
//    	chain->SetBranchAddress("Muon_GLB_Pz", &Muon_GLB_Pz);
//    	chain->SetBranchAddress("Muon_GLB_eta", &Muon_GLB_eta);
//    	chain->SetBranchAddress("Muon_GLB_phi", &Muon_GLB_phi);

//    	chain->SetBranchAddress("Muon_TuneP_pT", &Muon_TuneP_pT);
//    	chain->SetBranchAddress("Muon_TuneP_pTError", &Muon_TuneP_pTError);
//    	chain->SetBranchAddress("Muon_TuneP_Px", &Muon_TuneP_Px);
//    	chain->SetBranchAddress("Muon_TuneP_Py", &Muon_TuneP_Py);
//    	chain->SetBranchAddress("Muon_TuneP_Pz", &Muon_TuneP_Pz);
//    	chain->SetBranchAddress("Muon_TuneP_eta", &Muon_TuneP_eta);
//    	chain->SetBranchAddress("Muon_TuneP_phi", &Muon_TuneP_phi);

//        chain->SetBranchAddress("Muon_stationMask", &Muon_stationMask);
//        chain->SetBranchAddress("Muon_nMatchesRPCLayers", &Muon_nMatchesRPCLayers);

//        chain->AddBranchToCache("isPFmuon", 1);
//        chain->AddBranchToCache("isGLBmuon", 1);
//        chain->AddBranchToCache("isTRKmuon", 1);
//        chain->AddBranchToCache("CosAngle", 1);
//        chain->AddBranchToCache("vtxTrkChi2", 1);
//        chain->AddBranchToCache("vtxTrkProb", 1);
//        chain->AddBranchToCache("vtxTrkNdof", 1);
//        chain->AddBranchToCache("vtxTrkCkt1Pt", 1);
//        chain->AddBranchToCache("vtxTrkCkt2Pt", 1);
//        chain->AddBranchToCache("vtxTrkDiEChi2", 1);
//        chain->AddBranchToCache("vtxTrkDiEProb", 1);
//        chain->AddBranchToCache("vtxTrkDiENdof", 1);
//        chain->AddBranchToCache("vtxTrkDiE1Pt", 1);
//        chain->AddBranchToCache("vtxTrkDiE2Pt", 1);
//        chain->AddBranchToCache("vtxTrkEMuChi2", 1);
//        chain->AddBranchToCache("vtxTrkEMuProb", 1);
//        chain->AddBranchToCache("vtxTrkEMuNdof", 1);
//        chain->AddBranchToCache("vtxTrkEMu1Pt", 1);
//        chain->AddBranchToCache("vtxTrkEMu2Pt", 1);

//        chain->AddBranchToCache("CosAngle_TuneP", 1);
//        chain->AddBranchToCache("vtxTrk1Pt_TuneP", 1);
//        chain->AddBranchToCache("vtxTrk2Pt_TuneP", 1);
//        chain->AddBranchToCache("vtxTrkChi2_TuneP", 1);
//        chain->AddBranchToCache("vtxTrkNdof_TuneP", 1);
//        chain->AddBranchToCache("vtxTrkProb_TuneP", 1);

        chain->AddBranchToCache("nMuon", 1);
        chain->AddBranchToCache("Muon_pT", 1);
        chain->AddBranchToCache("Muon_eta", 1);
        chain->AddBranchToCache("Muon_phi", 1);
//        chain->AddBranchToCache("Muon_muonType", 1);
//        chain->AddBranchToCache("Muon_chi2dof", 1);
//        chain->AddBranchToCache("Muon_muonHits", 1);
//        chain->AddBranchToCache("Muon_nSegments", 1);
//        chain->AddBranchToCache("Muon_nMatches", 1);
//        chain->AddBranchToCache("Muon_trackerLayers", 1);

        // chain->AddBranchToCache("Muon_trackerHitsGLB", 1);
//        chain->AddBranchToCache("Muon_pixelHitsGLB", 1);
//        chain->AddBranchToCache("Muon_trackerLayersGLB", 1);

//        chain->AddBranchToCache("Muon_pixelHits", 1);
//        chain->AddBranchToCache("Muon_dxyVTX", 1);
//        chain->AddBranchToCache("Muon_dzVTX", 1);
//        chain->AddBranchToCache("Muon_trkiso", 1);
        
        chain->AddBranchToCache("Muon_Px", 1);
        chain->AddBranchToCache("Muon_Py", 1);
        chain->AddBranchToCache("Muon_Pz", 1);

//        chain->AddBranchToCache("Muon_dB", 1);
        
        chain->AddBranchToCache("Muon_charge", 1);
        
//        chain->AddBranchToCache("Muon_PfChargedHadronIsoR04", 1);
//        chain->AddBranchToCache("Muon_PfNeutralHadronIsoR04", 1);
//        chain->AddBranchToCache("Muon_PfGammaIsoR04", 1);
//        chain->AddBranchToCache("Muon_PFSumPUIsoR04", 1);

//        chain->AddBranchToCache("Muon_Best_pT", 1);
//        chain->AddBranchToCache("Muon_Best_pTError", 1);
//        chain->AddBranchToCache("Muon_Best_Px", 1);
//        chain->AddBranchToCache("Muon_Best_Py", 1);
//        chain->AddBranchToCache("Muon_Best_Pz", 1);
//        chain->AddBranchToCache("Muon_Best_eta", 1);
//        chain->AddBranchToCache("Muon_Best_phi", 1);

//        chain->AddBranchToCache("Muon_Inner_pT", 1);
//        chain->AddBranchToCache("Muon_Inner_pTError", 1);
//        chain->AddBranchToCache("Muon_Inner_eta", 1);
//        chain->AddBranchToCache("Muon_Inner_phi", 1);
//        chain->AddBranchToCache("Muon_Inner_Px", 1);
//        chain->AddBranchToCache("Muon_Inner_Py", 1);
//        chain->AddBranchToCache("Muon_Inner_Pz", 1);

//        chain->AddBranchToCache("Muon_Outer_pT", 1);
//        chain->AddBranchToCache("Muon_Outer_pTError", 1);
//        chain->AddBranchToCache("Muon_Outer_Px", 1);
//        chain->AddBranchToCache("Muon_Outer_Py", 1);
//        chain->AddBranchToCache("Muon_Outer_Pz", 1);
//        chain->AddBranchToCache("Muon_Outer_eta", 1);
//        chain->AddBranchToCache("Muon_Outer_phi", 1);

//        chain->AddBranchToCache("Muon_GLB_pT", 1);
//        chain->AddBranchToCache("Muon_GLB_pTError", 1);
//        chain->AddBranchToCache("Muon_GLB_Px", 1);
//        chain->AddBranchToCache("Muon_GLB_Py", 1);
//        chain->AddBranchToCache("Muon_GLB_Pz", 1);
//        chain->AddBranchToCache("Muon_GLB_eta", 1);
//        chain->AddBranchToCache("Muon_GLB_phi", 1);

//        chain->AddBranchToCache("Muon_TuneP_pT", 1);
//        chain->AddBranchToCache("Muon_TuneP_pTError", 1);
//        chain->AddBranchToCache("Muon_TuneP_eta", 1);
//        chain->AddBranchToCache("Muon_TuneP_phi", 1);
//        chain->AddBranchToCache("Muon_TuneP_Px", 1);
//        chain->AddBranchToCache("Muon_TuneP_Py", 1);
//        chain->AddBranchToCache("Muon_TuneP_Pz", 1);

//        chain->AddBranchToCache("Muon_stationMask", 1);
//        chain->AddBranchToCache("Muon_nMatchesRPCLayers", 1);

    }

    void TurnOnBranches_Electron()
    {
        chain->SetBranchStatus("Nelectrons", 1);
        chain->SetBranchStatus("Electron_Energy", 1);
        chain->SetBranchStatus("Electron_pT", 1);
        chain->SetBranchStatus("Electron_eta", 1);
        chain->SetBranchStatus("Electron_phi", 1);
        chain->SetBranchStatus("Electron_charge", 1);
//        chain->SetBranchStatus("Electron_gsfpT", 1);
//        chain->SetBranchStatus("Electron_gsfPx", 1);
//        chain->SetBranchStatus("Electron_gsfPy", 1);
//        chain->SetBranchStatus("Electron_gsfPz", 1);
//        chain->SetBranchStatus("Electron_gsfEta", 1);
//        chain->SetBranchStatus("Electron_gsfPhi", 1);
//        chain->SetBranchStatus("Electron_gsfCharge", 1);
//        chain->SetBranchStatus("Electron_etaSC", 1);
//        chain->SetBranchStatus("Electron_phiSC", 1);
//        chain->SetBranchStatus("Electron_etaWidth", 1);
//        chain->SetBranchStatus("Electron_phiWidth", 1);
//        chain->SetBranchStatus("Electron_dEtaIn", 1);
//        chain->SetBranchStatus("Electron_dPhiIn", 1);
//        chain->SetBranchStatus("Electron_sigmaIEtaIEta", 1);
//        chain->SetBranchStatus("Electron_HoverE", 1);
//        chain->SetBranchStatus("Electron_fbrem", 1);
//        chain->SetBranchStatus("Electron_eOverP", 1);
//        chain->SetBranchStatus("Electron_InvEminusInvP", 1);
//        chain->SetBranchStatus("Electron_dxyVTX", 1);
//        chain->SetBranchStatus("Electron_dzVTX", 1);
//        chain->SetBranchStatus("Electron_dxy", 1);
//        chain->SetBranchStatus("Electron_dz", 1);
//        chain->SetBranchStatus("Electron_dxyBS", 1);
//        chain->SetBranchStatus("Electron_dzBS", 1);
//        chain->SetBranchStatus("Electron_chIso03", 1);
//        chain->SetBranchStatus("Electron_nhIso03", 1);
//        chain->SetBranchStatus("Electron_phIso03", 1);
//        chain->SetBranchStatus("Electron_ChIso03FromPU", 1);

//        chain->SetBranchStatus("Electron_mHits", 1);
//        chain->SetBranchStatus("Electron_EnergySC", 1);
//        chain->SetBranchStatus("Electron_preEnergySC", 1);
//        chain->SetBranchStatus("Electron_rawEnergySC", 1);
//        chain->SetBranchStatus("Electron_etSC", 1);
//        chain->SetBranchStatus("Electron_E15", 1);
//        chain->SetBranchStatus("Electron_E25", 1);
//        chain->SetBranchStatus("Electron_E55", 1);
//        chain->SetBranchStatus("Electron_RelPFIso_dBeta", 1);
//        chain->SetBranchStatus("Electron_RelPFIso_Rho", 1);
//        chain->SetBranchStatus("Electron_r9", 1);
//        chain->SetBranchStatus("Electron_ecalDriven", 1);
//        chain->SetBranchStatus("Electron_passConvVeto", 1);

//        chain->SetBranchStatus("Electron_passLooseID", 1);
//        chain->SetBranchStatus("Electron_passMediumID", 1);
//        chain->SetBranchStatus("Electron_passTightID", 1);
//        chain->SetBranchStatus("Electron_passMVAID_WP80", 1);
//        chain->SetBranchStatus("Electron_passMVAID_WP90", 1);
//        chain->SetBranchStatus("Electron_passHEEPID", 1);

    	chain->SetBranchAddress("Nelectrons", &Nelectrons);
    	chain->SetBranchAddress("Electron_Energy", &Electron_Energy);
    	chain->SetBranchAddress("Electron_pT", &Electron_pT);
    	chain->SetBranchAddress("Electron_eta", &Electron_eta);
    	chain->SetBranchAddress("Electron_phi", &Electron_phi);
    	chain->SetBranchAddress("Electron_charge", &Electron_charge);
//    	chain->SetBranchAddress("Electron_gsfpT", &Electron_gsfpT);
//    	chain->SetBranchAddress("Electron_gsfPx", &Electron_gsfPx);
//    	chain->SetBranchAddress("Electron_gsfPy", &Electron_gsfPy);
//    	chain->SetBranchAddress("Electron_gsfPz", &Electron_gsfPz);
//    	chain->SetBranchAddress("Electron_gsfEta", &Electron_gsfEta);
//    	chain->SetBranchAddress("Electron_gsfPhi", &Electron_gsfPhi);
//    	chain->SetBranchAddress("Electron_gsfCharge", &Electron_gsfCharge);
//    	chain->SetBranchAddress("Electron_etaSC", &Electron_etaSC);
//    	chain->SetBranchAddress("Electron_phiSC", &Electron_phiSC);
//    	chain->SetBranchAddress("Electron_etaWidth", &Electron_etaWidth);
//    	chain->SetBranchAddress("Electron_phiWidth", &Electron_phiWidth);
//    	chain->SetBranchAddress("Electron_dEtaIn", &Electron_dEtaIn);
//    	chain->SetBranchAddress("Electron_dPhiIn", &Electron_dPhiIn);
//    	chain->SetBranchAddress("Electron_sigmaIEtaIEta", &Electron_sigmaIEtaIEta);
//    	chain->SetBranchAddress("Electron_HoverE", &Electron_HoverE);
//    	chain->SetBranchAddress("Electron_fbrem", &Electron_fbrem);
//    	chain->SetBranchAddress("Electron_eOverP", &Electron_eOverP);
//    	chain->SetBranchAddress("Electron_InvEminusInvP", &Electron_InvEminusInvP);
//    	chain->SetBranchAddress("Electron_dxyVTX", &Electron_dxyVTX);
//    	chain->SetBranchAddress("Electron_dzVTX", &Electron_dzVTX);
//    	chain->SetBranchAddress("Electron_dxy", &Electron_dxy);
//    	chain->SetBranchAddress("Electron_dz", &Electron_dz);
//    	chain->SetBranchAddress("Electron_dxyBS", &Electron_dxyBS);
//    	chain->SetBranchAddress("Electron_dzBS", &Electron_dzBS);
//    	chain->SetBranchAddress("Electron_chIso03", &Electron_chIso03);
//    	chain->SetBranchAddress("Electron_nhIso03", &Electron_nhIso03);
//    	chain->SetBranchAddress("Electron_phIso03", &Electron_phIso03);
//    	chain->SetBranchAddress("Electron_ChIso03FromPU", &Electron_ChIso03FromPU);

//    	chain->SetBranchAddress("Electron_mHits", &Electron_mHits);
//    	chain->SetBranchAddress("Electron_EnergySC", &Electron_EnergySC);
//    	chain->SetBranchAddress("Electron_preEnergySC", &Electron_preEnergySC);
//    	chain->SetBranchAddress("Electron_rawEnergySC", &Electron_rawEnergySC);
//    	chain->SetBranchAddress("Electron_etSC", &Electron_etSC);
//    	chain->SetBranchAddress("Electron_E15", &Electron_E15);
//    	chain->SetBranchAddress("Electron_E25", &Electron_E25);
//    	chain->SetBranchAddress("Electron_E55", &Electron_E55);
//    	chain->SetBranchAddress("Electron_RelPFIso_dBeta", &Electron_RelPFIso_dBeta);
//    	chain->SetBranchAddress("Electron_RelPFIso_Rho", &Electron_RelPFIso_Rho);
//    	chain->SetBranchAddress("Electron_r9", &Electron_r9);
//    	chain->SetBranchAddress("Electron_ecalDriven", &Electron_ecalDriven);
//        chain->SetBranchAddress("Electron_passConvVeto", &Electron_passConvVeto);

//        chain->SetBranchAddress("Electron_passLooseID", &Electron_passLooseID);
//        chain->SetBranchAddress("Electron_passMediumID", &Electron_passMediumID);
//        chain->SetBranchAddress("Electron_passTightID", &Electron_passTightID);
//        chain->SetBranchAddress("Electron_passMVAID_WP80", &Electron_passMVAID_WP80);
//        chain->SetBranchAddress("Electron_passMVAID_WP90", &Electron_passMVAID_WP90);
//        chain->SetBranchAddress("Electron_passHEEPID", &Electron_passHEEPID);

        chain->AddBranchToCache("Nelectrons", 1);
        chain->AddBranchToCache("Electron_Energy", 1);
        chain->AddBranchToCache("Electron_pT", 1);
        chain->AddBranchToCache("Electron_eta", 1);
        chain->AddBranchToCache("Electron_phi", 1);
        chain->AddBranchToCache("Electron_charge", 1);
//        chain->AddBranchToCache("Electron_gsfpT", 1);
//        chain->AddBranchToCache("Electron_gsfPx", 1);
//        chain->AddBranchToCache("Electron_gsfPy", 1);
//        chain->AddBranchToCache("Electron_gsfPz", 1);
//        chain->AddBranchToCache("Electron_gsfEta", 1);
//        chain->AddBranchToCache("Electron_gsfPhi", 1);
//        chain->AddBranchToCache("Electron_gsfCharge", 1);
//        chain->AddBranchToCache("Electron_etaSC", 1);
//        chain->AddBranchToCache("Electron_phiSC", 1);
//        chain->AddBranchToCache("Electron_etaWidth", 1);
//        chain->AddBranchToCache("Electron_phiWidth", 1);
//        chain->AddBranchToCache("Electron_dEtaIn", 1);
//        chain->AddBranchToCache("Electron_dPhiIn", 1);
//        chain->AddBranchToCache("Electron_sigmaIEtaIEta", 1);
//        chain->AddBranchToCache("Electron_HoverE", 1);
//        chain->AddBranchToCache("Electron_fbrem", 1);
//        chain->AddBranchToCache("Electron_eOverP", 1);
//        chain->AddBranchToCache("Electron_InvEminusInvP", 1);
//        chain->AddBranchToCache("Electron_dxyVTX", 1);
//        chain->AddBranchToCache("Electron_dzVTX", 1);
//        chain->AddBranchToCache("Electron_dxy", 1);
//        chain->AddBranchToCache("Electron_dz", 1);
//        chain->AddBranchToCache("Electron_dxyBS", 1);
//        chain->AddBranchToCache("Electron_dzBS", 1);
//        chain->AddBranchToCache("Electron_chIso03", 1);
//        chain->AddBranchToCache("Electron_nhIso03", 1);
//        chain->AddBranchToCache("Electron_phIso03", 1);
//        chain->AddBranchToCache("Electron_ChIso03FromPU", 1);

//        chain->AddBranchToCache("Electron_mHits", 1);
//        chain->AddBranchToCache("Electron_EnergySC", 1);
//        chain->AddBranchToCache("Electron_preEnergySC", 1);
//        chain->AddBranchToCache("Electron_rawEnergySC", 1);
//        chain->AddBranchToCache("Electron_etSC", 1);
//        chain->AddBranchToCache("Electron_E15", 1);
//        chain->AddBranchToCache("Electron_E25", 1);
//        chain->AddBranchToCache("Electron_E55", 1);
//        chain->AddBranchToCache("Electron_RelPFIso_dBeta", 1);
//        chain->AddBranchToCache("Electron_RelPFIso_Rho", 1);
//        chain->AddBranchToCache("Electron_r9", 1);
//        chain->AddBranchToCache("Electron_ecalDriven", 1);
//        chain->AddBranchToCache("Electron_passConvVeto", 1);

//        chain->AddBranchToCache("Electron_passLooseID", 1);
//        chain->AddBranchToCache("Electron_passMediumID", 1);
//        chain->AddBranchToCache("Electron_passTightID", 1);
//        chain->AddBranchToCache("Electron_passMVAID_WP80", 1);
//        chain->AddBranchToCache("Electron_passMVAID_WP90", 1);
//        chain->AddBranchToCache("Electron_passHEEPID", 1);
    }


    void GetEvent(Int_t i)
    {
        if(!chain) return;
        
        chain->GetEntry(i);
    }

    void ActivateBranches(TString brNames) {
      if (brNames.Length()==0) return;
      std::stringstream ss(brNames.Data());
      TString br;
      while (!ss.eof()) {
        ss >> br;
        chain->SetBranchStatus(br,1);
        std::cout << "activating branch <" << br << ">\n";
      }
    }

    Bool_t isTriggered(TString HLT)
    {
        Bool_t isTrigger = false;
        if( HLT == "HLT_IsoMu20_v* || HLT_IsoTkMu20_v*" )
        {
            for( Int_t k = 0; k < HLT_ntrig; k++ )
            {
                if( (HLT_trigName->at((unsigned int)k)) == "HLT_IsoMu20_v*" || (HLT_trigName->at((unsigned int)k)) == "HLT_IsoTkMu20_v*" )
                {
                    if( HLT_trigFired[k] == 1 )
                    {
                        isTrigger = true;
                        break;
                    }
                }
            }
        }
        else if( HLT == "HLT_IsoMu27_v* || HLT_IsoTkMu27_v*" )
        {
            for( Int_t k = 0; k < HLT_ntrig; k++ )
            {
                if( (HLT_trigName->at((unsigned int)k)) == "HLT_IsoMu27_v*" || (HLT_trigName->at((unsigned int)k)) == "HLT_IsoTkMu27_v*" )
                {
                    if( HLT_trigFired[k] == 1 )
                    {
                        isTrigger = true;
                        break;
                    }
                }
            }
        }
        else if( HLT == "HLT_IsoMu24_v* || HLT_IsoTkMu24_v*" )
        {
            for( Int_t k = 0; k < HLT_ntrig; k++ )
            {
                if( (HLT_trigName->at((unsigned int)k)) == "HLT_IsoMu24_v*" || (HLT_trigName->at((unsigned int)k)) == "HLT_IsoTkMu24_v*" )
                {
                    if( HLT_trigFired[k] == 1 )
                    {
                        isTrigger = true;
                        break;
                    }
                }
            }
        }
        else if( HLT == "HLT_Mu50_v* || HLT_TkMu50_v*" )
        {
            for( Int_t k = 0; k < HLT_ntrig; k++ )
            {
                if( (HLT_trigName->at((unsigned int)k)) == "HLT_Mu50_v*" || (HLT_trigName->at((unsigned int)k)) == "HLT_TkMu50_v*" )
                {
                    if( HLT_trigFired[k] == 1 )
                    {
                        isTrigger = true;
                        break;
                    }
                }
            }
        }
        else
        {
            for( Int_t k = 0; k < HLT_ntrig; k++ )
            {
                if( (HLT_trigName->at((unsigned int)k)) == HLT )
                {
                    if( HLT_trigFired[k] == 1 )
                    {
                        isTrigger = true;
                        break;
                    }
                }
            }
        }

        return isTrigger;
    }

};
