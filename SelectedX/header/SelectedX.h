///////////////////////////////////////////////////////////////////////////////////////
/// 2018.07.03: First version created by Marijus Ambrozas
/// 2018.07.27: Divided SelectedEE and SelectedMuMu into LongSelected* and Selected*
/// 2018.08.01: Created SelectedEMu and LongSelectedEMu classes, edited descriptions
///////////////////////////////////////////////////////////////////////////////////////
#pragma once

#include <TTree.h>
#include <TTreeCache.h>
#include <TChain.h>
#include <vector>

class SelectedX
{
public:
    TChain *chain;
    Bool_t File_Given = kFALSE;

    Double_t GENEvt_weight;
    Int_t nPileUp;
    Int_t nVertices;
};

///////////////////////////////////////////////////////////////////////////////////
///                                                                             ///
///     The LongSelectedMuMu_t class                                            ///
///                                                                             ///
///     Stores information about a pair of muons that passed the DY->MuMu       ///
///     selection. It is possible to repeat the full event selection using      ///
///     this class instead of regular ntuple (except for the IsHardProcess).    ///
///                                                                             ///
///     How to use (examples):                                                  ///
///                                                                             ///
///         To read from file:                                                  ///
///             > #include "SelectedX.h"                                        ///
///             > TChain *ch = new TChain("...");                               ///
///             > ch->Add("...");                                               ///
///             > LongSelectedMuMu_t mu;                                        ///
///             > mu.CreateFromChain(ch);                                       ///
///             > mu.GetEvent(...);                                             ///
///             > cout << mu.Muon_Px->at(...);                                  ///
///                                                                             ///
///         To create an empty class (e.g. to fill with values and write        ///
///         into a tree):                                                       ///
///             > #include "SelectedX.h"                                        ///
///             > LongSelectedMuMu_t mu;                                        ///
///             > mu.CreateNew();                                               ///
///             > mu.Muon_Py->push_back(...);                                   ///
///                                                                             ///
///////////////////////////////////////////////////////////////////////////////////
class LongSelectedMuMu_t : public SelectedX
{
public:
    //Event Information
    Int_t runNum;
    Int_t lumiBlock;
    Int_t evtNum;

    Bool_t isHardProcess;

    //Trigger variables
    Int_t HLT_ntrig;
    std::vector<int> *HLT_trigFired;
    std::vector<string> *HLT_trigName;
//    std::vector<double> *HLT_trigPt;
    std::vector<double> *HLT_trigEta;
    std::vector<double> *HLT_trigPhi;

    Double_t Muon_InvM;

    // -- Physical Variables -- //
    std::vector<double> *Muon_pT;
    std::vector<double> *Muon_eta;
    std::vector<double> *Muon_phi;
    std::vector<double> *Muon_Energy;
    std::vector<int> *Muon_charge;

    // -- Cut variables -- //
//    std::vector<int> *Muon_muonType;
    std::vector<double> *Muon_chi2dof;
    std::vector<int> *Muon_muonHits;
    std::vector<int> *Muon_nSegments;
    std::vector<int> *Muon_nMatches;
    std::vector<int> *Muon_trackerLayers;

//    std::vector<int> *Muon_trackerHitsGLB;
//    std::vector<int> *Muon_pixelHitsGLB;
//    std::vector<int> *Muon_trackerLayersGLB;

    std::vector<int> *Muon_pixelHits;
    std::vector<double> *Muon_dxyVTX;
    std::vector<double> *Muon_dzVTX;
    std::vector<double> *Muon_trkiso;
    std::vector<int> *isGLBmuon;
    std::vector<int> *isPFmuon;
    std::vector<int> *isTRKmuon;

    // -- for invariant mass calculation -- //
    std::vector<double> *Muon_Px;
    std::vector<double> *Muon_Py;
    std::vector<double> *Muon_Pz;

//    Double_t Muon_dB[MaxN];

    //PF information
    std::vector<double> *Muon_PfChargedHadronIsoR04;
    std::vector<double> *Muon_PfNeutralHadronIsoR04;
    std::vector<double> *Muon_PfGammaIsoR04;
    std::vector<double> *Muon_PFSumPUIsoR04;

    //Dimuon variables
    std::vector<double> *CosAngle;
    std::vector<double> *vtxTrkChi2;
    std::vector<double> *vtxTrkProb;
    std::vector<double> *vtxTrkNdof;
    std::vector<double> *vtxTrkCkt1Pt;
    std::vector<double> *vtxTrkCkt2Pt;
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

    // -- Various Track Information -- //
    std::vector<double> *Muon_Best_pT;
    std::vector<double> *Muon_Best_pTError;
    std::vector<double> *Muon_Best_Px;
    std::vector<double> *Muon_Best_Py;
    std::vector<double> *Muon_Best_Pz;
    std::vector<double> *Muon_Best_eta;
    std::vector<double> *Muon_Best_phi;

    std::vector<double> *Muon_Inner_pT;
    std::vector<double> *Muon_Inner_pTError;
    std::vector<double> *Muon_Inner_Px;
    std::vector<double> *Muon_Inner_Py;
    std::vector<double> *Muon_Inner_Pz;
    std::vector<double> *Muon_Inner_eta;
    std::vector<double> *Muon_Inner_phi;

    std::vector<double> *Muon_Outer_pT;
    std::vector<double> *Muon_Outer_pTError;
    std::vector<double> *Muon_Outer_Px;
    std::vector<double> *Muon_Outer_Py;
    std::vector<double> *Muon_Outer_Pz;
    std::vector<double> *Muon_Outer_eta;
    std::vector<double> *Muon_Outer_phi;

    std::vector<double> *Muon_GLB_pT;
    std::vector<double> *Muon_GLB_pTError;
    std::vector<double> *Muon_GLB_Px;
    std::vector<double> *Muon_GLB_Py;
    std::vector<double> *Muon_GLB_Pz;
    std::vector<double> *Muon_GLB_eta;
    std::vector<double> *Muon_GLB_phi;

    std::vector<double> *Muon_TuneP_pT;
    std::vector<double> *Muon_TuneP_pTError;
    std::vector<double> *Muon_TuneP_Px;
    std::vector<double> *Muon_TuneP_Py;
    std::vector<double> *Muon_TuneP_Pz;
    std::vector<double> *Muon_TuneP_eta;
    std::vector<double> *Muon_TuneP_phi;

//    std::vector<int> *Muon_stationMask;
//    std::vector<int> *Muon_nMatchesRPCLayers;

    // -- Default constructor -- //

    void CreateNew()
    {
        HLT_trigFired = new std::vector<int>;
        HLT_trigName = new std::vector<string>;
//        HLT_trigPt = new std::vector<double>;
        HLT_trigEta = new std::vector<double>;
        HLT_trigPhi = new std::vector<double>;

        Muon_pT = new std::vector<double>;
        Muon_eta = new std::vector<double>;
        Muon_phi = new std::vector<double>;
        Muon_Energy = new std::vector<double>;
        Muon_charge = new std::vector<int>;

//        Muon_muonType = new std::vector<int>;
        Muon_chi2dof = new std::vector<double>;
        Muon_muonHits = new std::vector<int>;
        Muon_nSegments = new std::vector<int>;
        Muon_nMatches = new std::vector<int>;
        Muon_trackerLayers = new std::vector<int>;

//        Muon_trackerHitsGLB = new std::vector<int>;
//        Muon_pixelHitsGLB = new std::vector<int>;
//        Muon_trackerLayersGLB = new std::vector<int>;

        Muon_pixelHits = new std::vector<int>;
        Muon_dxyVTX = new std::vector<double>;
        Muon_dzVTX = new std::vector<double>;
        Muon_trkiso = new std::vector<double>;
        isGLBmuon = new std::vector<int>;
        isPFmuon = new std::vector<int>;
        isTRKmuon = new std::vector<int>;

        Muon_Px = new std::vector<double>;
        Muon_Py = new std::vector<double>;
        Muon_Pz = new std::vector<double>;

    //    Muon_dB = new std::vector<double>;        

        Muon_PfChargedHadronIsoR04 = new std::vector<double>;
        Muon_PfNeutralHadronIsoR04 = new std::vector<double>;
        Muon_PfGammaIsoR04 = new std::vector<double>;
        Muon_PFSumPUIsoR04 = new std::vector<double>;

        CosAngle = new std::vector<double>;
        vtxTrkChi2 = new std::vector<double>;
        vtxTrkProb = new std::vector<double>;
        vtxTrkNdof = new std::vector<double>;
        vtxTrkCkt1Pt = new std::vector<double>;
        vtxTrkCkt2Pt = new std::vector<double>;
//        vtxTrkDiEChi2 = new std::vector<double>;
//        vtxTrkDiEProb = new std::vector<double>;
//        vtxTrkDiENdof = new std::vector<double>;
//        vtxTrkDiE1Pt = new std::vector<double>;
//        vtxTrkDiE2Pt = new std::vector<double>;
//        vtxTrkEMuChi2 = new std::vector<double>;
//        vtxTrkEMuProb = new std::vector<double>;
//        vtxTrkEMuNdof = new std::vector<double>;
//        vtxTrkEMu1Pt = new std::vector<double>;
//        vtxTrkEMu2Pt = new std::vector<double>;

        Muon_Best_pT = new std::vector<double>;
        Muon_Best_pTError = new std::vector<double>;
        Muon_Best_Px = new std::vector<double>;
        Muon_Best_Py = new std::vector<double>;
        Muon_Best_Pz = new std::vector<double>;
        Muon_Best_eta = new std::vector<double>;
        Muon_Best_phi = new std::vector<double>;

        Muon_Inner_pT = new std::vector<double>;
        Muon_Inner_pTError = new std::vector<double>;
        Muon_Inner_Px = new std::vector<double>;
        Muon_Inner_Py = new std::vector<double>;
        Muon_Inner_Pz = new std::vector<double>;
        Muon_Inner_eta = new std::vector<double>;
        Muon_Inner_phi = new std::vector<double>;

        Muon_Outer_pT = new std::vector<double>;
        Muon_Outer_pTError = new std::vector<double>;
        Muon_Outer_Px = new std::vector<double>;
        Muon_Outer_Py = new std::vector<double>;
        Muon_Outer_Pz = new std::vector<double>;
        Muon_Outer_eta = new std::vector<double>;
        Muon_Outer_phi = new std::vector<double>;

        Muon_GLB_pT = new std::vector<double>;
        Muon_GLB_pTError = new std::vector<double>;
        Muon_GLB_Px = new std::vector<double>;
        Muon_GLB_Py = new std::vector<double>;
        Muon_GLB_Pz = new std::vector<double>;
        Muon_GLB_eta = new std::vector<double>;
        Muon_GLB_phi = new std::vector<double>;

        Muon_TuneP_pT = new std::vector<double>;
        Muon_TuneP_pTError = new std::vector<double>;
        Muon_TuneP_Px = new std::vector<double>;
        Muon_TuneP_Py = new std::vector<double>;
        Muon_TuneP_Pz = new std::vector<double>;
        Muon_TuneP_eta = new std::vector<double>;
        Muon_TuneP_phi = new std::vector<double>;

    //    Muon_stationMask = new std::vector<int>;
    //    Muon_nMatchesRPCLayers = new std::vector<int>;
    }

    // -- Constructor with chain -- //
    void CreateFromChain(TChain *chainptr)
    {
        chain = chainptr;
        Int_t cachesize = 500000000; // -- 500MB -- //
        chain->SetCacheSize(cachesize);

        // -- Setting addresses -- //
        chain->SetBranchAddress("nVertices", &nVertices);
        chain->SetBranchAddress("runNum", &runNum);
        chain->SetBranchAddress("lumiBlock", &lumiBlock);
        chain->SetBranchAddress("evtNum", &evtNum);
        chain->SetBranchAddress("nPileUp", &nPileUp);
        chain->SetBranchAddress("isHardProcess", &isHardProcess);
        chain->SetBranchAddress("GENEvt_weight", &GENEvt_weight);

        chain->SetBranchAddress("HLT_trigName", &HLT_trigName);
        chain->SetBranchAddress("HLT_trigFired", &HLT_trigFired);
        chain->SetBranchAddress("HLT_ntrig", &HLT_ntrig);
//        chain->SetBranchAddress("HLT_trigPt", &HLT_trigPt);
        chain->SetBranchAddress("HLT_trigEta", &HLT_trigEta);
        chain->SetBranchAddress("HLT_trigPhi", &HLT_trigPhi);

        chain->SetBranchAddress("Muon_InvM", &Muon_InvM);

        chain->SetBranchAddress("Muon_pT", &Muon_pT);
        chain->SetBranchAddress("Muon_eta", &Muon_eta);
        chain->SetBranchAddress("Muon_phi", &Muon_phi);
        chain->SetBranchAddress("Muon_Energy", &Muon_Energy);
        chain->SetBranchAddress("Muon_charge", &Muon_charge);
//        chain->SetBranchAddress("Muon_muonType", &Muon_muonType);
        chain->SetBranchAddress("Muon_chi2dof", &Muon_chi2dof);
        chain->SetBranchAddress("Muon_muonHits", &Muon_muonHits);
        chain->SetBranchAddress("Muon_nSegments", &Muon_nSegments);
        chain->SetBranchAddress("Muon_nMatches", &Muon_nMatches);
        chain->SetBranchAddress("Muon_trackerLayers", &Muon_trackerLayers);
//        chain->SetBranchAddress("Muon_pixelHitsGLB", &Muon_pixelHitsGLB);
//        chain->SetBranchAddress("Muon_trackerLayersGLB", &Muon_trackerLayersGLB);
        chain->SetBranchAddress("Muon_pixelHits", &Muon_pixelHits);
        chain->SetBranchAddress("Muon_dxyVTX", &Muon_dxyVTX);
        chain->SetBranchAddress("Muon_dzVTX", &Muon_dzVTX);
        chain->SetBranchAddress("Muon_trkiso", &Muon_trkiso);
        chain->SetBranchAddress("isPFmuon", &isPFmuon);
        chain->SetBranchAddress("isGLBmuon", &isGLBmuon);
        chain->SetBranchAddress("isTRKmuon", &isTRKmuon);
        chain->SetBranchAddress("Muon_Px", &Muon_Px );
        chain->SetBranchAddress("Muon_Py", &Muon_Py );
        chain->SetBranchAddress("Muon_Pz", &Muon_Pz );

//        chain->SetBranchAddress("Muon_dB", Muon_dB );

        chain->SetBranchAddress("Muon_PfChargedHadronIsoR04", &Muon_PfChargedHadronIsoR04);
        chain->SetBranchAddress("Muon_PfNeutralHadronIsoR04", &Muon_PfNeutralHadronIsoR04);
        chain->SetBranchAddress("Muon_PfGammaIsoR04", &Muon_PfGammaIsoR04);
        chain->SetBranchAddress("Muon_PFSumPUIsoR04", &Muon_PFSumPUIsoR04);

        chain->SetBranchAddress("CosAngle", &CosAngle);
        chain->SetBranchAddress("vtxTrkChi2", &vtxTrkChi2);
        chain->SetBranchAddress("vtxTrkProb", &vtxTrkProb);
        chain->SetBranchAddress("vtxTrkNdof", &vtxTrkNdof);
        chain->SetBranchAddress("vtxTrkCkt1Pt", &vtxTrkCkt1Pt);
        chain->SetBranchAddress("vtxTrkCkt2Pt", &vtxTrkCkt2Pt);
//        chain->SetBranchAddress("vtxTrkDiEChi2", &vtxTrkDiEChi2);
//        chain->SetBranchAddress("vtxTrkDiEProb", &vtxTrkDiEProb);
//        chain->SetBranchAddress("vtxTrkDiENdof", &vtxTrkDiENdof);
//        chain->SetBranchAddress("vtxTrkDiE1Pt", &vtxTrkDiE1Pt);
//        chain->SetBranchAddress("vtxTrkDiE2Pt", &vtxTrkDiE2Pt);
//        chain->SetBranchAddress("vtxTrkEMuChi2", &vtxTrkEMuChi2);
//        chain->SetBranchAddress("vtxTrkEMuProb", &vtxTrkEMuProb);
//        chain->SetBranchAddress("vtxTrkEMuNdof", &vtxTrkEMuNdof);
//        chain->SetBranchAddress("vtxTrkEMu1Pt", &vtxTrkEMu1Pt);
//        chain->SetBranchAddress("vtxTrkEMu2Pt", &vtxTrkEMu2Pt);

        chain->SetBranchAddress("Muon_Best_pT", &Muon_Best_pT);
        chain->SetBranchAddress("Muon_Best_pTError", &Muon_Best_pTError);
        chain->SetBranchAddress("Muon_Best_Px", &Muon_Best_Px);
        chain->SetBranchAddress("Muon_Best_Py", &Muon_Best_Py);
        chain->SetBranchAddress("Muon_Best_Pz", &Muon_Best_Pz);
        chain->SetBranchAddress("Muon_Best_eta", &Muon_Best_eta);
        chain->SetBranchAddress("Muon_Best_phi", &Muon_Best_phi);

        chain->SetBranchAddress("Muon_Inner_pT", &Muon_Inner_pT);
        chain->SetBranchAddress("Muon_Inner_pTError", &Muon_Inner_pTError);
        chain->SetBranchAddress("Muon_Inner_Px", &Muon_Inner_Px);
        chain->SetBranchAddress("Muon_Inner_Py", &Muon_Inner_Py);
        chain->SetBranchAddress("Muon_Inner_Pz", &Muon_Inner_Pz);
        chain->SetBranchAddress("Muon_Inner_eta", &Muon_Inner_eta);
        chain->SetBranchAddress("Muon_Inner_phi", &Muon_Inner_phi);

        chain->SetBranchAddress("Muon_Outer_pT", &Muon_Outer_pT);
        chain->SetBranchAddress("Muon_Outer_pTError", &Muon_Outer_pTError);
        chain->SetBranchAddress("Muon_Outer_Px", &Muon_Outer_Px);
        chain->SetBranchAddress("Muon_Outer_Py", &Muon_Outer_Py);
        chain->SetBranchAddress("Muon_Outer_Pz", &Muon_Outer_Pz);
        chain->SetBranchAddress("Muon_Outer_eta", &Muon_Outer_eta);
        chain->SetBranchAddress("Muon_Outer_phi", &Muon_Outer_phi);

        chain->SetBranchAddress("Muon_GLB_pT", &Muon_GLB_pT);
        chain->SetBranchAddress("Muon_GLB_pTError", &Muon_GLB_pTError);
        chain->SetBranchAddress("Muon_GLB_Px", &Muon_GLB_Px);
        chain->SetBranchAddress("Muon_GLB_Py", &Muon_GLB_Py);
        chain->SetBranchAddress("Muon_GLB_Pz", &Muon_GLB_Pz);
        chain->SetBranchAddress("Muon_GLB_eta", &Muon_GLB_eta);
        chain->SetBranchAddress("Muon_GLB_phi", &Muon_GLB_phi);

        chain->SetBranchAddress("Muon_TuneP_pT", &Muon_TuneP_pT);
        chain->SetBranchAddress("Muon_TuneP_pTError", &Muon_TuneP_pTError);
        chain->SetBranchAddress("Muon_TuneP_Px", &Muon_TuneP_Px);
        chain->SetBranchAddress("Muon_TuneP_Py", &Muon_TuneP_Py);
        chain->SetBranchAddress("Muon_TuneP_Pz", &Muon_TuneP_Pz);
        chain->SetBranchAddress("Muon_TuneP_eta", &Muon_TuneP_eta);
        chain->SetBranchAddress("Muon_TuneP_phi", &Muon_TuneP_phi);

//        chain->SetBranchAddress("Muon_stationMask", &Muon_stationMask);
//        chain->SetBranchAddress("Muon_nMatchesRPCLayers", &Muon_nMatchesRPCLayers);

        // -- Adding to cache -- //
        chain->AddBranchToCache("nVertices", 1);
        chain->AddBranchToCache("runNum", 1);
        chain->AddBranchToCache("lumiBlock", 1);
        chain->AddBranchToCache("evtNum", 1);
        chain->AddBranchToCache("nPileUp", 1);
        chain->AddBranchToCache("isHardProcess", 1);
        chain->AddBranchToCache("GENEvt_weight", 1);

        chain->AddBranchToCache("HLT_trigName", 1);
        chain->AddBranchToCache("HLT_ntrig", 1);
        chain->AddBranchToCache("HLT_trigFired", 1);
//        chain->AddBranchToCache("HLT_trigPt", 1);
        chain->AddBranchToCache("HLT_trigEta", 1);
        chain->AddBranchToCache("HLT_trigPhi", 1);

        chain->AddBranchToCache("Muon_InvM", 1);

        chain->AddBranchToCache("isPFmuon", 1);
        chain->AddBranchToCache("isGLBmuon", 1);
        chain->AddBranchToCache("isTRKmuon", 1);

        chain->AddBranchToCache("CosAngle", 1);
        chain->AddBranchToCache("vtxTrkChi2", 1);
        chain->AddBranchToCache("vtxTrkProb", 1);
        chain->AddBranchToCache("vtxTrkNdof", 1);
        chain->AddBranchToCache("vtxTrkCkt1Pt", 1);
        chain->AddBranchToCache("vtxTrkCkt2Pt", 1);
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

        chain->AddBranchToCache("Muon_pT", 1);
        chain->AddBranchToCache("Muon_eta", 1);
        chain->AddBranchToCache("Muon_phi", 1);
        chain->AddBranchToCache("Muon_Energy", 1);
        chain->AddBranchToCache("Muon_charge", 1);
//        chain->AddBranchToCache("Muon_muonType", 1);
        chain->AddBranchToCache("Muon_chi2dof", 1);
        chain->AddBranchToCache("Muon_muonHits", 1);
        chain->AddBranchToCache("Muon_nSegments", 1);
        chain->AddBranchToCache("Muon_nMatches", 1);
        chain->AddBranchToCache("Muon_trackerLayers", 1);

//        chain->AddBranchToCache("Muon_pixelHitsGLB", 1);
//        chain->AddBranchToCache("Muon_trackerLayersGLB", 1);

        chain->AddBranchToCache("Muon_pixelHits", 1);
        chain->AddBranchToCache("Muon_dxyVTX", 1);
        chain->AddBranchToCache("Muon_dzVTX", 1);
        chain->AddBranchToCache("Muon_trkiso", 1);

        chain->AddBranchToCache("Muon_Px", 1);
        chain->AddBranchToCache("Muon_Py", 1);
        chain->AddBranchToCache("Muon_Pz", 1);

//        chain->AddBranchToCache("Muon_dB", 1);

        chain->AddBranchToCache("Muon_PfChargedHadronIsoR04", 1);
        chain->AddBranchToCache("Muon_PfNeutralHadronIsoR04", 1);
        chain->AddBranchToCache("Muon_PfGammaIsoR04", 1);
        chain->AddBranchToCache("Muon_PFSumPUIsoR04", 1);

        chain->AddBranchToCache("Muon_Best_pT", 1);
        chain->AddBranchToCache("Muon_Best_pTError", 1);
        chain->AddBranchToCache("Muon_Best_Px", 1);
        chain->AddBranchToCache("Muon_Best_Py", 1);
        chain->AddBranchToCache("Muon_Best_Pz", 1);
        chain->AddBranchToCache("Muon_Best_eta", 1);
        chain->AddBranchToCache("Muon_Best_phi", 1);

        chain->AddBranchToCache("Muon_Inner_pT", 1);
        chain->AddBranchToCache("Muon_Inner_pTError", 1);
        chain->AddBranchToCache("Muon_Inner_eta", 1);
        chain->AddBranchToCache("Muon_Inner_phi", 1);
        chain->AddBranchToCache("Muon_Inner_Px", 1);
        chain->AddBranchToCache("Muon_Inner_Py", 1);
        chain->AddBranchToCache("Muon_Inner_Pz", 1);

        chain->AddBranchToCache("Muon_Outer_pT", 1);
        chain->AddBranchToCache("Muon_Outer_pTError", 1);
        chain->AddBranchToCache("Muon_Outer_Px", 1);
        chain->AddBranchToCache("Muon_Outer_Py", 1);
        chain->AddBranchToCache("Muon_Outer_Pz", 1);
        chain->AddBranchToCache("Muon_Outer_eta", 1);
        chain->AddBranchToCache("Muon_Outer_phi", 1);

        chain->AddBranchToCache("Muon_GLB_pT", 1);
        chain->AddBranchToCache("Muon_GLB_pTError", 1);
        chain->AddBranchToCache("Muon_GLB_Px", 1);
        chain->AddBranchToCache("Muon_GLB_Py", 1);
        chain->AddBranchToCache("Muon_GLB_Pz", 1);
        chain->AddBranchToCache("Muon_GLB_eta", 1);
        chain->AddBranchToCache("Muon_GLB_phi", 1);

        chain->AddBranchToCache("Muon_TuneP_pT", 1);
        chain->AddBranchToCache("Muon_TuneP_pTError", 1);
        chain->AddBranchToCache("Muon_TuneP_eta", 1);
        chain->AddBranchToCache("Muon_TuneP_phi", 1);
        chain->AddBranchToCache("Muon_TuneP_Px", 1);
        chain->AddBranchToCache("Muon_TuneP_Py", 1);
        chain->AddBranchToCache("Muon_TuneP_Pz", 1);

//        chain->AddBranchToCache("Muon_stationMask", 1);
//        chain->AddBranchToCache("Muon_nMatchesRPCLayers", 1);

        File_Given = kTRUE;
    }

    void Ready()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->StopCacheLearningPhase();
    }

    void TurnOnAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("nVertices", 1);
        chain->SetBranchStatus("runNum", 1);
        chain->SetBranchStatus("lumiBlock", 1);
        chain->SetBranchStatus("evtNum", 1);
        chain->SetBranchStatus("nPileUp", 1);
        chain->SetBranchStatus("isHardProcess", 1);
        chain->SetBranchStatus("GENEvt_weight", 1);

        chain->SetBranchStatus("HLT_trigName", 1);
        chain->SetBranchStatus("HLT_ntrig", 1);
        chain->SetBranchStatus("HLT_trigFired", 1);
//        chain->SetBranchStatus("HLT_trigPt", 1);
        chain->SetBranchStatus("HLT_trigEta", 1);
        chain->SetBranchStatus("HLT_trigPhi", 1);

        chain->SetBranchStatus("Muon_InvM", 1);

        chain->SetBranchStatus("isPFmuon", 1);
        chain->SetBranchStatus("isGLBmuon", 1);
        chain->SetBranchStatus("isTRKmuon", 1);

        chain->SetBranchStatus("CosAngle", 1);
        chain->SetBranchStatus("vtxTrkChi2", 1);
        chain->SetBranchStatus("vtxTrkProb", 1);
        chain->SetBranchStatus("vtxTrkNdof", 1);
        chain->SetBranchStatus("vtxTrkCkt1Pt", 1);
        chain->SetBranchStatus("vtxTrkCkt2Pt", 1);
//        chain->SetBranchStatus("vtxTrkDiEChi2", 1);
//        chain->SetBranchStatus("vtxTrkDiEProb", 1);
//        chain->SetBranchStatus("vtxTrkDiENdof", 1);
//        chain->SetBranchStatus("vtxTrkDiE1Pt", 1);
//        chain->SetBranchStatus("vtxTrkDiE2Pt", 1);
//        chain->SetBranchStatus("vtxTrkEMuChi2", 1);
//        chain->SetBranchStatus("vtxTrkEMuProb", 1);
//        chain->SetBranchStatus("vtxTrkEMuNdof", 1);
//        chain->SetBranchStatus("vtxTrkEMu1Pt", 1);
//        chain->SetBranchStatus("vtxTrkEMu2Pt", 1);

        chain->SetBranchStatus("Muon_pT", 1);
        chain->SetBranchStatus("Muon_eta", 1);
        chain->SetBranchStatus("Muon_phi", 1);
        chain->SetBranchStatus("Muon_Energy", 1);
        chain->SetBranchStatus("Muon_charge", 1);
//        chain->SetBranchStatus("Muon_muonType", 1);
        chain->SetBranchStatus("Muon_chi2dof", 1);
        chain->SetBranchStatus("Muon_muonHits", 1);
        chain->SetBranchStatus("Muon_nSegments", 1);
        chain->SetBranchStatus("Muon_nMatches", 1);
        chain->SetBranchStatus("Muon_trackerLayers", 1);

//        chain->SetBranchStatus("Muon_pixelHitsGLB", 1);
//        chain->SetBranchStatus("Muon_trackerLayersGLB", 1);

        chain->SetBranchStatus("Muon_pixelHits", 1);
        chain->SetBranchStatus("Muon_dxyVTX", 1);
        chain->SetBranchStatus("Muon_dzVTX", 1);
        chain->SetBranchStatus("Muon_trkiso", 1);

        chain->SetBranchStatus("Muon_Px", 1);
        chain->SetBranchStatus("Muon_Py", 1);
        chain->SetBranchStatus("Muon_Pz", 1);
//        chain->SetBranchStatus("Muon_dB", 1);

        chain->SetBranchStatus("Muon_PfChargedHadronIsoR04", 1);
        chain->SetBranchStatus("Muon_PfNeutralHadronIsoR04" ,1);
        chain->SetBranchStatus("Muon_PfGammaIsoR04", 1);
        chain->SetBranchStatus("Muon_PFSumPUIsoR04", 1);

        chain->SetBranchStatus("Muon_Best_pT", 1);
        chain->SetBranchStatus("Muon_Best_pTError", 1);
        chain->SetBranchStatus("Muon_Best_Px", 1);
        chain->SetBranchStatus("Muon_Best_Py", 1);
        chain->SetBranchStatus("Muon_Best_Pz", 1);
        chain->SetBranchStatus("Muon_Best_eta", 1);
        chain->SetBranchStatus("Muon_Best_phi", 1);

        chain->SetBranchStatus("Muon_Inner_pT", 1);
        chain->SetBranchStatus("Muon_Inner_pTError", 1);
        chain->SetBranchStatus("Muon_Inner_eta", 1);
        chain->SetBranchStatus("Muon_Inner_phi", 1);
        chain->SetBranchStatus("Muon_Inner_Px", 1);
        chain->SetBranchStatus("Muon_Inner_Py", 1);
        chain->SetBranchStatus("Muon_Inner_Pz", 1);

        chain->SetBranchStatus("Muon_Outer_pT", 1);
        chain->SetBranchStatus("Muon_Outer_pTError", 1);
        chain->SetBranchStatus("Muon_Outer_Px", 1);
        chain->SetBranchStatus("Muon_Outer_Py", 1);
        chain->SetBranchStatus("Muon_Outer_Pz", 1);
        chain->SetBranchStatus("Muon_Outer_eta", 1);
        chain->SetBranchStatus("Muon_Outer_phi", 1);

        chain->SetBranchStatus("Muon_GLB_pT", 1);
        chain->SetBranchStatus("Muon_GLB_pTError", 1);
        chain->SetBranchStatus("Muon_GLB_Px", 1);
        chain->SetBranchStatus("Muon_GLB_Py", 1);
        chain->SetBranchStatus("Muon_GLB_Pz", 1);
        chain->SetBranchStatus("Muon_GLB_eta", 1);
        chain->SetBranchStatus("Muon_GLB_phi", 1);

        chain->SetBranchStatus("Muon_TuneP_pT", 1);
        chain->SetBranchStatus("Muon_TuneP_pTError", 1);
        chain->SetBranchStatus("Muon_TuneP_eta", 1);
        chain->SetBranchStatus("Muon_TuneP_phi", 1);
        chain->SetBranchStatus("Muon_TuneP_Px", 1);
        chain->SetBranchStatus("Muon_TuneP_Py", 1);
        chain->SetBranchStatus("Muon_TuneP_Pz", 1);

//        chain->SetBranchStatus("Muon_stationMask", 1);
//        chain->SetBranchStatus("Muon_nMatchesRPCLayers", 1);
    }

    void TurnOffAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("*", 0);
    }

    void TurnOnBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0)
        {
            std::cout << "Error: no branches to turn off" << endl;
            return;
        }
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof())
        {
            ss >> br;
            if (chain->GetBranchStatus(br)==kFALSE)
            {
                chain->SetBranchStatus(br,1);
                std::cout << "Activating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already active\n";
        }
    }

    void TurnOffBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0)
        {
            std::cout << "Error: no branches to turn off" << endl;
            return;
        }
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof())
        {
            ss >> br;
            if (chain->GetBranchStatus(br)==kTRUE)
            {
                chain->SetBranchStatus(br,0);
                std::cout << "Deactivating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already inactive\n";
        }
    }

    void GetEvent(Int_t i)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if(!chain) return;

        chain->GetEntry(i);
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
                    if( HLT_trigFired->at(k) == 1 )
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
                    if( HLT_trigFired->at(k) == 1 )
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
                    if( HLT_trigFired->at(k) == 1 )
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
                    if( HLT_trigFired->at(k) == 1 )
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
                    if( HLT_trigFired->at(k) == 1 )
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

///////////////////////////////////////////////////////////////////////////////////
///                                                                             ///
///     The SelectedMuMu_t class                                                ///
///                                                                             ///
///     Stores vectors (mostly) with only the most important information        ///
///     about a pair of muons that passed the DY->MuMu selection.               ///
///                                                                             ///
///     How to use (examples):                                                  ///
///                                                                             ///
///         To read from file:                                                  ///
///             > #include "SelectedX.h"                                        ///
///             > TChain *ch = new TChain("...");                               ///
///             > ch->Add("...");                                               ///
///             > SelectedMuMu_t mu;                                            ///
///             > mu.CreateFromChain(ch);                                       ///
///             > mu.GetEvent(...);                                             ///
///             > cout << mu.Muon_pT->at(...);                                  ///
///                                                                             ///
///         To create an empty class (e.g. to fill with values and write        ///
///         into a tree):                                                       ///
///             > #include "SelectedX.h"                                        ///
///             > SelectedMuMu_t mu;                                            ///
///             > mu.CreateNew();                                               ///
///             > mu.Muon_pT->push_back(...);                                   ///
///                                                                             ///
///////////////////////////////////////////////////////////////////////////////////
class SelectedMuMu_t : public SelectedX
{
public:
    Bool_t isSelPassed;

    Double_t Muon_InvM;
    std::vector<double> *Muon_pT;
    std::vector<double> *Muon_eta;
    std::vector<double> *Muon_phi;
    std::vector<double> *Muon_Energy;
    std::vector<int> *Muon_charge;
    std::vector<double> *Muon_TuneP_pT;
    std::vector<double> *Muon_TuneP_eta;
    std::vector<double> *Muon_TuneP_phi;
    std::vector<int> *Muon_trackerLayers;

    // -- Constructor -- //
    void CreateNew()
    {
        Muon_pT = new std::vector<double>;
        Muon_eta = new std::vector<double>;
        Muon_phi = new std::vector<double>;
        Muon_Energy = new std::vector<double>;
        Muon_charge = new std::vector<int>;

        Muon_TuneP_pT = new std::vector<double>;
        Muon_TuneP_eta = new std::vector<double>;
        Muon_TuneP_phi = new std::vector<double>;

        Muon_trackerLayers = new std::vector<int>;
    }

    // -- Constructor with chain -- //
    void CreateFromChain(TChain *chainptr)
    {
        chain = chainptr;
//        Int_t cachesize = 500000000; // -- 500MB -- //
//        chain->SetCacheSize(cachesize);

        // -- Setting addresses -- //
        chain->SetBranchAddress("nVertices", &nVertices);
        chain->SetBranchAddress("GENEvt_weight", &GENEvt_weight);
        chain->SetBranchAddress("nPileUp", &nPileUp);
        chain->SetBranchAddress("isSelPassed", &isSelPassed);
        chain->SetBranchAddress("Muon_InvM", &Muon_InvM);
        chain->SetBranchAddress("Muon_pT", &Muon_pT);
        chain->SetBranchAddress("Muon_eta", &Muon_eta);
        chain->SetBranchAddress("Muon_phi", &Muon_phi);
        chain->SetBranchAddress("Muon_Energy", &Muon_Energy);
        chain->SetBranchAddress("Muon_charge", &Muon_charge);
        chain->SetBranchAddress("Muon_TuneP_pT", &Muon_TuneP_pT);
        chain->SetBranchAddress("Muon_TuneP_eta", &Muon_TuneP_eta);
        chain->SetBranchAddress("Muon_TuneP_phi", &Muon_TuneP_phi);
        chain->SetBranchAddress("Muon_trackerLayers", &Muon_trackerLayers);

        // -- Adding to cache -- //
//        chain->AddBranchToCatche("nVertices", 1);
//        chain->AddBranchToCache("GENEvt_weight", 1);
//        chain->AddBranchToCache("nPileUp", 1);
//        chain->AddBranchToCache("isSelPassed", 1);
//        chain->AddBranchToCache("Muon_InvM", 1);
//        chain->AddBranchToCache("Muon_pT", 1);
//        chain->AddBranchToCache("Muon_eta", 1);
//        chain->AddBranchToCache("Muon_phi", 1);
//        chain->AddBranchToCache("Muon_Energy", 1);
//        chain->AddBranchToCache("Muon_charge", 1);
//        chain->AddBranchToCache("Muon_TuneP_pT", 1);
//        chain->AddBranchToCache("Muon_TuneP_eta", 1);
//        chain->AddBranchToCache("Muon_TuneP_phi", 1);
//        chain->AddBranchToCache("Muon_trackerLayers", 1);

        File_Given = kTRUE;
    }

    void CreateFromLongSelectedMuMu(LongSelectedMuMu_t *Mu, Bool_t SelPassed)
    {
        GENEvt_weight = Mu->GENEvt_weight;
        nVertices = Mu->nVertices;
        nPileUp = Mu->nPileUp;
        isSelPassed = SelPassed;
        Muon_InvM = Mu->Muon_InvM;

        Muon_pT = Mu->Muon_pT;
        Muon_eta = Mu->Muon_eta;
        Muon_phi = Mu->Muon_phi;
        Muon_Energy = Mu->Muon_Energy;
        Muon_charge = Mu->Muon_charge;
        Muon_TuneP_pT = Mu->Muon_TuneP_pT;
        Muon_TuneP_eta = Mu->Muon_TuneP_eta;
        Muon_TuneP_phi = Mu->Muon_TuneP_phi;
        Muon_trackerLayers = Mu->Muon_trackerLayers;
    }

    void Ready()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->StopCacheLearningPhase();
    }

    void TurnOnAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("nVertices", 1);
        chain->SetBranchStatus("GENEvt_weight", 1);
        chain->SetBranchStatus("nPileUp", 1);
        chain->SetBranchStatus("isSelPassed", 1);
        chain->SetBranchStatus("Muon_InvM", 1);
        chain->SetBranchStatus("Muon_pT", 1);
        chain->SetBranchStatus("Muon_eta", 1);
        chain->SetBranchStatus("Muon_phi", 1);
        chain->SetBranchStatus("Muon_Energy", 1);
        chain->SetBranchStatus("Muon_charge", 1);
        chain->SetBranchStatus("Muon_TuneP_pT", 1);
        chain->SetBranchStatus("Muon_TuneP_eta", 1);
        chain->SetBranchStatus("Muon_TuneP_phi", 1);
        chain->SetBranchStatus("Muon_trackerLayers", 1);
    }

    void TurnOffAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("*", 0);
    }

    void TurnOnBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0)
        {
            std::cout << "Error: no branches to turn off" << endl;
            return;
        }
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof())
        {
            ss >> br;
            if (chain->GetBranchStatus(br)==kFALSE)
            {
                chain->SetBranchStatus(br,1);
                std::cout << "Activating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already active\n";
        }
    }

    void TurnOffBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0)
        {
            std::cout << "Error: no branches to turn off" << endl;
            return;
        }
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof())
        {
            ss >> br;
            if (chain->GetBranchStatus(br)==kTRUE)
            {
                chain->SetBranchStatus(br,0);
                std::cout << "Deactivating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already inactive\n";
        }
    }

    void GetEvent(Int_t i)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if(!chain) return;

        chain->GetEntry(i);
    }

};

///////////////////////////////////////////////////////////////////////////////////
///                                                                             ///
///     The LongSelectedEE_t class                                              ///
///                                                                             ///
///     Stores information about a pair of electrons that passed the DY->ee     ///
///     selection. It is possible to repeat the full event selection using      ///
///     this class instead of regular ntuple (except for the IsHardProcess).    ///
///                                                                             ///
///     How to use:                                                             ///
///                                                                             ///
///         To read from file:                                                  ///
///             > #include <SelectedX.h>                                        ///
///             > TChain *ch = new TChain("...");                               ///
///             > ch->Add("...");                                               ///
///             > LongSelectedEE_t ele;                                         ///
///             > ele.CreateFromChain(ch);                                      ///
///             > ele.GetEvent(...);                                            ///
///             > cout << ele.Electron_Px->at(...);                             ///
///                                                                             ///
///         To create an empty class (e.g. to fill with values and write        ///
///         into a tree):                                                       ///
///             > #include <SelectedX.h>                                        ///
///             > LongSelectedEE_t ele;                                         ///
///             > ele.CreateNew();                                              ///
///             > ele.Electron_Py->push_back(...);                              ///
///                                                                             ///
///////////////////////////////////////////////////////////////////////////////////
class LongSelectedEE_t : public SelectedX
{
public:
    //Event Information
    Int_t runNum;
    Int_t lumiBlock;
    Int_t evtNum;

    //Trigger variables
    Int_t HLT_ntrig;
    std::vector<int> *HLT_trigFired;
    std::vector<string> *HLT_trigName;
//    std::vector<double> *HLT_trigPt;
    std::vector<double> *HLT_trigEta;
    std::vector<double> *HLT_trigPhi;

    Bool_t isHardProcess;

    Double_t Electron_InvM;

    // -- Electron Variables -- //  
    std::vector<double> *Electron_pT;
    std::vector<double> *Electron_eta;
    std::vector<double> *Electron_phi;
    std::vector<double> *Electron_Energy;
    std::vector<int> *Electron_charge;

    std::vector<double> *Electron_gsfpT;
    std::vector<double> *Electron_gsfPx;
    std::vector<double> *Electron_gsfPy;
    std::vector<double> *Electron_gsfPz;
    std::vector<double> *Electron_gsfEta;
    std::vector<double> *Electron_gsfPhi;
    std::vector<double> *Electron_gsfCharge;

    std::vector<double> *Electron_etaSC;
    std::vector<double> *Electron_phiSC;

    std::vector<double> *Electron_etaWidth;
    std::vector<double> *Electron_phiWidth;
    std::vector<double> *Electron_dEtaIn;
    std::vector<double> *Electron_dEtaInSeed;
    std::vector<double> *Electron_dPhiIn;
    std::vector<double> *Electron_sigmaIEtaIEta;
    std::vector<double> *Electron_Full5x5_SigmaIEtaIEta;
    std::vector<double> *Electron_HoverE;
    std::vector<double> *Electron_fbrem;
    std::vector<double> *Electron_eOverP;
    std::vector<double> *Electron_InvEminusInvP;
    std::vector<double> *Electron_dxyVTX;
    std::vector<double> *Electron_dzVTX;
    std::vector<double> *Electron_dxy;
    std::vector<double> *Electron_dz;
    std::vector<double> *Electron_dxyBS;
    std::vector<double> *Electron_dzBS;
    std::vector<double> *Electron_chIso03;
    std::vector<double> *Electron_nhIso03;
    std::vector<double> *Electron_phIso03;
    std::vector<double> *Electron_ChIso03FromPU;
    std::vector<int> *Electron_mHits;
    std::vector<double> *Electron_EnergySC;
    std::vector<double> *Electron_preEnergySC;
    std::vector<double> *Electron_rawEnergySC;
    std::vector<double> *Electron_etSC;
    std::vector<double> *Electron_E15;
    std::vector<double> *Electron_E25;
    std::vector<double> *Electron_E55;
    std::vector<double> *Electron_RelPFIso_dBeta;
    std::vector<double> *Electron_RelPFIso_Rho;
    std::vector<double> *Electron_r9;
    std::vector<double> *Electron_ecalDriven;
    std::vector<double> *Electron_passConvVeto;

//    std::vector<bool> *Electron_passLooseID;
    std::vector<bool> *Electron_passMediumID;
//    std::vector<bool> *Electron_passTightID;
//    std::vector<bool> *Electron_passMVAID_WP80;
//    std::vector<bool> *Electron_passMVAID_WP90;
//    std::vector<bool> *Electron_passHEEPID;

    // -- Default constructor -- //
    void CreateNew()
    {
        HLT_trigFired = new std::vector<int>;
        HLT_trigName = new std::vector<string>;
//        HLT_trigPt = new std::vector<double>;
        HLT_trigEta = new std::vector<double>;
        HLT_trigPhi = new std::vector<double>;

        Electron_pT = new std::vector<double>;
        Electron_eta = new std::vector<double>;
        Electron_phi = new std::vector<double>;
        Electron_Energy = new std::vector<double>;
        Electron_charge = new std::vector<int>;

        Electron_gsfpT = new std::vector<double>;
        Electron_gsfPx = new std::vector<double>;
        Electron_gsfPy = new std::vector<double>;
        Electron_gsfPz = new std::vector<double>;
        Electron_gsfEta = new std::vector<double>;
        Electron_gsfPhi = new std::vector<double>;
        Electron_gsfCharge = new std::vector<double>;

        Electron_etaSC = new std::vector<double>;
        Electron_phiSC = new std::vector<double>;

        Electron_etaWidth = new std::vector<double>;
        Electron_phiWidth = new std::vector<double>;
        Electron_dEtaIn = new std::vector<double>;
        Electron_dEtaInSeed = new std::vector<double>;
        Electron_dPhiIn = new std::vector<double>;
        Electron_sigmaIEtaIEta = new std::vector<double>;
        Electron_Full5x5_SigmaIEtaIEta = new std::vector<double>;
        Electron_HoverE = new std::vector<double>;
        Electron_fbrem = new std::vector<double>;
        Electron_eOverP = new std::vector<double>;
        Electron_InvEminusInvP = new std::vector<double>;
        Electron_dxyVTX = new std::vector<double>;
        Electron_dzVTX = new std::vector<double>;
        Electron_dxy = new std::vector<double>;
        Electron_dz = new std::vector<double>;
        Electron_dxyBS = new std::vector<double>;
        Electron_dzBS = new std::vector<double>;
        Electron_chIso03 = new std::vector<double>;
        Electron_nhIso03 = new std::vector<double>;
        Electron_phIso03 = new std::vector<double>;
        Electron_ChIso03FromPU = new std::vector<double>;
        Electron_mHits = new std::vector<int>;
        Electron_EnergySC = new std::vector<double>;
        Electron_preEnergySC = new std::vector<double>;
        Electron_rawEnergySC = new std::vector<double>;
        Electron_etSC = new std::vector<double>;
        Electron_E15 = new std::vector<double>;
        Electron_E25 = new std::vector<double>;
        Electron_E55 = new std::vector<double>;
        Electron_RelPFIso_dBeta = new std::vector<double>;
        Electron_RelPFIso_Rho = new std::vector<double>;
        Electron_r9 = new std::vector<double>;
        Electron_ecalDriven = new std::vector<double>;
        Electron_passConvVeto = new std::vector<double>;

//        Electron_passLooseID = new std::vector<bool>;
        Electron_passMediumID = new std::vector<bool>;
//        Electron_passTightID = new std::vector<bool>;
//        Electron_passMVAID_WP80 = new std::vector<bool>;
//        Electron_passMVAID_WP90 = new std::vector<bool>;
//        Electron_passHEEPID = new std::vector<bool>;
    }

    // -- Constructor with chain -- //
    void CreateFromChain(TChain *chainptr)
    {
        chain = chainptr;
        Int_t cachesize = 500000000; // -- 500MB -- //
        chain->SetCacheSize(cachesize);

        // -- Setting addresses -- //
        chain->SetBranchAddress("nVertices", &nVertices);
        chain->SetBranchAddress("runNum", &runNum);
        chain->SetBranchAddress("lumiBlock", &lumiBlock);
        chain->SetBranchAddress("evtNum", &evtNum);
        chain->SetBranchAddress("nPileUp", &nPileUp);
        chain->SetBranchAddress("HLT_trigName", &HLT_trigName);
        chain->SetBranchAddress("HLT_trigFired", &HLT_trigFired);
        chain->SetBranchAddress("HLT_ntrig", &HLT_ntrig);
//        chain->SetBranchAddress("HLT_trigPt", &HLT_trigPt);
        chain->SetBranchAddress("HLT_trigEta", &HLT_trigEta);
        chain->SetBranchAddress("HLT_trigPhi", &HLT_trigPhi);
        chain->SetBranchAddress("GENEvt_weight", &GENEvt_weight);
        chain->SetBranchAddress("isHardProcess", &isHardProcess);
        chain->SetBranchAddress("Electron_InvM", &Electron_InvM);

        chain->SetBranchAddress("Electron_pT", &Electron_pT);
        chain->SetBranchAddress("Electron_eta", &Electron_eta);
        chain->SetBranchAddress("Electron_phi", &Electron_phi);
        chain->SetBranchAddress("Electron_Energy", &Electron_Energy);
        chain->SetBranchAddress("Electron_charge", &Electron_charge);
        chain->SetBranchAddress("Electron_gsfpT", &Electron_gsfpT);
        chain->SetBranchAddress("Electron_gsfPx", &Electron_gsfPx);
        chain->SetBranchAddress("Electron_gsfPy", &Electron_gsfPy);
        chain->SetBranchAddress("Electron_gsfPz", &Electron_gsfPz);
        chain->SetBranchAddress("Electron_gsfEta", &Electron_gsfEta);
        chain->SetBranchAddress("Electron_gsfPhi", &Electron_gsfPhi);
        chain->SetBranchAddress("Electron_gsfCharge", &Electron_gsfCharge);
        chain->SetBranchAddress("Electron_etaSC", &Electron_etaSC);
        chain->SetBranchAddress("Electron_phiSC", &Electron_phiSC);
        chain->SetBranchAddress("Electron_etaWidth", &Electron_etaWidth);
        chain->SetBranchAddress("Electron_phiWidth", &Electron_phiWidth);
        chain->SetBranchAddress("Electron_dEtaIn", &Electron_dEtaIn);
        chain->SetBranchAddress("Electron_dEtaInSeed", &Electron_dEtaInSeed);
        chain->SetBranchAddress("Electron_dPhiIn", &Electron_dPhiIn);
        chain->SetBranchAddress("Electron_sigmaIEtaIEta", &Electron_sigmaIEtaIEta);
        chain->SetBranchAddress("Electron_Full5x5_SigmaIEtaIEta", &Electron_Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("Electron_HoverE", &Electron_HoverE);
        chain->SetBranchAddress("Electron_fbrem", &Electron_fbrem);
        chain->SetBranchAddress("Electron_eOverP", &Electron_eOverP);
        chain->SetBranchAddress("Electron_InvEminusInvP", &Electron_InvEminusInvP);
        chain->SetBranchAddress("Electron_dxyVTX", &Electron_dxyVTX);
        chain->SetBranchAddress("Electron_dzVTX", &Electron_dzVTX);
        chain->SetBranchAddress("Electron_dxy", &Electron_dxy);
        chain->SetBranchAddress("Electron_dz", &Electron_dz);
        chain->SetBranchAddress("Electron_dxyBS", &Electron_dxyBS);
        chain->SetBranchAddress("Electron_dzBS", &Electron_dzBS);
        chain->SetBranchAddress("Electron_chIso03", &Electron_chIso03);
        chain->SetBranchAddress("Electron_nhIso03", &Electron_nhIso03);
        chain->SetBranchAddress("Electron_phIso03", &Electron_phIso03);
        chain->SetBranchAddress("Electron_ChIso03FromPU", &Electron_ChIso03FromPU);
        chain->SetBranchAddress("Electron_mHits", &Electron_mHits);
        chain->SetBranchAddress("Electron_EnergySC", &Electron_EnergySC);
        chain->SetBranchAddress("Electron_preEnergySC", &Electron_preEnergySC);
        chain->SetBranchAddress("Electron_rawEnergySC", &Electron_rawEnergySC);
        chain->SetBranchAddress("Electron_etSC", &Electron_etSC);
        chain->SetBranchAddress("Electron_E15", &Electron_E15);
        chain->SetBranchAddress("Electron_E25", &Electron_E25);
        chain->SetBranchAddress("Electron_E55", &Electron_E55);
        chain->SetBranchAddress("Electron_RelPFIso_dBeta", &Electron_RelPFIso_dBeta);
        chain->SetBranchAddress("Electron_RelPFIso_Rho", &Electron_RelPFIso_Rho);
        chain->SetBranchAddress("Electron_r9", &Electron_r9);
        chain->SetBranchAddress("Electron_ecalDriven", &Electron_ecalDriven);
        chain->SetBranchAddress("Electron_passConvVeto", &Electron_passConvVeto);
//        chain->SetBranchAddress("Electron_passLooseID", &Electron_passLooseID);
        chain->SetBranchAddress("Electron_passMediumID", &Electron_passMediumID);
//        chain->SetBranchAddress("Electron_passTightID", &Electron_passTightID);
//        chain->SetBranchAddress("Electron_passMVAID_WP80", &Electron_passMVAID_WP80);
//        chain->SetBranchAddress("Electron_passMVAID_WP90", &Electron_passMVAID_WP90);
//        chain->SetBranchAddress("Electron_passHEEPID", &Electron_passHEEPID);

        // -- Adding to cache -- //
        chain->AddBranchToCache("nVertices", 1);
        chain->AddBranchToCache("runNum", 1);
        chain->AddBranchToCache("lumiBlock", 1);
        chain->AddBranchToCache("evtNum", 1);
        chain->AddBranchToCache("nPileUp", 1);
        chain->AddBranchToCache("HLT_trigName", 1);
        chain->AddBranchToCache("HLT_ntrig", 1);
        chain->AddBranchToCache("HLT_trigFired", 1);
//        chain->AddBranchToCache("HLT_trigPt", 1);
        chain->AddBranchToCache("HLT_trigEta", 1);
        chain->AddBranchToCache("HLT_trigPhi", 1);
        chain->AddBranchToCache("GENEvt_weight", 1);
        chain->AddBranchToCache("isHardProcess", 1);
        chain->AddBranchToCache("Electron_InvM", 1);

        chain->AddBranchToCache("Electron_pT", 1);
        chain->AddBranchToCache("Electron_eta", 1);
        chain->AddBranchToCache("Electron_phi", 1);
        chain->AddBranchToCache("Electron_Energy", 1);
        chain->AddBranchToCache("Electron_charge", 1);
        chain->AddBranchToCache("Electron_gsfpT", 1);
        chain->AddBranchToCache("Electron_gsfPx", 1);
        chain->AddBranchToCache("Electron_gsfPy", 1);
        chain->AddBranchToCache("Electron_gsfPz", 1);
        chain->AddBranchToCache("Electron_gsfEta", 1);
        chain->AddBranchToCache("Electron_gsfPhi", 1);
        chain->AddBranchToCache("Electron_gsfCharge", 1);
        chain->AddBranchToCache("Electron_etaSC", 1);
        chain->AddBranchToCache("Electron_phiSC", 1);
        chain->AddBranchToCache("Electron_etaWidth", 1);
        chain->AddBranchToCache("Electron_phiWidth", 1);
        chain->AddBranchToCache("Electron_dEtaIn", 1);
        chain->AddBranchToCache("Electron_dEtaInSeed", 1);
        chain->AddBranchToCache("Electron_dPhiIn", 1);
        chain->AddBranchToCache("Electron_sigmaIEtaIEta", 1);
        chain->AddBranchToCache("Electron_Full5x5_SigmaIEtaIEta", 1);
        chain->AddBranchToCache("Electron_HoverE", 1);
        chain->AddBranchToCache("Electron_fbrem", 1);
        chain->AddBranchToCache("Electron_eOverP", 1);
        chain->AddBranchToCache("Electron_InvEminusInvP", 1);
        chain->AddBranchToCache("Electron_dxyVTX", 1);
        chain->AddBranchToCache("Electron_dzVTX", 1);
        chain->AddBranchToCache("Electron_dxy", 1);
        chain->AddBranchToCache("Electron_dz", 1);
        chain->AddBranchToCache("Electron_dxyBS", 1);
        chain->AddBranchToCache("Electron_dzBS", 1);
        chain->AddBranchToCache("Electron_chIso03", 1);
        chain->AddBranchToCache("Electron_nhIso03", 1);
        chain->AddBranchToCache("Electron_phIso03", 1);
        chain->AddBranchToCache("Electron_ChIso03FromPU", 1);
        chain->AddBranchToCache("Electron_mHits", 1);
        chain->AddBranchToCache("Electron_EnergySC", 1);
        chain->AddBranchToCache("Electron_preEnergySC", 1);
        chain->AddBranchToCache("Electron_rawEnergySC", 1);
        chain->AddBranchToCache("Electron_etSC", 1);
        chain->AddBranchToCache("Electron_E15", 1);
        chain->AddBranchToCache("Electron_E25", 1);
        chain->AddBranchToCache("Electron_E55", 1);
        chain->AddBranchToCache("Electron_RelPFIso_dBeta", 1);
        chain->AddBranchToCache("Electron_RelPFIso_Rho", 1);
        chain->AddBranchToCache("Electron_r9", 1);
        chain->AddBranchToCache("Electron_ecalDriven", 1);
        chain->AddBranchToCache("Electron_passConvVeto", 1);
//        chain->AddBranchToCache("Electron_passLooseID", 1);
        chain->AddBranchToCache("Electron_passMediumID", 1);
//        chain->AddBranchToCache("Electron_passTightID", 1);
//        chain->AddBranchToCache("Electron_passMVAID_WP80", 1);
//        chain->AddBranchToCache("Electron_passMVAID_WP90", 1);
//        chain->AddBranchToCache("Electron_passHEEPID", 1);

        File_Given = kTRUE;
    }

    void Ready()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->StopCacheLearningPhase();
    }

    void TurnOnAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("nVertices", 1);
        chain->SetBranchStatus("runNum", 1);
        chain->SetBranchStatus("lumiBlock", 1);
        chain->SetBranchStatus("evtNum", 1);
        chain->SetBranchStatus("nPileUp", 1);

        chain->SetBranchStatus("HLT_trigName", 1);
        chain->SetBranchStatus("HLT_ntrig", 1);
        chain->SetBranchStatus("HLT_trigFired", 1);
//        chain->SetBranchStatus("HLT_trigPt", 1);
        chain->SetBranchStatus("HLT_trigEta", 1);
        chain->SetBranchStatus("HLT_trigPhi", 1);

        chain->SetBranchStatus("GENEvt_weight", 1);
        chain->SetBranchStatus("isHardProcess", 1);
        chain->SetBranchStatus("Electron_InvM", 1);

        chain->SetBranchStatus("Electron_pT", 1);
        chain->SetBranchStatus("Electron_eta", 1);
        chain->SetBranchStatus("Electron_phi", 1);
        chain->SetBranchStatus("Electron_Energy", 1);
        chain->SetBranchStatus("Electron_charge", 1);
        chain->SetBranchStatus("Electron_gsfpT", 1);
        chain->SetBranchStatus("Electron_gsfPx", 1);
        chain->SetBranchStatus("Electron_gsfPy", 1);
        chain->SetBranchStatus("Electron_gsfPz", 1);
        chain->SetBranchStatus("Electron_gsfEta", 1);
        chain->SetBranchStatus("Electron_gsfPhi", 1);
        chain->SetBranchStatus("Electron_gsfCharge", 1);
        chain->SetBranchStatus("Electron_etaSC", 1);
        chain->SetBranchStatus("Electron_phiSC", 1);
        chain->SetBranchStatus("Electron_etaWidth", 1);
        chain->SetBranchStatus("Electron_phiWidth", 1);
        chain->SetBranchStatus("Electron_dEtaIn", 1);
        chain->SetBranchStatus("Electron_dEtaInSeed", 1);
        chain->SetBranchStatus("Electron_dPhiIn", 1);
        chain->SetBranchStatus("Electron_sigmaIEtaIEta", 1);
        chain->SetBranchStatus("Electron_Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("Electron_HoverE", 1);
        chain->SetBranchStatus("Electron_fbrem", 1);
        chain->SetBranchStatus("Electron_eOverP", 1);
        chain->SetBranchStatus("Electron_InvEminusInvP", 1);
        chain->SetBranchStatus("Electron_dxyVTX", 1);
        chain->SetBranchStatus("Electron_dzVTX", 1);
        chain->SetBranchStatus("Electron_dxy", 1);
        chain->SetBranchStatus("Electron_dz", 1);
        chain->SetBranchStatus("Electron_dxyBS", 1);
        chain->SetBranchStatus("Electron_dzBS", 1);
        chain->SetBranchStatus("Electron_chIso03", 1);
        chain->SetBranchStatus("Electron_nhIso03", 1);
        chain->SetBranchStatus("Electron_phIso03", 1);
        chain->SetBranchStatus("Electron_ChIso03FromPU", 1);

        chain->SetBranchStatus("Electron_mHits", 1);
        chain->SetBranchStatus("Electron_EnergySC", 1);
        chain->SetBranchStatus("Electron_preEnergySC", 1);
        chain->SetBranchStatus("Electron_rawEnergySC", 1);
        chain->SetBranchStatus("Electron_etSC", 1);
        chain->SetBranchStatus("Electron_E15", 1);
        chain->SetBranchStatus("Electron_E25", 1);
        chain->SetBranchStatus("Electron_E55", 1);
        chain->SetBranchStatus("Electron_RelPFIso_dBeta", 1);
        chain->SetBranchStatus("Electron_RelPFIso_Rho", 1);
        chain->SetBranchStatus("Electron_r9", 1);
        chain->SetBranchStatus("Electron_ecalDriven", 1);
        chain->SetBranchStatus("Electron_passConvVeto", 1);

//        chain->SetBranchStatus("Electron_passLooseID", 1);
        chain->SetBranchStatus("Electron_passMediumID", 1);
//        chain->SetBranchStatus("Electron_passTightID", 1);
//        chain->SetBranchStatus("Electron_passMVAID_WP80", 1);
//        chain->SetBranchStatus("Electron_passMVAID_WP90", 1);
//        chain->SetBranchStatus("Electron_passHEEPID", 1);
    }

    void TurnOffAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("*", 0);
    }

    void TurnOnBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0) return;
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof()) {
            ss >> br;
            if (chain->GetBranchStatus(br)==kFALSE)
            {
                chain->SetBranchStatus(br,1);
                std::cout << "Activating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already active\n";
        }
    }

    void TurnOffBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0) return;
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof()) {
            ss >> br;
            if (chain->GetBranchStatus(br)==kTRUE)
            {
                chain->SetBranchStatus(br,0);
                std::cout << "Deactivating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already inactive\n";
        }
    }

    void GetEvent(Int_t i)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->GetEntry(i);
    }

    Bool_t isTriggered(TString HLT)
    {
        Bool_t isTrigger = false;
        for( Int_t k = 0; k < HLT_ntrig; k++ )
        {
            if( (HLT_trigName->at((unsigned int)k)) == HLT )
            {
                if( HLT_trigFired->at(k) == 1 )
                {
                    isTrigger = true;
                    break;
                }
            }
        }
        return isTrigger;
    }
};

///////////////////////////////////////////////////////////////////////////////////
///                                                                             ///
///     The SelectedEE_t class                                                  ///
///                                                                             ///
///     Stores vectors (mostly) with only the most important information        ///
///     abouta pair of electrons that passed the DY->ee selection               ///
///                                                                             ///
///     How to use:                                                             ///
///                                                                             ///
///         To read from file:                                                  ///
///             > #include <SelectedX.h>                                        ///
///             > TChain *ch = new TChain("...");                               ///
///             > ch->Add("...");                                               ///
///             > SelectedEE_t ele;                                             ///
///             > ele.CreateFromChain(ch);                                      ///
///             > ele.GetEvent(...);                                            ///
///             > cout << ele.Electron_pT->at(...);                             ///
///                                                                             ///
///         To create an empty class (e.g. to fill with values and write        ///
///         into a tree):                                                       ///
///             > #include <SelectedX.h>                                        ///
///             > SelectedEE_t ele;                                             ///
///             > ele.CreateNew();                                              ///
///             > ele.Electron_pT->push_back(...);                              ///
///                                                                             ///
///////////////////////////////////////////////////////////////////////////////////
class SelectedEE_t : public SelectedX
{
public:
    Bool_t isSelPassed;

    Double_t Electron_InvM;
    std::vector<double> *Electron_pT;
    std::vector<double> *Electron_eta;
    std::vector<double> *Electron_phi;
    std::vector<double> *Electron_Energy;
    std::vector<int> *Electron_charge;
    std::vector<double> *Electron_etaSC;
    std::vector<double> *Electron_phiSC;

    // -- Default constructor -- //
    void CreateNew()
    {
        Electron_pT = new std::vector<double>;
        Electron_eta = new std::vector<double>;
        Electron_phi = new std::vector<double>;
        Electron_Energy = new std::vector<double>;
        Electron_charge = new std::vector<int>;
        Electron_etaSC = new std::vector<double>;
        Electron_phiSC = new std::vector<double>;
    }

    // -- Constructor with chain -- //
    void CreateFromChain(TChain *chainptr)
    {
        chain = chainptr;
        Int_t cachesize = 500000000; // -- 500MB -- //
        chain->SetCacheSize(cachesize);

        // -- Setting addresses -- //
        chain->SetBranchAddress("nVertices", &nVertices);
        chain->SetBranchAddress("GENEvt_weight", &GENEvt_weight);
        chain->SetBranchAddress("nPileUp", &nPileUp);
        chain->SetBranchAddress("isSelPassed", &isSelPassed);
        chain->SetBranchAddress("Electron_InvM", &Electron_InvM);
        chain->SetBranchAddress("Electron_pT", &Electron_pT);
        chain->SetBranchAddress("Electron_eta", &Electron_eta);
        chain->SetBranchAddress("Electron_phi", &Electron_phi);
        chain->SetBranchAddress("Electron_Energy", &Electron_Energy);
        chain->SetBranchAddress("Electron_charge", &Electron_charge);
        chain->SetBranchAddress("Electron_etaSC", &Electron_etaSC);
        chain->SetBranchAddress("Electron_phiSC", &Electron_phiSC);

        // -- Adding to cache -- //
//        chain->AddBranchToCache("nVertices", 1);
//        chain->AddBranchToCache("GENEvt_weight", 1);
//        chain->AddBranchToCache("nPileUp", 1);
//        chain->AddBranchToCache("isSelPassed", 1);
//        chain->AddBranchToCache("Electron_InvM", 1);
//        chain->AddBranchToCache("Electron_pT", 1);
//        chain->AddBranchToCache("Electron_eta", 1);
//        chain->AddBranchToCache("Electron_phi", 1);
//        chain->AddBranchToCache("Electron_Energy", 1);
//        chain->AddBranchToCache("Electron_charge", 1);
//        chain->AddBranchToCache("Electron_etaSC", 1);
//        chain->AddBranchToCache("Electron_phiSC", 1);

        File_Given = kTRUE;
    }

    void CreateFromLongSelectedEE(LongSelectedEE_t *Ele, Bool_t SelPassed)
    {
        nVertices = Ele->nVertices;
        GENEvt_weight = Ele->GENEvt_weight;
        nPileUp = Ele->nPileUp;
        isSelPassed = SelPassed;
        Electron_InvM = Ele->Electron_InvM;

        Electron_pT = Ele->Electron_pT;
        Electron_eta = Ele->Electron_eta;
        Electron_phi = Ele->Electron_phi;
        Electron_Energy = Ele->Electron_Energy;
        Electron_charge = Ele->Electron_charge;
        Electron_etaSC = Ele->Electron_etaSC;
        Electron_phiSC = Ele->Electron_phiSC;
    }

    void Ready()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->StopCacheLearningPhase();
    }

    void TurnOnAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("nVertices", 1);
        chain->SetBranchStatus("GENEvt_weight", 1);
        chain->SetBranchStatus("nPileUp", 1);
        chain->SetBranchStatus("Electron_InvM", 1);
        chain->SetBranchStatus("Electron_pT", 1);
        chain->SetBranchStatus("Electron_eta", 1);
        chain->SetBranchStatus("Electron_phi", 1);
        chain->SetBranchStatus("Electron_Energy", 1);
        chain->SetBranchStatus("Electron_charge", 1);
        chain->SetBranchStatus("Electron_etaSC", 1);
        chain->SetBranchStatus("Electron_phiSC", 1);
    }

    void TurnOffAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("*", 0);
    }

    void TurnOnBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0) return;
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof()) {
            ss >> br;
            if (chain->GetBranchStatus(br)==kFALSE)
            {
                chain->SetBranchStatus(br,1);
                std::cout << "Activating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already active\n";
        }
    }

    void TurnOffBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0) return;
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof()) {
            ss >> br;
            if (chain->GetBranchStatus(br)==kTRUE)
            {
                chain->SetBranchStatus(br,0);
                std::cout << "Deactivating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already inactive\n";
        }
    }

    void GetEvent(Int_t i)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->GetEntry(i);
    }
};

///////////////////////////////////////////////////////////////////////////////////
///                                                                             ///
///     The LongSelectedEMu_t class                                             ///
///                                                                             ///
///     Stores information about an electron-muon pair that passed the EMu      ///
///     selection. It is possible to repeat the full event selection using      ///
///     this class instead of regular ntuple (except for the IsHardProcess).    ///
///                                                                             ///
///     How to use (examples):                                                  ///
///                                                                             ///
///         To read from file:                                                  ///
///             > #include "SelectedX.h"                                        ///
///             > TChain *ch = new TChain("...");                               ///
///             > ch->Add("...");                                               ///
///             > LongSelectedEMu_t emu;                                        ///
///             > emu.CreateFromChain(ch);                                      ///
///             > emu.GetEvent(...);                                            ///
///             > cout << emu.Muon_Px;                                          ///
///                                                                             ///
///         To create an empty class (e.g. to fill with values and write        ///
///         into a tree):                                                       ///
///             > #include "SelectedX.h"                                        ///
///             > LongSelectedEMu_t emu;                                        ///
///             > emu.Electron_Py=...;                                          ///
///                                                                             ///
///////////////////////////////////////////////////////////////////////////////////
class LongSelectedEMu_t : public SelectedX
{
public:
    //Event Information
    Int_t runNum;
    Int_t lumiBlock;
    Int_t evtNum;

    Bool_t isHardProcess;

    //Trigger variables
    Int_t HLT_ntrig;
    std::vector<int> *HLT_trigFired;
    std::vector<string> *HLT_trigName;
//    std::vector<double> *HLT_trigPt;
    std::vector<double> *HLT_trigEta;
    std::vector<double> *HLT_trigPhi;

    Double_t EMu_InvM;

    // ----------- MUON INFORMATION ------------ //

    // -- Physical Variables -- //
    Double_t Muon_pT;
    Double_t Muon_eta;
    Double_t Muon_phi;
    Double_t Muon_Energy;
    Int_t Muon_charge;

    // -- Cut variables -- //
//    Int_t Muon_muonType;
    Double_t Muon_chi2dof;
    Int_t Muon_muonHits;
    Int_t Muon_nSegments;
    Int_t Muon_nMatches;
    Int_t Muon_trackerLayers;

//    Int_t Muon_trackerHitsGLB;
//    Int_t Muon_pixelHitsGLB;
//    Int_t Muon_trackerLayersGLB;

    Int_t Muon_pixelHits;
    Double_t Muon_dxyVTX;
    Double_t Muon_dzVTX;
    Double_t Muon_trkiso;
    Int_t isGLBmuon;
    Int_t isPFmuon;
    Int_t isTRKmuon;

    // -- for invariant mass calculation -- //
    Double_t Muon_Px;
    Double_t Muon_Py;
    Double_t Muon_Pz;

//    Double_t Muon_dB;

    //PF information
    Double_t Muon_PfChargedHadronIsoR04;
    Double_t Muon_PfNeutralHadronIsoR04;
    Double_t Muon_PfGammaIsoR04;
    Double_t Muon_PFSumPUIsoR04;

    //Dimuon variables
//    std::vector<double> *CosAngle;
//    std::vector<double> *vtxTrkChi2;
//    std::vector<double> *vtxTrkProb;
//    std::vector<double> *vtxTrkNdof;
//    std::vector<double> *vtxTrkCkt1Pt;
//    std::vector<double> *vtxTrkCkt2Pt; // If there is just one muon, probably these don't have much use
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

    // -- Various Track Information -- //
    Double_t Muon_Best_pT;
    Double_t Muon_Best_pTError;
    Double_t Muon_Best_Px;
    Double_t Muon_Best_Py;
    Double_t Muon_Best_Pz;
    Double_t Muon_Best_eta;
    Double_t Muon_Best_phi;

    Double_t Muon_Inner_pT;
    Double_t Muon_Inner_pTError;
    Double_t Muon_Inner_Px;
    Double_t Muon_Inner_Py;
    Double_t Muon_Inner_Pz;
    Double_t Muon_Inner_eta;
    Double_t Muon_Inner_phi;

    Double_t Muon_Outer_pT;
    Double_t Muon_Outer_pTError;
    Double_t Muon_Outer_Px;
    Double_t Muon_Outer_Py;
    Double_t Muon_Outer_Pz;
    Double_t Muon_Outer_eta;
    Double_t Muon_Outer_phi;

    Double_t Muon_GLB_pT;
    Double_t Muon_GLB_pTError;
    Double_t Muon_GLB_Px;
    Double_t Muon_GLB_Py;
    Double_t Muon_GLB_Pz;
    Double_t Muon_GLB_eta;
    Double_t Muon_GLB_phi;

    Double_t Muon_TuneP_pT;
    Double_t Muon_TuneP_pTError;
    Double_t Muon_TuneP_Px;
    Double_t Muon_TuneP_Py;
    Double_t Muon_TuneP_Pz;
    Double_t Muon_TuneP_eta;
    Double_t Muon_TuneP_phi;

//    Int_t Muon_stationMask;
//    Int_t Muon_nMatchesRPCLayers;

    // ------------ ELECTRON INFORMATION ------------- //

    Double_t Electron_pT;
    Double_t Electron_eta;
    Double_t Electron_phi;
    Double_t Electron_Energy;
    Int_t Electron_charge;

    Double_t Electron_gsfpT;
    Double_t Electron_gsfPx;
    Double_t Electron_gsfPy;
    Double_t Electron_gsfPz;
    Double_t Electron_gsfEta;
    Double_t Electron_gsfPhi;
    Double_t Electron_gsfCharge;

    Double_t Electron_etaSC;
    Double_t Electron_phiSC;

    Double_t Electron_etaWidth;
    Double_t Electron_phiWidth;
    Double_t Electron_dEtaIn;
    Double_t Electron_dEtaInSeed;
    Double_t Electron_dPhiIn;
    Double_t Electron_sigmaIEtaIEta;
    Double_t Electron_Full5x5_SigmaIEtaIEta;
    Double_t Electron_HoverE;
    Double_t Electron_fbrem;
    Double_t Electron_eOverP;
    Double_t Electron_InvEminusInvP;
    Double_t Electron_dxyVTX;
    Double_t Electron_dzVTX;
    Double_t Electron_dxy;
    Double_t Electron_dz;
    Double_t Electron_dxyBS;
    Double_t Electron_dzBS;
    Double_t Electron_chIso03;
    Double_t Electron_nhIso03;
    Double_t Electron_phIso03;
    Double_t Electron_ChIso03FromPU;
    Int_t Electron_mHits;
    Double_t Electron_EnergySC;
    Double_t Electron_preEnergySC;
    Double_t Electron_rawEnergySC;
    Double_t Electron_etSC;
    Double_t Electron_E15;
    Double_t Electron_E25;
    Double_t Electron_E55;
    Double_t Electron_RelPFIso_dBeta;
    Double_t Electron_RelPFIso_Rho;
    Double_t Electron_r9;
    Double_t Electron_ecalDriven;
    Double_t Electron_passConvVeto;

//    Bool_t Electron_passLooseID;
    Bool_t Electron_passMediumID;
//    Bool_t Electron_passTightID;
//    Bool_t Electron_passMVAID_WP80;
//    Bool_t Electron_passMVAID_WP90;
//    Bool_t Electron_passHEEPID;

    // -- Constructor -- //
    void CreateNew()
    {
        HLT_trigFired = new std::vector<int>;
        HLT_trigName = new std::vector<string>;
        HLT_trigEta = new std::vector<double>;
        HLT_trigPhi = new std::vector<double>;
//        CosAngle = new std::vector<double>;
//        vtxTrkChi2 = new std::vector<double>;
//        vtxTrkProb = new std::vector<double>;
//        vtxTrkNdof = new std::vector<double>;
//        vtxTrkCkt1Pt = new std::vector<double>;
//        vtxTrkCkt2Pt = new std::vector<double>;
//        vtxTrkDiEChi2 = new std::vector<double>;
//        vtxTrkDiEProb = new std::vector<double>;
//        vtxTrkDiENdof = new std::vector<double>;
//        vtxTrkDiE1Pt = new std::vector<double>;
//        vtxTrkDiE2Pt = new std::vector<double>;
//        vtxTrkEMuChi2 = new std::vector<double>;
//        vtxTrkEMuProb = new std::vector<double>;
//        vtxTrkEMuNdof = new std::vector<double>;
//        vtxTrkEMu1Pt = new std::vector<double>;
//        vtxTrkEMu2Pt = new std::vector<double>;

        return;
    }

    // -- Constructor with chain -- //
    void CreateFromChain(TChain *chainptr)
    {
        chain = chainptr;
        Int_t cachesize = 500000000; // -- 500MB -- //
        chain->SetCacheSize(cachesize);

        // -- Setting addresses -- //
        chain->SetBranchAddress("nVertices", &nVertices);
        chain->SetBranchAddress("runNum", &runNum);
        chain->SetBranchAddress("lumiBlock", &lumiBlock);
        chain->SetBranchAddress("evtNum", &evtNum);
        chain->SetBranchAddress("nPileUp", &nPileUp);
        chain->SetBranchAddress("isHardProcess", &isHardProcess);
        chain->SetBranchAddress("GENEvt_weight", &GENEvt_weight);

        chain->SetBranchAddress("HLT_trigName", &HLT_trigName);
        chain->SetBranchAddress("HLT_trigFired", &HLT_trigFired);
        chain->SetBranchAddress("HLT_ntrig", &HLT_ntrig);
//        chain->SetBranchAddress("HLT_trigPt", &HLT_trigPt);
        chain->SetBranchAddress("HLT_trigEta", &HLT_trigEta);
        chain->SetBranchAddress("HLT_trigPhi", &HLT_trigPhi);

        chain->SetBranchAddress("EMu_InvM", &EMu_InvM);

        // -- Muon information -- //
        chain->SetBranchAddress("Muon_pT", &Muon_pT);
        chain->SetBranchAddress("Muon_eta", &Muon_eta);
        chain->SetBranchAddress("Muon_phi", &Muon_phi);
        chain->SetBranchAddress("Muon_Energy", &Muon_Energy);
        chain->SetBranchAddress("Muon_charge", &Muon_charge);
//        chain->SetBranchAddress("Muon_muonType", &Muon_muonType);
        chain->SetBranchAddress("Muon_chi2dof", &Muon_chi2dof);
        chain->SetBranchAddress("Muon_muonHits", &Muon_muonHits);
        chain->SetBranchAddress("Muon_nSegments", &Muon_nSegments);
        chain->SetBranchAddress("Muon_nMatches", &Muon_nMatches);
        chain->SetBranchAddress("Muon_trackerLayers", &Muon_trackerLayers);
//        chain->SetBranchAddress("Muon_pixelHitsGLB", &Muon_pixelHitsGLB);
//        chain->SetBranchAddress("Muon_trackerLayersGLB", &Muon_trackerLayersGLB);
        chain->SetBranchAddress("Muon_pixelHits", &Muon_pixelHits);
        chain->SetBranchAddress("Muon_dxyVTX", &Muon_dxyVTX);
        chain->SetBranchAddress("Muon_dzVTX", &Muon_dzVTX);
        chain->SetBranchAddress("Muon_trkiso", &Muon_trkiso);
        chain->SetBranchAddress("isPFmuon", &isPFmuon);
        chain->SetBranchAddress("isGLBmuon", &isGLBmuon);
        chain->SetBranchAddress("isTRKmuon", &isTRKmuon);
        chain->SetBranchAddress("Muon_Px", &Muon_Px );
        chain->SetBranchAddress("Muon_Py", &Muon_Py );
        chain->SetBranchAddress("Muon_Pz", &Muon_Pz );
//        chain->SetBranchAddress("Muon_dB", Muon_dB );

        chain->SetBranchAddress("Muon_PfChargedHadronIsoR04", &Muon_PfChargedHadronIsoR04);
        chain->SetBranchAddress("Muon_PfNeutralHadronIsoR04", &Muon_PfNeutralHadronIsoR04);
        chain->SetBranchAddress("Muon_PfGammaIsoR04", &Muon_PfGammaIsoR04);
        chain->SetBranchAddress("Muon_PFSumPUIsoR04", &Muon_PFSumPUIsoR04);

//        chain->SetBranchAddress("CosAngle", &CosAngle);
//        chain->SetBranchAddress("vtxTrkChi2", &vtxTrkChi2);
//        chain->SetBranchAddress("vtxTrkProb", &vtxTrkProb);
//        chain->SetBranchAddress("vtxTrkNdof", &vtxTrkNdof);
//        chain->SetBranchAddress("vtxTrkCkt1Pt", &vtxTrkCkt1Pt);
//        chain->SetBranchAddress("vtxTrkCkt2Pt", &vtxTrkCkt2Pt);
//        chain->SetBranchAddress("vtxTrkDiEChi2", &vtxTrkDiEChi2);
//        chain->SetBranchAddress("vtxTrkDiEProb", &vtxTrkDiEProb);
//        chain->SetBranchAddress("vtxTrkDiENdof", &vtxTrkDiENdof);
//        chain->SetBranchAddress("vtxTrkDiE1Pt", &vtxTrkDiE1Pt);
//        chain->SetBranchAddress("vtxTrkDiE2Pt", &vtxTrkDiE2Pt);
//        chain->SetBranchAddress("vtxTrkEMuChi2", &vtxTrkEMuChi2);
//        chain->SetBranchAddress("vtxTrkEMuProb", &vtxTrkEMuProb);
//        chain->SetBranchAddress("vtxTrkEMuNdof", &vtxTrkEMuNdof);
//        chain->SetBranchAddress("vtxTrkEMu1Pt", &vtxTrkEMu1Pt);
//        chain->SetBranchAddress("vtxTrkEMu2Pt", &vtxTrkEMu2Pt);

        chain->SetBranchAddress("Muon_Best_pT", &Muon_Best_pT);
        chain->SetBranchAddress("Muon_Best_pTError", &Muon_Best_pTError);
        chain->SetBranchAddress("Muon_Best_Px", &Muon_Best_Px);
        chain->SetBranchAddress("Muon_Best_Py", &Muon_Best_Py);
        chain->SetBranchAddress("Muon_Best_Pz", &Muon_Best_Pz);
        chain->SetBranchAddress("Muon_Best_eta", &Muon_Best_eta);
        chain->SetBranchAddress("Muon_Best_phi", &Muon_Best_phi);

        chain->SetBranchAddress("Muon_Inner_pT", &Muon_Inner_pT);
        chain->SetBranchAddress("Muon_Inner_pTError", &Muon_Inner_pTError);
        chain->SetBranchAddress("Muon_Inner_Px", &Muon_Inner_Px);
        chain->SetBranchAddress("Muon_Inner_Py", &Muon_Inner_Py);
        chain->SetBranchAddress("Muon_Inner_Pz", &Muon_Inner_Pz);
        chain->SetBranchAddress("Muon_Inner_eta", &Muon_Inner_eta);
        chain->SetBranchAddress("Muon_Inner_phi", &Muon_Inner_phi);

        chain->SetBranchAddress("Muon_Outer_pT", &Muon_Outer_pT);
        chain->SetBranchAddress("Muon_Outer_pTError", &Muon_Outer_pTError);
        chain->SetBranchAddress("Muon_Outer_Px", &Muon_Outer_Px);
        chain->SetBranchAddress("Muon_Outer_Py", &Muon_Outer_Py);
        chain->SetBranchAddress("Muon_Outer_Pz", &Muon_Outer_Pz);
        chain->SetBranchAddress("Muon_Outer_eta", &Muon_Outer_eta);
        chain->SetBranchAddress("Muon_Outer_phi", &Muon_Outer_phi);

        chain->SetBranchAddress("Muon_GLB_pT", &Muon_GLB_pT);
        chain->SetBranchAddress("Muon_GLB_pTError", &Muon_GLB_pTError);
        chain->SetBranchAddress("Muon_GLB_Px", &Muon_GLB_Px);
        chain->SetBranchAddress("Muon_GLB_Py", &Muon_GLB_Py);
        chain->SetBranchAddress("Muon_GLB_Pz", &Muon_GLB_Pz);
        chain->SetBranchAddress("Muon_GLB_eta", &Muon_GLB_eta);
        chain->SetBranchAddress("Muon_GLB_phi", &Muon_GLB_phi);

        chain->SetBranchAddress("Muon_TuneP_pT", &Muon_TuneP_pT);
        chain->SetBranchAddress("Muon_TuneP_pTError", &Muon_TuneP_pTError);
        chain->SetBranchAddress("Muon_TuneP_Px", &Muon_TuneP_Px);
        chain->SetBranchAddress("Muon_TuneP_Py", &Muon_TuneP_Py);
        chain->SetBranchAddress("Muon_TuneP_Pz", &Muon_TuneP_Pz);
        chain->SetBranchAddress("Muon_TuneP_eta", &Muon_TuneP_eta);
        chain->SetBranchAddress("Muon_TuneP_phi", &Muon_TuneP_phi);

//        chain->SetBranchAddress("Muon_stationMask", &Muon_stationMask);
//        chain->SetBranchAddress("Muon_nMatchesRPCLayers", &Muon_nMatchesRPCLayers);

        // -- Electron information -- //
        chain->SetBranchAddress("Electron_Energy", &Electron_Energy);
        chain->SetBranchAddress("Electron_pT", &Electron_pT);
        chain->SetBranchAddress("Electron_eta", &Electron_eta);
        chain->SetBranchAddress("Electron_phi", &Electron_phi);
        chain->SetBranchAddress("Electron_charge", &Electron_charge);
        chain->SetBranchAddress("Electron_gsfpT", &Electron_gsfpT);
        chain->SetBranchAddress("Electron_gsfPx", &Electron_gsfPx);
        chain->SetBranchAddress("Electron_gsfPy", &Electron_gsfPy);
        chain->SetBranchAddress("Electron_gsfPz", &Electron_gsfPz);
        chain->SetBranchAddress("Electron_gsfEta", &Electron_gsfEta);
        chain->SetBranchAddress("Electron_gsfPhi", &Electron_gsfPhi);
        chain->SetBranchAddress("Electron_gsfCharge", &Electron_gsfCharge);
        chain->SetBranchAddress("Electron_etaSC", &Electron_etaSC);
        chain->SetBranchAddress("Electron_phiSC", &Electron_phiSC);
        chain->SetBranchAddress("Electron_etaWidth", &Electron_etaWidth);
        chain->SetBranchAddress("Electron_phiWidth", &Electron_phiWidth);
        chain->SetBranchAddress("Electron_dEtaIn", &Electron_dEtaIn);
        chain->SetBranchAddress("Electron_dEtaInSeed", &Electron_dEtaInSeed);
        chain->SetBranchAddress("Electron_dPhiIn", &Electron_dPhiIn);
        chain->SetBranchAddress("Electron_sigmaIEtaIEta", &Electron_sigmaIEtaIEta);
        chain->SetBranchAddress("Electron_Full5x5_SigmaIEtaIEta", &Electron_Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("Electron_HoverE", &Electron_HoverE);
        chain->SetBranchAddress("Electron_fbrem", &Electron_fbrem);
        chain->SetBranchAddress("Electron_eOverP", &Electron_eOverP);
        chain->SetBranchAddress("Electron_InvEminusInvP", &Electron_InvEminusInvP);
        chain->SetBranchAddress("Electron_dxyVTX", &Electron_dxyVTX);
        chain->SetBranchAddress("Electron_dzVTX", &Electron_dzVTX);
        chain->SetBranchAddress("Electron_dxy", &Electron_dxy);
        chain->SetBranchAddress("Electron_dz", &Electron_dz);
        chain->SetBranchAddress("Electron_dxyBS", &Electron_dxyBS);
        chain->SetBranchAddress("Electron_dzBS", &Electron_dzBS);
        chain->SetBranchAddress("Electron_chIso03", &Electron_chIso03);
        chain->SetBranchAddress("Electron_nhIso03", &Electron_nhIso03);
        chain->SetBranchAddress("Electron_phIso03", &Electron_phIso03);
        chain->SetBranchAddress("Electron_ChIso03FromPU", &Electron_ChIso03FromPU);
        chain->SetBranchAddress("Electron_mHits", &Electron_mHits);
        chain->SetBranchAddress("Electron_EnergySC", &Electron_EnergySC);
        chain->SetBranchAddress("Electron_preEnergySC", &Electron_preEnergySC);
        chain->SetBranchAddress("Electron_rawEnergySC", &Electron_rawEnergySC);
        chain->SetBranchAddress("Electron_etSC", &Electron_etSC);
        chain->SetBranchAddress("Electron_E15", &Electron_E15);
        chain->SetBranchAddress("Electron_E25", &Electron_E25);
        chain->SetBranchAddress("Electron_E55", &Electron_E55);
        chain->SetBranchAddress("Electron_RelPFIso_dBeta", &Electron_RelPFIso_dBeta);
        chain->SetBranchAddress("Electron_RelPFIso_Rho", &Electron_RelPFIso_Rho);
        chain->SetBranchAddress("Electron_r9", &Electron_r9);
        chain->SetBranchAddress("Electron_ecalDriven", &Electron_ecalDriven);
        chain->SetBranchAddress("Electron_passConvVeto", &Electron_passConvVeto);
//        chain->SetBranchAddress("Electron_passLooseID", &Electron_passLooseID);
        chain->SetBranchAddress("Electron_passMediumID", &Electron_passMediumID);
//        chain->SetBranchAddress("Electron_passTightID", &Electron_passTightID);
//        chain->SetBranchAddress("Electron_passMVAID_WP80", &Electron_passMVAID_WP80);
//        chain->SetBranchAddress("Electron_passMVAID_WP90", &Electron_passMVAID_WP90);
//        chain->SetBranchAddress("Electron_passHEEPID", &Electron_passHEEPID);

        // -- Adding to cache -- //
        chain->AddBranchToCache("nVertices", 1);
        chain->AddBranchToCache("runNum", 1);
        chain->AddBranchToCache("lumiBlock", 1);
        chain->AddBranchToCache("evtNum", 1);
        chain->AddBranchToCache("nPileUp", 1);
        chain->AddBranchToCache("isHardProcess", 1);
        chain->AddBranchToCache("GENEvt_weight", 1);

        chain->AddBranchToCache("HLT_trigName", 1);
        chain->AddBranchToCache("HLT_ntrig", 1);
        chain->AddBranchToCache("HLT_trigFired", 1);
//        chain->AddBranchToCache("HLT_trigPt", 1);
        chain->AddBranchToCache("HLT_trigEta", 1);
        chain->AddBranchToCache("HLT_trigPhi", 1);

        chain->AddBranchToCache("EMu_InvM", 1);

        // -- Muon -- //
        chain->AddBranchToCache("isPFmuon", 1);
        chain->AddBranchToCache("isGLBmuon", 1);
        chain->AddBranchToCache("isTRKmuon", 1);

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

        chain->AddBranchToCache("Muon_pT", 1);
        chain->AddBranchToCache("Muon_eta", 1);
        chain->AddBranchToCache("Muon_phi", 1);
        chain->AddBranchToCache("Muon_Energy", 1);
        chain->AddBranchToCache("Muon_charge", 1);
//        chain->AddBranchToCache("Muon_muonType", 1);
        chain->AddBranchToCache("Muon_chi2dof", 1);
        chain->AddBranchToCache("Muon_muonHits", 1);
        chain->AddBranchToCache("Muon_nSegments", 1);
        chain->AddBranchToCache("Muon_nMatches", 1);
        chain->AddBranchToCache("Muon_trackerLayers", 1);
//        chain->AddBranchToCache("Muon_pixelHitsGLB", 1);
//        chain->AddBranchToCache("Muon_trackerLayersGLB", 1);
        chain->AddBranchToCache("Muon_pixelHits", 1);
        chain->AddBranchToCache("Muon_dxyVTX", 1);
        chain->AddBranchToCache("Muon_dzVTX", 1);
        chain->AddBranchToCache("Muon_trkiso", 1);
        chain->AddBranchToCache("Muon_Px", 1);
        chain->AddBranchToCache("Muon_Py", 1);
        chain->AddBranchToCache("Muon_Pz", 1);
//        chain->AddBranchToCache("Muon_dB", 1);

        chain->AddBranchToCache("Muon_PfChargedHadronIsoR04", 1);
        chain->AddBranchToCache("Muon_PfNeutralHadronIsoR04", 1);
        chain->AddBranchToCache("Muon_PfGammaIsoR04", 1);
        chain->AddBranchToCache("Muon_PFSumPUIsoR04", 1);

        chain->AddBranchToCache("Muon_Best_pT", 1);
        chain->AddBranchToCache("Muon_Best_pTError", 1);
        chain->AddBranchToCache("Muon_Best_Px", 1);
        chain->AddBranchToCache("Muon_Best_Py", 1);
        chain->AddBranchToCache("Muon_Best_Pz", 1);
        chain->AddBranchToCache("Muon_Best_eta", 1);
        chain->AddBranchToCache("Muon_Best_phi", 1);

        chain->AddBranchToCache("Muon_Inner_pT", 1);
        chain->AddBranchToCache("Muon_Inner_pTError", 1);
        chain->AddBranchToCache("Muon_Inner_eta", 1);
        chain->AddBranchToCache("Muon_Inner_phi", 1);
        chain->AddBranchToCache("Muon_Inner_Px", 1);
        chain->AddBranchToCache("Muon_Inner_Py", 1);
        chain->AddBranchToCache("Muon_Inner_Pz", 1);

        chain->AddBranchToCache("Muon_Outer_pT", 1);
        chain->AddBranchToCache("Muon_Outer_pTError", 1);
        chain->AddBranchToCache("Muon_Outer_Px", 1);
        chain->AddBranchToCache("Muon_Outer_Py", 1);
        chain->AddBranchToCache("Muon_Outer_Pz", 1);
        chain->AddBranchToCache("Muon_Outer_eta", 1);
        chain->AddBranchToCache("Muon_Outer_phi", 1);

        chain->AddBranchToCache("Muon_GLB_pT", 1);
        chain->AddBranchToCache("Muon_GLB_pTError", 1);
        chain->AddBranchToCache("Muon_GLB_Px", 1);
        chain->AddBranchToCache("Muon_GLB_Py", 1);
        chain->AddBranchToCache("Muon_GLB_Pz", 1);
        chain->AddBranchToCache("Muon_GLB_eta", 1);
        chain->AddBranchToCache("Muon_GLB_phi", 1);

        chain->AddBranchToCache("Muon_TuneP_pT", 1);
        chain->AddBranchToCache("Muon_TuneP_pTError", 1);
        chain->AddBranchToCache("Muon_TuneP_eta", 1);
        chain->AddBranchToCache("Muon_TuneP_phi", 1);
        chain->AddBranchToCache("Muon_TuneP_Px", 1);
        chain->AddBranchToCache("Muon_TuneP_Py", 1);
        chain->AddBranchToCache("Muon_TuneP_Pz", 1);

//        chain->AddBranchToCache("Muon_stationMask", 1);
//        chain->AddBranchToCache("Muon_nMatchesRPCLayers", 1);

        // -- Electron -- //
        chain->AddBranchToCache("Electron_pT", 1);
        chain->AddBranchToCache("Electron_eta", 1);
        chain->AddBranchToCache("Electron_phi", 1);
        chain->AddBranchToCache("Electron_Energy", 1);
        chain->AddBranchToCache("Electron_charge", 1);
        chain->AddBranchToCache("Electron_gsfpT", 1);
        chain->AddBranchToCache("Electron_gsfPx", 1);
        chain->AddBranchToCache("Electron_gsfPy", 1);
        chain->AddBranchToCache("Electron_gsfPz", 1);
        chain->AddBranchToCache("Electron_gsfEta", 1);
        chain->AddBranchToCache("Electron_gsfPhi", 1);
        chain->AddBranchToCache("Electron_gsfCharge", 1);
        chain->AddBranchToCache("Electron_etaSC", 1);
        chain->AddBranchToCache("Electron_phiSC", 1);
        chain->AddBranchToCache("Electron_etaWidth", 1);
        chain->AddBranchToCache("Electron_phiWidth", 1);
        chain->AddBranchToCache("Electron_dEtaIn", 1);
        chain->AddBranchToCache("Electron_dEtaInSeed", 1);
        chain->AddBranchToCache("Electron_dPhiIn", 1);
        chain->AddBranchToCache("Electron_sigmaIEtaIEta", 1);
        chain->AddBranchToCache("Electron_Full5x5_SigmaIEtaIEta", 1);
        chain->AddBranchToCache("Electron_HoverE", 1);
        chain->AddBranchToCache("Electron_fbrem", 1);
        chain->AddBranchToCache("Electron_eOverP", 1);
        chain->AddBranchToCache("Electron_InvEminusInvP", 1);
        chain->AddBranchToCache("Electron_dxyVTX", 1);
        chain->AddBranchToCache("Electron_dzVTX", 1);
        chain->AddBranchToCache("Electron_dxy", 1);
        chain->AddBranchToCache("Electron_dz", 1);
        chain->AddBranchToCache("Electron_dxyBS", 1);
        chain->AddBranchToCache("Electron_dzBS", 1);
        chain->AddBranchToCache("Electron_chIso03", 1);
        chain->AddBranchToCache("Electron_nhIso03", 1);
        chain->AddBranchToCache("Electron_phIso03", 1);
        chain->AddBranchToCache("Electron_ChIso03FromPU", 1);
        chain->AddBranchToCache("Electron_mHits", 1);
        chain->AddBranchToCache("Electron_EnergySC", 1);
        chain->AddBranchToCache("Electron_preEnergySC", 1);
        chain->AddBranchToCache("Electron_rawEnergySC", 1);
        chain->AddBranchToCache("Electron_etSC", 1);
        chain->AddBranchToCache("Electron_E15", 1);
        chain->AddBranchToCache("Electron_E25", 1);
        chain->AddBranchToCache("Electron_E55", 1);
        chain->AddBranchToCache("Electron_RelPFIso_dBeta", 1);
        chain->AddBranchToCache("Electron_RelPFIso_Rho", 1);
        chain->AddBranchToCache("Electron_r9", 1);
        chain->AddBranchToCache("Electron_ecalDriven", 1);
        chain->AddBranchToCache("Electron_passConvVeto", 1);
//        chain->AddBranchToCache("Electron_passLooseID", 1);
        chain->AddBranchToCache("Electron_passMediumID", 1);
//        chain->AddBranchToCache("Electron_passTightID", 1);
//        chain->AddBranchToCache("Electron_passMVAID_WP80", 1);
//        chain->AddBranchToCache("Electron_passMVAID_WP90", 1);
//        chain->AddBranchToCache("Electron_passHEEPID", 1);

        File_Given = kTRUE;
    }

    void Ready()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->StopCacheLearningPhase();
    }

    void TurnOnAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("nVertices", 1);
        chain->SetBranchStatus("runNum", 1);
        chain->SetBranchStatus("lumiBlock", 1);
        chain->SetBranchStatus("evtNum", 1);
        chain->SetBranchStatus("nPileUp", 1);
        chain->SetBranchStatus("isHardProcess", 1);
        chain->SetBranchStatus("GENEvt_weight", 1);
        chain->SetBranchStatus("HLT_trigName", 1);
        chain->SetBranchStatus("HLT_ntrig", 1);
        chain->SetBranchStatus("HLT_trigFired", 1);
//        chain->SetBranchStatus("HLT_trigPt", 1);
        chain->SetBranchStatus("HLT_trigEta", 1);
        chain->SetBranchStatus("HLT_trigPhi", 1);

        chain->SetBranchStatus("EMu_InvM", 1);

        // -- Muon -- //
        chain->SetBranchStatus("isPFmuon", 1);
        chain->SetBranchStatus("isGLBmuon", 1);
        chain->SetBranchStatus("isTRKmuon", 1);

//        chain->SetBranchStatus("CosAngle", 1);
//        chain->SetBranchStatus("vtxTrkChi2", 1);
//        chain->SetBranchStatus("vtxTrkProb", 1);
//        chain->SetBranchStatus("vtxTrkNdof", 1);
//        chain->SetBranchStatus("vtxTrkCkt1Pt", 1);
//        chain->SetBranchStatus("vtxTrkCkt2Pt", 1);
//        chain->SetBranchStatus("vtxTrkDiEChi2", 1);
//        chain->SetBranchStatus("vtxTrkDiEProb", 1);
//        chain->SetBranchStatus("vtxTrkDiENdof", 1);
//        chain->SetBranchStatus("vtxTrkDiE1Pt", 1);
//        chain->SetBranchStatus("vtxTrkDiE2Pt", 1);
//        chain->SetBranchStatus("vtxTrkEMuChi2", 1);
//        chain->SetBranchStatus("vtxTrkEMuProb", 1);
//        chain->SetBranchStatus("vtxTrkEMuNdof", 1);
//        chain->SetBranchStatus("vtxTrkEMu1Pt", 1);
//        chain->SetBranchStatus("vtxTrkEMu2Pt", 1);

        chain->SetBranchStatus("Muon_pT", 1);
        chain->SetBranchStatus("Muon_eta", 1);
        chain->SetBranchStatus("Muon_phi", 1);
        chain->SetBranchStatus("Muon_Energy", 1);
        chain->SetBranchStatus("Muon_charge", 1);
//        chain->SetBranchStatus("Muon_muonType", 1);
        chain->SetBranchStatus("Muon_chi2dof", 1);
        chain->SetBranchStatus("Muon_muonHits", 1);
        chain->SetBranchStatus("Muon_nSegments", 1);
        chain->SetBranchStatus("Muon_nMatches", 1);
        chain->SetBranchStatus("Muon_trackerLayers", 1);
//        chain->SetBranchStatus("Muon_pixelHitsGLB", 1);
//        chain->SetBranchStatus("Muon_trackerLayersGLB", 1);
        chain->SetBranchStatus("Muon_pixelHits", 1);
        chain->SetBranchStatus("Muon_dxyVTX", 1);
        chain->SetBranchStatus("Muon_dzVTX", 1);
        chain->SetBranchStatus("Muon_trkiso", 1);
        chain->SetBranchStatus("Muon_Px", 1);
        chain->SetBranchStatus("Muon_Py", 1);
        chain->SetBranchStatus("Muon_Pz", 1);
//        chain->SetBranchStatus("Muon_dB", 1);
        chain->SetBranchStatus("Muon_PfChargedHadronIsoR04", 1);
        chain->SetBranchStatus("Muon_PfNeutralHadronIsoR04" ,1);
        chain->SetBranchStatus("Muon_PfGammaIsoR04", 1);
        chain->SetBranchStatus("Muon_PFSumPUIsoR04", 1);
        chain->SetBranchStatus("Muon_Best_pT", 1);
        chain->SetBranchStatus("Muon_Best_pTError", 1);
        chain->SetBranchStatus("Muon_Best_Px", 1);
        chain->SetBranchStatus("Muon_Best_Py", 1);
        chain->SetBranchStatus("Muon_Best_Pz", 1);
        chain->SetBranchStatus("Muon_Best_eta", 1);
        chain->SetBranchStatus("Muon_Best_phi", 1);
        chain->SetBranchStatus("Muon_Inner_pT", 1);
        chain->SetBranchStatus("Muon_Inner_pTError", 1);
        chain->SetBranchStatus("Muon_Inner_eta", 1);
        chain->SetBranchStatus("Muon_Inner_phi", 1);
        chain->SetBranchStatus("Muon_Inner_Px", 1);
        chain->SetBranchStatus("Muon_Inner_Py", 1);
        chain->SetBranchStatus("Muon_Inner_Pz", 1);
        chain->SetBranchStatus("Muon_Outer_pT", 1);
        chain->SetBranchStatus("Muon_Outer_pTError", 1);
        chain->SetBranchStatus("Muon_Outer_Px", 1);
        chain->SetBranchStatus("Muon_Outer_Py", 1);
        chain->SetBranchStatus("Muon_Outer_Pz", 1);
        chain->SetBranchStatus("Muon_Outer_eta", 1);
        chain->SetBranchStatus("Muon_Outer_phi", 1);
        chain->SetBranchStatus("Muon_GLB_pT", 1);
        chain->SetBranchStatus("Muon_GLB_pTError", 1);
        chain->SetBranchStatus("Muon_GLB_Px", 1);
        chain->SetBranchStatus("Muon_GLB_Py", 1);
        chain->SetBranchStatus("Muon_GLB_Pz", 1);
        chain->SetBranchStatus("Muon_GLB_eta", 1);
        chain->SetBranchStatus("Muon_GLB_phi", 1);
        chain->SetBranchStatus("Muon_TuneP_pT", 1);
        chain->SetBranchStatus("Muon_TuneP_pTError", 1);
        chain->SetBranchStatus("Muon_TuneP_eta", 1);
        chain->SetBranchStatus("Muon_TuneP_phi", 1);
        chain->SetBranchStatus("Muon_TuneP_Px", 1);
        chain->SetBranchStatus("Muon_TuneP_Py", 1);
        chain->SetBranchStatus("Muon_TuneP_Pz", 1);
//        chain->SetBranchStatus("Muon_stationMask", 1);
//        chain->SetBranchStatus("Muon_nMatchesRPCLayers", 1);

        // -- Electron -- //
        chain->SetBranchStatus("Electron_pT", 1);
        chain->SetBranchStatus("Electron_eta", 1);
        chain->SetBranchStatus("Electron_phi", 1);
        chain->SetBranchStatus("Electron_Energy", 1);
        chain->SetBranchStatus("Electron_charge", 1);
        chain->SetBranchStatus("Electron_gsfpT", 1);
        chain->SetBranchStatus("Electron_gsfPx", 1);
        chain->SetBranchStatus("Electron_gsfPy", 1);
        chain->SetBranchStatus("Electron_gsfPz", 1);
        chain->SetBranchStatus("Electron_gsfEta", 1);
        chain->SetBranchStatus("Electron_gsfPhi", 1);
        chain->SetBranchStatus("Electron_gsfCharge", 1);
        chain->SetBranchStatus("Electron_etaSC", 1);
        chain->SetBranchStatus("Electron_phiSC", 1);
        chain->SetBranchStatus("Electron_etaWidth", 1);
        chain->SetBranchStatus("Electron_phiWidth", 1);
        chain->SetBranchStatus("Electron_dEtaIn", 1);
        chain->SetBranchStatus("Electron_dEtaInSeed", 1);
        chain->SetBranchStatus("Electron_dPhiIn", 1);
        chain->SetBranchStatus("Electron_sigmaIEtaIEta", 1);
        chain->SetBranchStatus("Electron_Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("Electron_HoverE", 1);
        chain->SetBranchStatus("Electron_fbrem", 1);
        chain->SetBranchStatus("Electron_eOverP", 1);
        chain->SetBranchStatus("Electron_InvEminusInvP", 1);
        chain->SetBranchStatus("Electron_dxyVTX", 1);
        chain->SetBranchStatus("Electron_dzVTX", 1);
        chain->SetBranchStatus("Electron_dxy", 1);
        chain->SetBranchStatus("Electron_dz", 1);
        chain->SetBranchStatus("Electron_dxyBS", 1);
        chain->SetBranchStatus("Electron_dzBS", 1);
        chain->SetBranchStatus("Electron_chIso03", 1);
        chain->SetBranchStatus("Electron_nhIso03", 1);
        chain->SetBranchStatus("Electron_phIso03", 1);
        chain->SetBranchStatus("Electron_ChIso03FromPU", 1);
        chain->SetBranchStatus("Electron_mHits", 1);
        chain->SetBranchStatus("Electron_EnergySC", 1);
        chain->SetBranchStatus("Electron_preEnergySC", 1);
        chain->SetBranchStatus("Electron_rawEnergySC", 1);
        chain->SetBranchStatus("Electron_etSC", 1);
        chain->SetBranchStatus("Electron_E15", 1);
        chain->SetBranchStatus("Electron_E25", 1);
        chain->SetBranchStatus("Electron_E55", 1);
        chain->SetBranchStatus("Electron_RelPFIso_dBeta", 1);
        chain->SetBranchStatus("Electron_RelPFIso_Rho", 1);
        chain->SetBranchStatus("Electron_r9", 1);
        chain->SetBranchStatus("Electron_ecalDriven", 1);
        chain->SetBranchStatus("Electron_passConvVeto", 1);
//        chain->SetBranchStatus("Electron_passLooseID", 1);
        chain->SetBranchStatus("Electron_passMediumID", 1);
//        chain->SetBranchStatus("Electron_passTightID", 1);
//        chain->SetBranchStatus("Electron_passMVAID_WP80", 1);
//        chain->SetBranchStatus("Electron_passMVAID_WP90", 1);
//        chain->SetBranchStatus("Electron_passHEEPID", 1);
    }

    void TurnOffAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("*", 0);
    }

    void TurnOnBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0)
        {
            std::cout << "Error: no branches to turn off" << endl;
            return;
        }
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof())
        {
            ss >> br;
            if (chain->GetBranchStatus(br)==kFALSE)
            {
                chain->SetBranchStatus(br,1);
                std::cout << "Activating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already active\n";
        }
    }

    void TurnOffBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0)
        {
            std::cout << "Error: no branches to turn off" << endl;
            return;
        }
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof())
        {
            ss >> br;
            if (chain->GetBranchStatus(br)==kTRUE)
            {
                chain->SetBranchStatus(br,0);
                std::cout << "Deactivating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already inactive\n";
        }
    }

    void GetEvent(Int_t i)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if(!chain) return;

        chain->GetEntry(i);
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
                    if( HLT_trigFired->at(k) == 1 )
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
                    if( HLT_trigFired->at(k) == 1 )
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
                    if( HLT_trigFired->at(k) == 1 )
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
                    if( HLT_trigFired->at(k) == 1 )
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
                    if( HLT_trigFired->at(k) == 1 )
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

///////////////////////////////////////////////////////////////////////////////////
///                                                                             ///
///     The SelectedEMu_t class                                                 ///
///                                                                             ///
///     Stores only the most important information about an electron-muon       ///
///     pair that passed the EMu selection.                                     ///
///                                                                             ///
///     How to use (examples):                                                  ///
///                                                                             ///
///         To read from file:                                                  ///
///             > #include "SelectedX.h"                                        ///
///             > TChain *ch = new TChain("...");                               ///
///             > ch->Add("...");                                               ///
///             > SelectedEMu_t emu;                                            ///
///             > emu.CreateFromChain(ch);                                      ///
///             > emu.GetEvent(...);                                            ///
///             > cout << emu.Electron_pT;                                      ///
///                                                                             ///
///         To create an empty class (e.g. to fill with values and write        ///
///         into a tree):                                                       ///
///             > #include "SelectedX.h"                                        ///
///             > SelectedEMu_t emu;                                            ///
///             > emu.Muon_pT=...;                                              ///
///                                                                             ///
///////////////////////////////////////////////////////////////////////////////////
class SelectedEMu_t : public SelectedX
{
public:
    Bool_t isSelPassed;

    Double_t EMu_InvM;

    Double_t Electron_pT;
    Double_t Electron_eta;
    Double_t Electron_phi;
    Double_t Electron_Energy;
    Int_t Electron_charge;
    Double_t Electron_etaSC;
    Double_t Electron_phiSC;

    Double_t Muon_pT;
    Double_t Muon_eta;
    Double_t Muon_phi;
    Double_t Muon_Energy;
    Int_t Muon_charge;
    Double_t Muon_TuneP_pT;
    Double_t Muon_TuneP_eta;
    Double_t Muon_TuneP_phi;
    Int_t Muon_trackerLayers;

    // -- Useless constructor -- //
    void CreateNew()
    {
        return;
    }

    // -- Constructor with chain -- //
    void CreateFromChain(TChain *chainptr)
    {
        chain = chainptr;
        Int_t cachesize = 500000000; // -- 500MB -- //
        chain->SetCacheSize(cachesize);

        // -- Setting addresses -- //
        chain->SetBranchAddress("nVertices", &nVertices);
        chain->SetBranchAddress("GENEvt_weight", &GENEvt_weight);
        chain->SetBranchAddress("nPileUp", &nPileUp);
        chain->SetBranchAddress("isSelPassed", &isSelPassed);
        chain->SetBranchAddress("EMu_InvM", &EMu_InvM);
        chain->SetBranchAddress("Electron_pT", &Electron_pT);
        chain->SetBranchAddress("Electron_eta", &Electron_eta);
        chain->SetBranchAddress("Electron_phi", &Electron_phi);
        chain->SetBranchAddress("Electron_Energy", &Electron_Energy);
        chain->SetBranchAddress("Electron_charge", &Electron_charge);
        chain->SetBranchAddress("Electron_etaSC", &Electron_etaSC);
        chain->SetBranchAddress("Electron_phiSC", &Electron_phiSC);
        chain->SetBranchAddress("Muon_pT", &Muon_pT);
        chain->SetBranchAddress("Muon_eta", &Muon_eta);
        chain->SetBranchAddress("Muon_phi", &Muon_phi);
        chain->SetBranchAddress("Muon_charge", &Muon_charge);
        chain->SetBranchAddress("Muon_Energy", &Muon_Energy);
        chain->SetBranchAddress("Muon_TuneP_pT", &Muon_TuneP_pT);
        chain->SetBranchAddress("Muon_TuneP_eta", &Muon_TuneP_eta);
        chain->SetBranchAddress("Muon_TuneP_phi", &Muon_TuneP_phi);
        chain->SetBranchAddress("Muon_trackerLayers", &Muon_trackerLayers);

        // -- Adding to cache -- //
//        chain->AddBranchToCache("nVertices", 1);
//        chain->AddBranchToCache("GENEvt_weight", 1);
//        chain->AddBranchToCache("nPileUp", 1);
//        chain->AddBranchToCache("isSelPassed", 1);
//        chain->AddBranchToCache("EMu_InvM", 1);
//        chain->AddBranchToCache("Electron_pT", 1);
//        chain->AddBranchToCache("Electron_eta", 1);
//        chain->AddBranchToCache("Electron_phi", 1);
//        chain->AddBranchToCache("Electron_Energy", 1);
//        chain->AddBranchToCache("Electron_charge", 1);
//        chain->AddBranchToCache("Electron_etaSC", 1);
//        chain->AddBranchToCache("Electron_phiSC", 1);
//        chain->AddBranchToCache("Muon_pT", 1);
//        chain->AddBranchToCache("Muon_eta", 1);
//        chain->AddBranchToCache("Muon_phi", 1);
//        chain->AddBranchToCache("Muon_charge", 1);
//        chain->AddBranchToCache("Muon_Energy", 1);
//        chain->AddBranchToCache("Muon_TuneP_pT", 1);
//        chain->AddBranchToCache("Muon_TuneP_eta", 1);
//        chain->AddBranchToCache("Muon_TuneP_phi", 1);
//        chain->AddBranchToCache("Muon_trackerLayers", 1);

        File_Given = kTRUE;
    }

    void CreateFromLongSelectedEMu(LongSelectedEMu_t *EMu, Bool_t SelPassed)
    {
        nVertices = EMu->nVertices;
        GENEvt_weight = EMu->GENEvt_weight;
        nPileUp = EMu->nPileUp;
        isSelPassed = SelPassed;
        EMu_InvM = EMu->EMu_InvM;

        Muon_pT = EMu->Muon_pT;
        Muon_eta = EMu->Muon_eta;
        Muon_phi = EMu->Muon_phi;
        Muon_Energy = EMu->Muon_Energy;
        Muon_charge = EMu->Muon_charge;
        Muon_TuneP_pT = EMu->Muon_TuneP_pT;
        Muon_TuneP_eta = EMu->Muon_TuneP_eta;
        Muon_TuneP_phi = EMu->Muon_TuneP_phi;
        Muon_trackerLayers = EMu->Muon_trackerLayers;

        Electron_pT = EMu->Electron_pT;
        Electron_eta = EMu->Electron_eta;
        Electron_phi = EMu->Electron_phi;
        Electron_Energy = EMu->Electron_Energy;
        Electron_charge = EMu->Electron_charge;
        Electron_etaSC = EMu->Electron_etaSC;
        Electron_phiSC = EMu->Electron_phiSC;
    }

    void Ready()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->StopCacheLearningPhase();
    }

    void TurnOnAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("nVertices", 1);
        chain->SetBranchStatus("GENEvt_weight", 1);
        chain->SetBranchStatus("nPileUp", 1);
        chain->SetBranchStatus("isSelPassed", 1);
        chain->SetBranchStatus("EMu_InvM", 1);
        chain->SetBranchStatus("Electron_pT", 1);
        chain->SetBranchStatus("Electron_eta", 1);
        chain->SetBranchStatus("Electron_phi", 1);
        chain->SetBranchStatus("Electron_Energy", 1);
        chain->SetBranchStatus("Electron_charge", 1);
        chain->SetBranchStatus("Electron_etaSC", 1);
        chain->SetBranchStatus("Electron_phiSC", 1);
        chain->SetBranchStatus("Muon_pT", 1);
        chain->SetBranchStatus("Muon_eta", 1);
        chain->SetBranchStatus("Muon_phi", 1);
        chain->SetBranchStatus("Muon_charge", 1);
        chain->SetBranchStatus("Muon_Energy", 1);
        chain->SetBranchStatus("Muon_TuneP_pT", 1);
        chain->SetBranchStatus("Muon_TuneP_eta", 1);
        chain->SetBranchStatus("Muon_TuneP_phi", 1);
        chain->SetBranchStatus("Muon_trackerLayers", 1);
    }

    void TurnOffAllBranches()
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        chain->SetBranchStatus("*", 0);
    }

    void TurnOnBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0)
        {
            std::cout << "Error: no branches to turn off" << endl;
            return;
        }
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof())
        {
            ss >> br;
            if (chain->GetBranchStatus(br)==kFALSE)
            {
                chain->SetBranchStatus(br,1);
                std::cout << "Activating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already active\n";
        }
    }

    void TurnOffBranches(TString brNames)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if (brNames.Length()==0)
        {
            std::cout << "Error: no branches to turn off" << endl;
            return;
        }
        std::stringstream ss(brNames.Data());
        TString br;
        while (!ss.eof())
        {
            ss >> br;
            if (chain->GetBranchStatus(br)==kTRUE)
            {
                chain->SetBranchStatus(br,0);
                std::cout << "Deactivating branch <" << br << ">\n";
            }
            else std::cout << "Branch <" << br << "> is already inactive\n";
        }
    }

    void GetEvent(Int_t i)
    {
        if (!chain || File_Given == kFALSE)
        {
            std::cout << "Error: no chain provided" << endl;
            return;
        }
        if(!chain) return;

        chain->GetEntry(i);
    }

};
