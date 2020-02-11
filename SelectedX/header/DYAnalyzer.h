// -- Class for common functions used in DY differential cross section measurement analysis @ 13 TeV -- //
// -- Author: KyoengPil Lee, 05 Dec. 2015 -- //
// -- Author: Dalmin Pai, 11 Sep. 2017 -- //
// -- Synchronize acceptance cuts in dilepton channels : 19 Jan. 2018 -- //
// -- Change sub-leading pT cut (17->28) and isolation (trkiso->RelPFIso_dBeta) in muon channel : 24 Jan. 2018 -- //
// -- Add N-1 selection of trkiso : 26 Jan. 2018 -- //
#pragma once

#include "Object.h"
#include "NtupleHandle.h"
#include "SelectedX.h"
#include <TGraphErrors.h>
#include <TH2.h>
#include <TGraphAsymmErrors.h>
#include <iostream>

#define Lumi 35867 // -- from Run2016B to Run2016H, JSON. unit: /pb, Updated at 2017.07.30 -- //
#define Lumi_BtoF 19721 // -- from Run2016B to Run2016F, JSON. unit: /pb, Updated at 2018.05.17 -- //
#define Lumi_GtoH 16146 // -- from Run2016G to Run2016H, JSON. unit: /pb, Updated at 2018.05.17 -- //
#define nMassBin 43
#define nMassBin2 86
#define nPtBinEndcap 9//8
#define nPtBinBarrel 18//16

class DYAnalyzer
{
public:

	TString HLT;
	Double_t LeadPtCut;
	Double_t SubPtCut;
	Double_t LeadEtaCut;
	Double_t SubEtaCut;

//	Double_t PileUpWeight[52];
	Double_t PileUpWeight[75];

	// -- For efficiency SF of BtoF -- //
        Double_t Eff_Reco_data_BtoF[15][1];
        Double_t Eff_Reco_MC_BtoF[15][1];

        Double_t Eff_ID_data_BtoF[4][6];
        Double_t Eff_ID_MC_BtoF[4][6];

        Double_t Eff_Iso_data_BtoF[4][6];
        Double_t Eff_Iso_MC_BtoF[4][6];

	Double_t Eff_HLT_data_BtoF[4][7];
	Double_t Eff_HLT_MC_BtoF[4][7];

        Double_t Eff_ID_data_BtoF_lead[4][6];
        Double_t Eff_ID_MC_BtoF_lead[4][6];
        Double_t Eff_ID_data_BtoF_sublead[4][6];
        Double_t Eff_ID_MC_BtoF_sublead[4][6];

        Double_t Eff_Iso_data_BtoF_lead[4][6];
        Double_t Eff_Iso_MC_BtoF_lead[4][6];
        Double_t Eff_Iso_data_BtoF_sublead[4][6];
        Double_t Eff_Iso_MC_BtoF_sublead[4][6];

        Double_t Eff_HLT_data_BtoF_lead[4][8];
        Double_t Eff_HLT_MC_BtoF_lead[4][8];
        Double_t Eff_HLT_data_BtoF_sublead[4][8];
        Double_t Eff_HLT_MC_BtoF_sublead[4][8];

	// -- For efficiency SF of GtoH -- //
        Double_t Eff_Reco_data_GtoH[15][1];
        Double_t Eff_Reco_MC_GtoH[15][1];

        Double_t Eff_ID_data_GtoH[4][6];
        Double_t Eff_ID_MC_GtoH[4][6];

        Double_t Eff_Iso_data_GtoH[4][6];
        Double_t Eff_Iso_MC_GtoH[4][6];

	Double_t Eff_HLT_data_GtoH[4][7];
	Double_t Eff_HLT_MC_GtoH[4][7];

        Double_t Eff_ID_data_GtoH_lead[4][6];
        Double_t Eff_ID_MC_GtoH_lead[4][6];
        Double_t Eff_ID_data_GtoH_sublead[4][6];
        Double_t Eff_ID_MC_GtoH_sublead[4][6];

        Double_t Eff_Iso_data_GtoH_lead[4][6];
        Double_t Eff_Iso_MC_GtoH_lead[4][6];
        Double_t Eff_Iso_data_GtoH_sublead[4][6];
        Double_t Eff_Iso_MC_GtoH_sublead[4][6];

        Double_t Eff_HLT_data_GtoH_lead[4][8];
        Double_t Eff_HLT_MC_GtoH_lead[4][8];
        Double_t Eff_HLT_data_GtoH_sublead[4][8];
        Double_t Eff_HLT_MC_GtoH_sublead[4][8];

	// -- For efficiency SF of electron -- //
	Double_t Eff_Reco_data[30][1];
	Double_t Eff_Reco_MC[30][1];

        Double_t Eff_ID_data[20][5];
        Double_t Eff_ID_MC[20][5];

        Double_t Eff_HLT_Leg2_data[20][7];
        Double_t Eff_HLT_Leg2_MC[20][7];

        // -- PVz weights -- //
        Double_t PVzWeight[80];

        // -- Fake rates -- //
        const double ptbin_barrel[nPtBinBarrel+1] = {52,60,70,80,90,100,120,140,160,180,200,250,300,350,400,450,500,700,1000};
        const double ptbin_endcap[nPtBinEndcap+1] = {52,60,70,80,90,100,150,200,500,1000};
        Double_t FR_barrel[nPtBinBarrel];
        Double_t FR_endcap[nPtBinEndcap];

	// -- Constructor -- //
	DYAnalyzer(TString HLTname);

	// -- Setup accetpance cuts -- //
	void AssignAccThreshold(TString HLTname, TString *HLT, Double_t *LeadPtCut, Double_t *SubPtCut, Double_t *LeadEtaCut, Double_t *SubEtaCut);

        Bool_t SeparateDYLLSample_isHardProcess(TString Tag, NtupleHandle *ntuple);
	Bool_t Separate_ttbarSample(TString Tag, NtupleHandle *ntuple, vector<GenOthers> *GenTopCollection);

	// -- outdated -- //
	Bool_t SeparateDYLLSample(TString Tag, NtupleHandle *ntuple);

	//////////////////////////////////
	// -- Setup pileup weighting -- //
	//////////////////////////////////
        void SetupPileUpReWeighting(Bool_t isMC);
	Double_t PileUpWeightValue(Int_t PileUp_MC);

	// -- for 80X -- //
        void SetupPileUpReWeighting_80X(Bool_t isMC, TString ROOTFileName);
	Double_t PileUpWeightValue_80X(Int_t PileUp_MC);

        // -- For PVz reweighting -- //
        void SetupPVzWeights(Bool_t isMC, TString whichX, TString fileName);
        Double_t PVzWeightValue(Double_t PVz);

        /////////////////////////////////////////
	// -- Setup Efficiency scale factor -- //
	/////////////////////////////////////////
        // MUONS
        void SetupEfficiencyScaleFactor_BtoF();
	void SetupEfficiencyScaleFactor_GtoH();
        void SetupEfficiencyScaleFactor_BtoF_new();
        void SetupEfficiencyScaleFactor_GtoH_new();
        Double_t EfficiencySF_EventWeight_HLT_BtoF(Muon mu1, Muon mu2);
        Double_t EfficiencySF_EventWeight_HLT_BtoF(SelectedMuMu_t *MuMu);
        Double_t EfficiencySF_EventWeight_HLT_BtoF(TLorentzVector mu1, TLorentzVector mu2);
        Double_t EfficiencySF_EventWeight_HLT_BtoF_new(SelectedMuMu_t *MuMu);
        Double_t EfficiencySF_EventWeight_HLT_GtoH(Muon mu1, Muon mu2);
        Double_t EfficiencySF_EventWeight_HLT_GtoH(SelectedMuMu_t *MuMu);
        Double_t EfficiencySF_EventWeight_HLT_GtoH(TLorentzVector mu1, TLorentzVector mu2);
        Double_t EfficiencySF_EventWeight_HLT_GtoH_new(SelectedMuMu_t *MuMu);
        Int_t Find_muon_PtBin_Reco(Double_t Pt);
        Int_t Find_muon_PtBin_ID(Double_t Pt);
        Int_t Find_muon_PtBin_Iso(Double_t Pt);
        Int_t Find_muon_PtBin_Trig(Double_t Pt);
        Int_t Find_muon_PtBin_Trig_new(Double_t Pt);
        Int_t Find_muon_EtaBin_Reco(Double_t eta);
        Int_t Find_muon_EtaBin_ID(Double_t eta);
        Int_t Find_muon_EtaBin_Iso(Double_t eta);
        Int_t Find_muon_EtaBin_Trig(Double_t eta);

        // ELECTRONS
        void SetupEfficiencyScaleFactor_electron();
        Double_t EfficiencySF_EventWeight_electron(Electron ele1, Electron ele2);
        Double_t EfficiencySF_EventWeight_electron(SelectedEE_t *EE);
        Int_t Find_electron_PtBin_Reco(Double_t Pt);
        Int_t Find_electron_PtBin_ID(Double_t Pt);
        Int_t Find_electron_PtBin_Trig(Double_t Pt);
        Int_t Find_electron_EtaBin_Reco(Double_t eta);
        Int_t Find_electron_EtaBin_ID(Double_t eta);
        Int_t Find_electron_EtaBin_Trig(Double_t eta);
	
	////////////////////////////
	// -- Event Selections -- //
	////////////////////////////
	Bool_t EventSelection(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
        Bool_t EventSelection(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection, vector< Int_t >* Index); // -- output: 2 muons passing event selection conditions and their indices -- //  //*Derived by Marijus Ambrozas 2018.07*//
        Bool_t EventSelection_Mu50(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
	Bool_t EventSelection_minusDimuonVtxCut(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
        Bool_t EventSelection_Zdiff_13TeV(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //

        Bool_t EventSelection_Zdiff_13TeV_HighPt(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection, //*Derived by Marijus Ambrozas 2018.07*//
                                                 vector< Int_t >* Index, Int_t &IndexDi); // -- output: 2 muons passing event selection conditions and their indices -- //
        Bool_t EventSelection_Zdiff_13TeV_HighPt(vector< Muon > MuonCollection, LongSelectedMuMu_t *ntuple, vector< Muon >* SelectedMuonCollection,  //*Derived by Marijus Ambrozas 2018.07*//
                                                 vector< Int_t >* Index, Int_t &IndexDi); // -- output: 2 muons passing event selection conditions and their indices -- //
        Bool_t EventSelection_Zdiff_13TeV_HighPt(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection,  //*Derived by Marijus Ambrozas 2018.07*//
                                                 vector< Int_t >* Index); // -- output: 2 muons passing event selection conditions and their indices -- //
        Bool_t EventSelection_Zdiff_13TeV_HighPt(vector< Muon > MuonCollection, LongSelectedMuMu_t *ntuple, vector< Muon >* SelectedMuonCollection,   //*Derived by Marijus Ambrozas 2018.07*//
                                                 vector< Int_t >* Index); // -- output: 2 muons passing event selection conditions and their indices -- //

	Bool_t EventSelection_Dijet(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
	Bool_t EventSelection_Wjet(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
	Bool_t EventSelection_CheckMoreThanOneDimuonCand(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection, Bool_t& isMoreThanOneCand); // -- output: 2 muons passing event selection conditions -- //

	// -- for N-1 cuts of muon channel -- //
	Bool_t EventSelection_Zdiff_13TeV_HighPt1(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt2(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt3(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt4(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt5(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt6(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt7(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt8(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt9(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt10(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt11(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);

        // -- FAKE RATE -- //
        Bool_t EventSelection_FR(vector<Muon> MuonCollection, NtupleHandle *ntuple, vector<Muon> *SelectedMuonCollection_nume, vector<Muon> *SelectedMuonCollection_deno); // -- output: muons passing numerator and denominator selection -- //
        Bool_t EventSelection_FRdijetEst(vector<Muon> MuonCollection, NtupleHandle *ntuple, vector<Muon> *SelectedMuonCollection_fail); // -- output: two muons passing regular selection but failing isolation requirements-- //
        Bool_t EventSelection_FRsingleJetEst(vector<Muon> MuonCollection, NtupleHandle *ntuple,  vector<Muon> *SelectedMuonCollection); // -- output: one muon passing full regular selection and another one passing the same selection but failing isolation requirements-- //
        Bool_t EventSelection_FakeMuons_Triggerless(vector<Muon> MuonCollection, NtupleHandle *ntuple, vector<Muon> *SelectedMuonCollection); // -- output: two muons passing regular selection but failing isolation requirements (at least one) -- //
        Bool_t EventSelection_FR(vector<Electron> ElectronCollection, NtupleHandle *ntuple, vector<Electron> *SelectedElectronCollection); // Electron selection
        void SetupFRvalues(TString filename, TString type="sigCtrl_template");
        Double_t FakeRate(Double_t p_T, Double_t eta);


	Bool_t isPassAccCondition_Muon(Muon Mu1, Muon Mu2);
	Bool_t isPassAccCondition_GenLepton(GenLepton genlep1, GenLepton genlep2);
	void CompareMuon(Muon *Mu1, Muon *Mu2, Muon *leadMu, Muon *subMu);
	void CompareGenLepton(GenLepton *genlep1, GenLepton *genlep2, GenLepton *leadgenlep, GenLepton *subgenlep);
	void DimuonVertexProbNormChi2(NtupleHandle *ntuple, Double_t Pt1, Double_t Pt2, Double_t *VtxProb, Double_t *VtxNormChi2);
        void DimuonVertexProbNormChi2(LongSelectedMuMu_t *ntuple, Double_t Pt1, Double_t Pt2, Double_t *VtxProb, Double_t *VtxNormChi2);

	// -- for electron channel -- //
	Bool_t EventSelection_Electron(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection); // -- output: 2 electrons passing event selection conditions -- //
	Bool_t EventSelection_ElectronChannel_NminusPFIso(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection); // -- output: 2 electrons passing event selection conditions -- //
	Bool_t EventSelection_ElectronChannel(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection); // -- output: 2 electrons passing event selection conditions -- //

        Bool_t EventSelection_ElectronChannel(vector< Electron > ElectronCollection, NtupleHandle *ntuple,          // -- output: 2 electrons passing event selection conditions
                                              vector< Electron >* SelectedElectronCollection, vector< Int_t >* Sel_Index);  // and their indices inside the ntuple vectors -- //*Derived by Marijus Ambrozas 2018.08.02*//
        Bool_t EventSelection_ElectronChannel(vector< Electron > ElectronCollection, LongSelectedEE_t *ntuple,          // -- output: 2 electrons passing event selection conditions
                                              vector< Electron >* SelectedElectronCollection, vector< Int_t >* Sel_Index);  // and their indices inside the ntuple vectors -- //*Derived by Marijus Ambrozas 2018.08.06*//

	// -- for N-1 cuts of electron channel -- //
	Bool_t EventSelection_ElectronChannel1(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel2(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel3(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel4(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel5(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel6(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel7(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel8(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t isPassAccCondition_Electron(Electron Elec1, Electron Elec2);
	Bool_t isPassAccCondition_GenLepton_ECALGAP(GenLepton genlep1, GenLepton genlep2);
	void CompareElectron(Electron *Elec1, Electron *Elec2, Electron *leadElec, Electron *subElec);

	// -- pre-FSR functions -- //
	void PostToPreFSR_byDressedLepton(NtupleHandle *ntuple, GenLepton *genlep_postFSR, Double_t dRCut, GenLepton *genlep_preFSR, vector< GenOthers >* GenPhotonCollection);
	void PostToPreFSR_byDressedLepton_AllPhotons(NtupleHandle *ntuple, GenLepton *genlep_postFSR, Double_t dRCut, GenLepton *genlep_preFSR, vector< GenOthers >* GenPhotonCollection);
	TString DecideFSRType(GenLepton preFSR1, GenLepton preFSR2, GenLepton postFSR1, GenLepton postFSR2);
        Double_t Calc_dR_GenLeptons(GenLepton genlep1, GenLepton genlep2);
        Double_t Calc_dR_GenLepton_GenOthers(GenLepton genlep1, GenOthers genlep2);

	// -- miscellaneous -- //
	void GenMatching(TString MuonType, NtupleHandle* ntuple, vector<Muon>* MuonCollection);
        void ConvertToTunePInfo(Muon &mu);
        void PrintOutDoubleMuInfo(Muon mu1, Muon mu2);
        Double_t GenMuonPt(TString MuonType, NtupleHandle* ntuple, Muon reco_mu);

	// -- emu method -- //
        Bool_t EventSelection_emu_method(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple,
                                              vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection); // -- output: 1 muon and 1 electron passing event selection conditions -- //
        Bool_t EventSelection_emu_method(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple, // Derived by Marijus Ambrozas 2018.08.07
                                              vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection,
                                              Int_t &Sel_Index_Mu, Int_t &Sel_Index_Ele); // -- output: 1 muon and 1 electron passing event selection conditions and their indices -- //
        Bool_t EventSelection_emu_method(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, LongSelectedEMu_t *ntuple, // Derived by Marijus Ambrozas 2018.08.07
                                              vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection,
                                              Int_t &Sel_Index_Mu, Int_t &Sel_Index_Ele); // -- output: 1 muon and 1 electron passing event selection conditions and their indices -- //
        void emuVertexProbNormChi2(NtupleHandle *ntuple, Double_t Pt1, Double_t Pt2, Double_t *VtxProb, Double_t *VtxNormChi2);
        void emuVertexProbNormChi2(LongSelectedEMu_t *ntuple, Double_t Pt1, Double_t Pt2, Double_t *VtxProb, Double_t *VtxNormChi2);
	Double_t EfficiencySF_EventWeight_emu_BtoF(Muon mu1, Electron ele2);
        Double_t EfficiencySF_EventWeight_emu_BtoF(SelectedEMu_t *EMu);
	Double_t EfficiencySF_EventWeight_emu_GtoH(Muon mu1, Electron ele2);
        Double_t EfficiencySF_EventWeight_emu_GtoH(SelectedEMu_t *EMu);
};

DYAnalyzer::DYAnalyzer(TString HLTname)
{
        if(HLTname == "None")
	{
                std::cout << "===================================================" << endl;
                std::cout << "[No specific trigger setting ... basic constructor]" << endl;
                std::cout << "===================================================" << endl;
		
		HLT = "None";
		LeadPtCut = 9999;
		SubPtCut = 9999;
		LeadEtaCut = 9999;
		SubEtaCut = 9999;
	}
	else
	{
		this->AssignAccThreshold(HLTname, &HLT, &LeadPtCut, &SubPtCut, &LeadEtaCut, &SubEtaCut);
                std::cout << "===========================================================" << endl;
                std::cout << "Trigger: " << HLT << endl;
                std::cout << "leading lepton pT Cut: " << LeadPtCut << endl;
                std::cout << "Sub-leading lepton pT Cut: " << SubPtCut << endl;
                std::cout << "leading lepton Eta Cut: " << LeadEtaCut << endl;
                std::cout << "sub-leading lepton Eta Cut: " << SubEtaCut << endl;
                std::cout << "===========================================================" << endl;
	}

}

void DYAnalyzer::AssignAccThreshold(TString HLTname, TString *HLT, Double_t *LeadPtCut, Double_t *SubPtCut, Double_t *LeadEtaCut, Double_t *SubEtaCut)
{
        if(HLTname == "IsoMu20")
	{
		*HLT = "HLT_IsoMu20_v*";
		*LeadPtCut = 22;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
        else if(HLTname == "IsoMu20_OR_IsoTkMu20")
	{
		*HLT = "HLT_IsoMu20_v* || HLT_IsoTkMu20_v*";
		*LeadPtCut = 22;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
        else if(HLTname == "IsoMu24")
	{
		*HLT = "HLT_IsoMu24_v*";
		*LeadPtCut = 26;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
        else if(HLTname == "IsoMu24_OR_IsoTkMu24") // added at 2017.08.01
	{
		*HLT = "HLT_IsoMu24_v* || HLT_IsoTkMu24_v*";
		//*LeadPtCut = 26;
		//*SubPtCut = 10;
		*LeadPtCut = 28;
		*SubPtCut = 17;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
        else if(HLTname == "Mu45_eta2p1")
	{
		*HLT = "HLT_Mu45_eta2p1_v*";
		*LeadPtCut = 46;
		*SubPtCut = 10;
		*LeadEtaCut = 2.1;
		*SubEtaCut = 2.4;
	}
        else if(HLTname == "Mu50")
	{
		*HLT = "HLT_Mu50_v*";
                *LeadPtCut = 52;
                *SubPtCut = 0;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
        else if(HLTname == "IsoMu20_SymmetricPt25")
	{
		*HLT = "HLT_IsoMu20_v*";
		*LeadPtCut = 25;
		*SubPtCut = 25;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
        else if(HLTname == "Ele17Ele12")
	{
		*HLT = "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
		*LeadPtCut = 25;
		*SubPtCut = 15;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
        else if(HLTname == "Ele22_eta2p1")
	{
		*HLT = "HLT_Ele22_eta2p1_WPLoose_Gsf_v*"; // -- Exist only for the data; "HLT_Ele22_eta2p1_WP75_Gsf_v*" should be used for MC
		*LeadPtCut = 25;
		*SubPtCut = 15;
		*LeadEtaCut = 2.1;
		*SubEtaCut = 2.1;
	}
        else if(HLTname == "Ele22_eta2p1_NoEtaCut")
	{
		*HLT = "HLT_Ele22_eta2p1_WPLoose_Gsf_v*"; // -- Exist only for the data; "HLT_Ele22_eta2p1_WP75_Gsf_v*" should be used for MC
		*LeadPtCut = 25;
		*SubPtCut = 15;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
        else if(HLTname == "Pt_30_10_eta_2p5")
	{
		*HLT = "None"; // -- just for acceptance test -- //
		*LeadPtCut = 30;
		*SubPtCut = 10;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
        else if(HLTname == "Ele23_WPLoose")
	{
		*HLT = "HLT_Ele23_WPLoose_Gsf_v*"; // -- Exist only for the data; "HLT_Ele22_eta2p1_WP75_Gsf_v*" should be used for MC
		*LeadPtCut = 30;
		*SubPtCut = 10;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
        else if(HLTname == "Ele23Ele12") // updated at 2017.02.21 by Dalmin Pai
	{
		*HLT = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
		*LeadPtCut = 28;
		*SubPtCut = 17;
		//*LeadEtaCut = 2.5; // -- Later it should exclude ECAL gap
		//*SubEtaCut = 2.5; // -- Later it should exclude ECAL gap
		*LeadEtaCut = 2.4; // -- Later it should exclude ECAL gap
		*SubEtaCut = 2.4; // -- Later it should exclude ECAL gap
	}
        else if (HLTname == "Photon_OR")
        {
            *HLT = "HLT_Photon*";
            *LeadPtCut = 28;
            *SubPtCut = 17;
            *LeadEtaCut = 2.4; // -- Later it should exclude ECAL gap
            *SubEtaCut = 2.4; // -- Later it should exclude ECAL gap
        }
	else
	{ 
                std::cout << "Wrong HLT name!: " << HLTname << endl;
		return; 
	}

}

Bool_t DYAnalyzer::Separate_ttbarSample(TString Tag, NtupleHandle *ntuple, vector<GenOthers> *GenTopCollection)
{
	Bool_t GenFlag = kFALSE;

	// -- Seperate ttbar events -- //
        if(Tag.Contains("ttbar"))
	{
		vector<GenOthers> GenOthersCollection;
		Int_t NGenOthers = ntuple->nGenOthers;
		for(Int_t i_gen=0; i_gen<NGenOthers; i_gen++)
		{
			GenOthers genothers;
			genothers.FillFromNtuple(ntuple, i_gen);
                        if(abs(genothers.ID) == 6 && genothers.isHardProcess)
                                GenOthersCollection.push_back(genothers);
		}

                if(GenOthersCollection.size() == 2) // -- Select the ttbar events from hard-process -- //
		{
			// -- Check top & anti-top pair -- //
                        if(GenOthersCollection[0].ID == GenOthersCollection[1].ID)
				printf("%d %d\n", GenOthersCollection[0].ID, GenOthersCollection[1].ID);

                        //if(Tag == "ttbar") // -- Select only evetns withtin M < 700 -- //
                        if(Tag == "ttbar" || Tag == "ttbarBackup") // -- Select only evetns withtin M < 700 -- //
			{
				TLorentzVector v1 = GenOthersCollection[0].Momentum;
				TLorentzVector v2 = GenOthersCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 700)
                                //if(reco_M > -999)
				{
					GenFlag = kTRUE;
                                        GenTopCollection->push_back(GenOthersCollection[0]);
                                        GenTopCollection->push_back(GenOthersCollection[1]);
				}
			}
			else // ex: ttbar_M700to1000, ttbar_M1000toInf
			{
				GenFlag = kTRUE;
                                GenTopCollection->push_back(GenOthersCollection[0]);
                                GenTopCollection->push_back(GenOthersCollection[1]);
			}
		}
		else
		{
			printf("Wrong? : more than two!!\n");
			printf("%d %d %d\n", GenOthersCollection[0].ID, GenOthersCollection[1].ID, GenOthersCollection[2].ID); //Check upto 3rd one
		}
	}
	// -- other cases(e.g. DY, WJets, Diboson...): pass
	else
		GenFlag = kTRUE; 

	return GenFlag;
}

Bool_t DYAnalyzer::SeparateDYLLSample_isHardProcess(TString Tag, NtupleHandle *ntuple)
{
	Bool_t GenFlag = kFALSE;

	// -- Seperate DYMuMu events from DYTauTau  -- //
        if(Tag.Contains("DYMuMu"))
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
                        if(genlep.isMuon() && genlep.isHardProcess)
                                GenLeptonCollection.push_back(genlep);
		}

                if(GenLeptonCollection.size() == 2) // -- Select the events containing 2 muons from hard-process -- //
		{
                        if(Tag == "DYMuMu_M50to200") // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 200)
					GenFlag = kTRUE;
			}
                        else if(Tag == "DYMuMu_M50to400") // -- Select only evetns withtin 50 < M < 400 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 400)
					GenFlag = kTRUE;
			}
                        else if(Tag == "DYMuMu_M50to100" || Tag == "DYMuMu_Photos_M50to100") // -- Select only evetns withtin 50 < M < 100 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 100)
					GenFlag = kTRUE;
			}
                        else if(Tag == "DYMuMu_M50to120")
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 120)
					GenFlag = kTRUE;
			}
                        else if(Tag == "DYMuMu_M120to200")
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M > 120)
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
        else if(Tag.Contains("DYEE"))
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
                        if(genlep.isElectron() && genlep.isHardProcess)
                                GenLeptonCollection.push_back(genlep);
		}

                if(GenLeptonCollection.size() == 2) // -- Select the events containing 2 electrons from hard-process -- //
		{
                        if(Tag == "DYEE_M50to200") // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 200)
					GenFlag = kTRUE;
			}
                        else if(Tag == "DYEE_M50to100") // -- Select only evetns withtin 50 < M < 100 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 100)
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
	// -- Separate DYTauTau events from MuMu events -- //
        else if(Tag.Contains("DYTauTau"))
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
                        if(abs(genlep.ID) == 15 && genlep.isHardProcess)
                                GenLeptonCollection.push_back(genlep);
		}

                if(GenLeptonCollection.size() == 2) // -- Select the events containing 2 taus from hard-process -- //
		{
			GenFlag = kTRUE;
		}
	}
	// -- Madgraph sample -- //
        else if(Tag.Contains("Madgraph"))
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
                        if(genlep.isMuon() && genlep.isHardProcess)
                                GenLeptonCollection.push_back(genlep);
		}

                if(GenLeptonCollection.size() == 2) // -- Select the events containing 2 muons from hard-process -- //
		{
                        if(Tag == "Madgraph_M50to150") // -- Select only evetns withtin 50 < M < 150 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 150)
					GenFlag = kTRUE;
			}
                        else if(Tag == "Madgraph_M10to50")
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M > 10)
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
        else if(Tag.Contains("DY"))
        {
                vector<GenLepton> GenLeptonCollection;
                Int_t NGenLeptons = ntuple->gnpair;
                for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
                {
                        GenLepton genlep;
                        genlep.FillFromNtuple(ntuple, i_gen);
                        if(genlep.isHardProcess)
                                GenLeptonCollection.push_back(genlep);
                }

                if(GenLeptonCollection.size() == 2) // -- Select the events containing 2 electrons from hard-process -- //
                {
                        if(Tag == "DY_M50to100") // -- Select only events withtin 50 < M < 100 -- //
                        {
                                TLorentzVector v1 = GenLeptonCollection[0].Momentum;
                                TLorentzVector v2 = GenLeptonCollection[1].Momentum;
                                Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 100)
                                        GenFlag = kTRUE;
                        }
                        else
                                GenFlag = kTRUE;
                }
        }
	// -- other cases(e.g. ttbar, WJets, Diboson...): pass
	else
		GenFlag = kTRUE; 

	return GenFlag;
}

Bool_t DYAnalyzer::SeparateDYLLSample(TString Tag, NtupleHandle *ntuple)
{
	Bool_t GenFlag = kFALSE;

	// -- Seperate DYMuMu events from DYTauTau  -- //
        if(Tag.Contains("DYMuMu"))
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
                        if(genlep.isMuon() && genlep.fromHardProcessFinalState)
                                GenLeptonCollection.push_back(genlep);
		}

                if(GenLeptonCollection.size() == 2) // -- Select the events containing 2 muons from hard-process -- //
		{
                        if(Tag == "DYMuMu_M50to200") // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 200)
					GenFlag = kTRUE;
			}
                        else if(Tag == "DYMuMu_M50to400") // -- Select only evetns withtin 50 < M < 400 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 400)
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
        else if(Tag.Contains("DYEE"))
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
                        if(genlep.isElectron() && genlep.fromHardProcessFinalState)
                                GenLeptonCollection.push_back(genlep);
		}

                if(GenLeptonCollection.size() == 2) // -- Select the events containing 2 electrons from hard-process -- //
		{
                        if(Tag == "DYEE_M50to200") // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
                                if(reco_M < 200)
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
	// -- Separate DYTauTau events from MuMu events -- //
        else if(Tag.Contains("DYTauTau"))
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
                        if(abs(genlep.ID) == 15 && genlep.fromHardProcessDecayed)
                                GenLeptonCollection.push_back(genlep);
		}

                if(GenLeptonCollection.size() == 2) // -- Select the events containing 2 taus from hard-process -- //
		{
			GenFlag = kTRUE;
		}
	}
	// -- other cases(e.g. ttbar, WJets, Diboson...): pass
	else
		GenFlag = kTRUE; 

	return GenFlag;
}

void DYAnalyzer::SetupPileUpReWeighting(Bool_t isMC)
{
        if(isMC == kFALSE) // -- for data -- //
	{
		for(Int_t i=0; i<52; i++)
			PileUpWeight[i] = 1;

		return;
	}
	
	// -- Only for the MC -- //
	TFile *f = new TFile("/home/kplee/CommonCodes/DrellYanAnalysis/ROOTFile_PUReWeight_v20160208_2nd_71mb.root");
	f->cd();
	TH1D *h_weight = (TH1D*)f->Get("h_PUReWeights");
        if(h_weight == NULL)
	{
                std::cout << "ERROR! ... No Weight histogram!"<< endl;
		return;
	}

	for(Int_t i=0; i<52; i++)
	{
		Int_t i_bin = i+1;
		PileUpWeight[i] = h_weight->GetBinContent(i_bin);
	}
}

Double_t DYAnalyzer::PileUpWeightValue(Int_t PileUp_MC)
{
        if(PileUp_MC < 0 || PileUp_MC > 51)
	{
                std::cout << "[PileUp_MC = " << PileUp_MC << "]: NO CORRESPONDING PU Weight! ... it returns 0" << endl;
		return 0;
	}
	return PileUpWeight[PileUp_MC];
}

void DYAnalyzer::SetupPileUpReWeighting_80X(Bool_t isMC, TString ROOTFileName)
{
        if(isMC == kFALSE) // -- for data -- //
	{
		for(Int_t i=0; i<75; i++)
			PileUpWeight[i] = 1;

		return;
	}
	
	// -- Only for the MC -- //
        TString FileLocation = "./etc/PileUp/80X/"+ROOTFileName;
	TFile *f = new TFile(FileLocation);
	TH1D *h_weight = (TH1D*)f->Get("h_PUReWeights");
        if(h_weight == NULL)
	{
                std::cout << "ERROR! ... No Weight histogram!"<< endl;
		return;
	}

	for(Int_t i=0; i<75; i++)
	{
		Int_t i_bin = i+1;
		PileUpWeight[i] = h_weight->GetBinContent(i_bin);
	}
}

Double_t DYAnalyzer::PileUpWeightValue_80X(Int_t PileUp_MC)
{
        if(PileUp_MC < 0 || PileUp_MC > 74)
	{
                std::cout << "[PileUp_MC = " << PileUp_MC << "]: NO CORRESPONDING PU Weight! ... it returns 0" << endl;
		return 0;
	}
	return PileUpWeight[PileUp_MC];
}


void DYAnalyzer::SetupEfficiencyScaleFactor_BtoF()
{
	TString Location = "./etc/effSF/effSF_muon/";
        std::cout << "[Tag&Probe efficiency is from " << Location+"*BtoF.root" << "]" << endl;

        // RECO
        TFile *f_reco = new TFile(Location+"Tracking_SF_RunBtoF.root");
        TGraphAsymmErrors *h_Reco_ratio = (TGraphAsymmErrors*)f_reco->Get("ratio_eff_eta3_dr030e030_corr");
        Int_t nEtaBins_reco = h_Reco_ratio->GetN();
        Int_t nPtBins_reco = 1;

        Double_t eta_reco[15]; Double_t eta_reco_ordered[15];
        Double_t SF_reco[15]; Double_t SF_reco_ordered[15];

        for(Int_t iter_y = 0; iter_y < nPtBins_reco; iter_y++) // in this case (2d distribution): x=eta, y=pT
        {
                for(Int_t i=0; i<nEtaBins_reco; i++)
                {
                        h_Reco_ratio->GetPoint(i, eta_reco[i], SF_reco[i]); // in this case: x=eta, y=SF
                }

                // -- Rearrangement in order of x (eta) values -- //
                Double_t etamin; // the minimum value in given iteration
                Double_t etalow = -9999; // the minimum value of the previous iteration
                for(Int_t j=0; j<nEtaBins_reco; j++)
                {
                        Int_t jj = -9999;

                        etamin = 9999;
                        for(Int_t k=0; k<nEtaBins_reco; k++)
                        {
                                if(etalow < eta_reco[k] && eta_reco[k] < etamin) // the lowest number but higher than the one from previos iteration
                                {
                                        jj = k;
                                        etamin = eta_reco[k];
                                }
                        }
                        eta_reco_ordered[j] = eta_reco[jj];
                        SF_reco_ordered[j] = SF_reco[jj];

                        etalow = etamin;
                } // End of rearrangement

                for(Int_t iter_x = 0; iter_x < nEtaBins_reco; iter_x++)
                {
                        Eff_Reco_data_BtoF[iter_x][iter_y] = SF_reco_ordered[iter_x]; // actually, it is the scale factor.
                        Eff_Reco_MC_BtoF[iter_x][iter_y] = 1;
                }
        }

        // ID
        TFile *f_ID = new TFile(Location+"ID_SF_RunBtoF.root");
//	TH2F *h_ID_data = (TH2F*)f_ID->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA");
//	TH2F *h_ID_MC = (TH2F*)f_ID->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesMC/abseta_pair_ne_MC");
        TH2F *h_ID_data = (TH2F*)f_ID->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/abseta_pt_DATA");
        TH2F *h_ID_MC = (TH2F*)f_ID->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");

        Int_t nEtaBins_ID = h_ID_data->GetNbinsX();
        Int_t nPtBins_ID = h_ID_data->GetNbinsY();

        for(Int_t iter_x = 0; iter_x < nEtaBins_ID; iter_x++)
        {
                for(Int_t iter_y = 0; iter_y < nPtBins_ID; iter_y++)
                {
                        Int_t i_etabin = iter_x + 1;
                        Int_t i_ptbin = iter_y + 1;

                        Double_t ID_data = h_ID_data->GetBinContent(i_etabin, i_ptbin);
                        Double_t ID_MC = h_ID_MC->GetBinContent(i_etabin, i_ptbin);

                        Eff_ID_data_BtoF[iter_x][iter_y] = ID_data;
                        Eff_ID_MC_BtoF[iter_x][iter_y] = ID_MC;
                }
        }

        // ISO
        TFile *f_iso = new TFile(Location+"ISO_SF_RunBtoF.root");
//	TH2F *h_Iso_data = (TH2F*)f_iso->Get("tkLooseISO_highptID_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA");
//	TH2F *h_Iso_MC = (TH2F*)f_iso->Get("tkLooseISO_highptID_newpt_eta/efficienciesMC/abseta_pair_ne_MC");
        TH2F *h_Iso_data = (TH2F*)f_iso->Get("TightISO_TightID_pt_eta/efficienciesDATA/abseta_pt_DATA");
        TH2F *h_Iso_MC = (TH2F*)f_iso->Get("TightISO_TightID_pt_eta/efficienciesMC/abseta_pt_MC");

        Int_t nEtaBins_iso = h_Iso_data->GetNbinsX();
        Int_t nPtBins_iso = h_Iso_data->GetNbinsY();

        for(Int_t iter_x = 0; iter_x < nEtaBins_iso; iter_x++)
        {
                for(Int_t iter_y = 0; iter_y < nPtBins_iso; iter_y++)
                {
                        Int_t i_etabin = iter_x + 1;
                        Int_t i_ptbin = iter_y + 1;

                        Double_t Iso_data = h_Iso_data->GetBinContent(i_etabin, i_ptbin);
                        Double_t Iso_MC = h_Iso_MC->GetBinContent(i_etabin, i_ptbin);

                        Eff_Iso_data_BtoF[iter_x][iter_y] = Iso_data;
                        Eff_Iso_MC_BtoF[iter_x][iter_y] = Iso_MC;
                }
        }

        // TRIGGER
        TFile *f_HLT = new TFile(Location+"Trigger_SF_RunBtoF.root");
        TH2F *h_HLT_data = (TH2F*)f_HLT->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
        TH2F *h_HLT_MC = (TH2F*)f_HLT->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/abseta_pt_MC");

        Int_t nEtaBins_HLT = h_HLT_data->GetNbinsX();
        Int_t nPtBins_HLT = h_HLT_data->GetNbinsY();

        for(Int_t iter_x = 0; iter_x < nEtaBins_HLT; iter_x++)
        {
                for(Int_t iter_y = 0; iter_y < nPtBins_HLT; iter_y++)
                {
                        Int_t i_etabin = iter_x + 1;
                        Int_t i_ptbin = iter_y + 1;

                        Double_t HLT_data = h_HLT_data->GetBinContent(i_etabin, i_ptbin);
                        Double_t HLT_MC = h_HLT_MC->GetBinContent(i_etabin, i_ptbin);

                        Eff_HLT_data_BtoF[iter_x][iter_y] = HLT_data;
                        Eff_HLT_MC_BtoF[iter_x][iter_y] = HLT_MC;
                }
        }
        std::cout << "Setting for efficiency correction factors (BtoF) is completed" << endl;
} // End of SetupEfficiencyScaleFactor_BtoF


void DYAnalyzer::SetupEfficiencyScaleFactor_BtoF_new()
{
    TString Location = "./etc/effSF/effSF_muon/";
    TFile *f1 = new TFile(Location+"New_SF_RunBtoF.root");
    std::cout << "[Tag&Probe efficiency is from " << Location+"New_SF_RunBtoF.root" << "]" << endl;


//    // RECO
//    TFile *f_reco = new TFile(Location+"Tracking_SF_RunBtoF.root");
//    TGraphAsymmErrors *h_Reco_ratio = (TGraphAsymmErrors*)f_reco->Get("ratio_eff_eta3_dr030e030_corr");
//    Int_t nEtaBins_reco = h_Reco_ratio->GetN();
//    Int_t nPtBins_reco = 1;

//    Double_t eta_reco[15]; Double_t eta_reco_ordered[15];
//    Double_t SF_reco[15]; Double_t SF_reco_ordered[15];

//    for(Int_t iter_y = 0; iter_y < nPtBins_reco; iter_y++) // in this case (2d distribution): x=eta, y=pT
//    {
//            for(Int_t i=0; i<nEtaBins_reco; i++)
//            {
//                    h_Reco_ratio->GetPoint(i, eta_reco[i], SF_reco[i]); // in this case: x=eta, y=SF
//            }

//            // -- Rearrangement in order of x (eta) values -- //
//            Double_t etamin; // the minimum value in given iteration
//            Double_t etalow = -9999; // the minimum value of the previous iteration
//            for(Int_t j=0; j<nEtaBins_reco; j++)
//            {
//                    Int_t jj = -9999;

//                    etamin = 9999;
//                    for(Int_t k=0; k<nEtaBins_reco; k++)
//                    {
//                            if(etalow < eta_reco[k] && eta_reco[k] < etamin) // the lowest number but higher than the one from previos iteration
//                            {
//                                    jj = k;
//                                    etamin = eta_reco[k];
//                            }
//                    }
//                    eta_reco_ordered[j] = eta_reco[jj];
//                    SF_reco_ordered[j] = SF_reco[jj];

//                    etalow = etamin;
//            } // End of rearrangement

//            for(Int_t iter_x = 0; iter_x < nEtaBins_reco; iter_x++)
//            {
//                    Eff_Reco_data_BtoF[iter_x][iter_y] = SF_reco_ordered[iter_x]; // actually, it is the scale factor.
//                    Eff_Reco_MC_BtoF[iter_x][iter_y] = 1;
//            }
//    }

    // ID
    TGraphAsymmErrors *h_data_lead_ID_eff[4];
    TGraphAsymmErrors *h_data_sublead_ID_eff[4];
    TGraphAsymmErrors *h_mc_lead_ID_eff[4];
    TGraphAsymmErrors *h_mc_sublead_ID_eff[4];

    h_data_lead_ID_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Leading_pteta_abseta0");
    h_data_lead_ID_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Leading_pteta_abseta1");
    h_data_lead_ID_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Leading_pteta_abseta2");
    h_data_lead_ID_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Leading_pteta_abseta3");

    h_data_sublead_ID_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Subleading_pteta_abseta0");
    h_data_sublead_ID_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Subleading_pteta_abseta1");
    h_data_sublead_ID_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Subleading_pteta_abseta2");
    h_data_sublead_ID_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Subleading_pteta_abseta3");

    h_mc_lead_ID_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Leading_pteta_abseta0");
    h_mc_lead_ID_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Leading_pteta_abseta1");
    h_mc_lead_ID_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Leading_pteta_abseta2");
    h_mc_lead_ID_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Leading_pteta_abseta3");

    h_mc_sublead_ID_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Subleading_pteta_abseta0");
    h_mc_sublead_ID_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Subleading_pteta_abseta1");
    h_mc_sublead_ID_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Subleading_pteta_abseta2");
    h_mc_sublead_ID_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Subleading_pteta_abseta3");

    Int_t nEtaBins_ID = 4;
    Int_t nPtBins_ID = 6;

    Double_t x_data_lead_ID[6];
    Double_t y_data_lead_ID[6];

    Double_t x_data_sublead_ID[6];
    Double_t y_data_sublead_ID[6];

    Double_t x_mc_lead_ID[6];
    Double_t y_mc_lead_ID[6];

    Double_t x_mc_sublead_ID[6];
    Double_t y_mc_sublead_ID[6];

    TGraphAsymmErrors *h_Data_lead_ID_eff;
    TGraphAsymmErrors *h_Data_sublead_ID_eff;
    TGraphAsymmErrors *h_MC_lead_ID_eff;
    TGraphAsymmErrors *h_MC_sublead_ID_eff;

    for(Int_t iter_x = 0; iter_x < nEtaBins_ID; iter_x++)
    {
        h_Data_lead_ID_eff = (TGraphAsymmErrors*)h_data_lead_ID_eff[iter_x]->Clone();
        h_Data_sublead_ID_eff = (TGraphAsymmErrors*)h_data_sublead_ID_eff[iter_x]->Clone();
        h_MC_lead_ID_eff = (TGraphAsymmErrors*)h_mc_lead_ID_eff[iter_x]->Clone();
        h_MC_sublead_ID_eff = (TGraphAsymmErrors*)h_mc_sublead_ID_eff[iter_x]->Clone();

        for(Int_t i = 0; i < nPtBins_ID; i++)
        {
            h_Data_sublead_ID_eff->GetPoint(i, x_data_sublead_ID[i], y_data_sublead_ID[i]);
            h_MC_sublead_ID_eff->GetPoint(i, x_mc_sublead_ID[i], y_mc_sublead_ID[i]);

            if(i == 0)
            {
                //It is just to initialize. They don't work in SF
                x_data_lead_ID[0] = 1;
                y_data_lead_ID[0] = 1;

                x_mc_lead_ID[0] = 1;
                y_mc_lead_ID[0] = 1;

                continue;
            }

            h_Data_lead_ID_eff->GetPoint(i-1, x_data_lead_ID[i], y_data_lead_ID[i]);
            h_MC_lead_ID_eff->GetPoint(i-1, x_mc_lead_ID[i], y_mc_lead_ID[i]);
        }

        for(Int_t iter_y = 0; iter_y < nPtBins_ID; iter_y++)
        {
            Eff_ID_data_BtoF_lead[iter_x][iter_y] = y_data_lead_ID[iter_y];
            Eff_ID_MC_BtoF_lead[iter_x][iter_y] = y_mc_lead_ID[iter_y];

            Eff_ID_data_BtoF_sublead[iter_x][iter_y] = y_data_sublead_ID[iter_y];
            Eff_ID_MC_BtoF_sublead[iter_x][iter_y] = y_mc_sublead_ID[iter_y];
        }
    }

    // ISO
    TGraphAsymmErrors *h_data_lead_iso_eff[4];
    TGraphAsymmErrors *h_data_sublead_iso_eff[4];
    TGraphAsymmErrors *h_mc_lead_iso_eff[4];
    TGraphAsymmErrors *h_mc_sublead_iso_eff[4];

    h_data_lead_iso_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Leading_pteta_abseta0");
    h_data_lead_iso_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Leading_pteta_abseta1");
    h_data_lead_iso_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Leading_pteta_abseta2");
    h_data_lead_iso_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Leading_pteta_abseta3");

    h_data_sublead_iso_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta0");
    h_data_sublead_iso_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta1");
    h_data_sublead_iso_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta2");
    h_data_sublead_iso_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta3");

    h_mc_lead_iso_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Leading_pteta_abseta0");
    h_mc_lead_iso_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Leading_pteta_abseta1");
    h_mc_lead_iso_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Leading_pteta_abseta2");
    h_mc_lead_iso_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Leading_pteta_abseta3");

    h_mc_sublead_iso_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta0");
    h_mc_sublead_iso_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta1");
    h_mc_sublead_iso_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta2");
    h_mc_sublead_iso_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta3");

    Int_t nEtaBins_iso = 4;
    Int_t nPtBins_iso = 6;

    Double_t x_data_lead_iso[6];
    Double_t y_data_lead_iso[6];

    Double_t x_data_sublead_iso[6];
    Double_t y_data_sublead_iso[6];

    Double_t x_mc_lead_iso[6];
    Double_t y_mc_lead_iso[6];

    Double_t x_mc_sublead_iso[6];
    Double_t y_mc_sublead_iso[6];

    TGraphAsymmErrors *h_Data_lead_iso_eff;
    TGraphAsymmErrors *h_Data_sublead_iso_eff;
    TGraphAsymmErrors *h_MC_lead_iso_eff;
    TGraphAsymmErrors *h_MC_sublead_iso_eff;

    for(Int_t iter_x = 0; iter_x < nEtaBins_iso; iter_x++)
    {
        h_Data_lead_iso_eff = (TGraphAsymmErrors*)h_data_lead_iso_eff[iter_x]->Clone();
        h_Data_sublead_iso_eff = (TGraphAsymmErrors*)h_data_sublead_iso_eff[iter_x]->Clone();
        h_MC_lead_iso_eff = (TGraphAsymmErrors*)h_mc_lead_iso_eff[iter_x]->Clone();
        h_MC_sublead_iso_eff = (TGraphAsymmErrors*)h_mc_sublead_iso_eff[iter_x]->Clone();

        for(Int_t i=0; i<nPtBins_iso; i++) //0 to 5
        {
            h_Data_sublead_iso_eff->GetPoint(i, x_data_sublead_iso[i], y_data_sublead_iso[i]);
            h_MC_sublead_iso_eff->GetPoint(i, x_mc_sublead_iso[i], y_mc_sublead_iso[i]);

            if(i == 0)
            {
                //It is just to initialize. They don't work in SF
                x_data_lead_iso[0] = 1;
                y_data_lead_iso[0] = 1;

                x_mc_lead_iso[0] = 1;
                y_mc_lead_iso[0] = 1;

                continue;
            }

            h_Data_lead_iso_eff->GetPoint(i-1, x_data_lead_iso[i], y_data_lead_iso[i]);
            h_MC_lead_iso_eff->GetPoint(i-1, x_mc_lead_iso[i], y_mc_lead_iso[i]);
        }

        for(Int_t iter_y = 0; iter_y < nPtBins_iso; iter_y++)
        {
            Eff_Iso_data_BtoF_lead[iter_x][iter_y] = y_data_lead_iso[iter_y];
            Eff_Iso_MC_BtoF_lead[iter_x][iter_y] = y_mc_lead_iso[iter_y];

            Eff_Iso_data_BtoF_sublead[iter_x][iter_y] = y_data_sublead_iso[iter_y];
            Eff_Iso_MC_BtoF_sublead[iter_x][iter_y] = y_mc_sublead_iso[iter_y];
        }
    }

    // TRIGGER
    TGraphAsymmErrors *h_data_lead_trig_eff[4];
    TGraphAsymmErrors *h_data_sublead_trig_eff[4];
    TGraphAsymmErrors *h_mc_lead_trig_eff[4];
    TGraphAsymmErrors *h_mc_sublead_trig_eff[4];

    h_data_lead_trig_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Leading_pteta_abseta0");
    h_data_lead_trig_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Leading_pteta_abseta1");
    h_data_lead_trig_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Leading_pteta_abseta2");
    h_data_lead_trig_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Leading_pteta_abseta3");

    h_data_sublead_trig_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Subleading_pteta_abseta0");
    h_data_sublead_trig_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Subleading_pteta_abseta1");
    h_data_sublead_trig_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Subleading_pteta_abseta2");
    h_data_sublead_trig_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Subleading_pteta_abseta3");

    h_mc_lead_trig_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Leading_pteta_abseta0");
    h_mc_lead_trig_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Leading_pteta_abseta1");
    h_mc_lead_trig_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Leading_pteta_abseta2");
    h_mc_lead_trig_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Leading_pteta_abseta3");

    h_mc_sublead_trig_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Subleading_pteta_abseta0");
    h_mc_sublead_trig_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Subleading_pteta_abseta1");
    h_mc_sublead_trig_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Subleading_pteta_abseta2");
    h_mc_sublead_trig_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_pair_dPhiPrimeDeg70_and_Subleading_pteta_abseta3");

    Int_t nEtaBins_HLT = 4;
    Int_t nPtBins_HLT = 8;

    Double_t x_data_lead_trig[8];
    Double_t y_data_lead_trig[8];

    Double_t x_data_sublead_trig[8];
    Double_t y_data_sublead_trig[8];

    Double_t x_mc_lead_trig[8];
    Double_t y_mc_lead_trig[8];

    Double_t x_mc_sublead_trig[8];
    Double_t y_mc_sublead_trig[8];

    TGraphAsymmErrors *h_Data_lead_trig_eff;
    TGraphAsymmErrors *h_Data_sublead_trig_eff;
    TGraphAsymmErrors *h_MC_lead_trig_eff;
    TGraphAsymmErrors *h_MC_sublead_trig_eff;

    for(Int_t iter_x = 0; iter_x < nEtaBins_HLT; iter_x++)
    {
        h_Data_lead_trig_eff = (TGraphAsymmErrors*)h_data_lead_trig_eff[iter_x]->Clone();
        h_Data_sublead_trig_eff = (TGraphAsymmErrors*)h_data_sublead_trig_eff[iter_x]->Clone();
        h_MC_lead_trig_eff = (TGraphAsymmErrors*)h_mc_lead_trig_eff[iter_x]->Clone();
        h_MC_sublead_trig_eff = (TGraphAsymmErrors*)h_mc_sublead_trig_eff[iter_x]->Clone();

        for(Int_t i=0; i<nPtBins_HLT; i++)
        {
            h_Data_lead_trig_eff->GetPoint(i, x_data_lead_trig[i], y_data_lead_trig[i]);
            h_MC_lead_trig_eff->GetPoint(i, x_mc_lead_trig[i], y_mc_lead_trig[i]);

            if(iter_x == 2 && i == 7)
            {
                h_Data_sublead_trig_eff->GetPoint(i, x_data_sublead_trig[i], y_data_sublead_trig[i]);

                x_mc_sublead_trig[i] = 1;
                y_mc_sublead_trig[i] = y_data_sublead_trig[i-1];

                continue;
            }
            else if(iter_x == 3 && i >= 6)
            {
                x_data_sublead_trig[i] = 1;
                y_data_sublead_trig[i] = y_data_sublead_trig[i-1];

                x_mc_sublead_trig[i] = 1;
                y_mc_sublead_trig[i] = y_data_sublead_trig[i-1];

                continue;
            }

            h_Data_sublead_trig_eff->GetPoint(i, x_data_sublead_trig[i], y_data_sublead_trig[i]);
            h_MC_sublead_trig_eff->GetPoint(i, x_mc_sublead_trig[i], y_mc_sublead_trig[i]);
        }

        for(Int_t iter_y = 0; iter_y < nPtBins_HLT; iter_y++)
        {
            Eff_HLT_data_BtoF_lead[iter_x][iter_y] = y_data_lead_trig[iter_y];
            Eff_HLT_MC_BtoF_lead[iter_x][iter_y] = y_mc_lead_trig[iter_y];

            Eff_HLT_data_BtoF_sublead[iter_x][iter_y] = y_data_sublead_trig[iter_y];
            Eff_HLT_MC_BtoF_sublead[iter_x][iter_y] = y_mc_sublead_trig[iter_y];
        }
    }
    std::cout << "Setting for efficiency correction factors (BtoF) is completed" << endl;
} // End of SetupEfficiencyScaleFactor_BtoF_new

void DYAnalyzer::SetupEfficiencyScaleFactor_GtoH()
{
	TString Location = "./etc/effSF/effSF_muon/";
        std::cout << "[Tag&Probe efficiency is from " << Location+"*GtoH.root" << "]" << endl;

        // RECO
        TFile *f_reco = new TFile(Location+"Tracking_SF_RunGtoH.root");
        TGraphAsymmErrors *h_Reco_ratio = (TGraphAsymmErrors*)f_reco->Get("ratio_eff_eta3_dr030e030_corr");
        Int_t nEtaBins_reco = h_Reco_ratio->GetN();
        Int_t nPtBins_reco = 1;

        Double_t eta_reco[15]; Double_t eta_reco_ordered[15];
        Double_t SF_reco[15]; Double_t SF_reco_ordered[15];

        for(Int_t iter_y = 0; iter_y < nPtBins_reco; iter_y++) // in this case (2d distribution): x=eta, y=pT
        {
                for(Int_t i=0; i<nEtaBins_reco; i++)
                {
                        h_Reco_ratio->GetPoint(i, eta_reco[i], SF_reco[i]); // in this case: x=eta, y=SF
                }

                // -- Rearrangement in order of x (eta) values -- //
                Double_t etamin; // the minimum value in given iteration
                Double_t etalow = -9999; // the minimum value of the previous iteration
                for(Int_t j=0; j<nEtaBins_reco; j++)
                {
                        Int_t jj = -9999;

                        etamin = 9999;
                        for(Int_t k=0; k<nEtaBins_reco; k++)
                        {
                                if(etalow < eta_reco[k] && eta_reco[k] < etamin) // the lowest number, but higher than the one from the previous iteration
                                {
                                        jj = k;
                                        etamin = eta_reco[k];
                                }
                        }
                        eta_reco_ordered[j] = eta_reco[jj];
                        SF_reco_ordered[j] = SF_reco[jj];

                        etalow = etamin;
                } // End of rearrangement

                for(Int_t iter_x = 0; iter_x < nEtaBins_reco; iter_x++)
                {
                        Eff_Reco_data_GtoH[iter_x][iter_y] = SF_reco_ordered[iter_x]; // actually, it is the scale factor.
                        Eff_Reco_MC_GtoH[iter_x][iter_y] = 1;
                }
        }

        // ID
        TFile *f_ID = new TFile(Location+"ID_SF_RunGtoH.root");
//	TH2F *h_ID_data = (TH2F*)f_ID->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA");
//	TH2F *h_ID_MC = (TH2F*)f_ID->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesMC/abseta_pair_ne_MC");
        TH2F *h_ID_data = (TH2F*)f_ID->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/abseta_pt_DATA");
        TH2F *h_ID_MC = (TH2F*)f_ID->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");

        Int_t nEtaBins_ID = h_ID_data->GetNbinsX();
        Int_t nPtBins_ID = h_ID_data->GetNbinsY();

        for(Int_t iter_x = 0; iter_x < nEtaBins_ID; iter_x++)
        {
                for(Int_t iter_y = 0; iter_y < nPtBins_ID; iter_y++)
                {
                        Int_t i_etabin = iter_x + 1;
                        Int_t i_ptbin = iter_y + 1;

                        Double_t ID_data = h_ID_data->GetBinContent(i_etabin, i_ptbin);
                        Double_t ID_MC = h_ID_MC->GetBinContent(i_etabin, i_ptbin);

                        Eff_ID_data_GtoH[iter_x][iter_y] = ID_data;
                        Eff_ID_MC_GtoH[iter_x][iter_y] = ID_MC;
                }
        }

        // ISO
        TFile *f_iso = new TFile(Location+"ISO_SF_RunGtoH.root");
//	TH2F *h_Iso_data = (TH2F*)f_iso->Get("tkLooseISO_highptID_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA");
//	TH2F *h_Iso_MC = (TH2F*)f_iso->Get("tkLooseISO_highptID_newpt_eta/efficienciesMC/abseta_pair_ne_MC");
        TH2F *h_Iso_data = (TH2F*)f_iso->Get("TightISO_TightID_pt_eta/efficienciesDATA/abseta_pt_DATA");
        TH2F *h_Iso_MC = (TH2F*)f_iso->Get("TightISO_TightID_pt_eta/efficienciesMC/abseta_pt_MC");

        Int_t nEtaBins_iso = h_Iso_data->GetNbinsX();
        Int_t nPtBins_iso = h_Iso_data->GetNbinsY();

        for(Int_t iter_x = 0; iter_x < nEtaBins_iso; iter_x++)
        {
                for(Int_t iter_y = 0; iter_y < nPtBins_iso; iter_y++)
                {
                        Int_t i_etabin = iter_x + 1;
                        Int_t i_ptbin = iter_y + 1;

                        Double_t Iso_data = h_Iso_data->GetBinContent(i_etabin, i_ptbin);
                        Double_t Iso_MC = h_Iso_MC->GetBinContent(i_etabin, i_ptbin);

                        Eff_Iso_data_GtoH[iter_x][iter_y] = Iso_data;
                        Eff_Iso_MC_GtoH[iter_x][iter_y] = Iso_MC;
                }
        }

        // TRIGGER
        TFile *f_HLT = new TFile(Location+"Trigger_SF_RunGtoH.root");
        TH2F *h_HLT_data = (TH2F*)f_HLT->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
        TH2F *h_HLT_MC = (TH2F*)f_HLT->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/abseta_pt_MC");

        Int_t nEtaBins_HLT = h_HLT_data->GetNbinsX();
        Int_t nPtBins_HLT = h_HLT_data->GetNbinsY();

        for(Int_t iter_x = 0; iter_x < nEtaBins_HLT; iter_x++)
        {
                for(Int_t iter_y = 0; iter_y < nPtBins_HLT; iter_y++)
                {
                        Int_t i_etabin = iter_x + 1;
                        Int_t i_ptbin = iter_y + 1;

                        Double_t HLT_data = h_HLT_data->GetBinContent(i_etabin, i_ptbin);
                        Double_t HLT_MC = h_HLT_MC->GetBinContent(i_etabin, i_ptbin);

                        Eff_HLT_data_GtoH[iter_x][iter_y] = HLT_data;
                        Eff_HLT_MC_GtoH[iter_x][iter_y] = HLT_MC;
                }
        }
        std::cout << "Setting for efficiency correction factors (GtoH) is completed" << endl;
} // End of SetupEfficiencyScaleFactor_GtoH

void DYAnalyzer::SetupEfficiencyScaleFactor_GtoH_new()
{
    TString Location = "./etc/effSF/effSF_muon/";
    TFile *f1 = new TFile(Location+"New_SF_RunGtoH.root");
    std::cout << "[Tag&Probe efficiency is from " << Location+"New_SF_RunGtoH.root" << "]" << endl;


//    // RECO
//    TFile *f_reco = new TFile(Location+"Tracking_SF_RunBtoF.root");
//    TGraphAsymmErrors *h_Reco_ratio = (TGraphAsymmErrors*)f_reco->Get("ratio_eff_eta3_dr030e030_corr");
//    Int_t nEtaBins_reco = h_Reco_ratio->GetN();
//    Int_t nPtBins_reco = 1;

//    Double_t eta_reco[15]; Double_t eta_reco_ordered[15];
//    Double_t SF_reco[15]; Double_t SF_reco_ordered[15];

//    for(Int_t iter_y = 0; iter_y < nPtBins_reco; iter_y++) // in this case (2d distribution): x=eta, y=pT
//    {
//            for(Int_t i=0; i<nEtaBins_reco; i++)
//            {
//                    h_Reco_ratio->GetPoint(i, eta_reco[i], SF_reco[i]); // in this case: x=eta, y=SF
//            }

//            // -- Rearrangement in order of x (eta) values -- //
//            Double_t etamin; // the minimum value in given iteration
//            Double_t etalow = -9999; // the minimum value of the previous iteration
//            for(Int_t j=0; j<nEtaBins_reco; j++)
//            {
//                    Int_t jj = -9999;

//                    etamin = 9999;
//                    for(Int_t k=0; k<nEtaBins_reco; k++)
//                    {
//                            if(etalow < eta_reco[k] && eta_reco[k] < etamin) // the lowest number but higher than the one from previos iteration
//                            {
//                                    jj = k;
//                                    etamin = eta_reco[k];
//                            }
//                    }
//                    eta_reco_ordered[j] = eta_reco[jj];
//                    SF_reco_ordered[j] = SF_reco[jj];

//                    etalow = etamin;
//            } // End of rearrangement

//            for(Int_t iter_x = 0; iter_x < nEtaBins_reco; iter_x++)
//            {
//                    Eff_Reco_data_BtoF[iter_x][iter_y] = SF_reco_ordered[iter_x]; // actually, it is the scale factor.
//                    Eff_Reco_MC_BtoF[iter_x][iter_y] = 1;
//            }
//    }

    // ID
    TGraphAsymmErrors *h_data_lead_ID_eff[4];
    TGraphAsymmErrors *h_data_sublead_ID_eff[4];
    TGraphAsymmErrors *h_mc_lead_ID_eff[4];
    TGraphAsymmErrors *h_mc_sublead_ID_eff[4];

    h_data_lead_ID_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Leading_pteta_abseta0");
    h_data_lead_ID_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Leading_pteta_abseta1");
    h_data_lead_ID_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Leading_pteta_abseta2");
    h_data_lead_ID_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Leading_pteta_abseta3");

    h_data_sublead_ID_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Subleading_pteta_abseta0");
    h_data_sublead_ID_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Subleading_pteta_abseta1");
    h_data_sublead_ID_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Subleading_pteta_abseta2");
    h_data_sublead_ID_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_Tight2012_from_Subleading_pteta_abseta3");

    h_mc_lead_ID_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Leading_pteta_abseta0");
    h_mc_lead_ID_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Leading_pteta_abseta1");
    h_mc_lead_ID_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Leading_pteta_abseta2");
    h_mc_lead_ID_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Leading_pteta_abseta3");

    h_mc_sublead_ID_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Subleading_pteta_abseta0");
    h_mc_sublead_ID_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Subleading_pteta_abseta1");
    h_mc_sublead_ID_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Subleading_pteta_abseta2");
    h_mc_sublead_ID_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_Tight2012_from_Subleading_pteta_abseta3");

    Int_t nEtaBins_ID = 4;
    Int_t nPtBins_ID = 6;

    Double_t x_data_lead_ID[6];
    Double_t y_data_lead_ID[6];

    Double_t x_data_sublead_ID[6];
    Double_t y_data_sublead_ID[6];

    Double_t x_mc_lead_ID[6];
    Double_t y_mc_lead_ID[6];

    Double_t x_mc_sublead_ID[6];
    Double_t y_mc_sublead_ID[6];

    TGraphAsymmErrors *h_Data_lead_ID_eff;
    TGraphAsymmErrors *h_Data_sublead_ID_eff;
    TGraphAsymmErrors *h_MC_lead_ID_eff;
    TGraphAsymmErrors *h_MC_sublead_ID_eff;

    for(Int_t iter_x = 0; iter_x < nEtaBins_ID; iter_x++)
    {
        h_Data_lead_ID_eff = (TGraphAsymmErrors*)h_data_lead_ID_eff[iter_x]->Clone();
        h_Data_sublead_ID_eff = (TGraphAsymmErrors*)h_data_sublead_ID_eff[iter_x]->Clone();
        h_MC_lead_ID_eff = (TGraphAsymmErrors*)h_mc_lead_ID_eff[iter_x]->Clone();
        h_MC_sublead_ID_eff = (TGraphAsymmErrors*)h_mc_sublead_ID_eff[iter_x]->Clone();

        for(Int_t i = 0; i < nPtBins_ID; i++)
        {
            h_Data_sublead_ID_eff->GetPoint(i, x_data_sublead_ID[i], y_data_sublead_ID[i]);
            h_MC_sublead_ID_eff->GetPoint(i, x_mc_sublead_ID[i], y_mc_sublead_ID[i]);

            if(i == 0)
            {
                //It is just to initialize. They don't work in SF
                x_data_lead_ID[0] = 1;
                y_data_lead_ID[0] = 1;

                x_mc_lead_ID[0] = 1;
                y_mc_lead_ID[0] = 1;

                continue;
            }

            h_Data_lead_ID_eff->GetPoint(i-1, x_data_lead_ID[i], y_data_lead_ID[i]);
            h_MC_lead_ID_eff->GetPoint(i-1, x_mc_lead_ID[i], y_mc_lead_ID[i]);
        }

        for(Int_t iter_y = 0; iter_y < nPtBins_ID; iter_y++)
        {
            Eff_ID_data_GtoH_lead[iter_x][iter_y] = y_data_lead_ID[iter_y];
            Eff_ID_MC_GtoH_lead[iter_x][iter_y] = y_mc_lead_ID[iter_y];

            Eff_ID_data_GtoH_sublead[iter_x][iter_y] = y_data_sublead_ID[iter_y];
            Eff_ID_MC_GtoH_sublead[iter_x][iter_y] = y_mc_sublead_ID[iter_y];
        }
    }

    // ISO
    TGraphAsymmErrors *h_data_lead_iso_eff[4];
    TGraphAsymmErrors *h_data_sublead_iso_eff[4];
    TGraphAsymmErrors *h_mc_lead_iso_eff[4];
    TGraphAsymmErrors *h_mc_sublead_iso_eff[4];

    h_data_lead_iso_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Leading_pteta_abseta0");
    h_data_lead_iso_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Leading_pteta_abseta1");
    h_data_lead_iso_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Leading_pteta_abseta2");
    h_data_lead_iso_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Leading_pteta_abseta3");

    h_data_sublead_iso_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta0");
    h_data_sublead_iso_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta1");
    h_data_sublead_iso_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta2");
    h_data_sublead_iso_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta3");

    h_mc_lead_iso_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Leading_pteta_abseta0");
    h_mc_lead_iso_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Leading_pteta_abseta1");
    h_mc_lead_iso_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Leading_pteta_abseta2");
    h_mc_lead_iso_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Leading_pteta_abseta3");

    h_mc_sublead_iso_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta0");
    h_mc_sublead_iso_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta1");
    h_mc_sublead_iso_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta2");
    h_mc_sublead_iso_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_dBeta_015_from_Tight2012_and_Subleading_pteta_abseta3");

    Int_t nEtaBins_iso = 4;
    Int_t nPtBins_iso = 6;

    Double_t x_data_lead_iso[6];
    Double_t y_data_lead_iso[6];

    Double_t x_data_sublead_iso[6];
    Double_t y_data_sublead_iso[6];

    Double_t x_mc_lead_iso[6];
    Double_t y_mc_lead_iso[6];

    Double_t x_mc_sublead_iso[6];
    Double_t y_mc_sublead_iso[6];

    TGraphAsymmErrors *h_Data_lead_iso_eff;
    TGraphAsymmErrors *h_Data_sublead_iso_eff;
    TGraphAsymmErrors *h_MC_lead_iso_eff;
    TGraphAsymmErrors *h_MC_sublead_iso_eff;

    for(Int_t iter_x = 0; iter_x < nEtaBins_iso; iter_x++)
    {
        h_Data_lead_iso_eff = (TGraphAsymmErrors*)h_data_lead_iso_eff[iter_x]->Clone();
        h_Data_sublead_iso_eff = (TGraphAsymmErrors*)h_data_sublead_iso_eff[iter_x]->Clone();
        h_MC_lead_iso_eff = (TGraphAsymmErrors*)h_mc_lead_iso_eff[iter_x]->Clone();
        h_MC_sublead_iso_eff = (TGraphAsymmErrors*)h_mc_sublead_iso_eff[iter_x]->Clone();

        for(Int_t i=0; i<nPtBins_iso; i++) //0 to 5
        {
            h_Data_sublead_iso_eff->GetPoint(i, x_data_sublead_iso[i], y_data_sublead_iso[i]);
            h_MC_sublead_iso_eff->GetPoint(i, x_mc_sublead_iso[i], y_mc_sublead_iso[i]);

            if(i == 0)
            {
                //It is just to initialize. They don't work in SF
                x_data_lead_iso[0] = 1;
                y_data_lead_iso[0] = 1;

                x_mc_lead_iso[0] = 1;
                y_mc_lead_iso[0] = 1;

                continue;
            }

            h_Data_lead_iso_eff->GetPoint(i-1, x_data_lead_iso[i], y_data_lead_iso[i]);
            h_MC_lead_iso_eff->GetPoint(i-1, x_mc_lead_iso[i], y_mc_lead_iso[i]);
        }

        for(Int_t iter_y = 0; iter_y < nPtBins_iso; iter_y++)
        {
            Eff_Iso_data_GtoH_lead[iter_x][iter_y] = y_data_lead_iso[iter_y];
            Eff_Iso_MC_GtoH_lead[iter_x][iter_y] = y_mc_lead_iso[iter_y];

            Eff_Iso_data_GtoH_sublead[iter_x][iter_y] = y_data_sublead_iso[iter_y];
            Eff_Iso_MC_GtoH_sublead[iter_x][iter_y] = y_mc_sublead_iso[iter_y];
        }
    }

    // TRIGGER
    TGraphAsymmErrors *h_data_lead_trig_eff[4];
    TGraphAsymmErrors *h_data_sublead_trig_eff[4];
    TGraphAsymmErrors *h_mc_lead_trig_eff[4];
    TGraphAsymmErrors *h_mc_sublead_trig_eff[4];

    h_data_lead_trig_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Leading_pteta_abseta0");
    h_data_lead_trig_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Leading_pteta_abseta1");
    h_data_lead_trig_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Leading_pteta_abseta2");
    h_data_lead_trig_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Leading_pteta_abseta3");

    h_data_sublead_trig_eff[0] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Subleading_pteta_abseta0");
    h_data_sublead_trig_eff[1] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Subleading_pteta_abseta1");
    h_data_sublead_trig_eff[2] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Subleading_pteta_abseta2");
    h_data_sublead_trig_eff[3] = (TGraphAsymmErrors*)f1->Get("Data_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Subleading_pteta_abseta3");

    h_mc_lead_trig_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Leading_pteta_abseta0");
    h_mc_lead_trig_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Leading_pteta_abseta1");
    h_mc_lead_trig_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Leading_pteta_abseta2");
    h_mc_lead_trig_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Leading_pteta_abseta3");

    h_mc_sublead_trig_eff[0] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Subleading_pteta_abseta0");
    h_mc_sublead_trig_eff[1] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Subleading_pteta_abseta1");
    h_mc_sublead_trig_eff[2] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Subleading_pteta_abseta2");
    h_mc_sublead_trig_eff[3] = (TGraphAsymmErrors*)f1->Get("MC_weight_IsoMu24_OR_IsoTkMu24_from_Tight2012_and_dBeta_015_and_Subleading_pteta_abseta3");

    Int_t nEtaBins_HLT = 4;
    Int_t nPtBins_HLT = 8;

    Double_t x_data_lead_trig[8];
    Double_t y_data_lead_trig[8];

    Double_t x_data_sublead_trig[8];
    Double_t y_data_sublead_trig[8];

    Double_t x_mc_lead_trig[8];
    Double_t y_mc_lead_trig[8];

    Double_t x_mc_sublead_trig[8];
    Double_t y_mc_sublead_trig[8];

    TGraphAsymmErrors *h_Data_lead_trig_eff;
    TGraphAsymmErrors *h_Data_sublead_trig_eff;
    TGraphAsymmErrors *h_MC_lead_trig_eff;
    TGraphAsymmErrors *h_MC_sublead_trig_eff;

    for(Int_t iter_x = 0; iter_x < nEtaBins_HLT; iter_x++)
    {
        h_Data_lead_trig_eff = (TGraphAsymmErrors*)h_data_lead_trig_eff[iter_x]->Clone();
        h_Data_sublead_trig_eff = (TGraphAsymmErrors*)h_data_sublead_trig_eff[iter_x]->Clone();
        h_MC_lead_trig_eff = (TGraphAsymmErrors*)h_mc_lead_trig_eff[iter_x]->Clone();
        h_MC_sublead_trig_eff = (TGraphAsymmErrors*)h_mc_sublead_trig_eff[iter_x]->Clone();

        for(Int_t i=0; i<nPtBins_HLT; i++)
        {
            h_Data_lead_trig_eff->GetPoint(i, x_data_lead_trig[i], y_data_lead_trig[i]);
            h_MC_lead_trig_eff->GetPoint(i, x_mc_lead_trig[i], y_mc_lead_trig[i]);

            if(iter_x == 2 && i == 7)
            {
                h_Data_sublead_trig_eff->GetPoint(i, x_data_sublead_trig[i], y_data_sublead_trig[i]);

                x_mc_sublead_trig[i] = 1;
                y_mc_sublead_trig[i] = y_data_sublead_trig[i-1];

                continue;
            }
            else if(iter_x == 3 && i >= 6)
            {
                x_data_sublead_trig[i] = 1;
                y_data_sublead_trig[i] = y_data_sublead_trig[i-1];

                x_mc_sublead_trig[i] = 1;
                y_mc_sublead_trig[i] = y_data_sublead_trig[i-1];

                continue;
            }

            h_Data_sublead_trig_eff->GetPoint(i, x_data_sublead_trig[i], y_data_sublead_trig[i]);
            h_MC_sublead_trig_eff->GetPoint(i, x_mc_sublead_trig[i], y_mc_sublead_trig[i]);
        }

        for(Int_t iter_y = 0; iter_y < nPtBins_HLT; iter_y++)
        {
            Eff_HLT_data_GtoH_lead[iter_x][iter_y] = y_data_lead_trig[iter_y];
            Eff_HLT_MC_GtoH_lead[iter_x][iter_y] = y_mc_lead_trig[iter_y];

            Eff_HLT_data_GtoH_sublead[iter_x][iter_y] = y_data_sublead_trig[iter_y];
            Eff_HLT_MC_GtoH_sublead[iter_x][iter_y] = y_mc_sublead_trig[iter_y];
        }
    }
    std::cout << "Setting for efficiency correction factors (BtoF) is completed" << endl;
} // End of SetupEfficiencyScaleFactor_GtoH_new

void DYAnalyzer::SetupEfficiencyScaleFactor_electron()
{
        TString Location = "./etc/effSF/effSF_electron/";
        std::cout << "[Tag&Probe efficiency is from " << Location+"*.root" << "]" << endl;

        // RECO
        TFile *f_reco = new TFile(Location+"Reco_SF.root");
        TGraphErrors *h_reco_sf = (TGraphErrors*)f_reco->Get("grSF1D_0");

        Int_t nEtaBins_reco = h_reco_sf->GetN();
        Int_t nPtBins_reco = 1;

        Double_t eta_reco[30]; Double_t eta_reco_ordered[30];
        Double_t SF_reco[30]; Double_t SF_reco_ordered[30];

        for(Int_t iter_y = 0; iter_y < nPtBins_reco; iter_y++) // in this case (2d distribution): x=eta, y=pT
        {
            for(Int_t i=0; i<nEtaBins_reco; i++)
            {
                h_reco_sf->GetPoint(i, eta_reco[i], SF_reco[i]); // in this case: x=eta, y=SF
            }

            // -- Rearrangement in order of x (eta) values -- //
            Double_t etamin; // the minimum value in given iteration
            Double_t etalow = -9999; // the minimum value of the previous iteration
            for(Int_t j=0; j<nEtaBins_reco; j++)
            {
                Int_t jj = -9999;

                etamin = 9999;
                for(Int_t k=0; k<nEtaBins_reco; k++)
                {
                    if(etalow < eta_reco[k] && eta_reco[k] < etamin) // the lowest number, but higher than the one from the previous iteration
                    {
                        jj = k;
                        etamin = eta_reco[k];
                    }
                }
                eta_reco_ordered[j] = eta_reco[jj];
                SF_reco_ordered[j] = SF_reco[jj];

                etalow = etamin;
//                std::cout << j << "  " << eta_reco_ordered[j] << "  " << SF_reco_ordered[j] << endl;
            } // End of rearrangement

            for(Int_t iter_x = 0; iter_x < nEtaBins_reco; iter_x++)
            {
                Eff_Reco_data[iter_x][iter_y] = SF_reco_ordered[iter_x]; // actually, it is the scale factor.
                Eff_Reco_MC[iter_x][iter_y] = 1;
//                cout << "Reco: bin# = [" << iter_x << ", " << iter_y << "]" << " eta = " << eta_reco_ordered[iter_x] << " sf = " << SF_reco_ordered[iter_x] << endl;
            }
        }

        // ID
//	TFile *f_ID = new TFile(Location+"MediumID_SF.root");
        TFile *f_ID = new TFile(Location+"Electron_MediumID_Run2016BtoH.root");
        TGraphErrors *h_id_sf_0 = (TGraphErrors*)f_ID->Get("grSF1D_0"); // Five graphs (functions of eta) for different pT bins
        TGraphErrors *h_id_sf_1 = (TGraphErrors*)f_ID->Get("grSF1D_1");
        TGraphErrors *h_id_sf_2 = (TGraphErrors*)f_ID->Get("grSF1D_2");
        TGraphErrors *h_id_sf_3 = (TGraphErrors*)f_ID->Get("grSF1D_3");
        TGraphErrors *h_id_sf_4 = (TGraphErrors*)f_ID->Get("grSF1D_4");

	Int_t nEtaBins_id = h_id_sf_0->GetN();
	Int_t nPtBins_id = 5;

//	Double_t eta_id[10]; Double_t eta_id_ordered[10];
//	Double_t SF_id[10]; Double_t SF_id_ordered[10];
        Double_t eta_id[20]; Double_t eta_id_ordered[20];
        Double_t SF_id[20]; Double_t SF_id_ordered[20];

	TGraphErrors *h_id_sf;
        for(Int_t iter_y = 0; iter_y < nPtBins_id; iter_y++) // In this case, x=eta, y=pT
	{
		if(iter_y == 0) h_id_sf = (TGraphErrors*)h_id_sf_0->Clone();
		else if(iter_y == 1) h_id_sf = (TGraphErrors*)h_id_sf_1->Clone();
		else if(iter_y == 2) h_id_sf = (TGraphErrors*)h_id_sf_2->Clone();
		else if(iter_y == 3) h_id_sf = (TGraphErrors*)h_id_sf_3->Clone();
		else if(iter_y == 4) h_id_sf = (TGraphErrors*)h_id_sf_4->Clone();

		for(Int_t i=0; i<nEtaBins_id; i++)
		{
                        h_id_sf->GetPoint(i, eta_id[i], SF_id[i]); // In this case, x=eta, y=SF
		}

                // -- Rearrangement in order of x (eta) values -- //
                Double_t etamin; // the minimum value in given iteration
                Double_t etalow = -9999; // the minimum value of the previous iteration
		for(Int_t j=0; j<nEtaBins_id; j++)
		{
			Int_t jj = -9999;

                        etamin = 9999;
			for(Int_t k=0; k<nEtaBins_id; k++)
			{
                                if(etalow < eta_id[k] && eta_id[k] < etamin) // the lowest number, but higher than the one from the previous iteration
				{
					jj = k;
                                        etamin = eta_id[k];
				}
			}
                        eta_id_ordered[j] = eta_id[jj];
                        SF_id_ordered[j] = SF_id[jj];

                        etalow = etamin;
//			std::cout << j << "  " << eta_id_ordered[j] << "  " << SF_id_ordered[j] << endl;
		} // End of rearrangement

		for(Int_t iter_x = 0; iter_x < nEtaBins_id; iter_x++)
		{
                        Eff_ID_data[iter_x][iter_y] = SF_id_ordered[iter_x]; // actually, it is the scale factor.
			Eff_ID_MC[iter_x][iter_y] = 1;
//			std::cout << "ID: bin# = [" << iter_x << ", " << iter_y << "]" << " eta = " << eta_id_ordered[iter_x] << " sf = " << SF_id_ordered[iter_x] << endl;
		}
	}

        // TRIGGER (leg2)
        //TFile *f_HLT = new TFile(Location+"Leg2_SF.root");
        TFile *f_HLT = new TFile(Location+"Electron_Leg2_SF.root");
        TGraphErrors *h_leg2_sf_0 = (TGraphErrors*)f_HLT->Get("grSF1D_0"); // 7 graphs (functions of eta) for different pT bins
        TGraphErrors *h_leg2_sf_1 = (TGraphErrors*)f_HLT->Get("grSF1D_1");
        TGraphErrors *h_leg2_sf_2 = (TGraphErrors*)f_HLT->Get("grSF1D_2");
        TGraphErrors *h_leg2_sf_3 = (TGraphErrors*)f_HLT->Get("grSF1D_3");
        TGraphErrors *h_leg2_sf_4 = (TGraphErrors*)f_HLT->Get("grSF1D_4");
        TGraphErrors *h_leg2_sf_5 = (TGraphErrors*)f_HLT->Get("grSF1D_5");
        TGraphErrors *h_leg2_sf_6 = (TGraphErrors*)f_HLT->Get("grSF1D_6");
        //TGraphErrors *h_leg2_sf_7 = (TGraphErrors*)f_HLT->Get("grSF1D_7");

        Int_t nEtaBins_leg2 = h_leg2_sf_0->GetN();
        //Int_t nPtBins_leg2 = 8;
        Int_t nPtBins_leg2 = 7;

        //Double_t eta_leg2[10]; Double_t eta_leg2_ordered[10];
        //Double_t SF_leg2[10]; Double_t SF_leg2_ordered[10];
        Double_t eta_leg2[20]; Double_t eta_leg2_ordered[20];
        Double_t SF_leg2[20]; Double_t SF_leg2_ordered[20];

        TGraphErrors *h_leg2_sf;
        for(Int_t iter_y = 0; iter_y < nPtBins_leg2; iter_y++) // In this case, x=eta, y=pT
        {
            if(iter_y == 0) h_leg2_sf = (TGraphErrors*)h_leg2_sf_0->Clone();
            else if(iter_y == 1) h_leg2_sf = (TGraphErrors*)h_leg2_sf_1->Clone();
            else if(iter_y == 2) h_leg2_sf = (TGraphErrors*)h_leg2_sf_2->Clone();
            else if(iter_y == 3) h_leg2_sf = (TGraphErrors*)h_leg2_sf_3->Clone();
            else if(iter_y == 4) h_leg2_sf = (TGraphErrors*)h_leg2_sf_4->Clone();
            else if(iter_y == 5) h_leg2_sf = (TGraphErrors*)h_leg2_sf_5->Clone();
            else if(iter_y == 6) h_leg2_sf = (TGraphErrors*)h_leg2_sf_6->Clone();
            //else if(iter_y == 7) h_leg2_sf = (TGraphErrors*)h_leg2_sf_7->Clone();

            for(Int_t i=0; i<nEtaBins_leg2; i++)
            {
                h_leg2_sf->GetPoint(i, eta_leg2[i], SF_leg2[i]); // In this case, x=eta, y=SF
            }

            // -- Rearrangement in order of x (eta) values -- //
            Double_t etamin; // the minimum value in given iteration
            Double_t etalow = -9999; // the minimum value of the previous iteration
            for(Int_t j=0; j<nEtaBins_leg2; j++)
            {
                Int_t jj = -9999;

                etamin = 9999;
                for(Int_t k=0; k<nEtaBins_leg2; k++)
                {
                    if(etalow < eta_leg2[k] && eta_leg2[k] < etamin) // the lowest number, but higher than the one from the previous iteration
                    {
                            jj = k;
                            etamin = eta_leg2[k];
                    }
                }
                eta_leg2_ordered[j] = eta_leg2[jj];
                SF_leg2_ordered[j] = SF_leg2[jj];

                etalow = etamin;
                //cout << j << "  " << eta_leg2_ordered[j] << "  " << SF_leg2_ordered[j] << endl;
            } // End of rearrangement

            for(Int_t iter_x = 0; iter_x < nEtaBins_leg2; iter_x++)
            {
                Eff_HLT_Leg2_data[iter_x][iter_y] = SF_leg2_ordered[iter_x]; // actually, it is the scale factor.
                Eff_HLT_Leg2_MC[iter_x][iter_y] = 1;
                //cout << "ID: bin# = [" << iter_x << ", " << iter_y << "]" << " eta = " << eta_leg2_ordered[iter_x] << " sf = " << SF_leg2_ordered[iter_x] << endl;
            }
        }

        std::cout << "Setting for efficiency correction factors is completed" << endl;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_BtoF(Muon mu1, Muon mu2)
{
    Double_t weight = -999;

    // -- Muon1 -- //
    Double_t Pt1 = mu1.Pt;
    Double_t eta1 = mu1.eta;

    Int_t ptbin1_Reco = Find_muon_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_muon_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    if(ptbin1_Reco == 9999 || etabin1_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Reco, etabin1_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Reco, etabin1_Reco, Pt1, eta1);
        return -999;
    }
    if(ptbin1_ID == 9999 || etabin1_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_ID, etabin1_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_ID, etabin1_ID, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Iso == 9999 || etabin1_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Iso, etabin1_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Iso, etabin1_Iso, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Trig == 9999 || etabin1_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Trig, etabin1_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Trig, etabin1_Trig, Pt1, eta1);
        return -999;
    }

    Double_t Eff_muon1_data = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso];

    // -- Muon2 -- //
    Double_t Pt2 = mu2.Pt;
    Double_t eta2 = mu2.eta;

    Int_t ptbin2_Reco = Find_muon_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_muon_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_muon_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_muon_EtaBin_ID(eta2);

    Int_t ptbin2_Iso = Find_muon_PtBin_Iso(Pt2);
    Int_t etabin2_Iso = Find_muon_EtaBin_Iso(eta2);

    Int_t ptbin2_Trig = Find_muon_PtBin_Trig(Pt2);
    Int_t etabin2_Trig = Find_muon_EtaBin_Trig(eta2);

    if(ptbin2_Reco == 9999 || etabin2_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Reco, etabin2_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Reco, etabin2_Reco, Pt2, eta2);
        return -999;
    }
    if(ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_ID, etabin2_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_ID, etabin2_ID, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Iso == 9999 || etabin2_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Iso, etabin2_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Iso, etabin2_Iso, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Trig == 9999 || etabin2_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Trig, etabin2_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Trig, etabin2_Trig, Pt2, eta2);
        return -999;
    }

    Double_t Eff_muon2_data = Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso];
    Double_t Eff_muon2_MC = Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso];

    // Trigger SF
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_data = Eff_HLT_data_BtoF[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_BtoF[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

    Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

    //std::cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2)
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

        printf("[Data]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID],
               Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID],
               Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

        printf("[MC]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID],
               Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID],
               Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

        printf("[SF] Weight = %.3lf\n", weight);
    }
    return weight;

}// End of EfficiencySF_EventWeight_HLT_BtoF(mu1, mu2)


Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_BtoF(TLorentzVector mu1, TLorentzVector mu2)
{
    Double_t weight = -999;

    // -- Muon1 -- //
    Double_t Pt1 = mu1.Pt();
    Double_t eta1 = mu1.Eta();

    Int_t ptbin1_Reco = Find_muon_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_muon_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    if(ptbin1_Reco == 9999 || etabin1_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Reco, etabin1_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Reco, etabin1_Reco, Pt1, eta1);
        return -999;
    }
    if(ptbin1_ID == 9999 || etabin1_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_ID, etabin1_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_ID, etabin1_ID, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Iso == 9999 || etabin1_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Iso, etabin1_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Iso, etabin1_Iso, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Trig == 9999 || etabin1_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Trig, etabin1_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Trig, etabin1_Trig, Pt1, eta1);
        return -999;
    }

    Double_t Eff_muon1_data = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso];

    // -- Muon2 -- //
    Double_t Pt2 = mu2.Pt();
    Double_t eta2 = mu2.Eta();

    Int_t ptbin2_Reco = Find_muon_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_muon_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_muon_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_muon_EtaBin_ID(eta2);

    Int_t ptbin2_Iso = Find_muon_PtBin_Iso(Pt2);
    Int_t etabin2_Iso = Find_muon_EtaBin_Iso(eta2);

    Int_t ptbin2_Trig = Find_muon_PtBin_Trig(Pt2);
    Int_t etabin2_Trig = Find_muon_EtaBin_Trig(eta2);

    if(ptbin2_Reco == 9999 || etabin2_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Reco, etabin2_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Reco, etabin2_Reco, Pt2, eta2);
        return -999;
    }
    if(ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_ID, etabin2_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_ID, etabin2_ID, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Iso == 9999 || etabin2_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Iso, etabin2_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Iso, etabin2_Iso, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Trig == 9999 || etabin2_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Trig, etabin2_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Trig, etabin2_Trig, Pt2, eta2);
        return -999;
    }

    Double_t Eff_muon2_data = Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso];
    Double_t Eff_muon2_MC = Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso];

    // Trigger SF
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_data = Eff_HLT_data_BtoF[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_BtoF[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

    Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

    //std::cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2)
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

        printf("[Data]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID],
               Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID],
               Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

        printf("[MC]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID],
               Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID],
               Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

        printf("[SF] Weight = %.3lf\n", weight);
    }
    return weight;

}// End of EfficiencySF_EventWeight_HLT_BtoF(mu1_4, mu2_4)


Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_BtoF_new(SelectedMuMu_t *MuMu)
{
    Double_t weight = -999;
    Double_t Pt1 = 0;
    Double_t eta1 = -9999;
    Double_t Pt2 = 0;
    Double_t eta2 = -9999;

    if (MuMu->Muon_pT->at(0) >= MuMu->Muon_pT->at(1))
    {
        Pt1 = MuMu->Muon_pT->at(0);
        eta1 = MuMu->Muon_eta->at(0);
        Pt2 = MuMu->Muon_pT->at(1);
        eta2 = MuMu->Muon_eta->at(1);
    }
    else
    {
        Pt1 = MuMu->Muon_pT->at(1);
        eta1 = MuMu->Muon_eta->at(1);
        Pt2 = MuMu->Muon_pT->at(0);
        eta2 = MuMu->Muon_eta->at(0);
    }

    // -- Muon1 LEADING -- //
    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig_new(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    if(ptbin1_ID == 9999 || etabin1_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_ID, etabin1_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_ID, etabin1_ID, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Iso == 9999 || etabin1_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Iso, etabin1_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Iso, etabin1_Iso, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Trig == 9999 || etabin1_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Trig, etabin1_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Trig, etabin1_Trig, Pt1, eta1);
        return -999;
    }

    Double_t Eff_muon1_data = Eff_ID_data_BtoF_lead[etabin1_ID][ptbin1_ID] * Eff_Iso_data_BtoF_lead[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_ID_MC_BtoF_lead[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_BtoF_lead[etabin1_Iso][ptbin1_Iso];

    // -- Muon2 SUBLEADING -- //
    Int_t ptbin2_ID = Find_muon_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_muon_EtaBin_ID(eta2);

    Int_t ptbin2_Iso = Find_muon_PtBin_Iso(Pt2);
    Int_t etabin2_Iso = Find_muon_EtaBin_Iso(eta2);

    Int_t ptbin2_Trig = Find_muon_PtBin_Trig_new(Pt2);
    Int_t etabin2_Trig = Find_muon_EtaBin_Trig(eta2);

    if(ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_ID, etabin2_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_ID, etabin2_ID, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Iso == 9999 || etabin2_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Iso, etabin2_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Iso, etabin2_Iso, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Trig == 9999 || etabin2_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Trig, etabin2_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Trig, etabin2_Trig, Pt2, eta2);
        return -999;
    }

    Double_t Eff_muon2_data = Eff_ID_data_BtoF_sublead[etabin2_ID][ptbin2_ID] * Eff_Iso_data_BtoF_sublead[etabin2_Iso][ptbin2_Iso];
    Double_t Eff_muon2_MC = Eff_ID_MC_BtoF_sublead[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_BtoF_sublead[etabin2_Iso][ptbin2_Iso];

    // Trigger SF
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF_lead[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_data = Eff_HLT_data_BtoF_sublead[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF_lead[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_BtoF_sublead[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

    Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

//    weight = (Eff_ID_data_BtoF_lead[etabin1_ID][ptbin1_ID] * Eff_ID_data_BtoF_sublead[etabin2_ID][ptbin2_ID]) /
//             (Eff_ID_MC_BtoF_lead[etabin1_ID][ptbin1_ID] * Eff_ID_MC_BtoF_sublead[etabin2_ID][ptbin2_ID]); // ONLY ID
//    weight = (Eff_muon1_data * Eff_muon2_data) / (Eff_muon1_MC * Eff_muon2_MC); // ID + ISO
    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2)
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

        printf("[Data]\n");
        printf("\t[Muon1] (ID, Iso): %.3lf, %.3lf)\n", Eff_ID_data_BtoF_lead[etabin1_ID][ptbin1_ID],
               Eff_Iso_data_BtoF_lead[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (ID, Iso): (%.3lf, %.3lf)\n", Eff_ID_data_BtoF_sublead[etabin2_ID][ptbin2_ID],
               Eff_Iso_data_BtoF_sublead[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

        printf("[MC]\n");
        printf("\t[Muon1] (ID, Iso): (%.3lf, %.3lf)\n", Eff_ID_MC_BtoF_lead[etabin1_ID][ptbin1_ID],
               Eff_Iso_MC_BtoF_lead[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (ID, Iso): (%.3lf, %.3lf)\n", Eff_ID_MC_BtoF_sublead[etabin2_ID][ptbin2_ID],
               Eff_Iso_MC_BtoF_sublead[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

        printf("[SF] Weight = %.3lf\n", weight);
    }
    return weight;

}// End of EfficiencySF_EventWeight_HLT_BtoF_new(mu1, mu2)


Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_BtoF(SelectedMuMu_t *MuMu)
{
    Double_t weight = -999;

    // -- Muon1 -- //
//    Double_t Pt1 = MuMu->Muon_TuneP_pT->at(0);
//    Double_t eta1 = MuMu->Muon_TuneP_eta->at(0);
    Double_t Pt1 = MuMu->Muon_pT->at(0);
    Double_t eta1 = MuMu->Muon_eta->at(0);

    Int_t ptbin1_Reco = Find_muon_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_muon_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    if(ptbin1_Reco == 9999 || etabin1_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Reco, etabin1_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Reco, etabin1_Reco, Pt1, eta1);
        cout << "pT1: " << Pt1 << "    eta1: " << eta1 << endl;
        return -999;
    }
    if(ptbin1_ID == 9999 || etabin1_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_ID, etabin1_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_ID, etabin1_ID, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Iso == 9999 || etabin1_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Iso, etabin1_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Iso, etabin1_Iso, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Trig == 9999 || etabin1_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Trig, etabin1_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Trig, etabin1_Trig, Pt1, eta1);
        return -999;
    }

    Double_t Eff_muon1_data = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso];

    // -- Muon2 -- //
//    Double_t Pt2 = MuMu->Muon_TuneP_pT->at(1);
//    Double_t eta2 = MuMu->Muon_TuneP_eta->at(1);
    Double_t Pt2 = MuMu->Muon_pT->at(1);
    Double_t eta2 = MuMu->Muon_eta->at(1);

    Int_t ptbin2_Reco = Find_muon_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_muon_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_muon_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_muon_EtaBin_ID(eta2);

    Int_t ptbin2_Iso = Find_muon_PtBin_Iso(Pt2);
    Int_t etabin2_Iso = Find_muon_EtaBin_Iso(eta2);

    Int_t ptbin2_Trig = Find_muon_PtBin_Trig(Pt2);
    Int_t etabin2_Trig = Find_muon_EtaBin_Trig(eta2);

    if(ptbin2_Reco == 9999 || etabin2_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Reco, etabin2_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Reco, etabin2_Reco, Pt2, eta2);
        return -999;
        cout << "pT2: " << Pt2 << "    eta2: " << eta2 << endl;
    }
    if(ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_ID, etabin2_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_ID, etabin2_ID, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Iso == 9999 || etabin2_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Iso, etabin2_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Iso, etabin2_Iso, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Trig == 9999 || etabin2_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Trig, etabin2_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Trig, etabin2_Trig, Pt2, eta2);
        return -999;
    }

    Double_t Eff_muon2_data = Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso];
    Double_t Eff_muon2_MC = Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso];

    // Trigger SF
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_data = Eff_HLT_data_BtoF[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_BtoF[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

    Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

    //std::cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2)
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

        printf("[Data]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID],
               Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID],
               Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

        printf("[MC]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID],
               Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID],
               Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

        printf("[SF] Weight = %.3lf\n", weight);
    }
    return weight;

}// End of EfficiencySF_EventWeight_HLT_BtoF(SelectedMuMu_t)


Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_GtoH(Muon mu1, Muon mu2)
{
    Double_t weight = -999;

    // -- Muon1 -- //
    Double_t Pt1 = mu1.Pt;
    Double_t eta1 = mu1.eta;

    Int_t ptbin1_Reco = Find_muon_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_muon_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    if(ptbin1_Reco == 9999 || etabin1_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Reco, etabin1_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Reco, etabin1_Reco, Pt1, eta1);
        return -999;
    }
    if(ptbin1_ID == 9999 || etabin1_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_ID, etabin1_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_ID, etabin1_ID, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Iso == 9999 || etabin1_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Iso, etabin1_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Iso, etabin1_Iso, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Trig == 9999 || etabin1_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Trig, etabin1_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Trig, etabin1_Trig, Pt1, eta1);
        return -999;
    }

    Double_t Eff_muon1_data = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso];

    // -- Muon2 -- //
    Double_t Pt2 = mu2.Pt;
    Double_t eta2 = mu2.eta;

    Int_t ptbin2_Reco = Find_muon_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_muon_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_muon_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_muon_EtaBin_ID(eta2);

    Int_t ptbin2_Iso = Find_muon_PtBin_Iso(Pt2);
    Int_t etabin2_Iso = Find_muon_EtaBin_Iso(eta2);

    Int_t ptbin2_Trig = Find_muon_PtBin_Trig(Pt2);
    Int_t etabin2_Trig = Find_muon_EtaBin_Trig(eta2);

    if(ptbin2_Reco == 9999 || etabin2_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Reco, etabin2_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Reco, etabin2_Reco, Pt2, eta2);
        return -999;
    }
    if(ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_ID, etabin2_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_ID, etabin2_ID, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Iso == 9999 || etabin2_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Iso, etabin2_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Iso, etabin2_Iso, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Trig == 9999 || etabin2_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Trig, etabin2_Trig) = (%d, %d, (pt, eta) = (%f, %f))\n", ptbin2_Trig, etabin2_Trig, Pt2, eta2);
        return -999;
    }
    Double_t Eff_muon2_data = Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso];
    Double_t Eff_muon2_MC = Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso];


    // -- Trigger part -- //
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_data = Eff_HLT_data_GtoH[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_GtoH[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

    Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

    //cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2)
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

        printf("[Data]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID],
               Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID],
               Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

        printf("[MC]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID],
               Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID],
               Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

        printf("[SF] Weight = %.3lf\n", weight);
    }
    return weight;

}// End of EfficiencySF_EventWeight_HLT_GtoH(mu1, mu2)


Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_GtoH(TLorentzVector mu1, TLorentzVector mu2)
{
    Double_t weight = -999;

    // -- Muon1 -- //
    Double_t Pt1 = mu1.Pt();
    Double_t eta1 = mu1.Eta();

    Int_t ptbin1_Reco = Find_muon_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_muon_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    if(ptbin1_Reco == 9999 || etabin1_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Reco, etabin1_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Reco, etabin1_Reco, Pt1, eta1);
        return -999;
    }
    if(ptbin1_ID == 9999 || etabin1_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_ID, etabin1_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_ID, etabin1_ID, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Iso == 9999 || etabin1_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Iso, etabin1_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Iso, etabin1_Iso, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Trig == 9999 || etabin1_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Trig, etabin1_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Trig, etabin1_Trig, Pt1, eta1);
        return -999;
    }

    Double_t Eff_muon1_data = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso];

    // -- Muon2 -- //
    Double_t Pt2 = mu2.Pt();
    Double_t eta2 = mu2.Eta();

    Int_t ptbin2_Reco = Find_muon_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_muon_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_muon_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_muon_EtaBin_ID(eta2);

    Int_t ptbin2_Iso = Find_muon_PtBin_Iso(Pt2);
    Int_t etabin2_Iso = Find_muon_EtaBin_Iso(eta2);

    Int_t ptbin2_Trig = Find_muon_PtBin_Trig(Pt2);
    Int_t etabin2_Trig = Find_muon_EtaBin_Trig(eta2);

    if(ptbin2_Reco == 9999 || etabin2_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Reco, etabin2_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Reco, etabin2_Reco, Pt2, eta2);
        return -999;
    }
    if(ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_ID, etabin2_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_ID, etabin2_ID, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Iso == 9999 || etabin2_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Iso, etabin2_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Iso, etabin2_Iso, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Trig == 9999 || etabin2_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Trig, etabin2_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Trig, etabin2_Trig, Pt2, eta2);
        return -999;
    }
    Double_t Eff_muon2_data = Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso];
    Double_t Eff_muon2_MC = Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso];


    // -- Trigger part -- //
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_data = Eff_HLT_data_GtoH[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_GtoH[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

    Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

    //cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2)
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

        printf("[Data]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID],
               Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID],
               Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

        printf("[MC]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID],
               Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID],
               Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

        printf("[SF] Weight = %.3lf\n", weight);
    }
    return weight;

}// End of EfficiencySF_EventWeight_HLT_GtoH(mu1_4, mu2_4)


Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_GtoH_new(SelectedMuMu_t *MuMu)
{
    Double_t weight = -999;
    Double_t Pt1 = 0;
    Double_t eta1 = -9999;
    Double_t Pt2 = 0;
    Double_t eta2 = -9999;

    if (MuMu->Muon_pT->at(0) >= MuMu->Muon_pT->at(1))
    {
        Pt1 = MuMu->Muon_pT->at(0);
        eta1 = MuMu->Muon_eta->at(0);
        Pt2 = MuMu->Muon_pT->at(1);
        eta2 = MuMu->Muon_eta->at(1);
    }
    else
    {
        Pt1 = MuMu->Muon_pT->at(1);
        eta1 = MuMu->Muon_eta->at(1);
        Pt2 = MuMu->Muon_pT->at(0);
        eta2 = MuMu->Muon_eta->at(0);
    }

    // -- Muon1 LEADING -- //
    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig_new(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    if(ptbin1_ID == 9999 || etabin1_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_ID, etabin1_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_ID, etabin1_ID, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Iso == 9999 || etabin1_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Iso, etabin1_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Iso, etabin1_Iso, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Trig == 9999 || etabin1_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Trig, etabin1_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Trig, etabin1_Trig, Pt1, eta1);
        return -999;
    }

    Double_t Eff_muon1_data = Eff_ID_data_GtoH_lead[etabin1_ID][ptbin1_ID] * Eff_Iso_data_GtoH_lead[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_ID_MC_GtoH_lead[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_GtoH_lead[etabin1_Iso][ptbin1_Iso];

    // -- Muon2 SUBLEADING -- //
    Int_t ptbin2_ID = Find_muon_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_muon_EtaBin_ID(eta2);

    Int_t ptbin2_Iso = Find_muon_PtBin_Iso(Pt2);
    Int_t etabin2_Iso = Find_muon_EtaBin_Iso(eta2);

    Int_t ptbin2_Trig = Find_muon_PtBin_Trig_new(Pt2);
    Int_t etabin2_Trig = Find_muon_EtaBin_Trig(eta2);

    if(ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_ID, etabin2_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_ID, etabin2_ID, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Iso == 9999 || etabin2_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Iso, etabin2_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Iso, etabin2_Iso, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Trig == 9999 || etabin2_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Trig, etabin2_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Trig, etabin2_Trig, Pt2, eta2);
        return -999;
    }

    Double_t Eff_muon2_data = Eff_ID_data_GtoH_sublead[etabin2_ID][ptbin2_ID] * Eff_Iso_data_GtoH_sublead[etabin2_Iso][ptbin2_Iso];
    Double_t Eff_muon2_MC = Eff_ID_MC_GtoH_sublead[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_GtoH_sublead[etabin2_Iso][ptbin2_Iso];

    // Trigger SF
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH_lead[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_data = Eff_HLT_data_GtoH_sublead[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH_lead[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_GtoH_sublead[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

    Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

//    weight = (Eff_ID_data_GtoH_lead[etabin1_ID][ptbin1_ID] * Eff_ID_data_GtoH_sublead[etabin2_ID][ptbin2_ID]) /
//             (Eff_ID_MC_GtoH_lead[etabin1_ID][ptbin1_ID] * Eff_ID_MC_GtoH_sublead[etabin2_ID][ptbin2_ID]); // ONLY ID
//    weight = (Eff_muon1_data * Eff_muon2_data) / (Eff_muon1_MC * Eff_muon2_MC); // ID + ISO
    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2)
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

        printf("[Data]\n");
        printf("\t[Muon1] (ID, Iso): %.3lf, %.3lf)\n", Eff_ID_data_GtoH_lead[etabin1_ID][ptbin1_ID],
               Eff_Iso_data_GtoH_lead[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (ID, Iso): (%.3lf, %.3lf)\n", Eff_ID_data_GtoH_sublead[etabin2_ID][ptbin2_ID],
               Eff_Iso_data_GtoH_sublead[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

        printf("[MC]\n");
        printf("\t[Muon1] (ID, Iso): (%.3lf, %.3lf)\n", Eff_ID_MC_GtoH_lead[etabin1_ID][ptbin1_ID],
               Eff_Iso_MC_GtoH_lead[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (ID, Iso): (%.3lf, %.3lf)\n", Eff_ID_MC_GtoH_sublead[etabin2_ID][ptbin2_ID],
               Eff_Iso_MC_GtoH_sublead[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

        printf("[SF] Weight = %.3lf\n", weight);
    }
    return weight;

}// End of EfficiencySF_EventWeight_HLT_GtoH_new(mu1, mu2)


Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_GtoH(SelectedMuMu_t *MuMu)
{
    Double_t weight = -999;

    // -- Muon1 -- //
//    Double_t Pt1 = MuMu->Muon_TuneP_pT->at(0);
//    Double_t eta1 = MuMu->Muon_TuneP_eta->at(0);
    Double_t Pt1 = MuMu->Muon_pT->at(0);
    Double_t eta1 = MuMu->Muon_eta->at(0);

    Int_t ptbin1_Reco = Find_muon_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_muon_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    if(ptbin1_Reco == 9999 || etabin1_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Reco, etabin1_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Reco, etabin1_Reco, Pt1, eta1);
        cout << "pT1: " << Pt1 << "    eta1: " << eta1 << endl;
        return -999;
    }
    if(ptbin1_ID == 9999 || etabin1_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_ID, etabin1_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_ID, etabin1_ID, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Iso == 9999 || etabin1_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Iso, etabin1_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Iso, etabin1_Iso, Pt1, eta1);
        return -999;
    }
    if(ptbin1_Trig == 9999 || etabin1_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin1_Trig, etabin1_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin1_Trig, etabin1_Trig, Pt1, eta1);
        return -999;
    }

    Double_t Eff_muon1_data = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso];

    // -- Muon2 -- //
//    Double_t Pt2 = MuMu->Muon_TuneP_pT->at(1);
//    Double_t eta2 = MuMu->Muon_TuneP_eta->at(1);
    Double_t Pt2 = MuMu->Muon_pT->at(1);
    Double_t eta2 = MuMu->Muon_eta->at(1);

    Int_t ptbin2_Reco = Find_muon_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_muon_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_muon_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_muon_EtaBin_ID(eta2);

    Int_t ptbin2_Iso = Find_muon_PtBin_Iso(Pt2);
    Int_t etabin2_Iso = Find_muon_EtaBin_Iso(eta2);

    Int_t ptbin2_Trig = Find_muon_PtBin_Trig(Pt2);
    Int_t etabin2_Trig = Find_muon_EtaBin_Trig(eta2);

    if(ptbin2_Reco == 9999 || etabin2_Reco == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Reco, etabin2_Reco) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Reco, etabin2_Reco, Pt2, eta2);
        cout << "pT2: " << Pt2 << "    eta2: " << eta2 << endl;
        return -999;
    }
    if(ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_ID, etabin2_ID) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_ID, etabin2_ID, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Iso == 9999 || etabin2_Iso == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Iso, etabin2_Iso) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Iso, etabin2_Iso, Pt2, eta2);
        return -999;
    }
    if(ptbin2_Trig == 9999 || etabin2_Trig == 9999)
    {
        printf("ERROR! Wrong assigned bin number ... (ptbin2_Trig, etabin2_Trig) = (%d, %d), (pt, eta) = (%f, %f)\n", ptbin2_Trig, etabin2_Trig, Pt2, eta2);
        return -999;
    }
    Double_t Eff_muon2_data = Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso];
    Double_t Eff_muon2_MC = Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso];


    // -- Trigger part -- //
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_data = Eff_HLT_data_GtoH[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_Trig][ptbin1_Trig];
    Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_GtoH[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

    Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

    //cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2)
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

        printf("[Data]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID],
               Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID],
               Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

        printf("[MC]\n");
        printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID],
               Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso]);
        printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID],
               Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso]);
        printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

        printf("[SF] Weight = %.3lf\n", weight);
    }
    return weight;

}// End of EfficiencySF_EventWeight_HLT_GtoH(SelectedMuMu_t)


Double_t DYAnalyzer::EfficiencySF_EventWeight_electron(Electron ele1, Electron ele2)
{
    Double_t weight = -999;

    // -- Electron1 -- //
    Double_t Pt1 = ele1.Pt;
    //Double_t eta1 = ele1.eta;
    Double_t eta1 = ele1.etaSC;

    Int_t ptbin1_Reco = Find_electron_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_electron_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_electron_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_electron_EtaBin_ID(eta1);

    Int_t ptbin1_Trig = Find_electron_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_electron_EtaBin_Trig(eta1);

    Double_t Eff_ele1_data = Eff_Reco_data[etabin1_Reco][ptbin1_Reco] * Eff_ID_data[etabin1_ID][ptbin1_ID];
    Double_t Eff_ele1_MC = Eff_Reco_MC[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC[etabin1_ID][ptbin1_ID];

    // -- Electron2 -- //
    Double_t Pt2 = ele2.Pt;
    //Double_t eta2 = ele2.eta;
    Double_t eta2 = ele2.etaSC;

    Int_t ptbin2_Reco = Find_electron_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_electron_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_electron_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_electron_EtaBin_ID(eta2);

    Int_t ptbin2_Trig = Find_electron_PtBin_Trig(Pt2);
    Int_t etabin2_Trig = Find_electron_EtaBin_Trig(eta2);

    Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
    Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC[etabin2_ID][ptbin2_ID];

    // TRIGGERS
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Eff_EventTrig_data = Eff_HLT_Leg2_data[etabin1_Trig][ptbin1_Trig] * Eff_HLT_Leg2_data[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_MC = Eff_HLT_Leg2_MC[etabin1_Trig][ptbin1_Trig] * Eff_HLT_Leg2_MC[etabin2_Trig][ptbin2_Trig];

    //cout << Eff_EventTrig_data << "\t" << Eff_EventTrig_MC << endl;

    Double_t Eff_data_all = Eff_ele1_data * Eff_ele2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_ele1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2)
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);
        printf("[SF] Weight = %.3lf\n", weight);
    }

    /*if((Pt1 < 25 && eta1 > 2.3) || (Pt2 < 25 && eta2 > 2.3))
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);
        cout << Eff_EventTrig_data << "\t" << Eff_EventTrig_MC << endl;
        cout << weight << endl;
    }*/

    return weight;

}// End of EfficiencySF_EventWeight_electron(ele1, ele2)


Double_t DYAnalyzer::EfficiencySF_EventWeight_electron(SelectedEE_t *EE)
{
    Double_t weight = -999;

    // -- Electron1 -- //
    Double_t Pt1 = EE->Electron_pT->at(0);
//    Double_t eta1 = EE->Electron_eta->at(0);
    Double_t eta1 = EE->Electron_etaSC->at(0);

    Int_t ptbin1_Reco = Find_electron_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_electron_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_electron_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_electron_EtaBin_ID(eta1);

    Int_t ptbin1_Trig = Find_electron_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_electron_EtaBin_Trig(eta1);

    Double_t Eff_ele1_data = Eff_Reco_data[etabin1_Reco][ptbin1_Reco] * Eff_ID_data[etabin1_ID][ptbin1_ID];
    Double_t Eff_ele1_MC = Eff_Reco_MC[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC[etabin1_ID][ptbin1_ID];

    // -- Electron2 -- //
    Double_t Pt2 = EE->Electron_pT->at(1);
//        Double_t eta2 = EE->Electron_eta->at(1);
    Double_t eta2 = EE->Electron_etaSC->at(1);

    Int_t ptbin2_Reco = Find_electron_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_electron_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_electron_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_electron_EtaBin_ID(eta2);

    Int_t ptbin2_Trig = Find_electron_PtBin_Trig(Pt2);
    Int_t etabin2_Trig = Find_electron_EtaBin_Trig(eta2);

    Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
    Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC[etabin2_ID][ptbin2_ID];

    // TRIGGERS
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Eff_EventTrig_data = Eff_HLT_Leg2_data[etabin1_Trig][ptbin1_Trig] * Eff_HLT_Leg2_data[etabin2_Trig][ptbin2_Trig];
    Eff_EventTrig_MC = Eff_HLT_Leg2_MC[etabin1_Trig][ptbin1_Trig] * Eff_HLT_Leg2_MC[etabin2_Trig][ptbin2_Trig];

    //cout << Eff_EventTrig_data << "\t" << Eff_EventTrig_MC << endl;

    Double_t Eff_data_all = Eff_ele1_data * Eff_ele2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_ele1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2)
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);
        printf("[SF] Weight = %.3lf\n", weight);
    }

    /*if((Pt1 < 25 && eta1 > 2.3) || (Pt2 < 25 && eta2 > 2.3))
    {
        printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);
        cout << Eff_EventTrig_data << "\t" << Eff_EventTrig_MC << endl;
        cout << weight << endl;
    }*/

    return weight;

}// End of EfficiencySF_EventWeight_electron(SelectedEE_t)


Double_t DYAnalyzer::EfficiencySF_EventWeight_emu_BtoF(Muon mu, Electron ele)
{
    Double_t weight = -999;

    // -- Muon1 -- //
    Double_t Pt1 = mu.Pt;
    Double_t eta1 = mu.eta;

    Int_t ptbin1_Reco = Find_muon_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_muon_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    // -- Electron2 -- //
    Double_t Pt2 = ele.Pt;
//    Double_t eta2 = ele.eta;
    Double_t eta2 = ele.etaSC;

    Int_t ptbin2_Reco = Find_electron_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_electron_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_electron_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_electron_EtaBin_ID(eta2);

    //Check about bin settings
    if(ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 ||
        ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ...");
        return -999;
    }

    //Muon1
    Double_t Eff_muon1_data = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso];

    //Electron2
    Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
    Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC[etabin2_ID][ptbin2_ID];

    //TRIGGER : SingleMuon trigger only
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_Trig][ptbin1_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_Trig][ptbin1_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC;

    // -- emu event SF -- //
    Double_t Eff_data_all = Eff_muon1_data * Eff_ele2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2) printf("[SF] Weight = %.3lf\n", weight);
    return weight;

}// End of EfficiencySF_EventWeight_emu_BtoF(mu, ele)


Double_t DYAnalyzer::EfficiencySF_EventWeight_emu_BtoF(SelectedEMu_t *EMu)
{
    Double_t weight = -999;

    // -- Muon1 -- //
    Double_t Pt1 = EMu->Muon_pT;
    Double_t eta1 = EMu->Muon_eta;

    Int_t ptbin1_Reco = Find_muon_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_muon_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    // -- Electron2 -- //
    Double_t Pt2 = EMu->Electron_pT;
//    Double_t eta2 = EMu->Electron_eta;
    Double_t eta2 = EMu->Electron_etaSC;

    Int_t ptbin2_Reco = Find_electron_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_electron_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_electron_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_electron_EtaBin_ID(eta2);

    //Check about bin settings
    if(ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 ||
        ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ...");
        return -999;
    }

    //Muon1
    Double_t Eff_muon1_data = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso];

    //Electron2
    Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
    Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC[etabin2_ID][ptbin2_ID];

    //TRIGGER : SingleMuon trigger only
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_Trig][ptbin1_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_Trig][ptbin1_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC;

    // -- emu event SF -- //
    Double_t Eff_data_all = Eff_muon1_data * Eff_ele2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2) printf("[SF] Weight = %.3lf\n", weight);
    return weight;

}// End of EfficiencySF_EventWeight_emu_BtoF(SelectedEMu_t)


Double_t DYAnalyzer::EfficiencySF_EventWeight_emu_GtoH(Muon mu, Electron ele)
{
    Double_t weight = -999;

    // -- Muon1 -- //
    Double_t Pt1 = mu.Pt;
    Double_t eta1 = mu.eta;

    Int_t ptbin1_Reco = Find_muon_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_muon_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    // -- Electron2 -- //
    Double_t Pt2 = ele.Pt;
//    Double_t eta2 = ele.eta;
    Double_t eta2 = ele.etaSC;

    Int_t ptbin2_Reco = Find_electron_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_electron_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_electron_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_electron_EtaBin_ID(eta2);

    //Check about bin settings
    if(ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 ||
        ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ...");
        return -999;
    }

    //Muon1
    Double_t Eff_muon1_data = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso];

    //Electron2
    Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
    Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC[etabin2_ID][ptbin2_ID];

    //TRIGGER : SingleMuon trigger only
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_Trig][ptbin1_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_Trig][ptbin1_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC;

    // -- emu event SF -- //
    Double_t Eff_data_all = Eff_muon1_data * Eff_ele2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2) printf("[SF] Weight = %.3lf\n", weight);
    return weight;

}// End of EfficiencySF_EventWeight_emu_GtoH(mu, ele)


Double_t DYAnalyzer::EfficiencySF_EventWeight_emu_GtoH(SelectedEMu_t *EMu)
{
    Double_t weight = -999;

    // -- Muon1 -- //
    Double_t Pt1 = EMu->Muon_pT;
    Double_t eta1 = EMu->Muon_eta;

    Int_t ptbin1_Reco = Find_muon_PtBin_Reco(Pt1);
    Int_t etabin1_Reco = Find_muon_EtaBin_Reco(eta1);

    Int_t ptbin1_ID = Find_muon_PtBin_ID(Pt1);
    Int_t etabin1_ID = Find_muon_EtaBin_ID(eta1);

    Int_t ptbin1_Iso = Find_muon_PtBin_Iso(Pt1);
    Int_t etabin1_Iso = Find_muon_EtaBin_Iso(eta1);

    Int_t ptbin1_Trig = Find_muon_PtBin_Trig(Pt1);
    Int_t etabin1_Trig = Find_muon_EtaBin_Trig(eta1);

    // -- Electron2 -- //
    Double_t Pt2 = EMu->Electron_pT;
//    Double_t eta2 = EMu->Electron_eta;
    Double_t eta2 = EMu->Electron_etaSC;

    Int_t ptbin2_Reco = Find_electron_PtBin_Reco(Pt2);
    Int_t etabin2_Reco = Find_electron_EtaBin_Reco(eta2);

    Int_t ptbin2_ID = Find_electron_PtBin_ID(Pt2);
    Int_t etabin2_ID = Find_electron_EtaBin_ID(eta2);

    //Check about bin settings
    if(ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 ||
        ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
    {
        printf("ERROR! Wrong assigned bin number ...");
        return -999;
    }

    //Muon1
    Double_t Eff_muon1_data = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso];
    Double_t Eff_muon1_MC = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso];

    //Electron2
    Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
    Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC[etabin2_ID][ptbin2_ID];

    //TRIGGER : SingleMuon trigger only
    Double_t Eff_EventTrig_data = 0;
    Double_t Eff_EventTrig_MC = 0;

    Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_Trig][ptbin1_Trig];
    Eff_EventTrig_data = Eff_Trig_muon1_data;

    Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_Trig][ptbin1_Trig];
    Eff_EventTrig_MC = Eff_Trig_muon1_MC;

    // -- emu event SF -- //
    Double_t Eff_data_all = Eff_muon1_data * Eff_ele2_data * Eff_EventTrig_data;
    Double_t Eff_MC_all = Eff_muon1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

    weight = Eff_data_all / Eff_MC_all;

    if(weight > 2) printf("[SF] Weight = %.3lf\n", weight);
    return weight;

} // End of EfficiencySF_EventWeight_emu_GtoH(SelectedEMu_t)


void DYAnalyzer::SetupPVzWeights(Bool_t isMC, TString whichX, TString fileName)
{
    if(isMC == kFALSE)
        {
            for(Int_t i=0; i<80; i++)
                PVzWeight[i] = 1;
            return;
        }
    // -- MC only -- //
    TFile *f = new TFile(fileName, "READ");
    TString WhichX = whichX;
    WhichX.ToLower();
    if (WhichX!="ee" && WhichX!="mumu" && WhichX!="emu" && WhichX!="combined")
    {
        cout << "ERROR! Wrong name for the final state!" << endl;
//        return;
    }
    TH1D *h_PVzWeights = (TH1D*)f->Get("h_PVzWeights_"+WhichX);
    if(h_PVzWeights == NULL)
    {
        cout << "ERROR! No Weight histogram!"<< endl;
        return;
    }
    for(Int_t i=0; i<80; i++)
    {
        Int_t i_bin = i+1;
        PVzWeight[i] = h_PVzWeights->GetBinContent(i_bin);
    }
    return;
}


Double_t DYAnalyzer::PVzWeightValue(Double_t PVz)
{
    if(PVz < -20 || PVz >= 20)
        return 0;

    Double_t binEdges[81];
    for(int i=0; i<=80; i++)
        binEdges[i] = -20 + 0.5 * i;

    Int_t PVz_bin = 9999;

    for(Int_t i=0; i<80; i++)
    {
        if(binEdges[i] <= PVz && PVz < binEdges[i+1])
        {
            PVz_bin = i;
            break;
        }
    }

    if(PVz_bin == 9999)
    {
        cout << "ERROR: Could not find the PVz bin! (PVz value: " << PVz << ")" << endl;
        return 0;
    }

    return PVzWeight[PVz_bin];
}


Bool_t DYAnalyzer::EventSelection(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].passTightID && MuonCollection[j].RelPFIso_dBeta < 0.15)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
                if(mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*"))
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

        if(isExistHLTMatchedMuon == kTRUE)
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
                if(nQMuons == 2)
		{
			Muon recolep1 = QMuonCollection[0];
			Muon recolep2 = QMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
                        Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

			Bool_t isOS = kFALSE;
                        if(recolep1.charge != recolep2.charge) isOS = kTRUE;

                        if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE)
			{
				isPassEventSelection = kTRUE;
                                SelectedMuonCollection->push_back(recolep1);
                                SelectedMuonCollection->push_back(recolep2);
			}
		}
                else if(nQMuons > 2)
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
                                if(Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*"))
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

                                                if(j_mu != i_mu) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

                                                        if(isPassAcc == kTRUE) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
                                                                if(VtxNormChi2_temp < VtxNormChi2_BestPair)
								{
									VtxNormChi2_BestPair = VtxNormChi2_temp;
									mu1_BestPair = Mu;
									mu2_BestPair = Mu_jth;
								}
							}
						}
					} // -- end of the loop for j_mu (finding for second muon)
				}
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

                        if(VtxNormChi2_BestPair < 999) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
                                Double_t Angle = reco_v1.Angle(reco_v2.Vect());

				Bool_t isOS = kFALSE;
                                if(mu1_BestPair.charge != mu2_BestPair.charge) isOS = kTRUE;

                                if(reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE)
				{
					isPassEventSelection = kTRUE;
                                        SelectedMuonCollection->push_back(mu1_BestPair);
                                        SelectedMuonCollection->push_back(mu2_BestPair);
				}
			}

                } // -- End of else if(nQMuons > 2) -- //

        } // -- End of if(isExistHLTMatchedMuon == kTRUE) -- //

	return isPassEventSelection;
}


Bool_t DYAnalyzer::EventSelection(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
                                  vector< Muon >* SelectedMuonCollection, // -- output: 2 muons passing event selection conditions -- //
                                  vector< Int_t >* Index) // -- output: 2 indices of muons that managed to pass the MuMu selection -- //
{
    Bool_t isPassEventSelection = kFALSE;

    //Collect qualified muons among muons
    vector< Muon > QMuonCollection;
    vector< Int_t > QIndex;
    for(Int_t j=0; j<(int)MuonCollection.size(); j++)
    {
        if(MuonCollection[j].passTightID && MuonCollection[j].RelPFIso_dBeta < 0.15)
        {
            QMuonCollection.push_back(MuonCollection[j]);
            QIndex.push_back(j);
        }
    }

    // -- Check the existence of at least one muon matched with HLT-object -- //
    Bool_t isExistHLTMatchedMuon = kFALSE;
    for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
    {
            Muon mu = QMuonCollection[i_mu];
            if(mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*"))
            {
                    isExistHLTMatchedMuon = kTRUE;
                    break;
            }
    }
    if(isExistHLTMatchedMuon == kTRUE)
    {
            Int_t nQMuons = (Int_t)QMuonCollection.size();
            Int_t nQIndices = (Int_t)QIndex.size();
            if(nQMuons == 2 && nQIndices == 2)
            {
                    Muon recolep1 = QMuonCollection[0];
                    Muon recolep2 = QMuonCollection[1];

                    // -- Check the Acceptance -- //
                    Bool_t isPassAcc = kFALSE;
                    isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

                    Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                    Double_t VtxProb = -999;
                    Double_t VtxNormChi2 = 999;
                    DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

                    TLorentzVector inner_v1 = recolep1.Momentum_Inner;
                    TLorentzVector inner_v2 = recolep2.Momentum_Inner;

                    // -- 3D open angle -- //
                    Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

                    Bool_t isOS = kFALSE;
                    if(recolep1.charge != recolep2.charge) isOS = kTRUE;

                    if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE)
                    {
                            isPassEventSelection = kTRUE;
                            SelectedMuonCollection->push_back(recolep1);
                            SelectedMuonCollection->push_back(recolep2);
                            Index->push_back(QIndex[0]);
                            Index->push_back(QIndex[1]);
                    }
            }
            else if(nQMuons > 2)
            {
                    Double_t VtxProb_BestPair = -1;
                    Double_t VtxNormChi2_BestPair = 999;
                    Muon mu1_BestPair;
                    Muon mu2_BestPair;
                    Int_t index_1, index_2;

                    for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
                    {
                            Muon Mu = QMuonCollection[i_mu];

                            // -- at least 1 muon should be matched with HLT objects in best pair -- //
                            if(Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*"))
                            {
                                    // -- Mu in this loop: QMuon Matched with HLT object -- //

                                    // -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
                                    for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
                                    {
                                            Muon Mu_jth = QMuonCollection[j_mu];

                                            if(j_mu != i_mu) // -- do not calculate vertex variables(prob, chi2). with itself -- //
                                            {
                                                    // -- Check that this pair is within acceptance -- //
                                                    Bool_t isPassAcc = kFALSE;
                                                    isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

                                                    if(isPassAcc == kTRUE) // -- Find best pair ONLY for the pairs within acceptance -- //
                                                    {
                                                            Double_t VtxProb_temp = -999;
                                                            Double_t VtxNormChi2_temp = 999;
                                                            DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

                                                            // -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- //
                                                            if(VtxNormChi2_temp < VtxNormChi2_BestPair)
                                                            {
                                                                    VtxNormChi2_BestPair = VtxNormChi2_temp;
                                                                    mu1_BestPair = Mu;
                                                                    mu2_BestPair = Mu_jth;
                                                                    index_1 = QIndex[i_mu];
                                                                    index_2 = QIndex[j_mu];
                                                            }
                                                    }
                                            }
                                    } // -- end of the loop for j_mu (finding for second muon)
                            }
                    } // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

                    if(VtxNormChi2_BestPair < 999) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
                    {
                            TLorentzVector reco_v1 = mu1_BestPair.Momentum;
                            TLorentzVector reco_v2 = mu2_BestPair.Momentum;
                            Double_t reco_M = (reco_v1 + reco_v2).M();

                            // -- 3D open angle is calculated using inner track information -- //
                            // -- 3D open angle -- //
                            Double_t Angle = reco_v1.Angle(reco_v2.Vect());

                            Bool_t isOS = kFALSE;
                            if(mu1_BestPair.charge != mu2_BestPair.charge) isOS = kTRUE;

                            if(reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE)
                            {
                                    isPassEventSelection = kTRUE;
                                    SelectedMuonCollection->push_back(mu1_BestPair);
                                    SelectedMuonCollection->push_back(mu2_BestPair);
                                    Index->push_back(index_1);
                                    Index->push_back(index_2);
                            }
                    }

            } // -- End of else if(nQMuons > 2) -- //

    } // -- End of if(isExistHLTMatchedMuon == kTRUE) -- //

    return isPassEventSelection;
} // End of EventSelection(SelectedMuon, index)


// -- Test using the trigger without isolation condition: HLT_Mu50_v* -- //
Bool_t DYAnalyzer::EventSelection_Mu50(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
                if(mu.isTrigMatched(ntuple, HLT))
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

        if(isExistHLTMatchedMuon == kTRUE)
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
                if(nQMuons == 2)
		{
			Muon recolep1 = QMuonCollection[0];
			Muon recolep2 = QMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
                        Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

                        if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005)
			{
				isPassEventSelection = kTRUE;
                                SelectedMuonCollection->push_back(recolep1);
                                SelectedMuonCollection->push_back(recolep2);
			}
		}
                else if(nQMuons > 2)
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
                                if(Mu.isTrigMatched(ntuple, HLT))
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

                                                if(j_mu != i_mu) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

                                                        if(isPassAcc == kTRUE) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
                                                                if(VtxNormChi2_temp < VtxNormChi2_BestPair)
								{
									VtxNormChi2_BestPair = VtxNormChi2_temp;
									mu1_BestPair = Mu;
									mu2_BestPair = Mu_jth;
								}
							}
						}
					} // -- end of the loop for j_mu (finding for second muon)
				}
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

                        if(VtxNormChi2_BestPair < 999) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
                                Double_t Angle = reco_v1.Angle(reco_v2.Vect());

                                if(reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005)
				{
					isPassEventSelection = kTRUE;
                                        SelectedMuonCollection->push_back(mu1_BestPair);
                                        SelectedMuonCollection->push_back(mu2_BestPair);
				}
			}

                } // -- End of else if(nQMuons > 2) -- //

        } // -- End of if(isExistHLTMatchedMuon == kTRUE) -- //

	return isPassEventSelection;
}


Bool_t DYAnalyzer::EventSelection_minusDimuonVtxCut(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
                if(mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*"))
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

        if(isExistHLTMatchedMuon == kTRUE)
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
                if(nQMuons == 2)
		{
			Muon recolep1 = QMuonCollection[0];
			Muon recolep2 = QMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
                        Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

                        // if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005)
                        if(reco_M > 10 && isPassAcc == kTRUE && Angle < TMath::Pi() - 0.005)
			{
				isPassEventSelection = kTRUE;
                                SelectedMuonCollection->push_back(recolep1);
                                SelectedMuonCollection->push_back(recolep2);
			}
		}
                else if(nQMuons > 2)
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
                                if(Mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*"))
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

                                                if(j_mu != i_mu) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

                                                        if(isPassAcc == kTRUE) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
                                                                if(VtxNormChi2_temp < VtxNormChi2_BestPair)
								{
									VtxNormChi2_BestPair = VtxNormChi2_temp;
									mu1_BestPair = Mu;
									mu2_BestPair = Mu_jth;
								}
							}
						}
					} // -- end of the loop for j_mu (finding for second muon)
				}
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

                        if(VtxNormChi2_BestPair < 999) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
                                Double_t Angle = reco_v1.Angle(reco_v2.Vect());

                                // if(reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005)
                                if(reco_M > 10 && Angle < TMath::Pi() - 0.005)
				{
					isPassEventSelection = kTRUE;
                                        SelectedMuonCollection->push_back(mu1_BestPair);
                                        SelectedMuonCollection->push_back(mu2_BestPair);
				}
			}

                } // -- End of else if(nQMuons > 2) -- //

        } // -- End of if(isExistHLTMatchedMuon == kTRUE) -- //

	return isPassEventSelection;
}

// -- Event selection used for differential Z cross section measurement @ 13TeV -- // For HighPt muon id with 60-120GeV mass bin
Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

// -- Event selection used for differential Z cross section measurement @ 13TeV -- // For all mass bin, and HighPt id
Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
                                                        vector< Muon >* SelectedMuonCollection, // -- output: 2 muons passing event selection conditions -- //
                                                        vector< Int_t >* Index, // -- output: 2 indices of muons that managed to pass the MuMu selection -- //
                                                        Int_t &IndexDi)  // -- output: Index of the dimuon vector -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
        vector< Int_t > QIndex;
        Int_t DiIndex = -1;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
            //if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].RelPFIso_dBeta < 0.10)
            {
                QMuonCollection.push_back(MuonCollection[j]);
                QIndex.push_back(j);
            }
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();
        if(nQMuons == 2 && nQIndices == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);
                for(UInt_t i=0; i<ntuple->vtxTrkProb->size(); i++)
                {
                    if(VtxProb == ntuple->vtxTrkProb->at(i))
                    {
                        DiIndex = i;
                    }
                }

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
                        Index->push_back(QIndex[0]);
                        Index->push_back(QIndex[1]);
                        IndexDi = DiIndex;
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
                Int_t Index_leading = -1;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
                                Index_leading = QIndex[i_mu1];
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
                Int_t Index_sub = -1;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
                                Index_sub = QIndex[i_mu2];
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

                for(UInt_t i=0; i<ntuple->vtxTrkProb->size(); i++)
                {
                    if(VtxProb == ntuple->vtxTrkProb->at(i))
                    {
                        DiIndex = i;
                    }
                }

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
                        Index->push_back(Index_leading);
                        Index->push_back(Index_sub);
                        IndexDi = DiIndex;
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

// -- Event selection used for differential Z cross section measurement @ 13TeV -- // For all mass bin, and HighPt id
Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
                                                        vector< Muon >* SelectedMuonCollection, // -- output: 2 muons passing event selection conditions -- //
                                                        vector< Int_t >* Index) // -- output: 2 indices of muons that managed to pass the MuMu selection -- //
{
        Bool_t isPassEventSelection = kFALSE;

        //Collect qualified muons among muons
        vector< Muon > QMuonCollection;
        vector< Int_t > QIndex;
        for(Int_t j=0; j<(int)MuonCollection.size(); j++)
        {
            if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
            //if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].RelPFIso_dBeta < 0.10)
            {
                QMuonCollection.push_back(MuonCollection[j]);
                QIndex.push_back(j);
            }
        }

        Int_t nQMuons = (Int_t)QMuonCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();
        if(nQMuons == 2 && nQIndices == 2)
        {
                Muon recolep1 = QMuonCollection[0];
                Muon recolep2 = QMuonCollection[1];

                // -- Check the Accpetance -- //
                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

                // -- Opposite sign condition -- //
                Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
                        isOppositeSign = kTRUE;

                Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

                TLorentzVector inner_v1 = recolep1.Momentum_Inner;
                TLorentzVector inner_v2 = recolep2.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
                        Index->push_back(QIndex[0]);
                        Index->push_back(QIndex[1]);
                }
        }
        else if(nQMuons > 2)
        {
                // -- More then 2 Qualified Muon: Select the muons with highest pT -- //
                Double_t Pt_leading = 0;
                Muon LeadingMuon;
                Double_t i_leading = 0;
                Int_t Index_leading = -1;
                for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
                {
                        Muon Mu = QMuonCollection[i_mu1];

                        // printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
                        {
                                Pt_leading = Mu.Pt;
                                LeadingMuon	= Mu;
                                i_leading = i_mu1;
                                Index_leading = QIndex[i_mu1];
                        }
                }

                Double_t Pt_sub = 0;
                Muon SubMuon;
                Int_t Index_sub = -1;
                for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
                {
                        if(i_mu2 == i_leading) continue;

                        Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
                        {
                                Pt_sub = Mu.Pt;
                                SubMuon	= Mu;
                                Index_sub = QIndex[i_mu2];
                        }
                }

                // printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

                // -- Check the Accpetance -- //
                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

                // -- Opposite sign condition -- //
                Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
                        isOppositeSign = kTRUE;

                Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

                TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
                TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
                        Index->push_back(Index_leading);
                        Index->push_back(Index_sub);
                }

        } // -- End of else if(nQMuons > 2) -- //

        return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt(vector< Muon > MuonCollection, LongSelectedMuMu_t *ntuple, // -- input: All muons in a event & NtupleHandle -- //
                                                        vector< Muon >* SelectedMuonCollection, // -- output: 2 muons passing event selection conditions -- //
                                                        vector< Int_t >* Index, // -- output: 2 indices of muons that managed to pass the MuMu selection -- //
                                                        Int_t &IndexDi)  // -- output: Index of the dimuon vector -- //
{
        Bool_t isPassEventSelection = kFALSE;

        //Collect qualified muons among muons
        vector< Muon > QMuonCollection;
        vector< Int_t > QIndex;
        Int_t DiIndex = -1;
        for(Int_t j=0; j<(int)MuonCollection.size(); j++)
        {
            if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
            //if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].RelPFIso_dBeta < 0.10)
            {
                QMuonCollection.push_back(MuonCollection[j]);
                QIndex.push_back(j);
            }
        }

        Int_t nQMuons = (Int_t)QMuonCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();
        if(nQMuons == 2 && nQIndices == 2)
        {
                Muon recolep1 = QMuonCollection[0];
                Muon recolep2 = QMuonCollection[1];

                // -- Check the Accpetance -- //
                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

                // -- Opposite sign condition -- //
                Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
                        isOppositeSign = kTRUE;

                Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);
                for(UInt_t i=0; i<ntuple->vtxTrkProb->size(); i++)
                {
                    if(VtxProb == ntuple->vtxTrkProb->at(i))
                    {
                        DiIndex = i;
                    }
                }

                TLorentzVector inner_v1 = recolep1.Momentum_Inner;
                TLorentzVector inner_v2 = recolep2.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
                        Index->push_back(QIndex[0]);
                        Index->push_back(QIndex[1]);
                        IndexDi = DiIndex;
                }
        }
        else if(nQMuons > 2)
        {
                // -- More then 2 Qualified Muon: Select the muons with highest pT -- //
                Double_t Pt_leading = 0;
                Muon LeadingMuon;
                Double_t i_leading = 0;
                Int_t Index_leading = -1;
                for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
                {
                        Muon Mu = QMuonCollection[i_mu1];

                        // printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
                        {
                                Pt_leading = Mu.Pt;
                                LeadingMuon	= Mu;
                                i_leading = i_mu1;
                                Index_leading = QIndex[i_mu1];
                        }
                }

                Double_t Pt_sub = 0;
                Muon SubMuon;
                Int_t Index_sub = -1;
                for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
                {
                        if(i_mu2 == i_leading) continue;

                        Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
                        {
                                Pt_sub = Mu.Pt;
                                SubMuon	= Mu;
                                Index_sub = QIndex[i_mu2];
                        }
                }

                // printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

                // -- Check the Accpetance -- //
                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

                // -- Opposite sign condition -- //
                Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
                        isOppositeSign = kTRUE;

                Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

                for(UInt_t i=0; i<ntuple->vtxTrkProb->size(); i++)
                {
                    if(VtxProb == ntuple->vtxTrkProb->at(i))
                    {
                        DiIndex = i;
                    }
                }

                TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
                TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
                        Index->push_back(Index_leading);
                        Index->push_back(Index_sub);
                        IndexDi = DiIndex;
                }

        } // -- End of else if(nQMuons > 2) -- //

        return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt(vector< Muon > MuonCollection, LongSelectedMuMu_t *ntuple, // -- input: All muons in a event & NtupleHandle -- //
                                                        vector< Muon >* SelectedMuonCollection, // -- output: 2 muons passing event selection conditions -- //
                                                        vector< Int_t >* Index) // -- output: 2 indices of muons that managed to pass the MuMu selection -- //
{
        Bool_t isPassEventSelection = kFALSE;

        //Collect qualified muons among muons
        vector< Muon > QMuonCollection;
        vector< Int_t > QIndex;
        for(Int_t j=0; j<(int)MuonCollection.size(); j++)
        {
            if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
            //if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].RelPFIso_dBeta < 0.10)
            {
                QMuonCollection.push_back(MuonCollection[j]);
                QIndex.push_back(j);
            }
        }

        Int_t nQMuons = (Int_t)QMuonCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();
        if(nQMuons == 2 && nQIndices == 2)
        {
                Muon recolep1 = QMuonCollection[0];
                Muon recolep2 = QMuonCollection[1];

                // -- Check the Accpetance -- //
                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

                // -- Opposite sign condition -- //
                Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
                        isOppositeSign = kTRUE;

                Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

                TLorentzVector inner_v1 = recolep1.Momentum_Inner;
                TLorentzVector inner_v2 = recolep2.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
                        Index->push_back(QIndex[0]);
                        Index->push_back(QIndex[1]);
                }
        }
        else if(nQMuons > 2)
        {
                // -- More then 2 Qualified Muon: Select the muons with highest pT -- //
                Double_t Pt_leading = 0;
                Muon LeadingMuon;
                Double_t i_leading = 0;
                Int_t Index_leading = -1;
                for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
                {
                        Muon Mu = QMuonCollection[i_mu1];

                        // printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
                        {
                                Pt_leading = Mu.Pt;
                                LeadingMuon	= Mu;
                                i_leading = i_mu1;
                                Index_leading = QIndex[i_mu1];
                        }
                }

                Double_t Pt_sub = 0;
                Muon SubMuon;
                Int_t Index_sub = -1;
                for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
                {
                        if(i_mu2 == i_leading) continue;

                        Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
                        {
                                Pt_sub = Mu.Pt;
                                SubMuon	= Mu;
                                Index_sub = QIndex[i_mu2];
                        }
                }

                // printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

                // -- Check the Accpetance -- //
                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

                // -- Opposite sign condition -- //
                Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
                        isOppositeSign = kTRUE;

                Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

                TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
                TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
                        Index->push_back(Index_leading);
                        Index->push_back(Index_sub);
                }

        } // -- End of else if(nQMuons > 2) -- //

        return isPassEventSelection;
}


// -- Event selection used for N-1 cuts -- // For all mass bin, and HighPt id
Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt1(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_isGLB() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt2(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_muonHits() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt3(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_nMatches() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt4(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_dpT_over_pT() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt5(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_dxyVTX() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt6(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt7(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_pixelHits() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt8(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_trackerLayers() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt9(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
//		if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt10(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
//		if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
//		if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt11(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            //if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
            if(MuonCollection[j].isHighPtMuon())
                QMuonCollection.push_back(MuonCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        if(nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(recolep1.charge != recolep2.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1);
                        SelectedMuonCollection->push_back(recolep2);
		}
	}
        else if(nQMuons > 2)
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

                        if(Mu.Pt > Pt_leading)
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
                        if(i_mu2 == i_leading) continue;

			Muon Mu = QMuonCollection[i_mu2];

                        if(Mu.Pt > Pt_sub)
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
                if(LeadingMuon.charge != SubMuon.charge)
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle(SubMuon.Momentum.Vect());

//		if(isPassAcc == kTRUE && isOppositeSign == kTRUE)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(LeadingMuon);
                        SelectedMuonCollection->push_back(SubMuon);
		}

        } // -- End of else if(nQMuons > 2) -- //

	return isPassEventSelection;
} // End of EventSelection_Zdiff_13TeV_HighPt11


Bool_t DYAnalyzer::EventSelection_FR(vector<Muon> MuonCollection, NtupleHandle *ntuple, vector<Muon>* SelectedMuonCollection_nume,
                                     vector<Muon>* SelectedMuonCollection_deno)
{
    Bool_t isPassEventSelection = kFALSE;
    Bool_t skip = kTRUE;
    for(Int_t j=0; j<(int)MuonCollection.size(); j++)
    { // Asking for only one tight muon to surpass trigger threshold
        if(MuonCollection[j].Pt > LeadPtCut && fabs(MuonCollection[j].eta) < SubEtaCut && MuonCollection[j].isTightMuon())
            skip = kFALSE;
    }
    if (skip == kTRUE)
        return isPassEventSelection;
    for(Int_t j=0; j<(int)MuonCollection.size(); j++)
    { // All other muons still have to pass these criteria
        if (fabs(MuonCollection[j].eta) < SubEtaCut && (MuonCollection[j].isTightMuon()))
        {
            isPassEventSelection = kTRUE;
            SelectedMuonCollection_deno->push_back(MuonCollection[j]);
            if (MuonCollection[j].RelPFIso_dBeta <= 0.15)
                SelectedMuonCollection_nume->push_back(MuonCollection[j]);
        }
    }
    return isPassEventSelection;
} // End of EventSelection_FR()



Bool_t DYAnalyzer::EventSelection_FRdijetEst(vector<Muon> MuonCollection, NtupleHandle *ntuple, vector<Muon>* SelectedMuonCollection_fail)
{
    Bool_t isPassEventSelection = kFALSE;
    vector< Muon > TempMuonCollection;
    for (Int_t j=0; j<(int)MuonCollection.size(); j++)
    {
        if (fabs(MuonCollection[j].eta) < SubEtaCut && MuonCollection[j].isTightMuon() && MuonCollection[j].RelPFIso_dBeta >= 0.15)
            TempMuonCollection.push_back(MuonCollection[j]);
    }
    if (TempMuonCollection.size() == 2)
    {
        if ((TempMuonCollection[0].Pt > LeadPtCut && TempMuonCollection[1].Pt > SubPtCut) || (TempMuonCollection[0].Pt > SubPtCut && TempMuonCollection[1].Pt > LeadPtCut))
        {
            isPassEventSelection = kTRUE;
            for (Int_t i=0; i<2; i++)
                SelectedMuonCollection_fail->push_back(TempMuonCollection[i]);
        }
    }
    return isPassEventSelection;
} // End of EventSelection_FRdijetEst()



Bool_t DYAnalyzer::EventSelection_FakeMuons_Triggerless(vector<Muon> MuonCollection, NtupleHandle *ntuple, vector<Muon>* SelectedMuonCollection)
{
    Bool_t isPassEventSelection = kFALSE;

    //Collect qualified muons among muons
    vector< Muon > QMuonCollection;
    for(Int_t j=0; j<(int)MuonCollection.size(); j++)
    {
        if(MuonCollection[j].isTightMuon())//passTightID)
            QMuonCollection.push_back(MuonCollection[j]);
    }

    Int_t nQMuons = (Int_t)QMuonCollection.size();
    if (nQMuons >= 2)
    {
        if(nQMuons == 2)
        {
            Muon recolep1 = QMuonCollection[0];
            Muon recolep2 = QMuonCollection[1];

            // -- Check the Accpetance -- //
            Bool_t isPassAcc = kFALSE;
            isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

            Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

            Double_t VtxProb = -999;
            Double_t VtxNormChi2 = 999;
            DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

            TLorentzVector inner_v1 = recolep1.Momentum_Inner;
            TLorentzVector inner_v2 = recolep2.Momentum_Inner;

            // -- 3D open angle -- //
            Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

            Bool_t isOS = kFALSE;
            if(recolep1.charge != recolep2.charge) isOS = kTRUE;

            if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005) // can have same sign
            {
                isPassEventSelection = kTRUE;
                SelectedMuonCollection->push_back(recolep1);
                SelectedMuonCollection->push_back(recolep2);
            }
        }
        else // nQMuons > 2
        {
            Double_t VtxProb_BestPair = -1;
            Double_t VtxNormChi2_BestPair = 999;
            Muon mu1_BestPair;
            Muon mu2_BestPair;

            for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
            {
                Muon Mu = QMuonCollection[i_mu];
                for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
                {
                    Muon Mu_jth = QMuonCollection[j_mu];
                    if(j_mu != i_mu) // do not calculate vertex variables(prob, chi2). with itself
                    {
                        // Check that this pair is within acceptance
                        Bool_t isPassAcc = kFALSE;
                        isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

                        if(isPassAcc == kTRUE) // Find best pair from pairs within acceptance
                        {
                            Double_t VtxProb_temp = -999;
                            Double_t VtxNormChi2_temp = 999;
                            DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

                            // Find best pair with smallest Chi2/dnof(VTX)
                            if(VtxNormChi2_temp < VtxNormChi2_BestPair)
                            {
                                VtxNormChi2_BestPair = VtxNormChi2_temp;
                                mu1_BestPair = Mu;
                                mu2_BestPair = Mu_jth;
                            }
                        }
                    }
                } // End of for(j_mu)
            } // End of for(i_mu)

            if(VtxNormChi2_BestPair < 999) // At least one pair within acceptance
            {
                TLorentzVector reco_v1 = mu1_BestPair.Momentum;
                TLorentzVector reco_v2 = mu2_BestPair.Momentum;
                Double_t reco_M = (reco_v1 + reco_v2).M();

                // 3D open angle
                Double_t Angle = reco_v1.Angle(reco_v2.Vect());

                Bool_t isOS = kFALSE;
                if(mu1_BestPair.charge != mu2_BestPair.charge) isOS = kTRUE;

                if(reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005) // can have same sign
                {
                    isPassEventSelection = kTRUE;
                    SelectedMuonCollection->push_back(mu1_BestPair);
                    SelectedMuonCollection->push_back(mu2_BestPair);
                }
            }
        } // End of if(nQMuons > 2)
    } // End of if(nQMuons >= 2)
    if (isPassEventSelection == kTRUE && (SelectedMuonCollection->at(0).RelPFIso_dBeta >= 0.15 || SelectedMuonCollection->at(1).RelPFIso_dBeta >= 0.15)) // At least one failing
        return isPassEventSelection;
    return kFALSE;
} // End of EventSelection_FakeMuons_Triggerless()



Bool_t DYAnalyzer::EventSelection_FRsingleJetEst(vector<Muon> MuonCollection, NtupleHandle *ntuple, vector<Muon> *SelectedMuonCollection)
{
    Bool_t isPassEventSelection = kFALSE;
    SelectedMuonCollection->clear();
    Muon TempMuon_lead, TempMuon_sublead;
    Double_t pT_lead=-999, pT_sublead=-999;
    Int_t nTightMuons = 0;
    for (Int_t j=0; j<(int)MuonCollection.size(); j++)
    {
        if (fabs(MuonCollection[j].eta) < SubEtaCut && MuonCollection[j].isTightMuon())
        {
            nTightMuons++;
            if (MuonCollection[j].Pt > pT_lead)
            {
                TempMuon_lead = MuonCollection[j];
                pT_lead = MuonCollection[j].Pt;
            }
            else if (MuonCollection[j].Pt > pT_sublead)
            {
                TempMuon_sublead = MuonCollection[j];
                pT_sublead = MuonCollection[j].Pt;
            }
        }
    }
    if (nTightMuons == 2 && pT_lead > LeadPtCut && pT_sublead > SubPtCut)
    {
        if (TempMuon_lead.RelPFIso_dBeta < 0.15 && TempMuon_sublead.RelPFIso_dBeta >= 0.15)
        {
            isPassEventSelection = kTRUE;
            SelectedMuonCollection->push_back(TempMuon_lead);
            SelectedMuonCollection->push_back(TempMuon_sublead);
        }
        else if (TempMuon_sublead.RelPFIso_dBeta < 0.15 && TempMuon_lead.RelPFIso_dBeta >= 0.15)
        {
            isPassEventSelection = kTRUE;
            SelectedMuonCollection->push_back(TempMuon_sublead);
            SelectedMuonCollection->push_back(TempMuon_lead);
        }
    }
    return isPassEventSelection;
} // End of EventSelection_FRsingleJetEst()


Bool_t DYAnalyzer::EventSelection_FR(vector<Electron> ElectronCollection, NtupleHandle *ntuple, vector<Electron> *SelectedElectronCollection) // Electron selection
{return kTRUE;}


void DYAnalyzer::SetupFRvalues(TString filename, TString type) // type can be "ratio", "template", "mixed", "sigCtrl_template" or "dalmin"
{
    // -- Setting up -- //
    std::cout << "Setting up fake rate values from " << filename << endl;
    std::cout << "Fake rate obtained using " << type << " method." << endl;
    TFile *f = new TFile(filename, "READ");
    TH1D *h_FR_barrel, *h_FR_endcap;
    f->GetObject("h_FR"+type+"_barrel", h_FR_barrel);  // +"up"  +"down"
    f->GetObject("h_FR"+type+"_endcap", h_FR_endcap);  // +"up"  +"down"

    // -- Getting values from histograms -- //
    for (Int_t i_bin=1; i_bin<=nPtBinBarrel; i_bin++)
    {
        FR_barrel[i_bin-1] = h_FR_barrel->GetBinContent(i_bin);
        if (i_bin <= nPtBinEndcap)
            FR_endcap[i_bin-1] = h_FR_endcap->GetBinContent(i_bin);
    }

    // -- Checking if everything has been done correctly -- //
    Int_t nProblem = 0;
    for (Int_t i=0; i<nPtBinBarrel; i++)
    {
        if (FR_barrel[i] >= 1 || FR_barrel[i] <= 0) {nProblem++; std::cout << i << "  " << FR_barrel[i] << endl;}
        if (i < nPtBinEndcap) { if (FR_endcap[i] >= 1 || FR_endcap[i] <= 0) {nProblem++; std::cout << i << "  " << FR_endcap[i] << endl;} }
    }
    if (nProblem)
        std::cout << "**************************************************\n" <<
                     "Problems were detected: fake rate values exceeding 'regular' probability boundaries have been found.\n" <<
                     "Please check the code.\n**************************************************" << endl;
    else std::cout << "Fake rates have been set up successfully." << endl;
    return;
} // End of SetupFRvalues()



Double_t DYAnalyzer::FakeRate(Double_t p_T, Double_t eta)
{
    Int_t i_bin = 0;
    Int_t stop = 0;
    if (fabs(eta) < 1.2) // barrel
    {
        while (!stop)
        {
            if (p_T < ptbin_barrel[i_bin + 1] || i_bin >= nPtBinBarrel-1) // Points exceeding boundaries are assigned last available FR value
                stop = 1;
            else
                i_bin++;
        }
        return FR_barrel[i_bin];
    }
    else // endcap
    {
        while (!stop)
        {
            if (p_T < ptbin_endcap[i_bin + 1] || i_bin >= nPtBinEndcap-1) // Points exceeding boundaries are assigned last available FR value
                stop = 1;
            else
                i_bin++;
        }
        return FR_endcap[i_bin];
    }
} // End of FakeRate()



Bool_t DYAnalyzer::EventSelection_emu_method(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple,
                                             vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection)
{
    Bool_t isPassEventSelection = kFALSE;

    //Collect qualified muons among muons
    vector< Muon > QMuonCollection;
    for(Int_t j=0; j<(int)MuonCollection.size(); j++)
    {
        if(/*MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10 &&*/
            MuonCollection[j].isTightMuon() && MuonCollection[j].relPFiso < 0.15 &&
            MuonCollection[j].Pt > SubPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut) // pT>17 && |eta|<2.4
                QMuonCollection.push_back(MuonCollection[j]);
    }

    //Collect qualified electrons among electrons
    vector< Electron > QElectronCollection;
    for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
    {
        Electron elec = ElectronCollection[j];
        if(elec.passMediumID == kTRUE && elec.Pt > SubPtCut && fabs(elec.etaSC) < LeadEtaCut &&
            !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566)) // pT>17 && |eta|<2.4
                QElectronCollection.push_back(ElectronCollection[j]);
    }

    Int_t nQMuons = (Int_t)QMuonCollection.size();
    Int_t nQElectrons = (Int_t)QElectronCollection.size();

    Double_t VtxProb_BestPair = -1;
    Double_t VtxNormChi2_BestPair = 999;
    Muon mu_BestPair;
    Electron el_BestPair;

    // -- Select muon with highest pT -- //
    for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
    {
        Muon Mu = QMuonCollection[i_mu];

        // -- muon should be matched with HLT objects in emu best pair -- //
        if(Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*"))
        {
            // -- Start another loop for finding electron (for electron, we don't need to check about trigger) -- //
            for(Int_t j_el=0; j_el<nQElectrons; j_el++)
            {
                Electron El = QElectronCollection[j_el];

                Double_t VtxProb_temp = -999;
                Double_t VtxNormChi2_temp = 999;
                emuVertexProbNormChi2(ntuple, El.gsfpT, Mu.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

                // -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- //
                if(VtxNormChi2_temp < VtxNormChi2_BestPair)
                {
                    VtxNormChi2_BestPair = VtxNormChi2_temp;
                    mu_BestPair = Mu;
                    el_BestPair = El;
                }
            } // -- end of the loop for j_el (finding for electron)
        }
    } // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

    if(VtxNormChi2_BestPair < 999)
    {
        TLorentzVector reco_v1 = mu_BestPair.Momentum;
        TLorentzVector reco_v2 = el_BestPair.Momentum;
        Double_t reco_M = (reco_v1 + reco_v2).M();

        Bool_t isPassAcc = kFALSE;
        if(mu_BestPair.Pt > LeadPtCut || el_BestPair.Pt > LeadPtCut) isPassAcc = kTRUE;

        // -- 3D open angle -- //
        Double_t Angle = reco_v1.Angle(reco_v2.Vect());

        if(reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isPassAcc == kTRUE)
        {
            isPassEventSelection = kTRUE;
            SelectedMuonCollection->push_back(mu_BestPair);
            SelectedElectronCollection->push_back(el_BestPair);
        }
    }

    return isPassEventSelection;
}

// Derived by Marijus Ambrozas 2018.08.07 to return indices of particles that passed the selection
Bool_t DYAnalyzer::EventSelection_emu_method(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple, // Input: electron, muon vectors and NtupleHandle
                                                vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection, // Output: Selected electron and muon
                                                  Int_t &Sel_Index_Mu, Int_t &Sel_Index_Ele)      // Output: Indices of selected electron and muon
{
    Bool_t isPassEventSelection = kFALSE;
    Sel_Index_Mu = -1; Sel_Index_Ele = -1;

    //Collect qualified muons among muons
    vector< Muon > QMuonCollection;
    vector< Int_t > QIndexMu;
    for(Int_t j=0; j<(int)MuonCollection.size(); j++)
    {
        if(MuonCollection[j].isTightMuon() && MuonCollection[j].relPFiso < 0.15 &&
            MuonCollection[j].Pt > SubPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut) // pT>17 && |eta|<2.4
        {
            QMuonCollection.push_back(MuonCollection[j]);
            QIndexMu.push_back(j);
        }
    }

    //Collect qualified electrons among electrons
    vector< Electron > QElectronCollection;
    vector< Int_t > QIndexEle;
    for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
    {
        Electron elec = ElectronCollection[j];
        if(elec.passMediumID == kTRUE && elec.Pt > SubPtCut && fabs(elec.etaSC) < LeadEtaCut &&
            !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566)) // pT>17 && |eta|<2.4
        {
            QElectronCollection.push_back(ElectronCollection[j]);
            QIndexEle.push_back(j);
        }
    }

    Int_t nQMuons = (Int_t)QMuonCollection.size();
    Int_t nQElectrons = (Int_t)QElectronCollection.size();
    Int_t nQIndicesMu = (Int_t)QIndexMu.size();
    Int_t nQIndicesEle = (Int_t)QIndexEle.size();

    Double_t VtxProb_BestPair = -1;
    Double_t VtxNormChi2_BestPair = 999;
    Int_t Index_mu = -1;
    Int_t Index_ele = -1;
    Muon mu_BestPair;
    Electron el_BestPair;

    if (nQMuons == nQIndicesMu && nQElectrons == nQIndicesEle)
    {
        for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
        {
            Muon Mu = QMuonCollection[i_mu];

            // -- muon should be matched with HLT objects in emu best pair -- //
            if(Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*"))
            {
                // -- Start another loop for finding electron (for electron, we don't need to check about trigger) -- //
                for(Int_t j_el=0; j_el<nQElectrons; j_el++)
                {
                    Electron El = QElectronCollection[j_el];

                    Double_t VtxProb_temp = -999;
                    Double_t VtxNormChi2_temp = 999;
                    emuVertexProbNormChi2(ntuple, El.gsfpT, Mu.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

                    // -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- //
                    if(VtxNormChi2_temp < VtxNormChi2_BestPair)
                    {
                        VtxNormChi2_BestPair = VtxNormChi2_temp;
                        mu_BestPair = Mu;
                        el_BestPair = El;
                        Index_mu = QIndexMu[i_mu];
                        Index_ele = QIndexEle[j_el];
                    }
                } // -- end of the loop for j_el (finding for electron)
            }
        }// -- end of the loop for i_mu (finding for the first muon matched with HLT matching)
    }

    if(VtxNormChi2_BestPair < 999)
    {
        TLorentzVector reco_v1 = mu_BestPair.Momentum;
        TLorentzVector reco_v2 = el_BestPair.Momentum;
        Double_t reco_M = (reco_v1 + reco_v2).M();

        Bool_t isPassAcc = kFALSE;
        if(mu_BestPair.Pt > LeadPtCut || el_BestPair.Pt > LeadPtCut) isPassAcc = kTRUE;

        // -- 3D open angle -- //
        Double_t Angle = reco_v1.Angle(reco_v2.Vect());

        if(reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isPassAcc == kTRUE)
        {
            isPassEventSelection = kTRUE;
            SelectedMuonCollection->push_back(mu_BestPair);
            SelectedElectronCollection->push_back(el_BestPair);
            Sel_Index_Mu = Index_mu;
            Sel_Index_Ele = Index_ele;
        }
    }

    return isPassEventSelection;
}

// test for emu event selection not using vtx cut
// Derived by Marijus Ambrozas 2018.08.07 to use LongSelectedEMu_t instead of NtupleHandle
Bool_t DYAnalyzer::EventSelection_emu_method(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, LongSelectedEMu_t *ntuple, // Input: electron, muon vectors and LongSelectedEMu_t
                                                vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection, // Output: Selected electron and muon
                                                  Int_t &Sel_Index_Mu, Int_t &Sel_Index_Ele)      // Output: Indices of selected electron and muon
{
    Bool_t isPassEventSelection = kFALSE;
    Sel_Index_Mu = -1; Sel_Index_Ele = -1;

    //Collect qualified muons among muons
    vector< Muon > QMuonCollection;
    vector< Int_t > QIndexMu;
    for(Int_t j=0; j<(int)MuonCollection.size(); j++)
    {
        if(/*MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10 &&*/
            MuonCollection[j].isTightMuon() && MuonCollection[j].relPFiso < 0.15 &&
            MuonCollection[j].Pt > SubPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut) // pT>17 && |eta|<2.4
        {
            QMuonCollection.push_back(MuonCollection[j]);
            QIndexMu.push_back(j);
        }
    }

    //Collect qualified electrons among electrons
    vector< Electron > QElectronCollection;
    vector< Int_t > QIndexEle;
    for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
    {
        Electron elec = ElectronCollection[j];
        if(elec.passMediumID == kTRUE && elec.Pt > SubPtCut && fabs(elec.etaSC) < LeadEtaCut &&
            !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566)) // pT>17 && |eta|<2.4
        {
            QElectronCollection.push_back(ElectronCollection[j]);
            QIndexEle.push_back(j);
        }
    }

    Int_t nQMuons = (Int_t)QMuonCollection.size();
    Int_t nQElectrons = (Int_t)QElectronCollection.size();
    Int_t nQIndicesMu = (Int_t)QIndexMu.size();
    Int_t nQIndicesEle = (Int_t)QIndexEle.size();

    Double_t VtxProb_BestPair = -1;
    Double_t VtxNormChi2_BestPair = 999;
    Int_t Index_mu = -1;
    Int_t Index_ele = -1;
    Muon mu_BestPair;
    Electron el_BestPair;

    if (nQMuons == nQIndicesMu && nQElectrons == nQIndicesEle)
    {
        for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
        {
            Muon Mu = QMuonCollection[i_mu];

            // -- muon should be matched with HLT objects in emu best pair -- //
            if(Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*"))
            {
                // -- Start another loop for finding electron (for electron, we don't need to check about trigger) -- //
                for(Int_t j_el=0; j_el<nQElectrons; j_el++)
                {
                    Electron El = QElectronCollection[j_el];

                    Double_t VtxProb_temp = -999;
                    Double_t VtxNormChi2_temp = 999;
                    emuVertexProbNormChi2(ntuple, El.gsfpT, Mu.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

                    // -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- //
                    if(VtxNormChi2_temp < VtxNormChi2_BestPair)
                    {
                        VtxNormChi2_BestPair = VtxNormChi2_temp;
                        mu_BestPair = Mu;
                        el_BestPair = El;
                        Index_mu = QIndexMu[i_mu];
                        Index_ele = QIndexEle[j_el];
                    }
                } // -- end of the loop for j_el (finding for electron)
            }
        }
    }// -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

    if(VtxNormChi2_BestPair < 999)
    {
        TLorentzVector reco_v1 = mu_BestPair.Momentum;
        TLorentzVector reco_v2 = el_BestPair.Momentum;
        Double_t reco_M = (reco_v1 + reco_v2).M();

        Bool_t isPassAcc = kFALSE;
        if(mu_BestPair.Pt > LeadPtCut || el_BestPair.Pt > LeadPtCut) isPassAcc = kTRUE;

        // -- 3D open angle -- //
        Double_t Angle = reco_v1.Angle(reco_v2.Vect());

        if(reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isPassAcc == kTRUE)
        {
            isPassEventSelection = kTRUE;
            SelectedMuonCollection->push_back(mu_BestPair);
            SelectedElectronCollection->push_back(el_BestPair);
            Sel_Index_Mu = Index_mu;
            Sel_Index_Ele = Index_ele;
        }
    }

    return isPassEventSelection;
}

/*
Bool_t DYAnalyzer::EventSelection_emu_method(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple,
						vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection)
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
                if(MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10
                        && MuonCollection[j].Pt > LeadPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut)
                        QMuonCollection.push_back(MuonCollection[j]);
	}

	//Collect qualified electrons among electrons
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                if(elec.passMediumID == kTRUE
                        //&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < 2.5 && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	Double_t VtxProb_BestPair = -1;
	Double_t VtxNormChi2_BestPair = 999;
	Muon mu_BestPair;
	Electron el_BestPair;

	for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
	{
		Muon Mu = QMuonCollection[i_mu];

		// -- muon should be matched with HLT objects in emu best pair -- //
                if(Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*"))
		{
			// -- Start another loop for finding electron (for electron, we don't need to check about trigger) -- //
			for(Int_t j_el=0; j_el<nQElectrons; j_el++)
			{
				Electron El = QElectronCollection[j_el];

				Double_t VtxProb_temp = -999;
				Double_t VtxNormChi2_temp = 999;
				emuVertexProbNormChi2(ntuple, Mu.Inner_pT, El.gsfpT, &VtxProb_temp, &VtxNormChi2_temp);

				// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
                                if(VtxNormChi2_temp < VtxNormChi2_BestPair)
				{
					VtxNormChi2_BestPair = VtxNormChi2_temp;
					mu_BestPair = Mu;
					el_BestPair = El;
				}
			} // -- end of the loop for j_el (finding for electron)
		}
	} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

        if(VtxNormChi2_BestPair < 999) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
	{
		TLorentzVector reco_v1 = mu_BestPair.Momentum;
		TLorentzVector reco_v2 = el_BestPair.Momentum;
		Double_t reco_M = (reco_v1 + reco_v2).M();

		// -- 3D open angle -- //
                Double_t Angle = reco_v1.Angle(reco_v2.Vect());

                if(reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(mu_BestPair);
                        SelectedElectronCollection->push_back(el_BestPair);
		}
	}

	return isPassEventSelection;
}
*/
Bool_t DYAnalyzer::isPassAccCondition_Muon(Muon Mu1, Muon Mu2)
{
	Bool_t isPassAcc = kFALSE;
	Muon leadMu, subMu;
	CompareMuon(&Mu1, &Mu2, &leadMu, &subMu);
        if(leadMu.Pt > LeadPtCut && fabs(leadMu.eta) < LeadEtaCut &&
                subMu.Pt  > SubPtCut  && fabs(subMu.eta)  < SubEtaCut)
		isPassAcc = kTRUE;

	return isPassAcc;
}

Bool_t DYAnalyzer::isPassAccCondition_GenLepton(GenLepton genlep1, GenLepton genlep2)
{
	Bool_t isPassAcc = kFALSE;

	GenLepton leadGenLep, subGenLep;
	CompareGenLepton(&genlep1, &genlep2, &leadGenLep, &subGenLep);
	
        if(leadGenLep.Pt > LeadPtCut && fabs(leadGenLep.eta) < LeadEtaCut &&
                subGenLep.Pt  > SubPtCut  && fabs(subGenLep.eta) < SubEtaCut)
		isPassAcc = 1;

	return isPassAcc;
}

void DYAnalyzer::CompareMuon(Muon *Mu1, Muon *Mu2, Muon *leadMu, Muon *subMu)
{
    if(Mu1->Pt > Mu2->Pt)
    {
        *leadMu = *Mu1;
        *subMu = *Mu2;
    }
    else
    {
        *leadMu = *Mu2;
        *subMu = *Mu1;
    }
}

void DYAnalyzer::CompareGenLepton(GenLepton *genlep1, GenLepton *genlep2, GenLepton *leadgenlep, GenLepton *subgenlep)
{
        if(genlep1->Pt > genlep2->Pt)
	{
		*leadgenlep = *genlep1;
		*subgenlep = *genlep2;
	}
	else
	{
		*leadgenlep = *genlep2;
		*subgenlep = *genlep1;
	}
}

void DYAnalyzer::DimuonVertexProbNormChi2(NtupleHandle *ntuple, Double_t Pt1, Double_t Pt2, Double_t *VtxProb, Double_t *VtxNormChi2)
{
	vector<double> *PtCollection1 = ntuple->vtxTrkCkt1Pt;
	vector<double> *PtCollection2 = ntuple->vtxTrkCkt2Pt;
	vector<double> *VtxProbCollection = ntuple->vtxTrkProb;

	Int_t NPt1 = (Int_t)PtCollection1->size();
	Int_t NPt2 = (Int_t)PtCollection2->size();
	Int_t NProb = (Int_t)VtxProbCollection->size();

        if(NPt1 != NPt2 || NPt2 != NProb || NPt1 != NProb)
                std::cout << "NPt1: " << NPt1 << " NPt2: " << NPt2 << " Nprob: " << NProb << endl;

        // std::cout << "Pt1: " << Pt1 << " Pt2: " << Pt2 << endl;

	for(Int_t i=0; i<NProb; i++)
	{
                // std::cout << "\tPtCollection1->at("<< i << "): " << PtCollection1->at(i) << " PtCollection2->at("<< i << "): " << PtCollection2->at(i) << endl;
                if((PtCollection1->at(i) == Pt1 && PtCollection2->at(i) == Pt2)  || (PtCollection1->at(i) == Pt2 && PtCollection2->at(i) == Pt1))
		{
			*VtxProb = VtxProbCollection->at(i);
			*VtxNormChi2 = ntuple->vtxTrkChi2->at(i) / ntuple->vtxTrkNdof->at(i);
			break;
		}
	}

	return;
}

void DYAnalyzer::DimuonVertexProbNormChi2(LongSelectedMuMu_t *ntuple, Double_t Pt1, Double_t Pt2, Double_t *VtxProb, Double_t *VtxNormChi2)
{
        vector<double> *PtCollection1 = ntuple->vtxTrkCkt1Pt;
        vector<double> *PtCollection2 = ntuple->vtxTrkCkt2Pt;
        vector<double> *VtxProbCollection = ntuple->vtxTrkProb;

        Int_t NPt1 = (Int_t)PtCollection1->size();
        Int_t NPt2 = (Int_t)PtCollection2->size();
        Int_t NProb = (Int_t)VtxProbCollection->size();

        if(NPt1 != NPt2 || NPt2 != NProb || NPt1 != NProb)
                std::cout << "NPt1: " << NPt1 << " NPt2: " << NPt2 << " Nprob: " << NProb << endl;

        // std::cout << "Pt1: " << Pt1 << " Pt2: " << Pt2 << endl;

        for(Int_t i=0; i<NProb; i++)
        {
                // std::cout << "\tPtCollection1->at("<< i << "): " << PtCollection1->at(i) << " PtCollection2->at("<< i << "): " << PtCollection2->at(i) << endl;
                if((PtCollection1->at(i) == Pt1 && PtCollection2->at(i) == Pt2)  || (PtCollection1->at(i) == Pt2 && PtCollection2->at(i) == Pt1))
                {
                        *VtxProb = VtxProbCollection->at(i);
                        *VtxNormChi2 = ntuple->vtxTrkChi2->at(i) / ntuple->vtxTrkNdof->at(i);
                        break;
                }
        }

        return;
}

// -- Event selecton for the electron channel (test) -- //
Bool_t DYAnalyzer::EventSelection_Electron(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
                if(elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 && elec.Pt > 15)
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if(nQElectrons == 2)
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
		}
	}
	return isPassEventSelection;

}

// -- Event selecton for the electron channel (2016.02.11) -- // modified at 17 May 2017 by Dalmin Pai
Bool_t DYAnalyzer::EventSelection_ElectronChannel(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                if(elec.passMediumID == kTRUE // modified by Dalmin Pai
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        if(nQElectrons == 2)
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
		}
	}
	return isPassEventSelection;

}

// -- Event selecton for the electron channel (to make SelectedEE_t) (2018.08.02) -- // derived by Marijus Ambrozas
Bool_t DYAnalyzer::EventSelection_ElectronChannel(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in an event & NtupleHandle -- //
                                                vector< Electron >* SelectedElectronCollection,  // -- output: 2 electrons passing event selection conditions -- //
                                                vector< Int_t >* Sel_Index)    // -- output: 2 indexes of electrons that passed the selection -- //
{
        Bool_t isPassEventSelection = kFALSE;

        // -- Electron ID -- //
        vector< Electron > QElectronCollection;
        vector< Int_t > QIndex;

        for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
        {
            Electron elec = ElectronCollection[j];
            if(elec.passMediumID == kTRUE // modified by Dalmin Pai
                    && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
            {
                QElectronCollection.push_back(ElectronCollection[j]);
                QIndex.push_back(j);
            }
        }

        Int_t nQElectrons = (Int_t)QElectronCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();

        if(nQElectrons == 2 && nQIndices == 2)
        {
                Electron recolep1 = QElectronCollection[0];
                Electron recolep2 = QElectronCollection[1];

                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

                Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
                {
                        isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
                        Sel_Index->push_back(QIndex[0]);
                        Sel_Index->push_back(QIndex[1]);
                }
        }
        return isPassEventSelection;

}

// -- Event re-selecton for the electron channel (from LongSelectedEE_t) (2018.08.06) -- // derived by Marijus Ambrozas
Bool_t DYAnalyzer::EventSelection_ElectronChannel(vector< Electron > ElectronCollection, LongSelectedEE_t *ntuple, // -- input: All electrons in an event & LongSelectedEE_t -- //
                                                vector< Electron >* SelectedElectronCollection,  // -- output: 2 electrons passing event selection conditions -- //
                                                vector< Int_t >* Sel_Index)    // -- output: 2 indexes of electrons that passed the selection -- //
{
        Bool_t isPassEventSelection = kFALSE;

        // -- Electron ID -- //
        vector< Electron > QElectronCollection;
        vector< Int_t > QIndex;

        for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
        {
            Electron elec = ElectronCollection[j];
            if(elec.passMediumID == kTRUE // modified by Dalmin Pai
                    && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
            {
                QElectronCollection.push_back(ElectronCollection[j]);
                QIndex.push_back(j);
            }
        }

        Int_t nQElectrons = (Int_t)QElectronCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();

        if(nQElectrons == 2 && nQIndices == 2)
        {
                Electron recolep1 = QElectronCollection[0];
                Electron recolep2 = QElectronCollection[1];

                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

                Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
                {
                        isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
                        Sel_Index->push_back(QIndex[0]);
                        Sel_Index->push_back(QIndex[1]);
                }
        }
        return isPassEventSelection;

}

// -- N-1 cuts Event selecton for the electron channel -- // modified at 17 May 2017 by Dalmin Pai //
Bool_t DYAnalyzer::EventSelection_ElectronChannel1(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
                //if(elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1
                //if(elec.passMediumID == kTRUE // modified by Dalmin Pai
                if(elec.isMediumElectron_2016dataFor80X_minus_Full5x5_SigmaIEtaIEta() == kTRUE
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if(nQElectrons == 2)
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
                //if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel2(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
                //if(elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1
                //if(elec.passMediumID == kTRUE // modified by Dalmin Pai
                if(elec.isMediumElectron_2016dataFor80X_minus_dEtaInSeed() == kTRUE
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if(nQElectrons == 2)
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
                //if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel3(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
                //if(elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1
                //if(elec.passMediumID == kTRUE // modified by Dalmin Pai
                if(elec.isMediumElectron_2016dataFor80X_minus_dPhiIn() == kTRUE
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if(nQElectrons == 2)
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
                //if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel4(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
                //if(elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1
                //if(elec.passMediumID == kTRUE // modified by Dalmin Pai
                if(elec.isMediumElectron_2016dataFor80X_minus_HoverE() == kTRUE
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if(nQElectrons == 2)
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
                //if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel5(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
                //if(elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1
                //if(elec.passMediumID == kTRUE // modified by Dalmin Pai
                if(elec.isMediumElectron_2016dataFor80X_minus_RelPFIso_Rho() == kTRUE
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if(nQElectrons == 2)
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
                //if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel6(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
                //if(elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1
                //if(elec.passMediumID == kTRUE // modified by Dalmin Pai
                if(elec.isMediumElectron_2016dataFor80X_minus_InvEminusInvP() == kTRUE
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if(nQElectrons == 2)
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
                //if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel7(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
                //if(elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1
                //if(elec.passMediumID == kTRUE // modified by Dalmin Pai
                if(elec.isMediumElectron_2016dataFor80X_minus_mHits() == kTRUE
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if(nQElectrons == 2)
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
                //if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel8(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
                //if(elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1
                //if(elec.passMediumID == kTRUE // modified by Dalmin Pai
                if(elec.isMediumElectron_2016dataFor80X_minus_passConvVeto() == kTRUE
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if(nQElectrons == 2)
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
                //if(reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
		}
	}
	return isPassEventSelection;

}
// End of N-1 cuts Event Selection

// -- Event selecton for the electron channel (2016.02.11) -- //
Bool_t DYAnalyzer::EventSelection_ElectronChannel_NminusPFIso(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
                // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
                if(elec.isMediumElectron_Spring25ns_minus_PFIso() && elec.ecalDriven == 1
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !(fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566))
                        QElectronCollection.push_back(ElectronCollection[j]);
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if(nQElectrons == 2)
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if(reco_M > 10 && isPassAcc == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back(recolep1);
                        SelectedElectronCollection->push_back(recolep2);
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::isPassAccCondition_Electron(Electron Elec1, Electron Elec2)
{
	Bool_t isPassAcc = kFALSE;
	Electron leadElec, subElec;
	CompareElectron(&Elec1, &Elec2, &leadElec, &subElec);
        if(leadElec.Pt > LeadPtCut && fabs(leadElec.etaSC) < LeadEtaCut && !(fabs(leadElec.etaSC) > 1.4442 && fabs(leadElec.etaSC) < 1.566) &&
                subElec.Pt  > SubPtCut  && fabs(subElec.etaSC)  < SubEtaCut && !(fabs(subElec.etaSC) > 1.4442 && fabs(subElec.etaSC) < 1.566))
		isPassAcc = kTRUE;

	return isPassAcc;
}


Bool_t DYAnalyzer::isPassAccCondition_GenLepton_ECALGAP(GenLepton genlep1, GenLepton genlep2)
{
	Bool_t isPassAcc = kFALSE;

	GenLepton leadGenLep, subGenLep;
	CompareGenLepton(&genlep1, &genlep2, &leadGenLep, &subGenLep);
	
        if(leadGenLep.Pt > LeadPtCut && fabs(leadGenLep.eta) < LeadEtaCut && !(fabs(leadGenLep.eta) > 1.4442 && fabs(leadGenLep.eta) < 1.566) &&
                subGenLep.Pt  > SubPtCut  && fabs(subGenLep.eta) < SubEtaCut && !(fabs(subGenLep.eta) > 1.4442 && fabs(subGenLep.eta) < 1.566))
		isPassAcc = 1;

	return isPassAcc;
}

void DYAnalyzer::CompareElectron(Electron *Elec1, Electron *Elec2, Electron *leadElec, Electron *subElec)
{
    if(Elec1->Pt > Elec2->Pt)
    {
        *leadElec = *Elec1;
        *subElec = *Elec2;
    }
    else
    {
        *leadElec = *Elec2;
        *subElec = *Elec1;
    }
}

void DYAnalyzer::PostToPreFSR_byDressedLepton(NtupleHandle *ntuple, GenLepton *genlep_postFSR, Double_t dRCut, GenLepton *genlep_preFSR, vector< GenOthers >* GenPhotonCollection)
{
	TLorentzVector genlep_Mom_postFSR = genlep_postFSR->Momentum;

	TLorentzVector SumPhotonMom;
	SumPhotonMom.SetPxPyPzE(0,0,0,0);

	Int_t NGenOthers = ntuple->nGenOthers;
	for(Int_t i_gen=0; i_gen<NGenOthers; i_gen++)
	{
		GenOthers genlep;
		genlep.FillFromNtuple(ntuple, i_gen);

		// -- Only for the photons whose mother is muon or anti-muon -- //
                if(fabs(genlep.ID) == 22 && fabs(genlep.Mother) == 13)
		{
			
			Double_t dR = Calc_dR_GenLepton_GenOthers(*genlep_postFSR, genlep);

			// -- Sum of all photon's momentum near the post-FSR muon -- //
                        if(dR < dRCut)
			{
				SumPhotonMom  = SumPhotonMom + genlep.Momentum;
                                GenPhotonCollection->push_back(genlep);
			}
		}
	}

	// -- Momentum(pre-FSR) = Momentum(post-FSR) + Sum of all Photon's momentum near the post-FSR muon -- //
	genlep_preFSR->Momentum = genlep_Mom_postFSR + SumPhotonMom;
	genlep_preFSR->Et = genlep_preFSR->Momentum.Et();
	genlep_preFSR->Pt = genlep_preFSR->Momentum.Pt();
	genlep_preFSR->eta = genlep_preFSR->Momentum.Eta();
	genlep_preFSR->phi = genlep_preFSR->Momentum.Phi();
	genlep_preFSR->Px = genlep_preFSR->Momentum.Px();
	genlep_preFSR->Py = genlep_preFSR->Momentum.Py();
	genlep_preFSR->Pz = genlep_preFSR->Momentum.Pz();
}

void DYAnalyzer::PostToPreFSR_byDressedLepton_AllPhotons(NtupleHandle *ntuple, GenLepton *genlep_postFSR, Double_t dRCut, GenLepton *genlep_preFSR, vector< GenOthers >* GenPhotonCollection)
{
	TLorentzVector genlep_Mom_postFSR = genlep_postFSR->Momentum;

	TLorentzVector SumPhotonMom;
	SumPhotonMom.SetPxPyPzE(0,0,0,0);

	Int_t NGenOthers = ntuple->nGenOthers;
	for(Int_t i_gen=0; i_gen<NGenOthers; i_gen++)
	{
		GenOthers genlep;
		genlep.FillFromNtuple(ntuple, i_gen);

		// -- all photons within dR < 0.1 -- //
                // if(fabs(genlep.ID) == 22 && fabs(genlep.Mother) == 13)
                if(fabs(genlep.ID) == 22)
		{
			
			Double_t dR = Calc_dR_GenLepton_GenOthers(*genlep_postFSR, genlep);

			// -- Sum of all photon's momentum near the post-FSR muon -- //
                        if(dR < dRCut)
			{
				SumPhotonMom  = SumPhotonMom + genlep.Momentum;
                                GenPhotonCollection->push_back(genlep);
			}
		}
	}

	// -- Momentum(pre-FSR) = Momentum(post-FSR) + Sum of all Photon's momentum near the post-FSR muon -- //
	genlep_preFSR->Momentum = genlep_Mom_postFSR + SumPhotonMom;
	genlep_preFSR->Et = genlep_preFSR->Momentum.Et();
	genlep_preFSR->Pt = genlep_preFSR->Momentum.Pt();
	genlep_preFSR->eta = genlep_preFSR->Momentum.Eta();
	genlep_preFSR->phi = genlep_preFSR->Momentum.Phi();
	genlep_preFSR->Px = genlep_preFSR->Momentum.Px();
	genlep_preFSR->Py = genlep_preFSR->Momentum.Py();
	genlep_preFSR->Pz = genlep_preFSR->Momentum.Pz();
}

TString DYAnalyzer::DecideFSRType(GenLepton preFSR1, GenLepton preFSR2, GenLepton postFSR1, GenLepton postFSR2)
{
	TString FSRType = "";

	Bool_t isPassAcc_preFSREvent = kFALSE;
	isPassAcc_preFSREvent = isPassAccCondition_GenLepton(preFSR1, preFSR2);

	Bool_t isPassAcc_postFSREvent = kFALSE;
	isPassAcc_postFSREvent = isPassAccCondition_GenLepton(postFSR1, postFSR2);


        if(isPassAcc_preFSREvent == kFALSE && isPassAcc_postFSREvent == kTRUE)
		FSRType = "A";

        else if(isPassAcc_preFSREvent == kTRUE && isPassAcc_postFSREvent == kTRUE)
		FSRType = "B";
	
        else if(isPassAcc_preFSREvent == kTRUE && isPassAcc_postFSREvent == kFALSE)
		FSRType = "C";

        else if(isPassAcc_preFSREvent == kFALSE && isPassAcc_postFSREvent == kFALSE)
		FSRType = "D";
	else
	{
                std::cout << "ERROR: NO FSR TYPE CORRESPONDING TO THIS EVENT" << endl;
		FSRType = "NOTAssigned";
	}

	return FSRType;
}

Double_t DYAnalyzer::Calc_dR_GenLeptons(GenLepton genlep1, GenLepton genlep2)
{
	Double_t eta1 = genlep1.eta;
	Double_t phi1 = genlep1.phi;

	Double_t eta2 = genlep2.eta;
	Double_t phi2 = genlep2.phi;

	Double_t diff_eta = eta1 - eta2;
	Double_t diff_phi = phi1 - phi2;

        Double_t dR = sqrt(diff_eta * diff_eta + diff_phi * diff_phi);
	return dR;
}

Double_t DYAnalyzer::Calc_dR_GenLepton_GenOthers(GenLepton genlep1, GenOthers genlep2)
{
	Double_t eta1 = genlep1.eta;
	Double_t phi1 = genlep1.phi;

	Double_t eta2 = genlep2.eta;
	Double_t phi2 = genlep2.phi;

	Double_t diff_eta = eta1 - eta2;
	Double_t diff_phi = phi1 - phi2;

        Double_t dR = sqrt(diff_eta * diff_eta + diff_phi * diff_phi);
	return dR;
}

void DYAnalyzer::GenMatching(TString MuonType, NtupleHandle* ntuple, vector<Muon>* MuonCollection)
{
	vector<GenLepton> GenLeptonCollection;
	Int_t NGenLeptons = ntuple->gnpair;

        if(MuonType == "PromptFinalState")
	{
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
                        if(genlep.isMuon() && genlep.isPromptFinalState)
                                GenLeptonCollection.push_back(genlep);
		}
	}
        else if(MuonType == "fromTau")
	{
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
                        if(genlep.isMuon() && genlep.isDirectPromptTauDecayProductFinalState)
                                GenLeptonCollection.push_back(genlep);
		}

	}
        else if(MuonType == "fromHardProcess")
	{
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
                        if(genlep.isMuon() && genlep.fromHardProcessFinalState)
                                GenLeptonCollection.push_back(genlep);
		}
	}
        else if (MuonType == "decayedHadron")
        {
            for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
            {
                    GenLepton genlep;
                    genlep.FillFromNtuple(ntuple, i_gen);
                    if(genlep.isMuon() && genlep.isDecayedLeptonHadron)
                            GenLeptonCollection.push_back(genlep);
            }
        }
	else
	{
                std::cout << "Incorrect MuonType!" << endl;
		return;
	}

	//Give Acceptance Cuts
        if(GenLeptonCollection.size() >= 2)
	{
		GenLepton leadGenLep, subGenLep;
		CompareGenLepton(&GenLeptonCollection[0], &GenLeptonCollection[1], &leadGenLep, &subGenLep);
                if(!(leadGenLep.Pt > LeadPtCut && subGenLep.Pt > SubPtCut && abs(leadGenLep.eta) < LeadEtaCut && abs(subGenLep.eta) < SubEtaCut))
			GenLeptonCollection.clear();
	}


	
	Int_t NMuons = (Int_t)MuonCollection->size();
	vector<Muon> RecoMuonCollection;
	//Copy all muons in MuonCollection into RecoMuonCollection
	for(Int_t i_mu=0; i_mu<NMuons; i_mu++)
                RecoMuonCollection.push_back(MuonCollection->at(i_mu));

	MuonCollection->clear();

	Int_t NGen = (Int_t)GenLeptonCollection.size();
	for(Int_t i_gen=0; i_gen<NGen; i_gen++)
	{
		GenLepton genlep = GenLeptonCollection[i_gen];
		Double_t gen_Pt = genlep.Pt;
		Double_t gen_eta = genlep.eta;
		Double_t gen_phi = genlep.phi;

		Int_t i_matched = -1;
		Double_t dPtMin = 1e10;
		for(Int_t i_reco=0; i_reco<NMuons; i_reco++)
		{
			Muon mu = RecoMuonCollection[i_reco];
			Double_t reco_Pt = mu.Pt;
			Double_t reco_eta = mu.eta;
			Double_t reco_phi = mu.phi;

                        Double_t dR = sqrt((gen_eta-reco_eta)*(gen_eta-reco_eta) + (gen_phi-reco_phi)*(gen_phi-reco_phi));
			Double_t dPt = fabs(gen_Pt - reco_Pt);
                        if(dR < 0.3)
			{
                                if(dPt < dPtMin)
				{
					i_matched = i_reco;
					dPtMin = dPt;
				}
			}
		}

                if(i_matched != -1)
                        MuonCollection->push_back(RecoMuonCollection[i_matched]);
	}

	return;
}

void DYAnalyzer::ConvertToTunePInfo(Muon &mu)
{
	// -- Use TuneP information -- //
	mu.Pt = mu.TuneP_pT;
	mu.eta = mu.TuneP_eta;
	mu.phi = mu.TuneP_phi;

        mu.Momentum.SetPtEtaPhiM(mu.Pt, mu.eta, mu.phi, M_Mu);
}

void DYAnalyzer::PrintOutDoubleMuInfo(Muon mu1, Muon mu2)
{
	printf("\t[Muon1] (pT, eta, phi, charge) = (%10.5lf, %10.5lf, %10.5lf, %.1d)\n", mu1.Pt, mu1.eta, mu1.phi, mu1.charge);
	printf("\t[Muon2] (pT, eta, phi, charge) = (%10.5lf, %10.5lf, %10.5lf, %.1d)\n", mu2.Pt, mu2.eta, mu2.phi, mu2.charge);
        Double_t reco_M = (mu1.Momentum + mu2.Momentum).M();
	printf("\t\tDilepton Mass: %10.5lf\n", reco_M);

}

Double_t DYAnalyzer::GenMuonPt(TString muonType, NtupleHandle* ntuple, Muon reco_mu)
{
        Double_t GenMuonPt = 0;

        TString MuonType = muonType;
        MuonType.ToLower();

        vector<GenLepton> GenLeptonCollection;
        Int_t NGenLeptons = ntuple->gnpair;

        if(MuonType == "fromhardprocessfinalstate")
        {
                for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
                {
                        GenLepton genlep;
                        genlep.FillFromNtuple(ntuple, i_gen);
                        if(genlep.isMuon() && genlep.fromHardProcessFinalState)
                                GenLeptonCollection.push_back(genlep);
                }
        }
        else if(MuonType == "finalstate_or_hadrondecay")
        {
                for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
                {
                        GenLepton genlep;
                        genlep.FillFromNtuple(ntuple, i_gen);
                        if(genlep.isMuon() && (genlep.fromHardProcessFinalState || genlep.isDecayedLeptonHadron || genlep.fromHardProcessDecayed))
                                GenLeptonCollection.push_back(genlep);
                }
        }
        else
        {
                cout << "Incorrect MuonType!" << endl;
                return GenMuonPt;
        }

        //Give Acceptance Cuts
        if(GenLeptonCollection.size() >= 2)
        {
                GenLepton leadGenLep, subGenLep;
                CompareGenLepton(&GenLeptonCollection[0], &GenLeptonCollection[1], &leadGenLep, &subGenLep);
                if(!(leadGenLep.Pt > LeadPtCut && subGenLep.Pt > SubPtCut && abs(leadGenLep.eta) < LeadEtaCut && abs(subGenLep.eta) < SubEtaCut))
                        GenLeptonCollection.clear();
        }


        Double_t reco_Pt = reco_mu.Pt;
        Double_t reco_eta = reco_mu.eta;
        Double_t reco_phi = reco_mu.phi;
        Double_t dPtMin = 1e10;

        Int_t NGen = (Int_t)GenLeptonCollection.size();
        for(Int_t i_gen=0; i_gen<NGen; i_gen++)
        {
                GenLepton genlep = GenLeptonCollection[i_gen];
                Double_t gen_Pt = genlep.Pt;
                Double_t gen_eta = genlep.eta;
                Double_t gen_phi = genlep.phi;

                Double_t dR = sqrt((gen_eta-reco_eta)*(gen_eta-reco_eta) + (gen_phi-reco_phi)*(gen_phi-reco_phi));
                Double_t dPt = fabs(gen_Pt - reco_Pt);

                if(dR < 0.2)
                {
                        if(dPt < dPtMin)
                        {
                                GenMuonPt = gen_Pt;
                                dPtMin = dPt;
                        }
                }

        }
        if(fabs(GenMuonPt-reco_Pt)/reco_Pt > 0.3) GenMuonPt = 0; // To avoid very bad mismatches

        return GenMuonPt;
}


Bool_t DYAnalyzer::EventSelection_Dijet(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > PassingMuonCollection;
	vector< Muon > FailingMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_dzVTX())
	    {
                if(MuonCollection[j].trkiso < 0.10)
                        PassingMuonCollection.push_back(MuonCollection[j]);
	    	else
                        FailingMuonCollection.push_back(MuonCollection[j]);
	    }
	}

	Int_t nFailMuon = (Int_t)FailingMuonCollection.size();

        if(nFailMuon >= 2) // -- Dijet events: contains more than 2 failing muons regardless of # passing muons -- //
	{
                if(nFailMuon == 2)
		{
			Muon recolep1 = FailingMuonCollection[0];
			Muon recolep2 = FailingMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
                        Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

			Bool_t isOS = kFALSE;
                        if(recolep1.charge != recolep2.charge) isOS = kTRUE;

                        // if(reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005)
                        if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE)
			{
				isPassEventSelection = kTRUE;
                                SelectedMuonCollection->push_back(recolep1);
                                SelectedMuonCollection->push_back(recolep2);
			}

                } // -- end of if(nFailMuon == 2) -- //
		else // -- # failing muons > 2 -- // 
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nFailMuon; i_mu++)
			{
				Muon Mu = FailingMuonCollection[i_mu];

				// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
				for(Int_t j_mu=0; j_mu<nFailMuon; j_mu++)
				{
					Muon Mu_jth = FailingMuonCollection[j_mu];

                                        if(j_mu != i_mu) // -- do not calculate vertex variables(prob, chi2). with itself -- //
					{
						// -- Check that this pair is within acceptance -- //
						Bool_t isPassAcc = kFALSE;
						isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

                                                if(isPassAcc == kTRUE) // -- Find best pair ONLY for the pairs within acceptance -- //
						{
							Double_t VtxProb_temp = -999;
							Double_t VtxNormChi2_temp = 999;
							DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

							// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
                                                        if(VtxNormChi2_temp < VtxNormChi2_BestPair)
							{
								VtxNormChi2_BestPair = VtxNormChi2_temp;
								mu1_BestPair = Mu;
								mu2_BestPair = Mu_jth;
							}
						}
					}
				} // -- end of the loop for j_mu (finding for second muon)
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

                        if(VtxNormChi2_BestPair < 999) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
                                Double_t Angle = reco_v1.Angle(reco_v2.Vect());

				Bool_t isOS = kFALSE;
                                if(mu1_BestPair.charge != mu2_BestPair.charge) isOS = kTRUE;

                                if(reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE)
				{
					isPassEventSelection = kTRUE;
                                        SelectedMuonCollection->push_back(mu1_BestPair);
                                        SelectedMuonCollection->push_back(mu2_BestPair);
				}
			}

		} // -- end of (# failing muons > 2) case -- //

        } // -- end of if(nFailMuon >= 2) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Wjet(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > PassingMuonCollection;
	vector< Muon > FailingMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_dzVTX())
	    {
                if(MuonCollection[j].trkiso < 0.10)
                        PassingMuonCollection.push_back(MuonCollection[j]);
	    	else
                        FailingMuonCollection.push_back(MuonCollection[j]);
	    }
	}

	Int_t nFailMuon = (Int_t)FailingMuonCollection.size();
	Int_t nPassMuon = (Int_t)PassingMuonCollection.size();

        if(nFailMuon == 1 && nPassMuon == 1) // -- W+Jets events: exactly (# pass muon , # fail muon) = (1, 1) -- //
	{
		Muon recolep1 = PassingMuonCollection[0]; // -- first one: passing muon -- //
		Muon recolep2 = FailingMuonCollection[1]; // -- second one: failing muon -- //

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

		Bool_t isOS = kFALSE;
                if(recolep1.charge != recolep2.charge) isOS = kTRUE;

                // if(reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005)
                if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE)
		{
			isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back(recolep1); // -- first one: passing muon -- //
                        SelectedMuonCollection->push_back(recolep2); // -- second one: failing muon -- //
		}
	}

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_CheckMoreThanOneDimuonCand(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection, Bool_t& isMoreThanOneCand) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;
	isMoreThanOneCand = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
            if(MuonCollection[j].isHighPtMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
                QMuonCollection.push_back(MuonCollection[j]);
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
                if(mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*"))
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

        if(isExistHLTMatchedMuon == kTRUE)
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
                if(nQMuons == 2)
		{
			Muon recolep1 = QMuonCollection[0];
			Muon recolep2 = QMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
                        Double_t Angle = recolep1.Momentum.Angle(recolep2.Momentum.Vect());

			Bool_t isOS = kFALSE;
                        if(recolep1.charge != recolep2.charge) isOS = kTRUE;

                        // if(reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005)
                        if(reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE)
			{
				isPassEventSelection = kTRUE;
                                SelectedMuonCollection->push_back(recolep1);
                                SelectedMuonCollection->push_back(recolep2);
			}
		}
                else if(nQMuons > 2)
		{
			isMoreThanOneCand = kTRUE;
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
                                if(Mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*"))
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

                                                if(j_mu != i_mu) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

                                                        if(isPassAcc == kTRUE) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
                                                                if(VtxNormChi2_temp < VtxNormChi2_BestPair)
								{
									VtxNormChi2_BestPair = VtxNormChi2_temp;
									mu1_BestPair = Mu;
									mu2_BestPair = Mu_jth;
								}
							}
						}
					} // -- end of the loop for j_mu (finding for second muon)
				}
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

                        if(VtxNormChi2_BestPair < 999) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
                                Double_t Angle = reco_v1.Angle(reco_v2.Vect());

				Bool_t isOS = kFALSE;
                                if(mu1_BestPair.charge != mu2_BestPair.charge) isOS = kTRUE;

                                if(reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE)
				{
					isPassEventSelection = kTRUE;
                                        SelectedMuonCollection->push_back(mu1_BestPair);
                                        SelectedMuonCollection->push_back(mu2_BestPair);
				}
			}

                } // -- End of else if(nQMuons > 2) -- //

        } // -- End of if(isExistHLTMatchedMuon == kTRUE) -- //

	return isPassEventSelection;
}


Int_t DYAnalyzer::Find_muon_PtBin_Reco(Double_t Pt)
{
        const Int_t nPtBins = 1;
        Double_t PtBinEdges[nPtBins+1] = {25, 500};

        Int_t ptbin = 9999;

        // -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- //
        if(Pt >= PtBinEdges[nPtBins])
                ptbin = nPtBins-1;
        // -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
        else if(Pt <= PtBinEdges[0])
                ptbin = 0;
        else
        {
                for(Int_t i=0; i<nPtBins; i++)
                {
                        if(Pt >= PtBinEdges[i] && Pt < PtBinEdges[i+1])
                        {
                                ptbin = i;
                                break;
                        }
                }
        }

        return ptbin;
}


Int_t DYAnalyzer::Find_muon_PtBin_ID(Double_t Pt)
{
    // -- HighPtID & TrkIso -- //
    //const Int_t nPtBins = 7;
    //Double_t PtBinEdges[nPtBins+1] = {20, 25, 30, 40, 50, 55, 60, 120};

    // -- TightID & PFIso -- //
    const Int_t nPtBins = 6;
    Double_t PtBinEdges[nPtBins+1] = {20, 25, 30, 40, 50, 60, 120};

    Int_t ptbin = 9999;

    // -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- //
    if(Pt > PtBinEdges[nPtBins])
            ptbin = nPtBins-1;
    // -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
    else if(Pt <= PtBinEdges[0])
            ptbin = 0;
    else
    {
        for(Int_t i=0; i<nPtBins; i++)
        {
            if(Pt > PtBinEdges[i] && Pt <= PtBinEdges[i+1])
            {
                ptbin = i;
                break;
            }
        }
    }

    return ptbin;
}


Int_t DYAnalyzer::Find_muon_PtBin_Iso(Double_t Pt)
{
        // -- HighPtID & TrkIso -- //
        //const Int_t nPtBins = 7;
        //Double_t PtBinEdges[nPtBins+1] = {20, 25, 30, 40, 50, 55, 60, 120};

        // -- TightID & PFIso -- //
        const Int_t nPtBins = 6;
        Double_t PtBinEdges[nPtBins+1] = {20, 25, 30, 40, 50, 60, 120};

        Int_t ptbin = 9999;

        // -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- //
        if(Pt > PtBinEdges[nPtBins])
                ptbin = nPtBins-1;
        // -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
        else if(Pt <= PtBinEdges[0])
                ptbin = 0;
        else
        {
                for(Int_t i=0; i<nPtBins; i++)
                {
                        if(Pt > PtBinEdges[i] && Pt <= PtBinEdges[i+1])
                        {
                                ptbin = i;
                                break;
                        }
                }
        }

        return ptbin;
}


Int_t DYAnalyzer::Find_muon_PtBin_Trig(Double_t Pt)
{
    // -- IsoMu24_OR_IsoTkMu24 -- //
    const Int_t nPtBins = 7;
    Double_t PtBinEdges[nPtBins+1] = {26, 30, 40, 50, 60, 120, 200, 500};

    Int_t ptbin = 9999;

    // -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- //
    if(Pt > PtBinEdges[nPtBins])
            ptbin = nPtBins-1;
    // -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
    else if(Pt <= PtBinEdges[0])
            ptbin = 0;
    else
    {
        for(Int_t i=0; i<nPtBins; i++)
        {
            if(Pt > PtBinEdges[i] && Pt <= PtBinEdges[i+1])
            {
                ptbin = i;
                break;
            }
        }
    }

    return ptbin;
}


Int_t DYAnalyzer::Find_muon_PtBin_Trig_new(Double_t Pt)
{
    // -- IsoMu24_OR_IsoTkMu24 -- //
    const Int_t nPtBins = 8;
    Double_t PtBinEdges[nPtBins+1] = {26, 30, 40, 50, 60, 80, 120, 200, 500};

    Int_t ptbin = 9999;

    // -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- //
    if(Pt > PtBinEdges[nPtBins])
            ptbin = nPtBins-1;
    // -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- //
    else if(Pt <= PtBinEdges[0])
            ptbin = 0;
    else
    {
        for(Int_t i=0; i<nPtBins; i++)
        {
            //if(Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1])
            if(Pt > PtBinEdges[i] && Pt <= PtBinEdges[i+1])
            {
                ptbin = i;
                break;
            }
        }
    }

    return ptbin;
}


Int_t DYAnalyzer::Find_muon_EtaBin_Reco(Double_t eta)
{
        const Int_t nEtaBins = 15;
        Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4};

        Int_t etabin = 9999;

        for(Int_t i=0; i<nEtaBins; i++)
        {
                if(eta >= EtaBinEdges[i] && eta < EtaBinEdges[i+1])
                //if(fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1])
                {
                        etabin = i;
                        break;
                }
        }

        return etabin;
}


Int_t DYAnalyzer::Find_muon_EtaBin_ID(Double_t eta)
{
        //const Int_t nEtaBins = 5;
        const Int_t nEtaBins = 4;
        //Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
        Double_t EtaBinEdges[nEtaBins+1] = {0, 0.9, 1.2, 2.1, 2.4};

        Int_t etabin = 9999;

        for(Int_t i=0; i<nEtaBins; i++)
        {
                //if(eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1])
                if(fabs(eta) >= EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1])
                {
                        etabin = i;
                        break;
                }
        }

        return etabin;
}


Int_t DYAnalyzer::Find_muon_EtaBin_Iso(Double_t eta)
{
        //const Int_t nEtaBins = 5;
        const Int_t nEtaBins = 4;
        //Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
        Double_t EtaBinEdges[nEtaBins+1] = {0, 0.9, 1.2, 2.1, 2.4};

        Int_t etabin = 9999;

        for(Int_t i=0; i<nEtaBins; i++)
        {
                //if(eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1])
                if(fabs(eta) >= EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1])
                {
                        etabin = i;
                        break;
                }
        }

        return etabin;
}


Int_t DYAnalyzer::Find_muon_EtaBin_Trig(Double_t eta)
{
        //const Int_t nEtaBins = 5;
        const Int_t nEtaBins = 4;
        //Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
        Double_t EtaBinEdges[nEtaBins+1] = {0, 0.9, 1.2, 2.1, 2.4};

        Int_t etabin = 9999;

        for(Int_t i=0; i<nEtaBins; i++)
        {
                //if(eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1])
                if(fabs(eta) >= EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1])
                {
                        etabin = i;
                        break;
                }
        }

        return etabin;
}


Int_t DYAnalyzer::Find_electron_PtBin_Reco(Double_t Pt)
{
        //const Int_t nPtBins = 4;
        const Int_t nPtBins = 1;
        //Double_t PtBinEdges[nPtBins+1] = {10, 22, 40, 70, 250};
        Double_t PtBinEdges[nPtBins+1] = {25, 500};

        Int_t ptbin = 9999;

        // -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- //
        if(Pt > PtBinEdges[nPtBins])
                ptbin = nPtBins-1;
        // -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
        else if(Pt < PtBinEdges[0])
                ptbin = 0;
        else
        {
                for(Int_t i=0; i<nPtBins; i++)
                {
                        if(Pt >= PtBinEdges[i] && Pt < PtBinEdges[i+1])
                        {
                                ptbin = i;
                                break;
                        }
                }
        }

        return ptbin;
}


Int_t DYAnalyzer::Find_electron_PtBin_ID(Double_t Pt)
{
        //const Int_t nPtBins = 4;
        const Int_t nPtBins = 5;
        //Double_t PtBinEdges[nPtBins+1] = {10, 22, 40, 70, 250};
        Double_t PtBinEdges[nPtBins+1] = {10, 20, 35, 50, 90, 150};

        Int_t ptbin = 9999;

        // -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- //
        if(Pt > PtBinEdges[nPtBins])
                ptbin = nPtBins-1;
        // -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
        else if(Pt < PtBinEdges[0])
                ptbin = 0;
        else
        {
                for(Int_t i=0; i<nPtBins; i++)
                {
                        if(Pt >= PtBinEdges[i] && Pt < PtBinEdges[i+1])
                        {
                                ptbin = i;
                                break;
                        }
                }
        }

        return ptbin;
}


Int_t DYAnalyzer::Find_electron_PtBin_Trig(Double_t Pt)
{
        //const Int_t nPtBins = 4;
        const Int_t nPtBins = 7;
        //Double_t PtBinEdges[nPtBins+1] = {10, 22, 40, 70, 250};
        Double_t PtBinEdges[nPtBins+1] = {10, 15, 20, 25, 30, 50, 90, 150};

        Int_t ptbin = 9999;

        // -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- //
        if(Pt > PtBinEdges[nPtBins])
                ptbin = nPtBins-1;
        // -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
        else if(Pt < PtBinEdges[0])
                ptbin = 0;
        else
        {
                for(Int_t i=0; i<nPtBins; i++)
                {
                        if(Pt >= PtBinEdges[i] && Pt < PtBinEdges[i+1])
                        {
                                ptbin = i;
                                break;
                        }
                }
        }

        return ptbin;
}


Int_t DYAnalyzer::Find_electron_EtaBin_Reco(Double_t eta)
{
        //const Int_t nEtaBins = 5;
        const Int_t nEtaBins = 30;
        //Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
        Double_t EtaBinEdges[nEtaBins+1] = {-2.5,-2.45,-2.4,-2.3,-2.2,-2.0,-1.8,-1.63,-1.566,-1.4442,-1.2,-1.0,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,1.0,1.2,1.4442,1.566,1.63,1.8,2.0,2.2,2.3,2.4,2.45,2.5};

        Int_t etabin = 9999;

        for(Int_t i=0; i<nEtaBins; i++)
        {
                if(eta >= EtaBinEdges[i] && eta < EtaBinEdges[i+1])
                //if(fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1])
                {
                        etabin = i;
                        break;
                }
        }

        return etabin;
}


Int_t DYAnalyzer::Find_electron_EtaBin_ID(Double_t eta)
{
        //const Int_t nEtaBins = 5;
        //const Int_t nEtaBins = 10;
        const Int_t nEtaBins = 22;
        //Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
        //Double_t EtaBinEdges[nEtaBins+1] = {-2.5, -2.0, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.0, 2.5};
        Double_t EtaBinEdges[nEtaBins+1] = {-2.4,-2.3,-2.2,-2.1,-2.0,-1.566,-1.4442,-1.1,-0.8,-0.4,-0.2,0.0,0.2,0.4,0.8,1.1,1.4442,1.566,2.0,2.1,2.2,2.3,2.4};

        Int_t etabin = 9999;

        for(Int_t i=0; i<nEtaBins; i++)
        {
                if(eta >= EtaBinEdges[i] && eta < EtaBinEdges[i+1])
                //if(fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1])
                {
                        //etabin = i;
                        if(eta < -1.566) etabin = i;
                        else if(-1.4442 < eta && eta < 1.4442) etabin = i-1;
                        else if(1.4442 < eta) etabin = i-2;

                        break;
                }
        }

        return etabin;
}


Int_t DYAnalyzer::Find_electron_EtaBin_Trig(Double_t eta)
{
        //const Int_t nEtaBins = 5;
        //const Int_t nEtaBins = 10;
        const Int_t nEtaBins = 22;
        //Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
        //Double_t EtaBinEdges[nEtaBins+1] = {-2.5, -2.0, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.0, 2.5};
        Double_t EtaBinEdges[nEtaBins+1] = {-2.4,-2.3,-2.2,-2.1,-2.0,-1.566,-1.4442,-1.1,-0.8,-0.4,-0.2,0.0,0.2,0.4,0.8,1.1,1.4442,1.566,2.0,2.1,2.2,2.3,2.4};

        Int_t etabin = 9999;

        for(Int_t i=0; i<nEtaBins; i++)
        {
                if(eta >= EtaBinEdges[i] && eta < EtaBinEdges[i+1])
                //if(fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1])
                {
                        //etabin = i;
                        if(eta < -1.566) etabin = i;
                        else if(-1.4442 < eta && eta < 1.4442) etabin = i-1;
                        else if(1.4442 < eta) etabin = i-2;

                        break;
                }
        }

        return etabin;
}


void DYAnalyzer::emuVertexProbNormChi2(NtupleHandle *ntuple, Double_t ele_Pt, Double_t mu_Pt, Double_t *VtxProb, Double_t *VtxNormChi2)
{
        vector<double> *PtCollection1 = ntuple->vtxTrkEMu1Pt; //electron
        vector<double> *PtCollection2 = ntuple->vtxTrkEMu2Pt; //muon
        vector<double> *VtxProbCollection = ntuple->vtxTrkEMuProb;

        Int_t NPt1 = (Int_t)PtCollection1->size();
        Int_t NPt2 = (Int_t)PtCollection2->size();
        Int_t NProb = (Int_t)VtxProbCollection->size();

        if(NPt1 != NPt2 || NPt2 != NProb || NPt1 != NProb)
                cout << "NPt1: " << NPt1 << " NPt2: " << NPt2 << " Nprob: " << NProb << endl;

        for(Int_t i=0; i<NProb; i++)
        {
                if(PtCollection1->at(i) == ele_Pt && PtCollection2->at(i) == mu_Pt)
                {
                        *VtxProb = VtxProbCollection->at(i);
                        *VtxNormChi2 = ntuple->vtxTrkEMuChi2->at(i) / ntuple->vtxTrkEMuNdof->at(i);
                        break;
                }
        }

        return;
}

void DYAnalyzer::emuVertexProbNormChi2(LongSelectedEMu_t *ntuple, Double_t ele_Pt, Double_t mu_Pt, Double_t *VtxProb, Double_t *VtxNormChi2)
{
        vector<double> *PtCollection1 = ntuple->vtxTrkEMu1Pt; //electron
        vector<double> *PtCollection2 = ntuple->vtxTrkEMu2Pt; //muon
        vector<double> *VtxProbCollection = ntuple->vtxTrkEMuProb;

        Int_t NPt1 = (Int_t)PtCollection1->size();
        Int_t NPt2 = (Int_t)PtCollection2->size();
        Int_t NProb = (Int_t)VtxProbCollection->size();

        if(NPt1 != NPt2 || NPt2 != NProb || NPt1 != NProb)
                cout << "NPt1: " << NPt1 << " NPt2: " << NPt2 << " Nprob: " << NProb << endl;

        for(Int_t i=0; i<NProb; i++)
        {
                if(PtCollection1->at(i) == ele_Pt && PtCollection2->at(i) == mu_Pt)
                {
                        *VtxProb = VtxProbCollection->at(i);
                        *VtxNormChi2 = ntuple->vtxTrkEMuChi2->at(i) / ntuple->vtxTrkEMuNdof->at(i);
                        break;
                }
        }

        return;
}
