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
#include <iostream>

#define Lumi 35867 // -- from Run2016B to Run2016H, JSON. unit: /pb, Updated at 2017.07.30 -- //
#define Lumi_HLTv4p2 865.919 // -- integrated luminosity before Run 257933 -- //
#define nMassBin 43

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

	Double_t Eff_RecoID_data[5][4];
	Double_t Eff_RecoID_MC[5][4];

	Double_t Eff_Iso_data[5][4];
	Double_t Eff_Iso_MC[5][4];

	Double_t Eff_HLTv4p2_data[5][4];
	Double_t Eff_HLTv4p2_MC[5][4];

	Double_t Eff_HLTv4p3_data[5][4];
	Double_t Eff_HLTv4p3_MC[5][4];

	// -- For efficiency SF of BtoF -- //
	Double_t Eff_RecoID_data_BtoF[4][7];
	Double_t Eff_RecoID_MC_BtoF[4][7];

	Double_t Eff_Iso_data_BtoF[4][7];
	Double_t Eff_Iso_MC_BtoF[4][7];

	Double_t Eff_HLT_data_BtoF[4][7];
	Double_t Eff_HLT_MC_BtoF[4][7];

	// -- For efficiency SF of GtoH -- //
	Double_t Eff_RecoID_data_GtoH[4][7];
	Double_t Eff_RecoID_MC_GtoH[4][7];

	Double_t Eff_Iso_data_GtoH[4][7];
	Double_t Eff_Iso_MC_GtoH[4][7];

	Double_t Eff_HLT_data_GtoH[4][7];
	Double_t Eff_HLT_MC_GtoH[4][7];

	// -- For efficiency SF of electron -- //
	Double_t Eff_Reco_data[30][1];
	Double_t Eff_Reco_MC[30][1];

	Double_t Eff_ID_data[10][5];
	Double_t Eff_ID_MC[10][5];

	// -- Constructor -- //
	DYAnalyzer(TString HLTname);

	// -- Setup accetpance cuts -- //
	void AssignAccThreshold(TString HLTname, TString *HLT, Double_t *LeadPtCut, Double_t *SubPtCut, Double_t *LeadEtaCut, Double_t *SubEtaCut);

	////////////////////////////
	// -- Setup MC samples -- //
	////////////////////////////
	void SetupMCsamples_v20160412_76X_MINIAODv2_CheckPremix( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents );
	void SetupMCsamples_Moriond17( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents );
	void SetupMCsamples_v20160309_76X_MiniAODv2( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents );
	void SetupMCsamples_v20160131_MiniAODv2( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents );
	void SetupMCsamples_v20160117_MiniAOD_JetMET( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents );
        void SetupDataSamples( TString Type, TString DataType, vector<TString> *ntupleDirectory, vector<TString> *Tag);
        Bool_t SeparateDYLLSample_isHardProcess(TString Tag, NtupleHandle *ntuple);
	Bool_t Separate_ttbarSample(TString Tag, NtupleHandle *ntuple, vector<GenOthers> *GenTopCollection);

	// -- outdated -- //
	Bool_t SeparateDYLLSample(TString Tag, NtupleHandle *ntuple);

	//////////////////////////////////
	// -- Setup pileup weighting -- //
	//////////////////////////////////
	void SetupPileUpReWeighting( Bool_t isMC );
	Double_t PileUpWeightValue(Int_t PileUp_MC);

	// -- for 76X -- //
	void SetupPileUpReWeighting_76X( Bool_t isMC );
	Double_t PileUpWeightValue_76X(Int_t PileUp_MC);

	// -- for 80X -- //
	void SetupPileUpReWeighting_80X( Bool_t isMC, TString ROOTFileName );
	Double_t PileUpWeightValue_80X(Int_t PileUp_MC);


	/////////////////////////////////////////
	// -- Setup Efficiency scale factor -- //
	/////////////////////////////////////////
	void SetupEfficiencyScaleFactor();
	void SetupEfficiencyScaleFactor(TString ROOTFileName);
	void SetupEfficiencyScaleFactor_BtoF();
	void SetupEfficiencyScaleFactor_GtoH();
	void SetupEfficiencyScaleFactor_electron();
	Double_t EfficiencySF_EventWeight_HLTv4p2(Muon mu1, Muon mu2);
	Double_t EfficiencySF_EventWeight_HLTv4p3(Muon mu1, Muon mu2);
	Double_t EfficiencySF_EventWeight_HLT_BtoF(Muon mu1, Muon mu2);
	Double_t EfficiencySF_EventWeight_HLT_GtoH(Muon mu1, Muon mu2);
	Double_t EfficiencySF_EventWeight_electron(Electron ele1, Electron ele2);
	Int_t FindPtBin(Double_t Pt);
	Int_t FindEtaBin(Double_t eta);
	Int_t FindPtBin_trig(Double_t Pt);
	Int_t FindEtaBin_trig(Double_t eta);
	Int_t FindPtBin_Reco(Double_t Pt);
	Int_t FindEtaBin_Reco(Double_t eta);
	Int_t FindPtBin_ID(Double_t Pt);
	Int_t FindEtaBin_ID(Double_t eta);

	// -- outdated -- //
	Double_t EfficiencySF_EventWeight(Muon mu1, Muon mu2, NtupleHandle *ntuple);
	Double_t EfficiencySF_EventWeight_RecoIdIso(Muon mu1, Muon mu2, NtupleHandle *ntuple);
	
	////////////////////////////
	// -- Event Selections -- //
	////////////////////////////
	Bool_t EventSelection(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
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
	Double_t Calc_dR_GenLeptons( GenLepton genlep1, GenLepton genlep2 );
	Double_t Calc_dR_GenLepton_GenOthers( GenLepton genlep1, GenOthers genlep2 );

	// -- miscellaneous -- //
	void GenMatching(TString MuonType, NtupleHandle* ntuple, vector<Muon>* MuonCollection);
	void ConvertToTunePInfo( Muon &mu );
	void PrintOutDoubleMuInfo( Muon mu1, Muon mu2 );

	// -- emu method -- //
        Bool_t EventSelection_emu_method_test(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple,
                                              vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection); // -- output: 1 muon and 1 electron passing event selection conditions -- //
        Bool_t EventSelection_emu_method_test(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple, // Derived by Marijus Ambrozas 2018.08.07
                                              vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection,
                                              Int_t &Sel_Index_Mu, Int_t &Sel_Index_Ele); // -- output: 1 muon and 1 electron passing event selection conditions and their indices -- //
        Bool_t EventSelection_emu_method_test(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, LongSelectedEMu_t *ntuple, // Derived by Marijus Ambrozas 2018.08.07
                                              vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection,
                                              Int_t &Sel_Index_Mu, Int_t &Sel_Index_Ele); // -- output: 1 muon and 1 electron passing event selection conditions and their indices -- //
//	Bool_t EventSelection_emu_method(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection); // -- output: 1 muon and 1 electron passing event selection conditions -- //
//	void emuVertexProbNormChi2(NtupleHandle *ntuple, Double_t Pt1, Double_t Pt2, Double_t *VtxProb, Double_t *VtxNormChi2);
	Double_t EfficiencySF_EventWeight_emu_BtoF(Muon mu1, Electron ele2);
	Double_t EfficiencySF_EventWeight_emu_GtoH(Muon mu1, Electron ele2);
};

DYAnalyzer::DYAnalyzer(TString HLTname)
{
	if( HLTname == "None" )
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
	if( HLTname == "IsoMu20" )
	{
		*HLT = "HLT_IsoMu20_v*";
		*LeadPtCut = 22;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "IsoMu20_OR_IsoTkMu20" )
	{
		*HLT = "HLT_IsoMu20_v* || HLT_IsoTkMu20_v*";
		*LeadPtCut = 22;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "IsoMu24" )
	{
		*HLT = "HLT_IsoMu24_v*";
		*LeadPtCut = 26;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "IsoMu24_OR_IsoTkMu24" ) // added at 2017.08.01
	{
		*HLT = "HLT_IsoMu24_v* || HLT_IsoTkMu24_v*";
		//*LeadPtCut = 26;
		//*SubPtCut = 10;
		*LeadPtCut = 28;
		*SubPtCut = 17;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "Mu45_eta2p1" )
	{
		*HLT = "HLT_Mu45_eta2p1_v*";
		*LeadPtCut = 46;
		*SubPtCut = 10;
		*LeadEtaCut = 2.1;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "Mu50" )
	{
		*HLT = "HLT_Mu50_v*";
		*LeadPtCut = 53;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "IsoMu20_SymmetricPt25" )
	{
		*HLT = "HLT_IsoMu20_v*";
		*LeadPtCut = 25;
		*SubPtCut = 25;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "Ele17Ele12" )
	{
		*HLT = "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
		*LeadPtCut = 25;
		*SubPtCut = 15;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
	else if( HLTname == "Ele22_eta2p1" )
	{
		*HLT = "HLT_Ele22_eta2p1_WPLoose_Gsf_v*"; // -- Exist only for the data; "HLT_Ele22_eta2p1_WP75_Gsf_v*" should be used for MC
		*LeadPtCut = 25;
		*SubPtCut = 15;
		*LeadEtaCut = 2.1;
		*SubEtaCut = 2.1;
	}
	else if( HLTname == "Ele22_eta2p1_NoEtaCut" )
	{
		*HLT = "HLT_Ele22_eta2p1_WPLoose_Gsf_v*"; // -- Exist only for the data; "HLT_Ele22_eta2p1_WP75_Gsf_v*" should be used for MC
		*LeadPtCut = 25;
		*SubPtCut = 15;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
	else if( HLTname == "Pt_30_10_eta_2p5" )
	{
		*HLT = "None"; // -- just for acceptance test -- //
		*LeadPtCut = 30;
		*SubPtCut = 10;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
	else if( HLTname == "Ele23_WPLoose" )
	{
		*HLT = "HLT_Ele23_WPLoose_Gsf_v*"; // -- Exist only for the data; "HLT_Ele22_eta2p1_WP75_Gsf_v*" should be used for MC
		*LeadPtCut = 30;
		*SubPtCut = 10;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
	else if( HLTname == "Ele23Ele12" ) // updated at 2017.02.21 by Dalmin Pai
	{
		*HLT = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
		*LeadPtCut = 28;
		*SubPtCut = 17;
		//*LeadEtaCut = 2.5; // -- Later it should exclude ECAL gap
		//*SubEtaCut = 2.5; // -- Later it should exclude ECAL gap
		*LeadEtaCut = 2.4; // -- Later it should exclude ECAL gap
		*SubEtaCut = 2.4; // -- Later it should exclude ECAL gap
	}
	else
	{ 
                std::cout << "Wrong HLT name!: " << HLTname << endl;
		return; 
	}

}

void DYAnalyzer::SetupMCsamples_v20160412_76X_MINIAODv2_CheckPremix( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents )
{
	if( Type == "DYMuMu_PU25" )
	{
		ntupleDirectory->push_back( "Premix/v20160412_76X_MINIAODv2_CheckPremix_CorrectDataSetName_DYJets_Classic_PU25" ); Tag->push_back( "DYMuMu_M50_PU25_Classic" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 8426438.0 );
		ntupleDirectory->push_back( "Premix/v20160412_76X_MINIAODv2_CheckPremix_CorrectDataSetName_DYJets_NonDeterministic_PU25" ); Tag->push_back( "DYMuMu_M50_PU25_NonDet" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 8286714.0 );
	}
	else if( Type == "DYEE_PU25" )
	{
		ntupleDirectory->push_back( "Premix/v20160412_76X_MINIAODv2_CheckPremix_CorrectDataSetName_DYJets_Classic_PU25" ); Tag->push_back( "DYEE_M50_PU25_Classic" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 8439605.0 );
		ntupleDirectory->push_back( "Premix/v20160412_76X_MINIAODv2_CheckPremix_CorrectDataSetName_DYJets_NonDeterministic_PU25" ); Tag->push_back( "DYEE_M50_PU25_NonDet" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 8300442.0 );
	}
}

void DYAnalyzer::SetupMCsamples_Moriond17( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents )
{
	if( Type == "DYMuMu_M10to50" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000" );
                Tag->push_back( "DYMuMu_M10to50_v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000" );
                Tag->push_back( "DYMuMu_M10to50_v2" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000" );
                Tag->push_back( "DYMuMu_M10to50_ext1v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "DYMuMu_M50to100" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
                ntupleDirectory->push_back( "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000" );
                Tag->push_back( "DYMuMu_M50to100" ); Xsec->push_back( 5869.58346/3.0 ); nEvents->push_back( 26175605.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "DYMuMu_M100toInf" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
                ntupleDirectory->push_back( "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000" );
                Tag->push_back( "DYMuMu_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3433295.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000" );
                Tag->push_back( "DYMuMu_M100to200_ext" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3433295.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000" );
                Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 56340.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000" );
                Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000" );
                Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 48188.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000" );
                Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44984.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000" );
                Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000" );
                Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40110.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000" );
                Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000" );
                Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 33360.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "DYMuMu_aMCNLO" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000" );
                Tag->push_back( "DYMuMu_M10to50_v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000" );
                Tag->push_back( "DYMuMu_M10to50_v2" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000" );
                Tag->push_back( "DYMuMu_M10to50_ext1v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000" );
                Tag->push_back( "DYMuMu_M50to100" ); Xsec->push_back( 5869.58346/3.0 ); nEvents->push_back( 26175605.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000" );
                Tag->push_back( "DYMuMu_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3433295.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000" );
                Tag->push_back( "DYMuMu_M100to200_ext" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3433295.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000" );
                Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 56340.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000" );
                Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000" );
                Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 48188.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000" );
                Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44984.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000" );
                Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000" );
                Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40110.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000" );
                Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights

                ntupleDirectory->push_back( "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000" );
                Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 33360.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "DYEE_M10to50" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000" );
                Tag->push_back( "DYEE_M10to50_v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000" );
                Tag->push_back( "DYEE_M10to50_v2" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000" );
                Tag->push_back( "DYEE_M10to50_ext1v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights
	}
	else if( Type == "DYEE_M50to100" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
                ntupleDirectory->push_back( "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000" );
                Tag->push_back( "DYEE_M50to100" ); Xsec->push_back( 5869.58346/3.0 ); nEvents->push_back( 26166194.0 ); //nEvents: sum of DYEE weights
	}
	else if( Type == "DYEE_M100toInf" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
                ntupleDirectory->push_back( "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000" );
                Tag->push_back( "DYEE_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3437885.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000" );
                Tag->push_back( "DYEE_M100to200_ext" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3437885.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000" );
                Tag->push_back( "DYEE_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 56144.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000" );
                Tag->push_back( "DYEE_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50420.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000" );
                Tag->push_back( "DYEE_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 48039.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000" );
                Tag->push_back( "DYEE_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 46114.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000" );
                Tag->push_back( "DYEE_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 44256.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000" );
                Tag->push_back( "DYEE_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 39712.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000" );
                Tag->push_back( "DYEE_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37287.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000" );
                Tag->push_back( "DYEE_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 34031.0 ); //nEvents: sum of DYEE weights
	}
	else if( Type == "DYEE_aMCNLO" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000" );
                Tag->push_back( "DYEE_M10to50_v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000" );
                Tag->push_back( "DYEE_M10to50_v2" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000" );
                Tag->push_back( "DYEE_M10to50_ext1v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000" );
                Tag->push_back( "DYEE_M50to100" ); Xsec->push_back( 5869.58346/3.0 ); nEvents->push_back( 26166194.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200/180326_143238/0000" );
                Tag->push_back( "DYEE_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3437885.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200_ext/180326_143324/0000" );
                Tag->push_back( "DYEE_M100to200_ext" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3437885.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M200to400/180326_143408/0000" );
                Tag->push_back( "DYEE_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 56144.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M400to500/180326_143512/0000" );
                Tag->push_back( "DYEE_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50420.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M500to700/180326_143600/0000" );
                Tag->push_back( "DYEE_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 48039.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M700to800/180326_143640/0000" );
                Tag->push_back( "DYEE_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 46114.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M800to1000/180326_143747/0000" );
                Tag->push_back( "DYEE_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 44256.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1000to1500/180326_143836/0000" );
                Tag->push_back( "DYEE_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 39712.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M1500to2000/180326_143921/0000" );
                Tag->push_back( "DYEE_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37287.0 ); //nEvents: sum of DYEE weights

                ntupleDirectory->push_back( "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M2000to3000/180326_144005/0000" );
                Tag->push_back( "DYEE_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 34031.0 ); //nEvents: sum of DYEE weights
	}
	else if( Type == "ttbar" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		//ntupleDirectory->push_back( "ttbar" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 154948878.0 ); //ttbar+ ttbarBackup
                ntupleDirectory->push_back( "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar/180326_142926/0000" );
                Tag->push_back( "ttbar" ); Xsec->push_back( 734.577 ); nEvents->push_back( 135949780.0 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
	}
	else if( Type == "ttbarBackup" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		//ntupleDirectory->push_back( "ttbarBackup" ); Tag->push_back( "ttbarBackup" ); Xsec->push_back( 831.76 ); nEvents->push_back( 154948878.0 ); //ttbar+ ttbarBackup
                ntupleDirectory->push_back( "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbarBackup/180326_143005/0000" );
                Tag->push_back( "ttbarBackup" ); Xsec->push_back( 734.577 ); nEvents->push_back( 135949780.0 ); //M(ttbar) < 700GeV, ttbar+ttbarBackup
	}
	else if( Type == "ttbar_M700toInf" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
                ntupleDirectory->push_back( "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M700to1000/180326_143059/0000" );
                Tag->push_back( "ttbar_M700to1000" ); Xsec->push_back( 76.605 ); nEvents->push_back( 38422582.0 ); //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)

                ntupleDirectory->push_back( "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar_M1000toInf/180326_143144/0000" );
                Tag->push_back( "ttbar_M1000toInf" ); Xsec->push_back( 20.578 ); nEvents->push_back( 24561630.0 ); //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
	}
	else if( Type == "DYTauTau_M10to50" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v1/180326_142925/0000" );
                Tag->push_back( "DYTauTau_M10to50_v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33080379.0 ); //nEvents: sum of DYTauTau weights

                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_v2/180326_143001/0000" );
                Tag->push_back( "DYTauTau_M10to50_v2" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33080379.0 ); //nEvents: sum of DYTauTau weights

                ntupleDirectory->push_back( "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_ext1v1/180326_143056/0000" );
                Tag->push_back( "DYTauTau_M10to50_ext1v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33080379.0 ); //nEvents: sum of DYTauTau weights
	}
	else if( Type == "DYTauTau_M50toInf" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		//ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104.0/3.0 ); nEvents->push_back( 27277866.0 ); //nEvents: sum of DYTauTau weights
                ntupleDirectory->push_back( "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M50toInf/180326_143143/0000" );
                Tag->push_back( "DYTauTau" ); Xsec->push_back( 1921.8 ); nEvents->push_back( 27277866.0 ); //nEvents: sum of DYTauTau weights, NNLO Xsec
	}
	else if( Type == "VVnST" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
                ntupleDirectory->push_back( "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW/180326_143800/0000" );
                Tag->push_back( "tW" ); Xsec->push_back( 35.85 ); nEvents->push_back( 6952830.0 );

                ntupleDirectory->push_back( "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tbarW/180326_143849/0000" );
                Tag->push_back( "tbarW" ); Xsec->push_back( 35.85 ); nEvents->push_back( 6933093.0 );

                ntupleDirectory->push_back( "ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/180326_143627/0000" ); // NOT SURE (there also is ZZTo4L)
                Tag->push_back( "ZZ" ); Xsec->push_back( 16.523 ); nEvents->push_back( 998034.0 );

                ntupleDirectory->push_back( "WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/180326_143414/0000" ); // NOT SURE (there also is WZTo3LNu)
                Tag->push_back( "WZ" ); Xsec->push_back( 47.13 ); nEvents->push_back( 2995828.0 );

                ntupleDirectory->push_back( "WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/180326_143237/0000" ); // NOT SURE (there also is WWTo2L2Nu)
                Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 6987123.0 );
	}
	else if( Type == "WJetsToLNu" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
                ntupleDirectory->push_back( "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo/180326_144617/0000" );
                Tag->push_back( "WJetsToLNu" ); Xsec->push_back( 61526.7 ); nEvents->push_back( 86731698.0 );

                ntupleDirectory->push_back( "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo_ext/180326_144652/0000" ); // There also is madgraph version
                Tag->push_back( "WJetsToLNu_ext" ); Xsec->push_back( 61526.7 ); nEvents->push_back( 86731698.0 );
	}
        else if( Type == "WJetsToLNu_test" )
        {
                // -- Background Samples -- //
                ntupleDirectory->push_back( "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_ext/180326_143105/0000" );
                Tag->push_back( "WJetsToLNu_test" ); Xsec->push_back( 61526.7 ); nEvents->push_back( 86731698.0 );
        }
	else if( Type == "QCDMuEnriched" )
	{
                std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
                ntupleDirectory->push_back( "QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt15to20/180326_143059/0000" );
                Tag->push_back( "QCDMuEnriched_Pt15to20" ); Xsec->push_back( 720648000*0.00042 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt20to30/180326_143144/0000" );
                Tag->push_back( "QCDMuEnriched_Pt20to30" ); Xsec->push_back( 1273190000*0.003 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt30to50/180326_143240/0000" );
                Tag->push_back( "QCDMuEnriched_Pt30to50" ); Xsec->push_back( 139803000*0.01182 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt50to80/180326_143340/0000" );
                Tag->push_back( "QCDMuEnriched_Pt50to80" ); Xsec->push_back( 19222500*0.02276 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120/180326_143419/0000" );
                Tag->push_back( "QCDMuEnriched_Pt80to120" ); Xsec->push_back( 2758420*0.03844 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt80to120_ext1/180326_143533/0000" );
                Tag->push_back( "QCDMuEnriched_Pt80to120_ext1" ); Xsec->push_back( 2758420*0.03844  ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170/180326_143612/0000" );
                Tag->push_back( "QCDMuEnriched_Pt120to170" ); Xsec->push_back( 469797*0.05362 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt120to170_backup/180326_143654/0000" );
                Tag->push_back( "QCDMuEnriched_Pt120to170_backup" ); Xsec->push_back( 469797*0.05362  ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300/180326_143750/0000" );
                Tag->push_back( "QCDMuEnriched_Pt170to300" ); Xsec->push_back( 117989*0.07335 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_ext1/180326_143849/0000" );
                Tag->push_back( "QCDMuEnriched_Pt170to300_ext1" ); Xsec->push_back( 117989*0.07335 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt170to300_backup/180326_143946/0000" );
                Tag->push_back( "QCDMuEnriched_Pt170to300_backup" ); Xsec->push_back( 117989*0.07335 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470/180326_144021/0000" );
                Tag->push_back( "QCDMuEnriched_Pt300to470" ); Xsec->push_back( 7820.25*0.10196 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext1/180326_144117/0000" );
                Tag->push_back( "QCDMuEnriched_Pt300to470_ext1" ); Xsec->push_back( 7820.25*0.10196 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt300to470_ext2/180326_144211/0000" );
                Tag->push_back( "QCDMuEnriched_Pt300to470_ext2" ); Xsec->push_back( 7820.25*0.10196 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600/180326_144301/0000" );
                Tag->push_back( "QCDMuEnriched_Pt470to600" ); Xsec->push_back( 645.528*0.12242 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt470to600_ext1/180326_144358/0000" );
                Tag->push_back( "QCDMuEnriched_Pt470to600_ext1" ); Xsec->push_back( 645.528*0.12242 ); nEvents->push_back( 1.0 );

//		ntupleDirectory->push_back( "QCDMuEnriched_Pt470to600_ext2" );      // DID NOT FIND THIS ONE
//                Tag->push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec->push_back( 645.528*0.12242 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800/180326_144534/0000" );
                Tag->push_back( "QCDMuEnriched_Pt600to800" ); Xsec->push_back( 187.109*0.13412 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_ext1/180326_144612/0000" );
                Tag->push_back( "QCDMuEnriched_Pt600to800_ext1" ); Xsec->push_back( 187.109*0.13412 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt600to800_backup/180326_144648/0000" );
                Tag->push_back( "QCDMuEnriched_Pt600to800_backup" ); Xsec->push_back( 187.109*0.13412 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000/180326_144736/0000" );
                Tag->push_back( "QCDMuEnriched_Pt800to1000" ); Xsec->push_back( 32.3486*0.14552 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext1/180326_144818/0000" );
                Tag->push_back( "QCDMuEnriched_Pt800to1000_ext1" ); Xsec->push_back( 32.3486*0.14552 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt800to1000_ext2/180326_144856/0000" );
                Tag->push_back( "QCDMuEnriched_Pt800to1000_ext2" ); Xsec->push_back( 32.3486*0.14552 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf/180326_144937/0000" );
                Tag->push_back( "QCDMuEnriched_Pt1000toInf" ); Xsec->push_back( 10.4305*0.15544 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/crab_QCDMuEnriched_Pt1000toInf_ext1/180326_145024/0000" );
                Tag->push_back( "QCDMuEnriched_Pt1000toInf_ext1" ); Xsec->push_back( 10.4305*0.15544 ); nEvents->push_back( 1.0 );
	}
	else if( Type == "QCDEMEnriched" )
	{
                std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
                ntupleDirectory->push_back( "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt20to30/180326_145104/0000" );
                Tag->push_back( "QCDEMEnriched_Pt20to30" ); Xsec->push_back( 557600000*0.0096 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50/180326_145144/0000" );
                Tag->push_back( "QCDEMEnriched_Pt30to50" ); Xsec->push_back( 136000000*0.073 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50_ext1/180326_145227/0000" );
                Tag->push_back( "QCDEMEnriched_Pt30to50_ext1" ); Xsec->push_back( 136000000*0.073 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80/180326_145308/0000" );
                Tag->push_back( "QCDEMEnriched_Pt50to80" ); Xsec->push_back( 19800000*0.146 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80_ext1/180326_145353/0000" );
                Tag->push_back( "QCDEMEnriched_Pt50to80_ext1" ); Xsec->push_back( 19800000*0.146 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120/180326_145437/0000" );
                Tag->push_back( "QCDEMEnriched_Pt80to120" ); Xsec->push_back( 2800000*0.125 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120_ext1/180326_145522/0000" );
                Tag->push_back( "QCDEMEnriched_Pt80to120_ext1" ); Xsec->push_back( 2800000*0.125 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170/180326_145602/0000" );
                Tag->push_back( "QCDEMEnriched_Pt120to170" ); Xsec->push_back( 477000*0.132 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170_ext1/180326_145701/0000" );
                Tag->push_back( "QCDEMEnriched_Pt120to170_ext1" ); Xsec->push_back( 477000*0.132 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt170to300/180326_145738/0000" );
                Tag->push_back( "QCDEMEnriched_Pt170to300" ); Xsec->push_back( 114000*0.165 ); nEvents->push_back( 1.0 );

                ntupleDirectory->push_back( "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt300toInf/180326_145836/0000" );
                Tag->push_back( "QCDEMEnriched_Pt300toInf" ); Xsec->push_back( 9000*0.15 ); nEvents->push_back( 1.0 );
	}
	else
                std::cout << "Wrong Type!" << endl;
}

void DYAnalyzer::SetupMCsamples_v20160309_76X_MiniAODv2( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents )
{
	if( Type == "Full" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 985598 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 999996 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 988416 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16520811.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7467514.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6309713.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 97994304 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6302525.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO" )
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6302525.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "M100to200" )
	{
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M100to200_25ns" ); Tag->push_back( "DYMuMu_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 227522.0 ); //nEvents: sum of weights within 10<M<50
	}
	else if( Type == "M50to200" )
	{
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6302525.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "M50toInf" )
	{
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50toInf" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6311695.0 );
	}
	else if( Type == "M10to50_M50toInf" )
	{
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50toInf" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6311695.0 );
	}
	else if( Type == "aMCNLO_M120Cut" )
	{
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to120" ); Xsec->push_back( 1975 ); nEvents->push_back( 6243307.0 );
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M100to200_25ns" ); Tag->push_back( "DYMuMu_M120to200" ); Xsec->push_back( 19.32 ); nEvents->push_back( 55554.0 );
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 2.731 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights
	}

	else if( Type == "Full_Include_M100to200" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 985598 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 999996 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 988416 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16520811.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7467514.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6309713.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 97994304 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to100" ); Xsec->push_back( 5869.58346/3.0 ); nEvents->push_back( 6061181.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M100to200_25ns" ); Tag->push_back( "DYMuMu_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 227522.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "Full_NoHighMass" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 985598 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 999996 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 988416 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16520811.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7467514.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6309713.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 97994304 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50toInf" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6311695.0 ); // -- sum of weight should be updated! -- //
	}
	else if( Type == "Full_Powheg" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 985598 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 999996 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 988416 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16520811.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7467514.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6309713.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 97994304 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160519_76X_MINIAODv2_Resubmit4_AdjustRunTime_ZMuMuPowheg_M50to120_25ns" ); Tag->push_back( "ZMuMu_M50to120" );  Xsec->push_back(1975);  nEvents->push_back(2971982.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M120to200_25ns" ); Tag->push_back( "ZMuMu_M120to200" );  Xsec->push_back(19.32);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M200to400_25ns" ); Tag->push_back( "ZMuMu_M200to400" );  Xsec->push_back(2.731);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M400to800_25ns" ); Tag->push_back( "ZMuMu_M400to800" );  Xsec->push_back(0.241);  nEvents->push_back(99600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M800to1400_25ns" );  Tag->push_back( "ZMuMu_M800to1400" );  Xsec->push_back(0.01678);  nEvents->push_back(97600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M1400to2300_25ns" );  Tag->push_back( "ZMuMu_M1400to2300" );  Xsec->push_back(0.00139);  nEvents->push_back(99200.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M2300to3500_25ns" );  Tag->push_back( "ZMuMu_M2300to3500" );  Xsec->push_back(0.00008948);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M3500to4500_25ns" );  Tag->push_back( "ZMuMu_M3500to4500" );  Xsec->push_back(0.000004135);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M4500to6000_25ns" );  Tag->push_back( "ZMuMu_M4500to6000" );  Xsec->push_back(4.56E-07);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M6000toInf_25ns" );  Tag->push_back( "ZMuMu_M6000toInf" );  Xsec->push_back(2.066E-08);  nEvents->push_back(99200.0);
	}
	else if( Type == "Full_M120Cut" )
	{
                // std::cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 985598 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 999996 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 988416 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16520811.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7467514.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6309713.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 97994304 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to120" ); Xsec->push_back( 1975 ); nEvents->push_back( 6243307.0 );
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M100to200_25ns" ); Tag->push_back( "DYMuMu_M120to200" ); Xsec->push_back( 19.32 ); nEvents->push_back( 55554.0 );
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 2.731 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO_Include_M100to200")
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to100" ); Xsec->push_back( 5869.58346/3.0 ); nEvents->push_back( 6061181.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M100to200_25ns" ); Tag->push_back( "DYMuMu_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 227522.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "Madgraph" )
	{
		ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M5to50_25ns" ); Tag->push_back( "Madgraph_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 631905.0 );
		ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50toInf" ); Xsec->push_back( 6014/3.0 ); nEvents->push_back( 3003455.0 );

		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M5to50_25ns" ); Tag->push_back( "Madgraph_M5to50" ); Xsec->push_back( 7160/3.0 ); nEvents->push_back( 2782834.0 );
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50toInf" ); Xsec->push_back( 4895/3.0 ); nEvents->push_back( 3003455.0 );

		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50to150" ); Xsec->push_back( (4895 - 6.58)/3.0 ); nEvents->push_back( 3003455.0 ); 
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M150toInf_25ns" ); Tag->push_back( "Madgraph_M150toInf" ); Xsec->push_back( 6.58/3.0 ); nEvents->push_back( 1.0 );

	}
	else if( Type == "MadgraphPowheg" ) // -- for estimation of syst. from unfolding -- //
	{
		ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M5to50_25ns" ); Tag->push_back( "Madgraph_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 631905.0 );
		ntupleDirectory->push_back( "76X/v20160519_76X_MINIAODv2_Resubmit4_AdjustRunTime_ZMuMuPowheg_M50to120_25ns" ); Tag->push_back( "ZMuMu_M50to120" );  Xsec->push_back(1975);  nEvents->push_back(2971982.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M120to200_25ns" ); Tag->push_back( "ZMuMu_M120to200" );  Xsec->push_back(19.32);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M200to400_25ns" ); Tag->push_back( "ZMuMu_M200to400" );  Xsec->push_back(2.731);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M400to800_25ns" ); Tag->push_back( "ZMuMu_M400to800" );  Xsec->push_back(0.241);  nEvents->push_back(99600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M800to1400_25ns" );  Tag->push_back( "ZMuMu_M800to1400" );  Xsec->push_back(0.01678);  nEvents->push_back(97600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M1400to2300_25ns" );  Tag->push_back( "ZMuMu_M1400to2300" );  Xsec->push_back(0.00139);  nEvents->push_back(99200.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M2300to3500_25ns" );  Tag->push_back( "ZMuMu_M2300to3500" );  Xsec->push_back(0.00008948);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M3500to4500_25ns" );  Tag->push_back( "ZMuMu_M3500to4500" );  Xsec->push_back(0.000004135);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M4500to6000_25ns" );  Tag->push_back( "ZMuMu_M4500to6000" );  Xsec->push_back(4.56E-07);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M6000toInf_25ns" );  Tag->push_back( "ZMuMu_M6000toInf" );  Xsec->push_back(2.066E-08);  nEvents->push_back(99200.0);


		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M5to50_25ns" ); Tag->push_back( "Madgraph_M5to50" ); Xsec->push_back( 7160/3.0 ); nEvents->push_back( 2782834.0 );
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50toInf" ); Xsec->push_back( 4895/3.0 ); nEvents->push_back( 3003455.0 );

		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50to150" ); Xsec->push_back( (4895 - 6.58)/3.0 ); nEvents->push_back( 3003455.0 ); 
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M150toInf_25ns" ); Tag->push_back( "Madgraph_M150toInf" ); Xsec->push_back( 6.58/3.0 ); nEvents->push_back( 1.0 );

	}
	else if( Type == "Powheg" ) // -- for estimation of syst. from unfolding -- //
	{
		ntupleDirectory->push_back( "76X/v20160519_76X_MINIAODv2_Resubmit4_AdjustRunTime_ZMuMuPowheg_M50to120_25ns" ); Tag->push_back( "ZMuMu_M50to120" );  Xsec->push_back(1975);  nEvents->push_back(2971982.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M120to200_25ns" ); Tag->push_back( "ZMuMu_M120to200" );  Xsec->push_back(19.32);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M200to400_25ns" ); Tag->push_back( "ZMuMu_M200to400" );  Xsec->push_back(2.731);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M400to800_25ns" ); Tag->push_back( "ZMuMu_M400to800" );  Xsec->push_back(0.241);  nEvents->push_back(99600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M800to1400_25ns" );  Tag->push_back( "ZMuMu_M800to1400" );  Xsec->push_back(0.01678);  nEvents->push_back(97600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M1400to2300_25ns" );  Tag->push_back( "ZMuMu_M1400to2300" );  Xsec->push_back(0.00139);  nEvents->push_back(99200.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M2300to3500_25ns" );  Tag->push_back( "ZMuMu_M2300to3500" );  Xsec->push_back(0.00008948);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M3500to4500_25ns" );  Tag->push_back( "ZMuMu_M3500to4500" );  Xsec->push_back(0.000004135);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M4500to6000_25ns" );  Tag->push_back( "ZMuMu_M4500to6000" );  Xsec->push_back(4.56E-07);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M6000toInf_25ns" );  Tag->push_back( "ZMuMu_M6000toInf" );  Xsec->push_back(2.066E-08);  nEvents->push_back(99200.0);


		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M5to50_25ns" ); Tag->push_back( "Madgraph_M5to50" ); Xsec->push_back( 7160/3.0 ); nEvents->push_back( 2782834.0 );
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50toInf" ); Xsec->push_back( 4895/3.0 ); nEvents->push_back( 3003455.0 );

		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50to150" ); Xsec->push_back( (4895 - 6.58)/3.0 ); nEvents->push_back( 3003455.0 ); 
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M150toInf_25ns" ); Tag->push_back( "Madgraph_M150toInf" ); Xsec->push_back( 6.58/3.0 ); nEvents->push_back( 1.0 );

	}
	else
                std::cout << "Wrong Type!" << endl;
}

void DYAnalyzer::SetupMCsamples_v20160131_MiniAODv2( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents )
{
	if( Type == "Full" )
	{
		// -- Background Samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 996944 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 978512 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 993640 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16541203.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7255646.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6419292.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 19757182 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6413327.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "Full_NoHighMass" )
	{
		// -- Background Samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 996944 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 978512 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 993640 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16541203.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7255646.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6419292.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 19757182 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6422093.0 ); //nEvents: sum of DYMuMu weights
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "M50_M200to400" )
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6422093.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO" )
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6413327.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "Powheg" )
	{
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M50to120_25ns" ); Tag->push_back( "ZMuMu_M50to120" );  Xsec->push_back(1975);  nEvents->push_back(2836871);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M120to200_25ns" ); Tag->push_back( "ZMuMu_M120to200" );  Xsec->push_back(19.32);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M200to400_25ns" ); Tag->push_back( "ZMuMu_M200to400" );  Xsec->push_back(2.731);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M400to800_25ns" ); Tag->push_back( "ZMuMu_M400to800" );  Xsec->push_back(0.241);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M800to1400_25ns" );  Tag->push_back( "ZMuMu_M800to1400" );  Xsec->push_back(0.01678);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M1400to2300_25ns" );  Tag->push_back( "ZMuMu_M1400to2300" );  Xsec->push_back(0.00139);  nEvents->push_back(99600);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M2300to3500_25ns" );  Tag->push_back( "ZMuMu_M2300to3500" );  Xsec->push_back(0.00008948);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M3500to4500_25ns" );  Tag->push_back( "ZMuMu_M3500to4500" );  Xsec->push_back(0.000004135);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M4500to6000_25ns" );  Tag->push_back( "ZMuMu_M4500to6000" );  Xsec->push_back(4.56E-07);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M6000toInf_25ns" );  Tag->push_back( "ZMuMu_M6000toInf" );  Xsec->push_back(2.066E-08);  nEvents->push_back(100000);
	}
	else if( Type == "Full_withoutM200to400" )
	{
		// -- Background Samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 996944 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 978512 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 993640 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16541203.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7255646.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6419292.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 19757182 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to400" ); Xsec->push_back( 6103.25346/3.0 ); nEvents->push_back( 6413327.0 ); //nEvents: sum of DYMuMu weights
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO_withoutM200to400" )
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to400" ); Xsec->push_back( 6103.25346/3.0 ); nEvents->push_back( 6413327.0 ); //nEvents: sum of DYMuMu weights
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO_M50toInf" )
	{
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50toInf" ); Xsec->push_back( 2008.4 ); nEvents->push_back( 6422093.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "aMCNLO_M200to400" )
	{
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO_DYEE" )
	{
                // std::cout << "Warning: # events should be adjusted using Sum weights of DYEE events (current one: DYMuMu SumWeights)" << endl;
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYEE_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7.29361e+06 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYEE_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6.40938e+06 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYEE_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18348 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYEE_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17410 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYEE_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17245 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYEE_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 16120 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYEE_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14397 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYEE_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 13857 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYEE_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13495 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYEE_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12859 ); //nEvents: sum of DYEE weights 
	}
	else if( Type == "aMCNLO_FEWZxSec" )
	{
		// xSec of M10-50 and M50 sample: aMC@NLO -- //
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6413327.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 2.59583 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.136235 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.0775862 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.0121251 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.010281 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.00546713 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.000735022 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.000176089 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else
                std::cout << "Wrong Type!" << endl;

	return;
}

void DYAnalyzer::SetupMCsamples_v20160117_MiniAOD_JetMET( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents )
{
	if( Type == "Full" )
	{
		// -- Background Samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 996168 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 991232 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 994416 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16518173 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7418362 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6430407 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 19899492 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7418362 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6430407 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18339 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 6951 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO" )
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7418362 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6430407 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18339 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 6951 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "Powheg" )
	{
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M50to120_25ns" ); Tag->push_back( "ZMuMu_M50to120" );  Xsec->push_back(1975);  nEvents->push_back(2848071);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M120to200_25ns" ); Tag->push_back( "ZMuMu_M120to200" );  Xsec->push_back(19.32);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M200to400_25ns" ); Tag->push_back( "ZMuMu_M200to400" );  Xsec->push_back(2.731);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M400to800_25ns" ); Tag->push_back( "ZMuMu_M400to800" );  Xsec->push_back(0.241);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M800to1400_25ns" );  Tag->push_back( "ZMuMu_M800to1400" );  Xsec->push_back(0.01678);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M1400to2300_25ns" );  Tag->push_back( "ZMuMu_M1400to2300" );  Xsec->push_back(0.00139);  nEvents->push_back(99600);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M2300to3500_25ns" );  Tag->push_back( "ZMuMu_M2300to3500" );  Xsec->push_back(0.00008948);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M3500to4500_25ns" );  Tag->push_back( "ZMuMu_M3500to4500" );  Xsec->push_back(0.000004135);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M4500to6000_25ns" );  Tag->push_back( "ZMuMu_M4500to6000" );  Xsec->push_back(4.56E-07);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M6000toInf_25ns" );  Tag->push_back( "ZMuMu_M6000toInf" );  Xsec->push_back(2.066E-08);  nEvents->push_back(100000);
	}
	else
                std::cout << "Wrong Type!" << endl;

	return;
}

void DYAnalyzer::SetupDataSamples(TString Type, TString DataType, vector<TString> *ntupleDirectory, vector<TString> *Tag)
{
    if ( Type == "Data" )
    {
        if( DataType == "DoubleEG_B" )
        {
            ntupleDirectory->push_back( "DoubleEG/crab_DoubleEG_RunB/180326_143532/0000" ); Tag->push_back( "DoubleEG_B_0000" );
            ntupleDirectory->push_back( "DoubleEG/crab_DoubleEG_RunB/180326_143532/0001" ); Tag->push_back( "DoubleEG_B_0001" );
        }
        else if( DataType == "DoubleEG_C" )
        {
            ntupleDirectory->push_back( "DoubleEG/crab_DoubleEG_RunC/180326_143612/0000" ); Tag->push_back( "DoubleEG_C" );
        }
        else if( DataType == "DoubleEG_D" )
        {
            ntupleDirectory->push_back( "DoubleEG/crab_DoubleEG_RunD/180326_143654/0000" ); Tag->push_back( "DoubleEG_D" );
        }
        else if( DataType == "DoubleEG_E" )
        {
            ntupleDirectory->push_back( "DoubleEG/crab_DoubleEG_RunE/180326_143750/0000" ); Tag->push_back( "DoubleEG_E" );
        }
        else if( DataType == "DoubleEG_F" )
        {
            ntupleDirectory->push_back( "DoubleEG/crab_DoubleEG_RunF/180326_143846/0000" ); Tag->push_back( "DoubleEG_F" );
        }
        else if( DataType == "DoubleEG_G" )
        {
            ntupleDirectory->push_back( "DoubleEG/crab_DoubleEG_RunG/180326_144559/0000" ); Tag->push_back( "DoubleEG_G_0000" );
            ntupleDirectory->push_back( "DoubleEG/crab_DoubleEG_RunG/180326_144559/0001" ); Tag->push_back( "DoubleEG_G_0001" );
        }
        else if( DataType == "DoubleEG_H" )
        {
            ntupleDirectory->push_back( "DoubleEG/crab_DoubleEG_RunHver2/180326_144638/0000" ); Tag->push_back( "DoubleEG_Hver2_0000" );
            ntupleDirectory->push_back( "DoubleEG/crab_DoubleEG_RunHver2/180326_144638/0001" ); Tag->push_back( "DoubleEG_Hver2_0001" );
            ntupleDirectory->push_back( "DoubleEG/crab_DoubleEG_RunHver3/180326_144719/0000" ); Tag->push_back( "DoubleEG_Hver3" );
        }
        else if( DataType == "SingleMuon_B" )
        {
            ntupleDirectory->push_back( "SingleMuon/crab_SingleMuon_RunB/180326_143105/0000" ); Tag->push_back( "SingleMuon_B" );
        }
        else if( DataType == "SingleMuon_C" )
        {
            ntupleDirectory->push_back( "SingleMuon/crab_SingleMuon_RunC/180326_143152/0000" ); Tag->push_back( "SingleMuon_C" );
        }
        else if( DataType == "SingleMuon_D" )
        {
            ntupleDirectory->push_back( "SingleMuon/crab_SingleMuon_RunD/180326_143257/0000" ); Tag->push_back( "SingleMuon_D" );
        }
        else if( DataType == "SingleMuon_E" )
        {
            ntupleDirectory->push_back( "SingleMuon/crab_SingleMuon_RunE/180326_143338/0000" ); Tag->push_back( "SingleMuon_E" );
        }
        else if( DataType == "SingleMuon_F" )
        {
            ntupleDirectory->push_back( "SingleMuon/crab_SingleMuon_RunF/180326_143419/0000" ); Tag->push_back( "SingleMuon_F" );
        }
        else if( DataType == "SingleMuon_G" )
        {
            ntupleDirectory->push_back( "SingleMuon/crab_SingleMuon_RunG/180326_144335/0000" ); Tag->push_back( "SingleMuon_G" );
        }
        else if( DataType == "SingleMuon_H" )
        {
            ntupleDirectory->push_back( "SingleMuon/crab_SingleMuon_RunHver2/180326_144412/0000" ); Tag->push_back( "SingleMuon_Hver2" );
            ntupleDirectory->push_back( "SingleMuon/crab_SingleMuon_RunHver3/180326_144454/0000" ); Tag->push_back( "SingleMuon_Hver3" );
        }
        else cout << "Wrong DataType!" << endl;
    }
    else cout << "Wrong Type!" << endl;
}

Bool_t DYAnalyzer::Separate_ttbarSample(TString Tag, NtupleHandle *ntuple, vector<GenOthers> *GenTopCollection)
{
	Bool_t GenFlag = kFALSE;

	// -- Seperate ttbar events -- //
	if( Tag.Contains("ttbar") )
	{
		vector<GenOthers> GenOthersCollection;
		Int_t NGenOthers = ntuple->nGenOthers;
		for(Int_t i_gen=0; i_gen<NGenOthers; i_gen++)
		{
			GenOthers genothers;
			genothers.FillFromNtuple(ntuple, i_gen);
			if( abs(genothers.ID) == 6 && genothers.isHardProcess )
				GenOthersCollection.push_back( genothers );
		}

		if( GenOthersCollection.size() == 2 ) // -- Select the ttbar events from hard-process -- //
		{
			// -- Check top & anti-top pair -- //
			if( GenOthersCollection[0].ID == GenOthersCollection[1].ID )
				printf("%d %d\n", GenOthersCollection[0].ID, GenOthersCollection[1].ID);

			//if( Tag == "ttbar" ) // -- Select only evetns withtin M < 700 -- //
			if( Tag == "ttbar" || Tag == "ttbarBackup" ) // -- Select only evetns withtin M < 700 -- //
			{
				TLorentzVector v1 = GenOthersCollection[0].Momentum;
				TLorentzVector v2 = GenOthersCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 700 )
				//if( reco_M > -999 )
				{
					GenFlag = kTRUE;
					GenTopCollection->push_back( GenOthersCollection[0] );
					GenTopCollection->push_back( GenOthersCollection[1] );
				}
			}
			else // ex: ttbar_M700to1000, ttbar_M1000toInf
			{
				GenFlag = kTRUE;
				GenTopCollection->push_back( GenOthersCollection[0] );
				GenTopCollection->push_back( GenOthersCollection[1] );
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
	if( Tag.Contains("DYMuMu") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isMuon() && genlep.isHardProcess )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 muons from hard-process -- //
		{
			if( Tag == "DYMuMu_M50to200" ) // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 200 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYMuMu_M50to400" ) // -- Select only evetns withtin 50 < M < 400 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 400 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYMuMu_M50to100" || Tag == "DYMuMu_Photos_M50to100" ) // -- Select only evetns withtin 50 < M < 100 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 100 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYMuMu_M50to120" )
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 120 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYMuMu_M120to200" )
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M > 120 )
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
	else if( Tag.Contains("DYEE") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isElectron() && genlep.isHardProcess )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 electrons from hard-process -- //
		{
			if( Tag == "DYEE_M50to200" ) // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 200 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYEE_M50to100" ) // -- Select only evetns withtin 50 < M < 100 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 100 )
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
	// -- Separate DYTauTau events from MuMu events -- //
	else if( Tag.Contains("DYTauTau") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( abs(genlep.ID) == 15 && genlep.isHardProcess )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 taus from hard-process -- //
		{
			GenFlag = kTRUE;
		}
	}
	// -- Madgraph sample -- //
	else if( Tag.Contains("Madgraph") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isMuon() && genlep.isHardProcess )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 muons from hard-process -- //
		{
			if( Tag == "Madgraph_M50to150" ) // -- Select only evetns withtin 50 < M < 150 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 150 )
					GenFlag = kTRUE;
			}
			else if( Tag == "Madgraph_M10to50" )
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M > 10 )
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
	if( Tag.Contains("DYMuMu") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isMuon() && genlep.fromHardProcessFinalState )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 muons from hard-process -- //
		{
			if( Tag == "DYMuMu_M50to200" ) // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 200 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYMuMu_M50to400" ) // -- Select only evetns withtin 50 < M < 400 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 400 )
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
	else if( Tag.Contains("DYEE") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isElectron() && genlep.fromHardProcessFinalState )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 electrons from hard-process -- //
		{
			if( Tag == "DYEE_M50to200" ) // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 200 )
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
	// -- Separate DYTauTau events from MuMu events -- //
	else if( Tag.Contains("DYTauTau") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( abs(genlep.ID) == 15 && genlep.fromHardProcessDecayed )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 taus from hard-process -- //
		{
			GenFlag = kTRUE;
		}
	}
	// -- other cases(e.g. ttbar, WJets, Diboson...): pass
	else
		GenFlag = kTRUE; 

	return GenFlag;
}

void DYAnalyzer::SetupPileUpReWeighting( Bool_t isMC )
{
	if( isMC == kFALSE ) // -- for data -- //
	{
		for(Int_t i=0; i<52; i++)
			PileUpWeight[i] = 1;

		return;
	}
	
	// -- Only for the MC -- //
	TFile *f = new TFile("/home/kplee/CommonCodes/DrellYanAnalysis/ROOTFile_PUReWeight_v20160208_2nd_71mb.root");
	f->cd();
	TH1D *h_weight = (TH1D*)f->Get("h_PUReWeights");
	if( h_weight == NULL )
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
	if( PileUp_MC < 0 || PileUp_MC > 51 )
	{
                std::cout << "[PileUp_MC = " << PileUp_MC << "]: NO CORRESPONDING PU Weight! ... it returns 0" << endl;
		return 0;
	}
	return PileUpWeight[PileUp_MC];
}

void DYAnalyzer::SetupPileUpReWeighting_76X( Bool_t isMC )
{
	if( isMC == kFALSE ) // -- for data -- //
	{
		for(Int_t i=0; i<50; i++)
			PileUpWeight[i] = 1;

		return;
	}
	
	// -- Only for the MC -- //
	TString FileLocation = "/home/kplee/CommonCodes/DrellYanAnalysis/ROOTFile_PUReWeight_76X_v20160404_71mb.root";
	TFile *f = new TFile(FileLocation);
	f->cd();
	TH1D *h_weight = (TH1D*)f->Get("h_PUReWeights");
	if( h_weight == NULL )
	{
                std::cout << "ERROR! ... No Weight histogram!"<< endl;
		return;
	}

	for(Int_t i=0; i<50; i++)
	{
		Int_t i_bin = i+1;
		PileUpWeight[i] = h_weight->GetBinContent(i_bin);
	}
}

Double_t DYAnalyzer::PileUpWeightValue_76X(Int_t PileUp_MC)
{
	if( PileUp_MC < 0 || PileUp_MC > 49 )
	{
                std::cout << "[PileUp_MC = " << PileUp_MC << "]: NO CORRESPONDING PU Weight! ... it returns 0" << endl;
		return 0;
	}
	return PileUpWeight[PileUp_MC];
}

void DYAnalyzer::SetupPileUpReWeighting_80X( Bool_t isMC, TString ROOTFileName )
{
	if( isMC == kFALSE ) // -- for data -- //
	{
		for(Int_t i=0; i<75; i++)
			PileUpWeight[i] = 1;

		return;
	}
	
	// -- Only for the MC -- //
	TString FileLocation = "./etc/PileUp/80X/"+ROOTFileName;
	TFile *f = new TFile(FileLocation);
	TH1D *h_weight = (TH1D*)f->Get("h_PUReWeights");
	if( h_weight == NULL )
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
	if( PileUp_MC < 0 || PileUp_MC > 74 )
	{
                std::cout << "[PileUp_MC = " << PileUp_MC << "]: NO CORRESPONDING PU Weight! ... it returns 0" << endl;
		return 0;
	}
	return PileUpWeight[PileUp_MC];
}

void DYAnalyzer::SetupEfficiencyScaleFactor()
{
	TString Location_TnP = "/home/kplee/CommonCodes/DrellYanAnalysis/ROOTFile_TagProbeEfficiency_v20160329.root";
        std::cout << "[Tag&Probe efficiency is from " << Location_TnP << " (Default, 74X)]" << endl;
	
	TFile *f = new TFile( Location_TnP );
	TH2D *h_RecoID_data = (TH2D*)f->Get("h_2D_Eff_RecoID_Data");
	TH2D *h_RecoID_MC = (TH2D*)f->Get("h_2D_Eff_RecoID_MC");

	TH2D *h_Iso_data = (TH2D*)f->Get("h_2D_Eff_Iso_Data");
	TH2D *h_Iso_MC = (TH2D*)f->Get("h_2D_Eff_Iso_MC");

	TH2D *h_HLTv4p2_data = (TH2D*)f->Get("h_2D_Eff_HLTv4p2_Data");
	TH2D *h_HLTv4p2_MC = (TH2D*)f->Get("h_2D_Eff_HLTv4p2_MC");

	TH2D *h_HLTv4p3_data = (TH2D*)f->Get("h_2D_Eff_HLTv4p3_Data");
	TH2D *h_HLTv4p3_MC = (TH2D*)f->Get("h_2D_Eff_HLTv4p3_MC");


	Int_t nEtaBins = h_RecoID_data->GetNbinsX();
	Int_t nPtBins = h_RecoID_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t RecoID_data = h_RecoID_data->GetBinContent(i_etabin, i_ptbin);
			Double_t RecoID_MC = h_RecoID_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t Iso_data = h_Iso_data->GetBinContent(i_etabin, i_ptbin);
			Double_t Iso_MC = h_Iso_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t HLTv4p2_data = h_HLTv4p2_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLTv4p2_MC = h_HLTv4p2_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t HLTv4p3_data = h_HLTv4p3_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLTv4p3_MC = h_HLTv4p3_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_RecoID_data[iter_x][iter_y] = RecoID_data;
			Eff_RecoID_MC[iter_x][iter_y] = RecoID_MC;

			Eff_Iso_data[iter_x][iter_y] = Iso_data;
			Eff_Iso_MC[iter_x][iter_y] = Iso_MC;

			Eff_HLTv4p2_data[iter_x][iter_y] = HLTv4p2_data;
			Eff_HLTv4p2_MC[iter_x][iter_y] = HLTv4p2_MC;

			Eff_HLTv4p3_data[iter_x][iter_y] = HLTv4p3_data;
			Eff_HLTv4p3_MC[iter_x][iter_y] = HLTv4p3_MC;
		}
	}
        std::cout << "Setting for efficiency correction factors is completed" << endl;
}

void DYAnalyzer::SetupEfficiencyScaleFactor(TString ROOTFileName)
{
	TString Location_TnP = "/home/kplee/CommonCodes/DrellYanAnalysis/"+ROOTFileName;
        std::cout << "[Tag&Probe efficiency is from " << Location_TnP << "]" << endl;

	TFile *f = new TFile( Location_TnP );

	TH2D *h_RecoID_data = (TH2D*)f->Get("h_2D_Eff_RecoID_Data");
	TH2D *h_RecoID_MC = (TH2D*)f->Get("h_2D_Eff_RecoID_MC");

	TH2D *h_Iso_data = (TH2D*)f->Get("h_2D_Eff_Iso_Data");
	TH2D *h_Iso_MC = (TH2D*)f->Get("h_2D_Eff_Iso_MC");

	TH2D *h_HLTv4p2_data = (TH2D*)f->Get("h_2D_Eff_HLTv4p2_Data");
	TH2D *h_HLTv4p2_MC = (TH2D*)f->Get("h_2D_Eff_HLTv4p2_MC");

	TH2D *h_HLTv4p3_data = (TH2D*)f->Get("h_2D_Eff_HLTv4p3_Data");
	TH2D *h_HLTv4p3_MC = (TH2D*)f->Get("h_2D_Eff_HLTv4p3_MC");


	Int_t nEtaBins = h_RecoID_data->GetNbinsX();
	Int_t nPtBins = h_RecoID_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t RecoID_data = h_RecoID_data->GetBinContent(i_etabin, i_ptbin);
			Double_t RecoID_MC = h_RecoID_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t Iso_data = h_Iso_data->GetBinContent(i_etabin, i_ptbin);
			Double_t Iso_MC = h_Iso_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t HLTv4p2_data = h_HLTv4p2_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLTv4p2_MC = h_HLTv4p2_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t HLTv4p3_data = h_HLTv4p3_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLTv4p3_MC = h_HLTv4p3_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_RecoID_data[iter_x][iter_y] = RecoID_data;
			Eff_RecoID_MC[iter_x][iter_y] = RecoID_MC;

			Eff_Iso_data[iter_x][iter_y] = Iso_data;
			Eff_Iso_MC[iter_x][iter_y] = Iso_MC;

			Eff_HLTv4p2_data[iter_x][iter_y] = HLTv4p2_data;
			Eff_HLTv4p2_MC[iter_x][iter_y] = HLTv4p2_MC;

			Eff_HLTv4p3_data[iter_x][iter_y] = HLTv4p3_data;
			Eff_HLTv4p3_MC[iter_x][iter_y] = HLTv4p3_MC;
		}
	}
        std::cout << "Setting for efficiency correction factors is completed" << endl;
}

void DYAnalyzer::SetupEfficiencyScaleFactor_BtoF()
{
	TString Location = "./etc/effSF/effSF_muon/";
        std::cout << "[Tag&Probe efficiency is from " << Location+"*BtoF.root" << "]" << endl;

	TFile *f1 = new TFile( Location+"ID_SF_RunBtoF.root" );
	TH2F *h_RecoID_data = (TH2F*)f1->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA");
	TH2F *h_RecoID_MC = (TH2F*)f1->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesMC/abseta_pair_ne_MC");

	TFile *f2 = new TFile( Location+"ISO_SF_RunBtoF.root" );
	TH2F *h_Iso_data = (TH2F*)f2->Get("tkLooseISO_highptID_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA");
	TH2F *h_Iso_MC = (TH2F*)f2->Get("tkLooseISO_highptID_newpt_eta/efficienciesMC/abseta_pair_ne_MC");

	TFile *f3 = new TFile( Location+"Trigger_SF_RunBtoF.root" );
	TH2F *h_HLT_data = (TH2F*)f3->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
	TH2F *h_HLT_MC = (TH2F*)f3->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/abseta_pt_MC");


	Int_t nEtaBins1 = h_RecoID_data->GetNbinsX();
	Int_t nPtBins1 = h_RecoID_data->GetNbinsY();

	Int_t nEtaBins2 = h_HLT_data->GetNbinsX();
	Int_t nPtBins2 = h_HLT_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins1; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins1; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t RecoID_data = h_RecoID_data->GetBinContent(i_etabin, i_ptbin);
			Double_t RecoID_MC = h_RecoID_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t Iso_data = h_Iso_data->GetBinContent(i_etabin, i_ptbin);
			Double_t Iso_MC = h_Iso_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_RecoID_data_BtoF[iter_x][iter_y] = RecoID_data;
			Eff_RecoID_MC_BtoF[iter_x][iter_y] = RecoID_MC;

			Eff_Iso_data_BtoF[iter_x][iter_y] = Iso_data;
			Eff_Iso_MC_BtoF[iter_x][iter_y] = Iso_MC;
		}
	}
	for(Int_t iter_x = 0; iter_x < nEtaBins2; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins2; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t HLT_data = h_HLT_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLT_MC = h_HLT_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_HLT_data_BtoF[iter_x][iter_y] = HLT_data;
			Eff_HLT_MC_BtoF[iter_x][iter_y] = HLT_MC;

                        //std::cout << "[iter_x, iter_y] : [" << iter_x << ", " << iter_y << "]" << endl;
                        //std::cout << "Eff_HLT_data_1 : " << Eff_HLT_data_BtoF[iter_x][iter_y] << endl;
                        //std::cout << "Eff_HLT_MC_1 : " << Eff_HLT_MC_BtoF[iter_x][iter_y] << endl;
		}
	}
        std::cout << "Setting for efficiency correction factors (BtoF) is completed" << endl;
}

void DYAnalyzer::SetupEfficiencyScaleFactor_GtoH()
{
	TString Location = "./etc/effSF/effSF_muon/";
        std::cout << "[Tag&Probe efficiency is from " << Location+"*GtoH.root" << "]" << endl;

	TFile *f1 = new TFile( Location+"ID_SF_RunGtoH.root" );
	TH2F *h_RecoID_data = (TH2F*)f1->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA");
	TH2F *h_RecoID_MC = (TH2F*)f1->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesMC/abseta_pair_ne_MC");

	TFile *f2 = new TFile( Location+"ISO_SF_RunGtoH.root" );
	TH2F *h_Iso_data = (TH2F*)f2->Get("tkLooseISO_highptID_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA");
	TH2F *h_Iso_MC = (TH2F*)f2->Get("tkLooseISO_highptID_newpt_eta/efficienciesMC/abseta_pair_ne_MC");

	TFile *f3 = new TFile( Location+"Trigger_SF_RunGtoH.root" );
	TH2F *h_HLT_data = (TH2F*)f3->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
	TH2F *h_HLT_MC = (TH2F*)f3->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/abseta_pt_MC");


	Int_t nEtaBins1 = h_RecoID_data->GetNbinsX();
	Int_t nPtBins1 = h_RecoID_data->GetNbinsY();

	Int_t nEtaBins2 = h_HLT_data->GetNbinsX();
	Int_t nPtBins2 = h_HLT_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins1; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins1; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t RecoID_data = h_RecoID_data->GetBinContent(i_etabin, i_ptbin);
			Double_t RecoID_MC = h_RecoID_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t Iso_data = h_Iso_data->GetBinContent(i_etabin, i_ptbin);
			Double_t Iso_MC = h_Iso_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_RecoID_data_GtoH[iter_x][iter_y] = RecoID_data;
			Eff_RecoID_MC_GtoH[iter_x][iter_y] = RecoID_MC;

			Eff_Iso_data_GtoH[iter_x][iter_y] = Iso_data;
			Eff_Iso_MC_GtoH[iter_x][iter_y] = Iso_MC;
		}
	}
	for(Int_t iter_x = 0; iter_x < nEtaBins2; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins2; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t HLT_data = h_HLT_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLT_MC = h_HLT_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_HLT_data_GtoH[iter_x][iter_y] = HLT_data;
			Eff_HLT_MC_GtoH[iter_x][iter_y] = HLT_MC;

                        //std::cout << "[iter_x, iter_y] : [" << iter_x << ", " << iter_y << "]" << endl;
                        //std::cout << "Eff_HLT_data_2 : " << Eff_HLT_data_GtoH[iter_x][iter_y] << endl;
                        //std::cout << "Eff_HLT_MC_2 : " << Eff_HLT_MC_GtoH[iter_x][iter_y] << endl;
		}
	}
        std::cout << "Setting for efficiency correction factors (GtoH) is completed" << endl;
}

void DYAnalyzer::SetupEfficiencyScaleFactor_electron()
{
	TString Location = "./etc/effSF/effSF_electron/";
        std::cout << "[Tag&Probe efficiency is from " << Location+"*.root" << "]" << endl;

	TFile *f1 = new TFile( Location+"Reco_SF.root" );
	TGraphErrors *h_reco_sf = (TGraphErrors*)f1->Get("grSF1D_0");

	TFile *f2 = new TFile( Location+"MediumID_SF.root" );
	TGraphErrors *h_id_sf_0 = (TGraphErrors*)f2->Get("grSF1D_0");
	TGraphErrors *h_id_sf_1 = (TGraphErrors*)f2->Get("grSF1D_1");
	TGraphErrors *h_id_sf_2 = (TGraphErrors*)f2->Get("grSF1D_2");
	TGraphErrors *h_id_sf_3 = (TGraphErrors*)f2->Get("grSF1D_3");
	TGraphErrors *h_id_sf_4 = (TGraphErrors*)f2->Get("grSF1D_4");

	Int_t nEtaBins_reco = h_reco_sf->GetN();
	Int_t nPtBins_reco = 1;

	Int_t nEtaBins_id = h_id_sf_0->GetN();
	Int_t nPtBins_id = 5;

	// -- For reco efficiency -- //
	Double_t x_reco[30]; Double_t xx_reco[30];
	Double_t y_reco[30]; Double_t yy_reco[30];

	for(Int_t iter_y = 0; iter_y < nPtBins_reco; iter_y++)
	{
		for(Int_t i=0; i<nEtaBins_reco; i++)
		{
			h_reco_sf->GetPoint(i, x_reco[i], y_reco[i]);
		}

		// -- Rearrangement in order of x values -- //
		Double_t xmin;
		Double_t xlow = -9999;
		for(Int_t j=0; j<nEtaBins_reco; j++)
		{
			Int_t jj = -9999;

			xmin = 9999;
			for(Int_t k=0; k<nEtaBins_reco; k++)
			{
				if( xlow < x_reco[k] && x_reco[k] < xmin )
				{
					jj = k;
					xmin = x_reco[k];
				}
			}
			xx_reco[j] = x_reco[jj];
			yy_reco[j] = y_reco[jj];

			xlow = xmin;
//			std::cout << j << "  " << xx_reco[j] << "  " << yy_reco[j] << endl;
		} // End of rearrangement

		for(Int_t iter_x = 0; iter_x < nEtaBins_reco; iter_x++)
		{
			Eff_Reco_data[iter_x][iter_y] = yy_reco[iter_x]; // actually, it is the scale factor.
			Eff_Reco_MC[iter_x][iter_y] = 1;
//			std::cout << "Reco: bin# = [" << iter_x << ", " << iter_y << "]" << " eta = " << xx_reco[iter_x] << " sf = " << yy_reco[iter_x] << endl;
		}
	}

	// -- For id efficiency -- //
	Double_t x_id[10]; Double_t xx_id[10];
	Double_t y_id[10]; Double_t yy_id[10];

	TGraphErrors *h_id_sf;
	for(Int_t iter_y = 0; iter_y < nPtBins_id; iter_y++)
	{
		if(iter_y == 0) h_id_sf = (TGraphErrors*)h_id_sf_0->Clone();
		else if(iter_y == 1) h_id_sf = (TGraphErrors*)h_id_sf_1->Clone();
		else if(iter_y == 2) h_id_sf = (TGraphErrors*)h_id_sf_2->Clone();
		else if(iter_y == 3) h_id_sf = (TGraphErrors*)h_id_sf_3->Clone();
		else if(iter_y == 4) h_id_sf = (TGraphErrors*)h_id_sf_4->Clone();

		for(Int_t i=0; i<nEtaBins_id; i++)
		{
			h_id_sf->GetPoint(i, x_id[i], y_id[i]);
		}

		// -- Rearrangement in order of x values -- //
		Double_t xmin;
		Double_t xlow = -9999;
		for(Int_t j=0; j<nEtaBins_id; j++)
		{
			Int_t jj = -9999;

			xmin = 9999;
			for(Int_t k=0; k<nEtaBins_id; k++)
			{
				if( xlow < x_id[k] && x_id[k] < xmin )
				{
					jj = k;
					xmin = x_id[k];
				}
			}
			xx_id[j] = x_id[jj];
			yy_id[j] = y_id[jj];

			xlow = xmin;
//			std::cout << j << "  " << xx_id[j] << "  " << yy_id[j] << endl;
		} // End of rearrangement

		for(Int_t iter_x = 0; iter_x < nEtaBins_id; iter_x++)
		{
			Eff_ID_data[iter_x][iter_y] = yy_id[iter_x]; // actually, it is the scale factor.
			Eff_ID_MC[iter_x][iter_y] = 1;
//			std::cout << "ID: bin# = [" << iter_x << ", " << iter_y << "]" << " eta = " << xx_id[iter_x] << " sf = " << yy_id[iter_x] << endl;
		}
	}

        std::cout << "Setting for efficiency correction factors is completed" << endl;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLTv4p2(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	if( ptbin1 == 9999 || etabin1 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin1, etabin1);
		return -999;
	}

	Double_t Eff_muon1_data = Eff_RecoID_data[etabin1][ptbin1] * Eff_Iso_data[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC[etabin1][ptbin1] * Eff_Iso_MC[etabin1][ptbin1];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2 = FindPtBin( Pt2 );
	Int_t etabin2 = FindEtaBin( eta2 );

	if( ptbin2 == 9999 || etabin2 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin2, etabin2);
		return -999;
	}
	Double_t Eff_muon2_data = Eff_RecoID_data[etabin2][ptbin2] * Eff_Iso_data[etabin2][ptbin2];
	Double_t Eff_muon2_MC = Eff_RecoID_MC[etabin2][ptbin2] * Eff_Iso_MC[etabin2][ptbin2];


	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLTv4p2_data[etabin1][ptbin1];
	Double_t Eff_Trig_muon2_data = Eff_HLTv4p2_data[etabin2][ptbin2];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLTv4p2_MC[etabin1][ptbin1];
	Double_t Eff_Trig_muon2_MC = Eff_HLTv4p2_MC[etabin2][ptbin2];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

        // std::cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("[Data]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin1][ptbin1], Eff_Iso_data[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin2][ptbin2], Eff_Iso_data[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin1][ptbin1], Eff_Iso_MC[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin2][ptbin2], Eff_Iso_MC[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("(ptbin1, etabin1, ptbin2, etabin2): (%d, %d, %d, %d)\n", ptbin1, etabin1, ptbin2, etabin2);
		
		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLTv4p3(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	if( ptbin1 == 9999 || etabin1 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin1, etabin1);
		return -999;
	}

	Double_t Eff_muon1_data = Eff_RecoID_data[etabin1][ptbin1] * Eff_Iso_data[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC[etabin1][ptbin1] * Eff_Iso_MC[etabin1][ptbin1];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2 = FindPtBin( Pt2 );
	Int_t etabin2 = FindEtaBin( eta2 );

	if( ptbin2 == 9999 || etabin2 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin2, etabin2);
		return -999;
	}
	Double_t Eff_muon2_data = Eff_RecoID_data[etabin2][ptbin2] * Eff_Iso_data[etabin2][ptbin2];
	Double_t Eff_muon2_MC = Eff_RecoID_MC[etabin2][ptbin2] * Eff_Iso_MC[etabin2][ptbin2];


	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLTv4p3_data[etabin1][ptbin1];
	Double_t Eff_Trig_muon2_data = Eff_HLTv4p3_data[etabin2][ptbin2];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLTv4p3_MC[etabin1][ptbin1];
	Double_t Eff_Trig_muon2_MC = Eff_HLTv4p3_MC[etabin2][ptbin2];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

        // std::cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("[Data]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin1][ptbin1], Eff_Iso_data[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin2][ptbin2], Eff_Iso_data[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin1][ptbin1], Eff_Iso_MC[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin2][ptbin2], Eff_Iso_MC[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("(ptbin1, etabin1, ptbin2, etabin2): (%d, %d, %d, %d)\n", ptbin1, etabin1, ptbin2, etabin2);
		
		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_BtoF(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	Int_t ptbin1_trig = FindPtBin_trig( Pt1 );
	Int_t etabin1_trig = FindEtaBin_trig( eta1 );

	if( ptbin1 == 9999 || etabin1 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin1, etabin1);
		return -999;
	}

	Double_t Eff_muon1_data = Eff_RecoID_data_BtoF[etabin1][ptbin1] * Eff_Iso_data_BtoF[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC_BtoF[etabin1][ptbin1] * Eff_Iso_MC_BtoF[etabin1][ptbin1];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2 = FindPtBin( Pt2 );
	Int_t etabin2 = FindEtaBin( eta2 );

	Int_t ptbin2_trig = FindPtBin_trig( Pt2 );
	Int_t etabin2_trig = FindEtaBin_trig( eta2 );

	if( ptbin2 == 9999 || etabin2 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin2, etabin2);
		return -999;
	}
	Double_t Eff_muon2_data = Eff_RecoID_data_BtoF[etabin2][ptbin2] * Eff_Iso_data_BtoF[etabin2][ptbin2];
	Double_t Eff_muon2_MC = Eff_RecoID_MC_BtoF[etabin2][ptbin2] * Eff_Iso_MC_BtoF[etabin2][ptbin2];


	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_trig][ptbin1_trig];
	Double_t Eff_Trig_muon2_data = Eff_HLT_data_BtoF[etabin2_trig][ptbin2_trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_trig][ptbin1_trig];
	Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_BtoF[etabin2_trig][ptbin2_trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

        //std::cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("[Data]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data_BtoF[etabin1][ptbin1], Eff_Iso_data_BtoF[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data_BtoF[etabin2][ptbin2], Eff_Iso_data_BtoF[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC_BtoF[etabin1][ptbin1], Eff_Iso_MC_BtoF[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC_BtoF[etabin2][ptbin2], Eff_Iso_MC_BtoF[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("(ptbin1, etabin1, ptbin2, etabin2): (%d, %d, %d, %d)\n", ptbin1, etabin1, ptbin2, etabin2);
		
		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_GtoH(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	Int_t ptbin1_trig = FindPtBin_trig( Pt1 );
	Int_t etabin1_trig = FindEtaBin_trig( eta1 );

	if( ptbin1 == 9999 || etabin1 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin1, etabin1);
		return -999;
	}

	Double_t Eff_muon1_data = Eff_RecoID_data_GtoH[etabin1][ptbin1] * Eff_Iso_data_GtoH[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC_GtoH[etabin1][ptbin1] * Eff_Iso_MC_GtoH[etabin1][ptbin1];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2 = FindPtBin( Pt2 );
	Int_t etabin2 = FindEtaBin( eta2 );

	Int_t ptbin2_trig = FindPtBin_trig( Pt2 );
	Int_t etabin2_trig = FindEtaBin_trig( eta2 );

	if( ptbin2 == 9999 || etabin2 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin2, etabin2);
		return -999;
	}
	Double_t Eff_muon2_data = Eff_RecoID_data_GtoH[etabin2][ptbin2] * Eff_Iso_data_GtoH[etabin2][ptbin2];
	Double_t Eff_muon2_MC = Eff_RecoID_MC_GtoH[etabin2][ptbin2] * Eff_Iso_MC_GtoH[etabin2][ptbin2];


	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_trig][ptbin1_trig];
	Double_t Eff_Trig_muon2_data = Eff_HLT_data_GtoH[etabin2_trig][ptbin2_trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_trig][ptbin1_trig];
	Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_GtoH[etabin2_trig][ptbin2_trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

        //std::cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("[Data]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data_GtoH[etabin1][ptbin1], Eff_Iso_data_GtoH[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data_GtoH[etabin2][ptbin2], Eff_Iso_data_GtoH[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC_GtoH[etabin1][ptbin1], Eff_Iso_MC_GtoH[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC_GtoH[etabin2][ptbin2], Eff_Iso_MC_GtoH[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("(ptbin1, etabin1, ptbin2, etabin2): (%d, %d, %d, %d)\n", ptbin1, etabin1, ptbin2, etabin2);
		
		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_electron(Electron ele1, Electron ele2)
{
	Double_t weight = -999;

	// -- Electron1 -- //
	Double_t Pt1 = ele1.Pt;
	//Double_t eta1 = ele1.eta;
	Double_t eta1 = ele1.etaSC;

	Int_t ptbin1_Reco = FindPtBin_Reco( Pt1 );
	Int_t etabin1_Reco = FindEtaBin_Reco( eta1 );

	Int_t ptbin1_ID = FindPtBin_ID( Pt1 );
	Int_t etabin1_ID = FindEtaBin_ID( eta1 );

	Double_t Eff_ele1_data = Eff_Reco_data[etabin1_Reco][ptbin1_Reco] * Eff_ID_data[etabin1_ID][ptbin1_ID];
	Double_t Eff_ele1_MC = 1;

	// -- Electron2 -- //
	Double_t Pt2 = ele2.Pt;
	//Double_t eta2 = ele2.eta;
	Double_t eta2 = ele2.etaSC;

	Int_t ptbin2_Reco = FindPtBin_Reco( Pt2 );
	Int_t etabin2_Reco = FindEtaBin_Reco( eta2 );

	Int_t ptbin2_ID = FindPtBin_ID( Pt2 );
	Int_t etabin2_ID = FindEtaBin_ID( eta2 );

	Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
	Double_t Eff_ele2_MC = 1;

	// -- This is trigger part -- // trigger SF is not yet.
/*	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_trig][ptbin1_trig];
	Double_t Eff_Trig_muon2_data = Eff_HLT_data_GtoH[etabin2_trig][ptbin2_trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_trig][ptbin1_trig];
	Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_GtoH[etabin2_trig][ptbin2_trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;*/

//	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
//	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;
	Double_t Eff_data_all = Eff_ele1_data * Eff_ele2_data;
	Double_t Eff_MC_all = Eff_ele1_MC * Eff_ele2_MC;

	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 ) printf("[SF] Weight = %.3lf\n", weight);

	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_emu_BtoF(Muon mu, Electron ele)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu.Pt;
	Double_t eta1 = mu.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	Int_t ptbin1_trig = FindPtBin_trig( Pt1 );
	Int_t etabin1_trig = FindEtaBin_trig( eta1 );

	// -- Electron2 -- //
	Double_t Pt2 = ele.Pt;
	//Double_t eta2 = ele.eta;
	Double_t eta2 = ele.etaSC;

	Int_t ptbin2_Reco = FindPtBin_Reco( Pt2 );
	Int_t etabin2_Reco = FindEtaBin_Reco( eta2 );

	Int_t ptbin2_ID = FindPtBin_ID( Pt2 );
	Int_t etabin2_ID = FindEtaBin_ID( eta2 );

	//Check about bin settings
	if( ptbin1 == 9999 || etabin1 == 9999 || ptbin1_trig == 9999 || etabin1_trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
	{
		printf("ERROR! Wrong assigned bin number ...");
		return -999;
	}

	//Muon1
	Double_t Eff_muon1_data = Eff_RecoID_data_BtoF[etabin1][ptbin1] * Eff_Iso_data_BtoF[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC_BtoF[etabin1][ptbin1] * Eff_Iso_MC_BtoF[etabin1][ptbin1];

	//Electron2
	Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
	Double_t Eff_ele2_MC = 1;

	//Trigger : We consider only SingleMuon trigger
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_trig][ptbin1_trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_trig][ptbin1_trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC;

	// -- emu event SF -- //
	Double_t Eff_data_all = Eff_muon1_data * Eff_ele2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 ) printf("[SF] Weight = %.3lf\n", weight);
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_emu_GtoH(Muon mu, Electron ele)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu.Pt;
	Double_t eta1 = mu.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	Int_t ptbin1_trig = FindPtBin_trig( Pt1 );
	Int_t etabin1_trig = FindEtaBin_trig( eta1 );

	// -- Electron2 -- //
	Double_t Pt2 = ele.Pt;
	//Double_t eta2 = ele.eta;
	Double_t eta2 = ele.etaSC;

	Int_t ptbin2_Reco = FindPtBin_Reco( Pt2 );
	Int_t etabin2_Reco = FindEtaBin_Reco( eta2 );

	Int_t ptbin2_ID = FindPtBin_ID( Pt2 );
	Int_t etabin2_ID = FindEtaBin_ID( eta2 );

	//Check about bin settings
	if( ptbin1 == 9999 || etabin1 == 9999 || ptbin1_trig == 9999 || etabin1_trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
	{
		printf("ERROR! Wrong assigned bin number ...");
		return -999;
	}

	//Muon1
	Double_t Eff_muon1_data = Eff_RecoID_data_GtoH[etabin1][ptbin1] * Eff_Iso_data_GtoH[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC_GtoH[etabin1][ptbin1] * Eff_Iso_MC_GtoH[etabin1][ptbin1];

	//Electron2
	Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
	Double_t Eff_ele2_MC = 1;

	//Trigger : We consider only SingleMuon trigger
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_trig][ptbin1_trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_trig][ptbin1_trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC;

	// -- emu event SF -- //
	Double_t Eff_data_all = Eff_muon1_data * Eff_ele2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 ) printf("[SF] Weight = %.3lf\n", weight);
	return weight;
}

Int_t DYAnalyzer::FindPtBin(Double_t Pt)
{
	//const Int_t nPtBins = 4;
	const Int_t nPtBins = 7;
	//Double_t PtBinEdges[nPtBins+1] = {10, 22, 40, 70, 250};
	Double_t PtBinEdges[nPtBins+1] = {20, 25, 30, 40, 50, 55, 60, 120};

	Int_t ptbin = 9999;

	// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
	if( Pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	// -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
	else if( Pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::FindPtBin_Reco(Double_t Pt)
{
	//const Int_t nPtBins = 4;
	const Int_t nPtBins = 1;
	//Double_t PtBinEdges[nPtBins+1] = {10, 22, 40, 70, 250};
	Double_t PtBinEdges[nPtBins+1] = {25, 500};

	Int_t ptbin = 9999;

	// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
	if( Pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	// -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
	else if( Pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::FindPtBin_ID(Double_t Pt)
{
	//const Int_t nPtBins = 4;
	const Int_t nPtBins = 5;
	//Double_t PtBinEdges[nPtBins+1] = {10, 22, 40, 70, 250};
	Double_t PtBinEdges[nPtBins+1] = {10, 20, 35, 50, 90, 150};

	Int_t ptbin = 9999;

	// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
	if( Pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	// -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
	else if( Pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::FindPtBin_trig(Double_t Pt)
{
	//const Int_t nPtBins = 4;
	const Int_t nPtBins = 7;
	//Double_t PtBinEdges[nPtBins+1] = {10, 22, 40, 70, 250};
	Double_t PtBinEdges[nPtBins+1] = {26, 30, 40, 50, 60, 120, 200, 500};

	Int_t ptbin = 9999;

	// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
	if( Pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	// -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
	else if( Pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::FindEtaBin(Double_t eta)
{
	//const Int_t nEtaBins = 5;
	const Int_t nEtaBins = 4;
	//Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
	Double_t EtaBinEdges[nEtaBins+1] = {0, 0.9, 1.2, 2.1, 2.4};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		//if( eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
		if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}

Int_t DYAnalyzer::FindEtaBin_Reco(Double_t eta)
{
	//const Int_t nEtaBins = 5;
	const Int_t nEtaBins = 30;
	//Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
	Double_t EtaBinEdges[nEtaBins+1] = {-2.5,-2.45,-2.4,-2.3,-2.2,-2.0,-1.8,-1.63,-1.566,-1.4442,-1.2,-1.0,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,1.0,1.2,1.4442,1.566,1.63,1.8,2.0,2.2,2.3,2.4,2.45,2.5};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		if( eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
		//if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}

Int_t DYAnalyzer::FindEtaBin_ID(Double_t eta)
{
	//const Int_t nEtaBins = 5;
	const Int_t nEtaBins = 10;
	//Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
	Double_t EtaBinEdges[nEtaBins+1] = {-2.5, -2.0, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.0, 2.5};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		if( eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
		//if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}

Int_t DYAnalyzer::FindEtaBin_trig(Double_t eta)
{
	//const Int_t nEtaBins = 5;
	const Int_t nEtaBins = 4;
	//Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
	Double_t EtaBinEdges[nEtaBins+1] = {0, 0.9, 1.2, 2.1, 2.4};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		//if( eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
		if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight(Muon mu1, Muon mu2, NtupleHandle *ntuple)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	if( ptbin1 == 9999 || etabin1 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin1, etabin1);
		return -999;
	}

	Double_t Eff_muon1_data = Eff_RecoID_data[etabin1][ptbin1] * Eff_Iso_data[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC[etabin1][ptbin1] * Eff_Iso_MC[etabin1][ptbin1];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2 = FindPtBin( Pt2 );
	Int_t etabin2 = FindEtaBin( eta2 );

	if( ptbin2 == 9999 || etabin2 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin2, etabin2);
		return -999;
	}
	Double_t Eff_muon2_data = Eff_RecoID_data[etabin2][ptbin2] * Eff_Iso_data[etabin2][ptbin2];
	Double_t Eff_muon2_MC = Eff_RecoID_MC[etabin2][ptbin2] * Eff_Iso_MC[etabin2][ptbin2];

	Bool_t isHLTv4p2 = kFALSE;
	if( ntuple->runNum < 257932.5 )
		isHLTv4p2 = kTRUE;

	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;
	if( isHLTv4p2 )
	{
		Double_t Eff_Trig_muon1_data = Eff_HLTv4p2_data[etabin1][ptbin1];
		Double_t Eff_Trig_muon2_data = Eff_HLTv4p2_data[etabin2][ptbin2];
		Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

		Double_t Eff_Trig_muon1_MC = Eff_HLTv4p2_MC[etabin1][ptbin1];
		Double_t Eff_Trig_muon2_MC = Eff_HLTv4p2_MC[etabin2][ptbin2];
		Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;
	}
	else
	{
		Double_t Eff_Trig_muon1_data = Eff_HLTv4p3_data[etabin1][ptbin1];
		Double_t Eff_Trig_muon2_data = Eff_HLTv4p3_data[etabin2][ptbin2];
		Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

		Double_t Eff_Trig_muon1_MC = Eff_HLTv4p3_MC[etabin1][ptbin1];
		Double_t Eff_Trig_muon2_MC = Eff_HLTv4p3_MC[etabin2][ptbin2];
		Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;
	}

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

        // std::cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("[Data]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin1][ptbin1], Eff_Iso_data[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin2][ptbin2], Eff_Iso_data[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin1][ptbin1], Eff_Iso_MC[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin2][ptbin2], Eff_Iso_MC[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("(ptbin1, etabin1, ptbin2, etabin2): (%d, %d, %d, %d)\n", ptbin1, etabin1, ptbin2, etabin2);
		
		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_RecoIdIso(Muon mu1, Muon mu2, NtupleHandle *ntuple)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	if( ptbin1 == 9999 || etabin1 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin1, etabin1);
		return -999;
	}

	Double_t Eff_muon1_data = Eff_RecoID_data[etabin1][ptbin1] * Eff_Iso_data[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC[etabin1][ptbin1] * Eff_Iso_MC[etabin1][ptbin1];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2 = FindPtBin( Pt2 );
	Int_t etabin2 = FindEtaBin( eta2 );

	if( ptbin2 == 9999 || etabin2 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin2, etabin2);
		return -999;
	}
	Double_t Eff_muon2_data = Eff_RecoID_data[etabin2][ptbin2] * Eff_Iso_data[etabin2][ptbin2];
	Double_t Eff_muon2_MC = Eff_RecoID_MC[etabin2][ptbin2] * Eff_Iso_MC[etabin2][ptbin2];

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC;

	weight = Eff_data_all / Eff_MC_all;

	return weight;
}

Bool_t DYAnalyzer::EventSelection(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    if( MuonCollection[j].isTightMuon() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
		if( mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

	if( isExistHLTMatchedMuon == kTRUE )
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
		if( nQMuons == 2)
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
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

			Bool_t isOS = kFALSE;
			if( recolep1.charge != recolep2.charge ) isOS = kTRUE;

			// if( reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005 )
//			if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
			if( reco_M > 60 && reco_M <120 && isPassAcc == kTRUE && isOS == kTRUE )
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}
		}
		else if( nQMuons > 2 )
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
				if( Mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

						if( j_mu != i_mu ) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

							if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
								if( VtxNormChi2_temp < VtxNormChi2_BestPair )
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

			if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

				Bool_t isOS = kFALSE;
				if( mu1_BestPair.charge != mu2_BestPair.charge ) isOS = kTRUE;

//				if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
				if( reco_M > 60 && reco_M < 120 && isOS == kTRUE )
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- End of else if( nQMuons > 2 ) -- //

	} // -- End of if( isExistHLTMatchedMuon == kTRUE ) -- //

	return isPassEventSelection;
}

// -- Test using the trigger without isolation condition: HLT_Mu50_v* -- //
Bool_t DYAnalyzer::EventSelection_Mu50(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    if( MuonCollection[j].isHighPtMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
		if( mu.isTrigMatched(ntuple, HLT) )
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

	if( isExistHLTMatchedMuon == kTRUE )
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
		if( nQMuons == 2)
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
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

			// if( reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005 )
			if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 )
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}
		}
		else if( nQMuons > 2 )
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
				if( Mu.isTrigMatched(ntuple, HLT) )
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

						if( j_mu != i_mu ) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

							if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
								if( VtxNormChi2_temp < VtxNormChi2_BestPair )
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

			if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

				if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 )
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- End of else if( nQMuons > 2 ) -- //

	} // -- End of if( isExistHLTMatchedMuon == kTRUE ) -- //

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
	    if( MuonCollection[j].isHighPtMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
		if( mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

	if( isExistHLTMatchedMuon == kTRUE )
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
		if( nQMuons == 2)
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
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

			// if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 )
			if( reco_M > 10 && isPassAcc == kTRUE && Angle < TMath::Pi() - 0.005 )
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}
		}
		else if( nQMuons > 2 )
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
				if( Mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

						if( j_mu != i_mu ) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

							if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
								if( VtxNormChi2_temp < VtxNormChi2_BestPair )
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

			if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

				// if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 )
				if( reco_M > 10 && Angle < TMath::Pi() - 0.005 )
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- End of else if( nQMuons > 2 ) -- //

	} // -- End of if( isExistHLTMatchedMuon == kTRUE ) -- //

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
	    if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
	    //if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].RelPFIso_dBeta < 0.10)
            {
                QMuonCollection.push_back( MuonCollection[j] );
                QIndex.push_back( j );
            }
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();
        if( nQMuons == 2 && nQIndices == 2 )
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);
                for( UInt_t i=0; i<ntuple->vtxTrkProb->size(); i++ )
                {
                    if( VtxProb == ntuple->vtxTrkProb->at(i) )
                    {
                        DiIndex = i;
                    }
                }

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
                        Index->push_back( QIndex[0] );
                        Index->push_back( QIndex[1] );
                        IndexDi = DiIndex;
		}
	}
	else if( nQMuons > 2 )
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

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

                for( UInt_t i=0; i<ntuple->vtxTrkProb->size(); i++ )
                {
                    if( VtxProb == ntuple->vtxTrkProb->at(i) )
                    {
                        DiIndex = i;
                    }
                }

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
                        Index->push_back( Index_leading );
                        Index->push_back( Index_sub );
                        IndexDi = DiIndex;
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
            if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
            //if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].RelPFIso_dBeta < 0.10)
            {
                QMuonCollection.push_back( MuonCollection[j] );
                QIndex.push_back( j );
            }
        }

        Int_t nQMuons = (Int_t)QMuonCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();
        if( nQMuons == 2 && nQIndices == 2 )
        {
                Muon recolep1 = QMuonCollection[0];
                Muon recolep2 = QMuonCollection[1];

                // -- Check the Accpetance -- //
                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

                // -- Opposite sign condition -- //
                Bool_t isOppositeSign = kFALSE;
                if( recolep1.charge != recolep2.charge )
                        isOppositeSign = kTRUE;

                Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

                TLorentzVector inner_v1 = recolep1.Momentum_Inner;
                TLorentzVector inner_v2 = recolep2.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
                if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back( recolep1 );
                        SelectedMuonCollection->push_back( recolep2 );
                        Index->push_back( QIndex[0] );
                        Index->push_back( QIndex[1] );
                }
        }
        else if( nQMuons > 2 )
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

                        if( Mu.Pt > Pt_leading )
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
                        if( i_mu2 == i_leading ) continue;

                        Muon Mu = QMuonCollection[i_mu2];

                        if( Mu.Pt > Pt_sub )
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
                if( LeadingMuon.charge != SubMuon.charge )
                        isOppositeSign = kTRUE;

                Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

                TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
                TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
                if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back( LeadingMuon );
                        SelectedMuonCollection->push_back( SubMuon );
                        Index->push_back( Index_leading );
                        Index->push_back( Index_sub );
                }

        } // -- End of else if( nQMuons > 2 ) -- //

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
            if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
            //if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].RelPFIso_dBeta < 0.10)
            {
                QMuonCollection.push_back( MuonCollection[j] );
                QIndex.push_back( j );
            }
        }

        Int_t nQMuons = (Int_t)QMuonCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();
        if( nQMuons == 2 && nQIndices == 2 )
        {
                Muon recolep1 = QMuonCollection[0];
                Muon recolep2 = QMuonCollection[1];

                // -- Check the Accpetance -- //
                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

                // -- Opposite sign condition -- //
                Bool_t isOppositeSign = kFALSE;
                if( recolep1.charge != recolep2.charge )
                        isOppositeSign = kTRUE;

                Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);
                for( UInt_t i=0; i<ntuple->vtxTrkProb->size(); i++ )
                {
                    if( VtxProb == ntuple->vtxTrkProb->at(i) )
                    {
                        DiIndex = i;
                    }
                }

                TLorentzVector inner_v1 = recolep1.Momentum_Inner;
                TLorentzVector inner_v2 = recolep2.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
                if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back( recolep1 );
                        SelectedMuonCollection->push_back( recolep2 );
                        Index->push_back( QIndex[0] );
                        Index->push_back( QIndex[1] );
                        IndexDi = DiIndex;
                }
        }
        else if( nQMuons > 2 )
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

                        if( Mu.Pt > Pt_leading )
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
                        if( i_mu2 == i_leading ) continue;

                        Muon Mu = QMuonCollection[i_mu2];

                        if( Mu.Pt > Pt_sub )
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
                if( LeadingMuon.charge != SubMuon.charge )
                        isOppositeSign = kTRUE;

                Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

                for( UInt_t i=0; i<ntuple->vtxTrkProb->size(); i++ )
                {
                    if( VtxProb == ntuple->vtxTrkProb->at(i) )
                    {
                        DiIndex = i;
                    }
                }

                TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
                TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
                if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back( LeadingMuon );
                        SelectedMuonCollection->push_back( SubMuon );
                        Index->push_back( Index_leading );
                        Index->push_back( Index_sub );
                        IndexDi = DiIndex;
                }

        } // -- End of else if( nQMuons > 2 ) -- //

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
            if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
            //if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].RelPFIso_dBeta < 0.10)
            {
                QMuonCollection.push_back( MuonCollection[j] );
                QIndex.push_back( j );
            }
        }

        Int_t nQMuons = (Int_t)QMuonCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();
        if( nQMuons == 2 && nQIndices == 2 )
        {
                Muon recolep1 = QMuonCollection[0];
                Muon recolep2 = QMuonCollection[1];

                // -- Check the Accpetance -- //
                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

                // -- Opposite sign condition -- //
                Bool_t isOppositeSign = kFALSE;
                if( recolep1.charge != recolep2.charge )
                        isOppositeSign = kTRUE;

                Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

                TLorentzVector inner_v1 = recolep1.Momentum_Inner;
                TLorentzVector inner_v2 = recolep2.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
                if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back( recolep1 );
                        SelectedMuonCollection->push_back( recolep2 );
                        Index->push_back( QIndex[0] );
                        Index->push_back( QIndex[1] );
                }
        }
        else if( nQMuons > 2 )
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

                        if( Mu.Pt > Pt_leading )
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
                        if( i_mu2 == i_leading ) continue;

                        Muon Mu = QMuonCollection[i_mu2];

                        if( Mu.Pt > Pt_sub )
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
                if( LeadingMuon.charge != SubMuon.charge )
                        isOppositeSign = kTRUE;

                Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

                // -- Dimuon Vtx cut -- //
                Double_t VtxProb = -999;
                Double_t VtxNormChi2 = 999;
                DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

                TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
                TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

                // -- 3D open angle -- //
                Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
                if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back( LeadingMuon );
                        SelectedMuonCollection->push_back( SubMuon );
                        Index->push_back( Index_leading );
                        Index->push_back( Index_sub );
                }

        } // -- End of else if( nQMuons > 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon_minus_isGLB() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon_minus_muonHits() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon_minus_nMatches() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon_minus_dpT_over_pT() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon_minus_dxyVTX() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon_minus_pixelHits() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon_minus_trackerLayers() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

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
	    //if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isHighPtMuon() )
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
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
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
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
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

// test for emu event selection not using vtx cut
// Updated to use synchronized acceptance cut : 19 Jan. 2018
Bool_t DYAnalyzer::EventSelection_emu_method_test(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple,
						vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection)
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
		if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10
			//&& MuonCollection[j].Pt > LeadPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut )
			&& MuonCollection[j].Pt > SubPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut ) // pT>17 && |eta|<2.4
			QMuonCollection.push_back( MuonCollection[j] );
	}

	//Collect qualified electrons among electrons
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.passMediumID == kTRUE
			//&& elec.Pt > SubPtCut && fabs(elec.etaSC) < 2.5 && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < LeadEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) ) // pT>17 && |eta|<2.4
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	Double_t Pt_mu = 0;
	Double_t Pt_el = 0;
	Muon mu_BestPair;
	Electron el_BestPair;

	// -- Select muon with highest pT -- //
	for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
	{
		Muon Mu = QMuonCollection[i_mu];

		// -- muon should be matched with HLT objects in emu best pair -- //
		if( Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") )
		{
			if( Mu.Pt > Pt_mu )
			{
				Pt_mu		= Mu.Pt;
				mu_BestPair	= Mu;
			}
		}
	}

	// -- Select electron with highest pT -- //
	for(Int_t j_el=0; j_el<nQElectrons; j_el++)
	{
		Electron El = QElectronCollection[j_el];

		if( El.Pt > Pt_el )
		{
			Pt_el		= El.Pt;
			el_BestPair	= El;
		}
	}

	//if( Pt_mu > 0 && Pt_el > 0 )
	if( Pt_mu > 0 && Pt_el > 0 && ( Pt_mu > LeadPtCut || Pt_el > LeadPtCut ) ) // At least one lepton has pT above 28 [GeV]
	{
		TLorentzVector reco_v1 = mu_BestPair.Momentum;
		TLorentzVector reco_v2 = el_BestPair.Momentum;
		Double_t reco_M = (reco_v1 + reco_v2).M();

		// -- 3D open angle -- //
		Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

		//if( reco_M > 10 && Angle < TMath::Pi() - 0.005 )
		if( reco_M > 10 )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( mu_BestPair );
			SelectedElectronCollection->push_back( el_BestPair );
		}
	}

	return isPassEventSelection;
}

// test for emu event selection not using vtx cut
// Derived by Marijus Ambrozas 2018.08.07 to return indices of particles that passed the selection
Bool_t DYAnalyzer::EventSelection_emu_method_test(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple, // Input: electron, muon vectors and NtupleHandle
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
                if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10
                        //&& MuonCollection[j].Pt > LeadPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut )
                        && MuonCollection[j].Pt > SubPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut ) // pT>17 && |eta|<2.4
                {
                        QMuonCollection.push_back( MuonCollection[j] );
                        QIndexMu.push_back( j );
                }
        }

        //Collect qualified electrons among electrons
        vector< Electron > QElectronCollection;
        vector< Int_t > QIndexEle;
        for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
        {
                Electron elec = ElectronCollection[j];
                if( elec.passMediumID == kTRUE
                        //&& elec.Pt > SubPtCut && fabs(elec.etaSC) < 2.5 && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < LeadEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) ) // pT>17 && |eta|<2.4
                {
                        QElectronCollection.push_back( ElectronCollection[j] );
                        QIndexEle.push_back( j );
                }
        }

        Int_t nQMuons = (Int_t)QMuonCollection.size();
        Int_t nQElectrons = (Int_t)QElectronCollection.size();
        Int_t nQIndicesMu = (Int_t)QIndexMu.size();
        Int_t nQIndicesEle = (Int_t)QIndexEle.size();

        Double_t Pt_mu = 0;
        Double_t Pt_el = 0;
        Int_t Index_mu = -1;
        Int_t Index_ele = -1;
        Muon mu_BestPair;
        Electron el_BestPair;

        // -- Select muon with highest pT -- //
        if ( nQMuons == nQIndicesMu )
        {
            for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
            {
                    Muon Mu = QMuonCollection[i_mu];

                    // -- muon should be matched with HLT objects in emu best pair -- //
                    if( Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") )
                    {
                            if( Mu.Pt > Pt_mu )
                            {
                                    Pt_mu           = Mu.Pt;
                                    mu_BestPair     = Mu;
                                    Index_mu        = QIndexMu[i_mu];
                            }
                    }
            }
        }

        // -- Select electron with highest pT -- //
        for(Int_t j_el=0; j_el<nQElectrons; j_el++)
        {
                Electron El = QElectronCollection[j_el];

                if( El.Pt > Pt_el )
                {
                        Pt_el		= El.Pt;
                        el_BestPair	= El;
                        Index_ele       = QIndexEle[j_el];
                }
        }

        //if( Pt_mu > 0 && Pt_el > 0 )
        if( Pt_mu > 0 && Pt_el > 0 && ( Pt_mu > LeadPtCut || Pt_el > LeadPtCut ) ) // At least one lepton has pT above 28 [GeV]
        {
                TLorentzVector reco_v1 = mu_BestPair.Momentum;
                TLorentzVector reco_v2 = el_BestPair.Momentum;
                Double_t reco_M = (reco_v1 + reco_v2).M();

                // -- 3D open angle -- //
                Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

                //if( reco_M > 10 && Angle < TMath::Pi() - 0.005 )
                if( reco_M > 10 )
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back( mu_BestPair );
                        SelectedElectronCollection->push_back( el_BestPair );
                        Sel_Index_Mu = Index_mu;
                        Sel_Index_Ele = Index_ele;
                }
        }

        return isPassEventSelection;
}

// test for emu event selection not using vtx cut
// Derived by Marijus Ambrozas 2018.08.07 to use LongSelectedEMu_t instead of NtupleHandle
Bool_t DYAnalyzer::EventSelection_emu_method_test(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, LongSelectedEMu_t *ntuple, // Input: electron, muon vectors and LongSelectedEMu_t
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
                if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10
                        //&& MuonCollection[j].Pt > LeadPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut )
                        && MuonCollection[j].Pt > SubPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut ) // pT>17 && |eta|<2.4
                {
                        QMuonCollection.push_back( MuonCollection[j] );
                        QIndexMu.push_back( j );
                }
        }

        //Collect qualified electrons among electrons
        vector< Electron > QElectronCollection;
        vector< Int_t > QIndexEle;
        for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
        {
                Electron elec = ElectronCollection[j];
                if( elec.passMediumID == kTRUE
                        //&& elec.Pt > SubPtCut && fabs(elec.etaSC) < 2.5 && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
                        && elec.Pt > SubPtCut && fabs(elec.etaSC) < LeadEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) ) // pT>17 && |eta|<2.4
                {
                        QElectronCollection.push_back( ElectronCollection[j] );
                        QIndexEle.push_back( j );
                }
        }

        Int_t nQMuons = (Int_t)QMuonCollection.size();
        Int_t nQElectrons = (Int_t)QElectronCollection.size();
        Int_t nQIndicesMu = (Int_t)QIndexMu.size();
        Int_t nQIndicesEle = (Int_t)QIndexEle.size();

        Double_t Pt_mu = 0;
        Double_t Pt_el = 0;
        Int_t Index_mu = -1;
        Int_t Index_ele = -1;
        Muon mu_BestPair;
        Electron el_BestPair;

        // -- Select muon with highest pT -- //
        if ( nQMuons == nQIndicesMu )
        {
            for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
            {
                    Muon Mu = QMuonCollection[i_mu];

                    // -- muon should be matched with HLT objects in emu best pair -- //
                    if( Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") )
                    {
                            if( Mu.Pt > Pt_mu )
                            {
                                    Pt_mu           = Mu.Pt;
                                    mu_BestPair     = Mu;
                                    Index_mu        = QIndexMu[i_mu];
                            }
                    }
            }
        }

        // -- Select electron with highest pT -- //
        for(Int_t j_el=0; j_el<nQElectrons; j_el++)
        {
                Electron El = QElectronCollection[j_el];

                if( El.Pt > Pt_el )
                {
                        Pt_el		= El.Pt;
                        el_BestPair	= El;
                        Index_ele       = QIndexEle[j_el];
                }
        }

        //if( Pt_mu > 0 && Pt_el > 0 )
        if( Pt_mu > 0 && Pt_el > 0 && ( Pt_mu > LeadPtCut || Pt_el > LeadPtCut ) ) // At least one lepton has pT above 28 [GeV]
        {
                TLorentzVector reco_v1 = mu_BestPair.Momentum;
                TLorentzVector reco_v2 = el_BestPair.Momentum;
                Double_t reco_M = (reco_v1 + reco_v2).M();

                // -- 3D open angle -- //
                Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

                //if( reco_M > 10 && Angle < TMath::Pi() - 0.005 )
                if( reco_M > 10 )
                {
                        isPassEventSelection = kTRUE;
                        SelectedMuonCollection->push_back( mu_BestPair );
                        SelectedElectronCollection->push_back( el_BestPair );
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
		if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10
			&& MuonCollection[j].Pt > LeadPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut )
			QMuonCollection.push_back( MuonCollection[j] );
	}

	//Collect qualified electrons among electrons
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.passMediumID == kTRUE
			//&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < 2.5 && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
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
		if( Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") )
		{
			// -- Start another loop for finding electron (for electron, we don't need to check about trigger) -- //
			for(Int_t j_el=0; j_el<nQElectrons; j_el++)
			{
				Electron El = QElectronCollection[j_el];

				Double_t VtxProb_temp = -999;
				Double_t VtxNormChi2_temp = 999;
				emuVertexProbNormChi2(ntuple, Mu.Inner_pT, El.gsfpT, &VtxProb_temp, &VtxNormChi2_temp);

				// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
				if( VtxNormChi2_temp < VtxNormChi2_BestPair )
				{
					VtxNormChi2_BestPair = VtxNormChi2_temp;
					mu_BestPair = Mu;
					el_BestPair = El;
				}
			} // -- end of the loop for j_el (finding for electron)
		}
	} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

	if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
	{
		TLorentzVector reco_v1 = mu_BestPair.Momentum;
		TLorentzVector reco_v2 = el_BestPair.Momentum;
		Double_t reco_M = (reco_v1 + reco_v2).M();

		// -- 3D open angle -- //
		Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

		if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( mu_BestPair );
			SelectedElectronCollection->push_back( el_BestPair );
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
	if( leadMu.Pt > LeadPtCut && fabs(leadMu.eta) < LeadEtaCut && 
		subMu.Pt  > SubPtCut  && fabs(subMu.eta)  < SubEtaCut )
		isPassAcc = kTRUE;

	return isPassAcc;
}

Bool_t DYAnalyzer::isPassAccCondition_GenLepton(GenLepton genlep1, GenLepton genlep2)
{
	Bool_t isPassAcc = kFALSE;

	GenLepton leadGenLep, subGenLep;
	CompareGenLepton(&genlep1, &genlep2, &leadGenLep, &subGenLep);
	
	if( leadGenLep.Pt > LeadPtCut && fabs(leadGenLep.eta) < LeadEtaCut &&
		subGenLep.Pt  > SubPtCut  && fabs(subGenLep.eta) < SubEtaCut )
		isPassAcc = 1;

	return isPassAcc;
}

void DYAnalyzer::CompareMuon(Muon *Mu1, Muon *Mu2, Muon *leadMu, Muon *subMu)
{
    if( Mu1->Pt > Mu2->Pt )
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
	if( genlep1->Pt > genlep2->Pt )
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

	if( NPt1 != NPt2 || NPt2 != NProb || NPt1 != NProb ) 
                std::cout << "NPt1: " << NPt1 << " NPt2: " << NPt2 << " Nprob: " << NProb << endl;

        // std::cout << "Pt1: " << Pt1 << " Pt2: " << Pt2 << endl;

	for(Int_t i=0; i<NProb; i++)
	{
                // std::cout << "\tPtCollection1->at("<< i << "): " << PtCollection1->at(i) << " PtCollection2->at("<< i << "): " << PtCollection2->at(i) << endl;
		if( ( PtCollection1->at(i) == Pt1 && PtCollection2->at(i) == Pt2 )  || ( PtCollection1->at(i) == Pt2 && PtCollection2->at(i) == Pt1 ) )
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

        if( NPt1 != NPt2 || NPt2 != NProb || NPt1 != NProb )
                std::cout << "NPt1: " << NPt1 << " NPt2: " << NPt2 << " Nprob: " << NProb << endl;

        // std::cout << "Pt1: " << Pt1 << " Pt2: " << Pt2 << endl;

        for(Int_t i=0; i<NProb; i++)
        {
                // std::cout << "\tPtCollection1->at("<< i << "): " << PtCollection1->at(i) << " PtCollection2->at("<< i << "): " << PtCollection2->at(i) << endl;
                if( ( PtCollection1->at(i) == Pt1 && PtCollection2->at(i) == Pt2 )  || ( PtCollection1->at(i) == Pt2 && PtCollection2->at(i) == Pt1 ) )
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
		if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 && elec.Pt > 15 )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
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
                // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
		//if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 
		if( elec.passMediumID == kTRUE // modified by Dalmin Pai
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		//if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

// -- Event selecton for the electron channel (to make SelectedEE_t) (2018.08.02) -- // derived by Marijus Ambrozas
Bool_t DYAnalyzer::EventSelection_ElectronChannel( vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in an event & NtupleHandle -- //
                                                vector< Electron >* SelectedElectronCollection,  // -- output: 2 electrons passing event selection conditions -- //
                                                vector< Int_t >* Sel_Index )    // -- output: 2 indexes of electrons that passed the selection -- //
{
        Bool_t isPassEventSelection = kFALSE;

        // -- Electron ID -- //
        vector< Electron > QElectronCollection;
        vector< Int_t > QIndex;

        for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
        {
            Electron elec = ElectronCollection[j];
            // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
            //if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1
            if( elec.passMediumID == kTRUE // modified by Dalmin Pai
                    && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
            {
                QElectronCollection.push_back( ElectronCollection[j] );
                QIndex.push_back( j );
            }
        }

        Int_t nQElectrons = (Int_t)QElectronCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if( nQElectrons == 2 && nQIndices == 2 )
        {
                Electron recolep1 = QElectronCollection[0];
                Electron recolep2 = QElectronCollection[1];

                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

                Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if( reco_M > 10 && isPassAcc == kTRUE )
                //if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
                {
                        isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back( recolep1 );
                        SelectedElectronCollection->push_back( recolep2 );
                        Sel_Index->push_back( QIndex[0] );
                        Sel_Index->push_back( QIndex[1] );
                }
        }
        return isPassEventSelection;

}

// -- Event re-selecton for the electron channel (from LongSelectedEE_t) (2018.08.06) -- // derived by Marijus Ambrozas
Bool_t DYAnalyzer::EventSelection_ElectronChannel( vector< Electron > ElectronCollection, LongSelectedEE_t *ntuple, // -- input: All electrons in an event & LongSelectedEE_t -- //
                                                vector< Electron >* SelectedElectronCollection,  // -- output: 2 electrons passing event selection conditions -- //
                                                vector< Int_t >* Sel_Index )    // -- output: 2 indexes of electrons that passed the selection -- //
{
        Bool_t isPassEventSelection = kFALSE;

        // -- Electron ID -- //
        vector< Electron > QElectronCollection;
        vector< Int_t > QIndex;

        for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
        {
            Electron elec = ElectronCollection[j];
            // std::cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
            //if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1
            if( elec.passMediumID == kTRUE // modified by Dalmin Pai
                    && elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
            {
                QElectronCollection.push_back( ElectronCollection[j] );
                QIndex.push_back( j );
            }
        }

        Int_t nQElectrons = (Int_t)QElectronCollection.size();
        Int_t nQIndices = (Int_t)QIndex.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

        if( nQElectrons == 2 && nQIndices == 2 )
        {
                Electron recolep1 = QElectronCollection[0];
                Electron recolep2 = QElectronCollection[1];

                Bool_t isPassAcc = kFALSE;
                isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

                Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

                if( reco_M > 10 && isPassAcc == kTRUE )
                //if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
                {
                        isPassEventSelection = kTRUE;
                        SelectedElectronCollection->push_back( recolep1 );
                        SelectedElectronCollection->push_back( recolep2 );
                        Sel_Index->push_back( QIndex[0] );
                        Sel_Index->push_back( QIndex[1] );
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
		//if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 
		//if( elec.passMediumID == kTRUE // modified by Dalmin Pai
		if( elec.isMediumElectron_2016dataFor80X_minus_Full5x5_SigmaIEtaIEta() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		//if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
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
		//if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 
		//if( elec.passMediumID == kTRUE // modified by Dalmin Pai
		if( elec.isMediumElectron_2016dataFor80X_minus_dEtaInSeed() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		//if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
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
		//if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 
		//if( elec.passMediumID == kTRUE // modified by Dalmin Pai
		if( elec.isMediumElectron_2016dataFor80X_minus_dPhiIn() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		//if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
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
		//if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 
		//if( elec.passMediumID == kTRUE // modified by Dalmin Pai
		if( elec.isMediumElectron_2016dataFor80X_minus_HoverE() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		//if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
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
		//if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 
		//if( elec.passMediumID == kTRUE // modified by Dalmin Pai
		if( elec.isMediumElectron_2016dataFor80X_minus_RelPFIso_Rho() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		//if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
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
		//if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 
		//if( elec.passMediumID == kTRUE // modified by Dalmin Pai
		if( elec.isMediumElectron_2016dataFor80X_minus_InvEminusInvP() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		//if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
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
		//if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 
		//if( elec.passMediumID == kTRUE // modified by Dalmin Pai
		if( elec.isMediumElectron_2016dataFor80X_minus_mHits() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		//if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
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
		//if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 
		//if( elec.passMediumID == kTRUE // modified by Dalmin Pai
		if( elec.isMediumElectron_2016dataFor80X_minus_passConvVeto() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		//if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
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
		if( elec.isMediumElectron_Spring25ns_minus_PFIso() && elec.ecalDriven == 1 
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

        // std::cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::isPassAccCondition_Electron(Electron Elec1, Electron Elec2)
{
	Bool_t isPassAcc = kFALSE;
	Electron leadElec, subElec;
	CompareElectron(&Elec1, &Elec2, &leadElec, &subElec);
	if( leadElec.Pt > LeadPtCut && fabs(leadElec.etaSC) < LeadEtaCut && !( fabs(leadElec.etaSC) > 1.4442 && fabs(leadElec.etaSC) < 1.566 ) &&
		subElec.Pt  > SubPtCut  && fabs(subElec.etaSC)  < SubEtaCut && !( fabs(subElec.etaSC) > 1.4442 && fabs(subElec.etaSC) < 1.566 ) )
		isPassAcc = kTRUE;

	return isPassAcc;
}


Bool_t DYAnalyzer::isPassAccCondition_GenLepton_ECALGAP(GenLepton genlep1, GenLepton genlep2)
{
	Bool_t isPassAcc = kFALSE;

	GenLepton leadGenLep, subGenLep;
	CompareGenLepton(&genlep1, &genlep2, &leadGenLep, &subGenLep);
	
	if( leadGenLep.Pt > LeadPtCut && fabs(leadGenLep.eta) < LeadEtaCut && !( fabs(leadGenLep.eta) > 1.4442 && fabs(leadGenLep.eta) < 1.566 ) &&
		subGenLep.Pt  > SubPtCut  && fabs(subGenLep.eta) < SubEtaCut && !( fabs(subGenLep.eta) > 1.4442 && fabs(subGenLep.eta) < 1.566 ) )
		isPassAcc = 1;

	return isPassAcc;
}

void DYAnalyzer::CompareElectron(Electron *Elec1, Electron *Elec2, Electron *leadElec, Electron *subElec)
{
    if( Elec1->Pt > Elec2->Pt )
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
		if( fabs(genlep.ID) == 22 && fabs(genlep.Mother) == 13)
		{
			
			Double_t dR = Calc_dR_GenLepton_GenOthers(*genlep_postFSR, genlep);

			// -- Sum of all photon's momentum near the post-FSR muon -- //
			if( dR < dRCut )
			{
				SumPhotonMom  = SumPhotonMom + genlep.Momentum;
				GenPhotonCollection->push_back( genlep );
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
		// if( fabs(genlep.ID) == 22 && fabs(genlep.Mother) == 13)
		if( fabs(genlep.ID) == 22 )
		{
			
			Double_t dR = Calc_dR_GenLepton_GenOthers(*genlep_postFSR, genlep);

			// -- Sum of all photon's momentum near the post-FSR muon -- //
			if( dR < dRCut )
			{
				SumPhotonMom  = SumPhotonMom + genlep.Momentum;
				GenPhotonCollection->push_back( genlep );
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


	if( isPassAcc_preFSREvent == kFALSE && isPassAcc_postFSREvent == kTRUE )
		FSRType = "A";

	else if( isPassAcc_preFSREvent == kTRUE && isPassAcc_postFSREvent == kTRUE)
		FSRType = "B";
	
	else if( isPassAcc_preFSREvent == kTRUE && isPassAcc_postFSREvent == kFALSE)
		FSRType = "C";

	else if( isPassAcc_preFSREvent == kFALSE && isPassAcc_postFSREvent == kFALSE)
		FSRType = "D";
	else
	{
                std::cout << "ERROR: NO FSR TYPE CORRESPONDING TO THIS EVENT" << endl;
		FSRType = "NOTAssigned";
	}

	return FSRType;
}

Double_t DYAnalyzer::Calc_dR_GenLeptons( GenLepton genlep1, GenLepton genlep2 )
{
	Double_t eta1 = genlep1.eta;
	Double_t phi1 = genlep1.phi;

	Double_t eta2 = genlep2.eta;
	Double_t phi2 = genlep2.phi;

	Double_t diff_eta = eta1 - eta2;
	Double_t diff_phi = phi1 - phi2;

	Double_t dR = sqrt( diff_eta * diff_eta + diff_phi * diff_phi );
	return dR;
}

Double_t DYAnalyzer::Calc_dR_GenLepton_GenOthers( GenLepton genlep1, GenOthers genlep2 )
{
	Double_t eta1 = genlep1.eta;
	Double_t phi1 = genlep1.phi;

	Double_t eta2 = genlep2.eta;
	Double_t phi2 = genlep2.phi;

	Double_t diff_eta = eta1 - eta2;
	Double_t diff_phi = phi1 - phi2;

	Double_t dR = sqrt( diff_eta * diff_eta + diff_phi * diff_phi );
	return dR;
}

void DYAnalyzer::GenMatching(TString MuonType, NtupleHandle* ntuple, vector<Muon>* MuonCollection)
{
	vector<GenLepton> GenLeptonCollection;
	Int_t NGenLeptons = ntuple->gnpair;

	if( MuonType == "PromptFinalState" )
	{
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isMuon() && genlep.isPromptFinalState )
				GenLeptonCollection.push_back( genlep );
		}
	}
	else if( MuonType == "fromTau")
	{
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isMuon() && genlep.isDirectPromptTauDecayProductFinalState )
				GenLeptonCollection.push_back( genlep );
		}

	}
	else if( MuonType == "fromHardProcess" )
	{
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isMuon() && genlep.fromHardProcessFinalState )
				GenLeptonCollection.push_back( genlep );
		}
	}
	else
	{
                std::cout << "Incorrect MuonType!" << endl;
		return;
	}

	//Give Acceptance Cuts
	if( GenLeptonCollection.size() >= 2 )
	{
		GenLepton leadGenLep, subGenLep;
		CompareGenLepton(&GenLeptonCollection[0], &GenLeptonCollection[1], &leadGenLep, &subGenLep);
		if( !(leadGenLep.Pt > LeadPtCut && subGenLep.Pt > SubPtCut && abs(leadGenLep.eta) < LeadEtaCut && abs(subGenLep.eta) < SubEtaCut) )
			GenLeptonCollection.clear();
	}


	
	Int_t NMuons = (Int_t)MuonCollection->size();
	vector<Muon> RecoMuonCollection;
	//Copy all muons in MuonCollection into RecoMuonCollection
	for(Int_t i_mu=0; i_mu<NMuons; i_mu++)
		RecoMuonCollection.push_back( MuonCollection->at(i_mu) );

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

			Double_t dR = sqrt( (gen_eta-reco_eta)*(gen_eta-reco_eta) + (gen_phi-reco_phi)*(gen_phi-reco_phi) );
			Double_t dPt = fabs(gen_Pt - reco_Pt);
			if( dR < 0.3 )
			{
				if( dPt < dPtMin )
				{
					i_matched = i_reco;
					dPtMin = dPt;
				}
			}
		}

		if( i_matched != -1 )
			MuonCollection->push_back( RecoMuonCollection[i_matched] );
	}

	return;
}

void DYAnalyzer::ConvertToTunePInfo( Muon &mu )
{
	// -- Use TuneP information -- //
	mu.Pt = mu.TuneP_pT;
	mu.eta = mu.TuneP_eta;
	mu.phi = mu.TuneP_phi;

	mu.Momentum.SetPtEtaPhiM( mu.Pt, mu.eta, mu.phi, M_Mu );
}

void DYAnalyzer::PrintOutDoubleMuInfo( Muon mu1, Muon mu2 )
{
	printf("\t[Muon1] (pT, eta, phi, charge) = (%10.5lf, %10.5lf, %10.5lf, %.1d)\n", mu1.Pt, mu1.eta, mu1.phi, mu1.charge);
	printf("\t[Muon2] (pT, eta, phi, charge) = (%10.5lf, %10.5lf, %10.5lf, %.1d)\n", mu2.Pt, mu2.eta, mu2.phi, mu2.charge);
	Double_t reco_M = ( mu1.Momentum + mu2.Momentum ).M();
	printf("\t\tDilepton Mass: %10.5lf\n", reco_M);

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
	    if( MuonCollection[j].isHighPtMuon_minus_dzVTX() )
	    {
	    	if( MuonCollection[j].trkiso < 0.10 )
	    		PassingMuonCollection.push_back( MuonCollection[j] );
	    	else
	    		FailingMuonCollection.push_back( MuonCollection[j] );
	    }
	}

	Int_t nFailMuon = (Int_t)FailingMuonCollection.size();

	if( nFailMuon >= 2 ) // -- Dijet events: contains more than 2 failing muons regardless of # passing muons -- // 
	{
		if( nFailMuon == 2 )
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
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

			Bool_t isOS = kFALSE;
			if( recolep1.charge != recolep2.charge ) isOS = kTRUE;

			// if( reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005 )
			if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}

		} // -- end of if( nFailMuon == 2 ) -- //
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

					if( j_mu != i_mu ) // -- do not calculate vertex variables(prob, chi2). with itself -- //
					{
						// -- Check that this pair is within acceptance -- //
						Bool_t isPassAcc = kFALSE;
						isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

						if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
						{
							Double_t VtxProb_temp = -999;
							Double_t VtxNormChi2_temp = 999;
							DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

							// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
							if( VtxNormChi2_temp < VtxNormChi2_BestPair )
							{
								VtxNormChi2_BestPair = VtxNormChi2_temp;
								mu1_BestPair = Mu;
								mu2_BestPair = Mu_jth;
							}
						}
					}
				} // -- end of the loop for j_mu (finding for second muon)
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

			if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

				Bool_t isOS = kFALSE;
				if( mu1_BestPair.charge != mu2_BestPair.charge ) isOS = kTRUE;

				if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- end of (# failing muons > 2) case -- //

	} // -- end of if( nFailMuon >= 2 ) -- //

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
	    if( MuonCollection[j].isHighPtMuon_minus_dzVTX() )
	    {
	    	if( MuonCollection[j].trkiso < 0.10 )
	    		PassingMuonCollection.push_back( MuonCollection[j] );
	    	else
	    		FailingMuonCollection.push_back( MuonCollection[j] );
	    }
	}

	Int_t nFailMuon = (Int_t)FailingMuonCollection.size();
	Int_t nPassMuon = (Int_t)PassingMuonCollection.size();

	if( nFailMuon == 1 && nPassMuon == 1) // -- W+Jets events: exactly (# pass muon , # fail muon ) = (1, 1) -- //
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
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

		Bool_t isOS = kFALSE;
		if( recolep1.charge != recolep2.charge ) isOS = kTRUE;

		// if( reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005 )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 ); // -- first one: passing muon -- //
			SelectedMuonCollection->push_back( recolep2 ); // -- second one: failing muon -- //
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
	    if( MuonCollection[j].isHighPtMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
		if( mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

	if( isExistHLTMatchedMuon == kTRUE )
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
		if( nQMuons == 2)
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
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

			Bool_t isOS = kFALSE;
			if( recolep1.charge != recolep2.charge ) isOS = kTRUE;

			// if( reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005 )
			if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}
		}
		else if( nQMuons > 2 )
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
				if( Mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

						if( j_mu != i_mu ) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

							if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
								if( VtxNormChi2_temp < VtxNormChi2_BestPair )
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

			if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

				Bool_t isOS = kFALSE;
				if( mu1_BestPair.charge != mu2_BestPair.charge ) isOS = kTRUE;

				if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- End of else if( nQMuons > 2 ) -- //

	} // -- End of if( isExistHLTMatchedMuon == kTRUE ) -- //

	return isPassEventSelection;
}

