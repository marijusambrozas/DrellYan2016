#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TROOT.h>
#include <TApplication.h>
#include <vector>
#include <TMath.h>
#include <TFormula.h>
#include <iostream>
#include <sstream>

// -- Customized Analyzer for EMu selection -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"

#define M_Mu 0.1056583715 // -- GeV -- //

// -- Muon Channel -- //
//void MakeReSelectedEMu(Int_t type, TString HLTname = "IsoMu24_OR_IsoTkMu24")
void MakeReSelectedEMu(TString HLTname = "IsoMu24_OR_IsoTkMu24")
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;
    TString OutputName = "ReSelectedEMu_WW_34.root";

//    TString DataType, DataLocation, DataLocation2, Type, OutputName;

//    if( type == 1 ) {
//            DataType = "B";
//            DataLocation = "SingleMuon_Run2016B";
//    }
//    else if( type == 2 ) {
//            DataType = "C";
//            DataLocation = "SingleMuon_Run2016C";
//    }
//    else if( type == 3 ) {
//            DataType = "D";
//            DataLocation = "SingleMuon_Run2016D";
//    }
//    else if( type == 4 ) {
//            DataType = "E";
//            DataLocation = "SingleMuon_Run2016E";
//    }
//    else if( type == 5 ) {
//            DataType = "F";
//            DataLocation = "SingleMuon_Run2016F";
//    }
//    else if( type == 6 ) {
//            DataType = "G";
//            DataLocation = "SingleMuon_Run2016G";
//    }
//    else if( type == 7 ) {
//            DataType = "H";
//            DataLocation = "SingleMuon_Run2016Hver2";
//            DataLocation2 = "SingleMuon_Run2016Hver3";
//    }

    Bool_t isMC = kTRUE;
//    if( type < 10  ) Type = "Data";
//    // -- Background MC samples -- //
//    else if( type == 21 ) Type = "ttbar";
//    else if( type == 22 ) Type = "ttbarBackup";
//    else if( type == 23 ) Type = "ttbar_M700toInf";
//    else if( type == 31 ) Type = "DYTauTau_M10to50";
//    else if( type == 32 ) Type = "DYTauTau_M50toInf";
//    else if( type == 41 ) Type = "VVnST";

    //Creating a file
    TFile* EMuFile = new TFile("/media/sf_DATA/"+OutputName, "RECREATE");
    TTree* EMuTree = new TTree("DYTree", "DYTree");

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

    // -- Each ntuple directory & corresponding Tags -- //
//    vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    Tag.push_back("WW");
    nEvents.push_back(6987123);
    Xsec.push_back(118.7);

    // -- Creating LongSelectedMuMu variables to assign branches -- //
    LongSelectedEMu_t EMu; EMu.CreateNew();

    EMuTree->Branch("nVertices", &EMu.nVertices);
    EMuTree->Branch("runNum", &EMu.runNum);
    EMuTree->Branch("lumiBlock", &EMu.lumiBlock);
    EMuTree->Branch("evtNum", &EMu.evtNum);
    EMuTree->Branch("nPileUp", &EMu.nPileUp);
    EMuTree->Branch("GENEvt_weight", &EMu.GENEvt_weight);
    EMuTree->Branch("HLT_ntrig", &EMu.HLT_ntrig);
    EMuTree->Branch("HLT_trigFired", &EMu.HLT_trigFired);
    EMuTree->Branch("HLT_trigName", &EMu.HLT_trigName);
//    EMuTree->Branch("HLT_trigPt", &EMu.HLT_trigPt);
    EMuTree->Branch("HLT_trigEta", &EMu.HLT_trigEta);
    EMuTree->Branch("HLT_trigPhi", &EMu.HLT_trigPhi);
    EMuTree->Branch("isHardProcess", &EMu.isHardProcess);
    EMuTree->Branch("EMu_InvM", &EMu.EMu_InvM);
    EMuTree->Branch("Muon_pT", &EMu.Muon_pT);
    EMuTree->Branch("Muon_eta", &EMu.Muon_eta);
    EMuTree->Branch("Muon_phi", &EMu.Muon_phi);
    EMuTree->Branch("isGLBmuon", &EMu.isGLBmuon);
    EMuTree->Branch("isPFmuon", &EMu.isPFmuon);
    EMuTree->Branch("isTRKmuon", &EMu.isTRKmuon);
    EMuTree->Branch("Muon_charge", &EMu.Muon_charge);
    EMuTree->Branch("Muon_chi2dof", &EMu.Muon_chi2dof);
    EMuTree->Branch("Muon_muonHits", &EMu.Muon_muonHits);
    EMuTree->Branch("Muon_nSegments", &EMu.Muon_nSegments);
    EMuTree->Branch("Muon_nMatches", &EMu.Muon_nMatches);
    EMuTree->Branch("Muon_trackerLayers", &EMu.Muon_trackerLayers);
    EMuTree->Branch("Muon_pixelHits", &EMu.Muon_pixelHits);
    EMuTree->Branch("Muon_dxyVTX", &EMu.Muon_dxyVTX);
    EMuTree->Branch("Muon_dzVTX", &EMu.Muon_dzVTX);
    EMuTree->Branch("Muon_trkiso", &EMu.Muon_trkiso);
    EMuTree->Branch("Muon_PfChargedHadronIsoR04", &EMu.Muon_PfChargedHadronIsoR04);
    EMuTree->Branch("Muon_PfNeutralHadronIsoR04", &EMu.Muon_PfNeutralHadronIsoR04);
    EMuTree->Branch("Muon_PfGammaIsoR04", &EMu.Muon_PfGammaIsoR04);
    EMuTree->Branch("Muon_PFSumPUIsoR04", &EMu.Muon_PFSumPUIsoR04);
    EMuTree->Branch("Muon_Px", &EMu.Muon_Px);
    EMuTree->Branch("Muon_Py", &EMu.Muon_Py);
    EMuTree->Branch("Muon_Pz", &EMu.Muon_Pz);
    EMuTree->Branch("Muon_Energy", &EMu.Muon_Energy);
    EMuTree->Branch("Muon_Best_pT", &EMu.Muon_Best_pT);
    EMuTree->Branch("Muon_Best_pTError", &EMu.Muon_Best_pTError);
    EMuTree->Branch("Muon_Best_Px", &EMu.Muon_Best_Px);
    EMuTree->Branch("Muon_Best_Py", &EMu.Muon_Best_Py);
    EMuTree->Branch("Muon_Best_Pz", &EMu.Muon_Best_Pz);
    EMuTree->Branch("Muon_Best_eta", &EMu.Muon_Best_eta);
    EMuTree->Branch("Muon_Best_phi", &EMu.Muon_Best_phi);
    EMuTree->Branch("Muon_Inner_pT", &EMu.Muon_Inner_pT);
    EMuTree->Branch("Muon_Inner_pTError", &EMu.Muon_Inner_pTError);
    EMuTree->Branch("Muon_Inner_Px", &EMu.Muon_Inner_Px);
    EMuTree->Branch("Muon_Inner_Py", &EMu.Muon_Inner_Py);
    EMuTree->Branch("Muon_Inner_Pz", &EMu.Muon_Inner_Pz);
    EMuTree->Branch("Muon_Inner_eta", &EMu.Muon_Inner_eta);
    EMuTree->Branch("Muon_Inner_phi", &EMu.Muon_Inner_phi);
    EMuTree->Branch("Muon_Outer_pT", &EMu.Muon_Outer_pT);
    EMuTree->Branch("Muon_Outer_pTError", &EMu.Muon_Outer_pTError);
    EMuTree->Branch("Muon_Outer_Px", &EMu.Muon_Outer_Px);
    EMuTree->Branch("Muon_Outer_Py", &EMu.Muon_Outer_Py);
    EMuTree->Branch("Muon_Outer_Pz", &EMu.Muon_Outer_Pz);
    EMuTree->Branch("Muon_Outer_eta", &EMu.Muon_Outer_eta);
    EMuTree->Branch("Muon_Outer_phi", &EMu.Muon_Outer_phi);
    EMuTree->Branch("Muon_GLB_pT", &EMu.Muon_GLB_pT);
    EMuTree->Branch("Muon_GLB_pTError", &EMu.Muon_GLB_pTError);
    EMuTree->Branch("Muon_GLB_Px", &EMu.Muon_GLB_Px);
    EMuTree->Branch("Muon_GLB_Py", &EMu.Muon_GLB_Py);
    EMuTree->Branch("Muon_GLB_Pz", &EMu.Muon_GLB_Pz);
    EMuTree->Branch("Muon_GLB_eta", &EMu.Muon_GLB_eta);
    EMuTree->Branch("Muon_GLB_phi", &EMu.Muon_GLB_phi);
    EMuTree->Branch("Muon_TuneP_pT", &EMu.Muon_TuneP_pT);
    EMuTree->Branch("Muon_TuneP_pTError", &EMu.Muon_TuneP_pTError);
    EMuTree->Branch("Muon_TuneP_Px", &EMu.Muon_TuneP_Px);
    EMuTree->Branch("Muon_TuneP_Py", &EMu.Muon_TuneP_Py);
    EMuTree->Branch("Muon_TuneP_Pz", &EMu.Muon_TuneP_Pz);
    EMuTree->Branch("Muon_TuneP_eta", &EMu.Muon_TuneP_eta);
    EMuTree->Branch("Muon_TuneP_phi", &EMu.Muon_TuneP_phi);
    EMuTree->Branch("Electron_pT", &EMu.Electron_pT);
    EMuTree->Branch("Electron_eta", &EMu.Electron_eta);
    EMuTree->Branch("Electron_phi", &EMu.Electron_phi);
    EMuTree->Branch("Electron_Energy", &EMu.Electron_Energy);
    EMuTree->Branch("Electron_charge", &EMu.Electron_charge);
    EMuTree->Branch("Electron_gsfpT", &EMu.Electron_gsfpT);
    EMuTree->Branch("Electron_gsfPx", &EMu.Electron_gsfPx);
    EMuTree->Branch("Electron_gsfPy", &EMu.Electron_gsfPy);
    EMuTree->Branch("Electron_gsfPz", &EMu.Electron_gsfPz);
    EMuTree->Branch("Electron_gsfEta", &EMu.Electron_gsfEta);
    EMuTree->Branch("Electron_gsfPhi", &EMu.Electron_gsfEta);
    EMuTree->Branch("Electron_gsfCharge", &EMu.Electron_gsfCharge);
    EMuTree->Branch("Electron_etaSC", &EMu.Electron_etaSC);
    EMuTree->Branch("Electron_phiSC", &EMu.Electron_phiSC);
    EMuTree->Branch("Electron_etaWidth", &EMu.Electron_etaWidth);
    EMuTree->Branch("Electron_phiWidth", &EMu.Electron_phiWidth);
    EMuTree->Branch("Electron_dEtaIn", &EMu.Electron_dEtaIn);
    EMuTree->Branch("Electron_dEtaInSeed", &EMu.Electron_dEtaInSeed);
    EMuTree->Branch("Electron_dPhiIn", &EMu.Electron_dPhiIn);
    EMuTree->Branch("Electron_sigmaIEtaIEta", &EMu.Electron_sigmaIEtaIEta);
    EMuTree->Branch("Electron_Full5x5_SigmaIEtaIEta", &EMu.Electron_Full5x5_SigmaIEtaIEta);
    EMuTree->Branch("Electron_HoverE", &EMu.Electron_HoverE);
    EMuTree->Branch("Electron_fbrem", &EMu.Electron_fbrem);
    EMuTree->Branch("Electron_eOverP", &EMu.Electron_eOverP);
    EMuTree->Branch("Electron_InvEminusInvP", &EMu.Electron_InvEminusInvP);
    EMuTree->Branch("Electron_dxyVTX", &EMu.Electron_dxyVTX);
    EMuTree->Branch("Electron_dzVTX", &EMu.Electron_dzVTX);
    EMuTree->Branch("Electron_dxy", &EMu.Electron_dxy);
    EMuTree->Branch("Electron_dz", &EMu.Electron_dz);
    EMuTree->Branch("Electron_dxyBS", &EMu.Electron_dxyBS);
    EMuTree->Branch("Electron_dzBS", &EMu.Electron_dxyBS);
    EMuTree->Branch("Electron_chIso03", &EMu.Electron_chIso03);
    EMuTree->Branch("Electron_nhIso03", &EMu.Electron_nhIso03);
    EMuTree->Branch("Electron_phIso03", &EMu.Electron_phIso03);
    EMuTree->Branch("Electron_ChIso03FromPU", &EMu.Electron_ChIso03FromPU);
    EMuTree->Branch("Electron_mHits", &EMu.Electron_mHits);
    EMuTree->Branch("Electron_EnergySC", &EMu.Electron_EnergySC);
    EMuTree->Branch("Electron_preEnergySC", &EMu.Electron_preEnergySC);
    EMuTree->Branch("Electron_rawEnergySC", &EMu.Electron_rawEnergySC);
    EMuTree->Branch("Electron_etSC", &EMu.Electron_etSC);
    EMuTree->Branch("Electron_E15", &EMu.Electron_E15);
    EMuTree->Branch("Electron_E25", &EMu.Electron_E25);
    EMuTree->Branch("Electron_E55", &EMu.Electron_E55);
    EMuTree->Branch("Electron_RelPFIso_dBeta", &EMu.Electron_RelPFIso_dBeta);
    EMuTree->Branch("Electron_RelPFIso_Rho", &EMu.Electron_RelPFIso_Rho);
    EMuTree->Branch("Electron_r9", &EMu.Electron_r9);
    EMuTree->Branch("Electron_ecalDriven", &EMu.Electron_ecalDriven);
    EMuTree->Branch("Electron_passConvVeto", &EMu.Electron_passConvVeto);
//    EMuTree->Branch("Electron_passLooseID", &EMu.Electron_passLooseID);
    EMuTree->Branch("Electron_passMediumID", &EMu.Electron_passMediumID);
//    EMuTree->Branch("Electron_passTightID", &EMu.Electron_passTightID);
//    EMuTree->Branch("Electron_passMVAID_WP80", &EMu.Electron_passMVAID_WP80);
//    EMuTree->Branch("Electron_passMVAID_WP90", &EMu.Electron_passMVAID_WP90);
//    EMuTree->Branch("Electron_passHEEPID", &EMu.Electron_passHEEPID);

    //Loop for all samples
//    const Int_t Ntup = ntupleDirectory.size();
    const Int_t Ntup = 1;       //Using just 1 ntuple for test
    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

//        cout << "\t<" << [i_tup] << ">" << endl;

        TChain *chain = new TChain("DYTree");
        chain->Add("/media/sf_DATA/LongSelectedEMu_WW_34.root");

        LongSelectedEMu_t *EMuInput = new LongSelectedEMu_t();
        EMuInput->CreateFromChain(chain);

        Double_t SumWeight = 0, SumWeight_Separated = 0;

        Int_t NEvents = chain->GetEntries();
//        Int_t NEvents = 10000;					// test using few events
        cout << "\t[Total Events: " << NEvents << "]" << endl;
        myProgressBar_t bar (NEvents);

        if ( EMuInput->File_Given == kFALSE ) cout << "Error: the chain was not assigned." << endl;
        else
        {
            for(Int_t i=0; i<NEvents; i++)
            {
                EMuInput->GetEvent(i);

                // -- Positive/Negative Gen-weights -- //
                EMuInput->GENEvt_weight < 0 ? EMu.GENEvt_weight = -1 : EMu.GENEvt_weight = 1;
                SumWeight += EMu.GENEvt_weight;

                // -- Normalization -- //
                Double_t TotWeight = EMu.GENEvt_weight;
                if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*EMu.GENEvt_weight;

                Bool_t TriggerFlag = kFALSE;
                TriggerFlag = EMuInput->isTriggered( analyzer->HLT );

                if( TriggerFlag == kTRUE && EMuInput->isHardProcess == kTRUE )
                {
                    EMu.isHardProcess = kTRUE;
                    // -- Reco level selection -- //
                    // Muon
                    vector< Muon > MuonCollection;
                    Muon mu;
                    mu.FillFromSelectedX(EMuInput);
                    analyzer->ConvertToTunePInfo( mu );
                    MuonCollection.push_back( mu );

                    // Electron
                    vector< Electron > ElectronCollection;
                    Electron ele;
                    ele.FillFromSelectedX(EMuInput);
                    ElectronCollection.push_back( ele );

                    // -- Event Selection -- //
                    vector< Muon > SelectedMuonCollection;
                    vector< Electron > SelectedElectronCollection;
                    Int_t Sel_Index_Mu, Sel_Index_Ele;  // Ntuple indexes of electron and muon that passed the selection
                    Bool_t isPassEventSelection = kFALSE;
                    isPassEventSelection = analyzer->EventSelection_emu_method_test(MuonCollection, ElectronCollection, EMuInput, &SelectedMuonCollection,
                                                                                    &SelectedElectronCollection, Sel_Index_Mu, Sel_Index_Ele);

                    if( isPassEventSelection == kTRUE && Sel_Index_Mu != -1 && Sel_Index_Ele != -1 )
                    {
                        Muon mu = SelectedMuonCollection[0];
                        Electron ele = SelectedElectronCollection[0];

                        if ( ele.Pt != EMuInput->Electron_pT || mu.Pt != EMuInput->Muon_TuneP_pT )
                        {
                            cout << "\n\tEntry no. " << i << ":  Error: pT values do not match!" << endl;
                            cout << "\t ele.Pt=" << ele.Pt << "\t Electron_pT=" << EMuInput->Electron_pT << endl;
                            cout << "\t mu.Pt=" << mu.Pt << "\t Muon_pT=" << EMuInput->Muon_pT << endl;
                            break;
                        }

                        EMu.nVertices = EMuInput->nVertices;
                        EMu.runNum = EMuInput->runNum;
                        EMu.lumiBlock = EMuInput->lumiBlock;
                        EMu.evtNum = EMuInput->evtNum;
                        EMu.nPileUp = EMuInput->nPileUp;
                        EMu.HLT_ntrig = EMuInput->HLT_ntrig;

                        EMu.EMu_InvM = (mu.Momentum + ele.Momentum).M();

                        Int_t zero_count = 0; // resolving if there is no more information in arrays

    //                    for (Int_t iter=0; iter<1000; iter++)
    //                    {
    //                        if(ntuple->HLT_trigEta[iter] && EMuInput->HLT_trigPhi[iter])
    //                        {
    //                            EMu.HLT_trigFired->push_back(ntuple->HLT_trigFired[iter]);
    //                            EMu.HLT_trigEta->push_back(ntuple->HLT_trigEta[iter]);
    //                            EMu.HLT_trigPhi->push_back(ntuple->HLT_trigPhi[iter]);
    //                            if(iter<ntuple->HLT_ntrig) EMu.HLT_trigName->push_back(ntuple->HLT_trigName->at(iter));
    //                        }
    //                        else
    //                        {
    //                            zero_count++;
    //                            if(zero_count>1) break;
    //                        }
    //                    }

                        for (Int_t iter=0; iter<EMuInput->HLT_ntrig; iter++)      // There are more trigEta and trigPhi values than nTrig suggests
                        {
                            EMu.HLT_trigFired->push_back(EMuInput->HLT_trigFired->at(iter));
                            EMu.HLT_trigEta->push_back(EMuInput->HLT_trigEta->at(iter));
                            EMu.HLT_trigPhi->push_back(EMuInput->HLT_trigPhi->at(iter));
                            EMu.HLT_trigName->push_back(EMuInput->HLT_trigName->at(iter));
                        }

                        EMu.Muon_pT = EMuInput->Muon_pT;
                        EMu.Muon_eta = EMuInput->Muon_eta;
                        EMu.Muon_phi = EMuInput->Muon_phi;
                        EMu.isGLBmuon = EMuInput->isGLBmuon;
                        EMu.isPFmuon = EMuInput->isPFmuon;
                        EMu.isTRKmuon = EMuInput->isTRKmuon;
                        EMu.Muon_charge = EMuInput->Muon_charge;
                        EMu.Muon_chi2dof = EMuInput->Muon_chi2dof;
                        EMu.Muon_muonHits = EMuInput->Muon_muonHits;
                        EMu.Muon_nSegments = EMuInput->Muon_nSegments;
                        EMu.Muon_nMatches = EMuInput->Muon_nMatches;
                        EMu.Muon_trackerLayers = EMuInput->Muon_trackerLayers;
                        EMu.Muon_pixelHits = EMuInput->Muon_pixelHits;
                        EMu.Muon_dxyVTX = EMuInput->Muon_dxyVTX;
                        EMu.Muon_dzVTX = EMuInput->Muon_dzVTX;
                        EMu.Muon_trkiso = EMuInput->Muon_trkiso;
                        EMu.Muon_PfChargedHadronIsoR04 = EMuInput->Muon_PfChargedHadronIsoR04;
                        EMu.Muon_PfNeutralHadronIsoR04 = EMuInput->Muon_PfNeutralHadronIsoR04;
                        EMu.Muon_PfGammaIsoR04 = EMuInput->Muon_PfGammaIsoR04;
                        EMu.Muon_PFSumPUIsoR04 = EMuInput->Muon_PFSumPUIsoR04;
                        EMu.Muon_Px = EMuInput->Muon_Px;
                        EMu.Muon_Py = EMuInput->Muon_Py;
                        EMu.Muon_Pz = EMuInput->Muon_Pz;
                        EMu.Muon_Energy = sqrt( EMuInput->Muon_Px*EMuInput->Muon_Px + EMuInput->Muon_Py
                                                *EMuInput->Muon_Py+ EMuInput->Muon_Pz*EMuInput->Muon_Pz + M_Mu*M_Mu );
                        EMu.Muon_Best_pT = EMuInput->Muon_Best_pT;
                        EMu.Muon_Best_pTError = EMuInput->Muon_Best_pTError;
                        EMu.Muon_Best_Px = EMuInput->Muon_Best_Px;
                        EMu.Muon_Best_Py = EMuInput->Muon_Best_Py;
                        EMu.Muon_Best_Pz = EMuInput->Muon_Best_Pz;
                        EMu.Muon_Best_eta = EMuInput->Muon_Best_eta;
                        EMu.Muon_Best_phi = EMuInput->Muon_Best_phi;
                        EMu.Muon_Inner_pT = EMuInput->Muon_Inner_pT;
                        EMu.Muon_Inner_pTError = EMuInput->Muon_Inner_pTError;
                        EMu.Muon_Inner_Px = EMuInput->Muon_Inner_Px;
                        EMu.Muon_Inner_Py = EMuInput->Muon_Inner_Py;
                        EMu.Muon_Inner_Pz = EMuInput->Muon_Inner_Pz;
                        EMu.Muon_Inner_eta = EMuInput->Muon_Inner_eta;
                        EMu.Muon_Inner_phi = EMuInput->Muon_Inner_phi;
                        EMu.Muon_Outer_pT = EMuInput->Muon_Outer_pT;
                        EMu.Muon_Outer_pTError = EMuInput->Muon_Outer_pTError;
                        EMu.Muon_Outer_Px = EMuInput->Muon_Outer_Px;
                        EMu.Muon_Outer_Py = EMuInput->Muon_Outer_Py;
                        EMu.Muon_Outer_Pz = EMuInput->Muon_Outer_Pz;
                        EMu.Muon_Outer_eta = EMuInput->Muon_Outer_eta;
                        EMu.Muon_Outer_phi = EMuInput->Muon_Outer_phi;
                        EMu.Muon_GLB_pT = EMuInput->Muon_GLB_pT;
                        EMu.Muon_GLB_pTError = EMuInput->Muon_GLB_pTError;
                        EMu.Muon_GLB_Px = EMuInput->Muon_GLB_Px;
                        EMu.Muon_GLB_Py = EMuInput->Muon_GLB_Py;
                        EMu.Muon_GLB_Pz = EMuInput->Muon_GLB_Pz;
                        EMu.Muon_GLB_eta = EMuInput->Muon_GLB_eta;
                        EMu.Muon_GLB_phi = EMuInput->Muon_GLB_phi;
                        EMu.Muon_TuneP_pT = EMuInput->Muon_TuneP_pT;
                        EMu.Muon_TuneP_pTError = EMuInput->Muon_TuneP_pTError;
                        EMu.Muon_TuneP_Px = EMuInput->Muon_TuneP_Px;
                        EMu.Muon_TuneP_Py = EMuInput->Muon_TuneP_Py;
                        EMu.Muon_TuneP_Pz = EMuInput->Muon_TuneP_Pz;
                        EMu.Muon_TuneP_eta = EMuInput->Muon_TuneP_eta;
                        EMu.Muon_TuneP_phi = EMuInput->Muon_TuneP_phi;
                        EMu.Electron_pT = EMuInput->Electron_pT;
                        EMu.Electron_eta = EMuInput->Electron_eta;
                        EMu.Electron_phi = EMuInput->Electron_phi;
                        EMu.Electron_Energy = EMuInput->Electron_Energy;
                        EMu.Electron_charge = EMuInput->Electron_charge;
                        EMu.Electron_gsfpT = EMuInput->Electron_gsfpT;
                        EMu.Electron_gsfPx = EMuInput->Electron_gsfPx;
                        EMu.Electron_gsfPy = EMuInput->Electron_gsfPy;
                        EMu.Electron_gsfPz = EMuInput->Electron_gsfPz;
                        EMu.Electron_gsfEta = EMuInput->Electron_gsfEta;
                        EMu.Electron_gsfPhi = EMuInput->Electron_gsfPhi;
                        EMu.Electron_gsfCharge = EMuInput->Electron_gsfCharge;
                        EMu.Electron_etaSC = EMuInput->Electron_etaSC;
                        EMu.Electron_phiSC = EMuInput->Electron_phiSC;
                        EMu.Electron_etaWidth = EMuInput->Electron_etaWidth;
                        EMu.Electron_phiWidth = EMuInput->Electron_phiWidth;
                        EMu.Electron_dEtaIn = EMuInput->Electron_dEtaIn;
                        EMu.Electron_dEtaInSeed = EMuInput->Electron_dEtaInSeed;
                        EMu.Electron_dPhiIn = EMuInput->Electron_dPhiIn;
                        EMu.Electron_sigmaIEtaIEta = EMuInput->Electron_sigmaIEtaIEta;
                        EMu.Electron_Full5x5_SigmaIEtaIEta = EMuInput->Electron_Full5x5_SigmaIEtaIEta;
                        EMu.Electron_HoverE = EMuInput->Electron_HoverE;
                        EMu.Electron_fbrem = EMuInput->Electron_fbrem;
                        EMu.Electron_eOverP = EMuInput->Electron_eOverP;
                        EMu.Electron_InvEminusInvP = EMuInput->Electron_InvEminusInvP;
                        EMu.Electron_dxyVTX = EMuInput->Electron_dxyVTX;
                        EMu.Electron_dzVTX = EMuInput->Electron_dzVTX;
                        EMu.Electron_dxy = EMuInput->Electron_dxy;
                        EMu.Electron_dz = EMuInput->Electron_dz;
                        EMu.Electron_dxyBS = EMuInput->Electron_dxyBS;
                        EMu.Electron_dzBS = EMuInput->Electron_dzBS;
                        EMu.Electron_chIso03 = EMuInput->Electron_chIso03;
                        EMu.Electron_nhIso03 = EMuInput->Electron_nhIso03;
                        EMu.Electron_phIso03 = EMuInput->Electron_phIso03;
                        EMu.Electron_ChIso03FromPU = EMuInput->Electron_ChIso03FromPU;
                        EMu.Electron_mHits = EMuInput->Electron_mHits;
                        EMu.Electron_EnergySC = EMuInput->Electron_EnergySC;
                        EMu.Electron_preEnergySC = EMuInput->Electron_preEnergySC;
                        EMu.Electron_rawEnergySC = EMuInput->Electron_rawEnergySC;
                        EMu.Electron_etSC = EMuInput->Electron_etSC;
                        EMu.Electron_E15 = EMuInput->Electron_E15;
                        EMu.Electron_E25 = EMuInput->Electron_E25;
                        EMu.Electron_E55 = EMuInput->Electron_E55;
                        EMu.Electron_RelPFIso_dBeta = EMuInput->Electron_RelPFIso_dBeta;
                        EMu.Electron_RelPFIso_Rho = EMuInput->Electron_RelPFIso_Rho;
                        EMu.Electron_r9 = EMuInput->Electron_r9;
                        EMu.Electron_ecalDriven = EMuInput->Electron_ecalDriven;
                        EMu.Electron_passConvVeto = EMuInput->Electron_passConvVeto;
    //                            EMu.Electron_passLooseID = EMuInput->Electron_passLooseID;
                        EMu.Electron_passMediumID = EMuInput->Electron_passMediumID;
    //                            EMu.Electron_passTightID = EMuInput->Electron_passTightID;
    //                            EMu.Electron_passMVAID_WP80 = EMuInput->Electron_passMVAID_WP80;
    //                            EMu.Electron_passMVAID_WP90 = EMuInput->Electron_passMVAID_WP90;
    //                            EMu.Electron_passHEEPID = EMuInput->Electron_passHEEPID;

                        EMuTree->Fill();

                        EMu.HLT_trigFired->clear();
                        EMu.HLT_trigName->clear();
                        EMu.HLT_trigEta->clear();
                        EMu.HLT_trigPhi->clear();


    //                    Double_t reco_Pt = (mu.Momentum + ele.Momentum).Pt();
    //                    Double_t reco_rapi = (mu.Momentum + ele.Momentum).Rapidity();

                    } // End of event selection

                } // End of if( isTriggered )

                bar.Draw(i);
            } // End of event iteration

        } // End of else()

        printf("\tTotal sum of weights: %.1lf\n", SumWeight);
        printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
        if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);
        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds" << endl;
        cout << "============================================================\n" << endl;

        cout << "Checking number of entries in output and input files:   ";
        if(chain->GetEntries()==EMuTree->GetEntries()) cout << "OK! (" << EMuTree->GetEntries() << ")" << endl;
        else cout << "Input: " << chain->GetEntries() << "     Output: " << EMuTree->GetEntries() << "\nCHECK SELECTION ALGORITHMS.\n" << endl;

    } //end of i_tup iteration

    EMuFile->cd();
    cout << "Writing into file...";
    EMuTree->Write();
    cout << "Finished." << endl << "Closing a file..." << endl;
    EMuFile->Close();
    if ( !EMuFile->IsOpen() ) cout << "File " << OutputName << " has been closed successfully." << endl;
    else cout << "FILE " << OutputName << " COULD NOT BE CLOSED!" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}
