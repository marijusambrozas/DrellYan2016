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
//void MakeLongSelectedEMu(Int_t type, TString HLTname = "IsoMu24_OR_IsoTkMu24")
void MakeLongSelectedEMu(TString HLTname = "IsoMu24_OR_IsoTkMu24")
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;
    TString OutputName = "LongSelectedEMu_WW_34.root";

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
//    cout << "Type: " << Type << endl;
//    if( type < 10 ) cout << "DataType: Run2016" << DataType << endl;

//    TString BaseLocation;
//    if( Type == "Data" ) BaseLocation = "/data9/DATA/DYntuple/v2.0";
//    else BaseLocation = "/data9/DATA/DYntuple/v2.1";
//    cout << "DATA location: " << BaseLocation << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

    // -- Each ntuple directory & corresponding Tags -- //
//    vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    Tag.push_back("WW");
    nEvents.push_back(6987123);
    Xsec.push_back(118.7);
//    if( Type == "Data" ) {
//        ntupleDirectory.push_back( "" ); Tag.push_back( "Data" ); // -- It will be filled later -- //
//    }
//    else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

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

        TChain *chain = new TChain("recoTree/DYTree");
//        //Set MC chain
//        if( isMC == kTRUE ) chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
//        //Set Data chain
//        else {
//            chain->Add(BaseLocation+"/"+DataLocation+"/*.root");
//            if(type==7) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
//        }
        chain->Add("/media/sf_DATA/WW_34.root/recoTree/DYTree;1"); // NEED A WAY TO TELL THE NUMBER OF CYCLES AND THEIR EXTENTION NAMES
        chain->Add("/media/sf_DATA/WW_34.root/recoTree/DYTree;2");

        NtupleHandle *ntuple = new NtupleHandle( chain );
        if( isMC == kTRUE ) {
            ntuple->TurnOnBranches_GenLepton(); // for all leptons
            ntuple->TurnOnBranches_GenOthers(); // for quarks
        }
        ntuple->TurnOnBranches_Electron();
        ntuple->TurnOnBranches_Muon();

        Double_t SumWeight = 0, SumWeight_Separated = 0;

        Int_t NEvents = chain->GetEntries();
//        Int_t NEvents = 10000;					// test using small events
        cout << "\t[Total Events: " << NEvents << "]" << endl;
        myProgressBar_t bar (NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            ntuple->GetEvent(i);

            // -- Positive/Negative Gen-weights -- //
            ntuple->GENEvt_weight < 0 ? EMu.GENEvt_weight = -1 : EMu.GENEvt_weight = 1;
            SumWeight += EMu.GENEvt_weight;

            // -- Separate DYLL samples -- //
            Bool_t GenFlag = kFALSE;
            GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

            // -- Separate ttbar samples -- //
            Bool_t GenFlag_top = kTRUE;

            // -- Normalization -- //
            Double_t TotWeight = EMu.GENEvt_weight;
            if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*EMu.GENEvt_weight;

            Bool_t TriggerFlag = kFALSE;
            TriggerFlag = ntuple->isTriggered( analyzer->HLT );

            if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
            {
                EMu.isHardProcess = kTRUE;
                // -- Reco level selection -- //
                // Muons
                vector< Muon > MuonCollection;
                Int_t N_Muons = ntuple->nMuon;
                for(Int_t i_reco=0; i_reco<N_Muons; i_reco++)
                {
                    Muon mu;
                    mu.FillFromNtuple(ntuple, i_reco);

                    // -- Convert to TuneP variables -- //
                    analyzer->ConvertToTunePInfo( mu );
                    MuonCollection.push_back( mu );
                }

                // Electrons
                vector< Electron > ElectronCollection;
                Int_t N_Electrons = ntuple->Nelectrons;
                for(Int_t i_reco=0; i_reco<N_Electrons; i_reco++)
                {
                        Electron ele;
                        ele.FillFromNtuple(ntuple, i_reco);

                        ElectronCollection.push_back( ele );
                }

                // -- Event Selection -- //
                vector< Muon > SelectedMuonCollection;
                vector< Electron > SelectedElectronCollection;
                Int_t Sel_Index_Mu, Sel_Index_Ele;  // Ntuple indexes of electron and muon that passed the selection
                Bool_t isPassEventSelection = kFALSE;
                isPassEventSelection = analyzer->EventSelection_emu_method_test(MuonCollection, ElectronCollection, ntuple, &SelectedMuonCollection,
                                                                                &SelectedElectronCollection, Sel_Index_Mu, Sel_Index_Ele);

                if( isPassEventSelection == kTRUE && Sel_Index_Mu != -1 && Sel_Index_Ele != -1 )
                {
                    Muon mu = SelectedMuonCollection[0];
                    Electron ele = SelectedElectronCollection[0];

                    EMu.nVertices = ntuple->nVertices;
                    EMu.runNum = ntuple->runNum;
                    EMu.lumiBlock = ntuple->lumiBlock;
                    EMu.evtNum = ntuple->evtNum;
                    EMu.nPileUp = ntuple->nPileUp;
                    EMu.HLT_ntrig = ntuple->HLT_ntrig;

                    EMu.EMu_InvM = (mu.Momentum + ele.Momentum).M();

                    Int_t zero_count = 0; // resolving if there is no more information in arrays

//                    for (Int_t iter=0; iter<1000; iter++)
//                    {
//                        if(ntuple->HLT_trigEta[iter] && ntuple->HLT_trigPhi[iter])
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

                    for (Int_t iter=0; iter<ntuple->HLT_ntrig; iter++)      // There are more trigEta and trigPhi values than nTrig suggests
                    {
                        EMu.HLT_trigFired->push_back(ntuple->HLT_trigFired[iter]);
                        EMu.HLT_trigEta->push_back(ntuple->HLT_trigEta[iter]);
                        EMu.HLT_trigPhi->push_back(ntuple->HLT_trigPhi[iter]);
                        EMu.HLT_trigName->push_back(ntuple->HLT_trigName->at(iter));
                    }

                    EMu.Muon_pT = ntuple->Muon_pT[Sel_Index_Mu];
                    EMu.Muon_eta = ntuple->Muon_eta[Sel_Index_Mu];
                    EMu.Muon_phi = ntuple->Muon_phi[Sel_Index_Mu];
                    EMu.isGLBmuon = ntuple->isGLBmuon[Sel_Index_Mu];
                    EMu.isPFmuon = ntuple->isPFmuon[Sel_Index_Mu];
                    EMu.isTRKmuon = ntuple->isTRKmuon[Sel_Index_Mu];
                    EMu.Muon_charge = ntuple->Muon_charge[Sel_Index_Mu];
                    EMu.Muon_chi2dof = ntuple->Muon_chi2dof[Sel_Index_Mu];
                    EMu.Muon_muonHits = ntuple->Muon_muonHits[Sel_Index_Mu];
                    EMu.Muon_nSegments = ntuple->Muon_nSegments[Sel_Index_Mu];
                    EMu.Muon_nMatches = ntuple->Muon_nMatches[Sel_Index_Mu];
                    EMu.Muon_trackerLayers = ntuple->Muon_trackerLayers[Sel_Index_Mu];
                    EMu.Muon_pixelHits = ntuple->Muon_pixelHits[Sel_Index_Mu];
                    EMu.Muon_dxyVTX = ntuple->Muon_dxyVTX[Sel_Index_Mu];
                    EMu.Muon_dzVTX = ntuple->Muon_dzVTX[Sel_Index_Mu];
                    EMu.Muon_trkiso = ntuple->Muon_trkiso[Sel_Index_Mu];
                    EMu.Muon_PfChargedHadronIsoR04 = ntuple->Muon_PfChargedHadronIsoR04[Sel_Index_Mu];
                    EMu.Muon_PfNeutralHadronIsoR04 = ntuple->Muon_PfNeutralHadronIsoR04[Sel_Index_Mu];
                    EMu.Muon_PfGammaIsoR04 = ntuple->Muon_PfGammaIsoR04[Sel_Index_Mu];
                    EMu.Muon_PFSumPUIsoR04 = ntuple->Muon_PFSumPUIsoR04[Sel_Index_Mu];
                    EMu.Muon_Px = ntuple->Muon_Px[Sel_Index_Mu];
                    EMu.Muon_Py = ntuple->Muon_Py[Sel_Index_Mu];
                    EMu.Muon_Pz = ntuple->Muon_Pz[Sel_Index_Mu];
                    EMu.Muon_Energy = sqrt( ntuple->Muon_Px[Sel_Index_Mu]*ntuple->Muon_Px[Sel_Index_Mu] + ntuple->Muon_Py[Sel_Index_Mu]*ntuple->Muon_Py[Sel_Index_Mu]
                                          + ntuple->Muon_Pz[Sel_Index_Mu]*ntuple->Muon_Pz[Sel_Index_Mu] + M_Mu*M_Mu );
                    EMu.Muon_Best_pT = ntuple->Muon_Best_pT[Sel_Index_Mu];
                    EMu.Muon_Best_pTError = ntuple->Muon_Best_pTError[Sel_Index_Mu];
                    EMu.Muon_Best_Px = ntuple->Muon_Best_Px[Sel_Index_Mu];
                    EMu.Muon_Best_Py = ntuple->Muon_Best_Py[Sel_Index_Mu];
                    EMu.Muon_Best_Pz = ntuple->Muon_Best_Pz[Sel_Index_Mu];
                    EMu.Muon_Best_eta = ntuple->Muon_Best_eta[Sel_Index_Mu];
                    EMu.Muon_Best_phi = ntuple->Muon_Best_phi[Sel_Index_Mu];
                    EMu.Muon_Inner_pT = ntuple->Muon_Inner_pT[Sel_Index_Mu];
                    EMu.Muon_Inner_pTError = ntuple->Muon_Inner_pTError[Sel_Index_Mu];
                    EMu.Muon_Inner_Px = ntuple->Muon_Inner_Px[Sel_Index_Mu];
                    EMu.Muon_Inner_Py = ntuple->Muon_Inner_Py[Sel_Index_Mu];
                    EMu.Muon_Inner_Pz = ntuple->Muon_Inner_Pz[Sel_Index_Mu];
                    EMu.Muon_Inner_eta = ntuple->Muon_Inner_eta[Sel_Index_Mu];
                    EMu.Muon_Inner_phi = ntuple->Muon_Inner_phi[Sel_Index_Mu];
                    EMu.Muon_Outer_pT = ntuple->Muon_Outer_pT[Sel_Index_Mu];
                    EMu.Muon_Outer_pTError = ntuple->Muon_Outer_pTError[Sel_Index_Mu];
                    EMu.Muon_Outer_Px = ntuple->Muon_Outer_Px[Sel_Index_Mu];
                    EMu.Muon_Outer_Py = ntuple->Muon_Outer_Py[Sel_Index_Mu];
                    EMu.Muon_Outer_Pz = ntuple->Muon_Outer_Pz[Sel_Index_Mu];
                    EMu.Muon_Outer_eta = ntuple->Muon_Outer_eta[Sel_Index_Mu];
                    EMu.Muon_Outer_phi = ntuple->Muon_Outer_phi[Sel_Index_Mu];
                    EMu.Muon_GLB_pT = ntuple->Muon_GLB_pT[Sel_Index_Mu];
                    EMu.Muon_GLB_pTError = ntuple->Muon_GLB_pTError[Sel_Index_Mu];
                    EMu.Muon_GLB_Px = ntuple->Muon_GLB_Px[Sel_Index_Mu];
                    EMu.Muon_GLB_Py = ntuple->Muon_GLB_Py[Sel_Index_Mu];
                    EMu.Muon_GLB_Pz = ntuple->Muon_GLB_Pz[Sel_Index_Mu];
                    EMu.Muon_GLB_eta = ntuple->Muon_GLB_eta[Sel_Index_Mu];
                    EMu.Muon_GLB_phi = ntuple->Muon_GLB_phi[Sel_Index_Mu];
                    EMu.Muon_TuneP_pT = ntuple->Muon_TuneP_pT[Sel_Index_Mu];
                    EMu.Muon_TuneP_pTError = ntuple->Muon_TuneP_pTError[Sel_Index_Mu];
                    EMu.Muon_TuneP_Px = ntuple->Muon_TuneP_Px[Sel_Index_Mu];
                    EMu.Muon_TuneP_Py = ntuple->Muon_TuneP_Py[Sel_Index_Mu];
                    EMu.Muon_TuneP_Pz = ntuple->Muon_TuneP_Pz[Sel_Index_Mu];
                    EMu.Muon_TuneP_eta = ntuple->Muon_TuneP_eta[Sel_Index_Mu];
                    EMu.Muon_TuneP_phi = ntuple->Muon_TuneP_phi[Sel_Index_Mu];
                    EMu.Electron_pT = ntuple->Electron_pT[Sel_Index_Ele];
                    EMu.Electron_eta = ntuple->Electron_eta[Sel_Index_Ele];
                    EMu.Electron_phi = ntuple->Electron_phi[Sel_Index_Ele];
                    EMu.Electron_Energy = ntuple->Electron_Energy[Sel_Index_Ele];
                    EMu.Electron_charge = ntuple->Electron_charge[Sel_Index_Ele];
                    EMu.Electron_gsfpT = ntuple->Electron_gsfpT[Sel_Index_Ele];
                    EMu.Electron_gsfPx = ntuple->Electron_gsfPx[Sel_Index_Ele];
                    EMu.Electron_gsfPy = ntuple->Electron_gsfPy[Sel_Index_Ele];
                    EMu.Electron_gsfPz = ntuple->Electron_gsfPz[Sel_Index_Ele];
                    EMu.Electron_gsfEta = ntuple->Electron_gsfEta[Sel_Index_Ele];
                    EMu.Electron_gsfPhi = ntuple->Electron_gsfPhi[Sel_Index_Ele];
                    EMu.Electron_gsfCharge = ntuple->Electron_gsfCharge[Sel_Index_Ele];
                    EMu.Electron_etaSC = ntuple->Electron_etaSC[Sel_Index_Ele];
                    EMu.Electron_phiSC = ntuple->Electron_phiSC[Sel_Index_Ele];
                    EMu.Electron_etaWidth = ntuple->Electron_etaWidth[Sel_Index_Ele];
                    EMu.Electron_phiWidth = ntuple->Electron_phiWidth[Sel_Index_Ele];
                    EMu.Electron_dEtaIn = ntuple->Electron_dEtaIn[Sel_Index_Ele];
                    EMu.Electron_dEtaInSeed = ntuple->Electron_dEtaInSeed[Sel_Index_Ele];
                    EMu.Electron_dPhiIn = ntuple->Electron_dPhiIn[Sel_Index_Ele];
                    EMu.Electron_sigmaIEtaIEta = ntuple->Electron_sigmaIEtaIEta[Sel_Index_Ele];
                    EMu.Electron_Full5x5_SigmaIEtaIEta = ntuple->Electron_Full5x5_SigmaIEtaIEta[Sel_Index_Ele];
                    EMu.Electron_HoverE = ntuple->Electron_HoverE[Sel_Index_Ele];
                    EMu.Electron_fbrem = ntuple->Electron_fbrem[Sel_Index_Ele];
                    EMu.Electron_eOverP = ntuple->Electron_eOverP[Sel_Index_Ele];
                    EMu.Electron_InvEminusInvP = ntuple->Electron_InvEminusInvP[Sel_Index_Ele];
                    EMu.Electron_dxyVTX = ntuple->Electron_dxyVTX[Sel_Index_Ele];
                    EMu.Electron_dzVTX = ntuple->Electron_dzVTX[Sel_Index_Ele];
                    EMu.Electron_dxy = ntuple->Electron_dxy[Sel_Index_Ele];
                    EMu.Electron_dz = ntuple->Electron_dz[Sel_Index_Ele];
                    EMu.Electron_dxyBS = ntuple->Electron_dxyBS[Sel_Index_Ele];
                    EMu.Electron_dzBS = ntuple->Electron_dzBS[Sel_Index_Ele];
                    EMu.Electron_chIso03 = ntuple->Electron_chIso03[Sel_Index_Ele];
                    EMu.Electron_nhIso03 = ntuple->Electron_nhIso03[Sel_Index_Ele];
                    EMu.Electron_phIso03 = ntuple->Electron_phIso03[Sel_Index_Ele];
                    EMu.Electron_ChIso03FromPU = ntuple->Electron_ChIso03FromPU[Sel_Index_Ele];
                    EMu.Electron_mHits = ntuple->Electron_mHits[Sel_Index_Ele];
                    EMu.Electron_EnergySC = ntuple->Electron_EnergySC[Sel_Index_Ele];
                    EMu.Electron_preEnergySC = ntuple->Electron_preEnergySC[Sel_Index_Ele];
                    EMu.Electron_rawEnergySC = ntuple->Electron_rawEnergySC[Sel_Index_Ele];
                    EMu.Electron_etSC = ntuple->Electron_etSC[Sel_Index_Ele];
                    EMu.Electron_E15 = ntuple->Electron_E15[Sel_Index_Ele];
                    EMu.Electron_E25 = ntuple->Electron_E25[Sel_Index_Ele];
                    EMu.Electron_E55 = ntuple->Electron_E55[Sel_Index_Ele];
                    EMu.Electron_RelPFIso_dBeta = ntuple->Electron_RelPFIso_dBeta[Sel_Index_Ele];
                    EMu.Electron_RelPFIso_Rho = ntuple->Electron_RelPFIso_Rho[Sel_Index_Ele];
                    EMu.Electron_r9 = ntuple->Electron_r9[Sel_Index_Ele];
                    EMu.Electron_ecalDriven = ntuple->Electron_ecalDriven[Sel_Index_Ele];
                    EMu.Electron_passConvVeto = ntuple->Electron_passConvVeto[Sel_Index_Ele];
//                            EMu.Electron_passLooseID = ntuple->Electron_passLooseID[Sel_Index_Ele];
                    EMu.Electron_passMediumID = ntuple->Electron_passMediumID[Sel_Index_Ele];
//                            EMu.Electron_passTightID = ntuple->Electron_passTightID[Sel_Index_Ele];
//                            EMu.Electron_passMVAID_WP80 = ntuple->Electron_passMVAID_WP80[Sel_Index_Ele];
//                            EMu.Electron_passMVAID_WP90 = ntuple->Electron_passMVAID_WP90[Sel_Index_Ele];
//                            EMu.Electron_passHEEPID = ntuple->Electron_passHEEPID[Sel_Index_Ele];

                    EMuTree->Fill();

                    EMu.HLT_trigFired->clear();
                    EMu.HLT_trigName->clear();
                    EMu.HLT_trigEta->clear();
                    EMu.HLT_trigPhi->clear();


//                    Double_t reco_Pt = (mu.Momentum + ele.Momentum).Pt();
//                    Double_t reco_rapi = (mu.Momentum + ele.Momentum).Rapidity();

                } // End of event selection

            } //End of if( isTriggered )

            bar.Draw(i);
        } //End of event iteration

        printf("\tTotal sum of weights: %.1lf\n", SumWeight);
        printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
        if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);

        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

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
