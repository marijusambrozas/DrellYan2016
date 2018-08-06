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

// -- Customized Analyzer for EE selection -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"

// -- Electron Channel -- //
//void MakeLongSelectedEE(Int_t type, Int_t Num = 100, Int_t isTopPtReweighting = 0, TString HLTname = "Ele23Ele12")
void MakeLongSelectedEE(TString HLTname = "Ele23Ele12")
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;
    TString OutputName = "LongSelectedEE_ZToEE_M4500to6000_2.root";

//    TString DataType, DataLocation, DataLocation2, Type;

//    if( type == 1 ) {
//            DataType = "B";
//            DataLocation = "DoubleEG_Run2016B";
//    }
//    if( type == 2 ) {
//            DataType = "C";
//            DataLocation = "DoubleEG_Run2016C";
//    }
//    if( type == 3 ) {
//            DataType = "D";
//            DataLocation = "DoubleEG_Run2016D";
//    }
//    if( type == 4 ) {
//            DataType = "E";
//            DataLocation = "DoubleEG_Run2016E";
//    }
//    if( type == 5 ) {
//            DataType = "F";
//            DataLocation = "DoubleEG_Run2016F";
//    }
//    if( type == 6 ) {
//            DataType = "G";
//            DataLocation = "DoubleEG_Run2016G";
//    }
//    if( type == 7 ) {
//            DataType = "H";
//            DataLocation = "DoubleEG_Run2016Hver2";
//            DataLocation2 = "DoubleEG_Run2016Hver3";
//    }

    Bool_t isMC = kTRUE;
//    if( type < 10  ) {Type = "Data"; isMC = kFALSE;}
//    // -- Signal MC samples -- //
//    if( type == 11 ) Type = "DYEE_M10to50";
//    if( type == 12 ) Type = "DYEE_M50to100";
//    if( type == 13 ) Type = "DYEE_M100toInf";
//    // -- Background MC samples -- //
//    if( type == 21 ) Type = "ttbar";
//    if( type == 22 ) Type = "ttbarBackup";
//    if( type == 23 ) Type = "ttbar_M700toInf";
//    if( type == 31 ) Type = "DYTauTau_M10to50";
//    if( type == 32 ) Type = "DYTauTau_M50toInf";
//    if( type == 41 ) Type = "VVnST";
//    if( type == 51 ) Type = "WJetsToLNu";
//    //if( type == 61 ) Type = "QCDEMEnriched";

    //Creating a file
    TFile* ElectronFile = new TFile("/media/sf_DATA/"+OutputName, "RECREATE");
    TTree* ElectronTree = new TTree("DYTree", "DYTree");


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
    Tag.push_back("ZToEE_M4500to6000"); // There is no such type mentioned in DYAnalyzer
    nEvents.push_back(74606);       // The actual number of events in file is taken
    Xsec.push_back(4.56E-07);       // Copied from ZToMuMu_M4500to6000
//    if( Type == "Data" ) {
//            ntupleDirectory.push_back( "" ); Tag.push_back( "Data" ); // -- It will be filled later -- //
//    }
//    else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

    // -- Creating LongSelectedEE variables to assign branches -- //
    LongSelectedEE_t EE; EE.CreateNew();

    ElectronTree->Branch("nVertices", &EE.nVertices);
    ElectronTree->Branch("runNum", &EE.runNum);
    ElectronTree->Branch("lumiBlock", &EE.lumiBlock);
    ElectronTree->Branch("evtNum", &EE.evtNum);
    ElectronTree->Branch("nPileUp", &EE.nPileUp);
    ElectronTree->Branch("GENEvt_weight", &EE.GENEvt_weight);
    ElectronTree->Branch("HLT_ntrig", &EE.HLT_ntrig);
    ElectronTree->Branch("HLT_trigFired", &EE.HLT_trigFired);
    ElectronTree->Branch("HLT_trigName", &EE.HLT_trigName);
//    ElectronTree->Branch("HLT_trigPt", &EE.HLT_trigPt);
    ElectronTree->Branch("HLT_trigEta", &EE.HLT_trigEta);
    ElectronTree->Branch("HLT_trigPhi", &EE.HLT_trigPhi);
    ElectronTree->Branch("isHardProcess", &EE.isHardProcess);
    ElectronTree->Branch("Electron_InvM", &EE.Electron_InvM);
    ElectronTree->Branch("Electron_pT", &EE.Electron_pT);
    ElectronTree->Branch("Electron_eta", &EE.Electron_eta);
    ElectronTree->Branch("Electron_phi", &EE.Electron_phi);
    ElectronTree->Branch("Electron_Energy", &EE.Electron_Energy);
    ElectronTree->Branch("Electron_charge", &EE.Electron_charge);
    ElectronTree->Branch("Electron_gsfpT", &EE.Electron_gsfpT);
    ElectronTree->Branch("Electron_gsfPx", &EE.Electron_gsfPx);
    ElectronTree->Branch("Electron_gsfPy", &EE.Electron_gsfPy);
    ElectronTree->Branch("Electron_gsfPz", &EE.Electron_gsfPz);
    ElectronTree->Branch("Electron_gsfEta", &EE.Electron_gsfEta);
    ElectronTree->Branch("Electron_gsfPhi", &EE.Electron_gsfEta);
    ElectronTree->Branch("Electron_gsfCharge", &EE.Electron_gsfCharge);
    ElectronTree->Branch("Electron_etaSC", &EE.Electron_etaSC);
    ElectronTree->Branch("Electron_phiSC", &EE.Electron_phiSC);
    ElectronTree->Branch("Electron_etaWidth", &EE.Electron_etaWidth);
    ElectronTree->Branch("Electron_phiWidth", &EE.Electron_phiWidth);
    ElectronTree->Branch("Electron_dEtaIn", &EE.Electron_dEtaIn);
    ElectronTree->Branch("Electron_dEtaInSeed", &EE.Electron_dEtaInSeed);
    ElectronTree->Branch("Electron_dPhiIn", &EE.Electron_dPhiIn);
    ElectronTree->Branch("Electron_sigmaIEtaIEta", &EE.Electron_sigmaIEtaIEta);
    ElectronTree->Branch("Electron_Full5x5_SigmaIEtaIEta", &EE.Electron_Full5x5_SigmaIEtaIEta);
    ElectronTree->Branch("Electron_HoverE", &EE.Electron_HoverE);
    ElectronTree->Branch("Electron_fbrem", &EE.Electron_fbrem);
    ElectronTree->Branch("Electron_eOverP", &EE.Electron_eOverP);
    ElectronTree->Branch("Electron_InvEminusInvP", &EE.Electron_InvEminusInvP);
    ElectronTree->Branch("Electron_dxyVTX", &EE.Electron_dxyVTX);
    ElectronTree->Branch("Electron_dzVTX", &EE.Electron_dzVTX);
    ElectronTree->Branch("Electron_dxy", &EE.Electron_dxy);
    ElectronTree->Branch("Electron_dz", &EE.Electron_dz);
    ElectronTree->Branch("Electron_dxyBS", &EE.Electron_dxyBS);
    ElectronTree->Branch("Electron_dzBS", &EE.Electron_dxyBS);
    ElectronTree->Branch("Electron_chIso03", &EE.Electron_chIso03);
    ElectronTree->Branch("Electron_nhIso03", &EE.Electron_nhIso03);
    ElectronTree->Branch("Electron_phIso03", &EE.Electron_phIso03);
    ElectronTree->Branch("Electron_ChIso03FromPU", &EE.Electron_ChIso03FromPU);
    ElectronTree->Branch("Electron_mHits", &EE.Electron_mHits);
    ElectronTree->Branch("Electron_EnergySC", &EE.Electron_EnergySC);
    ElectronTree->Branch("Electron_preEnergySC", &EE.Electron_preEnergySC);
    ElectronTree->Branch("Electron_rawEnergySC", &EE.Electron_rawEnergySC);
    ElectronTree->Branch("Electron_etSC", &EE.Electron_etSC);
    ElectronTree->Branch("Electron_E15", &EE.Electron_E15);
    ElectronTree->Branch("Electron_E25", &EE.Electron_E25);
    ElectronTree->Branch("Electron_E55", &EE.Electron_E55);
    ElectronTree->Branch("Electron_RelPFIso_dBeta", &EE.Electron_RelPFIso_dBeta);
    ElectronTree->Branch("Electron_RelPFIso_Rho", &EE.Electron_RelPFIso_Rho);
    ElectronTree->Branch("Electron_r9", &EE.Electron_r9);
    ElectronTree->Branch("Electron_ecalDriven", &EE.Electron_ecalDriven);
    ElectronTree->Branch("Electron_passConvVeto", &EE.Electron_passConvVeto);
//    ElectronTree->Branch("Electron_passLooseID", &EE.Electron_passLooseID);
    ElectronTree->Branch("Electron_passMediumID", &EE.Electron_passMediumID);
//    ElectronTree->Branch("Electron_passTightID", &EE.Electron_passTightID);
//    ElectronTree->Branch("Electron_passMVAID_WP80", &EE.Electron_passMVAID_WP80);
//    ElectronTree->Branch("Electron_passMVAID_WP90", &EE.Electron_passMVAID_WP90);
//    ElectronTree->Branch("Electron_passHEEPID", &EE.Electron_passHEEPID);

    //Loop for all samples
//	const Int_t Ntup = ntupleDirectory.size();
    const Int_t Ntup = 1;       //Using just 1 ntuple for test
    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

//        cout << "\t<" << Tag[i_tup] << ">" << endl;

        TChain *chain = new TChain("recoTree/DYTree");
        chain->Add("/media/sf_DATA/ZToEE_M4500to6000_2.root/recoTree/DYTree;7");
        chain->Add("/media/sf_DATA/ZToEE_M4500to6000_2.root/recoTree/DYTree;8");

//        //Set MC chain
//        if( isMC == kTRUE ) chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
//        //Set Data chain
//        else {
//                chain->Add(BaseLocation+"/"+DataLocation+"/*.root");
//                if( type == 7 ) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
//        }

        NtupleHandle *ntuple = new NtupleHandle( chain );
        if( isMC == kTRUE ) {
                ntuple->TurnOnBranches_GenLepton(); // for all leptons
                ntuple->TurnOnBranches_GenOthers(); // for quarks
        }
        ntuple->TurnOnBranches_Electron();

        Double_t SumWeight = 0, SumWeight_Separated = 0;

        Int_t NEvents = chain->GetEntries();
//        Int_t NEvents = 10000;		// test using few events
        cout << "\t[Total Events: " << NEvents << "]" << endl;
        myProgressBar_t bar (NEvents);
        Int_t timesPassed = 0;

        for(Int_t i=0; i<NEvents; i++)
        {
            ntuple->GetEvent(i);

            // -- Positive/Negative Gen-weights -- //
            ntuple->GENEvt_weight < 0 ? EE.GENEvt_weight = -1 : EE.GENEvt_weight = 1;
            SumWeight += EE.GENEvt_weight;

            // -- Separate DYLL samples -- //
            Bool_t GenFlag = kFALSE;
            GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

            // -- Separate ttbar samples -- //
            Bool_t GenFlag_top = kTRUE;
            //Bool_t GenFlag_top = kFALSE;
            //vector<GenOthers> GenTopCollection;
            //GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

            // -- Normalization -- //
            Double_t TotWeight = EE.GENEvt_weight;
            if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*EE.GENEvt_weight;

            Bool_t TriggerFlag = kFALSE;
            TriggerFlag = ntuple->isTriggered( analyzer->HLT );

            if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
            {
                EE.isHardProcess = kTRUE;
                // -- Reco level selection -- //
                vector< Electron > ElectronCollection;
                Int_t NLeptons = ntuple->Nelectrons;
                for(Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                {
                    Electron ele;
                    ele.FillFromNtuple(ntuple, i_reco);
                    ElectronCollection.push_back( ele );
                }

                // -- Event Selection -- //
                vector< Electron > SelectedElectronCollection;
                vector< Int_t > Sel_Index; // Ntuple indexes of electrons that passed the selection
                Bool_t isPassEventSelection = kFALSE;
                isPassEventSelection = analyzer->EventSelection_ElectronChannel(ElectronCollection, ntuple, &SelectedElectronCollection, &Sel_Index);

                if( isPassEventSelection == kTRUE )
                {
                    timesPassed++;
                    Electron ele1 = SelectedElectronCollection[0];
                    Electron ele2 = SelectedElectronCollection[1];

                    EE.nVertices = ntuple->nVertices;
                    EE.runNum = ntuple->runNum;
                    EE.lumiBlock = ntuple->lumiBlock;
                    EE.evtNum = ntuple->evtNum;
                    EE.nPileUp = ntuple->nPileUp;
                    EE.HLT_ntrig = ntuple->HLT_ntrig;

                    EE.Electron_InvM = (ele1.Momentum + ele2.Momentum).M();

                    Int_t zero_count = 0; // resolving if there is no more information in arrays

//                    for (Int_t iter=0; iter<1000; iter++)
//                    {
//                        if(ntuple->HLT_trigEta[iter] && ntuple->HLT_trigPhi[iter])
//                        {
//                            EE.HLT_trigFired->push_back(ntuple->HLT_trigFired[iter]);
//                            EE.HLT_trigEta->push_back(ntuple->HLT_trigEta[iter]);
//                            EE.HLT_trigPhi->push_back(ntuple->HLT_trigPhi[iter]);
//                            if(iter<ntuple->HLT_ntrig) EE.HLT_trigName->push_back(ntuple->HLT_trigName->at(iter));
//                        }
//                        else
//                        {
//                            zero_count++;
//                            if(zero_count>1) break;
//                        }
//                    }

                    for (Int_t iter=0; iter<ntuple->HLT_ntrig; iter++)
                    {
                        EE.HLT_trigFired->push_back(ntuple->HLT_trigFired[iter]);
                        EE.HLT_trigEta->push_back(ntuple->HLT_trigEta[iter]);
                        EE.HLT_trigPhi->push_back(ntuple->HLT_trigPhi[iter]);
                        EE.HLT_trigName->push_back(ntuple->HLT_trigName->at(iter));
                    }

                    if(Sel_Index.size()!=2) cout << "======== ERROR: The number of electrons saved is not 2 ========" << endl;
                    else
                    {
                        for (UInt_t iter=0; iter<Sel_Index.size(); iter++)
                        {
                            Int_t index = Sel_Index[iter];

                            EE.Electron_pT->push_back(ntuple->Electron_pT[index]);
                            EE.Electron_eta->push_back(ntuple->Electron_eta[index]);
                            EE.Electron_phi->push_back(ntuple->Electron_phi[index]);
                            EE.Electron_Energy->push_back(ntuple->Electron_Energy[index]);
                            EE.Electron_charge->push_back(ntuple->Electron_charge[index]);
                            EE.Electron_gsfpT->push_back(ntuple->Electron_gsfpT[index]);
                            EE.Electron_gsfPx->push_back(ntuple->Electron_gsfPx[index]);
                            EE.Electron_gsfPy->push_back(ntuple->Electron_gsfPy[index]);
                            EE.Electron_gsfPz->push_back(ntuple->Electron_gsfPz[index]);
                            EE.Electron_gsfEta->push_back(ntuple->Electron_gsfEta[index]);
                            EE.Electron_gsfPhi->push_back(ntuple->Electron_gsfPhi[index]);
                            EE.Electron_gsfCharge->push_back(ntuple->Electron_gsfCharge[index]);
                            EE.Electron_etaSC->push_back(ntuple->Electron_etaSC[index]);
                            EE.Electron_phiSC->push_back(ntuple->Electron_phiSC[index]);
                            EE.Electron_etaWidth->push_back(ntuple->Electron_etaWidth[index]);
                            EE.Electron_phiWidth->push_back(ntuple->Electron_phiWidth[index]);
                            EE.Electron_dEtaIn->push_back(ntuple->Electron_dEtaIn[index]);
                            EE.Electron_dEtaInSeed->push_back(ntuple->Electron_dEtaInSeed[index]);
                            EE.Electron_dPhiIn->push_back(ntuple->Electron_dPhiIn[index]);
                            EE.Electron_sigmaIEtaIEta->push_back(ntuple->Electron_sigmaIEtaIEta[index]);
                            EE.Electron_Full5x5_SigmaIEtaIEta->push_back(ntuple->Electron_Full5x5_SigmaIEtaIEta[index]);
                            EE.Electron_HoverE->push_back(ntuple->Electron_HoverE[index]);
                            EE.Electron_fbrem->push_back(ntuple->Electron_fbrem[index]);
                            EE.Electron_eOverP->push_back(ntuple->Electron_eOverP[index]);
                            EE.Electron_InvEminusInvP->push_back(ntuple->Electron_InvEminusInvP[index]);
                            EE.Electron_dxyVTX->push_back(ntuple->Electron_dxyVTX[index]);
                            EE.Electron_dzVTX->push_back(ntuple->Electron_dzVTX[index]);
                            EE.Electron_dxy->push_back(ntuple->Electron_dxy[index]);
                            EE.Electron_dz->push_back(ntuple->Electron_dz[index]);
                            EE.Electron_dxyBS->push_back(ntuple->Electron_dxyBS[index]);
                            EE.Electron_dzBS->push_back(ntuple->Electron_dzBS[index]);
                            EE.Electron_chIso03->push_back(ntuple->Electron_chIso03[index]);
                            EE.Electron_nhIso03->push_back(ntuple->Electron_nhIso03[index]);
                            EE.Electron_phIso03->push_back(ntuple->Electron_phIso03[index]);
                            EE.Electron_ChIso03FromPU->push_back(ntuple->Electron_ChIso03FromPU[index]);
                            EE.Electron_mHits->push_back(ntuple->Electron_mHits[index]);
                            EE.Electron_EnergySC->push_back(ntuple->Electron_EnergySC[index]);
                            EE.Electron_preEnergySC->push_back(ntuple->Electron_preEnergySC[index]);
                            EE.Electron_rawEnergySC->push_back(ntuple->Electron_rawEnergySC[index]);
                            EE.Electron_etSC->push_back(ntuple->Electron_etSC[index]);
                            EE.Electron_E15->push_back(ntuple->Electron_E15[index]);
                            EE.Electron_E25->push_back(ntuple->Electron_E25[index]);
                            EE.Electron_E55->push_back(ntuple->Electron_E55[index]);
                            EE.Electron_RelPFIso_dBeta->push_back(ntuple->Electron_RelPFIso_dBeta[index]);
                            EE.Electron_RelPFIso_Rho->push_back(ntuple->Electron_RelPFIso_Rho[index]);
                            EE.Electron_r9->push_back(ntuple->Electron_r9[index]);
                            EE.Electron_ecalDriven->push_back(ntuple->Electron_ecalDriven[index]);
                            EE.Electron_passConvVeto->push_back(ntuple->Electron_passConvVeto[index]);
//                            EE.Electron_passLooseID->push_back(ntuple->Electron_passLooseID[index]);
                            EE.Electron_passMediumID->push_back(ntuple->Electron_passMediumID[index]);
//                            EE.Electron_passTightID->push_back(ntuple->Electron_passTightID[index]);
//                            EE.Electron_passMVAID_WP80->push_back(ntuple->Electron_passMVAID_WP80[index]);
//                            EE.Electron_passMVAID_WP90->push_back(ntuple->Electron_passMVAID_WP90[index]);
//                            EE.Electron_passHEEPID->push_back(ntuple->Electron_passHEEPID[index]);

                        } // End of vector filling

                    } // End of else()

                    ElectronTree->Fill();

                    EE.HLT_trigFired->clear();
                    EE.HLT_trigEta->clear();
                    EE.HLT_trigPhi->clear();
                    EE.HLT_trigName->clear();
                    EE.Electron_pT->clear();
                    EE.Electron_eta->clear();
                    EE.Electron_phi->clear();
                    EE.Electron_Energy->clear();
                    EE.Electron_charge->clear();
                    EE.Electron_gsfpT->clear();
                    EE.Electron_gsfPx->clear();
                    EE.Electron_gsfPy->clear();
                    EE.Electron_gsfPz->clear();
                    EE.Electron_gsfEta->clear();
                    EE.Electron_gsfPhi->clear();
                    EE.Electron_gsfCharge->clear();
                    EE.Electron_etaSC->clear();
                    EE.Electron_phiSC->clear();
                    EE.Electron_etaWidth->clear();
                    EE.Electron_phiWidth->clear();
                    EE.Electron_dEtaIn->clear();
                    EE.Electron_dEtaInSeed->clear();
                    EE.Electron_dPhiIn->clear();
                    EE.Electron_sigmaIEtaIEta->clear();
                    EE.Electron_Full5x5_SigmaIEtaIEta->clear();
                    EE.Electron_HoverE->clear();
                    EE.Electron_fbrem->clear();
                    EE.Electron_eOverP->clear();
                    EE.Electron_InvEminusInvP->clear();
                    EE.Electron_dxyVTX->clear();
                    EE.Electron_dzVTX->clear();
                    EE.Electron_dxy->clear();
                    EE.Electron_dz->clear();
                    EE.Electron_dxyBS->clear();
                    EE.Electron_dzBS->clear();
                    EE.Electron_chIso03->clear();
                    EE.Electron_nhIso03->clear();
                    EE.Electron_phIso03->clear();
                    EE.Electron_ChIso03FromPU->clear();
                    EE.Electron_mHits->clear();
                    EE.Electron_EnergySC->clear();
                    EE.Electron_preEnergySC->clear();
                    EE.Electron_rawEnergySC->clear();
                    EE.Electron_etSC->clear();
                    EE.Electron_E15->clear();
                    EE.Electron_E25->clear();
                    EE.Electron_E55->clear();
                    EE.Electron_RelPFIso_dBeta->clear();
                    EE.Electron_RelPFIso_Rho->clear();
                    EE.Electron_r9->clear();
                    EE.Electron_ecalDriven->clear();
                    EE.Electron_passConvVeto->clear();
//                    EE.Electron_passLooseID->clear();
                    EE.Electron_passMediumID->clear();
//                    EE.Electron_passTightID->clear();
//                    EE.Electron_passMVAID_WP80->clear();
//                    EE.Electron_passMVAID_WP90->clear();
//                    EE.Electron_passHEEPID->clear();

//                    Double_t reco_Pt = (ele1.Momentum + ele2.Momentum).Pt();
//                    Double_t reco_rapi = (ele1.Momentum + ele2.Momentum).Rapidity();

                } // End of event selection

            } //End of if( isTriggered )

            bar.Draw(i);
        } //End of event iteration
        cout << "\t" << timesPassed << " events have passed the event selection." << endl;

        printf("\tTotal sum of weights: %.1lf\n", SumWeight);
        printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
        if ( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);

        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

    } //end of i_tup iteration

    ElectronFile->cd();
    cout << "Writing into file...";
    ElectronTree->Write();
    cout  << "Finished." << endl << "Closing a file..." << endl;
    ElectronFile->Close();
    if ( !ElectronFile->IsOpen() ) cout << "File " << OutputName << " has been closed successfully." << endl;
    else cout << "FILE " << OutputName << " COULD NOT BE CLOSED!" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}
