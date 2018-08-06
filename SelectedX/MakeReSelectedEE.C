#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TAttMarker.h>
#include <TROOT.h>
#include <TApplication.h>
#include <vector>
#include <TMath.h>
#include <TFormula.h>
#include <iostream>
#include <sstream>

// -- Customized Analyzer for MuMu reselection from selected MuMu file -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "header/myProgressBar_t.h"

#define M_Mu 0.1056583715 // -- GeV -- //

// -- Electron Channel -- //
//void MakeReSelectedEE(Int_t type, TString HLTname = "IsoMu24_OR_IsoTkMu24")
void MakeReSelectedEE(TString HLTname = "Ele23Ele12")
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;
    TString OutputName = "ReSelectedEE_ZToEE_M4500to6000_2.root";

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
//        ntupleDirectory.push_back( "" ); Tag.push_back( "Data" ); // -- It will be filled later -- //
//    }
//    else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

    // -- Creating SelectedMuMu variables to assign branches -- //
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
//    const Int_t Ntup = ntupleDirectory.size();
    const Int_t Ntup = 1;
    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

//        cout << "\t<" << [i_tup] << ">" << endl;

        TChain *chain = new TChain("DYTree");
        chain->Add("/media/sf_DATA/LongSelectedEE_ZToEE_M4500to6000_2.root");

//        //Set MC chain
//        if( isMC == kTRUE ) chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
//        //Set Data chain
//        else {
//            chain->Add(BaseLocation+"/"+DataLocation+"/*.root");
//            if(type==7) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
//        }       

        LongSelectedEE_t *EEInput = new LongSelectedEE_t();
        EEInput->CreateFromChain(chain);

        Double_t SumWeight = 0, SumWeight_Separated = 0;

        Int_t NEvents = chain->GetEntries();
//        Int_t NEvents = 10000;					// test using small events
        cout << "\t[Total Events: " << NEvents << "]" << endl;
        myProgressBar_t bar(NEvents);

        if ( EEInput->File_Given == kFALSE ) cout << "Error: the chain was not assigned." << endl;
        else
        {
            for(Int_t i=0; i<NEvents; i++)
            {
                EEInput->GetEvent(i);

                // -- Positive/Negative Gen-weights -- //
                EEInput->GENEvt_weight < 0 ? EE.GENEvt_weight = -1 : EE.GENEvt_weight = 1;
                SumWeight += EE.GENEvt_weight;

                // -- Normalization -- //
                Double_t TotWeight = EE.GENEvt_weight;
                if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*EE.GENEvt_weight;

                Bool_t TriggerFlag = kFALSE;
                TriggerFlag = EEInput->isTriggered( analyzer->HLT );

                if( TriggerFlag == kTRUE && EEInput->isHardProcess == kTRUE )
                {

                    EE.isHardProcess = kTRUE;
                    // -- Reco level selection -- //
                    vector< Electron > ElectronCollection;
                    Int_t NLeptons = EEInput->Electron_charge->size();
                    for(Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                    {
                        Electron ele;
                        ele.FillFromSelectedX(EEInput, i_reco);
                        ElectronCollection.push_back( ele );
                    }

                    // -- Event Selection -- //
                    vector< Electron > SelectedElectronCollection;
                    vector< Int_t > Sel_Index; // Ntuple indexes of electrons that passed the selection
                    Bool_t isPassEventSelection = kFALSE;
                    isPassEventSelection = analyzer->EventSelection_ElectronChannel(ElectronCollection, EEInput, &SelectedElectronCollection, &Sel_Index);

                    if( isPassEventSelection == kTRUE )
                    {
                        Electron ele1 = SelectedElectronCollection[0];
                        Electron ele2 = SelectedElectronCollection[1];

                        EE.nVertices = EEInput->nVertices;
                        EE.runNum = EEInput->runNum;
                        EE.lumiBlock = EEInput->lumiBlock;
                        EE.evtNum = EEInput->evtNum;
                        EE.nPileUp = EEInput->nPileUp;
                        EE.HLT_ntrig = EEInput->HLT_ntrig;

                        EE.Electron_InvM = (ele1.Momentum + ele2.Momentum).M();

                        Int_t zero_count = 0; // resolving if there is no more information in arrays

                        for (UInt_t iter=0; iter<EEInput->HLT_trigEta->size(); iter++)
                        {
                                EE.HLT_trigFired->push_back(EEInput->HLT_trigFired->at(iter));
                                EE.HLT_trigEta->push_back(EEInput->HLT_trigEta->at(iter));
                                EE.HLT_trigPhi->push_back(EEInput->HLT_trigPhi->at(iter));
                                if(((Int_t)(iter))<EEInput->HLT_ntrig) EE.HLT_trigName->push_back(EEInput->HLT_trigName->at(iter));
                        }

                        if(Sel_Index.size()!=2) cout << "======== ERROR: The number of electrons saved is not 2 ========" << endl;
                        else
                        {
                            for (UInt_t iter=0; iter<Sel_Index.size(); iter++)
                            {
                                Int_t index = Sel_Index[iter];

                                EE.Electron_pT->push_back(EEInput->Electron_pT->at(index));
                                EE.Electron_eta->push_back(EEInput->Electron_eta->at(index));
                                EE.Electron_phi->push_back(EEInput->Electron_phi->at(index));
                                EE.Electron_Energy->push_back(EEInput->Electron_Energy->at(index));
                                EE.Electron_charge->push_back(EEInput->Electron_charge->at(index));
                                EE.Electron_gsfpT->push_back(EEInput->Electron_gsfpT->at(index));
                                EE.Electron_gsfPx->push_back(EEInput->Electron_gsfPx->at(index));
                                EE.Electron_gsfPy->push_back(EEInput->Electron_gsfPy->at(index));
                                EE.Electron_gsfPz->push_back(EEInput->Electron_gsfPz->at(index));
                                EE.Electron_gsfEta->push_back(EEInput->Electron_gsfEta->at(index));
                                EE.Electron_gsfPhi->push_back(EEInput->Electron_gsfPhi->at(index));
                                EE.Electron_gsfCharge->push_back(EEInput->Electron_gsfCharge->at(index));
                                EE.Electron_etaSC->push_back(EEInput->Electron_etaSC->at(index));
                                EE.Electron_phiSC->push_back(EEInput->Electron_phiSC->at(index));
                                EE.Electron_etaWidth->push_back(EEInput->Electron_etaWidth->at(index));
                                EE.Electron_phiWidth->push_back(EEInput->Electron_phiWidth->at(index));
                                EE.Electron_dEtaIn->push_back(EEInput->Electron_dEtaIn->at(index));
                                EE.Electron_dEtaInSeed->push_back(EEInput->Electron_dEtaInSeed->at(index));
                                EE.Electron_dPhiIn->push_back(EEInput->Electron_dPhiIn->at(index));
                                EE.Electron_sigmaIEtaIEta->push_back(EEInput->Electron_sigmaIEtaIEta->at(index));
                                EE.Electron_Full5x5_SigmaIEtaIEta->push_back(EEInput->Electron_Full5x5_SigmaIEtaIEta->at(index));
                                EE.Electron_HoverE->push_back(EEInput->Electron_HoverE->at(index));
                                EE.Electron_fbrem->push_back(EEInput->Electron_fbrem->at(index));
                                EE.Electron_eOverP->push_back(EEInput->Electron_eOverP->at(index));
                                EE.Electron_InvEminusInvP->push_back(EEInput->Electron_InvEminusInvP->at(index));
                                EE.Electron_dxyVTX->push_back(EEInput->Electron_dxyVTX->at(index));
                                EE.Electron_dzVTX->push_back(EEInput->Electron_dzVTX->at(index));
                                EE.Electron_dxy->push_back(EEInput->Electron_dxy->at(index));
                                EE.Electron_dz->push_back(EEInput->Electron_dz->at(index));
                                EE.Electron_dxyBS->push_back(EEInput->Electron_dxyBS->at(index));
                                EE.Electron_dzBS->push_back(EEInput->Electron_dzBS->at(index));
                                EE.Electron_chIso03->push_back(EEInput->Electron_chIso03->at(index));
                                EE.Electron_nhIso03->push_back(EEInput->Electron_nhIso03->at(index));
                                EE.Electron_phIso03->push_back(EEInput->Electron_phIso03->at(index));
                                EE.Electron_ChIso03FromPU->push_back(EEInput->Electron_ChIso03FromPU->at(index));
                                EE.Electron_mHits->push_back(EEInput->Electron_mHits->at(index));
                                EE.Electron_EnergySC->push_back(EEInput->Electron_EnergySC->at(index));
                                EE.Electron_preEnergySC->push_back(EEInput->Electron_preEnergySC->at(index));
                                EE.Electron_rawEnergySC->push_back(EEInput->Electron_rawEnergySC->at(index));
                                EE.Electron_etSC->push_back(EEInput->Electron_etSC->at(index));
                                EE.Electron_E15->push_back(EEInput->Electron_E15->at(index));
                                EE.Electron_E25->push_back(EEInput->Electron_E25->at(index));
                                EE.Electron_E55->push_back(EEInput->Electron_E55->at(index));
                                EE.Electron_RelPFIso_dBeta->push_back(EEInput->Electron_RelPFIso_dBeta->at(index));
                                EE.Electron_RelPFIso_Rho->push_back(EEInput->Electron_RelPFIso_Rho->at(index));
                                EE.Electron_r9->push_back(EEInput->Electron_r9->at(index));
                                EE.Electron_ecalDriven->push_back(EEInput->Electron_ecalDriven->at(index));
                                EE.Electron_passConvVeto->push_back(EEInput->Electron_passConvVeto->at(index));
    //                            EE.Electron_passLooseID->push_back(EEInput->Electron_passLooseID->at(index));
                                EE.Electron_passMediumID->push_back(EEInput->Electron_passMediumID->at(index));
    //                            EE.Electron_passTightID->push_back(EEInput->Electron_passTightID->at(index));
    //                            EE.Electron_passMVAID_WP80->push_back(EEInput->Electron_passMVAID_WP80->at(index));
    //                            EE.Electron_passMVAID_WP90->push_back(EEInput->Electron_passMVAID_WP90->at(index));
    //                            EE.Electron_passHEEPID->push_back(EEInput->Electron_passHEEPID->at(index));

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
        }

        printf("\tTotal sum of weights: %.1lf\n", SumWeight);
        printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
        if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);
        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds" << endl;
        cout << "============================================================\n" << endl;

        cout << "Checking number of entries in output and input files:   ";
        if(chain->GetEntries()==ElectronTree->GetEntries()) cout << "OK! (" << ElectronTree->GetEntries() << ")" << endl;
        else cout << "Input: " << chain->GetEntries() << "     Output: " << ElectronTree->GetEntries() << "\nCHECK SELECTION ALGORITHMS.\n" << endl;

    } //end of i_tup iteration

    ElectronFile->cd();
    cout << "Writing into file...";
    ElectronTree->Write();
    cout << "Finished." << endl << "Closing a file..." << endl;
    ElectronFile->Close();
    if (!ElectronFile->IsOpen()) cout << "File " << OutputName << " has been closed successfully." << endl;
    else cout << "FILE " << OutputName << " COULD NOT BE CLOSED!" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}
