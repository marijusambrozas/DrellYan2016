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

void MakeReSelectedEE ( TString HLTname );
void MakeReSelectedMuMu ( TString HLTname );
void MakeReSelectedEMu ( TString HLTname );

void MakeReSelectedX ( TString whichX, TString HLTname = "DEFAULT" )
{
    TString HLT;
    Int_t Xselected = 0;
    if ( whichX.Contains("EE") || whichX.Contains("ee") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "Ele23Ele12";
        else HLT = HLTname;
        cout << "\n*******       MakeReSelectedEE ( " << HLT << " )       *******" << endl;
        MakeReSelectedEE(HLT);
    }
    if ( whichX.Contains("MuMu") || whichX.Contains("mumu") || whichX.Contains("MUMU") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****   MakeReSelectedMuMu ( " << HLT << " )   *****" << endl;
        MakeReSelectedMuMu(HLT);
    }
    if ( whichX.Contains("EMu") || whichX.Contains("emu") || whichX.Contains("Emu") || whichX.Contains("eMu") || whichX.Contains("EMU") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****    MakeReSelectedEMu ( " << HLT << " )   *****" << endl;
        MakeReSelectedEMu(HLT);
    }
    if ( Xselected == 0 ) cout << "Wrong arument!" << endl;
}

/// ----------------------------- Electron Channel ------------------------------ ///
//void MakeReSelectedEE ( Int_t type, TString HLTname )
void MakeReSelectedEE( TString HLTname )
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
    TFile* ElectronFile = new TFile("/media/sf_DATA/test/"+OutputName, "RECREATE");
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
        chain->Add("/media/sf_DATA/test/LongSelectedEE_ZToEE_M4500to6000_2.root");

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

//                        for (UInt_t iter=0; iter<EEInput->HLT_trigEta->size(); iter++)
//                        {
//                            EE.HLT_trigFired->push_back(EEInput->HLT_trigFired->at(iter));
//                            EE.HLT_trigEta->push_back(EEInput->HLT_trigEta->at(iter));
//                            EE.HLT_trigPhi->push_back(EEInput->HLT_trigPhi->at(iter));
//                            if(((Int_t)(iter))<EEInput->HLT_ntrig) EE.HLT_trigName->push_back(EEInput->HLT_trigName->at(iter));
//                        }

                        for (Int_t iter=0; iter<EEInput->HLT_ntrig; iter++)
                        {
                            EE.HLT_trigFired->push_back(EEInput->HLT_trigFired->at(iter));
                            EE.HLT_trigEta->push_back(EEInput->HLT_trigEta->at(iter));
                            EE.HLT_trigPhi->push_back(EEInput->HLT_trigPhi->at(iter));
                            EE.HLT_trigName->push_back(EEInput->HLT_trigName->at(iter));
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

        printf("\tSum of weights: %.1lf\n", SumWeight);
        if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);
        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds" << endl;
        cout << "===========================================================" << endl;

        cout << "Checking number of entries in output and input files:   ";
        if(chain->GetEntries()==ElectronTree->GetEntries()) cout << "OK! (" << ElectronTree->GetEntries() << ")" << endl;
        else cout << "Input: " << chain->GetEntries() << "     Output: " << ElectronTree->GetEntries() << "\nCHECK SELECTION ALGORITHMS.\n" << endl;

    } //end of i_tup iteration

    ElectronFile->cd();
    cout << "Writing into file...";
    Int_t write;
    write = ElectronTree->Write();
    if ( write )
    {
        cout << "Finished." << endl << "Closing a file..." << endl;
        ElectronFile->Close();
        if (!ElectronFile->IsOpen()) cout << "File " << OutputName << " has been closed successfully." << endl;
        else cout << "FILE " << OutputName << " COULD NOT BE CLOSED!" << endl;
    }
    else
    {
        cout << " Writing was NOT successful!" << endl;
        ElectronFile->Close();
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}

/// -------------------------------- Muon Channel ------------------------------------ ///
//void MakeReSelectedMuMu ( Int_t type, TString HLTname )
void MakeReSelectedMuMu ( TString HLTname )
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;
    TString OutputName = "ReSelectedMuMu_ZToMuMu_M4500to6000_4.root";

//    TString DataType, DataLocation, DataLocation2, Type, OutputName;

//    if( type == 1 ) {
//        DataType = "B";
//        DataLocation = "SingleMuon_Run2016B";
//        OutputName = "SelectedMuMu_SingleMuon_Run2016B.root";
//    }
//    if( type == 2 ) {
//        DataType = "C";
//        DataLocation = "SingleMuon_Run2016C";
//        OutputName = "SelectedMuMu_SingleMuon_Run2016C.root";
//    }
//    if( type == 3 ) {
//        DataType = "D";
//        DataLocation = "SingleMuon_Run2016D";
//        OutputName = "SelectedMuMu_SingleMuon_Run2016D.root";
//    }
//    if( type == 4 ) {
//        DataType = "E";
//        DataLocation = "SingleMuon_Run2016E";
//        OutputName = "SelectedMuMu_SingleMuon_Run2016E.root";
//    }
//    if( type == 5 ) {
//        DataType = "F";
//        DataLocation = "SingleMuon_Run2016F";
//        OutputName = "SelectedMuMu_SingleMuon_Run2016F.root";
//    }
//    if( type == 6 ) {
//        DataType = "G";
//        DataLocation = "SingleMuon_Run2016G";
//        OutputName = "SelectedMuMu_SingleMuon_Run2016G.root";
//    }
//    if( type == 7 ) {
//        DataType = "H";
//        DataLocation = "SingleMuon_Run2016Hver2";
//        DataLocation2 = "SingleMuon_Run2016Hver3";
//        OutputName = "SelectedMuMu_SingleMuon_Run2016H.root";
//    }

    Bool_t isMC = kTRUE;
//    if( type < 10  ) {Type = "Data"; isMC = kFALSE;}
//    // -- Signal MC samples -- //
//    if( type == 11 ) {
//        Type = "DYMuMu_M10to50";
//        OutputName = "SelectedMuMu_DYMuMu_M10to50.root";
//    }
//    if( type == 12 ) {
//        Type = "DYMuMu_M50to100";
//        OutputName = "SelectedMuMu_DYMuMu_M50to100.root";
//    }
//    if( type == 13 ) {
//        Type = "DYMuMu_M100toInf";
//        OutputName = "SelectedMuMu_DYMuMu_M100toInf.root";
//    }
//    // -- Background MC samples -- //
//    if( type == 21 ) {
//        Type = "ttbar";
//        OutputName = "SelectedMuMu_ttbar.root";
//    }
//    if( type == 22 ) {
//        Type = "ttbarBackup";
//        OutputName = "SelectedMuMu_ttbarBackup.root";
//    }
//    if( type == 23 ) {
//        Type = "ttbar_M700toInf";
//        OutputName = "SelectedMuMu_ttbar_M700toInf.root";
//    }
//    if( type == 31 ) {
//        Type = "DYTauTau_M10to50";
//        OutputName = "SelectedMuMu_DYTauTau_M10to50.root";
//    }
//    if( type == 32 ) {
//        Type = "DYTauTau_M50toInf";
//        OutputName = "SelectedMuMu_DYTauTau_M50toInf.root";
//    }
//    if( type == 41 ) {
//        Type = "VVnST";
//        OutputName = "SelectedMuMu_VVnST.root";
//    }
//    if( type == 51 ) {
//        Type = "WJetsToLNu";
//        OutputName = "SelectedMuMu_WJetsToLNu.root";
//    }
//if( type == 61 ) {
//        Type = "QCDMuEnriched";
//        OutputName = "SelectedMuMu_QCDMuEnriched.root";
//    }

    //Creating a file
    TFile* MuonFile = new TFile("/media/sf_DATA/test/"+OutputName, "RECREATE");
    TTree* MuonTree = new TTree("DYTree", "DYTree");

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
    Tag.push_back("ZMuMu_M4500to6000");
    nEvents.push_back(100000);
    Xsec.push_back(4.56E-07);
//    if( Type == "Data" ) {
//        ntupleDirectory.push_back( "" ); Tag.push_back( "Data" ); // -- It will be filled later -- //
//    }
//    else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

    // -- Creating SelectedMuMu variables to assign branches -- //
    LongSelectedMuMu_t MuMu; MuMu.CreateNew();

    MuonTree->Branch("nVertices", &MuMu.nVertices);
    MuonTree->Branch("runNum", &MuMu.runNum);
    MuonTree->Branch("lumiBlock", &MuMu.lumiBlock);
    MuonTree->Branch("evtNum", &MuMu.evtNum);
    MuonTree->Branch("nPileUp", &MuMu.nPileUp);
    MuonTree->Branch("GENEvt_weight", &MuMu.GENEvt_weight);
    MuonTree->Branch("HLT_ntrig", &MuMu.HLT_ntrig);
    MuonTree->Branch("HLT_trigFired", &MuMu.HLT_trigFired);
    MuonTree->Branch("HLT_trigName", &MuMu.HLT_trigName);
//    MuonTree->Branch("HLT_trigPt", &MuMu.HLT_trigPt);
    MuonTree->Branch("HLT_trigEta", &MuMu.HLT_trigEta);
    MuonTree->Branch("HLT_trigPhi", &MuMu.HLT_trigPhi);
    MuonTree->Branch("isHardProcess", &MuMu.isHardProcess);
    MuonTree->Branch("Muon_pT", &MuMu.Muon_pT);
    MuonTree->Branch("Muon_eta", &MuMu.Muon_eta);
    MuonTree->Branch("Muon_phi", &MuMu.Muon_phi);
    MuonTree->Branch("isGLBmuon", &MuMu.isGLBmuon);
    MuonTree->Branch("isPFmuon", &MuMu.isPFmuon);
    MuonTree->Branch("isTRKmuon", &MuMu.isTRKmuon);
    MuonTree->Branch("Muon_charge", &MuMu.Muon_charge);
    MuonTree->Branch("Muon_chi2dof", &MuMu.Muon_chi2dof);
    MuonTree->Branch("Muon_muonHits", &MuMu.Muon_muonHits);
    MuonTree->Branch("Muon_nSegments", &MuMu.Muon_nSegments);
    MuonTree->Branch("Muon_nMatches", &MuMu.Muon_nMatches);
    MuonTree->Branch("Muon_trackerLayers", &MuMu.Muon_trackerLayers);
    MuonTree->Branch("Muon_pixelHits", &MuMu.Muon_pixelHits);
    MuonTree->Branch("Muon_dxyVTX", &MuMu.Muon_dxyVTX);
    MuonTree->Branch("Muon_dzVTX", &MuMu.Muon_dzVTX);
    MuonTree->Branch("Muon_trkiso", &MuMu.Muon_trkiso);
    MuonTree->Branch("Muon_PfChargedHadronIsoR04", &MuMu.Muon_PfChargedHadronIsoR04);
    MuonTree->Branch("Muon_PfNeutralHadronIsoR04", &MuMu.Muon_PfNeutralHadronIsoR04);
    MuonTree->Branch("Muon_PfGammaIsoR04", &MuMu.Muon_PfGammaIsoR04);
    MuonTree->Branch("Muon_PFSumPUIsoR04", &MuMu.Muon_PFSumPUIsoR04);
    MuonTree->Branch("Muon_Px", &MuMu.Muon_Px);
    MuonTree->Branch("Muon_Py", &MuMu.Muon_Py);
    MuonTree->Branch("Muon_Pz", &MuMu.Muon_Pz);
    MuonTree->Branch("Muon_Energy", &MuMu.Muon_Energy);
    MuonTree->Branch("Muon_InvM", &MuMu.Muon_InvM);
    MuonTree->Branch("Muon_Best_pT", &MuMu.Muon_Best_pT);
    MuonTree->Branch("Muon_Best_pTError", &MuMu.Muon_Best_pTError);
    MuonTree->Branch("Muon_Best_Px", &MuMu.Muon_Best_Px);
    MuonTree->Branch("Muon_Best_Py", &MuMu.Muon_Best_Py);
    MuonTree->Branch("Muon_Best_Pz", &MuMu.Muon_Best_Pz);
    MuonTree->Branch("Muon_Best_eta", &MuMu.Muon_Best_eta);
    MuonTree->Branch("Muon_Best_phi", &MuMu.Muon_Best_phi);
    MuonTree->Branch("Muon_Inner_pT", &MuMu.Muon_Inner_pT);
    MuonTree->Branch("Muon_Inner_pTError", &MuMu.Muon_Inner_pTError);
    MuonTree->Branch("Muon_Inner_Px", &MuMu.Muon_Inner_Px);
    MuonTree->Branch("Muon_Inner_Py", &MuMu.Muon_Inner_Py);
    MuonTree->Branch("Muon_Inner_Pz", &MuMu.Muon_Inner_Pz);
    MuonTree->Branch("Muon_Inner_eta", &MuMu.Muon_Inner_eta);
    MuonTree->Branch("Muon_Inner_phi", &MuMu.Muon_Inner_phi);
    MuonTree->Branch("Muon_Outer_pT", &MuMu.Muon_Outer_pT);
    MuonTree->Branch("Muon_Outer_pTError", &MuMu.Muon_Outer_pTError);
    MuonTree->Branch("Muon_Outer_Px", &MuMu.Muon_Outer_Px);
    MuonTree->Branch("Muon_Outer_Py", &MuMu.Muon_Outer_Py);
    MuonTree->Branch("Muon_Outer_Pz", &MuMu.Muon_Outer_Pz);
    MuonTree->Branch("Muon_Outer_eta", &MuMu.Muon_Outer_eta);
    MuonTree->Branch("Muon_Outer_phi", &MuMu.Muon_Outer_phi);
    MuonTree->Branch("Muon_GLB_pT", &MuMu.Muon_GLB_pT);
    MuonTree->Branch("Muon_GLB_pTError", &MuMu.Muon_GLB_pTError);
    MuonTree->Branch("Muon_GLB_Px", &MuMu.Muon_GLB_Px);
    MuonTree->Branch("Muon_GLB_Py", &MuMu.Muon_GLB_Py);
    MuonTree->Branch("Muon_GLB_Pz", &MuMu.Muon_GLB_Pz);
    MuonTree->Branch("Muon_GLB_eta", &MuMu.Muon_GLB_eta);
    MuonTree->Branch("Muon_GLB_phi", &MuMu.Muon_GLB_phi);
    MuonTree->Branch("Muon_TuneP_pT", &MuMu.Muon_TuneP_pT);
    MuonTree->Branch("Muon_TuneP_pTError", &MuMu.Muon_TuneP_pTError);
    MuonTree->Branch("Muon_TuneP_Px", &MuMu.Muon_TuneP_Px);
    MuonTree->Branch("Muon_TuneP_Py", &MuMu.Muon_TuneP_Py);
    MuonTree->Branch("Muon_TuneP_Pz", &MuMu.Muon_TuneP_Pz);
    MuonTree->Branch("Muon_TuneP_eta", &MuMu.Muon_TuneP_eta);
    MuonTree->Branch("Muon_TuneP_phi", &MuMu.Muon_TuneP_phi);
    MuonTree->Branch("CosAngle", &MuMu.CosAngle);
    MuonTree->Branch("vtxTrkChi2", &MuMu.vtxTrkChi2);
    MuonTree->Branch("vtxTrkProb", &MuMu.vtxTrkProb);
    MuonTree->Branch("vtxTrkNdof", &MuMu.vtxTrkNdof);
    MuonTree->Branch("vtxTrkCkt1Pt", &MuMu.vtxTrkCkt1Pt);
    MuonTree->Branch("vtxTrkCkt2Pt", &MuMu.vtxTrkCkt2Pt);

    //Loop for all samples
//    const Int_t Ntup = ntupleDirectory.size();
    const Int_t Ntup = 1;
    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

//        cout << "\t<" << [i_tup] << ">" << endl;

        TChain *chain = new TChain("DYTree");
//        //Set MC chain
//        if( isMC == kTRUE ) chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
//        //Set Data chain
//        else {
//            chain->Add(BaseLocation+"/"+DataLocation+"/*.root");
//            if(type==7) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
//        }
        chain->Add("/media/sf_DATA/test/LongSelectedMuMu_ZToMuMu_M4500to6000_4.root");

        LongSelectedMuMu_t * MuMuInput = new LongSelectedMuMu_t();
        MuMuInput->CreateFromChain(chain);

        Double_t SumWeight = 0, SumWeight_Separated = 0;

        Int_t NEvents = chain->GetEntries();
//        Int_t NEvents = 10000;					// test using small events
        cout << "\t[Total Events: " << NEvents << "]" << endl;
        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            MuMuInput->GetEvent(i);

            // -- Positive/Negative Gen-weights -- //
            MuMuInput->GENEvt_weight < 0 ? MuMu.GENEvt_weight = -1 : MuMu.GENEvt_weight = 1;
            SumWeight += MuMu.GENEvt_weight;

            // -- Normalization -- //
            Double_t TotWeight = MuMu.GENEvt_weight;
            if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu.GENEvt_weight;

            Bool_t TriggerFlag = kFALSE;
            TriggerFlag = MuMuInput->isTriggered( analyzer->HLT );

            if( TriggerFlag == kTRUE && MuMuInput->isHardProcess == kTRUE )
            {
                MuMu.isHardProcess = kTRUE;
                // -- Reco level selection -- //
                vector< Muon > MuonCollection;
                Int_t NLeptons = MuMuInput->Muon_pT->size();
                for(Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                {
                    Muon mu;
                    mu.FillFromSelectedX(MuMuInput, i_reco);

                    // -- Convert to TuneP variables -- //
                    analyzer->ConvertToTunePInfo( mu );
                    MuonCollection.push_back( mu );
                }

                // -- Event Selection -- //
                vector< Muon > SelectedMuonCollection;
                vector< Int_t > Sel_Index;
                Int_t IndexDi;
                Bool_t isPassEventSelection = kFALSE;
                isPassEventSelection = analyzer->EventSelection_Zdiff_13TeV_HighPt(MuonCollection, MuMuInput, &SelectedMuonCollection, &Sel_Index, IndexDi);

                if( isPassEventSelection == kTRUE )
                {
                    Muon mu1 = SelectedMuonCollection[0];
                    Muon mu2 = SelectedMuonCollection[1];

                    MuMu.nVertices = MuMuInput->nVertices;
                    MuMu.runNum = MuMuInput->runNum;
                    MuMu.lumiBlock = MuMuInput->lumiBlock;
                    MuMu.evtNum = MuMuInput->evtNum;
                    MuMu.nPileUp = MuMuInput->nPileUp;
                    MuMu.HLT_ntrig = MuMuInput->HLT_ntrig;

                    MuMu.Muon_InvM = (mu1.Momentum + mu2.Momentum).M();

                    if( IndexDi!=-1 )
                    {
                        MuMu.CosAngle->push_back(MuMuInput->CosAngle->at(IndexDi));
                        MuMu.vtxTrkChi2->push_back(MuMuInput->vtxTrkChi2->at(IndexDi));
                        MuMu.vtxTrkProb->push_back(MuMuInput->vtxTrkProb->at(IndexDi));
                        MuMu.vtxTrkNdof->push_back(MuMuInput->vtxTrkNdof->at(IndexDi));
                        MuMu.vtxTrkCkt1Pt->push_back(MuMuInput->vtxTrkCkt1Pt->at(IndexDi));
                        MuMu.vtxTrkCkt2Pt->push_back(MuMuInput->vtxTrkCkt2Pt->at(IndexDi));
                    }
                    else cout << "== ERROR: Event selection was passed but no dimuon index registered. ==" << endl;

//                    for (UInt_t iter=0; iter<MuMuInput->HLT_trigEta->size(); iter++)
//                    {
//                        if(MuMuInput->HLT_trigEta->at(iter) && MuMuInput->HLT_trigPhi->at(iter))
//                        {
//                            MuMu.HLT_trigFired->push_back(MuMuInput->HLT_trigFired->at(iter));
//                            MuMu.HLT_trigEta->push_back(MuMuInput->HLT_trigEta->at(iter));
//                            MuMu.HLT_trigPhi->push_back(MuMuInput->HLT_trigPhi->at(iter));
//                            if(((Int_t)(iter))<MuMuInput->HLT_ntrig) MuMu.HLT_trigName->push_back(MuMuInput->HLT_trigName->at(iter));
//                        }
//                        else
//                        {
//                            cout << "HLT_trigEta and HLT_trigPhi do not exist. Breaking at " << iter << " iteration." << endl;
//                            break;
//                        }
//                    }

                    for (Int_t iter=0; iter<MuMuInput->HLT_ntrig; iter++)
                    {
                        MuMu.HLT_trigFired->push_back(MuMuInput->HLT_trigFired->at(iter));
                        MuMu.HLT_trigEta->push_back(MuMuInput->HLT_trigEta->at(iter));
                        MuMu.HLT_trigPhi->push_back(MuMuInput->HLT_trigPhi->at(iter));
                        MuMu.HLT_trigName->push_back(MuMuInput->HLT_trigName->at(iter));
                    }

                    if(Sel_Index.size()!=2) cout << "=========== ERROR: The number of muons saved is not 2 ===========" << endl;
                    else
                    {
                        for (UInt_t iter=0; iter<Sel_Index.size(); iter++)
                        {
                            Int_t index = Sel_Index[iter];

                            MuMu.Muon_pT->push_back(MuMuInput->Muon_pT->at(index));
                            MuMu.Muon_eta->push_back(MuMuInput->Muon_eta->at(index));
                            MuMu.Muon_phi->push_back(MuMuInput->Muon_phi->at(index));
                            MuMu.isGLBmuon->push_back(MuMuInput->isGLBmuon->at(index));
                            MuMu.isPFmuon->push_back(MuMuInput->isPFmuon->at(index));
                            MuMu.isTRKmuon->push_back(MuMuInput->isTRKmuon->at(index));
                            MuMu.Muon_charge->push_back(MuMuInput->Muon_charge->at(index));
                            MuMu.Muon_chi2dof->push_back(MuMuInput->Muon_chi2dof->at(index));
                            MuMu.Muon_muonHits->push_back(MuMuInput->Muon_muonHits->at(index));
                            MuMu.Muon_nSegments->push_back(MuMuInput->Muon_nSegments->at(index));
                            MuMu.Muon_nMatches->push_back(MuMuInput->Muon_nMatches->at(index));
                            MuMu.Muon_trackerLayers->push_back(MuMuInput->Muon_trackerLayers->at(index));
                            MuMu.Muon_pixelHits->push_back(MuMuInput->Muon_pixelHits->at(index));
                            MuMu.Muon_dxyVTX->push_back(MuMuInput->Muon_dxyVTX->at(index));
                            MuMu.Muon_dzVTX->push_back(MuMuInput->Muon_dzVTX->at(index));
                            MuMu.Muon_trkiso->push_back(MuMuInput->Muon_trkiso->at(index));
                            MuMu.Muon_PfChargedHadronIsoR04->push_back(MuMuInput->Muon_PfChargedHadronIsoR04->at(index));
                            MuMu.Muon_PfNeutralHadronIsoR04->push_back(MuMuInput->Muon_PfNeutralHadronIsoR04->at(index));
                            MuMu.Muon_PfGammaIsoR04->push_back(MuMuInput->Muon_PfGammaIsoR04->at(index));
                            MuMu.Muon_PFSumPUIsoR04->push_back(MuMuInput->Muon_PFSumPUIsoR04->at(index));
                            MuMu.Muon_Px->push_back(MuMuInput->Muon_Px->at(index));
                            MuMu.Muon_Py->push_back(MuMuInput->Muon_Py->at(index));
                            MuMu.Muon_Pz->push_back(MuMuInput->Muon_Pz->at(index));
                            MuMu.Muon_Energy->push_back(MuMuInput->Muon_Energy->at(index));

                            MuMu.Muon_Best_pT->push_back(MuMuInput->Muon_Best_pT->at(index));
                            MuMu.Muon_Best_pTError->push_back(MuMuInput->Muon_Best_pTError->at(index));
                            MuMu.Muon_Best_Px->push_back(MuMuInput->Muon_Best_Px->at(index));
                            MuMu.Muon_Best_Py->push_back(MuMuInput->Muon_Best_Py->at(index));
                            MuMu.Muon_Best_Pz->push_back(MuMuInput->Muon_Best_Pz->at(index));
                            MuMu.Muon_Best_eta->push_back(MuMuInput->Muon_Best_eta->at(index));
                            MuMu.Muon_Best_phi->push_back(MuMuInput->Muon_Best_phi->at(index));
                            MuMu.Muon_Inner_pT->push_back(MuMuInput->Muon_Inner_pT->at(index));
                            MuMu.Muon_Inner_pTError->push_back(MuMuInput->Muon_Inner_pTError->at(index));
                            MuMu.Muon_Inner_Px->push_back(MuMuInput->Muon_Inner_Px->at(index));
                            MuMu.Muon_Inner_Py->push_back(MuMuInput->Muon_Inner_Py->at(index));
                            MuMu.Muon_Inner_Pz->push_back(MuMuInput->Muon_Inner_Pz->at(index));
                            MuMu.Muon_Inner_eta->push_back(MuMuInput->Muon_Inner_eta->at(index));
                            MuMu.Muon_Inner_phi->push_back(MuMuInput->Muon_Inner_phi->at(index));
                            MuMu.Muon_Outer_pT->push_back(MuMuInput->Muon_Outer_pT->at(index));
                            MuMu.Muon_Outer_pTError->push_back(MuMuInput->Muon_Outer_pTError->at(index));
                            MuMu.Muon_Outer_Px->push_back(MuMuInput->Muon_Outer_Px->at(index));
                            MuMu.Muon_Outer_Py->push_back(MuMuInput->Muon_Outer_Py->at(index));
                            MuMu.Muon_Outer_Pz->push_back(MuMuInput->Muon_Outer_Pz->at(index));
                            MuMu.Muon_Outer_eta->push_back(MuMuInput->Muon_Outer_eta->at(index));
                            MuMu.Muon_Outer_phi->push_back(MuMuInput->Muon_Outer_phi->at(index));
                            MuMu.Muon_GLB_pT->push_back(MuMuInput->Muon_GLB_pT->at(index));
                            MuMu.Muon_GLB_pTError->push_back(MuMuInput->Muon_GLB_pTError->at(index));
                            MuMu.Muon_GLB_Px->push_back(MuMuInput->Muon_GLB_Px->at(index));
                            MuMu.Muon_GLB_Py->push_back(MuMuInput->Muon_GLB_Py->at(index));
                            MuMu.Muon_GLB_Pz->push_back(MuMuInput->Muon_GLB_Pz->at(index));
                            MuMu.Muon_GLB_eta->push_back(MuMuInput->Muon_GLB_eta->at(index));
                            MuMu.Muon_GLB_phi->push_back(MuMuInput->Muon_GLB_phi->at(index));
                            MuMu.Muon_TuneP_pT->push_back(MuMuInput->Muon_TuneP_pT->at(index));
                            MuMu.Muon_TuneP_pTError->push_back(MuMuInput->Muon_TuneP_pTError->at(index));
                            MuMu.Muon_TuneP_Px->push_back(MuMuInput->Muon_TuneP_Px->at(index));
                            MuMu.Muon_TuneP_Py->push_back(MuMuInput->Muon_TuneP_Py->at(index));
                            MuMu.Muon_TuneP_Pz->push_back(MuMuInput->Muon_TuneP_Pz->at(index));
                            MuMu.Muon_TuneP_eta->push_back(MuMuInput->Muon_TuneP_eta->at(index));
                            MuMu.Muon_TuneP_phi->push_back(MuMuInput->Muon_TuneP_phi->at(index));
                        } // End of vector filling
                    } // End of else()

                    MuonTree->Fill();

                    MuMu.HLT_trigFired->clear();
                    MuMu.HLT_trigName->clear();
                    MuMu.HLT_trigEta->clear();
                    MuMu.HLT_trigPhi->clear();
                    MuMu.Muon_pT->clear();
                    MuMu.Muon_eta->clear();
                    MuMu.Muon_phi->clear();
                    MuMu.isGLBmuon->clear();
                    MuMu.isPFmuon->clear();
                    MuMu.isTRKmuon->clear();
                    MuMu.Muon_charge->clear();
                    MuMu.Muon_chi2dof->clear();
                    MuMu.Muon_muonHits->clear();
                    MuMu.Muon_nSegments->clear();
                    MuMu.Muon_nMatches->clear();
                    MuMu.Muon_trackerLayers->clear();
                    MuMu.Muon_pixelHits->clear();
                    MuMu.Muon_dxyVTX->clear();
                    MuMu.Muon_dzVTX->clear();
                    MuMu.Muon_trkiso->clear();
                    MuMu.Muon_PfChargedHadronIsoR04->clear();
                    MuMu.Muon_PfNeutralHadronIsoR04->clear();
                    MuMu.Muon_PfGammaIsoR04->clear();
                    MuMu.Muon_PFSumPUIsoR04->clear();
                    MuMu.Muon_Px->clear();
                    MuMu.Muon_Py->clear();
                    MuMu.Muon_Pz->clear();
                    MuMu.Muon_Energy->clear();
                    MuMu.Muon_Best_pT->clear();
                    MuMu.Muon_Best_pTError->clear();
                    MuMu.Muon_Best_Px->clear();
                    MuMu.Muon_Best_Py->clear();
                    MuMu.Muon_Best_Pz->clear();
                    MuMu.Muon_Best_eta->clear();
                    MuMu.Muon_Best_phi->clear();
                    MuMu.Muon_Inner_pT->clear();
                    MuMu.Muon_Inner_pTError->clear();
                    MuMu.Muon_Inner_Px->clear();
                    MuMu.Muon_Inner_Py->clear();
                    MuMu.Muon_Inner_Pz->clear();
                    MuMu.Muon_Inner_eta->clear();
                    MuMu.Muon_Inner_phi->clear();
                    MuMu.Muon_Outer_pT->clear();
                    MuMu.Muon_Outer_pTError->clear();
                    MuMu.Muon_Outer_Px->clear();
                    MuMu.Muon_Outer_Py->clear();
                    MuMu.Muon_Outer_Pz->clear();
                    MuMu.Muon_Outer_eta->clear();
                    MuMu.Muon_Outer_phi->clear();
                    MuMu.Muon_GLB_pT->clear();
                    MuMu.Muon_GLB_pTError->clear();
                    MuMu.Muon_GLB_Px->clear();
                    MuMu.Muon_GLB_Py->clear();
                    MuMu.Muon_GLB_Pz->clear();
                    MuMu.Muon_GLB_eta->clear();
                    MuMu.Muon_GLB_phi->clear();
                    MuMu.Muon_TuneP_pT->clear();;
                    MuMu.Muon_TuneP_pTError->clear();
                    MuMu.Muon_TuneP_Px->clear();
                    MuMu.Muon_TuneP_Py->clear();
                    MuMu.Muon_TuneP_Pz->clear();
                    MuMu.Muon_TuneP_eta->clear();
                    MuMu.Muon_TuneP_phi->clear();
                    MuMu.CosAngle->clear();
                    MuMu.vtxTrkChi2->clear();
                    MuMu.vtxTrkProb->clear();
                    MuMu.vtxTrkNdof->clear();
                    MuMu.vtxTrkCkt1Pt->clear();
                    MuMu.vtxTrkCkt2Pt->clear();

//                    Double_t reco_Pt = (mu1.Momentum + mu2.Momentum).Pt();
//                    Double_t reco_rapi = (mu1.Momentum + mu2.Momentum).Rapidity();

                } // End of event selection

            } //End of if( isTriggered )

            bar.Draw(i);
        } //End of event iteration

        printf("\tSum of weights: %.1lf\n", SumWeight);
        if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);
        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds" << endl;
        cout << "===========================================================" << endl;

        cout << "Checking number of entries in output and input files:   ";
        if(chain->GetEntries()==MuonTree->GetEntries()) cout << "OK! (" << MuonTree->GetEntries() << ")" << endl;
        else cout << "Input: " << chain->GetEntries() << "     Output: " << MuonTree->GetEntries() << "\nCHECK SELECTION ALGORITHMS.\n" << endl;

    } //end of i_tup iteration

    MuonFile->cd();
    cout << "Writing into file...";
    Int_t write;
    write = MuonTree->Write();
    if ( write )
    {
        cout << "Finished." << endl << "Closing a file..." << endl;
        MuonFile->Close();
        if (!MuonFile->IsOpen()) cout << "File " << OutputName << " has been closed successfully." << endl;
        else cout << "FILE " << OutputName << " COULD NOT BE CLOSED!" << endl;
    }
    else
    {
        cout << " Writing was NOT successful!" << endl;
        MuonFile->Close();
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}

/// --------------------------------- EMu events --------------------------------- ///
//void MakeReSelectedEMu ( Int_t type, TString HLTname )
void MakeReSelectedEMu ( TString HLTname )
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
    TFile* EMuFile = new TFile("/media/sf_DATA/test/"+OutputName, "RECREATE");
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
        chain->Add("/media/sf_DATA/test/LongSelectedEMu_WW_34.root");

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

        printf("\tSum of weights: %.1lf\n", SumWeight);
        if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);
        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds" << endl;
        cout << "===========================================================" << endl;

        cout << "Checking number of entries in output and input files:   ";
        if(chain->GetEntries()==EMuTree->GetEntries()) cout << "OK! (" << EMuTree->GetEntries() << ")" << endl;
        else cout << "Input: " << chain->GetEntries() << "     Output: " << EMuTree->GetEntries() << "\nCHECK SELECTION ALGORITHMS.\n" << endl;

    } //end of i_tup iteration

    EMuFile->cd();
    cout << "Writing into file...";
    Int_t write;
    write = EMuTree->Write();
    if ( write )
    {
        cout << "Finished." << endl << "Closing a file..." << endl;
        EMuFile->Close();
        if ( !EMuFile->IsOpen() ) cout << "File " << OutputName << " has been closed successfully." << endl;
        else cout << "FILE " << OutputName << " COULD NOT BE CLOSED!" << endl;
    }
    else
    {
        cout << " Writing was NOT successful!" << endl;
        EMuFile->Close();
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}
