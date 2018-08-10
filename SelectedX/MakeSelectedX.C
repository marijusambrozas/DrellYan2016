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

// -- Macro for making new data files with only selection-passing events  -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"

void MakeSelectedEE ( Int_t type, TString HLTname );
void MakeSelectedMuMu ( Int_t type, TString HLTname );
void MakeSelectedEMu ( Int_t type, TString HLTname );


void MakeSelectedX ( TString whichX, Int_t type = -1, TString HLTname = "DEFAULT" )
{
    TString HLT;
    Int_t Xselected = 0;
    if ( whichX.Contains("EE") || whichX.Contains("ee") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "Ele23Ele12";
        else HLT = HLTname;
        cout << "\n*******      MakeSelectedEE ( " << type << ", " << HLT << " )      *******" << endl;
        MakeSelectedEE(type, HLT);
    }
    if ( whichX.Contains("MuMu") || whichX.Contains("mumu") || whichX.Contains("MUMU") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****  MakeSelectedMuMu ( " << type << ", " << HLT << " )  *****" << endl;
        MakeSelectedMuMu(type, HLT);
    }
    if ( whichX.Contains("EMu") || whichX.Contains("emu") || whichX.Contains("Emu") || whichX.Contains("eMu") || whichX.Contains("EMU") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****   MakeSelectedEMu ( " << type << ", " << HLT << " )  *****" << endl;
        MakeSelectedEMu(type, HLT);
    }
    if ( Xselected == 0 ) cout << "Wrong arument!" << endl;
}


/// ----------------------------- Electron Channel ------------------------------ ///
//void MakeSelectedEE(Int_t type, Int_t Num = 100, Int_t isTopPtReweighting = 0, TString HLTname = "Ele23Ele12")
void MakeSelectedEE ( Int_t type, TString HLTname )
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TString DataType, Type;

    if( type == -1 ) Type = "ZToEE_M4500to6000_2";
    if( type == 1 ) { DataType = "DoubleEG_B"; }
    if( type == 2 ) { DataType = "DoubleEG_C"; }
    if( type == 3 ) { DataType = "DoubleEG_D"; }
    if( type == 4 ) { DataType = "DoubleEG_E"; }
    if( type == 5 ) { DataType = "DoubleEG_F"; }
    if( type == 6 ) { DataType = "DoubleEG_G"; }
    if( type == 7 ) { DataType = "DoubleEG_H"; }
    Bool_t isMC = kTRUE;
    if( type < 10 && type > -1 ) {Type = "Data"; isMC = kFALSE;}
    // -- Signal MC samples -- //
    if( type == 11 ) Type = "DYEE_M10to50";
    if( type == 12 ) Type = "DYEE_M50to100";
    if( type == 13 ) Type = "DYEE_M100toInf";
    // -- Background MC samples -- //
    if( type == 21 ) Type = "ttbar";
    if( type == 22 ) Type = "ttbarBackup";
    if( type == 23 ) Type = "ttbar_M700toInf";
    if( type == 31 ) Type = "DYTauTau_M10to50";
    if( type == 32 ) Type = "DYTauTau_M50toInf";
    if( type == 41 ) Type = "VVnST";
    if( type == 51 ) Type = "WJetsToLNu";
    if( type == 61 ) Type = "QCDEMEnriched";

//    //Creating a file
//    TFile* ElectronFile;
//    if ( type == -1 ) ElectronFile = new TFile("/media/sf_DATA/test/SelectedEE_"+Type+".root", "RECREATE");
//    else if ( type < 10 ) ElectronFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedEE/Data/SelectedEE_"+Type+".root", "RECREATE");
//    else if (type < 20 )  ElectronFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedEE/MC_signal/SelectedEE_"+Type+".root", "RECREATE");
//    else ElectronFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedEE/MC_bkg/SelectedEE_"+Type+".root", "RECREATE");

//    TTree* ElectronTree = new TTree("DYTree", "DYTree");


    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    cout << "Type: " << Type << endl;
    if( type < 10 && type > -1 ) cout << "DataType: Run2016" << DataType << endl;

    TString BaseLocation;
    if( Type == "Data" ) BaseLocation = "/xrootd/store/user/dpai/_v2p3_";
    else BaseLocation = "/xrootd/store/user/dpai/_v2p3_";

    if ( type == -1 ) cout << "DATA location: /media/sf_DATA" << endl;
    else cout << "DATA location: " << BaseLocation << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

    // -- Each ntuple directory & corresponding Tags -- //
    vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    if ( type == -1 )
    {
        Tag.push_back("ZToEE_M4500to6000_2"); // There is no such type mentioned in DYAnalyzer
        nEvents.push_back(74606);       // The actual number of events in file is taken
        Xsec.push_back(4.56E-07);       // Copied from ZToMuMu_M4500to6000
    }
    else if( Type == "Data" ) analyzer->SetupDataSamples(Type, DataType, &ntupleDirectory, &Tag);
    else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

//    // -- Creating LongSelectedEE variables to assign branches -- //
//    SelectedEE_t EE; EE.CreateNew();

//    ElectronTree->Branch("isSelPassed", &EE.isSelPassed);
//    ElectronTree->Branch("nPileUp", &EE.nPileUp);
//    ElectronTree->Branch("GENEvt_weight", &EE.GENEvt_weight);
//    ElectronTree->Branch("Electron_InvM", &EE.Electron_InvM);
//    ElectronTree->Branch("Electron_pT", &EE.Electron_pT);
//    ElectronTree->Branch("Electron_eta", &EE.Electron_eta);
//    ElectronTree->Branch("Electron_phi", &EE.Electron_phi);
//    ElectronTree->Branch("Electron_Energy", &EE.Electron_Energy);
//    ElectronTree->Branch("Electron_charge", &EE.Electron_charge);

    //Loop for all samples
    Int_t Ntup;
    if ( type == -1 ) Ntup = 1;    // Using just 1 ntuple for test
    else Ntup = ntupleDirectory.size();

    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

        cout << "\t<" << Tag[i_tup] << ">" << endl;

        //Creating a file
        TFile* ElectronFile;
        if ( type == -1 ) ElectronFile = new TFile("/media/sf_DATA/test/SelectedEE_"+Tag[i_tup]+".root", "RECREATE");
        else if ( type < 10 ) ElectronFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedEE/Data/SelectedEE_"+Tag[i_tup]+".root", "RECREATE");
        else if (type < 20 )  ElectronFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedEE/MC_signal/SelectedEE_"+Tag[i_tup]+".root", "RECREATE");
        else ElectronFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedEE/MC_bkg/SelectedEE_"+Tag[i_tup]+".root", "RECREATE");

        TTree* ElectronTree = new TTree("DYTree", "DYTree");
        // -- Creating LongSelectedEE variables to assign branches -- //
        SelectedEE_t EE; EE.CreateNew();

        ElectronTree->Branch("isSelPassed", &EE.isSelPassed);
        ElectronTree->Branch("nPileUp", &EE.nPileUp);
        ElectronTree->Branch("GENEvt_weight", &EE.GENEvt_weight);
        ElectronTree->Branch("Electron_InvM", &EE.Electron_InvM);
        ElectronTree->Branch("Electron_pT", &EE.Electron_pT);
        ElectronTree->Branch("Electron_eta", &EE.Electron_eta);
        ElectronTree->Branch("Electron_phi", &EE.Electron_phi);
        ElectronTree->Branch("Electron_Energy", &EE.Electron_Energy);
        ElectronTree->Branch("Electron_charge", &EE.Electron_charge);

        TChain *chain = new TChain("recoTree/DYTree");
        if ( type == -1 )   // For testing
        {
            chain->Add("/media/sf_DATA/test/ZToEE_M4500to6000_2.root/recoTree/DYTree;7");
            chain->Add("/media/sf_DATA/test/ZToEE_M4500to6000_2.root/recoTree/DYTree;8");
//            chain->Add("/media/sf_DATA/test/ZToEE_M4500to6000_2.root");
        }
        else chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");

        NtupleHandle *ntuple = new NtupleHandle( chain );
        if( isMC == kTRUE )
        {
            ntuple->TurnOnBranches_GenLepton(); // for all leptons
            ntuple->TurnOnBranches_GenOthers(); // for quarks
        }
        ntuple->TurnOnBranches_Electron();

        Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

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
            SumWeightRaw += ntuple->GENEvt_weight;

            // -- Separate DYLL samples -- //
            Bool_t GenFlag = kFALSE;
            GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

            // -- Separate ttbar samples -- //
            Bool_t GenFlag_top = kFALSE;
            vector<GenOthers> GenTopCollection;
            GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

            if( GenFlag == kTRUE && GenFlag_top == kTRUE ) SumWeight_Separated += EE.GENEvt_weight;

            // -- Normalization -- //
            Double_t TotWeight = EE.GENEvt_weight;
            if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*EE.GENEvt_weight;

            Bool_t TriggerFlag = kFALSE;
            TriggerFlag = ntuple->isTriggered( analyzer->HLT );

            if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
            {
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

                    EE.isSelPassed = kTRUE;
                    EE.nPileUp = ntuple->nPileUp;
                    EE.Electron_InvM = (ele1.Momentum + ele2.Momentum).M();

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

                        } // End of vector filling

                    } // End of else()

                    ElectronTree->Fill();

                    EE.Electron_pT->clear();
                    EE.Electron_eta->clear();
                    EE.Electron_phi->clear();
                    EE.Electron_Energy->clear();
                    EE.Electron_charge->clear();

                } // End of event selection

            } //End of if( isTriggered )

            bar.Draw(i);
        } //End of event iteration
        cout << "\t" << timesPassed << " events have passed the event selection." << endl;

        if ( isMC == kTRUE )
        {
            printf("\tTotal sum of weights: %.1lf\n", SumWeight);
            printf("\tSum of weights of Separated events: %.1lf\n", SumWeight_Separated);
            printf("\tSum of unchanged (to 1 or -1) weights: %.1lf\n", SumWeightRaw);
            printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);
        }

        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        // Writing
        ElectronFile->cd();
        cout << "Writing into file...";
        Int_t write;
        write = ElectronTree->Write();
        if ( write )
        {
            cout << " Finished." << endl << "Closing a file..." << endl;
            ElectronFile->Close();
            if ( !ElectronFile->IsOpen() ) cout << "File SelectedEE_" << Tag[i_tup] << ".root has been closed successfully." << endl;
            else cout << "FILE SelectedEE_" << Type << ".root COULD NOT BE CLOSED!" << endl;
        }
        else
        {
            cout << " Writing was NOT successful!" << endl;
            ElectronFile->Close();
        }

    } //end of i_tup iteration

//    ElectronFile->cd();
//    cout << "Writing into file...";
//    Int_t write;
//    write = ElectronTree->Write();
//    if ( write )
//    {
//        cout << " Finished." << endl << "Closing a file..." << endl;
//        ElectronFile->Close();
//        if ( !ElectronFile->IsOpen() ) cout << "File SelectedEE_" << Type << ".root has been closed successfully." << endl;
//        else cout << "FILE SelectedEE_" << Type << ".root COULD NOT BE CLOSED!" << endl;
//    }
//    else
//    {
//        cout << " Writing was NOT successful!" << endl;
//        ElectronFile->Close();
//    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

}// end of MakeSelectedEE


/// -------------------------------- Muon Channel ------------------------------------ ///
void MakeSelectedMuMu ( Int_t type, TString HLTname )
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TString DataType, Type;

    if ( type == -1 ) Type = "ZToMuMu_M4500to6000_4";
    if ( type == 1 ) DataType = "SingleMuon_B";
    if ( type == 2 ) DataType = "SingleMuon_C";
    if ( type == 3 ) DataType = "SingleMuon_D";
    if ( type == 4 ) DataType = "SingleMuon_E";
    if ( type == 5 ) DataType = "SingleMuon_F";
    if ( type == 6 ) DataType = "SingleMuon_G";
    if ( type == 7 ) DataType = "SingleMuon_H";
    Bool_t isMC = kTRUE;
    if ( type < 10 && type > -1 ) {Type = "Data"; isMC = kFALSE;}
    // -- Signal MC samples -- //
    if ( type == 11 ) Type = "DYMuMu_M10to50";
    if ( type == 12 ) Type = "DYMuMu_M50to100";
    if ( type == 13 ) Type = "DYMuMu_M100toInf";
    // -- Background MC samples -- //
    if ( type == 21 ) Type = "ttbar";
    if ( type == 22 ) Type = "ttbarBackup";
    if ( type == 23 ) Type = "ttbar_M700toInf";
    if ( type == 31 ) Type = "DYTauTau_M10to50";
    if ( type == 32 ) Type = "DYTauTau_M50toInf";
    if ( type == 41 ) Type = "VVnST";
    if ( type == 51 ) Type = "WJetsToLNu";
    if ( type == 61 ) Type = "QCDMuEnriched";

    //Creating a file
    TFile* MuonFile;
    if ( type == -1 ) MuonFile = new TFile("/media/sf_DATA/test/SelectedMuMu_"+Type+".root", "RECREATE");
    else if ( type < 10 && type > -1 ) MuonFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedMuMu/Data/SelectedMuMu_"+Type+".root", "RECREATE");
    else if (type < 20 )  MuonFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedMuMu/MC_signal/SelectedMuMu_"+Type+".root", "RECREATE");
    else MuonFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedMuMu/MC_bkg/SelectedMuMu_"+Type+".root", "RECREATE");

    TTree* MuonTree = new TTree("DYTree", "DYTree");

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    cout << "Type: " << Type << endl;
    if( type < 10 && type > -1 ) cout << "DataType: Run2016" << DataType << endl;

    TString BaseLocation;
    if( Type == "Data" ) BaseLocation = "/xrootd/store/user/dpai/_v2p3_";
    else BaseLocation = "/xrootd/store/user/dpai/_v2p3_";
    cout << "DATA location: " << BaseLocation << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

    // -- Each ntuple directory & corresponding Tags -- //
    vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    if (type == -1 )
    {
        Tag.push_back("ZMuMu_M4500to6000");
        nEvents.push_back(100000);  // Does not match the number of events inside the MC file provided
        Xsec.push_back(4.56E-07);
    }
    else if( Type == "Data" ) analyzer->SetupDataSamples(Type, DataType, &ntupleDirectory, &Tag);
    else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

    // -- Creating SelectedMuMu variables to assign branches -- //
    SelectedMuMu_t MuMu; MuMu.CreateNew();

    MuonTree->Branch("isSelPassed", &MuMu.isSelPassed);
    MuonTree->Branch("nPileUp", &MuMu.nPileUp);
    MuonTree->Branch("GENEvt_weight", &MuMu.GENEvt_weight);
    MuonTree->Branch("Muon_pT", &MuMu.Muon_pT);
    MuonTree->Branch("Muon_eta", &MuMu.Muon_eta);
    MuonTree->Branch("Muon_phi", &MuMu.Muon_phi);
    MuonTree->Branch("Muon_charge", &MuMu.Muon_charge);
    MuonTree->Branch("Muon_Energy", &MuMu.Muon_Energy);
    MuonTree->Branch("Muon_InvM", &MuMu.Muon_InvM);
    MuonTree->Branch("Muon_TuneP_pT", &MuMu.Muon_TuneP_pT);
    MuonTree->Branch("Muon_TuneP_eta", &MuMu.Muon_TuneP_eta);
    MuonTree->Branch("Muon_TuneP_phi", &MuMu.Muon_TuneP_phi);

    //Loop for all samples
    Int_t Ntup;
    if ( type == -1 ) Ntup = 1;  // Using just 1 ntuple for test
    else Ntup = ntupleDirectory.size();

    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

        cout << "\t<" << Tag[i_tup] << ">" << endl;

        TChain *chain = new TChain("recoTree/DYTree");
        if ( type == -1 )   // For testing
        {
            chain->Add("/media/sf_DATA/test/ZToMuMu_M4500to6000_4.root/recoTree/DYTree;2"); // NEED A WAY TO TELL THE NUMBER OF CYCLES AND THEIR EXTENTION NAMES
            chain->Add("/media/sf_DATA/test/ZToMuMu_M4500to6000_4.root/recoTree/DYTree;3");
//            chain->Add("/media/sf_DATA/test/ZToMuMu_M4500to6000_4.root");
        }
        else chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");

        NtupleHandle *ntuple = new NtupleHandle( chain );
        if( isMC == kTRUE )
        {
            ntuple->TurnOnBranches_GenLepton(); // for all leptons
            ntuple->TurnOnBranches_GenOthers(); // for quarks
        }
        ntuple->TurnOnBranches_Muon();

        Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

        Int_t NEvents = chain->GetEntries();
//        Int_t NEvents = 10000;                    // test using few events
        cout << "\t[Total Events: " << NEvents << "]" << endl;
        myProgressBar_t bar (NEvents);
        Int_t timesPassed = 0;

        for(Int_t i=0; i<NEvents; i++)
        {
            ntuple->GetEvent(i);

            // -- Positive/Negative Gen-weights -- //           // IS THIS NECESSARY?
            ntuple->GENEvt_weight < 0 ? MuMu.GENEvt_weight = -1 : MuMu.GENEvt_weight = 1;
            SumWeight += MuMu.GENEvt_weight;
            SumWeightRaw += ntuple->GENEvt_weight;

            // -- Separate DYLL samples -- //
            Bool_t GenFlag = kFALSE;
            GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

            // -- Separate ttbar samples -- //
            Bool_t GenFlag_top = kFALSE;
            vector<GenOthers> GenTopCollection;
            GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

            if( GenFlag == kTRUE && GenFlag_top == kTRUE ) SumWeight_Separated += MuMu.GENEvt_weight;

            // -- Normalization -- //
            Double_t TotWeight = MuMu.GENEvt_weight;
            if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu.GENEvt_weight;

            Bool_t TriggerFlag = kFALSE;
            TriggerFlag = ntuple->isTriggered( analyzer->HLT );

            if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
            {
                // -- Reco level selection -- //
                vector< Muon > MuonCollection;
                Int_t NLeptons = ntuple->nMuon;
                for(Int_t i_reco=0; i_reco<NLeptons; i_reco++)
                {
                    Muon mu;
                    mu.FillFromNtuple(ntuple, i_reco);

                    // -- Convert to TuneP variables -- //
                    analyzer->ConvertToTunePInfo( mu );
                    MuonCollection.push_back( mu );
            }

                // -- Event Selection -- //
                vector< Muon > SelectedMuonCollection;
                vector< Int_t > Sel_Index;
                Bool_t isPassEventSelection = kFALSE;
                isPassEventSelection = analyzer->EventSelection_Zdiff_13TeV_HighPt(MuonCollection, ntuple, &SelectedMuonCollection, &Sel_Index);

                if( isPassEventSelection == kTRUE )
                {
                    timesPassed++;
                    Muon mu1 = SelectedMuonCollection[0];
                    Muon mu2 = SelectedMuonCollection[1];

                    MuMu.isSelPassed = kTRUE;
                    MuMu.nPileUp = ntuple->nPileUp;
                    MuMu.Muon_InvM = (mu1.Momentum + mu2.Momentum).M();

                    if(Sel_Index.size()!=2) cout << "======== ERROR: The number of muons saved is not 2 ========" << endl;
                    else
                    {
                        for (UInt_t iter=0; iter<Sel_Index.size(); iter++)
                        {
                            Int_t index = Sel_Index[iter];

                            MuMu.Muon_pT->push_back(ntuple->Muon_pT[index]);
                            MuMu.Muon_eta->push_back(ntuple->Muon_eta[index]);
                            MuMu.Muon_phi->push_back(ntuple->Muon_phi[index]);
                            MuMu.Muon_charge->push_back(ntuple->Muon_charge[index]);

                            Double_t Mu_E = sqrt( ntuple->Muon_Px[index]*ntuple->Muon_Px[index] + ntuple->Muon_Py[index]*ntuple->Muon_Py[index]
                                                  + ntuple->Muon_Pz[index]*ntuple->Muon_Pz[index] + M_Mu*M_Mu );

                            MuMu.Muon_Energy->push_back(Mu_E);

                            MuMu.Muon_TuneP_pT->push_back(ntuple->Muon_TuneP_pT[index]);
                            MuMu.Muon_TuneP_eta->push_back(ntuple->Muon_TuneP_eta[index]);
                            MuMu.Muon_TuneP_phi->push_back(ntuple->Muon_TuneP_phi[index]);
                        } // End of vector filling

                    } // End of else()

                    MuonTree->Fill();

                    MuMu.Muon_pT->clear();
                    MuMu.Muon_eta->clear();
                    MuMu.Muon_phi->clear();
                    MuMu.Muon_charge->clear();
                    MuMu.Muon_Energy->clear();
                    MuMu.Muon_TuneP_pT->clear();;
                    MuMu.Muon_TuneP_eta->clear();
                    MuMu.Muon_TuneP_phi->clear();

                } // End of event selection

            } //End of if( isTriggered )

            bar.Draw(i);

        } //End of event iteration
        cout << "\t" << timesPassed << " events have passed the event selection." << endl;

        if ( isMC == kTRUE )
        {
            printf("\tTotal sum of weights: %.1lf\n", SumWeight);
            printf("\tSum of weights of Separated events: %.1lf\n", SumWeight_Separated);
            printf("\tSum of unchanged (to 1 or -1) weights: %.1lf\n", SumWeightRaw);
            printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);
        }

        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

    } //end of i_tup iteration

    MuonFile->cd();
    cout << "Writing into file...";
    Int_t write;
    write = MuonTree->Write();
    if ( write )
    {
        cout << " Finished." << endl << "Closing a file..." << endl;
        MuonFile->Close();
        if ( !MuonFile->IsOpen() ) cout << "File SelectedMuMu_" << Type << ".root has been closed successfully." << endl;
        else cout << "FILE SelectedMuMu_" << Type << ".root COULD NOT BE CLOSED!" << endl;
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

}// end of MakeSelectedMuMu


/// --------------------------------- EMu events --------------------------------- ///
void MakeSelectedEMu ( Int_t type, TString HLTname )
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TString DataType, Type;

    if ( type == -1 ) Type = "WW_34";
    else if ( type == 1 ) DataType = "SingleMuon_B";
    else if ( type == 2 ) DataType = "SingleMuon_C";
    else if ( type == 3 ) DataType = "SingleMuon_D";
    else if ( type == 4 ) DataType = "SingleMuon_E";
    else if ( type == 5 ) DataType = "SingleMuon_F";
    else if ( type == 6 ) DataType = "SingleMuon_G";
    else if ( type == 7 ) DataType = "SingleMuon_H";
    Bool_t isMC = kTRUE;
    if ( type < 10  && type > -1 ) { Type = "Data"; isMC = kFALSE; }
    // -- Background MC samples -- //
    else if ( type == 21 ) Type = "ttbar";
    else if ( type == 22 ) Type = "ttbarBackup";
    else if ( type == 23 ) Type = "ttbar_M700toInf";
    else if ( type == 31 ) Type = "DYTauTau_M10to50";
    else if ( type == 32 ) Type = "DYTauTau_M50toInf";
    else if ( type == 41 ) Type = "VVnST";

    //Creating a file
    TFile* EMuFile;
    if ( type == -1 ) EMuFile = new TFile("/media/sf_DATA/test/SelectedEMu_"+Type+".root", "RECREATE");
    else if ( type < 10 && type > -1 ) EMuFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedEMu/Data/SelectedEMu_"+Type+".root", "RECREATE");
    else EMuFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedEMu/MC_bkg/SelectedEMu_"+Type+".root", "RECREATE");

    TTree* EMuTree = new TTree("DYTree", "DYTree");

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    cout << "Type: " << Type << endl;
    if( type < 10  && type > -1) cout << "DataType: Run2016 " << DataType << endl;

    TString BaseLocation;
    if( Type == "Data" ) BaseLocation = "/xrootd/store/user/dpai/_v2p3_";
    else BaseLocation = "/xrootd/store/user/dpai/_v2p3_";
    cout << "DATA location: " << BaseLocation << endl;

    TStopwatch totaltime;
    totaltime.Start();

    DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

    // -- Each ntuple directory & corresponding Tags -- //
    vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    if ( type == -1 )
    {
        Tag.push_back("WW");
        nEvents.push_back(6987123);
        Xsec.push_back(118.7);
    }
    else if( Type == "Data" ) analyzer->SetupDataSamples(Type, DataType, &ntupleDirectory, &Tag);
    else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

    // -- Creating LongSelectedMuMu variables to assign branches -- //
    SelectedEMu_t EMu; EMu.CreateNew();

    EMuTree->Branch("isSelPassed", &EMu.isSelPassed);
    EMuTree->Branch("nPileUp", &EMu.nPileUp);
    EMuTree->Branch("GENEvt_weight", &EMu.GENEvt_weight);
    EMuTree->Branch("EMu_InvM", &EMu.EMu_InvM);
    EMuTree->Branch("Muon_pT", &EMu.Muon_pT);
    EMuTree->Branch("Muon_eta", &EMu.Muon_eta);
    EMuTree->Branch("Muon_phi", &EMu.Muon_phi);
    EMuTree->Branch("Muon_charge", &EMu.Muon_charge);
    EMuTree->Branch("Muon_Energy", &EMu.Muon_Energy);
    EMuTree->Branch("Muon_TuneP_pT", &EMu.Muon_TuneP_pT);
    EMuTree->Branch("Muon_TuneP_eta", &EMu.Muon_TuneP_eta);
    EMuTree->Branch("Muon_TuneP_phi", &EMu.Muon_TuneP_phi);
    EMuTree->Branch("Electron_pT", &EMu.Electron_pT);
    EMuTree->Branch("Electron_eta", &EMu.Electron_eta);
    EMuTree->Branch("Electron_phi", &EMu.Electron_phi);
    EMuTree->Branch("Electron_Energy", &EMu.Electron_Energy);
    EMuTree->Branch("Electron_charge", &EMu.Electron_charge);

    //Loop for all samples
    Int_t Ntup;
    if ( type == -1 ) Ntup = 1;       //Using just 1 ntuple for test
    else Ntup = ntupleDirectory.size();

    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

        cout << "\t<" << Tag[i_tup] << ">" << endl;

        TChain *chain = new TChain("recoTree/DYTree");
        if ( type == -1 )
        {
            chain->Add("/media/sf_DATA/test/WW_34.root/recoTree/DYTree;1"); // NEED A WAY TO TELL THE NUMBER OF CYCLES AND THEIR EXTENTION NAMES
            chain->Add("/media/sf_DATA/test/WW_34.root/recoTree/DYTree;2");
//            chain->Add("/media/sf_DATA/test/WW_34.root");
        }
        else chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");

        NtupleHandle *ntuple = new NtupleHandle( chain );
        if( isMC == kTRUE ) {
            ntuple->TurnOnBranches_GenLepton(); // for all leptons
            ntuple->TurnOnBranches_GenOthers(); // for quarks
        }
        ntuple->TurnOnBranches_Electron();
        ntuple->TurnOnBranches_Muon();

        Double_t SumWeight = 0, SumWeight_Separated = 0, SumWeightRaw = 0;

        Int_t NEvents = chain->GetEntries();
//        Int_t NEvents = 10000;					// test using small events
        cout << "\t[Total Events: " << NEvents << "]" << endl;
        myProgressBar_t bar (NEvents);
        Int_t timesPassed = 0;

        for(Int_t i=0; i<NEvents; i++)
        {
            ntuple->GetEvent(i);

            // -- Positive/Negative Gen-weights -- //
            ntuple->GENEvt_weight < 0 ? EMu.GENEvt_weight = -1 : EMu.GENEvt_weight = 1;
            SumWeight += EMu.GENEvt_weight;
            SumWeightRaw += ntuple->GENEvt_weight;

            // -- Separate DYLL samples -- //
            Bool_t GenFlag = kFALSE;
            GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

            // -- Separate ttbar samples -- //
            Bool_t GenFlag_top = kFALSE;
            vector<GenOthers> GenTopCollection;
            GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

            if( GenFlag == kTRUE && GenFlag_top == kTRUE ) SumWeight_Separated += EMu.GENEvt_weight;

            // -- Normalization -- //
            Double_t TotWeight = EMu.GENEvt_weight;
            if( isMC == kTRUE ) TotWeight = (L*Xsec[i_tup]/nEvents[i_tup])*EMu.GENEvt_weight;

            Bool_t TriggerFlag = kFALSE;
            TriggerFlag = ntuple->isTriggered( analyzer->HLT );

            if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
            {
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
                    timesPassed++;
                    Muon mu = SelectedMuonCollection[0];
                    Electron ele = SelectedElectronCollection[0];

                    EMu.isSelPassed = kTRUE;
                    EMu.nPileUp = ntuple->nPileUp;
                    EMu.EMu_InvM = (mu.Momentum + ele.Momentum).M();
                    EMu.Muon_pT = ntuple->Muon_pT[Sel_Index_Mu];
                    EMu.Muon_eta = ntuple->Muon_eta[Sel_Index_Mu];
                    EMu.Muon_phi = ntuple->Muon_phi[Sel_Index_Mu];
                    EMu.Muon_charge = ntuple->Muon_charge[Sel_Index_Mu];
                    EMu.Muon_Energy = sqrt( ntuple->Muon_Px[Sel_Index_Mu]*ntuple->Muon_Px[Sel_Index_Mu] + ntuple->Muon_Py[Sel_Index_Mu]*ntuple->Muon_Py[Sel_Index_Mu]
                                          + ntuple->Muon_Pz[Sel_Index_Mu]*ntuple->Muon_Pz[Sel_Index_Mu] + M_Mu*M_Mu );
                    EMu.Muon_TuneP_pT = ntuple->Muon_TuneP_pT[Sel_Index_Mu];
                    EMu.Muon_TuneP_eta = ntuple->Muon_TuneP_eta[Sel_Index_Mu];
                    EMu.Muon_TuneP_phi = ntuple->Muon_TuneP_phi[Sel_Index_Mu];
                    EMu.Electron_pT = ntuple->Electron_pT[Sel_Index_Ele];
                    EMu.Electron_eta = ntuple->Electron_eta[Sel_Index_Ele];
                    EMu.Electron_phi = ntuple->Electron_phi[Sel_Index_Ele];
                    EMu.Electron_Energy = ntuple->Electron_Energy[Sel_Index_Ele];
                    EMu.Electron_charge = ntuple->Electron_charge[Sel_Index_Ele];

                    EMuTree->Fill();

                } // End of event selection

            } //End of if( isTriggered )

            bar.Draw(i);
        } //End of event iteration
        cout << "\t" << timesPassed << " events have passed the event selection." << endl;

        if ( isMC == kTRUE )
        {
            printf("\tTotal sum of weights: %.1lf\n", SumWeight);
            printf("\tSum of weights of Separated events: %.1lf\n", SumWeight_Separated);
            printf("\tSum of unchanged (to 1 or -1) weights: %.1lf\n", SumWeightRaw);
            printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);
        }

        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

    } //end of i_tup iteration

    EMuFile->cd();
    cout << "Writing into file...";
    Int_t write;
    write = EMuTree->Write();
    if ( write )
    {
        cout << " Finished." << endl << "Closing a file..." << endl;
        EMuFile->Close();
        if ( !EMuFile->IsOpen() ) cout << "File SelectedEMu_" << Type << ".root has been closed successfully." << endl;
        else cout << "FILE SelectedEMu_" << Type << ".root COULD NOT BE CLOSED!" << endl;
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

}// end of MakeSelectedEMu
