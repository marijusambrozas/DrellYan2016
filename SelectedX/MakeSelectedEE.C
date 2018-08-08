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
//void MakeSelectedEE(Int_t type, Int_t Num = 100, Int_t isTopPtReweighting = 0, TString HLTname = "Ele23Ele12")
void MakeSelectedEE(Int_t type = -1 , TString HLTname = "Ele23Ele12")
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TString DataType, DataLocation, DataLocation2, Type;

    if( type == -1 ) Type = "ZToEE_M4500to6000_2";

    if( type == 1 ) {
            DataType = "B";
            DataLocation = "DoubleEG_Run2016B";
    }
    if( type == 2 ) {
            DataType = "C";
            DataLocation = "DoubleEG_Run2016C";
    }
    if( type == 3 ) {
            DataType = "D";
            DataLocation = "DoubleEG_Run2016D";
    }
    if( type == 4 ) {
            DataType = "E";
            DataLocation = "DoubleEG_Run2016E";
    }
    if( type == 5 ) {
            DataType = "F";
            DataLocation = "DoubleEG_Run2016F";
    }
    if( type == 6 ) {
            DataType = "G";
            DataLocation = "DoubleEG_Run2016G";
    }
    if( type == 7 ) {
            DataType = "H";
            DataLocation = "DoubleEG_Run2016Hver2";
            DataLocation2 = "DoubleEG_Run2016Hver3";
    }

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

    //Creating a file
    TFile* ElectronFile;
    if ( type == -1 ) ElectronFile = new TFile("/media/sf_DATA/SelectedEE_"+Type+".root", "RECREATE");
    else if ( type < 10 ) ElectronFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedEE/Data/SelectedEE_"+Type+".root", "RECREATE");
    else if (type < 20 )  ElectronFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedEE/MC_signal/SelectedEE_"+Type+".root", "RECREATE");
    else ElectronFile = new TFile("/xrootd/store/user/mambroza/SelectedX_v1/SelectedEE/MC_bkg/SelectedEE_"+Type+".root", "RECREATE");

    TTree* ElectronTree = new TTree("DYTree", "DYTree");


    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    cout << "Type: " << Type << endl;
    if( type < 10 ) cout << "DataType: Run2016" << DataType << endl;

    TString BaseLocation;
    if( Type == "Data" ) BaseLocation = "/xrd/store/user/dpai/_v2p3_";
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
        Tag.push_back("ZToEE_M4500to6000"); // There is no such type mentioned in DYAnalyzer
        nEvents.push_back(74606);       // The actual number of events in file is taken
        Xsec.push_back(4.56E-07);       // Copied from ZToMuMu_M4500to6000
    }
    else if( Type == "Data" )
    {
        ntupleDirectory.push_back( "" ); Tag.push_back( "Data" ); // -- It will be filled later -- //
    }
    else analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);

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

    //Loop for all samples
    Int_t Ntup;
    if ( type == -1 ) Ntup = 1;    //Using just 1 ntuple for test
    else Ntup = ntupleDirectory.size();

    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

        cout << "\t<" << Tag[i_tup] << ">" << endl;

        TChain *chain = new TChain("recoTree/DYTree");
        if ( type == -1 )   // For testing
        {
            chain->Add("/media/sf_DATA/ZToEE_M4500to6000_2.root/recoTree/DYTree;7");
            chain->Add("/media/sf_DATA/ZToEE_M4500to6000_2.root/recoTree/DYTree;8");
//            chain->Add("/media/sf_DATA/ZToEE_M4500to6000_2.root/recoTree/DYTree");
        }
        else
        {
            //Set MC chain
            if( isMC == kTRUE ) chain->Add(BaseLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
            //Set Data chain
            else {
                    chain->Add(BaseLocation+"/"+DataLocation+"/*.root");
                    if( type == 7 ) chain->Add(BaseLocation+"/"+DataLocation2+"/*.root");
            }
        }

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
//            Bool_t GenFlag_top = kFALSE;
//            vector<GenOthers> GenTopCollection;
//            GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

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

        printf("\tTotal sum of weights: %.1lf\n", SumWeight);
        printf("\tSum of weights of Separated events: %.1lf\n", SumWeight_Separated);
        if ( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);

        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

    } //end of i_tup iteration

    ElectronFile->cd();
    cout << "Writing into file...";
    ElectronTree->Write();
    cout  << "Finished." << endl << "Closing a file..." << endl;
    ElectronFile->Close();
    if ( !ElectronFile->IsOpen() ) cout << "File SelectedEE_" << Type << ".root has been closed successfully." << endl;
    else cout << "FILE SelectedEE_" << Type << ".root COULD NOT BE CLOSED!" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}
