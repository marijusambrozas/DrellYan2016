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
//void MakeSelectedEMu(Int_t type, TString HLTname = "IsoMu24_OR_IsoTkMu24")
void MakeSelectedEMu(TString HLTname = "IsoMu24_OR_IsoTkMu24")
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;
    TString OutputName = "SelectedEMu_WW_34.root";

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
