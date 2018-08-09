#include <TChain.h>
#include <TTree.h>
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

// -- Customized Analyzer for MuMu selection -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "header/myProgressBar_t.h"

#define M_Mu 0.1056583715 // -- GeV -- //

// -- Muon Channel -- //
void MakeSelectedMuMu(Int_t type = -1, TString HLTname = "IsoMu24_OR_IsoTkMu24")
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
            Bool_t GenFlag_top = kTRUE;
//            Bool_t GenFlag_top = kFALSE;
//            vector<GenOthers> GenTopCollection;
//            GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

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

        printf("\tTotal sum of weights: %.1lf\n", SumWeight);
        printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
        printf("\tSum of unchanged (to 1 or -1) weights: %.1lf", SumWeightRaw);
        if ( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", L*Xsec[i_tup]/nEvents[i_tup]);

        Double_t LoopRunTime = looptime.CpuTime();
        cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

    } //end of i_tup iteration

    MuonFile->cd();
    cout << "Writing into file...";
    MuonTree->Write();
    cout << "Finished." << endl << "Closing a file..." << endl;
    MuonFile->Close();
    if ( !MuonFile->IsOpen() ) cout << "File SelectedMuMu_" << Type << ".root has been closed successfully." << endl;
    else cout << "FILE SelectedMuMu_" << Type << ".root COULD NOT BE CLOSED!" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}
