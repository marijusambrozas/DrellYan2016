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
#include <fstream>
#include <TH1.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>


// -- Macro for making new data files with only selection-passing events  -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.cc"
#include "./header/FileMgr.h"
#include "./etc/RoccoR/RoccoR.cc"

void CheckGammaJetsNormalization (Bool_t Debug = kFALSE)
{   
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;   

    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile* ElectronFile = TFile::Open("~/DrellYan2016/GammaJetsTest.root", "RECREATE");
    ElectronFile->cd();

    TH1D *h_g_pT = new TH1D("h_g_pT", "Gamma pT", 500, 0, 5000);

    // Loop for all processes
    for (Process_t pr=_GJets_20to100; pr<=_GJets_2000to5000; pr=next(pr))
    {
        Mgr.SetProc(pr, kTRUE);
        cout << "===========================================================" << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "BaseLocation: " << Mgr.BaseLocation << endl << endl;

        Int_t Ntup = Mgr.FullLocation.size();

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            Mgr.SetupChain(i_tup, chain);

            NtupleHandle *ntuple = new NtupleHandle(chain);
            ntuple->TurnOnBranches_GenOthers(); // for photons

            Int_t NEvents = chain->GetEntries();
            if (Debug == kTRUE) NEvents = 1000; // using few events for debugging

            cout << "\t[Total Events: " << NEvents << "]" << endl;
            myProgressBar_t bar(NEvents);

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                ntuple->GetEvent(i);

                // -- Positive/Negative Gen-weights -- //
                Double_t gen_weight = 0;
                ntuple->GENEvt_weight < 0 ? gen_weight = -1 : gen_weight = 1;

                for (Int_t i_gen=0; i_gen<ntuple->nGenOthers; i_gen++)
                {
                    if (ntuple->GenOthers_ID[i_gen] == 22 && ntuple->GenLepton_isHardProcess[i_gen])
                    {
                        h_g_pT->Fill(ntuple->GenOthers_pT[i_gen], gen_weight*Mgr.Xsec[i_tup]*Lumi/Mgr.Wsum[i_tup]);
                    }
                }

                if (!Debug) bar.Draw(i);
            } // End of event iteration

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        } // End of i_tup iteration

        // Writing
        cout << "Writing into file...";
        ElectronFile->cd();
        Int_t write;
        write = h_g_pT->Write();

        if (write)
        {
            cout << " Histogram writing finished." << endl << "Closing a file..." << endl;
            ElectronFile->Close();
            if (!ElectronFile->IsOpen()) cout << "File GammaJetsTest.root has been closed successfully.\n" << endl;
            else cout << "FILE GammaJetsTest.root COULD NOT BE CLOSED!\n" << endl;
        }
        else
        {
            cout << " Writing was NOT successful!\n" << endl;
            ElectronFile->Close();
        }
        cout << "===========================================================\n" << endl;

    } // End of i_proc iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of CheckGammaJetsNormalization
