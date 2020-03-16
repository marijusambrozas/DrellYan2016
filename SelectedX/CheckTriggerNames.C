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

void CheckTriggerNames ()
{
    FileMgr Mgr;
    vector<Process_t> Processes = Mgr.FindProc("SinglePhoton_C");

    vector<TString> TrigNames;

    // Loop for all processes
    for (Int_t i_proc=0; i_proc<Nproc; i_proc++)
    {
        Mgr.SetProc(Processes[i_proc]);
        Int_t Ntup = Mgr.FullLocation.size();

        // Loop for all samples in a process
        for (Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TChain *chain = new TChain(Mgr.TreeName[i_tup]);
            Mgr.SetupChain(i_tup, chain);

            NtupleHandle *ntuple = new NtupleHandle(chain);
            ntuple->TurnOnBranches_Electron();

            Int_t NEvents = 1000000;

            myProgressBar_t bar(NEvents);

            // Loop for all events in the chain
            for (Int_t i=0; i<NEvents; i++)
            {
                ntuple->GetEvent(i);

                for( Int_t k = 0; k < ntuple->HLT_ntrig; k++ )
                {

                    Int_t write = 1;
                    for (Int_t j=0; j<TrigNames.size(); j++)
                    {
                        if (HLT_trigName->at((unsigned int)k) == TrigNames[j]) write=0;
                    }
                    if (write == 1) TrigNames.push_back(HLT_trigName->at((unsigned int)k));
                }

                bar.Draw(i);
            } // End of event iteration

        } // End of i_tup iteration

    } // End of i_proc iteration
    cout << "Found trigger names: " << endl;
    for (Int_t j=0; j<TrigNames.size(); j++)
        cout << TrigNames[j] << endl;
    cout << "----------------------------------------------" << endl;

} // End of CheckTriggerNames()
