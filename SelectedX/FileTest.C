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

// -- Macro for checking if all the events are added to the chain  -- //
#include "./header/myProgressBar_t.h"
#include "./header/FileMgr.h"

void FileTest ( TString type = "" )
{
    FileMgr fm(_SingleMuon_B);
    TChain *ch_full = new TChain(fm.TreeName[0]);
    ch_full->Add(fm.FullLocation[0]);
    ch_full->Add(fm.FullLocation[1]);
    UInt_t EvtNo_full = ch->GetEntries();
    UInt_t EvtNo_separate = 0;

    myProgressBar_t bar(1317);
    for ( Int_t i=1; i<=1317; i++ )
    {
        stringstream ss;
        ss << i;
        TString Number = ss.str();
        TChain *ch_separate = new TChain(fm.TreeName[0]);
        if ( i < 1000 )
            ch_separate->Add(fm.BaseLocation+"SingleMuon/crab_SingleMuon_RunB/180326_143105/0000/ntuple_skim_"+Number+".root");
        else
            ch_separate->Add(fm.BaseLocation+"SingleMuon/crab_SingleMuon_RunB/180326_143105/0001/ntuple_skim_"+Number+".root");
        EvtNo_separate += ch_separate->GetEntries();
        bar.Draw(i);
    }

    cout << "Number of events in full chain: " << EvtNo_full << endl;
    cout << "Number of events when creating separate chains: " << EvtNo_separate << endl;
    if ( EvtNo_full == EvtNo_separate ) cout << "All good." << endl;
    else cout << "There are problems." << endl;

} // End of FileTest()
