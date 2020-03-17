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

#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/FileMgr.h"

void CheckDielectronTrigger ()
{
    TChain *chain = new TChain("recoTree/DYTree");
    chain->Add("~/Desktop/ntuple_skim_1.root");

    NtupleHandle *ntuple = new NtupleHandle(chain);
//    ntuple->TurnOnBranches_Electron();

    Int_t NEvents = 1000000;

    myProgressBar_t bar(NEvents);

    // Loop for all events in the chain
    for (Int_t i=0; i<NEvents; i++)
    {
        ntuple->GetEvent(i);
        Int_t brk = 0;

        for( Int_t k = 0; k < ntuple->HLT_ntrig; k++ )
        {
            if (ntuple->HLT_trigName->at((unsigned int)k) == "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")
            {
                brk = 1;
                cout << k << ": pT = " << ntuple->HLT_trigPt[k] << "   eta = " << ntuple->HLT_trigEta[k] << "   phi = " << ntuple->HLT_trigPhi[k] << endl;
            }
        }
        if (brk == 1) cout << endl;

//        bar.Draw(i);
    } // End of event iteration

} // End of CheckTriggerNames()
