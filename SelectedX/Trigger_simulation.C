#include <TROOT.h>
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
#include <TF1.h>

// -- Macro for simulating trigger prescale  -- //

void Trigger_simulation()
{
    TF1 *f_rand = new TF1("f_rand", "31949*x^(-4)", 22, 200);
    TH1D *h_pT = new TH1D("h_pT", "Simulated p_T", 178, 22, 200);
    TH1D *h_trig = new TH1D("h_trig", "Fired triggers in single event", 8, 0-0.5, 8-0.5);
    TH1D *h_2trig = new TH1D("h_2trig", "What fired: Photon_50, Photon_75 or both?", 3, 0-0.5, 3-0.5);

    Double_t thresholds[8] = {22, 30, 36, 50, 75, 90, 120, 175};
    Double_t prescales[8] = {0.0016/36.47, 0.0066/36.47, 0.0132/36.47, 0.0264/36.47, 0.13/36.47, 0.26/36.47, 0.54/36.47, 1};
    Double_t counters[8] = {99999, 99999, 99999, 99999, 99999, 99999, 99999, 99999};
    for (Int_t i=0; i<1e8; i++)
    {
        Double_t p_T = f_rand->GetRandom();
        Int_t isTriggered = 0;
        Int_t triggers[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        Int_t nTrig = 0;
        for (Int_t i_tr=0; i_tr<8; i_tr++)
        {
            if (p_T > thresholds[i_tr])
            {
                if (counters[i_tr] >= 1/prescales[i_tr])
                {
                    isTriggered = 1;
                    nTrig ++;
                    triggers[i_tr] = 1;
                    counters[i_tr]= 0;
                }
                counters[i_tr] += 1;
            }
        }
        h_trig->Fill(nTrig);
        if (triggers[3] == 1 && triggers[4] == 1) h_2trig->Fill(2); // Photon_50 AND Photon_75
        else if (triggers[3] == 1) h_2trig->Fill(0); // only Photon_50
        else if (triggers[4] == 1) h_2trig->Fill(1); // only Photon_75
        if (isTriggered) h_pT->Fill(p_T);
    }

    TCanvas *c_pT = new TCanvas("c_pT", "p_T", 800, 800);
    h_pT->SetStats(0);
    h_pT->Draw();
    c_pT->Update();

    TCanvas *c_trig = new TCanvas("c_trig", "Fired triggers in single event", 800, 800);
    h_trig->SetStats(0);
    h_trig->Draw();
    c_trig->Update();

    TCanvas *c_2trig = new TCanvas("c_2trig", "Fired triggers (out of 2)", 800, 800);
    h_2trig->SetStats(0);
    h_2trig->Draw();
    c_2trig->Update();

} // End of Trigger_simulation()
