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

void Trigger_simulation(Int_t type=1)
{
    TF1 *f_rand = new TF1("f_rand", "31949*x^(-4)", 28, 3000);
    TH1D *h_pT = new TH1D("h_pT", "Simulated p_T (real)", 478, 22, 500);
    TH1D *h_pT_pres = new TH1D("h_pT_pres", "Simulated p_T (prescaled)", 478, 22, 500);
    TH1D *h_pT_corr = new TH1D("h_pT_corr", "Simulated p_T (corrected)", 478, 22, 500);
    TH1D *h_trig = new TH1D("h_trig", "Fired triggers in single event", 9, 0-0.5, 9-0.5);
    TH1D *h_2trig = new TH1D("h_2trig", "What fired: Photon_50, Photon_75 or both?", 3, 0-0.5, 3-0.5);
    TH1D *h_trigcount = new TH1D("h_trigcount", "Firing counts for each trigger", 8, 0-0.5, 8-0.5);

    Double_t thresholds[8] = {22, 30, 36, 50, 75, 90, 120, 175};
    Double_t prescales[8] = {0.0016/36.47, 0.0066/36.47, 0.0132/36.47, 0.0264/36.47, 0.13/36.47, 0.26/36.47, 0.54/36.47, 1};
//    Double_t counters[8] = {99999, 99999, 99999, 99999, 99999, 99999, 99999, 99999};
    Int_t counters[8] = {1, 1, 1, 1, 1, 1, 1, 1};

    for (Int_t i=0; i<3.5e8; i++)
    {
        if (i % 7000000 == 0) cout << "+"; cout.flush();
        Double_t p_T = f_rand->GetRandom();
        h_pT->Fill(p_T);
        Double_t weight = 0;
        Int_t isTriggered = 0;
        Int_t triggers[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        Int_t nTrig = 0;
        for (Int_t i_tr=0; i_tr<8; i_tr++)
        {
            if (p_T > thresholds[i_tr])
            {
//                          if (counters[i_tr] >= 1/prescales[i_tr])
                if (counters[i_tr]%((int)(1/prescales[i_tr])) == 0)
                {
                    isTriggered = 1;
                    nTrig ++;
                    triggers[i_tr] = 1;
                    h_trigcount->Fill(i_tr);
//                    counters[i_tr]= 0;
                }
                if (type == 1) counters[i_tr] += 1;
            }
            else break;
        }
        if (type != 1)
        {
            for (Int_t i_tr=0; i_tr<8; i_tr++) counters[i_tr] += 1;
        }
        h_trig->Fill(nTrig);
        if (triggers[3] == 1 && triggers[4] == 1) h_2trig->Fill(2); // Photon_50 AND Photon_75
        else if (triggers[3] == 1) h_2trig->Fill(0); // only Photon_50
        else if (triggers[4] == 1) h_2trig->Fill(1); // only Photon_75
        if (isTriggered)
        {
            h_pT_pres->Fill(p_T);
            for (Int_t i_tr=0; i_tr<8; i_tr++)
            {
                if (p_T > thresholds[i_tr])
                    weight = weight + prescales[i_tr];
//                if (triggers[i_tr] == 1)
//                    weight = weight + prescales[i_tr];
            }
            if (p_T > 175) weight = 1;
            h_pT_corr->Fill(p_T, 1/weight);
        }
    }
    cout << endl;

    TCanvas *c_pT = new TCanvas("c_pT", "p_T", 800, 800);
    h_pT->SetStats(0);
    h_pT->SetTitle("");
    h_pT->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_pT->GetXaxis()->SetTitleOffset(1.15);
    h_pT->GetYaxis()->SetTitle("# entries");
    h_pT->GetYaxis()->SetTitleOffset(1.15);
    h_pT_pres->SetStats(0);
    h_pT_corr->SetStats(0);
    h_pT->SetLineColor(kGreen);
    h_pT_corr->SetLineColor(kRed);
    h_pT->Draw("hist");
    h_pT_corr->Draw("samehist");
    h_pT_pres->Draw("samehist");
    h_pT_corr->GetYaxis()->SetRangeUser(20, 2e7);
    c_pT->SetLogy();
    TLegend *legend = new TLegend(0.4, 0.7, 0.9, 0.9);
    legend->AddEntry(h_pT, "Simulated p_{T} (true)", "l");
    legend->AddEntry(h_pT_pres, "Simulated p_{T} (prescaled)", "l");
    legend->AddEntry(h_pT_corr, "Simulated p_{T} (corrected)", "l");
    legend->Draw();
    c_pT->Update();

    TCanvas *c_trig = new TCanvas("c_trig", "Fired triggers in single event", 800, 800);
    h_trig->SetStats(0);
    h_trig->Draw();
    c_trig->SetLogy();
    c_trig->Update();

    TCanvas *c_2trig = new TCanvas("c_2trig", "Fired triggers (out of 2)", 800, 800);
    h_2trig->SetStats(0);
    h_2trig->Draw();
    c_2trig->SetLogy();
    c_2trig->Update();

    TCanvas *c_trigcount = new TCanvas("c_trigcount", "Firing counts for each trigger", 800, 800);
    h_trigcount->SetStats(0);
    h_trigcount->Draw();
    c_trigcount->SetLogy();
    c_trigcount->Update();

    cout << "------------------------  SIMULATION  -------------------------------" << endl;
    cout << "Firing counts for each trigger:" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "22\t30\t36\t50\t75\t90\t120\t175" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    for (Int_t z=1; z<9; z++)
    {
        cout << h_trigcount->GetBinContent(z) << "\t";
    }
    cout << "\n---------------------------------------------------------------------\n" << endl;

    cout << "Fired photon trigger numbers (triggers/event):" << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    cout << "1\t\t2\t3\t4\t5\t6\t7\t8" << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    for (Int_t z=2; z<10; z++)
    {
        cout << h_trig->GetBinContent(z) << "\t";
    }
    cout << "\n-------------------------------------------------------------------------\n" << endl;

    cout << "How many events with Photon50 and Photon75 triggers?" << endl;
    cout << "HLT_Photon50\tHLT_Photon75\tBoth" << endl;
    cout << h_2trig->GetBinContent(1) << "\t\t" << h_2trig->GetBinContent(2) << "\t\t" << h_2trig->GetBinContent(3) << endl << endl;

} // End of Trigger_simulation()
