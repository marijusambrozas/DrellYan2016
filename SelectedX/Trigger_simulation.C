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
#include <TCanvas.h>
#include <TLegend.h>

// -- Macro for simulating trigger prescale  -- //

void Trigger_simulation(Int_t type=1)
{
    TF1 *f_rand = new TF1("f_rand", "31949*x^(-4)", 25, 3000);
    TH1D *h_pT = new TH1D("h_pT", "Simulated p_T (real)", 478, 22, 500);
    TH1D *h_pT_pres = new TH1D("h_pT_pres", "Simulated p_T (prescaled)", 478, 22, 500);
    TH1D *h_pT_corr = new TH1D("h_pT_corr", "Simulated p_T (corrected)", 478, 22, 500);
    TH1D *h_pT_corr_alt = new TH1D("h_pT_corr_alt", "Simulated p_T (alternatively corrected)", 478, 22, 500);

    Double_t thresholds[8] = {22, 30, 36, 50, 75, 90, 120, 175};
    Int_t prescales[8] = {20000, 5000, 2500, 1000, 300, 150, 70, 1};
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
                if (counters[i_tr]%prescales[i_tr] == 0)
                {
                    isTriggered = 1;
                    nTrig ++;
                    triggers[i_tr] = 1;
                }
                counters[i_tr] += 1;
            }
            else break;
        }

        Int_t f_22=0, f_30=0, f_36=0, f_50=0, f_75=0, f_90=0, f_120=0, f_175=0;
        if (isTriggered)
        {
            h_pT_pres->Fill(p_T);
            for (Int_t i_tr=0; i_tr<8; i_tr++)
            {
                if (p_T > thresholds[i_tr])
                    weight += + 1./((double)prescales[i_tr]);
            }
            if (p_T > 175) weight = 1;
            h_pT_corr->Fill(p_T, 1/weight);

            if (triggers[0] == 1 && p_T > 22 && p_T <= 30)
                f_22 = 1;
            if (triggers[1] == 1 && p_T > 30 && p_T <= 36)
                f_30 = 1;
            if (triggers[2] == 1 && p_T > 36 && p_T <= 50)
                f_36 = 1;
            if (triggers[3] == 1 && p_T > 50 && p_T <= 75)
                f_50 = 1;
            if (triggers[4] == 1 && p_T > 75 && p_T < 90)
                f_75 = 1;
            if (triggers[5] == 1 && p_T > 90 && p_T < 120)
                f_90 = 1;
            if (triggers[6] == 1 && p_T > 120 && p_T < 175)
                f_120 = 1;
            if (triggers[7] == 1 && p_T > 175)
                f_175 = 1;

            if (f_22)
                h_pT_corr_alt->Fill(p_T, prescales[0]);
            else if (f_30)
                h_pT_corr_alt->Fill(p_T, prescales[1]);
            else if (f_36)
                h_pT_corr_alt->Fill(p_T, prescales[2]);
            else if (f_50)
                h_pT_corr_alt->Fill(p_T, prescales[3]);
            else if (f_75)
                h_pT_corr_alt->Fill(p_T, prescales[4]);
            else if (f_90)
                h_pT_corr_alt->Fill(p_T, prescales[5]);
            else if (f_120)
                h_pT_corr_alt->Fill(p_T, prescales[6]);
            else if (f_175)
                h_pT_corr_alt->Fill(p_T, prescales[7]);
        }
    } // End of for(event)
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
    h_pT_corr_alt->SetLineColor(kOrange);
    h_pT->Draw("hist");
    h_pT_corr->Draw("samehist");
    h_pT_corr_alt->Draw("samehist");
    h_pT_pres->Draw("samehist");
    h_pT_corr->GetYaxis()->SetRangeUser(20, 2e7);
    c_pT->SetLogy();
    TLegend *legend = new TLegend(0.4, 0.7, 0.9, 0.9);
    legend->AddEntry(h_pT, "Simulated p_{T} (true)", "l");
    legend->AddEntry(h_pT_pres, "Simulated p_{T} (prescaled)", "l");
    legend->AddEntry(h_pT_corr, "Simulated p_{T} (corrected)", "l");
    legend->AddEntry(h_pT_corr_alt, "Simulated p_{T} (corrected alt)", "l");
    legend->Draw();
    c_pT->Update();

} // End of Trigger_simulation()
