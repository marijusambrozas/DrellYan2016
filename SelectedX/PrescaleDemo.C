#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TROOT.h>
#include <vector>
#include <TMath.h>
#include <THistPainter.h>
#include <iostream>

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/myProgressBar_t.h"
#include "./header/FileMgr.h"
#include "./header/PrescaleProvider.h"

PrescaleProvider pp("etc/prescale/triggerData2016");

Double_t PrescalePhoton (Int_t trig_fired, Double_t HLT_pT, Int_t runNo, Int_t lumiSec);

void PrescaleDemo (TString type = "SinglePhoton_full", Bool_t trigPtMatching = kFALSE, Bool_t DEBUG = kFALSE)
{
    TH1D *h_HLT_pT_full = new TH1D("h_HLT_pT_full", "", 500, 0, 500);
    TH1D *h_HLT_pT_uncorr_full = new TH1D("h_HLT_pT_uncorr_full", "", 500, 0, 500);

    FileMgr Mgr;
    vector<Process_t> Processes = Mgr.FindProc(type);
    Int_t Nproc = Processes.size();
    if (!Nproc)
    {
        cout << "No processes found." << endl;
        return;
    }

    TFile *f;
    TString Dir = "./";
    TString debug = "";
    if (DEBUG) debug = "_DEBUG";

    // Loop for all datasets
    for (Int_t i_proc=0; i_proc<Nproc; i_proc++)
    {
        Mgr.SetProc(Processes[i_proc], kTRUE);
        Process_t pr = Mgr.CurrentProc;

        // -- Output ROOTFile -- //
        f = new TFile(Dir+"FR_Hist_TEST_E_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root", "RECREATE");

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        // -- Creating Histograms -- //
        TH1D* h_HLT_pT_uncorr = new TH1D("h_HLT_pT_uncorr", "h_HLT_pT_uncorr", 500, 0, 500); h_HLT_pT_uncorr->Sumw2();
        TH1D* h_HLT_pT = new TH1D("h_HLT_pT", "h_HLT_pT", 500, 0, 500); h_HLT_pT->Sumw2();

        // Setting up variables
        std::vector<int> *trig_fired = new std::vector<int>;
        std::vector<int> *trig_matched = new std::vector<int>;
        std::vector<double> *trig_pT = new std::vector<double>;
        std::vector<double> *pT = new std::vector<double>;
        Int_t runNum;
        Int_t lumiBlock;
        Double_t prescale_alt;

        // Setting up chain
        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        chain->SetBranchStatus("trig_fired", 1);
        chain->SetBranchStatus("trig_matched", 1);
        chain->SetBranchStatus("trig_pT", 1);
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("runNum", 1);
        chain->SetBranchStatus("lumiBlock", 1);
        chain->SetBranchAddress("trig_fired", &trig_fired);
        chain->SetBranchAddress("trig_matched", &trig_matched);
        chain->SetBranchAddress("trig_pT", &trig_pT);
        chain->SetBranchAddress("p_T", &pT);
        chain->SetBranchAddress("runNum", &runNum);
        chain->SetBranchAddress("lumiBlock", &lumiBlock);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Number of events: " << NEvents << "]" << endl;
        if (DEBUG) NEvents = 100;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            if (DEBUG){
                cout << "\nEvt " << i << endl;
                cout << "nTrig = " << trig_fired->size() << endl;
                cout << "Triggers:" << endl;
                for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
                {
                    cout << "Photon" << trig_fired->at(i_tr) << "  p_T: " << trig_pT->at(i_tr) << endl;
                }
            }

            std::vector<Double_t> pTs; // Filling this to avoid double counting
            for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
            {
                if (trigPtMatching)
                {
                    Double_t dpT = fabs(trig_pT->at(i_tr) - pT->at(trig_matched->at(i_tr)))/trig_pT->at(i_tr);
                    if (dpT > 0.3) continue;
                }

                prescale_alt = PrescalePhoton(trig_fired->at(i_tr), trig_pT->at(i_tr), runNum, lumiBlock);
                Int_t match_found = 0;
                if (pTs.size())
                {
                    for (UInt_t i=0; i<pTs.size(); i++)
                    {
                        if (trig_pT->at(i_tr) == pTs[i])
                            match_found = 1;
                        else if (prescale_alt > 0) pTs.push_back(trig_pT->at(i_tr));
                    }
                }
                else if (prescale_alt > 0) pTs.push_back(trig_pT->at(i_tr));
                if (!match_found && prescale_alt > 0)
                {
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_tr));
                    h_HLT_pT->Fill(trig_pT->at(i_tr), prescale_alt);
                    h_HLT_pT_uncorr_full->Fill(trig_pT->at(i_tr));
                    h_HLT_pT_full->Fill(trig_pT->at(i_tr), prescale_alt);
                }
            }// End of for(triggers)

        }// End of event loop

        f->cd();
        cout << "\tWriting into file...";

        h_HLT_pT_uncorr->Write();
        h_HLT_pT->Write();
        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _SinglePhoton_B) break;

        f->Close();
        if (!f->IsOpen()) cout << "File " << Dir+"FR_Hist_TEST_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " has been closed successfully.\n" << endl;
        else cout << "FILE " << Dir+"FR_Hist_TEST_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    // Drawing
    TLegend *legend = new TLegend (0.55, 0.8, 0.95, 0.95);
    legend->AddEntry(h_HLT_pT_full, "Prescale weights applied", "l");
    legend->AddEntry(h_HLT_pT_uncorr_full, "Unweighted", "l");
    TCanvas *c_HLT_pT_full = new TCanvas("c", "", 800, 800);
    c_HLT_pT_full->SetTopMargin(0.05);
    c_HLT_pT_full->SetRightMargin(0.05);
    c_HLT_pT_full->SetBottomMargin(0.15);
    c_HLT_pT_full->SetLeftMargin(0.15);
    h_HLT_pT_full->SetStats(0);
    h_HLT_pT_uncorr_full->SetStats(0);
    h_HLT_pT_full->SetLineColor(kRed);
    h_HLT_pT_full->GetXaxis()->SetTitle("HLT object p_{#lower[-0.2]{T}} [GeV/c]");
    h_HLT_pT_full->GetXaxis()->SetTitleSize(0.062);
    h_HLT_pT_full->GetXaxis()->SetTitleOffset(0.95);
    h_HLT_pT_full->GetXaxis()->SetLabelSize(0.048);
    h_HLT_pT_full->GetXaxis()->SetNdivisions(6);
    h_HLT_pT_full->GetYaxis()->SetTitle("Number of events");
    h_HLT_pT_full->GetYaxis()->SetTitleSize(0.05);
    h_HLT_pT_full->GetYaxis()->SetTitleOffset(1.4);
    h_HLT_pT_full->GetYaxis()->SetLabelSize(0.043);
    h_HLT_pT_full->Draw("hist");
    h_HLT_pT_uncorr_full->Draw("samehist");
    c_HLT_pT_full->SetLogy();
    c_HLT_pT_full->SetGridy();
    c_HLT_pT_full->SetGridx();
    legend->Draw();
    c_HLT_pT_full->Update();

} // End of PrescaleDemo()


Double_t PrescalePhoton (Int_t trig_fired, Double_t HLT_pT, Int_t runNo, Int_t lumiSec)
{
    Double_t prescale = 0.0;
    if (trig_fired == 22 && HLT_pT > 22 && HLT_pT < 30)
    {
        Int_t prescale_int = pp.hltPrescale("HLT_Photon22_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG18", runNo, lumiSec);
        prescale = (Double_t)prescale_int;
    }
    else if (trig_fired == 30 && HLT_pT > 30 && HLT_pT < 36)
    {
        Int_t prescale_int = pp.hltPrescale("HLT_Photon30_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        prescale = (Double_t)prescale_int;
    }
    else if (trig_fired == 36 && HLT_pT > 36 && HLT_pT < 50)
    {
        Int_t prescale_int = pp.hltPrescale("HLT_Photon36_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        prescale = (Double_t)prescale_int;
    }
    else if (trig_fired == 50 && HLT_pT > 50 && HLT_pT < 75)
    {
        Int_t prescale_int = pp.hltPrescale("HLT_Photon50_v", runNo, lumiSec); // all possible l1 seeds are unprescaled
        prescale = (Double_t)prescale_int;
    }
    else if (trig_fired == 75 && HLT_pT > 75 && HLT_pT < 90)
    {
        Int_t prescale_int = pp.hltPrescale("HLT_Photon75_v", runNo, lumiSec); // all possible l1 seeds are unprescaled
        prescale = (Double_t)prescale_int;
    }
    else if (trig_fired == 90 && HLT_pT > 90 && HLT_pT < 120)
    {
        Int_t prescale_int = pp.hltPrescale("HLT_Photon90_v", runNo, lumiSec); // all possible l1 seeds are unprescaled
        prescale = (Double_t)prescale_int;
    }
    else if (trig_fired == 120 && HLT_pT > 120 && HLT_pT < 175)
    {
        Int_t prescale_int = pp.hltPrescale("HLT_Photon120_v", runNo, lumiSec); // all possible l1 seeds are unprescaled
        prescale = (Double_t)prescale_int;
    }
    else if (trig_fired == 175 && HLT_pT > 175)
        prescale = 1.0; // unprescaled

    return prescale;
}
