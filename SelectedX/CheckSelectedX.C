#include <TChain.h>
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
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <sstream>

// -- Macro to check if selected and reselected events are the same -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "header/myProgressBar_t.h"

static inline Bool_t AreVectorsTheSame(vector<Double_t> *vec1, vector<Double_t> *vec2, TString name="");
static inline Bool_t AreVectorsTheSame(vector<Int_t> *vec1, vector<Int_t> *vec2, TString name="");
static inline Bool_t AreValuesTheSame(Double_t val1, Double_t val2, TString name="");
static inline Bool_t AreValuesTheSame(Int_t val1, Int_t val2, TString name="");
static inline Bool_t AreValuesTheSame(Bool_t val1, Bool_t val2, TString name);
static inline Bool_t AreHistogramsTheSame(TH1D *hist1, TH1D *hist2, TString name="", Bool_t setBreak=kFALSE);

void CheckSelectedEE(Bool_t DrawHistos, Bool_t BreakIfProblem);
void CheckSelectedMuMu(Bool_t DrawHistos, Bool_t BreakIfProblem);
void CheckSelectedEMu(Bool_t DrawHistos, Bool_t BreakIfProblem);


void CheckSelectedX ( TString WhichX, Bool_t DrawHistos = kFALSE, Bool_t BreakIfProblem = kFALSE )
{
    TString whichX = WhichX;
    whichX.ToUpper();
    Int_t Xselected = 0;
    if ( whichX.Contains("EE") )
    {
        Xselected++;
        cout << "\n**********         CheckSelectedEE         **********" << endl;
        CheckSelectedEE(DrawHistos, BreakIfProblem);
    }
    if ( whichX.Contains("MUMU") )
    {
        Xselected++;
        cout << "\n**********        CheckSelectedMuMu        **********" << endl;
        CheckSelectedMuMu(DrawHistos, BreakIfProblem);
    }
    if ( whichX.Contains("EMU") )
    {
        Xselected++;
        cout << "\n**********        CheckSelectedEMu         **********" << endl;
        CheckSelectedEMu(DrawHistos, BreakIfProblem);
    }
    if ( Xselected == 0 ) cout << "Wrong arument!" << endl;
}


/// ----------------------------- Electron Channel ------------------------------ ///
void CheckSelectedEE ( Bool_t DrawHistos, Bool_t BreakIfProblem )
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0; // Need checking
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "\n[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

    TStopwatch totaltime;
    totaltime.Start();

    // -- Each ntuple directory & corresponding Tags -- //
//    vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    Tag.push_back("ZToEE_M4500to6000");
    nEvents.push_back(74606);       // The actual number of events in file is taken
    Xsec.push_back(4.56E-07);       // Copied from ZToMuMu_M4500to6000
    const Int_t Ntup = 1;
    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

//        cout << "\t<" << [i_tup] << ">" << endl;

        TChain *ch_s = new TChain("DYTree");
        TChain *ch_l = new TChain("DYTree");
        ch_s->Add("/media/sf_DATA/test/SelectedEE_ZToEE_M4500to6000_2.root");
        ch_l->Add("/media/sf_DATA/test/LongSelectedEE_ZToEE_M4500to6000_2.root");

        SelectedEE_t *EE_s = new SelectedEE_t();
        LongSelectedEE_t *EE_l = new LongSelectedEE_t();
        EE_s->CreateFromChain(ch_s);
        EE_l->CreateFromChain(ch_l);

        TH1D* h1_pT_s = new TH1D("h1pts", "h1_pT_s", 400, 0, 4000);
        TH1D* h1_pT_l = new TH1D("h1ptl", "h1_pT_l", 400, 0, 4000);
        TH1D* h1_eta_s = new TH1D("h1etas", "h1_eta_s", 120, -3, 3);
        TH1D* h1_eta_l = new TH1D("h1etal", "h1_eta_l", 120, -3, 3);
        TH1D* h1_phi_s = new TH1D("h1phis", "h1_phi_s", 400, -4, 4);
        TH1D* h1_phi_l = new TH1D("h1phil", "h1_phi_l", 400, -4, 4);
        TH1D* h1_Energy_s = new TH1D("h1energys", "h1_Energy_s", 400, 0, 8000);
        TH1D* h1_Energy_l = new TH1D("h1energyl", "h1_Energy_l", 400, 0, 8000);
        TH1D* h1_invm_s = new TH1D("h1invms", "h1_invm_s", 350, 0, 7000);
        TH1D* h1_invm_l = new TH1D("h1invml", "h1_invm_l", 350, 0, 7000);
        TH1D* h1_nPileUp_s = new TH1D("h1npileups", "h1_nPileUp_s", 60, 0, 60);
        TH1D* h1_nPileUp_l = new TH1D("h1npileupl", "h1_nPileUp_l", 60, 0, 60);

        Bool_t AllOkVec[8];
        Bool_t AllOkHist[6];
        Bool_t AllOkV = kTRUE, AllOkH = kTRUE, AllOk = kTRUE;

        Int_t NEvents_s = ch_s->GetEntries();
        Int_t NEvents_l = ch_l->GetEntries();
        cout << "\tTotal selected Events: " << NEvents_s << endl;
        cout << "\tTotal LongSelected Events: " << NEvents_l << endl;
        if(NEvents_s!=NEvents_l)
        {
            cout << "\tNUMBER IS NOT THE SAME!!" << endl;
            AllOk = kFALSE;
        }
        else
        {
            myProgressBar_t bar(NEvents_s);
            cout << "\tChecking separate events:" << endl;
            for(Int_t i=0; i<NEvents_s; i++)
            {
                EE_s->GetEvent(i);
                EE_l->GetEvent(i);

                AllOkVec[0] = AreVectorsTheSame(EE_s->Electron_pT, EE_l->Electron_pT, "pT"); if (AllOkVec[0]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[1] = AreVectorsTheSame(EE_s->Electron_eta, EE_l->Electron_eta, "eta"); if (AllOkVec[1]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[2] = AreVectorsTheSame(EE_s->Electron_phi, EE_l->Electron_phi, "phi"); if (AllOkVec[2]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[3] = AreVectorsTheSame(EE_s->Electron_Energy, EE_l->Electron_Energy, "Energy"); if (AllOkVec[3]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[4] = AreVectorsTheSame(EE_s->Electron_charge, EE_l->Electron_charge, "charge"); if (AllOkVec[4]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[5] = AreValuesTheSame(EE_s->Electron_InvM, EE_l->Electron_InvM, "InvM"); if (AllOkVec[5]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[6] = AreValuesTheSame(EE_s->nPileUp, EE_l->nPileUp, "nPileUp"); if (AllOkVec[6]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[7] = AreValuesTheSame(EE_s->GENEvt_weight, EE_l->GENEvt_weight, "event weight"); if (AllOkVec[7]==kFALSE && BreakIfProblem==kTRUE) break;

                for (Int_t j=0; j<8; j++)
                {
                    if(AllOkVec[j]==kFALSE)
                    {
                        AllOkV = kFALSE;
                        cout << "\tProblems in entry no." << i << endl;
                    }
                }

                // -- Normalization -- //
//                Double_t TotWeight_s = MuMu_s->GENEvt_weight;
//                Double_t TotWeight_l = MuMu_l->GENEvt_weight;
//                TotWeight_s = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_s->GENEvt_weight;
//                TotWeight_l = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_l->GENEvt_weight;

                // -- Histogram filling -- //
                for (Int_t j=0; j<2; j++)
                {
                    h1_pT_s->Fill(EE_s->Electron_pT->at(j), EE_s->GENEvt_weight);
                    h1_pT_l->Fill(EE_l->Electron_pT->at(j), EE_l->GENEvt_weight);
                    h1_eta_s->Fill(EE_s->Electron_eta->at(j), EE_s->GENEvt_weight);
                    h1_eta_l->Fill(EE_l->Electron_eta->at(j), EE_l->GENEvt_weight);
                    h1_phi_s->Fill(EE_s->Electron_phi->at(j), EE_s->GENEvt_weight);
                    h1_phi_l->Fill(EE_l->Electron_phi->at(j), EE_l->GENEvt_weight);
                    h1_Energy_s->Fill(EE_s->Electron_Energy->at(j), EE_s->GENEvt_weight);
                    h1_Energy_l->Fill(EE_l->Electron_Energy->at(j), EE_l->GENEvt_weight);
                }
                h1_nPileUp_s->Fill(EE_s->nPileUp, EE_s->GENEvt_weight);
                h1_nPileUp_l->Fill(EE_l->nPileUp, EE_l->GENEvt_weight);
                h1_invm_s->Fill(EE_s->Electron_InvM, EE_s->GENEvt_weight);
                h1_invm_l->Fill(EE_l->Electron_InvM, EE_l->GENEvt_weight);
                bar.Draw(i);
            } //End of event iteration
        } //End of if( same event numbers )

        if ( AllOk==kTRUE && AllOkV==kTRUE ) cout << "\tNo problems found so far." << endl;
        else cout << "\tProblems were detected." << endl;
        cout << "\tChecking histograms:  ";

        // -- Histogram checking -- //
        AllOkHist[0] = AreHistogramsTheSame(h1_pT_s, h1_pT_l, "pT", BreakIfProblem);
        AllOkHist[1] = AreHistogramsTheSame(h1_eta_s, h1_eta_l, "eta", BreakIfProblem);
        AllOkHist[2] = AreHistogramsTheSame(h1_phi_s, h1_phi_l, "phi", BreakIfProblem);
        AllOkHist[3] = AreHistogramsTheSame(h1_Energy_s, h1_Energy_l, "Energy", BreakIfProblem);
        AllOkHist[4] = AreHistogramsTheSame(h1_invm_s, h1_invm_l, "InvM", BreakIfProblem);
        AllOkHist[5] = AreHistogramsTheSame(h1_nPileUp_s, h1_nPileUp_l, "phi", BreakIfProblem);

        for ( Int_t j=0; j<6; j++ )
        {
            if ( AllOkHist[j]==kFALSE ) AllOkH = kFALSE;
        }

        // -- Drawing histos -- //
        if ( DrawHistos==kTRUE )
        {
            TCanvas* c_pT = new TCanvas ("pT", "pT", 1000, 1000);
            TCanvas* c_eta = new TCanvas ("eta", "eta", 1000, 1000);
            TCanvas* c_phi = new TCanvas ("phi", "phi", 1000, 1000);
            TCanvas* c_Energy = new TCanvas ("Energy", "Energy", 1000, 1000);
            TCanvas* c_invm = new TCanvas ("invm", "invm", 1000, 1000);
            TCanvas* c_nPileUp = new TCanvas ("nPileUp", "nPileUp", 1000, 1000);

            TLegend* leg = new TLegend (0.1, 0.8, 0.3, 0.9);

            h1_pT_s->SetLineColor(kBlue);
            h1_pT_s->SetLineWidth(2);
            h1_pT_l->SetLineColor(kRed);
            h1_eta_s->SetLineColor(kBlue);
            h1_eta_s->SetLineWidth(2);
            h1_eta_l->SetLineColor(kRed);
            h1_phi_s->SetLineColor(kBlue);
            h1_phi_s->SetLineWidth(2);
            h1_phi_l->SetLineColor(kRed);
            h1_Energy_s->SetLineColor(kBlue);
            h1_Energy_s->SetLineWidth(2);
            h1_Energy_l->SetLineColor(kRed);
            h1_invm_s->SetLineColor(kBlue);
            h1_invm_s->SetLineWidth(2);
            h1_invm_l->SetLineColor(kRed);
            h1_nPileUp_s->SetLineColor(kBlue);
            h1_nPileUp_s->SetLineWidth(2);
            h1_nPileUp_l->SetLineColor(kRed);

            leg->AddEntry(h1_pT_s, "Selected EE", "l");
            leg->AddEntry(h1_pT_l, "Long_Selected EE", "l");

            c_pT->cd(); h1_pT_s->Draw(); h1_pT_l->Draw("SAME"); leg->Draw(); c_pT->Update();
            c_eta->cd(); h1_eta_s->Draw(); h1_eta_l->Draw("SAME"); leg->Draw(); c_eta->Update();
            c_phi->cd(); h1_phi_s->Draw(); h1_phi_l->Draw("SAME"); leg->Draw(); c_phi->Update();
            c_Energy->cd(); h1_Energy_s->Draw(); h1_Energy_l->Draw("SAME"); leg->Draw(); c_Energy->Update();
            c_invm->cd(); h1_invm_s->Draw(); h1_invm_l->Draw("SAME"); leg->Draw(); c_invm->Update();
            c_nPileUp->cd(); h1_nPileUp_s->Draw(); h1_nPileUp_l->Draw("SAME"); leg->Draw(); c_nPileUp->Update();
        }

        if ( AllOkH==kTRUE ) cout << "All bin values match." << endl;
        else cout << "Problems were detected." << endl;

        if ( AllOk==kTRUE && AllOkV==kTRUE && AllOkH==kTRUE) cout << "Selected and LongSelected events match perfectly! Hooray!" << endl;
    } //end of i_tup iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds." << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

}// end of CheckSelectedEE


/// -------------------------------- Muon Channel ------------------------------------ ///
void CheckSelectedMuMu ( Bool_t DrawHistos, Bool_t BreakIfProblem )
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0; // Need checking
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "\n[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

    TStopwatch totaltime;
    totaltime.Start();

    // -- Each ntuple directory & corresponding Tags -- //
//    vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    Tag.push_back("ZMuMu_M4500to6000");
    nEvents.push_back(100000);  // Primary dataset actually contained 19672 events.
    Xsec.push_back(4.56E-07);   // Needs checking
    const Int_t Ntup = 1;
    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

//        cout << "\t<" << [i_tup] << ">" << endl;

        TChain *ch_s = new TChain("DYTree");
        TChain *ch_l = new TChain("DYTree");
        ch_s->Add("/media/sf_DATA/test/SelectedMuMu_ZToMuMu_M4500to6000_4.root");
        ch_l->Add("/media/sf_DATA/test/LongSelectedMuMu_ZToMuMu_M4500to6000_4.root");

        SelectedMuMu_t *MuMu_s = new SelectedMuMu_t();
        LongSelectedMuMu_t *MuMu_l = new LongSelectedMuMu_t();
        MuMu_s->CreateFromChain(ch_s);
        MuMu_l->CreateFromChain(ch_l);

        TH1D* h1_pT_s = new TH1D("h1pts", "h1_pT_s", 400, 0, 4000);
        TH1D* h1_pT_l = new TH1D("h1ptl", "h1_pT_l", 400, 0, 4000);
        TH1D* h1_TuneP_pT_s = new TH1D("h1tuneppts", "h1_TuneP_pT_s", 400, 0, 4000);
        TH1D* h1_TuneP_pT_l = new TH1D("h1tunepptl", "h1_TuneP_pT_l", 400, 0, 4000);
        TH1D* h1_eta_s = new TH1D("h1etas", "h1_eta_s", 120, -3, 3);
        TH1D* h1_eta_l = new TH1D("h1etal", "h1_eta_l", 120, -3, 3);
        TH1D* h1_TuneP_eta_s = new TH1D("h1tunepetas", "h1_TuneP_eta_s", 120, -3, 3);
        TH1D* h1_TuneP_eta_l = new TH1D("h1tunepetal", "h1_TuneP_eta_l", 120, -3, 3);
        TH1D* h1_phi_s = new TH1D("h1phis", "h1_phi_s", 400, -4, 4);
        TH1D* h1_phi_l = new TH1D("h1phil", "h1_phi_l", 400, -4, 4);
        TH1D* h1_TuneP_phi_s = new TH1D("h1tunepphis", "h1_TuneP_phi_s", 400, -4, 4);
        TH1D* h1_TuneP_phi_l = new TH1D("h1tunepphil", "h1_TuneP_phi_l", 400, -4, 4);
        TH1D* h1_Energy_s = new TH1D("h1energys", "h1_Energy_s", 400, 0, 8000);
        TH1D* h1_Energy_l = new TH1D("h1energyl", "h1_Energy_l", 400, 0, 8000);
        TH1D* h1_invm_s = new TH1D("h1invms", "h1_invm_s", 350, 0, 7000);
        TH1D* h1_invm_l = new TH1D("h1invml", "h1_invm_l", 350, 0, 7000);
        TH1D* h1_nPileUp_s = new TH1D("h1npileups", "h1_nPileUp_s", 60, 0, 60);
        TH1D* h1_nPileUp_l = new TH1D("h1npileupl", "h1_nPileUp_l", 60, 0, 60);

        Bool_t AllOkVec[11];
        Bool_t AllOkHist[9];
        Bool_t AllOkV = kTRUE, AllOkH = kTRUE, AllOk = kTRUE;

        Int_t NEvents_s = ch_s->GetEntries();
        Int_t NEvents_l = ch_l->GetEntries();
        cout << "\tTotal selected Events: " << NEvents_s << endl;
        cout << "\tTotal LongSelected Events: " << NEvents_l << endl;
        if(NEvents_s!=NEvents_l)
        {
            cout << "\tNUMBER IS NOT THE SAME!!" << endl;
            AllOk = kFALSE;
        }
        else
        {
            myProgressBar_t bar(NEvents_s);
            cout << "\tChecking separate events:" << endl;
            for(Int_t i=0; i<NEvents_s; i++)
            {
                MuMu_s->GetEvent(i);
                MuMu_l->GetEvent(i);

                AllOkVec[0] = AreVectorsTheSame(MuMu_s->Muon_pT, MuMu_l->Muon_pT, "pT"); if (AllOkVec[0]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[1] = AreVectorsTheSame(MuMu_s->Muon_TuneP_pT, MuMu_l->Muon_TuneP_pT, "TuneP_pT"); if (AllOkVec[1]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[2] = AreVectorsTheSame(MuMu_s->Muon_eta, MuMu_l->Muon_eta, "eta"); if (AllOkVec[2]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[3] = AreVectorsTheSame(MuMu_s->Muon_TuneP_eta, MuMu_l->Muon_TuneP_eta, "TuneP_eta"); if (AllOkVec[3]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[4] = AreVectorsTheSame(MuMu_s->Muon_phi, MuMu_l->Muon_phi, "phi"); if (AllOkVec[4]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[5] = AreVectorsTheSame(MuMu_s->Muon_TuneP_phi, MuMu_l->Muon_TuneP_phi, "TuneP_phi"); if (AllOkVec[5]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[6] = AreVectorsTheSame(MuMu_s->Muon_Energy, MuMu_l->Muon_Energy, "Energy"); if (AllOkVec[6]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[7] = AreVectorsTheSame(MuMu_s->Muon_charge, MuMu_l->Muon_charge, "charge"); if (AllOkVec[7]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[8] = AreValuesTheSame(MuMu_s->Muon_InvM, MuMu_l->Muon_InvM, "InvM"); if (AllOkVec[8]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[9] = AreValuesTheSame(MuMu_s->nPileUp, MuMu_l->nPileUp, "nPileUp"); if (AllOkVec[9]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[10] = AreValuesTheSame(MuMu_s->GENEvt_weight, MuMu_l->GENEvt_weight, "event weight"); if (AllOkVec[10]==kFALSE && BreakIfProblem==kTRUE) break;

                for (Int_t j=0; j<11; j++)
                {
                    if(AllOkVec[j]==kFALSE)
                    {
                        AllOkV = kFALSE;
                        cout << "\tProblems in entry no." << i << endl;
                    }
                }

                // -- Normalization -- //
//                Double_t TotWeight_s = MuMu_s->GENEvt_weight;
//                Double_t TotWeight_l = MuMu_l->GENEvt_weight;
//                TotWeight_s = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_s->GENEvt_weight;
//                TotWeight_l = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_l->GENEvt_weight;

                // -- Histogram filling -- //
                for (Int_t j=0; j<2; j++)
                {
                    h1_pT_s->Fill(MuMu_s->Muon_pT->at(j), MuMu_s->GENEvt_weight);
                    h1_pT_l->Fill(MuMu_l->Muon_pT->at(j), MuMu_l->GENEvt_weight);
                    h1_TuneP_pT_s->Fill(MuMu_s->Muon_TuneP_pT->at(j), MuMu_s->GENEvt_weight);
                    h1_TuneP_pT_l->Fill(MuMu_l->Muon_TuneP_pT->at(j), MuMu_l->GENEvt_weight);
                    h1_eta_s->Fill(MuMu_s->Muon_eta->at(j), MuMu_s->GENEvt_weight);
                    h1_eta_l->Fill(MuMu_l->Muon_eta->at(j), MuMu_l->GENEvt_weight);
                    h1_TuneP_eta_s->Fill(MuMu_s->Muon_TuneP_eta->at(j), MuMu_s->GENEvt_weight);
                    h1_TuneP_eta_l->Fill(MuMu_l->Muon_TuneP_eta->at(j), MuMu_l->GENEvt_weight);
                    h1_phi_s->Fill(MuMu_s->Muon_phi->at(j), MuMu_s->GENEvt_weight);
                    h1_phi_l->Fill(MuMu_l->Muon_phi->at(j), MuMu_l->GENEvt_weight);
                    h1_TuneP_phi_s->Fill(MuMu_s->Muon_TuneP_phi->at(j), MuMu_s->GENEvt_weight);
                    h1_TuneP_phi_l->Fill(MuMu_l->Muon_TuneP_phi->at(j), MuMu_l->GENEvt_weight);
                    h1_Energy_s->Fill(MuMu_s->Muon_Energy->at(j), MuMu_s->GENEvt_weight);
                    h1_Energy_l->Fill(MuMu_l->Muon_Energy->at(j), MuMu_l->GENEvt_weight);
                }
                h1_nPileUp_s->Fill(MuMu_s->nPileUp, MuMu_s->GENEvt_weight);
                h1_nPileUp_l->Fill(MuMu_l->nPileUp, MuMu_l->GENEvt_weight);
                h1_invm_s->Fill(MuMu_s->Muon_InvM, MuMu_s->GENEvt_weight);
                h1_invm_l->Fill(MuMu_l->Muon_InvM, MuMu_l->GENEvt_weight);
                bar.Draw(i);
            } //End of event iteration
        } //End of if( same event numbers )

        if ( AllOk==kTRUE && AllOkV==kTRUE ) cout << "\tNo problems found so far." << endl;
        else cout << "\tProblems were detected." << endl;
        cout << "\tChecking histograms:  ";

        // -- Histogram checking -- //
        AllOkHist[0] = AreHistogramsTheSame(h1_pT_s, h1_pT_l, "pT", BreakIfProblem);
        AllOkHist[1] = AreHistogramsTheSame(h1_TuneP_pT_s, h1_TuneP_pT_l, "TuneP_pT", BreakIfProblem);
        AllOkHist[2] = AreHistogramsTheSame(h1_eta_s, h1_eta_l, "eta", BreakIfProblem);
        AllOkHist[3] = AreHistogramsTheSame(h1_TuneP_eta_s, h1_TuneP_eta_l, "TuneP_eta", BreakIfProblem);
        AllOkHist[4] = AreHistogramsTheSame(h1_phi_s, h1_phi_l, "phi", BreakIfProblem);
        AllOkHist[5] = AreHistogramsTheSame(h1_TuneP_phi_s, h1_TuneP_phi_l, "TuneP_phi", BreakIfProblem);
        AllOkHist[6] = AreHistogramsTheSame(h1_Energy_s, h1_Energy_l, "Energy", BreakIfProblem);
        AllOkHist[7] = AreHistogramsTheSame(h1_invm_s, h1_invm_l, "InvM", BreakIfProblem);
        AllOkHist[8] = AreHistogramsTheSame(h1_nPileUp_s, h1_nPileUp_l, "phi", BreakIfProblem);

        for ( Int_t j=0; j<9; j++ )
        {
            if ( AllOkHist[j]==kFALSE ) AllOkH = kFALSE;
        }

        // -- Drawing histos -- //
        if ( DrawHistos==kTRUE )
        {
            TCanvas* c_pT = new TCanvas ("pT", "pT", 1000, 1000);
            TCanvas* c_TuneP_pT = new TCanvas ("TuneP_pT", "TuneP_pT", 1000, 1000);
            TCanvas* c_eta = new TCanvas ("eta", "eta", 1000, 1000);
            TCanvas* c_TuneP_eta = new TCanvas ("TuneP_eta", "TuneP_eta", 1000, 1000);
            TCanvas* c_phi = new TCanvas ("phi", "phi", 1000, 1000);
            TCanvas* c_TuneP_phi = new TCanvas ("TuneP_phi", "TuneP_phi", 1000, 1000);
            TCanvas* c_Energy = new TCanvas ("Energy", "Energy", 1000, 1000);
            TCanvas* c_invm = new TCanvas ("invm", "invm", 1000, 1000);
            TCanvas* c_nPileUp = new TCanvas ("nPileUp", "nPileUp", 1000, 1000);

            TLegend* leg = new TLegend (0.1, 0.8, 0.3, 0.9);

            h1_pT_s->SetLineColor(kBlue);
            h1_pT_s->SetLineWidth(2);
            h1_pT_l->SetLineColor(kRed);
            h1_TuneP_pT_s->SetLineColor(kBlue);
            h1_TuneP_pT_s->SetLineWidth(2);
            h1_TuneP_pT_l->SetLineColor(kRed);
            h1_eta_s->SetLineColor(kBlue);
            h1_eta_s->SetLineWidth(2);
            h1_eta_l->SetLineColor(kRed);
            h1_TuneP_eta_s->SetLineColor(kBlue);
            h1_TuneP_eta_s->SetLineWidth(2);
            h1_TuneP_eta_l->SetLineColor(kRed);
            h1_phi_s->SetLineColor(kBlue);
            h1_phi_s->SetLineWidth(2);
            h1_phi_l->SetLineColor(kRed);
            h1_TuneP_phi_s->SetLineColor(kBlue);
            h1_TuneP_phi_s->SetLineWidth(2);
            h1_TuneP_phi_l->SetLineColor(kRed);
            h1_Energy_s->SetLineColor(kBlue);
            h1_Energy_s->SetLineWidth(2);
            h1_Energy_l->SetLineColor(kRed);
            h1_invm_s->SetLineColor(kBlue);
            h1_invm_s->SetLineWidth(2);
            h1_invm_l->SetLineColor(kRed);
            h1_nPileUp_s->SetLineColor(kBlue);
            h1_nPileUp_s->SetLineWidth(2);
            h1_nPileUp_l->SetLineColor(kRed);

            leg->AddEntry(h1_pT_s, "Selected MuMu", "l");
            leg->AddEntry(h1_pT_l, "Long_Selected MuMu", "l");

            c_pT->cd(); h1_pT_s->Draw(); h1_pT_l->Draw("SAME"); leg->Draw(); c_pT->Update();
            c_TuneP_pT->cd(); h1_TuneP_pT_s->Draw(); h1_TuneP_pT_l->Draw("SAME"); leg->Draw(); c_TuneP_pT->Update();
            c_eta->cd(); h1_eta_s->Draw(); h1_eta_l->Draw("SAME"); leg->Draw(); c_eta->Update();
            c_TuneP_eta->cd(); h1_TuneP_eta_s->Draw(); h1_TuneP_eta_l->Draw("SAME"); leg->Draw(); c_TuneP_eta->Update();
            c_phi->cd(); h1_phi_s->Draw(); h1_phi_l->Draw("SAME"); leg->Draw(); c_phi->Update();
            c_TuneP_phi->cd(); h1_TuneP_phi_s->Draw(); h1_TuneP_phi_l->Draw("SAME"); leg->Draw(); c_TuneP_phi->Update();
            c_Energy->cd(); h1_Energy_s->Draw(); h1_Energy_l->Draw("SAME"); leg->Draw(); c_Energy->Update();
            c_invm->cd(); h1_invm_s->Draw(); h1_invm_l->Draw("SAME"); leg->Draw(); c_invm->Update();
            c_nPileUp->cd(); h1_nPileUp_s->Draw(); h1_nPileUp_l->Draw("SAME"); leg->Draw(); c_nPileUp->Update();
        }

        if ( AllOkH==kTRUE ) cout << "All bin values match." << endl;
        else cout << "Problems were detected." << endl;

        if ( AllOk==kTRUE && AllOkV==kTRUE && AllOkH==kTRUE) cout << "Selected and LongSelected events match perfectly! Hooray!" << endl;
    } //end of i_tup iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds." << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

}// end of CheckSelectedMuMu


/// --------------------------------- EMu events --------------------------------- ///
void CheckSelectedEMu(Bool_t DrawHistos, Bool_t BreakIfProblem )
{
    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0; // Need checking
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "\n[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

    TStopwatch totaltime;
    totaltime.Start();

    // -- Each ntuple directory & corresponding Tags -- //
//    vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
    Tag.push_back("WW");
    nEvents.push_back(6987123);
    Xsec.push_back(118.7);
    const Int_t Ntup = 1;
    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

//        cout << "\t<" << [i_tup] << ">" << endl;

        TChain *ch_s = new TChain("DYTree");
        TChain *ch_l = new TChain("DYTree");
        ch_s->Add("/media/sf_DATA/test/SelectedEMu_WW_34.root");
        ch_l->Add("/media/sf_DATA/test/LongSelectedEMu_WW_34.root");

        SelectedEMu_t *EMu_s = new SelectedEMu_t();
        LongSelectedEMu_t *EMu_l = new LongSelectedEMu_t();
        EMu_s->CreateFromChain(ch_s);
        EMu_l->CreateFromChain(ch_l);

        TH1D* h1_pT_ele_s = new TH1D("h1pteles", "h1_pT_ele_s", 400, 0, 4000);
        TH1D* h1_pT_ele_l = new TH1D("h1ptelel", "h1_pT_ele_l", 400, 0, 4000);
        TH1D* h1_eta_ele_s = new TH1D("h1etaeles", "h1_eta_ele_s", 120, -3, 3);
        TH1D* h1_eta_ele_l = new TH1D("h1etaelel", "h1_eta_ele_l", 120, -3, 3);
        TH1D* h1_phi_ele_s = new TH1D("h1phieles", "h1_phi_ele_s", 400, -4, 4);
        TH1D* h1_phi_ele_l = new TH1D("h1phielel", "h1_phi_ele_l", 400, -4, 4);
        TH1D* h1_Energy_ele_s = new TH1D("h1energyeles", "h1_Energy_ele_s", 400, 0, 8000);
        TH1D* h1_Energy_ele_l = new TH1D("h1energyelel", "h1_Energy_ele_l", 400, 0, 8000);
        TH1D* h1_pT_mu_s = new TH1D("h1ptmus", "h1_pT_mu_s", 400, 0, 4000);
        TH1D* h1_pT_mu_l = new TH1D("h1ptmul", "h1_pT_mu_l", 400, 0, 4000);
        TH1D* h1_eta_mu_s = new TH1D("h1etamus", "h1_eta_mu_s", 120, -3, 3);
        TH1D* h1_eta_mu_l = new TH1D("h1etamul", "h1_eta_mu_l", 120, -3, 3);
        TH1D* h1_phi_mu_s = new TH1D("h1phimus", "h1_phi_mu_s", 400, -4, 4);
        TH1D* h1_phi_mu_l = new TH1D("h1phimul", "h1_phi_mu_l", 400, -4, 4);
        TH1D* h1_Energy_mu_s = new TH1D("h1energymus", "h1_Energy_mu_s", 400, 0, 8000);
        TH1D* h1_Energy_mu_l = new TH1D("h1energymul", "h1_Energy_mu_l", 400, 0, 8000);
        TH1D* h1_invm_s = new TH1D("h1invms", "h1_invm_s", 350, 0, 7000);
        TH1D* h1_invm_l = new TH1D("h1invml", "h1_invm_l", 350, 0, 7000);
        TH1D* h1_nPileUp_s = new TH1D("h1npileups", "h1_nPileUp_s", 60, 0, 60);
        TH1D* h1_nPileUp_l = new TH1D("h1npileupl", "h1_nPileUp_l", 60, 0, 60);

        Bool_t AllOkVec[13];
        Bool_t AllOkHist[10];
        Bool_t AllOkV = kTRUE, AllOkH = kTRUE, AllOk = kTRUE;

        Int_t NEvents_s = ch_s->GetEntries();
        Int_t NEvents_l = ch_l->GetEntries();
        cout << "\tTotal selected Events: " << NEvents_s << endl;
        cout << "\tTotal LongSelected Events: " << NEvents_l << endl;
        if(NEvents_s!=NEvents_l)
        {
            cout << "\tNUMBER IS NOT THE SAME!!" << endl;
            AllOk = kFALSE;
        }
        else
        {
            myProgressBar_t bar(NEvents_s);
            cout << "\tChecking separate events:" << endl;
            for(Int_t i=0; i<NEvents_s; i++)
            {
                EMu_s->GetEvent(i);
                EMu_l->GetEvent(i);

                AllOkVec[0] = AreValuesTheSame(EMu_s->Electron_pT, EMu_l->Electron_pT, "Electron pT"); if (AllOkVec[0]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[1] = AreValuesTheSame(EMu_s->Electron_eta, EMu_l->Electron_eta, "Electron eta"); if (AllOkVec[1]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[2] = AreValuesTheSame(EMu_s->Electron_phi, EMu_l->Electron_phi, "Electron phi"); if (AllOkVec[2]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[3] = AreValuesTheSame(EMu_s->Electron_Energy, EMu_l->Electron_Energy, "Electron Energy"); if (AllOkVec[3]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[4] = AreValuesTheSame(EMu_s->Electron_charge, EMu_l->Electron_charge, "Electron charge"); if (AllOkVec[4]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[5] = AreValuesTheSame(EMu_s->Muon_pT, EMu_l->Muon_pT, "Muon pT"); if (AllOkVec[5]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[6] = AreValuesTheSame(EMu_s->Muon_eta, EMu_l->Muon_eta, "Muon eta"); if (AllOkVec[6]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[7] = AreValuesTheSame(EMu_s->Muon_phi, EMu_l->Muon_phi, "Muon phi"); if (AllOkVec[7]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[8] = AreValuesTheSame(EMu_s->Muon_Energy, EMu_l->Muon_Energy, "Muon Energy"); if (AllOkVec[8]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[9] = AreValuesTheSame(EMu_s->Muon_charge, EMu_l->Muon_charge, "Muon charge"); if (AllOkVec[9]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[10] = AreValuesTheSame(EMu_s->EMu_InvM, EMu_l->EMu_InvM, "InvM"); if (AllOkVec[10]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[11] = AreValuesTheSame(EMu_s->nPileUp, EMu_l->nPileUp, "nPileUp"); if (AllOkVec[11]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[12] = AreValuesTheSame(EMu_s->GENEvt_weight, EMu_l->GENEvt_weight, "event weight"); if (AllOkVec[12]==kFALSE && BreakIfProblem==kTRUE) break;

                for (Int_t j=0; j<12; j++)
                {
                    if(AllOkVec[j]==kFALSE)
                    {
                        AllOkV = kFALSE;
                        cout << "\tProblems in entry no." << i << endl;
                    }
                }

                // -- Normalization -- //
//                Double_t TotWeight_s = MuMu_s->GENEvt_weight;
//                Double_t TotWeight_l = MuMu_l->GENEvt_weight;
//                TotWeight_s = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_s->GENEvt_weight;
//                TotWeight_l = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_l->GENEvt_weight;

                // -- Histogram filling -- //
                h1_pT_ele_s->Fill(EMu_s->Electron_pT, EMu_s->GENEvt_weight);
                h1_pT_ele_l->Fill(EMu_l->Electron_pT, EMu_l->GENEvt_weight);
                h1_eta_ele_s->Fill(EMu_s->Electron_eta, EMu_s->GENEvt_weight);
                h1_eta_ele_l->Fill(EMu_l->Electron_eta, EMu_l->GENEvt_weight);
                h1_phi_ele_s->Fill(EMu_s->Electron_phi, EMu_s->GENEvt_weight);
                h1_phi_ele_l->Fill(EMu_l->Electron_phi, EMu_l->GENEvt_weight);
                h1_Energy_ele_s->Fill(EMu_s->Electron_Energy, EMu_s->GENEvt_weight);
                h1_Energy_ele_l->Fill(EMu_l->Electron_Energy, EMu_l->GENEvt_weight);
                h1_pT_mu_s->Fill(EMu_s->Muon_pT, EMu_s->GENEvt_weight);
                h1_pT_mu_l->Fill(EMu_l->Muon_pT, EMu_l->GENEvt_weight);
                h1_eta_mu_s->Fill(EMu_s->Muon_eta, EMu_s->GENEvt_weight);
                h1_eta_mu_l->Fill(EMu_l->Muon_eta, EMu_l->GENEvt_weight);
                h1_phi_mu_s->Fill(EMu_s->Muon_phi, EMu_s->GENEvt_weight);
                h1_phi_mu_l->Fill(EMu_l->Muon_phi, EMu_l->GENEvt_weight);
                h1_Energy_mu_s->Fill(EMu_s->Muon_Energy, EMu_s->GENEvt_weight);
                h1_Energy_mu_l->Fill(EMu_l->Muon_Energy, EMu_l->GENEvt_weight);
                h1_nPileUp_s->Fill(EMu_s->nPileUp, EMu_s->GENEvt_weight);
                h1_nPileUp_l->Fill(EMu_l->nPileUp, EMu_l->GENEvt_weight);
                h1_invm_s->Fill(EMu_s->EMu_InvM, EMu_s->GENEvt_weight);
                h1_invm_l->Fill(EMu_l->EMu_InvM, EMu_l->GENEvt_weight);
                bar.Draw(i);
            } //End of event iteration
        } //End of if( same event numbers )

        if ( AllOk==kTRUE && AllOkV==kTRUE ) cout << "\tNo problems found so far." << endl;
        else cout << "\tProblems were detected." << endl;
        cout << "\tChecking histograms:  ";

        // -- Histogram checking -- //
        AllOkHist[0] = AreHistogramsTheSame(h1_pT_ele_s, h1_pT_ele_l, "Electron pT", BreakIfProblem);
        AllOkHist[1] = AreHistogramsTheSame(h1_eta_ele_s, h1_eta_ele_l, "Electron eta", BreakIfProblem);
        AllOkHist[2] = AreHistogramsTheSame(h1_phi_ele_s, h1_phi_ele_l, "Electron phi", BreakIfProblem);
        AllOkHist[3] = AreHistogramsTheSame(h1_Energy_ele_s, h1_Energy_ele_l, "Electron Energy", BreakIfProblem);
        AllOkHist[4] = AreHistogramsTheSame(h1_pT_mu_s, h1_pT_mu_l, "Muon pT", BreakIfProblem);
        AllOkHist[5] = AreHistogramsTheSame(h1_eta_mu_s, h1_eta_mu_l, "Muon eta", BreakIfProblem);
        AllOkHist[6] = AreHistogramsTheSame(h1_phi_mu_s, h1_phi_mu_l, "Muon phi", BreakIfProblem);
        AllOkHist[7] = AreHistogramsTheSame(h1_Energy_mu_s, h1_Energy_mu_l, "Muon Energy", BreakIfProblem);
        AllOkHist[8] = AreHistogramsTheSame(h1_invm_s, h1_invm_l, "InvM", BreakIfProblem);
        AllOkHist[9] = AreHistogramsTheSame(h1_nPileUp_s, h1_nPileUp_l, "phi", BreakIfProblem);

        for ( Int_t j=0; j<9; j++ )
        {
            if ( AllOkHist[j]==kFALSE ) AllOkH = kFALSE;
        }

        // -- Drawing histos -- //
        if ( DrawHistos==kTRUE )
        {
            TCanvas* c_pT_ele = new TCanvas ("Electron pT", "Electron pT", 1000, 1000);
            TCanvas* c_eta_ele = new TCanvas ("Electron eta", "Electron eta", 1000, 1000);
            TCanvas* c_phi_ele = new TCanvas ("Electron phi", "Electron phi", 1000, 1000);
            TCanvas* c_Energy_ele = new TCanvas ("Electron Energy", "Electron Energy", 1000, 1000);
            TCanvas* c_pT_mu = new TCanvas ("Muon pT", "Muon pT", 1000, 1000);
            TCanvas* c_eta_mu = new TCanvas ("Muon eta", "Muon eta", 1000, 1000);
            TCanvas* c_phi_mu = new TCanvas ("Muon phi", "Muon phi", 1000, 1000);
            TCanvas* c_Energy_mu = new TCanvas ("Muon Energy", "Muon Energy", 1000, 1000);
            TCanvas* c_invm = new TCanvas ("invm", "invm", 1000, 1000);
            TCanvas* c_nPileUp = new TCanvas ("nPileUp", "nPileUp", 1000, 1000);

            TLegend* leg = new TLegend (0.1, 0.8, 0.3, 0.9);

            h1_pT_ele_s->SetLineColor(kBlue);
            h1_pT_ele_s->SetLineWidth(2);
            h1_pT_ele_l->SetLineColor(kRed);
            h1_eta_ele_s->SetLineColor(kBlue);
            h1_eta_ele_s->SetLineWidth(2);
            h1_eta_ele_l->SetLineColor(kRed);
            h1_phi_ele_s->SetLineColor(kBlue);
            h1_phi_ele_s->SetLineWidth(2);
            h1_phi_ele_l->SetLineColor(kRed);
            h1_Energy_ele_s->SetLineColor(kBlue);
            h1_Energy_ele_s->SetLineWidth(2);
            h1_Energy_ele_l->SetLineColor(kRed);
            h1_pT_mu_s->SetLineColor(kBlue);
            h1_pT_mu_s->SetLineWidth(2);
            h1_pT_mu_l->SetLineColor(kRed);
            h1_eta_mu_s->SetLineColor(kBlue);
            h1_eta_mu_s->SetLineWidth(2);
            h1_eta_mu_l->SetLineColor(kRed);
            h1_phi_mu_s->SetLineColor(kBlue);
            h1_phi_mu_s->SetLineWidth(2);
            h1_phi_mu_l->SetLineColor(kRed);
            h1_Energy_mu_s->SetLineColor(kBlue);
            h1_Energy_mu_s->SetLineWidth(2);
            h1_Energy_mu_l->SetLineColor(kRed);
            h1_invm_s->SetLineColor(kBlue);
            h1_invm_s->SetLineWidth(2);
            h1_invm_l->SetLineColor(kRed);
            h1_nPileUp_s->SetLineColor(kBlue);
            h1_nPileUp_s->SetLineWidth(2);
            h1_nPileUp_l->SetLineColor(kRed);

            leg->AddEntry(h1_pT_ele_s, "Selected EMu", "l");
            leg->AddEntry(h1_pT_ele_l, "Long_Selected EMu", "l");

            c_pT_ele->cd(); h1_pT_ele_s->Draw(); h1_pT_ele_l->Draw("SAME"); leg->Draw(); c_pT_ele->Update();
            c_eta_ele->cd(); h1_eta_ele_s->Draw(); h1_eta_ele_l->Draw("SAME"); leg->Draw(); c_eta_ele->Update();
            c_phi_ele->cd(); h1_phi_ele_s->Draw(); h1_phi_ele_l->Draw("SAME"); leg->Draw(); c_phi_ele->Update();
            c_Energy_ele->cd(); h1_Energy_ele_s->Draw(); h1_Energy_ele_l->Draw("SAME"); leg->Draw(); c_Energy_ele->Update();
            c_pT_mu->cd(); h1_pT_mu_s->Draw(); h1_pT_mu_l->Draw("SAME"); leg->Draw(); c_pT_mu->Update();
            c_eta_mu->cd(); h1_eta_mu_s->Draw(); h1_eta_mu_l->Draw("SAME"); leg->Draw(); c_eta_mu->Update();
            c_phi_mu->cd(); h1_phi_mu_s->Draw(); h1_phi_mu_l->Draw("SAME"); leg->Draw(); c_phi_mu->Update();
            c_Energy_mu->cd(); h1_Energy_mu_s->Draw(); h1_Energy_mu_l->Draw("SAME"); leg->Draw(); c_Energy_mu->Update();
            c_invm->cd(); h1_invm_s->Draw(); h1_invm_l->Draw("SAME"); leg->Draw(); c_invm->Update();
            c_nPileUp->cd(); h1_nPileUp_s->Draw(); h1_nPileUp_l->Draw("SAME"); leg->Draw(); c_nPileUp->Update();
        }

        if ( AllOkH==kTRUE ) cout << "All bin values match." << endl;
        else cout << "Problems were detected." << endl;

        if ( AllOk==kTRUE && AllOkV==kTRUE && AllOkH==kTRUE) cout << "Selected and LongSelected events match perfectly! Hooray!" << endl;
    } //end of i_tup iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds." << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

}// end of CheckSelectedEMu


Bool_t AreVectorsTheSame( vector<Double_t> *vec1, vector<Double_t> *vec2, TString name )
{
    Bool_t areThey = kTRUE;
    if ( vec1->size()!=vec2->size() )
    {
        cout << name << " vector sizes do not match!!" << endl;
        areThey = kFALSE;
    }
    else
    {
        for ( UInt_t i=0; i<vec1->size(); i++ )
        {
            if( vec1->at(i)!=vec2->at(i) )
            {
                cout << name << "[" << i << "] values do not match!!" << endl;
                areThey = kFALSE;
            }
        }
    }
    return areThey;
}

Bool_t AreVectorsTheSame( vector<Int_t> *vec1, vector<Int_t> *vec2, TString name )
{
    Bool_t areThey = kTRUE;
    if ( vec1->size()!=vec2->size() )
    {
        cout << name << " vector sizes do not match!!" << endl;
        areThey = kFALSE;
    }
    else
    {
        for ( UInt_t i=0; i<vec1->size(); i++ )
        {
            if( vec1->at(i)!=vec2->at(i) )
            {
                cout << name << "[" << i << "] values do not match!!" << endl;
                areThey = kFALSE;
            }
        }
    }
    return areThey;
}

Bool_t AreValuesTheSame( Double_t val1, Double_t val2, TString name )
{
    if( val1!=val2 )
    {
        cout << name << " values do not match!!" << endl;
        return kFALSE;
    }
    else return kTRUE;
}

Bool_t AreValuesTheSame( Int_t val1, Int_t val2, TString name )
{
    if( val1!=val2 )
    {
        cout << name << " values do not match!!" << endl;
        return kFALSE;
    }
    else return kTRUE;
}

Bool_t AreValuesTheSame( Bool_t val1, Bool_t val2, TString name )
{
    if( val1!=val2 )
    {
        cout << name << " values do not match!!" << endl;
        return kFALSE;
    }
    else return kTRUE;
}

Bool_t AreHistogramsTheSame( TH1D *hist1, TH1D *hist2, TString name, Bool_t setBreak )
{
    Bool_t areThey = kTRUE;
    if( hist1->GetSize()!=hist2->GetSize() )
    {
        cout << name << " histogram sizes do not match!!" << endl;
        areThey = kFALSE;
    }
    else
    {
        for ( Int_t i=1; i<hist1->GetSize()-1; i++ )
        {
            if( hist1->GetBinContent(i)!=hist2->GetBinContent(i) )
            {
                cout << "Histogram " << name << " bin " << i << ": Values do not match!!" << endl;
                areThey = kFALSE;
                if ( setBreak==kTRUE ) break;
            }
        }
    }
    return areThey;
}
