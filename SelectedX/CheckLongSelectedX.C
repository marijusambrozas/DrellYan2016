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

void CheckLongSelectedEE(Bool_t DrawHistos, Bool_t BreakIfProblem);
void CheckLongSelectedMuMu(Bool_t DrawHistos, Bool_t BreakIfProblem);
void CheckLongSelectedEMu(Bool_t DrawHistos, Bool_t BreakIfProblem);


void CheckLongSelectedX ( TString whichX, Bool_t DrawHistos = kFALSE, Bool_t BreakIfProblem = kFALSE )
{
    Int_t Xselected = 0;
    if ( whichX.Contains("EE") || whichX.Contains("ee") )
    {
        Xselected++;
        cout << "\n*********        CheckLongSelectedEE        *********" << endl;
        CheckLongSelectedEE(DrawHistos, BreakIfProblem);
    }
    if ( whichX.Contains("MuMu") || whichX.Contains("mumu") || whichX.Contains("MUMU") )
    {
        Xselected++;
        cout << "\n*********       CheckLongSelectedMuMu       *********" << endl;
        CheckLongSelectedMuMu(DrawHistos, BreakIfProblem);
    }
    if ( whichX.Contains("EMu") || whichX.Contains("emu") || whichX.Contains("Emu") || whichX.Contains("eMu") || whichX.Contains("EMU") )
    {
        Xselected++;
        cout << "\n*********       CheckLongSelectedEMu        *********" << endl;
        CheckLongSelectedEMu(DrawHistos, BreakIfProblem);
    }
    if ( Xselected == 0 ) cout << "Wrong arument!" << endl;
}


/// ----------------------------- Electron Channel ------------------------------ ///
void CheckLongSelectedEE ( Bool_t DrawHistos, Bool_t BreakIfProblem )
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
    Tag.push_back("ZToEE_M4500to6000");     // There is no such type mentioned in DYAnalyzer
    nEvents.push_back(74606);   // The actual number of events in file is taken
    Xsec.push_back(4.56E-07);   // Copied from ZToMuMu_M4500to6000
    const Int_t Ntup = 1;
    for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
    {
        TStopwatch looptime;
        looptime.Start();

//        cout << "\t<" << [i_tup] << ">" << endl;

        TChain *ch_s = new TChain("DYTree");
        TChain *ch_re = new TChain("DYTree");
        ch_s->Add("/media/sf_DATA/test/LongSelectedEE_ZToEE_M4500to6000_2.root");
        ch_re->Add("/media/sf_DATA/test/ReSelectedEE_ZToEE_M4500to6000_2.root");

        LongSelectedEE_t *EE_s = new LongSelectedEE_t();
        LongSelectedEE_t *EE_re = new LongSelectedEE_t();
        EE_s->CreateFromChain(ch_s);
        EE_re->CreateFromChain(ch_re);

        TH1D* h1_pT_s = new TH1D("h1pts", "h1_pT_s", 400, 0, 4000);
        TH1D* h1_pT_re = new TH1D("h1ptre", "h1_pT_re", 400, 0, 4000);
        TH1D* h1_eta_s = new TH1D("h1etas", "h1_eta_s", 120, -3, 3);
        TH1D* h1_eta_re = new TH1D("h1etare", "h1_eta_re", 120, -3, 3);
        TH1D* h1_invm_s = new TH1D("h1invms", "h1_invm_s", 350, 0, 7000);
        TH1D* h1_invm_re = new TH1D("h1invmre", "h1_invm_re", 350, 0, 7000);
        TH1D* h1_trigPhi_s = new TH1D("h1trigphis", "h1_trigPhi_s", 400, -4, 4);
        TH1D* h1_trigPhi_re = new TH1D("h1trigphire", "h1_trigPhi_re", 400, -4, 4);
        TH1D* h1_dzVTX_s = new TH1D("h1dzvtxs", "h1_dzVTX_s", 100, -0.05, 0.05);
        TH1D* h1_dzVTX_re = new TH1D("h1dzvtxre", "h1_dzVTX_re", 100, -0.05, 0.05);
        TH1D* h1_Full5x5sigma_s = new TH1D("h1full5x5sigmas", "h1_Full5x5sigma_s", 80, 0, 0.08);
        TH1D* h1_Full5x5sigma_re = new TH1D("h1full5x5sigmare", "h1_Full5x5sigma_re", 80, 0, 0.08);

        Bool_t AllOkVec[12];
        Bool_t AllOkHist[6];
        Bool_t AllOkV = kTRUE, AllOkH = kTRUE, AllOk = kTRUE;

        Int_t NEvents_s = ch_s->GetEntries();
        Int_t NEvents_re = ch_re->GetEntries();
        cout << "\tTotal selected Events: " << NEvents_s << endl;
        cout << "\tTotal reselected Events: " << NEvents_re << endl;
        if(NEvents_s!=NEvents_re)
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
                EE_re->GetEvent(i);

                AllOkVec[0] = AreValuesTheSame(EE_s->evtNum, EE_re->evtNum, "event number"); if (AllOkVec[0]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[1] = AreValuesTheSame(EE_s->isHardProcess, EE_re->isHardProcess, "isHardProcess"); if (AllOkVec[1]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[2] = AreValuesTheSame(EE_s->nPileUp, EE_re->nPileUp, "nPileUp"); if (AllOkVec[2]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[3] = AreValuesTheSame(EE_s->GENEvt_weight, EE_re->GENEvt_weight, "event weight"); if (AllOkVec[3]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[4] = AreValuesTheSame(EE_s->HLT_ntrig, EE_re->HLT_ntrig, "HLT_ntrig"); if (AllOkVec[4]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[5] = AreVectorsTheSame(EE_s->Electron_pT, EE_re->Electron_pT, "pT"); if (AllOkVec[5]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[6] = AreVectorsTheSame(EE_s->HLT_trigPhi, EE_re->HLT_trigPhi, "HLT_trigPhi"); if (AllOkVec[6]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[7] = AreVectorsTheSame(EE_s->Electron_charge, EE_re->Electron_charge, "charge"); if (AllOkVec[7]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[8] = AreVectorsTheSame(EE_s->Electron_gsfPx, EE_re->Electron_gsfPx, "gsfPx"); if (AllOkVec[8]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[9] = AreVectorsTheSame(EE_s->Electron_phiSC, EE_re->Electron_phiSC, "phiSC"); if (AllOkVec[9]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[10] = AreVectorsTheSame(EE_s->HLT_trigEta, EE_re->HLT_trigEta, "HLT_trigEta"); if (AllOkVec[10]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[11] = AreVectorsTheSame(EE_s->Electron_E55, EE_re->Electron_E55, "E55"); if (AllOkVec[11]==kFALSE && BreakIfProblem==kTRUE) break;

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
//                Double_t TotWeight_re = MuMu_re->GENEvt_weight;
//                TotWeight_s = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_s->GENEvt_weight;
//                TotWeight_re = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_re->GENEvt_weight;

                // -- Histogram filling -- //
                for (Int_t j=0; j<2; j++)
                {
                    h1_pT_s->Fill(EE_s->Electron_pT->at(j), EE_s->GENEvt_weight);
                    h1_pT_re->Fill(EE_re->Electron_pT->at(j), EE_re->GENEvt_weight);
                    h1_eta_s->Fill(EE_s->Electron_eta->at(j), EE_s->GENEvt_weight);
                    h1_eta_re->Fill(EE_re->Electron_eta->at(j), EE_re->GENEvt_weight);
                    h1_dzVTX_s->Fill(EE_s->Electron_dzVTX->at(j), EE_s->GENEvt_weight);
                    h1_dzVTX_re->Fill(EE_re->Electron_dzVTX->at(j), EE_re->GENEvt_weight);
                    h1_Full5x5sigma_s->Fill(EE_s->Electron_Full5x5_SigmaIEtaIEta->at(j), EE_s->GENEvt_weight);
                    h1_Full5x5sigma_re->Fill(EE_re->Electron_Full5x5_SigmaIEtaIEta->at(j), EE_re->GENEvt_weight);
                }
                h1_invm_s->Fill(EE_s->Electron_InvM, EE_s->GENEvt_weight);
                h1_invm_re->Fill(EE_re->Electron_InvM, EE_re->GENEvt_weight);
                for (UInt_t j=0; j<EE_s->HLT_trigPhi->size(); j++)
                {
                    h1_trigPhi_s->Fill(EE_s->HLT_trigPhi->at(j), EE_s->GENEvt_weight);
                    h1_trigPhi_re->Fill(EE_re->HLT_trigPhi->at(j), EE_re->GENEvt_weight);
                }

                bar.Draw(i);
            } //End of event iteration
        } //End of if( same event numbers )

        if ( AllOk==kTRUE && AllOkV==kTRUE ) cout << "\tNo problems found so far." << endl;
        else cout << "\tProblems were detected." << endl;
        cout << "\tChecking histograms:  ";

        // -- Histogram checking -- //
        AllOkHist[0] = AreHistogramsTheSame(h1_pT_s, h1_pT_re, "pT", BreakIfProblem);
        AllOkHist[1] = AreHistogramsTheSame(h1_eta_s, h1_eta_re, "eta", BreakIfProblem);
        AllOkHist[2] = AreHistogramsTheSame(h1_invm_s, h1_invm_re, "InvM", BreakIfProblem);
        AllOkHist[3] = AreHistogramsTheSame(h1_trigPhi_s, h1_trigPhi_re, "HLT_trigPhi", BreakIfProblem);
        AllOkHist[4] = AreHistogramsTheSame(h1_dzVTX_s, h1_dzVTX_re, "dzVTX", BreakIfProblem);
        AllOkHist[5] = AreHistogramsTheSame(h1_Full5x5sigma_s, h1_Full5x5sigma_re, "Full5x5_SigmaIEtaIEta", BreakIfProblem);

        for ( Int_t j=0; j<6; j++ )
        {
            if ( AllOkHist[j]==kFALSE ) AllOkH = kFALSE;
        }

        // -- Drawing histos -- //
        if ( DrawHistos==kTRUE )
        {
            TCanvas* c_pT = new TCanvas ("pT", "pT", 1000, 1000);
            TCanvas* c_eta = new TCanvas ("eta", "eta", 1000, 1000);
            TCanvas* c_invm = new TCanvas ("invm", "invm", 1000, 1000);
            TCanvas* c_trigPhi = new TCanvas ("trigPhi", "trigPhi", 1000, 1000);
            TCanvas* c_dzVTX = new TCanvas ("dzVTX", "dzVTX", 1000, 1000);
            TCanvas* c_Full5x5sigma = new TCanvas ("Full5x5_SigmaIEtaIEta", "Full5x5_SigmaIEtaIEta", 1000, 1000);
            TLegend* leg = new TLegend (0.1, 0.8, 0.3, 0.9);

            h1_pT_s->SetLineColor(kBlue);
            h1_pT_s->SetLineWidth(2);
            h1_pT_re->SetLineColor(kRed);
            h1_eta_s->SetLineColor(kBlue);
            h1_eta_s->SetLineWidth(2);
            h1_eta_re->SetLineColor(kRed);
            h1_invm_s->SetLineColor(kBlue);
            h1_invm_s->SetLineWidth(2);
            h1_invm_re->SetLineColor(kRed);
            h1_trigPhi_s->SetLineColor(kBlue);
            h1_trigPhi_s->SetLineWidth(2);
            h1_trigPhi_re->SetLineColor(kRed);
            h1_dzVTX_s->SetLineColor(kBlue);
            h1_dzVTX_s->SetLineWidth(2);
            h1_dzVTX_re->SetLineColor(kRed);
            h1_Full5x5sigma_s->SetLineColor(kBlue);
            h1_Full5x5sigma_s->SetLineWidth(2);
            h1_Full5x5sigma_re->SetLineColor(kRed);
            leg->AddEntry(h1_pT_s, "Selected EE", "l");
            leg->AddEntry(h1_pT_re, "Reselected EE", "l");

            c_pT->cd(); h1_pT_s->Draw(); h1_pT_re->Draw("SAME"); leg->Draw(); c_pT->Update();
            c_eta->cd(); h1_eta_s->Draw(); h1_eta_re->Draw("SAME"); leg->Draw(); c_eta->Update();
            c_invm->cd(); h1_invm_s->Draw(); h1_invm_re->Draw("SAME"); leg->Draw(); c_invm->Update();
            c_trigPhi->cd(); h1_trigPhi_s->Draw(); h1_trigPhi_re->Draw("SAME"); leg->Draw(); c_trigPhi->Update();
            c_dzVTX->cd(); h1_dzVTX_s->Draw(); h1_dzVTX_re->Draw("SAME"); leg->Draw(); c_dzVTX->Update();
            c_Full5x5sigma->cd(); h1_Full5x5sigma_s->Draw(); h1_Full5x5sigma_re->Draw("SAME"); leg->Draw(); c_Full5x5sigma->Update();
        }

        if ( AllOkH==kTRUE ) cout << "All bin values match." << endl;
        else cout << "Problems were detected." << endl;

        if ( AllOk==kTRUE && AllOkV==kTRUE && AllOkH==kTRUE ) cout << "Selected and reselected events match perfectly! Hooray!" << endl;
    } //end of i_tup iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds." << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

}// end of CheckLongSelectedMuMu

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
                cout << vec1->at(i) << "    " << vec2->at(i) << endl;
                areThey = kFALSE;
            }
        }
    }
    return areThey;
}// end of CheckLongSelectedEE


/// -------------------------------- Muon Channel ------------------------------------ ///
void CheckLongSelectedMuMu ( Bool_t DrawHistos, Bool_t BreakIfProblem )
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
        TChain *ch_re = new TChain("DYTree");
        ch_s->Add("/media/sf_DATA/test/LongSelectedMuMu_ZToMuMu_M4500to6000_4.root");
        ch_re->Add("/media/sf_DATA/test/ReSelectedMuMu_ZToMuMu_M4500to6000_4.root");

        LongSelectedMuMu_t *MuMu_s = new LongSelectedMuMu_t();
        LongSelectedMuMu_t *MuMu_re = new LongSelectedMuMu_t();
        MuMu_s->CreateFromChain(ch_s);
        MuMu_re->CreateFromChain(ch_re);

        TH1D* h1_pT_s = new TH1D("h1pts", "h1_pT_s", 400, 0, 4000);
        TH1D* h1_pT_re = new TH1D("h1ptre", "h1_pT_re", 400, 0, 4000);
        TH1D* h1_TuneP_pT_s = new TH1D("h1tuneppts", "h1_TuneP_pT_s", 400, 0, 4000);
        TH1D* h1_TuneP_pT_re = new TH1D("h1tunepptre", "h1_TuneP_pT_re", 400, 0, 4000);
        TH1D* h1_eta_s = new TH1D("h1etas", "h1_eta_s", 120, -3, 3);
        TH1D* h1_eta_re = new TH1D("h1etare", "h1_eta_re", 120, -3, 3);
        TH1D* h1_TuneP_eta_s = new TH1D("h1tunepetas", "h1_TuneP_eta_s", 120, -3, 3);
        TH1D* h1_TuneP_eta_re = new TH1D("h1tunepetare", "h1_TuneP_eta_re", 120, -3, 3);
        TH1D* h1_invm_s = new TH1D("h1invms", "h1_invm_s", 350, 0, 7000);
        TH1D* h1_invm_re = new TH1D("h1invmre", "h1_invm_re", 350, 0, 7000);
        TH1D* h1_trigPhi_s = new TH1D("h1trigphis", "h1_trigPhi_s", 400, -4, 4);
        TH1D* h1_trigPhi_re = new TH1D("h1trigphire", "h1_trigPhi_re", 400, -4, 4);
        TH1D* h1_dzVTX_s = new TH1D("h1dzvtxs", "h1_dzVTX_s", 100, -0.05, 0.05);
        TH1D* h1_dzVTX_re = new TH1D("h1dzvtxre", "h1_dzVTX_re", 100, -0.05, 0.05);
        TH1D* h1_vtxTrkProb_s = new TH1D("h1vtxtrkprobs", "h1_vtxTrkProb_s", 55, 0, 1.1);
        TH1D* h1_vtxTrkProb_re = new TH1D("h1vtxtrkprobre", "h1_vtxTrkProb_re", 55, 0, 1.1);

        Bool_t AllOkVec[12];
        Bool_t AllOkHist[8];
        Bool_t AllOkV = kTRUE, AllOkH = kTRUE, AllOk = kTRUE;

        Int_t NEvents_s = ch_s->GetEntries();
        Int_t NEvents_re = ch_re->GetEntries();
        cout << "\tTotal selected Events: " << NEvents_s << endl;
        cout << "\tTotal reselected Events: " << NEvents_re << endl;
        if(NEvents_s!=NEvents_re)
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
                MuMu_re->GetEvent(i);

                AllOkVec[0] = AreValuesTheSame(MuMu_s->evtNum, MuMu_re->evtNum, "event number"); if (AllOkVec[0]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[1] = AreValuesTheSame(MuMu_s->isHardProcess, MuMu_re->isHardProcess, "isHardProcess"); if (AllOkVec[1]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[2] = AreValuesTheSame(MuMu_s->nPileUp, MuMu_re->nPileUp, "nPileUp"); if (AllOkVec[2]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[3] = AreValuesTheSame(MuMu_s->GENEvt_weight, MuMu_re->GENEvt_weight, "event weight"); if (AllOkVec[3]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[4] = AreValuesTheSame(MuMu_s->HLT_ntrig, MuMu_re->HLT_ntrig, "HLT_ntrig"); if (AllOkVec[4]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[5] = AreVectorsTheSame(MuMu_s->Muon_pT, MuMu_re->Muon_pT, "pT"); if (AllOkVec[5]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[6] = AreVectorsTheSame(MuMu_s->HLT_trigPhi, MuMu_re->HLT_trigPhi, "HLT_trigPhi"); if (AllOkVec[6]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[7] = AreVectorsTheSame(MuMu_s->Muon_charge, MuMu_re->Muon_charge, "charge"); if (AllOkVec[7]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[8] = AreVectorsTheSame(MuMu_s->Muon_Px, MuMu_re->Muon_Px, "px"); if (AllOkVec[8]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[9] = AreVectorsTheSame(MuMu_s->Muon_TuneP_phi, MuMu_re->Muon_TuneP_phi, "TuneP_phi"); if (AllOkVec[9]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[10] = AreVectorsTheSame(MuMu_s->HLT_trigEta, MuMu_re->HLT_trigEta, "HLT_trigEta"); if (AllOkVec[10]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[11] = AreVectorsTheSame(MuMu_s->CosAngle, MuMu_re->CosAngle, "CosAngle"); if (AllOkVec[11]==kFALSE && BreakIfProblem==kTRUE) break;

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
//                Double_t TotWeight_re = MuMu_re->GENEvt_weight;
//                TotWeight_s = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_s->GENEvt_weight;
//                TotWeight_re = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_re->GENEvt_weight;

                // -- Histogram filling -- //
                for (Int_t j=0; j<2; j++)
                {
                    h1_pT_s->Fill(MuMu_s->Muon_pT->at(j), MuMu_s->GENEvt_weight);
                    h1_pT_re->Fill(MuMu_re->Muon_pT->at(j), MuMu_re->GENEvt_weight);
                    h1_TuneP_pT_s->Fill(MuMu_s->Muon_TuneP_pT->at(j), MuMu_s->GENEvt_weight);
                    h1_TuneP_pT_re->Fill(MuMu_re->Muon_TuneP_pT->at(j), MuMu_re->GENEvt_weight);
                    h1_eta_s->Fill(MuMu_s->Muon_eta->at(j), MuMu_s->GENEvt_weight);
                    h1_eta_re->Fill(MuMu_re->Muon_eta->at(j), MuMu_re->GENEvt_weight);
                    h1_TuneP_eta_s->Fill(MuMu_s->Muon_TuneP_eta->at(j), MuMu_s->GENEvt_weight);
                    h1_TuneP_eta_re->Fill(MuMu_re->Muon_TuneP_eta->at(j), MuMu_re->GENEvt_weight);
                    h1_dzVTX_s->Fill(MuMu_s->Muon_dzVTX->at(j), MuMu_s->GENEvt_weight);
                    h1_dzVTX_re->Fill(MuMu_re->Muon_dzVTX->at(j), MuMu_re->GENEvt_weight);
                }
                h1_invm_s->Fill(MuMu_s->Muon_InvM, MuMu_s->GENEvt_weight);
                h1_invm_re->Fill(MuMu_re->Muon_InvM, MuMu_re->GENEvt_weight);
                h1_vtxTrkProb_s->Fill(MuMu_s->vtxTrkProb->at(0), MuMu_s->GENEvt_weight);
                h1_vtxTrkProb_re->Fill(MuMu_re->vtxTrkProb->at(0), MuMu_re->GENEvt_weight);
                for (UInt_t j=0; j<MuMu_s->HLT_trigPhi->size(); j++)
                {
                    h1_trigPhi_s->Fill(MuMu_s->HLT_trigPhi->at(j), MuMu_s->GENEvt_weight);
                    h1_trigPhi_re->Fill(MuMu_re->HLT_trigPhi->at(j), MuMu_re->GENEvt_weight);
                }

                bar.Draw(i);
            } //End of event iteration
        } //End of if( same event numbers )

        if ( AllOk==kTRUE && AllOkV==kTRUE ) cout << "\tNo problems found so far." << endl;
        else cout << "\tProblems were detected." << endl;
        cout << "\tChecking histograms:  ";

        // -- Histogram checking -- //
        AllOkHist[0] = AreHistogramsTheSame(h1_pT_s, h1_pT_re, "pT", BreakIfProblem);
        AllOkHist[1] = AreHistogramsTheSame(h1_TuneP_pT_s, h1_TuneP_pT_re, "TuneP_pT", BreakIfProblem);
        AllOkHist[2] = AreHistogramsTheSame(h1_eta_s, h1_eta_re, "eta", BreakIfProblem);
        AllOkHist[3] = AreHistogramsTheSame(h1_TuneP_eta_s, h1_TuneP_eta_re, "TuneP_eta", BreakIfProblem);
        AllOkHist[4] = AreHistogramsTheSame(h1_invm_s, h1_invm_re, "InvM", BreakIfProblem);
        AllOkHist[5] = AreHistogramsTheSame(h1_trigPhi_s, h1_trigPhi_re, "HLT_trigPhi", BreakIfProblem);
        AllOkHist[6] = AreHistogramsTheSame(h1_dzVTX_s, h1_dzVTX_re, "dzVTX", BreakIfProblem);
        AllOkHist[7] = AreHistogramsTheSame(h1_vtxTrkProb_s, h1_vtxTrkProb_re, "vtxTrkProb", BreakIfProblem);

        for ( Int_t j=0; j<8; j++ )
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
            TCanvas* c_invm = new TCanvas ("invm", "invm", 1000, 1000);
            TCanvas* c_trigPhi = new TCanvas ("trigPhi", "trigPhi", 1000, 1000);
            TCanvas* c_dzVTX = new TCanvas ("dzVTX", "dzVTX", 1000, 1000);
            TCanvas* c_vtxTrkProb = new TCanvas ("vtxTrkProb", "vtxTrkProb", 1000, 1000);
            TLegend* leg = new TLegend (0.1, 0.8, 0.3, 0.9);

            h1_pT_s->SetLineColor(kBlue);
            h1_pT_s->SetLineWidth(2);
            h1_pT_re->SetLineColor(kRed);
            h1_TuneP_pT_s->SetLineColor(kBlue);
            h1_TuneP_pT_s->SetLineWidth(2);
            h1_TuneP_pT_re->SetLineColor(kRed);
            h1_eta_s->SetLineColor(kBlue);
            h1_eta_s->SetLineWidth(2);
            h1_eta_re->SetLineColor(kRed);
            h1_TuneP_eta_s->SetLineColor(kBlue);
            h1_TuneP_eta_s->SetLineWidth(2);
            h1_TuneP_eta_re->SetLineColor(kRed);
            h1_invm_s->SetLineColor(kBlue);
            h1_invm_s->SetLineWidth(2);
            h1_invm_re->SetLineColor(kRed);
            h1_trigPhi_s->SetLineColor(kBlue);
            h1_trigPhi_s->SetLineWidth(2);
            h1_trigPhi_re->SetLineColor(kRed);
            h1_dzVTX_s->SetLineColor(kBlue);
            h1_dzVTX_s->SetLineWidth(2);
            h1_dzVTX_re->SetLineColor(kRed);
            h1_vtxTrkProb_s->SetLineColor(kBlue);
            h1_vtxTrkProb_s->SetLineWidth(2);
            h1_vtxTrkProb_re->SetLineColor(kRed);
            leg->AddEntry(h1_pT_s, "Selected MuMu", "l");
            leg->AddEntry(h1_pT_re, "Reselected MuMu", "l");

            c_pT->cd(); h1_pT_s->Draw(); h1_pT_re->Draw("SAME"); leg->Draw(); c_pT->Update();
            c_TuneP_pT->cd(); h1_TuneP_pT_s->Draw(); h1_TuneP_pT_re->Draw("SAME"); leg->Draw(); c_TuneP_pT->Update();
            c_eta->cd(); h1_eta_s->Draw(); h1_eta_re->Draw("SAME"); leg->Draw(); c_eta->Update();
            c_TuneP_eta->cd(); h1_TuneP_eta_s->Draw(); h1_TuneP_eta_re->Draw("SAME"); leg->Draw(); c_TuneP_eta->Update();
            c_invm->cd(); h1_invm_s->Draw(); h1_invm_re->Draw("SAME"); leg->Draw(); c_invm->Update();
            c_trigPhi->cd(); h1_trigPhi_s->Draw(); h1_trigPhi_re->Draw("SAME"); leg->Draw(); c_trigPhi->Update();
            c_dzVTX->cd(); h1_dzVTX_s->Draw(); h1_dzVTX_re->Draw("SAME"); leg->Draw(); c_dzVTX->Update();
            c_vtxTrkProb->cd(); h1_vtxTrkProb_s->Draw(); h1_vtxTrkProb_re->Draw("SAME"); leg->Draw(); c_vtxTrkProb->Update();
        }

        if ( AllOkH==kTRUE ) cout << "All bin values match." << endl;
        else cout << "Problems were detected." << endl;

        if ( AllOk==kTRUE && AllOkV==kTRUE && AllOkH==kTRUE ) cout << "Selected and reselected events match perfectly! Hooray!" << endl;
    } //end of i_tup iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds." << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

}// end of CheckLongSelectedMuMu


/// --------------------------------- EMu events --------------------------------- ///
void CheckLongSelectedEMu ( Bool_t DrawHistos, Bool_t BreakIfProblem )
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
        TChain *ch_re = new TChain("DYTree");
        ch_s->Add("/media/sf_DATA/test/LongSelectedEMu_WW_34.root");
        ch_re->Add("/media/sf_DATA/test/ReSelectedEMu_WW_34.root");

        LongSelectedEMu_t *EMu_s = new LongSelectedEMu_t();
        LongSelectedEMu_t *EMu_re = new LongSelectedEMu_t();
        EMu_s->CreateFromChain(ch_s);
        EMu_re->CreateFromChain(ch_re);

        TH1D* h1_pT_ele_s = new TH1D("h1pteles", "h1_pT_ele_s", 400, 0, 4000);
        TH1D* h1_pT_ele_re = new TH1D("h1ptelere", "h1_pT_ele_re", 400, 0, 4000);
        TH1D* h1_eta_ele_s = new TH1D("h1etaeles", "h1_eta_ele_s", 120, -3, 3);
        TH1D* h1_eta_ele_re = new TH1D("h1etaelere", "h1_eta_ele_re", 120, -3, 3);
        TH1D* h1_pT_mu_s = new TH1D("h1ptmus", "h1_pT_mu_s", 400, 0, 4000);
        TH1D* h1_pT_mu_re = new TH1D("h1ptmure", "h1_pT_mu_re", 400, 0, 4000);
        TH1D* h1_eta_mu_s = new TH1D("h1etamus", "h1_eta_mu_s", 120, -3, 3);
        TH1D* h1_eta_mu_re = new TH1D("h1etamure", "h1_eta_mu_re", 120, -3, 3);
        TH1D* h1_invm_s = new TH1D("h1invms", "h1_invm_s", 350, 0, 7000);
        TH1D* h1_invm_re = new TH1D("h1invmre", "h1_invm_re", 350, 0, 7000);
        TH1D* h1_trigPhi_s = new TH1D("h1trigphis", "h1_trigPhi_s", 400, -4, 4);
        TH1D* h1_trigPhi_re = new TH1D("h1trigphire", "h1_trigPhi_re", 400, -4, 4);

        Bool_t AllOkVec[16];
        Bool_t AllOkHist[6];
        Bool_t AllOkV = kTRUE, AllOkH = kTRUE, AllOk = kTRUE;

        Int_t NEvents_s = ch_s->GetEntries();
        Int_t NEvents_re = ch_re->GetEntries();
        cout << "\tTotal selected Events: " << NEvents_s << endl;
        cout << "\tTotal reselected Events: " << NEvents_re << endl;
        if(NEvents_s!=NEvents_re)
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
                EMu_re->GetEvent(i);

                AllOkVec[0] = AreValuesTheSame(EMu_s->evtNum, EMu_re->evtNum, "event number"); if (AllOkVec[0]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[1] = AreValuesTheSame(EMu_s->isHardProcess, EMu_re->isHardProcess, "isHardProcess"); if (AllOkVec[1]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[2] = AreValuesTheSame(EMu_s->nPileUp, EMu_re->nPileUp, "nPileUp"); if (AllOkVec[2]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[3] = AreValuesTheSame(EMu_s->GENEvt_weight, EMu_re->GENEvt_weight, "event weight"); if (AllOkVec[3]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[4] = AreValuesTheSame(EMu_s->HLT_ntrig, EMu_re->HLT_ntrig, "HLT_ntrig"); if (AllOkVec[4]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[5] = AreValuesTheSame(EMu_s->Electron_pT, EMu_re->Electron_pT, "Electron pT"); if (AllOkVec[5]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[6] = AreVectorsTheSame(EMu_s->HLT_trigPhi, EMu_re->HLT_trigPhi, "HLT_trigPhi"); if (AllOkVec[6]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[7] = AreValuesTheSame(EMu_s->Electron_charge, EMu_re->Electron_charge, "Electron charge"); if (AllOkVec[7]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[8] = AreValuesTheSame(EMu_s->Electron_gsfPx, EMu_re->Electron_gsfPx, "Electron gsfPx"); if (AllOkVec[8]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[9] = AreValuesTheSame(EMu_s->Electron_phiSC, EMu_re->Electron_phiSC, "Electron phiSC"); if (AllOkVec[9]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[10] = AreVectorsTheSame(EMu_s->HLT_trigEta, EMu_re->HLT_trigEta, "HLT_trigEta"); if (AllOkVec[10]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[11] = AreValuesTheSame(EMu_s->Electron_E55, EMu_re->Electron_E55, "Electron E55"); if (AllOkVec[11]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[12] = AreValuesTheSame(EMu_s->Muon_pT, EMu_re->Muon_pT, "Muon pT"); if (AllOkVec[12]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[13] = AreValuesTheSame(EMu_s->Muon_charge, EMu_re->Muon_charge, "Muon charge"); if (AllOkVec[13]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[14] = AreValuesTheSame(EMu_s->Muon_Px, EMu_re->Muon_Px, "Muon px"); if (AllOkVec[14]==kFALSE && BreakIfProblem==kTRUE) break;
                AllOkVec[15] = AreValuesTheSame(EMu_s->Muon_TuneP_phi, EMu_re->Muon_TuneP_phi, "Muon TuneP_phi"); if (AllOkVec[15]==kFALSE && BreakIfProblem==kTRUE) break;

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
//                Double_t TotWeight_re = MuMu_re->GENEvt_weight;
//                TotWeight_s = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_s->GENEvt_weight;
//                TotWeight_re = (L*Xsec[i_tup]/nEvents[i_tup])*MuMu_re->GENEvt_weight;

                // -- Histogram filling -- //
                h1_pT_ele_s->Fill(EMu_s->Electron_pT, EMu_s->GENEvt_weight);
                h1_pT_ele_re->Fill(EMu_re->Electron_pT, EMu_re->GENEvt_weight);
                h1_eta_ele_s->Fill(EMu_s->Electron_eta, EMu_s->GENEvt_weight);
                h1_eta_ele_re->Fill(EMu_re->Electron_eta, EMu_re->GENEvt_weight);
                h1_pT_mu_s->Fill(EMu_s->Muon_pT, EMu_s->GENEvt_weight);
                h1_pT_mu_re->Fill(EMu_re->Muon_pT, EMu_re->GENEvt_weight);
                h1_eta_mu_s->Fill(EMu_s->Muon_eta, EMu_s->GENEvt_weight);
                h1_eta_mu_re->Fill(EMu_re->Muon_eta, EMu_re->GENEvt_weight);
                h1_invm_s->Fill(EMu_s->EMu_InvM, EMu_s->GENEvt_weight);
                h1_invm_re->Fill(EMu_re->EMu_InvM, EMu_re->GENEvt_weight);
                for (UInt_t j=0; j<EMu_s->HLT_trigPhi->size(); j++)
                {
                    h1_trigPhi_s->Fill(EMu_s->HLT_trigPhi->at(j), EMu_s->GENEvt_weight);
                    h1_trigPhi_re->Fill(EMu_re->HLT_trigPhi->at(j), EMu_re->GENEvt_weight);
                }

                bar.Draw(i);
            } //End of event iteration
        } //End of if( same event numbers )

        if ( AllOk==kTRUE && AllOkV==kTRUE ) cout << "\tNo problems found so far." << endl;
        else cout << "\tProblems were detected." << endl;
        cout << "\tChecking histograms:  ";

        // -- Histogram checking -- //
        AllOkHist[0] = AreHistogramsTheSame(h1_pT_ele_s, h1_pT_ele_re, "Electron pT", BreakIfProblem);
        AllOkHist[1] = AreHistogramsTheSame(h1_eta_ele_s, h1_eta_ele_re, "Electron eta", BreakIfProblem);
        AllOkHist[2] = AreHistogramsTheSame(h1_pT_mu_s, h1_pT_mu_re, "Muon pT", BreakIfProblem);
        AllOkHist[3] = AreHistogramsTheSame(h1_eta_mu_s, h1_eta_mu_re, "Muon eta", BreakIfProblem);
        AllOkHist[4] = AreHistogramsTheSame(h1_invm_s, h1_invm_re, "InvM", BreakIfProblem);
        AllOkHist[5] = AreHistogramsTheSame(h1_trigPhi_s, h1_trigPhi_re, "HLT_trigPhi", BreakIfProblem);

        for ( Int_t j=0; j<6; j++ )
        {
            if ( AllOkHist[j]==kFALSE ) AllOkH = kFALSE;
        }

        // -- Drawing histos -- //
        if ( DrawHistos==kTRUE )
        {
            TCanvas* c_pT_ele = new TCanvas ("Electron pT", "Electron pT", 1000, 1000);
            TCanvas* c_eta_ele = new TCanvas ("Electron eta", "Electron eta", 1000, 1000);
            TCanvas* c_pT_mu = new TCanvas ("Muon pT", "Muon pT", 1000, 1000);
            TCanvas* c_eta_mu = new TCanvas ("Muon eta", "Muon eta", 1000, 1000);
            TCanvas* c_invm = new TCanvas ("invm", "invm", 1000, 1000);
            TCanvas* c_trigPhi = new TCanvas ("trigPhi", "trigPhi", 1000, 1000);
            TLegend* leg = new TLegend (0.1, 0.8, 0.3, 0.9);

            h1_pT_ele_s->SetLineColor(kBlue);
            h1_pT_ele_s->SetLineWidth(2);
            h1_pT_ele_re->SetLineColor(kRed);
            h1_eta_ele_s->SetLineColor(kBlue);
            h1_eta_ele_s->SetLineWidth(2);
            h1_eta_ele_re->SetLineColor(kRed);
            h1_pT_mu_s->SetLineColor(kBlue);
            h1_pT_mu_s->SetLineWidth(2);
            h1_pT_mu_re->SetLineColor(kRed);
            h1_eta_mu_s->SetLineColor(kBlue);
            h1_eta_mu_s->SetLineWidth(2);
            h1_eta_mu_re->SetLineColor(kRed);
            h1_invm_s->SetLineColor(kBlue);
            h1_invm_s->SetLineWidth(2);
            h1_invm_re->SetLineColor(kRed);
            h1_trigPhi_s->SetLineColor(kBlue);
            h1_trigPhi_s->SetLineWidth(2);
            h1_trigPhi_re->SetLineColor(kRed);
            leg->AddEntry(h1_pT_ele_s, "Selected EMu", "l");
            leg->AddEntry(h1_pT_ele_re, "Reselected EMu", "l");

            c_pT_ele->cd(); h1_pT_ele_s->Draw(); h1_pT_ele_re->Draw("SAME"); leg->Draw(); c_pT_ele->Update();
            c_eta_ele->cd(); h1_eta_ele_s->Draw(); h1_eta_ele_re->Draw("SAME"); leg->Draw(); c_eta_ele->Update();
            c_pT_mu->cd(); h1_pT_mu_s->Draw(); h1_pT_mu_re->Draw("SAME"); leg->Draw(); c_pT_mu->Update();
            c_eta_mu->cd(); h1_eta_mu_s->Draw(); h1_eta_mu_re->Draw("SAME"); leg->Draw(); c_eta_mu->Update();
            c_invm->cd(); h1_invm_s->Draw(); h1_invm_re->Draw("SAME"); leg->Draw(); c_invm->Update();
            c_trigPhi->cd(); h1_trigPhi_s->Draw(); h1_trigPhi_re->Draw("SAME"); leg->Draw(); c_trigPhi->Update();
        }

        if ( AllOkH==kTRUE ) cout << "All bin values match." << endl;
        else cout << "Problems were detected." << endl;

        if ( AllOk==kTRUE && AllOkV==kTRUE && AllOkH==kTRUE ) cout << "Selected and reselected events match perfectly! Hooray!" << endl;
    } //end of i_tup iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds." << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

}// end of CheckLongSelectedEMu

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
