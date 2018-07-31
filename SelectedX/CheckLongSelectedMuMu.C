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
static inline Bool_t AreHistogramsTheSame(TH1D *hist1, TH1D *hist2, TString name="", Bool_t setBreak=kFALSE);

void CheckLongSelectedMuMu(Bool_t DrawHistos = kFALSE, Bool_t BreakIfProblem = kFALSE)
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
        ch_s->Add("/media/sf_DATA/LongSelectedMuMu_ZToMuMu_M4500to6000_4.root");
        ch_re->Add("/media/sf_DATA/ReSelectedMuMu_ZToMuMu_M4500to6000_4.root");

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

        Bool_t AllOk = kTRUE;

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

                AllOk = AreValuesTheSame(MuMu_s->evtNum, MuMu_re->evtNum, "event number"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;
                AllOk = AreValuesTheSame(MuMu_s->isHardProcess, MuMu_re->isHardProcess, "isHardProcess"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;
                AllOk = AreValuesTheSame(MuMu_s->nPileUp, MuMu_re->nPileUp, "nPileUp"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;
                AllOk = AreValuesTheSame(MuMu_s->GENEvt_weight, MuMu_re->GENEvt_weight, "event weight"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;
                AllOk = AreValuesTheSame(MuMu_s->HLT_ntrig, MuMu_re->HLT_ntrig, "HLT_ntrig"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;
                AllOk = AreVectorsTheSame(MuMu_s->Muon_pT, MuMu_re->Muon_pT, "pT"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;
                AllOk = AreVectorsTheSame(MuMu_s->HLT_trigPhi, MuMu_re->HLT_trigPhi, "HLT_trigPhi"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;
                AllOk = AreVectorsTheSame(MuMu_s->Muon_charge, MuMu_re->Muon_charge, "charge"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;
                AllOk = AreVectorsTheSame(MuMu_s->Muon_Px, MuMu_re->Muon_Px, "px"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;
                AllOk = AreVectorsTheSame(MuMu_s->Muon_TuneP_phi, MuMu_re->Muon_TuneP_phi, "TuneP_phi"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;
                AllOk = AreVectorsTheSame(MuMu_s->HLT_trigEta, MuMu_re->HLT_trigEta, "HLT_trigEta"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;
                AllOk = AreVectorsTheSame(MuMu_s->CosAngle, MuMu_re->CosAngle, "CosAngle"); if (AllOk==kFALSE && BreakIfProblem==kTRUE) break;

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

        if ( AllOk==kTRUE )
        {
            cout << "\tNo problems found so far." << endl;      
            cout << "\tChecking histograms:  ";
            AllOk = AreHistogramsTheSame(h1_pT_s, h1_pT_re, "pT", BreakIfProblem);
            AllOk = AreHistogramsTheSame(h1_TuneP_pT_s, h1_TuneP_pT_re, "TuneP_pT", BreakIfProblem);
            AllOk = AreHistogramsTheSame(h1_eta_s, h1_eta_re, "eta", BreakIfProblem);
            AllOk = AreHistogramsTheSame(h1_TuneP_eta_s, h1_TuneP_eta_re, "TuneP_eta", BreakIfProblem);
            AllOk = AreHistogramsTheSame(h1_invm_s, h1_invm_re, "InvM", BreakIfProblem);
            AllOk = AreHistogramsTheSame(h1_trigPhi_s, h1_trigPhi_re, "HLT_trigPhi", BreakIfProblem);
            AllOk = AreHistogramsTheSame(h1_dzVTX_s, h1_dzVTX_re, "dzVTX", BreakIfProblem);
            AllOk = AreHistogramsTheSame(h1_vtxTrkProb_s, h1_vtxTrkProb_re, "vtxTrkProb", BreakIfProblem);
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

        if ( AllOk==kTRUE )
        {
            cout << "All bin values match.\nSelected and reselected events match perfectly! Hooray!" << endl;
        }
        else cout << "Problems were detected." << endl;
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
