#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TAttMarker.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TApplication.h>
#include <vector>
#include <TMath.h>
#include <THistPainter.h>
#include <TFormula.h>
#include <iostream>
#include <sstream>

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/FileMgr.h"
#include "./header/PrescaleProvider.h"

PrescaleProvider pp("etc/prescale/triggerData2016");

void E_FR_PrescaleTest (Bool_t DEBUG = kFALSE);
void Mu_PrescaleTest (Bool_t DEBUG = kFALSE);
void EMu_PrescaleTest (Bool_t DEBUG = kFALSE);
Double_t PrescalePhoton (Double_t HLT_pT, Int_t runNo, Int_t lumiSec);
Double_t PrescalePhoton (Int_t trig_fired, Double_t HLT_pT, Int_t runNo, Int_t lumiSec);

void FR_PrescaleTest (TString WhichX = "")
{
    TString whichX = WhichX;
    whichX.ToUpper();
    Int_t Xselected = 0;
    Bool_t DEBUG = kFALSE;
    if (whichX.Contains("DEBUG"))
    {
        DEBUG = kTRUE;
        cout << "**** DEBUG MODE: Running with 100 events only. ****" << endl;
    }
    if (whichX.Contains("EMU"))
    {
        Xselected++;
        cout << "\n*****  EMu_PrescaleTest()  *****" << endl;
        EMu_PrescaleTest(DEBUG); // not ready yet
    }
    else if (whichX.Contains("MU"))
    {
        Xselected++;
        cout << "\n*****  Mu_PrescaleTest()  *****" << endl;
        Mu_PrescaleTest(DEBUG); // not ready yet
    }
    else if (whichX.Contains("E"))
    {
        Xselected++;
        cout << "\n*****  E_FR_PrescaleTest()  *****" << endl;
        E_FR_PrescaleTest(DEBUG); // not ready yet
    }

    if (Xselected == 0) cout << "Wrong arument! \nType in: >> .x FR_PrescaleTest.C+(\"whichX\")" << endl;

} // End of FR_PrescaleTest()

void E_FR_PrescaleTest (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    TH1D *h_HLT_pT_full = new TH1D("h_HLT_pT_full", "", 500, 0, 500);
    TH1D *h_HLT_pT_uncorr_full = new TH1D("h_HLT_pT_uncorr_full", "", 500, 0, 500);

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Electron/";
    TString debug = "";
    if (DEBUG) debug = "_DEBUG";

    for (Process_t pr=_SinglePhoton_B; pr<=_SinglePhoton_H; pr=next(pr))
    {
        Mgr.SetProc(pr);

        // -- Output ROOTFile -- //
        f = new TFile(Dir+"FR_Hist_TEST_E_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root", "RECREATE");

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        // -- Creating Histograms -- //
        TH1D* h_pT_uncorr = new TH1D("h_pT_uncorr", "h_pT_uncorr", 500, 0, 500); h_pT_uncorr->Sumw2();
        TH1D* h_pT = new TH1D("h_pT", "h_pT", 500, 0, 500); h_pT->Sumw2();
        TH1D* h_HLT_pT_uncorr = new TH1D("h_HLT_pT_uncorr", "h_HLT_pT_uncorr", 500, 0, 500); h_HLT_pT_uncorr->Sumw2();
        TH1D* h_HLT_pT = new TH1D("h_HLT_pT", "h_HLT_pT", 500, 0, 500); h_HLT_pT->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<int> *trig_fired = new std::vector<int>;
        std::vector<int> *trig_matched = new std::vector<int>;
        std::vector<double> *trig_pT = new std::vector<double>;
        std::vector<int> *prescale = new std::vector<int>;
        Int_t runNum;
        Int_t lumiBlock;
        Int_t prescale_alt_int;
        Double_t prescale_alt;

        TChain *chain = new TChain("FRTree");

        chain->Add(Dir+"SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("trig_fired", 1);
        chain->SetBranchStatus("trig_matched", 1);
        chain->SetBranchStatus("trig_pT", 1);
        chain->SetBranchStatus("runNum", 1);
        chain->SetBranchStatus("lumiBlock", 1);
        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("trig_fired", &trig_fired);
        chain->SetBranchAddress("trig_matched", &trig_matched);
        chain->SetBranchAddress("trig_pT", &trig_pT);
        chain->SetBranchAddress("runNum", &runNum);
        chain->SetBranchAddress("lumiBlock", &lumiBlock);
        chain->SetBranchAddress("prescale_factor", &prescale);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        if (DEBUG == kTRUE) NEvents = 100;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (DEBUG == kTRUE){
                cout << "\nEvt " << i << endl;
                cout << "nEle = " << p_T->size() << endl;
                cout << "p_T[1] = " << p_T->at(0) << endl;
            }

            if (DEBUG == kTRUE)
            {
                cout << "Triggers:" << endl;
                for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
                {
                    cout << "Photon" << trig_fired->at(i_tr) << "  p_T: " << trig_pT->at(i_tr) << "  matched to " << trig_matched->at(i_tr) << endl;
                }
            }

            // We add prescales here
            /*
            std::vector<Double_t> pTs;
//            if (trig_fired->size() > 1) continue;
            Int_t match_found = 0;
            for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
            {
                if (pTs.size())
                {
                    for (UInt_t i=0; i<pTs.size(); i++)
                    {
                        if (trig_pT->at(i_tr) == pTs[i])
                            match_found = 1;
                        else pTs.push_back(trig_pT->at(i_tr));
                    }
                }
                else pTs.push_back(trig_pT->at(i_tr));
            }
//            if (pTs.size() > 1) continue;
            for (UInt_t i=0; i<pTs.size(); i++)
            {
                prescale_alt = PrescalePhoton(pTs.at(i), runNum, lumiBlock);
                h_HLT_pT_uncorr->Fill(pTs.at(i));
                h_HLT_pT->Fill(pTs.at(i), prescale_alt);
                h_HLT_pT_full->Fill(pTs.at(i), prescale_alt);
            }
            */

            // We do not add prescales here
            std::vector<Double_t> pTs;
            for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
            {
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
            }



//            prescale_alt = PrescalePhoton(trig_pT->at(0), runNum, lumiBlock);
//            h_HLT_pT_uncorr->Fill(trig_pT->at(0));
//            h_HLT_pT->Fill(trig_pT->at(0), prescale_alt);
//            h_HLT_pT_full->Fill(trig_pT->at(0), prescale_alt);

//            for (UInt_t i_ele=0; i_ele<p_T->size(); i_ele++)
//            {
//                if (p_T->at(i_ele) != p_T->at(i_ele))
//                {
//                    cout << p_T->at(i_ele) << endl;
//                    continue;
//                }
//                if (p_T->at(i_ele) <= 25) continue;
//                if (DEBUG == kTRUE) cout << "i_ele = " << i_ele << endl;
/*
                Int_t matched22=0, matched30=0, matched36=0, matched50=0, matched75=0, matched90=0, matched120=0, matched175=0;
                Int_t matched22_wopT=0, matched30_wopT=0, matched36_wopT=0, matched50_wopT=0, matched75_wopT=0, matched90_wopT=0, matched120_wopT=0;
                Int_t i_22=-1, i_30=-1, i_36=-1, i_50=-1, i_75=-1, i_90=-1, i_120=-1, i_175=-1;
                prescale_alt = 1;

                for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
                {
//                    if (trig_pT->at(i_tr) < 22) continue;
//                    if (((UInt_t)(trig_matched->at(i_tr))) == i_ele)
//                    {
                        if (trig_fired->at(i_tr) == 22)
                        {
                            matched22_wopT = 1;
                            if (trig_pT->at(i_tr) > 22 && trig_pT->at(i_tr) <= 30)
                            {
                                i_22 = i_tr;
                                matched22 = 1;
                            }
                        }
                        if (trig_fired->at(i_tr) == 30)
                        {
                            matched30_wopT = 1;
                            if (trig_pT->at(i_tr) > 30 && trig_pT->at(i_tr) <= 175)
                            {
                                i_30 = i_tr;
                                matched30 = 1;
                            }
                        }

                        if (trig_fired->at(i_tr) == 36)
                        {
                            matched36_wopT = 1;
                            if (trig_pT->at(i_tr) > 36 && trig_pT->at(i_tr) <= 50)
                            {
                                i_36 = i_tr;
                                matched36 = 1;
                            }
                        }
                        if (trig_fired->at(i_tr) == 50)
                        {
                            matched50_wopT = 1;
                            if (trig_pT->at(i_tr) > 50 && trig_pT->at(i_tr) <= 75)
                            {
                                i_50 = i_tr;
                                matched50 = 1;
                            }
                        }
                        if (trig_fired->at(i_tr) == 75)
                        {
                            matched75_wopT = 1;
                            if (trig_pT->at(i_tr) > 75 && trig_pT->at(i_tr) <= 90)
                            {
                                i_75 = i_tr;
                                matched75 = 1;
                            }
                        }
                        if (trig_fired->at(i_tr) == 90)
                        {
                            matched90_wopT = 1;
                            if (trig_pT->at(i_tr) > 90 && trig_pT->at(i_tr) <= 120)
                            {
                                i_90 = i_tr;
                                matched90 = 1;
                            }
                        }
                        if (trig_fired->at(i_tr) == 120)
                        {
                            matched120_wopT = 1;
                            if (trig_pT->at(i_tr) > 120 && trig_pT->at(i_tr) <= 175)
                            {
                                i_120 = i_tr;
                                matched120 = 1;
                            }
                        }
                        if (trig_fired->at(i_tr) == 175 && trig_pT->at(i_tr) > 175)
                        {
                            i_175 = i_tr;
                            matched175 = 1;
                        }

//                    }
                }

                if (matched22) // Matched an electron to HLT_Photon22 with 22<HLT_pT<30
                { // get the prescale
                    prescale_alt_int = pp.hltPrescale("HLT_Photon22_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG18", runNum, lumiBlock);
                    prescale_alt = (double)prescale_alt_int;
//                    if (pr == _SinglePhoton_B) prescale_alt *= 1.37;
                    h_HLT_pT->Fill(trig_pT->at(i_22), prescale_alt); // fill HLT_pT histo with HLT object pT using the prescale weight
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_22)); // fill HLT_pT histo with HLT object pT without weights
//                    h_pT->Fill(p_T->at(i_ele), prescale_alt); // fill pT histo with matched electron pT using the prescale weight
//                    h_pT_uncorr->Fill(p_T->at(i_ele)); // fill pT histo with matched electron pT without weights
                    h_HLT_pT_full->Fill(trig_pT->at(i_22), prescale_alt);
                }
                if (matched30) // Matched an electron to HLT_Photon22 with 30<HLT_pT<36
                {
                    prescale_alt_int = pp.hltPrescale("HLT_Photon30_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG26", runNum, lumiBlock);
                    prescale_alt = (double)prescale_alt_int;
//                    if (pr == _SinglePhoton_B) prescale_alt *= 1.37;
                    h_HLT_pT->Fill(trig_pT->at(i_30), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_30));
//                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
//                    h_pT_uncorr->Fill(p_T->at(i_ele));
                    h_HLT_pT_full->Fill(trig_pT->at(i_30), prescale_alt);
                }
                if (matched36)
                {
                    prescale_alt_int = pp.hltPrescale("HLT_Photon36_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG26", runNum, lumiBlock);
                    prescale_alt = (double)prescale_alt_int;
//                    if (pr == _SinglePhoton_B) prescale_alt *= 1.37;
                    h_HLT_pT->Fill(trig_pT->at(i_36), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_36));
//                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
//                    h_pT_uncorr->Fill(p_T->at(i_ele));
                    h_HLT_pT_full->Fill(trig_pT->at(i_36), prescale_alt);
                }
                if (matched50)
                {
                    prescale_alt_int = 0;
                    if (pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon50_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon50_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock);
                    prescale_alt = (double)prescale_alt_int;
//                    if (pr == _SinglePhoton_B) prescale_alt *= 1.37;
                    h_HLT_pT->Fill(trig_pT->at(i_50), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_50));
//                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
//                    h_pT_uncorr->Fill(p_T->at(i_ele));
                    h_HLT_pT_full->Fill(trig_pT->at(i_50), prescale_alt);
                }
                if (matched75)
                {
                    prescale_alt_int = 0;
                    if (pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon75_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon75_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock);
                    prescale_alt = (double)prescale_alt_int;
//                    if (pr == _SinglePhoton_B) prescale_alt *= 1.37;
                    h_HLT_pT->Fill(trig_pT->at(i_75), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_75));
//                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
//                    h_pT_uncorr->Fill(p_T->at(i_ele));
                    h_HLT_pT_full->Fill(trig_pT->at(i_75), prescale_alt);
                }
                if (matched90)
                {
                    prescale_alt_int = 0;
                    if (pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon90_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon90_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock);
                    prescale_alt = (double)prescale_alt_int;
//                    if (pr == _SinglePhoton_B) prescale_alt *= 1.37;
                    h_HLT_pT->Fill(trig_pT->at(i_90), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_90));
//                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
//                    h_pT_uncorr->Fill(p_T->at(i_ele));
                    h_HLT_pT_full->Fill(trig_pT->at(i_90), prescale_alt);
                }
                if (matched120)
                {
                    prescale_alt_int = 0;
                    if (pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon120_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon120_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock);
                    prescale_alt = (double)prescale_alt_int;
//                    if (pr == _SinglePhoton_B) prescale_alt *= 1.37;
                    h_HLT_pT->Fill(trig_pT->at(i_120), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_120));
//                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
//                    h_pT_uncorr->Fill(p_T->at(i_ele));
                    h_HLT_pT_full->Fill(trig_pT->at(i_120), prescale_alt);
                }
                if (matched175)
                {
                    prescale_alt_int = 0;
                    if (pp.l1Prescale("L1_SingleEG30", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG30", runNum, lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG32", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG32", runNum, lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG34", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG34", runNum, lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG36", runNum, lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG38", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG38", runNum, lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock) != 0)
                        prescale_alt_int = pp.hltPrescale("HLT_Photon175_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock);
                    prescale_alt = (double)prescale_alt_int;
                    h_HLT_pT->Fill(trig_pT->at(i_175), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_175));
//                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
//                    h_pT_uncorr->Fill(p_T->at(i_ele));
                    h_HLT_pT_full->Fill(trig_pT->at(i_175), prescale_alt);
                }
//                else continue;
*/
//            }// End of i_ele iteration

            if (DEBUG == kFALSE) bar.Draw(i);

        }// End of event iteration

        f->cd();
        cout << "\tWriting into file...";

        h_pT_uncorr->Write();
        h_pT->Write();
        h_HLT_pT_uncorr->Write();
        h_HLT_pT->Write();
        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _SinglePhoton_B) break;

        f->Close();
        if (!f->IsOpen()) cout << "File " << Dir+"FR_Hist_TEST_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " has been closed successfully.\n" << endl;
        else cout << "FILE " << Dir+"FR_Hist_TEST_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

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

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of E_FR_PrescaleTest()


void Mu_PrescaleTest (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;
    PrescaleProvider pp("etc/prescale/triggerData2016");

    TFile *f;
//    TString Dir = "../";
    TString Dir = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    Int_t n50=0, n2050=0, n2750=0;

    cout << "Run " << 273158 << " LumiS " << 1 << " prescale " << pp.hltPrescale("HLT_Mu20_v", 273158, 1) << " " << pp.l1Prescale("L1_SingleMu22", 273158, 1) << endl;

    for (Process_t pr=_SingleMuon_F; pr<=_SingleMuon_F; pr=next(pr))
    {
        Mgr.SetProc(pr);

        // -- Output ROOTFile -- //
        f = new TFile(Dir+"FR_Hist_TEST_Mu_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root", "RECREATE");

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        // -- Creating Histograms -- //
        TH1D* h_pT_uncorr = new TH1D("h_pT_uncorr", "h_pT_uncorr", 500, 0, 500); h_pT_uncorr->Sumw2();
        TH1D* h_pT = new TH1D("h_pT", "h_pT", 500, 0, 500); h_pT->Sumw2();
        TH1D* h_HLT_pT_uncorr = new TH1D("h_HLT_pT_uncorr", "h_HLT_pT_uncorr", 500, 0, 500); h_HLT_pT_uncorr->Sumw2();
        TH1D* h_HLT_pT = new TH1D("h_HLT_pT", "h_HLT_pT", 500, 0, 500); h_HLT_pT->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<int> *trig_fired = new std::vector<int>;
        std::vector<int> *trig_matched = new std::vector<int>;
        std::vector<double> *trig_pT = new std::vector<double>;
        std::vector<int> *prescale_factor = new std::vector<int>;
        Int_t runNum;
        Int_t lumiBlock;
        Int_t prescale_alt;
        Double_t prescale;

        TChain *chain = new TChain("FRTree");

//        chain->Add(Dir+"SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        chain->Add("~/Desktop/SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("trig_fired", 1);
        chain->SetBranchStatus("trig_matched", 1);
        chain->SetBranchStatus("trig_pT", 1);
        chain->SetBranchStatus("runNum", 1);
        chain->SetBranchStatus("lumiBlock", 1);
        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("trig_fired", &trig_fired);
        chain->SetBranchAddress("trig_matched", &trig_matched);
        chain->SetBranchAddress("trig_pT", &trig_pT);
        chain->SetBranchAddress("prescale_factor", &prescale_factor);
        chain->SetBranchAddress("runNum", &runNum);
        chain->SetBranchAddress("lumiBlock", &lumiBlock);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        if (DEBUG == kTRUE) NEvents = 20;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (DEBUG == kTRUE)
            {
                cout << "\nEvt " << i << endl;
                cout << "nMu = " << p_T->size() << endl;
                cout << "RunNum = " << runNum << "  LumiSection = " << lumiBlock << endl;
            }

            if (p_T->at(0) > p_T->at(1)) h_pT->Fill(p_T->at(0));
            else h_pT->Fill(p_T->at(1));

//            for (UInt_t i_mu=0; i_mu<p_T->size(); i_mu++)
//            {
//                if (p_T->at(i_mu) != p_T->at(i_mu))
//                {
//                    cout << p_T->at(i_mu) << endl;
//                    continue;
//                }
//                if (p_T->at(i_mu) <= 27) continue;
//                if (DEBUG == kTRUE) cout << "i_mu = " << i_mu << endl;

//                Int_t matched20=0, matched27=0, matched0=0, high20=0, high27=0, high50=0;
//                Int_t i_20=-1, i_27=-1, i_0=-1;
//                prescale_alt = 1;

//                for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
//                {
//                    if (trig_pT->at(i_tr) < 20) continue;
//                    if (((UInt_t)(trig_matched->at(i_tr))) == i_mu)
//                    {
//                        if (trig_fired->at(i_tr) == 1/* && trig_pT->at(i_tr) > 50*/)
//                        {
//                            i_0 = i_tr;
//                            matched0 = 1;
//                        }
//                        if (trig_fired->at(i_tr) == 20 && trig_pT->at(i_tr) > 20 && trig_pT->at(i_tr) <= 27)
//                        {
//                            i_20 = i_tr;
//                            matched20 = 1;
//                        }
//                        if (trig_fired->at(i_tr) == 27 && trig_pT->at(i_tr) > 27 && trig_pT->at(i_tr) <= 50)
//                        {
//                            i_27 = i_tr;
//                            matched27 = 1;
//                        }
//                        if (trig_fired->at(i_tr) == 1 && trig_pT->at(i_tr) > 50)
//                            high50++;
//                        if (trig_fired->at(i_tr) == 20 && trig_pT->at(i_tr) > 50)
//                            high20++;
//                        if (trig_fired->at(i_tr) == 27 && trig_pT->at(i_tr) > 50)
//                            high27++;
//                        if (DEBUG == kTRUE)
//                        {
//                            cout << "HLT: " << trig_fired->at(i_tr) << " matched to " << trig_matched->at(i_tr) << " with HLT_pT=" << trig_pT->at(i_tr);
//                            cout << ",  pT=" << p_T->at(trig_matched->at(i_tr)) << endl;
//                        }
//                    }
//                }
//                if (DEBUG == kTRUE)
//                {
//                    cout << "i_0=" << i_0 << "  i_20=" << i_20 << "  i_27=" << i_27 << "\nmatched0=" << matched0 << "  matched20=" << matched20;
//                    cout << "  matched27=" << matched27 << "\nhigh50=" << high50 << "  high20=" << high20 << "  high27=" << high27 << endl;
//                    if (matched20 || high20) cout << "Mu20 prescale: " << pp.hltPrescale("HLT_Mu20_v", runNum, lumiBlock) << " * "
//                                                  << pp.l1Prescale("L1_SingleMu18", runNum, lumiBlock) << endl;
//                    if (matched27 || high27) cout << "Mu27 prescale: " << pp.hltPrescale("HLT_Mu27_v", runNum, lumiBlock) << " * "
//                                                  << pp.l1Prescale("L1_SingleMu22", runNum, lumiBlock) << endl;
//                }
//                if (high50) n50++;
//                if (high20 && high50) n2050++;
//                if (high27 && high50) n2750++;

//                if (matched0 == 1)
//                {
//                    prescale_alt = 1;
//                    h_HLT_pT->Fill(trig_pT->at(i_0), prescale_alt);
//                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_0));
//                    h_pT->Fill(p_T->at(i_mu), prescale_alt);
//                    h_pT_uncorr->Fill(p_T->at(i_mu));
//                }
//                else if (matched20 == 1 && matched0 == 0)
//                {
//                    prescale_alt = pp.hltPrescale("HLT_Mu20_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleMu18", runNum, lumiBlock)/10;
//                    prescale = 1000/*4.44*/;
//                    h_HLT_pT->Fill(trig_pT->at(i_20), prescale_alt);
//                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_20));
//                    h_pT->Fill(p_T->at(i_mu), prescale_alt);
//                    h_pT_uncorr->Fill(p_T->at(i_mu));
//                }
//                else if (matched27 == 1 && matched0 == 0 && matched20 == 0)
//                {
//                    prescale_alt = pp.hltPrescale("HLT_Mu27_v", runNum, lumiBlock)/10;
//                    prescale = 2/*1.061*/;
//                    h_HLT_pT->Fill(trig_pT->at(i_27), prescale_alt);
//                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_27));
//                    h_pT->Fill(p_T->at(i_mu), prescale_alt);
//                    h_pT_uncorr->Fill(p_T->at(i_mu));
//                }
//                else continue;
//
//            }// End of i_mu iteration

            if (DEBUG == kFALSE) bar.Draw(i);

        }// End of event iteration

        f->cd();
        cout << "\tWriting into file...";

        h_pT_uncorr->Write();
        h_pT->Write();
        h_HLT_pT_uncorr->Write();
        h_HLT_pT->Write();
        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _SingleMuon_B) break;

        f->Close();
        if (!f->IsOpen()) cout << "File " << Dir+"FR_Hist_TEST_Mu_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " has been closed successfully.\n" << endl;
        else cout << "FILE " << Dir+"FR_Hist_TEST_Mu_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    cout << "Mu50: " << n50 << endl;
    cout << "Mu50 && Mu20: " << n2050 << endl;
    cout << "Mu50 && Mu27: " << n2750 << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of Mu_PrescaleTest()


void EMu_PrescaleTest (Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;
    PrescaleProvider pp("etc/prescale/triggerData2016");

    TFile *f;
//    TString Dir = "../";
    TString Dir = "/media/sf_DATA/FR/EMu/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    Int_t n50=0, n2050=0, n2750=0;

    cout << "Run " << 273158 << " LumiS " << 1 << " prescale " << pp.hltPrescale("HLT_Mu20_v", 273158, 1) << " " << pp.l1Prescale("L1_SingleMu22", 273158, 1) << endl;

    for (Process_t pr=_SingleMuon_B; pr<=_SingleMuon_H; pr=next(pr))
    {
        Mgr.SetProc(pr);

        // -- Output ROOTFile -- //
        f = new TFile(Dir+"FR_Hist_TEST_EMu_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root", "RECREATE");

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        // -- Creating Histograms -- //
        TH1D* h_pT_uncorr = new TH1D("h_pT_uncorr", "h_pT_uncorr", 500, 0, 500); h_pT_uncorr->Sumw2();
        TH1D* h_pT = new TH1D("h_pT", "h_pT", 500, 0, 500); h_pT->Sumw2();
        TH1D* h_HLT_pT_uncorr = new TH1D("h_HLT_pT_uncorr", "h_HLT_pT_uncorr", 500, 0, 500); h_HLT_pT_uncorr->Sumw2();
        TH1D* h_HLT_pT = new TH1D("h_HLT_pT", "h_HLT_pT", 500, 0, 500); h_HLT_pT->Sumw2();

        Double_t mu_p_T;
        std::vector<int> *trig_fired = new std::vector<int>;
        std::vector<double> *trig_pT = new std::vector<double>;
        std::vector<int> *prescale_factor = new std::vector<int>;
        Int_t runNum;
        Int_t lumiBlock;
        Int_t prescale_alt;
        Double_t prescale;

        TChain *chain = new TChain("FRTree");

        chain->Add(Dir+"SelectedForBKGest_EMu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_EMu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        chain->SetBranchStatus("mu_p_T", 1);
        chain->SetBranchStatus("trig_fired", 1);
        chain->SetBranchStatus("trig_pT", 1);
        chain->SetBranchStatus("runNum", 1);
        chain->SetBranchStatus("lumiBlock", 1);
        chain->SetBranchAddress("mu_p_T", &mu_p_T);
        chain->SetBranchAddress("trig_fired", &trig_fired);
        chain->SetBranchAddress("trig_pT", &trig_pT);
        chain->SetBranchAddress("prescale_factor", &prescale_factor);
        chain->SetBranchAddress("runNum", &runNum);
        chain->SetBranchAddress("lumiBlock", &lumiBlock);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        if (DEBUG == kTRUE) NEvents = 20;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (DEBUG == kTRUE)
            {
                cout << "\nEvt " << i << endl;
                cout << "RunNum = " << runNum << "  LumiSection = " << lumiBlock << endl;
                cout << "nTrig = " << trig_fired->size() << endl;
            }

            if (mu_p_T != mu_p_T)
            {
                cout << mu_p_T << endl;
                continue;
            }
            if (mu_p_T <= 27) continue;

            Int_t matched20=0, matched27=0, matched0=0, high20=0, high27=0, high50=0;
            Int_t i_20=-1, i_27=-1, i_0=-1;
            prescale_alt = 1;

            for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
            {
                if (trig_pT->at(i_tr) < 20) continue;
                {
                    if (trig_fired->at(i_tr) == 1 && trig_pT->at(i_tr) > 50)
                    {
                        i_0 = i_tr;
                        matched0 = 1;
                    }
                    if (trig_fired->at(i_tr) == 20 && trig_pT->at(i_tr) > 20 && trig_pT->at(i_tr) <= 27)
                    {
                        i_20 = i_tr;
                        matched20 = 1;
                    }
                    if (trig_fired->at(i_tr) == 27 && trig_pT->at(i_tr) > 27 && trig_pT->at(i_tr) <= 50)
                    {
                        i_27 = i_tr;
                        matched27 = 1;
                    }
                    if (trig_fired->at(i_tr) == 1 && trig_pT->at(i_tr) > 50)
                        high50++;
                    if (trig_fired->at(i_tr) == 20 && trig_pT->at(i_tr) > 50)
                        high20++;
                    if (trig_fired->at(i_tr) == 27 && trig_pT->at(i_tr) > 50)
                        high27++;
                    if (DEBUG == kTRUE)
                    {
                        cout << "HLT: " << trig_fired->at(i_tr) << " with HLT_pT=" << trig_pT->at(i_tr) << ",  pT=" << mu_p_T << endl;
                    }
                }
            }
            if (DEBUG == kTRUE)
            {
                cout << "i_0=" << i_0 << "  i_20=" << i_20 << "  i_27=" << i_27 << "\nmatched0=" << matched0 << "  matched20=" << matched20;
                cout << "  matched27=" << matched27 << "\nhigh50=" << high50 << "  high20=" << high20 << "  high27=" << high27 << endl;
                if (matched20 || high20) cout << "Mu20 prescale: " << pp.hltPrescale("HLT_Mu20_v", runNum, lumiBlock) << " * "
                                              << pp.l1Prescale("L1_SingleMu18", runNum, lumiBlock) << endl;
                if (matched27 || high27) cout << "Mu27 prescale: " << pp.hltPrescale("HLT_Mu27_v", runNum, lumiBlock) << " * "
                                              << pp.l1Prescale("L1_SingleMu22", runNum, lumiBlock) << endl;
            }
            if (high50) n50++;
            if (high20 && high50) n2050++;
            if (high27 && high50) n2750++;

            if (matched0 == 1)
            {
                prescale_alt = 1;
                h_HLT_pT->Fill(trig_pT->at(i_0), prescale_alt);
                h_HLT_pT_uncorr->Fill(trig_pT->at(i_0));
                h_pT->Fill(mu_p_T, prescale_alt);
                h_pT_uncorr->Fill(mu_p_T);
            }
            else if (matched20 == 1)
            {
                prescale_alt = pp.hltPrescale("HLT_Mu20_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleMu18", runNum, lumiBlock);
                prescale = 1000/*4.44*/;
                h_HLT_pT->Fill(trig_pT->at(i_20), prescale_alt);
                h_HLT_pT_uncorr->Fill(mu_p_T);
                h_pT->Fill(mu_p_T, prescale_alt);
                h_pT_uncorr->Fill(mu_p_T);
            }
            else if (matched27 == 1)
            {
                prescale_alt = pp.hltPrescale("HLT_Mu27_v", runNum, lumiBlock);
                prescale = 2/*1.061*/;
                h_HLT_pT->Fill(trig_pT->at(i_27), prescale_alt);
                h_HLT_pT_uncorr->Fill(trig_pT->at(i_27));
                h_pT->Fill(mu_p_T, prescale_alt);
                h_pT_uncorr->Fill(mu_p_T);
            }
            else continue;

            if (DEBUG == kFALSE) bar.Draw(i);

        }// End of event iteration

        f->cd();
        cout << "\tWriting into file...";

        h_pT_uncorr->Write();
        h_pT->Write();
        h_HLT_pT_uncorr->Write();
        h_HLT_pT->Write();
        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _SingleMuon_B) break;

        f->Close();
        if (!f->IsOpen()) cout << "File " << Dir+"FR_Hist_TEST_EMu_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " has been closed successfully.\n" << endl;
        else cout << "FILE " << Dir+"FR_Hist_TEST_EMu_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    cout << "Mu50: " << n50 << endl;
    cout << "Mu50 && Mu20: " << n2050 << endl;
    cout << "Mu50 && Mu27: " << n2750 << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of EMu_PrescaleTest()

Double_t PrescalePhoton (Double_t HLT_pT, Int_t runNo, Int_t lumiSec)
{
    Double_t prescale = 0.0;
    if (HLT_pT > 22 && HLT_pT < 30)
    {
        Int_t prescale_int = pp.hltPrescale("HLT_Photon22_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG18", runNo, lumiSec);
        prescale = (Double_t)prescale_int;
    }
    else if (HLT_pT > 30 && HLT_pT < 36)
    {
        Int_t prescale_int1 = pp.hltPrescale("HLT_Photon22_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG18", runNo, lumiSec);
        Int_t prescale_int2 = pp.hltPrescale("HLT_Photon30_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        Double_t one_over_prescale = 1 / ((Double_t)prescale_int1) + 1 / ((Double_t)prescale_int2);
        prescale = 1 / one_over_prescale;
    }
    else if (HLT_pT > 36 && HLT_pT < 50)
    {
        Int_t prescale_int1 = pp.hltPrescale("HLT_Photon22_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG18", runNo, lumiSec);
        Int_t prescale_int2 = pp.hltPrescale("HLT_Photon30_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        Int_t prescale_int3 = pp.hltPrescale("HLT_Photon36_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        Double_t one_over_prescale = 1 / ((Double_t)prescale_int1) + 1 / ((Double_t)prescale_int2) + 1 / ((Double_t)prescale_int3);
        prescale = 1 / one_over_prescale;
    }
    else if (HLT_pT > 50 && HLT_pT < 75)
    {
        Int_t prescale_int1 = pp.hltPrescale("HLT_Photon22_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG18", runNo, lumiSec);
        Int_t prescale_int2 = pp.hltPrescale("HLT_Photon30_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        Int_t prescale_int3 = pp.hltPrescale("HLT_Photon36_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        Int_t prescale_int4 = pp.hltPrescale("HLT_Photon50_v", runNo, lumiSec);

        Double_t one_over_prescale = 1 / ((Double_t)prescale_int1) + 1 / ((Double_t)prescale_int2) + 1 / ((Double_t)prescale_int3) +
                                     1 / ((Double_t)prescale_int4);
        prescale = 1 / one_over_prescale;
    }
    else if (HLT_pT > 75 && HLT_pT < 90)
    {
        Int_t prescale_int1 = pp.hltPrescale("HLT_Photon22_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG18", runNo, lumiSec);
        Int_t prescale_int2 = pp.hltPrescale("HLT_Photon30_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        Int_t prescale_int3 = pp.hltPrescale("HLT_Photon36_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        Int_t prescale_int4 = pp.hltPrescale("HLT_Photon50_v", runNo, lumiSec);
        Int_t prescale_int5 = pp.hltPrescale("HLT_Photon75_v", runNo, lumiSec);
        Double_t one_over_prescale = 1 / ((Double_t)prescale_int1) + 1 / ((Double_t)prescale_int2) + 1 / ((Double_t)prescale_int3) +
                                     1 / ((Double_t)prescale_int4) + 1 / ((Double_t)prescale_int5);
        prescale = 1 / one_over_prescale;
    }
    else if (HLT_pT > 90 && HLT_pT < 120)
    {
        Int_t prescale_int1 = pp.hltPrescale("HLT_Photon22_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG18", runNo, lumiSec);
        Int_t prescale_int2 = pp.hltPrescale("HLT_Photon30_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        Int_t prescale_int3 = pp.hltPrescale("HLT_Photon36_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        Int_t prescale_int4 = pp.hltPrescale("HLT_Photon50_v", runNo, lumiSec);
        Int_t prescale_int5 = pp.hltPrescale("HLT_Photon75_v", runNo, lumiSec);
        Int_t prescale_int6 = pp.hltPrescale("HLT_Photon90_v", runNo, lumiSec);
        Double_t one_over_prescale = 1 / ((Double_t)prescale_int1) + 1 / ((Double_t)prescale_int2) + 1 / ((Double_t)prescale_int3) +
                                     1 / ((Double_t)prescale_int4) + 1 / ((Double_t)prescale_int5) + 1 / ((Double_t)prescale_int6);
        prescale = 1 / one_over_prescale;
    }
    else if (HLT_pT > 120 && HLT_pT < 175)
    {
        Int_t prescale_int1 = pp.hltPrescale("HLT_Photon22_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG18", runNo, lumiSec);
        Int_t prescale_int2 = pp.hltPrescale("HLT_Photon30_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        Int_t prescale_int3 = pp.hltPrescale("HLT_Photon36_v", runNo, lumiSec) * pp.l1Prescale("L1_SingleEG26", runNo, lumiSec);
        Int_t prescale_int4 = pp.hltPrescale("HLT_Photon50_v", runNo, lumiSec);
        Int_t prescale_int5 = pp.hltPrescale("HLT_Photon75_v", runNo, lumiSec);
        Int_t prescale_int6 = pp.hltPrescale("HLT_Photon90_v", runNo, lumiSec);
        Int_t prescale_int7 = pp.hltPrescale("HLT_Photon120_v", runNo, lumiSec);
        Double_t one_over_prescale = 1 / ((Double_t)prescale_int1) + 1 / ((Double_t)prescale_int2) + 1 / ((Double_t)prescale_int3) +
                                     1 / ((Double_t)prescale_int4) + 1 / ((Double_t)prescale_int5) + 1 / ((Double_t)prescale_int6) +
                                     1 / ((Double_t)prescale_int7);
        prescale = 1 / one_over_prescale;
    }
    else if (HLT_pT > 175)
        prescale =1;

    return prescale;
}


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
