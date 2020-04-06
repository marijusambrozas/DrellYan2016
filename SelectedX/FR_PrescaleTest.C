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

void FR_PrescaleTest (Bool_t DEBUG = kFALSE)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;
    PrescaleProvider pp("etc/prescale/triggerData2016");

    TFile *f;
//    TString Dir = "../";
    TString Dir = "/media/sf_DATA/FR/Electron/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

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
        Int_t runNum;
        Int_t lumiBlock;
        Int_t prescale_alt;

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

            for (UInt_t i_ele=0; i_ele<p_T->size(); i_ele++)
            {
                if (p_T->at(i_ele) != p_T->at(i_ele))
                {
                    cout << p_T->at(i_ele) << endl;
                    continue;
                }
                if (p_T->at(i_ele) <= 28) continue;
                if (DEBUG == kTRUE) cout << "i_ele = " << i_ele << endl;

                Int_t matched22=0, matched30=0, matched36=0, matched50=0, matched75=0, matched90=0, matched120=0, matched175=0;
                Int_t i_22=-1, i_30=-1, i_36=-1, i_50=-1, i_75=-1, i_90=-1, i_120=-1, i_175=-1;
                prescale_alt = 1;

                for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
                {
                    if (trig_fired->at(i_tr) < 22) continue;
                    if (((UInt_t)(trig_matched->at(i_tr))) == i_ele)
                    {
                        if (trig_fired->at(i_tr) == 22 && trig_pT->at(i_tr) > 22 && trig_pT->at(i_tr) < 30)
                        {
                            i_22 = i_tr;
                            matched22 = 1;
                        }
                        if (trig_fired->at(i_tr) == 30 && trig_pT->at(i_tr) > 30 && trig_pT->at(i_tr) < 36)
                        {
                            i_30 = i_tr;
                            matched30 = 1;
                        }
                        if (trig_fired->at(i_tr) == 36 && trig_pT->at(i_tr) > 36 && trig_pT->at(i_tr) < 50)
                        {
                            i_36 = i_tr;
                            matched36 = 1;
                        }
                        if (trig_fired->at(i_tr) == 50 && trig_pT->at(i_tr) > 50 && trig_pT->at(i_tr) < 75)
                        {
                            i_50 = i_tr;
                            matched50 = 1;
                        }
                        if (trig_fired->at(i_tr) == 75 && trig_pT->at(i_tr) > 75 && trig_pT->at(i_tr) < 90)
                        {
                            i_75 = i_tr;
                            matched75 = 1;
                        }
                        if (trig_fired->at(i_tr) == 90 && trig_pT->at(i_tr) > 90 && trig_pT->at(i_tr) < 120)
                        {
                            i_90 = i_tr;
                            matched90 = 1;
                        }
                        if (trig_fired->at(i_tr) == 120 && trig_pT->at(i_tr) > 120 && trig_pT->at(i_tr) < 175)
                        {
                            i_120 = i_tr;
                            matched120 = 1;
                        }
                        else if (trig_fired->at(i_tr) == 175 && trig_pT->at(i_tr) > 175)
                        {
                            i_175 = i_tr;
                            matched175 = 1;
                        }

                    }
                }
                if (matched22+matched30+matched36+matched50+matched75+matched90+matched120 > 1) continue;
                if (matched22==0 && matched30==0 && matched36==0 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0) continue;
                else if (matched22 == 1 && matched30 == 0 && matched36 == 0 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = pp.hltPrescale("HLT_Photon22_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG18", runNum, lumiBlock);
                    h_HLT_pT->Fill(trig_pT->at(i_22), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_22));
                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele));
                }
                else if (matched30 == 1 && matched36 == 0 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = pp.hltPrescale("HLT_Photon30_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG26", runNum, lumiBlock);
                    h_HLT_pT->Fill(trig_pT->at(i_30), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_30));
                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele));
                }
                else if (matched36 == 1 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = pp.hltPrescale("HLT_Photon36_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG26", runNum, lumiBlock);
                    h_HLT_pT->Fill(trig_pT->at(i_36), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_36));
                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele));
                }
                else if (matched50 == 1 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = pp.hltPrescale("HLT_Photon50_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock);
                    h_HLT_pT->Fill(trig_pT->at(i_50), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_50));
                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele));
                }
                else if (matched75 == 1 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = pp.hltPrescale("HLT_Photon75_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock);
                    h_HLT_pT->Fill(trig_pT->at(i_75), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_75));
                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele));
                }
                else if (matched90 == 1 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = pp.hltPrescale("HLT_Photon90_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock);
                    h_HLT_pT->Fill(trig_pT->at(i_90), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_90));
                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele));
                }
                else if (matched120 == 1 && matched175 == 0)
                {
                    prescale_alt = pp.hltPrescale("HLT_Photon120_v", runNum, lumiBlock) * pp.l1Prescale("L1_SingleEG40", runNum, lumiBlock);
                    h_HLT_pT->Fill(trig_pT->at(i_120), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_120));
                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele));
                }
                else if (matched175 == 1)
                {
                    prescale_alt = 1.;
                    h_HLT_pT->Fill(trig_pT->at(i_175), prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_175));
                    h_pT->Fill(p_T->at(i_ele), prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele));
                }
                else continue;

            }// End of i_ele iteration

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
        if (!f->IsOpen()) cout << "File " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " has been closed successfully.\n" << endl;
        else cout << "FILE " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of FR_PrescaleTest()
