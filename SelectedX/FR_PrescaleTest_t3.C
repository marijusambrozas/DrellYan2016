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
#include "./header/NtupleHandle.h"

void FR_PrescaleTest_t3 (Bool_t DEBUG = kFALSE)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;
    PrescaleProvider pp("etc/prescale/triggerData2016");

    TFile *f;
    TString Dir = "../";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    for (Process_t pr=_SinglePhoton_B; pr<=_SinglePhoton_C; pr=next(pr))
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
        TH1D* h_HLT_pT_uncorr = new TH1D("h_HLT_pT_uncorr", "h_HLT_pT_uncorr", 500, 0, 500); h_HLT_pT_uncorr->Sumw2();
        TH1D* h_HLT_pT = new TH1D("h_HLT_pT", "h_HLT_pT", 500, 0, 500); h_HLT_pT->Sumw2();

        TChain *chain = new TChain(Mgr.TreeName[0]);

        chain->Add(Mgr.FullLocation[0]+"*.root");

        NtupleHandle *ntuple = new NtupleHandle(chain);

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

                for (Int_t i_tr=0; i_tr<ntuple->HLT_ntrig; i_tr++)
                {
                    if (!(ntuple->HLT_trigFired->at(i_tr))) continue;
                    if (!(ntuple->HLT_trigName->at(i_tr)).Contains("Photon")) continue;
                    if (ntuple->HLT_trigPt->at(i_tr) <= 22) continue;
                    if (ntuple->HLT_trigName->at(i_tr) == "HLT_Photon22_v*" && ntuple->HLT_trigPt->at(i_tr) > 22 && ntuple->HLT_trigPt->at(i_tr) <= 30)
                    {
                        i_22 = i_tr;
                        matched22 = 1;
                    }
                    if (ntuple->HLT_trigName->at(i_tr) == "HLT_Photon30_v*" && ntuple->HLT_trigPt->at(i_tr) > 30 && ntuple->HLT_trigPt->at(i_tr) <= 36)
                    {
                        i_30 = i_tr;
                        matched30 = 1;
                    }
                    if (ntuple->HLT_trigName->at(i_tr) == "HLT_Photon36_v*" && ntuple->HLT_trigPt->at(i_tr) > 36 && ntuple->HLT_trigPt->at(i_tr) <= 50)
                    {
                        i_36 = i_tr;
                        matched36 = 1;
                    }
                    if (ntuple->HLT_trigName->at(i_tr) == "HLT_Photon50_v*" && ntuple->HLT_trigPt->at(i_tr) > 50 && ntuple->HLT_trigPt->at(i_tr) <= 75)
                    {
                        i_50 = i_tr;
                        matched50 = 1;
                    }
                    if (ntuple->HLT_trigName->at(i_tr) == "HLT_Photon75_v*" && ntuple->HLT_trigPt->at(i_tr) > 75 && ntuple->HLT_trigPt->at(i_tr) <= 90)
                    {
                        i_75 = i_tr;
                        matched75 = 1;
                    }
                    if (ntuple->HLT_trigName->at(i_tr) == "HLT_Photon90_v*" && ntuple->HLT_trigPt->at(i_tr) > 90 && ntuple->HLT_trigPt->at(i_tr) <= 120)
                    {
                        i_90 = i_tr;
                        matched90 = 1;
                    }
                    if (ntuple->HLT_trigName->at(i_tr) == "HLT_Photon120_v*" && ntuple->HLT_trigPt->at(i_tr) > 120 && ntuple->HLT_trigPt->at(i_tr) <= 175)
                    {
                        i_120 = i_tr;
                        matched120 = 1;
                    }
                    if (ntuple->HLT_trigName->at(i_tr) == "HLT_Photon175_v*" && ntuple->HLT_trigPt->at(i_tr) > 175)
                    {
                        i_175 = i_tr;
                        matched175 = 1;
                    }

                }
                Int_t prescale_alt = 1;
                if (matched22+matched30+matched36+matched50+matched75+matched90+matched120+matched175 != 1) continue;
                if (matched22 == 1)
                {
                    prescale_alt = pp.hltPrescale("HLT_Photon22_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG18", ntuple->runNum, ntuple->lumiBlock);
                    h_HLT_pT->Fill(ntuple->HLT_trigPt->at(i_22), prescale_alt);
                    h_HLT_pT_uncorr->Fill(ntuple->HLT_trigPt->at(i_22));
                }
                else if (matched30 == 1)
                {
                    prescale_alt = pp.hltPrescale("HLT_Photon30_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG26", ntuple->runNum, ntuple->lumiBlock);
                    h_HLT_pT->Fill(ntuple->HLT_trigPt->at(i_30), prescale_alt);
                    h_HLT_pT_uncorr->Fill(ntuple->HLT_trigPt->at(i_30));
                }
                else if (matched36 == 1)
                {
                    prescale_alt = pp.hltPrescale("HLT_Photon36_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG26", ntuple->runNum, ntuple->lumiBlock);
                    h_HLT_pT->Fill(ntuple->HLT_trigPt->at(i_36), prescale_alt);
                    h_HLT_pT_uncorr->Fill(ntuple->HLT_trigPt->at(i_36));
                }
                else if (matched50 == 1)
                {
                    prescale_alt = 0;
                    if (pp.l1Prescale("L1_SingleEG36", ntuple->runNum, ntuple->lumiBlock) != 0)
                        prescale_alt = pp.hltPrescale("HLT_Photon50_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG36", ntuple->runNum, ntuple->lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG40", ntuple->runNum, ntuple->lumiBlock) != 0)
                        prescale_alt = pp.hltPrescale("HLT_Photon50_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG40", ntuple->runNum, ntuple->lumiBlock);
                    h_HLT_pT->Fill(ntuple->HLT_trigPt->at(i_50), prescale_alt);
                    h_HLT_pT_uncorr->Fill(ntuple->HLT_trigPt->at(i_50));
                }
                else if (matched75 == 1)
                {
                    prescale_alt = 0;
                    if (pp.l1Prescale("L1_SingleEG36", ntuple->runNum, ntuple->lumiBlock) != 0)
                        prescale_alt = pp.hltPrescale("HLT_Photon75_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG36", ntuple->runNum, ntuple->lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG40", ntuple->runNum, ntuple->lumiBlock) != 0)
                        prescale_alt = pp.hltPrescale("HLT_Photon75_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG40", ntuple->runNum, ntuple->lumiBlock);
                    h_HLT_pT->Fill(ntuple->HLT_trigPt->at(i_75), prescale_alt);
                    h_HLT_pT_uncorr->Fill(ntuple->HLT_trigPt->at(i_75));
                }
                else if (matched90 == 1)
                {
                    prescale_alt = 0;
                    if (pp.l1Prescale("L1_SingleEG36", ntuple->runNum, ntuple->lumiBlock) != 0)
                        prescale_alt = pp.hltPrescale("HLT_Photon90_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG36", ntuple->runNum, ntuple->lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG40", ntuple->runNum, ntuple->lumiBlock) != 0)
                        prescale_alt = pp.hltPrescale("HLT_Photon90_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG40", ntuple->runNum, ntuple->lumiBlock);
                    h_HLT_pT->Fill(ntuple->HLT_trigPt->at(i_90), prescale_alt);
                    h_HLT_pT_uncorr->Fill(ntuple->HLT_trigPt->at(i_90));
                }
                else if (matched120 == 1)
                {
                    prescale_alt = 0;
                    if (pp.l1Prescale("L1_SingleEG36", ntuple->runNum, ntuple->lumiBlock) != 0)
                        prescale_alt = pp.hltPrescale("HLT_Photon120_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG36", ntuple->runNum, ntuple->lumiBlock);
                    if (pp.l1Prescale("L1_SingleEG40", ntuple->runNum, ntuple->lumiBlock) != 0)
                        prescale_alt = pp.hltPrescale("HLT_Photon120_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG40", ntuple->runNum, ntuple->lumiBlock);
                    h_HLT_pT->Fill(ntuple->HLT_trigPt->at(i_120), prescale_alt);
                    h_HLT_pT_uncorr->Fill(ntuple->HLT_trigPt->at(i_120));
                }
                else if (matched175 == 1)
                {
                     prescale_alt = 0;
                     if (pp.l1Prescale("L1_SingleEG30", ntuple->runNum, ntuple->lumiBlock) != 0)
                         prescale_alt = pp.hltPrescale("HLT_Photon175_v", ntuple->runNum, lntuple->umiBlock) * pp.l1Prescale("L1_SingleEG30", ntuple->runNum, ntuple->lumiBlock);
                     if (pp.l1Prescale("L1_SingleEG32", ntuple->runNum, ntuple->lumiBlock) != 0)
                         prescale_alt = pp.hltPrescale("HLT_Photon175_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG32", ntuple->runNum, ntuple->lumiBlock);
                     if (pp.l1Prescale("L1_SingleEG34", ntuple->runNum, ntuple->lumiBlock) != 0)
                         prescale_alt = pp.hltPrescale("HLT_Photon175_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG34", ntuple->runNum, ntuple->lumiBlock);
                     if (pp.l1Prescale("L1_SingleEG36", ntuple->runNum, ntuple->lumiBlock) != 0)
                         prescale_alt = pp.hltPrescale("HLT_Photon175_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG36", ntuple->runNum, ntuple->lumiBlock);
                     if (pp.l1Prescale("L1_SingleEG38", ntuple->runNum, ntuple->lumiBlock) != 0)
                         prescale_alt = pp.hltPrescale("HLT_Photon175_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG38", ntuple->runNum, ntuple->lumiBlock);
                     if (pp.l1Prescale("L1_SingleEG40", ntuple->runNum, ntuple->lumiBlock) != 0)
                         prescale_alt = pp.hltPrescale("HLT_Photon175_v", ntuple->runNum, ntuple->lumiBlock) * pp.l1Prescale("L1_SingleEG40", ntuple->runNum, ntuple->lumiBlock);                    h_HLT_pT->Fill(trig_pT->at(i_175), prescale_alt);
                    h_HLT_pT->Fill(ntuple->HLT_trigPt->at(i_175), prescale_alt);
                    h_HLT_pT_uncorr->Fill(ntuple->HLT_trigPt->at(i_175));
                }

            if (DEBUG == kFALSE) bar.Draw(i);

        }// End of event iteration

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

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of FR_PrescaleTest_t3()
