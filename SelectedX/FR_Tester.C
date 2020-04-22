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
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/FileMgr.h"

void E_QCD_Tester (TString FRtype, Bool_t DEBUG);
void E_WJET_Tester (TString FRtype, Bool_t DEBUG);
void Mu_QCD_Tester (TString FRtype, Bool_t DEBUG);
void Mu_WJET_Tester (TString FRtype, Bool_t DEBUG);
void EMu_WJET_Tester (TString FRtype, Bool_t DEBUG);

// -- Drell-Yan mass bins -- //
const Int_t binnum = 43;
const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126,
                               133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};


void FR_Tester (TString WhichX = "", TString FRtype = "MC")
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
        cout << "\n*****  EMu_WJET_Tester(" << FRtype << ")  *****" << endl;
        EMu_WJET_Tester(FRtype, DEBUG);
    }
    else if (whichX.Contains("MU"))
    {
        if (whichX.Contains("QCD"))
        {
            Xselected++;
            cout << "\n*****  Mu_QCD_Tester(" << FRtype << ")  *****" << endl;
            Mu_QCD_Tester(FRtype, DEBUG);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            Xselected++;
            cout << "\n*****  Mu_WJET_Tester(" << FRtype << ")  *****" << endl;
            Mu_WJET_Tester(FRtype, DEBUG);
        }
    }
    else if (whichX.Contains("E"))
    {
        if (whichX.Contains("QCD"))
        {
            Xselected++;
            cout << "\n*****  E_QCD_Tester(" << FRtype << ")  *****" << endl;
            E_QCD_Tester(FRtype, DEBUG);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            Xselected++;
            cout << "\n*****  E_WJET_Tester(" << FRtype << ")  *****" << endl;
            E_WJET_Tester(FRtype, DEBUG);
        }
    }

    if (Xselected == 0) cout << "Wrong arument! \nType in: >> .x FR_HistMaker.C+(\"whichX\")" << endl;

} // End of FR_Tester()


/// -------------------------------- Electron Channel ------------------------------------ ///
void E_QCD_Tester (TString FRtype, Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Electron/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    DYAnalyzer *analyzer = new DYAnalyzer("Ele23Ele12");

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues_ele(Dir+"FakeRate_electron.root", FRtype);

    TH1D *h_mass_2FR = new TH1D("h_mass_2FR", "", 10, 0, 1000);
    TH1D *h_mass_1FR = new TH1D("h_mass_1FR", "", 10, 0, 1000);

    for (Process_t pr=_QCDEMEnriched_20to30; pr<_EndOf_QCDEMEnriched_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<double> *etaSC = new std::vector<double>;
        std::vector<double> *Full5x5_SigmaIEtaIEta = new std::vector<double>;
        std::vector<double> *HoverE = new std::vector<double>;
        std::vector<double> *dEtaInSeed = new std::vector<double>;
        std::vector<double> *dPhiIn = new std::vector<double>;
        std::vector<int> *mHits = new std::vector<int>;
        std::vector<int> *passMediumID = new std::vector<int>;
        Int_t nPU;
        Double_t PVz;
        Double_t gen_weight;
        Double_t prefiring_weight;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("etaSC", 1);
        chain->SetBranchStatus("Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("HoverE", 1);
        chain->SetBranchStatus("dEtaInSeed", 1);
        chain->SetBranchStatus("dPhiIn", 1);
        chain->SetBranchStatus("mHits", 1);
        chain->SetBranchStatus("passMediumID", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);

        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("etaSC", &etaSC);
        chain->SetBranchAddress("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("HoverE", &HoverE);
        chain->SetBranchAddress("dEtaInSeed", &dEtaInSeed);
        chain->SetBranchAddress("dPhiIn", &dPhiIn);
        chain->SetBranchAddress("mHits", &mHits);
        chain->SetBranchAddress("passMediumID", &passMediumID);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;


        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // QCD selection
            if (p_T->size() != 2) continue;
            if ((fabs(etaSC->at(0)) > 1.4442 && fabs(etaSC->at(0)) < 1.566) || fabs((etaSC->at(1)) > 1.4442 && fabs(etaSC->at(1)) < 1.566)) continue;
            if (fabs(etaSC->at(0)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.013 || HoverE->at(0) >= 0.13 ||
                                                fabs(dEtaInSeed->at(0)) >= 0.01 || fabs(dPhiIn->at(0)) >= 0.07 || mHits->at(0) > 1)) continue;
            if (fabs(etaSC->at(1)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.013 || HoverE->at(1) >= 0.13 ||
                                                fabs(dEtaInSeed->at(1)) >= 0.01 || fabs(dPhiIn->at(1)) >= 0.07 || mHits->at(1) > 1)) continue;
            if (fabs(etaSC->at(0)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.035 || HoverE->at(0) >= 0.13 || mHits->at(0) > 1)) continue;
            if (fabs(etaSC->at(1)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.035 || HoverE->at(1) >= 0.13 || mHits->at(1) > 1)) continue;
            if (passMediumID->at(0) == 1 && passMediumID->at(1)  == 1) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << passMediumID->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << passMediumID->at(1) << endl;

            if (DEBUG == kTRUE)
            {
                cout << "Evt " << i << endl;
                cout << "nElectrons = " << p_T->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "\tpassMediumID[0] = " << passMediumID->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\tpassMediumID[1] = " << passMediumID->at(1) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector ele1, ele2;
            ele1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
            ele2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
            Double_t mass = (ele1+ele2).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- PVz weights -- //
            Double_t PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = prefiring_weight;

            if (DEBUG == kTRUE) cout << "PU weight: " << PUWeight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight << endl;

            // -- Normalization -- //
            Double_t TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            // -- FR WEIGHTS -- //
            if (passMediumID->at(0) == 1 && passMediumID->at(1) == 0) // First passing
            {
                Double_t FRweight, FR;
                FR = analyzer->FakeRate_ele(p_T->at(1), etaSC->at(1));
                FRweight = FR / (1 - FR);
                h_mass_1FR->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * FRweight);
            }
            else if (passMediumID->at(0) == 0 && passMediumID->at(1) == 1) // Second passing
            {
                Double_t FRweight, FR;
                FR = analyzer->FakeRate_ele(p_T->at(0), etaSC->at(0));
                FRweight = FR / (1 - FR);
                h_mass_1FR->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * FRweight);
            }
            else // Both failing
            {
                Double_t FRweight = 1;
                Double_t FR1, FR2;
                FR1 = analyzer->FakeRate_ele(p_T->at(0), etaSC->at(0));
                FR2 = analyzer->FakeRate_ele(p_T->at(1), etaSC->at(1));
                FRweight = 2 * (FR1 / (1 - FR1)) * (FR2 / (1 - FR2));
                if (DEBUG == kTRUE) cout << "FR1 = " << FR1 << "   FR2 = " << FR2 << "   FRweight = " << FRweight << endl;
                h_mass_2FR->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * FRweight);
            }

        }// End of event iteration

        cout << " Finished.\n" << endl;

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TCanvas *c_comparison = new TCanvas("FR test","FR test (QCD)", 800, 800);
    c_comparison->SetRightMargin(0.05);
    c_comparison->SetTopMargin(0.05);
    c_comparison->SetBottomMargin(0.12);
    c_comparison->SetLeftMargin(0.13);
    h_mass_1FR->SetDirectory(0);
    h_mass_2FR->SetDirectory(0);
//    for (Int_t i=1; i<=10; i++)
//    {
//        h_mass_2FR->SetBinError(i, sqrt(h_mass_2FR->GetBinContent(i)));
//    }
    h_mass_1FR->SetStats(0);
    h_mass_2FR->SetStats(0);
    h_mass_1FR->SetMarkerStyle(kFullDotLarge);
    h_mass_2FR->SetMarkerStyle(kFullSquare);
    h_mass_1FR->SetLineColor(kBlack);
    h_mass_2FR->SetLineColor(kRed);
    h_mass_2FR->SetMarkerColor(kRed);
    h_mass_1FR->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");
    h_mass_1FR->GetXaxis()->SetNoExponent(1);
    h_mass_1FR->GetXaxis()->SetMoreLogLabels(1);
    h_mass_1FR->GetXaxis()->SetTitleOffset(1);
    h_mass_1FR->GetXaxis()->SetTitleSize(0.05);
    h_mass_1FR->GetXaxis()->SetLabelSize(0.04);
    h_mass_1FR->GetYaxis()->SetTitle("Number of events");
    h_mass_1FR->GetYaxis()->SetTitleSize(0.05);
    h_mass_1FR->GetYaxis()->SetTitleOffset(1.25);
    h_mass_1FR->GetYaxis()->SetLabelSize(0.04);
    h_mass_1FR->GetYaxis()->SetRangeUser(10, 1e6);
    TLegend *legend = new TLegend(0.4, 0.8, 0.95, 0.95);
    legend->AddEntry(h_mass_1FR, "QCD MC (1 fail, FR applied)", "lp");
    legend->AddEntry(h_mass_2FR, "QCD MC (2 fail, FR applied x2)", "lp");
    h_mass_1FR->Draw();
    h_mass_2FR->Draw("same");
    legend->Draw();
    c_comparison->SetLogy();
    c_comparison->SetGridx();
    c_comparison->SetGridy();
    c_comparison->Update();

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of E_QCD_Tester()


void E_WJET_Tester (TString FRtype, Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Electron/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    DYAnalyzer *analyzer = new DYAnalyzer("Ele23Ele12");

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues_ele(Dir+"FakeRate_electron.root", FRtype, 1);

    TH1D *h_mass, *h_mass_compare;
    TH1D *h_mass_1FR = new TH1D("h_mass_1FR", "h_mass_1FR", binnum, massbins);

    for (Process_t pr=_WJets; pr<=_WJets_ext2v5; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<double> *etaSC = new std::vector<double>;
        std::vector<double> *Full5x5_SigmaIEtaIEta = new std::vector<double>;
        std::vector<double> *HoverE = new std::vector<double>;
        std::vector<double> *dEtaInSeed = new std::vector<double>;
        std::vector<double> *dPhiIn = new std::vector<double>;
        std::vector<int> *mHits = new std::vector<int>;
        std::vector<int> *passMediumID = new std::vector<int>;
        Int_t nPU;
        Double_t PVz;
        Double_t gen_weight;
        Double_t prefiring_weight;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("etaSC", 1);
        chain->SetBranchStatus("Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("HoverE", 1);
        chain->SetBranchStatus("dEtaInSeed", 1);
        chain->SetBranchStatus("dPhiIn", 1);
        chain->SetBranchStatus("mHits", 1);
        chain->SetBranchStatus("passMediumID", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);

        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("etaSC", &etaSC);
        chain->SetBranchAddress("Full5x5_SigmaIEtaIEta", &Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("HoverE", &HoverE);
        chain->SetBranchAddress("dEtaInSeed", &dEtaInSeed);
        chain->SetBranchAddress("dPhiIn", &dPhiIn);
        chain->SetBranchAddress("mHits", &mHits);
        chain->SetBranchAddress("passMediumID", &passMediumID);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;


        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // QCD selection
            if (p_T->size() != 2) continue;
            if (passMediumID->at(0) == 1 && passMediumID->at(1)  == 1) continue;
            if (passMediumID->at(0) == 0 && passMediumID->at(1)  == 0) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
            if ((fabs(etaSC->at(0)) > 1.4442 && fabs(etaSC->at(0)) < 1.566) || fabs((etaSC->at(1)) > 1.4442 && fabs(etaSC->at(1)) < 1.566)) continue;
            if (fabs(etaSC->at(0)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.013 || HoverE->at(0) >= 0.13 ||
                                                fabs(dEtaInSeed->at(0)) >= 0.01 || fabs(dPhiIn->at(0)) >= 0.07 || mHits->at(0) > 1)) continue;
            if (fabs(etaSC->at(1)) < 1.4442 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.013 || HoverE->at(1) >= 0.13 ||
                                                fabs(dEtaInSeed->at(1)) >= 0.01 || fabs(dPhiIn->at(1)) >= 0.07 || mHits->at(1) > 1)) continue;
            if (fabs(etaSC->at(0)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(0) >= 0.035 || HoverE->at(0) >= 0.13 || mHits->at(0) > 1)) continue;
            if (fabs(etaSC->at(1)) > 1.566 && (Full5x5_SigmaIEtaIEta->at(1) >= 0.035 || HoverE->at(1) >= 0.13 || mHits->at(1) > 1)) continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << passMediumID->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << passMediumID->at(1) << endl;

            if (DEBUG == kTRUE)
            {
                cout << "Evt " << i << endl;
                cout << "nElectrons = " << p_T->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "\tpassMediumID[0] = " << passMediumID->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\tpassMediumID[1] = " << passMediumID->at(1) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector ele1, ele2;
            ele1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Elec);
            ele2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Elec);
            Double_t mass = (ele1+ele2).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- PVz weights -- //
            Double_t PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = prefiring_weight;

            if (DEBUG == kTRUE) cout << "PU weight: " << PUWeight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight << endl;

            // -- Normalization -- //
            Double_t TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            // -- FR WEIGHTS -- //
            Double_t FRweight=1, FR=-1;
            if (passMediumID->at(0) == 1 && passMediumID->at(1) == 0) // First passing
                FR = analyzer->FakeRate_ele(p_T->at(1), etaSC->at(1));
            else if (passMediumID->at(0) == 0 && passMediumID->at(1) == 1) // First failing
                FR = analyzer->FakeRate_ele(p_T->at(0), etaSC->at(0));
            FRweight = FR / (1 - FR);
            h_mass_1FR->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * FRweight);

        }// End of event iteration

        cout << " Finished.\n" << endl;

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TFile *f_orig = new TFile("/media/sf_DATA/SelectedEE/Histos/Hist_SelectedEE_Bkg_Full.root", "READ");
    f_orig->GetObject("h_mass_SelectedEE_WJets_Full", h_mass);
    h_mass->SetDirectory(0);
    f_orig->Close();
    TFile *f_est = new TFile("/media/sf_DATA/SelectedEE/Histos/EstWJets_EE.root", "READ");
    f_est->GetObject("h_WJET_est_fit", h_mass_compare);
    h_mass_compare->SetDirectory(0);
    f_est->Close();

    TCanvas *c_comparison = new TCanvas("FR test","FR test (W+Jets)", 800, 800);
    c_comparison->SetRightMargin(0.05);
    c_comparison->SetTopMargin(0.05);
    c_comparison->SetBottomMargin(0.12);
    c_comparison->SetLeftMargin(0.13);
    h_mass_1FR->SetDirectory(0);
    h_mass->SetStats(0);
    h_mass_1FR->SetStats(0);
    h_mass_compare->SetStats(0);
    h_mass->SetMarkerStyle(kFullDotLarge);
    h_mass_1FR->SetMarkerStyle(kFullSquare);
    h_mass_compare->SetMarkerStyle(33);
    h_mass->SetLineColor(kBlack);
    h_mass_1FR->SetLineColor(kRed);
    h_mass_1FR->SetMarkerColor(kRed);
    h_mass_compare->SetLineColor(kGreen+2);
    h_mass_compare->SetMarkerColor(kGreen+2);
    h_mass_compare->SetMarkerSize(1.5);
    h_mass->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");
    h_mass->GetXaxis()->SetNoExponent(1);
    h_mass->GetXaxis()->SetMoreLogLabels(1);
    h_mass->GetXaxis()->SetTitleOffset(1);
    h_mass->GetXaxis()->SetTitleSize(0.05);
    h_mass->GetXaxis()->SetLabelSize(0.04);
    h_mass->GetYaxis()->SetTitle("Number of events");
    h_mass->GetYaxis()->SetTitleSize(0.05);
    h_mass->GetYaxis()->SetTitleOffset(1.25);
    h_mass->GetYaxis()->SetLabelSize(0.04);
    h_mass->GetYaxis()->SetRangeUser(1, 1e4);
    TLegend *legend = new TLegend(0.5, 0.7, 0.95, 0.95);
    legend->AddEntry(h_mass, "W+Jets MC (2 pass)", "lp");
    legend->AddEntry(h_mass_1FR, "W+Jets MC (1 fail, FR applied)", "lp");
    legend->AddEntry(h_mass_compare, "W+Jets est (data-driven template)", "lp");
    h_mass->Draw();
    h_mass_1FR->Draw("same");
    h_mass_compare->Draw("same");
    legend->Draw();
    c_comparison->SetLogx();
    c_comparison->SetLogy();
    c_comparison->SetGridx();
    c_comparison->SetGridy();
    c_comparison->Update();

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
} // End of E_WJET_Tester()


/// -------------------------------- Muon Channel ------------------------------------ ///
void Mu_QCD_Tester (TString FRtype, Bool_t DEBUG)
{
    return; // Not ready yet
    /*
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    f = new TFile(Dir+"QCDest_Mu"+debug+".root", "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Mu50");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_BtoF();
    analyzer->SetupEfficiencyScaleFactor_GtoH();

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues(Dir+"FakeRate_muon.root", "sigCtrl_template");

    TH1D *h_mass_test[_EndOf_SingleMuon_Normal];
    TH1D *h_mass_test_SS[_EndOf_SingleMuon_Normal];

    TH1D* h_FRweight = new TH1D("h_FRweight", "FR weights", 100, 0, 0.4);

    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();
        Int_t nPass = 0;
        Double_t avgFRweight = 0;

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_mass = new TH1D("h_mass_"+Mgr.Procname[pr], "h_mass_"+Mgr.Procname[pr], binnum, massbins); h_mass->Sumw2();
        TH1D* h_mass_SS = new TH1D("h_mass_SS_"+Mgr.Procname[pr], "h_mass_SS_"+Mgr.Procname[pr], binnum, massbins); h_mass_SS->Sumw2();
        TH1D* h_mass_forFit = new TH1D("h_mass_forFit_"+Mgr.Procname[pr], "h_mass_forFit_"+Mgr.Procname[pr], 37, 15, 200); h_mass_forFit->Sumw2();
        TH1D* h_mass_SS_forFit = new TH1D("h_mass_SS_forFit_"+Mgr.Procname[pr], "h_mass_SS_forFit_"+Mgr.Procname[pr], 37, 15, 200); h_mass_SS_forFit->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX_"+Mgr.Procname[pr], "h_nVTX_"+Mgr.Procname[pr], 50, 0, 50); h_nVTX->Sumw2();
        TH1D* h_pT_lead = new TH1D("h_pT_lead_"+Mgr.Procname[pr], "h_pT_lead_"+Mgr.Procname[pr], 100, 0, 1000); h_pT_lead->Sumw2();
        TH1D* h_pT_sublead = new TH1D("h_pT_sublead_"+Mgr.Procname[pr], "h_pT_sublead_"+Mgr.Procname[pr], 100, 0, 1000); h_pT_sublead->Sumw2();
        TH2D* h2_pT = new TH2D("h2_pT_"+Mgr.Procname[pr], "h2_pT_"+Mgr.Procname[pr], 49, 10, 500, 49, 10, 500);
        h_mass_test[pr] = new TH1D("h_test_"+Mgr.Procname[pr], "", binnum, massbins); h_mass_test[pr]->Sumw2();
        h_mass_test[pr]->SetDirectory(0);
        h_mass_test_SS[pr] = new TH1D("h_test_SS_"+Mgr.Procname[pr], "", binnum, massbins); h_mass_test_SS[pr]->Sumw2();
        h_mass_test_SS[pr]->SetDirectory(0);

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t evt_weight;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");

        if (type == 0) // Mu50 files
        {
            chain->Add(Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
            if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        }
        else if (type == 1) // Triggerless files
        {
            chain->Add(Dir+"SelectedForBKGest_Mu_Triggerless"+Mgr.Procname[Mgr.CurrentProc]+".root");
            if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_Mu_Triggerless"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        }
        else
        {
            cout << "Wrong type! Select 0 (Mu50) or 1 (Triggerless)" << endl;
            return;
        }
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("relPFiso", 1);
        chain->SetBranchStatus("TRKiso", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
//        chain->SetBranchStatus("evt_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);
        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("charge", &charge);
        chain->SetBranchAddress("relPFiso", &relPFiso);
        chain->SetBranchAddress("TRKiso", &TRKiso);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
//        chain->SetBranchAddress("evt_weight", &evt_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        UInt_t nPassMu=0, nFailMu=0;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            for (UInt_t m=0; m<relPFiso->size(); m++)
            {
                if (relPFiso->at(m) < 0.15) nPassMu++;
                else nFailMu++;
            }
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // QCD selection
            if (p_T->size() != 2) continue;
//            if (charge->at(0) == charge->at(1)) continue;
            if (relPFiso->at(0) < 0.15 || relPFiso->at(1) < 0.15) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
//            if (p_T->at(0) < 2 || p_T->at(1) < 2) continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso->at(1) << endl;

            nPass++;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 5) break;
                cout << "Evt " << i << endl;
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector mu1, mu2;
            mu1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Mu);
            mu2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Mu);
            Double_t mass = (mu1+mu2).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t weight1 = 0, weight2 = 0, effweight = 1;
            if (Mgr.isMC == kTRUE)
            {
                weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF(mu1, mu2);
                weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH(mu1, mu2);
                effweight = (Lumi_BtoF * weight1 + Lumi_GtoH * weight2) / Lumi;
            }

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "Eff weight: " << effweight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight <<
                                        "\tTop pT weight: " << TopPtWeight << endl;

            // -- FR WEIGHTS -- //
            Double_t FRweight = 1;
            Double_t FR1, FR2;
            FR1 = analyzer->FakeRate(p_T->at(0), eta->at(0));
            FR2 = analyzer->FakeRate(p_T->at(1), eta->at(1));
            FRweight = FR1 / (1 - FR1) * FR2 / (1 - FR2);
            if (DEBUG == kTRUE) cout << "FR1 = " << FR1 << "   FR2 = " << FR2 << "   FRweight = " << FRweight << endl;
            avgFRweight += FRweight;
            h_FRweight->Fill(FRweight);

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            if (charge->at(0) != charge->at(1))
            {
                h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_forFit->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_test[pr]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            else
            {
                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_SS_forFit->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_test_SS[pr]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            if (mu1.Pt() > mu2.Pt())
            {
                h2_pT->Fill(mu1.Pt(), mu2.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_pT_lead->Fill(mu1.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_pT_sublead->Fill(mu2.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            }
            else
            {
                h2_pT->Fill(mu2.Pt(), mu1.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_pT_lead->Fill(mu2.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_pT_sublead->Fill(mu1.Pt(), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            }

        }// End of event iteration

        cout << "\t " << nPass << " events have passed the selection." << endl;
        cout << "\t Average FR weight: " << avgFRweight / NEvents << endl;
        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }
        cout << "\t # passed muons: " << nPassMu << endl;
        cout << "\t # failed muons: " << nFailMu << endl;

        f->cd();
        cout << "\tWriting into file...";

        h_mass->Write();
        h_mass_forFit->Write();
        h_mass_SS->Write();
        h_mass_SS_forFit->Write();
        h_nVTX->Write();
        h_pT_lead->Write();
        h_pT_sublead->Write();
        h2_pT->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TCanvas *c_data = new TCanvas("data","data", 800, 800);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_C]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_D]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_E]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_F]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_G]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_H]);
    h_mass_test[_SingleMuon_B]->Draw("hist");
    c_data->SetLogx();
    c_data->Update();

    TCanvas *c_data_SS = new TCanvas("data_SS","data SS", 800, 800);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_C]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_D]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_E]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_F]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_G]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_H]);
    h_mass_test_SS[_SingleMuon_B]->Draw("hist");
    c_data_SS->SetLogx();
    c_data_SS->Update();

    TCanvas *c_DY = new TCanvas("dy","dy", 800, 800);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_50to100]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_100to200]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_200to400]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_400to500]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_500to700]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_700to800]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_800to1000]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_1000to1500]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_1500to2000]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_2000to3000]);
    h_mass_test[_DY_10to50]->Draw("hist");
    c_DY->SetLogx();
    c_DY->Update();

    TCanvas *c_DY_SS = new TCanvas("dy_SS","dy SS", 800, 800);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_50to100]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_100to200]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_200to400]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_400to500]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_500to700]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_700to800]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_800to1000]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_1000to1500]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_1500to2000]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_2000to3000]);
    h_mass_test_SS[_DY_10to50]->Draw("hist");
    c_DY_SS->SetLogx();
    c_DY_SS->Update();

    TCanvas *c_tt = new TCanvas("ttbar","ttbar", 800, 800);
    h_mass_test[_ttbar]->Add(h_mass_test[_ttbar_700to1000]);
    h_mass_test[_ttbar]->Add(h_mass_test[_ttbar_1000toInf]);
    h_mass_test[_ttbar]->Draw("hist");
    c_tt->SetLogx();
    c_tt->Update();

    TCanvas *c_tt_SS = new TCanvas("ttbar_SS","ttbar SS", 800, 800);
    h_mass_test_SS[_ttbar]->Add(h_mass_test_SS[_ttbar_700to1000]);
    h_mass_test_SS[_ttbar]->Add(h_mass_test_SS[_ttbar_1000toInf]);
    h_mass_test_SS[_ttbar]->Draw("hist");
    c_tt_SS->SetLogx();
    c_tt_SS->Update();

    TCanvas *c_wjets = new TCanvas("wjets","wjets", 800, 800);
    h_mass_test[_WJets]->Add(h_mass_test[_WJets_ext2v5]);
    h_mass_test[_WJets]->Draw("hist");
    c_wjets->SetLogx();
    c_wjets->Update();

    TCanvas *c_wjets_SS = new TCanvas("wjets_SS","wjets SS", 800, 800);
    h_mass_test_SS[_WJets]->Add(h_mass_test_SS[_WJets_ext2v5]);
    h_mass_test_SS[_WJets]->Draw("hist");
    c_wjets_SS->SetLogx();
    c_wjets_SS->Update();

    TH1D *h_wjets_SS_scaled = (TH1D*)h_mass_test_SS[_WJets]->Clone("h_wjets_SS_scaled");
    h_wjets_SS_scaled->Rebin(8);
    h_wjets_SS_scaled->Scale(3);
    h_wjets_SS_scaled->SetLineColor(kRed);
    h_wjets_SS_scaled->SetDirectory(0);
    TH1D *h_wjets_OS_compare = (TH1D*)h_mass_test[_WJets]->Clone("h_wjets_OS_compare");
    h_wjets_OS_compare->SetDirectory(0);
    h_wjets_OS_compare->Rebin(8);
    TCanvas *c_wjets_compare = new TCanvas("wjets_compare", "WJets OS vs SSx3", 800, 800);
    h_wjets_OS_compare->Draw("hist");
    h_wjets_SS_scaled->Draw("histsame");
    c_wjets_compare->SetLogx();
    c_wjets_compare->Update();
    cout << "Ratio of W+Jets opposite-sign and same-sign integrals: " << h_mass_test[_WJets]->Integral()/h_mass_test_SS[_WJets]->Integral() << endl;

    TCanvas *c_FRweight = new TCanvas("FRweight","FR weights", 800, 800);
    h_FRweight->Draw();
    c_FRweight->Update();

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"QCDest_Mu"+debug+".root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"QCDest_Mu"+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
    */
} // End of Mu_QCD_Tester()


void Mu_WJET_Tester (TString FRtype, Bool_t DEBUG)
{
    return; // Not ready yet
    /*
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    // -- Output ROOTFile -- //
    f = new TFile(Dir+"WJETest_Mu"+debug+".root", "RECREATE");

    DYAnalyzer *analyzer = new DYAnalyzer("Mu50");

    // -- For efficiency SF -- //
    analyzer->SetupEfficiencyScaleFactor_BtoF();
    analyzer->SetupEfficiencyScaleFactor_GtoH();

    // -- For W+Jets estimation from Fake Rate -- //
    analyzer->SetupFRvalues(Dir+"FakeRate_muon.root", "sigCtrl_template");

    TH1D *h_mass_test[_EndOf_SingleMuon_Normal];
    TH1D *h_mass_test_SS[_EndOf_SingleMuon_Normal];

    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();
        Int_t nPass = 0;
        Double_t avgFRweight = 0;

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_mass = new TH1D("h_mass_"+Mgr.Procname[pr], "h_mass"+Mgr.Procname[pr], binnum, massbins); h_mass->Sumw2();
        TH1D* h_mass_SS = new TH1D("h_mass_SS_"+Mgr.Procname[pr], "h_mass_SS_"+Mgr.Procname[pr], binnum, massbins); h_mass_SS->Sumw2();
        TH1D* h_mass_forFit = new TH1D("h_mass_forFit_"+Mgr.Procname[pr], "h_mass_forFit_"+Mgr.Procname[pr], 37,15,200); h_mass_forFit->Sumw2();
        TH1D* h_mass_SS_forFit = new TH1D("h_mass_SS_forFit_"+Mgr.Procname[pr], "h_mass_SS_forFit_"+Mgr.Procname[pr], 37,15,200); h_mass_SS_forFit->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX_"+Mgr.Procname[pr], "h_nVTX"+Mgr.Procname[pr], 50, 0, 50); h_nVTX->Sumw2();
        TH1D* h_eta_diff = new TH1D("h_eta_diff_"+Mgr.Procname[pr], "h_eta_diff_"+Mgr.Procname[pr], 50, 0, 5); h_eta_diff->Sumw2();
        TH1D* h_phi_diff = new TH1D("h_phi_diff_"+Mgr.Procname[pr], "h_phi_diff_"+Mgr.Procname[pr], 63, -3.15, 3.15); h_phi_diff->Sumw2();
        h_mass_test[pr] = new TH1D("h_test_"+Mgr.Procname[pr], "", binnum, massbins); h_mass_test[pr]->Sumw2();
        h_mass_test[pr]->SetDirectory(0);
        h_mass_test_SS[pr] = new TH1D("h_test_SS_"+Mgr.Procname[pr], "", binnum, massbins); h_mass_test_SS[pr]->Sumw2();
        h_mass_test_SS[pr]->SetDirectory(0);

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t evt_weight;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");

        if (type == 0) // Mu50 files
        {
            chain->Add(Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
            if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        }
        else if (type == 1) // IsoMu24 files
        {
            chain->Add(Dir+"SelectedForWJETest_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
            if (DEBUG == kTRUE) cout << Dir+"SelectedForWJETest_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        }
        else if (type == 2) // Triggerless files
        {
            chain->Add(Dir+"SelectedForBKGest_Mu_Triggerless"+Mgr.Procname[Mgr.CurrentProc]+".root");
            if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_Mu_Triggerless"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        }
        else
        {
            cout << "Wrong type! Select 0 (Mu50), 1 (IsoMu24) or 2 (Triggerless)" << endl;
            return;
        }
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("relPFiso", 1);
        chain->SetBranchStatus("TRKiso", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
//        chain->SetBranchStatus("evt_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);
        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("charge", &charge);
        chain->SetBranchAddress("relPFiso", &relPFiso);
        chain->SetBranchAddress("TRKiso", &TRKiso);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
//        chain->SetBranchAddress("evt_weight", &evt_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;


        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // W+Jets selection
            if (p_T->size() != 2) continue;
//            if (charge->at(0) == charge->at(1)) continue;
            if (relPFiso->at(0) < 0.15 && relPFiso->at(1) < 0.15) continue;
            if (relPFiso->at(0) > 0.15 && relPFiso->at(1) > 0.15) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
//            if (p_T->at(0) < 2 || p_T->at(1) < 2) continue;

            nPass++;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << charge->at(0) << " " << relPFiso->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << charge->at(1) << " " << relPFiso->at(1) << endl;

            if (DEBUG == kTRUE)
            {
                if (nPass >= 5) break;
                cout << "Evt " << i << endl;
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector mu1, mu2;
            mu1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Mu);
            mu2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Mu);

            Double_t mass = (mu1+mu2).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if (Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t weight1 = 0, weight2 = 0, effweight = 1;
            if (Mgr.isMC == kTRUE)
            {
                weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF(mu1, mu2);
                weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH(mu1, mu2);
                if (weight1<0 || weight2<0)
                    cout << "weights: " << weight1 << " " << weight2 << "   Pt: " << p_T->at(0) << " " << mu1.Pt() << " " << p_T->at(1)
                         << " " << mu2.Pt() << "   Eta: " << eta->at(0) << " " << mu1.Eta() << " " << eta->at(1) << " " << mu2.Eta() << endl;
                effweight = (Lumi_BtoF * weight1 + Lumi_GtoH * weight2) / Lumi;
            }

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            if (DEBUG == kTRUE) cout << "Eff weight: " << effweight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight <<
                                        "\tTop pT weight: " << TopPtWeight << endl;

            // -- FR WEIGHTS -- //
            Double_t FRweight = 1;
            if (relPFiso->at(0) >= 0.15 && relPFiso->at(1) >= 0.15) // Both failing -- skip
            {
                cout << "Both failed" << endl; // Comment this when not needed
                continue;
            }
            else // Only one failing
            {
                Double_t FR;
                if (relPFiso->at(0) > relPFiso->at(1)) // First fails, second passes
                    FR = analyzer->FakeRate(p_T->at(0), eta->at(0));
                else // Second fails, first passes
                    FR = analyzer->FakeRate(p_T->at(1), eta->at(1));
                FRweight = FR / (1 - FR);
                if (DEBUG == kTRUE) cout << "FR = " << FR << "   FRweight = " << FRweight << endl;
                avgFRweight += FRweight;
            }

            // -- Normalization -- //
            Double_t TotWeight = 1;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            // -- Histogram filling -- //
            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_eta_diff->Fill(fabs(eta->at(0)-eta->at(1)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            h_phi_diff->Fill(phi->at(0)-phi->at(1), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
            if (charge->at(0) != charge->at(1))
            {
                h_mass->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_forFit->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_test[pr]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            else
            {
                h_mass_SS->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_SS_forFit->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight * FRweight);
                h_mass_test_SS[pr]->Fill(mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }


        }// End of event iteration

        cout << "\t " << nPass << " events have passed the selection." << endl;
        cout << "\t Average FR weight: " << avgFRweight / NEvents << endl;
        if(Mgr.isMC == kTRUE)
        {
            cout << "\t *** Cross section: " << Mgr.Xsec[0] << endl;
            cout << "\t *** Sum of weights: " << Mgr.Wsum[0] << endl;
            printf("\t *** Normalization factor: %.8f\n\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);
        }

        f->cd();
        cout << "\tWriting into file...";

        h_mass->Write();
        h_mass_forFit->Write();
        h_mass_SS->Write();
        h_mass_SS_forFit->Write();
        h_nVTX->Write();
        h_eta_diff->Write();
        h_phi_diff->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TCanvas *c_data = new TCanvas("data","data", 800, 800);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_C]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_D]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_E]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_F]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_G]);
    h_mass_test[_SingleMuon_B]->Add(h_mass_test[_SingleMuon_H]);
    h_mass_test[_SingleMuon_B]->Draw("hist");
    c_data->SetLogx();
    c_data->Update();

    TCanvas *c_data_SS = new TCanvas("data_SS","data SS", 800, 800);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_C]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_D]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_E]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_F]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_G]);
    h_mass_test_SS[_SingleMuon_B]->Add(h_mass_test_SS[_SingleMuon_H]);
    h_mass_test_SS[_SingleMuon_B]->Draw("hist");
    c_data_SS->SetLogx();
    c_data_SS->Update();

    TCanvas *c_DY = new TCanvas("dy","dy", 800, 800);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_50to100]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_100to200]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_200to400]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_400to500]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_500to700]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_700to800]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_800to1000]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_1000to1500]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_1500to2000]);
    h_mass_test[_DY_10to50]->Add(h_mass_test[_DY_2000to3000]);
    h_mass_test[_DY_10to50]->Draw("hist");
    c_DY->SetLogx();
    c_DY->Update();

    TCanvas *c_DY_SS = new TCanvas("dy_SS","dy SS", 800, 800);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_50to100]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_100to200]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_200to400]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_400to500]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_500to700]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_700to800]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_800to1000]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_1000to1500]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_1500to2000]);
    h_mass_test_SS[_DY_10to50]->Add(h_mass_test_SS[_DY_2000to3000]);
    h_mass_test_SS[_DY_10to50]->Draw("hist");
    c_DY_SS->SetLogx();
    c_DY_SS->Update();

    TCanvas *c_tt = new TCanvas("ttbar","ttbar", 800, 800);
    h_mass_test[_ttbar]->Add(h_mass_test[_ttbar_700to1000]);
    h_mass_test[_ttbar]->Add(h_mass_test[_ttbar_1000toInf]);
    h_mass_test[_ttbar]->Draw("hist");
    c_tt->SetLogx();
    c_tt->Update();

    TCanvas *c_tt_SS = new TCanvas("ttbar_SS","ttbar SS", 800, 800);
    h_mass_test_SS[_ttbar]->Add(h_mass_test_SS[_ttbar_700to1000]);
    h_mass_test_SS[_ttbar]->Add(h_mass_test_SS[_ttbar_1000toInf]);
    h_mass_test_SS[_ttbar]->Draw("hist");
    c_tt_SS->SetLogx();
    c_tt_SS->Update();

    TCanvas *c_wjets = new TCanvas("wjets","wjets", 800, 800);
    h_mass_test[_WJets]->Add(h_mass_test[_WJets_ext2v5]);
    h_mass_test[_WJets]->Draw("hist");
    c_wjets->SetLogx();
    c_wjets->Update();

    TCanvas *c_wjets_SS = new TCanvas("wjets_SS","wjets SS", 800, 800);
    h_mass_test_SS[_WJets]->Add(h_mass_test_SS[_WJets_ext2v5]);
    h_mass_test_SS[_WJets]->Draw("hist");
    c_wjets_SS->SetLogx();
    c_wjets_SS->Update();

    TH1D *h_wjets_SS_scaled = (TH1D*)h_mass_test_SS[_WJets]->Clone("h_wjets_SS_scaled");
    h_wjets_SS_scaled->Rebin(8);
    h_wjets_SS_scaled->Scale(3);
    h_wjets_SS_scaled->SetLineColor(kRed);
    h_wjets_SS_scaled->SetDirectory(0);
    TH1D *h_wjets_OS_compare = (TH1D*)h_mass_test[_WJets]->Clone("h_wjets_OS_compare");
    h_wjets_OS_compare->SetDirectory(0);
    h_wjets_OS_compare->Rebin(8);
    TCanvas *c_wjets_compare = new TCanvas("wjets_compare", "WJets OS vs SSx3", 800, 800);
    h_wjets_OS_compare->Draw("hist");
    h_wjets_SS_scaled->Draw("histsame");
    c_wjets_compare->SetLogx();
    c_wjets_compare->Update();
    cout << "Ratio of W+Jets opposite-sign and same-sign integrals: " << h_mass_test[_WJets]->Integral()/h_mass_test_SS[_WJets]->Integral() << endl;

    f->Close();
    if (!f->IsOpen()) cout << "File " << Dir+"WJETest_Mu"+debug+".root" << " has been closed successfully.\n" << endl;
    else cout << "FILE " << Dir+"WJETest_Mu"+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
    */
} // End of Mu_WJET_Tester()


void EMu_WJET_Tester (TString FRtype, Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/EMu/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    DYAnalyzer *analyzer = new DYAnalyzer("IsoMu24_OR_IsoTkMu24");

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues_ele("/media/sf_DATA/FR/Electron/FakeRate_electron.root", "subtract");
    analyzer->SetupFRvalues("/media/sf_DATA/FR/Muon/FakeRate_muon.root", "sigCtrl_template");

    TH1D *h_mass, *h_mass_compare;
    TH1D *h_mass_1FR = new TH1D("h_mass_1FR", "h_mass_1FR", binnum, massbins);

    for (Process_t pr=_WJets; pr<=_WJets_ext2v5; pr=next(pr))
    {
        Mgr.SetProc(pr);

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

        std::vector<double> *e_p_T = new std::vector<double>;
        std::vector<double> *e_eta = new std::vector<double>;
        std::vector<double> *e_phi = new std::vector<double>;
        std::vector<double> *e_etaSC = new std::vector<double>;
        std::vector<int> *e_charge = new std::vector<int>;
        std::vector<double> *e_Full5x5_SigmaIEtaIEta = new std::vector<double>;
        std::vector<double> *e_dEtaInSeed = new std::vector<double>;
        std::vector<double> *e_dPhiIn = new std::vector<double>;
        std::vector<double> *e_HoverE = new std::vector<double>;
        std::vector<double> *e_InvEminusInvP = new std::vector<double>;
        std::vector<int> *e_mHits = new std::vector<int>;
        std::vector<int> *e_passConvVeto = new std::vector<int>;
        std::vector<double> *e_relPFiso_Rho = new std::vector<double>;
        std::vector<int> *e_passMediumID = new std::vector<int>;
        std::vector<double> *mu_p_T = new std::vector<double>;
        std::vector<double> *mu_eta = new std::vector<double>;
        std::vector<double> *mu_phi = new std::vector<double>;
        std::vector<int> *mu_charge = new std::vector<int>;
        std::vector<double> *mu_relPFiso_dBeta = new std::vector<double>;
//        Double_t e_p_T, e_eta, e_phi,*e_etaSC;
//        Int_t e_passMediumID;
//        Double_t mu_p_T, mu_eta, mu_phi;
//        Double_t mu_relPFiso_dBeta;
        Int_t nPU;
        Double_t PVz;
        Double_t gen_weight;
        Double_t prefiring_weight;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForBKGest_EMu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_EMu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("e_p_T", 1);
        chain->SetBranchStatus("e_eta", 1);
        chain->SetBranchStatus("e_phi", 1);
        chain->SetBranchStatus("e_etaSC", 1);
        chain->SetBranchStatus("e_charge", 1);
        chain->SetBranchStatus("e_Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("e_dEtaInSeed", 1);
        chain->SetBranchStatus("e_dPhiIn", 1);
        chain->SetBranchStatus("e_HoverE", 1);
        chain->SetBranchStatus("e_InvEminusInvP", 1);
        chain->SetBranchStatus("e_mHits", 1);
        chain->SetBranchStatus("e_passConvVeto", 1);
        chain->SetBranchStatus("e_relPFiso_Rho", 1);
        chain->SetBranchStatus("e_passMediumID", 1);
        chain->SetBranchStatus("mu_p_T", 1);
        chain->SetBranchStatus("mu_eta", 1);
        chain->SetBranchStatus("mu_phi", 1);
        chain->SetBranchStatus("mu_charge", 1);
        chain->SetBranchStatus("mu_relPFiso_dBeta", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchAddress("e_p_T", &e_p_T);
        chain->SetBranchAddress("e_eta", &e_eta);
        chain->SetBranchAddress("e_phi", &e_phi);
        chain->SetBranchAddress("e_etaSC", &e_etaSC);
        chain->SetBranchAddress("e_charge", &e_charge);
        chain->SetBranchAddress("e_Full5x5_SigmaIEtaIEta", &e_Full5x5_SigmaIEtaIEta);
        chain->SetBranchAddress("e_dEtaInSeed", &e_dEtaInSeed);
        chain->SetBranchAddress("e_dPhiIn", &e_dPhiIn);
        chain->SetBranchAddress("e_HoverE", &e_HoverE);
        chain->SetBranchAddress("e_InvEminusInvP", &e_InvEminusInvP);
        chain->SetBranchAddress("e_mHits", &e_mHits);
        chain->SetBranchAddress("e_passConvVeto", &e_passConvVeto);
        chain->SetBranchAddress("e_relPFiso_Rho", &e_relPFiso_Rho);
        chain->SetBranchAddress("e_passMediumID", &e_passMediumID);
        chain->SetBranchAddress("mu_p_T", &mu_p_T);
        chain->SetBranchAddress("mu_eta", &mu_eta);
        chain->SetBranchAddress("mu_phi", &mu_phi);
        chain->SetBranchAddress("mu_charge", &mu_charge);
        chain->SetBranchAddress("mu_relPFiso_dBeta", &mu_relPFiso_dBeta);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;


        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // WJets selection
            if (e_passMediumID->at(0) == 1 && mu_relPFiso_dBeta->at(0) < 0.15) continue;
            if (e_passMediumID->at(0) == 0 && mu_relPFiso_dBeta->at(0) >= 0.15) continue;
            if (e_p_T->at(0) < 17 || mu_p_T->at(0) < 17) continue;
            if (e_p_T->at(0) < 28 && mu_p_T->at(0) < 28) continue;
            if (fabs(e_etaSC->at(0)) > 1.4442 && fabs(e_etaSC->at(0)) < 1.566) continue;
            if (fabs(e_etaSC->at(0)) < 1.4442 && (e_Full5x5_SigmaIEtaIEta->at(0) >= 0.013 || e_HoverE->at(0) >= 0.13 ||
                                                  fabs(e_dEtaInSeed->at(0)) >= 0.01 || fabs(e_dPhiIn->at(0)) >= 0.07 || e_mHits->at(0) > 1)) continue;
            if (fabs(e_etaSC->at(0)) > 1.566 && (e_Full5x5_SigmaIEtaIEta->at(0) >= 0.035 || e_HoverE->at(0) >= 0.13 || e_mHits->at(0) > 1)) continue;

            if (e_p_T->at(0) != e_p_T->at(0)) cout << e_p_T->at(0) << " " << e_eta->at(0) << " " << e_phi->at(0) << " " << e_charge->at(0) << " " << e_passMediumID->at(0) << endl;
            if (mu_p_T->at(0) != mu_p_T->at(0)) cout << mu_p_T->at(0) << " " << mu_eta->at(0) << " " << mu_phi->at(0) << " " << mu_charge->at(0) << " " << mu_relPFiso_dBeta->at(0) << endl;

            if (DEBUG == kTRUE)
            {
                cout << "Evt " << i << endl;
                cout << "Electron p_T = " << e_p_T->at(0);
                cout << "\teta = " << e_eta->at(0);
                cout << "\tphi = " << e_phi->at(0) << endl;
                cout << "\tpassMediumID = " << e_passMediumID->at(0) << endl;
                cout << "Muon p_T = " << mu_p_T->at(0);
                cout << "\teta = " << mu_eta->at(0);
                cout << "\tphi = " << mu_phi->at(0) << endl;
                cout << "\trelPFiso = " << mu_relPFiso_dBeta->at(0) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector ele, mu;
            ele.SetPtEtaPhiM(e_p_T->at(0), e_eta->at(0), e_phi->at(0), M_Elec);
            mu.SetPtEtaPhiM(mu_p_T->at(0), mu_eta->at(0), mu_phi->at(0), M_Mu);
            Double_t mass = (ele+mu).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- PVz weights -- //
            Double_t PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = prefiring_weight;

            if (DEBUG == kTRUE) cout << "PU weight: " << PUWeight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight << endl;

            // -- Normalization -- //
            Double_t TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl << endl;

            // -- FR WEIGHTS -- //
            Double_t FRweight=1, FR=-1;
            if (e_passMediumID->at(0) == 0 && mu_relPFiso_dBeta->at(0) < 0.15) // Electron fails, muon passes
                FR = analyzer->FakeRate_ele(e_p_T->at(0), e_etaSC->at(0));
            else if (e_passMediumID->at(0) == 1 && mu_relPFiso_dBeta->at(0) >= 0.15) // Muon fails, electron passes
                FR = analyzer->FakeRate(mu_p_T->at(0), mu_eta->at(0));
            FRweight = FR / (1 - FR);
            h_mass_1FR->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * FRweight);

        }// End of event iteration

        cout << " Finished.\n" << endl;

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TFile *f_orig = new TFile("/media/sf_DATA/SelectedEMu/Histos/Hist_SelectedEMu_Bkg_Full.root", "READ");
    f_orig->GetObject("h_emu_mass_SelectedEMu_WJets_Full", h_mass);
    h_mass->SetDirectory(0);
    f_orig->Close();

    TCanvas *c_comparison = new TCanvas("FR test","FR test (W+Jets)", 800, 800);
    c_comparison->SetRightMargin(0.05);
    c_comparison->SetTopMargin(0.05);
    c_comparison->SetBottomMargin(0.12);
    c_comparison->SetLeftMargin(0.13);
    h_mass_1FR->SetDirectory(0);
    h_mass->SetStats(0);
    h_mass_1FR->SetStats(0);
    h_mass->SetMarkerStyle(kFullDotLarge);
    h_mass_1FR->SetMarkerStyle(kFullSquare);
    h_mass->SetLineColor(kBlack);
    h_mass_1FR->SetLineColor(kRed);
    h_mass_1FR->SetMarkerColor(kRed);
    h_mass->GetXaxis()->SetTitle("m_{e#mu} [GeV/c^{2}]");
    h_mass->GetXaxis()->SetNoExponent(1);
    h_mass->GetXaxis()->SetMoreLogLabels(1);
    h_mass->GetXaxis()->SetTitleOffset(1);
    h_mass->GetXaxis()->SetTitleSize(0.05);
    h_mass->GetXaxis()->SetLabelSize(0.04);
    h_mass->GetYaxis()->SetTitle("Number of events");
    h_mass->GetYaxis()->SetTitleSize(0.05);
    h_mass->GetYaxis()->SetTitleOffset(1.25);
    h_mass->GetYaxis()->SetLabelSize(0.04);
    h_mass->GetYaxis()->SetRangeUser(1, 1e4);
    TLegend *legend = new TLegend(0.5, 0.7, 0.95, 0.95);
    legend->AddEntry(h_mass, "W+Jets MC (2 pass)", "lp");
    legend->AddEntry(h_mass_1FR, "W+Jets MC (1 fail, FR applied)", "lp");
    h_mass->Draw();
    h_mass_1FR->Draw("same");
    legend->Draw();
    c_comparison->SetLogx();
    c_comparison->SetLogy();
    c_comparison->SetGridx();
    c_comparison->SetGridy();
    c_comparison->Update();

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
} // End of EMu_WJET_Tester()
