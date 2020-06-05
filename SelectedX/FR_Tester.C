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
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    DYAnalyzer *analyzer = new DYAnalyzer("Mu20_OR_Mu27_OR_Mu50");

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues(Dir+"FakeRate_muon.root", FRtype);

    TH1D *h_mass_2FR = new TH1D("h_mass_2FR", "", 10, 0, 1000);
    TH1D *h_mass_1FR = new TH1D("h_mass_1FR", "", 10, 0, 1000);

    for (Process_t pr=_QCDMuEnriched_20to30; pr<_EndOf_QCDMuEnriched_Normal; pr=next(pr))
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
        analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        Int_t nPU;
        Double_t PVz;
        Double_t gen_weight;
        Double_t prefiring_weight;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("charge", 1);
        chain->SetBranchStatus("relPFiso", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);

        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("charge", &charge);
        chain->SetBranchAddress("relPFiso", &relPFiso);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);

        Int_t NEvents = chain->GetEntries();
        if (DEBUG) NEvents = 15;
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;


        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (!DEBUG) bar.Draw(i);

            // QCD selection
            if (p_T->size() != 2) continue;
            if (relPFiso->at(0) < 0.15 && relPFiso->at(1) < 0.15) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;
            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << relPFiso->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << relPFiso->at(1) << endl;

            if (DEBUG == kTRUE)
            {
                cout << "\n===============================" << endl;
                cout << "Evt " << i << endl;
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "\tcharge[0] = " << charge->at(0) << endl;
                cout << "\tPFiso[0] = " << relPFiso->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\tcharge[1] = " << charge->at(1) << endl;
                cout << "\tPFiso[1] = " << relPFiso->at(1) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector mu1, mu2;
            mu1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Mu);
            mu2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Mu);
            Double_t mass = (mu1+mu2).M();

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- PVz weights -- //
            Double_t PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;

            if (DEBUG == kTRUE) cout << "PU weight: " << PUWeight << "\tPVz weight: " << PVzWeight << "\nL1 weight:" << L1weight << endl;

            // -- Normalization -- //
            Double_t TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Normalization weight " << TotWeight << endl;

            // -- FR WEIGHTS -- //
            if (relPFiso->at(0) < 0.15 && relPFiso->at(1) >= 0.15) // First passing
            {
                Double_t FRweight, FR;
                FR = analyzer->FakeRate(p_T->at(1), eta->at(1));
                FRweight = FR / (1 - FR);
                h_mass_1FR->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * FRweight);
                if (DEBUG == kTRUE) cout << "Total weight " << TotWeight * PUWeight * PVzWeight * L1weight * FRweight << endl;
            }
            else if (relPFiso->at(0) >= 0.15 && relPFiso->at(1) < 0.15) // Second passing
            {
                Double_t FRweight, FR;
                FR = analyzer->FakeRate(p_T->at(0), eta->at(0));
                FRweight = FR / (1 - FR);
                h_mass_1FR->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * FRweight);
                if (DEBUG == kTRUE) cout << "Total weight " << TotWeight * PUWeight * PVzWeight * L1weight * FRweight << endl;
            }
            else // Both failing
            {
                Double_t FRweight = 1;
                Double_t FR1, FR2;
                FR1 = analyzer->FakeRate(p_T->at(0), eta->at(0));
                FR2 = analyzer->FakeRate(p_T->at(1), eta->at(1));
                FRweight = (FR1 / (1 - FR1)) * (FR2 / (1 - FR2));
                if (DEBUG == kTRUE) cout << "FR1 = " << FR1 << "   FR2 = " << FR2 << "   FRweight = " << FRweight << endl;
                h_mass_2FR->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * FRweight);
                if (DEBUG == kTRUE) cout << "Total weight " << TotWeight * PUWeight * PVzWeight * L1weight * FRweight << endl;
            }

        }// End of event iteration

        cout << " Finished.\n" << endl;
        cout << "===========================================================\n" << endl;
        if (DEBUG && pr == _QCDMuEnriched_30to50) break;
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
    h_mass_1FR->GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
    h_mass_1FR->GetXaxis()->SetNoExponent(1);
    h_mass_1FR->GetXaxis()->SetMoreLogLabels(1);
    h_mass_1FR->GetXaxis()->SetTitleOffset(1);
    h_mass_1FR->GetXaxis()->SetTitleSize(0.05);
    h_mass_1FR->GetXaxis()->SetLabelSize(0.04);
    h_mass_1FR->GetYaxis()->SetTitle("Number of events");
    h_mass_1FR->GetYaxis()->SetTitleSize(0.05);
    h_mass_1FR->GetYaxis()->SetTitleOffset(1.25);
    h_mass_1FR->GetYaxis()->SetLabelSize(0.04);
    h_mass_1FR->GetYaxis()->SetRangeUser(1e-2, 1e5);
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

} // End of Mu_QCD_Tester()


void Mu_WJET_Tester (TString FRtype, Bool_t DEBUG)
{
    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    FileMgr Mgr;

    TFile *f;
    TString Dir = "/media/sf_DATA/FR/Muon/";
    TString debug = "";
    if (DEBUG == kTRUE) debug = "_DEBUG";

    DYAnalyzer *analyzer = new DYAnalyzer("Mu20_OR_Mu27_OR_Mu50");

    // -- For QCD estimation from Fake Rate -- //
    analyzer->SetupFRvalues(Dir+"FakeRate_muon.root", FRtype);

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
        analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        Int_t nPU;
        Double_t PVz;
        Double_t gen_weight;
        Double_t prefiring_weight;

        TChain *chain = new TChain("FRTree");
        chain->Add(Dir+"SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForBKGest_Mu_Triggerless_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;

        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("charge", 1);
        chain->SetBranchStatus("relPFiso", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);

        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("charge", &charge);
        chain->SetBranchAddress("relPFiso", &relPFiso);
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
            if (p_T->size() != 2) continue;
            if (relPFiso->at(0) < 0.15 && relPFiso->at(1)  < 0.15) continue;
            if (relPFiso->at(0) >= 0.15 && relPFiso->at(1) >= 0.15) continue;
            if (p_T->at(0) < 17 || p_T->at(1) < 17) continue;
            if (p_T->at(0) < 28 && p_T->at(1) < 28) continue;

            if (p_T->at(0) != p_T->at(0)) cout << p_T->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << relPFiso->at(0) << endl;
            if (p_T->at(1) != p_T->at(1)) cout << p_T->at(1) << " " << eta->at(1) << " " << phi->at(1) << " " << relPFiso->at(1) << endl;

            if (DEBUG == kTRUE)
            {
                cout << "Evt " << i << endl;
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[0] = " << p_T->at(0);
                cout << "\teta[0] = " << eta->at(0);
                cout << "\tphi[0] = " << phi->at(0) << endl;
                cout << "\tcharge[0] = " << charge->at(0) << endl;
                cout << "\tPFiso[0] = " << relPFiso->at(0) << endl;
                cout << "p_T[1] = " << p_T->at(1);
                cout << "\teta[1] = " << eta->at(1);
                cout << "\tphi[1] = " << phi->at(1) << endl;
                cout << "\tcharge[1] = " << charge->at(1) << endl;
                cout << "\tPFiso[1] = " << relPFiso->at(1) << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector mu1, mu2;
            mu1.SetPtEtaPhiM(p_T->at(0), eta->at(0), phi->at(0), M_Mu);
            mu2.SetPtEtaPhiM(p_T->at(1), eta->at(1), phi->at(1), M_Mu);
            Double_t mass = (mu1+mu2).M();

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
            if (DEBUG == kTRUE) cout << "Normalization weight " << TotWeight << endl << endl;

            // -- FR WEIGHTS -- //
            Double_t FRweight=1, FR=-1;
            if (relPFiso->at(0) < 0.15 && relPFiso->at(1) >= 0.15) // First passing
                FR = analyzer->FakeRate(p_T->at(1), eta->at(1));
            else if (relPFiso->at(0) >= 0.15 && relPFiso->at(1) < 0.15) // First failing
                FR = analyzer->FakeRate(p_T->at(0), eta->at(0));
            FRweight = FR / (1 - FR);
            h_mass_1FR->Fill(mass, TotWeight * PUWeight * PVzWeight * L1weight * FRweight);

        }// End of event iteration

        cout << " Finished.\n" << endl;

        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    TFile *f_orig = new TFile("/media/sf_DATA/SelectedMuMu/Histos/Hist_SelectedMuMu_Bkg_Full.root", "READ");
    f_orig->GetObject("h_mass_SelectedMuMu_WJets_Full", h_mass);
    h_mass->SetDirectory(0);
    f_orig->Close();
    TFile *f_est = new TFile("/media/sf_DATA/SelectedMuMu/Histos/EstWJets_MuMu.root", "READ");
    f_est->GetObject("h_WJET_fit_SS", h_mass_compare);
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
    h_mass->GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
    h_mass->GetXaxis()->SetNoExponent(1);
    h_mass->GetXaxis()->SetMoreLogLabels(1);
    h_mass->GetXaxis()->SetTitleOffset(1);
    h_mass->GetXaxis()->SetTitleSize(0.05);
    h_mass->GetXaxis()->SetLabelSize(0.04);
    h_mass->GetYaxis()->SetTitle("Number of events");
    h_mass->GetYaxis()->SetTitleSize(0.05);
    h_mass->GetYaxis()->SetTitleOffset(1.25);
    h_mass->GetYaxis()->SetLabelSize(0.04);
    h_mass->GetYaxis()->SetRangeUser(1e-2, 1e3);
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

        Double_t e_p_T, e_eta, e_etaSC, e_phi;
        Int_t e_charge;
        Double_t e_Full5x5_SigmaIEtaIEta, e_dEtaInSeed, e_dPhiIn, e_HoverE, e_InvEminusInvP;
        Double_t e_relPFiso_dBeta, e_relPFiso_Rho;
        Double_t e_chIso03, e_nhIso03, e_phIso03, e_ChIso03FromPU;
        Int_t e_mHits, e_passConvVeto, e_passMediumID;
        Double_t mu_p_T, mu_eta, mu_phi;
        Int_t mu_charge;
        Double_t mu_relPFiso_dBeta;
        std::vector<int> *trig_fired = new std::vector<int>;
        std::vector<double> *trig_pT = new std::vector<double>;
        std::vector<int> *prescale_factor = new std::vector<int>;
        Double_t MET_pT, MET_phi;
        Int_t nPU, nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;
        Int_t runNum;
        Int_t lumiBlock;

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
            if (e_passMediumID == 1 && mu_relPFiso_dBeta < 0.15) continue;
            if (e_passMediumID == 0 && mu_relPFiso_dBeta >= 0.15) continue;
            if (e_p_T < 17 || mu_p_T < 17) continue;
            if (e_p_T < 28 && mu_p_T < 28) continue;
            if (fabs(e_etaSC) > 1.4442 && fabs(e_etaSC) < 1.566) continue;
            if (fabs(e_etaSC) < 1.4442 && (e_Full5x5_SigmaIEtaIEta >= 0.013 || e_HoverE >= 0.13 ||
                                           fabs(e_dEtaInSeed) >= 0.01 || fabs(e_dPhiIn) >= 0.07 || e_mHits > 1)) continue;
            if (fabs(e_etaSC) > 1.566 && (e_Full5x5_SigmaIEtaIEta >= 0.035 || e_HoverE >= 0.13 || e_mHits > 1)) continue;

            if (e_p_T != e_p_T) cout << e_p_T << " " << e_eta << " " << e_phi << " " << e_charge << " " << e_passMediumID << endl;
            if (mu_p_T != mu_p_T) cout << mu_p_T << " " << mu_eta << " " << mu_phi << " " << mu_charge << " " << mu_relPFiso_dBeta << endl;

            if (DEBUG == kTRUE)
            {
                cout << "Evt " << i << endl;
                cout << "Electron p_T = " << e_p_T;
                cout << "\teta = " << e_eta;
                cout << "\tphi = " << e_phi << endl;
                cout << "\tpassMediumID = " << e_passMediumID << endl;
                cout << "Muon p_T = " << mu_p_T;
                cout << "\teta = " << mu_eta;
                cout << "\tphi = " << mu_phi << endl;
                cout << "\trelPFiso = " << mu_relPFiso_dBeta << endl;
                cout << "\nPVz = " << PVz << endl;
            }

            TLorentzVector ele, mu;
            ele.SetPtEtaPhiM(e_p_T, e_eta, e_phi, M_Elec);
            mu.SetPtEtaPhiM(mu_p_T, mu_eta, mu_phi, M_Mu);
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
            if (e_passMediumID == 0 && mu_relPFiso_dBeta < 0.15) // Electron fails, muon passes
                FR = analyzer->FakeRate_ele(e_p_T, e_etaSC);
            else if (e_passMediumID == 1 && mu_relPFiso_dBeta >= 0.15) // Muon fails, electron passes
                FR = analyzer->FakeRate(mu_p_T, mu_eta);
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
    TFile *f_est = new TFile("/media/sf_DATA/SelectedEMu/Histos/EstWJets_EMu.root", "READ");
    f_est->GetObject("h_WJET_est", h_mass_compare);
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
    h_mass->SetMarkerStyle(kFullDotLarge);
    h_mass_1FR->SetMarkerStyle(kFullSquare);
    h_mass->SetLineColor(kBlack);
    h_mass_1FR->SetLineColor(kRed);
    h_mass_1FR->SetMarkerColor(kRed);
    h_mass_compare->SetStats(0);
    h_mass_compare->SetMarkerStyle(33);
    h_mass_compare->SetLineColor(kGreen+2);
    h_mass_compare->SetMarkerColor(kGreen+2);
    h_mass_compare->SetMarkerSize(1.5);
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
    legend->AddEntry(h_mass_compare, "W+Jets est (data-driven)", "lp");
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
} // End of EMu_WJET_Tester()
