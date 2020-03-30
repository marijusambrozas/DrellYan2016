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

void E_FR_HistMaker (Bool_t DEBUG);
void Mu_FR_HistMaker (Bool_t DEBUG);

void Mu_QCD_HistMaker (Bool_t DEBUG, Int_t type);
void Mu_WJET_HistMaker (Bool_t DEBUG, Int_t type);

// -- Drell-Yan mass bins -- //
const Int_t binnum = 43;
const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126,
                               133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

// -- Mass bins for background estimation -- //
const Int_t binnum2 = 86;
const Double_t massbins2[87] = {15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62, 64,
                                66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5,
                                115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200,
                                210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830,
                                915, 1000, 1250, 1500, 2250, 3000};


void FR_HistMaker_alt (TString WhichX = "", Int_t type=1)
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
    if (whichX.Contains("MU"))
    {
        Xselected++;
        if (whichX.Contains("QCD"))
        {
            cout << "\n*****  Mu_QCD_HistMaker(" << type << ")  *****" << endl;
            Mu_QCD_HistMaker(DEBUG, type);
        }
        else if (whichX.Contains("W") && whichX.Contains("JET"))
        {
            cout << "\n*****  Mu_WJET_HistMaker(" << type << ")  *****" << endl;
            Mu_WJET_HistMaker(DEBUG, type);
        }
        else
        {
            cout << "\n*****  Mu_FR_HistMaker  *****" << endl;
            Mu_FR_HistMaker(DEBUG);
        }
    }
    else if (whichX.Contains("E"))
    {
        Xselected++;
        cout << "\n*****    E_FR_HistMaker    *****" << endl;
        E_FR_HistMaker(DEBUG);
    }

    if (Xselected == 0) cout << "Wrong arument! \nType in: >> .x FR_HistMaker.C+(\"whichX\", type)" << endl;

} // End of HistMaker()


/// ----------------------------- Electron Channel ------------------------------ ///
void E_FR_HistMaker (Bool_t DEBUG)
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

    UInt_t n2MC=0, n2Data=0;
    UInt_t n22175=0, n30175=0, n36175=0, n50175=0, n75175=0, n90175=0, n120175=0, n175=0;

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

        DYAnalyzer *analyzer = new DYAnalyzer("Photon_OR");

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For efficiency SF -- //
        analyzer->SetupEfficiencyScaleFactor_electron();

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_pT_uncorr = new TH1D("h_pT_uncorr", "h_pT_uncorr", 500, 0, 500); h_pT_uncorr->Sumw2();
        TH1D* h_pT = new TH1D("h_pT", "h_pT", 500, 0, 500); h_pT->Sumw2();
        TH1D* h_HLT_pT_uncorr = new TH1D("h_HLT_pT_uncorr", "h_HLT_pT_uncorr", 500, 0, 500); h_HLT_pT_uncorr->Sumw2();
        TH1D* h_HLT_pT = new TH1D("h_HLT_pT", "h_HLT_pT", 500, 0, 500); h_HLT_pT->Sumw2();
        TH1D* h_eta = new TH1D("h_eta", "h_eta", 48, -2.4, 2.4); h_eta->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *passMediumID = new std::vector<int>;
        std::vector<int> *trig_fired = new std::vector<int>;
        std::vector<int> *trig_matched = new std::vector<int>;
        std::vector<double> *trig_pT = new std::vector<double>;
        Double_t MET_pT, MET_phi;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;
//        Double_t prescale_factor;

        TChain *chain = new TChain("FRTree");

        chain->Add(Dir+"SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_E_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("passMediumID", 1);
        chain->SetBranchStatus("trig_fired", 1);
        chain->SetBranchStatus("trig_matched", 1);
        chain->SetBranchStatus("trig_pT", 1);
        chain->SetBranchStatus("MET_pT", 1);
        chain->SetBranchStatus("MET_phi", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
        chain->SetBranchStatus("gen_weight", 1);
        chain->SetBranchStatus("top_weight", 1);
        chain->SetBranchStatus("prefiring_weight", 1);
        chain->SetBranchStatus("prefiring_weight_up", 1);
        chain->SetBranchStatus("prefiring_weight_down", 1);
//        chain->SetBranchStatus("prescale_factor", 1);
        chain->SetBranchAddress("p_T", &p_T);
        chain->SetBranchAddress("eta", &eta);
        chain->SetBranchAddress("phi", &phi);
        chain->SetBranchAddress("passMediumID", &passMediumID);
        chain->SetBranchAddress("trig_fired", &trig_fired);
        chain->SetBranchAddress("trig_matched", &trig_matched);
        chain->SetBranchAddress("trig_pT", &trig_pT);
        chain->SetBranchAddress("MET_pT", &MET_pT);
        chain->SetBranchAddress("MET_phi", &MET_phi);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);
//        chain->SetBranchAddress("prescale_factor", &prescale_factor);

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
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[1] = " << p_T->at(0) << endl;
                cout << "eta[1] = " << eta->at(0) << endl;
                cout << "phi[1] = " << phi->at(0) << endl;
            }

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if(Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t effweight = 1;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE) TopPtWeight = top_weight;

            // -- Normalization -- //
            Double_t TotWeight = gen_weight;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl;

            if (Mgr.isMC == kTRUE && p_T->size() > 1) n2MC += TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight;
            if (Mgr.isMC == kFALSE && p_T->size() > 1) n2Data++;

            // For finding the leading electron
            TLorentzVector ele_lead;
            ele_lead.SetPtEtaPhiM(0, 0, 0, M_Elec);

            if (p_T->size() != passMediumID->size())
            {
                cout << "ERROR: vector sizes do not match!" << endl;
                break;
            }


            if (DEBUG == kTRUE)
            {
                cout << "Triggers:" << endl;
                for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
                {
                    cout << "Photon" << trig_fired->at(i_tr) << "  p_T: " << trig_pT->at(i_tr) << "  matched to " << trig_matched->at(i_tr) << endl;
                }
            }

            Double_t prescale_alt = 0;
//            Double_t tr_highest = -9999;
//            Int_t i_highest = -1;
//            for (Int_t i_tr=0; i_tr<(int)(trig_fired->size()); i_tr++)
//            {
//                if (trig_fired->at(i_tr) > tr_highest)
//                {
//                    tr_highest = trig_fired->at(i_tr);
//                    i_highest = i_tr;
//                }
//            }
//            if (i_highest < 0) continue;
//            if (i_highest > (int)(trig_fired->size()))
//            {
//                cout << "Highest pT trigger index is higher than a number of indices!  i_high=" << i_highest << "   size=" << trig_fired->size() << endl;
//                break;
//            }

            Int_t fired22=0, fired30=0, fired36=0, fired50=0, fired75=0, fired90=0, fired120=0, fired175=0;
            for (UInt_t i_tr=0; i_tr<trig_fired->size(); i_tr++)
            {
                if (trig_fired->at(i_tr) == 22) fired22++;
                if (trig_fired->at(i_tr) == 30) fired30++;
                if (trig_fired->at(i_tr) == 36) fired36++;
                if (trig_fired->at(i_tr) == 50) fired50++;
                if (trig_fired->at(i_tr) == 75) fired75++;
                if (trig_fired->at(i_tr) == 90) fired90++;
                if (trig_fired->at(i_tr) == 120) fired120++;
                if (trig_fired->at(i_tr) == 175) fired175++;
            }
            if (fired175 && !fired22 && !fired30 && !fired36 && !fired50 && !fired75 && !fired90 && !fired120) n175++;
            if (fired175 && fired120) n120175++;
            if (fired175 && fired90) n90175++;
            if (fired175 && fired75) n75175++;
            if (fired175 && fired50) n50175++;
            if (fired175 && fired36) n36175++;
            if (fired175 && fired30) n30175++;
            if (fired175 && fired22) n22175++;

            for (UInt_t i_ele=0; i_ele<p_T->size(); i_ele++)
            {
                if (p_T->at(i_ele) != p_T->at(i_ele))
                {
                    cout << p_T->at(i_ele) << " " << eta->at(i_ele) << " " << phi->at(i_ele) << endl;
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
//                        matched = 1;
//                        prescale_alt = analyzer->getPrescale_alt(trig_pT->at(i_tr));
//                        prescale_alt += analyzer->getPrescale(trig_fired->at(i_tr)+1); //BAD

                    }
                }
                if (matched22+matched30+matched36+matched50+matched75+matched90+matched120 > 1) continue;
                if (matched22==0 && matched30==0 && matched36==0 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0) continue;
                else if (matched22 == 1 && matched30 == 0 && matched36 == 0 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = 813./18621470.;
                    h_HLT_pT->Fill(trig_pT->at(i_22), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_22), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_pT->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_eta->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                }
                else if (matched30 == 1 && matched36 == 0 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = 3211./18621470.;
                    h_HLT_pT->Fill(trig_pT->at(i_30), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_30), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_pT->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_eta->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                }
                else if (matched36 == 1 && matched50 == 0 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = 6372./18621470.;
                    h_HLT_pT->Fill(trig_pT->at(i_36), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_36), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_pT->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_eta->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                }
                else if (matched50 == 1 && matched75 == 0 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = 12648./18621470.;
                    h_HLT_pT->Fill(trig_pT->at(i_50), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_50), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_pT->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_eta->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                }
                else if (matched75 == 1 && matched90 == 0 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = 63170./18621470.;
                    h_HLT_pT->Fill(trig_pT->at(i_75), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_75), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_pT->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_eta->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                }
                else if (matched90 == 1 && matched120 == 0 && matched175 == 0)
                {
                    prescale_alt = 126981./18621470.;
                    h_HLT_pT->Fill(trig_pT->at(i_90), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_90), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_pT->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_eta->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                }
                else if (matched120 == 1 && matched175 == 0)
                {
                    prescale_alt = 260278./18621470.;
                    h_HLT_pT->Fill(trig_pT->at(i_120), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_120), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_pT->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_eta->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                }
                else if (matched175 == 1)
                {
                    prescale_alt = 1.;
                    h_HLT_pT->Fill(trig_pT->at(i_175), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_HLT_pT_uncorr->Fill(trig_pT->at(i_175), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_pT->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
                    h_pT_uncorr->Fill(p_T->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_eta->Fill(eta->at(i_ele), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                }

            }// End of i_ele iteration

//            prescale_alt = analyzer->getPrescale(trig_fired->at(i_highest)+1);
//            h_HLT_pT->Fill(trig_pT->at(i_highest), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
//            h_HLT_pT_uncorr->Fill(trig_pT->at(i_highest), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//            h_pT->Fill(p_T->at(trig_matched->at(i_highest)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight / prescale_alt);
//            h_pT_uncorr->Fill(p_T->at(trig_matched->at(i_highest)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
//            h_eta->Fill(eta->at(trig_matched->at(i_highest)), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

            if (DEBUG == kFALSE) bar.Draw(i);

        }// End of event iteration

        if(Mgr.isMC == kTRUE) printf("\tNormalization factor: %.8f\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);

        f->cd();
        cout << "\tWriting into file...";

        h_pT_uncorr->Write();
        h_pT->Write();
        h_HLT_pT_uncorr->Write();
        h_HLT_pT->Write();
        h_eta->Write();
        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _SinglePhoton_B) break;

        f->Close();
        if (!f->IsOpen()) cout << "File " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " has been closed successfully.\n" << endl;
        else cout << "FILE " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    cout << "Total weighted events with 2 or more denominator electrons in MC: " << n2MC << endl;
    cout << "Total events with 2 or more denominator electrons in Data: " << n2Data << endl;

    cout << "HLT_Photon175 trigger was activated " << n175 << " times" << endl;
    cout << "HLT_Photon120+HLT_Photon175 trigger was activated " << n120175 << " times" << endl;
    cout << "HLT_Photon90+HLT_Photon175 trigger was activated " << n90175 << " times" << endl;
    cout << "HLT_Photon75+HLT_Photon175 trigger was activated " << n75175 << " times" << endl;
    cout << "HLT_Photon50+HLT_Photon175 trigger was activated " << n50175 << " times" << endl;
    cout << "HLT_Photon36+HLT_Photon175 trigger was activated " << n36175 << " times" << endl;
    cout << "HLT_Photon30+HLT_Photon175 trigger was activated " << n30175 << " times" << endl;
    cout << "HLT_Photon22+HLT_Photon175 trigger was activated " << n22175 << " times" << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of E_FR_HistMaker()


/// -------------------------------- Muon Channel ------------------------------------ ///
void Mu_FR_HistMaker (Bool_t DEBUG)
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

    UInt_t n2MC=0, n2Data=0;

    for (Process_t pr=_DY_10to50; pr<_EndOf_SingleMuon_Normal; pr=next(pr))
    {
        Mgr.SetProc(pr);

        // -- Output ROOTFile -- //
        f = new TFile(Dir+"FR_Hist_Mu_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root", "RECREATE");

        cout << "===========================================================" << endl;
        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Xsec: " << Mgr.Xsec[0] << endl;
        cout << "Wsum: " << Mgr.Wsum[0] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "Directory: " << Dir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        DYAnalyzer *analyzer = new DYAnalyzer("Mu50");

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X(Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root");

        // -- For efficiency SF -- //
        analyzer->SetupEfficiencyScaleFactor_BtoF();
        analyzer->SetupEfficiencyScaleFactor_GtoH();

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D* h_pT_barrel_nume = new TH1D("h_pT_barrel_nume", "h_pT_barrel_nume", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_nume->Sumw2();
        TH1D* h_pT_endcap_nume = new TH1D("h_pT_endcap_nume", "h_pT_endcap_nume", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_nume->Sumw2();
        TH1D* h_pT_barrel_deno = new TH1D("h_pT_barrel_deno", "h_pT_barrel_deno", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno->Sumw2();
        TH1D* h_pT_endcap_deno = new TH1D("h_pT_endcap_deno", "h_pT_endcap_deno", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno->Sumw2();
        TH1D* h_pT_barrel_ctrl = new TH1D("h_pT_barrel_ctrl", "h_pT_barrel_ctrl", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl->Sumw2();
        TH1D* h_pT_endcap_ctrl = new TH1D("h_pT_endcap_ctrl", "h_pT_endcap_ctrl", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl->Sumw2();
        TH1D* h_eta_nume = new TH1D("h_eta_nume", "h_eta_nume", 48, -2.4, 2.4); h_eta_nume->Sumw2();
        TH1D* h_eta_deno = new TH1D("h_eta_deno", "h_eta_deno", 48, -2.4, 2.4); h_eta_deno->Sumw2();
        TH1D* h_eta_ctrl = new TH1D("h_eta_ctrl", "h_eta_ctrl", 48, -2.4, 2.4); h_eta_ctrl->Sumw2();
        TH1D* h_PFiso_barrel_nume = new TH1D("h_PFiso_barrel_nume", "h_PFiso_barrel_nume", 30, 0, 0.15); h_PFiso_barrel_nume->Sumw2();
        TH1D* h_PFiso_endcap_nume = new TH1D("h_PFiso_endcap_nume", "h_PFiso_endcap_nume", 30, 0, 0.15); h_PFiso_endcap_nume->Sumw2();
        TH1D* h_PFiso_barrel_deno = new TH1D("h_PFiso_barrel_deno", "h_PFiso_barrel_deno", 100, 0, 5); h_PFiso_barrel_deno->Sumw2();
        TH1D* h_PFiso_endcap_deno = new TH1D("h_PFiso_endcap_deno", "h_PFiso_endcap_deno", 100, 0, 5); h_PFiso_endcap_deno->Sumw2();
        TH1D* h_PFiso_barrel_ctrl = new TH1D("h_PFiso_barrel_ctrl", "h_PFiso_barrel_ctrl", 100, 0, 5); h_PFiso_barrel_ctrl->Sumw2();
        TH1D* h_PFiso_endcap_ctrl = new TH1D("h_PFiso_endcap_ctrl", "h_PFiso_endcap_ctrl", 100, 0, 5); h_PFiso_endcap_ctrl->Sumw2();
        TH1D* h_TRKiso_barrel_nume = new TH1D("h_TRKiso_barrel_nume", "h_TRKiso_barrel_nume", 50, 0, 1); h_TRKiso_barrel_nume->Sumw2();
        TH1D* h_TRKiso_endcap_nume = new TH1D("h_TRKiso_endcap_nume", "h_TRKiso_endcap_nume", 50, 0, 1); h_TRKiso_endcap_nume->Sumw2();
        TH1D* h_TRKiso_barrel_deno = new TH1D("h_TRKiso_barrel_deno", "h_TRKiso_barrel_deno", 100, 0, 5); h_TRKiso_barrel_deno->Sumw2();
        TH1D* h_TRKiso_endcap_deno = new TH1D("h_TRKiso_endcap_deno", "h_TRKiso_endcap_deno", 100, 0, 5); h_TRKiso_endcap_deno->Sumw2();
        TH1D* h_TRKiso_barrel_ctrl = new TH1D("h_TRKiso_barrel_ctrl", "h_TRKiso_barrel_ctrl", 50, 0.15, 5); h_TRKiso_barrel_deno->Sumw2();
        TH1D* h_TRKiso_endcap_ctrl = new TH1D("h_TRKiso_endcap_ctrl", "h_TRKiso_endcap_ctrl", 50, 0.15, 5); h_TRKiso_endcap_deno->Sumw2();
        TH1D* h_MET = new TH1D("h_MET", "h_MET", 100, 0, 1000); h_MET->Sumw2();
        TH1D* h_MT_barrel_nume = new TH1D("h_MT_barrel_nume", "h_MT_barrel_nume", 500, 0, 1000); h_MT_barrel_nume->Sumw2();
        TH1D* h_MT_endcap_nume = new TH1D("h_MT_endcap_nume", "h_MT_endcap_nume", 500, 0, 1000); h_MT_endcap_nume->Sumw2();
        TH1D* h_MT_barrel_deno = new TH1D("h_MT_barrel_deno", "h_MT_barrel_deno", 500, 0, 1000); h_MT_barrel_deno->Sumw2();
        TH1D* h_MT_endcap_deno = new TH1D("h_MT_endcap_deno", "h_MT_endcap_deno", 500, 0, 1000); h_MT_endcap_deno->Sumw2();
        TH1D* h_MT_barrel_ctrl = new TH1D("h_MT_barrel_ctrl", "h_MT_barrel_ctrl", 500, 0, 1000); h_MT_barrel_ctrl->Sumw2();
        TH1D* h_MT_endcap_ctrl = new TH1D("h_MT_endcap_ctrl", "h_MT_endcap_ctrl", 500, 0, 1000); h_MT_endcap_ctrl->Sumw2();
        TH1D* h_nVTX = new TH1D("h_nVTX", "h_nVTX", 50, 0, 50); h_nVTX->Sumw2();

        TH1D* h_PFiso_barrel_nume_50to70   = new TH1D("h_PFiso_barrel_nume_50to70",   "h_PFiso_barrel_nume_50to70",   15, 0, 0.15); h_PFiso_barrel_nume_50to70  ->Sumw2();
        TH1D* h_PFiso_endcap_nume_50to70   = new TH1D("h_PFiso_endcap_nume_50to70",   "h_PFiso_endcap_nume_50to70",   15, 0, 0.15); h_PFiso_endcap_nume_50to70  ->Sumw2();
        TH1D* h_PFiso_barrel_deno_50to70   = new TH1D("h_PFiso_barrel_deno_50to70",   "h_PFiso_barrel_deno_50to70",   50, 0, 5);    h_PFiso_barrel_deno_50to70  ->Sumw2();
        TH1D* h_PFiso_endcap_deno_50to70   = new TH1D("h_PFiso_endcap_deno_50to70",   "h_PFiso_endcap_deno_50to70",   50, 0, 5);    h_PFiso_endcap_deno_50to70  ->Sumw2();
        TH1D* h_PFiso_barrel_ctrl_50to70   = new TH1D("h_PFiso_barrel_ctrl_50to70",   "h_PFiso_barrel_ctrl_50to70",   50, 0.15, 5);    h_PFiso_barrel_ctrl_50to70  ->Sumw2();
        TH1D* h_PFiso_endcap_ctrl_50to70   = new TH1D("h_PFiso_endcap_ctrl_50to70",   "h_PFiso_endcap_ctrl_50to70",   50, 0.15, 5);    h_PFiso_endcap_ctrl_50to70  ->Sumw2();
        TH1D* h_PFiso_barrel_nume_70to100  = new TH1D("h_PFiso_barrel_nume_70to100",  "h_PFiso_barrel_nume_70to100",  15, 0, 0.15); h_PFiso_barrel_nume_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap_nume_70to100  = new TH1D("h_PFiso_endcap_nume_70to100",  "h_PFiso_endcap_nume_70to100",  15, 0, 0.15); h_PFiso_endcap_nume_70to100 ->Sumw2();
        TH1D* h_PFiso_barrel_deno_70to100  = new TH1D("h_PFiso_barrel_deno_70to100",  "h_PFiso_barrel_deno_70to100",  50, 0, 5);    h_PFiso_barrel_deno_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap_deno_70to100  = new TH1D("h_PFiso_endcap_deno_70to100",  "h_PFiso_endcap_deno_70to100",  50, 0, 5);    h_PFiso_endcap_deno_70to100 ->Sumw2();
        TH1D* h_PFiso_barrel_ctrl_70to100  = new TH1D("h_PFiso_barrel_ctrl_70to100",  "h_PFiso_barrel_ctrl_70to100",  50, 0.15, 5);    h_PFiso_barrel_ctrl_70to100 ->Sumw2();
        TH1D* h_PFiso_endcap_ctrl_70to100  = new TH1D("h_PFiso_endcap_ctrl_70to100",  "h_PFiso_endcap_ctrl_70to100",  50, 0.15, 5);    h_PFiso_endcap_ctrl_70to100 ->Sumw2();
        TH1D* h_PFiso_barrel_nume_100to500 = new TH1D("h_PFiso_barrel_nume_100to500", "h_PFiso_barrel_nume_100to500", 15, 0, 0.15); h_PFiso_barrel_nume_100to500->Sumw2();
        TH1D* h_PFiso_endcap_nume_100to500 = new TH1D("h_PFiso_endcap_nume_100to500", "h_PFiso_endcap_nume_100to500", 15, 0, 0.15); h_PFiso_endcap_nume_100to500->Sumw2();
        TH1D* h_PFiso_barrel_deno_100to500 = new TH1D("h_PFiso_barrel_deno_100to500", "h_PFiso_barrel_deno_100to500", 50, 0, 5);    h_PFiso_barrel_deno_100to500->Sumw2();
        TH1D* h_PFiso_endcap_deno_100to500 = new TH1D("h_PFiso_endcap_deno_100to500", "h_PFiso_endcap_deno_100to500", 50, 0, 5);    h_PFiso_endcap_deno_100to500->Sumw2();
        TH1D* h_PFiso_barrel_ctrl_100to500 = new TH1D("h_PFiso_barrel_ctrl_100to500", "h_PFiso_barrel_ctrl_100to500", 50, 0.15, 5);    h_PFiso_barrel_ctrl_100to500->Sumw2();
        TH1D* h_PFiso_endcap_ctrl_100to500 = new TH1D("h_PFiso_endcap_ctrl_100to500", "h_PFiso_endcap_ctrl_100to500", 50, 0.15, 5);    h_PFiso_endcap_ctrl_100to500->Sumw2();

        TH1D* h_pT_barrel_nume_50to70   = new TH1D("h_pT_barrel_nume_50to70 ",  "h_pT_barrel_nume_50to70 ",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_nume_50to70  ->Sumw2();
        TH1D* h_pT_endcap_nume_50to70   = new TH1D("h_pT_endcap_nume_50to70 ",  "h_pT_endcap_nume_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_nume_50to70  ->Sumw2();
        TH1D* h_pT_barrel_deno_50to70   = new TH1D("h_pT_barrel_deno_50to70 ",  "h_pT_barrel_deno_50to70 ",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno_50to70  ->Sumw2();
        TH1D* h_pT_endcap_deno_50to70   = new TH1D("h_pT_endcap_deno_50to70 ",  "h_pT_endcap_deno_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno_50to70  ->Sumw2();
        TH1D* h_pT_barrel_ctrl_50to70   = new TH1D("h_pT_barrel_ctrl_50to70 ",  "h_pT_barrel_ctrl_50to70 ",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl_50to70  ->Sumw2();
        TH1D* h_pT_endcap_ctrl_50to70   = new TH1D("h_pT_endcap_ctrl_50to70 ",  "h_pT_endcap_ctrl_50to70 ",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl_50to70  ->Sumw2();
        TH1D* h_pT_barrel_nume_70to100  = new TH1D("h_pT_barrel_nume_70to100",  "h_pT_barrel_nume_70to100",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_nume_70to100 ->Sumw2();
        TH1D* h_pT_endcap_nume_70to100  = new TH1D("h_pT_endcap_nume_70to100",  "h_pT_endcap_nume_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_nume_70to100 ->Sumw2();
        TH1D* h_pT_barrel_deno_70to100  = new TH1D("h_pT_barrel_deno_70to100",  "h_pT_barrel_deno_70to100",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno_70to100 ->Sumw2();
        TH1D* h_pT_endcap_deno_70to100  = new TH1D("h_pT_endcap_deno_70to100",  "h_pT_endcap_deno_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno_70to100 ->Sumw2();
        TH1D* h_pT_barrel_ctrl_70to100  = new TH1D("h_pT_barrel_ctrl_70to100",  "h_pT_barrel_ctrl_70to100",  nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl_70to100 ->Sumw2();
        TH1D* h_pT_endcap_ctrl_70to100  = new TH1D("h_pT_endcap_ctrl_70to100",  "h_pT_endcap_ctrl_70to100",  nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl_70to100 ->Sumw2();
        TH1D* h_pT_barrel_nume_100to500 = new TH1D("h_pT_barrel_nume_100to500", "h_pT_barrel_nume_200to500", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_nume_100to500->Sumw2();
        TH1D* h_pT_endcap_nume_100to500 = new TH1D("h_pT_endcap_nume_100to500", "h_pT_endcap_nume_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_nume_100to500->Sumw2();
        TH1D* h_pT_barrel_deno_100to500 = new TH1D("h_pT_barrel_deno_100to500", "h_pT_barrel_deno_200to500", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_deno_100to500->Sumw2();
        TH1D* h_pT_endcap_deno_100to500 = new TH1D("h_pT_endcap_deno_100to500", "h_pT_endcap_deno_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_deno_100to500->Sumw2();
        TH1D* h_pT_barrel_ctrl_100to500 = new TH1D("h_pT_barrel_ctrl_100to500", "h_pT_barrel_ctrl_200to500", nPtBinBarrel, analyzer->ptbin_barrel); h_pT_barrel_ctrl_100to500->Sumw2();
        TH1D* h_pT_endcap_ctrl_100to500 = new TH1D("h_pT_endcap_ctrl_100to500", "h_pT_endcap_ctrl_200to500", nPtBinEndcap, analyzer->ptbin_endcap); h_pT_endcap_ctrl_100to500->Sumw2();

        std::vector<double> *p_T = new std::vector<double>;
        std::vector<double> *eta = new std::vector<double>;
        std::vector<double> *phi = new std::vector<double>;
        std::vector<int> *charge = new std::vector<int>;
        std::vector<double> *relPFiso = new std::vector<double>;
        std::vector<double> *TRKiso = new std::vector<double>;
        Double_t MET_pT, MET_phi, MET_sumEt;
        Int_t nPU;
        Int_t nVTX;
        Double_t PVz;
//        Double_t evt_weight;
        Double_t gen_weight, top_weight;
        Double_t prefiring_weight, prefiring_weight_up, prefiring_weight_down;

        TChain *chain = new TChain("FRTree");

        chain->Add(Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root");
        if (DEBUG == kTRUE) cout << Dir+"SelectedForFR_Mu_"+Mgr.Procname[Mgr.CurrentProc]+".root" << endl;
        chain->SetBranchStatus("p_T", 1);
        chain->SetBranchStatus("eta", 1);
        chain->SetBranchStatus("phi", 1);
        chain->SetBranchStatus("relPFiso", 1);
        chain->SetBranchStatus("TRKiso", 1);
        chain->SetBranchStatus("MET_pT", 1);
        chain->SetBranchStatus("MET_phi", 1);
        chain->SetBranchStatus("MET_sumEt", 1);
        chain->SetBranchStatus("nPU", 1);
        chain->SetBranchStatus("nVTX", 1);
        chain->SetBranchStatus("PVz", 1);
//        chain->SetBranchStatus("evt_weight", 1);
        chain->SetBranchStatus("gen_weight", 1);
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
        chain->SetBranchAddress("MET_pT", &MET_pT);
        chain->SetBranchAddress("MET_phi", &MET_phi);
        chain->SetBranchAddress("MET_sumEt", &MET_sumEt);
        chain->SetBranchAddress("nPU", &nPU);
        chain->SetBranchAddress("nVTX", &nVTX);
        chain->SetBranchAddress("PVz", &PVz);
//        chain->SetBranchAddress("evt_weight", &evt_weight);
        chain->SetBranchAddress("gen_weight", &gen_weight);
        chain->SetBranchAddress("top_weight", &top_weight);
        chain->SetBranchAddress("prefiring_weight", &prefiring_weight);
        chain->SetBranchAddress("prefiring_weight_up", &prefiring_weight_up);
        chain->SetBranchAddress("prefiring_weight_down", &prefiring_weight_down);

        Int_t NEvents = chain->GetEntries();
        cout << "\t[Sum of weights: " << Mgr.Wsum[0] << "]" << endl;
        cout << "\t[Number of events: " << NEvents << "]" << endl;

        if (DEBUG == kTRUE) NEvents = 100;

        myProgressBar_t bar(NEvents);

        for(Int_t i=0; i<NEvents; i++)
        {
            chain->GetEntry(i);
            if (DEBUG == kTRUE){
                cout << "Evt " << i << endl;
                cout << "nMuons = " << p_T->size() << endl;
                cout << "p_T[1] = " << p_T->at(0) << endl;
                cout << "eta[1] = " << eta->at(0) << endl;
                cout << "phi[1] = " << phi->at(0) << endl;
            }

            // -- Pileup-Reweighting -- //
            Double_t PUWeight = 1;
            if(Mgr.isMC == kTRUE) PUWeight = analyzer->PileUpWeightValue_80X(nPU);
            if (DEBUG == kTRUE) cout << "PU weight " << PUWeight << endl;

            // -- efficiency weights -- //
            Double_t weight1 = 0, weight2 = 0, effweight = 1;

            // -- PVz weights -- //
            Double_t PVzWeight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) PVzWeight = analyzer->PVzWeightValue(PVz);

            // -- L1 prefiring weights -- //
            Double_t L1weight = 1;
            if (Mgr.isMC == kTRUE && !Mgr.Tag[0].Contains("QCD")) L1weight = prefiring_weight;

            // -- Top pT weights -- //
            Double_t TopPtWeight = 1;
            if (Mgr.isMC == kTRUE && Mgr.Tag[0].Contains("ttbar")) TopPtWeight = top_weight;

            // -- Normalization -- //
            Double_t TotWeight = gen_weight;
            if (Mgr.isMC == kTRUE) TotWeight = (Lumi * Mgr.Xsec[0] / Mgr.Wsum[0]) * gen_weight;
            if (DEBUG == kTRUE) cout << "Total weight " << TotWeight << endl;

            if (Mgr.isMC == kTRUE && p_T->size() > 1) n2MC += TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight;
            if (Mgr.isMC == kFALSE && p_T->size() > 1) n2Data++;

            // For finding the leading muon
            TLorentzVector mu_lead;
            mu_lead.SetPtEtaPhiM(0, 0, 0, M_Mu);
            Double_t iso_lead = -9999;

            Double_t iso_count = 0;
            for (UInt_t i_mu=0; i_mu<p_T->size(); i_mu++)
            {
                if (p_T->at(i_mu) <= 52) continue;
                if (relPFiso->at(i_mu) < 0.15) iso_count++;

                // Selecting leading muon (could also try finding a muon with the best isolation)
                if (p_T->at(i_mu) > mu_lead.Pt())
                {
                    mu_lead.SetPtEtaPhiM(p_T->at(i_mu), eta->at(i_mu), phi->at(i_mu), M_Mu);
                    iso_lead = relPFiso->at(i_mu);
                }
            }
            if (iso_count > 1) continue;
            Double_t dTheta = mu_lead.Phi() - MET_phi;
            Double_t MT = sqrt(2 * mu_lead.Pt() * MET_pT * (1 - cos(dTheta)));
//            if (MT >= 60) continue;
//            if (MET_pT >= 50) continue;

            h_nVTX->Fill(nVTX, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            h_MET->Fill(MET_pT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

            for (UInt_t i_mu=0; i_mu<p_T->size(); i_mu++)
            {
                if (p_T->at(i_mu) != p_T->at(i_mu))
                {
                    cout << p_T->at(i_mu) << " " << eta->at(i_mu) << " " << phi->at(i_mu) << " " << charge->at(i_mu) << " " << relPFiso->at(i_mu) << endl;
                    continue;
                }
                if (p_T->at(i_mu) <= 52) continue;
                if (DEBUG == kTRUE) cout << "i_mu = " << i_mu << endl;

                // -- Efficiency scale factor -- //
//                if(Mgr.isMC == kTRUE)
//                {
//                    weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF_new(MuMu);
//                    weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH_new(MuMu);
//                    effweight = (Lumi_BtoF * weight1 + Lumi_GtoH * weight2) / Lumi;
//                }

                if (relPFiso->at(i_mu) < 0.15) // Signal/Numerator
                {
                    h_eta_nume->Fill(eta->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(eta->at(i_mu)) < 1.2) // Barrel
                    {
                        h_pT_barrel_nume->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_barrel_nume->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_TRKiso_barrel_nume->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                        if (p_T->at(i_mu) < 70)
                        {
                            h_pT_barrel_nume_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_nume_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (p_T->at(i_mu) < 100)
                        {
                            h_pT_barrel_nume_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_nume_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else //if (p_T->at(i_mu) < 500)
                        {
                            h_pT_barrel_nume_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_nume_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                    }
                    else // Endcap
                    {
                        h_pT_endcap_nume->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap_nume->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_TRKiso_endcap_nume->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                        if (p_T->at(i_mu) < 70)
                        {
                            h_pT_endcap_nume_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_nume_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (p_T->at(i_mu) < 100)
                        {
                            h_pT_endcap_nume_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_nume_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else //if (p_T->at(i_mu) < 500)
                        {
                            h_pT_endcap_nume_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_nume_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                    }
                } // End of if(Signal/Numerator)
                else // Control
                {
                    h_eta_ctrl->Fill(eta->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    if (fabs(eta->at(i_mu)) < 1.2) // Barrel
                    {
                        h_pT_barrel_ctrl->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_barrel_ctrl->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_TRKiso_barrel_ctrl->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                        if (p_T->at(i_mu) < 70)
                        {
                            h_pT_barrel_ctrl_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_ctrl_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (p_T->at(i_mu) < 100)
                        {
                            h_pT_barrel_ctrl_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_ctrl_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else //if (p_T->at(i_mu) < 500)
                        {
                            h_pT_barrel_ctrl_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_barrel_ctrl_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                    }
                    else // Endcap
                    {
                        h_pT_endcap_ctrl->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap_ctrl->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_TRKiso_endcap_ctrl->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                        if (p_T->at(i_mu) < 70)
                        {
                            h_pT_endcap_ctrl_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_ctrl_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else if (p_T->at(i_mu) < 100)
                        {
                            h_pT_endcap_ctrl_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_ctrl_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                        else //if (p_T->at(i_mu) < 500)
                        {
                            h_pT_endcap_ctrl_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                            h_PFiso_endcap_ctrl_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        }
                    }
                }// End of if(Control)
                // Denominator
                h_eta_deno->Fill(eta->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                if (fabs(eta->at(i_mu)) < 1.2) // Barrel
                {
                    h_pT_barrel_deno->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_PFiso_barrel_deno->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_TRKiso_barrel_deno->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                    if (p_T->at(i_mu) < 70)
                    {
                        h_pT_barrel_deno_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_barrel_deno_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else if (p_T->at(i_mu) < 100)
                    {
                        h_pT_barrel_deno_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_barrel_deno_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else //if (p_T->at(i_mu) < 500)
                    {
                        h_pT_barrel_deno_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_barrel_deno_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                }
                else // Endcap
                {
                    h_pT_endcap_deno->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_PFiso_endcap_deno->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    h_TRKiso_endcap_deno->Fill(TRKiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);

                    if (p_T->at(i_mu) < 70)
                    {
                        h_pT_endcap_deno_50to70->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap_deno_50to70->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else if (p_T->at(i_mu) < 100)
                    {
                        h_pT_endcap_deno_70to100->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap_deno_70to100->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                    else //if (p_T->at(i_mu) < 500)
                    {
                        h_pT_endcap_deno_100to500->Fill(p_T->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                        h_PFiso_endcap_deno_100to500->Fill(relPFiso->at(i_mu), TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                    }
                }
            }// End of i_mu iteration

            if (fabs(mu_lead.Eta()) < 1.2)
            {
                h_MT_barrel_deno->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                if (iso_lead < 0.15) h_MT_barrel_nume->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                else h_MT_barrel_ctrl->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }
            else
            {
                h_MT_endcap_deno->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                if (iso_lead < 0.15) h_MT_endcap_nume->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
                else h_MT_endcap_ctrl->Fill(MT, TotWeight * PUWeight * effweight * PVzWeight * L1weight * TopPtWeight);
            }

            bar.Draw(i);

        }// End of event iteration

        if(Mgr.isMC == kTRUE) printf("\tNormalization factor: %.8f\n", Lumi*Mgr.Xsec[0]/Mgr.Wsum[0]);

        f->cd();
        cout << "\tWriting into file...";

        h_pT_barrel_nume->Write();
        h_pT_endcap_nume->Write();
        h_pT_barrel_deno->Write();
        h_pT_endcap_deno->Write();
        h_pT_barrel_ctrl->Write();
        h_pT_endcap_ctrl->Write();
        h_eta_nume->Write();
        h_eta_deno->Write();
        h_eta_ctrl->Write();
        h_PFiso_barrel_nume->Write();
        h_PFiso_endcap_nume->Write();
        h_PFiso_barrel_deno->Write();
        h_PFiso_endcap_deno->Write();
        h_PFiso_barrel_ctrl->Write();
        h_PFiso_endcap_ctrl->Write();
        h_TRKiso_barrel_nume->Write();
        h_TRKiso_endcap_nume->Write();
        h_TRKiso_barrel_deno->Write();
        h_TRKiso_endcap_deno->Write();
        h_TRKiso_barrel_ctrl->Write();
        h_TRKiso_endcap_ctrl->Write();
        h_MET->Write();
        h_MT_barrel_nume->Write();
        h_MT_endcap_nume->Write();
        h_MT_barrel_deno->Write();
        h_MT_endcap_deno->Write();
        h_MT_barrel_ctrl->Write();
        h_MT_endcap_ctrl->Write();
        h_nVTX->Write();

        h_pT_barrel_nume_50to70->Write();
        h_pT_endcap_nume_50to70->Write();
        h_pT_barrel_deno_50to70->Write();
        h_pT_endcap_deno_50to70->Write();
        h_pT_barrel_ctrl_50to70->Write();
        h_pT_endcap_ctrl_50to70->Write();
        h_pT_barrel_nume_70to100->Write();
        h_pT_endcap_nume_70to100->Write();
        h_pT_barrel_deno_70to100->Write();
        h_pT_endcap_deno_70to100->Write();
        h_pT_barrel_ctrl_70to100->Write();
        h_pT_endcap_ctrl_70to100->Write();
        h_pT_barrel_nume_100to500->Write();
        h_pT_endcap_nume_100to500->Write();
        h_pT_barrel_deno_100to500->Write();
        h_pT_endcap_deno_100to500->Write();
        h_pT_barrel_ctrl_100to500->Write();
        h_pT_endcap_ctrl_100to500->Write();

        h_PFiso_barrel_nume_50to70->Write();
        h_PFiso_endcap_nume_50to70->Write();
        h_PFiso_barrel_deno_50to70->Write();
        h_PFiso_endcap_deno_50to70->Write();
        h_PFiso_barrel_ctrl_50to70->Write();
        h_PFiso_endcap_ctrl_50to70->Write();
        h_PFiso_barrel_nume_70to100->Write();
        h_PFiso_endcap_nume_70to100->Write();
        h_PFiso_barrel_deno_70to100->Write();
        h_PFiso_endcap_deno_70to100->Write();
        h_PFiso_barrel_ctrl_70to100->Write();
        h_PFiso_endcap_ctrl_70to100->Write();
        h_PFiso_barrel_nume_100to500->Write();
        h_PFiso_endcap_nume_100to500->Write();
        h_PFiso_barrel_deno_100to500->Write();
        h_PFiso_endcap_deno_100to500->Write();
        h_PFiso_barrel_ctrl_100to500->Write();
        h_PFiso_endcap_ctrl_100to500->Write();

        cout << " Finished.\n" << endl;

        if (DEBUG == kTRUE && pr == _DY_10to50) break;
        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _QCDMuEnriched_1000toInf) pr = _EndOf_DoubleEG_Normal; // next -- SingleMuon_B

        f->Close();
        if (!f->IsOpen()) cout << "File " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " has been closed successfully.\n" << endl;
        else cout << "FILE " << Dir+"FR_Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root" << " COULD NOT BE CLOSED!\n" << endl;
        cout << "===========================================================\n" << endl;
    } // End of pr iteration

    cout << "Total weighted events with 2 or more denominator muons in MC: " << n2MC << endl;
    cout << "Total events with 2 or more denominator muons in Data: " << n2Data << endl;

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of Mu_FR_HistMaker()


/// -------------------------------- Muon Channel ------------------------------------ ///
void Mu_QCD_HistMaker (Bool_t DEBUG, Int_t type=1)
// type=0 uses files with Mu50 trigger
// type=1 -- no trigger
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

} // End of Mu_QCD_HistMaker()


void Mu_WJET_HistMaker (Bool_t DEBUG, Int_t type=2)
// type=0 uses files with Mu50 trigger
// type=1 -- IsoMu24_OR_IsoTkMu24
// type=2 -- no trigger
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

} // End of Mu_WJET_HistMaker()
